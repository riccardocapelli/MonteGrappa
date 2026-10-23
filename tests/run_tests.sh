#!/bin/bash
#
# MonteGrappa test suite
#
#   tests/run_tests.sh [options] [group ...]
#
# groups (default: all the fast ones)
#   errors      wrong input must stop with a message, not crash
#   regression  fixed-seed runs compared with tests/reference (numerical tolerance,
#               so that references made on one platform hold on the others)
#   echeck      every move type on several systems with the keyword echeck: the
#               incremental energy, the backup copy and the cell list must always
#               agree with a full recalculation
#   restart     a run restarted from its checkpoint must give identical files
#   cells       cell list on and off must give identical trajectories
#   tools       grappino and mgp2pdb
#   mpi         MPI executable run on a single rank (if it was built)
#   stat        (slow, not in the default) detailed balance: the BGS local move
#               must sample the same mean energy as pivot moves
#
# options
#   --asan              use bin/montegrappa_asan and fail on sanitizer reports
#   --update-reference  rewrite tests/reference with the current executable
#   --keep              keep the working directory
#
# Exit status: number of failed tests (0 = all passed).

ROOT=$(cd "$(dirname "$0")/.." && pwd)
DATA=$ROOT/tests/data
REF=$ROOT/tests/reference
MG=$ROOT/bin/montegrappa
ASAN=0; UPDATE=0; KEEP=0; TGROUPS=""

for a in "$@"; do
	case $a in
		--asan) ASAN=1; MG=$ROOT/bin/montegrappa_asan ;;
		--update-reference) UPDATE=1 ;;
		--keep) KEEP=1 ;;
		-h|--help) sed -n '2,30p' "$0"; exit 0 ;;
		*) TGROUPS="$TGROUPS $a" ;;
	esac
done
[ -z "$TGROUPS" ] && TGROUPS="errors regression echeck restart cells tools mpi"

[ -x "$MG" ] || { echo "executable $MG not found: build it first"; exit 1; }

WORK=$(mktemp -d "${TMPDIR:-/tmp}/mgtests.XXXXXX")
NPASS=0; NFAIL=0; FAILED=""

pass() { NPASS=$((NPASS+1)); printf "  PASS  %s\n" "$1"; }
fail() { NFAIL=$((NFAIL+1)); FAILED="$FAILED\n  $1"; printf "  FAIL  %s: %s\n" "$1" "$2"; }

# write a parameter file: par <file> <lines with \n>
par() { printf "$2" > "$1"; }

# run montegrappa in a new directory: mg <dir> <pol> <pot> <par lines> ; returns its exit status
mg() {
	local d=$WORK/$1
	rm -rf "$d"; mkdir -p "$d"
	cp "$DATA/$2" "$DATA/$3" "$d/"
	[ -f "$DATA/polyala.op" ] && cp "$DATA/polyala.op" "$d/"
	par "$d/run.par" "$4"
	( cd "$d" && "$MG" "$2" "$3" run.par > stdout.log 2> stderr.log )
}

# sanitizer reports or a crash in a run directory
sanitizer_ok() {
	! grep -q -E "ERROR: AddressSanitizer|runtime error:" "$1/stderr.log"
}

# tracked and recomputed energy at the end of a run must agree
final_ok() {
	awk '/^Final:/ { split($2,a,"="); split($3,b,"="); n++; d=a[2]-b[2]; if (d<0) d=-d;
	                 if (d > 1e-6*(1+(b[2]<0?-b[2]:b[2]))) bad=1 }
	     END { exit (n==0 || bad) }' "$1/stderr.log"
}

# numerical comparison of two text files, token by token
same_numbers() {
	awk -v atol=1e-4 -v rtol=1e-6 '
		FNR==NR { a[NR]=$0; na=NR; next }
		{ nb=FNR; if (!(FNR in a)) { bad=1; exit }
		  n1=split(a[FNR],x); n2=split($0,y); if (n1!=n2) { bad=1; exit }
		  for (i=1;i<=n1;i++) {
			if (x[i] ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/ && y[i] ~ /^[-+]?[0-9.]+([eE][-+]?[0-9]+)?$/) {
				d=x[i]-y[i]; if (d<0) d=-d; s=(y[i]<0?-y[i]:y[i]);
				if (d > atol + rtol*s) { bad=1; exit }
			} else if (x[i]!=y[i]) { bad=1; exit }
		  } }
		END { exit (bad || na!=nb) }' "$1" "$2"
}

######################################################################
group_errors() {
	echo "errors"
	mg err_temp polyala.pol polyala.pot "nchains 1\nnstep 10\npivot 1\n"
	if [ $? -ne 0 ] && grep -q "Temp" $WORK/err_temp/stderr.log; then pass "missing Temp"; else fail "missing Temp" "not refused"; fi

	mg err_opminim polyala.pol polyala.pot "nchains 1\nTemp 1\nnstep 10\npivot 1\nop_minim steep\n"
	if [ $? -ne 0 ] && grep -q "op_minim" $WORK/err_opminim/stderr.log; then pass "unknown op_minim"; else fail "unknown op_minim" "not refused"; fi

	mg err_nchains polyala.pol polyala.pot "nchains 2\nTemp 1\nnstep 10\npivot 1\n"
	if [ $? -ne 0 ] && grep -q "chain" $WORK/err_nchains/stderr.log && sanitizer_ok $WORK/err_nchains; then pass "wrong nchains"; else fail "wrong nchains" "not refused"; fi

	local long=$(printf 'x%.0s' $(seq 1 300))
	( cd $WORK && "$MG" "$long.pol" polyala.pot run.par > /dev/null 2> long.log )
	if [ $? -ne 0 ] && grep -q "too long" $WORK/long.log; then pass "file name too long"; else fail "file name too long" "not refused"; fi

	mkdir -p $WORK/err_ckp; cp $DATA/polyala.pot $WORK/err_ckp/
	par $WORK/err_ckp/run.par "nchains 1\nTemp 1\nnstep 10\npivot 1\n"
	( cd $WORK/err_ckp && "$MG" checkpoint polyala.pot run.par > /dev/null 2> stderr.log )
	if [ $? -ne 0 ] && grep -q "checkpoint" $WORK/err_ckp/stderr.log; then pass "missing checkpoint"; else fail "missing checkpoint" "not refused"; fi
}

######################################################################
# name|pol|pot|parameter lines
REGRESSION_CASES="
ca|ca.pol|ca.pot|nchains 1\nTemp 1.0\nflip 1\npivot 4\nmpivot 1\nnstep 20000\n
polyala|polyala.pol|polyala.pot|nchains 1\nTemp 1.0\nflip 3\npivot 4\nmpivot 1\nsidechain 2\nnstep 20000\n
three|three.pol|homo.pot|nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\nflip 1\npivot 5\nmovecom 5\nmoverot 5\ncomcluster 10\nrotcluster 10\nnstep 20000\n
two|two.pol|two_hf.pot|nchains 2\nTemp 2.0\nflip 1\npivot 5\nmovecom 5\nmoverot 5\nsidechain 3\nbackrub 5\nbacksiderub 7\nmflip 7\nlpivot 7\nmovebias 4\nnmul_local 5\nnstep 10000\n
hb|polyala.pol|polyala_hb.pot|nchains 1\nTemp 1.0\nhb\nflip 1\npivot 4\nmpivot 2\nbackrub 5\nmovebias 5\nnmul_local 5\nnstep 10000\n
optim|polyala.pol|polyala.pot|nchains 1\nTemp 1.0\nmpivot 1\npivot 4\nnrun 2\nnstep 60000\nalways_restart\nrecord_native\nop_file polyala.op\nop_minim sample\nop_itermax 200\nop_step 0.2\nop_T 1.0\nop_deltat 5000\nop_wait 10000\nop_rw 3.5\nop_r0 1.0\n
"

group_regression() {
	echo "regression"
	echo "$REGRESSION_CASES" | while IFS="|" read -r name pol pot lines; do
		[ -z "$name" ] && continue
		mg reg_$name $pol $pot "seed 11\nnprinttrj 5000\nnprinte 1000\nnprintlog 100000\n$lines"
		rc=$?
		d=$WORK/reg_$name
		outputs=$(cd $d && ls energy*.ene last*.pol newpot_*.pot 2>/dev/null)
		if [ $UPDATE = 1 ]; then
			mkdir -p $REF/$name; rm -f $REF/$name/*
			for f in $outputs; do cp $d/$f $REF/$name/; done
			echo "  reference of $name updated"
			continue
		fi
		if [ $rc -ne 0 ] || ! sanitizer_ok $d; then echo "FAIL" > $d/result; printf "  FAIL  regression %s: run failed\n" $name; continue; fi
		if ! final_ok $d; then echo "FAIL" > $d/result; printf "  FAIL  regression %s: tracked energy differs from recomputed\n" $name; continue; fi
		res=PASS
		for f in $(cd $REF/$name 2>/dev/null && ls); do
			same_numbers $REF/$name/$f $d/$f || { res="FAIL"; printf "  FAIL  regression %s: %s differs from the reference\n" $name $f; }
		done
		[ -d $REF/$name ] || { res=FAIL; printf "  FAIL  regression %s: no reference\n" $name; }
		echo $res > $d/result
		[ $res = PASS ] && printf "  PASS  regression %s\n" $name
	done
	# the loop runs in a subshell: collect its results
	for r in $WORK/reg_*/result; do
		[ -f $r ] || continue
		if [ "$(cat $r)" = PASS ]; then NPASS=$((NPASS+1)); else NFAIL=$((NFAIL+1)); FAILED="$FAILED\n  regression $(basename $(dirname $r) | sed 's/^reg_//')"; fi
	done
}

######################################################################
MOVES="flip pivot mpivot sidechain lpivot mflip movecom movebias moverot comcluster rotcluster backrub backsiderub"

# system|pol|pot|header
ECHECK_SYSTEMS="
polyala|polyala.pol|polyala.pot|nchains 1\nTemp 2.0\n
angles|polyala.pol|polyala_ang.pot|nchains 1\nTemp 2.0\n
hbonds|polyala.pol|polyala_hb.pot|nchains 1\nTemp 1.0\nhb\n
three|three.pol|homo.pot|nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.5\n
two|two.pol|two_hf.pot|nchains 2\nTemp 2.0\n
"

group_echeck() {
	echo "echeck (every move type, also with the cell list forced on)"
	local sys pol pot hdr mv cells extra n bad
	for cells in -1 1; do
		echo "$ECHECK_SYSTEMS" | while IFS="|" read -r sys pol pot hdr; do
			[ -z "$sys" ] && continue
			n=0; bad=""
			for mv in $MOVES; do
				extra="flip 2\n$mv 1\n"; [ $mv = flip ] && extra="flip 1\n"
				mg ec_${sys}_$mv $pol $pot "${hdr}seed 42\nnstep 1500\ncells $cells\necheck\nnmul_local 4\nnprinttrj 100000\nnprinte 100000\nnprintlog 100000\n$extra"
				rc=$?; d=$WORK/ec_${sys}_$mv
				if [ $rc -ne 0 ] || grep -q ECHECK $d/stderr.log || ! sanitizer_ok $d; then bad="$bad $mv"; fi
			done
			if [ -z "$bad" ]; then echo PASS; else echo "FAIL $sys cells=$cells:$bad"; fi > $WORK/ec_result_${sys}_$cells
		done
	done
	for r in $WORK/ec_result_*; do
		name=$(basename $r | sed 's/^ec_result_//')
		if [ "$(cat $r)" = PASS ]; then pass "echeck $name"; else fail "echeck $name" "$(cat $r | sed 's/^FAIL [^:]*://')"; fi
	done
}

######################################################################
# a run interrupted at its checkpoint and restarted must give identical files
restart_case() {
	local name=$1 pol=$2 pot=$3 lines=$4
	mg rst_$name $pol $pot "seed 5\nnstep 4000\nncheckpoint 3000\nnprinte 200\nnprinttrj 1000\nnprintlog 100000\n$lines"
	local a=$WORK/rst_$name b=$WORK/rst_${name}_b
	rm -rf $b; cp -r $a $b
	awk '$1<3000' $a/energy.ene > $b/energy.ene
	awk '{print} /^ENDMDL/{n++; if(n==3) exit}' $a/traj.pdb > $b/traj.pdb
	rm -f $b/last.pol
	( cd $b && "$MG" checkpoint $pot run.par > stdout.log 2> stderr.log )
	local res=ok f
	for f in energy.ene traj.pdb last.pol; do cmp -s $a/$f $b/$f || res="$f differs"; done
	sanitizer_ok $b || res="sanitizer report"
	if [ "$res" = ok ]; then pass "restart $name"; else fail "restart $name" "$res"; fi
}

group_restart() {
	echo "restart"
	restart_case polyala polyala.pol polyala.pot "nchains 1\nTemp 1.0\nmpivot 1\npivot 4\nsidechain 2\nflip 3\nmovebias 5\nnmul_local 5\n"
	restart_case three three.pol homo.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\nflip 1\npivot 5\ncomcluster 5\nrotcluster 5\nmovecom 7\nmoverot 7\n"
	restart_case two two.pol two_hf.pot "nchains 2\nTemp 2.0\nflip 1\npivot 5\nmovecom 5\nmoverot 5\nsidechain 3\nbackrub 5\nmflip 7\nlpivot 7\nmovebias 4\nnmul_local 5\n"
}

######################################################################
cells_case() {
	local name=$1 pol=$2 pot=$3 lines=$4 c
	for c in 0 1; do
		mg cells_${name}_$c $pol $pot "seed 3\nnprinttrj 10000\nnprinte 500\nnprintlog 100000\ncells $c\n$lines"
	done
	if cmp -s $WORK/cells_${name}_0/energy.ene $WORK/cells_${name}_1/energy.ene && cmp -s $WORK/cells_${name}_0/last.pol $WORK/cells_${name}_1/last.pol \
		&& final_ok $WORK/cells_${name}_1 && sanitizer_ok $WORK/cells_${name}_1
	then pass "cells $name"; else fail "cells $name" "trajectories differ"; fi
}

group_cells() {
	echo "cells"
	cells_case chain1100 chain1100.pol chain1100.pot "nchains 1\nnstep 5000\nTemp 1.0\nflip 1\npivot 2\nnosidechain\nnodihpot\nnoangpot\nnohfields\n"
	cells_case three three.pol homo.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\nnstep 20000\nflip 1\npivot 5\ncomcluster 5\nrotcluster 5\nmovecom 7\nmoverot 7\n"
	cells_case two two.pol two_hf.pot "nchains 2\nTemp 2.0\nnstep 5000\nflip 1\npivot 5\nmovecom 5\nmoverot 5\nsidechain 3\nbackrub 5\nmflip 7\nmovebias 4\nnmul_local 5\n"
}

######################################################################
group_tools() {
	echo "tools"
	local d=$WORK/grappino
	if [ -x $ROOT/bin/grappino ]; then
		mkdir -p $d; cp $ROOT/tutorial/tutorial_1/1YPC.pdb $ROOT/tutorial/tutorial_1/1YPC-CA.in $d/
		( cd $d && $ROOT/bin/grappino 1YPC-CA.in > g.log 2>&1 )
		if cmp -s $d/1YPC-CA.pol $DATA/ca.pol && cmp -s $d/1YPC-CA.pot $DATA/ca.pot; then pass "grappino (tutorial 1)"; else fail "grappino (tutorial 1)" "output differs from tests/data/ca.*"; fi
	else
		echo "  SKIP  grappino not built"
	fi
	if [ -x $ROOT/bin/mgp2pdb ]; then
		d=$WORK/mgp2pdb; mkdir -p $d
		( cd $d && $ROOT/bin/mgp2pdb -i $DATA/two.pol -o two.pdb > log 2>&1 )
		if [ "$(grep -c '^ATOM' $d/two.pdb 2>/dev/null)" = 150 ]; then pass "mgp2pdb (two chains)"; else fail "mgp2pdb (two chains)" "wrong number of atoms"; fi
	else
		echo "  SKIP  mgp2pdb not built"
	fi
}

######################################################################
group_mpi() {
	echo "mpi (single rank)"
	local MPI=$ROOT/bin/montegrappa_mpi
	[ $ASAN = 1 ] && { echo "  SKIP  not with --asan"; return; }
	[ -x $MPI ] || { echo "  SKIP  montegrappa_mpi not built"; return; }
	local d=$WORK/mpi; mkdir -p $d; cp $DATA/polyala.pol $DATA/polyala.pot $d/
	par $d/run.par "nchains 1\nseed 5\nnstep 4000\nncheckpoint 3000\nnprinte 200\nnprinttrj 1000\nnprintlog 100000\nmpivot 1\npivot 4\nntemp 1\ntemperatures\n1.0\n"
	( cd $d && $MPI polyala.pol polyala.pot run.par > stdout.log 2> stderr.log )
	if [ $? -ne 0 ] || ! grep -q "^Final" $d/proc_0.log; then fail "mpi run" "did not complete (is MPI usable without mpirun?)"; return; fi
	awk '/^Final:/ { split($2,a,"="); split($3,b,"="); d=a[2]-b[2]; exit (d>1e-6||d<-1e-6) }' $d/proc_0.log && pass "mpi run" || fail "mpi run" "tracked energy differs from recomputed"
	cp -r $d $d.b; awk '$1<3000' $d/energy_proc0.ene > $d.b/energy_proc0.ene; rm -f $d.b/last_proc0.pol
	( cd $d.b && $MPI checkpoint polyala.pot run.par > stdout.log 2> stderr.log )
	if cmp -s $d/energy_proc0.ene $d.b/energy_proc0.ene && cmp -s $d/last_proc0.pol $d.b/last_proc0.pol; then pass "mpi restart"; else fail "mpi restart" "files differ"; fi
}

######################################################################
# Detailed balance: on a 30-residue Go fragment at T=0.6, the BGS local move (plus a
# few pivots) must give the same mean energy as pivot moves alone, within 3 sigma.
group_stat() {
	echo "stat (detailed balance, a few minutes)"
	local setup seed
	for setup in pivot bgs; do
		for seed in 11 22 33 44; do
			case $setup in
				pivot) mv="pivot 1\nmpivot 1\n" ;;
				bgs) mv="movebias 1\nnmul_local 6\npivot 10\n" ;;
			esac
			( mg stat_${setup}_$seed fragment.pol fragment.pot "nchains 1\nseed $seed\nnstep 8000000\nTemp 0.6\n${mv}nprinttrj 100000000\nnprinte 100\nnprintlog 100000000\n" ) &
		done
	done
	wait
	for setup in pivot bgs; do
		for seed in 11 22 33 44; do
			awk '{n++; e[n]=$2} END{s=0; for(i=int(n/2)+1;i<=n;i++) s+=e[i]; print s/(n-int(n/2))}' $WORK/stat_${setup}_$seed/energy.ene
		done | awk -v s=$setup '{m+=$1; m2+=$1*$1; n++} END{mu=m/n; printf "%s %f %f\n", s, mu, sqrt((m2/n-mu*mu)/(n-1))}'
	done > $WORK/stat_summary
	cat $WORK/stat_summary | sed 's/^/        /'
	if awk 'NR==1{m1=$2;s1=$3} NR==2{m2=$2;s2=$3} END{d=m1-m2; if(d<0)d=-d; exit !(d < 3*sqrt(s1*s1+s2*s2))}' $WORK/stat_summary
	then pass "BGS vs pivot mean energy"; else fail "BGS vs pivot mean energy" "difference above 3 sigma"; fi
}

######################################################################
echo "MonteGrappa tests with $MG"
echo "working directory $WORK"
for g in $TGROUPS; do
	case $g in
		errors|regression|echeck|restart|cells|tools|mpi|stat) group_$g ;;
		*) echo "unknown group $g"; NFAIL=$((NFAIL+1)) ;;
	esac
done

[ $KEEP = 1 ] || [ $NFAIL -gt 0 ] || rm -rf $WORK
echo
echo "passed $NPASS, failed $NFAIL"
[ $NFAIL -gt 0 ] && { printf "failed:$FAILED\n"; echo "(outputs kept in $WORK)"; }
exit $NFAIL
