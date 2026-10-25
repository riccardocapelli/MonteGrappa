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
#   observables observables of aggregation on sheets of known structure
#   avbmc       aggregation-volume-bias move: input checks, echeck on periodic
#               systems, and the exact distribution of non-interacting chains
#   mpi_multi   parallel tempering with 4 ranks (only where mpirun works): run,
#               exact restart, and the replica at the lowest temperature must
#               sample the same mean energy as plain runs at that temperature
#   stat        (slow, not in the default) detailed balance: the BGS local move
#               must sample the same mean energy as pivot moves, and so must
#               interacting chains with and without the avbmc move
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
GRAPPINO=$ROOT/bin/grappino
ASAN=0; UPDATE=0; KEEP=0; TGROUPS=""

for a in "$@"; do
	case $a in
		--asan) ASAN=1; MG=$ROOT/bin/montegrappa_asan; [ -x $ROOT/bin/grappino_asan ] && GRAPPINO=$ROOT/bin/grappino_asan
			export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0} ;;		# memory still allocated at exit is not an error here
		--update-reference) UPDATE=1 ;;
		--keep) KEEP=1 ;;
		-h|--help) sed -n '2,30p' "$0"; exit 0 ;;
		*) TGROUPS="$TGROUPS $a" ;;
	esac
done
[ -z "$TGROUPS" ] && TGROUPS="errors regression echeck restart cells tools mpi observables avbmc mpi_multi"

[ -x "$MG" ] || { echo "executable $MG not found: build it first"; exit 1; }

WORK=${MG_TEST_WORKDIR:-$(mktemp -d "${TMPDIR:-/tmp}/mgtests.XXXXXX")}; mkdir -p "$WORK"
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

	mg err_nmul polyala.pol polyala.pot "nchains 1\nTemp 1\nnstep 10\nmovebias 1\n"
	if [ $? -ne 0 ] && grep -q "nmul_local" $WORK/err_nmul/stderr.log; then pass "movebias without nmul_local"; else fail "movebias without nmul_local" "not refused"; fi

	# hydrogen bonds of atoms that are not in the polymer; hydrogen bonds in the pot file without hb
	d=$WORK/err_hbatom; rm -rf $d; mkdir -p $d; cp $DATA/polyala.pol $d/; { cat $DATA/polyala_hb.pot; printf "9999 \ta\t1\n"; } > $d/p.pot
	par $d/run.par "nchains 1\nTemp 1\nnstep 10\nhb\npivot 1\n"
	( cd $d && "$MG" polyala.pol p.pot run.par > stdout.log 2> stderr.log )
	if [ $? -ne 0 ] && grep -q "hydrogen bond of atoms 9999" $d/stderr.log && sanitizer_ok $d; then pass "hydrogen bond of a missing atom"; else fail "hydrogen bond of a missing atom" "not refused"; fi
	mg err_nohb polyala.pol polyala_hb.pot "nchains 1\nTemp 1\nnstep 10\npivot 1\n"
	if [ $? -eq 0 ] && grep -q "hb is not in the par file" $WORK/err_nohb/stderr.log && sanitizer_ok $WORK/err_nohb; then pass "hydrogen bonds without hb ignored"; else fail "hydrogen bonds without hb ignored" "see $WORK/err_nohb"; fi

	local long=$(printf 'x%.0s' $(seq 1 300))
	( cd $WORK && "$MG" "$long.pol" polyala.pot run.par > /dev/null 2> long.log )
	if [ $? -ne 0 ] && grep -q "too long" $WORK/long.log; then pass "file name too long"; else fail "file name too long" "not refused"; fi

	mkdir -p $WORK/err_box; cp $DATA/three.pol $WORK/err_box/
	printf "[ global ]\nboxtype p\nboxsize 3\n\n" | cat - $DATA/homo.pot > $WORK/err_box/small.pot
	par $WORK/err_box/run.par "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1\nnstep 10\nflip 1\n"
	( cd $WORK/err_box && "$MG" three.pol small.pot run.par > /dev/null 2> stderr.log )
	if [ $? -ne 0 ] && grep -q "periodic box must be larger" $WORK/err_box/stderr.log; then pass "periodic box too small"; else fail "periodic box too small" "not refused"; fi

	mg err_shell three.pol homo_pbc.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1\nnstep 10\nflip 1\nshell\n"
	if [ $? -ne 0 ] && grep -q "shells cannot be used" $WORK/err_shell/stderr.log; then pass "shells with periodic box"; else fail "shells with periodic box" "not refused"; fi

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
three_pbc|three.pol|homo_pbc.pot|nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.5\n
two_pbc|two.pol|two_hf_pbc.pot|nchains 2\nTemp 2.0\n
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
	# 40 chains (the limit used to be 25), whole-chain and cluster moves only
	for pot in homo.pot homo_pbc70.pot; do
		mg ec_many_$pot many.pol $pot "nchains 40\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\nseed 3\nnstep 3000\necheck\nflip 1\npivot 5\nmovecom 5\nmoverot 5\ncomcluster 3\nrotcluster 3\nnprinttrj 100000\nnprinte 100000\nnprintlog 100000\n"
		if [ $? -eq 0 ] && ! grep -q ECHECK $WORK/ec_many_$pot/stderr.log && sanitizer_ok $WORK/ec_many_$pot; then echo PASS; else echo "FAIL many: cluster moves"; fi > $WORK/ec_result_many_40chains_${pot%.pot}
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
	mg rst_$name $pol $pot "seed 5\nnstep 4000\nncheckpoint 3000\nnprinte 200\nnprintobs 200\nnprinttrj 1000\nnprintlog 100000\n$lines"
	local a=$WORK/rst_$name b=$WORK/rst_${name}_b
	rm -rf $b; cp -r $a $b
	awk '$1<3000' $a/energy.ene > $b/energy.ene
	awk '/^#/ || $1<3000' $a/obs.dat > $b/obs.dat
	awk '{print} /^ENDMDL/{n++; if(n==3) exit}' $a/traj.pdb > $b/traj.pdb
	rm -f $b/last.pol
	( cd $b && "$MG" checkpoint $pot run.par > stdout.log 2> stderr.log )
	local res=ok f
	for f in energy.ene traj.pdb last.pol obs.dat; do cmp -s $a/$f $b/$f || res="$f differs"; done
	sanitizer_ok $b || res="sanitizer report"
	if [ "$res" = ok ]; then pass "restart $name"; else fail "restart $name" "$res"; fi
}

group_restart() {
	echo "restart"
	restart_case polyala polyala.pol polyala.pot "nchains 1\nTemp 1.0\nmpivot 1\npivot 4\nsidechain 2\nflip 3\nmovebias 5\nnmul_local 5\n"
	restart_case three three.pol homo.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\nflip 1\npivot 5\ncomcluster 5\nrotcluster 5\nmovecom 7\nmoverot 7\n"
	restart_case two two.pol two_hf.pot "nchains 2\nTemp 2.0\nflip 1\npivot 5\nmovecom 5\nmoverot 5\nsidechain 3\nbackrub 5\nmflip 7\nlpivot 7\nmovebias 4\nnmul_local 5\n"
	restart_case three_pbc three.pol homo_pbc.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\nflip 1\npivot 5\ncomcluster 5\nrotcluster 5\nmovecom 7\nmoverot 7\n"
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
	cells_case many_pbc many.pol homo_pbc70.pot "nchains 40\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\nnstep 5000\nflip 1\npivot 5\nmovecom 5\nmoverot 5\ncomcluster 3\nrotcluster 3\n"
	cells_case two_pbc two.pol two_hf_pbc.pot "nchains 2\nTemp 2.0\nnstep 5000\nflip 1\npivot 5\nmovecom 5\nmoverot 5\nsidechain 3\nbackrub 5\nmflip 7\nmovebias 4\nnmul_local 5\n"
}

######################################################################
group_tools() {
	echo "tools"
	local d=$WORK/grappino
	if [ -x $GRAPPINO ]; then
		mkdir -p $d; cp $ROOT/tutorial/tutorial_1/1YPC.pdb $ROOT/tutorial/tutorial_1/1YPC-CA.in $d/
		( cd $d && $GRAPPINO 1YPC-CA.in > g.log 2>&1 )
		if cmp -s $d/1YPC-CA.pol $DATA/ca.pol && cmp -s $d/1YPC-CA.pot $DATA/ca.pot; then pass "grappino (tutorial 1)"; else fail "grappino (tutorial 1)" "output differs from tests/data/ca.*"; fi
	else
		echo "  SKIP  grappino not built"
	fi
	if [ -x $GRAPPINO ]; then
		# Go model of identical chains: the 16 chains of the TTR(105-115) protofilament share their atom types
		d=$WORK/grappino_homo; mkdir -p $d; cp $DATA/ttr_protofilament.pdb $d/
		printf "pdbfile ttr_protofilament.pdb\nbackbone_atoms 1 CA\nnosidechain\natomtypes go_homo\nimin 4\npotential go\nr_hardcore 4.0\nr_native 7.5\ngo-energy -1.0\ndebug 1\npolfile ttr.pol\npotfile ttr.pot\n" > $d/ttr.in
		( cd $d && $GRAPPINO ttr.in > g.log 2>&1 )
		if grep -q "Found 32 Go contacts" $d/g.log && [ "$(awk '$1~/^[0-9]+$/ && NF==11 {print $4}' $d/ttr.pol | sort -u | wc -l)" -eq 11 ] && [ "$(awk '$1~/^[0-9]+$/ && NF==11' $d/ttr.pol | wc -l)" -eq 176 ]
		then pass "grappino go_homo (TTR protofilament)"; else fail "grappino go_homo (TTR protofilament)" "wrong contacts or types"; fi

		# all-atom model with rotamers of 16 chains numbered alike, then sidechain moves on it
		d=$WORK/grappino_aa; mkdir -p $d; cp $DATA/ttr_protofilament.pdb $d/
		printf "pdbfile ttr_protofilament.pdb\nrotamers\nrotfile $ROOT/lib/rotamers.lib\npdb_rot\natomtypes go_homo\nbackbone_atoms 3 N CA C\nlocked_atoms N\nimin 4\nr_hardcore 2.0\nr_native 4.5\npotential go\ngo-energy -0.2\ngo_dihedrals\ne_dih1 1.2\ne_dih3 0.4\ndebug 1\npolfile ttr.pol\npotfile ttr.pot\n" > $d/ttr.in
		( cd $d && $GRAPPINO ttr.in > g.log 2>&1 )
		par $d/run.par "nchains 16\nseed 1\nnstep 300\nTemp 0.5\nsidechain 1\nflip 2\necheck\nnprinttrj 1000\nnprinte 1000\nnprintlog 1000\n"
		( cd $d && "$MG" ttr.pol ttr.pot run.par > o.log 2> e.log )
		if [ $? -eq 0 ] && grep -q "Found 358 Go contacts" $d/g.log && ! grep -q ECHECK $d/e.log && ! grep -q "WARNING: overlap" $d/e.log
		then pass "grappino all-atom rotamers, 16 chains"; else fail "grappino all-atom rotamers, 16 chains" "see $d"; fi

		# rotamers of a protein with glycines (not in the library): with the pdb rotamer and the CB of the pdb
		# (default), with library rotamers only and the pseudo-CB (cb_pdb 0), and sidechains built from scratch;
		# montegrappa runs on all of them with consistent energies
		for kind in pdb_rot lib scratch; do
			d=$WORK/grappino_rot_$kind; mkdir -p $d; cp $ROOT/tutorial/tutorial_3/hairpin.pdb $d/
			case $kind in pdb_rot) extra="pdb_rot\n" ;; lib) extra="cb_pdb 0\n" ;; scratch) extra="sidescratch\n" ;; esac
			printf "pdbfile hairpin.pdb\nrotamers\nrotfile $ROOT/lib/rotamers.lib\n${extra}atomtypes $ROOT/lib/atomtypes.1.lib\nbackbone_atoms 3 N CA C\nimin 3\npotential go\nr_hardcore 1.0\nr_native 3.0\ngo-energy -1.5\npolfile h.pol\npotfile h.pot\n" > $d/h.in
			( cd $d && $GRAPPINO h.in > g.log 2>&1 ) && \
			par $d/run.par "nchains 1\nseed 3\nnstep 3000\nTemp 1.0\nflip 1\npivot 1\nsidechain 2\nbackrub 1\nmovebias 1\nnmul_local 5\necheck\nnprinttrj 100000\nnprinte 1000\nnprintlog 100000\n" && \
			( cd $d && "$MG" h.pol h.pot run.par > o.log 2> e.log ) && grep -q "^Final" $d/e.log && ! grep -q ECHECK $d/e.log && \
			if [ $kind = lib ]; then grep -q "pseudo-CB" $d/g.log; else ! grep -q "pseudo-CB" $d/g.log; fi
			if [ $? -eq 0 ]; then pass "grappino rotamers with glycines ($kind)"; else fail "grappino rotamers with glycines ($kind)" "see $d"; fi
		done

		# backbone hydrogen bonds (amide H built on N) with the reduced types and a pair potential from a
		# file, statistical dihedrals from tables (phi on CA but not of proline, psi on C: 21 per chain);
		# montegrappa with hb and the dihedral tables on it is consistent
		d=$WORK/grappino_hb; mkdir -p $d; cp $DATA/ttr_protofilament.pdb $DATA/pairs_small.dat $DATA/dihedrals_small.dat $d/
		printf "pdbfile ttr_protofilament.pdb\nrotamers\nrotfile $ROOT/lib/rotamers.lib\npdb_rot\natomtypes $ROOT/lib/atomtypes.reduced.lib\nbackbone_atoms 3 N CA C\nlocked_atoms N\nimin 6\nr_hardcore 1.5\npotential pairs\npairsfile pairs_small.dat\nhb_backbone\ndihfile dihedrals_small.dat\npolfile h.pol\npotfile h.pot\n" > $d/h.in
		( cd $d && $GRAPPINO h.in > g.log 2>&1 )
		par $d/run.par "nchains 16\nseed 2\nnstep 300\nTemp 0.5\nhb\nflip 2\npivot 1\nsidechain 1\nbackrub 2\necheck\nnprinttrj 1000\nnprinte 1000\nnprintlog 1000\n"
		( cd $d && "$MG" h.pol h.pot run.par > o.log 2> e.log )
		if [ $? -eq 0 ] && grep -q "Added 144 amide hydrogens" $d/g.log && grep -q "144 donors, 176 acceptors" $d/g.log && \
		   [ "$(awk '/hydrogen_bonds/{f=1;next} /^\[/{f=0} f && NF==3' $d/h.pot | wc -l)" -eq 320 ] && ! grep -q ECHECK $d/e.log && ! grep -q "WARNING: overlap" $d/e.log && \
		   grep -q "Read 13 bins of dihedral tables" $d/g.log && [ "$(awk '/dihedral_weight/{f=1;next} /^\[/{f=0} f && NF==4' $d/h.pot | wc -l)" -eq 336 ] && \
		   awk 'NR==1 && $5>1 {ok=1} END{exit !ok}' $d/energy.ene
		then pass "grappino backbone hydrogen bonds, reduced types, pairs and dihedrals from files"; else fail "grappino backbone hydrogen bonds, reduced types, pairs and dihedrals from files" "see $d"; fi

		# splice written as given; CB of a CACB model of many chains taken from its own chain;
		# only the first model of a multi-model pdb
		d=$WORK/grappino_misc; mkdir -p $d
		awk '/^ATOM/' $DATA/ttr_protofilament.pdb > $d/one.pdb
		{ echo "MODEL        1"; cat $d/one.pdb; echo "ENDMDL"; echo "MODEL        2"; cat $d/one.pdb; echo "ENDMDL"; } > $d/two_models.pdb
		printf "pdbfile two_models.pdb\nmodel CACB\nbackbone_atoms 2 CA CB\nnosidechain\natomtypes go\nimin 4\npotential go\nr_hardcore 3.0\nr_native 6.5\ngo-energy -1.0\nsplice 0.7 0.5\npolfile m.pol\npotfile m.pot\n" > $d/m.in
		( cd $d && $GRAPPINO m.in > g.log 2>&1 )
		cb_pdb=$(awk '/^ATOM/ && substr($0,13,4)==" CB " && substr($0,22,1)=="P" && substr($0,23,4)+0==2 {printf "%.3f %.3f %.3f", substr($0,31,8), substr($0,39,8), substr($0,47,8)}' $d/one.pdb)
		cb_pol=$(awk '$1~/^[0-9]+$/ && NF==11 && $3=="CB" && $6==2 && $7==15 {printf "%.3f %.3f %.3f", $8, $9, $10}' $d/m.pol)
		if grep -q "in 16 chains" $d/g.log && grep -q "^splice 0.700000 0.500000" $d/m.pot && [ -n "$cb_pdb" ] && [ "$cb_pdb" = "$cb_pol" ]
		then pass "grappino splice, CACB of many chains, multi-model pdb"; else fail "grappino splice, CACB of many chains, multi-model pdb" "see $d (CB pdb '$cb_pdb', pol '$cb_pol')"; fi
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

	# three interacting rigid chains in a periodic box: whole-chain moves with and without avbmc
	for setup in ref avb; do
		for seed in 1 2 3 4; do
			mv="movecom 1\ndx_com 8\nmoverot 1\ndtheta 3\n"; [ $setup = avb ] && mv="${mv}avbmc 1\navb_rin 8\n"
			( mg stat_avb_${setup}_$seed three.pol homo_pbc.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.5\nseed $seed\nnstep 1000000\n${mv}nprinttrj 100000000\nnprinte 100\nnprintlog 100000000\n" ) &
		done
	done
	wait
	for setup in ref avb; do
		for seed in 1 2 3 4; do
			awk '{n++; e[n]=$2} END{s=0; for(i=int(n/5)+1;i<=n;i++) s+=e[i]; print s/(n-int(n/5))}' $WORK/stat_avb_${setup}_$seed/energy.ene
		done | awk -v s=$setup '{m+=$1; m2+=$1*$1; n++} END{mu=m/n; printf "%s %f %f\n", s, mu, sqrt((m2/n-mu*mu)/(n-1))}'
	done > $WORK/stat_avb_summary
	sed 's/^/        /' $WORK/stat_avb_summary
	if awk 'NR==1{m1=$2;s1=$3} NR==2{m2=$2;s2=$3} END{d=m1-m2; if(d<0)d=-d; exit !(d < 3*sqrt(s1*s1+s2*s2))}' $WORK/stat_avb_summary
	then pass "avbmc vs whole-chain moves mean energy"; else fail "avbmc vs whole-chain moves mean energy" "difference above 3 sigma"; fi
}

######################################################################
# fraction of pairs of chains whose centres (minimum image in a box of side L)
# are closer than R, over the frames of a pdb trajectory
PAIRS_IN_AWK='
function flush(   a,b,dx,dy,dz) { for (a=0;a<nc;a++) for (b=a+1;b<nc;b++) {
	dx=x[a]/n[a]-x[b]/n[b]; dy=y[a]/n[a]-y[b]/n[b]; dz=z[a]/n[a]-z[b]/n[b];
	dx-=L*int(dx/L+(dx>0?0.5:-0.5)); dy-=L*int(dy/L+(dy>0?0.5:-0.5)); dz-=L*int(dz/L+(dz>0?0.5:-0.5));
	np++; if (dx*dx+dy*dy+dz*dz<R*R) nin++ }
	split("",x); split("",y); split("",z); split("",n); nc=0 }
/^TITLE/ { flush(); next }
/^ATOM/ { i=index("ABCDEFGHIJKLMNOPQRSTUVWXYZ",substr($0,22,1))-1; if (i+1>nc) nc=i+1;
	x[i]+=substr($0,31,8); y[i]+=substr($0,39,8); z[i]+=substr($0,47,8); n[i]++ }
END { flush(); print nin/np }'

group_avbmc() {
	echo "avbmc (aggregation-volume-bias move)"
	local sys c seed

	mg avb_err_box three.pol homo.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1\nnstep 10\navbmc 1\navb_rin 5\n"
	if [ $? -ne 0 ] && grep -q "periodic box" $WORK/avb_err_box/stderr.log; then pass "avbmc without periodic box refused"; else fail "avbmc without periodic box refused" "not refused"; fi
	mg avb_err_rin three.pol homo_pbc.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1\nnstep 10\navbmc 1\n"
	if [ $? -ne 0 ] && grep -q "avb_rin" $WORK/avb_err_rin/stderr.log; then pass "avbmc without avb_rin refused"; else fail "avbmc without avb_rin refused" "not refused"; fi

	# energy, backup copy and cell list consistent after every move
	for c in -1 1; do
		for sys in "three|three.pol|homo_pbc.pot|nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.5\n" \
		           "two|two.pol|two_hf_pbc.pot|nchains 2\nTemp 2.0\nsidechain 3\n" \
		           "many|many.pol|homo_pbc70.pot|nchains 40\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\n"; do
			IFS="|" read -r name pol pot hdr <<< "$sys"
			mg avb_ec_${name}_$c $pol $pot "${hdr}seed 7\nnstep 3000\ncells $c\necheck\navbmc 1\navb_rin 8\nflip 2\nnprinttrj 100000\nnprinte 100000\nnprintlog 100000\n"
			if [ $? -eq 0 ] && ! grep -q ECHECK $WORK/avb_ec_${name}_$c/stderr.log && sanitizer_ok $WORK/avb_ec_${name}_$c
			then pass "avbmc echeck $name cells=$c"; else fail "avbmc echeck $name cells=$c" "see $WORK/avb_ec_${name}_$c"; fi
		done
	done

	# non-interacting chains: the centres are independent and uniform in the box, so a pair
	# is in the bonded region with probability V_in/V (0.03351 for R=6, L=30); a missing
	# or wrong bias factor gives several times that
	[ $ASAN = 1 ] && { echo "  SKIP  distribution test not with --asan"; return; }
	for seed in 1 2 3 4; do
		( mg avb_free_$seed three.pol free_pbc.pot "nchains 3\nnosidechain\nnodihpot\nnoangpot\nTemp 1.0\nseed $seed\nnstep 1000000\navbmc 1\navb_rin 6\navb_pin 0.3\nnprinttrj 100\nnprinte 1000000\nnprintlog 1000000\n" ) &
	done
	wait
	for seed in 1 2 3 4; do awk -v L=30 -v R=6 "$PAIRS_IN_AWK" $WORK/avb_free_$seed/traj.pdb; done | \
		awk '{m+=$1; m2+=$1*$1; n++} END{mu=m/n; printf "%f %f\n", mu, sqrt((m2/n-mu*mu)/(n-1))}' > $WORK/avb_free_summary
	sed 's/^/        pairs in the bonded region (exact 0.033510): /' $WORK/avb_free_summary
	if awk '{d=$1-0.033510; if(d<0)d=-d; exit !(d < 4*$2+0.002)}' $WORK/avb_free_summary
	then pass "avbmc distribution of non-interacting chains"; else fail "avbmc distribution of non-interacting chains" "$(cat $WORK/avb_free_summary)"; fi
}

######################################################################
# Three sheets of 4 chains of known structure (parallel in register, antiparallel
# in register, parallel shifted by two residues): the first line of obs.dat.
group_observables() {
	echo "observables"
	local f expect
	for f in parallel antiparallel shifted; do
		mg obs_$f obs_$f.pol obs.pot "nchains 4\nseed 1\nnstep 1\nTemp 1.0\nnosidechain\nnodihpot\nnoangpot\nflip 1\ndw_flip 0.00001\nnprintobs 1\nnprinttrj 10\nnprinte 10\nnprintlog 10\n"
		case $f in
			parallel)     expect='$2==1 && $3==4 && $5==3 && $6==3 && $7==0 && $8==3 && $10*$10<0.01 && $12>0.99' ;;
			antiparallel) expect='$2==1 && $3==4 && $5==3 && $6==0 && $7==3 && $9==3 && $11*$11<0.01 && $12>0.99' ;;
			shifted)      expect='$2==1 && $3==4 && $5==3 && $6==3 && $8==0 && $10>-2.1 && $10<-1.9 && $12>0.99' ;;
		esac
		if awk "NR==2 { exit !($expect) } END { if (NR<2) exit 1 }" $WORK/obs_$f/obs.dat; then pass "observables $f sheet"; else fail "observables $f sheet" "$(sed -n 2p $WORK/obs_$f/obs.dat)"; fi
	done
}

######################################################################
# Parallel tempering on 4 ranks, 30-residue Go fragment at T = 0.6 ... 0.9
PT_PAR="nchains 1\nnprinttrj 100000000\nnprintlog 100000000\npivot 1\nmpivot 1\nntemp 4\ntemperatures\n0.6\n0.7\n0.8\n0.9\nstep_exchange 1000\n"

# mpirun with a timeout (a deadlock must fail the test, not hang it)
mpirun4() {
	timeout 900 mpirun --oversubscribe -np 4 "$@"
}

group_mpi_multi() {
	echo "mpi_multi (parallel tempering, 4 ranks)"
	local MPI=$ROOT/bin/montegrappa_mpi d seed ok
	[ $ASAN = 1 ] && { echo "  SKIP  not with --asan"; return; }
	[ -x $MPI ] || { echo "  SKIP  montegrappa_mpi not built"; return; }
	if ! command -v mpirun > /dev/null || ! timeout 60 mpirun --oversubscribe -np 2 true > /dev/null 2>&1; then echo "  SKIP  mpirun does not work here"; return; fi

	# run: all ranks complete, consistent energies, exchanges accepted
	d=$WORK/pt_run; mkdir -p $d; cp $DATA/fragment.pol $DATA/fragment.pot $d/
	par $d/run.par "${PT_PAR}seed 7\nnstep 200000\nnprinte 1000\n"
	( cd $d && mpirun4 $MPI fragment.pol fragment.pot run.par > stdout.log 2> stderr.log )
	ok=1
	for r in 0 1 2 3; do
		[ -f $d/proc_$r.log ] && awk '/^Final:/ { split($2,a,"="); split($3,b,"="); x=a[2]-b[2]; f=1; exit (x>1e-6||x<-1e-6) } END { if (!f) exit 1 }' $d/proc_$r.log || ok=0
	done
	if [ $ok = 1 ]; then pass "pt run (4 ranks)"; else fail "pt run (4 ranks)" "did not complete or energies inconsistent (see $d)"; fi
	if grep -h "Exchanges with rank" $d/stderr.log | awk '{split($0,a,"= "); if (a[2]+0 <= 0) bad=1; n++} END { exit (n!=3 || bad) }'
	then pass "pt exchanges accepted"; else fail "pt exchanges accepted" "$(grep -h 'Exchanges' $d/stderr.log | tr '\n' ' ')"; fi

	# exact restart on 4 ranks
	d=$WORK/pt_rst; mkdir -p $d; cp $DATA/fragment.pol $DATA/fragment.pot $d/
	par $d/run.par "${PT_PAR}seed 5\nnstep 4000\nncheckpoint 3000\nnprinte 200\n"
	( cd $d && mpirun4 $MPI fragment.pol fragment.pot run.par > stdout.log 2> stderr.log )
	rm -rf $d.b; cp -r $d $d.b
	for r in 0 1 2 3; do awk '$1<3000' $d/energy_proc$r.ene > $d.b/energy_proc$r.ene; rm -f $d.b/last_proc$r.pol; done
	( cd $d.b && mpirun4 $MPI checkpoint fragment.pot run.par > stdout.log 2> stderr.log )
	ok=1
	for r in 0 1 2 3; do cmp -s $d/energy_proc$r.ene $d.b/energy_proc$r.ene && cmp -s $d/last_proc$r.pol $d.b/last_proc$r.pol || ok=0; done
	if [ $ok = 1 ]; then pass "pt restart (4 ranks)"; else fail "pt restart (4 ranks)" "files differ (see $d.b)"; fi

	# the lowest-temperature replica samples the canonical distribution at T=0.6
	for seed in 11 22 33 44; do
		d=$WORK/pt_stat_$seed; mkdir -p $d; cp $DATA/fragment.pol $DATA/fragment.pot $d/
		par $d/run.par "${PT_PAR}seed $seed\nnstep 2000000\nnprinte 100\n"
		( cd $d && mpirun4 $MPI fragment.pol fragment.pot run.par > stdout.log 2> stderr.log )
		mg pt_ref_$seed fragment.pol fragment.pot "nchains 1\nseed $seed\nnstep 2000000\nTemp 0.6\npivot 1\nmpivot 1\nnprinttrj 100000000\nnprinte 100\nnprintlog 100000000\n" &
	done
	wait
	for kind in pt ref; do
		for seed in 11 22 33 44; do
			if [ $kind = pt ]; then f=$WORK/pt_stat_$seed/energy_proc0.ene; else f=$WORK/pt_ref_$seed/energy.ene; fi
			awk '{n++; e[n]=$2} END{s=0; for(i=int(n/2)+1;i<=n;i++) s+=e[i]; print s/(n-int(n/2))}' $f
		done | awk -v k=$kind '{m+=$1; m2+=$1*$1; n++} END{mu=m/n; printf "%s %f %f\n", k, mu, sqrt((m2/n-mu*mu)/(n-1))}'
	done > $WORK/pt_summary
	sed 's/^/        /' $WORK/pt_summary
	if awk 'NR==1{m1=$2;s1=$3} NR==2{m2=$2;s2=$3} END{d=m1-m2; if(d<0)d=-d; exit !(NR==2 && d < 3*sqrt(s1*s1+s2*s2))}' $WORK/pt_summary
	then pass "pt lowest replica vs plain run"; else fail "pt lowest replica vs plain run" "difference above 3 sigma"; fi
}

######################################################################
echo "MonteGrappa tests with $MG"
echo "working directory $WORK"
for g in $TGROUPS; do
	case $g in
		errors|regression|echeck|restart|cells|tools|mpi|observables|avbmc|mpi_multi|stat) group_$g ;;
		*) echo "unknown group $g"; NFAIL=$((NFAIL+1)) ;;
	esac
done

[ $KEEP = 1 ] || [ $NFAIL -gt 0 ] || rm -rf $WORK
echo
echo "passed $NPASS, failed $NFAIL"
[ $NFAIL -gt 0 ] && { printf "failed:$FAILED\n"; echo "(outputs kept in $WORK)"; }
exit $NFAIL
