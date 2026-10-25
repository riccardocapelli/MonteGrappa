/*
 * dcd2pdb: convert a dcd trajectory of MonteGrappa (trjformat dcd) to pdb.
 *
 *    dcd2pdb traj_top.pdb traj.dcd [first [last [stride]]] > traj_frames.pdb
 *
 * The topology (the ATOM lines of traj_top.pdb, written by montegrappa
 * next to the dcd file) gives the names of the atoms; each frame of the dcd
 * file becomes a MODEL, with a CRYST1 record when it has a unit cell.
 * Frames are numbered from 0; last=-1 means up to the end.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINEMAX 200

static void Die(const char *msg)
{
	fprintf(stderr,"dcd2pdb: %s\n",msg);
	exit(1);
}

// read a fortran record of exactly n bytes (n<0: any length, returned)
static int Record(FILE *fp, void *buf, int n, int maxn)
{
	int a,b;

	if (fread(&a,sizeof(int),1,fp)!=1) return -1;
	if ((n>=0 && a!=n) || a<0 || a>maxn) Die("unexpected record length (not a dcd file of montegrappa?)");
	if (fread(buf,1,a,fp)!=(size_t)a || fread(&b,sizeof(int),1,fp)!=1 || a!=b) return -1;
	return a;
}

int main(int argc, char *argv[])
{
	FILE *ft,*fd;
	char (*lines)[LINEMAX],aux[LINEMAX],head[84],*title;
	int nat=0,nmax=1000,natdcd,icntrl[20],cell,ntitle,iframe,first=0,last=-1,stride=1,i;
	float *x;
	double uc[6];

	if (argc<3 || argc>6)
	{
		fprintf(stderr,"usage: dcd2pdb traj_top.pdb traj.dcd [first [last [stride]]] > traj_frames.pdb\n");
		exit(1);
	}
	if (argc>3) first = atoi(argv[3]);
	if (argc>4) last = atoi(argv[4]);
	if (argc>5) stride = atoi(argv[5]);
	if (first<0 || stride<1) Die("first must be >= 0 and stride >= 1");

	// topology
	ft = fopen(argv[1],"r");
	if (!ft) Die("cannot open topology file");
	lines = malloc(nmax*sizeof(*lines));
	if (!lines) Die("cannot allocate");
	while (fgets(aux,LINEMAX,ft))
		if (!strncmp(aux,"ATOM  ",6) || !strncmp(aux,"HETATM",6))
		{
			if (strlen(aux)<54) Die("short ATOM line in topology file");
			if (nat==nmax)
			{
				nmax *= 2;
				lines = realloc(lines,nmax*sizeof(*lines));
				if (!lines) Die("cannot allocate");
			}
			strncpy(lines[nat],aux,LINEMAX-1);
			lines[nat][LINEMAX-1] = '\0';
			nat++;
		}
	fclose(ft);

	// header
	fd = fopen(argv[2],"rb");
	if (!fd) Die("cannot open dcd file");
	if (Record(fd,head,84,84)!=84 || strncmp(head,"CORD",4)) Die("not a dcd file (or of the other endianness)");
	memcpy(icntrl,head+4,sizeof(icntrl));
	cell = icntrl[10];
	if (fread(&i,sizeof(int),1,fd)!=1 || i<4) Die("bad title record");
	title = malloc(i);
	if (!title || fread(title,1,i,fd)!=(size_t)i || fread(&ntitle,sizeof(int),1,fd)!=1) Die("bad title record");
	free(title);
	if (Record(fd,&natdcd,4,4)!=4) Die("bad number of atoms");
	if (natdcd!=nat) { fprintf(stderr,"dcd2pdb: %d atoms in the dcd file, %d in the topology\n",natdcd,nat); exit(1); }

	x = malloc(3*nat*sizeof(float));
	if (!x) Die("cannot allocate");
	for (iframe=0;last<0 || iframe<=last;++iframe)
	{
		if (cell && Record(fd,uc,sizeof(uc),sizeof(uc))<0) break;
		if (Record(fd,x,nat*sizeof(float),nat*sizeof(float))<0 || Record(fd,x+nat,nat*sizeof(float),nat*sizeof(float))<0 ||
			Record(fd,x+2*nat,nat*sizeof(float),nat*sizeof(float))<0) break;
		if (iframe<first || (iframe-first)%stride) continue;
		printf("MODEL %8d\n",iframe);
		if (cell) printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",uc[0],uc[2],uc[5],uc[4],uc[3],uc[1]);
		for (i=0;i<nat;++i)
			printf("%.30s%8.3f%8.3f%8.3f%s",lines[i],x[i],x[nat+i],x[2*nat+i],strlen(lines[i])>54 ? lines[i]+54 : "\n");
		printf("ENDMDL\n");
	}
	fclose(fd);
	free(x); free(lines);
	return 0;
}
