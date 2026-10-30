/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


    Copyright (C) 2014  G. Tiana, F. Villa, Y. Zhan, R. Capelli, C.Paissoni, P. Sormanni, R. Meloni

    MonteGrappa is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


    This program is free. We just ask you to cite in the publications you get
    making use of    MonteGrappa the following article:

    G. Tiana, F. Villa, Y. Zhan, R. Capelli, C.Paissoni, P. Sormanni, E. Heard,
    L. Giorgetti and R. Meloni, "MonteGrappa: an iterative Monte Carlo program
    to optimize biomolecular potentials in simplified models" Comp. Phys. Comm.
    Volume 186, January 2015, Pages 93–104.     DOI: 10.1016/j.cpc.2014.09.012

 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * trajectory.c
 *
 * Trajectory files, in the format chosen with trjformat:
 *  - pdb (default): one model per frame (PrintPDBStream);
 *  - dcd: binary CHARMM/NAMD format (single precision, read by VMD,
 *    MDAnalysis, mdtraj...), with the unit cell in a periodic box. The atoms
 *    are in the order of the pdb frames; their names are in a topology file
 *    <traj>_top.pdb, the initial structure in a single model.
 * A restart appends to the trajectory of the interrupted run, as for the pdb
 * format; a frame of the dcd file cut by the interruption is overwritten.
 */

#include "montegrappa.h"

#define DCD_NTITLE	2
#define DCD_HEADER	(4+84+4 + 4+4+80*DCD_NTITLE+4 + 4+4+4)	// bytes before the first frame
#define DCD_NSET	8					// offset of the number of frames

static void DCDRecord(FILE *fp, const void *data, int n)
{
	if (fwrite(&n,sizeof(int),1,fp)!=1 || fwrite(data,1,n,fp)!=(size_t)n || fwrite(&n,sizeof(int),1,fp)!=1)
		Error("Cannot write dcd trajectory");
}

static int CountAtoms(struct s_polymer *p, int npol)
{
	int ic,i,n=0;

	for (ic=0;ic<npol;++ic)
		for (i=0;i<(p+ic)->nback;++i)
			n += 1 + (((p+ic)->back)+i)->nside;
	return n;
}

static long DCDFrameSize(int natoms, int cell)
{
	return (cell ? 4+6*sizeof(double)+4 : 0) + 3*(4+4L*natoms+4);
}

/********************************************************************
 Topology of the dcd file: the structure in a single model, chains
 separated by TER, atoms in the order of the frames
 ********************************************************************/
static void PrintPDBTopology(char *fname, struct s_polymer *p, int npol)
{
	int ic,i,j,n=0;
	char cc,bx[9],by[9],bz[9];
	FILE *fp;
	struct s_back *b;
	struct s_side *s;

	fp = fopen(fname,"w");
	if (!fp) Error("Cannot open topology file of the dcd trajectory");
	fprintf(fp,"REMARK   atoms in the order of the frames of the dcd trajectory\n");
	for (ic=0;ic<npol;++ic)
	{
		cc = 'A' + ic%26;
		for (i=0;i<(p+ic)->nback;++i)
		{
			b = ((p+ic)->back)+i;
			fprintf(fp,"ATOM  %5d  %-4s%3s %c%4d    %s%s%s  1.00  0.00      %4d %c\n",(++n)%100000,b->type,
					b->aa,cc,b->iaa%10000,PDBCoord(b->pos.x,bx),PDBCoord(b->pos.y,by),PDBCoord(b->pos.z,bz),ic%10000,b->type[0]);
			for (j=0;j<b->nside;++j)
			{
				s = (b->side)+j;
				fprintf(fp,"ATOM  %5d  %-4s%3s %c%4d    %s%s%s  1.00  0.00      %4d %c\n",(++n)%100000,s->type,
						b->aa,cc,b->iaa%10000,PDBCoord(s->pos.x,bx),PDBCoord(s->pos.y,by),PDBCoord(s->pos.z,bz),ic%10000,s->type[0]);
			}
		}
		fprintf(fp,"TER\n");
	}
	fprintf(fp,"END\n");
	fclose(fp);
}

/********************************************************************
 Open the trajectory file <base>.pdb or <base>.dcd; a new dcd file gets
 its header and the topology <base>_top.pdb
 ********************************************************************/
FILE *OpenTrajectory(char *base, struct s_polymer *p, struct s_mc_parms *parms, struct s_potential *u, int restart)
{
	FILE *fp;
	char fname[FNAMEMAX+50],title[80*DCD_NTITLE+1];
	int icntrl[20],n,k,cell;
	long size,fsize;

	parms->trj_cell = (u->boxtype=='p');
	parms->trj_natoms = CountAtoms(p,parms->npol);
	parms->trj_nframes = 0;

	if (!parms->trjdcd)
	{
		snprintf(fname,sizeof(fname),"%s.pdb",base);
		return fopen(fname,restart ? "a" : "w");
	}

	snprintf(fname,sizeof(fname),"%s.dcd",base);
	if (restart && (fp=fopen(fname,"r+b")))
	{
		// check the header and go after the last complete frame
		if (fread(&n,sizeof(int),1,fp)!=1 || n!=84 || fread(title,1,4,fp)!=4 || strncmp(title,"CORD",4) ||
			fseek(fp,DCD_HEADER-8,SEEK_SET) || fread(&n,sizeof(int),1,fp)!=1 || n!=parms->trj_natoms ||
			fseek(fp,4+4+4*10,SEEK_SET) || fread(&cell,sizeof(int),1,fp)!=1 || (cell!=0)!=parms->trj_cell)
			Error("The dcd trajectory to continue was not written for this system");
		fseek(fp,0,SEEK_END);
		size = ftell(fp);
		fsize = DCDFrameSize(parms->trj_natoms,parms->trj_cell);
		parms->trj_nframes = (int) ((size-DCD_HEADER)/fsize);
		fseek(fp,DCD_HEADER+parms->trj_nframes*fsize,SEEK_SET);
		return fp;
	}

	fp = fopen(fname,"w+b");
	if (!fp) return NULL;

	memset(icntrl,0,sizeof(icntrl));
	icntrl[0] = 0;					// frames, updated at each frame
	icntrl[1] = 0;					// first step
	icntrl[2] = parms->nprinttrj;			// steps between frames
	{ float delta=1.; memcpy(icntrl+9,&delta,sizeof(float)); }
	icntrl[10] = parms->trj_cell;			// unit cell in each frame
	icntrl[19] = 24;				// CHARMM version
	n = 84;
	if (fwrite(&n,sizeof(int),1,fp)!=1 || fwrite("CORD",1,4,fp)!=4 || fwrite(icntrl,sizeof(int),20,fp)!=20 || fwrite(&n,sizeof(int),1,fp)!=1)
		Error("Cannot write dcd trajectory");

	memset(title,' ',sizeof(title));
	n = 4+80*DCD_NTITLE;
	snprintf(title,80,"REMARKS MonteGrappa trajectory");
	snprintf(title+80,80,"REMARKS atoms as in %s_top.pdb",base);
	for (k=0;k<80*DCD_NTITLE;++k) if (title[k]=='\0') title[k]=' ';
	k = DCD_NTITLE;
	if (fwrite(&n,sizeof(int),1,fp)!=1 || fwrite(&k,sizeof(int),1,fp)!=1 || fwrite(title,1,80*DCD_NTITLE,fp)!=80*DCD_NTITLE ||
		fwrite(&n,sizeof(int),1,fp)!=1)
		Error("Cannot write dcd trajectory");
	DCDRecord(fp,&(parms->trj_natoms),sizeof(int));

	snprintf(fname,sizeof(fname),"%s_top.pdb",base);
	PrintPDBTopology(fname,p,parms->npol);

	return fp;
}

/********************************************************************
 Write a frame of the trajectory
 ********************************************************************/
void PrintTrajectory(FILE *fp, struct s_polymer *p, struct s_mc_parms *parms, struct s_potential *u)
{
	int ic,i,j,n=0,nat;
	float *x;
	double cell[6];
	struct s_back *b;

	if (!parms->trjdcd)
	{
		PrintPDBStream(p,parms->npol,fp);
		return;
	}

	nat = parms->trj_natoms;
	if (CountAtoms(p,parms->npol)!=nat) Error("The number of atoms changed during the dcd trajectory");
	x = (float *) malloc(3*nat*sizeof(float));
	if (!x) Error("Cannot allocate frame of the dcd trajectory");
	for (ic=0;ic<parms->npol;++ic)
		for (i=0;i<(p+ic)->nback;++i)
		{
			b = ((p+ic)->back)+i;
			x[n] = (b->pos).x; x[nat+n] = (b->pos).y; x[2*nat+n] = (b->pos).z; n++;
			for (j=0;j<b->nside;++j)
			{
				x[n] = ((b->side)+j)->pos.x; x[nat+n] = ((b->side)+j)->pos.y; x[2*nat+n] = ((b->side)+j)->pos.z; n++;
			}
		}

	if (parms->trj_cell)
	{
		// a, gamma, b, beta, alpha, c
		cell[0] = cell[2] = cell[5] = u->lbox;
		cell[1] = cell[3] = cell[4] = 90.;
		DCDRecord(fp,cell,sizeof(cell));
	}
	DCDRecord(fp,x,nat*sizeof(float));
	DCDRecord(fp,x+nat,nat*sizeof(float));
	DCDRecord(fp,x+2*nat,nat*sizeof(float));
	free(x);

	// number of frames in the header
	parms->trj_nframes++;
	fseek(fp,DCD_NSET,SEEK_SET);
	if (fwrite(&(parms->trj_nframes),sizeof(int),1,fp)!=1) Error("Cannot write dcd trajectory");
	fseek(fp,0,SEEK_END);
}
