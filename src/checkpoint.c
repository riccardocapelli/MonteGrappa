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
 * checkpoint.c
 *
 * Exact restart: the checkpoint file contains the state of the MC loop at the
 * end of a step (step, run, counters, state of the random number generator),
 * the matrix of the potential when the optimization changes it, the state of
 * simulated tempering and the coordinates of all backbone and sidechain atoms,
 * written with 17 significant digits, so that they are read back bit by bit.
 * The file is written to a temporary name and then renamed, so an interrupted
 * write never destroys the previous checkpoint. The frames recorded for the
 * optimization of the potential, which can be many, are appended to a
 * separate file at each checkpoint (OP_SaveFrames), whose length is in the
 * checkpoint as that of the other output files.
 */

#include "montegrappa.h"
#include <errno.h>
#include <unistd.h>

#define CHECKPOINT_VERSION 4

/*****************************************************************************
 Name of the checkpoint file of a given rank
 *****************************************************************************/
void CheckpointName(char *fname, size_t size, int my_rank)
{
	#ifdef ACTIVE_MPI
	snprintf(fname,size,"checkpoint_proc%d.dat",my_rank);
	#else
	snprintf(fname,size,"checkpoint.dat");
	(void) my_rank;
	#endif
}

/*****************************************************************************
 Name of the file of the frames of the optimization of a given rank
 *****************************************************************************/
void CheckpointOPName(char *fname, size_t size, int my_rank)
{
	#ifdef ACTIVE_MPI
	snprintf(fname,size,"checkpoint_op_proc%d.dat",my_rank);
	#else
	snprintf(fname,size,"checkpoint_op.dat");
	(void) my_rank;
	#endif
}

/*****************************************************************************
 Write checkpoint c, the matrix of the potential if it is optimized, the
 state of simulated tempering and the coordinates of the chains of p
 *****************************************************************************/
void WriteCheckpoint(char *fname, struct s_checkpoint *c, struct s_polymer *p, struct s_mc_parms *parms, struct s_potential *pot)
{
	char tmpname[FNAMEMAX+10];
	FILE *fp;
	int ic,i,j,npol=parms->npol;
	struct s_back *b;

	snprintf(tmpname,sizeof(tmpname),"%s.tmp",fname);
	fp = fopen(tmpname,"w");
	if (!fp) Error("Cannot open checkpoint file for writing");

	fprintf(fp,"montegrappa_checkpoint %d\n",CHECKPOINT_VERSION);
	fprintf(fp,"polname %s\n",c->polname);
	fprintf(fp,"istep %llu\n",c->istep);
	fprintf(fp,"irun %d\n",c->irun);
	fprintf(fp,"rng %llu %llu %llu %llu\n",c->rng[0],c->rng[1],c->rng[2],c->rng[3]);
	fprintf(fp,"counters %d %d %d %d %d %lld %lld\n",c->iprinttrj,c->iprintlog,c->iprinte,c->icheckpoint,c->ptcount,c->acc,c->mov);
	fprintf(fp,"moves %d\n",NMOVES);
	for (i=0;i<NMOVES;++i) fprintf(fp,"%d %lld %lld\n",c->mcount[i],c->macc[i],c->mdone[i]);
	fprintf(fp,"files %d",NCKPFILES);
	for (i=0;i<NCKPFILES;++i) fprintf(fp," %lld",c->off[i]);
	fprintf(fp,"\nstate %.17g %.17g %d %d %d %d %d %d\n",c->t,c->etot,c->st_iprint,c->st_o,c->nconf,c->op_icount,c->exc,c->exa);

	// the optimization changes the matrix of the potential at the end of each run
	if (strcmp(parms->op_minim,"none"))
	{
		fprintf(fp,"potential %d\n",c->ntypes);
		for (i=0;i<c->ntypes;++i)
			for (j=0;j<c->ntypes;++j)
				fprintf(fp,"%.17g%c",pot->e[i][j],j==c->ntypes-1 ? '\n' : ' ');
	}

	#ifdef STEMPERING
	if (parms->stempering)
	{
		fprintf(fp,"stempering\n");
		STWriteState(fp,parms->p);
	}
	#endif

	fprintf(fp,"chains %d\n",npol);
	for (ic=0;ic<npol;++ic)
	{
		fprintf(fp,"chain %d %d\n",ic,(p+ic)->nback);
		for (i=0;i<(p+ic)->nback;++i)
		{
			b = ((p+ic)->back)+i;
			fprintf(fp,"%d %d %.17g %.17g %.17g\n",b->irot,b->nside,b->pos.x,b->pos.y,b->pos.z);
			for (j=0;j<b->nside;++j)
				fprintf(fp,"%.17g %.17g %.17g\n",((b->side)+j)->pos.x,((b->side)+j)->pos.y,((b->side)+j)->pos.z);
		}
	}

	if (fflush(fp) || ferror(fp)) Error("Cannot write checkpoint file");
	fclose(fp);
	if (rename(tmpname,fname)) Error("Cannot rename checkpoint file");
}

/*****************************************************************************
 Read the state part of a checkpoint (everything but the coordinates)
 *****************************************************************************/
void ReadCheckpointHeader(char *fname, struct s_checkpoint *c)
{
	FILE *fp;
	int version,nmoves,i;

	fp = fopen(fname,"r");
	if (!fp) Error("Cannot open checkpoint file (restart with 'checkpoint' as polymer file)");

	if (fscanf(fp," montegrappa_checkpoint %d",&version)!=1 || version!=CHECKPOINT_VERSION)
		Error("Checkpoint file not written by this version of MonteGrappa");
	if (fscanf(fp," polname %255s",c->polname)!=1) Error("Cannot read polname in checkpoint");
	if (fscanf(fp," istep %llu",&(c->istep))!=1) Error("Cannot read istep in checkpoint");
	if (fscanf(fp," irun %d",&(c->irun))!=1) Error("Cannot read irun in checkpoint");
	if (fscanf(fp," rng %llu %llu %llu %llu",&(c->rng[0]),&(c->rng[1]),&(c->rng[2]),&(c->rng[3]))!=4)
		Error("Cannot read random number generator state in checkpoint");
	if (fscanf(fp," counters %d %d %d %d %d %lld %lld",&(c->iprinttrj),&(c->iprintlog),&(c->iprinte),
			&(c->icheckpoint),&(c->ptcount),&(c->acc),&(c->mov))!=7)
		Error("Cannot read counters in checkpoint");
	if (fscanf(fp," moves %d",&nmoves)!=1 || nmoves!=NMOVES) Error("Wrong number of move types in checkpoint");
	for (i=0;i<NMOVES;++i)
		if (fscanf(fp,"%d %lld %lld",&(c->mcount[i]),&(c->macc[i]),&(c->mdone[i]))!=3)
			Error("Cannot read move counters in checkpoint");
	if (fscanf(fp," files %d",&i)!=1 || i!=NCKPFILES) Error("Cannot read lengths of the output files in checkpoint");
	for (i=0;i<NCKPFILES;++i)
		if (fscanf(fp,"%lld",&(c->off[i]))!=1) Error("Cannot read lengths of the output files in checkpoint");
	if (fscanf(fp," state %lf %lf %d %d %d %d %d %d",&(c->t),&(c->etot),&(c->st_iprint),&(c->st_o),&(c->nconf),&(c->op_icount),&(c->exc),&(c->exa))!=8)
		Error("Cannot read state of the MC loop in checkpoint");
	fclose(fp);
}

/*****************************************************************************
 Open the checkpoint after the line that starts with key, which is copied to
 line (NULL if there is no such line)
 *****************************************************************************/
static FILE *CheckpointSection(char *fname, char *key, char *line, int size)
{
	FILE *fp;
	size_t n=strlen(key);

	fp = fopen(fname,"r");
	if (!fp) Error("Cannot open checkpoint file");
	while (fgets(line,size,fp)!=NULL)
	{
		if (!strncmp(line,key,n) && (line[n]==' ' || line[n]=='\n')) return fp;
		if (!strncmp(line,"chains ",7)) break;			// the sections come before the coordinates
	}
	fclose(fp);
	return NULL;
}

/*****************************************************************************
 Overwrite the matrix of the potential with the one in the checkpoint, if
 there is one (the optimization changed it)
 *****************************************************************************/
void ReadCheckpointPotential(char *fname, struct s_potential *u, int ntypes)
{
	FILE *fp;
	char aux[500];
	int i,j,n;

	fp = CheckpointSection(fname,"potential",aux,500);
	if (!fp) return;
	if (sscanf(aux,"potential %d",&n)!=1 || n!=ntypes) Error("Matrix of the potential in checkpoint differs from the potential file");
	for (i=0;i<ntypes;++i)
		for (j=0;j<ntypes;++j)
			if (fscanf(fp,"%lf",&(u->e[i][j]))!=1) Error("Cannot read matrix of the potential in checkpoint");
	fclose(fp);
}

#ifdef STEMPERING
/*****************************************************************************
 Overwrite the state of simulated tempering with the one in the checkpoint
 *****************************************************************************/
void ReadCheckpointST(char *fname, struct st_stuff *st)
{
	FILE *fp;
	char aux[500];

	fp = CheckpointSection(fname,"stempering",aux,500);
	if (!fp) Error("No state of simulated tempering in checkpoint");
	STReadState(fp,st);
	fclose(fp);
}
#endif

/*****************************************************************************
 Cut an output file to the length it had at the checkpoint (off<0: the file
 was not written), so that what was written after it is not repeated
 *****************************************************************************/
void TruncateOutput(char *fname, long long off)
{
	if (off < 0) return;
	if (truncate(fname,(off_t) off))
	{
		if (errno==ENOENT && off==0) return;
		fprintf(stderr,"Cannot cut %s to %lld bytes: %s\n",fname,off,strerror(errno));
		Error("Cannot restart the output files");
	}
}

/*****************************************************************************
 Overwrite the coordinates and rotamers of p with those in the checkpoint;
 p must already have the topology of the system (read from the polymer file)
 *****************************************************************************/
void ReadCheckpointCoordinates(char *fname, struct s_polymer *p, int npol)
{
	FILE *fp;
	char aux[500];
	int ic,i,j,npolr,icr,nback,nside;
	struct s_back *b;

	fp = fopen(fname,"r");
	if (!fp) Error("Cannot open checkpoint file");

	// skip the header
	while (fgets(aux,500,fp)!=NULL)
		if (!strncmp(aux,"chains ",7)) break;
	if (sscanf(aux,"chains %d",&npolr)!=1 || npolr!=npol) Error("Number of chains in checkpoint differs from nchains");

	for (ic=0;ic<npol;++ic)
	{
		if (fscanf(fp," chain %d %d",&icr,&nback)!=2 || icr!=ic || nback!=(p+ic)->nback)
			Error("Chain in checkpoint differs from polymer file");
		for (i=0;i<nback;++i)
		{
			b = ((p+ic)->back)+i;
			if (fscanf(fp,"%d %d %lf %lf %lf",&(b->irot),&nside,&(b->pos.x),&(b->pos.y),&(b->pos.z))!=5 || nside!=b->nside ||
					(nside>0 && (b->irot<0 || b->irot>=b->nrot)))
				Error("Residue in checkpoint differs from polymer file");
			for (j=0;j<nside;++j)
				if (fscanf(fp,"%lf %lf %lf",&(((b->side)+j)->pos.x),&(((b->side)+j)->pos.y),&(((b->side)+j)->pos.z))!=3)
					Error("Cannot read sidechain in checkpoint");
		}
	}
	fclose(fp);
}
