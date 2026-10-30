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
 * optimizepot.h
 *
 *  Created on: Sep 19, 2011
 *      Author: guido
 */

#include "montegrappa.h"


/***************************************************
 Allocate op structure and read constrain from file
 ***************************************************/
struct s_optimizepot *InitializeOptimizePot(struct s_mc_parms *parms, int ntypes, struct s_potential *u, FILE *fproc, int iproc, FILE *fp, int *nrestr)
{
	int i,j,nframes;
	struct s_optimizepot *x;

	// if it does need to optimize potential, allocate the structure and set record=0
	if (!strcmp(parms->op_minim,"none"))
	{
		x = (struct s_optimizepot *) calloc(1,sizeof(struct s_optimizepot));
		if (x==NULL) Error("Cannot allocate s_optimizepot");
		x->record= 0;
		return x;
	}

	fprintf(fproc,"\nActivating module OPTIMIZE POTENTIAL\n\nInitialize structures\n");
	
	// read restrains from file (and allocates input structure)
	if(iproc==0)
	{
		parms->op_input = ReadOPRestrains(parms);
		fprintf(stderr," ntypes=%d\n OP_NCONTMAX=%d\n nrestr=%d\n nframes=%llu\n",ntypes,OP_NCONTMAX,parms->op_input->ndata,parms->nstep/parms->op_deltat+1);
	
	*nrestr=parms->op_input->ndata;
	}	


	

	// allocate structure for recording conformations
	#ifdef ACTIVE_MPI
	MPI_Bcast(nrestr,1,MPI_INT,0,MPI_COMM_WORLD);
	if(iproc!=0){ 
	parms->op_input = Allo_op_input(*nrestr);
	}
	parms->op_input = send_op_input((*nrestr),parms->op_input);
	fprintf(fproc," ntypes=%d\n OP_NCONTMAX=%d\n nrestr=%d\n nframes=%llu\n",ntypes,OP_NCONTMAX,parms->op_input->ndata,parms->nstep/parms->op_deltat+1);
	
	#endif
	// frames of this rank; with op_allreplicas rank 0 also keeps those of the others
	nframes = (parms->nstep / parms->op_deltat) + 2;
	#ifdef STEMPERING
	if (parms->stempering && (parms->p)->st_printpdb > 0 && (parms->p)->st_printpdb < parms->op_deltat)	// harvested every st_printpdb steps
		nframes = (parms->nstep / (parms->p)->st_printpdb) + 2;
	#endif
	#ifdef ACTIVE_MPI
	if (parms->op_allreplicas && iproc==0) nframes *= parms->ntemp;
	#else
	if (parms->op_allreplicas) Error("op_allreplicas needs parallel tempering (montegrappa_mpi)");
	#endif
	x = AlloOptimizePot(parms,ntypes,parms->op_input->ndata,nframes,fproc);

	#ifdef OP_DEBUG
	 fprintf(fp,"LOG of OPTIMIZEPOTENTIAL\n");
	#endif

	// set default hardcores and width
	if (parms->op_r>0)
		for (i=0;i<ntypes;++i)
			for (j=i;j<ntypes;++j)					// i,j are atom types, not positions along the chain
				if ( u->r_2[i][j] <0 )
	  			{
					 u->r_2[i][j] = parms->op_r * parms->op_r;
					 u->r_2[j][i] = parms->op_r * parms->op_r;
				  	 u->r0_2[i][j] = parms->op_r0 * parms->op_r0;
					 u->r0_2[j][i] = parms->op_r0 * parms->op_r0;
				}

	return x;
}

/***************************************************
 Allocate op structure
 The contacts of each frame are stored sparse; the
 arrays of contacts and entries grow when needed
 ***************************************************/
struct s_optimizepot *AlloOptimizePot(struct s_mc_parms *parms, int ntypes, int nrestr, int nframes, FILE *fproc)
{
	struct s_optimizepot *x;
	int i,j;

	x = (struct s_optimizepot *) calloc(1,sizeof(struct s_optimizepot));
	if (x==NULL) Error("Cannot allocate s_optimizepot");

	x->ntypes = ntypes;
	x->cindex = AlloIntMatrix(ntypes,ntypes);
	for (i=0;i<ntypes;++i)
		for (j=0;j<ntypes;++j) x->cindex[i][j] = -1;
	x->nallocont = 1;		// i.e., how many multiples of OP_NCONTMAX
	x->it1 = AlloInt(OP_NCONTMAX);
	x->it2 = AlloInt(OP_NCONTMAX);
	x->fmul = AlloDouble(OP_NCONTMAX);
	x->ftouched = AlloInt(OP_NCONTMAX);

	x->nnzmax = 16*(long long)OP_NCONTMAX;
	x->cidx = (int *) malloc(x->nnzmax*sizeof(int));
	x->cmul = (double *) malloc(x->nnzmax*sizeof(double));
	x->fstart = (long long *) calloc(nframes+1,sizeof(long long));
	if (!x->cidx || !x->cmul || !x->fstart) Error("Cannot allocate the contacts of the frames");

	x->eold = AlloDouble(nframes);
	x->enew = AlloDouble(nframes);	
	x->efix = AlloDouble(nframes);
	x->t = AlloDouble(nframes);
	x->lnden = AlloDouble(nframes);
	x->restrain = AlloDoubleMatrix(nframes,nrestr);
	x->rest_av = AlloDouble(nrestr);
	x->w = AlloDouble(nframes);
	x->wold = AlloDouble(nframes);
	x->sw = AlloDouble(nrestr);
	x->swold = AlloDouble(nrestr);

	fprintf(fproc,"Allocate OP structures for %d frames, %lf MB (the contacts of the frames grow as needed)\n",nframes,
		((double)ntypes*ntypes*sizeof(int)+7.*nframes*sizeof(double)+nframes*nrestr*sizeof(double))/1024/1024);
	
	x->nframesmax = nframes;
	x->nframes = 0;
	x->icount = 0;
	x->record= 0;
	x->ncontacts = 0;
	x->nnz = 0;
	x->nftouched = 0;

	return x;
}

/***************************************************
 Forget the recorded frames and contacts
 ***************************************************/
void OP_Reset(struct s_optimizepot *x)
{
	int i;

	if (x->cindex)
		for (i=0;i<x->ncontacts;++i)
			x->cindex[x->it1[i]][x->it2[i]] = x->cindex[x->it2[i]][x->it1[i]] = -1;
	x->ncontacts = 0;
	x->nnz = 0;
	x->nframes = 0;
	x->nftouched = 0;
	if (x->fstart) x->fstart[0] = 0;
	if (x->fsave) fclose(x->fsave);			// the next checkpoint starts a new file of frames
	x->fsave = NULL;
	x->nfsaved = x->ncsaved = 0;
	free(x->cstart); free(x->cfr); free(x->cmc);		// the copy by pair of types of the last reweighting
	x->cstart = NULL; x->cfr = NULL; x->cmc = NULL;
}

void FreeOpt(struct s_optimizepot *x,struct s_mc_parms *parms,int nrestr){

	int i;

	if (x->cindex)
	{
		for (i=0;i<x->ntypes;++i) free(x->cindex[i]);
		free(x->cindex);
	}
	free(x->it1);
	free(x->it2);
	free(x->fmul);
	free(x->ftouched);
	if (x->fsave) fclose(x->fsave);
	free(x->fstart);
	free(x->cidx);
	free(x->cmul);
	free(x->eold);
	free(x->enew);
	free(x->efix);
	free(x->t);
	free(x->lnden);
	free(x->cstart);
	free(x->cfr);
	free(x->cmc);
	free(x->rest_av);
	free(x->w);
	free(x->wold);
	free(x->sw);
	free(x->swold);
	if (x->restrain) FreeDoubleMatrix((x->restrain),x->nframesmax);
	free(x);
}

/***************************************************
 Read list of restrains from file.
 Format:  i1  i2  i3  i4  type  value  sigma 
 A negative sigma marks a validation restraint: its
 average is reported (chi2 with |sigma|), but it is
 not fitted
 ***************************************************/

int NFitted(struct s_optimizepot_input *in)
{
	int i,n=0;

	for (i=0;i<in->ndata;++i)
		if (in->sigma[i]>0) n++;
	return n;
}

struct s_optimizepot_input *ReadOPRestrains(struct s_mc_parms *parms)
{
	struct s_optimizepot_input *x;
	int i;
	char aux[500];
	FILE *fp;
        	
	// allocate input structure
	x = (struct s_optimizepot_input *) calloc(1,sizeof(struct s_optimizepot_input));
	if (x==NULL) Error("Cannot allocate s_optimizepot_input");

	// open restrain file
	fp = fopen(parms->fnop,"r");
	if (!fp) Error("Cannot open restrain file");	
	
	//read first lines
	//fscanf(fp,"%s",aux);
	fgets(aux,500,fp);
	if (sscanf(aux,"ndata %d",&(x->ndata)) != 1) Error("First line of restrain data should be: ndata <int>");

	// allocate arrays
	x->expdata = AlloDouble(x->ndata);
	x->sigma = AlloDouble(x->ndata);
	x->datatype = AlloInt(x->ndata);
	x->i1 = AlloInt(x->ndata);
	x->i2 = AlloInt(x->ndata);
	x->i3 = AlloInt(x->ndata);
	x->i4 = AlloInt(x->ndata);

	// read file
	for (i=0;i<x->ndata;++i)
	{
		fgets(aux,500,fp);

		if ( sscanf(aux,"%d %d %d %d %d %lf %lf",&(x->i1[i]),&(x->i2[i]),&(x->i3[i]),&(x->i4[i]),&(x->datatype[i]),&(x->expdata[i]),&(x->sigma[i])) != 7 )
		{
			if ( sscanf(aux,"%d %d %d %lf %lf",&(x->i1[i]),&(x->i2[i]),&(x->datatype[i]),&(x->expdata[i]),&(x->sigma[i])) != 5 )
			{
				fprintf(stderr,"line: %d string: %s\n",i,aux); 
				Error("Error in reading restrain file"); 
			}
		}
		else if ( x->datatype[i] != 3) Error("If it reads i3 and i4 in restrain file, datatype must be 3");
		if ( x->sigma[i] == 0 ) Error("sigma of a restraint cannot be 0 (negative: validation restraint)");
	}
	if ( NFitted(x) == 0 ) Error("No restraint to fit (all sigma are negative)");
	if ( NFitted(x) < x->ndata ) fprintf(stderr," %d restraints are for validation only\n",x->ndata-NFitted(x));

	fprintf(stderr," %d restrains read\n",x->ndata);
	fclose(fp);

	return x;
}

void FreeOptInput(struct s_optimizepot_input *x)
{
	free(x->expdata);
	free(x->sigma);
	free(x->datatype);
	free(x->i1);
	free(x->i2);
	free(x->i3);
	free(x->i4);
	free(x);

}
/***************************************************
 Get the actual values of the resrtain from chain ipol
 of polymer p
 ***************************************************/
void OP_GetRestrain(int it, struct s_polymer *p, int ipol, struct s_optimizepot_input *opi)
{
	int i,k,l,m,i1,i2,i3,i4;
	double xraw,d;

	// loop on all restrains 
	for (i=0;i<opi->ndata;++i)
	{
		
		i1 = opi->i1[i];
		i2 = opi->i2[i];
		xraw = 0.;
		
		// contact between amino acids
		if ( opi->datatype[i] == 0 )
		{
			for (k=0;k<(((p+ipol)->back)+i1)->ncontacts;++k)
				if ( *(((((p+ipol)->back)+i1)->contacts)+k) == i2 && *(((((p+ipol)->back)+i1)->contacts_p)+k) == ipol ) xraw = 1.;
		}
		// distance between atoms
		else if ( opi->datatype[i] == 1 )
		{
				
		xraw = FastSqrt( Dist2( *(*(((p+ipol)->vback)+i1)) , *(*(((p+ipol)->vback)+i2)) ), p->tables );
					
		}
		// 1/r^6 between atoms
	 	else if ( opi->datatype[i] == 2 )
		{
			d = Dist2( *(*(((p+ipol)->vback)+i1)) , *(*(((p+ipol)->vback)+i2)) );
			xraw = 1. / (d*d*d);
		}
		// contact between blocks of amino acids
		else if ( opi->datatype[i] == 3 )
		{
			i3 = opi->i3[i];		
			i4 = opi->i4[i];		
			for (k=i1;k<=i2;++k)
				for (l=i3;l<=i4;++l)
					for (m=0;m<(((p+ipol)->back)+k)->ncontacts;++m)
						if ( *(((((p+ipol)->back)+k)->contacts)+m) == l && *(((((p+ipol)->back)+k)->contacts_p)+m) == ipol ) xraw = 1.;
		}
		
		p->op->restrain[it][i] = xraw;
	}

	
	#ifdef OP_DEBUG
	 FILE *fp;
	 fp = fopen("op_log","a");
	 for (i=0;i<opi->ndata;++i)
	 	fprintf(fp,"restrains: snap=%5d restr=%4d  value=%lf  T\n",it,i,p->op->restrain[it][i]);
	 fclose(fp);
	#endif

}

/***************************************************
 Index of the contact between types a1 and a2,
 added to the list if new
 ***************************************************/
static int OP_ContactIndex(struct s_optimizepot *x, int a1, int a2)
{
	int c;

	c = x->cindex[a1][a2];
	if (c>=0) return c;

	c = x->ncontacts;
	if ( c >= OP_NCONTMAX * x->nallocont )			// enlarge the arrays of the contacts
	{
		#ifndef REALLOCATE
		Error("OP_NCONTMAX too small");
		#endif
		x->nallocont ++;
		x->it1 = (int *) realloc(x->it1, OP_NCONTMAX * x->nallocont * sizeof(int));
		x->it2 = (int *) realloc(x->it2, OP_NCONTMAX * x->nallocont * sizeof(int));
		x->fmul = (double *) realloc(x->fmul, OP_NCONTMAX * x->nallocont * sizeof(double));
		x->ftouched = (int *) realloc(x->ftouched, OP_NCONTMAX * x->nallocont * sizeof(int));
		if (!x->it1 || !x->it2 || !x->fmul || !x->ftouched) Error("Cannot reallocate the contacts of the optimization");
		memset(x->fmul+c,0,(OP_NCONTMAX * x->nallocont - c)*sizeof(double));
	}
	x->it1[c] = a1;
	x->it2[c] = a2;
	x->cindex[a1][a2] = x->cindex[a2][a1] = c;
	x->ncontacts ++;

	return c;
}

/***************************************************
 Add energy to energy table of the frame being
 recorded
 ***************************************************/
void OP_AddEnergy(struct s_polymer *p, int a1, int a2, double mul)
{
	struct s_optimizepot *x = p->op;
	int c;

	c = OP_ContactIndex(x,a1,a2);
	if (x->fmul[c] == 0) x->ftouched[x->nftouched++] = c;		// first time in this frame
	x->fmul[c] += mul;
}

static int CompareInt(const void *a, const void *b)
{
	return (*(const int *)a > *(const int *)b) - (*(const int *)a < *(const int *)b);
}

/***************************************************
 Store the contacts of the frame being recorded as
 frame x->nframes, sorted by contact index
 ***************************************************/
static void OP_EndFrame(struct s_optimizepot *x)
{
	int i,c;

	if (x->cidx==NULL || x->nnz + x->nftouched > x->nnzmax)
	{
		if (x->cidx==NULL) x->nnzmax = 16*(long long)OP_NCONTMAX;		// freed by OP_ByContact
		while (x->nnz + x->nftouched > x->nnzmax) x->nnzmax *= 2;
		x->cidx = (int *) realloc(x->cidx, x->nnzmax*sizeof(int));
		x->cmul = (double *) realloc(x->cmul, x->nnzmax*sizeof(double));
		if (!x->cidx || !x->cmul) Error("Cannot reallocate the contacts of the frames");
	}
	qsort(x->ftouched,x->nftouched,sizeof(int),CompareInt);
	for (i=0;i<x->nftouched;++i)
	{
		c = x->ftouched[i];
		x->cidx[x->nnz] = c;
		x->cmul[x->nnz] = x->fmul[c];
		x->nnz ++;
		x->fmul[c] = 0;
	}
	x->nftouched = 0;
	x->fstart[x->nframes+1] = x->nnz;
}

/***************************************************
 Turn the contacts stored by frame into contacts
 stored by pair of types (frames ascending), which
 is what the reweighting needs, and free the first
 ***************************************************/
static void OP_ByContact(struct s_optimizepot *x)
{
	int f,c;
	long long k,*next;

	x->cstart = (long long *) calloc(x->ncontacts+1,sizeof(long long));
	next = (long long *) malloc((x->ncontacts+1)*sizeof(long long));
	x->cfr = (int *) malloc((x->nnz+1)*sizeof(int));
	x->cmc = (double *) malloc((x->nnz+1)*sizeof(double));
	if (!x->cstart || !next || !x->cfr || !x->cmc) Error("Cannot allocate the contacts by pair of types");

	for (k=0;k<x->nnz;++k) x->cstart[x->cidx[k]+1] ++;
	for (c=0;c<x->ncontacts;++c) x->cstart[c+1] += x->cstart[c];
	memcpy(next,x->cstart,(x->ncontacts+1)*sizeof(long long));
	for (f=0;f<x->nframes;++f)
		for (k=x->fstart[f];k<x->fstart[f+1];++k)
		{
			x->cfr[next[x->cidx[k]]] = f;
			x->cmc[next[x->cidx[k]]] = x->cmul[k];
			next[x->cidx[k]] ++;
		}
	free(next);
	free(x->cidx); free(x->cmul);
	x->cidx = NULL; x->cmul = NULL;
	x->nnzmax = 0;
}

#ifdef OP_DEBUG
/***************************************************
 Geometric factor of contact iw in frame ifr (0 if
 absent), after OP_ByContact
 ***************************************************/
static double OP_Mul(struct s_optimizepot *x, int ifr, int iw)
{
	long long lo=x->cstart[iw], hi=x->cstart[iw+1]-1, mid;

	while (lo<=hi)
	{
		mid = (lo+hi)/2;
		if (x->cfr[mid]==ifr) return x->cmc[mid];
		if (x->cfr[mid]<ifr) lo = mid+1;
		else hi = mid-1;
	}
	return 0;
}
#endif

/***************************************************
 Record the present conformation as a frame
 ***************************************************/
void OP_RecordFrame(struct s_polymer *p, struct s_potential *pot, struct s_mc_parms *parms, int my_rank, double t)
{
	struct s_optimizepot *op = p->op;

	if (op->nframes >= op->nframesmax) Error("Too many frames recorded for the optimization of the potential");
	op->record = 1;         // fill it1,it2 and the contacts of the frame
	TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank);      // record energy
	op->record = 0;
	OP_EndFrame(op);
	#ifdef OP_DEBUG
	FILE *fp;
	long long k;
	fp = fopen("op_log","a");
	for (k=op->fstart[op->nframes];k<op->fstart[op->nframes+1];++k)
		fprintf(fp,"contacts: snap=%5d  %d-%d\tm=%lf\n",op->nframes,op->it1[op->cidx[k]],op->it2[op->cidx[k]],op->cmul[k]);
	fclose(fp);
	#endif
	OP_GetRestrain(op->nframes,p,0,parms->op_input);            // record restrains
	op->t[op->nframes] = t;						// record temperature (the current one with simulated tempering)
	op->nframes ++;
}

/***************************************************
 Frames for the checkpoint: append to fname the
 contacts and frames recorded since the last call
 (a contact as 'C' it1 it2, a frame as 'F' eold efix
 t, n entries (contact, factor), the restraints) and
 return its length. OP_Reset starts a new file
 ***************************************************/
long long OP_SaveFrames(struct s_optimizepot *x, int nrestr, char *fname)
{
	int i;
	long long k,n;
	char c;

	if (!x->fsave)
	{
		x->fsave = fopen(fname,"wb");
		if (!x->fsave) Error("Cannot open the file of the frames of the checkpoint");
	}
	for (i=x->ncsaved;i<x->ncontacts;++i)
	{
		c = 'C';
		fwrite(&c,1,1,x->fsave);
		fwrite(x->it1+i,sizeof(int),1,x->fsave);
		fwrite(x->it2+i,sizeof(int),1,x->fsave);
	}
	for (i=x->nfsaved;i<x->nframes;++i)
	{
		c = 'F';
		n = x->fstart[i+1]-x->fstart[i];
		fwrite(&c,1,1,x->fsave);
		fwrite(x->eold+i,sizeof(double),1,x->fsave);
		fwrite(x->efix+i,sizeof(double),1,x->fsave);
		fwrite(x->t+i,sizeof(double),1,x->fsave);
		fwrite(&n,sizeof(long long),1,x->fsave);
		for (k=x->fstart[i];k<x->fstart[i+1];++k)
		{
			fwrite(x->cidx+k,sizeof(int),1,x->fsave);
			fwrite(x->cmul+k,sizeof(double),1,x->fsave);
		}
		fwrite(x->restrain[i],sizeof(double),nrestr,x->fsave);
	}
	x->ncsaved = x->ncontacts;
	x->nfsaved = x->nframes;
	if (fflush(x->fsave) || ferror(x->fsave)) Error("Cannot write the file of the frames of the checkpoint");
	return (long long) ftell(x->fsave);
}

/***************************************************
 Restart: record again the contacts and frames of
 fname (cut to its length at the checkpoint), in the
 same order, so with the same indices, and continue
 the file
 ***************************************************/
void OP_LoadFrames(struct s_optimizepot *x, int nrestr, char *fname)
{
	FILE *fp;
	char c;
	int a1,a2,ci;
	long long k,n;
	double cm;

	fp = fopen(fname,"rb");
	if (!fp) Error("Cannot open the file of the frames of the checkpoint");
	while (fread(&c,1,1,fp)==1)
	{
		if (c=='C')
		{
			if (fread(&a1,sizeof(int),1,fp)!=1 || fread(&a2,sizeof(int),1,fp)!=1 || a1<0 || a2<0 || a1>=x->ntypes || a2>=x->ntypes)
				Error("Cannot read a contact in the file of the frames of the checkpoint");
			if (OP_ContactIndex(x,a1,a2) != x->ncontacts-1) Error("Contact repeated in the file of the frames of the checkpoint");
		}
		else if (c=='F')
		{
			if (x->nframes >= x->nframesmax) Error("Too many frames in the file of the frames of the checkpoint");
			if (fread(x->eold+x->nframes,sizeof(double),1,fp)!=1 || fread(x->efix+x->nframes,sizeof(double),1,fp)!=1 ||
					fread(x->t+x->nframes,sizeof(double),1,fp)!=1 || fread(&n,sizeof(long long),1,fp)!=1 || n<0 || n>x->ncontacts)
				Error("Cannot read a frame in the file of the frames of the checkpoint");
			for (k=0;k<n;++k)
			{
				if (fread(&ci,sizeof(int),1,fp)!=1 || fread(&cm,sizeof(double),1,fp)!=1 || ci<0 || ci>=x->ncontacts)
					Error("Cannot read a frame in the file of the frames of the checkpoint");
				x->ftouched[x->nftouched++] = ci;
				x->fmul[ci] = cm;
			}
			OP_EndFrame(x);
			if (fread(x->restrain[x->nframes],sizeof(double),nrestr,fp)!=(size_t)nrestr)
				Error("Cannot read a frame in the file of the frames of the checkpoint");
			x->nframes ++;
		}
		else Error("Wrong record in the file of the frames of the checkpoint");
	}
	fclose(fp);
	x->ncsaved = x->ncontacts;
	x->nfsaved = x->nframes;
	x->fsave = fopen(fname,"ab");
	if (!x->fsave) Error("Cannot open the file of the frames of the checkpoint");
}

#ifdef ACTIVE_MPI
/***************************************************
 Send/receive long arrays in pieces (MPI counts are
 int)
 ***************************************************/
#define OP_CHUNK (1<<26)
static void SendLong(void *buf, long long n, MPI_Datatype type, int size, int dest, int tag)
{
	long long i;
	for (i=0;i<n;i+=OP_CHUNK)
		MPI_Send((char *)buf+i*size, (int)(n-i<OP_CHUNK ? n-i : OP_CHUNK), type, dest, tag, MPI_COMM_WORLD);
}
static void RecvLong(void *buf, long long n, MPI_Datatype type, int size, int src, int tag)
{
	long long i;
	for (i=0;i<n;i+=OP_CHUNK)
		MPI_Recv((char *)buf+i*size, (int)(n-i<OP_CHUNK ? n-i : OP_CHUNK), type, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

/***************************************************
 op_allreplicas: rank 0 collects the frames of all
 the ranks, after its own; the contact indices of
 the other ranks are mapped to those of rank 0
 ***************************************************/
void OP_GatherFrames(struct s_optimizepot *x, struct s_mc_parms *parms, int my_rank, int nprocs)
{
	int r,i,nf,nc,nr=parms->op_input->ndata,*it1,*it2,*cidx,*map;
	long long nnz,*fstart,k;
	double *cmul,*eold,*efix,*t;

	if (my_rank != 0)
	{
		int head[2] = {x->nframes,x->ncontacts};

		MPI_Send(head,2,MPI_INT,0,100,MPI_COMM_WORLD);
		MPI_Send(&(x->nnz),1,MPI_LONG_LONG,0,101,MPI_COMM_WORLD);
		MPI_Send(x->it1,x->ncontacts,MPI_INT,0,102,MPI_COMM_WORLD);
		MPI_Send(x->it2,x->ncontacts,MPI_INT,0,103,MPI_COMM_WORLD);
		MPI_Send(x->fstart,x->nframes+1,MPI_LONG_LONG,0,104,MPI_COMM_WORLD);
		SendLong(x->cidx,x->nnz,MPI_INT,sizeof(int),0,105);
		SendLong(x->cmul,x->nnz,MPI_DOUBLE,sizeof(double),0,106);
		MPI_Send(x->eold,x->nframes,MPI_DOUBLE,0,107,MPI_COMM_WORLD);
		MPI_Send(x->efix,x->nframes,MPI_DOUBLE,0,108,MPI_COMM_WORLD);
		MPI_Send(x->t,x->nframes,MPI_DOUBLE,0,109,MPI_COMM_WORLD);
		for (i=0;i<x->nframes;++i)
			MPI_Send(x->restrain[i],nr,MPI_DOUBLE,0,110,MPI_COMM_WORLD);
		return;
	}

	for (r=1;r<nprocs;++r)
	{
		int head[2];

		MPI_Recv(head,2,MPI_INT,r,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		nf = head[0];
		nc = head[1];
		MPI_Recv(&nnz,1,MPI_LONG_LONG,r,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		if (x->nframes + nf > x->nframesmax) Error("Too many frames gathered for the optimization of the potential");

		it1 = AlloInt(nc+1);
		it2 = AlloInt(nc+1);
		map = AlloInt(nc+1);
		fstart = (long long *) malloc((nf+1)*sizeof(long long));
		cidx = (int *) malloc((nnz+1)*sizeof(int));
		cmul = (double *) malloc((nnz+1)*sizeof(double));
		eold = AlloDouble(nf+1);
		efix = AlloDouble(nf+1);
		t = AlloDouble(nf+1);
		if (!fstart || !cidx || !cmul) Error("Cannot allocate the frames of the other replicas");

		MPI_Recv(it1,nc,MPI_INT,r,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(it2,nc,MPI_INT,r,103,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(fstart,nf+1,MPI_LONG_LONG,r,104,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		RecvLong(cidx,nnz,MPI_INT,sizeof(int),r,105);
		RecvLong(cmul,nnz,MPI_DOUBLE,sizeof(double),r,106);
		MPI_Recv(eold,nf,MPI_DOUBLE,r,107,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(efix,nf,MPI_DOUBLE,r,108,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(t,nf,MPI_DOUBLE,r,109,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		for (i=0;i<nc;++i) map[i] = OP_ContactIndex(x,it1[i],it2[i]);

		for (i=0;i<nf;++i)
		{
			for (k=fstart[i];k<fstart[i+1];++k)
			{
				x->ftouched[x->nftouched++] = map[cidx[k]];
				x->fmul[map[cidx[k]]] = cmul[k];
			}
			OP_EndFrame(x);
			x->eold[x->nframes] = eold[i];
			x->efix[x->nframes] = efix[i];
			x->t[x->nframes] = t[i];
			MPI_Recv(x->restrain[x->nframes],nr,MPI_DOUBLE,r,110,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			x->nframes ++;
		}

		free(it1); free(it2); free(map); free(fstart); free(cidx); free(cmul);
		free(eold); free(efix); free(t);
	}
}
#endif

/***************************************************
 Denominators of the weights of the frames.
 Frames from one temperature: lnden = -eold/t.
 Frames from K temperatures (op_allreplicas): MBAR,
 lnden_n = ln sum_k N_k exp(f_k - eold_n/T_k), with
 the free energies f_k solved self-consistently
 (Shirts and Chodera, J. Chem. Phys. 129, 124105
 (2008)), so that exp(-enew/op_T - lnden) is the
 weight of the frame at op_T with the new energy
 ***************************************************/
#define OP_MAXTEMP 64
static void OP_Denominators(struct s_optimizepot *x, struct s_mc_parms *parms)
{
	int n,k,K=0,iter,kk[OP_MAXTEMP];
	double T[OP_MAXTEMP],lnN[OP_MAXTEMP],f[OP_MAXTEMP],fnew[OP_MAXTEMP],a,amax,sum,dmax;

	for (n=0;n<x->nframes;++n)
	{
		for (k=0;k<K;++k) if (fabs(x->t[n]-T[k])<EPSILON) break;
		if (k==K)
		{
			if (K==OP_MAXTEMP) Error("Too many temperatures in the frames of the optimization");
			T[K] = x->t[n];
			kk[K] = 0;
			K ++;
		}
		kk[k] ++;
	}

	if (K==1 || !parms->op_allreplicas)
	{
		for (n=0;n<x->nframes;++n) x->lnden[n] = -x->eold[n]/x->t[n];
		return;
	}

	for (k=0;k<K;++k) { lnN[k] = log((double)kk[k]); f[k] = 0; }
	for (iter=0;iter<100000;++iter)
	{
		// denominators with the present f
		for (n=0;n<x->nframes;++n)
		{
			amax = -9E19;
			for (k=0;k<K;++k)
			{
				a = lnN[k] + f[k] - x->eold[n]/T[k];
				if (a>amax) amax = a;
			}
			sum = 0;
			for (k=0;k<K;++k) sum += exp(lnN[k] + f[k] - x->eold[n]/T[k] - amax);
			x->lnden[n] = amax + log(sum);
		}
		// new f_k = -ln sum_n exp(-eold_n/T_k - lnden_n), with f_0 = 0
		for (k=0;k<K;++k)
		{
			amax = -9E19;
			for (n=0;n<x->nframes;++n)
				if (-x->eold[n]/T[k] - x->lnden[n] > amax) amax = -x->eold[n]/T[k] - x->lnden[n];
			sum = 0;
			for (n=0;n<x->nframes;++n) sum += exp(-x->eold[n]/T[k] - x->lnden[n] - amax);
			fnew[k] = -amax - log(sum);
		}
		dmax = 0;
		for (k=0;k<K;++k)
		{
			fnew[k] -= fnew[0];
			if (fabs(fnew[k]-f[k])>dmax) dmax = fabs(fnew[k]-f[k]);
			f[k] = fnew[k];
		}
		if (dmax<1E-9) break;
	}

	fprintf(stderr,"> MBAR on %d temperatures, %d iterations%s\n>   T\tframes\tf\n",K,iter+1,(dmax<1E-9 ? "" : " (NOT CONVERGED)"));
	for (k=0;k<K;++k) fprintf(stderr,">   %lf\t%d\t%lf\n",T[k],kk[k],f[k]);
}

/***************************************************
 Random search in energy space
 ***************************************************/
void OP_SamplePotential(struct s_polymer *p, struct s_mc_parms *parms, int ntypes, double **u, int irun, double **ematold)
{
	int istep,irunstart,iacc=0,a1,a2,iw,ifr,i,j,n=0,nstep=0,ntry,accept,nneff=0;
	double chi2,chi2old,chi2best,deltae,eav=0,eav2=0,de2=0,ndata,neff0,chi2val0;
	double **ebest;
	char aux[500];
	FILE *fp, *fc;
	
	irunstart=parms->chi2start;

	if ( p->op->nframes == 0 ) Error("Cannot optimize potential: no frame recorded for OP_SamplePotential");
	if ( p->op->ncontacts == 0 ) Error("Cannot optimize potential: no contact recorded for OP_SamplePotential");

	//keep old matrix for comparison
	CopyDoubleMatrix(u,ematold,ntypes,ntypes);

	fprintf(stderr,"> %d frames, %d pairs of types in contact, %.1f per frame (%.0f MB)\n",p->op->nframes,p->op->ncontacts,
		(double)p->op->nnz/p->op->nframes,(double)p->op->nnz*(sizeof(int)+sizeof(double))/1024/1024);
	OP_ByContact(p->op);
	OP_Denominators(p->op,parms);


	// initial chi2 (reduced with the number of fitted restraints)
	chi2old = OP_function(u,p->op,parms->op_input,parms);
	ndata = NFitted(parms->op_input);
	neff0 = p->op->neff;
	chi2val0 = Chi2Validation(p->op->rest_av,parms->op_input);
	{
		double dmax=0;				// with the original matrix, the energy rebuilt from the contacts is the recorded one
		for (ifr=0;ifr<p->op->nframes;++ifr)
			if (fabs(p->op->enew[ifr]-p->op->eold[ifr])>dmax) dmax = fabs(p->op->enew[ifr]-p->op->eold[ifr]);
		fprintf(stderr,"> largest |Enew-Eold| with the original matrix = %g\n",dmax);
	}

	// with op_tchi2>0 the search can go uphill: keep the best matrix found
	chi2best = chi2old;
	ebest = AlloDoubleMatrix(ntypes,ntypes);
	CopyDoubleMatrix(u,ebest,ntypes,ntypes);
	
//	fprintf(stderr,"Potential Optimization \n> Initial value of Chi2 = %lf | Reduced = %lf\n",chi2old,chi2old/(parms->op_input->ndata));
	fprintf(stderr,"Potential Optimization \n> Initial value of Chi2 %lf\n",chi2old/ndata);
	if (ndata < parms->op_input->ndata) fprintf(stderr,"> Initial value of Chi2 of validation restraints %lf\n",chi2val0);

	fprintf(stderr,"> History\t\t\tChi2");
	//print chi2
	if(parms->nrun>1)
	{
		snprintf(aux,sizeof(aux),"chi2.dat");
		fc = fopen(aux,"a");
		if (ndata < parms->op_input->ndata)			// run, chi2 of the fitted and of the validation restraints
			fprintf(fc,"%d\t%lf\t%lf\n",irun+irunstart, chi2old/ndata, chi2val0);
		else
			fprintf(fc,"%d\t%lf\n",irun+irunstart, chi2old/ndata);
		fclose(fc);
	}
	
	// minimize
	MGInfo_Optimize(parms,irun,0,parms->op_itermax,chi2old/ndata);
	for (istep=0;istep<parms->op_itermax;++istep)
	{
		iw = irand(p->op->ncontacts);
		deltae = (frand()-0.5) * parms->op_step;

		ntry = 0;
		do { deltae = (frand()-0.5) * parms->op_step; ntry++; }			// generate a random deltae
		while ( ( u[p->op->it1[iw]][p->op->it2[iw]]+deltae < parms->op_emin ||
			 u[p->op->it1[iw]][p->op->it2[iw]]+deltae > parms->op_emax ) && ntry < 1000 );		// check if matrix element within boundaries
		if ( u[p->op->it1[iw]][p->op->it2[iw]]+deltae < parms->op_emin ||
			 u[p->op->it1[iw]][p->op->it2[iw]]+deltae > parms->op_emax )	// element out of [op_emin,op_emax] by more than op_step/2
		{
			nstep ++;
			continue;
		}

		chi2 = OP_functionDiff(iw,deltae,p->op,parms->op_input,parms);

		// Metropolis on the reduced chi2 at temperature op_tchi2 (op_tchi2=0: only downhill steps)
		if (chi2 < chi2old) accept = 1;
		else if (parms->op_tchi2 > 0) accept = ( frand() < exp( -(chi2-chi2old)/ndata/parms->op_tchi2 ) );
		else accept = 0;
		// the reweighted averages are reliable only while enough frames keep a weight (relative to the start:
		// with frames from many temperatures the effective fraction is already below 1 with the original matrix)
		if (accept && p->op->neff < parms->op_neff * neff0) { accept = 0; nneff++; }

		// accept
		if (accept)
		{
			a1 = p->op->it1[iw];
			a2 = p->op->it2[iw];
			u[a1][a2] += deltae;
			if (a1 != a2) u[a2][a1] += deltae;			// a diagonal element must change only once
			chi2old = chi2;
			iacc ++;
			if (chi2 < chi2best)
			{
				chi2best = chi2;
				CopyDoubleMatrix(u,ebest,ntypes,ntypes);
			}
		}
		// reject
		else OP_UndoDiff(iw,deltae,p->op,parms->op_input,parms);	// return to previous value

		nstep ++;

		// print something
	//	if (!(istep%parms->op_print)) fprintf(stderr,"   %5d\t\t%lf\t%lf\n",istep,chi2old, chi2old/(parms->op_input->ndata));
		if (!(istep%parms->op_print)) fprintf(stderr,"   %5d\t\t%lf\n",istep,chi2old/ndata);
		MGInfo_Optimize(parms,irun,istep+1,parms->op_itermax,chi2old/ndata);
		// exit condition
		if (chi2old<parms->op_stop) break;
	}

	// return the best matrix, with its averages (rest_av of the last proposal otherwise)
	CopyDoubleMatrix(ebest,u,ntypes,ntypes);
	FreeDoubleMatrix(ebest,ntypes);
	chi2old = OP_function(u,p->op,parms->op_input,parms);

	// calculate the properties of the new matrix
	for (i=0;i<ntypes;++i)
		for (j=i+1;j<ntypes;++j) 
		{ 
			eav += u[i][j]; 
			eav2 += u[i][j]*u[i][j]; 
			de2 += (u[i][j]-ematold[i][j])*(u[i][j]-ematold[i][j]);
			n++; 
		}


	fprintf(stderr,"> Result of the sampling: Chi2=%lf\n",chi2old/ndata);
	if (ndata < parms->op_input->ndata) fprintf(stderr,"> Chi2 of validation restraints (reweighted) = %lf\n",Chi2Validation(p->op->rest_av,parms->op_input));
	fprintf(stderr,"> effective fraction of frames = %lf (at the start %lf)\n",p->op->neff,neff0);
	if (nstep) fprintf(stderr,"> accepted changes = %lf\n",(double)iacc/nstep);
	if (nstep && parms->op_neff>0) fprintf(stderr,"> changes rejected for op_neff = %lf\n",(double)nneff/nstep);
	fprintf(stderr,"> average E = %lf +/- %lf\n",eav/n,sqrt(eav2/n-eav*eav/n/n));
	fprintf(stderr,"> rms difference = %lf\n\n",sqrt(de2/n));
	

	// print to file
	snprintf(aux,sizeof(aux),"restraints_%d.dat",irun);
	fp = fopen(aux,"w");
	for (i=0;i<parms->op_input->ndata;++i)
		fprintf(fp,"%d\t%lf\t%lf\t%lf\n",i,p->op->rest_av[i],parms->op_input->expdata[i],parms->op_input->sigma[i]);
	fclose(fp);

	// reset structures
	for (ifr=0;ifr<p->op->nframes;++ifr)
	{
		p->op->eold[ifr] = 0;
		p->op->efix[ifr] = 0;
		p->op->enew[ifr] = 0;
		p->op->t[ifr] = -1;
		for (i=0;i<parms->op_input->ndata;++i)  p->op->restrain[ifr][i] = 0;
	}
	for (i=0;i<parms->op_input->ndata;++i) p->op->rest_av[i] = 0;
	OP_Reset(p->op);
	p->op->icount = 0;
}


/***************************************************
 Weights of all the frames from their energies enew,
 with the sums over the frames for the averages
 ***************************************************/
static void OP_Weights(struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms)
{
	double emax=-9E19,boltz;
	int ifr,ir;

	for (ir=0;ir<in->ndata;++ir) x->sw[ir]=0.;
	x->z = x->z2 = 0.;

	for (ifr=0;ifr<x->nframes;++ifr)
	{
		// calculate emax to avoid overflow in the exp
		if (-x->enew[ifr]/parms->op_T - x->lnden[ifr] > emax) 
						emax = -x->enew[ifr]/parms->op_T - x->lnden[ifr];
	}
	x->eshift = emax;

	// calculate thermal averages
	for (ifr=0;ifr<x->nframes;++ifr)
	{
		boltz = exp(-x->enew[ifr]/parms->op_T - x->lnden[ifr] - emax);
		x->w[ifr] = boltz;

		// make the average
		for (ir=0;ir<in->ndata;++ir) 
			x->sw[ir] +=  x->restrain[ifr][ir] * boltz;
                                                                                                                            
		x->z += boltz;
		x->z2 += boltz*boltz;
	}
}

/***************************************************
 Averages of the restraints, effective fraction of
 frames and chi2 from the sums of the weights
 ***************************************************/
static double OP_Averages(struct s_optimizepot *x, struct s_optimizepot_input *in)
{
	int ir;

	for (ir=0;ir<in->ndata;++ir) x->rest_av[ir] = x->sw[ir] / x->z;
	x->neff = x->z*x->z/x->z2/x->nframes;

	return Chi2(x->rest_av,in->expdata,in->sigma,in->ndata);
}

/***************************************************
 Calculate the chi2 with experimental data, updating
 the vector rest_av and e_new
 ***************************************************/
double OP_function(double **e, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms)
{
	double mul;
	int ifr,icont,it1,it2;
	long long k;

	// calculate energies
	for (ifr=0;ifr<x->nframes;++ifr)
		x->enew[ifr] = x->efix[ifr];			// one-body energy
	for (icont=0;icont<x->ncontacts;++icont)		// new two-body energy, contact by contact (in each frame in the order of the contacts)
	{
		it1 = x->it1[icont];
		it2 = x->it2[icont];
		for (k=x->cstart[icont];k<x->cstart[icont+1];++k)
		{
			mul = x->cmc[k];
			x->enew[x->cfr[k]] += e[it1][it2] * mul;
		}
	}
	x->ndiff = 0;
	OP_Weights(x,in,parms);

	return OP_Averages(x,in);
}

/***************************************************
 Calculate chi2 between two vectors
 ***************************************************/
double Chi2(double *x, double *xexp, double *sigma, int n)
{
	int i;
	double chi2=0.;

	for (i=0;i<n;++i)
		if (sigma[i]>0)						// validation restraints have sigma<0
			chi2 += (x[i] - xexp[i])*(x[i] - xexp[i]) / ( sigma[i] * sigma[i] );

	return chi2;
}

/***************************************************
 Reduced chi2 of the validation restraints (sigma<0)
 ***************************************************/
double Chi2Validation(double *x, struct s_optimizepot_input *in)
{
	int i,n=0;
	double chi2=0.;

	for (i=0;i<in->ndata;++i)
		if (in->sigma[i]<0)
		{
			chi2 += (x[i] - in->expdata[i])*(x[i] - in->expdata[i]) / ( in->sigma[i] * in->sigma[i] );
			n++;
		}

	return n ? chi2/n : 0.;
}

/***************************************************
 Calculate the chi2 with experimental data, given
 the energy difference in the contact iw
 it needs x->enew and changes x->enew[ifr]
 Only the frames that contain iw change their weight:
 the sums over the frames are updated with them.
 They are computed again from scratch every
 OP_NDIFFMAX steps (the rounding errors of the
 updates do not pile up) and when a weight would
 leave the range of the doubles.
 ***************************************************/
double OP_functionDiff(int iw, double deltae, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms)
{
	double boltz,dw,a,*r;
	int f,ir,full;
	long long k,m;

	full = ( ++x->ndiff >= OP_NDIFFMAX );
	for (k=x->cstart[iw];k<x->cstart[iw+1] && !full;++k)
	{
		f = x->cfr[k];
		a = -(x->enew[f] + deltae * x->cmc[k])/parms->op_T - x->lnden[f] - x->eshift;
		if (a > OP_EXPMAX) full = 1;
	}

	x->fullold = full;
	if (full)
	{
		for (k=x->cstart[iw];k<x->cstart[iw+1];++k)
			x->enew[x->cfr[k]] += deltae * x->cmc[k];
		x->ndiff = 0;
		OP_Weights(x,in,parms);
		return OP_Averages(x,in);
	}

	// keep the sums to undo the step
	x->zold = x->z;
	x->z2old = x->z2;
	for (ir=0;ir<in->ndata;++ir) x->swold[ir] = x->sw[ir];

	// the frames that contain the iw contact change their energy
	for (k=x->cstart[iw],m=0;k<x->cstart[iw+1];++k,++m)
	{
		f = x->cfr[k];
		x->enew[f] += deltae * x->cmc[k];
		boltz = exp(-x->enew[f]/parms->op_T - x->lnden[f] - x->eshift);
		x->wold[m] = x->w[f];
		dw = boltz - x->w[f];
		x->w[f] = boltz;
		x->z += dw;
		x->z2 += boltz*boltz - x->wold[m]*x->wold[m];
		r = x->restrain[f];
		for (ir=0;ir<in->ndata;++ir) x->sw[ir] += dw * r[ir];
	}

	// all the weights became too small to be summed with precision
	if (x->z < OP_ZMIN)
	{
		x->ndiff = 0;
		x->fullold = 1;
		OP_Weights(x,in,parms);
	}

	#ifdef OP_DEBUG
	 FILE *fp;
	 fp = fopen("op_log","a");
	 fprintf(fp,"minim:  iw=%d deltae=%lf\n",iw,deltae);
 	 for (f=0;f<x->nframes;++f)
	 	fprintf(fp,"\teold[%d]=%lf enew[%d]=%lf efix[%d]=%lf mul[%d]=%lf",f,x->eold[f],f,x->enew[f],f,x->efix[f],f,OP_Mul(x,f,iw));
         fprintf(fp,"\neshift=%lf\n",x->eshift);
	 fclose(fp);
	#endif

	return OP_Averages(x,in);
}

/***************************************************
 Undo the last OP_functionDiff (same iw and deltae)
 ***************************************************/
void OP_UndoDiff(int iw, double deltae, struct s_optimizepot *x, struct s_optimizepot_input *in, struct s_mc_parms *parms)
{
	int ir;
	long long k,m;

	for (k=x->cstart[iw];k<x->cstart[iw+1];++k)
		x->enew[x->cfr[k]] -= deltae * x->cmc[k];

	if (x->fullold) OP_Weights(x,in,parms);
	else
	{
		for (k=x->cstart[iw],m=0;k<x->cstart[iw+1];++k,++m) x->w[x->cfr[k]] = x->wold[m];
		x->z = x->zold;
		x->z2 = x->z2old;
		for (ir=0;ir<in->ndata;++ir) x->sw[ir] = x->swold[ir];
	}
	OP_Averages(x,in);
}

