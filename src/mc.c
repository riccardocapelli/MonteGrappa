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
 * mc.c
 *
 *  Created on: Sep 20, 2010
 *      Author: guido
 */



#include "montegrappa.h"
#include <sys/stat.h>
/*****************************************************************************
 Monte Carlo loop
 *****************************************************************************/
volatile sig_atomic_t softexit;

// current length of an output file (-1 if it is not open)
static long long FileLength(FILE *fp)
{
	if (!fp) return -1;
	if (fflush(fp)) Error("Cannot write output file");
	return (long long) ftell(fp);
}

/*****************************************************************************
 Write the checkpoint of the state at the end of step istep-1 of run irun,
 with the lengths of the output files the restart cuts them to
 *****************************************************************************/
static void SaveCheckpoint(struct s_checkpoint *chk, char *fname, struct s_polymer *p, struct s_potential *pot, struct s_mc_parms *parms, int irun,
		unsigned long long istep, int iprinttrj, int iprintlog, int iprinte, int *mcount, long long *macc, long long *mdone,
		int ptcount, FILE *ftrj, FILE *fe, FILE *fexchange, int my_rank)
{
	int i;

	chk->istep = istep;
	chk->irun = irun;
	GetRNGState(chk->rng);
	chk->iprinttrj = iprinttrj;
	chk->iprintlog = iprintlog;
	chk->iprinte = iprinte;
	chk->icheckpoint = 0;
	chk->acc = parms->acc;
	chk->mov = parms->mov;
	for(i=0;i<NMOVES;++i)
	{
		chk->mcount[i] = mcount[i];
		chk->macc[i] = macc[i];
		chk->mdone[i] = mdone[i];
	}
	chk->ptcount = ptcount;
	chk->etot = p->etot;
	chk->off[0] = FileLength(ftrj);
	chk->off[1] = strcmp(parms->fne,"") ? FileLength(fe) : -1;
	chk->off[2] = FileLength(parms->fobs);
	chk->off[3] = FileLength(fexchange);
	chk->off[4] = chk->off[5] = chk->off[6] = chk->off[7] = -1;
	#ifdef STEMPERING
	if (parms->stempering)
	{
		chk->off[4] = FileLength((parms->p)->st_ftout);
		if ((parms->p)->st_ttarget_harvest) chk->off[5] = FileLength((parms->p)->st_pdbf);
	}
	#endif
	if (strcmp(parms->op_minim,"none") && (my_rank==0 || parms->op_allreplicas))	// the frames first: the checkpoint refers to them
	{
		char fop[FNAMEMAX];

		CheckpointOPName(fop,sizeof(fop),my_rank);
		chk->off[6] = OP_SaveFrames(p->op,parms->op_input->ndata,fop);
		if (my_rank==0)							// chi2.dat is written at the end of each run
		{
			struct stat sb;
			chk->off[7] = stat("chi2.dat",&sb) ? 0 : (long long) sb.st_size;
		}
	}
	WriteCheckpoint(fname,chk,p,parms,pot);
}

/*****************************************************************************
 Variables of the MC loop that the checkpoint keeps
 *****************************************************************************/
static void SaveLoopState(struct s_checkpoint *chk, double t, int st_iprint, int st_o, int nconf, int op_icount,
		int *ex_count, int *ex_acc, int my_rank, int ntemp)
{
	chk->t = t;
	chk->st_iprint = st_iprint;
	chk->st_o = st_o;
	chk->nconf = nconf;
	chk->op_icount = op_icount;
	chk->exc = chk->exa = 0;
	if (ex_count && my_rank<ntemp-1)
	{
		chk->exc = ex_count[my_rank];
		chk->exa = ex_acc[my_rank];
	}
}

/*****************************************************************************
 Debug aid (keyword echeck): compare the incrementally updated energy with a
 full recalculation after a move, report the first mismatches of each move
 type and resynchronize p->etot so that later moves are tested independently
 *****************************************************************************/
static void CheckEnergyDrift(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, struct s_mc_parms *parms, int move, unsigned long long istep, int my_rank)
{
	static int nreported[NMOVES],nreported_old[NMOVES],nreported_side[NMOVES],nreported_rev[NMOVES];
	double etrue;
	int ic,i,k,q,found;
	struct s_back *a,*b;

	// after a move, oldp must be an exact copy of p (positions and contacts)
	for (ic=0;ic<parms->npol;++ic)
		for (i=0;i<(p+ic)->nback;++i)
		{
			a = ((p+ic)->back)+i;
			b = ((oldp+ic)->back)+i;
			found = (a->ncontacts == b->ncontacts && a->pos.x == b->pos.x && a->pos.y == b->pos.y && a->pos.z == b->pos.z);
			for (k=0;k<a->ncontacts && found;++k)
			{
				found = 0;
				for (q=0;q<b->ncontacts;++q)
					if (a->contacts[k]==b->contacts[q] && a->contacts_p[k]==b->contacts_p[q] && a->e[k]==b->e[q]) found = 1;
			}
			if (!found && nreported_old[move]++ < 5)
				fprintf(stderr,"ECHECK: oldp differs from p after move %d at step %llu (chain %d residue %d)\n",move,istep,ic,i);
		}

	// every contact entry knows where its partner keeps the reverse one (contacts_r)
	if ((k = ContactsReverseErrors(p,parms->npol)) + (q = ContactsReverseErrors(oldp,parms->npol)) > 0 && nreported_rev[move]++ < 5)
		fprintf(stderr,"ECHECK: %d (p) and %d (oldp) contact entries with a wrong reverse index after move %d at step %llu\n",k,q,move,istep);

	// every sidechain atom where its rotamer puts it (a move must rebuild all the sidechains built from what it moved)
	for (ic=0;ic<parms->npol && !parms->nosidechains;++ic)
		for (i=0;i<(p+ic)->nback;++i)
		{
			struct s_rotamers *r;
			struct vector v;
			int out;

			a = ((p+ic)->back)+i;
			for (k=0;k<a->nside;++k)
			{
				r = a->side[k].rot + a->irot;
				v = Spherical2Cartesian(*((p+ic)->vback[r->b1]),*((p+ic)->vback[r->b2]),*((p+ic)->vback[r->b3]),r->ang,p->tables,&out);
				if (out && Dist2(v,a->side[k].pos) > 1e-12 && nreported_side[move]++ < 5)
					fprintf(stderr,"ECHECK: sidechain of residue %d (chain %d) does not match its rotamer after move %d at step %llu\n",i,ic,move,istep);
			}
		}

	// before TotalEnergy, which rebuilds them
	if (p->cells && (k = CellsCheck(p)) > 0)
		fprintf(stderr,"ECHECK: cell list out of date for %d residues after move %d at step %llu\n",k,move,istep);

	etrue = TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,0,my_rank);
	if (fabs(etrue - p->etot) > 1e-6*(1.+fabs(etrue)))
	{
		if (nreported[move]++ < 5)
			fprintf(stderr,"ECHECK: energy drift after move %d at step %llu: Etot=%lf Etot(true)=%lf\n",move,istep,p->etot,etrue);
		p->etot = etrue;
	}
}


int Do_MC(struct s_polymer *p, struct s_polymer *fragment, struct s_polymer *replica, struct s_polymer *native, struct s_potential *pot, struct s_mc_parms *parms, FILE *ftrj, FILE *fe, struct s_polymer *oldp, FILE *fproc,int irun, struct s_checkpoint *chk, struct s_mpi_parms *mpiparms)
{
	int i,ok=0,iprinttrj=0,iprintlog=0,iprinte=0,icheckpoint=0,ntm=0,ci,mcount[NMOVES];
	long long macc[NMOVES],mdone[NMOVES];
//	int anneal_count=0, anneal_status=0;

	double t;
	unsigned long long istep=0;
	char fchkp[FNAMEMAX];
		
	if(parms->shell==1)
		parms->ishell=0;
	
	#ifdef OPTIMIZEPOT
	struct s_optimizepot *op = p->op;
	op->icount = (parms->record_native && !(chk->restart && chk->istep>0)) ? parms->op_deltat : 0;	// record_native: the starting structure is the first frame
	OP_Reset(op);
	op->record = 0;
	#endif

	int st_iprint=0;		// simulated tempering
	int st_o=0;
	int nconf=0;
       

	#ifdef ACTIVE_MPI
	int my_rank=mpiparms->my_rank;
	MPI_Datatype Backtype=mpiparms->Back_mpi;
    MPI_Datatype Sidetype=mpiparms->Side_mpi;
	MPI_Datatype Rottype=mpiparms->Rot_mpi;
	MPI_Status astatus=mpiparms->astatus;
	
	int ptempering_count=0;
	int ex_count[(parms->ntemp)-1],ex_acc[(parms->ntemp)-1];
	for (i=0;i<parms->ntemp-1;++i) { ex_count[i]=0; ex_acc[i]=0; }
	t=parms->T[my_rank];

	char fexname[100];
	// one file per rank (they would overwrite each other) and per run, as the energy files
	if (parms->nrun>1) snprintf(fexname,sizeof(fexname),"exchange_run%d_proc%d.dat",irun,my_rank);
	else snprintf(fexname,sizeof(fexname),"exchange_proc%d.dat",my_rank);
	if (chk->restart) TruncateOutput(fexname,chk->off[3]);
	FILE *fexchange=fopen(fexname,chk->restart ? "a" : "w");		// a restart continues the file of the run
	if (!fexchange) Error("Cannot open exchange file");
	#else
	int my_rank=0;
	FILE *fexchange=NULL;
	int *ex_count=NULL,*ex_acc=NULL;
	t=parms->T;
	#endif
	softexit =0;
	int stop=0;		// softexit as agreed by all ranks
	
	// print debug info
	if (parms->debug>0)
		fprintf(stderr,"Starting step loop....\n");
	if (parms->debug>1)
	{
		fprintf(stderr,"Moves = %d %d %d %d\n",parms->movetype[0],parms->movetype[1],
				parms->movetype[2],parms->movetype[3]);
		for (ci=0;ci<parms->npol;++ci)
			for (i=0;i<(p+ci)->nback;++i) ntm += (((p+ci)->back)+i)->move;
		fprintf(stderr,"Backbone atoms which can be moved = %d\n",ntm);
	}
	
	CopyAllPolymers(p,oldp,parms->npol,parms->nosidechains,parms->nosidechains);
	if (pot->pp) SetPairParameters(pot,pot->npptypes);		// the matrix of the pairs may have changed since SetPairCutoff
	parms->acc = 0;
	parms->mov = 0;
	
	for(i=0;i<NMOVES;++i)
	{
		mcount[i]=0;
		macc[i]=0;
		mdone[i]=0;
	}

	// continue from a checkpoint: restore the loop state saved at the end of a step
	CheckpointName(fchkp,sizeof(fchkp),my_rank);
	if (chk->restart)
	{
		istep = chk->istep;
		iprinttrj = chk->iprinttrj;
		iprintlog = chk->iprintlog;
		iprinte = chk->iprinte;
		icheckpoint = chk->icheckpoint;
		parms->acc = chk->acc;
		parms->mov = chk->mov;
		for(i=0;i<NMOVES;++i)
		{
			mcount[i] = chk->mcount[i];
			macc[i] = chk->macc[i];
			mdone[i] = chk->mdone[i];
		}
		#ifdef ACTIVE_MPI
		ptempering_count = chk->ptcount;
		if (my_rank<parms->ntemp-1) { ex_count[my_rank] = chk->exc; ex_acc[my_rank] = chk->exa; }
		#endif
		SetRNGState(chk->rng);
		t = chk->t;
		#ifdef STEMPERING
		st_iprint = chk->st_iprint;
		st_o = chk->st_o;
		nconf = chk->nconf;
		if (parms->stempering)
		{
			ReadCheckpointST(fchkp,parms->p);
			TruncateOutput((parms->p)->st_nftout,chk->off[4]);	// opened in append by ReadMcParms
			if ((parms->p)->st_ttarget_harvest)
			{
				TruncateOutput((parms->p)->st_pdbnf,chk->off[5]);
				(parms->p)->st_pdbf = fopen((parms->p)->st_pdbnf,"a");
				if (!(parms->p)->st_pdbf) Error("Cannot open pdbnf file for harvesting");
			}
		}
		#endif
		op->icount = chk->op_icount;
		if (strcmp(parms->op_minim,"none") && (my_rank==0 || parms->op_allreplicas))
		{
			char fop[FNAMEMAX];

			CheckpointOPName(fop,sizeof(fop),my_rank);
			TruncateOutput(fop,chk->off[6]);
			OP_LoadFrames(op,parms->op_input->ndata,fop);
		}
	}
	time_t tcheckpoint = time(NULL);
	#ifndef ACTIVE_MPI
	int ptempering_count=0;
	#endif
	if (!chk->restart)						// a restart of this run starts at least from here
	{
		SaveLoopState(chk,t,st_iprint,st_o,nconf,op->icount,ex_count,ex_acc,my_rank,parms->ntemp);
		SaveCheckpoint(chk,fchkp,p,pot,parms,irun,istep,iprinttrj,iprintlog,iprinte,mcount,macc,mdone,ptempering_count,ftrj,fe,fexchange,my_rank);
	}
	if (my_rank==0)
	{
		MGInfo_StartRun(parms,istep);
		MGInfo_MC(parms,irun,istep,p->etot,t,1);
	}

 
	#ifdef OPTIMIZEPOT
	if(my_rank==0 || parms->op_allreplicas)				// op_allreplicas: every replica records its frames
	{
		if (strcmp(parms->op_minim,"none")  && ok>-1 && !chk->restart)		// the frame of the start (a restart continues the loop)
		{
		#ifdef STEMPERING
			if (parms->stempering)
			{
				if ((parms->p)->st_ttarget_harvest==1 && st_iprint>=(parms->p)->st_printpdb && t <= (parms->p)->st_ttarget + EPSILON && (parms->p)->st_nm==1)
				{
					OP_RecordFrame(p,pot,parms,my_rank,t);
					st_iprint = 0;
				}
			if ((parms->p)->st_ttarget_harvest>0 && ok>-1 && t <= (parms->p)->st_ttarget + EPSILON )
				st_iprint ++;
            if (st_o == 1)
				t = (parms->p)->st_temp[(parms->p)->st_itemp];
			}
			else if (op->icount == parms->op_deltat)
			#else
			if (op->icount == parms->op_deltat)
			#endif
			{
				#ifdef ACTIVE_MPI
				if (my_rank==0) fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
						(parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T[my_rank]);
				#else
				fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
						(parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T);
				#endif
				OP_RecordFrame(p,pot,parms,my_rank,t);
				op->icount=0;
			}
			if (istep >= parms->op_wait)
				op->icount ++;
		}
	}
	#endif


	////////////////////
	// MC LOOP        //
	////////////////////
	while (istep<parms->nstep && stop==0)			// a restart from the last step does nothing
	{
	

		#ifdef DEBUG
		if (parms->debug>1)
			fprintf(stderr,"-----------------\nStep %llu\n",istep);
		#endif
		
		// make a move (one per allowed type)
		if (mcount[0] == parms->movetype[0])						// flip of backbone
		{
			ok = MoveBackboneFlip(p,oldp,parms,pot,istep,parms->debug,t);
			if (ok>-1) mcount[0] = 0;
			if (ok==1) macc[0]++;
			mdone[0]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,0,istep,my_rank);
		}

	        if (mcount[3] == parms->movetype[3])                                            // move of sidechain
                {
                        ok = MoveSidechain(p,oldp,parms,pot,istep,parms->debug,t);
                        if (ok==-1) ok = 0;						// no residue with rotamers found: just skip this move
                        if (ok>-1) mcount[3] = 0;
                        if (ok==1) macc[3]++;
                        mdone[3]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,3,istep,my_rank);
                }

               if (mcount[4] == parms->movetype[4])                                            // loose pivot
                {
                        ok = MoveLoosePivot(p,oldp,pot,parms->nmul_lpivot,parms,t);
                        if (ok>-1) mcount[4] = 0;
                        if (ok==1) macc[4]++;
                        mdone[4]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,4,istep,my_rank);
                }

		if (mcount[5] == parms->movetype[5])                                            // multiple flip
                {
                        ok = MoveMultipleFlip(p,oldp,parms,pot,istep,parms->debug,t);
                        if (ok>-1) mcount[5] = 0;
                        if (ok==1) macc[5]++;
                        mdone[5]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,5,istep,my_rank);
                }

		if (mcount[7] == parms->movetype[7])
                {
                        ok=LocalMove(p,oldp,fragment,pot,parms->nmul_local,parms,t);
                        if(ok>-1) mcount[7]=0;
                        if(ok==1) macc[7]++;
                        mdone[7]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,7,istep,my_rank);
                }

		if( mcount[11] == parms->movetype[11])
		{
			ok=BackRub(p,oldp,pot,parms,t);
			if(ok>-1) mcount[11]=0;
			if(ok==1) macc[11]++;
			mdone[11]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,11,istep,my_rank);
		}
	    	if ( mcount[12] == parms->movetype[12] )
		{
			ok=BackSideRub(p,oldp,pot,parms,t);
			if(ok>-1) mcount[12]=0;
			if(ok==1) macc[12]++;
			mdone[12]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,12,istep,my_rank);

		}
		

		
		if (mcount[1] == parms->movetype[1])						// pivot of backbone
		{
			ok = MoveBackbonePivot(p,oldp,pot,parms,t);
			if (ok>-1) mcount[1] = 0;
			if (ok==1) macc[1]++;
			mdone[1]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,1,istep,my_rank);
		}
	
		if (mcount[2] == parms->movetype[2])						// multiple pivot of backbone
		{
			ok = MoveMultiplePivot(p,oldp,pot,parms->nmul_mpivot,parms,t);
			if (ok>-1) mcount[2] = 0;
			if (ok==1) macc[2]++;
			mdone[2]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,2,istep,my_rank);
		}
		
		if (mcount[6] == parms->movetype[6])						// move of center of mass
		{
			ok = MoveCoM(p,oldp,parms,pot,istep,parms->debug,t);
			if (ok>-1) mcount[6] = 0;
			if (ok==1) macc[6]++;
			mdone[6]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,6,istep,my_rank);
		}

		if (mcount[8] == parms->movetype[8])
		{
			ok=MoveRotation(p,oldp,parms,pot,istep,parms->debug,t);
			if(ok>-1) mcount[8]=0;
			if(ok==1) macc[8]++;
			mdone[8]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,8,istep,my_rank);
		}
				
		// if you want to anneal
		//if (parms->anneal) Anneal(parms,&t,&anneal_count,&anneal_status,&ok,&ishell,mcount);


		if (mcount[9] == parms->movetype[9])
		{
			ok=MoveClusterCoM(p,oldp,parms,pot,istep,parms->debug,t,fproc,my_rank);
			if(ok>-1) mcount[9]=0;
			if(ok==1) macc[9]++;
			mdone[9]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,9,istep,my_rank);
		}

		if (mcount[10] == parms->movetype[10])
		{
			ok=MoveClusterRot(p,oldp,parms,pot,istep,parms->npol,parms->debug,t,fproc,my_rank);
			if(ok>-1) mcount[10]=0;
			if(ok==1) macc[10]++;
			mdone[10]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,10,istep,my_rank);
		}

		if (mcount[13] == parms->movetype[13])						// aggregation-volume-bias move
		{
			ok=MoveAVB(p,oldp,parms,pot,t);
			if(ok>-1) mcount[13]=0;
			if(ok==1) macc[13]++;
			mdone[13]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,13,istep,my_rank);
		}

		if (mcount[14] == parms->movetype[14])						// configurational-bias regrowth of an end
		{
			ok=MoveCBMC(p,oldp,parms,pot,t);
			if(ok>-1) mcount[14]=0;
			if(ok==1) macc[14]++;
			mdone[14]++; if (parms->echeck) CheckEnergyDrift(p,oldp,pot,parms,14,istep,my_rank);
		}

/************************************************\
 * 	print di varie ed eventuali		*
\************************************************/

		// print log
		if (iprintlog == parms->nprintlog && ok>-2)
		{
			if(my_rank==0)		/* stampa solo per rank=0 */
			fprintf(parms->flog,"Step = %llu\tE=%lf\tok=%d\tT=%lf\n",istep,p->etot,ok,t);
			PrintEnergies_Parallel(stderr,parms->npol,istep,p,my_rank);
			//CountContacts(stderr,p,parms,parms->mov);
	
			#ifdef DEBUG
			fprintf(stderr,"%llu Etot=%lf\tEtot(true)=%lf ok=%d\n",istep,p->etot,
			TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank),ok);	// energy update in Move...
			
			if (parms->debug>4 && ok==1)
				PrintStructure(p,parms->npol,stderr,parms->shell);
			#endif
			iprintlog=0;
		}

		// print trajectory
		if (ftrj != NULL && iprinttrj == parms->nprinttrj && ok>-1)
		{
//DEBUG			#ifdef ACTIVE_MPI
//			fprintf(stderr,"printing trajectory MPI: rank %d, nback = %d\n",my_rank,(p)->nback);
//			#endif
			//CountContacts(stderr,p,parms,parms->mov);			

			//fprintf(stderr,"DEBUG: ok=%d\tTotal Energy %lf\n",ok,TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank));
			//fpirintf(stderr,"DEBUG: \t\t\tetot=%lf\n",p->etot);
			snprintf(p->title,sizeof(p->title),"step %llu\tE=%lf",istep,p->etot);
			PrintTrajectory(ftrj,p,parms,pot);
			fflush(ftrj);
			iprinttrj = 0;
			
			

			//fprintf(stderr,"Total Energy:\t%lf\n\n",TotalEnergy(p,pot,parms,parms->npol,0,parms->nosidechains,parms->debug,0));
		}

		// print the observables of aggregation
		if (parms->fobs && parms->nprintobs>0 && istep%parms->nprintobs==0)
			PrintObservables(parms->fobs,p,parms->npol,istep);

		// print energies
		if (fe != NULL && iprinte == parms->nprinte && ok>-1)
		{
			PrintEnergies(fe,parms->npol,istep,p);
			fflush(fe);
			iprinte = 0;
		}
		
		
//		fprintf(stderr,"ETOT\t%lf\n",p->etot);
		// update shell
		if (parms->shell==1)
		{//DEBUG SHELL
                        if( istep%parms->nshell==0)
			{
                          SetTrueShellFlags(p,parms);
			  UpdateShell(p,parms);
			  SetFalseShellFlags(p,parms);
                          CopyShell(p,oldp,parms);		
			  parms->ishell=0;		
			}
			if ( ( parms->ishell == 1 ) && ok>-1)
			{
				UpdateShell(p,parms);
	                        CopyShell(p,oldp,parms);
				parms->ishell=0;
			}
		}		
//ASTEMPERING 	
		#ifdef STEMPERING	
		if (parms->stempering)
		{
			st_o = STempering(p->etot,istep,(parms->p));

			if(!strcmp(parms->op_minim,"none"))
			{
				if ((parms->p)->st_ttarget_harvest==1 && st_iprint>=(parms->p)->st_printpdb && t <= (parms->p)->st_ttarget + EPSILON && (parms->p)->st_nm==1)
                {
					snprintf(p->title,sizeof(p->title),"step %llu\tE=%lf",istep,p->etot);
					PrintPDBStream(p,parms->npol,(parms->p)->st_pdbf);
					nconf++;

					if(parms->nconf!=-1)
						if(nconf>parms->nconf-1)
						{
							fprintf(stderr,"Reached target temperature (Ttarget=%lf).\nENDING SIMULATION\n",(parms->p)->st_ttarget);

							exit(0);
					
						}

					st_iprint = 0;
				}
				if ((parms->p)->st_ttarget_harvest>0 && ok>-1 && t <= (parms->p)->st_ttarget + EPSILON )
				st_iprint ++;
				
				if (st_o == 1) t = (parms->p)->st_temp[(parms->p)->st_itemp];
			}

		}
		#endif



		#ifdef OPTIMIZEPOT
        if(my_rank==0 || parms->op_allreplicas)				// op_allreplicas: every replica records its frames
		{
			if (strcmp(parms->op_minim,"none")  && ok>-1)
			{
				#ifdef STEMPERING
                if (parms->stempering)
				{
					if ((parms->p)->st_ttarget_harvest==1 && st_iprint>=(parms->p)->st_printpdb && t <= (parms->p)->st_ttarget + EPSILON && (parms->p)->st_nm==1)
					{

                                        
                                        //#ifdef ACTIVE_MPI

                                        //fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                        //       (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T[my_rank]);
                                        //#else
                                        //fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                        //       (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T);
                                        //#endif
                                        

                                        OP_RecordFrame(p,pot,parms,my_rank,t);

					 st_iprint = 0;

                                       
                                                        }                               

					if ((parms->p)->st_ttarget_harvest>0 && ok>-1 && t <= (parms->p)->st_ttarget + EPSILON ) st_iprint ++;

                                         if (st_o == 1) t = (parms->p)->st_temp[(parms->p)->st_itemp];

                                        }

                                        
                                else if (op->icount == parms->op_deltat)
                                
				#else
				if (op->icount == parms->op_deltat)
                                

				#endif
                                {        
                                        #ifdef ACTIVE_MPI
                                        
                                        if (my_rank==0) fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                               (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T[my_rank]);
                                        #else
                                        fprintf(stderr,"> recording polymer %d/%llu\tE=%lf\tT=%lf\n",op->nframes,
                                               (parms->nstep-parms->op_wait)/parms->op_deltat,p->etot,parms->T);
                                        #endif
                                        

                                        OP_RecordFrame(p,pot,parms,my_rank,t);
                                        op->icount=0;
                                }
                                if (istep >= parms->op_wait)
                                        op->icount ++;

                	}
                }

                #endif

                 #ifdef ACTIVE_MPI
                 if(ptempering_count==(parms->nstep_exchange)/2 && parms->ntemp>1)
                 {
			
                 	//Process syncing before exchange

                        MPI_Barrier(MPI_COMM_WORLD);
			
                        if(parms->debug>2) 
				if(my_rank==0)
					fprintf(stderr,"\nstep=%llu\n",istep);
                	//if(my_rank==0)fprintf(stderr,"~MPI: PT ...");
			ExchangePol(p,replica,oldp,parms,pot,my_rank,parms->ntemp,0,ex_count,ex_acc,Backtype,Sidetype,Rottype,astatus,istep,fexchange);	
			//if(my_rank==0)fprintf(stderr,"\n");
			fflush(fexchange);	

		}	
                
		if(ptempering_count==(parms->nstep_exchange) && parms->ntemp>1)
                {

		
                	MPI_Barrier(MPI_COMM_WORLD);
                        if(parms->debug>2) if(my_rank==0)
				fprintf(stderr,"\nstep=%llu\n",istep); 
                        ExchangePol(p,replica,oldp,parms,pot,my_rank,parms->ntemp,1,ex_count,ex_acc,Backtype,Sidetype,Rottype,astatus,istep,fexchange);
			
			ptempering_count=0;                        
                }

		// a soft exit (SIGUSR1) must stop all ranks at the same step, or the others hang in the next exchange
		if (parms->ntemp>1)
		{
			if (ptempering_count==1) { int se=softexit; MPI_Allreduce(&se,&stop,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD); }
		}
		else stop = softexit;
                
		ptempering_count++;
		#else
		stop = softexit;
		#endif

		
                




		// advance counters
                for (i=0;i<NMOVES;++i){ 
			mcount[i]++;
			if(mcount[i]>200000000) mcount[i]=0; //reset in order to prevent too big number
			}
                if (ok>-1)
                {
                	istep ++;
					iprinttrj ++;
                    iprintlog ++;
					iprinte ++;
					icheckpoint ++;
                 //       parms->ishell ++;
                }

		// rebuild the cell list from time to time (the shape of the system changes, see cells.c)
		if (p->cells && istep%10000==0) CellsBuild(p,pot->rcut,pot->boxtype=='p' ? pot->lbox : 0.);

		// checkpoint at the end of a step, so that a restart continues exactly from here
		// (by default on wall-clock time: writing a checkpoint does not change the trajectory);
		// a soft exit writes one before stopping
		if (stop || (parms->ncheckpoint>0 ? icheckpoint >= parms->ncheckpoint :
				(icheckpoint % 1000 == 0 && difftime(time(NULL),tcheckpoint) >= parms->tcheckpoint)))
		{
			tcheckpoint = time(NULL);
			icheckpoint = 0;
			SaveLoopState(chk,t,st_iprint,st_o,nconf,op->icount,ex_count,ex_acc,my_rank,parms->ntemp);
			SaveCheckpoint(chk,fchkp,p,pot,parms,irun,istep,iprinttrj,iprintlog,iprinte,mcount,macc,mdone,ptempering_count,ftrj,fe,fexchange,my_rank);
		}

		// progress in the file mginfo (every tmginfo seconds)
		if (my_rank==0 && istep%1000==0) MGInfo_MC(parms,irun,istep,p->etot,t,0);

	}
	///////////////////////////////////////// end MC LOOP
	if (my_rank==0) MGInfo_EndMC(parms,irun);
		
	// print summary
	#ifdef ACTIVE_MPI	
	fclose(fexchange);	
	fprintf(stderr,"\nRank %d Acceptance:\t%lld / %lld = %lf\n",my_rank,parms->acc,parms->mov,parms->mov ? (double)parms->acc/parms->mov : 0.);
	if (my_rank<parms->ntemp-1 && ex_count[my_rank]>0)
		fprintf(stderr,"Rank %d Exchanges with rank %d:\t%d / %d = %lf\n",my_rank,my_rank+1,ex_acc[my_rank],ex_count[my_rank],(double)ex_acc[my_rank]/ex_count[my_rank]);
	#else
	fprintf(stderr,"\nAcceptance:\t%lld / %lld = %lf\n",parms->acc,parms->mov,parms->mov ? (double)parms->acc/parms->mov : 0.);
	#endif
	{
		static const char *movename[NMOVES] = {"Flip","Pivot","MPivot","Sidechain","LPivot","MFlip","CoM shift",
							"Local","Rotation","Cluster CoM","Cluster Rot","BackRub","BackSideRub","AVB","CBMC"};
		for (i=0;i<NMOVES;++i)
			if (parms->movetype[i]!=-1)
				fprintf(stderr,"\t%s:\t%lld / %lld = %lf\n",movename[i],macc[i],mdone[i],mdone[i] ? (double)macc[i]/mdone[i] : 0.);
	}

	fprintf(stderr,"\n");
	if (stop) return 1;
	snprintf(p->title,sizeof(p->title),"final");
	PrintTrajectory(ftrj,p,parms,pot);

	return 0;
}

/*****************************************************************************
 Make a flip in a random backbone atom
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveBackboneFlip(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t)
{
	int iw,ip,ok=0,out,ra,rb;
	double dw,deltaE=0.;

	ip = irand(parms->npol);						// which chain to move
	iw = irand( (p+ip)->nback );						// which backbone to move

	#ifdef DEBUG
	 if (debug>2) fprintf(stderr,"step=%d flip ip = %d iw=%d\n",istep,ip,iw); fflush(stderr);
	#endif

	// the residues whose sidechain is built from iw (iw-1..iw+1, and the last oxygen of the chain from iw+2)
	ra = iw-1; rb = iw+1;
	SideRange(parms,ip,(p+ip)->nback,iw,iw,&ra,&rb);

	// old energies
	deltaE = - GetEnergyMonomerRange(p,ra,rb,ip);			// old energy of iw with the others (and of the residues whose sidechains move)

	if (!parms->noangpot)
	{
		if (iw>0) deltaE -= (((p+ip)->back)+iw-1)->e_ang ;		// old angular energy
		if (iw<(p+ip)->nback-1) deltaE -= (((p+ip)->back)+iw+1)->e_ang;
	}
	if (!parms->nodihpot)
	{
		if (iw>0) deltaE -= (((p+ip)->back)+iw-1)->e_dih;		// old dihedral energy
		deltaE -= (((p+ip)->back)+iw)->e_dih;
		if (iw<(p+ip)->nback-1) deltaE -= (((p+ip)->back)+iw+1)->e_dih;
		if (iw<(p+ip)->nback-2) deltaE -= (((p+ip)->back)+iw+2)->e_dih;
	}
	
	if (!parms->nohfields)
	{
		//We store ALL the H fields energy for the change of configuration
		deltaE -= GetHFields(p,parms->npol);
	}

	// make the move in p
	if (iw==0)
	{	
		ok = MoveHead((p+ip),parms);						// move the head
		ok *= AddSidechain(p,ra,rb,ip);
	}
	else if (iw==(p+ip)->nback-1)
	{
		ok = MoveTail((p+ip),parms);						// move the tail
		ok *= AddSidechain(p,ra,rb,ip);
	}
	else
	{
		dw = 2. * parms->dw_flip * (0.5 - frand());				// angle to flip
		ok = Flip((p+ip),iw,dw);						// make the flip
		ok *= AddSidechain(p,ra,rb,ip);
	}

	// check if constrains are violated

	if (parms->a_cloose>0)
	{
		if (iw-2>=0)
			if ( DAbs( Angle( (((p+ip)->back)+iw-2)->pos, (((p+ip)->back)+iw-1)->pos, (((p+ip)->back)+iw)->pos, (p+ip)->tables, &out)
			- (((p+ip)->back)+iw-1)->a_next) > parms->a_cloose || out==0) ok=0;
		if (ok && iw+2<(p+ip)->nback)
			if ( DAbs( Angle( (((p+ip)->back)+iw)->pos, (((p+ip)->back)+iw+1)->pos, (((p+ip)->back)+iw+2)->pos, (p+ip)->tables, &out)
					- (((p+ip)->back)+iw+1)->a_next) > parms->a_cloose || out==0) ok=0;
	}
	if (ok) ok = CheckLockedDihedralRange(p,ip,iw-1,iw+2,parms->d_cloose);			// a flip of iw changes the dihedrals iw-1..iw+2

	// calculate new energy
	if (ok)
	{
		deltaE += EnergyMonomerRange(p,pot,ra,rb,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);		// new energy of iw with the others
		if (!parms->noangpot)
		{
			if (iw>0) deltaE += EnergyAngles(p,pot,iw-1,ip,1);
			if (iw<(p+ip)->nback-1) deltaE += EnergyAngles(p,pot,iw+1,ip,1);
		}
		if (!parms->nodihpot)
		{
			if (iw>0) deltaE += EnergyDihedrals(p,pot,iw-1,ip,1);
			deltaE += EnergyDihedrals(p,pot,iw,ip,1);
			if (iw<(p+ip)->nback-1) deltaE += EnergyDihedrals(p,pot,iw+1,ip,1);
			if (iw<(p+ip)->nback-2) deltaE += EnergyDihedrals(p,pot,iw+2,ip,1);
		}
	
		if (!parms->nohfields)
		{
			deltaE += EnergyHFieldsAll(p,pot,parms);
		}
	}

	// Metropolis acceptance
	if (ok) ok = Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{
		UpdateMonomerRange(oldp,p,ra,rb,ip,parms->shell);							// return to the old position, contacts, etc.
	}
	else					// move accepted
	{
		UpdateMonomerRange(p,oldp,ra,rb,ip,parms->shell);							// update oldp
		p->etot += deltaE;
		parms->acc ++;
	}

	#ifdef DEBUG
	if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
	#endif
	parms->mov ++;

	if (ok == 0) return 0;
	else return 1;
}

/*****************************************************************************
 Copy an element of a polymer structure
 Watch the copy of the pointer to the rotamers --> do not free.
 *****************************************************************************/
static void CopyResidue(struct s_back *from, struct s_back *to, int contacts);

void CopyResiduePositions(struct s_back *from, struct s_back *to)
{
	CopyResidue(from,to,1);
}

/*****************************************************************************
 The same, with its contacts only if contacts=1
 *****************************************************************************/
static void CopyResidue(struct s_back *from, struct s_back *to, int contacts)
{
	int k;

	to->irot = from->irot;
	(to->pos).x = (from->pos).x;
	(to->pos).y = (from->pos).y;
	(to->pos).z = (from->pos).z;
	(to->sph).ang = (from->sph).ang;
	(to->sph).dih = (from->sph).dih;
	(to->sph).r = (from->sph).r;
	to->ext = from->ext;

	//sidechain
	for (k=0;k<from->nside;++k)
	{
		(((to->side)+k)->pos).x = (((from->side)+k)->pos).x;
		(((to->side)+k)->pos).y = (((from->side)+k)->pos).y;
		(((to->side)+k)->pos).z = (((from->side)+k)->pos).z;
	}

	// contacts
	if (contacts)
	{
		for (k=0;k<from->ncontacts;++k)
		{
			(to->contacts)[k] = (from->contacts)[k];
			(to->contacts_p)[k] = (from->contacts_p)[k];
			(to->contacts_r)[k] = (from->contacts_r)[k];
			(to->e)[k] = (from->e)[k];
		}
		to->ncontacts = from->ncontacts;
	}

	//shell
	for (k=0;k<from->nshell;++k)
	{
	
		(to->shell)[k] = (from->shell)[k];
		(to->shell_p)[k] = (from->shell_p)[k];
	}
	to->nshell = from->nshell;

}



void CopyResiduePositions_NOCONT(struct s_back *from, struct s_back *to)
{
        int k;

        to->irot = from->irot;
        (to->pos).x = (from->pos).x;
        (to->pos).y = (from->pos).y;
        (to->pos).z = (from->pos).z;
        (to->sph).ang = (from->sph).ang;
        (to->sph).dih = (from->sph).dih;
        (to->sph).r = (from->sph).r;
	to->ext = from->ext;
	
	 for (k=0;k<from->nside;++k)
        {
                (((to->side)+k)->pos).x = (((from->side)+k)->pos).x;
                (((to->side)+k)->pos).y = (((from->side)+k)->pos).y;
                (((to->side)+k)->pos).z = (((from->side)+k)->pos).z;
        }


	for (k=0;k<to->nshell;++k)
        {
                (to->shell)[k] = (from->shell)[k];
                (to->shell_p)[k] = (from->shell_p)[k];
        }
        to->nshell = from->nshell;

}

/*****************************************************************************
 Orthonormal frame of three positions (origin a, first axis along b-a)
 *****************************************************************************/
static void RigidFrame(struct vector a, struct vector b, struct vector c, double f[3][3])
{
	double n,d;
	int k;

	f[0][0] = b.x-a.x; f[0][1] = b.y-a.y; f[0][2] = b.z-a.z;
	n = sqrt(f[0][0]*f[0][0]+f[0][1]*f[0][1]+f[0][2]*f[0][2]);
	for (k=0;k<3;++k) f[0][k] /= n;
	f[1][0] = c.x-a.x; f[1][1] = c.y-a.y; f[1][2] = c.z-a.z;
	d = f[1][0]*f[0][0]+f[1][1]*f[0][1]+f[1][2]*f[0][2];
	for (k=0;k<3;++k) f[1][k] -= d*f[0][k];
	n = sqrt(f[1][0]*f[1][0]+f[1][1]*f[1][1]+f[1][2]*f[1][2]);
	for (k=0;k<3;++k) f[1][k] /= n;
	f[2][0] = f[0][1]*f[1][2]-f[0][2]*f[1][1];
	f[2][1] = f[0][2]*f[1][0]-f[0][0]*f[1][2];
	f[2][2] = f[0][0]*f[1][1]-f[0][1]*f[1][0];
}

// the position x (before the move) carried by the rigid motion of the frames fo (before) and fn (after)
static int RigidSame(struct vector x, struct vector ao, double fo[3][3], struct vector an, double fn[3][3], struct vector y)
{
	double c[3],z[3];
	int k;

	for (k=0;k<3;++k) c[k] = fo[k][0]*(x.x-ao.x) + fo[k][1]*(x.y-ao.y) + fo[k][2]*(x.z-ao.z);
	for (k=0;k<3;++k) z[k] = c[0]*fn[0][k] + c[1]*fn[1][k] + c[2]*fn[2][k];
	return (an.x+z[0]-y.x)*(an.x+z[0]-y.x) + (an.y+z[1]-y.y)*(an.y+z[1]-y.y) + (an.z+z[2]-y.z)*(an.z+z[2]-y.z) < RIGID_TOL2;
}

// widen the bounding box lo,hi to contain x
static void RigidGrow(struct vector x, double lo[3], double hi[3])
{
	double v[3]={x.x,x.y,x.z};
	int k;

	for (k=0;k<3;++k)
	{
		if (v[k]<lo[k]) lo[k] = v[k];
		if (v[k]>hi[k]) hi[k] = v[k];
	}
}

/*****************************************************************************
 After a pivot move, the residues [from,to] of chain ip were rotated
 rigidly with their backbone (from and to farther than parms->rigid from
 the pivot): flag those whose sidechain atoms (and the atoms their hydrogen
 bonds refer to) moved rigidly too, comparing with oldp. A sidechain that
 did not match its rotamer (a move rebuilt only the residues next to the
 atom it moved) was rebuilt by AddSidechain and is not flagged. Returns
 the flags of the residues of the chain, or NULL if none is used.
 *****************************************************************************/
static char *RigidFlags(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, struct s_mc_parms *parms, int ip, int from, int to)
{
	static char *rig=NULL;
	static int nrig=0;
	struct s_back *b,*bo;
	struct vector ao,an;
	double fo[3][3],fn[3][3];
	int i,is,ia,n=(p+ip)->nback;

	if (parms->rigid<0) return NULL;
	if (from<0) from = 0;
	if (to>n-1) to = n-1;
	if (to-from < 2) return NULL;
	if (n > nrig)
	{
		rig = (char *) realloc(rig,n);
		if (!rig) Error("Cannot allocate in RigidFlags");
		nrig = n;
	}
	memset(rig,0,n);

	b = (p+ip)->back;
	bo = (oldp+ip)->back;
	ao = bo[from].pos;
	an = b[from].pos;
	RigidFrame(bo[from].pos,bo[from+1].pos,bo[from+2].pos,fo);
	RigidFrame(b[from].pos,b[from+1].pos,b[from+2].pos,fn);

	for (i=from;i<=to;++i)
	{
		rig[i] = 1;
		if (!parms->nosidechains)
			for (is=0;is<b[i].nside && rig[i];++is)
				rig[i] = RigidSame(bo[i].side[is].pos,ao,fo,an,fn,b[i].side[is].pos);
		if (parms->hb && rig[i])
		{
			if (pot->hb[b[i].ia]) rig[i] = RigidSame(*((oldp+ip)->vback[pot->hb_iam[b[i].ia]]),ao,fo,an,fn,*((p+ip)->vback[pot->hb_iam[b[i].ia]]));
			for (is=0;is<b[i].nside && rig[i];++is)
			{
				ia = b[i].side[is].ia;
				if (pot->hb[ia]) rig[i] = RigidSame(*((oldp+ip)->vback[pot->hb_iam[ia]]),ao,fo,an,fn,*((p+ip)->vback[pot->hb_iam[ia]]));
			}
		}
	}

	// in a periodic box the image of a pair can change with the rotation, unless
	// the rigid atoms fit in a box whose diagonal is shorter than half a side
	if (pot->boxtype=='p')
	{
		double lo[3]={1e30,1e30,1e30},hi[3]={-1e30,-1e30,-1e30},d2=0;
		int k;

		for (i=from;i<=to;++i)
			if (rig[i])
			{
				RigidGrow(bo[i].pos,lo,hi);
				if (!parms->nosidechains)
					for (is=0;is<bo[i].nside;++is) RigidGrow(bo[i].side[is].pos,lo,hi);
				if (parms->hb)
				{
					if (pot->hb[b[i].ia]) RigidGrow(*((oldp+ip)->vback[pot->hb_iam[b[i].ia]]),lo,hi);
					for (is=0;is<b[i].nside;++is)
					{
						ia = b[i].side[is].ia;
						if (pot->hb[ia]) RigidGrow(*((oldp+ip)->vback[pot->hb_iam[ia]]),lo,hi);
					}
				}
			}
		for (k=0;k<3;++k) if (hi[k]>lo[k]) d2 += (hi[k]-lo[k])*(hi[k]-lo[k]);
		if (4.*d2 >= pot->lbox*pot->lbox) return NULL;
	}
	return rig;
}

/*****************************************************************************
 Make a pivot move in a random backbone atom (very inefficient!)
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveBackbonePivot(struct s_polymer *p,struct s_polymer *oldp, struct s_potential *pot, struct s_mc_parms *parms, double t)
{
	int iw,ok=0,half,ip,ra,rb;
	double dw=0,deltaE;
	char *rig;
		
	ip = irand(parms->npol);

													// which chain to move
	half = (p+ip)->nback / 2;
	iw = 1 + irand((p+ip)->nback-2);

	if (parms->randdw==1) dw = parms->dw_pivot * (0.5 - frand());						// dihedral to pivot
	else if (parms->randdw==2) dw = parms->dw_pivot * gauss(0.,1.);

	#ifdef DEBUG
		if (parms->debug>2) fprintf(stderr,"pivot ip = %d iw=%d dw=%lf half=%d\n",ip,iw,dw,half);
	#endif

	
	// if iw belongs to the first half, pivot backward (moving dih of iw)
	if (iw < half)
	{

		if ( (((p+ip)->back)+iw+1)->move == 0 ) { parms->mov++; return 0;}

		// calculate old energy (the two-body one after the move)
		deltaE = 0;
		if (!parms->nodihpot)
		{
			deltaE -= (((p+ip)->back)+iw+1)->e_dih;							// old dihedral energy
		}
		if (!parms->nohfields)
		{
			deltaE -= GetHFields(p,parms->npol);										// old H fields energy
		}
		
		// move
		ok = PivotBackward((p+ip),iw+1,dw,iw,parms);							// moved atoms are in [0,iw-1]; changed dihedral is iw+1
		ra = 0; rb = iw;									// and the sidechains built from [0,iw-1]
		SideRange(parms,ip,(p+ip)->nback,0,iw-1,&ra,&rb);
		ok *= AddSidechain(p,ra,rb,ip);

		// the pairs of residues that moved rigidly together do not change: old two-body energy without them (from oldp)
		rig = RigidFlags(p,oldp,pot,parms,ip,0,iw-1-parms->rigid);
		deltaE -= GetEnergyMonomerRangeRigid(oldp,ra,rb,ip,rig);
		deltaE += EnergyMonomerRangeRigid(p,pot,ra,rb,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb,rig);			// new energy of iw with the others
		if (!parms->nodihpot)
			deltaE += EnergyDihedrals(p,pot,iw+1,ip,1);
		if (!parms->nohfields)
		{
			deltaE += EnergyHFieldsAll(p,pot,parms);
		}
	}
	// if it belong to the second half, pivot forward (moving dih of iw)
	else
	{
	
		if ( (((p+ip)->back)+iw)->move == 0 ) { parms->mov++; return 0; }

		// calculate old energy (the two-body one after the move)
		deltaE = 0;
		if (!parms->nodihpot)
		{
			deltaE -= (((p+ip)->back)+iw)->e_dih;							// old dihedral energy
		}
		if (!parms->nohfields)
		{
			deltaE -= GetHFields(p,parms->npol);										// old H fields energy
		}
		// move
		ok = PivotForward((p+ip),iw-1,dw,(p+ip)->nback-iw-1,parms);			// moved atoms are in [iw+1,nback-1]; changed dihedral of iw
		ra = iw; rb = (p+ip)->nback-1;							// and the sidechains built from [iw+1,nback-1]
		SideRange(parms,ip,(p+ip)->nback,iw+1,(p+ip)->nback-1,&ra,&rb);
		ok *= AddSidechain(p,ra,rb,ip);

		rig = RigidFlags(p,oldp,pot,parms,ip,iw+1+parms->rigid,(p+ip)->nback-1);
		deltaE -= GetEnergyMonomerRangeRigid(oldp,ra,rb,ip,rig);
		deltaE += EnergyMonomerRangeRigid(p,pot,ra,rb,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb,rig);
		if (!parms->nodihpot)
			deltaE += EnergyDihedrals(p,pot,iw,ip,1);
		if (!parms->nohfields)
		{
			deltaE += EnergyHFieldsAll(p,pot,parms);
		}
	}

	#ifdef DEBUG
	 int dir = (iw-half)<0 ? -1 : 1;
	 if (parms->debug>2) fprintf(stderr," ip=%d iw=%d dw=%lf direction=%d\n deltaE=%lf\n",ip,iw,dw,dir,deltaE); fflush(stderr);
	#endif

	ok *= Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{

		UpdateMonomerRange(oldp,p,ra,rb,ip,0);                               // return to the old position, contacts, etc.
		
	}
	else					// move accepted
	{
	        UpdateMonomerRange(p,oldp,ra,rb,ip,0);                               // update oldp

		if(parms->shell==1 )		// if shells are active, update them 
		{
//			UpdateShell(p,parms);
//			CopyShell(p,oldp,parms);
			parms->ishell=1;
			if (iw < half )
			  UpdateShellFlags(p,ip,0,iw);
			else
                          UpdateShellFlags(p,ip,iw,(p+ip)->nback-1);
		}
		p->etot += deltaE;
		parms->acc ++;
	}

	#ifdef DEBUG
	if (parms->debug>2) fprintf(stderr," accept=%d\n",ok); fflush(stderr);
	#endif
	parms->mov++;

	return ok;
}

/*****************************************************************************
 Make nmul consecutive pivot moves at random
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveMultiplePivot(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, int nmul, struct s_mc_parms *parms, double t)
{

	int half,ok=1,m,iw,ip,idir,ra,rb;
	double dw=0,deltaE;
	char *rig;

	if(nmul<0) nmul=2+irand(-nmul-1);				//generate random nmul in [2,-mul] if required

	ip = irand(parms->npol);
	half = (p+ip)->nback / 2;
	iw = 1 + irand((p+ip)->nback-2);

	if (iw<half-nmul-1) idir = -1;				// which direction to go
	else if (iw>half+nmul+1) idir=1;
	else idir = 2 * irand(2) - 1;

	#ifdef DEBUG
		if (parms->debug>2) fprintf(stderr,"mpivot ip = %d iw=%d half=%d\n",ip,iw,half);
	#endif

	// if iw belongs to the first half, pivot backward
	if (idir == -1)
	{
		if ( iw+2-nmul < 2 ) nmul = iw;										// if nmul beyond beginning of the chain, decrease it

		deltaE = 0;											// old energy (the two-body one after the move)
		if (!parms->nodihpot)
			for (m=0;m<nmul;++m) deltaE -= (((p+ip)->back)+iw-m+1)->e_dih;	// old dihedral energy
		if (!parms->nohfields){
			deltaE -= GetHFields(p,parms->npol);										// old H fields energy
		}
			
		
		for (m=0;m<nmul;++m)												// apply nmul pivot moves
			if ( (((p+ip)->back)+iw-m+1)->move == 1 )						// check if it you can move it
			{
				if (parms->randdw==1) dw = parms->dw_mpivot * (0.5 - frand());				// dihedral to pivot
				else if (parms->randdw==2) dw = parms->dw_mpivot *gauss(0.,1.);

 				ok *= PivotBackward((p+ip),iw+1-m,dw,iw-m,parms);				// moves in [0,iw-1]
			}
		ra = 0; rb = iw;									// and the sidechains built from [0,iw-1]
		SideRange(parms,ip,(p+ip)->nback,0,iw-1,&ra,&rb);
		ok *= AddSidechain(p,ra,rb,ip);

//		deltaE += EnergyMonomerRange(p,pot,0,iw,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);				// new energy of iw with the others
		rig = RigidFlags(p,oldp,pot,parms,ip,0,iw-nmul-parms->rigid);					// [0,iw-nmul] gets all the rotations
		deltaE -= GetEnergyMonomerRangeRigid(oldp,ra,rb,ip,rig);
		deltaE += EnergyMonomerRangeRigid(p,pot,ra,rb,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb,rig);    


		if (!parms->nodihpot)													// new dihedral energy
			for (m=0;m<nmul;++m) deltaE += EnergyDihedrals(p,pot,iw-m+1,ip,1);
		if (!parms->nohfields)
		{
			deltaE += EnergyHFieldsAll(p,pot,parms);
		}
	}
	// if it belong to the second half, pivot forward
	else
	{
		if ( iw-2+nmul > (p+ip)->nback - 3 ) nmul = (p+ip)->nback - 1 - iw;			// if nmul is beyond the end, decrease it

		deltaE = 0;											// old energy (the two-body one after the move)
		if (!parms->nodihpot)
			for (m=0;m<nmul;++m)deltaE -= (((p+ip)->back)+iw+m)->e_dih;				// old dihedral energy
		if (!parms->nohfields)
			deltaE -= GetHFields(p,parms->npol);												// old H fields energy

		for (m=0;m<nmul;++m)														// apply nmul pivot moves
			if ( (((p+ip)->back)+iw+m)->move == 1 )								// check if it you can move it
			{
				if (parms->randdw==1) dw = parms->dw_mpivot * (0.5 - frand());					// dihedral to pivot
				else if (parms->randdw==2) dw = parms->dw_mpivot * gauss(0.,1.);

				ok *= PivotForward((p+ip),iw-1+m,dw,(p+ip)->nback-iw-1-m,parms);	// moves in [iw+1,nback-1]
			}
		ra = iw; rb = (p+ip)->nback-1;							// and the sidechains built from [iw+1,nback-1]
		SideRange(parms,ip,(p+ip)->nback,iw+1,(p+ip)->nback-1,&ra,&rb);
		ok *= AddSidechain(p,ra,rb,ip);

		//shell disattivate: mossa globale
		rig = RigidFlags(p,oldp,pot,parms,ip,iw+nmul+parms->rigid,(p+ip)->nback-1);			// [iw+nmul,nback-1] gets all the rotations
		deltaE -= GetEnergyMonomerRangeRigid(oldp,ra,rb,ip,rig);
		deltaE += EnergyMonomerRangeRigid(p,pot,ra,rb,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb,rig);   
	//	deltaE += EnergyMonomerRange(p,pot,iw,p->nback-1,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);			// new energy of iw with the others
		if (!parms->nodihpot)														// new dihedral energy
			for (m=0;m<nmul;++m) deltaE += EnergyDihedrals(p,pot,iw+m,ip,1);
		if (!parms->nohfields)
		{
			deltaE += EnergyHFieldsAll(p,pot,parms);
		}
	}

//	fprintf(stderr,"\ndeltaE=%lf\n",deltaE);
	#ifdef DEBUG
	if (parms->debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
	#endif

	ok *= Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{	

	//	if (idir == -1) UpdateMonomerRange(oldp,p,0,iw,ip,parms->shell);				// return to the old position, contacts, etc. (if moved first half)
		UpdateMonomerRange(oldp,p,ra,rb,ip,0);


			
	}
	else					// move accepted
	{
//		if (idir == -1) UpdateMonomerRange(p,oldp,0,iw,ip,parms->shell);				// update oldp (if moved first half)
		UpdateMonomerRange(p,oldp,ra,rb,ip,0);



		if(parms->shell==1 )            // if shells are active, update them 
                {
                    //    UpdateShell(p,parms);
                      //  CopyShell(p,oldp,parms);
                        parms->ishell=1;
                }
 
		p->etot += deltaE;
		parms->acc ++;
	}

	#ifdef DEBUG
	if (parms->debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
	#endif

	parms->mov++;

	return ok;
}



/*****************************************************************************
 Metropolis acceptance rule
 1 = accept, 0 = reject
 *****************************************************************************/
int Metropolis(double deltaE, double T, struct s_tables *t)
{
	// the random number is drawn also when deltaE<=0, so that the stream of random
	// numbers does not depend on the sign of a deltaE which is zero up to rounding
	double r = frand();

	if (deltaE<=0)
		return 1;
	return ( r < FastExp(-deltaE/T,t) );
}


/*****************************************************************************
 Copy the position of the monomer w to another structure, updating also
 contacts, partial energies and shell
 The other structure should be identical to the former except for residue w
 *****************************************************************************/
void UpdateMonomer(struct s_polymer *from, struct s_polymer *to, int w, int n, int shell)
{
	UpdateChainRanges(from,to,1,&n,&w,&w);
}

void UpdateMonomerRange(struct s_polymer *from, struct s_polymer *to, int wfrom, int wto, int p, int shell)
{
	UpdateChainRanges(from,to,1,&p,&wfrom,&wto);
}

// residue i of chain c is in one of the ranges of UpdateChainRanges
static inline int InRanges(int c, int i, int nr, int *ch, int *wfrom, int *wto)
{
	int r;

	for (r=0;r<nr;++r)
		if (ch[r]==c && i>=wfrom[r] && i<=wto[r]) return 1;
	return 0;
}

/*****************************************************************************
 Copy residues [wfrom[r],wto[r]] of chains ch[r] (r=0..nr-1), and the contact
 lists that the move changed (in p: see ContactsTouched), from "from" to
 "to". The lists it did not change are the same in both.
 *****************************************************************************/
void UpdateChainRanges(struct s_polymer *from, struct s_polymer *to, int nr, int *ch, int *wfrom, int *wto)
{
	int r,i,w,j,p,wf,wt,nt,*tc,*ti;
	struct s_back *bf,*bt;

	for (r=0;r<nr;++r)
	{
		p = ch[r];
		wf = (wfrom[r]<0) ? 0 : wfrom[r];
		wt = (wto[r]>(from+p)->nback-1) ? (from+p)->nback-1 : wto[r];

		// Update angle/dihedral energies from w-1 to w+2
		for (w=wf;w<=wt;++w)
		{
			(((to+p)->back)+w)->e_ang = (((from+p)->back)+w)->e_ang;
			if (w>0) (((to+p)->back)+w-1)->e_ang = (((from+p)->back)+w-1)->e_ang;
			if (w<(from+p)->nback-1) (((to+p)->back)+w+1)->e_ang = (((from+p)->back)+w+1)->e_ang;

			if (w>0) (((to+p)->back)+w-1)->e_dih = (((from+p)->back)+w-1)->e_dih;
			(((to+p)->back)+w)->e_dih = (((from+p)->back)+w)->e_dih;
			if (w<((from+p)->nback)-1) (((to+p)->back)+w+1)->e_dih = (((from+p)->back)+w+1)->e_dih;
			if (w<((from+p)->nback)-2) (((to+p)->back)+w+2)->e_dih = (((from+p)->back)+w+2)->e_dih;

			(((to+p)->back)+w)->irot = (((from+p)->back)+w)->irot;
			(((to+p)->back)+w)->e_hfs = (((from+p)->back)+w)->e_hfs;
		}
	}

	// Copy the contacts changed by the move (the moved residues are copied below)
	nt = ContactsTouched(&tc,&ti);
	for (i=0;i<nt;++i)
	{
		if (InRanges(tc[i],ti[i],nr,ch,wfrom,wto)) continue;
		bf = ((from+tc[i])->back)+ti[i];
		bt = ((to+tc[i])->back)+ti[i];
		bt->e_hfs = bf->e_hfs;							// H fields depend on the number of contacts
		bt->ncontacts = bf->ncontacts;
		for (j=0;j<bf->ncontacts;++j)
		{
			bt->contacts[j] = bf->contacts[j];
			bt->contacts_p[j] = bf->contacts_p[j];
			bt->contacts_r[j] = bf->contacts_r[j];
			bt->e[j] = bf->e[j];
		}
	}

	// Finally, copy the moved residues themselves
	for (r=0;r<nr;++r)
	{
		p = ch[r];
		wf = (wfrom[r]<0) ? 0 : wfrom[r];
		wt = (wto[r]>(from+p)->nback-1) ? (from+p)->nback-1 : wto[r];
		for (w=wf;w<=wt;++w)
			CopyResidue( ((from+p)->back)+w, ((to+p)->back)+w,
				(((from+p)->back)+w)->touched == contacts_stamp || (((to+p)->back)+w)->touched == contacts_stamp );
		if (to->cells) CellsUpdateRange(to,p,wf,wt);
	}

	// the next move starts a new stamp
	ContactsNewMove();
}






/*****************************************************************************
 Move the sidechain among the possible rotamers
	 returns 1 if it could make the move, 0 if nothing was done,
	 -1 if it couldn't find a sidechain to move
 *****************************************************************************/
int MoveSidechain(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *mc_parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t)
{
	int iw,ok=0,chk=0,ir,ip;
	double deltaE;


	// which chain to move
	ip = irand(mc_parms->npol);

	// select among sidechains which have more than 1 rotamer
	do
	{
		iw = irand((p+ip)->nback);
		++chk;
		if (chk>(p+ip)->nback*3) {mc_parms->mov++; return -1;}
	}
	while ( (((p+ip)->back)+iw)->nside < 1 || (((p+ip)->back)+iw)->nrot ==1 );


	#ifdef DEBUG
	 if (debug>2) fprintf(stderr,"step=%d sidechain iw=%d\n",istep,iw); fflush(stderr);
	#endif

	deltaE = - GetEnergyMonomer(p,ip,iw);								// old energy of iw with the others

	if (!mc_parms->nohfields){
		deltaE -= GetHFields(p,mc_parms->npol);										// old H fields energy
	}
	
	ir = irand( (((p+ip)->back)+iw)->nrot );
	(((p+ip)->back)+iw)->irot = ir;
	ok = AddSidechain(p,iw,iw,ip);

	deltaE += EnergyMonomer(p,pot,iw,ip,mc_parms->npol,1,mc_parms->shell,mc_parms->nosidechains,mc_parms->disentangle,mc_parms->hb);			// new energy of iw with the others
	//deltaE += EnergyMonomer(p,pot,iw,ip,mc_parms->npol,1,0,mc_parms->nosidechains,mc_parms->disentangle,mc_parms->hb); 
	
	if (!mc_parms->nohfields)
	{
		deltaE += EnergyHFieldsAll(p,pot,mc_parms);
	}
	
	#ifdef DEBUG
	if (debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
	#endif

	ok *= Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{
		// return to the old position, contacts, etc.
		UpdateMonomer(oldp,p,iw,ip,mc_parms->shell);
		//UpdateMonomer(oldp,p,iw,ip,0);
	}
	else					// move accepted
	{
		UpdateMonomer(p,oldp,iw,ip,mc_parms->shell);
		//update oldp
		//UpdateMonomer(p,oldp,iw,ip,0);
		p->etot += deltaE;
		mc_parms->acc ++;

	}

	#ifdef DEBUG
	if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
	#endif

	mc_parms->mov ++;

	if (ok == 0) return 0;
	else return 1;
}

void CompareStructures(struct s_polymer *a, struct s_polymer *b, int nc, int natoms)
{
	int ic,i,j,irot;


	for (ic=0;ic<nc;++ic)
		for (i=0;i<(a+ic)->nback;++i)
		{
			irot = (((a+ic)->back)+i)->irot;

			if ( (((a+ic)->back)+i)->ia != (((b+ic)->back)+i)->ia ) fprintf(stderr,"*** i=%d ic=%d ia=%d %d\n",i,ic,(((a+ic)->back)+i)->ia,(((b+ic)->back)+i)->ia);
			if ( (((a+ic)->back)+i)->iaa != (((b+ic)->back)+i)->iaa ) fprintf(stderr,"*** i=%d ic=%d iaa=%d %d\n",i,ic,(((a+ic)->back)+i)->iaa,(((b+ic)->back)+i)->iaa);
			if ( (((a+ic)->back)+i)->itype != (((b+ic)->back)+i)->itype ) fprintf(stderr,"*** i=%d ic=%d itype=%d %d\n",i,ic,(((a+ic)->back)+i)->itype,(((b+ic)->back)+i)->itype);
			if ( (((a+ic)->back)+i)->nside != (((b+ic)->back)+i)->nside ) fprintf(stderr,"*** i=%d ic=%d nside=%d %d\n",i,ic,(((a+ic)->back)+i)->nside,(((b+ic)->back)+i)->nside);
			if ( (((a+ic)->back)+i)->move != (((b+ic)->back)+i)->move ) fprintf(stderr,"*** i=%d ic=%d move=%d %d\n",i,ic,(((a+ic)->back)+i)->move,(((b+ic)->back)+i)->move);
			if ( (((a+ic)->back)+i)->irot != (((b+ic)->back)+i)->irot ) fprintf(stderr,"*** i=%d ic=%d irot=%d %d\n",i,ic,(((a+ic)->back)+i)->irot,(((b+ic)->back)+i)->irot);
			if ( (((a+ic)->back)+i)->nrot != (((b+ic)->back)+i)->nrot ) fprintf(stderr,"*** i=%d ic=%d nrot=%d %d\n",i,ic,(((a+ic)->back)+i)->nrot,(((b+ic)->back)+i)->nrot);
			if ( (((a+ic)->back)+i)->ncontacts != (((b+ic)->back)+i)->ncontacts ) fprintf(stderr,"*** i=%d ic=%d ncotacts=%d %d\n",i,ic,(((a+ic)->back)+i)->ncontacts,(((b+ic)->back)+i)->ncontacts);
			if ( (((a+ic)->back)+i)->pos.x != (((b+ic)->back)+i)->pos.x ) fprintf(stderr,"*** i=%d ic=%d pos.x=%lf %lf\n",i,ic,(((a+ic)->back)+i)->pos.x,(((b+ic)->back)+i)->pos.x);
			if ( (((a+ic)->back)+i)->pos.y != (((b+ic)->back)+i)->pos.y ) fprintf(stderr,"*** i=%d ic=%d pos.y=%lf %lf\n",i,ic,(((a+ic)->back)+i)->pos.y,(((b+ic)->back)+i)->pos.y);
			if ( (((a+ic)->back)+i)->pos.z != (((b+ic)->back)+i)->pos.z ) fprintf(stderr,"*** i=%d ic=%d pos.z=%lf %lf\n",i,ic,(((a+ic)->back)+i)->pos.z,(((b+ic)->back)+i)->pos.z);
			if ( (((a+ic)->back)+i)->sph.ang!= (((b+ic)->back)+i)->sph.ang ) fprintf(stderr,"*** i=%d ic=%d ang=%lf %lf\n",i,ic,(((a+ic)->back)+i)->sph.ang,(((b+ic)->back)+i)->sph.ang);
			if ( (((a+ic)->back)+i)->sph.dih!= (((b+ic)->back)+i)->sph.dih ) fprintf(stderr,"*** i=%d ic=%d dih=%lf %lf\n",i,ic,(((a+ic)->back)+i)->sph.dih,(((b+ic)->back)+i)->sph.dih);
			if ( (((a+ic)->back)+i)->sph.r!= (((b+ic)->back)+i)->sph.r ) fprintf(stderr,"*** i=%d ic=%d r=%lf %lf\n",i,ic,(((a+ic)->back)+i)->sph.r,(((b+ic)->back)+i)->sph.r);
			for (j=0;j<(((a+ic)->back)+i)->nside;++j)
			{
				if ( (((((a+ic)->back)+i)->side)+j)->ia != (((((b+ic)->back)+i)->side)+j)->ia ) fprintf(stderr,"*** i=%d ic=%d side %d ia=%d %d\n",i,ic,j,
					(((((a+ic)->back)+i)->side)+j)->ia,(((((b+ic)->back)+i)->side)+j)->ia);
				if ( (((((a+ic)->back)+i)->side)+j)->pos.x != (((((b+ic)->back)+i)->side)+j)->pos.x ) fprintf(stderr,"*** i=%d ic=%d side %d x=%lf %lf\n",i,ic,j,
					(((((a+ic)->back)+i)->side)+j)->pos.x,(((((b+ic)->back)+i)->side)+j)->pos.x);
				if ( (((((a+ic)->back)+i)->side)+j)->pos.y != (((((b+ic)->back)+i)->side)+j)->pos.y ) fprintf(stderr,"*** i=%d ic=%d side %d y=%lf %lf\n",i,ic,j,
					(((((a+ic)->back)+i)->side)+j)->pos.y,(((((b+ic)->back)+i)->side)+j)->pos.y);
				if ( (((((a+ic)->back)+i)->side)+j)->pos.z != (((((b+ic)->back)+i)->side)+j)->pos.z ) fprintf(stderr,"*** i=%d ic=%d side %d z=%lf %lf\n",i,ic,j,
					(((((a+ic)->back)+i)->side)+j)->pos.z,(((((b+ic)->back)+i)->side)+j)->pos.z);
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b1 != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b1 ) fprintf(stderr,"*** i=%d ic=%d side %d b1=%d %d\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b1,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b1 );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b2 != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b2 ) fprintf(stderr,"*** i=%d ic=%d side %d b2=%d %d\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b2,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b2 );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b3 != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b3 ) fprintf(stderr,"*** i=%d ic=%d side %d b3=%d %d\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->b3,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->b3 );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.ang != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.ang ) fprintf(stderr,"*** i=%d ic=%d side %d ang=%lf %lf\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.ang,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.ang );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.dih != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.dih ) fprintf(stderr,"*** i=%d ic=%d side %d dih=%lf %lf\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.dih,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.dih );
				if ( (((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.r != (((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.r ) fprintf(stderr,"*** i=%d ic=%d side %d r=%lf %lf\n",
					i,ic,j,(((((((a+ic)->back)+i)->side)+j)->rot)+irot)->ang.r,(((((((b+ic)->back)+i)->side)+j)->rot)+irot)->ang.r );



			}

		}

	for (ic=0;ic<nc;++ic)
		for (i=0;i<natoms;++i)
				{
					if ( (*((a+ic)->vback[i])).x != (*((b+ic)->vback[i])).x ) fprintf(stderr,"*** ia=%d x=%lf %lf\n",i,(*((a+ic)->vback[i])).x,(*((b+ic)->vback[i])).x);
					if ( (*((a+ic)->vback[i])).y != (*((b+ic)->vback[i])).y ) fprintf(stderr,"*** ia=%d y=%lf %lf\n",i,(*((a+ic)->vback[i])).y,(*((b+ic)->vback[i])).y);
					if ( (*((a+ic)->vback[i])).z != (*((b+ic)->vback[i])).z ) fprintf(stderr,"*** ia=%d z=%lf %lf\n",i,(*((a+ic)->vback[i])).z,(*((b+ic)->vback[i])).z);
				}
}



/*****************************************************************************
 Make nmul consecutive pivot moves at random and keep fixed  the nmul+1,
 chaniging the bond length
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/

int MoveLoosePivot(struct s_polymer *p, struct s_polymer *oldp, struct s_potential *pot, int nmul, struct s_mc_parms *parms, double t)
{
      int half,ok=1,m,iw,ip,idir,ra,rb;
      double dw=0,deltaE,rc2;

      if(nmul<0) nmul=2+irand(-nmul-1);                     //generate random nmul in [2,-mul] if required

      rc2 = parms->r_cloose * parms->r_cloose;

      ip = irand(parms->npol);
      half = (p+ip)->nback / 2;
      iw = 1 + irand((p+ip)->nback-2);

      if (iw<half-nmul-1) idir = -1;                        // which direction to go
      else if (iw>half+nmul+1) idir=1;
      else idir = 2 * irand(2) - 1;

      #ifdef DEBUG
            if (parms->debug>2) fprintf(stderr,"lpivot ip = %d iw=%d half=%d  nmul=%d\n",ip,iw,half,nmul);
      #endif

      if (idir == -1)
      {
            if ( iw-nmul-1 < 2 ) nmul = iw - 1;                                                       // if nmul beyond beginning of the chain, decrease it
            if (nmul < 1) { parms->mov ++; return 0; }                                                // nothing to move

            ra = iw-nmul-1; rb = iw;									// and the sidechains built from [iw-nmul,iw-1]
            SideRange(parms,ip,(p+ip)->nback,iw-nmul,iw-1,&ra,&rb);
            deltaE = -GetEnergyMonomerRange(p,ra,rb,ip);                                 // calculate old energy
            if (!parms->nodihpot)
				for (m=0;m<nmul+3;++m) deltaE -= (((p+ip)->back)+iw-m+1)->e_dih;  // old dihedral energy in [iw-nmul-1,iw+1]
            if (!parms->noangpot)
            {
				deltaE -= (((p+ip)->back)+iw-nmul-1)->e_ang ;                                 // old angular energy
				deltaE -= (((p+ip)->back)+iw-nmul)->e_ang;
			}
			if (!parms->nohfields)
				deltaE -= GetHFields(p,parms->npol);													  // old H fields energy
		  
		  
            for (m=0;m<nmul;++m)                                                                            // apply nmul pivot moves
                  if ( (((p+ip)->back)+iw+1-m)->move == 1 )                                           // PivotBackward(k) changes dihedral k
                  {
                        if (parms->randdw==1) dw = parms->dw_lpivot * (0.5 - frand());                            // dihedral to pivot
                        else if (parms->randdw==2) dw = parms->dw_lpivot * gauss(0.,1.);

                        ok *= PivotBackward((p+ip),iw+1-m,dw,nmul-m,parms);                     // moves iw-1-m => overall moves in [iw-nmul,iw-1] (and sidechains of iw-imul-1 and iw)
                  }
            if (iw-nmul-1>=0)
                  if ( DAbs(Dist2( (((p+ip)->back)+iw-nmul)->pos, (((p+ip)->back)+iw-nmul-1)->pos ) - (((p+ip)->back)+iw-nmul-1)->d2_next) > rc2 ) ok=0;           // if covalent bond is broken
            // the stretched bond changes the angles at iw-nmul-1, iw-nmul and the dihedrals iw-nmul-1..iw-nmul+1
            if (ok) ok = CheckAngleRange(p,ip,iw-nmul-1,iw-nmul,parms->a_cloose) && CheckLockedDihedralRange(p,ip,iw-nmul-1,iw-nmul+1,parms->d_cloose);

            if (ok==1)
            {
				ok = AddSidechain(p,ra,rb,ip);

				deltaE += EnergyMonomerRange(p,pot,ra,rb,ip,parms->npol,parms->shell,1,
                                                                        parms->nosidechains,parms->disentangle,parms->hb);                      // new energy of iw with the others
				if (!parms->nodihpot)                                                                                                         // new dihedral energy
                        for (m=0;m<nmul+3;++m) deltaE += EnergyDihedrals(p,pot,iw-m+1,ip,1);
				if (!parms->noangpot)                                                                                                         // new angle energy
				{
                        deltaE += EnergyAngles(p,pot,iw-nmul-1,ip,1);
                        deltaE += EnergyAngles(p,pot,iw-nmul,ip,1);
				}
				
				if (!parms->nohfields)
				{
					deltaE += EnergyHFieldsAll(p,pot,parms);
				}
            }
      }
      else
      {
            if ( iw+2+nmul > (p+ip)->nback - 1 ) nmul = (p+ip)->nback - 3 - iw;                // if nmul is beyond the end, decrease it
            if (nmul < 1) { parms->mov ++; return 0; }                                                // nothing to move

            ra = iw; rb = iw+nmul+1;									// and the sidechains built from [iw+1,iw+nmul]
            SideRange(parms,ip,(p+ip)->nback,iw+1,iw+nmul,&ra,&rb);
            deltaE = -GetEnergyMonomerRange(p,ra,rb,ip);                                // calculate old energy
            if (!parms->nodihpot)
                  for (m=0;m<nmul+3;++m) deltaE -= (((p+ip)->back)+iw+m)->e_dih;               // old dihedral energy
            if (!parms->noangpot)
            {
				deltaE -= (((p+ip)->back)+iw+nmul+1)->e_ang ;                                  // old angular energy
				deltaE -= (((p+ip)->back)+iw+nmul)->e_ang;
			}
		    if (!parms->nohfields)
				deltaE -= GetHFields(p,parms->npol);													  // old H fields energy

            for (m=0;m<nmul;++m)                                                               // apply nmul pivot moves
				if ( (((p+ip)->back)+iw+m)->move == 1 )                                        // PivotForward(k) changes dihedral k+1
				{
					if (parms->randdw==1) dw = parms->dw_lpivot * (0.5 - frand());                      // dihedral to pivot
					else if (parms->randdw==2) dw = parms->dw_lpivot * gauss(0.,1.);

					ok *= PivotForward((p+ip),iw-1+m,dw,nmul-m,parms);                            // moves iw+1+m => overall moves in [iw+1,iw+nmul]
				}
            if (iw+nmul+1<(p+ip)->nback)
                  if ( DAbs(Dist2( (((p+ip)->back)+iw+nmul)->pos, (((p+ip)->back)+iw+nmul+1)->pos ) - (((p+ip)->back)+iw+nmul)->d2_next) > rc2 ) ok=0;           // if covalent bond is broken
            // the stretched bond changes the angles at iw+nmul, iw+nmul+1 and the dihedrals iw+nmul..iw+nmul+2
            if (ok) ok = CheckAngleRange(p,ip,iw+nmul,iw+nmul+1,parms->a_cloose) && CheckLockedDihedralRange(p,ip,iw+nmul,iw+nmul+2,parms->d_cloose);

            if (ok==1)
            {
				ok = AddSidechain(p,ra,rb,ip);

				deltaE += EnergyMonomerRange(p,pot,ra,rb,ip,parms->npol,parms->shell,1,
                                                                                    parms->nosidechains,parms->disentangle,parms->hb);                // new energy of iw with the others
				if (!parms->nodihpot)                                                                                                               // new dihedral energy
					for (m=0;m<nmul+3;++m) deltaE += EnergyDihedrals(p,pot,iw+m,ip,1);
				if (!parms->noangpot)
				{
					deltaE += EnergyAngles(p,pot,iw+nmul+1,ip,1);                                                                     // new angular energy
					deltaE += EnergyAngles(p,pot,iw+nmul,ip,1);
				}
				
				if (!parms->nohfields)
				{
					deltaE += EnergyHFieldsAll(p,pot,parms);
				}
				
            }
      }

      #ifdef DEBUG
      if (parms->debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
      #endif

      if (ok==1) ok = Metropolis(deltaE,t,p->tables);

      if (ok == 0)                  // move rejected
      {
            if (idir == -1) UpdateMonomerRange(oldp,p,ra,rb,ip,parms->shell);            // return to the old position, contacts, etc. (if moved first half)
            else UpdateMonomerRange(oldp,p,ra,rb,ip,parms->shell);                       // ... if moved second half
      }
      else                          // move accepted
      {
            if (idir == -1) UpdateMonomerRange(p,oldp,ra,rb,ip,parms->shell);                        // update oldp (if moved first half)
            else UpdateMonomerRange(p,oldp,ra,rb,ip,parms->shell);                       // ... if moved second half
            p->etot += deltaE;
            parms->acc ++;
      }

      #ifdef DEBUG
      if (parms->debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
      #endif

      parms->mov ++;

      return ok;
}






void SoftExit(int sig)
{
	(void) sig;
	softexit = 1;
}

/*****************************************************************************
 Make a flip of a whole piece of the chain
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveMultipleFlip(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t)
{
	int iw1,iw2,ip,ok=0,i,out,ra,rb;
	double dw,deltaE=0;

	ip = irand(parms->npol);						// which chain to move
	if ((p+ip)->nback < 10) { parms->mov++; return 0; }			// no room for a multiple flip
	iw1 = irand( (p+ip)->nback-6 ) + 3;					// which backbone to move as first
	// the loop below needs some i in [-nmul_mflip,nmul_mflip-1] with |i|>=3 and iw1+i in [3,nback-4]
	if ( (iw1-3 < 3 || parms->nmul_mflip < 3) && (iw1+3 > (p+ip)->nback-4 || parms->nmul_mflip < 4) ) { parms->mov++; return 0; }
	do { i = irand(2*parms->nmul_mflip) - parms->nmul_mflip; } while ( iw1+i > (p+ip)->nback-4 || iw1+i < 3 || (i>-3 && i<3) );
	iw2 = iw1 + i;								// ... which to move as last
	if (iw2<iw1) { i = iw1; iw1=iw2; iw2=i; }				// iw1 should be lower than iw2 for geometry

	#ifdef DEBUG
	 if (debug>2) fprintf(stderr,"step=%d mflip ip=%d iw1=%d iw2=%d\n",istep,ip,iw1,iw2); fflush(stderr);
	#endif

	// old energies
	ra = iw1-1; rb = iw2+1;						// and the sidechains built from [iw1,iw2]
	SideRange(parms,ip,(p+ip)->nback,iw1,iw2,&ra,&rb);
	deltaE = - GetEnergyMonomerRange(p,ra,rb,ip);			// old energy  (iw1-1 & iw2+1 because I move their sidechains)

	if (!parms->noangpot)
	{
		if (iw1>0) deltaE -= (((p+ip)->back)+iw1-1)->e_ang ;		// old angular energy
		if (iw2<(p+ip)->nback-1) deltaE -= (((p+ip)->back)+iw2+1)->e_ang;
	}
	if (!parms->nodihpot)
	{
		if (iw1>0) deltaE -= (((p+ip)->back)+iw1-1)->e_dih;		// old dihedral energy
		deltaE -= (((p+ip)->back)+iw1)->e_dih;
		if (iw2<(p+ip)->nback-1) deltaE -= (((p+ip)->back)+iw2+1)->e_dih;
		if (iw2<(p+ip)->nback-2) deltaE -= (((p+ip)->back)+iw2+2)->e_dih;
	}
	if (!parms->nohfields)
		deltaE -= GetHFields(p,parms->npol);								// old H fields energy

	// make the move in p

	dw = 2. * parms->dw_mflip * (0.5 - frand());				// angle to flip
	ok = FlipFragment((p+ip),iw1,iw2,dw);					// make the flip

	if (ok==0)
	{
		UpdateMonomerRange(oldp,p,ra,rb,ip,parms->shell);		// return to the old position, contacts, etc.
		parms->mov ++;
		return 0;
	}

	// check if constrains are violated
	if (parms->a_cloose>0)
	{
		if (iw1-2>=0)
			if ( DAbs( Angle( (((p+ip)->back)+iw1-2)->pos, (((p+ip)->back)+iw1-1)->pos, (((p+ip)->back)+iw1)->pos, (p+ip)->tables, &out)
			- (((p+ip)->back)+iw1-1)->a_next) > parms->a_cloose  || out==0 ) ok=0;
		if (ok && iw2+2<(p+ip)->nback)
			if ( DAbs( Angle( (((p+ip)->back)+iw2)->pos, (((p+ip)->back)+iw2+1)->pos, (((p+ip)->back)+iw2+2)->pos, (p+ip)->tables, &out)
					- (((p+ip)->back)+iw2+1)->a_next) > parms->a_cloose  || out==0 ) ok=0;
	}
	// flipping [iw1,iw2] changes the dihedrals iw1-1, iw1, iw2+1 and iw2+2
	if (ok) ok = CheckLockedDihedralRange(p,ip,iw1-1,iw1,parms->d_cloose) && CheckLockedDihedralRange(p,ip,iw2+1,iw2+2,parms->d_cloose);

	if (ok) ok = AddSidechain(p,ra,rb,ip);


	// calculate new energy
	if (ok)
	{
		deltaE += EnergyMonomerRange(p,pot,ra,rb,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);		// new energy of iw with the others
		if (!parms->noangpot)
		{
			if (iw1>0) deltaE += EnergyAngles(p,pot,iw1-1,ip,1);
			if (iw2<(p+ip)->nback-1) deltaE += EnergyAngles(p,pot,iw2+1,ip,1);
		}
		if (!parms->nodihpot)
		{
			if (iw1>0) deltaE += EnergyDihedrals(p,pot,iw1-1,ip,1);
			deltaE += EnergyDihedrals(p,pot,iw1,ip,1);
			if (iw2<(p+ip)->nback-1) deltaE += EnergyDihedrals(p,pot,iw2+1,ip,1);
			if (iw2<(p+ip)->nback-2) deltaE += EnergyDihedrals(p,pot,iw2+2,ip,1);
		}
		
		if (!parms->nohfields)
		{
			deltaE += EnergyHFieldsAll(p,pot,parms);
		}

		#ifdef DEBUG
		if (debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
		#endif
	}
	// Metropolis acceptance
	if (ok) ok = Metropolis(deltaE,t,p->tables);

	if (ok == 0) 			// move rejected
	{
		UpdateMonomerRange(oldp,p,ra,rb,ip,parms->shell);							// return to the old position, contacts, etc.
	}
	else					// move accepted
	{
		UpdateMonomerRange(p,oldp,ra,rb,ip,parms->shell);							// update oldp
		p->etot += deltaE;
		parms->acc ++;

		#ifdef DEBUG
			if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
		#endif
	}

	parms->mov ++;

	if (ok == 0) return 0;
	else return 1;
}

/*****************************************************************************
 Moves center of mass of a chain
	 returns 1 if it could make the move, 0 if nothing was done
 *****************************************************************************/
int MoveCoM(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
		unsigned long long istep, int debug, double t)
{
	int ip,ok;
	double dx,dy,dz,deltaE;

	ip = irand(parms->npol);
	dx = parms->dx_com * frand() - parms->dx_com / 2;
	dy = parms->dx_com * frand() - parms->dx_com / 2;
	dz = parms->dx_com * frand() - parms->dx_com / 2;

	#ifdef DEBUG
	 if (debug>2) fprintf(stderr,"step=%d com ip=%d\n",istep,ip); fflush(stderr);
	#endif

	// old energies
	deltaE=-GetEnergyMonomerRange(p,0,(p+ip)->nback-1,ip);						// old energy  (iw1-1 & iw2+1 because I move their sidechains)
	
	if (!parms->nohfields)
	{
		deltaE -= GetHFields(p,parms->npol);										// old H fields energy
	}
	
	
	// move
	DisplaceCoM(p,ip,dx,dy,dz);

	// new energy
	deltaE+=EnergyMonomerRange(p,pot,0,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
	
	// new hfields
	if (!parms->nohfields)
	{
		deltaE += EnergyHFieldsAll(p,pot,parms);
	}
	
#ifdef DEBUG
  	   if (debug>2) fprintf(stderr,"deltaE=%lf\n",deltaE); fflush(stderr);
  	#endif

	// Metropolis acceptance
	ok=Metropolis(deltaE,t,p->tables);

	if (ok == 0) 				// move rejected
	{
		UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,0);							// return to the old position, contacts, etc.
	}
	else					// move accepted
	{
		UpdateMonomerRange(p,oldp,0,(p+ip)->nback-1,ip,0);							// update oldp
	
                if(parms->shell==1 )            // if shells are active, update them 
                {
	//		UpdateShell(p,parms);
          //              CopyShell(p,oldp,parms);
                        parms->ishell=1;
                }

		p->etot+=deltaE;
		parms->acc ++;

		#ifdef DEBUG
			if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
		#endif
	}

	parms->mov++;

	return ok;
}

int MoveRotation(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
                unsigned long long istep, int debug, double t)
{
	int ip,ok,idir;
	double dtheta,deltaE;
	ip=irand(parms->npol);
	
	dtheta= parms->dtheta * (0.5 - frand());


	deltaE = - GetEnergyMonomerRange(p,0,(p+ip)->nback-1,ip);

	if (!parms->nohfields)
		deltaE -= GetHFields(p,parms->npol);								// old H fields energy

	idir=irand(3);
	
	if(idir==0)
		RotationX(p,ip,dtheta);
	else if(idir==1)
		RotationY(p,ip,dtheta);
	else
		RotationZ(p,ip,dtheta);

	
	
	deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
	
	if (!parms->nohfields)
		deltaE += EnergyHFieldsAll(p,pot,parms);

	ok = Metropolis(deltaE,t,p->tables);

        if (ok == 0)                            // move rejected
        {
                UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,0);                                                     // return to the old position, contacts, etc.
        }
        else                                    // move accepted
        {

                UpdateMonomerRange(p,oldp,0,(p+ip)->nback-1,ip,0); 
       		if(parms->shell==1 )            // if shells are active, update them 
        	{
           //     	UpdateShell(p,parms);
         //       	CopyShell(p,oldp,parms);
                        parms->ishell=1;
        	}

		p->etot += deltaE;
        	parms->acc ++;

        	#ifdef DEBUG
	        	if (debug>2) fprintf(stderr,"accept=%d\n",ok); fflush(stderr);
        	#endif
        }

        parms->mov++;
	if (ok == 0) return 0;
        else return 1;

}



/*****************************************************************************
 Returns 1 if the chains of cluster ic have no contact with other chains
 *****************************************************************************/
static int ClusterIsIsolated(struct s_polymer *p, int *members, int n, int *cluster_of, int ic)
{
	int i,ires,k;
	struct s_back *b;

	for (i=0;i<n;++i)
		for (ires=0;ires<(p+members[i])->nback;++ires)
		{
			b = ((p+members[i])->back)+ires;
			for (k=0;k<b->ncontacts;++k)
				if (cluster_of[(b->contacts_p)[k]] != ic) return 0;
		}
	return 1;
}

/*****************************************************************************
 Periodic box: shift the chains of a cluster by box vectors (which does not
 change any energy) so that each one is at the image nearest to the chain
 through which it is reached. Returns 0 if the cluster wraps around the box
 (e.g. a fibril through the periodic boundary), which then cannot be rotated.
 *****************************************************************************/
static int UnwrapCluster(struct s_polymer *p, int *members, int n, int *cluster_of, int ic, double L)
{
	int head,i,k,a,c2,placed[n],idx;
	double dx,dy,dz;
	struct s_back *b,*b2;

	for (i=0;i<n;++i) placed[i] = (i==0);
	for (head=0;head<n;++head)
	{
		a = members[head];
		if (!placed[head]) continue;
		for (i=0;i<(p+a)->nback;++i)
		{
			b = ((p+a)->back)+i;
			for (k=0;k<b->ncontacts;++k)
			{
				c2 = (b->contacts_p)[k];
				if (cluster_of[c2] != ic || c2 == a) continue;
				for (idx=0;members[idx]!=c2;++idx) ;
				b2 = ((p+c2)->back)+(b->contacts)[k];
				dx = b2->pos.x - b->pos.x;
				dy = b2->pos.y - b->pos.y;
				dz = b2->pos.z - b->pos.z;
				if (!placed[idx])
				{
					DisplaceCoM(p,c2,-L*round(dx/L),-L*round(dy/L),-L*round(dz/L));
					placed[idx] = 1;
				}
				else if (fabs(dx) > 0.5*L || fabs(dy) > 0.5*L || fabs(dz) > 0.5*L) return 0;
			}
		}
	}
	for (i=0;i<n;++i) if (!placed[i]) return 0;				// the visit of FindClusters reaches all of them
	return 1;
}

/*****************************************************************************
 Rotate the chains members[0..n-1] by dtheta (radians) around the axis idir
 (0=x, 1=y, 2=z) through the centre of mass of their backbone atoms
 *****************************************************************************/
static void RotateChains(struct s_polymer *p, int *members, int n, int idir, double dtheta)
{
	int k,i,j,nat=0;
	double c[3]={0,0,0},cs=cos(dtheta),sn=sin(dtheta),u,v;
	double *x,*y;
	struct s_back *b;
	struct vector *r;

	for (k=0;k<n;++k)
		for (i=0;i<(p+members[k])->nback;++i)
		{
			b = ((p+members[k])->back)+i;
			c[0] += b->pos.x; c[1] += b->pos.y; c[2] += b->pos.z;
			nat++;
		}
	for (i=0;i<3;++i) c[i] /= nat;

	for (k=0;k<n;++k)
		for (i=0;i<(p+members[k])->nback;++i)
			for (j=-1;j<(((p+members[k])->back)+i)->nside;++j)	// j=-1 is the backbone atom
			{
				r = (j<0) ? &((((p+members[k])->back)+i)->pos) : &((((((p+members[k])->back)+i)->side)+j)->pos);
				if (idir==0)      { x = &(r->y); y = &(r->z); u = r->y - c[1]; v = r->z - c[2]; *x = c[1] + cs*u - sn*v; *y = c[2] + sn*u + cs*v; }
				else if (idir==1) { x = &(r->z); y = &(r->x); u = r->z - c[2]; v = r->x - c[0]; *x = c[2] + cs*u - sn*v; *y = c[0] + sn*u + cs*v; }
				else              { x = &(r->x); y = &(r->y); u = r->x - c[0]; v = r->y - c[1]; *x = c[0] + cs*u - sn*v; *y = c[1] + sn*u + cs*v; }
			}
}

/*****************************************************************************
 Rigid translation (rot=0) or rotation (rot=1) of a cluster of chains.
 A move after which the cluster touches other chains is rejected: it would
 change the decomposition into clusters, making the reverse move impossible.
 *****************************************************************************/
static int MoveCluster(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot, double t, int rot)
{
	int npol=parms->npol,member[npol],start[npol+1],cluster_of[npol],wfrom[npol],wto[npol];
	int ncluster,ic,n,i,ok,idir=0,*m;
	double dx=0,dy=0,dz=0,dtheta=0,deltaE=0.;

	if (rot) dtheta = parms->dtheta * (0.5 - frand());

	ncluster = FindClusters(p,npol,member,start,cluster_of);
	ic = irand(ncluster);
	m = member+start[ic];
	n = start[ic+1]-start[ic];
	for(i=0;i<n;++i)
	{
		wfrom[i] = 0;
		wto[i] = (p+m[i])->nback-1;
	}

	if (rot)
	{
		idir = irand(3);
		if (pot->boxtype=='p' && !UnwrapCluster(p,m,n,cluster_of,ic,pot->lbox))
		{
			UpdateChainRanges(oldp,p,n,m,wfrom,wto);			// undo the shifts
			parms->mov ++;
			return 0;
		}
		for(i=0;i<n;++i) deltaE -= GetEnergyMonomerRange(p,0,(p+m[i])->nback-1,m[i]);
		RotateChains(p,m,n,idir,dtheta);
		if (p->cells)							// all chains moved before their energies are computed
			for(i=0;i<n;++i) CellsUpdateRange(p,m[i],0,(p+m[i])->nback-1);
		for(i=0;i<n;++i)
			deltaE += EnergyMonomerRange(p,pot,0,(p+m[i])->nback-1,m[i],npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
	}
	else
	{
		// all chains are displaced before any energy is computed: moving them one at a
		// time would make them overlap with their neighbours not moved yet
		dx = parms->dx_clm * frand() - parms->dx_clm / 2;
		dy = parms->dx_clm * frand() - parms->dx_clm / 2;
		dz = parms->dx_clm * frand() - parms->dx_clm / 2;
		for(i=0;i<n;++i) deltaE -= GetEnergyMonomerRange(p,0,(p+m[i])->nback-1,m[i]);
		for(i=0;i<n;++i) DisplaceCoM(p,m[i],dx,dy,dz);
		if (p->cells)
			for(i=0;i<n;++i) CellsUpdateRange(p,m[i],0,(p+m[i])->nback-1);
		for(i=0;i<n;++i)
			deltaE += EnergyMonomerRange(p,pot,0,(p+m[i])->nback-1,m[i],npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
	}

	if (ClusterIsIsolated(p,m,n,cluster_of,ic))
		ok = Metropolis(deltaE,t,p->tables);
	else
		ok = 0;

	if (ok == 0)                        // move rejected: return to the old position, contacts, etc.
		UpdateChainRanges(oldp,p,n,m,wfrom,wto);
	else                                // move accepted
	{
		UpdateChainRanges(p,oldp,n,m,wfrom,wto);
		if(parms->shell==1) parms->ishell=1;
		p->etot += deltaE;
		parms->acc ++;
	}
	parms->mov ++;

	return ok;
}

int MoveClusterCoM(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
                  unsigned long long istep, int debug, double t, FILE *fproc, int iproc)
{
	return MoveCluster(p,oldp,parms,pot,t,0);
}

int MoveClusterRot(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot,
                  unsigned long long istep, int npol,int debug, double t, FILE *fproc, int iproc)
{
	return MoveCluster(p,oldp,parms,pot,t,1);
}


/*****************************************************************************
 Aggregation-volume-bias move (AVBMC, Chen & Siepmann, J. Phys. Chem. B 105,
 11275 (2001)): a target chain j is chosen; with probability avb_pin a chain
 out of the bonded region of j (centre of mass farther than avb_rin) is put
 at a random point of that region, otherwise a chain in the region is put at
 a random point of the rest of the box. The moved chain keeps its internal
 conformation and takes a random orientation. Periodic box only, so that the
 volumes of the two regions are exact (V_in + V_out = L^3).
 *****************************************************************************/
static struct vector ChainCentre(struct s_polymer *p, int ic)
{
	struct vector c={0,0,0};
	int i;

	for (i=0;i<(p+ic)->nback;++i)
	{
		c.x += (((p+ic)->back)+i)->pos.x;
		c.y += (((p+ic)->back)+i)->pos.y;
		c.z += (((p+ic)->back)+i)->pos.z;
	}
	c.x /= (p+ic)->nback; c.y /= (p+ic)->nback; c.z /= (p+ic)->nback;
	return c;
}

// rotate chain ic by the matrix R around its centre c, then put the centre at c2
static void PlaceChain(struct s_polymer *p, int ic, double R[3][3], struct vector c, struct vector c2)
{
	int i,j;
	double u,v,w;
	struct vector *r;

	for (i=0;i<(p+ic)->nback;++i)
		for (j=-1;j<(((p+ic)->back)+i)->nside;++j)			// j=-1 is the backbone atom
		{
			r = (j<0) ? &((((p+ic)->back)+i)->pos) : &((((((p+ic)->back)+i)->side)+j)->pos);
			u = r->x - c.x; v = r->y - c.y; w = r->z - c.z;
			r->x = c2.x + R[0][0]*u + R[0][1]*v + R[0][2]*w;
			r->y = c2.y + R[1][0]*u + R[1][1]*v + R[1][2]*w;
			r->z = c2.z + R[2][0]*u + R[2][1]*v + R[2][2]*w;
		}
}

// uniformly distributed rotation matrix (random unit quaternion, Shoemake)
static void RandomRotation(double R[3][3])
{
	double u1=frand(),u2=2*PI*frand(),u3=2*PI*frand();
	double a=sqrt(1-u1)*sin(u2),b=sqrt(1-u1)*cos(u2),c=sqrt(u1)*sin(u3),d=sqrt(u1)*cos(u3);	// q = d + a i + b j + c k

	R[0][0] = 1-2*(b*b+c*c); R[0][1] = 2*(a*b-c*d);   R[0][2] = 2*(a*c+b*d);
	R[1][0] = 2*(a*b+c*d);   R[1][1] = 1-2*(a*a+c*c); R[1][2] = 2*(b*c-a*d);
	R[2][0] = 2*(a*c-b*d);   R[2][1] = 2*(b*c+a*d);   R[2][2] = 1-2*(a*a+b*b);
}

int MoveAVB(struct s_polymer *p,struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot, double t)
{
	int npol=parms->npol,j,k,ip,in_move,nin=0,nout=0,cand[npol],ncand=0,ok;
	double L=pot->lbox,rin=parms->avb_rin,r2,vin,vout,factor,deltaE,R[3][3];
	struct vector cj,ci,c2,d;

	j = irand(npol);
	in_move = (frand() < parms->avb_pin);
	cj = ChainCentre(p,j);

	// chains in and out of the bonded region of j; the candidates to move
	for (k=0;k<npol;++k)
		if (k!=j)
		{
			ci = ChainCentre(p,k);
			d.x = ci.x-cj.x; d.y = ci.y-cj.y; d.z = ci.z-cj.z;
			d.x -= L*round(d.x/L); d.y -= L*round(d.y/L); d.z -= L*round(d.z/L);
			if (d.x*d.x+d.y*d.y+d.z*d.z < rin*rin) { nin++; if (in_move==0) cand[ncand++] = k; }
			else { nout++; if (in_move==1) cand[ncand++] = k; }
		}

	// the random numbers of the proposal are always drawn, also when nothing can be moved
	ip = ncand ? cand[irand(ncand)] : -1;
	RandomRotation(R);
	do
	{
		d.x = (frand()-0.5)*L; d.y = (frand()-0.5)*L; d.z = (frand()-0.5)*L;
		r2 = d.x*d.x+d.y*d.y+d.z*d.z;
	}
	while ( in_move ? (r2 >= rin*rin) : (r2 < rin*rin) );			// uniform in the sphere, or in the box out of it

	if (ip<0) { parms->mov++; return 0; }

	vin = 4./3.*PI*rin*rin*rin;
	vout = L*L*L - vin;
	if (in_move) factor = (1-parms->avb_pin)*vin*nout / (parms->avb_pin*vout*(nin+1));
	else factor = parms->avb_pin*vout*nin / ((1-parms->avb_pin)*vin*(nout+1));

	deltaE = - GetEnergyMonomerRange(p,0,(p+ip)->nback-1,ip);
	if (!parms->nohfields) deltaE -= GetHFields(p,npol);

	ci = ChainCentre(p,ip);
	c2.x = cj.x+d.x; c2.y = cj.y+d.y; c2.z = cj.z+d.z;
	PlaceChain(p,ip,R,ci,c2);

	deltaE += EnergyMonomerRange(p,pot,0,(p+ip)->nback-1,ip,npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
	if (!parms->nohfields) deltaE += EnergyHFieldsAll(p,pot,parms);

	ok = Metropolis(deltaE - t*log(factor),t,p->tables);		// acceptance min[1, factor exp(-deltaE/T)]

	if (ok == 0)
		UpdateMonomerRange(oldp,p,0,(p+ip)->nback-1,ip,0);
	else
	{
		UpdateMonomerRange(p,oldp,0,(p+ip)->nback-1,ip,0);
		if(parms->shell==1) parms->ishell=1;
		p->etot += deltaE;
		parms->acc ++;
	}
	parms->mov ++;

	return ok;
}


/*****************************************************************************
 Fast energy of backbone atom j of chain ip with the backbone atoms already in
 place during a regrowth (dir>0: those before j in its chain, dir<0: those
 after it; all the other chains) and of the dihedral it completes
 *****************************************************************************/
static double CBMCEnergy(struct s_polymer *p, struct s_potential *pot, struct s_mc_parms *parms, int ip, int j, int dir)
{
	int cq,q,q0,q1,tooclose;
	double e,etot=0;

	for (cq=0;cq<parms->npol;++cq)
	{
		q0 = 0; q1 = (p+cq)->nback-1;
		if (cq==ip) { if (dir>0) q1 = j-1; else q0 = j+1; }
		for (q=q0;q<=q1;++q)
		{
			tooclose = (cq==ip && Abs(j-q) <= pot->g_imin);
			e = EnergyPair(p,pot,j,q,ip,cq,0,1,parms->disentangle,tooclose,0);	// backbone only, contacts untouched
			if (e>=LARGE && !parms->disentangle) return LARGE;
			etot += e;
		}
	}
	if (!parms->nodihpot) etot += EnergyDihedrals(p,pot,dir>0 ? j-1 : j+2,ip,0);
	return etot;
}

// log of the sum of exp(-u_k/T) over the trials with u_k < LARGE, summed relative to
// the lowest of them (umin, LARGE if none) so that it neither overflows nor underflows
static double CBMCLogWeight(double *u, int nt, double t, double *umin)
{
	double w=0;
	int k;

	*umin = LARGE;
	for (k=0;k<nt;++k) if (u[k] < *umin) *umin = u[k];
	if (*umin >= LARGE) return 0;
	for (k=0;k<nt;++k) if (u[k] < LARGE) w += exp(-(u[k]-*umin)/t);
	return log(w) - *umin/t;
}

/*****************************************************************************
 Configurational-bias regrowth of an end of a chain (Siepmann and Frenkel,
 Mol. Phys. 75, 59 (1992)): the last n backbone atoms of one end (n uniform
 in 1..nmax_cbmc) are grown again one at a time, keeping their bond lengths and
 angles. The dihedral of each one (if its bond can rotate) is chosen among
 ntrial_cbmc random trials with the Boltzmann weight of a fast energy (CBMCEnergy);
 the old conformation is weighted in the same way, with its own dihedral as
 one of the trials. The move is accepted with
	min[1, W_new/W_old exp(-(dE - (u_new - u_old))/T)]
 where W are the Rosenbluth weights (products over the atoms of the sums of
 the weights of the trials), dE the true change of energy and u the sums of
 the fast energies of the chosen trials.
 *****************************************************************************/
int MoveCBMC(struct s_polymer *p, struct s_polymer *oldp, struct s_mc_parms *parms, struct s_potential *pot, double t)
{
	int ip,nb,n,dir,m0,m1,ra,rb,kd0,kd1,k,i,j,nt,ok=1,out,K=parms->cbmc_k;
	double deltaE,usel_old=0,usel_new=0,lnw_old=0,lnw_new=0,lnw,umin,r,u[K];
	struct vector x[K],xold;
	struct s_back *bk;

	ip = irand(parms->npol);
	nb = (p+ip)->nback;
	bk = (p+ip)->back;
	dir = irand(2) ? 1 : -1;						// 1: regrow the end of the chain, -1: its beginning
	n = 1 + irand(parms->cbmc_n);
	if (n > nb-3) { parms->mov++; return 0; }				// three atoms must stay to place the first one
	if (dir>0) { m0 = nb-n; m1 = nb-1; kd0 = m0-1; kd1 = nb-2; }
	else { m0 = 0; m1 = n-1; kd0 = 2; kd1 = m1+2; }			// the dihedrals k (atoms k-2..k+1) that change
	if (kd0 < 2) kd0 = 2;
	if (kd1 > nb-2) kd1 = nb-2;

	struct angles g[n];						// bond length and angle of the regrown atoms (kept)

	ra = m0; rb = m1;						// and the sidechains built from them
	SideRange(parms,ip,nb,m0,m1,&ra,&rb);

	// old energy
	deltaE = -GetEnergyMonomerRange(p,ra,rb,ip);
	if (!parms->nodihpot) for (k=kd0;k<=kd1;++k) deltaE -= (bk+k)->e_dih;
	if (!parms->nohfields) deltaE -= GetHFields(p,parms->npol);

	// Rosenbluth weight of the old conformation, retracing its growth
	for (i=0;i<n;++i)
	{
		j = (dir>0) ? m0+i : m1-i;
		if (dir>0) g[i] = Cartesian2Spherical((bk+j-3)->pos,(bk+j-2)->pos,(bk+j-1)->pos,(bk+j)->pos,p->tables,&out);
		else g[i] = Cartesian2Spherical((bk+j+3)->pos,(bk+j+2)->pos,(bk+j+1)->pos,(bk+j)->pos,p->tables,&out);
		if (!out) ok = 0;
		nt = ((bk+(dir>0 ? j-1 : j+2))->move == 1) ? K : 1;		// a dihedral which cannot rotate is kept
		xold = (bk+j)->pos;
		u[0] = CBMCEnergy(p,pot,parms,ip,j,dir);
		for (k=1;k<nt;++k)
		{
			struct angles gt = g[i];
			gt.dih = 360.*frand() - 180.;
			if (dir>0) (bk+j)->pos = Spherical2Cartesian((bk+j-3)->pos,(bk+j-2)->pos,(bk+j-1)->pos,gt,p->tables,&out);
			else (bk+j)->pos = Spherical2Cartesian((bk+j+3)->pos,(bk+j+2)->pos,(bk+j+1)->pos,gt,p->tables,&out);
			u[k] = out ? CBMCEnergy(p,pot,parms,ip,j,dir) : LARGE;
		}
		(bk+j)->pos = xold;
		if (u[0] >= LARGE) ok = 0;
		else { lnw_old += CBMCLogWeight(u,nt,t,&umin); usel_old += u[0]; }
	}

	// grow the new conformation
	for (i=0;i<n && ok;++i)
	{
		j = (dir>0) ? m0+i : m1-i;
		nt = ((bk+(dir>0 ? j-1 : j+2))->move == 1) ? K : 1;
		for (k=0;k<nt;++k)
		{
			struct angles gt = g[i];
			if (nt > 1) gt.dih = 360.*frand() - 180.;
			if (dir>0) x[k] = Spherical2Cartesian((bk+j-3)->pos,(bk+j-2)->pos,(bk+j-1)->pos,gt,p->tables,&out);
			else x[k] = Spherical2Cartesian((bk+j+3)->pos,(bk+j+2)->pos,(bk+j+1)->pos,gt,p->tables,&out);
			(bk+j)->pos = x[k];
			u[k] = out ? CBMCEnergy(p,pot,parms,ip,j,dir) : LARGE;
		}
		lnw = CBMCLogWeight(u,nt,t,&umin);
		if (umin >= LARGE) { ok = 0; break; }
		r = frand() * exp(lnw + umin/t);				// trial k with probability exp(-u_k/T)/w (weights relative to the lowest u)
		for (k=0;k<nt-1;++k)
		{
			if (u[k] < LARGE) r -= exp(-(u[k]-umin)/t);
			if (r < 0) break;
		}
		while (u[k] >= LARGE) k--;					// (only if rounding left r>=0 past the last allowed trial)
		(bk+j)->pos = x[k];
		lnw_new += lnw;
		usel_new += u[k];
	}

	// true energy of the new conformation
	if (ok) ok = AddSidechain(p,ra,rb,ip);
	if (ok)
	{
		deltaE += EnergyMonomerRange(p,pot,ra,rb,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
		if (!parms->nodihpot) for (k=kd0;k<=kd1;++k) deltaE += EnergyDihedrals(p,pot,k,ip,1);
		if (!parms->nohfields) deltaE += EnergyHFieldsAll(p,pot,parms);
		ok = Metropolis(deltaE - (usel_new-usel_old) - t*(lnw_new-lnw_old),t,p->tables);
	}

	if (ok == 0)
		UpdateMonomerRange(oldp,p,ra,rb,ip,0);
	else
	{
		UpdateMonomerRange(p,oldp,ra,rb,ip,0);
		if(parms->shell==1) parms->ishell=1;
		p->etot += deltaE;
		parms->acc ++;
	}
	parms->mov ++;

	return ok;
}
