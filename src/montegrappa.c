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
 * montegrappa.c
 *
 *  Created on: Sep 15, 2010
 *      Author: guido
 */

#include "montegrappa.h"
#include <unistd.h>

int main(int argc, char *argv[])
{
	struct s_polymer *polymer, *startp, *oldp,*fragment=NULL,*replica;
	struct s_mc_parms *parms=0;
	struct s_potential *u;
	char polname[FNAMEMAX],potname[FNAMEMAX],fname[2*FNAMEMAX];
	int i,nat,ntypes,irun,ci,out,nrestr;
	struct s_checkpoint ckp;
	char ckpname[FNAMEMAX];
	int restart=0;
	int nbackmax,nsidemax,nrotmax,natommax;

	FILE *ftrj=NULL,*fe=NULL,*fproc=NULL,*fp=NULL;


	double **ematold;

	if (argc<4)
	{
		Welcome(stderr);
		fprintf(stderr,"\nusage:\n$ montegrappa polymer.pol potential.pot parameters.par\n\n");
		exit(1);
	}


	int my_rank=0;

	#ifdef ACTIVE_MPI

	int nprocs,j;
	time_t rawtime;
	MPI_Init(&argc,&argv);
	MPI_Status astatus = {0};			// passed by value to the functions that send the structures
	struct s_mpi_parms *mpiparms=malloc(sizeof(struct s_mpi_parms));

	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);	

	Create_back_datatype(&(mpiparms->Back_mpi));
	Create_side_datatype(&(mpiparms->Side_mpi));
	Create_rot_datatype(&(mpiparms->Rot_mpi));
	Create_pot_datatype(&(mpiparms->Pot_mpi));


	mpiparms->my_rank=my_rank;
	mpiparms->nprocs=nprocs;
	#endif

	if (argc<4) Error("Usage: montegrappa <polymer file | checkpoint> <potential file> <parameter file>");
	if (strlen(argv[1])>=FNAMEMAX || strlen(argv[2])>=FNAMEMAX) Error("File name too long (FNAMEMAX)");
	strcpy(polname,argv[1]);
	
	memset(&ckp,0,sizeof(ckp));
	CheckpointName(ckpname,sizeof(ckpname),my_rank);
	if(!strcmp(polname,"checkpoint"))
	{
		// restart: the checkpoint of each rank contains the name of the original polymer file
		restart = 1;
		ReadCheckpointHeader(ckpname,&ckp);
		strcpy(polname,ckp.polname);
		ckp.restart = 1;
		if (my_rank==0) fprintf(stderr,"Restart from %s at run %d, step %llu (polymer file %s)\n",ckpname,ckp.irun,ckp.istep,polname);
	}
	else
		strcpy(ckp.polname,polname);
	
	
	// if kill -SIGUSR1 given, soft exit (on every rank: they must stop together)
	signal(SIGUSR1,SoftExit);

	#ifdef ACTIVE_MPI
	if(my_rank==0)
	{
	#endif
		Welcome(stderr);
		strcpy(potname,argv[2]);


		fprintf(stderr,"Simulation parameters:\n");
		parms = ReadMcParms(argv[3]);
		if (parms->seed==-1) parms->seed = TimeSeed();			// decided once, before it is sent to the other ranks
		//speedup calculation
		parms->r2shell=(parms->r2shell)*(parms->r2shell);
	
		// square the r_contact (hfields)
		parms->r_contact=(parms->r_contact)*(parms->r_contact);
		
	

		#ifndef DEBUG
		if (parms->debug>0) Error("To debug, the program must be compiled with #define DEBUG");
		#else
		fprintf(stderr,"WARNING: executable compiled with DEBUG option\n");
		#endif

		#ifdef ACTIVE_MPI
       	//printf(stderr,"~MPI: %d Procs\n",nprocs);
         	if(nprocs!= parms->ntemp)
           	{
               	Error("WARNING! ntemp must be equal to nprocs");
                 	return -1;
     		}
        	if(parms->nreplicas!=nprocs && parms->nreplicas!=1)
		{
		Error("WARNING! nreplicas must be equal to nprocs");
			return -1;
		}   	
		if(parms->nreplicas!=parms->ntemp&& parms->nreplicas!=1)
		{
		Error("WARNING! nreplicas must be equal to ntemp");
			return -1;
		}
	#endif

	#ifdef ACTIVE_MPI
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	#endif

	#ifdef ACTIVE_MPI
	parms=send_parms(mpiparms->my_rank,nprocs,parms);
	if (strcmp(parms->fnproc,""))
	{
            snprintf(fname,sizeof(fname),"%s_%d.log",parms->fnproc,my_rank);
            fproc = fopen(fname,"w");
            if (!fproc)
		      Error("Cannot open log-process file for writing");
	}

	
	if(parms->debug>1) 
		fprintf(fproc, "iproc=%d\tnrun=%d\tnstep=%llu\tnpol=%d\ttemp=%f\tnstep_exchange=%d\n", my_rank, parms->nrun, parms->nstep, parms->npol, parms->T[mpiparms->my_rank],parms->nstep_exchange);
	#else
	fproc = stderr;
	#endif
	//fprintf(stderr,"proc %d , seed = %d\n",my_rank,parms->seed);
	parms->seed = Randomize(parms->seed,my_rank);		// one independent stream per replica
	
	// allocate only what the polymer file needs (all ranks read the same file)
	ScanPolymerSizes(polname,parms->npol,&nbackmax,&nsidemax,&nrotmax,&natommax);
	startp = AlloPolymer(parms->npol,nbackmax,nsidemax,nrotmax,natommax,parms->shell,parms->nosidechains,parms->flog);
	polymer = AlloPolymer(parms->npol,nbackmax,nsidemax,nrotmax,natommax,parms->shell,parms->nosidechains,parms->flog);
	oldp = AlloPolymer(parms->npol,nbackmax,nsidemax,nrotmax,natommax,parms->shell,parms->nosidechains,parms->flog);
	replica = AlloPolymer(parms->npol,nbackmax,nsidemax,nrotmax,natommax,parms->shell,parms->nosidechains,fproc);

	if( parms->nmul_local !=0 && parms->movetype[7]!=-1)
	{
//      fprintf(stderr,"Fragment length: %d(+1)\n",parms->nmul_local);
//		fprintf(stderr,"npol %d \n",parms->npol);	      
		fragment=Allo_Fragment(parms->npol,parms->nmul_local,parms->nmul_local,parms->flog);
	}
	
	//fprintf(stderr,"\n");
	startp->tables = InitTables(fproc);                         // keep tables only in (polymer+0)

	for (i=0;i<parms->npol;++i) (polymer+i)->tables = startp->tables; // tables of all polymers point to same address
	for (i=0;i<parms->npol;++i) (replica+i)->tables = startp->tables; // tables of all polymers point to same address
	for (i=1;i<parms->npol;++i) (startp+i)->tables = startp->tables;  // tables of all polymers point to same address

	// Read polymer stuff from file
	//if(parms->nreplicas==1)
	//{
	{
		#ifdef ACTIVE_MPI
		if(my_rank==0)
		{
			ReadPolymer(polname,startp,parms->flog,parms->npol, parms->debug, &nat, &ntypes);
			ComputeIAPDB(startp,parms);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		for (i=0;i<parms->npol;++i)
		{
			send_struct(&((startp+i)->nback),mpiparms->my_rank, nprocs, &nat, &ntypes, astatus);
			startp = send_pol(mpiparms->my_rank, nprocs, (startp+i)->nback, mpiparms->Back_mpi, mpiparms->Side_mpi, mpiparms->Rot_mpi, startp, astatus, i, parms->nshell, parms->nosidechains);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(parms->debug>1)
			fprintf(fproc,"iproc=%d\tnat=%d\tntypes=%d\tnback=%d\n",mpiparms->my_rank,nat,ntypes, startp->nback);
		#else
		fprintf(stderr,"POLNAME: %s\n",polname);
		ReadPolymer(polname,startp,parms->flog,parms->npol, parms->debug, &nat, &ntypes);
		ComputeIAPDB(startp,parms);
		#endif
	
	}
	//}//end of if nreplicas==1
        

/*
	else if(parms->nreplicas==parms->ntemp)
	{
		fprintf(stderr,"BBBBBB");
		#ifdef ACTIVE_MPI

			int poly;
			for(poly=parms->ntemp-1;poly>0;poly--)
			{
                		if(my_rank==0)
				{
		        		ReadPolymer(parms->input_names[poly],startp,parms->flog,parms->npol, parms->debug, &nat, &ntypes);
                        		ComputeIAPDB(startp,parms);
		                }
			

				for(i=0;i<parms->npol;++i)
				{
					send_struct_to_proc(&((startp+i)->nback),mpiparms-> my_rank,poly, nprocs, &nat, &ntypes, astatus);
        			        startp = send_pol_to_proc(mpiparms->my_rank,poly,nprocs, (startp+i)->nback, mpiparms->Back_mpi, mpiparms->Side_mpi, mpiparms->Rot_mpi, startp, astatus, i, parms->nshell, parms->nosidechains);


				}	
		                MPI_Barrier(MPI_COMM_WORLD);

			}

		ReadPolymer(parms->input_names[0],startp,parms->flog,parms->npol, parms->debug, &nat, &ntypes);
                ComputeIAPDB(startp,parms);
	
		

		#endif
	}
	MPI_Barrier(MPI_COMM_WORLD);	

*/

	for (i=0;i<parms->npol;++i)
		if(!AddSidechain(startp,0,(startp+i)->nback-1,i))
			Error("Cannot add sidechain in starting chain");

	// Read potential (needs the number of atoms from ReadPolymer)
	u = AlloPotential(nat,ntypes,parms->noangpot,parms->nodihpot,parms->hb,parms->nohfields);

	if(my_rank==0)
	{
		ReadPotential(potname,u,parms,nat,ntypes);
	}

	#ifdef ACTIVE_MPI
	MPI_Bcast(&(parms->nohfields),1,MPI_INT,0,MPI_COMM_WORLD);	// ReadPotential may switch H fields off
	u = send_pot(nat,ntypes,parms->noangpot,parms->nodihpot,parms->nohfields,parms->hb, my_rank, nprocs,mpiparms-> Pot_mpi, u, astatus);
	//POTENZIALE SUI DIEDRI 
	if(u->dih_ram)
	{
        	for(i=0;i<2;++i)
                  for(j=0;j<2;++j)
                        if(u->sigma[i][j]<=0)
                              if(my_rank==0)
                                    Error("Sigma <=0 !!");

        	if(u->dih0[0][0]!= -57)
                  if(my_rank==0)
                        fprintf(stderr,"WARNING! phi0_a is usually equal to -57!!!!!!!!\n");

        	if(u->dih0[1][0]!= -129)
                  if(my_rank==0)
                        fprintf(stderr,"WARNING! phi0_b is usually equal to -129!!!!!!!!\n");
        	if(u->dih0[0][1]!= -47)
                  if(my_rank==0)
                        fprintf(stderr,"WARNING! psi0_a is usually equal to -47!!!!!!!!\n");

        	if(u->dih0[1][1]!= 124)
                  if(my_rank==0)
                        fprintf(stderr,"WARNING! psi0_b is usually equal to 124!!!!!!!!\n");
     	}
	#endif
	   
	#ifdef OPTIMIZEPOT
	char fname_op[50];
	#ifdef ACTIVE_MPI
	snprintf(fname_op,sizeof(fname_op),"op_log_proc%d.pdb",my_rank);
	#else
	snprintf(fname_op,sizeof(fname_op),"op_log");
	#endif
	fp = fopen(fname_op,"w");
	startp->op = InitializeOptimizePot(parms,ntypes,u,fproc,my_rank,fp,&nrestr);	

	ematold = AlloDoubleMatrix(ntypes,ntypes);
	#endif
	// the shell list is rebuilt only every nshell steps: its radius must exceed the
	// interaction range by a skin, otherwise pairs entering the range are missed
	if (parms->shell)
	{
		double r2max=0., rshell;
		int it,jt;
		for (it=0;it<ntypes;++it)
			for (jt=0;jt<ntypes;++jt)
				if ((u->r_2)[it][jt] > r2max) r2max = (u->r_2)[it][jt];
		rshell = sqrt(r2max) + SHELL_SKIN;
		if (parms->r2shell < rshell*rshell)
		{
			if (my_rank==0) fprintf(stderr,"WARNING: r2shell=%lf is too small for interactions up to %lf A: using %lf\n",
						sqrt(parms->r2shell),sqrt(r2max),rshell);
			parms->r2shell = rshell*rshell;
		}
	}

	SetPairCutoff(u,ntypes);
	if (u->boxtype=='p' && parms->shell) Error("shells cannot be used with a periodic box (the cell list does the same job)");
	if (parms->movetype[13]!=-1 && (u->boxtype!='p' || 2*parms->avb_rin >= u->lbox)) Error("avbmc needs a periodic box of side larger than 2*avb_rin");
	if (parms->hb && my_rank==0 && CheckHBReferences(startp,u,parms->npol,nat))
		fprintf(stderr,"WARNING: the energy of these hydrogen bonds is not updated when only the reference atom moves\n");
	SetSidechainExtents(startp,parms->npol,parms->nosidechains);
	// pivot moves skip the pairs within the part they rotate rigidly (not with disentangle, where overlaps count)
	parms->rigid = parms->disentangle ? -1 : RigidMargin(startp,u,parms->npol,nat,parms->hb,parms->nosidechains);
	fprintf(fproc,"Rigid margin of pivot moves: %d\n",parms->rigid);

	// cell list for the energy of moved residues (not with shells, which play the same role)
	{
		int nres=0;
		for (i=0;i<parms->npol;++i) nres += (startp+i)->nback;
		if (!parms->shell)
			CellsAlloc(polymer,startp,parms->npol,parms->cells==1,parms->cells==-1 && nres>=CELLS_AUTO_MIN);
	}
	CheckOverlaps(startp,u,parms,parms->npol,parms->nosidechains,1,fproc);

	// calculate initial distances, angles and dihedrals to constrain loose moves
	if(parms->debug>1)
            fprintf(stderr,"Initial distances, angles and dihedrals: ");
	
	for(ci=0;ci<parms->npol;++ci)
		for(i=0;i<(startp+ci)->nback-1;++i)
		{
		      (((startp+ci)->back)+i)->d2_next = Dist2( (((startp+ci)->back)+i)->pos, (((startp+ci)->back)+i+1)->pos );

			if (i>0) 
			{
			      (((startp+ci)->back)+i)->a_next = Angle( (((startp+ci)->back)+i-1)->pos , (((startp+ci)->back)+i)->pos,
					(((startp+ci)->back)+i+1)->pos, startp->tables, &out);
				if (!out) Error("Zero angle in initial structure");
			}

			if (i>1)
			{ 
				(((startp+ci)->back)+i)->d_next = Dihedral( (((startp+ci)->back)+i-2)->pos, (((startp+ci)->back)+i-1)->pos ,
					(((startp+ci)->back)+i)->pos, (((startp+ci)->back)+i+1)->pos, startp->tables, &out);
				if (!out) Error("Flat angle in initial structure");
			}

			if(parms->debug>1)
                        fprintf(stderr,"d2(%d %d)=%lf ",ci,i,(((startp+ci)->back)+i)->d2_next);

			if(parms->debug>1 && i>0)
                        fprintf(stderr,"ang(%d %d)=%lf ",ci,i,(((startp+ci)->back)+i)->a_next);

			if (parms->debug>1 && i>1)
                        fprintf(stderr,"dih(%d %d)=%lf ",ci,i,(((startp+ci)->back)+i)->d_next);
		}
	
        fprintf(stderr,"\n");
	#ifdef ACTIVE_MPI
	fprintf(stderr,"> rank %d Starting simulation...\n",my_rank);
	#else
	fprintf(stderr,"Starting simulation...\n");
	#endif
	// do Monte Carlo
	#ifdef OPTIMIZEPOT
	if(parms->chi2start==0 && !restart)
		remove("chi2.dat");
	#endif
	for (irun=(restart ? ckp.irun : parms->chi2start);irun<parms->nrun;++irun)
	{
		if (ckp.restart)						// continue the interrupted run from the checkpoint
		{
			CopyAllPolymers(startp,polymer,parms->npol,parms->nosidechains,parms->nosidechains);
			ReadCheckpointCoordinates(ckpname,polymer,parms->npol);
			#ifdef ACTIVE_MPI
			CopyAllPolymers(polymer,replica,parms->npol,parms->nosidechains,parms->nosidechains);
			#endif
		}
		else if (irun==0 || parms->always_restart)
		{
			CopyAllPolymers(startp,polymer,parms->npol,parms->nosidechains,parms->nosidechains);
			#ifdef ACTIVE_MPI
			CopyAllPolymers(startp,replica,parms->npol,parms->nosidechains,parms->nosidechains);
			#endif
		}
		polymer->etot = TotalEnergy(polymer,u,parms,parms->npol,1,parms->nosidechains,parms->debug,my_rank);
		#ifdef ACTIVE_MPI
		replica->etot = TotalEnergy(replica,u,parms,parms->npol,1,parms->nosidechains,parms->debug,my_rank);
		#endif
		if (polymer->etot >= LARGE)
			fprintf(stderr,"WARNING: the structure has overlapping atoms (E=%g); the energy reported during the run keeps\n"
				"         their hard-core terms also after they are resolved (consider the keyword disentangle)\n",polymer->etot);
		if(my_rank==0)
		{
			fprintf(stderr,"\nPolymer energy=%lf\tstep=%d\n",polymer->etot,irun);
			#ifdef ACTIVE_MPI
			fprintf(stderr,"Replica energy=%lf\tstep=%d\n\n",replica->etot,irun);
			#endif
		}

		if (parms->shell==1)
		{
			
			SetTrueShellFlags(polymer,parms);
			UpdateShell(polymer,parms);
			SetFalseShellFlags(polymer,parms);

		}
		#ifdef ACTIVE_MPI
		if (parms->shell==1) 	UpdateShell(replica,parms);
		#endif

		#ifdef ACTIVE_MPI
		if (parms->nrun>1) fprintf(stderr,"> rank %d Starting run %d\n",my_rank,irun);
		#else
		if (parms->nrun>1) fprintf(stderr,"Starting run %d\n",irun);
		#endif
		#ifdef ACTIVE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		if(my_rank==0)
			fprintf(stderr,"Initial energy = %lf\n",polymer->etot);

		#endif
		
		// open trajectory file (name without the extension, which depends on trjformat)
		#ifdef ACTIVE_MPI
		if (parms->nrun>1)
			snprintf(fname,sizeof(fname),"%s_run%d_proc%d",parms->fntrj,irun,my_rank);
		else
			snprintf(fname,sizeof(fname),"%s_proc%d",parms->fntrj,my_rank);
		#else
		if (parms->nrun>1)
			snprintf(fname,sizeof(fname),"%s_run%d",parms->fntrj,irun);
		else
			snprintf(fname,sizeof(fname),"%s",parms->fntrj);
		#endif
		ftrj = OpenTrajectory(fname,polymer,parms,u,ckp.restart);	// a restart continues the files of the run
		if (!ftrj)
			Error("Cannot open trajectory file for writing");
		if (!ckp.restart)
		{
			snprintf(polymer->title,sizeof(polymer->title),"initial\tE=%lf",polymer->etot);
			PrintTrajectory(ftrj,polymer,parms,u);
		}

		// open energy file
		if (strcmp(parms->fne,""))
		{
			#ifdef ACTIVE_MPI
			if (parms->nrun>1)
				snprintf(fname,sizeof(fname),"%s_run%d_proc%d.ene",parms->fne,irun,my_rank);
			else
				snprintf(fname,sizeof(fname),"%s_proc%d.ene",parms->fne,my_rank);
			#else
			if (parms->nrun>1)
				snprintf(fname,sizeof(fname),"%s_run%d.ene",parms->fne,irun);
			else
				snprintf(fname,sizeof(fname),"%s.ene",parms->fne);
			#endif
			fe = fopen(fname,ckp.restart ? "a" : "w");
			if (!fe) 
				Error("Cannot open energy file for writing");
			if (!ckp.restart) PrintEnergies(fe,parms->npol,0,polymer);
		}

		// open the file of the observables of aggregation
		parms->fobs = NULL;
		if (parms->nprintobs > 0)
		{
			#ifdef ACTIVE_MPI
			if (parms->nrun>1) snprintf(fname,sizeof(fname),"obs_run%d_proc%d.dat",irun,my_rank);
			else snprintf(fname,sizeof(fname),"obs_proc%d.dat",my_rank);
			#else
			if (parms->nrun>1) snprintf(fname,sizeof(fname),"obs_run%d.dat",irun);
			else snprintf(fname,sizeof(fname),"obs.dat");
			#endif
			parms->fobs = fopen(fname,ckp.restart ? "a" : "w");
			if (!parms->fobs) Error("Cannot open the file of observables");
			if (!ckp.restart) PrintObservablesHeader(parms->fobs);
		}

		//if(parms->npol>1)
		//{
			//fprintf(stderr,"\nChecking number of Contacts\tstep=0\n");
			//CountContacts(stderr,polymer,parms,0);
		//}
		fprintf(stderr,"\n");

		#ifdef ACTIVE_MPI
		if(parms->shell==1)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			//fprintf(stderr,"rank %d\tnshell=%d\n",my_rank,parms->nshell);
			//SET LOCAL MOVE LIMIT ACCORDING TO iT_bias
			MPI_Barrier(MPI_COMM_WORLD);

			if(parms->movetype[7]>0)
			if(my_rank>=parms->iT_bias)
			{
				parms->movetype[7]=0;
			}

			//fprintf(stderr,"rank %d\tbgs frequency= %d\n",my_rank,parms->movetype[7]);
			MPI_Barrier(MPI_COMM_WORLD);
		
		}

		#endif
	
		// MAKE MONTE CARLO
		#ifdef ACTIVE_MPI		
		Do_MC(polymer,fragment,replica,startp,u,parms,ftrj,fe,oldp,fproc,irun,&ckp,mpiparms);
		#else
		Do_MC(polymer,fragment,replica,startp,u,parms,ftrj,fe,oldp,fproc,irun,&ckp,0);
		#endif
		ckp.restart = 0;						// the following runs start normally

		// check the energy before the optimisation of the potential changes it
		CheckOverlaps(polymer,u,parms,parms->npol,parms->nosidechains,1,fproc);
		fprintf(fproc,"Final: Etot=%lf\tEtot(true)=%lf\n",polymer->etot,TotalEnergy(polymer,u,parms,parms->npol,0,parms->nosidechains,parms->debug,my_rank));
		#ifdef OPTIMIZEPOT
		if (strcmp(parms->op_minim,"none"))
		{
			#ifdef ACTIVE_MPI
            			
			if (parms->op_allreplicas) OP_GatherFrames(polymer->op,parms,my_rank,nprocs);
			if(my_rank==0)
			OP_SamplePotential(polymer,parms,ntypes,u->e,irun,ematold);
			u = send_pot(nat,ntypes,parms->noangpot,parms->nodihpot,parms->nohfields,parms->hb, my_rank, nprocs, mpiparms->Pot_mpi, u, astatus);
			MPI_Barrier(MPI_COMM_WORLD);

			#else
                  		
			OP_SamplePotential(polymer,parms,ntypes,u->e,irun,ematold);
            #endif
			if(my_rank==0)
			{
				snprintf(fname,sizeof(fname),"newpot_%d.pot",irun);
				PrintPotential(u,fname,nat,ntypes,parms->noangpot,parms->nodihpot,parms->hb,parms->nohfields);
			}
		}
		#endif



		// open last conformation file
		#ifdef ACTIVE_MPI
             	
		if (parms->nrun>1) 
			snprintf(fname,sizeof(fname),"%s_run%d_proc%d.pol",parms->flastp,irun,my_rank);
		else
			snprintf(fname,sizeof(fname),"%s_proc%d.pol",parms->flastp,my_rank);
            	
		#else
             
		if (parms->nrun>1)
			snprintf(fname,sizeof(fname),"%s_run%d.pol",parms->flastp,irun);
		else
		{
			snprintf(fname,sizeof(fname),"%s.pol",parms->flastp);
		}
		#endif




		PrintPolymer(fname,polymer,parms->npol);

		// close files
		if (strcmp(parms->fntrj,"")) fclose(ftrj);
		if (strcmp(parms->fne,"")) fclose(fe);
		if (parms->fobs) { fclose(parms->fobs); parms->fobs = NULL; }

	
		if (my_rank==0) MGInfo_EndRun(parms,irun);
		#ifdef ACTIVE_MPI
            	if (my_rank==0) fprintf(stderr,"End of iteration %d\n",irun);
		#else
		fprintf(stderr,"End of iteration %d\n",irun);
            	#endif

		#ifdef STEMPERING
        if (parms->stempering)
        FreeStempering((parms->p));
        else free(parms->p);

		if(irun<((parms->nrun)-1))
			ResetStStuff(parms,argv[3]);
		#endif
	}
	
	#ifdef OPTIMIZEPOT
	fclose(fp);
	#endif

	
	CellsFree(polymer);
	FreePolymer(polymer,parms->npol, nbackmax, nsidemax, parms->shell, parms->nosidechains);
	FreePolymer(oldp,parms->npol, nbackmax, nsidemax, parms->shell, parms->nosidechains);	
	
	//FreeOpt(startp->op,parms,startp->op->ndata);	// commentato da bob
	if (strcmp(parms->op_minim,"none"))
	FreeOpt(startp->op,parms,parms->op_input->ndata);
	else free(startp->op);
	
	FreeTables(startp->tables);	
	FreePolymer(startp,parms->npol, nbackmax, nsidemax, parms->shell, parms->nosidechains);
	FreePolymer(replica,parms->npol, nbackmax, nsidemax, parms->shell, parms->nosidechains);
	FreeDoubleMatrix(ematold,ntypes);
	FreePotential(u,ntypes,parms->noangpot,parms->nodihpot,parms->hb);
     	
	if( parms->movetype[7]!=-1)
	{
	  FreeFragment(fragment,parms);
	}

	if(my_rank==0 && strcmp(parms->op_minim,"none"))	
	FreeOptInput(parms->op_input);	

	free(parms);
	
	if (my_rank==0) fprintf(stderr,"\n\n\"Monte Grappa, tu sei la mia patria,\nsovra te il nostro sole risplende\"\n\n");

	
	#ifdef ACTIVE_MPI
	time(&rawtime);
	fprintf(fproc,"Run finished at %s\n\n",asctime(localtime(&rawtime)));
	MPI_Type_free(&mpiparms->Pot_mpi);
	MPI_Type_free(&mpiparms->Rot_mpi);
	MPI_Type_free(&mpiparms->Side_mpi);
	MPI_Type_free(&mpiparms->Back_mpi);


	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	#endif


	return 0;
}

void Welcome(FILE *fp)
{
	time_t t = time(NULL);
	struct tm *tm = localtime(&t);


	fprintf(fp,"\n\n");
	fprintf(fp,"             ,/k.\n");
	fprintf(fp,"            /  ih,\t Monte^Grappa \n");     
	fprintf(fp,"       ,-' ,  `:7b \t\t  %d.%d \n",NVER,NSUBVER);
	fprintf(fp,"     _.-/   '  /b.`.4p,\n");
	fprintf(fp,"  --   ,    ,-' ^6x, `.'^=._\n");
	fprintf(fp,"\n");
	fprintf(fp,"\nG. Tiana, 2016\n");
	fprintf(fp,"\n");
	fprintf(fp,"%s", asctime(tm));
	fprintf(fp,"pid = %d\n",getpid());
	fprintf(fp,"\nExtensions:\n");
	#ifdef ACTIVE_MPI
	fprintf(fp,"+ [MPI]\t\t is ON\n");
	#else
	fprintf(fp,"- [MPI]\t\t is OFF\n");
	#endif
	#ifdef ACTIVE_OPENMP
	fprintf(fp,"+ [OPENMP]\t is ON\t(experimental!)\n");
	#else
	fprintf(fp,"- [OPENMP]\t is OFF\n");
	#endif
	#ifdef STEMPERING
	fprintf(fp,"+ [STEMPERING]\t is ON\n");
	#else
	fprintf(fp,"- [STEMPERING]\t is OFF\n");
	#endif

	fprintf(fp,"\n\n");



}

