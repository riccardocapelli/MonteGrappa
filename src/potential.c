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
    making use of MonteGrappa the following article:

    G. Tiana, F. Villa, Y. Zhan, R. Capelli, C.Paissoni, P. Sormanni, E. Heard,
    L. Giorgetti and R. Meloni, "MonteGrappa: an iterative Monte Carlo program
    to optimize biomolecular potentials in simplified models" Comp. Phys. Comm.
    Volume 186, January 2015, Pages 93–104.     DOI: 10.1016/j.cpc.2014.09.012

   
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * potential.c
 *
 *  Created on: Sep 23, 2010
 *      Author: guido
 */

#include "montegrappa.h"

unsigned long long contacts_stamp = 1;	// stamp of the current move for the contact lists it changes

// residues (chain, index) whose contact list changed since the last ContactsNewMove, each once
static int *touched_c=NULL,*touched_i=NULL,ntouched=0,touchedmax=0;

static inline void Touch(struct s_back *b, int c, int i)
{
	if (b->touched == contacts_stamp) return;
	b->touched = contacts_stamp;
	if (ntouched >= touchedmax)
	{
		touchedmax = 2*touchedmax + 1024;
		touched_c = (int *) realloc(touched_c,touchedmax*sizeof(int));
		touched_i = (int *) realloc(touched_i,touchedmax*sizeof(int));
		if (!touched_c || !touched_i) Error("Cannot allocate the list of changed contacts");
	}
	touched_c[ntouched] = c;
	touched_i[ntouched] = i;
	ntouched ++;
}

/********************************************************************
 Residues whose contact list changed in the current move (AddContact,
 ResetContactsMonomer and TotalEnergy record them): the moves copy
 only these between p and oldp (UpdateChainRanges)
 ********************************************************************/
int ContactsTouched(int **c, int **i)
{
	*c = touched_c;
	*i = touched_i;
	return ntouched;
}

void ContactsNewMove(void)
{
	ntouched = 0;
	contacts_stamp ++;
}

/********************************************************************
 Largest distance at which two atoms interact (well or hard core):
 EnergyPair skips pairs of residues farther than this plus the
 extent of their sidechains
 ********************************************************************/
void SetPairCutoff(struct s_potential *u, int ntypes)
{
	int i,j;
	double r2max = u->g_r0hard * u->g_r0hard;

	for (i=0;i<ntypes;++i)
		for (j=0;j<ntypes;++j)
		{
			if ((u->r_2)[i][j] > r2max) r2max = (u->r_2)[i][j];
			if ((u->r0_2)[i][j] > r2max) r2max = (u->r0_2)[i][j];
		}
	u->rcut = sqrt(r2max);
	SetPairParameters(u,ntypes);
}

/********************************************************************
 Copy of r0_2, r_2 and e of each pair of types in a single array, read
 by the pair energy of atoms: to be called again whenever the matrices
 change (optimization of the potential)
 ********************************************************************/
void SetPairParameters(struct s_potential *u, int ntypes)
{
	int i,j;

	if (!u->pp || u->npptypes != ntypes)
	{
		free(u->pp);
		u->pp = (struct s_pairpar *) malloc((size_t)ntypes*ntypes*sizeof(struct s_pairpar));
		if (!u->pp) Error("Cannot allocate the parameters of the pairs");
		u->npptypes = ntypes;
	}
	for (i=0;i<ntypes;++i)
		for (j=0;j<ntypes;++j)
		{
			u->pp[i*ntypes+j].r02 = (u->r0_2)[i][j];
			u->pp[i*ntypes+j].rw2 = (u->r_2)[i][j];
			u->pp[i*ntypes+j].e = (u->e)[i][j];
		}
}

/********************************************************************
 The angular factor of a hydrogen bond uses the atom hb_iam preceding
 the donor/acceptor: if it belongs to another backbone unit, a move of
 that unit changes the bond energy without recomputing it. Returns the
 number of such donors/acceptors.
 ********************************************************************/
int CheckHBReferences(struct s_polymer *p, struct s_potential *u, int npol, int nat)
{
	int ic,i,is,k,g=0,nbad=0,*owner;

	owner = (int *) malloc(nat*sizeof(int));
	if (!owner) Error("Cannot allocate in CheckHBReferences");
	for (k=0;k<nat;++k) owner[k] = -1;
	for (ic=0;ic<npol;++ic)
		for (i=0;i<(p+ic)->nback;++i,++g)
		{
			owner[(((p+ic)->back)+i)->ia] = g;
			for (is=0;is<(((p+ic)->back)+i)->nside;++is)
				owner[(((((p+ic)->back)+i)->side)+is)->ia] = g;
		}
	for (k=0;k<nat;++k)
		if ((u->hb)[k] && (u->hb_iam[k]<0 || u->hb_iam[k]>=nat || owner[u->hb_iam[k]] != owner[k]))
		{
			if (nbad < 5) fprintf(stderr,"WARNING: hydrogen bond atom %d refers to atom %d, which is not in its backbone unit\n",k,u->hb_iam[k]);
			nbad++;
		}
	free(owner);
	return nbad;
}

/********************************************************************
 Largest distance of the sidechain atoms of each residue from its
 backbone atom (also kept up to date by AddSidechain)
 ********************************************************************/
void SetSidechainExtents(struct s_polymer *p, int npol, int nosidechains)
{
	int ci,i,is;
	double e2,e2max;

	for (ci=0;ci<npol;++ci)
		for (i=0;i<(p+ci)->nback;++i)
		{
			e2max = 0;
			if (!nosidechains)
				for (is=0;is<(((p+ci)->back)+i)->nside;++is)
				{
					e2 = Dist2( (((((p+ci)->back)+i)->side)+is)->pos, (((p+ci)->back)+i)->pos );
					if (e2 > e2max) e2max = e2;
				}
			(((p+ci)->back)+i)->ext = sqrt(e2max);
		}
}

/********************************************************************
 Pair energy of the whole system with the cell list: only residue pairs
 in neighbouring cells, summed in the same order as the plain loops of
 TotalEnergy (chain ci, chain cj>=ci, residue i, residue j), so that the
 result is identical. Returns 0 (nothing done) if some residue has too
 many candidates for the cell list to be useful.
 ********************************************************************/
static int TotalPairEnergyCells(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int update, int nosidechains, double *etot)
{
	struct s_cells *c = p->cells;
	int g,n,k,ci,cj,i,j,tooclose,*start,*list,*ptr,nlist=0,size;

	start = (int *) malloc((c->nres+1)*sizeof(int));
	ptr = (int *) malloc(c->nres*sizeof(int));
	size = 32*c->nres;
	list = (int *) malloc(size*sizeof(int));
	if (!start || !ptr || !list) Error("Cannot allocate in TotalPairEnergyCells");

	// sorted candidates of every residue
	for (g=0;g<c->nres;++g)
	{
		start[g] = nlist;
		n = CellsCandidates(p,c->chain[g],c->res[g]);
		if (n < 0) { free(start); free(ptr); free(list); return 0; }
		if (nlist+n > size)
		{
			size = 2*(nlist+n);
			list = (int *) realloc(list,size*sizeof(int));
			if (!list) Error("Cannot allocate in TotalPairEnergyCells");
		}
		for (k=0;k<n;++k) list[nlist++] = c->cand[k];
	}
	start[c->nres] = nlist;

	for (ci=0;ci<npol;++ci)
	{
		for (i=0;i<(p+ci)->nback;++i) ptr[c->offset[ci]+i] = start[c->offset[ci]+i];
		for (cj=ci;cj<npol;++cj)
			for (i=0;i<(p+ci)->nback;++i)
			{
				g = c->offset[ci]+i;
				// candidates are sorted by global index: those of chain cj are a contiguous block
				while (ptr[g] < start[g+1] && list[ptr[g]] < c->offset[cj]) ptr[g]++;
				for (k=ptr[g]; k<start[g+1] && list[k] < c->offset[cj+1]; ++k)
				{
					j = list[k] - c->offset[cj];
					if (ci!=cj || i<j)					// in the same polymer, do not count twice a contact
					{
						tooclose = 0;
						if (ci==cj && j-i <= u->g_imin ) tooclose = 1;
						*etot += EnergyPair(p,u,i,j,ci,cj,update,nosidechains,parms->disentangle,tooclose,parms->hb);
					}
				}
			}
	}
	free(start); free(ptr); free(list);
	return 1;
}

double TotalEnergy(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int update, int nosidechains, int debug, int iproc)
{
	int i,j,ci,cj,tooclose;
	double etot=0,eang=0,edih=0,ehfs=0,ebox=0;

	SetSidechainExtents(p,npol,nosidechains);			// coordinates may have been loaded in bulk (start, exchange, restart)
	if (p->cells) CellsBuild(p,u->rcut,u->boxtype=='p' ? u->lbox : 0.);

	if (update)
		#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
		#endif
		  for (ci=0;ci<npol;++ci)
			for (i=0;i<(p+ci)->nback;++i)
                  {
                        for (j=0;j<(((p+ci)->back)+i)->ncontacts;++j)  ((((p+ci)->back)+i)->e)[j] = 0;
                        (((p+ci)->back)+i)->ncontacts = 0;
                  }
	if (update)								// (not in the parallel loop: the list is shared)
		for (ci=0;ci<npol;++ci)
			for (i=0;i<(p+ci)->nback;++i) Touch(((p+ci)->back)+i,ci,i);

	if (p->cells && p->cells->active && TotalPairEnergyCells(p,u,parms,npol,update,nosidechains,&etot))
		;							// pairs from the cell list
	else if(update)
	{
	for (ci=0;ci<npol;++ci)                                     // loops on polymers
            for (cj=ci;cj<npol;++cj)
                  for (i=0;i<(p+ci)->nback;++i)             // loops on backbone atoms of each polymer
                        for (j=0;j<(p+cj)->nback;++j)
                              if (ci!=cj || i<j)                        // in the same polymer, do not count twice a contact
                              {
                                    tooclose = 0;
                                    if (ci==cj && j-i <= u->g_imin ) tooclose = 1;
                                    etot += EnergyPair(p,u,i,j,ci,cj,update,nosidechains,parms->disentangle,tooclose,parms->hb); 
                             }
	}
	else
	{
	#ifdef _OPENMP
	#pragma omp parallel for private(cj,i,j,tooclose) reduction(+:etot)
	#endif
	for (ci=0;ci<npol;++ci)                                     // loops on polymers
            for (cj=ci;cj<npol;++cj)
                  for (i=0;i<(p+ci)->nback;++i)             // loops on backbone atoms of each polymer
                        for (j=0;j<(p+cj)->nback;++j)
                              if (ci!=cj || i<j)                        // in the same polymer, do not count twice a contact
                              {
                                    tooclose = 0;
                                    if (ci==cj && j-i <= u->g_imin ) tooclose = 1;
                                    etot += EnergyPair(p,u,i,j,ci,cj,update,nosidechains,parms->disentangle,tooclose,parms->hb);
                             }
	}
	if (debug>0 && iproc==0) fprintf(parms->flog," Epairs = %lf\n",etot);

	if (!parms->noangpot || !parms->nodihpot || !parms->nohfields || u->boxtype != 'n')
	{
		if(update)
		{
			for (ci=0;ci<npol;++ci)
			{
				for (i=0;i<(p+ci)->nback;++i)
				{
					if (!parms->noangpot){
   
					eang += EnergyAngles(p,u,i,ci,update);                      // angular potential

					}
					if (!parms->nodihpot){
						edih += EnergyDihedrals(p,u,i,ci,update);					// dihedral potential
 					}
					if ( (!parms->nohfields) ){
						ehfs += EnergyHFields(p,u,parms,i,ci,update);					// H fields potential
					}
					if (u->boxtype != 'n')                                                  // box potential
						if (EnergyBox(p,u,i,ci) == 1) ebox = LARGE;
				}
			}
		}
		else
		{
		#ifdef _OPENMP
		#pragma omp parallel for shared(ebox) private(i) reduction(+:eang,edih,ehfs)
		#endif
			for (ci=0;ci<npol;++ci)
			{
				for (i=0;i<(p+ci)->nback;++i)
				{
					if (!parms->noangpot)
						eang += EnergyAngles(p,u,i,ci,update);                  // angular potential
					if (!parms->nodihpot)
						edih += EnergyDihedrals(p,u,i,ci,update);				// dihedral potential
					if ( (!parms->nohfields) )
						ehfs += EnergyHFields(p,u,parms,i,ci,update);					// H fields potential
					if (u->boxtype != 'n')
						if (EnergyBox(p,u,i,ci) == 1) ebox = LARGE;				// box potential
				}
			}
		}
	
	}
		
	etot += (eang + edih + ebox + ehfs);

	if (debug>0 && iproc==0) fprintf(parms->flog," Eang  = %lf\n Edih  = %lf\n Ehfs  = %lf\n Etot  = %lf\n",eang,edih,ehfs,etot);


      #ifdef OPTIMIZEPOT
       if (strcmp(parms->op_minim,"none"))
       {
            p->op->eold[ p->op->nframes ] = etot;
            p->op->efix[ p->op->nframes ] = eang + edih + ebox + ehfs;
       }
      #endif


      return etot;
}




/********************************************************************
 Hard core and square well of two atoms of types a1, a2 at square
 distance r2 (a hydrogen bond multiplies the well by hbmul, -1 if none).
 Returns 1 if the hard core is hit without disentangle (the pair
 energy is then LARGE).
 ********************************************************************/
static inline int AtomPairEnergy(struct s_polymer *p, struct s_potential *u, int a1, int a2, double r2, int tooclose, int disentangle,
		int record, double *etot, double *econt, int *incontact)
{
	double e,mul,rw2;
	const struct s_pairpar *q = u->pp + a1*u->npptypes + a2;

	if ( r2 < (tooclose ? u->g_r0hard * u->g_r0hard : q->r02) )
	{
		if (!disentangle) return 1;
		else *etot += LARGE;
	}
	if (tooclose) return 0;
	rw2 = q->rw2;
	if ( r2 < rw2 )							// close in space but not in sequence: the element is defined (rw2>0)
	{
		e = q->e;
		mul = 1.;
		if (u->splice==1)
			if (rw2 * u->kr2_splice < r2) mul = u->ke_splice;

		#ifdef OPTIMIZEPOT
		if (record && mul>0) OP_AddEnergy(p,a1,a2,mul);
		#endif

		e *= mul;
		*etot += e;
		*econt += e;
		*incontact = 1;
	}
	return 0;
}

/********************************************************************
 Calculate the energy of a pair of residues,
 updating contact lists and partial energies if update==1.
 ********************************************************************/
static double EnergyPairAtoms(struct s_polymer *p, struct s_potential *u, int i, int j, int ci, int cj, int nosidechains, int disentangle, int tooclose, int hb, double *econt, int *incontact)
{
	int is,js,a1,a2,ia1,ia2,out,record=0;
	double r2,etot=0,e,r02,cos1,cos2,mul,rmax,r2bb;
	struct vector sh = {0.,0.,0.},v2,v2m,pj,psj;			// periodic box: shift of the atoms of j to the image nearest to i
	struct s_back *bi = ((p+ci)->back)+i, *bj = ((p+cj)->back)+j;
	struct s_side *si,*sj;

	// within same aa do not interact
	if (ci==cj && bi->iaa == bj->iaa) return 0;

	// no atom of the two residues can be within the interaction range
	if (u->boxtype=='p') sh = ImageShift( bi->pos, bj->pos, u->lbox );
	r2bb = Dist2S( bi->pos, bj->pos, sh );
	rmax = u->rcut + bi->ext + bj->ext;
	if (r2bb > rmax*rmax) return 0;

	#ifdef OPTIMIZEPOT
	record = p->op->record;
	#endif
	// the atoms of j in the image of i (as in Dist2S: b+s-a)
	pj.x = bj->pos.x + sh.x; pj.y = bj->pos.y + sh.y; pj.z = bj->pos.z + sh.z;

	// interaction backbone-backbone (i-j)
	if ( ci!=cj || (i!=(j+1) && i!=(j-1)) )						// if not consecutive (i.e., no covalent interaction)
		if (AtomPairEnergy(p,u,bi->itype,bj->itype,r2bb,tooclose,disentangle,record,&etot,econt,incontact)) return LARGE;

  	if (nosidechains == 1) return etot;

	for (js=0;js<bj->nside;++js)
	{
		sj = bj->side+js;
		psj.x = sj->pos.x + sh.x; psj.y = sj->pos.y + sh.y; psj.z = sj->pos.z + sh.z;
		a2 = sj->itype;

		// interaction backbone-sidechain (i-js)
		r2 = (psj.x-bi->pos.x)*(psj.x-bi->pos.x) + (psj.y-bi->pos.y)*(psj.y-bi->pos.y) + (psj.z-bi->pos.z)*(psj.z-bi->pos.z);
		if (AtomPairEnergy(p,u,bi->itype,a2,r2,tooclose,disentangle,record,&etot,econt,incontact)) return LARGE;

		for (is=0;is<bi->nside;++is)
		{
			// interaction sidechain-sidechain (is-js)
			si = bi->side+is;
			a1 = si->itype;
			r2 = (psj.x-si->pos.x)*(psj.x-si->pos.x) + (psj.y-si->pos.y)*(psj.y-si->pos.y) + (psj.z-si->pos.z)*(psj.z-si->pos.z);

			if (!hb)
			{
				if (AtomPairEnergy(p,u,a1,a2,r2,tooclose,disentangle,record,&etot,econt,incontact)) return LARGE;
				continue;
			}

			if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
			else r02 = (u->r0_2)[a1][a2];

			if ( r2 < r02 )
			{
				if (!disentangle) return LARGE;
				else etot += LARGE;
			}
			if ( r2 < (u->r_2)[a1][a2] && !tooclose)				// if they are close in space anbd not in sequence (the element is defined)
			{
				e = (u->e)[a1][a2];
				mul = 1;
				if (u->splice==1)
					if ((u->r_2)[a1][a2] * u->kr2_splice < r2) mul =  u->ke_splice;

				// calculate hydrogen bond (if any)
				ia1 = si->ia;
				ia2 = sj->ia;
				if ((u->hb)[ia1]*(u->hb)[ia2]==2)		// donor-acceptor pair
				{
				 // a hydrogen bond acts only in the first half of splice (other pairs keep ke_splice)
				 if (u->splice==1 && (u->r_2)[a1][a2] * u->kr2_splice < r2) mul = 0;
				 else
				 {
					v2 = *(*(((p+cj)->vback)+ia2));
					v2.x += sh.x; v2.y += sh.y; v2.z += sh.z;
					v2m = *(*(((p+cj)->vback)+(u->hb_iam)[ia2]));
					v2m.x += sh.x; v2m.y += sh.y; v2m.z += sh.z;
					cos1 = CosAngle( (*(*(((p+ci)->vback)+(u->hb_iam)[ia1]))), (*(*(((p+ci)->vback)+ia1))), 
							v2, p->tables, &out );
					if (out==0) return LARGE;
					cos2 = CosAngle( (*(*(((p+ci)->vback)+ia1))), v2, v2m, p->tables , &out);
					if (out==0) return LARGE;
					mul = FastSqrt( DAbs(cos1*cos2),p->tables);
				 }
				}

				#ifdef OPTIMIZEPOT
				if (record && mul>0) OP_AddEnergy(p,a1,a2,mul);
				#endif
				e *= mul;
				etot += e;
				{ *econt += e; *incontact = 1; }
			}
		}
	}

	// interaction backbone-sidechain (is-j)
	for (is=0;is<bi->nside;++is)
	{
		si = bi->side+is;
		r2 = (pj.x-si->pos.x)*(pj.x-si->pos.x) + (pj.y-si->pos.y)*(pj.y-si->pos.y) + (pj.z-si->pos.z)*(pj.z-si->pos.z);
		if (AtomPairEnergy(p,u,bj->itype,si->itype,r2,tooclose,disentangle,record,&etot,econt,incontact)) return LARGE;
	}
	return etot;
}

/********************************************************************
 Interaction energy of residues i (chain ci) and j (chain cj); if
 update, the pair is recorded once in their contact lists (also when
 a hard core is found, as the move is then rejected and restored)
 ********************************************************************/
double EnergyPair(struct s_polymer *p, struct s_potential *u, int i, int j, int ci, int cj, int update, int nosidechains, int disentangle, int tooclose, int hb)
{
	double e,econt=0;
	int incontact=0;

	e = EnergyPairAtoms(p,u,i,j,ci,cj,nosidechains,disentangle,tooclose,hb,&econt,&incontact);
	if (update && incontact) AddContact(p,i,j,ci,cj,econt);
	return e;
}

/********************************************************************
 Calculate the energy of ith residue with all the others,
 updating contact lists and partial energies. (single loop)
 ********************************************************************/
double EnergyMonomer(struct s_polymer *p, struct s_potential *u, int i, int ci, int npol, int update, int shell, int nosidechains,int disentangle,int hb)
{
	int j,cj,k,tooclose;
	double e,etot=0;

/*	// box energy
	if (u->boxtype != 'n')
		if (EnergyBox(p,u,i,ci) == 1) return LARGE;
*/
	// delete contact list and partial energies of i and of its contacting monomers
	if ( update ) ResetContactsMonomer(p,i,ci);


	// box energy
	if (u->boxtype != 'n')
	if (EnergyBox(p,u,i,ci) == 1) return LARGE;


	if (shell == 0)
	{
		for (cj=0;cj<npol;++cj)
		{
			for (j=0;j<(p+cj)->nback;++j)
			{
				if ( ci!= cj || i!=j )
				{
					tooclose = 0;
					if (ci == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
					e = EnergyPair(p,u,i,j,ci,cj,update,nosidechains,disentangle,tooclose,hb);

					if (e>=LARGE && !disentangle) return LARGE;
					etot += e;
				}
			}
		}
	}

	// use shell
	else
	{
		for (k=0;k<(((p+ci)->back)+i)->nshell;++k)
		{
			j = ((((p+ci)->back)+i)->shell)[k];									// what is the monomer in the shell
			cj = ((((p+ci)->back)+i)->shell_p)[k];	
	//		if(i==26 && ci == 0 && j == 37 && cj == 0) 							// ... what polymers it belongs to
	//		fprintf(stderr,"interaction atom %d of %d // %d of %d \n",i,ci,j,cj);	
			//SSSS
			if( ci!=cj || i!=j)
			{
         			//if(i==26 && ci == 0 && j == 37 && cj == 0)     fprintf(stderr,"interaction atom %d of %d // %d of %d \n",i,ci,j,cj);

				tooclose=0;
				if (ci == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
				e = EnergyPair(p,u,i,j,ci,cj,update,nosidechains,disentangle,tooclose,hb);
				//if(i==136 && ci == 0 && j == 151 && cj == 0)  
				//fprintf(stderr,"SISHELL interaction atom %d of %d // %d of %d E=%lf\tirot1=%d\tirot2=\%d\n",i,ci,j,cj,e,(((p+ci)->back)+i)->irot,(((p+cj)->back)+j)->irot);
 

				if (e>=LARGE && !disentangle) return LARGE;
				etot += e;
			}
		}
	}

	return etot;
}

/********************************************************************
 Residues g0..g1-1 of the compact copy whose sphere (CellsSphere) is
 within reach of the sphere of centre (xi,yi,zi) and radius plus range ri
 (the others cannot interact): appended to c->cand from position nhit
 (branchless, in their order). Returns the new number of candidates.
 ********************************************************************/
static inline int CandidatesRange(struct s_cells *c, int g0, int g1, double xi, double yi, double zi, double ri, struct s_potential *u, int nhit)
{
	double dx,dy,dz,rmax;
	int g;

	if (u->boxtype=='p')
	{
		double box=u->lbox;
		for (g=g0;g<g1;++g)
		{
			dx = c->sx[g]-xi; dy = c->sy[g]-yi; dz = c->sz[g]-zi;
			dx -= box*round(dx/box); dy -= box*round(dy/box); dz -= box*round(dz/box);
			rmax = ri + c->srad[g];
			c->cand[nhit] = g;
			nhit += ( dx*dx + dy*dy + dz*dz <= rmax*rmax );
		}
	}
	else
		for (g=g0;g<g1;++g)
		{
			dx = c->sx[g]-xi; dy = c->sy[g]-yi; dz = c->sz[g]-zi;
			rmax = ri + c->srad[g];
			c->cand[nhit] = g;
			nhit += ( dx*dx + dy*dy + dz*dz <= rmax*rmax );
		}
	return nhit;
}

/********************************************************************
 Calculate the energy of energies of monomers in a range, updating
 contacts and partial energies
 ********************************************************************/
double EnergyMonomerRange(struct s_polymer *p, struct s_potential *u, int from, int to, int ip, int npol, int shell, int update, int nosidechains, int disentangle, int hb)
{
	return EnergyMonomerRangeRigid(p,u,from,to,ip,npol,shell,update,nosidechains,disentangle,hb,NULL);
}

/********************************************************************
 The same, when the residues of chain ip (within [from,to]) with
 rig[i]=1 moved rigidly together: their pairs keep energy and contacts
 (rig=NULL computes all the pairs)
 ********************************************************************/
double EnergyMonomerRangeRigid(struct s_polymer *p, struct s_potential *u, int from, int to, int ip, int npol, int shell, int update, int nosidechains, int disentangle, int hb, const char *rig)
{
	int i,j,cj,k,tooclose;
	double etot=0,e;

	if (from<0) from=0;
	if (to>(p+ip)->nback - 1) to=(p+ip)->nback - 1;



	// delete contact list and partial energies of i and of its contacting monomers
	if ( update )
		//#pragma omp parallel for 
		for (i=from;i<=to;++i)
		{
			if (rig && rig[i]) ResetContactsMonomerOutside(p,i,ip,rig);
			else ResetContactsMonomer(p,i,ip);
		}

	// box energy
	if (u->boxtype != 'n')
        	for (i=from;i<=to;++i)
			if (EnergyBox(p,u,i,ip) == 1) return LARGE;

	if (shell == 0 && p->cells)
	{
		CellsUpdateRange(p,ip,from,to);
		if (p->cells->dirty) CellsBuild(p,u->rcut,u->boxtype=='p' ? u->lbox : 0.);
	}

	// if not shell, with cell list
	if (shell == 0 && p->cells && p->cells->active)
	{
		struct s_cells *c = p->cells;
		double rmax,dx,dy,dz;
		int ncand,k,g,gi;

		for (i=from;i<=to;++i)
		{
			gi = c->offset[ip]+i;
			ncand = CellsCandidates(p,ip,i);
			if (ncand<0)						// too many: all residues, in the same order
			{
				ncand = c->nres;
				for (k=0;k<ncand;++k) c->cand[k] = k;
			}
			for (k=0;k<ncand;++k)
			{
				g = c->cand[k];
				cj = c->chain[g];
				j = c->res[g];
				rmax = u->rcut + c->srad[gi] + c->srad[g];				// spheres of the residues (see CellsSphere)
				dx = c->sx[g]-c->sx[gi]; dy = c->sy[g]-c->sy[gi]; dz = c->sz[g]-c->sz[gi];
				if (c->box > 0) { dx -= c->box*round(dx/c->box); dy -= c->box*round(dy/c->box); dz -= c->box*round(dz/c->box); }
				if ( dx*dx + dy*dy + dz*dz > rmax*rmax ) continue;

				if ( ( ip!=cj || j<from || j>to || i<j ) && !( rig && ip==cj && rig[i] && rig[j] ) )	// j in [from,to] of the same chain: i<j to count once; skip rigid pairs
				{
					tooclose = 0;
					if (ip == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
					e = EnergyPair(p,u,i,j,ip,cj,update,nosidechains,disentangle,tooclose,hb);
					if (e>=LARGE && !disentangle) return LARGE;
					etot += e;
				}
			}
		}
	}

	// if not shell, all residues from the compact copy of the positions
	else if (shell == 0 && p->cells)
	{
		struct s_cells *c = p->cells;
		double xi,yi,zi,ri;
		int nhit,k,g,gi,r0=to+1,r1=from-1;

		if (rig)							// first and last rigid residue
			for (i=from;i<=to;++i)
				if (rig[i]) { if (i<r0) r0 = i; r1 = i; }

		for (i=from;i<=to;++i)
		{
			gi = c->offset[ip]+i;
			xi = c->sx[gi]; yi = c->sy[gi]; zi = c->sz[gi]; ri = u->rcut + c->srad[gi];

			// spheres of the residues (tighter than the test at the beginning of EnergyPair); branchless, the hits keep their order
			// (a rigid residue skips the rigid ones, whose pairs are kept)
			if (rig && rig[i])
			{
				nhit = CandidatesRange(c,0,c->offset[ip]+r0,xi,yi,zi,ri,u,0);
				for (g=c->offset[ip]+r0;g<=c->offset[ip]+r1;++g)
					if (!rig[g-c->offset[ip]]) nhit = CandidatesRange(c,g,g+1,xi,yi,zi,ri,u,nhit);
				nhit = CandidatesRange(c,c->offset[ip]+r1+1,c->nres,xi,yi,zi,ri,u,nhit);
			}
			else nhit = CandidatesRange(c,0,c->nres,xi,yi,zi,ri,u,0);

			for (k=0;k<nhit;++k)
			{
				cj = c->chain[c->cand[k]];
				j = c->res[c->cand[k]];
				if ( ( ip!=cj || j<from || j>to || i<j ) && !( rig && ip==cj && rig[i] && rig[j] ) )	// j in [from,to] of the same chain: i<j to count once; skip rigid pairs
				{
					tooclose = 0;
					if (ip == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
					e = EnergyPair(p,u,i,j,ip,cj,update,nosidechains,disentangle,tooclose,hb);
					if (e>=LARGE && !disentangle) return LARGE;
					etot += e;
				}
			}
		}
	}

	// if not shell
	else if (shell == 0)
	{
		struct s_back *bi,*bj;
		double rmax;

		for (i=from;i<=to;++i)
		{
			bi = ((p+ip)->back)+i;
			for (cj=0;cj<npol;++cj)
			{
				for (j=0;j<(p+cj)->nback;++j)
				{
					// same test as at the beginning of EnergyPair, to skip the call for far pairs
					bj = ((p+cj)->back)+j;
					rmax = u->rcut + bi->ext + bj->ext;
					if ( (bi->pos.x-bj->pos.x)*(bi->pos.x-bj->pos.x) + (bi->pos.y-bj->pos.y)*(bi->pos.y-bj->pos.y)
						+ (bi->pos.z-bj->pos.z)*(bi->pos.z-bj->pos.z) > rmax*rmax ) continue;

					if ( ( ip!=cj || j<from || j>to || i<j ) && !( rig && ip==cj && rig[i] && rig[j] ) )	// j in [from,to] of the same chain: i<j to count once; skip rigid pairs
					{
						//fprintf(stderr,"energypair atomi %d di %d / %d di %d \n",i,ip,j,cj);
						tooclose = 0;
						if (ip == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
						e = EnergyPair(p,u,i,j,ip,cj,update,nosidechains,disentangle,tooclose,hb);
						if (e>=LARGE && !disentangle) return LARGE;
						etot += e;
					}
				}
			}
		}
	}

	// if shell
	else
	{
		//slow #pragma omp parallel for private(k,i,j,cj,tooclose,e) reduction(+:etot)		
		for (i=from;i<=to;++i)
		{
			for (k=0;k<(((p+ip)->back)+i)->nshell;++k)
			{
				j = ((((p+ip)->back)+i)->shell)[k];							// what is the monomer in the shell
				cj = ((((p+ip)->back)+i)->shell_p)[k];						// ... what polymers it belongs to
				if ( ip!=cj || j<from || j>to || i<j )					// if j in [from,to] of the same chain, require i<j to count only once
				{
					
					tooclose = 0;
					if (ip == cj && Abs(i-j) <= u->g_imin) tooclose = 1;
					e = EnergyPair(p,u,i,j,ip,cj,update,nosidechains,disentangle,tooclose,hb);
					if (e>=LARGE && !disentangle) e= LARGE;
					etot += e;
				}
			}
		}
	}
	return etot;
}

/********************************************************************
 Tell backbones i and j that they are in contact; if they are already in
 contact, update the partial energy
 ********************************************************************/
void AddContact(struct s_polymer *p, int i, int j, int ci, int cj, double e)
{
	struct s_back *bi = ((p+ci)->back)+i;
	struct s_back *bj = ((p+cj)->back)+j;

	Touch(bi,ci,i);
	Touch(bj,cj,j);

	// each residue pair is evaluated once after the contacts of the moved residues
	// have been reset (ResetContactsMonomer, TotalEnergy), so it is not in the lists yet
	#ifdef DEBUG
	int k;
	for (k=0;k<bi->ncontacts;++k)
		if (bi->contacts[k]==j && bi->contacts_p[k]==cj) Error("AddContact: contact already present");
	#endif

	if (bi->ncontacts+1 >= NCONTMAX || bj->ncontacts+1 >= NCONTMAX) Error("NCONTMAX too small");

	bi->contacts[bi->ncontacts] = j;
	bi->contacts_p[bi->ncontacts] = cj;
	bi->contacts_r[bi->ncontacts] = bj->ncontacts;		// each entry knows where its partner is
	bi->e[bi->ncontacts] = e;

	bj->contacts[bj->ncontacts] = i;
	bj->contacts_p[bj->ncontacts] = ci;
	bj->contacts_r[bj->ncontacts] = bi->ncontacts;
	bj->e[bj->ncontacts] = e;

	bi->ncontacts ++;
	bj->ncontacts ++;
}

/********************************************************************
 Number of contact entries whose reverse index (contacts_r) does not
 point back to them (for echeck)
 ********************************************************************/
int ContactsReverseErrors(struct s_polymer *p, int npol)
{
	int ci,i,k,nerr=0;
	struct s_back *b,*m;

	for (ci=0;ci<npol;++ci)
		for (i=0;i<(p+ci)->nback;++i)
		{
			b = ((p+ci)->back)+i;
			for (k=0;k<b->ncontacts;++k)
			{
				m = ((p+b->contacts_p[k])->back)+b->contacts[k];
				if (b->contacts_r[k] < 0 || b->contacts_r[k] >= m->ncontacts || m->contacts[b->contacts_r[k]] != i ||
				    m->contacts_p[b->contacts_r[k]] != ci || m->contacts_r[b->contacts_r[k]] != k) nerr++;
			}
		}
	return nerr;
}

/********************************************************************
 Remove entry l of the contact list of b (residue i of chain c),
 moving the last one in its place (and telling its partner)
 ********************************************************************/
static inline void RemoveContactEntry(struct s_polymer *p, struct s_back *b, int c, int i, int l)
{
	int last = b->ncontacts-1;
	struct s_back *m;

	if (l != last)
	{
		b->contacts[l] = b->contacts[last];
		b->contacts_p[l] = b->contacts_p[last];
		b->contacts_r[l] = b->contacts_r[last];
		b->e[l] = b->e[last];
		m = ((p+b->contacts_p[l])->back)+b->contacts[l];
		m->contacts_r[b->contacts_r[l]] = l;
		Touch(m,b->contacts_p[l],b->contacts[l]);
	}
	b->ncontacts --;
	Touch(b,c,i);
}

/********************************************************************
 Tell i and j that they are in the same shell
 ********************************************************************/
void AddShell(struct s_polymer *p, int i, int j, int ci, int cj)
{
	int k,ok;

	// check if j is already among the shell of i
	ok=1;
	for (k=0;k<(((p+ci)->back)+i)->nshell;++k)
		if ( ((((p+ci)->back)+i)->shell)[k] == j && ((((p+ci)->back)+i)->shell_p)[k] == cj) ok=0;
	if (ok==1) // if not, add it
	{
		if ( (((p+ci)->back)+i)->nshell >= NSHELLMAX) Error("NSHELLMAX too small");
		((((p+ci)->back)+i)->shell)[ (((p+ci)->back)+i)->nshell ] = j;
		((((p+ci)->back)+i)->shell_p)[ (((p+ci)->back)+i)->nshell ] = cj;
		(((p+ci)->back)+i)->nshell ++;
	}

	// check if i is already among the shell of j
	ok=1;
	for (k=0;k<(((p+cj)->back)+j)->nshell;++k)
		if ( ((((p+cj)->back)+j)->shell)[k] == i && ((((p+cj)->back)+j)->shell_p)[k] == ci) ok=0;
	if (ok==1) // if not, add it
	{
		if ( (((p+cj)->back)+j)->nshell >= NSHELLMAX) Error("NSHELLMAX too small");
		((((p+cj)->back)+j)->shell)[ (((p+cj)->back)+j)->nshell ] = i;
		((((p+cj)->back)+j)->shell_p)[ (((p+cj)->back)+j)->nshell ] = ci;
		(((p+cj)->back)+j)->nshell ++;
	}
}

/********************************************************************
 Reset the contacts of i-th monomer, and all references to it of
 its contacting monomers (contact list and energies)
 ********************************************************************/
void ResetContactsMonomer(struct s_polymer *p, int i, int ci)
{
	int k,neigh,cn;
	struct s_back *bi = ((p+ci)->back)+i;

	// remove i from its neighbours (the entry of i in the list of the k-th is contacts_r[k])
	for (k=0;k<bi->ncontacts;++k)
	{
		neigh = bi->contacts[k];
		cn = bi->contacts_p[k];
		RemoveContactEntry(p,((p+cn)->back)+neigh,cn,neigh,bi->contacts_r[k]);
		bi->e[k] = 0;								// reset its energy
	}
	// reset the neighbour list of i
	if ((((p+ci)->back)+i)->ncontacts) Touch(((p+ci)->back)+i,ci,i);
	(((p+ci)->back)+i)->ncontacts = 0;
}

/********************************************************************
 As ResetContactsMonomer, but keeping the contacts of i with the
 residues of its own chain with rig[j]=1
 ********************************************************************/
void ResetContactsMonomerOutside(struct s_polymer *p, int i, int ci, const char *rig)
{
	struct s_back *bi = ((p+ci)->back)+i, *bn;
	int k,neigh,cn,nk=0;

	for (k=0;k<bi->ncontacts;++k)
	{
		neigh = bi->contacts[k];
		cn = bi->contacts_p[k];
		bn = ((p+cn)->back)+neigh;
		if (cn==ci && rig[neigh])						// kept, in the same order
		{
			if (nk != k)
			{
				bi->contacts[nk] = neigh;
				bi->contacts_p[nk] = cn;
				bi->contacts_r[nk] = bi->contacts_r[k];
				bi->e[nk] = bi->e[k];
				bn->contacts_r[bi->contacts_r[k]] = nk;			// i moved in its own list
				Touch(bn,cn,neigh);
			}
			nk++;
			continue;
		}
		RemoveContactEntry(p,bn,cn,neigh,bi->contacts_r[k]);		// remove i from the list of neigh
	}
	for (k=nk;k<bi->ncontacts;++k) bi->e[k] = 0;
	if (nk < bi->ncontacts) Touch(bi,ci,i);
	bi->ncontacts = nk;
}

/********************************************************************
 The residues whose backbone or sidechain atoms the rotamers of the
 sidechain of each residue are built from (AddSidechain): a move that
 changes the position of one of them must rebuild that sidechain too
 (see SideRange). Usually the residue itself and its neighbours; the
 last oxygen of a chain refers to the N two residues before.
 ********************************************************************/
void SetSidechainDependencies(struct s_polymer *p, struct s_mc_parms *parms, int nat)
{
	int ic,i,is,ir,k,ref[3],*owner,nwarn=0;
	struct s_back *b;
	struct s_rotamers *r;

	owner = (int *) malloc(nat*sizeof(int));
	parms->sidelo = (int **) malloc(parms->npol*sizeof(int *));
	parms->sidehi = (int **) malloc(parms->npol*sizeof(int *));
	if (!owner || !parms->sidelo || !parms->sidehi) Error("Cannot allocate in SetSidechainDependencies");
	parms->sidereach = 0;
	for (ic=0;ic<parms->npol;++ic)
	{
		parms->sidelo[ic] = (int *) malloc((p+ic)->nback*sizeof(int));
		parms->sidehi[ic] = (int *) malloc((p+ic)->nback*sizeof(int));
		if (!parms->sidelo[ic] || !parms->sidehi[ic]) Error("Cannot allocate in SetSidechainDependencies");
		for (k=0;k<nat;++k) owner[k] = -1;
		for (i=0;i<(p+ic)->nback;++i)
		{
			b = ((p+ic)->back)+i;
			if (b->ia>=0 && b->ia<nat) owner[b->ia] = i;
			for (is=0;is<b->nside;++is)
				if (b->side[is].ia>=0 && b->side[is].ia<nat) owner[b->side[is].ia] = i;
		}
		for (i=0;i<(p+ic)->nback;++i)
		{
			b = ((p+ic)->back)+i;
			parms->sidelo[ic][i] = parms->sidehi[ic][i] = i;
			for (is=0;is<b->nside;++is)
				for (ir=0;ir<b->nrot;++ir)
				{
					r = b->side[is].rot+ir;
					ref[0] = r->b1; ref[1] = r->b2; ref[2] = r->b3;
					for (k=0;k<3;++k)
					{
						if (ref[k]<0 || ref[k]>=nat || owner[ref[k]]<0)
							Error("A rotamer refers to an atom which is not in its chain");
						if (owner[ref[k]] != i && (((p+ic)->back)+owner[ref[k]])->ia != ref[k] && !nwarn++)
							fprintf(stderr,"WARNING: a rotamer of residue %d refers to a sidechain atom of residue %d: a change of rotamer\n"
								"         of the latter does not rebuild it\n",i,owner[ref[k]]);
						if (owner[ref[k]] < parms->sidelo[ic][i]) parms->sidelo[ic][i] = owner[ref[k]];
						if (owner[ref[k]] > parms->sidehi[ic][i]) parms->sidehi[ic][i] = owner[ref[k]];
					}
				}
			if (i-parms->sidelo[ic][i] > parms->sidereach) parms->sidereach = i-parms->sidelo[ic][i];
			if (parms->sidehi[ic][i]-i > parms->sidereach) parms->sidereach = parms->sidehi[ic][i]-i;
		}
	}
	free(owner);
}

/********************************************************************
 Extend the range [a,b] of residues of chain ip (of nback residues) that
 a move rebuilds and recomputes, to all the residues whose sidechain is
 built from residues in [m0,m1], those whose backbone the move changed
 ********************************************************************/
void SideRange(struct s_mc_parms *parms, int ip, int nback, int m0, int m1, int *a, int *b)
{
	int k,k0=m0-parms->sidereach,k1=m1+parms->sidereach;

	if (k0<0) k0 = 0;
	if (k1>nback-1) k1 = nback-1;
	for (k=k0;k<=k1;++k)
		if (parms->sidelo[ip][k] <= m1 && parms->sidehi[ip][k] >= m0)
		{
			if (k < *a) *a = k;
			if (k > *b) *b = k;
		}
}

/********************************************************************
 How far (in backbone units) the positions that the energy of a unit
 depends on can be from it: the reference atoms of the rotamers of its
 sidechain and, with hydrogen bonds, the atoms preceding its donors and
 acceptors. A pivot move rotates rigidly the units farther than this
 from the end of the rotated part, if their sidechains were consistent
 with their rotamers (see RigidFlags in mc.c). -1 if some reference is
 not found (pivot moves then compute all the pairs).
 ********************************************************************/
int RigidMargin(struct s_polymer *p, struct s_potential *u, int npol, int nat, int hb, int nosidechains)
{
	int ic,i,is,ir,k,d=0,ref[4],nref,*owner;
	struct s_back *b;
	struct s_rotamers *r;

	owner = (int *) malloc(nat*sizeof(int));
	if (!owner) Error("Cannot allocate in RigidMargin");
	for (ic=0;ic<npol;++ic)
	{
		for (k=0;k<nat;++k) owner[k] = -1;
		for (i=0;i<(p+ic)->nback;++i)
		{
			b = ((p+ic)->back)+i;
			if (b->ia>=0 && b->ia<nat) owner[b->ia] = i;
			for (is=0;is<b->nside;++is)
				if (b->side[is].ia>=0 && b->side[is].ia<nat) owner[b->side[is].ia] = i;
		}
		for (i=0;i<(p+ic)->nback;++i)
		{
			b = ((p+ic)->back)+i;
			for (is=0;is<b->nside && !nosidechains;++is)
				for (ir=0;ir<b->nrot;++ir)
				{
					r = b->side[is].rot+ir;
					ref[0] = r->b1; ref[1] = r->b2; ref[2] = r->b3; nref = 3;
					if (hb && b->side[is].ia>=0 && b->side[is].ia<nat && u->hb[b->side[is].ia]) ref[nref++] = u->hb_iam[b->side[is].ia];
					for (k=0;k<nref;++k)
					{
						if (ref[k]<0 || ref[k]>=nat || owner[ref[k]]<0) { free(owner); return -1; }
						if (Abs(owner[ref[k]]-i) > d) d = Abs(owner[ref[k]]-i);
					}
				}
			if (hb && b->ia>=0 && b->ia<nat && u->hb[b->ia])
			{
				k = u->hb_iam[b->ia];
				if (k<0 || k>=nat || owner[k]<0) { free(owner); return -1; }
				if (Abs(owner[k]-i) > d) d = Abs(owner[k]-i);
			}
		}
	}
	free(owner);
	return d;
}

/********************************************************************
 Read the total interaction energy of a monomer with the others
 from its energy structure
 ********************************************************************/
double GetEnergyMonomer(struct s_polymer *p, int ip, int iw)
{
	double e=0;
	int i;

	for (i=0;i<(((p+ip)->back)+iw)->ncontacts;++i)
	{
		e += ((((p+ip)->back)+iw)->e)[i];
	}
	return e;
}

/********************************************************************
 Read the total interaction energy of a set of monomers
 from its energy structure (ip is the id of the chain)
 ********************************************************************/
double GetEnergyMonomerRange(struct s_polymer *p, int from, int to, int ip)
{
	return GetEnergyMonomerRangeRigid(p,from,to,ip,NULL);
}

/********************************************************************
 The same without the pairs of residues of chain ip with rig[i]=1
 ********************************************************************/
double GetEnergyMonomerRangeRigid(struct s_polymer *p, int from, int to, int ip, const char *rig)
{
	double e=0;
	int i,iw,j,cj;
	
	for (iw=from;iw<=to;++iw)
	{
		if (iw>=0 && iw<(p+ip)->nback)
		{
			for (i=0;i<(((p+ip)->back)+iw)->ncontacts;++i)
			{
				j = ((((p+ip)->back)+iw)->contacts)[i];
				cj = ((((p+ip)->back)+iw)->contacts_p)[i];
				if ( (ip!=cj || j<from || j>to || iw<j)				// count only once contacts in [from,to] of the same chain
					&& !(rig && ip==cj && rig[iw] && rig[j]) )
					e += ((((p+ip)->back)+iw)->e)[i];
			}
		}
	}
	return e;
}

/********************************************************************
 Read the H Fields interaction energy
 from its energy structure
 ********************************************************************/
double GetHFields(struct s_polymer *p, int npol)
{
	double e=0.;
	int ic,iw;

	// a move changes the number of contacts also of residues of other chains
	for (ic=0;ic<npol;++ic)
		for (iw=0;iw<(p+ic)->nback;++iw)
			e += (((p+ic)->back)+iw)->e_hfs;
	return e;
}

/********************************************************************
 Calculate (and store) the H fields energy of all the chains
 ********************************************************************/
double EnergyHFieldsAll(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms)
{
	double e=0.;
	int ic,iw;

	for (ic=0;ic<parms->npol;++ic)
		for (iw=0;iw<(p+ic)->nback;++iw)
			e += EnergyHFields(p,u,parms,iw,ic,1);
	return e;
}

/********************************************************************
 Print contacts and partial energies for debugging purposes
 ********************************************************************/
void PrintContacts(FILE *fp, struct s_polymer *p,int ip, unsigned long long step)
{
	int i,j;

	fprintf(fp,"\n\n------------------------\nSTEP = %llu\n",step);

	for (i=0;i<(p+ip)->nback;++i)
	{
		fprintf(fp,"Back = %d\n\tside =",i);
		for (j=0;j<(((p+ip)->back)+i)->nside;++j)
			fprintf(fp,"%d ",(((((p+ip)->back)+i)->side)+j)->ia);
		fprintf(fp,"\nContacts = ");
		for (j=0;j<(((p+ip)->back)+i)->ncontacts;++j)
			fprintf(fp,"%d ",((((p+ip)->back)+i)->contacts)[j]);
		fprintf(fp,"\nContactsP= ");
	        for (j=0;j<(((p+ip)->back)+i)->ncontacts;++j)
                fprintf(fp,"%d ",((((p+ip)->back)+i)->contacts_p)[j]);
		fprintf(fp,"\n");
		for (j=0;j<(((p+ip)->back)+i)->ncontacts;++j)
			fprintf(fp,"%lf ",((((p+ip)->back)+i)->e)[j]);
		fprintf(fp,"\n");
	}
}


void CountContacts(FILE *fp,struct s_polymer *polymer,struct s_mc_parms *parms,unsigned long long step)
{
	int i,tempncontacts,iback,ipol,icont,matcont[parms->npol];
    for(i=0;i<parms->npol;++i)
	{
		for(ipol=0;ipol<parms->npol;++ipol)
			matcont[ipol]=0;
		tempncontacts=0;
		fprintf(fp,"> Chain %d\n  %d backbone atoms\n",i,(polymer+i)->nback);
		for(iback=0;iback<(polymer+i)->nback;++iback)
		{
			tempncontacts+=((polymer+i)->back+iback)->ncontacts;
			for(icont=0;icont<((polymer+i)->back+iback)->ncontacts;++icont)
				matcont[ ((polymer+i)->back+iback)->contacts_p[icont] ]++;
                            
                                
		}
     
		fprintf(fp,"  %d contacts ( ",tempncontacts);
		for(ipol=0;ipol<parms->npol;++ipol)
			fprintf(fp," %d ",matcont[ipol]);
		fprintf(fp,")\n");
	}
}



void SetTrueShellFlags(struct s_polymer *p,struct s_mc_parms *parms){
 int ip,i;
 #ifdef _OPENMP
 #pragma omp parallel for
 #endif
 for(ip=0;ip<parms->npol;++ip)
   for(i=0;i<(p+ip)->nback;++i)
      (((p+ip)->back)+i)->shellflag=true;

}

void SetFalseShellFlags(struct s_polymer *p,struct s_mc_parms *parms){
 int ip,i;
 #ifdef _OPENMP
 #pragma omp parallel for
 #endif
 for(ip=0;ip<parms->npol;++ip)
   for(i=0;i<(p+ip)->nback;++i)
      (((p+ip)->back)+i)->shellflag=false;

}


void UpdateShellFlags(struct s_polymer *p,int ip,int ifrom,int ito){

  int i;
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(i=ifrom;i<ito;i++)
	(((p+ip)->back)+i)->shellflag=true;


}


/********************************************************************
 Update shells of neighbours
 ********************************************************************/
void UpdateShell(struct s_polymer *p, struct s_mc_parms *parms)
{
 
	int i,j,is,js,ci,cj;
	double r2;
        bool tmpflag;
	//for (ci=0;ci<parms->npol;++ci)
         // for (i=0;i<(p+ci)->nback;++i)
          //  (((p+ci)->back)+i)->nshell = 0;


	// update contacts shells
	for (ci=0;ci<parms->npol;++ci)
	{
		for (i=0;i<(p+ci)->nback;++i)
		{
			if( (((p+ci)->back)+i)->shellflag==true)
			{ // update shell of i and its neighbours
				(((p+ci)->back)+i)->nshell=0;
				for(cj=0;cj<parms->npol;++cj)
				{
					for(j=0;j<(p+cj)->nback;++j)
					{
						tmpflag=false;
						if( (((p+cj)->back)+j)->shellflag==false ) // if j is not a moved residue
						{ // interaction backbone i - backbone j
							r2=Dist2( (((p+ci)->back)+i)->pos, (((p+cj)->back)+j)->pos );
							if( r2<parms->r2shell) 
							{
							  AddShell(p,i,j,ci,cj);
						          tmpflag=true;
							}
							else
							{
								// interaction backbone i - sidechain j
		                                                for (js=0;js< (((p+cj)->back)+j)->nside;++js)
								{ // interaction backbone i - sidechain j
									if(tmpflag==true) break;
		                                                        r2 = Dist2( (((p+ci)->back)+i)->pos, (((((p+cj)->back)+j)->side)+js)->pos );
									if ( r2 < parms->r2shell )
									{
										AddShell(p,i,j,ci,cj);
										tmpflag=true;
										break; // if shell is added, then check the next backbone i+1 of polymer ci
									}
								}
								if(tmpflag==false) 
								  //interaction backbone j - sidechain i
								  for (is=0;is< (((p+ci)->back)+i)->nside;++is) 
                                                       		  {			
									
                		              	                         	r2 = Dist2( (((((p+ci)->back)+i)->side)+is)->pos, ((((p+cj)->back)+j)->pos ));
		                                                        	if ( r2 < parms->r2shell )
										{
											AddShell(p,i,j,ci,cj);                       // add to shell
											tmpflag=true;
											break;
										}
                                                       		} // end of for loop backbone j - sidechain i
							} // end of else
														
							if(tmpflag==false)	
							  // interaction sidechain i - sidechain j
							  for (is=0;is< (((p+ci)->back)+i)->nside;++is)
	                                                  {
								if(tmpflag==true) break;
								for(js=0;js< (((p+cj)->back)+j)->nside;++js)
								{
                                                        	  r2 = Dist2( (((((p+cj)->back)+j)->side)+js)->pos, (((((p+ci)->back)+i)->side)+is)->pos );
                                                                  if ( r2 < parms->r2shell )
								  {
									AddShell(p,i,j,ci,cj);                       // add to shell
									tmpflag=true;
									break;
                                                		  } 
								} 
							  } // end of loop on sidechain atoms of i
						} // end of if j is not moved
		
				
					

					} // end of loop on backbone atoms of polymer cj
	
				} // end of loop on polymers cj
			
				// now set the flag of backbone i in polymer ci to zero
                               (((p+ci)->back)+i)->shellflag=false;
                   	      
			} //end of if backbone i of ci is moved
		        	
                        


		} // end of loop on backbone atoms i of ci

	} // end of loop on polymers ci
	


	
/* old code starts here
	 for (ci=0;ci<parms->npol;++ci)
                for (i=0;i<(p+ci)->nback;++i)
			(((p+ci)->back)+i)->nshell = 0;


	//  double loop
	#ifdef _OPENMP
	#pragma omp parallel for private(cj,i,j,r2,js,is)
	#endif
	for (ci=0;ci<parms->npol;++ci)
	{
		for (cj=ci;cj<parms->npol;++cj)
		{
			for (i=0;i<(p+ci)->nback;++i)
			{
				for (j=0;j<(p+cj)->nback;++j)
				{
					if (ci!=cj || i<j)
					{
						// interaction backbone-backbone (i-j)
						r2 = Dist2( (((p+ci)->back)+i)->pos, (((p+cj)->back)+j)->pos );
						if ( r2 < parms->r2shell ) AddShell(p,i,j,ci,cj);					// add to shell

						for (js=0;js< (((p+cj)->back)+j)->nside;++js)
						{
							// interaction backbone-sidechain (i-js)
							r2 = Dist2( (((p+ci)->back)+i)->pos, (((((p+cj)->back)+j)->side)+js)->pos );
							if ( r2 < parms->r2shell ) AddShell(p,i,j,ci,cj);			// add to shell

							for (is=0;is< (((p+ci)->back)+i)->nside;++is)
							{
								// interaction sidechain-sidechain (is-js)
								r2 = Dist2( (((((p+ci)->back)+i)->side)+is)->pos, (((((p+cj)->back)+j)->side)+js)->pos );
//								if(i==4&&j==148) fprintf(stderr,"updateshell: i=4,j=148;\tdist2=%lf\n",r2); 
								if ( r2 < parms->r2shell ) AddShell(p,i,j,ci,cj);			// add to shell
							}

						}

						// interaction backbone-sidechain (is-j)
						for (is=0;is< (((p+ci)->back)+i)->nside;++is)
						{
								r2 = Dist2( (((p+cj)->back)+j)->pos, (((((p+ci)->back)+i)->side)+is)->pos );
								if ( r2 < parms->r2shell ) AddShell(p,i,j,ci,cj);			// add to shell
						}
					}
				}
			}
		}
	}
*/ // old code ends here
}



void CopyShell(struct s_polymer *from,struct s_polymer *to,struct s_mc_parms *parms)
{

	int i,ci,cs;
	for (ci=0;ci<parms->npol;++ci)
	{
		for (i=0;i<(from+ci)->nback;++i)
		{
                        (((to+ci)->back)+i)->nshell = (((from+ci)->back)+i)->nshell ;
			for(cs=0;cs<(((from+ci)->back)+i)->nshell;++cs)
			{
				(((to+ci)->back)+i)->shell[cs]=(((from+ci)->back)+i)->shell[cs];		
			        (((to+ci)->back)+i)->shell_p[cs]=(((from+ci)->back)+i)->shell_p[cs];	
			}
		}
	}
}



/********************************************************************
 Angular potential of backbone iw of chain ic
 ********************************************************************/
double EnergyAngles(struct s_polymer *p, struct s_potential *u, int iw, int ic, int update)
{
	int ia,out;
	double ang,e;

	if (iw<1 || iw>(p+ic)->nback-2) return 0;			// angles are not defined for the first and last backbone atom

	ia = (((p+ic)->back)+iw)->ia;

	ang = Angle( (((p+ic)->back)+iw-1)->pos, (((p+ic)->back)+iw)->pos, (((p+ic)->back)+iw+1)->pos, p->tables, &out );
	if (out==0) return LARGE;

	e = (u->e_ang)[ia] * ( ang - (u->ang0)[ia] ) * ( ang - (u->ang0)[ia] );

	if (update)  (((p+ic)->back)+iw)->e_ang = e;

	return e;
}

/********************************************************************
 Dihedral potential of backbone iw of chain ic
 (defined as dihedral among iw-2,iw-1,iw,iw+1)
 ********************************************************************/
double EnergyDihedrals(struct s_polymer *p, struct s_potential *u, int iw, int ic, int update)
{
	int ia,iaa,i, out;
	double e=0,dih01,dih03,dih;

	if (iw<2 || iw>(p+ic)->nback-2) return 0;			// dihedrals are not defined for the first and last two backbone atoms

	ia = (((p+ic)->back)+iw)->ia;
        iaa = (((p+ic)->back)+iw)->iaa - (((p+ic)->back))->iaa;

	dih01 = (u->dih01)[ia];
	dih03 = (u->dih03)[ia];

	dih = Dihedral( (((p+ic)->back)+iw-2)->pos, (((p+ic)->back)+iw-1)->pos, (((p+ic)->back)+iw)->pos, (((p+ic)->back)+iw+1)->pos, p->tables, &out );
	if (out==0) { fprintf(stderr,"WARNING: Error in evaluating dihedral energy"); return LARGE; }

	// periodic kind of potential, i.e. e1*cos(d-d0) + e3*cos(3(d-d0))
	if (u->dih_periodic)
		e = (u->e_dih1)[ia] * ( 1 - FastCos(dih-dih01,p->tables) ) + (u->e_dih3)[ia] * ( 1 - FastCos(3.*(dih-dih03),p->tables) );

	// tabled kind of potential
	if (u->dih_tabled)
	{
		i = (int) ((dih+180.)/u->g_dihbin);
        #ifdef DEBUG
	      if (i<0 || i>NDIHFMAX) { fprintf(stderr,"WARNING: dih=%lf i=%d in EnergyDihedrals\n",dih,i); return LARGE; }
        #endif
		if ((u->dih_which)[ia] == 0)  // is a phi
			e += u->g_dihe * (u->dih_pa)[ia] * (u->dih_f_phi_a)[i] + u->g_dihe * (u->dih_pb)[ia] * (u->dih_f_phi_b)[i];
		else						  // is a psi
			e += u->g_dihe * (u->dih_pa)[ia] * (u->dih_f_psi_a)[i] + u->g_dihe * (u->dih_pb)[ia] * (u->dih_f_psi_b)[i];
	}

    // Gaussian potential on dihedrals (Propensity)
	
    if(u->dih_ram)
    {
        if((iw+1)%3==0) // psi dihedral
		{
            for(i=0;i<2;++i)
            { 
                e += - (u->e_dihram) * (u->ab_propensity[i][iaa]/u->sigma[i][1]) * exp(-0.5*(((dih - u->dih0[i][1])/u->sigma[i][1])*((dih - u->dih0[i][1])/u->sigma[i][1])));
            }
		}
        if((iw-1)%3==0)   // phi dihedral
		{
            for(i=0;i<2;++i)
             {
                e += - (u->e_dihram) * ( (u->ab_propensity[i][iaa]/u->sigma[i][0]) * exp(-0.5*(((dih - u->dih0[i][0])/u->sigma[i][0])*((dih - u->dih0[i][0])/u->sigma[i][0]))) );
             }
		}
    }
    /*
    //buca infinita sui diedri
    if(u->dih_ram)
    {
        if((iw+1)%3==0) {   //di tipo psi
            for(i=0;i<2;++i) {
                if ( u->ab_propensity[i][iaa] != 0) {
                    //fprintf(stderr,"iaa: %d\t propensity: %lf\tdihpsi: %d\tsigmapsi: %lf\tdih %lf  \n",iaa+1,u->ab_propensity[i][iaa],u->dih0[i][1],u->sigma[i][1],dih);
                    if ( dih >= (u->dih0[i][1]-u->sigma[i][1]) && dih <= (u->dih0[i][1]+u->sigma[i][1])) e+=0;
                    else e+= LARGE;
                }
            }
        }
        if((iw-1)%3==0) {   //di tipo phi
            for(i=0;i<2;++i) {
                if ( u->ab_propensity[i][iaa] != 0) {
                    //fprintf(stderr,"iaa: %d\t propensity: %lf\tdihphi: %d\tsigmaphi: %lf\tdih %lf  \n",iaa+1,u->ab_propensity[i][iaa],u->dih0[i][0],u->sigma[i][0],dih);
                    if ( dih >= (u->dih0[i][0]-u->sigma[i][0]) && dih <= (u->dih0[i][0]+u->sigma[i][0])) e+=0;
                    else e+= LARGE;
                }
            }
        }
        
    }*/



	if (update)  (((p+ic)->back)+iw)->e_dih = e;

	return e;
}

/***********************************************
 H fields potential
 ***********************************************/

double EnergyHFields(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int iw, int ic, int update){

	double e=0.;
	
	// Loop over all the sidechain atoms
	e += (double)(((p+ic)->back)+iw)->ncontacts * u->h_values[((((p+ic)->back)+iw)->side)->itype];
	
	// Store the energy value on the backbone atom
	if (update)
	{
		(((p+ic)->back)+iw)->e_hfs = e;
	}
		
	return e;
}

/********************************************************************
 Print energy terms
 ********************************************************************/
void PrintEnergies(FILE *fp, int nc, unsigned long long step, struct s_polymer *p)
{
	int i,j;
	double epair=0,eang=0,edih=0,ehfs=0;
	
	for (i=0;i<nc;++i)
		for (j=0;j<(p+i)->nback;++j)
		{
			epair += GetEnergyMonomer(p,i,j) / 2.;
			eang += (((p+i)->back)+j)->e_ang;
			edih += (((p+i)->back)+j)->e_dih;
			ehfs += (((p+i)->back)+j)->e_hfs;
		}

	fprintf(fp,"%llu\t%9lf\t%9lf\t%9lf\t%9lf\t%9lf\n",step,epair+eang+edih+ehfs,epair,eang,edih,ehfs);
}


/********************************************************************
 Print energy terms (MPI routine)
 ********************************************************************/
void PrintEnergies_Parallel(FILE *fp, int nc, unsigned long long step, struct s_polymer *p,int my_rank)
{
	int i,j;
	double epair=0,eang=0,edih=0,ehfs=0;

	for (i=0;i<nc;++i)
	{
		for (j=0;j<(p+i)->nback;++j)
		{
			epair += GetEnergyMonomer(p,i,j)/2.;
			eang += (((p+i)->back)+j)->e_ang;
			edih += (((p+i)->back)+j)->e_dih;
			ehfs += (((p+i)->back)+j)->e_hfs;
		}
	}
	
	fprintf(fp,"%d\t%llu\t%9lf\t%9lf\t%9lf\t%9lf\t%9lf\n",my_rank,step,epair+eang+edih+ehfs,epair,eang,edih,ehfs);
}

/********************************************************************
 Potential copy
 ********************************************************************/
void CopyPotential(struct s_potential *from, struct s_potential *to, int nat, int ntypes)
{
	int i,j;

	for (i=0;i<ntypes;++i)
		for (j=0;j<ntypes;++j)
		{
			to->r0_2[i][j] = from->r0_2[i][j];
			to->r_2[i][j] = from->r_2[i][j];
			to->e[i][j] = from->e[i][j];
		}

	for (i=0;i<nat;++i)
	{
		to->ang0[i] = from->ang0[i];
		to->e_ang[i] = from->e_ang[i];
		to->e_dih1[i] = from->e_dih1[i];
		to->dih01[i] = from->dih01[i];
		to->e_dih3[i] = from->e_dih3[i];
		to->dih03[i] = from->dih03[i];
	}

	to->g_r0hard = from->g_r0hard;
	to->g_ehomo = from->g_ehomo;
	to->g_rhomo = from->g_rhomo;
	to->g_anglek = from->g_anglek;
	to->g_angle0 = from->g_angle0;
	to->g_dihk1 = from->g_dihk1;
	to->g_dih01 = from->g_dih01;
	to->g_dihk3 = from->g_dihk3;
	to->g_dih03 = from->g_dih03;
}

/********************************************************************
 Check overlaps
 	 pr=1	print overlaps
 	 pr=0   only return 1 if something overlaps
 ********************************************************************/

int CheckOverlaps(struct s_polymer *p, struct s_potential *u, struct s_mc_parms *parms, int npol, int nosidechains, int pr, FILE *fproc)
{
	int o=0,tooclose;
	int i,j,ci,cj,a1,a2,is,js,ia1,ia2,iaa1,iaa2;
	char atom1[5],am1[5];
	char atom2[5],am2[5];
	double r2,r02;
	struct vector sh={0.,0.,0.};		// periodic box: image of j nearest to i

	if (pr) fprintf(fproc,"\nChecking overlaps... r0hard = %lf\n",u->g_r0hard);
	// backbone loop
	for (ci=0;ci<npol;++ci)                                     // loops on polymers
	{
		for (cj=ci;cj<npol;++cj)
		{
			for (i=0;i<(p+ci)->nback;++i)             // loops on backbone atoms of each polymer
			{
				for (j=0;j<(p+cj)->nback;++j)
				{
					if (ci!=cj || (i<j && (((p+ci)->back)+i)->iaa != (((p+cj)->back)+j)->iaa))                      // in the same polymer, do not count twice a contact
					{
						if (u->boxtype=='p') sh = ImageShift( (((p+ci)->back)+i)->pos, (((p+cj)->back)+j)->pos, u->lbox );
						tooclose = 0;
						if (ci==cj && j-i <= u->g_imin ) tooclose = 1;

						// interaction backbone-backbone (i-j)
						if (ci!=cj || (i!=j+1 && i!=j-1)) // if not consecutive (i.e., no covalent interaction)
						{
							a1 = (((p+ci)->back)+i)->itype;
							a2 = (((p+cj)->back)+j)->itype;
							r2 = Dist2S( (((p+ci)->back)+i)->pos, (((p+cj)->back)+j)->pos, sh );
							if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
							else r02 = (u->r0_2)[a1][a2];
							if ( r2 < r02 )
							{
								strcpy(atom1, (((p+ci)->back)+i)->type);
								ia1= (((p+ci)->back)+i)->ia;
								strcpy(am1, (((p+ci)->back)+i)->aa);
								iaa1= (((p+ci)->back)+i)->iaa;
								strcpy(atom2,(((p+cj)->back)+j)->type);
								ia2= (((p+cj)->back)+j)->ia;
								strcpy(am2,(((p+cj)->back)+j)->aa);
								iaa2= (((p+cj)->back)+j)->iaa;

								o=1;
								if (pr) fprintf(fproc,"WARNING: overlap between backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d) and backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d): d=%lf d0=%lf\n",
									i,ci,a1,atom1,ia1,am1,iaa1,j,cj,a2,atom2,ia2,am2,iaa2,sqrt(r2),sqrt((u->r0_2)[a1][a2]));
							}
						}

						for (js=0;js< (((p+cj)->back)+j)->nside;++js)
						{
							// interaction backbone-sidechain (i-js)
							a1 = (((p+ci)->back)+i)->itype;
							a2 = (((((p+cj)->back)+j)->side)+js)->itype;
							r2 = Dist2S( (((p+ci)->back)+i)->pos, (((((p+cj)->back)+j)->side)+js)->pos, sh );
							if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
							else r02 = (u->r0_2)[a1][a2];
							if ( r2 < r02 )
							{
								strcpy(atom1, (((p+ci)->back)+i)->type);
								ia1= (((p+ci)->back)+i)->ia;
								strcpy(am1, (((p+ci)->back)+i)->aa);
								iaa1= (((p+ci)->back)+i)->iaa;
								strcpy(atom2, (((((p+cj)->back)+j)->side)+js)->type);
								ia2= (((((p+cj)->back)+j)->side)+js)->ia;
								strcpy(am2, (((p+cj)->back)+j)->aa);
								iaa2= (((p+cj)->back)+j)->iaa;

								o=1;
								if (pr) fprintf(fproc,"WARNING: overlap between backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d) and sidechain %d of backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d): d=%lf d0=%lf\n",
                                                                                          i,ci,a1,atom1,ia1,am1,iaa1,js,j,cj,a2,atom2,ia2,am2,iaa2,sqrt(r2),sqrt((u->r0_2)[a1][a2]));
							}

							for (is=0;is< (((p+ci)->back)+i)->nside;++is)
							{
								// interaction sidechain-sidechain (is-js)
								a1 = (((((p+ci)->back)+i)->side)+is)->itype;
								a2 = (((((p+cj)->back)+j)->side)+js)->itype;
								r2 = Dist2S( (((((p+ci)->back)+i)->side)+is)->pos, (((((p+cj)->back)+j)->side)+js)->pos, sh );
								if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
								else r02 = (u->r0_2)[a1][a2];
								if ( r2 < r02 )
								{
									strcpy(atom1, (((((p+ci)->back)+i)->side)+is)->type);
									ia1= (((((p+ci)->back)+i)->side)+is)->ia;
									strcpy(am1, (((p+ci)->back)+i)->aa);
									iaa1= (((p+ci)->back)+i)->iaa;
									strcpy(atom2, (((((p+cj)->back)+j)->side)+js)->type);
									ia2= (((((p+cj)->back)+j)->side)+js)->ia;
									strcpy(am2, (((p+cj)->back)+j)->aa);
									iaa2= (((p+cj)->back)+j)->iaa;

									o=1;
									if (pr)
									{
										fprintf(fproc,"WARNING: overlap between sidechain %d of backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d) and sidechain %d of backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d): d=%lf d0=%lf\n",is,i,ci,a1,atom1,ia1,am1,iaa1,js,j,cj,a2,atom2,ia2,am2,iaa2,sqrt(r2),sqrt((u->r0_2)[a1][a2]));
									}
								}
							}

						}

                                    // interaction backbone-sidechain (is-j)
						for (is=0;is< (((p+ci)->back)+i)->nside;++is)
						{
							a1 = (((p+cj)->back)+j)->itype;
							a2 = (((((p+ci)->back)+i)->side)+is)->itype;
							r2 = Dist2S( (((((p+ci)->back)+i)->side)+is)->pos, (((p+cj)->back)+j)->pos, sh );
							if (tooclose) r02 = u->g_r0hard * u->g_r0hard;
							else r02 = (u->r0_2)[a1][a2];
							if ( r2 < r02 )
							{
								strcpy(atom1, (((p+cj)->back)+j)->type);
								ia1= (((p+cj)->back)+j)->ia;
								strcpy(am1, (((p+cj)->back)+j)->aa);
								iaa1= (((p+cj)->back)+j)->iaa;
								strcpy(atom2, (((((p+ci)->back)+i)->side)+is)->type);
								ia2= (((((p+ci)->back)+i)->side)+is)->ia;
								strcpy(am2, (((p+ci)->back)+i)->aa);
								iaa2= (((p+ci)->back)+i)->iaa;

								o=1;
								if (pr) fprintf(fproc,"WARNING: overlap between sidechain %d of backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d) and backbone %d (chain %d, type %d, atom %s, id %d of amino acid %s %d): d=%lf d0=%lf\n",
									is,i,ci,a1,atom1,ia1,am1,iaa1,j,cj,a2,atom2,ia2,am2,iaa2,sqrt(r2),sqrt((u->r0_2)[a1][a2]));
							}
						}

					}
				}
			}
		}
	}
	return o;
}



/********************************************************************
 Energy of chains in a box, returns 1 if any backbone atom is outside
 the box, 0 iw all right.
 'c'=cubic box, 's'=spherical box
 ********************************************************************/
int EnergyBox(struct s_polymer *p, struct s_potential *u, int iw, int ic)
{
	if (u->boxtype == 'p') return 0;				// periodic box: no walls

	double d2;

	// cubic box
	if (u->boxtype == 'c')
	{
		if ( ((((p+ic)->back)+iw)->pos).x < -u->boxsize || ((((p+ic)->back)+iw)->pos).x > u->boxsize ) return 1;
		if ( ((((p+ic)->back)+iw)->pos).y < -u->boxsize || ((((p+ic)->back)+iw)->pos).y > u->boxsize ) return 1;
		if ( ((((p+ic)->back)+iw)->pos).z < -u->boxsize || ((((p+ic)->back)+iw)->pos).z > u->boxsize ) return 1;
	}

	// spherical box
	else if (u->boxtype == 's')
	{
		d2 = (((((p+ic)->back)+iw)->pos).x * ((((p+ic)->back)+iw)->pos).x) + (((((p+ic)->back)+iw)->pos).y * ((((p+ic)->back)+iw)->pos).y) +
			(((((p+ic)->back)+iw)->pos).z * ((((p+ic)->back)+iw)->pos).z);
		if (d2 > u->boxsize * u->boxsize) return 1;
	}

	return 0;
}

int EnergyBoxPolymer(struct s_polymer *p, struct s_potential *u, int ic)
{
	int i;

	for (i=0;i<(p+ic)->nback;++i)
		if (EnergyBox(p,u,i,ic)==1) return 1;

	return 0;

}
