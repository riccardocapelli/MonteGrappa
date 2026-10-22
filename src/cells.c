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
 * cells.c
 *
 * Cell list of the backbone atoms of the system, used by EnergyMonomerRange
 * to find the residues which can interact with a moved one without scanning
 * the whole system.
 *
 * The space is divided in cubic cells of side L >= rcut + 2*extcap, where
 * rcut is the interaction range of the potential and extcap bounds the extent
 * of the sidechains (s_back.ext), so that two residues in non-neighbouring
 * cells cannot interact. Cells are hashed into a table of buckets (the system
 * has no box); a bucket can contain residues of far cells, which are then
 * discarded by the distance test of the caller. Each bucket is a doubly
 * linked list of residues, indexed globally as offset[chain]+residue.
 *
 * The list follows the coordinates lazily: EnergyMonomerRange updates the
 * residues it is going to evaluate, UpdateChainRanges the residues restored
 * after a rejected move and TotalEnergy rebuilds it after bulk changes. A
 * sidechain whose extent exceeds extcap marks the list for a rebuild.
 *
 * Candidates are visited in the same order as the plain loop over all
 * residues, so the energies (and trajectories) are identical with and
 * without cell list; at each rebuild the list is switched off when the
 * system spans too few cells for it to pay off (unless forced).
 */

#include "montegrappa.h"

#define CELLS_EXT_MARGIN 0.5		// slack added to the largest sidechain extent (A)

/*****************************************************************************
 Allocate the cell list of system p (attached to p->cells), whose npol chains
 have the same number of residues as those of topo
 *****************************************************************************/
void CellsAlloc(struct s_polymer *p, struct s_polymer *topo, int npol, int forced, int allowed)
{
	struct s_cells *c;
	int ic,i,g;

	c = (struct s_cells *) calloc(1,sizeof(struct s_cells));
	if (!c) Error("Cannot allocate cell list");

	c->npol = npol;
	c->offset = (int *) calloc(npol+1,sizeof(int));
	if (!c->offset) Error("Cannot allocate cell list");
	for (ic=0;ic<npol;++ic) c->offset[ic+1] = c->offset[ic] + (topo+ic)->nback;
	c->nres = c->offset[npol];

	c->hsize = 1024;
	while (c->hsize < 2*c->nres) c->hsize *= 2;
	c->mask = c->hsize - 1;

	c->head = (int *) malloc(c->hsize*sizeof(int));
	c->next = (int *) malloc(c->nres*sizeof(int));
	c->prev = (int *) malloc(c->nres*sizeof(int));
	c->bucket = (int *) malloc(c->nres*sizeof(int));
	c->chain = (int *) malloc(c->nres*sizeof(int));
	c->res = (int *) malloc(c->nres*sizeof(int));
	c->cand = (int *) malloc(c->nres*sizeof(int));
	c->x = (double *) malloc(c->nres*sizeof(double));
	c->y = (double *) malloc(c->nres*sizeof(double));
	c->z = (double *) malloc(c->nres*sizeof(double));
	c->ext = (double *) malloc(c->nres*sizeof(double));
	if (!c->head || !c->next || !c->prev || !c->bucket || !c->chain || !c->res || !c->cand) Error("Cannot allocate cell list");
	if (!c->x || !c->y || !c->z || !c->ext) Error("Cannot allocate cell list");

	for (ic=0;ic<npol;++ic)
		for (i=0;i<(topo+ic)->nback;++i)
		{
			g = c->offset[ic]+i;
			c->chain[g] = ic;
			c->res[g] = i;
		}

	c->dirty = 1;
	c->forced = forced;
	c->allowed = allowed || forced;
	c->active = 0;
	p->cells = c;
}

void CellsFree(struct s_polymer *p)
{
	struct s_cells *c = p->cells;

	if (!c) return;
	free(c->offset); free(c->head); free(c->next); free(c->prev);
	free(c->bucket); free(c->chain); free(c->res); free(c->cand);
	free(c->x); free(c->y); free(c->z); free(c->ext);
	free(c);
	p->cells = NULL;
}

static inline int CellCoord(double x, double invL)
{
	return (int) floor(x*invL);
}

static inline int CellBucket(struct s_cells *c, int ix, int iy, int iz)
{
	return (int) ( ( (unsigned) ix * 73856093u ^ (unsigned) iy * 19349663u ^ (unsigned) iz * 83492791u ) & (unsigned) c->mask );
}

static inline void CellsUnlink(struct s_cells *c, int g)
{
	if (c->prev[g] >= 0) c->next[c->prev[g]] = c->next[g];
	else c->head[c->bucket[g]] = c->next[g];
	if (c->next[g] >= 0) c->prev[c->next[g]] = c->prev[g];
}

static inline void CellsLink(struct s_cells *c, int g, int h)
{
	c->bucket[g] = h;
	c->prev[g] = -1;
	c->next[g] = c->head[h];
	if (c->head[h] >= 0) c->prev[c->head[h]] = g;
	c->head[h] = g;
}

/*****************************************************************************
 Rebuild the whole list from the current coordinates
 *****************************************************************************/
void CellsBuild(struct s_polymer *p, double rcut)
{
	struct s_cells *c = p->cells;
	struct s_back *b;
	int g,h,d;
	double extmax=0,lo[3]={1e30,1e30,1e30},hi[3]={-1e30,-1e30,-1e30},v[3],ncells;

	for (g=0;g<c->nres;++g)
	{
		b = ((p+c->chain[g])->back)+c->res[g];
		c->x[g] = b->pos.x; c->y[g] = b->pos.y; c->z[g] = b->pos.z; c->ext[g] = b->ext;
		if (b->ext > extmax) extmax = b->ext;
		v[0] = b->pos.x; v[1] = b->pos.y; v[2] = b->pos.z;
		for (d=0;d<3;++d)
		{
			if (v[d] < lo[d]) lo[d] = v[d];
			if (v[d] > hi[d]) hi[d] = v[d];
		}
	}
	c->rcut = rcut;
	c->extcap = extmax + CELLS_EXT_MARGIN;
	c->L = rcut + 2.*c->extcap;
	c->invL = 1./c->L;

	// worth it if the expected candidates (27 cells out of those spanned by the system)
	// plus the cost of visiting the cells are well below scanning all residues
	ncells = 1.;
	for (d=0;d<3;++d) ncells *= (hi[d]-lo[d])*c->invL + 1.;
	c->active = c->forced || ( c->allowed && c->nres*27./ncells + 60. < 0.5*c->nres );
	c->dirty = 0;
	if (!c->active) return;

	for (h=0;h<c->hsize;++h) c->head[h] = -1;
	for (g=0;g<c->nres;++g)
	{
		b = ((p+c->chain[g])->back)+c->res[g];
		h = CellBucket(c,CellCoord(b->pos.x,c->invL),CellCoord(b->pos.y,c->invL),CellCoord(b->pos.z,c->invL));
		CellsLink(c,g,h);
	}
}

/*****************************************************************************
 Move residues [from,to] of chain ic to the cells of their current position
 *****************************************************************************/
void CellsUpdateRange(struct s_polymer *p, int ic, int from, int to)
{
	struct s_cells *c = p->cells;
	struct s_back *b;
	int i,g,h;

	if (c->dirty) return;					// will be rebuilt before use
	if (from<0) from=0;
	if (to>(p+ic)->nback-1) to=(p+ic)->nback-1;

	for (i=from;i<=to;++i)
	{
		b = ((p+ic)->back)+i;
		g = c->offset[ic]+i;
		c->x[g] = b->pos.x; c->y[g] = b->pos.y; c->z[g] = b->pos.z; c->ext[g] = b->ext;
		if (!c->active) continue;					// only the compact copy
		if (b->ext > c->extcap) { c->dirty = 1; return; }
		h = CellBucket(c,CellCoord(b->pos.x,c->invL),CellCoord(b->pos.y,c->invL),CellCoord(b->pos.z,c->invL));
		if (h != c->bucket[g])
		{
			CellsUnlink(c,g);
			CellsLink(c,g,h);
		}
	}
}

/*****************************************************************************
 Residues (global indices, in increasing order) which can interact with
 residue i of chain ic; returns their number, the list is in c->cand.
 Returns -1 when they are so many that scanning all residues is cheaper.
 *****************************************************************************/
int CellsCandidates(struct s_polymer *p, int ic, int i)
{
	struct s_cells *c = p->cells;
	struct s_back *b = ((p+ic)->back)+i;
	int ix,iy,iz,dx,dy,dz,h,k,n=0,nb=0,seen,buckets[27],g,t;

	ix = CellCoord(b->pos.x,c->invL);
	iy = CellCoord(b->pos.y,c->invL);
	iz = CellCoord(b->pos.z,c->invL);

	for (dx=-1;dx<=1;++dx)
		for (dy=-1;dy<=1;++dy)
			for (dz=-1;dz<=1;++dz)
			{
				h = CellBucket(c,ix+dx,iy+dy,iz+dz);
				seen = 0;					// two neighbouring cells may share a bucket
				for (k=0;k<nb;++k) if (buckets[k]==h) seen = 1;
				if (seen) continue;
				buckets[nb++] = h;
				for (g=c->head[h];g>=0;g=c->next[g]) c->cand[n++] = g;
			}

	if (4*n > c->nres) return -1;

	// increasing order, as the plain loop over chains and residues (same sums)
	for (k=1;k<n;++k)
	{
		t = c->cand[k];
		for (g=k-1;g>=0 && c->cand[g]>t;--g) c->cand[g+1] = c->cand[g];
		c->cand[g+1] = t;
	}
	return n;
}

/*****************************************************************************
 Debug aid (echeck): number of residues whose compact copy or cell does not
 match their current coordinates
 *****************************************************************************/
int CellsCheck(struct s_polymer *p)
{
	struct s_cells *c = p->cells;
	struct s_back *b;
	int g,h,nbad=0;

	if (!c || c->dirty) return 0;
	for (g=0;g<c->nres;++g)
	{
		b = ((p+c->chain[g])->back)+c->res[g];
		if (c->x[g]!=b->pos.x || c->y[g]!=b->pos.y || c->z[g]!=b->pos.z || c->ext[g]!=b->ext) { nbad++; continue; }
		if (c->active)
		{
			h = CellBucket(c,CellCoord(b->pos.x,c->invL),CellCoord(b->pos.y,c->invL),CellCoord(b->pos.z,c->invL));
			if (h != c->bucket[g]) nbad++;
		}
	}
	return nbad;
}

