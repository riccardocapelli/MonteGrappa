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
 * observables.c
 *
 * Observables of aggregation, written every nprintobs steps:
 *  - clusters of chains (connected components of the graph of inter-chain
 *    contacts): number, largest and mean size
 *  - pairs of chains in contact, classified by the cosine of their end-to-end
 *    vectors as parallel (> OBS_COSMIN), antiparallel (< -OBS_COSMIN) or other
 *  - register of parallel and antiparallel pairs: from the residues in
 *    contact, the mean shift between the residue numbers (0 = in register)
 *  - nematic order parameter of the end-to-end vectors of the chains
 * They only use the contact lists and vectors within the chains, so they
 * hold also in a periodic box.
 */

#include "montegrappa.h"

#define OBS_COSMIN 0.8		// |cos| above which two chains are (anti)parallel
#define OBS_INREG  0.5		// |register| below which a pair is in register

/*****************************************************************************
 Clusters of chains: connected components of the graph whose edges are the
 contacts between different chains (a chain without such contacts is a
 cluster by itself). Returns the number of clusters; the chains of cluster
 k are member[start[k]..start[k+1]-1].
 *****************************************************************************/
int FindClusters(struct s_polymer *p, int npol, int *member, int *start, int *cluster_of)
{
	int ic,q,i,k,c2,head,tail=0,nc=0;
	struct s_back *b;

	for (ic=0;ic<npol;++ic) cluster_of[ic] = -1;

	for (ic=0;ic<npol;++ic)
		if (cluster_of[ic] < 0)
		{
			start[nc] = tail;
			cluster_of[ic] = nc;
			member[tail++] = ic;
			for (head=start[nc];head<tail;++head)				// breadth-first visit
			{
				q = member[head];
				for (i=0;i<(p+q)->nback;++i)
				{
					b = ((p+q)->back)+i;
					for (k=0;k<b->ncontacts;++k)
					{
						c2 = (b->contacts_p)[k];
						if (cluster_of[c2] < 0)
						{
							cluster_of[c2] = nc;
							member[tail++] = c2;
						}
					}
				}
			}
			nc++;
		}
	start[nc] = tail;
	return nc;
}

// unit vector from the first to the last backbone atom of chain ic
static struct vector ChainDirection(struct s_polymer *p, int ic)
{
	struct vector u;
	struct s_back *b0 = (p+ic)->back, *b1 = ((p+ic)->back)+(p+ic)->nback-1;
	double n;

	u.x = b1->pos.x - b0->pos.x;
	u.y = b1->pos.y - b0->pos.y;
	u.z = b1->pos.z - b0->pos.z;
	n = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
	if (n > 0) { u.x /= n; u.y /= n; u.z /= n; }
	return u;
}

void PrintObservablesHeader(FILE *fp)
{
	fprintf(fp,"# step  nclusters  largest  mean_size  pairs_in_contact  parallel  antiparallel  "
		"parallel_in_register  antiparallel_in_register  mean_register_parallel  mean_register_antiparallel  nematic_S\n");
}

void PrintObservables(FILE *fp, struct s_polymer *p, int npol, unsigned long long istep)
{
	int member[npol],start[npol+1],cluster_of[npol],sumreg_n[npol];
	double sumreg[npol],cosab;
	struct vector u[npol];
	int nc,k,largest=0,a,b,i,q,first_a,first_b,len_b;
	int npairs=0,npar=0,nanti=0,nparreg=0,nantireg=0;
	double regpar=0,reganti=0,reg,Q[3][3]={{0}},S,ev,m,qq,pp,phi,r;
	struct s_back *ba;

	// clusters
	nc = FindClusters(p,npol,member,start,cluster_of);
	for (k=0;k<nc;++k) if (start[k+1]-start[k] > largest) largest = start[k+1]-start[k];

	for (a=0;a<npol;++a) u[a] = ChainDirection(p,a);

	// pairs of chains in contact: orientation and register (from the residue numbers in contact)
	for (a=0;a<npol;++a)
	{
		for (b=0;b<npol;++b) { sumreg[b] = 0; sumreg_n[b] = 0; }
		first_a = ((p+a)->back)->iaa;
		for (i=0;i<(p+a)->nback;++i)
		{
			ba = ((p+a)->back)+i;
			for (k=0;k<ba->ncontacts;++k)
			{
				b = (ba->contacts_p)[k];
				if (b <= a) continue;
				first_b = ((p+b)->back)->iaa;
				len_b = (((p+b)->back)+(p+b)->nback-1)->iaa - first_b + 1;
				q = (((p+b)->back)+(ba->contacts)[k])->iaa - first_b;
				cosab = u[a].x*u[b].x + u[a].y*u[b].y + u[a].z*u[b].z;
				if (cosab >= 0) sumreg[b] += q - (ba->iaa - first_a);					// parallel: residue k faces k
				else sumreg[b] += q - ((len_b-1) - (ba->iaa - first_a));				// antiparallel: k faces len-1-k
				sumreg_n[b]++;
			}
		}
		for (b=a+1;b<npol;++b)
			if (sumreg_n[b] > 0)
			{
				npairs++;
				cosab = u[a].x*u[b].x + u[a].y*u[b].y + u[a].z*u[b].z;
				reg = sumreg[b]/sumreg_n[b];
				if (cosab > OBS_COSMIN)
				{
					npar++; regpar += reg;
					if (fabs(reg) < OBS_INREG) nparreg++;
				}
				else if (cosab < -OBS_COSMIN)
				{
					nanti++; reganti += reg;
					if (fabs(reg) < OBS_INREG) nantireg++;
				}
			}
	}

	// nematic order parameter: largest eigenvalue of Q = < (3 u u - 1)/2 >
	for (a=0;a<npol;++a)
	{
		double v[3] = {u[a].x,u[a].y,u[a].z};
		for (i=0;i<3;++i)
			for (k=0;k<3;++k)
				Q[i][k] += (1.5*v[i]*v[k] - (i==k ? 0.5 : 0.)) / npol;
	}
	// eigenvalues of the symmetric traceless 3x3 matrix (trigonometric solution)
	pp = Q[0][1]*Q[0][1] + Q[0][2]*Q[0][2] + Q[1][2]*Q[1][2];
	m = (Q[0][0]+Q[1][1]+Q[2][2])/3.;
	qq = ( (Q[0][0]-m)*(Q[0][0]-m) + (Q[1][1]-m)*(Q[1][1]-m) + (Q[2][2]-m)*(Q[2][2]-m) + 2.*pp ) / 6.;
	if (qq <= 0) S = m;
	else
	{
		double B[3][3],det;
		for (i=0;i<3;++i) for (k=0;k<3;++k) B[i][k] = (Q[i][k] - (i==k ? m : 0.)) / sqrt(qq);
		det = B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1]) - B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0]) + B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]);
		r = det/2.;
		if (r <= -1) phi = PI/3.;
		else if (r >= 1) phi = 0;
		else phi = acos(r)/3.;
		ev = m + 2.*sqrt(qq)*cos(phi);
		S = ev;
	}

	fprintf(fp,"%llu\t%d\t%d\t%lf\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",istep,nc,largest,(double)npol/nc,npairs,npar,nanti,
		nparreg,nantireg,npar ? regpar/npar : 0.,nanti ? reganti/nanti : 0.,S);
	fflush(fp);
}
