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

#include "montegrappa.h"



struct s_polymer *Allo_Fragment(int npol, int nback,int nangmax, FILE *flog)
{
	
	int ipol;
	struct s_polymer *f;

	f=(struct s_polymer *)calloc(npol,sizeof(struct s_polymer));
	if(!f)      Error("Cannot Allocate fragment");


	for(ipol=0;ipol<npol;++ipol)
	{
		
		(f+ipol)->nback=nback;
		(f+ipol)->back=(struct s_back *)calloc((nback+3),sizeof(struct s_back));
		if(!(f+ipol)->back)     Error("\tAlloFragment(): Cannot Allocate backbone");

                (f+ipol)->A=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->G=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->L=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->Y=AlloDoubleMatrix(nangmax,nangmax);
                (f+ipol)->g_ang=AlloDouble(nangmax);      
                (f+ipol)->d_ang=AlloDouble(nangmax);
	
	}

	return f;
	
}

double DetTriang( double **Mat, int n )
{

	int i;
	double det=1.0;
	
	for(i=0; i<n; ++i)
   		det = det * Mat[i][i];

	return det;

}


void FreeFragment(struct s_polymer *fragment,struct s_mc_parms *parms)
{

	int i;
	if (!fragment) return;
	for(i=0;i<parms->npol;++i)
        {
	        FreeDoubleMatrix((fragment+i)->A,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->G,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->L,parms->nmul_local);
                FreeDoubleMatrix((fragment+i)->Y,parms->nmul_local);
                free((fragment+i)->d_ang );
                free((fragment+i)->g_ang);
		free((fragment+i)->back);	
	}

	free(fragment);

}


void InvertTriang( double **Inv, double **Mat, int n )
{
	int i,k,a;
	
	for(i=0; i<n; ++i)
   		for(k=0; k<n; ++k)
     			Inv[i][k]=0.0;

	for(i=0; i<n; ++i)
   		Inv[i][i] = 1./Mat[i][i];

	for(i=0; i<n; ++i)
  		for(a=1; a<=i; a++)
  		{
     			for( k= i-a; k<i; ++k)
       				Inv[i][i-a] += (Mat[i][k] * Inv[k][i-a]);
       			Inv[i][i-a] /= -Mat[i][i];
  		}

	return;

}




double Squared_n_Norma( double *vect, int dim )
{
	int i;
	double norm=0.0;

	for(i=0; i<dim; ++i)
    		norm += vect[i]*vect[i];
	
	return norm;
}

void TransposedMatOnVect( double **Mat, double *vect, double *risu, int n)
{
       int i,j;

      for(i=0; i<n; ++i)
      {
            risu[i] = 0.0;
            for(j=0; j<n; ++j)
                  risu[i] += Mat[j][i] * vect[j];
      }

      return;
}




void MatA( double **A, double **G, int dim,double a,double b)
{
     int i, j;
      
      for(i=0; i<dim; i++)
            for(j=0; j<dim; j++)
            { 
		if(i==j)
                	A[i][i]= 0.5*a*(1.+ b*G[i][j]);
                else
                	A[i][j]= 0.5*a*b*G[i][j]; 
              }
     return;
		   
		    
}


void Cholesky_2( double **L, double **A, int dim)
{
      int i,j,k;
      double sum=0.0;
      
      for( i=0; i<dim; i++)
            for( j=i; j<dim; j++)
            {

	          sum=A[i][j];
                  for(k=0;k<i;k++)
                        sum += -L[i][k]*L[j][k];
                  if(i == j)
                  {
                        if(sum<=0.0)
                              Error("Cholesky decomposition failed (check the BGS parameters!)");
                        L[i][i]= sqrt(sum);
                  }
                  else
                        L[j][i]= sum/L[i][i];
            }
            
      return;
}
      











      















	  
/*****************************************************************************
 Biased Gaussian steps (Favrin, Irbaeck, Sjunnesson, J Chem Phys 114 (2001) 8154)

 The free dihedrals kd[0..nang-1] of the segment [iw,iw+nmul-1] are changed
 by dphi ~ det(L) exp(-1/2 dphi^T A dphi), A = L L^T = a/2 (1 + b G), where
 G_ij = sum_e dr_e/dphi_i . dr_e/dphi_j over the three last moved atoms
 iw+nmul-3..iw+nmul-1, which are bound to the fixed rest of the chain.
 Dihedral k is the rotation around the bond (k-1,k), which moves the atoms
 k+1..iw+nmul-1 (PivotForward(k-1)); angles are in degrees.
 *****************************************************************************/
static void BGS_GMatrix(struct s_polymer *p, int ip, int iw, int nmul, int *kd, int nang, double **G)
{
	int m,n,e,c;
	double der[nang][3][3],norm;
	struct vector ax,r;
	struct s_back *b = (p+ip)->back;

	for (m=0;m<nang;++m)
	{
		ax.x = (b+kd[m])->pos.x - (b+kd[m]-1)->pos.x;
		ax.y = (b+kd[m])->pos.y - (b+kd[m]-1)->pos.y;
		ax.z = (b+kd[m])->pos.z - (b+kd[m]-1)->pos.z;
		norm = sqrt(Norm2(ax)) * 180. / PI;			// derivative per degree
		ax.x /= norm; ax.y /= norm; ax.z /= norm;

		for (e=0;e<3;++e)
		{
			if (iw+nmul-3+e < kd[m]+1)				// atom not moved by this dihedral
			{
				for (c=0;c<3;++c) der[m][e][c] = 0.;
				continue;
			}
			r.x = (b+iw+nmul-3+e)->pos.x - (b+kd[m])->pos.x;
			r.y = (b+iw+nmul-3+e)->pos.y - (b+kd[m])->pos.y;
			r.z = (b+iw+nmul-3+e)->pos.z - (b+kd[m])->pos.z;
			der[m][e][0] = ax.y * r.z - ax.z * r.y;
			der[m][e][1] = ax.z * r.x - ax.x * r.z;
			der[m][e][2] = ax.x * r.y - ax.y * r.x;
		}
	}

	for (m=0;m<nang;++m)
		for (n=0;n<nang;++n)
		{
			G[m][n] = 0.;
			for (e=0;e<3;++e)
				for (c=0;c<3;++c)
					G[m][n] += der[m][e][c] * der[n][e][c];
		}
}

/*****************************************************************************
 Weight of the proposal, given g = L^T dphi: det(L) exp(-|g|^2/2)
 *****************************************************************************/
static double BGS_Weight(double **L, double *g, int nang)
{
	return DetTriang(L,nang) * exp(-0.5*Squared_n_Norma(g,nang));
}

int localpivot(int ip,int iw,struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
  
  
  double deltaE,W1=0.,W2=0.;

  int nang,i,m,ok=1;

  if (iw < nmul+3)
  {
    if ( (((p+ip)->back)+iw+1)->move == 0 ) { return 0;}

    deltaE = -GetEnergyMonomerRange(p,0,iw,ip);							// two-body energy
  
    if (!parms->nodihpot)
      deltaE -= (((p+ip)->back)+iw+1)->e_dih;							// old dihedral energy
    if (!parms->nohfields)
      deltaE -= GetHFields(p,parms->npol);										// old H fields energy
  
    ok = PivotBackward((p+ip),iw+1,0.2*(0.5-frand()),iw,parms);							// moved atoms are in [0,iw-1]; changed dihedral is iw+1

    ok *= AddSidechain(p,0,iw,ip);
    
    deltaE += EnergyMonomerRange(p,pot,0,iw,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);			// new energy of iw with the others
    
    if (!parms->nodihpot)
      deltaE += EnergyDihedrals(p,pot,iw+1,ip,1);
    if (!parms->nohfields)
	deltaE += EnergyHFieldsAll(p,pot,parms);

    ok *= Metropolis(deltaE,t,p->tables);

    if (ok == 0) 			// move rejected
    {
      UpdateMonomerRange(oldp,p,0,iw,ip,0);                         
      return 0;	
    }
    else
    {
      UpdateMonomerRange(p,oldp,0,iw,ip,0); 
      parms->acc++;
      p->etot+=deltaE;
      return 1;
    }
	
    return 0;
    
  }
  
  
  else if ( iw > ( (p+ip)->nback -3 -nmul   )  ) // do small forward
  { 
     if ( (((p+ip)->back)+iw)->move == 0 ) { return 0; }

     deltaE=-GetEnergyMonomerRange(p,iw,(p+ip)->nback-1,ip);
     
     if (!parms->nodihpot)
        deltaE -= (((p+ip)->back)+iw)->e_dih;		
    if (!parms->nohfields)
       deltaE -= GetHFields(p,parms->npol);

    ok = PivotForward((p+ip),iw-1,0.2*(0.5-frand()),(p+ip)->nback-iw-1,parms);			// moved atoms are in [iw+1,nback-1]; changed dihedral of iw
   
    ok *= AddSidechain(p,iw,(p+ip)->nback-1,ip);
    
    deltaE += EnergyMonomerRange(p,pot,iw,(p+ip)->nback-1,ip,parms->npol,0,1,parms->nosidechains,parms->disentangle,parms->hb);
   
    if (!parms->nodihpot)
	deltaE += EnergyDihedrals(p,pot,iw,ip,1);
   
    if (!parms->nohfields)
    {
      deltaE += EnergyHFieldsAll(p,pot,parms);
      
    }

    ok *= Metropolis(deltaE,t,p->tables);
    
    if (ok == 0) 			// move rejected
    {
      UpdateMonomerRange(oldp,p,iw,(p+ip)->nback-1,ip,0);                           // ... if moved second half
      return 0;	
    }
    else
    {
      UpdateMonomerRange(p,oldp,iw,(p+ip)->nback-1,ip,0); 
      parms->acc++;
      p->etot+=deltaE;
      return 1;
    }

    
  } // end of small forward
  
  
  else {

    struct s_back *b = (p+ip)->back;
    struct s_polymer *f = fragment+ip;
    int kd[nmul],k;
    double rnd;

    // free dihedrals of the segment (see BGS_GMatrix)
    nang=0;
    for (k=iw+1;k<=iw+nmul-2;++k)
      if ((b+k)->move) kd[nang++] = k;
    if (nang==0) return 0;

    deltaE=-GetEnergyMonomerRange(p,iw,iw+nmul,ip);

    if (!parms->nodihpot)
      for(i=0;i<nmul+3;i++)
        deltaE -= (b+iw+i)->e_dih;							// old dihedral energy
    if (!parms->noangpot)
    {
      deltaE -= (b+iw+nmul-1)->e_ang;							// old angular energy at the closure
      deltaE -= (b+iw+nmul)->e_ang;
    }
    if (!parms->nohfields)
      deltaE -= GetHFields(p,parms->npol);

    // draw g ~ N(0,1) and set dphi = L^-T g, so that g = L^T dphi
    for (m=0;m<nang;++m) f->g_ang[m] = gauss(0.,1.);
    BGS_GMatrix(p,ip,iw,nmul,kd,nang,f->G);
    MatA(f->A,f->G,nang,parms->bgs_a,parms->bgs_b);
    Cholesky_2(f->L,f->A,nang);
    W1 = BGS_Weight(f->L,f->g_ang,nang);
    InvertTriang(f->Y,f->L,nang);
    TransposedMatOnVect(f->Y,f->g_ang,f->d_ang,nang);

    for (m=0;m<nang;++m)
      ok *= PivotForward((p+ip),kd[m]-1,f->d_ang[m],iw+nmul-1-kd[m],parms);

    // the bond, angles and locked dihedrals that join the segment to the rest of the chain
    if (DAbs( Dist2( (b+iw+nmul-1)->pos, (b+iw+nmul)->pos ) - (b+iw+nmul-1)->d2_next ) > 0.05) ok=0;
    if (ok) ok = CheckAngleRange(p,ip,iw+nmul-1,iw+nmul,2.0);
    if (ok) ok = CheckLockedDihedralRange(p,ip,iw+nmul-1,iw+nmul+1,2.0);

    if (ok && !parms->nosidechains)
      ok = AddSidechain(p,iw,iw+nmul,ip);

    if (ok)
    {
      deltaE += EnergyMonomerRange(p,pot,iw,iw+nmul,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);

      if (!parms->nodihpot)
        for(i=0;i<nmul+3;i++)
          deltaE+=EnergyDihedrals(p,pot,iw+i,ip,1);
      if (!parms->noangpot)
      {
        deltaE += EnergyAngles(p,pot,iw+nmul-1,ip,1);
        deltaE += EnergyAngles(p,pot,iw+nmul,ip,1);
      }
      if (!parms->nohfields)
        deltaE += EnergyHFieldsAll(p,pot,parms);

      // weight of the reverse step (-dphi) from the new conformation
      BGS_GMatrix(p,ip,iw,nmul,kd,nang,f->G);
      MatA(f->A,f->G,nang,parms->bgs_a,parms->bgs_b);
      Cholesky_2(f->L,f->A,nang);
      TransposedMatOnVect(f->L,f->d_ang,f->g_ang,nang);
      W2 = BGS_Weight(f->L,f->g_ang,nang);

      // the proposal is not symmetric: the ratio of weights enters also when deltaE<=0
      rnd = (W2/W1) * exp(-deltaE/t);
      if (rnd < 1. && frand() >= rnd) ok = 0;
    }

    if (ok)
    {
      UpdateMonomerRange(p,oldp,iw,iw+nmul,ip,parms->shell);
      parms->acc++;
      p->etot += deltaE;
      return 1;
    }

    UpdateMonomerRange(oldp,p,iw,iw+nmul,ip,parms->shell);
    return 0;

  }// end of local forward
  
  
  return 0; 
  
}
	  

int LocalMove(struct s_polymer *p, struct s_polymer *oldp,struct s_polymer *fragment,struct s_potential *pot,int nmul,struct s_mc_parms *parms, double t)
{
        int ip,iw;
	int ok=1;
	
	
	ip   = irand(parms->npol);
	iw   = irand( (p+ip)->nback );
		
	ok=localpivot(ip,iw,p,oldp,fragment,pot,nmul,parms,t);  	

            	
	parms->mov++;
        
	return ok;

}	


/*****************************************************************************
 Backrub move: rotate the five backbone atoms iw-2..iw+2 (iw is a CA) around
 the axis joining the CAs iw-3 and iw+3. The bond angles at iw-3 and iw+3 and
 the dihedrals iw-3, iw-2, iw+3, iw+4 change. If newrot, a new rotamer is also
 chosen for residue iw (backsiderub).
 *****************************************************************************/
static int BackRubMove(struct s_polymer *p, struct s_polymer *oldp,struct s_potential *pot,struct s_mc_parms *parms,double t,int newrot)
{
  int ok=1;
  int iw,ip;
  int out;

  double dw,deltaE;

  ip=irand( parms->npol);
  iw=3*irand((p+ip)->nback/3)+1;
  dw = 2.0 * (0.5 - frand());

  // the angle checks need the backbone atoms iw-4 and iw+4
  if (iw<4 || iw+4>(p+ip)->nback-1)
    return 0;

  parms->mov++;

  deltaE = -GetEnergyMonomerRange(p,iw-3,iw+3,ip);                                                 // two-body energy

  if (!parms->nodihpot)
  {
      deltaE -= (((p+ip)->back)+iw-3)->e_dih;                                                   // old dihedral energy
      deltaE -= (((p+ip)->back)+iw-2)->e_dih;
      deltaE -= (((p+ip)->back)+iw+3)->e_dih;
      deltaE -= (((p+ip)->back)+iw+4)->e_dih;
  }
  if (!parms->noangpot)
  {
      deltaE -= (((p+ip)->back)+iw-3)->e_ang;                                                   // old angular energy
      deltaE -= (((p+ip)->back)+iw+3)->e_ang;
  }
  if (!parms->nohfields)
      deltaE -= GetHFields(p,parms->npol);

  ok=BackFlip(p+ip,iw,dw);

  if (ok)
    if( DAbs( Angle( (((p+ip)->back)+iw-4)->pos, (((p+ip)->back)+iw-3)->pos, (((p+ip)->back)+iw-2)->pos,p->tables,&out)-(((p+ip)->back)+iw-3)->a_next) > 5.0)
      ok=0;
  if (ok)
    if( DAbs( Angle( (((p+ip)->back)+iw+2)->pos, (((p+ip)->back)+iw+3)->pos, (((p+ip)->back)+iw+4)->pos,p->tables,&out)-(((p+ip)->back)+iw+3)->a_next) > 5.0)
      ok=0;

  if (ok)
  {
    if (newrot)
      (((p+ip)->back)+iw)->irot = irand( (((p+ip)->back)+iw)->nrot );

    if(!parms->nosidechains)
      ok *= AddSidechain(p,iw-3,iw+3,ip);
  }

  if (ok)
  {
    deltaE += EnergyMonomerRange(p,pot,iw-3,iw+3,ip,parms->npol,parms->shell,1,parms->nosidechains,parms->disentangle,parms->hb);

    if (!parms->nodihpot)
    {
      deltaE += EnergyDihedrals(p,pot,iw-3,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw-2,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw+3,ip,1);
      deltaE += EnergyDihedrals(p,pot,iw+4,ip,1);
    }
    if (!parms->noangpot)
    {
      deltaE += EnergyAngles(p,pot,iw-3,ip,1);
      deltaE += EnergyAngles(p,pot,iw+3,ip,1);
    }
    if (!parms->nohfields)
      deltaE += EnergyHFieldsAll(p,pot,parms);

    ok *= Metropolis(deltaE,t,p->tables);
  }

  if (ok == 0)                        // move rejected
  {
    UpdateMonomerRange(oldp,p,iw-3,iw+3,ip,0);
    return 0;
  }

  UpdateMonomerRange(p,oldp,iw-3,iw+3,ip,0);
  parms->acc++;
  p->etot+=deltaE;
  return 1;
}


int BackRub(struct s_polymer *p, struct s_polymer *oldp,struct s_potential *pot,struct s_mc_parms *parms,double t)
{
  return BackRubMove(p,oldp,pot,parms,t,0);
}


int BackSideRub(struct s_polymer *p, struct s_polymer *oldp,struct s_potential *pot,struct s_mc_parms *parms,double t)
{
  return BackRubMove(p,oldp,pot,parms,t,1);
}
