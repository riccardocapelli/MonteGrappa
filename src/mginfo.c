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
 * mginfo.c
 *
 * Progress of the simulation in the file mginfo (as the mdinfo file of
 * Amber), rewritten every tmginfo seconds of wall-clock time by rank 0:
 * run, step, energy, acceptance, speed, time elapsed and expected to the end
 * of the run and of the whole simulation. During the optimization of the
 * potential it shows the step of the optimization and the chi2.
 * The time to the end of the simulation comes from the runs already completed
 * (Monte Carlo and optimization); before the first one ends it counts the
 * Monte Carlo only. The file is written aside and renamed, so that it is never
 * read half written.
 */

#include "montegrappa.h"

static struct
{
	int on;
	double tstart;			// start of the simulation (of this execution, after a restart)
	double tlast;			// last time the file was written
	double trun;			// start of this run
	double tmc;			// end of the Monte Carlo of this run (0 while it goes on)
	double top;			// start of the optimization of the potential
	unsigned long long step0;	// first step of this run (a restart does not start from 0)
	int ndone;			// runs completed and their mean durations
	double dtrun,dtmc;
} mgi;

static double Now(void)
{
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC,&ts);
	return ts.tv_sec + 1e-9*ts.tv_nsec;
}

static void HMS(char *s, size_t size, double t)
{
	long n = t>0 ? (long)(t+0.5) : 0;

	if (n >= 86400) snprintf(s,size,"%ldd %02ld:%02ld:%02ld",n/86400,(n/3600)%24,(n/60)%60,n%60);
	else snprintf(s,size,"%02ld:%02ld:%02ld",n/3600,(n/60)%60,n%60);
}

// header, time to the end and expected end; body writes the lines of the phase
static void Write(struct s_mc_parms *parms, int irun, double left, const char *note, void (*body)(FILE *,void *), void *data)
{
	FILE *fp;
	char s1[32],s2[32],date[64];
	double t = Now(),total;
	time_t end;

	fp = fopen("mginfo.tmp","w");
	if (!fp) { mgi.on = 0; fprintf(stderr,"WARNING: cannot write mginfo, no more progress reports\n"); return; }

	end = time(NULL);
	strftime(date,sizeof(date),"%a %b %d %H:%M:%S %Y",localtime(&end));
	fprintf(fp,"MonteGrappa progress, %s\n\n",date);
	fprintf(fp,"run                  %d of %d\n",irun+1,parms->nrun);
	body(fp,data);

	// rest of this run, then the following ones as long as the completed ones
	total = left + (parms->nrun-irun-1) * (mgi.ndone ? mgi.dtrun : (mgi.tmc>0 ? mgi.tmc-mgi.trun : t-mgi.trun+left));
	HMS(s1,sizeof(s1),t-mgi.trun);
	HMS(s2,sizeof(s2),left);
	fprintf(fp,"elapsed (run)        %s\n",s1);
	fprintf(fp,"remaining (run)      %s%s\n",s2,note);
	HMS(s1,sizeof(s1),t-mgi.tstart);
	HMS(s2,sizeof(s2),total);
	fprintf(fp,"elapsed (all)        %s\n",s1);
	fprintf(fp,"remaining (all)      %s%s\n",s2,mgi.ndone ? "" : "   (Monte Carlo only until a run is completed)");
	end = time(NULL) + (time_t)total;
	strftime(date,sizeof(date),"%a %b %d %H:%M:%S %Y",localtime(&end));
	fprintf(fp,"expected end         %s\n",date);
	fclose(fp);
	rename("mginfo.tmp","mginfo");
	mgi.tlast = t;
}

void MGInfo_StartRun(struct s_mc_parms *parms, unsigned long long step0)
{
	mgi.on = parms->tmginfo > 0;
	mgi.trun = mgi.tlast = Now();
	if (mgi.tstart==0) mgi.tstart = mgi.trun;
	mgi.tmc = 0;
	mgi.step0 = step0;
}

struct mc_data { struct s_mc_parms *parms; unsigned long long istep; double e,T,speed; };

static void MCBody(FILE *fp, void *data)
{
	struct mc_data *d = data;
	struct s_mc_parms *parms = d->parms;

	fprintf(fp,"phase                Monte Carlo\n");
	fprintf(fp,"step                 %llu of %llu (%.1f%%)\n",d->istep,parms->nstep,100.*d->istep/parms->nstep);
	fprintf(fp,"temperature          %lf\n",d->T);
	fprintf(fp,"energy               %lf\n",d->e);
	fprintf(fp,"acceptance           %lf\n",parms->mov ? (double)parms->acc/parms->mov : 0.);
	fprintf(fp,"speed                %.0f steps/s\n",d->speed);
}

// from the step loop: writes when tmginfo seconds have passed (or always if force)
void MGInfo_MC(struct s_mc_parms *parms, int irun, unsigned long long istep, double e, double T, int force)
{
	struct mc_data d = { parms, istep, e, T, 0. };
	double t;

	if (!mgi.on) return;
	t = Now();
	if (!force && t-mgi.tlast < parms->tmginfo) return;
	if (t > mgi.trun && istep > mgi.step0) d.speed = (istep-mgi.step0) / (t-mgi.trun);
	if (d.speed > 0)
		Write(parms,irun,(parms->nstep-istep)/d.speed + (mgi.ndone ? mgi.dtrun-mgi.dtmc : 0.),
			mgi.ndone ? "" : "   (Monte Carlo only)",MCBody,&d);
	else							// start of the run: as long as the completed ones
		Write(parms,irun,mgi.ndone ? mgi.dtrun : 0.,mgi.ndone ? "" : "   (not known yet)",MCBody,&d);
}

static void EndMCBody(FILE *fp, void *data)
{
	fprintf(fp,"phase                end of the Monte Carlo, frames for the optimization\n");
}

// end of the Monte Carlo of the run; with the optimization, the frames are prepared next
void MGInfo_EndMC(struct s_mc_parms *parms, int irun)
{
	mgi.tmc = Now();
	if (mgi.on && strcmp(parms->op_minim,"none"))
		Write(parms,irun,mgi.ndone ? mgi.dtrun-mgi.dtmc : 0.,mgi.ndone ? "" : "   (not known yet)",EndMCBody,NULL);
}

struct op_data { int istep,nstep; double chi2,speed; };

static void OPBody(FILE *fp, void *data)
{
	struct op_data *d = data;

	fprintf(fp,"phase                optimization of the potential\n");
	fprintf(fp,"step                 %d of %d (%.1f%%)\n",d->istep,d->nstep,d->nstep ? 100.*d->istep/d->nstep : 100.);
	fprintf(fp,"chi2                 %lf\n",d->chi2);
	if (d->speed>0) fprintf(fp,"speed                %.1f steps/s\n",d->speed);
}

// from the optimization of the potential (rank 0), with its reduced chi2
// (step 0 is its start and is always written)
void MGInfo_Optimize(struct s_mc_parms *parms, int irun, int istep, int nstep, double chi2)
{
	struct op_data d = { istep, nstep, chi2, 0. };
	double t;

	if (!mgi.on) return;
	t = Now();
	if (istep==0) mgi.top = t;
	else if (t-mgi.tlast < parms->tmginfo) return;
	if (istep>0 && t>mgi.top) d.speed = istep / (t-mgi.top);
	Write(parms,irun,d.speed>0 ? (nstep-istep)/d.speed : 0.,d.speed>0 ? "" : "   (not known yet)",OPBody,&d);
}

struct end_data { int last; };

static void EndBody(FILE *fp, void *data)
{
	fprintf(fp,"phase                %s\n",((struct end_data *)data)->last ? "finished" : "end of run");
}

// at the end of each run: its durations enter the estimate of the following ones
void MGInfo_EndRun(struct s_mc_parms *parms, int irun)
{
	struct end_data d = { irun==parms->nrun-1 };
	double t = Now();

	if (!mgi.on) return;
	if (mgi.tmc==0) mgi.tmc = t;
	mgi.dtrun = (mgi.dtrun*mgi.ndone + (t-mgi.trun)) / (mgi.ndone+1);
	mgi.dtmc = (mgi.dtmc*mgi.ndone + (mgi.tmc-mgi.trun)) / (mgi.ndone+1);
	mgi.ndone ++;
	Write(parms,irun,0.,"",EndBody,&d);
}
