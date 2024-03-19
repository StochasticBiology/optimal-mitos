#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND drand48()

// parameter set for a simulation
typedef struct {
  double D, kon, koff, rhoon, rhooff, activep, V, dmito, kmito, inter;
} Params;

// produce gaussian random number
double gsl_ran_gaussian(const double sigma)
{
  double x, y, r2;

  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      x = -1 + 2 * RND;
      y = -1 + 2 * RND;

      /* see if it is in the unit circle */
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  /* Box-Muller transform */
  return sigma * y * sqrt (-2.0 * log (r2) / r2);
}

// create an adjacency matrix matching given n, e statistics from a diffusive physical simulation
void AMFromSimulation(Params P, int nsample, double *meanmin, double *meanedges, int output, char *fname)
{
  // physical parameters -- somewhat arbitrary
  double cx = 100, cy = 30; // cell size
  double thresh = 1;        // encounter threshold
  int n = 140;              // number of mitos
  double idx, idy, ix, iy;
  double *x, *y;
  double *dx, *dy;
  int *active;
  int *am;
  int i, j;
  long int timer;
  double *scale;
  FILE *fp, *fp2;
  double r2;
  int sample;
  double mindist;
  int edges;
  char fname2[200];

  // initialise output files
  sprintf(fname2, "%s-am.csv", fname);
  if(output != 0)
    {
      fp = fopen(fname, "w");
      fprintf(fp, "frame,mito,x,y\n");
      fp2 = fopen(fname2, "w");
      fprintf(fp2, "frame,mito1,mito2\n");
    }

  timer = 0;
  // allocate memory for positions and convenient matrix
  x = (double*)malloc(sizeof(double)*n);
  y = (double*)malloc(sizeof(double)*n);
  dx = (double*)malloc(sizeof(double)*n);
  dy = (double*)malloc(sizeof(double)*n);
  active = (int*)malloc(sizeof(int)*n);
  am = (int*)malloc(sizeof(int)*n*n);
  scale = (double*)malloc(sizeof(double)*n);

  // loop over different random instances
  for(sample = 0; sample < nsample; sample++)
    {
      if(sample % 10 == 0)
	printf("%s sample %i\n", fname, sample);
      // initialise empty encounter matrix and random positions
      for(i = 0; i < n*n; i++)
	am[i] = 0;
      for(i = 0; i < n; i++)
	{
	  x[i] = RND*cx;
	  y[i] = RND*cy;
	  dx[i] = dy[i] = 0;
	  active[i] = (i < P.activep*n);
	}
      // loop while we have fewer than desired encounters
      timer = 0;
      while(timer < 1e2)
	{
	  timer++;

	  for(i = 0; i < n; i++)
	    scale[i] = 1;
      
	  // look for encounters
	  for(i = 0; i < n; i++)
	    {
	      for(j = i+1; j < n; j++)
		{
		  if(active[i] && active[j])
		    {
		      // distance between i and j
		      r2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
		      // if we're below a threshold, scale motion accordingly
		      if(r2 < P.dmito*P.dmito)
			{
			  scale[i] = P.kmito;
			  scale[j] = P.kmito;
			}
		      // if we're below a threshold and haven't scored this encounter yet, do so
		      if(am[i*n+j] == 0)
			{
			  if(r2 < thresh*thresh)
			    {
			      // populate encounter matrix and adj mat
			      am[i*n+j] = am[j*n+i] = 1;
			      if(output != 0 && sample == 0)
				{
				  fprintf(fp2, "%li,%i,%i\n", timer, i, j);
				}
			   
			    }
			}
		    }
		}
	    }
	  // put a uniform diffusion kernel on each position
	  for(i = 0; i < n; i++)
	    {
	      // activate / inactivate mito (modelling motion in and out of plane)
	      if(active[i])
		{
		  if(RND < P.rhooff) active[i] = 0;
		}
	      else
		{
		  if(RND < P.rhoon) active[i] = 1;
		}
	      if(active[i])
		{
		  // apply "hydrodynamic" interactions if appropriate
		  if(P.inter)
		    {
		      idx = 0; idy = 0;
		      for(j = 0; j < n; j++)
			{
			  if((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) < fabs(P.inter))
			    {
			      idx += (x[i]-x[j]);
			      idy += (y[i]-y[j]);
			    }
			}
		      if(idx < 0) ix = -1*(P.inter < 0 ? -1 : 1); else ix = 1*(P.inter < 0 ? -1 : 1);
		      if(idy < 0) iy = -1*(P.inter < 0 ? -1 : 1); else iy = 1*(P.inter < 0 ? -1 : 1);
		    }
		  // attach/detach from cytoskeleton
		  if(dx[i] == 0 && dy[i] == 0 && RND < P.kon)
		    {
		      // choose horizontal or vertical motion randomly
		      if(RND < 0.5) { dx[i] = (P.inter ? ix : (RND < 0.5 ? -1 : 1))*P.V; dy[i] = 0; }
		      else { dx[i] = 0; dy[i] = (P.inter ? iy : (RND < 0.5 ? -1 : 1))*P.V; }
		    }
		  else if((dx[i] != 0 || dy[i] != 0) && RND < P.koff)
		    {
		      // back to diffusion
		      dx[i] = 0; dy[i] = 0;
		    }

		  // apply diffusion kernel if we're not ballistic
		  if(dx[i] == 0 && dy[i] == 0)
		    {
		      x[i] += (gsl_ran_gaussian(2.*P.D)+ix*P.D)*scale[i];
		      y[i] += (gsl_ran_gaussian(2.*P.D)+iy*P.D)*scale[i];
		    }
		  // otherwise ballistic motion
		  else
		    {
		      x[i] += dx[i]*scale[i];
		      y[i] += dy[i]*scale[i];
		    }
		  // reflecting boundaries
		  if(x[i] < 0) x[i] = 0;
		  if(x[i] > cx) x[i] = cx;
		  if(y[i] < 0) y[i] = 0;
		  if(y[i] > cy) y[i] = cy;
		}
	    }
	  // output state to file if desired
	  if(output != 0 && sample == 0)
	    {
	      for(i = 0; i < n; i++)
		fprintf(fp, "%li,%i,%.3f,%.3f\n", timer, i, x[i], y[i]);
	    }
	}

      // record statistics of mean minimum distance and edges in adj mat
      meanmin[sample] = meanedges[sample] = 0;
      for(i = 0; i < n; i++)
	{
	  mindist = -1; edges = 0;
	  for(j = 0; j < n; j++)
	    {
	      if(i != j)
		{
		  r2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]);
		  if(r2 < mindist || mindist == -1)
		    mindist = r2;
		}
	      if(am[i*n+j]) edges++;
	    }
	  // add to growing mean score
	  meanmin[sample] += sqrt(mindist);
	  meanedges[sample] += (edges/2);
	}
      // normalise mean score
      meanmin[sample] /= n;
      meanedges[sample] /= n;
    }
  
  // free up memory and close up
  free(x); free(y); free(dx); free(dy); free(am); free(scale); free(active);
  if(output != 0)
    {
      fclose(fp);						    
      fclose(fp2);
    }
}

int main(int argc, char *argv[])
{
  Params P;
  int i;
  FILE *fp;
  double *meanmin, *meanedges;
  int ns = 100;
  int expt = 0;
  int fullsim;

  // from command line, are we running the full parameter sweep?
  fillsim = 0;
  if(argc == 2)
    {
      if(strcmp(argv[1], "--full") == 0)
	fullsim = 1;
    }
  
  meanmin = (double*)malloc(sizeof(double)*ns*100);
  meanedges = (double*)malloc(sizeof(double)*ns*100);

  // diffusion is 0.1um2 s-1, frame is 2s, so 0.2um2 frame-1

  double scaled = 0.2, scalev = 0.5;

  // here we go through a collection of specifically-chosen parameterisations, simulating behaviour and outputting trajectories to files
  P.inter = 0;
  // 0
  P.D = scaled*1; P.kon = 0; P.koff = 0; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 0; P.kmito = 1;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test0.csv"); expt++;
  // 1
  P.D = scaled*1; P.kon = 0.1; P.koff = 0.1; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 0; P.kmito = 1;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test1.csv"); expt++;
  // 2
  P.D = scaled*1; P.kon = 0.1; P.koff = 0.1; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*2; P.dmito = 0; P.kmito = 1;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test2.csv"); expt++;
  // 3
  P.D = scaled*1; P.kon = 0; P.koff = 0; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 5; P.kmito = 0.5;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test3.csv"); expt++;
  // 4
  P.D = scaled*1; P.kon = 0; P.koff = 0; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 5; P.kmito = 2;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test4.csv"); expt++;
  // 5
  P.D = scaled*0.25; P.kon = 0; P.koff = 0; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 5; P.kmito = 4;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test5.csv"); expt++;
  // 6
  P.D = scaled*0.5; P.kon = 0; P.koff = 0; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 5; P.kmito = 2;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test6.csv"); expt++;
  // 7
  P.D = scaled*1; P.kon = 1; P.koff = 0.1; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 2; P.kmito = 0.5;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test7.csv"); expt++;
  
  P.inter = 30;
  // 8
  P.D = scaled*1; P.kon = 0; P.koff = 0; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 0; P.kmito = 1;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test8.csv"); expt++;

  P.inter = -1000;
  // 9
  P.D = scaled*1; P.kon = 0; P.koff = 0; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 0; P.kmito = 1;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test9.csv"); expt++;

  P.inter = -20;
  // 10
  P.D = scaled*1; P.kon = 0; P.koff = 0; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*1; P.dmito = 0; P.kmito = 1;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test10.csv"); expt++;

  P.inter = 0;
  // 11
  P.D = scaled*1; P.kon = 0.1; P.koff = 0.1; P.rhoon = 1; P.rhooff = 0; P.activep = 1; P.V = scalev*5; P.dmito = 0; P.kmito = 1;
  AMFromSimulation(P, ns, &(meanmin[ns*expt]), &(meanedges[ns*expt]), 1, "test11.csv"); expt++;
 
  fp= fopen("outstats.csv", "w");
  fprintf(fp, "params,sample,meanmin,meanedges\n");
  for(i = 0; i < ns*expt; i++)
    fprintf(fp, "%i,%i,%f,%f\n", i/ns, i%ns, meanmin[i], meanedges[i]);

  double scaled = 0.2, scalev = 0.5;

  if(fullsim == 1)
    {
      fp= fopen("outstatsscan.csv", "w");
      fprintf(fp, "D,kon,koff,V,dmito,kmito,expt,sample,meanmin,meanedges\n");
      fclose(fp);
  
      P.rhoon = 1; P.rhooff = 0; P.activep = 1;
      for(P.D = 0.05; P.D <= 0.21; P.D *= 2) {
	for(P.kon = 0; P.kon <= 1.01; P.kon += 0.25) {
	  for(P.koff = 0; P.koff <= 1.01; P.koff += 0.25) {
	    for(P.V = 0.5; P.V <= 4.01; P.V *= 2) {
	      for(P.dmito = 1; P.dmito <= 4.01; P.dmito *= 2) {
		for(P.kmito = 0.25; P.kmito <= 4.01; P.kmito *= 2) {

		  if(P.kmito * P.D <= 1.01)
		    {
		      sprintf(str, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f", P.D, P.kon, P.koff, P.V, P.dmito, P.kmito);
		      AMFromSimulation(P, ns, meanmin, meanedges, 0, str); expt++;
		      fp= fopen("outstatsscan.csv", "a");
  
		      for(i = 0; i < ns; i++)
			fprintf(fp, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%i,%i,%f,%f\n", P.D, P.kon, P.koff, P.V, P.dmito, P.kmito, expt, i, meanmin[i], meanedges[i]);
		      fclose(fp);
		    }
		}
	      }
	    }
	  }
	}
      }
    }
 
  return 0;
}
