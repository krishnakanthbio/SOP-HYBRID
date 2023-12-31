// clang-format off

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Paolo Raiteri (Curtin University)
   
   Modified by : Krishnakanth B, Theoretical Biophysics Laboratory, 
   Molecular Biophysics Unit, Indian Institute of Science, Bangalore - 560012  
------------------------------------------------------------------------- */

#include "pair_lj_sop.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include <iostream>


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJ_SOP::PairLJ_SOP(LAMMPS *lmp) : Pair(lmp) {
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJ_SOP::~PairLJ_SOP()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_inner);
    memory->destroy(cut_inner_sq);
    memory->destroy(epsilon_lj);
    memory->destroy(sigma);
    
    memory->destroy(lj1); // LJ force factor 1
    memory->destroy(lj2); // LJ force factor 2
    memory->destroy(lj3); // LJ energy factor 3
    memory->destroy(lj4); // LJ energy factor 4

  }
}

/* ---------------------------------------------------------------------- */

void PairLJ_SOP::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *domain_id = atom->domain_id;       // domain tag of the atom
  tagint *residue_id = atom->residue_id;     // residue tag of the atom
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double rr, d, dd, tt, dt, dp, philj;
  tagint i_residue_id, j_residue_id, i_domain_id, j_domain_id;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  tagint *tag = atom->tag;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    i_residue_id = residue_id[i];
    i_domain_id = domain_id[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j_residue_id = residue_id[j];
      j_domain_id = domain_id[j];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        if ((i_domain_id != j_domain_id) || (i_domain_id == 0) || (j_domain_id == 0)){
//        std::cout << "acting\n";
        forcelj = r6inv*(lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        } else{ 
        forcelj = 0.0;
        }
        
        if (rsq > cut_inner_sq[itype][jtype]) {
          philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);

          rr = sqrt(rsq);
          dp = (cut[itype][jtype] - cut_inner[itype][jtype]);
          d = (rr-cut_inner[itype][jtype]) / dp;
          dd = 1.-d;
	// tapering function - mdf style
          tt = (1. + 3.*d + 6.*d*d)*dd*dd*dd;
	// minus the derivative of the tapering function
          dt = 30.* d*d * dd*dd * rr / dp;

          forcelj = forcelj*tt + philj*dt;
        } else {
          tt = 1;
        }
        fpair = factor_lj*forcelj*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) { // begin of eflag calculation
           if((i_domain_id != j_domain_id) || (i_domain_id == 0) || (j_domain_id == 0)){
		evdwl = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);  // SOP lj potential        
          } else{
          evdwl = 0.0;
          }
        
          
          if (rsq > cut_inner_sq[itype][jtype]) evdwl *= tt;

          evdwl *= factor_lj;

          if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
        } // end of eflag claculation
      } // end of loop over neighbour list
    } // end of loop over j atoms
  } // end of loop over i atoms

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJ_SOP::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_inner,n+1,n+1,"pair:cut_inner");
  memory->create(cut_inner_sq,n+1,n+1,"pair:cut_inner_sq");
  memory->create(epsilon_lj,n+1,n+1,"pair:epsilon_lj");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
 
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJ_SOP::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  cut_inner_global = utils::numeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);
  if (cut_inner_global <= 0.0 || cut_inner_global > cut_global)
    error->all(FLERR,"Illegal pair_style command");

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJ_SOP::coeff(int narg, char **arg)
{
  if (narg != 4 && narg != 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_lj_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);


  
  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;
  if (narg == 6) {
    cut_inner_one = utils::numeric(FLERR,arg[6],false,lmp);
    cut_one = utils::numeric(FLERR,arg[7],false,lmp);
  }
  if (cut_inner_global <= 0.0 || cut_inner_global > cut_global)
    error->all(FLERR,"Illegal pair_style command");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_lj[i][j] = epsilon_lj_one;
      sigma[i][j] = sigma_one;
//      epsilon_rep[i][j] = epsilon_rep_one;

      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJ_SOP::init_one(int i, int j)
{
//  if (setflag[i][j] == 0) {
//    epsilon_lj[i][j] = mix_energy(epsilon_lj[i][i],epsilon_lj[j][j],
//                               sigma[i][i],sigma[j][j]);
//    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
//    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
//    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
//  }

  cut_inner_sq[i][j] = cut_inner[i][j]*cut_inner[i][j];

  lj1[i][j] = 12.0 * epsilon_lj[i][j] * pow(sigma[i][j],12.0); // LJ SOP force related
  lj2[i][j] = 12.0 * epsilon_lj[i][j] * pow(sigma[i][j],6.0); // LJ SOP force related
  lj3[i][j] = epsilon_lj[i][j] * pow(sigma[i][j],12.0); // LJ SOP energy related
  lj4[i][j] = 2.0 * epsilon_lj[i][j] * pow(sigma[i][j],6.0); // LJ SOP energy related



  cut[j][i] = cut[i][j];  // BUG FIX

  cut_inner[j][i] = cut_inner[i][j];
  cut_inner_sq[j][i] = cut_inner_sq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];


  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJ_SOP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon_lj[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut_inner[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJ_SOP::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&epsilon_lj[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_inner[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&epsilon_lj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_inner[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJ_SOP::write_restart_settings(FILE *fp)
{
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJ_SOP::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJ_SOP::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon_lj[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJ_SOP::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",
              i,j,epsilon_lj[i][j],sigma[i][j],
              cut_inner[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJ_SOP::single(int i, int j, int itype, int jtype,
                             double rsq,
                             double /*factor_coul*/, double factor_lj,
                             double &fforce)
{
  double r2inv,r6inv,forcelj,philj;
  double rr, dp, d, tt, dt, dd;
  tagint i_residue_id, j_residue_id, i_domain_id, j_domain_id;
  i_residue_id = atom->residue_id[i];
  j_residue_id = atom->residue_id[j];
  i_domain_id = atom->domain_id[i]; 
  j_domain_id = atom->domain_id[j];
  
  
  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;

 if((i_domain_id != j_domain_id) || (i_domain_id == 0) || (j_domain_id == 0)){
         forcelj = r6inv*(lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
         philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);
  }else{
  	forcelj = 0.0;
  	philj = 0.0;
  }  
  
  
  if (rsq > cut_inner_sq[itype][jtype]) {

    rr = sqrt(rsq);
    dp = (cut[itype][jtype] - cut_inner[itype][jtype]);
    d = (rr - cut_inner[itype][jtype]) / dp;
    dd = 1-d;
    tt = (1. + 3.*d + 6.*d*d)* dd*dd*dd;
    dt = 30.* d*d * dd*dd * rr / dp;

    forcelj = forcelj*tt + philj*dt;
    philj *= tt;
  }

  fforce = factor_lj*forcelj*r2inv;

  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJ_SOP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon_lj") == 0) return (void *) epsilon_lj;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return nullptr;
}
