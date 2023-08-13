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

#include "pair_rep_sop.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include <iostream>


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairREP_SOP::PairREP_SOP(LAMMPS *lmp) : Pair(lmp) {
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairREP_SOP::~PairREP_SOP()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon_rep); // excluded volume energy scale 
    memory->destroy(sigma); // Sigma between given bead pairs
    
    memory->destroy(rep1); // excluded volume force factor 1
    memory->destroy(rep2); // excluded volume force factor 2


  }
}

/* ---------------------------------------------------------------------- */

void PairREP_SOP::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,ev_rep,fpair;
  double rsq,r2inv,r6inv,force_rep;
  int *ilist;

  ev_rep = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special; //  Array containing local neighbours
  tagint *tag = atom->tag;  // Mapping local id to  global id of the atoms
  int my_i;

  // loop over neighbors of my atoms

  for (i = 0; i < nlocal; i++) {

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

       for (my_i= 0; my_i < nspecial[i][0]; my_i++){
        j = atom->map(special[i][my_i]);
        if(i<j){
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;       
        jtype = type[j];
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        force_rep =  r6inv*r2inv*rep1[itype][jtype];
//     std::cout << r2inv<< "\n";
//             std::cout << i << "\t"<< j <<" "<< force_rep<<"\n";
        f[i][0] += delx*force_rep;
        f[i][1] += dely*force_rep;
        f[i][2] += delz*force_rep;
//        if (newton_pair || j < nlocal) {
//       std::cout << i <<" "<< j <<"\n";
          f[j][0] -= delx*force_rep;
          f[j][1] -= dely*force_rep;
          f[j][2] -= delz*force_rep;
  //      }  
        if(eflag) {
        ev_rep = r6inv*rep2[itype][jtype];
        }
        if(evflag) ev_tally(i,j,nlocal,newton_pair,ev_rep,0.0,fpair,delx,dely,delz);      
        
        }
        } // end of 1-2 exv interactions

       for (my_i= nspecial[i][0]; my_i < nspecial[i][1]; my_i++){
        j = atom->map(special[i][my_i]);
        if(i<j){
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;       
        jtype = type[j];
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        force_rep =  r6inv*r2inv*rep1[itype][jtype];

        
        f[i][0] += delx*force_rep;
        f[i][1] += dely*force_rep;
        f[i][2] += delz*force_rep;
    //    if (newton_pair || j < nlocal) {
          f[j][0] -= delx*force_rep;
          f[j][1] -= dely*force_rep;
          f[j][2] -= delz*force_rep;
      //  }  
        if(eflag) {
        ev_rep = r6inv*rep2[itype][jtype];
        }
        if(evflag) ev_tally(i,j,nlocal,newton_pair,ev_rep, 0.0,fpair,delx,dely,delz);      
        
        }
        } // end of 1-2 exv interactions


  } // end of loop over i atoms

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairREP_SOP::allocate()
{        
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon_rep,n+1,n+1,"pair:epsilon_rep");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(rep1,n+1,n+1,"pair:rep1");
  memory->create(rep2,n+1,n+1,"pair:rep2");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairREP_SOP::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
  if (cut_global <= 0.0)
    error->all(FLERR,"Illegal pair_style command");

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairREP_SOP::coeff(int narg, char **arg)
{
  if (narg != 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_rep_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);
  double cut_one = cut_global;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_rep[i][j] = epsilon_rep_one;
      sigma[i][j] = sigma_one;

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

double PairREP_SOP::init_one(int i, int j)
{

//  if (setflag[i][j] == 0) {
//    epsilon_lj[i][j] = mix_energy(epsilon_lj[i][i],epsilon_lj[j][j],
//                               sigma[i][i],sigma[j][j]);
//    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
//    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
//    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
//  }

//  cut_inner_sq[i][j] = cut_inner[i][j]*cut_inner[i][j];

  rep1[i][j] = 6.0*epsilon_rep[i][j] * pow(sigma[i][j],6.0); // excluded volume force related
  rep2[i][j] = epsilon_rep[i][j] * pow(sigma[i][j],6.0); // excluded volume energy related

  cut[j][i] = cut[i][j];  // BUG FIX
  rep1[j][i] = rep1[i][j];
  rep2[j][i] = rep2[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairREP_SOP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon_rep[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairREP_SOP::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&epsilon_rep[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&epsilon_rep[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairREP_SOP::write_restart_settings(FILE *fp)
{
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairREP_SOP::read_restart_settings(FILE *fp)
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

void PairREP_SOP::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon_rep[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairREP_SOP::write_data_all(FILE *fp)
{ 
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",
              i,j,epsilon_rep[i][j],sigma[i][j],cut[i][j]);
              
}

/* ---------------------------------------------------------------------- */

double PairREP_SOP::single(int i, int j, int itype, int jtype,
                             double rsq,
                             double /*factor_coul*/, double factor_lj,
                             double &fforce)
{
  double r2inv,r6inv,force_rep,ev_rep;
  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  force_rep =  r6inv*rep1[itype][jtype];
  ev_rep = r6inv*rep2[itype][jtype];
  fforce = force_rep;

  return ev_rep;
}

/* ---------------------------------------------------------------------- */

void *PairREP_SOP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon_rep") == 0) return (void *) epsilon_rep;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return nullptr;
}
