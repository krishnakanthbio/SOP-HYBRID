/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Modified by : Krishnakanth B, Theoretical Biophysics Laboratory, 
   Molecular Biophysics Unit, Indian Institute of Science, Bangalore - 560012
------------------------------------------------------------------------- */


#include "pair_debye_sop.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include <iostream>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairDebye_SOP::PairDebye_SOP(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairDebye_SOP::~PairDebye_SOP()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    memory->destroy(epsilon_coul);
    memory->destroy(kappa);
  }
}

/* ---------------------------------------------------------------------- */

void PairDebye_SOP::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul, fpair, r, rinv, screening;
  double rsq, r2inv, r6inv, forcecoul,  factor_coul, dielectric_prefactor;
  int *ilist, *jlist, *numneigh, **firstneigh;

  ecoul = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqr2e = force->qqr2e;  
  tagint *domain_id = atom->domain_id;
  tagint i_domain_id, j_domain_id;
  tagint *tag = atom->tag;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    i_domain_id = domain_id[i];
    
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j_domain_id = domain_id[j];

      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      

     if( (i_domain_id == 0) || (j_domain_id == 0) || (i_domain_id != j_domain_id) ){
         dielectric_prefactor = 1.0/dielectric_uf;
        }else{
//        dielectric_prefactor = 1.0/dielectric_f;
        dielectric_prefactor = 0.0;//1.0/dielectric_f;
//        std::cout<<qqr2e<<"\n";
        }      

      if (rsq < cutsq[itype][jtype]){
        r2inv = 1.0 / rsq;

        if (rsq < cut_coulsq[itype][jtype]){
          r = sqrt(rsq);
          rinv = 1.0/r;
          screening = exp(-kappa[itype][jtype]*r);
          forcecoul = qqr2e*dielectric_prefactor *qtmp * q[j] * screening * (kappa[itype][jtype] + rinv);
          }
        else{
          forcecoul = 0.0;
          }

        fpair = (epsilon_coul[itype][jtype] *factor_coul *forcecoul ) * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }
//       if(abs(fpair) > 1.0){
//     std::cout  << tag[i] << " " <<  tag[j] <<" " <<fpair << " "<< dielectric_prefactor<<"\n";
//   }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype]){
            ecoul = factor_coul*qqr2e*dielectric_prefactor*epsilon_coul[itype][jtype]  * qtmp * q[j] * rinv * screening;
            }
          else{
            ecoul = 0.0;
            }

        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, 0.0, ecoul, fpair, delx, dely, delz);

      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDebye_SOP::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(cut_coul, np1, np1, "pair:cut_coul");
  memory->create(cut_coulsq, np1, np1, "pair:cut_coulsq");
  memory->create(epsilon_coul, np1, np1, "pair:epsilon_coul");
  memory->create(kappa, np1, np1, "pair:kappa");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDebye_SOP::settings(int narg, char **arg)
{
  if (narg < 2 || narg > 4) error->all(FLERR, "Illegal pair_style command");
  
    dielectric_uf = utils::numeric(FLERR, arg[0], false, lmp);
  dielectric_f = utils::numeric(FLERR, arg[1], false, lmp);  
  if (narg >= 3){
    kappa_global = utils::numeric(FLERR, arg[2], false, lmp);
  }
  if (narg == 4){
    cut_coul_global = utils::numeric(FLERR, arg[3], false, lmp);  
    }

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          kappa[i][j] = kappa_global;
          cut_coul[i][j] = cut_coul_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDebye_SOP::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 6) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_coul_one = utils::numeric(FLERR, arg[2], false, lmp);
  double kappa_one = kappa_global;
  double cut_coul_one = cut_coul_global;
  if (narg >= 3) kappa_one = utils::numeric(FLERR, arg[3], false, lmp);
  if (narg == 4) cut_coul_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      kappa[i][j] = kappa_one;
      epsilon_coul[i][j] = epsilon_coul_one;
      cut_coul[i][j] = cut_coul_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDebye_SOP::init_style()
{
  if (!atom->q_flag) error->all(FLERR, "Pair style debye/sop requires atom attribute q");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDebye_SOP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    error->all(FLERR, "Pair coefficients not set by the user");
  }

  double cut =  cut_coul[i][j];

  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDebye_SOP::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon_coul[i][j], sizeof(double), 1, fp);
        fwrite(&kappa[i][j], sizeof(double), 1, fp);
        fwrite(&cut_coul[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDebye_SOP::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &epsilon_coul[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &kappa[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut_coul[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon_coul[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&kappa[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_coul[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDebye_SOP::write_restart_settings(FILE *fp)
{
  fwrite(&dielectric_uf, sizeof(double), 1, fp);
  fwrite(&dielectric_f, sizeof(double), 1, fp);
  fwrite(&cut_coul_global, sizeof(double), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDebye_SOP::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &dielectric_uf, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &dielectric_f, sizeof(double), 1, fp, nullptr, error);    
    utils::sfread(FLERR, &cut_coul_global, sizeof(double), 1, fp, nullptr, error);
  }
  MPI_Bcast(&dielectric_uf, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dielectric_f, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_coul_global, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairDebye_SOP::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g\n", i, epsilon_coul[i][i], kappa[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairDebye_SOP::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j, epsilon_coul[i][j], kappa[i][j], cut_coul[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairDebye_SOP::single(int i, int j, int itype, int jtype, double rsq, double factor_coul,
                                double /*factor_lj*/, double &fforce)
{
  double r2inv, r6inv, forcecoul, phicoul, screening, rinv, r;
  double dielectric_prefactor;
  tagint i_domain_id, j_domain_id;
  i_domain_id = atom->domain_id[i];
  j_domain_id = atom->domain_id[j];  
  
     if( (i_domain_id == 0) || (j_domain_id == 0) || (i_domain_id != j_domain_id)){
         dielectric_prefactor = 1.0/dielectric_uf;
        }else{
//        dielectric_prefactor = 1.0/dielectric_f;        
        dielectric_prefactor = 0.0;//1.0/dielectric_f;
//        std::cout<<qqr2e<<"\n";
        }       

  r2inv = 1.0 / rsq;
  if (rsq < cut_coulsq[itype][jtype]){
	r = sqrt(rsq);
	rinv = 1.0/r;
	screening = exp(-kappa[itype][jtype]*r);
	forcecoul = force->qqr2e *dielectric_prefactor * epsilon_coul[itype][jtype]*atom->q[i]*atom->q[j]*screening*(kappa[itype][jtype] + rinv);
  }else{
  forcecoul = 0.0;
  }
  fforce = factor_coul*forcecoul * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul =  force->qqr2e *epsilon_coul[itype][jtype]* dielectric_prefactor * atom->q[i]*atom->q[j] * rinv * screening;
    eng += factor_coul * phicoul;
  }


  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairDebye_SOP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon_coul") == 0) return (void *) epsilon_coul;
  if (strcmp(str, "kappa") == 0) return (void *) kappa;
  if (strcmp(str, "cut_coul") == 0) return (void *) cut_coul;
  return nullptr;
}
