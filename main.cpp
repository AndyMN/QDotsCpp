#include <iostream>
#include <fstream>

#include <armadillo>
#include "gnuplot-iostream.h"

#include "Compound.h"
#include "PotWell.h"
#include "HamiltonianMaker.h"
#include "MatrixSolver.h"

int main(int argc, char** argv){
  Gnuplot gp;

  Compound GaAs(gaas_params);
  arma::vec k_vector = arma::linspace<arma::vec>(-3,3);
  int num_elem_kvec = k_vector.n_elem;

  PotWell potential_well_x("x");
  HamiltonianMaker hamiltonian_maker_x(&GaAs, &potential_well_x);
  MatrixSolver matrix_solver_x(&hamiltonian_maker_x);


  arma::vec heavy_hole_fractions_x = matrix_solver_x.heavy_hole_mixing(k_vector);
  arma::vec light_hole_fractions_x = matrix_solver_x.light_hole_mixing(k_vector);
  arma::vec rotated_heavy_hole_fractions_x = matrix_solver_x.rotated_heavy_hole_mixing(k_vector);


  arma::mat hh_x = arma::zeros<arma::mat>(num_elem_kvec, 2);
  arma::mat lh_x = arma::zeros<arma::mat>(num_elem_kvec, 2);
  arma::mat rot_hh_x = arma::zeros<arma::mat>(num_elem_kvec, 2);

  hh_x.col(0) = k_vector;
  hh_x.col(1) = heavy_hole_fractions_x;

  lh_x.col(0) = k_vector;
  lh_x.col(1) = light_hole_fractions_x;

  rot_hh_x.col(0) = k_vector;
  rot_hh_x.col(1) = rotated_heavy_hole_fractions_x;

  PotWell potential_well_z("z");
  HamiltonianMaker hamiltonian_maker_z(&GaAs, &potential_well_z);
  MatrixSolver matrix_solver_z(&hamiltonian_maker_z);

  arma::vec heavy_hole_fractions_z = matrix_solver_z.heavy_hole_mixing(k_vector);
  arma::vec light_hole_fractions_z = matrix_solver_z.light_hole_mixing(k_vector);

  arma::mat hh_z = arma::zeros<arma::mat>(num_elem_kvec, 2);
  arma::mat lh_z = arma::zeros<arma::mat>(num_elem_kvec, 2);

  hh_z.col(0) = k_vector;
  hh_z.col(1) = heavy_hole_fractions_z;

  lh_z.col(0) = k_vector;
  lh_z.col(1) = light_hole_fractions_z;



  gp << "set xrange [-3:3]\nset yrange [0:1]\n";
  gp << "plot ";
  gp << gp.binFile1d(hh_z, "record") << "with lines title 'heavy hole_z'";
  gp << ", ";
  gp << gp.binFile1d(rot_hh_x, "record") << "with lines title 'rot heavy hole_x'";
  gp << std::endl;




  return 0;
}
