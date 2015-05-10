#ifndef MATRIXSOLVER_H
#define MATRIXSOLVER_H

#include "HamiltonianMaker.h"


class MatrixSolver{
 public:
  MatrixSolver(HamiltonianMaker* hamiltonian_maker);
  arma::vec heavy_hole_mixing(arma::vec k_vector);
  arma::vec light_hole_mixing(arma::vec k_vector);
  arma::vec rotated_heavy_hole_mixing(arma::vec k_vector);
  arma::vec rotated_light_hole_mixing(arma::vec k_vector);

  HamiltonianMaker* hamiltonian_maker();


 private:
  HamiltonianMaker* hamiltonian_maker_;
  arma::vec mixing_calculator(std::string mixing_result, arma::vec k_vector);


};

#endif // MATRIXSOLVER_H
