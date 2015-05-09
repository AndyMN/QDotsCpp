#include "MatrixSolver.h"

MatrixSolver::MatrixSolver(HamiltonianMaker* hamiltonian_maker) : hamiltonian_maker_(hamiltonian_maker){
}

arma::vec MatrixSolver::heavy_hole_mixing(arma::vec k_vector){
  int num_elements = k_vector.n_elem;
  int num_gridpoints = hamiltonian_maker_->num_gridpoints();

  arma::vec heavy_hole_fractions(num_elements, arma::fill::zeros);

  for (int i = 0; i < num_elements; i++) {
    arma::cx_mat hamiltonian = hamiltonian_maker_->hamiltonian_luttinger_kohn(k_vector[i]);
    arma::vec eigvals;
    arma::cx_mat eigvecs;
    arma::eig_sym(eigvals, eigvecs, hamiltonian);

    arma::cx_vec ground_state_eigenvector = eigvecs.col(0);

    arma::cx_vec heavy_hole1 = ground_state_eigenvector.subvec(0,num_gridpoints-1);
    arma::cx_vec light_hole1 = ground_state_eigenvector.subvec(num_gridpoints, 2*num_gridpoints-1);
    arma::cx_vec light_hole2 = ground_state_eigenvector.subvec(2*num_gridpoints, 3*num_gridpoints-1);
    arma::cx_vec heavy_hole2 = ground_state_eigenvector.subvec(3*num_gridpoints, 4*num_gridpoints-1);

    double normsq_heavy_hole1 = pow(arma::norm(heavy_hole1),2);
    double normsq_light_hole1 = pow(arma::norm(light_hole1),2);
    double normsq_light_hole2 = pow(arma::norm(light_hole2),2);
    double normsq_heavy_hole2 = pow(arma::norm(heavy_hole2),2);

    double normsq_total_wavefunction = normsq_heavy_hole1 + normsq_heavy_hole2 + normsq_light_hole1 + normsq_light_hole2;

    heavy_hole_fractions[i] = (normsq_heavy_hole1 + normsq_heavy_hole2)/normsq_total_wavefunction;
  }
  return heavy_hole_fractions;
}

arma::vec MatrixSolver::light_hole_mixing(arma::vec k_vector){
  int num_elements = k_vector.n_elem;
  int num_gridpoints = hamiltonian_maker_->num_gridpoints();

  arma::vec light_hole_fractions(num_elements, arma::fill::zeros);

  for (int i = 0; i < num_elements; i++) {
    arma::cx_mat hamiltonian = hamiltonian_maker_->hamiltonian_luttinger_kohn(k_vector[i]);
    arma::vec eigvals;
    arma::cx_mat eigvecs;
    arma::eig_sym(eigvals, eigvecs, hamiltonian);

    arma::cx_vec ground_state_eigenvector = eigvecs.col(0);

    arma::cx_vec heavy_hole1 = ground_state_eigenvector.subvec(0,num_gridpoints-1);
    arma::cx_vec light_hole1 = ground_state_eigenvector.subvec(num_gridpoints, 2*num_gridpoints-1);
    arma::cx_vec light_hole2 = ground_state_eigenvector.subvec(2*num_gridpoints, 3*num_gridpoints-1);
    arma::cx_vec heavy_hole2 = ground_state_eigenvector.subvec(3*num_gridpoints, 4*num_gridpoints-1);

    double normsq_heavy_hole1 = pow(arma::norm(heavy_hole1),2);
    double normsq_light_hole1 = pow(arma::norm(light_hole1),2);
    double normsq_light_hole2 = pow(arma::norm(light_hole2),2);
    double normsq_heavy_hole2 = pow(arma::norm(heavy_hole2),2);

    double normsq_total_wavefunction = normsq_heavy_hole1 + normsq_heavy_hole2 + normsq_light_hole1 + normsq_light_hole2;

    light_hole_fractions[i] = (normsq_light_hole1 + normsq_light_hole2)/normsq_total_wavefunction;
  }
  return light_hole_fractions;
}

