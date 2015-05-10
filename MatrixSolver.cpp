#include "MatrixSolver.h"
#include <cmath>
#include <fstream>

MatrixSolver::MatrixSolver(HamiltonianMaker* hamiltonian_maker) : hamiltonian_maker_(hamiltonian_maker){
}

arma::vec MatrixSolver::heavy_hole_mixing(arma::vec k_vector){
  return mixing_calculator("heavy_hole", k_vector);
}

arma::vec MatrixSolver::light_hole_mixing(arma::vec k_vector){
  return mixing_calculator("light_hole", k_vector);
}


arma::vec MatrixSolver::rotated_heavy_hole_mixing(arma::vec k_vector){
  return mixing_calculator("rot_heavy_hole", k_vector);
}

arma::vec MatrixSolver::rotated_light_hole_mixing(arma::vec k_vector){
  return mixing_calculator("rot_light_hole", k_vector);
}

arma::vec MatrixSolver::mixing_calculator(std::string mixing_result, arma::vec k_vector){
  bool light_hole_rotation = false;
  bool heavy_hole_rotation = false;
  bool light_hole = false;
  bool heavy_hole = false;

  if (mixing_result == "rot_light_hole"){
    light_hole_rotation = true;
  }else if (mixing_result == "rot_heavy_hole") {
    heavy_hole_rotation = true;
  } else if (mixing_result == "light_hole") {
    light_hole = true;
  } else if (mixing_result == "heavy_hole") {
    heavy_hole = true;
  } else {
    std::cerr << "Not recognized method" << std::endl;
  }

  int num_direction_potwell = hamiltonian_maker_->potential_well()->num_direction();
  int num_elements = k_vector.n_elem;
  int num_gridpoints = hamiltonian_maker_->num_gridpoints();
  int num_matrixdim = hamiltonian_maker_->num_bandmodel();

  arma::vec hole_fractions(num_elements, arma::fill::zeros);
  arma::vec eigvals(num_gridpoints, arma::fill::zeros);
  arma::cx_mat eigvecs = arma::zeros<arma::cx_mat>(num_gridpoints, num_gridpoints);

  arma::cx_vec ground_state_eigenvector(num_gridpoints*num_matrixdim, arma::fill::zeros);
  arma::cx_vec heavy_hole1(num_gridpoints, arma::fill::zeros);
  arma::cx_vec light_hole1(num_gridpoints, arma::fill::zeros);
  arma::cx_vec light_hole2(num_gridpoints, arma::fill::zeros);
  arma::cx_vec heavy_hole2(num_gridpoints, arma::fill::zeros);

  arma::cx_vec rotated_heavy_hole1;
  arma::cx_vec rotated_light_hole1;
  arma::cx_vec rotated_light_hole2;
  arma::cx_vec rotated_heavy_hole2;

  if (light_hole_rotation || heavy_hole_rotation) {
    arma::cx_vec rotated_heavy_hole1(num_gridpoints, arma::fill::zeros);
    arma::cx_vec rotated_light_hole1(num_gridpoints, arma::fill::zeros);
    arma::cx_vec rotated_light_hole2(num_gridpoints, arma::fill::zeros);
    arma::cx_vec rotated_heavy_hole2(num_gridpoints, arma::fill::zeros);
  }


  double normsq_heavy_hole1 = 0;
  double normsq_heavy_hole2 = 0;
  double normsq_light_hole1 = 0;
  double normsq_light_hole2 = 0;
  double normsq_total_wavefunction = 0;

  for (int i = 0; i < num_elements; i++) {
    arma::cx_mat hamiltonian = hamiltonian_maker_->hamiltonian_luttinger_kohn(k_vector[i]);
    arma::eig_sym(eigvals, eigvecs, hamiltonian);

    ground_state_eigenvector = eigvecs.col(0);

    heavy_hole1 = ground_state_eigenvector.subvec(0, num_gridpoints-1);
    light_hole1 = ground_state_eigenvector.subvec(num_gridpoints, 2*num_gridpoints-1);
    light_hole2 = ground_state_eigenvector.subvec(2*num_gridpoints, 3*num_gridpoints-1);
    heavy_hole2 = ground_state_eigenvector.subvec(3*num_gridpoints, 4*num_gridpoints-1);

	if (light_hole_rotation || heavy_hole_rotation) {
	  if (num_direction_potwell == 1){
        //POTWELL IN X DIRECTION SO ROTATION TO THE Z DIRECTION
        rotated_heavy_hole1 = sqrt(2)/4 * (heavy_hole1 + heavy_hole2) + sqrt(6)/4 * (light_hole1 + light_hole2);
        rotated_light_hole1 = sqrt(6)/4 * (heavy_hole2 - heavy_hole1) + sqrt(2)/4 * (light_hole2 - light_hole1);
        rotated_light_hole2 = sqrt(6)/4 * (heavy_hole1 + heavy_hole2) + sqrt(2)/4 * (-light_hole1 - light_hole2);
        rotated_heavy_hole2 = sqrt(2)/4 * (heavy_hole2 - heavy_hole1) + sqrt(6)/4 * (light_hole1 - light_hole2);
      }
	}
    if (light_hole_rotation || heavy_hole_rotation) {
      normsq_heavy_hole1 = pow(arma::norm(rotated_heavy_hole1), 2);
      normsq_light_hole1 = pow(arma::norm(rotated_light_hole1), 2);
      normsq_light_hole2 = pow(arma::norm(rotated_light_hole2), 2);
      normsq_heavy_hole2 = pow(arma::norm(rotated_heavy_hole2), 2);
	} else if (light_hole || heavy_hole) {
	  normsq_heavy_hole1 = pow(arma::norm(heavy_hole1), 2);
      normsq_light_hole1 = pow(arma::norm(light_hole1), 2);
      normsq_light_hole2 = pow(arma::norm(light_hole2), 2);
      normsq_heavy_hole2 = pow(arma::norm(heavy_hole2), 2);
	}

	normsq_total_wavefunction = normsq_heavy_hole1 + normsq_heavy_hole2 + normsq_light_hole1 + normsq_light_hole2;
    if (light_hole || light_hole_rotation) {
      hole_fractions[i] = (normsq_light_hole1 + normsq_light_hole2)/normsq_total_wavefunction;
    } else if (heavy_hole || heavy_hole_rotation) {
      hole_fractions[i] = (normsq_heavy_hole1 + normsq_heavy_hole2)/normsq_total_wavefunction;
    }
    std::cout << mixing_result << " : " << normsq_total_wavefunction << std::endl;
  }
  return hole_fractions;
}


HamiltonianMaker* MatrixSolver::hamiltonian_maker(){
  return hamiltonian_maker_;
}
