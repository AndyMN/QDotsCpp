#include "HamiltonianMaker.h"


HamiltonianMaker::HamiltonianMaker(Compound* compound, PotWell* potential_well,
							 int num_gridpoints, int x_min, int x_max,
							 int num_bandmodel)
    : compound_(compound),
	  potential_well_(potential_well),
	  num_bandmodel_(num_bandmodel){

  set_grid_params(num_gridpoints, x_min, x_max);
  set_units();
}


void HamiltonianMaker::set_grid_params(int num_gridpoints, int x_min, int x_max){
  num_gridpoints_ = num_gridpoints;
  x_max_ = x_max;
  x_min_ = x_min;
  x_axis_ = arma::linspace(x_min, x_max, num_gridpoints);
  step_size_ = ((double)x_max-(double)x_min)/x_axis_.n_elem;
  left_potential_well_boundary_ = floor(x_axis_.n_elem/2) - floor(potential_well_->width()/(2*step_size_));
  right_potential_well_boundary_ = floor(x_axis_.n_elem/2) + ceil(potential_well_->width()/(2*step_size_));
}

void HamiltonianMaker::set_units(){
  double hbar = 6.58211928*exp10(-13); // meV.s
  double mass_electron = 5.6778*exp10(-13); // meV.s^2/cm^2
  unit_energy_ = (pow(hbar,2)*pow(pi,2))/(2*mass_electron*pow(potential_well_->width_in_cm(),2));
  unit_potential_ = potential_well_->depth()/unit_energy_;
  unit_delta_ = compound_->delta()/unit_energy_;
}

arma::cx_mat HamiltonianMaker::hamiltonian_luttinger_kohn(double k){
  double y1 = compound_->y1();
  double y2 = compound_->y2();
  double y3 = compound_->y3();

  arma::cx_mat P = arma::zeros<arma::cx_mat>(num_gridpoints_, num_gridpoints_);
  arma::cx_mat Q = arma::zeros<arma::cx_mat>(num_gridpoints_, num_gridpoints_);
  arma::cx_mat R = arma::zeros<arma::cx_mat>(num_gridpoints_, num_gridpoints_);
  arma::cx_mat S = arma::zeros<arma::cx_mat>(num_gridpoints_, num_gridpoints_);
  arma::cx_mat V = arma::zeros<arma::cx_mat>(num_gridpoints_, num_gridpoints_);

  arma::cx_double diagonal_P_element;
  arma::cx_double subdiagonal_P_element;
  arma::cx_double superdiagonal_P_element;

  arma::cx_double diagonal_Q_element;
  arma::cx_double subdiagonal_Q_element;
  arma::cx_double superdiagonal_Q_element;

  arma::cx_double diagonal_R_element;
  arma::cx_double subdiagonal_R_element;
  arma::cx_double superdiagonal_R_element;

  arma::cx_double diagonal_S_element;
  arma::cx_double subdiagonal_S_element;
  arma::cx_double superdiagonal_S_element;


  if (potential_well_->num_direction() == 1) {
    diagonal_P_element.real((y1/pow(pi,2)) * (pow(k,2) + 2.0/pow(step_size_,2)));
    subdiagonal_P_element.real((-y1/(pow(step_size_,2)) * pow(pi,2)));
    superdiagonal_P_element = subdiagonal_P_element;

    diagonal_Q_element.real((y2/pow(pi,2)) * (-2*pow(k,2) + 2.0/pow(step_size_,2)));
    subdiagonal_Q_element.real(-y2/(pow(step_size_,2) * pow(pi,2)));
    superdiagonal_Q_element = subdiagonal_Q_element;

    diagonal_R_element.real((sqrt(3)*y2/pow(pi,2)) * (-2.0/pow(step_size_,2)));
    subdiagonal_R_element.real((sqrt(3)/pow(pi,2)) * (y2/pow(step_size_,2)));
    superdiagonal_R_element.real((sqrt(3)/pow(pi,2)) * (y2/pow(step_size_,2)));

	subdiagonal_S_element.imag(sqrt(3)*y3*k/(pow(pi,2)*step_size_));
    superdiagonal_S_element = -subdiagonal_S_element;
  } else if (potential_well_->num_direction() == 3){
    diagonal_P_element.real((y1/pow(pi,2)) * (pow(k,2) + 2.0/pow(step_size_,2)));
    subdiagonal_P_element.real((-y1/(pow(step_size_,2)*pow(pi,2))));
    superdiagonal_P_element = subdiagonal_P_element;

    diagonal_Q_element.real((y2/pow(pi,2)) * (pow(k,2) - 4.0/pow(step_size_,2)));
    subdiagonal_Q_element.real((2*y2/(pow(step_size_,2)*pow(pi,2))));
    superdiagonal_Q_element = subdiagonal_Q_element;

    diagonal_R_element.real((1.0/pow(pi,2)) * (-sqrt(3)*y2*pow(k,2)));

    subdiagonal_S_element.imag(y3*sqrt(3)*k/(pow(pi,2)*step_size_));
    superdiagonal_S_element = -subdiagonal_S_element;
  }

  arma::cx_vec diagonal_P(num_gridpoints_);
  arma::cx_vec subdiagonal_P(num_gridpoints_-1);
  arma::cx_vec superdiagonal_P(num_gridpoints_-1);

  diagonal_P.fill(diagonal_P_element);
  subdiagonal_P.fill(subdiagonal_P_element);
  superdiagonal_P.fill(superdiagonal_P_element);

  arma::cx_vec diagonal_Q(num_gridpoints_);
  arma::cx_vec subdiagonal_Q(num_gridpoints_-1);
  arma::cx_vec superdiagonal_Q(num_gridpoints_-1);

  diagonal_Q.fill(diagonal_Q_element);
  subdiagonal_Q.fill(subdiagonal_Q_element);
  superdiagonal_Q.fill(superdiagonal_Q_element);

  arma::cx_vec diagonal_R(num_gridpoints_);
  arma::cx_vec subdiagonal_R(num_gridpoints_-1);
  arma::cx_vec superdiagonal_R(num_gridpoints_-1);

  diagonal_R.fill(diagonal_R_element);
  subdiagonal_R.fill(subdiagonal_R_element);
  superdiagonal_R.fill(superdiagonal_R_element);

  arma::cx_vec diagonal_S(num_gridpoints_);
  arma::cx_vec subdiagonal_S(num_gridpoints_-1);
  arma::cx_vec superdiagonal_S(num_gridpoints_-1);

  diagonal_S.fill(diagonal_S_element);
  subdiagonal_S.fill(subdiagonal_S_element);
  superdiagonal_S.fill(superdiagonal_S_element);

  arma::cx_vec potential(num_gridpoints_);
  arma::cx_double complex_unit_potential(unit_potential_,0);
  potential.fill(complex_unit_potential);
  V = diagmat(potential);

  P.diag() = diagonal_P;
  P.diag(-1) = subdiagonal_P;
  P.diag(1) = superdiagonal_P;

  Q.diag() = diagonal_Q;
  Q.diag(-1) = subdiagonal_Q;
  Q.diag(1) = superdiagonal_Q;

  R.diag() = diagonal_R;
  R.diag(-1) = subdiagonal_R;
  R.diag(1) = superdiagonal_R;

  S.diag() = diagonal_S;
  S.diag(-1) = subdiagonal_S;
  S.diag(1) = superdiagonal_S;

  arma::cx_mat luttinger_kohn_hamiltonian = arma::zeros<arma::cx_mat>(num_gridpoints_*num_bandmodel_,num_gridpoints_*num_bandmodel_);

  // FIRST ROW
  luttinger_kohn_hamiltonian.submat(0, 0, num_gridpoints_-1, num_gridpoints_-1) = P+Q+V;
  luttinger_kohn_hamiltonian.submat(0, num_gridpoints_, num_gridpoints_-1, 2*num_gridpoints_-1) = -1*S;
  luttinger_kohn_hamiltonian.submat(0, 2*num_gridpoints_, num_gridpoints_-1, 3*num_gridpoints_-1) = R;
  //luttinger_kohn_hamiltonian.submat(0, 3*num_gridpoints_, num_gridpoints_-1, 4*num_gridpoints_-1) = THIS IS A ZERO BLOCK

  // SECOND ROW
  luttinger_kohn_hamiltonian.submat(num_gridpoints_, 0, 2*num_gridpoints_-1, num_gridpoints_-1) = -1*S.t();
  luttinger_kohn_hamiltonian.submat(num_gridpoints_, num_gridpoints_, 2*num_gridpoints_-1, 2*num_gridpoints_-1) = P-Q+V;
  //luttinger_kohn_hamiltonian.submat(num_gridpoints_, 2*num_gridpoints_, 2*num_gridpoints_-1, 3*num_gridpoints_-1) = THIS IS A ZERO BLOCK
  luttinger_kohn_hamiltonian.submat(num_gridpoints_, 3*num_gridpoints_, 2*num_gridpoints_-1, 4*num_gridpoints_-1) = R;

  // THIRD ROW
  luttinger_kohn_hamiltonian.submat(2*num_gridpoints_, 0, 3*num_gridpoints_-1, num_gridpoints_-1) = R.t();
  //luttinger_kohn_hamiltonian.submat(2*num_gridpoints_, num_gridpoints_, 3*num_gridpoints_-1, 2*num_gridpoints_-1) = THIS IS A ZERO BLOCK
  luttinger_kohn_hamiltonian.submat(2*num_gridpoints_, 2*num_gridpoints_, 3*num_gridpoints_-1, 3*num_gridpoints_-1) = P-Q+V;
  luttinger_kohn_hamiltonian.submat(2*num_gridpoints_, 3*num_gridpoints_, 3*num_gridpoints_-1, 4*num_gridpoints_-1) = S;

  // FOURTH ROW
  //luttinger_kohn_hamiltonian.submat(3*num_gridpoints_, 0, 4*num_gridpoints_-1, num_gridpoints_-1) = THIS IS A ZERO BLOCK
  luttinger_kohn_hamiltonian.submat(3*num_gridpoints_, num_gridpoints_, 4*num_gridpoints_-1, 2*num_gridpoints_-1) = R.t();
  luttinger_kohn_hamiltonian.submat(3*num_gridpoints_, 2*num_gridpoints_, 4*num_gridpoints_-1, 3*num_gridpoints_-1) = S.t();
  luttinger_kohn_hamiltonian.submat(3*num_gridpoints_, 3*num_gridpoints_, 4*num_gridpoints_-1, 4*num_gridpoints_-1) = P+Q+V;


  return luttinger_kohn_hamiltonian;

}

arma::vec HamiltonianMaker::x_axis(){
  return x_axis_;
}

double HamiltonianMaker::step_size(){
  return step_size_;
}

int HamiltonianMaker::left_potential_well_boundary(){
  return left_potential_well_boundary_;
}

int HamiltonianMaker::right_potential_well_boundary(){
  return right_potential_well_boundary_;
}

int HamiltonianMaker::x_max(){
  return x_max_;
}

int HamiltonianMaker::x_min(){
  return x_min_;
}

int HamiltonianMaker::num_bandmodel(){
  return num_bandmodel_;
}

int HamiltonianMaker::num_gridpoints(){
  return num_gridpoints_;
}

double HamiltonianMaker::unit_delta(){
  return unit_delta_;
}

double HamiltonianMaker::unit_energy(){
  return unit_energy_;
}

double HamiltonianMaker::unit_potential(){
  return unit_potential_;
}

