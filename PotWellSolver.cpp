

#include "PotWellSolver.h"


PotWellSolver::PotWellSolver(Compound* compound, PotWell* potential_well,
							 int num_gridpoints, int x_min, int x_max,
							 int num_bandmodel)
    : compound_(compound),
	  potential_well_(potential_well),
	  num_bandmodel_(num_bandmodel){

  set_grid_params(num_gridpoints, x_min, x_max);
  set_units();
}


void PotWellSolver::set_grid_params(int num_gridpoints, int x_min, int x_max){
  num_gridpoints_ = num_gridpoints;
  x_max_ = x_max;
  x_min_ = x_min;
  x_axis_ = arma::linspace(x_min, x_max, num_gridpoints);
  step_size_ = ((double)x_max-(double)x_min)/x_axis_.n_elem;
  left_potential_well_boundary_ = floor(x_axis_.n_elem/2) - floor(potential_well_->width()/(2*step_size_));
  right_potential_well_boundary_ = floor(x_axis_.n_elem/2) + ceil(potential_well_->width()/(2*step_size_));
}

void PotWellSolver::set_units(){
  double hbar = 6.58211928*exp10(-13); // meV.s
  double mass_electron = 5.6778*exp10(-13); // meV.s^2/cm^2
  unit_energy_ = (pow(hbar,2)*pow(pi,2))/(2*mass_electron*pow(potential_well_->width_in_cm(),2));
  unit_potential_ = potential_well_->depth()/unit_energy_;
  unit_delta_ = compound_->delta()/unit_energy_;
}

arma::vec PotWellSolver::x_axis(){
  return x_axis_;
}

double PotWellSolver::step_size(){
  return step_size_;
}

int PotWellSolver::left_potential_well_boundary(){
  return left_potential_well_boundary_;
}

int PotWellSolver::right_potential_well_boundary(){
  return right_potential_well_boundary_;
}

int PotWellSolver::x_max(){
  return x_max_;
}

int PotWellSolver::x_min(){
  return x_min_;
}

int PotWellSolver::num_bandmodel(){
  return num_bandmodel_;
}

int PotWellSolver::num_gridpoints(){
  return num_gridpoints_;
}

double PotWellSolver::unit_delta(){
  return unit_delta_;
}

double PotWellSolver::unit_energy(){
  return unit_energy_;
}

double PotWellSolver::unit_potential(){
  return unit_potential_;
}

