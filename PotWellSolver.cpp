

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

}
