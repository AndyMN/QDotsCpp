#ifndef POTWELLSOLVER_H
#define POTWELLSOLVER_H

#include <armadillo>

#include "Compound.h"
#include "PotWell.h"

class PotWellSolver{
 public:
  PotWellSolver(Compound* compound, PotWell* potential_well,
                int num_gridpoints = 100, int x_min = -3,
                int x_max = 3, int num_bandmodel = 4);

  void set_grid_params(int num_gridpoints, int x_min = -3, int x_max = 3);

  Compound* compound();
  PotWell* potential_well();

  int num_gridpoints();
  int x_max();
  int x_min();
  arma::vec x_axis();
  double step_size();
  int left_potential_well_boundary();
  int right_potential_well_boundary();

  double unit_energy();
  double unit_potential();
  double unit_delta();
  int num_bandmodel();

 private:
  Compound* compound_;
  PotWell* potential_well_;

  int num_gridpoints_;
  int x_max_;
  int x_min_;
  arma::vec x_axis_;
  double step_size_;
  int left_potential_well_boundary_;
  int right_potential_well_boundary_;

  void set_units();
  double unit_energy_;
  double unit_potential_;
  double unit_delta_;
  int num_bandmodel_ = 4;

  static constexpr double pi = 3.14159265359;






};

#endif
