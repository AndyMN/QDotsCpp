#include <iostream>

#include <armadillo>

#include "Compound.h"
#include "PotWell.h"
#include "PotWellSolver.h"

int main(int argc, char** argv)
  {
  arma::mat A = arma::randu<arma::mat>(4,5);
  arma::mat B = arma::randu<arma::mat>(4,5);

  Compound GaAs(gaas_params);
  PotWell potential_well_x("z");


  PotWellSolver potential_well_solver(&GaAs, &potential_well_x);


  return 0;
  }
