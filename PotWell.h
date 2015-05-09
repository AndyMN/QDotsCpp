#ifndef POTWELL_H
#define POTWELL_H

#include <string>
#include "Compound.h"

class PotWell{
 public:
  PotWell(std::string direction = "x", double depth = 130, double width = 1);

  void set_depth(double depth);
  void set_width(double width);
  void set_direction(std::string direction);
  void set_width_in_cm(double width_in_cm);

  double depth();
  double width();
  double width_in_cm();
  int num_direction();
  std::string direction();

 private:
  double depth_;
  double width_;
  double width_in_cm_ = 100*exp10(-8);
  int num_direction_;
  std::string direction_;


};

#endif
