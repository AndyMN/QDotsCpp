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

  double depth();
  double width();
  double num_direction();
  std::string direction();

 private:
  double depth_;
  double width_;
  double num_direction_;
  std::string direction_;


};

#endif
