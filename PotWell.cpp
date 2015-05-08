#include <algorithm>
#include <iostream>

#include "PotWell.h"


PotWell::PotWell(std::string direction, double depth, double width){
  set_direction(direction);
  set_depth(depth);
  set_width(width);
}

void PotWell::set_direction(std::string direction){
  direction_ = direction;
  std::transform(direction_.begin(), direction_.end(), direction_.begin(), ::tolower);

  if (direction == "x"){
    num_direction_ = 1;
  } else if (direction == "y"){
    num_direction_ = 2;
  } else if (direction == "z"){
    num_direction_ = 3;
  } else{
    std::cerr << "Direction not recognized. Default direction (x) initialized." << std::endl;
    num_direction_ = 1;
    direction_ = "x";
  }
}

void PotWell::set_depth(double depth){
  depth_ = depth;
}

void PotWell::set_width(double width){
  width_ = width;
}

void PotWell::set_width_in_cm(double width_in_cm){
  width_in_cm_ = width_in_cm;
}

double PotWell::width_in_cm(){
  return width_in_cm_;
}

double PotWell::depth(){
  return depth_;
}

double PotWell::width(){
  return width_;
}

std::string PotWell::direction(){
  return direction_;
}

double PotWell::num_direction(){
  return num_direction_;
}


