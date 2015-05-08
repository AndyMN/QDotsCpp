#include "Compound.h"


const std::vector<double> gaas_params = {6.85, 2.1, 2.9, 341};
const std::vector<double> si_params = {4.26, 0.38, 1.56, 44};

Compound::Compound(std::vector<double> params, std::string compound_name){
  set_params(params);
  set_compound_name(compound_name);
}

void Compound::set_compound_name(std::string compound_name){
  if (params_ == gaas_params){
    compound_name_ = "GaAs";
  } else if (params_ == si_params){
    compound_name_ = "Si";
  } else{
    compound_name_ = compound_name;
  }

}

void Compound::set_params(std::vector<double> params){
  params_ = params;
  y1_ = params[0];
  y2_ = params[1];
  y3_ = params[2];
  delta_ = params[3];
}

double Compound::delta(){
  return delta_;
}

double Compound::y1(){
  return y1_;
}

double Compound::y2(){
  return y2_;
}

double Compound::y3(){
  return y3_;
}

std::vector<double> Compound::params(){
  return params_;
}

std::string Compound::compound_name(){
  return compound_name_;
}


