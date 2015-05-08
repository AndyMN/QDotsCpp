#ifndef COMPOUND_H
#define COMPOUND_H

#include <vector>
#include <string>

extern const std::vector<double> gaas_params;
extern const std::vector<double> si_params;

class Compound{
 public:
  Compound(std::vector<double> params, std::string compound_name = "\0");
  std::string compound_name();
  std::vector<double> params();
  double y1();
  double y2();
  double y3();
  double delta();

 private:
  void set_params(std::vector<double> params);
  void set_compound_name(std::string compound_name);

  std::vector<double> params_;
  double y1_;
  double y2_;
  double y3_;
  double delta_;
  std::string compound_name_;

};


#endif // COMPOUND_H
