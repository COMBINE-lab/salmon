#ifndef SIMPLE_POS_BIAS_HPP
#define SIMPLE_POS_BIAS_HPP

class SimplePosBias {
public:
  SimplePosBias(int32_t numBins);
  void addMass(int32_t bin, double mass);
  void addMass(int32_t pos, int32_t length, double mass);
private:
  std::vector<double> masses_;
  std::vector<double> cdf_;
  bool isLogged_{true};
  bool haveCDF_{false};
};

#endif // SIMPLE_POS_BIAS_HPP
