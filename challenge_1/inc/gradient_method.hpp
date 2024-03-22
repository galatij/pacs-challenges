#ifndef _GRADIENT_METHOD_
#define _GRADIENT_METHOD_

#include "parameters.hpp"
#include <vector>
#include <functional>


double norm2_square(const std::vector<double> &);

double compute_minimum(
    const std::function<double(const std::vector<double> &)> &,
    const std::function<std::vector<double>(const std::vector<double> &)> &,
    const params &,
    bool verbose);

double compute_minimum_gradient(
    const std::function<double(const std::vector<double> &)> &,
    const std::function<std::vector<double>(const std::vector<double> &)> &,
    const params &,
    bool verbose);

double compute_minimum_momentum_1(
    const std::function<double(const std::vector<double> &)> &,
    const std::function<std::vector<double>(const std::vector<double> &)> &,
    const params &,
    bool verbose);

double compute_minimum_momentum_2(
    const std::function<double(const std::vector<double> &)> &,
    const std::function<std::vector<double>(const std::vector<double> &)> &,
    const params &,
    bool verbose);

double compute_minimum_momentum_nesterov(
    const std::function<double(const std::vector<double> &)> &,
    const std::function<std::vector<double>(const std::vector<double> &)> &,
    const params &,
    bool verbose);


template<step_type T>
double compute_step(double a0, double mu, unsigned int k, double sigma,
                    std::vector<double> x,
                    const std::function<double(std::vector<double>)> &f,
                    const std::function<std::vector<double>(std::vector<double>)> &df);

bool converged(const std::vector<double> &x1,
               const std::vector<double> &x0,
               const std::function<double (std::vector<double>)> &f, 
               const double tol_res,
               const double tol_x,
               bool verbose
               );


//------------------
void print_vector(const std::vector<double> &vec);

#endif // _GRADIENT_METHOD_