#ifndef _PARAMETERS_
#define _PARAMETERS_

#include <vector>
#include <iostream>

enum step_type{
    Expo,
    Inverse,
    Line
};

enum heavy_ball_variant{
    variant_1,
    variant_2,
    variant_Nesterov
};

struct params{
    unsigned method = 0;            /* 0 = gradient method (default)
                                       1 = hevay_ball (variant_1)
                                       2 = hevay_ball (variant_2)
                                       3 = hevay_ball (variant_Nesterov)
                                    */

    unsigned int max_it = 2000;
    double tol_res = 1e-08;
    double tol_x = 1e-08;
    std::vector<double> x0{};
    double a0 = 0.2;
    double mu = 0.2;
    double sigma = 0.2;
    step_type step = Expo; 
    heavy_ball_variant var = variant_1;
};

std::ostream &operator<<(std::ostream &, const params &);

#endif // _PARAMETERS_