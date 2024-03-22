#include "parameters.hpp"

std::ostream & operator<<(std::ostream &os, const params & p){
    os << "PARAMETERS:\n";
    os << "max_iter: " << p.max_it << "\n";
    os << "tol_res: " << p.tol_res << "\n";
    os << "tol_x: " << p.tol_x << "\n";
    os << "x0: [ ";
    for (const double &v: p.x0)
        os << v << " ";
    os << "]\n";
    os << "a0: " << p.a0 << "\n";
    os << "mu: " << p.mu << "\n";

    if (p.method == 0)
    {
        os << "step type: ";
        switch(p.step){
            case (step_type::Expo):
                os << "Exponential decay" << std::endl;
                break;
            case (step_type::Inverse):
                os << "Inverse decay" << std::endl;
                break;
            case (step_type::Line):
                os << "Approximate Line Search (sigma = " <<
                p.sigma << ")" << std::endl;
                break;
        }
    }
    else
    {
        os << "variant: ";
        switch(p.var){
            case (heavy_ball_variant::variant_1):
                os << "First variant" << std::endl;
                break;
            case (heavy_ball_variant::variant_2):
                os << "Second variant" << std::endl;
                break;
            case (heavy_ball_variant::variant_Nesterov):
                os << "Nesterov" << std::endl;
                break;
        }
    }
    return os;
}