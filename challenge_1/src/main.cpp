#include "read_parameters.hpp"
#include "gradient_method.hpp"
#include "GetPot"
#include <iostream>
double f(const std::vector<double> &);
std::vector<double> df(const std::vector<double> &x);

int main(int argc, char *argv[]){
    
    // read the input json file
    GetPot cl(argc, argv);
    std::string filename = cl.follow("data.json", "-f");
    auto pos = filename.find(".json");
    // check I found the input parameter file
    if(pos == std::string::npos)
    {
        std::cerr << "missing input json file" << std::endl;
        return 1;
    }

    // read the method to be used, gradient_method by default
    unsigned method = cl.follow(0, "-m");

    bool verbose = cl.search(1, "-v");

    // read the parameters
    params p;
    p = read_parameters(filename, method);
    std::cout << p << std::endl;
    
    double minimum = compute_minimum(f, df, p, verbose);
    std::cout << "Minimum: " << minimum << std::endl;

    return 0;
}

double f(const std::vector<double> &x){
    return x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0];
}

std::vector<double> df(const std::vector<double> &x){
    double df1 = x[1] + 12*x[0]*x[0]*x[0] + 3;
    double df2 = x[0] + 2*x[1];
    return std::vector<double>{df1, df2};
}

