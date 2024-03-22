#include "read_parameters.hpp"

#include <fstream>
#include <vector>
#include "json.hpp"

using json = nlohmann::json;
params read_parameters(const std::string &filename, const unsigned method){
    std::ifstream check(filename);
    params p;
    // check i opened the file
    if(!check)
    {
        std::cerr << "ERROR: Parameter file " << filename
                  << " does not exist" << std::endl;
        std::cerr << "Reverting to default values (gradient method)" << std::endl;
        p.x0 = std::vector<double> {0,0};
        check.close();

        // return default values
        return p;
    }
    else
        check.close();
    
    p.method = method;
    std::ifstream f(filename);
    json data = json::parse(f);
    p.max_it = data.value("max_it", 10);
    p.tol_res = data.value("tol_res", 0.);
    p.tol_x = data.value("tol_x", 0.);
    p.x0 = data["x0"].get<std::vector<double>>();
    
    if (method >= 1 && method <= 3)
    {
        p.a0 = data["momentum"].value("a0", 0.2);
        p.mu = data["momentum"].value("mu", 0.9);
        switch(method){
            case (1):
                p.var = data["momentum"].value("var", heavy_ball_variant::variant_1);
                break;
            case (2):
                p.var = data["momentum"].value("var", heavy_ball_variant::variant_2);
                break;
            case (3):
                p.var = data["momentum"].value("var", heavy_ball_variant::variant_Nesterov);
                break;
        }
        return p;
    }
    if (method != 0){
        std::cout << "Invalid method number, using default gradient method..." << std::endl;
    }
    p.a0 = data["gradient_method"].value("a0", 0.2);
    p.mu = data["gradient_method"].value("mu", 0.2);
    p.sigma = data["gradient_method"].value("sigma", 0.2);
    p.step = data["gradient_method"].value("step", step_type::Expo);

    return p;
}
