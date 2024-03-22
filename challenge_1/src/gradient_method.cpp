#include "gradient_method.hpp"
#include <cmath>
#include <iomanip>

// general function that calls the specific method
double compute_minimum(
    const std::function<double(const std::vector<double> &)> &f,
    const std::function<std::vector<double>(const std::vector<double> &)> &df,
    const params &p,
    bool verbose)
{
    if (p.method == 1)
    {
        return compute_minimum_momentum_1(f,df,p,verbose);
    }
    else if (p.method == 2)
    {
        return compute_minimum_momentum_2(f,df,p,verbose);
    }
    else if (p.method == 3)
    {
        return compute_minimum_momentum_nesterov(f,df,p,verbose);
    }
    // other possible methods...
    return compute_minimum_gradient(f,df,p,verbose);

}


// iterative Gradient Method (default)
double compute_minimum_gradient(
    const std::function<double(const std::vector<double> &)> &f,
    const std::function<std::vector<double> (const std::vector<double> &)> &df,
    const params &p,
    bool verbose)
{   
    // first guess
    std::vector<double> x0 (p.x0);
    std::vector<double> x1 (p.x0);
    unsigned iter{};
    std::cout << std::setprecision(10);
    do{
        // move x0 to the next value computed during the previous iteration
        // RMK: the value owned by x1 after the swap should be discarded, no useful anymore
        std::swap(x0,x1);
        if (verbose){
            std::cout << "\nit = " << iter << "\n";
            std::cout << "x" << iter << " = "; print_vector(x0);
            std::cout << "f(x" << iter << ") = " << f(x0) << "\n";
        }

        /* computing alpha_k
           RMK: THE TEMPLATE PARAMETER SHOULD BE MODIFIED TO CHOOSE
                THE METHOD FOR COMPUTING THE STEP
        */
        double alpha = compute_step<step_type::Line>(p.a0, p.mu, iter, p.sigma,
                                    x0, f, df);
        // std::cout << "alpha" << iter << " = " << alpha << "\n";

        // assigning df(x0) to x1 
        x1 = df(x0);

        // multiply by alpha_k by using std::transform and a lambda
        std::transform(x1.begin(), x1.end(), x1.begin(),
                       [&alpha](auto& c){return c*alpha;});

        // perform x_k - alpha_k*df(x_k) by using std::tranform
        std::transform(x0.cbegin(), x0.cend(), x1.cbegin(),
                       x1.begin(), std::minus<double>());
        
        ++iter;
    } while(!converged(x1,x0,f,p.tol_res,p.tol_x,verbose) && iter <= p.max_it);
    return f(x1);
}


template<step_type T>
double compute_step(double a0, double mu, unsigned int k, double sigma,
                    std::vector<double> x,
                    const std::function<double(std::vector<double>)> &f,
                    const std::function<std::vector<double>(std::vector<double>)> &df)
{
    if constexpr (T == step_type::Expo){
        return a0*exp(-mu*k);
    }
    else if constexpr (T == step_type::Inverse){
        return a0 / (1 + mu*k);
    }
    else if constexpr (T == step_type::Line) {  // armijo rule
        double fx = f(x);
        std::vector<double> new_x(x.size());

        // values needed in the convergence criterium
        std::vector<double> dfx = df(x);
        double dfx_square = norm2_square(df(x));

        // perform first iteration with a0
        a0*=2;
        unsigned count{};
        do{
            a0 /=2;
            new_x = dfx;
            
            // multiply by alpha_k by using std::transform and a lambda
            std::transform(new_x.begin(), new_x.end(), new_x.begin(),
                           [&a0](auto& c){return c*a0;});

            // compute the difference between x and alpha*grad(f)(x)
            std::transform(x.cbegin(), x.cend(), new_x.begin(),
                           new_x.begin(), std::minus<double>());
            ++count;
        } while ( (f(x) - f(new_x) ) < sigma*a0*dfx_square && count < 10);
        return a0;
    } 
    return a0;
}

//  check convergence
bool converged(const std::vector<double> &x1,
               const std::vector<double> &x0,
               const std::function<double (std::vector<double>)> &f, 
               const double tol_res,
               const double tol_x,
               bool verbose 
               )
{
    std::vector<double> x_diff(x0.size());
    std::transform(x1.cbegin(), x1.cend(), x0.cbegin(), x_diff.begin(), std::minus<double>() );
    if (norm2_square(x_diff) < tol_x){
        std::cout << "\nConverged on the step length (= "
                  << norm2_square(x_diff) << ")" << std::endl;
        return true;
    }
    else if (std::abs(f(x1) - f(x0)) < tol_res){
        std::cout << "\nConverged on the residual (=  "
                  << std::abs(f(x1) - f(x0)) << ")" << std::endl;
        return true;
    }
    if (verbose)
        std::cout << "residual = " << std::abs(f(x1) - f(x0)) << std::endl;
    return false;
}


// compute the euclidean norm (squared)
double norm2_square(const std::vector<double> &vec)
{
    double sum {};
    for (double v: vec)
        sum+= v*v;
    return std::sqrt(sum);
};


void print_vector(const std::vector<double> &vec){
    std::cout << "[ ";
    for (auto v: vec){
        std::cout << v << " ";
    }
    std::cout << "]" << std::endl;
}


// EXTRAS -----------------------------------

// first verion of heavy_ball method
double compute_minimum_momentum_1(
    const std::function<double(const std::vector<double> &)> &f,
    const std::function<std::vector<double> (const std::vector<double> &)> &df,
    const params &p,
    bool verbose)
{   
    double mu = p.mu;
    double alpha = p.a0;
    // initial guess
    std::vector<double> x0 (p.x0);
    std::vector<double> x1 (p.x0);

    // initialize d0
    std::vector<double> d0 (df(x0));
    std::transform(d0.begin(),d0.end(),d0.begin(),[&alpha](auto&c){return -c*alpha;});
    std::vector<double> d1(d0);

    unsigned iter{};
    std::cout << std::setprecision(10);
    do{
        // move x0 and d0 to next iteration
        std::swap(x0,x1);
        std::swap(d0,d1);
        if (verbose){
            std::cout << "\nit = " << iter << "\n";
            std::cout << "x" << iter << " = "; print_vector(x0);
            std::cout << "f(x" << iter << ") = " << f(x0) << "\n";
        }

        // computing x_(k+1)
        std::transform(x0.cbegin(),x0.cend(),d0.cbegin(),x1.begin(),std::plus<double>());

        // computind d_(k+1) (mu=0.9, alpha = const)
        std::vector<double> d_temp (df(x1));
        std::transform(d0.begin(),d0.end(),d1.begin(),[&mu](auto&c){return c*mu;});
        std::transform(d_temp.begin(),d_temp.end(),d_temp.begin(),[&alpha](auto&c){return alpha*c;});
        std::transform(d1.begin(),d1.end(),d_temp.begin(),d1.begin(), std::minus<double>());
        
        ++iter;
    } while(!converged(x1,x0,f,p.tol_res,p.tol_x, verbose) && iter <= p.max_it);
    return f(x1);
}

// second version of heavy ball method
double compute_minimum_momentum_2(
    const std::function<double(const std::vector<double> &)> &f,
    const std::function<std::vector<double> (const std::vector<double> &)> &df,
    const params &p,
    bool verbose)
{   
    double mu = p.mu;
    double alpha = p.a0;

    // initialization
    std::vector<double> x0(p.x0);       // initialized with x0
    std::vector<double> df_x0(df(x0));
    std::vector<double> x1 (x0);
    for (unsigned i{}; i < x0.size(); ++i){
        x1[i] -= alpha*df_x0[i];        // computing x1
    }
    std::vector<double> df_x1(df(x1));
    std::vector<double> x2(x1);

    // swap to perform correctly the first iteration
    std::swap(x0,x1);
    unsigned iter{1};
    std::cout << std::setprecision(10);
    do{
        std::swap(x2,x1);   // now x1 contains previous x2
        std::swap(x2,x0);   // now x0 contains previous x1
        
        if (verbose){
            std::cout << "\nit = " << iter << "\n";
            std::cout << "x" << iter << " = "; print_vector(x1);
            std::cout << "f(x" << iter << ") = " << f(x1) << "\n";
        }

        // computing x_(k+1)
        df_x1= df(x1);
        for (unsigned i{}; i < x0.size(); ++i){
            x2[i] = x1[i] - alpha*df_x1[i] + mu*(x1[i] - x0[i]);
        }    
        ++iter;
    } while(!converged(x2,x1,f,p.tol_res,p.tol_x, verbose) && iter <= p.max_it);
    return f(x1);
}

// nersterov for heavy ball method
double compute_minimum_momentum_nesterov(
    const std::function<double(const std::vector<double> &)> &f,
    const std::function<std::vector<double> (const std::vector<double> &)> &df,
    const params &p,
    bool verbose)
{   
    double mu = p.mu;
    double alpha = p.a0;

    // initial guess (unordered for accomplishing the first iteration)
    std::vector<double> x0(p.x0);       // initialized with x0
    std::vector<double> df_x0(df(x0));
    std::vector<double> x1 (x0);        // initialized with x0
    for (unsigned i{}; i < x0.size(); ++i){
        x1[i] -= alpha*df_x0[i];        // computing x1
    }
    std::vector<double> df_x1(df(x1));
    std::vector<double> x2(x1);

    std::swap(x0,x1);
    unsigned iter{1};
    std::cout << std::setprecision(10);
    do{
        std::swap(x2,x1);   // now x1 contains previous x2
        std::swap(x2,x0);   // now x0 contains previous x1
        
        if (verbose){
            std::cout << "\nit = " << iter << "\n";
            std::cout << "x" << iter << " = "; print_vector(x1);
            std::cout << "f(x" << iter << ") = " << f(x1) << "\n";
        }
        
        // computing x_(k+1)
        df_x1= df(x1);
        std::vector<double> y(x0.size());
        for (unsigned i{}; i < x0.size(); ++i)
            y[i] = x1[i] + mu*(x1[i] - x0[i]);
        df_x1= df(y);
        for (unsigned i{}; i < x0.size(); ++i)
            x2[i] = y[i] - alpha*df_x1[i];    
        ++iter;
    } while(!converged(x2,x1,f,p.tol_res,p.tol_x, verbose) && iter <= p.max_it);
    return f(x1);
}