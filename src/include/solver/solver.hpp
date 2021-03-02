#ifndef SOLVER
#define SOLVER

#include <iostream>
#include <stdexcept>
#include <math.h>

#include <solver/differentiate.hpp>

namespace solver { 
    //
// Solver for N equations using a real type R
template <
    typename C, // complex number type
    typename R  // real number type (for tolerance parameters, etc.)
>
class Solver {
    private:
        // Function to solve
        const std::function<C(C&)> f_;
        //Eigen::Matrix<C,N,1> x0;
        // Accepted difference between two iterations
        R tol_; 
        // Numerical value for the differentiation step (only relevant when
        // using numerical differentiation)
        R h_;
        // Maximum number of iterations
        int maxiter_ = 100;
        // Print each iteration?
        bool log_iters_ = false;
        int  log_precision_ = 15;

    public:
        //----------------------------------------------------------------------
        // Construct a Solver object from a vector of functions. 
        // Differentiation is performed numerically.
        Solver(const std::function<C(C&)> &f) :
            
            f_(f) {
                //h_   = std::sqrt(std::numeric_limits<R>::epsilon());
                h_   = 1e-12;
                tol_ = 1e-8;
            };

        //----------------------------------------------------------------------
        // Setters
        void set_h(R h) {h_ = h;};
        void set_tol(R tol) {tol_ = tol;};
        void set_maxiter(int maxiter) {maxiter_ = maxiter;};
        void set_log(int precision) {
            log_iters_ = true;
            log_precision_ = precision;
        };
        void unset_log() {
            log_iters_ = false;
        }
        
        //----------------------------------------------------------------------
        // Getters
        int maxiter() {return maxiter_;};

        //----------------------------------------------------------------------
        // Solve for f using x0 as initial value
        C solve(C x0) 
        {

            using solver::differentiate;
            C jacobian, inv_jacobian;
            C x(x0), xold;
            R desv = tol_ + 1;

            int niter = 0;

            C val;

            while ( desv > tol_ ) {
                val = differentiate<C>(f_, x, h_);
                jacobian = std::move(val);

                inv_jacobian = 1/jacobian;
                
                C F;

                F = f_(x); 

                xold = x;
                x = x - inv_jacobian * F;

                desv = abs(x - xold);

                if ( log_iters_ ) {
                    std::cout << std::setw(8) << "( NR: " << niter << " )";
                    std::cout 
                        << std::setprecision(log_precision_) 
                        << std::setw(log_precision_ + 10) << std::left
                        << x;
                    std::cout << std::endl;
                }


                if ( niter++ > maxiter_ ) {
                    throw std::runtime_error(
                            "Maximum number of iterations reached."
                            );
                    return x;
                }
            }

            return x;
        };

};

}; // namespace solver

#endif
