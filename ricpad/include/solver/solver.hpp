#ifndef SOLVER
#define SOLVER

#include <iostream>
#include <vector>
#include <stdexcept>
#include <math.h>

#include <Eigen/Dense>

#include <solver/differentiate.hpp>

using std::cout;
using std::endl;

namespace solver { 
    //
// Solver for N equations using a real type R
template <
    typename C, // complex number type
    typename R, // real number type (for tolerance parameters, etc.)
    int N       // number of equations and unknowns
>
class Solver {
    private:
        // Vector with the functions
        const std::vector<std::function<C(Eigen::Matrix<C,N,1>&)>> f_;
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
        Solver(
            const std::vector<std::function<C(Eigen::Matrix<C,N,1>&)>> &f
            ) :
            
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
        Eigen::Matrix<C,N,1> solve(
            Eigen::Matrix<C,N,1> x0
            ) 
        {

            using solver::differentiate;
            Eigen::Matrix<C, N, N> jacobian, inv_jacobian;
            Eigen::Matrix<C, N, 1> x(x0), xold;
            R desv = tol_ + 1;

            int niter = 0;

            C val;

            while ( desv > tol_ ) {
                for ( int i = 0; i < x.size(); i++ ) {
                    for ( int j = 0; j < x.size(); j++ ) {
                        val = differentiate<C, N>(f_[i], x, j, h_);
                        jacobian(i, j) = std::move(val);
                    }
                }

                inv_jacobian = jacobian.inverse();
                
                Eigen::Matrix<C, N, 1> F;

                for ( int i = 0; i < N; i++ ) F(i) = f_[i](x); 

                xold = x;
                x = x - inv_jacobian * F;

                desv = (x - xold).norm();

                if ( log_iters_ ) {
                    std::cout << std::setw(8) << "( NR: " << niter << " )";
                    for ( int i = 0; i < N; i++ ) {
                        std::cout 
                            << std::setprecision(log_precision_) 
                            << std::setw(log_precision_ + 10) << std::left
                            << x(i);
                    }
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

        /* Construct a Solver object from a function and its derivative
        Solver(
            const std::function<R(R)> f,
            const std::function<R(R)> df
            )
        */
};

}; // namespace solver

#endif
