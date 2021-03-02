#include <vector>
#include <functional>
#include <algorithm>

#include <Eigen/Dense>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/eigen.hpp>
#include <boost/program_options.hpp>

#include <ricpad/hankdet.hpp>
#include <solver/solver.hpp>
#include <tf.hpp>

namespace mp = boost::multiprecision;
using mp::mpfr_float;
using mp::mpfr_float;

namespace po = boost::program_options;

po::options_description 
    //mandatory("Required options"), 
    optional("Non-mandatory options"), 
    hidden, opts;
po::variables_map vm;

const auto help_message = "Usage: tf-ricpad OPTIONS\n\n";

void print_help_message(std::string const &a = "") {
    cout << help_message;
    //cout << mandatory << endl;
    cout << optional << endl;
}

int main(int argc, char* argv[]) {
    optional.add_options()
        ("help", po::value<std::string>()
         ->zero_tokens()
         ->implicit_value("")
         ->notifier(&print_help_message)
         , "Print this message")
        ("x0", po::value<std::string>()->default_value("-1.6"), 
         "Initial value of the second coefficient in the expansion of f(x).")
        ("Dmin", po::value<int>()->default_value(3), "Starting D value.")
        ("Dmax", po::value<int>()->default_value(-1), 
         "Maximum D value. If Dmax == -1, Dmax is assumed to be infinite.")
        ("Dstep", po::value<int>()->default_value(1), 
         "Distance between D values.")
        ("d", po::value<int>()->default_value(3), "Value of d.")
        ("no-auto-precision", po::bool_switch(),
         "This program automatically sets the number of digits used in its "
         "computations, the Newton-Raphson step size and its tolerance. "
         "Setting this option uses a fixed value for each. ")
        ("ndigits", po::value<int>()->default_value(40),
         "Number of digits for the numerical calculations. "
         "If auto-precision is selected, this accounts for the starting number"
         " of digits.")
        ("tol", po::value<std::string>()->default_value("1E-10"),
         "Starting value of the tolerance for the Newton-Raphson method. "
         "If --no-auto-precision is set, this value is used throughout the "
         "whole computation." )
        ("h", po::value<std::string>()->default_value("1E-20"),
         "Starting value of the step size for the Newton-Raphson method. "
         "If --no-auto-precision is set, this value is used throughout the "
         "whole computation." )
        ("log-nr", po::bool_switch()->default_value(false), 
         "Set this option to print out each Newton-Raphson iteration")
        ("nr-max-iter", po::value<int>()->default_value(20), 
         "Maximum number of Newton-Raphson iterations")
        ("strong-field", po::bool_switch(),
         "Solve the Thomas-Fermi equation for atoms in a strong magnetic field"
         " instead")
        ;

    opts./*add(mandatory).*/add(optional).add(hidden);

    po::store(
        po::command_line_parser(argc, argv)
        .options(opts).run(), vm
        );

    po::notify(vm);

    // ------------------------------------------------------------------------
    // Here we read the options
    // ------------------------------------------------------------------------
   
    // Did the user ask for help?
    if ( vm.count("help") ) {
        return 1;
    }

    // Dmin
    int Dmin = vm["Dmin"].as<int>();

    // Dmax
    int Dmax = vm["Dmax"].as<int>();

    // Step between D values
    int Dstep = vm["Dstep"].as<int>();

    // d
    int d = vm["d"].as<int>();

    // Number of digits for numerical computations
    int ndigits = vm["ndigits"].as<int>();
    mpfr_float::default_precision(ndigits);
    mpfr_float::default_precision(ndigits);

    // Tolerance for the Newton-Raphson method
    mpfr_float tol;
    tol = mpfr_float(vm["tol"].as<std::string>());

    // Step size for the Newton-Raphson method
    mpfr_float h;
    h = mpfr_float(vm["h"].as<std::string>());

    // Maximum number of Newton-Raphson iterations
    int maxiter = vm["nr-max-iter"].as<int>();

    // ------------------------------------------------------------------------
    // Here we accomodate to the selected options and/or exit the program if 
    // some of them are wrong
    // ------------------------------------------------------------------------
    
    // Dmin and Dmax should be at least 3
    if ( Dmin < 3 ) {
        std::cout << "Dmin should be at least 3" << std::endl;
        return 1;
    }
    if ( Dmax > -1 && Dmax < 3 ) {
        std::cout << "Dmax should be at least 3" << std::endl;
        return 1;
    }

    // Handle the numerical precision of our computations
    if ( ndigits < 15 ) {
        std::cout << "ndigits should be at least 15" << std::endl;
        return 1;
    } else {
        mpfr_float::default_precision(ndigits);
    }

    // Initial x0 value
    Eigen::Matrix<mpfr_float,1,1> x0;
    x0(0) = mpfr_float(vm["x0"].as<std::string>());

    int D;

    // Define the lambda function for the Hankel determinants
    std::function<
        mpfr_float(Eigen::Matrix<mpfr_float, 1, 1>&)
        > f;
    std::vector<mpfr_float> v;

    if ( vm["strong-field"].as<bool>() ) {
        f = [&D, &d] ( Eigen::Matrix<mpfr_float,1,1> &x ) -> mpfr_float {
            std::vector<mpfr_float> v;

            v = coefs_strong<mpfr_float>(2*D+d, x(0));
            v.erase(v.begin(), v.begin()+d+1);
            return ricpad::hankdet::hankdet<mpfr_float>(D, v);
        };
    } else {
        f = [&D, &d] ( Eigen::Matrix<mpfr_float,1,1> &x ) -> mpfr_float {
            std::vector<mpfr_float> v;

            v = coefs<mpfr_float>(2*D+d, x(0));
            v.erase(v.begin(), v.begin()+d+1);
            return ricpad::hankdet::hankdet<mpfr_float>(D, v);
        };
    }

    std::vector<decltype(f)> vf;
    vf.push_back(f);

    // Here we define the Solver object that will solve the H[D,d] = 0 equation.
    solver::Solver<mpfr_float, mpfr_float, 1> s(vf);
    s.set_tol(tol);
    s.set_h(h);
    s.set_maxiter(maxiter);

    // ------------------------------------------------------------------------
    // Here starts the actual computation
    // ------------------------------------------------------------------------

    Eigen::Matrix<mpfr_float,1,1> x, xtry, xold;
    x = x0/2;
    mpfr_float dE;

    int nfailed = 0;

    for ( D = Dmin; Dmax < 0 || D<=Dmax ; D = D + Dstep ) {
        Eigen::Matrix<mpfr_float,1,1> ans;

        if ( vm["log-nr"].as<bool>() ) {
            s.set_log(ndigits);
        }

        try { 
            xtry = s.solve(x);
            xold = x;
            x = xtry;

            dE = abs(x(0) - xold(0));


            if ( ! vm["no-auto-precision"].as<bool>() ) {
                int curr_ndigits = -int(floor(log10(dE)));
                tol = mp::min(tol, dE/1e10);
                h = tol*tol;
                ndigits = std::max(4*curr_ndigits, ndigits);
                ndigits = std::max(-2*int(floor(log10(h))), ndigits);
            }

            std::cout 
                << "D = " << std::setw(3) << D 
                << " " << std::setw(ndigits+5) 
                       << std::setprecision(ndigits) << std::left << 2*x(0,0) 
                << " " << std::setw(10) << std::setprecision(4) << dE
                << " digits: " << ndigits 
                << " tol: " << tol 
                << " h: " << h; 

            std::cout << std::endl;

            if ( ! vm["no-auto-precision"].as<bool>() ) {
                mpfr_float::default_precision(ndigits);
                s.set_tol(tol);
                s.set_h(h);
            }

        } catch ( const std::runtime_error& e ) {
            nfailed += 1;
            std::cout 
                << "Newton-Raphson failed after " << s.maxiter() 
                << " iterations for D = " << D << "." << std::endl;

            if ( nfailed >= 3 ) {
                throw std::runtime_error(
                    "The Newton-Raphson method failed to converge for " + 
                    std::to_string(nfailed) + " consecutive D values. "
                    );
            }
        }
    }

    return 0;
}
