#include <vector>
#include <functional>
#include <algorithm>

#include <boost/multiprecision/mpfr.hpp>
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
po::positional_options_description positional;
po::variables_map vm;

const auto help_message = 
        "This program uses the Hankel-Pade method to compute the slope of "
        "the solution of the Thomas-Fermi equation at origin.\n\n"
        "Usage: tf-ricpad MODE OPTIONS\n\n"
        "MODE can be either 'isolated', for the isolated atom, "
        "or 'strong-field', for an atom in a strong field.\n\n";

void print_help_message(std::string const &a = "") {
    std::cout << help_message;
    std::cout << optional << std::endl;
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
        ("d", po::value<int>()->default_value(-1), "Value of d."
         "Set automatically to 3 for the isolated atom equation and to 4 for "
         "the strong field one.")
        /*("no-auto-precision", po::bool_switch(),
         "This program automatically sets the number of digits used in its "
         "computations, the Newton-Raphson step size and its tolerance. "
         "Setting this option uses a fixed value for each. ")*/
        ("tol", po::value<std::string>()->default_value("1E-10"),
         "Starting value of the tolerance for the Newton-Raphson method. "
         "This number gets automatically increased so that it always has "
         "ten more digits than the difference between two consecutive "
         "Hankel determinant roots.")
        ("h", po::value<std::string>()->default_value("1E-20"),
         "Starting value of the step size for the Newton-Raphson method. "
         "This number gets updated automatically to be at least tol^2.")
        ("ndigits", po::value<int>()->default_value(40),
         "Number of digits for the numerical calculations. "
         "The program increases this number automatically when ndigits is less"
         " than -2*log10(h).")
        ("log-nr", po::bool_switch()->default_value(false), 
         "Set this option to print out each Newton-Raphson iteration.")
        ("nr-max-iter", po::value<int>()->default_value(20), 
         "Maximum number of Newton-Raphson iterations.")
        ;

    hidden.add_options()
        ("mode", po::value<std::string>(), "" ) 
        ;

    positional.add("mode", -1);

    opts./*add(mandatory).*/add(optional).add(hidden);

    po::store(
        po::command_line_parser(argc, argv)
        .options(opts).positional(positional).run(), vm
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
    
    // Check if the user asked for help or if they didn't set the 
    // mode correctly.
    const std::vector<std::string> available_modes = {
        "isolated", "strong-field"
    };

    if ( vm.count("mode") ) {
        std::string mode = vm["mode"].as<std::string>();
        if ( std::none_of(
                available_modes.begin(),
                available_modes.end(),
                [&mode](const std::string &a) -> bool {return a == mode;}
                )
           ) {
            std::cout << "Mode " << mode << " not available." 
                << std::endl << std::endl;
            print_help_message("");
        }
    } else {
        print_help_message();
        return 1;
    }

    if ( d == -1 ) {
        if ( vm["mode"].as<std::string>() == "isolated" ) {
            d = 3;
        } else { 
            d = 4;
        }
    }

    // ------------------------------------------------------------------------
    // Here we accomodate to the selected options and/or exit the program if 
    // some of them are wrong
    // ------------------------------------------------------------------------
    
    // Dmin and Dmax should be at least 3
    if ( Dmin < 3 ) {
        std::cout << "Dmin should be at least 3." << std::endl;
        return 1;
    }
    if ( Dmax > -1 && Dmax < 3 ) {
        std::cout << "Dmax should be at least 3." << std::endl;
        return 1;
    }

    // Handle the numerical precision of our computations
    if ( ndigits < 15 ) {
        std::cout << "ndigits should be at least 15." << std::endl;
        return 1;
    } else {
        mpfr_float::default_precision(ndigits);
    }

    // Initial x0 value
    mpfr_float x0;
    x0 = mpfr_float(vm["x0"].as<std::string>());

    int D;

    // Define the lambda function for the Hankel determinants
    std::function<mpfr_float(mpfr_float&)> f;

    if ( vm["mode"].as<std::string>() == "strong-field" ) {
        f = [&D, &d] ( mpfr_float &x ) -> mpfr_float {
            std::vector<mpfr_float> v;

            v = coefs_strong<mpfr_float>(2*D+d, x/2);
            v.erase(v.begin(), v.begin()+d+1);
            return ricpad::hankdet::hankdet<mpfr_float>(D, v);
        };
    } else if ( vm["mode"].as<std::string>() == "isolated" ) {
        f = [&D, &d] ( mpfr_float &x ) -> mpfr_float {
            std::vector<mpfr_float> v;

            v = coefs<mpfr_float>(2*D+d, x/2);
            v.erase(v.begin(), v.begin()+d+1);
            return ricpad::hankdet::hankdet<mpfr_float>(D, v);
        };
    }

    // Here we define the Solver object that will solve the H[D,d] = 0 equation.
    solver::Solver<mpfr_float, mpfr_float> s(f);
    s.set_tol(tol);
    s.set_h(h);
    s.set_maxiter(maxiter);

    // ------------------------------------------------------------------------
    // Here starts the actual computation
    // ------------------------------------------------------------------------

    mpfr_float x, xtry, xold;
    x = x0;
    mpfr_float dE;

    int nfailed = 0;

    for ( D = Dmin; Dmax < 0 || D<=Dmax ; D = D + Dstep ) {
        mpfr_float ans;

        if ( vm["log-nr"].as<bool>() ) {
            s.set_log(ndigits);
        }

        try { 
            xtry = s.solve(x);
            xold = x;
            x = xtry;

            dE = abs(x - xold);

            //if ( ! vm["no-auto-precision"].as<bool>() ) {
                int curr_ndigits = -int(floor(log10(dE)));
                tol = mp::min(tol, dE/1e10);
                h = tol*tol;
                ndigits = std::max(4*curr_ndigits, ndigits);
                ndigits = std::max(-2*int(floor(log10(h))), ndigits);
            //}

            std::cout 
                << "D = " << std::setw(3) << D 
                << " " << std::setw(ndigits+5) 
                       << std::setprecision(ndigits) << std::left << x 
                << " " << std::setw(10) << std::setprecision(4) << dE
                << " digits: " << ndigits 
                << " tol: " << tol 
                << " h: " << h; 

            std::cout << std::endl;

            //if ( ! vm["no-auto-precision"].as<bool>() ) {
                mpfr_float::default_precision(ndigits);
                s.set_tol(tol);
                s.set_h(h);
            //}

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
