#include <tf.hpp>

namespace mp = boost::multiprecision;
using mp::mpfr_float;

std::vector<mpfr_float> coefs_strong(int N, mpfr_float f2) {
    std::vector<mpfr_float> fj;

    fj.push_back(mpfr_float(1));
    fj.push_back(mpfr_float(0));
    fj.push_back(mpfr_float(f2));
    fj.push_back(mpfr_float(0));

    for ( int j = 4; j <= N; j++ ) {
        mpfr_float A(0), t3;


        for ( int k = 0; k < j-2; k++ ) {
            t3 = 0;

            A +=  (k+2)*(k+1)*fj[k+2]*fj[j-k-2]
                + (k+1)*(j-k-2)*fj[k+1]*fj[j-k-1]
                - 2*t3;
        }

        if ( j > 4 ) 
            A += -2 * fj[j-5];

        fj.push_back(-A/(j*(j-2)));
    }

    return fj;
}

std::vector<mpfr_float> coefs(int N, mpfr_float f2) {
    std::vector<mpfr_float> fj;

    fj.push_back(mpfr_float(1));
    fj.push_back(mpfr_float(0));
    fj.push_back(mpfr_float(f2));

    for ( int j = 3; j <= N; j++ ) {
        mpfr_float A(0), t3;

        for ( int k = 0; k < j-2; k++ ) {
            t3 = 0;
            for ( int l = 0; l <= k; l++ ) {
                t3 += fj[l]*fj[k-l]*fj[j-k-3];
            }
            A +=  (k+2)*(k+1)*fj[k+2]*fj[j-k-2]
                + (k+1)*(j-k-2)*fj[k+1]*fj[j-k-1]
                - 2*t3;
        }

        fj.push_back(-A/(j*(j-2)));
    }

    return fj;
}
