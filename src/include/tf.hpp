#ifndef TF_HPP
#define TF_HPP

#include <vector>
#include <boost/multiprecision/mpfr.hpp>

namespace mp = boost::multiprecision;
using mp::mpfr_float;

template <typename num_t>
std::vector<num_t> coefs_strong(int N, num_t f2) {
    std::vector<num_t> fj;

    fj.push_back(num_t(1));
    fj.push_back(num_t(0));
    fj.push_back(num_t(f2));
    fj.push_back(num_t(0));

    for ( int j = 4; j <= N; j++ ) {
        num_t A(0), t3;


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

template <typename num_t>
std::vector<num_t> coefs(int N, num_t f2) {
    std::vector<num_t> fj;

    fj.push_back(num_t(1));
    fj.push_back(num_t(0));
    fj.push_back(num_t(f2));

    for ( int j = 3; j <= N; j++ ) {
        num_t A(0), t3;

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
#endif
