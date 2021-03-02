#ifndef RICPAD_DIFFERENTIATE
#define RICPAD_DIFFERENTIATE

namespace solver {

template <typename T, int N>
T differentiate(
        // A function which takes an Eigen Matrix
        const std::function<T(T&)>& f,
        // The variables contained in an Eigen matrix
        const T& x,
        // With respect to which variable we differentiate
        const int k,
        // Step size
        const T& h
      )
{
    Eigen::Matrix<T, N, 1> xp(x), xm(x);

    xp(k) += h;
    xm(k) -= h;

    T ans = f(xp) - f(xm);
    ans /= (h*T(2));

    return ans;
};

}; // namespace differentiate
#endif
