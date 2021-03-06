#ifndef RICPAD_DIFFERENTIATE
#define RICPAD_DIFFERENTIATE

namespace solver {

template <typename T>
T differentiate(
        // A function which takes an Eigen Matrix
        const std::function<T(T&)> &f,
        // The variables contained in an Eigen matrix
        const T& x,
        // Step size
        const T& h
      )
{
    T xp(x), xm(x);

    xp += h;
    xm -= h;

    T ans = f(xp) - f(xm);
    ans /= (h*T(2));

    return ans;
};

}; // namespace differentiate
#endif
