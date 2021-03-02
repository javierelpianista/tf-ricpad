# tf-ricpad

A tool for solving the Thomas-Fermi equations using the Hankel-Padé method

# Usage

Usage: tf-ricpad OPTIONS

Non-mandatory options:
  --help                  Print this message
  --x0 arg (=-1.6)        Initial value of the second coefficient in the
                          expansion of f(x).
  --Dmin arg (=3)         Starting D value.
  --Dmax arg (=-1)        Maximum D value. If Dmax == -1, Dmax is assumed to be
                          infinite.
  --Dstep arg (=1)        Distance between D values.
  --d arg (=3)            Value of d.
  --no-auto-precision     This program automatically sets the number of digits
                          used in its computations, the Newton-Raphson step
                          size and its tolerance. Setting this option uses a
                          fixed value for each.
  --ndigits arg (=40)     Number of digits for the numerical calculations. If
                          auto-precision is selected, this accounts for the
                          starting number of digits.
  --tol arg (=1E-10)      Starting value of the tolerance for the
                          Newton-Raphson method. If --no-auto-precision is set,
                          this value is used throughout the whole computation.
  --h arg (=1E-20)        Starting value of the step size for the
                          Newton-Raphson method. If --no-auto-precision is set,
                          this value is used throughout the whole computation.
  --log-nr                Set this option to print out each Newton-Raphson
                          iteration
  --nr-max-iter arg (=20) Maximum number of Newton-Raphson iterations
  --strong-field          Solve the Thomas-Fermi equation for atoms in a strong
                          magnetic field instead
