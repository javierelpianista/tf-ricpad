# tf-ricpad

A tool for solving the Thomas-Fermi equations using the Hankel-Padé method

# Installation 

## Requirements

This program uses the following libraries:
* boost
* Eigen
* GMP
* MPFR

Additionally, it requires cmake for installation

## Instructions

Clone over this repo to your pc, then create a directory where the program will
be compiled, `cd` to that directory, and build using cmake. For example:

```
git clone https://github.com/javierelpianista/tf-ricpad
cd tf-ricpad
mkdir build
cd build
cmake ..
cmake --build .
```

# Usage

Run `tf-ricpad` from the build directory with the `--help` option to see 
instructions on how to use it.
