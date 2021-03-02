# tf-ricpad

A tool for solving the Thomas-Fermi equations using the Hankel-Padé method.

# Installation 

## Requirements

This program is intended to be run on GNU/Linux systes. It uses the following libraries:

* `Boost::multiprecision` (header only, included with regular Boost)
* `Boost::program_options` (compiled library, some distros require separate installation)
* `GMP`
* `MPFR`

Additionally, it requires cmake for installation.

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
