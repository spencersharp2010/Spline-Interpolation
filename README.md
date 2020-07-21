# Spline Interpolation
The user provides a set of interpolation points, as well as the polynomial degree via a Python interface. The code interactively calculates the control points needed to draw a curve through the provided points using functions in C++ scripts.

This repo also contains various other Python files which plot B-splines in 1D and 2D using shape functions and De Boor's algorithm.

pybind11 is used as the interface between Python and C++. Functions and classes from C++ are exposed using "bindings", and these are contained in the folder python/bindings. Once exposed, these C++ functions and classes can be compiled and accessed from Python. Thus, one can take advantage of the efficiency of C++ by writing the computationally demanding functions in C++ while using the user-friendly aspects of Python to input data, as well as plotting and printing.

**Note**: when compiling on Linux, after running the `make` command, you must also run `make install` to copy Python files over to the build folder. On Windows, you must build the `install` folder for the same aforementioned reasons.
