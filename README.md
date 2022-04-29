# matrix-fun-cpp
**matrix-fun-cpp** is a small matrix library containing a basic matrix data structure and some basic matrix operations and algorithms written in c++.  
The project was inspired by the *Linear Algebra* course I had at the university, I wanted to see how I could implement some of the discussed matrix operations and algorithms in code.  

The project currently is not well tested, use with caution.

## Building & Running
1. `git clone https://github.com/sbencoding/matrix-fun-cpp.git && cd matrix-fun-cpp`
2. `make`
3. `./main`

This will download the project, build it and run some basic calculation showcasing the functions the library is able to do.

If you wish to use this in your own project, you will need to include some headers.

## Files
* `Makefile` - Responsible for building the project
* `main.cpp` - Small driver program using the matrix data structure and algorithms
* `matrix.hpp` - Defines the `Matrix` data structure and most basic matrix operations (inverse, transpose, exponentiation, multiplication, etc.)
* `eigen.hpp` - Defines functionality to calculate eigenvalues and eigenvectors.
* `approx.hpp` - Defines functionality to project and orthogonalize vectors + calculate least squares solutions of systems
* `solver.hpp` - Defines some solver functions for polynomials (lacks proper polynomial solver)
* `polynomial.hpp` - Defines a basic polynomial, mainly used by eigenvalue calculations internally

## Future tasks
- [ ] Extensive unit testing
- [ ] Look into possible optimizations
- [ ] Generic polynomial solver
- [ ] Documentation
