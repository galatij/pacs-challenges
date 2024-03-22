# FIRST CHALLENGE PACS
### Stefano Galati
### A.A. 2023.2024

## Code organization
- parameters.hpp: contains the parameters needed by the method and the overload of the streaming operator for the class params

- read_parameters.hpp: allows to parse the parameters from a json file

- gradient_method.hpp: contains the following functions:
  1. norm2_square: given a vector of double, computes its Euclidean norm (squared)
  2. multiply_scalar_vec: performs multiplication by a scalar
  3. converged: to check convergence of the gradient method
  4. compute_step: depending of the method chosen, computes the step to be used in each ieration of the gradient method
  5. compute_minimum: generic function to be called in main.cpp that calls the specific function depending on the method to be used for the iterative process
  6. compute_minimum_gradient: implements the first scheme (used by default if either the method asked is not valid or if the input parameter file is not found)
  7. compute_minimum_momentum_1: implements the first variant of the heavy_ball scheme
  8. compute_minimum_momentum_2: implements the second variant of the heavy_ball scheme
  9. compute_minimum_momentum_nesterov: implements the Nesterov iteration scheme for the heavy_ball scheme


## Compiling
To compile, you just need to type `make`


## Notes for running the code
In the Makefile, you should change the path of the included repository where the utilities for parsing are stored, i.e.
```bash
.../pacs-examples/Examples/Include'
```

The user may run the code with the following options:
- `-m[0-3]`: to choose the iterative scheme 
- - `0` --> gradient method (default)
- - `1` --> heavy ball method (first variant)
- - `2` --> heavy ball method (second variant)
- - `3` --> heavy ball method (Nesterov iteration)
- `-v`: verbosity (prints all the iterations)
- `-f`: to input a different json file for parameters (default = data.json)