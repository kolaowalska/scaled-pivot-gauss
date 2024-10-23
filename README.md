This code implements a numerical method to solve a system of linear equations 𝐴𝑥=𝑏using Gaussian elimination with scaled partial pivoting. 
The function `elimination` performs Gaussian elimination, transforming the matrix 𝐴 into an upper triangular form and solving for the unknown vector 𝑥 using back-substitution.
The function `solveEquations` uses the result of elimination in an iterative refinement process to reduce the residual (error) until it is smaller than a given tolerance 𝜖.
The `scale` vector helps in scaled partial pivoting by normalizing the rows for better numerical stability.
