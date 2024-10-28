This program implements a numerical method to solve a system of linear equations 𝐴𝑥=𝑏 using Gaussian elimination with scaled partial pivoting. In terms of numerical stability, scaled partial pivoting reduces errors that arise when matrix rows differ greatly in magnitude. The iterative residual vector refinement minimizes the error, keeping the solution more accurate, especially in ill-conditioned systems.

For a $n \times n$ matrix, the $k$-th step of the algorithm (where $k \in \\{1, \ldots, n-1\\}$), the program operates on a submatrix produced after the deletion of $k-1$ columns and $k-1$ rows. At the beginning of the algorithm, $s_i$ norm is calculated for each matrix row (in this case the maximum norm), and in the $k$-th step we choose the $a_{ik}$-th element from the first column of the submatrix, whose $|a_{ik}|/s_i$ is the largest.

The program consists of the following key components: 
* `elimination` function - performs Gaussian elimination, transforming the matrix 𝐴 into an upper triangular form and solving for the unknown vector 𝑥 using back-substitution
* `solveEquations` function - uses the result of elimination in an iterative refinement process to reduce the residual (error) until it is smaller than a given tolerance 𝜖
* `scale` vector - helps in scaled partial pivoting by normalizing the rows for better numerical stability

