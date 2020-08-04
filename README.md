# epi-proj
A collection of code written to compute epigraphical/level set projections numerically.

These implementations are based on the algorithms and results of Hoheisel, Burke, and Quang, found in [Link to GenProxNew preprint/whatever](https://math.mcgill.ca/hoheisel/) and [Link to K-epigraph projection paper](https://www.math.mcgill.ca/hoheisel/composite_conjugate.pdf) (currently this points to the Convex Convex-Composite paper but the work in Quang's paper that generalizes the epigraph projection problem is more relevant).

## L1 Ball Projection
Code for computing the projection of a vector onto the <img src="https://latex.codecogs.com/svg.latex?l_1" title="l_1" /> unit ball. 

## Negative Sum of Logs Projection
Code for computing the projection of a vector onto the epigraph of the negative sum of logs function.

## Nuclear Norm Projection
Code for computing the projection of a matrix onto the nuclear norm unit ball, a spectral generalization of the <img src="https://latex.codecogs.com/svg.latex?l_1" title="l_1" /> projection based on the the vector of singular values of a given matrix.

## Negative LogDet Projection
Code for computing the projection of a matrix onto the level sets of the function <img src="https://latex.codecogs.com/svg.latex?X\mapsto&space;-\log&space;\det&space;X" title="X\mapsto -\log \det X" />, which is a spectral generalization of the negative sum of logs function based on the vector of eigenvalues of a given matrix.

## Gram Matrix Epigraph Projection
Code for computing the projection of a matrix onto the K-epigraph of the function <img src="https://latex.codecogs.com/svg.latex?X\mapsto&space;X^TX" title="X\mapsto X^TX" />. Here K represents the positive semidefinite cone.
