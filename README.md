# epi-proj
A collection of code written to compute epigraphical/level set projections numerically.

These implementations are based on the algorithms and results of Hoheisel, Burke, and Quang, found in [Link to GenProxNew preprint/whatever](https://math.mcgill.ca/hoheisel/) and [Link to K-epigraph projection paper](https://www.math.mcgill.ca/hoheisel/composite_conjugate.pdf) (currently this points to the Convex Convex-Composite paper but the work in Quang's paper that generalizes the epigraph projection problem is more relevant).

## L1 Ball Projection
Code for computing the projection of a vector onto the <img src="https://render.githubusercontent.com/render/math?math=l_1"> unit ball using a semismooth Newton method. 

## Negative Sum of Logs Projection
Code for computing the projection of a vector onto the epigraph of the negative sum of logs function using a semismooth Newton method.

## Nuclear Norm Projection
Code for computing the projection of a matrix onto the nuclear norm unit ball, a spectral generalization of the <img src="https://render.githubusercontent.com/render/math?math=l_1"> projection based on the the vector of singular values of a given matrix.

## Negative LogDet Projection
Code for computing the projection of a matrix onto the level sets of the function <img src="https://render.githubusercontent.com/render/math?math=X \mapsto -\log\detX">, which is a spectral generalization of the negative sum of logs function based on the vector of eigenvalues of a given matrix.

## Gram Matrix Epigraph Projection
Code for computing the projection of a matrix onto the K-epigraph of the function <img src="https://render.githubusercontent.com/render/math?math=X \mapsto X^TX">
