# epi-proj
A collection of code written to compute epigraphical/level set projections numerically.

These implementations are based on the algorithms and results of Hoheisel, Friedlander, and Goodwin, which can be found at https://arxiv.org/abs/2102.06809.

Note that the code written in C requires the [GNU Scientific Library](https://www.gnu.org/software/gsl/) (GSL) for generating random numbers.
The code in this repository is free and distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html).

## L1 Ball Projection
Code for computing the projection of a vector onto the <img src="https://latex.codecogs.com/svg.latex?l_1" title="l_1" /> unit ball. Found in projl1.c  

## Negative Sum of Logs Projection
Code for computing the projection of a vector onto the epigraph of the negative sum of logs function. Found in projminuslog.c
