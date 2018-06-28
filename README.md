# How to solve an LQ problem with time-varying targets by Riccati's theory

In this short tutorial, we explain how to use Riccati's theory to solve an LQ control problem with targets. The related MATLAB code is downloadable freely.

We consider the optimal control problem:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cmin_%7Bu%5Cin%20L%5E2%280%2CT%29%7DJ%28u%29%3D%5Cfrac12%20%5Cleft%5B%20%5Cint_0%5ET%20%5C%7Cu%28t%29-q%28t%29%5C%7C%5E2%20dt&plus;%5Cbeta%5Cint_0%5ET%20%5C%7CC%28x%28t%29-z%28t%29%29%5C%7C%5E2%20dt&plus;%5Cgamma%20%5C%7CD%28x%28T%29-z%28T%29%29%5C%7C%5E2%5Cright%5D%2C"></p>

where

<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20%5Cfrac%7Bd%7D%7Bdt%7Dx%28t%29&plus;Ax%28t%29%3DBu%28t%29%5Chspace%7B0.6%20cm%7D%20%26%20t%5Cin%20%280%2CT%29%5C%5C%20x%280%29%3Dx_0.%20%5Cend%7Bcases%7D"></p>

In the above control problem, <img src="https://latex.codecogs.com/gif.latex?A%5Cin%5Cmathcal%7BM%7D_%7Bn%5Ctimes%20n%7D">, <img src="https://latex.codecogs.com/gif.latex?B%5Cin%20%5Cmathcal%7BM%7D_%7Bn%5Ctimes%20m%7D">, <img src="https://latex.codecogs.com/gif.latex?C%5Cin%20%5Cmathcal%7BM%7D_%7Br%5Ctimes%20n%7D"> and <img src="https://latex.codecogs.com/gif.latex?D%5Cin%5Cmathcal%7BM%7D_%7Br%5Ctimes%20n%7D">. The control <img src="https://latex.codecogs.com/gif.latex?u%3A%5B0%2CT%5D%5Clongrightarrow%20%5Cmathbb%7BR%7D%5Em">, while the state <img src="https://latex.codecogs.com/gif.latex?x%3A%5B0%2CT%5D%5Clongrightarrow%20%5Cmathbb%7BR%7D%5En">. The control target is <img src="https://latex.codecogs.com/gif.latex?q%5Cin%20C%5E1%28%5B0%2CT%5D%3B%5Cmathbb%7BR%7D%5Em%29"> and the state target is <img src="https://latex.codecogs.com/gif.latex?z%5Cin%20C%5E1%28%5B0%2CT%5D%3B%5Cmathbb%7BR%7D%5En%29">. <img src="https://latex.codecogs.com/gif.latex?%5Cbeta%5Cgeq%200"> and <img src="https://latex.codecogs.com/gif.latex?%5Cgamma%5Cgeq%200"> are positive parameters.

By the Direct Methods in the Calculus of Variations and strict convexity, the above problem admits an unique optimal control.

We compute the optimal pair (optimal control, optimal state) by using the well-known Riccati's theory (see, for instance, [[1](https://epubs.siam.org/doi/pdf/10.1137/130907239), Lemma 2.6] and [[2](https://www.ljll.math.upmc.fr/trelat/fichiers/livreopt2.pdf), section 4.3]).

For further details regarding the algorithm, we refer to [RiccatiAlgorithm.pdf](https://github.com/ChairOfComputationalMathematics/RiccatiLQ/blob/master/RiccatiAlgorithm.pdf).

## Example

Take
<p align="center"><img src="https://latex.codecogs.com/gif.latex?A%3D%20%5Cbegin%7Bpmatrix%7D%202%26-1%5C%5C%20-1%262%20%5Cend%7Bpmatrix%7D%2C%5Chspace%7B0.2%20cm%7DB%3D%20%5Cbegin%7Bpmatrix%7D%201%5C%5C%200%20%5Cend%7Bpmatrix%7D%2C%5Chspace%7B0.2%20cm%7DC%3D%20%5Cbegin%7Bpmatrix%7D%201%260%5C%5C%200%261%20%5Cend%7Bpmatrix%7D%2C%5Chspace%7B0.2%20cm%7D%5Cmbox%7Band%7D%5Chspace%7B0.2%20cm%7DD%3D%20%5Cbegin%7Bpmatrix%7D%200%260%5C%5C%200%260%20%5Cend%7Bpmatrix%7D."></p>

Choose <img src="https://latex.codecogs.com/gif.latex?%5Cbeta%3D26">, <img src="https://latex.codecogs.com/gif.latex?%5Cgamma%3D0">, <img src="https://latex.codecogs.com/gif.latex?x_0%3D%5B1.4%3B1.4%5D">, <img src="https://latex.codecogs.com/gif.latex?q%5Cequiv%200">, <img src="https://latex.codecogs.com/gif.latex?z%28t%29%3D%5B%5Csin%28t%29%3B%5Csin%28t%29%5D"> and T=10. We obtain the following figures:

<p align="center">
  <img src="state_1.png">
</p>

<p align="center">Click
<a href="https://github.com/ChairOfComputationalMathematics/RiccatiLQ/blob/master/state_1.png" target="_blank">here</a> to open state_1.png.
</p>

<p align="center">
  <a href="https://github.com/ChairOfComputationalMathematics/RiccatiLQ/blob/master/state_2.png" target="_blank"><img src="state_2.png">
</p>

<p align="center">
  <a href="https://github.com/ChairOfComputationalMathematics/RiccatiLQ/blob/master/state_2.png" target="_blank">Click here to open state_2.png.</a>
</p>

<p align="center">
  <a href="https://github.com/ChairOfComputationalMathematics/RiccatiLQ/blob/master/control.png" target="_blank"><img src="control.png">
</p>

<p align="center">
  <a href="https://github.com/ChairOfComputationalMathematics/RiccatiLQ/blob/master/control.png" target="_blank">Click here to open control.png.</a>
</p>



Since the parameter <img src="https://latex.codecogs.com/gif.latex?%5Cbeta"> is large enough and the control acts only on the first component of the state equation
* the first component of the state is close to the target;
* the second component of the state is less close to the target;
* the control is far from its target.

The algorithm described in this guide can be employed to test the fulfillment of the turnpike property (see, e.g., [[1](https://epubs.siam.org/doi/pdf/10.1137/130907239)] and [[3](https://arxiv.org/abs/1402.3263)]). In agreement with the theory, the turnpike effect is evident if:
* the targets are constants;
* (A,B) is controllable;
* (A,C) is observable, <img src="https://latex.codecogs.com/gif.latex?%5Cbeta%3E0"> and <img src="https://latex.codecogs.com/gif.latex?%5Cgamma%3D0">;
* the time horizon T is large enough.





## Author

* **Dario Pighin**

## References

[[1](https://epubs.siam.org/doi/pdf/10.1137/130907239)] A. PORRETTA and E. ZUAZUA, _Long time versus steady state optimal control_, SIAM Journal
on Control and Optimization, 51 (2013), pp. 4242–4273.

[[2](https://www.ljll.math.upmc.fr/trelat/fichiers/livreopt2.pdf)] E. TRÉLAT, _Contrôle optimal: théorie & applications_, Vuibert, 2008.

[[3](https://arxiv.org/abs/1402.3263)] E. TRÉLAT and E. ZUAZUA, _The turnpike property in finite-dimensional nonlinear optimal control_, Journal of Differential Equations, 258 (2015), pp. 81–114.

## Acknowledgments

This project has received funding from the European Research Council (ERC) under the European  Union’s Horizon 2020 research and innovation programme (grant agreement No. 694126-DyCon).
 
[DyCon Webpage](http://cmc.deusto.eus/dycon/)
