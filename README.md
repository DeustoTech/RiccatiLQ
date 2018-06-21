# How to solve an LQ problem with nonzero targets by Riccati's theory

In this short tutorial, we explain how to use Riccati's theory to solve an LQ control problem with targets. The related MATLAB code is downloadable freely.

We consider the optimal control problem:

<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cmin_%7Bu%5Cin%20L%5E2%280%2CT%29%7DJ%28u%29%3D%5Cfrac12%20%5Cleft%5B%20%5Cint_0%5ET%20%5C%7Cu%28t%29-q%28t%29%5C%7C%5E2%20dt&plus;%5Cbeta%5Cint_0%5ET%20%5C%7CC%28x%28t%29-z%28t%29%29%5C%7C%5E2%20dt&plus;%5Cgamma%20%5C%7CD%28x%28T%29-z%28T%29%29%5C%7C%5E2%5Cright%5D%2C"></p>

where <img src="https://latex.codecogs.com/gif.latex?%5CGamma%20%3D%20%5CGamma_D%20%5Ccup%20%5CGamma_N"> and <img src ="https://latex.codecogs.com/gif.latex?%5CGamma_D%20%5Ccap%20%5CGamma_N%3D%20%5Cemptyset">.

<p align="center"><img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bcases%7D%20%5Cfrac%7Bd%7D%7Bdt%7Dx%28t%29&plus;Ax%28t%29%3DBu%28t%29%5Chspace%7B0.6%20cm%7D%20%26%20t%5Cin%20%280%2CT%29%5C%5C%20x%280%29%3Dx_0.%20%5Cend%7Bcases%7D"></p>

In the above control problem, <img src="https://latex.codecogs.com/gif.latex?A%5Cin%5Cmathcal%7BM%7D_%7Bn%5Ctimes%20n%7D">,
<img src="https://latex.codecogs.com/gif.latex?C%5Cin%20%5Cmathcal%7BM%7D_%7Br%5Ctimes%20n%7D"> and <img src="https://latex.codecogs.com/gif.latex?D%5Cin%5Cmathcal%7BM%7D_%7Br%5Ctimes%20n%7D">. The control <img src="https://latex.codecogs.com/gif.latex?u%3A%5B0%2CT%5D%5Clongrightarrow%20%5Cmathbb%7BR%7D%5Em">, while the state <img src="https://latex.codecogs.com/gif.latex?x%3A%5B0%2CT%5D%5Clongrightarrow%20%5Cmathbb%7BR%7D%5En">. The control target is <img src="https://latex.codecogs.com/gif.latex?q%5Cin%20C%5E1%28%5B0%2CT%5D%3B%5Cmathbb%7BR%7D%5Em%29"> and the state target is <img src="https://latex.codecogs.com/gif.latex?z%5Cin%20C%5E1%28%5B0%2CT%5D%3B%5Cmathbb%7BR%7D%5En%29">. <img src="https://latex.codecogs.com/gif.latex?%5Cbeta%5Cgeq%200"> and <img src="https://latex.codecogs.com/gif.latex?%5Cgamma%5Cgeq%200"> are positive parameters.

By the Direct Methods in the Calculus of Variations and strict convexity, the above problem admits an unique optimal control.

We compute the optimal pair (optimal control, optimal state) by using the well-known Riccati's theory (see, for instance, \cite[Lemma 2.6]{porretta2013long} and \cite[section 4.3]{TAT}).

For more details and references, see the pdf.




## Author

* **Dario Pighin**

## Acknowledgments

This project has received funding from the European Research Council (ERC) under the European  Union’s Horizon 2020 research and innovation programme (grant agreement No. 694126-DyCon).
 
[DyCon Webpage](http://cmc.deusto.eus/dycon/)
