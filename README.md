# NLbalancing
Software that uses polynomials to approximately solve the nonlinear (NL) balancing problem for systems with polynomial nonlinearities.  The description of the NL balancing algorithms are provided in papers [1-6].

[1] J. Borggaard and L. Zietsman, “The quadratic-quadratic regulator problem: approximating feedback controls for quadratic-in-state nonlinear systems,” in 2020 American Control Conference (ACC), Jul. 2020, pp. 818–823. doi: 10.23919/ACC45564.2020.9147286

[2] J. Borggaard and L. Zietsman, “On approximating polynomial-quadratic regulator problems,” IFAC-PapersOnLine, vol. 54, no. 9, pp. 329–334, 2021, doi: 10.1016/j.ifacol.2021.06.090

[3] B. Kramer, S. Gugercin, and J. Borggaard, “Nonlinear balanced truncation: Part 2—model reduction on manifolds,” arXiv, Feb. 2023. doi: 10.48550/ARXIV.2302.02036

[4] B. Kramer, S. Gugercin, J. Borggaard, and L. Balicki, “Nonlinear balanced truncation: Part 1—computing energy functions,” arXiv, Dec. 2022. doi: 10.48550/ARXIV.2209.07645

[5] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞ energy functions for polynomial drift nonlinear systems,” ACC 2024.

[6] N. A. Corbin and B. Kramer, “Computing Solutions to the Polynomial-Polynomial Regulator Problem,” CDC 2024.

[7] N. A. Corbin and B. Kramer, “Scalable computation of 𝓗_∞ energy functions for polynomial control-affine systems,” IEEE TAC, May 2025.

## Installation Notes
Clone these repository: 
```
  git clone https://www.github.com/cnick1/PPR.git
  git clone https://www.github.com/cnick1/KroneckerTools.git
  git clone https://www.github.com/cnick1/NLbalancing.git
```
then modify the path in **setKroneckerToolsPath.m**

The installation can be tested in Matlab by typing
```
>> examplesForPaper3
```
which produces the results for [7]. 

The details of some of our functions and test examples are provided below.  


## How to use this software
We consider polynomial control-affine dynamical systems of the form 
$$\dot{\mathbf{x}}  = \mathbf{A} \mathbf{x} + \sum^\ell_{i=2}\mathbf{F}_ i {\mathbf{x}}^{ⓘ} + \mathbf{B} \mathbf{u}  + \sum^\ell_{i=1} \mathbf{G} _ i({\mathbf{x}}^{ⓘ} \otimes \mathbf{u}),$$
$$\mathbf{y}        = \mathbf{C} \mathbf{x} + \sum^\ell_{i=2} \mathbf{H} _ i {\mathbf{x}}^{ⓘ},$$

where $\mathbf{A} \in \mathbb{R}^{n\times n}$, $\mathbf{F}_ i \in \mathbb{R}^{n\times n^i}$, $\mathbf{B} \in \mathbb{R}^{n\times m}$, $\mathbf{G}_ i \in \mathbb{R}^{n \times mn^i}$, $\mathbf{C} \in \mathbb{R}^{p\times n}$, and $\mathbf{H}_ i \in \mathbb{R}^{p\times n^i}$.
We also assume $[\mathbf{A},\mathbf{B}]$ a controllable pair and $[\mathbf{A},\mathbf{C}]$ a detectable pair.
Letting $\eta = 1-\gamma^{-2}$, where $\gamma$ is the $\mathcal{H}_ \infty$ parameter, the past and future energy functions are defined as 
$$\mathcal{E}_ \gamma^{-}(\mathbf{x}_ 0)  \coloneqq  \min_{\substack{\mathbf{u} \in L_{2}(-\infty, 0] \\ \mathbf{x}(-\infty) = 0,  \mathbf{x}(0) = \mathbf{x}_ 0}}  \frac{1}{2} \int_{-\infty}^{0} \eta \Vert \mathbf{y}(t) \Vert^2  +  \Vert \mathbf{u}(t) \Vert^2 {\rm{d}}t,$$

$$\mathcal{E}_ \gamma^{+}(\mathbf{x}_ 0)  \coloneqq \max_{\substack{\mathbf{u} \in L_{2}[0,\infty) \\ \mathbf{x}(0) = \mathbf{x}_ 0,   \mathbf{x}(\infty) = 0}}  \frac{1}{2} \int_{0}^{\infty} \Vert \mathbf{y}(t) \Vert^2  +  \frac{\Vert \mathbf{u}(t) \Vert^2}{\eta} {\rm{d}}t.$$

and they solve the Hamilton-Jacobi-Bellman Partial Differential Equations (HJB PDEs)
$$0 =  \frac{\partial \mathcal{E}_ \gamma^{-}(\mathbf{x})}{\partial \mathbf{x}} \mathbf{f}(\mathbf{x}) + \frac{1}{2}  \frac{\partial \mathcal{E}_ \gamma^{-}(\mathbf{x})}{\partial \mathbf{x}} \mathbf{g}(\mathbf{x}) \mathbf{g}(\mathbf{x})^\top \frac{\partial^\top \mathcal{E}_ \gamma^{-}(\mathbf{x})}{\partial \mathbf{x}} - \frac{\eta}{2}  \mathbf{h}(\mathbf{x})^\top  \mathbf{h}(\mathbf{x})$$
$$0 =  \frac{\partial \mathcal{E}_ \gamma^{+}(\mathbf{x})}{\partial \mathbf{x}} \mathbf{f}(\mathbf{x})   - \frac{\eta}{2} \frac{\partial \mathcal{E}_ \gamma^{+}(\mathbf{x})}{\partial \mathbf{x}} \mathbf{g}(\mathbf{x}) \mathbf{g}(\mathbf{x})^\top \frac{\partial^\top \mathcal{E}_\gamma^{+}(\mathbf{x})}{\partial \mathbf{x}} + \frac{1}{2}\mathbf{h}(\mathbf{x})^\top \mathbf{h}(\mathbf{x})$$

As $\eta$ goes to one, the closed-loop HJB past and future energy functions are recovered, whereas as $\eta$ goes to zero, the open-loop nonlinear controllability and observability energy functions are recovered.
We use Al'brekht's method to solve the HJB PDEs for polynomial expansions of the energy functions 
$$\mathcal{E}_ \gamma^-(\mathbf{x}) = \frac{1}{2} \sum_{i=2}^d \mathbf{v}_ i^\top {\mathbf{x}}^{i} \qquad \text{and} \qquad \mathcal{E}_ \gamma^+(\mathbf{x}) = \frac{1}{2} \sum_{i=2}^d \mathbf{w}_i^\top {\mathbf{x}}^{i}, $$
though a scalable implementation has not been provided until this work.

For a given set of polynomial dynamics defined by the cell arrays `f,g,h` and a permissible value of `eta`, the functions `approxFutureEnergy()` and `approxPastEnergy()` will return the energy function polynomial coefficients $\mathbf{v}_i$ and $\mathbf{w}_i$ up to degree $d=$`degree`:
```
>>  [w] = approxFutureEnergy(f,g,h,eta,degree);
>>  [v] = approxPastEnergy(f,g,h,eta,degree);
```
`approxFutureEnergy()` and `approxPastEnergy()` correspond to Algorithm 1 in reference [1] or [6].
The variables `f,g,h` are cell arrays containing the polynomial coefficients for the dynamics, i.e. 
``` 
>>   f = {A, F2,...};
>>   g = {B, G1, G2,...};
>>   h = {C,H2,...};
```
To aid in computing these polynomial coefficients, we provide the function `utils/approxPolynomialDynamics.m`; given symbolic expressions for $\mathbf{f}(\mathbf{x}),\mathbf{g}(\mathbf{x}),\mathbf{h}(\mathbf{x})$, the function will return the polynomial coefficients to degree $d$ in Kronecker product form using multivariate Taylor series expansions.

The returned variables `v` and `w` are cell arrays with `v{2}` being a vector of dimension $n^2 \times 1$, up to `v{degree+1}` which is a vector of dimension $n^{d+1} \times 1$. 
These can be thought of as vectorized tensors; for example `v{2}=V2(:)`. 
Alternatively, `v{k}` are often reshaped as $n \times n^{k}$ matrices for efficient computations. 

From an initial `x0`, we can compute the approximation to the energy function as
```
>>  E = (1/2)*( v{2}*kron(x0,x0) + ... + v{degree+1}*kron(kron(... ,x0),x0) );
```
or, using the utility function,
```
>>  E = (1/2)*kronPolyEval(v,x0,degree);
```
TODO: Document input-normal and output-diagonal transformations.


## Description of Files
#### setKroneckerToolsPath

Defines the path to the KroneckerTools directory containing functions for working with Kronecker product expressions.
KroneckerTools can be downloaded from github.com/cnick1/KroneckerTools.
The default assumes that NLbalancing and KroneckerTools lie in the same directory and uses relative pathnames.
This should be changed if you use different locations. 

#### LyapProduct

Efficiently computes the product of a special Kronecker sum matrix (aka an N-Way Lyapunov matrix) with a vector.  
This is done by reshaping the vector, performing matrix-matrix products, then reshaping the answer.  
We could also utilize the matrization of the associated tensor.

## Examples

### runExample1.m
Approximates the future and past energy functions for a one-dimensional model problem motivated by the literature.  
A quadratic approximation to this problem appears as Example 1 in [1], and the full polynomial problem appears as Example 1 in [6].

### runExample2.m
The example is based on a two-dimensional problem found in Kawano and Scherpen, IEEE Transactions on Automatic Control, 2016.
This example approximates the future and past energy functions, then computes an approximation to the input-normal transformation.  
An quadratic approximation to the original model is considered as Example 2 in [1], and the full quadratic-bilinear model is considered as Example 2 in [6].

### runExample3
This example demonstrates the scalability and convergence of the proposed approach on a finite-element discretization of Burgers' equation.
This is Example 3 in [1].

### runExample4
This example demonstrates the scalability and convergence of the proposed approach on a finite-element discretization of the Kuramoto-Sivashinsky equation.
This is Example 4 in [1].

### runExample6
This example demonstrates the scalability and convergence of the proposed approach on a finite-element discretization of a nonlinear beam.
This is Example 3 in [6].

### runExample7
This example demonstrates controllers based on the energy functions on a 3D aircraft stall stabilization model from Garrard 1977.

### runExample8
This example demonstrates the scalability and convergence of the proposed approach on a finite-element discretization of a nonlinear heat equation (reaction-diffusion problem).
This is Example 1 in [5].

<!-- ### runExample12
This example is for testing the output-diagonalization transformation.
The system is a polynomial approximation of the 2D model from Fujimoto and Scherpen 2001, 2005, 2010. -->


## Algorithms from Kramer, Gugercin, and Borggaard, Part 2:

### Algorithm 1 is implemented in _inputNormalTransformation.m_

### Algorithm 2 is implemented in _approximateSingularValueFunctions.m_

## References
```
@inproceedings{Borggaard2020,
  author           = {Borggaard, Jeff and Zietsman, Lizette},
  booktitle        = {2020 American Control Conference (ACC)},
  doi              = {10.23919/ACC45564.2020.9147286},
  month            = jul,
  pages            = {818--823},
  title            = {The quadratic-quadratic regulator problem: approximating feedback controls for quadratic-in-state nonlinear systems},
  year             = {2020}
}
```
```
@article{Borggaard2021,
  author           = {Borggaard, Jeff and Zietsman, Lizette},
  doi              = {10.1016/j.ifacol.2021.06.090},
  journal          = {{IFAC}-{PapersOnLine}},
  number           = {9},
  pages            = {329--334},
  publisher        = {Elsevier {BV}},
  title            = {On approximating polynomial-quadratic regulator problems},
  volume           = {54},
  year             = {2021}
}
```
```
@unpublished{Kramer2023,
  archiveprefix    = {arXiv},
  author           = {Kramer, Boris and Gugercin, Serkan and Borggaard, Jeff},
  doi              = {10.48550/ARXIV.2302.02036},
  eprint           = {2302.02036},
  month            = feb,
  note             = {{\em arXiv:2302.02036}},
  primaryclass     = {math.OC},
  publisher        = {arXiv},
  title            = {Nonlinear balanced truncation: {P}art 2---model reduction on manifolds},
  year             = {2023}
}
```
```
@unpublished{Kramer2022,
  archiveprefix    = {arXiv},
  author           = {Kramer, Boris and Gugercin, Serkan and Borggaard, Jeff and Balicki, Linus},
  doi              = {10.48550/ARXIV.2209.07645},
  eprint           = {arXiv:2209.07645v2},
  month            = dec,
  note             = {{\em arXiv:2209.07645v2}},
  primaryclass     = {math.OC},
  publisher        = {arXiv},
  title            = {Nonlinear balanced truncation: {P}art 1---computing energy functions},
  year             = {2022}
}
```
```
@Unpublished{Corbin2023,
  author           = {Corbin, Nicholas A. and Kramer, Boris},
  note             = {Submitted to ACC 2024},
  title            = {Scalable computation of {$\mathcal{H}_\infty$} energy functions for polynomial drift nonlinear systems},
  year             = {2023},
}
```
```
@Unpublished{Corbin2023a,
  author           = {Corbin, Nicholas A. and Kramer, Boris},
  note             = {Submitted to IEEE Transactions on Automatic Control},
  title            = {Scalable computation of {$\mathcal{H}_\infty$} energy functions for polynomial control-affine systems},
  year             = {2023},
}
```
