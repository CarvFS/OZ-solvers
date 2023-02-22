# Solvers for integral equations in liquid state physics[^1]
[^1]: In development.

In this repository it will be find codes in different languages for 

* Solving Ornstein-Zernike equation using: 
  - Picard iteration method;
    - [x] MATLAB
    - [x] Python
    - [x] C
  - Modified Direct Inversion in the Iterative Subspace (MDIIS) method.
    - [x] MATLAB
    - [x] Python
    - [ ] C

* Solving 1D Reference Interaction Site Model (1D-RISM) using the MDIIS method.
  - [x] MATLAB
  - [ ] Python
  - [ ] C

## Tecnical details

The Fourier transforms are computed according to [1]. For more details of how those methods works, one may check the references [2-4].

## References

[1]: <a href="https://www.sciencedirect.com/science/article/pii/0021999171900210" >Numerical fourier transforms in one, two, and three dimensions for liquid state calculations</a>

[2]: <a href="https://www.scielo.br/j/jbchs/a/Rrhx8PT4FbwSzqctPTX4xky/abstract/?lang=en" >Radial Distribution Function for a Hard Sphere Liquid: A Modified Percus-Yevick and Hypernetted-Chain Closure Relations</a>

[3]: <a href="https://www.sciencedirect.com/science/article/abs/pii/S0378437121003381" >Thermodynamic consistency by a modified Perkus–Yevick theory using the Mittag-Leffler function</a>

[4]: <a href="https://onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1096-987X(19990715)20:9%3C928::AID-JCC4%3E3.0.CO;2-X" >Solution of three-dimensional reference interaction site model and hypernetted chain equations for simple point charge water by modified method of direct inversion in iterative subspace</a>
