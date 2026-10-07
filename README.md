# ExpHR: Exponential Hermite--Rosenbrock Methods

#### Nguyen Van Hoang<sup>1</sup> and Vu Thai Luan<sup>1</sup>

<sup>1</sup>Department of Mathematics & Statistics, Texas Tech University,  
1108 Memorial Circle, Lubbock, Texas 79409, USA

This repository contains MATLAB implementations and numerical experiments associated with the manuscript:

**Vu Thai Luan and Nguyen Van Hoang, “Exponential Hermite--Rosenbrock Methods for Semilinear PDEs with Strongly Stiff Nonlinearities,” 2026.**

## Overview

Exponential Hermite--Rosenbrock (`ExpHR`) methods are high-order exponential integrators for stiff semilinear parabolic problems of the form

$$
u'(t)=Lu(t)+N(u(t))=:F(u(t)),
\qquad
u(t_0)=u_0,
$$

where \(L\) represents the stiff linear differential operator and \(N\) is a sufficiently smooth nonlinear function that may itself contribute substantially to the stiffness.

The main idea of ExpHR methods is to combine the continuous local linearization used in exponential Rosenbrock methods with the Hermite-type derivative information used in exponential Hermite methods.

At each time step, the full vector field is linearized around the current numerical solution \(u_n\). We define

$$
J_n=F'(u_n)=L+N'(u_n),
$$

and introduce the nonlinear remainder

$$
N_n(u)=F(u)-J_nu.
$$

Its derivative is

$$
N_n'(u)
=
F'(u)-J_n
=
N'(u)-N'(u_n),
$$

and therefore

$$
N_n'(u_n)=0.
$$

The Hermite-type information is introduced through the directional derivative of the nonlinear remainder along the solution,

$$
\frac{d}{dt}N_n(u(t))
=
N_n'(u(t))F(u(t)).
$$

Accordingly, the ExpHR correction is defined by

$$
H_{ni}=N_n'(U_{ni})F(U_{ni}).
$$

The general \(s\)-stage ExpHR method is

$$
\begin{aligned}
U_{ni}
&=
u_n+c_i h_n\varphi_1(c_i h_nJ_n)F(u_n)
+h_n^2\sum_{j=2}^{i-1}a_{ij}(h_nJ_n)H_{nj},
\qquad i=2,\ldots,s,
\\
u_{n+1}
&=
u_n+h_n\varphi_1(h_nJ_n)F(u_n)
+h_n^2\sum_{i=2}^{s}b_i(h_nJ_n)H_{ni}.
\end{aligned}
$$

The matrix functions are defined by

$$
\varphi_0(z)=e^z,
\qquad
\varphi_k(z)
=
\int_0^1
e^{(1-\theta)z}
\frac{\theta^{k-1}}{(k-1)!}\,d\theta,
\qquad k\geq1.
$$

## Relation to ExpH and ExpRB

The main distinction among ExpH, ExpRB, and ExpHR methods is:

- `ExpH`: fixed operator \(L\) + Hermite-type derivative information;
- `ExpRB`: step-dependent Jacobian \(J_n\) + value-based nonlinear remainder corrections;
- `ExpHR`: step-dependent Jacobian \(J_n\) + Hermite-type derivative corrections.

For ExpH methods, the fixed splitting

$$
F(u)=Lu+N(u)
$$

is retained throughout the integration. The matrix functions therefore depend only on \(L\), while Hermite information is introduced through

$$
N'(u)F(u).
$$

Thus, ExpH methods enrich the approximation of the nonlinear term through derivative information, but the nonlinear stiffness remains outside the exponential propagation.

Exponential Rosenbrock methods instead use the full step-dependent Jacobian

$$
J_n=F'(u_n)=L+N'(u_n),
$$

so that stiffness associated with the nonlinear term is incorporated directly into the exponential operator.

Their nonlinear corrections are typically based on value differences of the form

$$
N_n(U_{ni})-N_n(u_n).
$$

ExpHR methods use the same full Jacobian \(J_n\), but replace the value-based correction by the Hermite-type derivative correction

$$
N_n'(U_{ni})F(U_{ni}).
$$

Thus, ExpHR combines the treatment of nonlinear stiffness characteristic of exponential Rosenbrock methods with Hermite-type information about the variation of the nonlinear remainder.

## Stiff order conditions

The stiff order conditions for ExpHR methods are derived up to order five.

The property

$$
N_n'(u_n)=0
$$

plays an important role in the local error expansion and in the derivation of the stiff order conditions.

The resulting framework allows the construction of practical ExpHR schemes of orders three through five.

## Implemented ExpHR methods

The repository contains ExpHR methods for both constant and adaptive time stepping.

### Constant time stepping

For constant step sizes, the method name indicates the order of the scheme and the number of stages.

The implemented constant-step methods are:

- third-order:
  - `expHR3s2`;

- fourth-order:
  - `expHR4s2`,
  - `expHR4s3`;

- fifth-order:
  - `expHR5s3`,
  - `expHR5s3a`,
  - `expHR5s4`,
  - `expHR5s4a`,
  - `expHR5s5`.

For example, `expHR5s4` denotes a fifth-order ExpHR method with four stages.

### Embedded and adaptive time stepping

For adaptive time stepping, embedded ExpHR pairs are used.

The two order digits indicate the orders of the primary and embedded formulas.

The implemented embedded pairs are:

- `expHR32s2`;
- `expHR42s2`;
- `expHR43s3`;
- `expHR53s3`;
- `expHR53s3a`;
- `expHR54s4`;
- `expHR54s4a`;
- `expHR54s5`.

For example, `expHR54s4` denotes a four-stage embedded \(5(4)\) ExpHR pair.

Similarly, `expHR32s2` denotes a two-stage embedded \(3(2)\) ExpHR pair.

## Embedded formulas

The ExpHR schemes are constructed together with embedded formulas for adaptive time stepping.

For an \(s\)-stage method, the internal stages are

$$
U_{ni}
=
u_n+c_i h_n\varphi_1(c_i h_nJ_n)F(u_n)
+h_n^2\sum_{j=2}^{i-1}a_{ij}(h_nJ_n)H_{nj},
\qquad i=2,\ldots,s.
$$

The primary approximation is

$$
u_{n+1}
=
u_n+h_n\varphi_1(h_nJ_n)F(u_n)
+h_n^2\sum_{i=2}^{s}b_i(h_nJ_n)H_{ni},
$$

while the embedded approximation is

$$
\bar u_{n+1}
=
u_n+h_n\varphi_1(h_nJ_n)F(u_n)
+h_n^2\sum_{i=2}^{s}\bar b_i(h_nJ_n)H_{ni}.
$$

Hence, the complete embedded ExpHR pair is

$$
\begin{aligned}
U_{ni}
&=
u_n+c_i h_n\varphi_1(c_i h_nJ_n)F(u_n)
+h_n^2\sum_{j=2}^{i-1}a_{ij}(h_nJ_n)H_{nj},
\\
u_{n+1}
&=
u_n+h_n\varphi_1(h_nJ_n)F(u_n)
+h_n^2\sum_{i=2}^{s}b_i(h_nJ_n)H_{ni},
\\
\bar u_{n+1}
&=
u_n+h_n\varphi_1(h_nJ_n)F(u_n)
+h_n^2\sum_{i=2}^{s}\bar b_i(h_nJ_n)H_{ni}.
\end{aligned}
$$

The primary and embedded formulas share the same internal stages and Hermite corrections.

Therefore, the embedded approximation introduces no additional internal stages.

The local error estimator is

$$
e_{n+1}
=
u_{n+1}-\bar u_{n+1},
$$

or equivalently,

$$
e_{n+1}
=
h_n^2
\sum_{i=2}^{s}
\left(
b_i(h_nJ_n)-\bar b_i(h_nJ_n)
\right)H_{ni}.
$$

This error estimate is used to accept or reject the current step and to determine the next time step.

## Efficient implementation

A key feature of the proposed ExpHR schemes is their low matrix-function implementation cost.

Although ExpHR uses the additional correction

$$
H_{ni}
=
N_n'(U_{ni})F(U_{ni}),
$$

this quantity is inexpensive for a broad class of local nonlinearities.

Since

$$
N_n'(U_{ni})
=
N'(U_{ni})-N'(u_n),
$$

the correction can be written as

$$
H_{ni}
=
\left(
N'(U_{ni})-N'(u_n)
\right)F(U_{ni}).
$$

For many reaction--diffusion problems, \(N'(u)\) is diagonal or acts pointwise. Hence, evaluating \(H_{ni}\) requires only inexpensive vector operations.

No additional nonlinear solve is required.

A central implementation feature of the fifth-order ExpHR methods developed in this work is that the required matrix-function actions can be organized using only two calls to a matrix-function routine per time step.

Thus, despite using additional derivative information, the dominant computational cost remains the evaluation of matrix-function actions.

## Numerical experiments

The numerical experiments assess the accuracy and computational efficiency of the proposed ExpHR methods using both constant and adaptive time stepping.

For constant step sizes, the values of the time step are selected so that the compared methods attain approximately the same error levels. This allows computational cost to be compared at comparable accuracy.

For adaptive time stepping, absolute and relative tolerances are selected so that the methods also attain comparable final-time error levels.

The proposed ExpHR methods are compared with representative exponential Rosenbrock methods, including

- `exprb43`;
- `exprb4s2`;
- `pexprb43`;
- `exprb53s3`;
- `pexprb54s4`;
- `pexprb54s5`.

Representative ExpH methods of comparable orders are also included.

The numerical experiments examine:

1. convergence order;
2. error versus time step;
3. error versus CPU time;
4. constant-step efficiency;
5. adaptive-step efficiency;
6. performance relative to ExpH methods;
7. performance relative to ExpRB methods;
8. behavior for problems with strongly stiff nonlinearities.

Unless otherwise stated, the numerical error is measured at the final time using the discrete maximum norm.

## Test problems

The numerical experiments focus on stiff semilinear parabolic and reaction--diffusion problems in which the nonlinear contribution may itself be strongly stiff.

The test problems include reaction--diffusion systems such as the Schnakenberg model and phase-field-type problems such as the Allen--Cahn equation.

For a general reaction--diffusion problem

$$
u_t=D\Delta u+f(u),
$$

the linear and nonlinear parts are

$$
L=D\Delta,
\qquad
N(u)=f(u),
$$

and the ExpHR Jacobian is

$$
J_n
=
D\Delta+f'(u_n).
$$

Thus, in contrast with methods based only on the fixed operator \(L\), the nonlinear stiffness represented by \(f'(u_n)\) is incorporated directly into the exponential propagation.

For the Allen--Cahn equation

$$
u_t
=
\Delta u
+
\frac{1}{\varepsilon^2}(u-u^3),
$$

we have

$$
L=\Delta,
\qquad
N(u)=\frac{1}{\varepsilon^2}(u-u^3),
$$

and

$$
N'(u)
=
\frac{1}{\varepsilon^2}(1-3u^2).
$$

As \(\varepsilon\) decreases, the nonlinear contribution becomes increasingly stiff, making this problem suitable for testing the performance of ExpHR methods in the strongly nonlinear stiff regime.

## Matrix-function evaluations

Linear combinations of \(\varphi\)-functions applied to vectors are evaluated using the adaptive Krylov-subspace routine `phipm_simul_iom`.

The routine can simultaneously compute expressions of the form

$$
\sum_{k=0}^{q}
\varphi_k(c_\ell hJ_n)V_k,
\qquad
\ell=1,\ldots,i,
$$

for several scaling factors \(c_1,\ldots,c_i\) within a single Krylov process.

This simultaneous evaluation is particularly useful for ExpHR methods because matrix-function actions associated with several stage abscissae can be combined within a single routine call.

Unless otherwise stated, the Krylov tolerance is set to

$$
10^{-10}.
$$

## Constant-step experiments

The constant-step experiments are used to verify the theoretical convergence orders and compare computational efficiency.

For each method, several numbers of time steps are considered. The corresponding final-time errors and CPU times are recorded.

The observed convergence order is computed from successive refinements.

Efficiency is assessed primarily through:

- error versus time step;
- error versus CPU time.

These experiments provide a direct comparison of ExpHR with ExpH and ExpRB methods of the same order.

## Adaptive-step experiments

For adaptive time stepping, the embedded ExpHR pairs are used to estimate the local temporal error.

Both absolute and relative tolerances are prescribed, and the time step is adjusted dynamically using the embedded error estimate.

The adaptive experiments compare:

- achieved final-time error;
- number of accepted time steps;
- evolution of the adaptive time step;
- CPU time;
- efficiency relative to adaptive ExpRB methods.

Particular attention is given to problems with strongly varying nonlinear stiffness.

## Computational environment

All numerical experiments were performed in MATLAB on a MacBook Pro equipped with:

- Apple M1 processor;
- 8 GB of memory.

The implementation and source codes used in the numerical experiments are provided in this repository.

## Dependencies

This repository uses external routines for Krylov-based evaluations of matrix-function actions and, where required, tensor operations.

The routine `phipm_simul_iom` is used to compute linear combinations of \(\varphi\)-functions acting on vectors. This routine is distributed under the BSD 3-Clause License.

Some tensor operations use KronPACK, which is distributed under the MIT License.

The original copyright and license notices of third-party routines are retained in the corresponding source files. See `ThirdPartyNotices.md` for details.

## License

The original code developed in this repository is distributed under the MIT License.

Third-party components remain under their respective licenses.
