# Codes of "A rigorous integrator and global existence for higher-dimensional semilinear parabolic PDEs via semigroup theory"

This repository contains the MATLAB codes associated with the paper:
"A rigorous integrator and global existence for higher-dimensional semilinear parabolic PDEs via semigroup theory"
by G W Duchesne , J-P Lessard and A Takayasu.

**Abstract**  In this paper, we introduce a general constructive method to compute solutions of initial value problems of semilinear parabolic partial differential equations via semigroup theory and computer-assisted proofs (CAPs). Once a numerical candidate for the solution is obtained via a finite dimensional projection, Chebyshev series expansions are used to solve the linearized equations about the approximation from which a solution map operator is constructed. Using the solution operator (which exists from semigroup theory), we define an infinite dimensional contraction operator whose unique fixed point together with its rigorous bounds provide the local inclusion of the solution. Applying this technique for multiple time steps leads to constructive proofs of existence of solutions over long time intervals. As applications, we study the 3D/2D Swift-Hohenberg, where we combine our method with explicit constructions of trapping regions to prove global existence of solutions of initial value problems converging asymptotically to nontrivial equilibria. A second application consists of the 2D Ohta-Kawasaki equation, providing a framework for handling derivatives in nonlinear terms.

These codes require *MATLAB* with [*INTLAB* - INTerval LABoratory](http://www.ti3.tu-harburg.de/rump/intlab/) (MATLAB toolbox for interval arithmetic) version 11.

---

A rough correspondence for some of the files & computational procedures in the paper are as follows:

### CAPs for global existence of solutions to 3D Swift-Hohenberg equation (sec 5.3.1)

```
>> cd 3D-SH
>> script_proof_GE_3DSH
```

### CAPs for global existence of solutions to 2D Swift-Hohenberg equation (sec 5.3.2)

Stripe pattern equilibrium:

```
>> cd ../2D-SH/
>> script_proof_GE1_2DSH
```

Spot pattern equilibrium:

```
>> script_proof_GE2_2DSH
```

Figs 2-4 are plotted by

```
>> script_plot_equilibria
```

### Rigorous integration of solution to 2D Ohta-Kawasaki equation (sec 6.1)

Stripe pattern state:

```
>> cd ../2D-OK/
>> script_integrate_2DOK_stripe
```

Spot pattern state:

```
>> script_integrate_2DOK_spot
```

Figs 5 & 6 are plotted by

```
>> script_plot_pattern_2DOK % Fig 5
>> script_plot_data_2DOK % Fig 6
```


Copyright (C) 2024  G W Duchesne, J-P Lessard and A Takayasu.
