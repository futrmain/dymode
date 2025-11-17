# Dymode

A parallel dynamic mode decomposition software.

## Overview

Dynamic Mode Decomposition (DMD) is an analysis technique that decomposes nonlinear systems into 
summable spatio-temporal modes, 
revealing patterns and dynamics that characterize the system's behavior. 

This reposirory contains an implementation of DMD for Matlab, as well as a C++ version that leverages 
parallel computing with MPI to enable analysis of massive datasets.
Related work and official documentation is availbale at [Dymode: A parallel dynamic mode decomposition software](https://kth.diva-portal.org/smash/record.jsf?pid=diva2%3A786623&dswid=3187).

## Theory

DMD is based on the Koopman operator $A$ that relates the state of a system from one time step to the next:

$$ \mathbf{x}_{k+1} \approx \mathbf{A}\mathbf{x}_k $$

Even though $A$ is a _linear_ operator, it is able to describe _non-linear_ systems when working in infinite dimensions. 
The eigendecomposition of $A$ then provides modes that are
1. tied to a unique frequency
2. coherent in space and time
3. may grow or decay over time
Points 2 and 3 make the decomposition different and more valuable than standard Fourier based decompositions.

Dynamic Mode Decomposition works by approximating the eigendecomposition of $A$ by that of a finite-dimensional operator which best advances the observed state of a system from one snapshot to the next. 
This can be done in several ways, but `dymode` uses singular value decomposition for robustness and accuracy. 

We seek to find the eigendecomposition of the operator $A$ that satisfies

$$A\mathbf{X}_{k} = \mathbf{X}_{k+1}$$

where $\mathbf{X}_{k}$ is a matrix of snapshots from time-step 0 to $k$, and $\mathbf{X}_{k+1}$ from 1 to $k+1$.
Using the singular value decomposition of $\mathbf{X}_{k} = U\Sigma V^{t}$, 

\begin{aligned}
A\mathbf{X}_{k} & = \mathbf{X}_{k+1} \
AU\Sigma V^{t} & = \mathbf{X}_{k+1}
\end{aligned}


The eigendecomposition of \(\mathbf{A}\) provides the DMD modes and eigenvalues, which characterize the temporal dynamics and spatial structures of the system.
