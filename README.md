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

Dynamic Mode Decomposition approximates a linear operator that best advances the system from one snapshot to the next. For a sequence of snapshots \(\mathbf{x}_1, \mathbf{x}_2, \ldots, \mathbf{x}_m\), DMD seeks an operator \(\mathbf{A}\) such that:


The eigendecomposition of \(\mathbf{A}\) provides the DMD modes and eigenvalues, which characterize the temporal dynamics and spatial structures of the system.
