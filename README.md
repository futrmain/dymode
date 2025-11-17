# Dymode

A parallel dynamic mode decomposition software.

## Overview

Dynamic Mode Decomposition (DMD) is an analysis technique that decomposes complex, nonlinear systems into summable modes, 
revealing patterns and dynamics that characterize the system's behavior. 
It provides an efficient framework for extracting spatial-temporal coherent modes from complex flow data.
This reposirory contains an implementation of DMD for Matlab, as well as a C++ version that leverages 
parallel computing with MPI to enable analysis of massive datasets that would be impractical with serial implementations.
Related work and official documentation is availbale at [Dymode: A parallel dynamic mode decomposition software](https://kth.diva-portal.org/smash/record.jsf?pid=diva2%3A786623&dswid=3187).

