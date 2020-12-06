## LIF network with white Gaussian noise

* This code solves the equation for average firing-rate of a LIF neuron
 driven by white Gaussian noise as given by Brunel (2000) (and others).

* Assumptions: Network has sparse random connectivity (C << N )
Thus, correlations of the fluctuating part of the synaptic inputs 
of different neurons are neglected in the limit C/N -> 0.

* How it solves: fsolve tries to find the roots from Eq.(21) from Brunel
(2000) iteratively. The conditions such as number of iterations or step 
are set by optimset.

by: rodrigo pena
pena@njit.edu / rfdop20@gmail.com                             
