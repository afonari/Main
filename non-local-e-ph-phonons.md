## Nonlocal electron-phonon coupling
Nonlocal el-ph coupling results from the interaction of the charge with lattice vibration (phonon). As decribed in literature [1] el-ph coupling can be related to the changes in couplings of the frontier orbitals of neighbouring molecules due to changes in geometry of the molecules. In particular, coupling constant that depends on the derivative of the transfer integral (site energy is another option) over normal mode displacement can be defined:  
![g_{i}](https://raw.github.com/alexandr-fonari/Main/master/pics/gi.png)  
where, *Qi* is dimensionless displacement, *t* is transfer integral computed for a pair of neighbouring molecules (units of energy), thus *gi* has units of energy.

Also, small polaron binding energy is defined as:  
![L](https://raw.github.com/alexandr-fonari/Main/master/pics/L.png)  
*L* has units of energy.

Because method of the finite differences will be used to obtain *gi*, for reproducibility purposes, step-size should be defined in a unique way. A natural choice is normal mode characteristic length defined as [2]:  
![q0](https://raw.github.com/alexandr-fonari/Main/master/pics/q0.png)  
*q0* has units of length.

#### Acknowledgments and references
1. V. Coropceanu, *et al.*, *Chem. Rev.*, **107**, 926 (2007): [10.1021/cr050140x](http://pubs.acs.org/doi/abs/10.1021/cr050140x).
1. [Quantum harmonic oscillator](http://en.wikiversity.org/wiki/Quantum_harmonic_oscillator) @ Wikiversity.