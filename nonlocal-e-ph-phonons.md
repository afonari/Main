## Nonlocal electron-phonon coupling through phonons calculation
Nonlocal el-ph coupling results from the interaction of the charge with lattice vibration (phonon). As decribed in literature [1] el-ph coupling can be related to the changes in couplings of the frontier orbitals of neighbouring molecules due to changes in geometry of the molecules. In particular, coupling constant that depends on the derivative of the transfer integral (site energy is another option) over normal mode displacement can be defined:  
![g_{i}](https://raw.github.com/alexandr-fonari/Main/master/pics/gi.png)  
where, *Qi* is dimensionless displacement, *t* is transfer integral computed for a pair of neighbouring molecules (units of energy), thus *gi* has units of energy.

Also, small polaron binding energy is defined as:  
![L](https://raw.github.com/alexandr-fonari/Main/master/pics/L.png)  
*L* has units of energy.

Derivative in (1) can be evaluated numerically, (using finite difference method), thus for reproducibility purposes, step-size should be defined in a unique way. A natural choice is normal mode characteristic length defined as [2]:  
![q0](https://raw.github.com/alexandr-fonari/Main/master/pics/q0.png)  
*q0* has units of length.

#### Example
For a choosen normal (e.g. ```ħω = 0.28069053 eV```) four displacements along a normal mode are made ```q0 = {-2, -1, 1, 2}```. The plot of ground state energy (SCF calculation) with respect to displacement looks like this:  
![Energy vs q0](https://raw.github.com/alexandr-fonari/Main/master/pics/e_vs_q0.png)  
If all calculations are done correctly, than  
```ħω/2 = k``` (4)  
where *k* is curvature of the parabola ```y=kx²/2```. Comparing left side (0.14034 eV) with the right side (0.14108 eV) of (4) it could be concluded that both are approximately equal.

#### Acknowledgments and references
1. V. Coropceanu, *et al.*, *Chem. Rev.*, **107**, 926 (2007): [10.1021/cr050140x](http://dx.doi.org/10.1021/cr050140x).
1. [Quantum harmonic oscillator](http://en.wikiversity.org/wiki/Quantum_harmonic_oscillator) @ Wikiversity.