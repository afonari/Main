# VASP implementation
Currently, phonons can be obtained in a straightforward way only at Γ-point (there are other, *less* straight ways, to access phonons at different k-points)[1].

## vasprun.xml
As the result of a phonon calculation (```IBRION = 5/6/7/8```), ```vasprun.xml``` will contain hessian matrix.

#### Acknowledgments and references
1. ```IBRION``` tag in [VASP Manual](http://cms.mpi.univie.ac.at/vasp/vasp/IBRION_tag_NFREE_tag.html).
1. [Quantum harmonic oscillator](http://en.wikiversity.org/wiki/Quantum_harmonic_oscillator) @ Wikiversity.