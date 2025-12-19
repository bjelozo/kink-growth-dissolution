# kink-growth-dissolution
Collective variables for Plumed and Gromacs to simulate kink growth and dissolution of molecules and ions

The code consists of following collective variables (CVs):
- Surface structure CV: ZrowC8C9GXY.cpp
- Gaussian biased CV kernel: ZsphereG2.cpp
- Wall CVs for growth prevention: ZpyXY.cpp

If you use these CVs, please cite following papers:
- Solubility prediction of organic molecules with molecular dynamics simulations; Z Bjelobrk, D Mendels, T Karmakar, M Parrinello, M Mazzotti; Crystal Growth & Design 21 (9), 5198-5205. https://doi.org/10.1021/acs.cgd.1c00546
- Solubility of Organic Salts in Solvent–Antisolvent Mixtures: A Combined Experimental and Molecular Dynamics Simulations Approach; Z Bjelobrk, AK Rajagopalan, D Mendels, T Karmakar, M Parrinello, M. Mazzotti; Journal of Chemical Theory and Computation 18 (8), 4952-4959 https://doi.org/10.1021/acs.jctc.2c00304

Please cite Plumed as described here: https://www.plumed.org/cite;
Also cite any other methodology such as the constant chemical potential algorithm or well-tempered Metadynamics if used in combination with the presented CVs:
- Molecular dynamics simulations of solutions at constant chemical potential; C Perego, M Salvalaglio, M Parrinello; The Journal of chemical physics 142 (14). https://doi.org/10.1063/1.4917200
- Well-tempered metadynamics: a smoothly converging and tunable free-energy method; A Barducci, G Bussi, M Parrinello; Physical review letters 100 (2), 020603. https://doi.org/10.1103/PhysRevLett.100.020603
