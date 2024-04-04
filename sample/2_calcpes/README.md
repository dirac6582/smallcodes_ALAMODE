# Examples

The purpose of the example is to reproduce Fig. 4 in 10.1103/PhysRevB.107.094305, where we calculated the anharmonic phonon energy of the rutile TiO2.


## Necessary files

- : Only need for generating displacement patterns
- POSCAR: equilibrium positions of the rutile TiO2
- anharm.xml: anharmonic force constants generated from `alm`


## Step1: generate displacement patterns

You can generate the displacement patterns using in the ALAMODE as 

```bash

```

In this case, we move atoms along the A2u phonon, which is the most anharmonic mode in the rutile TiO2.

The calculated DFT energies are given in .


## Step2: calculate anharmonic phonon energies

Using the displacement paterns made in Step1, we calculate the potential energy as 

```bash
python plot_pes_to_tadano.py -i TiO2224_anharm.xml -p disp_POSCAR/A2u_0_4/POSCAR disp_POSCAR/A2u_0_4/disp{01..10}.POSCAR
```

The results are collected in `calc_pes_result.txt`,

```bash
# max displacement in Ang and energy in eV.
  0.00000000   0.00000000    0.00000000    0.00000000    0.00000000
  0.03850026   0.00156511    0.00156511    0.00161142    0.00161142
  0.07700052   0.00626046    0.00626046    0.00700138    0.00700138
```

where the first column represents the maximum displacement of the given configurations, followed by the potential energy up to the second to 6-th order in eV unit.



## Step3: Make a figure to compare anharmonic energies with DFT ones

To see if the phonon model works fine, we compare the potential energy from the phonon model and DFT ones.
