# smallcodes_ALAMODE


## Introduction

`smallcodes_ALAMODE` includes several usefull python scripts for ALAMODE, which is anharmonic phonon calculator. Below are the list of codes.

- plot_ifc.py : plotting atomic distances vs IFCs from ALAMODE ifc files(.fcs) to check if IFCs are converged with respect to atomic distances.
- calc_pes.py : calculating potential energy surface from ALAMODE ifc files(.xml).
- plot_ltc.py : plotting Lattice Thermal Conductivity from ALAMODE *.kl files.

## Requirements

Following python modules are required.

- ase
- xml
- matplotlib

## Examples (see also sample/)

### plot_ifc.py

![ifc_2th](image/2th_order_ifc.jpeg)

![ifc_6th](image/6th_order_ifc.jpeg)


### calc_pes.py

`calc_pes.py` calculates the potential energy from an anharmonic phonon model for given atomic displacements. This calculation is very important when you want to check the validity of your model compared to the DFT energies.

The typical workflow is as follows. After making the anharmonic phonon model of `*.xml` from , you prepare both equilibrium and displaced atomic configurations in VASP POSCAR format for your primitive cell. Then execute the code with

```bash
./calc_pes.py -p <POSCAR_equilibrium> -i <anharm.xml> -T <temperature> -f <maxorder> POSCAR_disp1 POSCAR_disp2 ...
```

Here, `<POSCAR_equilibrium>` contains the equilibrium positions, `<anharm.xml>` is your model file, `<POSCAR_disp>` are displaced atomic configurations. `<maxorder>` specifies the maximum order you want to calculate, which is useful when you only want to calculate harmonic energy (`maxorder=2`) or quartic energy (`maxorder=4`) even though your model contains up to 6-th order IFCs.

In practice, you want to calculate the potential energy surface for a certain phonon mode at the Gamma point. In this case, `` in the ALAMODE is useful to generate this kind of atomic coordinates.

Examples are given in the `sample/calc_pes` directory.

### plot_ltc.py

