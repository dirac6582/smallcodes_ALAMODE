# Examples

The purpose of the example is to reproduce Fig. 4 in 10.1103/PhysRevB.107.094305, where we calculated the anharmonic phonon energy of the rutile TiO2.


## Necessary files

- TiO2224_gamma.evec: Only need for generating displacement patterns
- POSCAR: equilibrium positions of the rutile TiO2
- anharm.xml: anharmonic force constants generated from `alm`


## Step1: generate displacement patterns

You can generate the displacement patterns using `displace.py` in the ALAMODE as 

```bash
displace.py --evec TiO2224_gamma.evec --VASP POSCAR --pes='0 6' --Qrange "0 0.5" --num_disp 41
```

In this case, we move atoms along the A2u phonon, which is the most anharmonic mode in the rutile TiO2. `--pes` specifies the 1st Q-point and 6th branch given in the `evec` file. `Qrange` gives the minimum and maximum normal coordinate amplitude in units of `amu^{1/2}*Angstrom` with the number of displacement being given by `--num_disp` option.  For details, please refer to `alamode`. The resultant `disp*.POSCAR` files and the calculated DFT energies are given in the `output/` directory.


## Step2: calculate anharmonic phonon energies

Using the displacement paterns made in Step1, we calculate the potential energy of the anharmonic phonon model as 

```bash
python ../../plot_pes_to_tadano.py -i TiO2224_anharm.xml -p POSCAR output/disp{01..41}.POSCAR
```

The results are collected in `calc_pes_result.txt`,

```bash
# max displacement in Ang and energy in eV.
  0.00000000   0.00000000    0.00000000    0.00000000    0.00000000
  0.00457191   0.00002207    0.00002207    0.00002208    0.00002208
  0.00914381   0.00008828    0.00008828    0.00008843    0.00008843
  0.01371572   0.00019864    0.00019864    0.00019938    0.00019938
```

where the first column represents the maximum displacement of the given configurations, followed by the potential energy up to the second to 6-th order in eV unit. The sample output is in the `output/` directory.


## Step3: Make a figure to compare anharmonic energies with DFT ones

To see if the phonon model works fine, we compare the potential energy from the phonon model and DFT ones. Here is a simple gnuplot script to make a graph.

```gnuplot
# output setting
set term pdfcairo enhanced # enhanced for latin
set output "rutileTiO2_A2u_pes.pdf"


# title
set title "A2u"

# axis label
set xlabel "Displacement [Ang]" offset 0,-1
set ylabel "Energy [meV]" offset -2,0
# set y2label "y2"

# margin setting (if num is large, margin becomes large)
set lmargin 15 #left
set bmargin 5 #bottom


set xrange [0:0.1]
set yrange [0:40]
set key left top

# plot
plot \
  "output/vasp20230117.txt"    u 1:(($2+69.82869692)*1000) title "vasp" ,\
  "calc_pes_result.txt" u ($1*0.529177249):($2*1000) title "harm",\
  "calc_pes_result.txt" u ($1*0.529177249):($4*1000) title "quartic"
```

We plot the harmonic and quartic potential energy compared to the VASP values. Be sure to use the relative potential energy from the equilibrium configuration for DFT results. The anharmonic phonon model only calculates the energy differences from the equilibrium. The sample plot is given in `output/rutileTiO2_A2u_pes.pdf`.