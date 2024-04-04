# ========================
# %file%
#
# ========================
# Time-stamp: <2024-04-05 01:04:10 amano>


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

