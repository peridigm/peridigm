
set terminal pdf enhanced font "Times,32" size 12in, 8in
set output "initial_displacement.pdf"
set title "Initial Displacement" font "Times,42"
set xlabel "r" font "Times,32"
set ylabel "Displacement" font "Times,32"
plot "WaveInCube.txt" using 1:2 with points pt 7 ps 1 lc 2 notitle

set terminal pdf enhanced font "Times,32" size 12in, 8in
set output "final_displacement.pdf"
set title "Final Displacement" font "Times,42"
set xlabel "r" font "Times,32"
set ylabel "Displacement" font "Times,32"
plot "WaveInCube.txt" using 1:3 with points pt 7 ps 1 lc 2 notitle
