set title "Resultats de MC2 en funció de la temperatura"
set xlabel "T"
set grid
set terminal png size 1700,2000 enhanced font "Helvetica,20"
set output 'MC2-grafic.png'

plot 'mc2-resultats.res' using 1:2 with lines title "<E/N>",'' using 1:3 with lines title "<E²|N>",'' using 1:4 with lines title "<M|N>",''using 1:5 with lines title "<M²|N>",'' using 1:6 with lines title "<|M|>",'' using 1:7 with lines title "VAR(E)",'' using 1:8 with lines title "VAR(M)"
