# Salve como 'plot_populacoes.gp'
# Para rodar no terminal: gnuplot plot_populacoes.gp

# Configuração da saída de imagem
set terminal pngcairo size 1000,700 enhanced font 'Helvetica,14'
set output 'populacoes_particulas.png'

# Títulos e Eixos
set title "Frações das Partículas na Equação de Estado (EoS)"
set xlabel "Densidade Bariônica Relativa ({/Symbol r}_B / {/Symbol r}_0)"
set ylabel "Fração de Partícula (Y_i = n_i / n_B)"

# Escala Logarítmica para ver partículas que surgem depois
set logscale y
set yrange [1e-5:2]
set format y "10^{%T}"
set xrange [0:*]

# Grades visuais
set grid ytics lc rgb "#e0e0e0"
set grid xtics lc rgb "#e0e0e0"

# Posição da legenda
set key outside right top spacing 1.5

# A coluna 1 é nb/n0. Para obtermos a densidade real nb (para dividir as frações),
# multiplicamos pela densidade de saturação nuclear
n0 = 0.153

# Plotagem
# O Gnuplot automaticamente ignorará avisos de "log(0)" nas baixas densidades, 
# traçando as linhas exatamente a partir de onde a partícula surge.
plot 'eos.dat' u 1:($6/($1/n0)) w l lw 2 lc rgb "black"      t 'n', \
     ''        u 1:($7/($1/n0)) w l lw 2 lc rgb "red"        t 'p', \
     ''        u 1:($4/($1/n0)) w l lw 2 lc rgb "blue"       t 'e^-', \
     ''        u 1:($5/($1/n0)) w l lw 2 lc rgb "cyan"       t '{/Symbol m}^-', \
     ''        u 1:($8/($1/n0)) w l lw 2 dt 2 lc rgb "green" t '{/Symbol L}^0', \
     ''        u 1:($9/($1/n0)) w l lw 2 dt 2 lc rgb "orange" t '{/Symbol S}^-', \
     ''        u 1:($10/($1/n0)) w l lw 2 dt 2 lc rgb "brown" t '{/Symbol S}^0', \
     ''        u 1:($11/($1/n0)) w l lw 2 dt 2 lc rgb "pink"  t '{/Symbol S}^+', \
     ''        u 1:($12/($1/n0)) w l lw 2 dt 4 lc rgb "purple" t '{/Symbol X}^-', \
     ''        u 1:($13/($1/n0)) w l lw 2 dt 4 lc rgb "dark-magenta" t '{/Symbol X}^0'