# ================================================
# Configuração Global
# ================================================
set datafile separator ','

# --- DEFINA SEUS PARÂMETROS AQUI ---
MODELS = "GM1 GM3"
# Coloque aqui os 3 valores exatos de B como aparecem no seu CSV (após o Rust converter 'd' para 'e')
B_VALS = "1.0e15 1.0e16 1.0e17" 

# Função segura de log10 numérico
safe_log(x) = (x > 0) ? log10(x) : NaN

# Estilo Visual
set style line 101 lc rgb '#e0e0e0' lt 1 lw 1
set grid back ls 101
set border 3 front lw 1.5
set tics nomirror

COLOR_MASS = "#1f77b4"
COLOR_RADIUS = "#d62728"

# ================================================
# Loop Externo: Uma Figura 2x2 para cada Campo B
# ================================================
do for [b_val in B_VALS] {

    # Gera um arquivo diferente para cada B (ex: grid_mr_B_1.0e15.png)
    set terminal pngcairo size 1200,1000 enhanced font 'Arial,12'
    set output sprintf('results/grid_mr_B_%s.png', b_val)
    
    # Adiciona o valor de B no título principal da imagem
    set multiplot layout 2,2 rowsfirst \
        title sprintf("Análise de Parâmetros: GM1 e GM3 (B = %s G)", b_val) font ",18" offset 0,-0.5

    # ------------------------------------------------
    # 1ª LINHA (SUPERIOR): APENAS MASSA
    # ------------------------------------------------
    do for [model in MODELS] {
        set title model font ",14"
        set xlabel "log_{10}({/Symbol x})"
        set ylabel "Massa Máxima [M_{/Symbol \304}]" textcolor rgb COLOR_MASS
        set key bottom left box opaque font ",10"

        # O filtro agora exige que col 1 == model E col 2 == b_val
        plot 'results/summary.csv' \
             using (strcol(1) eq model && strcol(2) eq b_val ? safe_log($3) : NaN):4 \
             with points pt 7 ps 1.2 lc rgb COLOR_MASS title 'Massa Máxima'
    }

    # ------------------------------------------------
    # 2ª LINHA (INFERIOR): APENAS RAIO
    # ------------------------------------------------
    do for [model in MODELS] {
        set title ""
        set xlabel "log_{10}({/Symbol x})"
        set ylabel "Raio Máximo [km]" textcolor rgb COLOR_RADIUS
        set key bottom left box opaque font ",10"

        plot 'results/summary.csv' \
             using (strcol(1) eq model && strcol(2) eq b_val ? safe_log($3) : NaN):5 \
             with points pt 7 ps 1.2 lc rgb COLOR_RADIUS title 'Raio Máximo'
    }

    unset multiplot
}