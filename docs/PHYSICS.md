# Física e Modelos

Este arquivo reúne os modelos físicos e as equações utilizadas no NSRS.

## Matéria Bariônica (Teoria RMF)

A matéria bariônica é descrita pela Hadrodinâmica Quântica (QHD), dentro da aproximação de Campo Médio Relativístico (RMF), onde os graus de liberdade fundamentais são os férmions que compõem o octeto completo de bárions, $b = (n, p, \Lambda^0, \Sigma^-, \Sigma^0, \Sigma^+, \Xi^-, \Xi^0)$, e os léptons, $\ell = (e^-, \mu^-)$, necessários para satisfazer o equilíbrio beta e a neutralidade de carga. Estes bárions interagem via acoplamento mínimo com campos mesônicos clássicos: o escalar $\sigma$, o vetorial $\omega^\mu$ e o isovetorial $\rho^\mu$.

A Lagrangiana da QHD é construída como a soma das Lagrangianas individuais para bárions, mésons, léptons e o campo eletromagnético:

$$
\mathcal{L}_{QHD} = \mathcal{L}_{\text{baryons}} + \mathcal{L}_{\text{mesons}} + \mathcal{L}_{\text{leptons}} + \mathcal{L}_{\text{EM}}
$$

### Setor Bariônico

As contribuições bariônicas ditam a dinâmica do octeto completo de bárions acoplado aos campos mesônicos e eletromagnético. A densidade Lagrangiana é dada por:

$$
\mathcal{L}_{\text{baryons}} = \sum_{b} \bar{\psi}_b \left[ \gamma_\mu \left(i \partial^\mu + e_b A^\mu - g_{\omega,b} \omega^\mu - g_{\rho,b} I_{3} \rho^\mu\right) - M^*_b \right] \psi_b
$$

onde $\psi_b$ e $\bar{\psi_b}$ representam o spinor de Dirac e seu adjunto para um dado bárion $b$, enquanto $\gamma_\mu$ denota as matrizes de Dirac padrão e $A^\mu$ é o potencial vetor eletromagnético. A carga elétrica do bárion é dada por $e_b$, e $I_3$ representa a terceira componente do operador de isospin para cada bárion específico. A massa efetiva do bárion, $M_{b}^{\ast}$, é dinamicamente modificada pelo campo escalar $\sigma$ e definida como $M_{b}^{\ast} \equiv m_{b} - g_{\sigma,b}\sigma$, onde $m_{b}$ representa a massa nua do bárion. Finalmente, os parâmetros $g_{\sigma,b}$, $g_{\omega,b}$ e $g_{\rho,b}$ denotam as constantes de acoplamento específicas do bárion com os campos escalar ($\sigma$), vetorial ($\omega$) e isovetorial ($\rho$), respectivamente.

### Setor Mesônico

O setor mesônico abrange os termos cinéticos e de massa dos campos mesônicos livres, bem como as auto-interações não-lineares do campo escalar parametrizadas por $\kappa$ e $\lambda$. Estes termos cúbicos e quárticos são estritamente necessários para reproduzir as propriedades empíricas da matéria nuclear simétrica na densidade de saturação dentro do modelo de Walecka estendido. A densidade de Lagrangiana mesônica é dada por:

$$
\begin{aligned}
    \mathcal{L}_{\text{mesons}} &= \frac{1}{2}\partial_{\mu}\sigma\partial^{\mu}\sigma - \frac{1}{2}m_{\sigma}^{2}\sigma^{2} - \frac{1}{4}\Omega_{\mu\nu}\Omega^{\mu\nu} + \frac{1}{2}m_{\omega}^{2}\omega_{\mu}\omega^{\mu} \\
    &\quad - \frac{1}{4} \mathbf{P}_{\mu\nu} \cdot \mathbf{P}^{\mu\nu} + \frac{1}{2}m_{\rho}^{2} \boldsymbol{\rho}_{\mu} \cdot \boldsymbol{\rho}^{\mu} - \frac{1}{3!}\kappa\sigma^{3} - \frac{1}{4!}\lambda\sigma^{4}
\end{aligned}
$$

onde o tensor mesônico simétrico é definido como:

$$
\Omega^{\mu \nu} = \partial^\mu \omega^\nu - \partial^\nu \omega^\mu
$$

e a força de campo $SU(2)$ é definida como:

$$
\mathbf{P}^{\mu \nu} = \partial^\mu \boldsymbol{\rho}^\nu - \partial^\nu \boldsymbol{\rho}^\mu + g_{\rho} \boldsymbol{\rho}^\mu \times \boldsymbol{\rho}^\nu
$$

### Setor Leptônico

Para garantir o equilíbrio químico e a neutralidade de carga global dentro da matéria da estrela de nêutrons, um setor léptônico deve ser incluído. Os léptons, especificamente elétrons ($e^-$) e múons ($\mu^-$), não participam da interação forte e são acoplados apenas ao campo eletromagnético. Sua dinâmica é governada pela seguinte densidade de Lagrangiana:

$$
\mathcal{L}_{\text{leptons}} = \sum_{\ell} \bar{\psi}_\ell \left[  \gamma^\mu (i \partial_\mu + e_\ell A_\mu) - m_\ell\right]\psi_\ell
$$

onde $\psi_\ell$ e $\bar{\psi}_\ell$ representam o spinor de Dirac e seu adjunto para o lépton $\ell$, e $m_\ell$ é a respectiva massa do lépton. O termo $e_\ell A_\mu$ representa o acoplamento mínimo dos léptons com o campo eletromagnético, onde $e_\ell$ é a carga elétrica do lépton. Além desta interação eletromagnética, os léptons são tratados como um gás de Fermi livre.

### Campo Eletromagnético

A dinâmica do próprio campo eletromagnético livre é descrita pela Lagrangiana de Maxwell padrão:

$$
\mathcal{L}_{\text{EM}} = -\frac{1}{4}F_{\mu\nu}F^{\mu\nu}
$$

Aqui, $F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu$ representa o tensor de força do campo eletromagnético de Maxwell padrão, e $A_\mu$ é o campo de fótons visíveis acoplado à corrente eletromagnética $J^\mu_{\text{EM}}$.

## Campo magnético: perfis e tensões

Todas as intensidades de campo são em Gauss (inclusive o parâmetro $\xi$ da
eletrodinâmica logarítmica). O código está em `src/core/magnetic.rs`.

### Perfis do campo local (`FieldProfile`)

| Perfil | Níveis de Landau | Energia magnética | Referência |
|---|---|---|---|
| `Constant` (padrão) | $B$ central `bg` em todas as densidades | BDD com $B_{\rm surf}=10^{15}$ G, $B_0=$ `bg` | comportamento legado |
| `Bdd` | $B(n_B)$ local | mesmo $B(n_B)$ | Bandyopadhyay, Chakrabarty & Pal, PRL 79, 2176 (1997) |
| `Dexheimer2017` | $B(\mu_B)$ local | não entra (padrão) | Dexheimer et al., PLB 773, 487 (2017) |

**BDD.** $B(n_B)=B_{\rm surf}+B_0\left[1-e^{-\beta(n_B/n_0)^\gamma}\right]$, com
$\beta=0.01$, $\gamma=3$ e $B_{\rm surf}=10^{15}$ G (intensidade máxima de superfície
observada em magnetares, adotada nos trabalhos recentes). Como o solver avança em
$\mu_n$, cada ponto usa uma iteração de ponto fixo $B\leftrightarrow n_B$ até
$|\Delta n_B/n_B|<10^{-10}$.

**Dexheimer et al. (2017), Eq. (1).** Ajuste a soluções de Einstein–Maxwell
(direção polar, campo poloidal):

$$
B(\mu_B)=\frac{(a+b\,\mu_B+c\,\mu_B^2)\,\mu}{B_c},\qquad B_c=4.414\times10^{13}\ {\rm G},
$$

com $\mu_B$ em MeV ($\mu_B=\mu_n$ em equilíbrio β), $\mu$ o momento de dipolo em
A m² e $(a,b,c)$ da Tabela 2 para $M_B=2.2\,M_\odot$ ou $1.6\,M_\odot$. O ajuste
cobre $\mu_B\approx 939$–$1500$ MeV. Abaixo de $m_N$ usa-se $B(m_N)$. Acima de
1500 MeV o polinômio é extrapolado até o vértice ($-b/2c$: 1734 MeV para
$M_B=2.2$, 1629 MeV para $M_B=1.6$) e congelado a partir dele. A extrapolação é
necessária: o centro da estrela de massa máxima tem $\mu_n\approx1584$ MeV (GM1) e
1525 MeV (GM3). Entre 1500 e 1600 MeV o campo cresce ~5%.

Limitações: o perfil é o da direção polar de uma estrela de massa bariônica e
dipolo fixos, e é usado aqui numa TOV esférica e isotrópica.

**Energia e tensões do campo na TOV** (`with_field_stress`). Por padrão entram para
`Constant` e `Bdd` (prática da literatura com o perfil BDD) e não entram para
`Dexheimer2017`: o ajuste vem de soluções de Einstein–Maxwell, em que o campo é
tratado na estrutura, e é destinado à EoS microscópica. Como $B(m_N)\approx4\times10^{17}$ G
para $\mu=3\times10^{32}$ A m², somar $B^2/8\pi\approx3$ MeV/fm³ à EoS criaria um envelope
sem matéria ($M_{\max}=4.4\,M_\odot$, $R_{1.4}=38$ km para GM1). Sem as tensões, o campo
na matéria muda $M_{\max}$ em −0.2% (GM1) e −0.1% (GM3). Note que a NLEM só atua pelas
tensões do campo; para estudá-la com este perfil use `with_field_stress(true)`.

**Aproximação termodinâmica.** Em cada ponto a EoS é resolvida com o campo local
como parâmetro externo, como nas duas referências. Com $B=B(\mu_B)$,

$$
\frac{dP_m}{d\mu_n}=n_B+\mathcal M\,\frac{dB}{d\mu_n},\qquad
\mathcal M=\left.\frac{\partial P_m}{\partial B}\right|_{\mu},
$$

e o termo de magnetização não é incluído em $n_B$. Para GM1 com
$\mu=3\times10^{32}$ A m² ele vale $5\times10^{-4}$–$2.4\times10^{-3}\,n_B$ entre
$n_0$ e $6n_0$ (verificado em `tests/magnetic_profiles.rs`). Para o perfil BDD, os
termos em $dB/dn_B$ desprezados valem $\lesssim10^{-3}$ até $B_0=10^{18}$ G e
~2% em $B_0=5\times10^{18}$ G.

### Momentos magnéticos anômalos (AMM)

Desligados por padrão; `HadronsMatter::with_anomalous_moments(true)` (ou `--amm` em
`nsrs study b`) acrescenta o termo de Pauli
$\tfrac12\kappa_b\mu_N\bar\psi_b\sigma_{\mu\nu}F^{\mu\nu}\psi_b$, com
$\kappa_b=\mu_b/\mu_N-q_b\,m_p/m_b$ (PDG; `KAPPA_B` em `constants.rs`) e
$\mu_N=e/2m_p$. Com $s=\pm1$ e $a_s=s\,\kappa_b\mu_N B$:

- bárions carregados: $E=\sqrt{k_z^2+\big(\sqrt{M^{*2}+2\nu|q_b|B}-a_s\big)^2}$, e
  $\nu=0$ tem um único estado de spin;
- bárions neutros: $E=\sqrt{k_z^2+\big(\sqrt{M^{*2}+k_\perp^2}-a_s\big)^2}$.

Para os neutros, com $\bar m=M^*-a_s$, $k_F=\sqrt{E_F^2-\bar m^2}$,
$A=\arcsin(\bar m/E_F)-\pi/2$ e $L=\ln[(E_F+k_F)/|\bar m|]$, cada estado de spin dá

$$
n=\frac{1}{2\pi^2}\Big[\frac{k_F^3}{3}-\frac{a_s}{2}\big(\bar m k_F+E_F^2A\big)\Big],\qquad
n_s=\frac{M^*}{4\pi^2}\big[E_Fk_F-\bar m^2L\big],
$$

$$
\epsilon=\frac{1}{4\pi^2}\Big[\frac{E_F^3k_F}{2}-\frac{\bar m}{4}\big(\bar mk_FE_F+\bar m^3L\big)
-\frac{a_s}{3}\big(E_F\bar mk_F+\bar m^3L\big)-\frac{2}{3}a_sE_F^3A\Big]
$$

(cf. Broderick, Prakash & Lattimer, ApJ 537, 351 (2000)). As formas fechadas foram
derivadas e conferidas contra quadratura numérica, $dP/dE_F=n$ e
$\partial(\epsilon-E_Fn)/\partial M^*=n_s$ (testes em `particles.rs`); com $a_s=0$
reduzem-se ao gás isotrópico. Com AMM a $10^{18}$ G, $dP/d\mu_n=n_B$ continua valendo
(teste `anomalous_moments_keep_pressure_consistent`). Para o nêutron
($\kappa_n=-1.913$) o desdobramento é $|a_s|\approx6$ MeV a $10^{18}$ G; o efeito na
estrutura só aparece acima de ~$10^{18}$ G (GM1, campo constante: $R_{1.4}$ cai de 13.01
para 12.64 km a $3\times10^{18}$ G). O valor de $\kappa_{\Sigma^0}$ é a média de
$\Sigma^\pm$ (não medido); léptons ficam sem AMM.

### Magnetização e pressão anisotrópica da matéria

A pressão termodinâmica da matéria, $P_\parallel=-\Omega$, é a pressão ao longo do
campo. Perpendicularmente às linhas de campo,

$$
P_\perp=P_\parallel-\mathcal M B,\qquad
\mathcal M B=B\left.\frac{\partial P_\parallel}{\partial B}\right|_{\mu}
$$

(Ferrer et al., PRC 82, 065802 (2010); Strickland, Dexheimer & Menezes, PRD 86,
125032 (2012)). $\mathcal MB$ é calculado em cada ponto por diferença central com
duas soluções em $B(1\pm10^{-5})$ a $\mu_n$ fixo (teste
`magnetization_is_the_field_derivative_of_the_parallel_pressure`). A coluna 2 da EoS
exporta a pressão usada na TOV: topologia anisotrópica,
$P_\perp^{\rm matéria}+P_\perp^{\rm campo}$; isotrópica (campo emaranhado),
$P_\parallel-\tfrac23\mathcal MB+(P_\parallel^{\rm campo}+2P_\perp^{\rm campo})/3$.
$\mathcal MB$ por linha vai para `<saída>_diag.txt`; é $P_\parallel$ que obedece
$dP/d\mu_n=n_B$.

Com campo constante de $10^{18}$ G (perfil `Constant`), $\mathcal MB$ chega a ~40%
da pressão da matéria em $n_B\sim0.04\,n_0$ e $P_\perp$ deixa de crescer (e pode ficar
negativa) com a densidade. Os critérios de validade da varredura usam a pressão sem o
termo de magnetização, que é monótona em $\mu$; a TOV ordena a EoS por $\epsilon$ e descarta
os trechos em que $P_\perp$ não cresce (matéria uniforme instável, região da crosta). Com os perfis dependentes de
densidade a matéria diluída vê $\sim B_{\rm surf}$ e o efeito desaparece
($|\mathcal MB/P|\lesssim2\%$ para BDD até $B_0=5\times10^{18}$ G).

### Tensões do campo e eletrodinâmica não linear

As partículas carregadas acoplam ao potencial vetor $A_\mu$ (acoplamento mínimo), portanto
o espectro de Landau depende de $B=\nabla\times A$ em todos os modelos. A eletrodinâmica
não linear altera apenas a energia e as tensões do próprio campo; no mesmo $B$, a matéria
é idêntica à do caso de Maxwell (teste `nlem_changes_only_the_field_stress`). O campo
prescrito (`bg` ou o perfil) é interpretado como $B$.

Para um campo magnético estático puro com Lagrangiana $L(B)$, a densidade de
energia é $\epsilon_B=-L$ e $H=d\epsilon_B/dB$. O tensor de tensões
$\sigma_{ij}=H_iB_j-\delta_{ij}(HB-\epsilon_B)$ dá

$$
P_\parallel=-\epsilon_B,\qquad P_\perp=HB-\epsilon_B
$$

(Soleng, PRD 52, 6178 (1995), Eq. 3, no caso logarítmico). A topologia
anisotrópica usa $P_{\rm mag}=P_\perp$; a isotrópica (campo emaranhado) usa
$(P_\parallel+2P_\perp)/3$.

| Modelo | $\epsilon_B/\epsilon_{\rm Maxwell}$ | $HB/\epsilon_{\rm Maxwell}$ | $P_\perp$ |
|---|---|---|---|
| Maxwell | 1 | 2 | $\epsilon_B$ |
| ModMax($\gamma$) | $e^{-\gamma}$ | $2e^{-\gamma}$ | $\epsilon_B$ |
| Log($\xi$) | $\ln(1+x)/x$ | $2/(1+x)$ | $\epsilon_{\rm Maxwell}\left[\frac{2}{1+x}-\frac{\ln(1+x)}{x}\right]$ |

com $\epsilon_{\rm Maxwell}=B^2/8\pi$ e $x=B^2/(2\xi^2)$. Para Maxwell e ModMax
($\epsilon_B\propto B^2$) as relações $P=\epsilon_B$ e $P=\epsilon_B/3$ continuam
valendo. Para o modelo logarítmico elas não valem: $P_\perp$ fica menor que
$\epsilon_B$ e se torna negativa para $x\gtrsim3.9$ ($B\gtrsim2.8\,\xi$).

O modelo Log é a Lagrangiana de Gaete e Helayël-Neto (EPJC 74, 3182 (2014)),
$L=-\beta^2\ln\left(1-\mathcal F/\beta^2-\mathcal G^2/2\beta^4\right)$ com
$\mathcal F=(E^2-B^2)/2$ e $\mathcal G=\mathbf E\cdot\mathbf B$, e $\xi\equiv\beta$ em
Gauss. Com $E=0$, $\mathcal G=0$ e ela coincide com a forma logarítmica de Soleng
(1995), de onde vêm as tensões acima; o termo $\mathcal G^2$ só importaria com campo
elétrico ou na propagação de fótons (birrefringência), não na estrutura estelar.

## Propriedades estelares: maré, momento de inércia e redshift

Junto com a TOV ($P$, $m$, $m_B$) o integrador resolve, com os mesmos passos (o
controle de erro usa só $P$, $m$, $m_B$; a curva M-R não muda), em unidades
geometrizadas:

**Maré** (Hinderer, ApJ 677, 1216 (2008); Postnikov, Prakash & Lattimer, PRD 82,
024016 (2010)): $r\,y'=-y^2-yF-r^2Q$, $y(0)=2$, com

$$
F=\frac{1-4\pi r^2(\epsilon-P)}{1-2m/r},\quad
Q=\frac{4\pi\left[5\epsilon+9P+(\epsilon+P)\,d\epsilon/dP\right]}{1-2m/r}
-\frac{6}{r^2(1-2m/r)}-4\left[\frac{m+4\pi r^3P}{r^2(1-2m/r)}\right]^2 .
$$

Na superfície, $y_R\to y_R-3\epsilon_s/\bar\epsilon$ ($\bar\epsilon=3M/4\pi R^3$) pela
descontinuidade de densidade (Damour & Nagar 2009). $k_2$ segue da fórmula fechada em
$C=M/R$ e $y_R$ (limite newtoniano $(2-y)/2(3+y)$ para $C<5\times10^{-3}$) e
$\Lambda=\tfrac23k_2C^{-5}$.

**Momento de inércia** (Hartle, ApJ 150, 1005 (1967)):
$\frac{1}{r^4}(r^4j\bar\omega')'+\frac{4j'}{r}\bar\omega=0$, $j=e^{-\nu/2}\sqrt{1-2m/r}$,
$\bar\omega(0)=1$. Fora da estrela $\bar\omega=\Omega-2J/r^3$, logo $J=R^4\bar\omega'(R)/6$,
$\Omega=\bar\omega(R)+2J/R^3$ e $I=J/\Omega$; $\bar I=I/M^3$.

**Redshift** de superfície: $z=(1-2C)^{-1/2}-1$.

Validação (`tests/stellar_properties.rs`): polítropo $n=1$ com $C\sim10^{-4}$,
$k_2=(15-\pi^2)/2\pi^2$ e $I=\tfrac23(1-6/\pi^2)MR^2$; densidade uniforme, $k_2=3/4$ e
$I=\tfrac25MR^2$; relação universal I-Love de Yagi & Yunes, Science 341, 365 (2013),
dentro de 1.5% para GM1, GM3 e FSU2 entre $1\,M_\odot$ e $M_{\max}$ (medido: <0.8%).
Com crosta, GM1 dá $\Lambda_{1.4}\approx850$.

Com `with_eos_output("x.dat")` o solver grava também `x_stars.txt` (EoS do núcleo +
crosta BPS; colunas M, R, $M_B$, $P_c$, $C$, $z$, $k_2$, $\Lambda$, $I$ [$10^{45}$ g cm²],
$\bar I$) e `x_diag.txt` (diagnósticos por linha da EoS). As colunas M-R anexadas a
`x.dat` continuam sem crosta, como antes.

## Propriedades de saturação e diagnósticos da EoS

`core::nuclear::saturation_properties` resolve matéria simétrica ($\mu_e=0$) com o
mesmo motor e localiza $P(\mu)=0$ no ramo denso: $E/A=\mu_n-\bar m_N$,
$K=9n_0/(dn/d\mu)$ (pois $P=n^2\,d(E/A)/dn$),
$J=k_F^2/6E_F^*+C_{\rho,\rm ef}^2n/8$ com $C_{\rho,\rm ef}^2=C_\rho^2/(1+2\Lambda_vC_\rho^2v_\omega^2)$
e $L=3n_0\,dJ/dn$. Reproduz Chen & Piekarewicz (2014) para FSU2 (K = 237.5, J = 37.56,
L = 112.6 MeV; artigo: 238.0, 37.62, 112.8) e Glendenning & Moszkowski (1991) para
GM1/GM3.

`io_utils::derived_diagnostics` (gravado em `x_diag.txt`) dá por linha $c_s^2=dP/d\epsilon$,
$\Gamma=(\epsilon+P)/P\,c_s^2$, frações $Y_p$, $Y_e$, $Y_\mu$, $Y_{\rm hyp}$ e o critério de URCA
direto nucleônico $k_{Fn}\le k_{Fp}+k_{F\ell}$ (Lattimer et al., PRL 66, 2701 (1991); momentos
de Fermi isotrópicos). O comando `nsrs report properties` resume saturação, estrelas, limiar de URCA e
início dos hyperons para GM1, GM3 e FSU2.

## Validação contra a literatura (Nível 2)

`tests/literature.rs` e `core::nuclear` reproduzem, com as parametrizações do NSRS:

| Modelo | Grandeza | NSRS | Referência |
|---|---|---|---|
| GM1 | $K$, $J$, $L$ (MeV) | 299.7, 32.48, 93.9 | 300.50, 32.52, 94.04 (Nam, Lim & Holt, arXiv:2510.15356, Tab. III) |
| GM3 | $K$, $J$, $L$ (MeV) | 239.8, 32.47, 89.6 | 240.04, 32.51, 89.75 (idem) |
| FSU2 | $n_0$, $E/A$, $M^*/M$, $K$, $J$, $L$ | 0.1503, −16.26, 0.593, 237.5, 37.56, 112.6 | 0.1505, −16.28, 0.593, 238.0, 37.62, 112.8 (Chen & Piekarewicz 2014) |
| GM1 | $M_{\max}$ só núcleons | 2.359 $M_\odot$ | 2.363 (Nam, Lim & Holt) |
| GM3 | $M_{\max}$ só núcleons | 2.015 $M_\odot$ | 2.018 (Nam, Lim & Holt) |
| FSU2 | $M_{\max}$ só núcleons | 2.071 $M_\odot$ | 2.07 ± 0.02 (Chen & Piekarewicz) |

Estrelas só com núcleons usam `with_hyperons(false)`. Para EoS rígidas a malha em
$\mu_n$ deve ir além do padrão (1.8 $M_N$): o GM1 só com núcleons termina em $4.9\,n_0$ e
a massa máxima sai truncada (2.342 em vez de 2.359). Os testes e o comando
`nsrs report properties` usam $\mu_n\le3M_N$ e verificam que o máximo não está no fim da sequência.

Diferença conhecida: para FSU2, Chen & Piekarewicz obtêm $R_{1.4}=14.42\pm0.26$ km com
uma interpolação politrópica entre a crosta externa BPS e o núcleo; o NSRS junta a
tabela BPS diretamente à EoS uniforme e obtém 13.95 km. Uma crosta interna unificada
é um passo pendente.

## Confronto com observações (Nível 3)

Vínculos em `input/observations/constraints.csv`, cada um com referência, DOI e arXiv
conferidos na fonte: massas de PSR J0348+0432 (Antoniadis et al. 2013) e PSR J0740+6620
(Fonseca et al. 2021); pontos M-R de NICER para J0030+0451 (Riley et al. 2019; Miller et
al. 2019) e J0740+6620 (Riley et al. 2021; Miller et al. 2021); $R_{1.4}=12.45\pm0.65$ km
(Miller et al. 2021); $\Lambda_{1.4}=190^{+390}_{-120}$ a 90% (LVC, PRL 121, 161101 (2018));
$n_0$, $E/A$, $K$ (Margueron, Hoffmann & Casali, PRC 97, 025805 (2018)); $J$, $L$ (Oertel et
al., RMP 89, 015007 (2017)).

`core::observations` converte cada vínculo numa distância $d$ em desvios-padrão (barras
assimétricas; intervalos de 90% divididos por 1.645): massa máxima,
$d=\max(0,(M_{\rm obs}-M_{\max})/\sigma_-)$; ponto M-R, menor distância normalizada ao ramo
estável; $R_{1.4}$ e $\Lambda_{1.4}$ interpolados na curva; propriedades nucleares, desvio
simples. $d\le1$ compatível, $1<d\le2$ tensão, $d>2$ excluído. O comando `nsrs report observations`
avalia GM1, GM3 e FSU2 (B = 0, crosta BPS, com e sem hyperons) e grava
`results/observations_report.csv`.

Resultado (excluídos, $d>2$):

| Modelo | com hyperons | só núcleons |
|---|---|---|
| GM1 | $\Lambda_{1.4}$ (890), $K$ (300) | $\Lambda_{1.4}$, $K$ |
| GM3 | massas de J0348 e J0740 ($M_{\max}=1.70$), pontos NICER de J0740 | nenhum (6 em tensão) |
| FSU2 | massas ($M_{\max}=1.60$), NICER J0740, $\Lambda_{1.4}$ (738) | $R_{1.4}$ (13.95), $\Lambda_{1.4}$ (866) |

Com hyperons (acoplamentos universais $x_\sigma=0.7$, $x_\omega=x_\rho=0.783$) GM3 e FSU2 não
sustentam $2\,M_\odot$ (problema dos hyperons). GM1 e FSU2 são rígidos demais para o
GW170817. Os raios dependem da junção crosta-núcleo (ver a nota sobre FSU2 acima).

## Setor escuro fermiônico

O setor escuro é opcional em `HadronsMatter` (builders `with_y_chi`, `with_m_chi`,
`with_m_x`, `with_g_d`, `with_epsilon`; `DarkPhotonsMatter` é um apelido) e acrescenta
um férmion de Dirac eletricamente neutro $\chi$
e um fóton escuro físico massivo $X_\mu$. Após diagonalizar a mistura cinética,
a convenção usada é

$$
\mathcal L_{\rm int}=\frac{g_Dj_D^\mu+\epsilon J_{\rm EM}^\mu}
{\sqrt{1-\epsilon^2}}X_\mu.
$$

Assim, o deslocamento de uma partícula visível de carga $q_i$ é
$\Delta_{X,i}=\epsilon e q_iX_0/\sqrt{1-\epsilon^2}$. O código usa
$E_{F,b}^*=\mu_b-V_b-\Delta_{X,b}$ e
$E_{F,\ell}=\mu_e-\Delta_{X,\ell}$; para elétrons e múons isto resulta em
$E_{F,\ell}=\mu_e+\epsilon eX_0/\sqrt{1-\epsilon^2}$.

A abundância escura é prescrita por ponto como $n_\chi=Y_\chi n_B$. O solver
determina simultaneamente $(\mu_e,S,W,R,X_0)$ e inclui a equação completa

$$
m_X^2X_0=\frac{g_Dn_\chi+\epsilon e n_Q}{\sqrt{1-\epsilon^2}},
$$

sem substituir $n_Q=0$ antes da solução. O gás escuro não recebe Landau, AMM
ou acoplamentos mesônicos. Seu potencial químico é
$\mu_\chi=\sqrt{k_{F\chi}^2+m_\chi^2}+g_DX_0/\sqrt{1-\epsilon^2}$, com
$k_{F\chi}=(3\pi^2n_\chi)^{1/3}$. Energia e pressão totais incluem uma única
contribuição vetorial $m_X^2X_0^2/2$.

Para $E_{F\chi}=\sqrt{k_{F\chi}^2+m_\chi^2}$, as contribuições cinéticas são

$$
\epsilon_\chi^{\rm kin}=\frac{k_FE_F(2k_F^2+m_\chi^2)
-m_\chi^4\ln[(k_F+E_F)/m_\chi]}{8\pi^2},
$$

$$
P_\chi^{\rm kin}=\frac{k_FE_F(2k_F^2-3m_\chi^2)
+3m_\chi^4\ln[(k_F+E_F)/m_\chi]}{24\pi^2},
$$

e satisfazem $\epsilon_\chi^{\rm kin}+P_\chi^{\rm kin}=n_\chi E_{F\chi}$.
A pressão total é calculada uma única vez pela relação de Gibbs

$$
P=\sum_b\mu_bn_b+\mu_e(n_e+n_\mu)+\mu_\chi n_\chi-\epsilon,
$$

que, sob neutralidade e a equação de Proca, recupera
$P_{\rm dark}=P_\chi^{\rm kin}+m_X^2X_0^2/2$ sem dupla contagem.

### Normalização dos termos vetoriais não lineares (FSU2)

Internamente, `gv` e `gr` representam $C_i=g_iM_N/m_i$ e os campos são
escalados como $v_\omega=g_\omega\omega_0/M_N$ e $v_\rho=g_\rho b_0/M_N$.
Para a Lagrangiana de Chen & Piekarewicz,
$\mathcal L\supset\frac{\zeta}{24}(g_\omega^2\omega_\mu\omega^\mu)^2
+\Lambda_v(g_\rho^2\mathbf b_\mu\cdot\mathbf b^\mu)(g_\omega^2\omega_\nu\omega^\nu)$,
com `rxi` $=\zeta/6$ e `lambda_v` $=\Lambda_v$, as equações de campo são

$$
v_\omega=C_\omega^2\left(\sum_b x_{\omega b}n_b-\texttt{rxi}\,v_\omega^3
-2\Lambda_v v_\omega v_\rho^2\right),\qquad
v_\rho=C_\rho^2\left(\sum_b x_{\rho b}I_{3b}n_b-2\Lambda_v v_\rho v_\omega^2\right),
$$

com o mesmo fator $C^2$ que multiplica `rb` e `rc` na equação do $\sigma$.
Como $\epsilon=\sum_b\epsilon_b^{\rm kin}+g_\omega\omega_0n_B+g_\rho b_0n_3-\mathcal L_{\rm mésons}$,
os termos vetoriais entram na densidade de energia como

$$
\epsilon\supset\frac{v_\omega^2}{2C_\omega^2}+\frac{v_\rho^2}{2C_\rho^2}
+\frac{\zeta}{8}v_\omega^4+3\Lambda_v v_\omega^2v_\rho^2 .
$$

Com essa forma, FSU2 reproduz em matéria simétrica $n_0\simeq0.1505$ fm$^{-3}$,
$E/A\simeq-16.28$ MeV e $M^*/M\simeq0.593$, e a pressão de Gibbs coincide com
$n^2\,\partial(\epsilon/n)/\partial n$ (testes em `physics.rs`). GM1 e GM3 não
são afetados porque `rxi` $=$ `lambda_v` $=0$. A matéria visível tem
uma única implementação (`physics.rs`, `particles.rs`, `eos.rs`), de modo que $Y_\chi=0$ continua
recuperando o caminho hadrônico.

### Limites desta etapa de validação

Os testes novos isolam primeiro o setor microscópico em `B=0`, isto é, sem
Landau nem AMM. Para preservar o perfil magnético legado, a exportação ainda
adiciona o piso macroscópico `bsurf = 10^11 T` mesmo quando `bg = 0`; portanto o
limite de vácuo automatizado demonstra
$n_\chi,X_0,\epsilon_{\rm dark},P_{\rm dark}\to0$, mas a energia e a pressão
totais exportadas conservam o pequeno piso eletromagnético legado.

Os comandos de campanha `nsrs dark scan` e `nsrs dark single` continuam
configurados com $B=10^{17}$ G. Essas campanhas magnetizadas são uma etapa
posterior à validação limpa do setor escuro e ainda exigem testes dedicados da
interação com Landau, AMM, perfil de campo e anisotropia antes de produção em
larga escala. O scan padrão usa somente GM1/GM3 e registra explicitamente
$m_\chi$ e $B$ no `summary.csv`.

---

## Outros Modelos

- Campo magnético dependente da densidade
- Topologias magnéticas (isotrópica/anisotrópica)
- Eletromagnetismo não linear (NLEM: Maxwell, ModMax, Log)
- Sistema TOV para sequências massa-raio

Para a documentação técnica operacional, consulte [DOCUMENTATION.md](DOCUMENTATION.md).
