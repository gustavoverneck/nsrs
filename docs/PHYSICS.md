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
| `Dexheimer2017` | $B(\mu_B)$ local | mesmo $B(\mu_B)$ | Dexheimer et al., PLB 773, 487 (2017) |

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
$\mathcal MB$ por linha vai para `<saída>_diag.dat`; é $P_\parallel$ que obedece
$dP/d\mu_n=n_B$.

Com campo constante de $10^{18}$ G (perfil `Constant`), $\mathcal MB$ chega a ~40%
da pressão da matéria em $n_B\sim0.04\,n_0$ e $P_\perp$ deixa de crescer com a
densidade; a EoS é encerrada com `NonMonotonic`. Com os perfis dependentes de
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

Os executáveis de campanha `darkphotons` e `single_darkphotons` continuam
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
