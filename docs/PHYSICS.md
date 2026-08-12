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

## Setor escuro fermiônico

`DarkPhotonsMatter` acrescenta um férmion de Dirac eletricamente neutro $\chi$
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

### Nota de normalização RMF legada

Internamente, `gv` e `gr` representam $C_i=g_iM_N/m_i$, enquanto os termos
não lineares de `equation_omega()` e `equation_rho()` mantêm a normalização
histórica do projeto. Para FSU2 (`rxi` ou `lambda_v` não nulos), essas equações
não são a derivada exata do funcional atualmente usado em `compute()`; na forma
escalada faltariam fatores $C_\omega^2$ ou $C_\rho^2$ nos termos não lineares.
GM1 e GM3 não são afetados porque ambos os coeficientes são zero. A presente
alteração preserva deliberadamente essa normalização para que $Y_\chi=0$
recupere o caminho hadrônico existente; uma correção deve modificar e validar
em conjunto `physics.rs`, `eos.rs`, `darkphotons.rs` e as parametrizações.

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
