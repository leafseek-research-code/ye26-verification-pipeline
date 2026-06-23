# Section 1: Introduction

The Kronecker product of two matrices $A\in\operatorname{Mat}_{n_0\times m_0}$ and $B\in\operatorname{Mat}_{n_1\times m_1}$, denoted by $A\otimes B$, is defined as
$$
(A\otimes B)[(i_0,i_1),(j_0,j_1)]:=A[i_0,j_0]\cdot B[i_1,j_1].
$$
For $k\geq1$, the $k$-th Kronecker power of a matrix $A$ is defined as $A^{\otimes k}:=\underbrace{A\otimes A\otimes\cdots\otimes A}_{k\text{ times}}$.

We study linear circuits: circuits in which each gate computes a fixed linear combination of its input gates, so that the whole circuit computes a linear transform. In particular, we focus on the following question:

- Fix a base matrix $A\in\operatorname{Mat}_{\mathfrak{c}_0\times\mathfrak{c}_1}$. As $k\to+\infty$, how can one design depth-2 linear circuits that compute the linear transform represented by $A^{\otimes k}$, with asymptotically small *size* or *maximum input-output degree*?

In the case where $\mathfrak{c}_0=\mathfrak{c}_1$ and $A$ has no zero entries, this question was partially answered by Alman, Guan, and Padaki [AGP23] through the lens of *matrix rigidity upper bounds*. In this paper, however, we provide another viewpoint on this question.

In particular, our work improves both main results of Alman and Li [AL25] concerning the disjointness matrix $R_1:=\begin{pmatrix}1&1\\1&0\end{pmatrix}$, and can be generalized to other matrices, such as the $c$-ary partial-match incidence matrix $P_1^{(c)}:=\begin{pmatrix} I_c&\mathbf 1_c\end{pmatrix}$, which will be discussed later in the paper.

Our results on $R_1$ and $P_1^{(c)}$ also have applications in algorithm design. For example:

1. Since $R_1$ is the orthogonality kernel, we derive algorithms for counting orthogonal vectors (#OV). This connection is central in the recent literature on Orthogonal Vectors and #OV [NW21, Wil24, AL25].
2. Since $P_1^{(c)}$ is the incidence matrix of $c$-ary partial match, we derive data structures for the $c$-ary Partial Match problem, a generalization of the binary Partial Match problem in [CIP02].
3. By the reductions established in [AL25], we derive smaller low-depth linear circuits for disjointness and for a certain class of non-Kronecker-power matrices.

## Section 1.1: Problem, Results, and Related Work

### Section 1.1.1: Depth-2 Linear Circuits

The main object of study in this paper is depth-2 linear circuits. A depth-2 linear circuit is a linear circuit with $2+1=3$ layers: an *input layer*, a *middle layer*, and an *output layer*. Each middle gate computes a fixed linear combination of the input gates, and each output gate computes a fixed linear combination of the middle gates.

In a depth-2 linear circuit $C$, if the three layers have $n,r,m$ gates, respectively, then $C$ can be written as an ordered pair of matrices $C=(A,B)$, where the matrices $A\in\operatorname{Mat}_{n\times r}$ and $B\in\operatorname{Mat}_{r\times m}$ represent the weights from the input gates to the middle gates and from the middle gates to the output gates, respectively.

Interpreting the input as a row vector, an input $x$ is transformed as $x\mapsto xA\mapsto xAB$. Therefore, the input-output behavior of this circuit is exactly the linear transformation represented by the matrix $M=AB$; hence we say that $C$ *computes* $M$.

For a circuit $C=(A,B)$, its size [JS13] is defined as
$$
\mathrm{S}(A,B):=\operatorname{nnz}(A)+\operatorname{nnz}(B),
$$
and its maximum input-output degree [AL25, Section 3.2], abbreviated as “degree,” is defined as
$$
\mathrm{D}(A,B):=
\max\left\{
\max_x\operatorname{nnz}(A_{x,*}),
\max_y\operatorname{nnz}(B_{*,y})
\right\}.
$$
Here $\operatorname{nnz}(\cdot)$ denotes the number of nonzero elements in a given matrix or vector, since only edges with nonzero weight need to be present in the circuit. For a detailed discussion of depth-2 linear circuits, see Appendix A.1.

### Section 1.1.2: Asymptotic Size and Degree Constants

Let $\mathcal{R}$ be a commutative semiring containing $0,1$, and let $\Gamma\subseteq\mathcal{R}$ be a multiplicatively closed subset containing $0,1$. Let $\textsf{Circ}_{\mathcal{R},\Gamma}^{(2)}$ be the set of all depth-$2$ linear circuits whose weights belong to $\Gamma$.

For a matrix $M\in\mathcal{R}^{n\times m}$, we define its *size cost* $\mathrm{S}(M)$ and *degree cost* $\mathrm{D}(M)$ as the minimum size and the minimum degree, respectively, among all circuits $C\in\textsf{Circ}_{\mathcal{R},\Gamma}^{(2)}$ that compute $M$. Formally:
$$
\mathrm{S}_{\mathcal{R},\Gamma}(M):=\min_{\substack{C=(A,B)\in\textsf{Circ}_{\mathcal{R},\Gamma}^{(2)}\\ AB=M}}\mathrm{S}(C),
\qquad
\mathrm{D}_{\mathcal{R},\Gamma}(M):=\min_{\substack{C=(A,B)\in\textsf{Circ}_{\mathcal{R},\Gamma}^{(2)}\\ AB=M}}\mathrm{D}(C).
$$
We also define its *asymptotic size constant* and *asymptotic degree constant* [AL25, Section 3.2] as
$$
\sigma_{\mathcal{R},\Gamma}(M):=\lim_{n\to\infty}\mathrm S_{\mathcal{R},\Gamma}(M^{\otimes n})^{1/n},
\qquad
\delta_{\mathcal{R},\Gamma}(M):=\lim_{n\to\infty}\mathrm D_{\mathcal{R},\Gamma}(M^{\otimes n})^{1/n}.
$$

Throughout this paper, we deliberately use the restrictive model where $\mathcal{R}=\mathbb{Z}$ and $\Gamma=\{0,\pm1\}$. Since $\mathcal{R}$ and $\Gamma$ are fixed, we suppress them from the notation, and use the abbreviations $\mathrm{S}(M), \mathrm{D}(M),\sigma(M),\delta(M)$ for $\mathrm{S}_{\mathcal{R},\Gamma}(M),\mathrm{D}_{\mathcal{R},\Gamma}(M),\sigma_{\mathcal{R},\Gamma}(M),\delta_{\mathcal{R},\Gamma}(M)$, respectively [[footnote 1]].

In this convention, every integer matrix $M$ admits at least one depth-$2$ linear circuit, and the limits defining $\sigma(M)$ and $\delta(M)$ indeed exist because both depth-2 linear circuit costs are submultiplicative under tensor products. Therefore, if we define
$$
s_n:=\log_2\mathrm{S}(M^{\otimes n}),
\qquad
d_n:=\log_2\mathrm{D}(M^{\otimes n}),
$$
then $\{s_n\}_{n\ge1}$ and $\{d_n\}_{n\ge1}$ are both subadditive sequences, so Fekete's lemma [Fek23] applies to them [[footnote 2]]. For more details on these matters, see Appendix A.1.

[[footnote 1]]: Because we can identify the elements $0,1$ and $-1$ in any field, it follows that $\sigma_{\mathbb{F},\mathbb{F}}(M)\leq\sigma_{\mathbb{Z},\{0,\pm1\}}(M)$ and $\delta_{\mathbb{F},\mathbb{F}}(M)\leq\delta_{\mathbb{Z},\{0,\pm1\}}(M)$ for every integer matrix $M$. Thus, our upper bounds on $R_1$ and on $P_1^{(c)}$ hold verbatim, with the same numerical constants, over an arbitrary field $\mathbb{F}$.

[[footnote 2]]: Computing these limits is usually much harder than proving that they exist. This is a familiar theme for asymptotic constants defined by powers: tensor rank and the matrix-multiplication exponent led Strassen to the asymptotic spectrum [Str86, Str87, Str88], while graph powers lead to Shannon capacity [Sha56], Lovász’s theta-function upper bounds [Lov79], and Zuiddam’s asymptotic spectrum of graphs [Zui19]. In most cases these asymptotic constants are not even known to be computable, let alone to have closed forms.

### Section 1.1.3: Main Theorem & Related Work

Our main theorem is the following.

**Theorem 1.1 (Corollary 2.8, Corollary 2.9, Theorem 2.10).** For $R_1=\begin{pmatrix}1&1\\1&0\end{pmatrix}$, we have
$$
\sigma(R_1)<2^{1.244844}<2.37,\qquad\delta(R_1)<2^{0.319866}<1.25.
$$
Moreover, for every integer $c\ge2$, let
$$
P_1^{(c)}:=\begin{pmatrix}I_c&\mathbf{1}_c\end{pmatrix}=
\begin{pmatrix}
1&0&\cdots&0&1\\
0&1&\cdots&0&1\\
\vdots&\vdots&\ddots&\vdots&\vdots\\
0&0&\cdots&1&1
\end{pmatrix}
$$
be the $c$-ary partial-match incidence matrix. Then
$$
\sigma(P_1^{(c)})\le
c\cdot\psi(\log_2c)=c\left(1+\Theta\left(\frac{\log\log c}{\log c}\right)\right),
$$
where $\psi(r)$ is the unique root $x>1$ of $x=1+\frac{1}{x^r}$, and
$$
\delta(P_1^{(c)})\le
c^{\frac{c-1}{(c-1)+c\log_2 c}}=2-\Theta\left(\frac{1}{\log c}\right).
$$
In particular, for $P_1:=P_1^{(2)}$,
$$
\sigma(P_1)\le 1+\sqrt{5}\approx3.236,
\qquad
\delta(P_1)\le 2^{1/3}<1.26.
$$
There is a rich history of work studying upper bounds on $\sigma(R_1)$ and $\delta(R_1)$. The starting point is Yates’ algorithm [Yat37], which in our language gives $\sigma(R_1)\leq\sqrt{6}$ and $\delta(R_1)\leq\sqrt{2}$. Subsequently, for $\sigma(R_1)$, Jukna and Sergeev [JS13], Alman, Guan, and Padaki [AGP23], Sergeev [Ser22], and Alman and Li [AL25] gave increasingly sharp upper bounds. For $\delta(R_1)$, Williams [Wil24] gave an upper bound on the asymptotic equality rank $\theta(R_1)$ [[footnote 3]], while Alman and Li [AL25] gave a much sharper upper bound on $\delta(R_1)$.

| Upper Bound             | Source / Formula    | Reference                                |
| :---------------------- | ------------------- | ---------------------------------------- |
| $2^{1.293}$             | $\sqrt{6}$          | [Yat37]                                  |
| $2^{1.272}$             | $1+\sqrt{2}$        | [JS13]                                   |
| $2^{1.258}$             | order-$3$ Sergeev   | [AGP23]                                  |
| $2^{1.2504}$            | order-$15$ Sergeev  | [Ser22]                                  |
| $2^{1.2503}$            | order-$18$ Sergeev  | [Ser22], ignoring "imbalanced" condition |
| $2^{1.2495}$            | asymptotic-spectrum | [AL25]                                   |
| $\color{red}2^{1.2449}$ | Lyapunov-exponent   | **This paper**                           |

Table 1.1.1: Historical progression of upper bounds on $\sigma(R_1)$.

| Upper Bound             | Source / Formula                     | Reference                      |
| ----------------------- | ------------------------------------ | ------------------------------ |
| $2^{0.5}$               | $\sqrt{2}$                           | [Yat37]                        |
| $2^{0.4644}$            | $5^{1/5}$                            | [Wil24]                        |
| $2^{0.4151}$            | $\frac{4}{3}$                        | [AL25], if without hole-fixing |
| $2^{0.333\overline{3}}$ | $2^{1/3}$                            | [AL25]                         |
| $\color{red}2^{0.3199}$ | Lyapunov-exponent + degree-landscape | **This paper**                 |

Table 1.1.2: Historical progression of upper bounds on $\delta(R_1)$.

Write $R_k:=R_1^{\otimes k}$ for the $k$-th Kronecker power of $R_1$, equivalently the disjointness matrix on $k$-bit vectors.

Since the result of [AL25] is the sharpest among previous works and is the culmination of all previous techniques, we describe only how [AL25] proves their upper bounds: $\sigma(R_1)<2^{1.249424}$ and $\delta(R_1)\leq2^{1/3}$. At a high level, their approach can be divided into the following two steps, although the implementations differ for size and degree:

1. Start with some fixed depth-2 linear circuits for a constant Kronecker power, then
2. Give an approach for converting these circuits for a fixed power into a circuit family for all powers.

For their size result, for any fixed decomposition $R_d=\sum_i U_iV_i^\textsf{T}$, they measure its $\alpha$-volume $\rho(\alpha):=\sum_i \mathrm{nnz}(U_i)^\alpha\mathrm{nnz}(V_i)^{1-\alpha}$. Their duality theorem (component 2) says, in effect, that these $\alpha$-volumes are the right complete barriers: to prove $\sigma(R_1)\le b$, it is enough to have a collection of decompositions such that for every $\alpha\in[0,1]$, at least one decomposition has per-coordinate $\alpha$-volume at most $b$. They proved this in two ways: one nonconstructively using Strassen's asymptotic spectrum, and one constructively by building a computation tree using a bottom-up “rebalancing” dynamic program, which recovers a circuit family whose size matches the above $\alpha$-volume lower bound exactly. With this clean objective in hand, they build a better partition of $R_{18}$ (component 1) using a generalized row/column-adding strategy together with a packing/constant-weight-code construction, yielding an exponent below $1.25$.

However, our paper is able to go beyond this clean lower bound because their lower bound is tight only for the size obtained directly from a computation tree, while in our paper, after forming circuits defined by computation trees, we further feed these circuits into our size transfer theorem.

For their degree result, their path is to use a degree-uniformization theorem (component 2) for Kronecker powers. Given a fixed depth-2 linear circuit $C$ for $R_d$, they first construct the tensor power $C^{\otimes N/d}$, which computes $R_N$. Grouping rows and columns by Hamming weight, they show that within each class, typical degrees follow an expected log-degree density function. By Hoeffding concentration, only a negligible fraction of rows and columns exceed this typical value. After deleting these exceptional lines to obtain a low-degree circuit, a hole-fixing lemma exploits coordinate-permutation symmetry to cover the deleted entries with subexponentially many shifted copies, so the overhead contributes zero to the asymptotic degree exponent. They apply this theorem to a specific circuit for $R_4$ (component 1), which is formed by the tensor product of a small partition circuit for $R_2$ and its transpose. This circuit has density function at most $4/3$ everywhere, and since $(1/4)(4/3)=1/3$, they conclude $\delta(R_1)\le 2^{1/3}$.

In fact, in the language of our paper, what they do can be viewed as precisely equivalent to applying our degree transfer theorem to the singleton set $\{C_3\otimes C_4\}$.

[[footnote 3]]: The asymptotic equality rank $\theta(M)$ is another asymptotic constant, defined as $\lim\limits_{n\to\infty}\mathrm Q(M^{\otimes n})^{1/n}$, where $\mathrm{Q}(M)$ denotes the equality rank of $M$, namely the minimum $q$ such that there exists a decomposition $M=\sum_{t=1}^qa_iE^{(t)}$, where each $E^{(t)}_{ij}=\mathbf{1}[f_t(i)=g_t(j)]$. It is easy to prove that $\mathrm{D}(M)\leq\mathrm{Q}(M)$, and hence $\delta(M)\leq\theta(M)$; see [AL25, Propositions 3.1 and 3.2]. The quantity $\theta(R_1)$ defined in [Wil24] is used not only to obtain algorithms for #OV, but also to connect equality-rank lower bounds with exact-threshold circuit lower bounds. However, [AL25] points out that if the goal is only to obtain algorithms for #OV, then one can instead use the smaller quantity $\delta(R_1)$.

## Section 1.2: Technical Overview

**Size & Degree Transfer Theorems.** In this paper, we give a unified approach to the size and degree questions via the pair of *size and degree transfer theorems* (component 2). Given a matrix $M\in\operatorname{Mat}_{\mathfrak{c}_0\times\mathfrak{c}_1}$ with no identically zero row or column, and a collection $\mathcal{C}$ of depth-2 linear circuits $C$, each computing some positive Kronecker power of $M$, these theorems output upper bounds on $\sigma(M)$ and $\delta(M)$. Specifically, for a circuit $C=(A,B)$ computing $M^{\otimes k}$, let $L_C(x):=\operatorname{nnz}(A_{x,*})$ denote the degree of input gate $x$, and let $R_C(y):=\operatorname{nnz}(B_{*,y})$ denote the degree of output gate $y$. Define
$$
\begin{aligned}
f_C(p)
&:=\frac1k\,\mathbb E_{X\sim p^k}\bigl[\log_2 L_C(X)\bigr],&p\in\text{probability distribution on }[\mathfrak{c}_0].\\[4pt]
g_C(q)
&:=\frac1k\,\mathbb E_{Y\sim q^k}\bigl[\log_2 R_C(Y)\bigr],&q\in\text{probability distribution on }[\mathfrak{c}_1].
\end{aligned}
$$
For a pair of probability distributions $(p,q)\in\Delta_{[\mathfrak{c}_0]}\times\Delta_{[\mathfrak{c}_1]}$, we define the manifestation of $C$ at $(p,q)$ as $P_C(p,q):=(f_C(p),g_C(q))$, and view $P_C(p,q)$ as a point in the first quadrant of the Cartesian plane. Fixing a pair $(p,q)$, let $\operatorname{conv}\{P_C(p,q):C\in\mathcal C\}$ denote the convex hull of all points $P_C(p,q)$ over $C\in\mathcal C$. Then the size transfer theorem states that
$$
\sigma(M)\le2^\xi,\qquad\xi:=
\max_{(p,q)\in\Delta_{[\mathfrak{c}_0]}\times\Delta_{[\mathfrak{c}_1]}}\left\{
\min_{(a,b)\in\operatorname{conv}\{P_C(p,q):C\in\mathcal C\}}
\max\{\text{H}_{\mathfrak{c}_0}(p)+a,\text{H}_{\mathfrak{c}_1}(q)+b\}\right\},
$$
where $\text{H}_{\mathfrak{c}_0}(p)$ and $\text{H}_{\mathfrak{c}_1}(q)$ denote the entropies of $p$ and $q$ in bits, respectively. The degree transfer theorem states that
$$
\delta(M)\le2^\phi,\qquad\phi:=
\max_{(p,q)\in\Delta_{[\mathfrak{c}_0]}\times\Delta_{[\mathfrak{c}_1]}}\left\{
\min_{(a,b)\in\operatorname{conv}\{P_C(p,q):C\in\mathcal C\}}\max\{a,b\}\right\}.
$$

These two theorems are illustrated in detail in Section 2.3 (Theorems 2.6 and 2.7).

**Nine Circuits for $R_1$, $c+2$ Circuits for $P_1^{(c)}$.** All our size and degree results are obtained by the size and degree transfer theorems. Specifically, for $R_1$, we construct nine circuits (component 1): $C_0,C_1,C_2,C_3,C_4,C_{\mathrm J},C_{\mathrm S},C_{\mathrm R},C_{\mathrm{R}^\textsf{T}}$. Plugging the singleton set $\{C_{\mathrm J}\}$ into the size transfer theorem gives $\sigma(R_1)<2^{1.244844}$, and plugging the set $\mathcal C_8:=\{C_0,C_1,C_2,C_3,C_4,C_{\mathrm S},C_{\mathrm R},C_{\mathrm{R}^\textsf{T}}\}$ into the degree transfer theorem gives $\delta(R_1)<2^{0.319866}$. For $P_1^{(c)}$, we construct $c+2$ circuits $C'_{\textsf{pull}},C'_{\textsf{push}},C'_{(0)},\ldots,C'_{(c-1)}$. Plugging the set $\{C'_{\textsf{pull}},C'_{\textsf{push}}\}$ into the size transfer theorem gives $\sigma(P_1^{(c)})\le c\cdot\psi(\log_2 c)$, and plugging the set $\{C'_{\textsf{pull}},C'_{(0)},\ldots,C'_{(c-1)}\}$ into the degree transfer theorem gives $\delta(P_1^{(c)})\le c^{\frac{c-1}{(c-1)+c\log_2 c}}$. The fact that these two sets yield our claimed bounds for $P_1^{(c)}$ is proved in Theorem 2.10. By contrast, proving that the two sets above yield our claimed bounds for $R_1$ is much more difficult. Among the circuits above, $C_{\mathrm J},C_{\mathrm S},C_{\mathrm R},C_{\mathrm{R}^\textsf{T}}$ have especially complicated constructions: they compute $R_{7\times10^{11}}$, while $C_0,C_1,C_2$ compute $R_1$, $C_3,C_4$ compute $R_2$, and $C'_{\textsf{pull}},C'_{\textsf{push}},C'_{(0)},\ldots,C'_{(c-1)}$ all compute $P_1^{(c)}$. For more details, see Sections 2.1 and 2.2.

**Computation Trees and Scripted Rebalancing.** Computation trees [AL25, Appendix B] are rooted trees that define depth-2 linear circuits. The circuits $C_{\mathrm J},C_{\mathrm S},C_{\mathrm R},C_{\mathrm{R}^\textsf{T}}$ are all defined by computation trees that have $10^{11}+1$ layers. In those trees, each internal node $x$ has attached to it a depth-2 linear circuit $C(x)$ computing $R_k$, for $k=7$. For each rank-1 term of $C(x)$, there is one child of $x$, and the term is attached to the edge from $x$ to that child. The whole tree defines a circuit by the sum, over all leaves in the tree, of the Kronecker products of the rank-1 terms attached to the edges along the corresponding root-to-leaf paths.

To construct these trees, we use scripted rebalancing, inspired by the rebalancing idea in [AL25, Appendix B].  Here, each node has a state $h$, and its attached circuit depends only on this state. A child's state is obtained by adding the state offset of the rank-1 term on its incoming edge to its parent's state, then clamping the result by $\max\{0,\min\{H-1,\cdot\}\}$. For a rank-1 term represented by vectors $(U,V^\textsf{T})$, its state offset is the nearest-integer rounding of the base-$Z$ logarithm of the ratio between the sum of $(1-\beta,\beta)$-Bernoulli weights of the indices of nonzero entries in $V^\textsf{T}$ and the sum of $(1-\alpha,\alpha)$-Bernoulli weights of the indices of nonzero entries in $U$. Letting $I$ be the prescribed map from $[H]$ to circuits computing $R_k$, $(k,H,I,Z,\alpha,\beta)$ is called a parameter set. In this paper, we handpicked four parameter sets with $k=7$, named Jessica, Sonetto, Regulus, and Regulus Transpose. The circuits $C_{\mathrm J},C_{\mathrm S},C_{\mathrm R},C_{\mathrm{R}^\textsf{T}}$ are obtained by applying scripted rebalancing to these four parameter sets, respectively, starting from root state $0$ and growing each tree to $10^{11}+1$ layers. Each resulting circuit computes $R_{7\times10^{11}}$. For more details, see Sections 2.1 and 2.2.

**Reducing Mean Log-Degree to Lyapunov Exponent.** It remains to prove that plugging $\{C_\mathrm{J}\}$ and $\mathcal{C}_8=\{\cdots\}$ into the respective transfer theorems indeed gives the claimed upper bounds on $\sigma(R_1)$ and on $\delta(R_1)$.

Fix $p\in[0,1]$. To compute $f_{C_{\mathrm J}}(\text{Ber}(p))$, recall what this quantity means: choose a left-gate index of $C_{\mathrm J}$ whose $7\times10^{11}$ bits are i.i.d. $\text{Ber}(p)$, take the base-2 logarithm of its degree, and then take the expectation. Recall that $C_{\mathrm J}$ is defined by a computation tree with $10^{11}+1$ layers, each step consuming a block of $7$ bits. We reveal the random left-gate index one block at a time. After revealing the first several blocks, we count the nodes in the corresponding layer that are compatible with this revealed prefix: namely, those nodes for which the product of the rank-1 terms along the path from the root to the node has a left vector that is nonzero at the revealed prefix index. Since the subtree below any node is determined solely by the node's state, we only need to keep a vector indexed by $[H]$, recording, for each state, the number of compatible nodes in the current layer with that state.

The layer-to-layer update is therefore linear: the next vector is the current vector multiplied on the right by a matrix. Moreover, this matrix depends only on how many of the newly revealed $7$ bits equal $1$. Under the i.i.d. $\text{Ber}(p)$ input distribution, this number has distribution $\text{B}(7,p)$. Hence $f_{C_{\mathrm J}}(\text{Ber}(p))$ is obtained as follows: start from $[1,0,\dots,0]$, at each of the $10^{11}$ steps independently choose a matrix from $A_0,\dots,A_7$ according to $\text{B}(7,p)$, multiply the current vector on the right by that matrix, then take the expected base-2 logarithm of the $\ell_1$-norm of the final vector and divide by $7\times10^{11}$. The matrices $A_0,\dots,A_7$ are determined only by the parameter set and can be computed from it. For more details, see Section 3.1.

**Toll-and-Potential Certificates & Machine Verification.** Consider a collection of matrices $A_0,\dots,A_7$ and a distribution $\text{B}(7,p)$ on them. For the corresponding random matrix product, choose a large number of steps, randomly select one matrix at each step to multiply the current vector, and take the logarithm of the resulting vector's $\ell_1$-norm divided by the number of steps. As the number of steps tends to infinity, this quantity converges almost surely to a fixed value: the *quenched Lyapunov exponent*. This concept is studied in random dynamical systems, especially in multiplicative ergodic theory for random matrix products.

Because $10^{11}$ is so large, the numerical value of $f_{C_{\mathrm{J}}}(\mathrm{Ber}(p))$ is extremely close to this exponent divided by $7$. In this paper, we find that any pair of arrays $(\mathfrak{a},\mathfrak{b})\in\mathbb{R}^{(k+1)\times(k+1)}\times\mathbb{R}^{(k+1)\times H}$ satisfying certain conditions determined by $A_0,\dots,A_7$, called a *valid certificate*, gives a polynomial $\widetilde\Phi_\mathfrak{a}(p)$ that upper-bounds the Lyapunov exponent for every distribution $\text{B}(7,p)$. This certificate theorem is derived from the Lagrangian dual of Gharavi and Anantharam's [GA05] convex-optimization-based upper bound.

Thus, for each circuit family $\mathrm{X}$, we choose a lattice on $[0,1]$, obtain optimal certificates at the lattice points, and construct a global upper bound $\ell_{\mathrm X}:[0,1]\to\mathbb R_{\ge0}$. For each $\text{X}$, the inequality $f_{C_{\text{X}}(7m,h)}(\mathrm{Ber}(p))\le\ell_{\text{X}}(p)$ holds for all $m\geq10^{11}$, all $h\in[H]$, and all $p\in[0,1]$; moreover, $\ell_\text{X}(p)$ is $13$-Lipschitz in $p$. These two facts allow us to write two computer programs whose accepting is a sufficient, but not necessary, condition for the desired numerical bounds: $\xi<1.244844$ after plugging $\{C_\mathrm{J}\}$ into the corresponding formula, and $\phi<0.319866$ after plugging $\mathcal{C}_8$ into the corresponding formula. The implementation is available at https://github.com/leafseek-research-code/ye26-verification-pipeline. For more details, see Sections 3.2 and 3.3.

## Section 1.3: Applications in Algorithm Design

In this section, we use our upper bounds on $\sigma(R_1),\delta(R_1),\delta(P_1^{(c)})$ to derive algorithms and data structures for a variety of problems. We assume a word-RAM model, but only for basic arithmetic and memory access in $O(1)$ time; we do not use bit-parallelism.

The first example is a deterministic offline algorithm for #OV.

**Theorem 1.2 (Deterministic #OV From $\delta(R_1)$).** For every $\varepsilon>0$, there is a deterministic offline algorithm that counts orthogonal pairs between two sets $U,V\subseteq\{0,1\}^d$, each of size at most $n$, in $O(n\cdot(\delta(R_1)+\varepsilon)^d)$ time using only $O(nd)$ space. In particular, by Theorem 1.1, #OV can be solved in $O(n\cdot2^{0.319866d})=o(n\cdot2^{d/3.126})=o(n\cdot1.25^d)$ time and $O(nd)$ space.

*Proof of Theorem 1.2.* See Appendix B.1. $\square$

This improves on previous algorithms when $1.471\log_2n\leq d\leq3.125\log_2n$.

| Dimension regime                          | Descriptive label           | Time complexity, up to $n^{o(1)}$ factors | Space complexity, up to $n^{o(1)}$ factors | Main method                                                  |
| ----------------------------------------- | --------------------------- | ----------------------------------------: | -----------------------------------------: | ------------------------------------------------------------ |
| $d=o(\log n)$                             | Trivial regime              |                                    $O(n)$ |                                     $O(n)$ | Zeta transform                                               |
| $d=c\log_2n$ with $c\leq 1.47$            | Logarithmic-dimensional I   |                                $O(n+2^d)$ |                                 $O(n+2^d)$ | Zeta transform                                               |
| $d=c\log_2n$ with $1.471\leq c\leq 3.125$ | Logarithmic-dimensional II  |                    $O(n\cdot2^{d/3.126})$ |                                     $O(n)$ | Depth-2 linear circuits (this paper)                         |
| $d=c\log_2n$ with $c\geq3.126$            | Logarithmic-dimensional III |                       $n^{2-1/O(\log c)}$ |                        $n^{2-1/O(\log c)}$ | Polynomial method + rectangle matrix multiplication [AWY14, CW20] |
| $d=\omega(\log n)\cap n^{o(1)}$           | OVC barrier regime          |                                  $O(n^2)$ |                                     $O(n)$ | Brute-force pair enumeration                                 |

Table 1.3.1: Complexity landscape for #OV, using algorithms known to the authors as of June 2026.

The reduction from #OV to $\delta(R_1)$ without the linear-space requirement is given in [AL25, Theorem 7.2]. [Wil24] gives a linear-space reduction from #OV to the asymptotic equality rank $\theta(R_1)$. Because our reduction targets degree rather than equality rank, ensuring linear space requires substantially more technical work.

As a second example, we give deterministic online data structures for the Subset Query and Partial Match primitives considered by Charikar, Indyk, and Panigrahy [CIP02], as well as for a generalized $c$-ary Partial Match primitive.

In the online **Subset Query** problem, a database stores binary strings in $\{0,1\}^d$ and supports a stream of insertion and query operations. An insertion adds a string to the database. A query gives a pattern in $\{0,*\}^d$ and requires the database to return the number of currently stored strings matching the pattern. Here, a query symbol $0$ matches only a stored $0$, whereas $*$ matches either $0$ or $1$.

In the online **$c$-ary Partial Match** problem, a database stores strings in $\{0,1,\ldots,c-1\}^d$ and supports a stream of insertion and query operations. The setup is the same as for Subset Query, except that a query gives a pattern in $\{0,1,\ldots,c-1,*\}^d$. The binary case $c=2$ is the Partial Match problem defined in [CIP02].

Indeed, if we label the $c$ rows of the $c\times(c+1)$ matrix $P_1^{(c)}=\begin{pmatrix} I_c&\mathbf 1_c\end{pmatrix}$ by the database symbols $0,1,\ldots,c-1$ and its $c+1$ columns by the query symbols $0,1,\ldots,c-1,*$, then $P_1^{(c)}$ is precisely the incidence matrix for the problem.

**Theorem 1.3 (Deterministic Subset Query and Partial Match From $\delta(R_1)$ and $\delta(P_1^{(c)})$).** For every $\varepsilon>0$:

1. Subset Query in $d$ dimensions admits an online deterministic data structure in which each insertion or query takes time
   $$
   O((\delta(R_1)+\varepsilon)^d);
   $$

2. For every fixed integer $c\ge2$, $c$-ary Partial Match in $d$ dimensions admits a deterministic data structure in which each insertion or query takes time
   $$
   O((\delta(P_1^{(c)})+\varepsilon)^d)
   \le
   O\!\left(\left(c^{\frac{c-1}{(c-1)+c\log_2 c}}+\varepsilon\right)^d\right).
   $$
   In particular, by Theorem 1.1, these upper bounds are $o(1.25^d)$ for Subset Query, $o(1.26^d)$ for Partial Match, $o(1.39^d)$ for $3$-ary Partial Match, and $o(1.46^d)$ for $4$-ary Partial Match.

*Proof of Theorem 1.3.* See Appendix B.1. $\square$

Finally, by the reductions established in [AL25], we derive smaller low-depth linear circuits for disjointness and for a certain class of non-Kronecker-power matrices.

For disjointness, Theorem 1.1 implies that the $N\times N$ matrix $R_n$ has a depth-2 linear circuit of size $O(N^{1.244844})$, where $N=2^n$. More generally, applying the reduction from [AL25, Corollary 5.2], for every integer $h\ge 1$, $R_n$ has a depth-$2h$ linear circuit of size $O(N^{1+0.244844/h})$.

Moreover, for every function $f:\{0,1\}^n\to\mathbb{R}$, define
$$
M_f[x,y]:=f(x_1\vee y_1,\ldots,x_n\vee y_n),
\qquad x,y\in\{0,1\}^n .
$$
Then, applying the reduction from [AL25, Corollary 5.3], for every integer $h\ge 1$, $M_f$ has a depth-$4h$ linear circuit over $\mathbb{R}$, allowing arbitrary real weights, of size $O(N^{1+0.244844/h})$ [[footnote 4]], where $N=2^n$.

[[footnote 4]]: Indeed, let $g$ be the upper Möbius transform of $f$, defined by $g(z):=\sum\limits_{w\ge z}(-1)^{\lVert w-z\rVert_1}f(w)$. Then
$$
M_f=A\Lambda_g A^{\mathsf T},
\qquad
A:=[\mathbf 1_{x\le z}]_{x,z}
=
\begin{pmatrix}
1 & 1\\
0 & 1
\end{pmatrix}^{\otimes n},
$$
where $\Lambda_g$ is the diagonal matrix with entries $g(z)$. The matrix $A$ is, up to a permutation of rows and columns, the same Kronecker-power matrix as $R_n$. Thus one may compute the two factors $A$ and $A^{\mathsf T}$ using the circuits for $R_n$ and absorb the diagonal factor $\Lambda_g$ into an adjacent linear layer. This yields the claimed depth-$4h$ real-weight linear circuit for $M_f$ with the same size bound.
