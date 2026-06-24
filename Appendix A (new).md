# Appendix A: Preliminaries

This appendix records the circuit model, Sergeev's decompositions, and computation-tree language used throughout the paper.

## Appendix A.1: Depth-2 Linear Circuits

In this paper, we define a "depth-2 linear circuit" as an ordered pair $C=(A,B)$ of matrices $A\in\operatorname{Mat}_{n\times r}$ and $B\in\operatorname{Mat}_{r\times m}$. We say that $C$ computes the $n\times m$ matrix $M=AB$. This definition corresponds to depth-2 linear circuits in the following sense:

- Consider a circuit with three layers: an input layer, a middle layer, and an output layer, with $n,r,m$ gates, respectively.
- Each middle gate computes a fixed linear combination of the input gates.
- Each output gate computes a fixed linear combination of the middle gates.
- The matrices $A$ and $B$ represent the weights from the input gates to the middle gates and from the middle gates to the output gates, respectively.
- Because the input-output behavior of the circuit is exactly the linear transformation represented by the matrix $M=AB$, saying that $C$ computes $M$ is justified.

In this convention, inputs are row vectors: an input $x$ is transformed as $x\mapsto xA\mapsto xAB$.

For a depth-2 linear circuit $C=(A,B)$, its size is defined as
$$
\mathrm{S}(C):=\operatorname{nnz}(A)+\operatorname{nnz}(B).
$$
Its maximum input-output degree is defined as
$$
\mathrm{D}(C):=
\max\left\{
\max_x\operatorname{nnz}(A_{x,*}),
\max_y\operatorname{nnz}(B_{*,y})
\right\}.
$$
These are the conventions used in [JS13] and [AL25, Section 3.2].

Let $\mathcal R$ be a commutative semiring containing $0,1$, and let $\Gamma\subseteq\mathcal R$ be a multiplicatively closed subset containing $0,1$. Let $\textsf{Circ}_{\mathcal R,\Gamma}^{(2)}$ be the set of all depth-$2$ linear circuits whose weights belong to $\Gamma$.

For a matrix $M$ whose entries belong to $\mathcal R$, define
$$
\mathrm{S}_{\mathcal R,\Gamma}(M):=\min_{\substack{C=(A,B)\in\textsf{Circ}_{\mathcal R,\Gamma}^{(2)}\\ AB=M}}\mathrm{S}(C),
\qquad
\mathrm{D}_{\mathcal R,\Gamma}(M):=\min_{\substack{C=(A,B)\in\textsf{Circ}_{\mathcal R,\Gamma}^{(2)}\\ AB=M}}\mathrm{D}(C),
$$
with the convention that the minimum is $+\infty$ if no such circuit exists. Also define
$$
\sigma_{\mathcal R,\Gamma}(M):=\lim_{n\to\infty}\mathrm S_{\mathcal R,\Gamma}(M^{\otimes n})^{1/n},
\qquad
\delta_{\mathcal R,\Gamma}(M):=\lim_{n\to\infty}\mathrm D_{\mathcal R,\Gamma}(M^{\otimes n})^{1/n},
$$
whenever these limits exist.

By [AL25, Theorem 3.1], if $\mathbb{F}$ and $\mathbb{K}$ are fields and $\mathbb{K}$ is a field extension of $\mathbb{F}$, then
$$
\sigma_{\mathbb{F},\mathbb{F}}(M)=\sigma_{\mathbb{K},\mathbb{K}}(M),
\qquad
\delta_{\mathbb{F},\mathbb{F}}(M)=\delta_{\mathbb{K},\mathbb{K}}(M).
$$

In this paper, we use the restrictive model $\mathcal R=\mathbb Z$ and $\Gamma=\{0,\pm1\}$. Since our constructions use only coefficients in $\{0,\pm1\}$, the corresponding upper bounds on
$$
\sigma_{\mathbb F,\mathbb F}(R_1),\qquad
\delta_{\mathbb F,\mathbb F}(R_1),\qquad
\sigma_{\mathbb F,\mathbb F}(P_1^{(c)}),\qquad
\delta_{\mathbb F,\mathbb F}(P_1^{(c)})
$$
also hold, with the same numerical constants, over any field $\mathbb F$, because all integers, in particular $0$, $1$, and $-1$, can be identified in $\mathbb F$, although it is possible that $1_{\mathbb F}=(-1)_{\mathbb F}$ in fields of characteristic $2$.

We define an equivalence relation $\cong$ corresponding to relabeling the middle gates. Formally, two circuits $(A_0,B_0)$ and $(A_1,B_1)$ with $r$ middle gates are equivalent, written
$$
(A_0,B_0)\cong(A_1,B_1),
$$
if and only if there exists an $r\times r$ permutation matrix $\Pi$ such that
$$
A_1=A_0\Pi,\qquad B_1=\Pi^{-1}B_0.
$$
It is straightforward to check that $\cong$ is an equivalence relation. Let $\widetilde{\textsf{Circ}}_{\mathcal R,\Gamma}^{(2)}$ be the set of equivalence classes in $\textsf{Circ}_{\mathcal R,\Gamma}^{(2)}$ under $\cong$. From now on, we work in $\widetilde{\textsf{Circ}}_{\mathcal R,\Gamma}^{(2)}$ unless otherwise stated.

We define two binary operations on depth-2 linear circuits: concatenation sum and Kronecker product.

**Concatenation sum.** If $C=(A,B)$ and $C'=(A',B')$, where $A$ and $A'$ have the same number of rows and $B$ and $B'$ have the same number of columns, then their concatenation sum $C+C'$ is defined as
$$
C+C':=\left(\begin{bmatrix}A&A'\end{bmatrix},\begin{bmatrix}B\\B'\end{bmatrix}\right).
$$
If $C$ computes $M$ and $C'$ computes $M'$, then $C+C'$ computes $M+M'$, because
$$
\begin{bmatrix}A&A'\end{bmatrix}\begin{bmatrix}B\\B'\end{bmatrix}=AB+A'B'=M+M'.
$$
**Kronecker product.** If $C=(A,B)$ and $C'=(A',B')$, then their Kronecker product $C\otimes C'$ is defined as
$$
C\otimes C':=(A\otimes A',B\otimes B').
$$
If $C$ computes $M$ and $C'$ computes $M'$, then $C\otimes C'$ computes $M\otimes M'$, because
$$
(A\otimes A')(B\otimes B')=(AB)\otimes(A'B')=M\otimes M'.
$$
This equality follows from the mixed-product property of the Kronecker product on matrices.

Modulo the middle-gate equivalence relation $\cong$, the concatenation sum is associative and commutative. The Kronecker product on circuits is associative, but not commutative. Finally, $\otimes$ distributes over $+$:
$$
C\otimes(C'+C'')\cong(C\otimes C')+(C\otimes C'').
$$
Indeed,
$$
\left(A\otimes\begin{bmatrix}A'&A''\end{bmatrix},B\otimes\begin{bmatrix}B'\\B''\end{bmatrix}\right)\cong\left(\begin{bmatrix}A\otimes A'&A\otimes A''\end{bmatrix},\begin{bmatrix}B\otimes B'\\B\otimes B''\end{bmatrix}\right).
$$
We define one unary operation on depth-2 linear circuits: transposition.

**Transposition.** If $C=(A,B)$, then its transposition $C^\textsf{T}$ is defined as
$$
C^\textsf{T}:=(B^\textsf{T},A^\textsf{T}).
$$
If $C$ computes $M$, then $C^\textsf{T}$ computes $M^\textsf{T}$, because
$$
B^\textsf{T}A^\textsf{T}=(AB)^\textsf{T}=M^\textsf{T}.
$$
Each middle gate of a depth-$2$ linear circuit contributes one rank-$1$ term. If the middle layer has $r$ gates, then
$$
AB=\sum_{\ell=0}^{r-1}A_{*,\ell}B_{\ell,*}.
$$
We call $(A_{*,\ell},B_{\ell,*})$ the $\ell$-th rank-$1$ term of the circuit. Depending on the context, we also regard it as a circuit with one middle gate, and we call $A_{*,\ell}B_{\ell,*}$ the corresponding rank-$1$ matrix. Thus the circuit $(A,B)$ and the concatenation sum of its $r$ rank-$1$ terms are two equivalent viewpoints on the same object:
$$
(A,B)=\sum_{\ell=0}^{r-1}\,(A_{*,\ell},B_{\ell,*}).
$$
We switch between these viewpoints when discussing Sergeev's constructions and computation trees.

From now on, we fix $\mathcal R,\Gamma$ and suppress them from the notation, writing $\mathrm{S}(M), \mathrm{D}(M),\sigma(M),\delta(M)$ for $\mathrm{S}_{\mathcal R,\Gamma}(M),\mathrm{D}_{\mathcal R,\Gamma}(M),\sigma_{\mathcal R,\Gamma}(M),\delta_{\mathcal R,\Gamma}(M)$, respectively. Furthermore, we discuss only cases in which $\mathcal R$ has no zero divisors and $M$ admits at least one depth-$2$ linear circuit over $\mathcal R,\Gamma$; that is, there exists $C\in\textsf{Circ}_{\mathcal R,\Gamma}^{(2)}$ such that $C$ computes $M^{\otimes1}=M$. In particular, the matrices considered in this paper, namely $R_1$ and $P_1^{(c)}$, with $\mathcal R=\mathbb Z$ and $\Gamma=\{0,\pm1\}$, indeed satisfy these assumptions.

Under these assumptions, the size and maximum input-output degree are indeed submultiplicative under Kronecker products of matrices.

**Lemma A.1.** Under the assumptions above, for every matrix $M$ and all integers $a,b\geq1$,
$$
\mathrm{S}(M^{\otimes(a+b)})
\le
\mathrm{S}(M^{\otimes a})\mathrm{S}(M^{\otimes b}),
$$
and
$$
\mathrm{D}(M^{\otimes(a+b)})
\le
\mathrm{D}(M^{\otimes a})\mathrm{D}(M^{\otimes b}).
$$

*Proof of Lemma A.1.* For size, let $(A_a,B_a)$ and $(A_b,B_b)$ be size-optimal circuits for $M^{\otimes a}$ and $M^{\otimes b}$, respectively. Their Kronecker product $(A_a\otimes A_b,B_a\otimes B_b)$ computes $M^{\otimes(a+b)}$. Therefore
$$
\begin{aligned}
\mathrm{S}(M^{\otimes(a+b)})
&\leq
\mathrm{S}((A_a\otimes A_b,B_a\otimes B_b))\\
&=
\operatorname{nnz}(A_a\otimes A_b)+\operatorname{nnz}(B_a\otimes B_b)\\
&=
\operatorname{nnz}(A_a)\operatorname{nnz}(A_b)
+
\operatorname{nnz}(B_a)\operatorname{nnz}(B_b)\\
&\le
(\operatorname{nnz}(A_a)+\operatorname{nnz}(B_a))
(\operatorname{nnz}(A_b)+\operatorname{nnz}(B_b))\\
&=
\mathrm{S}((A_a,B_a))\mathrm{S}((A_b,B_b))\\
&=
\mathrm{S}(M^{\otimes a})\mathrm{S}(M^{\otimes b}).
\end{aligned}
$$
For maximum input-output degree, let $(A_a,B_a)$ and $(A_b,B_b)$ be degree-optimal circuits for $M^{\otimes a}$ and $M^{\otimes b}$, respectively. Their Kronecker product $(A_a\otimes A_b,B_a\otimes B_b)$ computes $M^{\otimes(a+b)}$. Therefore
$$
\begin{aligned}
\mathrm{D}(M^{\otimes(a+b)})
&\leq
\mathrm{D}((A_a\otimes A_b,B_a\otimes B_b))\\
&=
\max\{
\max\{\operatorname{nnz}((A_a\otimes A_b)_{(i,i'),*})\},
\max\{\operatorname{nnz}((B_a\otimes B_b)_{*,(j,j')})\}
\}\\
&=
\max\{
\max\{\operatorname{nnz}((A_a)_{i,*})\operatorname{nnz}((A_b)_{i',*})\},
\max\{\operatorname{nnz}((B_a)_{*,j})\operatorname{nnz}((B_b)_{*,j'})\}
\}\\
&\leq
\max\{
\max\{\mathrm{D}((A_a,B_a))\mathrm{D}((A_b,B_b))\},
\max\{\mathrm{D}((A_a,B_a))\mathrm{D}((A_b,B_b))\}
\}\\
&=
\mathrm{D}((A_a,B_a))\mathrm{D}((A_b,B_b))\\
&=
\mathrm{D}(M^{\otimes a})\mathrm{D}(M^{\otimes b}).
\end{aligned}
$$
Hence the claimed inequalities hold. $\square$

If $M$ is the zero matrix, then we set $\sigma(M)=\delta(M)=0$.

Now assume that $M$ is nonzero. Then $M^{\otimes n}$ is nonzero for every $n\ge1$, so $\mathrm{S}(M^{\otimes n})$ and $\mathrm{D}(M^{\otimes n})$ are positive integers. By Lemma A.1, the sequences
$$
\log_2\mathrm{S}(M^{\otimes n})
\qquad\text{and}\qquad
\log_2\mathrm{D}(M^{\otimes n})
$$
are subadditive, since the logarithm turns the multiplicative inequalities in Lemma A.1 into additive ones. By Fekete's lemma [Fek23],
$$
\lim_{n\to\infty}\frac1n\log_2\mathrm{S}(M^{\otimes n})
=
\inf_{n\ge1}\frac1n\log_2\mathrm{S}(M^{\otimes n}),
$$
and
$$
\lim_{n\to\infty}\frac1n\log_2\mathrm{D}(M^{\otimes n})
=
\inf_{n\ge1}\frac1n\log_2\mathrm{D}(M^{\otimes n}).
$$
Equivalently, the limits
$$
\sigma(M):=\lim_{n\to\infty}\mathrm{S}(M^{\otimes n})^{1/n},
\qquad
\delta(M):=\lim_{n\to\infty}\mathrm{D}(M^{\otimes n})^{1/n}
$$
must exist. For every nonzero matrix $M$, we have $\sigma(M),\delta(M)\ge1$.

## Appendix A.2: Sergeev's Decompositions

Sergeev's constructions [Ser22; AL25, Section 5.1] are useful rank-$1$ decompositions of the disjointness matrices, and in this paper they serve as the building blocks of the computation trees defined in Section 2.1.

In this section, for convenience, we index the rows and columns of $R_k=(R_1)^{\otimes k}$ by subsets of $[k]$, so
$$
R_k(S,T)=\mathbf 1[S\cap T=\varnothing].
$$

Fix $k\ge1$ and an identifier $c\in\{0,\ldots,2^k-1\}$. Define
$$
w(c):=(w_0,\ldots,w_{k-1},\texttt{R}),\qquad
w_i:=\begin{cases}
\texttt{R},&\text{the bit of }c\text{ at weight }2^{k-1-i}\text{ is }0,\\
\texttt{C},&\text{the bit of }c\text{ at weight }2^{k-1-i}\text{ is }1.
\end{cases}
$$
For instance, if $k=7$ and $c=43=(0101011)_2$, then $w(c)=\texttt{RCRCRCCR}$.

We now define $S_k(c)$ as the depth-$2$ linear circuit obtained by applying the following algorithm to the word $w(c)$.

1. Initialize counters $a=0$ and $b=0$, and initialize an empty list of rank-$1$ terms.

2. Read the word $w(c)$ from left to right. For each symbol:

   - If the symbol is $\texttt{R}$, perform a $(\texttt{R},a,b)$ step, and then replace $a$ by $a+1$.
   - If the symbol is $\texttt{C}$, perform a $(\texttt{C},a,b)$ step, and then replace $b$ by $b+1$.

3. The algorithm uses only the following two kinds of steps.

   - A $(\texttt{R},a,b)$ step appends, for every subset $S\subseteq[k]$ satisfying $|S|=a$, the rank-$1$ term $(U,V^\textsf{T})$ to the list, where $U$ is the indicator vector of the row $S$, and $V^\textsf{T}$ is the indicator row vector of the set of columns
     $$
     \{T\subseteq[k]\setminus S:\ |T|\ge b\}.
     $$

   - A $(\texttt{C},a,b)$ step appends, for every subset $T\subseteq[k]$ satisfying $|T|=b$, the rank-$1$ term $(U,V^\textsf{T})$ to the list, where $V^\textsf{T}$ is the indicator row vector of the column $T$, and $U$ is the indicator vector of the set of rows
     $$
     \{S\subseteq[k]\setminus T:\ |S|\ge a\}.
     $$

Finally, $S_k(c)$ is defined as a depth-$2$ linear circuit as the concatenation sum of the $2^k$ rank-$1$ terms appended.

**Fact A.2.1.** Every $S_k(c)$ is a depth-$2$ linear circuit computing $R_k$, with exactly $2^k$ middle gates.

*Proof of Fact A.2.1.* Let $k_0$ and $k_1$ denote the number of non-sentinel $\texttt{R}$ and $\texttt{C}$ symbols in $w(c)$, respectively. Then $k_0+k_1=k$. Therefore the number of appended rank-$1$ terms is
$$
\sum_{i=0}^{k_0}\binom ki+\sum_{j=0}^{k_1-1}\binom kj=
\sum_{i=0}^{k_0}\binom ki+\sum_{i=k_0+1}^{k}\binom ki=2^k.
$$


Fix a pair $(S,T)$ with $S,T\subseteq[k]$, and write $s:=|S|$ and $t:=|T|$.

If $S\cap T\neq\varnothing$, then no appended rank-$1$ term contributes to the entry $(S,T)$, since in every appended term each row set is disjoint from each column set. Hence the sum is $0$ at this entry.

If $S\cap T=\varnothing$, then a rank-$1$ term containing the entry $(S,T)$ can be appended only in one of the following two types of steps:

1. A $(\texttt{R},a,b)$ step with $a=s$ and $b\leq t$.
2. A $(\texttt{C},a,b)$ step with $a\leq s$ and $b=t$.

Indeed, for a $(\texttt{R},a,b)$ step to append a rank-$1$ term containing $(S,T)$, the term must be the one indexed by row $S$, so $a=s$ and $b\le t$; for a $(\texttt{C},a,b)$ step, it must be the one indexed by column $T$, so $a\le s$ and $b=t$.

Consider the monotone path of counters $(a,b)$ while scanning $w(c)$: a $\texttt{R}$ step moves right and a $\texttt{C}$ step moves up. The path starts at $(0,0)$ and ends at $(k_0+1,k_1)$, where $k_0+k_1=k$. Since $S\cap T=\varnothing$, we have $s+t\le k$, so the endpoint lies outside
$$
\{0,\ldots,s\}\times\{0,\ldots,t\}.
$$
Thus the path leaves this rectangle exactly once. If it leaves by a $\texttt{R}$ step, then the step has type $(\texttt{R},s,b)$ with $b\le t$, and the term indexed by row $S$ contains $(S,T)$. If it leaves by a $\texttt{C}$ step, then the step has type $(\texttt{C},a,t)$ with $a\le s$, and the term indexed by column $T$ contains $(S,T)$. These are exactly the possible contributing steps, so exactly one appended rank-$1$ term contributes to $(S,T)$. Hence the sum is $1$ at this entry. $\square$

In fact, the transpose of a Sergeev block is the Sergeev block of the complemented identifier. Set
$$
\overline c:=2^k-1-c.
$$
Since $\overline c$ flips every bit of $c$, one can obtain $w(\overline c)$ from $w(c)$ by swapping every non-sentinel $\texttt{R}$ and $\texttt{C}$.

**Fact A.2.2.** For every $k\ge1$ and every $c\in\{0,\ldots,2^k-1\}$,
$$
S_k(c)^\textsf{T}\cong S_k(\overline c).
$$

*Proof of Fact A.2.2.* We compare the appended rank-$1$ terms. In the non-sentinel prefix, passing from $c$ to $\overline c$ exchanges $\texttt{R}$ steps and $\texttt{C}$ steps. Consequently, if a non-sentinel step of $S_k(c)$ has counters $(a,b)$, then the corresponding step of $S_k(\overline c)$ has counters $(b,a)$.

A $(\texttt{R},a,b)$ step of $S_k(c)$ appends, for each $|S|=a$, the term
$$
\left(
\mathbf 1_{{S}},
\mathbf 1_{{T\subseteq[k]\setminus S:\ |T|\ge b}}
\right).
$$
Its transpose is exactly the term appended by the corresponding $(\texttt{C},b,a)$ step of $S_k(\overline c)$. Similarly, every non-sentinel $(\texttt{C},a,b)$ step of $S_k(c)$ transposes to the corresponding $(\texttt{R},b,a)$ step of $S_k(\overline c)$.

The only remaining step is the final sentinel $\texttt{R}$. If the non-sentinel prefix of $w(c)$ contains $k_0$ symbols $\texttt{R}$ and $k_1$ symbols $\texttt{C}$, then this sentinel appends
$$
\left(
\mathbf 1_{{S}},
\mathbf 1_{{[k]\setminus S}}
\right)
\qquad (|S|=k_0).
$$
In $S_k(\overline c)$, the sentinel counters are $(k_1,k_0)$, so its sentinel terms are
$$
\left(
\mathbf 1_{{T}},
\mathbf 1_{{[k]\setminus T}}
\right)
\qquad (|T|=k_1).
$$
Taking $T=[k]\setminus S$ gives exactly the transposed sentinel terms. Hence the two concatenation sums agree after transposition, up to relabeling middle gates. $\square$

## Appendix A.3: Computation Tree Constructions

In this section, we define computation trees [AL25, Appendix B] and use them to justify the recursively defined circuit families in Section 2.1.

Fix a base matrix $M$. A computation tree for $M^{\otimes N}$ is a finite rooted tree with the following data:

1. Each internal node $x$ has attached to it a finite depth-2 linear circuit $C_x$ computing $M^{\otimes n_x}$, where $n_x\ge1$.
2. The children of $x$ are in bijection with the rank-1 terms of $C_x$.
3. If a child corresponds to the rank-1 term $r$, then the edge to that child has $r$ attached to it and has weight $n_x$.
4. Every root-to-leaf path has total edge weight $N$.

The weighted depth of a node, denoted by $\operatorname{depth}(x)$, is the total edge weight on the path from the root to that node. The remaining weighted depth at a node $x$ is
$$
N_x:=N-\operatorname{depth}(x).
$$
The leaves are therefore exactly the nodes with $N_x=0$.

If $r=(U,V^\textsf{T})$ is a rank-1 term, write $\operatorname{Mat}(r):=UV^\textsf{T}$. Define a depth-2 linear circuit $\mathcal C(x)$ for each node $x$ as follows. If $x$ is a leaf, let
$$
\mathcal C(x)=([1],[1]),
$$
the trivial circuit computing $M^{\otimes0}$. If $x$ is an internal node, let
$$
\mathcal C(x)=
\sum_{y\in\operatorname{Ch}(x)}
r_{x,y}\otimes\mathcal C(y),
$$
where $r_{x,y}$ is the rank-1 term attached to the edge from $x$ to $y$.

By the distributivity of $\otimes$ over $+$, $\mathcal C(x)$ is the concatenation sum, over all leaves below $x$, of the Kronecker products of the rank-1 terms attached to the edges along the corresponding paths.

**Lemma A.3.1.** For every node $x$ in a computation tree, the subtree circuit $\mathcal C(x)$ computes $M^{\otimes N_x}$. In particular, if $\rho$ is the root, then $\mathcal C(\rho)$ computes $M^{\otimes N}$.

*Proof of Lemma A.3.1.* We prove the claim by induction from the leaves upward. If $x$ is a leaf, then $N_x=0$, and $\mathcal C(x)=([1],[1])$ computes $M^{\otimes0}$.

Now let $x$ be an internal node and assume the claim for its children. Since the children of $x$ are indexed by the rank-1 terms of $C_x$, we may write
$$
C_x=\sum_{y\in\operatorname{Ch}(x)}r_{x,y}.
$$
Because $C_x$ computes $M^{\otimes n_x}$,
$$
\sum_{y\in\operatorname{Ch}(x)}
\operatorname{Mat}(r_{x,y})
=
M^{\otimes n_x}.
$$
For every child $y$, the edge from $x$ to $y$ has weight $n_x$, so $N_y=N_x-n_x$. By induction, $\mathcal C(y)$ computes $M^{\otimes(N_x-n_x)}$. Therefore $\mathcal C(x)$ computes
$$
\begin{aligned}
\sum_{y\in\operatorname{Ch}(x)}
\operatorname{Mat}(r_{x,y})\otimes M^{\otimes(N_x-n_x)}
&=
\left(
\sum_{y\in\operatorname{Ch}(x)}
\operatorname{Mat}(r_{x,y})
\right)
\otimes M^{\otimes(N_x-n_x)}\\
&=
M^{\otimes n_x}\otimes M^{\otimes(N_x-n_x)}
=
M^{\otimes N_x}.
\end{aligned}
$$
This proves the induction step, and the statement for the root follows by taking $x=\rho$. $\square$

Finally, we connect computation trees with the time-homogeneous finite-state construction in Section 2.1.

Fix a parameter set $\mathrm{X}=(k,H,I,Z,\alpha,\beta)$ with $I=(\mathbf c,\boldsymbol{\iota})$, a multiple $n$ of $k$, and an initial state $h\in[H]$. Define a computation tree for $R_n$ with the following data:

1. Each node has a state in $[H]$, and the root has state $h$.
2. Each node of weighted depth less than $n$ is internal. If such a node has state $u$, then the Sergeev decomposition $S_k(c_{\boldsymbol{\iota}(u)})$ is attached to it.
3. The children of a node with state $u$ are in bijection with the rank-1 terms $(U_{c_{\boldsymbol{\iota}(u)},r},V_{c_{\boldsymbol{\iota}(u)},r}^\textsf{T})$ of this Sergeev decomposition. The edge corresponding to the $r$-th term has this term attached to it and has weight $k$.
4. The child corresponding to the $r$-th term has state $\operatorname{clamp}_H(u+\Delta_{c_{\boldsymbol{\iota}(u)},r})$. Nodes of weighted depth $n$ are leaves.

Since every edge has weight $k$ and $n$ is a multiple of $k$, every root-to-leaf path has total edge weight $n$.

**Proof of Proposition 2.1.** In this computation tree, every internal node has a Sergeev decomposition $S_k(c_{\boldsymbol{\iota}(u)})$ attached to it. By Fact A.2.1, every such $S_k(c_{\boldsymbol{\iota}(u)})$ computes $R_k$. Since every edge has weight $k$ and every root-to-leaf path has total weight $n$, Lemma A.3.1 proves that the root circuit computes $R_n$.

Unfolding the definition of the root circuit gives exactly the recurrence from Section 2.1:
$$
C_{\mathrm{X}}(0,u)=([1],[1]),
$$
and, for every multiple $n$ of $k$ with $n\ge k$,
$$
C_{\mathrm{X}}(n,u)
=
\sum_{r=0}^{2^k-1}
(U_{c_i,r},V_{c_i,r}^\textsf{T})
\otimes
C_{\mathrm{X}}\!\left(n-k,\operatorname{clamp}_H(u+\Delta_{c_i,r})\right),
\qquad i=\boldsymbol{\iota}(u).
$$
Therefore $C_{\mathrm{X}}(n,h)$ is the root circuit, and hence computes $R_n$. $\square$
