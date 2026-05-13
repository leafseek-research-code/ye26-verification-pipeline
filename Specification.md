File structure:

```
01_write_circuit_parameters.py
02_circuit_parameters/
	01_jessica.txt
	02_sonetto.txt
	03_regulus.txt
	03_regulus_transposed.txt
03_build_transition_matrices.py
04_transition_matrices/
	01_jessica.txt
	02_sonetto.txt
	03_regulus.txt
	03_regulus_transposed.txt
05_generate_lattice_points.py
06_lattice_points/
	01_jessica.txt
	02_sonetto.txt
	03_regulus.txt
	03_regulus_transposed.txt
07_compute_dual_certificates.py
08_dual_certificates/
	01_jessica.txt
	02_sonetto.txt
	03_regulus.txt
	03_regulus_transposed.txt
09_verify_size_bound.py
10_verify_degree_bound.py
11_figures/
	01_certificate_curves.png
	02_value_heatmap_full.png
	03_value_heatmap_zoom.png
	04_strategy_heatmap_full.png
	05_strategy_heatmap_zoom.png
```

# `01_write_circuit_parameters.py`

This script writes the parameters into files in the directory `02_circuit_parameters/`. Each file contains exactly seven lines, with no trailing newline after the final line. The format is:

1. $K$: the number of bits in each block.
2. $H$: the number of tilt states in the finite window.
3. $L$: the number of available Sergeev building blocks.
4. The $L$ block identifiers $c_1,...,c_L$.
5. The $L$ multiplicities $m_1,...,m_L$.
6. $Z$: the logarithmic base used in the tilt rounding rule.
7. $\alpha,\beta$: the left and right Bernoulli weights.

The file does not store the whole index map $I$ explicitly.

Instead, lines 4 and 5 are interpreted as a piecewise-constant partition of $[H]=\{0,\dots,H-1\}$.

If the cumulative sums of the multiplicities are $M_a=m_1+\cdots+m_a$, then $I(h)=c_a$ precisely when $M_{a-1}\leq h<M_a$, with $M_0=0$.

In other words, the first $m_1$ states use block $c_1$, the next $m_2$ states use block $c_2$, and so on.

(This $I$ is exactly the $I$ inside the definition of the matrices $A_i$ in chapter `03_build_transition_matrices.py`. This needs emphasis; otherwise, in chapter 03, one could forget.)

The script must write the following four files.

### `01_jessica.txt`

The file content is

```text
7
612
6
0 20 42 85 107 127
14 22 270 270 22 14
1.2599210498948732
0.5 0.5
```

### `02_sonetto.txt`

The file content is

```text
7
360
6
0 21 42 85 106 127
8 22 150 150 22 8
1.2599210498948732
0.38 0.38
```

### `03_regulus.txt`

The file content is

```text
7
360
8
0 42 43 45 77 85 86 127
8 36 20 20 160 64 44 8
1.2599210498948732
0.336 0.414
```

### `03_regulus_transposed.txt`

The file content is

```text
7
360
8
0 41 42 50 82 84 85 127
8 44 64 160 20 20 36 8
1.2599210498948732
0.414 0.336
```

# `03_build_transition_matrices.py`

This script converts each of the parameter files in `02_circuit_parameters/` to a file in `04_transition_matrices/`.

For a set of parameters $(K,H,I,Z,\alpha,\beta)$, it creates the matrices
$$
A_0,\dots,A_K\in(\mathbb Z_{\geq 0})^{H\times H}.
$$

### 1. Reconstructing the building blocks

For each block identifier $c$ appearing in the parameter file, construct the Sergeev decomposition $\mathcal{S}(c)$ of order $K$ corresponding to $c$ as follows.

First write $c$ in binary using exactly $K$ bits. Read bit $0$ as $\texttt{R}$ and bit $1$ as $\texttt{C}$. This produces a word
$$
w(c)\in\{\texttt{R},\texttt{C}\}^K.
$$
For example, for $K=7$ and $c=43$, as $43=(0101011)_2$, it maps to the string $\texttt{RCRCRCC}$.

Append one final sentinel $\texttt{R}$, obtaining a word of length $K+1$.

Now run Sergeev's row/column procedure on $w(c)$:

1. Start with integers $r=0,c=0$, and an empty list of rank-1 terms.

2. For each position $k=1,2,\dots,K+1$:

   - If the current symbol is $\texttt{R}$, then for every subset $S\subseteq[K]$ with $|S|=r$, append the rank-1 matrix
     $$
     UV^\textsf{T}
     $$
     into the list, where $U$ is the indicator vector of the single row $S$, and $V$ is the indicator vector of all columns $T\subseteq[K]\setminus S$ with $|T|\geq c$.

     After appending all such terms, increase $r$ by $1$.

   - If the current symbol is $\texttt{C}$, then for every subset $T\subseteq[K]$ with $|T|=c$, append the rank-1 matrix
     $$
     UV^\textsf{T}
     $$
     into the list, where $V$ is the indicator vector of the single column $T$, and $U$ is the indicator vector of all rows $S\subseteq[K]\setminus T$ with $|S|\geq r$.

     After appending all such terms, increase $c$ by $1$.

The resulting list is a decomposition of $\begin{psmallmatrix}1&1\\1&0\end{psmallmatrix}^{\otimes K}$ into exactly $2^K$ rank-1 $0/1$ matrices. Here rows and columns are indexed by subsets of $[K]$, or equivalently by $K$-bit strings.

(Using bitmasks is beneficial for the implementation here.)

### 2. Weighted nonzero counts and tilts

For a vector $\vec{x}$ indexed by $K$-bit strings and parameters $(\gamma_0,\gamma_1)$, define the weighted nonzero count
$$
\operatorname{wnnz}_{\gamma_0,\gamma_1}(\vec{x})
=
\sum_{i:\,\vec{x}_i\neq 0}
(\gamma_0)^{K-\|i\|_1}
(\gamma_1)^{\|i\|_1},
$$
where $\|i\|_1$ is the Hamming weight of $i$.

Since every Sergeev term is a $0/1$ indicator vector, this is just the weighted measure of its support.

For each rank-1 term $(U_j,V_j)$ in $\mathcal{S}(c)$, compute
$$
D_j
=
\operatorname{round}\!\left(
\log_Z
\frac{\operatorname{wnnz}_{1-\beta,\beta}(V_j)}
     {\operatorname{wnnz}_{1-\alpha,\alpha}(U_j)}
\right),
$$
where $\operatorname{round}$ means nearest integer.

Also, print to the console some information about the worst rounding: that is, print information on the occurred instance where the relevant value is within the smallest absolute difference from a half-integer, making the rounding suffer the worst.

### 3. Defining the matrices $A_i$

For each $i\in\{0,1,\dots,K\}$, define the canonical row
$$
0^{K-i}1^i,
$$
viewed either as a $K$-bit string or as the subset $\{0,\dots,i-1\}$. Its integer index is $2^i-1$.

(Any fixed row of weight $i$ would give the same count.)

For every old state label $h_\text{old}\in\{0,\dots,H-1\}$:

1. Determine the active block identifier $c = I(h_\text{old})$ from the piecewise-constant map encoded in the parameter file.

2. Enumerate all rank-1 terms $(U_j,V_j)$ in the Sergeev decomposition $\mathcal{S}(c)$.

3. Keep exactly those terms for which the canonical row $0^{K-i}1^i$ lies in the support of $U_j$.

4. For each such term, compute
   $$
   h_\text{new}=\operatorname{clamp}_H(h_\text{old}+D_j).
   $$

5. Increment the matrix entry $A_i[h_\text{old},h_\text{new}]$ by $1$.

Equivalently,
$$
A_i[h_\text{old},h_\text{new}]
=
\#\Bigl\{
j:\,
0^{K-i}1^i\in\operatorname{supp}(U_j),\ 
\operatorname{clamp}_H(h_\text{old}+D_j)=h_\text{new}
\Bigr\},
$$
where the terms $(U_j,V_j)$ range over the decomposition $\mathcal{S}(I(h_\text{old}))$, and the $\text{clamp}$ function is
$$
\operatorname{clamp}_H(t)=
\begin{cases}
0,& t<0,\\
t,& 0\le t\le H-1,\\
H-1,& t>H-1.
\end{cases}
$$
Every entry is therefore a nonnegative integer. There is no need to store zero entries.

### 4. Output format

For each parameter set $(K,H,I,Z,\alpha,\beta)$, write one text file in `04_transition_matrices/`.

- The first line is

  ```text
  K H (# of entries)
  ```

  in our examples, the first two integers are just `7 360` or `7 612`, while the third integer is the number of nonzero entries total across the $K+1$ matrices $A_0,\dots,A_K$.

- Then write one nonzero entry per line in the form

  ```text
  i h_old h_new value
  ```

  meaning that $A_i[h_{\text{old}},h_{\text{new}}]=\text{value}$.

Only nonzero entries should be written, and each nonzero entry should appear exactly once.

The nonzero entries should appear first in ascending order of $i$, second in ascending order of $h_\text{old}$, and third in ascending order of $h_\text{new}$.

A typical file therefore looks like

```text
7 360 6726
0 0 14 1
...
7 359 345 1
```

# `05_generate_lattice_points.py`

This script generates the lattice files in the directory `06_lattice_points/`, writing one file per parameter set.

Each file contains a strictly increasing list of lattice points $p\in(0,1)$, one value per line.

### `01_jessica.txt`

The file is

```text
0.1
... (step 0.1)
0.3
... (step 0.01)
0.34
... (step 0.001)
0.36
... (step 0.0001)
0.3708
... (step 0.00003)
0.372
... (step 0.0001)
0.38
... (step 0.001)
0.4
... (step 0.01)
0.5
... (step 0.1)
0.9
```

### `02_sonetto.txt`

The file is

```text
0.01
... (step 0.01)
0.2
... (step 0.001)
0.5
... (step 0.01)
0.99
```

### `03_regulus.txt`

The file is

```text
0.01
... (step 0.01)
0.2
... (step 0.001)
0.32
... (step 0.01)
0.99
```

### `03_regulus_transposed.txt`

The file is

```text
0.01
... (step 0.01)
0.32
... (step 0.001)
0.44
... (step 0.01)
0.99
```

# `07_compute_dual_certificates.py`

This script is run separately for each parameter set.

Its job is to read the corresponding lattice and matrix data, solve one convex optimization problem for each lattice point, and write the resulting dual certificates to the directory `08_dual_certificates/`.

Let the lattice be $p_0 < p_1 < \cdots < p_L$.

For each lattice point $p$, define
$$
P_i=\binom{K}{i}p^i(1-p)^{K-i},\qquad i=0,\dots,K.
$$
For that fixed $p$, introduce free variables
$$
\alpha_{i,j}\qquad (0\le i,j\le K)
$$
and
$$
\beta_{i,u}\qquad (0\le i\le K,\ 0\le u\le H-1).
$$
The optimization problem is:

$$
\text{minimize}\qquad
\sum_{i=0}^{K}\sum_{j=0}^{K} P_iP_j\,\alpha_{i,j}
$$
subject to
$$
\sum_{j=0}^{K}\sum_{v=0}^{H-1}
A_j[u,v]\,2^{-\alpha_{i,j}-\beta_{i,u}+\beta_{j,v}}
\le 1,
\qquad
\text{for every } i\in\{0,\dots,K\},\ u\in\{0,\dots,H-1\}.
$$

### Implementation requirements

Use **CVXPY** to build the optimization model and **MOSEK** to solve it.

When constructing the model, do **not** try to compress the problem with clever vectorized expressions or matrix tricks. Build it directly and transparently, using the explicit sums described above.

Each constraint should be modeled as a `log_sum_exp` constraint in the natural-log domain, exactly as one would normally do for base-$e$ exponentials. In particular, use the ordinary CVXPY `log_sum_exp` machinery rather than trying to encode base-$2$ exponentials directly.

On top of everything, create $K+1$ additional constraints that
$$
\beta_{i,\frac{H}{2}-1}+\beta_{i,\frac{H}{2}}=0\qquad(i=0,\dots,K).
$$
This is to remove the $K+1$ translational gauge freedoms between the $\alpha,\beta$ array values, for sake of numerical stability.

### Post-processing the solver output

MOSEK returns only an approximate solution, so the raw output must be corrected before it is written as a certificate.

After solving:

1. Divide every $\alpha_{i,j}$ and every $\beta_{i,u}$ by $\ln 2$.

   This converts the solution back to the base-$2$ language used in the mathematical statement.

2. Manually verify every constraint by evaluating
   $$
   \sum_{j=0}^{K}\sum_{v=0}^{H-1}
   A_j[u,v]\,2^{-\alpha_{i,j}-\beta_{i,u}+\beta_{j,v}}
   $$
   for each pair $(i,u)$.

3. Let $S_0$ be the largest of these evaluated sums over all constraints.

4. Add $\log_2(S_0)+10^{-8}$ to **every** $\alpha_{i,j}$ simultaneously.

   This correction must be applied regardless of whether $S_0\le 1$ or $S_0>1$.

### Output format

For each parameter set, write one file in `08_dual_certificates/`. The file contains one block for each lattice point $p$, in lattice order.

A block has the following form:

1. One line containing the value of $p$.
2. $K+1$ lines containing the matrix $\alpha$, one row per line, with $K+1$ floating-point values on each line.
3. $K+1$ lines containing the array $\beta$, one row per line, with $H$ floating-point values on each line.

For example, the structure of `02_sonetto.txt` is

```text
0.01
(K+1 lines, each with K+1 values: the α matrix)
(K+1 lines, each with H values: the β rows)
0.02
(K+1 lines, each with K+1 values: the α matrix)
(K+1 lines, each with H values: the β rows)
...
0.99
(K+1 lines, each with K+1 values: the α matrix)
(K+1 lines, each with H values: the β rows)
```

A similar layout is used for the other families.

Formatting rules:

- Do not place a trailing newline at the end of the file.
- On lines containing multiple numbers, do not place a trailing space after the final entry.
- Print every floating-point number using the same format as `%14.10lf`, i.e. width $14$ and exactly $10$ digits after the decimal point.

# `09_verify_size_bound.py`

This program verifies the size exponent bound
$$
\boldsymbol{\sigma}:=\max_{p\in[0,1]}\bigl\{\ell_\text{Jessica}(p)+\mathrm h(p)\bigr\}\le 1.244844,
$$
where
$$
\mathrm h(p):=-p\log_2 p-(1-p)\log_2(1-p).
$$

The corrected definition of the certificate envelope is
$$
\ell_X(p)=\frac1K\Bigl(\text{best certificate objective at }p-\operatorname{Ent}(P)\Bigr),
$$
where
$$
P_i=\binom{K}{i}p^i(1-p)^{K-i},\qquad i=0,\dots,K,
$$
and
$$
\operatorname{Ent}(P):=-\sum_{i=0}^{K}P_i\log_2 P_i
$$
is the entropy, in bits, of the binomial law $P=(P_0,\dots,P_K)$.

Because $\operatorname{Ent}(P)$ has sharp slopes near $0$ and $1$, the verifier does **not** use it directly in the rigorous interval estimate.

Instead it uses a continuous piecewise-linear lower envelope $\underline E$ of $\operatorname{Ent}(P)$.
Since subtracting a smaller entropy produces a larger value of $\ell_X$, this only makes the proof harder and is therefore rigorous.

## 1. Inputs

The program reads:

- `06_lattice_points/01_jessica.txt`
- `06_lattice_points/02_sonetto.txt`
- `06_lattice_points/03_regulus.txt`
- `06_lattice_points/03_regulus_transposed.txt`
- `08_dual_certificates/01_jessica.txt`
- `08_dual_certificates/02_sonetto.txt`
- `08_dual_certificates/03_regulus.txt`
- `08_dual_certificates/03_regulus_transposed.txt`

Only Jessica is needed for the numerical proof of the size bound, but all four families are needed for the plot.

## 2. Evaluating one stored certificate at an arbitrary parameter

Fix one stored certificate, meaning one matrix
$$
\alpha=(\alpha_{i,j})_{0\le i,j\le 7}.
$$
At a fresh parameter value $p\in[0,1]$, its raw objective value is
$$
\Phi_\alpha(p)
:=
\sum_{i=0}^{7}\sum_{j=0}^{7}
\binom{7}{i}\binom{7}{j}
p^{i+j}(1-p)^{14-i-j}\,\alpha_{i,j}.
$$

This is exactly the old certificate re-evaluated against the new binomial law.
No feasibility check is needed here: feasibility was already enforced when the certificate file was created, and the constraints do not depend on $p$.

For efficiency, precompute the degree-$14$ Bernstein coefficients
$$
b_t(\alpha)
:=
\frac{1}{\binom{14}{t}}
\sum_{\substack{0\le i,j\le 7\\ i+j=t}}
\binom{7}{i}\binom{7}{j}\alpha_{i,j},
\qquad 0\le t\le 14.
$$
Then
$$
\Phi_\alpha(p)=\sum_{t=0}^{14} b_t(\alpha)\binom{14}{t}p^t(1-p)^{14-t}.
$$

## 3. The entropy correction and its lower envelope

For $K=7$, define
$$
P_i(p)=\binom{7}{i}p^i(1-p)^{7-i},
\qquad
E(p):=\operatorname{Ent}(P(p))
=
-\sum_{i=0}^{7}P_i(p)\log_2 P_i(p).
$$

The function $E$ is symmetric under $p\leftrightarrow 1-p$, satisfies $E(0)=E(1)=0$, and is concave on $[0,1]$.
Define the anchor value
$$
E_{0.01}:=E(0.01),
$$
and the lower envelope
$$
\underline E(p)=
\begin{cases}
\frac{E_{0.01}}{0.01}\,p, & 0\le p\le 0.01,\\[1ex]
E(p), & 0.01\le p\le 0.99,\\[1ex]
\frac{E_{0.01}}{0.01}\,(1-p), & 0.99\le p\le 1.
\end{cases}
$$

This is a valid lower envelope:

- on $[0.01,0.99]$ it is equal to $E$ by definition;
- on $[0,0.01]$ it is the secant from $(0,0)$ to $(0.01,E(0.01))$, which lies below the concave graph of $E$;
- on $[0.99,1]$ the same follows by symmetry.

Therefore
$$
\underline E(p)\le E(p)
\qquad\text{for all }p\in[0,1].
$$

For one stored certificate $\alpha$, define its corrected evaluation curve by
$$
\widetilde\Phi_\alpha(p):=\frac{1}{7}\bigl(\Phi_\alpha(p)-\underline E(p)\bigr).
$$

## 4. Definition of the four corrected envelope curves

For any one of the four families
$$
X\in\{\text{Jessica},\text{Sonetto},\text{Regulus},\text{Regulus}^\mathsf T\},
$$
let its lattice be
$$
p_0<p_1<\cdots<p_L,
$$
and let $\alpha^{(0)},\dots,\alpha^{(L)}$ be the stored certificates in that order.

Define the corrected certificate-evaluation curves
$$
\widetilde\Phi^{(r)}_X(p):=\frac{1}{7}\bigl(\Phi^{(r)}_X(p)-\underline E(p)\bigr),
\qquad 0\le r\le L.
$$

Then define $\ell_X:[0,1]\to\mathbb R$ by using the best adjacent certificate:

- if $p\in[p_r,p_{r+1}]$, set
  $$
  \ell_X(p):=\min\bigl\{\widetilde\Phi^{(r)}_X(p),\widetilde\Phi^{(r+1)}_X(p)\bigr\};
  $$

- if $p\in[0,p_0]$, set
  $$
  \ell_X(p):=\widetilde\Phi^{(0)}_X(p);
  $$

- if $p\in[p_L,1]$, set
  $$
  \ell_X(p):=\widetilde\Phi^{(L)}_X(p).
  $$

Equivalently,
$$
\ell_X(p)=\frac{1}{7}\Bigl(\min\{\text{usable raw certificate values at }p\}-\underline E(p)\Bigr),
$$
because the same entropy correction is subtracted from every certificate.

## 5. A computable Lipschitz bound for the corrected curves

For a degree-$m$ Bernstein polynomial
$$
F(p)=\sum_{t=0}^{m} b_t\binom{m}{t}p^t(1-p)^{m-t},
$$
one has
$$
F'(p)=
m\sum_{t=0}^{m-1}(b_{t+1}-b_t)\binom{m-1}{t}p^t(1-p)^{m-1-t}.
$$
Hence
$$
|F'(p)|
\le
m\max_{0\le t\le m-1}|b_{t+1}-b_t|
\qquad\text{for all }p\in[0,1].
$$

Apply this with $m=14$. For one stored certificate $\alpha$, define
$$
L(\alpha):=
14\max_{0\le t\le 13}
\bigl|b_{t+1}(\alpha)-b_t(\alpha)\bigr|.
$$
Then $\Phi_\alpha$ is $L(\alpha)$-Lipschitz on $[0,1]$.

Now define
$$
S_E:=\frac{E(0.01)}{0.01}.
$$
Because $E$ is concave and symmetric, the lower envelope $\underline E$ is globally $S_E$-Lipschitz on $[0,1]$.

Therefore every corrected curve $\widetilde\Phi_\alpha$ is
$$
\frac{L(\alpha)+S_E}{7}
$$
-Lipschitz on $[0,1]$.

Define the family-wide constants
$$
L_\text{Jessica}:=\max_r \frac{L\bigl(\alpha^{(r)}_\text{Jessica}\bigr)+S_E}{7},
$$
and similarly
$$
L_\text{Sonetto},\quad L_\text{Regulus},\quad L_{\text{Regulus}^\mathsf T}.
$$

Because the minimum of finitely many $L$-Lipschitz functions is again $L$-Lipschitz, each envelope $\ell_X$ is $L_X$-Lipschitz.

Note: Instead of directly using those as the Lipschitz constant in the verification process (of both the size bound and the degree bound), the programs should check if all four constants are $\leq13$, then use $13$ as the Lipschitz constant in the verification process. When reporting this effective verification Lipschitz constant to the console, both programs should print it in the format `L_*=`.

## 6. Certification of the bound $\boldsymbol{\sigma}\le 1.244843$

Use the stricter threshold $1.244843$ in place of $1.244844$.

Take the Jessica lattice, adjoin the two endpoint intervals $[0,p_0]$ and $[p_L,1]$, and subdivide **every** resulting interval into $1000$ equal subintervals.
This yields a fine partition of $[0,1]$ into closed subintervals $[a,b]$.

For each such $[a,b]$, compute the rigorous upper bound
$$
U[a,b]
:=
\ell_\text{Jessica}(b)
+
L_\text{Jessica}(b-a)
+
\max_{p\in[a,b]}\mathrm h(p).
$$

The binary-entropy term is evaluated exactly:
$$
\max_{p\in[a,b]}\mathrm h(p)
=
\mathrm h\bigl(\min\{\max\{1/2,a\},\,b\}\bigr),
$$
because $\mathrm h$ is concave and maximized at $p=1/2$.

If
$$
U[a,b]\le 1.244843
$$
for every fine subinterval, then the whole proof is complete.

If any interval fails this test, the program must raise an error and print the offending interval, the computed upper bound, and the endpoint values used to obtain it.

## 7. Plot `01_certificate_curves.png`

The program also draws six curves on the common domain $p\in[0,1]$:

1. $\ell_\text{Jessica}(p)$
2. $\max\{-0.1,\ell_\text{Jessica}(p)+\mathrm h(p)-1\}$
3. $\ell_\text{Sonetto}(p)$
4. $\ell_\text{Regulus}(p)$
5. the horizontal line at $1.244843-1=0.244843$
6. $\ell_{\text{Regulus}^\mathsf T}(p)$

### Sampling rule

Evaluate all six curves at
$$
p=\frac{m}{10000},\qquad m=0,1,\dots,10000,
$$
and connect consecutive sample points by straight line segments.

### Fixed colors

Use the following fixed character-inspired colors:

- Jessica:
  - `#D7E0E1`
- Jessica excess minus one:
  - `#C4CED0`
- Target minus one:
  - `#B4ABDA`
- Sonetto:
  - `#D9782B`
- Regulus:
  - `#D94B3A`
- Regulus$^\mathsf T$, darker:
  - `#9E2F2A`

### Plot requirements

Use one informative plot with:

- x-axis label `p`
- y-axis label `value`
- legend enabled
- a light grid
- x-limits exactly `[0, 1]`

A good default is to draw all lines with line width $2$.

The horizontal reference line must be labeled by the decimal value of `VERIFICATION_THRESHOLD - 1`, so with the current target it appears in the legend as `0.244843`. The program should also print the target bound to the console in the format `Target Bound: 1.244843`.

Save the image as

```text
11_figures/01_certificate_curves.png
```

# `10_verify_degree_bound.py`

This program verifies
$$
\boldsymbol{\delta}:=\max_{p,q\in[0,1]}\delta(p,q)\le 0.3199,
$$
using a recursive rectangle certifier, and also generates the four heatmaps

1. `11_figures/02_value_heatmap_full.png`
2. `11_figures/03_value_heatmap_zoom.png`
3. `11_figures/04_strategy_heatmap_full.png`
4. `11_figures/05_strategy_heatmap_zoom.png`

## 1. The eight candidate points

For each $(p,q)\in[0,1]^2$, define eight points
$$
P_s(p,q)=(f_s(p),g_s(q)),\qquad s=0,\dots,7,
$$
by
$$
\begin{aligned}
P_0(p,q) &= (0,\;1-q),\\
P_1(p,q) &= (1-p,\;0),\\
P_2(p,q) &= (p,\;q),\\
P_3(p,q) &= (p(1-p),\;\tfrac12(1-q^2)),\\
P_4(p,q) &= (\tfrac12(1-p^2),\;q(1-q)),\\
P_5(p,q) &= (\ell_\text{Sonetto}(p),\;\ell_\text{Sonetto}(q)),\\
P_6(p,q) &= (\ell_\text{Regulus}(p),\;\ell_{\text{Regulus}^\mathsf T}(q)),\\
P_7(p,q) &= (\ell_{\text{Regulus}^\mathsf T}(p),\;\ell_\text{Regulus}(q)).
\end{aligned}
$$

The three envelope curves $\ell_\text{Sonetto},\ell_\text{Regulus},\ell_{\text{Regulus}^\mathsf T}$ are defined exactly as in `09_verify_size_bound.py`, from the corresponding certificate files.

## 2. The 20 strategies

There are two kinds of strategies.

### 2.1 The two “one” strategies

Only the two one-strategies
$$
2\qquad\text{and}\qquad 5
$$
are considered.

For each such $i\in\{2,5\}$, the one-strategy using only the point
$$
P_i(p,q)=(x_i,y_i)
$$
has value
$$
\delta_i(p,q):=\max\{x_i,y_i\}.
$$

This is the whole rule for one-strategies: no segment is involved.

### 2.2 The eighteen “two” strategies

Only the following unordered pairs are considered:
$$
\begin{aligned}
&(0,1),(0,2),(0,4),(0,5),(0,6),(0,7),\\
&(1,2),(1,3),(1,5),(1,6),(1,7),\\
&(2,5),(2,6),(2,7),\\
&(3,6),(4,7),(5,6),(5,7).
\end{aligned}
$$

For each such pair $(i,j)$, write
$$
A=P_i(p,q)=(x_\mathrm I,y_\mathrm I),
\qquad
B=P_j(p,q)=(x_\mathrm{II},y_\mathrm{II}),
$$
and define
$$
d_\mathrm I:=x_\mathrm I-y_\mathrm I,
\qquad
d_\mathrm{II}:=x_\mathrm{II}-y_\mathrm{II}.
$$

Fix
$$
\mathrm{EPS}:=10^{-11}.
$$

A two-strategy is declared **valid** only when the segment $AB$ has one numerically certified real intersection with the diagonal $y=x$, in the following deliberately unfavorable sense:

- valid if
  $$
  d_\mathrm I<-\mathrm{EPS}\ \text{ and }\ d_\mathrm{II}>\mathrm{EPS},
  $$
  or
  $$
  d_\mathrm I>\mathrm{EPS}\ \text{ and }\ d_\mathrm{II}<-\mathrm{EPS};
  $$

- invalid otherwise.

Thus endpoint contacts, nearly-diagonal segments, and all floating-point ambiguous cases are treated as **no certified real intersection**.

If the pair is valid, the intersection parameter is
$$
t_*=\frac{d_\mathrm I}{d_\mathrm I-d_\mathrm{II}},
$$
which automatically lies in $(0,1)$, and the intersection point is
$$
(z,z)=(1-t_*)A+t_*B.
$$
The two-strategy value is then
$$
\delta_{i,j}(p,q):=z.
$$

If the pair is invalid, define
$$
\delta_{i,j}(p,q):=+\infty.
$$

So a two-strategy contributes only when there is one robustly certified interior crossing of the diagonal.

## 3. The value landscape and the strategy landscape

There are exactly
$$
20
$$
candidate strategies.

Define the value landscape
$$
\delta(p,q):=
\min\Bigl(
\delta_{0,1}(p,q),\delta_{0,2}(p,q),\delta_{0,4}(p,q),\delta_{0,5}(p,q),\delta_{0,6}(p,q),\delta_{0,7}(p,q),
\delta_{1,2}(p,q),\delta_{1,3}(p,q),\delta_{1,5}(p,q),\delta_{1,6}(p,q),\delta_{1,7}(p,q),
\delta_2(p,q),\delta_{2,5}(p,q),\delta_{2,6}(p,q),\delta_{2,7}(p,q),
\delta_{3,6}(p,q),\delta_{4,7}(p,q),\delta_5(p,q),\delta_{5,6}(p,q),\delta_{5,7}(p,q)
\Bigr).
$$

Use the following fixed deterministic strategy order for computation and storage:
$$
\begin{aligned}
&2,\ 5,\ 0{+}1,\ 0{+}2,\ 0{+}4,\ 0{+}5,\ 0{+}6,\ 0{+}7,\\
&1{+}2,\ 1{+}3,\ 1{+}5,\ 1{+}6,\ 1{+}7,\ 2{+}5,\ 2{+}6,\ 2{+}7,\\
&3{+}6,\ 4{+}7,\ 5{+}6,\ 5{+}7.
\end{aligned}
$$
This fixed order governs computation and storage. It is not the same as the swatch-bar display order used in the final plots.

Then define the strategy landscape
$$
s^*(p,q)\in\{0,\dots,19\}
$$
to be the index of the first strategy in that list attaining the minimum.

Ties must be broken by first occurrence, with tolerance $10^{-10}$.  
Concretely: if a newly computed candidate is smaller than the current best by more than $10^{-10}$, replace the incumbent; otherwise keep the earlier one.

Because invalid two-strategies evaluate to $+\infty$, they simply never win the minimum.

This fixed order governs the computation. The color bar in the final plots uses a different fixed display order, specified in Section 9.

## 4. A computable global Lipschitz constant

### 4.1 The basic coordinate curves

The first five candidate coordinates come from
$$
0,\ 1-p,\ p,\ p(1-p),\ \tfrac12(1-p^2)
$$
and their $q$-counterparts.  
Each of these has derivative magnitude at most $1$ on $[0,1]$.

### 4.2 The certificate curves

For $\ell_\text{Sonetto},\ell_\text{Regulus},\ell_{\text{Regulus}^\mathsf T}$, use exactly the same Bernstein-coefficient method as in `09_verify_size_bound.py`.

For a stored certificate $\alpha$, define
$$
L(\alpha):=
14\max_{0\le t\le 13}\bigl|b_{t+1}(\alpha)-b_t(\alpha)\bigr|.
$$
Then set
$$
L_\text{Sonetto}:=\max_r L\bigl(\alpha^{(r)}_\text{Sonetto}\bigr),
$$

$$
L_\text{Regulus}:=\max_r L\bigl(\alpha^{(r)}_\text{Regulus}\bigr),
$$

$$
L_{\text{Regulus}^\mathsf T}:=\max_r L\bigl(\alpha^{(r)}_{\text{Regulus}^\mathsf T}\bigr).
$$

Because each envelope is the minimum of adjacent certificate curves, it inherits the same family-wide Lipschitz constant.

### 4.3 The effective bound used by the verifier

Define
$$
L_*:=\max\{1,\ L_\text{Sonetto},\ L_\text{Regulus},\ L_{\text{Regulus}^\mathsf T}\}.
$$

In the actual implementation, after those family constants are computed, the verifier should use `13` instead whenever all three of them are at most `13`; otherwise it should use the largest computed value.

Then every coordinate appearing in every point $P_s(p,q)$ is $L_*$-Lipschitz in its own variable.

For the eight one-strategies this is immediate:
$$
\delta_i(p,q)=\max\{x_i,y_i\}
$$
is the maximum of two coordinate functions, so it is also $L_*$-Lipschitz in each coordinate.

Because the two-strategy intersection test uses the nonzero tolerance $\mathrm{EPS}=10^{-11}$, the exact global Lipschitz statement from the idealized intersection model is no longer literally true.  
What remains true, and is sufficient for the proof, is the following safe comparison bound:

if
$$
|p-p'|\le \Delta_p,\qquad |q-q'|\le \Delta_q,
$$
then
$$
\delta(p',q')
\le
\delta(p,q)+L_*(\Delta_p+\Delta_q)+4\,\mathrm{EPS}.
$$

This is the inequality that should guide the certification argument.

Since
$$
4\,\mathrm{EPS}=4\cdot 10^{-11},
$$
it is enough in the implementation to **ignore** this additive term and instead use the slightly stricter numerical acceptance threshold
$$
0.3198989999.
$$
Indeed,
$$
0.3198989999+4\cdot 10^{-11}<0.319899.
$$

So the script should certify against `0.3198989999` internally, and then conclude the stated theorem
$$
\delta(p,q)\le 0.319899
\qquad\text{for all }(p,q)\in[0,1]^2.
$$

## 5. Recursive proof of $ \delta(p,q)\le 0.319899 $

The certifier works on rectangles
$$
R=[p_0,p_1]\times[q_0,q_1].
$$

### Step 1: sample the $101\times 101$ grid

For
$$
0\le a,b\le 100,
$$
set
$$
p_a=\frac{(100-a)p_0+ap_1}{100},
\qquad
q_b=\frac{(100-b)q_0+bq_1}{100}.
$$

Evaluate $\delta(p_a,q_b)$ on all $10201$ grid points.

During the whole recursive validation process, keep track of the global maximum sampled value
$$
M_\mathrm{sample}:=\max \delta(p_a,q_b)
$$
over all sampled grid points encountered, together with one location and one winning strategy where it is attained.

This information must be printed at the end of a successful run.

### Step 2: acceptance test

Every point of $R$ is within
$$
\frac{p_1-p_0}{200}
$$
in the $p$-coordinate and within
$$
\frac{q_1-q_0}{200}
$$
in the $q$-coordinate of some sampled grid point.

So define, for each sampled point,
$$
U_R(p_a,q_b):=
\delta(p_a,q_b)
+
L_*\left(\frac{p_1-p_0}{200}+\frac{q_1-q_0}{200}\right).
$$

Because of the bound from Section 4.3, the whole rectangle is certified once
$$
U_R(p_a,q_b)\le 0.3198989999
$$
holds for **every** sampled grid point.

This stricter threshold absorbs the omitted $4\,\mathrm{EPS}$ term.

### Step 3: immediate error reporting for an explicit violation

If during the validation process the script ever encounters a sampled grid point with
$$
\delta(p_a,q_b)>0.319899,
$$
then it has found an actual violation of the theorem, not merely a failure of the certification method.

In that case the program must raise an error immediately, in the same spirit as `09_verify_size_bound.py`, and print at least:

- the rectangle $R=[p_0,p_1]\times[q_0,q_1]$ in which it happened,
- the sampled point $(p_a,q_b)$,
- the sampled value $\delta(p_a,q_b)$,
- the winning strategy index,
- and the corresponding human-readable label (for example `2` or `0+7`).

### Step 4: subdivision

If the rectangle is not yet certified, subdivide it into two children:

- at even recursion depth, split in the $p$-direction at $(p_0+p_1)/2$;
- at odd recursion depth, split in the $q$-direction at $(q_0+q_1)/2$.

Then recurse on both children.

### Step 5: termination requirement

The program succeeds only if the recursion terminates with every leaf rectangle certified.

If not, it must print the first failing rectangle together with:

- its depth,

- its bounds $[p_0,p_1]\times[q_0,q_1]$,

- the largest sampled upper bound
  $$
  \max U_R(p_a,q_b)
  $$
  on that rectangle,

- the sampled point where this largest upper bound occurs,

- the sampled value $\delta(p_a,q_b)$ there,

- the winning strategy index in $\{0,\dots,19\}$,

- and the corresponding human-readable label, written as `i` for a one-strategy and `i+j` for a two-strategy.

At the end of a successful run, the script must print:

```
Starting Verification Algorithm V...
Target Bound: 0.3198989999
Using Threads: 16
Certificates loaded successfully.
------------------------------------------------------------
Verification Complete.
Time Elapsed: 2152.(2 decimal digits) seconds.
Total Recursive Calls: 496815
Total Grid Points Checked: 5068009815
Max Recursion Depth: 35
Max Delta_C Found: 0.(14 decimal digits)
Location of Max: (p=0.(14 decimal digits), q=(14 decimal digits))
RESULT: CERTIFIED (Max Found <= Target)
```

(Something like this.)

## 6. Heatmap data

The program evaluates the landscapes on two cell-centered grids.

### Full square

Domain:
$$
(p,q)\in[0,1]\times[0,1].
$$
Grid:
$$
p_i=\frac{i+1/2}{1600},
\qquad
q_j=\frac{j+1/2}{1600},
\qquad
0\le i,j<1600.
$$

Store:

- `best_full[i,j] = \delta(p_i,q_j)`
- `idx_full[i,j] = s^*(p_i,q_j)`

### Zoom square

Domain:
$$
(p,q)\in[0.2,0.4]\times[0.2,0.4].
$$
Grid:
$$
p_i=0.2+\frac{i+1/2}{1600}\cdot 0.2,
\qquad
q_j=0.2+\frac{j+1/2}{1600}\cdot 0.2,
\qquad
0\le i,j<1600.
$$

Store:

- `best_zoom[i,j] = \delta(p_i,q_j)`
- `idx_zoom[i,j] = s^*(p_i,q_j)`

So the four PNGs are built from the four arrays

- `best_full`
- `best_zoom`
- `idx_full`
- `idx_zoom`

## 7. Common figure geometry

All four figures use the same canvas:

- figure size `10 x 10` inches
- DPI `600`

So every saved image is exactly `6000 x 6000` pixels.

Use one fixed plotting box for the square heatmap,

```text
AX_BBOX = [0.08, 0.10, 0.80, 0.80]
```

and do **not** call `tight_layout()` or save with `bbox_inches='tight'`.

This keeps the full and zoom pictures perfectly aligned.

## 8. Value heatmaps

The first two figures visualize the scalar field $(p,q)\mapsto \delta(p,q)$.

### Array orientation

The stored arrays use index order `[i, j] = [p-index, q-index]`, but image plotting expects vertical index first.  
Therefore transpose before plotting:
$$
F^\mathsf T.
$$

So the final display has:

- horizontal axis = $p$
- vertical axis = $q$

### Colormap and normalization

Use the fixed continuous colormap built from

```text
#03051a, #0b3d91, #006d9c, #24b3b3, #9fffcf
```

with fixed normalization
$$
[0.31,\ 0.32].
$$

Values below `0.31` clip to the darkest color, and values above `0.32` clip to the lightest color.

Use

```text
interpolation = 'nearest'
```

so the picture shows the sampled grid directly.

### Full value heatmap

Save:

```text
11_figures/02_value_heatmap_full.png
```

Use:

- extent `[0, 1, 0, 1]`
- title `Cost landscape: δ(p,q) on [0,1] × [0,1]`
- x-label `p`
- y-label `q`
- ticks `0.0, 0.2, 0.4, 0.6, 0.8, 1.0`

Add a vertical colorbar at

```text
COLORBAR_BBOX = [0.90, 0.10, 0.03, 0.80]
```

with label `δ (cost)` and ticks

```text
0.31, 0.311, 0.312, ..., 0.32
```

### Zoomed value heatmap

Save:

```text
11_figures/03_value_heatmap_zoom.png
```

Use:

- extent `[0.25, 0.4, 0.25, 0.4]`
- title `Cost landscape: δ(p,q) on [0.25,0.4] × [0.25,0.4]`
- x-label `p`
- y-label `q`
- ticks `0.25, 0.28, 0.31, 0.34, 0.37, 0.40`

Use the **same** colormap, the **same** normalization, and the **same** colorbar placement as in the full plot. The zoomed heatmaps must actually be computed on the midpoint grid of `[0.25,0.4]×[0.25,0.4]`, namely with sample coordinates `0.25 + (m+1/2)·0.15/1600` in each axis for `m=0,1,\dots,1599`.

## 9. Strategy heatmaps

The last two figures visualize the categorical map $(p,q)\mapsto s^*(p,q)$.

### Strategy indexing and labels

Use the fixed list of 20 strategies from Section 3.

For display, use the following text labels:

- one-strategy `i` gets label

  ```text
  i
  ```

- two-strategy `(i,j)` gets label

  ```text
  i+j
  ```

### Color palette

Use the following fixed categorical palette:

```text
0+1 -> #F8F4EA
0+2 -> #F1E8D8
1+2 -> #E6DAC0
2   -> #F4BE18
5   -> #E69A00
2+5 -> #69B2D8
2+6 -> #2D6C9B
2+7 -> #143454
0+4 -> #F7C8C1
4+7 -> #F2A08F
1+7 -> #E97866
0+7 -> #D85A4B
5+7 -> #B83A31
1+5 -> #8F241D
1+3 -> #E1E9BD
3+6 -> #C8DF8B
0+6 -> #A9CE59
1+6 -> #7FB64A
5+6 -> #569632
0+5 -> #2F7428
```

The mapping

```text
strategy label -> color
```

must be fixed once and then reused unchanged.

### Full strategy heatmap

Save:

```text
11_figures/04_strategy_heatmap_full.png
```

Plot the RGB image obtained from `idx_full`, again transposed so that horizontal = $p$ and vertical = $q$.

Use:

- extent `[0, 1, 0, 1]`
- title `Strategy landscape: best strategy on [0,1] × [0,1]`
- x-label `p`
- y-label `q`
- ticks `0.0, 0.2, 0.4, 0.6, 0.8, 1.0`

Do **not** use a numeric colorbar.  
Instead place a one-column swatch bar on the right and draw all 20 colored swatches with their labels.

The swatches must appear in the following precise display order:
$$
\begin{aligned}
&0{+}1,\ 0{+}2,\ 1{+}2,\ 2,\ 5,\ 2{+}5,\ 2{+}6,\ 2{+}7,\\
&0{+}4,\ 4{+}7,\ 1{+}7,\ 0{+}7,\ 5{+}7,\ 1{+}5,\ 1{+}3,\\
&3{+}6,\ 0{+}6,\ 1{+}6,\ 5{+}6,\ 0{+}5.
\end{aligned}
$$
This display order is intentionally different from the computation order used for certification and storage.

### Zoomed strategy heatmap

Save:

```text
11_figures/05_strategy_heatmap_zoom.png
```

This is the same categorical plot, but on the zoom domain $[0.2,0.4]^2$, using `idx_zoom`.

Use:

- extent `[0.25, 0.4, 0.25, 0.4]`
- title `Strategy landscape: best strategy on [0.25,0.4] × [0.25,0.4]`
- x-label `p`
- y-label `q`
- ticks `0.25, 0.28, 0.31, 0.34, 0.37, 0.40`

Use the **same** category colors, the **same** label texts, and the **same** one-column swatch-bar order as in the full strategy plot.

The purpose of `05_strategy_heatmap_zoom.png` is to make the fine combinatorial partition inside the near-extremal window readable: it should be understood as a magnified view of the full strategy landscape, not as a separately recolored or separately indexed figure.

## 10. One-sentence summary of the script

The script evaluates all $20$ strategies—two one-point strategies and eighteen robustly filtered two-point strategies—on a recursive proof grid to certify
$$
\delta(p,q)\le 0.319899
\quad\text{for all }(p,q)\in[0,1]^2,
$$
prints the maximum sampled value encountered during the validation process, and then turns the resulting scalar field $\delta$ and categorical field $s^*$ into two value heatmaps and two strategy heatmaps with shared geometry and shared color conventions.
