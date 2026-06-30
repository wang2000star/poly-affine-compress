## Minimizing the number of ANF monomials via affine input transformations

### Problem statement

Let $f: \mathbb{F}_2^n \to \mathbb{F}_2$ be a Boolean function given in Algebraic Normal Form (ANF):

$$f(x) = \bigoplus_{s \in \{0,1\}^n} a_s \cdot x^s, \quad a_s \in \mathbb{F}_2$$

where $x^s = \prod_{i: s_i=1} x_i$ and multiplication is in the Boolean ring ($x_i^2 = x_i$). The complexity measure $T(f)$ is the number of monomials with $a_s = 1$.

We consider **pure algebraic substitution**: find an affine transformation $z = Mx \oplus b$ with $M \in \mathbb{F}_2^{m \times n}$, $b \in \mathbb{F}_2^m$ (no restriction on $m$, $M$ need not be square or full-rank) and a Boolean function $g: \mathbb{F}_2^m \to \mathbb{F}_2$ such that

$$f(x) = g(Mx \oplus b)$$

holds as a formal polynomial identity in $\mathbb{F}_2[x]/(x_i^2 - x_i)$.

The goal is to minimize $T(g)$ relative to $T(f)$. Note that $m$ can differ from $n$, and $\deg(g) = \deg(f)$ for all admissible $M,b$.

### What we know

For quadratic functions ($\deg(f) = 2$), the problem is well-understood: the rank $r$ of the associated alternating bilinear form $B = A + A^{\mathsf{T}}$ (where $A$ is the strictly upper-triangular coefficient matrix) is an affine invariant, and $T_{\min} \le \lfloor r/2 \rfloor + \delta$ for some $\delta \in \{0,1\}$.

For $\deg(f) \ge 3$, no such complete characterization exists.

### Heuristic approaches tested

1. **Greedy XOR merge**: Iteratively replace $x_i \to x_i \oplus x_j$ when it reduces $T(f)$. This is equivalent to elementary row operations on $M$. It works remarkably well for functions with low intrinsic dimension (e.g., $f(x) = G(Ax)$ where $G$ depends on only $k \ll n$ linear forms), achieving $>99\%$ reduction in $T$. However, it fails for random cubic functions.

2. **Variable complement search**: Try $z_i = x_i$ or $z_i = x_i \oplus 1$ for each variable. This gives modest improvements (5–35%) for small structured functions but is negligible for large random functions.

3. **Walsh-Hadamard spectrum guidance**: Directions with large spectral coefficients are candidates for dimension-reducing transformations. This helps for some structured functions but rarely finds improvements the greedy merge misses.

### Open questions

1. **Polynomial-time approximation**: Is there a polynomial-time algorithm (in $n$ and $T(f)$) that, for any $f$, finds $M,b$ achieving $T(g)$ within a factor of $O(\text{poly}(n))$ of the optimum?

2. **Intrinsic complexity measure**: For $\deg(f) \ge 3$, does there exist an efficiently computable invariant (analogous to the quadratic rank) that characterizes the minimum possible $T(g)$ under affine input transformations?

3. **Random functions**: For a random $n$-variable Boolean function (where $\mathbb{E}[T(f)] = 2^{n-1}$), what is the expected minimum $T(g)$ achievable via affine input transformations? Is it $\Theta(2^{n-1})$ or significantly smaller?

4. **Connection to coding theory**: The ANF coefficient vector $a \in \mathbb{F}_2^{2^n}$ lives in an orbit under the action of the affine group $\text{AGL}(n,\mathbb{F}_2)$. $T(f)$ is the Hamming weight of $a$. The problem asks for the minimum weight in this orbit — can coding-theoretic bounds (e.g., Delsarte's linear programming bound) be applied?

Any pointers to relevant literature or algorithmic ideas would be greatly appreciated.
