## Minimizing ANF monomials via affine input substitution

Let $f: \mathbb{F}_2^n \to \mathbb{F}_2$ be a Boolean function in Algebraic Normal Form:

$$f(x) = \bigoplus_{s \in \{0,1\}^n} a_s \cdot x^s, \quad a_s \in \mathbb{F}_2$$

where $x^s = \prod_{i: s_i=1} x_i$ in the Boolean ring ($x_i^2 = x_i$). $T(f) = |\{s : a_s = 1\}|$ is the number of monomials.

I want to find an affine transformation $z = Mx \oplus b$ with $M \in \mathbb{F}_2^{m \times n}$, $b \in \mathbb{F}_2^m$ ($m$ can differ from $n$, $M$ need not be square or full rank) and a Boolean function $g: \mathbb{F}_2^m \to \mathbb{F}_2$ such that

$$f(x) = g(Mx \oplus b)$$

holds as a **formal polynomial identity** in $\mathbb{F}_2[x]/(x_i^2 - x_i)$, and $T(g) < T(f)$.

Note this is a pure algebraic substitution — no function evaluation, no information loss, and $\deg(g) = \deg(f)$ for any $M,b$.

For quadratic functions ($\deg f = 2$), the problem is well understood: the rank of $B = A + A^\mathsf{T}$ (the alternating bilinear form) determines the minimum $T(g)$. For $\deg f \ge 3$, I'm not aware of any similar invariant or efficient algorithm.

**Question**: Is there any known result or algorithm for this problem when $\deg f \ge 3$? Any pointers to relevant literature would be appreciated.
