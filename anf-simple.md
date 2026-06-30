$x\in \mathbb{F}_2^n$，$f(x)$是$n$元布尔函数，$T(f)$是$f(x)$的ANF形式的单项式个数。

对于一个已知的$f(x)$，显然存在$M\in \mathbb{F}_2^{m \times n},b\in\mathbb{F}_2^m$，令$z=Mx\oplus b\in \mathbb{F}_2^m$，使得$g(z)=g(Mx\oplus b)=f(x)$。

例如，$f(x)=x_0x_1+x_0x_2$，令$z_0=x_0,z_1=x_1\oplus x_2$，则$g(z)=z_0z_1$

问题1：$g(z)$的代数次数与$f(x)$的是否相同？

问题2：如何构造出$M,b$，使得$\frac{T(g)}{T(f)}$尽可能地小 ？有没有高效的算法？