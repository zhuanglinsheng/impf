# Polynomial Approximation

In an implicit function approximation problem $F(y(x), x)=0$ where $F:\mathbb{R}^m\times\mathbb{R}^n\to\mathbb{R}^m$ and $y:\mathbb{R}^n\to\mathbb{R}^m$, it is natural to approximate the implicit function $y(x_1, x_2, ..., x_n)$ with polynomials. As an example, a third order polynomial approximation is 

$$
y(x_1, x_2, ..., x_n) \approx c_0 + \sum^n_{i=1}c_ix_i + \sum_{i=1}^n\sum_{j\ge i}^nc_{i,j}x_ix_j + \sum_{i=1}^n\sum_{j\ge i}^n\sum_{k\ge j}^nc_{i,j,k}x_ix_jx_k
$$

We consider the approximation on a small cube $[l_i, u_i]^n_{i=1}$, where $l_i$ and $u_i$ are the $i$-th lower and upper bound of $x_i$, in other words, $l_i \le x_i \le u_i$ for $i=1, 2, ..., n$.

## Sampling the Polynomial Bases

The function `impf_polybases_sampling_on_cube` will sample the polynomial bases on a cube, whose prototype is 

```c
impf_matrix * impf_polybases_sampling_on_cube(
    const unsigned int m,           /* Number of implicit functions */
    const unsigned int n,           /* Number of the parameters of each i.f. */
    const unsigned int poly_order,  /* Order of the polynomial */
    const double * ls,              /* Lower bounds of the cube */
    const double * us,              /* Upper bounds of the cube */
);
```

Based on the above inputs, the function returns a matrix of sampled polynomial bases 
$X$, which is a matrix of the shape $N\times D$, where $N$ is the sample size (number of rows) and 

$$
D \equiv 1+{n \choose 1}+{n+1\choose 2}+...+{n-1+p\choose p}
$$

is the number of polynomial bases, and $p$ is the order of the polynomial `poly_order`. The matrix row number (sample size) $N$ is greater or equal to the column number $D$. 

Take the 2nd order polynomial as an example, the sampled base is:

|   | 1 | $x_1$      | ... | $x_n$      | $x_1x_1$       | ... | $x_nx_n$       | 
|:--|:-:|:----------:|:---:|:----------:|:--------------:|:---:|:--------------:|
| 1 | 1 | $l_1$      |     | $l_n$      | $l_1^2$        |     | $l_n^2$        |
| 2 | 1 | $l_1+h_1$  |     | $l_n+h_n$  | $(l_1+h_1)^2$  |     | $(l_n+h_n)^2$  |
| 3 | 1 | $l_1+2h_1$ |     | $l_n+2h_n$ | $(l_1+2h_1)^2$ |     | $(l_n+2h_n)^2$ |
|   ...                                                                         |
| N | 1 | $u_1$      |     | $u_n$      | $u_1^2$        |     | $u_n^2$        |

## Selecting the Coefficients

The polynomial coefficients, denoted by $c$, is a vector of the length $D$. The form of $c$ is like (2nd order polynomial as an example)

$$
c = (c_0, c_1, ..., c_n, c_{11}, ..., c_{nn})
$$

We need to select the optimal coefficients so that the approximated implicit function zeros the outer nonlinear function $F$ as much as possible. Commonly, people choose the Euclidean norm to measure the distance of a vector to the origin. Hence, the approximation problem can be described as 

$$
\min_c\ \left\|\begin{pmatrix}
F(cX_1, X_1)\\
F(cX_2, X_2)\\
...\\
F(cX_N, X_N)
\end{pmatrix}\right\|
$$

where $X_i$ is the $i$-th row of the matrix $X$. If $F$ is a linear function, then the above problem is in essence a least square problem. Generally, problems of the above form can be solved by nonlinear least squared algorithms. 

TO DO... (nonlinear least square)

## Error Estimation and Meshing

Suppose the shape of a cube $[l_i,u_i]_{i=1}^n$ is $h \equiv (h_1,...,h_n)$, i.e., $h_i = u_i-l_i$ for $i=1, 2, ..., n$. Then, the polynomial approximation error, i.e., the Taylor residual term, is $\mathcal{O}(\|h\|^{p+1})$. More specifically, the Lagrangian residual term is upper bounded by 

$$
\sup_{0\le\alpha\le h} D^{(p+1)}y(l+\alpha)(\alpha/\|\alpha\|) \cdot \frac{1}{(p+1)!} 
$$

where $l = (l_1, l_2, ..., l_n)$ and $D^{(p+1)}y(l+\alpha)(\alpha)$ is the $(p+1)$-th order derivative of $y$ at $l+\alpha$ in the direction $\alpha/\|\alpha\|$, and

$$
D^{(p+1)}y(l+\alpha)(\alpha/\|\alpha\|) = \sum_{k_1+...+k_n=p+1}\frac{\partial^{p+1}y(l+\alpha)}{\partial x_1^{k_1}...\partial x_n^{k_n}}\frac{\alpha_1^{k_1}...\alpha_n^{k_n}}{\|\alpha\|^{p+1}}
$$

