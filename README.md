# matmulexp
a calculator for [matrix multiplication exponent](https://en.wikipedia.org/wiki/Coppersmith%E2%80%93Winograd_algorithm), translated to Python from [Jan van den Brand's calculator in Javascript](https://people.kth.se/~janvdb/matrix.html).

**UPD:** Bounds in this program has been updated to reflect recent progress by [Alman, Duan, Williams, Xu, Xu, and Zhou](https://arxiv.org/abs/2404.16349).

main code: `rectmm.py`. In which `omega(a,b,c)` is the exponent of multiplying a `n^a * n^b` matrix and a `n^b * n^c` matrix.

We hope that with Python libraries such as `scipy.optimize`, it'll become less painful to analyze rectangular-matrix-multiplication-based algorithms.

## demo
We include `rectmm` as well as some optimizers. Here we use optimizers in `scipy.optimize`.

    from scipy.optimize import fsolve, minimize
    from rectmm import *

Below is *the* exponent of matrix multiplication, currently `2.371339`.

    w = omega(1, 1, 1)

We analyze the complexity of several published algorithms as demonstration:

[APSP in directed unweighted graph](https://arxiv.org/pdf/cs/0008011.pdf) is in `O(n^{2+mu})` time where `omega(1,mu,1) = 1+2mu`. The following code shows current `mu<=0.5274994`:

    >>> fsolve(lambda mu: omega(1, mu, 1) - (1 + 2 * mu), 0.5)
    array([0.52749936])

[Dominance product](https://pdfs.semanticscholar.org/84cf/60d1bb1ab7e8b22734066119549a5bd003a8.pdf) is in `O(n^rho)` time where `rho = omega(1,4-rho,1)`. The following code shows current `rho<=2.65805`:

    >>> fsolve(lambda rho: omega(1, 4 - rho, 1) - rho, 2.4)
    array([2.65804365])

An [algorithm](http://drops.dagstuhl.de/opus/volltexte/2018/9048/pdf/LIPIcs-ICALP-2018-44.pdf) for all-pair nondecreasing paths run in `n^(t+omega) + n^(3-t+s) + n^(3-s-q) + n^(t+omega(1,1-s,1)+2q)` time. We use `scipy.optimize.minimize` to find the best tuple of parameters `(t,s,q)`:

    >>> apnp = minimize(lambda x: max(x[0]+w,3-x[0]+x[1],3-x[1]-x[2],x[0]+omega(1,1-x[1],1)+2*x[2]), [0.1,0.1,0.1], method = 'Nelder-Mead')
    >>> print(apnp.fun, apnp.x, apnp.success)
    2.7693405337969965 [0.39796647 0.167307   0.06347337] True

That is, when `t=0.39796647, s=0.167307` and `q=0.06347337`, this algorithm runs in `O(n^2.769341)` time.

Theorem 4.2 of [this paper](https://arxiv.org/pdf/1905.05067.pdf) states an operation of a data structure in `n^(e2+e1) + n^(omega(1,e1,e2)-e1) + n^(omega(1,1,e2)-e2)` time, which we find the optimal values:

    >>> matinv = minimize(lambda x: max(x[1]+x[0], omega(1,x[0],x[1])-x[0], omega(1,1,x[1])-x[1]), [0.2,0.2], method = 'Nelder-Mead')
    >>> print(matinv.fun, matinv.x, matinv.success)
    1.4056447265523806 [0.55041856 0.85513407] True

That is, when `e1=0.55041856` and `e2=0.85513407`, the operation is in `O(n^1.405645)` time.
