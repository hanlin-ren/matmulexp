# matmulexp
a calculator for [matrix multiplication exponent](https://en.wikipedia.org/wiki/Coppersmith%E2%80%93Winograd_algorithm), translated to Python from [Jan van den Brand's calculator in Javascript](https://people.kth.se/~janvdb/matrix.html).

**UPD:** Bounds in this program has been updated to reflect recent progress by [Duan, Wu, and Zhou](https://arxiv.org/abs/2210.10173) and [Le Gall](https://arxiv.org/abs/2307.06535). We also included some old bounds by [Le Gall and Urrutia](https://arxiv.org/abs/1708.05622). We thank Martin Schirneck for informing us about the updates.

main code: `rectmm.py`. In which `omega(a,b,c)` is the exponent of multiplying a `n^a * n^b` matrix and a `n^b * n^c` matrix.

We hope that with Python libraries such as `scipy.optimize`, it'll become less painful to analyze rectangular-matrix-multiplication-based algorithms.

## demo
We include `rectmm` as well as some optimizers. Here we use optimizers in `scipy.optimize`.

    from scipy.optimize import fsolve, minimize
    from rectmm import *

Below is *the* exponent of matrix multiplication, currently `2.371866`.

    w = omega(1, 1, 1)

We analyze the complexity of several published algorithms as demonstration:

[APSP in directed unweighted graph](https://arxiv.org/pdf/cs/0008011.pdf) is in `O(n^{2+mu})` time where `omega(1,mu,1) = 1+2mu`. The following code shows current `mu<=0.52853`:

    >>> fsolve(lambda mu: omega(1, mu, 1) - (1 + 2 * mu), 0.5)
    array([0.52852557])

[Dominance product](https://pdfs.semanticscholar.org/84cf/60d1bb1ab7e8b22734066119549a5bd003a8.pdf) is in `O(n^rho)` time where `rho = omega(1,4-rho,1)`. The following code shows current `rho<=2.65805`:

    >>> fsolve(lambda rho: omega(1, 4 - rho, 1) - rho, 2.4)
    array([2.65804365])

An [algorithm](http://drops.dagstuhl.de/opus/volltexte/2018/9048/pdf/LIPIcs-ICALP-2018-44.pdf) for all-pair nondecreasing paths run in `n^(t+omega) + n^(3-t+s) + n^(3-s-q) + n^(t+omega(1,1-s,1)+2q)` time. We use `scipy.optimize.minimize` to find the best tuple of parameters `(t,s,q)`:

    >>> minimize(lambda x: max(x[0]+w,3-x[0]+x[1],3-x[1]-x[2],x[0]+omega(1,1-x[1],1)+2*x[2]), [0.1,0.1,0.1], method = 'Nelder-Mead')
           message: Optimization terminated successfully.
           success: True
            status: 0
               fun: 2.769618825006436
                 x: [ 3.978e-01  1.674e-01  6.305e-02]
               nit: 97
              nfev: 171
     final_simplex: (array([[ 3.978e-01,  1.674e-01,  6.305e-02],
                           [ 3.977e-01,  1.674e-01,  6.306e-02],
                           [ 3.978e-01,  1.674e-01,  6.298e-02],
                           [ 3.978e-01,  1.674e-01,  6.303e-02]]), array([ 2.770e+00,  2.770e+00,  2.770e+00,  2.770e+00]))

That is, when `t=0.3978, s=0.1674` and `q=0.06305`, this algorithm runs in `O(n^2.76962)` time.

Theorem 4.2 of [this paper](https://arxiv.org/pdf/1905.05067.pdf) states an operation of a data structure in `n^(e2+e1) + n^(omega(1,e1,e2)-e1) + n^(omega(1,1,e2)-e2)` time, which we find the optimal values:

    >>> minimize(lambda x: max(x[1]+x[0], omega(1,x[0],x[1])-x[0], omega(1,1,x[1])-x[1]), [0.2,0.2], method = 'Nelder-Mead')
           message: Optimization terminated successfully.
           success: True
            status: 0
               fun: 1.406852908543138
                 x: [ 5.512e-01  8.556e-01]
               nit: 67
              nfev: 123
     final_simplex: (array([[ 5.512e-01,  8.556e-01],
                           [ 5.512e-01,  8.556e-01],
                           [ 5.512e-01,  8.556e-01]]), array([ 1.407e+00,  1.407e+00,  1.407e+00]))

That is, when `e1=0.5512` and `e2=0.8556`, the operation is in `O(n^1.40686)` time.
