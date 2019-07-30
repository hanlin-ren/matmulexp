# matmulexp
a calculator for [matrix multiplication exponent](https://en.wikipedia.org/wiki/Coppersmith%E2%80%93Winograd_algorithm), translated to Python from [Jan van den Brand's calculator in Javascript](https://people.kth.se/~janvdb/matrix.html). (Please refer to papers there for the correctness of this calculator.)

main code: `rectmm.py`. In which `omega(a,b,c)` is the exponent of multiplying a `n^a * n^b` matrix and a `n^b * n^c` matrix.

We hope that with Python libraries such as `scipy.optimize`, it'll become less painful to analyze rectangular-matrix-multiplication-based algorithms.

## demo
We include `rectmm` as well as some optimizers. Here we use optimizers in `scipy.optimize`.

    from scipy.optimize import fsolve, minimize
    from rectmm import *

Below is *the* exponent of matrix multiplication, currently `2.3728639`.

    w = omega(1, 1, 1)

We analyze the complexity of several published algorithms as demonstration:

[APSP in directed unweighted graph](https://arxiv.org/pdf/cs/0008011.pdf) is in `O(n^{2+mu})` time where `omega(1,mu,1) = 1+2mu`. The following code shows current `mu<=0.52853`:

    >>> fsolve(lambda mu: omega(1, mu, 1) - (1 + 2 * mu), 0.5)
    array([0.52852575])

[Dominance product](https://pdfs.semanticscholar.org/84cf/60d1bb1ab7e8b22734066119549a5bd003a8.pdf) is in `O(n^rho)` time where `rho = omega(1,4-rho,1)`. The following code shows current `rho<=2.65805`:

    >>> fsolve(lambda rho: omega(1, 4 - rho, 1) - rho, 2.4)
    array([2.65804365])

An [algorithm](http://drops.dagstuhl.de/opus/volltexte/2018/9048/pdf/LIPIcs-ICALP-2018-44.pdf) for all-pair nondecreasing paths run in `n^(t+omega) + n^(3-t+s) + n^(3-s-q) + n^(t+omega(1,1-s,1)+2q)` time. We use `scipy.optimize.minimize` to find the best tuple of parameters `(t,s,q)`:

    >>> minimize(lambda x: max(x[0]+w,3-x[0]+x[1],3-x[1]-x[2],x[0]+omega(1,1-x[1],1)+2*x[2]), [0.1,0.1,0.1], method = 'Nelder-Mead')
     final_simplex: (array([[0.397063  , 0.16686704, 0.06319079],
           [0.39709779, 0.16685612, 0.06325981],
           [0.39708387, 0.16687988, 0.06313003],
           [0.39715275, 0.1669564 , 0.06310318]]), array([2.76994217, 2.76996169, 2.76999009, 2.77001665]))
               fun: 2.7699421717571737
           message: 'Optimization terminated successfully.'
              nfev: 161
               nit: 88
            status: 0
           success: True
                 x: array([0.397063  , 0.16686704, 0.06319079])
That is, when `t=0.397063, s=0.16686704` and `q=0.06319079`, this algorithm runs in `O(n^2.76995)` time.

Theorem 4.2 of [this paper](https://arxiv.org/pdf/1905.05067.pdf) states an operation of a data structure in `n^(e2+e1) + n^(omega(1,e1,e2)-e1) + n^(omega(1,1,e2)-e2)` time, which we find the optimal values:

    >>> minimize(lambda x: max(x[1]+x[0], omega(1,x[0],x[1])-x[0], omega(1,1,x[1])-x[1]), [0.2,0.2], method = 'Nelder-Mead')
     final_simplex: (array([[0.55125019, 0.85560732],
           [0.55127205, 0.85559352],
           [0.55116682, 0.85559014]]), array([1.40688456, 1.40688808, 1.40688894]))
               fun: 1.4068845552536864
           message: 'Optimization terminated successfully.'
              nfev: 122
               nit: 66
            status: 0
           success: True
                 x: array([0.55125019, 0.85560732])
That is, when `e1=0.55125019` and `e2=0.85560732`, the operation is in `O(n^1.40689)` time.
