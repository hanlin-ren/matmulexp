from scipy.optimize import fsolve, minimize

Ks =     [0, 0.321334,     0.33,     0.34,     0.35,     0.40,     0.45,     0.50]
Omegas = [2,        2, 2.000100, 2.000600, 2.001363, 2.009541, 2.023788, 2.042994]
Ks +=       [0.527661,     0.55,     0.60,     0.65,     0.70,     0.75,     0.80,     0.85]
Omegas +=   [2.055322, 2.066134, 2.092631, 2.121734, 2.153048, 2.186210, 2.220929, 2.256984]
Ks +=           [0.90,     0.95,     1.00,     1.10,     1.20,     1.30,     1.40,     1.50]
Omegas +=   [2.294209, 2.332440, 2.371552, 2.452056, 2.535063, 2.621644, 2.708400, 2.794941]
Ks +=           [1.75,     2.00,     2.50,     3.00,     4.00,     5.00]
Omegas +=   [3.021591, 3.250385, 3.720468, 4.198809, 5.171210, 6.156708]

def interpolated_omega(k) :
	for i in range(1, len(Ks)) :
		if Ks[i] >= k :
			intervalLength = Ks[i] - Ks[i - 1]
			offset = (k - Ks[i - 1]) / intervalLength
			return Omegas[i - 1] * (1 - offset) + Omegas[i] * offset
	return k - Ks[-1] + Omegas[-1]

def omega(a, b, c) :
	c, b, a = sorted([1.0 * a, 1.0 * b, 1.0 * c])
	if c < 0.0001 : return a + b + c
	if a - b < 0.0001 : b = a
	if b - c < 0.0001 : c = b
	if a == b :
		return a * interpolated_omega(c / a)
	if b == c :
		return b * interpolated_omega(a / b)
	return min((a - b) + b * interpolated_omega(c / b), (b - c) + c * interpolated_omega(a / c))

if __name__ == '__main__' :
	w = omega(1, 1, 1)
	mu = fsolve(lambda mu: omega(1, mu, 1) - (1 + 2 * mu), 0.5)
	rho = fsolve(lambda rho: omega(1, 4 - rho, 1) - rho, 2.4)
	apnp = minimize(lambda x: max(x[0]+w,3-x[0]+x[1],3-x[1]-x[2],x[0]+omega(1,1-x[1],1)+2*x[2]), [0.1,0.1,0.1], method = 'Nelder-Mead')
	matinv = minimize(lambda x: max(x[1]+x[0], omega(1,x[0],x[1])-x[0], omega(1,1,x[1])-x[1]), [0.2,0.2], method = 'Nelder-Mead')
	print(w)
	print(2+mu[0])
	print(rho[0])
	print(apnp.fun, apnp.x, apnp.success)
	print(matinv.fun, matinv.x, matinv.success)

