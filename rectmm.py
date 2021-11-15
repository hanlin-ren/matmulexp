from scipy.optimize import fsolve, minimize

Ks = [0, 0.31389, 0.32, 0.33, 0.34, 0.35, 0.40, 0.45, 0.50, 0.5286, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0]
Omegas = [2, 2, 2.000064, 2.000448, 2.001118, 2.001957, 2.010314, 2.024801, 2.044183, 2.057085, 2.067488, 2.093981, 2.123097, 2.154399, 2.187543, 2.222256, 2.258317, 2.295544, 2.333789, 2.3728596, 2.453481, 2.536550, 2.621644, 2.708400, 2.796537, 3.021591, 3.251640, 3.721503, 4.199712, 5.171210, 6.157233];

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
	print w
	print 2+mu[0]
	print rho[0]
	print apnp.fun, apnp.x, apnp.success
	print matinv.fun, matinv.x, matinv.success

