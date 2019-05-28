import numpy as np

a = 1e14
# s = 1e12
frac = 1e-4
EPS = 1e4

for s in np.linspace(0.01 * a, 0.99 * a):
	print("s = ", s)
	ok = False
	MAXIT = 500
	it = 0
	while ok != True and it < MAXIT:
		dx = frac * a
		tmp = ((s + dx) - s) / s
		print("fractional change = ", tmp)

		if tmp < 1e-4:
			frac *= 5
		else:
			ok = True

		if dx < EPS:
			dx = EPS
			break
		it+=1
print(s, dx)
