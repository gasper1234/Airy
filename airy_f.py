import decimal
from decimal import Decimal
from scipy.special import airy as airy
from scipy.special import gamma as gamma
import numpy as np
import matplotlib.pyplot as plt

#majhne

g_1_3 = Decimal(gamma(1/3))
g_2_3 = Decimal(gamma(2/3))
alpha = Decimal(0.355028053887817239)
beta = Decimal(0.258819403792806798)

koef_f = ['s' for _ in range(500)]
koef_f[0] = Decimal('1')

def f(x, n):
	sum = 0
	for i in range(n):
		if koef_f[i] == 's':
			koef_f[i] = koef_f[i-1] * Decimal((1/3+i-1) * 3 / ((3*i) * (3*i-1) * (3*i-2)))
		sum += Decimal(x ** (3*i)) * koef_f[i]
	return sum

koef_g = ['s' for _ in range(500)]
koef_g[0] = Decimal('1')

def g(x, n):
	sum = 0
	for i in range(n):
		if koef_g[i] == 's':
			koef_g[i] = koef_g[i-1] * Decimal((2/3+i-1) * 3 / ((3*i+1) * (3*i) * (3*i-1)))
		sum += Decimal(x ** (3*i+1)) * koef_g[i]
	return sum

def Ai_T(x, n):
	return alpha*f(x,n) - beta*g(x,n)

#velike

def zeta(x):
	return 2/3 * abs(x)**(3/2)

koef_L = ['s' for _ in range(1000)]
koef_L[0] = Decimal('1')

def L(z):
	sum = Decimal('0')
	p_term = np.inf
	for i in range(1000):
		if koef_L[i] == 's':
			koef_L[i] = koef_L[i-1] * ((3*i-Decimal('0.5')) * (3*i+Decimal('0.5')-2) * (3*(i-1)+Decimal('0.5'))) / (54 * i * (i-Decimal('0.5')))
		c_term = koef_L[i] * Decimal(z ** (-1*i))
		if c_term > p_term:
			break
		sum += c_term
		p_term = c_term
	return sum

koef_P = ['s' for _ in range(1000)]
koef_P[0] = Decimal('1')

def P(z):
	sum = Decimal('0')
	p_term = np.inf
	for i in range(1000):
		if koef_P[i] == 's':
			gamma_fa = 1
			for j in range(1, 7):
				gamma_fa *= (6*j+Decimal('0.5')-i)
			koef_P[i] = koef_P[i-1] * Decimal('-1') * gamma_fa / (54 ** 2 * 2*i * (2*i-1) * (2*i-Decimal('0.5')) * (2*(i-1)+Decimal('0.5')))
		c_term = koef_P[i] * Decimal(z ** (-2*i))
		if abs(c_term) > abs(p_term):
			break
		sum += c_term
		p_term = c_term
	return sum

koef_Q = ['s' for _ in range(1000)]
koef_Q[0] = Decimal('1')/Decimal('54')*Decimal('2.5')*Decimal('1.5')

def Q(z):
	sum = Decimal('0')
	p_term = np.inf
	for i in range(1000):
		if koef_Q[i] == 's':
			gamma_fa = 1
			for j in range(1, 7):
				gamma_fa *= (6*j+Decimal('0.5')+3-i)
			koef_Q[i] = koef_Q[i-1] * Decimal('-1') * gamma_fa / (54 ** 2 * 2*i * (2*i+1) * (2*i+Decimal('0.5')) * (2*i-Decimal('0.5')))#todo preveri koef!!!!!!!!!!!!!!
		c_term = koef_Q[i] * Decimal(z ** (-2*i-1))
		if abs(c_term) > abs(p_term):
			break
		sum += c_term
		p_term = c_term
	return sum

def Ai_A(x):
	if x > 0:
		return Decimal(np.exp(-1*zeta(x)) / (2*np.sqrt(np.pi)*x**0.25)) * L(zeta(x))
	else:
		return Decimal(1 / (2*np.sqrt(np.pi) * (-1*x)**(1/4))) * ( Decimal(np.sin(zeta(x)-np.pi/4)) * Q(zeta(x))  +  Decimal(np.cos(zeta(x)-np.pi/4)) * P(zeta(x)))

data = [i for i in range(-30, 30, 7)]

for j in range(len(data)):
	print(data[j], 'true', airy(data[j])[0], 'asimpt', Ai_A(data[j]), 'dif', abs(Decimal(airy(data[j])[0])-Ai_A(data[j])), 'rel dif', abs(Decimal(airy(data[j])[0])-Ai_A(data[j]))/Decimal(airy(data[j])[0]))

#print(koef_L)
'''
N = 400
for j in range(len(data)):
		print(data[j], 'taylor', abs(Decimal(airy(data[j])[0])-Ai_T(data[j], N)), 'rel dif', abs(Decimal(airy(data[j])[0])-Ai_T(data[j], N))/Decimal(airy(data[j])[0]))
		print('asimpt', abs(Decimal(airy(data[j])[0])-Ai_A(data[j])), 'rel dif', abs(Decimal(airy(data[j])[0])-Ai_A(data[j]))/Decimal(airy(data[j])[0]))
'''
#print(koef_Q)
#print(koef_P)