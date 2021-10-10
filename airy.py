from airy_f import *
'''
N = 2
points = [10**i for i in range(-N-1, N)]
x = np.array(points)
print(x)
y_true = airy(x)[0]
y_T = abs(Ai_T(x)-y_true)
y_A = abs(Ai_A(x)-y_true)
n = 2
y_A_2 = abs(Ai_A(x)-y_true)
print(y_T, y_A)

plt.plot(x, y_T, 'x', label='taylor')
plt.plot(x, y_A, '+', label='asimptotska_1')
plt.plot(x, y_A_2, 'o',label='asimptotska_2')
#plt.plot(x, y_true, label='true')
plt.legend()
plt.xscale('log')
plt.show()
'''
