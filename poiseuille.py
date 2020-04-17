import numpy as np
import matplotlib.pyplot as plt

rd = 2
raio = 40

dy = 1.0
dx = 1.0
dt = 1.0

tau = 0.56000000

ni = (2*tau - 1)*(dx*dx) / (6*dt)

v_max = 0.003

dP = v_max * 4 * ni * dx / (raio*raio)

vec = np.zeros([(rd*raio) + 1])

print('ni: ',ni)
print('dP: ', dP)

for i in range(0,len(vec),1):
    vec[i] = (dP*(raio*raio - np.abs((raio - i)*(raio - i)))) / (4*ni*dx)

plt.plot(vec,'-',c='red')
plt.grid()
plt.title('Escoamento de Poiseuille')
plt.show()
