import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

name1 = "Matriz da Velocidade em X - refinada.txt"
name2 = "Matriz da Velocidade em X.txt"

data1 = np.loadtxt(name1)
data2 = np.loadtxt(name2)

print(len(data1[900][:]))
print(len(data2[900][:]))

vec1 = (data1[900][:])
vec2 = (data2[900][:])

namelist = ['Refinado','Grosseiro']

plt.plot(vec1,c='red')
plt.plot(vec2,c='blue')
plt.legend(namelist,loc="lower right")
plt.title('Perfil de velocidade em X')
plt.grid()
plt.show()

###################################################################
#Comparando as duas velocidades
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

namelist = ['Analitico','Numérico']

plt.plot(vec,'-',c='red')
plt.plot(vec1,'*',c='blue')
plt.title('Comparação entre a velocidade analítica e a numérica para 80 pontos - Refinada')
plt.legend(namelist,loc="upper right")
plt.grid()
plt.show()

###############################################################################

fig5 = plt.figure()

erro = np.zeros(len(vec1))
x    = np.arange(len(vec1))

#não sei se está certinho assim
for i in range(len(vec1)):
    erro[i] = np.absolute(vec[i] - vec1[i])

plt.bar(x, erro, align='center', alpha=0.5)
plt.xticks(x, x)
plt.ylabel('Diferença')
plt.xlabel('Posição')
plt.title('Erro em diferença absoluta')
plt.show()
