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
plt.legend(namelist,loc="upper right")
plt.title('Perfil de velocidade em X')
plt.grid()
plt.show()
