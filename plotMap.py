import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

name1 = "Matriz da Densidade - refinada.txt"
name2 = "Matriz da Velocidade em X - refinada.txt"
name3 = "Matriz da Velocidade em Y - refinada.txt"

#name1 = "Matriz da Densidade.txt"
#name2 = "Matriz da Velocidade em X.txt"
#name3 = "Matriz da Velocidade em Y.txt"

data1 = np.loadtxt(name1)
data2 = np.loadtxt(name2)
data3 = np.loadtxt(name3)

media = sum(data1)
#print(sum(media)/(2000.0*80.0))
print(sum(media)/(1000.0*40.0))

fig1 = plt.figure()

cb1 = plt.imshow(np.transpose(data1))
plt.colorbar(cb1)
plt.title("Mapa de densidade")

plt.show()

fig2 = plt.figure()

cb2 = plt.imshow(np.transpose(data2))
plt.colorbar(cb2)
plt.title("Mapa de velocidade em X")

plt.show()

fig3 = plt.figure()

cb3 = plt.imshow(np.transpose(data3))
plt.colorbar(cb3)
plt.title("Mapa de velocidade em Y")

plt.show()

fig4 = plt.figure()
vec = (data2[900][:])
plt.plot(vec)
plt.title("Escoamento ao longo da linha")
#fig4.savefig("Escoamento ao longo da linha.png")

plt.show()
