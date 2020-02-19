import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#arquivos a serem usados
name1 = "results/direções grosseira.txt"
name2 = "results/direções Pós-Processamento.txt"


#abre os arquivos na forma de uma lista
dataF = np.loadtxt(name2)
dataC = np.loadtxt(name1)

#razao de divisao
rd = 2

#tamanho da matriz a ser alocada (igual ao da matriz grosseira vezes o fator de divisão em Y)
tamX = 1000
tamY = 20

#criando a nova matriz
rho = np.zeros([tamX,rd*tamY])
ux  = np.zeros([tamX,rd*tamY])
uy  = np.zeros([tamX,rd*tamY])

ux_fine = np.zeros([tamY])
#-----------------------------------------------------------------CALCULANDO A DENSIDADE-----------------------------------------------------#
rhoMed = 0
it = 0
for i in range(0,tamX,1):
    for j in range(0,tamY,1):
        for c in range(9):
            rho[i][j] = rho[i][j] + dataC[it][c]
        rhoMed = rhoMed + rho[i][j]
        it = it+1

it = 0
for i in range(0,tamX,1):
    for j in range(0,tamY,1):
        for c in range(9):
            rho[i][tamY + j] = rho[i][tamY + j] + dataF[it][c]
        rhoMed = rhoMed + rho[i][j]
        it = it+1

#--------------------------------------------------------------CALCULANDO A VELOCIDADE EM X-----------------------------------------------------#
it = 0
for i in range(0,tamX,1):
    for j in range(0,tamY,1):
        ux[i][tamY + j] = dataC[it][1] + dataC[it][5] + dataC[it][8] - (dataC[it][3] + dataC[it][6] + dataC[it][7])
        ux[i][tamY + j] = ux[i][tamY + j]/(rho[i][tamY + j] + 0.0000001)
        it = it+1

it = 0
for i in range(0,tamX,1):
    for j in range(0,tamY,1):
        ux[i][j] = dataF[it][1] + dataF[it][5] + dataF[it][8] - (dataF[it][3] + dataF[it][6] + dataF[it][7])
        ux[i][j] = ux[i][j]/(rho[i][j] + 0.0000001)
        if(i == 900):
            ux_fine[j] = ux[i][j]
        it = it+1

#--------------------------------------------------------------CALCULANDO A VELOCIDADE EM Y-----------------------------------------------------#
it = 0
for i in range(0,tamX,1):
    for j in range(0,tamY,1):
        uy[i][j] = dataC[it][2] + dataC[it][5] + dataC[it][6] - (dataC[it][4] + dataC[it][7] + dataC[it][8])
        uy[i][j] = uy[i][j]/(rho[i][j]+0.0000001)
        it = it+1

it = 0
for i in range(0,tamX,1):
    for j in range(0,tamY,1):
        uy[i][tamY + j] = dataF[it][2] + dataF[it][5] + dataF[it][6] - (dataF[it][4] + dataF[it][7] + dataF[it][8])
        uy[i][tamY + j] = uy[i][tamY + j]/(rho[i][tamY + j] + 0.0000001)
        it = it+1



###############################################################################
raio = 19.5

dy = 1.0
dx = 1.0
dt = 1.0

tau = 0.56000000

ni = (2*tau - 1)*(dx*dx) / (6*dt)

v_max = 0.003

dP = v_max * 4 * ni * dx / (raio*raio)

vec = np.zeros([(rd*20)])

print('ni: ',ni)
print('dP: ', dP)

for i in range(0,len(vec),1):
    vec[i] = (dP*(raio*raio - np.abs((raio - i)*(raio - i)))) / (4*ni*dx)
###############################################################################

print(rhoMed / (tamX*(rd*tamY)))

fig = plt.figure()
cbRho = plt.imshow(rho)
plt.colorbar(cbRho)
plt.ylabel('y')
plt.xlabel('x')
plt.title('Densidade no domínio')

###############################################################################

fig2 = plt.figure()
cbUx = plt.imshow(ux)
plt.colorbar(cbUx)
plt.ylabel('y')
plt.xlabel('x')
plt.title('Velocidade na direção X no domínio')

###############################################################################

fig3 = plt.figure()
cbUy = plt.imshow(uy)
plt.colorbar(cbUy)
plt.ylabel('y')
plt.xlabel('x')
plt.title('Velocidade na direção Y no domínio')

###############################################################################

fig4 = plt.figure()

namelist = ['Analítico','Numérico']

plt.plot(vec,'-')
#plt.plot(ux_fine,'*',c='red')
plt.plot(ux[900][:],'o')
plt.legend(namelist,loc="upper right")
plt.ylabel('y')
plt.xlabel('x')
plt.title('Perfil de velocidade em X')

plt.show()
###############################################################################
fig5 = plt.figure()

erro = np.zeros(len(vec))
x    = np.arange(len(vec))

#não sei se está certinho assim
for i in range(len(vec)):
    print(i, ' | analitoco: ', vec[i], ' | numérico: ', ux[900][i])
    #erro[i] = (np.absolute(vec[i] - ux[900][i]) / (vec[i] + 0.0001))*100
    erro[i] = np.absolute(vec[i] - ux[900][i])

plt.bar(x, erro, align='center', alpha=0.5)
plt.xticks(x, x)
plt.ylabel('Diferença')
plt.xlabel('Posição')
plt.title('Erro em diferença absoluta')

plt.show()
