/*
Implementação de um método de refino de malha incorporada usando LBM. Baseado no artigo "A generic, mass conservative local grid refinement technique for lattice-Boltzmann schemes"
compilado usando a versão gcc (Ubuntu 7.4.0-1ubuntu1~18.04.1) 7.4.0
comando para compilação:
  gcc -o main2D main2D.c q2d9_library.c -lm

  ou

  make

O problema consiste basicamente em resolver duas matrizes separadas, cada uma referente a uma das malhas e lidar com os pontos de intersecção entre as mesmas.
*/

#include "q2d9_library.h"

void EscreveDirecoes(Celula **funcDistribuicao, char *nome, int tamX, int tamY)
{
  FILE *arq;
  arq = fopen(nome,"w");
  for(int i = 0; i < tamX; i++)
  {
    for(int j = 0; j < tamY; j++)
    {
      for(int c = 0; c < 9; c++)
      {
        fprintf(arq, " %f ", funcDistribuicao[i][j].f[c]);
      }
      fprintf(arq, "\n");
    }
  }

}

/*
Função usada para escrever as direções de cada célula da matriz refinada na ordem que elas coincidem com a matriz grosseira para o pós processamento
*/

void EscreveDirecoes_refinada(Celula **funcDistribuicao, char *nome,int tamX, int tamY, int razaoDeDivisao)
{
  FILE *arq;
  arq = fopen(nome,"w");
  for(int i = 0; i < tamX; i++)
  {
    for(int j = 0; j < tamY; j++)
    {
      for(int i2 = 0; i2 < razaoDeDivisao; i2++)
      {
        for(int j2 = 0; j2 < razaoDeDivisao; j2++)
        {
          for(int c = 0; c < 9; c++)
          {
            fprintf(arq, " %f ", funcDistribuicao[i*razaoDeDivisao + i2][j*razaoDeDivisao + j2].f[c]);
          }
          fprintf(arq, "\n");
        }
      }
    }
  }
}

void EscreveLinha(Celula **funcDistribuicao, char *nome, int tamX, int linha)
{
  FILE *arq;
  arq = fopen(nome,"w");
  //printf("nome: %s\nlinha: %d\n", nome, linha);
  for(int i = 0; i < tamX; i++)
  {
    fprintf(arq, " |");
    for(int c = 0; c < 9; c++)
    {
      fprintf(arq, " %f ", funcDistribuicao[i][linha].f[c]);
    }
    fprintf(arq, "|");
  }
}

void EscreveTudo(Celula **funcDistribuicao, char *nome, int tamX, int tamY)
{
    FILE *arq;
    arq = fopen(nome,"w");
    fprintf(arq, "#ux     uy      rho\n");
    for(int i = 0; i < tamX; i++)
    {
        for(int j = 0; j < tamY; j ++)
        {
            fprintf(arq,"%f    %f    %f\n", funcDistribuicao[i][j].ux, funcDistribuicao[i][j].uy, funcDistribuicao[i][j].rho);
        }
    }
    fclose(arq);
}

void escreveMatriz(Celula **funcDistribuicao, int tamX, int tamY)
{
    FILE *arqX, *arqY, *arqRho;
    arqX   = fopen("Matriz da Velocidade em X.txt","w");
    arqY   = fopen("Matriz da Velocidade em Y.txt","w");
    arqRho = fopen("Matriz da Densidade.txt","w");

    for(int i = 0; i < tamX; i++)
    {
        for(int j = 0; j < tamY; j++)
        {
            fprintf(arqX,   " %f ", funcDistribuicao[i][j].ux);
            fprintf(arqY,   " %f ", funcDistribuicao[i][j].uy);
            fprintf(arqRho, " %f ", funcDistribuicao[i][j].rho);
        }
        fprintf(arqX,   "\n");
        fprintf(arqY,   "\n");
        fprintf(arqRho, "\n");
    }
    fclose(arqX);
    fclose(arqY);
    fclose(arqRho);
}

void escreveMatriz_refinada(Celula **funcDistribuicao, int tamX, int tamY)
{
    FILE *arqX, *arqY, *arqRho;
    arqX   = fopen("Matriz da Velocidade em X - refinada.txt","w");
    arqY   = fopen("Matriz da Velocidade em Y - refinada.txt","w");
    arqRho = fopen("Matriz da Densidade - refinada.txt","w");

    for(int i = 0; i < tamX; i++)
    {
        for(int j = 0; j < tamY; j++)
        {
            fprintf(arqX,   " %f ", funcDistribuicao[i][j].ux);
            fprintf(arqY,   " %f ", funcDistribuicao[i][j].uy);
            fprintf(arqRho, " %f ", funcDistribuicao[i][j].rho);
        }
        fprintf(arqX,   "\n");
        fprintf(arqY,   "\n");
        fprintf(arqRho, "\n");
    }
    fclose(arqX);
    fclose(arqY);
    fclose(arqRho);
}

void condContornoPoiseuille_coarse(Celula **funcDistribuicao, int tamX, int tamY, double velX, double velY)
{
    double rhoAux = 0;
    for(int j = 0; j < tamY; j++)
    {
        rhoAux = 0;
        //known velocity na face da entrada
        rhoAux = funcDistribuicao[0][j].f[0] + funcDistribuicao[0][j].f[2] + funcDistribuicao[0][j].f[4];
        rhoAux += 2*(funcDistribuicao[0][j].f[3] + funcDistribuicao[0][j].f[6] + funcDistribuicao[0][j].f[7]);
        rhoAux = (double)(rhoAux / (1.0 - velX));

        funcDistribuicao[0][j].f[1] = funcDistribuicao[0][j].f[3] + (2.0*rhoAux*velX/3.0);
        funcDistribuicao[0][j].f[5] = funcDistribuicao[0][j].f[7] + (rhoAux*velX/6.0);
        funcDistribuicao[0][j].f[8] = funcDistribuicao[0][j].f[6] + (rhoAux*velX/6.0);

        //--------------------------------------------------------------------------------------------------------//

        //open boundary na face da saída
        if(j>0)//cuidar
        {
          funcDistribuicao[tamX-1][j].f[1] = 2.0*funcDistribuicao[tamX-2][j].f[1] - funcDistribuicao[tamX-3][j].f[1];
          funcDistribuicao[tamX-1][j].f[5] = 2.0*funcDistribuicao[tamX-2][j].f[5] - funcDistribuicao[tamX-3][j].f[5];
          funcDistribuicao[tamX-1][j].f[8] = 2.0*funcDistribuicao[tamX-2][j].f[8] - funcDistribuicao[tamX-3][j].f[8];
        }
    }

    for(int i = 0; i < tamX; i++)
    {
      //bounce back na face do topo
      funcDistribuicao[i][tamY-1].f[4] = funcDistribuicao[i][tamY-1].f[2];
      funcDistribuicao[i][tamY-1].f[7] = funcDistribuicao[i][tamY-1].f[5];
      funcDistribuicao[i][tamY-1].f[8] = funcDistribuicao[i][tamY-1].f[6];
    }
}

void condContornoPoiseuille_fine(Celula **funcDistribuicao, int tamX, int tamY, double velX, double velY)
{
    double rhoAux = 0;
    for(int j = 0; j < tamY; j++)
    {
        rhoAux = 0;
        //known velocity na face da entrada
        rhoAux = funcDistribuicao[0][j].f[0] + funcDistribuicao[0][j].f[2] + funcDistribuicao[0][j].f[4];
        rhoAux += 2*(funcDistribuicao[0][j].f[3] + funcDistribuicao[0][j].f[6] + funcDistribuicao[0][j].f[7]);
        rhoAux = (double)(rhoAux / (1.0 - velX));

        funcDistribuicao[0][j].f[1] = funcDistribuicao[0][j].f[3] + (2.0*rhoAux*velX/3.0);
        funcDistribuicao[0][j].f[5] = funcDistribuicao[0][j].f[7] + (rhoAux*velX/6.0);
        funcDistribuicao[0][j].f[8] = funcDistribuicao[0][j].f[6] + (rhoAux*velX/6.0);

        //--------------------------------------------------------------------------------------------------------//

        //open boundary na face da saída
        if(j>0)
        {
          funcDistribuicao[tamX-1][j].f[1] = 2.0*funcDistribuicao[tamX-2][j].f[1] - funcDistribuicao[tamX-3][j].f[1];
          funcDistribuicao[tamX-1][j].f[5] = 2.0*funcDistribuicao[tamX-2][j].f[5] - funcDistribuicao[tamX-3][j].f[5];
          funcDistribuicao[tamX-1][j].f[8] = 2.0*funcDistribuicao[tamX-2][j].f[8] - funcDistribuicao[tamX-3][j].f[8];
        }
    }

    for(int i = 0; i < tamX; i++)
    {
      //bounce back no fundo
      funcDistribuicao[i][0].f[2] = funcDistribuicao[i][0].f[4];
      funcDistribuicao[i][0].f[5] = funcDistribuicao[i][0].f[7];
      funcDistribuicao[i][0].f[6] = funcDistribuicao[i][0].f[8];
    }
}


void condContornoBordas(Celula **funcDistribuicao_c, Celula **funcEquilibrio_c, int tamX_c, int tamY_c, Celula **funcDistribuicao_f, Celula **funcEquilibrio_f, int tamX_f, int tamY_f)
{
    double rho = 0.0;
    double ux  = 0.0;
    //superior esquerdo
    rho = funcDistribuicao_c[0][tamY_c-1].rho;
    ux  = funcDistribuicao_c[0][tamY_c-1].ux;
    funcDistribuicao_c[0][tamY_c-1].f[1] = funcDistribuicao_c[0][tamY_c-1].f[3] + (funcEquilibrio_c[0][tamY_c-1].f[1] - funcEquilibrio_c[0][tamY_c-1].f[3]);
    funcDistribuicao_c[0][tamY_c-1].f[4] = funcDistribuicao_c[0][tamY_c-1].f[2] + (funcEquilibrio_c[0][tamY_c-1].f[4] - funcEquilibrio_c[0][tamY_c-1].f[2]);
    funcDistribuicao_c[0][tamY_c-1].f[8] = funcDistribuicao_c[0][tamY_c-1].f[6] + (funcEquilibrio_c[0][tamY_c-1].f[8] - funcEquilibrio_c[0][tamY_c-1].f[6]);
    funcDistribuicao_c[0][tamY_c-1].f[7] = 0.5*(rho - rho*ux - 2*funcDistribuicao_c[0][tamY_c-1].f[6] - 2*funcDistribuicao_c[0][tamY_c-1].f[3] - funcDistribuicao_c[0][tamY_c-1].f[0] - funcDistribuicao_c[0][tamY_c-1].f[2] - funcDistribuicao_c[0][tamY_c-1].f[4]);
    funcDistribuicao_c[0][tamY_c-1].f[5] = rho*ux + funcDistribuicao_c[0][tamY_c-1].f[7] + funcDistribuicao_c[0][tamY_c-1].f[6] + funcDistribuicao_c[0][tamY_c-1].f[3] - funcDistribuicao_c[0][tamY_c-1].f[8] - funcDistribuicao_c[0][tamY_c-1].f[1];

    //superior direito
    rho = funcDistribuicao_c[tamX_c-1][tamY_c-1].rho;
    ux  = funcDistribuicao_c[tamX_c-1][tamY_c-1].ux;
    funcDistribuicao_c[tamX_c-1][tamY_c-1].f[3] = funcDistribuicao_c[tamX_c-1][tamY_c-1].f[1] + (funcEquilibrio_c[tamX_c-1][tamY_c-1].f[3] -funcEquilibrio_c[tamX_c-1][tamY_c-1].f[1]);
    funcDistribuicao_c[tamX_c-1][tamY_c-1].f[4] = funcDistribuicao_c[tamX_c-1][tamY_c-1].f[2] + (funcEquilibrio_c[tamX_c-1][tamY_c-1].f[4] -funcEquilibrio_c[tamX_c-1][tamY_c-1].f[2]);
    funcDistribuicao_c[tamX_c-1][tamY_c-1].f[7] = funcDistribuicao_c[tamX_c-1][tamY_c-1].f[5] + (funcEquilibrio_c[tamX_c-1][tamY_c-1].f[7] -funcEquilibrio_c[tamX_c-1][tamY_c-1].f[5]);
    funcDistribuicao_c[tamX_c-1][tamY_c-1].f[6] = 0.5*(rho - rho*ux - funcDistribuicao_c[tamX_c-1][tamY_c-1].f[0] - funcDistribuicao_c[tamX_c-1][tamY_c-1].f[2] - funcDistribuicao_c[tamX_c-1][tamY_c-1].f[4] - 2*funcDistribuicao_c[tamX_c-1][tamY_c-1].f[3] - 2*funcDistribuicao_c[tamX_c-1][tamY_c-1].f[7]);
    funcDistribuicao_c[tamX_c-1][tamY_c-1].f[8] = rho*ux + funcDistribuicao_c[tamX_c-1][tamY_c-1].f[6] + funcDistribuicao_c[tamX_c-1][tamY_c-1].f[7] + funcDistribuicao_c[tamX_c-1][tamY_c-1].f[3] - funcDistribuicao_c[tamX_c-1][tamY_c-1].f[1] - funcDistribuicao_c[tamX_c-1][tamY_c-1].f[5];

    //infeior direito
    rho = funcDistribuicao_f[tamX_f-1][0].rho;
    ux  = funcDistribuicao_f[tamX_f-1][0].ux;
    funcDistribuicao_f[tamX_f-1][0].f[3] = funcDistribuicao_f[tamX_f-1][0].f[1] + (funcEquilibrio_f[tamX_f-1][0].f[3] - funcEquilibrio_f[tamX_f-1][0].f[1]);
    funcDistribuicao_f[tamX_f-1][0].f[2] = funcDistribuicao_f[tamX_f-1][0].f[4] + (funcEquilibrio_f[tamX_f-1][0].f[2] - funcEquilibrio_f[tamX_f-1][0].f[4]);
    funcDistribuicao_f[tamX_f-1][0].f[6] = funcDistribuicao_f[tamX_f-1][0].f[8] + (funcEquilibrio_f[tamX_f-1][0].f[6] - funcEquilibrio_f[tamX_f-1][0].f[8]);
    funcDistribuicao_f[tamX_f-1][0].f[7] = 0.5*(rho - rho*ux - 2*funcDistribuicao_f[tamX_f-1][0].f[6] - 2*funcDistribuicao_f[tamX_f-1][0].f[3] - funcDistribuicao_f[tamX_f-1][0].f[0] - funcDistribuicao_f[tamX_f-1][0].f[2] - funcDistribuicao_f[tamX_f-1][0].f[4]);
    funcDistribuicao_f[tamX_f-1][0].f[5] = rho*ux + funcDistribuicao_f[tamX_f-1][0].f[7] + funcDistribuicao_f[tamX_f-1][0].f[6] + funcDistribuicao_f[tamX_f-1][0].f[3] - funcDistribuicao_f[tamX_f-1][0].f[8] - funcDistribuicao_f[tamX_f-1][0].f[1];

    //inferior esquerdo
    rho = funcDistribuicao_f[0][0].rho;
    ux  = funcDistribuicao_f[0][0].ux;
    funcDistribuicao_f[0][0].f[2] = funcDistribuicao_f[0][0].f[4] + (funcDistribuicao_f[0][0].f[2] - funcDistribuicao_f[0][0].f[4]);
    funcDistribuicao_f[0][0].f[1] = funcDistribuicao_f[0][0].f[3] + (funcDistribuicao_f[0][0].f[1] - funcDistribuicao_f[0][0].f[3]);
    funcDistribuicao_f[0][0].f[5] = funcDistribuicao_f[0][0].f[7] + (funcDistribuicao_f[0][0].f[5] - funcDistribuicao_f[0][0].f[7]);
    funcDistribuicao_f[0][0].f[6] = 0.5*(rho - rho*ux - funcDistribuicao_f[0][0].f[0] - funcDistribuicao_f[0][0].f[2] - funcDistribuicao_f[0][0].f[4] - 2*funcDistribuicao_f[0][0].f[3] - 2*funcDistribuicao_f[0][0].f[7]);
    funcDistribuicao_f[0][0].f[8] = rho*ux + funcDistribuicao_f[0][0].f[6] + funcDistribuicao_f[0][0].f[7] + funcDistribuicao_f[0][0].f[3] - funcDistribuicao_f[0][0].f[1] - funcDistribuicao_f[0][0].f[5];
}

int main()
{
    //pesos
    double centro = 4.0/9.0;
    double retas = 1.0/9.0;
    double diagonais = 1.0/36.0;
    double pesos[9]  = {centro, retas, retas, retas, retas, diagonais, diagonais, diagonais, diagonais};

    //razão de divisão
    int razaoDeDivisao = 2;

    int tamX_c = 1000, tamY_c = 20;//malha grosseira
    int tamX_f = tamX_c*razaoDeDivisao, tamY_f = tamY_c*razaoDeDivisao;//malha refinada

    //tamanho do domínio
    double dx = 1.0;
    double dy = 1.0;


    //características dos passos
    double dt    = 1.0;
    double alpha = 0.01;
    double a     = 3*alpha;
    double omega = (double)(1.0 / ( a + 0.5 ));//da onde vem esse termo exatamente

    //tempo de relaxação
    double tau_c =  1/omega;
    double tau_f = tau_c*(razaoDeDivisao + 0.5) + 0.5;//no artigo falta um parenteses no denominador
    printf("tau_c: %f\ntau_f: %f\n",tau_c,tau_f);

    int numSteps = 1000;

    //caracteristicas do fluido
    double ux0   = 0.002;
    double uy0   = 0.0;
    double rho   = 1.0;
    double Re    = ux0*tamX_c / alpha;

    printf("Reynolds: %f\n", Re);

    //alocação das matrizes que representam o domínio
    Celula **funcDistribuicao_c, **funcEquilibrio_c;
    Celula **funcDistribuicao_f, **funcEquilibrio_f;

    funcDistribuicao_c = (Celula **)malloc(tamX_c*sizeof(Celula *));
    funcEquilibrio_c   = (Celula **)malloc(tamX_c*sizeof(Celula *));

    funcDistribuicao_f = (Celula **)malloc(tamX_f*sizeof(Celula *));
    funcEquilibrio_f   = (Celula **)malloc(tamX_f*sizeof(Celula *));

    for(int i = 0; i < tamX_c; i++)
    {
        funcDistribuicao_c[i] = (Celula *)malloc(tamY_c*sizeof(Celula));
        funcEquilibrio_c[i]   = (Celula *)malloc(tamY_c*sizeof(Celula));
    }

    //Aqui tenho que alocar a quantidade de colunas refinadas em Y mais a razão de divisão que representa as células da interface
    for(int i = 0; i < tamX_f; i++)
    {
      funcDistribuicao_f[i] = (Celula *)malloc((tamY_f+razaoDeDivisao)*sizeof(Celula));
      funcEquilibrio_f[i]   = (Celula *)malloc((tamY_f+razaoDeDivisao)*sizeof(Celula));
    }

    //inicializo o problema para a malha grosseira
    for(int i = 0; i < tamX_c; i++)
    {
        for(int j = 0; j < tamY_c; j ++)
        {
            funcDistribuicao_c[i][j].ux  = ux0;
            funcDistribuicao_c[i][j].uy  = uy0;
            funcDistribuicao_c[i][j].rho = rho;
        }
    }

    //inicializo o problema para a malha refinada
    for(int i = 0; i < tamX_f; i++)
    {
        for(int j = 0; j < tamY_f; j ++)
        {
            funcDistribuicao_f[i][j].ux  = ux0;
            funcDistribuicao_f[i][j].uy  = uy0;
            funcDistribuicao_f[i][j].rho = rho/(razaoDeDivisao*razaoDeDivisao);
        }
    }

    //função distribuição original calculada como a função de equilíbrio no estado inicial
    calculaFuncEquilibrio2D(funcDistribuicao_c,funcDistribuicao_c,pesos,tamX_c,tamY_c);
    calculaFuncEquilibrio2D(funcDistribuicao_f,funcDistribuicao_f,pesos,tamX_f,tamY_f);

    //laço principal
    for(int k = 0; k < numSteps; k++)
    {
        printf("laço: %d\n",k );

        //calcula a nova função de equilíbrio
        calculaFuncEquilibrio2D(funcDistribuicao_c,funcEquilibrio_c,pesos,tamX_c,tamY_c);
        calculaFuncEquilibrio2D(funcDistribuicao_f,funcEquilibrio_f,pesos,tamX_f,tamY_f);

        //calcula a colisão | passo 1
        calculaColisao2D(funcDistribuicao_c,funcEquilibrio_c,tau_c,tamX_c,tamY_c);
        calculaColisao2D(funcDistribuicao_f,funcEquilibrio_f,tau_f,tamX_f,tamY_f);

        //calcula a explosão | passo 2
        calculaExplosao2D(funcDistribuicao_c,funcDistribuicao_f,tamX_c,tamY_c,tamX_f,tamY_f,razaoDeDivisao);

        //propagação | passo 3
        calculaPropagacao2D(funcDistribuicao_c,tamX_c,tamY_c);
        calculaPropagacao2D(funcDistribuicao_f,tamX_f,tamY_f+razaoDeDivisao);//propaga na matriz refinada e na interface

        //passo 4
        for(int i = 0; i < (razaoDeDivisao)-1; i++)
        {
          calculaFuncEquilibrio2D(funcDistribuicao_f,funcEquilibrio_f,pesos,tamX_f,tamY_f);
          calculaColisao2D(funcDistribuicao_f,funcEquilibrio_f,tau_f,tamX_f,tamY_f);//colide duas vezes somente na matriz refinada
          //calculaColisao2D(funcDistribuicao_f,funcEquilibrio_f,tau_f,tamX_f,tamY_f+razaoDeDivisao-1);//colide na intersecção (para ajustar os tempos) | parece estar dando problemas
          calculaPropagacao2D(funcDistribuicao_f,tamX_f,tamY_f+razaoDeDivisao);//propaga na matriz refinada e na interface
        }

        //passo 5 | tá bem ruim isso daqui
        for(int i = 0; i < tamX_c; i++)
        {
          funcDistribuicao_c[i][0].f[2] = 0;
          funcDistribuicao_c[i][0].f[5] = 0;
          funcDistribuicao_c[i][0].f[6] = 0;

          for(int i2 = 0; i2 < razaoDeDivisao; i2++)
          {
            for(int j2 = 0; j2 < razaoDeDivisao; j2++)
            {
              funcDistribuicao_c[i][0].f[2] += funcDistribuicao_f[i*razaoDeDivisao + i2][tamY_f + j2].f[2];
              funcDistribuicao_c[i][0].f[5] += funcDistribuicao_f[i*razaoDeDivisao + i2][tamY_f + j2].f[5];
              funcDistribuicao_c[i][0].f[6] += funcDistribuicao_f[i*razaoDeDivisao + i2][tamY_f + j2].f[6];
            }
          }
        }

        //Condições de contorno
        condContornoPoiseuille_coarse(funcDistribuicao_c,tamX_c,tamY_c,ux0,uy0);
        condContornoPoiseuille_fine(funcDistribuicao_f,tamX_f,tamY_f,ux0,uy0);
        condContornoBordas(funcDistribuicao_c, funcEquilibrio_c, tamX_c, tamY_c, funcDistribuicao_f, funcEquilibrio_f, tamX_f, tamY_f);

        //Recalcula a densidade e a velocidade
        calculaDensidade2D(funcDistribuicao_c,tamX_c,tamY_c);
        calculaDensidade2D(funcDistribuicao_f,tamX_f,tamY_f);

        calculaVelocidade2D(funcDistribuicao_c,tamX_c,tamY_c);
        calculaVelocidade2D(funcDistribuicao_f,tamX_f,tamY_f);

        for(int j = 1; j < tamY_c; j++)
        {
          funcDistribuicao_c[0][j].uy = 0.0;
          funcDistribuicao_c[tamX_c-1][j].uy = 0.0;
        }

        for(int j = 1; j < tamY_f; j++)
        {
          funcDistribuicao_f[0][j].uy = 0.0;
          funcDistribuicao_f[tamX_f-1][j].uy = 0.0;
        }
    }

    EscreveDirecoes(funcDistribuicao_c,"direções grosseira.txt",tamX_c,tamY_c);
    //EscreveDirecoes_refinada(funcDistribuicao_f,"direções refinadas.txt",tamX_c,tamY_c,razaoDeDivisao);
    escreveMatriz(funcDistribuicao_c,tamX_c,tamY_c);
    //escreveMatriz_refinada(funcDistribuicao_f,tamX_f,tamY_f);

    for(int i = 0; i < tamX_c; i++)
    {
      for(int j = 0; j < tamY_c; j++)
      {
        for(int c = 0; c < 9; c++)
        {
          funcDistribuicao_c[i][j].f[c] = funcDistribuicao_f[i*razaoDeDivisao][j*razaoDeDivisao].f[c] + funcDistribuicao_f[i*razaoDeDivisao+1][j*razaoDeDivisao].f[c] + funcDistribuicao_f[i*razaoDeDivisao][j*razaoDeDivisao+1].f[c] + funcDistribuicao_f[i*razaoDeDivisao+1][j*razaoDeDivisao+1].f[c];
        }
      }
    }
    calculaDensidade2D(funcDistribuicao_c,tamX_c,tamY_c);
    calculaVelocidade2D(funcDistribuicao_c,tamX_c,tamY_c);
    escreveMatriz_refinada(funcDistribuicao_c,tamX_c,tamY_c);
    EscreveDirecoes(funcDistribuicao_c,"direções Pós-Processamento.txt",tamX_c,tamY_c);

    return 0;
}
