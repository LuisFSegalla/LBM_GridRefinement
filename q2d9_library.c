#include "q2d9_library.h"

//Função usada para calcular a densidade na célula
//Retorna a densidade da célula no parâmetro rho
void calculaDensidade2D(Celula **funcDistribuicao, int tamX, int tamY)
{
    for(int i = 0; i < tamX; i++)
    {
        for(int j = 0; j < tamY; j++)
        {
            funcDistribuicao[i][j].rho = 0;
            for(int q = 0; q < 9; q++)
            {
                funcDistribuicao[i][j].rho += funcDistribuicao[i][j].f[q];
            }
        }
    }
}

//Função usada para calcular a velocidade de uma célula
//Retorna a velocidade no eixo X e no eixo Y nas variáveis ux e uy respectivamente
void calculaVelocidade2D(Celula **funcDistribuicao, int tamX, int tamY)
{
    for(int i = 0; i < tamX; i++)
    {
        for(int j = 0; j < tamY; j++)
        {
            funcDistribuicao[i][j].ux = 0;
            funcDistribuicao[i][j].uy = 0;

            funcDistribuicao[i][j].ux = funcDistribuicao[i][j].f[1] + funcDistribuicao[i][j].f[5] + funcDistribuicao[i][j].f[8]
                                      - (funcDistribuicao[i][j].f[3] + funcDistribuicao[i][j].f[6] + funcDistribuicao[i][j].f[7]);

            funcDistribuicao[i][j].uy = funcDistribuicao[i][j].f[2] + funcDistribuicao[i][j].f[5] + funcDistribuicao[i][j].f[6]
                                      - (funcDistribuicao[i][j].f[4] + funcDistribuicao[i][j].f[7] + funcDistribuicao[i][j].f[8]);

            funcDistribuicao[i][j].ux = funcDistribuicao[i][j].ux / funcDistribuicao[i][j].rho;
            funcDistribuicao[i][j].uy = funcDistribuicao[i][j].uy / funcDistribuicao[i][j].rho;
        }
    }
}

//Função usada para encontrar a função de equilíbrio baseada na função distribuição e nos pesos para o caso q2d9
//Retorna uma matriz da struct Celula com a função de equilíbrio para todas as direções em cada ponto da matriz dentro da struct Celula funcEquilibrio
void calculaFuncEquilibrio2D(Celula **funcDistribuicao, Celula **funcEquilibrio, double *pesos, int tamX, int tamY)
{
    double rho, ux, uy;
    for(int i = 0; i < tamX; i++)
    {
        for(int j = 0; j < tamY; j++)
        {
            rho = funcDistribuicao[i][j].rho;
            ux  = funcDistribuicao[i][j].ux;
            uy  = funcDistribuicao[i][j].uy;

            funcEquilibrio[i][j].f[0] = rho*pesos[0]*( 1                                         - 1.5*(ux*ux + uy*uy) );

            funcEquilibrio[i][j].f[1] = rho*pesos[1]*( 1 + 3*ux        + 4.5*ux*ux               - 1.5*(ux*ux + uy*uy) );

            funcEquilibrio[i][j].f[2] = rho*pesos[2]*( 1 + 3*uy        + 4.5*uy*uy               - 1.5*(ux*ux + uy*uy) );

            funcEquilibrio[i][j].f[3] = rho*pesos[3]*( 1 - 3*ux        + 4.5*ux*ux               - 1.5*(ux*ux + uy*uy) );

            funcEquilibrio[i][j].f[4] = rho*pesos[4]*( 1 - 3*uy        + 4.5*uy*uy               - 1.5*(ux*ux + uy*uy) );

            funcEquilibrio[i][j].f[5] = rho*pesos[5]*( 1 + 3*(ux + uy) + 4.5*(ux + uy)*(ux + uy) - 1.5*(ux*ux + uy*uy) );

            funcEquilibrio[i][j].f[6] = rho*pesos[6]*( 1 + 3*(uy - ux) + 4.5*(uy - ux)*(uy - ux) - 1.5*(ux*ux + uy*uy) );

            funcEquilibrio[i][j].f[7] = rho*pesos[7]*( 1 - 3*(uy + ux) + 4.5*(ux + uy)*(ux + uy) - 1.5*(ux*ux + uy*uy) );

            funcEquilibrio[i][j].f[8] = rho*pesos[8]*( 1 + 3*(ux - uy) + 4.5*(ux - uy)*(ux - uy) - 1.5*(ux*ux + uy*uy) );
        }
    }
}


//Função usada para calcular a colisão interna das células
//Retorna a nova função distribuição depois da colisão dentro da struct Celula funcDistribuicao
void calculaColisao2D(Celula **funcDistribuicao, Celula **funcEquilibrio, double tau, int tamX, int tamY)
{
    for(int i = 0; i < tamX; i++)
    {
        for(int j = 0; j < tamY; j++)
        {
            funcDistribuicao[i][j].f[0] = funcDistribuicao[i][j].f[0] - ( (funcDistribuicao[i][j].f[0] - funcEquilibrio[i][j].f[0]) / tau );

            funcDistribuicao[i][j].f[1] = funcDistribuicao[i][j].f[1] - ( (funcDistribuicao[i][j].f[1] - funcEquilibrio[i][j].f[1]) / tau );

            funcDistribuicao[i][j].f[2] = funcDistribuicao[i][j].f[2] - ( (funcDistribuicao[i][j].f[2] - funcEquilibrio[i][j].f[2]) / tau );

            funcDistribuicao[i][j].f[3] = funcDistribuicao[i][j].f[3] - ( (funcDistribuicao[i][j].f[3] - funcEquilibrio[i][j].f[3]) / tau );

            funcDistribuicao[i][j].f[4] = funcDistribuicao[i][j].f[4] - ( (funcDistribuicao[i][j].f[4] - funcEquilibrio[i][j].f[4]) / tau );

            funcDistribuicao[i][j].f[5] = funcDistribuicao[i][j].f[5] - ( (funcDistribuicao[i][j].f[5] - funcEquilibrio[i][j].f[5]) / tau );

            funcDistribuicao[i][j].f[6] = funcDistribuicao[i][j].f[6] - ( (funcDistribuicao[i][j].f[6] - funcEquilibrio[i][j].f[6]) / tau );

            funcDistribuicao[i][j].f[7] = funcDistribuicao[i][j].f[7] - ( (funcDistribuicao[i][j].f[7] - funcEquilibrio[i][j].f[7]) / tau );

            funcDistribuicao[i][j].f[8] = funcDistribuicao[i][j].f[8] - ( (funcDistribuicao[i][j].f[8] - funcEquilibrio[i][j].f[8]) / tau );
        }
    }
}


//Função usada para calcular a propagação das particulas dentro do domínio
//Retorna a função pós propagação dentro da matriz de struct Celula funcDistribuicao
void calculaPropagacao2D(Celula **funcDistribuicao, int tamX, int tamY)
{
    for(int i = 0; i < tamX-1; i++)
    {
        for(int j = 0; j < tamY-1; j++)
        {
            funcDistribuicao[i][j].f[7] = funcDistribuicao[i+1][j+1].f[7];
            funcDistribuicao[i][j].f[3] = funcDistribuicao[i+1][j].f[3];
        }
    }

    for(int i = 0; i < tamX-1; i++)
    {
        for(int j = tamY-1; j > 0; j--)
        {
            funcDistribuicao[i][j].f[6] = funcDistribuicao[i+1][j-1].f[6];
            funcDistribuicao[i][j].f[2] = funcDistribuicao[i][j-1].f[2];
        }
    }

    for(int i = tamX-1; i > 0; i--)
    {
        for(int j = tamY-1; j > 0; j--)
        {
            funcDistribuicao[i][j].f[5] = funcDistribuicao[i-1][j-1].f[5];
            funcDistribuicao[i][j].f[1] = funcDistribuicao[i-1][j].f[1];
        }
    }

    for(int i = tamX-1; i > 0; i--)
    {
        for(int j = 0; j < tamY-1; j++)
        {
            funcDistribuicao[i][j].f[8] = funcDistribuicao[i-1][j+1].f[8];
            funcDistribuicao[i][j].f[4] = funcDistribuicao[i][j+1].f[4];
        }
    }



    for(int j = 0; j < tamY-1; j++)
    {
        funcDistribuicao[0][j].f[4] = funcDistribuicao[0][j+1].f[4];
    }

    for(int j = tamY-1; j > 0; j--)
    {
        funcDistribuicao[tamX-1][j].f[2] = funcDistribuicao[tamX-1][j-1].f[2];
    }

    //--------------------------------------------------------------------------ADICIONEI MAS NÃO TENHO CERTEZA
    for(int i = 0; i < tamX-1; i++)
    {
      funcDistribuicao[i][tamY-1].f[3] = funcDistribuicao[i+1][tamY-1].f[3];
    }
    for(int i = tamX-1; i > 0; i--)
    {
      funcDistribuicao[i][0].f[1] = funcDistribuicao[i-1][0].f[1];
    }

}

//------------------------------------------------------------------------------funções em relação a malha refinada-----------------------------------------------------//

void calculaExplosao2D(Celula **funcDistribuicao_c, Celula **funcDistribuicao_f, int tamX_c, int tamY_c, int tamX_f, int tamY_f, int razaoDeDivisao)//funcionando
{//não está mt bom, mas funciona
  double divValue = 0;
  for(int i = 0; i < tamX_c; i++)
  {
    for(int c = 0; c < 9; c++)
    {
      divValue = funcDistribuicao_c[i][0].f[c]/(razaoDeDivisao*razaoDeDivisao);
      for(int j = 0; j < razaoDeDivisao; j++)
      {
        for(int i2 = 0; i2 < razaoDeDivisao; i2++)
        {
          funcDistribuicao_f[i*razaoDeDivisao + i2][tamY_f + j].f[c] = divValue;
        }
      }
    }
  }

}

void calculaDensidade2D_fronteira(Celula **funcDistribuicao, int tamX_f, int tamY_f, int razaoDeDivisao)
{
    for(int i = 0; i < tamX_f; i++)
    {
        for(int j = 0; j < tamY_f+razaoDeDivisao; j++)
        {
          funcDistribuicao[i][j].rho = 0;
          for(int q = 0; q < 9; q++)
          {
            funcDistribuicao[i][j].rho += funcDistribuicao[i][j].f[q];
          }
        }
    }
}
