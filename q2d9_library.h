#ifndef Q2D9_LIBRARY_H
#define Q2D9_LIBRARY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

/*
 * distribuição das direcoes | baseada no livro pg. 70
 * |6 2 5|
 * |3 0 1|
 * |7 4 8|
 */

/*
  Estrutura que representa as células no domínio para a configuração Q2D9.
  Cada campo representa uma das caracterísiticas que estamos interessados.
*/
struct celula_t
{
    double rho;
    double ux;
    double uy;

    double f[9];//função de distribuição
}typedef Celula;

//Função usada para calcular a densidade na célula
//Retorna a densidade da célula no parâmetro rho
void calculaDensidade2D(Celula **funcDistribuicao, int tamX, int tamY);

//Função usada para calcular a velocidade de uma célula
//Retorna a velocidade no eixo X e no eixo Y nas variáveis ux e uy respectivamente
void calculaVelocidade2D(Celula **funcDistribuicao, int tamX, int tamY);

//Função usada para encontrar a função de equilíbrio baseada na função distribuição e nos pesos para o caso q2d9
//Retorna uma matriz da struct Celula com a função de equilíbrio para todas as direções em cada ponto da matriz dentro da struct Celula funcEquilibrio
void calculaFuncEquilibrio2D(Celula **funcDistribuicao, Celula **funcEquilibrio, double *pesos, int tamX, int tamY);

//Função usada para calcular a colisão interna das células
//Retorna a nova função distribuição depois da colisão dentro da struct Celula funcDistribuicao
void calculaColisao2D(Celula **funcDistribuicao, Celula **funcEquilibrio, double tau, int tamX, int tamY);

//Função usada para calcular a propagação das particulas dentro do domínio
//Retorna a função pós propagação dentro da matriz de struct Celula funcDistribuicao
void calculaPropagacao2D(Celula **funcDistribuicao, int tamX, int tamY);

//Calcula o processo de explosão, neste caso "explode" a última linha da matriz grosseira para a interface (primeiras linhas) da matriz refinada.
void calculaExplosao2D(Celula **funcDistribuicao_c, Celula **funcDistribuicao_f, int tamX_c, int tamY_c, int tamX_f, int tamY_f, int razaoDeDivisao);

//Calcula a densidade para todo o domínio refinado, incluindo a região de interface. Usado principalmente para depurar o código
void calculaDensidade2D_fronteira(Celula **funcDistribuicao, int tamX_f, int tamY_f, int razaoDeDivisao);

/*****************************************************************************************************************************************************/

#endif
