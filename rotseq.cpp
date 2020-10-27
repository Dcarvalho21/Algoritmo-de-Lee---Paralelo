//Grupo: Daniel Carvalho de Oliveira
//Compilação: g++ rotpar.cpp -o rotpar -fopenmp -Wall

#include <iostream>
#include <fstream>
#include <deque>
#include <climits>
#include <omp.h>

#define inf INT_MAX
using namespace std;

//Define a estrutura celula, cada celula contem a posição de uma celula do grid
class cel
{
    public:
        int i;
        int j;
};

//Verifica se o vizinho de uma celula existe e retorna a celula vizinha caso sim,
//caso não retorna uma celula com -1 em suas posições
cel verifica_viz(int i, cel cl, int m, int n)
{
    cel vz;
    vz.i = -1;
    vz.j = -1;
    
    if(cl.i < m && cl.j < n)
    {
        switch(i)
        {
            case 0:
                if (cl.i > 0)
                    {
                        vz.i = cl.i -1;
                        vz.j = cl.j;
                    }
            break;

            case 1:
                if(cl.j < n-1)
                    {
                        vz.i = cl.i;
                        vz.j = cl.j + 1;
                    }
            break;
        
            case 2:
                if (cl.i < m-1)
                    {
                        vz.i = cl.i + 1;
                        vz.j = cl.j;
                    }
            break;

            case 3:
                if (cl.j > 0)
                    {
                        vz.i = cl.i;
                        vz.j = cl.j - 1;
                    }
            break;
        }
    }
    
    return vz;
}

int main(int argc, char* argv[])
{
    int m, n, ** grid = NULL; // Declaração da matriz e suas variaveis
    int obstaculos;
    bool achou = false;
    double tini, tfin;
    ifstream entrada;
    ofstream saida;
    cel origem, destino, cl, viz;
    deque <cel> fila, caminho;

    entrada.open(argv[1]);
    
    //Aloca o espaço na memoria da matriz dinamicamente
    entrada >> m;
    entrada >> n;
    grid = new int*[m];
    if (m)
    {
        grid[0] = new int[m * n];
        for (int i = 0; i < m; i++)
            grid[i] = grid[0] + i * n;
    }
    else
    {
        cout << "Insira um valor valido para as linhas ou colunas" << endl;
        return 0;
    }
    
    //Inicializa a matriz com o valor de INT_MAX em suas celulas
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            grid[i][j] = inf;
    
    entrada >> origem.i;
    entrada >> origem.j;
    grid[origem.i][origem.j] = 0;
    
    entrada >> destino.i;
    entrada >> destino.j;
    
    //Le e coloca os obstaculos no grid
    entrada >> obstaculos;
    for (int k = 0; k < obstaculos; k++)
    {
        int i, j;
        cel obst_origem;

        entrada >> obst_origem.i;
        entrada >> obst_origem.j;
        entrada >> i;
        entrada >> j;
        
        for (int l = obst_origem.i; l <= obst_origem.i + i -1; l++)
            for (int c = obst_origem.j; c <= obst_origem.j + j -1; c++)
                grid[l][c] = -1;
    }
    
    entrada.close();

    tini = omp_get_wtime(); //Começa a contar o tempo de execução

    //Algoritmo (fase de expansão)
    fila.emplace_back(origem);

    while (!fila.empty() && !achou)
    {
        cl = fila.front();
        fila.pop_front();
        
        if (cl.i == destino.i && cl.j == destino.j)
            achou = true;
        
        else
        {
            for (int i = 0; i < 4; i++) // Verifica os vizinhos da celula atual em sentido horario
            {
                viz = verifica_viz(i, cl, m, n);
                
                if (viz.i != -1)
                    if (grid[viz.i][viz.j] == inf)
                    {
                        grid[viz.i][viz.j] = grid[cl.i][cl.j] + 1;
                        fila.push_back(viz);
                    }
            }
        }
    }
    
    fila.clear(); // Limpa a fila usada na primeira parte do algoritmo
    fila = deque<cel>(); // Desaloca o espaço usado por essa mesma fila 

    // Algoritmo (fase de backtracking)
    if (achou)
    {
        cl = destino;
        caminho.push_back(cl);
        while(cl.i != origem.i || cl.j != origem.j)
        {
            for (int i = 0; i < 4; i++) // Verifica os vizinhos da celula atual em sentido horario
            {
                viz = verifica_viz(i, cl, m, n);
                
                if (viz.i != -1)
                    if (grid[viz.i][viz.j] == grid[cl.i][cl.j] - 1)
                    {
                        cl = viz;
                        caminho.push_back(viz);
                        i = 4; //Para de checar os vizinhos da celula atual quando achar o certo
                    }
            }
        }
        
        //Termina de contar o tempo de execução e imprime ele na saída padrão
        tfin = omp_get_wtime();
        cout << "Tempo de execução (sequencial): " << tfin - tini << "s" << endl;
        
        //Coloca o resultado do algortimo no arquivo de saida informado
        saida.open(argv[2]);
        saida << grid[destino.i][destino.j] << endl;
        while (!caminho.empty())
        {
            
            saida << caminho.back().i << " " << caminho.back().j << endl;
            caminho.pop_back();
        }
        saida.close();
    }
    else
        cout << "Caminho não encontrado!" << endl;

    // Desaloca o espaço usado pela matriz
    if (m) delete [] grid[0];
    delete [] grid;

    return 0;
}