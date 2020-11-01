// Grupo: Daniel Carvalho de Oliveira
// SO utilizado: MX Linux 19.2 (debian-based)
// Compilador utilizado: g++ 8.3.0
// Compilação: g++ rotpar.cpp -o rotpar -fopenmp -Wall

#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <climits>
#include <omp.h>

#define inf INT_MAX
using namespace std;

// Define a estrutura celula, cada celula contem a posição de uma celula do grid
class cel
{
    public:
        int i;
        int j;
};

// Verifica se o vizinho de uma celula existe e retorna a celula vizinha caso sim,
// caso não retorna uma celula com -1 em suas posições
cel verifica_viz(int i, cel cl, const int m, const int n)
{
    cel vz;
    vz.i = -1;
    vz.j = -1;
    
    if (cl.i < m && cl.j < n)
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

// Declaração de redução do openmp pelo usuário, nesse caso uma união de vetores 
#pragma omp declare reduction (merge : std::vector<cel> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

int main(int argc, char* argv[])
{
    int m, n, ** grid = NULL; // Declaração da matriz e suas variaveis
    int obstaculos;
    bool achou = false;
    double tini, tfin;
    ifstream entrada;
    ofstream saida;
    cel origem, destino;
    deque <cel> caminho;
    vector <cel> fila;

    entrada.open(argv[1]);
    
    // Aloca o espaço na memoria da matriz
    entrada >> m;
    entrada >> n;
    grid = new int*[m];
    for (int i = 0; i < m; ++i)
    {
        grid[i] = new int[n];
    }
    
    // Inicializa a matriz com o valor de INT_MAX em suas celulas
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            grid[i][j] = inf;
    
    entrada >> origem.i;
    entrada >> origem.j;
    grid[origem.i][origem.j] = 0;
    
    entrada >> destino.i;
    entrada >> destino.j;
    
    // Lê e coloca os obstáculos no grid
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

    tini = omp_get_wtime(); // Começa a contar o tempo de execução

    // Algoritmo (fase de expansão)
    fila.push_back(origem);

    while (!fila.empty() && !achou)
    {     
        long unsigned int f = fila.size();
        vector<cel> fila_aux;
        
        #pragma omp parallel for reduction(merge:fila_aux)
        for (long unsigned int c = 0; c < f; c++)
        {
            cel cl; 
            cl = fila.at(c);
            
            if (cl.i == destino.i && cl.j == destino.j)
                achou = true;
                
            else
            {
                for (int i = 0; i < 4; i++) // Verifica os vizinhos da celula
                {
                    cel viz;
                    viz = verifica_viz(i, cl, m, n);

                    if (viz.i > -1 && viz.j > -1)
                        if (grid[viz.i][viz.j] == inf)
                        {
                            grid[viz.i][viz.j] = grid[cl.i][cl.j] + 1; // Atualiza o valor da distância da célula no grid
                            fila_aux.push_back(viz); // Insere o vizinho na lista de celulas a serem visitadas
                        }
                }
            }
        }
        fila = fila_aux; // Passa para o próximo nível de expansão
    }
    
    if (!fila.empty())
        fila.clear(); // Limpa a fila usada na primeira parte do algoritmo
    fila = vector<cel>(); // Desaloca o espaço usado por essa mesma fila 
    

    // Algoritmo (fase de backtracking)
    if (achou)
    {
        cel cl;
        cl = destino;
        caminho.push_back(cl);
        
        while(cl.i != origem.i || cl.j != origem.j)
        {
            for (int i = 0; i < 4; i++) // Verifica os vizinhos da celula
            {
                cel viz;
                viz = verifica_viz(i, cl, m, n);
                
                if (viz.i != -1)
                    if (grid[viz.i][viz.j] == grid[cl.i][cl.j] - 1)
                    {
                        cl = viz;
                        caminho.push_back(viz); // Monta o caminho da celula destino à celula origem
                        i = 4; // Para de checar os vizinhos da celula atual quando achar o certo
                    }
            }
        }
        
        // Termina de contar o tempo de execução e imprime ele na saída padrão
        tfin = omp_get_wtime();
        cout << "Tempo de execução (paralelo): " << tfin - tini << "s" << endl;
        
        // Coloca o resultado do algortimo no arquivo de saída informado
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
    if (grid != NULL) 
        delete [] grid[0];
    delete [] grid;
    grid = NULL;
    
    return 0;
}
