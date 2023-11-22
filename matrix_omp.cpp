#include <bits/stdc++.h>
#include <cfloat>
#include <omp.h>
using namespace std;

const int N = 15;
int final_path[N + 1];
bool visited[N];
//int final_res = INT_MAX;
double final_res = DBL_MAX;


// Definición de la función para calcular la distancia entre dos puntos
double calcularDistancia(std::pair<double, double> punto1, std::pair<double, double> punto2) {
    double dx = punto2.first - punto1.first;
    double dy = punto2.second - punto1.second;
    return std::sqrt(dx * dx + dy * dy);
}

// Función para crear la matriz de distancias
std::vector<std::vector<double>> crearMatrizDistancias(const std::vector<std::pair<double, double>>& puntos) {
    int n = puntos.size();

    // Inicializar la matriz con distancias
    std::vector<std::vector<double>> matrizDistancias(n, std::vector<double>(n, 0.0));

    // Calcular las distancias y almacenar en la matriz
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double distancia = calcularDistancia(puntos[i], puntos[j]);
            matrizDistancias[i][j] = distancia;
            matrizDistancias[j][i] = distancia; // La matriz es simétrica
        }
    }

    return matrizDistancias;
}

// Función para imprimir la matriz de distancias
void imprimirMatrizDistancias(const std::vector<std::vector<double>>& matrizDistancias) {
    for (const auto& fila : matrizDistancias) {
        for (double distancia : fila) {
            std::cout << distancia << " ";
        }
        std::cout << std::endl;
    }
}

void copyToFinal(int curr_path[]) {
    //#pragma omp for ordered
    for (int i = 0; i < N; i++)
        final_path[i] = curr_path[i];
    #pragma omp critical
    final_path[N] = curr_path[0];
}

int firstMin(const vector<vector<double>>& matrizDistancias, int i) {
    double min = DBL_MAX;
    int n = matrizDistancias.size();
    for (int k = 0; k < n; k++)
        if (matrizDistancias[i][k] < min && i != k)
            min = matrizDistancias[i][k];
    return min;
}

int secondMin(const vector<vector<double>>& matrizDistancias, int i) {
    double first = DBL_MAX, second = DBL_MAX;
    int n = matrizDistancias.size();
    
    for (int j = 0; j < n; j++) {
        if (i == j)
            continue;

        if (matrizDistancias[i][j] <= first) {
            second = first;
            first = matrizDistancias[i][j];
        } else if (matrizDistancias[i][j] <= second &&
                   matrizDistancias[i][j] != first)
            second = matrizDistancias[i][j];
    }
    return second;
}

void TSPRec(const vector<vector<double>>& matrizDistancias, double curr_bound, double curr_weight,
            int level, int curr_path[]) {
    if (level == N) {
        if (matrizDistancias[curr_path[level - 1]][curr_path[0]] != 0) {
            double curr_res = curr_weight + matrizDistancias[curr_path[level - 1]][curr_path[0]];

            if (curr_res < final_res) {
                
                copyToFinal(curr_path);
                final_res = curr_res;
            }
        }
        return;
    }

    int n = matrizDistancias.size();

    #pragma omp task
    for (int i = 0; i < n; i++) {
        if (matrizDistancias[curr_path[level - 1]][i] != 0 && !visited[i]) {
            double temp = curr_bound;
            curr_weight += matrizDistancias[curr_path[level - 1]][i];

            if (level == 1)
                curr_bound -= ((firstMin(matrizDistancias, curr_path[level - 1]) +
                                firstMin(matrizDistancias, i)) / 2);
            else
                curr_bound -= ((secondMin(matrizDistancias, curr_path[level - 1]) +
                                firstMin(matrizDistancias, i)) / 2);

            if (curr_bound + curr_weight < final_res) {
                curr_path[level] = i;
                visited[i] = true;
                TSPRec(matrizDistancias, curr_bound, curr_weight, level + 1, curr_path);
            }

            curr_weight -= matrizDistancias[curr_path[level - 1]][i];
            curr_bound = temp;

            memset(visited, false, sizeof(visited));
            for (int j = 0; j <= level - 1; j++)
                visited[curr_path[j]] = true;
        }
    }
    #pragma omp task wait
}

void TSP(const vector<vector<double>>& matrizDistancias) {
    int n = matrizDistancias.size();
    int curr_path[n + 1];

    int curr_bound = 0;
    memset(curr_path, -1, sizeof(curr_path));
    memset(visited, 0, sizeof(visited));

    //#pragma omp for
    for (int i = 0; i < n; i++)
        curr_bound += (firstMin(matrizDistancias, i) + secondMin(matrizDistancias, i));

    curr_bound = (curr_bound & 1) ? curr_bound / 2 + 1 : curr_bound / 2;

    visited[0] = true;
    curr_path[0] = 0;

    //#pragma omp single
    TSPRec(matrizDistancias, curr_bound, 0, 1, curr_path);
}

int main() {
    ifstream archivo("xqf10.txt");

    omp_set_num_threads(7);

    if (!archivo.is_open()) {
        cerr << "No se pudo abrir el archivo." << endl;
        return 1;
    }

    string linea;
    vector<pair<double, double>> puntos;

    while (getline(archivo, linea)) {
        if (linea.find("NODE_COORD_SECTION") != string::npos) {
            while (getline(archivo, linea) && linea != "EOF") {
                istringstream iss(linea);
                int id;
                double x, y;
                iss >> id >> x >> y;
                puntos.emplace_back(x, y);
            }
            break;
        }
    }

    archivo.close();

    vector<vector<double>> matrizDistancias = crearMatrizDistancias(puntos);

    imprimirMatrizDistancias(matrizDistancias);

    // Obtener el tiempo de inicio
    double start_time = omp_get_wtime();


     #pragma omp parallel
    {
        #pragma omp single
        //#pragma omp task
        TSP(matrizDistancias);
    }


    // Obtener el tiempo de finalización
     double end_time = omp_get_wtime();

    // Calcular la duración en milisegundos
    double duration = end_time - start_time;

    // Imprimir el tiempo transcurrido con printf
    printf("Tiempo transcurrido: %lf segundos\n", duration);


    //printf("Minimum cost: %d\n", final_res);
    printf("Minimum cost: %lf\n", final_res);

    printf("Path Taken: ");
    for (int i = 0; i <= N; i++)
        printf("%d ", final_path[i]);

    return 0;
}
