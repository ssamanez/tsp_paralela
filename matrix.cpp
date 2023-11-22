#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

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

int main() {
    std::ifstream archivo("xqf5.txt");

    if (!archivo.is_open()) {
        std::cerr << "No se pudo abrir el archivo." << std::endl;
        return 1;
    }

    std::string linea;
    std::vector<std::pair<double, double>> puntos;

    // Leer el archivo línea por línea
    while (std::getline(archivo, linea)) {
        if (linea.find("NODE_COORD_SECTION") != std::string::npos) {
            // Comenzar a leer las coordenadas
            while (std::getline(archivo, linea) && linea != "EOF") {
                std::istringstream iss(linea);
                int id;
                double x, y;
                iss >> id >> x >> y;
                puntos.emplace_back(x, y);
            }
            break;
        }
    }

    archivo.close();

    // Crear la matriz de distancias
    std::vector<std::vector<double>> matrizDistancias = crearMatrizDistancias(puntos);

    // Imprimir la matriz de distancias
    imprimirMatrizDistancias(matrizDistancias);

    return 0;
}
