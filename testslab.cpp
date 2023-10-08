#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <chrono>

// Open header named "testslab.h"
//#include "testslab.h"

using namespace std;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // Semilla para el generador de numeros aleatorios
mt19937 gen(seed);                                                           // Generador de numeros aleatorios, variable global

double rand_uniforme(double a, double b)
{
    uniform_real_distribution<double> dis(a, b);
    return dis(gen);
}

double modulo(vector<double> v)
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

vector<double> normalizar(vector<double> v)
{
    double mod = modulo(v);
    vector<double> v_norm;
    v_norm.push_back(v[0] / mod);
    v_norm.push_back(v[1] / mod);
    v_norm.push_back(v[2] / mod);
    return v_norm;
}

class Neutron
{
public:
    vector<double> posicion(void)
    {
        vector<double> pos = {x, y, z};
        return pos;
    }
    vector<double> direccion(void)
    {
        vector<double> dir = {mu, 0.0, 0.0};
        dir.at(1) = sqrt(1 - mu * mu) * cos(phi);
        dir.at(2) = sqrt(1 - mu * mu) * sin(phi);
        return normalizar(dir);
    }
    double peso(void)
    {
        return w;
    }
    void transportar(double Sigma_s, double Sigma_a)
    {
        double camino = path(Sigma_s, Sigma_a);
        change_dir();
        vector<double> dir = direccion();
        x = x + camino * direccion()[0];
        y = y + camino * direccion()[1];
        z = z + camino * direccion()[2];
    }
    void set_peso(double w_nuevo)
    {
        w = w_nuevo;
    }
    int interaccion(double Sigma_s, double Sigma_a) // 1 scattering, 2 absorcion
    {
        double r = rand_uniforme(0, 1);
        if (r < Sigma_s / (Sigma_s + Sigma_a))
        {
            return 1; // scattering
        }

        return 2; // absorsion
    }
    Neutron(double a, double b) : x(rand_uniforme(a, b)), y(rand_uniforme(a, b)), z(rand_uniforme(a, b)),
                                  mu(rand_uniforme(-1, 1)), phi(rand_uniforme(0, 2 * (3.1415926))), w(1.0) {}
    ~Neutron() {}

private:
    double path(double Sigma_s, double Sigma_a)
    {
        return -log(rand_uniforme(0, 1)) / (Sigma_s + Sigma_a);
    }
    void change_dir(void)
    {
        mu = rand_uniforme(-1, 1);
        phi = rand_uniforme(0, 2 * (3.141592));
    }
    double x, y, z, mu, phi, w; // coordenadas cartesianas, dirección, peso,se podría incluir la energía
};

// funcion que verifica si el neutron esta dentro del slab
bool dentro_slab(Neutron n, double a)
{
    if (abs(n.posicion()[0]) <= a)
    {
        return true;
    }
    return false;
}

// funcion que asigna la posicion del neutron en la grilla (solo en x), -a<x<a
int posicion_grilla(double x, double a, int n)
{
    double delta = 2 * a / n;
    float i = 0.0;
    while (x > -a + i * delta)
    {
        i++;
    }
    return i - 1;
}

vector<double> norm_tali(vector<double> v, int n, double delta)
{
    for (int i = 0; i < v.size(); i++)
    {
        v.at(i) = v.at(i) / (n * delta); // normalizo el tali por el numero de neutrones y el ancho de la celda
    }
    return v;
}

int main()
{
    double a = 15.0;       // [cm] slab de ancho 2a
    double Sigma_s = 0.5;  // [cm^-1]
    double Sigma_a = 0.05; // [cm^-1]
    double w_c = 0.5;      // peso de corte

    int N = 10000; // numero de neutrones
    // discretizacion del slab
    int n = 100; // numero de celdas
    double delta = (2 * a) / n;

    vector<double> tali_flux(n, 0.0); // tali de flujo
    vector<double> tali_int(n, 0.0);      // tali de interaccion
    vector<double> tali_error(n, 0.0); // tali de error

    auto start_time = std::chrono::high_resolution_clock::now(); // Start time

    for (int i = 0; i < N; i++)
    {
        // Estimate time to finish
        if (i == N / 1000)
        {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto time_elapsed = std::chrono::duration_cast<std::chrono::microseconds>(current_time - start_time).count();
            int estimate_time_second = round(static_cast<double>(time_elapsed) * (N / (i + 1)) - 1) / 1e6;
            std::cout << "Estimated time to finish: " << estimate_time_second << " seconds" << std::endl;
        }

        // Show progress in a bar of 10 steps growing one step each 10 %
        if (i % (N / 10) == 0 && i != 0)
        {
            cout << "Progress: ";
            for (int j = 0; j < i / (N / 10); j++)
            {
                cout << "#";
            }
            for (int j = i / (N / 10); j < 10; j++)
            {
                cout << "-";
            }
            cout << endl;
        }

        Neutron n1(-a, a);
        int j = 0;
        while (j < 1000)
        {
            n1.transportar(Sigma_s, Sigma_a);
            if (!dentro_slab(n1, a))
            {
                break;
            } 
            tali_flux.at(posicion_grilla(n1.posicion()[0], a, n)) += n1.peso() / (Sigma_s + Sigma_a);
            tali_int.at(posicion_grilla(n1.posicion()[0], a, n)) += 1.0;
            tali_error.at(posicion_grilla(n1.posicion()[0], a, n)) += n1.peso() / (Sigma_s + Sigma_a) * n1.peso() / (Sigma_s + Sigma_a);
            n1.set_peso(n1.peso() * Sigma_s / (Sigma_s + Sigma_a)); // método de absorción implicita
            if (n1.peso() < w_c) // método ruleta rusa
            {
                if (rand_uniforme(0, 1) < 1 - w_c)
                {
                    break;
                }
                n1.set_peso(n1.peso() / (1 - w_c));
            }
            j++;
        }
    }

    // normalizo los talis
    tali_flux = norm_tali(tali_flux, n, delta);
    tali_int = norm_tali(tali_int, n, delta);
    for (int i = 0; i < n; i++)  
    {
        tali_error.at(i) = sqrt(tali_error.at(i));
    }
    tali_error = norm_tali(tali_error, n, delta);
     
    // importo el tali a un archivo de texto
    ofstream datos;
    datos.open("datos.txt");
    if (!datos.is_open())
    {
        cout << "No se pudo abrir el archivo" << endl;
        return 1;
    }
    // escribo los datos en el archivo
    for (int i = 0; i < n; i++)
    {
        datos << -a + i * delta << "," << tali_flux.at(i) << "," << tali_int.at(i) << "," << tali_error.at(i) << endl;
    }
    datos.close();

    return 0;
}
