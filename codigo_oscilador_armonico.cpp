#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

// --- FUNCIÓN DE POTENCIAL ---

double PotencialArmonico(int j, int N, double h, double L) {
    const double omega = 200.0;    // Frecuencia del oscilador armónico
    const double x_centro = L / 2.0;   // Centro del potencial 
    
    double x = j * h; // Calculamos la posición física real dado j
    double dx = x - x_centro;
    
    // V(x) = (omega^2 / 4) * (x - x_centro)^2 
    double V = (omega * omega / 4.0) * dx * dx;
    
    // Reescalamos el potencial para tener V_tilde    
    return V * h * h;
}


//---- EFECTO STARK -----
double PotencialStark(int j, int N, double h, double L, double lambda_stark) {
    const double omega = 200.0;
    const double x_centro = L / 2.0;
    
    double x = j * h; // Calculamos la posición física real dado j
    double dx = x - x_centro;
    
    // V(x) = (omega^2 / 4) * dx^2 - lambda * dx
    double V = (omega * omega / 4.0) * dx * dx - lambda_stark * dx;
    
    return V * h * h;
}


// --- FUNCIONES DE INICIALIZACIÓN DE LA FUNCIÓN DE ONDA ---

// ----Inicializar con autofuncion ------- 

void InicializarConAutofuncion(vector<complex<double>>& phi, int N, double h, int n_cuantico, double omega, double L) {
    double factor_escala = sqrt(omega / 2.0);

    for (int j = 0; j <= N; ++j) {
        double x = j * h; // Calculamos la posición física real dado j
        double x_centro = L / 2.0; 
        double dx = x - x_centro;

        // Hacemos un cambio de variable y= raiz(m omega / hbar) * dx  para trabajar con H_n(y)   
        double y = factor_escala * dx;
        double gauss = exp(- y * y / 2.0);
        
        double H_n = 0.0; // Inicializamos a cero
        
        // Definimos las funciones de hermite para n=0 y n=1 y establecemos la relación de recurrencia
        if (n_cuantico == 0) {
            H_n = 1.0;
        }
        else if (n_cuantico == 1) {
            H_n = 2.0 * y;
        }
        else {
            // Necesitamos H0 y H1 para iniciar la relación de recurrencia
            double H_menos2 = 1.0;      
            double H_menos1 = 2.0 * y;  

            // Relación de recurrencia: H_n(y) = 2y*H_{n-1}(y) - 2(n-1)*H_{n-2}(y)

            for (int k = 2; k <= n_cuantico; ++k) {
                H_n = 2.0 * y * H_menos1 - 2.0 * (k - 1) * H_menos2;
                H_menos2 = H_menos1;
                H_menos1 = H_n;
            }
        }
        
        // La funcion de onda es Hn * gaussiana  
        phi[j] = H_n * gauss;
        
        // Condiciones de contorno
        if (j == 0 || j == N) phi[j] = 0.0;
    }
}


// ---- Inicializar con gaussiana -------
void InicializarConGaussiana(vector<complex<double>>& phi, int N, double h, double x0, double sigma) {
    for (int j = 0; j <= N; ++j) {
        double x = j * h; // Calculamos la posición física real dado j
        
        // forma de la gaussiana: exp(- (x - x0)^2 / (2 * sigma^2))
        phi[j] = exp(-pow(x - x0, 2) / (2.0 * sigma * sigma));
        
        // Condiciones de contorno
        if (j == 0 || j == N) {
            phi[j] = 0.0;
        }
    }
}


// --- CÁLCULO DE LA NORMA ----
double CalcularNorma(const vector<complex<double>>& phi, int N, double h) {
    double integral_norma = 0.0;
    for (int j = 0; j <= N; ++j) {
        integral_norma += norm(phi[j]);
    }
    return integral_norma * h;
}


//---- OPERADOR MOMENTO -----
complex<double> OperadorMomento(int j, const vector<complex<double>>& phi, double h) {
    complex<double> i(0, 1);
    return -i * (phi[j + 1] - phi[j - 1]) / (2.0 * h); // Usamos la derivada central 
}

//---- OPERADOR MOMENTO CUADRADO -----

complex<double> OperadorMomentoCuadrado(int j, const vector<complex<double>>& phi, double h) {
    return -(phi[j + 1] - 2.0 * phi[j] + phi[j - 1]) / (h * h); // Usamos la segunda derivada central 
}

// --- CÁLCULO DE LOS VALORES ESPERADOS ---

//--- Valor esperado de x ----  
double CalcularEsperadoX(const vector<complex<double>>& phi, int N, double h) {
    double esp_x = 0.0;
    for (int j = 0; j <= N; ++j) {
        double x = j * h; // Calculamos la posición física real dado j
        esp_x += x * norm(phi[j]) ;
    }
    return esp_x * h;
}


//--- Valor esperado de x^2----
double CalcularEsperadoX2(const vector<complex<double>>& phi, int N, double h) {
    double esp_x2 = 0.0;
    for (int j = 0; j <= N; ++j) {
        double x = j * h;
        esp_x2 += x * x * norm(phi[j]) ;
    }
    return esp_x2 * h;
}


//--- Valor esperado de p----
double CalcularEsperadoP(const vector<complex<double>>& phi, int N, double h) {
    complex<double> esp_p = 0.0;
    for (int j = 1; j < N; ++j) {
        esp_p += conj(phi[j]) * OperadorMomento(j, phi, h) * h; //  (Phi*) * P * Phi  
    }
    return esp_p.real(); 
}

//--- Valor esperado de p^2----
double CalcularEsperadoP2(const vector<complex<double>>& phi, int N, double h) {
    complex<double> esp_p2 = 0.0;
    for (int j = 1; j < N; ++j) {
        esp_p2 += conj(phi[j]) * OperadorMomentoCuadrado(j, phi, h) * h; //  (Phi*) * P^2 * Phi  
    }
    return esp_p2.real();
}

//--- Valor esperado de la energía sin el efecto Stark -----
double CalcularEsperadoEnergia(const vector<complex<double>>& phi, int N, double h, double omega, double L) {
    
    // Calculamos la parte de la energia potencial    
    double esp_V = 0.0;
    for (int j = 0; j <= N; ++j) {
        double x = j * h; // Calculamos la posición física real dado j
        double x_centro = L / 2.0;
        double V = (omega * omega / 4.0) * pow(x - x_centro, 2); // Necesitamos el potencial    
        esp_V += V * norm(phi[j]) * h;
    }
    // Calculamos la parte de la energia cinetica
    double esp_T = CalcularEsperadoP2(phi, N, h); 
    return esp_T + esp_V;
}

// --- Valor esperado de la energia con el efecto Stark -----
double CalcularEsperadoEnergiaStark(const vector<complex<double>>& phi, int N, double h, double omega, double lambda_stark, double L) {
    double esp_V = 0.0;
    double x_centro = L / 2.0;
    for (int j = 0; j <= N; ++j) {
        double x = j * h;
        double dx = x - x_centro;
        
        // Añadimos el efecto de la perturbacion al potencial (-lambda * x)
        double V_cont = (omega * omega / 4.0) * dx * dx - lambda_stark * dx; 
        esp_V += V_cont * norm(phi[j]) * h;
    }
    double esp_T = CalcularEsperadoP2(phi, N, h); 
    return esp_T + esp_V;
}


int main() {
    
    // --- PARÁMETROS DE LA RED ---
    const double L = 0.8;           // Longitud de la caja 
    const int N = 1000;             // Intervalos en los que la vamos a dividir  
    const double h = L / N;         // Paso espacial 
    const double s = 0.0001;        // Paso temporal   
    
    const double PI = acos(-1.0); 
    double s_tilde = s / (h * h);   // Reescalamos el paso temporal

    complex<double> i(0, 1);

    // Vectores del sistema
    vector<complex<double>> phi(N + 1);
    vector<double> V_tilde(N + 1);
    vector<complex<double>> alpha(N);
    vector<complex<double>> beta(N);
    vector<complex<double>> chi(N + 1);
    vector<complex<double>> A0(N, 0.0);

    const double omega = 200.0; 
   
    //------ Parámetros del Efecto Stark -----
    const double lambda_stark = 50.0; // Intensidad del campo electrico   
    const int t_encendido = 1000;     // Paso temporal donde se activa el campo



    // -----Elegimos paquete de ondas inicial ------
    
    // Inicializamos con una autofuncion 
    InicializarConAutofuncion(phi, N, h, 0, omega, L); 


    // Inicializamos con una Gaussiana  
    //InicializarConGaussiana(phi, N, h, 0.3, 1.0 / 16.0); // Podemos elegir x0 y sigma


    // Inicializaciamos el potencial    
    for (int j = 0; j <= N; ++j) {
        V_tilde[j] = PotencialArmonico(j, N, h, L);
    }


    // Normalización de la función de onda elegida
    double suma_norma = 0.0;
    for (int j = 0; j <= N; ++j) {
        suma_norma += norm(phi[j]);
    }
    double norma = sqrt(suma_norma * h);
    for (int j = 0; j <= N; ++j) {
        phi[j] = phi[j] / norma;
    }

    // CÁLCULO DE A0 Y ALPHA
    for (int j = 1; j < N; ++j) {
        A0[j] = complex<double>(-2.0 - V_tilde[j], 2.0 / s_tilde);
    }

    // Calculamos alpha usando A0
    alpha[N - 1] = 0.0; 
    for (int j = N - 1; j >= 1; --j) {
        alpha[j - 1] = -1.0 / (A0[j] + alpha[j]); 
    }


    // Definimos los archivos que vamos a usar para guardar los datos    
    ofstream file_animacion("evolucion_probabilidad_stark50.dat"); // Archivo para guardar la evolución temporal de la probabilidad
    ofstream file_fases("evolucion_fases_stark50.dat");  // Archivo para guardar la evolución temporal de las fases
    ofstream file_observables("observables_stark50.dat"); // Archivo para guardar los valores esperados de los observables

    file_observables << "# t  Norma  <x>  <p>  <H>  DeltaX*DeltaP\n";

    // --- BUCLE DE TIEMPO ---
    for (int t = 0; t <= 5000; ++t) { 

       
        // Encendido del campo electrico del efecto stark
            /*
            if (t == t_encendido) {
                // Actualizamos el potencial para todo el espacio
                for (int j = 0; j <= N; ++j) {
                    V_tilde[j] = PotencialStark(j, N, h, L, lambda_stark);
                }
                
                // Al cambiar V_tilde, tenemos que recalcular A0 y alpha    
                for (int j = 1; j < N; ++j) {
                    A0[j] = complex<double>(-2.0 - V_tilde[j], 2.0 / s_tilde);
                }
                alpha[N - 1] = 0.0; 
                for (int j = N - 1; j >= 1; --j) {
                    alpha[j - 1] = -1.0 / (A0[j] + alpha[j]); 
                }
                cout << "Campo electrico encendido en el paso t = " << t << " (Tiempo: " << t * s << ")" << endl;
            }
            */
        

        double integral_norma = CalcularNorma(phi, N, h);
        double esp_x = CalcularEsperadoX(phi, N, h);
        double esp_x2 = CalcularEsperadoX2(phi, N, h);
        double esp_p = CalcularEsperadoP(phi, N, h);
        double esp_p2 = CalcularEsperadoP2(phi, N, h);
        
        //Si no hay efecto Stark calculamos la energia asi:
        double esp_H = CalcularEsperadoEnergia(phi, N, h, omega, L);

        
        /*
        // Si hay efecto Stark calculamos la energia asi:
        double esp_H = 0.0;
            if (t >= t_encendido) {
                esp_H = CalcularEsperadoEnergiaStark(phi, N, h, omega, lambda_stark, L);
            } else {
                esp_H = CalcularEsperadoEnergia(phi, N, h, omega, L);
            }
        */
        

        // Calculamos la incertidumbre usando los valores esperados
        double delta_x = sqrt(esp_x2 - esp_x * esp_x); 
        double delta_p = sqrt(esp_p2 - esp_p * esp_p);

        // Guardamos los datos de 50 en 50 pasos para la animacion (densidad de probabilidad y fases)
        if (t % 50 == 0) {
            for (int j = 0; j <= N; ++j) {
                file_animacion << norm(phi[j]) << " ";
                file_fases << phi[j].real() << " " << phi[j].imag() << " ";
            }
            file_animacion << "\n";
            file_fases << "\n";
        }

        file_observables << fixed << setprecision(15) << t * s << " " << integral_norma << " " 
                         << esp_x << " "  << esp_p << " " 
                         << esp_H << " "  << delta_x * delta_p << "\n";


        // --- ALGORITMO CRANK-NICOLSON ---
        beta[N - 1] = 0.0;
        for (int j = N - 1; j >= 1; --j) {
            complex<double> b = ((4.0 * i) / s_tilde) * phi[j]; 
            beta[j - 1] = (b - beta[j]) / (A0[j] + alpha[j]);   
        }

        chi[0] = 0.0; 
        for (int j = 0; j < N; ++j) {
            chi[j + 1] = alpha[j] * chi[j] + beta[j]; 
        }

        for (int j = 0; j <= N; ++j) {
            phi[j] = chi[j] - phi[j];
        }
    }

    file_animacion.close();
    file_fases.close();
    file_observables.close();
    
    cout << "Simulacion finalizada." << endl;
    return 0;
}
