/*
 * Programación orientada al rendimiento C++
 * Grupo 82, Equipo 5
 * Marco Antonio Merola Nia:100413665
 * Adrian Perez Ruiz Nia: 100429044
 * Victor Herranz Sanchez Nia: 100429052
 * Alejandro Gonzalez Nuñez Nia: 10042913
 */

/*Librerias incluidas*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <vector>
#include <chrono>
#include <omp.h>

/*Constantes*/
const double MEAN = 1E21; //media de la distribución normal
const double STDDEV = 1E15;//desviación típica de la distribución normal
const double GRAVITY = 6.674E-11; //Constante de gravitación

//Estructura que contiene los elementos necesarios para caracterizar un objeto
struct vectores
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    std::vector<double> ax;
    std::vector<double> ay;
    std::vector<double> az;
    std::vector<double> fx;
    std::vector<double> fy;
    std::vector<double> fz;
    std::vector<double> masa;
    std::vector<bool> existe;
};

int main(int argc, char *argv[])
{
    auto begin = std::chrono::high_resolution_clock::now();
    /*Comprobar el número de parametros introducidos*/
    if (argc < 6)
    {
        std::cerr << "Error: Wrong number of parameters.\n "<< argv[0]<< " invoked with "<< argc-1<< " parameters. \n";
        return -1;
    }

    int num_objects = std::stoi(argv[1]); //Número de objetos
    int num_iterations = std::stoi(argv[2]); //Número de iteraciones (pasos de tiempo) a simular
    int random_seed = std::stoi(argv[3]); //Semilla para las funciones generadoras de números aleatorios
    double size_enclosure = std::stod(argv[4]); //Es el tamaño de cada lado del cubo perfecto con origen en 0
    double time_step = std::stod(argv[5]); //indica el incremento del tiempo de cada interación

    /*Comprobar si los parametros introducidos son correctos*/
    if (num_objects <= 0 )
    {
        std::cerr << "Error: Invalid number of objects.\n "<< argv[0]<< " invoked with "<< argc-1<< " parameters. \n";
        return -2;
    }
    else if (num_iterations <= 0)
    {
        std::cerr << "Error: Invalid number of iterations.\n "<< argv[0]<< " invoked with "<< argc-1<< " parameters. \n";
        return -2;
    }

    else if (random_seed <= 0)
    {
        std::cerr << "Error: Invalid random seed.\n "<< argv[0]<< " invoked with "<< argc-1<< " parameters. \n";
        return -2;
    }

    else if (size_enclosure <= 0)
    {
        std::cerr << "Error: Invalid size enclosure.\n "<< argv[0]<< " invoked with "<< argc-1<< " parameters. \n";
        return -2;
    }

    else if (time_step <= 0)
    {
        std::cerr << "Error: Invalid delta time.\n "<< argv[0]<< " invoked with "<< argc-1<< " parameters. \n";
        return -2;
    }

    /*Generador de números aleatorios a partir de una semilla*/
    std::mt19937_64 generator(random_seed);

    /*Distribucion uniforme*/
    std::uniform_real_distribution <> uni_dis(0.0,size_enclosure);

    /*Distribucion normal*/
    std::normal_distribution <> normal_dis(MEAN,STDDEV);

    /*Creación de un TXT*/
    std::ofstream init_config;
    init_config.open("init_config2.txt");
    init_config << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << "\n";
    struct vectores objetos;
    int contador=0;

    /*Creacion de cada uno de los objetos*/
    for (int i = 0; i < num_objects ; ++i) //creamos cada objeto en función al numero de objetos que hemos recibido como parámetro
    {
        double x = uni_dis(generator);
        double y = uni_dis(generator);
        double z = uni_dis(generator);
        double  masa = normal_dis(generator);
        objetos.x.push_back(x);//valores de posiciones generados con una distribución uniforme
        objetos.y.push_back(y);
        objetos.z.push_back(z);
        objetos.vx.push_back(0);
        objetos.vy.push_back(0);
        objetos.vz.push_back(0);
        objetos.ax.push_back(0);
        objetos.ay.push_back(0);
        objetos.az.push_back(0);
        objetos.fx.push_back(0);
        objetos.fy.push_back(0);
        objetos.fz.push_back(0);
        objetos.existe.push_back(true);
        objetos.masa.push_back(masa); //valores de las masas que se calculan usando una districbución normal
        init_config << std::fixed << std::setprecision(3) << objetos.x[i] << " " << objetos.y[i] << " " << objetos.z[i] << " " << objetos.vx[i] << " " << objetos.vy[i] << " " << objetos.vz[i] << " " << objetos.masa[i] << "\n";

    }
    init_config.close();

    //Comprobamos que no haya ninguna colisión antes de comenzar con la simulación
    for (int i = 0; i < num_objects; i++)
    {
        //Comprobamos que el objeto exista, así nos evitamos calculos innecesarios
        if (objetos.existe[i])
        {
            for (int j = 0; j < num_objects; j++)
            {
                //Solo queremos calcular una colisión, ya que si el objeto 0 colisiona con el 5, no hace falta comprobar que el objeto 5 colisiona con el 0
                if ( j>i )
                {
                    //Volvemos a comprobar si el otro objeto existe
                    if (objetos.existe[j])
                    {
                        //Calculamos la distancia entre ambos objetos
                        double norm = std::sqrt((((objetos.x[j]) - (objetos.x[i])) * ((objetos.x[j]) - (objetos.x[i]))) + (((objetos.y[j]) - (objetos.y[i])) * ((objetos.y[j]) - (objetos.y[i]))) + (((objetos.z[j]) - (objetos.z[i])) * ((objetos.z[j]) - (objetos.z[i]))));
                        if (norm <= 1)
                        {
                            //Fusionamos ambos objetos en uno solo, y nos quedamos siempre con el objeto de menor indice
                            objetos.vx[i] += objetos.vx[j];
                            objetos.vy[i] += objetos.vy[j];
                            objetos.vz[i] += objetos.vz[j];
                            objetos.masa[i] += objetos.masa[j];
                            //Eliminamos el objeto de mayor índice
                            objetos.existe[j] = false;
                            contador+=1;
                        }
                    }
                }
                else
                {
                    //Para evitar comprobaciones innecesarias
                    j=i;
                }
            }
        }
    }


    int nthreads; //Numero de hilos
    int threads_id; //Identificador de cada hilo

    omp_set_num_threads(16); //Creamos tantos hilos como indique el numero que esta entre parentesis. Tiene que ser potencia de 2

    //Creamos los hilos
    #pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }

    #pragma omp parallel
    {
        threads_id = omp_get_thread_num();
    }

    // Vector auxiliar para sincronizar las fuerzas de cada hilo
    std::vector<std::vector<double>> x_thread(nthreads, std::vector<double>(num_objects,0.0));
    std::vector<std::vector<double>> y_thread(nthreads, std::vector<double>(num_objects,0.0));
    std::vector<std::vector<double>> z_thread(nthreads, std::vector<double>(num_objects,0.0));

    // Vector auxiliar para sincronizar las aceleraciones de cada hilo
    std::vector<std::vector<double>> ax_thread(nthreads, std::vector<double>(num_objects,0.0));
    std::vector<std::vector<double>> ay_thread(nthreads, std::vector<double>(num_objects,0.0));
    std::vector<std::vector<double>> az_thread(nthreads, std::vector<double>(num_objects,0.0));

    // Vector auxiliar para sincronizar las velocidades de cada hilo
    std::vector<std::vector<double>> vx_thread(nthreads, std::vector<double>(num_objects,0.0));
    std::vector<std::vector<double>> vy_thread(nthreads, std::vector<double>(num_objects,0.0));
    std::vector<std::vector<double>> vz_thread(nthreads, std::vector<double>(num_objects,0.0));

    // Vector auxiliar para sincronizar las posiciones de cada hilo
    std::vector<std::vector<double>> px_thread(nthreads, std::vector<double>(num_objects,0.0));
    std::vector<std::vector<double>> py_thread(nthreads, std::vector<double>(num_objects,0.0));
    std::vector<std::vector<double>> pz_thread(nthreads, std::vector<double>(num_objects,0.0));

    //Comienzo de la simulación de los objetos
    for (int k = 0; k < num_iterations; k++)
    {
        //Calculo de fuerzas
        /*Iniciamos la paralelización, para ello utilizamos hilos que se repartiran las fuerzas a calcular entre
        los diferentes objetos, es decir, se crearan x hilos y cada hilo calculará los objetos que se les asignen*/
        #pragma omp parallel for
        for (int i = 0; i < num_objects; i++)
        {
            //Comprobamos que el objeto exista
            if (objetos.existe[i])
            {
                //Iniciamos la simulación, poniendo las fuerzas a 0, ya que después de cada iteración, las fuerzas se reinician, no son acumulativas en cada iteración
                objetos.fx[i] = 0;
                objetos.fy[i]  = 0;
                objetos.fz[i]  = 0;
                //Reiniciamos también la de los hilos
                x_thread[threads_id][i] = 0;
                y_thread[threads_id][i] = 0;
                z_thread[threads_id][i] = 0;

                for (int j = 0; j < num_objects; j++)
                {
                    //Comprobamos que el objeto exista
                    if (objetos.existe[j])
                    {
                        //Calculamos la distancia entre ambos objetos
                        double norm = std::sqrt((((objetos.x[j]) - (objetos.x[i])) * ((objetos.x[j]) - (objetos.x[i]))) + (((objetos.y[j]) - (objetos.y[i])) * ((objetos.y[j]) - (objetos.y[i]))) + (((objetos.z[j]) - (objetos.z[i])) * ((objetos.z[j]) - (objetos.z[i]))));
                        //Elevamos la norma al cubo para más adelante evitarnos repetir operaciones innecesariamente
                        norm = (norm*norm*norm);
                        //Guardamos una variable con estos calculos, ya que son constantes y nos evitan repetir operaciones innecesariamente
                        double gmm = (GRAVITY * (objetos.masa[i]) * (objetos.masa[j]));
                        if (norm != 0)
                        {
                            //Se guarda la fuerza calculada a cada objeto de su respectivo hilo, para después sumarle el valor a los diferentes objetos
                            x_thread[threads_id][i] += (gmm * ((objetos.x[j]) - (objetos.x[i]))) / (norm);
                            y_thread[threads_id][i] += (gmm * ((objetos.y[j]) - (objetos.y[i]))) / (norm);
                            z_thread[threads_id][i] += (gmm * ((objetos.z[j]) - (objetos.z[i]))) / (norm);
                        }
                    }
                }
            }
        }

        //Sincronizacion de los resultados
        for (int i = 0; i < nthreads; ++i)
        {
            for (int j = 0; j < num_objects; ++j)
            {
                //Comprobamos que el objeto exista
                if (objetos.existe[j])
                {
                    //A cada objeto se le suma la fuerza que haya calculado el respectivo hilo
                    objetos.fx[j] += x_thread[i][j];
                    objetos.fy[j] += y_thread[i][j];
                    objetos.fz[j] += z_thread[i][j];
                }
            }
        }

        /*Iniciamos la paralelización, para ello utilizamos hilos que se repartiran los calculos entre
        los diferentes objetos, es decir, se crearan x hilos y cada hilo calculará los objetos que se les asignen*/
        #pragma omp parallel for
        for (int i = 0; i < num_objects; i++)
        {
            //Comprobamos que el objeto exista
            if (objetos.existe[i])
            {
                //Al igual que en el calculo de fuerzas, leemos los datos de memoria del objeto en si, y guardamos el valor en los vectores
                //Para que después se sincronicen los valores de los objetos con lo calculado en cada hilo
                //Guardamos el valor de la masa en una variable, para evitarnos más accesos a la estructura
                double masa = objetos.masa[i];
                ax_thread[threads_id][i] = ((objetos.fx[i]) / (masa));
                vx_thread[threads_id][i] = ((objetos.vx[i]) + ((ax_thread[threads_id][i]) * (time_step)));
                px_thread[threads_id][i] = (objetos.x[i] + ((vx_thread[threads_id][i]) * (time_step)));
                ay_thread[threads_id][i] = ((objetos.fy[i]) / (masa));
                vy_thread[threads_id][i] = ((objetos.vy[i]) + ((ay_thread[threads_id][i]) * (time_step)));
                py_thread[threads_id][i] = (objetos.y[i] + ((vy_thread[threads_id][i]) * (time_step)));
                az_thread[threads_id][i] = ((objetos.fz[i]) / (masa));
                vz_thread[threads_id][i] = ((objetos.vz[i]) + ((az_thread[threads_id][i]) * (time_step)));
                pz_thread[threads_id][i] = (objetos.z[i] + ((vz_thread[threads_id][i]) * (time_step)));
            }
        }

        //Sincronizamos los calculos obtenidos por cada hilo
        for (int i = 0; i < nthreads; ++i)
        {
            for (int j = 0; j < num_objects; ++j)
            {
                //Comprobamos que el objeto exista
                if (objetos.existe[j])
                {
                    //Comprobamos que el hilo i alberga los datos calculados para el objeto j (Así nos ahorramos asignaciones o sumas de valores nulos)
                    if (px_thread[i][j]!=0)
                    {
                        //Asignamos los valores calculados por cada hilo a los objetos
                        objetos.ax[j] = ax_thread[i][j];
                        objetos.vx[j] = vx_thread[i][j];
                        objetos.x[j] = px_thread[i][j];
                        //Hay un rebote debido a una posición negativa
                        if (objetos.x[j]<=0)
                        {
                            objetos.x[j] = 0;
                            objetos.vx[j] = ((-1)*(objetos.vx[j]));
                        }
                        //Hay un rebote debido a una posición mayor a la del recinto recibido por parametro
                        else if(objetos.x[j]>=size_enclosure)
                        {
                            objetos.x[j] = (size_enclosure);
                            objetos.vx[j] = ((-1)*(objetos.vx[j]));
                        }
                        //Asignamos los valores calculados por cada hilo a los objetos
                        objetos.ay[j] = ay_thread[i][j];
                        objetos.vy[j] = vy_thread[i][j];
                        objetos.y[j] = py_thread[i][j];
                        //Hay un rebote debido a una posición negativa
                        if (objetos.y[j]<=0)
                        {
                            objetos.y[j] = 0;
                            objetos.vy[j] = ((-1)*(objetos.vy[j]));
                        }
                        //Hay un rebote debido a una posición mayor a la del recinto recibido por parametro
                        else if(objetos.y[j]>=size_enclosure)
                        {
                            objetos.y[j] = (size_enclosure);
                            objetos.vy[j] = ((-1)*(objetos.vy[j]));
                        }
                        //Asignamos los valores calculados por cada hilo a los objetos
                        objetos.az[j] = az_thread[i][j];
                        objetos.vz[j] = vz_thread[i][j];
                        objetos.z[j] = pz_thread[i][j];
                        //Hay un rebote debido a una posición negativa
                        if (objetos.z[j]<=0)
                        {
                            objetos.z[j] = 0;
                            objetos.vz[j] = ((-1)*(objetos.vz[j]));
                        }
                        //Hay un rebote debido a una posición mayor a la del recinto recibido por parametro
                        else if(objetos.z[j]>=size_enclosure)
                        {
                            objetos.z[j] = (size_enclosure);
                            objetos.vz[j] = ((-1)*(objetos.vz[j]));
                        }
                    }
                }
            }
        }

        //Fase de colisiones
        for (int i = 0; i < num_objects; i++)
        {
            //Comprobamos que el objeto exista, así nos evitamos calculos innecesarios
            if (objetos.existe[i])
            {
                for (int j = 0; j < num_objects; j++)
                {
                    //Solo queremos calcular una colisión, ya que si el objeto 0 colisiona con el 5, no hace falta comprobar que el objeto 5 colisiona con el 0
                    if ( j>i )
                    {
                        //Volvemos a comprobar si el otro objeto existe
                        if (objetos.existe[j])
                        {
                            //Calculamos la distancia entre ambos objetos
                            double norm = std::sqrt((((objetos.x[j]) - (objetos.x[i])) * ((objetos.x[j]) - (objetos.x[i]))) + (((objetos.y[j]) - (objetos.y[i])) * ((objetos.y[j]) - (objetos.y[i]))) + (((objetos.z[j]) - (objetos.z[i])) * ((objetos.z[j]) - (objetos.z[i]))));
                            if (norm <= 1)
                            {
                                //Fusionamos ambos objetos en uno solo, y nos quedamos siempre con el objeto de menor indice
                                objetos.vx[i] += objetos.vx[j];
                                objetos.vy[i] += objetos.vy[j];
                                objetos.vz[i] += objetos.vz[j];
                                objetos.masa[i] += objetos.masa[j];
                                //Eliminamos el objeto de mayor índice
                                objetos.existe[j] = false;
                                contador+=1;
                            }
                        }
                    }
                    else
                    {
                        //Para evitar comprobaciones innecesarias
                        j=i;
                    }
                }
            }
        }
    }

    //Escritura del txt correspondiente a los diferentes valores de cada objeto
    std::ofstream final_config;
    final_config.open("final_config2.txt");
    final_config << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << num_objects-contador << "\n";
    //Escribimos los valores de cada objeto, siempre y cuando el objeto exista
    for (int i = 0; i < num_objects; i++)
    {
        if (objetos.existe[i] == true)
        {

            final_config << std::fixed << std::setprecision(3) << objetos.x[i] << " " << objetos.y[i] << " "
                         << objetos.z[i] << " " << objetos.vx[i] << " " << objetos.vy[i] << " " << objetos.vz[i] << " "
                         << objetos.masa[i] << "\n";
        }
    }

    final_config.close();
    //Cerramos el fichero
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
    std::cout << "Time: "<< elapsed.count()<< " milliseconds\n";
    return 0;
}

