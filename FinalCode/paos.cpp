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
struct objeto
{
    double fx = 0;
    double fy = 0;
    double fz = 0;
    double ax = 0;
    double vx = 0;
    double x;
    double ay = 0;
    double vy = 0;
    double y;
    double az = 0;
    double vz = 0;
    double z;
    double masa;
    bool existe = true;
};

int main(int argc, char *argv[])
{
    auto begin = std::chrono::high_resolution_clock::now();
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
    init_config.open("init_config4.txt");
    init_config << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << "\n";
    objeto* objetos = new objeto[num_objects];
    int contador=0;
    //register int i;
    //register int j;
    //register int k;

    /*Creacion de cada uno de los objetos*/
    for (int i = 0; i < num_objects ; ++i) //creamos cada objeto en función al numero de objetos que hemos recibido como parámetro
    {
        objetos[i].x = uni_dis(generator);//valores de posiciones generados con una distribución uniforme
        objetos[i].y = uni_dis(generator);
        objetos[i].z = uni_dis(generator);
        objetos[i].masa = normal_dis(generator); //valores de las masas que se calculan usando una districbución normal
        init_config << std::fixed << std::setprecision(3) << objetos[i].x << " " << objetos[i].y << " " << objetos[i].z << " " << objetos[i].vx << " " << objetos[i].vy << " " << objetos[i].vz << " " << objetos[i].masa << "\n";
    }

    init_config.close();

    //Comprobamos que no haya ninguna colisión antes de comenzar con la simulación
    for (int i = 0; i < num_objects; i++)
    {
        if (objetos[i].existe)
        {
            //Comprobamos que el objeto exista, así nos evitamos calculos innecesarios
            for (int j = 0; j < num_objects; j++)
            {
                //Solo queremos calcular una colisión, ya que si el objeto 0 colisiona con el 5, no hace falta comprobar que el objeto 5 colisiona con el 0
                if ( j>i )
                {
                    //Volvemos a comprobar si el otro objeto existe
                    if (objetos[j].existe)
                    {
                    //Calculamos la distancia entre ambos objetos
                    double norm = std::sqrt((((objetos[j].x) - (objetos[i].x)) * ((objetos[j].x) - (objetos[i].x))) + (((objetos[j].y) - (objetos[i].y)) * ((objetos[j].y) - (objetos[i].y))) + (((objetos[j].z) - (objetos[i].z)) * ((objetos[j].z) - (objetos[i].z))));
                    if (norm <= 1)
                        {
                            //Fusionamos ambos objetos en uno solo, y nos quedamos siempre con el objeto de menor indice
                            objetos[i].vx += objetos[j].vx;
                            objetos[i].vy += objetos[j].vy;
                            objetos[i].vz += objetos[j].vz;
                            objetos[i].masa += objetos[j].masa;
                            //Eliminamos el objeto de mayor índice
                            objetos[j].existe = false;
                            contador+=1;
                        }
                    }
                }
                else {
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
            if (objetos[i].existe)
            {
                //Iniciamos la simulación, poniendo las fuerzas a 0, ya que después de cada iteración, las fuerzas se reinician, no son acumulativas en cada iteración
                objetos[i].fx = 0;
                objetos[i].fy = 0;
                objetos[i].fz = 0;
                //Reiniciamos también la de los hilos
                x_thread[threads_id][i] = 0;
                y_thread[threads_id][i] = 0;
                z_thread[threads_id][i] = 0;

                for (int j = 0; j < num_objects; j++)
                {
                    //Comprobamos que el objeto exista
                    if (objetos[j].existe)
                    {
                        //Calculamos la distancia entre ambos objetos
                        double norm = std::sqrt((((objetos[j].x) - (objetos[i].x)) * ((objetos[j].x) - (objetos[i].x))) + (((objetos[j].y) - (objetos[i].y)) * ((objetos[j].y) - (objetos[i].y))) + (((objetos[j].z) - (objetos[i].z)) * ((objetos[j].z) - (objetos[i].z))));
                        //Elevamos la norma al cubo para más adelante evitarnos repetir operaciones innecesariamente
                        norm = (norm*norm*norm);
                        //Guardamos una variable con estos calculos, ya que son constantes y nos evitan repetir operaciones innecesariamente
                        double gmm = (GRAVITY * (objetos[i].masa) * (objetos[j].masa));
                        if (norm != 0) 
                        {
                            //Se guarda la fuerza calculada a cada objeto de su respectivo hilo, para después sumarle el valor a los diferentes objetos
                            x_thread[threads_id][i] += (gmm * ((objetos[j].x) - (objetos[i].x))) / (norm);
                            y_thread[threads_id][i] += (gmm * ((objetos[j].y) - (objetos[i].y))) / (norm);
                            z_thread[threads_id][i] += (gmm * ((objetos[j].z) - (objetos[i].z))) / (norm);
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
                if (objetos[j].existe)
                {
                    //A cada objeto se le suma la fuerza que haya calculado el respectivo hilo
                    objetos[j].fx += x_thread[i][j];
                    objetos[j].fy += y_thread[i][j];
                    objetos[j].fz += z_thread[i][j];
                }
            }
        }

        /*Iniciamos la paralelización, para ello utilizamos hilos que se repartiran los calculos entre
        los diferentes objetos, es decir, se crearan x hilos y cada hilo calculará los objetos que se les asignen*/
        #pragma omp parallel for
        for (int i = 0; i < num_objects; i++)
        {

            //Comprobamos que el objeto exista
            if (objetos[i].existe)
            {
                //Al igual que en el calculo de fuerzas, leemos los datos de memoria del objeto en si, y guardamos el valor en los vectores
                //Para que después se sincronicen los valores de los objetos con lo calculado en cada hilo
                //Guardamos el valor de la masa en una variable, para evitarnos más accesos a la estructura
                double masa = objetos[i].masa;
                ax_thread[threads_id][i] = ((objetos[i].fx) / (masa));
                vx_thread[threads_id][i] = ((objetos[i].vx) + ((ax_thread[threads_id][i]) * (time_step)));
                px_thread[threads_id][i] = (objetos[i].x + ((vx_thread[threads_id][i]) * (time_step)));
                ay_thread[threads_id][i] = ((objetos[i].fy) / (masa));
                vy_thread[threads_id][i] = ((objetos[i].vy) + ((ay_thread[threads_id][i]) * (time_step)));
                py_thread[threads_id][i] = (objetos[i].y + ((vy_thread[threads_id][i]) * (time_step)));
                az_thread[threads_id][i] = ((objetos[i].fz) / (masa));
                vz_thread[threads_id][i] = ((objetos[i].vz) + ((az_thread[threads_id][i]) * (time_step)));
                pz_thread[threads_id][i] = (objetos[i].z + ((vz_thread[threads_id][i]) * (time_step)));
            }
        }

        //Sincronizamos los calculos obtenidos por cada hilo
        for (int i = 0; i < nthreads; ++i) 
        {
            for (int j = 0; j < num_objects; ++j) 
            {
                //Comprobamos que el objeto exista
                if (objetos[j].existe)
                {
                    //Comprobamos que el hilo i alberga los datos calculados para el objeto j (Así nos ahorramos asignaciones o sumas de valores nulos)
                    if (px_thread[i][j]!=0)
                    {
                        //Asignamos los valores calculados por cada hilo a los objetos
                        objetos[j].ax = ax_thread[i][j];
                        objetos[j].vx = vx_thread[i][j];
                        objetos[j].x = px_thread[i][j];
                        //Hay un rebote debido a una posición negativa
                        if (objetos[j].x<=0)
                        {
                            objetos[j].x = 0;
                            objetos[j].vx = ((-1)*(objetos[j].vx));
                        }
                        //Hay un rebote debido a una posición mayor a la del recinto recibido por parametro
                        else if(objetos[j].x>=size_enclosure)
                        {
                            objetos[j].x = (size_enclosure);
                            objetos[j].vx = ((-1)*(objetos[j].vx));
                        }
                        //Asignamos los valores calculados por cada hilo a los objetos
                        objetos[j].ay = ay_thread[i][j];
                        objetos[j].vy = vy_thread[i][j];
                        objetos[j].y = py_thread[i][j];
                        //Hay un rebote debido a una posición negativa
                        if (objetos[j].y<=0)
                        {
                            objetos[j].y = 0;
                            objetos[j].vy = ((-1)*(objetos[j].vy));
                        }
                        //Hay un rebote debido a una posición mayor a la del recinto recibido por parametro
                        else if(objetos[j].y>=size_enclosure)
                        {
                            objetos[j].y = (size_enclosure);
                            objetos[j].vy = ((-1)*(objetos[j].vy));
                        }
                        //Asignamos los valores calculados por cada hilo a los objetos
                        objetos[j].az = az_thread[i][j];
                        objetos[j].vz = vz_thread[i][j];
                        objetos[j].z = pz_thread[i][j];
                        //Hay un rebote debido a una posición negativa
                        if (objetos[j].z<=0)
                        {
                            objetos[j].z = 0;
                            objetos[j].vz = ((-1)*(objetos[j].vz));
                        }
                        //Hay un rebote debido a una posición mayor a la del recinto recibido por parametro
                        else if(objetos[j].z>=size_enclosure)
                        {
                            objetos[j].z = (size_enclosure);
                            objetos[j].vz = ((-1)*(objetos[j].vz));
                        }
                    }
                }
            } 
        }

        //Fase de colisiones
        for (int i = 0; i < num_objects; i++)
        {
            if (objetos[i].existe)
            {
                //Comprobamos que el objeto exista, así nos evitamos calculos innecesarios
                for (int j = 0; j < num_objects; j++)
                {
                    //Solo queremos calcular una colisión, ya que si el objeto 0 colisiona con el 5, no hace falta comprobar que el objeto 5 colisiona con el 0
                    if ( j>i )
                    {
                        //Volvemos a comprobar si el otro objeto existe
                        if (objetos[j].existe)
                        {
                            //Calculamos la distancia entre ambos objetos
                            double norm = std::sqrt((((objetos[j].x) - (objetos[i].x)) * ((objetos[j].x) - (objetos[i].x))) + (((objetos[j].y) - (objetos[i].y)) * ((objetos[j].y) - (objetos[i].y))) + (((objetos[j].z) - (objetos[i].z)) * ((objetos[j].z) - (objetos[i].z))));
                            if (norm <= 1)
                            {
                                //Fusionamos ambos objetos en uno solo, y nos quedamos siempre con el objeto de menor indice
                                objetos[i].vx += objetos[j].vx;
                                objetos[i].vy += objetos[j].vy;
                                objetos[i].vz += objetos[j].vz;
                                objetos[i].masa += objetos[j].masa;
                                //Eliminamos el objeto de mayor índice
                                objetos[j].existe = false;
                                contador+=1;
                            }
                        }
                    }
                    else {
                        //Para evitar comprobaciones innecesarias
                        j=i;
                    }
                }
            }
        }
    }

    //Escritura del txt correspondiente a los diferentes valores de cada objeto
    std::ofstream final_config;
    final_config.open("final_config4.txt");
    final_config << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << num_objects-contador << "\n";
    //Escribimos los valores de cada objeto, siempre y cuando el objeto exista
    for (int i = 0; i < num_objects; i++)
    {   if (objetos[i].existe)
        {

            final_config << std::fixed << std::setprecision(3) << objetos[i].x << " " << objetos[i].y << " "
                         << objetos[i].z << " " << objetos[i].vx << " " << objetos[i].vy << " " << objetos[i].vz << " "
                         << objetos[i].masa << "\n";
        }
    }

    final_config.close();
    //Cerramos el fichero
    //Liberamos memoria
    delete [] objetos;
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
    std::cout << "Time: "<< elapsed.count()<< " milliseconds\n";
    return 0;
}
