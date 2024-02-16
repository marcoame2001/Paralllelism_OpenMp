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
#include <omp.h>


/*Constantes*/
const double MEAN = 1E21; //media de la distribución normal
const double STDDEV = 1E15;//desviación típica de la distribución normal
const double GRAVITY = 6.674E-11; //Constante de gravitación

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

    for (int i = 0; i < num_objects; i++)
    {
        if (objetos[i].existe)
        {
            for (int j = 0; j < num_objects; j++)
            {
                if ( j>i )
                {
                    if (objetos[j].existe)
                    {
                    double norm = std::sqrt((((objetos[j].x) - (objetos[i].x)) * ((objetos[j].x) - (objetos[i].x))) + (((objetos[j].y) - (objetos[i].y)) * ((objetos[j].y) - (objetos[i].y))) + (((objetos[j].z) - (objetos[i].z)) * ((objetos[j].z) - (objetos[i].z))));
                    if (norm <= 1)
                        {
                            objetos[i].vx += objetos[j].vx;
                            objetos[i].vy += objetos[j].vy;
                            objetos[i].vz += objetos[j].vz;
                            objetos[i].masa += objetos[j].masa;
                            objetos[j].existe = false;
                            contador+=1;
                        }
                    }
                }
                else {
                    j=i;
                }
            }
        }
    }
    int nthreads;
    #pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    int threads_id;
    
    // Paralelizamos
    #pragma omp parallel
    {
        threads_id = omp_get_thread_num();
    }

    /*Calculo de fuerzas entre objetos*/
    for (int k = 0; k < num_iterations; k++)
    {
        //std::cout << "ITERACION: " << k << "\n";
        for (int i = 0; i < num_objects; i++)
        {
            if (objetos[i].existe)
            {
                objetos[i].fx = 0;
                objetos[i].fy = 0;
                objetos[i].fz = 0;
                double sum_fx = 0.0;
                double norm = 0.0;
                double gmm = 0.0;

                omp_set_num_threads(4);
                #pragma omp parallel for reduction(+:sum_fx) private(norm,gmm)
                for (int j = 0; j < num_objects; j++)
                {
                    if (objetos[j].existe)
                    {   
                        norm = std::sqrt((((objetos[j].x) - (objetos[i].x)) * ((objetos[j].x) - (objetos[i].x))) + (((objetos[j].y) - (objetos[i].y)) * ((objetos[j].y) - (objetos[i].y))) + (((objetos[j].z) - (objetos[i].z)) * ((objetos[j].z) - (objetos[i].z))));
                        norm = (norm*norm*norm);
                        gmm = (GRAVITY * (objetos[i].masa) * (objetos[j].masa));
                        if (norm != 0) 
                        { 
                            sum_fx += (gmm * ((objetos[j].x) - (objetos[i].x))) / (norm);
                            objetos[i].fy += (gmm * ((objetos[j].y) - (objetos[i].y))) / (norm);
                            objetos[i].fz += (gmm * ((objetos[j].z) - (objetos[i].z))) / (norm);
                        }
                    }
                    //gravitatory_force(&objetos[i].masa, &objetos[j].masa, &objetos[i].x, &objetos[j].x, &objetos[i].y, &objetos[j].y, &objetos[i].z, &objetos[j].z, &objetos[i].fx, &objetos[i].fy, &objetos[i].fz);
                }
                objetos[i].fx = sum_fx;
            }
        }

        for (int i = 0; i < num_objects; i++)
        {
            if (objetos[i].existe)
            {   
                double masa = objetos[i].masa;
                objetos[i].ax = ((objetos[i].fx) / (masa));
                objetos[i].vx = ((objetos[i].vx) + ((objetos[i].ax) * (time_step)));
                objetos[i].x = (objetos[i].x +((objetos[i].vx) * (time_step)));
                objetos[i].ay = ((objetos[i].fy) / (masa));
                objetos[i].vy = ((objetos[i].vy) + ((objetos[i].ay) * (time_step)));
                objetos[i].y = (objetos[i].y + ((objetos[i].vy) * (time_step)));
                objetos[i].az = ((objetos[i].fz) / (masa));
                objetos[i].vz = ((objetos[i].vz) + ((objetos[i].az) * (time_step)));
                objetos[i].z = (objetos[i].z + ((objetos[i].vz) * (time_step)));
                //aceleration_vector(&objetos[i].ax, &objetos[i].ay, &objetos[i].az, &objetos[i].fx, &objetos[i].fy, &objetos[i].fz, &objetos[i].masa);
                //speed_vector(&objetos[i].vx, &objetos[i].vy, &objetos[i].vz, &objetos[i].ax, &objetos[i].ay, &objetos[i].az, &time_step);
                //position_vector(&objetos[i].x, &objetos[i].y, &objetos[i].z, &objetos[i].vx, &objetos[i].vy, &objetos[i].vz, &time_step);
                if (objetos[i].x<=0)
                {
                    objetos[i].x = 0;
                    objetos[i].vx = ((-1)*(objetos[i].vx));
                }
                else if(objetos[i].x>=size_enclosure)
                {
                    objetos[i].x = (size_enclosure);
                    objetos[i].vx = ((-1)*(objetos[i].vx));
                }
                if (objetos[i].y<=0)
                {
                    objetos[i].y = 0;
                    objetos[i].vy = ((-1)*(objetos[i].vy));
                }
                else if(objetos[i].y>=size_enclosure)
                {
                    objetos[i].y = (size_enclosure);
                    objetos[i].vy = ((-1)*(objetos[i].vy));
                }
                if (objetos[i].z<=0)
                {
                    objetos[i].z = 0;
                    objetos[i].vz = ((-1)*(objetos[i].vz));
                }
                else if(objetos[i].z>=size_enclosure)
                {
                    objetos[i].z = (size_enclosure);
                    objetos[i].vz = ((-1)*(objetos[i].vz));
                }
                //bounce(&objetos[i].x, &objetos[i].y, &objetos[i].z, &objetos[i].vx, &objetos[i].vy, &objetos[i].vz, &size_enclosure);
                /*
                std::cout << "Time step: " << time_step << "\n";
                std::cout << "Masa objeto " << i << ": " << objetos[i].masa << "\n";
                std::cout << "Fuerza X del objeto " << i << ": " << objetos[i].fx << "\n";
                std::cout << "Fuerza Y del objeto " << i << ": " << objetos[i].fy << "\n";
                std::cout << "Fuerza Z del objeto " << i << ": " << objetos[i].fz << "\n";
                std::cout << "Aceleracion X del objeto " << i << ": " << objetos[i].ax << "\n";
                std::cout << "Aceleracion Y del objeto " << i << ": " << objetos[i].ay << "\n";
                std::cout << "Aceleracion Z del objeto " << i << ": " << objetos[i].az << "\n";
                std::cout << "Velocidad X del objeto " << i << ": " << objetos[i].vx << "\n";
                std::cout << "Velocidad Y del objeto " << i << ": " << objetos[i].vy << "\n";
                std::cout << "Velocidad Z del objeto " << i << ": " << objetos[i].vz << "\n";
                std::cout << "Posicion X del objeto " << i << ": " << objetos[i].x << "\n";
                std::cout << "Posicion Y del objeto " << i << ": " << objetos[i].y << "\n";
                std::cout << "Posicion Z del objeto " << i << ": " << objetos[i].z << "\n";
                std::cout << "-------------------------------------------------------------\n";
                */
            }
        }
        
        for (int i = 0; i < num_objects; i++)
        {
            if (objetos[i].existe)
            {
                for (int j = 0; j < num_objects; j++)
                {
                    if ( j>i )
                    {
                        if (objetos[j].existe)
                        {
                            double norm = std::sqrt((((objetos[j].x) - (objetos[i].x)) * ((objetos[j].x) - (objetos[i].x))) + (((objetos[j].y) - (objetos[i].y)) * ((objetos[j].y) - (objetos[i].y))) + (((objetos[j].z) - (objetos[i].z)) * ((objetos[j].z) - (objetos[i].z))));
                            if (norm <= 1 && i!=j)
                            {
                                objetos[i].vx = objetos[j].vx + objetos[i].vx;
                                objetos[i].vy =  objetos[j].vy + objetos[i].vy;
                                objetos[i].vz = objetos[j].vz + objetos[i].vz;
                                objetos[i].masa = objetos[j].masa + objetos[i].masa;
                                objetos[j].existe = false;
                                contador+=1;
                            }
                        }
                    }
                    else
                    {
                        j=i;
                    }
                }
            }
        }
    }

    std::ofstream final_config;
    final_config.open("final_config4.txt");
    final_config << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << num_objects-contador << "\n";
    for (int i = 0; i < num_objects; i++)
    {   if (objetos[i].existe)
        {

            final_config << std::fixed << std::setprecision(3) << objetos[i].x << " " << objetos[i].y << " "
                         << objetos[i].z << " " << objetos[i].vx << " " << objetos[i].vy << " " << objetos[i].vz << " "
                         << objetos[i].masa << "\n";
        }
    }

    final_config.close();
    delete [] objetos;
    return 0;
}
