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
    register int i;
    register int j;
    register int k;

    /*Creacion de cada uno de los objetos*/
    for (i = 0; i < num_objects ; ++i) //creamos cada objeto en función al numero de objetos que hemos recibido como parámetro
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

    for (i = 0; i < num_objects; i++)
    {
        if (objetos.existe[i])
        {
            for (j = 0; j < num_objects; j++)
            {
                if ( j>i )
                {
                    if (objetos.existe[j])
                    {
                        double norm = std::sqrt((((objetos.x[j]) - (objetos.x[i])) * ((objetos.x[j]) - (objetos.x[i]))) + (((objetos.y[j]) - (objetos.y[i])) * ((objetos.y[j]) - (objetos.y[i]))) + (((objetos.z[j]) - (objetos.z[i])) * ((objetos.z[j]) - (objetos.z[i]))));
                        if (norm <= 1)
                        {
                            objetos.vx[i] += objetos.vx[j];
                            objetos.vy[i] += objetos.vy[j];
                            objetos.vz[i] += objetos.vz[j];
                            objetos.masa[i] += objetos.masa[j];
                            objetos.existe[j] = false;
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



    for (k = 0; k < num_iterations; k++)
    {
        //std::cout << "ITERACION: " << k << "\n";

        for (i = 0; i < num_objects; i++)
        {
            if (objetos.existe[i])
            {
                objetos.fx[i] = 0;
                objetos.fy[i]  = 0;
                objetos.fz[i]  = 0;

                for (j = 0; j < num_objects; j++)
                {
                    if (objetos.existe[j])
                    {
                        double norm = std::sqrt((((objetos.x[j]) - (objetos.x[i])) * ((objetos.x[j]) - (objetos.x[i]))) + (((objetos.y[j]) - (objetos.y[i])) * ((objetos.y[j]) - (objetos.y[i]))) + (((objetos.z[j]) - (objetos.z[i])) * ((objetos.z[j]) - (objetos.z[i]))));
                        norm = (norm*norm*norm);
                        double gmm = (GRAVITY * (objetos.masa[i]) * (objetos.masa[j]));
                        if (norm != 0) {
                            objetos.fx[i] += (gmm * ((objetos.x[j]) - (objetos.x[i]))) / (norm);
                            objetos.fy[i] += (gmm * ((objetos.y[j]) - (objetos.y[i]))) / (norm);
                            objetos.fz[i] += (gmm * ((objetos.z[j]) - (objetos.z[i]))) / (norm);
                        }
                    }
                    //gravitatory_force(&objetos[i].masa, &objetos[j].masa, &objetos[i].x, &objetos[j].x, &objetos[i].y, &objetos[j].y, &objetos[i].z, &objetos[j].z, &objetos[i].fx, &objetos[i].fy, &objetos[i].fz);
                }
            }
        }

        for (i = 0; i < num_objects; i++)
        {
            if (objetos.existe[i])
            {
                double masa = objetos.masa[i];
                objetos.ax[i] = ((objetos.fx[i]) / masa);
                objetos.ay[i] = ((objetos.fy[i]) / masa);
                objetos.az[i] = ((objetos.fz[i]) / masa);
                //aceleration_vector(&objetos[i].ax, &objetos[i].ay, &objetos[i].az, &objetos[i].fx, &objetos[i].fy, &objetos[i].fz, &objetos[i].masa);
                objetos.vx[i] = ((objetos.vx[i]) + ((objetos.ax[i]) * (time_step)));
                objetos.vy[i] = ((objetos.vy[i]) + ((objetos.ay[i]) * (time_step)));
                objetos.vz[i] = ((objetos.vz[i]) + ((objetos.az[i]) * (time_step)));
                //speed_vector(&objetos[i].vx, &objetos[i].vy, &objetos[i].vz, &objetos[i].ax, &objetos[i].ay, &objetos[i].az, &time_step);
                objetos.x[i] += ((objetos.vx[i]) * (time_step));
                objetos.y[i] += ((objetos.vy[i]) * (time_step));
                objetos.z[i] += ((objetos.vz[i]) * (time_step));
                //position_vector(&objetos[i].x, &objetos[i].y, &objetos[i].z, &objetos[i].vx, &objetos[i].vy, &objetos[i].vz, &time_step);
                if (objetos.x[i]<=0)
                {
                    objetos.x[i] = 0;
                    objetos.vx[i] = ((-1)*(objetos.vx[i]));
                }
                if(objetos.x[i]>=size_enclosure)
                {
                    objetos.x[i] = (size_enclosure);
                    objetos.vx[i] = ((-1)*(objetos.vx[i]));
                }
                if (objetos.y[i]<=0)
                {
                    objetos.y[i] = 0;
                    objetos.vy[i] = ((-1)*(objetos.vy[i]));
                }
                if(objetos.y[i]>=size_enclosure)
                {
                    objetos.y[i] = (size_enclosure);
                    objetos.vy[i] = ((-1)*(objetos.vy[i]));
                }
                if (objetos.z[i]<=0)
                {
                    objetos.z[i] = 0;
                    objetos.vz[i] = ((-1)*(objetos.vz[i]));
                }
                if(objetos.z[i]>=size_enclosure)
                {
                    objetos.z[i] = (size_enclosure);
                    objetos.vz[i] = ((-1)*(objetos.vz[i]));
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
        for (i = 0; i < num_objects; i++)
        {
            if (objetos.existe[i])
            {
                for (j = 0; j < num_objects; j++)
                {
                    if ( j>i )
                    {
                        if (objetos.existe[j])
                        {
                            double norm = std::sqrt((((objetos.x[j]) - (objetos.x[i])) * ((objetos.x[j]) - (objetos.x[i]))) + (((objetos.y[j]) - (objetos.y[i])) * ((objetos.y[j]) - (objetos.y[i]))) + (((objetos.z[j]) - (objetos.z[i])) * ((objetos.z[j]) - (objetos.z[i]))));
                            if (norm <= 1)
                            {
                                objetos.vx[i] += objetos.vx[j];
                                objetos.vy[i] += objetos.vy[j];
                                objetos.vz[i] += objetos.vz[j];
                                objetos.masa[i] += objetos.masa[j];
                                objetos.existe[j] = false;
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
    final_config.open("final_config2.txt");
    final_config << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << num_objects-contador << "\n";
    for (i = 0; i < num_objects; i++)
    {   if (objetos.existe[i] == true)
        {

            final_config << std::fixed << std::setprecision(3) << objetos.x[i] << " " << objetos.y[i] << " "
                         << objetos.z[i] << " " << objetos.vx[i] << " " << objetos.vy[i] << " " << objetos.vz[i] << " "
                         << objetos.masa[i] << "\n";
        }
    }

    final_config.close();
    return 0;
}

