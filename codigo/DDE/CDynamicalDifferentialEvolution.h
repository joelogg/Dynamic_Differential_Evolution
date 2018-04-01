#ifndef C_DYNAMICAL_DIFFERENTIAL_EVOLUTION_H
#define C_DYNAMICAL_DIFFERENTIAL_EVOLUTION_H

#include <vector>
#include <iostream>
#include <omp.h>
#include <limits>
#include "RandomUniform.h"

#define ACKLEY		 0
#define SCHWEFEL	 1
#define FUNCION3    	 2
#define PI           3.14159265358979323846  

using namespace std;


int thread_count = 8;


class CDynamicalDifferentialEvolution 
{

public:    
    CDynamicalDifferentialEvolution()
    { }

public:
    
    vector<double> resta(vector<double> p1, vector<double> p2) 
    {
        vector<double> tmp;
        for (int i = 0; i < p1.size(); i++)
            tmp.push_back(p1[i] - p2[i]);
        return tmp;
    }

    vector<double> suma(vector<double> p1, vector<double> p2) 
    {
        vector<double> tmp;
        for (int i = 0; i < p1.size(); i++)
            tmp.push_back(p1[i] + p2[i]);
        return tmp;
    }

    vector<double> multiplicacion(vector<double> p, double factor) 
    {
        vector<double> tmp;
        for (int i = 0; i < p.size(); i++)
            tmp.push_back(p[i] * factor);
        return tmp;
    }
    
    double calVal(vector<double> pattern, int funcion)
    {
        if(funcion== ACKLEY)
        {
            double a = 20;
            double b = .2;
            double c = 2*PI;
            double tmp1 = 0.0;
            double tmp2 = 0.0;
            int size = pattern.size();
            
            for (int i = 0; i < pattern.size(); i++)
            {
                double val1 = pattern[i] *  pattern[i];
                tmp1 +=  val1;
		
                double val2 = cos(c *  pattern[i]);
                tmp2 +=  val2;
            }
            return -a * exp(-b * sqrt(tmp1/size)) - exp(tmp2/size) + a + exp(1);
        }
        else if(funcion == SCHWEFEL)
        {
            double tmp1 = 0.0;
            int size = pattern.size();

            for (int i = 0; i < size; i++)
            {
                double val1 = pattern[i] *  sin(sqrt(fabs(pattern[i])));
                tmp1 +=  val1;
            }
            return (418.9829 * size) - tmp1;
        }
        else if(funcion == FUNCION3)
        {
            double tmp1 = 0.0;
            
            int size = pattern.size();
            for (int i = 0; i < size; i++)
            {
                tmp1 +=  pattern[i] *  pattern[i];
            }
            return .5 - ((pow(sin(tmp1),2)-.5) / pow((1.0 + 0.001 * tmp1),2));
	
        }
        else 
            return 0.5;
        
    }
    
    void run(int dimension1, int NP, int funcion, int TYPE,
            double pmin, double pmax, double f, double gr, int iterations)
    {
        
        int dimension = dimension1;
        vector<vector<double>>  population; 
        vector<double> pOptimal = vector<double>(dimension, 0.0);

        //#pragma omp parallel for num_threads(thread_count) 
        for (int i = 0; i < NP; i++)
        {            
            vector<double> aux = vector<double>(dimension, 0.0);
            for (int i = 0; i < dimension; i++) 
            {
                aux[i] = RandomUniform(pmin, pmax); 
            }
            //#pragma omp critical
            population.push_back(aux);

        }

        get_optimal(NP, population, pOptimal, TYPE, funcion);
        cout << get_optimal_value(pOptimal, funcion) << "; ";
        
        #pragma omp parallel num_threads(thread_count) 
	for(int iter=0; iter<iterations; iter++)
        {
            
        
            #pragma omp for schedule(static,1)
            for (int i = 0; i <  NP; i++)
            {
                vector<double> tmp;
                //---------mutacion------------
                int a = (int)RandomUniform(0.0, NP-1.0);
		int b = (int)RandomUniform(0.0, (double)NP-1.0);
		tmp = suma(pOptimal , multiplicacion(resta(population[a],population[b]), f)  );
                
              
                //-------------recombinacion-----------
		for (int k = 0; k < dimension; k++)
                {
                    double rand_value = RandomUniform(0.0, 1.0);
                    if (rand_value < gr)
                    	tmp[k] = tmp[k];
                    else	
                        tmp[k] = population[i][k];
		}
                
              
                
                //----------Seleccion---------------
                if (TYPE == 0) 
                { // minimun
                    if (calVal(tmp, funcion) < calVal(population[i], funcion))
                    {
                        population[i] = tmp;
                    }
                    if (calVal(population[i], funcion) < calVal(pOptimal, funcion))
                    {
                        #pragma omp critical
                        pOptimal = population[i];					
                    }
                }
		else if (TYPE == 1) 
                { // maximum 
                    if (calVal(tmp, funcion) > calVal(population[i], funcion))
                    {
                        population[i] = tmp;
                    }
                    if (calVal(population[i], funcion) > calVal(pOptimal, funcion))
                    {
                        #pragma omp critical
                        pOptimal = population[i];
                    }
		}
            }
	
       
            //if ((iter) % 5 == 0)
                //cout << get_optimal_value(pOptimal, funcion) << "; ";
                
           
	}
	cout << get_optimal_value(pOptimal, funcion) << "; ";
        cout << endl;	
    }

    void get_optimal (int NP, vector<vector<double>>  population, vector<double> &pOptimal, int TYPE, int funcion)
    {
        double value = 0.0;
        double aux;
	if (TYPE == 0) 
        { // minimun
            value = numeric_limits<double>::max();
            for (int i = 0; i <  NP; i++)
            {
                aux = calVal(population[i], funcion);
                if (aux < value)
                {
                    value = aux;
                    pOptimal = population[i];
		}
            }
	}
        else if (TYPE == 1) 
        { // maximum
            value = numeric_limits<double>::min();
            for (int i = 0; i <  NP; i++)
            {
                aux = calVal(population[i], funcion);
                if (aux> value)
                {
                    value = aux;
                    pOptimal = population[i];
		}
            }
	}
    }

	
    double get_optimal_value(vector<double> pOptimal, int funcion)
    {
        return calVal(pOptimal, funcion);
        //return pOptimal[0];
    }
};


#endif