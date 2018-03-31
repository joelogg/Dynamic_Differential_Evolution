#ifndef C_DYNAMICAL_DIFFERENTIAL_EVOLUTION_H
#define C_DYNAMICAL_DIFFERENTIAL_EVOLUTION_H

#include <vector>
#include <iostream>
#include <omp.h>
#include <limits>
#include "RandomUniform.h"

using namespace std;


int thread_count = 8;

template <typename pattern_t, int NP, typename function_t, int TYPE> 
class CDynamicalDifferentialEvolution 
{

public:
    vector<pattern_t>  population; 
    vector<pattern_t>  noise_population; 
    vector<pattern_t>  trial_population; 
    function_t 	function;
    pattern_t 	pOptimal;

public:
    CDynamicalDifferentialEvolution(double pmin, double pmax, double f, double gr, int iterations)
    {
        #pragma omp parallel for num_threads(thread_count) private(NP, pmin, pmax) shared(population)
        for (int i = 0; i < NP; i++)
        {
            population.push_back(pattern_t(pmin, pmax));
	}
	
                
        get_optimal();
        
        #pragma omp parallel num_threads(thread_count) shared(pOptimal, iterations, iter) \
        private(TYPE, NP, gr, f, population) schedule(dynamic)
	for(int iter=0; iter<iterations; iter++)
        {
            #pragma omp for
            for (int i = 0; i <  NP; i++)
            {
                pattern_t tmp;
                //---------mutacion------------
                int a = (int)RandomUniform(0.0, (double)NP-1.0);
		int b = (int)RandomUniform(0.0, (double)NP-1.0);
		tmp = pOptimal + ((population[a] - population[b]) * f)  ;
                
                
                //-------------recombinacion-----------
		for (int k = 0; k < population[0].features.size(); k++)
                {
                    double rand_value = RandomUniform(0.0, 1.0);
                    if (rand_value < gr)
                    	tmp.features[k] = tmp.features[k];
                    else	
                        tmp.features[k] = population[i].features[k];
		}
                
                
                
                //----------Seleccion---------------
                if (TYPE == 0) 
                { // minimun
                    if (function(tmp) < function(population[i]))
                    {
                        population[i] = tmp;
                    }
                    if (function(population[i]) < function(pOptimal))
                    {
                        //#pragma omp critical
                        pOptimal = population[i];					
                    }
                }
		else if (TYPE == 1) 
                { // maximum 
                    if (function(tmp) > function(population[i]))
                    {
                        population[i] = tmp;
                    }
                    if (function(population[i]) > function(pOptimal))
                    {
                        //#pragma omp critical
                        pOptimal = population[i];
                    }
		}
            }
	
       
            if ((iter) % 5 == 0)
                cout << get_optimal_value() << "; ";
                
            
	}
	
        cout << endl;		
    }

    void get_optimal (void)
    {
        double value = 0.0;
	if (TYPE == 0) 
        { // minimun
            value = numeric_limits<double>::max();
                for (int i = 0; i <  NP; i++)
                {
                    if (function(population[i]) < value)
                    {
                        value = function(population[i]);
			pOptimal = population[i];
			}
		}
	}
        else if (TYPE == 1) 
        { // maximum
            value = numeric_limits<double>::min();
            for (int i = 0; i <  NP; i++)
            {
                if (function(population[i]) > value)
                {
                    value = function(population[i]);
                    pOptimal = population[i];
		}
            }
	}
    }

	
    double get_optimal_value(void)
    {
        return function(pOptimal);
    }
};


#endif