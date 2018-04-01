#include <iostream>
#include <math.h>
#include <time.h>

using namespace std;

#include "CPattern.h"
#include "CDynamicalDifferentialEvolution.h"



#define PI           3.14159265358979323846  
#define MINIMUM 	 0
#define MAXIMUM		 1
#define ACKLEY		 0
#define SCHWEFEL	 1
#define FUNCION3    	 2




int main (int, char*[])
{
        srand(time(NULL));
	int nRep = 100;
        
        CDynamicalDifferentialEvolution evolution_t;
	while (nRep-- > 0)
	{
            cout<<"Repeticion:  "<< (100-nRep)<<endl;
		
		
            //----------Ackley-------
            double  max_value = 32.7680;
            double  min_value = -32.7680;
            // d=2
            //evolution_t.run(2,200, ACKLEY, MINIMUM,
            //            min_value, max_value, .2, .7, 50); 
            // d=8
            evolution_t.run(8,2000, ACKLEY, MINIMUM,
                        min_value, max_value, .2, .7, 50); 
            
            
            
            /*/----------Schwefel----------
            double  max_value = 500;
            double  min_value = -500;
            // d=8
            //evolution_t.run(2,200, SCHWEFEL, MINIMUM,
            //            min_value, max_value, .2, .6, 50); // for Schwefel d8
            // d=8
            //evolution_t.run(8,20000, SCHWEFEL, MINIMUM,
            //            min_value, max_value, .01, .4, 50); // for Schwefel d8
            */
	
            
            /*/----------Funcion 3----------
            double  max_value= 100;
            double  min_value= -100;
            // d=8
            //evolution_t.run(2,200, FUNCION3, MAXIMUM,
            //            min_value, max_value, .2, .8, 50); // for Schwefel d8
            // d=8
            evolution_t.run(8,20000, FUNCION3, MAXIMUM,
                        min_value, max_value, .2, .8, 50); // for Schwefel d8
            */
            
            
            
	}
	cout << "terminado";
	
	//getchar();
	return true;
}
