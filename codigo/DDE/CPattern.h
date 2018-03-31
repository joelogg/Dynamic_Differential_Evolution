#ifndef CPATTERN_H
#define CPATTERN_H

#include <iostream>
#include <vector>
#include "RandomUniform.h"

using namespace std;

template <typename T, int SIZE>
class CPattern {
public:
    vector<T> features;

public:

    CPattern(const CPattern& pattern) 
    {
        features = pattern.features;
    }

    CPattern(double MIN = 0.0, double MAX = 0.0) 
    {
        features = vector<T>(SIZE, 0.0);
        for (int i = 0; i < features.size(); i++) 
        {
            features[i] = RandomUniform(MIN, MAX); //(double)MIN + RandomUniform(0.0, SCOPE) /*rand()%(int)SCOPE*/ /(double)SCOPE * (double)(MAX  - MIN);
        }
    }

    CPattern operator-(CPattern p) 
    {
        CPattern tmp;
        for (int i = 0; i < features.size(); i++)
            tmp.features[i] = features[i] - p.features[i];
        return tmp;
    }

    CPattern operator+(CPattern p) {
        CPattern tmp;
        for (int i = 0; i < features.size(); i++)
            tmp.features[i] = features[i] + p.features[i];
        return tmp;
    }

    CPattern operator*(double _factor) 
    {
        CPattern tmp;
        for (int i = 0; i < features.size(); i++)
            tmp.features[i] = features[i] * _factor;
        return tmp;
    }

    /*template <typename _T, int _SIZE, int _SCOPE>
    friend ostream& operator<< (ostream& out, CPattern<_T,_SIZE,_SCOPE> _pattern){
            for (int i = 0; i < _pattern._features.size(); i++){
                    out << _pattern._features[i] << " ";
            }
            return out;
    }*/

};

#endif