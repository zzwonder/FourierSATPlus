#ifndef OPTIMIZE_H
#define OPTIMIZE_H
#include<vector>
#include<string>
#include<dlib/optimization.h>
#include<iostream>
#include "bdd_gradient.h"
typedef dlib::matrix<double,0,1> column_vector;

class Optimizer{
    public:
    int num_of_vars;
    BDD *bdd; 
        std::vector<double> *x;    

    Optimizer();
    Optimizer(int n, BDD *bdd);
    double minimize();
    double fval(const column_vector& m);
    column_vector a_grad (const column_vector& m);
};


#endif
