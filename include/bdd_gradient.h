#ifndef BDD_GRADIENT_H
#define BDD_GRADIENT_H
#include "../include/read_formula.h"
#include<string>
#include "../cudd/"
vector<float> *grad(Formula formula, vector<float> *clause_weights, bdd, vector<float> *x);
float fval(Formula formula,vector<float> *clause_weights, bdd, vector<float> *x);


#endif
