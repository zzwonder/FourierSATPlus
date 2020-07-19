#ifndef BDD_GRADIENT_H
#define BDD_GRADIENT_H
#include "../include/read_formula.h"
#include<string>
#include "cudd.h"
#include <sys/types.h>
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <iostream>


class BDD{
    public:
    DdManager *gbm;
    std::vector<float> *clause_weights;
    std::map<DdNode*,float> *forward_message;
    std::map<DdNode*, float> *backward_message;  
    std::map<DdNode*, std::vector<DdNode*> * > *hipar;  
    std::map<DdNode*, std::vector<DdNode*> * > *lopar;  
    std::vector<DdNode*> *roots;

    BDD();
    BDD(Formula *formula);
    void change_weights_in_BDD(std::vector<float> *clause_weights);
    
    std::vector<float> *grad(std::vector<float> *x);
    float fval(std::vector<float> *x);
    void print_dd (DdManager *gbm, DdNode *dd, int n, int pr );   
 
    private:
    int generate_parents;
    void message_clean();
    void backward_message_clean();
    void build_BDD_for_clause(Formula *formula, int ci);
    void forward_pass(std::vector<float> *x);
    void backward_pass(std::vector<float> *x, std::vector<float> *grad);
    
    DdNode *BDD_for_PB(DdManager *gbm, std::vector<int> *literals, std::vector<int> *coefs, int rhs, int size, int sum, int material_left, int comparator, std::map<std::pair<int,int>, DdNode *> *nodes_storage);
    DdNode *BDD_for_XOR(DdManager *gbm, std::vector<int> *literals, int size, int product, std::map<std::pair<int,int>, DdNode *> *nodes_storage);
};
#endif
