#include "../include/bdd_gradient.h"
#include "../include/read_formula.h"
#include <stdio.h>
#include <iostream>
#include <queue>
#include <set>
#include "cudd.h"
BDD::BDD(){
    this->gbm = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0); /* Initialize a new BDD manager. */
    this->roots = new std::vector<DdNode*>; 
}

BDD::BDD(Formula *formula){
    this->roots = new std::vector<DdNode*>;
    this->forward_message = new std::map<DdNode*,float>; 
    this->backward_message = new std::map<DdNode*,float>; 
    this->hipar = new std::map<DdNode*, std::vector<int>* >; 
    this->lopar = new std::map<DdNode*, std::vector<int>* >; 
    this->gbm = Cudd_Init(formula->num_of_vars, 0 ,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 0); /* Initialize a new BDD manager. */
    this->clause_weights = formula->clause_weights;
    for ( int ci = 0; ci < formula->clauses->size(); ci++){
        this->build_BDD_for_clause(formula, ci); 
    }
    //this->print_dd(this->gbm, (*roots)[16],24,4);
}

void BDD::print_dd (DdManager *gbm, DdNode *dd, int n, int pr )
{
    printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(gbm)); /*Reports the number of live nodes in BDDs and ADDs*/
    printf("DdManager vars: %d | ", Cudd_ReadSize(gbm) ); /*Returns the number of BDD variables in existence*/
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(gbm) ); /*Returns the number of times reordering has occurred*/
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(gbm) ); /*Returns the memory in use by the manager measured in bytes*/
    Cudd_PrintDebug(gbm, dd, n, pr);  // Prints to the standard output a DD and its statistics: number of nodes, number of leaves, number of minterms.
}



static int sum(std::vector<int> *v){
    int s = 0;
    for (int i = 0; i < v->size(); i++){
        s += (*v)[i];
    }
    return s;
}

void BDD::build_BDD_for_clause(Formula *formula, int ci){
    char ctype = (*formula->clause_type)[ci];
    DdNode *root;
    if ( ctype == 'c'){
        root = this->BDD_for_PB(gbm, (*formula->clauses)[ci], (*formula->coefs)[ci], (*formula->klist)[ci], 0, 0, (*formula->clauses)[ci]->size(), GE);
    }
    else if ( ctype == 'x'){
        root = this->BDD_for_XOR(gbm, (*formula->clauses)[ci], 0, 1);
    }
    else if ( ctype == 'p'){
        if ( (*formula->comparators)[ci] == GE){       
            root = this->BDD_for_PB(gbm, (*formula->clauses)[ci], (*formula->coefs)[ci], (*formula->klist)[ci], 0, 0, sum((*formula->coefs)[ci]), GE);
        }
        else{
            root = this->BDD_for_PB(gbm, (*formula->clauses)[ci], (*formula->coefs)[ci], (*formula->klist)[ci], 0, 0, sum((*formula->coefs)[ci]), EQ);
        }
    }
    this->roots->push_back(root);
}

DdNode *BDD::BDD_for_PB(DdManager *gbm, std::vector<int> *literals, std::vector<int> *coefs, int rhs, int size, int sum, int material_left, int comparator){
    DdNode *res;
    if (comparator == GE){
        if ( sum >= rhs){
            res = Cudd_ReadOne(gbm);
            Cudd_Ref(res);
            return res;
        } 
        else if ( sum + material_left < rhs){
            res = Cudd_ReadLogicZero(gbm);
            Cudd_Ref(res);
            return res;
        }
    }
    else if (comparator == EQ){
        if ((sum > rhs) || (sum + material_left < rhs) ){
            res = Cudd_ReadLogicZero(gbm);
            Cudd_Ref(res);
            return res;
        } 
        else if ( (material_left == 0) && (sum==rhs) ){
            res = Cudd_ReadOne(gbm);
            Cudd_Ref(res);
           return res;
        }
    }    

    int current_literal = (*literals)[size];
    int current_variable = abs(current_literal);

    DdNode *current_node = Cudd_bddIthVar(gbm, current_variable - 1);

    DdNode *true_child = this->BDD_for_PB(gbm, literals, coefs, rhs, size + 1, sum + (*coefs)[size], material_left - (*coefs)[size], comparator);
    DdNode *false_child = this->BDD_for_PB(gbm, literals, coefs, rhs, size + 1, sum, material_left - (*coefs)[size], comparator);
    
    if (current_literal > 0){
        res = Cudd_bddIte(gbm, current_node, true_child, false_child);
    }
    else{
        res = Cudd_bddIte(gbm, current_node, false_child, true_child);
    }
    Cudd_Ref(res);       
    return res;
}


DdNode *BDD::BDD_for_XOR(DdManager *gbm, std::vector<int> *literals, int size, int product){
    DdNode *res;
    if ( size == literals->size()){
        if (product == -1){
            res = Cudd_ReadOne(gbm);
            Cudd_Ref(res);
            return res;
        }
        else{
            res = Cudd_ReadLogicZero(gbm);
            Cudd_Ref(res);
            return res;
        }
    }
    int current_literal = (*literals)[size];
    int current_variable = abs(current_literal);
    DdNode *current_node = Cudd_bddIthVar(gbm, current_variable - 1);
    
    DdNode *true_child = this->BDD_for_XOR(gbm, literals, size + 1, -product);
    
    DdNode *false_child = this->BDD_for_XOR(gbm, literals, size + 1, product);
    if (current_literal > 0){
        res = Cudd_bddIte(gbm, current_node, true_child, false_child);
    }
    else{
        res = Cudd_bddIte(gbm, current_node, false_child, true_child);
    }
    Cudd_Ref(res);
    return res;
}
void BDD::message_clean(){

}


void BDD::forward_pass(std::vector<float> *x){
    this->message_clean();
// TO DO: 1. fill in the forward message
//        2. generate the parent map
    std::queue<DdNode*> Q;
    std::set<DdNode*> S;
    for( int i=0; i < this->roots->size(); i++){
        DdNode *root = (*roots)[i];
        Q.push(root);
        float weight = (*this->clause_weights)[i];
        (*this->forward_message)[root] = weight;
    }
    while ( !Q.empty()){
        DdNode *node = Q.front();
        Q.pop();
        int current_variable = Cudd_NodeReadIndex(node);
        DdNode *hi_child = Cudd_T(node);
        DdNode *lo_child = Cudd_E(node);
        if ( hi_child != NULL){
            Q.push(hi_child);
            if ( !this->forward_message->count(hi_child)){
                (*this->forward_message)[hi_child] = 0;
            }
            (*this->forward_message)[hi_child] += (*this->forward_message)[node] * (*x)[current_variable];
        }
        if ( lo_child != NULL){
            if ( !this->forward_message->count(lo_child)){
                (*this->forward_message)[lo_child] = 0;
            }
            (*this->forward_message)[lo_child] += (*this->forward_message)[node] * (1 - (*x)[current_variable]);
            Q.push(lo_child);
        }
   }

            
            
}

void BDD::backward_pass(std::vector<float> *x, std::vector<float> *grad){
// to be finished
}

void BDD::change_weights_in_BDD(std::vector<float> *clause_weights){
    this->clause_weights = clause_weights;
}

std::vector<float> *BDD::grad(std::vector<float> *x){
    int n = x->size();
    std::vector<float> *grad = new std::vector<float>(n);
    for( int i=0;i<n;i++ ) (*grad)[i] = 0;
    return grad;
}

float BDD::fval(std::vector<float> *x){
    int fval = 0;
    
    return fval;
}
