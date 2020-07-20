#include "../include/bdd_gradient.h"
#include "../include/read_formula.h"
#include <stdio.h>
#include <iostream>
#include <queue>
#include <set>
#include "cudd.h"

class Comp_forward
{
public:
    bool operator()(DdNode *node1, DdNode *node2)
    {
        return Cudd_NodeReadIndex(node1) > Cudd_NodeReadIndex(node2);
    }
};

class Comp_backward
{
public:
    bool operator()(DdNode *node1, DdNode *node2)
    {
        return Cudd_NodeReadIndex(node1) < Cudd_NodeReadIndex(node2);
    }
};

static int sum(std::vector<int> *v){
    int s = 0;
    for (int i = 0; i < v->size(); i++){
        s += (*v)[i];
    }
    return s;
}

static int sum(std::vector<double> *v){
    double s = 0;
    for (int i = 0; i < v->size(); i++){
        s += (*v)[i];
    }
    return s;
}


BDD::BDD(){
    this->gbm = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0); /* Initialize a new BDD manager. */
    this->roots = new std::vector<DdNode*>;
    this->generate_parents = 0; 
}

BDD::BDD(Formula *formula){
    this->roots = new std::vector<DdNode*>;
    this->forward_message = new std::map<DdNode*,double>; 
    this->backward_message = new std::map<DdNode*,double>; 
    this->hipar = new std::map<DdNode*, std::vector<DdNode*>* >; 
    this->lopar = new std::map<DdNode*, std::vector<DdNode*>* >; 
    this->gbm = Cudd_Init(formula->num_of_vars, 0 ,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS, 0); /* Initialize a new BDD manager. */
    this->clause_weights = formula->clause_weights;
    this->num_of_vars = formula->num_of_vars;
    for ( int ci = 0; ci < formula->clauses->size(); ci++){
        this->build_BDD_for_clause(formula, ci); 
    }
    this->generate_parents = 0;
    this->sum_of_clause_weights = sum(clause_weights);
    std::cout<<"sum of weights = "<< this->sum_of_clause_weights<<std::endl;
    std::cout<<"building BDDs: finished."<<std::endl;
}

void BDD::print_dd (DdManager *gbm, DdNode *dd, int n, int pr )
{
    printf("DdManager nodes: %ld | ", Cudd_ReadNodeCount(gbm)); /*Reports the number of live nodes in BDDs and ADDs*/
    printf("DdManager vars: %d | ", Cudd_ReadSize(gbm) ); /*Returns the number of BDD variables in existence*/
    printf("DdManager reorderings: %d | ", Cudd_ReadReorderings(gbm) ); /*Returns the number of times reordering has occurred*/
    printf("DdManager memory: %ld \n", Cudd_ReadMemoryInUse(gbm) ); /*Returns the memory in use by the manager measured in bytes*/
    Cudd_PrintDebug(gbm, dd, n, pr);  // Prints to the standard output a DD and its statistics: number of nodes, number of leaves, number of minterms.
}


void BDD::build_BDD_for_clause(Formula *formula, int ci){
    std::map<std::pair<int,int>, DdNode *> nodes_storage;
    char ctype = (*formula->clause_type)[ci];
    DdNode *root;
    if ( ctype == 'c'){
        root = this->BDD_for_PB(gbm, (*formula->clauses)[ci], (*formula->coefs)[ci], (*formula->klist)[ci], 0, 0, (*formula->clauses)[ci]->size(), GE, &nodes_storage );
    }
    else if ( ctype == 'x'){
        root = this->BDD_for_XOR(gbm, (*formula->clauses)[ci], 0, 1, &nodes_storage);
    }
    else if ( ctype == 'p'){
        if ( (*formula->comparators)[ci] == GE){       
            root = this->BDD_for_PB(gbm, (*formula->clauses)[ci], (*formula->coefs)[ci], (*formula->klist)[ci], 0, 0, sum((*formula->coefs)[ci]), GE, &nodes_storage);
        }
        else{
            root = this->BDD_for_PB(gbm, (*formula->clauses)[ci], (*formula->coefs)[ci], (*formula->klist)[ci], 0, 0, sum((*formula->coefs)[ci]), EQ, &nodes_storage);
        }
    }
    this->roots->push_back(root);
}

DdNode *BDD::BDD_for_PB(DdManager *gbm, std::vector<int> *literals, std::vector<int> *coefs, int rhs, int size, int sum, int material_left, int comparator, std::map<std::pair<int,int>, DdNode *> *nodes_storage){
    DdNode *res;
    std::pair<int,int> size_sum_pair = std::make_pair(size,sum);
    if ( nodes_storage->count(size_sum_pair)){
        res = nodes_storage->find(size_sum_pair)->second;
        Cudd_Ref(res);
        return res;
    }  
    if (comparator == GE){
        if ( sum >= rhs){
            res = Cudd_ReadOne(gbm);
            Cudd_Ref(res);
            return res;
        } 
        else if ( sum + material_left < rhs){
            //res = Cudd_ReadLogicZero(gbm);
            res = Cudd_ReadZero(gbm);
            Cudd_Ref(res);
            return res;
        }
    }
    else if (comparator == EQ){
        if ((sum > rhs) || (sum + material_left < rhs) ){
            //res = Cudd_ReadLogicZero(gbm);
            res = Cudd_ReadZero(gbm);
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

    DdNode *current_node = Cudd_addIthVar(gbm, current_variable - 1);
    DdNode *true_child = this->BDD_for_PB(gbm, literals, coefs, rhs, size + 1, sum + (*coefs)[size], material_left - (*coefs)[size], comparator, nodes_storage);
    DdNode *false_child = this->BDD_for_PB(gbm, literals, coefs, rhs, size + 1, sum, material_left - (*coefs)[size], comparator, nodes_storage);
    Cudd_Ref(true_child); 
    Cudd_Ref(current_node); 
    Cudd_Ref(false_child); 
    if (current_literal > 0){
        res = Cudd_addIte(gbm, current_node, true_child, false_child);
    }
    else{
        res = Cudd_addIte(gbm, current_node, false_child, true_child);
    }
    Cudd_Ref(res);
    (*nodes_storage)[size_sum_pair] = res;      
    return res;
}


DdNode *BDD::BDD_for_XOR(DdManager *gbm, std::vector<int> *literals, int size, int product,  std::map<std::pair<int,int>, DdNode *> *nodes_storage){
    DdNode *res;
    std::pair<int,int> size_product_pair = std::make_pair(size,product);
    if ( nodes_storage->count(size_product_pair)){
        res = nodes_storage->find(size_product_pair)->second;
        return res;
    }

    if ( size == literals->size()){
        if (product == -1){
            res = Cudd_ReadOne(gbm);
            Cudd_Ref(res);
            return res;
        }
        else{
            res = Cudd_ReadZero(gbm);
            Cudd_Ref(res);
            return res;
        }
    }
    int current_literal = (*literals)[size];
    int current_variable = abs(current_literal);
    DdNode *current_node = Cudd_addIthVar(gbm, current_variable - 1);
    
    DdNode *true_child = this->BDD_for_XOR(gbm, literals, size + 1, -product, nodes_storage);
    
    DdNode *false_child = this->BDD_for_XOR(gbm, literals, size + 1, product, nodes_storage);
    if (current_literal > 0){
        res = Cudd_addIte(gbm, current_node, true_child, false_child);
    }
    else{
        res = Cudd_addIte(gbm, current_node, false_child, true_child);
    }
    Cudd_Ref(res);
    (*nodes_storage)[size_product_pair] = res;      
    return res;
}
void BDD::message_clean(){
    this->forward_message->clear();
}

void BDD::backward_message_clean(){
    this->backward_message->clear();
}

void BDD::forward_pass(std::vector<double> *x){
// TO DO: 1. fill in the forward message
//        2. generate the parent map
    this->message_clean();
    std::priority_queue<DdNode*, std::vector<DdNode*>, Comp_forward> Q;
    std::set<DdNode*> S;
    for( int i=0; i < this->roots->size(); i++){
        DdNode *root = (*roots)[i];
        if (!S.count(root)){
            Q.push(root);
            S.insert(root);
        }
        double weight = (*this->clause_weights)[i];
        if ( !this->forward_message->count(root)){
                (*this->forward_message)[root] = 0;
        }
        (*this->forward_message)[root] += weight;
    }
    while ( !Q.empty()){
        DdNode *node = Q.top();
        Q.pop();
        int current_variable = Cudd_NodeReadIndex(node);
        if ( Cudd_IsConstant(node)) continue;
        DdNode *hi_child = Cudd_T(node);
        DdNode *lo_child = Cudd_E(node);
        
        if ( !this->generate_parents){
            if ( !this->hipar->count(hi_child)){
                (*this->hipar)[hi_child] = new std::vector<DdNode*>;
            }
            if ( !this->lopar->count(lo_child)){
                (*this->lopar)[lo_child] = new std::vector<DdNode*>;
            }  
            (*this->hipar)[hi_child]->push_back(node);
            (*this->lopar)[lo_child]->push_back(node);
        } 
     
        double hi_message = (*this->forward_message)[node] * (*x)[current_variable];
        double lo_message = (*this->forward_message)[node] - hi_message;
        if ( hi_child != NULL){
            if ( !this->forward_message->count(hi_child)){
                (*this->forward_message)[hi_child] = 0;
            }
            (*this->forward_message)[hi_child] += hi_message;
            if (!S.count(hi_child)){
                Q.push(hi_child);
                S.insert(hi_child);
            }
        }
        if ( lo_child != NULL){
            if ( !this->forward_message->count(lo_child)){
                (*this->forward_message)[lo_child] = 0;
            }
            (*this->forward_message)[lo_child] += lo_message;
            if (!S.count(lo_child)){
                Q.push(lo_child);
                S.insert(lo_child);
            }
        }
   }
   this->generate_parents = 1;
}

void BDD::backward_pass(std::vector<double> *x, std::vector<double> *grad){
// to be finished
    this->backward_message->clear(); 
    DdNode *one = Cudd_ReadOne(this->gbm);
    (*this->backward_message)[one] = 1;
    std::priority_queue<DdNode*, std::vector<DdNode*>, Comp_backward> Q;
    std::set<DdNode*> S;
    S.insert(one);
    Q.push(one);
    while ( !Q.empty()){
        DdNode *node = Q.top();
        Q.pop();
        int current_variable = Cudd_NodeReadIndex(node);
        if ( (this->hipar)->count(node)){
            for ( int i = 0; i < (*this->hipar)[node]->size(); i++){
                DdNode *hipar = (*(*this->hipar)[node])[i];
                int hipar_variable = Cudd_NodeReadIndex(hipar);
            	if (! S.count(hipar)){
                	S.insert(hipar);
                	Q.push(hipar);
            	}
            	if ( !(this->backward_message)->count(hipar)){
                	(*this->backward_message)[hipar] = 0;
            	}
            	(*this->backward_message)[hipar] += (*this->backward_message)[node] * (*x)[hipar_variable];
            }
        }
        if ( this->lopar->count(node)){
        	for ( int i = 0; i < (*this->lopar)[node]->size(); i++){
            		DdNode *lopar = (*(*this->lopar)[node])[i];
            		int lopar_variable = Cudd_NodeReadIndex(lopar);
             		if ( !S.count(lopar)){
                	S.insert(lopar);
                	Q.push(lopar);
            	}
            	if ( !(this->backward_message)->count(lopar)){
                	(*this->backward_message)[lopar] = 0;
           	 }
            	(*this->backward_message)[lopar] += (*this->backward_message)[node] * (1 - (*x)[lopar_variable]);
        //   	 std::cout<< "lopar receives" << (*this->backward_message)[node] * (*x)[lopar_variable];
        	}
        }
        if ( Cudd_IsConstant(node)) continue;
        (*grad)[current_variable] -= (*this->forward_message)[node] * ( (*this->backward_message)[Cudd_E(node)] - (*this->backward_message)[Cudd_T(node)]);
    }
}

void BDD::change_weights_in_BDD(std::vector<double> *clause_weights){
    this->clause_weights = clause_weights;
}

std::vector<double> *BDD::grad(std::vector<double> *x){
    int n = x->size();
    std::vector<double> *grad = new std::vector<double>(n);
    std::vector<double> a;
    for( int i=0; i<x->size(); i++) a.push_back( (1 - (*x)[i]) / 2);
    this->forward_pass(&a);
    this->backward_pass(&a, grad);
    return grad;
}

double BDD::fval(std::vector<double> *x){
    DdNode *one = Cudd_ReadOne(this->gbm);
    std::vector<double> a;
    for( int i=0; i<x->size(); i++) a.push_back( (1 - (*x)[i]) / 2);
    this->forward_pass(&a);
    double fval = sum(this->clause_weights) - 2 * (*this->forward_message)[one];
    std::cout<<"fval = " <<fval << std::endl;
    return fval;
}
