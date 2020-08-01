#include <dlib/optimization.h>
#include <iostream>
#include "../include/optimize.h"
#include <random>
#include <algorithm>
#define CHANGE_WEIGHTS
#define MAX_TRIALS_BEFORE_RESTART 15
#include "omp.h"

#define CORES 4

using namespace std;
using namespace dlib;

Optimizer::Optimizer(){
}

Optimizer::Optimizer(int n, BDD *bdd, std::string name){
    srand(time(NULL));
    this->bdd = bdd;
    this->num_of_vars = n;
    this->solved_flag = 0;
    this->x = new std::vector<double>;
    this->unsat_clauses = new std::vector<int>;
    for(int i=0; i< n; i++){
        double xi = (double) std::rand();
        this->x->push_back( 1 - 2 * (xi/(double)RAND_MAX));
    }
    this->optimizer_name = name;
    this->number_of_trials_this_start = 0;
}


double Optimizer::fval_for_optimizer(const column_vector& m)
{
    std::vector<double> x;
    for (int i=0; i < this->bdd->num_of_vars; i++){
        x.push_back((double)m(i));
    }
    return this->bdd->fval(&x);
}

const column_vector Optimizer::grad_for_optimizer (const column_vector& m)
{
    std::vector<double> *grad;
    std::vector<double> x;
    for (int i=0; i < this->bdd->num_of_vars; i++){
        x.push_back((double)m(i));
    }
    column_vector res(this->bdd->num_of_vars);
    grad = this->bdd->grad(&x);
    for (int i=0; i < this->bdd->num_of_vars; i++){
        res(i) = (*grad)[i];
    }
    return res;
}

void Optimizer::update_weights(){
    for ( int i = 0; i < this->unsat_clauses->size(); i++){
        (*this->bdd->clause_weights) [ (*this->unsat_clauses)[i] ] *= 2;
    }
}

static std::vector<double> *rounding( column_vector *m, int n){
    std::vector<double> *res = new std::vector<double>;
    for ( int i = 0; i < n; i++){
        if ( (*m)(i) > 0 )
            res->push_back(1);
        else res->push_back(-1);
    }
    return res;  
}

static void random_restart(std::vector<double> *x){
    int n = x->size();
    for(int i=0; i< n; i++){
        double xi = (double) std::rand();
        (*x)[i] =  1 - 2 * (xi/(double)RAND_MAX);
    }

}

double Optimizer::minimize()
{
    std::vector<double> *x = this->x;
    int n = this->num_of_vars;
    random_restart(this->x);
    this->number_of_trials_this_start ++;
    if (  number_of_trials_this_start  > MAX_TRIALS_BEFORE_RESTART ){
         std::cout<<"random restart"<<std::endl;
         random_restart(this->x);
         this->bdd->restart_weights();
         this->number_of_trials_this_start = 0;
    }
    column_vector starting_point(n);
    for(int i = 0; i < n; i++){
        starting_point(i) = (double) (*x)[i];
    }
    find_min_box_constrained(cg_search_strategy(), 
                            objective_delta_stop_strategy(1e-9),  
                           [this](const column_vector& a){ return this->fval_for_optimizer(a);}, 
                           [this](const column_vector& a){ return this->grad_for_optimizer(a);},
                           starting_point, -1.0, 1.0);
    std::vector<double> *roundedx = rounding(&starting_point,n);
    int num_unsat_clause = this->bdd->verify_solution(roundedx,this->unsat_clauses);
    double fval = fval_for_optimizer(starting_point);
    #ifdef CHANGE_WEIGHTS
        this->update_weights();
    #endif
    std::cout<<"num of unsat: "<<num_unsat_clause<<std::endl;
    if (num_unsat_clause==0){ 
        cout<<"solved"<<endl; 
        cout <<"solution: \n" << starting_point << endl;
//        cout <<"number of trials: \n" << number_of_total_trials << endl;
        this->solved_flag = 1;
    }
    return 0;
}

Optimizer_Portfolio::Optimizer_Portfolio(){
}


Optimizer_Portfolio::Optimizer_Portfolio(int ncores, int max_trials_per_start, int num_of_vars, BDD *original_bdd){
    this->ncores = ncores;
    this->max_trials_per_start = max_trials_per_start;
    this->num_of_vars = num_of_vars;
    this->bdd = original_bdd;
}

void Optimizer_Portfolio::solve(){
    BDD bdd_group[this->ncores];
    Optimizer optimizer_group[this->ncores];
    for ( int i = 0; i < this->ncores; i++){
        bdd_group[i] = BDD(this->bdd);
        optimizer_group[i] = Optimizer(this->num_of_vars, &bdd_group[i], "BFGS");
    }
    int num_trials = 0;
    bool solved_flag = 0;
    cout<<"ncores "<<this->ncores<<endl;
    while(!solved_flag){
        #pragma omp parallel for num_threads(this->ncores)
        for ( int i = 0; i < this->ncores; i++){
            optimizer_group[i].minimize();
            int tid = omp_get_thread_num();
            std::cout<<"thread id "<<tid<<std::endl;
            solved_flag |= optimizer_group[i].solved_flag;
        }
        num_trials += this->ncores;
    }
    cout<<"number of trials: "<<num_trials<<endl;
}
