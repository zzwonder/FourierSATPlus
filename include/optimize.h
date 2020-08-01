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
    void update_weights();    
    std::vector<int> *unsat_clauses;   
    std::string optimizer_name;
    int number_of_trials_this_start;
    Optimizer();
    Optimizer(int n, BDD *bdd,std::string name="BFGS");
    double minimize();
    double fval(const column_vector& m);
    column_vector a_grad (const column_vector& m);
    bool solved_flag;
    
    private:
    double fval_for_optimizer(const column_vector& m);
    const column_vector grad_for_optimizer (const column_vector& m);

};


class Optimizer_Portfolio{
    public:
        int ncores;
        int max_trials_per_start;
        int num_of_vars;
        BDD *bdd;
        Optimizer_Portfolio();
        Optimizer_Portfolio(int ncores, int max_trials_per_start, int num_of_vars, BDD *original_bdd);
        void solve();
};
#endif
