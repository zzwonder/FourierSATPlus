#include <dlib/optimization.h>
#include <iostream>
#include "../include/optimize.h"
#include <random>
#include <algorithm>


using namespace std;
using namespace dlib;

BDD *global_bdd;

// ----------------------------------------------------------------------------------------
// Below we create a few functions.  When you get down into main() you will see that
// we can use the optimization algorithms to find the minimums of these functions.
// ----------------------------------------------------------------------------------------

Optimizer::Optimizer(){
}

Optimizer::Optimizer(int n, BDD *bdd){
    srand(time(NULL));
    this->bdd = bdd;
    this->num_of_vars = n;
    this->x = new std::vector<double>;
    for(int i=0; i< n; i++){
        double xi = (double) std::rand();
        this->x->push_back( 1 - 2 * (xi/(double)RAND_MAX));
    }
}


double fval_for_optimizer(const column_vector& m)
{
    std::vector<double> x;
    for (int i=0; i < global_bdd->num_of_vars; i++){
        x.push_back((double)m(i));
    }
    return global_bdd->fval(&x);
}

// This is a helper function used while optimizing the rosen() function.  
const column_vector grad_for_optimizer (const column_vector& m)
{
    std::vector<double> *grad;
    std::vector<double> x;
    for (int i=0; i < global_bdd->num_of_vars; i++){
        x.push_back((double)m(i));
    }
    column_vector res(global_bdd->num_of_vars);
    grad = global_bdd->grad(&x);
    for (int i=0; i < global_bdd->num_of_vars; i++){
        res(i) = (*grad)[i];
    }
    return res;
}

double Optimizer::minimize()
{
    global_bdd = this->bdd;
    std::vector<double> *x = this->x;
    int n = this->num_of_vars;
    while(1){
        try
        {
            for(int i=0; i< n; i++){
                double xi = (double) std::rand();
                (*this->x)[i] =  1 - 2 * (xi/(double)RAND_MAX);
            }

            column_vector starting_point(n);
            for(int i = 0; i < n; i++){
                starting_point(i) = (double) (*x)[i];
            }
        //find_min_box_constrained(lbfgs_search_strategy(10),  
            find_min_box_constrained(bfgs_search_strategy(), 
                                 objective_delta_stop_strategy(1e-9),  
                                fval_for_optimizer, grad_for_optimizer, starting_point, -1.0, 1.0);
            double fval = fval_for_optimizer(starting_point);
            cout << endl << "fval: \n" << fval << endl;
            if ( (fval + global_bdd->sum_of_clause_weights) < 1e-2){ 
                cout<<"solved"<<endl; 
                cout <<"solution: \n" << starting_point << endl;
                break;}
        }
        catch (std::exception& e)
        {
            cout << e.what() << endl;
        }
    }
    return 0;
}


