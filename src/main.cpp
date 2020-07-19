#include "../include/read_formula.h"
#include "../include/bdd_gradient.h"
#include "../include/main.h"
#include <chrono>

void sign(std::vector<float> *x){
    for ( int i = 0; i<x->size();i++){
        if ((*x)[i] > 0) (*x)[i] = 1;
        else (*x)[i] = -1;
    }
}

int main(){
    Formula formula;
    formula.read_DIMACS("benchmarks/testset/cubic_vc_50_39.cnf");
//    formula.read_DIMACS("benchmarks/testset/rand3_v400_c1680_0.cnf");
 //   formula.read_DIMACS("benchmarks/testset/cubic_vc_300_99.cnf");
    BDD bdd = BDD(&formula);
    std::vector<float> x;
    for ( int i = 0; i < formula.num_of_vars; i++){
        //float xi = (float) rand();
        //x.push_back( 1 - 2 * (xi/(float)RAND_MAX));
        float xi = 0.6;
        x.push_back( 1 - 2 * xi);
    }
    //sign(&x);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<float> *grad;

    for (int j = 0; j<1;j++){
       float fval = bdd.fval(&x);
       std::cout<<"fval = "<<fval <<std::endl;
        grad = bdd.grad(&x);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << duration / 1e6 <<std::endl;
    for (int i = 0; i < formula.num_of_vars; i++){
        std::cout<<(*grad)[i]<<" ";
    }
    std::cout<<std::endl;
    return 0;
}

