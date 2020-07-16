#include "../include/read_formula.h"
#include "../include/bdd_gradient.h"

int main(){
    Formula formula;
//   formula.read_DIMACS("benchmarks/testset/cubic_vc_50_0.cnf");
    formula.read_DIMACS("benchmarks/testset/8_0.cnf");
    formula.print(); 
    BDD bdd = BDD(&formula);
    
    return 0;
}
