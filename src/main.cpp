#include "../include/read_formula.h"


int main(){
    Formula formula;
    formula.read_DIMACS("benchmarks/cubic_vc_50_0.cnf");
    formula.print(); 
    return 0;
}
