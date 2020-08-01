#ifndef READ_FORMULA_H
#define READ_FORMULA_H
#include<vector>
#include<string>
#define GE 1
#define EQ 2

class Formula{
    public:
    int num_of_vars;
    std::vector<std::vector<int> *> *clauses;
    std::vector<double> *clause_weights;
    std::vector<double> *clause_weights_original;
    std::vector<char> *clause_type;
    std::vector<int> *klist;
    std::vector<int> *comparators;
    std::vector<std::vector<int> *> *coefs;
    
    Formula();
    void add_clause(std::vector<int> *literals, int k, char ctype, double weight, std::vector<int> *coefsL, int comparator);
    void read_DIMACS(std::string file);
    void print();
    double compute_clause_weight(int n, int k, char ctype);
    std::vector<int> unsat_clauses(std::vector<int> *x); 
};


#endif
