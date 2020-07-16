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
    std::vector<float> *clause_weights;
    std::vector<char> *clause_type;
    std::vector<int> *klist;
    std::vector<int> *comparators;
    std::vector<std::vector<int> *> *coefs;
    
    Formula();
    void add_clause(std::vector<int> *literals, int k, char ctype, float weight, std::vector<int> *coefsL, int comparator);
    void read_DIMACS(std::string file);
    void print();
    float compute_clause_weight(int n, int k, char ctype); 
};


#endif
