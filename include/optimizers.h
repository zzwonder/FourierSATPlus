#ifndef OPTIMIZER_H
#define OPTIMIZER_H
#include<vector>
#include<string>
class Formula{
    int num_of_vars;
    vector<vector<int> *> *clauses;
    vector<float> *clause_weights;
    vector<char> *clause_type;
    vector<int> *klist;
    vector<int> *comparators;
    vector<vector<int> *> *coefs;
    
    Formula();
    void add_clause(vector<int> *literals, int k, char ctype, float weight, vector<int> *coefs = NULL, int comparator = NULL);
    void read_DIMACS(string file);
     
};


#endif
