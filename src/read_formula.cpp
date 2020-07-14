#include "../include/read_formula.h"
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<stdlib.h>
#include <iterator>

Formula::Formula(){
    this->num_of_vars = 0;
    this->clauses = new std::vector<std::vector<int> *>;
    this->coefs = new std::vector<std::vector<int> *>;
    this->clause_weights = new std::vector<float>;
    this->clause_type = new std::vector<char>;
    this->klist = new std::vector<int>;
    this->comparators = new std::vector<int>;
}

void Formula::add_clause(std::vector<int>* literals, int k, char ctype, float weight, std::vector<int> *coefs = NULL, int comparator =0){
   this->clauses->push_back(literals);
   this->klist->push_back(k);
   this->clause_type->push_back(ctype);
   this->clause_weights->push_back(weight);
   this->coefs->push_back(coefs);
   this->comparators->push_back(comparator);
}

void Formula::print(){
    for (int i=0;i<this->clauses->size();i++){
        std::cout<<(*this->clause_type)[i]<<" ";
        for (int j=0;j<(*this->clauses)[i]->size();j++){
            std::cout<<(*(*this->clauses)[i])[j]<<" ";
        }
        std::cout<<std::endl;
    }
}

void Formula::read_DIMACS(std::string file){
    std::ifstream input_file;
    input_file.open(file.c_str());
    if(!input_file.is_open()){
        std::perror("error open");
        exit(1);
    }
    std::string line;
    while(std::getline(input_file,line)){
        std::stringstream ss(line);
    	std::istream_iterator<std::string> begin(ss);
    	std::istream_iterator<std::string> end;
    	std::vector<std::string> split(begin, end);
        std::vector<int> *literals = new std::vector<int>;
        int weight;
        char ctype;
        if ((split.size() >= 2) && (split[1] == "#variable=") ){
            this->num_of_vars = std::stoi(split[2]);
        }
        if ( (split.size()==0) || (split[0]=="c") or (split[0])=="*") continue;
        if ( split[0] == "p"){
            this->num_of_vars = std::stoi(split[2]);   
        }
        else if (split[0] == "g"){
            int k = std::stoi(split[1]);
            weight = this->compute_clause_weight(this->num_of_vars,k,'g');
            if (k>=0){
                for(int i=1; i<=this->num_of_vars;i++){
                    literals->push_back(i);}
                this->add_clause(literals,k,'c',weight);
	    }
            else{
                for(int i=1; i<=this->num_of_vars;i++){
                    literals->push_back(-i);}
                this->add_clause(literals, this->num_of_vars + k, 'c', weight);
            }
        }
        else if (split[0] == "x"){
            for(int i=1; i<split.size()-1;i++){
                literals->push_back(std::stoi(split[i]));
            }
            weight = this->compute_clause_weight(literals->size(),1,'x');
            this->add_clause(literals, 1 , 'x', weight);
        }
        else if (split[0] == "d" ){
            int k = std::stoi(split[1]);
            if (k>=0){
            	for(int i=2; i<split.size()-1;i++){
                	literals->push_back(std::stoi(split[i]));
            	}
                weight = this->compute_clause_weight(literals->size(),k,'d');
                this->add_clause(literals, k , 'c', weight);
            }
             else{
                for(int i=2; i<split.size()-1;i++){
                        literals->push_back(-std::stoi(split[i]));
                }
                weight = this->compute_clause_weight(literals->size(),k,'d');
                this->add_clause(literals, literals->size() + k , 'c', weight);
            }

            
        }
        else if (split[0] == "n" ){
           for(int i=1; i<split.size()-1;i++){
                literals->push_back(std::stoi(split[i]));
            }
            weight = this->compute_clause_weight(literals->size(),1,'n');
            this->add_clause(literals, 1 , 'n', weight);            
        }
        else if (split[1][0] == 'x' ){  // pseudo-Boolean constraints
            std::vector<int> *coefs = new std::vector<int>;
            for(int i=0; i < (split.size() - 3) / 2;i++){
                coefs->push_back(std::stoi(split[i*2]));
                literals->push_back(std::stoi(split[i*2+1]));
            }
            std::string comparator = split[split.size()-3];
            int k = std::stoi(split[split.size() -2];
            canonicalize(literals, &k, &comparator); // comparator: >= (1), =(2)  
            weight = this->compute_clause_weight(literals->size(),k,'p');
            this->add_clause(literals, k , 'p', weight);
        }
        else{
            for(int i=0; i<split.size()-1;i++){
                literals->push_back(stoi(split[i]));
    	    }
            weight = this->compute_clause_weight(literals->size(),1,'c');
            this->add_clause(literals, 1 , 'c', weight);
        }
    }
}


float Formula::compute_clause_weight(int n, int k, char ctype){
    return n;
}
