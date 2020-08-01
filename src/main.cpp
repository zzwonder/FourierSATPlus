#include "../include/read_formula.h"
#include "../include/bdd_gradient.h"
#include "../include/optimize.h"
#include "../include/main.h"
#include <chrono>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <stdint.h>
#include <errno.h>
#include <limits.h>

struct Parameter {
    int solver;
    int k;
    double eps;
    int max_iter;
    int is_unspeficied_wcnf;
    int verbose;
    int n_trial;
    FILE *fin;
    char *fin_name;
    double beta;
    int adapt;
    int ncores;
    int max_trial_per_start;
};

enum {MAXCUT=0, MAXSAT};

void print_usage(char* prog_name, Parameter *param)
{
    printf( "%s [OPTIONS] INPUT: Mixing method for SDP\n", prog_name); 
    printf( "OPTIONS:\n");
    printf( "\t-s SOLVER: type of solver\n");
    printf( "\t           \"-s maxcut\" for maximum cut\n");
    printf( "\t           \"-s maxsat\" for maximum SAT (default)\n");
    printf( "\t-k RANK: rank of solution (default auto)\n");
    printf( "\t         use \"-k /2\" to divide the rank by 2\n");
    printf( "\t-e EPS: stopping threshold (default %.4e)\n", param->eps);
    printf( "\t-t MAX_ITER_PER_START: maximum trial per restart (default %d)\n", param->max_trial_per_start);
    printf( "\t             use \"-t max\" for INT_MAX\n");
    printf( "\t-r N_TRIAL: number of trial in evaluation (default %d)\n", param->n_trial);
    printf( "\t-u: use unspeficied wcnf format\n");
    printf( "\t-v: verbose\n");
    printf( "\t-b beta: momentum scalar (default %f)\n", param->beta);
    printf( "\t-n ncores: number of cores (default %f)\n", param->ncores);
}

void get_parameter(int argc, char **argv, Parameter *param)
{
    Parameter _param = {
        MAXCUT, // solver
         -1, // k
	 1e-3, //eps
         1000, // maxiter
         0, // is unspecified cnf
         0, // verbose
        1, // n trials
        NULL, //fin
        NULL, // fin name
        0, //beta
        0, // adapt
        1, // ncores
        5 //max trial per start
    };

    if(argc <= 1){
        print_usage(argv[0], &_param);
        exit(0);
    }

    char **p = argv+1;
    int i;
    for(i=1; i<argc; i++, p++){
        if(!strcmp(*p, "-s")){
            if(i+1 >= argc) break;
            if(!strcmp(p[1], "maxcut")){
                _param.solver = MAXCUT;
            }else if(!strcmp(p[1], "maxsat")){
                _param.solver = MAXSAT;
            }else {
                int ret = sscanf(p[1], "%d", &_param.solver);
                if(ret != 1 || !(_param.solver >=0 && _param.solver <= 1)) break;
            }
            i++, p++;
        }else if(!strcmp(*p, "-k")){
            if(i+1  >= argc) break;
            int ret = sscanf(p[1], "/%d", &_param.k);
            if(ret==1){
                _param.k *= -1;
            }else{
                ret = sscanf(p[1], "%d", &_param.k);
                if(ret != 1 || _param.k <= 0) break;
            }
            i++, p++;
        }else if(!strcmp(*p, "-e")){
            if(i+1 >= argc) break; 
            int ret = sscanf(p[1], "%lf", &_param.eps);
            if(ret != 1) break; 
            i++, p++;
        }else if(!strcmp(*p, "-t")){
            if(i+1 >= argc) break; 
            if(!strcmp(p[1], "max")){
                _param.max_trial_per_start = INT_MAX;
            }else{
                int ret = sscanf(p[1], "%d", &_param.max_trial_per_start);
                if(ret != 1) break;
            }
            i++, p++;
        }else if(!strcmp(*p, "-n")){
            if(i+1 >= argc) break;
             int ret = sscanf(p[1], "%d", &_param.ncores);
            if(ret != 1) break;
            i++, p++;
        }
        else if(!strcmp(*p, "-r")){
            if(i+1 >= argc) break;
            int ret = sscanf(p[1], "%d", &_param.n_trial);
            if(ret != 1) break; 
            if(_param.n_trial < 1)
                _param.n_trial = 1;
            i++, p++;
        }else if(!strcmp(*p, "-u")){
            _param.is_unspeficied_wcnf = 1;
        }else if(!strcmp(*p, "-v")){
            _param.verbose = 1;
        }else if(!strcmp(*p,"-b")){
            if(i+1 >= argc) break;
            int ret = sscanf(p[1], "%f", &_param.beta);
            if(ret != 1) break;
            i++; p++;
        }else if(!strcmp(*p,"-a")){
            if(i+1 >= argc) break;
            int ret = sscanf(p[1], "%d", &_param.adapt);
            if(ret != 1) break;
            i++; p++;
        }  
        else if(i+1 == argc){
            _param.fin = fopen(*p, "r");
            if(!_param.fin){
                fprintf(stderr, "%s\n", strerror(errno));
                exit(1);
            }
            _param.fin_name = strdup(*p);
        }else{
            printf("Error: no such parameter\n");
            break;
        }
    }
    if(i != argc || !_param.fin){
        print_usage(argv[0], &_param);
        exit(0);
    }
    *param = _param;
}


int main(int argc, char **argv){
    Parameter param;
    get_parameter(argc, argv, &param);
    srand48(0);
    Formula formula;
    formula.read_DIMACS(param.fin_name);
    BDD bdd = BDD(&formula);
    auto t1 = std::chrono::high_resolution_clock::now();
    Optimizer_Portfolio op = Optimizer_Portfolio(param.ncores, param.max_trial_per_start, formula.num_of_vars, &bdd);
    op.solve();    
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    std::cout << duration / 1e6<<"s"<<std::endl;
    std::cout<<std::endl;
    return 0;
}

