//
// Created by yirao zhang on 2023-06-26.
//

#ifndef COMPOSITEILM_GENERATE_INDIVIDUALS_H
#define COMPOSITEILM_GENERATE_INDIVIDUALS_H
#include "ILM.h"
void generate(Individual* pop);
void assign_int_inf(Individual* pop, int id);
void generate_clustered(Individual pop[], int K);
#endif //COMPOSITEILM_GENERATE_INDIVIDUALS_H
