#include "libsvm.h"

const char *svm_check_parameter(const svm_parameter *param)
{
    // svm_type

    int svm_type = param->svm_type;
    if(svm_type != C_SVC &&
       svm_type != NU_SVC &&
       svm_type != ONE_CLASS &&
       svm_type != EPSILON_SVR &&
       svm_type != NU_SVR)
        return "unknown svm type";
    
    // kernel_type
    
    int kernel_type = param->kernel_type;
    if(kernel_type != GKM &&
       kernel_type != EST_FULL &&
       kernel_type != EST_TRUNC &&
       kernel_type != EST_TRUNC_RBF &&
       kernel_type != EST_TRUNC_PW &&
       kernel_type != EST_TRUNC_PW_RBF)
        return "unknown kernel type";

    if(param->L < 2)
        return "L < 2";

    if(param->L > 12)
        return "L > 12";

    if(param->k > param->L)
        return "k > L";

    if(param->d > (param->L - param->k))
        return "d > L - k";

    // cache_size,eps,C,nu,p,shrinking
    if(param->cache_size <= 0)
        return "cache_size <= 0";

    if(param->eps <= 0)
        return "eps <= 0";

    if(svm_type == C_SVC ||
       svm_type == EPSILON_SVR ||
       svm_type == NU_SVR)
        if(param->C <= 0)
            return "C <= 0";

    if(svm_type == NU_SVC ||
       svm_type == ONE_CLASS ||
       svm_type == NU_SVR)
        if(param->nu <= 0 || param->nu > 1)
            return "nu <= 0 or nu > 1";

    if(svm_type == EPSILON_SVR)
        if(param->p < 0)
            return "p < 0";

    if(param->shrinking != 0 &&
       param->shrinking != 1)
        return "shrinking != 0 and shrinking != 1";

    if(param->probability != 0 &&
       param->probability != 1)
        return "probability != 0 and probability != 1";

    if(param->probability == 1 &&
       svm_type == ONE_CLASS)
        return "one-class SVM probability output not supported yet";

    return NULL;
}


void svm_destroy_param(svm_parameter* param)
{
    free(param->weight_label);
    free(param->weight);
}