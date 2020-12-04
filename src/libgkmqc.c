#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "libgkm.h"
#include "libsvm.h"

#define CLOG_MAIN
#include "clog.h"

int count_sequences(const char *filename)
{
    FILE *fp = fopen(filename,"r");
    int nseqs = 0;

    if(fp == NULL) {
        clog_error(CLOG(LOGGER_ID), "can't open file");
        exit(1);
    }

    //count the number of sequences for memory allocation
    while(readline(fp)!=NULL) {
        if (line[0] == '>') {
            ++nseqs;
        }
    }
    fclose(fp);
    
    return nseqs;
}

void read_fasta_file(const char *filename, int offset, int label)
{
    FILE *fp = fopen(filename, "r");

    if(fp == NULL) {
        clog_error(CLOG(LOGGER_ID), "can't open file");
        exit(1);
    }

    int iseq = -1;
    char seq[MAX_SEQ_LENGTH];
    char sid[MAX_SEQ_LENGTH];
    int seqlen = 0;
    sid[0] = '\0';
    while (readline(fp)) {
        if (iseq >= prob.l) {
            clog_error(CLOG(LOGGER_ID), "error occured while reading sequence file (%d >= %d).\n", iseq, prob.l);
            exit(1);
        }

        if (line[0] == '>') {
            if (((iseq % 1000) == 0)) {
                clog_info(CLOG(LOGGER_ID), "reading... %d", iseq);
            }

            if (iseq >= 0) {
                prob.y[offset + iseq] = label;
                prob.x[offset + iseq].d = gkmkernel_new_object(seq, sid, offset + iseq);
            }
            ++iseq;

            seq[0] = '\0'; //reset sequence
            seqlen = 0;
            char *ptr = strtok(line," \t\r\n");
            if (strlen(ptr) >= MAX_SEQ_LENGTH) {
                clog_error(CLOG(LOGGER_ID), "maximum sequence id length is %d.\n", MAX_SEQ_LENGTH-1);
                exit(1);
            }
            strcpy(sid, ptr+1);
        } else {
            if (seqlen < MAX_SEQ_LENGTH-1) {
                if ((((size_t) seqlen) + strlen(line)) >= MAX_SEQ_LENGTH) {
                    clog_warn(CLOG(LOGGER_ID), "maximum sequence length allowed is %d. The first %d nucleotides of %s will only be used (Note: You can increase the MAX_SEQ_LENGTH parameter in libsvm_gkm.h and recompile).", MAX_SEQ_LENGTH-1, MAX_SEQ_LENGTH-1, sid);
                    int remaining_len = MAX_SEQ_LENGTH - seqlen - 1;
                    line[remaining_len] = '\0';
                }
                strcat(seq, line);
                seqlen += (int) strlen(line);
            }
        }
    }

    //last one
    prob.y[offset + iseq] = label;
    prob.x[offset + iseq].d = gkmkernel_new_object(seq, sid, offset + iseq);

    clog_info(CLOG(LOGGER_ID), "reading... done");

    fclose(fp);
}

void read_problem(struct svm_problem *prob, const char *posfile, const char *negfile)
{
    int n1 = count_sequences(posfile);
    int n2 = count_sequences(negfile);

    prob->l = n1+n2;

    prob->y = (double *) malloc (sizeof(double) * ((size_t) prob->l));
    prob->x = (union svm_data *) malloc(sizeof(union svm_data) * ((size_t) prob->l));

    clog_info(CLOG(LOGGER_ID), "reading %d sequences from %s", n1, posfile);
    read_fasta_file(posfile, 0, 1);

    clog_info(CLOG(LOGGER_ID), "reading %d sequences from %s", n2, negfile);
    read_fasta_file(negfile, n1, -1);
}

/*  
* Main Wrapper Function
*/

typedef struct _gkmOpt {
    int kernel_type;
    int L;
    int k;
    int d;
    uint8_t M;
    double H;
    double gamma;
    char *posfile;
    char *negfile;
    int nthreads;
    int verbosity;
} gkmOpt;

int gkmKernelPyWrapper(gkmOpt *opts, double **kmat, int *nseqs)
{
    int i, c;
    char* error_msg;
    struct svm_problem prob;
    struct svm_parameter param;
    struct svm_model *model;

    // Set up options from gkmOpt
    param.svm_type = C_SVC;
    param.kernel_type = opts->kernel_type;
    param.L = opts->L;
    param.k = opts->k;
    param.d = opts->d;
    param.M = opts->M;
    param.H = opts->H;

    // belows are related to training step
    param.gamma = opts->gamma;
    param.cache_size = 100; // dummy
    param.C = 1; // dummy
    param.eps = 1e-3; // dummy
    param.shrinking = 0; // dummy
    param.nr_weight = 0; // dummy
    param.weight_label = (int *) malloc(sizeof(int)*1); // ??
    param.weight = (double *) malloc(sizeof(double)*1); // ??

    switch(opts->verbosity) 
    {
        case 0:
            clog_set_level(LOGGER_ID, CLOG_ERROR);
            break;
        case 1:
            clog_set_level(LOGGER_ID, CLOG_WARN);
            break;
        case 2:
            clog_set_level(LOGGER_ID, CLOG_INFO);
            break;
        case 3:
            clog_set_level(LOGGER_ID, CLOG_DEBUG);
            break;
        case 4:
            clog_set_level(LOGGER_ID, CLOG_TRACE);
            break;
        default:
            fprintf(stderr, "Unknown verbosity: %d\n", verbosity);
            print_usage_and_exit();
    }

    gkmkernel_set_num_threads(opts->nthreads);

    // build-up kernel and check integrity of options
    gkmkernel_init(&param);
    read_problem(&prob, opts->posfile, opts->negfile);
    error_msg = svm_check_parameter(&param);

    if(error_msg) {
        clog_error(CLOG(LOGGER_ID), error_msg);
        return 1; // error exit 

    // build kernel
    Kernel::Kernel(int l, svm_data const * x_, const svm_parameter& param):kernel_type(param.kernel_type)
    // spec svm_data <- prob.x, param <- param
    // TODO: Kernel::k_function
    //model = svm_train(&prob, &param);

    // remove allocated memory
    for (i = 0; i < prob.l; i++) {
        gkmkernel_delete_object(prob.x[i].d);
    }
    svm_destroy_param(&param);
    free(prob.y);
    free(prob.x);

    return 0;
}