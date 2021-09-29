#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <pthread.h>

#include "libgkm.h"

#define CLOG_MAIN
#include "clog.h"

typedef struct _gkmkernel_pthread_data_t {
    gkm_kernel *kernel;
    int nthreads;
    int task_ind;
    double **res;
} gkmkernel_pthread_data_t;

const char *gkm_check_parameter(const gkm_parameter *param)
{
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

    return NULL;
}

/*  
* Main Wrapper Function
*/

static void *pthread_gkmkernel_kernelfunc_batch_all(void *ptr)
{
    int i;
    gkmkernel_pthread_data_t *td = (gkmkernel_pthread_data_t *) ptr;
    int NTHREADS = td->nthreads;

    for (i = 0; i < td->kernel->prob_num/NTHREADS; i++) {
        int a = (i*NTHREADS) + td->task_ind;
        gkmkernel_kernelfunc_batch_all(td->kernel, a, 0, a, td->res[i]); 
        if (i % 10 == 0) {
            clog_info(CLOG(LOGGER_ID), "  Thread %d, i = %d", td->task_ind, i);
        }
    }

    return 0;
}

int main(int argc, char** argv)
{
    int i, j, k;
    gkm_parameter param;
    svm_problem prob;
    gkmOpt _opts, *opts;
    opts = &_opts;

    /* Initialize the logger */
    if (clog_init_fd(LOGGER_ID, 1) != 0) {
        fprintf(stderr, "Logger initialization failed.\n");
        return 1;
    }

    clog_set_fmt(LOGGER_ID, LOGGER_FORMAT);
    clog_set_level(LOGGER_ID, CLOG_INFO);

    if (argc != 4) {
        fprintf(stderr, "Usage: gkmkern [posfile] [negfile] [outfile]\n");
        fprintf(stderr,"Wrong number of arguments [%d].\n", argc-1);
        exit(0);
    }

	int index = 1;
	char *posfile = argv[index++];
	char *negfile = argv[index++];
	char *outfile = argv[index];

    // Set default value ...
    // opts->kernel_type = EST_TRUNC_PW;
    opts->kernel_type = EST_TRUNC; // XXX: This is for comparison with gkmSVM 2.0
    opts->L = 10;
    opts->k = 6;
    opts->d = 3;
    opts->M = 50;
    opts->H = 50;
    opts->gamma = 1.0;
    opts->verbosity = 2;
    opts->nthreads = 4;
    opts->posfile = posfile;
    opts->negfile = negfile;

    // Set up options from gkmOpt
    param.kernel_type = opts->kernel_type;
    param.L = opts->L;
    param.k = opts->k;
    param.d = opts->d;
    param.M = opts->M;
    param.H = opts->H;
    param.gamma = opts->gamma;
    param.nthreads = 1; // no longer used...

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
            fprintf(stderr, "Unknown verbosity: %d\n", opts->verbosity);
            exit(0);
    }

    clog_info(CLOG(LOGGER_ID), "Arguments:");
    clog_info(CLOG(LOGGER_ID), "  posfile = %s", posfile);
    clog_info(CLOG(LOGGER_ID), "  negfile = %s", negfile);
    clog_info(CLOG(LOGGER_ID), "  outfile = %s", outfile);

    clog_info(CLOG(LOGGER_ID), "Parameters:");
    clog_info(CLOG(LOGGER_ID), "  kernel-type = %d", param.kernel_type);
    clog_info(CLOG(LOGGER_ID), "  L = %d", param.L);
    clog_info(CLOG(LOGGER_ID), "  k = %d", param.k);
    clog_info(CLOG(LOGGER_ID), "  d = %d", param.d);
    if (param.kernel_type == EST_TRUNC_RBF || param.kernel_type == EST_TRUNC_PW_RBF) {
    clog_info(CLOG(LOGGER_ID), "  gamma = %g", param.gamma);
    }
    if (param.kernel_type == EST_TRUNC_PW || param.kernel_type == EST_TRUNC_PW_RBF) {
    clog_info(CLOG(LOGGER_ID), "  M = %d", param.M);
    clog_info(CLOG(LOGGER_ID), "  H = %g", param.H);
    }

    const char *error_msg = gkm_check_parameter(&param);
    if(error_msg) {
        clog_error(CLOG(LOGGER_ID), error_msg);
        return 1; // error exit 
    }

    // build-up kernel and check integrity of options
    gkm_kernel *kernel = gkmkernel_init(&param);
    gkmkernel_read_problems(kernel, &prob, opts->posfile, opts->negfile);
    gkmkernel_build_tree(kernel, prob.x, prob.l); //build kmertree using the entire problem set, which will then be used by gkmkernel_kernelfunc_batch_all

    // multithreads by line
    int NTHREADS = opts->nthreads;
    gkmkernel_pthread_data_t *td;
    pthread_t *threads;

    td = (gkmkernel_pthread_data_t *) malloc(sizeof(gkmkernel_pthread_data_t) * ((size_t) NTHREADS));
    threads= (pthread_t *) malloc(sizeof(pthread_t) * ((size_t) NTHREADS));
    
    // initialize data for threads
    for (i=0; i<NTHREADS; i++) {
        td[i].kernel = kernel;
        td[i].task_ind = i;
        td[i].nthreads = NTHREADS;
        td[i].res = (double **) malloc(sizeof(double *) * ((size_t) prob.l/(unsigned int)NTHREADS));
        for(j = 0; j < prob.l/NTHREADS; j++) {
            td[i].res[j] = (double *) malloc(sizeof(double) * 10000);
        }
    }

    int rc[NTHREADS];
    for (i=1; i<NTHREADS; i++) {
        rc[i] = pthread_create(&threads[i], NULL, pthread_gkmkernel_kernelfunc_batch_all, (void *) &td[i]);
        if (rc[i]) {
            clog_error(CLOG(LOGGER_ID), "failed to create thread. pthread_create() returned %d", rc[i]);
        } else {
            clog_trace(CLOG(LOGGER_ID), "thread %d was created.", i);
        }
    }

    for (i=0; i<NTHREADS; i++) {
        if (i == 0) {
            pthread_gkmkernel_kernelfunc_batch_all(&td[i]);
        } else {
            if (rc[i] == 0) {
                //wait thread return
                pthread_join(threads[i], NULL);
            } else {
                //if failed to run thread, execute the function in the main process
                pthread_gkmkernel_kernelfunc_batch_all(&td[i]);
            }
        }
    }

    FILE *fo = fopen(outfile, "w");
    if (fo == NULL) {
        perror ("error occurred while opening a file");
        return 1;
    }

    for(i = 0; i < prob.l/NTHREADS; i++) {
        for(j = 0; j < NTHREADS; j++) {
            for(k = 0; k < (i*NTHREADS) + j; k++) {
				fprintf(fo, "%e\t", td[j].res[i][k]);
            }
            fprintf(fo, "1.0\t\n");
        }
    }

    // free memory
    for(j = 0; j < NTHREADS; j++) {
        for(i = 0; i < prob.l/NTHREADS; i++) {
            free(td[j].res[i]);
        }
        free(td[j].res);
    }
    free(td);

    // remove allocated memory
    for (i = 0; i < prob.l; i++) {
        gkmkernel_delete_object(prob.x[i]);
    }
    free(prob.y);
    free(prob.x);

    gkmkernel_destroy(kernel);
    clog_free(LOGGER_ID);
    return 0;
}
