#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "libgkm.h"

#define CLOG_MAIN
#include "clog.h"

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

int main(int argc, char** argv)
{
    int i, j;
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
    param.nthreads = opts->nthreads;

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

    FILE *fo = fopen(outfile, "w");
    if (fo == NULL) {
        perror ("error occurred while opening a file");
        return 1;
    }

    double res[10000];
    for(i = 0; i < prob.l; i++) {
        gkmkernel_kernelfunc_batch_all(kernel, i, 0, prob.l, res);
        if (i % 10 == 0) {
            clog_info(CLOG(LOGGER_ID), "  i = %d", i);
        }
		for(j = 0; j <= i; j++) {
			if(j < i) {
				fprintf(fo, "%e\t", res[j]);
			}
			else if (i == j) {
				fprintf(fo, "1.0\t");
			}
		}
		fprintf(fo, "\n");
    }
	fclose(fo);

    // remove allocated memory
    for (i = 0; i < prob.l; i++) {
        gkmkernel_delete_object(prob.x[i]);
    }
    free(prob.y);
    free(prob.x);

    gkmkernel_destroy(kernel);

    return 0;
}
