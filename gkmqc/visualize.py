import numpy as np
import logging

# return gkmQC score and curve stats
def gkmQC_stat(eval_file):

    # Parse eval file
    auc_list = []
    f = open(eval_file)
    for line in f.readlines():
        txt, _, num, avg, std = line.split()
        if float(num) >= 5000:
            auc_list.append([int(txt.split('.')[-2][3:]), float(avg), float(std)])
    f.close()
    auc_list.sort(key=lambda x: x[0])
    auc_scores = list(list(zip(*auc_list))[1])
    
    # limit to top 20 peak subsets
    n = len(auc_scores)
    if n > 20:
        auc_scores = auc_scores[:20]
        n = 20
    
    # calculate gkmQC score
    auc_max = max(auc_scores)
    auc_min = min(auc_scores)
    score = sum(auc_scores) / (auc_max - auc_min)  
    logging.info("gkmQC score = %.3f", score)

    # Visualize gkmQC curve
    try:
        import matplotlib.pyplot as plt
    except:
        logging.info("Matplotlib is not installed in the conda environment. Curve PDF file will not be created.")
        return
    
    plt.figure(figsize=(10, 10))
    rank_l, avg_l, std_l = zip(*auc_list[:20])
    plt.errorbar(rank_l, avg_l, yerr=std_l, label=eval_file)
    plt.ylim(0.5, 1.0)
    plt.xlim(0, 21)

    # create curve pdf file
    cpdf_file = eval_file.replace(".eval.out", ".curve.pdf")
    plt.title("%s\ngkmQC score = %.3f" % (eval_file, score))
    plt.xlabel("The rank of peak subsets")
    plt.ylabel("Peak predictability (AUC)")
    plt.savefig(cpdf_file)
    logging.info("Curve PDF file has been created: %s", cpdf_file)
    return

def report(args):
    gkmQC_stat(args.eval_file)