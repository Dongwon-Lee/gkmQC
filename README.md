## gkmQC: gapped k-mer-SVM quality check


gkmQC is a sequence-based quality assessment and refinement of
chromatin accessibility data using gkm-SVM.

It trains a support vector regression (SVR) using gapped-kmer kernels
(Ghandi et al., 2014; Lee, 2016), and learns sequence features that modulate
gene expressions. We use LIBSVM (Chang & Lin 2011) for implementing SVR.

requires sklearn, numpy working on Python 3.x

* Ghandi M†, Lee D†, Mohammad-Noori M, & Beer MA. Enhanced Regulatory Sequence Prediction Using Gapped k-mer Features. PLoS Comput Biol 10, e1003711 (2014). doi:10.1371/journal.pcbi.1003711 *† Co-first authors*

* Lee D. LS-GKM: A new gkm-SVM for large-scale Datasets. Bioinformatics btw142 (2016). doi:10.1093/bioinformatics/btw142

* Chang C.-C and Lin C.-J. LIBSVM : a library for support vector machines. ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011.


TBA...  

Please email Seong Kyu Han (seongkyu.han AT childrens DOT harvard DOT edu) if you have any questions.
