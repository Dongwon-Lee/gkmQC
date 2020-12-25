## gkmQC: gapped k-mer-SVM quality check

gkmQC is a sequence-based quality assessment and refinement of
chromatin accessibility data using gkm-SVM.
It trains a support vector classifier (SVC) using gapped-kmer kernels
(Ghandi et al., 2014; Lee, 2016), and learns sequence features that modulate
gene expressions. We use LIBSVM (Chang & Lin 2011) for implementing SVR.

requires 

* Python >=3
* numpy
* sklearn
* bitarray
* pyfasta

Set conda virtual environment
```bash
$ conda env create -f environment.yml
$ conda activate gkmqc
```

Please compile C library for gkm-kernel
```bash
$ cd ./src
$ make && make install
```

Download Null-seq index (hg38: 5.8GB)
```bash
$ cd ./data
$ wget https://www.dropbox.com/s/wtjylew5ybim29x/gkmqc.idx.hg38.tar.xz?dl=0
$ tar xvfJ gkmqc.idx.hg38.tar.xz
```

Build your own null-seq index
(takes 15 mins with 10-threads)
```bash
$ cd ./data
# Download zipped fa files split by chromosome
$ wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
$ ../bin/gkmqc.py buildidx -i hg38.chromFa.tar.gz -g hg38 -@ [threads]
```

Evaluate your called peaks and check your gkmQC curve.
(takes 1 ~ 2hrs with 10-threads)
```bash
$ ../bin/gkmqc.py evaluate -i test.narrowPeak -g hg38 -n test -@ [threads]
$ cat ./test/test.gkmqc.eval.out

```
You can check the options with -h arg of gkmqc.py
```bash
$ cd ./bin
$ ./gkmqc.py -h
$ ./gkmqc.py buildidx -h # Building null-seq index
$ ./gkmqc.py evaluate -h # run gkm-SVM to evaluate peaks
```

Please cite below papers
* Han SK, Sampson MG, Lee D. Quality assessment and refinement of non-coding regulatory map using a sequence-based predictive model.
* Ghandi M†, Lee D†, Mohammad-Noori M, & Beer MA. Enhanced Regulatory Sequence Prediction Using Gapped k-mer Features. PLoS Comput Biol 10, e1003711 (2014). doi:10.1371/journal.pcbi.1003711 *† Co-first authors*
* Lee D. LS-GKM: A new gkm-SVM for large-scale Datasets. Bioinformatics btw142 (2016). doi:10.1093/bioinformatics/btw142
* Chang C.-C and Lin C.-J. LIBSVM : a library for support vector machines. ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011.

