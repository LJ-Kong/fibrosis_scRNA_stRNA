for x in `cat ../samples.txt`; do echo Rscript 2.run_bayes_prism.r $x $x.ted.rds; done
