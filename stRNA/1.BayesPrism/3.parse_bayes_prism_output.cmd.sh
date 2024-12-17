cat ../cell_types.txt | xargs -I {} echo Rscript 3.parse_bayes_prism_output.r \'{}\' all
cat ../cell_groups.txt | xargs -I {} echo Rscript 3.parse_bayes_prism_output.r \'{}\' merged
