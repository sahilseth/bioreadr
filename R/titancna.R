

# https://github.com/gavinha/TitanCNA/tree/master/scripts/R_scripts
# numClusters=3
# numCores=4
# for ploidy in [2, 3, 4]:
# for num_clusters in [1, 2, 3]:
# ## run TITAN for each ploidy (2,3,4) and clusters (1 to numClusters)
# echo "Maximum number of clusters: $numClusters";
# for ploidy in $(seq 2 4)
# do
# echo "Running TITAN for $i clusters.";
# outDir=run_ploidy$ploidy
# mkdir $outDir
# for numClust in $(seq 1 $numClusters)
# do
# echo "Running for ploidy=$ploidy";
# Rscript titanCNA_v1.10.1.R --id test --hetFile test.het.txt --cnFile test.cn.txt \
# --numClusters $numClust --numCores $numCores --normal_0 0.5 --ploidy_0 $ploidy \
# --chrs "c(1:22, \"X\")" --estimatePloidy TRUE --outDir $outDir
# done
# echo "Completed job for $numClust clusters."
# done
# 
# ## select optimal solution
# Rscript selectSolution.R --ploidyRun2=run_ploidy2 --ploidyRun3=run_ploidy3 --ploidyRun4=run_ploidy4 --threshold=0.05 --outFile optimalClusters.txt