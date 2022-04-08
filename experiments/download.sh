DATA=$1

mkdir -p results_hpc

if [[ $DATA == "all" ]]; then
  rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/cms/code/experiments/cms-conformal/results/ results_hpc/ --exclude "*detailed*"
elif [[ $DATA == "detailed" ]]; then
  mkdir -p results_hpc/detailed/
#  rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/cms/code/experiments/cms-conformal/results/covid/detailed/exp1_cms_covid_d3_w50000_n1000000_s1_mcmc_*_95.txt results_hpc/detailed/
#  rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/cms/code/experiments/cms-conformal/results/words/detailed/exp1_cms_words_d3_w50000_n1000000_s1_mcmc_*_95.txt results_hpc/detailed/
#  rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/cms/code/experiments/cms-conformal/results/covid/detailed/exp1_cms_covid_d3_w5000_n1000000_s1_mcmc_*_95.txt results_hpc/detailed/
#  rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/cms/code/experiments/cms-conformal/results/words/detailed/exp1_cms_words_d3_w5000_n1000000_s1_mcmc_*_95.txt results_hpc/detailed/
else
  rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/cms/code/experiments/cms-conformal/results/$DATA/marginal/* results_hpc/$DATA/marginal/
  rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/cms/code/experiments/cms-conformal/results/$DATA/conditional/* results_hpc/$DATA/conditional/
  mkdir -p results_hpc/$DATA
fi

#rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/cms/code/experiments/cms-conformal/results_conditional/* results_hpc/results_conditional/
