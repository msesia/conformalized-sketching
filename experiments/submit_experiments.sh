#!/bin/bash


# Parameters
CONF=102

#SKETCH_LIST=("cms" "cms-cu")
SKETCH_LIST=("cms" "cms-cu")

if [[ $CONF == 1 ]]; then
  DATA_LIST=("dp-1" "dp-100" "pyp-1-0.5" "pyp-100-0.5")
  N_LIST=(10000 20000 50000 100000 200000)
  CONFIDENCE_LIST=(50 95)
  METHOD_LIST=("classical" "bootstrap" "conformal-adaptive" "conformal-constant")
  W_LIST=1000

elif [[ $CONF == 2 ]]; then
  DATA_LIST=("pyp-5000-0.0" "pyp-5000-0.1" "pyp-5000-0.2" "pyp-5000-0.3" "pyp-5000-0.4" "pyp-5000-0.5" "pyp-5000-0.6" "pyp-5000-0.7" "pyp-5000-0.8" "pyp-5000-0.9" "pyp-5000-0.95")
  N_LIST=(100000) # 5000 10000 50000)
  CONFIDENCE_LIST=(95)
  METHOD_LIST=("classical" "bayesian-dp-?" "bootstrap" "conformal-constant" "conformal-adaptive")
  W_LIST=1000

elif [[ $CONF == 3 ]]; then
  DATA_LIST=("zipf-1.05" "zipf-1.1" "zipf-1.2" "zipf-1.3" "zipf-1.4" "zipf-1.5" "zipf-1.6" "zipf-1.7" "zipf-1.8" "zipf-1.9" "zipf-2.0")
  N_LIST=(1000000)
  CONFIDENCE_LIST=(95)
  #METHOD_LIST=("classical" "bootstrap" "conformal-adaptive" "conformal-constant")
  METHOD_LIST=("classical" "bayesian-dp-?" "bootstrap" "conformal-constant" "conformal-adaptive")
  W_LIST=(1000)

elif [[ $CONF == 101 ]]; then
  DATA_LIST=("covid")
  N_LIST=(1000000)
  CONFIDENCE_LIST=(50 95)
  METHOD_LIST=("classical" "bayesian-dp-?" "bootstrap" "conformal-constant" "conformal-adaptive")
  W_LIST=(1000 2000 5000 10000 20000 50000)
  SKETCH_LIST=("cms" "cms-cu")

elif [[ $CONF == 102 ]]; then
  DATA_LIST=("words")
  N_LIST=(1000000)
  CONFIDENCE_LIST=(50 95)
  METHOD_LIST=("classical" "bayesian-dp-?" "bootstrap" "conformal-constant" "conformal-adaptive")
  W_LIST=(1000 2000 5000 10000 20000 50000)
  SKETCH_LIST=("cms" "cms-cu")

fi

D_LIST=3
SEED_LIST=$(seq 1 10)
POSTERIOR_LIST=("mcmc") #approximate")


# Slurm parameters
MEMO=20G                             # Memory required (5 GB) # Note: 20 GB are needed for covid data
TIME=00-00:10:00                    # Time required (20m) # 2 hours for Bayesian model on large dataset
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

OUT_DIR="results"


# Loop over configurations and chromosomes
for DATA in "${DATA_LIST[@]}"; do
  OUT_DIR_DATA=$OUT_DIR"/"$DATA
  mkdir -p $OUT_DIR_DATA
  mkdir -p $OUT_DIR_DATA"/detailed"
  mkdir -p $OUT_DIR_DATA"/marginal"
  mkdir -p $OUT_DIR_DATA"/conditional"

  for SKETCH in "${SKETCH_LIST[@]}"; do
    for CONFIDENCE in "${CONFIDENCE_LIST[@]}"; do
      for POSTERIOR in "${POSTERIOR_LIST[@]}"; do
        for SEED in $SEED_LIST; do
          for D in "${D_LIST[@]}"; do
            for W in "${W_LIST[@]}"; do
              for N in "${N_LIST[@]}"; do
                for METHOD in "${METHOD_LIST[@]}"; do

                  if [[ $METHOD != conformal* ]]; then
                    if [[ $CONFIDENCE = 50 ]]; then
                      continue
                    fi
                  fi

                  if [[ $METHOD = conformal* ]]; then
                    N_BINS_LIST=(5)
                    N_TRACK_LIST=(5000)
                    UNIQUE_LIST=(0 1)
                  else
                    N_BINS_LIST=("NA")
                    N_TRACK_LIST=("NA")
                    UNIQUE_LIST=("NA")
                  fi

                  for UNIQUE in "${UNIQUE_LIST[@]}"; do
                    for N_BINS in "${N_BINS_LIST[@]}"; do
                      for N_TRACK in "${N_TRACK_LIST[@]}"; do

                        if [[ $METHOD = conformal* ]]; then
                          METHOD_NAME=$METHOD"_unique"$UNIQUE"_bins"$N_BINS"_track"$N_TRACK
                        else
                          METHOD_NAME=$METHOD
                        fi

                        JOBN=$SKETCH"_"$DATA"_d"$D"_w"$W"_n"$N"_s"$SEED"_"$POSTERIOR"_"$METHOD_NAME"_"$CONFIDENCE
                        OUT_FILE=$OUT_DIR_DATA"/detailed/exp1_"$JOBN".txt"
                        COMPLETE=0
                        if [[ -f $OUT_FILE ]]; then
                          COMPLETE=1
                        fi

                        if [[ $COMPLETE -eq 0 ]]; then
                          # Script to be run
                          SCRIPT="experiment_1.sh $SKETCH $DATA $D $W $N $METHOD $UNIQUE $N_BINS $N_TRACK $SEED $POSTERIOR $CONFIDENCE"
                          # Define job name for this chromosome
                          OUTF=$LOGS"/"$JOBN".out"
                          ERRF=$LOGS"/"$JOBN".err"
                          # Assemble slurm order for this job
                          ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
                          # Print order
                          echo $ORD
                          # Submit order
                          #$ORD
                          # Run command now
                          #./$SCRIPT
                        fi
                      done
                    done
                  done
                done
              done
            done
          done
        done
      done
    done
  done
done
