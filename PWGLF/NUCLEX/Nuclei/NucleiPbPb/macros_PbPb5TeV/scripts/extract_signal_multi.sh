#!/bin/bash
# Will use all cores at a time

cat scripts/parallel.txt | xargs -P6 -L 1 scripts/extract_signal.sh
rm ../results/signal.root
hadd ../results/signal.root ../results/signal_*.root