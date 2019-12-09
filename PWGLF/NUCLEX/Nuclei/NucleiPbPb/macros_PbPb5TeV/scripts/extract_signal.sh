#!/bin/bash
if [[ $# < 1 ]]
then
multirun=0
echo "set run for all directory"
else
multirun=$1
echo "set run for $multirun directory"
fi
rootver=$(root-config --version)
if (( ${rootver:0:1} < 6 )); then
  echo "ROOT version 6 is required for this analysis."
else
  root -b -l << EOF
.L src/RooGausExp.cxx+
.L src/RooGausDExp.cxx+
.L src/RooDSCBShape.cxx+
.L src/FitModules.cxx+
.x Signal.cc($multirun)
EOF

fi

