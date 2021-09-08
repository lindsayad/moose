#!/bin/bash

export METHOD=dbg
export HOME=/home/lindad
cd $HOME/projects/moose
rm -r HTML
rm Coverage.*
lcov -c -i -d . -o Coverage.baseline
cd test
./run_tests -j64 --re fvkernels -t --longest-jobs 20
./run_tests -j64 --re functor_properties -t --longest-jobs 20
cd ..
lcov -c -d . -o Coverage.out
lcov -a Coverage.baseline -a Coverage.out -o Coverage.combined
genhtml Coverage.combined -o HTML
