#! /bin/bash

echo Compiling C++ Code ...
make 2>"compiler_error.txt"
if [ $? != 0 ]
then
    printf "C++ compilation failed!, Not running the R script.\n"
    printf "less compiler_error.txt to see the compiler output.\n"
else
    ./BayesCRE.R -e ./data/KB/BEL_LargeCorpus.ents -r ./data/KB/BEL_LargeCorpus.rels -d ./data/Simulations/simulated_evidence.txt -m ./data/KB/MeSH.txt -n 200000 -u 20000 -k ./output/sim_mesh_prob.txt -y ./output/sim_marg.txt -t ./output/sim_hyp.txt -o ./output/sim_mesh_app.txt
fi
