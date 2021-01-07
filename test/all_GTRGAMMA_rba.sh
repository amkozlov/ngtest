$RAXNG --parse --msa $DATADIR/dna8.fa --model $DATADIR/dna8.part --prefix ${PREFIX}_parse --seed 1 --threads $THREADS --redo --brlen linked $NGARGS

rba=${PREFIX}_parse.raxml.rba

$RAXNG --all --msa $rba --prefix $PREFIX --seed 1 --threads $THREADS --brlen linked --tree pars{5},rand{5} --bs-trees 20 --redo $NGARGS
