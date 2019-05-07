for constype in "" "STRICT" "mr80" "MRE";
do	

  $RAXNG --consense $constype --tree $DATADIR/dna2311.bsrep --prefix $PREFIX --seed 1 --threads 1 --redo

done

