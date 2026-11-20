#!/bin/bash

for S in human rat chicken cobra sturgeon stingray; do
	echo "Submitted $S"
  qsub -v SPECIES=$S "/PATH/TO/pbmc.mt/SNP/batch.simulation.pbs"
done

