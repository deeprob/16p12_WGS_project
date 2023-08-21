#!/bin/bash





# combos of length2
combos=2
for phenotype in Child_ID_DD Child_behav Child_psych Child_nervous_system Child_congenital Child_craniofacial
do
Rscript 3_rarecomb.R $phenotype $combos
done


combos=3
for phenotype in Child_ID_DD Child_behav Child_psych Child_nervous_system Child_congenital Child_craniofacial
do
Rscript 3_rarecomb.R $phenotype $combos
done




