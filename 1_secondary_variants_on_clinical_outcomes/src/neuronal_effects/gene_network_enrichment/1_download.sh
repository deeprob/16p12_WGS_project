2_download.sh


#=====================
# download networks
#=====================

# on cluster
# cp /data5/anastasia/network_analysis_16p12/hub_genes/networks/* ~/mv_files


# # locally
# scp awt5304@uniport.bx.psu.edu:~/mv_files/* networks



#=====================
# download networks
#=====================

# on cluster
rm ~/mv_files/*
cp /data/3q29_network/RNASeq_network/brain.degnorm-ge2.prob-gept02.dat ~/mv_files
cp /data/3q29_network/RNASeq_network/brain.genes.protein-coding.txt ~/mv_files


# locally
scp awt5304@uniport.bx.psu.edu:~/mv_files/* networks





