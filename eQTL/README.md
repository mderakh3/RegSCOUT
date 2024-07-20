Script 1: transfer_files.sh
Purpose: This script is designed to run on your local machine. It transfers your input files to the computeCanada to the specified directory. 

Script 2: run_eqtl.sh
Purpose: This script is run on your computeCanada. It downloads all the cell subtype specific eQTL indexed files and uses tabix function to find gene associated with the input lists of SNPs.

Script 3: download_files.sh
Purpose: This script is run on your local machine. It downloads the output files from computeCanada to your local directory.

Transfer run_eqtl.sh to computeCanada through: 
scp "/Users/m.hossein_drn/Documents/aim_1/eQTL_gene_linkage/run_eqtl.sh" mderakh3@gra-login1.computecanada.ca:~/moho/}
