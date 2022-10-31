<p align="center">
  <img src="https://github.com/felipehcoutinho/RaFAH/blob/main/RaFAh_Logo_From_craiyon_215839_watercolour_of_bacteriophage_T4_infecting_a_bacterium.png" width="400" height="400" alt="RaFAH logo generated with Craiyon"/>
</p>

# RaFAH: Random Forest Assignment of Hosts (v0.3)

## Introduction

RaFAH is a tool for performing in silico host predictions for viruses of Archaea and Bacteria.
The rationale behind RaFAH is to use machine learning, specifically Random Forests, to assign host to a viruses based on its genomic content. For a detailed description of how RaFAH was developed and reports on its performance across multiple datasets please refer to [our publication](https://doi.org/10.1101/2020.09.25.313155)

***

## Dependencies
RaFAH is written in Perl and R and uses multiple external depencies.

In order to perform host predictions through RaFAH the following dependencies are necessary:

- BioPerl
- Prodigal (v2.60)
- HMMER (v3.1b2)
- R (v3.6.3)
- R libraries: ranger

In case you are interested in training custom models some additional dependencies are necessary:

- MMSeqs2 (Release 12-113e3)
- MUSCLE (v3.8.1551)

A docker container with all the necessary dependencies, files, and scripts is available in [the docker hub](https://hub.docker.com/repository/docker/fhcoutinho/rafah/)
In the docker container the RaFAH.pl script is located in /usr/bin/

***

## Setup

These steps only need to be performed once:
 
To setup RaFAH download and decompress all the files listed in [the repository](https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/)


These include the RaFAH scripts and the necessary associated data to run them with the default pre-computed models:
- RaFAH.pl
- HP_Ranger_Model_3_Valid_Cols.txt
- HP_Ranger_Model_3_Filtered_0.9_Valids.hmm;
- RaFAH_Predict_Host.R
- RaFAH_train_Model.R;
- MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData

Everytime you run RaFAH you will need to specify the full path of these files in your system. This is to ensure the correct files are used in the case there are multiple models in the same machine. If you are always going to use the default model or if you want to avoid having to type the full path to these files everytime you run RaFAH you can modify the script itself to always point to a given set of files by default. To do so, open the file RaFAH.pl and edit the following lines to point to the location of the previously downloaded files in your machine.

**Note: edit only the text between double quotes and nothing else.**

- Line 26: Full path to HP_Ranger_Model_3_Valid_Cols.txt
- Line 27: Full path to HP_Ranger_Model_3_Filtered_0.9_Valids.hmm
- Line 28: Full path to RaFAH_Predict_Host.R
- Line 29: Full path to RaFAH_train_Model.R
- Line 30: Full path to MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData

Save the edited version of RaFAH.pl

***

## Running RaFAH

For instruction on how to run RaFAH:

	perl RaFAH.pl --help

There are two major functions in RaFAH: train and predict.
predict will perform host predictions using a set of pre-computed model files
train will perform the necessary analysis to train a custom model and generate the necessary files to be later used with the predict function. By default all output files are written to the current working directory. You can change this by specifying the desired output directory with --output_dir.

### Predict function

We recommend that you use the toy set of 10 genomes for this tutorial. Download and decompress the files available [here](https://sourceforge.net/projects/rafah/files/Data/Toy_Set.tgz/download)

Let's start with the predict function. Starting from a set of viral genomes this function will perform gene calling with Prodigal, query those genes against the specificed Hmmer database with hmmscan, parse the output of hmmscan to generate a Genome X OG score table, and use said table to perform host prediction with pre-computed Random Forest models built with the Ranger function in R. 

Assuming you have a set of viral genomes for which you want to predict the hosts using the default pre-computed models provided with RaFAH. Put the fasta format file(s) containing genomes in a directory of your choice.  Let's call this directory Genomes for the purpose of this tutorial and assume that the extension of the files is .fasta. Also, we will include a prefix to the files so we can keep track of the results from different runs. If you have altered the script as described above you would simply run:

	perl RaFAH.pl --predict --genomes_dir Genomes/ --extension .fasta --file_prefix RaFAH_1

If you have not altered the script as described above you would need to specify the location of the files used by RaFAH so the command would be:

	perl RaFAH.pl --predict --genomes_dir Genomes/ --extension .fasta --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData

Once this runs the following files should be generated:

- RaFAH_1_CDS_Prediction.faa
- RaFAH_1_CDS_Prediction.fna
- RaFAH_1_CDS_Prediction.gff
- RaFAH_1_CDSxClusters_Prediction
- RaFAH_1_Genomes_Prediction.fasta
- RaFAH_1_Genome_to_OGs_Score_Min_Score_50-Max_evalue_1e-05_Prediction.tsv
- RaFAH_1_Host_Predictions.tsv
- RaFAH_1_Seq_Info_Prediction.tsv

These respectively correspond to:

- Protein sequences of the CDS identified in the genomes
- DNA sequences of the CDS identified in the genomes
- Gff file with the CDS coordinates
- File gereated by querying the proteins sequences from RaFAH_1_CDS_Prediction.faa against the HP_Ranger_Model_3_Filtered_0.9_Valids.hmm database with hmmscan
- Fasta file with the merged sequences from Genomes/
- Genome x OG table in tsv format generated by parsing RaFAH_1_CDSxClusters_Prediction
- Table containing all the host predictions and their host prediction scores
- Table containing the information collected for the genomic sequences during prediction (i.e. Length, CDS count, source file, description, filtered host prediction, and host prediction score)


The RaFAH_1_Seq_Info_Prediction.tsv table will include some basic information about each input sequence (i.e. Length, CDS count, source file, and description) and the host prediction generated by RaFAH, specifically the taxonomic genus (Predicted_Host column) and the score of the prediction (Predicted_Host_Score). Scores range from 0 to 1 and the closer to 1 the more confident RaFAH is of the prediction. By default, RaFAH only outputs predictions with a minimum score of 0.14 (which provides approximately 90% precision at the phylum level). To change this behaviour we can change the --min_cutoff parameter. For example, if you want RaFAH to output all predictions regardless of their score you could run:

	perl RaFAH.pl --min_cutoff 0 --predict --genomes_dir Genomes/ --extension .fasta --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData

**Note: Processing is performed per sequence, sequences in the same file are not treated as fragments of a larger genome.
 	
RaFAH supports multi-threading. You can set up the number os threads to be used through the --threads parameter. So if you want to use 24 threads run:

	perl RaFAH.pl --threads 24 --predict --genomes_dir Genomes/ --extension .fasta --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData
	
RaFAH uses Prodigal to identify protein sequences among viral genomes. If you prefer, you can provide your own pre-computed protein sequences through the --merged_cds_file_name variable. This will skip gene calling with Prodigal and directly query the protein sequences provided against the Hmmer database. Assuming your file is named Proteins.faa run:

	perl RaFAH.pl --predict --merged_cds_file_name Proteins.faa --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData

**Note: the protein sequences in the file must follow the Prodigal naming scheme, i.e. protein sequences must be named as the Id of the original genomic sequence followed by _SEQNUM. Example: Scaffold_1_1, Scaffold_1_2, Scaffold_1_3, GenomeA_1, GenomeA_2...**

**Note: --merged_cds_file_name takes as input a single file with the protein sequences from all viral genomes as opposed to a directory and an extension**

Likewise, you can also skip the hmmscan search by directly providing the output of this search with --hmmer_hits_file_name. For example:

	perl RaFAH.pl --predict --hmmer_hits_file_name Proteins_vs_HP_Ranger_Model_3_Filtered_0.9_Valids --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData
	
In the next step RaFAH parses the hmmscan output file to generate a .tsv table of Genomes (rows) x OGs (columns). This parsing step can also be skipped by directly providing this .tsv file in which the cells are filled with the Bitscores of the best hits of proteins in a given genome against each of the OGs used by the model, or 0 if there are no hits or the hits fall below the 50 bitscore threshold.

	perl RaFAH.pl --predict --genomexog_table_file_name Scores_Proteins_vs_HP_Ranger_Model_3_Filtered_0.9_Valids.tsv --valid_ogs_file HP_Ranger_Model_3_Valid_Cols.txt --hmmer_db_file_name HP_Ranger_Model_3_Filtered_0.9_Valids.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData

**Note: The columns of this table must be named according to the protein OGs listed in the _Valid_Cols.txt file specifc to the model. The table must include all columns and nothing else, otherwise RaFAH will be unable to perform the predictions**

The goal of being able to start RaFAH at different stages is to speed things up in case something goes wrong during the analysis (e.g. lack of a dependency). This way RaFAH can pick up where it left off. We do not recommend providing data generated with other tools (such as BLAST or HHsuite) because RaFAH was not designed to operate with the file formats of other tools and the results might not be comparable (e.g. different gene calling stratetgies or bitscore calculations)

This concludes the tutorial of the predict function.

***

### Train function

This function was designed for users interested in training their own Random Forest Host Prediction models. Starting from a set of viral genomes this function will perform gene calling with Prodigal, cluster the genes into orthologous groups with MMSeqs2, align the genes of each orthologous group with MUSCLE, build HMMs for each aligned orthologous group with hmmbuild, build a single Hmmer database with all HMMs with hmmpress, query the predicted genes against the generated Hmmer database with hmmscan, parse the output of hmmscan to generate a Genome X OG score table, and use said table to train a host prediction model with Ranger function in R.

Like the predict function, the train function needs a set of fasta file(s) with the genomic sequences of the viruses which will be used to train a custom model. In addition a .tsv file stating the known hosts of the genomic sequences must be provided. We can add a prefix to the output files so we known exactly which output files to use when using the model with the predict function. Assuming you have altered the script as described above and you want to use 72 threads you can run:

	perl RaFAH.pl --train --genomes_dir Genomes/ --extension .fasta --true_host Genomes_Hosts.tsv --file_prefix Custom_Model_1 --threads 72

**Note: The true_host file should have no header, just one sequence ID per line followed by the host separated by a tab**

**Note: RaFAH default models were trained by providing host information at the level of genus. You can alter that behavior by providing hosts at a different taxonomic level such as species or phylum**

If you have not altered the script as described above you would need to specify the location of the Train Script used by R, so the command would be:

	perl RaFAH.pl --train --genomes_dir Genomes/ --extension .fasta --true_host Genomes_Hosts.tsv --file_prefix Custom_Model_1 --r_script_train_file_name RaFAH_train_Model.R --threads 72
	
The following files and directories should be generated:

- Custom_Model_1_Aligned_OGs_Training (Directory)
- Custom_Model_1_CDS_Training.faa
- Custom_Model_1_CDS_Training.fna
- Custom_Model_1_CDS_Training.gff
- Custom_Model_1_CDSxClusters_Training
- Custom_Model_1_Genomes_Training.fasta
- Custom_Model_1_Genome_to_OGs_Score_Min_Score_50-Max_evalue_1e-05_Training.tsv
- Custom_Model_1_hmm_OGs_Training (Directory)
- Custom_Model_1_Merged_Training.hmm
- Custom_Model_1_Merged_Training.hmm.h3f
- Custom_Model_1_Merged_Training.hmm.h3i
- Custom_Model_1_Merged_Training.hmm.h3m
- Custom_Model_1_Merged_Training.hmm.h3p
- Custom_Model_1_MMSeqs_Training_all_seqs.fasta
- Custom_Model_1_MMSeqs_Training_cluster.tsv
- Custom_Model_1_MMSeqs_Training_rep_seq.fasta
- Custom_Model_1_Model_Training.RData
- Custom_Model_1_OGs_Training (Directory)
- Custom_Model_1_Seq_Info_Training.tsv
- Custom_Model_1_Valid_OGs_Training.tsv
- tmp_Cluster (Directory)

These respectively correspond to:

- Directory containing the files of the aligned orthologous groups
- Protein sequences of the CDS identified in the genomes
- DNA sequences of the CDS identified in the genomes
- Gff file with the CDS coordinates
- File gereated by querying the proteins sequences from Custom_Model_1_CDS_Training.faa against the Custom_Model_1_Merged_Training.hmm database with hmmscan
- Fasta file with the merged sequences from Genomes/
- Directory of MMSeqs2 intermediate files
- Genome x OG table in tsv format generated by parsing Custom_Model_1_CDSxClusters_Training with an additional column for the Host of each genome
- Directory containing the HMMs generated with hmmbuild
- HMM file generated by merging all the .hmm files from Custom_Model_1_hmm_OGs_Training/
- .h3f,.h3i, .h3m and .h3p files generated by running hmmpress on Custom_Model_1_Merged_Training.hmm
- Protein sequences of the CDS used by MMSeqs2
- Table with the CDS to Orthologous Group (OG) affiliation generated by MMSeqs2
- Protein sequences of the CDS representative of each OG generated by MMSeqs2
- The Random Forest model generated with R
- Directory containing the CDS sequences split into files according to their OG affiliation
- Table containing the information collected for the genomic sequences during training (i.e. Length, CDS count, source file, description, and host)
- List of OGs used by the custom model and that must be included as columns in the Genomes x OG table later used for prediction
- tmp directory used by MMSeqs2

Much like the predict function, the train function can also pick up from different stages, which is done by specifying different input files.

To provide a set of precomputed protein sequences run:

	perl RaFAH.pl --train --merged_cds_file_name Custom_Model_1_CDS_Training.faa --true_host Genomes_Hosts.tsv --file_prefix Custom_Model_1 --r_script_train_file_name RaFAH_train_Model.R --threads 72
	
To provide a precomputed file with the result from hmmscan run:

	perl RaFAH.pl --train --hmmer_hits_file_name Custom_Model_1_CDSxClusters_Training --true_host Genomes_Hosts.tsv --file_prefix Custom_Model_1 --r_script_train_file_name RaFAH_train_Model.R --threads 72
	
To provide a precomputed table of Genome x OG bitscores run:

	perl RaFAH.pl --train --genomexog_table_file_name Custom_Model_1_Genome_to_OGs_Score_Min_Score_50-Max_evalue_1e-05_Training.tsv --true_host Genomes_Hosts.tsv --file_prefix Custom_Model_1 --r_script_train_file_name RaFAH_train_Model.R --threads 72
	
**Note: Another usage of this feature is the possibility to edit the Genome x OG bitscores table before it is used to train the models. For example, you might want to exclude redundant OGs (i.e. those that have strong correlations with another OG) or to remove sparse OGs (i.e. those that have 0 values in most cells). This will make your models more simple and your computations run faster**

Once the train pipeline has completed you should have all the files necessary to perform predictions using the custom model. To do that use the predict function and point to the files that have previously been generated by the train function so that RaFAH knows to use the model you have created, which you can do by directly altering the RaFAH.pl script as described above, or by indicating these files on the command line. An example command for an unaltered script would be:

	perl RaFAH.pl --predict --genomes_dir Toy_Set/ --extension .fasta --valid_ogs_file Custom_Model_1_Valid_OGs_Training.tsv --hmmer_db_file_name Custom_Model_1_Merged_Training.hmm --r_script_predict_file_name RaFAH_Predict_Host.R --r_model_file_name Custom_Model_1_Model_Training.RData --threads 72 --file_prefix RaFAH

This concludes the tutorial of the train function.

***

### Fetch function

This function simply automates the process of downloading and decompressing the necessary files to run RaFAH. If you have already downloaded the files from the repository you do not need to use this function. Run it with:

	perl RaFAH.pl --fetch

After running this the following files should be in your working directory:

- HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3f
- HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3m
- HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3i
- HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3p
- HP_Ranger_Model_3_Valid_Cols.txt
- RaFAH_Predict_Host.R
- RaFAH_Train_Model.R
- MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData


