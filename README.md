<p align="center">
  <img src="https://github.com/felipehcoutinho/FinHo/blob/main/FinHo_Logo.png" width="400" height="400" alt="Virathon logo generated with Craiyon"/>
</p>

# FinHo
Host-Virus matches for Bacteria and Archaea

## Citation
[Coutinho FH, Cabello-Yeves PJ, Gonzalez-Serrano R, Rosselli R, López-Pérez M, Zemskaya TI, Zakharenko AS, Ivanov VG, and
Rodriguez-Valera F. New viral biogeochemical roles revealed through metagenomic analysis of Lake Baikal. Microbiome . Microbiome; 2020;8:163.](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00936-4)

## Dependencies (must be included in the system path):
- Bioperl
- ncbi blast+
- tRNAscan-SE
- CRISPRDetect
- WIsH (Who Is the Host?)

## Usage:

perl FinHo.pl - --hosts Hosts_dir/ - --phages Phage_dir/

### General options:
- --format | Extension of the files containing genomic sequences in the Host and Phage directories (default = fasta)
- --threads | Number of threads to use for blast search (default = 12)
- --use_ref | Use host sequences from RefSeq  instead of custom genomes
- --parse_only | Skip BLAST searches and only parse the Results (ALL the original files must be in the folder and with the same names)
- --help | Prints this message and exits.

### Homology Matches options:
- --min_bitscore_hml | Minimum bitscore for considering a match (default = 0)
- --min_alignment_hml | Minimum alignment length for considering a match (default = 1000nt)
- --min_identity_hml | Minimum % identity for considering a match (default = 85%)
- --max_evalue_hml | Maximum e-value for considering a match (default = 0.001)
- --min_qcov_hml | Minimum query % coverage for considering a match (default = 0)
- --skip_hml | Flag to skip Homology Match step

### RNA Matches options:
- --min_bitscore_trna | Minimum bitscore for considering a match (default = 0)
- --min_alignment_trna | Minimum alignment length for considering a match (default = 60nt)
- --min_identity_trna | Minimum % identity for considering a match (default = 95%)
- --max_evalue_trna | Maximum e-value for considering a match (default = 0.001)
- --min_qcov_trna | Minimum query % coverage for considering a match (default = 95%)
- --skip_trna | Flag to skip tRNA Match step

### CRISPR Matches options:
- --min_bitscore_CRISPR | Minimum bitscore for considering a match (default = 40)
- --min_alignment_CRISPR | Minimum alignment length for considering a match (default = 0nt)
- --min_identity_CRISPR | Minimum % identity for considering a match (default = 100%)
- --max_evalue_CRISPR | Maximum e-value for considering a match (default = 1)
- --min_qcov_CRISPR | Minimum query % coverage for considering a match (default = 100%)
- --max_mismatch_CRISPR | Maximum number of mismatches for considering a match (default = 0nt)
- --skip_CRISPR | Flag to skip CRISPR Match step

### WIsH options:
- --min_kmer_score | Minimum score for considering a kmer match (default = -10)
- --skip_kmer | Flag to skip Kmer search step
