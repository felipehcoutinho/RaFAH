#!/usr/bin/env perl
#This script takes as input a folder with fasta files of phage genomes and second folder of bacteiral genomas
#and finds phage-host associations based on homology matches, tRNAs and CRISPRs
#Requires ncbi blast+, Bioperl, tRNAscan-SE, cd-hit to be installed and included in the system path

use warnings;
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Tools::GFF;

#Define default values of program parameters 
my $format = "fasta"; #extension of files in the phage and host folders to be processed
my $host_folder;
my $phage_folder;
my $prok_info_file;
my $seq_to_genome_file ;
my $prok_info_ref_file;
my $taxonomy_db = "NCBI";

#Default values of cutoffs used in the blastn searches to consider a match as valid for the homology matches, trna and crisprs
my $hml_ref_db_file;
my $min_score_hml = 0;
my $min_algn_hml = 1000;
my $min_id_hml = 85;
my $max_evalue_hml = 0.001;
my $min_qcov_hml = 0;
my $max_mismatch_hml = 1000;
my $skip_hml;

my $min_score_trna = 0;
my $min_algn_trna = 60;
my $min_id_trna = 100;
my $max_evalue_trna = 0.001;
my $min_qcov_trna = 95;
my $min_tscan_score = 30;
my $max_mismatch_trna = 0;
my $skip_trna;

my $file_CRISPR;
my $CRISPR_spacer_ref_file;
my $min_score_CRISPR = 40;
my $min_algn_CRISPR = 0;
my $min_id_CRISPR = 100; 
my $max_evalue_CRISPR = 1; 
my $min_qcov_CRISPR = 100; 
my $max_mismatch_CRISPR = 0; 
my $skip_CRISPR;

my $skip_kmer;
my $wish_models_ref_folder;
my $min_kmer_score = -10;

my $parse_only;
my $use_ref;
my $threads = "1";
my $help;


#Receive input parameters from user

GetOptions(
'hosts=s' => \$host_folder,
'phages=s' => \$phage_folder,
'prok_info_file=s' => \$prok_info_file,
'seq_to_genome_file=s' => \$seq_to_genome_file,
'prok_info_ref_file=s' => \$prok_info_ref_file,
'taxonomy_db=s' => \$taxonomy_db,
'min_bitscore_hml=s' => \$min_score_hml,
'min_alignment_hml=s' => \$min_algn_hml,
'min_identity_hml=s' => \$min_id_hml, 
'max_evalue_hml=s' => \$max_evalue_hml,
'min_qcov_hml=s' => \$min_qcov_hml,
'skip_hml' => \$skip_hml,
'min_tscan_score=s' => \$min_tscan_score,
'min_bitscore_trna=s' => \$min_score_trna,
'min_alignment_trna=s' => \$min_algn_trna,
'min_identity_trna=s' => \$min_id_trna, 
'max_evalue_trna=s' => \$max_evalue_trna,
'min_qcov_trna=s' => \$min_qcov_trna,
'skip_trna' => \$skip_trna,
'CRISPR_spacer_ref_file=s' => \$CRISPR_spacer_ref_file,
'min_bitscore_CRISPR=s' => \$min_score_CRISPR,
'min_alignment_CRISPR=s' => \$min_algn_CRISPR,
'min_identity_CRISPR=s' => \$min_id_CRISPR, 
'max_evalue_CRISPR=s' => \$max_evalue_CRISPR,
'max_mismatch_CRISPR=s' => \$max_mismatch_CRISPR,
'skip_CRISPR' => \$skip_CRISPR,
'min_kmer_score=s' => \$min_kmer_score,
'skip_kmer' => \$skip_kmer,
'threads=s' => \$threads,
'format=s' => \$format,
'parse_only' => \$parse_only,
'use_ref' => \$use_ref,
'help' => \$help
);

#help message:

die
"Usage:

perl FinHo.pl --hosts Hosts_dir/ --phages Phage_dir/

\t\tOptional parameters

\tGeneral options:
--format | Extension of the files containing genomic sequences in the Host and Phage directories (default = fasta)
--threads | Number of threads to use for blast search (default = 12)
--use_ref | Use host sequences from RefSeq  instead of custom genomes
--parse_only | Skip BLAST searches and only parse the Results (ALL the original files must be in the folder and with the same names)
--help | Prints this message and exits.

\tHomology Matches options:
--min_bitscore_hml | Minimum bitscore for considering a match (default = 0)
--min_alignment_hml | Minimum alignment length for considering a match (default = 1000nt)
--min_identity_hml | Minimum % identity for considering a match (default = 85%)
--max_evalue_hml | Maximum e-value for considering a match (default = 0.001)
--min_qcov_hml | Minimum query % coverage for considering a match (default = 0)
--skip_hml | Flag to skip Homology Match step

\ttRNA Matches options:
--min_bitscore_trna | Minimum bitscore for considering a match (default = 0)
--min_alignment_trna | Minimum alignment length for considering a match (default = 60nt)
--min_identity_trna | Minimum % identity for considering a match (default = 95%)
--max_evalue_trna | Maximum e-value for considering a match (default = 0.001)
--min_qcov_trna | Minimum query % coverage for considering a match (default = 95%)
--skip_trna | Flag to skip tRNA Match step

\tCRISPR Matches options:
--min_bitscore_CRISPR | Minimum bitscore for considering a match (default = 40)
--min_alignment_CRISPR | Minimum alignment length for considering a match (default = 0nt)
--min_identity_CRISPR | Minimum % identity for considering a match (default = 100%)
--max_evalue_CRISPR | Maximum e-value for considering a match (default = 1)
--min_qcov_CRISPR | Minimum query % coverage for considering a match (default = 100%)
--max_mismatch_CRISPR | Maximum number of mismatches for considering a match (default = 0nt)
--skip_CRISPR | Flag to skip CRISPR Match step

\tWiSH options:
--min_kmer_score | Minimum score for considering a kmer match (default = -10)
--skip_kmer | Flag to skip Kmer search step

Dependencies (must be included in the system path):
ncbi blast+
Bioperl
tRNAscan-SE
WiSH (Who Is the Host?)
" if ($help);

die "Folder of host genomes and phage genomes are required. Use --help for instructions\n" unless ($phage_folder);

my $score_index = 0;
my %seen_combos;
my %index;
my %revindex;
my %scores;
my %prok_info;
my %sequence_to_genome;
my $full_tax_var_id = $taxonomy_db."_Full_taxonomy";
my %tax_dict = qw(d Domain p Phylum c Class o Order f Family g Genus s Species);
central();

#This subroutine calls all other subroutines of the program
sub central {
	if ($use_ref) {
		print "Will look for taxonomy in $full_tax_var_id column from $prok_info_ref_file\n";
		my %params = ("file" => $prok_info_ref_file);
		my $info_ref = read_table(\%params);
		%prok_info = %$info_ref;
	} elsif ($prok_info_file) {
		print "Will look for taxonomy in $full_tax_var_id column from $prok_info_file\n";
		my %params = ("file" => $prok_info_file);
		my $info_ref = read_table(\%params);
		%prok_info = %$info_ref;
	}
	
	if ($seq_to_genome_file) {
		my %params = ("file" => $seq_to_genome_file);
		my $info_ref = read_table(\%params);
		%sequence_to_genome = %$info_ref;
	}
	
	open INDEX, "> FinHo_Seq_Index.tsv" or die "$!";
	merge_and_index($host_folder,"FinHo_Host_Sequences.fasta","Host") unless ($use_ref); 
	merge_and_index($phage_folder,"FinHo_Phage_Sequences.fasta","Phage");
	close INDEX;
	call_CRISPRDetect() unless ($skip_CRISPR);
	call_blast() unless ($skip_hml);
	call_trna_scan() unless ($skip_trna);
	#call_wish() unless ($skip_kmer);
	my %params = ("matrix" => \%scores, "file" => "Find_Host_Results.tsv");
	print_matrix(\%params);
}

sub merge_and_index {
	print "Indexing sequences from $_[0]\n";
#Merges and indexes all sequences from the phage and host folders into two distinct files. Keep track of sequence IDs and renames sequences with replicated names
	#Iterate over all fasta file
	my $dir = $_[0];
	my @files = glob("$dir*$format");
	#Count number of sequences seen
	my $counter = 0;
	#Run through all fasta files and print all sequences to $_[0]
	my $seq_out = Bio::SeqIO->new(-file => "> $_[1]", -format => "fasta");

	foreach my $file (@files) {
		my $seq_in = Bio::SeqIO->new(-file => "< $file", -format => "fasta");

		while (my $seq_obj = $seq_in->next_seq) {
			$counter++;
			my $new_id = $_[2]."_Seq_".$counter;
			$index{$new_id}{"File"} = $file; #Keep track of sequences ids in files in %index
			$index{$new_id}{"O_Id"} = $seq_obj->id;
			print INDEX $new_id."\t".$seq_obj->id."\t".$seq_obj->desc."\t".$file."\n";
			die "Repeated id $index{$new_id}{O_Id} in $file\n" if (defined $revindex{$seq_obj->id});
			$revindex{$seq_obj->id} = $new_id;
			$seq_obj->id($new_id);
			$seq_out->write_seq($seq_obj);
			
		} 
	}	

}

sub call_wish {
	print "Searching for Kmer matches\n";
	my $wish_models_folder = "FinHo_CustomWIsHModels/";
	
	if ($use_ref) {
		$wish_models_folder = $wish_models_ref_folder;
	} 
	
	unless ($parse_only) {
		unless ($use_ref) {
			#Call  WISH
			print "Building WIsH models\n";
			system("/home/felipe/felipe/WIsH-master/WIsH -c build -g $host_folder -m $wish_models_folder -t $threads");
		}
		print "Running WIsH with models from $wish_models_folder\n";
		system("mkdir FinHo_WIsHOutput");
		system("/home/felipe/felipe/WIsH-master/WIsH -c predict -g $phage_folder -m $wish_models_folder -r FinHo_WIsHOutput/ -b -t $threads")
	}
	collect_wish("Kmer","FinHo_WIsHOutput/prediction.list");
}

sub collect_wish {
	my ($etype,$input_file) = @_;
	print "Parsing $etype data from $input_file\n";

	open IN, "< $input_file" or die "$!";
	my ($phage_id,$host_id);
	
	my $header = <IN>;
	while (my $line = <IN>) {
		chomp $line;
		my @values = split /\t/, $line;
		$phage_id = $values[0];
		$host_id = $values[1];
		
		#Ignore any hits that dont satisfy the criteria for min score;
		print "Line seems to be incomplete:\n$line\n" unless (@values == 4);
		next if ($values[2] < $min_kmer_score);

		next if (defined $seen_combos{$values[0]}{$values[1]});
		$seen_combos{$host_id}{$phage_id} = 1;
		
		$score_index++;
		$scores{"Evidence_Type"}{$score_index} = $etype;
		$scores{"Phage_Sequence"}{$score_index} = $phage_id;
		
		my $phage_nid = $revindex{$phage_id};
		$scores{"Phage_Sequence_File"}{$score_index} = $index{$phage_nid}{"File"};
		
		
		$scores{"Host_Sequence"}{$score_index} = "NA";
		$scores{"Host_Sequence_File"}{$score_index} = $host_id;
		$scores{"Score"}{$score_index} = $values[2];

		#print "Host Genome: $host_id\n\tGTDB: $sequence_to_genome{GTDB_Full_taxonomy}{$host_id}\n\tNCBI: $sequence_to_genome{NCBI_Full_taxonomy}{$host_id}\n"; sleep 1;
		
		
		#print "Host Genome: $host_id\n\tFull Taxonomy: $full_taxonomy\n";
		if (defined $prok_info{$full_tax_var_id}{$host_id}) {
			my $full_taxonomy = $prok_info{$full_tax_var_id}{$host_id};
			if (($full_taxonomy ne "NA") and ($full_taxonomy ne "Unknown")) {
				my @levels = split /;/, $full_taxonomy;
				my $level_count = @levels;
				unless ($level_count == 7) {
					print "Warning $full_tax_var_id for $host_id does not have 7 levels!\n$full_taxonomy\n";
				}
				foreach my $posit (1..$level_count) {
					$posit--;
					my ($rank,$taxon) = split /__/, $levels[$posit];
					$posit++;
					my $out_var_name = "Host_Taxonomy_$posit"."_$tax_dict{$rank}";
					$scores{$out_var_name}{$score_index} = $taxon;
				}
			}
		} else {
			print "WARNING! No $full_tax_var_id defined for $host_id\n";
		}	
	}
	close IN;	
}

sub call_CRISPRDetect {
	print "Searching for CRISPR spacer matches\n";
	
	my $file_CRISPR = "FinHo_CRISPR_SpacersxVir.blastn";
	my $query_CRISPR = "FinHo_Host_CRISPR_Spacers.fasta";
	
	if ($use_ref) {
		$file_CRISPR = "FinHo_CRISPR_Spacers_RefxVir.blastn";
		$query_CRISPR = $CRISPR_spacer_ref_file;
	} 
	
	unless ($parse_only) {
		unless ($use_ref) {
			#Call  CRISPERdetect using default parameters and using maximumm number of threads available in the system. This command must be changed according to the location of the cirpsr detect script 
			system("perl /home/felipe/felipe/CRISPRDetect_2.2/CRISPRDetect.pl -f FinHo_Host_Sequences.fasta -o FinHo_Host_CRISPRDetect_Output -T 0 -array_quality_score_cutoff 3");
			print "Parsing Host_CRISPRDetect_Output.gff\n"; #Parse the gff output of CRISPR detect, fetch the sequencens anotated with a primary tag of binding_site, and print them to the fasta file Host_CRISPR_Spacers.fasta
			my $gffio = Bio::Tools::GFF->new(-file => "< FinHo_Host_CRISPRDetect_Output.gff", -gff_version => 3);
			my $out_obj = Bio::SeqIO->new(-file => "> FinHo_Host_CRISPR_Spacers.fasta", -format => "fasta");
			while(my $feature = $gffio->next_feature()) {
				my $type = $feature->primary_tag();
				next unless ($type eq "binding_site");
				my $id = $feature->seq_id();
				my $start = $feature->start();
				my $end = $feature->end();
				my ($seq) = $feature->get_tag_values('Note');
				my ($repeat_id) = $feature->get_tag_values('ID');
				my $new_seq = Bio::Seq->new(-seq => $seq, -id => "$id|$start|$end|$repeat_id");
				$out_obj->write_seq($new_seq);
				}
		}
		#Query sequences from Host_CRISPR_Spacers.fasta to the file containing viral sequences using blastn with special parameters set for CRISPR searches
		system("makeblastdb -in FinHo_Phage_Sequences.fasta -dbtype  nucl -title FinHo_VirDB -out FinHo_VirDB");
		print "Performing blast search for CRISPR Spacers\n";
		system("blastn -task blastn-short -db FinHo_VirDB -query $query_CRISPR -out $file_CRISPR -outfmt '6 std qcovs' -evalue 1 -word_size 7 -gapopen 10 -gapextend 2 -penalty -1 -dust no -max_target_seqs 100000 -num_threads $threads");
	}
	
	#Parse the output of the blast search and keep track of the matches within the specified criteria
	collect("CRISPR",$file_CRISPR,$min_id_CRISPR,$min_algn_CRISPR,$max_mismatch_CRISPR,$max_evalue_CRISPR,$min_score_CRISPR,$min_qcov_CRISPR);
}


sub call_trna_scan {
	print "Searching for tRNA matches\n";
	my %coords;
	
	my $file_trna = "FinHo_PhagexHost_tRNAs.blastn";
	my $db_hml = "FinHo_HostDB";
		
	if ($use_ref) {
		$file_trna = "FinHo_PhagexHost_tRNAs_Ref.blastn";
		$db_hml = $hml_ref_db_file;
	} 
		
	unless ($parse_only) {
		#Identify tRNA sequences in the viral sequences with tRNAscan-SE, parse the output and print the tRNA sequences to Viral_tRNAs.fasta
		system("/home/felipe/felipe/tRNAscan-SE-2.0/tRNAscan-SE -B -Q -c /home/felipe/felipe/tRNAscan-SE-2.0/tRNAscan-SE.conf --thread $threads -o FinHo_Vir_tRNAs_Report FinHo_Phage_Sequences.fasta");
	}
		
	print "Parsing tRNA_Scan Output\n";
	open IN, "< FinHo_Vir_tRNAs_Report" or die "$!";

	my $header1 = <IN>;
	my $header2 = <IN>;
	my $header3 = <IN>;

	my $rep_count = 0;
	my $valid_count = 0;
	my %valid_ids;

	while (my $line = <IN>) {
		$rep_count++;
		chomp $line;
		my @values = split /\s+/, $line;

		if (($values[4] =~ /Pseudo/) or ($values[4] =~ /Undet/)) {
			#print "Skipping: [$values[0]] tRNA [$values[1]] in:\n$line\n" 
			} else {
			$values[0] =~ s/\s//g;
			$values[1] =~ s/\D//g;
			$valid_ids{$values[0]} = 1;
			$valid_count++;
			#print "Passing: [$values[0]] in $line\n";
			$coords{$values[0]}{".trna$values[1]"}{"Start"} = $values[2];
			$coords{$values[0]}{".trna$values[1]"}{"End"} = $values[3];
			$coords{$values[0]}{".trna$values[1]"}{"Reversed"} = 1 if ($values[2] > $values[3]);
			$coords{$values[0]}{".trna$values[1]"}{"Score"} = $values[8];
		}
		
	}
	close IN;

	print "Obtained $valid_count valid tRNA reports out of $rep_count\n";

	unless ($parse_only) {
		print "Fetching tRNA sequences from Phage_Sequences.fasta\n";
		my $seq_in = Bio::SeqIO->new(-file => "< FinHo_Phage_Sequences.fasta", -format => "fasta");
		my $seq_out = Bio::SeqIO->new(-file => "> FinHo_Viral_tRNAs.fasta", -format => "fasta");
			
		while (my $in_seq_obj = $seq_in->next_seq) {
			my $id = $in_seq_obj->id();
			my $base_count = $in_seq_obj->length();
			if ($base_count == 0) {print "Sequence $id is empty and will be Skipped!\n"; next} 
			unless ($in_seq_obj->alphabet eq 'dna') {print "Skiping Weird sequence:$id\n"; next}
			if (defined $valid_ids{$id}) {
				#print "Fetching tRNAs for sequence $id\n";
				my $trna_ids_ref = $coords{$id};
				my %trna_ids_list = %$trna_ids_ref;
				foreach my $trna_id (sort keys %trna_ids_list) {
					my $trna_seq;
					my $out_seq_obj;
					if (defined $coords{$id}{$trna_id}{"Reversed"}) {
						$trna_seq = $in_seq_obj->subseq($coords{$id}{$trna_id}{"End"},$coords{$id}{$trna_id}{"Start"});
						$out_seq_obj = Bio::Seq->new(-id => $id.$trna_id, -seq => $trna_seq);
						$out_seq_obj = $out_seq_obj->revcom();
						#$out_seq_obj->seq($revcom);
					} else {
						$trna_seq = $in_seq_obj->subseq($coords{$id}{$trna_id}{"Start"},$coords{$id}{$trna_id}{"End"});
						$out_seq_obj = Bio::Seq->new(-id => $id.$trna_id, -seq => $trna_seq);
					} 
					$seq_out->write_seq($out_seq_obj);
				}
			} else {
				#print "No tRNA coords defined for: $id Skipping\n"
			}
		}
		
		
		#query Host_tRNAs.fasta against a database made from Host_tRNAs.fasta using blastn 
		print "Performing blast search\n";
		system("blastn -db $db_hml -query FinHo_Viral_tRNAs.fasta -out $file_trna -outfmt '6 std qcovs' -evalue 0.1 -max_target_seqs 100000 -num_threads $threads");
	}
	
	#Parse the output of the blast search and keep track of the matches within the specified criteria
	collect("tRNA",$file_trna,$min_id_trna,$min_algn_trna,$max_mismatch_trna,$max_evalue_trna,$min_score_trna,$min_qcov_trna);
}

sub call_blast {
	print "Searching for Homology matches\n";	
	my $file_hml = "FinHo_PhagexHost.blastn";
	my $db_hml = "FinHo_HostDB";
	
	if ($use_ref) {
		$file_hml = "FinHo_PhagexHost_Ref.blastn";
		$db_hml = $hml_ref_db_file;
	} 
	
	unless ($parse_only) {
		#Query the phage sequences from Phage_Sequences.fasta host sequences against Host_Sequences.fasta with blastn
		system("makeblastdb -in FinHo_Host_Sequences.fasta -dbtype nucl -title FinHo_HostDB -out FinHo_HostDB") unless ($use_ref);
		system("blastn -db $db_hml -query FinHo_Phage_Sequences.fasta -out $file_hml -outfmt '6 std qcovs' -evalue 0.1 -max_target_seqs 100000 -num_threads $threads");
	}
	
	#Parse the output of the blast search and keep track of the matches within the specified criteria
	collect("Homology",$file_hml,$min_id_hml,$min_algn_hml,$max_mismatch_hml,$max_evalue_hml,$min_score_hml,$min_qcov_hml);
}

sub collect {
	my ($etype,$input_file,$min_id,$min_algn,$max_mismatch,$max_evalue,$min_score,$min_qcov) = @_;
	print "Parsing $etype data from $input_file\n";

	open IN, "< $input_file" or die "$!";
	my ($phage_id,$host_id);
	
	while (my $line = <IN>) {
		chomp $line;
		my @values = split /\t/, $line;
		if ($etype eq "CRISPR") {
			$phage_id = $values[1];
			$host_id = $values[0];
			if ($use_ref) { $host_id =~ s/\|Spacer(.)*// } else {$host_id =~ s/\|(.)*$//;};
		} elsif ($etype eq "tRNA") {
			$phage_id = $values[0];
			$host_id = $values[1];
			$phage_id =~ s/\.trna(.)*$//;
		} elsif ($etype eq "Homology") {
			$phage_id = $values[0];
			$host_id = $values[1];
		}
		#Ignore any hits that dont satisfy the criteiria for alignment length, identity, bitscore and evalue;
		print "Line seems to be incomplete:\n$line\n" unless (@values == 13);
		next if ($values[2] < $min_id);
		next if ($values[3] < $min_algn);
		next if ($values[4] > $max_mismatch);
		next if ($values[10] > $max_evalue);
		next if ($values[11] < $min_score);
		next if ($values[12] < $min_qcov);
		next if (defined $seen_combos{$values[0]}{$values[1]});
		$seen_combos{$host_id}{$phage_id} = 1;
		
		$score_index++;
		$scores{"Evidence_Type"}{$score_index} = $etype;
		$scores{"Phage_Sequence"}{$score_index} = $index{$phage_id}{"O_Id"};
		$scores{"Phage_Sequence_File"}{$score_index} = $index{$phage_id}{"File"};
		#if ($use_ref) {$scores{"Host_Sequence"}{$score_index} = $prok_info{"NCBI_Access_Number"}{$host_id}} else {$scores{"Host_Sequence"}{$score_index} = $index{$host_id}{"O_Id"}}; #$host_id;}
		$scores{"Host_Sequence"}{$score_index} = $host_id;
		if ($use_ref) { $scores{"Host_Sequence_File"}{$score_index} = "Reference_File" } else {$scores{"Host_Sequence_File"}{$score_index} = $index{$host_id}{"File"} }
		$scores{"Score"}{$score_index} = $values[11];
		$scores{"Identity"}{$score_index} = $values[2];
		$scores{"Alignment_Length"}{$score_index} = $values[3];
		$scores{"Mismatches"}{$score_index} = $values[4];
		$scores{"Query_Coverage"}{$score_index} = $values[12];
		unless ($use_ref) {
			#print "Changing $host_id to ".$index{$host_id}{"O_Id"}."\n";
			$host_id = $index{$host_id}{"O_Id"};
		}
		
		my $host_nid = $host_id;
		$host_nid =~ s/\|(.)*$// if (($etype eq "CRISPR") or ($etype eq "tRNA"));
		my $genome = "NA";
		$genome = $sequence_to_genome{"Genome"}{$host_nid} if (defined $sequence_to_genome{"Genome"}{$host_nid});
		if (defined $prok_info{$full_tax_var_id}{$genome}) {
			my $full_taxonomy = $prok_info{$full_tax_var_id}{$genome};
			if (($full_taxonomy ne "NA") and ($full_taxonomy ne "Unknown")) {
				my @levels = split /;/, $full_taxonomy;
				my $level_count = @levels;
				unless ($level_count == 7) {
					print "Warning $full_tax_var_id for $host_id $host_nid from $genome does not have 7 levels!\n$full_taxonomy\n";
				}
				foreach my $posit (1..$level_count) {
					$posit--;
					my ($rank,$taxon) = split /__/, $levels[$posit];
					$posit++;
					my $out_var_name = "Host_Taxonomy_$posit"."_$tax_dict{$rank}";
					$scores{$out_var_name}{$score_index} = $taxon;
					}
				}
		} else {
			print "WARNING! No $full_tax_var_id defined for $host_id $host_nid from $genome\n";
		}
		
		
		
	}
	close IN;	
}

 
sub read_table {
#Reads a table as a 2D hash. Columns are the first dimension and rows the second. Receives a reference of a hash as the parameters to run the subroutine  
	my %parameters	= %{$_[0]};
	die "No Input File Defined.\n" unless (defined $parameters{"file"});
	my $file = $parameters{"file"};
	my $sep = "\t";
	my %blistR;
	my %blistC; 
	$sep = $parameters{"separator"} if (defined $parameters{"separator"});	
	%blistR = %{$parameters{"blistR"}} if (defined $parameters{"blistR"});
	%blistC = %{$parameters{"blistC"}} if (defined $parameters{"blistC"});
	
	my %matrix;

	print "Retrieving data from $file\n";
	open INPUT, "< $file" or die "$!";
	my $header = <INPUT>;
	$header =~ s/\r//g;
	#$header =~ s/_\d//g;
	chomp $header;
	my @cols = split /$sep/, $header;
	
	my $values_count;
	while (my $line = <INPUT>) {
		#Read each line and assign it to the 2D hash
		$line =~ s/\r//g;
		chomp $line;
		my @rows = split /$sep/, $line;
		next if (defined $blistR{$rows[0]});
		foreach my $posit (1..(@cols - 1)) {
			next if (defined $blistC{$cols[$posit]}); 			
			$matrix{$cols[$posit]}{$rows[0]} = $rows[$posit];
			$values_count++;
		}
	}
	close INPUT;
	print "Obtained $values_count values from $file.\n";
	return \%matrix;
}


sub print_matrix {
    my %parameters  = %{$_[0]};
         
    die "No Matrix (2D Hash) Defined.\n" unless (defined $parameters{"matrix"});
    my %matrix = %{$parameters{"matrix"}};  
    my $output_file = "Output_Matrix";
    my $sep = "\t"; 
     
    $output_file = $parameters{"file"} if (defined $parameters{"file"});
    $sep = $parameters{"separator"} if (defined $parameters{"separator"});      
 
     
    my $rows_ref = list_level2(\%matrix);
    my %rows = %{$rows_ref};
 
    open OUTPUT, "> $output_file" or die "$!";   
    print OUTPUT "Hit_#";
    #print the Header   
    foreach my $col (sort keys %matrix) {
        print OUTPUT "$sep"."$col";
        }
    print OUTPUT "\n";
    #print the row names followed by the values and a line break
    foreach my $row (sort keys %rows) {
        print OUTPUT "$row";
        foreach my $col (sort keys %matrix) {
            if (defined $matrix{$col}{$row}) {print OUTPUT "$sep"."$matrix{$col}{$row}"} else {print OUTPUT "$sep"."NA"};
        }
    print OUTPUT "\n";
    }
    close OUTPUT;
}

sub list_level2 {
    #Receives a reference to a 2D hash. returns a reference to a 1D hash in which the keys are all the keys from the 2nd dimesion of the input hash 
    my %hash = %{$_[0]};
    my %list;   
    foreach my $param1 (keys %hash) {
		my $hash2ref =  $hash{$param1};
		my %hash2 = %$hash2ref;
        foreach my $param2 (keys %hash2) {
            $list{$param2} = 1;
        }
    }
    return \%list;
}