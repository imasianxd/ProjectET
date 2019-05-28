#! /usr/bin/perl -w

# Author: Edita Karosiene, 07-01-2013

use strict;

use Getopt::Long;
use Env;
use Sys::Hostname;

#set environment variable for $NETMHCIIpan
my $NETMHCIIpan = $ENV{'NETMHCIIpan'};
# determine platform
my $UNIX = `uname -s`;
my $AR = `uname -m`;
chomp $UNIX;
chomp $AR;
my $platform_name = "$UNIX" . "_" . "$AR";
# set environment variable for $PLATFORM
my $PLATFORM = $ENV{'PLATFORM'}; 

# specify version of the method
my $version = "3.0";

# variables used to save options
my ($rdir, $inptype, $length, $tdir, $choose_chain, $chain_a, $chain_b, $alleles, $thr_aff_S, $thr_aff_W, $thr_rank_S, $thr_rank_W,
    $filter, $aff_f, $rank_f, $file, $dirty, $hlaseq, $hlaseq_a, $sort, $fast, $xls, $xlsfile, $unique, $w, $list, $help);

# associating options values with variables
&GetOptions (
'rdir:s' => \$rdir,
'tdir:s' => \$tdir,
'inptype:i' => \$inptype,
'length:s' => \$length,
'affS=s' => \$thr_aff_S,
'affW=s' => \$thr_aff_W,
'rankS=s' => \$thr_rank_S,
'rankW=s' => \$thr_rank_W,
'choose' => \$choose_chain,
'cha=s' => \$chain_a,
'chb=s' => \$chain_b,
'filter:s' => \$filter,
'rankF=s' => \$rank_f,
'affF=s' => \$aff_f,
'xlsfile:s' => \$xlsfile,
'a:s' => \$alleles,
'f:s'=>   \$file,
'hlaseq:s' => \$hlaseq,
'hlaseqA:s' => \$hlaseq_a,
's' => \$sort,
'fast' => \$fast,
'u' => \$unique,
'list' => \$list,
'xls' => \$xls,
'w' => \$w,
'dirty' => \$dirty,
'h' => \$help
);

# Define session Id
my $sId = $$;

##############################################################################################################
# Define default values for optional variables and check if provided options and files are in the right format
##############################################################################################################

# When the help option is defined
if (defined $help) {
	&usage();
	exit;
}

# working directory
if (!defined $rdir) {
	$rdir = $NETMHCIIpan; 
}
elsif (defined $rdir and $rdir eq "") {
	print "ERROR: insufficient arguments for option -rdir\nCheck ussage of the program using -h option\n";
	exit;
}
unless (-e $rdir) {
	print "ERROR: directory $rdir doesn't exist\n" ;
	exit;
}

# temporary directory
if (!defined $tdir) {
	$tdir = $ENV{'TMPDIR'} . "/tmp_$sId";
}
elsif (defined $tdir and $tdir eq ""){
	print "ERROR: insufficient arguments for option -tdir\nCheck ussage of the program using -h option\n";
	exit;
}

#When the list option was defined
if (defined $list) {
	system ("cat $rdir/data/convert_pseudo.dat | cut -f2 -d ' '");
	exit;
}

# input file format
my $input_format = "";
if (!defined $inptype) {
	$inptype = 0;
	$input_format = "FASTA";	
}
elsif ($inptype == 0) {
	$input_format = "FASTA";
}
elsif ($inptype == 1) {
	$input_format = "PEPTIDE";
}
unless ($inptype == 0 or $inptype == 1) {
	print "ERROR: input type should be 1 for peptide format and 0 for FASTA format\n";
	exit;
}

#peptide length
if (!defined $length and $inptype == 0) {
	$length = 15;
}
elsif (defined $length and $length eq "") {
	print "ERROR: insufficient arguments for option -length\nCheck ussage of the program using -h option\n";
	exit;
}
elsif (defined $length and $inptype == 0 and $length ne "") {
	if ($length < 9 ) {
		print "ERROR: peptide length should be 9 or more amino acids\n";
		exit;
	}
	
}

# affinity threshold for strong binders
if(!defined $thr_aff_S ) {
	$thr_aff_S = 50;
}
unless (&isAnumber($thr_aff_S)) {
	print "ERROR: threshold affinity for strong binders should be a number\n";
	exit;
}

# rank threshold for weak binders
if (!defined $thr_rank_W) {
	$thr_rank_W = 2;
}
unless (&isAnumber($thr_rank_W)) {
	print "ERROR: threshold rank value for weak binders should be a number\n";
	exit;
}

# affinity threshold for weak binders
if(!defined $thr_aff_W ) {
	$thr_aff_W = 500;
}
unless (&isAnumber($thr_aff_W)) {
	print "ERROR: threshold affinity for weak binders should be a number\n";
	exit;
}

# rank threshold for strong binders
if (!defined $thr_rank_S) {
	$thr_rank_S = 0.5;
}
unless (&isAnumber($thr_rank_S)) {
	print "ERROR: threshold rank value for strong binders should be a number\n";
	exit;
}

# if filter for the output is defined
if (!defined $filter) {
	$filter = 0;
}
elsif (defined $filter and $filter eq "") {
	print "ERROR:insufficient arguments for option -filter\nCheck ussage of the program using -h option\n";
	exit;
}
elsif (defined $filter and ($filter !=0 and $filter != 1)) {
	print "ERROR:wrong arguments for option -filter\nCheck ussage of the program using -h option\n";
	exit;
}

# Checking if hla sequence was defined
my $hlaseq_option=0;
if (defined $hlaseq and $hlaseq ne "") {
	$hlaseq_option = 1;
	my $first = 0;
	open (IN, "<", $hlaseq) or die "Can not open the file with MHC sequence $hlaseq $!\n";
	while (defined (my $line =<IN>)) {
		chomp $line;
		if ($line =~ m/^>/) {
			$first++;
		}
		elsif ($line !~ m/^>/ and $line !~ m/^[a-zA-Z]+$/) {
			print "ERROR: provided MHC sequence of beta chain is in a wrong format\n";
			exit;
		}
	}
	close IN;
	if ($first == 0) { 
		print "ERROR: provided MHC sequence of beta chain is not in FASTA format\n";
		exit;
	}
	# Checking if hla sequence was defined for alpha chain
	if (defined $hlaseq_a and $hlaseq_a ne "") {
		$hlaseq_option = 2;
		my $first = 0;
		open (IN, "<", $hlaseq_a) or die "Can not open the file with MHC sequence $hlaseq_a $!\n";
		while (defined (my $line =<IN>)) {
			chomp $line;
			if ($line =~ m/^>/) {
				$first++;
			}
			elsif ($line !~ m/^>/ and $line !~ m/^[a-zA-Z]+$/) {
				print "ERROR: provided MHC sequence of alpha chain is in a wrong format\n";
				exit;
			}
		}
		close IN;
		if ($first == 0) { 
			print "ERROR: provided MHC sequence of alpha chain is not in FASTA format\n";
			exit;
		}
	}
	elsif (defined $hlaseq_a and $hlaseq_a eq "") {
		print "ERROR: insufficient arguments for option -hlaseqA\nCheck ussage of the program using -h option\n";
		exit;
	}
}
elsif (defined $hlaseq_a and !defined $hlaseq) {
	print "ERROR: FASTA files for both MHC sequence chains should be specified, please enter MHC beta chain sequence in FASTA format\n";
	exit;
}		
elsif (defined $hlaseq and $hlaseq eq "") {
	print "ERROR: insufficient arguments for option -hlaseq\nCheck ussage of the program using -h option\n";
	exit;
}

# rank threshold for filtering output
unless (&isAnumber($rank_f)) {
	print "ERROR: threshold rank value for filtering should be a number\n";
	exit;
}
# affinity threshold for filtering output
unless (&isAnumber($aff_f)) {
	print "ERROR: threshold affinity value for filtering should be a number\n";
	exit;
}

# if filter is set to "Yes" or 1 (for the command line), give the default values for rank and affinity
my $filter_message ="";
if ($filter == 1) {
	if (!defined $rank_f) {
		$rank_f = 2;
	}
	if (!defined $aff_f) {
		$aff_f = 500;
	}
}
elsif ($filter == 0 and (defined $rank_f or defined $aff_f)) {
	if (!defined $w) {
		$filter_message = "# Filter was set to 0, input will not be filtered. Threshold value(s) will be ignored\n";
	}
}

# File name for xls output
my $xlsfilename;
if (!defined $xlsfile) {
	$xlsfilename = "NetMHCIIpan_out.xls";
}
elsif (defined $xlsfile and $xlsfile eq "") {
	print "ERROR: insufficient arguments for option -xlsfile\nCheck ussage of the program using -h option\n";
	exit;
}
else {
	$xlsfilename = $xlsfile;
}

# Saving the alleles in question (one allele or a comma separated file with many alleles)
my @alleles_input = ();

if (! defined $alleles and ! defined $choose_chain and ! defined $hlaseq) {
	$alleles = ("DRB1_0101");
}
elsif (! defined $alleles and defined $choose_chain) {
	if (defined $chain_a and defined $chain_b) {
		$alleles = ("HLA-" . $chain_a . "-" . $chain_b);
	}
	elsif ((defined $chain_a) and (!defined $chain_b)) {
		print "ERROR: beta chain is not defined\n";
		exit;
	}
	elsif ((!defined $chain_a) and (defined $chain_b)) {
		print "ERROR: alpha chain is not defined\n";
		exit;
	}
	else {
		print "ERROR: alpha and beta chains are not defined \n";
		exit;
	}
}
elsif (defined $alleles and defined $choose_chain) {
	print "ERROR: you should only choose one option from '-a' and '-choose'\n";
	exit;
}		
	
elsif (defined $alleles and $alleles eq "") {
	print "ERROR: insufficient arguments for option -a\nCheck ussage of the program using -h option\n";
	exit;
}

# converting allele names to the new nomenclature
my @alleles = ();
my %allele_names = ();

open(IN, "<", "$rdir/data/convert_pseudo.dat") or die "Can't open the file for converting allele names $rdir/data/convert_pseudo.dat $! \n";
while (defined (my $line = <IN>)) {
	chomp $line;
	my @tmp = split (" ", $line);
	my $old_name = $tmp[1];
	my $new_name = $tmp[0];
	$allele_names{$old_name} = $new_name;
}
close IN;

# splitting the list of alleles separated by commas into array of alleles
if (! defined $hlaseq) {
	@alleles_input = split (",", $alleles);
	foreach my $a (@alleles_input) {
		my $new_a;
		if (exists $allele_names{$a}) {
			$new_a = $allele_names{$a};
		}
		else {
			print "ERROR: Can not provide predictions for the chosen molecule, see the list of possible molecules\n";
			exit;
		}	
		push @alleles, $new_a;
	}
}
else {
	@alleles = ("USER_DEF");
}

# Testing if the input file is in FASTA format when no type specified and in peptide format when peptide type is specified
my %expected_affinities = ();
my $flag_expected_aff = 0;

if ($inptype == 0) {
	open (IN, "<", $file) or die "Can not open input file $! \n";
	my $first = 0;
	while ($first == 0 and defined (my $input_line =<IN>)) {
		if ($input_line !~ m/^>/) {
			print "ERROR: Input file is not in FASTA format\n";
			exit;
		}
		$first = 1;
	}
	close IN;
}

# If the input is in PEPTIDE format, check if all the peptides are of the same length
elsif ($inptype == 1) {
	my @lengths_peptides = ();
	open (IN, "<", $file) or die "Can not open input file $! \n";
	while (defined (my $input_line =<IN>)) {
		chomp $input_line;
		if ($input_line !~/^\s/) {
			my @tmp = split (" ", $input_line);
			my $pep = $tmp[0];
			if (defined ($tmp[1]) and $tmp[1] ne "") {
				$expected_affinities{$pep} = $tmp[1];
				$flag_expected_aff = 1;
			}		
					
			if ($flag_expected_aff == 1 and defined $pep and !defined $tmp[1]) {
				print "ERROR: Each line in the file should have the same format: peptide and binding affinity or just peptide sequence\n";
				exit;
			}	
			if ($pep !~ m/^[A-Za-z]+$/) {
				print "ERROR: Wrong format of the input file - check if the file contains peptides\n";
				exit;
			}
		}		
	}
	close IN;
}

# input file	
if (! defined $file and ! defined $help) {
	print "ERROR: no input data\n";
	print "Usage: ./NetMHCIIpan-3.0.pl [-h] [args] -f [fastafile/peptidefile]\n";
	exit;
}
elsif (defined $file and $file eq ""){
	print "ERROR: insufficient arguments for option -f\nCheck ussage of the program using -h option\n";
	exit;
}
elsif (defined $file) {
	unless (-e $file) {
		print "ERROR: input file $file doesn't exist\n";
		exit;
	}
}

# Check if fast mode was specified and choose the correct synapse list file
my $synlist_file = "";
if (defined $fast) {
	$synlist_file = "$rdir/data/synlist.all.short";
}
else {
	$synlist_file = "$rdir/data/synlist.all";
}

# creating temporary directory
$tdir = "$tdir";
mkdir "$tdir";
unless (-e $tdir) {
	print "ERROR: directory $tdir doesn't exist\n";
	exit;
}

# Set the print of the results to 0
my $print_count = 0;

# Putting all the pseudosequences into the hash so that we have a look up table for each allele (used for ranks and etc.)
my %pseudoseqs = ();
open (IN, "<", "$rdir/data/pseudosequences.dat") or die "Can not open the file $rdir/data/pseudosequences.dat $! \n";
while (defined (my $line =<IN>)) {
	chomp $line;
	my @tmp = split (" ", $line);
	my $allele = $tmp[0];
	my $seq = $tmp[1];
	$pseudoseqs{$allele} = $seq;
}
close IN;

# Defining hashes to be used for saving the results for making xls file
my %pos = ();
my %pep = ();
my %prot_id = ();
my %log = ();
my %nm = ();
my %bl = ();
my %rank = ();

############################################################
### If MHC sequence has been uploaded - hlaseq option(s) ###
############################################################
if ($hlaseq_option != 0) {
	my $allele = "USER_DEF";
	my $rank_flag = 0;
	my $RESULT = "";
	
	# Getting pseudosequence for running the predictions
	my $pseudosequence = "";
	my $pseudoseqA = "";
	my $pseudoseqB = "";
	if ($hlaseq_option == 1) {
		$pseudoseqB = `cat $hlaseq | $PLATFORM/bin/mhcfsa2psseq -p $rdir/data/full_beta.positions -m $rdir/data/BLOSUM30 -gf 14 -r $rdir/data/DRB10101.fsa -n -- | grep -v '#' | cut -f2 -d ' '`;
		$pseudoseqA = `cat $rdir/data/DRA.pseudo`;
		chomp $pseudoseqA;
		chomp $pseudoseqB;		
	}
	elsif ($hlaseq_option == 2) {
		$pseudoseqA = `cat $hlaseq_a | $PLATFORM/bin/mhcfsa2psseq -p $rdir/data/full_alpha.positions -m $rdir/data/BLOSUM30 -gf 14 -r $rdir/data/DRA10101.fsa -n -- | grep -v '#' | cut -f2 -d ' '`;
		$pseudoseqB = `cat $hlaseq | $PLATFORM/bin/mhcfsa2psseq -p $rdir/data/full_beta.positions -m $rdir/data/BLOSUM30 -gf 14 -r $rdir/data/DRB10101.fsa -n -- | grep -v '#' | cut -f2 -d ' '`;
		chomp $pseudoseqA;
		chomp $pseudoseqB;
		
	}
	
	if ($pseudoseqA eq "" or $pseudoseqB eq "") {
		print "ERROR: an error occured when getting pseudosequence from MHC fasta sequence\n";
		exit;
	}
	else {
		$pseudosequence = $pseudoseqA . $pseudoseqB;
	}
	if (length($pseudosequence) != 34) {
		print "ERROR: pseudosequence was not obtained correctly: $pseudosequence\n";
		exit;
	}
	
	my $new_pseudosequence = &changeGap($pseudosequence);	
	
	open (OUT, ">", "$tdir/pseudosequence_$$") or die "Can not open the file for writing '$tdir/pseudosequence_$$' $! \n";
	print OUT "USER_DEF $new_pseudosequence\n";
	close OUT;
	
	# Split FASTA files
	if ($inptype == 0) {
		my $FASTA_files = &SplitFiles($file);
		
		foreach my $fasta_file (@{$FASTA_files}) {
		
			
			# Preparing data to submit for predictions
			my ($identity, $input_file) = &PrepareInput($allele, $inptype, $flag_expected_aff, $fasta_file, $tdir);
	
			# Running the prediction script
			my $output_file = "$tdir/output_$$.out";
			system ("cat  $input_file | $PLATFORM/bin/smm_align_flank_gd_nn_play_pssm_pan2 -aX -mhc $tdir/pseudosequence_$$ -rsc -elpfr -eplen -fl 3 -blr -bl -afs $synlist_file -- | grep -v '#' > $output_file");
			
		
			# Accesing the results and preparing them to print as an output
			$RESULT = &AccessResults($allele, $output_file, $flag_expected_aff, \%expected_affinities, $identity, $rdir, $thr_aff_S, $thr_aff_W, $thr_rank_S, $thr_rank_W, $rank_flag, 0);
			
			# Saving the results into hashes for creation of xls file
			if (defined $xls) {
				my @result_lines = split ("\n", $RESULT);
				foreach my $result_line (@result_lines) {
					chomp $result_line;
					my @tmp = split (" ", $result_line);
					push (@{$pos{"USER_DEF"}}, $tmp[0]);
					push (@{$pep{"USER_DEF"}}, $tmp[2]);
					push (@{$prot_id{"USER_DEF"}}, $tmp[3]);
					push (@{$log{"USER_DEF"}}, $tmp[6]);
					push (@{$nm{"USER_DEF"}}, $tmp[7]);
		
					if (defined $tmp[8] and $flag_expected_aff !=1) {
						push (@{$bl{"USER_DEF"}}, 1);
					}
					elsif (defined $tmp[9]) {
						 push (@{$bl{"USER_DEF"}}, 1);
					}
					else {
						push (@{$bl{"USER_DEF"}}, 0);
					}
				}
			}
	
			# Filtering the results if the option for printing only the best binding core was specified
			if (defined $unique) {
				($RESULT) = &FindBestCore($RESULT, $flag_expected_aff, $rank_flag);
			}
			
			# Modifying the results if the filtering or the sorting was used
			my ($FINAL_RESULT, $count_strong, $count_weak) = &ModifyResults($RESULT, $filter, $rank_f, $aff_f, $tdir, $sort, $rank_flag);

			# Printing the results
			&Print($allele, $FINAL_RESULT, $count_strong, $count_weak, $input_format, $pseudosequence, $thr_rank_S, $thr_rank_W, $rank_f, $aff_f, $filter_message, $rank_flag);	
		}
	}
	else{ 

		# Preparing data to submit for predictions
		my ($identity, $input_file) = &PrepareInput($allele, $inptype, $flag_expected_aff, $file, $tdir);
	
		# Running the prediction script
		my $output_file = "$tdir/output_$$.out";
		system ("cat  $input_file | $PLATFORM/bin/smm_align_flank_gd_nn_play_pssm_pan2 -aX -mhc $tdir/pseudosequence_$$ -rsc -elpfr -eplen -fl 3 -blr -bl -afs $synlist_file -- | grep -v '#' > $output_file");
		#unlink "$tdir/pseudosequence_$$";
	
		# Accesing the results and preparing them to print as an output
		$RESULT = &AccessResults($allele, $output_file, $flag_expected_aff, \%expected_affinities, $identity, $rdir, $thr_aff_S, $thr_aff_W, $thr_rank_S, $thr_rank_W, $rank_flag, 0);
	
		# Saving the results into hashes for creation of xls file
		if (defined $xls) {
			my @result_lines = split ("\n", $RESULT);
			foreach my $result_line (@result_lines) {
				chomp $result_line;
				my @tmp = split (" ", $result_line);
				push (@{$pos{"USER_DEF"}}, $tmp[0]);
				push (@{$pep{"USER_DEF"}}, $tmp[2]);
				push (@{$prot_id{"USER_DEF"}}, $tmp[3]);
	  			push (@{$log{"USER_DEF"}}, $tmp[6]);
				push (@{$nm{"USER_DEF"}}, $tmp[7]);
		
				if (defined $tmp[8] and $flag_expected_aff !=1) {
					push (@{$bl{"USER_DEF"}}, 1);
				}
				elsif (defined $tmp[9]) {
				 	push (@{$bl{"USER_DEF"}}, 1);
				}
				else {
					push (@{$bl{"USER_DEF"}}, 0);
				}
			}
		}
		# Filtering the results if the option for printing only the best binding core was specified
		if (defined $unique) {
			($RESULT) = &FindBestCore($RESULT, $flag_expected_aff, $rank_flag);
		}
		
		# Modifying the results if the filtering or the sorting was used
		my ($FINAL_RESULT, $count_strong, $count_weak) = &ModifyResults($RESULT, $filter, $rank_f, $aff_f, $tdir, $sort, $rank_flag);

		# Printing the results
		&Print($allele, $FINAL_RESULT, $count_strong, $count_weak, $input_format, $pseudosequence, $thr_rank_S, $thr_rank_W, $rank_f, $aff_f, $filter_message, $rank_flag);
	}
	unlink "$tdir/pseudosequence_$$";
}	
	
##########################################
### If allele or allele list was given ###
##########################################
if ($hlaseq_option == 0){
	foreach my $allele (@alleles) {
		my $rank_flag = 1;
		my $RESULT = "";
		
		# Split FASTA files
		if ($inptype == 0) {
			my $FASTA_files = &SplitFiles($file);
		
			foreach my $fasta_file (@{$FASTA_files}) {	
				# Preparing the data to submit for predictions
				my ($identity, $input_file) = &PrepareInput($allele, $inptype, $flag_expected_aff, $fasta_file, $tdir);
			
				# Running the prediction script
				my $output_file = "$tdir/output_$$.out";
				system ("cat  $input_file | $PLATFORM/bin/smm_align_flank_gd_nn_play_pssm_pan2 -aX -mhc $rdir/data/pseudosequences.dat -rsc -elpfr -eplen -fl 3 -blr -bl -afs $synlist_file -- | grep -v '#' > $output_file");
		
				# Accesing the results and preparing them to print as an output
				$RESULT = &AccessResults($allele, $output_file, $flag_expected_aff, \%expected_affinities, $identity, $rdir, $thr_aff_S, $thr_aff_W, $thr_rank_S, $thr_rank_W, $rank_flag, \%pseudoseqs);
			
				# Saving the results into hashes for creation of xls file
				if (defined $xls) {		
					my @lines = split ("\n", $RESULT);
		
					foreach my $result_line (@lines) {
						chomp $result_line;
						my @tmp = split (" ", $result_line);			
						push (@{$pos{$allele}}, $tmp[0]);
						push (@{$pep{$allele}}, $tmp[2]);
						push (@{$prot_id{$allele}}, $tmp[3]);
						push (@{$log{$allele}}, $tmp[6]);
						push (@{$nm{$allele}}, $tmp[7]);
						push (@{$rank{$allele}}, $tmp[8]);
						if (defined $tmp[9] and $flag_expected_aff !=1) {
							push (@{$bl{$allele}}, 1);
						}
						elsif (defined $tmp[10]) {
							push (@{$bl{$allele}}, 1);
						}
						else {
							push (@{$bl{$allele}}, 0);
						}			
					}
				}
	
				# Filtering the results if the option for printing only the best binding core was specified
				if (defined $unique) {
					($RESULT) = &FindBestCore($RESULT, $flag_expected_aff, $rank_flag);
				}
				
				# Modifying the results if the filtering or the sorting was used
				my ($FINAL_RESULT, $count_strong, $count_weak) = &ModifyResults($RESULT, $filter, $rank_f, $aff_f, $tdir, $sort, $rank_flag);
		
				# Printing the results
				&Print($allele, $FINAL_RESULT, $count_strong, $count_weak, $input_format, 0, $thr_rank_S, $thr_rank_W, $rank_f, $aff_f, $filter_message, $rank_flag);
			}
		}
		else {
			# Preparing the data to submit for predictions
			my ($identity, $input_file) = &PrepareInput($allele, $inptype, $flag_expected_aff, $file, $tdir);
		
			# Running the prediction script
			my $output_file = "$tdir/output_$$.out";
			system ("cat  $input_file | $PLATFORM/bin/smm_align_flank_gd_nn_play_pssm_pan2 -aX -mhc $rdir/data/pseudosequences.dat -rsc -elpfr -eplen -fl 3 -blr -bl -afs $synlist_file -- | grep -v '#' > $output_file");
		
			# Accesing the results and preparing them to print as an output
			$RESULT = &AccessResults($allele, $output_file, $flag_expected_aff, \%expected_affinities, $identity, $rdir, $thr_aff_S, $thr_aff_W, $thr_rank_S, $thr_rank_W, $rank_flag, \%pseudoseqs);
						
			# Saving the results into hashes for creation of xls file
			if (defined $xls) {		
				my @lines = split ("\n", $RESULT);
			
				foreach my $result_line (@lines) {
					chomp $result_line;
					my @tmp = split (" ", $result_line);			
					push (@{$pos{$allele}}, $tmp[0]);
					push (@{$pep{$allele}}, $tmp[2]);
					push (@{$prot_id{$allele}}, $tmp[3]);
					push (@{$log{$allele}}, $tmp[6]);
					push (@{$nm{$allele}}, $tmp[7]);
					push (@{$rank{$allele}}, $tmp[8]);
					if (defined $tmp[9] and $flag_expected_aff !=1) {
						push (@{$bl{$allele}}, 1);
					}
					elsif (defined $tmp[10]) {
						push (@{$bl{$allele}}, 1);
					}
					else {
						push (@{$bl{$allele}}, 0);
					}			
				}
			}
			# Filtering the results if the option for printing only the best binding core was specified
			if (defined $unique) {
				($RESULT) = &FindBestCore($RESULT, $flag_expected_aff, $rank_flag);
			}		
			# Modifying the results if the filtering or the sorting was used
			my ($FINAL_RESULT, $count_strong, $count_weak) = &ModifyResults($RESULT, $filter, $rank_f, $aff_f, $tdir, $sort, $rank_flag);
		
			# Printing the results
			&Print($allele, $FINAL_RESULT, $count_strong, $count_weak, $input_format, 0, $thr_rank_S, $thr_rank_W, $rank_f, $aff_f, $filter_message, $rank_flag);
		}
	}
}

#######################################
### Preparing XLS file if specified ###
#######################################
# Acessing and printing  the results into the tab separated file for .xls output if the option "xls" option was specified
if (defined $xls) {
	foreach my $allele1 (@alleles) {
		foreach my $allele2 (@alleles) {
			if ($#{$pos{$allele1}} != $#{$pos{$allele2}}) {
				print "ERROR occured when creating Excel file: $allele1 and $allele2 resulted into different number of peptides\n";
				exit;
			}
		}
	}		
	my $first = "FALSE";
	my $file_name;
	if (defined $w) {
		$file_name = "/usr/opt/www/pub/CBS/services/NetMHCIIpan-3.0/tmp/$$"."_NetMHCIIpan.xls";
	}
	else {
		$file_name = $xlsfilename;
	}
	open (OUT, ">", $file_name) or die "Can not open the file $!\n";
	foreach my $query_allele (@alleles) {
		if ($first eq "FALSE") {
			print OUT "\t\t\t\t$query_allele";
			$first = "TRUE";
		}
		else {
			print OUT "\t\t\t$query_allele";
		}
	}
	$first = "FALSE";
	foreach my $query_allele (@alleles) {
		if ($first eq "FALSE") {
			if (!exists $rank{$alleles[0]}) {
				print OUT "\nPos\tPeptide\tID\t1-log50k\tnM";
			}
			else {
				print OUT "\nPos\tPeptide\tID\t1-log50k\tnM\tRank";
			}			
			$first = "TRUE";
		}
		else {
			print OUT "\t1-log50k\tnM\tRank";
		}
	}
	print OUT "\tAve\tNB\n";

	for (my $i = 0; $i <= $#{$pos{$alleles[0]}}; $i++) {
		my $nb = $bl{$alleles[0]}[$i];
		my $log_sum = $log{$alleles[0]}[$i];
		if (!exists $rank{$alleles[0]}) {
			print OUT "$pos{$alleles[0]}[$i]\t$pep{$alleles[0]}[$i]\t$prot_id{$alleles[0]}[$i]\t$log{$alleles[0]}[$i]\t$nm{$alleles[0]}[$i]\t";
		}
		else {
			print OUT "$pos{$alleles[0]}[$i]\t$pep{$alleles[0]}[$i]\t$prot_id{$alleles[0]}[$i]\t$log{$alleles[0]}[$i]\t$nm{$alleles[0]}[$i]\t$rank{$alleles[0]}[$i]\t";
		}
		for (my $n = 1; $n <= $#alleles; $n++){
			$nb += $bl{$alleles[$n]}[$i];
			$log_sum += $log{$alleles[$n]}[$i];
			print OUT "$log{$alleles[$n]}[$i]\t$nm{$alleles[$n]}[$i]\t$rank{$alleles[$n]}[$i]\t";
		}
		my $avg = $log_sum/scalar(@alleles);
		$avg = sprintf("%.4f", $avg);
		print OUT "$avg\t$nb\n";
	}
	close OUT;
	my $short_name = "/services/NetMHCIIpan-3.0/tmp/$$"."_NetMHCIIpan.xls";
	if (defined $w) {
		print "Link to output xls file <a href='$short_name'>NetMHCIIpan_out.xls</a>\n";
	}	
}
#################### Finished creating .xls file #################################################
	
# Deleting temporary directory if the dirty mode was not chosen
if (!defined $dirty) {
	system ("rm -r $tdir");
}
else {
	print "\n\nTemporary files have been saved here: $tdir\n";
}
###################################
###   S U B R O U T I N E S	###
###################################

## check for integer
sub isInt {
    my $test=shift;
    
    if ($test =~ m/^\d+$/ && $test>=0) {
	return 1; }
    else {
	return 0; }
}

## check format of a number
sub isAnumber {
    my $test = shift;

    eval {
        local $SIG{__WARN__} = sub {die $_[0]};
        $test += 0;
    };
    if ($@) {
	return 0;}
    else {
	return 1;} 
}

## if there is a gap in a pseudosequence, change the position of it
sub changeGap {
    my $seq = shift;
    my $new_seq = "";
    my $no_gaps = 0;
    for (my $i = 0; $i < length($seq); $i++) {
    	my $symbol = substr ($seq, $i, 1);
	if ($symbol eq "-") {
		$no_gaps = 1;
	}
	
    }
    if ($no_gaps == 1) {
    	for (my $i = 0; $i < length($seq); $i++) {
    		my $symbol = substr ($seq, $i, 1);
		if ($i < 6 or $i > 6 ){
			if ($symbol ne "-") {
				$new_seq .= $symbol;
			}
		}
		elsif ($i == 6) {
			$new_seq .= "$symbol-";			
		}
	}
    }
    else {	
    	$new_seq = $seq;
}	
    return $new_seq;
}

sub SplitFiles {
	my ($file) = @_;
	my $flag = 0;
	my @fasta_files_array = ();
	my $fasta_file;
	open (IN, "<", $file) or die "Can not open input file $file for reading $! \n";
	my $unique_id_count = 0;
	while (defined (my $line = <IN>)) {
		chomp $line;
		
		if ($line =~ m/^>(\S+)\s*/) {## MN has changed  m/^>(\S+)\s+/ to  m/^>(\S+)\s*/	04092013
			$unique_id_count++;
			if ($flag == 0) {
				$fasta_file = $tdir ."/". $1 . "_" . $unique_id_count;
				open (OUT, ">", $fasta_file) or die "Can not open input file $fasta_file for writing $! \n";
				print OUT "$line\n";
				$flag++;
			}
			else {
				close OUT;
				push (@fasta_files_array, $fasta_file);
				$fasta_file = $tdir . "/". $1 . "_" . $unique_id_count;
				open (OUT, ">", $fasta_file) or die "Can not open input file $fasta_file for writing $! \n";
				print OUT "$line\n";	
			}
						
		}
		if ($line !~ m/^>/) {
			print OUT "$line\n";
		}		
	}
	close OUT;
	push (@fasta_files_array, $fasta_file);
	return (\@fasta_files_array);		
}
sub PrepareInput {
	my ($allele, $inptype, $flag_expected_aff, $file, $tdir) = @_;
	my $identity = "Sequence";
	# If the input was in PEPTIDE format
	if ($inptype == 1) {
		open (OUT, ">", "$tdir/input_file.dat") or die "Can not open the file '$tdir/input_file.dat' for writing $!\n";
		open (IN, "<", $file) or die "Can not open the file $file for reading $! \n";
		while (defined (my $line = <IN>)) {
			chomp $line;
			if ($flag_expected_aff == 1 and $line =~ m/^(\w+)\s+(\d+\.*\d*)/) {
				my $pep = $1;
				my $aff = $2;
				if (length($pep) < 9 ) {
					print "ERROR: peptide length should be 9 or more amino acids\n";
					exit;
				}
				print OUT "$pep\t$aff\t$allele\n";				
			}
			elsif ($flag_expected_aff == 0) {
				if (length($line) < 9 ) {
					print "ERROR: peptide length should be 9 or more amino acids\n";
					exit;
				}
				print OUT "$line\t2.000\t$allele\n";
			}
		}
		close IN;
		close OUT;
	}	
	if ($inptype == 0) {
		my $seq = "";
		open (IN, "<", $file) or die "Can not open input file $file for reading $! \n";
		while (defined (my $line = <IN>)) {
			chomp $line;
			if ($line !~ m/^>/) {
				$seq .= $line;
			}
			if ($line =~ m/^>(\S+)\s*/) {
				$identity = $1;
			}
		}
		close IN;
		open (OUT, ">", "$tdir/input_file.dat") or die "Can not open the file '$tdir/input_file.dat' for writing $!\n";
		for (my $i=0; $i < length($seq) - $length + 1; $i++) {
			my $pep = substr ($seq, $i, $length);
			print OUT "$pep\t2.000\t$allele\n";
		}
		close OUT;
	}
	return ($identity, "$tdir/input_file.dat");
}

sub AccessResults {
	my ($allele, $output_file, $flag_expected_aff,  $expected_affinities, $identity, $rdir, $thr_aff_S, $thr_aff_W, $thr_rank_S, $thr_rank_W , $rank_flag, $pseudoseqs) = @_;
	my $RESULT = "";
	my $rank = "";
		
	open (IN, "<", $output_file) or die "Can not open the file $output_file for reading $! \n";
	while (defined (my $line = <IN>) ) {
		chomp $line;		
		if ($line =~ m/Error\. Cannot find MHC name/) {
			print "ERROR: Can not provide predictions for the chosen molecule, see the list of possible molecules\n";
			exit;
		}
		if ($line =~ m/Error\. Wrong line format/) {
			print "ERROR: Input contains lines in a wrong format. Check if you have empty lines\n";
			exit;
		}
		if ($line =~ m/Error\. Not elements in trainlist/) {
			print "ERROR: An error occurred creating input file.\n";
			exit;
		}
		my @tmp = split (" ", $line);
		my $pos = $tmp[0];
		my $core = $tmp[1];
		my $score = $tmp[4];
		my $peptide = $tmp[5];
		my $output_allele = $tmp[6];
		
		my $length = length($peptide);
		if ($length < 9 ) {
			print "ERROR: peptide length should be 9 or more amino acids\n";
			exit;
		}					
		# Finding the position of the core in the peptide
		my $pos_core;
		for (my $i = 0; $i < length($peptide) - 9 + 1; $i++) {
			my $new_core = substr ($peptide, $i, 9);
			if ($new_core eq $core) {
				$pos_core = $i;
			}
		}
		
		# Check if expected affinities for peptides were specified and prepare them for printing if they were
		my $expected = "";			
		if ($flag_expected_aff == 1) {
			$expected = $expected_affinities->{$peptide};
			$expected = sprintf("%12.6f", $expected);
		}	
			
		# Calculating affinity from the log score
		my $aff = exp((1-$score)*log(50000));		
		
		# Changing format for printing
		$aff = sprintf("%12.2f", $aff);
		$score = sprintf("%4.3f", $score);
				
		# Finding the rank if alleele name was given and not MHC sequence
		if ($rank_flag == 1) {
			
			my @RANKS = ();
			my %SCORES = ();
			
			my $allele_sequence = $pseudoseqs->{$allele};
			open (IN1, "<", "$rdir/data/thresholds/$allele_sequence.thr_cmb") or die "Can not open the file $rdir/data/thresholds/$allele_sequence.thr_cmb $! \n";
			while (defined (my $line =<IN1>)) {
				chomp $line;
				my @tmp = split (" ", $line);
				push (@RANKS, $tmp[1]);
				if ($length <= 19) {
					push (@{$SCORES{$length}}, $tmp[$length-7]);
				}
				else {
					push (@{$SCORES{$length}}, $tmp[12]);
				}				
			}
			close IN1;
			my $flag = 0;
			for (my $i = 0; $i <= $#RANKS; $i++) {
				if ($score >= $SCORES{$length}[$i] and $flag == 0) {
					$flag = 1;
					$rank = $RANKS[$i];
				}
				if ($i == $#RANKS and $score < $SCORES{$length}[$i]) {
					$rank = $RANKS[$#RANKS];
				}
			}
			$rank = sprintf("%6.2f", $rank);
		}
		
		## Finding the level of binding
		my $level ="";
		if ($rank_flag == 1) {
			if (($aff <= $thr_aff_S) or ($rank <= $thr_rank_S)) {
				$level = "<=SB";
			}
			elsif (($aff <= $thr_aff_W and $aff > $thr_aff_S) or ($rank <= $thr_rank_W and $rank > $thr_rank_S)) {
				$level = "<=WB";
			}
		}
		else {
			if ($aff <= $thr_aff_S) {
				$level = "<=SB";			
			}
			elsif ($aff <= $thr_aff_W and $aff > $thr_aff_S) {
				$level = "<=WB";			
			}
		}
		
		# Saving the results into one variable
		if ($rank_flag == 1) {
			$RESULT .= "$pos $output_allele $peptide $identity $pos_core $core $score $aff $rank $expected $level\n";
		}
		else {
			$RESULT .= "$pos $output_allele $peptide $identity $pos_core $core $score $aff $expected $level\n";
		}
	}
	return $RESULT;
}

sub ModifyResults {
	my ($RESULT, $filter, $rank_f, $aff_f, $tdir, $sort, $rank_flag) = @_;
	my $RESULT_MOD = "";
	my $count_strong = 0;
	my $count_weak = 0;
	my $FINAL_RESULT = "";
	
	my @result_lines = split ("\n", $RESULT);
	foreach my $result_line (@result_lines) {
		my $rank = "";
		my $level = "";
		my $expected = "";
		my @scores = split (" ", $result_line);
		my $pos = sprintf("%6s", $scores[0]);
		my $pos_core = sprintf("%4s", $scores[4]);
		my $output_allele = sprintf("%22s", $scores[1]);
		my $peptide = sprintf("%19s", $scores[2]);
		my $identity = sprintf("%11s", $scores[3]);
		my $core = sprintf("%12s", $scores[5]);
		my $score = sprintf("%4.3f", $scores[6]);
		my $aff = sprintf("%12.2f", $scores[7]);
		if ($rank_flag == 1) {
			$rank = sprintf("%6.2f", $scores[8]);
		}
		
		$aff = sprintf("%13s", $aff);
		$score = sprintf("%13s", $score);
		
		if ($rank_flag == 1) {
			if (defined $scores[10]) {
				$expected = sprintf("%12.3f", $scores[9]);
				$level = $scores[10];
			}
			elsif (defined $scores[9] and !defined $scores[10] and $flag_expected_aff != 1) {
				$level = $scores[9];
			}
			elsif (defined $scores[9] and !defined $scores[10] and $flag_expected_aff == 1) {
				$expected = sprintf("%12.3f", $scores[9]);
			}	
		}
		else {
			if (defined $scores[9]) {
				$expected = sprintf("%12.3f", $scores[8]);
				$level = $scores[9];
			}
			elsif (defined $scores[8] and !defined $scores[9] and $flag_expected_aff != 1) {
				$level = $scores[8];
			}
			elsif (defined $scores[8] and !defined $scores[9] and $flag_expected_aff == 1) {
				$expected = sprintf("%12.3f", $scores[8]);
			}
		}
		
		# Finding the number of strong and weak binders
		if ($level eq "<=SB") {
			$count_strong++;			
		}
		elsif ($level eq "<=WB") {
			$count_weak++;			
		}
		
		#If the filter for filtering output was defined
		if ($filter == 1) {
			if ($rank_flag == 1) {
				if ($rank <= $rank_f or $aff <= $aff_f) {
					$RESULT_MOD .= "$pos $output_allele $peptide $identity $pos_core $core $score $aff $rank $expected $level\n";
				}
				if (($level eq "<=SB") and ($rank > $rank_f and $aff > $aff_f)) {
					$count_strong--;
				}
				if (($level eq "<=WB") and ($rank > $rank_f and $aff > $aff_f)) {
					$count_weak--;
				}
			}
			else {
				if ($aff <= $aff_f) {
					$RESULT_MOD .= "$pos $output_allele $peptide $identity $pos_core $core $score $aff $expected $level\n";
				}
				if (($level eq "<=SB") and ($aff > $aff_f)) {
					$count_strong--;
				}
				if (($level eq "<=WB") and ($aff > $aff_f)) {
					$count_weak--;
				}
			}			
		}
		else {
			if ($rank_flag == 1) {
				$RESULT_MOD .= "$pos $output_allele $peptide $identity $pos_core $core $score $aff $rank $expected $level\n";
			}
			else {
				$RESULT_MOD .= "$pos $output_allele $peptide $identity $pos_core $core $score $aff $expected $level\n";
			}
		}
	}
		
	# Printing result lines into the file
	open (OUT, ">", "$tdir/results.out") or die "Can not open the file $tdir/results.out for writing $!\n";
	print OUT $RESULT_MOD;
	close OUT;

	# If the sort option was specified, sort the results based on affinity	
	if (defined $sort)  {
		system("cat $tdir/results.out | sort -nrk7 > $tdir/final_results.out");
		open (IN, "<", "$tdir/final_results.out") or die "Can not open the file $tdir/final_results.out for reading $!\n";
		while (defined (my $line = <IN>)) {
			$FINAL_RESULT .= $line;
		}
		close IN;
	}
	else {
		$FINAL_RESULT = $RESULT_MOD;
	}
	return ($FINAL_RESULT, $count_strong, $count_weak);
}

sub Print {
	my ($allele, $FINAL_RESULT, $count_strong, $count_weak, $input_format, $pseudosequence, $thr_rank_S, $thr_rank_W, $rank_f, $aff_f, $filter_message, $rank_flag) = @_;
	# defining the names for the initial lines to print 
	my $pos_print = sprintf("%6s", "pos");
	my $pos_core_print = sprintf("%4s", "Pos");
	my $allele_print = sprintf("%22s", "Allele");
	my $peptide_print = sprintf("%17s", "peptide");
	my $core_print = sprintf("%12s", "Core");
	my $identity_print = sprintf("%11s", "Identity");
	my $score_print = sprintf("%13s", "1-log50k(aff)");
	my $affinity_print = sprintf ("%13s", "Affinity(nM)");
	my $rank_print = sprintf("%6s", "\%Rank");
	my $level_print = " BindingLevel";
	my $new_thr_aff_S = sprintf("%.3f", $thr_aff_S);
	my $new_thr_aff_W = sprintf("%.3f", $thr_aff_W);
	my $new_aff_f = "";
	if (defined $aff_f) {
		$new_aff_f = sprintf("%.3f", $aff_f);
	}
	my $expected_affinity_print = "";
	if ($flag_expected_aff == 1) {
		$expected_affinity_print = sprintf("%11s", "Exp_Binding");
	}
	if ($print_count == 0) {
		print "# NetMHCIIpan version 3.0\n\n" , 
      		"# Input is in $input_format format\n\n";
		
		if ($inptype ==0) {
			print "# Peptide length $length\n\n";
		}      		
		if ($rank_flag == 0) {
			print "# Use user MHC pseudo sequence $pseudosequence USER_DEF\n\n";
		}
		print  "# Threshold for Strong binding peptides (IC50)\t$new_thr_aff_S nM\n",
      		"# Threshold for Weak binding peptides (IC50)\t$new_thr_aff_W nM\n\n";
		if ($rank_flag == 1) {
			print "# Threshold for Strong binding peptides (\%Rank)\t$thr_rank_S%\n",
			"# Threshold for Weak binding peptides (\%Rank)\t$thr_rank_W%\n";
		}
		if ($filter == 1) {
			print "\n# Threshold for filtering output (\%Rank)\t$rank_f%\n",
			"# Threshold for filtering output (IC50)\t\t$aff_f nM\n";
		}
		if (defined $filter_message) {
			print "\n\n$filter_message\n";
		}
      		
		$print_count++;
	}
		
	# Printing the results for each allele
	print "\n# Allele: $allele\n";	
	if ($rank_flag == 1) {
		print "-------------------------------------------------------------------------------------------------------------------------------------------\n",
		"$pos_print $allele_print $peptide_print $identity_print $pos_core_print $core_print $score_print $affinity_print $rank_print $expected_affinity_print $level_print\n",
		"-------------------------------------------------------------------------------------------------------------------------------------------\n";	
	}
	else {
		print "-------------------------------------------------------------------------------------------------------------------------------------------\n",
		"$pos_print $allele_print $peptide_print $identity_print $pos_core_print $core_print $score_print $affinity_print $expected_affinity_print $level_print\n",
		"-------------------------------------------------------------------------------------------------------------------------------------------\n";
	}	
	print $FINAL_RESULT;
	print "-------------------------------------------------------------------------------------------------------------------------------------------\n",
	      "Number of strong binders: $count_strong Number of weak binders: $count_weak\n",
	      "-------------------------------------------------------------------------------------------------------------------------------------------\n";
}

sub FindBestCore {
	my ($FINAL_RESULT, $flag_expected_aff, $rank_flag) = @_;
	my $RESULT = "";
	my @result_lines = split("\n", $FINAL_RESULT);
	my @core_pos = ();
	my @cores = ();
	my @scores =();
	my %same_core_lines = ();
	my @unique_lines = ();
	my $count_strong = 0;
	my $count_weak = 0;
	foreach my $line (@result_lines) {
		my @tmp = split( " ", $line);
		my $core_pos = $tmp[4];
		my $core = $tmp[5];
		my $score = $tmp[6];
		push (@core_pos, $core_pos);
		push (@cores, $core);
		push (@scores, $score);	
	}
	my $no = 0;
	for (my $n = 0; $n <= $#result_lines; $n++) {
		$no++;
		if ($cores[$n] eq $cores[$n-1] and $core_pos[$n] == ($core_pos[$n-1] - 1))  {
			$same_core_lines{$result_lines[$n]} = $scores[$n];
		}
		if ($n == 0 and $#result_lines != 0 and $cores[$n] ne $cores[$n+1]) {
			push (@unique_lines, $result_lines[$n]);
		}
		elsif ($n == 0 and $#result_lines == 0){
			push (@unique_lines, $result_lines[$n]);
		}
		
		if ($no == 11 or ($cores[$n] ne $cores[$n-1]) or $n == $#result_lines) {
			my $first = 0;
			foreach my $key (sort { $same_core_lines{$b} <=> $same_core_lines{$a} } keys %same_core_lines) {
				if ($first == 0 and exists $same_core_lines{$key}) {
					push (@unique_lines, $key);
					$first++;
				}
			}
			%same_core_lines = ();
			$no = 0;				
		}
		if ($n != $#result_lines and ($cores[$n] ne $cores[$n-1] and $cores[$n] ne $cores[$n+1]) and $n != 0) {
			push (@unique_lines, $result_lines[$n]);
		}
		
		if ($n == $#result_lines and $cores[$n] ne $cores[$n-1]) {
			push (@unique_lines, $result_lines[$n]);
		}
	}
	foreach my $unique_line (@unique_lines) {
		my $level = "";
		my @tmp = split (" ", $unique_line);
			if ($rank_flag == 1) {
				if (defined $tmp[10]) {
					$level = $tmp[10];
				}
				elsif (defined $tmp[9] and !defined $tmp[10] and $flag_expected_aff != 1) {
					$level = $tmp[9];
				}	
			}
			else {
				if (defined $tmp[9]) {
					$level = $tmp[9];
				}
				elsif (defined $tmp[8] and !defined $tmp[9] and $flag_expected_aff != 1) {
					$level = $tmp[8];
				}				
			}
		if ($level eq "<=SB") {
			$count_strong++;			
		}
		elsif ($level eq "<=WB") {
			$count_weak++;			
		}		
	}
	$RESULT = join("\n", @unique_lines);
	$RESULT .= "\n";
	return ($RESULT, $count_strong, $count_weak);
}

# program usage	
sub usage {
	print "\nUsage: ./NetMHCIIpan-3.0.pl [-h] [args] -f [fastafile/peptidefile]\n";
	print "Command line options:\n\n";
	printf ("%-16s\t%-32s\t%-12s\n",  "PARAMETER", "DEFAULT VALUE", "DESCRIPTION");
	printf ("%-16s\t%-33s\t%-12s\n",  "[-rdir dirname]", "$NETMHCIIpan", "Home directory for NetMHCIIpan");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-tdir dirname]", "$TMPDIR/tmp_\$\$", "Temporary directory");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-a name]", "DRB1_0101", "HLA allele");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-choose]", " ", "Choose alpha and beta chains separately");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-cha name]", "", "Alpha chain name");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-chb name]", "", "Beta chain name");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-affS float]", "50.000", "Threshold for strong binders (IC50)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-affW float]", "500.000", "Threshold for weak binders (IC50)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-rankS float]", "0.5", "Threshold for strong binders (\%Rank)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-rankW float]", "2", "Threshold for weak binders (\%Rank)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-filter int]", "0", "Filter output [1]");	
	printf ("%-16s\t%-32s\t%-12s\n",  "[-affF float]", "500", "Threshold for filtering output (IC50). Used only when filter is 1");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-rankF float]", "2", "Threshold for filtering output (\%Rank). Used only when filter is 1. Ignored if -hlaseq option was used.");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-hlaseq filename]", " ", "File with full length MHC beta chain sequence (used alone for HLA-DR)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-hlaseqA filename]", " ", "File with full length MHC alpha chain sequence (used with -hlaseq option)");	
	printf ("%-16s\t%-32s\t%-12s\n",  "[-inptype int]", "0", "Input type [0] FASTA [1] Peptide");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-length int_array]", "15", "Peptide length. Necessary for FASTA input only.");	
	printf ("%-16s\t%-32s\t%-12s\n",  "[-s]", "0", "Sort output on descending affinity");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-u]", "0", "Print unique binding core only");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-fast]", "0", "Use fast mode (10 best networks)");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-f filename]", " ", "File with the input data");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-xls]", "0", "Save output into xls file");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-xlsfile filename]", "NetMHCIIpan_out.xls", "File name for xls output");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-dirty]", "0", "Dirty mode, leave tmp dir+files");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-w]", "0", "w option for webface");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-list]", "0", "Print the list of possible alleles and exit");
	printf ("%-16s\t%-32s\t%-12s\n",  "[-h]", "0", "Print this message and exit");
}
