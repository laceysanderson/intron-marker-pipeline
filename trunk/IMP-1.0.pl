#!/usr/bin/perl
use Getopt::Long;
require File::Spec;

$RepeatMasker_Path = "/usr/local/share/applications/BinfoTools/RepeatMasker/RepeatMasker";
$BLAT_Path = "/usr/local/share/applications/BinfoTools/BLAT/blat";
$GeneSeqer_Path = "/usr/local/share/applications/BinfoTools/GeneSeqer/bin";
$GeneSeqer_Path_makearray = join("/", $GeneSeqer_Path, "MakeArray");
$GeneSeqer_Path_main = join("/", $GeneSeqer_Path, "GeneSeqer");
$Primer3_Path = "/usr/local/share/applications/BinfoTools/Primer3/src/primer3_core";
$extract_fasta_Path = "/usr/local/share/applications/BinfoTools/extract_fasta.pl";

########################################################################################################################################

#COPYRIGHT 2009 Lacey Sanderson, Perumal Vijayan, Kirstin Bett

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


########################################################################################################################################

#This program is designed to take two datasets: query sequences (a single Batch FASTA) & genome sequences (a single Batch FASTA) and
#1) Mask the query sequences using RepeatMasker
#2) Align the query ESTs against the genome sequences using BLAT or GeneSeqer
#3) Use the output from this process to determine the number of introns a given EST might have
#and their start and end locations
#4) Use this information to design primers (using Primer3) which flank each putative intron

#Synopsis:

#Usage:
#IMP [options] database query

#Where:
#Database and query are each in a single batch FASTA file.

#Options:
#       -r|noMask                  Does not run repeat masker. Thus BLAT will be run on query file directly.
#       -a|Alignment=s             Choose the program to be used to align the ESTs to the Genome where s = b for blat and g for geneseqer. Default: Geneseqer
#       -d|rmDash=s                Removes any dashes (-) from the query sequence and replaces them with a space (" "). Also ensures that the query is the proper sequence for input into the pipeline
#                                  s = q (process query only), t (process template only) or b (process both)
#		    -5extensionFwd=s		       Allows the user to enter s, which will be added to the 5' end of the forward primer
#		    -5extensionRvs=s		       Allows the user to enter s, which will be added to the 5' end of the reverse primer

#		    -autocurate=s			         Will pick the best primer pair per intron (P), best intron per gene model (I) and best gene model per EST (G). s = P, I, G or any combination

#    Following options only applicable with -a=b
#       -o|ooc                     Use only in combination with -b. Allows the generation of the ooc required for BLAT (5.ooc for cross-species DNA to DNA)
#       -q|useQuality              Adds quality information for designing primers. File containing the quality information
#                                  must be of the form queryfilename.qual

#    Following options only applicable with -a=g
#       -p|preProcess              Use with Geneseqer. Allows preprocessing of the query by running MakeArray distributed with Geneseqer
#       -s|species=s               Use with Geneseqer. Specifies splice site model to use where the default is Medicago. Other splice site models available include "human", "mouse", "rat", "chicken", "Drosophila", 
#                                  "Daphnia", "nematode", "yeast", "Aspergillus", "Arabidopsis", "maize", "rice", "generic"

#    Following options for configuring Primer3
#       -PoptSize=i                Optimum length (in bases) of a primer oligo. Primer3 will attempt to pick primers close to this length.
#       -PminSize=i                Minimum acceptable length of a primer.  Must be greater than 0 and less than or equal to PRIMER_MAX_SIZE.
#       -PmaxSize=i                Maximum acceptable length (in bases) of a primer.  Currently this parameter cannot be larger than 35.  This limit is governed by
#                                  maximum oligo size for which primer3's melting-temperature is valid.
#       -PoptTm=f                  Optimum melting temperature(Celsius) for a primer oligo. Primer3 will try to pick primers with melting temperatures are close to this temperature.
#       -PminTm=f                  Minimum acceptable melting temperature(Celsius) for a primeroligo.
#       -PmaxTm=f                  Maximum acceptable melting temperature(Celsius) for a primer oligo.
#       -PMask=i                   This option allows for intelligent design of primers in sequence in which masked regions (for example repeat-masked regions) are lower-cased.
#                                  A value of 1 directs primer3 to reject primers overlapping lowercase a base exactly at the 3' end; whereas, a value of 0 inactivates this feature.
#       -PminGC=f                  Minimum allowable percentage of Gs and Cs in any primer.

#######################################################################################################################################
#          MAIN PROGRAM
#######################################################################################################################################

##subroutine declarations
## all subroutines are at the bottom
sub print_synopsis();
sub print_header($);
sub get_primerset($$$$$);

sub blat_parser($$$$);
sub get_output_BLAT($$$$);

sub geneseqer_parser($$$);
sub get_output_GENESEQER($$$$);

sub parse_Primer3_out ($$);

sub rmDash($$$);
sub geneseqer_autocurate ($$$$);
sub blat_autocurate ($$$$);

my $devnull = File::Spec->devnull();

##PARSE THE COMMAND LINE OPTIONS AND CHECK FOR VALIDITY-------------------------------------------------
if (($#ARGV + 1) == 0) {
  &print_synopsis;
  exit;
}

##option defaults
%options = ('a' => 'g', 'Bt' => 'dnax', 'Bq' => 'dnax', 's' => 'Medicago', 'PoptSize' => 20, 'PminSize' => 18, 'PmaxSize' => 27, 'PoptTm' => 60.0, 'PminTm' => 57.0, 'PmaxTm' => 63.0, 'PMask' => 1, 'PminGC' => 20.0, 'd' => 'q');

##Processing of options, 
GetOptions(\%options, 'a|Alignment=s', 'o|ooc', 'r|noMask', 'p|preProcess', 'q|useQuality', 'Bt=s', 'Bq=s', 's|species=s', 'PoptSize=i', 'PminSize=i', 'PmaxSize=i', 'PoptTm=f', 'PminTm=f', 'PmaxTm=f', 'PMask=i', 'PminGC=f', 'd|rmDash=s', 'v|verbose', '5extensionFwd=s', '5extensionRvs=s', 'c|autocurate=s');

##get arguements
$database = $ARGV[0];
$queryPATH =  $ARGV[1];
$queryPATH =~ /(?:\.)*\/(?:\S*\/)*(\S*)/;
$query = $1;

print "The command line arguements are: \ndatabase-", $database, ".\nquery-", $queryPATH, "($query)", ".\nOptions are ";

while( my ($k, $v) = each %options ) {
    print "($k, $v), ";
}

print ".\n\n";

##check files exist
if (-e $database) {
	print "File: $database found...\n";
} else {
	print "ERROR! File: $database not found!\n";
	exit 1;
}
if (-e $queryPATH) {
	print "File: $queryPATH found...\n";
} else {
	print "ERROR: File: $queryPATH not found!\n";
	exit 2;
}

##Validate file contents (optional)----------------------------------------------------------------------------------
if (defined $options{d}) {
    ##remove the dashes from the sequence and ensure the rest is of a recognized format [GACTgactNn ]
    &rmDash($queryPATH, $database, \%options);
    if ($options{d} eq 'q') {$query = "$query.rmDash"; $queryPATH = "$queryPATH.rmDash";}
    elsif ($options{d} eq 't') {$database = "$database.rmDash";}
    elsif ($options{d} eq 'b') {$query = "$query.rmDash"; $queryPATH = "$queryPATH.rmDash"; $database = "$database.rmDash";} 
    elsif ($options{d} eq 'n') {}
    else { print "ERRPR: Unrecorgnized option $option{d} for -d\n"; exit;}
}

##Program Flow proper------------------------------------------------------------------------------------------------
##use BLAT for phase 2 of pipeline...................................................................................
if ($options{a} eq "b") {

  ##run blat without masking query
  if ($options{r} == 1) {
    
    ##make overused tile file
    if ($options{o} == 1) {
  	  print "\n\nGenerating the 5.ooc file used by BLAT...\n";
  	  $out = `$BLAT_Path -t=$options{Bt} $database -q=$options{Bq} $queryPATH  -qMask=lower -makeOoc=5.ooc ooc.out`;
  	  if ($options{v} ==1) {
        print "$out\n";
  	  }
    } #end of make overused tile file
    
	  print "\n\nAligning ESTs in $query with $database using BLAT...\n";
 
	  $out = `$BLAT_Path -t=$options{Bt} $database -q=$options{Bq} $queryPATH -ooc=5.ooc -noHead -qMask=lower $query.blat.out.psl`;
  	if ($options{v} ==1) {
  	  print "The output of Blat is: \n $out\n";
  	  print "BLAT completed with an exit code of $?.\n";
  	}

	  if ($? != 0) {
	    print "ERROR: BLAT completed with a non-zero exit code of $?\n";
	    exit (200+$?);
	  }
	  
  ##run blat after masking query
  } else {
	  print "\nMasking repeats in $query using Repeat Masker...\n";
	  
	  ##run RepeatMasker
	  system $RepeatMasker_Path, '-xsmall', '-nocut', $queryPATH;
	  $masked = "$queryPATH.masked";
	  if ($options{v} ==1) {
	    print "RepeatMasker completed with an exit code of $?\n";
	    print "The masked sequences is in $masked\n";
	  }

	  if ($? != 0) {
	    print "ERROR: RepeatMasker completed with a non-zero exit code of $?\n";
	    exit (100 + $?);
	  }

    ##make overused tile file
	  if ($options{o} == 1) {
	    print "\n\nGenerating the 5.ooc file used by BLAT...\n";
	    $out = `$BLAT_Path -t=$options{Bt} $database -q=$options{Bq} $masked  -qMask=lower -makeOoc=5.ooc ooc.out`;
	  }

     ##run BLAT
	   print "\n\nAligning ESTs in $masked with $database using BLAT...\n";
	   $out = `$BLAT_Path -t=$options{Bt} $database -q=$options{Bq} $masked -ooc=5.ooc -noHead -qMask=lower $query.blat.out.psl`;
	   if ($options{v} ==1) {
       print "The output of Blat is: \n $out\n";
       print "BLAT completed with an exit code of $?\n";
     }

	   if ($? != 0 ) {
       print "ERROR: BLAT completed with a non-zero exit code of $?";
       exit (200+$?);
	   }
  }##end of if else mask queries and run blat

  ##parse BLAT output
  $blatRef = &blat_parser ($query, $queryPATH, $database, \%options);
  @blatRecords = @{$blatRef};

  ##run primer3
  print "\n\nDesigning primers using primer3...\n";
  $out = `$Primer3_Path < $query.primer3.in.boulderIO > $query.primer3.out.boulderIO | > $query.primer3.error.boulderIO`;

  $ref = &parse_Primer3_out($query, $database);
  @primer3_out = @{$ref};
  
  ##Compile full output
  &get_output_BLAT($query, $queryPATH, $blatRef, $ref);

##use Geneseqer for phase 2 of pipeline...................................................................................
} elsif ($options{a} eq "g") {

  ##Run GeneSeqer w/out masking query
  if ($options{r} == 1) {

    ##run Geneseqer after preprocessing
    if ($options{p} == 1) {
      ##preprocessing
      print "\n\nPreprocessing $queryPATH for Geneseqer Alignments...\n";
      $out = `$GeneSeqer_Path_makearray $queryPATH`;
	    if ($options{v} ==1) {
		    print $out;
	    }

      ##GeneSeqer
      $out = "";
	    print "Aligning $queryPATH with $database using Geneseqer...\n";
		  $out = `$GeneSeqer_Path_main -s $options{s} -D $queryPATH -R -L $database -o $query.GeneSequer.out 2> $devnull`;
	    if ($options{v} ==1) {
		    print $out;
	    }

    ##run GeneSeqer w/out preprocessing
    } else {

      $out = "";
      print "\n\nAligning $queryPATH with $database using Geneseqer...\n";
      $out = `$GeneSeqer_Path_main -s $options{s} -E $queryPATH -R -L $database -o $query.GeneSequer.out 2> $devnull`;
      if ($options{v} ==1) {
		    print $out;
	    }

    }##end of if else preprocessing

  ##Run GeneSeqer after masking query
  } else {

    ##run RepeatMasker
    print "\nMasking repeats in $query using Repeat Masker...\n";
    system $RepeatMasker_Path, '-xsmall', '-nocut', $queryPATH;
    $masked = "$queryPATH.masked";
    
    if ($options{v} ==1) {
      print "RepeatMasker completed with an exit code of $?\n";
	    print "The masked sequences is in $masked\n";
	  }

    if ($? != 0) {
      print "ERROR: RepeatMasker completed with a non-zero exit code of $?\n";
      exit (100 + $?);
    }

	  ##run Geneseqer w/ preprocessing
	  if ($options{p} == 1) {
      print "\n\nPreprocessing $masked for Geneseqer Alignments...\n";
	    $out = `$GeneSeqer_Path_makearray $masked`;
	    
	    if ($options{v} ==1) {
		    print $out;
	    }

      $out = "";
      print "Aligning $masked with $database using Geneseqer...\n";
	    $out = `$GeneSeqer_Path_main -s $options{s} -D $masked -R -L $database -o $query.GeneSequer.out 2> $devnull`;
	    if ($options{v} ==1) {
			  print $out;
	    }

    ##run GeneSeqer w/out preprocessing
	  } else {

	    $out = "";
      print "\n\nAligning $masked with $database using Geneseqer...\n";
      $out = `$GeneSeqer_Path_main -s $options{s} -E $masked -R -L $database -o $query.GeneSequer.out 2> $devnull`;
      if ($options{v} ==1) {
		    print $out;
	    }
	    
	  } ##end of if else preprocessing
  } ##end of if else mask query

  ##parse GeneSeqer output
  $geneseqerRef = &geneseqer_parser($query, $queryPATH, \%options);
  @Genes = @{ $geneseqerRef };

  ##run primer3
  print "\n\nDesigning primers using primer3...\n";
  $out = `$Primer3_Path < $query.primer3.in.boulderIO > $query.primer3.out.boulderIO | > $query.primer3.error.boulderIO`;
	
  $ref = &parse_Primer3_out($query, $database);
  @primer3_out = @{$ref};

  ##generate complete output
  $refOUT = &get_output_GENESEQER($query, $queryPATH, $geneseqerRef, $ref);

} ##end of if else blat or geneseqer

#######################################################################################################################################
#          SUBROUTINES
#######################################################################################################################################

## Prints the synopsis to standard output where the synopsis is above and includes usage and option descriptions
## Returns: nothing
## Arguements: none
sub print_synopsis () {

  print "Usage:\n";
  print "IMP [options] database query\n";
  print "Where:\n";
  print "Database and query are each in a single batch FASTA file\n\n";

  print "Options:\n";
  print "\t-r|noMask\t\tDoes not run repeat masker. Thus BLAT will be run on query file directly.\n";
  print "\t-a|Alignment=s\t\tChoose the program to be used to align the ESTs to the Genome where s = b for blat and g for geneseqer. Default: Geneseqer.\n";
  print "\t-d|rmDash=s\t\tRemoves any dashes (-) from the query sequence and replaces them with a space (\" \"). Also ensures that the query is the proper sequence for input into the pipeline. s = q (process query only), t (process template only) or b (process both).\n";
  print "  Following options only applicable with -a=b.\n";
  print "\t-o|ooc\t\t\tUse only in combination with -b. Allows the generation of the ooc required for BLAT (5.ooc for cross-species DNA to DNA).\n";
  print "\t-q|useQuality\t\tAdds quality information for designing primers. File containing the quality information must be of the form queryfilename.qual.\n";
  print "  Following options only applicable with -a=g.\n";
  print "\t-p|preProcess\t\tUse with Geneseqer. Allows preprocessing of the query by running MakeArray distributed with Geneseqer.\n";
  print "\t-s|species=s\t\tUse with Geneseqer. Specifies splice site model to use where the default is Medicago. Other splice site models available include human, mouse, rat, chicken, Drosophila, Daphnia, nematode, yeast, Aspergillus, Arabidopsis, maize, rice, generic.\n";

  print "  Following options for configuring Primer3.\n";
  print "\t-PoptSize=i\t\tOptimum length (in bases) of a primer oligo. Primer3 will attempt to pick primers close to this length.\n";
  print "\t-PminSize=i\t\tMinimum acceptable length of a primer.  Must be greater than 0 and less than or equal to PRIMER_MAX_SIZE.\n";
  print "\t-PmaxSize=i\t\tMaximum acceptable length (in bases) of a primer.  Currently this parameter cannot be larger than 35.  This limit is governed by maximum oligo size for which primer3's melting-temperature is valid.\n";
  print "\t-PoptTm=f\t\tOptimum melting temperature(Celsius) for a primer oligo. Primer3 will try to pick primers with melting temperatures are close to this temperature.\n";
  print "\t-PminTm=f\t\tMinimum acceptable melting temperature(Celsius) for a primeroligo.\n";
  print "\t-PmaxTm=f\t\tMaximum acceptable melting temperature(Celsius) for a primer oligo.\n";
  print "\t-PMask=i\t\tThis option allows for intelligent design of primers in sequence in which masked regions (for example repeat-masked regions) are lower-cased. A value of 1 directs primer3 to reject primers overlapping lowercase a base exactly at the 3' end; whereas, a value of 0 inactivates this feature.\n";
  print "\t-PminGC=f\t\tMinimum allowable percentage of Gs and Cs in any primer.\n";

} ##Sub End of print_synopsis

## compiles the header for full output
## returns a string containing the full header
## arguement is a string of either "blat" or "geneseqer depending which header is needed"
sub print_header($) {
  $mode = $_[0]; ##either blat or geneseqer
  if ($mode eq "blat") {
    $header = "Query\t\t\t\t\t\tTemplate\t\t\t\tIntron\t\t\t\t\t\t\t\tPrimers\n"."Name\tNo. Repeats\tStrand\tSize\tStart\tEnd\t|Name\tSize\tStart\tEnd\t|Match/Mismatch\tNo.Introns\tCurr. Intron\tStart\tEnd\tLength\tTrg Length\tIncluded Reg\t|\tLeft\t\t\t\t\tRight\t\t\t\t\tProduct\t\tDiagnostic\n"."\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tPrimerset\tStart\tLength\tSequence\tTm\tGC\t|Start\tLength\tSequence\tTm\tGC\t|EST Product Size\tGenomic Product Size\t|Pair Penalty\tCompl-Any\tCompl-End\tF/R Penalty\tF/R Self-Any\tF/R Self-End\t F/R End Stability\n";
  } elsif ($mode eq "geneseqer") {
    $header = "Query\t\t\t\t\tTemplate\t\t\t\tIntron\t\t\t\t\t\t\t\t\t\tPrimers\n"."Name\tNo. Repeats\tStrand\tStart\tEnd\t|Name\tStrand\tStart\tEnd\t|Similarity\tNo.Introns\tCurr. Intron\tStart\tEnd\tLength\tTrg Length\tIncluded Reg\tD,A P-value\tD,A Score\t|Left\t\t\t\t\t\tRight\t\t\t\t\tProduct\t\tDiagnostic\n"."\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tPrimerSet\tStart\tLength\tSequence\tTm\tGC\t|Start\tLength\tSequence\tTm\tGC\t|EST Product Size\tGenomic Product Size\t|Pair Penalty\tCompl-Any\tCompl-End\tF/R Penalty\tF/R Self-Any\tF/R Self-End\t F/R End Stability\n";
  } else {
    print "ERROR: Couldn't find header type.\n";
  }
   return $header;
}

## gets all the information for a given primerset (reduces repition for 5 primersets)
## returns an array containing 1 value per column pertaining to primers in final output
## arguements supplied include the current primerset (1-5), the genomic product size, and three arrays containing details pertaining to the left, right and pair values of a given primer3 record
sub get_primerset($$$$$) {
  $curr_primerset = $_[0]; 
  $Gnm_prod_size = $_[1]; 
  @recP3_left = @{$_[2]};
  @recP3_right = @{$_[3]};
  @recP3_pair = @{$_[4]};
	
	#Adds tail if one is provided
  if (exists $options{'5extensionFwd'}) {
	  $Fwd = "$options{'5extensionFwd'}$recP3_left[4]";
  } else {
	  $Fwd = $recP3_left[4];
  }
  if (exists $options{'5extensionRvs'}) {
    $Rvs = "$options{'5extensionRvs'}$recP3_right[4]";
  } else {
    $Rvs = $recP3_right[4];
  }
  
	@primerset = ($curr_primerset, $recP3_left[7], $recP3_left[5], $Fwd, $recP3_left[2], $recP3_left[3], $recP3_right[7], $recP3_right[5], $Rvs, $recP3_right[2], $recP3_right[3], $recP3_pair[2], $Gnm_prod_size, $recP3_pair[5], $recP3_pair[3], $recP3_pair[4], "$recP3_left[8]/$recP3_right[8]", "$recP3_left[9]/$recP3_right[9]", "$recP3_left[10]/$recP3_right[10]", "$recP3_left[11]/$recP3_right[11]");
  
  return \@primerset;
} ##Sub end get_primerset

##-----GeneSeqer related subroutines-----------------------------------------------------------------------------------------------

## Parses the GeneSeqer output and creates Primer3 input
## Returns a Genes hash descrbing each Gene Model
## Arguements supplied include $query, $queryPATH, %options
sub geneseqer_parser ($$$) {
  $query = $_[0];
  $queryPATH = $_[1];
  %options = %{ $_[2] };

  ##parse GeneSeqer ouput in file GeneSequer.out
  open(GSQR, "$query.GeneSequer.out")  || die "Can't open $query.GeneSequer.out: $!";

  ##remove header
  $line = <GSQR>;
  while ($line) {
    if ($line =~ /_+[^A-Z][^a-z]/) {
      print "Parsing $query.GeneSeqer.out...\n";
      last;
    }
    $line = <GSQR>;
  }

  ##separate into records delinated by ----
  $line = <GSQR>;
  $i=0;

  ##REMOVE RECORDS NOT CONTAINING ALIGNMENTS
  $a = 0; #for progess ... so that you only get a dot every 10 records processed
  $TotalUnmatched = 0;
  $TotalMatched = 0;
  print "Removing the records that do not contain alignments from the Geneseqer output";
  while ($line) {                                                                   #for each record in the output
    if ($a % 10 == 0) { print ".";}
	  $a ++;

	  @rec = ();
	  while ($line) {                                                               #for each line in a given record
	    if ($line =~ /^_+$/) {
		    $line = <GSQR>;
		    last;
	    } else {
		    ##add to records 2D array
		    @rec = (@rec, $line);
		    $line = <GSQR>;
	    }
	  }##end of GeneSeqer record (delineated by ----)

    ##check the last read record contains an alignment
	  $noMatch = FALSE;
	  foreach $ln (@rec) {
	    ##removes all records with No significant EST matches were found.
	    if ($ln =~ /^.*No.*$/ ) {
		    $noMatch = TRUE;
		    $TotalUnmatched = $TotalUnmatched + 1;
		    last;
	    }
	    ##removes all records with Total number of EST alignments reported: 0
	    if ($ln =~ /^.*reported: 0.*$/) {
		    $TotalUnmatched = $TotalUnmatched + 1;
		    $noMatch = TRUE;
		    last;
	    }
	    ##removes all 3-phase translations of preceeding alignment from records
	    if ($ln =~ /^.*translation.*$/) {
		    $noMatch = TRUE;
		    last;
	    }
	  } ##end of foreach Geneseqer record

    ##checks that first line is not empty (for true matches with alignments 1st line will contain Sequence ##: etc.)
	  if ($rec[1] =~ /^.*finished.*$/ ) {
	    $noMatch = TRUE;
	  }

	  if ($noMatch eq FALSE) {
	    $gsqrRecords[$#gsqrRecords+1] = [ @rec ];
	    $TotalMatched ++;
	  }
  } ##while there are still lines in GeneSeqer file

  print "\n$TotalMatched($#gsqrRecords) query FASTA records were aligned ($TotalUnmatched were not).\n";
  
  ##GENERATES PRIMER3 INPUT WHILE PARSING GENESEQER ALIGNMENT INFORMATION
  open(P3IN, ">$query.primer3.in.boulderIO");

  ##put global parameters into the first record of the Primer3 input file
  print P3IN "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0\n";
  print P3IN "PRIMER_PRODUCT_SIZE_RANGE=37-1000\n";
  print P3IN "PRIMER_TM_SANTALUCIA=1\n";
  print P3IN "PRIMER_OPT_SIZE=$options{'PoptSize'}\n";
  print P3IN "PRIMER_MIN_SIZE=$options{'PminSize'}\n";
  print P3IN "PRIMER_MAX_SIZE=$options{'PmaxSize'}\n";
  print P3IN "PRIMER_OPT_TM=$options{'PoptTm'}\n";
  print P3IN "PRIMER_MIN_TM=$options{'PminTm'}\n";
  print P3IN "PRIMER_MAX_TM=$options{'PmaxTm'}\n";
  print P3IN "PRIMER_LOWERCASE_MASKING=$options{'PMask'}\n";
  print P3IN "PRIMER_MIN_GC=$options{'PminGC'}\n";

  print "\n\nProcessing information from the alignments & generating Primer3 input file";

  ## Extract Sequence, Exon/Intron & alignment information for each alignment
  @Genes = ();
  $a=0; #so that a . only gets printed out for every 10 records

  for($i=0;$i<=$#gsqrRecords;$i++) {                                                              #loops through each record in the output
    if ($a % 10 == 0) { print ".";}
    $a ++;

	  @ei = @{$gsqrRecords[$i]};                                                                  #gets all lines of a given record & store in array
	  for($j=0;$j<=$#ei;$j++) {                                                                   #for each line of a given record does the following
	    @g = ();                                                                      #holds the arrays describing a given gene prediction

	    if ($ei[$j] =~ /^.*EST sequence.*$/) {
		    $PROCESSING = "SEQ";
		    $j = $j + 2;
		    $seq = "";
		    ##get the sequence and remove line numbers
		    $CONT = TRUE;
		    while ($CONT eq TRUE) {                                                           #builds up the EST sequence from the record
		      $GI = $ei[$j];
		      if ($GI =~ /^\s*\d+\s*([GACTgactNn]+)\s*([GACTgactNn]*)\s*([GACTgactNn]*)\s*([GACTgactNn]*)\s*([GACTgactNn]*)\s*([GACTgactNn]*).*$/) {
			      $seq = "$seq$1$2$3$4$5$6";
	        } else {
		        $CONT = FALSE;
	        }
		      $j++;
		    }##end of sequence build up loop
		    @g[$#g+1] = $seq;                                                   #adds the EST sequence as the first element of the gene

        #print "\nProcessing the Exon/Intron information...";
		    $GI = $ei[$j];
		    $CONT = TRUE;
		    while ( $CONT eq TRUE) {                                                           #builds up the exon intron information (each is its own array in g)
		      if ($GI =~ /^.*Exon.*$/ ) {
			      $PROCESSING = "EXON";
			      #Extract info from this line
			      if ($GI =~ /^.*Exon\s*(\S+)\s*(\S+)\s*(\S+)\s*\(\s*(\S*)\s*n\)\;\s*\S*\s*(\S*)\s*(\S*)\s*\(\s*\S*\s*n\)\;\s*score:\s*(\S*).*$/) {
			        @exon = (E, $1, $2, $3, $4, $5, $6, $7);
			        @g[$#g+1] = [ @exon ];
			      } else {
			        print "ERROR: Exon description doesn't match expected format!\n";
			      }
		      } elsif ($GI =~ /^.*Intron.*$/) {
			      $PROCESSING = "INTRON";
			      if ($GI =~ /^.*Intron\s*(\S+)\s*(\S+)\s*(\S+)\s*\(\s*(\S*)\s*n\)\;\s*Pd\:\s*(\S*)\s*\(s\:\s*(\S*)\)\,\s*Pa\:\s*(\S*)\s*\(s\:\s*(\S*)\).*$/) {
			        @intron = (I, $1, $2, $3, $4, $5, $6, $7, $8);
			        @g[$#g+1] = [ @intron ];
			      } else {
			        print "ERROR: Intron description doesn't match expected format!\n";
			      }
		      } else {
			      if ($PROCESSING eq "SEQ") {} else {$CONT = FALSE};
		      }##end of if exon, intron, other
		      $j ++;
		      $GI = $ei[$j];
		    }#end of building up exon/intron information

        #print "\nProcessing alignment information.\n";
		    $GI= $ei[$j];
		    if ($GI =~ /^.*MATCH\s*(\S*)([+|-])\s*(\S*)([+|-])\s*(\S*)\s*(\S*)\s*(\S*).*$/) {    #collects info about the EST and BAC and alignment
		      @info = ($#Genes+1, BAC,$1,$2, EST, $3, $4, SCORES, $5, $6, $7);            #NOTE: $#Genes+1 is the index in @Genes this gene description will be stored at & is unique
		      @g[$#g+1] = [ @info ];
		    }
		    
		    $Genes[$#Genes+1] = [ @g ];                                                          #adds g (seq, exons/introns) to the set of genes

        #print "Generating Primer3 input";
		    for($k=1; $k<=$#g;$k++) {
		      if ($g[$k][0] eq "I") {
      			print P3IN "PRIMER_SEQUENCE_ID=", $g[$#g][2], "_", $g[$#g][3], "_", $g[$#g][5], "_", $g[$#g][6], "_", $g[$#g][0], "_", $g[$k][1], "\n";
      			print P3IN "SEQUENCE=", $g[0], "\n";
      			print P3IN "TARGET=$g[$k-1][6],2\n";  #end on EST of exon before the current intron

      			$inclRegS = $g[$k-1][6] - $g[$k-1][4];
      			if ($inclRegS < 0) {$inclRegS=1;}
      			$e1 = $g[$k-1][4];
      			$e2 =$g[$k+1][4];
      			if ($e1 > 100) {$e1 = 100;}
      			if ($e2 > 100) {$e2 = 100;}			
      			$inclRegLen = $e1 + $e2;
    			
  			    print P3IN "INCLUDED_REGION=$inclRegS,$inclRegLen\n";
  			    print P3IN "=\n";
		      } ##end of if I
        } ##end of generating primer3input for
        
      } #end of if ($ei[$j] =~ /^.*EST sequence.*$/)
    } #end of for each line of a given record
  } #end of for each record

  print "\n";

  close GSQR;
  close P3IN;

  print "$#Genes geneseqer alignments were recorded.\n";

  #returns a reference for the root array containing all information about the predicted genes from the alignments
  return \@Genes;

}##Sub End of geneseqer_parser

## Compiles all output (before autocuration) from RepeatMasker, GeneSeqer & Primer3
## Returns an array containing each record printed to out.IMP for use in curation
## Arguements supplied include $query, $queryPATH, the reference to the Genes array & the reference to the primer3 output
sub get_output_GENESEQER($$$$) {
  $query = $_[0];
  $queryPATH = $_[1];
  @Genes = @{ $_[2] };
  @primer3_out = @{ $_[3] };

  print "\nGenerating program output...\n";

  ##open output files
  open(PO, ">$query.primerOrder.IMP");
  open(FO, ">$query.out.IMP");

  ##get repeat masker information
  open(RM, "$queryPATH.out");
  <RM>; <RM>; <RM>; ##remove header
  %Smasked = ();
  while (<RM>) {
  	$_ =~ /^\s*\S*\s*\S*\s*\S*\s*\S*\s*(\S*).*$/;
  	if (exists $Smasked{$1}) {
  		$val = $Smasked{$1};
  		$Smasked{$1} = $val ++;
  	} else {
  		$Smasked{$1} = 1;
  	}
  }

	##array containing all output record arrays
	@output = ();
	
  ##print headers
  print PO "Name\t%GC\tTm\tLength\tSequence\tPrimerSet\n";

  print FO &print_header("geneseqer");
  
  for ($j=0; $j<= $#primer3_out; $j++) {
		@recP3 = @{ $primer3_out[$j] };
  	@recP3_left = @{$recP3[1]};
  	$lp_seq = $recP3_left[4];
  	
    ##if there is at least 1 primer
  	if ($lp_seq =~ /[GACTgactNn]+/) {
  	  
  	  ##get geneseqer portion
      $seq_id = $recP3[0];
      @Gnsqr_id = split('_', $seq_id);
      $index = $Gnsqr_id[$#Gnsqr_id-2];
      
      @recG = @{ $Genes[$index]  };
      $curr_intron = @Gnsqr_id[$#Gnsqr_id-1]; 	  
  	  
      #Ex: @rec->seq|E|I|E|I|E|I|E|Info = ((8+1) -3) /2 =6/2=3
      $noIntrons = (($#recG+1) - 3)/2;
      ##start/end of query look at start of 1st exon and end of last exon
      $ESTstart = @{$recG[1]}[5];
      $ESTend = @{$recG[$#recG-1]}[6];
      $GnmStart = @{$recG[1]}[2];
      $GnmEnd = @{$recG[$#recG-1]}[3];

      @infoG = @{$recG[$#recG]};
      @infoP = @{$recP3[$#recP3]};
      if (exists $Smasked{$infoG[5]}) {$rmNo = $Smasked{$infoG[5]};} else {$rmNo = 0;}
      ##Ex: @recG->seq|E|I|E|I|E|I|E|Info = 2*2=4
      $index = $curr_intron * 2;
      @intron = @{$recG[$index]};
   	  
  	  @geneseqerOut = ($infoG[5], $rmNo, $infoG[6], $ESTstart, $ESTend, $infoG[2], $infoG[3], $GnmStart, $GnmEnd, $infoG[8]); ##the Geneseqer Query/Template portion
  	  @geneseqerOut = (@geneseqerOut, $noIntrons, $intron[1], $intron[2], $intron[3], $intron[4], 2, $infoP[6], "$intron[5],$intron[7]", "$intron[6],$intron[8]"); ##the Geneseqer intron portion
	    
      ##For each primerset
      for($i=0; $i<5; $i++) {
        if ($lp_seq =~ /[GACTgactNn]+/) {
          $iS = ($i * 3) + 1;
          @recP3_left = @{$recP3[$iS]};
        	@recP3_right = @{$recP3[$iS+1]};
        	@recP3_pair = @{$recP3[$iS+2]};
      	
        	$lp_seq = $recP3_left[4];
    	
          $Gnm_prod_size = 0;
    	    $Gnm_prod_size = $recP3_pair[2] + $recB[23];
	    
          $ref = &get_primerset($i+1, $Gnm_prod_size, \@recP3_left, \@recP3_right, \@recP3_pair);
          @primerRec = @{$ref};
      
    	    print FO join(" \t ", @geneseqerOut), " \t "; ##prints the blat portion
    	    print FO join(" \t ", @primerRec), "\n";  ##prints the primer portion
    
          ##prints to the simple primer output
    	    print PO $infoG[5], "_", $primerRec[0],"_Fwd\t ", $primerRec[5], " \t ", $primerRec[4], " \t ", $primerRec[2], " \t ", $primerRec[3], " \t ", $primerRec[0], "\n";
    	    print PO $infoG[5], "_", $primerRec[0],"_Rvs\t ", $primerRec[10], " \t ", $primerRec[9], " \t ", $primerRec[7], " \t ", $primerRec[8], " \t ", $primerRec[0], "\n";
	    
    	    #adds into array
    	    @record = (@geneseqerOut, @primerRec, $j, $index); #last two are primer3 id and blat id
    	    $output[$#output+1] = [ @record ];
    	  } ##end of if sequence for current primerset
      } ##end of for each primerset

	  } ##end of if there is @ least 1 primerset
	  
	} ##End of for each primer3 record
  
  close PO;
  close FO;
  close RM;
	
	if ($options{c} =~ /G/) {@opt = (1);} else {@opt = (0);}
	if ($options{c} =~ /I/) {@opt = (@opt, 1);} else {@opt = (@opt, 0);}
	if ($options{c} =~ /P/) {@opt = (@opt, 1);} else {@opt = (@opt, 0);}
	&geneseqer_autocurate($query, $queryPATH, \@output, \@opt);
	
}##Sub End of get_output_GENESEQER

##-----Blat related subroutines-----------------------------------------------------------------------------------------------

## Parses output from BLAT and generated Primer3 input file
## Returns an array containing each blat record
## Arguements supplied include $query, $queryPATH, $database, $options
sub blat_parser ($$$$) {
  $query = $_[0]; 
  $queryPATH = $_[1];
  $database = $_[2];
  %options = %{ $_[3]};

  ##parse blat output
  print "Generating primer3 input file...\n";
  open(BOUT, "$query.blat.out.psl") || die "Can't open $query.blat.out.psl: $!";
  open(P3IN, ">$query.primer3.in.boulderIO") || die "Can't open $query.primer3.in.boulderIO: $!";

  ##also keeps each line of the blat output in a 2D array @blat

  ##put global parameters into the first record of the Primer3 input file
  print P3IN "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0\n";
  print P3IN "PRIMER_PRODUCT_SIZE_RANGE=37-1000\n";
  print P3IN "PRIMER_TM_SANTALUCIA=1\n";
  print P3IN "PRIMER_LOWERCASE_MASKING=1\n";
  print P3IN "PRIMER_OPT_SIZE=$options{'PoptSize'}\n";
  print P3IN "PRIMER_MIN_SIZE=$options{'PminSize'}\n";
  print P3IN "PRIMER_MAX_SIZE=$options{'PmaxSize'}\n";
  print P3IN "PRIMER_OPT_TM=$options{'PoptTm'}\n";
  print P3IN "PRIMER_MIN_TM=$options{'PminTm'}\n";
  print P3IN "PRIMER_MAX_TM=$options{'PmaxTm'}\n";
  print P3IN "PRIMER_LOWERCASE_MASKING=$options{'PMask'}\n";
  print P3IN "PRIMER_MIN_GC=$options{'PminGC'}\n";

  @blat = ();
  $unqBLid = 100;
  while (<BOUT>) {
    chomp;
    @blatOutput = split(/\t/, $_);
    
    ##processes the line
    @qSTART = split(',', $blatOutput[19]);
    @blkSizes = split(',', $blatOutput[18]);
    @tSTART = split(',', $blatOutput[20]);
    $numIntrons = ($#qSTART+1) -1;

    for ($i=0; $i<($numIntrons-1); $i++) { ##for each intron
      $iStart = $qSTART[$i] + $blkSizes[$i];
      $SpliceSiteLen = $qSTART[($i+1)] - $iStart;
      $iLength = $tSTART[($i+1)] - ($tSTART[$i] + $blkSizes[$i]);
      $includeRegStart = $qSTART[$i];
      $includeRegLength = ($blkSizes[$i] + $SpliceSiteLen + $blkSizes[($i+1)]);

      ##put the line into the blat array    
      @b = ($blatOutput[0], $blatOutput[1], $blatOutput[2], $blatOutput[3], $blatOutput[4], $blatOutput[5], $blatOutput[6], $blatOutput[7], $blatOutput[8], $blatOutput[9], $blatOutput[10], $blatOutput[11], $blatOutput[12]);
      @b = (@b, $blatOutput[13], $blatOutput[14], $blatOutput[15], $blatOutput[16], $blatOutput[17], $blatOutput[18], $blatOutput[19], $blatOutput[20], $iStart, $SpliceSiteLen, $iLength, $includeRegStart, $includeRegLength, $i+1, $unqBLid);
      $blat[$#blat+1] = [ @b ];

      ##find the corresponding sequence
      ##search file for qname $blatOutput[9]
      ##-q option tells the program to use quality information when designing primers
      if ($options{q} == 1) {
        $out = `perl $extract_fasta_Path -q -c $blatOutput[9] $queryPATH $query.$blatOutput[9].fasta`;
        open(SEQ, "$query.$blatOutput[9].fasta");
        $line = <SEQ>; #remove header
        $seq = " ";
        while (<SEQ>) {
          chomp;
          $seq = "$seq$_";
        }
  
        open(QUAL, "$query.$blatOutput[9].fasta.qual");
        $line = <QUAL>; #remove header
        $qual = " ";
        while (<QUAL>) {
          chomp;
          $qual = "$qual $_";
        }
        
        close SEQ;
        $out = `rm "$query.$blatOutput[9].fasta"`;
        close QUAL;
        $out = `rm $query.$blatOutput[9].fasta.qual`; 

      } else {
        $out = `perl $extract_fasta_Path -c $blatOutput[9] $queryPATH $query.$blatOutput[9].fasta`;
        open(SEQ, "$query.$blatOutput[9].fasta");
        $line = <SEQ>; #remove header
        $seq = " ";
        while (<SEQ>) {
          chomp;
          $seq = "$seq$_";
        }

        close SEQ;
        $out = `rm "$query.$blatOutput[9].fasta"`;
      }##end of if else qual
      
      ##put in file for input to primer3
      print P3IN "PRIMER_SEQUENCE_ID=", $unqBLid, "_", $blatOutput[9], "_", $blatOutput[8], "_", $blatOutput[0], "/", $blatOutput[1], "_Intron", $i, "\n";
      print P3IN "SEQUENCE=$seq\n";
      if ($options{q} == 1) {print P3IN "PRIMER_SEQUENCE_QUALITY=$qual\n";}
      print P3IN "TARGET=$iStart,$SpliceSiteLen\n";
      print P3IN "INCLUDED_REGION=$includeRegStart,$includeRegLength\n";
      print P3IN "=\n";

      $unqBLid++;
    }##end of for each intron
  }##end of while blat output

  close BOUT;
  close P3IN;

  ##return the reference to @blat
  return \@blat;
}##Sub End of blat_parser

## Compiles output from RepeatMasker, BLAt and primer3 to create out.IMP
## Returns an array containing every record printed to out.IMP
## Arguemens include $query, $queryPATH, array containing blat records and one containing primer3 records
sub get_output_BLAT($$$$) {
  $query = $_[0];
  $queryPATH = $_[1];
  @blat = @{$_[2]};
  @primer3 = @{$_[3]};

  print "Generating program output...\n";
  
  ##parse primer3 output
  open(PO, ">$query.primerOrder.IMP");
  open(FO, ">$query.out.IMP");
  
	##array containing all output record arrays
	@output = ();

  ##get repeat masker information
  open(RM, "$queryPATH.out");
  <RM>; <RM>; <RM>; ##remove header
  %Smasked = ();
  while (<RM>) {
  	$_ =~ /^\s*\S*\s*\S*\s*\S*\s*\S*\s*(\S*).*$/;
  	if (exists $Smasked{$1}) {
  		$val = $Smasked{$1};
  		$Smasked{$1} = $val ++;
  	} else {
  		$Smasked{$1} = 1;
  	}
  }
  
  ##print headers
  print PO "Name\tLength\t%GC\tTm\tSequence\t\n";
  print FO &print_header("blat");
  
  for ($j=0; $j<= $#primer3_out; $j++) {
  	@recP3 = @{ $primer3_out[$j] };
  	@recP3_left = @{$recP3[1]};
  	$lp_seq = $recP3_left[4];
  	
    ##if there is at least 1 primer
  	if ($lp_seq =~ /[GACTgactNn]+/) {
  	  
  	  ##get blat portion
	    $seq_id = $recP3[0];
	    @blat_id = split('_', $seq_id); ##unique blat id is the 1st # of the Seq id

	    $index = $blat_id[0] - 100;  ##unique blat id = 100 + blat index
	    @recB = @{ $blat[$index]  };  ##blat record
	    @blat_id[$#blat_id] =~ /Intron(\d+)/; ##last part of seq id is Intron# where # is the current intron
	    $curr_intron = $1;

      if (exists $Smasked{$recB[9]}) {$rmNo = $Smasked{$recB[9]};} else {$rmNo = 0;}
      $iend = $recB[21] + $recB[23];
      $noI = $recB[17] -1;
      @blatOutput = ($recB[9], $rmNo, $recB[8], $recB[10], $recB[11], $recB[12], $recB[13], $recB[14], $recB[15], $recB[16], "$recB[0]/$recB[1]", $noI, $recB[($#recB-1)], $recB[21], $iend, $recB[23], $recB[22], "$recB[24],$recB[25]");

      ##For each primerset
      for($i=0; $i<5; $i++) {
        if ($lp_seq =~ /[GACTgactNn]+/) {
          $iS = ($i * 3) + 1;
          @recP3_left = @{$recP3[$iS]};
        	@recP3_right = @{$recP3[$iS+1]};
        	@recP3_pair = @{$recP3[$iS+2]};
      	
        	$lp_seq = $recP3_left[4];
    	
          $Gnm_prod_size = 0;
    	    $Gnm_prod_size = $recP3_pair[2] + $recB[23];
	    
          $ref = &get_primerset($i+1, $Gnm_prod_size, \@recP3_left, \@recP3_right, \@recP3_pair);
          @primerRec = @{$ref};
      
    	    print FO join(" \t ", @blatOutput), " \t "; ##prints the blat portion
    	    print FO join(" \t ", @primerRec), "\n";  ##prints the primer portion
    
          ##prints to the simple primer output
    	    print PO $recB[9], "_", $primerRec[0],"_Fwd\t ", $primerRec[5], " \t ", $primerRec[4], " \t ", $primerRec[2], " \t ", $primerRec[3], " \t ", $primerRec[0], "\n";
    	    print PO $recB[9], "_", $primerRec[0],"_Rvs\t ", $primerRec[10], " \t ", $primerRec[9], " \t ", $primerRec[7], " \t ", $primerRec[8], " \t ", $primerRec[0], "\n";
	    
    	    #adds into array
    	    @record = (@blatOutput, @primerRec, $j, $index); #last two are primer3 id and blat id
    	    $output[$#output+1] = [ @record ];
    	  } ##end of if sequence for current primerset
      } ##end of for each primerset

	  } ##end of if there is @ least 1 primerset

  } #end of for each primer3 output record
  
  
  close PO;
  close FO;
  close P3OUT;

	if ($options{c} =~ /G/) {@opt = (1);} else {@opt = (0);}
	if ($options{c} =~ /I/) {@opt = (@opt, 1);} else {@opt = (@opt, 0);}
	if ($options{c} =~ /P/) {@opt = (@opt, 1);} else {@opt = (@opt, 0);}
	&blat_autocurate($query, $queryPATH, \@output, \@opt);
	
}

##-----Primer3 related subroutines-----------------------------------------------------------------------------------------------

## Parses primer3 input for addition to final output
## Returns an array containing each primer3 record
## Arguements include query and database
sub parse_Primer3_out ($$) {
  $query = $_[0];
  $genome = $_[1];

  open(P3OUT, "$query.primer3.out.boulderIO");

  @p3Records = ();
  $index = 0;

  use Boulder::Stream;
  $stream = Boulder::Stream->new(*P3OUT,);
  while ($record = $stream->get()) {
    ##Primer set 1
    ##get the information for the primer portion of output
  	$seq_id = $record->PRIMER_SEQUENCE_ID;
  	$sd = "_$index";
  	$s = "$seq_id$sd";
  	@curRec = ($s);

  	$left_tm = $record->PRIMER_LEFT_TM;
  	$left_GC = $record->PRIMER_LEFT_GC_PERCENT;
  	$left_seq = $record->PRIMER_LEFT_SEQUENCE;
  	$left_len = length($left_seq);
  	@left_target = split(',', $record->PRIMER_LEFT);
  	$left_start = $left_target[0];
  	$left_penalty = $record->PRIMER_LEFT_PENALTY;
  	$left_selfAny = $record->PRIMER_LEFT_SELF_ANY;
  	$left_selfEnd = $record->PRIMER_LEFT_SELF_END;
  	$left_endStable = $record->PRIMER_LEFT_END_STABILITY;

  	$right_tm = $record->PRIMER_RIGHT_TM;
  	$right_GC = $record->PRIMER_RIGHT_GC_PERCENT;
  	$right_seq = $record->PRIMER_RIGHT_SEQUENCE;

  	$right_len = length($right_seq);
  	@right_target = split(',', $record->PRIMER_RIGHT);
  	$right_start = $right_target[0];
  	$right_penalty = $record->PRIMER_RIGHT_PENALTY;
  	$right_selfAny = $record->PRIMER_RIGHT_SELF_ANY;
  	$right_selfEnd = $record->PRIMER_RIGHT_SELF_END;
  	$right_endStable = $record->PRIMER_RIGHT_END_STABILITY;

  	$EST_prod_size = $record->PRIMER_PRODUCT_SIZE;
  	$pair_compl_any = $record->PRIMER_PAIR_COMPL_ANY;
  	$pair_compl_end = $record->PRIMER_PAIR_COMPL_END;
  	$pair_penalty = $record->PRIMER_PAIR_PENALTY;
  	$inclReg = $record->INCLUDED_REGION;

  	if ($left_seq =~ /[GACTgactNn]+/) {
	    @left = ("L", 1, $left_tm, $left_GC, $left_seq, $left_len, $record->PRIMER_LEFT, $left_start, $left_penalty, $left_selfAny, $left_selfEnd, $left_endStable);
	    @curRec[$#curRec+1] = [ @left ];
	    @right = ("R", 1, $right_tm, $right_GC, $right_seq, $right_len, $record->PRIMER_RIGHT, $right_start, $right_penalty, $right_selfAny, $right_selfEnd, $right_endStable);
	    @curRec[$#curRec+1] = [ @right ];
	    @pair = ("I", 1, $EST_prod_size, $pair_compl_any, $pair_compl_end, $pair_penalty, $inclReg);
	    @curRec[$#curRec+1] = [ @pair ];
	  } #end of if this record contains a primer description

    ##Primerset 2
    ##get the information for the primer portion of output
  	$left_tm = $record->PRIMER_LEFT_1_TM;
  	$left_GC = $record->PRIMER_LEFT_1_GC_PERCENT;
  	$left_seq = $record->PRIMER_LEFT_1_SEQUENCE;
  	$left_len = length($left_seq);
  	@left_target = split(',', $record->PRIMER_LEFT_1);
  	$left_start = $left_target[0];
  	$left_penalty = $record->PRIMER_LEFT_1_PENALTY;
  	$left_selfAny = $record->PRIMER_LEFT_1_SELF_ANY;
  	$left_selfEnd = $record->PRIMER_LEFT_1_SELF_END;
  	$left_endStable = $record->PRIMER_LEFT_1_END_STABILITY;
    
  	$right_tm = $record->PRIMER_RIGHT_1_TM;
  	$right_GC = $record->PRIMER_RIGHT_1_GC_PERCENT;
  	$right_seq = $record->PRIMER_RIGHT_1_SEQUENCE;        
  	$right_len = length($right_seq);
  	@right_target = split(',', $record->PRIMER_RIGHT_1);
  	$right_start = $right_target[0];
  	$right_penalty = $record->PRIMER_RIGHT_1_PENALTY;
  	$right_selfAny = $record->PRIMER_RIGHT_1_SELF_ANY;
  	$right_selfEnd = $record->PRIMER_RIGHT_1_SELF_END;
  	$right_endStable = $record->PRIMER_RIGHT_1_END_STABILITY;

  	$EST_prod_size = $record->PRIMER_PRODUCT_SIZE_1;
  	$pair_compl_any = $record->PRIMER_PAIR_1_COMPL_ANY;
  	$pair_compl_end = $record->PRIMER_PAIR_1_COMPL_END;
  	$pair_penalty = $record->PRIMER_PAIR_PENALTY_1;
  	$inclReg = $record->INCLUDED_REGION;

  	if ($left_seq =~ /[GACTgactNn]+/) {
	    @left = ("L", 2, $left_tm, $left_GC, $left_seq, $left_len, $record->PRIMER_LEFT_1, $left_start, $left_penalty, $left_selfAny, $left_selfEnd, $left_endStable);
	    @curRec[$#curRec+1] = [ @left ];
	    @right = ("R", 2, $right_tm, $right_GC, $right_seq, $right_len, $record->PRIMER_RIGHT_1, $right_start, $right_penalty, $right_selfAny, $right_selfEnd, $right_endStable);
	    @curRec[$#curRec+1] = [ @right ];
	    @pair = ("I", 2, $EST_prod_size, $pair_compl_any, $pair_compl_end, $pair_penalty, $inclReg);
	    @curRec[$#curRec+1] = [ @pair ];
	    #print "Left: @left\n Right: @right\n Pair: @pair\n";
	  } #end of if this record contains a primer description

    ##Primerset 3
    ##get the information for the primer portion of output
  	$left_tm = $record->PRIMER_LEFT_2_TM;
  	$left_GC = $record->PRIMER_LEFT_2_GC_PERCENT;
  	$left_seq = $record->PRIMER_LEFT_2_SEQUENCE;    
  	$left_len = length($left_seq);
  	@left_target = split(',', $record->PRIMER_LEFT_2);
  	$left_start = $left_target[0];
  	$left_penalty = $record->PRIMER_LEFT_2_PENALTY;
  	$left_selfAny = $record->PRIMER_LEFT_2_SELF_ANY;
  	$left_selfEnd = $record->PRIMER_LEFT_2_SELF_END;
  	$left_endStable = $record->PRIMER_LEFT_2_END_STABILITY;
    
  	$right_tm = $record->PRIMER_RIGHT_2_TM;
  	$right_GC = $record->PRIMER_RIGHT_2_GC_PERCENT;
  	$right_seq = $record->PRIMER_RIGHT_2_SEQUENCE;
  	$right_len = length($right_seq);
  	@right_target = split(',', $record->PRIMER_RIGHT_2);
  	$right_start = $right_target[0];
  	$right_penalty = $record->PRIMER_RIGHT_2_PENALTY;
  	$right_selfAny = $record->PRIMER_RIGHT_2_SELF_ANY;
  	$right_selfEnd = $record->PRIMER_RIGHT_2_SELF_END;
  	$right_endStable = $record->PRIMER_RIGHT_2_END_STABILITY;

  	$EST_prod_size = $record->PRIMER_PRODUCT_SIZE_2;
  	$pair_compl_any = $record->PRIMER_PAIR_2_COMPL_ANY;
  	$pair_compl_end = $record->PRIMER_PAIR_2_COMPL_END;
  	$pair_penalty = $record->PRIMER_PAIR_PENALTY_2;
  	$inclReg = $record->INCLUDED_REGION;

  	if ($left_seq =~ /[GACTgactNn]+/) {
	    @left = ("L", 3, $left_tm, $left_GC, $left_seq, $left_len, $record->PRIMER_LEFT_2, $left_start, $left_penalty, $left_selfAny, $left_selfEnd, $left_endStable);
	    @curRec[$#curRec+1] = [ @left ];
	    @right = ("R", 3, $right_tm, $right_GC, $right_seq, $right_len, $record->PRIMER_RIGHT_2, $right_start, $right_penalty, $right_selfAny, $right_selfEnd, $right_endStable);
	    @curRec[$#curRec+1] = [ @right ];
	    @pair = ("I", 3, $EST_prod_size, $pair_compl_any, $pair_compl_end, $pair_penalty, $inclReg);
	    @curRec[$#curRec+1] = [ @pair ];
	  } #end of if this record contains a primer description

    ##Primerset 4
    ##get the information for the primer portion of output
  	$left_tm = $record->PRIMER_LEFT_3_TM;
  	$left_GC = $record->PRIMER_LEFT_3_GC_PERCENT;
  	$left_seq = $record->PRIMER_LEFT_3_SEQUENCE;
  	$left_len = length($left_seq);
  	@left_target = split(',', $record->PRIMER_LEFT_3);
  	$left_start = $left_target[0];
  	$left_penalty = $record->PRIMER_LEFT_3_PENALTY;
  	$left_selfAny = $record->PRIMER_LEFT_3_SELF_ANY;
  	$left_selfEnd = $record->PRIMER_LEFT_3_SELF_END;
  	$left_endStable = $record->PRIMER_LEFT_3_END_STABILITY;
    
  	$right_tm = $record->PRIMER_RIGHT_3_TM;
  	$right_GC = $record->PRIMER_RIGHT_3_GC_PERCENT;
  	$right_seq = $record->PRIMER_RIGHT_3_SEQUENCE;
  	$right_len = length($right_seq);
  	@right_target = split(',', $record->PRIMER_RIGHT_3);
  	$right_start = $right_target[0];
  	$right_penalty = $record->PRIMER_RIGHT_3_PENALTY;
  	$right_selfAny = $record->PRIMER_RIGHT_3_SELF_ANY;
  	$right_selfEnd = $record->PRIMER_RIGHT_3_SELF_END;
  	$right_endStable = $record->PRIMER_RIGHT_3_END_STABILITY;

  	$EST_prod_size = $record->PRIMER_PRODUCT_SIZE_3;
  	$pair_compl_any = $record->PRIMER_PAIR_3_COMPL_ANY;
  	$pair_compl_end = $record->PRIMER_PAIR_3_COMPL_END;
  	$pair_penalty = $record->PRIMER_PAIR_PENALTY_3;
  	$inclReg = $record->INCLUDED_REGION;

  	if ($left_seq =~ /[GACTgactNn]+/) {
	    @left = ("L", 4, $left_tm, $left_GC, $left_seq, $left_len, $record->PRIMER_LEFT_3, $left_start, $left_penalty, $left_selfAny, $left_selfEnd, $left_endStable);
	    @curRec[$#curRec+1] = [ @left ];
	    @right = ("R", 4, $right_tm, $right_GC, $right_seq, $right_len, $record->PRIMER_RIGHT_3, $right_start, $right_penalty, $right_selfAny, $right_selfEnd, $right_endStable);
	    @curRec[$#curRec+1] = [ @right ];
	    @pair = ("I", 4, $EST_prod_size, $pair_compl_any, $pair_compl_end, $pair_penalty, $inclReg);
	    @curRec[$#curRec+1] = [ @pair ];
	  } #end of if this record contains a primer description

    ##Primerset 5
    ##get the information for the primer portion of output
  	$left_tm = $record->PRIMER_LEFT_4_TM;
  	$left_GC = $record->PRIMER_LEFT_4_GC_PERCENT;
  	$left_seq = $record->PRIMER_LEFT_4_SEQUENCE;
  	$left_len = length($left_seq);
  	@left_target = split(',', $record->PRIMER_LEFT_4);
  	$left_start = $left_target[0];
  	$left_penalty = $record->PRIMER_LEFT_4_PENALTY;
  	$left_selfAny = $record->PRIMER_LEFT_4_SELF_ANY;
  	$left_selfEnd = $record->PRIMER_LEFT_4_SELF_END;
  	$left_endStable = $record->PRIMER_LEFT_4_END_STABILITY;

  	$right_tm = $record->PRIMER_RIGHT_4_TM;
  	$right_GC = $record->PRIMER_RIGHT_4_GC_PERCENT;
  	$right_seq = $record->PRIMER_RIGHT_4_SEQUENCE;
  	$right_len = length($right_seq);
  	@right_target = split(',', $record->PRIMER_RIGHT_4);
  	$right_start = $right_target[0];
  	$right_penalty = $record->PRIMER_RIGHT_4_PENALTY;
  	$right_selfAny = $record->PRIMER_RIGHT_4_SELF_ANY;
  	$right_selfEnd = $record->PRIMER_RIGHT_4_SELF_END;
  	$right_endStable = $record->PRIMER_RIGHT_4_END_STABILITY;

  	$EST_prod_size = $record->PRIMER_PRODUCT_SIZE_4;
  	$pair_compl_any = $record->PRIMER_PAIR_4_COMPL_ANY;
  	$pair_compl_end = $record->PRIMER_PAIR_4_COMPL_END;
  	$pair_penalty = $record->PRIMER_PAIR_PENALTY_4;
  	$inclReg = $record->INCLUDED_REGION;

  	if ($left_seq =~ /[GACTgactNn]+/) {
	    @left = ("L", 5, $left_tm, $left_GC, $left_seq, $left_len, $record->PRIMER_LEFT_4, $left_start, $left_penalty, $left_selfAny, $left_selfEnd, $left_endStable);
	    @curRec[$#curRec+1] = [ @left ];
	    @right = ("R", 5, $right_tm, $right_GC, $right_seq, $right_len, $record->PRIMER_RIGHT_4, $right_start, $right_penalty, $right_selfAny, $right_selfEnd, $right_endStable);
	    @curRec[$#curRec+1] = [ @right ];
	    @pair = ("I", 5, $EST_prod_size, $pair_compl_any, $pair_compl_end, $pair_penalty, $inclReg);
	    @curRec[$#curRec+1] = [ @pair ];
	  } #end of if this record contains a primer description

  	@r = @{$curRec[1]};
  	if ($r[4] =~ /[GACTgactNn]+/) {
  	    @p3Records[$#p3Records+1] = [@curRec];
  	}
  	
	  $index ++;
  } #end of for each primer3 record

  close P3OUT;
  close STS;

  print "Primer sets were found for $#p3Records of the GeneSeqer Alignments...\n";
  return \@p3Records;

} ##end of parse_Primer3_out

##-----Seq Curation related subroutines-----------------------------------------------------------------------------------------------

## Removes any dashes and compresses sequence data into a single line for input into RepeatMasker ie:clenses the sequence files
## Returns nothing
## Arguements include $queryPATH, $database & reference to options hash
sub rmDash($$$) {

  my $query = $_[0];
  my $database = $_[1];
  my %options = %{$_[2]};

  if ($options{d} eq "q") { @files = ($query);}
  elsif ($options{d} eq "t") { @files = ($database);}
  elsif ($options{d} eq "b") { @files = ($query, $database);}
  elsif ($options{d} eq "n") { return;}
  else { print "ERROR: unrecognized option $option{d} for -d\n";}

  foreach $file (@files) {
  	open(Q, "$file");
  	open(NEWQ, ">$file.rmDash");

  	$count = 0;
  	$currRec = "";
  	print "\nVerification of format of $file...\n";
	
  	while (<Q>) {
	    if ($_ =~ /^>(\S*)\s*$/) {
	    	if ($currRec =~ /^>(\S*)\s*\n(?:\s*[GACTgactNn]+\s*)+/) {
	    		print NEWQ $currRec, "\n";
	    	}
  			$currRec = "$_";
  			$PriorityIDs = "$PriorityIDs\t$_";
  			$PROGRESS = "HEADER";
	    } elsif ($_ =~ /^(?:\s*[GACTgactNn]+\s*)+$/) {
	  		chomp;
  			if($PROGRESS eq "HEADER") {$currRec = $currRec.$_;}
  			if($PROGRESS eq "SEQ") {$currRec = $currRec.$_;}
  			$PROGRESS = "SEQ";
	    } elsif ($_ =~ /^(?:\s*[GACTgactNn-]+\s*)+$/) {
	  		chomp;
  			$_ =~ tr/-/ /;
  			if($PROGRESS eq "HEADER") {$currRec = $currRec.$_;}
  			if($PROGRESS eq "SEQ") {$currRec = $currRec.$_;}
  			$PROGRESS = "SEQ";
	    } else {
  			print "ERROR: Sequence in an unrecognized format: $_\n";
  			$PROGRESS = "ERROR";
  			$count ++;
	    } ##end of if header or sequence
  	} #end of while input from Query

  	#print last record
  	if ($currRec =~ /^>(\S*)\s*\n(?:\s*[GACTgactNn]+\s*)+/) {
  		print NEWQ $currRec, "\n";
  	}
	    	
  	print "$count Errors found in $file.\n";	
  	close Q;
  	close NEWQ;
  }##end of foreach file
}##Sub End rmDash

## Curates the out.IMP dataset by a number of user specified criteria (ie: only output 1 intron per Gene Model)
## Returns nothing
## Arguements include query, query path, out.IMP array & options
sub geneseqer_autocurate ($$$$) {
  $query = $_[0];
  $queryPATH = $_[1];
  @output = @{ $_[2] };
  @opt = @{ $_[3] };	
	
	print "Autocurating $query.out.IMP for:\n"; 
	##implement each level of autocuration independantly allowing any combination
	##after each level have a hash containing arrays with all output feilds	
	
	#Create hash for autocuration
	#Order of keys ESTname -> GenemodelNo. -> IntronNo. -> P3recordNo. -> PrimersetNo.
	%curated = ();
	foreach $ref (@output) { #for each line of output
		@rec = @{$ref};
		$EST = $rec[0];
		$GM = $rec[$#rec];
		$I = $rec[11];
		$P3 = $rec[$#rec-1];
		$PNo = $rec[19];
		pop(@rec); pop(@rec); ##remove last two elements (don't want to print p3 and blat id)
		##adds each line from the out.IMP to the hash
		$curated{$EST}{$GM}{$I}{$P3}{$PNo} = [@rec];
	}

	##prints the hash to the file
	open(T,">$query.sorted.IMP");
  print T &print_header("geneseqer");
  
	while( my ($ik, $iv) = each %curated ) {
		print T "\nEST: $ik\n";
		while( my ($jk, $jv) = each %{$iv} ) {
			print T "\tGene Model: $jk\n";
			while( my ($kk, $kv) = each %{$jv} ) {
				print T "\t\tIntron $kk\n";
				while( my ($lk, $lv) = each %{$kv} ) {
					print T "\t\t\tP3 Record: $lk\n";
					while( my ($mk, $mv) = each %{$lv} ) {
						print T "Primerset:$mk = ", join("\t", @{$mv}), "\n";	
					}
				}	
			}
		}
	}
	close T;
	
	open(O,">$query.curated.details");
	
	##curate: pick only 1 gene model per query EST (the one with the highest similarity)
	if ($opt[0]  == 1) {
		print "G)\tGene model with highest similarity per query EST\n";
		print O "Below are the details for determining which is the best gene model (highest similarity between query and template) per EST. If two gene models have the exact same similarity then both are kept.\n";
		%curated1 = ();
		##for each EST
		while( my ($ik, $iv) = each %curated ) {
			print O "\nEST: $ik\n";
			%GeneModels = %{$iv};
			##for each Genemodel of a given EST
			$bestSim = 0;
			$curSim = 0;
			@bestRec = ();
			@introns = ();
			while( my ($jk, $jv) = each %GeneModels) {
				print O "\tGene Model: $jk\n";
				##pick the first element
				while (my ($kk, $kv) = each(%{$jv}) ) {
					while (my ($lk, $lv) = each(%{$kv})) {
						while (my ($mk, $mv) = each(%{$lv})) {
							$curSim = @{$mv}[9];
							print O "Similarity: $curSim\tIntron:$kk, Primerset:$mk\n";
							if ($curSim > $bestSim) {
								@bestRec = ("$jk,$kk,$lk,$mk");
								$bestSim = $curSim;
							} elsif ($curSim == $bestSim) {
								@bestRec = (@bestRec, "$jk,$kk,$lk,$mk");
							} else {
								#keep last one as bestRec
							}
							$curSim = 0;
						}
					}
				}
			}
			
			for ($a=0; $a <= $#bestRec; $a++) {
				@indices = split(",", $bestRec[$a]);
				print O "Best] GeneModel:$indices[0], Intron:$indices[1], Primerset:$indices[3]\n";
				$curated1{$ik}{$indices[0]}{$indices[1]}{$indices[2]}{$indices[3]} = [ @{ $curated{$ik}{$indices[0]}{$indices[1]}{$indices[2]}{$indices[3]} } ];
			}
		}
		%curated = %curated1;
	}

  ##curate: pick only 1 intron per gene model (the one with the highest simimlarity score)
	if ($opt[1]  == 1) {
		print "I)\tintron with highest similarity score per query EST Gene model(EST/BAC alignment)\n";
		print O "\n\n", "-"x200;
		print O "\nBelow are the details outlining the process for chosing the best intron (Highest average donor/acceptor similarity) per Gene model. If two introns have the same average similarity then both are kept.\n";
		%curated2 = ();
		#for each EST
		while( my ($ik, $iv) = each %curated ) {
			print O "\nEST: $ik\n";
			#for each Genemodel of a given EST
			while( my ($jk, $jv) = each %{$iv}) {
				print O "\tGene Model: $jk\n";
				$bestSim = 0;
				$curSim = 0;
				@bestRec = ();
				#for each intron
				while( my ($kk, $kv) = each %{$jv}) {
					#pick the first element
					while (my ($lk, $lv) = each(%{$kv})) {
						while (my ($mk, $mv) = each(%{$lv})) {
							$curSim = (@{$mv}[19]+@{$mv}[20])/2;
							print O "Similarity: $curSim\tIntron:$kk, Primerset:$mk\n";
							if ($curSim >= $bestSim) {
								if ($curSim > $bestSim) {@bestRec = ();}
								@bestRec = (@bestRec, "$kk,$lk,$mk");
								$bestSim = $curSim;
							} else {
								#keep last one as bestRec
							}
							$curSim = 0;
						}
					}
				}
				
				for ($a=0; $a <= $#bestRec; $a++) {
					@indices = split(",", $bestRec[$a]);
					print O "Best] Intron:$indices[0], Primerset:$indices[2]\n";
					$curated2{$ik}{$jk}{$indices[0]}{$indices[1]}{$indices[2]} = [ @{ $curated{$ik}{$jk}{$indices[0]}{$indices[1]}{$indices[2]}} ];
				}
			}
		}
		%curated = %curated2;
	}
	
	##curate: pick only 1 primer pair per intron (the one with the lowest primer penalty)
	if ($opt[2]  == 1) {
		print "P)\tPrimer pair with lowest primer pair penalty per intron\n";
		print O "\n\n", "-"x200;
		print O "\nBelow are the details outlining the process of choosing the best primer pair (the lowest primer pair penalty) per intron. If two primer pairs have the same primer pair penalty both will be kept.\n";
		%curated3 = ();
		#for each EST
		while( my ($ik, $iv) = each %curated ) {
			print O "\nEST: $ik\n";
			#for each Genemodel of a given EST
			while( my ($jk, $jv) = each %{$iv}) {
				print O "\tGene Model: $jk\n";
				#for each intron
				while( my ($kk, $kv) = each %{$jv}) {
					print O "\t\tIntron $kk\n";
					$bestPen = 100;
					$curPen = 100;
					@bestRec = ();
					#for each primer3 record
					while( my ($lk, $lv) = each %{$kv} ) {
						print O "\t\t\tP3 Record: $lk\n";
						#for each primer pair
						while( my ($mk, $mv) = each %{$lv} ) {
							$curPen = @{$mv}[34];
							print O "Pair Penalty: $curPen\tPrimerset:$mk\n";
							if ($curPen <= $bestPen) {
								if ($curPen < $bestPen) {@bestRec = ();}
								@bestRec = (@bestRec, "$kk,$lk,$mk");
								$bestPen = $curPen;
							} else {
							#keep last one as bestRec
							}
						$curPen = 100;
						}
					}
					
					for ($a=0; $a <= $#bestRec; $a++) {
						@indices = split(",", $bestRec[$a]);
						print O "Best] Primerset:$indices[2]\n";
						$curated3{$ik}{$jk}{$indices[0]}{$indices[1]}{$indices[2]} = [ @{ $curated{$ik}{$jk}{$indices[0]}{$indices[1]}{$indices[2]}} ];
					}
				}
			}
		}
		%curated = %curated3;
	}
	
	close O;
	
	##print final curated results
	open (CURATED, ">$query.curated.IMP");
  print CURATED &print_header("geneseqer");
  
	while( my ($ik, $iv) = each %curated ) {
		print CURATED "\nEST: $ik\n";
		while( my ($jk, $jv) = each %{$iv} ) {
			print CURATED "\tGene Model: $jk\n";
			while( my ($kk, $kv) = each %{$jv} ) {
				print CURATED "\t\tIntron $kk\n";
				while( my ($lk, $lv) = each %{$kv} ) {
					while( my ($mk, $mv) = each %{$lv} ) {
						print CURATED "Primerset: $mk = ", join("\t", @{$mv}), "\n";	
					}
				}	
			}
		}
	} ##while curated records to print
	close CURATED;

}##Sub End autocurate

## Curates the out.IMP dataset by a number of user specified criteria (ie: only output 1 intron per Gene Model)
## Returns nothing
## Arguements include query, query path, out.IMP array & options
sub blat_autocurate ($$$$) {
  $query = $_[0];
  $queryPATH = $_[1];
  @output = @{ $_[2] };
  @opt = @{ $_[3] };	
	
	print "Autocurating $query.out.IMP for:\n"; 
	##implement each level of autocuration independantly allowing any combination
	##after each level have a hash containing arrays with all output feilds	
	
	##Create hash for autocuration
	##Order of keys ESTname -> GenemodelNo. -> IntronNo. -> P3recordNo. -> PrimersetNo.
	%curated = ();
	foreach $ref (@output) { ##for each line of output
		@rec = @{$ref};
		$EST = $rec[0];
		$GM = $rec[$#rec];
		$I = $rec[12];
		$P3 = $rec[$#rec-1];
		$PNo = $rec[18];
		pop(@rec); pop(@rec); ##remove last two elements (don't want to print p3 and blat id)
		##adds each line from the out.IMP to the hash
		$curated{$EST}{$GM}{$I}{$P3}{$PNo} = [@rec];
	}

	##prints the hash to the file
	open(T,">$query.sorted.IMP");
  print T &print_header("blat");
  
	while( my ($ik, $iv) = each %curated ) {
		print T "\nEST: $ik\n";
		while( my ($jk, $jv) = each %{$iv} ) {
			print T "\tGene Model: $jk\n";
			while( my ($kk, $kv) = each %{$jv} ) {
				print T "\t\tIntron $kk\n";
				while( my ($lk, $lv) = each %{$kv} ) {
					print T "\t\t\tP3 Record: $lk\n";
					while( my ($mk, $mv) = each %{$lv} ) {
						print T "Primerset:$mk = ", join("\t", @{$mv}), "\n";	
					}
				}	
			}
		}
	}
	close T;
	
	open(O,">$query.curated.details");
	
	##curate: pick only 1 gene model per query EST (the one with the highest similarity)
	if ($opt[0]  == 1) {
		print "G)\tGene model with highest match/mismatch ratio per query EST\n";
		print O "Below are the details for determining which is the best gene model (highest similarity between query and template) per EST. If two gene models have the exact same similarity then both are kept.\n";
		%curated1 = ();
		##for each EST
		while( my ($ik, $iv) = each %curated ) {
			print O "\nEST: $ik\n";
			%GeneModels = %{$iv};
			##for each Genemodel of a given EST
			$bestSim = 0;
			$curSim = 0;
			@bestRec = ();
			@introns = ();
			while( my ($jk, $jv) = each %GeneModels) {
				print O "\tGene Model: $jk\n";
				##pick the first element
				while (my ($kk, $kv) = each(%{$jv}) ) {
					while (my ($lk, $lv) = each(%{$kv})) {
						while (my ($mk, $mv) = each(%{$lv})) {
							$curSim = eval(@{$mv}[10]);
							print O "Similarity: $curSim\tIntron:$kk, Primerset:$mk\n";
							if ($curSim > $bestSim) {
								@bestRec = ("$jk,$kk,$lk,$mk");
								$bestSim = $curSim;
							} elsif ($curSim == $bestSim) {
								@bestRec = (@bestRec, "$jk,$kk,$lk,$mk");
							} else {
								#keep last one as bestRec
							}
							$curSim = 0;
						}
					}
				}
			}
			
			for ($a=0; $a <= $#bestRec; $a++) {
				@indices = split(",", $bestRec[$a]);
				print O "Best] GeneModel:$indices[0], Intron:$indices[1], Primerset:$indices[3]\n";
				$curated1{$ik}{$indices[0]}{$indices[1]}{$indices[2]}{$indices[3]} = [ @{ $curated{$ik}{$indices[0]}{$indices[1]}{$indices[2]}{$indices[3]} } ];
			}
		}
		%curated = %curated1;
	}

  ##curate: pick only 1 intron per gene model (the one with the highest simimlarity score)
	if ($opt[1]  == 1) {
    print "I)\tintron with highest similarity score per query EST Gene model(EST/BAC alignment)\n";
    print "\tSorry but results cannot be filtered based on intron predictions when using BLAT\n";
    
    print O "Results cannot be filtered based on intron predictions when using BLAT\n";
    # print O "\n\n", "-"x200;
    # print O "\nBelow are the details outlining the process for chosing the best intron (Highest average donor/acceptor similarity) per Gene model. If two introns have the same average similarity then both are kept.\n";
    # %curated2 = ();
    # #for each EST
    # while( my ($ik, $iv) = each %curated ) {
    #   print O "\nEST: $ik\n";
    #   #for each Genemodel of a given EST
    #   while( my ($jk, $jv) = each %{$iv}) {
    #     print O "\tGene Model: $jk\n";
    #     $bestSim = 0;
    #     $curSim = 0;
    #     @bestRec = ();
    #     #for each intron
    #     while( my ($kk, $kv) = each %{$jv}) {
    #       #pick the first element
    #       while (my ($lk, $lv) = each(%{$kv})) {
    #         while (my ($mk, $mv) = each(%{$lv})) {
    #           $curSim = (@{$mv}[19]+@{$mv}[20])/2;
    #           print O "Similarity: $curSim\tIntron:$kk, Primerset:$mk\n";
    #           if ($curSim >= $bestSim) {
    #             if ($curSim > $bestSim) {@bestRec = ();}
    #             @bestRec = (@bestRec, "$kk,$lk,$mk");
    #             $bestSim = $curSim;
    #           } else {
    #             #keep last one as bestRec
    #           }
    #           $curSim = 0;
    #         }
    #       }
    #     }
    #     
    #     for ($a=0; $a <= $#bestRec; $a++) {
    #       @indices = split(",", $bestRec[$a]);
    #       print O "Best] Intron:$indices[0], Primerset:$indices[2]\n";
    #       $curated2{$ik}{$jk}{$indices[0]}{$indices[1]}{$indices[2]} = [ @{ $curated{$ik}{$jk}{$indices[0]}{$indices[1]}{$indices[2]}} ];
    #     }
    #   }
    # }
    # %curated = %curated2;
	}
	
	##curate: pick only 1 primer pair per intron (the one with the lowest primer penalty)
	if ($opt[2]  == 1) {
		print "P)\tPrimer pair with lowest primer pair penalty per intron\n";
		print O "\n\n", "-"x200;
		print O "\nBelow are the details outlining the process of choosing the best primer pair (the lowest primer pair penalty) per intron. If two primer pairs have the same primer pair penalty both will be kept.\n";
		%curated3 = ();
		#for each EST
		while( my ($ik, $iv) = each %curated ) {
			print O "\nEST: $ik\n";
			#for each Genemodel of a given EST
			while( my ($jk, $jv) = each %{$iv}) {
				print O "\tGene Model: $jk\n";
				#for each intron
				while( my ($kk, $kv) = each %{$jv}) {
					print O "\t\tIntron $kk\n";
					$bestPen = 100;
					$curPen = 100;
					@bestRec = ();
					#for each primer3 record
					while( my ($lk, $lv) = each %{$kv} ) {
						print O "\t\t\tP3 Record: $lk\n";
						#for each primer pair
						while( my ($mk, $mv) = each %{$lv} ) {
							$curPen = @{$mv}[31];
							print O "Pair Penalty: $curPen\tPrimerset:$mk\n";
							if ($curPen <= $bestPen) {
								if ($curPen < $bestPen) {@bestRec = ();}
								@bestRec = (@bestRec, "$kk,$lk,$mk");
								$bestPen = $curPen;
							} else {
							#keep last one as bestRec
							}
						$curPen = 100;
						}
					}
					
					for ($a=0; $a <= $#bestRec; $a++) {
						@indices = split(",", $bestRec[$a]);
						print O "Best] Primerset:$indices[2]\n";
						$curated3{$ik}{$jk}{$indices[0]}{$indices[1]}{$indices[2]} = [ @{ $curated{$ik}{$jk}{$indices[0]}{$indices[1]}{$indices[2]}} ];
					}
				}
			}
		}
		%curated = %curated3;
	}
	
	close O;
	
	##print final curated results
	open (CURATED, ">$query.curated.IMP");
  print CURATED &print_header("blat");
  
	while( my ($ik, $iv) = each %curated ) {
		print CURATED "\nEST: $ik\n";
		while( my ($jk, $jv) = each %{$iv} ) {
			print CURATED "\tGene Model: $jk\n";
			while( my ($kk, $kv) = each %{$jv} ) {
				print CURATED "\t\tIntron $kk\n";
				while( my ($lk, $lv) = each %{$kv} ) {
					while( my ($mk, $mv) = each %{$lv} ) {
						print CURATED "Primerset: $mk = ", join("\t", @{$mv}), "\n";	
					}
				}	
			}
		}
	} ##while curated records to print
	close CURATED;

}##Sub End blat_autocurate
