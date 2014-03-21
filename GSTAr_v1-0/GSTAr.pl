#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my $version_number = "1.0";
my $help = help_message($version_number);

# if there are no arguments, return the help message and quit
unless($ARGV[0]) {
    die "$help";
}

# declare global variables
our($opt_q,$opt_r,$opt_h,$opt_v,$opt_t,$opt_a);

# Initialize option defaults
$opt_r = 0.65;


# get options
getopts('r:hvqta');

# If user wants help or version, give it, and die
if($opt_h) {
    die "$help";
}
if($opt_v) {
    die "GSTAr.pl version $version_number\n";
}

# ensure transcripts and queries files are readable
my $transcripts_file = pop @ARGV;
unless(-r $transcripts_file) {
    die "FATAL: transcripts.fasta file was not readable\n\n$help";
}
my $queries_file = pop @ARGV;
unless(-r $queries_file) {
    die "FATAL: queries.fasta file was not readable\n\n$help";
}

# Start speaking to user, unless quiet mode is on
unless($opt_q) {
    print STDERR "\nGSTAr version $version_number\n";
    print STDERR `date`;
    print STDERR "\nDependency checks:\n";
}

# Check dependencies
my $RNAplex_check = check_RNAplex();
unless($RNAplex_check) {
    die "FAIL\n\n$help";
}

unless($opt_q) {
    print STDERR "\n";
}

# Validate options

unless(($opt_r > 0) and ($opt_r <= 1)) {
    die "FATAL: option r must be greater than 0 and less than or equal to 1\n\n$help";
}

# Hash the transcriptome into memory
my %txome = hash_transcriptome($transcripts_file);

# Begin work, one query at a time
unless($opt_q) {
    print STDERR "\n";
}

# Print headers for output
print "\# GSTAr version $version_number\n";
print "\# ";
print `date`;
print "\# Queries: $queries_file\n";
print "\# Transcripts: $transcripts_file\n";
print "\# Minimum Free Energy Ratio cutoff (option -r): $opt_r\n";
print "\# Sorted by: ";
if($opt_a) {
    print "AllenScore\n";
} else {
    print "MFEratio\n";
}
print "\# Output Format: ";
if($opt_t) {
    print "Tabular\n";
    # print a header
    print "Query\tTranscript\tTStart\tTStop\tTSlice\tMFEperfect\tMFEsite\tMFEratio\tAllenScore\tPaired\tUnpaired\tStructure\tSequence\n";
} else {
    print "Pretty\n";
}

(open(QUERIES, "$queries_file")) || die "FATAL: Failed to open queries file $queries_file\n\n$help";
while (<QUERIES>) {
    my $header = $_;
    my $full_qseq = uc <QUERIES>;

    # remove newlines and the > from header
    $header =~ s/\r//g;
    $header =~ s/\n//g;
    $header =~ s/>//g;
    $full_qseq =~ s/\r//g;
    $full_qseq =~ s/\n//g;
    
    # Report to user
    unless($opt_q) {
	print STDERR "Working on query $header ...";
    }
    
    # Crude cleanup and validation of sequence
    $full_qseq =~ s/U/T/g; ## U's to T's if present
    unless($full_qseq =~ /^[ATGC]+$/) {
	## not suppressed by quiet mode, because it is a warning
	print STDERR "\n\tWARNING: query sequence $header has non AUTGC characters. It is not being analyzed\!\n";
	next;
    }
    if((length $full_qseq) > 26) {
	print STDERR "\n\tWARNING: query sequence $header is too large. The limit is 26nts. It is not being analyzed\!\n";
	next;
    }
    if((length $full_qseq) < 15) {
	print STDERR "\n\tWARNING: query sequence $header is too small. It must be larger than 14 nts. It is not being analyzed\!\n";
	next;
    }
    
    $full_qseq = uc $full_qseq;
    $full_qseq =~ s/T/U/g;
    
    # Write the query fasta file and ge the tx lengths in the process.
    my $q_file = write_q(\$full_qseq,\%txome);  ## The file name is GSTAr_temp.fasta
    
    # Call the master sub-routine for RNAplex analyses
    my @first_pass = RNAplex_analysis(\$full_qseq,\$header,\%txome, \$q_file);
    
    # Sort, and in the process eliminate redundant alignments based on slice site within the query set
    my @winners = sort_first_pass(@first_pass);

    # Output
    unless($opt_q) {
	my $final_count = scalar @winners;
	if($final_count) {
	    print STDERR " Found $final_count valid alignments\n";
	} else {
	    print STDERR " NO valid alignments found\n";
	}
    }
    if($opt_t) {
	tabular_output(@winners);
    } else {
	pretty_output(@winners);
    }
}
close QUERIES;

# remove temp file
system "rm -f GSTAr_temp.fasta";

###########################
# Here be sub-routines

sub help_message {
    my($version) = @_;
    my $message = "GSTAr.pl : Generic Small RNA-Transcriptome Aligner
Usage: GSTAr.pl \[options\] queries.fasta transcripts.fasta

Version: $version

Options:
-h Print help message and quit
-v Print version and quit
-q Quiet mode .. no log/progress information to STDERR
-t Tabular output format ... More suitable for automated parsing
-a Sort by Allen et al. score instead of the default MFEratio
-r [float >0..1] Minimum Free Energy Ratio cutoff. Default: 0.65


Dependencies (must be in PATH):
  RNAplex (from Vienna RNA Package)

Documentation: perldoc GSTAr.pl
";
    return $message;
}

sub check_RNAplex {
    unless($opt_q) {
	print STDERR "\tRNAplex: ";
    }
    (open(RD, "RNAplex --version 2>&1 |")) || return 0;
    my $one = <RD>;
    close RD;
    if($one =~ /^RNAplex \d/) {
	unless($opt_q) {
	    print STDERR "PASS $one";
	}
	return 1;
    } else {
	return 0;
    }
}

sub RNAplex_analysis {
    my($qseq,$query_name,$txome,$qfile) = @_;  ## passed by reference
    
    # paranoia
    my $u_qseq = uc $$qseq;
    $u_qseq =~ s/T/U/g;
    
    # get the perfect MFE
        
    my $perfect_MFE = get_perfect_MFE($u_qseq);
    unless($perfect_MFE) {
	die "ABORT: Failed to get perfect MFE in sub-routine RNAplex_analysis for query sequence $u_qseq\n";
    }
    
    # Interaction length is length query + 10. 
    my $int_length = (length $u_qseq) + 10;
    
    # spacing left on default zero .. return maximal number of sites regardles of spacing. Redundant sites will be sorted out later
    
    # option -e .. 
    my $option_e = sprintf("%.2f",($opt_r * $perfect_MFE));
    
    # Call it
    (open(PLEX, "RNAplex -f 2 -e $option_e -z $int_length < $$qfile |")) || die "ABORT: Failed to open RNAplex main job in sub-routine RNAplex_analysis for query $$query_name\n";
    
    my $plex_brax;
    my $plex_mfe;
    my @local_pos = ();
    my $mfe_ratio;

    my $padded_brax;
    my @tx_pos = ();
    my $tx_ungapped_seq;
    my $ungapped;
    my @to_be_gapped;
    my @gapped = (); 
    my @out = ();
    my $outstring;
    
    my $tx;
    
    my $pline;
    
    while (<PLEX>) {
	
	# test
	#print STDERR "Working on plex line $_";
	
	chomp;
	$pline = $_;
	if($pline =~ /^>/) {
	    $pline =~ s/>//g;
	    unless($pline eq "query") {
		$tx = $pline;
	    }
	} elsif ($pline =~ /^[\.\(]+\&/) {
	    
	    ## TEST
	    #print STDERR "\tworking on bracketer pline $pline\n";
	    
	    ## all other data.
	    ## ensure that the perfect_MFE has been set
	    unless($perfect_MFE) {
		die "ABORT: Failed to set MFE of perfect match in sub-routine RNAplex_analysis\n";
	    }
	    
            # check mfe ratio. If it fails op_r, cease analysis of this alignment
	    $plex_mfe = get_plex_MFE($pline);
	    
	    #print STDERR "\tplex_mfe: $plex_mfe\n";

	    $mfe_ratio = $plex_mfe / $perfect_MFE;
	    
	    #print STDERR "\tmfe_ratio is $mfe_ratio ";
	    
	    if($mfe_ratio > 1) {
		$mfe_ratio = 1;
	    }
	    if($mfe_ratio < $opt_r) {
		next;  ## enclosing loop is the <PLEX> loop. So this goes to next alignment
	    }
	    
	    #print STDERR "  proceeding\n";
	    
	    # get the plex_brax
	    $plex_brax = get_plex_brax($pline);
	    
	    #print STDERR "\tplex_brax is $plex_brax\n";
	    
	    
	    # If you are here, then the alignment qualifies for reporting.

	    # Process the brax
	    # First, pad with dots so that the length of the query is represented
	    @local_pos = get_local_pos($pline);

	    $padded_brax = pad_plex_brax(\$plex_brax,\@local_pos,\$u_qseq);
	    
	    # After pad_plex_brax above, the @local_pos array will be modified (if needed) to reflect the local
	    #  coordinates within the transcript.

	    # Ensure that the adjusted sequence location is still valid.  If not, abort analysis of this alignment
	    if(($local_pos[0] < 1) or
	       ((length $$txome{$tx}) < $local_pos[1])) {
		next;
	    }
	    
	    $tx_ungapped_seq = substr($$txome{$tx}, ($local_pos[0] - 1), ($local_pos[1] - $local_pos[0] + 1));
	    
	    # Ensure that this sequence is all AUCG
	    unless($tx_ungapped_seq =~ /^[AUCG]+$/) {
		next;
	    }
	    
	    # Generate the initial, ungapped sequence string
	    $ungapped = "$tx_ungapped_seq" . "&" . "$u_qseq";
	    
	    # Process the alignments to remove trailing nts on the transcript
	    @to_be_gapped = no_trailing(\$padded_brax,\$ungapped,\@local_pos);
	    
	    # Gap-ify both the nt and structure alignment strings
	    @gapped = gapify($to_be_gapped[0],$to_be_gapped[1]);
	    
	    # Quality control
	    my $qc = quality_control(@gapped);
	    unless($qc) {
		print STDERR "\nWARNING: Skipping an alignment for $$query_name at transcript $tx $local_pos[0] to $local_pos[1] because of a failed alignment parse\n\t@gapped\n";
		next;
	    }
	    
	    # Begin to build output string
	    $outstring = "$$query_name"; ## reset ..  [0] Query name
	    $outstring .= "\t$tx";  ## [1] Transcript name
	    $outstring .= "\t$local_pos[0]"; ## [2] Transcript start position
	    $outstring .= "\t$local_pos[1]"; ## [3] Transcript stop position
	    
	    # Determine transcript position corresponding to position 10 of the query .. the 'slicing site'
	    my $slice_site = compute_slice_site(\@gapped,\@local_pos);
	    
	    $outstring .= "\t$slice_site";  ## [4] Transcript slice position (pos 10 of query)
	    $outstring .= "\t$perfect_MFE"; ## [5] MFE of perfect match
	    $outstring .= "\t$plex_mfe";  ## [6] MFE of this alignment
	    $outstring .= "\t$mfe_ratio"; ## [7] MFE Ratio of this alignment to perfect match
	    
	    # Calculate the Allen et al. score
	    my $allen = allen_score(@gapped);
	    $outstring .= "\t$allen";  ## [8] Allen et al. score

	    # Assess paired, unpaired, 5' SS, 3' SS, SIL, AIL, and Bulges.
	    my @struc = assess_pairing($gapped[0],$local_pos[0],$u_qseq);
	    $outstring .= "\t$struc[0]";  ## [9] Paired string
	    $outstring .= "\t$struc[1]"; ## [10] Annotated unpaired string
	    
	    # Add on the alignments themselves
	    $outstring .= "\t$gapped[0]";  ## [11] Bracket alignment
	    $outstring .= "\t$gapped[1]"; ## [12] Sequence alignment
	    
            # test
	    #print "$outstring\n";
	    #my @txb = get_al_array($gapped[0],0);
	    #my @qxb = get_al_array($gapped[0],1);
	    #my @txs = get_al_array($gapped[1],0);
	    #my @qxs = get_al_array($gapped[1],1);
	    #my $tb_pretty = join('',@txb);
	    #my $ts_pretty = join('',@txs);
	    #my $qb_pretty = reverse(join('',@qxb));
	    #my $qs_pretty = reverse(join('',@qxs));
	    
            #print "$ts_pretty\n$tb_pretty\n$qb_pretty\n$qs_pretty\n";
	    
	    # test
	    #print "$outstring\n";
	    
	    push(@out,$outstring);
	}
    }
    close PLEX;
    return @out;
}
		
sub hash_transcriptome {
    my($tx_file) = @_;
    unless($opt_q) {
	print STDERR "\nHashing the transcriptome data...";
    }
    my %hash = ();
    (open(TX, "$tx_file")) || die "ABORT: Failed to open transcripts_file $tx_file in sub-routine hash_transcriptome\n";
    my $key;
    my $seq;
    my $line;
    while (<TX>) {
	chomp;
	if($_ =~ /^>(\S+)$/) {
	    if($seq) {
		$hash{$key} = $seq;
	    }
	    $seq = '';
	    $key = $1;
	} elsif ($_ =~ /^>/) {
	    die "ABORT: Invalid transcript name $_ \n\tTranscript names must contain no whitespace! Please edit your transcriptome file\n";
	} else {
	    $line = uc $_; ## force uppercase
	    $line =~ s/T/U/g;  ## RNA-i-fy
	    $seq .= $line;
	}
    }
    close TX;
    # enter last one
    $hash{$key} = $seq;
    unless($opt_q) {
	print STDERR " Done\n";
    }
    return %hash;
}
	   
sub get_plex_brax {
    my($line) = @_;
    my $brax;
    if($line =~ /^([\(\.\)\&]+)\s/) {
	$brax = $1;
    }
    if($brax) {
	return $brax;
    } else {
	die "ABORT: Failed to parse out RNAplex structure from input line $line in sub-routine get_plex_brax\n";
    }
}

sub get_plex_MFE {
    my($line) = @_;
    my $mfe;
    if($line =~ /\(([\s\-\d\.]+)\)/){
	$mfe = $1;
    }
    if($mfe) {
	return $mfe;
    } else {
	#print STDERR "\nWARNING: Failed to parse out MFE from input line $line in sub-routine get_plex_MFE .. reported as 0\n";
	return 0;
    }
}

sub get_local_pos {
    my($line) = @_;
    my @local_pos = ();
    if($line =~ /(\d+),(\d+).*:.*(\d+),(\d+)/) {
	push(@local_pos,$1);
	push(@local_pos,$2);
	push(@local_pos,$3);
	push(@local_pos,$4);
    }
    if(@local_pos) {
	return @local_pos;
    } else {
	die "ABORT: Failed to parse out local coordinates from input line $line in sub-routine get_local_pos\n";
    }
}

sub pad_plex_brax {
    my($inbrax,$pos,$qseq) = @_; # passed by reference .. scalar, array, scalar
    
    # test
    #print "\nInput: $$inbrax @$pos\n";
    
    my $q_len = length $$qseq;
    my $t_part;
    my $q_part;
    if($$inbrax =~ /^([\(\.]+)\&([\)\.]+)$/) {
	$t_part = $1;
	$q_part = $2;
    } else {
	die "ABORT: Failed to split the incoming brax in sub-routine pad_plex_brax\n";
    }
    
    my $i;
    my $left_added = 0;  ## left of query, right of transcript
    for($i = $$pos[2]; $i > 1; --$i) {
	$q_part = "." . "$q_part";
	$t_part = "$t_part" . ".";
	++$left_added;
    }
    my $right_added = 0; ## right of query, left of transcript
    for($i = $$pos[3]; $i < $q_len; ++$i) {
	$q_part = "$q_part" . ".";
	$t_part = "." . "$t_part";
	++$right_added;
    }

    # modify the local position array
    $$pos[0] -= $right_added;
    $$pos[1] += $left_added;
    $$pos[2] -= $left_added;
    $$pos[3] += $right_added;
    
    
    
    
    # Ensure against rare cases when the transcript does not fully extend to the limits of the query, likely because
    #  the orignial longer transcript sequence did not fully capture the site.
    my $t_5_dots = 0;
    my $t_3_dots = 0;
    my $q_5_dots = 0;
    my $q_3_dots = 0;
    if($t_part =~ /^\.+/) {
	$t_5_dots = length $&;
    }
    if($t_part =~ /\.+$/) {
	$t_3_dots = length $&;
    }
    if($q_part =~ /^\.+/) {
	$q_5_dots = length $&;
    }
    if($q_part =~ /\.+$/) {
	$q_3_dots = length $&;
    }
    
    # Add dots to t_5 if needed
    for($i = ($q_3_dots - $t_5_dots); $i > 0; --$i) {
	$t_part = "." . "$t_part";
	--$$pos[0];
    }
    # Add dots to t_3 if needed
    for($i = ($q_5_dots - $t_3_dots); $i > 0; --$i) {
	$t_part .= ".";
	++$$pos[1];
    }

    # rejoin
    my $outbrax = "$t_part" . "&" . "$q_part";
    
    return $outbrax;
    
}

sub get_tx_final {
    my($lpos,$old_tx) = @_; ## passed by reference .. array and scalar
    my $old_start;
    my $tx_key;
    if($$old_tx =~ /^(\S+):(\d+)-\d+$/) {
	$tx_key = $1;
	$old_start = $2;
    } else {
	die "ABORT: Failed to parse transcript location string $$old_tx in sub-routine get_tx_final\n";
    }
    my @out = ();
    $out[0] = $tx_key;
    $out[1] = $old_start + $$lpos[0] - 1;
    $out[2] = $old_start + $$lpos[1] - 1;
    
    return @out;
}

sub gapify {
    my($inbrax,$inseq) = @_;

    # test
    #print "\ngapify inbrax: $inbrax inseq: $inseq\n";
    #exit;

    my @out = ();
    my %left_right = get_left_right($inbrax);
    my $last_left = 0;
    my $last_right = (length $inbrax) + 1;
    my @inbrax_chars = split ('',$inbrax);
    my @inseq_chars = split ('', $inseq);
    my $left_delta;
    my $right_delta;
    my $delta;
    my $index;
    my $outbrax_left;
    my $outbrax_right;
    my $outseq_left;
    my $outseq_right;
    my $x;
    for (my $i = 1; $i <= (length $inbrax); ++$i) {
	if(exists($left_right{$i})) {
	    
	    $left_delta = $i - $last_left + 1;
	    $right_delta = $last_right - $left_right{$i} + 1;
	    $delta = $left_delta - $right_delta;
	    
	    # add trailing characters, should be all dots
	    if($left_delta) {
		for($x = $last_left; $x < ($i - 1); ++$x) {
		    $outbrax_left .= "$inbrax_chars[$x]";
		    $outseq_left .= "$inseq_chars[$x]";
		}
	    }
	    if($right_delta) {
		for($x = ($last_right - 2); $x >= $left_right{$i}; --$x) {
		    if($outbrax_right) {
			$outbrax_right = "$inbrax_chars[$x]" . "$outbrax_right";
			$outseq_right = "$inseq_chars[$x]" . "$outseq_right";
		    } else {
			$outbrax_right = "$inbrax_chars[$x]";
			$outseq_right = "$inseq_chars[$x]";
		    }
		}
	    }
	    # add gaps if need be
	    if($delta > 0) {
		# gap(s) added on the right side
		for($x = $delta; $x > 0; --$x) {
		    if($outbrax_right) {
			$outbrax_right = "-" . "$outbrax_right";
			$outseq_right = "-" . "$outseq_right";
		    } else {
			$outbrax_right = "-";
			$outseq_right = "-";
		    }
		}
	    } elsif ($delta < 0) {
		# gap(s) added on the left side
		for($x = $delta; $x < 0; ++$x) {
		    $outbrax_left .= "-";
		    $outseq_left .= "-";
		}
	    }
	    
	    # Add the current pair
	    $outbrax_left .= "$inbrax_chars[($i - 1)]";
	    $outseq_left .= "$inseq_chars[($i - 1)]";
	    if($outbrax_right) {
		$outbrax_right = "$inbrax_chars[($left_right{$i} - 1)]" . "$outbrax_right";
		$outseq_right = "$inseq_chars[($left_right{$i} - 1)]" . "$outseq_right";
	    } else {
		$outbrax_right = "$inbrax_chars[($left_right{$i} - 1)]";
		$outseq_right = "$inseq_chars[($left_right{$i} - 1)]" ;
	    }
	    
	    # reset
	    $last_left = $i;
	    $last_right = $left_right{$i};
		
	} elsif ($inbrax_chars[($i - 1)] eq "\&") {
	    # Have reach the middle. Just report any unreported dots since the last pair
	    $left_delta = $i - $last_left + 1;
	    $right_delta = $last_right - ($i + 1);
	    
	    # add trailing characters, should be all dots
	    if($left_delta) {
		for($x = $last_left; $x < ($i - 1); ++$x) {
		    $outbrax_left .= "$inbrax_chars[$x]";
		    $outseq_left .= "$inseq_chars[$x]";
		}
	    }
	    if($right_delta) {
		for($x = ($last_right - 2); $x >= $i; --$x) {
		    $outbrax_right = "$inbrax_chars[$x]" . "$outbrax_right";
		    $outseq_right = "$inseq_chars[$x]" . "$outseq_right";
		}
	    }
	    last;
	}
    }

    $out[0] = "$outbrax_left" . "\&" . "$outbrax_right";
    $out[1] = "$outseq_left" . "\&" . "$outseq_right";
    
    # TEST
    #print "\nINPUT:\n", "\t$inbrax\n", "\t$inseq\n";
    #print "OUTPUT:\n", "\t$out[0]\n", "\t$out[1]\n";
    #my $rev_right_brax = reverse $outbrax_right;
    #my $rev_right_seq = reverse $outseq_right;
    #print "$outseq_left\n$outbrax_left\n$rev_right_brax\n$rev_right_seq\n";
    
    return @out;
    
}

sub get_left_right {
    my($brax)  = @_;
    my %hash = ();
    my @chars = split('',$brax);
    ## Note that the & is included in the above array
    my $i = 0;
    my @lefts = ();
    my $left;
    foreach my $char (@chars) {
	++$i;
	if($char eq "\(") {
	    push(@lefts, $i);
	} elsif ($char eq "\)") {
	    $left = pop @lefts;
	    $hash{$left} = $i;
	}
    }
    return %hash;
}

sub no_trailing {
    my($brax,$ungapped,$pos) = @_;  ## passed by reference .. scalar, scalar, array
    my @out = ();
    my @bs = split("\&",$$brax);
    my @ss = split("\&",$$ungapped);
    my $n_left_end_dots = 0;
    my $n_right_end_dots = 0;
    my $n_left_mid_dots = 0;
    my $n_right_mid_dots = 0;
    
    my $n_end_trim = 0;
    my $n_mid_trim = 0;
 
    ## test
    #print "\nno_trailing INPUT $$brax $$ungapped @$pos\n";
    

    if($bs[0] =~ /^\.+/) {
	$n_left_end_dots += (length $&);
    }
    if($bs[0] =~ /\.+$/) {
	$n_left_mid_dots += (length $&);
    }
    if($bs[1] =~ /^\.+/) {
	$n_right_mid_dots += (length $&);
    }
    if($bs[1] =~ /\.+$/) {
	$n_right_end_dots += (length $&);
    }
    
    if($n_left_end_dots > $n_right_end_dots) {
	$n_end_trim = $n_left_end_dots - $n_right_end_dots;
    }
    if($n_left_mid_dots > $n_right_mid_dots) {
	$n_mid_trim = $n_left_mid_dots - $n_right_mid_dots;
    }
    
    my $junk;
    my @bleft = split ('', $bs[0]);
    my @sleft = split ('', $ss[0]);
    my $i;
    for($i = 0; $i < $n_end_trim; ++$i) {
	$junk = shift @bleft;
	$junk = shift @sleft;
	++$$pos[0];
    }
    for($i = 0; $i < $n_mid_trim; ++$i) {
	$junk = pop @bleft;
	$junk = pop @sleft;
	--$$pos[1];
    }
    
    my $new_brax_left = join('',@bleft);
    my $new_seq_left = join('',@sleft);
    $out[0] = "$new_brax_left" . "\&" . "$bs[1]";
    $out[1] = "$new_seq_left" . "\&" . "$ss[1]";
    
    ## test
    #print "no_trailing OUTPUT $out[0] $out[1] @$pos\n";
    
    #exit;
    
    return @out;
}
	
sub compute_slice_site {
    my($gapped,$local_pos) = @_;  ## passed by reference.  Two arrays.
    my $q_seq;
    my $t_seq;
    if($$gapped[0] =~ /^(\S+)\&(\S+)$/) {
	$t_seq = $1;
	$q_seq = $2;
    } else {
	die "ABORT: In sub-routine compute_slice_site, failed to parse alignment $$gapped[0]\n";
    }
    my @q_char = split ('', $q_seq);
    my @t_char = split ('', $t_seq);
    
    my $real_t_pos = $$local_pos[1] + 1;
    my $real_q_pos = 0;
    my $tch;
    
    foreach my $qch (@q_char) {
	$tch = pop @t_char;
	unless($tch eq "-") {
	    --$real_t_pos;
	}
	unless($qch eq "-") {
	    ++$real_q_pos;
	}
	if($real_q_pos == 10) {
	    return $real_t_pos;
	}
    }
    # Should not get here
    die "ABORT: Failure in sub-routine compute_slice_site\n";
}

sub allen_score {
    my(@al) = @_;
    my $score = 0;
    
    my @tbrax = get_al_array($al[0],0);
    my @qbrax = get_al_array($al[0],1);
    my @tseq = get_al_array($al[1],0);
    my @qseq = get_al_array($al[1],1);
    
    my $q_pos = 0;
    my $ts;
    my $qb;
    my $tb;
    my $i = -1;

    #print "INPUT: @al\n";


    foreach my $qs (@qseq) {
	++$i;
	unless($qs eq "-") {
	    ++$q_pos;
	}
	$qb = $qbrax[$i];
	$ts = pop @tseq;
	$tb = pop @tbrax;
	
	#print "q_pos: $q_pos ts $ts qs $qs tb $tb qb $qb\n";
	
	# trap error
	#unless($tb) {
	    #print "\nError trap in sub-routine allen_acore INPUT: @al\n";
	    #print "\ttbrax: @tbrax\n";
	    #print "\tqbrax: @qbrax\n";
	    #print "\ttseq: @tseq\n";
	    #print "\tqseq: @qseq\n";
	    #exit;
	#}
	
	if(($tb eq ".") or ($qb eq ".")) {  ## acounts for bulges, assymmetric internal loops as well
	    $score += 1;
	    #print "\tpenalty 1\n";
	    if(($q_pos >= 2) and ($q_pos <= 12)) {
		#print "\tpenalty doubled\n";
		$score += 1;
	    }
	} elsif ((($ts eq "G") and ($qs eq "U") and ($tb eq "\(") and ($qb eq "\)")) or
		 (($ts eq "U") and ($qs eq "G") and ($tb eq "\(") and ($qb eq "\)"))) {
	    $score += 0.5;
	    #print "\tGU wobble penalty 0.5\n";
	    if(($q_pos >= 2) and ($q_pos <= 12)) {
		#print "\tGU woblle penalty doubled\n";
		$score += 0.5;
	    }
	}
    }
    #print "final score: $score\n";
    #exit;
    return $score;
}

sub get_al_array {
    my($string,$index) = @_;
    my @two = split ("\&", $string);
    my $chunk = $two[$index];
    my @out = split ('', $chunk);
    return @out;
}
    
sub assess_pairing {
    my($brax,$tx_start,$u_qseq) = @_;
    my @out = ();
    $out[0] = check_pairs($brax,$tx_start,$u_qseq);
    my $unpaired = check_unpaired($brax,$tx_start,$u_qseq);
    if($unpaired) {
	$out[1] = annotate_unpaired($unpaired,$tx_start,$u_qseq);
    } else {
	$out[1] = "NA"; ## All are paired
    }

    # test
    #print "sub assess_pairing INPUT $brax $tx_start $u_qseq\n";
    #print "\tpairs: $out[0]\n";
    #print "\tUNpairs: $out[1]\n";
    
    return @out;
}

sub check_pairs {
    my($brax,$tx_start,$u_qseq) = @_;
    my $out;
    
    my @tbrax = get_al_array($brax,0);
    my @qbrax = reverse(get_al_array($brax,1));
    
    my $i = -1;
    my $t_pos = $tx_start - 1;

    my $q_pos = (length $u_qseq) + 1;
    
    my $last_t_start;
    my $last_t_paired;
    my $last_q_start;
    my $last_q_paired;
    
    my $q;
    
    my $output;
    
    foreach my $t (@tbrax) {
	++$i;
	$q = $qbrax[$i];
	
	# check for closure
	if(($t ne "\(") or ($q ne "\)")) {
	    if($last_t_start) {
		if($output) {
		    $output = "$last_q_paired" . "-" . "$last_q_start" . "," . "$last_t_paired" . "-" . "$last_t_start" . "\;" . "$output";
		} else {
		    $output = "$last_q_paired" . "-" . "$last_q_start" . "," . "$last_t_paired" . "-" . "$last_t_start";
		}
	    }
	    # reset
	    $last_t_start = '';
	    $last_t_paired = '';
	    $last_q_start = '';
	    $last_q_paired = '';
	}
	
	# Update positions
	unless($t eq "-") {
	    ++$t_pos;
	}
	unless($q eq "-") {
	    --$q_pos;
	}
	
	# If the current position is a pair, track it
	if(($t eq "\(") and ($q eq "\)")) {
	    $last_t_paired = $t_pos;
	    $last_q_paired = $q_pos;
	    unless($last_t_start) {
		$last_t_start = $t_pos;
	    }
	    unless($last_q_start) {
		$last_q_start = $q_pos;
	    }
	}
    }
    # check for closure at end
    if($last_t_start) {
	if($output) {
	    $output = "$last_q_paired" . "-" . "$last_q_start" . "," . "$last_t_paired" . "-" . "$last_t_start" . "\;" . "$output";
	} else {
	    $output = "$last_q_paired" . "-" . "$last_q_start" . "," . "$last_t_paired" . "-" . "$last_t_start";
	}
    }
    # reset
    $last_t_start = '';
    $last_t_paired = '';
    $last_q_start = '';
    $last_q_paired = '';
    
    return $output;
}
	    
sub check_unpaired {
    my($brax,$tx_start,$u_qseq) = @_;
    my $out;
    
    my @tbrax = get_al_array($brax,0);
    my @qbrax = reverse(get_al_array($brax,1));
    
    my $i = -1;
    my $t_pos = $tx_start - 1;

    my $q_pos = (length $u_qseq) + 1;
    
    my $last_t_start;
    my $last_t_unpaired;
    my $last_q_start;
    my $last_q_unpaired;
    
    my $q;
    
    my $output;
    
    foreach my $t (@tbrax) {
	++$i;
	$q = $qbrax[$i];
	
	# check for closure
	if(($t eq "\(") and ($q eq "\)")) {
	    # test
	    #print "Hello world\n";
	    
	    if($last_t_start) {
		
		# test
		#print "In reporting loop 1\n";
		
		if($output) {
		    $output = "$last_q_unpaired" . "-" . "$last_q_start" . "," . "$last_t_unpaired" . "-" . "$last_t_start" . "\;" . "$output";
		} else {
		    $output = "$last_q_unpaired" . "-" . "$last_q_start" . "," . "$last_t_unpaired" . "-" . "$last_t_start";
		}
	    }
	    # reset
	    $last_t_start = '';
	    $last_t_unpaired = '';
	    $last_q_start = '';
	    $last_q_unpaired = '';
	}
	
	# Update positions
	unless($t eq "-") {
	    ++$t_pos;
	}
	unless($q eq "-") {
	    --$q_pos;
	}
	
	# If the current position is not a pair, track it
	if(($t ne "\(") or ($q ne "\)")) {
	    if($t eq "-") {
		if($last_t_unpaired) {
		    unless($last_t_unpaired =~ /\d+/) {
			$last_t_unpaired = "x";
		    }
		} else {
		    $last_t_unpaired = "x";
		}
		unless($last_t_start) {
		    $last_t_start = "x"; 
		}
	    } else {
		$last_t_unpaired = $t_pos;
		if($last_t_start) {
		    if($last_t_start eq "x") {
			$last_t_start = $t_pos;
		    }
		} else {
		    $last_t_start = $t_pos;
		}
	    }

	    if($q eq "-") {
		if($last_q_unpaired) {
		    unless($last_q_unpaired =~ /\d+/) {
			$last_q_unpaired = "x";
		    }
		} else {
		    $last_q_unpaired = "x";
		}
		unless($last_q_start) {
		    $last_q_start = "x"; 
		}
	    } else {
		$last_q_unpaired = $q_pos;
		if($last_q_start) {
		    if($last_q_start eq "x") {
			$last_q_start = $q_pos;
		    }
		} else {
		    $last_q_start = $q_pos;
		}
	    }
	}
    }
    # check for closure at end
    if($last_t_start) {
	# test
	#print "In reporting loop 2\n";
	if($output) {
	    $output = "$last_q_unpaired" . "-" . "$last_q_start" . "," . "$last_t_unpaired" . "-" . "$last_t_start" . "\;" . "$output";
	} else {
	    $output = "$last_q_unpaired" . "-" . "$last_q_start" . "," . "$last_t_unpaired" . "-" . "$last_t_start";
	}
    }
    # reset
    $last_t_start = '';
    $last_t_unpaired = '';
    $last_q_start = '';
    $last_q_unpaired = '';
    
    return $output;
}

sub annotate_unpaired {
    my($in,$tx_start,$qseq) = @_;
    my @unp = split ("\;", $in);
    my $out;
    my $q_start;
    my $q_stop;
    my $t_start;
    my $t_stop;
    my $max_q = length $qseq;
    foreach my $entry (@unp) {
	if($out) {
	    $out .= "\;" . "$entry";
	} else {
	    $out = $entry;
	}
	
	if($entry =~ /^(\S+)-(\S+)\,(\S+)-(\S+)$/) {
	    $q_start = $1;
	    $q_stop = $2;
	    $t_start = $4;
	    $t_stop = $3;
	} else {
	    die "ABORT: Failed to parse entry $entry in sub-routine annotate_unpaired\n";
	}
	
	# check for any x's first
	if(($q_start =~ /^\d+$/) and
	   ($q_stop =~ /^\d+$/) and
	   ($t_start =~ /^\d+$/) and
	   ($t_stop =~ /^\d+$/)) {
	    
	    ## All numeric .. cannot be a bulge
	    my $q_delta = $q_stop - $q_start + 1;
	    my $t_delta = $t_stop - $t_start + 1;
	    
	    if($q_start == 1) {
		# UP5
		$out .= "\[UP5\]";
	    } elsif ($q_stop == $max_q) {
		# UP3
		$out .= "\[UP3\]";
	    } elsif ($q_delta == $t_delta) {
		# SIL
		$out .= "\[SIL\]";
	    } elsif ($q_delta > $t_delta) {
		# AILq
		$out .= "\[AILq\]";
	    } elsif ($q_delta < $t_delta) {
		# AILt
		$out .= "\[AILt\]";
	    }
	} else {
	    ## Somebody is non-numeric. Should be either query with x-x or transcript with x-x. In either case a bulge
	    ## Verify
	    if(($q_start =~ /^x$/) and ($q_stop =~ /^x$/)) {
		## BULt
		$out .= "\[BULt\]";
	    } elsif (($t_start =~ /^x$/) and ($t_stop =~ /^x$/)) {
		## BULq
		$out .= "\[BULq\]";
	    } else {
		## not suppressed by quiet mode because it is a warning
		print STDERR "\n\tWARNING: Failed to parse entry $entry in sub-routine annotate_unpaired .. reporting type \?\n";
		$out .= "\[\?\]";
	    }
	}
    }
    return $out;
}

sub sort_first_pass {
    my(@in) = @_;
    my %sort = ();
    my %all = ();
    my @fields = ();
    my $key;
    foreach my $entry (@in) {
	@fields = split ("\t", $entry);
	$key = "$fields[1]" . ":" . "$fields[4]";  ## slice site
	if(exists($sort{$key})) {
	    if($opt_a) {
		if($fields[8] < $sort{$key}) {
		    # keep the new one
		    $sort{$key} = $fields[8];
		    $all{$key} = $entry;
		}
	    } else {
		if($fields[7] > $sort{$key}) {
		    # keep the new one
		    $sort{$key} = $fields[7];
		    $all{$key} = $entry;
		}
	    }
	} else {
	    if($opt_a) {
		$sort{$key} = $fields[8];
		$all{$key} = $entry;
	    } else {
		$sort{$key} = $fields[7];
		$all{$key} = $entry;
	    }
	}
    }
    
    # get an array of keys sorted by value in the sort hash
    my @sorted = ();
    if($opt_a) {
	# Want to sort ascending if sorting by Allen et al. scores
	@sorted = sort{$sort{$a} <=> $sort{$b}} keys(%sort);
    } else {
	# Otherwise want to sort descending for MFEratio
	@sorted = sort{$sort{$b} <=> $sort{$a}} keys(%sort);
    }
    
    # and to finish
    my @out = ();
    foreach my $skey (@sorted) {
	push(@out, $all{$skey});
    }
    return @out;
}
    
sub quality_control {
    my(@gapped) = @_;
    my @tbrax = get_al_array($gapped[0],0);
    my @qbrax = get_al_array($gapped[0],1);
    my @tseq = get_al_array($gapped[1],0);
    my @qseq = get_al_array($gapped[1],1);
    
    my $fail;

    # All four alignment subsections need to have the same size after gapping

    unless(((scalar @tbrax) == (scalar @qbrax)) and
	   ((scalar @tbrax) == (scalar @tseq)) and
	   ((scalar @tbrax) == (scalar @qseq))) {
	$fail = 1;
	
	# test
	#print STDERR "QC Failed equivalent array size check for @gapped\n";
	#print STDERR "tbrax @tbrax\n";
	#print STDERR "qbrax @qbrax\n";
	#print STDERR "tseq  @tseq\n";
	#print STDERR "qseq  @qseq\n";
	
    }
    
    # No unexpected characters in sequences or brackets
    my $char;
    my $t_pairs = 0;
    my $q_pairs = 0;
    foreach $char (@tseq) {
	unless($char =~ /^[AUCG-]$/) {
	    $fail = 1;
	    
	    # test
	    #print STDERR "QC Failed character test for tseq entry $char from array @tseq\n";
	    
	}
    }
    foreach $char (@qseq) {
	unless($char =~ /^[AUCG-]$/) {
	    $fail = 1;
	    
	    # test
	    #print STDERR "QC Failed character test for qseq entry $char from array @qseq\n";
	    
	}
    }
    foreach $char (@tbrax) {
	unless($char =~ /^[\(\.-]$/) {
	    $fail = 1;
	    
	    #test
	    #print STDERR "QC Failed character test for tbrax entry $char from array @tbrax\n";
	    
	}
	if($char eq "\(") {
	    ++$t_pairs;
	}
    }
    foreach $char (@qbrax) {
	unless($char =~ /^[\)\.-]$/) {
	    $fail = 1;
	    
	    #test
	    #print STDERR "QC Failed character test for qbrax entry $char from array @qbrax\n";
	    
	}
	if($char eq "\)") {
	    ++$q_pairs;
	}
    }
    
    # Number of pairs must be the same in tbrax and qbrax
    unless($t_pairs == $q_pairs) {
	$fail = 1;
	
	# test
	#print STDERR "QC Failed equal pairs test with t of $t_pairs and q of $q_pairs for entries @tbrax and @qbrax\n";
	
    }
    
    if($fail) {
	return 0;
    } else {
	return 1;
    }
}
    
sub tabular_output {
    my(@out) = @_;
    foreach my $string (@out) {
	print "$string\n";
    }
}

sub pretty_output {
    my(@out) = @_;
    my @fields = ();
    my $i;
    print "----------------------------------------------------------------\n";
    foreach my $string (@out) {
	@fields = split ("\t", $string);
	print "5\' ";
	my @tseq = get_al_array($fields[12],0);
	print (@tseq," 3\'");
	print " Transcript: $fields[1]",":","$fields[2]","-","$fields[3]"," Slice Site:$fields[4]\n";
	my @tbrax = get_al_array($fields[11],0);
	my @qbrax = reverse(get_al_array($fields[11],1));
	my @qseq = reverse(get_al_array($fields[12],1));
	$i = -1;
	print "   ";
	foreach my $tb (@tbrax) {
	    ++$i;
	    my $qb = $qbrax[$i];
	    my $ts = $tseq[$i];
	    my $qs = $qseq[$i];
	    if(($tb eq "\(") and ($qb eq "\)")) {
		# is it a G-U wobble?
		if((($ts eq "G") and ($qs eq "U")) or
		   (($ts eq "U") and ($qs eq "G"))) {
		    print "o";
		} else {
		    print "\|";
		}
	    } else {
		print " ";
	    }
	}
	print "\n";
	print "3\' ";
	print (@qseq, " 5\'");
	print " Query: $fields[0]\n";
	print "MFE of perfect match: $fields[5]\n";
	print "MFE of this site: $fields[6]\n";
	print "MFEratio: $fields[7]\n";
	print "Allen et al. score: $fields[8]\n";
	my @paired = split ("\;", $fields[9]);
	print "Paired Regions \(query5\'-query3\',transcript3\'-transcript5\'\)\n";
	foreach my $pair (@paired) {
	    print "\t$pair\n";
	}
	print "Unpaired Regions \(query5\'-query3\',transcript3\'-transcript5\'\)\n";
	my $coord;
	my $unp_type;
	if($fields[10] eq "NA") {
	    print "\tNA\n";
	} else {
	    my @unpaired = split ("\;", $fields[10]);
	    foreach my $unp (@unpaired) {
		if($unp =~ /^(\S+-\S+,\S+-\S+)\[(\S+)\]$/) {
		    print "\t$1\t$2: ";
		    if($2 eq "UP5") {
			print "Unpaired region at 5\' of query\n";
		    } elsif ($2 eq "UP3") {
			print "Unpaired region at 3\' of query\n";
		    } elsif ($2 eq "SIL") {
			print "Symmetric internal loop\n";
		    } elsif ($2 eq "AILq") {
			print "Asymmetric internal loop with more nts on the query side\n";
		    } elsif ($2 eq "AILt") {
			print "Asymmetric internal loop with more nts on the transcript side\n";
		    } elsif ($2 eq "BULt") {
			print "Bulge on transcript side\n";
		    } elsif ($2 eq "BULq") {
			print "Bulge on query side\n";
		    } else {
			print "Unknown - Likely a bug in the parsing of this alignment\n";
		    }
		} else {
		    die "ABORT: in sub-routine pretty_output : failed to parse unpaired string $unp\n";
		}
	    }
	}
	print "----------------------------------------------------------------\n";
    }
}

sub write_q {  ## passed by reference
    my($qseq,$txome) = @_;
    my $q_file = "GSTAr_temp.fasta";
    (open(Q, ">$q_file")) || die "ABORT: Failed to open temp file for writing in sub-routine write_q\n";
    
    my $t_seq;
    my $t_name;
    while (($t_name,$t_seq) = each %$txome) {
	print Q ">$t_name\n$t_seq\n>query\n$$qseq\n";
    }
    close Q;
    ## test

    return $q_file;
}

sub get_perfect_MFE {
    my($seq) = @_;
    my $perfect = reverse $seq;
    $perfect =~ tr/AUCG/UAGC/;
    (open(TEMP, ">GSTAr_perfect_tmp.fasta")) || die "ABORT: Failed to create a temp file in sub-routine get_perfect_MFE\n";
    print TEMP ">perfect\n$perfect\n>Query\n$seq\n";
    close TEMP;
    
    # Call RNAplex under default parameters
    (open(PLEX, "RNAplex < GSTAr_perfect_tmp.fasta |")) || die "ABORT: RNAplex call failed in sub-routine get_perfect_MFE\n";
    
    my $perfect_MFE;
    
    while (<PLEX>) {
	chomp;
	if($_ =~ /^[\.\(]+\&/) {
	    $perfect_MFE = get_plex_MFE($_);
	    last;
	}
    }
    close PLEX;
    #print "PERFECT: $perfect_MFE\n";
    #exit;
    system "rm -f GSTAr_perfect_tmp.fasta";
    return $perfect_MFE;
}
	    
	
__END__

=head1 LICENSE

GSTAr.pl

Copyright (c) 2013 Michael J. Axtell

This program is free software: you can redistribute it and/or modify                                                                                                                
it under the terms of the GNU General Public License as published by                                                                                                                
the Free Software Foundation, either version 3 of the License, or                                                                                                                   
(at your option) any later version.                                                                                                                                                 
                                                                                                                                                                                    
This program is distributed in the hope that it will be useful,                                                                                                                     
but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                                                                      
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                                                                                       
GNU General Public License for more details.                                                                                                                                        
                                                                                                                                                                                    
You should have received a copy of the GNU General Public License                                                                                                                   
along with this program.  If not, see <http://www.gnu.org/licenses/>. 

=head1 SYNOPSIS

GSTAr : Generic Small RNA-Transcriptome Aligner.

Flexible RNAplex-based alignment of miRNAs and siRNAs (15-26 nts) to a transcriptome.

=head1 AUTHOR

Michael J. Axtell, Penn State University, mja18@psu.edu

=head1 VERSION

1.0 : September 17, 2013

=head1 INSTALL

=head2 Dependencies

	perl
	RNAplex (from Vienna RNA package)

RNAplex must be executable from your PATH. GSTAr.pl was developed using RNAplex version RNAplex 2.1.3. It has not been tested on other versoins of RNAplex.

=head2 Installation

No "real" installation. If the script is in your working directory, you can call it with 

	./GSTAr.pl

For convenience, you can add it to your PATH.  e.g.

	sudo mv GSTAr.pl /usr/bin/

GSTAr.pl expects to find perl in /usr/bin/perl .. if not, edit line 1 (the hashbang) accordingly.

=head1 USAGE

	GSTAr.pl [options] queries.fasta transcriptome.fasta

Output alignments go to STDOUT, and can be redirected to a file with > or piped to another process with |

Log and progress information goes to STDERR, and can be suppressed with option -q (quiet mode).

=head2 Options

Options:  
                                                                                                                                                                          
-h Print help message and quit
                                                                                                                                                    
-v Print version and quit
                                                                                                                                                           
-q Quiet mode .. no log/progress information to STDERR 
                                                                                                                             
-t Tabular output format ... More suitable for parsing   
                                                                                                                           
-a Sort by Allen et al. score instead of the default MFEratio 
                                                                                            
-r [float >0..1] Minimum Free Energy Ratio cutoff. Default: 0.65

=head1 METHODS

GSTAr.pl is essentially a wrapper and parser for RNAplex desgined for aligning short queries (15-26nts) against the rev-comp. strand of a eukaryotic transcriptome. For each query sequence, the minimum free energy of a perfectly complelementary sequence is calculated using RNAplex under all default parameters. Following this, the same query is then analyzed against the entire transcriptome. Hits where the MFEratio (i.e. MFE / MFE-perfect) is >= the cutoff established by option r are retained and parsed.

The detailed RNAplex parameters (see RNAplex man page for details) for each query analysis are 

	-f 2 : Fast mode .. structure based on approximated plex model.
	-e [minMFEratio * perfectMFE] : Minimum acceptable MFE value to keep a hit
	-z : [10 + query_length] : Acceptable alignments can span no more than length of query + 10nts.
	

Slice Site is the transcript nt opposite nt 10 of the query. This is where one should look to find evidence of AGO-catalyzed slicing in the event that a) the transcript was really a target of the query at that site, and b) it really was sliced. GSTAr makes no judgements on the likelihood of either of those events, and the recording of the Slice Site position should NOT be taken as evidence that slicing exists or is even possible at that site.

For a given query, all returned sites are non-redundant. Redunancy is based upon the putative slice site only. By default, the output is sorted in descending order (best to worst) according to the MFE ratio. In the alternative option -a mode, the results are instead sorted in ascending order (best to worst) according to the Allen et al. score.

Allen et al. score: This is a score for plant miRNA/siRNA-target interactions based on the position-specific penalties described by Allen et al. (2005) Cell, 121:207-221 [PMID: 15851028].  Specifically, mismatched query bases or target-bulged bases, are penalized 1. G-U wobbles are penalized 0.5. These penalties are double within positions 2-13 of the query.

=head1 WARNINGS

=head2 NOT a target predictor

GSTAr is very explicitly NOT a target predictor for miRNAs or siRNAs. It is only an aligner based on RNA-RNA hybridization thermodynamic predictions. Users should make no claims as to whether the identified alignments are actually targets of the query without independent data of some sort. 

=head2 Slice Sites are NOT predictions of slicing

Although GSTAr reports a "Slicing Site" position for each alignment, this is merely for conveneince when using GSTAr alignments to guide subsequent experiments searching for AGO-catalyzed slicing evidence. No claim is made that any alignment is actually AGO-cleaved or even theoretically AGO-cleavable.

=head2 Not for whole genomes

GSTAr holds the entire contents of the transcripts.fasta file in memory to speed the isolation of sub-sequences. This will be impractical in terms of memory usage if a user attempts to load a whole genome.  Similarly, GSTAr will only search for pairing between the top strand of the transcripts.fasta file, making it also impractical for a genome analysis, where sites might be on either strand.

=head2 Temp files

GSTAr writes temp files to the working directory. Their contents change dynamically during a run, and they will be deleted at the end of a run. So, don't mess with them during a run. In addition, it is a very bad idea to have two GSTAr runs operating concurrently from the same working directory because there will be clashes and overwrites for these temp files.

=head2 Not too fast

GSTAr uses RNAplex (Tafer and Hofacker, 2008. Bioinformatics 24:2657-63, PMID: 18434344, doi:10.1093/bioinformatics/btn193), which is exceptionally fast for an inter-molecular RNA-RNA hybridization calculator. However, when applied to entire eukaryotic transcriptomes the CPU time per query is still significant. Run time is only slightly affected (much less than 2-fold) by the setting of -r. Setting tabular mode (option -t) also increases speed just a tiny bit for runs with a low option -r. In tests with the Arabidopsis transcriptome (33,602 mRNAs, total nts=51,074,197), a single 21nt miRNA query typically takes about 90-110 seconds to complete.

=head2 No ambiguity codes

Query sequences with characters other than A, T, U, C, or G (case-insensitive) will not be analyzed, and a warning will be sent to the user. Transcript sub-sequences for potential alignments will be *silently* ignored if they contain any characters other than A, T, U, C,or G (case-insensitive).

=head2 Small queries

Query sequences must be small (between 15 and 26nts). Queries that don't meet these size requirements will not be analyzed and a warning sent to the user.

=head2 Redundant output

GSTAr guarantees that, FOR A GIVEN QUERY, the returned alignments will be unique in terms of their PUTATIVE SLICING SITE POSITION. However, the same query could generate multiple overlapping aligments that each have different putative slice sites.  Furthermore, if different queries in a multi-query analysis have similar (or identical !) sequences, the same alignment position (based on putative slicing site position) could be returned multiple times, once for each of the similar/identical queries. Therefore, for multi-query result files, there is no guarantee of non-redundancy among the returned sites.

=head1 OUTPUT

Both output formats begin with a series of commented lines that provide details about the run.

=head2 Pretty output

The default output "pretty" mode is easily parsed by humans and should be self-explanatory

=head2 Tabular output

Tabular output (which occurs when option -t is specified) is a tab-delimited text file. The first non-commented line is the column headings, with meanings as follows:

1: Query: Name of query

2: Transcript: Name of transcript

3: TStart: One-based start position of the alignment within the transcript

4: TStop: One-based stop position of the alignment within the transcript

5: TSLice: One-based position of the alignment opposite position 10 of the query

6: MFEperfect: Minimum free energy of a perfectly matched site (approximate)

7: MFEsite: Minimum free energy of the alignment in question

8: MFEratio: MFEsite / MFEperfect

9: AllenScore: Penalty score calculated per Allen et al. (2005) Cell, 121:207-221 [PMID: 15851028].

10: Paired: String representing paired positions in the query and transcript. The format is Query5'-Query3',Transcript3'-Transcript5'.  Positions are one-based. Discrete blocks of pairing are separated by ;

11: Unpaired: String representing unpaired positions in the query and transcript. The format is Query5'-Query3',Transcript3'-Transcript5'[code].  Possible codes are "UP5" (Unpaired region at 5' end of query), "UP3" (Unpaired region at 3' end of query), "SIL" (symmetric internal loop), "AILt" (asymmetric internal loop with more unpaired nts on the transcript side), "AILq" (asymmetric internal loop with more unpaired nts on the query side), "BULt" (Bulged on the transcript side), or "BULq" (bulged on the query side). Positions are one-based. Discrete blocks of pairing are separated by ;

12: Structure: Aligned secondary structure. The region before the "&" represents the transcript, 5'-3', while the region after the "&" represents the query, 5'-3'.  "(" represents a transcript base that is paired, ")" represents a query based that is paired, "." represents an unpaired base, and "-" represents a gap inserted to facilitate alignment.

13: Sequence: Aligned sequence. The region before the "&" represents the transcript, 5'-3', while the region after the "&" represents the query, 5'-3'. 

