#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use Math::CDF 'pbinom';

my $version_number = "4.4";
my $help = help_message($version_number);

# if there are no arguments, return the help message and quit
unless($ARGV[0]) {
    die "$help";
}

# declare global variables
our($opt_h,$opt_v,$opt_q,$opt_d,$opt_e,$opt_p,$opt_o,$opt_g,$opt_u,$opt_n,$opt_a,$opt_r,$opt_t,$opt_c);

# Initialize opt_p and opt_c defaults
$opt_p = 1;
$opt_c = 4;

# get options
getopts('o:d:e:g:u:n:p:r:c:hvqlat');

# If user wants help or version, give it, and die
if($opt_h) {
    die "$help";
}
if($opt_v) {
    die "CleaveLand4.pl version $version_number\n";
}

# Start speaking to user, unless quiet mode is on
unless($opt_q) {
    print STDERR "\nCleaveLand4 version $version_number\n";
    print STDERR `date`;
}

# Parse options to determine the mode
our($mode);
$mode = determine_mode();
unless($mode) {
    die "$help";
}
unless($opt_q) {
    print STDERR "\nMode: $mode\n";
    print STDERR "Checking Dependencies\n";
}


# Check dependencies
unless($mode == 4) {
    my $bowtie_check = check_bowtie_version();
    unless($bowtie_check) {
	die "FAIL\n\n$help";
    }
    my $bowtie_build_check = check_bowtie_build_version();
    unless($bowtie_build_check) {
	die "FAIL\n\n$help";
    }
}

if(($mode == 1) or ($mode == 2)) {
    my $RNAplex_check = check_RNAplex();
    unless($RNAplex_check) {
	die "FAIL\n\n$help";
    }
    my $GSTAr_check = check_GSTAr_version();
    unless($GSTAr_check) {
	die "FAIL\n\n$help";
    }
}
my $samtools_version;
if(($mode == 1) or ($mode == 3)) {
    $samtools_version = check_samtools();
    unless($samtools_version) {
	die "FAIL\n\n$help";
    }
}

if($opt_o) {
    my $R_check = check_R_version();
    unless($R_check) {
	die "FAIL\n\n$help";
    }
}
unless($opt_q) {
    print STDERR "\n";
}

# Option p must be a number between 0 and 1
unless($opt_p =~ /^[\.\d]+$/) {
    die "FATAL: Invalid entry for option p. Must be a number >0 and <= 1\n\n$help";
}
unless(($opt_p > 0) and ($opt_p <= 1)) {
    die "FATAL: Invalid entry for option p. Must be a number >0 and <= 1\n\n$help";
}

# Option c must be an integer between 0 and 4
unless($opt_c =~ /^\d$/) {
    die "FATAL: Invlaid entry for option c. Must be an integer between 0 and 4\n\n$help";
}
unless(($opt_c >= 0) and ($opt_c <= 4)) {
    die "FATAL: Invlaid entry for option c. Must be an integer between 0 and 4\n\n$help";
}    

# If option o is on, ensure the directory is created, and write the R script
my $R_script;
if($opt_o) {
    if(-d $opt_o) {
	unless($opt_q) {
	    print STDERR "T-plot directory $opt_o already exists, new T-plots will be placed there.\n";
	}
    } else {
	system "mkdir $opt_o";
	unless($opt_q) {
	    print STDERR "Created directory $opt_o for T-plots.\n";
	}
    }
    $R_script = write_R_script();
}

# if $opt_a is on for modes 3 or 4 warn the user
if((($mode == 3) or ($mode == 4)) and $opt_a) {
    # not suppressed because it is a warning
    print STDERR "\nWARNING: option a is irrelevant when running in mode $mode \!\n";
}

# Do simple verification that any file paths passed in are readable
my $readable = verify_readable();
unless($readable) {
    die "$help";
}

# Get a degradome density file path
my $dd_file;
if(($mode == 1) or ($mode == 3)) {
    $dd_file = make_deg_density($samtools_version);
} else {
    $dd_file = $opt_d;
}

# Get info from degradome file
my %deg_header = get_dd_header($dd_file);
unless(%deg_header) {
    die "$help";
}
my %deg_data = get_dd_data($dd_file);  ## key is tx id, value is all entries for the locus (inc newlines)
my %deg_tx_lens = get_dd_tx_lens($dd_file);  ## key is tx id, value is length of the tx
unless($opt_q) {
    print STDERR "Loaded degradome density data from file $dd_file\n";
    print STDERR "\tTranscriptome: $deg_header{'tx_name'}\n";
}

# get a GSTAr file path
my $GSTAr_file;
if(($mode == 1) or ($mode == 2)) {
    $GSTAr_file = make_GSTAr_file();
} else {
    $GSTAr_file = $opt_g;
}

# get info from GSTAr file
my %G_header = get_G_header($GSTAr_file);
unless(%G_header) {
    die "$help";
}
unless($opt_q) {
    print STDERR "Loaded GSTAr alignment file $GSTAr_file\n";
    print STDERR "\tTranscriptome: $G_header{'transcripts'}\n";
    print STDERR "\tQueries: $G_header{'queries'}\n";
    print STDERR "\tRanked by: $G_header{'sort'}\n";
}

# Transcriptomes must match!
unless($G_header{'transcripts'} eq $deg_header{'tx_name'}) {
    die "\nFATAL: Transcriptome has to be the same for the degradome density file and the GSTAr alignment file\!\n$help";
}


# Print a header to STDOUT
print "\# CleaveLand4 $version_number\n";
print "\# ";
print `date`;
print "\# Degradome Density File: $dd_file\n";
print "\# Query Alignment file: $GSTAr_file\n";
print "\# Transcriptome: $G_header{'transcripts'}\n";
print "\# P-value cutoff: $opt_p\n";
print "\# Category cutoff: $opt_c\n";
print "\# Output Format: ";
if($opt_t) {
    print "Tabular\n";
} else {
    print "Pretty\n";
}

# Work through all of the alignments
unless($opt_q) {
    print STDERR "\nBeginning analysis phase...\n\n";
}

(open(G, "$GSTAr_file")) || die "\nFATAL: Failed to open GSTAr file $GSTAr_file in main\n";
my $gline;
my @gfields = ();
my $cat;
my $last_q_name = "NULL";
my $rank;
my $chance;
my $p_val;
my $t_file;
my %nr_tabular = ();
my $slice_site;
my $ext_entry;
my $cat_digit;
while (<G>) {
    chomp;
    if($_ =~ /^\#/) {
	next;
    }
    if($_ =~ /^Query\tTranscript\tTStart/) {
	if($opt_t) {
	    print "SiteID\t$_\tDegradomeCategory\tDegradomePval\tTplot_file_path\n";
	}
	next; ## column names
    }
    @gfields = split ("\t", $_);
    $gline = $_;
    # track the ranks
    if($gfields[0] ne $last_q_name) {
	$rank = 0;
	unless($last_q_name eq "NULL") {
	    unless($opt_q) {
		print STDERR "Finished with Query: $last_q_name\n"
	    }
	}
    }

    $last_q_name = $gfields[0];
    ++$rank;
    
    
    $cat = check_cat(\%deg_data,\$gfields[1],\$gfields[4]);
    
    if($cat =~ /cat(\d+)$/) {
	$cat_digit = $1;
    } else {
	next;
    }
    $chance = $deg_header{$cat} / $deg_header{'tx_ef_size'};
    $p_val = 1 - (pbinom(0,$rank,$chance));
    $slice_site = "$gfields[1]" . ":" . "$gfields[4]";
    if(($p_val <= $opt_p) and ($cat_digit <= $opt_c)) {
	# check for redundancy
	if(exists($nr_tabular{$slice_site})) {
	    my @old_fields = split ("\t", $nr_tabular{$slice_site});
	    if($p_val < $old_fields[-1]) {
		# the new one is better than the old one
		$ext_entry = "$gline\t$cat_digit\t$p_val";
		$nr_tabular{$slice_site} = $ext_entry;
	    }
	} else {
	    $ext_entry = "$gline\t$cat_digit\t$p_val";
	    $nr_tabular{$slice_site} = $ext_entry;
	}
    }
}
close G;
unless($opt_q) {
    print STDERR "Finished with Query: $last_q_name\n"
}

my $n_ok = scalar(keys %nr_tabular);
unless($opt_q) {
    print STDERR "\nFound a total of $n_ok non-redundant putative slicing sites with p-values <= $opt_p and categories <= $opt_c\. Outputting results.\n";
}
my @sorted_slice_sites = sort(keys %nr_tabular);
foreach my $sss (@sorted_slice_sites) {
    @gfields = split ("\t", $nr_tabular{$sss});
    # Make T-plot if user wanted one
    if($opt_o) {
	$t_file = make_t_plot($opt_o,$deg_data{$gfields[1]}, $deg_tx_lens{$gfields[1]}, $nr_tabular{$sss}, $gfields[-2], $gfields[-1]);
    } else {
	$t_file = "NOT_CREATED";
    }
    # Report to STDOUT
    if($opt_t) {
	tabular_report_hit($nr_tabular{$sss},$t_file,$sss);
    } else {
	pretty_report_hit($deg_data{$gfields[1]}, $nr_tabular{$sss}, $gfields[-2], $gfields[-1], $dd_file, $t_file, $sss);
    }
}


# Delete the T-plot script if you wrote it and it's still there
if(($opt_o) and ($R_script)) {
    system "rm -f $R_script";
}

# Have an A-1 day!
exit;

####################
# Here be sub-routines

sub help_message {
    my($version) = @_;
    my $message = "
CleaveLand4.pl : Finding sliced targets of small RNAs from degradome data

Version: $version

Usage: CleaveLand4.pl \[options\] > \[out.txt\]

Options:
-h Print help message and quit
-v Print version and quit
-q Quiet mode .. no log/progress information to STDERR
-a Sort small RNA / transcript alignments by Allen et al. score instead of default MFEratio -- for GSTAr
-t Output in tabular format instead of the default verbose format
-r [float >0..1] Minimum Free Energy Ratio cutoff. Default: 0.65 -- for GSTAr
-o [string] : Produce T-plots in the directory indicated by the string. If the dir does not exist, it will be created
-d [string] : Path to degradome density file.
-e [string] : Path to FASTA-formatted degradome reads.
-g [string] : Path to GSTAr-created tabular formatted query-transcript alignments.
-u [string] : Path to FASTA-formatted small RNA queries 
-n [string] : Path to FASTA-formatted transcriptome
-p [float >0..1] : p-value for reporting. Default is 1 (no p-value filtering).
-c [integer 0..4] : Maximum category for reporting. Default is 4 (all categories reported).

Modes:
1. Align degradome data, align small RNA queries, and analyze. 
  REQUIRED OPTIONS: -e, -u, -n
  DISALLOWED OPTIONS: -d, -g
2. Use existing degradome density file, align small RNA queries, and analyze.
  REQUIRED OPTIONS: -d, -u, -n
  DISALLOWED OPTIONS: -e, -g
3. Align degradome data, use existing small RNA query alignments, and analyze.
  REQUIRED OPTIONS: -e, -n, -g
  DISALLOWED OPTIONS: -d, -u
  IRRELEVANT OPTIONS: -a, -r
4. Use existing degradome density file and existing small RNA query alignments, and analyze.
  REQUIRED OPTIONS: -d, -g
  DISALLOWED OPTIONS: -e, -u
  IRRELEVANT OPTIONS: -a, -r

Dependencies (must be in PATH):
  R \[only if making T-plots\]]
  GSTAr.pl \[modes 1 and 2\ .. Version 1.0 or higher]
  bowtie (0.12.x OR 1.x) \[modes 1, 2, and 3\]
  bowtie-build \[modes 1, 2, and 3\]
  RNAplex (from Vienna RNA Package) \[modes 1 and 2\]
  samtools \[modes 1 and 3\]

Documentation: perldoc CleaveLand4.pl
";
    return $message;
}

sub check_bowtie_version {
    unless($opt_q) {
	print STDERR "\tbowtie: ";
    }
    (open(BTV, "bowtie --version |")) || return 0;
    my $vline = <BTV>;
    close BTV;
    my $version;
    if($vline =~ /^bowtie version (\S+)/) {
	$version = $1;
	if(($version =~ /^0\.12/) or ($version =~ /^1\./)) {
	    unless($opt_q) {
		print STDERR "PASS version $version\n";
	    }
	    return 1;
	} else {
	    return 0;
	}
    } else {
	return 0;
    }
    # should never be here
    return 0;
}

sub check_bowtie_build_version {
    unless($opt_q) {
	print STDERR "\tbowtie-build: ";
    }
    (open(BBV, "bowtie-build --version |")) || return 0;
    my $vline = <BBV>;
    close BBV;
    if($vline =~ /^bowtie-build version (\S+)/) {
	unless($opt_q) {
	    print STDERR "PASS version $1\n";
	}
	return 1;
    } else { 
	return 0;
    }
}

sub check_RNAplex {
    unless($opt_q) {
	print STDERR "\tRNAplex: ";
    }
    (open(RD, "RNAplex --version |")) || return 0;
    my $v = <RD>;
    close RD;
    unless($opt_q) {
	print STDERR "PASS $v";
    }
    return 1;
}

sub check_samtools {
    unless($opt_q) {
	print STDERR "\tsamtools: ";
    }
    (open(SAM, "samtools 2>&1 |")) || return 0;
    my $version;
    while (<SAM>) {
	chomp;
	if($_ =~ /^Version/) {
	    $version = $_;
	}
    }
    close SAM;
    if($version) {
	unless($opt_q) {
	    print STDERR "PASS $version\n";
	}
	return $version;
    } else {
	return 0;
    }
}

sub check_GSTAr_version {
    unless($opt_q) {
	print STDERR "\tGSTAr: ";
    }
    (open(GS, "GSTAr.pl -v 2>&1 |")) || return 0;
    my $version = <GS>;
    close GS;
    if($version =~ /version (\d+)\.\d/) { ## must be higher than 0.x as of CleaveLand 4.1
	if($1 < 1) {
	    return 0;
	} else {
	    unless($opt_q) {
		print STDERR "PASS $version";
	    }
	    return 1;
	}
    } else {
	return 0;
    }
}
    
sub check_R_version {
    unless($opt_q) {
	print STDERR "\tR: ";
    }
    (open(R, "R --version |")) || return 0;
    my $version = <R>;
    close R;
    if($version =~ /version/) {
	unless($opt_q) {
	    print STDERR "PASS $version";
	}
	return 1;
    } else {
	return 0;
    }
}

sub check_ebwt {
    my($base) = @_;
    unless($opt_q) {
	print STDERR "Bowtie index files for $base ";
    }
    my @ebwt_suffixes = qw(.1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt);
    my $ebwt_count = 0;
    my $ebwt_file;
    foreach my $ebwt_suffix (@ebwt_suffixes) {
	$ebwt_file = "$base" . "$ebwt_suffix";
	if(-r $ebwt_file) {
	    ++$ebwt_count;
	}
    }
    if($ebwt_count == 6) {
	unless($opt_q) {
	    print STDERR "PRESENT\n";
	}
	return 1;
    } elsif ($ebwt_count > 0) {
	unless($opt_q) {
	    print STDERR "INCOMPLETE only $ebwt_count of the 6 expeted index files were found.\n";
	}
	die "\tABORTING .. please clean up old incompleted ebwt files for transcript $base and try again\n";
    } else {
	unless($opt_q) {
	    print STDERR "ABSENT\n";
	}
	return 0;
    }
}

sub bowtie_build {
    my($fasta) = @_;
    unless($opt_q) {
	print STDERR "Building bowtie index for file $fasta ...\n";
    }
    my $log = "$fasta" . "_bowtie_build_log.txt";
    system("bowtie-build $fasta $fasta > $log");
    unless($opt_q) {
	print STDERR " Done\n";
	print STDERR "Log for bowtie-build job is at $log\n";
    }
    return 1;
}

sub determine_mode {
    if(($opt_e) and ($opt_u) and ($opt_n)) {
	# Looks like mode 1.  Check for conflicts
	# e is incompatible with d, u is incompatible with g
	if($opt_d) {
	    print STDERR "\nFATAL: Incompatible options e and d.\n";
	    return 0;
	} elsif ($opt_g) {
	    print STDERR "\nFATAL: Incompatible options u and g.\n";
	    return 0;
	} else {
	    return 1;
	}
    } elsif (($opt_d) and ($opt_u) and ($opt_n)) {
	# Looks like mode 2. Check for conflicts
	# e is incompatible with d, u is incompatible with g
	if($opt_e) {
	    print STDERR "\nFATAL: Incompatible options e and d.\n";
	    return 0;
	} elsif ($opt_g) {
	    print STDERR "\nFATAL: Incompatible options u and g.\n";
	    return 0;
	} else {
	    return 2;
	}
    } elsif (($opt_e) and ($opt_n) and ($opt_g)) {
	# Looks like mode 3. Check for conflicts.
	# e is incompatible with d, u is incompatible with g
	if($opt_d) {
	    print STDERR "\nFATAL: Incompatible options e and d.\n";
	    return 0;
	} elsif ($opt_u) {
	    print STDERR "\nFATAL: Incompatible options u and g.\n";
	    return 0;
	} else {
	    return 3;
	}
    } elsif (($opt_d) and ($opt_g)) {
	# Looks like mode 4. Check for conflicts.
	# e is incompatible with d, u is incompatible with g
	if($opt_e) {
	    print STDERR "\nFATAL: Incompatible options e and d.\n";
	    return 0;
	} elsif ($opt_u) {
	    print STDERR "\nFATAL: Incompatible options u and g.\n";
	    return 0;
	} else {
	    return 4;
	}
    } else {
	print STDERR "\nFATAL: Option combination does not conform to any of the four modes.\n";
	return 0;
    }
}
	
sub verify_readable {
    my $fail = 0;
    if($opt_d) {
	unless(-r $opt_d) {
	    print STDERR "FATAL: File $opt_d from option d was not readable\n";
	    ++$fail;
	}
    }
    if($opt_e) {
	unless(-r $opt_e) {
	    print STDERR "FATAL: File $opt_e from option e was not readable\n";
	    ++$fail;
	}
    }
    if($opt_g) {
	unless(-r $opt_g) {
	    print STDERR "FATAL: File $opt_g from option g was not readable\n";
	    ++$fail;
	}
    }
    if($opt_u) {
	unless(-r $opt_u) {
	    print STDERR "FATAL: File $opt_u from option u was not readable\n";
	    ++$fail;
	}
    }
    if($opt_n) {
	unless(-r $opt_n) {
	    print STDERR "FATAL: File $opt_n from option n was not readable\n";
	    ++$fail;
	}
    }
    if($fail) {
	return 0;
    } else {
	return 1;
    }
}

sub make_deg_density {
    my($stv) = @_; ## samtools version
    unless($opt_q) {
	print STDERR "\nInitiating alignment of degradome reads to transcriptome.\n";
    }

    # Search for the bowtie indices
    my $ebwt_check = check_ebwt($opt_n);
    # If absent, call bowtie-build
    unless($ebwt_check) {
	bowtie_build($opt_n);
    }
    
    my $bam_name = "$opt_n" . "_sorted";
    
    # call the bowtie --> samtools pipeline to make a sorted bam file
    # Keep it quiet if needed
    # samtools sort syntax depends on samtools version.
    my $sort_call;
    if($stv =~ /^Version: 1\./) {
	my $tmp_sort = "$bam_name" . "_sorttemp";
	my $full_bam_name = "$bam_name" . '.bam';
	$sort_call = "samtools sort -T $tmp_sort -O bam - > $full_bam_name";
    } else {
	$sort_call = "samtools sort - $bam_name";
    }
    if($opt_q) {
	system "bowtie -f -v 1 --best -k 1 --norc -S $opt_n $opt_e 2> /dev/null \| sed -e 's/SO:unsorted/SO:coordinate/' 2> /dev/null \| samtools view -S -b -u - 2> /dev/null \| $sort_call 2> /dev/null";
    } else {
	system "bowtie -f -v 1 --best -k 1 --norc -S $opt_n $opt_e \| sed -e 's/SO:unsorted/SO:coordinate/' \| samtools view -S -b -u - 2> /dev/null \| $sort_call 2> /dev/null";
    }
    # Verify the bam file is there
    my $bam_file = "$bam_name" . ".bam";
    unless(-r $bam_file) {
	die "ABORT: No bam file $bam_file found after bowtie call in sub-routine make_deg_density\n";
    }
    
    # report to user 
    unless($opt_q) {
	print STDERR "Alignments complete. Generating degradome density file.\n";
    }
    
    # parse the header to get all of the tx lengths
    # also track the total number of characters
    my %tx_lens = ();
    my $tx_chars = 0;
    (open(HEADER, "samtools view -H $bam_file \|")) || die "ABORT: Failed to open bam file $bam_file in sub-routine make_deg_density\n";
    while (<HEADER>) {
	chomp;
	if($_ =~ /^\@SQ.*SN:(\S+).*LN:(\d+)/) {
	    $tx_lens{$1} = $2;
	    $tx_chars += $2;
	} elsif ($_ =~ /^\@SQ/) {
	    die "ABORT: Failure to parse bam SQ header line at $_\n";
	}
    }
    close HEADER;
    
    # create the degradome density file and write header
    my $dd_file_tmp = "$opt_e" . "_dd.txt" . ".tmp";
    (open(DD, ">$dd_file_tmp")) || die "ABORT: Failed to create initial tmp degradome density file $dd_file_tmp in sub-routine make_deg_density\n";
    print DD "\# CleaveLand4 degradome density\n";
    print DD "\# ";
    print DD `date`;
    print DD "\# Degradome Reads:$opt_e\n";
    print DD "\# Transcriptome:$opt_n\n";
    print DD "\# TranscriptomeCharacters:$tx_chars\n";
    
    # go through the alignments
    (open(SAM, "samtools view $bam_file \|")) || die "ABORT: Failed to open bam file $bam_file in sub-routine make_deg_density\n";
    my @samfields = ();
    my $last_tx = "NULL";
    my %freq = ();
    my $outstring;
    my %read_sizes = ();
    my %categorized = (
	0 => 0,
	1 => 0,
	2 => 0,
	3 => 0,
	4 => 0,
	);
    while (<SAM>) {
	chomp;
	# no headers
	if($_ =~ /^@/) {
	    next;
	}
	@samfields = split ("\t", $_);
	# no unmapped reads
	if($samfields[1] & 4) {
	    next;
	}
	# paranoid but just in case no rev comps either
	if($samfields[1] & 16) {
	    next;
	}
	
	# if it is a new tx, report the last one
	if(($samfields[2] ne $last_tx) and ($last_tx ne "NULL")) {
	    
	    # Report
	    $outstring = report_dd(\%freq,\%tx_lens,\$last_tx,\%categorized);
	    print DD "$outstring\n";
	    
	    # clear
	    %freq = ();
	}
	
	# record
	++$freq{$samfields[3]};
	++$read_sizes{(length $samfields[9])};
	
	# reset
	$last_tx = $samfields[2];
	
    }
    close SAM;
    # Report the last one
    $outstring = report_dd(\%freq,\%tx_lens,\$last_tx,\%categorized);
    print DD "$outstring\n";
    close DD;
    
    # remove the bam file
    system "rm -f $bam_file";
    
    # calculate the mean read size
    my $read_sum;
    my $read_n;
    foreach my $size (keys %read_sizes) {
	$read_n += $read_sizes{$size};
	$read_sum += ($size * $read_sizes{$size});
    }
    my $mean_read = int($read_sum / $read_n);
    
    # estimate number of tx positions that could have possibly had a read
    my $est_tx_avail = $tx_chars - ((scalar(keys %tx_lens)) * $mean_read);
    
    # re-open the DD file to add the final information to the header
    (open(DDTMP, "$dd_file_tmp")) || die "ABORT: Failed to re-open temp dd file in sub-routine make_deg_density\n";
    my $dd_file = $dd_file_tmp;
    $dd_file =~ s/\.tmp//g;
    (open(DD, ">$dd_file")) || die "ABORT: Failed to open final dd file in sub-routine make_deg_density\n";
    my $done = 0;
    while (<DDTMP>) {
	if($done) {
	    print DD "$_";
	} elsif ($_ =~ /^\@/) {
	    print DD "\# Mean Degradome Read Size:$mean_read\n";
	    print DD "\# Estimated effective Transcriptome Size:$est_tx_avail\n";
	    print DD "\# Category 0:$categorized{'0'}\n";
	    print DD "\# Category 1:$categorized{'1'}\n";
	    print DD "\# Category 2:$categorized{'2'}\n";
	    print DD "\# Category 3:$categorized{'3'}\n";
	    print DD "\# Category 4:$categorized{'4'}\n";
	    $done = 1;
	    print DD "$_";
	} else {
	    print DD "$_";
	}
    }
    close DDTMP;
    close DD;
    
    # remove DD tmp file
    system "rm -f $dd_file_tmp";
    
    # report to user
    unless($opt_q) {
	print STDERR "Degradome density file completed ... $dd_file\n";
    }
    return $dd_file;
}

sub report_dd {
    my($freq,$lens,$id,$cat) = @_;  ## passed by reference .. hash, hash, scalar
    my $out;
    # print the info for the transcript
    $out .= "\@ID:$$id\n";
    if(exists($$lens{$$id})) {
	$out .= "\@LN:$$lens{$$id}\n";
    } else {
	die "ABORT: Failed to find length of transcript $$id in hash in sub-routine report_dd\n";
    }
    
    # get the max and mean 
    my $max = 0;
    my $mean;
    my $pos;
    my $sum;
    my $occupied = 0;
    foreach $pos (keys %$freq) {
	++$occupied;
	$sum += $$freq{$pos};
	if($$freq{$pos} > $max) {
	    $max = $$freq{$pos};
	}
    }
    $mean = $sum / $occupied;
    
    # how many positions are maximal?
    my $n_max = 0;
    foreach $pos (keys %$freq) {
	if($$freq{$pos} == $max) {
	    ++$n_max;
	}
    }
    
    # now write it
    my @sorted = sort {$a <=> $b} (keys %$freq);
    foreach $pos (@sorted) {
	$out .= "$pos\t";
	$out .= "$$freq{$pos}\t";
	if($$freq{$pos} == 1) {
	    $out .= "4\n";
	    ++$$cat{'4'};
	} elsif ($$freq{$pos} <= $mean) {
	    $out .= "3\n";
	    ++$$cat{'3'};
	} elsif ($$freq{$pos} == $max) {
	    if($n_max > 1) {
		$out .= "1\n";
		++$$cat{'1'};
	    } else {
		$out .= "0\n";
		++$$cat{'0'};
	    }
	} else {
	    $out .= "2\n";
	    ++$$cat{'2'};
	}
    }
    return $out;
}

sub get_dd_header {
    my($file) = @_;
    (open(DD, "$file")) || die "ABORT: Failed to open degradome density file $file in sub-routine get_dd_header\n";
    my $one = <DD>;
    chomp $one;
    unless($one eq "\# CleaveLand4 degradome density") {
	print STDERR "FATAL: File $file does not have the expected header for a CleaveLand4 degradome density file\n";
	return 0;
    }
    my %out = ();
    while (<DD>) {
	chomp;
	unless($_ =~ /^\#/) {
	    last;
	}
	if($_ =~ /Transcriptome:(\S+)/) {
	    $out{'tx_name'} = $1;
	} elsif ($_ =~ /Estimated effective Transcriptome Size:(\d+)/) {
	    $out{'tx_ef_size'} = $1;
	} elsif ($_ =~ /Category 0:(\d+)/) {
	    $out{'cat0'} = $1;
	} elsif ($_ =~ /Category 1:(\d+)/) {
	    $out{'cat1'} = $1;
	} elsif ($_ =~ /Category 2:(\d+)/) {
	    $out{'cat2'} = $1;
	} elsif ($_ =~ /Category 3:(\d+)/) {
	    $out{'cat3'} = $1;
	} elsif ($_ =~ /Category 4:(\d+)/) {
	    $out{'cat4'} = $1;
	}
    }
    close DD;
    if((exists($out{'tx_name'})) and
	   (exists($out{'tx_ef_size'})) and
	   (exists($out{'cat0'})) and
	   (exists($out{'cat1'})) and
	   (exists($out{'cat2'})) and
	   (exists($out{'cat3'})) and
	   (exists($out{'cat4'}))) {
	return %out;
    } else {
	print STDERR "FATAL: File $file does not have the expected header for a CleaveLand4 degradome density file\n";
	return 0;
    }
}

sub get_dd_data {
    my($file) = @_;
    (open(DD, "$file")) || die "ABORT: Failed to open dd file $file in sub-routine get_dd_data\n";
    my %out = ();
    my $string;
    my $last_id;
    while (<DD>) {
	if($_ =~ /^\@ID:(\S+)/) {
	    if($last_id) {
		$out{$last_id} = $string;
	    }
	    $string = '';
	    $last_id = $1;
	} elsif ($_ =~ /^\d/) {
	    $string .= $_;
	}
    }
    close DD;
    ## last one
    $out{$last_id} = $string;
    return %out;
}
    
sub get_dd_tx_lens {
    my($file) = @_;
    (open(DD, "$file")) || die "ABORT: Failed to open dd file $file in sub-routine get_dd_data\n";
    my %out = ();
    my $id;
    my $length;
    while (<DD>) {
	if($_ =~ /^\@ID:(\S+)/) {
	    $id = $1;
	} elsif ($_ =~ /^\@LN:(\d+)/) {
	    $length = $1;
	    $out{$id} = $length;
	}
    }
    close DD;
    return %out;
}

sub make_GSTAr_file {
    unless($opt_q) {
	print STDERR "Initiating GSTAr run.\n";
    }
    my $command = "GSTAr.pl -t";
    if($opt_q) {
	$command .= " -q";
    }
    if($opt_a) {
	$command .= " -a";
    }
    if($opt_r) {
	$command .= " -r $opt_r";
    }
    $command .= " $opt_u $opt_n";
    
    my $file = "$opt_u" . "_GSTAr.txt";
    $command .= " > $file";
    
    unless($opt_q) {
	print STDERR "GSTAr command is $command\n\n";
    }
    system "$command";
    
    unless($opt_q) {
	print STDERR "\n\nGSTAr run complete\n";
    }

    return $file;
}

sub get_G_header {
    my($file) = @_;
    my %out = ();
    (open(G, "$file")) || die "ABORT: Failed to open GSTAr file $file in sub-routine get_G_header\n";
    my $one = <G>;
    chomp $one;
    unless($one =~ /^\# GSTAr/) {
	print STDERR "FATAL: File $file does not appear to be a GSTAr alignment file\n";
	return 0;
    }
    while (<G>) {
	chomp;
	unless($_ =~ /^\#/) {
	    last;
	}
	if($_ =~ /Queries: (\S+)/) {
	    $out{'queries'} = $1;
	} elsif ($_ =~ /Transcripts: (\S+)/) {
	    $out{'transcripts'} = $1;
	} elsif ($_ =~ /Sorted by: (\S+)/) {
	    $out{'sort'} = $1;
	} elsif ($_ =~ /Output Format: (\S+)/) {
	    $out{'format'} = $1;
	}
    }
    close G;
    if(($out{'queries'}) and
       ($out{'transcripts'}) and
       ($out{'sort'}) and
       ($out{'format'})) {
	
	if($out{'format'} eq "Tabular") {
	    return %out;
	} else {
	    print STDERR "FATAL: GSTAr file $file is not in tabular format\n";
	    return 0;
	}
    } else {
	print STDERR "FATAL: File $file does not appear to be a GSTAr alignment file\n";
	return 0;
    }
}

	
sub check_cat {
    my($deg_data,$id,$slice) = @_;  ## passed by reference .. hash, scalar, scalar
    my $cat = "cat";
    # Is there any data at all for the transcript?
    unless(exists($$deg_data{$$id})) {
	return $cat;
    }
    
    # get the lines
    my @lines = split ("\n", $$deg_data{$$id});
    my @fields = ();
    foreach my $line (@lines) {
	@fields = split ("\t", $line);
	if($fields[0] == $$slice) {
	    $cat .= "$fields[2]";
	    last;
	}
    }
    return $cat;
}

sub pretty_report_hit {
    my($deg_data,$gline,$cat,$p_val,$dd_file,$tfile,$key) = @_;
    print "---------------------------------------------------\n";
    # First, print the gline in pretty format
    my @fields = split ("\t", $gline);
    print "\n5\' ";
    my @tseq = get_al_array($fields[12],0);
    print (@tseq," 3\'");
    print " Transcript: $fields[1]",":","$fields[2]","-","$fields[3]"," Slice Site:$fields[4]\n";
    my @tbrax = get_al_array($fields[11],0);
    my @qbrax = reverse(get_al_array($fields[11],1));
    my @qseq = reverse(get_al_array($fields[12],1));
    my $i = -1;
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
    print " Query: $fields[0]\n\n";
    print "SiteID: $key\n";
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
		die "ABORT: in sub-routine report_hit : failed to parse unpaired string $unp\n";
	    }
	}
    }
    
    # Now, print the hit information
    print "\nDegradome data file: $dd_file\n";
    print "Degardome Category: $cat\n";
    print "Degradome p-value: $p_val\n";
    print "T-Plot file: $tfile\n\n";
    # And the details
    print "Position\tReads\tCategory\n";
    my @dlines = split ("\n", $deg_data);
    my @dfields = ();
    foreach my $dline (@dlines) {
	print "$dline";
	@dfields = split ("\t", $dline);
	if($dfields[0] == $fields[4]) {
	    print "\t<<<<<<<<<<\n";
	} else {
	    print "\n";
	}
    }
    
    print "---------------------------------------------------\n";
    return 1;
}
    
sub get_al_array {
    my($string,$index) = @_;
    my @two = split ("\&", $string);
    my $chunk = $two[$index];
    my @out = split ('', $chunk);
    return @out;
}
    
    
sub tabular_report_hit {
    my($base,$tfile,$key) = @_;
    print "$key\t$base\t$tfile\n";
}

sub make_t_plot {
    my($dir_name,$deg_data,$tx_len,$gline,$cat_digit,$pval) = @_;

    # open tmp file for data
    my $args1 = "Tplot_tmp.txt";
    (open(TMP, ">$args1")) || die "ABORT: Failed to open a temp file for Tplot data in sub-routine make_t_plot\n";
    
    # get required arguments
    my @gfields = split ("\t", $gline);
    my $args2 = "$dir_name" . "\/" . "$gfields[0]" . "_" . "$gfields[1]" . "_" . "$gfields[4]" . "_TPlot.pdf";
    my $args3 = "T\=$gfields[1]" . "_" . "Q\=$gfields[0]" . "_" . "S\=$gfields[4]";
    my $args4 = "category\=$cat_digit" . "_" . "p\=$pval";
    
    # print out the data to temp file
    print TMP "Position\tAll\tSite\n";
    my $i;
    my $junk;
    my @deg_lines = split ("\n", $deg_data);
    for($i = 1; $i <= $tx_len; ++$i) {
	print TMP "$i\t";
	if(@deg_lines) {
	    if($deg_lines[0] =~ /^$i\t(\d+)/) {
		print TMP "$1\t";
		if($i == $gfields[4]) {
		    print TMP "$1\n";
		} else {
		    print TMP "NA\n";
		}
		$junk = shift @deg_lines;
	    } else {
		print TMP "0\tNA\n";
	    }
	} else {
	    print TMP "0\tNA\n";
	}
    }
    close TMP;
    
    # Call R
    system "R --vanilla --quiet --slave --args $args1 $args2 $args3 $args4 < CleaveLand4_Tplotter.R 2> /dev/null";
    
    # clean up TMP
    system "rm -f $args1";
    
    # return the file path
    return $args2;
}

sub write_R_script {
    my $file = "CleaveLand4_Tplotter.R";
    (open(R, ">$file")) || die "ABORT: Failed to open R script file in sub-routine write_R_script\n";
    print R "args <- commandArgs(trailingOnly = TRUE)\n";
    print R "data <- read.table(args[1], header=TRUE, sep=\"\\t\")\n";
    #### args[1] = file of data
    print R "pdf(file=args[2])\n";
    #### args[2] = name of output pdf file to make
    print R "plot(data\$Position, data\$All, type=\"l\", xlab=\"Transcript Position\", ylab=\"Degradome 5\' end Frequency\", main=args[3])\n";
    #### args[3] = Title to print on the plot
    print R "points(data\$Position, data\$Site, pch=16, cex=2, col=\"red\")\n";
    print R "mtext(args[4])\n";
    #### args[4] = mtext
    print R "crap <- dev.off()\n";
    close R;
    return $file;
}

__END__

=head1 LICENSE

CleaveLand4.pl

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

CleaveLand4 : Finding evidence of sliced targets of small RNAs from degradome data

=head1 AUTHOR

Michael J. Axtell, Penn State University, mja18@psu.edu

=head1 VERSION

4.3 : November 7, 2013

=head1 INSTALL

=head2 Dependencies - Required Perl Modules

	Getopt::Std
	Math::CDF

CleaveLand4.pl is a perl program, so it needs perl installed on your system. It was developed on perl version 5.10.0, and hasn't been tested on other versions (but there is no reason to suspect problems with other perl 5.x versions). CleaveLand4.pl will not compile unless the above two Perl modules are also installed in your perl's @INC. Getopt::Std is pre-loaded into most (all?) Perl distros. But you may need to install Math::CDF from CPAN. Only one Math::CDF function is used by CleaveLand4 ('pbinom').

=head2 Dependencies - PATH executables

	bowtie (version 0.12.x or 1.x)
	bowtie-build
	RNAplex (from Vienna RNA package)
	GSTAr.pl (Version 1.0 or higher -- distrubuted with CleaveLand4)
	R
	samtools
	
	
All of the above must be executable from your PATH. Depending on the mode of the CleaveLand4 run (see below), only a subset of these programs may be required for a given run.

=head2 Installation

Except for the above dependencies, there is no "real" installation. If the script is in your working directory, you can call it with 

	./CleaveLand4.pl

For convenience, you can add it to your PATH.  e.g.

	sudo mv CleaveLand4.pl /usr/bin/

GSTAr.pl expects to find perl in /usr/bin/perl .. if not, edit line 1 (the hashbang) accordingly.

=head1 USAGE

	CleaveLand4.pl [options] > [out.txt]

Log and progress information goes to STDERR, and can be suppressed with option -q (quiet mode).

=head2 Options  
                                                                                                                                                                          
-h Print help message and quit

-v Print version and quit

-q Quiet mode .. no log/progress information to STDERR

-a Sort small RNA / transcript alignments by Allen et al. score* instead of default MFEratio -- for GSTAr

-t Output in tabular format instead of the default verbose format

-r [float >0..1] Minimum Free Energy Ratio cutoff. Default: 0.65 -- for GSTAr

-o [string] : Produce T-plots in the directory indicated by the string. If the dir does not exist, it will be created

-d [string] : Path to degradome density file.

-e [string] : Path to FASTA-formatted degradome reads.

-g [string] : Path to GSTAr-created tabular formatted query-transcript alignments.

-u [string] : Path to FASTA-formatted small RNA queries 

-n [string] : Path to FASTA-formatted transcriptome

-p [float >0..1] : p-value for reporting. Default is 1 (no p-value filtering).

-c [integer 0..4] : Maximum category for reporting. Default is 4 (all categories reported).

*Allen et al. score: This is a score based on the position-specific penalties described by Allen et al. (2005) Cell, 121:207-221 [PMID: 15851028].  Specifically, mismatched query bases or target-bulged bases, are penalized 1. G-U wobbles are penalized 0.5. These penalties are double within positions 2-13 of the query.

=head2 Modes

CleaveLand4 runs in one of four different modes. Each mode has a required set of options and a disallowed set of options, as described below:

Mode 1: Align degradome data and create degradome density file, perform new small RNA query/transcriptome alignment with GSTAr, and analyze. Required options: -e, -u, -n. Disallowed options: -d, -g.

Mode 2: Use existing degradome density file, perform new small RNA query/transcriptome alignment with GSTAr, and analyze. Required options: -d, -u, -n. Disallowed options: -e, -g.

Mode 3: Align degradome data and create degradome density file, use existing GSTAr alignments, and analyze. Required options: -e, -n, -g. Disallowed options: -d, -u.

Mode 4: Use existing degradome density file and existing GSTAr alignments, and analyze. Required options: -d, -g. Disallowed options: -e, -u.

=head1 METHODS

=head2 Degradome data --> transcriptome alignments --> degradome density file creation (modes 1 and 3)

Degradome data is aligned to the reference transcriptome using bowtie. If needed, the bowtie indices for the transcript are built with bowtie-build using default parameters. This results in the creation of six files, each including "ebwt" in their suffix.  Alignment parameters allow zero or one mismatch, and are only allowed to the forward strand of the transcriptome. In the case of multiple valid alignments only one is randomly selected and reported. (The specific bowtie command used is "bowtie -f -v 1 --best -k 1 --norc -S"). The alignment process uses samtools to generate a sorted BAM alignment file from the bowtie output stream.

After creation of the sorted BAM alignment file, the alignments are parsed to quantify the density of observed 5' ends at each nt of the transcriptome. The results are written to a 'degradome density' file in the working directory. The BAM alignment file is deleted at the completion of this process. The degradome density files contain the position of the transcript, the number of 5' ends at that position, and the degradome peak 'category'. Categories are determined as follows:

Category 4: Just one read at that position

Category 3: >1 read, but below or equal to the average* depth of coverage on the transcript

Category 2: >1 read, above the average* depth, but not the maximum on the transcript

Cateogry 1: >1 read, equal to the maximum on the transcript, when there is >1 position at maximum value

Cateogry 0: >1 read, equal to the maximum on the transcript, when there is just 1 position at the the maximum value

* Note that the average does not include all of the 'zeroes' for non-occupied positions within a transcript. Instead, it is the average of all positions that have at least one read.

=head2 Small RNA --> transcriptome alignments with GSTAr (modes 1 and 2)

Potential target sites are generated with GSTAr.pl, which ships with CleaveLand4.  Options -r and -a are passed to GSTAr. By default, potential target sites are sorted in descending order by MFE ratio. If the -a switch is present, this is changed to ascending order based on Allen et al. score. Note that GSTAr can be expected to take 90-120 seconds per query when analyzing a typically sized transcritpome. GSTAr.pl is always called to output in tabular format by CleaveLand4.pl. Only alignments to the reverse-comp strand of the transcriptome are considered. Upon completion, a GSTAr alignment file is written to the working directory. See the GSTAr documentation for more details on this program.

=head2 Analysis (all modes)

After loading valid degradome density and GSTAr alignment files, CleaveLand4 first checks to ensure that the transcriptome (as noted in the headers) is the same. If so, analysis progresses.  For each alignment in the GSTAr alignment file, CleaveLand4 searches the degradome density file to see if there are any degradome reads at the predicted slicing site. If there are, a p-value is calculated. The p-value takes into account both the noise in the degradome density file and the quality of the small RNA-transcriptome alignment. First, the chances of observing a deagrdome 'peak' of the given category by random chance is calculated. The chance is the total number of peaks of the given category divided by the effective transcriptome size*.  Then, the quality of the alignment is simply the rank of the alignment in the GSTAr alignment file (which is either sorted by MFE ratio [default] or by Allen et al. score). The p-value is calculated as the binomial probability of observing one or more 'hits' in 'x' trials given probability 'c', where 'x' is the rank of the alignment, and 'c' is the chance described above.

* The effective transcriptome size is the total number of bases in the transcriptome - (n * mean_read_size), where 'n' is the number of transcripts.  This adjustment accounts for the fact that the very ends of the transcripts could not possibly have any mapped 5' ends.

Any hits with a p-value <= the cutoff specified by option -p AND a category <= the cutoff specified by option -c are output to STDOUT. By default, all hits are reported (option p default is 1, option c default is 4).

=head1 INPUT FILE FORMAT REQUIREMENTS

=head2 Newlines

All files are assumed to have "\n" as newline characters. Files with MS-DOS text encoding, or others, that do not conform to this assumption will cause unexpected behavior and likely meaningless results.

=head2 Transcriptome (option -n)

This must be a multiline FASTA file. The headers should be short and simple and devoid of whitespace (e.g. ">AT1G12345" is good, ">AT1G12345 | this is my favorite gene | it is awesome" is not. The filename of the transcriptome file should also be devoid of whitespace.

=head2 Degradome reads (option -e)

This must be a multiline FASTA file. The reads are assumed to have already clipped to remove adapters. Furthermore, the reads must not have been collapsed in any way. In other words, each read off the sequencer should have an entry. Sequences that appear 50 times in the raw data from 50 different reads should each have a line.

Finally, CleaveLand4 assumes that each degradome read represents the 5'-3' sequence of a transcript, and that the first nt of each read represents the 5' end of an RNA.

=head2 Small RNA Queries (option -u)

This must be a multiline FASTA file with the full sequence of a given small RNA on one line (e.g. each line is either a header beginning with ">" or the full-length sequence of the small RNA). The headers should be short and simple and devoid of whitespace (e.g. ">ath-miR169a" is good, ">ath-miR169a MIMAT0000200 Arabidopsis thaliana miR169a" is not. Sequences can have either T's or U's.

Note: miRBase often has several mature miRNAs with exactly the same sequence, relfecting paralogous miRNA genes within a species. There is no use querying the same sequence multiple times, so it is a good idea to collapse the redundancy by query sequence when making a file of small RNA queries.

=head2 Degradome density files (option -d)

Most of the time, these will be files created by previous runs of CleaveLand4 that will have the suffix "_dd.txt". If you don't like the alignment parameters that CleaveLand4 uses, you could create your own degradome denstiy files. The format specification is:

Header region: Lines begin with "#".  Here is an example:

[line1]# CleaveLand4 degradome density

[line2]# Fri Sep 13 09:21:38 EDT 2013

[line3]# Degradome Reads:GSM278335.fasta

[line4]# Transcriptome:TAIR10_cdna_20110103_rgmupdated_cleand.fasta

[line5]# TranscriptomeCharacters:51074197

[line6]# Mean Degradome Read Size:20

[line7]# Estimated effective Transcriptome Size:50402157

[line8]# Category 0:16430

[line9]# Category 1:3456

[line10]# Category 2:95747

[line11]# Category 3:207062

[line12]# Category 4:78279

CleaveLand4 demands that the first line of a valid degradome density file is "# CleaveLand4 degradome density". It also requires all other header lines, except the date on line 2, to be present. All of this information is required for analysis.

Data Region: Each transcript begins with two lines as follows:

[line1]@ID:AT1G50920.1

[line2]@LN:2394

The @ID: gives the name of the transcript, while the @LN: gives the length of the transcript.
After this, each line gives a one-based position, the number of 5' ends at that positions, and the degradome category. The data lines are tab-delimited. Note that positons with zero reads are NOT shown.

=head2 GSTAr query-transcriptome alignments (option -g)

These are files created by GSTAr. If they were created as part of a CleaveLand4 run, they will have the suffix "_GSTAr.txt". They must be in the 'tabular' format, and have a proper header as shown below:

[line1]# GSTAr version 1.0

[line2]# Thu Sep 12 13:56:58 EDT 2013

[line3]# Queries: test_mir.fasta

[line4]# Transcripts: TAIR10_cdna_20110103_rgmupdated_cleand.fasta

[line5]# Hit seed length required to initiate RNAplex analysis (option -s): 7

[line6]# Minimum Free Energy Ratio cutoff (option -r): 0.65

[line7]# Sorted by: MFEratio

[line8]# Output Format: Tabular

Is is strongly recommended NOT to try to produce these files by means not involving CleaveLand4/GSTAr.

=head1 OUTPUT

=head2 Pretty format

By default, CleaveLand prints hits that pass the p-value and category filters to STDOUT in a human-readble, verbose format that is self-explanatory. A header (lines beginning with "#") is printed giving basic information on the analysis.

=head2 Tabular format

If option -t is specified, any hits passing the p-value filter are printed in a tab-delimited format. First, a header (lines beginning with "#") is printed giving basic information on the analysis. After that, a line giving the names of the columns is printed. Each subsequent line gives information on a single hit. The format is similar to that of the GSTAr alignments. Column information is:

1: SiteID: A unique name (within the scope of the output of a particular run) for the putative slicing site. In the form "[transcript]:[slice_site]". The output is sorted by SiteID.

2: Query: Name of query

3: Transcript: Name of transcript

4: TStart: One-based start position of the alignment within the transcript

5: TStop: One-based stop position of the alignment within the transcript

6: TSLice: One-based position of the alignment opposite position 10 of the query

7: MFEperfect: Minimum free energy of a perfectly matched site (approximate)

8: MFEsite: Minimum free energy of the alignment in question

9: MFEratio: MFEsite / MFEperfect

10: AllenScore: Penalty score calculated per Allen et al. (2005) Cell, 121:207-221 [PMID: 15851028].

11: Paired: String representing paired positions in the query and transcript. The format is Query5'-Query3',Transcript3'-Transcript5'.  Positions are one-based. Discrete blocks of pairing are separated by ;

12: Unpaired: String representing unpaired positions in the query and transcript. The format is Query5'-Query3',Transcript3'-Transcript5'[code].  Possible codes are "UP5" (Unpaired region at 5' end of query), "UP3" (Unpaired region at 3' end of query), "SIL" (symmetric internal loop), "AILt" (asymmetric internal loop with more unpaired nts on the transcript side), "AILq" (asymmetric internal loop with more unpaired nts on the query side), "BULt" (Bulged on the transcript side), or "BULq" (bulged on the query side). Positions are one-based. Discrete blocks of pairing are separated by ;

13: Structure: Aligned secondary structure. The region before the "&" represents the transcript, 5'-3', while the region after the "&" represents the query, 5'-3'.  "(" represents a transcript base that is paired, ")" represents a query based that is paired, "." represents an unpaired base, and "-" represents a gap inserted to facilitate alignment.

14: Sequence: Aligned sequence. The region before the "&" represents the transcript, 5'-3', while the region after the "&" represents the query, 5'-3'. 

15: DegradomeCategory:

	Category 4: Just one read at that position

	Category 3: >1 read, but below or equal to the average* depth of coverage on the transcript

	Category 2: >1 read, above the average* depth, but not the maximum on the transcript

	Cateogry 1: >1 read, equal to the maximum on the transcript, when there is >1 position at maximum value

	Cateogry 0: >1 read, equal to the maximum on the transcript, when there is just 1 position at the the maximum value

16: DegradomePval: p-value for this degradome hit.

17: Tplot_file_path: File path for the T-plot of this hit.

* Note that the average does not include all of the 'zeroes' for non-occupied positions within a transcript. Instead, it is the average of all positions that have at least one read.
=head2 T-plots

If the user requests them by including the -o option, T-plots for each hit that passes the p-value cutoff are created and written to the directory specificied by option -o. The black line on the plot shows all of the degradome data, and the red dot shows the putative slicing site. The title of each T-plot indicates the transcript ("T="), the query ("Q="), and the putative slicing site ("S="), as well as the category and p-value.

Existing T-plot files in the -o directory with the same name will be over-written without warning.

=head1 WARNINGS

=head2 Don't believe the hype - part 1

Under default settings, CleaveLand4 reports ALL putative slicing sites with ANY degradome reads at all, regardless of the liklihood of a given hit of being due to random chance. Without any filtering, most of your hits probably ARE due to random chance. This means that there will be many many hits of Category 4 (just one read) and/or at very high p-values. Exercise skepticism when interpreting these results. More confidence in the reality of a given slicing event can come from restricting analysis to hits with low p-values and/or high categories, and, even better, observing the slicing event in multiple degradome libraries.

=head2 Don't believe the hype - part 2

The p-value calculation is built around the ASSUMPTION that the rank order of alignments for a given query reflects their liklihood of being functional. Under default settings, GSTAr will sort the alignments for each query based on descending MFE ratio. Alternatively, when option a is specified for the GSTAr run, they will be sorted in ascending order according the Allen et al. score.  The extent to which the p-values are trustworthy is dictated directly by the extent to which you believe that MFEratio or Allen et al. scores are predictive of function. If you don't believe that, you should treat the p-values with due skepticism, and focus instead on high category hits and reproducibility between libraries.

=head2 Not for whole genomes

Degradome alignment by CleaveLand4 only searches the top strand of the transcriptome. Also, GSTAr holds the entire contents of the transcripts.fasta file in memory to speed the isolation of sub-sequences, and CleaveLand4 holds the entire contents of the degradome density file in memory. This will be impractical in terms of memory usage if a user attempts to load a whole genome.  Similarly, GSTAr will only search for pairing between the top strand of the transcripts.fasta file, making it also impractical for a genome analysis, where sites might be on either strand.

=head2 Temp files

CleaveLand4 writes several temp files during the course of a run.  So, don't mess with them during a run. In addition, it is a very bad idea to have two CleaveLand4 runs operating concurrently from the same working directory. CleaveLand4 will clean up all temp files at the conclusion of a run.

=head2 Not too fast in modes 1 and 2 (and maybe 3).

GSTAr is a very fast intermolecular RNA-RNA hybridization calculator. But when applied to whole transcriptomes, it is still very time-consuming. When running in modes 1 or 2, plan on about 90-120 seconds per query during the GSTAr phase.  Additionally, bowtie alignments and index bulding (modes 1 or 3) can also be time-consuming. Finally, requesting T-plots can also slow things down, especially if a great number of hits are being returned.

=head2 No ambiguity codes

Query sequences with characters other than A, T, U, C, or G (case-insensitive) will not be analyzed, and a warning will be sent to the user. Transcript sub-sequences for potential query alignments will be *silently* ignored if they contain any characters other than A, T, U, C,or G (case-insensitive).

=head2 Small queries

GSTAR demands that query sequences must be small (15-26 nts). Queries that don't meet these size requirements will not be analyzed and a warning sent to the user.

=head2 No redundancy

For a GIVEN QUERY, GSTAr alignments are non-redundant in terms of the slicing site of the alignment. However, a single query can have multiple overlapping alignment patterns that have differing predicted slicing sites. In addition, if multiple queries are similar in sequence, the same alignment (in terms of putative slicing site) can be present multiple times for different queries. If CleaveLand4 identifies more than one alignment at the same putative slicing site, it will only report the one with the best (lowest) p-value (subject to the maximum allowed p-value, option -p). Therefore there should be no redundancy in the putative slice sites returned by a given CleaveLand4 run.

=head2 No reverse-compatibility

Degradome density files created by versions of CleaveLand prior to 4.0 are NOT compatible with CleaveLand 4. Sorry.

=head2 Change in category definitions

The categories used by CleaveLand4 differ slightly from those used in CleaveLand3 and earlier. Specifically, categories 3 and 2 now rely upon calculating the mean, not the median, level of coverage in the transcript. In addition, transcript positions with zero coverage are no longer included in the calculation. The effect of this is to make category 2 hits much more rare, and category 3 hits much more common.

=head2 Slicing at position 10

CleaveLand4 only looks for evidence of slicing at position 10 relative to the aligned small RNA. There is no ambiguity -- data at position 11 or 9 is not relevant to CleaveLand4. This is because, as far as I know, there is no direct evidence showing Argonaute proteins cut anywhere besides position 10. However, there IS clear evidence for isomirs: lower-abundance variants of miRNAs with alternative 5' or 3' ends. Isomirs with alternative 5' ends could certainly cause offset slicing. If you wish to search for slicing at slightly 'off' locations with CleaveLand4, you will need to explictly query with the isomirs of interest.

