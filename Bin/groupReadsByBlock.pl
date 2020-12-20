#! /usr/bin/perl
use strict;
use FindBin qw($Bin);
no strict "refs";

my $root_path = $Bin;
my $outDir = $ARGV[3];
mkdir $outDir unless -d $outDir;

my $asmDir = $ARGV[4];
my $shellDir = $ARGV[5];
my $refGenome = $ARGV[6];

my %hash;
my %num;
my %FH;
my %BL;
my %CHR;
open VCF, $ARGV[0] or die $!;
while(<VCF>){
	chomp;
	next if /^#/;
	next if /HOMVAR/; #HOMVAR
	next if /HETVAR/; #unphased HETVAR
	my $info = (split /\s+/, $_)[-1];
	my ($block, $BA1, $BA2);
#	$_ =~ /(.+):(.*):(.*):(.*):(.*):(.*)/;
#	($block, $BA1, $BA2) = ($2, $5, $6);
	$info =~ /(.+):(.+):(.*):(.*):(.*):(.*):(.*)/;
	($block, $BA1, $BA2) = ($3, $6, $7);
	
	my @group1 = split /;/, $BA1;
	my @group2 = split /;/, $BA2;
	for (my $i=0; $i<@group1; $i++){
		my $barcode = $group1[$i];
		$num{$barcode}{$block}{0}++;
	}
	for (my $i=0; $i<@group2; $i++){
        my $barcode = $group2[$i];
		$num{$barcode}{$block}{1}++;
    }
	$block =~ /(.+)_(\d+)_(\d+)/;
	my $chr = $1;
	$BL{$block} =$chr;
	$CHR{$chr} =1;
}

foreach my $barcode (keys %num){
	foreach my $block (keys %{$num{$barcode}}){
		if ($num{$barcode}{$block}{0} > $num{$barcode}{$block}{1}){
			$hash{$barcode}{$block} =0;
		}elsif ($num{$barcode}{$block}{0} < $num{$barcode}{$block}{1}){
			$hash{$barcode}{$block} =1;
		}else{
			next;
		}
	}
}
close VCF;
foreach my $block (keys %BL){
	my $hap0fq1 = "$block.hap0.read.1";
	my $hap0fq2 = "$block.hap0.read.2";
	my $hap1fq1 = "$block.hap1.read.1";
	my $hap1fq2 = "$block.hap1.read.2";
	my $fh1fq1 = \*{$hap0fq1};
	my $fh1fq2 = \*{$hap0fq2};
	my $fh2fq1 = \*{$hap1fq1};
	my $fh2fq2 = \*{$hap1fq2};
	open $fh1fq1, ">$outDir/$block.hap0.read.1.fq" or die $!;
	open $fh1fq2, ">$outDir/$block.hap0.read.2.fq" or die $!;
	open $fh2fq1, ">$outDir/$block.hap1.read.1.fq" or die $!;
	open $fh2fq2, ">$outDir/$block.hap1.read.2.fq" or die $!;
	$FH{$hap0fq1}=1;
	$FH{$hap0fq2}=1;
	$FH{$hap1fq1}=1;
	$FH{$hap1fq2}=1;
}

open FQ1, $ARGV[1] or die $1;
open FQ2, $ARGV[2] or die $1;
my $num;
while(<FQ1>){
#	$num++;
#	if ($num%4==1){
		my $read1_id = $_;
		my $read1_base = <FQ1>;
		my $read1_symb = <FQ1>;
		my $read1_qual = <FQ1>;
		my $read2_id = <FQ2>;
		my $read2_base = <FQ2>;
        my $read2_symb = <FQ2>;
        my $read2_qual = <FQ2>;
#		$read1_id =~ /(.+)\#(.+)\/1/;
		$read1_id =~ /(.+)\/(.+)\/(.+)\/1/;		
		my $barcode = $2;
		foreach my $block (keys %{$hash{$barcode}}){
			my $hap = $hash{$barcode}{$block};
			my $fq1 = "$block.hap$hap.read.1";
			my $fq2 = "$block.hap$hap.read.2";
			my $fhfq1 = \*{$fq1};	
			my $fhfq2 = \*{$fq2};
			print $fhfq1 "$read1_id"."$read1_base"."$read1_symb"."$read1_qual";
			print $fhfq2 "$read2_id"."$read2_base"."$read2_symb"."$read2_qual";
		}
#	}
}
close FQ1;
close FQ2;
foreach my $hapfq (keys %FH){
	my $fh = \*{$hapfq};
	close $fh;
}
foreach my $chr (keys %CHR){
	open ASM, "> $shellDir/soapdenovo.$chr.sh" or die $!;
	foreach my $block (%BL){
		next unless $BL{$block} eq $chr;
		open HAP0, "> $asmDir/$block.hap0.config" or die $1;
		print HAP0 "#maximal read length
max_rd_len=100
[LIB]
#average insert size
avg_ins=200
##if sequence needs to be reversed
reverse_seq=0
##in which part(s) the reads are used
asm_flags=3
##use only first 100 bps of each read
rd_len_cutoff=100
##in which order the reads are used while scaffolding
rank=1
## cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=5
##minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=35
##a pair of fastq file, read 1 file should always be followed by read 2 file
q1=$outDir/$block.hap0.read.1.fq
q2=$outDir/$block.hap0.read.2.fq
";
		close HAP0;
		open HAP1, "> $asmDir/$block.hap1.config" or die $1;
	    print HAP1 "#maximal read length
max_rd_len=100
[LIB]
#average insert size
avg_ins=200
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#use only first 100 bps of each read
rd_len_cutoff=100
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=5
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=35
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=$outDir/$block.hap1.read.1.fq
q2=$outDir/$block.hap1.read.2.fq
";
		close HAP1;
		print ASM "$root_path/SOAPdenovo2/SOAPdenovo-63mer all -s $asmDir/$block.hap0.config -K 35 -M 0 -o $asmDir/$block.hap0 1>$asmDir/$block.hap0.log 2>$asmDir/$block.hap0.error\n";
		print ASM "$root_path/SOAPdenovo2/SOAPdenovo-63mer all -s $asmDir/$block.hap1.config -K 35 -M 0 -o $asmDir/$block.hap1 1>$asmDir/$block.hap1.log 2>$asmDir/$block.hap1.error\n";
		print ASM "$root_path/MUMmer-3.23/dnadiff $refGenome -p $asmDir/$block.hap0 $asmDir/$block.hap0.scafSeq\n";
		print ASM "$root_path/MUMmer-3.23/dnadiff $refGenome -p $asmDir/$block.hap1 $asmDir/$block.hap1.scafSeq\n";
	}
	close ASM;
}

