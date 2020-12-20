#! /usr/bin/perl
use strict;
use FindBin qw($Bin);
use Cwd qw(abs_path getcwd);

my $root_path = $Bin;
my $abs = abs_path (getcwd());
my $len_file = $ARGV[0];
my $abs = abs_path (getcwd());

die "Usage:perl $0 <*.fa.fai> <indexed sorted bam> <vcf> <window size> <max LFR> <min SupportBarcodes> <min Link> <shellDir> <outDir_tmp> <outDir_final> <fq1> <fq2> <outDir_reads> <asmDir> <ref Genome>
Example:perl $0 ref.fa.fai input.bam input.vcf 10000000 300000 1 1 shDir outDir_tmp outDir_final fq1 fq2 outDir_reads asmDir refGenome

This script is writen NOT only for individual genome phasing based on hg19 reference genome, but also for any diploid species phasing based on any given reference genome. Besides of phasing function, users can also achieve the haplotype-resolved assembly as well as the corresponding variant calls, based on the phased reults.
<window size> is set for each phasing window, 10MB is suggested.
<max LFR> is decided by the maximum length of DNA fragment
<min SupportBarcodes> is minimum count of supporting barcode for each variant
<min Link> is minimum count of linked info

Author: Sun Yuhui
Email: sunyuhui\@genomics.cn\n" unless @ARGV==15;


my $bam = $ARGV[1];
unless ($bam =~ /^\//){
	$bam = "$abs"."/"."$bam";
}
#$bam = abs_path ($bam);
#

my $hg19 = 1; 
my $vcf = $ARGV[2];
unless ($vcf =~ /^\//){
	$vcf = "$abs"."/"."$vcf";
}
open VCF, $vcf or die $!;
while(<VCF>){
	chomp;
	next if /^#/;
	my $line = $_;
	my @info = split /\s+/, $line;
	my $chr = $info[0];
	if ($chr =~ /^chr/){
		$hg19 = 1;
	}else{
		$hg19 = 0;
	}
	last;
}
close VCF;

#$vcf = abs_path ($vcf);
my $win_size = $ARGV[3];
my $max_LFR_len = $ARGV[4];
my $min_barcode_num = $ARGV[5];
my $min_link = $ARGV[6];


my $shellDir = $ARGV[7];
$shellDir = abs_path ($shellDir);
mkdir $shellDir unless -d $shellDir;

my $outDir = $ARGV[8];
$outDir = abs_path ($outDir);
mkdir $outDir unless -d $outDir;

my $mergeDir = $ARGV[9];
$mergeDir = abs_path ($mergeDir);
mkdir $mergeDir unless -d $mergeDir;

my ($fq1, $fq2, $groupReadsDir, $asmDir) = @ARGV[10,11,12,13];
$groupReadsDir = abs_path $groupReadsDir;
$asmDir = abs_path $asmDir;
mkdir $groupReadsDir unless -d $groupReadsDir;
mkdir $asmDir unless -d $asmDir;

my $refGenome = $ARGV[14];
$refGenome = abs_path $refGenome;

my $size = 10000000;
open LEN, $len_file or die $!;
while(<LEN>){
	chomp;
	my ($chr, $len) = (split /\s+/, $_)[0,1];
#	next if /_/;
#	next unless /chr[\d+|X]/; #TODO
	my $num = int($len/$size)+1;
	my ($beg, $end);
	my $rename_chr = $chr;
#	if ($hg19 eq 0){ #TODO
#		$rename_chr =~ s/chr//;
#	}
	my $asmCom;
	open SH, ">$shellDir/run.$chr.sh" or die $!;
	for my $i (1..$num){
		if ($i eq 1){
			$beg = 1;
			$end = $size;
		}elsif ($i >1 and $i <$num){
			$beg = ($i-1)*$size+1-500000;
			$end = $i *$size;
		}elsif ($i eq $num){
			$beg = ($i-1)*$size+1-500000;
			$end = $len;
		}
		print SH "perl $root_path/Bin/LongHap.pl $bam $vcf $rename_chr:$beg-$end $max_LFR_len $min_barcode_num $min_link $outDir/$chr\_$beg\_$end.log 2>$outDir/$chr\_$beg\_$end.error\n";

	}
	print SH "perl $root_path/Bin/merge.pl $outDir $chr $win_size $mergeDir/$chr.log\n";
	print SH "cat $outDir/$rename_chr\_*.hete.barcodes |sort |uniq >$mergeDir/$chr.hete.barcodes\n";
	print SH "cat $outDir/$rename_chr\_*.homo.barcodes |sort |uniq >$mergeDir/$chr.homo.barcodes\n";
	print SH "cat $outDir/$rename_chr\_*.unknown.barcodes |sort |uniq >$mergeDir/$chr.unknown.barcodes\n";
	print SH "perl $root_path/Bin/log2vcf.pl $mergeDir/$chr.log $mergeDir/$chr.hete.barcodes $mergeDir/$chr.homo.barcodes $vcf >$mergeDir/$chr.vcf\n";
	print SH "perl $root_path/Bin/groupReadsByBlock.pl $mergeDir/$chr.vcf $fq1 $fq2 $groupReadsDir $asmDir $shellDir $refGenome\n";
	print SH "sh $shellDir/soapdenovo.$chr.sh";
	close SH;
}
close LEN;

