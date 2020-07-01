#!/usr/bin/perl
# /data/projects/phillippy/projects/rDNA/getCompleteFrequenciesRibosomeForNewVariants.pl

use Getopt::Long;
use strict;
use Data::Dumper;

my $samtools_bin = 'samtools';
my $BAM;
my $outputFile;
my $minMappingQuality = 30;
my $minBaseQuality = 35;
my $referenceGenome;

GetOptions (
	'BAM:s' => \$BAM,
	'outputFile:s' => \$outputFile,
	'minMappingQuality:s' => \$minMappingQuality,
	'minBaseQuality:s' => \$minBaseQuality,
	'referenceGenome:s' => \$referenceGenome,
);

unless(-e $BAM)
{
	die "--BAM file $BAM not found";
}

unless($outputFile)
{
	die "Please specify --outputFile"
}
my $BAMi = $BAM.'.bai';
unless(-e $BAMi)
{
	die "BAM index file $BAMi not found";
}

print "Parameters:\n";
print "\t", "BAM", ": ", $BAM, "\n";
print "\t", "minMappingQuality", ": ", $minMappingQuality, "\n";
print "\t", "minBaseQuality", ": ", $minBaseQuality, "\n";

my %BAM_regions_lengths;
my $idxstats_command = qq($samtools_bin idxstats $BAM);
my $idxstats_output = `$idxstats_command`;
foreach my $line (split(/\n/, $idxstats_output))
{
	die unless($line =~ /^(\S+)\t(\d+)/);
	my $id = $1;
	next if($id eq '*');
	$BAM_regions_lengths{$id} = $2;
}	
die unless(scalar(keys %BAM_regions_lengths) == 1);
my $regionName = (keys %BAM_regions_lengths)[0];

# die "Weird BAM -- RIBOSOME region is not length 44838" unless($BAM_regions_lengths{"RIBOSOME"} == 44838);

my $outputFile3 = $outputFile;
open(OUTPUT3, '>', $outputFile3) or die "Cannot open $outputFile3";
# print OUTPUT3 join("\t", qw/ChromosomeID Position N A C G T OtherAlleles/), "\n";

my $samtools_command = qq(samtools mpileup -d 100000 -r $regionName --fasta-ref $referenceGenome -q$minMappingQuality -Q$minBaseQuality $BAM);
# die $samtools_command;
# $samtools_command = qq($samtools_bin mpileup -d 50000 -r RIBOSOME -q$minMappingQuality -Q$minBaseQuality $BAM);
print "Now executing:\n", $samtools_command, "\n";

my %alleles_nextColumn;
my $alleles_nextColumn_position;

my %deletionsCarryForward;
my %interestingPositions;
my %supplementaryPositions; 
my %allele_frequencies;
open(FREQ, "$samtools_command |") or die "Cannot open pipe $samtools_command";
while(<FREQ>)
{
	my $line = $_;
	chomp($line);
	
	my @fields = split(/\t/, $line);
	
	my $chr = $fields[0];
	my $position = $fields[1];
	my $refbase = $fields[2];
	my $coverage = $fields[3];
	my $alleles = $fields[4];	

	my $lookupKey = key_for_interestingPosition($chr, $position);

	
	my @alleles;
	for(my $i = 0; $i < length($alleles); $i++)
	{
		my $C = substr($alleles, $i, 1);
		if($C eq '^')
		{
			$i++;
			next;
		}
		next if($C eq '$');
		
		if(($C eq '+') or ($C eq '-'))
		{
			my $j = 1;
			while(substr($alleles, $i + $j + 1, 1) =~ /^\d$/)
			{
				$j++;
			}
			my $lengthInsertion = substr($alleles, $i+1, $j);
			die unless($lengthInsertion =~ /^\d+$/);
			my $allele = substr($alleles, $i + 1 + $j, $lengthInsertion);
			if(substr($alleles, $i, 1) eq '+')
			{
				$alleles[$#alleles] .= $allele;
			}
			else
			{
				# push(@alleles, '-');
				my $deletionAllele = '-' x $lengthInsertion;
				die unless(length($deletionAllele) == $lengthInsertion);
				$alleles[$#alleles] .= $deletionAllele;	
				my $isPlusStrand;
				if($allele =~ /^[ACGTN]+$/)
				{
					$isPlusStrand = 1;
				}	
				else
				{
					die Dumper("Weird deletionAllele", $allele, $alleles) unless ($allele =~ /^[acgtn]+$/);
					$isPlusStrand = 0;
				}
				die unless(defined $isPlusStrand);
				
				die unless($lengthInsertion >= 1);
				for(my $refPosII = $position+1; $refPosII <= $position + $lengthInsertion; $refPosII++)
				{
					my $lookupkey_carryForward = key_for_interestingPosition($chr, $refPosII);
					my $allele_strand = ($isPlusStrand) ? '(' : ')';
					$deletionsCarryForward{$lookupkey_carryForward}{$allele_strand}++;
				}
				#$next_alleles_nextColumn{$deletionAllele}++;
			}
			$i = ($i + 1 + $j + $lengthInsertion) - 1;
		}
		elsif(($C eq '>') or ($C eq '<') or ($C eq '*'))
		{
			if(($C eq '>') or ($C eq '<'))
			{
				die Dumper("This is unexpected, position $position of BAM $BAM");
			}
			push(@alleles, '-');		
		}
		else
		{
			my $a = substr($alleles, $i, 1);
			if($a eq '.')
			{
				$a = uc($refbase);
			}
			elsif($a eq ',')
			{
				$a = lc($refbase);
			}
			
			die Dumper("Weird allele $a") unless($a =~ /^[ACGTN]$/i);
			push(@alleles, $a);
		}	
	}	
	
	my %alleles_counts;
	foreach my $allele (@alleles)
	{
		next if($allele eq '-');
		$alleles_counts{$allele}++;
	}
	if(exists $deletionsCarryForward{$lookupKey})
	{
		foreach my $allele (keys %{$deletionsCarryForward{$lookupKey}})
		{
			$alleles_counts{$allele} += $deletionsCarryForward{$lookupKey}{$allele};			
		}
	}
	
	$allele_frequencies{$lookupKey} = \%alleles_counts;

	print OUTPUT3 join("\t", $position, $refbase, join(';', map {$_ . '=' . $alleles_counts{$_}} keys %alleles_counts)), "\n";
}
close(FREQ);	

print "\nDone. Output file: $outputFile3 \n";

sub key_for_interestingPosition
{
	my $chr = shift;
	my $position = shift;
	return $position;
}

sub interestingPositions_from_key
{
	my $key = shift;
	return ('RIBOSOME', $key);
}

sub unique
{
	my @input = @_;
	my %_input = map {$_ => 1} @input;
	return keys %_input;
}
