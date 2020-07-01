# original location: /gpfs/project/dilthey/projects/COVID-19/generate_IGV.pl

use strict;
use Data::Dumper;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename;
use File::Spec;

my $COVID_ref;
my $outputFile;
my $BAM;
my $VCF_Medaka;
my $VCF_Nanopolish;

GetOptions (
	'covidReference:s' => \$COVID_ref,
	'outputFile:s' => \$outputFile,
	'BAM:s' => \$BAM,
	'VCF_Medaka:s' => \$VCF_Medaka,
	'VCF_Nanopolish:s' => \$VCF_Nanopolish,
);

die "Please specify --outputFile" unless($outputFile);
die "Please specify --BAM" unless($BAM);
die "Please specify --VCF_Medaka" unless($VCF_Medaka);
die "Please specify --VCF_Nanopolish" unless($VCF_Nanopolish);

die "File --BAM not existing" unless(-e $BAM);
die "File --VCF_Medaka not existing" unless(-e $VCF_Medaka);
die "File --VCF_Nanopolish not existing" unless(-e $VCF_Nanopolish);


#die "Path of --BAM ($BAM_path) not equal to path of --outputFile ($outputFile_path)" unless($BAM_path eq $outputFile_path);
#die "Path of --VCF_Medaka not equal to path of --outputFile" unless($VCF_Medaka_path eq $outputFile_path);
#die "Path of --VCF_Nanopolish not equal to path of --outputFile" unless($VCF_Nanopolish_path eq $outputFile_path);

my $outputDir = dirname($outputFile);

my $BAM_N = File::Spec->abs2rel(abs_path($BAM), $outputDir);
my $VCFM_N = File::Spec->abs2rel(abs_path($VCF_Nanopolish), $outputDir);
my $VCFN_N = File::Spec->abs2rel(abs_path($VCF_Medaka), $outputDir);
my $REF_N = File::Spec->abs2rel(abs_path($COVID_ref), $outputDir);

open(XML, '>', $outputFile) or die "Cannot open $outputFile";
print XML qq(<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="${REF_N}" hasGeneTrack="false" hasSequenceTrack="true" locus="MN908947.3:1-29903" version="8">
    <Resources>
        <Resource path="${VCFM_N}"/>
        <Resource path="${VCFN_N}"/>
        <Resource path="${BAM_N}"/>
    </Resources>
    <Panel height="84" name="DataPanel" width="1262">
        <Track clazz="org.broad.igv.variant.VariantTrack" color="0,0,178" displayMode="EXPANDED" fontSize="10" id="${VCFN_N}" name="${VCFN_N}" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>
        <Track clazz="org.broad.igv.variant.VariantTrack" color="0,0,178" displayMode="EXPANDED" fontSize="10" id="${VCFM_N}" name="${VCFM_N}" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>
    </Panel>
    <Panel height="13302" name="Panel1583780401168" width="1262">
        <Track autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" color="175,175,175" colorScale="ContinuousColorScale;0.0;937.0;255,255,255;175,175,175" fontSize="10" id="${BAM_N}_coverage" name="${BAM_N} Coverage" snpThreshold="0.2" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="937.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track clazz="org.broad.igv.sam.SpliceJunctionTrack" fontSize="10" height="60" id="${BAM_N}_junctions" name="${BAM_N} Junctions" visible="false"/>
        <Track clazz="org.broad.igv.sam.AlignmentTrack" displayMode="EXPANDED" experimentType="OTHER" fontSize="10" id="${BAM_N}" name="${BAM_N}" visible="true">
            <RenderOptions/>
        </Track>
    </Panel>
    <Panel height="44" name="FeaturePanel" width="1262">
        <Track clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" visible="true"/>
    </Panel>
    <PanelLayout dividerFractions="0.177319587628866,0.8989690721649485"/>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>
);
close(XML);
