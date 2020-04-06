# HTS-Metabarcoding-Primer-Trimming-Scripts
Perl scripts to trim the Copepod, Insect, or Mollusc primer and spacer regions from HTS metabarcoding data.  Marshall, N, &amp; Stepien, C. (2020). Macroinvertebrate community diversity and habitat quality relationships along a large river from targeted eDNA metabarcoding assays.

***TrimSpacePrimer*** is a stand alone perl script that can be called from a perl interpreter. 

## Using *TrimSpacePrimer* ##

***TrimSpacePrimer*** trims paired end metabarcoding reads returned from Illumina sequencing for subsequent merging in other programs such as Dada2, Unoise, or OBITools. ***TrimSpacePrimer*** works on the parent directory level. It should be run from the parent directory in which all of your subdirectories with Illumina sequencing results are stored. This allows you to trim individual sets of paired sequencing results files. ***TrimSpacePrimer*** handles zipped fastq.gz files. It will deposit your trimmed meta-barcoding results as unzipped fastq files into a new subdirectory named with your parent directory name. If you have used primers with identical spacer inserts to those published in Klymus *et al.* 2017, *Plos One*, **12**(5): e0177643, you can use those inserts to remove instances of index hops. ***All inputs are case INsensitive! fastq.gz files must use the demultiplexed illumina naming convention (forward read file ends in R1.001.fastq.gz, reverse read file ends in R2_001.fastq.gz)***


Which script to use is based on the primerset used for metabarcoding:
Copepod 16S (Clarke et al., 2017) – CopTrimSpacePrimer.pl
Insect 16S (Epp et al., 2012) – InsTrimSpacePrimer.pl
Mollusk 16S (Klymus et al., 2017) – MolTrimSpacePrimer.pl

***TrimSpacePrimer*** identifies which spacer region was used for each sample by counting the times each spacer appears within reads for each sample, and the spacer with the highest occurrence is considered the used spacer (see Klymus et al., 2017). Any sequences without the considered appropriate spacer region are removed from the dataset, which allows the removal of index hops! 
SpacerE: TCCTATG + CGTACTAGATGTACGA
SpacerF: ATGCTACAGT + TCACTAGCTGACGC
SpacerG: CGAGGCTACAACTC + GAGTAGCTGA
SpacerH: GATACGATCTCGCACTC + ATCGGCT

Next the reverse primer is identified and trimmed from the reverse read file and the forward primer is trimmed from the forward read file.  If the forward or reverse primers are not found, both the forward and reverse reads are removed from the dataset.
Copepod 16S: TAAGGTAGCATARTAATTWG + TAATTCAACATCGAGGTC
Insect 16S: TGCAAAGGTAGCATAATMATTAG + TCCATAGGGTCTTCTCGTC
Mollusk 16S: RRWRGACRAGAAGACCCT + ARTCCAACATCGAGGT

Finally, the reads are trimmed to a specified length for each marker.  The ***InsTrimSpacePrimer*** removes both the forward and reverse primers from both the forward and reverse reads.  The ***CopTrimSpacePrimer*** and ***MolTrimSpacePrimer*** yields amplicons of varying length, resulting in some reads only having the forward or reverse primer.  Thus these scripts trim reads to a length of 185bp or 125bp based on their length following primer removal.
