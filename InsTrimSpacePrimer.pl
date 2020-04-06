#!user/local/perl
#Created by C. Pells, M. R. Snyder, and N. T. Marshall 2017

#Script trims spacers and primers from the fastq files of high throughput sequencing for a specific primer set

#You must run this script in the directory with your sample files.  The sample files need to be unzipped fastq files in their own directories.  For example, if I was running this script to trim the primers from sampels mSite01 and mSite02, I would have a directory with three items in it.  1) this script, 2) a directory labeled 'mSite01' which has the associated unzipped forward and reverse fastq files, and 3) a directory labeled 'mSite02' which has the associated unzipped forward and reverse fastq files


use Cwd;
use File::Basename;
#use warnings;
use List::Util "max";

$StartTime= localtime;

$MasterDir = getcwd; #obtains the current directory
$DirSmall = basename ($MasterDir);

opendir (DIR, $MasterDir);
@objects = readdir (DIR);
closedir (DIR);
foreach (@objects){
	print $_,"\n";
}

$DadaDir = $MasterDir."/"."$DirSmall"."Dada2";
mkdir $DadaDir unless -d $DadaDir;

$OutSumName= "$DirSmall"."DadaTrimSummary.txt";
open (OUTS, ">", "$OutSumName" || die "cannot open Summary.txt\n");
print OUTS "Sample\tSpacer\tTotal\tForward w/ correct spacer\tReverse w/ correct spacer\tF w/ FPrimer\tR w/ RPrimer\tF w/ RPrimer\tR w/ FPrimer\tF correct len\tR correct len\tF incorrect len\tR incorrect len\tReads w/ both primers correct len\n";

@Dirs = ();
@DirsSmall = ();
foreach $O (0..$#objects){
	$CurrDir = "";
	if ($objects[$O] =~ /(^[a-zA-Z]{4}[0-9]{2})/){ #searches directories for samples, they need to be labeled as marker then site (i.e., mSite01).  You can change this code to search for a different naming system if you would like.
        	$CurrDir = $MasterDir."/".$objects[$O]; #appends directory name to full path
        	push (@Dirs, $CurrDir);
        	push (@DirsSmall, substr($objects[$O], 0, 6));
    	}
}

foreach (@Dirs){
	print $_,"\n";#checks that all directories were read in
}

foreach $S (0..$#Dirs){
	@files = ();
        opendir (DIR, $Dirs[$S]) || die "cannot open $Dirs[$S]: $!";
        @files = readdir DIR; #reads in all files in a directory
        closedir DIR;
        @AbsFiles = ();
        foreach $F (0..$#files){
        	$AbsFileName = $Dirs[$S]."/".$files[$F]; #appends file name to full path
        	push (@AbsFiles, $AbsFileName);
    	}
        %RSeqHash;
        %RQualHash;
        %RNameHash;
        $eSpacer = "ATGTACGAA|[ACGTN]TGTACGAA|A[ACGTN]GTACGAA|AT[ACGTN]TACGAA|ATG[ACGTN]ACGAA|ATGT[ACGTN]CGAA|ATGTA[ACGTN]GAA|ATGTAC[ACGTN]AA|ATGTACG[ACGTN]A|ATGTACGA[ACGTN]";
        $fSpacer = "GCTGACGCA|[ACGTN]CTGACGCA|G[ACGTN]TGACGCA|GC[ACGTN]GACGCA|GCT[ACGTN]ACGCA|GCTG[ACGTN]CGCA|GCTGA[ACGTN]GCA|GCTGAA[ACGTN]CA|GCTGACG[ACGTN]A|GCTGACGC[ACGTN]";
    	$gSpacer = "GTAGCTGAA|[ACGTN]TAGCTGAA|G[ACGTN]AGCTGAA|GT[ACGTN]GCTGAA|GTA[ACGTN]CTGAA|GTAG[ACGTN]TGAA|GTAGC[ACGTN]GAA|GTAGCT[ACGTN]AA|GTAGCTG[ACGTN]A|GTAGCTGA[ACGTN]";
        $hSpacer = "ATCGGCTA|[ACGTN]TCGGCTA|A[ACGTN]CGGCTA|AT[ACGTN]GGCTA|ATC[ACGTN]GCTA|ATCG[ACGTN]CTA|ATCGG[ACGTN]TA|ATCGGC[ACGTN]A|ATCGGCT[ACGTN]";
        $e = 0;
        $f = 0;
        $g = 0;
        $h = 0;
	$eSp = "TCCTATG|[ACGTN]CCTATG|T[ACGTN]CTATG|TC[ACGTN]TATG|TCC[ACGTN]ATG|TCCT[ACGTN]TG|TCCTA[ACGTN]G|TCCTAT[ACGTN]"; 
        $fSp = "GCTACAGT|[ACGTN]CTACAGT|G[ACGTN]TACAGT|GC[ACGTN]ACAGT|GCT[ACGTN]CAGT|GCTA[ACGTN]AGT|GCTAC[ACGTN]GT|GCTACA[ACGTN]T|GCTACAG[ACGTN]";
        $gSp = "TACAACTC|[ACGTN]ACAACTC|T[ACGTN]CAACTC|TA[ACGTN]AACTC|TAC[ACGTN]ACTC|TACA[ACGTN]CTC|TACAA[ACGTN]TC|TACAAC[ACGTN]C|TACAACT[ACGTN]";
        $hSp = "TCGCACTC|[ACGTN]CGCACTC|T[ACGTN]GCACTC|TC[ACGTN]CACTC|TCG[ACGTN]ACTC|TCGC[ACGTN]CTC|TCGCA[ACGTN]TC|TCGCAC[ACGTN]C|TCGCACT[ACGTN]";
     	foreach $AF (0..$#AbsFiles){
         	if (($AbsFiles[$AF] =~ m/_R2_001\.fastq$/) and ($AbsFiles[$AF] !~ m/\.gz/)) { #finds reverse fastq file
            		@readbuffer=();
            		#read in reverse fastq
            		open (INPUT2, $AbsFiles[$AF]) || die "Can't open file: $!\n";
            		$RTrimFastq = "$DadaDir"."/"."$DirsSmall[$S]"."_XXXX_R2_001.fastq";
            		open (OUTR, ">", "$RTrimFastq") || die "Can't open file: $!\n";
            		$x=0;
            		$xB=0;
            		$xC=0;
            		$xD=0;
            		$xE=0;
			$CS=0;
            		$multi=0;
            		$None=0;
            		while (<INPUT2>) {
               			chomp ($_);
                		push(@readbuffer, $_);
                    		if (@readbuffer == 4) {
                        		$cc=0;
                        		@splitSeq=split(//, $readbuffer[1]);
                        		@push=();
                        		for $xx(0..49){
                                		push(@push, $splitSeq[$xx]);
                        		}
                        		$FindSpace=join('', @push);
                        		if ($FindSpace =~ /^$gSpacer/){
                            			$g++;
						$cc++;
                        		}
                        		if ($FindSpace =~ /^$fSpacer/) {
                            			$f++;
						$cc++;
                        		}
                        		if ($FindSpace =~ /^$eSpacer/) {
                            			$e++;
						$cc++;
                        		}
                        		if ($FindSpace =~ /^$hSpacer/) {
                            			$h++;
						$cc++;
                        		}	
                        		@readbuffer = ();
                        		if ($cc==0){
                            			$None++;
                        		}
                        		if ($cc>1){
                            			$multi++;
                        		}
                   		}
            		}
            		print "E = $e\nF = $f\nG = $g\nH = $h\nNone = $None\nMultiple Spacer = $multi\n";
            		$CorrSpacer = "";
            		$maxSpacer = max ($g, $e, $f, $h);
            		print "Max spacer = $maxSpacer\n";
            		if ($maxSpacer == $e) {
                		$CorrSpacer = $eSpacer;
				$CorrFSpacer = $eSp;
				$Space="E";
                		print "$DirsSmall[$S] has e spacer\n"
            		}
            		elsif ($maxSpacer == $f) {
                		$CorrSpacer = $fSpacer;
				$CorrFSpacer = $fSp;
                                $Space="F";
                		print "$DirsSmall[$S] has f spacer\n"
            		}
            		elsif ($maxSpacer == $g) {
                		$CorrSpacer = $gSpacer;
				$CorrFSpacer = $gSp;
                                $Space="G";
                		print "$DirsSmall[$S] has g spacer\n"
            		}
            		elsif ($maxSpacer == $h) {
                		$CorrSpacer = $hSpacer;
				$CorrFSpacer = $hSp;
                                $Space="H";
                		print "$DirsSmall[$S] has h spacer\n"
            		}
            		print "Spacer is $CorrSpacer\n";
            		seek INPUT2, 0, 0;
           		while (<INPUT2>){
                		chomp ($_);
                		push(@readbuffer, $_);
                		if (@readbuffer == 4) {
                    			$x++;
                    			@splitSeq=split(//, $readbuffer[1]);
                    			for $xx(0..49){
                        			push(@push, $splitSeq[$xx]);
                    			}
                    			$FindSpace=join('', @push);
					@push=();
                    			if ($FindSpace =~ /^$CorrSpacer/){
						$CS++;
                        			if ($readbuffer[1] =~ /TCCATAGGGTCTTCTCGTC|[NACGT]CCATAGGGTCTTCTCGTC|T[NACGT]CATAGGGTCTTCTCGTC|TC[NACGT]ATAGGGTCTTCTCGTC|TCC[NACGT]TAGGGTCTTCTCGTC|TCCA[NACGT]AGGGTCTTCTCGTC|TCCAT[NACGT]GGGTCTTCTCGTC|TCCATA[NACGT]GGTCTTCTCGTC|TCCATAG[NACGT]GTCTTCTCGTC|TCCATAGG[NACGT]TCTTCTCGTC|TCCATAGGG[NACGT]CTTCTCGTC|TCCATAGGGT[NACGT]TTCTCGTC|TCCATAGGGTC[NACGT]TCTCGTC|TCCATAGGGTCT[NACGT]CTCGTC|TCCATAGGGTCTT[NACGT]TCGTC|TCCATAGGGTCTTC[NACGT]CGTC|TCCATAGGGTCTTCT[NACGT]GTC|TCCATAGGGTCTTCTC[NACGT]TC|TCCATAGGGTCTTCTCG[NACGT]C|TCCATAGGGTCTTCTCGT[NACGT]/){	
							$xB++;
                        				$SeqEnd=length($readbuffer[1]);
                        				$SeqStart=($-[0]) + 19;
							if ($readbuffer[1] =~ /CTAAT[TG]ATTATGCTACCTTTGCA|[NACGT]TAAT[TG]ATTATGCTACCTTTGCA|C[NACGT]AAT[TG]ATTATGCTACCTTTGCA|CT[NACGT]AT[TG]ATTATGCTACCTTTGCA|CTA[NACGT]T[TG]ATTATGCTACCTTTGCA|CTAA[NACGT][TG]ATTATGCTACCTTTGCA|CTAAT[NACGT]ATTATGCTACCTTTGCA|CTAAT[TG][NACGT]TTATGCTACCTTTGCA|CTAAT[TG]A[NACGT]TATGCTACCTTTGCA|CTAAT[TG]AT[NACGT]ATGCTACCTTTGCA|CTAAT[TG]ATT[NACGT]TGCTACCTTTGCA|CTAAT[TG]ATTA[NACGT]GCTACCTTTGCA|CTAAT[TG]ATTAT[NACGT]CTACCTTTGCA|CTAAT[TG]ATTATG[NACGT]TACCTTTGCA|CTAAT[TG]ATTATGC[NACGT]ACCTTTGCA|CTAAT[TG]ATTATGCT[NACGT]CCTTTGCA|CTAAT[TG]ATTATGCTA[NACGT]CTTTGCA|CTAAT[TG]ATTATGCTAC[NACGT]TTTGCA|CTAAT[TG]ATTATGCTACC[NACGT]TTGCA|CTAAT[TG]ATTATGCTACCT[NACGT]TGCA|CTAAT[TG]ATTATGCTACCTT[NACGT]GCA|CTAAT[TG]ATTATGCTACCTTT[NACGT]CA|CTAAT[TG]ATTATGCTACCTTTG[NACGT]A|CTAAT[TG]ATTATGCTACCTTTGC[NACGT]/){ #Compliment forward primer	
								$SeqEnd = ($-[0]) - $SeqStart;
                            					$xC++;
							}
                            				@Name = split(/\s/,$readbuffer[0]);
                            				$rsn = $Name[0];
                            				$RRCS = substr($readbuffer[1], $SeqStart, $SeqEnd);
                            				$lenR = length ($RRCS);
                            				if ($lenR > 50){
                                				$xD++;
                                				$RNameHash{$rsn} = $readbuffer[0];
                                				$RSeqHash{$rsn} = $RRCS;
                                				$RQualHash{$rsn} = substr($readbuffer[3], $SeqStart, $SeqEnd);
                            				}
                    					else {
                        					$xE++;
                    					}
                				}
					}
        				@readbuffer = ();
					#last;
				}
    			}
		}
    		foreach $AFx (0..$#AbsFiles){
        		if (($AbsFiles[$AFx] =~ m/_R1_001\.fastq$/) and ($AbsFiles[$AFx] !~ m/\.gz/)) { #finds forward fastq file
				open (INPUT1, $AbsFiles[$AFx]) || die "Can't open file: $!\n";
            			$FTrimFastq = "$DadaDir"."/"."$DirsSmall[$S]"."_XXXX_R1_001.fastq";
            			open (OUTF, ">", "$FTrimFastq") || die "Can't open file: $!\n";
            			@readbuffer = ();
            			$y=0;
            			$yB=0;
            			$yC=0;
				$yD=0;
				$yE=0;
            			$q=0;
				$FCS=0;
            			while (<INPUT1>){
                			chomp ($_);
                			push(@readbuffer, $_);
                			if (@readbuffer == 4) {
                    				$y++;
						@splitSeq=split(//, $readbuffer[1]);
                                                for $xx(0..49){
                                                        push(@push, $splitSeq[$xx]);
                                                }
                                                $FindSpace=join(//, @push);
                                                @push=();
                                                if ($FindSpace =~ /^$CorrFSpacer/){
							$FCS++;
                                                        $SeqEnd=length($readbuffer[1]);
							if ($readbuffer[1] =~ /TGCAAAGGTAGCATAAT[AC]ATTAG|[ACGTN]GCAAAGGTAGCATAAT[AC]ATTAG|T[ACGTN]CAAAGGTAGCATAAT[AC]ATTAG|TG[ACGTN]AAAGGTAGCATAAT[AC]ATTAG|TGC[ACGTN]AAGGTAGCATAAT[AC]ATTAG|TGCA[ACGTN]AGGTAGCATAAT[AC]ATTAG|TGCAA[ACGTN]GGTAGCATAAT[AC]ATTAG|TGCAAA[ACGTN]GTAGCATAAT[AC]ATTAG|TGCAAAG[ACGTN]TAGCATAAT[AC]ATTAG|TGCAAAGG[ACGTN]AGCATAAT[AC]ATTAG|TGCAAAGGT[ACGTN]GCATAAT[AC]ATTAG|TGCAAAGGTA[ACGTN]CATAAT[AC]ATTAG|TGCAAAGGTAG[ACGTN]ATAAT[AC]ATTAG|TGCAAAGGTAGC[ACGTN]TAAT[AC]ATTAG|TGCAAAGGTAGCA[ACGTN]AAT[AC]ATTAG|TGCAAAGGTAGCAT[ACGTN]AT[AC]ATTAG|TGCAAAGGTAGCATA[ACGTN]T[AC]ATTAG|TGCAAAGGTAGCATAA[ACGTN][AC]ATTAG|TGCAAAGGTAGCATAAT[ACGTN]ATTAG|TGCAAAGGTAGCATAAT[AC][ACGTN]TTAG|TGCAAAGGTAGCATAAT[AC]A[ACGTN]TAG|TGCAAAGGTAGCATAAT[AC]AT[ACGTN]AG|TGCAAAGGTAGCATAAT[AC]ATT[ACGTN]G|TGCAAAGGTAGCATAAT[AC]ATTA[ACGTN]/){ #Forward primer
								$SeqStart=($-[0]) + 23;
								$yB++;
								$SeqEnd=length($readbuffer[1]);
								if ($readbuffer[1] =~ /GACGAGAAGACCCTATGGA|[NACGT]ACGAGAAGACCCTATGGA|G[NACGT]CGAGAAGACCCTATGGA|GA[NACGT]GAGAAGACCCTATGGA|GAC[NACGT]AGAAGACCCTATGGA|GACG[NACGT]GAAGACCCTATGGA|GACGA[NACGT]AAGACCCTATGGA|GACGAG[NACGT]AGACCCTATGGA|GACGAGA[NACGT]GACCCTATGGA|GACGAGAA[NACGT]ACCCTATGGA|[NACGT]CCCTATGGA|GACGAGAAGA[NACGT]CCTATGGA|GACGAGAAGAC[NACGT]CTATGGA|GACGAGAAGACC[NACGT]TATGGA|GACGAGAAGACCC[NACGT]ATGGA|GACGAGAAGACCCT[NACGT]TGGA|GACGAGAAGACCCTA[NACGT]GGA|GACGAGAAGACCCTAT[NACGT]GA|GACGAGAAGACCCTATG[NACGT]A|GACGAGAAGACCCTATGG[NACGT]/){ #Compliment Reverse primer
									$SeqEnd = ($-[0]) - $SeqStart;
                        						$yC++;
								}
                        					@Name = split(/\s/,$readbuffer[0]);
                        					$fsn = $Name[0]; #trims forward seq name
                        					$SeqFTrim = substr($readbuffer[1], $SeqStart, $SeqEnd);
                        					$lenF = length ($SeqFTrim);
                        					$QualFTrim = substr($readbuffer[3], $SeqStart, $SeqEnd);
                        					if ($lenF > 50){
                            						$yD++;
                            						if (exists($RSeqHash{$fsn})){ #checks to see if forward seq name is present in reverse seq hash
                                						print OUTR "$RNameHash{$fsn}\n$RSeqHash{$fsn}\n+\n$RQualHash{$fsn}\n";
                                						print OUTF "$readbuffer[0]\n$SeqFTrim\n+\n$QualFTrim\n";
                                						$q++
                            						}
                        					}
								else{
									$yE++;
								}
                    					}
						}
                				@readbuffer = ();
                			}
            			}
            			close INPUT1;
            			close INPUT2;
            			close OUTF;
            			close OUTR;
            			last;
        		}
    		}
	}
    	$CurrSamp = $DirsSmall[$S];
    	print OUTS "$CurrSamp\t$Space\t$y\t$FCS\t$CS\t$yB\t$xB\t$yC\t$xC\t$yD\t$xD\t$yE\t$xE\t$q\n";
}






