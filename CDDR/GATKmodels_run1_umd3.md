# Variant filtration using GATK VQSR #


## Trial 1 ##

The training and truth set is the high quality filtered dataset with 218714 SNPs at 

/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz

The training false dataset is the dataset which includes all the above variants along with false positive variants that are heterozygous in the non-target animals (/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/gatk_datasets/falseset/falseset.vcf.gz) 


Here is the SNP count per chromosome in each of the datasets:


**CHROM**|**falseset**|**selected.snps**
:-----:|:-----:|:-----:
chr1|157627|13814
chr2|131249|12177
chr3|92082|8411
chr4|56090|3991
chr5|90899|7391
chr6|66420|5032
chr7|139162|11316
chr8|67618|4369
chr9|35417|2342
chr10|99406|7589
chr11|174043|13191
chr12|37167|2031
chr13|107162|8266
chr14|112913|8140
chr15|11327|1097
chr16|68929|5494
chr17|79303|5740
chr18|32119|1741
chr19|76676|4600
chr20|101277|9768
chr21|140743|11270
chr22|28347|1437
chr23|7959|447
chr24|62277|4564
chr25|6538|333
chr26|60735|6243
chr27|15941|804
chr28|8111|592
chr29|52832|2850
chrX|83725|53674
Total|2204094|218714


This is the first successful trial.

working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/gatk_datasets/

GATK VQSR has failed to run successfully with running the command without the option -maxGaussian 4 


    GenomeAnalysisTK -T VariantRecalibrator \ 
					 -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa \
					 -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz \
					 -resource:HQF,known=false,training=true,truth=true,prior=15.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz \
					 -resource:HQFHET,known=false,training=true,truth=false,prior=10.0 falseset.vcf.gz \
					 -an DP -an MQB -an MQSB -an MQ0F \
					 -mode SNP --maxGaussians 4 \
					-recalFile output.recal \
					-tranchesFile output.tranches \
					-rscriptFile output.plots.R

    INFO  09:48:04,325 HelpFormatter - Executing as kiranmayee.bakshy@assembler2.agil.barc.ba.ars.usda.gov on Linux 4.8.13-100.fc23.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0_131-b11.
    INFO  09:48:04,326 HelpFormatter - Date/Time: 2018/10/11 09:48:04
    INFO  09:48:04,327 HelpFormatter - --------------------------------------------------------------------------------
    INFO  09:48:04,327 HelpFormatter - --------------------------------------------------------------------------------
    INFO  09:48:04,395 GenomeAnalysisEngine - Strictness is SILENT
    INFO  09:48:04,629 GenomeAnalysisEngine - Downsampling Settings: Method: BY_SAMPLE, Target Coverage: 1000
    WARN  09:48:05,455 IndexDictionaryUtils - Track input doesn't have a sequence dictionary built in, skipping dictionary validation
    WARN  09:48:05,456 IndexDictionaryUtils - Track HQF doesn't have a sequence dictionary built in, skipping dictionary validation
    WARN  09:48:05,458 IndexDictionaryUtils - Track HQFHET doesn't have a sequence dictionary built in, skipping dictionary validation
    INFO  09:48:05,655 GenomeAnalysisEngine - Preparing for traversal
    INFO  09:48:05,665 GenomeAnalysisEngine - Done preparing for traversal
    INFO  09:48:05,665 ProgressMeter - [INITIALIZATION COMPLETE; STARTING PROCESSING]
    INFO  09:48:05,666 ProgressMeter - | processed |time |per 1M |   |   total | remaining
    INFO  09:48:05,666 ProgressMeter -Location | sites | elapsed | sites | completed | runtime |   runtime
    INFO  09:48:05,682 TrainingSet - Found HQF track:   Known = false   Training = true Truth = truePrior = Q15.0
    INFO  09:48:05,683 TrainingSet - Found HQFHET track:Known = false   Training = true Truth = false   Prior = Q10.0
    INFO  09:48:35,677 ProgressMeter -   chr1:13948737110830.030.0 s   4.5 m0.5%95.4 m  94.9 m
    INFO  09:49:05,682 ProgressMeter -   chr1:32280220269925.060.0 s   3.7 m1.2%82.4 m  81.4 m
    INFO  09:49:35,691 ProgressMeter -   chr1:49832502409561.090.0 s   3.7 m1.9%80.1 m  78.6 m
    INFO  09:50:05,697 ProgressMeter -   chr1:76943880556657.0   120.0 s   3.6 m2.9%69.2 m  67.2 m
    INFO  09:50:35,700 ProgressMeter -  chr1:100657338702045.0 2.5 m   3.6 m3.8%66.1 m  63.6 m
    INFO  09:51:05,704 ProgressMeter -  chr1:124344115870110.0 3.0 m   3.4 m4.7%64.2 m  61.2 m
    INFO  09:51:35,707 ProgressMeter -  chr1:151740161   1051284.0 3.5 m   3.3 m5.7%61.4 m  57.9 m
    INFO  09:52:05,711 ProgressMeter -   chr2:15677648   1232572.0 4.0 m   3.2 m6.5%61.2 m  57.2 m
    INFO  09:52:35,715 ProgressMeter -   chr2:42676625   1400524.0 4.5 m   3.2 m7.6%59.6 m  55.1 m
    INFO  09:53:05,720 ProgressMeter -   chr2:74324334   1536100.0 5.0 m   3.3 m8.7%57.2 m  52.2 m
    INFO  09:53:35,951 ProgressMeter -   chr2:99196452   1681713.0 5.5 m   3.3 m9.7%56.8 m  51.3 m
    INFO  09:54:05,957 ProgressMeter -  chr2:121778049   1837674.0 6.0 m   3.3 m   10.5%57.0 m  51.0 m
    INFO  09:54:36,009 ProgressMeter -   chr3:11441308   2032465.0 6.5 m   3.2 m   11.5%56.4 m  49.9 m
    INFO  09:55:06,025 ProgressMeter -   chr3:36774467   2209266.0 7.0 m   3.2 m   12.5%56.1 m  49.1 m
    INFO  09:55:36,032 ProgressMeter -   chr3:55288361   2359736.0 7.5 m   3.2 m   13.2%56.9 m  49.4 m
    INFO  09:56:06,043 ProgressMeter -   chr3:82520542   2541489.0 8.0 m   3.2 m   14.2%56.3 m  48.3 m
    INFO  09:56:36,449 ProgressMeter -  chr3:115952058   2722295.0 8.5 m   3.1 m   15.5%55.0 m  46.5 m
    INFO  09:57:06,453 ProgressMeter -   chr4:16829216   2905564.0 9.0 m   3.1 m   16.3%55.2 m  46.2 m
    INFO  09:57:36,457 ProgressMeter -   chr4:44996136   3097617.0 9.5 m   3.1 m   17.4%54.7 m  45.2 m
    INFO  09:58:06,460 ProgressMeter -   chr4:71850471   3245747.010.0 m   3.1 m   18.4%54.5 m  44.5 m
    INFO  09:58:36,462 ProgressMeter -  chr4:101455013   3436332.010.5 m   3.1 m   19.5%53.9 m  43.4 m
    INFO  09:59:06,471 ProgressMeter -  chr4:120139358   3632779.011.0 m   3.0 m   20.2%54.5 m  43.5 m
    INFO  09:59:36,477 ProgressMeter -   chr5:19116645   3791751.011.5 m   3.0 m   20.9%55.0 m  43.5 m
    INFO  10:00:06,481 ProgressMeter -   chr5:46996436   3963376.012.0 m   3.0 m   22.0%54.6 m  42.6 m
    INFO  10:00:36,487 ProgressMeter -   chr5:72317058   4146484.012.5 m   3.0 m   22.9%54.5 m  42.0 m
    INFO  10:01:06,490 ProgressMeter -   chr5:95183817   4321890.013.0 m   3.0 m   23.8%54.7 m  41.7 m
    INFO  10:01:36,493 ProgressMeter -  chr5:117346674   4490931.013.5 m   3.0 m   24.6%54.8 m  41.3 m
    INFO  10:02:06,497 ProgressMeter -   chr6:19822950   4667765.014.0 m   3.0 m   25.5%54.9 m  40.9 m
    INFO  10:02:36,500 ProgressMeter -   chr6:44626088   4847107.014.5 m   3.0 m   26.4%54.8 m  40.3 m
    INFO  10:03:06,503 ProgressMeter -   chr6:68506173   5019150.015.0 m   3.0 m   27.3%54.9 m  39.9 m
    INFO  10:03:36,507 ProgressMeter -   chr6:92273364   5204979.015.5 m   3.0 m   28.2%54.9 m  39.4 m
    INFO  10:04:07,090 ProgressMeter -chr7:4062403   5386561.016.0 m   3.0 m   29.4%54.5 m  38.5 m
    INFO  10:04:37,104 ProgressMeter -   chr7:27217466   5542330.016.5 m   3.0 m   30.3%54.6 m  38.0 m
    INFO  10:05:07,117 ProgressMeter -   chr7:47014012   5678656.017.0 m   3.0 m   31.0%54.9 m  37.8 m
    INFO  10:05:37,120 ProgressMeter -   chr7:76414489   5827812.017.5 m   3.0 m   32.1%54.5 m  37.0 m
    INFO  10:06:07,193 ProgressMeter -  chr7:104787591   5998847.018.0 m   3.0 m   33.2%54.3 m  36.3 m
    INFO  10:06:37,537 ProgressMeter -   chr8:16633037   6190269.018.5 m   3.0 m   34.1%54.3 m  35.8 m
    INFO  10:07:07,540 ProgressMeter -   chr8:45219917   6379223.019.0 m   3.0 m   35.2%54.1 m  35.0 m
    INFO  10:07:44,206 ProgressMeter -   chr8:65828257   6493379.019.6 m   3.0 m   36.0%54.6 m  35.0 m
    INFO  10:08:14,209 ProgressMeter -   chr8:94606598   6670856.020.1 m   3.0 m   37.0%54.4 m  34.2 m
    INFO  10:08:44,213 ProgressMeter -chr9:6338137   6833654.020.6 m   3.0 m   38.0%54.3 m  33.7 m
    INFO  10:09:14,693 ProgressMeter -   chr9:29048331   7014820.021.2 m   3.0 m   38.8%54.5 m  33.3 m
    INFO  10:09:44,698 ProgressMeter -   chr9:55856026   7179522.021.7 m   3.0 m   39.8%54.3 m  32.7 m
    INFO  10:10:14,702 ProgressMeter -   chr9:84065688   7370885.022.2 m   3.0 m   40.9%54.2 m  32.0 m
    INFO  10:10:44,704 ProgressMeter -   chr10:3523196   7557272.022.7 m   3.0 m   41.8%54.1 m  31.5 m
    INFO  10:11:14,706 ProgressMeter -  chr10:35946264   7747185.023.2 m   3.0 m   43.1%53.8 m  30.6 m
    INFO  10:11:44,708 ProgressMeter -  chr10:62462554   7910339.023.7 m   3.0 m   44.1%53.7 m  30.0 m
    INFO  10:12:14,711 ProgressMeter -  chr10:94639044   8055143.024.2 m   3.0 m   45.3%53.3 m  29.2 m
    INFO  10:12:44,715 ProgressMeter -  chr11:10687868   8215455.024.7 m   3.0 m   46.0%53.5 m  28.9 m
    INFO  10:13:14,717 ProgressMeter -  chr11:40206787   8396442.025.2 m   3.0 m   47.1%53.3 m  28.2 m
    INFO  10:13:44,721 ProgressMeter -  chr11:61346961   8568671.025.7 m   3.0 m   47.9%53.5 m  27.9 m
    INFO  10:14:14,724 ProgressMeter -  chr11:85102994   8714006.026.2 m   3.0 m   48.8%53.5 m  27.4 m
    INFO  10:14:44,727 ProgressMeter -   chr12:5847479   8857447.026.7 m   3.0 m   49.9%53.4 m  26.8 m
    INFO  10:15:14,731 ProgressMeter -  chr12:27505982   9033714.027.2 m   3.0 m   50.7%53.5 m  26.4 m
    INFO  10:15:44,734 ProgressMeter -  chr12:55735204   9218520.027.7 m   3.0 m   51.8%53.4 m  25.8 m
    INFO  10:16:14,737 ProgressMeter -  chr12:72715871   9374567.028.2 m   3.0 m   52.4%53.7 m  25.6 m
    INFO  10:16:52,722 ProgressMeter -chr13:822385   9558825.028.8 m   3.0 m   53.1%54.2 m  25.4 m
    INFO  10:17:22,724 ProgressMeter -  chr13:25410051   9732213.029.3 m   3.0 m   54.1%54.2 m  24.9 m
    INFO  10:17:52,726 ProgressMeter -  chr13:48076242   9854453.029.8 m   3.0 m   54.9%54.2 m  24.5 m
    INFO  10:18:22,729 ProgressMeter -  chr13:82098092   1.003561E730.3 m   3.0 m   56.2%53.9 m  23.6 m
    INFO  10:18:52,732 ProgressMeter -  chr14:19070710   1.0166752E730.8 m   3.0 m   57.0%54.0 m  23.2 m
    INFO  10:19:22,735 ProgressMeter -  chr14:47329282   1.0317357E731.3 m   3.0 m   58.0%53.9 m  22.6 m
    INFO  10:19:52,738 ProgressMeter -  chr14:73978186   1.0492405E731.8 m   3.0 m   59.0%53.8 m  22.0 m
    INFO  10:20:22,742 ProgressMeter -  chr15:10423632   1.0688588E732.3 m   3.0 m   59.8%54.0 m  21.7 m
    INFO  10:20:52,745 ProgressMeter -  chr15:34525200   1.0858519E732.8 m   3.0 m   60.7%54.0 m  21.2 m
    INFO  10:21:22,747 ProgressMeter -  chr15:60173263   1.1055702E733.3 m   3.0 m   61.7%53.9 m  20.7 m
    INFO  10:21:52,750 ProgressMeter -   chr16:3027107   1.1261381E733.8 m   3.0 m   62.8%53.8 m  20.0 m
    INFO  10:22:22,752 ProgressMeter -  chr16:21076235   1.142114E734.3 m   3.0 m   63.4%54.0 m  19.8 m
    INFO  10:22:52,755 ProgressMeter -  chr16:60244124   1.1617982E734.8 m   3.0 m   64.9%53.6 m  18.8 m
    INFO  10:23:22,757 ProgressMeter -   chr17:4028604   1.1794728E735.3 m   3.0 m   65.9%53.6 m  18.3 m
    INFO  10:23:52,760 ProgressMeter -  chr17:21682282   1.1922915E735.8 m   3.0 m   66.5%53.8 m  18.0 m
    INFO  10:24:22,762 ProgressMeter -  chr17:41250371   1.211141E736.3 m   3.0 m   67.3%53.9 m  17.7 m
    INFO  10:24:52,765 ProgressMeter -  chr17:73018225   1.2298316E736.8 m   3.0 m   68.5%53.7 m  16.9 m
    INFO  10:25:22,768 ProgressMeter -  chr18:29370357   1.2473773E737.3 m   3.0 m   69.6%53.5 m  16.2 m
    INFO  10:25:52,771 ProgressMeter -  chr18:59733206   1.2644592E737.8 m   3.0 m   70.8%53.4 m  15.6 m
    INFO  10:26:22,773 ProgressMeter -  chr19:24568788   1.2805606E738.3 m   3.0 m   71.9%53.2 m  14.9 m
    INFO  10:26:52,776 ProgressMeter -  chr19:56538677   1.2958822E738.8 m   3.0 m   73.1%53.0 m  14.2 m
    INFO  10:27:22,780 ProgressMeter -  chr20:21182655   1.3124443E739.3 m   3.0 m   74.2%52.9 m  13.6 m
    INFO  10:27:52,783 ProgressMeter -  chr20:45195250   1.3277996E739.8 m   3.0 m   75.1%53.0 m  13.2 m
    INFO  10:28:22,786 ProgressMeter -  chr20:66372665   1.3421048E740.3 m   3.0 m   75.9%53.1 m  12.8 m
    INFO  10:28:52,789 ProgressMeter -  chr21:14594563   1.3554474E740.8 m   3.0 m   76.7%53.2 m  12.4 m
    INFO  10:29:22,791 ProgressMeter -  chr21:34753882   1.3680171E741.3 m   3.0 m   77.4%53.3 m  12.0 m
    INFO  10:29:52,794 ProgressMeter -  chr21:66266195   1.3851222E741.8 m   3.0 m   78.6%53.1 m  11.4 m
    INFO  10:30:23,136 ProgressMeter -chr22:174560   1.3882587E742.3 m   3.0 m   78.8%53.6 m  11.4 m
    INFO  10:30:53,139 ProgressMeter -  chr22:23980757   1.4032327E742.8 m   3.0 m   79.7%53.7 m  10.9 m
    INFO  10:31:23,212 ProgressMeter -  chr22:57294712   1.4222358E743.3 m   3.0 m   81.0%53.4 m  10.2 m
    INFO  10:31:53,221 ProgressMeter -  chr23:23018713   1.4428663E743.8 m   3.0 m   82.0%53.4 m   9.6 m
    INFO  10:32:23,232 ProgressMeter -  chr23:32250655   1.4607336E744.3 m   3.0 m   82.4%53.8 m   9.5 m
    INFO  10:32:53,817 ProgressMeter -   chr24:6421695   1.4785483E744.8 m   3.0 m   83.4%53.7 m   8.9 m
    INFO  10:33:23,826 ProgressMeter -  chr24:24309518   1.4954764E745.3 m   3.0 m   84.0%53.9 m   8.6 m
    INFO  10:33:53,835 ProgressMeter -  chr24:52054836   1.5123684E745.8 m   3.0 m   85.1%53.8 m   8.0 m
    INFO  10:34:23,855 ProgressMeter -  chr25:16322365   1.5307702E746.3 m   3.0 m   86.1%53.8 m   7.5 m
    INFO  10:34:53,862 ProgressMeter -   chr26:2401974   1.5490876E746.8 m   3.0 m   87.2%53.7 m   6.9 m
    INFO  10:35:23,959 ProgressMeter -  chr26:25863575   1.5647882E747.3 m   3.0 m   88.1%53.7 m   6.4 m
    INFO  10:35:53,969 ProgressMeter -  chr26:51370037   1.5816768E747.8 m   3.0 m   89.0%53.7 m   5.9 m
    INFO  10:36:23,983 ProgressMeter -  chr27:23213053   1.6013565E748.3 m   3.0 m   89.9%53.7 m   5.4 m
    INFO  10:36:53,985 ProgressMeter -   chr28:9114638   1.6196229E748.8 m   3.0 m   91.1%53.6 m   4.8 m
    INFO  10:37:23,989 ProgressMeter -  chr28:35988450   1.6383529E749.3 m   3.0 m   92.1%53.5 m   4.2 m
    INFO  10:37:54,174 ProgressMeter -  chr29:11000123   1.6562354E749.8 m   3.0 m   92.9%53.6 m   3.8 m
    INFO  10:38:24,177 ProgressMeter -  chr29:34232511   1.6746015E750.3 m   3.0 m   93.8%53.6 m   3.3 m
    INFO  10:38:54,181 ProgressMeter -   chrX:71933804   1.688864E750.8 m   3.0 m   97.1%52.3 m  90.0 s
    INFO  10:39:24,184 ProgressMeter -   chrX:88112728   1.7013446E751.3 m   3.0 m   97.7%52.5 m  71.0 s
    INFO  10:39:54,187 ProgressMeter -  chrX:148803211   1.7146573E751.8 m   3.0 m  100.0%51.8 m   0.0 s
    INFO  10:40:11,708 VariantDataManager - DP:  mean = 4099.67  standard deviation = 906.40
    INFO  10:40:13,178 VariantDataManager - MQB: mean = 0.90 standard deviation = 0.28
    INFO  10:40:14,895 VariantDataManager - MQSB:mean = 1.00 standard deviation = 0.01
    INFO  10:40:16,561 VariantDataManager - MQ0F:mean = 0.00 standard deviation = 0.00
    INFO  10:40:24,190 ProgressMeter -  chrX:148803211   1.7146573E752.3 m   3.1 m  100.0%52.3 m   0.0 s
    INFO  10:40:31,098 VariantDataManager - Annotations are now ordered by their information content: [DP, MQSB, MQ0F, MQB]
    INFO  10:40:31,641 VariantDataManager - Training with 2200673 variants after standard deviation thresholding.
    INFO  10:40:31,659 GaussianMixtureModel - Initializing model with 100 k-means iterations...
    INFO  10:40:54,193 ProgressMeter -  chrX:148803211   1.7146573E752.8 m   3.1 m  100.0%52.8 m   0.0 s
    INFO  10:41:24,196 ProgressMeter -  chrX:148803211   1.7146573E753.3 m   3.1 m  100.0%53.3 m   0.0 s
    INFO  10:41:54,198 ProgressMeter -  chrX:148803211   1.7146573E753.8 m   3.1 m  100.0%53.8 m   0.0 s
    INFO  10:42:09,760 VariantRecalibratorEngine - Finished iteration 0.
    INFO  10:42:24,201 ProgressMeter -  chrX:148803211   1.7146573E754.3 m   3.2 m  100.0%54.3 m   0.0 s
    INFO  10:43:05,759 ProgressMeter -  chrX:148803211   1.7146573E755.0 m   3.2 m  100.0%55.0 m   0.0 s
    INFO  10:43:16,486 VariantRecalibratorEngine - Finished iteration 5.Current change in mixture coefficients = 0.28551
    INFO  10:43:35,761 ProgressMeter -  chrX:148803211   1.7146573E755.5 m   3.2 m  100.0%55.5 m   0.0 s
    INFO  10:43:52,121 VariantRecalibratorEngine - Finished iteration 10.   Current change in mixture coefficients = 0.02388
    INFO  10:44:05,764 ProgressMeter -  chrX:148803211   1.7146573E756.0 m   3.3 m  100.0%56.0 m   0.0 s
    INFO  10:44:28,072 VariantRecalibratorEngine - Finished iteration 15.   Current change in mixture coefficients = 0.00543
    INFO  10:44:36,195 ProgressMeter -  chrX:148803211   1.7146573E756.5 m   3.3 m  100.0%56.5 m   0.0 s
    INFO  10:44:50,641 VariantRecalibratorEngine - Convergence after 18 iterations!
    INFO  10:44:54,289 VariantRecalibratorEngine - Evaluating full set of 17156523 variants...
    INFO  10:45:06,197 ProgressMeter -  chrX:148803211   1.7146573E757.0 m   3.3 m  100.0%57.0 m   0.0 s
    INFO  10:45:20,833 VariantDataManager - Training with worst 507538 scoring variants --> variants with LOD <= -5.0000.
    INFO  10:45:20,833 GaussianMixtureModel - Initializing model with 100 k-means iterations...
    INFO  10:45:36,199 ProgressMeter -  chrX:148803211   1.7146573E757.5 m   3.4 m  100.0%57.5 m   0.0 s
    INFO  10:45:46,207 VariantRecalibratorEngine - Finished iteration 0.
    INFO  10:45:51,187 VariantRecalibratorEngine - Finished iteration 5.Current change in mixture coefficients = 0.04370
    INFO  10:45:55,875 VariantRecalibratorEngine - Finished iteration 10.   Current change in mixture coefficients = 0.01443
    INFO  10:46:00,490 VariantRecalibratorEngine - Finished iteration 15.   Current change in mixture coefficients = 0.00275
    INFO  10:46:01,403 VariantRecalibratorEngine - Convergence after 16 iterations!
    INFO  10:46:02,508 VariantRecalibratorEngine - Evaluating full set of 17156523 variants...
    INFO  10:46:06,202 ProgressMeter -  chrX:148803211   1.7146573E758.0 m   3.4 m  100.0%58.0 m   0.0 s
    INFO  10:46:36,204 ProgressMeter -  chrX:148803211   1.7146573E758.5 m   3.4 m  100.0%58.5 m   0.0 s
    INFO  10:47:06,207 ProgressMeter -  chrX:148803211   1.7146573E759.0 m   3.4 m  100.0%59.0 m   0.0 s
    INFO  10:47:36,210 ProgressMeter -  chrX:148803211   1.7146573E759.5 m   3.5 m  100.0%59.5 m   0.0 s
    INFO  10:48:06,214 ProgressMeter -  chrX:148803211   1.7146573E760.0 m   3.5 m  100.0%60.0 m   0.0 s
    INFO  10:48:35,889 TrancheManager - Finding 4 tranches for 17156523 variants
    INFO  10:48:36,217 ProgressMeter -  chrX:148803211   1.7146573E760.5 m   3.5 m  100.0%60.5 m   0.0 s
    INFO  10:48:52,005 TrancheManager -   Tranche threshold 100.00 => selection metric threshold 0.000
    INFO  10:48:54,723 TrancheManager -   Found tranche for 100.000: 0.000 threshold starting with variant 0; running score is 0.000
    INFO  10:48:54,723 TrancheManager -   Tranche is Tranche ts=100.00 minVQSLod=-39996.3475 known=(0 @ 0.0000) novel=(17156523 @ 2.1309) truthSites(218714 accessible, 218714 called), name=anonymous]
    INFO  10:48:54,724 TrancheManager -   Tranche threshold 99.90 => selection metric threshold 0.001
    INFO  10:48:57,344 TrancheManager -   Found tranche for 99.900: 0.001 threshold starting with variant 3780286; running score is 0.001
    INFO  10:48:57,345 TrancheManager -   Tranche is Tranche ts=99.90 minVQSLod=-58.4241 known=(0 @ 0.0000) novel=(13376237 @ 2.2001) truthSites(218714 accessible, 218495 called), name=anonymous]
    INFO  10:48:57,345 TrancheManager -   Tranche threshold 99.00 => selection metric threshold 0.010
    INFO  10:48:59,879 TrancheManager -   Found tranche for 99.000: 0.010 threshold starting with variant 5294345; running score is 0.010
    INFO  10:48:59,880 TrancheManager -   Tranche is Tranche ts=99.00 minVQSLod=-4.0027 known=(0 @ 0.0000) novel=(11862178 @ 2.2194) truthSites(218714 accessible, 216526 called), name=anonymous]
    INFO  10:48:59,880 TrancheManager -   Tranche threshold 90.00 => selection metric threshold 0.100
    INFO  10:49:02,476 TrancheManager -   Found tranche for 90.000: 0.100 threshold starting with variant 7211419; running score is 0.100
    INFO  10:49:02,476 TrancheManager -   Tranche is Tranche ts=90.00 minVQSLod=1.8899 known=(0 @ 0.0000) novel=(9945104 @ 2.2243) truthSites(218714 accessible, 196842 called), name=anonymous]
    INFO  10:49:02,483 VariantRecalibrator - Writing out recalibration table...
    INFO  10:49:06,218 ProgressMeter -  chrX:148803211   1.7146573E761.0 m   3.6 m  100.0%61.0 m   0.0 s
    INFO  10:49:36,221 ProgressMeter -  chrX:148803211   1.7146573E761.5 m   3.6 m  100.0%61.5 m   0.0 s
    INFO  10:50:06,225 ProgressMeter -  chrX:148803211   1.7146573E762.0 m   3.6 m  100.0%62.0 m   0.0 s
    INFO  10:50:36,228 ProgressMeter -  chrX:148803211   1.7146573E762.5 m   3.6 m  100.0%62.5 m   0.0 s
    INFO  10:51:03,426 VariantRecalibrator - Writing out visualization Rscript file...
    INFO  10:51:06,132 VariantRecalibrator - Building DP x MQSB plot...
    INFO  10:51:06,143 VariantRecalibratorEngine - Evaluating full set of 3660 variants...
    INFO  10:51:06,230 ProgressMeter -  chrX:148803211   1.7146573E763.0 m   3.7 m  100.0%63.0 m   0.0 s
    INFO  10:51:06,402 VariantRecalibratorEngine - Evaluating full set of 3660 variants...
    INFO  10:51:06,618 VariantRecalibrator - Building DP x MQ0F plot...
    INFO  10:51:06,624 VariantRecalibratorEngine - Evaluating full set of 3660 variants...
    INFO  10:51:06,859 VariantRecalibratorEngine - Evaluating full set of 3660 variants...
    INFO  10:51:07,057 VariantRecalibrator - Building DP x MQB plot...
    INFO  10:51:07,062 VariantRecalibratorEngine - Evaluating full set of 3660 variants...
    INFO  10:51:07,296 VariantRecalibratorEngine - Evaluating full set of 3660 variants...
    INFO  10:51:07,495 VariantRecalibrator - Building MQSB x MQ0F plot...
    INFO  10:51:07,498 VariantRecalibratorEngine - Evaluating full set of 3721 variants...
    INFO  10:51:07,764 VariantRecalibratorEngine - Evaluating full set of 3721 variants...
    INFO  10:51:07,945 VariantRecalibrator - Building MQSB x MQB plot...
    INFO  10:51:07,947 VariantRecalibratorEngine - Evaluating full set of 3721 variants...
    INFO  10:51:08,176 VariantRecalibratorEngine - Evaluating full set of 3721 variants...
    INFO  10:51:08,358 VariantRecalibrator - Building MQ0F x MQB plot...
    INFO  10:51:08,360 VariantRecalibratorEngine - Evaluating full set of 3721 variants...
    INFO  10:51:08,602 VariantRecalibratorEngine - Evaluating full set of 3721 variants...
    INFO  10:51:08,800 VariantRecalibrator - Executing: Rscript /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/gatk_datasets/output.plots.R
    INFO  10:51:20,936 VariantRecalibrator - Executing: Rscript (resource)org/broadinstitute/gatk/tools/walkers/variantrecalibration/plot_Tranches.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/gatk_datasets/output.tranches 2.15
    INFO  10:51:21,887 ProgressMeter -done   1.7156496E763.3 m   3.7 m  100.0%63.3 m   0.0 s
    INFO  10:51:21,888 ProgressMeter - Total runtime 3796.22 secs, 63.27 min, 1.05 hours
    ------------------------------------------------------------------------------------------
    Done. There were 3 WARN messages, the first 3 are repeated below.
    WARN  09:48:05,455 IndexDictionaryUtils - Track input doesn't have a sequence dictionary built in, skipping dictionary validation
    WARN  09:48:05,456 IndexDictionaryUtils - Track HQF doesn't have a sequence dictionary built in, skipping dictionary validation
    WARN  09:48:05,458 IndexDictionaryUtils - Track HQFHET doesn't have a sequence dictionary built in, skipping dictionary validation
    ------------------------------------------------------------------------------------------
    

## Apply recalibration ##


    GenomeAnalysisTK -T ApplyRecalibration \
					 -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa \
					 -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz \
					-mode SNP \
					--ts_filter_level 99.0 \
					-recalFile output.recal \
					-tranchesFile output.tranches \
					-o recalibrated.vcf


    [kiranmayee.bakshy@assembler2 trial_1]$ gunzip -c recalibrated.vcf.gz | grep -v '#' | grep -c 'PASS'
    11862181
    [kiranmayee.bakshy@assembler2 trial_1]$ gunzip -c recalibrated.vcf.gz | grep -v '#' | grep -c 'VQSRTrancheSNP99.00to99.90'
    1514050
    [kiranmayee.bakshy@assembler2 trial_1]$ gunzip -c recalibrated.vcf.gz | grep -v '#' | grep -c 'VQSRTrancheSNP99.90to100.00'
    3780292



## Observations

There are many true positives and very less false positives. 

I think I should change the falseset to only high quality variants QUAL >= 999, 5 > MAF < 50, HWE p-valuE > 1E-05 which are homozygous variants in the haplotype groups and any state in the non-target animals.

## Trial 2 ##

Trial - 2 is similar to the above but with few more annotations added to the model: MQ and AN

This trial looks even worse and I need to change the training false dataset and include BQB annotation to the model.

So far i.e. in trial 1 DP and MQSB seems to be important parameters for training a VQSR model and the plot looks better than the others with a clear distinction to the good and bad variants.

I also should give more tranches and try to see if the model does better.


## Trial 3 ##


    GenomeAnalysisTK -T VariantRecalibrator \
					 -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa \
					 -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz \
					 -resource:HQFiltered,known=false,training=true,truth=true,prior=12.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz \
					-resource:HQFalse,known=false,training=true,truth=false,prior=10.0 training_falseset.vcf.gz \
					-an DP -an MQB -an MQSB -an MQ0F -an BQB -an RPB \ 
					-mode SNP --maxGaussians 8 \
					-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -tranche 50.0 \
					-recalFile trial_3/output.recal \
					-tranchesFile trial_3/output.tranches \
					-rscriptFile trial_3/output.plots.R`

 
This trial looks even more worse...so go back to the trial 1 settings and change the parameters instead of changing the falseset. 

## Trial 4 ##


    GenomeAnalysisTK -T VariantRecalibrator \
					 -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa \ 
					 -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz \
					 -resource:HQF,known=false,training=true,truth=true,prior=15.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz \
					 -resource:HQFHET,known=false,training=true,truth=false,prior=10.0 falseset.vcf.gz \
					 -an DP -an MQSB -an BQB -an RPB 
					 -mode SNP --maxGaussians 4 \
					 -recalFile trial_4/output.recal 
					 -tranchesFile trial_4/output.tranches 
					 -rscriptFile trial_4/output.plots.R

  

 	GenomeAnalysisTK -T ApplyRecalibration \
					 -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa \
					 -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz \
					-mode SNP \
					--ts_filter_level 90.0 \
					-recalFile output.recal \
					-tranchesFile output.tranches \
					-o recalibrated.vcf


Trial 4 seems to be better than all other models as there is a formation of a full green spot that indicates that all the true variants clustered at that spot.

This green cluster is missing all the other models. 

One more last trial: Just add one more annotation MQ0F which may be important...

## Trial 5 ##

 	GenomeAnalysisTK -T VariantRecalibrator \
					 -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa \ 
					 -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz \
					 -resource:HQF,known=false,training=true,truth=true,prior=15.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz \
					 -resource:HQFHET,known=false,training=true,truth=false,prior=10.0 falseset.vcf.gz \
					 -an DP -an MQSB -an BQB -an RPB -an MQ0F
					 -mode SNP --maxGaussians 4 \
					 -recalFile trial_5/output.recal 
					 -tranchesFile trial_5/output.tranches 
					 -rscriptFile trial_5/output.plots.R

OK, this is even worse than trial 4. 

## Trial 6 ##

    GenomeAnalysisTK \ 
    -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa \ 
    -T VariantRecalibrator \ 
    -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz \ 
    -resource:HQF,known=false,training=true,truth=true,prior=15.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz \ 
    -resource:HQFHET,known=false,training=true,truth=false,prior=10.0 falseset/falseset.vcf.gz \ 
    -an DP \
    -an MQSB \
    -an BQB \
    -an RPB \
    -an MQB \
    -an ICB \ 
    -mode SNP \ 
    --maxGaussians 4 \ 
    -recalFile Run1/trial_6/output.recal \ 
    -tranchesFile Run1/trial_6/output.tranches \ 
    -rscriptFile Run1/trial_6/output.plots.R \

This also is no better than trial 4.

Ok, now lets proceed with the trial 4 model and plot some density profiles of all, passed and failed sites.

The attributes that can be plotted are QUAL, DP, VQSLOD, MQSB, RPB, BQB, MQB, ICB, MAF, HWE

All these attributes should be extracted from the recalibrated.vcf.gz file

bcftools to extract the fields and Rscript to plot 

    # First fill in the tags MAF and HWE p-values in the recalibrated vcf file
	bcftools plugin fill-tags Run1/trial_4/recalibrated.vcf.gz -O z -o Run1/trial_4/recalibrated_maf.vcf.gz -- -t MAF,HWE
    
	# Tabix index for future
	tabix -p vcf recalibrated_maf.vcf.gz 
	
	# Now extract all the required fields from the vcf 
    bcftools query -f '%ID\t%QUAL\t%FILTER\t%INFO/DP\t%INFO/MQSB\t%INFO/RPB\t%INFO/MQB\t%INFO/MQ\t%INFO/ICB\t%INFO/BQB\t%INFO/VQSLOD\t%INFO/MAF\t%INFO/HWE\n' -o recalibrated.att recalibrated_maf.vcf.gz 
    
 
Now plot the attributes using the following Rscript:

    #!/usr/bin/env Rscript
    #Date:10/22/2018
    #Author:Kiranmayee Bakshy
    
    # A program to plot VCF attributes after recalibration by GATK
    library(ggplot2)
    library(ggridges)
    library(data.table)
    
    df<-fread("recalibrated.att", sep="\t", header=F, stringsAsFactors=F)
    colnames(df)<-c("id","qual","status", "dp", "mqsb", "rpb", "mqb", "mq", "icb", "bqb", "vqslod", "maf", "hwe")
    
	# There are few NA's in the following attributes which are probably missing because of the nature of the attribute.
    cols<-c("mqsb", "rpb", "mqb", "icb", "bqb")
    df[cols]<-sapply(df[cols], as.numeric)
    df<-as.data.frame(df)
    summary(df)
    
    pdf(file="recalibrated.att.Rplots.pdf", onefile = T)
    ggplot() + geom_density(data = df, aes(qual, fill = "Alldata")) + geom_density(data=subset(df, status=="PASS"), aes(qual, fill = "Recalibrated")) + labs(x = "QUAL", y = "Density") + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data = df, aes(mq, fill = "Alldata")) + geom_density(data=subset(df, status=="PASS"), aes(mq, fill = "Recalibrated")) + labs(x = "MQ", y = "Density") + coord_cartesian(xlim = c(45, 60)) + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data = df, aes(dp, fill = "Alldata")) + geom_density(data=subset(df, status=="PASS"), aes(dp, fill = "Recalibrated")) + labs(x = "DP", y = "Density") + coord_cartesian(xlim = c(0, 15000)) + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data = df, aes(vqslod, fill = "Alldata")) + geom_density(data=subset(df, status=="PASS"), aes(vqslod, fill = "Recalibrated")) + labs(x = "VQSLOD", y = "Density") + coord_cartesian(xlim = c(-250, 7)) + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data = df, aes(maf, fill = "Alldata")) + geom_density(data=subset(df, status=="PASS"), aes(maf, fill = "Recalibrated")) + labs(x = "MAF", y = "Density") + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data = df, aes(hwe, fill = "Alldata")) + geom_density(data=subset(df, status=="PASS"), aes(hwe, fill = "Recalibrated")) + labs(x = "HWE (p-value)", y = "Density") + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data= subset(df, !is.na(mqsb)), aes(mqsb, fill = "Alldata")) + geom_density(data=subset(df, mqsb!="NA" & status=="PASS"), aes(mqsb, fill = "Recalibrated")) + labs(x = "MQSB", y = "Density") + coord_cartesian(xlim = c(0.95, 1.01)) + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data= subset(df, !is.na(mqb)), aes(mqb, fill = "Alldata")) + geom_density(data=subset(df, mqb!="NA" & status=="PASS"), aes(mqb, fill = "Recalibrated")) + labs(x = "MQB", y = "Density") +scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black")) 
    ggplot() + geom_density(data= subset(df, !is.na(bqb)), aes(bqb, fill = "Alldata")) + geom_density(data=subset(df, bqb!="NA" & status=="PASS"), aes(bqb, fill = "Recalibrated")) + labs(x = "BQB", y = "Density") + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data= subset(df, !is.na(rpb)), aes(rpb, fill = "Alldata")) + geom_density(data=subset(df, rpb!="NA" & status=="PASS"), aes(rpb, fill = "Recalibrated")) + labs(x = "RPB", y = "Density") + scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    ggplot() + geom_density(data= subset(df, !is.na(icb)), aes(icb, fill = "Alldata")) + geom_density(data=subset(df, icb!="NA" & status=="PASS"), aes(icb, fill = "Recalibrated")) + labs(x = "ICB", y = "Density")+ scale_fill_manual("Dataset", values=c("Alldata" = alpha("green", 0.3), "Recalibrated" = alpha("red",0.3))) + scale_color_manual(values=c("Alldata" = "black", "Recalibrated" = "black"))
    dev.off()
    
    
    pdf(file="recalibrated.att.Rplots_ridges.pdf", onefile = T)
    ggplot(df, aes(qual, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "QUAL", y = "Tranche") + theme(legend.position="none")
    ggplot(df, aes(dp, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "DP", y = "Tranche") +   xlim(0, 12500) + theme(legend.position="none")
    ggplot(df, aes(mq, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "MQ", y = "Tranche") +   xlim(45, 61) + theme(legend.position="none")
    ggplot(df, aes(vqslod, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "VQSLOD", y = "Tranche") +   xlim(-1200, 60) + theme(legend.position="none")
    ggplot(df, aes(maf, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "MAF", y = "Tranche") +  theme(legend.position="none")
    ggplot(df, aes(hwe, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "HWE (p-value)", y = "Tranche") + theme(legend.position="none")
    ggplot(data=subset(df, !is.na(mqsb)), aes(mqsb, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "MQSB", y = "Tranche") + theme(legend.position="none")
    ggplot(data=subset(df, !is.na(mqb)), aes(mqb, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "MQB", y = "Tranche") + theme(legend.position="none")
    ggplot(data=subset(df, !is.na(rpb)), aes(rpb, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "RPB", y = "Tranche") + theme(legend.position="none")
    ggplot(data=subset(df, !is.na(icb)), aes(icb, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "ICB", y = "Tranche") + theme(legend.position="none")
    ggplot(data=subset(df, !is.na(bqb)), aes(bqb, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "BQB", y = "Tranche") + theme(legend.position="none")
    dev.off()


All the plots are density distributions of the vcf properties for all the variants recalibrated using the trial_4 GATK model.

Tranche 99.90 to 100.00: 10% are false positives
Tranche 99.00 to 99.90: 1.0 % are false-positives
Tranche 90.00 to 99.00: 0.1 % are false-positives 
PASS: Tranche 90.0 -> no false positives allowed 

 
Looks like this model has efficiently filtered away the false-positives as observed from some of the attributes for ex: HWE, DP, MQ, QUAL, MQSB and RPB