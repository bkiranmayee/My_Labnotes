VQSR trials results analysis
============================

There are currently 2 models under study:

- model1: Training sets used => bovineHD, 1000bulls, dbSNP
- model2: 					 => model1 + high filtered ARS custom dataset				

The input data set is calibrated using both the models and I need to analyze the results...

**The some stats of model building are summarized in the tables below:**

**MODEL**|**Convergence (iterations)**|**Positive training set**|**Negative training set**
:-----:|:-----:|:-----:|:-----:
1|Good: 82 |4,134,582|144,946
 |Bad: 8| | 
2|Good: 83|4,370,029|167,442
 |Bad: 3| | 

**MODEL**|**DP**| |**MQSB**| |**BQB**| |**RPB**| |**VQSLOD**| | 
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
 |Mean|SD|Mean|SD|Mean|SD|Mean|SD|Mean|Min|Max
1|4302.6|1024.6|0.84|0.32|0.5|0.37|0.61|0.34|3.955|-545.4|9.05
2|4295.9|1013.5|0.85|0.31|0.5|0.37|0.61|0.34|4.016|-3.289|165.9


**Summary of 99.0 tranche**  

**MODEL**|**TS:99.0 minVQSLod**|**known**|**novel**|**Accesible Truth (584625)**|**Positive training sites that passed**|**Positive training sites that failed**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
1|3.28|4,932,596|7,178,883|578,778|3,408,874|728,721
2|3.12|5,075,681|7,602,879|578,778|3,769,288|603,850

 

**Objectives:**


- Find the overlap between the sites of both models i.e. passed and failed sites
- Next plot the density distribution of all the attributes including MAF


So, first step is separate the passed and failed sites and determine the amount of overlap...

This can be done using the bcftools...

 	# Separating the passed and failed sites
    bcftools filter -i 'FILTER="PASS"' recalibrated.3.2.vcf -O z -o recalibrated_pass.3.2.vcf.gz 
    bcftools filter -e 'FILTER="PASS"' recalibrated.3.2.vcf -O z -o recalibrated_fail.3.2.vcf.gz 
    bcftools filter -i 'FILTER="PASS"' recalibrated.1.2.vcf -O z -o recalibrated_pass.1.2.vcf.gz 
    bcftools filter -e 'FILTER="PASS"' recalibrated.1.2.vcf -O z -o recalibrated_fail.1.2.vcf.gz
	bgzip  run1/recalibrated.1.2.vcf
    tabix -p vcf run1/recalibrated_pass.1.2.vcf.gz
    tabix -p vcf run3/recalibrated_pass.3.2.vcf.gz
    tabix -p vcf run1/recalibrated_fail.1.2.vcf.gz
    tabix -p vcf run3/recalibrated_fail.3.2.vcf.gz
    mkdir analysis
 

  	# Determine the overlap between sites

	bcftools isec -p analysis/pass run1/recalibrated_pass.1.2.vcf.gz run3/recalibrated_pass.3.2.vcf.gz
	bcftools isec -p analysis/fail run1/recalibrated_fail.1.2.vcf.gz run3/recalibrated_fail.3.2.vcf.gz
    grep -v -E '^#' 0000.vcf | cut -f 1 | sort | uniq -c > uniq_1.2_per_chrom
    grep -v -E '^#' 0001.vcf | cut -f 1 | sort | uniq -c > uniq_3.2_per_chrom


 |**PASSED**| | |**FAILED**| | 
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
**CHROM**|**MODEL1**|**MODEL2**|**Common to both**|**MODEL1**|**MODEL2**|**Common**
1|49779|81163|751029|81163|49779|322997
2|33277|60368|561153|60368|33277|221072
3|31932|59582|553905|59582|31932|224130
4|34135|61115|596093|61115|34135|262341
5|36090|63628|572366|63628|36090|254974
6|38939|54253|541846|54253|38939|232377
7|28617|52759|465637|52759|28617|198826
8|26545|51342|489710|51342|26545|191929
9|32565|51352|511470|51352|32565|206186
10|19571|42443|424196|42443|19571|147633
11|24039|51925|466774|51925|24039|171174
12|40996|45533|458285|45533|40996|284266
13|11903|36873|340171|36873|11903|110618
14|17924|42125|385389|42125|17924|138965
15|31265|44581|425603|44581|31265|188148
16|20768|38682|365112|38682|20768|138843
17|23465|45822|378622|45822|23465|142064
18|11372|28939|265071|28939|11372|93024
19|6586|25031|212789|25031|6586|58595
20|19178|37166|333563|37166|19178|131746
21|11639|31005|256051|31005|11639|95771
22|9545|26961|272464|26961|9545|80741
23|21696|32994|330856|32994|21696|156345
24|14502|30732|307014|30732|14502|98267
25|6906|22779|212245|22779|6906|54434
26|13352|27420|254007|27420|13352|88880
27|10404|23294|237884|23294|10404|84495
28|9631|21099|213673|21099|9631|67701
29|16807|29640|275129|29640|16807|104488
**total**|**653428**|**1220606**|**11458107**|**1220606**|**653428**|**4551030**


It is interesting to note that the exact same passed sites which were unique to MODEL1 have been failed by MODEL2 and vice versa. This is perhaps because of the different VQSLOD scores given to each of the site by each model.

This may indicate that the Holstein breed specific sites that are picked up by MODEL2 have been ignored by MODEL1. I need to analyze further.

    library(VennDiagram)
    Loading required package: grid
    Loading required package: futile.logger
    > x<-list(Model1=c(1:12111535), Model2=c(653429:13332141))
    > venn.diagram(x,"gatk_models_passed.png", lty= "blank", fill = c("gold1", "firebrick")
    + )
    [1] 1
    > x<-list(Model1=c(1:5771636), Model2=c(5771637:6425065))
    > venn.diagram(x,"gatk_models_failed.png", lty= "blank", fill = c("gold1", "firebrick"))
    [1] 1
    > x<-list(Model1=c(1:5771636), Model2=c(653429:6425065))
    > venn.diagram(x,"gatk_models_failed.png", lty= "blank", fill = c("gold1", "firebrick"))
    [1] 1
    > x<-list(Model1=c(1:5771636), Model2=c(1220606:6425064))
    > venn.diagram(x,"gatk_models_failed.png", lty= "blank", fill = c("gold1", "firebrick"))
    [1] 1

 


![Passed_sites](https://i.imgur.com/5pJFcWc.png)

![failed_sites](https://i.imgur.com/Ir9LcUR.png)
	
### Now plot density distributions of the attributes ###

The attributes that can be plotted are QUAL, AF, AC, AN, DP, VQSLOD, MQSB, RPB, BQB, MQB, ICB, MAF, HWE

All these attributes should be extracted from the recalibrated.vcf.gz files

bcftools to extract the fields and Rscript to plot 

    bcftools plugin fill-tags run1/recalibrated.1.2.vcf.gz -O z -o run1/recalibrated_tagged.1.2.vcf.gz -- -t AF,MAF,HWE
    
	tabix -p vcf run1/recalibrated_tagged.1.2.vcf.gz 

    bcftools query -f '%ID\t%QUAL\t%FILTER\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/DP\t%INFO/MQSB\t%INFO/RPB\t%INFO/MQB\t%INFO/MQ\t%INFO/ICB\t%INFO/BQB\t%INFO/VQSLOD\t%INFO/MAF\t%INFO/HWE\n' -o run1/recalibrated_1.att run1/recalibrated_tagged.1.2.vcf.gz 
    
	bcftools plugin fill-tags run3/recalibrated.3.2.vcf.gz -O z -o run3/recalibrated_tagged.3.2.vcf.gz -- -t AF,MAF,HWE
    
	tabix -p vcf run3/recalibrated_tagged.3.2.vcf.gz 

    bcftools query -f '%ID\t%QUAL\t%FILTER\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/DP\t%INFO/MQSB\t%INFO/RPB\t%INFO/MQB\t%INFO/MQ\t%INFO/ICB\t%INFO/BQB\t%INFO/VQSLOD\t%INFO/MAF\t%INFO/HWE\n' -o run3/recalibrated_3.att run3/recalibrated_tagged.3.2.vcf.gz 
 
Now plot the attributes using the following Rscript:

    #!/usr/bin/env Rscript
    #Date:10/22/2018
    #Author:Kiranmayee Bakshy
    
    # A program to plot VCF attributes after recalibration by GATK
    library(ggplot2)
    library(ggridges)
    library(data.table)
    
	args = commandArgs(trailingOnly=TRUE)

	# test if there is at least one argument: if not, return an error
	if (length(args)==0) {
	  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
	} else if (length(args)==1) {
	  # default output file
	  args[2] = "out.pdf"
	}

    df<-fread(args[1], sep="\t", header=F, stringsAsFactors=F)
    colnames(df)<-c("id","qual","status", "ac", "af", "an", "dp", "mqsb", "rpb", "mqb", "mq", "icb", "bqb", "vqslod", "maf", "hwe")
    
    cols<-c("mqsb", "rpb", "mqb", "icb", "bqb")
    df[cols]<-sapply(df[cols], as.numeric)
    
    summary(df)
    
    pdf(file=args[2], onefile = T)
    ggplot(df, aes(qual, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "QUAL", y = "Tranche") + theme(legend.position="none")
	ggplot(df, aes(ac, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "Alt Allele count", y = "Tranche") + theme(legend.position="none")
	ggplot(df, aes(af, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "Alt Allele Freq", y = "Tranche") + theme(legend.position="none")
	ggplot(df, aes(an, status)) + geom_density_ridges(aes(fill = as.factor(status))) + labs(x = "Total Allele NUmber", y = "Tranche") + theme(legend.position="none")
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


All the plots are density distributions of the vcf properties for all the variants re calibrated using both the models separately.


Now plot the attributes of the passed and failed variants separately...

### Summary of the attributes ###

    **MODEL 1:**
    
     id 			qual			 status			ac
     Length:17883171	Min.   :  3.01   Length:17883171	Min.   :  1.00  
     Class :character   1st Qu.:422.00   Class :character   1st Qu.:  5.00  
     Mode  :character   Median :999.00   Mode  :character   Median : 40.00  
    Mean   :755.91  Mean   : 79.42  
    3rd Qu.:999.00  3rd Qu.:127.00  
    Max.   :999.00  Max.   :344.00  
    
     af 				an  			dp 			mqsb   
     Min.   :0.002907   Min.   :  2.0   Min.   :1   Min.   :0.0000  
     1st Qu.:0.014535   1st Qu.:344.0   1st Qu.: 3679   1st Qu.:0.8226  
     Median :0.116279   Median :344.0   Median : 4236   Median :0.9974  
     Mean   :0.232313   Mean   :342.2   Mean   : 4348   Mean   :0.7982  
     3rd Qu.:0.372093   3rd Qu.:344.0   3rd Qu.: 4828   3rd Qu.:0.9996  
     Max.   :1.000000   Max.   :344.0   Max.   :66018   Max.   :1.0128  
    NA's   :1637

     rpb 			mqb  			mq 				icb
     Min.   :0.00	Min.   :0.00	Min.   : 0.00   Min.   :0.00
     1st Qu.:0.171st Qu.:0.001st Qu.:59.00   1st Qu.:0.00
     Median :0.62Median :1.00Median :59.00   Median :0.67
     Mean   :0.55Mean   :0.65Mean   :56.99   Mean   :0.52
     3rd Qu.:0.903rd Qu.:1.003rd Qu.:59.00   3rd Qu.:0.95
     Max.   :1.01Max.   :1.01Max.   :60.00   Max.   :1.00
     NA's   :47150   NA's   :47150   NA's   :145325 
 
     bqb			vqslod  			 maf   				hwe
     Min.   :0.00**	Min.   :-545.400**   Min.   :0.00000   Min.   :0.0000  
     1st Qu.:0.091st Qu.:   2.840   1st Qu.:0.01453   1st Qu.:0.3406  
     Median :0.50Median :   4.250   Median :0.09593   Median :0.8620  
     Mean   :0.49**Mean   :   3.955**   Mean   :0.15186   Mean   :0.6783  
     3rd Qu.:0.863rd Qu.:   5.440   3rd Qu.:0.26608   3rd Qu.:1.0000  
     Max.   :1.01**Max.   :   9.050**   Max.   :0.50000   Max.   :1.0000  
     NA's   :47150



    **MODEL2:**
    
      
     bqb        	vqslod 				maf			   	 hwe
     Min.   :0.00**	Min.   : -3.289**   Min.   :0.00000   Min.   :0.0000  
     1st Qu.:0.091st Qu.:  2.880   1st Qu.:0.01453   1st Qu.:0.3406  
     Median :0.50Median :  4.220   Median :0.09593   Median :0.8620  
     Mean   :0.49**Mean   :  4.016**   Mean   :0.15186   Mean   :0.6783  
     3rd Qu.:0.863rd Qu.:  5.410   3rd Qu.:0.26608   3rd Qu.:1.0000  
     Max.   :1.01**Max.   :165.970**   Max.   :0.50000   Max.   :1.0000  
     NA's   :47150                                               

All the attributes summary of the models are the same except vqslod which is what we expect...

It is interesting to note that the min, max vqslod of model2 are much higher than that of model1 although the mean is nearly the same...which is good.

I guess this is an improvement to our model.

 


**How many positive training sites are present among the passed and failed sites of the 2 models?**

	#MODEL1:
    gunzip -c recalibrated_pass.1.2.vcf.gz | grep -c 'POSITIVE'
    3408874
    gunzip -c recalibrated_fail.1.2.vcf.gz | grep -c 'POSITIVE'
    728721
	Total = 4,137,595
	#MODEL2:
	gunzip -c recalibrated_pass.3.2.vcf.gz | grep -c 'POSITIVE'
	3769288
	gunzip -c recalibrated_fail.3.2.vcf.gz | grep -c 'POSITIVE'
	603850
	Total = 4,373,138

Difference between total positive training sites between the models is **235,543** and these must belong to the highly filtered Holstein specific data that was added to the MODEL2.


Make boxplots of the VQSLOD scores and QUAL scores of all, passed and failed sites separately to compare both the models...

#### Summary of some statistics of the attributes of passed and failed sites to compare between models: ####

    df3<-readRDS("merged_attributes.rds")
    p1<-do.call(cbind, lapply(df3[which(df3$status1=="PASS" & df3$model=="model1"), c(2,4:16)],summary))
    p2<-do.call(cbind, lapply(df3[which(df3$status1=="PASS" & df3$model=="model2"), c(2,4:16)],summary))
    f1<-do.call(cbind, lapply(df3[which(df3$status1=="Fail" & df3$model=="model1"), c(2,4:16)],summary))
    f2<-do.call(cbind, lapply(df3[which(df3$status1=="Fail" & df3$model=="model2"), c(2,4:16)],summary))
    write.table(t(f1), "f1.txt", sep="\t", quote=F, row.names=T)
    write.table(t(f2), "f2.txt", sep="\t", quote=F, row.names=T)
    write.table(t(p1), "p1.txt", sep="\t", quote=F, row.names=T)
    write.table(t(p2), "p2.txt", sep="\t", quote=F, row.names=T)



**PASSED SITES**|**Min.**| |**1st Qu.**| |**Median**| |**Mean**| |**3rd Qu.**| |**Max.**| |**NA's**| 
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
 |**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**
**qual**|3.01|3.01|999|999|999|999|829|819.6|999|999|999|999| | 
**ac**|1|1|10|9|49|47|86.32|84.41|137|134|344|344| | 
**af**|0.002907|0.002907|0.02907|0.02616|0.1433|0.1366|0.2513|0.2456|0.3988|0.3895|1|1| | 
**an**|2|4|344|344|344|344|343.5|343.7|344|344|344|344| | 
**dp**|1508|1294|3640|3743|4122|4227|4177|4304|4672|4749|8331|66020| | 
**mqsb**|0|0|0.9962|0.9943|0.9995|0.9994|0.8687|0.8813|0.9999|0.9998|1|1|1| 
**rpb**|0|0|0.2314|0.2555|0.6608|0.6676|0.5728|0.5805|0.9182|0.9181|1|1|38408|38965
**mqb**|0|0|0.964|0.9207|0.9995|0.9994|0.7961|0.784|1|1|1.013|1.013|38408|38965
**mq**|0|0|59|59|59|59|57.83|57.96|59|59|60|60| | 
**icb**|4.47E-38|4.47E-38|0.002947|0.002307|0.7782|0.7663|0.5844|0.5755|0.9599|0.9578|1|1|107776|110611
**bqb**|0|0|0.09539|0.09978|0.4983|0.4937|0.4889|0.4875|0.8634|0.8586|1|1|38408|38965
**vqslod**|3.28|3.13|4.22|4.07|5.12|4.86|5.086|4.968|5.73|5.69|**9.05**|**166**| | 
**maf**|0|0|0.02339|0.02047|0.1192|0.1134|0.1641|0.1611|0.282|0.2762|0.5|0.5| | 
**hwe**|0|0|0.3294|0.3379|0.7588|0.7698|0.6588|0.6651|1|1|1|1| | 




**FAILED SITES**|**Min.**||**1st Qu.**||**Median**| |**Mean**| |**3rd Qu.**| |**Max.**| |**NA's**| 
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
 |**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**|**MODEL1**|**MODEL2**
**qual**|3.01|3.01|54|56|999|999|602.6|600.8|999|999|999|999| | 
**ac**|1|1|3|3|18|19|64.94|67.26|103|109|344|344| | 
**af**|0.002907|0.002907|0.008721|0.008721|0.05294|0.05814|0.1926|0.1999|0.3052|0.3227|1|1| | 
**an**|2|2|344|344|344|344|339.3|338.5|344|344|344|344| | 
**dp**|1|1|3848|3464|4518|4276|4706|4454|5179|5064|66020|23500| | 
**mqsb**|0|0|0.2836|0.1375|0.8493|0.7732|0.6503|0.5957|0.9756|0.9592|1.013|1.013|1636|1637
**rpb**|0|0|0.09363|0.05189|0.5064|0.447|0.4911|0.4636|0.8652|0.8542|1.013|1.013|8742|8185
**mqb**|0|0|6.60E-14|2.02E-16|0.009106|0.003233|0.3582|0.3399|0.9379|0.9135|1.013|1.013|8742|8185
**mq**|0|0|56|55|59|59|55.22|54.62|59|59|60|60| | 
**icb**|4.47E-38|4.47E-38|0.0001368|0.0001368|0.003761|0.008208|0.3839|0.3836|0.8912|0.8907|1|1|37549|34714
**bqb**|0|0|0.08898|0.07667|0.4946|0.5063|0.4864|0.4895|0.8632|0.8742|1.013|1.013|8742|8185
**vqslod**|**-545.4**|**-3.289**|1.52|1.32|2.17|2|1.581|1.695|2.8|2.62|3.28|3.13| | 
**maf**|0|0|0.006098|0.008721|0.0436|0.04942|0.1261|0.1293|0.2238|0.2297|0.5|0.5| | 
**hwe**|0|0|0.3878|0.3666|1|1|0.7192|0.7104|1|1|1|1| | 






    #!/usr/bin/env Rscript
    #Date:11/16/2018
    #Author:Kiranmayee Bakshy
    
    # A program to plot comparisons of all the VCF attributes after recalibration by GATK VQSR models
    library(ggplot2)
    library(data.table)
    
    args = commandArgs(trailingOnly=TRUE)
    
    # test if there are at least two argument: if not, return an error
    if (length(args)==0) {
      stop("At least two arguments must be supplied (vcf attributes of 2 files to be compared).\n", call.=FALSE)
    } else if (length(args)==2) {
      # default output file
      args[3] = "out_passed.pdf"
      args[4] = "out_failed.pdf"
    }
    
    
    df_1<-fread(args[1], sep="\t", header=F, stringsAsFactors=F)
    df_2<-fread(args[2], sep="\t", header=F, stringsAsFactors=F)
    
    message("loaded attributes")
    
    colnames(df_1)<-c("id","qual","status", "ac", "af", "an", "dp", "mqsb", "rpb", "mqb", "mq", "icb", "bqb", "vqslod", "maf", "hwe")
    colnames(df_2)<-c("id","qual","status", "ac", "af", "an", "dp", "mqsb", "rpb", "mqb", "mq", "icb", "bqb", "vqslod", "maf", "hwe")
    
    df1<-as.data.frame(df_1)
    df2<-as.data.frame(df_2)
    rm(df_1,df_2)
      
    cols<-c("mqsb", "rpb", "mqb", "icb", "bqb")
    df1[cols]<-sapply(df1[cols], as.numeric)
    df2[cols]<-sapply(df2[cols], as.numeric)
    df1$model<-"model1"
    df2$model<-"model2"
    
    df1$status1<-gsub("VQSRTranche.*", "Fail", df1$status)
    df2$status1<-gsub("VQSRTranche.*", "Fail", df2$status)
    df3<-rbind(df1,df2)
    
    saveRDS(df3, "merged_attribures.rds")
    rm(df1,df2)
    
    message("summarizing the subsets...")
    print("Model1 failed sites summary:")
    summary(subset(df3, status1=="Fail" & model=="model1"))
    print("Model1 passed sites summary:") 
    summary(subset(df3, status1=="PASS" & model=="model1"))
    print("Model2 failed sites summary:") 
    summary(subset(df3, status1=="Fail" & model=="model2"))
    print("Model2 passed sites summary:") 
    summary(subset(df3, status1=="PASS" & model=="model2"))

    message("plotting attributes of passed sites...")
    pdf(file=args[2], onefile = T)
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(qual), fill = as.factor(model), alpha=0.3) + labs(x = "QUAL", y = "Density")
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(af), fill = as.factor(model), alpha=0.3) + labs(x = "Alt Allele Freq", y = "Density")
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(maf), fill = as.factor(model), alpha=0.3) + labs(x = "Minor allele freq", y = "Density")
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(ac), fill = as.factor(model), alpha=0.3) + labs(x = "Alt allele count", y = "Density")
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(an), fill = as.factor(model), alpha=0.3) + labs(x = "Total allele number", y = "Density")
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(mq), fill = as.factor(model), alpha=0.3) + labs(x = "MQ", y = "Density") + coord_cartesian(xlim = c(45, 60))
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(dp), fill = as.factor(model), alpha=0.3) + labs(x = "DP", y = "Density") + coord_cartesian(xlim = c(0, 12500))
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(vqslod), fill = as.factor(model), alpha=0.3) + labs(x = "VQSLOD", y = "Density")
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(maf), fill = as.factor(model), alpha=0.3) + labs(x = "MAF", y = "Density") 
    ggplot() + geom_density(data = subset(df3, status=="PASS"), aes(hwe), fill = as.factor(model), alpha=0.3) + labs(x = "HWE (p-value)", y = "Density")
    ggplot() + geom_density(data= subset(df3, !is.na(mqsb) & status=="PASS"), aes(mqsb), fill = as.factor(model), alpha=0.3) + labs(x = "MQSB", y = "Density") 
    ggplot() + geom_density(data= subset(df3, !is.na(mqb) & status=="PASS"), aes(mqb), fill = as.factor(model), alpha=0.3) + labs(x = "MQB", y = "Density")
    ggplot() + geom_density(data= subset(df3, !is.na(bqb) & status=="PASS"), aes(bqb), fill = as.factor(model), alpha=0.3) + labs(x = "BQB", y = "Density")
    ggplot() + geom_density(data= subset(df3, !is.na(rpb) & status=="PASS"), aes(rpb), fill = as.factor(model), alpha=0.3) + labs(x = "RPB", y = "Density")
    ggplot() + geom_density(data= subset(df3, !is.na(icb) & status=="PASS"), aes(icb), fill = as.factor(model), alpha=0.3) + labs(x = "ICB", y = "Density")
    dev.off()
    
    message("plotting attributes of failed sites...")
    pdf(file=args[3], onefile = T)
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(qual), fill = as.factor(model), alpha=0.3) + labs(x = "QUAL", y = "Density")
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(af), fill = as.factor(model), alpha=0.3) + labs(x = "Alt Allele Freq", y = "Density")
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(maf), fill = as.factor(model), alpha=0.3) + labs(x = "Minor allele freq", y = "Density")
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(ac), fill = as.factor(model), alpha=0.3) + labs(x = "Alt allele count", y = "Density")
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(an), fill = as.factor(model), alpha=0.3) + labs(x = "Total allele number", y = "Density")
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(mq), fill = as.factor(model), alpha=0.3) + labs(x = "MQ", y = "Density") + coord_cartesian(xlim = c(45, 60))
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(dp), fill = as.factor(model), alpha=0.3) + labs(x = "DP", y = "Density") + coord_cartesian(xlim = c(0, 12500))
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(vqslod), fill = as.factor(model), alpha=0.3) + labs(x = "VQSLOD", y = "Density")
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(maf), fill = as.factor(model), alpha=0.3) + labs(x = "MAF", y = "Density") 
    ggplot() + geom_density(data = subset(df3, status1=="Fail"), aes(hwe), fill = as.factor(model), alpha=0.3) + labs(x = "HWE (p-value)", y = "Density")
    ggplot() + geom_density(data= subset(df3, !is.na(mqsb) & status1=="Fail"), aes(mqsb), fill = as.factor(model), alpha=0.3) + labs(x = "MQSB", y = "Density") 
    ggplot() + geom_density(data= subset(df3, !is.na(mqb) & status1=="Fail"), aes(mqb), fill = as.factor(model), alpha=0.3) + labs(x = "MQB", y = "Density")
    ggplot() + geom_density(data= subset(df3, !is.na(bqb) & status1=="Fail"), aes(bqb), fill = as.factor(model), alpha=0.3) + labs(x = "BQB", y = "Density")
    ggplot() + geom_density(data= subset(df3, !is.na(rpb) & status1=="Fail"), aes(rpb), fill = as.factor(model), alpha=0.3) + labs(x = "RPB", y = "Density")
    ggplot() + geom_density(data= subset(df3, !is.na(icb) & status1=="Fail"), aes(icb), fill = as.factor(model), alpha=0.3) + labs(x = "ICB", y = "Density")
    dev.off()


**How much of the rare variants were captured by each model?**

    df3<-readRDS("merged_attributes.rds")
    dim(subset(df3, status1=="PASS" & model=="model1"))
    [1] 12111535   18
    dim(subset(df3, status1=="PASS" & model=="model1" & maf <= 0.05))
    [1] 4128878  18
    dim(subset(df3, status1=="PASS" & model=="model2"))
    [1] 12678713   18
    dim(subset(df3, status1=="PASS" & model=="model2" & maf <= 0.05))
    [1] 4464481  18

% of rare variants captured by: 
MODEL1: 34.09 %
MODEL2: 35.2 %




**Examine the culprit column categories for failed sites between both the models...**

    bcftools query -f '%ID\t%INFO/culprit\n' recalibrated_fail.1.2.vcf.gz > model1_failed_culprit.txt
    bcftools query -f '%ID\t%INFO/culprit\n' recalibrated_fail.3.2.vcf.gz >  model2_failed_culprit.txt
    
    perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f model1_failed_culprit.txt -c 1 -o model1_failed_culprit_count
    perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f model2_failed_culprit.txt -c 1 -o model2_failed_culprit_count


 |**Count**| |**Proportion**| 
:-----:|:-----:|:-----:|:-----:|:-----:
Culprit|Model1|Model2 |Model1|Model2 
BQB|140285|49344|2.430593336|0.948110255
DP|462278|672482|8.009479461|12.92126865
MQSB|4717650|3957503|81.73852266|76.0406367
RPB|451423|525129|7.821404538|10.08998439
total|5771636|5204458|100|100


Now perform a paired t-test to determine if the difference in the effect of culprits on failed sites between the models is significant...

    R
	culprit <- read.delim("U:/GATK/ars-ucdv1.2/culprit.txt", stringsAsFactors=FALSE)
    t.test(culprit$Model1,culprit$Model2.,paired=T)
    
    	Paired t-test
    
    data:  culprit$Model1 and culprit$Model2.
    t = 8.6577e-16, df = 3, p-value = 1
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -7.345857  7.345857
    sample estimates:
    mean of the differences 
       1.998401e-15 
    
    plot(culprit$Model1, culprit$Model2., pch = 16, xlab="Model1", ylab = "Model2")
    abline(0,1, col ="blue", lwd=2)
    cor(culprit$Model1, culprit$Model2.)
	[1] 0.9965416

p-value > 0.05, so we fail to reject the null hypothesis that the difference in means between the models is equal to zero. 

The culprit categories do not have a statistically significant effect on the failed sites between the two models.


#### No. of rare variants per chrom (MAF<0.05) that passed the TS99:  ####

    R	
    df3<-readRDS("merged_attributes.rds")
    df<-subset(df3, status1=="PASS" & model=="model1" & maf <= 0.05)
    write.table(df[,1], "model1_rare_passed", sep="\t", row.names=F, quote=F)
    df<-subset(df3, status1=="PASS" & model=="model2" & maf <= 0.05)
    write.table(df[,1], "model2_rare_passed", sep="\t", row.names=F, quote=F)
    

    gunzip -c run1/recalibrated.1.2.vcf.gz | grep -Fwf model1_rare_passed | perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -c 0 -o model1_rare_passed_per_chrom -e '#'
    gunzip -c run3/recalibrated.3.2.vcf.gz | grep -Fwf model2_rare_passed | perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -c 0 -o model2_rare_passed_per_chrom -e '#'


**chrom**|**Model1**|**Model2**
:-----:|:-----:|:-----:
1|285887|307812
2|228631|246842
3|196533|212904
4|212607|229033
5|215386|233261
6|196317|208889
7|178919|194067
8|175507|190200
9|184903|197781
10|164057|177615
11|164689|180124
12|165963|172083
13|118017|130213
14|145719|159930
15|140760|148971
16|130360|141127
17|137676|149306
18|84508|92455
19|63262|70765
20|129696|141201
21|89425|99743
22|100960|110192
23|113848|120356
24|106016|115033
25|66883|74307
26|89176|97726
27|75449|81813
28|73950|80077
29|93774|100655


**Now perform paired t-test to understand if the results are statistically significant...**
	
	R
	df <- read.delim("U:/GATK/ars-ucdv1.2/rare_var_per_chrom.txt", stringsAsFactors=FALSE)
    t.test(df$Model1,df$Model2,paired=T)
    
    	Paired t-test
    
    data:  df$Model1 and df$Model2
    t = -14.678, df = 28, **p-value = 1.121e-14**
    alternative hypothesis: true difference in means is not equal to 0
    95 percent confidence interval:
     -13187.523  -9957.512
    sample estimates:
    mean of the differences 
      -11572.52 
	
	plot(df$Model1,df$Model2, pch = 16, xlab="Model1", ylab = "Model2", main = "Rare variants (MAF < 0.05) per Chromosome")
	abline((lm(df$Model2~df$Model1), col ="blue", lwd=2))


The number of rare variants per chromosome that passed the TS99 under trained Model2 are significantly higher than that of Model1.

There seems to be no chromosome bias with rare variants per chromosome analysis...

 
![](https://i.imgur.com/9Y809C2.png)



**How many variants were picked up +/- 1kb of each of the truth training site in both the models? Is it picking more sites for model2 than model1?**

Its better do this analyses with bed coordinates rather than full vcf files...
	
	# make bovineHD bed
	bcftools query -f'%CHROM\t%POS0\t%END\t%ID\n' bovineHD.vcf.gz > modelling/analysis/bovineHD.bed
	# Add 1kb to each site left and right	
	/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedops --range 1000 -u modelling/analysis/bovineHD.bed > modelling/analysis/bovineHD.plus.1kb.bed

Now using a bash script count the number of passed variants that were picked up by each model 

	#!/usr/bin/sh
	while read i; 
	do bcftools view -r $i $1 | bcftools +counts | grep 'SNP' >> $2_counts_per_site
	done <bovineHD.plus.1kb.bed

Now running the script to get the count...

	sh counts_per_interval.sh /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/gatk/modelling/run1/recalibrated_pass.1.2.vcf.gz model1

	sh counts_per_interval.sh /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/gatk/modelling/run3/recalibrated_pass.3.2.vcf.gz model2

Still running...

**Are there any functional difference between the passed sites of both the models?**

Lets annotate the passed sites vcf files of each model using the custom ARSUCD1.2 SnpEff annotation database and check the stats...


