    
#Metagenomics technology comparison

##EggNOG Annotations Data Processing 

This is an attempt to find out if there is a significant difference in sequencing technology Pacbio(p) vs Illumina(i) in the field of metagenomics.


- The samples from cattle rumen were sequenced and assembled using both the technologies.
 
- The major constituents of the microbiome includes Bacteria(b), Archaea(a) and Eukaryota(e)...
 
- Both the metagenome assemblies were then annotated using eggnog-mapper/bin/diamond blastp
 
- The fields that we are interested in the annotation files are COG, KEGG and GO terms. 


#### working directory: 
    > getwd()
    [1] "/scinet01/gov/usda/ars/scinet/project/rumen_longread_metagenome_assembly/analysis/eggnog"

    ia<-read.table("illumina_megahit.Archaea.eggnog.emapper.annotations",sep="\t",header=F,na.strings="",stringsAsFactors=F,quote="")
    
    header<-c("query_name","seed_eggNOG_ortholog","seed_ortholog_evalue","seed_ortholog_score","predicted_gene_name","GO_terms","KEGG_KOs","BiGG_reactions","Annotation_tax_scope","OGs","bestOG|evalue|score","COG_cat","eggNOG_annot")
    
    ie<-read.table("illumina_megahit.Eukaryota.eggnog.emapper.annotations",sep="\t",header=F,na.strings="",stringsAsFactors=F,quote="")

    pa<-read.table("pacbio_final.Archaea.eggnog.emapper.annotations",sep="\t",header=F,na.strings="",stringsAsFactors=F,quote="")
    names(pa)<-header
    
	pb<-read.table("pacbio_final.Bacteria.eggnog.emapper.annotations",sep="\t",header=F,na.strings="",stringsAsFactors=F,quote="")
    names(pb)<-header
    
	pe<-read.table("pacbio_final.Eukaryota.eggnog.emapper.annotations",sep="\t",header=F,na.strings="",stringsAsFactors=F,quote="")
    names(pe)<-header
    library(dplyr)
    cog_ia<-ia %>% group_by(COG_cat) %>% summarise(ia_cog=length(query_name))
    cog_pa<-pa %>% group_by(COG_cat) %>% summarise(pa_cog=length(query_name))
    cog_ib<-ib %>% group_by(COG_cat) %>% summarise(ib_cog=length(query_name))
    cog_pb<-pb %>% group_by(COG_cat) %>% summarise(pb_cog=length(query_name))
    cog_ie<-ie %>% group_by(COG_cat) %>% summarise(ie_cog=length(query_name))
    cog_pe<-pe %>% group_by(COG_cat) %>% summarise(pe_cog=length(query_name))

    cog_all<-Reduce(function(x,y) merge(x, y, all=T), list(cog_ia,cog_pa,cog_ib,cog_pb,cog_ie,cog_pe))
    
	cog<-cog_all[complete.cases(cog_all),]

    cog_Archea<-merge(cog_ia,cog_pa,by="COG_cat")
    cog_Bacteria<-merge(cog_ib,cog_pb,by="COG_cat")
    cog_Eukaryote<-merge(cog_ie,cog_pe,by="COG_cat")

    save(ia,ib,ie,pa,pb,pe, file="annotations.RData")
	save(cog_all,cog_Archea,cog_Bacteria,cog_Eukaryote,cog,file="COG_counts.RData")


    
The COG data is now ready for analysis


####BACTERIA####

**Wilcoxon signed rank test:**
 
Wilcoxon signed-rank test is a non-parametric statistical hypothesis test used to compare repeated measurements on the same sample to assess whether their population mean ranks differ. 

Assumptions: 
1. Data are paired and come from the same population
2. each pair is chosen randomly and independently
3. the data are measured on at least an interval scale 

H<sub>0</sub>: The median difference between the sequence assemblies (observations) is zero.

H<sub>1</sub>: median difference between the sequence assemblies is different

*Bonferroni correction = 0.05/97 = 5.1e-4*



	>  wilcox.test(cog_Bacteria$ib_cog,cog_Bacteria$pb_cog, paired=T)

        Wilcoxon signed rank test with continuity correction

	data:  cog_Bacteria$ib_cog and cog_Bacteria$pb_cog
	V = 4256.5, p-value = 2.124e-15
	alternative hypothesis: true location shift is not equal to 0

**Conclusion:** After adjusting for multiple correction, there is a significant evidence (P-value=2.124e-15) that the median difference between the sequence assemblies is not zero.


**Chi-square test for Independence:**

H<sub>0</sub>: The sequencing technology (Illumina, Pacbio) and COG category are independent

H<sub>1</sub>: The sequencing technology and COG category are associated

    > chisq.test(cog_all$ib_cog,cog_all$pb_cog)
    
    Pearson's Chi-squared test
    
    data:  cog_all$ib_cog and cog_all$pb_cog
    X-squared = 6621.5, df = 5772, p-value = 2.091e-14
    
    Warning message:
    In chisq.test(cog_all$ib_cog, cog_all$pb_cog) :
      Chi-squared approximation may be incorrect
   

    > chisq.test(cog_all$ib_cog,cog_all$pb_cog,simulate.p.value=T, B=10000)
    
    Pearson's Chi-squared test with simulated p-value (based on 10000
    replicates)
    
    data:  cog_all$ib_cog and cog_all$pb_cog
    X-squared = 6621.5, df = NA, p-value = 9.999e-05


**Conclusion**: After multiple correction (Bonferroni correction=5.1e-4), there is sufficient evidence that the sequencing technology and COG category are associated (p-value=9.999e-05). 

However, it should be noted that the association between variables can also exist as a result of large number of datapoints.

Now to compare the proportions of COG categories across the two groups using Chisq test for homogeneity, the data set should be transposed...

    cog_Bacteria<-cog_Bacteria[1:nrow(cog_Bacteria)-1,] ## remove the NA group
    cog_Bacteria$COG_cat<-gsub(", ","", cog_Bacteria$COG_cat)
    row.names(cog_Bacteria)<-cog_Bacteria$COG_cat
    cog_Bacteria_t <- as.data.frame(t(cog_Bacteria),stringsAsFactors=F) # Transpose the dataset
    library(magrittr)
	cog_Bacteria_t<-cog_Bacteria_t[2:3,]
    cog_Bacteria_t %<>% mutate_if(is.character,as.numeric) # change all the columns to numeric



**Chi-square test for Homogeneity:**

Is the distribution of COG categories same across the two sequence assemblies (Illumina, Pacbio)?

H<sub>0</sub>: Distribution of COG categories are the same for both Illumina and Pacbio assemblies

H<sub>1</sub>: Distributions are not same

    > chisq.test(cog_Bacteria_t)
    
    Pearson's Chi-squared test
    
    data:  cog_Bacteria_t
    X-squared = 4127.8, df = 96, p-value < 2.2e-16
    
    Warning message:
    In chisq.test(cog_Bacteria_t) : Chi-squared approximation may be incorrect

When the validity of the asymptotic approximation cannot be trusted, the permutation approach can be implemented in the chisq.test() function...

    > chisq.test(cog_Bacteria_t, simulate.p.value=T, B=10000)
    
    Pearson's Chi-squared test with simulated p-value (based on 10000
    replicates)
    
    data:  cog_Bacteria_t
    X-squared = 4127.8, df = NA, p-value = 9.999e-05

**Conclusion**: There is significant evidence after Bonferroni correction (5.1e-4) that the distribution of COG categories are not the same for both the sequence assemblies (P-value=9.999e-05)
 

**ARCHAEA**

*Bonferroni correction = 0.05/39 = 0.0012*

    > wilcox.test(cog_Archea$ia_cog,cog_Archea$pa_cog, paired=T)
    
    Wilcoxon signed rank test with continuity correction
    
    data:  cog_Archea$ia_cog and cog_Archea$pa_cog
    V = 774, p-value = 8.646e-08
    alternative hypothesis: true location shift is not equal to 0
    
    Warning messages:
    1: In wilcox.test.default(cog_Archea$ia_cog, cog_Archea$pa_cog, paired = T) :
      cannot compute exact p-value with ties
    2: In wilcox.test.default(cog_Archea$ia_cog, cog_Archea$pa_cog, paired = T) :
      cannot compute exact p-value with zeroes


Chisq test:    

    > chisq.test(cog_Archea_t)
    
    Pearson's Chi-squared test
    
    data:  cog_Archea_t
    X-squared = 297.9, df = 38, p-value < 2.2e-16
    
    Warning message:
    In chisq.test(cog_Archea_t) : Chi-squared approximation may be incorrect
    
    
    > chisq.test(cog_Archea_t, simulate.p.value=T, B=10000)
    
    Pearson's Chi-squared test with simulated p-value (based on 10000
    replicates)
    
    data:  cog_Archea_t
    X-squared = 297.9, df = NA, p-value = 9.999e-05



**EUKARYOTE**

*Bonferroni correction = 0.05/48 = 0.001*

    > wilcox.test(cog_Eukaryote$ie_cog,cog_Eukaryote$pe_cog, paired=T)
    
    Wilcoxon signed rank test with continuity correction
    
    data:  cog_Eukaryote$ie_cog and cog_Eukaryote$pe_cog
    V = 1014.5, p-value = 2.085e-08
    alternative hypothesis: true location shift is not equal to 0
    
    Warning messages:
    1: In wilcox.test.default(cog_Eukaryote$ie_cog, cog_Eukaryote$pe_cog,  :
      cannot compute exact p-value with ties
    2: In wilcox.test.default(cog_Eukaryote$ie_cog, cog_Eukaryote$pe_cog,  :
      cannot compute exact p-value with zeroes


###Graphical display of the available data


    library(ggplot2)
    library(tidyr)
    library(dplyr)
	load("COG_counts.RData")

    names(cog_Bacteria)<-c("COG_cat","Illumina", "Pacbio")
    
    cog_Bacteria <- cog_Bacteria %>% mutate(Illumina.prop = (Illumina*100)/sum(Illumina))
    cog_Bacteria <- cog_Bacteria %>% mutate(Pacbio.prop = (Pacbio*100)/sum(Pacbio))

	names(cog_Archea)<-c("COG_cat","Illumina", "Pacbio")
    cog_Archea <- cog_Archea %>% mutate(Illumina.prop = (Illumina*100)/sum(Illumina))
    cog_Archea <- cog_Archea %>% mutate(Pacbio.prop = (Pacbio*100)/sum(Pacbio))
    
	names(cog_Eukaryote)<-c("COG_cat","Illumina", "Pacbio")
    cog_Eukaryote <- cog_Eukaryote %>% mutate(Illumina.prop = (Illumina*100)/sum(Illumina))
    cog_Eukaryote <- cog_Eukaryote %>% mutate(Pacbio.prop = (Pacbio*100)/sum(Pacbio))

    cBl<-gather(cog_Bacteria, "category", "count", 4:5)
    
    bp<-ggplot(cBl, aes(x=category, y=count, fill=COG_cat))+
    geom_bar(stat = "identity")
    bp + coord_polar("y", start=0)
    
	## Too many categories and the plot does not look good
	## Try to condense the categories
	## I got the following condensed categories from [CloVR](http://clovr.org/docs/clusters-of-orthologous-groups-cogs/)
    poor<-c("R","S")
    cell<-c("D","M","N","O","T","U","V","W","Y","Z","DZ","MOT","MU","NO","NOU","NT","NU","OT","OU","OW","TU","TV","TW","TZ","UW","UY","VW")
    info_sto<-c("A","B","J","K","L","AJ","BK","KL")
    met<-c("C","E","F","G","H","I","P","Q","CE","CG","CH","CP","EF","EFP", "EG", "EGP", "EH", "EI","EP", "EQ","FG","FH","FP","HI","HP","IQ","PQ")
    
	## all categories in the dataset
    cog_cat<-cog_Bacteria$COG_cat 
    
	## What about the remaining overlapping categories?
	## below is a function to subtract character vectors
    cNew <- unlist(sapply(cog_cat[!duplicated(cog_cat)],
    function(item, tab1, tab2) {
    rep(item,
    tab1[item] - ifelse(item %in% names(tab2), tab2[item], 0))
    }, tab1=table(cog_cat), tab2=table(cell)))
    
    cNew <- unlist(sapply(cNew[!duplicated(cNew)],
    function(item, tab1, tab2) {
    rep(item,
    tab1[item] - ifelse(item %in% names(tab2), tab2[item], 0))
    }, tab1=table(cNew), tab2=table(met)))
    
    cNew <- unlist(sapply(cNew[!duplicated(cNew)],
    function(item, tab1, tab2) {
    rep(item,
    tab1[item] - ifelse(item %in% names(tab2), tab2[item], 0))
    }, tab1=table(cNew), tab2=table(info_sto)))
    
    cNew <- unlist(sapply(cNew[!duplicated(cNew)],
    function(item, tab1, tab2) {
    rep(item,
    tab1[item] - ifelse(item %in% names(tab2), tab2[item], 0))
    }, tab1=table(cNew), tab2=table(poor)))
    
	## create a dataframe to populate with the counts of condensed categories
    func_cat<-data.frame(matrix(vector(), 5,5, dimnames=list(c(),c("Category", "Illumina", "Pacbio", "Illumina.prop", "Pacbio.prop"))), stringsAsFactors = F)
	## Start populating it manually
    func_cat$Category<-c("Cellular.processes", "Metabolism", "Information.storage","Poorly.characterized", "Combination")
    subset(cog_Bacteria, COG_cat %in% poor) %>% summarise(Illumina=sum(Illumina),Pacbio=sum(Pacbio))
      Illumina Pacbio
    1   390184 195710
    subset(cog_Bacteria, COG_cat %in% cell) %>% summarise(Illumina=sum(Illumina),Pacbio=sum(Pacbio))
      Illumina Pacbio
    1   357969 164688
    subset(cog_Bacteria, COG_cat %in% info_sto) %>% summarise(Illumina=sum(Illumina),Pacbio=sum(Pacbio))
      Illumina Pacbio
    1   317437 150515
    subset(cog_Bacteria, COG_cat %in% met) %>% summarise(Illumina=sum(Illumina),Pacbio=sum(Pacbio))
      Illumina Pacbio
    1   453408 245850
    subset(cog_Bacteria, COG_cat %in% cNew) %>% summarise(Illumina=sum(Illumina),Pacbio=sum(Pacbio))
      Illumina Pacbio
    1 6621   3047

    ## Looks like there are NAs in the categories which should be added to the poorly characterized group
	func_cat$Pacbio<-c(164688,245850,150515,283507,3047)
	func_cat$Illumina<-c(357969,453408,317437,563814,6621)
	func_cat <- func_cat %>% mutate(Illumina.prop = (Illumina*100)/sum(Illumina))
    func_cat <- func_cat %>% mutate(Pacbio.prop = (Pacbio*100)/sum(Pacbio))
    fc<-gather(func_cat, "Technology", "Proportion", 4:5)
	fc$Technology<-gsub(".prop","",fc$Technology)
	ggplot(fc, aes(x=Technology, y=Proportion, fill=Category))+ geom_bar(stat = "identity") + ggtitle("Annotation:COG categories:Bacteria")
	#save the plot as COG_categories_Bacteria.jpeg

	## Scatter plot to show the correlation between assemblies and COG annotations 
	library(gridExtra)
	b<-ggplot(cog_Bacteria, aes(Illumina.prop, Pacbio.prop, label=COG_cat)) + geom_point(color="#33CC33") + geom_smooth(method=lm, se=F) + ggtitle("COG_categories:Bacteria") + xlab("Illumina") + ylab("Pacbio") + geom_text_repel(data = subset(cog_Bacteria, Illumina.prop > 2.5 | Pacbio.prop > 2.5))
    a<-ggplot(cog_Archea, aes(Illumina.prop, Pacbio.prop, label=COG_cat)) + geom_point(color="#FF6666") + geom_smooth(method=lm, se=F) + ggtitle("COG_categories:Archaea") + xlab("Illumina") + ylab("Pacbio") + geom_text_repel(data = subset(cog_Archea, Illumina.prop > 2.5 | Pacbio.prop > 2.5)) 
    e<-ggplot(cog_Eukaryote, aes(Illumina.prop, Pacbio.prop, label=COG_cat)) + geom_point(color="#0099FF") + geom_smooth(method=lm, se=F) + ggtitle("COG_categories:Eukaryota") + xlab("Illumina") + ylab("Pacbio") + geom_text_repel(data = subset(cog_Eukaryote, Illumina.prop > 2.5 | Pacbio.prop > 2.5)) 
    
    grid.arrange(b,a,e,ncol=3) ## plot saved as Scatterplot_COG.jpeg


	## Scatterplot of assemblies COG annotation categories showing all the three groups of organisms 
    cog<-cog_all[complete.cases(cog_all),]	## select the COG categories for which we have counts in all 3 groups
    names(cog)<-c("COG","Illumina.Archaea","Pacbio.Archaea","Illumina.Bacteria","Pacbio.Bacteria","Illumina.Eukaryota","Pacbio.Eukaryota")
    cog_prop <- cog %>% mutate_each(funs((.)*100/sum(.)), -COG) ## calculate their relative proportions
    cogl<-gather(cog_prop, "Key", "Value", 2:7)
    cogl$Technology <- gsub("\\..*","",cogl$Key)
    cogl$Group <- gsub(".*\\.","",cogl$Key)
    cogl$Key<-NULL
    cogl_org<-spread(cogl,key = Technology, value=Value, convert=T)
    ggplot(cogl_org, aes(Illumina,Pacbio, color=Group,shape=Group)) + geom_point() + geom_smooth(method=lm, se=F) + ggtitle("COG categories") # plot saved as Scatterplot_COG_all.jpeg
    ## Add labels on the plot
	library(ggrepel)
    ggplot(cogl_org, aes(Illumina,Pacbio, color=Group,shape=Group, label=COG)) + geom_point() + geom_smooth(method=lm, se=F) + ggtitle("COG categories")+geom_text_repel( data = subset(cogl_org, Illumina > 5 | Pacbio > 5)) # plot saved as Scatterplot_COG_all.labs.jpeg


Updating the document:

The grouped COG categories were resolved by using a custom script to split the grouped data into individual categories.

The script is as follows: 

    #!/usr/bin/perl
    # This is a script to split grouped COG category counts to individual category counts 
    
    use strict;
    use warnings;
    use Data::Dumper qw(Dumper);
    
    if (scalar @ARGV != 2 ) {
    print STDERR "Usage: splitCOG.pl <tab-delimited COG counts file> <output filename>\n";
    exit(1);
    }
    
    my ($input, $output) = @ARGV;
    
    #my @main_cat = ("D","M","N","O","T","U","V","W","Y","Z","A","B","J","K","L","C","E","F","G","H","I","P","Q","R","S");
    open(my $IN, "< $input") || die "Could not open input file: $input!\n";
    open(my $OUTFILE, "> $output" );	
    print {$OUTFILE} "COG_category\tIA\tPA\tIB\tPB\tIE\tPE\n";
    
    my @cogIdx;
    my %cogCls;
    
    my $head = <$IN>;
    chomp $head;
    my @hsegs = split(/\t/, $head);
    shift @hsegs;
    @cogIdx = @hsegs;
    #print Dumper \@cogIdx;
    
    while(my $line = <$IN>){
    	chomp $line;
    	my @cols = split(/\t/, $line);
    	my $cog = shift @cols;
    	for(my $x = 0; $x < scalar(@cols); $x++){
    			my $cls = $cogIdx[$x];
    			$cogCls{$cog}->{$cls} = $cols[$x];
    	}			
    }
    #print Dumper \%cogCls;
    
    my %newcogCls = map { $_ => $cogCls{$_} } grep { /^[A-Z]$/ } keys %cogCls;
    #print Dumper \%newcogCls;
    
    foreach my $cat (sort keys %cogCls)  {
    if ($cat =~ /.*\,.*/) {
    		my @subcat = split(/, /, $cat);
    		foreach my $i(@subcat){
    			foreach my $cls(@cogIdx){
    			if ($cogCls{$cat}->{$cls} ne "NA"){
    				$newcogCls{$i}->{$cls} += sprintf('%.0f', ($cogCls{$cat}->{$cls})/scalar(@subcat));
    			} else {next;
    		 }
    		}
    	}
    	} else {next;} 
    }
    #print Dumper \%cogCls;
    #print Dumper \%newcogCls;
    
    
    
    foreach my $cat (sort keys %newcogCls) {
    	my @val;
    	foreach my $cls (@cogIdx) {
    	push (@val, $newcogCls{$cat}{$cls});
    	}
    	printf {$OUTFILE} join("\t", $cat, @val) . "\n";	
    }
    print "done\n";
    close $IN;
    close $OUTFILE;
    exit;


The output file is loaded in to R as final_cog_counts
    
    final_cog_counts <- read.delim("U:/metagenomics/final_cog_counts", stringsAsFactors=FALSE)
    names(final_cog_counts)<-c("COG","Illumina.Archaea","Pacbio.Archaea","Illumina.Bacteria","Pacbio.Bacteria","Illumina.Eukaryota","Pacbio.Eukaryota")
    cog_prop_final <- final_cog_counts %>% mutate_each(funs((.)*100/sum(.)), -COG)
    
    coglf$Technology <- gsub("\\..*","",coglf$Key)
    > coglf$Group <- gsub(".*\\.","",coglf$Key)
    > coglf$Key<-NULL
    > coglf_org<-spread(coglf,key = Technology, value=Value, convert=T)
    > ggplot(coglf_org, aes(Illumina,Pacbio, color=Group,shape=Group)) + geom_point() + geom_smooth(method=lm, se=F) + ggtitle("COG categories")
    > library(ggrepel)
    > ggplot(coglf_org, aes(Illumina,Pacbio, color=Group,shape=Group, label=COG)) + geom_point() + geom_smooth(method=lm, se=F) + ggtitle("COG categories")+geom_text_repel( data = subset(coglf_org, Illumina > 5 | Pacbio > 5))
    > ggplot(coglf_org, aes(Illumina,Pacbio, color=Group,shape=Group, label=COG)) + geom_point() + geom_smooth(method=lm) + ggtitle("COG categories")+geom_text_repel( data = subset(coglf_org, Illumina > 5 | Pacbio > 5))
    

Now try to plot the individual plots to better show the confidence intervals.

	library(gridExtra)
	b<-ggplot(cog_prop_final, aes(Illumina.Bacteria, Pacbio.Bacteria, label=COG)) + geom_point(color="#33CC33") + geom_smooth(method=lm) + ggtitle("COG_categories:Bacteria") + xlab("Illumina") + ylab("Pacbio") + geom_text_repel(data = subset(cog_prop_final, Illumina.Bacteria > 2.5 | Pacbio.Bacteria > 2.5))
    a<-ggplot(cog_prop_final, aes(Illumina.Archaea, Pacbio.Archaea, label=COG)) + geom_point(color="#FF6666") + geom_smooth(method=lm) + ggtitle("COG_categories:Archaea") + xlab("Illumina") + ylab("Pacbio") + geom_text_repel(data = subset(cog_prop_final, Illumina.Archaea > 2.5 | Pacbio.Archaea > 2.5)) 
    e<-ggplot(cog_prop_final, aes(Illumina.Eukaryota, Pacbio.Eukaryota, label=COG)) + geom_point(color="#0099FF") + geom_smooth(method=lm) + ggtitle("COG_categories:Eukaryota") + xlab("Illumina") + ylab("Pacbio") + geom_text_repel(data = subset(cog_prop_final, Illumina.Eukaryota > 2.5 | Pacbio.Eukaryota > 2.5)) 
    
    grid.arrange(b,a,e,ncol=3) ## plot saved as Scatterplot_COG.jpeg


Calculate correlation between the technologies for each category and taxa:

cor(cog_prop_final$Illumina.Archaea, cog_prop_final$Pacbio.Archaea)
[1] 0.9825274

cor(cog_prop_final$Illumina.Bacteria, cog_prop_final$Pacbio.Bacteria)
[1] 0.9968469

cor(cog_prop_final$Illumina.Eukaryota, cog_prop_final$Pacbio.Eukaryota)
[1] 0.5875916

> t.test(cog_prop_final$Illumina.Archaea, cog_prop_final$Pacbio.Archaea, paired = TRUE) 

	Paired t-test

data:  cog_prop_final$Illumina.Archaea and cog_prop_final$Pacbio.Archaea
t = 3.1718e-19, df = 23, p-value = 1
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.3977562  0.3977562
sample estimates:
mean of the differences 
           6.098637e-20 

> wilcox.test(cog_prop_final$Illumina.Archaea, cog_prop_final$Pacbio.Archaea, paired = TRUE)

	Wilcoxon signed rank test

data:  cog_prop_final$Illumina.Archaea and cog_prop_final$Pacbio.Archaea
V = 149, p-value = 0.9888
alternative hypothesis: true location shift is not equal to 0

> t.test(cog_prop_final$Illumina.Bacteria, cog_prop_final$Pacbio.Bacteria, paired = TRUE) 

	Paired t-test

data:  cog_prop_final$Illumina.Bacteria and cog_prop_final$Pacbio.Bacteria
t = -1.958e-15, df = 23, p-value = 1
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1810611  0.1810611
sample estimates:
mean of the differences 
          -1.713774e-16 

> wilcox.test(cog_prop_final$Illumina.Bacteria, cog_prop_final$Pacbio.Bacteria, paired = TRUE)

	Wilcoxon signed rank test

data:  cog_prop_final$Illumina.Bacteria and cog_prop_final$Pacbio.Bacteria
V = 160, p-value = 0.7898
alternative hypothesis: true location shift is not equal to 0

> t.test(cog_prop_final$Illumina.Eukaryota, cog_prop_final$Pacbio.Eukaryota, paired = TRUE) 

	Paired t-test

data:  cog_prop_final$Illumina.Eukaryota and cog_prop_final$Pacbio.Eukaryota
t = 1.184e-16, df = 23, p-value = 1
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1.871407  1.871407
sample estimates:
mean of the differences 
           1.071101e-16 

> wilcox.test(cog_prop_final$Illumina.Eukaryota, cog_prop_final$Pacbio.Eukaryota, paired = TRUE)

	Wilcoxon signed rank test

data:  cog_prop_final$Illumina.Eukaryota and cog_prop_final$Pacbio.Eukaryota
V = 130, p-value = 0.5838
alternative hypothesis: true location shift is not equal to 0
