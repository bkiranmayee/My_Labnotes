IGC filtered variants haplotype block analysis
==============================================

IGC alternate haplotype variants filtered by their QUAL, AAF and total allele number were lifted over to new assembly ARS-UCDv1.2.

The liftover failed coordinates were arbitrarily given new assembly coordinates. The vcf file is at /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/igc_filtered_all_ars-ucdv1.2.vcf.gz

working directory: */mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/igc_haplotype_analysis*   

 plink was used to calculate pairwise LD between all the 5002 SNPs...

    module load plink/1.90b4.4-2017-05-21
    
    plink --vcf igc_filtered_all_ars-ucdv1.2.vcf.gz --cow --keep-allele-order --make-bed --out igc_haplotype_analysis/igc
    
    plink --bfile igc_haplotype_analysis/igc --cow --r2 --ld-window 999999 --ld-window-kb 160000 --ld-window-r2 0 --out igc_haplotype_analysis/igc-ld
    
    cat igc-ld.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > igc.ld.summary



#### Plotting in R: ####

    dfr <- read.delim("igc.ld.summary",sep="",header=T,check.names=F,stringsAsFactors=F)
    colnames(dfr) <- c("dist","rsq")
    dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=100000))
    dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
    dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
    + end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
    + mid=start+((end-start)/2))
    p<-ggplot()+
    +   geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour="grey20")+
    +   geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
    +   labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+ theme_bw()
    pdf("ld.scatter.pdf")
    p
    Warning messages:
    1: Removed 1 rows containing missing values (geom_point).
    2: Removed 1 rows containing missing values (geom_path).
    dev.off()
    pdf
      2
    summary(dfr)
      dist   rsqdistc
     Min.   :   0   Min.   :0.00000   (-1,4999]: 325344
     1st Qu.:   43373   1st Qu.:0.01929   (4999,9999]  : 268514
     Median :  112988   Median :0.07717   (9999,14999] : 251512
     Mean   :  131067   Mean   :0.14648   (14999,19999]: 225293
     3rd Qu.:  201450   3rd Qu.:0.20459   (19999,24999]: 207646
     Max.   :54174401   Max.   :1.00000   (Other)  :6258038
      NA's :  7


#### Haplotype block analysis: ####

    plink --bfile igc no-pheno-req --blocks --blocks-max-kb 200 --out igc-bl
    R
    dfr <- read.delim("igc-bl.blocks.det",sep="",header=T,check.names=F,stringsAsFactors=F)
    colnames(dfr) <- tolower(colnames(dfr))
    
    # ld block density
    p <- ggplot(dfr,aes(x=kb))+
      geom_density(size=0.5,colour="grey40")+
      labs(x="LD block length (Kb)",y="Density")+
      theme_bw()
    
    ggsave("igc-ld-blocks.png",p,height=8,width=8,units="cm",dpi=250)

![IGC_haploblock_density](https://i.imgur.com/6kTHFPS.png)