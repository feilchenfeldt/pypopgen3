#!/usr/bin/awk -f
# This script takes a samtools mpileup as input and
# and outputs CHROM, POS, cov, rmsmq, mq0
#USAGE:
#samtools mpileup --output-MQ --count-orphans sample1.bam [sample2.bam [sample3.bam [...]]] | mq_stats.awk > mapping_stats.tsv

# array to convert ASCI to PHRED
BEGIN {print "CHROM\tPOS\tcov\trmsmq\tmq0";
       for (n=0; n<256; n++) ord[sprintf("%c",n)]=n;
       OFS="\t"
    }
{   
    cov=0; mq0=0; sm=0; 
    n_samples = (NF-3)/4
    for (j=7; j<=NF; j+=4)
    {
     #alternative way to calculate coverage
     #use this OR cov++ below; they give virtually the same results
     #but can be off by 1 or 2 probably due to some reads being excluded
     #in the calculation of samtools
     #cov = cov + $(j-3)
        if (length($j) > 0 && $j != "*" && $(j-3)>0) 
        {
             split($j, a, ""); 
             #print j;

             for (i=1; i <= length(a); i++) 
               {
                s=ord[a[i]]-33;
                #this would be an alternative to determine coverage, but we now take the
                #info directly from the column above
                cov++ 
                if (s>0) {sm=sm+s^2}
                else {mq0++} 
                
                };
        };
    };
    if (cov-mq0 == 0) {rmsmq="NA"} else {rmsmq=sqrt(sm*1.0/(cov-mq0))};
    print $1, $2, cov, rmsmq, mq0
        
}



