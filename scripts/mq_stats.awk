#!/usr/bin/awk -f
#array to convert ASCI to PHRED
BEGIN {for (n=0; n<256; n++) ord[sprintf("%c",n)]=n;
       OFS="\t"
    }
{   
    cov=0; mq0=0; sm=0; 
    n_samples = (NF-3)/4
    for (j=7; j<=NF; j+=4) 
    {
         split($j, a, ""); 
         #print j;
         for (i=1; i <= length(a); i++) 
           {
            s=ord[a[i]]-33;
            cov++ 
            if (s>0) {sm=sm+s^2}
            else {mq0++} 
            
            };
    };
    if (cov-mq0 == 0) {rmsmq="NA"} else {rmsmq=sqrt(sm*1.0/(cov-mq0))};
    print $1, $2, cov, rmsmq, mq0
        
}



