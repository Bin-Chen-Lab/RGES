
## RGES

This pipeline is to compute Reverse Gene Expression Score (RGES) published by Chen B. et al (Nature Communications, 2017). In this work, we found that the drug's ability to reverse cancer gene expression (RGES) correlates to its efficacy (IC50) in cancer cell lines.

Since the code for this publication is  complicated,  we would recommend to play with the example at the begining. This example includes computing and summarizing RGES, and correlating sRGES with drug efficacy. You can simply replace the disease signature with yours to compute sRGES for the disease of interest.

To run the example, you need 
1) download example data from <https://www.synapse.org/#!Synapse:syn6182429/files/>
2) change workspace in runRGESExample.R
3) run runRGESExample.R. 
```{r cars}
Rscript("runRGESExample")
```

To run the entire code, you need 
1) download data from <https://www.synapse.org/#!Synapse:syn6182429/files/>
2) change workspace in workflow.R
3) run workflow.R. It is a very time-consuming and complicated pipeline, so we recommend to run the specific component at one time.




