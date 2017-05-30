
## RGES

This pipeline is to compute Reverse Gene Expression Score (RGES) published by Chen B. et al (Nature Communications, 2017). In this work, we found that the drug's ability to reverse cancer gene expression (RGES) correlates to its efficacy (IC50) in cancer cell lines.

To run the code, you need 
1) download data from <https://www.synapse.org/#!Synapse:syn6182429/files/>
2) change workspace in workflow.R
3) run workflow.R. The entire pipeline is rather complicated, so we would recommend to play with the example at the begining.
```{r cars}
Rscript("runRGESExample")
```



