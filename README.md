# 2023Autumn_Ginger

The analysis is consist of 2 parts:  
1. Snakemake pipeline consist of rules for:  
    1-1. QC fastq  
    1-2. Download reference (genome, gtf from NCBI)  
    1-3. Perform Ginger to Arabidopsis blastp  
    1-4. Quantification using RSEM

2. R Markdowns by Batch   
    2-1. Exploratory data analysis (visualize expression profiles)  
    EDA_batch1.Rmd  
    EDA_batch2.Rmd

    2-2. Differential expression analysis  
    DEA_batch1.Rmd  
    DEA_batch2.Rmd

    2-3. GO enrichment analysis (with PlantGSEA)  
    GO_enrichment_analysis_batch1.Rmd   
    GO_enrichment_analysis_batch2.Rmd
    