These scripts allow cross-species comparison of DEGs using gene family information from Phytozome. 

First, [we created](https://github.com/plant-plasticity/Evolutionary-flexibility-in-flooding-response-2019/blob/master/Gene-families-Phytozome/GeneFamiliesPhytozome.R) [a gene family file](https://github.com/plant-plasticity/Evolutionary-flexibility-in-flooding-response-2019/blob/master/Gene-families-Phytozome/phytozomev11_gene_families_addSl3_1.csv) that contains the gene family information for our target species available in Phytozome and adds gene models introduced in the newer tomato annotation.

To link S. pennellii (not present in Phytozome) we use a [separate Sl-to-Sp list](https://github.com/plant-plasticity/Evolutionary-flexibility-in-flooding-response-2019/blob/master/Gene-families-Phytozome/ITAG3.10_Spenn.xlsx) based on synteny and orthology. For more details, see Supplementary Methods.

Then, [we compare DEGs across species](https://github.com/plant-plasticity/Evolutionary-flexibility-in-flooding-response-2019/blob/master/Gene-families-Phytozome/DEGs_Comparisons_species.R). The input is the DEG_AllContrasts.csv file from limma-voom analysis, and the output in this example where all submergence-upregulated genes are compared in all species, is a Venn diagram of the overlapping upregulated genes.

