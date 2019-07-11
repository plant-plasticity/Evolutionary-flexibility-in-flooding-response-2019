The ATAC-seq workflow was composed by Marko Bajic (marko.bajic@emory.edu or marko.bajic25@gmail.com) and consists of the following files


"Workflow-Overview.pdf"
This is a visual representation of the steps between obtaining sequencing reads and identifying replicated peaks between different cell types or conditions among which to look for accessibility level changes


"ATAC-Seq-Pipeline-Overview.pdf"
This is a more descriptive visual representation of a large amount of the ATAC-seq Data Processing


"ATACseq-Bioinformatics-Notes.txt"
This text file goes over the very initial setup you would go through to get most of the required pieces onto your computer in order to do ATAC-seq data analysis
Notably, my approach is to do everything in Terminal and R, with some utilization of public websites and apps for visualization, such as SeqPlots


"ATAC-seq-Descriptive-Workflow.txt"
This is a fairly thorough outline of most of the steps (I say most because it's always possible that I may have missed describing one or two steps) involved with collecting sequencing reads, aligning, filtering for quality, visualizing, calling THSs, evaluating levels, annotating nearest genes...
This is the first file that really starts to go into the code used in the ATAC-seq analysis, but in a meta way, talking about the steps instead of just showing you what was done


"Annotating-Genomic-Locations-Of-Coordainte-File.txt"
I included this file as a standout because it is a very useful few steps to annotate coordinates when PAVIS (website) is either down or does not contain the genome annotation for your species of interest


"Example-R-Script-For-Evaluating-Counts-In-THSs.txt"
In order to analyze the level of accessibility between a set of THSs you can use this sample script as a starting off point
For this study, we did a 1:1 comparison (control to submergence), however I have since applied DESeq2 for multi level comparisons (several different cell types) and have used clustering to separate out and identify THSs that are most accessible for one specific level even if there are 5 or more present
Because this was not the case I did not include this script and workflow, but if you are interested in more complicated evaluations such as this please reach out to me and I can send you such information, or I will upload it to GitHub once that work is finished or I have rewritten the script not to be too indicative of ongoing work


In the subfolder "Actual-Ran-Code" there is a comprehensive collection of all of the different commands ran in all the species if it pertains to ATAC-seq data analysis
Most notably, the most recently worked on "LOG-THSs-LOG.txt" file has a significant (~18k) amount of lines that makes up the latter parts of the analysis, such as localizing motif sequencing and annotating where they are relative to TSS, which THSs have these motifs, and which genes have a motif found in their promoter and if this motif is accessible or not
Additionally, the "Actual-Ran-Code" subfolder contains the R scripts used to evaluate the count levels in the THSs of each species





