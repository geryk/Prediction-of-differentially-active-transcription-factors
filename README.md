# Prediction-of-differential-activity-of-transcription-factors
# DESCRIPTION
Set of R function for prediction of differentially active transcription factors (TFs). The presented method analyzes only TF-TG interactions for which it is known whether this interaction is inhibitory or activating. It is designed to compare the state of target genes (TGs) associated with some TF between two groups of samples representing two different conditions. It is applicable in global mode - where the average expression of TGs is compared between two conditions and in single-pair mode, where TGs expression is compared between two paired samples. The activity score ATF is computed for every TF and its statistical significance estimated by bootstrap technique. Definition of ATF is based on the idea that in case of increased activation of TF in one group of samples we can expect increased expression of its positively regulated TGs and decreased expression of its negatively regulated TGs within the mentioned sample group in comparison with the second group. In case of decreased TF activity, we expect exactly the opposite behavior of TGs. The precise definition of ATF is presented within the publication XXX.

# USAGE
1. source file predictActiveTFs.R into the R workspace.
   
2. optionally load inputs.RData into the R workspace. The data contains  gene regulatory network (grn_Biochimie_regEffect) and 2 topTables (topTable_paired_proteomics and topTable_rnaSeq_RIN7)
   
3. run the following function: results=predictActiveTFs(**grn**, **exprTable**, **nRand**, **isGlobal**, **nThreads**), with all input variables properly set.

**grn** - gene regulatory network in the form of matrix or data.frame, where 1. st column contain names of TFs, 2. st column contain theirs TGs  and
3. st column contain "Activation" if TF increases expression of the TG from the same row or "Inhibition" if the TF decreases its expression.

**exprTable** - 2 possibilities: a) top table produced by limma package, where 1.st column, (named "gene") contains genes and 2.st column (named "logFC") contains logFC values
OR b) expression matrix with normalised expression values, where rows represents genes and columns represents samples. The row names contain gene names and column names contains sample names.
It is expected that expression matrix contains exclusively paired samples, corresponding to two conditions. 
One sample from each pair must start with "N" and the other one with "T", the remaining part of sample name is optional but must be same for paired samples.
For example: T01 and N01 represents correct names for paired samples. 

**nRand** - number of permutations to obtain random estimate of TF activities (ATF).

**isGlobal** - TRUE if **exprTable** is top table, FASLE otherwise

**nThreads** - number of threads for parallelisation
