# Prediction-of-differential-activity-of-transcription-factors
# DESCRIPTION
Set of R function for prediction of differentially active transcription factors. The presented method analyzes only TF-TG interactions for which it is known whether this interaction is inhibitory or activating. It is designed to compare the state of target genes associated with some TF between two groups of samples representing two different conditions. It is applicable in global mode - where the average expression of TGs is compared between two conditions and in single-pair mode, where TGs expression is compared between two paired samples. The activity score ATF is computed for every TF and its statistical significance estimated by bootstrap technique. Definition of ATF is based on the idea that in case of increased activation of TF in one group of samples we can expect increased expression of its positively regulated TGs and decreased expression of its negatively regulated TGs within the mentioned sample group in comparison with the second group. In case of decreased TF activity, we expect exactly the opposite behavior of TGs. The precise definition of ATF is presented within the publication XXX.

# USAGE
1. source file predictActiveTFs.R into the R workspace.
2. run the following function: results=predictActiveTFs(#grn, exprTable, nRand, isGlobal, nThreads), with all input variables properly set.
3. grn - 
