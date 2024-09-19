reformat<-function(validTable)
{
  validTable=as.data.frame(validTable)
  cols_char=c("genes_posT_pos", "genes_negT_neg", "genes_posT_neg", "genes_negT_pos",
              "hugo", "samples")
  colNames=colnames(validTable)
  for (i in 1: length(colNames))
      {
       if (colNames[i] %in% cols_char)
          {
           validTable[[i]]=as.character(validTable[[i]])
          }
       else
          {
           validTable[[i]]=as.numeric(validTable[[i]])
          }
      }
  return(validTable)
}  

correctPvals<-function(validationTable, isGlobal=T)
{
  library(qvalue)  
  
  get2tailedP<-function(p1,p2)
  {
   twoTiledP=2*pmin(validationTable[[p1]], validationTable[[p2]])
   twoTiledP[twoTiledP>1]=1
   return(twoTiledP)
  }
      
  getCorrP<-function(corrPval, Lsample)
  {
   L=Lsample & (!is.na(validationTable$pval_fractConsist)) 
   corrPval$pval_fract_consist_q[L]=qvalue(validationTable$pval_fractConsist[L])$qvalues
   L=Lsample & (!is.na(validationTable$pval_fractInconsist)) 
   corrPval$pval_fract_inconsist_q[L]=qvalue(validationTable$pval_fractInconsist[L])$qvalues
   
   L=Lsample & (!is.na(corrPval$pval_consist_inconsist)) 
   corrPval$pval_consist_inconsist_q[L]=qvalue(corrPval$pval_consist_inconsist[L])$qvalues
   L=Lsample & (!is.na(corrPval$pval_mixed)) 
   corrPval$pval_mixed_q[L]=qvalue(corrPval$pval_mixed[L])$qvalues
   L=Lsample & (!is.na(corrPval$pval_diffMean))
   corrPval$pval_diffMean_q[L]=qvalue(corrPval$pval_diffMean[L])$qvalues
   return(corrPval)
  }
      
  validationTable=reformat(validationTable)
  
  L=validationTable$nExprTargets_pos==0 & validationTable$nExprTargets_neg==0
  L_pval=grepl("pval_",colnames(validationTable))
  validationTable[L,c(L_pval)]=NA
  
  nRows=nrow(validationTable)
  nCols=ncol(validationTable)
  vect=numeric(nRows)
  corrPval=data.frame("pval_fract_consist_q"=vect,
                      "pval_fract_inconsist_q"=vect,
                      "pval_consist_inconsist"=vect,
                      "pval_consist_inconsist_q"=vect,
                      "pval_mixed"=vect,
                      "pval_mixed_q"=vect,
                      "pval_diffMean"=vect,
                      "pval_diffMean_q"=vect)
  corrPval$pval_consist_inconsist=get2tailedP("pval_consist_inconsist_higher","pval_consist_inconsist_lower")
  corrPval$pval_mixed=get2tailedP("pval_mixed_higher","pval_mixed_lower")
  corrPval$pval_diffMean=get2tailedP("pval_diffMean_higher","pval_diffMean_lower")
  
  if (isGlobal==F)
     {
      samples=unique(validationTable$samples)
      for (sample in samples)
          {
           L_sample=validationTable$samples==sample
           corrPval=getCorrP(corrPval, L_sample)
          }
     }
  else
     {
      L_sample=rep(T,nRows)
      corrPval=getCorrP(corrPval, L_sample)
     }
  
  validationTable=cbind(validationTable, corrPval)
  
  ## replace corr pval to NA when no targets gene expressed
  validationTable[L,c((nCols+1):ncol(validationTable))]=NA
  
  return(validationTable)  
}

getMutTable_genesOnly<-function(TFs_unique, samples_T_unique)
{
  nTFs=length(TFs_unique)
  nSamples=length(samples_T_unique)
  nRows=nTFs*nSamples
  mutTable_genesOnly=data.frame("hugo"=character(nRows),
                                "samples"=character(nRows), stringsAsFactors=F) 
  start=1
  for (i in 1:nSamples)
  {
    mutTable_genesOnly$hugo[start:(start+nTFs-1)]=TFs_unique
    mutTable_genesOnly$samples[start:(start+nTFs-1)]=samples_T_unique[i]
    start=start+nTFs
  }
  return(mutTable_genesOnly)
}

getStatistics<-function(diffExpr_pos, diffExpr_neg)
{
  nExprTargets_pos=sum(!is.na(diffExpr_pos))
  nExprTargets_neg=sum(!is.na(diffExpr_neg))
  nExprTargets=nExprTargets_pos+nExprTargets_neg
  
  ## statistics for pvalue computation 
  if (nExprTargets>0)
    {
     L_pos_pos=diffExpr_pos>0 # positive TGs with positive expression difference
     L_neg_neg=diffExpr_neg<0 # negative TGs with negative expression difference
     L_pos_pos[is.na(L_pos_pos)]=F
     L_neg_neg[is.na(L_neg_neg)]=F
    
     L_pos_neg=diffExpr_pos<0 # positive TGs with negative expression difference
     L_neg_pos=diffExpr_neg>0 # negative TGS with positive expression difference
     L_pos_neg[is.na(L_pos_neg)]=F
     L_neg_pos[is.na(L_neg_pos)]=F
     
     nConsist_pos=sum(L_pos_pos)
     nConsist_neg=sum(L_neg_neg)
     nInconsist_pos=sum(L_pos_neg)
     nInconsist_neg=sum(L_neg_pos)
    
     Pp=diffExpr_pos[L_pos_pos] # expression values of positive TGs with positive expression difference
     Pn=diffExpr_pos[L_pos_neg] # expression values of positive TGs with negative expression difference
     Nn=diffExpr_neg[L_neg_neg] # expression values of negative TGs with negative expression difference
     Np=diffExpr_neg[L_neg_pos] # expression values of negative TGs with positive expression difference
    
     meanConsist=mean(c(-Nn,Pp))
     meanNotConsist=mean(c(-Pn,Np))
     if (length(meanConsist)==0 | sum(is.na(meanConsist)==F)==0)
        {
         meanConsist=0
        }
     if (length(meanNotConsist)==0 | sum(is.na(meanNotConsist)==F)==0)
        {
         meanNotConsist=0
        }
   
     genes_posT_pos=names(diffExpr_pos)[L_pos_pos]
     genes_negT_neg=names(diffExpr_neg)[L_neg_neg]
     genes_posT_neg=names(diffExpr_pos)[L_pos_neg]
     genes_negT_pos=names(diffExpr_neg)[L_neg_pos]
      
     sumDiff_posUp_negDown=sum(c(-Nn,-Np,Pp,Pn)) ## consist scenario...
     sumDiff_posDown_negUp=-sumDiff_posUp_negDown ## inconsist scenario... =sum(c(Nn,Np,-Pp,-Pn))
     sumDiff_posUp_negUp=sum(c(Nn,Np,Pp,Pn)) ## consistent effect on positive tgs and inconsistent effect on negative tgs 
     sumDiff_posDown_negDown=-sumDiff_posUp_negUp ## inconsistent effect on positive tgs and consistent effect on negative tgs... =sum(c(-Nn,-Np,-Pp,-Pn))
      
     ## test if consistent differences are significantly higher than inconsistent and vice versa
     diffMean=meanConsist-meanNotConsist
  
     ## fraction of consistent and inconsistent differences
     fract_consist=(nConsist_pos+nConsist_neg)/nExprTargets
     fract_inconsist=(nInconsist_pos+nInconsist_neg)/nExprTargets
    
     ## compute TF activity (ATF)
     meanEffect_consist_inconsist=sumDiff_posUp_negDown/nExprTargets
     meanEffect_mixed=sumDiff_posUp_negUp/nExprTargets
    }
  else
  {
    nConsist_pos=NA
    nConsist_neg=NA
    nInconsist_pos=NA
    nInconsist_neg=NA
    fract_consist=NA
    fract_inconsist=NA
    meanEffect_consist_inconsist=NA
    meanEffect_mixed=NA
    diffMean=NA
    genes_posT_pos=NA
    genes_negT_neg=NA
    genes_posT_neg=NA
    genes_negT_pos=NA
  }
  
  stat=list("nExprTargets"=nExprTargets,"nExprTargets_pos"=nExprTargets_pos, "nExprTargets_neg"=nExprTargets_neg,
            "nConsist_pos"=nConsist_pos, "nConsist_neg"=nConsist_neg, "nInconsist_pos"=nInconsist_pos, "nInconsist_neg"=nInconsist_neg,
            "fract_consist"=fract_consist, "fract_inconsist"=fract_inconsist,
            "meanEffect_consist_inconsist"=meanEffect_consist_inconsist, "meanEffect_mixed"=meanEffect_mixed,
            "diffMean"=diffMean,
            "genes_posT_pos"=genes_posT_pos, "genes_negT_neg"=genes_negT_neg,
            "genes_posT_neg"=genes_posT_neg, "genes_negT_pos"=genes_negT_pos)
  return(stat)
}

permuteTest_regEffect<-function(diffExpr, stat, nTargets_pos, nTargets_neg, nRand)
{
  n_higher_fractConsist=0
  n_higher_fractInconsist=0
  n_higher_consist_inconsist=0
  n_lower_consist_inconsist=0
  n_higher_mixed=0
  n_lower_mixed=0
  n_higher_diff=0
  n_lower_diff=0
  
  nTargets=nTargets_pos+nTargets_neg
  diffExprPositive=diffExpr[diffExpr>0]
  diffExprNegative=diffExpr[diffExpr<0]
  
  for (i in 1:nRand)
  {
    ## get statistics on random sample #########################################
    diffExpr_pos_rand=sample(diffExpr, nTargets_pos, replace=TRUE)
    diffExpr_neg_rand=sample(diffExpr, nTargets_neg, replace=TRUE)
    stat_rand=getStatistics(diffExpr_pos_rand, diffExpr_neg_rand)
    
    ## more strict random model ################################################
    diffExpr_pos_rand2=c(sample(diffExprPositive, stat$nConsist_pos, replace=TRUE), sample(diffExprNegative, (nTargets_pos-stat$nConsist_pos), replace=TRUE))
    diffExpr_neg_rand2=c(sample(diffExprNegative, stat$nConsist_neg, replace=TRUE), sample(diffExprPositive, (nTargets_neg-stat$nConsist_neg), replace=TRUE))
    stat_rand2=getStatistics(diffExpr_pos_rand2, diffExpr_neg_rand2)
    
    ## fract stat ##
    if (stat_rand$fract_consist>=stat$fract_consist)
       {
        n_higher_fractConsist=n_higher_fractConsist+1
       }
    if (stat_rand$fract_inconsist>=stat$fract_inconsist)
       {
        n_higher_fractInconsist=n_higher_fractInconsist+1
       }
    
    ## TF activity (ATF) stats ##
    if (stat_rand$meanEffect_consist_inconsist>=stat$meanEffect_consist_inconsist)
       {
        n_higher_consist_inconsist=n_higher_consist_inconsist+1
       }
    if (stat_rand$meanEffect_consist_inconsist<=stat$meanEffect_consist_inconsist)
       {
        n_lower_consist_inconsist=n_lower_consist_inconsist+1
       }
    
    ## TF activity mixed stat ##
    if (stat_rand$meanEffect_mixed>=stat$meanEffect_mixed)
       {
        n_higher_mixed=n_higher_mixed+1
       }
    if (stat_rand$meanEffect_mixed<=stat$meanEffect_mixed)
       {
        n_lower_mixed=n_lower_mixed+1
       }
    
    ## diff means stat ##
    if (stat_rand2$diffMean>=stat$diffMean)
       {
        n_higher_diff=n_higher_diff+1
       }
    if (stat_rand2$diffMean<=stat$diffMean)
       {
        n_lower_diff=n_lower_diff+1
       }
  }
  
  pval_fractConsist=max(n_higher_fractConsist,1)/nRand
  pval_fractInconsist=max(n_higher_fractInconsist,1)/nRand
  
  pval_consist_inconsist_higher=max(n_higher_consist_inconsist,1)/nRand
  pval_consist_inconsist_lower=max(n_lower_consist_inconsist,1)/nRand
  
  pval_mixed_higher=max(n_higher_mixed,1)/nRand
  pval_mixed_lower=max(n_lower_mixed,1)/nRand
  
  pval_diffMean_higher=max(n_higher_diff,1)/nRand
  pval_diffMean_lower=max(n_lower_diff,1)/nRand

  return(list("pval_fractConsist"=pval_fractConsist, 
              "pval_fractInconsist"=pval_fractInconsist, 
              "pval_consist_inconsist_higher"=pval_consist_inconsist_higher, 
              "pval_consist_inconsist_lower"=pval_consist_inconsist_lower, 
              "pval_mixed_higher"=pval_mixed_higher,
              "pval_mixed_lower"=pval_mixed_lower,
              "pval_diffMean_higher"=pval_diffMean_higher,
              "pval_diffMean_lower"=pval_diffMean_lower))
}

computeStat_tfSample<-function(validationTable_regEffect, i, exprTable, L_tg_pos, L_tg_neg, nRand, isImputeMissing, isGlobal)
{
  if (isGlobal==F)
     {
      ## compute difference in expresion between two (T and N) paired samples
      sampleID_T=validationTable_regEffect$samples[i]
      sampleID_N=sub("T","N",sampleID_T)
      expr_T=exprTable[,sampleID_T]
      expr_N=exprTable[,sampleID_N]
      if (isImputeMissing==T)
         {
          expr_T[is.na(expr_T)]=0
          expr_N[is.na(expr_N)]=0
         }
      diffExpr=expr_T-expr_N
     }
  else
     {
      ######## computation based on FOLDCHANGE solely 
      diffExpr=exprTable$logFC
      names(diffExpr)=exprTable$gene
      diffExpr_pos=diffExpr[L_tg_pos]
      diffExpr_neg=diffExpr[L_tg_neg]
     }
 
  ## compute statistic #########################################################
  diffExpr_pos=diffExpr[L_tg_pos]
  diffExpr_neg=diffExpr[L_tg_neg]
  stat=getStatistics(diffExpr_pos, diffExpr_neg)
  
  if (stat$nExprTargets>0)
     {
      diffExpr=diffExpr[!is.na(diffExpr)]
      n_pos=stat$nExprTargets_pos
      n_neg=stat$nExprTargets_neg  
   
    
      ## compute pvalues by permutation simulation
      ### ORIGINAL VERSION - COMPUTATIONALLY INTESIVE ####
      pvals=permuteTest_regEffect(diffExpr, stat, n_pos, n_neg, nRand)
    
      ## fill results variable ###################################################
      tfName=validationTable_regEffect$hugo[i]
      validationTable_regEffect$tfDiffExpr[i]=diffExpr[tfName]
      validationTable_regEffect$nExprTargets_pos[i]=stat$nExprTargets_pos
      validationTable_regEffect$nExprTargets_neg[i]=stat$nExprTargets_neg
      if (stat$nExprTargets_pos>0)
         {
          validationTable_regEffect$fract_posT_pos[i]=length(stat$genes_posT_pos)/stat$nExprTargets_pos
          validationTable_regEffect$fract_posT_neg[i]=length(stat$genes_posT_neg)/stat$nExprTargets_pos
         }  
      if (stat$nExprTargets_neg>0)
         {
          validationTable_regEffect$fract_negT_neg[i]=length(stat$genes_negT_neg)/stat$nExprTargets_neg
          validationTable_regEffect$fract_negT_pos[i]=length(stat$genes_negT_pos)/stat$nExprTargets_neg
         }
      validationTable_regEffect$genes_posT_pos[i]=paste(stat$genes_posT_pos, collapse = ",")
      validationTable_regEffect$genes_negT_neg[i]=paste(stat$genes_negT_neg, collapse = ",")
      validationTable_regEffect$genes_posT_neg[i]=paste(stat$genes_posT_neg, collapse = ",")
      validationTable_regEffect$genes_negT_pos[i]=paste(stat$genes_negT_pos, collapse = ",")
      
      validationTable_regEffect$meanEffect_consist_inconsist[i]=stat$meanEffect_consist_inconsist
      validationTable_regEffect$meanEffect_mixed[i]=stat$meanEffect_mixed
      validationTable_regEffect$diffMean[i]=stat$diffMean
      
      validationTable_regEffect$pval_fractConsist[i]=pvals$pval_fractConsist
      validationTable_regEffect$pval_fractInconsist[i]=pvals$pval_fractInconsist
      validationTable_regEffect$pval_consist_inconsist_higher[i]=pvals$pval_consist_inconsist_higher
      validationTable_regEffect$pval_consist_inconsist_lower[i]=pvals$pval_consist_inconsist_lower
      validationTable_regEffect$pval_mixed_higher[i]=pvals$pval_mixed_higher
      validationTable_regEffect$pval_mixed_lower[i]=pvals$pval_mixed_lower
      validationTable_regEffect$pval_diffMean_higher[i]=pvals$pval_diffMean_higher
      validationTable_regEffect$pval_diffMean_lower[i]=pvals$pval_diffMean_lower
  }
  
  return(validationTable_regEffect)
}

predictActiveTFs<-function(grn, exprTable, nRand, isGlobal, nThreads)
{
  computeStat<-function(indexes, validationTable_regEffect, exprTable, L_tg_pos, L_tg_neg, nRand)
  {
    ## get subset of validationTable
    validationTable_regEffect=validationTable_regEffect[indexes,]
    nRows=nrow(validationTable_regEffect)
    for (i in 1:nRows)
    {print(i)
      if (validationTable_regEffect$nTargets_pos[i]+validationTable_regEffect$nTargets_neg[i]>0)
      {
        ## paired test
        tfName=validationTable_regEffect$hugo[i]
        validationTable_regEffect=computeStat_tfSample(validationTable_regEffect, 
                                                       i, 
                                                       exprTable, 
                                                       L_tg_pos[,tfName], 
                                                       L_tg_neg[,tfName], 
                                                       nRand, 
                                                       isImputeMissing=F,
                                                       isGlobal=isGlobal)
      }
    }
    return(validationTable_regEffect)
  }
  
  tfs_unique=unique(grn[,1])  
  if (isGlobal==F)
     {
      ## generate TF-sample table
      L_T=grepl("T", colnames(exprTable), fixed=T)
      samples_T_unique=colnames(exprTable)[L_T]
      tfSampleTable=getMutTable_genesOnly(tfs_unique, samples_T_unique)
      geneNames=rownames(exprTable)
     }
  else
     {
      ## get unique TFs
      tfSampleTable=as.data.frame(tfs_unique)
      colnames(tfSampleTable)="hugo"
      geneNames=exprTable$gene
     }
  
  ## generate results variable #################################################
  nRows=nrow(tfSampleTable)
  evn=numeric(nRows)
  evch=character(nRows)
  validationTable_regEffect=data.frame("tfDiffExpr"=evn,"allelicFraction"=evn, "nTargets_pos"=evn, "nExprTargets_pos"=evn, "nTargets_neg"=evn, "nExprTargets_neg"=evn,
                                       "fract_posT_pos"=evn, "fract_negT_neg"=evn, "fract_posT_neg"=evn, "fract_negT_pos"=evn, 
                                       "genes_posT_pos"=evch,"genes_negT_neg"=evch, "genes_posT_neg"=evch,"genes_negT_pos"=evch,
                                       "meanEffect_consist_inconsist"=evn, "meanEffect_mixed"=evn, "diffMean"=evn,
                                       "pval_fractConsist"=evn,"pval_fractInconsist"=evn,
                                       "pval_consist_inconsist_higher"=evn,"pval_consist_inconsist_lower"=evn,
                                       "pval_mixed_higher"=evn,"pval_mixed_lower"=evn,
                                       "pval_diffMean_higher"=evn,"pval_diffMean_lower"=evn)
                                       
  validationTable_regEffect=cbind(tfSampleTable, validationTable_regEffect) 
  
  ## precompute output variable ###
  nTfs=length(tfs_unique)
  L_tg_pos=array(FALSE,c(nrow(exprTable),nTfs))
  L_tg_neg=array(FALSE,c(nrow(exprTable),nTfs))
  colnames(L_tg_pos)=tfs_unique
  colnames(L_tg_neg)=tfs_unique
  for (i in 1:nTfs)
      {
       tfName=validationTable_regEffect$hugo[i]
       L_tg_grn_pos=grn[,1]==tfName & grn[,3]=="Activation"
       L_tg_grn_neg=grn[,1]==tfName & grn[,3]=="Inhibition"
       targetGenes_pos=grn[L_tg_grn_pos,2]
       targetGenes_neg=grn[L_tg_grn_neg,2]
       nTargets_pos=length(targetGenes_pos)
       nTargets_neg=length(targetGenes_neg)
       L=validationTable_regEffect$hugo==tfName
       validationTable_regEffect$nTargets_pos[L]=nTargets_pos
       validationTable_regEffect$nTargets_neg[L]=nTargets_neg
       L_tg_pos[,i]=geneNames %in% targetGenes_pos
       L_tg_neg[,i]=geneNames %in% targetGenes_neg
      }
  
 
  ## parallel computing ########################################################
  library(parallel)
  nCores=nThreads
  cl=makeCluster(nCores)
  clusterExport(cl,
                c("computeStat_tfSample",
                  "permuteTest_regEffect",
                  "getStatistics"),
                envir = .GlobalEnv)
  indexesAll=seq(1:nrow(validationTable_regEffect))
  
  ## NOT PARALELISED VERSION --FOR DEBUGING
  #validationTable_regEffect=computeStat(indexesAll, validationTable_regEffect, exprTable, L_tg_pos, L_tg_neg, nRand)
  
  validationTable_regEffect=parSapply(cl,indexesAll,computeStat, validationTable_regEffect, exprTable, L_tg_pos, L_tg_neg, nRand)
  validationTable_regEffect=t(validationTable_regEffect)
  stopCluster(cl)
  
  ## pvalues correction
  validationTable_regEffect=correctPvals(validationTable_regEffect, isGlobal=isGlobal)
  
  return(validationTable_regEffect)
} 
