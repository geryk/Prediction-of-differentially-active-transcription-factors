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
  
  getCorrP<-function(corrPval, L)
  {
   corrPval$pval_fract_BH[L]=p.adjust(validationTable$pval_fract[L], method="BH")
   corrPval$pval_mean_consist_BH[L]=p.adjust(validationTable$pval_mean_consist[L], method="BH")
   corrPval$pval_mean_inconsist_BH[L]=p.adjust(validationTable$pval_mean_inconsist[L], method="BH")
   corrPval$pval_meanRatio_BH[L]=p.adjust(validationTable$pval_meanRatio[L], method="BH")
   corrPval$pval_mean_consist_inconsist_q[L]=qvalue(corrPval$pval_mean_consist_inconsist[L])$qvalues
   return(corrPval)
  }
      
  validationTable=reformat(validationTable)
  
  L=validationTable$nExprTargets_pos==0 & validationTable$nExprTargets_neg==0
  L_pval=grepl("pval_",colnames(validationTable))
  validationTable[L,c(L_pval)]=NA
  
  nRows=nrow(validationTable)
  nCols=ncol(validationTable)
  vect=numeric(nRows)
  corrPval=data.frame("pval_fract_BH"=vect,
                      "pval_mean_consist_BH"=vect,
                      "pval_mean_inconsist_BH"=vect,
                      "pval_meanRatio_BH"=vect,
                      "pval_mean_consist_inconsist"=vect,
                      "pval_mean_consist_inconsist_q"=vect)
  
  corrPval$pval_mean_consist_inconsist=2*pmin(validationTable$pval_mean_consist, validationTable$pval_mean_inconsist)
  corrPval$pval_mean_consist_inconsist[corrPval$pval_mean_consist_inconsist>1]=1
  
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
      Ltrue=rep(T,nRows)
      corrPval=getCorrP(corrPval, Ltrue)
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

getStatistics<-function(diffExpr_pos, diffExpr_neg, isConsist)
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
    
    Pp=diffExpr_pos[L_pos_pos]
    Pn=diffExpr_pos[L_pos_neg]
    Nn=diffExpr_neg[L_neg_neg]
    Np=diffExpr_neg[L_neg_pos]
    
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
    
    if (isConsist==T)
    {
      nConsist_pos=sum(L_pos_pos)
      nConsist_neg=sum(L_neg_neg)
      genes_posT_pos=names(diffExpr_pos)[L_pos_pos]
      genes_negT_neg=names(diffExpr_neg)[L_neg_neg]
      genes_posT_neg=names(diffExpr_pos)[L_pos_neg]
      genes_negT_pos=names(diffExpr_neg)[L_neg_pos]
      #sumNegDiff=sum((-1*diffExpr_neg),na.rm=T)
      #sumPosDiff=sum(diffExpr_pos, na.rm=T)
      
      #sumDiff_pos=sum(c(Pp,Pn))
      #sumDiff_neg=sum(c(-Nn,-Np))
      
      sumDiff_posUp_negDown=sum(c(-Nn,-Np,Pp,Pn)) ## consist scenario...
      sumDiff_posUp_negUp=sum(c(Nn,Np,Pp,Pn)) ## consistent effect on positive tgs and inconsistent effect on negative tgs 
      sumDiff_posDown_negDown=sum(c(-Nn,-Np,-Pp,-Pn)) ## inconsistent effect on positive tgs and consistent effect on negative tgs
      sumDiff_posDown_negUp=sum(c(Nn,Np,-Pp,-Pn)) ## inconsist scenario...
      
      meanRatio=meanConsist-meanNotConsist
    }
    else if (isConsist==F) ## DEPRECATED - NOT USED
    {
      nConsist_pos=sum(L_pos_neg)
      nConsist_neg=sum(L_neg_pos)
      genes_consist_pos=names(diffExpr_pos)[L_pos_neg]
      genes_consist_neg=names(diffExpr_neg)[L_neg_pos]
      #sumNegDiff=sum(diffExpr_neg,na.rm=T)
      #sumPosDiff=sum((-1*diffExpr_pos), na.rm=T)
      sumDiff=sum(c(Nn,Np,-Pp,-Pn))
      sumDiff_pos=sum(c(-Pp,-Pn))
      sumDiff_neg=sum(c(Nn,Np))
      meanRatio=meanNotConsist-meanConsist
    }
    fract_consist=(nConsist_pos+nConsist_neg)/nExprTargets
    
    meanEffect_posUp_negDown=sumDiff_posUp_negDown/nExprTargets
    meanEffect_posUp_negUp=sumDiff_posUp_negUp/nExprTargets
    meanEffect_posDown_negDown=sumDiff_posDown_negDown/nExprTargets
    meanEffect_posDown_negUp=sumDiff_posDown_negUp/nExprTargets
    
    #meanEffectDiff_pos=sumDiff_pos/nExprTargets_pos
    #meanEffectDiff_neg=sumDiff_neg/nExprTargets_neg
  }
  else
  {
    nConsist_pos=NA
    nConsist_neg=NA
    sumNegDiff=NA
    sumPosDiff=NA
    fract_consist=NA
    meanEffect_posUp_negDown=NA
    meanEffect_posUp_negUp=NA
    meanEffect_posDown_negDown=NA
    meanEffect_posDown_negUp=NA
    meanRatio=NA
    genes_posT_pos=NA
    genes_negT_neg=NA
    genes_posT_neg=NA
    genes_negT_pos=NA
  }
  
  stat=list("nExprTargets"=nExprTargets,"nExprTargets_pos"=nExprTargets_pos, "nExprTargets_neg"=nExprTargets_neg,
            "nConsist_pos"=nConsist_pos, "nConsist_neg"=nConsist_neg, 
            "fract_consist"=fract_consist, 
            "meanEffect_consist"=meanEffect_posUp_negDown, 
            "meanEffect_posUp_negUp"=meanEffect_posUp_negUp,
            "meanEffect_posDown_negDown"=meanEffect_posDown_negDown,
            "meanEffect_inconsist"=meanEffect_posDown_negUp,
            "meanRatio"=meanRatio,
            "genes_posT_pos"=genes_posT_pos, "genes_negT_neg"=genes_negT_neg,
            "genes_posT_neg"=genes_posT_neg, "genes_negT_pos"=genes_negT_pos)
  return(stat)
}

permuteTest_regEffect<-function(diffExpr, stat, nTargets_pos, nTargets_neg, nRand, isConsist)
{
  n_higher_fractConsist=0
  n_higher_mean_consist=0
  n_higher_mean_posUp_negUp=0
  n_higher_mean_posDown_negDown=0
  n_higher_mean_inconsist=0
  n_higher_ratio=0
  nTargets=nTargets_pos+nTargets_neg
  
  if (isConsist==T)
     {
      diffExprPositive=diffExpr[diffExpr>0]
      diffExprNegative=diffExpr[diffExpr<0]
     }
  else
     {
      diffExprPositive=diffExpr[diffExpr<0]
      diffExprNegative=diffExpr[diffExpr>0]
     }
  
  print("nTargets pos/neg:")
  print(nTargets_pos)
  print(nTargets_neg)
  
  for (i in 1:nRand)
  {
    ## get statistics on random sample 
    diffExpr_pos_rand=sample(diffExpr, nTargets_pos, replace=TRUE)
    diffExpr_neg_rand=sample(diffExpr, nTargets_neg, replace=TRUE)
    stat_rand=getStatistics(diffExpr_pos_rand, diffExpr_neg_rand, isConsist)
    
    ## more strict random model
    diffExpr_pos_rand2=c(sample(diffExprPositive, stat$nConsist_pos, replace=TRUE), sample(diffExprNegative, (nTargets_pos-stat$nConsist_pos), replace=TRUE))
    diffExpr_neg_rand2=c(sample(diffExprNegative, stat$nConsist_neg, replace=TRUE), sample(diffExprPositive, (nTargets_neg-stat$nConsist_neg), replace=TRUE))
    stat_rand2=getStatistics(diffExpr_pos_rand2, diffExpr_neg_rand2, isConsist)
    
    ## fract stat
    if (stat_rand$fract_consist>=stat$fract_consist)
    {
      n_higher_fractConsist=n_higher_fractConsist+1
    }
    
    ## mean stats
    if (stat_rand$meanEffect_consist>=stat$meanEffect_consist)
    {
      n_higher_mean_consist=n_higher_mean_consist+1
    }
    if (stat_rand$meanEffect_posUp_negUp>=stat$meanEffect_posUp_negUp)
    {
      n_higher_mean_posUp_negUp=n_higher_mean_posUp_negUp+1
    }
    if (stat_rand$meanEffect_posDown_negDown>=stat$meanEffect_posDown_negDown)
    {
      n_higher_mean_posDown_negDown=n_higher_mean_posDown_negDown+1
    }
    if (stat_rand$meanEffect_inconsist>=stat$meanEffect_inconsist)
    {
      n_higher_mean_inconsist=n_higher_mean_inconsist+1
    }
    
    ## ratio stat
    if (stat_rand2$meanRatio>=stat$meanRatio)
    {
      n_higher_ratio=n_higher_ratio+1
    }
  }
  
  pval_fract=max(n_higher_fractConsist,1)/nRand
  
  pval_mean_consist=max(n_higher_mean_consist,1)/nRand
  pval_mean_posUp_negUp=max(n_higher_mean_posUp_negUp,1)/nRand
  pval_mean_posDown_negDown=max(n_higher_mean_posDown_negDown,1)/nRand
  pval_mean_inconsist=max(n_higher_mean_inconsist,1)/nRand
  
  pval_meanRatio=max(n_higher_ratio,1)/nRand
  
  return(list("pval_fract"=pval_fract, 
              "pval_mean_consist"=pval_mean_consist, 
              "pval_mean_posUp_negUp"=pval_mean_posUp_negUp, 
              "pval_mean_posDown_negDown"=pval_mean_posDown_negDown, 
              "pval_mean_inconsist"=pval_mean_inconsist,
              "pval_meanRatio"=pval_meanRatio))
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
      ######## computation basen on FOLDCHANGE solely 
      diffExpr=exprTable$logFC
      names(diffExpr)=exprTable$gene
      diffExpr_pos=diffExpr[L_tg_pos]
      diffExpr_neg=diffExpr[L_tg_neg]
     }
 
  ## compute statistic #########################################################
  diffExpr_pos=diffExpr[L_tg_pos]
  diffExpr_neg=diffExpr[L_tg_neg]
  stat=getStatistics(diffExpr_pos, diffExpr_neg, isConsist=T)
  
  if (stat$nExprTargets>0)
     {
      diffExpr=diffExpr[!is.na(diffExpr)]
      n_pos=stat$nExprTargets_pos
      n_neg=stat$nExprTargets_neg  
   
    
      ## compute pvalues by permutation simulation
      ### ORIGINAL VERSION - COMPUTATIONALLY INTESIVE ####
      pvals=permuteTest_regEffect(diffExpr, stat, n_pos, n_neg, nRand, isConsist=T)
    
      ## fill results variable ###################################################
      tfName=validationTable_regEffect$hugo[i]
      validationTable_regEffect$tfDiffExpr[i]=diffExpr[tfName]
      validationTable_regEffect$nExprTargets_pos[i]=stat$nExprTargets_pos
      validationTable_regEffect$nExprTargets_neg[i]=stat$nExprTargets_neg
      validationTable_regEffect$meanEffect_consist[i]=stat$meanEffect_consist
      validationTable_regEffect$meanEffect_posUp_negUp[i]=stat$meanEffect_posUp_negUp
      validationTable_regEffect$meanEffect_posDown_negDown[i]=stat$meanEffect_posDown_negDown
      validationTable_regEffect$meanEffect_inconsist[i]=stat$meanEffect_inconsist
      validationTable_regEffect$pval_fract[i]=pvals$pval_fract
      validationTable_regEffect$pval_mean_consist[i]=pvals$pval_mean_consist
      validationTable_regEffect$pval_mean_posUp_negUp[i]=pvals$pval_mean_posUp_negUp
      validationTable_regEffect$pval_mean_posDown_negDown[i]=pvals$pval_mean_posDown_negDown
      validationTable_regEffect$pval_mean_inconsist[i]=pvals$pval_mean_inconsist
      validationTable_regEffect$pval_meanRatio[i]=pvals$pval_meanRatio
      validationTable_regEffect$genes_posT_pos[i]=paste(stat$genes_posT_pos, collapse = ",")
      validationTable_regEffect$genes_negT_neg[i]=paste(stat$genes_negT_neg, collapse = ",")
      validationTable_regEffect$genes_posT_neg[i]=paste(stat$genes_posT_neg, collapse = ",")
      validationTable_regEffect$genes_negT_pos[i]=paste(stat$genes_negT_pos, collapse = ",")
    
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
  }
  
  return(validationTable_regEffect)
}

predictActiveTFs<-function(grn, exprTable, nRand, isGlobal, nThreads)
{
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
                                       "meanEffect_consist"=evn, "meanEffect_posUp_negUp"=evn, "meanEffect_posDown_negDown"=evn,"meanEffect_inconsist"=evn,
                                       "pval_fract"=evn, 
                                       "pval_mean_consist"=evn, "pval_mean_posUp_negUp"=evn, "pval_mean_posDown_negDown"=evn, "pval_mean_inconsist"=evn, 
                                       "pval_meanRatio"=evn,
                                       "genes_posT_pos"=evch,"genes_negT_neg"=evch, "genes_posT_neg"=evch,"genes_negT_pos"=evch)
  validationTable_regEffect=cbind(tfSampleTable, validationTable_regEffect) 
  
  ## precompute inputs #########################################################
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
  #validationTable_regEffect_new=computeStat(indexesAll, validationTable_regEffect, exprTable, L_tg_pos, L_tg_neg, nRand)
  
  validationTable_regEffect_new=parSapply(cl,indexesAll,computeStat, validationTable_regEffect, exprTable, L_tg_pos, L_tg_neg, nRand)
  validationTable_regEffect_new=t(validationTable_regEffect_new)
  stopCluster(cl)
  
  return(validationTable_regEffect_new)
} 
