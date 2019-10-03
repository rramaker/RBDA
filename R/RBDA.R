#' Rank Based Differential Expression Analysis
#' This function identifies differentially expressed genes from a raw count table
#' @param countData A raw count table (dataframe or matrix) with samples as columns and genes as rows.
#' @param colData A dataframe with samples as rows and columns factor variables indicating experimental groups.
#' @param testVariable Character indicating the column name in colData that reresents the variable of interest for differential expression analysis.
#' @param batch_family_variable This paramater is used in N=1 cases only. Character that specifies a column name in colData that indicates other samples sequenced in the same batch as the sample of interest or related family members. It can be used as a downstream filter when comparing one sample to a large reference set. Defaults to NULL.
#' @param randomSeed A numeric indicating a random seed for reproducible analysis. Defaults to 1990
#' @param minP A numeric between 0 and one indicating the minimum possible P-value computed via random sampling from the countTable. The larger this number is the faster the compute time will be. Defaults to 0.000005.
#' @param numCores A numeric indicating the number of cores used for computing. Most laptops can readily handle up to 4. Defaults to 4.
#' @return Returns a data frame containing differential expression results. The first column, "variance_rank," indicates the percentile rank of a gene's rank variance. Low numbers indicate a gene exhibits very low variance across samples. The second column "test_statistic" is the RBDA test statistic. The third column "p_value" indicates the probability of observing the test statistic by chance alone. A fourth column, "min_batch_fam_p_value," will only result if the batch_family_variable parameter is specified and will provide the minimum p-value computed for any samples indicated by this factor variable in colData.
#' @export
#' @import parallel
#' @examples
#'data("BRCA_Counts")
#'colData<-data.frame(testVariable=as.factor(c(rep(0,5),rep(1,5))), row.names=colnames(BRCA_Counts))
#'RBDA(countData = BRCA_Counts,colData,testVariable = "testVariable", randomSeed = 1990, numCores=1)

RBDA<-function(countData, colData, testVariable, batch_family_variable=NULL, randomSeed=1990, minP=0.000005, numCores=4){
  set.seed(randomSeed)
  countRanks<-apply(countData, 2, function(x) rank(-x, ties.method = "random"))
  countRanks<-apply(countRanks,2,as.numeric)
  if(length(which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[2]))==1){
    countVar_list<-list(apply(countRanks[,which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[1])],1,function(x) sum(abs(x-mean(x))^2)),0)
  }
  if(length(which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[2]))>1){
    countVar_list<-lapply(1:length(levels(as.factor(colData[,testVariable]))), function(var) apply(countRanks[,which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[var])],1,function(x) sum(abs(x-mean(x))^2)))
  }
  totVar<-apply(countRanks,1,function(x) sum(abs(x-mean(x))^2))/apply(do.call(rbind, countVar_list), 2, sum)
  totVar[which(totVar=="Inf")]<-NA

  samp_func<-function(countRanks, seed){
    set.seed(seed)
    SampFrame<-countRanks[,sample(1:ncol(countRanks),ncol(countRanks))]
    Samp_countVar1<-apply(SampFrame[,which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[1])],1,function(x) sum(abs(x-mean(x))^2))
    if(length(which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[2]))==1){
      Samp_countVar2<-0
    }
    if(length(which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[2]))>1){
      Samp_countVar2<-apply(SampFrame[,which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[2])],1,function(x) sum(abs(x-mean(x))^2))
    }
    Samp_totVar<-apply(SampFrame,1,function(x) sum(abs(x-mean(x))^2))/(Samp_countVar1+Samp_countVar2)
    Samp_totVar[which(Samp_totVar=="Inf")]<-NA
    return(Samp_totVar)
  }

  if(!is.null(batch_family_variable)){
    batch_fam_samps<-which(colData[,batch_family_variable]==levels(as.factor(colData[,batch_family_variable]))[2])
    batch_fam_statistic_list<-list()
    for(samp in batch_fam_samps[batch_fam_samps!=which(colData[,testVariable]==levels(as.factor(colData[,testVariable]))[2])]){
      batch_fam_countVar1<-apply(countRanks[,-samp],1,function(x) sum(abs(x-mean(x))^2))
      batch_fam_countVar2<-0
      batch_fam_totVar<-apply(countRanks,1,function(x) sum(abs(x-mean(x))^2))/(batch_fam_countVar1+batch_fam_countVar2)
      batch_fam_totVar[which(batch_fam_totVar=="Inf")]<-NA
      batch_fam_statistic_list[[as.character(samp)]]<-batch_fam_totVar
    }
    max_batch_fam_statistic<-apply(do.call(cbind, batch_fam_statistic_list),1,max)
    Samp_totVar_tot<-unlist(parallel::mclapply(1:ceiling((1/minP)/nrow(countData)), function(x) samp_func(countRanks, x),mc.cores = numCores))
    batch_fam_Ps<-((rank(c(-max_batch_fam_statistic,-Samp_totVar_tot),ties.method="max")[1:length(max_batch_fam_statistic)]) - rank(-max_batch_fam_statistic,ties.method="max"))/length(Samp_totVar_tot)
    Ps<-((rank(c(-totVar,-Samp_totVar_tot),ties.method="max")[1:length(totVar)]) - rank(-totVar,ties.method="max"))/length(Samp_totVar_tot)
    resultFrame<-do.call(cbind,list((rank(countVar_list[[which.max(table(colData[,1]))]],ties.method = "max")/length(countVar_list[[which.max(table(colData[,1]))]]))*100,totVar,Ps,batch_fam_Ps))
    colnames(resultFrame)<-c("variance_rank", "test_statistic", "p_value","min_batch_fam_p_value")
    row.names(resultFrame)<-row.names(countData)
    return(resultFrame)
  }

  if(is.null(batch_family_variable)){
    Samp_totVar_tot<-unlist(parallel::mclapply(1:ceiling((1/minP)/nrow(countData)), function(x) samp_func(countRanks, x),mc.cores = numCores))
    Ps<-((rank(c(-totVar,-Samp_totVar_tot),ties.method="max")[1:length(totVar)]) - rank(-totVar,ties.method="max"))/length(Samp_totVar_tot)
    resultFrame<-do.call(cbind,list((rank(countVar_list[[which.max(table(colData[,1]))]],ties.method = "max")/length(countVar_list[[which.max(table(colData[,1]))]]))*100,totVar,Ps))
    colnames(resultFrame)<-c("variance_rank", "test_statistic", "p_value")
    row.names(resultFrame)<-row.names(countData)
    return(resultFrame)
  }
}
