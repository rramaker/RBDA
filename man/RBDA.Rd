% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RBDA.R
\name{RBDA}
\alias{RBDA}
\title{Rank Based Differential Expression Analysis
This function identifies differentially expressed genes from a raw count table}
\usage{
RBDA(countData, colData, testVariable, batch_family_variable = NULL,
  randomSeed = 1990, minP = 5e-06, numCores = 4)
}
\arguments{
\item{countData}{A raw count table (dataframe or matrix) with samples as columns and genes as rows.}

\item{colData}{A dataframe with samples as rows and columns factor variables indicating experimental groups.}

\item{testVariable}{Character indicating the column name in colData that reresents the variable of interest for differential expression analysis.}

\item{batch_family_variable}{This paramater is used in N=1 cases only. Character that specifies a column name in colData that indicates other samples sequenced in the same batch as the sample of interest or related family members. It can be used as a downstream filter when comparing one sample to a large reference set. Defaults to NULL.}

\item{randomSeed}{A numeric indicating a random seed for reproducible analysis. Defaults to 1990}

\item{minP}{A numeric between 0 and one indicating the minimum possible P-value computed via random sampling from the countTable. The larger this number is the faster the compute time will be. Defaults to 0.000005.}

\item{numCores}{A numeric indicating the number of cores used for computing. Most laptops can readily handle up to 4. Defaults to 4.}
}
\value{
Returns a data frame containing differential expression results. The first column, "variance_rank," indicates the percentile rank of a gene's rank variance. Low numbers indicate a gene exhibits very low variance across samples. The second column "test_statistic" is the RBDA test statistic. The third column "p_value" indicates the probability of observing the test statistic by chance alone. A fourth column, "min_batch_fam_p_value," will only result if the batch_family_variable parameter is specified and will provide the minimum p-value computed for any samples indicated by this factor variable in colData.
}
\description{
Rank Based Differential Expression Analysis
This function identifies differentially expressed genes from a raw count table
}
\examples{
data("BRCA_Counts")
colData<-data.frame(testVariable=as.factor(c(rep(0,5),rep(1,5))), row.names=colnames(BRCA_Counts))
RBDA(countData = BRCA_Counts,colData,testVariable = "testVariable", randomSeed = 1990, numCores=1)
}

