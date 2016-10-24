# This code reads the kaggle data for OLR30 and performs HD and MICE multiple imputations and outputs AUC
# confidence intervals for each of the K tasks

rm(list =ls())

source('functionLibrary.R')

#---------------------------------------------------------------------
#---------------- Read and preprocess the data

XY = read.table('data/OLR30_data_wNAs.tsv',sep='\t',header = TRUE)
#XY = read.table('data/OLR30_data_meanImputed.tsv',sep='\t',header = TRUE)
#XY = read.table('data/STL_data_wNAs.tsv',sep='\t',header = TRUE)
#XY = read.table('data/STL_data_meanImputed.tsv',sep='\t',header = TRUE)



# XY is a data.frame with each row a patient and p columns that are features (biomarkers) and
# its last K columns are the binary responses for each disease (1 = yes, 0 = no)
# For STL the p features are not filtered and it contains all possible biomarkers in Kaggle data
# For OLR-30 a different version of the data should be read and the p features should be only the biomarkers 
# from Table 3 of the paper are used. 

# For this code XY is not pre-imputed and it contains NAs among the p features

# Separate features from responses
nc = ncol(XY)
X = XY[,1:(nc-5)]
Y = XY[,(nc-4):nc]

K = ncol(Y)

auc_CI = matrix(rep(0,4*K),nrow=K,ncol=4)
colnames(auc_CI)=c('HD_L','HD_U','MICE_L','MICE_U')
rownames(auc_CI) = colnames(Y)

stm = proc.time()

for (t in 1:K)
{
  XYt = as.data.frame(cbind(X,Y[,t]))
  result = MultipleIMputationCI(XYt, nMI=1000, impMethod = "sample", miceMaxIt = 30)
  auc_CI[t,c(1,2)] = result
  
  result = MultipleIMputationCI(XYt, nMI=1000, impMethod = "pmm", miceMaxIt = 30)
  auc_CI[t,c(3,4)] = result
  
  etm = proc.time()
  cat('Task',t,'finished in',etm[3]-stm[3],'seconds\n')
  stm = etm
}

print(auc_CI)
