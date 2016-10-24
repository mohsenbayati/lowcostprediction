rm(list =ls())

# This is the code we used to create a synthetic version of Kaggle data.
# Responses are not changed but lab values are randomly sampled from the distribution of all patients
# missing patterns are also preserved


# Uncomment each line and run the code once.
fileName = "OLR30_data_wNAs.tsv"
#fileName = "OLR30_data_meanImputed.tsv"
#fileName = "STL_data_wNAs.tsv"
#fileName = "STL_data_meanImputed.tsv"

XY = read.table(paste("development files/data/",fileName,sep=""),sep='\t',header = TRUE)

# Separate features from responses
nc = ncol(XY)
K = 5
p = nc-K
n = nrow(XY)


X = XY[,1:p]
Y = XY[,(p+1):nc]

missingPattern = is.na(X)

X=XY[,1:p]
Xsyn = X

for (c in 1:p){
  inds = !is.na(X[,c])
  Xsyn[,c]=sample(X[inds,c],n,replace=T)
}

Xsyn[missingPattern] = NA

XY = cbind(Xsyn,Y)

destination = paste(paste(getwd(),"/published/data/",sep=""),fileName,sep="")

XY = write.table(XY,file=destination,sep='\t')