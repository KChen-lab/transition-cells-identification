# For performance evaluation, we compared our method with CellRank and Mutrans. Our method 
# uses transition index to reflect the probability. While for CellRank, we calculated the 
# probability of each cell not being in a macrostate. We used transition entropy provided 
# by Mutrans to correspond the probability of a cell to be a transition cell. The performance 
# of these methods was compared using receiver operating characteristic curve.

library(dplyr)
library(Seurat)
library(pROC)

compute_roc_cellrank<-function(file){
	cellrank<-read.csv(file)
	cellrank$X<-sapply(cellrank$X,function(x) gsub('\\.','-',x))
	cellrank$X[grep('^X[0-9].*',cellrank$X)]<-sub('X','',cellrank$X[grep('^X[0-9].*',cellrank$X)])
	cellrank<-cellrank %>% mutate('non_terminal'=1-terminal_states_probs)
	data<-AddMetaData(data,cellrank %>% select(X,non_terminal) %>% tibble::column_to_rownames('X'))
	res<-roc(data$celltype,data$non_terminal)
	print(res)
	return(data)
}

compute_roc_mutran<-function(file){
    mutrans<-read.csv(file,row.names=1)
    data<-AddMetaData(data,mutrans)
    res<-roc(data$celltype,data$entropy)
    print(res)
    return(data)
}

compute_roc_pearson<-function(file){
	data<-readRDS(file)
	res<-roc(data$celltype,data$transition_index)
	print(res)
	return(data)
}
