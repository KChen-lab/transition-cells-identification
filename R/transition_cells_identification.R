# To separate cells in transition from cells in stable states, we developed an analytical 
# framework based on the conclusion we drew from the SDEs model: gene pair-wise correlation 
# coefficients for transition cells are close to -1 and 1

library(dplyr)
library(Seurat)
library(scLink)

# We divided our methodology into 4 steps:
# Step 1: Data pre-processing
# Step 2: Selecting cells for Pearson’s correlation coefficients calculation
# Step 3: Top most variable genes selection
# Step 4: Transition index calculation

gene_pearson<-function(data,highly_variable_gene=NULL,n_neighbor=200){
        res<-NULL
        exprs<-as.matrix(GetAssayData(data))
        # selecting top most variable genes according to expression variation
        if(is.null(highly_variable_gene)){
                highly_variable_gene<-order(apply(as.matrix(exprs),1,function(x) var(x)),decreasing=T)[1:100]
                highly_variable_gene<-rownames(exprs)[highly_variable_gene]
        }
        # finding nearest neighbors of each cell in PCA dimensions
        pca_emb<-data@reductions$pca@cell.embeddings[,1:20]
        cell_dist<-as.matrix(dist(pca_emb,diag=T,upper=T))
        for(c in data %>% colnames()){
                if(ncol(data)<n_neighbor)
                next
                expr<-t(exprs[highly_variable_gene,order(cell_dist[c,])[1:n_neighbor]])
                # calculating gene pair-wise Pearson's correlation coefficients after filtering dropouts
                pearson_cor<-sclink_cor(expr = expr, ncores = 5)
                diag(pearson_cor)<-NA
                res<-rbind(res,as.numeric(pearson_cor))
        }
        rownames(res)<-colnames(data)
        gene_name<-highly_variable_gene
        name_tmp<-character()
        for(i in gene_name){
                for(j in gene_name)
                        name_tmp<-c(name_tmp,paste(i,j,sep='_'))
         }
        colnames(res)<-name_tmp
        return(res)
}

find_ks_d<-function(x,y){ # x:control
        n.x <- length(x)
        n.y <- length(y)
        w <- c(x, y)
        z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
        w_sort<-sort(w)
        w_max<-max(w_sort[z==max(z)])
        z_sub<-z[w_sort>w_max]
        w_min<-min(w_sort[w_sort>w_max][z_sub==min(z_sub)])
        return(list('w_min'=w_min,'w_max'=w_max))
}

# main
n_neighbors<-300
res<-list()
# Since we only use cells and their nearest neighbors to calculate gene pair-
# wise Pearson's correlation coefficients, to save computational cost, we do 
# this calculation cluster by cluster.
for(i in data$seurat_clusters %>% unique()){
        print(i)
        data_sub<-subset(data,seurat_clusters == i)
        if(ncol(data_sub)<n_neighbors)
                next
        data_sub<-FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = 2000)
        data_sub<- ScaleData(data_sub,vars.to.regress="nCount_RNA")
        data_sub<- RunPCA(data_sub, features = VariableFeatures(object = data_sub))
        data_sub<- FindNeighbors(data_sub, dims = 1: 20,annoy.metric='cosine')
        # calculating gene pair-wise Pearson's correlation coefficients using gene_pearson function
        pearson<-gene_pearson(data_sub,n_neighbor=n_neighbors)
        # fine tuning  (optional)
        pearson<-abs(pearson)
        res_1<-apply(pearson,1,function(x) ks.test(x,pearson[1,],alternative='greater')$statistic)
        res_2<-apply(pearson,1,function(x) ks.test(x,pearson[1,],alternative='less')$statistic)
        res_tmp<-find_ks_d(pearson[which.max(res_1),],pearson[which.max(res_2),])
        tmp<-apply(pearson,1,function(x) sum(abs(x)>min(unlist(res_tmp)) & abs(x)<max(unlist(res_tmp)),na.rm=T)/sum(abs(x)>=0,na.rm=T))
        data_sub<-AddMetaData(data_sub,data.frame('pearson'=tmp))
        cell<-data_sub@meta.data %>% filter(pearson>quantile(data_sub$pearson,0.8,na.rm=T)) %>% rownames()
        hvg<-sort(apply(as.matrix(GetAssayData(data_sub[,cell])),1,function(x) var(x)),decreasing=T)[1:100] %>% names
        pearson<-gene_pearson(data_sub,n_neighbor=n_neighbors,highly_variable_gene=hvg)
        # calculating transition index based on gene pair-wise Pearson's correlation coefficients
        res_1<-apply(pearson,1,function(x) ks.test(x,pearson[1,],alternative='greater')$statistic)
        res_2<-apply(pearson,1,function(x) ks.test(x,pearson[1,],alternative='less')$statistic)
        res_tmp<-find_ks_d(pearson[which.max(res_1),],pearson[which.max(res_2),])
        res[[as.character(i)]]<-apply(pearson,1,function(x) sum(abs(x)>min(unlist(res_tmp)) & abs(x)<max(unlist(res_tmp)),na.rm=T)/sum(abs(x)>=0,na.rm=T))
}
res<-Reduce(function(x,y) rbind(data.frame('transition_index'=x),data.frame('transition_index'=y)), res)
data<-AddMetaData(data,res)
saveRDS(res,'data_with_transition_index.rds')
