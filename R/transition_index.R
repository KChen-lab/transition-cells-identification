#' Calculate GPPCCs
#'
#' This is the function calculating gene pair-wise Pearson's correlation coefficient (GPPCCs). 
#' The GPPCCs calculated using our method are assumed from the cells that exhibit similar 
#' regulatory profiles and around equilibrium points. Thus, for each cell, we select the cell 
#' together with its nearest neighbors in Principal Component Analysis (PCA) dimensions based 
#' on cosine distance to calculate GPPCCs. We also use scLink in the calculation to deal with
#' the gene expression values affected by dropout. 
#' @param data Seurat object containing both count matrix and metadata such as neighboring 
#' results for a single-cell dataset.
#' @param highly_variable_gene A gene list containing top most variable genes. Calculated 
#' by \code{var} function by default.
#' @param n_neighbor The number of neighboring cells used to calculate GPPCCs. 300 by default.
#' @return A matrix of GPPCCs. Rows are cells and columns are gene pairs.
#' @import Seurat
#' @import scLink
#' @import dplyr
#' @export
gene_pearson<-function(data,highly_variable_gene=NULL,n_neighbor=300,n_gene=50,...){
        res<-NULL
        exprs<-as.matrix(GetAssayData(data))
        if(is.null(highly_variable_gene)){
                highly_variable_gene<-order(apply(as.matrix(exprs),1,function(x) var(x)),decreasing=T)[1:n_gene]
                highly_variable_gene<-rownames(exprs)[highly_variable_gene]
        }
        for(c in data %>% colnames()){
                expr<-t(exprs[highly_variable_gene,data@neighbors$RNA.nn@nn.idx[match(c,data@neighbors$RNA.nn@cell.names),]])
                pearson_cor<-scLink::sclink_cor(expr = expr, ncores = 5)
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


find_ks_d<-function(x,y){
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

#' Calculate transition index
#'
#' This is the function calculating transition index that reflects the transition
#' probability of a cell. The transition index we defined is inspired by Kolmogorov-Smirnov 
#' statistic. We assume there are both transition cells and stable cells captured by  
#' scRNA-seq in a developmental process. We first find the archetype of transition 
#' cell and stable cell by using modified Kolmogorov-Smirnov statistic. As proved 
#' mathematically, there are more gene pairs whose absolute value of Pearson’s 
#' correlation coefficients are close to 1 in transition cells than stable cells. 
#' Thus, based on the distribution difference, by finding the maximum difference 
#' among all cells, we identify archetypal stable cells and transition cells. To 
#' count and compare the number of gene pairs exceeding such threshold in each cell, 
#' we define transition index by summarizing the percentage of gene pairs whose 
#' absolute value of Pearson’s correlation coefficients are between argmax_x(D_max)    
#' and argmin_x(D_max). 
#'
#' @param data Seurat object containing both count matrix and metadata such as neighboring 
#' results for a single-cell dataset.
#' @param highly_variable_gene A gene list containing top most variable genes. Calculated 
#' by \code{var} function by default.
#' @param group The cell group highly_variable_genes selected based on. It should be a column
#' name of data@meta.data, such as 'seurat_clusters','cell_type','time_point'. Using 
#' 'seurat_clusters' by default. 
#' @param n_neighbor The number of neighboring cells used to calculate GPPCCs. 300 by default.
#' @param n_gene The number of top most variable genes used in calculating GPPCCs. 50 by default.
#' @return A Seurat object. Calculated transition index is store as a column in the metadata.
#' @import Seurat
#' @import scLink
#' @import dplyr
#' @export
transition_index<-function(data,highly_variable_gene=NULL,group='seurat_clusters',n_neighbors=300,n_gene=50,...){
    res<-list()
    for(i in data@meta.data[,group] %>% unique()){
            print(i)
            data_sub<-data[,data@meta.data[,group]==i]
            if(ncol(data_sub)<=n_neighbors)
                    next
            data_sub<-Seurat::FindVariableFeatures(data_sub, selection.method = "vst", nfeatures = 2000)
            data_sub<- Seurat::ScaleData(data_sub,vars.to.regress="nCount_RNA")
            data_sub<- Seurat::RunPCA(data_sub, features = VariableFeatures(object = data_sub),npcs=min(20,ncol(data_sub)))
            data_sub<- Seurat::FindNeighbors(data_sub, dims = 1: min(20,ncol(data_sub)),return.neighbor = TRUE,k.param=n_neighbors,annoy.metric='cosine')
            pearson<-gene_pearson(data_sub,highly_variable_gene=highly_variable_gene,n_neighbor=n_neighbors,n_gene=n_gene)
            # fine tuning  (optional)
            pearson<-abs(pearson)
            res_1<-apply(pearson,1,function(x) ks.test(x,pearson[1,],alternative='greater')$statistic)
            res_2<-apply(pearson,1,function(x) ks.test(x,pearson[1,],alternative='less')$statistic)
            res_tmp<-find_ks_d(pearson[which.max(res_1),],pearson[which.max(res_2),])
            tmp<-apply(pearson,1,function(x) sum(abs(x)>min(unlist(res_tmp)) & abs(x)<max(unlist(res_tmp)),na.rm=T)/sum(abs(x)>=0,na.rm=T))
            data_sub<-Seurat::AddMetaData(data_sub,data.frame('pearson'=tmp))
            if(min(tmp)==max(tmp)){
                    res[[as.character(i)]]<-pearson
                    next
            }
            cell<-data_sub@meta.data %>% dplyr::filter(pearson>quantile(data_sub$pearson,0.8,na.rm=T)) %>% rownames()
            hvg<-sort(apply(as.matrix(GetAssayData(data_sub[,cell])),1,function(x) var(x)),decreasing=T)[1:n_gene] %>% names
            pearson<-gene_pearson(data_sub,n_neighbor=n_neighbors,highly_variable_gene=hvg)
            res[[as.character(i)]]<-pearson
    }
    res<-Reduce(function(x,y) rbind(x,y),res)
    res<-abs(res)
    res_1<-apply(res,1,function(x) ks.test(x,res[1,],alternative='greater')$statistic)
    res_2<-apply(res,1,function(x) ks.test(x,res[1,],alternative='less')$statistic)
    res_tmp<-find_ks_d(res[which.max(res_1),],res[which.max(res_2),])
    tmp<-apply(res,1,function(x) sum(abs(x)>min(unlist(res_tmp)) & abs(x)<max(unlist(res_tmp)),na.rm=T)/sum(abs(x)>=0,na.rm=T))
    data<-Seurat::AddMetaData(data,data.frame('transition_index'=tmp))
    return(data)
}


