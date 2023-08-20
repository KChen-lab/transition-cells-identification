library(Seurat)
library(dplyr)

# find transition cells
# B cells
stats=data@meta.data %>% filter(celltype=='B',transition_index>0.22) %>% group_by(treatment,hamster) %>% summarize(n=n()) %>% left_join(data@meta.data %>%
    group_by(treatment,hamster) %>% summarize(n_t=n())) %>% mutate(per=n/n_t) %>% group_by(treatment) %>% 
    mutate(mean_n=mean(n),mean_per=mean(per),sd_n=sd(n),sd_per=sd(per))
# A tibble: 16 x 9
# Groups:   treatment [5]
   treatment hamster     n   n_t     per mean_n mean_per  sd_n  sd_per
   <chr>     <chr>   <int> <int>   <dbl>  <dbl>    <dbl> <dbl>   <dbl>
 1 aaUntr    Ha1        24   841 0.0285    12.5   0.0175  8.02 0.00895
 2 aaUntr    Ha2         7   957 0.00731   12.5   0.0175  8.02 0.00895
 3 aaUntr    Ha3         7   485 0.0144    12.5   0.0175  8.02 0.00895
 4 aaUntr    Ha4        12   604 0.0199    12.5   0.0175  8.02 0.00895
 5 adeno2x   Ha1        28  1117 0.0251    21     0.0226  8.89 0.00571
 6 adeno2x   Ha2        11   686 0.0160    21     0.0226  8.89 0.00571
 7 adeno2x   Ha3        24   902 0.0266    21     0.0226  8.89 0.00571
 8 att2x     Ha1        45   996 0.0452    37.3   0.0481 13.3  0.00978
 9 att2x     Ha2        45  1123 0.0401    37.3   0.0481 13.3  0.00978
10 att2x     Ha3        22   373 0.0590    37.3   0.0481 13.3  0.00978
11 mRNA2x    Ha1        15  1043 0.0144    15.3   0.0136  2.52 0.00104
12 mRNA2x    Ha2        13  1046 0.0124    15.3   0.0136  2.52 0.00104
13 mRNA2x    Ha3        18  1285 0.0140    15.3   0.0136  2.52 0.00104
14 mRNAatt   Ha1        29  1090 0.0266    23     0.0315 11.3  0.00514
15 mRNAatt   Ha2        30   814 0.0369    23     0.0315 11.3  0.00514
16 mRNAatt   Ha3        10   322 0.0311    23     0.0315 11.3  0.00514

p1<-ggplot(stats %>% select(treatment,mean_per,sd_per) %>% unique()) +
    geom_bar(aes(x=treatment, y=mean_per,fill=factor(treatment)), stat="identity",alpha=0.7) +
    geom_errorbar( aes(x=treatment, ymin=mean_per-sd_per, ymax=mean_per+sd_per), width=0.4, alpha=0.9) +
    theme_bw()+theme(legend.position="none")+ylab('transition B cells % of all blood cells')

transition_cells<-subset(data,celltype=='B'& transition_index>0.22) %>% colnames()
ctrl<-subset(data,celltype=='B' & transition_index<0.22) %>% colnames()
res<-FindMarkers(data,ident.1=transition_cells,ident.2=ctrl)
res %>% filter(p_val_adj<0.01) %>% write.csv('deg_B_cell_transition_vs_ctrl.csv')

# classical monocytes
pdf('classical_monocyte_transition_index_dist.pdf')
transition_index_sub<-subset(data,celltype=='Classical monocyte')$transition_index
fit <- Mclust(transition_index_sub,G=2)
plot(density(transition_index_sub),main="Transition index distribution for classical monocyte")
rug(jitter(transition_index_sub[fit$classification==2 & fit$z[,2]>0.8]),col='red')
rug(jitter(transition_index_sub[fit$classification==1 & fit$z[,1]>0.8]),col='blue')
dev.off()
monocyte<-subset(data,celltype=='Classical monocyte')
monocyte$transition<-NA
monocyte$transition[names(transition_index_sub[fit$classification==2 & fit$z[,2]>0.8])]<-'transition'
monocyte$transition[names(transition_index_sub[fit$classification==1 & fit$z[,1]>0.8])]<-'stable'
stats=monocyte@meta.data %>% filter(celltype=='Classical monocyte',transition=='transition') %>% group_by(treatment,hamster) %>% summarize(n=n()) %>% left_join(data@meta.data %>%
    group_by(treatment,hamster) %>% summarize(n_t=n())) %>% mutate(per=n/n_t) %>% group_by(treatment) %>% 
    mutate(mean_n=mean(n),mean_per=mean(per),sd_n=sd(n),sd_per=sd(per))
stats
# A tibble: 16 x 9
# Groups:   treatment [5]
   treatment hamster     n   n_t     per mean_n mean_per  sd_n  sd_per
   <chr>     <chr>   <int> <int>   <dbl>  <dbl>    <dbl> <dbl>   <dbl>
 1 aaUntr    Ha1        24   841 0.0285    27.8   0.0430 10.7  0.0238
 2 aaUntr    Ha2        18   957 0.0188    27.8   0.0430 10.7  0.0238
 3 aaUntr    Ha3        26   485 0.0536    27.8   0.0430 10.7  0.0238
 4 aaUntr    Ha4        43   604 0.0712    27.8   0.0430 10.7  0.0238
 5 adeno2x   Ha1        36  1117 0.0322    33.3   0.0374  6.43 0.00497
 6 adeno2x   Ha2        26   686 0.0379    33.3   0.0374  6.43 0.00497
 7 adeno2x   Ha3        38   902 0.0421    33.3   0.0374  6.43 0.00497
 8 att2x     Ha1        22   996 0.0221    21     0.0285  5.57 0.0102
 9 att2x     Ha2        26  1123 0.0232    21     0.0285  5.57 0.0102
10 att2x     Ha3        15   373 0.0402    21     0.0285  5.57 0.0102
11 mRNA2x    Ha1        27  1043 0.0259    28.7   0.0259  3.79 0.00566
12 mRNA2x    Ha2        33  1046 0.0315    28.7   0.0259  3.79 0.00566
13 mRNA2x    Ha3        26  1285 0.0202    28.7   0.0259  3.79 0.00566
14 mRNAatt   Ha1        25  1090 0.0229    19.7   0.0234 14.7  0.0144
15 mRNAatt   Ha2        31   814 0.0381    19.7   0.0234 14.7  0.0144
16 mRNAatt   Ha3         3   322 0.00932   19.7   0.0234 14.7  0.0144

p2<-ggplot(stats %>% select(treatment,mean_per,sd_per) %>% unique()) +
    geom_bar(aes(x=treatment, y=mean_per,fill=factor(treatment)), stat="identity",alpha=0.7) +
    geom_errorbar( aes(x=treatment, ymin=mean_per-sd_per, ymax=mean_per+sd_per), width=0.4, alpha=0.9) +
    theme_bw()+theme(legend.position="none")+ylab('transition classical monocytes % of all blood cells')

transition_cells<-subset(monocyte,transition=='transition') %>% colnames()
ctrl<-subset(monocyte,transition=='stable') %>% colnames()
res<-FindMarkers(data,ident.1=transition_cells,ident.2=ctrl)
res %>% filter(p_val_adj<0.01) %>% write.csv('deg_classical_monocytes_transition_vs_ctrl.csv')

# platelets
stats=data@meta.data %>% filter(celltype=='Platelet',transition_index>0.5) %>% group_by(treatment,hamster) %>% summarize(n=n()) %>% left_join(data@meta.data %>%
    group_by(treatment,hamster) %>% summarize(n_t=n())) %>% mutate(per=n/n_t) %>% group_by(treatment) %>% 
    mutate(mean_n=mean(n),mean_per=mean(per),sd_n=sd(n),sd_per=sd(per))
# A tibble: 9 x 9
# Groups:   treatment [5]
  treatment hamster     n   n_t      per mean_n mean_per  sd_n   sd_per
  <chr>     <chr>   <int> <int>    <dbl>  <dbl>    <dbl> <dbl>    <dbl>
1 aaUntr    Ha4         1   604 0.00166    1    0.00166  NA    NA
2 adeno2x   Ha2         1   686 0.00146    1    0.00146  NA    NA
3 att2x     Ha1         5   996 0.00502    7.33 0.0103    3.21  0.00555
4 att2x     Ha2        11  1123 0.00980    7.33 0.0103    3.21  0.00555
5 att2x     Ha3         6   373 0.0161     7.33 0.0103    3.21  0.00555
6 mRNA2x    Ha2         1  1046 0.000956   1    0.000956 NA    NA
7 mRNAatt   Ha1         4  1090 0.00367    3.33 0.00618   1.15  0.00544
8 mRNAatt   Ha2         2   814 0.00246    3.33 0.00618   1.15  0.00544
9 mRNAatt   Ha3         4   322 0.0124     3.33 0.00618   1.15  0.00544

p3<-ggplot(stats %>% select(treatment,mean_per,sd_per) %>% unique()) +
    geom_bar(aes(x=treatment, y=mean_per,fill=factor(treatment)), stat="identity",alpha=0.7) +
    geom_errorbar( aes(x=treatment, ymin=mean_per-sd_per, ymax=mean_per+sd_per), width=0.4, alpha=0.9) +
    theme_bw()+theme(legend.position="none")+ylab('transition platelet % of all blood cells')
ggarrange(p1,p2,p3,ncol=3)

transition_cells<-subset(data,celltype=='Platelet' & pearson>0.5) %>% colnames()
ctrl<-subset(data,celltype=='Platelet' & pearson<0.5) %>% colnames()
res<-FindMarkers(data,ident.1=transition_cells,ident.2=ctrl)
res %>% filter(p_val_adj<0.01) %>% write.csv('deg_platelets_transition_vs_ctrl.csv')
