library(data.table)
library(stringr)
library(dplyr)
library(tidyverse)

library(AID)
library(EnvStats)
library(cytominer)
library(nortest)
library(ggplot2)
library(corrplot)
library(stringr)
library(ComplexHeatmap)

# loading the filtering data: 
df_clean_filtered<-fread("comparison/df_clean_filtered.csv",header = T, sep = ",")
variables<-fread("comparison/variables.csv",header = T,sep = ",")
variables<-variables$variables
Genes<-unique(df_clean_filtered$Gene)
Controls<-sort(Genes[which(str_detect(Genes,"Human"))])
Genes<-sort(Genes[which(!str_detect(Genes,"Human"))])

# do a t-test of each gene knock-out groups vs. controls of each features
t.test.func<-function(df,gene,variable) {
  df<-dplyr::select(df,"Gene",variable)
  df_c<-dplyr::filter(df,str_detect(Gene,"Human"))
  df_g<-dplyr::filter(df,Gene==gene)
  res<-t.test(df_g[,2],df_c[,2], var.equal = TRUE)
  temp<-data.frame(feature=variable,mean_gene=res$estimate["mean of x"],mean_control=res$estimate["mean of y"],pval=res$p.value)
  rownames(temp)<-variable
  temp<-temp %>% mutate(diff= mean_gene - mean_control)
  return(temp)
}

df_t_test<-list()
for (j in 1:length(Genes)) {
  temp<-list()
  for (i in 1:length(variables)) {
    temp[[i]]<-t.test.func(df_clean_filtered,Genes[j],variables[i])
  }
  df_t_test[[j]]<-do.call(rbind,temp)
}
names(df_t_test)<-Genes
saveRDS(df_t_test,"df_t_test_groups.rds")

# load the t_test result
df_t_test<-readRDS("df_t_test_groups.rds")

# extract pval of each comparison to the controls
df_t_test_pval<-list()
for (i in 1:length(Genes)) {
  temp<-df_t_test[[i]]
  temp<-data.frame(pval=temp$pval)
  rownames(temp)<-variables
  temp<-as.data.frame(t(temp))
  rownames(temp)<-Genes[i]
  temp$Gene<-Genes[i]
  df_t_test_pval[[i]]<-temp
}
df_t_test_pval<-do.call(rbind,df_t_test_pval)
write.csv(df_t_test_pval,"df_clean_t_test_filtered.csv",row.names = F)

# do plot pval of the t_test results 
df_t_test_pval<-fread("df_clean_t_test_filtered.csv",header = T,sep = ",")
df <- df_t_test_pval %>% melt(id="Gene",measure.vars = variables)
df <- df %>% mutate(log=-log10(value))

## Color xlabs and reorder: 
df_v<-data.frame(feat=variables)
df_v<-mutate(df_v,color = case_when(!str_detect(feat,"Correlation") & str_detect(feat,"corNucleus") ~ 1,
                                    !str_detect(feat,"Correlation") & str_detect(feat,"corER") ~ 2,
                                    !str_detect(feat,"Correlation") & str_detect(feat,"corRNA") ~ 3,
                                    !str_detect(feat,"Correlation") & str_detect(feat,"corAGP") ~ 4,
                                    !str_detect(feat,"Correlation") & str_detect(feat,"corMito") ~ 5,
                                    str_detect(feat,"Correlation") ~ 6,
                                    TRUE ~ 0))
df_v<-df_v[order(df_v$color),]
color <- case_when(!str_detect(df_v$feat,"Correlation") & str_detect(df_v$feat,"corNucleus") ~ "blue",
                   !str_detect(df_v$feat,"Correlation") & str_detect(df_v$feat,"corER") ~ "green",
                   !str_detect(df_v$feat,"Correlation") & str_detect(df_v$feat,"corRNA") ~ "cyan",
                   !str_detect(df_v$feat,"Correlation") & str_detect(df_v$feat,"corAGP") ~ "red",
                   !str_detect(df_v$feat,"Correlation") & str_detect(df_v$feat,"corMito") ~ "pink",
                   str_detect(df_v$feat,"Correlation") ~ "black",
                   TRUE ~ "grey")

df<-left_join(df,df_v,by=c("variable"="feat"))
df_order<-df_v$feat
## plot:
p<-ggplot(df,aes(factor(df$variable),df$log,color=factor(df$Gene)))
p<-p+geom_point(size=0.5)+theme_bw()+theme(panel.grid=element_blank(),
                                           panel.border=element_blank(),
                                           axis.line=element_line(size=1,colour="black"),
                                           axis.text.x = element_text(size = 2,angle = 90, vjust = 0.5, hjust = 1,face = "bold",colour = color),
                                           axis.title.y = element_text(size = 10, face = "bold"),axis.title.x = element_text(size = 10, face = "bold"),
                                           legend.position = "right"
)
p<-p+scale_x_discrete(limits=df_order)
p<-p+ylab("-log10(pvalue)")+xlab("Morphological Features")+scale_color_manual("Group",values = rainbow(62))

pdf("pval_t_test_groups.pdf",width = 40,height = 10,bg="transparent")
p
dev.off()


# draw volcano plot for each groups

v_plots <- function(df,gene) {
  df<-mutate(df,color=case_when(pval<1e-40 & diff > 1.5 ~ "Upregulated",
                                pval<1e-40 & diff < (-1.5) ~ "Downregulated",
                                TRUE ~ "Black"))
  n_up<-nrow(filter(df,color=="Upregulated"))
  n_down<-nrow(filter(df,color=="Downregulated"))
  n_outlier<-nrow(filter(df,abs(diff) >50))
  
  p<-ggplot(df,aes(y=-log10(pval),x=diff)) + geom_point(aes(color=as.factor(color)),alpha=0.44)
  p<-p+geom_hline(yintercept = -log10(1e-40),linetype="dashed")+
    geom_vline(xintercept = -1.5,linetype="dashed")+
    geom_vline(xintercept = 1.5,linetype="dashed")+xlim(-50,50)
  p<-p+theme(panel.grid.major = element_line(size=0.5,colour = "grey96"),
             panel.grid.minor = element_line(size=0.2,colour = "grey96"),
             panel.border = element_blank(),
             axis.line = element_line(size=1,colour = "black"),
             legend.position = "right",
             axis.title.y = element_text(size = 12),
             axis.title.x = element_text(size=12))
  p<-p+scale_color_manual(name=NULL,
                          breaks=c("Upregulated","Downregulated","Black"),
                          limits=c("Upregulated","Downregulated","Black"),
                          values=c("green","red","grey40"),
                          labels=c("Up_features","Down_features","Not Significant"))
  p<-p+labs(title = gene,subtitle = paste0("n_upfeat: ",n_up,"\n","n_downfeat: ",n_down,"\n","n_outlier: ",n_outlier))
  p<-p+scale_y_continuous(breaks = seq(0,320,by=20),limits = c(0,320))
  print(p)
  
}

pdf("vocalno_plot.pdf",height = 6,width=18,bg="transparent")
for (i in 1:length(df_t_test)) {v_plots(df_t_test[[i]],names(df_t_test)[i])} 
dev.off()

for (i in 2) {v_plots(df_t_test[[i]],names(df_t_test)[i])} 



