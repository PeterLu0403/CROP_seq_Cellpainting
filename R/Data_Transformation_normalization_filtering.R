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

setwd("~/Project/Dataset_20191216/")
df<-fread("Raw_merged.csv",header = T,sep = ",")

## remove location features:
loc<-c("Center","Location","Parent","Children")
cn<-colnames(df)
cn<-cn[which(!str_detect(cn, paste(loc,collapse = "|")))]
## remove Object_Number from features: 
cn<-cn[which(!str_detect(cn,"Object_Number"))]
## remove eulernumber
cn<-cn[which(!str_detect(cn,"EulerNumber"))]
## extract metadata and variables
meta<-cn[1056:1062]
variables<-setdiff(cn,meta)[-(1:2)]

object_m<-dplyr::select(df,meta) %>% separate(Gene,into=(c("Number","Gene","seq"))) %>% mutate(Gene = ifelse(str_detect(Gene,"Human"),paste(Gene,seq,sep = "_"),Gene))
object_m<-dplyr::select(object_m,meta) %>% mutate(Batch_Plate_Gene=paste(Batch,Plate,Gene,sep = "_"))

meta<-colnames(object_m)
object_v<-dplyr::select(df,variables)
object<-cbind(object_m,object_v)
Genes<-unique(object$Gene)
Controls<-sort(Genes[which(str_detect(Genes,"Human"))])
Genes<-sort(Genes[which(!str_detect(Genes,"Human"))])

# Transformation:
## shift data to get rid of negative values
df_shift<-apply(object %>% dplyr::select(one_of(variables)), 2, function(x) x+abs(min(x, na.rm = T))+1)
df_shift<-cbind(object_m,df_shift)

# select control groups
df_control<-dplyr::filter(df_shift,str_detect(Gene,"Human"))

## Calculate lambda value of boxcox transformation in control groups and do BoxCox transformation
boxcox_lam<-function(x){
  temp<-boxcox(x,lambda = c(-10,10),objective.name = "Log-Likelihood",optimize = T)
  lambda<-temp[[1]]
  return(lambda)
}
lambda<-apply(df_control %>% select(one_of(variables)),2,boxcox_lam)
write.csv(lambda,"csv_files/boxcox.csv",row.names = F)
lambda<-fread("csv_files/boxcox.csv",header = T,sep = ",") %>% as.matrix() %>% as.vector()

temp<-df_control %>% select(one_of(variables))
df_boxcox<-NULL
for (i in 1:length(variables)) {
  df_boxcox[[i]]<-boxcoxTransform(temp[,..i] %>% as.matrix() %>% as.vector(),lambda[i])
}
names(df_boxcox)<-variables
df_boxcox<-as.data.frame(df_boxcox)
df_boxcox<-cbind(df_control %>% select(meta),df_boxcox)


## Transform data using generalized_log:
df_transformed<-cytominer::transform(df_control,variables,operation="generalized_log")

## Distribution testing
df_ad<-data.frame(Variables=variables, raw_data=numeric(length(variables)),BoxCox=numeric(length(variables)),Generalized_log=numeric(length(variables)))
for (i in 1:length(variables)){
  ad<-ad.test(select(df_control,variables)[,..i] %>% as.matrix() %>% as.vector())
  df_ad$raw_data[i]<-ad$statistic
  ad<-ad.test(select(df_boxcox,variables)[,..i] %>% as.matrix() %>% as.vector())
  df_ad$BoxCox[i]<-ad$statistic
  ad<-ad.test(select(df_transformed,variables)[,..i] %>% as.matrix() %>% as.vector())
  df_ad$Generalized_log[i]<-ad$statistic
  }
df_ad<-gather(df_ad,"group","value",2:4)
## Plot ad test
p<-ggplot(df_ad,aes(x=reorder(Variables,-value),y=value,group=group))+geom_line(aes(color=group),size=0.1)
p<-p+theme_bw()+ theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(size=0.8,colour="black"),
        plot.title = element_text(size= 16,face = "bold"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p<-p+ylab("Statistic Value A")+xlab("Features")
p<-p+scale_color_discrete(name="Methods:")
p<-p+ggtitle("Aderson-Darling normality test")
p
pdf("Figures/Normality test.pdf",width = 13,height = 8,bg="transparent")
p
dev.off()

# Do boxcox transformation in all groups:
temp<- df_shift %>% select(one_of(variables))
df_boxcox<-NULL
for (i in 1:length(variables)) {
  df_boxcox[[i]]<-boxcoxTransform(temp[,..i] %>% as.matrix() %>% as.vector(),lambda[i])
}
names(df_boxcox)<-variables
df_boxcox<-as.data.frame(df_boxcox)

df_boxcox<-cbind(df_shift %>% select(meta),df_boxcox)
df_boxcox_control<-df_boxcox %>% filter(str_detect(Gene,"Human"))
# Normalization, drop features with NAs
df_normalized<-cytominer::normalize(df_boxcox,variables,strata=c("Batch","Plate"),sample=df_boxcox_control,operation = "robustize")
variables<-setdiff(colnames(df_normalized),meta)
df_clean<-cytominer::variable_select(df_normalized,variables,operation = "drop_na_columns",cutoff=0)

# Batch effects evaluation 
rm(df,df_boxcox,df_boxcox_control,df_control,df_normalized,df_shift,object,object_m,object_v,temp)
Spear_cor<-function(df,group) {
  temp <- dplyr::filter(df,Gene==group)
  temp_t <- temp %>% select(one_of(variables)) %>% as.matrix()
  rownames(temp_t)<-paste(temp$Batch,temp$Plate,temp$Well,"Position",temp$Position,temp$ObjectNumber,sep="_")
  temp_t<-t(temp_t)
  Spear_cor<-cor(temp_t,method = "spearman")
  order<-corrMatOrder(Spear_cor,order = "FPC")
  Spear_cor<-Spear_cor[order,order]
  return(Spear_cor)
}
## Calculate Spearman correlation value for each group:
df_spear<-NULL
for (i in 1:length(Genes)) {
  df_spear[[i]]<-Spear_cor(df_clean,Genes[i])
}
names(df_spear)<-Genes

df_spear_control<-NULL
for (i in 1:length(Controls)) {
  df_spear_control[[i]]<-Spear_cor(df_clean,Controls[i])
}
names(df_spear_control)<-Controls

saveRDS(df_spear,"spear_cor_obj.rds")
saveRDS(df_spear_control,"spear_cor_control_obj.rds")

## Do plot of the Spearman matrix of each group: 
col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))

for (i in 1:length(Genes)) {
  png(file = paste0("sc_Spearman_Cor/",Genes[i],".png"),height=5000,width=7000,res = 300)
  print(Heatmap(df_spear[[i]],name = Genes[i],show_row_names = F,show_column_names = F, cluster_rows = F, cluster_columns = F, col = col_fun))
  dev.off()
}
for (i in 1:length(Controls)) {
  png(file = paste0("sc_Spearman_Cor/",Controls[i],".png"),height=5000,width=7000,res = 300)
  print(Heatmap(df_spear_control[[i]],name = Controls[i],show_row_names = F,show_column_names = F, cluster_rows = F, cluster_columns = F, col = col_fun))
  dev.off()
}

## calculate mean value of absolute Spearman correlation value of each group
mean_spear_abs<-function(df, group) {
  temp<-abs(df)
  avg<-apply(temp,2,mean)
  temp<-data.frame(obj=names(avg),spear=avg,Group=group)
  return(temp)
}

df_spear_mean_genes<-NULL
for (i in 1:length(Genes)) {
  df_spear_mean_genes[[i]]<-mean_spear_abs(df_spear[[i]],Genes[i])  
}
df_spear_mean_genes<-do.call(rbind,df_spear_mean_genes)

df_spear_mean_controls<-NULL
for (i in 1:length(Controls)) {
  df_spear_mean_controls[[i]]<-mean_spear_abs(df_spear_control[[i]],Controls[i])
}
df_spear_mean_controls<-do.call(rbind,df_spear_mean_controls)

df_spear_mean<-rbind(df_spear_mean_genes,df_spear_mean_controls)
write.csv(df_spear_mean,"csv_files/spear_cor_mean_obj.csv",row.names = F)

## draw a plot of the mean value of absolute Spearman correlation value of each group
df_spear_mean<-fread("csv_files/spear_cor_mean_obj.csv",header = T,sep = ",")
p<-ggplot(df_spear_mean,aes(Group,spear,color=Group))
p<-p+ geom_violin()+geom_point(position = position_jitter(width = 0.15),size=0.1)+theme_bw()+
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(size=0.8,colour="black"),
        plot.title = element_text(size= 16,face = "bold"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.text.x = element_text(size = 12,angle = 90,hjust = 1,vjust = 0.5, face = "plain",color="black"))
p<-p+ylab("mean absolute Spearman correlation value")+xlab("")+ylim(0,0.5)
p<-p+scale_color_discrete(name="Group:")
p
png(file = "Figures/mean absolute Spearmans cor.png",height=5000,width=7000,res = 300)
print(p)
dev.off()

## Filter the data at the cutoff of 0.15:
filter_spear<-function(df_matrix,df_spear_filter,group) {
  temp<-filter(df_spear_filter,Group==group)
  df_matrix_filter<-df_matrix[temp$obj,temp$obj]
  return(df_matrix_filter)
}
# 0.15: 
df_spear_filter_0.15<-filter(df_spear_mean,spear>=0.15)
df_spear_genes_filter_0.15<-NULL
for (i in 1:length(Genes)) {
  df_spear_genes_filter_0.15[[i]]<-filter_spear(df_spear[[i]],df_spear_filter_0.15,Genes[i])
}
df_spear_control_filter_0.15<-NULL
for (i in 1:length(Controls)) {
  df_spear_control_filter_0.15[[i]]<-filter_spear(df_spear_control[[i]],df_spear_filter_0.15,Controls[i])
}
# re-draw the filtered Spearman correlation matrix of each group
#0.15:
for (i in 1:length(Genes)) {
  png(file = paste0("0.15/",Genes[i],".png"),height=5000,width=7000,res = 300)
  print(Heatmap(df_spear_genes_filter_0.15[[i]],name = Genes[i],show_row_names = F,show_column_names = F, cluster_rows = F, cluster_columns = F, col = col_fun))
  dev.off()
}
for (i in 1:length(Controls)) {
  png(file = paste0("0.15/",Controls[i],".png"),height=5000,width=7000,res = 300)
  print(Heatmap(df_spear_control_filter_0.15[[i]],name = Controls[i],show_row_names = F,show_column_names = F, cluster_rows = F, cluster_columns = F, col = col_fun))
  dev.off()
}



