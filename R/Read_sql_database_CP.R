## This script is to extract features from individual sqlite database and merge them together
library(RSQLite)
library(dplyr)
library(stringr)
library(data.table)
library(tidyverse)

# The function of reading object file from sqlite database:
object_extraction<- function(sql.path) {
  db <- dbConnect(RSQLite::SQLite(), sql.path)
  tables<-dbListTables(db)
  lDataFrames <- vector("list", length=length(tables))
  for (i in seq(along=tables)) {
    lDataFrames[[i]] <- dbGetQuery(conn=db, statement=paste("SELECT * FROM '", tables[[i]], "'", sep=""))
  }
  image<-lDataFrames[[which(str_detect(tables,"Image"))]]
  object<-lDataFrames[[which(str_detect(tables,"Object"))]]
  meta_image<-c("ImageNumber",
                "Metadata_Cell",
                "Metadata_Plate",
                "Metadata_Position",
                "Metadata_Scene",
                "Metadata_Well")
  image<-image[which(str_detect(colnames(image),paste(meta_image,collapse = "|")))]
  names(image)<-c("ImageNumber","Batch","Plate","Position","Scene","Well")
  image$ImageNumber<-as.integer(image$ImageNumber)
  object<- inner_join(object,image,by = "ImageNumber")
  object<- object %>% dplyr::mutate(Batch_Plate_Well=paste(paste(Batch,Plate,sep = "_"),Well,sep = "_"))
  return(object)
}

# The experiment contains 3 batches of plates (B1, D2, and D3) with each batch contains 5 plates (P1-P5):
setwd("~/Project/Dataset_20191216/")
folders<-list.files(pattern = "(B|D)\\d{1}_P\\d{1}")
# First, extract the features from each plate: 
for (f in 1:length(folders)) {
  folder <- folders[f]
  temp<-list.files(path=paste0("./",folder),pattern = "batch_\\d+_out")
  path<-vector()
  for (i in 1:length(temp)) {
    path[i]<-paste(folder,temp[i],"hTMC_analysis.db",sep = "/")
  }
  object<-list()
  for (i in 1:length(temp)) {
    object[[i]]<-object_extraction(path[i])
  }
  df<-do.call(rbind,object)
  write.csv(df,paste0("DB_CSV/",folder,".csv"),row.names = F)
}

# Add metadata of which gene was knockout in each well to each plate:
files<-list.files(path="./DB_CSV", pattern = "(B|D)\\d{1}_P\\d{1}.csv")
df<-list()
df_m<-list()
df_f<-list()
for (i in 1:length(files)){
  df[[i]]<-fread(paste0("DB_CSV/",files[i]),header=T,sep = ",")
  df_m[[i]]<-fread(paste0("metadata/",folders[i],"_metadata.csv"),header = T,sep = ",")
  df_f[[i]]<-inner_join(df[[i]],df_m[[i]],by="Well")
  write.csv(df_f[[i]],paste("DB_CSV/Merged",files[i],sep = "/"),row.names = F)
}
# Merged all batches of plates:
df<- do.call(rbind,df_f)
write.csv(df,"Raw_merged.csv",row.names = F)


