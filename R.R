##BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(ComplexHeatmap)

library(ggplot2)
library(showtext)
library(Hmisc)
windowsFonts(myFont = windowsFont("Arial"))

library(RColorBrewer)
library(circlize)
library(arulesViz) 


library(tidyverse)

library(data.table)



normalizepath <- function(xma,ycls) {
  
  dffdrmetma=data.frame()
  
  xmahe=cbind(t(xma),ycls)
  xmahe=data.frame(xmahe)
  

  
  
  xmaheg=split(xmahe,xmahe$ycls)
  
  print(dim(xmaheg$`1`))###2型
  print(dim(xmaheg$`2`))###1型
  print(dim(xmaheg$`3`))###3型
  # # [1] 26 40
  # # [1] 13 40
  # # [1] 21 40
  

  
  
  print(c(xmahe$ycls))
  clst=unique(xmahe$ycls, fromLast = TRUE)
  print(clst)
  
  # for ( i in clst){
  #   print(i)
  #   i=str(i)
  #   print(rbind(apply(xmaheg$`i`, 2, mean),dffdrmetma,make.row.names=T,make.names=T))
  #   # dffdrmetma=rbind(dffdrmetma,apply(xmaheg$`i`, 2, mean))
  #   
  # }
  dffdrmetma=rbind(apply(xmaheg$`2`, 2, mean),apply(xmaheg$`1`, 2, mean),apply(xmaheg$`3`, 2, mean))
  print(dffdrmetma)
  
  

  return(data.frame(dffdrmetma))
  

  
}









# fr="sigheatmap.txt"
# fr="sigheatmap260.txt"
fr="sigheatmap25.txt"
fr="sigheatmap23.txt"
wr=paste(strsplit(fr,split='.gct'),".png")
wr=paste(strsplit(fr,split='.gct'),".pdf")





dxy= read.table(fr, header = TRUE, sep = "\t",row.names = 1,quote = "")

print(dxy)

# genelist=c("P12004", "P23246", "P11387", "Q15185", "Q14978", "P11413", "P78527", "P06400", "P14866", "Q99873", "Q9NPL8", "Q9NY12", "P05386", "P38159")

hang=dim(dxy)[1]-2

dtt=dxy[c(1:hang),]
print(dtt)

dxyt=data.frame(t(dxy))
print(dxyt)
tisu=dxyt$sample
subtypetr=dxyt$subtype



print(tisu)

print(subtypetr)


# nsam=60
# 
# 
# nc=c(t(dxy[1,]))
# 
# samt=nc[(length(nc)-60+1):length(nc)]
# print(length(samt))
# print(samt)
# samt=as.integer(samt)
# 
# 
# class=rep(samt, each =1, times = 1)

# datasub3=data[,grep("1",data[2,])]##预后最差的组22
#
#
# datasub3$V2293=gsub("3","III",datasub3$V2293)
#
# print(dim(datasub3))
#
#
# datasub2=subset(data,V2293 ==1)##预后较中的组26
#
# datasub2$V2293=gsub("1","II",datasub2$V2293)
#
# print(dim(datasub2))
#
#
#
# datasub1=subset(data,V2293 ==2)##预后最好的组12
#
# datasub1$V2293=gsub("2","I",datasub1$V2293)
#
# print(dim(datasub1))
# newbind=rbind(datasub1,datasub2,datasub3)









# dx=dxy[c(grep("^HALLMARK",row.names(dxy))),]
# 
# # print(dx)
# 
# fdr=dx[,c(grep("^fdr",names(dx)))]
# 
# 
# fdr=as.data.frame(lapply(fdr,as.numeric),row.names=row.names(fdr))
# 
# 
# print(dim(fdr))
# 
# 
# 
# 
# 
# nes=dx[,c(grep("^S[0-9]",names(dx)))]
# # nes=as.data.frame(lapply(nes,as.numeric),optional = FALSE)
# 
# 
# nes=as.data.frame(lapply(nes,as.numeric),row.names=row.names(nes))
# 
# print(dim(nes))
# 
# # write.table(data.frame(nes),"nes.txt",col.names=TRUE, row.names=TRUE,sep = "\t",quote=FALSE)##CHEK 前716列是T的原始表达值，中716列是FC的原始表达值，后716列是FC的中位数归一化表达值
# # 
# 
# 
# nesguoyi = scale(nes, center = TRUE, scale = TRUE)##按每列归一化 zscore
# # nesguoyi=nes
# 
# print(nesguoyi)
# 
# nesguoyiaddhang =sweep(nesguoyi,1, apply(nesguoyi,1,mean,na.rm=T))#减去每行的中位数（na.rm=T表示求中位数时不计入缺失值）
# 
# 
# write.table(data.frame(nesguoyiaddhang),"nesguoyiaddhang.txt",col.names=TRUE, row.names=TRUE,sep = "\t",quote=FALSE)##CHEK 前716列是T的原始表达值，中716列是FC的原始表达值，后716列是FC的中位数归一化表达值
# 
# 
# print(nesguoyiaddhang)
# 
# print(dim(nesguoyiaddhang))



# fdrhe=as.data.table(cbind(t(fdr),class))
# 
# print(fdrhe)
# 
# # fdrhe[,lapply(.SD,mean), .SDcols=c("class","V39")]
# 
# fdrhe[,lapply(.SD,sum),by=V2]





# 
# fdrhe=cbind(t(fdr),class)
# fdrhe=data.frame(fdrhe)
# 
# print(fdrhe)
# 
# 
# fdrheg=split(fdrhe,fdrhe$class)
# 
# print(dim(fdrheg$`1`))###2型
# print(dim(fdrheg$`2`))###1型
# print(dim(fdrheg$`3`))###3型
# # # [1] 26 40
# # # [1] 13 40
# # # [1] 21 40
# 
# print(apply(fdrheg$`2`, 2, mean))
# print(apply(fdrheg$`1`, 2, mean))
# print(apply(fdrheg$`3`, 2, mean))
# 
# 
# dffdrmean=rbind(apply(fdrheg$`2`, 2, mean),apply(fdrheg$`1`, 2, mean),apply(fdrheg$`3`, 2, mean))
# print(dffdrmean)








# 
# fdrhe1=fdrhe %>% select(contains("1"))##starts_with
# print(names(fdrhe1))
# 
# fdrhe11=fdrhe1[c(-1),]
# 
# print(fdrhe11)
# fdrhe1mean=rowMeans(fdrhe11,na.rm = T)
# print(fdrhe1mean)
# 
# 
# print(11111111111)
# 
# 
# 
# 
# fdrhe2=fdrhe %>% select(contains("2"))##starts_with
# print(names(fdrhe2))
# 
# fdrhe3=fdrhe %>% select(contains("3"))##starts_with
# print(names(fdrhe3))
# 
# 
# neshe1=neshe %>% select(contains("1"))##starts_with
# print(names(neshe1))
# 
# 
# neshe2=neshe %>% select(contains("2"))##starts_with
# print(names(neshe2))
# 
# neshe3=neshe %>% select(contains("3"))##starts_with
# print(names(neshe3))
# 
# 
# 






#
# # d = sweep(dratio,1, apply(dratio,1,median,na.rm=T))#
#
#
#


# dt<-cbind(fdr,nes)
# print(dim(dt))
# print(dt)

# write.table(data.frame(dt),"dt.txt",col.names=TRUE, row.names=TRUE,sep = "\t",quote=FALSE)##CHEK 前716列是T的原始表达值，中716列是FC的原始表达值，后716列是FC的中位数归一化表达值
# cl=dt[c("class"),]
# print(cl)
#


# fdrheg=normalizepath(fdr,class)
# print(fdrheg)
# 
# nesguoyiaddhangg=normalizepath(nesguoyiaddhang,class)
# print(nesguoyiaddhangg)




# print(subset (nesguoyiaddhangg, select = ycls))


# classtr=nesguoyiaddhangg$ycls
# print(classtr)
# dtt<-as.matrix(t(subset (nesguoyiaddhangg, select = -ycls)))
# 
# fdrtt=as.matrix(t(subset (fdrheg, select = -ycls)))
# print(dtt)
# 
# print(classtr)
# 
# print(fdrtt)
#
#
# ha1 = HeatmapAnnotation(Subtype=nesguoyiaddhang$class,col = list( Subtype = c("0" = "#13FC00","1"="#0037FE","2"="#FE0000","3"="LightSalmon","4"="grey","5"="red")),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11))##绿；蓝；红
# ha1 = HeatmapAnnotation(Subtype=classtr,col = list( Subtype = c("1"="#0037FE","2"="#FE0000","3"="LightSalmon")),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11))##绿；蓝；红
# 
# 
# 
# # ha2 = HeatmapAnnotation("Lymphatic metastasis"=newbind$Details.of.lymph.node.involvement, col = list("Lymphatic metastasis" = c("Non metastasis"="Gainsboro","Metastasis"="#FE0000")), annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11),gp = gpar(col='white'),show_legend=T  )
# # ha2 = HeatmapAnnotation(Subtype2=class,col = list( Subtype2 = c("0" = "#13FC00","1"="#0037FE","2"="#FE0000","3"="LightSalmon","4"="grey","5"="red")),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11))##绿；蓝；红
# ha2 = HeatmapAnnotation(Subtype=classtr,col = list( Subtype = c("1"="#0037FE","2"="#FE0000","3"="LightSalmon")),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11))##绿；蓝；红

# ha = HeatmapAnnotation(Subtype=classtr,col = list( Subtype = c("2"="#0037FE","1"="LightSalmon","3"="#FE0000")),annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 11))##绿；蓝；红
###0037FE蓝色  #FE0000 红色

#
#

# ha1 = HeatmapAnnotation(Sample=tisu,col = list( Subtype = c("P" = "#6A7BBB","T"="#F05E6B")),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11))##绿；蓝；红
ha1 = HeatmapAnnotation(Sample=tisu,col = list( Sample =  c("P" = "#84A5D6","T"="#F48B8C")),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11))##绿；蓝；红

ha2 = HeatmapAnnotation(Subtype=subtypetr,col = list( Subtype = c("S-I" = "#6A7BBB","S-II"="#FFAC1D","S-III"="#F05E6B")),annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11))##绿；蓝；红


ha = c(ha1, ha2,gap = unit(0.5, "mm"))
print(ha)

#
# col_fun = colorRamp2(c(-1, 0, 1), c("#0037FE", "white", "#FE0000"))
#
# col_fun(seq(-1,  1))

# col_fun = colorRamp2(c(0.5, 0, 2), c("#4184BB", "white", "red"))
col_fun = colorRamp2(c(0.5, 1, 1.5), c("skyblue", "white", "red"))

col_fun(seq(0.5,  1.5))

print(222222222221)
# pa=cluster::pam(dtt,k=3)

# small_mat = mat[1:50, 1:18]
# print(small_mat)




dtt<-as.matrix(as.data.frame(lapply(dtt,as.numeric),row.names=row.names(dtt)))
print(dtt)

ht_list2=Heatmap(dtt, name = "Z-score(Protein FC-value)",col=col_fun,rect_gp = gpar(col = "white", lwd = 2),
                 # annotation_legend_param  = list(
                 #   name=list(
                 #   at = c(-2, 0, 2),
                 #   labels = c("low", "zero", "high"),
                 #   title = "Z-score(Protein FC-value)",
                 #   legend_height = unit(4, "cm"),
                 #   title_position = "lefttop-rot"
                 # )),
                 heatmap_legend_param = list(
                   at = c(0.5,  1.5),
                   labels = c("Low",  "High"),
                   title = "Relative Protein abundance",
                   legend_height = unit(4, "cm"),
                   title_position = "lefttop-rot"
                 ),

                 # row_split = paste("pam",pa$clustering),

                 top_annotation =ha,
                 # right_annotation = anomaker,






                                  # cell_fun = function(j, i, x, y, width, height, fill) {
                                  #   if( fdrtt< 0.01)
                                  # 
                                  #     grid.points( x, y, pch = "*", size = unit(4, "mm"))
                                  # },       ###可以
                                  # 


                 row_names_max_width = max_text_width(
                   rownames(dtt),
                   gp = gpar(fontsize = 12)),



                 column_names_rot = 75,column_names_gp=gpar(fontsize=8),cluster_rows=F,cluster_columns = F,show_row_names=T,show_column_names =F,row_dend_reorder = TRUE,

                 # column_names_rot = 75,column_names_gp=gpar(fontsize=8),cluster_rows=row_dend,cluster_columns = F,show_row_names=F,show_column_names = T,row_dend_reorder = T,





                 clustering_distance_rows = "spearman",show_row_dend=F, show_column_dend=F, row_dend_width = unit(13, "cm"),column_dend_height = unit(7, "cm")
                 )

print(111111111)

# draw(ht_list2)

# pdf(wr,width = 30.56/2.54,height = 90.56/2.54)
pdf(wr,width = 22.56/2.54,height = 10.56/2.54)
draw(ht_list2,merge_legend = TRUE,heatmap_legend_side = "left",
     annotation_legend_side = "left")
print(111111111)
dev.off()
