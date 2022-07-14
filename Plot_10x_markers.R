options(repos = list(CRAN="http://cran.rstudio.com/"))


list.of.packages <- c(c("Seurat","cowplot","ggplot2","ggpubr"))
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}

library(Seurat)
library(cowplot)
library(ggplot2)
library(ggpubr)


features <- unlist(strsplit(commandArgs(trailingOnly = T),","))

#plotting function
set_theme_featureplot <- function(plot,markers,numcol = 3,length = F,names,dir,fill = '#9d9c9c',reverse = F){
  numrow <- ceiling(length(markers) / numcol)
  plot_list <- list()
  for(i in 1:length(markers)){
    plot2 <- plot[[i]]+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))+ggtitle(names[i])+labs(x= "UMAP 1", y = "UMAP 2")
    plot2 <- plot2+theme(panel.background = element_rect(fill = fill),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if(reverse){
      plot2 <- plot2 + scale_x_reverse()
    }
    plot_list[[i]]<- plot2 
    print(plot2)
  } 
  return(ggarrange(plotlist = plot_list,ncol = numcol,nrow = numrow,align = "hv"))
}


#Plot without N2B27 group
load("Seurat_merged_no_N2b27_norm_cc_count_mito.Robj")

pdf("UMAP_wo_N2B27.pdf",width = 7,height = 5)
DimPlot(Seurat_merged,reduction="umap",label=F)+labs(x="UMAP 1",y = "UMAP 2")+theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+scale_x_reverse()
dev.off()

set_theme_featureplot(FeaturePlot(Seurat_merged,order=T,features = features ,cols = c("black","yellow")),markers = seq(length(features)),names = features,numcol = 2,fill = "#ffffff",reverse = T )
ggsave("UMAP_markers_wo_N2B27.pdf",width = 9, height = 2.5 * length(features))


#Plot with N2B27 group
load("Seurat_cc_mito_count_all_samples_no12h.Robj")

pdf("UMAP_with_N2B27.pdf",width = 7,height = 5)
DimPlot(Seurat_merged,reduction="umap",label=F)+labs(x="UMAP 1",y = "UMAP 2")+theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dev.off()

set_theme_featureplot(FeaturePlot(Seurat_merged,order=T,features = features ,cols = c("black","yellow")),markers = seq(length(features)),names = features,numcol = 2,fill = "#ffffff")
ggsave("UMAP_markers_with_N2B27.pdf",width = 9, height = 2.5 * length(features))

