library(pheatmap)
setwd("C:/Users/David/Documents/GitHub/biosynthetic_network_robustness/Figures/Figure5")

x = read.csv("Figure5_Data_PM.csv",check.names=FALSE)

C = c(312,309,310,311,361,253,419,418,252,416,417,420,421,47,48,308,457)+1

x1 = x[,C]
rownames(x1)=x$metabolite[]

y = read.csv("Figure5_Data_Metadata.csv",check.names=FALSE)
y1 = y[C-1,c(4,3,1)]
RN = y$`GENRE Filename`
rownames(y1)=RN[C-1]

L = list("Class" = c("Methanobacteria"="orangered",
                     "Gammaproteobacteria"="navyblue",
                     "Betaproteobacteria"="royalblue4",
                     "Alphaproteobacteria"="royalblue3",
                     "Epsilonproteobacteria"="royalblue2",
                     "Deltaproteobacteria"="royalblue1",
                     "Bacteroidia"="honeydew3",
                     "Flavobacteriia"="honeydew2",
                     "Spirochaetia"="violetred2",
                     "Fusobacteriia"="mediumorchid3",
                     "Clostridia"="darkslategray",
                     "Negativicutes"="darkslategray4",
                     "Bacilli"="darkslategray3",
                     "Mollicutes"="darkslategray2",
                     "Erysipelotrichia"="darkslategray1",
                     "Actinobacteria"="green4",
                     "Coriobacteriia"="green3",
                     "Synergistia"="red2",
                     "Anaerolineae"="indianred4",
                     "TM7 [C-1]"="yellow2",
                     "SR1 [C-1]"="orange2",
                     "Chlamydiia"="hotpink2",
                     "Avg"="white"),
         
         "Phylum" = c("Euryarchaeota"="orangered",
                      "Proteobacteria"="royalblue1",
                      "Bacteroidetes"="honeydew3",
                      "Spirochaetes"="violetred2",
                      "Fusobacteria"="mediumorchid3",
                      "Firmicutes"="darkslategray1",
                      "Actinobacteria"="green4",
                      "Synergistetes"="red2",
                      "Chloroflexi"="indianred4",
                      "TM7 Saccharibacteria"="yellow2",
                      "SR1"="orange2",
                      "Chlamydiae"="hotpink2",
                      "Avg"="white"),
         
         "Gram Stain"=c("Negative"="gray33",
                        "Positive"="purple2",
                        "Avg"="white"))


CL = paste(y$Genus,y$Species,y$`Culture Collection (Strain)`)
CL[457] = "Average oral MB"
ColumnLabels = CL[C-1]

colfunc <- colorRampPalette(c("white","firebrick4"))
x2 = pheatmap(x1,
              cellwidth = 12,
              cellheight = 6,
              annotation_colors = L,
              annotation_col = y1,
              scale="none",
              color = colfunc(100),
              filename="figure_A.pdf",
              show_colnames = TRUE,
              labels_col = ColumnLabels,
              labels_row = rownames(x1),
              fontsize_col = 6,
              fontsize_row = 6,
              clustering_method = "average",
              clustering_distance_rows = "euclidean",
              clustering_distance_cols = "euclidean",
              width=17,
              height=22,
              cluster_cols=T,
              cluster_rows=T)

#colnames(x1)[x2$tree_col$order]

