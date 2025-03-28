library(rgl) 
library(geomorph)
library(Morpho)
library(RColorBrewer)
library(vegan)
library(lattice)
library(RRphylo)
library(phytools)
library(ggplot2)
library(dplyr)
library(coin)
library(rcompanion)
library(maps)
library(mapdata)
library(patchwork)
library(MetBrewer)
library(caper)
library(convevol)
library(cluster)
library(factoextra)
library(cli)
library(dendextend) 
library(car)
library(amap)

setwd("C:/Users/...")

###DATA GEOMORPHO

#Coordinates
Data <- read.csv("Coords.csv", sep=";")
Data <- Data[,-1] 
Data <- arrayspecs(Data, 24, 3)

#Select missing landmarks
Data[1,,12] <- Data[2,,12] <- Data[3,,12]<- Data[4,,12]<- Data[5,,12]<- Data[7,,12] <- Data[9,,12]<- Data[10,,12]<- Data[13,,12] <- NA  

Data[1,,13]<- Data[2,,13]<- Data[3,,13]<- Data[4,,13]<- Data[5,,13]<- Data[6,,13]<- Data[7,,13]<- Data[20,,13]<- Data[21,,13]<- Data[22,,13]<- Data[23,,13]<- Data[24,,13] <- NA

Data[7,,24]<- Data[13,,24]<- NA 

Data[1,,26]<- Data[2,,26]<- Data[3,,26]<- Data[4,,26]<- Data[5,,26]<- Data[7,,26]<- Data[15,,26]<- Data[16,,26]<- Data[17,,26]<- Data[18,,26]<- Data[19,,26]<- Data[20,,26]<- Data[21,,26]<- Data[22,,26]<- Data[23,,26]<- Data[24,,26]<- NA
#Estimate missing landmarks
estimate.missing(Data, method = "TPS")
coordData <- estimate.missing(Data,method = "TPS")

#sliding landmarks file
curve<-as.matrix(read.table(file = "Curve.csv", header = F, sep=";", skip = 0, fill = TRUE)) 

#Names
Indivlist <- read.table(file = "ID.csv", header = F, sep=";", skip = 0, fill = TRUE)
dimnames(coordData)[3] <- Indivlist	
#Data creation
data<- list("land" = coordData,"curve"=curve,"ID"=Indivlist) 
str(data) 


###PROCRUSTES
R.gpa <- gpagen(coordData,curves=curve, ProcD = T)
gpa_2d <- two.d.array(R.gpa$coords) #GPA coords
row.names(gpa_2d) <- Indivlist[,1]


###MEAN VALUES PER TAXA
#Taxa
Taxa <- read.csv("Taxa.csv", sep=";", header=F)
rum2d_sp <- cbind.data.frame(Taxa[,1], gpa_2d)
split_sp <- split(rum2d_sp, rum2d_sp[,1])
row.names(rum2d_sp) <- Indivlist[,1]

#mean
for (i in 1:25)
{
  split_sp[[i]] <- split_sp[[i]][,2:73]
  df <- list(unique(Taxa[,1]))
  assign(paste("sp",i, sep=""), as.matrix(split_sp[[i]]))
  
  if (nrow(get(paste("sp",i, sep=""))) > 1) {
    assign(paste("df",i, sep=""), arrayspecs(get(paste0("sp", i)), 24, 3))
    assign(paste("df",i, sep=""), mshape(get(paste0("df", i))))
  }		
  
  else	{	
    assign(paste("df",i, sep=""), matrix(get(paste0("sp", i)), ncol=72))
    assign(paste("df",i, sep=""), arrayspecs(get(paste0("df", i)), 24, 3))
  }		
  
}

#New dataset with one specimen per species
data_25 <- bindArr(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,
                   df15,df16,df17,df18,df19,df20,df21,df22,df23,df24,df25, along =3)
#Clades attribution
clades <- read.table(file = "Family.csv", header = F, sep=";", skip = 0, fill = TRUE)  
clades <- unlist(clades, use.names=FALSE)

#data mean creation
rm_names <- as.data.frame(ls(split_sp))
dimnames(data_25)[3] <- rm_names
data_mean <- list("land" = data_25, "curve"=curve,"species"=rm_names,"clades"=clades)



###PROCRUSTES MEAN
R.gpa2 <- gpagen(data_25,curves=curve, ProcD = T)
names(R.gpa2$Csize) <- rm_names[,1]
gpa_m<- two.d.array(R.gpa2$coords)
row.names(gpa_m) <- rm_names[,1]



###PGLS SHAPE ~ CENTROID SIZE

#Tree
tree <- read.tree("tree_GM.TRE") #import tree
tree <- multi2di(tree)
tree$tip.label <- gsub("'","",tree$tip.label)
tree
gdf <- geomorph.data.frame(gpa_m, phy = tree)
match(tree$tip.label,rownames(gpa_m))
#Test
Size.pgls <- procD.pgls(R.gpa2$coords ~ R.gpa2$Csize, phy = tree, data = R.gpa2, iter = 999)
summary(Size.pgls)



###PHYLO SIGNAL

#Pagel's_lambda
phylosig(tree,gpa_m, method="lambda",test=TRUE)
#K_Blomberg
phylosig(tree,gpa_m, method="K",test=TRUE)



###PCA
pca <- gm.prcomp(R.gpa2$coord)
row.names(pca$x) <- rm_names[,1]
orpdata_exp <- pca$x

eigenvalues <- pca$sdev^2  
# ploteigenvalues
plot(eigenvalues, type = "b", pch = 19, col = "blue", 
     xlab = "Component number", 
     ylab = "Eigenvalues",
     main = "Extant")



###PCA Visualisation

##Plot
#Plot Legend 
gp<-as.factor(unlist(clades))
names(gp) <-row.names(clades)
col.gp <- col.gp1 <- c("#99610a","#67322e","#175449", "#c38f16","#6e948c", "red")
names(col.gp) <- levels(gp)
col.gp <- col.gp[match(gp, names(col.gp))]
##Legend shape
shp.gp <- c(24,21,23,24,22,21)
names(shp.gp) <- levels(gp)
shp.gp <- shp.gp[match(gp, names(shp.gp))]

#Phylomorphospace PC1 PC2
phylomorphospace(tree,pca$x[,1:2],node.col=col.gp,bg=col.gp,bty="l",ftype="off",node.by.map=T,
                 node.size=c(0,1.2),xlab="PC1 (26.16%)",ylab="PC2 (15.37%)")
text(pca$x[,1], pca$x[,2], labels = rownames(pca$x), cex = 0.7, pos = 4, col = col.gp)


#Phylomorphospace PC3 PC4
phylomorphospace(tree,pca$x[,3:4],node.col=col.gp,bg=col.gp,bty="l",ftype="off",node.by.map=T,
                 node.size=c(0,1.2),xlab="PC3 (10.49%)",ylab="PC4 (8.23%)")
text(pca$x[,3], pca$x[,4], labels = rownames(pca$x), cex = 0.7, pos = 4, col = col.gp)


##Shape visualisation

#PC1
plotRefToTarget(pca$shapes$shapes.comp1$min, pca$shapes$shapes.comp1$max,
                gridPars = gridPar(link.lwd=0.5,pt.bg=0.001,tar.pt.size=1.5,pt.size=1,tar.pt.bg = "black", tar.link.col="black",
                                   tar.link.lwd = 10), method = "points", links = curve)
rgl.snapshot('C:/Users/.../PC1max_.png')

#PC2
plotRefToTarget(pca$shapes$shapes.comp2$min, pca$shapes$shapes.comp2$max,
                gridPars = gridPar(link.lwd=0.5,pt.bg=0.001,tar.pt.size=1.5,pt.size=1,tar.pt.bg = "black", tar.link.col="black",
                                   tar.link.lwd = 10), method = "points", links = curve)
rgl.snapshot('C:/Users/.../PC2max_.png')







##DATA TABLE
BVBM <- read.table(file = "BVBM.csv", header = T, sep=";", skip = 0, fill = TRUE)


## EQ_Orliac

# plot
plot_base <- ggplot() +
  theme_classic() +
  labs(y = "EQ Orliac", x = "Family", title = "Encephalization quotient values per Family (Extant + Extinct)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

hoplitomeryx_mean_EQ <- mean(BVBM %>% filter(Family == "Hoplitomeryx") %>% pull(EQ_Orliac), na.rm = TRUE)

plot_combined <- plot_base +
  geom_boxplot(data = BVBM %>% filter(!Family %in% c("Hoplitomeryx")), 
               aes(x = Family, y = EQ_Orliac, fill = Family, color = Family),  # Associer `color` à `Family`
               outlier.shape = NA, alpha = 0.5, size=0.8,show.legend = FALSE) +  
  geom_point(data = BVBM %>% filter(Extant_Fossil == "Extant" & Family != "Hoplitomeryx"), 
             aes(x = Family, y = EQ_Orliac, color = Family, shape = "Extant"), 
             size = 2.5, show.legend = FALSE) +  
  geom_point(data = BVBM %>% filter(Extant_Fossil == "Fossil" & Family != "Hoplitomeryx"), 
             aes(x = Family, y = EQ_Orliac, color = Family, shape = "Fossil"), 
             size = 2.5, show.legend = FALSE) +  
  geom_text(data = BVBM %>% filter(Species %in% c("Candiacervus_ropalophorus", "Dicrocerus_sp", "Samotherium_sp", 
                                                  "Antifer_ensenadensis", "Myotragus_balearicus")),
            aes(x = Family, y = EQ_Orliac, label = Species), 
            hjust = 0.5, vjust = 1.5, size = 3) +  
  geom_hline(yintercept = hoplitomeryx_mean_EQ, 
             color = "black", linetype = "dashed", size = 1) + 
  scale_fill_manual(values = c("#99610a", "#67322e", "#175449", "#c38f16", "#6e948c", "#122c43","red")) +  
  scale_color_manual(values = c("#99610a", "#67322e", "#175449", "#c38f16", "#6e948c", "#122c43","red")) +  
  scale_shape_manual(name = "Type", values = c("Extant" = 16, "Fossil" = 15), 
                     labels = c("Extant" = "Extant", "Fossil" = "Extinct")) + 
  theme(legend.position = "right") +  
  guides(color = "none") +  
  scale_x_discrete(limits = c("Bovidae", "Moschidae", "Cervidae", "Giraffidae", "Antilocapridae","Tragulidae"))  

print(plot_combined)

#One-way test
BVBM$Family <- as.factor(BVBM$Family)
oneway_test(EQ_Orliac~Family,data=BVBM) 

#paiwise Permutation test
PT_EQ<-pairwisePermutationTest(EQ_Orliac~Family,data=BVBM,method="fdr")
PT_EQ



##RAEC

plot_base4 <- ggplot() +
  theme_classic() +
  labs(y = "RAEC (%)", x = "Family", title = "RAEC values per Family (Extant + Extinct)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_combined4 <- plot_base4 +
  geom_boxplot(data = BVBM %>% filter(!Family %in% c("Hoplitomeryx")), 
               aes(x = Family, y = RAEC, fill = Family, color = Family),  # Bordure = remplissage
               outlier.shape = NA, alpha = 0.5, size = 0.8, show.legend = FALSE) +  # Taille des bordures ajustée
  geom_point(data = BVBM %>% filter(Extant_Fossil == "Extant"), 
             aes(x = Family, y = RAEC, color = Family, shape = "Extant"), 
             size = 2.5, show.legend = FALSE) +  
  geom_point(data = BVBM %>% filter(Extant_Fossil == "Fossil"), 
             aes(x = Family, y = RAEC, color = Family, shape = "Fossil"), 
             size = 2.5, show.legend = FALSE) +  
  geom_text(data = BVBM %>% filter(Species %in% c("Antifer_ensenadensis", "Myotragus_balearicus")),
            aes(x = Family, y = RAEC, label = Species), 
            hjust = 0.5, vjust = 1.5, size = 3) +  
  scale_fill_manual(values = c("#99610a", "#67322e", "#175449","#c38f16","red","#6e948c","#122c43")) +  
  scale_color_manual(values = c("#99610a", "#67322e", "#175449","#c38f16","red","#6e948c","#122c43")) +   
  scale_shape_manual(name = "Type", values = c("Extant" = 16, "Fossil" = 15), 
                     labels = c("Extant" = "Extant", "Fossil" = "Extinct")) +  
  theme(legend.position = "right") +  
  guides(color = "none") +  
  scale_x_discrete(limits = c("Bovidae", "Moschidae", "Cervidae","Antilocapridae","Tragulidae")) +
  geom_hline(data = BVBM %>% filter(Family == "Hoplitomeryx"), 
             aes(yintercept = mean(RAEC, na.rm = TRUE), color = "Hoplitomeryx"),  # Ligne associée à Hoplitomeryx
             linetype = "dashed", size = 1)  

print(plot_combined4)

#One-way test
BVBM$Family<-as.factor(BVBM$Family)
oneway_test(RAEC~Family,data=BVBM) 

#paiwise Permutation test
PT_RAEC<-pairwisePermutationTest(RAEC~Family,data=BVBM,method="fdr")
PT_RAEC




## Olfactory bulbs Volume Ratio

plot_base5 <- ggplot() +
  theme_classic() +
  labs(y = "Olfcatory bulbs volume vs Endocast volume (%)", x = "Family", title = "Olfcatory bulbs volume ratio per Family (Extant + Extinct)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_combined5 <- plot_base5 +
  geom_boxplot(data = BVBM %>% filter(!Family %in% c("Hoplitomeryx")), 
               aes(x = Family, y = ObV_EV, fill = Family, color = Family),  # Bordure = remplissage
               outlier.shape = NA, alpha = 0.5, size = 0.8, show.legend = FALSE) +  # Taille des bordures ajustée
  geom_point(data = BVBM %>% filter(Extant_Fossil == "Extant"), 
             aes(x = Family, y = ObV_EV, color = Family, shape = "Extant"), 
             size = 2.5, show.legend = FALSE) +  
  geom_point(data = BVBM %>% filter(Extant_Fossil == "Fossil"), 
             aes(x = Family, y = ObV_EV, color = Family, shape = "Fossil"), 
             size = 2.5, show.legend = FALSE) +  
  geom_text(data = BVBM %>% filter(Species %in% c("Antifer_ensenadensis", "Myotragus_balearicus","Kobus_ellippsiprymnus")),
            aes(x = Family, y = ObV_EV, label = Species), 
            hjust = 0.5, vjust = 1.5, size = 3) +  
  scale_fill_manual(values = c("#99610a", "#67322e", "#175449","#c38f16","#122c43","#6e948c","red")) +  
  scale_color_manual(values = c("#99610a", "#67322e", "#175449","#c38f16","red","#6e948c","#122c43")) +   
  scale_shape_manual(name = "Type", values = c("Extant" = 16, "Fossil" = 15), 
                     labels = c("Extant" = "Extant", "Fossil" = "Extinct")) +  
  theme(legend.position = "right") +  
  guides(color = "none") +  
  scale_x_discrete(limits = c("Bovidae", "Moschidae", "Cervidae","Antilocapridae","Tragulidae")) +
  geom_hline(data = BVBM %>% filter(Family == "Hoplitomeryx"), 
             aes(yintercept = mean(ObV_EV, na.rm = TRUE), color = "Hoplitomeryx"),  # Ligne associée à Hoplitomeryx
             linetype = "dashed", size = 1)  

print(plot_combined5)

#One-way test
oneway_test(ObV_EV~Family,data=BVBM) #not significant




##BRM VS BM: Régression LM

model <- lm(Log_BrM ~ Log_BM, data=BVBM)
summary(model)



##BRM VS BM: PGLS 

#data cleaning
BVBM_clean <- BVBM %>%
  filter(!is.na(Log_BM) & !is.na(Log_BrM)) %>%
  dplyr::select(-Neo_Ratio, -ObS_ES)

#tree
tree2 <- read.tree("tree_BVBM.TRE")
tree2 <- multi2di(tree2)
tree2$tip.label <- gsub("'", "", tree2$tip.label)
plot(tree2)
# Prune the tree to match cleaned dataset
tree2_clean <- drop.tip(tree2, setdiff(tree2$tip.label, BVBM_clean$Species))

# Row names of BVBM_clean match the tip labels of tree2_clean
BVBM_clean <- BVBM_clean[BVBM_clean$Species %in% tree2_clean$tip.label, ]

# BrM VS BM
PGLS_BVBM <- comparative.data(tree2_clean, BVBM_clean, Species, vcv = TRUE, vcv.dim = 3)
PGLSmod1 <- pgls(Log_BrM ~ Log_BM, data = PGLS_BVBM, lambda = 'ML', delta = 'ML')
summary(PGLSmod1)

# coeff models
intercept1 <- coef(model)[1]  # Intercept model LM
slope1 <- coef(model)[2]      # Slope model LM
intercept2 <- coef(PGLSmod1)[1]  # Intercept model PGLS
slope2 <- coef(PGLSmod1)[2]      # Slope medel PGLS

# R² 
r_squared_LM <- summary(model)$r.squared
r_squared_PGLS <- summary(PGLSmod1)$r.squared



##BRM VS BM: plot

ggplot(BVBM_clean, aes(x = Log_BM, y = Log_BrM)) +
  geom_point(aes(color = Family, shape = Extant_Fossil), size = 3) + 
  geom_abline(intercept = intercept1, slope = slope1, color = "black", linetype = "dashed", size = 1) +  # Ligne LM
  geom_abline(intercept = intercept2, slope = slope2, color = "darkgray", linetype = "dashed", size = 1) +    # Ligne PGLS
  scale_shape_manual(values = c(19, 15)) + 
  scale_fill_manual(values = c("#99610a", "#67322e", "#175449", "#c38f16", "red", "#6e948c", "#122c43")) +
  scale_color_manual(values = c("#99610a", "#67322e", "#175449", "#c38f16", "red", "#6e948c", "#122c43")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = 'bold')) +
  labs(x = "Log Body Mass (g)", y = "Log Brain Mass (g)", title = "Brain Mass vs Body Mass with LM and PGLS Fits") +
  annotate("text", x = 5.5, y = 3, 
           label = paste("LM Slope:", round(slope1, 3), "| R²:", round(r_squared_LM, 3), "\n",
                         "PGLS Slope:", round(slope2, 3), "| R²:", round(r_squared_PGLS, 3)), 
           hjust = 1, vjust = 1, size = 4, color = "black", fontface = "bold")


