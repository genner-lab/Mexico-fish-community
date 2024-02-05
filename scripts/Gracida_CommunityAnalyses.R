#Packages required

library(vegan)
library(ape)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggokabeito)
library(jtools)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(magrittr)
library(indicspecies)
library(reshape2)
library(data.table)

#Analysis 1: Rarefaction of species diversity

Gracida_Species96_t <- read.table("Gracida_Species96_t.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#Analysis 1: Make the six accumulation curves, using raw capture data

pdf("Noh-Bec_plot.pdf")
sp_raw_NohBec <- specaccum(Gracida_Species96_t[2:17], "random")
NohBec_plot <- plot(sp_raw_NohBec, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw_NohBec, col="brown3", add=TRUE)
dev.off()

pdf("Petcacab_plot.pdf")
sp_raw_Petcacab <- specaccum(Gracida_Species96_t[18:33], "random")
Petcacab_plot <- plot(sp_raw_Petcacab, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw_Petcacab, col="brown3", add=TRUE)
dev.off()

pdf("Santa_Teresita_plot.pdf")
sp_raw_Santa_Teresita <- specaccum(Gracida_Species96_t[34:49], "random")
Santa_Teresita_plot <- plot(sp_raw_Santa_Teresita, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw_Santa_Teresita, col="brown3", add=TRUE)
dev.off()

pdf("San_Felipe_plot.pdf")
sp_raw_San_Felipe <- specaccum(Gracida_Species96_t[50:65], "random")
San_Felipe_plot <- plot(sp_raw_San_Felipe, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw_San_Felipe, col="brown3", add=TRUE)
dev.off()

pdf("Coabas_plot.pdf")
sp_raw_Coabas <- specaccum(Gracida_Species96_t[66:80], "random")
Coabas_plot <- plot(sp_raw_Coabas, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw_Coabas, col="brown3", add=TRUE)
dev.off()

pdf("Chacambacab_plot.pdf")
sp_raw_Chacambacab <- specaccum(Gracida_Species96_t[81:94], "random")
Chacambacab_plot <- plot(sp_raw_Chacambacab, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw_Chacambacab, col="brown3", add=TRUE)
dev.off()

#Analysis 2: Principal Coordinates Analysis

Gracida_Species96 <- read.table("Gracida_Species96.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#Analysis 2: Remove samples that yielded no catch (all zero values), and log10 transform the data

Gracida_Species93 <- Gracida_Species96[rowSums(Gracida_Species96[,3:20])>0,] 
log_Gracida_Species93 <- Gracida_Species93
log_Gracida_Species93[,3:20] <- log10((Gracida_Species93[,3:20])+1)

#Analysis 2: Generation distance matrix and use for Principal Coordinates Analysis 

SpeciesMatrixHell.D <- vegdist(log_Gracida_Species93[,3:20], method="hellinger")
PCOA_SpeciesMatrixHell<- pcoa(SpeciesMatrixHell.D, correction="none", rn=NULL)
PCOA_SpeciesMatrixHell

PCOAscoresHell <- PCOA_SpeciesMatrixHell$vectors
PCOAscoresHell <- as.data.frame(PCOAscoresHell)
PCOAscoresHell <- cbind(log_Gracida_Species93[,1:2],PCOAscoresHell)

#Analysis 2: Plot the PCOA, here using a star plot.

colours = c("#999999", "#E69F00", "#56B4E9", "black", "red","#CCCC00")

PlotHell <- ggscatter(PCOAscoresHell, x = "Axis.1", y = "Axis.2", color = "Lake",
                      ellipse = FALSE, ellipse.type = "convex", mean.point = TRUE,
                      star.plot = TRUE)  +
  scale_colour_manual(values = colours) +
  labs(x ="PCOA1 (27.4% of variation)", y = "PCOA2 (14.7% of variation)") +
  theme(legend.position = "right")
PlotHell

#Analysis 3: Permanova (permutation Anova)

GlobalTest <- adonis2(SpeciesMatrixHell.D ~ Lake, permutations = 100000, data = Gracida_Species93)
GlobalTest

PairwiseTest <- pairwise.adonis(SpeciesMatrixHell.D,Gracida_Species93$Lake, perm=100000)
PairwiseTest

#Analysis 4: Indicator Species analysis

Abund <- Gracida_Species93[,3:20]
Gracida_EnvPredic93 <- read.table("Gracida_EnvPredic93.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
Group <- Gracida_EnvPredic93$Lake
indval <- multipatt(Abund, Group, control = how(nperm=999)) 

#Analysis 4: Make table of IndVal scores

IndVal_Index <- indval$A*indval$B*100
IndVal_Index <- IndVal_Index[,1:6]
IndVal_Index <- as.data.frame(IndVal_Index)
IndVal_Index <- tibble::rownames_to_column(IndVal_Index, "Species")
IndVal_Index_Long = melt(IndVal_Index , id = c("Species"))
colnames(IndVal_Index_Long)[2] <- "Habitat"
colnames(IndVal_Index_Long)[3] <- "IndVal"

level_order <- c('Gambusia yucatana',
                 'Mayaheros urophthalmus',
                 'Thorichthys meeki',
                 'Trichromis salvini',
                 'Atherinella alvarezi',
                 'Rhamdia guatemalensis',
                 'Poecilia kykesis',
                 'Vieja melanurus',
                 'Belonesox belizanus',
                 'Cribroheros robertsoni',
                 'Dorosoma petenense',
                 'Parachromis multifasciatus',                
                 'Petenia splendida',
                 'Poecilia mexicana',
                 'Astyanax bacalarensis',
                 'Cryptoheros chetumalensis',
                 'Hyphessobrycon compressus',
                 'Oreochromis sp'
)

level_order2 <- c('Noh-Bec',
                  'Santa Teresita',
                  'San Felipe',
                  'Chacambacab',
                  'Caobas',
                  'Petcacab')

IndValPlot = ggplot(IndVal_Index_Long, aes(x = factor(Habitat, level = level_order2), y = factor(Species, level = level_order))) + 
  geom_point(aes(size = IndVal, fill = Habitat), alpha = 0.7, shape = 21)  +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "black", "red","#CCCC00")) +
  scale_size_continuous(limits = c(0.000000001, 55)) + scale_y_discrete(limits = rev) +
  labs( x= "", y = "", size = "Indicator Value", fill = "") +
  theme_classic() + theme(legend.position='right') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
IndValPlot 

#Analysis 5: Redundancy analysis

#Analysis 5: Load environmental data, for 93 samples with catches

Gracida_EnvPredic93$LogDepth <- log10(Gracida_EnvPredic93$Depth)
Gracida_EnvPredic93$LogArea <- log10(Gracida_EnvPredic93$Area)
Gracida_EnvPredic93$Oreo <- as.factor(Gracida_EnvPredic93$Oreo)
log_Gracida_Species93_Native <- cbind(log_Gracida_Species93[,1:11],log_Gracida_Species93[,13:20])

RDA=rda(log_Gracida_Species93_Native[,3:19] ~ LogDepth+Trans+Shore+pH+DO+Temp+LogArea+Subs+Oreo, Gracida_EnvPredic93)
plot(RDA)
anova(RDA, by="axis",nperms=10000)
anova(RDA, by="margin",nperms=10000)
summary(RDA)

#Analysis 5: Plot Site RDA values

GracidaRDA_Sites <- RDA$CCA$wa
GracidaRDA_Sites <- cbind(Gracida_Species93[1:2], GracidaRDA_Sites)
GracidaRDA_Sites

RDA_SitePlot <- ggscatter(GracidaRDA_Sites, x = "RDA1", y = "RDA2", color = "Lake",
                          ellipse = TRUE, ellipse.type = "confidence", mean.point = FALSE, ellipse.level = 0.95,
                          star.plot = FALSE, ellipse.alpha = 0) +
  scale_x_continuous(limits = c(-0.4, 0.4), breaks =c(-0.4,0,0.4)) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks =c(-0.6,0,0.6)) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue") +
  geom_vline(xintercept=0, linetype="dashed", color = "blue") +
  theme(legend.position = "right") +
  scale_colour_manual(values = colours) +  xlab("RDA Axis 1 (18.1% of variation)") + ylab("RDA Axis 2 (5.7% of variation)")
RDA_SitePlot

#Analysis 5: Place Environmental RDA values in dataframes for plotting

GracidaRDA_Env <- as.data.frame(row.names(RDA$CCA$biplot))
GracidaRDA_Env <- cbind(GracidaRDA_Env, as.data.frame(RDA$CCA$biplot))
names(GracidaRDA_Env)[1] <- "Env_Variable"
Variables =c('Depth','Transparency','Shore proximity','pH','Oxygen','Temperature','Lake area','Substrate','Tilapia P/A')
GracidaRDA_Env$Env_Variable_Full <- Variables

RDA_EnvPlot <- ggscatter(GracidaRDA_Env, x = "RDA1", y = "RDA2",
                          ellipse = TRUE, ellipse.type = "FALSE", mean.point = FALSE,
                          star.plot = FALSE) +
  scale_x_continuous(limits = c(-0.8, .8), breaks =c(-0.8,0,0.8)) +
  scale_y_continuous(limits = c(-0.8, .8), breaks =c(-0.8,0,0.8)) +
  geom_label_repel(aes(label = Env_Variable_Full),
                   max.overlaps = Inf, point.padding = 0,
                   label.size = NA,segment.color = 'transparent') +
  geom_hline(yintercept=0, linetype="dashed", color = "blue") +
  geom_vline(xintercept=0, linetype="dashed", color = "blue") +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = RDA1,
                   yend = RDA2),colour='grey',
               arrow = arrow(length = unit(0.0, "cm"))) +
  xlab("RDA Axis 1 (18.1% of variation)") + ylab("RDA Axis 2 (5.7% of variation)")
RDA_EnvPlot

#Analysis 5: The combined figure

Combined_Figure_VII <- ggarrange(RDA_SitePlot,RDA_EnvPlot, common.legend = TRUE, 
                labels = c("A", "B"),
                ncol = 1, nrow = 2, legend="right")
Combined_Figure_VII 

#Analysis 6: Diversity vs. environmental variables

Diversity <- read.table("Diversity.txt",header=TRUE, fill=TRUE,sep="\t",check.names=FALSE)
Diversity$log10Depth <- log10(Diversity$Depth)

#Analysis 6: Associations between diversity and environmental variables using linear models

RichnessModel <-  lm(SpeciesRichness ~ Lake + log10Depth + Trans + Shore + pH + DO + Temp, data = Diversity)
summary(RichnessModel)
anova(RichnessModel)

colours = c("#999999", "#E69F00", "#56B4E9", "black", "red","#CCCC00")

SpeciesxDO <- ggplot(Diversity, aes(x = DO, y = SpeciesRichness)) +
  geom_point(aes(color = Lake)) + scale_colour_manual(values = colours) +
  geom_smooth(method = "lm",  se = TRUE) + theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Dissolved oxygen (mg/l)") + ylab("Species richness")
SpeciesxDO

SpeciesxDepth <- ggplot(Diversity, aes(x = log10(Depth), y = SpeciesRichness)) +
  geom_point(aes(color = Lake)) + scale_colour_manual(values = colours) +
  geom_smooth(method = "lm",  se = TRUE) + theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y ="Species richness", x =  bquote('Depth (m, '*log[10]*'transformed)'),)
SpeciesxDepth

#Analysis 6:Combined Figure

Combined_Figure_IV <- ggarrange(SpeciesxDO,SpeciesxDepth, common.legend = TRUE, 
                             labels = c("A", "B"),
                             ncol = 1, nrow = 2, legend="right")
Combined_Figure_IV

#Analysis 7:Barchart of mean catch per unit effort

Cichlids <- read.table("CichlidSummary.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
NonCichlids <- read.table("NonCichlidSummary.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#Analysis 7: Melt the dataframes and make plots

Cichlids_long <- melt(setDT(Cichlids), id.vars = c("Lake"), variable.name = "Taxon")

Cichlids_plot <- ggplot(data=Cichlids_long, aes(x=Lake, y=value, fill = Taxon)) +
  geom_bar(stat="identity") + theme_classic() + scale_fill_okabe_ito() +
  xlab("Species") + ylab("Mean number of individuals per sample") + guides(fill = guide_legend(reverse = FALSE)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
Cichlids_plot

NonCichlids_long <- melt(setDT(NonCichlids), id.vars = c("Lake"), variable.name = "Taxon")

NonCichlids_plot <- ggplot(data=NonCichlids_long, aes(x=Lake, y=value, fill = Taxon)) +
  geom_bar(stat="identity") + theme_classic() + scale_fill_okabe_ito() +
  xlab("Species") + ylab("Mean number of individuals per sample") + guides(fill = guide_legend(reverse = TRUE)) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
NonCichlids_plot

#Analysis 7: The combined figure of both plots

Combined_Figure_II <- ggarrange(Cichlids_plot,NonCichlids_plot, common.legend = FALSE, 
                              labels = c("A", "B"),
                              ncol = 1, nrow = 2, legend="right", align = c("hv"))
Combined_Figure_II

#End of analyses
