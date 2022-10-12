#Ensure the following packgages are installed and live

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(pairwiseAdonis)
library(ggpubr)
library(data.table)
library(ggokabeito)

#Analysis 1: A non-metric multidimensional scaling analysis, ordination plot and associated Permanovas.

#Analysis 1: Read in species abundance dataframe, which is log10(x+1) transformed

Gracida_Species93 <- read.table("Gracida_Species93.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#Analysis 1: Conduct a non-metric multidimensional scaling analysis (nMDS)

Gracida_Species93_MDS <- metaMDS(Gracida_Species93[,2:19], labels=Gracida_Species93$Site, trymax = 100)

#Analysis 1: Construct a dataframe of the nMDS values of each sample

NMDS = data.frame(MDS1 = Gracida_Species93_MDS$points[,1], MDS2=Gracida_Species93_MDS$points[,2],group=Gracida_Species93$Site)

#Analysis 1: Plot the nMDS, here using a star plot, save as 6 x 8

ggscatter(NMDS, x = "MDS1", y = "MDS2", color = "group",
                       ellipse = FALSE, ellipse.type = "convex", mean.point = TRUE,
                       star.plot = TRUE, legend = "right") +
  scale_x_continuous(limits = c(-1.5, 1.5), breaks =c(-1.5,-1,-0.5,0,0.5,1,1.5)) +
  scale_y_continuous(limits = c(-1.5, 1.5), breaks =c(-1.5,-1,-0.5,0,0.5,1,1.5)) +
  xlab("nMDS Axis 1") + ylab("nMDS Axis 2") + scale_color_okabe_ito()

#Analysis 1: Conduct a Permanova (permutation Anova)

Gracida_Species93$Site <- as.factor(Gracida_Species93$Site)

adonis(Gracida_Species93[,2:19] ~ Site, data=Gracida_Species93, method ="bray", permutations=10000)

#Analysis 1: Conduct a Permanova (permutation Anova)

pairwise.adonis(Gracida_Species93[,2:19], Gracida_Species93$Site)

#Analysis 2: Canonical Correspondence Analysis, associated Anova and various plots on all species except Oreochromis

#Analysis 2: Read in environmental and species abundance dataframe

Gracida_EnvPredic93 <- read.table("Gracida_EnvPredic93.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
Gracida_Native93 <- cbind(Gracida_Species93[,1:10],Gracida_Species93[,12:19])

#Analysis 2: Canonical Correspondence Analysis across all sites

CCA_All<- cca (Gracida_Native93 [,2:18] ~ Area + Depth + Trans + Shore + pH + DO + Temp + Oreo + Subs, data = Gracida_EnvPredic93)
summary(CCA_All)

#Analysis 2: Test of the statistical signficance of environmental predictors in the Canonical Correspondence Analysis

anova(CCA_All, by ="margin", permutations=10000)

#Analysis 2: Place CCA values in dataframes for plotting

GracidaCCA_Sites <- CCA_All$CCA$wa
GracidaCCA_Sites <- cbind(GracidaCCA_Sites, Gracida_Native93[1:1])
GracidaCCA_Sites

GracidaCCA_Env <- read.table("VariableNames.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
GracidaCCA_Env <- cbind(GracidaCCA_Env,CCA_All$CCA$biplot)
GracidaCCA_Env

GracidaCCA_Species <- read.table("SpeciesNames.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
GracidaCCA_Species <- cbind(GracidaCCA_Species,CCA_All$CCA$v)
GracidaCCA_Species

#Analysis 2: Place CCA values in dataframes for plotting

#Analysis 2: The plot of species values

CCA_SpeciesPlot <- ggscatter(GracidaCCA_Species, x = "CCA1", y = "CCA2",
                             ellipse = FALSE, ellipse.type = "convex", mean.point = FALSE,
                             star.plot = FALSE) +
  scale_x_continuous(limits = c(-4, 4), breaks =c(-4,-2,0,2,4)) +
  scale_y_continuous(limits = c(-4, 4), breaks =c(-4,-2,0,2,4)) +
  geom_label_repel(aes(label = Species),
                   max.overlaps = Inf,
                   segment.color = 'grey',point.padding = 0.3) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue") +
  geom_vline(xintercept=0, linetype="dashed", color = "blue") +
  xlab("CCA Axis 1 (6.7% of variation)") + ylab("CCA Axis 2 (5.6% of variation)") +
  scale_color_okabe_ito()
CCA_SpeciesPlot

#Analysis 2: The plot of site values

CCA_SitePlot <- ggscatter(GracidaCCA_Sites, x = "CCA1", y = "CCA2", color = "Site",
          ellipse = FALSE, ellipse.type = "convex", mean.point = TRUE,
          star.plot = TRUE) +
  geom_hline(yintercept=0, linetype="dashed", color = "blue") +
  geom_vline(xintercept=0, linetype="dashed", color = "blue") +
  scale_x_continuous(limits = c(-5, 5), breaks =c(-5,-2.5,0,2.5,5)) +
  scale_y_continuous(limits = c(-5, 5), breaks =c(-5,-2.5,0,2.5,5)) +
  xlab("CCA Axis 1 (6.7% of variation)") + ylab("CCA Axis 2 (5.6% of variation)") +
  scale_color_okabe_ito()
CCA_SitePlot

#Analysis 2: The plot of environmental values

CCA_EnvPlot <- ggscatter(GracidaCCA_Env, x = "CCA1", y = "CCA2",
          ellipse = FALSE, ellipse.type = "convex", mean.point = FALSE,
          star.plot = FALSE) +
  scale_x_continuous(limits = c(-1, 1), breaks =c(-1,-0.5,0,0.5,1)) +
  scale_y_continuous(limits = c(-1, 1), breaks =c(-1,-0.5,0,0.5,1)) +
  geom_label_repel(aes(label = Variable),
                   max.overlaps = Inf, point.padding = 0,
                   label.size = NA,segment.color = 'transparent') +
  geom_hline(yintercept=0, linetype="dashed", color = "blue") +
  geom_vline(xintercept=0, linetype="dashed", color = "blue") +
  geom_segment(aes(x = 0,
                 y = 0,
                 xend = CCA1,
                 yend = CCA2),colour='grey',
             arrow = arrow(length = unit(0.0, "cm"))) +
  xlab("CCA Axis 1 (6.7% of variation)") + ylab("CCA Axis 2 (5.6% of variation)") +
  scale_color_okabe_ito()
CCA_EnvPlot

#Analysis 2: The combined figure of all three (export as 10 height x 12 width, note this will need to be edited to move the legend)

Combined_Figure1 <- ggarrange(CCA_SpeciesPlot,CCA_SitePlot,CCA_EnvPlot, common.legend = TRUE, 
                labels = c("A", "B", "C"),
                ncol = 2, nrow = 2, legend="right")
Combined_Figure1 

#Analysis 3: Diversity vs. environmental variables

#Analysis 3: Read in diversity and abundance dataframe

Diversity <- read.table("Diversity.txt",header=TRUE, fill=TRUE,sep="\t",check.names=FALSE)

#Analysis 3: Associations between abundance, diversity and environmental variables using linear models

AbundanceModel <-  lm(Abundance ~ Site + Depth + Secchi + Shore + pH + DO + Temperature, data = Diversity)
summary(AbundanceModel)
anova(AbundanceModel)

RichnessModel <-  lm(Species ~ Site + Depth + Secchi + Shore + pH + DO + Temperature, data = Diversity)
summary(RichnessModel)
anova(RichnessModel)

#Analysis 3:The combined figure of all three (export as 6 height x 8 width, note this will need to be edited to move the legend)

SpeciesxDO <- ggplot(Diversity, aes(x = DO, y = Species)) +
  geom_point(aes(color = Site)) +
  geom_smooth(aes(color = Site, fill = Site), method = "lm",  se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Dissolved oxygen (mg/l))") + ylab("Species richness") +
  scale_color_okabe_ito()
SpeciesxDO

SpeciesxpH <- ggplot(Diversity, aes(x = pH, y = Species)) +
  geom_point(aes(color = Site)) +
  geom_smooth(aes(color = Site, fill = Site), method = "lm",  se = FALSE) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("pH") + ylab("Species richness") +
  scale_color_okabe_ito()
SpeciesxpH

DOxpH <- ggplot(Diversity, aes(x = DO, y = pH)) +
  geom_point(aes(color = Site)) +
  geom_smooth(aes(color = Site, fill = Site), method = "lm",  se = FALSE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Dissolved oxygen (mg/l))") + ylab("pH") +
  scale_color_okabe_ito()
DOxpH

#Analysis 3:Combined Figure, save as 8 x 10

Combined_Figure2 <- ggarrange(SpeciesxDO,SpeciesxpH,DOxpH, common.legend = TRUE, 
                             labels = c("A", "B", "C"),
                             ncol = 2, nrow = 2, legend="right")
Combined_Figure2 

#Analysis 4:Barchart of mean catch per unit effort

#Analysis 4: Read in mean catch per sample per location dataframe

Cichlids <- read.table("CichlidSummary.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#Analysis 4: Melt the dataframes and make plots

Cichlids_long <- melt(setDT(Cichlids), id.vars = c("Lake"), variable.name = "Taxon")
Cichlids_plot <- ggplot(data=Cichlids_long, aes(x=Lake, y=value, fill = Taxon)) +
  geom_bar(stat="identity") + theme_classic() + scale_fill_okabe_ito() +
  xlab("Species") + ylab("Mean number of individuals per sample")
Cichlids_plot

NonCichlids <- read.table("NonCichlidSummary.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)
NonCichlids_long <- melt(setDT(NonCichlids), id.vars = c("Lake"), variable.name = "Taxon")
NonCichlids_plot <- ggplot(data=NonCichlids_long, aes(x=Lake, y=value, fill = Taxon)) +
  geom_bar(stat="identity") + theme_classic() + scale_fill_okabe_ito() +
  xlab("Species") + ylab("Mean number of individuals per sample")
NonCichlids_plot

#Analysis 4: The combined figure of both plots (export as 10 height x 10 width)

Combined_Figure3 <- ggarrange(Cichlids_plot,NonCichlids_plot, common.legend = FALSE, 
                              labels = c("A", "B"),
                              ncol = 1, nrow = 2, legend="right")
Combined_Figure3 

#Analysis 5: Rarefaction of species diversity

#Analysis 5: Read in the species abundance dataframe

Gracida_Species93_t <- read.table("Gracida_Species93_transpose.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

#Analysis 5: Make the six accumulation curves, using raw capture data, and print each one to a pdf

pdf(file = 'Noh_Bec.pdf')
sp_raw_NohBec <- specaccum(Gracida_Species93_t[2:17], "random")
NohBec_plot <- plot(sp_raw, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw, col="brown3", add=TRUE)
dev.off()

pdf(file = 'Petcacab.pdf')
sp_raw_Petcacab <- specaccum(Gracida_Species93_t[18:33], "random")
Petcacab_plot <- plot(sp_raw, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw, col="brown3", add=TRUE)
dev.off()

pdf(file = 'Santa_Teresita.pdf')
sp_raw_Santa_Teresita <- specaccum(Gracida_Species93_t[34:49], "random")
Santa_Teresita_plot <- plot(sp_raw, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw, col="brown3", add=TRUE)
dev.off()

pdf(file = 'San_Felipe.pdf')
sp_raw_San_Felipe <- specaccum(Gracida_Species93_t[50:65], "random")
San_Felipe_plot <- plot(sp_raw, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw, col="brown3", add=TRUE)
dev.off()

pdf(file = 'Coabas.pdf')
sp_raw_Coabas <- specaccum(Gracida_Species93_t[66:80], "random")
Coabas_plot <- plot(sp_raw, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw, col="brown3", add=TRUE)
dev.off()

pdf(file = 'Chacanbacan.pdf')
sp_raw_Chacanbacan <- specaccum(Gracida_Species93_t[81:94], "random")
Chacanbacan_plot <- plot(sp_raw, ci.type="poly", col="blue", ci.lty=0, ci.col="lightblue", las=1, xlim=c(1, 16), ylim=c(0,20))
boxplot(sp_raw, col="brown3", add=TRUE)
dev.off()

#End of analyses
