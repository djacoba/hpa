library(tidyverse)

#prepare HPA IHC data
hpa.kidney <- read_csv("hpa_kidney.csv")
hpa.kidney %>% dplyr::filter(Reliability != "Uncertain") -> hpa.kidney.filt

#separate glomerular and tubular data
hpa.kidney.filt %>% dplyr::filter(`Cell type` == "cells in glomeruli") -> hpa.glom.filt
hpa.kidney.filt %>% dplyr::filter(`Cell type` == "cells in tubules") -> hpa.tub.filt

#prepare HPA RNA-seq data
rna.kidney <- read_csv("rna_kidney.csv")
rna.kidney %>% dplyr::mutate(Abundance = case_when(nTPM >= 13.3 ~ "High",
                                                   nTPM >= 3.0 & nTPM < 13.3 ~ "Medium",
                                                   nTPM < 3.0 & nTPM > 0 ~ "Low",
                                                   nTPM == 0 ~ "Not detected")) -> rna.kidney

rna.kidney %>% dplyr::select(Gene, nTPM, Abundance) -> rna.kidney

#merge IHC and RNA-seq data
left_join(hpa.glom.filt, rna.kidney) -> hpa.glom
left_join(hpa.tub.filt, rna.kidney) -> hpa.tub

#make sure expression levels are set as factors
hpa.glom$Level <- factor(hpa.glom$Level, levels=c('Not detected', 'Low', 'Medium', 'High'))
hpa.glom$Abundance <- factor(hpa.glom$Abundance, levels=c('Not detected', 'Low', 'Medium', 'High'))

hpa.tub$Level <- factor(hpa.tub$Level, levels=c('Not detected', 'Low', 'Medium', 'High'))
hpa.tub$Abundance <- factor(hpa.tub$Abundance, levels=c('Not detected', 'Low', 'Medium', 'High'))

#chi-square test
chi.glom <- chisq.test(hpa.glom$Level, hpa.glom$Abundance)
chi.glom
chi.glom$observed

chi.tub <- chisq.test(hpa.tub$Level, hpa.tub$Abundance)
chi.tub
chi.tub$observed

#filter for mRNA-protein pairs with discordant expression levels
hpa.glom %>% dplyr::filter(Level != Abundance) -> glom.uneq
hpa.tub %>% dplyr::filter(Level != Abundance) -> tub.uneq

#figure for glom (Fig 1A)
col.cor <- ifelse(hpa.glom$Level==hpa.glom$Abundance, "#5785C1", 
                  ifelse(hpa.glom$Level == "Not detected" & hpa.glom$Abundance != "Not detected", "#C93312",
                         ifelse(hpa.glom$Level != "Not detected" & hpa.glom$Abundance == "Not detected", "#C93312","#899DA4")))
ggplot(hpa.glom, aes(Level,Abundance)) + geom_jitter(position = position_jitter(height = .5, width=.5), shape=21, 
                                                     size=5, alpha=0.7, colour="black", fill=col.cor) +
  theme_bw() + theme(axis.text=element_text(size=30), axis.title=element_text(size=30), plot.title=element_text(size=30)) +
  labs(x = "Glomerular Protein Level", y = "Cortical mRNA Level")

dev.copy(png,"new_glom2.png", width=800, height=800)
dev.off()

#figure for tub (Fig 1B)
col.cor <- ifelse(hpa.tub$Level==hpa.tub$Abundance, "#5785C1", 
                  ifelse(hpa.tub$Level == "Not detected" & hpa.tub$Abundance != "Not detected", "#C93312",
                         ifelse(hpa.tub$Level != "Not detected" & hpa.tub$Abundance == "Not detected", "#C93312","#899DA4")))
ggplot(hpa.tub, aes(Level,Abundance)) + geom_jitter(position = position_jitter(height = .5, width=.5), shape=21, 
                                                    size=5, alpha=0.7, colour="black", fill=col.cor) +
  theme_bw() + theme(axis.text=element_text(size=30), axis.title=element_text(size=30), plot.title=element_text(size=30)) +
  labs(x = "Tubular Protein Level", y = "Cortical mRNA Level")

dev.copy(png,"new_tub2.png", width=800, height=800)
dev.off()

#Table 1

#counting how many mRNA-protein pairs are concordant and not
sum(hpa.glom$Level == hpa.glom$Abundance)
sum(hpa.glom$Level != hpa.glom$Abundance)
sum(hpa.tub$Level == hpa.tub$Abundance)
sum(hpa.tub$Level != hpa.tub$Abundance)

#discordant pairs
glom.uneq %>% dplyr::filter(Level=="Not detected" & Abundance != "Not detected") %>% write_csv("glom_uneq_noprot_new.csv") #4953
tub.uneq %>% dplyr::filter(Level=="Not detected" & Abundance != "Not detected") %>% write_csv("tub_uneq_noprot_new.csv") #2122
glom.uneq %>% dplyr::filter(Level!="Not detected" & Abundance == "Not detected") %>% write_csv("glom_uneq_norna_new.csv") #76
tub.uneq %>% dplyr::filter(Level!="Not detected" & Abundance == "Not detected") %>% write_csv("tub_uneq_norna_new.csv") #141

#read RNA-seq validation dataset
levin.glom <- read_csv("levin_glom.csv")
levin.tub <- read_csv("levin_tub.csv")

#pairs with no mRNA but with IHC signals (detectable in Levin RNA-seq) SuppTab1
glom.norna <- read_csv("glom_uneq_norna_new.csv")
inner_join(glom.norna, levin.glom)  %>% dplyr::filter(mean>0) #62
inner_join(glom.norna, levin.glom)  %>% dplyr::filter(mean>=1) #5
tub.norna <- read_csv("tub_uneq_norna_new.csv")
inner_join(tub.norna, levin.tub)  %>% dplyr::filter(mean>0) #120
inner_join(tub.norna, levin.tub)  %>% dplyr::filter(mean>=1) #3


#pairs with no IHC signals but with RNA (nTPM >=1)
glom.noprot <- read_csv("glom_uneq_noprot_new.csv")
glom.noprot %>% dplyr::filter(nTPM >=1) #3773
tub.noprot <- read_csv("tub_uneq_noprot_new.csv")
tub.noprot %>% dplyr::filter(nTPM >=1) #1172

#read mass spec validation dataset
kpmp.glom <- read_csv("kpmp_glom2.csv")
kpmp.tub <- read_csv("kpmp_tub2.csv")

#pairs with no IHC signals but with RNA (detectable in mass spec) SuppTab1
inner_join(glom.noprot, kpmp.glom, by = c("Gene" = "Ensembl gene ID")) %>% dplyr::filter(mean >= 1) %>% View() #637
inner_join(tub.noprot, kpmp.tub, by = c("Gene" = "Ensembl gene ID")) %>% dplyr::filter(mean >= 1) %>% View() #171
