library(limma)
library(data.table)
library(dplyr)
library(matrixStats)

#The toxgroup samples
samples <- read.csv('/project/bf528/project_3/groups/group_3_mic_info.csv')

# The normalization table.
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)

# Subsets of the normalization dat that correspond to Leflunomide, fluconazole, and ifosfamide, respectively. 

rma.subset.leflu <-rma[paste0('X', samples$array_id[samples$chemical=='LEFLUNOMIDE' | samples$chemical == 'Control'])]

rma.subset.fluc <-rma[paste0('X', samples$array_id[samples$chemical=='FLUCONAZOLE' | samples$chemical == 'Control'])]

rma.subset.ifo <-rma[paste0('X', samples$array_id[samples$chemical=='IFOSFAMIDE' | samples$chemical == 'Control'])]


# The design matrices for our respective chemicals

design.leflu <- model.matrix(~factor(samples$chemical, levels=c("Control", "LEFLUNOMIDE"))
colnames(design.leflu) <- c("Intercept", "LEFLUNOMIDE")

design.fluc <- model.matrix(~factor(samples$chemical, levels=c("Control", "FLUCONAZOLE"))
colnames(design.fluc) <- c("Intercept", "FLUCONAZOLE")

design.ifo <- model.matrix(~factor(samples$chemical, levels=c("Control", "IFOSFAMIDE"))
colnames(design.ifo) <- c("Intercept", "IFOSFAMIDE")

# Leflunomide Results

fit.leflu <- lmFit(rma.subset.leflu, design.leflu)
fit.leflu <- eBayes(fit.leflu)

t.leflu <- topTable(fit.leflu, coef=2, n=nrow(rma.subset.leflu), adjust='BH')

write.csv(t.leflu,'limma_results_leflu.csv')

# Floconazole Results 

fit.fluc <- lmFit(rma.subset.fluc, design.fluc)
fit.fluc <- eBayes(fit.fluc)

t.fluc <- topTable(fit.fluc, coef=2, n=nrow(rma.subset.fluc), adjust='BH')

write.csv(t.fluc,'limma_results_fluc.csv')

#Ifosfamide Results

fit.ifo <- lmFit(rma.subset.ifo, design.ifo)
fit.ifo <- eBayes(fit.ifo)

t.ifo <- topTable(fit.ifo, coef=2, n=nrow(rma.subset.ifo), adjust='BH')

write.csv(t.ifo,'limma_results_ifo.csv')
