# This code covers all of the things necessary to do the parts in number 5, granted you 
# already have the sorted limma results. 

leflu <- read.csv("limma_results_leflu_sorted.csv")
ifo <- read.csv("limma_results_ifo_sorted.csv")
fluc <- read.csv("limma_results_fluc_sorted.csv")

n_rows = nrow(leflu)

N_leflu <- 0
N_ifo <- 0
N_fluc <- 0

# These steps allow you to get the differentially expressed genes. 

for (i in 1:n_rows) {
  if (leflu$adj.P.Val[i] < 0.05) {
    N_leflu <- N_leflu + 1
  }
  if (ifo$adj.P.Val[i] < 0.05) {
    N_ifo <- N_ifo + 1
  }
  if (fluc$adj.P.Val[i] < 0.05) {
    N_fluc <- N_fluc + 1
  }
}

# The top ten most expressed genes. The csv files will be used in the final report. 

fluc_tt <- fluc[1:10, c(2,7)]
colnames(fluc_tt) <- c("Probe Name", "Adj. P-value")
ifo_tt <- ifo[1:10, c(2,7)]
colnames(ifo_tt) <- c("Probe Name", "Adj. P-value")
leflu_tt <- leflu[1:10, c(2,7)]
colnames(leflu_tt) <- c("Probe Name", "Adj. P-value")

write.csv(fluc_tt, "fluc_tt.csv")
write.csv(ifo_tt, "ifo_tt.csv")
write.csv(leflu_tt, "leflu_tt.csv")

# These are the matrices that contain the data that we really care about. 
fluc_de <- fluc[1:N_fluc, ]
ifo_de <- ifo[1:N_ifo, ]
leflu_de <- leflu[1:N_leflu, ]

# Histogram of the Log fold changes in fluc. 
hist(fluc_de$logFC,
     main = "Log Fold Change of DE Genes in Fluconazole", 
     xlab = "Log Fold Change", 
     xlim = c(-4, 4),
     ylim = c(0, 700),
     breaks = 80)

# Histogram of the log fold changes in IFO. 
hist(ifo_de$logFC,
     main = "Log Fold Change of DE Genes in Ifosfamide", 
     xlab = "Log Fold Change", 
     xlim = c(-2, 2),
     breaks = 16)

# Histogram of the log fold changes in LEFLU. 
hist(leflu_de$logFC,
     main = "Log Fold Change of DE Genes in Leflunomide",
     xlab = "Log Fold Change",
     xlim = c(-4, 8),
     ylim= c(0,500),
     breaks = 48)

# P-val vs. log fold change. log plot so that way it 
# doesn't look the same at the previous thing. All log plots except 
# for IFo because ifo sucks and doesn't have enough stuff in it. 

plot( fluc_de$logFC, fluc_de$adj.P.Val, type = "p", pch = 20,
     log = "y",
     xlab = "Log Fold Change", ylab = "Nominal p-value", 
     main = "Fluconazole: Log Fold Change vs. Nominal p-value")

plot(ifo_de$logFC, ifo_de$adj.P.Val, type = "p", pch = 20,
     xlab = "Log Fold Change", ylab = "Nominal p-value", 
     xlim = c(-2, 2),
     ylim = c(0,0.05),
     main = "Ifosfamide: Log Fold Change vs. Nominal p-value")

plot(leflu_de$logFC, leflu_de$adj.P.Val, type = "p", pch = 20,
     log = "y", 
     xlab = "Log Fold Change", ylab = "Nominal p-value", 
     xlim = c(-7.5, 7.5),
     ylim = c(0,0.05),
     main = "Leflunomide: Log Fold Change vs. Nominal p-value")
