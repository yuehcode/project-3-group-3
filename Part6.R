# This code covers everything that needed to be done
# in part 6 of Project 3 for BF528, Spring 2020. 

# The RNA-seq data. 
ifo_seq <- read.csv('deseq_ifo.csv', as.is = TRUE)
fluc_seq <- read.csv('deseq_fluc.csv', as.is = TRUE)
leflu_seq <- read.csv('deseq_leflu.csv', as.is = TRUE)

# THis is the matrix that we will use to map the IDs for the 
# RNa to the probes in the microarray data later. 
names <- read.csv('refseq_affy_map.csv', as.is = TRUE)

n_names <- nrow(names)

# Getting the matrices to have only the markers that were 
# at adj p-val of 0.05. 

# Sorting by p-adj is fine since we've only got a couple thousand data points.
# If this were much larger then I'd have to do this
# another way.

ifo_seq <- ifo_seq[order(ifo_seq$padj), ]
fluc_seq <- fluc_seq[order(fluc_seq$padj), ]
leflu_seq <- leflu_seq[order(leflu_seq$padj), ]

n_ifo <- 0
n_fluc <- 0
n_leflu <- 0

# Figure out the count of DE genes in each matrix. 
# This is a really dumb way of doing it,
# but it works so I don't care. 

for (i in 1:length(ifo_seq$X)) {
  if (ifo_seq$padj[i] < 0.05) {
    n_ifo <- n_ifo + 1
  }
}

for (i in 1:length(fluc_seq$X)) {
  if (fluc_seq$padj[i] < 0.05) {
    n_fluc <- n_fluc + 1
  }
}

for (i in 1:length(leflu_seq$X)) {
  if (leflu_seq$padj[i] < 0.05) {
    n_leflu <- n_leflu + 1
  }
}

# Now these matrices are only the DE genes, sorted
# by adjusted p-value. 

ifo_seq <- ifo_seq[1:n_ifo, ]
fluc_seq <- fluc_seq[1:n_fluc, ]
leflu_seq <- leflu_seq[1:n_leflu, ]

# End of creating relevant matrices

# Mapping refseqIDs to probeset IDs using the names matrix.
# Some refseqIDs don't have a corresponding probeset ID,
# so for those we just label them "NOPROBE" so that 
# way it's easier to catch later on and I don't 
# have to deal with any weird stuff with nulls. 

# I'm sure there's a way I could do all this without 
# repeating the code three times, but it wasn't really worth 
# it for me here. 

for (i in 1:n_ifo) {
  seq_name <- ifo_seq$X[i]
  if (seq_name %in% names$REFSEQ) {
    probe <- names$PROBEID[which(names$REFSEQ == seq_name)][1]
    if (is.na(probe)) {
      ifo_seq$X[i] <- "NOPROBE"
    }
    else {
      ifo_seq$X[i] <- names$PROBEID[which(names$REFSEQ == seq_name)]
    }
  }
  else{
    ifo_seq$X[i] <- "NONAME"
  }
}

for (i in 1:n_fluc) {
  seq_name <- fluc_seq$X[i]
  if (seq_name %in% names$REFSEQ) {
    probe <- names$PROBEID[which(names$REFSEQ == seq_name)][1]
    if (is.na(probe)) {
      fluc_seq$X[i] <- "NOPROBE"
    }
    else {
      fluc_seq$X[i] <- names$PROBEID[which(names$REFSEQ == seq_name)]
    }
  }
  else{
    fluc_seq$X[i] <- "NONAME"
  }
}

for (i in 1:n_leflu) {
  seq_name <- leflu_seq$X[i]
  if (seq_name %in% names$REFSEQ) {
    probe <- names$PROBEID[which(names$REFSEQ == seq_name)][1]
    if (is.na(probe)) {
      leflu_seq$X[i] <- "NOPROBE"
    }
    else {
      leflu_seq$X[i] <- names$PROBEID[which(names$REFSEQ == seq_name)]
    }
  }
  else{
    leflu_seq$X[i] <- "NONAME"
  }
}

# End of matching probes to stuff.

# The microarray data from earlier, so that way we can calculate concordance. 

leflu_mic <- read.csv("limma_results_leflu_sorted.csv")
ifo_mic <- read.csv("limma_results_ifo_sorted.csv")
fluc_mic <- read.csv("limma_results_fluc_sorted.csv")

leflu_mic <- leflu_mic[order(leflu_mic$adj.P.Val), ]
leflu_mic <- leflu_mic[1:which(leflu_mic$adj.P.Val > 0.05)[1], ]

ifo_mic <- ifo_mic[order(ifo_mic$adj.P.Val), ]
ifo_mic <- ifo_mic[1:which(ifo_mic$adj.P.Val > 0.05)[1], ]

fluc_mic <- fluc_mic[order(fluc_mic$adj.P.Val), ]
fluc_mic <- fluc_mic[1:which(fluc_mic$adj.P.Val > 0.05)[1], ]

# All the microarray data has been taken to only have DE genes at p_adj < 0.05
# and they are all ordered on adj.P.Val. 

#Ifo concordance

# If a thing's in both sets and its logFC is in the same direction
# for both, gets considered an intersection. 
# I tried doing the 
# background corrected, but I got really small negative numbers for IFO, and
# none of the other values were between 0 and 1. 

n_0 <- 0 # Number of observed items in the intersection
n_1 <- length(ifo_mic$X) # Number of things in subset 1
n_2 <- length(ifo_seq$X) # Number of things in subset 2
N <- n_1 + n_2 # Total number of things in the whole set. 

for (i in 1:n_1) {
  if (ifo_mic$X[i] %in% ifo_seq$X) {
    
    idx <- which(ifo_seq$X == ifo_mic$X[i])[1]
    if ((ifo_mic$logFC[i] > 0) & (ifo_seq$log2FoldChange[idx] > 0)) {
      n_0 <- n_0 + 1
    }
    if (ifo_mic$logFC[i] < 0 & ifo_seq$log2FoldChange[idx] < 0) {
      n_0 <- n_0 + 1
    }
  }
}

conc_ifo <- (2 * (n_0 * N - n_1 * n_2)/(n_0 + N -n_1 - n_2)) / (n_1 + n_2)
conc_ifo_unadj <- 2*n_0 / N

#Fluc concordance

n_0 <- 0 # Number of observed items in the intersection
n_1 <- length(ifo_seq$X) # Number of things in subset 1
n_2 <- length(ifo_mic$X) # Number of things in subset 2

N <- n_1 + n_2 # Total number of things in the whole set. 

for (i in 1:n_1) {
  if (fluc_mic$X[i] %in% fluc_seq$X) {
    
    idx <- which(fluc_seq$X == fluc_mic$X[i])[1]
    if ((fluc_mic$logFC[i] > 0) & (fluc_seq$log2FoldChange[idx] > 0)) {
      n_0 <- n_0 + 1
    }
    if (fluc_mic$logFC[i] < 0 & fluc_seq$log2FoldChange[idx] < 0) {
      n_0 <- n_0 + 1
    }
  }
}

conc_fluc <- (2 * (n_0 * N - n_1 * n_2)/(n_0 + N -n_1 - n_2)) / (n_1 + n_2)
conc_fluc_unadj <- 2*n_0 / N


#Leflu Concordance 

n_0 <- 0 # Number of observed items in the intersection
n_1 <- length(leflu_seq$X) # Number of things in subset 1
n_2 <- length(leflu_mic$X) # Number of things in subset 2

N <- n_1 + n_2 # Total number of things in the whole set. 

for (i in 1:n_1) {
  if (leflu_mic$X[i] %in% leflu_seq$X) {
    
    idx <- which(leflu_seq$X == leflu_mic$X[i])[1]
    if ((leflu_mic$logFC[i] > 0) & (leflu_seq$log2FoldChange[idx] > 0)) {
      n_0 <- n_0 + 1
    }
    if (leflu_mic$logFC[i] < 0 & leflu_seq$log2FoldChange[idx] < 0) {
      n_0 <- n_0 + 1
    }
  }
}

conc_leflu <- (2 * (n_0 * N - n_1 * n_2)/(n_0 + N -n_1 - n_2)) / (n_1 + n_2)
conc_leflu_unadj <- 2*n_0 / N

# End of calculating concordance

# Begin plotting conc. vs. # of DEGs (RNA-seq) and plotting conc. vs. DEGs (microarray)

plt_pnts <- matrix(c("IFO", length(ifo_seq$X), conc_ifo_unadj, "FLU", length(fluc_seq$X), conc_fluc_unadj, 
                     "LEF", length(leflu_seq$X), conc_leflu_unadj)
                   , nrow = 3, ncol = 3, byrow=TRUE)

# I just added the labels later in Photoshop. 

plot(plt_pnts[,2:3], xlab = "Treatment Effect \n (number of DEGs from RNA-seq)",
     ylab = "Concordance of DEG")


plt_pnts2 <- matrix(c("IFO", length(ifo_mic$X), conc_ifo_unadj, "FLU", length(fluc_mic$X), conc_fluc_unadj, 
                     "LEF", length(leflu_mic$X), conc_leflu_unadj)
                   , nrow = 3, ncol = 3, byrow=TRUE)

plot(plt_pnts2[,2:3], xlab = "Treatment Effect \n (number of DEGs from Microarray)",
     ylab = "Concordance of DEG")

# End plotting DEGs vs. stuff

# Begin plotting of above and below median genes of concordance. 

# Separate into above and below median groups. 

fluc_bm <- fluc_seq[1:ceiling(length(fluc_seq$X) / 2), ]
fluc_am <- fluc_seq[(ceiling(length(fluc_seq$X) / 2) + 1):length(fluc_seq$X), ]

ifo_bm <- ifo_seq[1:ceiling(length(ifo_seq$X) / 2), ]
ifo_am <- ifo_seq[(ceiling(length(ifo_seq$X) / 2) + 1):length(ifo_seq$X), ]

leflu_bm <- leflu_seq[1:ceiling(length(leflu_seq$X) / 2), ]
leflu_am <- leflu_seq[(ceiling(length(leflu_seq$X) / 2) + 1):length(leflu_seq$X), ]

# Calculate the concordance for every point in below-median FLUC points. 

# There is probably an easier way to do this, but I'm not sure how to do 
# something like that in R, and I was running out of time. 

#FLU BM calculations. 

n_0 <- 0
N <- 0

conc_fluc_bm <- list(0)

for (i in 1:length(fluc_bm$X)) {
  N <- N + 1
  if (fluc_bm$X[i] %in% fluc_mic$X) {
    idx <- which (fluc_mic$X == fluc_bm$X[i])[1]
    
    if ((fluc_bm$log2FoldChange[i] > 0 & fluc_mic$logFC[idx] > 0) | 
        (fluc_bm$log2FoldChange[i] < 0 & fluc_mic$logFC[idx] < 0)) {
      n_0 <- n_0 + 1
      N <- N + 1
    }
  }
  
  
  if (N == 0) {
    new_conc <- 0
  }
  else {
    new_conc <- 2 * n_0 / N * 100
  }
  conc_fluc_bm <- c(conc_fluc_bm, new_conc)
}

# Same thing but for IFO. 

n_0 <- 0
N <- 0

conc_ifo_bm <- list(0)

for (i in 1:length(ifo_bm$X)) {
  N <- N + 1
  if (ifo_bm$X[i] %in% ifo_mic$X) {
    idx <- which (ifo_mic$X == ifo_bm$X[i])[1]
    
    if ((ifo_bm$log2FoldChange[i] > 0 & ifo_mic$logFC[idx] > 0) | 
        (ifo_bm$log2FoldChange[i] < 0 & ifo_mic$logFC[idx] < 0)) {
      n_0 <- n_0 + 1
      N <- N + 1
    }
  }
  
  if (N == 0) {
    new_conc <- 0
  }
  else {
    new_conc <- 2 * n_0 / N * 100
  }
  conc_ifo_bm <- c(conc_ifo_bm, new_conc)
}

#Aaaaand same thing for LEFLU

n_0 <- 0
N <- 0

conc_leflu_bm <- list(0)

for (i in 1:length(leflu_bm$X)) {
  N <- N + 1
  if (leflu_bm$X[i] %in% leflu_mic$X) {
    idx <- which (leflu_mic$X == leflu_bm$X[i])[1]
    
    if ((leflu_bm$log2FoldChange[i] > 0 & leflu_mic$logFC[idx] > 0) | 
        (leflu_bm$log2FoldChange[i] < 0 & leflu_mic$logFC[idx] < 0)) {
      n_0 <- n_0 + 1
      N <- N + 1
    }
  }
  
  if (N == 0) {
    new_conc <- 0
  }
  else {
    new_conc <- 2 * n_0 / N * 100
  }
  conc_leflu_bm <- c(conc_leflu_bm, new_conc)
}

# Plot the results
plot(1:length(conc_fluc_bm), conc_fluc_bm, type = "l", col = "red",
     main = "Below-median Expression",
     xlab = "Number of DEGS Ordered by Fold Change", 
     ylab = "Cross Platform DEG Concordance (%)",
     #xlim = c(0, 600),
     ylim = c(0,100))

lines(1:length(conc_ifo_bm), conc_ifo_bm, col = "blue")
lines(1:length(conc_leflu_bm), conc_leflu_bm, col = "green")

legend(100, 100, legend=c("FLU", "IFO", "LEF"), col = c("red", 'blue', 'green'), lty = 1)

#Do the same thing for Above-Median Expression

#FLU stuff
# It's the exact same code, but with am instead of bm. 

n_0 <- 0
N <- 0

conc_fluc_am <- list(0)

for (i in 1:length(fluc_am$X)) {
  N <- N + 1
  if (fluc_am$X[i] %in% fluc_mic$X) {
    idx <- which (fluc_mic$X == fluc_am$X[i])[1]
    
    if ((fluc_am$log2FoldChange[i] > 0 & fluc_mic$logFC[idx] > 0) | 
        (fluc_am$log2FoldChange[i] < 0 & fluc_mic$logFC[idx] < 0)) {
      n_0 <- n_0 + 1
      N <- N + 1
    }
  }
  
  if (N == 0) {
    new_conc <- 0
  }
  else {
    new_conc <- 2 * n_0 / N * 100
  }
  conc_fluc_am <- c(conc_fluc_am, new_conc)
}

# Same thing but for IFO. 

n_0 <- 0
N <- 0

conc_ifo_am <- list(0)

for (i in 1:length(ifo_am$X)) {
  N <- N + 1
  if (ifo_am$X[i] %in% ifo_mic$X) {
    idx <- which (ifo_mic$X == ifo_am$X[i])[1]
    
    if ((ifo_am$log2FoldChange[i] > 0 & ifo_mic$logFC[idx] > 0) | 
        (ifo_am$log2FoldChange[i] < 0 & ifo_mic$logFC[idx] < 0)) {
      n_0 <- n_0 + 1
      N <- N + 1
    }
  }
  
  if (N == 0) {
    new_conc <- 0
  }
  else {
    new_conc <- 2 * n_0 / N * 100
  }
  conc_ifo_am <- c(conc_ifo_am, new_conc)
}

#Aaaaand same thing for LEFLU

n_0 <- 0
N <- 0

conc_leflu_am <- list(0)

for (i in 1:length(leflu_am$X)) {
  N <- N + 1
  if (leflu_am$X[i] %in% leflu_mic$X) {
    idx <- which (leflu_mic$X == leflu_am$X[i])[1]
    
    if ((leflu_am$log2FoldChange[i] > 0 & leflu_mic$logFC[idx] > 0) | 
        (leflu_am$log2FoldChange[i] < 0 & leflu_mic$logFC[idx] < 0)) {
      n_0 <- n_0 + 1
      N <- N + 1 
    }
  }
  
  if (N == 0) {
    new_conc <- 0
  }
  else {
    new_conc <- 2 * n_0 / N * 100
  }
  conc_leflu_am <- c(conc_leflu_am, new_conc)
}

# Plot the results
plot(1:length(conc_fluc_am), conc_fluc_am, type = "l", col = "red",
     main = "Above-median Expression",
     xlab = "Number of DEGS Ordered by Fold Change", 
     ylab = "Cross Platform DEG Concordance (%)",
     #xlim = c(0,900),
     ylim = c(0,100))

lines(1:length(conc_ifo_am), conc_ifo_am, col = "blue")
lines(1:length(conc_leflu_am), conc_leflu_am, col = "green")

legend(1000, 50, legend=c("FLU", "IFO", "LEF"), col = c("red", 'blue', 'green'), lty = 1)

# End of the Median Expression Plotting.
# Beginning of the bar plot stuff, and the end of this whole thing. 

bar_plt_pts <- c(conc_ifo_unadj, 
                        conc_ifo_bm[[length(conc_ifo_bm)]], 
                        conc_ifo_am[[length(conc_ifo_am)]],
                        conc_leflu_unadj,
                        conc_leflu_bm[[length(conc_leflu_bm)]],
                        conc_leflu_am[[length(conc_leflu_am)]],
                        conc_fluc_unadj,
                        conc_fluc_bm[[length(conc_fluc_bm)]],
                        conc_fluc_am[[length(conc_fluc_am)]])

bar_plt_pts <- bar_plt_pts * 100

bar_names = c("IFO All", "IFO BM", "IFO AM",
              "LEF All", "LEF BM", "LEF AM",
              "FLU All", "FLU BM", "FLU AM")

barplot(bar_plt_pts, names = bar_names, 
        main = "Concordance by Chemical and Relation to the Median",
        xlab = "Chemical and Relation to the Median",
        ylab = "Concordance (%)",
        ylim = c(0,100),
        col = c("blue", "blue", "blue", 
                  "green", "green", "green", 
                  "red", "red", "red"))
               
       
#Fluc correlation chart between log fold changes in microarray and RNA-seq

mtx_corr_fluc <- matrix(nrow=1, ncol=2)

for (i in 1:length(fluc_mic$X)) {
  if (fluc_mic$X[i] %in% fluc_seq$X) {
    idx <- which(fluc_seq$X == fluc_mic$X[i])[1]
    mtx_corr_fluc <- rbind(mtx_corr_fluc, 
                           c(fluc_mic$logFC[i],fluc_seq$log2FoldChange[idx]))
  }
}

plot(mtx_corr_fluc, pch = 20, xlab = "Microarray Log Fold Change of DEGs",
     ylab = "RNA-seq Log Fold Change of DEGs", 
     main = "Correlation of Log Fold Changes for Fluconazole")

abline(fit <-lm(mtx_corr_fluc[,2] ~ mtx_corr_fluc[,1]))

r2_stg <- paste("R^2= ", format(summary(fit)$adj.r.squared, digits=4))

legend("topleft", bty="n", legend=r2_stg)


# Same thing but for LEF

mtx_corr_leflu <- matrix(nrow=1, ncol=2)

for (i in 1:length(leflu_mic$X)) {
  if (leflu_mic$X[i] %in% leflu_seq$X) {
    idx <- which(leflu_seq$X == leflu_mic$X[i])[1]
    mtx_corr_leflu <- rbind(mtx_corr_leflu, 
                           c(leflu_mic$logFC[i],leflu_seq$log2FoldChange[idx]))
  }
}

plot(mtx_corr_leflu, pch = 20, xlab = "Microarray Log Fold Change of DEGs",
     ylab = "RNA-seq Log Fold Change of DEGs", 
     main = "Correlation of Log Fold Changes for Leflunomide")

abline(fit <-lm(mtx_corr_leflu[,2] ~ mtx_corr_leflu[,1]))

r2_stg <- paste("R^2= ", format(summary(fit)$adj.r.squared, digits=4))

legend("topleft", bty="n", legend=r2_stg)
# The end. 