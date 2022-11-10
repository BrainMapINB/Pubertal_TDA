#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())
# Save original parameters
op <- par()

############################################################################
# Sample description
############################################################################

# Read phenotypic data
info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Phenotypic/phenotypic.csv")
# Factorize subjects' ID
info$BIDS.ID <- factor(info$BIDS.ID)

# Longitudinal scatter-plot
# Sort original variables
info <- info[order(info$Age),]
edad <- info$Age
id <- as.character(info$BIDS.ID)
# Extract labels and frequencies
id_lab <- unique(id)
id_num <- sapply(id_lab,function(x)sum(id==x))
n <- length(id_lab)
sexo <- info$Sex
qc <- abs(info$QC-1)
# Set color
if(!require(scales)) install.packages("scales"); library(scales)
huecol <- hue_pal()(2)
# Set axis frequencies
seq1y <- seq(floor(min(edad)), ceiling(max(edad)), 1)
seq2y <- seq(floor(min(edad)),ceiling(max(edad)),2)
# Plot canvas
#svg(filename = "Results/LongPlot.svg", width = 5, height = 6)
plot(edad, 1:length(edad), type = "n", axes = F,
     ylim = c(0,n), xlim = c(floor(min(edad)), ceiling(max(edad))),
     main = "Longitudinal plot", ylab = "Subject", xlab = "Age (years)")
axis(1, seq1y, F); axis(1,seq2y)
axis(2, seq(0, 10*ceiling(n/10), 10),las=1)
# Draw points
for(ii in 1:n){
  # Find sessions
  pos <- which(id==id_lab[ii])
  # Extract sex character
  ifelse(sexo[pos[1]] == "F", sexchar <- huecol[1], sexchar <- huecol[2])
  # Draw points
  if(length(pos)==1) points(edad[pos], ii, pch = 21, bg = ifelse(qc[pos],"black",sexchar), cex = 0.8, col = "black") else{
    lines(c(edad[min(pos)],edad[max(pos)]),c(ii,ii), lwd = 1, col = sexchar)
    for(jj in pos) points(edad[jj], ii, pch = 21, bg = ifelse(qc[jj],"black",sexchar), cex = 0.8, col = "black")
  }
}
# Legend
legend(min(edad), n, c("F","M","QC"), pch = 21, pt.bg = c(huecol,"#000000"), bty = "n")
#dev.off()

# Plot Pubertal Development Scale (PDS)
# Load 'ggplot2' package
if(require(ggplot2)==0) install.packages("ggplot2"); library(ggplot2)
# Load 'lme4' package
if(require(lme4)==0) install.packages("lme4"); library(lme4)

# Regress out intra-individual trends
ex_lme <- lmer(PDSe~1 + (1|BIDS.ID), data = info, REML = "F")
info$PDSr <- resid(ex_lme) + ex_lme@beta

# Plot PSD vs Age
g1 <- ggplot(data = info, mapping = aes(x = Age, y = PDSe)) +
    geom_line(aes(group = BIDS.ID, color = Sex), size=.25) +
    geom_point(aes(color = Sex)) +
    geom_smooth(aes(y=PDSe, color = Sex), size=1.5) +
    theme_bw() +
    labs(title = "Pubertal Development Scale", y="PDS")
show(g1)
#ggsave(filename = "Results/agePDS.svg", plot = g1,
#       device = "svg", width = 4, height = 3)

############################################################################
# Compute Topological Data Analysis features
############################################################################

# Clear workspace
rm(list = ls())

# Read phenotypic data
info <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Phenotypic/phenotypic.csv")
# Remove those sessions with excesive motion artifact
info <- info[which(as.logical(info$QC)),]
# Refactorize subjects' ID
info$BIDS.ID <- factor(info$BIDS.ID)

# Generate output name
tda_dir <- file.path(getwd(),"Features")
if(!dir.exists(tda_dir)) dir.create(tda_dir)
outfile <- file.path(tda_dir,"betti.csv")

# Run Rips filtration for each FC matrix
if(!file.exists(outfile)){
  
  # 1) Compute TDA features
  
  # Extract time-series files (those with GSR)
  ts_files <- paste0("https://raw.githubusercontent.com/BrainMapINB/Pubertal_Connectome/main/Brain/GSR/sub-",
                     info$BIDS.ID,"/ses-",info$Session,"/pp150v_wGSR_NIHPD_MNI_P264_ts.txt")
  
  # Compute connectivity matrices
  ConnMx <- sapply(ts_files,
                   function(x){
                     ts <- read.table(x)
                     cmx <- cor(ts)
                     if(sum(is.na(cmx))>0) cmx[is.na(cmx)] <- 0
                     return(cmx)
                   }
  )
  # Reshape object
  nroi <- sqrt(dim(ConnMx)[1])
  ConnMx_dim <- c(nroi, nroi, nrow(info))
  ConnMx <- array(ConnMx, dim = ConnMx_dim)
  
  # Generate output name
  outfile <- file.path(tda_dir, "rips.rds")
  
  # Run Rips filtration for each FC matrix
  if(!file.exists(outfile)){
    
    # Load Betti-Rips function
    source("https://raw.githubusercontent.com/BrainMapINB/Pubertal_TDA/main/Scripts/Functions/rips_betti.R")
    
    # Create empty list
    all_rips <- vector("list", nrow(info))
    
    # Apply Rips filtrations
    for(ii in 1:nrow(info)){
      print(ii)
      all_rips[[ii]] <- rips_betti(ConnMx[,,ii], 1)
    }
    # Export output
    saveRDS(all_rips, outfile)
  }
  
  # Generate 1000 permuted TDA
  null_file <- file.path(tda_dir, "rips_1000perm.rds")
  if(!file.exists(null_file)){
    # Load Betti-Rips function
    source("https://raw.githubusercontent.com/BrainMapINB/Pubertal_TDA/main/Scripts/Functions/rips_betti.R")
    # Set empty object
    null_rips <- vector("list", 1000)
    for(pp in 1:1000){
      # Print progress
      print(pp)
      # Extract one connectivity matrix
      pmx <- array(as.numeric(0), dim = ConnMx_dim[1:2])
      ptri <- (ConnMx[,,sample(1:nrow(info),1)])[upper.tri(pmx)]
      # Permute
      pmx[upper.tri(pmx)] <- ptri[sample(1:length(ptri))]
      pmx <- pmx + t(pmx)
      diag(pmx) <- 1
      # Apply TDA
      null_rips[[pp]] <- rips_betti(pmx, 1)
    }
    # Export output
    saveRDS(null_rips, null_file)
  }
  
  # Read filtration list
  rips_list <- readRDS(outfile)
  
  # Set columns for Betti# areas
  info$b1_auc <- info$b0_auc <- as.numeric(0)
  
  # Extract Betti-0/1 areas
  for(ii in 1:nrow(info)){
    # Extract diagram
    rips_diag <- rips_list[[ii]][[1]][-1,]
    # Stepwise intervals
    rips_int <- rips_diag[,3]-rips_diag[,2]
    # Store by dimension
    info$b0_auc[ii] <- sum(rips_int[which(rips_diag[,1]==0)])
    info$b1_auc[ii] <- sum(rips_int[which(rips_diag[,1]==1)])
  }
  
  # 3) Compute TDA features at the functional network level
  
  # Extract Functional Network features
  atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_TDA/main/Atlas/P264/power264.csv")
  rsn_idx <- (1:nlevels(atlas$Network))[-12]
  
  # Load Betti-Rips function
  source("https://raw.githubusercontent.com/BrainMapINB/Pubertal_TDA/main/Scripts/Functions/rips_betti.R")
  
  # Create empty objects
  FN_b0_auc <- FN_b1_auc <- matrix(0, nrow = nrow(info), ncol = length(rsn_idx))
  
  # Compute TDA features for Functional Networks
  for(pp in 1:length(rsn_idx)){
    # Print progress
    print(pp)
    
    # Extract network positions
    net_idx <- which(atlas$Network==levels(atlas$Network)[rsn_idx[pp]])
    
    # Create sub-connectivity matrix
    ConnMx_sub <- ConnMx[net_idx,net_idx,]
    
    # Compute Rips diagram
    for(ii in 1:nrow(info)){
      rips_diag <- rips_betti(ConnMx_sub[,,ii],1)[[1]][-1,]
      # Stepwise intervals
      rips_int <- rips_diag[,3]-rips_diag[,2]
      # Store by dimension
      FN_b0_auc[ii,pp] <- sum(rips_int[which(rips_diag[,1]==0)])
      if(sum(rips_diag[,1])>0){
        FN_b1_auc[ii,pp] <- sum(rips_int[which(rips_diag[,1]==1)])
      }
    }
  }
  
  # Set column names
  colnames(FN_b0_auc) <- paste0("b0auc_",levels(atlas$Net)[rsn_idx])
  colnames(FN_b1_auc) <- paste0("b1auc_",levels(atlas$Net)[rsn_idx])
  
  # Concatenate outputs
  info <- cbind(info, FN_b0_auc, FN_b1_auc)
  # Store results
  outfile <- file.path(tda_dir,"betti.csv")
  write.csv(info, outfile, row.names = F)
}

# Plot sample averages
# Read data
datos <- read.csv(file.path(tda_dir,"betti.csv"))
rips_list <- readRDS(file.path(tda_dir, "rips.rds"))
null_list <- readRDS(file.path(tda_dir, "rips_1000perm.rds"))

# Extract Betti-0
outfile <- file.path(tda_dir,"b0avg.rds")
if(!file.exists(outfile)){
  nnodes <- sum(rips_list[[1]][[1]][,1]==0)
  b0 <- t(sapply(1:nrow(datos), function(x) rips_list[[x]][[1]][1:nnodes,3]))
  # Betti-0 intercept and 95% CI
  b0_avg <- array(as.numeric(0), dim = c(nnodes,4))
  b0_avg[1,1] <- 1
  for(ii in 2:nnodes){
    datos$y <- b0[,ii]
    fit <- summary(lmer(y ~ 1 + (1|BIDS.ID), data = datos))
    b0_avg[ii,1:2] <- fit$coefficients[1:2]
  }
  # Extract null distribution too
  b0 <- t(sapply(1:length(null_list), function(x) null_list[[x]][[1]][1:nnodes,3]))
  b0_avg[,3] <- colMeans(b0)
  b0_avg[,4] <- apply(b0,2,sd)
  # Add initial filtration
  b0_avg <- rbind(b0_avg,rep(0,ncol(b0_avg)))
  b0_avg <- b0_avg[-1,]
  # Save sample averages
  saveRDS(b0_avg, outfile)
}
# Read Betti-0 sample intercept
b0_avg <- readRDS(outfile)

# Plot
#pdf("Results/TDA_b0_avg.pdf", 5,3)
plot(b0_avg[,1], 1:nrow(b0_avg), type = "n", xlim = c(0,1), las = 1, frame.plot = F,
     ylab = "Betti-0", xlab = "Filtration Value", main = "Betti-0 sample intercept")
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2, las = 1)
# Null dist
polygon(c(b0_avg[,3]-b0_avg[,4]*1.96, rev(b0_avg[,3]+b0_avg[,4]*1.96)),
        c(1:nrow(b0_avg),nrow(b0_avg):1),
        col = scales::alpha("gray",0.35), border = F)
lines(b0_avg[,3], 1:nrow(b0_avg), col = "gray", lwd = 1)
# Observed
polygon(c(b0_avg[,1]-b0_avg[,2]*1.96, rev(b0_avg[,1]+b0_avg[,2]*1.96)),
        c(1:nrow(b0_avg),nrow(b0_avg):1),
        col = scales::alpha("blue",0.35), border = F)
lines(b0_avg[,1], 1:nrow(b0_avg), col = "darkblue", lwd = 1)
# Legend
legend("topright", legend = c("Observed","Random"), col = c("darkblue","gray"), lwd = 2)
#dev.off()

# Print the sample Betti-0 equal to one
print(b0_avg[1,1]) # Betti-0 LME intercept (mean de facto)
print(b0_avg[1,2]) # Betti-0 standard error
print(b0_avg[1,2]*sqrt(nrow(info))) # Betti-0 standard deviation
print(b0_avg[1,2]*1.96) # Betti-0 95%-CI

# Add lines
abline(v=b0_avg[1,1], col="blue",lwd=2)
abline(v=b0_avg[1,1]-b0_avg[1,2]*1.96, col="red",lwd=2,lty=2)
abline(v=b0_avg[1,1]+b0_avg[1,2]*1.96, col="red",lwd=2,lty=2)
#dev.off()

# Betti-1
tseq <- seq(0, 1, length = 1000)
outfile <- file.path(tda_dir,"b1avg.rds")
if(!file.exists(outfile)){
  if(!require(TDA)) install.packages("TDA"); library(TDA)
  b1 <- t(sapply(1:nrow(datos), function(x) landscape(rips_list[[x]][[1]], tseq = tseq)))
  # Betti-0 intercept and 95% CI
  b1_avg <- array(as.numeric(0), dim = c(1000,4))
  for(ii in 1:1000){
    if(sum(b1[,ii])>0){
      datos$y <- b1[,ii]
      fit <- tryCatch(lmer(y ~ 1 + (1|BIDS.ID), data = datos), error=function(e) NA)
      if(!is.na(fit)){
        fit <- summary(lmer(y ~ 1 + (1|BIDS.ID), data = datos))
        b1_avg[ii,1:2] <- fit$coefficients[1:2]
      }
    }
  }
  # Null distribution Betti-1
  b1 <- t(sapply(1:length(null_list), function(x) landscape(null_list[[x]][[1]], tseq = tseq)))
  b1_avg[,3] <- colMeans(b1)
  b1_avg[,4] <- apply(b1,2,sd)
  b1_avg <- b1_avg*length(tseq)
  # Save sample averages
  saveRDS(b1_avg, outfile)
}
b1_avg <- readRDS(outfile)

# Plot
#pdf("Results/TDA_b1_avg.pdf", 5,3)
plot(x = tseq, y =  b1_avg[,3], type = "n", xlim = c(0,1),
     ylim = c(0,250), las = 1, frame.plot = F,
     ylab = "Betti-1", xlab = "Filtration Value",
     main = "Betti-1 sample intercept")
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2, las = 1)
# Null dist
polygon(c(tseq,rev(tseq)),
        c(b1_avg[,3]-b1_avg[,4]*1.96, rev(b1_avg[,3]+b1_avg[,4]*1.96)),
        col = scales::alpha("gray",0.35), border = F)
lines(tseq, b1_avg[,3], col = "gray", lwd = 1)
# Observed
polygon(c(tseq,rev(tseq)),
        c(b1_avg[,1]-b1_avg[,2]*1.96, rev(b1_avg[,1]+b1_avg[,2]*1.96)),
        col = scales::alpha("blue",0.35), border = F)
lines(tseq, b1_avg[,1], col = "darkblue", lwd = 1)
# Legend
#legend("topright", legend = c("OBS","PERM"), col = c("darkblue","gray"), lwd = 2)
#dev.off()

# Compute the sample intercept for the first and last hole of Betti-1
datos$b1_last <- datos$b1_one <- as.numeric(NA)
for(rr in 1:length(rips_list)){
  aux_diag <- rips_list[[rr]]$diagram
  aux_b1 <- which(aux_diag[,1]==1)
  datos$b1_one[rr] <- min(aux_diag[aux_b1,2])
  datos$b1_last[rr] <- max(aux_diag[aux_b1,3])
}
# Extract intercepts
# First hole
fit <- summary(lmer(b1_one ~ 1 + (1|BIDS.ID), data = datos))
fit$coefficients[1:2] # Intercept and s.e.
fit$coefficients[2]*sqrt(nrow(datos)) #s.d.
fit$coefficients[2]*1.96 # 95%-CI
# Add lines
abline(v=fit$coefficients[1], col="blue",lwd=2)
abline(v=fit$coefficients[1]-fit$coefficients[2]*1.96, col="red",lwd=2,lty=2)
abline(v=fit$coefficients[1]+fit$coefficients[2]*1.96, col="red",lwd=2,lty=2)

# Last hole
fit <- summary(lmer(b1_last ~ 1 + (1|BIDS.ID), data = datos))
fit$coefficients[1:2] # Intercept and s.e.
fit$coefficients[2]*sqrt(nrow(datos)) #s.d.
fit$coefficients[2]*1.96 # 95%-CI
# Add lines
abline(v=fit$coefficients[1], col="blue",lwd=2)
abline(v=fit$coefficients[1]-fit$coefficients[2]*1.96, col="red",lwd=2,lty=2)
abline(v=fit$coefficients[1]+fit$coefficients[2]*1.96, col="red",lwd=2,lty=2)
#dev.off()

# Extended material
# Compute the relationship between B1-AUC and number of holes
datos$b1_holes <- sapply(1:length(rips_list), function(x) sum(rips_list[[x]]$diagram[,1]==1))
# Regress out intra-individual trends
# B1-AUC
ex_lme <- lmer(b1_auc~1 + (1|BIDS.ID), data = datos, REML = "F")
datos$b1_auc_e <- resid(ex_lme) + ex_lme@beta
# B1-Holes
ex_lme <- lmer(b1_holes~1 + (1|BIDS.ID), data = datos, REML = "F")
datos$b1_holes_e <- resid(ex_lme) + ex_lme@beta
# Plot
summary(lm(b1_holes_e~b1_auc_e, datos))
# linear trend + confidence interval
p3 <- ggplot(datos, aes(x=b1_auc_e, y=b1_holes_e)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_bw() +
  labs(title = "Betti-1 AUC vs Holes",
       y="# Holes", x="B1-AUC")
show(p3)
#ggsave(filename = "B1_AUC-Holes.pdf", plot = p3, device = "pdf",
#       path = "Results/", scale = 1, width = 3, height = 3)
