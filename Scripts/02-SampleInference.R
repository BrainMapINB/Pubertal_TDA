#!/usr/bin/env Rscript

# Clear workspace
rm(list = ls())
# Save original parameters
op <- par()

# Read phenotypic data
datos <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_TDA/main/Features/betti.csv")

############################################################################
# Sample inference
############################################################################

# Remove columns with NA's
datos <- datos[,which(apply(datos,2,function(x) sum(is.na(x)))==0)]
# Factorize character variables
datos$BIDS.ID <- factor(datos$BIDS.ID)
datos$Sex <- factor(datos$Sex)
datos$Coil <- factor(datos$Coil)

# 1) Find best model (AIC)

# Compute whole-connectome AIC
tda_var <- c("b0_auc","b1_auc")

# Load 'lme4', 'gamm4', & 'ggplot2' packages
if(require(lme4)==0) install.packages("lme4"); library(lme4)
if(require(gamm4)==0) install.packages("gamm4"); library(gamm4)
if(require(ggplot2)==0) install.packages("ggplot2"); library(ggplot2)

# Create summary table
stab <- matrix(data = as.numeric(NA), nrow = 8, ncol = 10)
colnames(stab) <- paste0(rep(tda_var, each = 5),"_", c("AIC","resid-SW","resid-p","RE-SW","RE-p"))

# Apply models to TDA features
for(tt in 1:length(tda_var)){
  print(tda_var[tt])
 
  # LME model (age slope)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ Age + FDRMS + Coil + (1|BIDS.ID)"))
  fitLIN <- lmer(formula = grp_form, data = datos, REML = F)
  stab[1,1+5*(tt-1)] <- AIC(fitLIN)
  swt <- shapiro.test(resid(fitLIN))
  stab[1,2+5*(tt-1)] <- swt$statistic
  stab[1,3+5*(tt-1)] <- swt$p.value
  swt <- shapiro.test(ranef(fitLIN)[[1]][[1]])
  stab[1,4+5*(tt-1)] <- swt$statistic
  stab[1,5+5*(tt-1)] <- swt$p.value
  
  # LME model (age-sex slope)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ Age*Sex + FDRMS + Coil + (1|BIDS.ID)"))
  fitLINi <- lmer(formula = grp_form, data = datos, REML = F)
  stab[2,1+5*(tt-1)] <- AIC(fitLINi)
  swt <- shapiro.test(resid(fitLINi))
  stab[2,2+5*(tt-1)] <- swt$statistic
  stab[2,3+5*(tt-1)] <- swt$p.value
  swt <- shapiro.test(ranef(fitLINi)[[1]][[1]])
  stab[2,4+5*(tt-1)] <- swt$statistic
  stab[2,5+5*(tt-1)] <- swt$p.value
  
  # GAMM model (age spline)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(Age) + FDRMS + Coil"))
  fitGAMM <- gamm4(formula = grp_form, random = ~ (1|BIDS.ID),
                   data = datos, REML = F)
  stab[3,1+5*(tt-1)] <- AIC(fitGAMM$mer)
  swt <- shapiro.test(resid(fitGAMM$mer))
  stab[3,2+5*(tt-1)] <- swt$statistic
  stab[3,3+5*(tt-1)] <- swt$p.value
  swt <- shapiro.test(ranef(fitGAMM$mer)[[1]][[1]])
  stab[3,4+5*(tt-1)] <- swt$statistic
  stab[3,5+5*(tt-1)] <- swt$p.value
  
  # GAMM model (age*sex splines)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(Age, by=Sex) + Sex + FDRMS + Coil"))
  fitGAMMi <- gamm4(formula = grp_form, random = ~ (1|BIDS.ID),
                    data = datos, REML = F)
  stab[4,1+5*(tt-1)] <- AIC(fitGAMMi$mer)
  swt <- shapiro.test(resid(fitGAMMi$mer))
  stab[4,2+5*(tt-1)] <- swt$statistic
  stab[4,3+5*(tt-1)] <- swt$p.value
  swt <- shapiro.test(ranef(fitGAMMi$mer)[[1]][[1]])
  stab[4,4+5*(tt-1)] <- swt$statistic
  stab[4,5+5*(tt-1)] <- swt$p.value
  
  # GAMM model (PDS spline)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(PDSe) + FDRMS + Coil"))
  fitGAMM <- gamm4(formula = grp_form, random = ~ (1|BIDS.ID),
                   data = datos, REML = F)
  stab[5,1+5*(tt-1)] <- AIC(fitGAMM$mer)
  swt <- shapiro.test(resid(fitGAMM$mer))
  stab[5,2+5*(tt-1)] <- swt$statistic
  stab[5,3+5*(tt-1)] <- swt$p.value
  swt <- shapiro.test(ranef(fitGAMM$mer)[[1]][[1]])
  stab[5,4+5*(tt-1)] <- swt$statistic
  stab[5,5+5*(tt-1)] <- swt$p.value
  
  # GAMM model (PDS*sex splines)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(PDSe, by=Sex) + Sex + FDRMS + Coil"))
  fitGAMMi <- gamm4(formula = grp_form, random = ~ (1|BIDS.ID),
                    data = datos, REML = F)
  stab[6,1+5*(tt-1)] <- AIC(fitGAMMi$mer)
  swt <- shapiro.test(resid(fitGAMMi$mer))
  stab[6,2+5*(tt-1)] <- swt$statistic
  stab[6,3+5*(tt-1)] <- swt$p.value
  swt <- shapiro.test(ranef(fitGAMMi$mer)[[1]][[1]])
  stab[6,4+5*(tt-1)] <- swt$statistic
  stab[6,5+5*(tt-1)] <- swt$p.value
  
  # GAMM model (Age+PDS spline)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(Age,PDSe) + FDRMS + Coil"))
  fitGAMM <- gamm4(formula = grp_form, random = ~ (1|BIDS.ID),
                   data = datos, REML = F)
  stab[7,1+5*(tt-1)] <- AIC(fitGAMM$mer)
  swt <- shapiro.test(resid(fitGAMM$mer))
  stab[7,2+5*(tt-1)] <- swt$statistic
  stab[7,3+5*(tt-1)] <- swt$p.value
  swt <- shapiro.test(ranef(fitGAMM$mer)[[1]][[1]])
  stab[7,4+5*(tt-1)] <- swt$statistic
  stab[7,5+5*(tt-1)] <- swt$p.value
  
  # GAMM model (Age+PDS spline)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(Age,PDSe, by=Sex) + Sex + FDRMS + Coil"))
  fitGAMMi <- gamm4(formula = grp_form, random = ~ (1|BIDS.ID),
                    data = datos, REML = F)
  stab[8,1+5*(tt-1)] <- AIC(fitGAMMi$mer)
  swt <- shapiro.test(resid(fitGAMMi$mer))
  stab[8,2+5*(tt-1)] <- swt$statistic
  stab[8,3+5*(tt-1)] <- swt$p.value
  swt <- shapiro.test(ranef(fitGAMMi$mer)[[1]][[1]])
  stab[8,4+5*(tt-1)] <- swt$statistic
  stab[8,5+5*(tt-1)] <- swt$p.value
  
}

# Export results
stab <- as.data.frame(stab)
stab <- cbind(c("LME-Age","LME-Age.Sex",
                "GAMM-Age","GAMM-Age.Sex",
                "GAMM-PDS","GAMM-PDS.Sex",
                "GAMM-AgePDS","GAMM-AgePDS.Sex"),
              stab)
names(stab)[1] <- "Model"
#write.csv(stab, "Results/ModelStats.csv",row.names = F)               

# 2) Plot results - whole-connectome
# Apply PDS-GAMM
for(tt in 1:length(tda_var)){
 
  # Compute model
  # GAMM model (PDS spline)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(PDSe) + FDRMS + Coil"))
  fitPDS <- gamm4(formula = grp_form, random = ~ (1|BIDS.ID),
                   data = datos, REML = F)
  print(tryCatch(AIC(fitPDS$mer), error=function(e) NA))
  stab <- summary(fitPDS$gam)$s.table
  show(stab)
  
  # Calculate TDA variable residuals
  datos$res <- datos[[tda_var[tt]]]
  plt_form <- as.formula("res ~ s(PDSe)")
  fit_plt <- gamm4(plt_form, data = datos)
  datos$preds <- predict(fit_plt$gam)
  datos$preds_se <- predict(fit_plt$gam, se = T)[[2]]
  # Generate plot
  (gSC <- ggplot(data = datos, mapping = aes(x = PDSe)) +
      geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) +
      geom_point(aes(y = res, color = Sex)) +
      geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                      ymin = preds - 1.96 * preds_se,
                      linetype = NA), alpha = 0.2) +
      geom_line(aes(y = preds), size = 1.5) +
      ylab(paste0(tda_var[tt]," (a.u.)")) + xlab("PDS") +
      ggtitle(paste0("F(",round(stab[2],2),")=",round(stab[3],2),"; p=",signif(stab[4],2))) +
      theme_light())
  # Save plot
  outfile <- paste0("Results/GAMM_TDA_",tda_var[tt],"_raw_sPDS.pdf")
  #ggsave(filename = outfile, plot = gSC, device = "pdf", width = 4, height = 3)
  
  # Calculate TDA variable residuals (for extended figure)
  res_frm <- as.formula(paste0(tda_var[tt], "~ FDRMS + Coil + (1|BIDS.ID)"))
  datos$res <- scale(residuals(lmer(res_frm, data = datos, REML = F)))
  plt_form <- as.formula("res ~ s(PDSe)")
  fit_plt <- gamm4(plt_form, data = datos)
  datos$preds <- predict(fit_plt$gam)
  datos$preds_se <- predict(fit_plt$gam, se = T)[[2]]
  # Generate plot
  (gSC <- ggplot(data = datos, mapping = aes(x = PDSe)) +
    geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) +
    geom_point(aes(y = res, color = Sex)) +
    geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                    ymin = preds - 1.96 * preds_se,
                    linetype = NA), alpha = 0.2) +
    geom_line(aes(y = preds), size = 1.5) +
    ylab(paste0(tda_var[tt]," (e)")) + xlab("PDS") +
    ggtitle(paste0("F(",round(stab[2],2),")=",round(stab[3],2),"; p=",signif(stab[4],2))) +
    theme_light())
  # Save plot
  outfile <- paste0("Results/GAMM_TDA_",tda_var[tt],"_resid_sPDS.pdf")
  #ggsave(filename = outfile, plot = gSC, device = "pdf", width = 4, height = 3)
 
  # Plot against Age as well, as a supportive info
  # GAMM model (Age spline)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(Age) + FDRMS + Coil"))
  fitAge <- gamm4(formula = grp_form, random = ~ (1|BIDS.ID),
                  data = datos, REML = F)
  print(tryCatch(AIC(fitAge$mer), error=function(e) NA))
  stab <- summary(fitAge$gam)$s.table
  show(stab)
  
  # Calculate TDA variable raw scores
  datos$res <- datos[[tda_var[tt]]]
  plt_form <- as.formula("res ~ + s(Age)")
  fit_plt <- gamm4(plt_form, data = datos)
  datos$preds <- predict(fit_plt$gam)
  datos$preds_se <- predict(fit_plt$gam, se = T)[[2]]
  # Generate plot
  (gSC <- ggplot(data = datos, mapping = aes(x = Age)) +
      geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) +
      geom_point(aes(y = res, color = Sex)) +
      geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                      ymin = preds - 1.96 * preds_se,
                      linetype = NA), alpha = 0.2) +
      geom_line(aes(y = preds), size = 1.5) +
      ylab(paste0(tda_var[tt]," (a.u.)")) + xlab("Age") +
      ggtitle(paste0("F(",round(stab[2],2),")=",round(stab[3],2),"; p=",signif(stab[4],2))) +
      theme_light())
  # Save plot
  outfile <- paste0("Results/GAMM_TDA_",tda_var[tt],"_raw_sAge.pdf")
  #ggsave(filename = outfile, plot = gSC, device = "pdf", width = 4, height = 3)
  
  # Calculate TDA variable residuals (for extended figure)
  plt_form <- as.formula("res ~ + s(Age)")
  fit_plt <- gamm4(plt_form, data = datos)
  datos$preds <- predict(fit_plt$gam)
  datos$preds_se <- predict(fit_plt$gam, se = T)[[2]]
  # Generate plot
  (gSC <- ggplot(data = datos, mapping = aes(x = Age)) +
    geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) +
    geom_point(aes(y = res, color = Sex)) +
    geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                    ymin = preds - 1.96 * preds_se,
                    linetype = NA), alpha = 0.2) +
    geom_line(aes(y = preds), size = 1.5) +
      ylab(paste0(tda_var[tt]," (e)")) + xlab("Age") +
      ggtitle(paste0("F(",round(stab[2],2),")=",round(stab[3],2),"; p=",signif(stab[4],2))) +
      theme_light())
  # Save plot
  outfile <- paste0("Results/GAMM_TDA_",tda_var[tt],"_resid_sAge.pdf")
  #ggsave(filename = outfile, plot = gSC, device = "pdf", width = 4, height = 3)
}

# 3) Apply PDS-GAMM at the intra-network level
# Compute FDR correction and volume assignation at the Functional Network level
if(!require(neurobase)) install.packages("neurobase"); library(neurobase)
# Read atlas info
atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_TDA/main/Atlas/P264/power264.csv")
n_fn <- length(table(atlas$Network))
# Set column names
net_lab <- names(table(atlas$Net)); net_lab[4:5] <- net_lab[5:4]

# Apply PDS-GAMM
tda_var <- c("b0auc","b1auc")
for(tt in 1:length(tda_var)){
 
  # Create empty matrix to store results
  DF <- matrix(data = as.numeric(NA), nrow = n_fn, ncol = 7)
  rownames(DF) <- net_lab
  colnames(DF) <- c("IPnum","IPval","slope",
                    "PDS_edf","PDS_F","PDS_p","PDS_pFDR")
 
  # For each functional network
  for(ff in (1:n_fn)[-12]){
   
    # Select column
    fn_frm <- as.formula(paste0(tda_var[tt],"_",net_lab[ff],
                                " ~ s(PDSe) + FDRMS + Coil"))
    # Apply GAMM model
    fit <- tryCatch(gamm4(formula = fn_frm,
                          random = ~ (1|BIDS.ID),
                          data = datos,
                          REML = F),
                    error=function(e) NA)
    
    if(!all(is.na(fit))){
      # Store GAMM results
      fit_s <- summary(fit$gam)$s.table
      DF[ff,4] <- fit_s[1] # Estimated degrees of freedom
      DF[ff,5] <- fit_s[3] # F-statistic
      DF[ff,6] <- fit_s[4] # p-value
     
      # Extract standardized residuals
      res_frm <- as.formula(paste0(tda_var[tt],"_",net_lab[ff],
                                   "~ FDRMS + Coil + (1|BIDS.ID)"))
      res_lme <- lmer(formula = res_frm, data = datos, REML = F)
      datos$res <- scale(residuals(res_lme))
      plt_form <- as.formula("res ~ s(PDSe)")
      fit_plt <- gamm4(plt_form, data = datos)
     
      # Store slope and point of inflection
      datos$preds <- predict(fit_plt$gam)
      lo <- loess(preds~PDSe, data = datos)
      x <- seq(1,4,0.01)
      y <- predict(lo,x)
      #plot(x,y, type = "l")
      infl <- c(FALSE, diff(diff(y)>0)!=0)
      DF[ff,1] <- sum(infl)
      DF[ff,2] <- 0
      if(DF[ff,1]>0){
        # Store the last inflection point
        DF[ff,2] <- x[infl][DF[ff,1]]
        # Compute slope after inflection point
        y <- y[x>=DF[ff,2]]
        x <- x[x>=DF[ff,2]]
      }
      # Compute slope
      DF[ff,3] <- coefficients(lm(y~x))[2]
    }
  }# for(ff in (1:n_fn)
 
  # Apply FDR
  DF <- DF[which(!is.na(DF[,1])),]
  DF <- as.data.frame(DF)
  DF$PDS_pFDR <- p.adjust(p = DF$PDS_p, method = "fdr")
  
  # Lastly, insert F values into a brain volume
  # It is necessary to download this file first
  # https://github.com/BrainMapINB/Pubertal_TDA/blob/main/Atlas/P264/FuntionalNetworks_P264.nii.gz
  # PDS term
  nii <- readnii("Atlas/P264/FuntionalNetworks_P264.nii.gz")
  for(ff in (1:n_fn)[-12]){
    # Avoid 'uncertain' label
    if(ff < 12){
      nii[nii==ff] <- DF$PDS_F[ff]
    } else {
      nii[nii==ff] <- DF$PDS_F[ff-1]
    }
  }
  nii[nii==12] <- 0
  outfile <- paste0("RSN_",tda_var[tt],"_PDS-F")
  #writenii(nim = nii, filename = file.path("Results",outfile))
 
  # Save results
  outfile <- paste0("RSN_",tda_var[tt],".csv")
  #write.csv(x = DF, file = file.path("Results",outfile), quote = F)
 
  # If any effect surpass the FDR correction, plot it
  if(any(DF$PDS_pFDR < 0.05)){
    # Find those
    for(ss in which(DF$PDS_pFDR < 0.05)){
     
      # Calculate TDA variable residuals
      res_frm <- as.formula(paste0(tda_var[tt],"_",(net_lab[-12])[ss], "~ FDRMS + Coil + (1|BIDS.ID)"))
      datos$res <- scale(residuals(lmer(res_frm, data = datos, REML = F)))
      plt_form <- as.formula("res ~ s(PDSe)")
      fit_plt <- gamm4(plt_form, data = datos)
      datos$preds <- predict(fit_plt$gam)
      datos$preds_se <- predict(fit_plt$gam, se = T)[[2]]
      # Generate plot
      # First, for the residuals (extended figure)
      ttl <- paste0("(F=",round(DF$PDS_F[ss],2),
                    ", EDF=",round(DF$PDS_edf[ss],2),
                    ", p=",signif(DF$PDS_p[ss],2),")")
      gSC <- ggplot(data = datos, mapping = aes(x = PDSe)) +
        geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) +
        geom_point(aes(y = res, color = Sex)) +
        geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                        ymin = preds - 1.96 * preds_se,
                        linetype = NA), alpha = 0.2) +
        geom_line(aes(y = preds), size = 1.5) +
        ylab(paste0(tda_var[tt]," (e)")) + 
        xlab("PDS") + 
        ggtitle(names(table(atlas$Network))[ss], ttl) +
        theme_light()
      show(gSC)
      # Save plot
      outfile <- paste0("Results/RSN_",tda_var[tt],"_res-",(net_lab[-12])[ss],"_sPDS.pdf")
      #ggsave(filename = outfile, plot = gSC, device = "pdf", width = 4, height = 3)
      
      # Then for the raw scores
      datos$res <- datos[[paste0(tda_var[tt],"_",(net_lab[-12])[ss])]]
      plt_form <- as.formula("res ~ + s(PDSe)")
      fit_plt <- gamm4(plt_form, data = datos)
      datos$preds <- predict(fit_plt$gam)
      datos$preds_se <- predict(fit_plt$gam, se = T)[[2]]
      # Plot
      gSC <- ggplot(data = datos, mapping = aes(x = PDSe)) +
        geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) +
        geom_point(aes(y = res, color = Sex)) +
        geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                        ymin = preds - 1.96 * preds_se,
                        linetype = NA), alpha = 0.2) +
        geom_line(aes(y = preds), size = 1.5) +
        ylab(paste0(tda_var[tt]," (a.u.)")) + 
        xlab("PDS") + 
        ggtitle(names(table(atlas$Network))[ss], ttl) +
        theme_light()
      show(gSC)
      # Save plot
      outfile <- paste0("Results/RSN_",tda_var[tt],"_raw-",(net_lab[-12])[ss],"_sPDS.pdf")
      #ggsave(filename = outfile, plot = gSC, device = "pdf", width = 4, height = 3)
      
    }
  }
}

############################################################################
# Extra commands
############################################################################

# Heat palette for BrainNet Viewer custom color palette
aux <- t(round(col2rgb(heat.colors(64))/255,3))
#write.table(aux, "Results/colpalette_heat64.txt", quote = F, sep = " ", row.names = F, col.names = F)
