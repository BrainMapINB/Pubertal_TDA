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
# Generate one sex columns
# Compute PDS-Sex interaction manually
datos$AgeM <- datos$Age; datos$AgeM[datos$Sex=="F"] <- 0
datos$PDSeM <- datos$PDSe; datos$PDSeM[datos$Sex=="F"] <- 0

# 1) Find best model (AIC)

# Compute whole-connectome AIC
tda_var <- c("b0_auc","b1_auc")

# Load 'gamlss', 'gamm4', & ' package for GAMM analyses
if(require(gamlss)==0) install.packages("gamlss"); library(gamlss)
if(require(gamm4)==0) install.packages("gamm4"); library(gamm4)
if(require(lme4)==0) install.packages("lme4"); library(lme4)

# Apply models to TDA features
for(tt in 1:length(tda_var)){
  print(tda_var[tt])
 
  # LME model (age slope)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ Age + FDRMS + Coil + (1|BIDS.ID)"))
  fitLIN <- tryCatch(lmer(grp_form,
                          data = datos,
                          REML = F),
                     error=function(e) NA)
  print(tryCatch(AIC(fitLIN), error=function(e) NA))
 
  # LME model (age-sex slope)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ Age*Sex + FDRMS + Coil + (1|BIDS.ID)"))
  fitLINi <- tryCatch(lmer(grp_form,
                           data = datos,
                           REML = F),
                      error=function(e) NA)
  print(tryCatch(AIC(fitLINi), error=function(e) NA))
 
  # GAMM model (age spline)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(Age) + FDRMS + Coil"))
  fitGAMM <- tryCatch(gamm4(grp_form,
                            random = ~ (1|BIDS.ID),
                            data = datos,
                            REML = F),
                      error=function(e) NA)
  print(tryCatch(AIC(fitGAMM$mer), error=function(e) NA))
 
  # GAMM model (age*sex splines)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ s(Age, by=Sex) + Sex + FDRMS + Coil"))
  fitGAMMi <- tryCatch(gamm4(grp_form,
                             random = ~ (1|BIDS.ID),
                             data = datos,
                             REML = F),
                       error=function(e) NA)
  print(tryCatch(AIC(fitGAMMi$mer), error=function(e) NA))
 
  # GAMM model (PDS LOESS)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ lo(~PDSe) + FDRMS + Coil + re(random = ~1|BIDS.ID)"))
  fitPDS <- tryCatch(gamlss(grp_form, data = datos,
                            control = gamlss.control(n.cyc = 30)),
                     error=function(e) NA)
  print(tryCatch(AIC(fitPDS), error=function(e) NA))
 
  # GAMM model (PDS*sex splines)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ lo(~PDSe)+lo(~PDSeM)+Sex + FDRMS + Coil + re(random = ~1|BIDS.ID)"))
  fitPDSi <- tryCatch(gamlss(grp_form, data = datos,
                             control = gamlss.control(n.cyc = 30)),
                      error=function(e) NA)
  print(tryCatch(AIC(fitPDSi), error=function(e) NA))
}

# 2) Plot results - whole-connectome
# Apply PDS-GAMM
for(tt in 1:length(tda_var)){
 
  # Compute model
  # GAMM model (PDS LOESS)
  grp_form <- as.formula(paste0(tda_var[tt], " ~ lo(~PDSe) + FDRMS + Coil + re(random = ~1|BIDS.ID)"))
  fitPDS <- tryCatch(gamlss(grp_form, data = datos,
                            control = gamlss.control(n.cyc = 30)),
                     error=function(e) NA)
  print(tryCatch(AIC(fitPDS), error=function(e) NA))
  drop1(fitPDS)
 
  # Calculate TDA variable residuals
  res_frm <- as.formula(paste0(tda_var[tt], "~ FDRMS + Coil + (1|BIDS.ID)"))
  datos$res <- scale(residuals(lmer(res_frm, data = datos, REML = F)))
  plt_form <- as.formula("res ~ lo(~PDSe)")
  fit_plt <- getSmo(gamlss(plt_form, data = datos))
  datos$preds <- predict(fit_plt)
  datos$preds_se <- predict(fit_plt, se = T)[[2]]
  # Generate plot
  (gSC <- ggplot(data = datos, mapping = aes(x = PDSe)) +
    geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) +
    geom_point(aes(y = res, color = Sex)) +
    geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                    ymin = preds - 1.96 * preds_se,
                    linetype = NA), alpha = 0.2) +
    geom_line(aes(y = preds), size = 1.5) +
    ylab("Residuals") + xlab("PDS") + ggtitle(tda_var[tt]) +
    theme_light())
  # Save plot
  outfile <- paste0("GAMM_TDA_",tda_var[tt],"_loPDS.svg")
  #ggsave(filename = outfile, plot = gSC, device = "svg", width = 4, height = 3)
 
  # Plot against Age as well, as a supportive info
  # Calculate TDA variable residuals
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
    ylab("Residuals") + xlab("Age (y.o.)") + ggtitle(tda_var[tt]) +
    theme_light())
  # Save plot
  outfile <- paste0("GAMM_TDA_",tda_var[tt],"_sAge.svg")
  #ggsave(filename = outfile, plot = gSC, device = "svg", width = 4, height = 3)
}

# 3) Apply PDS-GAMM at the intra-network level
# Compute FDR correction and volume assignation at the Functional Network level
if(!require(neurobase)) install.packages("neurobase"); library(neurobase)
# Read atlas info
atlas <- read.csv("https://raw.githubusercontent.com/BrainMapINB/Pubertal_TDA/main/Atlas/P264/power264.csv")
n_fn <- nlevels(atlas$Network)
# Set column names
net_lab <- levels(atlas$Net); net_lab[4:5] <- net_lab[5:4]

# Apply PDS-GAMM
tda_var <- c("b0auc","b1auc")
for(tt in 1:length(tda_var)){
 
  # Create empty matrix to store results
  DF <- matrix(data = as.numeric(NA), nrow = n_fn, ncol = 14)
  rownames(DF) <- net_lab
  colnames(DF) <- c("IPnum","IPval","slope","full_AIC",
                    "PDS_df","PDS_AIC","PDS_LRT","PDS_p","PDS_pFDR",
                    "re_df","re_AIC","re_LRT","re_p","re_pFDR")
 
  # For each functional network
  for(ff in (1:n_fn)[-12]){
   
    # Select column
    fn_frm <- as.formula(paste0(tda_var[tt],"_",net_lab[ff],
                                " ~ lo(~PDSe) + FDRMS + Coil + re(random = ~1|BIDS.ID)"))
    # Apply GAMM model
    fit <- tryCatch(gamlss(fn_frm, data = datos),
                    error=function(e) NA)
   
    if(!all(is.na(fit))){
      # Store GAMLSS results via dropping single terms
      fit_drop <- drop1(fit)
      DF[ff,4] <- fit_drop[1,2] # Full Model AIC
      DF[ff,5:8] <- unlist(fit_drop[2,1:4]) # PDS term
      DF[ff,10:13] <- unlist(fit_drop[5,1:4]) # Random effects term
     
      # Extract standardized residuals
      res_frm <- as.formula(paste0(tda_var[tt],"_",net_lab[ff],
                                   "~ FDRMS + Coil + (1|BIDS.ID)"))
      res_lme <- lmer(formula = res_frm, data = datos, REML = F)
      datos$res <- scale(residuals(res_lme))
      plt_form <- as.formula("res ~ lo(~PDSe)")
      fit_plt <- getSmo(gamlss(plt_form, data = datos))
     
      # Store slope and point of inflection
      x <- seq(1,4,0.01)
      y <- predict(fit_plt, newdata = x)
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
 
 
  # Compute manual p-values for low LRT
  DF <- DF[which(!is.na(DF[,1])),]
  DF <- as.data.frame(DF)
  # p-values for PDS term
  for(ii in which(is.na(DF$PDS_p))){
    # Positive degrees of freedom
    DF$PDS_df[ii] <- abs(DF$PDS_df[ii])
    # Compute LRT-pvalue
    DF$PDS_p[ii] <- pchisq(q = DF$PDS_LRT[ii], df = DF$PDS_df[ii], lower.tail = F)
  }
  # p-values for Random-effects term
  for(ii in which(is.na(DF$re_p))){
    # Positive degrees of freedom
    DF$re_df[ii] <- abs(DF$re_df[ii])
    # Compute LRT-pvalue
    DF$re_p[ii] <- pchisq(q = DF$re_LRT[ii], df = DF$re_df[ii], lower.tail = F)
  }
  # Apply FDR
  DF$PDS_pFDR <- p.adjust(p = DF$PDS_p, method = "fdr")
  DF$re_pFDR <- p.adjust(p = DF$re_p, method = "fdr")

  # Lastly, insert LRT values into a brain volume
  # It is necesary to download this file first
  # https://github.com/BrainMapINB/Pubertal_TDA/blob/main/Atlas/P264/FuntionalNetworks_P264.nii.gz
  # PDS term
  nii <- readnii("FuntionalNetworks_P264.nii.gz")
  for(ff in (1:n_fn)[-12]){
    # Avoid 'uncertain' label
    if(ff < 12){
      nii[nii==ff] <- DF$PDS_LRT[ff]
    } else {
      nii[nii==ff] <- DF$PDS_LRT[ff-1]
    }
  }
  nii[nii==12] <- 0
  outfile <- paste0("RSN_",tda_var[tt],"_PDS-LRT")
  writenii(nim = nii, filename = outfile)
 
  # Save results
  outfile <- paste0("RSN_",tda_var[tt],".csv")
  #write.csv(x = DF, file = outfile, quote = F)
 
  # If any effect surpass the FDR correction, plot it
  if(any(DF$PDS_pFDR < 0.05)){
    # Find those
    for(ss in which(DF$PDS_pFDR < 0.05)){
     
      # Calculate TDA variable residuals
      res_frm <- as.formula(paste0(tda_var[tt],"_",(net_lab[-12])[ss], "~ FDRMS + Coil + (1|BIDS.ID)"))
      datos$res <- scale(residuals(lmer(res_frm, data = datos, REML = F)))
      plt_form <- as.formula("res ~ lo(~PDSe)")
      fit_plt <- getSmo(gamlss(plt_form, data = datos))
      datos$preds <- predict(fit_plt)
      datos$preds_se <- predict(fit_plt, se = T)[[2]]
      # Generate plot
      gSC <- ggplot(data = datos, mapping = aes(x = PDSe)) +
        #geom_line(aes(y = res, group = BIDS.ID, color = Sex), size=.25) +
        #geom_point(aes(y = res, color = Sex)) +
        geom_ribbon(aes(ymax = preds + 1.96 * preds_se,
                        ymin = preds - 1.96 * preds_se,
                        linetype = NA), alpha = 0.2) +
        geom_line(aes(y = preds), size = 1.5) +
        ylab("Residuals") + xlab("PDS") + ggtitle(tda_var[tt]) +
        theme_light()
      show(gSC)
      # Save plot
      outfile <- paste0("RSN_",tda_var[tt],"-",(net_lab[-12])[ss],"_loPDS_sct.svg")
      #ggsave(filename = outfile, plot = gSC, device = "svg", width = 4, height = 3)
    }
  }
}


############################################################################
# Extra commands
############################################################################

# Heat palette for BrainNet Viewer custom color palette
aux <- t(round(col2rgb(heat.colors(64))/255,3))
#write.table(aux, "colpalette_heat64.txt", quote = F, sep = " ", row.names = F, col.names = F)
aux <- aux[nrow(aux):1,]
#write.table(aux, "colpalette_heat64_rev.txt", quote = F, sep = " ", row.names = F, col.names = F)


