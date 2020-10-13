# ======================================================================================================
# Size-spectra analysis
# Paul Carvalho - paulcarvalho@uri.edu or pgcarvalh@gmail.com
# Updated: OCtober 13, 2020
#
# ======================================================================================================

# Set up ===============================================================================================
rm(list=ls())

# Directories
setwd("C:/Users/pgcar/Google Drive/Paul Carvalho/dissertation/chapter 2/analysis/size_spectra")

# Functions
source("functions.r")

# Libraries
library(plotrix)
library(ggplot2)
library(splitstackshape)
library(data.table)
library(tidyverse)
library(mgcv)
library(mgcViz)
library(MuMIn)
library(readxl)
library(varhandle)
library(visreg)
library(devtools)
devtools::install_github("andrew-edwards/sizeSpectra")
devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
library(sizeSpectra)
library(ggpubr)
library(rstatix)
library(fitdistrplus)
library(stringr)
library(dplyr)

# Load data
fish.df <- read.csv("../fish_spectra_data.csv")
fish.df <- fish.df[,-1] 
covariates.df <- read.csv("../covariate_spectra_data.csv")

# Subset data by region
ra <- fish.df %>% filter(region == "raja_ampat") 
wa <- fish.df %>% filter(region == "wakatobi")
lo <- fish.df %>% filter(region == "lombok")

# Subset by trophic group
carn.df <- fish.df %>% filter(tp == "Carnivore")
herb.df <- fish.df %>% filter(tp == "Herbivore")

# Table S1
table_s1 <- fish.df %>%
  dplyr::select(genus_species, functional_group, tp) %>%
  distinct(.) %>%
  mutate(functional_group = str_to_sentence(functional_group))
# write.csv(table_s1,"table_s1.csv")

# Check differences in slope between divers ============================================================
# Plot size spectra slopes for carnivores and herbivores at each regions for each observer
observerFunc_plot <- slope_regAndFunc(carn.df, herb.df)

# Plot size spectra slopes (all fishes) at each region for each observer
observerReg_plot <- slope_reg()

# Remove observer "ch" (initials PS) from analysis due to different overall size spectrum slope with other divers in Lombok
fish.df <- fish.df %>% filter(observer != "ch")    
lo <- lo %>% filter(observer != "ch")    
carn.df <- carn.df %>% filter(observer != "ch")
herb.df <- herb.df %>% filter(observer != "ch")

# Set MLE parameters ===================================================================================
# raja ampat (ra) biomass
ra.input <- set.params(ra$biomass_kg)
# wakatobi (wa) biomass
wa.input <- set.params(wa$biomass_kg)
# lombok (lo) biomass
lo.input <- set.params(lo$biomass_kg)

mgpVals <- c(1.6,0.5,0) # mgp values 2.0, 0.5, 0
xLim <- 10^par("usr")[1:2]
yLim <- 10^par("usr")[3:4]

# MLE Raja Ampat biomass ===============================================================================
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.ra <- mle_b(region="raja_ampat", x=ra.input$biomass, log_x=ra.input$log.biomass, sum_log_x=ra.input$sum.log.biomass,
                 x_min=ra.input$min.biomass, x_max=ra.input$max.biomass)
PLB.bMLE.ra.b <- PLB.return.ra[[1]] 
PLB.minLL.ra.b <- PLB.return.ra[[2]]

# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.ra.b <- PLB.minLL.ra.b$minimum
x <- ra.input$biomass
rax.PLB = seq(min(ra.input$biomass), max(ra.input$biomass), length=1000) # x values to plot PLB. Note
                                                                         # that these encompass the data, and are not based
                                                                         # on the binning (in MEE Figure 6 the line starts as
                                                                         # min(x), not the first bin.
ray.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.ra.b, xmin = min(x.PLB),
    xmax = max(x.PLB))) * length(ra.input$biomass)
spectra.text <- as.character(round(PLB.bMLE.ra.b, 2))
rab_plot <- ggplot() +
  geom_point(aes(x = (sort(ra.input$biomass, decreasing=TRUE)), y = (1:length(ra.input$biomass))),
             color = "#E69F00", size = 2, alpha = 0.3) +
  xlab(expression(paste("Body sizes, ", italic("x"), " (kg)"))) +
  ylab(expression(paste("Number of body sizes", " ">=" ", italic("x")))) +
  xlim((c(ra.input$min.biomass, ra.input$max.biomass))) +
  ylim((c(1, length(ra.input$biomass)))) +
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10', breaks = c(0, 1, 5, 10)) +
  geom_line(aes(x = rax.PLB, y = ray.PLB), col = 'black', lwd = 1) +
  annotate("text", x=0.0722, y=15, label="Raja Ampat") +
  annotate("text", x=0.07, y=5, label = expression(paste(italic("b = "), -1.58))) +
  theme_classic()

# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.ra.b - 0.5, PLB.bMLE.ra.b + 0.5, 0.00001) 
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
  PLB.LLvals[i] = negLL.PLB(bvec[i], x=ra.input$biomass, n=length(ra.input$biomass), xmin=ra.input$min.biomass,
  xmax=ra.input$max.biomass, sumlogx=ra.input$sum.log.biomass)   
}
critVal = PLB.minNegLL.ra.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
rabIn95 <- c(min(bIn95), max(bIn95))

# MLE Wakatobi biomass ===============================================================================

# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.wa <- mle_b(region="wakatobi", x=wa.input$biomass, log_x=wa.input$log.biomass, sum_log_x=wa.input$sum.log.biomass,
                 x_min=wa.input$min.biomass, x_max=wa.input$max.biomass)
PLB.bMLE.wa.b <- PLB.return.wa[[1]] 
PLB.minLL.wa.b <- PLB.return.wa[[2]]

# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.wa.b <- PLB.minLL.wa.b$minimum
x <- wa.input$biomass
wax.PLB = seq(min(wa.input$biomass), max(wa.input$biomass), length=1000) # x values to plot PLB. Note
                                                                         # that these encompass the data, and are not based
                                                                         # on the binning (in MEE Figure 6 the line starts as
                                                                         # min(x), not the first bin.
way.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.wa.b, xmin = min(x.PLB),
    xmax = max(x.PLB))) * length(wa.input$biomass)
spectra.text <- as.character(round(PLB.bMLE.wa.b, 2))
wab_plot <- ggplot() +
  geom_point(aes(x = (sort(wa.input$biomass, decreasing=TRUE)), y = (1:length(wa.input$biomass))), 
             color = "#56B4E9", size = 2, alpha = 0.3) +
  xlab(expression(paste("Body sizes, ", italic("x"), " (kg)"))) +
  ylab(expression(paste("Number of body sizes", " ">=" ", italic("x")))) +
  xlim((c(wa.input$min.biomass, wa.input$max.biomass))) +
  ylim((c(1, length(wa.input$biomass)))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  scale_x_continuous(trans = 'log10', breaks = c(0, 1, 5, 10)) +
  geom_line(aes(x = wax.PLB, y = way.PLB), col = 'black', lwd = 1) +
  annotate("text", x=0.069, y=5.5, label="Wakatobi") +
  annotate("text", x=0.07, y=1.5, label = expression(paste(italic("b = "), -1.71))) +
  theme_classic()

# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.wa.b - 0.5, PLB.bMLE.wa.b + 0.5, 0.00001) 
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
  PLB.LLvals[i] = negLL.PLB(bvec[i], x=wa.input$biomass, n=length(wa.input$biomass), xmin=wa.input$min.biomass,
  xmax=wa.input$max.biomass, sumlogx=wa.input$sum.log.biomass)   
}
critVal = PLB.minNegLL.wa.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
wabIn95 <- c(min(bIn95), max(bIn95))

# MLE Lombok biomass ===============================================================================

# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.lo <- mle_b(region="lombok", x=lo.input$biomass, log_x=lo.input$log.biomass, sum_log_x=lo.input$sum.log.biomass,
                 x_min=lo.input$min.biomass, x_max=lo.input$max.biomass)
PLB.bMLE.lo.b <- PLB.return.lo[[1]] 
PLB.minLL.lo.b <- PLB.return.lo[[2]]

# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.lo.b <- PLB.minLL.lo.b$minimum
x <- lo.input$biomass
lox.PLB = seq(min(lo.input$biomass), max(lo.input$biomass), length=1000) # x values to plot PLB. Note
                                                                         # that these encompass the data, and are not based
                                                                         # on the binning (in MEE Figure 6 the line starts as
                                                                         # min(x), not the first bin.
loy.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.lo.b, xmin = min(x.PLB),
    xmax = max(x.PLB))) * length(lo.input$biomass)
spectra.text <- as.character(round(PLB.bMLE.lo.b, 2))
lob_plot <- ggplot() +
  geom_point(aes(x = (sort(lo.input$biomass, decreasing=TRUE)), y = (1:length(lo.input$biomass))), 
             color = "#009E73", size = 2, alpha = 0.3) +
  xlab(expression(paste("Body sizes, ", italic("x"), " (kg)"))) +
  ylab(expression(paste("Number of body sizes", " ">=" ", italic("x")))) +
  xlim((c(lo.input$min.biomass, lo.input$max.biomass))) +
  ylim((c(1, length(lo.input$biomass)))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000)) +
  scale_x_continuous(trans = 'log10', breaks = c(0, 1, 5, 10)) +
  geom_line(aes(x = lox.PLB, y = loy.PLB), col = 'black', lwd = 1) +
  annotate("text", x=0.069, y=15, label="Lombok") +
  annotate("text", x=0.07, y=5, label = expression(paste(italic("b = "), -2.06))) +
  theme_classic()

# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.lo.b - 0.5, PLB.bMLE.lo.b + 0.5, 0.00001) 
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
  PLB.LLvals[i] = negLL.PLB(bvec[i], x=lo.input$biomass, n=length(lo.input$biomass), xmin=lo.input$min.biomass,
  xmax=lo.input$max.biomass, sumlogx=lo.input$sum.log.biomass)   
}
critVal = PLB.minNegLL.lo.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
lobIn95 <- c(min(bIn95), max(bIn95))

# Plot estimated size spectra slopes with confidence intervals =========================================
slopes_ci_df <- data.frame(region = c("Raja Ampat", "Wakatobi", "Lombok"),
                           b = c(PLB.bMLE.ra.b, PLB.bMLE.wa.b, PLB.bMLE.lo.b),
                           b_lo = c(rabIn95[1], wabIn95[1],lobIn95[1]),
                           b_up = c(rabIn95[2], wabIn95[2],lobIn95[2]))
slopes_ci_df$region <- as.factor(slopes_ci_df$region)
slopes_ci_df$region <- factor(slopes_ci_df$region, levels = c("Raja Ampat", "Wakatobi", "Lombok"))
slopes_ci <- ggplot() +
  geom_point(data = slopes_ci_df, aes(x = c(1,1,1), y = b, color = region)) +
  geom_errorbar(data = slopes_ci_df, aes(x = c(1,1,1), ymin = b_lo, ymax = b_up, color = region), width = 0.05) +
  ylab(expression(italic("b"))) +
  xlab("") +
  scale_x_continuous(limits = c(0.9,1.5)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 14),
        legend.position = c(0.5,0.5))

# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggarrange(rab_plot, wab_plot, lob_plot, slopes_ci, ncol=2, nrow=2)

# calculate size spectra slope for each study site =====================================================
site.names <- as.character(unique(fish.df$site_name)) # save site names as a vector
covariates.df$b <- NA
covariates.df$b.weight <- NA
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
for(i in 1:length(site.names)){
  site.i <- site.names[i] # get name of first site
  site.i.df <- fish.df %>% filter(site_name == site.i) # save site data frame
  # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
  # as a starting point for nlm for MLE of b for PLB model.
  tmp.input <- set.params(site.i.df$biomass_kg)
  PLB.return <- mle_b(region="na", x=tmp.input$biomass, log_x=tmp.input$log.biomass, sum_log_x=tmp.input$sum.log.biomass,
    x_min=tmp.input$min.biomass, x_max=tmp.input$max.biomass)
  stderr <- sqrt(abs(diag(solve(PLB.return[[2]]$hessian))))
  stdev <- stderr*(sqrt(length(site.i.df$biomass_kg)))
  inv.var <- 1/(stdev^2)
  
  PLB.bMLE.b <- PLB.return[[1]]
  PLB.minLL.b <- PLB.return[[2]]$minimum
  covariates.df$b[which(covariates.df$site_name == site.i)] <- PLB.bMLE.b
  covariates.df$b.weight[which(covariates.df$site_name == site.i)] <- inv.var
}
 
# slope (b) in relation to human population gravity ====================================================
ggplot() +
  geom_point(data=covariates.df, aes(x=log(Grav_tot), y=b, color=region)) +
  geom_smooth(data=covariates.df, aes(x=log(Grav_tot), y=b), color="black", method="lm") +
  theme_classic() +
  labs(x="Human population gravity", y="Size spectra slope (b)") +
  scale_color_discrete(name=NULL, labels=c("Lombok","Raja Ampat","Wakatobi")) +
  # scale_x_continuous(limits = c(0,900), expand = c(0,0)) +
  theme(legend.position = c(0.9,0.9))

lm.1 <- lm(covariates.df$b ~ log(covariates.df$Grav_tot))
summary(lm.1)
plot(lm.1)
cooksd <- cooks.distance(lm.1)
plot(cooksd, pch="*")
abline(h = 4 * mean(cooksd, na.rm=T), col = "red")
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

min(covariates.df$b)

# remove outlier
tmp.covariates.df <- covariates.df[-c(5,16),] 
ggplot() +
  geom_point(data=tmp.covariates.df, aes(x=log(Grav_tot), y=b, color=region)) +
  geom_smooth(data=tmp.covariates.df, aes(x=log(Grav_tot), y=b), color="black", method="lm") +
  theme_classic() +
  labs(x="log(Human population gravity)", y="Size spectra slope (b)") +
  scale_color_discrete(name=NULL, labels=c("Lombok","Raja Ampat","Wakatobi")) +
  # scale_x_continuous(limits = c(0,900), expand = c(0,0)) +
  theme(legend.position = c(0.9,0.9))

lm.2 <- lm(tmp.covariates.df$b ~ log(tmp.covariates.df$Grav_tot))
summary(lm.2)

# slope (b) in relation to biomass (kg/ha) =============================================================

# calculate the mean kg/ha for each site
tmp <- fish.df %>%
  dplyr::group_by(site_name, transect, observer) %>%
  dplyr::summarize(biomass_250m_sq = sum(biomass_kg)) %>%
  dplyr::group_by(site_name) %>%
  dplyr::summarize(mean_bio_density = mean(biomass_250m_sq)) %>%
  mutate(mean_bio_hectare = mean_bio_density * 40) %>%
  dplyr::select(site_name, mean_bio_hectare)

covariates.df <- covariates.df %>%
    left_join(., tmp, by = "site_name")

ggplot() +
  geom_point(data=covariates.df, aes(x=mean_bio_hectare, y=b, color=region)) +
  geom_smooth(data=covariates.df, aes(x=mean_bio_hectare, y=b), color="black", method="lm") +
  theme_classic() +
  labs(x="Mean biomass (kg/ha)", y="Size spectra slope (b)") +
  scale_color_discrete(name=NULL, labels=c("Lombok","Raja Ampat","Wakatobi")) +
  scale_x_continuous(limits = c(0,3000), expand = c(0,0)) +
  theme(legend.position = c(0.9,0.9))
  
lm.3 <- lm(covariates.df$b ~ covariates.df$mean_bio_hectare)
summary(lm.3)

# GAMs =================================================================================================

# Pairs plot for covariates
panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
pairs(covariates.df[,c(14,9,10,11)], 
      labels = c("Biomass", "Hard coral cover", "Algal cover", "Structural complexity"),
      lower.panel = panel.cor,
      upper.panel = panel.smooth,
      pch = 19)
cor(covariates.df$mean_bio_hectare, covariates.df$hard_coral, method = "pearson")
cor(covariates.df$mean_bio_hectare, covariates.df$algae, method = "pearson")
cor(covariates.df$mean_bio_hectare, covariates.df$mean_complexity, method = "pearson")

covariates.df$region <- as.factor(covariates.df$region)

# Model with all covariates and region as a fixed effect
gam1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) + region,
          data = covariates.df, family = Gamma(link=log), na.action = "na.fail", weights = covariates.df$b.weight)
plot(gam1)
summary(gam1)
gam.check(gam1)
# Test for concurvity
concurvity(gam1, full = TRUE)
concurvity(gam1, full = FALSE)
dd <- dredge(gam1, rank = "AICc")

# Remove structural complexity
gam2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + region,
          data = covariates.df, family = Gamma(link=log), na.action = "na.fail", weights = covariates.df$b.weight)
summary(gam2)
concurvity(gam2, full = FALSE)
dredge(gam2, rank = "AICc")

# Remove hard coral
gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(algae, k=3) + s(mean_complexity, k=3) + region,
          data = covariates.df, family = Gamma(link=log), na.action = "na.fail", weights = covariates.df$b.weight)
summary(gam3)
concurvity(gam3, full = FALSE)
dredge(gam3, rank = "AICc")


# Plot sum of AICc weights
dd.gam2 <- dredge(gam2)
dd.gam2.df <- data.frame(Biomass = dd.gam2$`s(mean_bio_hectare, k = 3)`,
                         Algae = dd.gam2$`s(algae, k = 3)`,
                         Hardcoral = dd.gam2$`s(hard_coral, k = 3)`,
                         weight = dd.gam2$weight)
dd.gam2.df <- dd.gam2.df %>% 
  gather(., "covariate", "included", -weight) %>%
  na.omit() %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarize(sum_weight = sum(weight)) %>%
  mutate(model = "Hard coral cover")

dd.gam3 <- dredge(gam3)
dd.gam3.df <- data.frame(Biomass = dd.gam3$`s(mean_bio_hectare, k = 3)`,
                         Algae = dd.gam3$`s(algae, k = 3)`,
                         Complexity = dd.gam3$`s(mean_complexity, k = 3)`,
                         weight = dd.gam3$weight)
dd.gam3.df <- dd.gam3.df %>% 
  gather(., "covariate", "included", -weight) %>%
  na.omit() %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarize(sum_weight = sum(weight)) %>%
  mutate(model = "Structural complexity")

dd.gam.df <- rbind(dd.gam2.df, dd.gam3.df)
dd.gam.df$covariate <- c("Algal cover", "Biomass", "Hard coral cover", "Algal cover", "Biomass", "Structural complexity")

# Plot
plot_base <- ggplot(data=dd.gam.df, aes(x = covariate, y = sum_weight, fill = model)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.5, lty = "dashed") +
  scale_fill_manual(values = c("#F1CE75", "#B2D1E8")) +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  facet_grid(rows = vars(model), scales = "free_y", space = "free_y", switch = "y") + #
  theme(strip.background = element_rect(fill="white"),
        strip.placement = "outside") +
  guides(fill = FALSE) +
  xlab("") +
  ylab("Sum of AICc weights")
  
# Plot model averages
gam2.modavg <- model.avg(dd.gam2)
gam2.modavg.coef <- coef(gam2.modavg)
gam2.modavg.confint <- confint(gam2.modavg)
plot(gam2.modavg.coef)

# # Load data
# fish.df <- read.csv("indonesia_uvc_fish.csv")
# benthos_df <- read.csv("indonesia_benthos.csv")
# gravity.df <- read.csv("../gravity/gravity.csv")
# strcomp_df <- read.csv("structural_complexity.csv")
# benthos_key <- read_xlsx("benthos_key.xlsx")

# # Fix site names in gravity dataframe
# gravity.df$site_name <- as.character(gravity.df$site_name)
# gravity.df$site_name[which(gravity.df$site_name == "Yanbuba 1")] <- "Yenbuba 1"
# gravity.df$site_name[which(gravity.df$site_name == "Saonek")] <- "Saonek Besar"
# gravity.df$site_name[which(gravity.df$site_name == "Yenbeser")] <- "Yenbesar"
# gravity.df$site_name[which(gravity.df$site_name == "Yansawai")] <- "Yensawai"
# gravity.df$site_name[which(gravity.df$site_name == "Mambetron site")] <- "Mambetron"
# gravity.df$site_name[which(gravity.df$site_name == "Gof Island")] <- "Gof"
# gravity.df$site_name[which(gravity.df$site_name == "Paniki Besar Island")] <- "Paniki Besar"

# # data cleaning for fish dataframe =====================================================================
# 
# # remove chaetodontids and pomacentrids larger than 20cm because the max sizes are 20cm
# fish.df <- fish.df %>%
#     filter(!(family == "chaetodontidae" & size_cm > 20)) %>%
#     filter(!(family == "pomacentridae" & size_cm > 20)) %>%
#     filter(!(site_name == "Sombano"))
# 
# # Expand dataframe such that abundance is 1 for each row
# fish.df <- fish.df %>%
#      dplyr::select(site_name, region, transect, genus_species, genus, species, family, size_cm, size_5cm_bin, abundance, 
#             functional_group, a, b, observer)
# df <- setDT(expandRows(fish.df,"abundance"))[, abundance := sprintf("1")][]
# fish.df <- as.data.frame(df) 
# 
# # Convert abundance to biomass
# fish.df <- fish.df %>%
#      mutate(abundance = as.numeric(abundance)) %>%
#      mutate(biomass_g = ((a * (size_cm ^ b)) * as.numeric(abundance))) %>%
#      mutate(biomass_kg = round((biomass_g/1000), digits = 2)) %>%
#      filter(size_cm > 10 & size_cm < 67) %>% # underwater visual census is best for fishes in this size range (Kulbicki 1998)
#      filter(family != "caesionidae") %>% # caesionids found in large school and observers have difficultly estimating abundance (Samoilys and Carlos 2000)
#      filter(biomass_kg != 0) %>%
#      filter(biomass_kg > 0.03) # fish smaller than this size are likely inadequately sampled
      
# # data cleaning for benthic dataframe ==================================================================
# 
# benthos_df <- benthos_df %>%
#     mutate(hard_soft = gsub("abiotik","abiotic", .$hard_soft)) %>%
#     mutate(hard_soft = gsub("makro algae", "macro algae", .$hard_soft)) %>%
#     filter(site_name != "Sombano")

# # trophic group assignement ============================================================================
# 
# fish.df <- fish.df %>% 
#   mutate(tp = ifelse(functional_group == "Obligate and facultative coral feeder" | functional_group == "benthic invertivore" |
#                      functional_group == "omnivore" | functional_group == "carnivore" | functional_group == "Sessile invertebrate feeder" |
#                      functional_group == "micro-invertivore" | functional_group == "macro-invertivore" |
#                      functional_group == "pisci-invertivore" | functional_group == "piscivore" |
#                      functional_group == "corallivore", "Carnivore", 
#                      ifelse(functional_group == "browser" | functional_group == "detritivore" |
#                             functional_group == "grazer" | functional_group == "scraper" | 
#                             functional_group == "large excavator" | functional_group == "small excavator" |
#                             functional_group == "grazer/detritivore" | functional_group == "excavator/scraper", "Herbivore", NA)))

# # calculate percent cover algae and hard coral for each site ===========================================
# 
# benthos_df1 <- benthos_df %>%
#     mutate(hard_coral = ifelse(.$hard_soft == "hard coral", 1, 0)) %>%
#     mutate(algae = ifelse(.$hard_soft == "macro algae" | .$hard_soft == "algae", 1, 0)) %>%
#     group_by(site_name, transect, observer) %>%
#     summarise_at(c("hard_coral","algae"), sum, na.rm = TRUE) %>%
#     group_by(site_name) %>%
#     summarise_at(c("hard_coral","algae"), mean, na.rm = TRUE) # hard coral and algae values are percent coverage
# 
# benthos_df2 <- benthos_df %>%
#     left_join(., benthos_key, by = "lifeform") %>%
#     select(site_name, transect, transect_point, observer, lifeform = new_id)
# benthos_df2_tmp <- as.data.frame(to.dummy(benthos_df2$lifeform, "lf"))
# benthos_df2 <- cbind(benthos_df2, benthos_df2_tmp)
# names(benthos_df2) <- c("site_name","transect","transect_point","observer","lifeform","DC",
#                         "HC_B","HC_D","HC_E","HC_F","HC_L","HC_M","HC_S",
#                         "HC_T","MA","NA","SC")              
# benthos_df2 <- benthos_df2 %>%
#     group_by(site_name, transect, observer) %>%
#     summarise_at(c("DC","HC_B","HC_D","HC_E","HC_F","HC_L","HC_M","HC_S","HC_T","MA","SC"), sum, na.rm = TRUE) %>%
#     group_by(site_name) %>%
#     summarise_at(c("DC","HC_B","HC_D","HC_E","HC_F","HC_L","HC_M","HC_S","HC_T","MA","SC"), mean, na.rm = TRUE)
# 
# gravity.df <- gravity.df %>%
#     left_join(., benthos_df1, by = "site_name")

# # calculate mean structural complexity for each site ===================================================
# 
# strcomp_df1 <- strcomp_df %>%
#     dplyr::group_by(site_name) %>%
#     dplyr::summarize(mean_complexity = mean(complexity, na.rm = TRUE))
# 
# gravity.df <- gravity.df %>%
#     left_join(., strcomp_df1, by = "site_name")

  
  


# Plot MLE for each site ====================================================================================

# site.names <- as.character(unique(fish.df$site_name)) # save site names as a vector
# gravity.df$b <- NA
# for(i in 1:length(site.names)){
#     site.i <- site.names[i] # set the site name
#     site.i.df <- fish.df %>% filter(site_name == site.i) # save site data frame
#     # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
#     # as a starting point for nlm for MLE of b for PLB model.
#     tmp.input <- set.params(site.i.df)
#     PLB.return <- mle_b(region=NA, x=tmp.input$biomass, log_x=tmp.input$log.biomass, sum_log_x=tmp.input$sum.log.biomass,
#                             x_min=tmp.input$min.biomass, x_max=tmp.input$max.biomass)
#     PLB.bMLE <- PLB.return[[1]] 
#     PLB.minLL <- PLB.return[[2]]
#         
#     # plot and find 95% confidence intervals for MLE method.
#     PLB.minNegLL <- PLB.minLL$minimum
#     x <- tmp.input$biomass # for plotting purposes
#     plot(sort(tmp.input$biomass, decreasing=TRUE), 1:length(tmp.input$biomass), log="xy",
#       xlab=expression(paste("Body sizes, ", italic(x), " (kg)")),
#       ylab=expression( paste("Number of Body ", sizes >= x)), mgp=mgpVals,
#       xlim = c(tmp.input$min.biomass, tmp.input$max.biomass), ylim = c(1, length(tmp.input$biomass)), axes=FALSE)
#       logTicks(xLim, yLim, xLabelBig = c(0, 1, 10, 100)) # Tick marks.
#       x.PLB = seq(min(tmp.input$biomass), max(tmp.input$biomass), length=1000) # x values to plot PLB. Note
#                                                                                # that these encompass the data, and are not based
#                                                                                # on the binning (in MEE Figure 6 the line starts as
#                                                                                # min(x), not the first bin.
#       y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE, xmin = min(x.PLB),
#         xmax = max(x.PLB))) * length(tmp.input$biomass)
#       lines(x.PLB, y.PLB, col="red", lwd=2)
#       text(x=0.065, y=5, labels=paste("b = ", as.character(round(PLB.bMLE, 2)), sep=""),
#         cex=1.5, pos=1, col="black")
#       title(main=site.i)
# }


# # GAMs =================================================================================================
# 
# # PCA
# grav_tmp <- gravity.df %>% select(site_name, region, mean_complexity)
# benthos_df3 <- benthos_df2 %>% left_join(., grav_tmp, by = "site_name")
# benthos.pca <- prcomp(gravity.df[,c(8:10)], center = TRUE, scale. = TRUE)
# summary(benthos.pca)
# ggbiplot(benthos.pca, groups = gravity.df$region)
# 
# 
# benthos.pca <- prcomp(benthos_df3[,c(2:12,14)], center = TRUE, scale. = TRUE)
# summary(benthos.pca)
# ggbiplot(benthos.pca, groups = benthos_df3$region)
# 
# 
# 
# # Models with Gamma family and log link
# gam1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3),
#           data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
# summary(gam1)
# concurvity(gam1, full=FALSE)
# dd <- dredge(gam1, rank = "AICc")

# concurvity shows potential covariance between multiple variables and structural complexity

gam1.1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) +
              te(mean_bio_hectare,mean_complexity, k=3) + te(hard_coral,mean_complexity, k=3),
              data = gravity.df, family = Gamma(link = log), na.action = "na.fail", weights = gravity.df$b.weight)
summary(gam1.1)
dd.1 <- dredge(gam1.1, rank = "AICc")

gam1.2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3),
              data = gravity.df, family = Gamma(link = log), na.action = "na.fail", weights = gravity.df$b.weight)
summary(gam1.2)

gam1.3 <- gam(abs(b) ~ s(mean_bio_hectare, k = 3),
              data = gravity.df, family = Gamma(link = log), na.action = "na.fail", weights = gravity.df$b.weight)
par(mfrow=c(2,2))
gam.check(gam1.3)
summary(gam1.3)
AICc(gam1.2, gam1.3)
anova(gam1, gam1.2, gam1.3)

# remove structural complexity
gam2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
summary(gam2)

# remove biomass
gam3 <- gam(abs(b) ~ s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
AICc(gam1, gam2, gam3) 

# removing structural complexity was a better model fit

# try interactions with structural complexity
gam4 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) +
              te(mean_bio_hectare, mean_complexity),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
gam5 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) +
              te(hard_coral, mean_complexity),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
aic <- AICc(gam1, gam2, gam3, gam4, gam5)
aic[order(aic$AICc),]

# removing structural complexity was a better model fit

# try interactions with hard coral
gam6 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) +
              te(mean_bio_hectare, hard_coral) + te(mean_bio_hectare, algae),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
summary(gam6)
dd1 <- dredge(gam6)
aic <- AICc(gam1, gam2, gam6)
aic[order(aic$AICc),]

# try interactions with algae
gam7 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) +
              te(mean_bio_hectare, algae) + te(algae, hard_coral),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
aic <- AICc(gam1, gam2, gam6, gam7)
aic[order(aic$AICc),]

# interactions with algae does not help

# try removing algae
gam8 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
summary(gam8)
aic <- AICc(gam1, gam2, gam6, gam8)
aic[order(aic$AICc),]
 
# removing algae improved model fit

# try biomass and hard coral cover with interaction
gam9 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + te(mean_bio_hectare, hard_coral),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail", weights = gravity.df$b.weight)
summary(gam9)
aic <- AICc(gam1, gam2, gam6, gam9)
aic[order(aic$AICc),]

# try biomass and interaction only with hard_coral
# gam10 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + te(mean_bio_hectare, hard_coral),
          # data = gravity.df, family = Gamma(link=log), na.action = "na.fail")
gam10 <- gam(abs(b) ~ s(mean_bio_hectare, k = 3), data = gravity.df, weights = gravity.df$b.weight,
             family = Gamma(link = log), na.action = "na.fail")
summary(gam10)
aic <- AICc(gam1, gam2, gam6, gam9, gam10)
aic[order(aic$AICc),]

# biomass as response and structural complexity as a predictor
gamb1 <- gam(mean_bio_hectare ~ s(hard_coral, k = 3),
             data = gravity.df, family = Gamma(link=log), na.action = "na.fail") 
summary(gamb1)
ddb1 <- dredge(gamb1)

gamb2 <- gam(hard_coral ~ s(mean_complexity, k = 3),
             data = gravity.df, family = Gamma(link=log), na.action = "na.fail") 
summary(gamb2)
newdb <- data.frame(mean_complexity = seq(min(gravity.df$mean_complexity), max(gravity.df$mean_complexity), length.out = 2000))
predb <- predict.gam(gamb2, newdb, se.fit = TRUE)
upr <- predb$fit + (2 * predb$se.fit)
lwr <- predb$fit - (2 * predb$se.fit)
newdb$b_pred <- exp(predb$fit)
newdb$b_pred_upr <- exp(upr)
newdb$b_pred_lwr <- exp(lwr)
tmp_grav <- gravity.df
gravity.df <- tmp_grav
gravity.df$region <- factor(gravity.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newdb) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newdb) +
  geom_point(aes(x = mean_complexity, y = hard_coral, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Hard coral cover") +
  xlab("Structural Complexity")

gamb3 <- gam(mean_bio_hectare ~ s(mean_complexity, k = 3), 
             data = gravity.df, family = Gamma(link=log), na.action = "na.fail")
summary(gamb3)
newdb <- data.frame(mean_complexity = seq(min(gravity.df$mean_complexity), max(gravity.df$mean_complexity), length.out = 2000))
predb <- predict.gam(gamb3, newdb, se.fit = TRUE)
upr <- predb$fit + (2 * predb$se.fit)
lwr <- predb$fit - (2 * predb$se.fit)
newdb$b_pred <- exp(predb$fit)
newdb$b_pred_upr <- exp(upr)
newdb$b_pred_lwr <- exp(lwr)
tmp_grav <- gravity.df
gravity.df <- tmp_grav
gravity.df$region <- factor(gravity.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newdb) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newdb) +
  geom_point(aes(x = mean_complexity, y = mean_bio_hectare, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Hard coral cover") +
  xlab("Structural Complexity")

AICc(gamb1, gamb2)


newdb <- data.frame(mean_complexity = seq(min(gravity.df$mean_complexity), max(gravity.df$mean_complexity), length.out = 2000))
predb <- predict.gam(gamb2, newdb, se.fit = TRUE)
upr <- predb$fit + (2 * predb$se.fit)
lwr <- predb$fit - (2 * predb$se.fit)
newdb$b_pred <- exp(predb$fit)
newdb$b_pred_upr <- exp(upr)
newdb$b_pred_lwr <- exp(lwr)
tmp_grav <- gravity.df
gravity.df <- tmp_grav
gravity.df$region <- factor(gravity.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newdb) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newdb) +
  geom_point(aes(x = mean_complexity, y = mean_bio_hectare, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Fishable Biomass") +
  xlab("Structural Complexity")


# Plot predicted values ================================================================================

newd <- data.frame(mean_bio_hectare = seq(min(gravity.df$mean_bio_hectare), max(gravity.df$mean_bio_hectare), length.out = 2000))
pred <- predict.gam(gam1.3, newd, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd$b_pred <- exp(pred$fit) * -1
newd$b_pred_upr <- exp(upr) * -1
newd$b_pred_lwr <- exp(lwr) * -1
tmp_grav <- gravity.df
gravity.df <- tmp_grav
gravity.df$region <- factor(gravity.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = mean_bio_hectare, y = b_pred), size = 1.5, data = newd) +
  geom_ribbon(aes(x = mean_bio_hectare, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd) +
  geom_point(aes(x = mean_bio_hectare, y = b, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.2)) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab(expression(paste("Size spectra (", italic(b), ")"))) +
  xlab("Biomass (kg/ha)")

# Check slope in relation to complexity and hard coral cover --------------

# Check residuals to make sure the model meets assumptions
gam3 <- gam(abs(b) ~ s(hard_coral, k=3), data = gravity.df, family = Gamma(link=log))
summary(gam3)
par(mfrow = c(2,2))
gam.check(gam3)
newd <- data.frame(hard_coral = seq(min(gravity.df$hard_coral), max(gravity.df$hard_coral), length.out = 2000))
pred <- predict.gam(gam3, newd, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd$b_pred <- exp(pred$fit) * -1
newd$b_pred_upr <- exp(upr) * -1
newd$b_pred_lwr <- exp(lwr) * -1
tmp_grav <- gravity.df
gravity.df <- tmp_grav
gravity.df$region <- factor(gravity.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = hard_coral, y = b_pred), size = 1.5, data = newd) +
  geom_ribbon(aes(x = hard_coral, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd) +
  geom_point(aes(x = hard_coral, y = b, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Size spectrum slope (b)") +
  xlab("hard_coral")



# Density of small, medium and large fishes at each site ==================================================

fish.sizes.df <- fish.df %>%
  mutate(size_cat = ifelse(biomass_kg < 0.2, "small",
                    ifelse(biomass_kg >= 0.2 & biomass_kg < 1.2, "medium", "large"))) %>%
  group_by(region, site_name, transect, observer, size_cat) %>%
  summarize(biomass_kg = sum(biomass_kg)) %>%
  group_by(region, size_cat) %>%
  summarize(bio = mean(biomass_kg),
            stderr = std.error(biomass_kg)) %>%
  mutate(bio_ha = bio * 40) %>%
  mutate(stderr_ha = stderr * 40) %>%
  ungroup()

fish.sizes.df <- fish.sizes.df %>% mutate(region = fct_relevel(region, "raja_ampat", "wakatobi", "lombok"))

ggplot(data = fish.sizes.df, aes(x = size_cat, y = bio_ha, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = bio_ha - stderr_ha, ymax = bio_ha + stderr_ha), stat = "identity",
                width = 0.2, position = position_dodge(0.9)) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.title = element_blank()) +
  scale_fill_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  scale_x_discrete(labels = c("Large", "Medium", "Small")) +
  labs(x = "Size category", y = "Biomass (kg/ha)")

# Species abundance for carnivores and herbivores ======================================================
species.carn <- fish.df %>%
  filter(tp == "Carnivore") %>%
  group_by(genus_species) %>%
  summarise(total = sum(abundance)) %>%
  filter(total >= 10)

ggplot() +
  geom_bar(data = species.carn, aes(x = reorder(genus_species, total), y = total), stat = "identity") +
  coord_flip()

species.herb <- fish.df %>%
  filter(tp == "Herbivore") %>%
  group_by(genus_species) %>%
  summarise(total = sum(abundance)) %>%
  filter(total >= 10)

ggplot() +
  geom_bar(data = species.herb, aes(x = reorder(genus_species, total), y = total), stat = "identity") +
  coord_flip()

# Trophic position analysis for entire dataset =========================================================

# Set MLE parameters

# Carnivore biomass
carn.input <- set.params(carn_df)
# Herbivore biomass
herb.input <- set.params(herb_df)

# MLE CARNIVORE
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.carn <- mle_b(region=NA, x=carn.input$biomass, log_x=carn.input$log.biomass, sum_log_x=carn.input$sum.log.biomass,
                 x_min=carn.input$min.biomass, x_max=carn.input$max.biomass)
PLB.bMLE.carn.b <- PLB.return.carn[[1]] 
PLB.minLL.carn.b <- PLB.return.carn[[2]]

trophic.plots(PLB.return.carn, PLB.bMLE.carn.b, PLB.minLL.carn.b, carn.input, mgpVals, troph_id = "Carnivores", panel = "a")
# title(main = "Carnivore")

test <- slope.conf.int(PLB.bMLE.carn.b, PLB.minLL.carn.b$minimum, carn.input)

# MLE HERBIVORE
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.herb <- mle_b(region=NA, x=herb.input$biomass, log_x=herb.input$log.biomass, sum_log_x=herb.input$sum.log.biomass,
                 x_min=herb.input$min.biomass, x_max=herb.input$max.biomass)
PLB.bMLE.herb.b <- PLB.return.herb[[1]] 
PLB.minLL.herb.b <- PLB.return.herb[[2]]

trophic.plots(PLB.return.herb, PLB.bMLE.herb.b, PLB.minLL.herb.b, herb.input, mgpVals, troph_id = "Herbivores", panel = "b")
# title(main="Herbivore")

test <- slope.conf.int(PLB.bMLE.herb.b, PLB.minLL.herb.b$minimum, herb.input)

# Trophic position analysis for each region =========================================================

# Create separate dataframes
ra.carn <- carn_df %>% filter(region == "raja_ampat")
ra.herb <- herb_df %>% filter(region == "raja_ampat")
wa.carn <- carn_df %>% filter(region == "wakatobi")
wa.herb <- herb_df %>% filter(region == "wakatobi")
lo.carn <- carn_df %>% filter(region == "lombok")
lo.herb <- herb_df %>% filter(region == "lombok")
        
# Set MLE parameters

# Raja Ampat 
ra.carn.input <- set.params(ra.carn)
ra.herb.input <- set.params(ra.herb)

# Wakatobi
wa.carn.input <- set.params(wa.carn)
wa.herb.input <- set.params(wa.herb)

# Lombok
lo.carn.input <- set.params(lo.carn)
lo.herb.input <- set.params(lo.herb)

# MLE CARNIVORE - Raja Ampat
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.ra.carn <- mle_b(region=NA, x=ra.carn.input$biomass, log_x=ra.carn.input$log.biomass, sum_log_x=ra.carn.input$sum.log.biomass,
                 x_min=ra.carn.input$min.biomass, x_max=ra.carn.input$max.biomass)
PLB.bMLE.ra.carn.b <- PLB.return.ra.carn[[1]] 
PLB.minLL.ra.carn.b <- PLB.return.ra.carn[[2]]

trophic.plots(PLB.return.ra.carn, PLB.bMLE.ra.carn.b, PLB.minLL.ra.carn.b, ra.carn.input, mgpVals)
title(main="Raja Ampat - Carnivore")

# MLE HERBIVORE - Raja Ampat
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.ra.herb <- mle_b(region=NA, x=ra.herb.input$biomass, log_x=ra.herb.input$log.biomass, sum_log_x=ra.herb.input$sum.log.biomass,
                 x_min=ra.herb.input$min.biomass, x_max=ra.herb.input$max.biomass)
PLB.bMLE.ra.herb.b <- PLB.return.ra.herb[[1]] 
PLB.minLL.ra.herb.b <- PLB.return.ra.herb[[2]]

trophic.plots(PLB.return.ra.herb, PLB.bMLE.ra.herb.b, PLB.minLL.ra.herb.b, ra.herb.input, mgpVals)
title(main="Raja Ampat - Herbivore")

# MLE CARNIVORE - Wakatobi
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.wa.carn <- mle_b(region=NA, x=wa.carn.input$biomass, log_x=wa.carn.input$log.biomass, sum_log_x=wa.carn.input$sum.log.biomass,
                 x_min=wa.carn.input$min.biomass, x_max=wa.carn.input$max.biomass)
PLB.bMLE.wa.carn.b <- PLB.return.wa.carn[[1]] 
PLB.minLL.wa.carn.b <- PLB.return.wa.carn[[2]]

trophic.plots(PLB.return.wa.carn, PLB.bMLE.wa.carn.b, PLB.minLL.wa.carn.b, wa.carn.input, mgpVals)
title(main="Wakatobi - Carnivore")

# MLE HERBIVORE - Wakatobi
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.wa.herb <- mle_b(region=NA, x=wa.herb.input$biomass, log_x=wa.herb.input$log.biomass, sum_log_x=wa.herb.input$sum.log.biomass,
                 x_min=wa.herb.input$min.biomass, x_max=wa.herb.input$max.biomass)
PLB.bMLE.wa.herb.b <- PLB.return.wa.herb[[1]] 
PLB.minLL.wa.herb.b <- PLB.return.wa.herb[[2]]

trophic.plots(PLB.return.wa.herb, PLB.bMLE.wa.herb.b, PLB.minLL.wa.herb.b, wa.herb.input, mgpVals)
title(main="Wakatobi - Herbivore")

# MLE CARNIVORE - Lombok
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.lo.carn <- mle_b(region=NA, x=lo.carn.input$biomass, log_x=lo.carn.input$log.biomass, sum_log_x=lo.carn.input$sum.log.biomass,
                 x_min=lo.carn.input$min.biomass, x_max=lo.carn.input$max.biomass)
PLB.bMLE.lo.carn.b <- PLB.return.lo.carn[[1]] 
PLB.minLL.lo.carn.b <- PLB.return.lo.carn[[2]]

trophic.plots(PLB.return.lo.carn, PLB.bMLE.lo.carn.b, PLB.minLL.lo.carn.b, lo.carn.input, mgpVals)
title(main="Lombok - Carnivore")

# MLE HERBIVORE - Lombok
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.lo.herb <- mle_b(region=NA, x=lo.herb.input$biomass, log_x=lo.herb.input$log.biomass, sum_log_x=lo.herb.input$sum.log.biomass,
                 x_min=lo.herb.input$min.biomass, x_max=lo.herb.input$max.biomass)
PLB.bMLE.lo.herb.b <- PLB.return.lo.herb[[1]] 
PLB.minLL.lo.herb.b <- PLB.return.lo.herb[[2]]

trophic.plots(PLB.return.lo.herb, PLB.bMLE.lo.herb.b, PLB.minLL.lo.herb.b, lo.herb.input, mgpVals)
title(main="Lombok - Herbivore")

# slope (b) for each trophic group in relation to biomass (kg/ha) ======================================

sites.trophic <- fish.df %>%
  dplyr::select(site_name, tp) %>%
  distinct() %>%
  filter(!(is.na(tp)))
sites.trophic$b <- NA
sites.trophic$b.weight <- NA
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
for(i in 1:length(sites.trophic$site_name)){
  site.i <- sites.trophic$site_name[i]
  tp.i <- sites.trophic$tp[i]
  site.tp.df <- fish.df %>% filter(site_name == site.i) %>% filter(tp == tp.i)
  # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
  # as a starting point for nlm for MLE of b for PLB model.
  tmp.input <- set.params(site.tp.df)
  PLB.return <- mle_b(region="na", x=tmp.input$biomass, log_x=tmp.input$log.biomass, sum_log_x=tmp.input$sum.log.biomass,
                x_min=tmp.input$min.biomass, x_max=tmp.input$max.biomass)
  stderr <- sqrt(abs(diag(solve(PLB.return[[2]]$hessian))))
  stdev <- stderr*(sqrt(length(site.tp.df$biomass_kg)))
  inv.var <- 1/(stdev^2)
  PLB.bMLE.b <- PLB.return[[1]]
  PLB.minLL.b <- PLB.return[[2]]
  sites.trophic$b[i] <- PLB.bMLE.b
  sites.trophic$b.weight[i] <- inv.var
}

# calculate the mean kg/ha for each site
tmp <- fish.df %>%
  dplyr::group_by(site_name, transect, observer, tp) %>%
  dplyr::summarize(biomass_250m_sq = sum(biomass_kg)) %>%
  dplyr::group_by(site_name, tp) %>%
  dplyr::summarize(mean_bio_density = mean(biomass_250m_sq)) %>%
  filter(!is.na(tp)) %>%
  mutate(mean_bio_hectare = mean_bio_density * 40) %>%
  dplyr::select(site_name, tp, mean_bio_hectare)

sites.trophic <- sites.trophic %>%
  left_join(., tmp, by = c("site_name" = "site_name", "tp" = "tp")) %>%
  left_join(., gravity.df, by = "site_name") %>%
  dplyr::select(site_name, tp, b = b.x, b.weight = b.weight.x, mean_bio_hectare = mean_bio_hectare.x, region, Grav_tot, hard_coral, 
         algae, mean_complexity, total_bio = mean_bio_hectare.y)

ggplot() +
  geom_point(data=sites.trophic, aes(x=total_bio, y=b, shape=region, color=tp), size = 2) +
  theme_classic() +
  labs(x="Biomass density (kg/ha)", y="Size spectra slope (b)") +
  scale_shape_discrete(name = NULL, labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  scale_color_discrete(name = NULL) 
  scale_x_continuous(limits = c(0,2500), expand = c(0,0))

# Abundance of small, medium and large herbs  and carns fishes at each site =================================

herb.sizes.df <- fish.df %>%
  filter(tp == "Herbivore") %>%
  mutate(size_cat = ifelse(biomass_kg < 0.2, "small",
                    ifelse(biomass_kg >= 0.2 & biomass_kg < 1.2, "medium", "large"))) %>%
  group_by(size_cat) %>%
  summarize(abundance = sum(abundance))
  # arrange(biomass_kg)

ggplot() +
  geom_bar(data = herb.sizes.df, aes(x = size_cat, y = abundance), stat = "identity") +
  # scale_y_continuous(lim = c(0, 1100)) +
  labs(x = "Body size (kg)", y = "Abundance", title = "Herbivores")



carn.sizes.df <- fish.df %>%
  filter(tp == "Carnivore") %>%
  mutate(size_cat = ifelse(biomass_kg < 0.2, "small",
                    ifelse(biomass_kg >= 0.2 & biomass_kg < 1.2, "medium", "large"))) %>%
  group_by(size_cat) %>%
  summarize(abundance = sum(abundance))
  # arrange(biomass_kg)

ggplot() +
  geom_bar(data = carn.sizes.df, aes(x = biomass_kg, y = abundance), stat = "identity") +
  # scale_y_continuous(lim = c(0, 1100)) +
  labs(x = "Body size (kg)", y = "Abundance", title = "Carnivores")


# Trophic group GAMs ===================================================================================

carn.gam.df <- sites.trophic %>% filter(tp == "Carnivore")
herb.gam.df <- sites.trophic %>% filter(tp == "Herbivore")
# merge to get the other functional groups biomass
carn.gam.df <- carn.gam.df %>%
  left_join(., herb.gam.df, by = "site_name") %>%
  dplyr::select(site_name, tp = tp.x, b = b.x, b.weight = b.weight.x, mean_bio_hectare = mean_bio_hectare.x, region = region.x,
         Grav_tot = Grav_tot.x, hard_coral = hard_coral.x, algae = algae.x, mean_complexity = mean_complexity.x,
         herb_bio = mean_bio_hectare.y, total_bio = total_bio.x) %>%
  mutate(carn.herb = mean_bio_hectare/herb_bio)

herb.gam.df <- herb.gam.df %>%
  left_join(., carn.gam.df, by = "site_name") %>%
  dplyr::select(site_name, tp = tp.x, b = b.x, b.weight = b.weight.x, mean_bio_hectare = mean_bio_hectare.x, region = region.x,
         Grav_tot = Grav_tot.x, hard_coral = hard_coral.x, algae = algae.x, mean_complexity = mean_complexity.x,
         carn_bio = mean_bio_hectare.y, total_bio = total_bio.x) %>%
  mutate(carn.herb = carn_bio/mean_bio_hectare)

# Carnivore global model
carn.gam1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) +
                   te(mean_bio_hectare, mean_complexity, k = 3) + te(hard_coral, mean_complexity, k = 3), 
                 data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = carn.gam.df$b.weight)
summary(carn.gam1)
plot(influence.gam(carn.gam1))
carn.dd <- dredge(carn.gam1, rank = "AICc")
concurvity(carn.gam1, full = FALSE)

carn.gam1.1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3, fx=TRUE) + s(mean_complexity, k=3, fx=TRUE) + te(mean_bio_hectare, mean_complexity, k = 3), 
                 data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = carn.gam.df$b.weight)
hist(influence.gam(carn.gam1.1))
summary(carn.gam1.1)
gam.check(carn.gam1.1)

carn.gam1.2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(mean_complexity, k=3), 
                 data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = carn.gam.df$b.weight)
summary(carn.gam1.2)


# Carnivore best fit model according to AICc
carn.gam2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3), weights = carn.gam.df$b.weight,
                 data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail")
summary(carn.gam2)
gam.check(carn.gam2)
concurvity(carn.gam2, full = FALSE)

carn.gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + te(mean_bio_hectare, hard_coral), 
                 data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = carn.gam.df$b.weight)
summary(carn.gam3)
AICc(carn.gam2, carn.gam3)


carn.gam.t1 <- gam(abs(b) ~ s(mean_complexity, k=3), data = carn.gam.df, Gamma(link = log), weights = carn.gam.df$b.weight)
summary(carn.gam.t1)
plot(carn.gam.t1)

# test covariance between biomass, coral, and complexity
carn.gamb1 <- gam(mean_bio_hectare ~ s(mean_complexity, k = 3), data = carn.gam.df, family = Gamma(link = log))
summary(carn.gamb1)

# Carn SI model
carn.gam.t1 <- gam(mean_bio_hectare ~ s(hard_coral, k = 3), data = carn.gam.df, Gamma(link = log))
summary(carn.gam.t1)
gam.check(carn.gam.t1)
plot(carn.gam.t1)

# Herbivore global model
herb.gam1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) + 
                   te(hard_coral, mean_complexity, k = 3), 
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = herb.gam.df$b.weight)
summary(herb.gam1)
concurvity(herb.gam1, full = FALSE)
herb.dd <- dredge(herb.gam1, rank = "AICc")

herb.gam1.1 <- gam(abs(b) ~ s(mean_bio_hectare, k = 3) + s(mean_complexity, k=3),
                   data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = herb.gam.df$b.weight)
gam.check(herb.gam1.1)
summary(herb.gam1.1)

herb.gamb1 <- gam(mean_bio_hectare ~ s(mean_complexity, k=3),
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = herb.gam.df$b.weight)
summary(herb.gamb1)
herb.gamb2 <- gam(abs(b) ~ s(mean_complexity, k=3),
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = herb.gam.df$b.weight)
summary(herb.gamb2)
herb.gamb3 <- gam(hard_coral ~ s(mean_complexity, k=3), data = herb.gam.df, 
                  family = Gamma(link = log))
summary(herb.gamb3)

# Herbivore best fit model accoring to AICc
herb.gam2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(algae, k = 3) + s(hard_coral, k=3), 
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail")
summary(herb.gam2)
gam.check(herb.gam2)
concurvity(herb.gam2, full = FALSE)

herb.gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(algae, k = 3) + s(mean_complexity, k=3), 
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail")
summary(herb.gam3)
gam.check(herb.gam3)
concurvity(herb.gam3, full = FALSE)
dredge(herb.gam3)
aic <- AICc(herb.gam2, herb.gam3)
aic[order(aic$AICc),]

herb.gam4 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(mean_complexity, k=3), data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail")
summary(herb.gam4)
concurvity(herb.gam4, full = FALSE)
aic <- AICc(herb.gam2, herb.gam3, herb.gam4)
aic[order(aic$AICc),]

herb.gam5 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(mean_complexity, k=3) + te(mean_complexity, hard_coral), 
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail")
summary(herb.gam5)
dredge(herb.gam5)
concurvity(herb.gam4, full = FALSE)
aic <- AICc(herb.gam2, herb.gam3, herb.gam4, herb.gam5)
aic[order(aic$AICc),]


# covariance
plot(hard_coral ~ mean_complexity, data = herb.gam.df)
herb.gamb1 <- gam(hard_coral ~ s(mean_complexity, k = 3), data = herb.gam.df, family = gaussian(link = log))
summary(herb.gamb1) # hard coral cover and complexity covary such that hard coral cover increases with structural complexity

herb.gam.t1 <- gam(mean_bio_hectare ~ s(hard_coral, k = 3), data = herb.gam.df, Gamma(link = log))
plot(herb.gam.t1)

summary(lm(log(gravity.df$mean_bio_hectare) ~ log(gravity.df$mean_complexity)))
ggplot() +
  geom_point(data = gravity.df, aes(x = log(mean_complexity), y = log(mean_bio_hectare))) +
  geom_smooth(data = gravity.df, aes(x = log(mean_complexity), y = log(mean_bio_hectare)), method = "lm")

# Plot predicted values for trophic groups =============================================================

# Influence of biomass on slope - CARNIVORE
# carn.gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3), data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail")
newd.carn.gam3 <- data.frame(mean_bio_hectare = seq(min(carn.gam.df$mean_bio_hectare), max(carn.gam.df$mean_bio_hectare), length.out = 2000),
                             mean_complexity = seq(min(carn.gam.df$mean_complexity), max(carn.gam.df$mean_complexity), length.out = 2000))
pred <- predict.gam(carn.gam1.1, newdata = newd.carn.gam3, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.carn.gam3$b_pred <- exp(pred$fit) * -1
newd.carn.gam3$b_pred_upr <- exp(upr) * -1
newd.carn.gam3$b_pred_lwr <- exp(lwr) * -1
carn.plot1 <- ggplot() +
  geom_line(aes(x = mean_bio_hectare, y = b_pred), size = 1.5, data = newd.carn.gam3) +
  geom_ribbon(aes(x = mean_bio_hectare, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.carn.gam3) +
  geom_point(aes(x = mean_bio_hectare, y = b, shape = region), data = carn.gam.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.2)) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  labs(x = "Biomass (kg/ha)", y = expression(paste("Size spectra (", italic("b"), ")")), tag = "a")
  # labs(x = "Biomass kg/ha", y = "Size spectra (b)")
  
# Influence of structtural complexity on slope - CARNIVORE
# carn.gam4 <- gam(abs(b) ~ s(hard_coral, k=3), data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail")
# newd.carn.gam4 <- data.frame(hard_coral = seq(min(carn.gam.df$hard_coral), max(carn.gam.df$hard_coral), length.out = 2000))
# pred <- predict.gam(carn.gam4, newdata = newd.carn.gam4, se.fit = TRUE)
# upr <- pred$fit + (2 * pred$se.fit)
# lwr <- pred$fit - (2 * pred$se.fit)
# newd.carn.gam4$b_pred <- exp(pred$fit) * -1
# newd.carn.gam4$b_pred_upr <- exp(upr) * -1
# newd.carn.gam4$b_pred_lwr <- exp(lwr) * -1
carn.plot2 <- ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newd.carn.gam3) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.carn.gam3) +
  geom_point(aes(x = mean_complexity, y = b, shape = region), data = carn.gam.df, size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  labs(x = "Structural complexity", y = "", tag = "b")

carn.gam4 <- gam(mean_bio_hectare ~ s(mean_complexity, k=3), data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail")
newd.carn.gam4 <- data.frame(mean_complexity = seq(min(carn.gam.df$mean_complexity), max(carn.gam.df$mean_complexity), length.out = 2000))
pred <- predict.gam(carn.gam4, newdata = newd.carn.gam4, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.carn.gam4$pred <- exp(pred$fit)
newd.carn.gam4$pred_upr <- exp(upr)
newd.carn.gam4$pred_lwr <- exp(lwr)
carn.plot3 <- ggplot() +
  geom_line(aes(x = mean_complexity, y = pred), size = 1.5, data = newd.carn.gam4) +
  geom_ribbon(aes(x = mean_complexity, ymin = pred_lwr, ymax = pred_upr), alpha = 0.5, fill = "gray", data = newd.carn.gam4) +
  geom_point(aes(x = mean_complexity, y = mean_bio_hectare, shape = region), data = carn.gam.df, size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  labs(x = "Structural complexity", y = "Biomass (kg/ha)", tag = "c")

library(ggpubr)
ggarrange(carn.plot1,carn.plot2,carn.plot3,)

carn.gam4 <- gam(abs(b) ~ s(mean_complexity, k=3), data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail")
newd.carn.gam4 <- data.frame(mean_complexity = seq(min(carn.gam.df$mean_complexity), max(carn.gam.df$mean_complexity), length.out = 2000))
pred <- predict.gam(carn.gam4, newdata = newd.carn.gam4, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.carn.gam4$b_pred <- exp(pred$fit) * -1
newd.carn.gam4$b_pred_upr <- exp(upr) * -1
newd.carn.gam4$b_pred_lwr <- exp(lwr) * -1
ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newd.carn.gam4) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.carn.gam4) +
  geom_point(aes(x = mean_complexity, y = b, shape = region), data = carn.gam.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Size spectrum slope") +
  xlab("mean_complexity")

# Influence of biomass on slope - HERBIVORE
# herb.gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3), data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail")
newd.herb.gam3 <- data.frame(mean_bio_hectare = seq(min(herb.gam.df$mean_bio_hectare), max(herb.gam.df$mean_bio_hectare), length.out = 2000),
                            mean_complexity = seq(min(herb.gam.df$mean_complexity), max(herb.gam.df$mean_complexity), length.out = 2000))
pred <- predict.gam(herb.gam1.1, newdata = newd.herb.gam3, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.herb.gam3$b_pred <- exp(pred$fit) * -1
newd.herb.gam3$b_pred_upr <- exp(upr) * -1
newd.herb.gam3$b_pred_lwr <- exp(lwr) * -1
herb.plot1 <- ggplot() +
  geom_line(aes(x = mean_bio_hectare, y = b_pred), size = 1.5, data = newd.herb.gam3) +
  geom_ribbon(aes(x = mean_bio_hectare, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.herb.gam3) +
  geom_point(aes(x = mean_bio_hectare, y = b, shape = region), data = herb.gam.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.2)) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  labs(x = "Biomass (kg/ha)", y = expression(paste("Size spectra (", italic("b"), ")")), tag = "a")

# Influence of mean complexity on slope - HERBIVORE
# herb.gam4 <- gam(abs(b) ~ s(mean_complexity, k=3), data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail")
# newd.herb.gam4 <- data.frame(mean_complexity = seq(min(herb.gam.df$mean_complexity), max(herb.gam.df$mean_complexity), length.out = 2000))
# pred <- predict.gam(herb.gam4, newdata = newd.herb.gam4, se.fit = TRUE)
# upr <- pred$fit + (2 * pred$se.fit)
# lwr <- pred$fit - (2 * pred$se.fit)
# newd.herb.gam4$b_pred <- exp(pred$fit) * -1
# newd.herb.gam4$b_pred_upr <- exp(upr) * -1
# newd.herb.gam4$b_pred_lwr <- exp(lwr) * -1
herb.plot2 <- ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newd.herb.gam3) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.herb.gam3) +
  geom_point(aes(x = mean_complexity, y = b, shape = region), data = herb.gam.df, size = 2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  labs(x = "Structural complexity", y = expression(paste("Size spectra (", italic("b"), ")")), tag = "b")

ggarrange(herb.plot1, herb.plot2, nrow = 2)

# SUPPLEMENTAL FIGURES =================================================================================

# Gravity against biomass -------------------------------------------------
gravity.df <- gravity.df %>%
  mutate(Grav_tot = as.numeric(Grav_tot)) %>%
  mutate(mean_bio_hectare = as.numeric(mean_bio_hectare))

lm.1 <- lm(log(mean_bio_hectare) ~ log(Grav_tot), data = gravity.df)
summary(lm.1)

ggplot() +
  geom_point(data = gravity.df, aes(x = log(Grav_tot), y = log(mean_bio_hectare), color = region)) +
  geom_smooth(data = gravity.df, aes(x = log(Grav_tot), y = log(mean_bio_hectare)), method = "lm") +
  theme_classic() +
  labs(x = "log(population gravity)", y = "log(biomass density [kg/ha])") +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok"))

# Fish biomass per region -------------------------------------------------
bio_region <- fish.df %>%
  group_by(region, site_name, transect, observer) %>%
  summarize(biomass_kg = sum(biomass_kg)) %>%
  group_by(region, site_name) %>%
  summarize(biomass_kg = mean(biomass_kg)) %>%
  group_by(region) %>%
  summarize(mean_bio = mean(biomass_kg),
            stderr_bio = std.error(biomass_kg)) %>%
  mutate(mean_bio_ha = mean_bio * 40) %>%
  mutate(stderr_bio_ha = stderr_bio * 40)

bio_region <- bio_region %>% mutate(region = fct_relevel(region, "raja_ampat", "wakatobi", "lombok"))

ggplot() +
  geom_pointrange(data = bio_region, aes(x = region, y = mean_bio_ha, ymin = mean_bio_ha - stderr_bio_ha, ymax = mean_bio_ha + stderr_bio_ha)) +
  theme_classic() +
  labs(x = "Region", y = "Mean biomass (kg/ha)") +
  scale_x_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok"))

# Herbivore biomass against hard coral cover ------------------------------

# hard coral cover and herbivore biomass
herb.gam.t1 <- gam(mean_bio_hectare ~ s(hard_coral, k = 3), data = herb.gam.df, Gamma(link = log))
summary(herb.gam.t1)
newd.herb.gam.t1 <- data.frame(hard_coral = seq(min(herb.gam.df$hard_coral), max(herb.gam.df$hard_coral), length.out = 2000))
pred <- predict.gam(herb.gam.t1, newdata = newd.herb.gam.t1, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.herb.gam.t1$b_pred <- exp(pred$fit)
newd.herb.gam.t1$b_pred_upr <- exp(upr)
newd.herb.gam.t1$b_pred_lwr <- exp(lwr)
ggplot() +
  geom_line(aes(x = hard_coral, y = b_pred), size = 1.5, data = newd.herb.gam.t1) +
  geom_ribbon(aes(x = hard_coral, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.herb.gam.t1) +
  geom_point(aes(x = hard_coral, y = mean_bio_hectare, shape = region), data = herb.gam.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Biomass density (kg/ha)") +
  xlab("Hard coral cover (%)") +
  scale_x_continuous(expand = c(0,0))

plot(herb.gam.t1)


# Overall fish biomass against hard coral cover ---------------------------

gam.t1 <- gam(mean_bio_hectare ~ s(hard_coral, k = 3), data = gravity.df, Gamma(link = log))
newd.herb.gam.t1 <- data.frame(hard_coral = seq(min(herb.gam.df$hard_coral), max(herb.gam.df$hard_coral), length.out = 2000))
pred <- predict.gam(herb.gam.t1, newdata = newd.herb.gam.t1, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.herb.gam.t1$b_pred <- exp(pred$fit)
newd.herb.gam.t1$b_pred_upr <- exp(upr)
newd.herb.gam.t1$b_pred_lwr <- exp(lwr)
ggplot() +
  geom_line(aes(x = hard_coral, y = b_pred), size = 1.5, data = newd.herb.gam.t1) +
  geom_ribbon(aes(x = hard_coral, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.herb.gam.t1) +
  geom_point(aes(x = hard_coral, y = mean_bio_hectare, shape = region), data = herb.gam.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Biomass density (kg/ha)") +
  xlab("Hard coral cover (%)") +
  scale_x_continuous(expand = c(0,0))

summary(gam.t1)


# Fish biomass and structural complexity ----------------------------------

gam.t2 <- gam(mean_complexity ~ s(mean_bio_hectare, k = 3), data = gravity.df, Gamma(link = log))
lm.t2 <- lm(mean_complexity ~ mean_bio_hectare, data = gravity.df)
newd.gam.t2 <- data.frame(mean_bio_hectare = seq(min(gravity.df$mean_bio_hectare), max(gravity.df$mean_bio_hectare), length.out = 2000))
pred <- predict.gam(gam.t2, newdata = newd.gam.t2, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.gam.t2$b_pred <- exp(pred$fit)
newd.gam.t2$b_pred_upr <- exp(upr)
newd.gam.t2$b_pred_lwr <- exp(lwr)
ggplot() +
  geom_line(aes(x = mean_bio_hectare, y = b_pred), size = 1.5, data = newd.gam.t2) +
  geom_ribbon(aes(x = mean_bio_hectare, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.gam.t2) +
  geom_point(aes(x = mean_bio_hectare, y = mean_complexity, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  xlab("Biomass density (kg/ha)") +
  ylab("Structural complexity") +
  scale_x_continuous(expand = c(0,0))

AIC(gam.t2, lm.t2)


# Fish biomass and algal cover --------------------------------------------

gam.t3 <- gam(mean_bio_hectare ~ s(algae, k = 3), data = gravity.df, Gamma(link = log))
newd.gam.t3 <- data.frame(algae = seq(min(gravity.df$algae), max(gravity.df$algae), length.out = 2000))
pred <- predict.gam(gam.t3, newdata = newd.gam.t3, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.gam.t3$b_pred <- exp(pred$fit)
newd.gam.t3$b_pred_upr <- exp(upr)
newd.gam.t3$b_pred_lwr <- exp(lwr)
ggplot() +
  geom_line(aes(x = algae, y = b_pred), size = 1.5, data = newd.gam.t3) +
  geom_ribbon(aes(x = algae, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.gam.t3) +
  geom_point(aes(x = algae, y = mean_bio_hectare, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Biomass density (kg/ha)") +
  xlab("Algae cover") +
  scale_x_continuous(expand = c(0,0))

summary(gam.t3)

# Hard coral and algal cover --------------------------------------------

gam.t4 <- gam(hard_coral ~ s(algae, k = 3), data = gravity.df, Gamma(link = log))
newd.gam.t4 <- data.frame(algae = seq(min(gravity.df$algae), max(gravity.df$algae), length.out = 2000))
pred <- predict.gam(gam.t4, newdata = newd.gam.t4, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.gam.t4$b_pred <- exp(pred$fit)
newd.gam.t4$b_pred_upr <- exp(upr)
newd.gam.t4$b_pred_lwr <- exp(lwr)
ggplot() +
  geom_line(aes(x = algae, y = b_pred), size = 1.5, data = newd.gam.t4) +
  geom_ribbon(aes(x = algae, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.gam.t4) +
  geom_point(aes(x = algae, y = hard_coral, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Hard coral cover") +
  xlab("Algae cover") +
  scale_x_continuous(expand = c(0,0))

summary(gam.t4)

# Hard coral and structural complexity --------------------------------------------

gam.t5 <- gam(hard_coral ~ s(mean_complexity, k = 3), data = gravity.df, Gamma(link = log))
newd.gam.t5 <- data.frame(mean_complexity = seq(min(gravity.df$mean_complexity), max(gravity.df$mean_complexity), length.out = 2000))
pred <- predict.gam(gam.t5, newdata = newd.gam.t5, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.gam.t5$b_pred <- exp(pred$fit)
newd.gam.t5$b_pred_upr <- exp(upr)
newd.gam.t5$b_pred_lwr <- exp(lwr)
ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newd.gam.t5) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.gam.t5) +
  geom_point(aes(x = mean_complexity, y = hard_coral, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Hard coral cover") +
  xlab("Structural complexity") +
  scale_x_continuous(expand = c(0,0))

summary(gam.t5)

# Algae and structural complexity --------------------------------------------

gam.t6 <- gam(algae ~ s(mean_complexity, k = 10), data = gravity.df)
newd.gam.t6 <- data.frame(mean_complexity = seq(min(gravity.df$mean_complexity), max(gravity.df$mean_complexity), length.out = 2000))
pred <- predict.gam(gam.t6, newdata = newd.gam.t6, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.gam.t6$b_pred <- exp(pred$fit)
newd.gam.t6$b_pred_upr <- exp(upr)
newd.gam.t6$b_pred_lwr <- exp(lwr)
ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newd.gam.t6) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.gam.t6) +
  geom_point(aes(x = mean_complexity, y = algae, shape = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Algae") +
  xlab("Structural complexity") +
  scale_x_continuous(expand = c(0,0))

summary(gam.t6)


# Carnivore biomass against structural complexity -------------------------

# structural complexity and carnivore biomass
carn.gam.t2 <- gam(mean_bio_hectare ~ s(mean_complexity, k = 3), data = carn.gam.df, Gamma(link = log))
newd.carn.gam.t2 <- data.frame(mean_complexity = seq(min(carn.gam.df$mean_complexity), max(carn.gam.df$mean_complexity), length.out = 2000))
pred <- predict.gam(carn.gam.t2, newdata = newd.carn.gam.t2, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd.carn.gam.t2$b_pred <- exp(pred$fit)
newd.carn.gam.t2$b_pred_upr <- exp(upr)
newd.carn.gam.t2$b_pred_lwr <- exp(lwr)
ggplot() +
  geom_line(aes(x = hard_coral, y = b_pred), size = 1.5, data = newd.herb.gam.t1) +
  geom_ribbon(aes(x = hard_coral, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd.herb.gam.t1) +
  geom_point(aes(x = hard_coral, y = mean_bio_hectare, shape = region), data = herb.gam.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_shape_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Biomass density (kg/ha)") +
  xlab("Hard coral cover (%)") +
  scale_x_continuous(expand = c(0,0))

plot(herb.gam.t1)


xmax = 65
xmin = 10
b = -1.5
n = rev(seq(1,1000,1))
x = seq(10,65,length.out = 1000)

nx = (n * x^b) * ((b + 1) / (xmax^(b+1) - xmin^(b+1)))

plot(x,nx)


# MLE bin method from Edwards et al 2020 ----------------------------------

# need to fix a and b params for some species
fish.df$a[which(fish.df$genus_species == "Zanclus cornutus")] <- 0.0147
fish.df$b[which(fish.df$genus_species == "Zanclus cornutus")] <- 3.37
fish.df$a[which(fish.df$genus_species == "Halichoeres melanurus")] <- 0.01091
fish.df$b[which(fish.df$genus_species == "Halichoeres melanurus")] <- 3.000
fish.df$a[which(fish.df$genus_species == "Hemigymnus fasciatus")] <- 0.01713
fish.df$b[which(fish.df$genus_species == "Hemigymnus fasciatus")] <- 3.000
fish.df$a[which(fish.df$genus_species == "Hemigymnus melapterus")] <- 0.01816
fish.df$b[which(fish.df$genus_species == "Hemigymnus melapterus")] <- 3.000
fish.df$a[which(fish.df$genus_species == "Thalassoma spp.")] <- 0.01830
fish.df$b[which(fish.df$genus_species == "Thalassoma spp.")] <- 2.862


# MLE bin - REGIONS -------------------------------------------------------

MLEbin.data.region <- fish.df %>%
  dplyr::select(site_name, region, genus_species, trophic.position = tp, abundance, size_cm, a, b, biomass_g) %>%
  dplyr::group_by(region, genus_species, a, b, size_cm) %>%
  dplyr::summarise(abundance = sum(abundance)) %>%
  dplyr::mutate(body.mass = a * size_cm ^ b) %>%
  dplyr::ungroup()

dataBin <- dplyr::mutate(MLEbin.data.region,
                         LngtMax = size_cm + 1)

dataBin$SpecCode <- as.integer(dataBin$genus_species)

dataBin <- dataBin %>%
  mutate(regionCode = if_else(region == "raja_ampat", 1,
                      if_else(region == "wakatobi", 2, 3)))

dataBin <- dplyr::mutate(dataBin, wmax = a * LngtMax ^ b)
dataBin <- dplyr::rename(dataBin, LngtMin = size_cm)
dataBin <- dplyr::rename(dataBin, wmin = body.mass)
dataBin <- dataBin %>% dplyr::select(Year = regionCode, SpecCode, LngtMin, LngtMax, LWa = a, 
                                     LWb= b, wmin, wmax, Number = abundance)

fullYears = sort(unique(dataBin$Year))
for(iii in 1:length(fullYears)){
  dataBinForLike <- dplyr::filter(dataBin,
                                  Year == fullYears[iii])
  dataBinForLike <- dplyr::select(dataBinForLike, SpecCode, wmin, wmax, Number)
  n = sum(dataBinForLike$Number)
  xmin = min(dataBinForLike$wmin)
  xmax = max(dataBinForLike$wmax)
  
  MLEbins.region.new = calcLike(negLL.fn = negLL.PLB.bins.species,
                                p = -1.9,
                                suppress.warnings = TRUE,
                                dataBinForLike = dataBinForLike,
                                n = n,
                                xmin = xmin,
                                xmax = xmax)
  if(iii == 1){
    MLEbins.region = data.frame(Year = fullYears[iii],
                                    xmin = xmin,
                                    xmax = xmax,
                                    n = n,
                                    b = MLEbins.region.new$MLE,
                                    confMin = MLEbins.region.new$conf[1],
                                    confMax = MLEbins.region.new$conf[2])
  } else {
    MLEbins.region = rbind(MLEbins.region,
                           c(fullYears[iii],
                             xmin, xmax, n, MLEbins.region.new$MLE, MLEbins.region.new$conf[1], MLEbins.region.new$conf[2]))
  }
}

MLEbins.region <- dplyr::rename(MLEbins.region, Region = Year) 
MLEbins.region <- MLEbins.region %>%
  mutate(Region = if_else(Region == 1, "Raja Ampat", if_else(Region == 2, "Wakatobi", "Lombok")))



ggplot() +
  geom_point(data = MLEbins.region, aes(x = Region, y = b)) +
  geom_errorbar(data = MLEbins.region, aes(x = Region, ymin = confMin, ymax = confMax), width= 0.1) +
  theme_classic() +
  labs(y = "Estimate of b", x = "")


# MLE bin - STUDY SITES ---------------------------------------------------

MLEbin.data.site <- fish.df %>%
  dplyr::select(site_name, genus_species, trophic.position = tp, abundance, size_cm, a, b, biomass_g) %>%
  dplyr::group_by(site_name, genus_species, a, b, size_cm) %>%
  dplyr::summarise(abundance = sum(abundance)) %>%
  dplyr::mutate(body.mass = a * size_cm ^ b) %>%
  dplyr::ungroup()

dataBin.site <- dplyr::mutate(MLEbin.data.site,
                         LngtMax = size_cm + 1)

dataBin.site$SpecCode <- as.integer(dataBin.site$genus_species)
dataBin.site$siteCode <- as.integer(dataBin.site$site_name)

# need key to match site codes with site names and site names with regions
region.site.key <- fish.df %>%
  dplyr::select(site_name, region) %>%
  dplyr::distinct(site_name, region)
site.code.key <- dataBin.site %>%
  dplyr::select(site_name, siteCode) %>%
  dplyr::distinct(site_name, siteCode) %>%
  dplyr::left_join(., region.site.key, by = "site_name")

dataBin.site <- dplyr::mutate(dataBin.site, wmax = a * LngtMax ^ b)
dataBin.site <- dplyr::rename(dataBin.site, LngtMin = size_cm)
dataBin.site <- dplyr::rename(dataBin.site, wmin = body.mass)
dataBin.site <- dataBin.site %>% dplyr::select(Year = siteCode, SpecCode, LngtMin, LngtMax, LWa = a, 
                                     LWb= b, wmin, wmax, Number = abundance)

fullYears = sort(unique(dataBin.site$Year))
for(iii in 1:length(fullYears)){
  dataBinForLike <- dplyr::filter(dataBin.site,
                                  Year == fullYears[iii])
  dataBinForLike <- dplyr::select(dataBinForLike, SpecCode, wmin, wmax, Number)
  n = sum(dataBinForLike$Number)
  xmin = min(dataBinForLike$wmin)
  xmax = max(dataBinForLike$wmax)
  
  MLEbins.site.new = calcLike(negLL.fn = negLL.PLB.bins.species,
                                p = -1.9,
                                suppress.warnings = TRUE,
                                dataBinForLike = dataBinForLike,
                                n = n,
                                xmin = xmin,
                                xmax = xmax)
  if(iii == 1){
    MLEbins.site = data.frame(Year = fullYears[iii],
                                    xmin = xmin,
                                    xmax = xmax,
                                    n = n,
                                    b = MLEbins.site.new$MLE,
                                    confMin = MLEbins.site.new$conf[1],
                                    confMax = MLEbins.site.new$conf[2])
  } else {
    MLEbins.site = rbind(MLEbins.site,
                           c(fullYears[iii],
                             xmin, xmax, n, MLEbins.site.new$MLE, MLEbins.site.new$conf[1], MLEbins.site.new$conf[2]))
  }
}

MLEbins.site <- dplyr::rename(MLEbins.site, site = Year)
MLEbins.site <- dplyr::rename(MLEbins.site, siteCode = site)
MLEbins.site <- MLEbins.site %>%
  dplyr::left_join(., site.code.key, by = "siteCode") %>%
  dplyr::mutate(region = fct_relevel(region, "raja_ampat", "wakatobi", "lombok"))
mean.b <- mean(MLEbins.site$b)

# find which size spectra are not equal to the average b
MLEbins.site.arrange <- MLEbins.site %>% dplyr::arrange(b)
MLEbins.site.arrange <- MLEbins.site.arrange %>%
  dplyr::mutate(test = if_else(confMin <= mean.b & confMax >= mean.b, "red", "black"))
prop.below <- MLEbins.site.arrange %>% 
  dplyr::filter(b < mean.b & test == "black") %>%
  dplyr::mutate(num = 1) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(num = sum(num)) %>%
  dplyr::mutate(perc = round(num/sum(num),2) * 100)
prop.above <- MLEbins.site.arrange %>% 
  dplyr::filter(b > mean.b & test == "black") %>%
  dplyr::mutate(num = 1) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(num = sum(num)) %>%
  dplyr::mutate(perc = round(num/sum(num),2) * 100)

ggplot() +
  geom_point(data = MLEbins.site.arrange, aes(x = reorder(site_name, -b), y = b, color = region), size = 2) +
  geom_errorbar(data = MLEbins.site.arrange, aes(x = reorder(site_name, -b), ymin = confMin, ymax = confMax), 
                color = MLEbins.site.arrange$test, width= 0.1) +
  theme_classic() +
  labs(y = "Estimate of b", x = "") +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = mean.b, linetype = "dashed", color = "grey", size = 1) +
  coord_flip()

b.sites.df <- MLEbins.site %>% dplyr::select(site_name, b.bins = b)

gravity.df <- gravity.df %>%
  dplyr::left_join(., b.sites.df, by = "site_name")

# rerun GAM models

# Models with Gamma family and log link
bin.gam1 <- gam(abs(b.bins) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3),
          data = gravity.df, family = Gamma(link=log), na.action = "na.fail")
dd <- dredge(bin.gam1, rank = "AICc")

# Check residuals to make sure the model meets assumptions
bin.gam2 <- gam(abs(b.bins) ~ s(mean_bio_hectare, k=3) + s(mean_complexity, k=3), data = gravity.df, family = Gamma(link=log))
summary(bin.gam2)
par(mfrow = c(2,2))
gam.check(bin.gam2)

bin.gam3 <- gam(abs(b.bins) ~ s(mean_bio_hectare, k=3), data = gravity.df, family = Gamma(link=log))
summary(bin.gam3)
par(mfrow = c(2,2))
gam.check(bin.gam3)

# Plot predicted values

newd <- data.frame(mean_bio_hectare = seq(min(gravity.df$mean_bio_hectare), max(gravity.df$mean_bio_hectare), length.out = 2000))
pred <- predict.gam(bin.gam3, newd, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd$b_pred <- exp(pred$fit) * -1
newd$b_pred_upr <- exp(upr) * -1
newd$b_pred_lwr <- exp(lwr) * -1
tmp_grav <- gravity.df
gravity.df <- tmp_grav
gravity.df$region <- factor(gravity.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = mean_bio_hectare, y = b_pred), size = 1.5, data = newd) +
  geom_ribbon(aes(x = mean_bio_hectare, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd) +
  geom_point(aes(x = mean_bio_hectare, y = b, color = region), data = gravity.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Size spectra (b)") +
  xlab("Biomass (kg/ha)")

# MLE bin - CARN AND HERB -------------------------------------------------

MLEbin.data.trophic <- fish.df %>%
  dplyr::select(site_name, region, genus_species, trophic.position = tp, abundance, size_cm, a, b, biomass_g) %>%
  dplyr::group_by(trophic.position, genus_species, a, b, size_cm) %>%
  dplyr::summarise(abundance = sum(abundance)) %>%
  dplyr::mutate(body.mass = a * size_cm ^ b) %>%
  dplyr::ungroup()

dataBin.trophic <- dplyr::mutate(MLEbin.data.trophic,
                         LngtMax = size_cm + 1)

dataBin.trophic$SpecCode <- as.integer(dataBin.trophic$genus_species)

dataBin.trophic <- dataBin.trophic %>%
  mutate(trophicCode = if_else(trophic.position == "Carnivore", 1,
                       if_else(trophic.position == "Herbivore", 2, 3)))

dataBin.trophic <- dplyr::mutate(dataBin.trophic, wmax = a * LngtMax ^ b)
dataBin.trophic <- dplyr::rename(dataBin.trophic, LngtMin = size_cm)
dataBin.trophic <- dplyr::rename(dataBin.trophic, wmin = body.mass)
dataBin.trophic <- dataBin.trophic %>% dplyr::select(Year = trophicCode, SpecCode, LngtMin, LngtMax, LWa = a, 
                                     LWb= b, wmin, wmax, Number = abundance)

fullYears = sort(unique(dataBin.trophic$Year))
for(iii in 1:length(fullYears)){
  dataBinForLike <- dplyr::filter(dataBin.trophic,
                                  Year == fullYears[iii])
  dataBinForLike <- dplyr::select(dataBinForLike, SpecCode, wmin, wmax, Number)
  n = sum(dataBinForLike$Number)
  xmin = min(dataBinForLike$wmin)
  xmax = max(dataBinForLike$wmax)
  
  MLEbins.trophic.new = calcLike(negLL.fn = negLL.PLB.bins.species,
                                p = -1.9,
                                suppress.warnings = TRUE,
                                dataBinForLike = dataBinForLike,
                                n = n,
                                xmin = xmin,
                                xmax = xmax)
  if(iii == 1){
    MLEbins.trophic = data.frame(Year = fullYears[iii],
                                    xmin = xmin,
                                    xmax = xmax,
                                    n = n,
                                    b = MLEbins.trophic.new$MLE,
                                    confMin = MLEbins.trophic.new$conf[1],
                                    confMax = MLEbins.trophic.new$conf[2])
  } else {
    MLEbins.trophic = rbind(MLEbins.trophic,
                           c(fullYears[iii],
                             xmin, xmax, n, MLEbins.trophic.new$MLE, MLEbins.trophic.new$conf[1], MLEbins.trophic.new$conf[2]))
  }
}

MLEbins.trophic <- dplyr::rename(MLEbins.trophic, trophic.position = Year) 
MLEbins.trophic <- MLEbins.trophic %>%
  mutate(trophic.position = if_else(trophic.position == 1, "Carnivore", if_else(trophic.position == 2, "Herbivore", "NA")))

ggplot() +
  geom_point(data = MLEbins.trophic, aes(x = trophic.position, y = b)) +
  geom_errorbar(data = MLEbins.trophic, aes(x = trophic.position, ymin = confMin, ymax = confMax), width= 0.1) +
  theme_classic() +
  labs(y = "Estimate of b", x = "")


# MLE bin - GAMS CARNS ------------------------------------------

MLEbin.data.site.carn <- fish.df %>%
  dplyr::filter(tp == "Carnivore") %>%
  dplyr::select(site_name, genus_species, trophic.position = tp, abundance, size_cm, a, b, biomass_g) %>%
  dplyr::group_by(site_name, genus_species, a, b, size_cm) %>%
  dplyr::summarise(abundance = sum(abundance)) %>%
  dplyr::mutate(body.mass = a * size_cm ^ b) %>%
  dplyr::ungroup()

dataBin.site.carn <- dplyr::mutate(MLEbin.data.site.carn,
                         LngtMax = size_cm + 1)

dataBin.site.carn$SpecCode <- as.integer(dataBin.site.carn$genus_species)
dataBin.site.carn$siteCode <- as.integer(dataBin.site.carn$site_name)

# need key to match site codes with site names and site names with regions
region.site.key <- fish.df %>%
  dplyr::select(site_name, region) %>%
  dplyr::distinct(site_name, region)
site.code.key <- dataBin.site.carn %>%
  dplyr::select(site_name, siteCode) %>%
  dplyr::distinct(site_name, siteCode) %>%
  dplyr::left_join(., region.site.key, by = "site_name")

dataBin.site.carn <- dplyr::mutate(dataBin.site.carn, wmax = a * LngtMax ^ b)
dataBin.site.carn <- dplyr::rename(dataBin.site.carn, LngtMin = size_cm)
dataBin.site.carn <- dplyr::rename(dataBin.site.carn, wmin = body.mass)
dataBin.site.carn <- dataBin.site.carn %>% dplyr::select(Year = siteCode, SpecCode, LngtMin, LngtMax, LWa = a, 
                                     LWb= b, wmin, wmax, Number = abundance)

fullYears = sort(unique(dataBin.site.carn$Year))
for(iii in 1:length(fullYears)){
  dataBinForLike <- dplyr::filter(dataBin.site.carn,
                                  Year == fullYears[iii])
  dataBinForLike <- dplyr::select(dataBinForLike, SpecCode, wmin, wmax, Number)
  n = sum(dataBinForLike$Number)
  xmin = min(dataBinForLike$wmin)
  xmax = max(dataBinForLike$wmax)
  
  MLEbins.site.carn.new = calcLike(negLL.fn = negLL.PLB.bins.species,
                                p = -1.9,
                                suppress.warnings = TRUE,
                                dataBinForLike = dataBinForLike,
                                n = n,
                                xmin = xmin,
                                xmax = xmax,
                                vecDiff = 0.6)
  if(iii == 1){
    MLEbins.site.carn = data.frame(Year = fullYears[iii],
                                    xmin = xmin,
                                    xmax = xmax,
                                    n = n,
                                    b = MLEbins.site.carn.new$MLE,
                                    confMin = MLEbins.site.carn.new$conf[1],
                                    confMax = MLEbins.site.carn.new$conf[2])
  } else {
    MLEbins.site.carn = rbind(MLEbins.site.carn,
                           c(fullYears[iii],
                             xmin, xmax, n, MLEbins.site.carn.new$MLE, MLEbins.site.carn.new$conf[1], MLEbins.site.carn.new$conf[2]))
  }
}

MLEbins.site.carn <- dplyr::rename(MLEbins.site.carn, site = Year)
MLEbins.site.carn <- dplyr::rename(MLEbins.site.carn, siteCode = site)
MLEbins.site.carn <- MLEbins.site.carn %>%
  dplyr::left_join(., site.code.key, by = "siteCode") %>%
  dplyr::mutate(region = fct_relevel(region, "raja_ampat", "wakatobi", "lombok"))
mean.carn.b <- mean(MLEbins.site.carn$b)

# find which size spectra are not equal to the average b
MLEbins.site.arrange.carn <- MLEbins.site.carn %>% dplyr::arrange(b)
MLEbins.site.arrange.carn <- MLEbins.site.arrange.carn %>%
  dplyr::mutate(test = if_else(confMin <= mean.carn.b & confMax >= mean.carn.b, "red", "black"))
prop.below <- MLEbins.site.arrange.carn %>% 
  dplyr::filter(b < mean.carn.b & test == "black") %>%
  dplyr::mutate(num = 1) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(num = sum(num)) %>%
  dplyr::mutate(perc = round(num/sum(num),2) * 100)
prop.above <- MLEbins.site.arrange.carn %>% 
  dplyr::filter(b > mean.carn.b & test == "black") %>%
  dplyr::mutate(num = 1) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(num = sum(num)) %>%
  dplyr::mutate(perc = round(num/sum(num),2) * 100)

ggplot() +
  geom_point(data = MLEbins.site.arrange.carn, aes(x = reorder(site_name, -b), y = b, color = region), size = 2) +
  geom_errorbar(data = MLEbins.site.arrange.carn, aes(x = reorder(site_name, -b), ymin = confMin, ymax = confMax), 
                color = MLEbins.site.arrange.carn$test, width= 0.1) +
  theme_classic() +
  labs(y = "Estimate of b", x = "", title = "Carnivores") +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = mean.carn.b, linetype = "dashed", color = "grey", size = 1) +
  coord_flip()

b.sites.carn.df <- MLEbins.site.carn %>% dplyr::select(site_name, b.bins.carn = b)

# get carnivore biomass from carn.gam.df
carn_bio <- carn.gam.df %>% dplyr::select(site_name, mean_bio_carn = mean_bio_hectare)

gravity.carn.df <- gravity.df %>%
  dplyr::left_join(., b.sites.carn.df, by = "site_name") %>%
  dplyr::left_join(., carn_bio, by = "site_name")

# rerun GAM models

# Models with Gamma family and log link
bin.carn.gam1 <- gam(abs(b.bins.carn) ~ s(mean_bio_carn, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3),
                    data = gravity.carn.df, family = Gamma(link=log), na.action = "na.fail")
dd <- dredge(bin.carn.gam1, rank = "AICc")

# Check residuals to make sure the model meets assumptions
bin.carn.gam2 <- gam(abs(b.bins.carn) ~ s(mean_bio_carn, k=3) + s(hard_coral, k=3), data = gravity.carn.df, 
                    family = Gamma(link=log))
summary(bin.carn.gam2)
par(mfrow = c(2,2))
gam.check(bin.carn.gam2)

bin.carn.gam3 <- gam(abs(b.bins.carn) ~ s(mean_bio_carn, k=3), data = gravity.carn.df, 
                    family = Gamma(link=log))
summary(bin.carn.gam3)
par(mfrow = c(2,2))
gam.check(bin.carn.gam3)

bin.carn.gam4 <- gam(abs(b.bins.carn) ~ s(hard_coral, k=3), data = gravity.carn.df, 
                    family = Gamma(link=log))
summary(bin.carn.gam4)
par(mfrow = c(2,2))
gam.check(bin.carn.gam4)

# Plot predicted values

newd <- data.frame(mean_bio_carn = seq(min(gravity.carn.df$mean_bio_carn), max(gravity.carn.df$mean_bio_carn), length.out = 2000))
pred <- predict.gam(bin.carn.gam3, newd, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd$b_pred <- exp(pred$fit) * -1
newd$b_pred_upr <- exp(upr) * -1
newd$b_pred_lwr <- exp(lwr) * -1
tmp_grav <- gravity.carn.df
gravity.carn.df$region <- factor(gravity.carn.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = mean_bio_carn, y = b_pred), size = 1.5, data = newd) +
  geom_ribbon(aes(x = mean_bio_carn, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd) +
  geom_point(aes(x = mean_bio_carn, y = b.bins.carn, color = region), data = gravity.carn.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Size spectra (b)") +
  xlab("Biomass (kg/ha)")

newd <- data.frame(hard_coral = seq(min(gravity.carn.df$hard_coral), max(gravity.carn.df$hard_coral), length.out = 2000))
pred <- predict.gam(bin.carn.gam4, newd, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd$b_pred <- exp(pred$fit) * -1
newd$b_pred_upr <- exp(upr) * -1
newd$b_pred_lwr <- exp(lwr) * -1
tmp_grav <- gravity.carn.df
gravity.carn.df$region <- factor(gravity.carn.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = hard_coral, y = b_pred), size = 1.5, data = newd) +
  geom_ribbon(aes(x = hard_coral, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd) +
  geom_point(aes(x = hard_coral, y = b.bins.carn, color = region), data = gravity.carn.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Size spectra (b)") +
  xlab("Hard coral cover (%)")

# MLE bin - GAMS HERBS ------------------------------------------

MLEbin.data.site.herb <- fish.df %>%
  dplyr::filter(tp == "Herbivore") %>%
  dplyr::select(site_name, genus_species, trophic.position = tp, abundance, size_cm, a, b, biomass_g) %>%
  dplyr::group_by(site_name, genus_species, a, b, size_cm) %>%
  dplyr::summarise(abundance = sum(abundance)) %>%
  dplyr::mutate(body.mass = a * size_cm ^ b) %>%
  dplyr::ungroup()

dataBin.site.herb <- dplyr::mutate(MLEbin.data.site.herb,
                         LngtMax = size_cm + 1)

dataBin.site.herb$SpecCode <- as.integer(dataBin.site.herb$genus_species)
dataBin.site.herb$siteCode <- as.integer(dataBin.site.herb$site_name)

# need key to match site codes with site names and site names with regions
region.site.key <- fish.df %>%
  dplyr::select(site_name, region) %>%
  dplyr::distinct(site_name, region)
site.code.key <- dataBin.site.herb %>%
  dplyr::select(site_name, siteCode) %>%
  dplyr::distinct(site_name, siteCode) %>%
  dplyr::left_join(., region.site.key, by = "site_name")

dataBin.site.herb <- dplyr::mutate(dataBin.site.herb, wmax = a * LngtMax ^ b)
dataBin.site.herb <- dplyr::rename(dataBin.site.herb, LngtMin = size_cm)
dataBin.site.herb <- dplyr::rename(dataBin.site.herb, wmin = body.mass)
dataBin.site.herb <- dataBin.site.herb %>% dplyr::select(Year = siteCode, SpecCode, LngtMin, LngtMax, LWa = a, 
                                     LWb= b, wmin, wmax, Number = abundance)

fullYears = sort(unique(dataBin.site.herb$Year))
for(iii in 1:length(fullYears)){
  dataBinForLike <- dplyr::filter(dataBin.site.herb,
                                  Year == fullYears[iii])
  dataBinForLike <- dplyr::select(dataBinForLike, SpecCode, wmin, wmax, Number)
  n = sum(dataBinForLike$Number)
  xmin = min(dataBinForLike$wmin)
  xmax = max(dataBinForLike$wmax)
  
  MLEbins.site.herb.new = calcLike(negLL.fn = negLL.PLB.bins.species,
                                p = -1.9,
                                suppress.warnings = TRUE,
                                dataBinForLike = dataBinForLike,
                                n = n,
                                xmin = xmin,
                                xmax = xmax,
                                vecDiff = 0.6)
  if(iii == 1){
    MLEbins.site.herb = data.frame(Year = fullYears[iii],
                                    xmin = xmin,
                                    xmax = xmax,
                                    n = n,
                                    b = MLEbins.site.herb.new$MLE,
                                    confMin = MLEbins.site.herb.new$conf[1],
                                    confMax = MLEbins.site.herb.new$conf[2])
  } else {
    MLEbins.site.herb = rbind(MLEbins.site.herb,
                           c(fullYears[iii],
                             xmin, xmax, n, MLEbins.site.herb.new$MLE, MLEbins.site.herb.new$conf[1], MLEbins.site.herb.new$conf[2]))
  }
}

MLEbins.site.herb <- dplyr::rename(MLEbins.site.herb, site = Year)
MLEbins.site.herb <- dplyr::rename(MLEbins.site.herb, siteCode = site)
MLEbins.site.herb <- MLEbins.site.herb %>%
  dplyr::left_join(., site.code.key, by = "siteCode") %>%
  dplyr::mutate(region = fct_relevel(region, "raja_ampat", "wakatobi", "lombok"))
mean.herb.b <- mean(MLEbins.site.herb$b)

# find which size spectra are not equal to the average b
MLEbins.site.arrange.herb <- MLEbins.site.herb %>% dplyr::arrange(b)
MLEbins.site.arrange.herb <- MLEbins.site.arrange.herb %>%
  dplyr::mutate(test = if_else(confMin <= mean.herb.b & confMax >= mean.herb.b, "red", "black"))
prop.below <- MLEbins.site.arrange.herb %>% 
  dplyr::filter(b < mean.herb.b & test == "black") %>%
  dplyr::mutate(num = 1) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(num = sum(num)) %>%
  dplyr::mutate(perc = round(num/sum(num),2) * 100)
prop.above <- MLEbins.site.arrange.herb %>% 
  dplyr::filter(b > mean.herb.b & test == "black") %>%
  dplyr::mutate(num = 1) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(num = sum(num)) %>%
  dplyr::mutate(perc = round(num/sum(num),2) * 100)

ggplot() +
  geom_point(data = MLEbins.site.arrange.herb, aes(x = reorder(site_name, -b), y = b, color = region), size = 2) +
  geom_errorbar(data = MLEbins.site.arrange.herb, aes(x = reorder(site_name, -b), ymin = confMin, ymax = confMax), 
                color = MLEbins.site.arrange.herb$test, width= 0.1) +
  theme_classic() +
  labs(y = "Estimate of b", x = "", title = "Herbivores") +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = mean.herb.b, linetype = "dashed", color = "grey", size = 1) +
  coord_flip()

b.sites.herb.df <- MLEbins.site.herb %>% dplyr::select(site_name, b.bins.herb = b)

# get carnivore biomass from carn.gam.df
herb_bio <- herb.gam.df %>% dplyr::select(site_name, mean_bio_herb = mean_bio_hectare)

gravity.herb.df <- gravity.df %>%
  dplyr::left_join(., b.sites.herb.df, by = "site_name") %>%
  dplyr::left_join(., herb_bio, by = "site_name")

# rerun GAM models

# Models with Gamma family and log link
bin.herb.gam1 <- gam(abs(b.bins.herb) ~ s(mean_bio_herb, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3),
                    data = gravity.herb.df, family = Gamma(link=log), na.action = "na.fail")
dd <- dredge(bin.herb.gam1, rank = "AICc")

# Check residuals to make sure the model meets assumptions
bin.herb.gam2 <- gam(abs(b.bins.herb) ~ s(mean_bio_herb, k=3) + s(mean_complexity, k=3) + s(hard_coral, k=3), 
                    data = gravity.herb.df, family = Gamma(link=log))
summary(bin.herb.gam2)
par(mfrow = c(2,2))
gam.check(bin.herb.gam2)

bin.herb.gam3 <- gam(abs(b.bins.herb) ~ s(mean_bio_herb, k=3), data = gravity.herb.df, 
                    family = Gamma(link=log))
summary(bin.herb.gam3)
par(mfrow = c(2,2))
gam.check(bin.herb.gam3)

bin.herb.gam4 <- gam(abs(b.bins.herb) ~ s(mean_complexity, k=3), data = gravity.herb.df, 
                    family = Gamma(link=log))
summary(bin.herb.gam4)
par(mfrow = c(2,2))
gam.check(bin.herb.gam4)

bin.herb.gam5 <- gam(abs(b.bins.herb) ~ s(hard_coral, k=3), data = gravity.herb.df, 
                    family = Gamma(link=log))
summary(bin.herb.gam5)
par(mfrow = c(2,2))
gam.check(bin.herb.gam5)

# Plot predicted values

newd <- data.frame(mean_bio_herb = seq(min(gravity.herb.df$mean_bio_herb), max(gravity.herb.df$mean_bio_herb), length.out = 2000))
pred <- predict.gam(bin.herb.gam3, newd, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd$b_pred <- exp(pred$fit) * -1
newd$b_pred_upr <- exp(upr) * -1
newd$b_pred_lwr <- exp(lwr) * -1
tmp_grav <- gravity.herb.df
gravity.herb.df$region <- factor(gravity.herb.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = mean_bio_herb, y = b_pred), size = 1.5, data = newd) +
  geom_ribbon(aes(x = mean_bio_herb, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd) +
  geom_point(aes(x = mean_bio_herb, y = b.bins.herb, color = region), data = gravity.herb.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Size spectra (b)") +
  xlab("Biomass (kg/ha)")

newd <- data.frame(mean_complexity = seq(min(gravity.herb.df$mean_complexity), max(gravity.herb.df$mean_complexity), length.out = 2000))
pred <- predict.gam(bin.herb.gam4, newd, se.fit = TRUE)
upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
newd$b_pred <- exp(pred$fit) * -1
newd$b_pred_upr <- exp(upr) * -1
newd$b_pred_lwr <- exp(lwr) * -1
gravity.herb.df$region <- factor(gravity.herb.df$region, levels = c("raja_ampat", "wakatobi", "lombok", 
                                        labels = c("Raja Ampat", "Wakatobi", "Lombok")))
ggplot() +
  geom_line(aes(x = mean_complexity, y = b_pred), size = 1.5, data = newd) +
  geom_ribbon(aes(x = mean_complexity, ymin = b_pred_lwr, ymax = b_pred_upr), alpha = 0.5, fill = "gray", data = newd) +
  geom_point(aes(x = mean_complexity, y = b.bins.herb, color = region), data = gravity.herb.df, size = 2) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_color_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  ylab("Size spectra (b)") +
  xlab("Structural complexity")



