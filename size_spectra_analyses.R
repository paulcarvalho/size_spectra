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
devtools::install_github("andrew-edwards/sizeSpectra")
devtools::install_github("vqv/ggbiplot")
devtools::install_github("m-clark/visibly")
library(sizeSpectra)
library(ggbiplot)
library(visibly)
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
ray.PLB = (1 - pPLB(x = rax.PLB, b = PLB.bMLE.ra.b, xmin = min(rax.PLB),
    xmax = max(rax.PLB))) * length(ra.input$biomass)
spectra.text <- as.character(round(PLB.bMLE.ra.b, 2))
rab_plot <- ggplot() +
  geom_point(aes(x = (sort(ra.input$biomass, decreasing=TRUE)), y = (1:length(ra.input$biomass))),
             color = "#E69F00", size = 2, alpha = 0.3) +
  xlab(expression(paste("Body sizes, ", italic("x"), " (kg)"))) +
  ylab(expression(paste("Number of body sizes", " ">=" ", italic("x"), "    "))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000),
                     limits = c(0.25, max(length(wa.input$biomass), length(ra.input$biomass), length(lo.input$biomass)))) +
  scale_x_continuous(trans = 'log10', breaks = c(0, 1, 5, 10),
                     limits = c(min(ra.input$min.biomass,wa.input$min.biomass,lo.input$min.biomass), 
                                max(ra.input$max.biomass,wa.input$max.biomass,lo.input$max.biomass))) +
  geom_line(aes(x = rax.PLB, y = ray.PLB), col = 'black', lwd = 1) +
  annotate("text", x=0.08, y=10, label="  Raja Ampat") +
  annotate("text", x=0.08, y=3, label = expression(paste(italic("b = "), -1.58))) +
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
way.PLB = (1 - pPLB(x = wax.PLB, b = PLB.bMLE.wa.b, xmin = min(wax.PLB),
    xmax = max(wax.PLB))) * length(wa.input$biomass)
spectra.text <- as.character(round(PLB.bMLE.wa.b, 2))
wab_plot <- ggplot() +
  geom_point(aes(x = (sort(wa.input$biomass, decreasing=TRUE)), y = (1:length(wa.input$biomass))), 
             color = "#56B4E9", size = 2, alpha = 0.3) +
  xlab(expression(paste("Body sizes, ", italic("x"), " (kg)"))) +
  ylab(expression(paste("Number of body sizes", " ">=" ", italic("x"), "    "))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000),
                     limits = c(0.25, max(length(wa.input$biomass), length(ra.input$biomass), length(lo.input$biomass)))) +
  scale_x_continuous(trans = 'log10', breaks = c(0, 1, 5, 10),
                     limits = c(min(ra.input$min.biomass,wa.input$min.biomass,lo.input$min.biomass), 
                                max(ra.input$max.biomass,wa.input$max.biomass,lo.input$max.biomass))) +
  geom_line(aes(x = wax.PLB, y = way.PLB), col = 'black', lwd = 1) +
  annotate("text", x=0.08, y=10, label="Wakatobi") +
  annotate("text", x=0.08, y=3, label = expression(paste(italic("  b = "), -1.71))) +
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
loy.PLB = (1 - pPLB(x = lox.PLB, b = PLB.bMLE.lo.b, xmin = min(lox.PLB),
    xmax = max(lox.PLB))) * length(lo.input$biomass)
spectra.text <- as.character(round(PLB.bMLE.lo.b, 2))
lob_plot <- ggplot() +
  geom_point(aes(x = (sort(lo.input$biomass, decreasing=TRUE)), y = (1:length(lo.input$biomass))), 
             color = "#009E73", size = 2, alpha = 0.3) +
  xlab(expression(paste("Body sizes, ", italic("x"), " (kg)"))) +
  ylab(expression(paste("Number of body sizes", " ">=" ", italic("x"), "    "))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000),
                     limits = c(0.25, max(length(wa.input$biomass), length(ra.input$biomass), length(lo.input$biomass)))) +
  scale_x_continuous(trans = 'log10', breaks = c(0, 1, 5, 10),
                     limits = c(min(ra.input$min.biomass,wa.input$min.biomass,lo.input$min.biomass), 
                                max(ra.input$max.biomass,wa.input$max.biomass,lo.input$max.biomass))) +
  geom_line(aes(x = lox.PLB, y = loy.PLB), col = 'black', lwd = 1) +
  annotate("text", x=0.07, y=10, label="Lombok") +
  annotate("text", x=0.08, y=3, label = expression(paste(italic("b = "), -2.06))) +
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
        legend.position = c(0.75,0.85))

# Arrange plots
ggarrange(rab_plot, wab_plot, lob_plot, slopes_ci, ncol=2, nrow=2, labels = c("a","b","c","d"))

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

# Biomass in relation to human population gravity ======================================================
ggplot() +
  geom_point(data=covariates.df, aes(x=log(Grav_tot), y=log(mean_bio_hectare), color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                     labels = c("Raja Ampat","Wakatobi","Lombok")) +
  geom_smooth(data=covariates.df, aes(x=log(Grav_tot), y=log(mean_bio_hectare)), color = "black", method = "lm") +
  labs(x="log(Human population gravity)", y="log(Biomass [kg/ha])") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2,0.2))
  

# slope (b) in relation to human population gravity ====================================================
ggplot() +
  geom_point(data=covariates.df, aes(x=log(Grav_tot), y=-log(abs(b)), color=region)) +
  geom_smooth(data=covariates.df, aes(x=log(Grav_tot), y=-log(abs(b))), color="black", method="lm") +
  theme_classic() +
  labs(x="Human population gravity", y="Size spectra slope (b)") +
  scale_color_discrete(name=NULL, labels=c("Lombok","Raja Ampat","Wakatobi")) +
  # scale_x_continuous(limits = c(0,900), expand = c(0,0)) +
  theme(legend.position = c(0.9,0.9))

lm.1 <- lm(log(abs(covariates.df$b)) ~ log(covariates.df$Grav_tot))
summary(lm.1)
plot(lm.1)
cooksd <- cooks.distance(lm.1)
plot(cooksd, pch="*")
abline(h = 4 * mean(cooksd, na.rm=T), col = "red")
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
min(covariates.df$b)

# remove outliers
tmp.covariates.df <- covariates.df[-c(5,16),] 
ggplot() +
  geom_point(data=tmp.covariates.df, aes(x=log(Grav_tot), y=-log(abs(b)), color=region)) +
  geom_smooth(data=tmp.covariates.df, aes(x=log(Grav_tot), y=-log(abs(b))), color="black", method="lm") +
  theme_classic() +
  labs(x="log(Human population gravity)", y="Size spectra slope (b)") +
  scale_color_discrete(name=NULL, labels=c("Lombok","Raja Ampat","Wakatobi")) +
  # scale_x_continuous(limits = c(0,900), expand = c(0,0)) +
  theme(legend.position = c(0.9,0.9))
lm.2 <- lm(-log(abs(tmp.covariates.df$b)) ~ log(tmp.covariates.df$Grav_tot))
summary(lm.2)

# GAMs =================================================================================================

# Pairs plot for covariates
panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}

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

pairs(covariates.df[,c(12,7,8,9)], 
      labels = c("Biomass", "Hard coral cover", "Algal cover", "Structural complexity"),
      lower.panel = panel.cor,
      upper.panel = panel.smooth,
      pch = 19)

covariates.df$region <- as.factor(covariates.df$region)

# Model 1 ==============================================================================================
# Model with all predictor variables and region as a fixed effect
gam1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) + region,
          data = covariates.df, family = Gamma(link=log), na.action = "na.fail", weights = covariates.df$b.weight)
plot(gam1, page = 1)
summary(gam1)
# Concurvity
concurvity(gam1, full = TRUE)
concurvity(gam1, full = FALSE)

# Model 2 ==============================================================================================
# Remove structural complexity due to covariance with hard coral cover (concurvity > 0.3)
gam2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + region,
          data = covariates.df, family = Gamma(link=log), na.action = "na.fail", weights = covariates.df$b.weight)
summary(gam2)
concurvity(gam2, full = FALSE)

# Model 3 ==============================================================================================
# Remove hard coral due to covariance with structural complexity (concurvity > 0.3)
gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(algae, k=3) + s(mean_complexity, k=3) + region,
          data = covariates.df, family = Gamma(link=log), na.action = "na.fail", weights = covariates.df$b.weight)
summary(gam3)
concurvity(gam3, full = FALSE)

# Plot sum of AICc weights from GAMs ===================================================================
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
  mutate(model = "GAMs with hard coral cover")

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
  mutate(model = "GAMs with structural complexity")

dd.gam.df <- rbind(dd.gam2.df, dd.gam3.df)
dd.gam.df$covariate <- c("Algal cover", "Biomass", "Hard coral cover", "Algal cover", "Biomass", "Structural complexity")
dd.gam.df$facet_fill_color <- c(rep("grey80",3),rep("grey40",3))

# Plot
aicc_sums_plot <- ggplot(data=dd.gam.df, aes(x = covariate, y = sum_weight, fill = model)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.5, lty = "dashed") +
  scale_fill_manual(values = c("grey80", "grey50")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.05)) +
  coord_flip() +
  facet_grid(rows = vars(model), scales = "free_y", space = "free_y") + #
  theme(strip.background = element_rect(fill=c("grey80", "grey50")),
        strip.placement = "outside",
        strip.text.y = element_text(angle = 90, face = "bold")) +
  guides(fill = FALSE) +
  xlab("") +
  ylab("Sum of AICc weights")

g <- ggplot_gtable(ggplot_build(aicc_sums_plot))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("grey80", "grey40")
j2 <- which(grepl('rect', g$grobs[[11]]$grobs[[1]]$childrenOrder))
g$grobs[[11]]$grobs[[1]]$children[[1]]$gp$fill <- "grey50"
grid.draw(g)
require(ggplotify)
aicc_sums_plot <- as_ggplot(g)

# GAM 2 partial effects ================================================================================
covariates.df$region <- factor(covariates.df$region, levels = c("raja_ampat", "wakatobi", "lombok"))

# Biomass
parEff_gam2 <- plot_gam(gam2)
parEff_gam2_bio <- parEff_gam2[1]$data %>% filter(term == "mean_bio_hectare")
gam2_bio <- ggplot() +
  geom_line(data = parEff_gam2_bio, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_gam2_bio, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = covariates.df, aes(x = mean_bio_hectare, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.25)) +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("")

# Hard coral
parEff_gam2_hardcoral <- parEff_gam2[1]$data %>% filter(term == "hard_coral")
gam2_hardcoral <- ggplot() +
  geom_line(data = parEff_gam2_hardcoral, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_gam2_hardcoral, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = covariates.df, aes(x = hard_coral, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Hard coral cover (%)")
  
# Algae
parEff_gam2_algae <- parEff_gam2[1]$data %>% filter(term == "algae")
gam2_algae <- ggplot() +
  geom_line(data = parEff_gam2_algae, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_gam2_algae, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = covariates.df, aes(x = algae, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Algal cover (%)")

# Region
df_gam2 <- covariates.df %>% 
  dplyr::select(mean_bio_hectare,hard_coral,algae,region)
newdf_gam2 <- data.frame(mean_bio_hectare=646.57,
                         hard_coral=32.371,
                         algae=3.3713,
                         region=levels(df_gam2$region))

data_list_gam2 <- newdf_gam2 %>% bind_cols(
  tibble::as_tibble(
  predict(gam2, ., se=TRUE))) %>%
  mutate(ll = gam2$family$linkinv(fit - 1.96*se.fit),
         ul = gam2$family$linkinv(fit + 1.96*se.fit),
         fit = gam2$family$linkinv(fit)) %>%
  dplyr::select(region, fit, ll, ul) %>%
  mutate(region = as.factor(region))
data_list_gam2

data_list_gam2$region <- factor(data_list_gam2$region, levels = c("raja_ampat", "wakatobi", "lombok"))

parEff_gam2_region <- data_list_gam2 %>%
  ggplot(aes(x = region, y = -fit, group = 1)) +
  geom_errorbar(aes(ymin = -ll, ymax = -ul), colour = c("#E69F00", "#56B4E9", "#009E73"), width = 0.1) +
  geom_line(color = "light blue4") +
  scale_x_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("")

ggarrange(gam2_bio, gam2_hardcoral, gam2_algae, parEff_gam2_region, nrow = 2, ncol = 2)

# GAM 3 partial effects ================================================================================
parEff_gam3 <- plot_gam(gam3)

# Biomass
parEff_gam3_bio <- parEff_gam3[1]$data %>% filter(term == "mean_bio_hectare")
gam3_bio <- ggplot() +
  geom_line(data = parEff_gam3_bio, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_gam3_bio, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = covariates.df, aes(x = mean_bio_hectare, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Biomass (kg/ha)")

# Structural complexity
parEff_gam3_comp <- parEff_gam3[1]$data %>% filter(term == "mean_complexity")
gam3_comp <- ggplot() +
  geom_line(data = parEff_gam3_comp, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_gam3_comp, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = covariates.df, aes(x = mean_complexity, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Structural complexity")

# Algae
parEff_gam3_algae <- parEff_gam3[1]$data %>% filter(term == "algae")
gam3_algae <- ggplot() +
  geom_line(data = parEff_gam3_algae, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_gam3_algae, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = covariates.df, aes(x = algae, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Algal cover (%)")

# Region
df_gam3 <- covariates.df %>% 
  dplyr::select(mean_bio_hectare,mean_complexity,algae,region)
newdf_gam3 <- data.frame(mean_bio_hectare=646.57,
                         mean_complexity=2.888,
                         algae=3.3713,
                         region=levels(df_gam3$region))

data_list_gam3 <- newdf_gam3 %>% bind_cols(
  tibble::as_tibble(
  predict(gam3, ., se=TRUE))) %>%
  mutate(ll = gam3$family$linkinv(fit - 1.96*se.fit),
         ul = gam3$family$linkinv(fit + 1.96*se.fit),
         fit = gam3$family$linkinv(fit)) %>%
  dplyr::select(region, fit, ll, ul) %>%
  mutate(region = as.factor(region))
data_list_gam3

data_list_gam3$region <- factor(data_list_gam3$region, levels = c("raja_ampat", "wakatobi", "lombok"))

parEff_gam3_region <- data_list_gam3 %>%
  ggplot(aes(x = region, y = -fit, group = 1)) +
  geom_errorbar(aes(ymin = -ll, ymax = -ul), colour = c("#E69F00", "#56B4E9", "#009E73"), width = 0.1) +
  geom_line(color = "light blue4") +
  scale_x_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("")

ggarrange(gam3_bio, gam3_comp, gam3_algae, parEff_gam3_region, nrow = 2, ncol = 2)

# Plot with AIC sum and partial effects ===================================================================
ggarrange(aicc_sums_plot, ggarrange(gam2_bio, gam3_bio, ncol=1, labels = c("b","c")), ncol = 2, labels = "a")

# Density of small, medium and large fishes at each site ==================================================
test <- fish.df %>%
  mutate(size_cat = ifelse(biomass_kg < 0.2, "small",
                    ifelse(biomass_kg >= 0.2 & biomass_kg < 1.2, "medium", "large"))) %>%
  dplyr::group_by(region, site_name, transect, observer, size_cat) %>%
  dplyr::summarize(biomass_kg = sum(biomass_kg)) %>%
  dplyr::group_by(region, site_name, size_cat) %>%
  dplyr::summarize(bio = mean(biomass_kg),
            stderr = std.error(biomass_kg)) %>%
  mutate(bio_ha = bio * 40) %>%
  mutate(stderr_ha = stderr * 40) %>%
  mutate(region_size = paste(region,size_cat, sep="")) %>%
  ungroup()

a1 <- aov(test$bio_ha ~ test$region_size)
summary(a1)
TukeyHSD(a1)
a2 <- aov(test$bio_ha ~ test$region + test$size_cat)
summary(a2)
TukeyHSD(a2)

fish.sizes.df <- fish.df %>%
  mutate(size_cat = ifelse(biomass_kg < 0.2, "small",
                    ifelse(biomass_kg >= 0.2 & biomass_kg < 1.2, "medium", "large"))) %>%
  dplyr::group_by(region, site_name, transect, observer, size_cat) %>%
  dplyr::summarize(biomass_kg = sum(biomass_kg)) %>%
  dplyr::group_by(region, size_cat) %>%
  dplyr::summarize(bio = mean(biomass_kg),
            stderr = std.error(biomass_kg)) %>%
  mutate(bio_ha = bio * 40) %>%
  mutate(stderr_ha = stderr * 40) %>%
  ungroup()

fish.sizes.df <- fish.sizes.df %>% 
  mutate(region = fct_relevel(region, "raja_ampat", "wakatobi", "lombok")) %>%
  mutate(size_cat = as.factor(size_cat)) %>%
  mutate(size_cat = fct_relevel(size_cat, "small","medium","large"))

ggplot(data = fish.sizes.df, aes(x = size_cat, y = bio_ha, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = bio_ha - stderr_ha, ymax = bio_ha + stderr_ha), stat = "identity",
                width = 0.2, position = position_dodge(0.9)) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = c("Small", "Medium", "Large")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                    labels = c("Raja Ampat", "Wakatobi", "Lombok"))+
  labs(x = "Size category", y = "Biomass (kg/ha)") +
  theme(legend.title = element_blank(),
        legend.position = c(0.15,0.85))
  
# Species abundance for carnivores and herbivores ======================================================
species.carn <- fish.df %>%
  filter(tp == "Carnivore") %>%
  dplyr::group_by(genus_species) %>%
  dplyr::summarise(total = sum(abundance)) %>%
  filter(total >= 10)

ggplot() +
  geom_bar(data = species.carn, aes(x = reorder(genus_species, total), y = total), stat = "identity") +
  coord_flip()

species.herb <- fish.df %>%
  filter(tp == "Herbivore") %>%
  dplyr::group_by(genus_species) %>%
  dplyr::summarise(total = sum(abundance)) %>%
  filter(total >= 10)

ggplot() +
  geom_bar(data = species.herb, aes(x = reorder(genus_species, total), y = total), stat = "identity") +
  coord_flip()

# MLE Carnivore and Herbivore ==========================================================================

# CARNIVORE MLE
carn.input <- set.params(carn.df$biomass_kg)
# MLE CARNIVORE
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.carn <- mle_b(region=NA, x=carn.input$biomass, log_x=carn.input$log.biomass, sum_log_x=carn.input$sum.log.biomass,
                 x_min=carn.input$min.biomass, x_max=carn.input$max.biomass)
PLB.bMLE.carn.b <- PLB.return.carn[[1]] 
PLB.minLL.carn.b <- PLB.return.carn[[2]]
# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.carn.b <- PLB.minLL.carn.b$minimum
x <- carn.input$biomass
carnx.PLB = seq(min(carn.input$biomass), max(carn.input$biomass), length=1000) # x values to plot PLB. Note
                                                                         # that these encompass the data, and are not based
                                                                         # on the binning (in MEE Figure 6 the line starts as
                                                                         # min(x), not the first bin.
carny.PLB = (1 - pPLB(x = carnx.PLB, b = PLB.bMLE.carn.b, xmin = min(carnx.PLB),
    xmax = max(carnx.PLB))) * length(carn.input$biomass)
spectra.text <- as.character(round(PLB.bMLE.carn.b, 2))
# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.carn.b - 0.5, PLB.bMLE.carn.b + 0.5, 0.00001) 
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
  PLB.LLvals[i] = negLL.PLB(bvec[i], x=carn.input$biomass, n=length(carn.input$biomass), xmin=carn.input$min.biomass,
  xmax=carn.input$max.biomass, sumlogx=carn.input$sum.log.biomass)   
}
critVal = PLB.minNegLL.carn.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
carnbIn95 <- c(min(bIn95), max(bIn95))

# HERBIVORE MLE
# Herbivore biomass
herb.input <- set.params(herb.df$biomass_kg)
# MLE HERBIVORE
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.herb <- mle_b(region=NA, x=herb.input$biomass, log_x=herb.input$log.biomass, sum_log_x=herb.input$sum.log.biomass,
                 x_min=herb.input$min.biomass, x_max=herb.input$max.biomass)
PLB.bMLE.herb.b <- PLB.return.herb[[1]] 
PLB.minLL.herb.b <- PLB.return.herb[[2]]
# plot and find 95% confidence intervals for MLE method.
PLB.minNegLL.herb.b <- PLB.minLL.herb.b$minimum
x <- herb.input$biomass
herbx.PLB = seq(min(herb.input$biomass), max(herb.input$biomass), length=1000) # x values to plot PLB. Note
                                                                         # that these encompass the data, and are not based
                                                                         # on the binning (in MEE Figure 6 the line starts as
                                                                         # min(x), not the first bin.
herby.PLB = (1 - pPLB(x = herbx.PLB, b = PLB.bMLE.herb.b, xmin = min(herbx.PLB),
    xmax = max(herbx.PLB))) * length(herb.input$biomass)
spectra.text <- as.character(round(PLB.bMLE.herb.b, 2))
# Values of b to test to obtain confidence interval. For the real movement data
# sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# symmetric interval here.
bvec = seq(PLB.bMLE.herb.b - 0.5, PLB.bMLE.herb.b + 0.5, 0.00001) 
PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
for(i in 1:length(bvec)){
  PLB.LLvals[i] = negLL.PLB(bvec[i], x=herb.input$biomass, n=length(herb.input$biomass), xmin=herb.input$min.biomass,
  xmax=herb.input$max.biomass, sumlogx=herb.input$sum.log.biomass)   
}
critVal = PLB.minNegLL.herb.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
bIn95 = bvec[ PLB.LLvals < critVal ]
herbbIn95 <- c(min(bIn95), max(bIn95))

# Plot Carn and Herb rank frequency plots
trophb_plot <- ggplot() +
  geom_point(aes(x = (sort(carn.input$biomass, decreasing=TRUE)), y = (1:length(carn.input$biomass))), 
             color = "#D55E00", size = 2, alpha = 0.25) +
  geom_point(aes(x = (sort(herb.input$biomass, decreasing=TRUE)), y = (1:length(herb.input$biomass))), 
             color = "#0072B2", size = 2, alpha = 0.25) +
  xlab(expression(paste("Body sizes, ", italic("x"), " (kg)"))) +
  ylab(expression(paste("Number of body sizes", " ">=" ", italic("x")))) +
  scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000),
                     limits = c(0.25, max(length(carn.input$biomass),length(herb.input$biomass)))) +
  scale_x_continuous(trans = 'log10', breaks = c(0, 1, 5, 10),
                     limits = c(min(carn.input$min.biomass, herb.input$min.biomass),
                                max(carn.input$max.biomass, herb.input$max.biomass))) +
  geom_line(aes(x = carnx.PLB, y = carny.PLB), col = "#D55E00", lwd = 1) +
  geom_line(aes(x = herbx.PLB, y = herby.PLB), col = "#0072B2", lwd = 1) +
  annotate("text", x=0.1, y=10, label = expression(italic("b")[Carnivore]*" = "* -1.97), color = "#D55E00") +
  annotate("text", x=0.1, y=3, label = expression(italic("b")[Herbivore]*" = "* -1.54), color = "#0072B2") +
  theme_classic()
trophb_plot

# Plot estimated size spectra slopes with confidence intervals (Carn/Herb) =============================
slopesTroph_ci_df <- data.frame(tp = c("Carnivore", "Herbivore"),
                           b = c(PLB.bMLE.carn.b, PLB.bMLE.herb.b),
                           b_lo = c(carnbIn95[1], herbbIn95[1]),
                           b_up = c(carnbIn95[2], herbbIn95[2]))
slopesTroph_ci_df$tp <- as.factor(slopesTroph_ci_df$tp)
slopesTroph_ci_df$tp <- factor(slopesTroph_ci_df$tp, levels = c("Carnivore", "Herbivore"))
slopesTroph_ci <- ggplot() +
  geom_point(data = slopesTroph_ci_df, aes(x = c(1,1), y = b, color = tp)) +
  geom_errorbar(data = slopesTroph_ci_df, aes(x = c(1,1), ymin = b_lo, ymax = b_up, color = tp), width = 0.05) +
  ylab(expression(italic("b"))) +
  xlab("") +
  scale_x_continuous(limits = c(0.9,1.5)) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
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
        legend.position = c(0.55,0.95))

ggarrange(trophb_plot, slopesTroph_ci, ncol=2, nrow=1, labels = c("a", "b"))

# Trophic position analysis for each region =========================================================
# Create separate dataframes
ra.carn <- carn.df %>% filter(region == "raja_ampat")
ra.herb <- herb.df %>% filter(region == "raja_ampat")
wa.carn <- carn.df %>% filter(region == "wakatobi")
wa.herb <- herb.df %>% filter(region == "wakatobi")
lo.carn <- carn.df %>% filter(region == "lombok")
lo.herb <- herb.df %>% filter(region == "lombok")
        
# Set MLE parameters

# Raja Ampat 
ra.carn.input <- set.params(ra.carn$biomass_kg)
ra.herb.input <- set.params(ra.herb$biomass_kg)

# Wakatobi
wa.carn.input <- set.params(wa.carn$biomass_kg)
wa.herb.input <- set.params(wa.herb$biomass_kg)

# Lombok
lo.carn.input <- set.params(lo.carn$biomass_kg)
lo.herb.input <- set.params(lo.herb$biomass_kg)

# MLE CARNIVORE - Raja Ampat
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.ra.carn <- mle_b(region=NA, x=ra.carn.input$biomass, log_x=ra.carn.input$log.biomass, sum_log_x=ra.carn.input$sum.log.biomass,
                 x_min=ra.carn.input$min.biomass, x_max=ra.carn.input$max.biomass)
PLB.bMLE.ra.carn.b <- PLB.return.ra.carn[[1]] 
PLB.minLL.ra.carn.b <- PLB.return.ra.carn[[2]]
ra.carnbIn95 <- slope.conf.int(PLB.bMLE.ra.carn.b, PLB.minLL.ra.carn.b$minimum, ra.carn.input)
  
# MLE HERBIVORE - Raja Ampat
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.ra.herb <- mle_b(region=NA, x=ra.herb.input$biomass, log_x=ra.herb.input$log.biomass, sum_log_x=ra.herb.input$sum.log.biomass,
                 x_min=ra.herb.input$min.biomass, x_max=ra.herb.input$max.biomass)
PLB.bMLE.ra.herb.b <- PLB.return.ra.herb[[1]] 
PLB.minLL.ra.herb.b <- PLB.return.ra.herb[[2]]
ra.herbbIn95 <- slope.conf.int(PLB.bMLE.ra.herb.b, PLB.minLL.ra.herb.b$minimum, ra.herb.input)

# MLE CARNIVORE - Wakatobi
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.wa.carn <- mle_b(region=NA, x=wa.carn.input$biomass, log_x=wa.carn.input$log.biomass, sum_log_x=wa.carn.input$sum.log.biomass,
                 x_min=wa.carn.input$min.biomass, x_max=wa.carn.input$max.biomass)
PLB.bMLE.wa.carn.b <- PLB.return.wa.carn[[1]] 
PLB.minLL.wa.carn.b <- PLB.return.wa.carn[[2]]
wa.carnbIn95 <- slope.conf.int(PLB.bMLE.wa.carn.b, PLB.minLL.wa.carn.b$minimum, wa.carn.input)

# MLE HERBIVORE - Wakatobi
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.wa.herb <- mle_b(region=NA, x=wa.herb.input$biomass, log_x=wa.herb.input$log.biomass, sum_log_x=wa.herb.input$sum.log.biomass,
                 x_min=wa.herb.input$min.biomass, x_max=wa.herb.input$max.biomass)
PLB.bMLE.wa.herb.b <- PLB.return.wa.herb[[1]] 
PLB.minLL.wa.herb.b <- PLB.return.wa.herb[[2]]
wa.herbbIn95 <- slope.conf.int(PLB.bMLE.wa.herb.b, PLB.minLL.wa.herb.b$minimum, wa.herb.input)

# MLE CARNIVORE - Lombok
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.lo.carn <- mle_b(region=NA, x=lo.carn.input$biomass, log_x=lo.carn.input$log.biomass, sum_log_x=lo.carn.input$sum.log.biomass,
                 x_min=lo.carn.input$min.biomass, x_max=lo.carn.input$max.biomass)
PLB.bMLE.lo.carn.b <- PLB.return.lo.carn[[1]] 
PLB.minLL.lo.carn.b <- PLB.return.lo.carn[[2]]
lo.carnbIn95 <- slope.conf.int(PLB.bMLE.lo.carn.b, PLB.minLL.lo.carn.b$minimum, lo.carn.input)

# MLE HERBIVORE - Lombok
# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# as a starting point for nlm for MLE of b for PLB model.
PLB.return.lo.herb <- mle_b(region=NA, x=lo.herb.input$biomass, log_x=lo.herb.input$log.biomass, sum_log_x=lo.herb.input$sum.log.biomass,
                 x_min=lo.herb.input$min.biomass, x_max=lo.herb.input$max.biomass)
PLB.bMLE.lo.herb.b <- PLB.return.lo.herb[[1]] 
PLB.minLL.lo.herb.b <- PLB.return.lo.herb[[2]]
lo.herbbIn95 <- slope.conf.int(PLB.bMLE.lo.herb.b, PLB.minLL.lo.herb.b$minimum, lo.herb.input)


slopesTrophReg_ci_df <- data.frame(region = c("Raja Ampat","Raja Ampat",
                                              "Wakatobi", "Wakatobi",
                                              "Lombok", "Lombok"),
                                   tp = rep(c("Carnivore", "Herbivore"), 3),
                                   b = c(PLB.bMLE.ra.carn.b, PLB.bMLE.ra.herb.b,
                                         PLB.bMLE.wa.carn.b, PLB.bMLE.wa.herb.b,
                                         PLB.bMLE.lo.carn.b, PLB.bMLE.lo.herb.b),
                                   b_lo = c(ra.carnbIn95[1], ra.herbbIn95[1],
                                            wa.carnbIn95[1], wa.herbbIn95[1],
                                            lo.carnbIn95[1], lo.herbbIn95[1]),
                                   b_up = c(ra.carnbIn95[2], ra.herbbIn95[2],
                                            wa.carnbIn95[2], wa.herbbIn95[2],
                                            lo.carnbIn95[2], lo.herbbIn95[2]))
slopesTrophReg_ci_df$region <- as.factor(slopesTrophReg_ci_df$region)
slopesTrophReg_ci_df$region <- factor(slopesTrophReg_ci_df$region, levels = c("Raja Ampat", "Wakatobi", "Lombok"))
slopesTrophReg_ci <- ggplot() +
  geom_point(data = slopesTrophReg_ci_df, aes(x = region, y = b, color = tp)) +
  geom_errorbar(data = slopesTrophReg_ci_df, aes(x = region, ymin = b_lo, ymax = b_up, color = tp), width = 0.05) +
  ylab(expression(italic("b"))) +
  xlab("") +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.title.y = element_text(angle = 0, vjust = 0.5, size = 14),
        legend.position = c(0.15,0.15))  


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
  tmp.input <- set.params(site.tp.df$biomass_kg)
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
  left_join(., covariates.df, by = "site_name") %>%
  dplyr::select(site_name, tp, b = b.x, b.weight = b.weight.x, mean_bio_hectare = mean_bio_hectare.x, region, Grav_tot, hard_coral, 
         algae, mean_complexity, total_bio = mean_bio_hectare.y)

ggplot() +
  geom_point(data=sites.trophic, aes(x=total_bio, y=b, shape=region, color=tp), size = 2) +
  theme_classic() +
  labs(x="Biomass density (kg/ha)", y="Size spectra slope (b)") +
  scale_shape_discrete(name = NULL, labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  scale_color_manual(name = NULL, values = c("#D55E00", "#0072B2")) +
  scale_x_continuous(limits = c(0,2500), expand = c(0,0))

# Density of small, medium and large herbs and carns fishes at each site =================================

carn.sizes.df <- fish.df %>%
  filter(tp == "Carnivore") %>%
  mutate(size_cat = ifelse(biomass_kg < 0.2, "small",
                    ifelse(biomass_kg >= 0.2 & biomass_kg < 1.2, "medium", "large"))) %>%
  dplyr::group_by(region, site_name, transect, observer, size_cat) %>%
  dplyr::summarise(abundance = sum(abundance)) %>%
  dplyr::group_by(region, size_cat) %>%
  dplyr::summarise(ab = mean(abundance, na.rm = TRUE),
                   stderr = sd(abundance, na.rm = TRUE)) %>%
  mutate(ab_ha = ab * 40) %>%
  mutate(stderr_ha = stderr * 40) %>%
  mutate(region = as.factor(region)) %>%
  mutate(region = fct_relevel(region, "raja_ampat", "wakatobi", "lombok"))

ggplot(data = carn.sizes.df, aes(x = size_cat, y = ab_ha, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_ha - stderr_ha, ymax = ab_ha + stderr_ha), stat = "identity",
                width = 0.2, position = position_dodge(0.9)) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = c("Large", "Medium", "Small")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                    labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  labs(title = "Carnivore", x = "Size category", y = "Biomass (kg/ha)") +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.85))
  
herb.sizes.df <- fish.df %>%
  filter(tp == "Herbivore") %>%
  mutate(size_cat = ifelse(biomass_kg < 0.2, "small",
                    ifelse(biomass_kg >= 0.2 & biomass_kg < 1.2, "medium", "large"))) %>%
  dplyr::group_by(region, site_name, transect, observer, size_cat) %>%
  dplyr::summarize(biomass_kg = sum(biomass_kg)) %>%
  dplyr::group_by(region, size_cat) %>%
  dplyr::summarize(bio = mean(biomass_kg),
                   stderr = std.error(biomass_kg)) %>%
  mutate(bio_ha = bio * 40) %>%
  mutate(stderr_ha = stderr * 40) %>%
  mutate(region = as.factor(region)) %>%
  mutate(region = fct_relevel(region, "raja_ampat", "wakatobi", "lombok"))

ggplot(data = herb.sizes.df, aes(x = size_cat, y = bio_ha, fill = region)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = bio_ha - stderr_ha, ymax = bio_ha + stderr_ha), stat = "identity",
                width = 0.2, position = position_dodge(0.9)) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = c("Large", "Medium", "Small")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"),
                    labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  labs(title = "Herbivore", x = "Size category", y = "Biomass (kg/ha)") +
  theme(legend.title = element_blank(),
        legend.position = c(0.85,0.85))


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

# Carvnivore GAM 1 =====================================================================================
carn.gam1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) + region,
                 data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = carn.gam.df$b.weight)
summary(carn.gam1)
plot(carn.gam1, pages = 1)
concurvity(carn.gam1, full = TRUE)
concurvity(carn.gam1, full = FALSE)

# Carvnivore GAM 2 =====================================================================================
# Remove structural complexity due to concurvity with hard coral > 0.3
carn.gam2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + region,
                 data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = carn.gam.df$b.weight)
summary(carn.gam2)
plot(carn.gam2, pages = 1)
concurvity(carn.gam2, full = FALSE)

# Carvnivore GAM 3 =====================================================================================
# Remove hard coral due to concurvity with strucural complexity > 0.3
carn.gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(mean_complexity, k=3) + s(algae, k=3) + region,
                 data = carn.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = carn.gam.df$b.weight)
summary(carn.gam3)
plot(carn.gam3, pages = 1)
concurvity(carn.gam3, full = FALSE)

# Herbivore GAM 1 =====================================================================================
herb.gam1 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + s(mean_complexity, k=3) + region,
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = herb.gam.df$b.weight)
summary(herb.gam1)
plot(herb.gam1, pages = 1)
concurvity(herb.gam1, full = TRUE)
concurvity(herb.gam1, full = FALSE)

# Herbivore GAM 2 =====================================================================================
# Remove structural complexity due to concurvity with hard coral > 0.3
herb.gam2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(hard_coral, k=3) + s(algae, k=3) + region,
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = herb.gam.df$b.weight)
summary(herb.gam2)
plot(herb.gam2, pages = 1)
concurvity(herb.gam2, full = TRUE)
concurvity(herb.gam2, full = FALSE)

# Herbivore GAM 3 =====================================================================================
# Remove hard coral due to concurvity with strucural complexity > 0.3
herb.gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(mean_complexity, k=3) + s(algae, k=3) + region,
                 data = herb.gam.df, family = Gamma(link = log), na.action = "na.fail", weights = herb.gam.df$b.weight)
summary(herb.gam3)
plot(herb.gam3, pages = 1)
concurvity(herb.gam3, full = FALSE)

# Plot sum of AICc weights ============================================================================
dd.carn.gam2 <- dredge(carn.gam2)
dd.carn.gam3 <- dredge(carn.gam3)
dd.herb.gam2 <- dredge(herb.gam2)
dd.herb.gam3 <- dredge(herb.gam3)

# Carn gams
dd.carn.gam2.df <- data.frame(Biomass = dd.carn.gam2$`s(mean_bio_hectare, k = 3)`,
                              Algae = dd.carn.gam2$`s(algae, k = 3)`,
                              Hardcoral = dd.carn.gam2$`s(hard_coral, k = 3)`,
                              weight = dd.carn.gam2$weight)
dd.carn.gam2.df <- dd.carn.gam2.df %>% 
  gather(., "covariate", "included", -weight) %>%
  na.omit() %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarize(sum_weight = sum(weight)) %>%
  mutate(tp = "Carnivore") %>%
  mutate(model = "GAMs with hard coral cover")
dd.carn.gam3.df <- data.frame(Biomass = dd.carn.gam3$`s(mean_bio_hectare, k = 3)`,
                              Algae = dd.carn.gam3$`s(algae, k = 3)`,
                              Complexity = dd.carn.gam3$`s(mean_complexity, k = 3)`,
                              weight = dd.carn.gam3$weight)
dd.carn.gam3.df <- dd.carn.gam3.df %>% 
  gather(., "covariate", "included", -weight) %>%
  na.omit() %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarize(sum_weight = sum(weight)) %>%
  mutate(tp = "Carnivore") %>%
  mutate(model = "GAMs with structural complexity")

# Herb gams
dd.herb.gam2.df <- data.frame(Biomass = dd.herb.gam2$`s(mean_bio_hectare, k = 3)`,
                              Algae = dd.herb.gam2$`s(algae, k = 3)`,
                              Hardcoral = dd.herb.gam2$`s(hard_coral, k = 3)`,
                              weight = dd.herb.gam2$weight)
dd.herb.gam2.df <- dd.herb.gam2.df %>% 
  gather(., "covariate", "included", -weight) %>%
  na.omit() %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarize(sum_weight = sum(weight)) %>%
  mutate(tp = "Herbivore") %>%
  mutate(model = "GAMs with hard coral cover")
dd.herb.gam3.df <- data.frame(Biomass = dd.herb.gam3$`s(mean_bio_hectare, k = 3)`,
                              Algae = dd.herb.gam3$`s(algae, k = 3)`,
                              Complexity = dd.herb.gam3$`s(mean_complexity, k = 3)`,
                              weight = dd.herb.gam3$weight)
dd.herb.gam3.df <- dd.herb.gam3.df %>% 
  gather(., "covariate", "included", -weight) %>%
  na.omit() %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarize(sum_weight = sum(weight)) %>%
  mutate(tp = "Herbivore") %>%
  mutate(model = "GAMs with structural complexity")

# Merge data frames
dd.gam.tp.df <- rbind(dd.carn.gam2.df, dd.carn.gam3.df, dd.herb.gam2.df, dd.herb.gam3.df)
dd.gam.tp.df$model <- as.factor(dd.gam.tp.df$model)
dd.gam.tp.df$covariate <- c("Algal cover", "Biomass", "Hard coral cover", 
                            "Algal cover", "Biomass", "Structural complexity",
                            "Algal cover", "Biomass", "Hard coral cover", 
                            "Algal cover", "Biomass", "Structural complexity")

# Plot
aicc_sums_tp_plot <- ggplot(data=dd.gam.tp.df, aes(x = covariate, y = sum_weight, fill = tp)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.5, lty = "dashed") +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  facet_grid(rows = vars(model), scales = "free_y", space = "free_y", switch = "y") +
  theme(strip.background = element_rect(fill="white"),
        strip.placement = "outside",
        legend.title = element_blank(),
        legend.position = c(0.92, 0.95)) +
  xlab("") +
  ylab("Sum of AICc weights") +
  scale_fill_manual(values = c("#D55E00", "#0072B2"))
aicc_sums_tp_plot

# Carn GAM 2 partial effects ================================================================================
# Dummy plot for legend
# dummyplot <- ggplot() +
#   geom_line(data = parEff_carngam2_bio, aes(x = value, y=-fit, color = "#D55E00"), size = 1) +
#   geom_line(data = parEff_carngam2_bio, aes(x = value, y=-fit, color = "#0072B2"), size = 1) +
#   scale_color_manual(values = c("#D55E00","#0072B2"),
#                      labels = c("Carnivore models", "Herbivore models")) +
#   theme_classic()

# Biomass
parEff_carngam2 <- plot_gam(carn.gam2)
parEff_carngam2_bio <- parEff_carngam2[1]$data %>% filter(term == "mean_bio_hectare")
carngam2_bio <- ggplot() +
  geom_point(data = carn.gam.df, aes(x = mean_bio_hectare, y = b, color = region)) +
  geom_line(data = parEff_carngam2_bio, aes(x = value, y=-fit), color = "#D55E00", size = 1) +
  geom_ribbon(data = parEff_carngam2_bio, aes(x = value, ymin = -ll, ymax = -ul), fill = "#D55E00", alpha = 0.2) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.25)) +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Carnivore biomass (kg/ha)")

# Hard coral
parEff_carngam2_hardcoral <- parEff_carngam2[1]$data %>% filter(term == "hard_coral")
carngam2_hardcoral <- ggplot() +
  geom_line(data = parEff_carngam2_hardcoral, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_carngam2_hardcoral, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = carn.gam.df, aes(x = hard_coral, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Hard coral cover (%)")
  
# Algae
parEff_carngam2_algae <- parEff_carngam2[1]$data %>% filter(term == "algae")
carngam2_algae <- ggplot() +
  geom_line(data = parEff_carngam2_algae, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_carngam2_algae, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = carn.gam.df, aes(x = algae, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Algal cover (%)")

# Region
df_carngam2 <- carn.gam.df %>% 
  dplyr::select(mean_bio_hectare,hard_coral,algae,region)
newdf_carngam2 <- data.frame(mean_bio_hectare=205.906,
                             hard_coral=32.371,
                             algae=3.371,
                             region=levels(df_gam2$region))

data_list_carngam2 <- newdf_carngam2 %>% bind_cols(
  tibble::as_tibble(
  predict(carn.gam2, ., se=TRUE))) %>%
  mutate(ll = carn.gam2$family$linkinv(fit - 1.96*se.fit),
         ul = carn.gam2$family$linkinv(fit + 1.96*se.fit),
         fit = carn.gam2$family$linkinv(fit)) %>%
  dplyr::select(region, fit, ll, ul) %>%
  mutate(region = as.factor(region))
data_list_carngam2
data_list_carngam2$region <- factor(data_list_carngam2$region, levels = c("raja_ampat", "wakatobi", "lombok"))
parEff_carngam2_region <- data_list_carngam2 %>%
  ggplot(aes(x = region, y = -fit, group = 1)) +
  geom_errorbar(aes(ymin = -ll, ymax = -ul), colour = c("#E69F00", "#56B4E9", "#009E73"), width = 0.1) +
  geom_line(color = "light blue4") +
  scale_x_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("")

ggarrange(carngam2_bio, carngam2_hardcoral, carngam2_algae, parEff_carngam2_region, nrow = 2, ncol = 2)

# Carn GAM 3 partial effects ================================================================================
# Biomass
parEff_carngam3 <- plot_gam(carn.gam3)
parEff_carngam3_bio <- parEff_carngam3[1]$data %>% filter(term == "mean_bio_hectare")
carngam3_bio <- ggplot() +
  geom_line(data = parEff_carngam3_bio, aes(x = value, y=-fit), color = "#D55E00", size = 1) +
  geom_ribbon(data = parEff_carngam3_bio, aes(x = value, ymin = -ll, ymax = -ul), fill = "#D55E00", alpha = 0.2) +
  geom_point(data = carn.gam.df, aes(x = mean_bio_hectare, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  # ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  ylab("") +
  xlab("Carnivore Biomass (kg/ha)")

# Structural complexity
parEff_carngam3_comp <- parEff_carngam3[1]$data %>% filter(term == "mean_complexity")
carngam3_comp <- ggplot() +
  geom_line(data = parEff_carngam3_comp, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_carngam3_comp, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = carn.gam.df, aes(x = mean_complexity, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Structural complexity")
  
# Algae
parEff_carngam3_algae <- parEff_carngam3[1]$data %>% filter(term == "algae")
carngam3_algae <- ggplot() +
  geom_line(data = parEff_carngam3_algae, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_carngam3_algae, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = carn.gam.df, aes(x = algae, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Algal cover (%)")

# Region
df_carngam3 <- carn.gam.df %>% 
  dplyr::select(mean_bio_hectare,mean_complexity,algae,region)
newdf_carngam3 <- data.frame(mean_bio_hectare=205.906,
                             mean_complexity=2.888,
                             algae=3.371,
                             region=levels(df_gam2$region))

data_list_carngam3 <- newdf_carngam3 %>% bind_cols(
  tibble::as_tibble(
  predict(carn.gam3, ., se=TRUE))) %>%
  mutate(ll = carn.gam3$family$linkinv(fit - 1.96*se.fit),
         ul = carn.gam3$family$linkinv(fit + 1.96*se.fit),
         fit = carn.gam3$family$linkinv(fit)) %>%
  dplyr::select(region, fit, ll, ul) %>%
  mutate(region = as.factor(region))
data_list_carngam3
data_list_carngam3$region <- factor(data_list_carngam3$region, levels = c("raja_ampat", "wakatobi", "lombok"))
parEff_carngam3_region <- data_list_carngam3 %>%
  ggplot(aes(x = region, y = -fit, group = 1)) +
  geom_errorbar(aes(ymin = -ll, ymax = -ul), colour = c("#E69F00", "#56B4E9", "#009E73"), width = 0.1) +
  geom_line(color = "light blue4") +
  scale_x_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("")

ggarrange(carngam3_bio, carngam3_comp, carngam3_algae, parEff_carngam3_region, nrow = 2, ncol = 2)

# Herb GAM 2 partial effects ================================================================================
# Biomass
parEff_herbgam2 <- plot_gam(herb.gam2)
parEff_herbgam2_bio <- parEff_herbgam2[1]$data %>% filter(term == "mean_bio_hectare")
herbgam2_bio <- ggplot() +
  geom_line(data = parEff_herbgam2_bio, aes(x = value, y=-fit), color = "#0072B2", size = 1) +
  geom_ribbon(data = parEff_herbgam2_bio, aes(x = value, ymin = -ll, ymax = -ul), fill = "#0072B2", alpha = 0.3) +
  geom_point(data = herb.gam.df, aes(x = mean_bio_hectare, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Herbivore biomass (kg/ha)")

# Hard coral
parEff_herbgam2_hardcoral <- parEff_herbgam2[1]$data %>% filter(term == "hard_coral")
herbgam2_hardcoral <- ggplot() +
  geom_line(data = parEff_herbgam2_hardcoral, aes(x = value, y=-fit), color = "#0072B2", size = 1) +
  geom_ribbon(data = parEff_herbgam2_hardcoral, aes(x = value, ymin = -ll, ymax = -ul), fill = "#0072B2", alpha = 0.3) +
  geom_point(data = herb.gam.df, aes(x = hard_coral, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Hard coral cover (%)")
  
# Algae
parEff_herbgam2_algae <- parEff_herbgam2[1]$data %>% filter(term == "algae")
herbgam2_algae <- ggplot() +
  geom_line(data = parEff_herbgam2_algae, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_herbgam2_algae, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = herb.gam.df, aes(x = algae, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Algal cover (%)")

# Region
df_herbgam2 <- herb.gam.df %>% 
  dplyr::select(mean_bio_hectare,hard_coral,algae,region)
newdf_herbgam2 <- data.frame(mean_bio_hectare=327.149,
                             hard_coral=32.371,
                             algae=3.371,
                             region=levels(df_gam2$region))

data_list_herbgam2 <- newdf_herbgam2 %>% bind_cols(
  tibble::as_tibble(
  predict(herb.gam2, ., se=TRUE))) %>%
  mutate(ll = herb.gam2$family$linkinv(fit - 1.96*se.fit),
         ul = herb.gam2$family$linkinv(fit + 1.96*se.fit),
         fit = herb.gam2$family$linkinv(fit)) %>%
  dplyr::select(region, fit, ll, ul) %>%
  mutate(region = as.factor(region))
data_list_herbgam2
data_list_herbgam2$region <- factor(data_list_herbgam2$region, levels = c("raja_ampat", "wakatobi", "lombok"))
parEff_herbgam2_region <- data_list_herbgam2 %>%
  ggplot(aes(x = region, y = -fit, group = 1)) +
  geom_errorbar(aes(ymin = -ll, ymax = -ul), colour = c("#E69F00", "#56B4E9", "#009E73"), width = 0.1) +
  geom_line(color = "light blue4") +
  scale_x_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("")

ggarrange(herbgam2_bio, herbgam2_hardcoral, herbgam2_algae, parEff_herbgam2_region, nrow = 2, ncol = 2)

# Herb GAM 3 partial effects ================================================================================
# Biomass
parEff_herbgam3 <- plot_gam(herb.gam3)
parEff_herbgam3_bio <- parEff_herbgam3[1]$data %>% filter(term == "mean_bio_hectare")
herbgam3_bio <- ggplot() +
  geom_line(data = parEff_herbgam3_bio, aes(x = value, y=-fit), color = "#0072B2", size = 1) +
  geom_ribbon(data = parEff_herbgam3_bio, aes(x = value, ymin = -ll, ymax = -ul), fill = "#0072B2", alpha = 0.3) +
  geom_point(data = herb.gam.df, aes(x = mean_bio_hectare, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  # ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  ylab("") + 
  xlab("Herbivore biomass (kg/ha)")

# Structural complexity
parEff_herbgam3_comp <- parEff_herbgam3[1]$data %>% filter(term == "mean_complexity")
herbgam3_comp <- ggplot() +
  geom_line(data = parEff_herbgam3_comp, aes(x = value, y=-fit), color = "#0072B2", size = 1) +
  geom_ribbon(data = parEff_herbgam3_comp, aes(x = value, ymin = -ll, ymax = -ul), fill = "#0072B2", alpha = 0.3) +
  geom_point(data = herb.gam.df, aes(x = mean_complexity, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  # ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  ylab("") +
  xlab("Structural complexity")
  
# Algae
parEff_herbgam3_algae <- parEff_herbgam3[1]$data %>% filter(term == "algae")
herbgam3_algae <- ggplot() +
  geom_line(data = parEff_herbgam3_algae, aes(x = value, y=-fit), size = 1) +
  geom_ribbon(data = parEff_herbgam3_algae, aes(x = value, ymin = -ll, ymax = -ul), fill = "light blue4", alpha = 0.3) +
  geom_point(data = herb.gam.df, aes(x = algae, y = b, color = region)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("Algal cover (%)")

# Region
df_herbgam3 <- herb.gam.df %>% 
  dplyr::select(mean_bio_hectare,mean_complexity,algae,region)
newdf_herbgam3 <- data.frame(mean_bio_hectare=327.149,
                             mean_complexity=2.888,
                             algae=3.371,
                             region=levels(df_gam2$region))

data_list_herbgam3 <- newdf_herbgam3 %>% bind_cols(
  tibble::as_tibble(
  predict(herb.gam3, ., se=TRUE))) %>%
  mutate(ll = herb.gam3$family$linkinv(fit - 1.96*se.fit),
         ul = herb.gam3$family$linkinv(fit + 1.96*se.fit),
         fit = herb.gam3$family$linkinv(fit)) %>%
  dplyr::select(region, fit, ll, ul) %>%
  mutate(region = as.factor(region))
data_list_herbgam3
data_list_herbgam3$region <- factor(data_list_herbgam3$region, levels = c("raja_ampat", "wakatobi", "lombok"))
parEff_herbgam3_region <- data_list_herbgam3 %>%
  ggplot(aes(x = region, y = -fit, group = 1)) +
  geom_errorbar(aes(ymin = -ll, ymax = -ul), colour = c("#E69F00", "#56B4E9", "#009E73"), width = 0.1) +
  geom_line(color = "light blue4") +
  scale_x_discrete(labels = c("Raja Ampat", "Wakatobi", "Lombok")) +
  theme_classic() +
  ylab(expression(paste("Size spectrum slope (", italic(b), ")"))) +
  xlab("")

ggarrange(herbgam3_bio, herbgam3_comp, herbgam3_algae, parEff_herbgam3_region, nrow = 2, ncol = 2)

# Figure 6 =============================================================================================
ggarrange(carngam2_bio, carngam3_bio,
          herbgam2_bio, herbgam3_bio,
          herbgam2_hardcoral, herbgam3_comp,
          ncol = 2, nrow = 3,
          labels = c("a","b","c","d","e","f"))

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



