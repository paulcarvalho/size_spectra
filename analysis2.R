#------  initial
load(file = 'covariates.df.Rdata')
library(mgcv)
library(MuMIn)
library(tidyverse)
library(dplyr)
#devtools::install_github("m-clark/visibly")
library(visibly)
#------
# Remove structural complexity
gam2 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + 
              s(hard_coral, k=3) + 
              s(algae, k=3) +
              region,
          data = covariates.df, family = Gamma(link=log), na.action = "na.fail", weights = covariates.df$b.weight)

plot(gam2,pages=1)
summary(gam2)
concurvity(gam2, full = FALSE)
dredge(gam2, rank = "AICc")

# Remove hard coral
gam3 <- gam(abs(b) ~ s(mean_bio_hectare, k=3) + s(algae, k=3) + s(mean_complexity, k=3) + region,
          data = covariates.df, family = Gamma(link=log), na.action = "na.fail", weights = covariates.df$b.weight)
summary(gam3)
concurvity(gam3, full = FALSE)
dredge(gam3, rank = "AICc")

#-----
# Plot sum of AICc weights
dd.gam2 <- dredge(gam2)
dd.gam2.dff <- data.frame(Biomass = dd.gam2$`s(mean_bio_hectare, k = 3)`,
                         Algae = dd.gam2$`s(algae, k = 3)`,
                         Hardcoral = dd.gam2$`s(hard_coral, k = 3)`,
                         weight = dd.gam2$weight)
dd.gam2.dff
dd.gam2.df <- dd.gam2.dff %>% 
  gather(., "covariate", "included", -weight) %>%
  na.omit() %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarize(sum_weight = sum(weight)) %>%
  mutate(model = "Hard coral cover")
dd.gam2.df

#-----
dd.gam3 <- dredge(gam3)
dd.gam3.dff <- data.frame(Biomass = dd.gam3$`s(mean_bio_hectare, k = 3)`,
                         Algae = dd.gam3$`s(algae, k = 3)`,
                         Complexity = dd.gam3$`s(mean_complexity, k = 3)`,
                         weight = dd.gam3$weight)
dd.gam3.df <- dd.gam3.dff %>% 
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
plot_base 
#-----  plot partial effect
df=covariates.df%>% dplyr::select(mean_bio_hectare,hard_coral,algae,region)
summary(df)
ss=plot_gam(gam2,main_var = 'algae')
ss
summary(ss$data)

newdf=data.frame(mean_bio_hectare=646.57,hard_coral=32.371,algae=seq(0,28.2500,length.out = 100),region='raja_ampat')
newdf
predict(gam2,newdata=newdf,type='response')
ss$data$fit
#------

line_color = '#7B321C'
ribbon_color = '#28688640'
model=gam2
df=covariates.df%>% dplyr::select(mean_bio_hectare,hard_coral,algae,region)
newdf<- data.frame(mean_bio_hectare=646.57,
                   hard_coral=32.371,
                   algae=3.3713,
                   region=levels(df$region))

data_list<- newdf%>%bind_cols(
    tibble::as_tibble(
      predict(model, ., se=TRUE))) %>%
  mutate(ll = model$family$linkinv(fit - 1.96*se.fit),
         ul = model$family$linkinv(fit + 1.96*se.fit),
         fit = model$family$linkinv(fit)) %>%
  select(region, fit, ll, ul) %>%
  rename(value = region) 
data_list

g<-data_list %>%
  ggplot(aes(x=value, y=fit,group=1)) +
  geom_errorbar(aes(ymin=ll, ymax=ul),width = 0.1) +
  geom_line(color=line_color) +theme_clean()
g
ggsave('Region_Effect.pdf')
