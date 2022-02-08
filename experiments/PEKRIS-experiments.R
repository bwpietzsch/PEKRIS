# data processing ---------------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population.csv",skip=6)

# convert krill structural length to real length in mm
df$mean_length_krill <- df$mean_length_krill / 0.02

# change steps to factor
df$X.step. <- as.factor(df$X.step.)

# rename levels to match year number
levels(df$X.step.) <- (c(1:nlevels(df$X.step.)))

# convert to numeric
df$X.step. <- as.numeric(levels(df$X.step.))[df$X.step.]

# calculate year number
df$year <- floor(df$X.step./365)

# calculate mean of repetitions and years
df2 <- aggregate(cbind(abundance_krill,mean_length_krill,sum_eggs_krill,max_chla_density,median_abundace_salp_overall)~species+chla_supply+year,df,mean)

# calculate max of salp abundance
df3 <- aggregate(max_abundance_salp_season~species+chla_supply+year,df,max)

# combine data frames
df2$max_abundance_salp_season <- df3$max_abundance_salp_season

# remove first 12 years (transient phase)
df2 <- df2[df2$year>12,]

# safe data frame 
write.csv(df2,"output/PEKRIS-population-trimmed.csv",row.names=F)

# delete everything
rm(list=ls())

# salps: time series ------------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import trimmed data frame
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# load library for plotting
library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position="bottom")

# load tidyr
library("tidyr")

# change table from wide to long format
ddf <- gather(df,key="measure",value="value",c("max_abundance_salp_season","median_abundance_salp_overall"))

# add maximum of each output as own column
ddf$maximum <- c(rep(max(df$max_abundance_salp_season),nrow(ddf)/2),rep(max(df$median_abundance_salp_overall),nrow(ddf)/2))

# calculate max chl a densities scaled to the maximum of each output
ddf$max_chla_density <- ddf$max_chla_density * ddf$maximum / max(df$max_chla_density)

# rename parameter values for plotting
ddf[ddf$chla_supply=="Const",2] <- "constant max chlorophyll a density"
ddf[ddf$chla_supply=="Lognorm",2] <- "varying max chlorophyll a density"
ddf[ddf$measure=="max_abundance_salp_season",8] <- "max peak abundance [n] salp"
ddf[ddf$measure=="median_abundance_salp_overall",8] <- "median abundance [n] salp"
ddf[ddf$species=="both",1] <- "salps & krill"
ddf[ddf$species=="krill",1] <- "krill only"

library("scales")

# plot all together
ggplot(ddf[ddf$species=="salps & krill",],aes(x=year,value)) +
  geom_ribbon(aes(x=year,ymax=max_chla_density,ymin=0,alpha=0.5),col="green4",fill="green4") +
  geom_line() +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous("model output",
                     labels=comma,
                     sec.axis=sec_axis(~./max(.)*max(df$max_chla_density),
                                       name="max chl a density [mg / m3]")) +
  facet_grid(measure~chla_supply,scales="free_y",switch="y") +
  labs(x="year") +
  scale_alpha(labels="chlorophyll a density") +
  guides(alpha=guide_legend(title=NULL),
         linetype=guide_legend(title=NULL))

# export this plot as PDF file
ggsave("figures/experiments-both-salps.pdf",width=6,height=5)  

# delete obsolete data
rm(list=ls())
dev.off()

# salps: test for differences ----------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# remove obsolete data
df <- df[df$species=="both",]

# test for normal distribution of data
aggregate(cbind(max_abundance_salp_season,median_abundance_salp_overall)~chla_supply,
          df,
          function(x) shapiro.test(x)$p.value)

# test of significant differences for salp max abundance
df$rk1 <- rank(df$max_abundance_salp_season) # rank transformation
dd <- anova(lm(rk1~chla_supply,df)) # anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # sum of Sum Squares
dd$`F value` = dd$`Sum Sq` / ddsum # divide sum square by proportions
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # calculate p-values
dd

# test of significant differences for salp max abundance
df$rk1 <- rank(df$median_abundance_salp_overall) # rank transformation
dd <- anova(lm(rk1~chla_supply,df)) # anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # sum of Sum Squares
dd$`F value` = dd$`Sum Sq` / ddsum # divide sum square by proportions
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # calculate p-values
dd

# delete obsolete data
rm(dd,ddsum)
df <- df[,-c(4,5,6,7,10)]

# rename levels of chl a supply
df$chla_supply <- as.factor(df$chla_supply)
levels(df$chla_supply) <- c("Constant Chl a","Lognorm Chl a")

# load ggplot
library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position=c(0.8,0.2))

library("tidyr") # load tidyr
library("scales") # load scales

# change table from wide to long format
ddf <- gather(df,key="measure",value="value",c("max_abundance_salp_season","median_abundance_salp_overall"))

ddf$measure <- c(rep("max abundance salp [n]",nrow(ddf)/2),
                 rep("med abundance salp [mm]",nrow(ddf)/2))

ggplot(ddf[ddf$value<8000,],aes(x=chla_supply,y=value)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E69F00","#56B4E9"),
                    labels=c("yes","no"),
                    name="salps present?") +
  scale_y_continuous(labels=comma) +
  facet_wrap(measure~.,
             scales="free_y",
             ncol=2) +
  labs(x="",y="") +
  theme(legend.text=element_text(size=12))

# export this plot as PDF file
ggsave("figures/experiments-differences-salp.pdf",width=6,height=4)

# delete everything
rm(list=ls())
dev.off()

# krill: time series ---------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import trimmed data frame
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# load library for plotting
library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position="bottom")

# load tidyr
library("tidyr")

# change table from wide to long format
ddf <- gather(df,key="measure",value="value",c("abundance_krill","mean_length_krill","sum_eggs_krill"))

# add maximum of each output as own column
ddf$maximum <- c(rep(max(df$abundance_krill),nrow(ddf)/3),rep(max(df$mean_length_krill),nrow(ddf)/3),rep(max(df$sum_eggs_krill),nrow(ddf)/3))

# calculate max chla densities scaled to the maximum of each output
ddf$max_chla_density <- ddf$max_chla_density * ddf$maximum / max(df$max_chla_density)

# rename parameter values for plotting
ddf[ddf$chla_supply=="Const",2] <- "constant max chlorophyll a density"
ddf[ddf$chla_supply=="Lognorm",2] <- "varying max chlorophyll a density"
ddf[ddf$measure=="abundance_krill",7] <- "abundance [n] krill"
ddf[ddf$measure=="mean_length_krill",7] <- "mean length [mm] krill"
ddf[ddf$measure=="sum_eggs_krill",7] <- "sum of eggs layed [n]"
ddf[ddf$species=="both",1] <- "salps & krill"
ddf[ddf$species=="krill",1] <- "krill only"

library("scales")
# plot all together
ggplot(ddf,aes(x=year,value)) +
  geom_ribbon(aes(x=year,ymax=max_chla_density,ymin=0,alpha=0.5),col="green4",fill="green4") +
  geom_line(aes(linetype=species)) +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous("model output",
                     labels=comma,
                     sec.axis=sec_axis(~./max(.)*max(df$max_chla_density),
                     name="max chl a density [mg / m3]")) +
  facet_grid(measure~chla_supply,scales="free_y",switch="y") +
  # facet_wrap(chla_supply~measure,scales="free_y",strip.position = "top") +
  labs(x="year") +
  scale_alpha(labels="chlorophyll a density") +
  guides(alpha=guide_legend(title=NULL),
         linetype=guide_legend(title=NULL))

# export this plot as PDF file
ggsave("figures/experiments-both.pdf",width=6,height=8)  

# delete everything
rm(list=ls())
dev.off()

# krill: test for differences ----------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# test for normal distribution of krill data
aggregate(cbind(abundance_krill,mean_length_krill,sum_eggs_krill)~species+chla_supply,
          df,
          function(x) shapiro.test(x)$p.value)

# test of significant differences for mean krill length
anova(lm(mean_length_krill~species*chla_supply,df))

# test of significant differences for krill abundance
df$rk1 <- rank(df$abundance_krill) # rank transformation
dd <- anova(lm(rk1~species*chla_supply,df)) # anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # sum of Sum Squares
dd$`F value` = dd$`Sum Sq` / ddsum # divide sum square by proportions
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # calculate p-values
dd

# test of significant differences for egg amount
df$rk2 <- rank(df$sum_eggs_krill) # rank transformation
dd <- anova(lm(rk2~species*chla_supply,df)) # anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # sum of Sum Squares
dd$`F value` = dd$`Sum Sq` / ddsum # divide sum square by proportions
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # calculate p-values
dd

# delete obsolete data
rm(dd,ddsum)
df <- df[,-c(8,9)]

# rename levels of chl a supply
df$chla_supply <- as.factor(df$chla_supply)
levels(df$chla_supply) <- c("Constant Chl a","Lognorm Chl a")

# load ggplot
library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position="bottom")

library("tidyr") # load tidyr
library("scales") # load scales

# change table from wide to long format
ddf <- gather(df,key="measure",value="value",c("abundance_krill","mean_length_krill","sum_eggs_krill"))

ddf$measure <- c(rep("abundance krill [n]",nrow(ddf)/3),
                 rep("mean length of krill [mm]",nrow(ddf)/3),
                 rep("sum of eggs laid [n]",nrow(ddf)/3))

ggplot(ddf,aes(x=chla_supply,y=value,fill=species)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E69F00","#56B4E9"),
                    labels=c("yes","no"),
                    name="salps present?") +
  scale_y_continuous(labels=comma) +
  facet_wrap(measure~.,
             scales="free_y",
             ncol=3) +
  labs(x="",y="") +
  theme(legend.text=element_text(size=12))

# export this plot as PDF file
# ggsave("figures/experiments-differences.pdf",width=7,height=6)
ggsave("figures/experiments-differences.pdf",width=9,height=4)

# delete everything
rm(list=ls())
dev.off()
