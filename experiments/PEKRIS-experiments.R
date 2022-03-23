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

# retrieve model outputs for salp (calculated only on the 180th day of each year)
df1 <- df[df$X.step. %in% seq(180,65*180,365),]

# calculate year number
df$year <- floor(df$X.step./365) + 1
df1$year <- floor(df1$X.step./365)

# calculate mean of repetitions and years
df2 <- aggregate(cbind(abundance_krill,mean_length_krill,sum_eggs_krill,max_chla_density)~species+chla_supply+year,df,mean)

# calculate mean of salp abundances
df1 <- aggregate(cbind(max_abundance_salp_season,median_abundance_salp_overall)~species+chla_supply+year,df1,mean)

# remove incomplete 31st year
df2 <- df2[df2$year<31,]

# combine data frames inserting NA for incomplete last year
df2$median_abundance_salp_overall <- c(df1$median_abundance_salp_overall[5:120],NA,NA,NA,NA)
df2$max_abundance_salp_season <- c(df1$max_abundance_salp_season[5:120],NA,NA,NA,NA)

# safe data frame 
write.csv(df2,"output/PEKRIS-population-trimmed.csv",row.names=F)

# delete everything
rm(list=ls())

# time series ---------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import trimmed data frame
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# delete first 12 years as transient phase
df <- df[df$year>12,]

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
# ddf <- gather(df,key="measure",value="value",c("abundance_krill","mean_length_krill","sum_eggs_krill"))
ddf <- gather(df,key="measure",value="value",c("abundance_krill","mean_length_krill","sum_eggs_krill","max_abundance_salp_season"))

# add maximum of each output as own column
# ddf$maximum <- c(rep(max(df$abundance_krill),nrow(ddf)/3),rep(max(df$mean_length_krill),nrow(ddf)/3),rep(max(df$sum_eggs_krill),nrow(ddf)/3))
ddf$maximum <- c(rep(max(df$abundance_krill),nrow(ddf)/4),
                 rep(max(df$mean_length_krill),nrow(ddf)/4),
                 rep(max(df$sum_eggs_krill),nrow(ddf)/4),
                 rep(max(df$max_abundance_salp_season,na.rm=T),nrow(ddf)/4))

# calculate max chl a densities scaled to the maximum of each output
ddf$max_chla_density <- ddf$max_chla_density * ddf$maximum / max(df$max_chla_density)

# rename parameter values for plotting
# ddf[ddf$chla_supply=="Const",2] <- "constant max chlorophyll a density"
# ddf[ddf$chla_supply=="Lognorm",2] <- "varying max chlorophyll a density"
# ddf[ddf$measure=="abundance_krill",7] <- "abundance [n] krill"
# ddf[ddf$measure=="mean_length_krill",7] <- "mean length [mm] krill"
# ddf[ddf$measure=="sum_eggs_krill",7] <- "sum of eggs layed [n]"
# ddf[ddf$species=="both",1] <- "salps & krill"
# ddf[ddf$species=="krill",1] <- "krill only"

library("scales")
str(ddf)
ddf[ddf$species=="krill"&ddf$measure=="max_abundance_salp_season",7] <- NA

ddf$measure <- as.factor(ddf$measure)
ddf$measure <- factor(ddf$measure,levels=c("abundance_krill","mean_length_krill","sum_eggs_krill","max_abundance_salp_season"))
levels(ddf$measure)

# plot all together
ggplot(ddf,aes(x=year,value)) +
  geom_ribbon(aes(x=year,ymax=max_chla_density,ymin=0,alpha=0.5),col="green4",fill="green4") +
  geom_line(aes(linetype=species),na.rm=T) +
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
# ggsave("figures/experiments-both.pdf",width=6,height=8)
ggsave("figures/experiments-both.pdf",width=7,height=10)  

# delete everything
rm(list=ls())
dev.off()

# test for differences ----------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# remove first 12 years as transient phase
df <- df[df$year>12,]

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

# all together -----------------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import trimmed data frame
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# filter runs of interest (after 12 years of transient phase, 
# variable chl a and both species present)
df <- df[df$year>12 & df$chla_supply!="Const" & df$species=="both",]

# load ggplot
library("ggplot2") 

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position="top")

# load tidyr + scales
library("tidyr") 
library("scales") 

# classify chl a and max abundance salps into three levels each
b <- df$max_chla_density
a <- quantile(b,probs=c(1/3,2/3,1))
df$max_chla_density <- sapply(b, function(x){
  ifelse(x<=a[1],a[1],
         ifelse(x<=a[2],a[2],a[3]))
})

b <- df$max_abundance_salp_season
a <- quantile(b,probs=c(1/3,2/3,1),na.rm=T)

df$max_abundance_salp_season <- sapply(b, function(x){
  ifelse(x<=a[1],a[1],
         ifelse(x<=a[2],a[2],a[3]))
})

# round levels
df$max_chla_density <- round(df$max_chla_density,digits=3)
df$max_abundance_salp_season <- round(df$max_abundance_salp_season,digits=0)

# calculate mean of krill outputs for each level combination
d2 <- aggregate(cbind(abundance_krill,mean_length_krill,sum_eggs_krill)~max_chla_density+max_abundance_salp_season,df,mean)

# convert to factors
d2$max_chla_density <- as.factor(d2$max_chla_density)
d2$max_abundance_salp_season <- as.factor(d2$max_abundance_salp_season)

# change table from wide to long format
d2 <- gather(d2,key="measure",value="value",c("abundance_krill","mean_length_krill","sum_eggs_krill"))

# create column with levels low, moderate and high for discrete color levels
b <- d2[d2$measure=="abundance_krill",4]
a <- quantile(b,probs=c(1/3,2/3,1))
b1 <- sapply(b, function(x){
  ifelse(x<=a[1],"low",
         ifelse(x<=a[2],"moderate","high"))
})

b <- d2[d2$measure=="mean_length_krill",4]
a <- quantile(b,probs=c(1/3,2/3,1))
b2 <- sapply(b, function(x){
  ifelse(x<=a[1],"low",
         ifelse(x<=a[2],"moderate","high"))
})

b <- d2[d2$measure=="sum_eggs_krill",4]
a <- quantile(b,probs=c(1/3,2/3,1))
b3 <- sapply(b, function(x){
  ifelse(x<=a[1],"low",
         ifelse(x<=a[2],"moderate","high"))
})

d2$levels <- c(b1,b2,b3)
d2$levels <- as.factor(d2$levels)
d2$levels <- factor(d2$levels,levels=c("low", "moderate", "high"))

# round values
d2[d2$measure=="abundance_krill",4] <- round(d2[d2$measure=="abundance_krill",4],digits=0)
d2[d2$measure=="mean_length_krill",4] <- round(d2[d2$measure=="mean_length_krill",4],digits=1)
d2[d2$measure=="sum_eggs_krill",4] <- round(d2[d2$measure=="sum_eggs_krill",4],digits=0)

# create tile plot
ggplot(d2,aes(x=max_chla_density,y=max_abundance_salp_season)) +
  geom_tile(aes(fill=levels)) +
  geom_text(aes(label=round(value,1),col=levels)) +
  scale_fill_viridis_d(option="D") +
  scale_color_manual(values=c("white","black","black"),name=NULL,labels=NULL,breaks=NULL) +
  facet_wrap(.~measure)

ggsave("figures/experiments-both-togethter.pdf",width=9,height=4)

# delete everything
rm(list=ls())
dev.off()

# OLD: salps: time series ------------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import trimmed data frame
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# remove years without data
df <- df[!is.na(df$max_abundance_salp_season),]

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
ddf$maximum <- c(rep(max(df$max_abundance_salp_season),nrow(ddf)/2),
                 rep(max(df$median_abundance_salp_overall),nrow(ddf)/2))

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

# OLD: krill: time series ---------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import trimmed data frame
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# delete first 12 years as transient phase
df <- df[df$year>12,]

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

# calculate max chl a densities scaled to the maximum of each output
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
  geom_line(aes(linetype=species),na.rm=T) +
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
ggsave("figures/experiments-both.pdf",width=6,height=8)

# delete everything
rm(list=ls())
dev.off()
