# data processing ---------------------------------------------------------

# set working directory accordingly
setwd("/home/bruno/Nextcloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population.csv",skip=6)

# convert krill structural length to real length in mm
df$l_k_mean <- df$l_k_mean / 0.02

# change steps to factor
df$X.step. <- as.factor(df$X.step.)

# rename levels to match year number
levels(df$X.step.) <- (c(1:nlevels(df$X.step.)))

# convert to numeric
df$X.step. <- as.numeric(levels(df$X.step.))[df$X.step.]

# calculate year number
df$year <- floor(df$X.step./365)

# calculate mean of repetitions and years
df <- aggregate(cbind(n_k_now,l_k_mean,n_k_eggs,chla_max)~species+chla_supply+year,df,mean)

# safe data frame 
write.csv(df,"output/PEKRIS-population-trimmed.csv",row.names=F)

# delete everything
rm(list=ls())

# time series ---------------------------------------

# set working directory accordingly
setwd("/home/bruno/Nextcloud/PEKRIS/013-experiments/")

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
ddf <- gather(df,key="measure",value="value",c("n_k_now","l_k_mean","n_k_eggs"))

# add maximum of each output as own column
ddf$maximum <- c(rep(max(df$n_k_now),nrow(ddf)/3),rep(max(df$l_k_mean),nrow(ddf)/3),rep(max(df$n_k_eggs),nrow(ddf)/3))

# calculate max chl a densities scaled to the maximum of each output
ddf$chla_max <- ddf$chla_max * ddf$maximum / max(df$chla_max)

# rename parameter values for plotting
ddf[ddf$chla_supply=="Const",2] <- "constant max chlorophyll a density"
ddf[ddf$chla_supply=="Lognorm",2] <- "varying max chlorophyll a density"
ddf[ddf$measure=="n_k_now",5] <- "abundance [n] krill"
ddf[ddf$measure=="l_k_mean",5] <- "mean length [mm] krill"
ddf[ddf$measure=="n_k_eggs",5] <- "sum of eggs layed [n]"
ddf[ddf$species=="both",1] <- "salps & krill"
ddf[ddf$species=="krill",1] <- "krill only"

library("scales")
# plot all together
ggplot(ddf,aes(x=year,value)) +
  geom_ribbon(aes(x=year,ymax=chla_max,ymin=0,alpha=0.5),col="green4",fill="green4") +
  geom_line(aes(linetype=species)) +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous("model output",
                     labels=comma,
                     sec.axis=sec_axis(~./max(.)*max(df$chla_max),
                     name="max chl a density [mg / m3]")) +
  facet_grid(measure~chla_supply,scales="free_y",switch="y") +
  labs(x="year") +
  scale_alpha(labels="chlorophyll a density") +
  guides(alpha=guide_legend(title=NULL),
         linetype=guide_legend(title=NULL))

# export this plot as PDF file
ggsave("figures/experiments-both.pdf",width=7,height=8)  

# delete everything
rm(list=ls())
dev.off()

# test for differences ----------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# test for normal distribution
aggregate(cbind(n_k_now,l_k_mean,n_k_eggs)~species+chla_supply,
          df,
          function(x) shapiro.test(x)$p.value)

# test of significant differences for mean length
anova(lm(l_k_mean~species*chla_supply,df))

# test of significant differences for abundance
df$rk1 <- rank(df$n_k_now) # rank transformation
dd <- anova(lm(rk1~species*chla_supply,df)) # anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # sum of Sum Squares
dd$`F value` = dd$`Sum Sq` / ddsum # divide sum square by proportions
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # calculate p-values
dd

# test of significant differences for egg amount
df$rk2 <- rank(df$n_k_eggs) # rank transformation
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
             legend.position=c(0.8,0.2))

library("tidyr") # load tidyr
library("scales") # load scales

# change table from wide to long format
ddf <- gather(df,key="measure",value="value",c("n_k_now","l_k_mean","n_k_eggs"))

ddf$measure <- c(rep("abundance krill [n]",nrow(ddf)/3),
                 rep("mean length of krill [mm]",nrow(ddf)/3),
                 rep("sum of eggs layed [n]",nrow(ddf)/3))

ggplot(ddf,aes(x=chla_supply,y=value,fill=species)) +
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
ggsave("figures/experiments-differences.pdf",width=7,height=6)

# delete everything
rm(list=ls())
dev.off()
