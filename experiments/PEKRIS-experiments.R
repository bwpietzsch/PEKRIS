# data processing ---------------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df1 <- read.csv("output/PEKRIS-population1.csv",skip=6)
df2 <- read.csv("output/PEKRIS-population2.csv",skip=6)
df3 <- read.csv("output/PEKRIS-population3.csv",skip=6)

# combine data frames
df <- rbind(df1,df2,df3)

# delete obsolete data
rm(df1,df2,df3)

# convert krill structural length to real length in mm
df$mean_length_krill <- df$mean_length_krill / 0.2

# change steps to factor
df$X.step. <- as.factor(df$X.step.)

# rename levels to match year number
levels(df$X.step.) <- (c(1:nlevels(df$X.step.)))

# convert to numeric
df$X.step. <- as.numeric(levels(df$X.step.))[df$X.step.]

# retrieve model outputs for salp (calculated only on the 180th day of each year)
df1 <- df[df$X.step. %in% seq(180,220*180,365),]

# calculate year number
df$year <- floor(df$X.step./365) + 1
df1$year <- floor(df1$X.step./365)

# calculate mean of repetitions and years
df2 <- aggregate(cbind(abundance_krill,mean_length_krill,sum_eggs_krill,max_chla_density)~species+chla_supply+year,df,mean)

# calculate mean and max of salp abundances
df1 <- aggregate(cbind(max_abundance_salp_season,median_abundance_salp_overall)~species+chla_supply+year,df1,mean)

# remove incomplete last year
df2 <- df2[df2$year<101,]

# combine data frames inserting NA for incomplete last year
df2$median_abundance_salp_overall <- c(df1$median_abundance_salp_overall[7:606])
df2$max_abundance_salp_season <- c(df1$max_abundance_salp_season[7:606])

# calculate number of produced eggs by krill per year
dl <- split(df2,f=list(df2$species,df2$chla_supply))

for (i in 1:length(dl)) {
  a <- dl[[i]][,6]
  b <- c(0,a)
  a <- c(a,0)
  ab <- a - b
  ab <- ab[1:50]
  dl[[i]][,6] <- ab
}

df2 <- dl[[1]]

for (i in 2:length(dl)) {df2 <- rbind(df2,dl[[i]])}

# safe data frame 
write.csv(df2,"output/PEKRIS-population-trimmed.csv",row.names=F)

# delete everything
rm(list=ls())
dev.off()

# time series ---------------------------------------

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

library("scales")

ddf[ddf$species=="krill"&ddf$measure=="max_abundance_salp_season",c(7,8)] <- NA

ddf$measure <- as.factor(ddf$measure)
levels(ddf$measure) <- c("abundance_krill [n]","max_abundance_salp [n]",
                         "mean_length_krill [mm]","sum_eggs_krill [n]")

ddf$measure <- factor(ddf$measure,levels=c("abundance_krill [n]",
                                           "mean_length_krill [mm]",
                                           "sum_eggs_krill [n]",
                                           "max_abundance_salp [n]"))

ddf$species <- as.factor(ddf$species)
levels(ddf$species) <- c("krill & salps","krill only","salps only")

# plot all together
ggplot(ddf[ddf$year>=30,],aes(x=year,value)) +
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
# ggsave("figures/experiments-both.tiff",width=7,height=10,dpi=800)
ggsave("C:/Users/Bruno/ownCloud/PEKRIS/015-manuscript1/Fig2.tiff",width=16,height=22,units="cm",dpi=800)

# delete everything
rm(list=ls())
dev.off()

# test for differences ----------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# remove transient phase (first 30 years)
df <- df[df$year>=30,]

# split data
df1 <- df[df$species!="salps",]
df2 <- df[df$species!="krill",]

# test for normal distribution of data
aggregate(cbind(abundance_krill,mean_length_krill,sum_eggs_krill)~species+chla_supply,
          df1,
          function(x) shapiro.test(x)$p.value)

# test for normal distribution of data
aggregate(max_abundance_salp_season~species+chla_supply,
          df2,
          function(x) shapiro.test(x)$p.value)

# Scheirer-Ray-Hare test for krill abundance
df1$rk1 <- rank(df1$abundance_krill) # transform data into ranks
dd = anova(lm(rk1~species*chla_supply,df1)) # perform anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # calculate total sq
dd$`F value` = c(dd$`Sum Sq`[c(1:3)] / ddsum,NA) # divide total sq by sum sq
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # chi sq test
dd

# Scheirer-Ray-Hare test for mean length of krill
df1$rk2 <- rank(df1$mean_length_krill) # transform data into ranks
dd = anova(lm(rk2~species*chla_supply,df1)) # perform anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # calculate total sq
dd$`F value` = c(dd$`Sum Sq`[c(1:3)] / ddsum,NA) # divide total sq by sum sq
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # chi sq test
dd

# Scheirer-Ray-Hare test for number of eggs released by krill
df1$rk3 <- rank(df1$sum_eggs_krill) # transform data into ranks
dd = anova(lm(rk3~species*chla_supply,df1)) # perform anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # calculate total sq
dd$`F value` = c(dd$`Sum Sq`[c(1:3)] / ddsum,NA) # divide total sq by sum sq
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # chi sq test
dd

# Scheirer-Ray-Hare test for salp abundance
df2$rk1 <- rank(df2$max_abundance_salp_season) # transform data into ranks
dd = anova(lm(rk1~species*chla_supply,df2)) # perform anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # calculate total sq
dd$`F value` = c(dd$`Sum Sq`[c(1:3)] / ddsum,NA) # divide total sq by sum sq
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # chi sq test
dd

# delete obsolete data
rm(df1,df2,ddsum,dd)

# rename levels of chl a supply
df$chla_supply <- as.factor(df$chla_supply)
levels(df$chla_supply) <- c("Constant","Lognorm")

df1 <- df[df$species!="salps",]
df2 <- df[df$species!="krill",]

# change from wide to long format
ddf1 <- gather(df1,key="measure",value="value",c("abundance_krill","mean_length_krill","sum_eggs_krill"))
ddf2 <- gather(df2,key="measure",value="value","max_abundance_salp_season")

# change species to factor
ddf1$species <- as.factor(ddf1$species)
ddf2$species <- as.factor(ddf2$species)

# change factor levels
levels(ddf1$species) <- c("both","single")
levels(ddf2$species) <- c("both","single")

# combine data frames
ddf <- rbind(ddf1[,-c(5,6)],ddf2[,-c(4:6,8)])

ddf$measure <- as.factor(ddf$measure)

levels(ddf$measure) <- c("abundance krill [n]","abundance salp [n]",
                         "mean length krill [mm]","sum of eggs released [n]")

ddf$measure <- factor(ddf$measure,levels=c("abundance krill [n]","mean length krill [mm]",
                                           "sum of eggs released [n]","abundance salp [n]"))

# load ggplot
library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             # text=element_text(size=9),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position="bottom")

library("tidyr") # load tidyr
library("scales") # load scales

ggplot(ddf,aes(x=chla_supply,y=value,fill=species)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E69F00","#56B4E9"),
                    labels=c("yes","no"),
                    name="competition present?") +
  scale_y_continuous(labels=comma) +
  facet_wrap(measure~.,
             scales="free_y",
             ncol=2) +
  labs(x="chorophyll a",y="") +
  theme(legend.text=element_text(size=12))

# export this plot as PDF file
# ggsave("figures/experiments-differences.tiff",width=9,height=4,dpi=800)
ggsave("C:/Users/Bruno/ownCloud/PEKRIS/015-manuscript1/Fig3.tiff",width=16,height=14,units="cm",dpi=800)

# delete everything
rm(list=ls())
dev.off()

# calculate medians -------------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# remove first 30 years of transient phase
df <- df[df$year>=30,]

# split data frame
dl <- split(df,f=list(df$species,df$chla_supply))

# get overview
names(dl)

# create column with chlorophyll supply scenario
df0 <- data.frame("chla_supply"=c(rep("Const",3),rep("Lognorm",3)))

# create column with species scenario
df0$species <- rep(c("both","krill","salps"),2)

# add median of model outputs
df0$abundance_krill <- unlist(lapply(dl,function(x){median(x[,4])}))
df0$mean_length_krill <- unlist(lapply(dl,function(x){median(x[,5])}))
df0$sum_eggs_krill <- unlist(lapply(dl,function(x){median(x[,6])}))
df0$max_abundance_salp_season <- unlist(lapply(dl,function(x){median(x[,9])}))

# view results
df0

# export results
write.csv2(df0,"C:/Users/bruno/ownCloud/PEKRIS/013-experiments/output/medians.csv",row.names=F)

# delete everything
rm(list=ls())
dev.off()

# plot impact salps ----------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/figures/")

# load packages
library("terra")
library("viridis")

# import PEKRIS snapshot as raster file
# get names of all available files
temp <- list.files(pattern=".asc")

# load all raster files into one list
datalist = lapply(temp, rast) 

# remove file extension
temp <- substr(temp,1,nchar(temp)-4)

# give the data its corresponding names
names(datalist) <- temp

rm(temp) # delete names

# calculate chl a reduction based on max chl a value
for (i in 1:length(datalist)) {
  m1 <- max(values(datalist[[i]]))
  datalist[[i]] <- (datalist[[i]] - m1) / m1 * 100
  rm(m1)
}

# convert to data frame for plotting
datalist <- lapply(datalist,as.data.frame,xy=T)

# create column with salp abundance
for (i in 1:length(datalist)) {
  datalist[[i]][,4] <- as.numeric(names(datalist)[i])
}

# rename columns of data frame 1
colnames(datalist[[1]]) <- c("x","y","chla","abundance")

# rename column names
for (i in 2:length(datalist)) {
  colnames(datalist[[i]]) <- colnames(datalist[[1]])
}

# delete i
rm(i)

# create empty data frame
df <- data.frame(NULL)

# combine all into one data frame
for (i in 1:length(datalist)) {
  if(i==1){df <- datalist[[i]]}
  if(i>=2){df <- rbind(df,datalist[[i]])}
}

# delete obsolete data
rm(datalist)

# create factor
df$abundance <- as.factor(df$abundance)
levels(df$abundance) <- c("1,000 salps","3,000 salps","7,000 salps","20,000 salps")

library("ggplot2")

# set black white theme
theme_set(theme_bw())

windowsFonts()

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             text=element_text(family="sans"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position="bottom")

ggplot(df,aes(x=x,y=y,fill=chla)) +
  geom_raster() +
  facet_wrap(.~abundance,ncol=1) +
  scale_fill_gradient(low="black",high="green3") +
  labs(fill="chl a reduction [%]")

ggsave("C:/Users/Bruno/ownCloud/PEKRIS/015-manuscript1/Fig1.tiff",width=8,height=21,units="cm",dpi=800)

# delete everything
rm(list=ls())
dev.off()
