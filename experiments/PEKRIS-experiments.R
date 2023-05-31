# data processing ---------------------------------------------------------

# set working directory accordingly
# setwd("C:/Users/Bruno/Nextcloud/PEKRIS/013-experiments/")
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population.csv",skip=6)

# convert krill structural length to real length in mm
df$mean_length_krill <- df$mean_length_krill / 0.2

# calculate densities of krill and salps
df$density_krill <- df$abundance_krill / (51 * 51 * 16) * 1000
df$density_salps <- df$max_abundance_salp_season / (51 * 51 * 16) * 1000

max(df$density_salps)

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

# calculate median of repetitions and years
df2 <- aggregate(cbind(density_krill,mean_length_krill,sum_eggs_krill,max_chla_density)~species+chla_supply+year,df,mean)

# calculate median and max of salp abundances
df1 <- aggregate(cbind(density_salps,median_abundance_salp_overall)~species+chla_supply+year,df1,mean)

# remove incomplete last year
df2 <- df2[df2$year<101,]

# combine data frames inserting NA for incomplete last year
df2$median_abundance_salp_overall <- c(df1$median_abundance_salp_overall[7:606])
df2$density_salps <- c(df1$density_salps[7:606])
# calculate number of produced eggs by krill per year

dl <- split(df2,f=list(df2$species,df2$chla_supply))

for (i in 1:length(dl)) {
  a <- dl[[i]][,6]
  b <- c(0,a)
  a <- c(a,0)
  ab <- a - b
  ab <- ab[1:100]
  dl[[i]][,6] <- ab
}

df2 <- dl[[1]]

for (i in 2:length(dl)) {df2 <- rbind(df2,dl[[i]])}

# remove transient phase (first 30 years)
df2 <- df2[df2$year>=30,]

# calculate krill egg density
df2$sum_eggs_krill <- df2$sum_eggs_krill / (51 * 51 * 16) * 1000

# safe data frame 
write.csv(df2,"output/PEKRIS-population-trimmed.csv",row.names=F)

# delete everything
rm(list=ls())
dev.off()

# time series ---------------------------------------

# set working directory accordingly
# setwd("C:/Users/Bruno/nextcloud/PEKRIS/013-experiments/")
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
             panel.border = element_rect(colour = "black", fill=NA),
             legend.position="bottom")

# load tidyr
library("tidyr")

# change table from wide to long format
ddf <- gather(df,
              key="measure",
              value="value",
              c("density_krill","mean_length_krill","sum_eggs_krill","density_salps"))

# add maximum of each output as own column
ddf$maximum <- c(rep(max(df$density_krill),nrow(ddf)/4),
                 rep(max(df$mean_length_krill),nrow(ddf)/4),
                 rep(max(df$sum_eggs_krill),nrow(ddf)/4),
                 rep(max(df$density_salps),nrow(ddf)/4))

# calculate max chl a densities scaled to the maximum of each output
ddf$max_chla_density <- ddf$max_chla_density * ddf$maximum / max(df$max_chla_density)

library("scales")

ddf[ddf$species=="krill"&ddf$measure=="density_salps",c(7,8)] <- NA

ddf$measure <- as.factor(ddf$measure)
levels(ddf$measure) <- c("krill density [n / 1000 m3]","salp density [n / 1000 m3]",
                         "krill mean length [mm]","krill eggs density [n / 1000 m3]")

ddf$measure <- factor(ddf$measure,levels=c("krill density [n / 1000 m3]",
                                           "krill mean length [mm]",
                                           "krill eggs density [n / 1000 m3]",
                                           "salp density [n / 1000 m3]"))

ddf$species <- as.factor(ddf$species)
levels(ddf$species) <- c("krill & salps","krill only","salps only")

# plot all together
ggplot(ddf,aes(x=year,value)) +
  geom_ribbon(aes(x=year,ymax=max_chla_density,ymin=0,alpha=0.5),col="green4",fill="green4") +
  geom_line(aes(linetype=species),na.rm=T) +
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  scale_y_continuous("model output",
                     labels=comma,
                     sec.axis=sec_axis(~./max(.,na.rm=T)*max(df$max_chla_density),
                                       name=expression(paste("max chl ",italic("a")," concentration [mg / m3]")))) +
  facet_grid(measure~chla_supply,scales="free_y",switch="y") +
  labs(x="year") +
  scale_alpha(labels="chlorophyll a concentration") +
  guides(alpha=guide_legend(title=NULL),
         linetype=guide_legend(title=NULL))

# export this plot as PDF file
# ggsave("figures/experiments-both.tiff",width=7,height=10,dpi=800)
ggsave("C:/Users/Bruno/ownCloud/PEKRIS/015-manuscript1/revision1/Fig3.tiff",width=16,height=22,units="cm",dpi=250)

# delete everything
rm(list=ls())
dev.off()

# test for differences ----------------------------------------------------

# set working directory accordingly
# setwd("C:/Users/Bruno/nextcloud/PEKRIS/013-experiments/")
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# load libraries
library("tidyr")

# import simulation results
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# split data
df1 <- df[df$species!="salps",]
df2 <- df[df$species!="krill",]

# test for normal distribution of data
aggregate(cbind(density_krill,mean_length_krill,sum_eggs_krill)~species+chla_supply,
          df1,
          function(x) shapiro.test(x)$p.value)

# test for normal distribution of data
aggregate(cbind(density_salps)~species+chla_supply,
          df2,
          function(x) shapiro.test(x)$p.value)

# Scheirer-Ray-Hare test for krill density
df1$rk1 <- rank(df1$density_krill) # transform data into ranks
dd = anova(lm(rk1~species*chla_supply,df1)) # perform anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # calculate total sq
dd$`F value` = c(dd$`Sum Sq`[c(1:3)] / ddsum,NA) # divide total sq by sum sq
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # chi sq test
dd

# ANOVA for mean length of krill
dd = anova(lm(mean_length_krill~species*chla_supply,df1)) # perform anova
dd

# Scheirer-Ray-Hare test for number of eggs released by krill
df1$rk3 <- rank(df1$sum_eggs_krill) # transform data into ranks
dd = anova(lm(rk3~species*chla_supply,df1)) # perform anova of ranks
ddsum = sum(dd$`Sum Sq`) / sum(dd$Df) # calculate total sq
dd$`F value` = c(dd$`Sum Sq`[c(1:3)] / ddsum,NA) # divide total sq by sum sq
dd$`Pr(>F)` = 1 - pchisq(dd$`F value`,1) # chi sq test
dd

# Scheirer-Ray-Hare test for salp abundance
df2$rk1 <- rank(df2$density_salps) # transform data into ranks
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
ddf1 <- gather(df1,
               key="measure",
               value="value",
               c("density_krill","mean_length_krill","sum_eggs_krill"))

ddf2 <- gather(df2,
               key="measure",
               value="value",
               c("density_salps"))

# change species to factor
ddf1$species <- as.factor(ddf1$species)
ddf2$species <- as.factor(ddf2$species)

# change factor levels
levels(ddf1$species) <- c("both","single")
levels(ddf2$species) <- c("both","single")

# combine data frames
ddf <- rbind(ddf1[,-c(5,6)],ddf2[,-c(4:6,8)])

ddf$measure <- as.factor(ddf$measure)

levels(ddf$measure) <- c("density krill [n/1000m3]","density salp [n/1000m3]",
                         "mean length krill [mm]","density krill eggs [n/1000m3]")

ddf$measure <- factor(ddf$measure,
                      levels=c("density krill [n/1000m3]","mean length krill [mm]",
                               "density krill eggs [n/1000m3]","density salp [n/1000m3]"))

# load ggplot
library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             # text=element_text(size=9),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA),
             legend.position="bottom")

library("scales") # load scales

ggplot(ddf,aes(x=chla_supply,y=value,fill=species)) +
  geom_boxplot() +
  # stat_summary(fun=mean, geom="point",size=2,shape=4,position = position_dodge2(width=0.75, preserve="single")) +
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
ggsave("C:/Users/Bruno/ownCloud/PEKRIS/015-manuscript1/revision1/Fig4.tiff",width=16,height=14,units="cm",dpi=300)

# delete everything
rm(list=ls())
dev.off()

# calculate medians -------------------------------------------------------

# set working directory accordingly
# setwd("C:/Users/Bruno/nextcloud/PEKRIS/013-experiments/")
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population-trimmed.csv")

# split data frame
dl <- split(df,f=list(df$species,df$chla_supply))

# get overview
names(dl)

# create column with chlorophyll supply scenario
df0 <- data.frame("chla_supply"=c(rep("Const",3),rep("Lognorm",3)))

# create column with species scenario
df0$species <- rep(c("both","krill","salps"),2)

str(dl[[1]])

# add median of model outputs
df0$density_krill <- unlist(lapply(dl,function(x){median(x[,4])}))
df0$mean_length_krill <- unlist(lapply(dl,function(x){median(x[,5])}))
df0$sum_eggs_krill <- unlist(lapply(dl,function(x){median(x[,6])}))
df0$density_salps <- unlist(lapply(dl,function(x){median(x[,9])}))

# view results
df0

# export results
# write.csv2(df0,"C:/Users/bruno/nextcloud/PEKRIS/015-manuscript1/revision1/medians.csv",row.names=F)
write.csv2(df0,"C:/Users/bruno/ownCloud/PEKRIS/015-manuscript1/revision1/medians.csv",row.names=F)

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
a <- substr(temp,1,nchar(temp)-3)
b <- substr(temp,nchar(temp)-1,nchar(temp))
  
rm(temp) # delete names

# give the data its corresponding names
names(datalist) <- a

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

# normalize maximum chl a to 0 % reduction
df$chla <- df$chla + 20

# create factor
df$abundance <- as.factor(df$abundance)

library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             text=element_text(family="sans"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position="right")

levels(df$abundance) <- c(paste("6,735 salps per 1,000 m3 on day",b[4]),
                          paste("10,668 salps per 1,000 m3 on day",b[1]),
                          paste("17,099 salps per 1,000 m3 on day",b[2]),
                          paste("20,015 salps per 1,000 m3 on day",b[3]))

ggplot(df,aes(x=x,y=y,fill=chla)) +
  geom_raster() +
  facet_wrap(.~abundance,ncol=2) +
  scale_fill_gradient(low="black",high="green3") +
  labs(fill="chl a \nreduction\n[%]",
       x="x coordinate",
       y="y coordinate")

ggsave("C:/Users/Bruno/ownCloud/PEKRIS/015-manuscript1/revision1/Fig2.tiff",
       width=16,
       height=12,
       units="cm",
       dpi=800)

# delete everything
rm(list=ls())
dev.off()
