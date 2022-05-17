# data processing ---------------------------------------------------------

# set working directory accordingly
setwd("C:/Users/Bruno/ownCloud/PEKRIS/013-experiments/")

# import simulation results
df <- read.csv("output/PEKRIS-population.csv",skip=6)

# convert krill structural length to real length in mm
df$mean_length_krill <- df$mean_length_krill / 0.2

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

# calculate mean and max of salp abundances
df1 <- aggregate(cbind(max_abundance_salp_season,median_abundance_salp_overall)~species+chla_supply+year,df1,mean)

# remove incomplete 31st year
df2 <- df2[df2$year<31,]

# combine data frames inserting NA for incomplete last year
# df2$median_abundance_salp_overall <- c(df1$median_abundance_salp_overall[5:120],NA,NA,NA,NA)
df2$median_abundance_salp_overall <- c(df1$median_abundance_salp_overall[5:124])
# df2$max_abundance_salp_season <- c(df1$max_abundance_salp_season[5:120],NA,NA,NA,NA)
df2$max_abundance_salp_season <- c(df1$max_abundance_salp_season[5:124])

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

# export for calculation of year 30 differences
# df1 <- df[df$year==30,]
# write.csv2(df1,"figures/year30-raw.csv",row.names = F)

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
ddf[ddf$species=="krill"&ddf$measure=="max_abundance_salp_season",c(7,8)] <- NA

ddf$measure <- as.factor(ddf$measure)
ddf$measure <- factor(ddf$measure,levels=c("abundance_krill","mean_length_krill","sum_eggs_krill","max_abundance_salp_season"))

ddf$species <- as.factor(ddf$species)
levels(ddf$species) <- c("krill & salps","krill only")

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
# ggsave("figures/experiments-both.tiff",width=7,height=10,dpi=800)
ggsave("figures/experiments-both.tiff",width=16,height=22,units="cm",dpi=800)

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
levels(df$chla_supply) <- c("Constant","Lognorm")

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

# change table from wide to long format
ddf <- gather(df,key="measure",value="value",c("abundance_krill","mean_length_krill","sum_eggs_krill"))

ddf$measure <- c(rep("abundance krill [n]",nrow(ddf)/3),
                 rep("mean length of krill [mm]",nrow(ddf)/3),
                 rep("sum of eggs released [n]",nrow(ddf)/3))

ggplot(ddf,aes(x=chla_supply,y=value,fill=species)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#E69F00","#56B4E9"),
                    labels=c("yes","no"),
                    name="salps present?") +
  scale_y_continuous(labels=comma) +
  facet_wrap(measure~.,
             scales="free_y",
             ncol=3) +
  labs(x="chorophyll a",y="") +
  theme(legend.text=element_text(size=12))

# export this plot as PDF file
# ggsave("figures/experiments-differences.tiff",width=9,height=4,dpi=800)
ggsave("figures/experiments-differences.tiff",width=16,height=7.5,units="cm",dpi=800)

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
levels(df$abundance) <- c("1,000 salps","2,000 salps","7,000 salps","15,000 salps","30,000 salps","50,000 salps")

library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6),
             legend.position="bottom")

ggplot(df,aes(x=x,y=y,fill=chla)) +
  geom_raster() +
  facet_wrap(.~abundance,ncol=2) +
  scale_fill_continuous(type="viridis") +
  labs(fill="chl a reduction [%]")

ggsave("salp-impact.tiff",width=16,height=21,units="cm",dpi=800)

# delete everything
rm(list=ls())
dev.off()
