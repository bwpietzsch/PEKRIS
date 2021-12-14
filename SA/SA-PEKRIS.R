# variance proportions ----------------------------------------------------

# set working directory accordingly
setwd("/home/bruno/Nextcloud/PEKRIS/009-SA/")

# import raw results
df <- read.csv("3-results/SA-PEKRIS.csv",skip=6)

# remove combinations without egg production
df <- df[df$n_k_max_eggs!=0,]

# test fo normal distributions
lm1 <- lm(n_s_max_total~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)
shapiro.test(resid(lm1))

# ANOVA proportions for n_s_max_total
df$n_s_max_total <- rank(df$n_s_max_total)

lm1 <- lm(n_s_max_total~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

# check variance homogeneity
plot(lm1,which=3)

g <- anova(lm1)
g$prop = round(g$`Sum Sq` / sum(g$`Sum Sq`) * 100,2)
g$parameter <- row.names(g)
row.names(g) <- NULL
g <- g[with(g,order(-g$prop)),]
g <- subset(g,prop>=1)
g <- g[,c(7,6)]
g$response <- "n_s_max_total"
oresult <- g

# trst for normal distribution
lm1 <- lm(n_s_med~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

shapiro.test(resid(lm1))

# ANOVA proportions for n_s_med
df$n_s_med <- rank(df$n_s_med)

lm1 <- lm(n_s_med~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

# check variance homogeneity
plot(lm1,which=3)

g <- anova(lm1)
g$prop = round(g$`Sum Sq` / sum(g$`Sum Sq`) * 100,2)
g$parameter <- row.names(g)
row.names(g) <- NULL
g <- g[with(g,order(-g$prop)),]
g <- subset(g,prop>=1)
g <- g[,c(7,6)]
g$response <- "n_s_med"
oresult <- rbind(oresult,g)

# test for normal distribution
lm1 <- lm(n_k_max_eggs~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

shapiro.test(resid(lm1))

# ANOVA proportions for n_k_max_eggs
df$n_k_max_eggs <- rank(df$n_k_max_eggs)

lm1 <- lm(n_k_max_eggs~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

# check variance homogeneity
plot(lm1,which=3)

g <- anova(lm1)
g$prop = round(g$`Sum Sq` / sum(g$`Sum Sq`) * 100,2)
g$parameter <- row.names(g)
row.names(g) <- NULL
g <- g[with(g,order(-g$prop)),]
g <- subset(g,prop>=1)
g <- g[,c(7,6)]
g$response <- "n_k_max_eggs"
oresult <- rbind(oresult,g)

# trst for normal distribution
lm1 <- lm(l_k_max~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

shapiro.test(resid(lm1))

# ANOVA proportions for l_k_max
df$l_k_max <- rank(df$l_k_max)

lm1 <- lm(l_k_max~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

# check variance homogeneity
plot(lm1,which=3)

g <- anova(lm1)
g$prop = round(g$`Sum Sq` / sum(g$`Sum Sq`) * 100,2)
g$parameter <- row.names(g)
row.names(g) <- NULL
g <- g[with(g,order(-g$prop)),]
g <- subset(g,prop>=1)
g <- g[,c(7,6)]
g$response <- "l_k_max"
oresult <- rbind(oresult,g)

# test for normal distribution
lm1 <- lm(d_k_firstrepro~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

shapiro.test(resid(lm1))

# ANOVA proportions for d_k_firstrepro
df$d_k_firstrepro <- rank(df$d_k_firstrepro)

lm1 <- lm(d_k_firstrepro~(vegetation_delay+oozoid_resp+krill_amount+chla_growth+salp_immiprob+salp_mortality+krill_mortality+blasto_resp+salp_amount+chla_decay+salp_length+salp_starvation+krill_hibernation+salp_halfsat+krill_halfsat)^2,df)

# check variance homogeneity
plot(lm1,which=3)

g <- anova(lm1)
g$prop = round(g$`Sum Sq` / sum(g$`Sum Sq`) * 100,2)
g$parameter <- row.names(g)
row.names(g) <- NULL
g <- g[with(g,order(-g$prop)),]
g <- subset(g,prop>=1)
g <- g[,c(7,6)]
g$response <- "d_k_firstrepro"
oresult <- rbind(oresult,g)

oresult <- as.data.frame(oresult)

# round proportion for plotting
oresult$prop <- round(oresult$prop,digits=0)

# remove parameters with less than 5 %
oresult <- oresult[oresult$prop>=5,]

# calculate % sum
ddf <- aggregate(prop~parameter,oresult,sum)

# store order for factor levels of plot
a <- ddf[order(-ddf$prop),1]

oresult <- oresult[oresult$parameter!="Residuals",]

# convert characters to factor
oresult$parameter <- as.factor(oresult$parameter)
oresult$response <- as.factor(oresult$response)

# rearrange levels
oresult$parameter <- factor(oresult$parameter,levels=a)

oresult$parameter <- factor(oresult$parameter,levels=c("chla_growth","chla_decay","krill_halfsat",
                                                       "salp_halfsat","oozoid_resp","salp_mortality",
                                                       "salp_length","Residuals"))

oresult$response <- factor(oresult$response,levels=c("n_s_max_total","n_s_med","l_k_max","n_k_max_eggs","d_k_firstrepro"))

# rename response variables
# levels(oresult$response) <- c("max abundance Salps [n]","median abundance Salps [n]","max size Krill [mm]",
#                               "max amount eggs per Krill [n]","first reproduction Krill [d]")


# plot data
library("ggplot2")

ggplot(oresult,aes(parameter, response)) +
  geom_tile(aes(fill = prop),col="black",na.rm = T) + 
  geom_text(aes(label = prop),na.rm=T) +
  scale_fill_gradient(low = "#A9E2F3", high = "#3A7EAB",name="variance [%]") +
  labs(x="",y="") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black",angle = 45, hjust = 1),
        axis.text.y = element_text(colour="black"),
        legend.position = "top")

ggsave("4-plots/SA-anova.pdf",width=5.5,height=4)

# delete everything
rm(list=ls())
dev.off()

# main effect plots --------------------------------------------------

# set working directory accordingly
setwd("/home/bruno/Nextcloud/PEKRIS/009-SA/")

# import raw results
df <- read.csv("3-results/SA-PEKRIS.csv",skip=6)

# convert krill length to mm
df$l_k_max <- df$l_k_max / 0.2

# load ggplot library
library("ggplot2")

# define function for multiple plots in ggplot
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"))

# plot main effects for d_k_firstrepro
p1 <- ggplot(df,aes(krill_halfsat,d_k_firstrepro)) +
  geom_point(size=0.5)
  
p2 <- ggplot(df,aes(chla_decay,d_k_firstrepro)) +
  geom_point(size=0.5)

p3 <- ggplot(df,aes(chla_growth,d_k_firstrepro)) +
  geom_point(size=0.5)

pdf("4-plots/SA-ME-firstrepro.pdf",width=9,height=3)
multiplot(p1,p2,p3,cols=3)
dev.off()

# plot main effects for n_k_max_eggs
p4 <- ggplot(df,aes(krill_halfsat,n_k_max_eggs)) +
  geom_point(size=0.5)

p5 <- ggplot(df,aes(chla_decay,n_k_max_eggs)) +
  geom_point(size=0.5)

p6 <- ggplot(df,aes(chla_growth,n_k_max_eggs)) +
  geom_point(size=0.5)

pdf("4-plots/SA-ME-maxeggs.pdf",width = 9,height = 3)
multiplot(p4,p5,p6,cols=3)
dev.off()

# plot main effect for l_k_max
p7 <- ggplot(df,aes(krill_halfsat,l_k_max)) +
  geom_point(size=0.5)

p8 <- ggplot(df,aes(chla_decay,l_k_max)) +
  geom_point(size=0.5)

p9 <- ggplot(df,aes(chla_growth,l_k_max)) +
  geom_point(size=0.5)

pdf("4-plots/SA-ME-lkmax.pdf",width = 9,height = 3)
multiplot(p7,p8,p9,cols=3)
dev.off()

# plot main effect for n_s_max_total
p10 <- ggplot(df,aes(oozoid_resp,n_s_max_total)) +
  geom_point(size=0.5)

p11 <- ggplot(df,aes(salp_halfsat,n_s_max_total)) +
  geom_point(size=0.5)

p12 <- ggplot(df,aes(salp_mortality,n_s_max_total)) +
  geom_point(size=0.5)

pdf("4-plots/SA-ME-nstotal.pdf",width = 9,height = 3)
multiplot(p10,p11,p12,cols=3)
dev.off()

# plot main effect for n_s_med
p13 <- ggplot(df,aes(oozoid_resp,n_s_med)) +
  geom_point(size=0.5)

p14 <- ggplot(df,aes(salp_mortality,n_s_med)) +
  geom_point(size=0.5)

p15 <- ggplot(df,aes(salp_halfsat,n_s_med)) +
  geom_point(size=0.5)

pdf("4-plots/SA-ME-nsmed.pdf",width = 9,height = 3)
multiplot(p13,p14,p15,cols=3)
dev.off()

dev.off()
rm(list=ls())
