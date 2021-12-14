# krill halfsat -----------------------------------------------------------

# set working directory accordingly
setwd("/home/bruno/Nextcloud/PEKRIS/012-calibration/")

# import results
df <- read.csv("results/CAL-PEKRIS-k-halfsat.csv",skip=6)

# convert length to mm
df$l_k_max <- df$l_k_max / 0.2

# load ggplot library
library("ggplot2")

# define multiplot function
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

# retrieve maximum possible krill_halfsat (above which no eggs are produced)
a <- max(df[df$n_k_max_eggs>0,21])

# plot results
p1 <- ggplot(df,aes(krill_halfsat,d_k_firstrepro)) +
  geom_point() + 
  geom_hline(yintercept = 3*365,col="blue") +
  geom_vline(xintercept = c(0.082,a),col="red")

p2 <- ggplot(df,aes(krill_halfsat,l_k_max)) +
  geom_point() +
  geom_vline(xintercept = c(0.082,a),col="red")

p3 <- ggplot(df,aes(krill_halfsat,n_k_max_eggs)) +
  geom_point() +
  geom_vline(xintercept = c(0.082,a),col="red")

multiplot(p1,p2,p3,cols=3)

pdf("plots/CAL-krill-halfsat.pdf",width=9,height=3)
multiplot(p1,p2,p3,cols=3)
dev.off()
dev.off()

# subset rows with first reproduction during third summer as reported
# by Bahlburg et al. (2021)
df <- df[df$d_k_firstrepro>1000 & df$d_k_firstrepro<1250,]

# plot results
p1 <- ggplot(df,aes(krill_halfsat,d_k_firstrepro)) +
  geom_point(size=0.5) + 
  geom_hline(yintercept = 3*365)

p2 <- ggplot(df,aes(krill_halfsat,l_k_max)) +
  geom_point(size=0.5)

p3 <- ggplot(df,aes(krill_halfsat,n_k_max_eggs)) +
  geom_point(size=0.5)

# plot multiplot
multiplot(p1,p2,p3,cols=3)

# save multiplot as PDF
pdf("plots/CAL-krill-halfsat-reduced.pdf",width=9,height=3)
multiplot(p1,p2,p3,cols=3)
dev.off()

# calculate means for every krill halfsat
df <- aggregate(cbind(l_k_max,n_k_max_eggs,d_k_firstrepro)~krill_halfsat,df,mean)

# extract best fit with Bahlburg et al. (2021)
df[df$l_k_max==max(df$l_k_max),]

# remove and close everything
dev.off()
rm(list=ls())

# blastozoids resp --------------------------------------------------------

# set working directory accordingly
setwd("/home/bruno/Nextcloud/PEKRIS/012-calibration/")

# import results
df <- read.csv("results/CAL-PEKRIS-blasto-resp.csv",skip=6)

# load ggplot library
library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6))#,
             #legend.position="none")

# define function for multi plotting with ggpplot2
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

# remove rows without salps individuals
df <- df[df$max..l..of.blastozoids!="<RuntimePrimitiveException>",]

# format results as numeric
df$max..l..of.blastozoids <- as.numeric(df$max..l..of.blastozoids)
df$min..l..of.blastozoids <- as.numeric(df$min..l..of.blastozoids)

# find maximum length for each blsato_resp value simulated
df <- aggregate(max..l..of.blastozoids~chla_supply+blasto_resp,df,max)

# find lowest blasto_resp where maximum possible length is not exceeded
a <- min(df[df$chla_supply=="Lognorm"&df$max..l..of.blastozoids<=6,2])

# plot result
p1 <- ggplot(df[df$chla_supply=="Lognorm",],aes(y=max..l..of.blastozoids,x=blasto_resp)) +
  geom_point() +
  geom_smooth() +
  geom_hline(yintercept=6,col="red") +
  geom_vline(xintercept=a,col="black") +
  geom_vline(xintercept=7.5,col="green4") +
  labs(y="max l blastozooids",title="Lognorm")

p2 <- ggplot(df[df$chla_supply=="Const",],aes(y=max..l..of.blastozoids,x=blasto_resp)) +
  geom_point() +
  geom_smooth() +
  geom_hline(yintercept=6,col="red") +
  geom_vline(xintercept=a,col="black") +
  geom_vline(xintercept=7.5,col="green4") +
  labs(y="max l of blastozooids",title="Const")

multiplot(p1,p2,cols=2)

pdf("plots/CAL-blasto-resp.pdf",width=8,height=4)
multiplot(p1,p2,cols=2)
dev.off()

rm(list=ls())
dev.off()

# oozoid resp --------------------------------------------------------

# set working directory accordingly
setwd("/home/bruno/Nextcloud/PEKRIS/012-calibration/")

# import results
df <- read.csv("results/CAL-PEKRIS-oozoid-resp.csv",skip=6)

# load ggplot library
library("ggplot2")

# set black white theme
theme_set(theme_bw())

# set black labels
theme_update(axis.text.x = element_text(colour="black"),
             axis.text.y = element_text(colour="black"),
             plot.title = element_text(hjust = 0.5),
             panel.border = element_rect(colour = "black", fill=NA, size=0.6))#,
#legend.position="none")

# define function for multi plotting with ggpplot2
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

# remove rows without salps individuals
df <- df[df$max..l..of.oozoids!="<RuntimePrimitiveException>",]

# format results as numeric
df$max..l..of.oozoids <- as.numeric(df$max..l..of.oozoids)
df$min..l..of.oozoids <- as.numeric(df$min..l..of.oozoids)

# find maximum length for each oozoid_resp value simulated
df <- aggregate(max..l..of.oozoids~chla_supply+oozoid_resp,df,max)

# plot result
p1 <- ggplot(df[df$chla_supply=="Lognorm",],aes(y=max..l..of.oozoids,x=oozoid_resp)) +
  geom_point() +
  geom_smooth() +
  geom_hline(yintercept=16,col="red") +
  labs(y="max l oozoids",title="Lognorm")

p2 <- ggplot(df[df$chla_supply=="Const",],aes(y=max..l..of.oozoids,x=oozoid_resp)) +
  geom_point() +
  geom_smooth() +
  geom_hline(yintercept=16,col="red") +
  labs(y="max l of oozoids",title="Const")

multiplot(p1,p2,cols=2)

pdf("plots/CAL-oozoid-resp.pdf",width=8,height=4)
multiplot(p1,p2,cols=2)
dev.off()

rm(list=ls())
dev.off()
