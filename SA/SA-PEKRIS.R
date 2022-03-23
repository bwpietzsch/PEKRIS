# set working directory accordingly
setwd("C:/Users/bruno/ownCloud/PEKRIS/009-SA")

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

levels(oresult$response) <- c("max_abundance_salp_day","median_abundance_salp_overall","max_length_krill",
                              "max_eggs_krill","age_of_first_reproduction_krill")

levels(oresult$parameter)[5] <- "oozoid_respiration"

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
