setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(gridExtra)
library(grid)
library(latex2exp)
library(mgcv)

nice_green <- rgb(0, 200, 50, maxColorValue = 255)

##################### Chapter on MSR

# 1: uniform vs Gaussian MSR, 1000 points
n <- 1e3
maxg_bins <- 200
u <- runif(n)
g <- rnorm(n)
rr_curve_comparison(List(u,g), c("uniform", "gaussian"), add_msr = T, discrete = F, max_gbins = maxg_bins)
##########################################


# 2: random unform samples of several sizes, M vs MSR

#df_approx <- make_random_msr_data_frame2(1000, lw = 1e1, hp = 1e2, verbose = T, perfect_method = F)
#df_perfect <- make_random_msr_data_frame2(1000, lw = 1e1, hp = 1e2, verbose = T, perfect_method = T)
df2 <- readRDS("Scrivania/Tesi/Rexperiments/unstructured_1e4_sample_MSR_positions.rda")

plot1 <- ggplot(df2, aes(x=log(M,10), y=msr)) + 
  geom_point(alpha = 0.2, color='black') +
  #geom_smooth(level = 0.9) +
  geom_smooth(method=lm) +
  xlab(TeX("$\\log_{10}M$")) +
  ylab("MSR") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

plot2 <- ggplot(undersample_data_frame((df2[df2$M>100,]),1000), aes(x=log(M,10), y=msr)) + 
  geom_point(alpha = 0.5) +
  #geom_smooth(level = 0.9) +
  geom_smooth() +
  xlab(TeX("$\\log_{10}M$")) +
  ylab("MSR") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))


#grid.arrange(plot1, plot2, ncol=2)
##########################################



##################### Chapter on applization of MSR to methylation data

# 1: WGBS stomach data
wgbs <- "Scrivania/Tesi/MethylationCode/MethylationData/wgbs/stomach.Rda"
wgbs <- readRDS_if_needed(wgbs)
colnames(wgbs) <- c("chr", "pos", "strand", "reads", "prop")
wgbs$prop <- wgbs$prop/max(wgbs$prop)

n <- 6
l <- length(wgbs$chr)
r1 <- 1:n
r2 <- (l-n+1):l
h <- (wgbs)[r1,]
rownames(h)<-r1
t <- (wgbs)[r2, ]
rownames(t)<-r2
p1<-grid.table(h)
p2<-grid.table(t)
haligned <- gtable_combine(tableGrob(h),tableGrob(t),along=1)
grid.arrange(haligned)
#########
nswgbs <- sum_strands(wgbs)
nswgbs$prop <- round(nswgbs$prop ,2)
n <- 6
l <- length(nswgbs$chr)
r1 <- 1:n
r2 <- (l-n+1):l
h <- (nswgbs)[r1,]
rownames(h)<-r1
t <- (nswgbs)[r2, ]
rownames(t)<-r2
p1<-grid.table(h)
p2<-grid.table(t)
haligned <- gtable_combine(tableGrob(h),tableGrob(t),along=1)
grid.arrange(haligned)
##########

binaryzed <- data.frame(chr=nswgbs$chr,
                        pos=nswgbs$pos,
                        methylated=(nswgbs$prop>=0.5)*1)

n <- 6
l <- length(binaryzed$chr)
r1 <- 1:n
r2 <- (l-n+1):l
h <- (binaryzed)[r1,]
rownames(h)<-r1
t <- (binaryzed)[r2, ]
rownames(t)<-r2
p1<-grid.table(h)
p2<-grid.table(t)
haligned <- gtable_combine(tableGrob(h),tableGrob(t),along=1)
grid.arrange(haligned)

##########################################

# 2: WGBS stomach data summed strands
wgbs <- readRDS_if_needed(wgbs)
wgbs <- sum_strands(wgbs)
##########################################

# 3: WGBS stomach binary
wgbs <- readRDS_if_needed(wgbs)
wgbs <- standard_binaryzer(wgbs)
wgbs$methylation <- wgbs$prop*1
wgbs[,c("chr", "Cpos", "methylation")]
##########################################

# 4: random bernoulli samples of several proportion, p vs MSR
msr_ecdf_ref <- readRDS("Scrivania/Tesi/MethylationCode/MethylationData/msr_ecdf_1e3.Rda")
length <- 1e3
sample_size <- 2e3
df <- make_random_msr_data_frame(sample_size, length = length, lw = 0.02, hp = 0.98)
colnames(df) <- c("x","y")
df = df[complete.cases(df),]
df = df[df$x>=0.02 & df$x<=0.98, ]

confidence <- 0.99
median_function <- extract_ecdf_function(msr_ecdf_ref, 0.5)
up_conf_function <- extract_ecdf_function(msr_ecdf_ref, 1-(1-confidence)/2)
lw_conf_function <- extract_ecdf_function(msr_ecdf_ref, (1-confidence)/2)

px <- seq(from = 0.02, to = 0.98, by = 0.01)
fit = median_function(px)
pred.int <- data.frame(x = px, upr = up_conf_function(px), lwr = lw_conf_function(px))

ggplot(pred.int, aes(x = x, y = fit), alpha = 0.9) +
  #ggtitle("Prediction interval for future observations from predict()") +
  geom_point(data = df, aes(x = x, y = y), alpha = 0.3) +
  geom_smooth(data = pred.int, aes(ymin = lwr, ymax = upr), stat = "identity") +
  xlab("methylation proportion")+ 
  ylab("MSR")+
  labs(colour="")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

##########################################

# 4: random bernoulli samples of several proportion, p vs MSR, 10^4
msr_ecdf_ref <- readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e4.Rda")
length <- 1e4
sample_size <- 2000
df <- make_random_msr_data_frame(sample_size, length = length, lw = 0.02, hp = 0.98)

colnames(df) <- c("x","y")
df = df[complete.cases(df),]
df = df[df$x>=0.02 & df$x<=0.98, ]

confidence <- 0.99
median_function <- extract_ecdf_function(msr_ecdf_ref, 0.5)
up_conf_function <- extract_ecdf_function(msr_ecdf_ref, 1-(1-confidence)/2)
lw_conf_function <- extract_ecdf_function(msr_ecdf_ref, (1-confidence)/2)

px <- seq(from = 0.05, to = 0.95, by = 0.01)
fit = median_function(px)
pred.int <- data.frame(x = px, upr = up_conf_function(px), lwr = lw_conf_function(px))

ggplot(pred.int, aes(x = x, y = fit), alpha = 0.9) +
  #ggtitle("Prediction interval for future observations from predict()") +
  geom_point(data = df, aes(x = x, y = y), alpha = 0.3) +
  geom_smooth(data = pred.int, aes(ymin = lwr, ymax = upr), stat = "identity") +
  xlab("methylation proportion")+ 
  ylab("MSR")+
  labs(colour="")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

##########################################


# 4: random bernoulli samples of several proportion, p vs MSR, 10^5
msr_ecdf_ref <- readRDS("../../MethylationCode/MethylationData/msr_ecdf_1e4.Rda")
length <- 1e5
sample_size <- 500
df <- make_random_msr_data_frame(sample_size, length = length, lw = 0, hp = 1)
colnames(df) <- c("x","y")
df = df[complete.cases(df),]
df = df[df$x>=0.02 & df$x<=0.98, ]

ggplot(data = df, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(formula = y ~ s(x), method = "gam", alpha = 0.2) +
  xlab("methylation proportion")+ 
  ylab("MSR")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

##########################################





# 5: experiment p vs MSR
msr_ecdf_ref <- readRDS("Scrivania/Tesi/MethylationCode/MethylationData/msr_ecdf_1e3.Rda")
wgbs <- "Scrivania/Tesi/MethylationCode/MethylationData/wgbs/stomach.Rda"
wgbs <- sum_strands(readRDS_if_needed(wgbs))
size <- 1e3
df <- readRDS("Scrivania/Tesi/Rexperiments/stomach_msr_complete_experimental_table_1000.Rda")
df <- readRDS("Scrivania/Tesi/Rexperiments/GM23248_msr_complete_experimental_table_1000.Rda")
df <- df[df$missing_prop<0.1 & df$msr_density>0.02 & df$msr_density<0.98 & df$nucleotides<3e5,]
# colnames(df) <- c("x","y")
# df = df[complete.cases(df),]
# df = df[df$x>=0.02 & df$x<=0.98, ]

confidence <- 0.9999
median_function <- extract_ecdf_function(msr_ecdf_ref, 0.5)
up_conf_function <- extract_ecdf_function(msr_ecdf_ref, 1-(1-confidence)/2)
lw_conf_function <- extract_ecdf_function(msr_ecdf_ref, (1-confidence)/2)


px <- seq(from = 0.02, to = 0.98, by = 0.01)
fit = median_function(px)
pred.int <- data.frame(x = px, upr = up_conf_function(px), lwr = lw_conf_function(px), x_rev = 1-px)


color1 = "dark blue"
color2 = "dark red"

p1<-ggplot(data = pred.int, aes(x = x, y = fit)) +
  #ggtitle("Prediction interval for future observations from predict()") +
  geom_smooth(data = pred.int, aes(ymin = lwr, ymax = upr), stat = "identity", linetype=0) +
  geom_point(data = undersample_data_frame(df, n = 1e3), aes(x = msr_density, y = msr), alpha = 0.333, col=color1) +
  #geom_density2d(data = df, aes(x = msr_density, y = msr))+
  xlab("methylation rate")+ 
  ylab(TeX("$MSR_1$"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2<-ggplot(data = pred.int, aes(x = x_rev, y = fit)) +
  #ggtitle("Prediction interval for future observations from predict()") +
  geom_smooth(data = pred.int, aes(ymin = lwr, ymax = upr), stat = "identity", linetype=0) +
  geom_point(data = undersample_data_frame(df, n = 1e3), aes(x = msr_density, y = inverted_msr), alpha = 0.333, col=color2) +
  xlab("methylation rate")+
  ylab(TeX("$MSR_0$"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))


grid.arrange(p1, p2, nrow=2)

##############################################


########################## ecdf0 and ecdf1 HIST 

# basic histogram
h1 <- ggplot(df, aes(x=ecdf)) + 
  geom_histogram(color = "white", fill = color1, alpha = 0.6) +
  xlab(TeX("ecdf_1"))


# basic histogram
h2 <- ggplot(df, aes(x=inverted_ecdf)) + 
  geom_histogram(color = "white", fill = color2, alpha = 0.6) +
  xlab(TeX("ecdf_0")) +
  ylab("") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


grid.arrange(h1, h2, ncol=2)

############################################


####################### MSR_0 and MSR_1 HIST

yl <- 6000
# basic histogram
h1 <- ggplot(df, aes(x=msr)) + 
  geom_histogram(color = "white", fill = color1, alpha = 0.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  xlab(TeX("$MSR_1$")) +
  ylim(NA,yl)


# basic histogram
h2 <- ggplot(df, aes(x=inverted_msr)) + 
  geom_histogram(color = "white", fill = color2, alpha = 0.6) +
  #xlim(0.06,NA)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  ylim(NA,yl)+
  xlab(TeX("$MSR_0$")) +
  ylab("") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


grid.arrange(h1, h2, ncol=2)

############################################

####################### MSR_0 and MSR_1 HIST

# basic histogram
h1 <- ggplot(df, aes(x=msr)) + 
  geom_density(fill = color1, alpha=0.6, adjust = 0.5)+
  ylab("") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(0,0.31))+
  xlab(TeX("$MSR_1$"))+
  theme(axis.title.y=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank())


# basic histogram
h2 <- ggplot(df, aes(x=inverted_msr)) + 
  geom_density(fill = color2, alpha=0.6, adjust = 0.5)+
  ylab("") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(0,0.31))+
  xlab(TeX("$MSR_0$"))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


grid.arrange(h1, h2, nrow=2)

############################################




####################### METH RATE AND MSR DENSITY HIST

l<-length(df$meth_rate)
data <- data.frame(
  type = c(rep("methylation rate", l), rep("proportion of ones", l) ),
  value = c( df$meth_rate, df$msr_density)
)

ggplot(data, aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme(legend.position = c(0.3, 0.5))+
  #theme_ipsum() +
  labs(fill="")

##################################################



####################### BASIC FEATURES CORRELATIONS
response_variable <- fragments_response_variable
basic_predictors <- c("nucleotides", "meth_rate")
advanced_predictors <- c("meth_autocorrelation", "drift", "meth_sd")
msr_related_predictors <- c("CGsites_msr","msr", "inverted_msr", "ecdf", "inverted_ecdf", "residual", "inverted_residual")
features <- c(basic_predictors, msr_related_predictors)
M <- cor((df[,features]), method = "pearson")
fnames<-c("nucleotides","meth rate", "CpG sites MSR","MSR(1)", "MSR(0)", "ecdf(1)", "ecdf(0)", "residual(1)", "residual(0)")
rownames(M)<-fnames
colnames(M)<-fnames
corrplot.mixed(M, lower = "pie", upper = "number")
corrplot(M,type = "upper", method = "number")
##########################################

####################### FEATURES CORRELATIONS WITH AUTOC
response_variable <- fragments_response_variable
basic_predictors <- c("nucleotides", "meth_rate")
msr_related_predictors <- c("msr", "inverted_msr", "ecdf", "inverted_ecdf", "residual", "inverted_residual")
features <- c(basic_predictors, msr_related_predictors, "meth_autocorrelation")
M <- cor(df[,features])["meth_autocorrelation", c(basic_predictors, msr_related_predictors), drop=FALSE]
fnames<-c("nucleotides","meth rate","MSR(1)", "MSR(0)", "ecdf(1)", "ecdf(0)", "residual(1)", "residual(0)")
rownames(M)<-"meth autocorrelation"
colnames(M)<-fnames
corrplot(M, cl.pos='n',  "number")
##########################################


####################### 1000 unif msr vs CG sites msr hists

random_unif_1000_msr <- readRDS("Scrivania/Tesi/Rexperiments/random_unifrom_1000_msr.rda")

l<-length(df$CGsites_msr)
data <- data.frame(
  type = c(rep("1000 random uniform points MSR", length(random_unif_1000_msr)), rep("1000 CpG sites MSR", length(df$CGsites_msr)) ),
  value = c( random_unif_1000_msr, df$CGsites_msr)
)


ggplot(data, aes(x=value, group=type, fill=type)) +
  geom_density(adjust=1.5, alpha=.4)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c("red", "navy")) +
  theme(legend.position = c(0.3, 0.5))+
  #theme_ipsum() +
  labs(fill="")+
  xlab("MSR")

##################################################

####################### 1000 unif msr vs CG sites msr hists

random_unif_10000_msr <- readRDS("Scrivania/Tesi/Rexperiments/random_unifrom_10000_msr.rda")
df <- readRDS("Scrivania/Tesi/Rexperiments/final/GM12878_msr_complete_experimental_table_10000.Rda")

l<-length(df$CGsites_msr)
data <- data.frame(
  type = c(rep("10,000 random uniform points MSR", length(random_unif_10000_msr)), rep("10,000 CpG sites MSR", length(df$CGsites_msr)) ),
  value = c( random_unif_10000_msr, df$CGsites_msr)
)

ggplot(data, aes(x=value, group=type, fill=type)) +
  geom_density(adjust=1.5, alpha=.4)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(NA,0.3))+
  scale_fill_manual(values=c("red", "navy")) +
  theme(legend.position = c(0.25, 0.5))+
  #xlim(NA,0.3)+
  #theme_ipsum() +
  labs(fill="")+
  xlab("MSR")

##################################################


####################### RANDOM VS 

w <- wgbs[139001:(139001+1000),]
fragment_bs <- standard_binaryzer2(w)
random_bs <- rbinom(length(fragment_bs),1,mean(fragment_bs, na.rm=T))==1
rr_curve_comparison(List(which(fragment_bs),which(random_bs)),c("1000 CpGs fragment","random"), add_msr = T, discrete = F, max_gbins = 100)

# rr_calc <- function(v){calculate_relevance_resolution_vector(v, na_tolerance = 0.2,na_values_handler = replace_nas_with_bin_prop, verbose = F)} 
# rr1 <- rr_calc(fragment_bs)
# rr2 <- rr_calc(random_bs)
# 
# compare_resolution_relevance_plot(List(rr1,rr2),c("1000 CpGs", "random"), title = " ")

##################################################


##################################################



# w <- wgbs[1e6:(1e6+1000),]
# Z <- standard_binaryzer2(w)*1
# #Z <- w$prop
# X <- 1:length(Z)
# Y <- 1
# data <- expand.grid(X, Y)
# data$Z <- Z
# 
# # Give extreme colors:
# ggplot(data, aes(X, Y, fill= Z)) + 
#   geom_tile() +
#   scale_fill_gradient(low="white", high="black") +
#   theme_void() +
#   theme(legend.position="none")
#   #+
#   #theme_ipsum()


#########

par(mfrow = c(3,1), mai = c(0.3, 0.3, 0.1, 0.3))

data <- df[df$missing_prop<0.01, ]

indexes <- c(which(data$ecdf==0 & data$inverted_ecdf==0 & abs(data$msr_density-0.5)<0.1)[1],
             which(abs(data$ecdf-0.5)<0.3 & abs(data$inverted_ecdf-0.5)<0.3 & abs(data$msr_density-0.5)<0.1)[1],
             which(data$inverted_ecdf==1 & abs(data$msr_density-0.5)<0.2)[1])

for(i in indexes)
{
  title <- paste("MSR(1):",round(data$msr[i],2),"  MSR(0):",round(data$inverted_msr[i],2),
                 "\necdf(1):",round(data$ecdf[i],2),"  ecdf(0):",round(data$inverted_ecdf[i],2),
                 "\nmeth rate: ", round(data$msr_density[i],2))
  show_fragment_meth(wgbs, data,i,F,F,title)
}

#########

indexes <- c(which(data$ecdf==0 & data$inverted_ecdf==0 & abs(data$msr_density-0.5)<0.1)[2],
             which(data$ecdf==0 & data$inverted_ecdf==0 & abs(data$msr_density-0.5)>0.2)[1],
             which(abs(data$ecdf-0.5)<0.2 & abs(data$inverted_ecdf-0.5)<0.2 & abs(data$msr_density-0.5)<0.1)[2],
             which(abs(data$ecdf-0.5)<0.2 & abs(data$inverted_ecdf-0.5)<0.2 & abs(data$msr_density-0.5)>0.2)[1],
             which(data$inverted_ecdf==1 & abs(data$msr_density-0.5)>0.2)[1],
             which(data$ecdf==1 & abs(data$msr_density-0.5)<0.2)[1]
             )

any(is.na(indexes))

par(mfrow = c(length(indexes),1), mai = c(0.3, 0.3, 0.1, 0.3))
for(i in indexes)
{
  title <- paste("\necdf(1):",round(data$ecdf[i],2),"  ecdf(0):",round(data$inverted_ecdf[i],2))
  show_fragment_meth(wgbs, data,i,F,F,title,-2)
}

#########

data <- df[df$missing_prop<0.01, ]
size <- 1000

indexes <- c(data$i_start[data$ecdf==0 & data$inverted_ecdf==0 & abs(data$msr_density-0.5)<0.1][1],
             data$i_start[abs(data$ecdf-0.5)<0.2 & abs(data$inverted_ecdf-0.5)<0.2 & abs(data$msr_density-0.5)<0.1][1],
             data$i_start[data$ecdf==1 & abs(data$msr_density-0.5)>0.2][1]
)

indexes <- c(which(data$ecdf==0 & data$inverted_ecdf==0 & abs(data$msr_density-0.5)<0.1)[1],
             which(abs(data$ecdf-0.5)<0.3 & abs(data$inverted_ecdf-0.5)<0.3 & abs(data$msr_density-0.5)<0.1)[1],
             which(data$inverted_ecdf==1 & abs(data$msr_density-0.5)<0.2)[1])

indexes <- data$i_start[indexes]

minimum_reads <- 3
genomewide = F
loess_span <- 0.2

pp <- lapply(indexes, function(i_start)
{
  ggplot_meth_fragment(wgbs, i_start, size, T,genomewide, minimum_reads, loess_span)
})

do.call("grid.arrange", c(pp, nrow=length(pp)))



##################################################









# LAST CHAPTER
all_cell_names <- c("H1", "K562", "GM12878", "GM23248", "Hela", "endodermal", "lung", "stomach")


# 1: overall meth levels

cell_names <- c("H1", "endodermal", "lung", "stomach", "Hela", "GM23248", "GM12878", "K562")
color1 <- "#69b3a2"
pp <- lapply(cell_names, function(cn)
{
  df <- data_tables[[cn]]
  df <- data.frame(x=df$meth_rate[!is.na(df$meth_rate)])
  yl <-  NA
  
  ggplot(df, aes(x=x)) + 
    geom_histogram(color = "white", fill = color1, alpha = 0.8) +
    ggtitle(cn) +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    annotate("text", x = 0.5, y = 2500, label = (paste("mean: ", round(mean(df$x),2)))) #+ 
    #ylim(NA,yl)

})

do.call("grid.arrange", c(pp, ncol=3, nrow = 3))

#####################################



# 2: overall expression levels

cell_names <- c("H1")
color1 <- "green"

pp <- lapply(cell_names, function(cn)
{
  v <- response_variable
  df <- data_tables[[cn]]
  x = df[,v]
  df <- data.frame(x=x[!is.na(x)])

  ggplot(df, aes(x=x)) + 
    geom_histogram(color = "white", fill = color1, alpha = 0.6) +
    ggtitle(cn) +
    xlab(TeX("$\\log_{2}(pmeTPM)$"))
    

})

do.call("grid.arrange", c(pp, ncol=1, nrow = 1))

#####################################


# 3: features vs expr matrix

cell_names <- c("H1", "endodermal", "lung", "stomach", "Hela", "GM23248", "GM12878", "K562")
rnames <- c("log(nucleotides)", "CpG density", "meth rate", "meth autocorrelation", "mean entropy", "meth sd",
            "CpG sites msr", "msr(1)","msr(0)","ecdf(1)", "ecdf(0)", "residual(1)","residual(0)")
M <- features_vs_expr_corr_matrix
rownames(M) <- rnames 
corrplot(M, "number")

#####################################



# 4: methylation rate + inverted_msr + expr

cell_names <- c("H1", "endodermal", "lung", "stomach", "Hela", "GM23248", "GM12878", "K562")
rnames <- c("log(nucleotides)", "CpG density", "meth rate", "meth autocorrelation", "mean entropy", "meth sd",
            "CpG sites msr", "msr(1)","msr(0)","ecdf(1)", "ecdf(0)", "residual(1)","residual(0)")



cell_names <- c("lung", "K562")
cell_names <- c("H1", "GM12878")
undersample <- 5000
c1 <- '#666699'
c2 <- '#00b300'
points_alpha <- 0.8

pp <- lapply(cell_names, function(cn)
{
  df <- data_tables[[cn]]
  xname <- "meth_rate"
  yname <- "inverted_msr"
  maxc <- 8
  
  if(!is.na(undersample))
    df <- undersample_data_frame(df, undersample)
  
  df <- df[df[,response_variable]<maxc,]
  
  df <- data.frame(x=df[,xname], y=df[,yname], ltpm=df[,response_variable])
  df <- df[complete.cases(df),]
  #data_table <- data_table[,c(xname, yname, response_variable)]
  p <- ggplot(df, aes(x, y, color=ltpm)) +
    geom_point(alpha = points_alpha) +
    scale_colour_gradient(low = c1, high = "green", na.value = NA, name = TeX("$\\log_2(pmeTPM)$"), limits=c(-4,maxc)) +
    #theme(legend.position=c(0,1), legend.justification=c(0,1)) +
    labs(x = "meth rate", y=TeX("MSR_0") , title = cn)+
    xlim(0,1)
  
  #if(cn==cell_names[1])
   # return(p + theme(legend.position = "none") + ylab(TeX("MSR_0")))
    
  if(cn!=cell_names[length(cell_names)])
    return(p + theme(legend.position = "none"))

  return(p)
  
})

do.call("grid.arrange", c(pp, ncol=1, nrow = 2))

#####################################



# 5: autoc vs expr



cn <- "H1"

  undersample <- 1000
  points_alpha <- 0.3
  df <- data_tables[[cn]]
  xname <- "meth_autocorrelation"
  yname <- response_variable
  maxc <- 8
  limits <- ylim(-4,15)
  
  if(!is.na(undersample))
    df <- undersample_data_frame(df, undersample)
  
  df <- data.frame(x=df[,xname], y=df[,yname])
  df <- df[complete.cases(df),]

  p1 <- ggplot(df, aes(x, y)) +
    geom_point(alpha = points_alpha) +
    xlab("methylation autocorrelation")+
    ylab(TeX("$\\log_2(pmeTPM)$"))+
    ggtitle(cn)+
    geom_smooth() 
  
  ######

  undersample <- 2000
  points_alpha <- 0.1
  
  df <- data_tables[[cn]]
  xname <- "CGsites_msr"
  yname <- response_variable
  maxc <- 8
  
  if(!is.na(undersample))
    df <- undersample_data_frame(df, undersample)
  
  df <- data.frame(x=df[,xname], y=df[,yname])
  df <- df[complete.cases(df),]
  
  p2<- ggplot(df, aes(x, y)) +
    geom_point(alpha = points_alpha) +
    xlab("CpG sites msr")+
    ylab(TeX("$\\log_2(pmeTPM)$"))+
    #ylab("")+
    #geom_smooth()+
    ggtitle(cn) + limits
  

grid.arrange(p1,p2,ncol=2)

#####################################








####### Models data manipulation

# graphic
plot_table <- grid.table

# data clean
max_miss <- 0.1
train_prop <- 0.75
n_splits <- 100
evaluation_metric <- prediction_correlation_score
id <- identifier

l <- length(merged_data_tables[,1])
ba <- mean(complete.cases(merged_data_tables[,basic_predictors]))
ad <- mean(complete.cases(merged_data_tables[,basic_and_advanced_predictors]))
al <- mean(complete.cases(merged_data_tables[,all_predictors]))

cat("\nbasic predictors data: ", ba)
cat("\nbasic and asv predictors data: ", ad)
cat("\nall predictors data: ", al)
cat("\n", max_miss, "max missing data: ", mean(merged_data_tables$missing_prop<max_miss))

merged_df <- merged_data_tables[complete.cases(merged_data_tables[,all_predictors]),]
merged_df <- merged_df[merged_df$missing_prop<max_miss, ]
cat("\nfinal fraction of data:", length(merged_df[,1])/l)

# model function

# linear model results

model_name <- "linear_with_basic_predictors"
#model_name <- "basic"
predictors <- basic_predictors
model_formula <- model_this(response_variable,predictors)
model_lambda <- function(d){{lm(data=d,formula = model_formula)}}
d1 <- validation_on_subsets_data_frame(n_splits, train_prop, merged_df,model_lambda,
                                      evaluation_metric, id, label = model_name)

model_name <- "linear_with_basic_and_advanced_predictors"
#model_name <- "advanced"
predictors <- basic_and_advanced_predictors
model_formula <- model_this(response_variable,predictors)
model_lambda <- function(d){{lm(data=d,formula = model_formula)}}
d2 <- validation_on_subsets_data_frame(n_splits, train_prop, merged_df,model_lambda,
                                       evaluation_metric, id, label = model_name)

model_name <- "linear_with_all_predictors"
#model_name <- "all"
predictors <- all_predictors
model_formula <- model_this(response_variable,predictors)
model_lambda <- function(d){{lm(data=d,formula = model_formula)}}
d3 <- validation_on_subsets_data_frame(n_splits, train_prop, merged_df,model_lambda,
                                       evaluation_metric, id, label = model_name)

linar_models_results <- rbind(d1,d2,d3)


model_name <- "lasso_with_all_predictors"
predictors <- all_predictors
lasso_lambda <- 0.1
model_formula <- model_this(response_variable,predictors)
eval_f <- function(m,d) {prediction_correlation_score.glmnet(m, d, response_variable, predictors)}
model_lambda <- function(d){lasso(response_variable, d[,c(response_variable, predictors)], lasso_lambda)}
d_lasso <- validation_on_subsets_data_frame(n_splits, train_prop, merged_df,model_lambda,
                                            eval_f, id, label = model_name)



# 6: linear models performances
data <- linar_models_results
data$predictors <- as.character(data$model)
data$predictors[data$model=="linear_with_basic_predictors"]<-("ba")
data$predictors[data$model=="linear_with_basic_and_advanced_predictors"]<-("adv")
data$predictors[data$model=="linear_with_all_predictors"]<-("all")
data$predictors = as.factor(data$predictors)
data$model <- data$type

data$r <- data$result
data$r2 <- data$result^2

data$predictors <- relevel(data$predictors, "ba")
cell_names <- c("H1","K562","GM12878")
cell_names <- c("H1", "endodermal", "lung", "stomach", "Hela", "GM23248", "GM12878", "K562")

ggplot(data[data$name %in% cell_names,], aes(x=predictors, y=r2, fill=model)) + 
  #geom_boxplot() +
  geom_violin() +
  facet_wrap(~name,nrow = 2)+
  ylim(0,0.7) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  ylab(TeX("$r^2$"))


#####################################


# genes lasso coefficients
# load the library
library(forcats)

colnames(coeff_matrix) <- c(cell_names, "Shared")
M <- round(coeff_matrix,2)

coefficients <- (M[,"Shared"])
merged_lasso_normalized_coefficients<-lasso_normalized_coefficients(response_variable, merged_model_data$train[,c(response_variable, basic_and_advanced_predictors)], 0.1,alpha)
coefficients <- as.numeric(merged_lasso_normalized_coefficients)

#rnames <- c("log(nucleotides)", "CpG density", "log(CpG count)","meth rate", "meth autocorrelation", "mean entropy", "meth sd",
            #"CpG sites msr", "msr(1)","msr(0)","ecdf(1)", "ecdf(0)", "residual(1)","residual(0)")

rnames <- c("log(nucleotides)", "CpG density", "log(CpG count)", "meth rate", "meth autocorrelation", "mean entropy", "meth sd")

data <- data.frame(name=as.factor(rnames), coefficient=coefficients)
data <- data[order(abs(data$coefficient)),]
data$name <- factor(data$name, levels = data$name)

ggplot(data, aes(x=name, y=coefficient)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  ggtitle("Lasso coefficients")

#####################################








###### GENEBODIES EXPERIMENT




####### Models data manipulation

# graphic
plot_table <- grid.table

# data clean
max_miss <- 0.1
train_prop <- 0.8
n_splits <- 100
evaluation_metric <- prediction_correlation_score
id <- identifier

l <- length(merged_data_tables[,1])
ba <- mean(complete.cases(merged_data_tables[,basic_predictors]))
ad <- mean(complete.cases(merged_data_tables[,basic_and_advanced_predictors]))
al <- mean(complete.cases(merged_data_tables[,all_predictors]))

cat("\nbasic predictors data: ", ba)
cat("\nbasic and asv predictors data: ", ad)
cat("\nall predictors data: ", al)
cat("\n", max_miss, "max missing data: ", mean(merged_data_tables$missing_prop<max_miss))

merged_df <- merged_data_tables[complete.cases(merged_data_tables[,basic_and_advanced_predictors]),]
merged_df <- merged_df[merged_df$missing_prop<max_miss, ]
#merged_df <- merged_df[merged_df$meth_rate>0.1 & merged_df$meth_rate<0.9, ]
cat("\nfinal fraction of data:", length(merged_df[,1])/l)

# model function

# linear model results

model_name <- "linear_with_basic_predictors"
#model_name <- "basic"
predictors <- basic_predictors
model_formula <- model_this(response_variable,predictors)
model_lambda <- function(d){{lm(data=d,formula = model_formula)}}
d1 <- validation_on_subsets_data_frame(n_splits, train_prop, merged_df,model_lambda,
                                       evaluation_metric, id, label = model_name)

model_name <- "linear_with_basic_and_advanced_predictors"
#model_name <- "advanced"
predictors <- basic_and_advanced_predictors
model_formula <- model_this(response_variable,predictors)
model_lambda <- function(d){{lm(data=d,formula = model_formula)}}
d2 <- validation_on_subsets_data_frame(n_splits, train_prop, merged_df,model_lambda,
                                       evaluation_metric, id, label = model_name)


linar_models_results <- rbind(d1,d2)


model_name <- "Lasso"
predictors <- basic_and_advanced_predictors
lasso_lambda <- 0.1
model_formula <- model_this(response_variable,predictors)
eval_f <- function(m,d) {prediction_correlation_score.glmnet(m, d, response_variable, predictors)}
model_lambda <- function(d){lasso(response_variable, d[,c(response_variable, predictors)], lasso_lambda)}
d_lasso <- validation_on_subsets_data_frame(n_splits, train_prop, merged_df,model_lambda,
                                            eval_f, id, label = model_name)



# linear models performances for genebodies
data <- linar_models_results
data$predictors <- as.character(data$model)
data$predictors[data$model=="linear_with_basic_predictors"]<-("basic")
data$predictors[data$model=="linear_with_basic_and_advanced_predictors"]<-("adv")
data$predictors = as.factor(data$predictors)
data$model <- data$type

data$r <- data$result
data$r2 <- data$result^2

data$predictors <- relevel(data$predictors, "basic")
cell_names <- c("H1","K562","GM12878")
cell_names <- c("H1", "endodermal", "Hela", "GM23248", "GM12878", "K562")

ggplot(data[data$name %in% cell_names,], aes(x=predictors, y=r2, fill=model)) + 
  #geom_boxplot() +
  geom_violin() +
  facet_wrap(~name,nrow = 2)+
  ylim(0,0.7) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  ylab(TeX("$R^2$"))


#####################################


# 7: lasso coefficients
# load the library
library(forcats)

colnames(coeff_matrix) <- c(cell_names, "Shared")
M <- round(coeff_matrix,2)

coefficients <- (M[,"Shared"])
merged_lasso_normalized_coefficients<-lasso_normalized_coefficients(response_variable, merged_model_data$train[,c(response_variable, all_predictors)], 0.1,alpha)
coefficients <- as.numeric(merged_lasso_normalized_coefficients)
#rnames <- c("log(nucleotides)", "CpG density", "log(CpG count)", "meth rate", "meth autocorrelation", "mean entropy", "meth sd")
rnames <- c("log(nucleotides)", "CpG density", "meth rate", "meth autocorrelation", "mean entropy", "meth sd",
            "CpG sites MSR", "MSR(1)", "MSR(0)","ecdf(0)","ecdf(1)","residual(1)","residual(0)")

data <- data.frame(name=as.factor(rnames), coefficient=coefficients)
data <- data[order(abs(data$coefficient)),]
data$name <- factor(data$name, levels = data$name)

ggplot(data, aes(x=name, y=coefficient)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  ggtitle("Lasso coefficients")

#####################################


# genes features vs expr matrix

cell_names <- c("H1", "endodermal", "lung", "stomach", "Hela", "GM23248", "GM12878", "K562")
rnames <- c("log(nucleotides)", "CpG density", "meth rate", "meth autocorrelation", "mean entropy", "meth sd",
            "CpG sites msr", "msr(1)","msr(0)","ecdf(1)", "ecdf(0)", "residual(1)","residual(0)")

rnames <- c("log(nucleotides)", "CpG density", "log(CpG count)", "meth rate", "meth autocorrelation", "mean entropy", "meth sd")

M <- features_vs_expr_corr_matrix[1:7,]
rownames(M) <- rnames 
corrplot(M, "number")

#####################################


# genes methylation rate + something + expr

cell_names <- c("H1", "endodermal", "Hela", "GM23248", "GM12878", "K562")
cell_names <- c("H1", "Hela", "GM12878", "K562")
undersample <- 5000
c1 <- '#666699'
c2 <- '#00b300'
points_alpha <- 0.5

sumbsampled_genes <- sample(data_tables[["H1"]]$gene_id, undersample, replace = F)

pp <- lapply(cell_names, function(cn)
{
  df <- data_tables[[cn]]
  xname <- "meth_rate"
  yname <- "mean_bernoulli_entropy"
  maxc <- 8
  
  if(!is.na(undersample))
    df <- df[df$gene_id %in% sumbsampled_genes, ]
  
  df <- df[df[,response_variable]<maxc,]
  
  df <- data.frame(x=df[,xname], y=df[,yname], ltpm=df[,response_variable])
  df <- df[complete.cases(df),]
  #data_table <- data_table[,c(xname, yname, response_variable)]
  p <- ggplot(df, aes(x, y, color=ltpm)) +
    geom_point(alpha = points_alpha) +
    scale_colour_gradient(low = c1, high = "green", na.value = NA, name = TeX("$\\log_2(pmeTPM)$"), limits=c(-4,maxc)) +
    #theme(legend.position=c(0,1), legend.justification=c(0,1)) +
    labs(x = "meth rate", y="mean entropy" , title = cn)+
    xlim(0,1)+
    ylim(0,0.9)
  
  #if(cn==cell_names[1])
  # return(p + theme(legend.position = "none") + ylab(TeX("MSR_0")))
  
  if(cn!=cell_names[length(cell_names)])
    return(p + theme(legend.position = "none")+
             theme(axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank()) )
  
  return(p + theme(legend.position = "none"))
  
})

do.call("grid.arrange", c(pp, ncol=1, nrow = 4))



# genes overall meth levels

cell_names <- c("H1", "endodermal", "Hela", "GM23248", "GM12878", "K562")
color1 <- "#69b3a2"
pp <- lapply(cell_names, function(cn)
{
  df <- data_tables[[cn]]
  df <- data.frame(x=df$meth_rate[!is.na(df$meth_rate)])
  yl <-  NA
  
  ggplot(df, aes(x=x)) + 
    geom_histogram(color = "white", fill = color1, alpha = 0.8) +
    ggtitle(cn) +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) + 
    annotate("text", x = 0.5, y = 2500, label = (paste("mean: ", round(mean(df$x),2)))) #+ 
  #ylim(NA,yl)
  
})

do.call("grid.arrange", c(pp, ncol=3, nrow = 3))

#####################################


## xgboost

M <- t(matrix(data = c(single_boost_results, shared_boost_result), ncol = 2))
colnames(M) <- cell_names
rownames(M) <- c("single", "shared")
cnames <- c("value","model","cell")
lname = "squared"

Mdf <- table_to_frame(M,cnames)
Mdf[,lname]<-"r"

M2df <- table_to_frame(M^2,cnames)
M2df[,lname]<-"R2"

results_table <- rbind(Mdf,M2df)
#####

Mdf$r <- Mdf$value
Mdf$r2 <- Mdf$r^2

Mdf$model <- factor(Mdf$model, levels = c("single","shared"))

ggplot(Mdf, aes(x=cell, y=r2)) +
  geom_bar(stat="identity", position = "identity", fill = "white", alpha=1, width=.4, aes(x=cell, y=r2, color=model)) +
  coord_flip() +
  geom_point(data = Mdf, aes(x=cell, y=r, color=model)) +
  ylab("")+
  xlab("") +
  ggtitle("Gradient Boosting performances")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits = c(0,1))# +
  #theme_bw()


shared_boost_result
single_boost_results

## importance (5x4)
predictors_names <- c("log(nucleotides)","CpG density","meth rate","meth autocorrelation","mean entropy", "meth sd", "MSR(1)", "MSR(0)", "CpG sites MSR")

shared_boost_importance_matrix <- xgb.importance(predictors_names, model = bstDense)
xgb.plot.importance(shared_boost_importance_matrix, rel_to_first = TRUE, xlab = "Relative importance")

##############################





# Error analysis


######### MODEL ANALYSIS VISUALIZATION
formula <- model_this(response_variable, basic_and_advanced_predictors)
df <- merged_data_tables[complete.cases(merged_data_tables[,c(response_variable, basic_and_advanced_predictors)]),]

merged_model <- lm(data=df, formula = formula)
df$predicted <- merged_model$fitted.values




######### MODEL ANALYSIS VISUALIZATION
H1_wgbs <- sum_strands(readRDS("Scrivania/Tesi/MethylationCode/MethylationData/wgbs/ENCFF601NBW_H1.rda"))

GM12878_wgbs <- sum_strands(readRDS("../../MethylationCode/MethylationData/wgbs/ENCFF279HCL_GM12878.Rda"))



cn <- "H1"
wgbs <- H1_wgbs
t <- df[df$cell_name==cn,]
prediction_error <- t$predicted - t[,response_variable]

plotter <- function(x)
{
  ggplot_meth_fragment2(wgbs, x, prop=T,genomewide=F, minimum_reads=1, loess_span = NA)+
    ggtitle(paste(cn, "\ntrue: ",round(x$log_pme_TPM,2), ", predicted: ", round(x$predicted,2), sep = ""))
}

i <- 1

# High Predicted - High True
x <- t[t$log_pme_TPM > 5,]
x <- (x[order(-x$predicted),])[i,]
x <- x[x$CG_count>400, ]
cat("true:",x$log_pme_TPM, "\npred:", x$predicted)

hh <- plotter(x)

# High Predicted - Low True
x <- t[t$log_pme_TPM < -3,]
x <- (x[order(-x$predicted),])[i,]
x <- x[x$CG_count>400, ]
cat("true:",x$log_pme_TPM, "\npred:", x$predicted)

hl <- plotter(x)

# Low Predicted - Low True
x <- t[t$log_pme_TPM < -3,]
x <- x[x$CG_count>400, ]
x <- (x[order(x$predicted),])[i,]
cat("true:",x$log_pme_TPM, "\npred:", x$predicted)

ll <- plotter(x)

# Low Predicted - High True
x <- t[t$log_pme_TPM > 5,]
x <- x[x$CG_count>400, ]
x <- (x[order(x$predicted),])[i,]
cat("true:",x$log_pme_TPM, "\npred:", x$predicted)

lh <- plotter(x)

grid.arrange(hh,hl,lh,ll)




# predicted - true scatter
cell_names <- c("H1", "endodermal", "Hela", "GM23248", "GM12878", "K562")
cn <- "GM12878"
ggplot(data = df[df$cell_name==cn,], aes(x=predicted,y=log_pme_TPM))+
  geom_point(alpha=0.1)+
  geom_smooth(method = lm)+
  ylim(-4,NA)+
  ggtitle(paste(cn, ",  r=",round(cor(df[df$cell_name==cn,]$predicted, df[df$cell_name==cn,response_variable]),2), sep=""))

library(ggplot2)
library(ggExtra)

# predicted - true scatter
cn <- "GM12878"
p <- ggplot(data = df[df$cell_name==cn,], aes(x=predicted,y=log_pme_TPM))+
  geom_point(alpha=0.03)+
  geom_smooth(method = lm)+
  ylim(-4,NA)+
  ylab("true")+
  ggtitle(paste("\n",cn, ",  r=",round(cor(df[df$cell_name==cn,]$predicted, df[df$cell_name==cn,response_variable]),2), sep=""))

ggMarginal(p, type="histogram", margins = "y")

# predicted - true scatter
cn <- "H1"
p <- ggplot(data = df[df$cell_name==cn,], aes(x=predicted,y=log_pme_TPM))+
  geom_point(alpha=0.04)+
  geom_smooth(method = lm)+
  ylim(-4,NA)+
  ylab("true")+
  ggtitle(paste("\n",cn, ",  r=",round(cor(df[df$cell_name==cn,]$predicted, df[df$cell_name==cn,response_variable]),2), sep=""))

ggMarginal(p, type="histogram", margins = "y")

# predicted - true scatter
cn1 <- "H1"
cn2 <- "GM12878"
ids = intersect(data_tables[[cn1]]$gene_id, data_tables[[cn2]]$gene_id)
data = data.frame(x = data_tables[[cn1]][data_tables[[cn1]]$gene_id %in% ids,response_variable],
                  y = data_tables[[cn2]][data_tables[[cn2]]$gene_id %in% ids,response_variable])
ggplot(data=data, aes(x=x,y=y))+
  geom_point(alpha=0.1)+
  geom_abline(col = "red")+
  ylim(-4,NA)+
  xlab(cn1)+
  ylab(cn2)+
  ggtitle(paste("r=",round(cor(data$x, data$y),2), sep=""))

###############




# bonus scatterplots

# predicted - true scatter
cn <- "H1"
under = 1000
df <- data.frame(x = data_tables[[cn]]$CGsites_msr, y= data_tables[[cn]]$log_total_pme_TPM)
df <- undersample_data_frame(df,under)
p1 <- ggplot(data = df, aes(x=x,y=y))+
  geom_point(alpha=0.2)+
  geom_smooth(method = lm)+
  xlim(0.26,0.3)+
  xlab("CpG sites MSR")+
  ylab(TeX("$\\log_2$(pmeTPM)"))+
  ggtitle(paste(cn,":  \nfragments of 10,000 CpGs", "    \nr = ", round(cor(df$x,df$y),2) ,sep=""))
  ylim(-4,16)
  #ylab("true")+
  #ggtitle(paste("\n",cn, ",  r=",round(cor(df[df$cell_name==cn,]$predicted, df[df$cell_name==cn,response_variable]),2), sep=""))


cn <- "GM12878"
df <- data.frame(x = data_tables[[cn]]$meth_autocorrelation, y= data_tables[[cn]]$log_total_pme_TPM)
df <- undersample_data_frame(df,under)
p2 <- ggplot(data = df, aes(x=x,y=y))+
  geom_point(alpha=0.2)+
  #geom_smooth()+
  xlab("meth autocorrelation")+
  ylab(TeX(""))+
  ylim(-4,16)+
  ggtitle(paste(cn,":  \nfragments of 10,000 CpGs", "    \nr = ", round(cor(df$x,df$y),2) ,sep=""))

grid.arrange(p1,p2, ncol=2)


#######################

# predicted - true scatter
cn <- "H1"
under = 1000
df <- data.frame(x = data_tables[[cn]]$meth_autocorrelation, y= data_tables[[cn]]$log_pme_TPM)
df <- df[complete.cases(df),]
df <- undersample_data_frame(df,under)
p1 <- ggplot(data = df, aes(x=x,y=y))+
  geom_point(alpha=0.4)+
  geom_smooth(method = lm)+
  #xlim(0.26,0.3)+
  #xlab("CpG sites MSR")+
  ylab(TeX("$\\log_2$(pmeTPM)"))+
  ggtitle(paste(cn, "    \nr = ", round(cor(df$x,df$y),2) ,sep="")) +
  ylim(-4,NA)
#ylab("true")+
#ggtitle(paste("\n",cn, ",  r=",round(cor(df[df$cell_name==cn,]$predicted, df[df$cell_name==cn,response_variable]),2), sep=""))

cn <- "H1"
under = 1000
df <- data.frame(x = data_tables[[cn]]$meth_rate, y= data_tables[[cn]]$log_pme_TPM)
df <- df[complete.cases(df),]
df <- undersample_data_frame(df,under)
p2 <- ggplot(data = df, aes(x=x,y=y))+
  geom_point(alpha=0.4)+
  geom_smooth(method = lm)+
  #xlim(0.26,0.3)+
  #xlab("CpG sites MSR")+
  ylab(TeX("$\\log_2$(pmeTPM)"))+
  ggtitle(paste(cn, "    \nr = ", round(cor(df$x,df$y),2) ,sep="")) +
  ylim(-4,NA)

#############

