# Load Libraries ####
library(texmex)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(extrafont)
library(MASS)
library(latex2exp)
loadfonts(device = "win")
palette(c("black","purple","cyan","orange"))

# load own functions ####
gpd_fun <- function(marg,par,par_value,family){
  xi <- marg$models[[par]]$coefficients[2]
  sigma <-  exp(marg$models[[par]]$coefficients[1]) 
  
  pdf <- 1 - (1 + xi * par_value / sigma)^(-1/xi)
  
  u <- marg$models[[par]]$threshold
  
  T_ <- c(1:1000,1:1000*1000,1:1000*1000000)
  zN <- u + sigma/xi*((T_*(1-marg$mqu[par]))^xi - 1)
  
  # return_period <- (xi/sigma * (par_value - u) + 1)^1/xi / (1-marg$mqu[par])
  
  return_period <- as.numeric(-(-(xi* (u - par_value - sigma/xi))/sigma)^(1/xi)/(-1 + marg$mqu[par]))
  
  return(return_period)
}

Scenario <- function(index,var.nr,return_graph, family = FALSE){
  if (family == FALSE){
    graph <- var.nr
    j <- return_graph
    j.list <- list(j1,j2,j3,j4,j5,j6)
    T.j <- c(T.1,T.2,T.3)
    
    cat(
      paste(paste('Extreme var:',list(var.1,var.2)[graph]),
            paste("Return period:",T.j[j]/365,"y"),
            '',
            paste(var.1,'Value:', list(j.list[j][[1]][[1]][[index]],j.list[j+3][[1]][[1]][[index]])[graph]),
            paste('Individual rp:',gpd_fun(marg, as.character(list(var.1,var.2)[graph]), as.numeric(list(j.list[j][[1]][[1]][[index]],j.list[j+3][[1]][[1]][[index]])[graph]))/365,'y'),
            '',
            paste(var.2,'Value:', list(j.list[j][[1]][[2]][[index]],j.list[j+3][[1]][[2]][[index]])[graph]),
            paste('Individual rp:',gpd_fun(marg, as.character(list(var.1,var.2)[3-graph]), as.numeric(list(j.list[j][[1]][[2]][[index]],j.list[j+3][[1]][[2]][[index]])[graph]))/365,'y')
            ,sep="\n"))
  }
  else{
    j <- return_graph
    j.list <- list(j10,j11,j12)
    T.j <- c(T.1,T.2,T.3)
    cat(
      paste(paste('Conditional models:',var.1,'&',var.2),
            paste("Return period:",T.j[j]/365,"y"),
            '',
            paste(var.1,'Value:', j.list[j][[1]][[1]][[index]]),
            paste('Individual rp:',gpd_fun(mAll[[var.1]]$margins, as.character(var.1), as.numeric(j.list[j][[1]][[1]][[index]]))/365,'y'),
            '',
            paste(var.2,'value:', j.list[j][[1]][[2]][[index]]),
            paste('Individual rp:',gpd_fun(mAll[[var.1]]$margins, as.character(var.2), as.numeric(j.list[j][[1]][[2]][[index]]))/365,'y'),
            '',
            paste('rp1/rp2:',(gpd_fun(mAll[[var.1]]$margins, as.character(var.1), as.numeric(j.list[j][[1]][[1]][[index]]))/gpd_fun(mAll[[var.1]]$margins, as.character(var.2), as.numeric(j.list[j][[1]][[2]][[index]])))),
            paste('rp2/rp1:',(gpd_fun(mAll[[var.1]]$margins, as.character(var.2), as.numeric(j.list[j][[1]][[2]][[index]]))/gpd_fun(mAll[[var.1]]$margins, as.character(var.1), as.numeric(j.list[j][[1]][[1]][[index]]))))
            ,sep="\n"))
  }
}

likiest.values <- function(value, marges, var, p.var, plot=FALSE){
  l.values <- matrix(0,length(p.var$data$real),1)
  for (i in 1:(length(p.var$data$real))){
    data.var <- p.var$data$simulated
    spec.data.var <- data.var[data.var[colnames(data.var) == var] > (value + marges[1]) & data.var[,colnames(data.var) == var] < (value + marges[2]),]
    id.max <- as.integer(which.max(kde2d(spec.data.var[,var], spec.data.var[,i], n=n)[3][[1]])/n)
    l.values[i,] <- kde2d(spec.data.var[,var], spec.data.var[,i], n=n)[2][[1]][id.max]
    if (plot == TRUE){
      # min_ <- kde2d(p.var.1$data$simulated[,var], p.var.1$data$simulated[,i], n=n)[[2]][1]
      # max_ <- kde2d(p.var.1$data$simulated[,var], p.var.1$data$simulated[,i], n=n)[[2]][200]
      # seq <- sequence(n,min_*10000,((max_-min_)/n*10000))/10000
      # Hs_lim <- kde2d(spec.data.var[,var], spec.data.var[,i], n=n)[[1]][1]
      # plot(seq,kde2d(spec.data.var[,var], spec.data.var[,i], n=n)[[3]][sequence(200,18,200)]/n,
      #      main = paste('T distribution of wind waves given Hs is in', round(Hs_lim,3),'+/- 0.01m'),
      #      xlab="T wind waves [s]",
      #      ylab="Density",
      #      type='line')
      
      image(kde2d(spec.data.var$Hs_ww, spec.data.var[,i], n=n))
      # image(kde2d(data.var$Hs_ww, data.var[,i], n=n),
      #       main="Density plot: 99 percentile Hs vs T",
      #       xlab="Hs wind waves in [m]",
      #       ylab="T wind waves [s]",
      #       xlim = c(4.8,7.5))
    }
  }
  rownames(l.values) <- colnames(spec.data.var)
  cat(paste('particles used for maximum estimation:',length(spec.data.var[,1]),'+ range:', marges[1],',',marges[2]))
  cat('\n\n')
  cat(paste('Likeliest values given',var,'=',value))
  return(l.values)
}




# load data ####
setwd('C:/Users/ianmu/OneDrive - Imperial College London/Thesis/Code/R/Data interpretation')
path <- "../../../../Thesis Data/ERA5/Daily data/"
File_names <- list.files(path = path, pattern="*.csv")
for (i in 1:length(File_names)){
  assign(File_names[i], read.csv(paste(path,File_names[i],sep = ''), header = FALSE))
}
for (i in 1:length(File_names)){
  assign(File_names[i],get(File_names[i])[get(File_names[i])$V1>=692504,])
  assign(File_names[i],get(File_names[i])[get(File_names[i])$V1<=1043088,])
}

# Section on DPF Threshold selection ####

c_finder <- function(index_var, umax_quant = .995){
  name <- names(index_var)
  index_var <- index_var[[1]]
  
  umax <- quantile(index_var, umax_quant)
  grfVar1 <- gpdRangeFit(index_var,umax=umax)
  mrlVar1 <- mrl(index_var)
  g1 <- ggplot(grfVar1)
  g2 <- ggplot(mrlVar1)
  grid.arrange(g1[[1]] + ggtitle(paste("Stability plot, scale parameter",name)),
               g1[[2]] + ggtitle(paste("Stability plot, shape parameter",name)),
               g2 + ggtitle("MRL plot"),ncol=2)
}

data_tba <- na.omit(data.frame(Hs_daily.csv$V2,Hs_ww_daily.csv$V2, dw_daily.csv$V2, Hs_s_m_daily.csv$V2, mod10_daily.csv$V2, P_daily.csv$V2,T_s_m_daily.csv$V2, T_ww_m_daily.csv$V2, surge_daily.csv$V2,dww_m_daily.csv$V2))
colnames(data_tba) <- c('Hs','Hs_ww','dw', 'Hs_s_m','mod10', 'P','T_s_m', 'T_ww_m', 'surge', 'd_ww')
index_var <- data_tba[1]
# Deleting Westerly data 
data_tba = data_tba[!((data_tba$d_ww<(-1+(pi/2)))*(data_tba$d_ww>(-1-(pi/2)))),]
c_finder(index_var,.99)
ecdf(index_var[[1]])(3)
quantile(index_var[[1]],.9)

# GPD_c <- c(3, .95, 2.25, 14, 3.5, 10, 6, .6)
# GPD_q <- c(0.9161418, 0.9896701, 0.9940134, 0.9067966, 0.9900223, 0.9641977, 0.862073, 0.981688)

GPD_c <- c(2.790131, 0.9542967, 1.385238, 13.83927, 0.5, 8.197759, 5.599755, 0.5932663, 2.910943)
GPD_q <- c(0.90, 0.99, 0.90, 0.90, 0.70, 0.80, 0.80, 0.98, 0.90)

# Section on preliminary dependency analyses ####
index_name <- 'Hs'
data_tba <- na.omit(data.frame(Hs_daily.csv$V2, Hs_ww_daily.csv$V2))
GGally::ggpairs(data_tba)
# chi loop for all vars on index var
for (i in 1:(length(data_tba))){
  chi <- chi(data_tba[, c(index_name, colnames(data_tba[i]))])
  assign(paste('chi.',index_name,'.',colnames(data_tba[i]),sep=''), chi)
  ggplot(chi, main=c("Chi"=paste("Chi: ", index_name,'.', colnames(data_tba[i]), sep=''), "ChiBar"=paste("Chi-bar: ",index_name,'.' ,colnames(data_tba[i]), sep='')))
}

# spearman loop for all vars on index var
for (i in 1:(length(data_tba))){
  spear <- bootMCS(data_tba[, c(index_name, colnames(data_tba[i]))],trace=1000,R=10)
  assign(paste('spear.', index_name,'.' ,colnames(data_tba[i]), sep=''), ggplot(spear, main=paste("MCS: ",index_name,'.', colnames(data_tba[i]), sep='')))
  gridExtra::grid.arrange(ggplot(spear, main=paste("MCS: ",index_name,'.', colnames(data_tba[i]), sep='')))
}


# Start scenario creation: Fitting ####
data_tba <- data.frame(Hs_ww_daily.csv$V2, dw_daily.csv$V2, Hs_s_m_daily.csv$V2, mod10_daily.csv$V2, P_daily.csv$V2,T_s_m_daily.csv$V2, T_ww_m_daily.csv$V2, surge_daily.csv$V2,Hs_daily.csv$V2)
colnames(data_tba) <- c('Hs_ww','dw', 'Hs_s_m','mod10', 'P','T_s_m', 'T_ww_m', 'surge','Hs')

# # Optional deletion of Westerly data 
# data_tba = data_tba[!((dww_m_daily.csv$V2<(-1+(pi/2)))*(dww_m_daily.csv$V2>(-1-(pi/2)))),]
# data_tba = na.omit(data_tba)

# mqu <- 0.95 # Threshold for pareto
# dqu.1 <- 0.8 # Thresholf for dependence var 1
# dqu.2 <- 0.8 # Thresholf for dependence var 2

mqu <- GPD_q # Threshold for pareto
dqu.1 <- .90 # Thresholf for dependence var 1
dqu.2 <- .90 # Thresholf for dependence var 2

marg <- migpd(data_tba, mqu=mqu, penalty = 'gaussian', family = gpd)

## marginals
g <- ggplot(marg)
for (i in 1:length(g)){
  do.call("grid.arrange", c(g[[i]], list(ncol=2, nrow=2)))
  ggplot(gpdRangeFit(data_tba[,i]))
  ggplot(mrl(data_tba[,i]))
}

## Defining dependence
# First choose 2 mean variables for return contour selection
var.1 <- 'Hs_ww'
var.2 <- 'T_ww_m'

start_mat.1 <- matrix(c(0.1, 0.1,# 
                        0.1, 0.1,
                        0.95, 0.0, # 0.5, 0.6
                        0.1, 0.1,
                        0.1, 0.1,
                        0.1, 0.1, # 0.5, 0.4
                        0.4, 0.2,
                        0.0, 0.0 ), # 0.4, 0.5
                      2,(length(data_tba)-1))

start_mat.2 <- matrix(c(0.1, 0.1,# checked for P with mqu = dqu = .95
                        0.1, 0.1,
                        0.1, 0.1,
                        0.1, 0.1,
                        0.1, 0.1,
                        0.1, 0.1,
                        0.1, 0.1,
                        0.1, 0.1),
                      2,(length(data_tba)-1))

# calculate dependence on var.1 & 2
Plot_parameter_contours <- TRUE
mex.var.1 <- mexDependence(marg, dqu=dqu.1, which = var.1, PlotLikDo=Plot_parameter_contours, 
                           start = start_mat.1, PlotLikRange=list(a=c(-1,1),b=c(-0.8,0.8)),v=12)
mex.var.2 <- mexDependence(marg, dqu=dqu.2, which = var.2, PlotLikDo=Plot_parameter_contours, 
                           start = start_mat.2, PlotLikRange=list(a=c(-1,1),b=c(-1,1)))

ggplot(mex.var.1)
ggplot(mex.var.2)

# bootstrap threshhold & start value check
R <- 1
par_to_check <- 1
range <- seq(0.8, 0.99, length = 10)
mrf <- mexRangeFit(marg, c(var.1,var.2)[par_to_check], trace=5, R=R, quantiles= range, start = list(start_mat.1, start_mat.2)[[par_to_check]])
ggplot(mrf)

# dependence interpretation using 2 variables sepparately (using a family later on!) ####
# Create predictions based on distributions
nsim.1 <- 100000
nsim.2 <- 100000

pqu.1 <- .99
pqu.2 <- .99

p.var.1 <- predict(mex.var.1, pqu=pqu.1, nsim=nsim.1)
p.var.2 <- predict(mex.var.2, pqu=pqu.2, nsim=nsim.2)

T.1 <- 10*365
T.2 <- 25*365
T.3 <- 100*365
n <- 1000/2

j1 <- JointExceedanceCurve(p.var.1,1/(T.1),which=c(var.1,var.2),n=n)
j2 <- JointExceedanceCurve(p.var.1,1/(T.2),which=c(var.1,var.2),n=n)
j3 <- JointExceedanceCurve(p.var.1,1/(T.3),which=c(var.1,var.2),n=n)
j4 <- JointExceedanceCurve(p.var.2,1/(T.1),which=c(var.2,var.1),n=n)
j5 <- JointExceedanceCurve(p.var.2,1/(T.2),which=c(var.2,var.1),n=n)
j6 <- JointExceedanceCurve(p.var.2,1/(T.3),which=c(var.2,var.1),n=n)

# adding predicted data from var.2 to g for nicer plotting
p.var.1$data[2]$simulated <- rbind(p.var.1$data[2]$simulated, p.var.2$data[2]$simulated)
p.var.1$data[2]$simulated[colnames(p.var.1$data[2]$simulated) == var.2] <- c(p.var.1$data[2]$simulated[,colnames(p.var.1$data[2]$simulated)==var.2][1:nsim.1],p.var.2$data[2]$simulated[,colnames(p.var.2$data[2]$simulated)==var.2])
p.var.1$data[2]$simulated[colnames(p.var.1$data[2]$simulated) == var.1] <- c(p.var.1$data[2]$simulated[,colnames(p.var.1$data[2]$simulated)==var.1][1:nsim.1],p.var.2$data[2]$simulated[,colnames(p.var.2$data[2]$simulated)==var.1])
p.var.1$data$CondLargest <- 0

# p.var.1$data[2]$simulated <- matrix(0,nsim.1,7) # delete all simulations
g <- ggplot(p.var.1,plot.=FALSE)

pl <- g[[(1:(length(names(p.var.1$mth))-1))[names(p.var.1$mth)[names(p.var.1$mth)!=var.1] == var.2]]] + # change this number to the right y variable
  geom_jointExcCurve(j1,aes(!!as.name(var.1),!!as.name(var.2)),col="red") +
  geom_jointExcCurve(j2,aes(!!as.name(var.1),!!as.name(var.2)),col="red") +
  geom_jointExcCurve(j3,aes(!!as.name(var.1),!!as.name(var.2)),col="red") + 
  geom_jointExcCurve(j4,aes(!!as.name(var.1),!!as.name(var.2)),col="red") +
  geom_jointExcCurve(j5,aes(!!as.name(var.1),!!as.name(var.2)),col="red") +
  geom_jointExcCurve(j6,aes(!!as.name(var.1),!!as.name(var.2)),col="red")
pl$layers[[5]]$aes_params$colour <- 'black'
pl

# Joint exceedance curve from bi family of conditional models ####
data_All <- data.frame(data_tba[var.1], data_tba[var.2])
mAll <- mexAll(data_All,mqu=0.9,dqu=rep(0.9,2))

p_all <- mexMonteCarlo(nSample=(10000*365),mexList=mAll)
j10 <- JointExceedanceCurve(p_all,1/(T.1),which=c(var.1,var.2),n = n)
j11 <- JointExceedanceCurve(p_all,1/(T.2),which=c(var.1,var.2),n = n)
j12 <- JointExceedanceCurve(p_all,1/(T.3),which=c(var.1,var.2),n = n)

Scenario(index = 437,var.nr = 1,return_graph = 3 ,family = TRUE) 

# dots <- data.frame(c(7.98559858540739,7.75409195018658,7.11857486512342,6.56106886100701,5.80894633303791,4.44383398590048,-4.44089209850063e-16),
#                    c(0.105061787402089,1.84222157076779,2.61100717702325,2.96780046747323,3.25412590533311,3.76905128974528,3.79431257962869))



T_ratio <- c(100,50,10,1,1/10,1/50,1/100)
Curve <- j12

dots <- data.frame(t(data.frame(c(Curve[[1]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[1]))], Curve[[2]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[1]))]),
           c(Curve[[1]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[2]))], Curve[[2]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[2]))]),
           c(Curve[[1]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[3]))], Curve[[2]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[3]))]),
           c(Curve[[1]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[4]))], Curve[[2]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[4]))]),
           c(Curve[[1]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[5]))], Curve[[2]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[5]))]),
           c(Curve[[1]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[6]))], Curve[[2]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[6]))]),
           c(Curve[[1]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[7]))], Curve[[2]][which.min(abs(gpd_fun(marg,var.1,Curve[[1]][c(1:(2*n))],TRUE)/gpd_fun(marg,var.2,Curve[[2]][c(1:(2*n))],TRUE) - T_ratio[7]))]))))
colnames(dots) <- c(var.1,var.2)
rownames(dots) <- c(1,2,3,4,5,6,7)

# dots <- data.frame(c(1.78076171702584e-05,8.01888854631365,5.98789256130688,7.8139006654471,3.30278275549826,6.77320304849972,5.11528548717413),
                   # c(8.77606727334638,0,5.15635080854053,2.92985154538487,7.4334858455457,4.40696759693255,6.02877173079173))

# dots_names <- c('T1/T2 = 0','T1/T2 = inf','T1/T2 = 1','T1/T2 = 500','T1/T2 = 1/500','T1/T2 = 10','T1/T2 = 1/10')
dots_names <- c('T1/T2 = 100','T1/T2 = 50','T1/T2 = 10','T1/T2 = 1','T1/T2 = 1/10','T1/T2 = 1/50','T1/T2 = 1/100')

pl <- ggplot(as.data.frame(p_all$MCsample[,c(var.1,var.2)]),aes(!!as.name(var.1),!!as.name(var.2))) +
  geom_point(data = as.data.frame(p_all$MCsample[,c(var.1,var.2)]), col="orange",alpha=0.3)  +
  geom_point(data = p.var.1$data$real[colnames(p.var.1$data$real) == var.1 | colnames(p.var.1$data$real) == var.2], col = 'red',alpha=0.5) + 
  geom_jointExcCurve(j12) + 
  geom_point(data = dots, aes(col = dots_names), size=3.5) + 
  ggtitle(expression('Joint exceedance scenarios: T'['c']*'= 100y: H'['s,wind']*' & T'['mean,wind'])) +
  labs(x = expression('H'['s,wind']*"[m]"),
       y = expression('T'['mean,wind']*"[s]"),
       color = "Legend") +
  theme_minimal() +
  theme(text         = element_text(size=11, family="LM Roman 10"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(face="bold"))
pl

# print 2 var scenario ####
Scenario(index = 1000,var.nr = 1,return_graph = 3) #189
Scenario(index = 400,var.nr = 1,return_graph = 3 ,family = TRUE) #189

# Find likeliest var values given scenario occurs
chosen.var.val <- dots[1,]
marges <- c(-0.1, 0.1)
n <- 200

likiest.values(value = as.numeric(chosen.var.val[1]), marges = marges, var = var.1, p.var = p.var.1, plot = FALSE)
likiest.values(value = as.numeric(chosen.var.val[2]), marges = marges, var = var.2, p.var = p.var.2, plot = FALSE)







# Joint exceedance curve 3 model adaptation (unused by the thesis) ####
var.1 <- 'Hs'
var.2 <- 'Hs_s_m'
var.3 <- 'Hs_ww'

svar.1 <- 'Hs_ww'
svar.2 <- 'Hs_s_m'
getvar <- 'Hs'

val.1 <- c(5.80894633303791*1/5,5.80894633303791*2/5,5.80894633303791*3/5,5.80894633303791*4/5,5.80894633303791*5/5)#dots[2:6,1]
val.2 <- c(3.25412590533311*1/5,3.25412590533311*2/5,3.25412590533311*3/5,3.25412590533311*4/5,3.25412590533311*5/5)#dots[1:6,1]

marges <- c(-0.25, 0.25)

data_All <- data.frame(data_tba[var.1], data_tba[var.2], data_tba[var.3])
mAll <- mexAll(data_All,mqu=0.9,dqu=rep(0.9,3))

dis3 <- 0
dis3_list <- list(dis3,dis3,dis3,dis3,dis3)
for (j in 1:50){
  p_all3 <- mexMonteCarlo(nSample=(10*365),mexList=mAll)
  
  for (i in 1:5){
    c1 <- (p_all3$MCsample[names(p_all3$MCsample)==svar.1] > val.1[i] + marges[1])*(p_all3$MCsample[names(p_all3$MCsample)==svar.1] < val.1[i] + marges[2])
    c2 <- (p_all3$MCsample[names(p_all3$MCsample)==svar.2] > val.2[i] + marges[1])*(p_all3$MCsample[names(p_all3$MCsample)==svar.2] < val.2[i] + marges[2])
    # dis3 <- c(dis3, p_all3$MCsample[(p_all3$MCsample*(c1*c2))[1]!=0,][getvar][,1])
    if (i == 1){
      dis3_list <- list(c(dis3_list[[1]],p_all3$MCsample[(p_all3$MCsample*(c1*c2))[1]!=0,][getvar][,1]),dis3_list[[2]],dis3_list[[3]],dis3_list[[4]],dis3_list[[5]])
    }
    if (i == 2){
      dis3_list <- list(dis3_list[[1]],c(dis3_list[[2]],p_all3$MCsample[(p_all3$MCsample*(c1*c2))[1]!=0,][getvar][,1]),dis3_list[[3]],dis3_list[[4]],dis3_list[[5]])
    }
    if (i == 3){
      dis3_list <- list(dis3_list[[1]],dis3_list[[2]],c(dis3_list[[3]],p_all3$MCsample[(p_all3$MCsample*(c1*c2))[1]!=0,][getvar][,1]),dis3_list[[4]],dis3_list[[5]])
    }
    if (i == 4){
      dis3_list <- list(dis3_list[[1]],dis3_list[[2]],dis3_list[[3]],c(dis3_list[[4]],p_all3$MCsample[(p_all3$MCsample*(c1*c2))[1]!=0,][getvar][,1]),dis3_list[[5]])
    }
    if (i == 5){
      dis3_list <- list(dis3_list[[1]],dis3_list[[2]],dis3_list[[3]],dis3_list[[4]],c(dis3_list[[5]],p_all3$MCsample[(p_all3$MCsample*(c1*c2))[1]!=0,][getvar][,1]))
    }
    }
  # print(paste('datapoints obtained:',length(dis3)-1))
  # if (length(dis3)-1 >= 100){
  #   break
  # }
  print(j)
  print(dis3_list)
}
# dis3 <- dis3[2:length(dis3)]
dis3

A <- hist(dis3,n = 10,)
B <- (A$breaks[which(A$counts==max(A$counts))] + A$breaks[which(A$counts==max(A$counts))+1])/2
# C <- c(val.1, val.2, B)
C <- data.frame(C, c(val.1, val.2, B))
