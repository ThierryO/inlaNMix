# code for Meehan, Michel, and Rue #############################################
# demonstration of R-INLA for analyzing N-mixture models
# 04/19/2017
# ##############################################################################





# directories, libraries, and helpers ##########################################
setwd("~/GitHub/Quantitative_Metrics/INLA_NMix")
library(runjags)
library(INLA)
library(unmarked)
library(ggplot2)
# multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# plotting theme
theme_acbs <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey85", colour = "grey20"),
          legend.key = element_rect(fill = "white", colour = NA),
          complete = TRUE)
}
# data generating function
data4sim <- function(n.sites = 72,   # number study sites
                     n.surveys = 3,  # short term replicates
                     n.years = 9,    # number years, MAKE ODD NUMBER
                     lam.b0 = 2.0,   # intercept for log lambda
                     lam.b1 = 2.0,   # slope for log lambda, covariate 1
                     lam.b2 = -3.0,  # slope for log lambda, covariate 2
                     lam.b4 = 1.0,   # slope for log lambda, year
                     p.b0 = 1.0,     # intercept for logit p
                     p.b1 = -2.0,    # slope for logit p, covariate 1
                     p.b3 = 1.0,     # slope for logit p covariate 3
                     disp.size = 3.0 # size of the overdisperison
                     ){
  # setup
  if(n.years %% 2 == 0) {n.years <- n.years + 1}     # make years odd
  N.tr <- array(dim = c(n.sites, n.years))           # array for true abund
  y1 <- array(dim = c(n.sites, n.surveys, n.years))  # array for eg1 counts
  y2 <- array(dim = c(n.sites, n.surveys, n.years))  # array for eg2 counts
  # abundance covariate values
  x.lam.1 <- array(as.numeric(scale(runif(n=n.sites, -0.5, 0.5), scale=F)),
                   dim=c(n.sites, n.years)) # site-level covariate 1
  x.lam.2 <- array(as.numeric(scale(runif(n=n.sites, -0.5, 0.5), scale=F)),
                   dim=c(n.sites, n.years)) # site-level covariate 1
  yrs <- 1:n.years; yrs <- (yrs - mean(yrs)) / (max(yrs - mean(yrs))) / 2
  yr <- array(rep(yrs, each=n.sites), dim=c(n.sites, n.years)) # std years
  # fill abundance array
  lam.tr <- exp(lam.b0 + lam.b1*x.lam.1 + lam.b2*x.lam.2 + lam.b4*yr)	# true lam
  for(i in 1:n.sites){
    for(k in 1:n.years){
    N.tr[i, k] <- rnbinom(n = 1, mu = lam.tr[i, k], size = disp.size) # true N
  }}
  # detection covariate values
  x.p.1 <- array(x.lam.1[,1], dim=c(n.sites,n.surveys,n.years))
  x.p.3 <- array(as.numeric(scale(runif(n=n.sites*n.surveys*n.years, -0.5, 0.5),
                                  scale=F)), dim=c(n.sites,n.surveys,n.years))
  # average x.p.3 per site-year
  x.p.3.mean <- apply(x.p.3, c(1,3), mean, na.rm=F)
  out1 <- c()
  for(k in 1:n.years){
   chunk1 <- x.p.3.mean[,k]
   chunk2 <- rep(chunk1, n.surveys)
   out1 <- c(out1, chunk2)
  }
  x.p.3.arr <- array(out1, dim=c(n.sites, n.surveys, n.years))
  # fill count array with site-yr x.p.3
  p.tr1 <- plogis(p.b0 + p.b1*x.p.1 + p.b3*x.p.3.arr) # true p
  for (i in 1:n.sites){
    for (k in 1:n.years){
      for (j in 1:n.surveys){
        y1[i,j,k] <- rbinom(1, size=N.tr[i,k], prob=p.tr1[i,j,k])
  }}}
  # fill count array with site-surv-yr x.p.3
  p.tr2 <- plogis(p.b0 + p.b1*x.p.1 + p.b3*x.p.3)     # true p
  for (i in 1:n.sites){
    for (k in 1:n.years){
      for (j in 1:n.surveys){
        y2[i,j,k] <- rbinom(1, size=N.tr[i,k], prob=p.tr2[i,j,k])
  }}}
  # return data
  return(list(n.sites=n.sites, n.surveys=n.surveys, n.years=n.years,
              x.p.1=x.p.1[,1,1], x.p.3=x.p.3, x.p.3.mean=x.p.3.mean,
              x.p.3.arr=x.p.3.arr, x.lam.1=x.lam.1[,1],
              x.lam.2=x.lam.2[,1], yr=yr[1,], y1=y1, y2=y2,
              lam.tr=lam.tr, N.tr=N.tr))
}
# ##############################################################################





# create dataset for examples 1 and 2 ##########################################
set.seed(12345)
sim.data <- data4sim()
# ##############################################################################





# jags analysis, example 1 #####################################################
jags.model.string <- "
model {
  # priors
  intP ~ dnorm(0, 0.01)       # detection intercept
  bCov1P ~ dnorm(0, 0.01)     # detection cov 1 effect
  bCov3P ~ dnorm(0, 0.01)     # detection cov 3 effect
  intLam ~ dnorm(0, 0.01)     # lambda intercept
  bCov1Lam ~ dnorm(0, 0.01)   # lambda cov 1 effect
  bCov2Lam ~ dnorm(0, 0.01)   # lambda cov 2 effect
  bYr ~ dnorm(0, 0.01)        # year effect
  overDisEst ~ dunif(0, 5)    # overdispersion size
  # abundance component
  for (k in 1:nYears){
    for (i in 1:nSites){
      N[i, k] ~ dnegbin(prob[i, k], overDisEst) # negative binomial
      prob[i, k] <- overDisEst / (overDisEst + lambda[i, k]) # overdispersion
      log(lambda[i, k]) <- intLam + (bCov1Lam * x.lam.1[i]) +
                            (bCov2Lam * x.lam.2[i]) + (bYr * yr[k])
  # detection component
      for (j in 1:nSurveys){
        y[i, j, k] ~ dbin(p[i,j,k], N[i,k])
        p[i, j, k] <- exp(lp[i,j,k]) / (1 + exp(lp[i,j,k]))
        lp[i, j, k] <- intP + (bCov1P * x.p.1[i]) + (bCov3P * x.p.3[i, k])
      } # close j loop
    } # close i loop
  } # close k loop
} # close model loop
"
# parameters to monitor
params <- c("intP", "bCov1P", "bCov3P", "intLam", "bCov1Lam","bCov2Lam",
            "bYr", "overDisEst")
# jags data
jags.data <- list(y = sim.data$y1, x.lam.1 = sim.data$x.lam.1,
             x.lam.2 = sim.data$x.lam.2, yr = sim.data$yr,
             x.p.1 = sim.data$x.p.1, x.p.3 = sim.data$x.p.3.mean,
             nSites = sim.data$n.sites, nSurveys = sim.data$n.surveys,
             nYears = sim.data$n.years)
# initial values
N.init <- sim.data$y1 # initial count values
N.init[is.na(N.init)] <- 1 # clean up NA's
N.init <- apply(N.init, c(1, 3), max) + 1 # zero values cause trouble
inits <- function() list(N = N.init, intLam = rnorm(1, 0, 0.01),
                     intP = rnorm(1, 0, 0.01), bCov1P = rnorm(1, 0, 0.01),
                     bCov2Lam = rnorm(1,0,0.01), bCov1Lam = rnorm(1, 0, 0.01),
                     bCov3P = rnorm(1, 0, 0.01), bYr = rnorm(1, 0, 0.01),
                     overDisEst = runif(1, 0.5, 2.5))
# set run parameters
nc <- 3; na <- 2500; nb <- 2500; ni <- 5000; nt <- 10
# run jags
ptm <- proc.time()
out.jags <- run.jags(model = jags.model.string, data = jags.data,
                     monitor = params, n.chains = nc, inits = inits,
                     burnin = nb, adapt = na, sample = ni, thin = nt,
                     modules = "glm on", method = "parallel")
# view summary
summary(out.jags)
# computing time
round(jags.time.1 <- proc.time() - ptm, 2)[3]
# ##############################################################################





# unmarked analysis, example 1##################################################
# format count data
y.unmk <- sim.data$y1[ , , 1]
for(i in 2:sim.data$n.years){
  y.chunk <- sim.data$y1[ , , i]
  y.unmk <- rbind(y.unmk, y.chunk)
}
# format covariates
x.lam.1.unmk <- rep(sim.data$x.lam.1, sim.data$n.years)
x.lam.2.unmk <- rep(sim.data$x.lam.2, sim.data$n.years)
yr.unmk <- rep(sim.data$yr, each = sim.data$n.sites)
x.p.1.unmk <- rep(sim.data$x.p.1, sim.data$n.years)
x.p.3.unmk <- c(sim.data$x.p.3.mean)
site.covs.unmk <- data.frame(x.lam.1.unmk, x.lam.2.unmk, yr.unmk,
                             x.p.1.unmk, x.p.3.unmk)
# make unmarked pcount frame
unmk.data <- unmarkedFramePCount(y = y.unmk, siteCovs = site.covs.unmk)
# run unmk model
ptm <- proc.time()
out.unmk <- pcount(~ x.p.1.unmk + x.p.3.unmk
                   ~ x.lam.1.unmk + x.lam.2.unmk + yr.unmk,
                    data = unmk.data, mixture = "NB")
# view summary
summary(out.unmk)
# computing time
round(unmk.time.1 <- proc.time() - ptm, 2)[3]
# ##############################################################################





# inla analysis, example 1 #####################################################
# format count data
y.inla <- sim.data$y1[ , , 1]
for(i in 2:sim.data$n.years){
  y.chunk <- sim.data$y1[ , , i]
  y.inla <- rbind(y.inla, y.chunk)
}
# format covariates
x.lam.1.inla <- rep(sim.data$x.lam.1, sim.data$n.years)
x.lam.2.inla <- rep(sim.data$x.lam.2, sim.data$n.years)
yr.inla <- rep(sim.data$yr, each = sim.data$n.sites)
x.p.1.inla <- rep(sim.data$x.p.1, sim.data$n.years)
x.p.3.inla <- c(sim.data$x.p.3.mean)
# make inla.mdata object
counts.and.count.covs <- inla.mdata(y.inla, 1, x.lam.1.inla, x.lam.2.inla,
                                    yr.inla)
# run inla model
ptm <- proc.time()
out.inla <- inla(counts.and.count.covs ~ 1 + x.p.1.inla + x.p.3.inla,
         data = list(counts.and.count.covs = counts.and.count.covs,
                     x.p.1.inla = x.p.1.inla, x.p.3.inla = x.p.3.inla),
         family = "nmixnb",
         control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.01,
                              prec.intercept = 0.01),
         control.family = list(hyper = list(theta1 = list(param = c(0, 0.01)),
                                          theta2 = list(param = c(0, 0.01)),
                                          theta3 = list(param = c(0, 0.01)),
                                          theta4 = list(param = c(0, 0.01)),
                                          theta5 = list(prior = "flat",
                                                        param = numeric()))))
# view summary
summary(out.inla)
# computing time
round(inla.time.1 <- proc.time() - ptm, 2)[3]
# ##############################################################################





# compile results, example 1 ###################################################
# get jags out
jags.out.tab <- round(as.data.frame(summary(out.jags)[,c(2, 1, 3)]), 3)
# get unmk out
e1 <- as.matrix(out.unmk@estimates@estimates$det@estimates)
ci1 <- as.matrix(confint(out.unmk, type="det"))
e2 <- as.matrix(out.unmk@estimates@estimates$state@estimates)
ci2 <- as.matrix(confint(out.unmk, type="state"))
e3 <- as.matrix(out.unmk@estimates@estimates$alpha@estimates)
ci3 <- as.matrix(confint(out.unmk, type="alpha"))
e4 <- rbind(e1,e2,e3)
ci4 <- rbind(ci1,ci2,ci3)
unmk.out.tab <- cbind(e4,ci4)
# get inla out
in1 <- out.inla$summary.fixed[,c(4, 3, 5)]
in2 <- out.inla$summary.hyperpar[,c(4, 3, 5)][-5,]
in3 <- 1/out.inla$summary.hyperpar[,c(4, 3, 5)][5,]
in3 <- in3[,c(1,3,2)]
inla.out.tab <- rbind(in1,in2,in3)
# combine out
all.out.tab <- round(cbind(jags.out.tab, unmk.out.tab, inla.out.tab),2)
colnames(all.out.tab) <- c("JAGS Median","JAGS Lower CrL", "JAGS Upper CrL",
                           "UNMK Estimate","UNMK Lower CL", "UNMK Upper CL",
                           "INLA Median","INLA Lower CrL", "INLA Upper CrL")
rownames(all.out.tab) <- c("P Int", "P Cov 1", "P Cov 3", "Lambda Int",
                           "Lambda Cov 1", "Lambda Cov 2", "Lambda Year",
                           "Overdispersion")
# get output from jags
jags.mcmc <- combine.mcmc(out.jags, thin=3, return.samples=5000)
jags.df <- as.data.frame(jags.mcmc)
# make density for intLam
d1 <- as.data.frame(out.inla$marginals.hyperpar$`beta[1] for NMix observations`)
d1$x <- (d1$x)
d2 <- as.data.frame(density(jags.df$intLam)[c(1,2)])
d2$x <- (d2$x)
col1 <- 4
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- (c(d3,d4,d5))
d7 <- data.frame(y=c(0,0,0), x=d6)
p1 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab("log(lambda) intercept") +
  ylab(" ") + scale_y_continuous(breaks=c(3,6,9)) + theme_acbs() +
  geom_vline(xintercept=(2), col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) + theme(axis.title.x = element_text(size=9)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))
# make density for x.1.lam
d1 <- as.data.frame(out.inla$marginals.hyperpar$`beta[2] for NMix observations`)
d2 <- as.data.frame(density(jags.df$bCov1Lam)[c(1,2)])
col1 <- 5
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p2 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab("log(lambda) covariate 1")  +
  ylab("Posterior density") + theme_acbs() +
  geom_vline(xintercept=2, col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) + theme(axis.title.x = element_text(size=9)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))
# make density for x.2.lam
d1 <- as.data.frame(out.inla$marginals.hyperpar$`beta[3] for NMix observations`)
d2 <- as.data.frame(density(jags.df$bCov2Lam)[c(1,2)])
col1 <- 6
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p3 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab("log(lambda) covariate 2")  +
  ylab(" ") + theme_acbs() + theme(axis.title.x = element_text(size=9)) +
  geom_vline(xintercept=-3, col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))
# make density for yr.lam
d1 <- as.data.frame(out.inla$marginals.hyperpar$`beta[4] for NMix observations`)
d2 <- as.data.frame(density(jags.df$bYr)[c(1,2)])
col1 <- 7
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p4 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab("log(lambda) year covariate")  +
  ylab(NULL) + theme_acbs() + theme(axis.title.x = element_text(size=9)) +
  geom_vline(xintercept=1, col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))
# make density for intP
d1 <- as.data.frame(out.inla$marginals.fixed$`(Intercept)`)
d1$x <- (d1$x)
d2 <- as.data.frame(density(jags.df$intP)[c(1,2)])
d2$x <- (d2$x)
col1 <- 1
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- (c(d3,d4,d5))
d7 <- data.frame(y=c(0,0,0), x=d6)
p5 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab("logit(p) intercept")  +
  ylab(NULL) + theme_acbs() + geom_vline(xintercept=(1), col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) + theme(axis.title.x = element_text(size=9)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))
# make density for x.1.p
d1 <- as.data.frame(out.inla$marginals.fixed$x.p.1.inla)
d2 <- as.data.frame(density(jags.df$bCov1P)[c(1,2)])
col1 <- 2
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p6 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab("logit(p) covariate 1")  + ylab(NULL) + theme_acbs() +
  geom_vline(xintercept=-2, col="gray20", lty=1)+
  geom_line(data=d7, aes(y=y, x=x)) + theme(axis.title.x = element_text(size=9)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))
# make density for x.3.p
d1 <- as.data.frame(out.inla$marginals.fixed$x.p.3.inla)
d2 <- as.data.frame(density(jags.df$bCov3P)[c(1,2)])
col1 <- 3
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p7 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab("logit(p) covariate 3")  + ylab(NULL) + theme_acbs() +
  geom_vline(xintercept=1, col="gray20", lty=1)+
  geom_line(data=d7, aes(y=y, x=x)) + theme(axis.title.x = element_text(size=9)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))
# make density for overdisp
d1 <- as.data.frame(out.inla$marginals.hyperpar$`overdispersion for NMix observations`)
d1$x <- 1 / d1$x; d1$y <- d1$y / 8
d2 <- as.data.frame(density(jags.df$overDisEst)[c(1,2)])
col1 <- 8
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- exp(c(d3,d4,d5))
d7 <- data.frame(y=c(0,0,0), x=d6)
p8 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab("Overdispersion parameter") + ylab(NULL) + theme_acbs() +
  geom_vline(xintercept=3, col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) + theme(axis.title.x = element_text(size=9)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))
# plot all
png("fig1.png", width = 6, height = 5, units = 'in', res = 600)
multiplot(p1,p2,p3,p4,p5,p6,p7,p8,cols=3)
dev.off()
# ##############################################################################





# simulation code, example 2 ###################################################
sim.fun <- function(){
  # same inputs as before for most parameters
  n.sites = 72; n.surveys = 3; n.years = 9
  lam.b0 = 2.0; lam.b1 = 2.0; lam.b2 = -3.0; lam.b4 = 1.0
  p.b0 = 1.0; p.b1 = -2.0
  disp.size = 3.0
  # now vary the effect size for detection covariate 3
  p.b3 = runif(1, -3.0, 3.0)
  # keep track of input values for later analysis
  real.vals <- c(p.b0, p.b1, p.b3, lam.b0, lam.b1, lam.b2, lam.b4, disp.size)
  # set jags run parameters
  nc <- 3; na<- 500; nb <- 100; ni <- 1000; nt <- 1
  # parameters to monitor
  params <- c("intP", "bCov1P", "bCov3P", "intLam", "bCov1Lam","bCov2Lam",
            "bYr", "overDisEst")
  # make data
  sim.data <- data4sim(n.sites = n.sites, n.surveys = n.surveys,
                       n.years = n.years,
                       lam.b0 = lam.b0, lam.b1 = lam.b1, lam.b2 = lam.b2,
                       lam.b4 = lam.b4,
                       p.b0 = p.b0, p.b1 = p.b1, p.b3 = p.b3,
                       disp.size = disp.size)
  # bundle jags with site-survey-year x.p.3
  jags.data.big <- list(y = sim.data$y2, x.lam.1 = sim.data$x.lam.1,
             x.lam.2 = sim.data$x.lam.2, yr = sim.data$yr,
             x.p.1 = sim.data$x.p.1, x.p.3 = sim.data$x.p.3,
             nSites = sim.data$n.sites, nSurveys = sim.data$n.surveys,
             nYears = sim.data$n.years)
  # initial values for first jags run
  N.init <- sim.data$y2 # initial count values
  N.init[is.na(N.init)] <- 1 # clean up NA's
  N.init <- apply(N.init, c(1, 3), max) + 1 # zero values cause trouble
  inits <- function() list(N = N.init, intLam = rnorm(1, 0, 0.01),
                     intP = rnorm(1, 0, 0.01), bCov1P = rnorm(1, 0, 0.01),
                     bCov2Lam = rnorm(1,0,0.01), bCov1Lam = rnorm(1, 0, 0.01),
                     bCov3P = rnorm(1, 0, 0.01), bYr = rnorm(1, 0, 0.01),
                     overDisEst = runif(1, 0.5, 2.5))
  # new jags model with expanded x.p.3
  jags.model.string.expand <- "
  model {
    # priors
    intP ~ dnorm(0, 0.01)       # detection intercept
    bCov1P ~ dnorm(0, 0.01)     # detection cov 1 effect
    bCov3P ~ dnorm(0, 0.01)     # detection cov 3 effect
    intLam ~ dnorm(0, 0.01)     # lambda intercept
    bCov1Lam ~ dnorm(0, 0.01)   # lambda cov 1 effect
    bCov2Lam ~ dnorm(0, 0.01)   # lambda cov 2 effect
    bYr ~ dnorm(0, 0.01)        # year effect
    overDisEst ~ dunif(0, 5)    # overdispersion size
    # abundance component
    for (k in 1:nYears){
      for (i in 1:nSites){
        N[i, k] ~ dnegbin(prob[i, k], overDisEst) # negative binomial
        prob[i, k] <- overDisEst / (overDisEst + lambda[i, k]) # overdispersion
        log(lambda[i, k]) <- intLam + (bCov1Lam * x.lam.1[i]) +
                              (bCov2Lam * x.lam.2[i]) + (bYr * yr[k])
    # detection component, note new i, j, k, subscript for x.p.3
        for (j in 1:nSurveys){
          y[i, j, k] ~ dbin(p[i,j,k], N[i,k])
          p[i, j, k] <- exp(lp[i,j,k]) / (1 + exp(lp[i,j,k]))
          lp[i, j, k] <- intP + (bCov1P * x.p.1[i]) + (bCov3P * x.p.3[i, j, k])
        } # close j loop
      } # close i loop
    } # close k loop
  } # close model loop
  "
  # call jags for site-survey-year x.p.3
  out.jags.big <- run.jags(model=jags.model.string.expand, data=jags.data.big,
                           monitor=params, n.chains=nc, inits=inits,
                           burnin=nb, adapt=na, sample=ni,
                           thin=nt, modules="glm on", method="parallel")
  # collect results
  bo <- as.data.frame(round(summary(out.jags.big),2))[,c(1:5,9,11)]
  bo$simtype <- "big"; bo$p.b3 <- p.b3; bo$par <- row.names(bo)
  bo$true.vals <- real.vals
  # bundle jags data with site-year x.p.3
  jags.data.small <- list(y = sim.data$y2, x.lam.1 = sim.data$x.lam.1,
             x.lam.2 = sim.data$x.lam.2, yr = sim.data$yr,
             x.p.1 = sim.data$x.p.1, x.p.3 = sim.data$x.p.3.arr,
             nSites = sim.data$n.sites, nSurveys = sim.data$n.surveys,
             nYears = sim.data$n.years)
  # call jags for site-year x.p.3
  out.jags.small <- run.jags(model=jags.model.string.expand,
                             data=jags.data.small,
                             monitor=params, n.chains=nc, inits=inits,
                             burnin=nb, adapt=na, sample=ni,
                             thin=nt, modules="glm on", method="parallel")
  so <- as.data.frame(round(summary(out.jags.small),2))[,c(1:5,9,11)]
  so$simtype <- "small"; so$p.b3 <- p.b3; so$par <- row.names(so)
  so$true.vals <- real.vals
  # combine results
  all.out <- rbind(bo,so)
  all.out$diffs <- all.out$Median - all.out$true.vals
  all.out$range <- all.out$Lower95 - all.out$Upper95
  row.names(all.out) <- NULL
  return(all.out)
}
# first round
ptm <- proc.time()
sim.out <- sim.fun()
# subsequent rounds
nsims <- 49
for(w in 2:nsims){
  iter.out <- sim.fun()
  sim.out <- rbind(sim.out, iter.out)
  cat(paste("\n\n\n this is iteration", w, "\n\n\n"))
}
proc.time() - ptm
round(simu.time.1 <- proc.time() - ptm, 2)[3]
# ##############################################################################





# explore results, example 2 ###################################################
par_names <- c(
  'intLam'="log(lambda) intercept",
  'bCov1Lam'="log(lambda) covariate 1",
  'bCov2Lam'="log(lambda) covariate 2",
  'bYr'="log(lambda) year covariate",
  'intP'="logit(p) intercept",
  'bCov1P'="logit(p) covariate 1",
  'bCov3P'="logit(p) covariate 3",
  'overDisEst'="Overdispersion parameter"
)
to.string <- as_labeller(par_names)
# plot
png("fig2.png", width = 6, height = 5, units = 'in', res = 600)
ggplot(data=sim.out, aes(x=p.b3, y=diffs, colour=simtype, linetype=simtype)) +
  geom_point(pch=1, size=1.5) +
  facet_wrap(~factor(par), scales="fixed",labeller=to.string) +
  xlab("Slope coefficient for logit(p) covariate 3") +
  ylab("Difference between posterior mean and true parameter value") +
  geom_smooth(span=3, se=F, size=0.5) +
  theme_acbs() + scale_color_manual(values=c("black","gray50")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        legend.position="none") + scale_linetype_manual(values=c(1,1))
dev.off()
# ##############################################################################





# unmarked analysis, example 3 #################################################
# get real mallard data from unmarked package
data(mallard)
# format unmarked data
mallard.y.unmk <- mallard.y
mallard.site.unmk <- mallard.site
mallard.obs.unmk <- mallard.obs
mallardUMF <- unmarkedFramePCount(mallard.y.unmk, siteCovs = mallard.site.unmk,
                                  obsCovs = mallard.obs.unmk)
# run unmarked model
mallard.out.unmk <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest,
                           mixture = "NB", mallardUMF, K=30)
summary(mallard.out.unmk)
# ##############################################################################





# inla analysis, example 3 #####################################################
mallard.y.inla <- mallard.y
mallard.length.inla <- mallard.site[,2]
mallard.elev.inla <- mallard.site[,1]
mallard.forest.inla <- mallard.site[,3]
mallard.ivel.inla <- rowMeans(mallard.obs$ivel, na.rm=T) # average per site
mallard.ivel.inla[is.na(mallard.ivel.inla)] <- mean(mallard.ivel.inla, na.rm=T)
mallard.date.inla <- rowMeans(mallard.obs$date, na.rm=T) # average per site
mallard.date2.inla <- mallard.date.inla^2
# make inla.mdata object
counts.and.count.covs <- inla.mdata(mallard.y.inla, 1, mallard.length.inla,
                                    mallard.elev.inla, mallard.forest.inla)
# run inla model
mallard.out.inla <- inla(counts.and.count.covs ~ 1 + mallard.ivel.inla +
                           mallard.date.inla + mallard.date2.inla,
         data = list(counts.and.count.covs = counts.and.count.covs,
                     mallard.ivel.inla = mallard.ivel.inla,
                     mallard.date.inla = mallard.date.inla,
                     mallard.date2.inla = mallard.date2.inla),
         family = "nmixnb",
         control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.01,
                              prec.intercept = 0.01),
         control.family = list(hyper = list(theta1 = list(param = c(0, 0.01)),
                                          theta2 = list(param = c(0, 0.01)),
                                          theta3 = list(param = c(0, 0.01)),
                                          theta4 = list(param = c(0, 0.01)),
                                          theta5 = list(prior = "flat",
                                                        param = numeric()))))
# ##############################################################################





# explore results, example 3 ###################################################
# get unmk out
e1 <- as.matrix(mallard.out.unmk@estimates@estimates$det@estimates)
ci1 <- as.matrix(confint(mallard.out.unmk, type="det"))
e2 <- as.matrix(mallard.out.unmk@estimates@estimates$state@estimates)
ci2 <- as.matrix(confint(mallard.out.unmk, type="state"))
e3 <- exp(as.matrix(mallard.out.unmk@estimates@estimates$alpha@estimates))
ci3 <- exp(as.matrix(confint(mallard.out.unmk, type="alpha")))
e4 <- rbind(e1,e2,e3)
ci4 <- rbind(ci1,ci2,ci3)
unmk.out.tab.2 <- cbind(e4,ci4); colnames(unmk.out.tab.2) <- c("X1", "X2", "X3")
# get inla out
in1 <- mallard.out.inla$summary.fixed[,c(4, 3, 5)]
in2 <- mallard.out.inla$summary.hyperpar[,c(4, 3, 5)][-5,]
in3 <- 1/mallard.out.inla$summary.hyperpar[,c(4, 3, 5)][5,]
in3 <- in3[,c(1,3,2)]
inla.out.tab.2 <- rbind(in1,in2,in3); colnames(inla.out.tab.2) <- c("X1", "X2",
                                                                    "X3")
# combine out
all.out.tab.2 <- round(rbind(unmk.out.tab.2, inla.out.tab.2),2)
colnames(all.out.tab.2) <- c("Estimate","Lower", "Upper")
all.out.tab.2$Parameter <- factor(rep(c("logit(p) intercept",
                                        "Survey intensity",
                                        "Survey date (SD)", "SD squared",
                                        "log(lambda) intercept",
                                        "Transect length", "Transect elevation",
                                        "Forest cover",
                                        "Overdispersion"), 2))
all.out.tab.2$Technique <- factor(rep(c("unmarked", "R-INLA"), each=9))
rownames(all.out.tab.2) <- NULL
pd <- position_dodge(0.2)
# plot
png("fig3.png", width = 6, height = 3, units = 'in', res = 600)
ggplot(data=all.out.tab.2, aes(y=Estimate, x=Parameter, col=Technique)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position=pd) +
  theme_acbs() + scale_color_manual(values=c("black","gray50")) +
  coord_flip() + geom_hline(yintercept=0, col="gray30") +
  ylab("Estimated value (95% CI or CrI)") +
  xlab("Model parameter")
dev.off()
# ##############################################################################





# wrap up ######################################################################
save.image("~/GitHub/Quantitative_Metrics/INLA_NMix/meehan_et_al_inla_nmix.RData")
# ##############################################################################


