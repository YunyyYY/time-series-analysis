library(pomp)
library(ggplot2)
library(plyr)
library(reshape2)
library(foreach)
library(doParallel)
registerDoParallel()
cores <- 20  # The number of cores on this machine 
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

theme_set(theme_bw())

set.seed(987654321L)

#### load part 1 model, set parameters ####
bsflu_statenames <- c("S","I","R1","R2")
bsflu_paramnames <- c("Beta","mu_IR","rho","mu_R1","mu_R2")

bsflu_dmeasure <- "
  lik = dpois(B,rho*R1+1e-6,give_log);
"
bsflu_rmeasure <- "
  B = rpois(rho*R1+1e-10);
  C = rpois(rho*R2);
"

bsflu_rprocess <- "
  double t1 = rbinom(S,1-exp(-Beta*I*dt));
  double t2 = rbinom(I,1-exp(-dt*mu_IR));
  double t3 = rbinom(R1,1-exp(-dt*mu_R1));
  double t4 = rbinom(R2,1-exp(-dt*mu_R2));
  S -= t1;
  I += t1 - t2;
  R1 += t2 - t3;
  R2 += t3 - t4;
"

bsflu_rinit <- "
 S=762;
 I=1;
 R1=0;
 R2=0;
"

bsflu_data <- read.table("bsflu_data.txt")
bsflu2 <- pomp(
  data=bsflu_data,
  times="day",
  t0=0,
  rprocess=euler(
    step.fun=Csnippet(bsflu_rprocess),
    delta.t=1/12
  ),
  rmeasure=Csnippet(bsflu_rmeasure),
  dmeasure=Csnippet(bsflu_dmeasure),
  partrans=parameter_trans(log=c("Beta","mu_IR","mu_R1","mu_R2"),logit="rho"),
  statenames=bsflu_statenames,
  paramnames=bsflu_paramnames,
  rinit=Csnippet(bsflu_rinit)
)

#### set run_level ####

run_level <- 1
switch(run_level, {  # Np is number of particles
  bsflu_Np=2000; bsflu_Nmif=100; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10
},{
  bsflu_Np=20000; bsflu_Nmif=200; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10
},{
  bsflu_Np=60000; bsflu_Nmif=300; bsflu_Neval=10; bsflu_Nglobal=100; bsflu_Nlocal=20}
)

bsflu_params <- data.matrix(
  read.table("mif_bsflu_params.csv",
             row.names=NULL,header=TRUE))
which_mle <- which.max(bsflu_params[,"logLik"])
bsflu_mle <- bsflu_params[which_mle,][bsflu_paramnames]
bsflu_fixed_params <- c(mu_R1=1/(sum(bsflu_data$B)/512),
                        mu_R2=1/(sum(bsflu_data$C)/512))

#### get `mifs_local` for global likelihood search use ####
# local search of the likelihood surface
bsflu_rw.sd <- 0.02;
bsflu_cooling.fraction.50 <- 0.5
stew(file=sprintf("out/local_search-%d.rda",run_level),{
  t_local <- system.time({
    mifs_local <- foreach(i=1:bsflu_Nlocal, .packages='pomp', .combine=c) %dopar% {
      mif2(bsflu2,
           params=bsflu_mle,
           Np=bsflu_Np,
           Nmif=bsflu_Nmif,
           cooling.fraction.50=bsflu_cooling.fraction.50,
           rw.sd=rw.sd(Beta=bsflu_rw.sd, mu_IR=bsflu_rw.sd, rho=bsflu_rw.sd)
      )
    }
  })
},seed=900242057,kind="L'Ecuyer")
stew(file=sprintf("out/lik_local-%d.rda",run_level),{
  t_local_eval <- system.time({
    liks_local <- foreach(i=1:bsflu_Nlocal,.combine=rbind) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu2,params=coef(mifs_local[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")
results_local <- data.frame(logLik=liks_local[,1],
                            logLik_se=liks_local[,2],
                            t(sapply(mifs_local,coef)))  # t() is transpose


#### A global likelihood search using random starting values ####
bsflu_box <- rbind(
  Beta=c(0.001,0.01),
  mu_IR=c(1.5,2),
  rho = c(0.8,0.85)
)
stew(file=sprintf("out/box_eval-%d.rda",run_level),{
  t_global <- system.time({
    mifs_global <- foreach(i=1:bsflu_Nglobal,.combine=c) %dopar% {
      mif2(
        mifs_local[[1]],  # object of class `mif2d_pomp`
        params=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),bsflu_fixed_params)
      )}
  })
},seed=1270401374,kind="L'Ecuyer")

#### get global ####

stew(file=sprintf("out/box_eval_log-%d.rda",run_level),{
  t_global.3 <- system.time({
    mifs_global.3 <- foreach(i=1:bsflu_Nglobal,.packages='pomp',.combine=c,.options.multicore=mcopts) %dopar% 
      mif2(
        mifs_global[[1]],
        start=c(Beta=exp(runif(1,log(bsflu_box[1,1]),log(bsflu_box[1,2])))
                ,apply(bsflu_box[2:3,],1,function(x) runif(1,x[1],x[2])),bsflu_fixed_params),
        Np=bsflu_Np,
        Nmif=bsflu_Nmif
      )
  })
},seed=931129,kind="L'Ecuyer")

stew(file=sprintf("out/lik_global_eval_log-%d.rda",run_level),{
  t_global_eval.3 <- system.time({
    liks_global.3 <- foreach(i=1:bsflu_Nglobal,.packages='pomp',.combine=rbind,.options.multicore=mcopts) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu2,params=coef(mifs_global.3[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=920804,kind="L'Ecuyer")

results_global.3 <- data.frame(logLik=liks_global.3[,1],logLik_se=liks_global.3[,2],t(sapply(mifs_global.3,coef)))
print(summary(results_global.3$logLik,digits=5))

par(mfrow=c(1,2))
hist(runif(5000,0.001,0.1),ylim=c(0,2600),xlab=expression(beta),
     main="Uniformly drawing \n from orginal scale")
hist(exp(runif(5000,log(0.001),log(0.1))),xlab=expression(beta),
     ylim=c(0,2600),
     main="Uniformly drawing \n from log scale")


