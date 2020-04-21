library(pomp)
library(ggplot2)
library(plyr)
library(reshape2)
library(foreach)
library(doParallel)
registerDoParallel()
theme_set(theme_bw())

set.seed(987654321L)

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
stew(file=sprintf("out/local_search-%d.rda",run_level, 1),{
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
stew(file=sprintf("out/lik_local-%d.rda",run_level, 1),{
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

print("Local search results:")
print(summary(results_local$logLik,digits=5))
# png(filename="out/q8-1-local-pair.png")
pairs(~logLik+Beta+mu_IR+rho, data=subset(results_local,logLik>max(logLik)-50))
# dev.off()

#### A global likelihood search using random starting values ####
bsflu_box <- rbind(
  Beta=c(0.004,0.005),
  mu_IR=c(1.5,2),
  rho = c(0.8,0.85)
)
stew(file=sprintf("out/box_eval-%d.rda",run_level, 1),{
  t_global <- system.time({
    mifs_global <- foreach(i=1:bsflu_Nglobal,.combine=c) %dopar% {
      mif2(
        mifs_local[[1]],  # object of class `mif2d_pomp`
        params=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),bsflu_fixed_params)
      )}
  })
},seed=1270401374,kind="L'Ecuyer")

stew(file=sprintf("out/lik_global_eval-%d.rda",run_level, 1),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:bsflu_Nglobal, .combine=rbind) %dopar% {
      evals <- replicate(bsflu_Neval,
                         logLik(pfilter(bsflu2,
                                        params=coef(mifs_global[[i]]),
                                        Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(logLik=liks_global[,1],
                             logLik_se=liks_global[,2],
                             t(sapply(mifs_global,coef)))

if (run_level > 1)
  write.table(rbind(results_local,results_global),
              file="out/mif_bsflu_params.csv",
              append=TRUE,col.names=TRUE,row.names=TRUE)

print("Global search results")
print(summary(results_global$logLik,digits=5))
# png(filename="out/q8-1-global-pair.png")
pairs(~logLik+Beta+mu_IR+rho, data=subset(results_global,logLik>max(logLik)-250))
# dev.off()
# png(filename="out/q8-1.png")
plot(mifs_global)
# dev.off()
