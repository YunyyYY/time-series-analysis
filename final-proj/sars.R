# run but not shown code
options(warn=-1)  # disable the warnings
library(knitr)
library(dplyr)
library(plyr)
library(reshape2)
library(pomp)
library(ggplot2)
theme_set(theme_bw())

library(doParallel)
registerDoParallel()
library(doRNG)
registerDoRNG(3899882)

df <- read.table("sars_HK.csv",header=TRUE, sep=",",colClasses=c("Date",rep("numeric",4)))
# typeof(df$Date)

df$day = seq(1, length(df$Date))

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt)); 
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")
sir_skel <- Csnippet("
  double dN_SIdt = Beta*S*I/N;
  double dN_IRdt = mu_IR*I;
  DS = -dN_SIdt;
  DI = dN_SIdt - dN_IRdt;
  DR = dN_IRdt;
  DH = dN_IRdt;
")
# At day zero, we’ll assume that I = 50 and R = 0, 
# but we don’t know the actual population, 
# so we treat N as a parameter to be estimated and let S(0) = N −50. 

sir_rinit <- Csnippet("
  S = nearbyint(N)-50;
  I = 50;
  R = 0;
  H = 0;
")

sir_dmeas <- Csnippet("lik = dbinom(Current,H,rho,give_log);")
sir_rmeas <- Csnippet("Current = rbinom(H,rho);")

sir <- pomp(
  subset(df,select=c(day,Current)),
  time="day",t0=0, rprocess=euler(sir_step,delta.t=1/6),
  rmeasure=sir_rmeas, dmeasure=sir_dmeas, rinit=sir_rinit,
  skeleton = vectorfield(sir_skel),
  paramnames=c("Beta","mu_IR","N","rho"), statenames=c("S","I","R","H"),
  accumvars="H",
  partrans=parameter_trans(log=c("Beta","mu_IR","N"),logit="rho")
) 

sars_rprocess <- "
double dN_SI = rbinom(S,1-exp(-Beta*I*dt)); 
double dN_IR1 = rbinom(I,1-exp(-dt*mu_IR)); 
double dN_R1R2 = rbinom(R1,1-exp(-dt*mu_R1)); 
double dN_R2R3 = rbinom(R2,1-exp(-dt*mu_R2)); 
S -= dN_SI;
I += dN_SI - dN_IR1;
R1 += dN_IR1 - dN_R1R2;
R2 += dN_R1R2 - dN_R2R3;
"
sars_dmeasure <- "
lik = dpois(Current,rho*R1+1e-10,give_log);
"
sars_rmeasure <- "
Current = rpois(rho*R1+1e-10);
"
# 1e-10 tolerance value prevents the code from crashing when all particles have R1=0.
sars_rinit <- " 
S=100000-50;
I=50;
R1=0;
R2=0; "
sars_statenames <- c("S","I","R1","R2")
sars_paramnames <- c("Beta","mu_IR","rho","mu_R1","mu_R2")

sars2 <- pomp( 
  data=subset(df,select=c(day,Current)), times="day", t0=0,
  rprocess=euler(step.fun=Csnippet(sars_rprocess),delta.t=1/12),
  rmeasure=Csnippet(sars_rmeasure), dmeasure=Csnippet(sars_dmeasure),
  partrans=parameter_trans(log=c("Beta","mu_IR","mu_R1","mu_R2"), logit="rho"),
  statenames=sars_statenames, paramnames=sars_paramnames, rinit=Csnippet(sars_rinit)
)

run_level <- 2
switch(run_level, {
  sars_Np=100; sars_Nmif=10; sars_Neval=10; sars_Nglobal=10; sars_Nlocal=10
},{
  sars_Np=20000; sars_Nmif=100; sars_Neval=10; sars_Nglobal=10; sars_Nlocal=10
},{
  sars_Np=60000; sars_Nmif=300; sars_Neval=10; sars_Nglobal=100; sars_Nlocal=20}
)
sars_params <- data.matrix( read.table("mif_sars_params.csv", row.names=NULL,header=TRUE))
which_mle <- which.max(sars_params[,"logLik"]) 
sars_mle <- sars_params[which_mle,][sars_paramnames]

sars_fixed_params <- c(mu_R1=0.99, mu_R2=0.8) # mu_confined, mu_recovered

sars_rw.sd <- 0.02; 
sars_cooling.fraction.50 <- 0.5 
stew(file=sprintf("local_search-%d.rda",run_level),{
  t_local <- system.time({
    mifs_local <- foreach(i=1:sars_Nlocal, .packages='pomp', .combine=c) %dopar% {
      mif2(sars2, params=sars_mle, Np=sars_Np,
           Nmif=sars_Nmif, cooling.fraction.50=sars_cooling.fraction.50,
           rw.sd=rw.sd(Beta=sars_rw.sd, mu_IR=sars_rw.sd, rho=sars_rw.sd))
    }
  })
},seed=900242057,kind="L'Ecuyer")

stew(file=sprintf("lik_local-%d.rda",run_level),{
  t_local_eval <- system.time({
    liks_local <- foreach(i=1:sars_Nlocal,.combine=rbind) %dopar% {
      evals <- replicate(sars_Neval, 
                         logLik(pfilter(sars2,params=coef(mifs_local[[i]]),Np=sars_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")

results_local <- data.frame(logLik=liks_local[,1], logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))

sars_box <- rbind(Beta=c(0,1), mu_IR=c(0,1),rho = c(0.7,1))

stew(file=sprintf("box_eval-%d.rda",run_level),{ 
  t_global <- system.time({
    mifs_global <- foreach(i=1:sars_Nglobal,.combine=c) %dopar% { 
      mif2(mifs_local[[1]], params=c(
        apply(sars_box,1,function(x)runif(1,x[1],x[2])),
        sars_fixed_params) 
      )}
  }) 
},seed=1270401374,kind="L'Ecuyer")

stew(file=sprintf("lik_global_eval-%d.rda",run_level),{ 
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:sars_Nglobal, .combine=rbind) %dopar% {
      evals <- replicate(sars_Neval, 
                         logLik(pfilter(sars2,params=coef(mifs_global[[i]]),Np=sars_Np))) 
      logmeanexp(evals, se=TRUE)
    } })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(
  logLik=liks_global[,1], logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))


if (run_level>1) 
  write.table(rbind(results_local,results_global),
              file="mif_sars_params.csv", append=TRUE,
              col.names=FALSE,row.names=FALSE)

It=20
nprof=20
profile.box <- profileDesign(  
  Beta=seq(0.1,1.5,length.out=It),
  lower=c(mu_IR=0,rho=0.5),
  upper=c(mu_IR=1.5,rho=1.0),
  nprof=nprof
)

sars_rw.sd <- 0.02;
stew(file=sprintf("profile_phi-%d.rda",It),{
  t_pro <- system.time({
    prof.llh <- foreach(i=1:dim(profile.box)[1], .packages='pomp', .combine=rbind) %dopar%{
      mif2(
        mifs_local[[1]],
        params=c(unlist(profile.box[i,]),sars_fixed_params),
        Np=sars_Np,
        Nmif=sars_Nmif,
        rw.sd=rw.sd(
          mu_IR=sars_rw.sd,
          rho=sars_rw.sd
        )
      ) -> mifs_pro
      evals = replicate(10, logLik(pfilter(mifs_pro,Np=sars_Np)))
      ll=logmeanexp(evals, se=TRUE)
      data.frame(as.list(coef(mifs_pro)),
                 loglik = ll[1],
                 loglik.se = ll[2])
    }
  })
},seed=5556129,kind="L'Ecuyer")






















