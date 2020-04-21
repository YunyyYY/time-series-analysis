require(plyr)
require(ggplot2)
require(dplyr)
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

#### profile likelihood ####

It=20
nprof=20
profile.box <- profileDesign(  
  Beta=exp(seq(log(0.0035),log(0.0045),length.out=It)),
  lower=c(mu_IR=1.5,rho=0.8),
  upper=c(mu_IR=2.0,rho=0.9),
  nprof=nprof
)

bsflu_rw.sd <- 0.02;
stew(file=sprintf("out/profile_phi-%d.rda",It),{
  t_pro <- system.time({
    prof.llh <- foreach(i=1:dim(profile.box)[1], .packages='pomp', .combine=rbind) %dopar%{
      mif2(
        mifs_local[[1]],
        params=c(unlist(profile.box[i,]),bsflu_fixed_params),
        Np=bsflu_Np,
        Nmif=bsflu_Nmif,
        rw.sd=rw.sd(
          mu_IR=bsflu_rw.sd,
          rho=bsflu_rw.sd
        )
      ) -> mifs_pro
      evals = replicate(10, logLik(pfilter(mifs_pro,Np=bsflu_Np)))
      ll=logmeanexp(evals, se=TRUE)
      data.frame(as.list(coef(mifs_pro)),
                 loglik = ll[1],
                 loglik.se = ll[2])
    }
  })
},seed=5556129,kind="L'Ecuyer")

prof.llh %<>%
  mutate(Beta=exp(signif(log(Beta),5))) %>%
  ddply(~Beta,subset,rank(-loglik)<=1)

a=max(prof.llh$loglik)
b=a-1.92
CI=which(prof.llh$loglik>=b)
c=prof.llh$Beta[min(CI)]
d=prof.llh$Beta[max(CI)]

# png(filename="q8-3-prorfile.png")
as.data.frame(prof.llh) %>%
  ggplot(aes(x=Beta,y=loglik))+
  geom_point()+
  geom_smooth(method="loess")+
  geom_hline(aes(yintercept=a),linetype="dashed")+
  geom_hline(aes(yintercept=b),linetype="dashed")+
  geom_vline(aes(xintercept=c),linetype="dashed")+
  geom_vline(aes(xintercept=d),linetype="dashed") -> img
plot(img)
# dev.off()

print(c(lower=c,upper=d))
