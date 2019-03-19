### Gibbs sampling for mean and variance estimation of unknown normal distribution

## Gibbs sampling with uninformative priors
set.seed(1237)                     # set randomizer seed
m=50000                            # set num iterations
MU=numeric(m);  THETA=numeric(m);  # sampled values
THETA[1]=1;                        # inital value
n=41; x.bar=9.6; x.var=2.73^2;     # data
mu.0=0; th.0=400;                  # mu prors
alp.0=1/2; kap.0=1/5;              # theta priors

for (i in 2:m)
{
  th.up=1/(n/THETA[i-1]+1/th.0)
  mu.up=(n*x.bar/THETA[i-1]+mu.0/th.0)*th.up
  MU[i]=rnorm(1,mu.up,sqrt(th.up))
  
  alp.up=n/2+alp.0
  kap.up=kap.0+((n-1)*x.var+n*(x.bar-MU[i])^2)/2
  THETA[i]=1/rgamma(1,alp.up,kap.up)
}

# Bayesian point and probability interval estimates
aft.brn=(m/2+1):m                                         # discard first half of simulations
mean(MU[aft.brn])                                         # point estimate mu
bi.MU=quantile(MU[aft.brn],c(.025,.975)); bi.MU           # 95% confidence interval mu
mean(THETA[aft.brn])                                      # point estimate theta
bi.THETA=quantile(THETA[aft.brn],c(.025,.975)); bi.THETA  # 95% confidence interval theta
SIGMA=sqrt(THETA)
mean(SIGMA[aft.brn])                                      # point estimate of sigma
bi.SIGMA=sqrt(bi.THETA); bi.SIGMA                         # 95% confidence interval sigma

# Plot diagnostic graphs for correlation and distribution of samples
par(mfrow=c(2,2))
  plot(aft.brn,MU[aft.brn],type="l")
  plot(aft.brn,SIGMA[aft.brn],type="l")
  hist(MU[aft.brn],prob=T); abline(v=bi.MU, col="red")
  hist(SIGMA[aft.brn],prob=T); abline(v=bi.SIGMA, col="red")
par(mfrow=c(1,1))

# Plot diagnostic graphs for ACF and cumulative average of samples
par(mfrow=c(2,2))
  acf(MU); acf(THETA);
  plot(1:m,cumsum(MU)/(1:m),type="l", ylim=c(9.4,9.8))
  plot(1:m,cumsum(SIGMA)/(1:m),type="l",ylim=c(2.6,2.8))
par(mfrow=c(1,1))


## Gibbs sampling with informative priors
set.seed(1237)                     # set randomizer seed
m=50000                            # num iterations
MU=numeric(m);  THETA=numeric(m);  # sampled values
THETA[1]=1;                        # inital value
n=41; x.bar=9.6; x.var=2.73^2;     # data
mu.0=10; th.0=1;                   # mu priors
alp.0=20; kap.0=200;               # theta priors

for (i in 2:m)
{
  th.up=1/(n/THETA[i-1]+1/th.0)
  mu.up=(n*x.bar/THETA[i-1]+mu.0/th.0)*th.up
  MU[i]=rnorm(1,mu.up,sqrt(th.up))
  
  alp.up=n/2+alp.0
  kap.up=kap.0+((n-1)*x.var+n*(x.bar-MU[i])^2)/2
  THETA[i]=1/rgamma(1,alp.up,kap.up)
}

# Bayesian point and probability interval estimates
aft.brn=(m/2+1):m                                         # discard first half of simulations
mean(MU[aft.brn])                                         # point estimate mu
bi.MU=quantile(MU[aft.brn],c(.025,.975)); bi.MU           # 95% confidence interval mu
mean(THETA[aft.brn])                                      # point estimate theta
bi.THETA=quantile(THETA[aft.brn],c(.025,.975)); bi.THETA  # 95% confidence interval theta
SIGMA=sqrt(THETA)
mean(SIGMA[aft.brn])                                      # point estimate of sigma
bi.SIGMA=sqrt(bi.THETA); bi.SIGMA                         # 95% confidence interval sigma

# Plot diagnostic graphs for correlation and distribution of samples
par(mfrow=c(2,2))
plot(aft.brn,MU[aft.brn],type="l")
plot(aft.brn,SIGMA[aft.brn],type="l")
hist(MU[aft.brn],prob=T); abline(v=bi.MU, col="red")
hist(SIGMA[aft.brn],prob=T); abline(v=bi.SIGMA, col="red")
par(mfrow=c(1,1))

# Plot diagnostic graphs for ACF and cumulative average of samples
par(mfrow=c(2,2))
acf(MU); acf(THETA);
plot(1:m,cumsum(MU)/(1:m),type="l", ylim=c(9.4,9.8))
plot(1:m,cumsum(SIGMA)/(1:m),type="l",ylim=c(2.6,2.8))
par(mfrow=c(1,1))



## Frequentist parameter estimation
x.bar                                       # point estimate mu
x.bar+c(-1,1)*qt(.975,40)*sqrt(x.var/n)     # 95% confidence interval mu
sqrt(x.var)                                 # point estimate sigma
sqrt((n-1)*x.var/qchisq(c(.975,.025),n-1))  # 95% confidence interval sigma

