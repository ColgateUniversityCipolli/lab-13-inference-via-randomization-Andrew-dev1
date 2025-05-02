## importing libraries
library(tidyverse)
library(pwr)
library(e1071)
library(boot.pval)
library(boot)

## 1a
zebra.finches.dat <- read_csv("zebrafinches.csv")

further.vec <- zebra.finches.dat$further
skew <- skewness(further.vec)
n <- length(further.vec)

further.test <- t.test(x=further.vec, mu = 0, alternative = "less")
t.far <- as.numeric(further.test$statistic) #finding the t value 
fz <- dnorm(t.far)
Fz <- pnorm(t.far)
(error = (skew / sqrt(n)) * ((2 * (t.far^2) + 1) / 6) * fz)

prob <- error + Fz



t.vecs <- seq(-10,10,0.1)
results <- tibble(
  t = t.vecs,
  fz = dnorm(t.vecs), 
  error = (skew / sqrt(n)) * ((2 * (t.vecs^2) + 1) / 6) * dnorm(t.vecs)
)

ggplot(results) +
  theme_bw()+
  geom_line(aes(x = t, y = error))

# 1c 



alpha = 0.05
# when close to rejection zone, there is high error rate
# but t of alpha, n-1 is approx qnorm of alpha
t.val <- qnorm(alpha) 
pdf <- dnorm(t.val)
(target.size = ((skew / (6*(0.1*alpha))) * (2 * (t.val^2) + 1) * pdf)^2)

###########################
#  Part 2:Boot strapping  #
###########################
R <- 10000
mu0 = 0 
fur.sd <- sd(zebra.finches.dat$further)
close.sd <- sd(zebra.finches.dat$closer)
diff.sd <- sd(zebra.finches.dat$diff)

resample <- tibble(close.t=numeric(R),
                   far.t=numeric(R),
                   diff.t=numeric(R),
                   close.mean =numeric(R),
                   far.mean =numeric(R),
                   diff.mean =numeric(R)
)

# resamples.null.closer <- tibble()
# resamples.null.further <- tibble()
# resamples.null.diff <- tibble()
resamples.shifted <- tibble()


for(i in 1:R){
  close.sample <- sample(x=zebra.finches.dat$closer,
                         size=n,
                         replace=T)
  far.sample <- sample(x=zebra.finches.dat$further,
                       size=n,
                       replace=T)
  diff.sample <- sample(x=zebra.finches.dat$diff,
                        size=n,
                        replace=T)

  resample$close.t[i] = (mean(close.sample)-mu0)/(close.sd/sqrt(n))
  resample$far.t[i] = (mean(far.sample)-mu0)/(fur.sd/sqrt(n))
  resample$diff.t[i] = (mean(diff.sample)-mu0)/(diff.sd/sqrt(n))
  resample$close.mean[i] = mean(close.sample)
  resample$far.mean[i] = mean(far.sample)
  resample$diff.mean[i] = mean(diff.sample)
  
}
view(resample)


## shifting the means
close.change <- mean(resample$close.t) 
far.change <- mean(resample$far.t) 
diff.change <- mean(resample$diff.t) 

# resamples.shifted <- resample$close.t - change
# mean(resamples.shifted)

resamples.shifted <- resample |>
  mutate(close.shifted = close.t - close.change)  |>
  mutate(far.shifted = far.t - far.change) |>
  mutate(diff.shifted = diff.t - diff.change) |>
  select(c(close.shifted, far.shifted, diff.shifted ))

## checking the sampling distribution
ggplot(resamples.shifted)+
  geom_histogram(aes(x = close.shifted, y =after_stat(density)))+
  theme_bw()
ggplot(resamples.shifted)+
  geom_histogram(aes(x = far.shifted, y =after_stat(density)))+
  theme_bw()
ggplot(resamples.shifted)+
  geom_histogram(aes(x = diff.shifted, y =after_stat(density)))+
  theme_bw()


## 2b 
## p-value = proportion of times data is against the null

# Bootstrap P-Value
boot.p.closer <- mean(resamples.shifted$close.shifted >= close.change)
boot.p.far <- mean(resamples.shifted$far.shifted <= far.change)
low = -diff.change
high = diff.change
p.low <- mean(resamples.shifted$diff.shifted <= low)
p.high <- mean(resamples.shifted$diff.shifted >= high)
boot.p.diff <- p.low +p.high

# T-Test P-Value
test.close <- t.test(zebra.finches.dat$closer,
                     mu = 0,
                     alternative = "greater")
close.ttest.p <- test.close$p.value

test.far <- t.test(zebra.finches.dat$further,
                     mu = 0,
                     alternative = "less")
far.ttest.p <- test.far$p.value

test.diff <- t.test(zebra.finches.dat$diff,
                   mu = 0,
                   alternative = "two.sided")
diff.ttest.p <- test.diff$p.value

(comparisons.p <- tibble(" "= c("close", "far", "diff"), Bootstrapped = 
                        c(boot.p.closer, boot.p.far, boot.p.diff),
                      t.test = c(close.ttest.p, far.ttest.p, diff.ttest.p)))


## Part c
close.5th <- quantile(resamples.shifted$close.shifted, 0.05)
far.5th <- quantile(resamples.shifted$far.shifted, 0.05)
diff.5th <- quantile(resamples.shifted$diff.shifted, 0.05)

close.percentile.t <- qt(0.05, df = n-1)
far.percentile.t <- qt(0.05, df = n-1)
diff.percentile.t <- qt(0.05, df = n-1)

(comparisons.percentiles <- tibble(" "= c("close", "far", "diff"), sampled = 
                         c(close.5th, far.5th, diff.5th),
                       t.val =c(close.percentile.t, far.percentile.t, 
                                 diff.percentile.t)))

## Part d
lower.close <- quantile(resample$close.mean, 0.025)
upper.close <- quantile(resample$close.mean, 0.925)
Bootstrap.CI.close <- c(lower.close, upper.close)

lower.far <- quantile(resample$far.mean, 0.025)
upper.far <- quantile(resample$far.mean, 0.925)
Bootstrap.CI.far <- c(lower.far, upper.far)

lower.diff <- quantile(resample$diff.mean, 0.025)
upper.diff <- quantile(resample$diff.mean, 0.925)
Bootstrap.CI.diff <- c(lower.diff, upper.diff)

# t tests confidence intervals
test.close <- t.test(zebra.finches.dat$closer,
                    mu = 0, conf.level = 0.95,
                    alternative = "two.sided")
close.CI <- test.close$conf.int

test.far <- t.test(zebra.finches.dat$further,
                    mu = 0, conf.level = 0.95,
                    alternative = "two.sided")
far.CI <- test.far$conf.int

test.diff <- t.test(zebra.finches.dat$diff,
                    mu = 0, conf.level = 0.95,
                    alternative = "two.sided")
diff.CI <- test.diff$conf.int

comparisons.CI <- tibble(
  "Condition" = c("close", "far", "diff"),
  "Bootstrapped Lower" = c(lower.close, lower.far, lower.diff),
  "Bootstrapped Upper" = c(upper.close, upper.far, upper.diff),
  "t-test Lower" = c(close.CI[1], far.CI[1], diff.CI[1]),
  "t-test Upper" = c(close.CI[2], far.CI[2], diff.CI[2])
)
view(comparisons.CI)

###########################
#   Part 3:Randomization  #
###########################
mu0 <- 0
R <- 10000
randomized <- tibble(close.means = numeric(R), far.means = numeric(R),
                  diff.means = numeric(R))

shifted.x.close <- zebra.finches.dat$closer - mu0
shifted.x.far <- zebra.finches.dat$further - mu0
shifted.x.diff <- zebra.finches.dat$diff - mu0

for(i in 1:R){
  curr.rand1 <- shifted.x.close *
    sample(x = c(-1, 1),
           size = length(shifted.x.close),
           replace = T)
  curr.rand2 <- shifted.x.far *
    sample(x = c(-1, 1),
           size = length(shifted.x.far),
           replace = T)
  curr.rand3 <- shifted.x.diff *
    sample(x = c(-1, 1),
           size = length(shifted.x.diff),
           replace = T)
  
  randomized$close.means[i] <- mean(curr.rand1)
  randomized$far.means[i] <- mean(curr.rand2)
  randomized$diff.means[i] <- mean(curr.rand3)
}

randomized <- randomized |>
  mutate(close.means = close.means + mu0,
         far.means = far.means + mu0, 
         diff.means = diff.means + mu0 ) # shifting data back

############################################
## 3b
## found observed data 
observed.close.mean <- mean(zebra.finches.dat$closer -mu0)
observed.far.mean <- mean(zebra.finches.dat$further -mu0)
observed.diff.mean <- mean(zebra.finches.dat$diff -mu0)

## p-values for each value type 
rand.p.close <- mean(randomized$close.means >= observed.close.mean)
rand.p.far <- mean(randomized$far.means <= observed.far.mean)

(delta <- abs(mean(zebra.finches.dat$diff) - mu0))
(low <- mu0 - delta) # mirror
(high<- mu0 + delta)   # xbar
rand.p.diff <- mean(randomized$diff.means <= low) +
  mean(randomized$diff.means >= high)

############################################
## 3c
## creating randomized Confidence Intervals 

R <- 1000
mu0.iterate <- 0.0001
starting.point <- mean(zebra.finches.dat$closer)
mu.lower.close <- starting.point

repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zebra.finches.dat$closer - mu.lower.close
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }

  rand <- rand |>
    mutate(xbars = xbars + mu.lower.close) # shifting back
  
  # p-value 
  delta <- abs(mean(zebra.finches.dat$closer) - mu.lower.close)
  low <- mu.lower.close - delta # mirror
  high<- mu.lower.close + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  # delta <- mean(zebra.finches.dat$closer)
  # p.val <- mean(rand$xbars >= delta)

  
  if(p.val < 0.05){
    break
  }else{
    mu.lower.close <- mu.lower.close - mu0.iterate
  }
}


mu.upper.close <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zebra.finches.dat$closer - mu.upper.close
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }

  rand <- rand |>
    mutate(xbars = xbars + mu.upper.close) # shifting back
  
  # p-value 
  delta <- abs(mean(zebra.finches.dat$closer) - mu.upper.close)
  (low <- mu.upper.close - delta) # mirror
  (high<- mu.upper.close + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  # delta <- mean(zebra.finches.dat$closer)
  # p.val <- mean(rand$xbars <= delta)
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper.close <- mu.upper.close + mu0.iterate
  }
}


## CI for further 
starting.point.far <- mean(zebra.finches.dat$further)
mu.lower.far <- starting.point.far

repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zebra.finches.dat$further - mu.lower.far
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  
  rand <- rand |>
    mutate(xbars = xbars + mu.lower.far) # shifting back
  
  # p-value 
  delta <- abs(mean(zebra.finches.dat$further) - mu.lower.far)
  low <- mu.lower.far - delta # mirror
  high<- mu.lower.far + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high)
  # delta <- mean(zebra.finches.dat$further)
  # p.val <- mean(rand$xbars <= delta)
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower.far <- mu.lower.far - mu0.iterate
  }
}


mu.upper.far <- starting.point.far
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zebra.finches.dat$further - mu.upper.far
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  
  rand <- rand |>
    mutate(xbars = xbars + mu.upper.far) # shifting back
  
  # p-value 
  delta <- abs(mean(zebra.finches.dat$further) - mu.upper.far)
  low <- mu.upper.far - delta # mirror
  high<- mu.upper.far + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  # delta <- mean(zebra.finches.dat$further)
  # p.val <- mean(rand$xbars >= delta)
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper.far <- mu.upper.far + mu0.iterate
  }
}


## CI for diff
starting.point.diff <- mean(zebra.finches.dat$diff)
mu.lower.diff <- starting.point.diff

repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zebra.finches.dat$diff - mu.lower.diff
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  
  rand <- rand |>
    mutate(xbars = xbars + mu.lower.diff) # shifting back
  
  # p-value 
  delta <- abs(mean(zebra.finches.dat$diff) - mu.lower.diff)
  low <- mu.lower.diff - delta # mirror
  high<- mu.lower.diff + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  
  if(p.val < 0.05){
    break
  }else{
    mu.lower.diff <- mu.lower.diff - mu0.iterate
  }
}


mu.upper.diff <- starting.point.diff
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- zebra.finches.dat$diff - mu.upper.diff
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  
  rand <- rand |>
    mutate(xbars = xbars + mu.upper.diff) # shifting back
  
  # p-value 
  delta <- abs(mean(zebra.finches.dat$diff) - mu.upper.diff)
  low <- mu.upper.diff - delta # mirror
  high<- mu.upper.diff + delta   # xbar
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper.diff <- mu.upper.diff + mu0.iterate
  }
}


rand.CI <- tibble(Condition = c("close", "far", "diff"), 
                  lower.limit = c(mu.lower.close, mu.lower.far, mu.lower.diff), 
                  upper.limit = c(mu.upper.close, mu.upper.far, mu.upper.diff))

