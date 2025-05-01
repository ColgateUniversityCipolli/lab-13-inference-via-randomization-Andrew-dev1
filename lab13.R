## importing libraries
library(tidyverse)
library(pwr)
library(e1071)
library(boot.pval)
library(boot)


zebra.finches.dat <- read_csv("zebrafinches.csv")

further.vec <- zebra.finches.dat$further
skew <- skewness(further.vec)
n <- length(further.vec)

further.test <- t.test(x=further.vec, mu = 0)
t.far <- as.numeric(further.test$statistic) #finding the t value 

# calculating cdf and pdf
pdf <- dnorm(t.far)
error = (skew / sqrt(n)) * ((2 * (t.far^2) + 1) / 6) * pdf


t.vecs <- seq(-10,10,0.1)
results <- tibble(
  t = t.vecs,
  pdf = dnorm(t.vecs), 
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
R <- 1000
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

resamples.shifted <- resample$close.t - change
mean(resamples.shifted)

resamples.shifted <- resample |>
  mutate(close.shifted = close.t - close.change)  |>
  mutate(far.shifted = far.t - far.change) |>
  mutate(diff.shifted = diff.t - diff.change) |>
  select(c(close.shifted, far.shifted, diff.shifted ))

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

(comparisons <- tibble(" "= c("close", "far", "diff"), Bootstrapped = 
                        c(boot.p.closer, boot.p.far, boot.p.diff),
                      t.test = c(close.ttest.p, far.ttest.p, diff.ttest.p)))


## Part c
close.5th <- quantile(resamples.shifted$close.shifted, 0.05)
far.5th <- quantile(resamples.shifted$far.shifted, 0.05)
diff.5th <- quantile(resamples.shifted$diff.shifted, 0.05)

close.percentile.t <- qt(0.05, df = n-1)
far.percentile.t <- qt(0.05, df = n-1)
diff.percentile.t <- qt(0.05, df = n-1)

(comparisons <- tibble(" "= c("close", "far", "diff"), sampled = 
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


