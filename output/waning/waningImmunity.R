library(distributions3)
library(ggplot2)

#Erlang(k, lambda)
#When k = 1, we have the exponential function
dur <- 180 #durability of mAb
erL_1 <- Erlang(1, 1/dur)
erL_2 <- Erlang(2, 2/dur)
erL_3 <- Erlang(3, 3/dur)
erL_4 <- Erlang(4, 4/dur)
erL_7 <- Erlang(7, 7/dur)
erL_10 <- Erlang(10, 10/dur)
#erL_15 <- Erlang(15, 15/dur)

cdf_1 <- 1- cdf(erL_1, 0:300)
cdf_2 <- 1- cdf(erL_2, 0:300)
cdf_3 <- 1 - cdf(erL_3, 0:300)
cdf_4 <- 1 - cdf(erL_4, 0:300)
cdf_7 <- 1 - cdf(erL_7, 0:300)
cdf_10 <- 1 - cdf(erL_10, 0:300)
#cdf_15 <- 1 - cdf(erL_15, 0:300)

er_data <- cbind(cdf_1, cdf_2, cdf_3, cdf_4, cdf_7, cdf_10)
#er_data <- cbind(cdf_1, cdf_2, cdf_3, cdf_4, cdf_7, cdf_10, cdf_15)

melt_er <- reshape2::melt(as.matrix(er_data))
melt_er$Var1 <- rep(0:300, 6)
#melt_er$Var1 <- rep(0:300, 7)
colnames(melt_er) <- c("Days", "Distribution", "Proportion")

#Find average value of 1-CDF for different values of k
#of Erlang distribution over set number of days
av_erl <- rep(NA, 6)
#av_erl <- rep(NA, 7)
days <- 150

k = 1
rate = k/dur
integrand <- function(x) {exp(-rate*x)}
result <- integrate(integrand, lower = 0, upper = days)
av_erl[1] <- result$value/days

k = 2
rate = k/dur
integrand <- function(x) {exp(-rate*x) + rate*x*exp(-rate*x)}
result <- integrate(integrand, lower = 0, upper = days)
av_erl[2] <- result$value/days

k = 3
rate = k/dur
integrand <- function(x) {exp(-rate*x) + rate*x*exp(-rate*x) + (1/2)*rate^2*x^2*exp(-rate*x)}
result <- integrate(integrand, lower = 0, upper = days)
av_erl[3] <- result$value/days

k = 4
rate = k/dur
integrand <- function(x) {exp(-rate*x) + rate*x*exp(-rate*x) +
                            (1/2)*rate^2*x^2*exp(-rate*x) +
                             (1/6)*rate^3*x^3*exp(-rate*x)}
result <- integrate(integrand, lower = 0, upper = days)
av_erl[4] <- result$value/days

k = 7
rate = k/dur
integrand <- function(x) {exp(-rate*x) + rate*x*exp(-rate*x) +
    (1/2)*rate^2*x^2*exp(-rate*x) +
    (1/6)*rate^3*x^3*exp(-rate*x) +
    (1/factorial(4))*rate^4*x^4*exp(-rate*x) +
    (1/factorial(5))*rate^5*x^5*exp(-rate*x) +
    (1/factorial(6))*rate^6*x^6*exp(-rate*x)}
result <- integrate(integrand, lower = 0, upper = days)
av_erl[5] <- result$value/days

k = 10
rate = k/dur
integrand <- function(x) {exp(-rate*x) + rate*x*exp(-rate*x) +
    (1/2)*rate^2*x^2*exp(-rate*x) +
    (1/6)*rate^3*x^3*exp(-rate*x) +
    (1/factorial(4))*rate^4*x^4*exp(-rate*x) +
    (1/factorial(5))*rate^5*x^5*exp(-rate*x) +
    (1/factorial(6))*rate^6*x^6*exp(-rate*x) +
    (1/factorial(7))*rate^7*x^7*exp(-rate*x) +
    (1/factorial(8))*rate^8*x^8*exp(-rate*x) +
    (1/factorial(9))*rate^9*x^9*exp(-rate*x)}
result <- integrate(integrand, lower = 0, upper = days)
av_erl[6] <- result$value/days

# k = 15
# rate = k/dur
# integrand <- function(x) {exp(-rate*x) + rate*x*exp(-rate*x) +
#     (1/2)*rate^2*x^2*exp(-rate*x) +
#     (1/6)*rate^3*x^3*exp(-rate*x) +
#     (1/factorial(4))*rate^4*x^4*exp(-rate*x) +
#     (1/factorial(5))*rate^5*x^5*exp(-rate*x) +
#     (1/factorial(6))*rate^6*x^6*exp(-rate*x) +
#     (1/factorial(7))*rate^7*x^7*exp(-rate*x) +
#     (1/factorial(8))*rate^8*x^8*exp(-rate*x) +
#     (1/factorial(9))*rate^9*x^9*exp(-rate*x) +
#     (1/factorial(10))*rate^10*x^10*exp(-rate*x) +
#     (1/factorial(11))*rate^11*x^11*exp(-rate*x) +
#     (1/factorial(12))*rate^12*x^12*exp(-rate*x) +
#     (1/factorial(13))*rate^13*x^13*exp(-rate*x) +
#     (1/factorial(14))*rate^14*x^14*exp(-rate*x)}
# result <- integrate(integrand, lower = 0, upper = days)
# av_erl[7] <- result$value/days


#Plot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

ggplot2::ggplot() +
  ggplot2::geom_line(data = melt_er, ggplot2::aes(x = Days, y = Proportion , col = Distribution), size = 1.1) +
  ggplot2::theme_bw()+
  ggplot2::scale_y_continuous(name = "Proportion of vaccinees still protected", breaks = seq(0, 1, 0.1))+
  ggplot2::scale_x_continuous(name = "Days post vaccination", breaks = seq(0, 300, 50)) +
  ggplot2::expand_limits(x = 0, y = 0) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::scale_color_discrete(name = "Erlang",
                       labels = c("k = 1", "k = 2", "k = 3", "k = 4", "k = 7", "k = 10")) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl[1], xend = days, yend = av_erl[1]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[1]) +
  ggplot2::geom_text(aes(x = days + 10, y= av_erl[1], label= round(av_erl[1], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl[2], xend = days, yend = av_erl[2]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[2]) +
  ggplot2::geom_text(aes(x = days + 10, y= av_erl[2], label= round(av_erl[2], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl[3], xend = days, yend = av_erl[3]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[3]) +
  ggplot2::geom_text(aes(x = days + 10, y= av_erl[3], label= round(av_erl[3], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl[4], xend = days, yend = av_erl[4]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[4]) +
  ggplot2::geom_text(aes(x = days + 10, y= av_erl[4], label= round(av_erl[4], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl[5], xend = days, yend = av_erl[5]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[5]) +
  ggplot2::geom_text(aes(x = days + 10, y= av_erl[5], label= round(av_erl[5], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl[6], xend = days, yend = av_erl[6]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[6]) +
  ggplot2::geom_text(aes(x = days + 10, y= av_erl[6], label= round(av_erl[6], digits = 2)), size = 4)


#Keep k = 3, change durability
dur <- c(100, 150, 180, 200, 210, 220) #durability of mAb
erL3_1 <- Erlang(3, 3/dur[1])
erL3_2 <- Erlang(3, 3/dur[2])
erL3_3 <- Erlang(3, 3/dur[3])
erL3_4 <- Erlang(3, 3/dur[4])
erL3_5 <- Erlang(3, 3/dur[5])
erL3_6 <- Erlang(3, 3/dur[6])

cdf3_1 <- 1- cdf(erL3_1, 0:300)
cdf3_2 <- 1- cdf(erL3_2, 0:300)
cdf3_3 <- 1 - cdf(erL3_3, 0:300)
cdf3_4 <- 1 - cdf(erL3_4, 0:300)
cdf3_5 <- 1 - cdf(erL3_5, 0:300)
cdf3_6 <- 1 - cdf(erL3_6, 0:300)

er3_data <- cbind(cdf3_1, cdf3_2, cdf3_3, cdf3_4, cdf3_5, cdf3_6)

melt_er3 <- reshape2::melt(as.matrix(er3_data))
melt_er3$Var1 <- rep(0:300, 6)
colnames(melt_er3) <- c("Days", "Distribution", "Proportion")

#Find average value of 1-CDF of Erlang distribution over set number of days
av_erl <- rep(NA, 6)
days <- 150
k = 3

for(i in 1:6)
{
  rate = k/dur[i]
  integrand <- function(x) {exp(-rate*x) + rate*x*exp(-rate*x) + (1/2)*rate^2*x^2*exp(-rate*x)}
  result <- integrate(integrand, lower = 0, upper = days)
  av_erl[i] <- result$value/days
}
av_erl_150 <- av_erl

av_erl <- rep(NA, 6)
days <- 180

for(i in 1:6)
{
  rate = k/dur[i]
  integrand <- function(x) {exp(-rate*x) + rate*x*exp(-rate*x) + (1/2)*rate^2*x^2*exp(-rate*x)}
  result <- integrate(integrand, lower = 0, upper = days)
  av_erl[i] <- result$value/days
}
av_erl_180 <- av_erl

#Plot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

ggplot2::ggplot() +
  ggplot2::geom_line(data = melt_er3, ggplot2::aes(x = Days, y = Proportion , col = Distribution), size = 1.1) +
  ggplot2::theme_bw()+
  ggplot2::scale_y_continuous(name = "Proportion of vaccinees still protected", breaks = seq(0, 1, 0.1))+
  ggplot2::scale_x_continuous(name = "Days post vaccination", breaks = seq(0, 300, 50)) +
  ggplot2::expand_limits(x = 0, y = 0) +
  ggplot2::coord_cartesian(expand = FALSE) +
  ggplot2::scale_color_discrete(name = "Average durability",
                               # labels = c("100 days", "150 days", "180 days", "200 days", "250 days", "300 days")) +
                               labels = c("100 days", "150 days", "180 days", "200 days", "210 days", "220 days")) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_150[1], xend = 150, yend = av_erl_150[1]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[1]) +
  ggplot2::geom_text(aes(x = 150 + 10, y= av_erl_150[1], label= round(av_erl_150[1], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_150[2], xend = 150, yend = av_erl_150[2]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[2]) +
  ggplot2::geom_text(aes(x = 150 + 10, y= av_erl_150[2], label= round(av_erl_150[2], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_150[3], xend = 150, yend = av_erl_150[3]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[3]) +
  ggplot2::geom_text(aes(x = 150 + 10, y= av_erl_150[3], label= round(av_erl_150[3], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_150[4], xend = 150, yend = av_erl_150[4]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[4]) +
  ggplot2::geom_text(aes(x = 150 + 10, y= av_erl_150[4], label= round(av_erl_150[4], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_150[5], xend = 150, yend = av_erl_150[5]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[5]) +
  ggplot2::geom_text(aes(x = 150 + 10, y= av_erl_150[5], label= round(av_erl_150[5], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_150[6], xend = 150, yend = av_erl_150[6]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[6]) +
  ggplot2::geom_text(aes(x = 150 + 10, y= av_erl_150[6], label= round(av_erl_150[6], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_180[1], xend = 180, yend = av_erl_180[1]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[1]) +
  ggplot2::geom_text(aes(x = 180 + 10, y= av_erl_180[1], label= round(av_erl_180[1], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_180[2], xend = 180, yend = av_erl_180[2]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[2]) +
  ggplot2::geom_text(aes(x = 180 + 10, y= av_erl_180[2], label= round(av_erl_180[2], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_180[3], xend = 180, yend = av_erl_180[3]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[3]) +
  ggplot2::geom_text(aes(x = 180 + 10, y= av_erl_180[3], label= round(av_erl_180[3], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_180[4], xend = 180, yend = av_erl_180[4]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[4]) +
  ggplot2::geom_text(aes(x = 180 + 10, y= av_erl_180[4], label= round(av_erl_180[4], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_180[5], xend = 180, yend = av_erl_180[5]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[5]) +
  ggplot2::geom_text(aes(x = 180 + 10, y= av_erl_180[5], label= round(av_erl_180[5], digits = 2)), size = 4) +
  ggplot2::geom_segment(aes(x = 0, y = av_erl_180[6], xend = 180, yend = av_erl_180[6]),
                        linetype = "dashed", size = 1.1, col = ggplotColours()[6]) +
  ggplot2::geom_text(aes(x = 180 + 10, y= av_erl_180[6], label= round(av_erl_180[6], digits = 2)), size = 4)

