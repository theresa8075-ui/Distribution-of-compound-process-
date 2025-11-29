# Distribution of Compound Poisson process
This repository consists of an R Shiny project that offers an opportunity to visualize and study an interactive compound Poisson process.The inter-arrival time of the events follows the exponential distribution and jump sizes associated with those events are also exponentially distributed. Users can modify the most important parameters of the model, such as the inter-arrival rate( λ) and the jump-size rate(β), and immediately see how these changes affect both the distribution and growth of the aggregate process S(t) over time. Using dynamic histograms and real-time computations, the Shiny app conveys an understanding of the sensitivity of S(t) with respect to different parameter settings. This becomes a powerful tool for learning, simulation, and stochastic process exploration.

# Objectives of the project
* Derive the distribution of compound poissoon process
* Plot the histogram at t=10,100,1000,10000
* Creating an R shiny app to visualize the histograms.
* Understanding the sensitivity of the parameters(inter-arrival rate( λ) and the jump-size rate(β)) on S(t) vs time

# Impact of inter arrival parameters on s(t)
*λ controls how often jumps happen.

When λ is large, events happen more frequently, so S(t) becomes bigger and smoother.

When λ is small, events are rare, so S(t) often stays near zero and the histogram is very skewed.

*β controls how big each jump is.

When β is large, each jump is small, so S(t) grows slowly.

When β is small, each jump is large, so S(t) becomes bigger and much more spread out. 

*Together, they decide the overall size and shape of S(t).

Thus both the mean and variability of S(t) increase as λ increases or β decreases since

E[S(t)]=λt\β​,Var(S(t))=2λt\β^2​.

# How to run the application 
It can be done by installing the required package called "shiny".

install.packages("shiny")


