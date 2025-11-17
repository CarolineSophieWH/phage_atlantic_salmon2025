## Deriving differences between treatments 

setwd("/path/")
## Read metadata information for the samples
metadata = read.csv("metadata.csv", header=T, as.is=T)
## Read OTU and viral OTU tables
OTUs = read.csv("Bozzi_OTUs.csv", sep=";", header=T, as.is=T)
votu = read.csv("vOTU_table.csv", header=T, as.is=T)
## adjust rownames and clean up a little, incl. transposing
row.names(votu) = votu$vOTU
votu$vOTU = NULL
sample = colnames(votu)
votu=as.data.frame(t(votu))
votu$Sample_name = sample

## Merge the data to make a usable dataset 
d = merge(x = metadata, y = OTUs, by.x="Sample_name", by.y="SampleID")
keep_cols = c("Sample_name", "Health_phenotype", "Bacterial_profile", "OTU1", "OTU2")
d = d[, keep_cols]

## Merge votu data into the merged metadata and OTU data frame
d = merge(d, votu, by="Sample_name")

## There are two models to consider here - 
## first, we need to look at the graph where we have, 
## Model 1. HealthStatus -> vOTUs <- OTUs and HealthStatus -> OTUs
## HealthStatus -> vOTUS, and HealthStatus -> OTUs
## Two things to also consider is whether to dichotomize OTUs and vOTUs

## Load the libraries needed
#BiocManager::install("rethinking")
library("rethinking")
library("tidyverse")
library("dagitty")
library("svglite")

dag1 <- dagitty("dag{ Phage <- Health -> BacterialProfile -> Phage }")
impliedConditionalIndependencies(dag1)
plot(dag1)

dag2 <- dagitty("dag{ viralOTUs <- Health -> BacterialProfile }")
impliedConditionalIndependencies(dag2)

# First, we are going to dichotomize the variables
# for vOTUs as having viral load or not. 
d$hasVOTU = as.integer(rowSums(d[,-c(1:5)]) != 0)

## Now let us look at the correlation, and conditional correlation. 
cor(as.integer(d$Bacterial_profile == "Aliivibrio sp."), d$hasVOTU )
# [1] 0.5653769

## conditional on health status
healthy = d$Health_phenotype == "Healthy"
cor(as.integer(d$Bacterial_profile[healthy] == "Aliivibrio sp."), d$hasVOTU[healthy])
cor(as.integer(d$Bacterial_profile[!healthy] == "Aliivibrio sp."), d$hasVOTU[!healthy])
# [1] 0.4844814
# [1] 0.5070926

## Make helper variables we will use in our analyses
d$BP = as.integer(d$Bacterial_profile == "Aliivibrio sp.")
d$hf = as.integer(healthy)

## now let us model hasVOTU as a binomial
## with only BP as being explanatory
m1.hasVOTU <- quap ( alist(
  hasVOTU ~ dbinom(1, p),
  logit(p) <- a + b * BP, 
  a ~ dnorm (0, 1.5),
  b ~ dnorm (0, 0.5)
), data = d)
  
## let us look at the priors to make sure that we are not 
## introducing any differences in the effects through the prior
prior <- extract.prior( m1.hasVOTU , n=1e4)
p <- cbind(inv_logit( prior$a + prior$b ), inv_logit(prior$a))
mean( ( p[,2] - p[,1] ) )
## remember that this is based on a sample so you might see some other values
## [1] 0.0004032817

## Let us now look at the posterior differences
post <- extract.samples(m1.hasVOTU)
p <- cbind(inv_logit(post$a + post$b), inv_logit(post$a))
mean( ( p[,2] - p[,1] ) )
## Again, from sampling
## [1] -0.156087
plot(precis(m1.hasVOTU))

## now include health status in the model 
m2.hasVOTU <- quap ( alist(
  hasVOTU ~ dbinom(1, p),
  logit(p) <- a + b * BP + c * hf, 
  a ~ dnorm (0, 1.5),
  b ~ dnorm (0, 0.5),
  c ~ dnorm (0, 0.5)
), data = d)

## Now let us again check the priors and posteriors, 
## we will use the inv_logit function to get the values on the original 
## probability scale
prior <- extract.prior( m2.hasVOTU , n=1e4)
p_prior <- cbind(inv_logit( prior$a + prior$b + prior$c),
           inv_logit( prior$a + prior$b ),
           inv_logit( prior$a + prior$c ), 
           inv_logit(prior$a))

post <- extract.samples(m2.hasVOTU, n=1e4)
p_post <- cbind(inv_logit( post$a + post$b + post$c),
           inv_logit( post$a + post$b ),
           inv_logit( post$a + post$c ), 
           inv_logit(post$a))


svglite("fig_S5.svg", width = 12.5, height = 8)

## Plot the posterior densities of having a viral otu based on health status 
## and bacterial profile. 
plot(density(p_post[,1]), xlim=c(0,1), ylim=c(0,5), 
     xlab="Viral OTUs", col="goldenrod1", main=" ",
     lwd=3, las=1)
lines(density(p_post[,2]), col="goldenrod1", lwd=3, lty=2)
lines(density(p_post[,3]), col="salmon2", lwd=3)
lines(density(p_post[,4]), col="salmon2", lwd=3, lty=2)
legend(
  "topright",
  legend = c(expression(italic(Aliivibrio)~sp.), expression(italic(Mycoplasma)~sp.)),
  col = c("goldenrod1", "salmon2"),
  lwd = 2,
  bty = "n"
)

dev.off()

legend(x=0.8, y=4, legend = c("Healthy", "Sick"), 
       title="Phenotype", title.font = 2, 
       col=c("black"), lwd=3, lty=c(1,2), bty="n")

fig_S5


plot(precis(m2.hasVOTU))

## Now let us try the other model dag1, and 
## see if the conditional independence holds - actually there is no conditional 
## independence implied by dag1

m3.hasVOTU <- quap ( alist(
  ##  BP <- Health -> hasVOTU <- BP
  hasVOTU ~ dbinom(1, p),
  logit(p) <- a + b * BP + c*hf, 
  a ~ dnorm (0, 1.5),
  b ~ dnorm (0, 0.5),
  c ~ dnorm(0, 0.5),
  ## now BP
  BP ~ dbinom(1, q), 
  logit(q) <- d + e*hf, 
  d ~ dnorm (0, 1.5),
  e ~ dnorm (0, 0.5)
), data = d)

prior <- extract.prior( m3.hasVOTU , n=1e4)
p <- cbind(inv_logit( prior$a + prior$b + prior$c),
           inv_logit( prior$a + prior$b ),
           inv_logit( prior$a + prior$c ), 
           inv_logit(prior$a))
mean( abs( p[,1] - p[,3] ) )

post <- extract.samples(m3.hasVOTU, n=1e4)
p <- cbind(inv_logit( post$a + post$b + post$c),
           inv_logit( post$a + post$b ),
           inv_logit( post$a + post$c ), 
           inv_logit(post$a))
mean( ( p[,1] - p[,3] ) )

## Now let us try dag1
m4.hasVOTU <- quap ( alist(
  ##  BP <- Health -> hasVOTU
  ## first hasVOTU
  hasVOTU ~ dbinom(1, p),
  logit(p) <- a + b * hf,
  a ~ dnorm (0, 1.5),
  b ~ dnorm (0, 0.5),
  ## now BP
  BP ~ dbinom(1, q), 
  logit(q) <- d + e*hf, 
  d ~ dnorm (0, 1.5),
  e ~ dnorm (0, 0.5)
), data = d)


prior <- extract.prior( m4.hasVOTU , n=1e4)
p <- cbind(inv_logit( prior$a + prior$b),
           inv_logit( prior$a),
           inv_logit( prior$d + prior$e ), 
           inv_logit(prior$d))
mean( abs( p[,1] - p[,2] ) )

post <- extract.samples(m4.hasVOTU, n=1e4)
p <- cbind(inv_logit( post$a + post$b),
           inv_logit( post$a ),
           inv_logit( post$d + post$e ), 
           inv_logit(post$d))
mean( abs( p[,1] - p[,2] ) )
precis(m4.hasVOTU)

## let us move on with M2, and model the viral OTU1
## Basically use a poisson model to look at the number of viral otus
## but only when it has any OTUs, we ignore the ones which have none. 
d$numVOTUs = rowSums(d[,c(6:27)] != 0)
d$numVOTUs[d$hasVOTU == 0] = NA

##Fit using a poisson outcome
m1.numVOTU <- quap ( alist(
  numVOTUs ~ dpois(lam),
  log(lam) <- a[BP], 
  a[BP] ~ dnorm (1, 0.75)
), data = list(numVOTUs = d$numVOTUs[!is.na(d$numVOTUs)], 
               BP = (d$BP[!is.na(d$numVOTUs)]+1))
)

## Get priors and posteriors
prior <- extract.prior( m1.numVOTU , n = 1e5 )
counts <- cbind(exp(prior$a[,1]), exp(prior$a[,2]))
mean((counts[,1] - counts[,2]))
## [1] -0.001693049

post <- extract.samples(m1.numVOTU, n=1e5)
counts <- cbind(exp(post$a[,1]), exp(post$a[,2]))
mean((counts[,1] - counts[,2]))
## [1] -8.5312

svglite("fig_S6.svg", width = 12.5, height = 8)

## plot the posterior density 
plot(density(exp(post$a[,2])), xlim=c(0,22), ylim=c(0,0.5), col="goldenrod1", 
     lwd=2, xlab="Viral OTUs", main=" ", las=1)
lines(density(exp(post$a[,1])), col="salmon2", lwd=2)
legend(
  "topright",
  legend = c(expression(italic(Aliivibrio)~sp.), expression(italic(Mycoplasma)~sp.)),
  col = c("goldenrod1", "salmon2"),
  lwd = 2,
  bty = "n"
)

dev.off()
