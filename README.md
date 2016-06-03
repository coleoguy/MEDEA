MEDEA
====

R package for the analysis of MEDEA systems


This package is in the early stages of development and should not be used for any analysis at this point.

```
install.packages("devtools")
library(devtools)
install_github('coleoguy/MEDEA')
library(Medea)


install.packages("devtools")
library(devtools)

install_github('coleoguy/MEDEA', build_vignettes = T)
library(Medea)

# here is an example of running it:

result <- modelMedea(max = 1000,
                     N = 1000000,
                     sims = 100,
                     loss = 0,
                     gain = 0,
                     bottleneck = 0,
                     bottleneck.freq = 7,
                     mon.both = F,
                     fate = T,
                     pop.init = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 950),
                     s = 0,
                     mutation = F,
                     verbose = F)



```

if you have questions or problems please let me know [coleoguy@gmail.com](mailto:coleoguy@gmail.com).
