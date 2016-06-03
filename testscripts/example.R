
install.packages("devtools")
library(devtools)

install_github('coleoguy/MEDEA', build_vignettes = T)
library(Medea)


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

