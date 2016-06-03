max = 1000
N = 1000000
sims = 1
loss = 0
gain = 0
bottleneck = 0
bottleneck.freq = 7
mon.both = T
fate = F
pop.init = c(0, 0, 0, 0, 0, 0, 0, 0, 
             0, 0, 0, 0, 0, 0, 50, 950)
s <- 0 # cost of a Madea allele additive fitness model and gene action
#ABAB ABAb ABaB ABab AbAB AbAb AbaB Abab 
#aBAB aBAb aBaB aBab abAB abAb abaB abab

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

