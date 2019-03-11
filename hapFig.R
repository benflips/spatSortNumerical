##  To plot the basic haploid model under a set of parameters causing fitness to decrease on the invasion front.

# the recursion equation tracking frequency of A over time
#   Wvec and Vvec are length=2 vectors containing fitness values of A and a haplotypes respectively
recHap<-function(p.init, Wvec, Vvec, time=30){
  p<-vector(mode = "numeric", length = time)
  p[1]<-p.init
  for (ii in 2:time){
    p[ii]<-Vvec[1]*Wvec[1]*p[ii-1]/(Vvec[1]*Wvec[1]*p[ii-1]+Vvec[2]*Wvec[2]*(1-p[ii-1]))
  }
  p
}

# generates mean temporal, spatial, and spatiotemporal fitness
meanFitness<-function(pVec, Wvec, Vvec){
  temporal<-pVec*Wvec[1]+(1-pVec)*Wvec[2]
  spatial<-pVec*Vvec[1]+(1-pVec)*Vvec[2]
  spatTemp<-pVec*Wvec[1]*Vvec[1]+(1-pVec)*Wvec[2]*Vvec[2]
  list(t=temporal, s=spatial, st=spatTemp)
}


pA0<-0.1 # initial frequency of A

W<-c(1, 1.4) #Temporal fitness, A, a
V<-c(1, 0.5)  #Spatial fitness, A, a
gens<-35

pFig<-recHap(pA0, W, V, time=gens)
fitFig<-meanFitness(pFig, W, V)
tFig<-1:gens


par(mfrow=c(1, 2))
plot(pFig~tFig, type = "l", xlab="Time (generations)", ylab="Frequency of allele A", bty="l")
plot(fitFig$t~tFig, type = "l", xlab="Time (generations)", ylab="Mean fitness", ylim=range(unlist(fitFig)), bty="l")
lines(fitFig$s~tFig, lty=2)
lines(fitFig$st~tFig, lty=3)
legend('topright', legend = c("Temporal", "Spatial", "Spatiotemporal"), lty = 1:3, bty = "n", title = "Aspect of fitness")
