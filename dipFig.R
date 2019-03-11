##  To plot the basic diploid model under a set of parameters causing fitness to decrease on the invasion front.

# the recursion equation tracking frequency of A over time
#   Wvec and Vvec are length=3 vectors containing fitness values of AA, Aa and aa genotypes respectively
recDip<-function(p.init, Wvec, Vvec, time=30){
  p<-vector(mode = "numeric", length = time)
  p[1]<-p.init
  for (ii in 2:time){
    numerator<-Vvec[1]*Wvec[1]*p[ii-1]^2+Vvec[2]*Wvec[2]*p[ii-1]*(1-p[ii-1])
    denominator<-Vvec[1]*Wvec[1]*p[ii-1]^2+2*Vvec[2]*Wvec[2]*p[ii-1]*(1-p[ii-1])+Vvec[3]*Wvec[3]*(1-p[ii-1])^2
    p[ii]<-numerator/denominator
  }
  p
}

# generates mean temporal, spatial, and spatiotemporal fitness
#   Assumes HWE...
meanFitnessDip<-function(pVec, Wvec, Vvec){
  temporal<-pVec^2*Wvec[1]+pVec*(1-pVec)*Wvec[2]+(1-pVec)^2*Wvec[3]
  spatial<-pVec^2*Vvec[1]+pVec*(1-pVec)*Vvec[2]+(1-pVec)^2*Vvec[3]
  spatTemp<-pVec^2*Wvec[1]*Vvec[1]+pVec*(1-pVec)*Wvec[2]*Vvec[2]+(1-pVec)*Wvec[3]*Vvec[3]
  list(t=temporal, s=spatial, st=spatTemp)
}


pA0<-0.1 # initial frequency of A

W<-c(1, 1, 1.1) #Temporal fitness, AA, Aa, aa
V<-c(1, 1, 0.5)  #Spatial fitness, AA, Aa, aa
gens<-60

pFig<-recDip(pA0, W, V, time=gens)
fitFig<-meanFitnessDip(pFig, W, V)
tFig<-1:gens


par(mfrow=c(1, 2))
plot(pFig~tFig, type = "l", xlab="Time (generations)", ylab="Frequency of allele A", bty="l")
plot(fitFig$t~tFig, type = "l", xlab="Time (generations)", ylab="Mean fitness", ylim=range(unlist(fitFig)), bty="l")
lines(fitFig$s~tFig, lty=2)
lines(fitFig$st~tFig, lty=3)
legend('bottomright', legend = c("Temporal", "Spatial", "Spatiotemporal"), lty = 1:3, bty = "n", title = "Aspect of fitness")
