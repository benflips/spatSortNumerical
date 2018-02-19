pA0<-0.1 # initial frequency of A

W<-c(0.9, 0.9, 1) #Temporal fitness, AA, Aa, aa
V<-c(0.4,  0.4, 0.1)  #Spatial fitness, AA, Aa, aa

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

example.p<-recDip(pA0, W, V)

plot(example.p, type = "l", xlab="Time", ylab="Allele frequency")
