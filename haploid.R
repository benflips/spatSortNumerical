pA0<-0.1 # initial frequency of A

W<-c(0.9, 1) #Temporal fitness, A, a
V<-c(0.4, 0.1)  #Spatial fitness, A, a

recHap<-function(p.init, Wvec, Vvec, time=30){
  p<-vector(mode = "numeric", length = time)
  p[1]<-p.init
  for (ii in 2:time){
    p[ii]<-Vvec[1]*Wvec[1]*p[ii-1]/(Vvec[1]*Wvec[1]*p[ii-1]+Vvec[2]*Wvec[2]*(1-p[ii-1]))
  }
  p
}

example.p<-recHap(pA0, W, V)

plot(example.p, type = "l", xlab="Time", ylab="Allele frequency")