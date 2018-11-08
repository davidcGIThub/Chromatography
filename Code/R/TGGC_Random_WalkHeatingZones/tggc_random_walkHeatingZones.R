rm(list=ls())
########################

###   NEW EQUATION  

####################

j <- 2 # number of compounds
n.mol <- 100
cols <- c('red','green','magenta','black','blue','grey10')
sig.sq <- 4.073e-6 #3.792e-4  # was 1e-6
taylor <- 1.647e-1 
#C <- c(-12.5075, -12.7202) 
#h <- c(11416, 11521)
X <- c(10,11) #alkane length
C <- -7.22173 - 0.47406*X
h <- 2284 + 819*X

delta.t <- 0.0001 # time step: number of n's per second
R <- 1.987 # Boltzman's constant
col.length <- 200 # column length [cm]

T.0 <- 300 # initial temperature [K]  ******Samuel T.0 should be about 310***********
#T.0 <- 270
u.m <- 55 # velocity of mobile phase [cm/sec]
w <- 6 # velocity of gradient [cm/sec]
tp.rate <- 0.5 #degrees C per second
#s <- 0.1
#c <- 0.0005


col.length <- 200

ramp.length <- 200 # cm
ramp.height <- 80 # degree C

#Temp.vec <- c( ramp.height * seq(from=1, to=0, length = 1e4*ramp.length),
 #              rep(0,times=1e4*(col.length-ramp.length) ) )
 
zone.len <- 1 #cm
window <- 2000
factor <- 1.08^(1/40)
kern <- factor^(-abs(-window:window))
kern <- kern/sum(kern)
holder <- seq(from=ramp.height, to=0, length=ramp.length/zone.len)
holder <- c(holder, rep(0,times=(col.length-ramp.length)/zone.len)-1)
Temp.vec <- rep(holder, each=1e4*zone.len)
for(m in (1+window):(ramp.length*1e4)){
     Temp.vec[m] <- sum(Temp.vec[(m-window):(m+window)]*kern)
}

#Temp.vec <- ramp.height * exp(-(seq(from=0,to=5, length=col.length*1e4)))
#plot(Temp.vec[(1:col.length)*1e4])
     
# pdf('26JuneFig15.pdf', width=7, height=4)
# plot(Temp.vec[seq(1,length(Temp.vec),length=200)], type='l', xlab='Column Length [cm]', ylab='Temperature') 
# dev.off() 
#expo.decay <- exp(-(1:(col.length*10))/10*decay.const/20*3)*T.ramp
#plot(expo.decay[(1:col.length)*10],type='l', xlab='Column Length [cm]', ylab='Temperature')

run.time <- 15000000

##### Start Standard Execution

x <- detector <- matrix(0, nrow=j, ncol=n.mol)
for(i in 1:j){   x[i,] <- sample(n.mol, size=n.mol, replace=F)/10 }
res.end <- matrix(0, nrow=j, ncol=2)
eluted.1 <- eluted.2 <- 1
 detected.1 <- detected.2 <- 0

 # pdf(file='/home/samuel/R/TGGC/lastrun%02d.pdf', height=10, width=12, onefile=F) 
    
for(n in 0:run.time){

T.all <- Temp.vec[pmin(x*1e4+1,col.length*1e4)]+tp.rate*n*delta.t+ T.0
   
 C.all <- C 
  k.all <- exp(h/(R*T.all)+C.all)
   D.t <- ( T.all^1.75*sig.sq + ((1+6*k.all+11*k.all^2)/((1+k.all)^2))*u.m^2*taylor/(T.all^1.75) )*delta.t
 
 
  x <- x + ( u.m*delta.t + rnorm(j*n.mol, mean = 0, sd = sqrt(2*D.t) ) ) * 
    as.numeric( runif(j*n.mol) < 1/(1+k.all ) ) 

#  1/(1+exp( h/(R*( T.all )) + C.all )
 
 if(n%%500==0){
      if(max(x[1,])>col.length){
          all.detected.1 <- which(x[1,] > col.length)
          new.detected.1 <- all.detected.1[!all.detected.1 %in% detected.1]
          detected.1 <- c(detected.1, new.detected.1)
          detector[1, new.detected.1] <- n*delta.t
          
     }
     if(max(x[2,])>col.length){
          all.detected.2 <- which(x[2,] > col.length)
          new.detected.2 <- all.detected.2[!all.detected.2 %in% detected.2]
          detected.2 <- c(detected.2, new.detected.2)
          detector[2, new.detected.2] <- n*delta.t
          if(length(detected.2)+length(detected.1)-2 >= 2*n.mol){ break() }
     }     
 }
 

  if(n%%50000==0){
    if(min(x)>col.length) break()
    par(mfrow=c(2,1), mar=c(2,1,1,6), las=1)
    plot(0, 0, type='n', xlim=c(0,col.length), ylim=c(0, 1), yaxt='n',ylab='')#
    mtext(text=paste(' Time:\n ', n*delta.t,' s'), side=4, cex=2)
    for(i in 1:j){
      points(x[i,], (1:n.mol)/n.mol, col=cols[i], pch=i, cex=0.5) 
      text(mean(x[i,]), 0.2+i*0.2, labels=round(4*sd(x[i,]),2 ), col=cols[i],pos=2,cex=2)
    }
    
    lines(xx <-seq(1, col.length, 0.1),  (Temp.vec[xx*1e4]+tp.rate*delta.t*n)/200,col='red')
    #abline(h=( (delta.t*n*b) )/200,col='red') #programmed temperature
     
    plot(density(x[1,]), type='n', xlim=c(0,col.length), yaxt='n', xlab='',ylab='',main='')
    for(i in 1:j){ lines(density(x[i,]), col=cols[i] ) }
    
    top.pks <- pk.head <- pk.tail <- numeric(j)
    peakx <- numeric(1)
    for(i in 1:j){
      temp <- density(x[i,])
      data <- temp$y
      peakx <- which.max(data)
      top.pks[i] <- temp$x[peakx]
      pk.head[i] <- temp$x[ max(which(data[1:peakx]<data[peakx]/2)) ]
      pk.tail[i] <- temp$x[ min(which(data[peakx:length(data)]< data[peakx]/2))+peakx ]
    }
    
    res <- numeric(j-1)
    for(i in 2:j){
      res[i-1] <- abs( (top.pks[i]-top.pks[i-1]) / 
                            ( ( 1.7*(pk.head[i-1]-pk.tail[i-1]) + 1.7*(pk.head[i]-pk.tail[i]) )/2) ) 
    }
    mtext(text=paste(' Ave \n Res:\n ', round(mean(res),2)), side=4, cex=2)
  } 
}
par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)

#### End Standard Execution

dev.off()

retention.time <- rowMeans(detector)
stand.dev <- apply(detector,1,sd)
resolutions <- rep(0,j)
      for(i in 2:j){resolutions[i] <- (retention.time[i]-retention.time[i-1])/(2*(stand.dev[i]+stand.dev[i-1]))} 
(rbind(stand.dev,retention.time, resolutions))


### Continuous ramp

#30 cm, 80 degree
#stand.dev        0.3566578   0.4373644
#retention.time 220.6815000 217.8245000
#resolutions      0.0000000  -1.7990680

# 30 cm, 60 degree
#stand.dev        0.3872775   0.3983207
#retention.time 221.1960000 218.2985000
#resolutions      0.0000000  -1.8441363

# 60 cm, 20 degree
#stand.dev        0.5101126   0.5066913
#retention.time 220.3285000 217.5375000
#resolutions      0.0000000  -1.3724376


### Exponetial decay, 200cm, 5 exp.fac, 270 T.0, continuous
# stand.dev        0.4279255   0.4471038
#retention.time 204.1810000 201.3835000
#resolutions      0.0000000  -1.5985181


## real alkane values: C10,C11
#####################

## 30 cm, 20 degree
#stand.dev        0.563449   0.4527344
#retention.time 130.640000 161.0175000
#resolutions      0.000000  14.9468594

# 40 degree ramp
# stand.dev        0.4287885   0.4234407
# retention.time 131.3270000 161.7290000
# resolutions      0.0000000  17.8367507

# 10 degree ramp
# stand.dev        0.6378901   0.6330021
#retention.time 131.4905000 161.9905000
#resolutions      0.0000000  11.9994440

# 30 degree ramp on 15 cm
# stand.dev        0.4899132   0.5231017
#retention.time 131.8605000 162.3990000
#resolutions      0.0000000  15.0730766


