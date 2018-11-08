rm(list=ls())
setwd("~/work/chromatography/Rcode/annotated/ExpSerpentine")
########################

###   NEW EQUATION:  This is a step by step program to generate a random walk model
###     for a varietly of moving thermal gradient profiles.
###		This code was UPDATED in Riva with Samuel in May 2016
####################

j <- 2 # number of compounds
n.mol <- 100 #****Number of molecules per compound used in simulation
cols <- c('red','green','magenta','black','blue','grey10')
sigma.not <- 10000*sqrt(1/12)

gamma.1 <-5.113e-3 # This is the molecular diffusion coefficient when using both temp and pressure.

taylor <- 1.647e-1  #  This is the most recent based on solving Anzi's data, using R=1.987.  (141105)
#                     This is the dispersion coefficient due to resistence to transport and taylor dispersion.
gamma.3 <-7.676e-8 # This is resistence to flow for adsorption/deporption (no pressure variable)

gamma.2 <- 1.217e-9 # This is Golay's formula for resistence to transport when using both temp and pressure in diffusion expression.
chk.sm <- 0
#
#
#C <- c(-12.5075, -12.7202) 
#h <- c(11416, 11521)
X <- c(9,10) #alkane length
C <- -7.22173 - 0.47406*X  #This is the formula for calculating C where X is the number of carbons.  Includes phase ratio
h <- 2284 + 819*X    #  calculation of enthalpy (need to check this)
h <- h* 4.18442979  # This is to change enthalpy from kcal units to Joules in the units.
C <- c(-11.1277,-12.4829)
h <- c(34926.27,40845.71)

C <- c(-11.1277,-11.2829)
h <- c(34926.27,35245.71)

delta.t <- 0.0001 # 0.001 ***# time step: Reciprocol is the number of time steps per second
R <- 8.3144621 #  Boltzman's constant in Joules per Kelvin mole.  I used to use R = 1.987 kcal per Kelvin mole.
col.length <- 20#0 #*** column length [cm]
col.length <- 1000 # measured in cm

T.0 <- 323 # initial temperature [Kelvin, K]
T.0 <- 313
T.0 <- 290
diameter <- 0.0100 #inner diameter of the column in centimeters.
#diameter <- 0.0001 #measured in meters for test.


p.i <- 57.16  #This is the column inlet pressure by gage measured in psi.  For true inlet must add atmosphere pressure.
#p.i <- 30
p.o <- 14.69595   #This is the outlet pressure or Atmospheric pressure in psi.
p.o <- 13.8



p.i <- p.i*6894.757 # This converts pressure from psi units to Pascals  (Anzi solution)
p.o <- p.o*6894.757 # This converts pressure from psi units to Pascals

#p.i <- p.i*10000/1.450377  #converting psi to Pa (pascals)  (Wiki Solution)
#p.o <- p.o*10000/1.450377  #converting psi to Pa

p.i <- p.i + p.o  #This is the true inlet pressure.



w <- 6 # velocity of gradient [cm/sec]
tp.rate <- 0.5 #degrees C per second. This is the rate at which the column is heated after a gradient is established.
#tp.rate <- 1.7
tp.rate <- .25 #This is the more reasonable rate of raising temperature under TPGC

#             rep(0,times=1e4*(col.length-ramp.length) )

Temp.granularity <- 1e2#4 # granularity of temperatures across column. data points/cm

power.ace <- 100


###################
#
#		Here is the temperature profile function.
#		It produces the initial heat gradient profile that we heat.
#		We need to make this a function call to either a data file of a function.
#
#		The values needed are:
#			Temp.granularity
#			ramp.length
#			ramp.height
#			zone.len
#			window
#			col.length
#
#		Return a vector  Temp.vec that has the profile of temperature along
#		the column.  We will then heat this profile up.
#
##################


plate.len <-  10 #60 #cm  Note that the plate width is calculated by subtracting the plate.len
#from the column length and then dividing the difference by the number of
# turns.  This are called isotherms.
#
decay.const <- 0.5 #0.05  # THis is the exponential decay factor, giving the rate of decay
T.ramp <- 50 #height of the temperature ramp; delta T across plate.  This is not T.0
expo.decay <- exp(-(1:(plate.len*Temp.granularity))/Temp.granularity*decay.const)*T.ramp
#This is determines the exponential decay
#
#plot(expo.decay[(1:plate.len)*Temp.granularity],type='l', xlab='Temperature', main=paste("PlateLength: ",plate.len, " DecayConst: ", decay.const, "NumbCorners: ", n.corners))     
# Execute the above line to look at exponential temperature decay down plate.

n.corners <- 200# number of switchbacks, and rungs  #This turns out to be the number of isotherms.
# n.corners <- col.length-plate.len
len.cross <- (col.length-plate.len)/n.corners  #This gives the length of each isotherm.

corner <- round(seq(1,length(expo.decay), length=n.corners+1))  #indexes the position of each of the corners
# on the exponential decay curve.
Temp.vec <- NULL
for(jj in 1:n.corners)
{
  Temp.vec <- c(Temp.vec, expo.decay[corner[jj]:(corner[jj+1])]) #This gives the temperature at each
  # of the corners
  Temp.vec <- c(Temp.vec, rep(expo.decay[corner[jj+1]],times=len.cross*Temp.granularity))
  #  This determines the temperature at each point along the isotherm.
}
Temp.vec <- Temp.vec[1:(col.length*Temp.granularity)]

#plot(Temp.vec[1:(col.length*10)*Temp.granularity/10],type='l', xlab='Temperature', main=expression(paste("PlateLength: ",plate.len, " DecayConst: ", decay.const, "NumbCorners: ", n.corners))) # Execute to look at temperature along column


window <- 100   #This is for the ksmooth function, below.  It should be proportional to the
# temperature granularity above.  Change either and you must change the other.
print(Temp.vec)
Temp.vec <- ksmooth(1:length(Temp.vec),Temp.vec,'normal',window)$y  #This smooths the abrupt step function produced above.
#  NOTE!!!  ksmooth takes a long time.  Samuel is not sure it is necessary.   I am thinking that if
#  it turns out to be necessary, then we could archive these gradient curves and recall as necessary.
#
Temp.vec <- unlist(Temp.vec)  #Here we make the 'list' structure from the object above into an array.
#plot(Temp.vec[(1:(length(Temp.vec)/1e3))*1e3],type='l')
#******************
#  Here I set the gradient to zero, giving programmed temperature or isothermal
#
#******************

#Temp.vec <-rep(0,length(Temp.vec))
#tp.rate <- 0
###################
#
#		This is the end of the temperature profile function
#
###################

Temp2Pow.vec <- seq(from=0, to=1000,by=1/power.ace)^1.646 #address the temperature /100, get the temp raised to the power
TempNoPow.vec <- seq(from=0, to=1000,by=1/power.ace) #address the temperature /100, get the temp. No Power here (HDT)

#################
#
#	The vector 'TempNoPow.vec' containing the basic temperature profile or template 
#		is returned from temperature profile function call.
#	The vector 'Temp2Pow' is sequence of temperatures raised to the 1.646 power
#		to adjust for the viscosity of helium.
#    
#################


run.time <- 1500000  #Number of potential steps to get compounds out of column.
#run.time <- 6000

pause.count <- 10000  #number of frames to skip in plotting results.  The 'big.matrix'
# takes a snapshot at every 'pause.count' time for illustrating the
# progress of the molecules.

big.matrix <- matrix(0,nrow=run.time/pause.count+1,ncol=n.mol*j)
#  This has the location of every molecule at each recorded time point.  Only
#   time points recorded are every "pause.count" apart.

##### Start Standard Execution

x <- detector <- matrix(0, nrow=j, ncol=n.mol)	# This matrix keeps track of x-axis location of each
# molecule of each analyte. Initially each is at point
# '0' as indicated here.  The y-axis location is simply the
# scaled value of the molecule's index number in the matrix x.

injection.width <- 2  #  Width of injection plug at time zero in cm.											

for(i in 1:j){   x[i,] <- sample(n.mol, size=n.mol, replace=F)*injection.width/n.mol } 	

# The above loop returns an order from 1 to n.mol for each molecule.  Since it is divided
# by n.mol, in essence we get each molecule uniquely with a value between 0 and 1.
# Multiplying by injection.width gives a plug that has the n.mol molecules uniformly distributed 
# from 0 cm to injection.width cm at the beginning of the column.
#
#   Note that the standard deviation below is calculated from s^2 and not from the var
#   of the uniform. HDT
#
#   Remember, x[i,] is the vector of locations in x-axis of i-th analyte.  This location is used 
#	at each step for each molecule to determine temperature and velocity. 
#	This is also used to plot location of the molecules.  The y-axis values for plot are determined
#	from the matrix column index.  (Rows indicate analyte and column indicate individual molecules.)
#

move.count <- 0

res.end <- matrix(0, nrow=j, ncol=2)
eluted.1 <- eluted.2 <- 1
detected.1 <- detected.2 <- 0

##########  Temporary pseudo-do-loop for debugging.

# n <- 5

##########  This is where the do-loop for n starts.

total.time.1 <- 0
total.time.2 <- 0
spread.1 <- 0
spread.2 <- 0

vel1_x = c()
vel2_x = c()
pos1_x = c()
pos2_x = c()
std1 = c()
std2 = c()

for(n in 0:run.time){
  Temp.plus <- Temp.vec + tp.rate*n*delta.t + T.0	
  # Here we heat up the whole column an amount 'tp.rate*delta.t 
  # from the last time increment.  To avoid numerical problems
  # we simply multiply the temperature increment amount by the
  # number of time steps, n, and add to initial temperature.
  # Initial temperature is T.0 plus the temperature gradient profile.
  
  #************
  # Now I hold the temperature constant.  This is for isothermal of gradient
  #  NOTE: Comment this out if allow heating of column over time.
  #Temp.plus <- Temp.vec + T.0
  #
  #*************
  
  T.all <- Temp.plus[pmin(abs(x*Temp.granularity+1),col.length*Temp.granularity)]
  
  C.all <- C 
  k.all <- exp(h/(R*T.all)+C.all)
  
  #   Here I adjust the diffusion calculation according to the position of each molecule and the temperature
  #   at that location, taking adjusting the velocity resulting from the different pressure and viscosity values.  (HDT)
  # 	To do this there are several steps.
  
  #	First step is to determine a numerical value for the integral of temperature times viscosity.
  #
  T.sum.pow <- cumsum(Temp2Pow.vec[Temp.plus*power.ace]) / Temp.granularity 	
  # integral of Temp^1.646
  # Note baseline viscosity cancels out.
  # T.sum.pow[length(T.sum.pow)]
  # is sum or integral from 0 to end.
  
  #***********
  tot.int <-T.sum.pow[col.length*Temp.granularity]  #This is the integral from 0 to L of T^1.646. 
  He.factor <- 0.000018662/(273.15^.646)  #Mult factor viscosity of He solvent.  This is delta.not in notes.
  
  
  C.star <- (p.i^2-p.o^2)/(2*He.factor*T.sum.pow[col.length*Temp.granularity]) 
  
  ## C.star is =to (pi^2-po^2)/(2*integral over L of T^1.646)
  #	I am doing the integral as simply the sum of the histogram pieces.  I could
  #	improve this step by using Simpson's rule or the Newton binomial trick if the 
  #	Temperature can be modeled as a simple function of two variates. (HDT)
  
  
  
  
  velocity.x <- (p.i^2-p.o^2)*Temp.plus[pmin(abs(x*Temp.granularity +1),col.length*Temp.granularity)] * (diameter/2)^2 /
    ((16*He.factor*tot.int )*sqrt(abs(p.i^2-(p.i^2-p.o^2)*T.sum.pow[pmin(abs(x*Temp.granularity+1),col.length*Temp.granularity)]/tot.int )))
  
  #  Velocity updated by HDT on 31 March 2015
  #
  #  
  
  
  p.x <- sqrt(abs((p.i^2-(p.i^2-p.o^2)*T.sum.pow[pmin(abs(x*Temp.granularity+1),col.length*Temp.granularity)]/tot.int)) ) 
  
  #This is the pressure at x
  #
  # Pressure updated by HDT on 31 March 2015
  
  
  Temp.prss.rat <- (T.all^1.75)/p.x
  
  
  
  D.t <- 2*(( Temp.prss.rat*gamma.1 + ((1+6*k.all+11*k.all^2)/((1+k.all)^2))*velocity.x^2*gamma.2/(Temp.prss.rat)+
                gamma.3*velocity.x^2*k.all/(1+k.all) )/(1+k.all))
  
  sigma.delta <- sqrt(abs(D.t*delta.t))
  #  
  #   
  x <- x + ( velocity.x*delta.t/(1 + k.all) + rnorm(j*n.mol, mean = 0, sd = sigma.delta ) ) 
  
  #
  #
  #	For the Milshtein correction term we have
  #  W.Lang <- rnorm(j*n.mol, mean = 0, sd = sigma.delta )
  # x <- x + ( velocity.x*delta.t +W.Lang +sigma.delta^2*((W.lang-delta.t)^2)/2 )/(1 + k.all)
  #
  #
  mu.1 <- mean(x[1,])
  mu.2 <- mean(x[2,])
  if(mu.1<= col.length){
    total.time.1 <- n*delta.t
    speed.1 <- velocity.x[1]/(1+k.all[1])
    spread.1 <- sd(x[1,])*(1+k.all[1])/(60*velocity.x[1])
  }
  
  if(mu.2<= col.length){
    total.time.2 <- n*delta.t
    speed.2 <- velocity.x[2]/(1+k.all[2])
    spread.2 <- sd(x[2,])*(1+k.all[2])/(60*velocity.x[2])
  }
  
  if(mu.1<= 20){
    time.1.20 <- n*delta.t
    speed.1.20 <- velocity.x[1]/(1+k.all[1])
    spread.1.20 <- sd(x[1,])*(1+k.all[1])/(60*velocity.x[1])
  }
  
  if(mu.2<= 20){
    time.2.20 <- n*delta.t
    speed.2.20 <- velocity.x[2]/(1+k.all[2])
    spread.2.20 <- sd(x[2,])*(1+k.all[2])/(60*velocity.x[2])
  }
  
  if(mu.1<= col.length/4){
    time.1.fourth <- n*delta.t
    speed.1.fourth <- velocity.x[1]/(1+k.all[1])
    spread.1.fourth <- sd(x[1,])*(1+k.all[1])/(60*velocity.x[1])
  }
  
  if(mu.2<= col.length/4){
    time.2.fourth <- n*delta.t
    speed.2.fourth <- velocity.x[2]/(1+k.all[2])
    spread.2.fourth <- sd(x[2,])*(1+k.all[2])/(60*velocity.x[2])
  }
  
  if(mu.1<= col.length/2){
    time.1.half <- n*delta.t
    speed.1.half <- velocity.x[1]/(1+k.all[1])
    spread.1.half <- sd(x[1,])*(1+k.all[1])/(60*velocity.x[1])
  }
  
  if(mu.2<= col.length/2){
    time.2.half <- n*delta.t
    speed.2.half <- velocity.x[2]/(1+k.all[2])
    spread.2.half <- sd(x[2,])*(1+k.all[2])/(60*velocity.x[2])
  }
  
  #  1/(1+exp( h/(R*( T.all )) + C.all )
  
  if(n%%500==0){
    if(max(x[1,],na.rm=TRUE)>col.length){
      all.detected.1 <- which(x[1,] > col.length)
      new.detected.1 <- all.detected.1[!all.detected.1 %in% detected.1]
      detected.1 <- c(detected.1, new.detected.1)
      detector[1, new.detected.1] <- n*delta.t          
    }
    if(max(x[2,],na.rm=TRUE)>col.length){
      all.detected.2 <- which(x[2,] > col.length)
      new.detected.2 <- all.detected.2[!all.detected.2 %in% detected.2]
      detected.2 <- c(detected.2, new.detected.2)
      detector[2, new.detected.2] <- n*delta.t
      if(length(detected.2)+length(detected.1)-2 >= 2*n.mol){ break() }
    }     
  }
  
  if(n%%pause.count==0){
    if(min(x)>col.length) break()
    vel1_x <- c(vel1_x , mean(velocity.x[1:100]) )
    vel2_x <- c(vel2_x , mean(velocity.x[101:200]) )
    pos1_x <-c(pos1_x , mean(x[1,]) )
    pos2_x <-c(pos1_x , mean(x[2,]) )
    std1 <- c(std1 , sd(x[1,]))
    std2 <- c(std2 , sd(x[2,]))
    move.count <- move.count + 1
    big.matrix[move.count,] <- x
    top.pks <- pk.head <- pk.tail <- numeric(j)
    peakx <- numeric(1)
  } 
}  #End of big.matrix for loop.

##big.matrix <- big.matrix[1:move.count,]





#### End Standard Execution










### STart picture drawing

cols <- c('red','green','blue','grey30')
Tmax <- delta.t * n * tp.rate + T.ramp + T.0
Tmin <- T.0
Tdelta <- Tmax-Tmin
peak_width <-rep(0,move.count)
plot_res <- rep(0,move.count)
time = rep(0,move.count)

for(mm in 1:move.count){
  png(paste(mm,"sec.png"))
  par(mfrow=c(2,1), mar=c(2,5,1,6), las=1)
  plot(0, 0, type='n', xlim=c(0,col.length), ylim=c(Tmin, Tmax), ylab='Temperature (K)')#
  mtext(text=sprintf("Time:\n %3.2f s",(mm-1)*pause.count*delta.t), side=4, cex=2)
  x <- matrix(big.matrix[mm,], nrow=j,ncol=n.mol)
  
  for(i in 1:j){
    points(x[i,], (1:n.mol)/n.mol*Tdelta+Tmin, col=cols[i], pch=i, cex=0.5) 
    text(col.length*0.33*i, Tmin+Tdelta*0.2, labels=round(4*sd(x[i,]),2 ), col=cols[i],pos=2,cex=2)
    
  }
  tmp <- (sd(x[1,],2)+sd(x[2,],2))/2
  peak_width[mm]<-4*tmp
  plot_res[mm]<-(mean(x[1,])-mean(x[2,]))/(4*tmp)
  time[mm]<-mm
  lines(xx <-seq(1, col.length, 0.1),  (Temp.vec[xx*Temp.granularity]+tp.rate*delta.t*(mm-1)*pause.count+T.0),col='red')
  #abline(h=( (delta.t*mm*pause.count*b) )/200,col='red') #programmed temperature     
  plot(density(x[1,]), type='n', xlim=c(0,col.length), yaxt='n', xlab='',ylab='Density',main='')
  for(i in 1:j){ lines(density(x[i,]), col=cols[i] ) }
  dev.off()
} 

##write.csv(peak_width,file="peak_width.csv",row.names=FALSE)
write.csv(plot_res,file="resolutionR.csv",row.names=FALSE)
write.csv(time,file="timeR.csv",row.names=FALSE)
write.csv(vel1_x,file="vel1_xR.csv",row.names=FALSE)
write.csv(vel2_x,file="vel2_xR.csv",row.names=FALSE)
write.csv(pos1_x,file="pos1_xR.csv",row.names=FALSE)
write.csv(pos2_x,file="pos2_xR.csv",row.names=FALSE)
write.csv(std1,file="std1R.csv",row.names=FALSE)
write.csv(std2,file="std2R.csv",row.names=FALSE)

##### End picture drawing



retention.time <- rowMeans(detector) 
stand.dev <- apply(detector,1,sd) 
resolutions <- rep(0,j) 
for(i in 2:j){resolutions[i] <- (retention.time[i]-retention.time[i-1])/(2*(stand.dev[i]+stand.dev[i-1]))} 
(rbind(stand.dev,retention.time, resolutions)) 




##### Simulation Results

plate.len <-  30 #cm
decay.const <- 0.05
T.ramp <- 90 #height of the temperature ramp; delta T across plate
n.corners <- 20 
#stand.dev       0.02365578   0.02487849
#retention.time 40.46400000  51.47850000
# resolutions     0.00000000 113.47136009


plate.len <-  30 #cm
decay.const <- 0.1
# stand.dev       0.02480225   0.01305582
# retention.time 52.17100000  70.30250000
# resolutions     0.00000000 239.46676175



## Double the length, but same number of turns shold make the gradient portions longer and isothermal shorter
plate.len <-  60 #cm   
decay.const <- 0.1  #first half same as plate.len=30, decay.const=0.1, but extends to 60 cm
# stand.dev       0.0223098   0.01507557
# retention.time 64.5135000  88.05500000
# resolutions     0.0000000 314.84910248



plate.len <-  60 #cm
decay.const <- 0.05  # Same shape as plate.len= 30, decay.const= 0.1
#stand.dev       0.02480225 0.008689364
#retention.time 52.12900000 70.24950
#resolutions     0.00000000 270.5230

plate.len <- 120 #cm
decay.const <- 0.025 # same shape as 30;0.1, and 60,0.05.
#stand.dev       0.01882938   0.01930615
# retention.time 52.05700000  70.14100000
# resolutions     0.00000000 237.10176732

