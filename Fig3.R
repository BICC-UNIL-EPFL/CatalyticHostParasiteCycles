# R code to generate the figures of the manuscript
# Energy-harnessing problem solving of primordial life:
#  modeling the emergence of cata-lytic host-nested parasite life cycles
# Copyright (c) 2022, Christian Iseli, Bernard Conrad and Magnus Pirovino
#
# LVhpr two energies model general program code:
pdf(file='_3_Fig3.pdf',width=15,height=24)

# 4 plots on the same page
layout(matrix(c(1,2,3,4),nrow=4,ncol=1))

# these values should no change between steps
dt = 1
H = 2.0
q = 10.0
kappa = 3
sigma = 1.1
theta = 0.1
alpha1 = 0.001
alpha2 = 0.000001
gamma1 = 0.0004
gamma2 = 0.015
gamma3 = 0.0045
a_1 = 0.001
a_0 = 0.001
b_1 = 0.001
b_0 = 0.0011
c_1 = 0.0013
c_0 = 0.0013
d_1 = 0.001
d_0 = 0.00143
e_1 = 0.0013
e_0 = 0.0013
f_1 = 0.001
f_0 = 0.00143
E1_lower = 20
E1_upper = 30
E2_lower = 3
E2_upper = 3.5

for (totaltime in list(90000, 500000)) {
  # begin ### LVhpr two energies model
  ### for Figs 3 (dt=3.0 for Figs 3’)
  m = floor(totaltime/dt) ### 500'000 for Figs 3 (m=90’500 for Figs 3’)
  ### In R, indices start at 1...
  X <- matrix(rep(0, 3*(m+1)), ncol=3, nrow=m+1)
  Y <- matrix(rep(0, 3*(m+1)), ncol=3, nrow=m+1)
  X[1,1] = 1.0  ### Parasite X1 initial population
  Y[1,1] = 0.2  ###  Hyperparasite Y1 initial population
  X[1,2] = 0.5  ### Parasite X2 initial population
  Y[1,2] = 3.0  ###  Hyperparasite Y2 initial population
  X[1,3] = 0.7  ### Parasite X3 initial population
  Y[1,3] = 3.2  ###  Hyperparasite Y3 initial population
  E1 <- rep(0, m+1)
  E1[1] = 21.0 ###  initial Energy E1
  vE1 <- rep(0, m+1)
  vE1[1] = 0.0 ###  initial virulence vE1
  from_below1 = FALSE
  from_above1 = TRUE
  E2 <- rep(0, m+1)
  E2[1] = 2.6  ###  initial Energy E2
  vE2 <- rep(0, m+1)
  vE2[1] = 1.0 ###  initial virulence vE2
  from_below2 = TRUE
  from_above2 = FALSE

  t = 0
  # so we go now from 2 to (m+1) ... yes, needs parenthesis otherwise +1 is applied to the vector 2:m...
  for (k in 2:(m+1)) {
    t = t+dt
    habitatrestrictionterm = 1-(X[k-1,1]+X[k-1,2]+X[k-1,3] +1/q*(Y[k-1,1]+Y[k-1,2]+Y[k-1,3]))/(kappa*H)
    X[k,1] = X[k-1,1] + dt*vE1[k-1]*X[k-1,1]*(a_1*habitatrestrictionterm-a_0*sigma*Y[k-1,1]/H)
    Y[k,1] = Y[k-1,1] + dt*vE1[k-1]*Y[k-1,1]*(-b_1+b_0*sigma* X[k-1,1]/H)
    X[k,2] = X[k-1,2] + dt*vE2[k-1]*X[k-1,2]*(c_1*habitatrestrictionterm-c_0*sigma*Y[k-1,2]/H)
    Y[k,2] = Y[k-1,2] + dt*vE2[k-1]*Y[k-1,2]*(-d_1+d_0*sigma* X[k-1,2]/H)
    X[k,3] = X[k-1,3] + dt*vE2[k-1]*X[k-1,3]*(e_1*habitatrestrictionterm-e_0*sigma*Y[k-1,3]/H)
    Y[k,3] = Y[k-1,3] + dt*vE2[k-1]*Y[k-1,3]*(-f_1+f_0*sigma* X[k-1,3]/H)

    E1[k] = E1[k-1]+dt*(-alpha1*H + gamma1/theta*vE1[k-1]*Y[k-1,1]+gamma2*vE2[k-1]*Y[k-1,2]-gamma3*vE2[k-1]*Y[k-1,3])
    E2[k] = E2[k-1]+dt*(-alpha2*E2[k-1] - gamma1*vE1[k-1]*Y[k-1,1]+gamma3*theta*vE2[k-1]*Y[k-1,3])

    below1 = (E1[k] < E1_lower)
    in_between1 = (E1[k] > E1_lower) && (E1[k] < E1_upper)
    if (from_below1 && (below1 || in_between1)) {
      vE1[k] = 1
    } else {
      vE1[k] = 0
    }
    if (from_below1 && E1[k] < E1_upper) {
      from_below1 = TRUE
    } else if (E1[k] < E1_lower) {
      from_below1 = TRUE
    } else {
      from_below1 = FALSE
    }
    from_above1 = ! from_below1

    below2 = (E2[k] < E2_lower)
    in_between2 = (E2[k] > E2_lower) && (E2[k] < E2_upper)
    if (from_below2 && (below2 || in_between2)) {
      vE2[k] = 1
    } else {
      vE2[k] = 0
    }
    if (from_below2 && E2[k] < E2_upper) {
      from_below2 = TRUE
    } else if (E2[k] < E2_lower) {
      from_below2 = TRUE
    } else {
      from_below2 = FALSE
    }
    from_above2 = ! from_below2
  } ### for k
  ### LVhpr two energies model

  # If possible, we would suggested the following:

  yE2_min = 0
  yE2_max = 4.5
  y1_min = min(c(min(X),min(Y)))
  y1_max = max(c(max(X),max(Y)))
  y2_min = -.25
  y2_max = 1.15

  # Fig 3.1: Virulence1, Virulence 2, E1, E2, as well as the dotted ones E1_upper, E1_lower, E2_upper, E2_lower
  # -	Virulence … different reds
  # -	E… different greens

    par(mar=c(5, 12, 4, 5) + 0.1)

    ## Plot first set of data and draw its axis
    plot(E1, ylim=c(0,max(E1)), axes=FALSE, xlab="", ylab="", type="l", col="green", main=paste("Energies and virulences using step",dt,"for total time of",totaltime))
    axis(2, ylim=c(0,max(E1)), col="green", col.axis="green", las=1)  ## las=1 makes horizontal labels
    mtext(2, text="Energy 1", line=2, col="green")
    box()
    abline(h=E1_lower, col="green", lty=3, lwd=2)
    abline(h=E1_upper, col="green", lty=2, lwd=2)

    ## Allow a second plot on the same graph
    par(new=T)
    plot(E2, ylim=c(yE2_min,yE2_max), axes=FALSE, xlab="", ylab="", type="l", col="aquamarine3", main="")
    axis(2, ylim=c(yE2_min,yE2_max), col="aquamarine3", col.axis="aquamarine3", las=1, line=3.5)
    mtext(2, text="Energy 2", line=5.5, col="aquamarine3")
    abline(h=E2_lower, col="aquamarine3", lty=3, lwd=2)
    abline(h=E2_upper, col="aquamarine3", lty=2, lwd=2)

    par(new=TRUE)
    ## Plot the virulence and put axis scale on right
    plot(vE1, xlab="", ylab="", ylim=c(-2,18), axes=FALSE, type="l", col="red")
    ## a little farther out (line=4) to make room for labels
    mtext("Step-Virulence 1",side=4,col="red",line=1)
    axis(4, at=c(0,1), ylim=c(-2,18), col="red",col.axis="red",las=1)

    par(new=TRUE)
    ## Plot the virulence and put axis scale on right
    plot(vE2, xlab="", ylab="", ylim=c(-3.5,16.5), axes=FALSE, type="l", col="pink")
    ## a little farther out (line=4) to make room for labels
    mtext("Step-Virulence 2",side=4,col="pink",line=2.5)
    axis(4, at=c(0,1), ylim=c(-3.5,16.5), col="pink",col.axis="pink",las=1)

    ## Draw the time axis
    axis(1,pretty(c(1,length(E1)),10))
    mtext("Time (Steps)",side=1,col="black",line=2.5)

    ## Add Legend
    legend("bottomleft",horiz=T,
      legend=c("Energy 1","Energy 2","Virulence 1","Virulence 2","E1_lower","E1_upper","E2_lower","E2_upper"),
      text.col=c("green","aquamarine3","red","pink","green","green","aquamarine3","aquamarine3"),
      col=c("green","aquamarine3","red","pink","green","green","aquamarine3","aquamarine3"),lty=c(1,1,1,1,3,2,3,2),lwd=c(1,1,1,1,2,2,2,2))

  # Fig 3.2: Parasite1, Hyperparasite1, Virulence1, Virulence 2, E1, E2 (without E1-2_upper-lower) 
  # -	Parasite, Hyperparasite… different blues
  # -	Virulence … different reds
  # -	E… different greens
  # Fig 3.3: Parasite2, Hyperparasite2, Virulence1, Virulence 2, E1, E2 (without E1-2_upper-lower) 
  # -	Parasite, Hyperparasite… different blues
  # -	Virulence … different reds
  # -	E… different greens
  # Fig 3.4: Parasite3, Hyperparasite3, Virulence1, Virulence 2, E1, E2 (without E1-2_upper-lower) 
  # -	Parasite, Hyperparasite… different blues
  # -	Virulence … different reds
  # -	E… different greens

  for (i in 1:3) {
    ## add extra space to right margin of plot within frame
    par(mar=c(5, 12, 4, 5) + 0.1)

    ## Plot first set of data and draw its axis
    plot(E1, ylim=c(0,max(E1)), axes=FALSE, xlab="", ylab="", type="l", col="green", main=paste("Parasite group",i,"using step",dt,"for total time of",totaltime))
    axis(2, ylim=c(0,max(E1)), col="green", col.axis="green", las=1)  ## las=1 makes horizontal labels
    mtext(2, text="Energy 1", line=2, col="green")
    box()

    ## Allow a second plot on the same graph
    par(new=T)
    plot(E2, ylim=c(yE2_min,yE2_max), axes=FALSE, xlab="", ylab="", type="l", col="aquamarine3", main="")
    axis(2, ylim=c(yE2_min,yE2_max), col="aquamarine3", col.axis="aquamarine3", las=1, line=3.5)
    mtext(2, text="Energy 2", line=5.5, col="aquamarine3")

    par(new=T)
    plot(X[,i], ylim=c(y1_min,y1_max), axes=FALSE, xlab="", ylab="", type="l", col="blue", main="")
    points(Y[,i], type="l", lwd=2, col="blue2")
    axis(2, ylim=c(y1_min,y1_max), col="blue", col.axis="blue", las=1, line=7.5)
    mtext(2, text="Populations", line=10, col="blue4")

    par(new=TRUE)
    ## Plot the virulence and put axis scale on right
    plot(vE1, xlab="", ylab="", ylim=c(-2,18), axes=FALSE, type="l", col="red")
    ## a little farther out (line=4) to make room for labels
    mtext("Step-Virulence 1",side=4,col="red",line=1)
    axis(4, at=c(0,1), ylim=c(-2,18), col="red",col.axis="red",las=1)

    par(new=TRUE)
    ## Plot the virulence and put axis scale on right
    plot(vE2, xlab="", ylab="", ylim=c(-3.5,16.5), axes=FALSE, type="l", col="pink")
    ## a little farther out (line=4) to make room for labels
    mtext("Step-Virulence 2",side=4,col="pink",line=2.5)
    axis(4, at=c(0,1), ylim=c(-3.5,16.5), col="pink",col.axis="pink",las=1)

    ## Draw the time axis
    axis(1,pretty(c(1,length(E1)),10))
    mtext("Time (Steps)",side=1,col="black",line=2.5)

    ## Add Legend
    legend("bottomleft",horiz=T,
      legend=c("Energy 1","Energy 2","Virulence 1","Virulence 2",paste("Parasite",i),paste("Hyperparasite",i)),
      text.col=c("green","aquamarine3","red","pink","blue","royalblue"),
      col=c("green","aquamarine3","red","pink","blue","blue2"),lty=c(1,1,1,1,1,1),lwd=c(1,1,1,1,1,2))
  }

}

dev.off()

###
