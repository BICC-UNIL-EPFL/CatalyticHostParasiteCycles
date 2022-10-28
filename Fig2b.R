# R code to generate the figures of the manuscript
# Energy-harnessing problem solving of primordial life:
#  modeling the emergence of cata-lytic host-nested parasite life cycles
# Copyright (c) 2022, Christian Iseli, Bernard Conrad and Magnus Pirovino
#
# LVhpr one energy model with instability general program code:
#begin ***LVhpr one energy model with instability
dt = 0.1
m = 5997
alpha = 0.0125
beta = 0.05
gamma = 0.2
a_1 = 0.05
a_0 = 0.2
b_1 = 0.12
b_0 = 0.4
H = 4
q = 10
kappa = 1
sigma = 1
E_lower = 0.5
E_upper = 2.5
from_below = TRUE
from_above = FALSE
Parasite_eqilibrium = b_1/(b_0*sigma)*H  ### equilibrium X
Hyperparasite_equilibrium = a_1*q/(b_0*sigma)*(b_0*sigma-b_1/kappa)/(q*a_0*sigma+a_1/kappa)*H  ### equilibrium Y
### In R, indices start at 1...
X <- rep(0, m+1)
X[1] = 0.9 ###  Parasite initial population
Y <- rep(0, m+1)
Y[1] = 0.7 ###  Hyperparasite initial population
E <- rep(0, m+1)
E[1]=0.6
Virulence <- rep(0, m+1)
Virulence[1] = 1
t_gamma_eq_0 = 300 ### starting time of fitness erosion
gamma_erosion_rate = 0.01
erosion_halftime = -log(1/2)/gamma_erosion_rate
null_energy_parasite_fitness = -0.01
t = 0

# so we go now from 2 to (m+1) ... yes, needs parenthesis otherwise +1 is applied to the vector 2:m...
for (k in 2:(m+1)) {
  t = t+dt
  E[k] = max(E[k-1]+dt*(-alpha*H-beta*Virulence[k-1]*X[k-1]+gamma*Virulence[k-1]*Y[k-1]),0)
  habitatrestrictionterm = 1-(X[k-1]+Y[k-1]/q)/(kappa*H)
  below = (E[k] < E_lower)
  in_between = (E[k] > E_lower) && (E[k] < E_upper)
  if (from_below && (below || in_between)) {
    Virulence[k] = 1
  } else {
    Virulence[k] = 0
  }
  if (from_below && E[k] < E_upper) {
    from_below = TRUE
  } else if (E[k] < E_lower) {
    from_below = TRUE
  } else {
    from_below = FALSE
  }
  from_above = ! from_below

  if (t < t_gamma_eq_0) {
    X[k] = X[k-1] + dt*Virulence[k-1]*X[k-1]*(a_1*habitatrestrictionterm-a_0*sigma*Y[k-1]/H)
    Y[k] = Y[k-1] + dt*Virulence[k-1]*Y[k-1]*(-b_1+b_0*sigma*X[k-1]/H)
  } else  {
    gamma = gamma*(1-dt*gamma_erosion_rate)
    if (E[k] > 0) {
      t_last_energy = t  # t_last_energy is the last time when E was still positive
      eroded_fitness = a_1*habitatrestrictionterm
    } else {
      eroded_fitness = ((t-t_last_energy)*null_energy_parasite_fitness
                        +(erosion_halftime-(t-t_last_energy))*a_1*habitatrestrictionterm)/erosion_halftime
    }
    X[k] = X[k-1] + dt*Virulence[k-1]*X[k-1]*(eroded_fitness - a_0*sigma*Y[k-1]/H)
    Y[k] = Y[k-1] + dt*Virulence[k-1]*Y[k-1]*(-b_1+b_0*sigma*X[k-1]/H)
  }
} ### for k
### end LVhpr one energy model

# -	For Figs 2a&b: Could you find a solution with the y-axis in the same way you did for Figs 3? Left axis is Energy and right axis is Virulence, but we should also have an axis-scale for Population Size? 
# -	For Fig 2b can we have text inside the graph as recommended in the word document attached?: 
# o	“start of hyperparasite degeneration” and
# o	“Erosion of hyperparasite energy catalysis” 
# Or should or can this be done by Bernard later on top of his document?

y1_min = -.2
y1_max = 2.6
y2_min = -.25
y2_max = 1.15
y3_min = -.1
y3_max = 1.5

pdf(file='_2_Fig2b.pdf',width=14,height=6)

## add extra space to right margin of plot within frame
par(mar=c(5, 9, 4, 6) + 0.1)

## Plot first set of data and draw its axis
plot(E, axes=FALSE, ylim=c(y1_min,y1_max), xlab="", ylab="", type="b", col="green", main=paste("Figure 2b using step",dt))
box()
abline(h=E_lower, col="darkgreen", lty=3, lwd=2)
abline(h=E_upper, col="darkgreen", lty=2, lwd=2)
abline(v=t_gamma_eq_0/dt, col="darkgrey", lty=2, lwd=2)
arrows(x0=t_gamma_eq_0/dt,y0=2.6,x1=t_gamma_eq_0/dt+1500,y1=2.6,col="darkgrey", length=.125)
axis(2, ylim=c(y1_min,y1_max),col="green",col.axis="green",las=1)  ## las=1 makes horizontal labels
mtext("Energy",side=2,line=2.5,col="green")

## Allow a second plot on the same graph
par(new=T)
plot(X, ylim=c(y3_min,y3_max), axes=FALSE, xlab="", ylab="", type="l", col="blue", main="")
points(Y, type="l", lwd=2, col="blue2")
axis(2, ylim=c(y3_min,y3_max), col="blue", col.axis="blue", las=1, line=4)
mtext(2, text="Populations", line=6.5, col="blue4")

par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(Virulence, pch=15,  xlab="", ylab="", ylim=c(y2_min,y2_max), axes=FALSE, type="b", col="red")
## a little farther out (line=4) to make room for labels
mtext("Step-Virulence",side=4,col="red",line=3)
axis(4, at=c(0,1), ylim=c(y2_min,y2_max), col="red",col.axis="red",las=1)

## Draw the time axis
axis(1,pretty(c(1,length(E)),10))
mtext("Time (Steps)",side=1,col="black",line=2.5)

## Add Legend
legend("bottomleft",horiz=T,legend=c("Energy","Virulence","Parasite","Hyperparasite","E_lower","E_upper"),
  text.col=c("green","red","blue", "blue2","darkgreen","darkgreen"),
  col=c("green","red","blue", "blue2","darkgreen","darkgreen"),lty=c(1,1,1,1,3,2),lwd=c(4,4,1,2,2,2))

dev.off()

###
