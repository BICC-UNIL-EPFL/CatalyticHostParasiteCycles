# R code to generate the figures of the manuscript
# Energy-harnessing problem solving of primordial life:
#  modeling the emergence of cata-lytic host-nested parasite life cycles
# Copyright (c) 2022, Christian Iseli, Bernard Conrad and Magnus Pirovino
#
### setup
dt = 0.01
m = 5997
E_lower = 0.5
E_upper = 0.9
###Take values of E from Excel file “_1 Figs 1a-b Virulence step & smoothstep 050521”,
###Sheet “Calculations”, Column E, Row 9:6006****
data <- read.table('_1_E_values.txt',header=TRUE,sep='\t',skip=7)
E <- data$E..lhs.
Step_virulence <- rep(0, m)
from_below = FALSE
for (k in 1:m) {
  below = E[k] < E_lower
  in_between = (E[k] > E_lower) && (E[k] < E_upper)
  if (from_below && (below || in_between)) {
    Step_virulence[k] = 1
  } else {
    # already set by default
    Step_virulence[k] = 0
  }
  if (from_below && E[k] < E_upper) {
    from_below = TRUE
  } else if (E[k] < E_lower) {
    from_below = TRUE
  } else {
    from_below = FALSE
  }
}

y1_min = -.2
y1_max = 2.3
y2_min = -.25
y2_max = 1.15

pdf(file='_1_Fig1a.pdf',width=8.3,height=5)

## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 6) + 0.1)

## Plot first set of data and draw its axis
plot(E, axes=FALSE, ylim=c(y1_min,y1_max), xlab="", ylab="", type="b", col="green", main="Figure 1a")
box()
abline(h=E_lower, col="darkgreen", lty=3, lwd=2)
abline(h=E_upper, col="darkgreen", lty=2, lwd=2)
axis(2, ylim=c(y1_min,y1_max),col="green",col.axis="green",las=1)  ## las=1 makes horizontal labels
mtext("Energy",side=2,line=2.5,col="green")

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(Step_virulence, pch=15,  xlab="", ylab="", ylim=c(y2_min,y2_max), axes=FALSE, type="b", col="red")
## a little farther out (line=4) to make room for labels
mtext("Step-Virulence",side=4,col="red",line=4) 
axis(4, at=c(0,1), ylim=c(y2_min,y2_max), col="red",col.axis="red",las=1)

## Draw the time axis
axis(1,pretty(c(1,length(E)),10))
mtext("Time (Steps)",side=1,col="black",line=2.5)  

## Add Legend
legend("bottomleft",horiz=T,legend=c("Energy","Virulence","E_lower","E_upper"),
  text.col=c("green","red","darkgreen","darkgreen"),
  col=c("green","red","darkgreen","darkgreen"),lty=c(1,1,3,2),lwd=c(4,4,2,2))

dev.off()

###
