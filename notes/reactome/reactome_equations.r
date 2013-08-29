


x = (-500:500) / 100
yp = 1 + 1 / (1 + exp(-x))
yn = 1 - 1 / (1 + exp(-x))
yr = 1 / (1 + exp(-x))
yc = 2 / (1 + exp(-x))

#par(mfrow=c(2,1))
plot(x,y, type='l', ylim=c(0,2), col='green', xlab='input', ylab='output', main="Reaction", sub="GREEN: Positive reg.    RED:Negative reg.    BLUE: Catalyst    BLACK: Requirement" )
lines(x,yn, type='l', col='red')
lines(x,yr, type='l', col='black')
lines(x,yc, type='l', col='blue')

