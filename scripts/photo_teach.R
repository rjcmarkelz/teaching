library(plantecophys)

?Photosyn()
r <- Photosyn(VPD=seq(0.5, 4, length=25), Vcmax=50, Jmax=100)



with(r, plot(VPD, ALEAF, type='l'))
p <- Aci(Ci = 100, Ca=400)
p
with(p, plot(Ci, ALEAF, type='l'))

head(iris)

pl2 <- ggplot(data = photo_out(), aes(x = Ci, y = ALEAF))
pl2 <- pl2 + geom_line(size = 2) + ylim(-3, 25)
pl2