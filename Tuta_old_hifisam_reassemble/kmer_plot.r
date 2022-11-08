
path <- getwd()

k21 <- read.table(paste0(path, "/tuab_21mers.histo"), header=F, sep="\t")
k27 <- read.table(paste0(path, "/tuab_27mers.histo"), header=F, sep="\t")
k31 <- read.table(paste0(path, "/tuab_31mers.histo"), header=F, sep="\t")
k31$V2 <- k31$V2/2 
# make k21 count (x axis) compatible to the other two
k21_trim <- k21[trimws(k21$V1) %in% trimws(k27$V1), ]

x11() # 
plot(k21_trim[2:100,],
     type="l",
     xlab="kmer count",
     ylab="frequency",
     ylim=range(1,max(k31$V2))) # plots the line from the data points 
lines(k27[2:100,], col="red")
lines(k31[2:100,], col="blue")





# replot the histogram
highest <- myhisto[which(myhisto$V2 == max(myhisto[10:200,2])),]

x11() # this only works for windows R I guess
plot(myhisto[10:200,],
     type="l",
     xlab="kmer count",
     ylab="frequency") # plots the line from the data points 
points(myhisto[10:200,]) # plot the data points
text(highest[,1], highest[,2], highest$V1,
     cex=1, pos=2,col="black")