#install.packages("scatterplot3d") 	# uncomment this if not already installed...
library("scatterplot3d") 		# load graphics library

fname <- "..\\data\\henonmap.dat"		# data file name
data <- read.table(fname)				# load data

seq <- data$V1		# extract sequence
n <- length(seq)	# length
w <- 6				# word size
qmin <- -10.0
qmax <- 10.0
dq <- 0.1

qseq <- seq(qmin, qmax, by = dq)	# sequence of q's
H <- numeric(length(qseq))			# entropy
C <- numeric(length(qseq))			# complexity

dllname <- "gwpe64.dll"			# dll name
dyn.load(dllname)				# load the dll

ans <- .C("gwpe", as.numeric(H), as.numeric(C), as.numeric(seq), 
	as.integer(n), as.integer(w), as.numeric(qmin), as.numeric(qmax), as.numeric(dq))
H <- ans[[1]]
C <- ans[[2]]

dyn.unload(dllname)

scatterplot3d(H,C,qseq, zlab="q", pch = 16, color="steelblue")


