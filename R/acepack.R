
"ace"<- function(x, y, wt = rep(1, nrow(x)), cat=NULL, mon=NULL, lin=NULL, circ=NULL, delrsq = 0.01)
{
	x <- as.matrix(x)
	if(delrsq <= 0) {
		cat("delrsq must be positive")
		return()
	}
	iy <- ncol(x) + 1
	l <- matrix(1, ncol = iy)
	if(!is.null(circ)) {
		for(i in 1:length(circ)) {
			if(circ[i] < 0 || circ[i] > nrow(x)) {
				cat("bad circ= specification")
				return()
			}
			if(circ[i] == 0) {
				cat(
					"response spec can only be lin or ordered (default)"
					)
				return()
			}
			else {
				nncol <- circ[i]
				if(l[nncol]!=2 & l[nncol]!=1) {
					cat(
						"conflicting transformation specifications"
						)
					return()
				}
				l[nncol] <- 2
			}
		}
	}
	if(!is.null(mon)) {
		for(i in 1:length(mon)) {
			if(mon[i] < 0 || mon[i] > nrow(x)) {
				cat("bad mon= specification")
				return()
			}
			if(mon[i] == 0) {
				cat(
					"response spec can only be lin or ordered (default)"
					)
				return()
			}
			else {
				nncol <- mon[i]
				if(l[nncol]!=3 && l[nncol]!=1) {
					cat(
						"conflicting transformation specifications"
						)
					return()
				}
				l[nncol] <- 3
			}
		}
	}
	if(!is.null(lin)) {
		for(i in 1:length(lin)) {
			if(lin[i] < 0 || lin[i] > nrow(x)) {
				cat("bad lin= specification")
				return()
			}
			if(lin[i] == 0) {
				nncol <- iy
			}
			else {
				nncol <- lin[i]
			}
			if(l[nncol]!=4 && l[nncol]!=1) {
				cat("conflicting transformation specifications"
					)
				return()
			}
			l[nncol] <- 4
		}
	}
	if(!is.null(cat)) {
		for(i in 1:length(cat)) {
			if(cat[i] < 0 || cat[i] > nrow(x)) {
				cat("bad cat= specification")
				return()
			}
			if(cat[i] == 0) {
				cat(
					"response spec can only be lin or ordered (default)"
					)
				return()
			}
			else {
				nncol <- cat[i]
				if(l[nncol]!=4 && l[nncol]!=1) {
					cat(
						"conflicting transformation specifications"
						)
					return()
				}
				l[nncol] <- 4
			}
		}
	}
	tx <- x
	ty <- y
	m <- matrix(0, nrow = nrow(x), ncol =  iy)
	z <- matrix(0, nrow = nrow(x), ncol = 12)
	z <- as.matrix(z)
	ns<-1
	mode(ns)<-"integer"
	mode(x) <- "double"
	mode(y) <- "double"
	mode(tx) <- "double"
	mode(ty) <- "double"
	mode(wt) <- "double"
	mode(m) <- "integer"
	mode(l) <- "integer"
	mode(delrsq) <- "double"
	mode(z) <- "double"
	junk <- .Fortran("mace",
			 p=ncol(x),
			 n=nrow(x),
			 x=t(x),
			 y=y,
			 w=wt,
			 l=l,
			 delrsq=delrsq,
			 ns=ns,
			 tx = tx,
			 ty = ty,
			 rsq = double(1),
			 ierr = integer(1),
			 m=m,
			 z=z
			 )
#	return(list(x=x, y=y, tx = junk$tx, ty = junk$ty, rsq = junk$rsq, l, m))
return(junk)
}


avas.formula <-function(formula, data=sys.frame(sys.parent()) ,subset=NULL,weights=NULL,cat=NULL,mon=NULL,lin=NULL,circ=NULL,delrsq=0.01,yspan=0){
  mf<-match.call()
  mf$cat<-mf$mon<-lin<-mf$circ<-mf$delrsq<-mf$span<-NULL
  mf[[1]]<-as.name("model.frame")
  mf$use.data<-TRUE
  mf<-eval(mf,sys.frame(sys.parent()))
  y<-model.response(mf)
  mm<-model.matrix(formula,mf)[,-1,drop=F]
  w<-if (is.null(weights)) rep(1,nrow(mm)) else model.extract(mf,"weights")
  avas(mm,y,w,cat,mon,lin,circ,delrsq,yspan)
}

"avas"<- function(x, y, wt = rep(1, nrow(x)), cat=NULL, mon=NULL, lin=NULL, circ=NULL, delrsq = 0.01, yspan
	 = 0)
{
	x <- as.matrix(x)
	if(delrsq <= 0) {
		cat("delrsq must be positive")
		return()
	}
	iy <- ncol(x) + 1
	l <- matrix(1, ncol = iy)
	if(!is.null(circ)) {
		for(i in 1:length(circ)) {
			if(circ[i] < 0 || circ[i] > nrow(x)) {
				cat("bad circ= specification")
				return()
			}
			if(circ[i] == 0) {
				cat(
					"response spec can only be lin or ordered (default)"
					)
				return()
			}
			else {
				nncol <- circ[i]
				if(l[nncol]!=2 & l[nncol]!=1) {
					cat(
						"conflicting transformation specifications"
						)
					return()
				}
				l[nncol] <- 2
			}
		}
	}
	if(!is.null(mon)) {
		for(i in 1:length(mon)) {
			if(mon[i] < 0 || mon[i] > nrow(x)) {
				cat("bad mon= specification")
				return()
			}
			if(mon[i] == 0) {
				cat(
					"response spec can only be lin or ordered (default)"
					)
				return()
			}
			else {
				nncol <- mon[i]
				if(l[nncol]!=3 && l[nncol]!=1) {
					cat(
						"conflicting transformation specifications"
						)
					return()
				}
				l[nncol] <- 3
			}
		}
	}
	if(!is.null(lin)) {
		for(i in 1:length(lin)) {
			if(lin[i] < 0 || lin[i] > nrow(x)) {
				cat("bad lin= specification")
				return()
			}
			if(lin[i] == 0) {
				nncol <- iy
			}
			else {
				nncol <- lin[i]
			}
			if(l[nncol]!=4 && l[nncol]!=1) {
				cat("conflicting transformation specifications"
					)
				return()
			}
			l[nncol] <- 4
		}
	}
	if(!is.null(cat)) {
		for(i in 1:length(cat)) {
			if(cat[i] < 0 || cat[i] > nrow(x)) {
				cat("bad cat= specification")
				return()
			}
			if(cat[i] == 0) {
				cat(
					"response spec can only be lin or ordered (default)"
					)
				return()
			}
			else {
				nncol <- cat[i]
				if(l[nncol]!=4 && l[nncol]!=1) {
					cat(
						"conflicting transformation specifications"
						)
					return()
				}
				l[nncol] <- 4
			}
		}
	}
	tx <- x
	ty <- y
	m <- matrix(0, nrow = nrow(x), ncol = ncol(x) + 2)
	z <- matrix(0, nrow = nrow(x), ncol = 17)
	z <- as.matrix(z)
	iters <- matrix(0, nrow = 100, ncol = 2)
	mode(x) <- "double"
	mode(y) <- "double"
	mode(tx) <- "double"
	mode(ty) <- "double"
	mode(wt) <- "double"
	mode(m) <- "integer"
	mode(l) <- "integer"
	mode(delrsq) <- "double"
	mode(z) <- "double"
	mode(yspan) <- "double"
	mode(iters) <- "double"
	junk <- .Fortran("avas",
		ncol(x),
		nrow(x),
		x,
		y,
		wt,
		l,
		delrsq,
		tx = tx,
		ty = ty,
		rsq = double(1),
		ierr = integer(1),
		m,
		z,
		yspan = yspan,
		niter = integer(1),
		iters = iters)
	junk$iters <- junk$iters[1:junk$niter,  ]
	return(list(x=x, y=y, tx = junk$tx, ty = junk$ty, rsq = junk$rsq, l, m, yspan = 
		junk$yspan, iters = junk$iters, niters = junk$niter))
}














