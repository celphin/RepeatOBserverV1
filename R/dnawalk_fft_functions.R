# RepeatOBserverV1
#-----------------------------
# functions

#' imagenan
#'
#' Displays a false colour image of two dimensional data set
#'
#' @param x 2D data to display
#' @param yline Number of lines in (left side) to where image begins
#' @param yma Number of characters in (left side)
#' @param col colour table
#' @param outside.below.color Colour below threshold
#' @param outside.above.color Colour above threshold
#' @param na.color Colour of NAs (see \code{\link[graphics]{par}})
#' @param ... other plot parameters eg zlim list with min and max range of color table data e.g. zlim=base::c(2,5)
#'
#' @return None
#'
#' @examples
#' d=base::data.frame(base::cbind(base::c(1:4),base::c(2,5,NA,NA) ))
#' imagenan(d)
#' @export
imagenan <- function(x,yline=3,yma=5,topyma=4, xline=3,xma=6,lnumr=10,lnumc=10,
                     lasval=1,cex.axis=0.7,cex.lab=0.7,
                     cex.main=0.7,widths=base::c(5,1.25), heights=base::c(1,0.5),
                     col = grDevices::topo.colors(255),outside.below.color='black',
                     outside.above.color='white',na.color='grey28',
                     ...){
  x<-base::as.matrix(x)
  .pardefault <-graphics::par(no.readonly = TRUE)
  reverse <- base::nrow(x) : 1
  x <- x[reverse,]
  zlim=base::c(base::min(x,na.rm=TRUE),base::max(x,na.rm=TRUE))
  if(zlim[2]<= zlim[1]+1e-10)zlim<-base::c(zlim[1],zlim[1]+1/100)  #  zlim[2]<= zlim[1]+1e-10
  if(!base::is.null(base::rownames(x)))rLabels <- base::rownames(x) else {rLabels<-base::c(reverse); base::rownames(x)<-base::c(reverse)}
  cLabels <- base::colnames(x)
  er<-base::nrow(x)
  ec<- base::ncol(x)
  main <-" "
  row_unit <-"Rows"
  col_unit<-"Columns"
  zunit<-"Intensity"
  rtick=1          # for plotting use start at 30 days for xtick marks
  ctick=1
  rtickinc=base::round(er/lnumr)    # for plotting use every 1 unit for "space" tick marks
  ctickinc=base::round(ec/lnumc)    # for plotting use every 1 unit for "space" tick marks
  if(rtickinc==0)rtickinc=1
  if(ctickinc==0)ctickinc=1
  rnames<-base::rownames(x)
  cnames<-base::colnames(x)
  attr=base::c(seq(rtick,er,by= rtickinc))
  attc=base::c(seq(ctick,ec,by= ctickinc))
  rlbls=rnames[seq(rtick,er,by= rtickinc)]
  clbls=cnames[seq(ctick,ec,by=ctickinc)]

  # check for additional function arguments
  if( base::length(base::list(...)) ){
    Lst <- base::list(...)
    if( !base::is.null(Lst$zlim) ){
      zlim<-base::c(Lst$zlim)
      if(zlim[2]<= zlim[1])zlim<-base::c(zlim[1],zlim[1]+1/100)
    }
    if( !base::is.null(Lst$yLabels) ){
      cLabels <- base::c(Lst$cLabels)
    }
    if( !base::is.null(Lst$xLabels) ){
      rLabels <- base::c(Lst$rLabels)
    }
    if( !base::is.null(Lst$main) ){
      main <- Lst$main
    }
    if( !base::is.null(Lst$cex.main) ){
      cex.main <- Lst$cex.main
    }
    if( !base::is.null(Lst$cex.axis) ){
      cex.axis <- Lst$cex.axis
    }
    if( !base::is.null(Lst$row_unit) ){
      row_unit <-  Lst$row_unit
    }
    if( !base::is.null(Lst$col_unit) ){
      col_unit <-  Lst$col_unit
    }
    if( !base::is.null(Lst$rlbls) ){
      rlbls <-  Lst$rlbls
    }
    if( !base::is.null(Lst$clbls) ){
      clbls <-  Lst$clbls
    }
    if( !base::is.null(Lst$zunit) ){
      zunit <-  Lst$zunit
    }
  }

  # check for null values
  if( base::is.null(rLabels) ){
    rLabels <- base::c(reverse)
  }
  if( base::is.null(cLabels) ){
    cLabels <- base::c(1: base::ncol(x))
  }

  graphics::layout(base::matrix(data=base::c(1,1,2,3), nrow=2, ncol=2), widths=widths, heights=heights)

  zstep <- (zlim[2] - zlim[1]) / base::length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA

  x[base::which(x<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  x[base::which(x>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  x[base::which(base::is.na(x>zlim[2]))] <- newz.na # same for newz.na

  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na

  col <- base::c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range
  ColorRamp <-col
  if(zlim[1]>=zlim[2] )zlim[2]<-zlim[1]+0.01*zlim[1]
  ColorLevels <- seq(zlim[1], zlim[2], length=100)

  graphics::par(mar =  base::c(5,4,4,2))
  graphics::par(mar =  base::c(6,yma,4,2))
  graphics::par(mar =  base::c(xma,yma,topyma,2))

  graphics::image(1:base::length(cLabels),1:base::length(rLabels), base::t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=zlim)
  graphics::mtext(row_unit,side=2,line=yline, cex=cex.lab)
  if( !base::is.null(main) ){
    graphics::title(main=main,cex.main=cex.main)
  }
  graphics::mtext(col_unit,side=1,line=xline, cex=cex.lab)
  graphics::axis(BELOW<-1, at=attc, labels=clbls, cex.axis=cex.axis,las= lasval, xlab=col_unit)
  graphics::axis(LEFT <-2, at=attr, labels=rlbls, las= HORIZONTAL<-1,cex.axis=cex.axis)
  graphics::par(mar = base::c(xma/5,yma/5,1,4))
  if(xma!=6){
    graphics::image(1, ColorLevels,
          base::matrix(data=ColorLevels[1:base::length(ColorLevels)-1], ncol=base::length(ColorLevels)-1),
          col=ColorRamp,cex.axis=cex.axis,
          xlab="",ylab="",
          xaxt="n") #cex.lab=cex.lab,ylab=zunit,
  } else {
    graphics::image(1, ColorLevels,
          base::matrix(data=ColorLevels[1:base::length(ColorLevels)-1], ncol=base::length(ColorLevels)-1),
          col=ColorRamp,cex.axis=cex.axis,
          xlab="",cex.lab=cex.lab,ylab=zunit,
          xaxt="n") #
  }

  graphics::par(.pardefault)
  graphics::layout(1)
  graphics::par(mfrow=base::c(1,1))

}

# plot with errors (either std dev or std error) original(orig) and/or
# signal, Equitable transform (zw) and least squares transform (lsx)
plotdata_with_errors<- function( dataset,data_std,rnames=1:base::length(dataset),cex=1,
                                 main="data",ylim=base::c(0,1),xlim=NULL,xlab="ROW",
                                 ylab="DATA VALUE",
                                 pch="O",type="p",col="black",lty=1,
                                 lineonly=FALSE,cex.main=1,
                                 x=NULL){

  #dataset and data_std are the data vector and error bars respectively
  numrows <- base::length(dataset)
  if(base::is.null(xlim))xlim<-base::c(1,numrows)
  if(base::is.null(x)){
    d = base::data.frame(
      x  = base::c(1:numrows)
      , y  = dataset
      , xsd = data_std
    )
  } else {
    d = base::data.frame(
      x  = x
      , y  = dataset
      , xsd = data_std
    )
    xlim<-base::c(base::min(x,na.rm=TRUE),base::max(x,na.rm=TRUE))
  }
  if(base::is.null(ylab)) ylab<-"Data Value"
  if(base::is.null(xlab)) xlab<-"Index"
  if(!lineonly){
    base::plot(d$x, d$y ,pch=pch, ylim= ylim,xlim= xlim,xlab=xlab, ylab=ylab,xaxt='n',yaxt='n',cex=cex,
         cex.lab=1.5, cex.axis=1.5,  cex.sub=1.5)
    base::with(
      data = d
      , expr = Hmisc::errbar(x, y, y+xsd, y-xsd, add=TRUE, pch=pch,type=type,col=col,lty=lty, cap=.01)
    )
    graphics::title(main=main,cex.main=cex.main)

    rtick<-xlim[1]
    er<-xlim[2]
    rtickinc<-((er-rtick)/10)

    if(xlim[1]==1 && xlim[2]==numrows) {
      rtick<-xlim[1]
      er<-xlim[2]
      rtickinc<-base::round((er-rtick)/10)
    }
    if(rtickinc==0)rtickinc=1
    graphics::axis(1,at=base::round(base::c(seq(rtick,er,by= rtickinc)),digits=4),labels=rnames[seq(rtick,er,by=rtickinc)],lwd=2,cex.lab=1.5, cex.axis=1.5,  cex.sub=1.5)

    graphics::axis(2,ylim= ylim,lwd=2,cex.lab=1.5, cex.axis=1.5,  cex.sub=1.5)
    graphics::box(lwd=2)
  } else {
    stats::line(d$x, d$y )
    base::with(
      data = d
      , expr = Hmisc::errbar(x, y, y+xsd, y-xsd, add=TRUE, pch=pch,type=type,col=col,lty=lty, cap=.01)
    )
  }

  base::return()
}

plot.frequency.spectrum <- function(fr=NULL, X.k, xlimits=NULL,main="",plotflag=TRUE) {
  if(base::is.null(fr)){
    plot.data  <- base::cbind(0:(base::length(X.k)-1), base::Mod(X.k))
    if(base::is.null(xlimits)) xlimits<-base::c(0,base::length(X.k))
  } else {
    plot.data  <- base::cbind(fr, base::Mod(X.k))
    if(base::is.null(xlimits))xlimits<-base::c(0,fr[base::length(fr)])
  }

  # TODO: why this scaling is necessary?
  plot.data[2:base::length(X.k),2] <- 2*plot.data[2:base::length(X.k),2]

  if(plotflag) base::plot(plot.data,  lwd=2, main=main, type="b",
                    xlab="Frequency (1/bp or Hz)", ylab="Strength",
                    xlim=xlimits, ylim=base::c(0,base::max(base::Mod(plot.data[,2]))))

  base::return(plot.data)
}

plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
  Xk.h <-  base::rbind(0,base::length(Xk))
  Xk.h[i+1] <- Xk[i+1] # i-th harmonic
  harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
  graphics::points(ts, harmonic.trajectory, type="l", col=color)
}

f <- function(t,w) {
  dc.component +
    base::sum( component.strength * base::sin(component.freqs*w*t + component.delay))
}

plot.fourier <- function(fourier.series, f.0, ts) {
  w <- 2*pi*f.0
  trajectory <- base::sapply(ts, function(t) fourier.series(t,w))
  base::plot(ts, trajectory, type="l", xlab="time", ylab="f(t)"); graphics::abline(h=0,lty=3)
  base::return(base::cbind(ts,trajectory))
}

#' dna_1_to_wax
#'
#' Makes DNA Walks
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
dna_1_to_wax<-function(dna_1,linerng=NULL){ #fastq<-TRUE   dna_1[1:10,1]   dn[1:10]
  mindna<-10
  if(base::is.null(linerng))linerng<-base::c(2,base::nrow(dna_1))
  waxy1<-dna_1[linerng[1]:linerng[2],1]   #was dn
  base::cat("\nlines in waxy1 ", base::length(waxy1),"\n")

  waxy<- base::paste0(waxy1,collapse="\n") #nchar(waxy)
  waxy1<-NULL;base::gc()
  base::cat("\nlength is ",nchar(waxy),"\n")

  if(mindna<nchar(waxy)){
    waxy<-base::casefold(waxy, upper=F)
    waxy<-stringr::str_replace_all(waxy, "([\n])", "")
    base::cat("\neliminate newline char : length is ",nchar(waxy),"\n")
    waxy<-(base::unlist(base::strsplit(waxy,split="")))
    base::cat("\nlength is ",base::length(waxy),"\n")
  }
  base::return(waxy)
}

tdna_to_wax<-function(tdna,main=""){
  mindna<-10
  waxy1<-(base::unlist(base::strsplit(tdna,split=" ")))
  waxy1<-(base::unlist(base::strsplit(waxy1,split="\n")))
  waxy1<-waxy1[(base::is.na(base::as.numeric(base::gsub("([0-9]+).*$", "\\1", waxy1))))]
  waxy<- base::paste0(waxy1,collapse="\n") #nchar(waxy)  #changed may 3

  base::cat("\nlength is ",nchar(waxy),"\n")
  if(mindna<nchar(waxy)){
    wax<-base::casefold(waxy, upper=F)
    wax<-stringr::str_replace_all(waxy, "([\n])", "")
    wax<-(base::unlist(base::strsplit(wax,split="")))
    graphics::barplot( base::prop.table(base::table(wax)),main= base::paste(main,"\nSense"))
  }
  base::cat("\nlength without newline char is ",nchar(waxy),"\n")
  base::return(wax)
}

fracD<-function(CG,AT,pflag=TRUE,main="",plotflag=TRUE ){
  dat <- base::data.frame(CG = CG,
                    AT = AT )
  uniquedat<-base::unique(dat) # base::plot(uniquedat)    base::length(dat$x)  base::plot(dat)
  fCGlim<-base::c(base::min(uniquedat$CG,na.rm = TRUE ), base::max(uniquedat$CG,na.rm = TRUE ))
  fATlim<-base::c(base::min(uniquedat$AT,na.rm = TRUE ), base::max(uniquedat$AT,na.rm = TRUE ))
  dCG<-fCGlim[2]-fCGlim[1]
  dAT<-fATlim[2]-fATlim[1]
  rangeboth<- base::max(dCG,dAT)+1
  total<-rangeboth*rangeboth
  fracdim<- 2*base::log(base::length(uniquedat$CG))/base::log(total)

  fracdim_rectangle<- 2*base::log(base::length(uniquedat$CG))/base::log((dCG+1)*(dAT+1))  #one of these 0 gives dimension 2
  m1<- base::paste("\nbase::length(CG)",base::length(CG),"|unique|",base::length(uniquedat$CG),
            "ranges(x,y)",dCG,dAT)

  if(pflag){
    xlim<-base::c(fCGlim[1],fCGlim[1]+rangeboth-1);ylim<-base::c(fATlim[1],fATlim[1]+rangeboth-1)
    base::plot(dat,main= base::paste(main,m1,"D",base::round(fracdim,digits=2) ),
         type="l",xlim=xlim,ylim=ylim,cex.main=0.6,las=2)
    graphics::points(dat$CG[1],dat$AT[1],type="p",pch=83,col="red",cex=2)   #S
    graphics::points(dat$CG[base::length(dat$CG)],dat$AT[base::length(dat$AT)],type="p",pch=69,col="red",cex=2)
  }
  if(plotflag){
    xlim<-base::c(fCGlim[1],fCGlim[1]+dCG);ylim<-base::c(fATlim[1],fATlim[1]+dAT)
    base::plot(dat,main= base::paste(main,m1,"D(rect)",base::round(fracdim_rectangle,digits=2) ),
         type="l",xlim=xlim,ylim=ylim,cex.main=0.6,las=2)
    graphics::points(dat$CG[1],dat$AT[1],type="p",pch=83,col="red",cex=2)
    graphics::points(dat$CG[base::length(dat$CG)],dat$AT[base::length(dat$AT)],type="p",pch=69,col="red",cex=2)
  }
  base::return(fracdim_rectangle)
}

seq_to_dnawalk<-function(x){

  ATval<- base::rep(0,base::length(base::c(x)))
  CGval<- base::rep(0,base::length(base::c(x)))
  MAGval<- base::rep(0,base::length(base::c(x)))

  Cindices<-base::grep("c",x=x)
  Gindices<-base::grep("g",x=x)
  Tindices<-base::grep("t",x=x)
  Aindices<-base::grep("a",x=x)

  ATval[Cindices]<-  0
  ATval[Gindices]<-  0
  ATval[Aindices]<-  1;  MAGval[Aindices]<-  1
  ATval[Tindices]<- (-1);  MAGval[Tindices]<- -1

  CGval[Cindices]<-  1;  MAGval[Cindices]<-  1
  CGval[Gindices]<- (-1);  MAGval[Gindices]<- -1
  CGval[Aindices]<-  0
  CGval[Tindices]<-  0

  Walklist<-base::list(atwalk=ATval,cgwalk=CGval,dnawalk=MAGval)
  base::return(Walklist)
}

#' walk_and_plot
#'
#' Plots DNA Walks
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
walk_and_plot<-function(tdna,main="",ptset=NULL,rng=NULL,seqlen=NULL,
                        sample_every=1,fullname=NULL,All_walk=NULL,
                        walkflag=TRUE,fracflag=FALSE, AT_flag=TRUE,
                        atflag=TRUE, spar=0.6,startval=1,endval=NULL,
                        printlong=FALSE,waxflag=FALSE,
                        freqerrorflag=TRUE,pflag=TRUE,plotflag=TRUE,
                        writeflag=TRUE){

  presampleflag<-FALSE
  mindna<-10
  if(base::is.null(All_walk)){   #All_walk

    if(!waxflag) wax<-tdna_to_wax(tdna) else wax<-tdna #

    if(mindna<base::length(wax)|| (!base::is.null(seqlen) && (seqlen <= base::length(wax))) ){
      if(!base::is.null(rng)){
        wax<-base::matrix(wax[rng[1]:rng[2]],nrow=1,ncol=base::length(rng[1]:rng[2]))
        if(!base::is.null(ptset)){
          main<- base::paste(main,"Rng (",rng[1],rng[2],")",startval+rng[1]-1,"_",startval+rng[2]-1)
        }
      } else {
        wax<-base::matrix(wax,nrow=1,ncol=base::length(wax)) # wax[,rng[1]:rng[2]]
      }
      # convert to cumulative walk
      Walklist<-seq_to_dnawalk(wax)     #nucleotide to numeric single values A+ve T-ve C+ve G-ve
      atwalk<-base::cumsum(Walklist$atwalk)   #cumulative sum of AT single walk values
      cgwalk<-base::cumsum(Walklist$cgwalk)    #cumulative sum of CG single walk values
      dnawalk<-base::cumsum(Walklist$dnawalk)  #cumulative sum of Net single walk values
    }
  } else{   # case where dna walk is previously found   startval and endval used

    #use  startval and rng to find values in All_walk
    starting<-base::as.numeric(base::rownames(All_walk)[1])
    if(base::is.null(endval))endval<-base::as.numeric(base::rownames(All_walk)[base::nrow(All_walk)])
    if(base::is.null(startval))startval<-base::as.numeric(base::rownames(All_walk)[1])
    if(sample_every<0){
      sval<-(startval-starting)/(-1*sample_every)+1
      eval<-(endval-starting)/(-1*sample_every)+1
    } else{
      sval<-startval-starting+1
      eval<-endval-starting+1
    }
    atwalk<-wax<-All_walk[sval:eval,"AT"]
    cgwalk<-All_walk[sval:eval,"CG"]

    dnawalk<-base::rep.int(0, base::length(sval:eval))   #base::rep.int(0,10)
    main<- base::paste(main,"\nsval eval (",sval,eval,") startval endval",startval,"_",endval)
  }

  if(mindna<base::length(wax)|| (!base::is.null(seqlen) && (seqlen <= base::length(wax))) ){
    len<-base::length(atwalk)
    if(!AT_flag){
      foo<-atwalk
      atwalk<-cgwalk
      cgwalk<-foo    #this reversal could cause mixup when plotting cg-AT plots
      foo<-NULL
      main<- base::paste("CG Walk Spar=",spar,"",main)
      walk<-base::cbind(atwalk,cgwalk)
      base::rownames(walk)<-startval+(0:(len-1));base::colnames(walk)<-base::c("CG","AT")
    } else {
      walk<-base::cbind(atwalk,cgwalk)
      base::rownames(walk)<-startval+(0:(len-1));base::colnames(walk)<-base::c("AT","CG")
    }

    if(!base::is.null(fullname)&& sample_every==1){
      ofile<- base::paste0(fullname,"_Table.txt")
    } else {
      ofile<-NULL
    }

    if(sample_every!=1){
      if(sample_every>1){
        whichval<-seq (1,len,by=sample_every)  #seq (1,20,by=2)   #whichval<-seq (1,140384,by=11)
      } else {
        presampleflag<-TRUE
        sample_every=sample_every*(-1);
        whichval<-seq (1,len,by=1);
      }

      dnawalk<-dnawalk[whichval]
      cgwalk<-cgwalk[whichval]
      atwalk<-atwalk[whichval]

      len<-base::length(atwalk)
      base::cat("\nnew length is",len,"\n")

      walk<-base::cbind(atwalk,cgwalk)
      if(atflag) base::colnames(walk)<-base::c("AT","CG") else base::colnames(walk)<-base::c("CG","AT")
      base::rownames(walk)<-startval+ base::seq(1,(len*sample_every),by=sample_every)-1    #base::rownames(dnawalkimaginaryat[1:len])
    }

    acq.freq <- 1/sample_every                    # data acquisition (sample) frequency (Hz)
    time     <- base::length(atwalk)* sample_every                    # measuring time interval (seconds)
    tim       <- base::seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s)

    if(base::is.null(seqlen)){
      seqlen<-base::length(atwalk)/2
      if(base::length(atwalk)> 1200)seqlen<-400
      if(base::length(atwalk)> 4000) seqlen<-base::round(base::length(atwalk)/5); #seqlen<-3000   4000/5
      if(base::length(atwalk)> 10000) seqlen<-base::round(base::length(atwalk)/5);  # 10000/10
      if(base::length(atwalk)> 20000) seqlen<-base::round(base::length(atwalk)/10);  # 20000/10
      if(base::length(atwalk)> 40000) seqlen<-base::round(base::length(atwalk)/20); # 40000/20
      if(base::length(atwalk)> 70000) seqlen<-base::round(base::length(atwalk)/30); # 70000/30
      if(base::length(atwalk)> 100000) seqlen<-base::round(base::length(atwalk)/40); # 100000/30
    }

    spectrumimage<-NULL;peaks<-NULL
    fracdim_list<-NULL
    sslist<- base::seq(1, (len-seqlen/sample_every+1),by= (seqlen/sample_every))
    # }
    for ( j in sslist){
      nextseq<-  atwalk[j:(j+seqlen/sample_every-1)];
      nextseqcg<-cgwalk[j:(j+seqlen/sample_every-1)]
      tim1<-tim[j:(j+seqlen/sample_every-1)]

      m1<- base::paste(main,"\n",sample_every*j,"to",sample_every*(j+seqlen/sample_every),":  ",
                sample_every*j+startval,sample_every*(j+seqlen/sample_every)+startval)

      fracDim<-fracD(CG=nextseqcg,AT=nextseq,main=m1,pflag=FALSE,plotflag=plotflag )

      x1<-base::c(base::min(tim1,na.rm=TRUE),base::max(tim1,na.rm=TRUE))
      nti<-20

      if(plotflag)base::plot(tim1,nextseq,type="l",xlab="",cex.main=0.75,
                       main= base::paste(m1,"D",base::round(fracDim,digits=2)),las=2,xaxp = base::c(x1[1],x1[2],nti))

      lowpass1.spline <- stats::smooth.spline(tim1,nextseq, spar = spar) ## Control spar for amount of smoothing
      if(plotflag)graphics::lines(stats::predict(lowpass1.spline, tim1), col = "red", lwd = 1)
      highpass1 <- nextseq - stats::predict(lowpass1.spline, tim1)$y
      if(plotflag)graphics::lines(tim1, highpass1, lwd =  2)
      if(plotflag)base::plot(tim1, highpass1,type="l",pch=15, lwd =  1,xlab="",cex.main=0.75,
                       main= base::paste(main,"\n",sample_every*j,"to",sample_every*(j+seqlen/sample_every),":  ",
                                  sample_every*j+startval,sample_every*(j+seqlen/sample_every)+startval,
                                  "D",base::round(fracDim,digits=2)),las=2,xaxp = base::c(x1[1],x1[2],nti) )


      # base::plot(nextseq,type="b",main= base::paste(j,"to",(j+seqlen)))
      samplefreq<-1/( seqlen/sample_every)
      fr<-(samplefreq*(1:(seqlen/sample_every))  )/sample_every   # seqlen<-1e5; sample_every<-4000
      if(writeflag){
        base::cat("\nsamplefreq ",samplefreq)
        if(freqerrorflag & j==1){
          upperfreqval<-(samplefreq*(1.5:(seqlen/sample_every+0.5))  )/sample_every
          lowerfreqval<-(samplefreq*(0.5:(seqlen/sample_every-0.5))  )/sample_every
          base::cat("\nlower repeat upper lowerlimit fequency upperlimit\n")
          for(f in 1:(base::length(upperfreqval)/2)){
            base::cat( base::round(1/upperfreqval[f] ,digits=4),
                 base::round(1/fr[f],digits=4),base::round(1/lowerfreqval[f],digits=4),
                 upperfreqval[f],fr[f],lowerfreqval[f],"\n")
          }
        }
      }
      X.k <- stats::fft(highpass1)                   # find all harmonics with stats::fft()

      plotdata1<-plot.frequency.spectrum(fr,X.k, xlimits=base::c(0,base::max(fr)/2),
                                         main= base::paste(sample_every*j,"to",sample_every*(j+seqlen/sample_every),":  ",
                                                    sample_every*j+startval,sample_every*(j+seqlen/sample_every)+startval),
                                         plotflag=plotflag)    # base::c(0,20)

      windowpeak<-pracma::findpeaks(plotdata1[,2])
      if(!base::is.null(windowpeak))

        if(writeflag){
          if(!base::is.null(windowpeak))      base::cat("\npeaks for window ",sample_every*j,sample_every*(j+seqlen/sample_every-1),
                                            "     ",startval+sample_every*(j-1),startval+sample_every*(j+seqlen/sample_every-1),"\n")
          windowpeak<-base::cbind(fr[windowpeak[,2]],base::round(1/fr,digits=4)[windowpeak[,2]],windowpeak)

          base::colnames(windowpeak)<-base::c("freq","repeat length","height", "index", "start", "end")
          if(!base::is.null(windowpeak)){
            if(!base::is.null(windowpeak)) if(base::nrow(windowpeak)>=200) base::print(windowpeak[1:200,]) else if(base::nrow(windowpeak)>=100) base::print(windowpeak[1:100,]) else if(base::nrow(windowpeak)>=50) base::print(windowpeak[1:50,]) else if(base::nrow(windowpeak)>=25) base::print(windowpeak[1:25,]) else if(base::nrow(windowpeak)>=10) base::print(windowpeak[1:10,])
          }
        } else{
          windowpeak<-base::matrix(base::c(0,0,0, 0, 0, 0),ncol=6,nrow=1)
          base::colnames(windowpeak)<-base::c("freq","repeat length","height", "index", "start", "end")
        }

      spectrumimage<-base::cbind(spectrumimage,plotdata1[,2])
      fracdim_list<-base::c(fracdim_list,fracDim)

    }

    base::rownames(spectrumimage)<-  base::paste0("1/",base::round(1/fr,digits=4))
    base::colnames(spectrumimage)<-startval+sample_every*
      base::seq(seqlen/(2*sample_every), (len-seqlen/(2*sample_every)+1),by= seqlen/sample_every)# base::rownames(spectrumimage)<-ts
    if( base::ncol(spectrumimage)<2){
      spectrumimage<-base::cbind(spectrumimage,spectrumimage)
      fracdim_list<-base::c(fracdim_list,fracdim_list)
    }

    base::names(fracdim_list)<-base::colnames(spectrumimage)

    minf<-base::as.numeric(names(fracdim_list))[1]
    maxf<-base::as.numeric(names(fracdim_list))[base::length(fracdim_list)]

    numrange<- base::round((maxf-minf)/(15*seqlen))
    if(numrange<1)numrange<-1
    if(plotflag)run_sum_Fractal(fracdim_list=fracdim_list, numrange=numrange ,pflag=pflag,
                                main= base::paste("",main))

    x<-plotdata1[,1];
    freqlim<-1:(base::nrow(spectrumimage)/2);

    rmean<-base::rowMeans(spectrumimage[freqlim,],na.rm=TRUE)   #this is average over range
    allrmean<-rmean; base::names(allrmean)<-fr[freqlim]
    spectrumimage1<-spectrumimage[freqlim,]
    rstd<-matrixStats::rowSds(spectrumimage[freqlim,],na.rm=TRUE)  #set ylim?
    stderr<-rstd/base::sqrt( base::ncol(spectrumimage))

    peaks<-pracma::findpeaks(rmean)
    if(!base::is.null(peaks)){
      peaks<-base::cbind(fr[peaks[,2]],base::round(1/fr,digits=4)[peaks[,2]],rstd[peaks[,2]],stderr[peaks[,2]],peaks)
      base::colnames(peaks)<-base::c("freq","repeat length","Std dev","std err","height", "index", "start", "end")

    } else{
      peaks<-base::matrix(base::c(0,0,0, 0, 0, 0,0,0),ncol=8,nrow=1)
      base::colnames(peaks)<-base::c("freq","repeat length","Std dev","std err","height", "index", "start", "end")
    }


    if(plotflag){
      if(base::nrow(spectrumimage)>200){freqlim<-1:200;
      rmean<-base::rowMeans(spectrumimage[freqlim,],na.rm=TRUE)

      rstd<-matrixStats::rowSds(spectrumimage[freqlim,],na.rm=TRUE)  #set ylim?
      stderr<-rstd/base::sqrt( base::ncol(spectrumimage))
      plotdata_with_errors( rmean,x=fr[freqlim],stderr, main= base::paste("Average and std err over all Windows\n",main),
                            ylim=NULL,xlim=NULL,xlab="1/repeatlength",ylab="Power",
                            pch="O",type="b",col="black",lty=1,lineonly=FALSE,cex.main=1)
      }

      if(base::nrow(spectrumimage)>80){freqlim<-1:80;
      rmean<-base::rowMeans(spectrumimage[freqlim,],na.rm=TRUE)
      rstd<-matrixStats::rowSds(spectrumimage[freqlim,],na.rm=TRUE)  #set ylim?
      stderr<-rstd/base::sqrt( base::ncol(spectrumimage))
      plotdata_with_errors( rmean,x=fr[freqlim],stderr, main= base::paste("Average and std err over all Windows\n",main),
                            ylim=NULL,xlim=NULL,xlab="1/repeatlength",ylab="Power",
                            pch="O",type="b",col="black",lty=1,lineonly=FALSE,cex.main=1)
      }
    }

    if(base::is.null(ptset)){
      #set the points to subdivide the fractal graph two sets one for the legend (ptset1) and one for the graph(ptset)
      numseg<-10;
      ptset1<-startval+ sample_every*base::seq((base::length(cgwalk )/numseg+1),base::length(cgwalk ),
                                         by=(base::length(cgwalk )/numseg))
      ptset<- base::seq((base::length(cgwalk)/numseg+1),base::length(cgwalk ),
                  by=(base::length(cgwalk )/numseg))

    } else{
      numseg<-base::length(ptset)
      if(ptset[1]>=startval){
        ptset1<-ptset
        ptset<-ptset-startval+1
        if(sample_every>1)ptset<-(ptset-startval)/sample_every+1
      } else {  # ptset are low values so increase ptset1
        ptset1<-startval+ptset
        if(sample_every>1) ptset<-base::ceiling((ptset/sample_every)) #base::ceiling(1)
      }
    }

    if(plotflag) {
      col<-grDevices::rainbow(base::length(ptset)) #grDevices::rainbow(10)
      base::plot(sample_every*(1:(base::length(dnawalk))),dnawalk[1:(base::length(dnawalk))],type="l",cex.main=0.75,
           main= base::paste("NET DNAWALK",main))
      base::plot(sample_every*(1:base::length(cgwalk)),cgwalk[1:base::length(cgwalk)],cex.main=0.75,
           type="l",main= base::paste("CG DNAWALK",main));
      base::plot(sample_every*(1:base::length(atwalk)),atwalk[1:base::length(atwalk)],
           type="l",cex.main=0.75,main= base::paste("AT DNAWALK",main))
      base::plot(cgwalk,atwalk,lty=1,lwd=1,
           type="l",cex.main=0.75,
           main= base::paste("CG-AT DNA WALK Sense only: increment ",
                      base::round((base::length(cgwalk)/numseg),digits=0),"\n", main));
      offsetlist<-base::c(base::round((base::max(cgwalk[1:(base::length(cgwalk))],na.rm=TRUE)-
                             base::min(cgwalk[1:base::length(cgwalk)],na.rm=TRUE))/5))
      if(!base::is.numeric(offsetlist[1])|| base::is.na(offsetlist[1])) {
        base::cat("\n reset offst to 0")
        offsetlist[1]<-0
      }

      if(base::length(ptset)!=0){
        ptset<-base::c(1, ptset,base::length(cgwalk))
        ptset1<-base::c(startval+1, ptset1,startval+sample_every*base::length(cgwalk))

      }

      for (offst in offsetlist ){  #offst<- 0 changes the plot scale on x to see the legend better
        xli<-base::c(base::min(cgwalk[1:base::length(cgwalk)]-offst,na.rm = TRUE),
               base::max(cgwalk[1:base::length(cgwalk)],na.rm = TRUE))

        base::plot(cgwalk[1:base::length(cgwalk)],atwalk[1:base::length(cgwalk)],lty=1,lwd=1,
             type="l",cex.main=0.75,
             xlim=xli,
             main= base::paste("Sense only: increment ",
                        base::round((sample_every*(base::length(cgwalk) )/numseg),digits=0),"\n", main));
        graphics::lines(xli, -xli); graphics::lines(xli, xli)

        if(base::length(ptset)!=0 && base::min(ptset,rm.na=TRUE)>0){
          for( jpt in 2:(base::length(ptset)-1)){
            lwd=1
            graphics::lines(cgwalk[base::round(ptset[jpt-1]):base::round(ptset[jpt])], atwalk[base::round(ptset[jpt-1]):base::round(ptset[jpt])],
                  lty=1,lwd=lwd,col=col[jpt-1])
            graphics::points(cgwalk[base::round(ptset[jpt])],atwalk[base::round(ptset[jpt])],
                   pch=base::as.character(jpt-1),cex=2, col="black")
          }
          for( jpt in 2:(base::length(ptset)-1)){
            graphics::points(cgwalk[base::round(ptset[jpt])],atwalk[base::round(ptset[jpt])],
                   pch=base::as.character(jpt-1),cex=2, col="black")
          }

          graphics::legend("topleft",legend=base::c(base::round(ptset1[2:(base::length(ptset1))],digits=0)  ),
                 pch=base::c(base::as.character(1:(base::length(ptset)-2)),"E"))
        }
        graphics::points(cgwalk[1],atwalk[1], pch=83,cex=2)
        graphics::points(cgwalk[(base::length(cgwalk))],atwalk[base::length(atwalk)], pch=69,cex=2)

      }

    }  #end of plotflag check

  } else {
    peaks<-base::matrix(base::c(0,0,0, 0, 0, 0,0,0),ncol=8,nrow=1);spectrumimage1<-NA; allrmean<-NA; fracdim_list<-NA
  }
  spectinfo<-base::list(peaks, spectrumimage1, allrmean, walk,ofile, fracdim_list)  #ofile is NULL if sample_every=1
  base::return(spectinfo)
} #end of walk and plot spar near has no effective smoothing

#' run_one_repeat
#'
#' Finds peaks in one repeat length
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_one_repeat<-function(repeat_val,All_spec,chromnam,numstd=5,colrange=NULL,sortstr=FALSE){
  if (!base::is.null(colrange)){
    All_spec<-All_spec[,((base::which(base::colnames(All_spec)==colrange[1]):base::which(base::colnames(All_spec)==colrange[2])))]
    base::print(base::colnames(All_spec))
  }
  All_spec1<-All_spec
  All_spec1[base::which(All_spec1==0.0)]<-NA
  maxspec<-base::max(All_spec1[repeat_val, ],na.rm=TRUE)
  meanspec<-base::mean(All_spec1[repeat_val, ],na.rm=TRUE)
  stdspec<-stats::sd(All_spec1[repeat_val, ],na.rm=TRUE)

  threshspec<-meanspec+numstd*stdspec    #base::length(peak_along_repeat)

  peak_along_repeat<-pracma::findpeaks(x=All_spec[repeat_val, ], minpeakheight = threshspec,
                               minpeakdistance = 1, threshold = 0, npeaks = 0, sortstr = sortstr)  #change to false Feb 5 2023
  if(base::length(peak_along_repeat)!=0){  #change June 17 2022
    peak_along_repeat<-base::cbind(base::as.numeric(peak_along_repeat[,1])/maxspec,peak_along_repeat)
    peak_along_repeat<-base::cbind(peak_along_repeat,base::as.numeric(base::colnames(All_spec)[peak_along_repeat[,3]]),
                             base::as.numeric(base::colnames(All_spec)[peak_along_repeat[,4]]),
                             base::as.numeric(base::colnames(All_spec)[peak_along_repeat[,5]]) )

    base::colnames(peak_along_repeat)<-base::c("Percent of peak","Power","peak index","sindex","eindex","peak_bp","start_peak_bp", "end_peak_bp")
    peak_along_repeat
  } #end of check to ensure peak value found
  if(!base::is.null(colrange) | base::length(peak_along_repeat)!=0){
    base::plot(base::as.numeric(base::colnames(All_spec1)),All_spec1[repeat_val ,],cex.main=0.75,
         main= base::paste("Chromosome",chromnam,"\n1/repeat length (bp)",repeat_val,"\nmean + ",numstd,
                    "*sigma= ",base::round(threshspec,digits=1)),
         ylim=base::c(0,maxspec),type="o",las=2,xlab="",ylab="",
         xaxp=base::c( base::min(base::as.numeric(base::colnames(All_spec1)),na.rm=TRUE),
                 base::max(base::as.numeric(base::colnames(All_spec1)),na.rm=TRUE),20) )       #best
  }
  if(base::length(peak_along_repeat)!=0){
    graphics::lines(base::as.numeric(base::colnames(All_spec1)[base::as.numeric(peak_along_repeat[,"peak index"])]),
          All_spec1[repeat_val ,base::as.numeric(peak_along_repeat[,"peak index"])],type="p",pch=15,col="red",cex=1.5)
    base::cat("track type=bedGraph")
  }
  base::return(peak_along_repeat)
} #end of function

CG_AT_content<-function(dn1,nam, winbox=5000){
  wax<-tdna_to_wax(dn1)
  wax<-base::matrix(wax,nrow=1,ncol=base::length(wax))  # wax[1,1:10]
  walklist<-seq_to_dnawalk(wax)     #nucleotide to numeric single values A+ve T-ve C+ve G-ve
  at<-walklist$atwalk
  cg<-walklist$cgwalk
  atwalk<-NULL;cgwalk<-NULL;CG_ATwalk<-NULL
  bpseq<-base::seq(1,base::length(at),by= winbox)
  for(s in base::seq(1,base::length(at),by= winbox)){
    e=s+winbox-1
    if(e>base::length(dn1))e<-base::length(dn1)
    vat<-base::sum(base::abs(at[s:e]))
    atwalk<-base::c(atwalk,vat)   #cumulative sum of AT single walk values
    vcg<-base::sum(base::abs(cg[s:e]) )
    cgwalk<-base::c(cgwalk,vcg )  #cumulative sum of CG single walk values
    if(vat>0){CG_ATrat<-vcg/vat } else { CG_ATrat<-NA}
    CG_ATwalk<-base::c(CG_ATwalk,CG_ATrat)  #cumulative sum of Net single walk values

    # atwalkcum<-base::cumsum(base::abs(Walklist$atwalk[s:e]))   #cumulative sum of AT single walk values
    # cgwalkcum<-base::cumsum(base::abs(Walklist$cgwalk[s:e]) )   #cumulative sum of CG single walk values
    # dnawalkcum<-base::cumsum(base::abs(Walklist$dnawalk[s:e]))  #cumulative sum of Net single walk values
  }

  peaks<-pracma::findpeaks(CG_ATwalk,sortstr=TRUE)   #height, index, start end of peak
  if(!base::is.null(peaks)){
    peaks<-base::cbind(bpseq[peaks[,2]],peaks)
    base::colnames(peaks)<-base::c("bpseq","CG_AT ratio", "index", "start", "end")

  } else{
    peaks<-base::matrix(base::c(0,0,0,0,0),ncol=5,nrow=1)
    base::colnames(peaks)<-base::c("bpseq","CG_AT ratio", "index", "start", "end")
  }
  base::print(peaks)
  #base::plot(walklist$cgwalk,walklist$atwalk,pch=15,type="o")
  base::plot(bpseq,atwalk,pch=15,type="o",main= base::paste(nam,winbox),las=2)
  base::plot(bpseq,cgwalk,pch=15,type="o",main= base::paste(nam,winbox),las=2)
  base::plot(bpseq,CG_ATwalk,pch=15,type="o",main= base::paste(nam,winbox,"\nmax=",peaks[1,1]),las=2)
  walklist<-base::list(cgwalk=cgwalk,atwalk=atwalk,CG_ATwalk=CG_ATwalk)
  base::return(walklist)
}

run_sum_Fractal<-function(fracdim_list, numrange= (80*1e6/(100000)),pflag=TRUE,main=""){
  if(base::length(fracdim_list)>20)nti<-20 else nti<-base::length(fracdim_list)
  x1<-base::c(base::min(base::as.numeric(names(fracdim_list)),na.rm=TRUE),base::max(base::as.numeric(names(fracdim_list)),na.rm=TRUE))
  nti<-20
  nstd<-3.5
  if(numrange<=2000)nstd<-3
  if(numrange<=300)nstd<-2
  if(numrange<=60)nstd<-1
  if(numrange<=30)nstd<-0.75
  meanD<-base::mean(fracdim_list,na.rm=TRUE)
  stdD<-stats::sd(fracdim_list,na.rm=TRUE)
  high_frac_thresh<-meanD+nstd*stdD
  low_frac_thresh<-meanD-nstd*stdD
  if(pflag){
    base::plot(base::as.numeric(names(fracdim_list)),fracdim_list,main= base::paste("Box Filling Dimension (Complexity) mean",
                                                                 base::round(meanD,digits=2),"sigma=",
                                                                 base::round(stdD,digits=2),"numrange",numrange,"\n",main),xlab="",
         ylim=base::c(0.8,2),type="b",pch=1,xaxp = base::c(x1[1],x1[2],nti),las=2,cex.main=0.6)



    highpts<-base::which(fracdim_list>high_frac_thresh)
    lowpts<-base::which(fracdim_list<low_frac_thresh)


    graphics::lines(base::as.numeric(names(fracdim_list)), base::rep(meanD,base::length(fracdim_list)),type="l",col="green",lwd=3)

    graphics::lines(base::as.numeric(names(fracdim_list[highpts])),fracdim_list[highpts],type="p",pch=16,col="red")
    graphics::lines(base::as.numeric(names(fracdim_list[lowpts])),fracdim_list[lowpts],type="p",pch=16,col="blue")

    graphics::legend("topleft",legend=base::c( base::paste0("+",nstd," Sigma"), base::paste0("-",nstd," Sigma")),pch=base::c(16,16),col=base::c("red","blue"))
    base::plot(base::as.numeric(names(fracdim_list)),fracdim_list,main= base::paste("Box Filling Dimension (Complexity)  mean",
                                                                 base::round(meanD,digits=2),"sigma=",
                                                                 base::round(stdD,digits=2),"numrange",numrange,"\n",main),xlab="",
         type="b",pch=1,xaxp = base::c(x1[1],x1[2],nti),las=2,cex.main=0.6)

    graphics::lines(base::as.numeric(names(fracdim_list)), base::rep(meanD,base::length(fracdim_list)),type="l",col="green",lwd=3)

    graphics::lines(base::as.numeric(names(fracdim_list[highpts])),fracdim_list[highpts],type="p",pch=16,col="red")
    graphics::lines(base::as.numeric(names(fracdim_list[lowpts])),fracdim_list[lowpts],type="p",pch=16,col="blue")
    graphics::legend("topleft",legend=base::c( base::paste0("+",nstd," Sigma"), base::paste0("-",nstd," Sigma")),pch=base::c(16,16),col=base::c("red","blue"))

  }
  base::cat("\n Box Filling Dimension (Complexity):        Max summary\n");base::print(base::summary(fracdim_list))

  meanfracD<-base::mean(fracdim_list,na.rm=TRUE)
  stdfracD<-stats::sd(fracdim_list,na.rm=TRUE)
  threshfracD<-meanfracD+nstd*stdfracD    #base::length(peak_along_repeat)    numstd<-3
  threshfracD_N<-meanfracD
  threshfracD<-0   #this is the threshold used

  bpval<-base::as.numeric(names(fracdim_list))  #

  deltabp<-(bpval[2]-bpval[1])
  deltabp_onebin<-deltabp*numrange
  base::cat("\n deltabp_onebin",deltabp_onebin)
  start_bpseqval<-base::seq(bpval[1],bpval[base::length(bpval)],by= deltabp_onebin)  #changed Sep 4 2022
  end_bpseqval<-start_bpseqval+deltabp_onebin
  end_bpseqval[base::length(end_bpseqval)]<-bpval[base::length(bpval)]
  which_ones<-base::which(fracdim_list >=threshfracD)
  which_ones_N_above<-base::which(fracdim_list >=threshfracD_N)  #for N count check values above and below mean
  which_ones_N_below<-base::which(fracdim_list <threshfracD_N)

  power_mean<-NULL
  power_sum<-power_N_above<-power_N_below<-NULL
  for(j in 1:base::length(start_bpseqval)){
    sbpval<-start_bpseqval[j]
    ebpval<-end_bpseqval[j]

    whichval<-base::which(bpval>=sbpval & bpval<=ebpval)
    sval<-whichval[1]
    eval<-whichval[base::length(whichval)]

    which_notNA<-base::which(!base::is.na(fracdim_list[sval: eval]))

    bplength_onebin<-base::length(which_notNA)*deltabp
    bp_notNA<-NULL
    for (k in 2: base::length(which_notNA)){ #change Nov 16
      bp_notNA<-base::c(bp_notNA,bpval[which_notNA[k]]-bpval[which_notNA[k-1]]) #  bp_notNA[7000:7200]
    }

    deltabpnotNA<-(bplength_onebin)/1.e6    # normalizes N and power_sum to per Mbp
    which_val<-base::which(which_ones<eval & which_ones>=sval)
    if(base::length(which_notNA)!=0){
      meanabove<-base::mean(fracdim_list[which_ones[which_val]] ,na.rm=TRUE)
      power_mean<-base::c(power_mean,meanabove)
      sumabove<-base::sum(fracdim_list[which_ones[which_val]] ,na.rm=TRUE)/deltabpnotNA
      power_sum<-base::c(power_sum,sumabove)
      which_val_above<-base::which(which_ones_N_above<eval & which_ones_N_above>=sval)
      Nval<-fracdim_list[which_ones_N_above[which_val_above]]
      Nabove<-base::length(base::which(!base::is.na(Nval))) /deltabpnotNA
      power_N_above<-base::c(power_N_above,Nabove)
      which_val_below<-base::which(which_ones_N_below<eval & which_ones_N_below>=sval)
      Nval<-fracdim_list[which_ones_N_below[which_val_below]]
      Nbelow<-base::length(base::which(!base::is.na(Nval))) /deltabpnotNA
      power_N_below<-base::c(power_N_below,Nbelow)
      base::names(power_mean)[base::length(power_mean)]<-sbpval+bplength_onebin/2
      base::names(power_sum)<-names(power_N_above)<-names(power_N_below)<-names(power_mean)
    } else{
      meanabove<-NA;power_mean<-base::c(power_mean,NA);power_sum<-base::c(power_sum,NA) ;
      Nbelow<-NA;Nabove<-NA;
      power_N_above<-base::c(power_N_above,NA);power_N_below<-base::c(power_N_below,NA)
      base::names(power_mean)[base::length(power_mean)]<-sbpval+bplength_onebin/2
      base::names(power_sum)<-names(power_N_above)<-names(power_N_below)<-names(power_mean)
      base::cat("\n no data : setting value to NA","\n")
    }
  }

  if( !base::is.null(power_mean)){

    maxspec<-base::max(power_mean,na.rm=TRUE)
    minspec<-base::min(power_mean,na.rm=TRUE)

    mincol<-base::which(minspec==power_mean)[1]
    maxcol<-base::which(maxspec==power_mean)[1]      #base::as.numeric
    min_mean_seqval<-base::as.numeric(names(power_mean[mincol]) )#+deltabp_onebin/2;
    max_mean_seqval<-base::as.numeric(names(power_mean[maxcol]) )#+deltabp_onebin/2
    if(pflag){

      graphics::barplot(power_mean,names.arg = base::names(power_mean),ylim=base::c(minspec,maxspec),
              main= base::paste("Box Filling Dimension (Complexity) sigma=",
                         base::round(stdD,digits=2),"numrange",numrange,"\n",
                         "Mean",
                         main,"\nMid-Points: min at:",min_mean_seqval,"max at:",max_mean_seqval),
              las=2,cex.axis = 0.65,cex.names=0.65,cex.main=0.6)
    }

    maxspec<-base::max(power_sum,na.rm=TRUE)
    minspec<-base::min(power_sum,na.rm=TRUE)

    if(base::is.na(minspec)|base::is.infinite(minspec)|base::is.null(minspec))minspec<-0
    if(base::is.na(maxspec)|base::is.infinite(maxspec)|base::is.null(maxspec))maxspec<-1
    if((maxspec-minspec)==0){ minspec<-0;maxspec<-1}
    mincol<-base::which(minspec==power_sum)[1]
    maxcol<-base::which(maxspec==power_sum)[1]      #base::as.numeric
    min_sum_seqval<-base::as.numeric(names(power_sum[mincol]) )#+deltabp_onebin/2;
    max_sum_seqval<-base::as.numeric(names(power_sum[maxcol]) )#+deltabp_onebin/2
    if(pflag){

      graphics::barplot(power_sum,names.arg = base::names(power_sum),ylim=base::c(minspec,maxspec),
              main= base::paste("Box Filling Dimension (Complexity) sigma=",
                         base::round(stdD,digits=2),"numrange",numrange,"\n",
                         "Sum/Mbp above",base::round(threshfracD,digits=2),
                         main,"\nMid-Points: min at:",min_sum_seqval,"max at:",max_sum_seqval),
              las=2,cex.axis = 0.65,cex.names=0.65,cex.main=0.6) #"#std above mean", nstd,"\n",
    }

    minspec<-maxspec<-NULL
    if(base::length(base::which(!base::is.na(power_N_below)))>0)maxspec<-base::max(power_N_below,na.rm=TRUE)  #x<-NULL; base::length(x); x<-base::c(NA,NA); base::length(base::which(!base::is.na(x)))
    if(base::length(base::which(!base::is.na(power_N_below)))>0)minspec<-base::min(power_N_below,na.rm=TRUE)

    if(!base::is.null(minspec)){

      if(base::is.na(minspec)|base::is.infinite(minspec)|base::is.null(minspec))minspec<-0
      if(base::is.na(maxspec)|base::is.infinite(maxspec)|base::is.null(maxspec))maxspec<-1
      if((maxspec-minspec)==0){ minspec<-0;maxspec<-1}
      mincol<-base::which(minspec==power_N_above)[1]
      maxcol<-base::which(maxspec==power_N_above)[1]      #base::as.numeric
      min_N_above_seqval<-base::as.numeric(names(power_N_above[mincol]) )#+deltabp_onebin/2;
      max_N_above_seqval<-base::as.numeric(names(power_N_above[maxcol]) )#+deltabp_onebin/2
      if(pflag){

        graphics::barplot(power_N_above,names.arg = base::names(power_N_above),ylim=base::c(minspec,maxspec),
                main= base::paste("Box Filling Dimension (Complexity) sigma=",
                           base::round(stdD,digits=2),"numrange",numrange,"\n",
                           "N/Mbp above mean",base::round(threshfracD_N,digits=2),
                           main,"\nMid-Points: min at:",min_N_above_seqval,"max at:",max_N_above_seqval),
                las=2,cex.axis = 0.65,cex.names=0.65,cex.main=0.6)
      }
    } else {
      min_N_above_seqval<-max_N_above_seqval<-NA  #
    }

    minspec<-maxspec<-NULL
    if(base::length(base::which(!base::is.na(power_N_below)))>0)maxspec<-base::max(power_N_below,na.rm=TRUE)  #x<-NULL; base::length(x); x<-base::c(NA,NA); base::length(base::which(!base::is.na(x)))
    if(base::length(base::which(!base::is.na(power_N_below)))>0)minspec<-base::min(power_N_below,na.rm=TRUE)

    if(!base::is.null(minspec)){
      if(base::is.na(minspec)|base::is.infinite(minspec))minspec<-0
      if(base::is.na(maxspec)|base::is.infinite(maxspec))maxspec<-1
      if((maxspec-minspec)==0){ minspec<-0;maxspec<-1}


      mincol<-base::which(minspec==power_N_below)[1]
      maxcol<-base::which(maxspec==power_N_below)[1]      #base::as.numeric
      min_N_below_seqval<-base::as.numeric(names(power_N_below[mincol]) )#+deltabp_onebin/2;
      max_N_below_seqval<-base::as.numeric(names(power_N_below[maxcol]) )#+deltabp_onebin/2
      if(pflag){

        graphics::barplot(power_N_below,names.arg = base::names(power_N_below),ylim=base::c(minspec,maxspec),
                main= base::paste("Box Filling Dimension (Complexity) sigma=",
                           base::round(stdD,digits=2),"numrange",numrange,"\n",
                           "N/Mbp below mean",base::round(threshfracD_N,digits=2),
                           main,"\nMid-Points: min at:",min_N_below_seqval,"max at:",max_N_below_seqval),
                las=2,cex.axis = 0.65,cex.names=0.65,cex.main=0.6)
      }
    } else {
      min_N_below_seqval <-max_N_below_seqval<-NA  #min_N_above_seqval<-max_N_above_seqval<-
    }

  } else {
    min_mean_seqval<-max_mean_seqval<-min_sum_seqval<-max_sum_seqval<- NA

    min_N_above_seqval<-max_N_above_seqval<- min_N_below_seqval <-max_N_below_seqval<-NA

  }

  pow_list<-base::list(power_mean=power_mean, min_mean_seqval=min_mean_seqval,max_mean_seqval=max_mean_seqval,
                 power_sum=power_sum, min_sum_seqval=min_sum_seqval,max_sum_seqval=max_sum_seqval,
                 power_N_above=power_N_above, min_N_above_seqval=min_N_above_seqval,max_N_above_seqval=max_N_above_seqval,
                 power_N_below=power_N_below, min_N_below_seqval=min_N_below_seqval,max_N_below_seqval=max_N_below_seqval
  )

  base::return(pow_list)

}

Ave_spectra<-function(nam=nam,All_spec_long5Mbp=All_spec_long5Mbp,fftlen=fftlen,dseq=(15000),numbp=1,
                      seqlim=NULL, repeatlim=NULL ){
  if(base::is.null(seqlim)){
    seqlim<-base::c(base::as.numeric(base::colnames(All_spec_long5Mbp)[1]),
              base::as.numeric(base::colnames(All_spec_long5Mbp)[ base::ncol(All_spec_long5Mbp)]))
  }
  Ave_spectra<-NULL
  for(j in base::seq(seqlim[1],(seqlim[2]-(dseq/2)),by=dseq)){
    if((j+dseq)<=seqlim[2]) e1<-j+dseq else  e1<-seqlim[2]
    bplim<-base::c(j,e1)
    spec_profile<- Average_spectral_profile(nam=nam,All_spec_long5Mbp=All_spec_long5Mbp,numbp=numbp,
                                            fftlen=fftlen,repeatlim=repeatlim,bplim=bplim,pspectra=FALSE)
    Ave_spectra<-base::cbind(Ave_spectra,spec_profile)
    if(!base::is.null(Ave_spectra)) base::colnames(Ave_spectra)[ base::ncol(Ave_spectra)]<-j+base::round(dseq/2)
  }
  base::return(Ave_spectra)
}

Average_spectral_profile<-function(nam=NULL,All_spec_long5Mbp=NULL,numbp=1,
                                   fftlen=5000,repeatlim=NULL,bplim=NULL,pspectra=TRUE){

  if((bplim[2]-bplim[1])>=(2*fftlen)){
    if(base::is.null(bplim)){
      bplim<-base::c(base::as.numeric(base::colnames(All_spec_long5Mbp)[1]),
               base::as.numeric(base::colnames(All_spec_long5Mbp)[ base::ncol(All_spec_long5Mbp)]))
    }


    main<- base::paste(nam,fftlen,"\nAve over",numbp,"repeats\n",repeatlim[1],repeatlim[2],bplim[1],bplim[2])
    repeats<-base::unlist(base::strsplit(base::rownames(All_spec_long5Mbp),"/"))
    if(base::is.null(repeatlim)){
      repeatlim<-base::c(repeats[1],
                   repeats[base::length(repeats)])
    }
    repeats<-base::as.numeric(repeats[base::seq(2,base::length(repeats), by=2)])
    whi1<-  base::which(repeats>=repeatlim[1])
    whi1<-whi1[base::length(whi1)]
    whi2<-  base::which(repeats<=repeatlim[2])
    whi2<-whi2[1]

    repeatindex<-base::c(whi2, whi1)
    freqlim<-(whi2:whi1)   #sets repeat range

    whi2<-base::which(base::as.numeric(base::colnames(All_spec_long5Mbp))<=bplim[2])
    bpindex<-  base::which(base::as.numeric(base::colnames(All_spec_long5Mbp))>=bplim[1])[1]:whi2[base::length(whi2)]
    repeats<- repeats[freqlim]


    rmean<-base::rowMeans(All_spec_long5Mbp[freqlim,bpindex],na.rm=TRUE)

    base::names(rmean)<-base::rownames(All_spec_long5Mbp[freqlim,])  # has average repeat profile over sequence range
    rmean_repeat_ave<-NULL; repeats_ave<-NULL
    if(numbp>1){
      for (j in base::seq(1,(base::length(rmean)-numbp),by = numbp)){
        rmean_repeat_ave<-base::c(rmean_repeat_ave, base::mean(rmean[j:(j+numbp-1)], rm.na=TRUE))
        base::names(rmean_repeat_ave)[base::length(rmean_repeat_ave)]<-names(rmean)[j+((numbp-1)/2)]
        repeats_ave<-base::c(repeats_ave, repeats[j+((numbp-1)/2)])
      }
    } else { rmean_repeat_ave<-rmean; repeats_ave<-repeats}
    if(pspectra){
      powval<-1.25
      for (r in 1:base::length(rmean_repeat_ave)){
        rmean_repeat_ave[r]<-rmean_repeat_ave[r]/((repeats_ave[r]-2+2 )*1)^powval
      }
      windowpeak<-pracma::findpeaks(rmean_repeat_ave)
      windowpeak<-base::cbind(repeats_ave[windowpeak[,2]],base::round(1/repeats_ave,digits=7)[windowpeak[,2]],windowpeak)

      base::colnames(windowpeak)<-base::c("repeat length","freq","height", "index", "start", "end")
      if(!base::is.null(windowpeak)&& (pspectra)){
        if(base::nrow(windowpeak)>=200) base::print(windowpeak[1:200,]) else if(base::nrow(windowpeak)>=100) base::print(windowpeak[1:100,]) else if(base::nrow(windowpeak)>=50) base::print(windowpeak[1:50,]) else if(base::nrow(windowpeak)>=25) base::print(windowpeak[1:25,]) else if(base::nrow(windowpeak)>=10) base::print(windowpeak[1:10,])
      } else{
        windowpeak<-base::matrix(base::c(0,0,0, 0, 0, 0),ncol=6,nrow=1)
        base::colnames(windowpeak)<-base::c("repeat length","freq","height", "index", "start", "end")
      }
    }

    if(pspectra){
      if(base::length(rmean_repeat_ave)>500)type<-"l" else type<-"o"
      base::plot((rmean_repeat_ave) , ylab="Ave Normalized Power",xlab="", xaxt="n",cex.axis=1.1,cex.lab=1.5,type=type,pch=15)
      graphics::title(main=main,cex.main=1)   #repeat_col_ave<-NULL
      graphics::mtext("1/Repeat length",side=1,line=4, cex=1)
      if(base::min(repeats_ave,rm.na=TRUE)>500)digitr<-0 else if(base::min(repeats_ave,rm.na=TRUE)>50) digitr<-2 else digitsr<-4
      graphics::axis(BELOW<-1, labels= base::paste0("1/",base::round(repeats_ave,digits=digitsr)),at=1:base::length(rmean_repeat_ave), cex.axis=1,las=2, xlab="1/Repeat #")   #    at=attc,
      for(j in 1:(base::nrow(windowpeak)-1)){
        fundamental<-1.0/(1/windowpeak[j+1,1]-1/windowpeak[j,1])
        base::cat("\nj",j,"peak2",windowpeak[j+1,1],"peak1",windowpeak[j,1],"Fundamental?=",fundamental)
      }
      base::cat("\n")
    }
  } else {base::cat("\n Invalid bplim <=fftlen: Return NULL");rmean_repeat_ave<-NULL }
  base::return(rmean_repeat_ave)
}

#' run_sum_bp
#'
#' Finds peaks in one repeat length
#'
#' @param  input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_sum_bp<-function(repeat_val,All_spec,chromnam,numstd2, numrange= 1500,pflag=FALSE){
  All_spec1<-All_spec
  #All_spec1[base::which(All_spec1<=1e-5)]<-NA

  meanspec<-base::mean(All_spec1[repeat_val, ],na.rm=TRUE)
  stdspec<-stats::sd(All_spec1[repeat_val, ],na.rm=TRUE)
  threshspec<-0 #meanspec+numstd2*stdspec
  bpval<-base::as.numeric(base::colnames(All_spec))
  deltabp<-(bpval[2]-bpval[1])
  deltabp_onebin<-deltabp*numrange
  start_bpseqval<-base::seq(bpval[1],bpval[base::length(bpval)],by= deltabp_onebin)
  end_bpseqval<-start_bpseqval+deltabp_onebin
  end_bpseqval[base::length(end_bpseqval)]<-bpval[base::length(bpval)]
  which_ones<-base::which(All_spec1[repeat_val, ] >=threshspec)
  power_sum<-NULL;power_mean<-NULL; N<-NULL
  for(j in 1:base::length(start_bpseqval)){
    sbpval<-start_bpseqval[j]
    ebpval<-end_bpseqval[j]
    whichval<-base::which(bpval>=sbpval & bpval<=ebpval)
    if(!base::is.na(whichval[1])){
      sval<-whichval[1]
      eval<-whichval[base::length(whichval)]
      which_notNA<-base::which(!base::is.na(All_spec1[repeat_val,sval: eval]))
      bplength_onebin<-base::length(which_notNA)*deltabp
      bp_notNA<-NULL
      for (k in 2: base::length(which_notNA)){
        bp_notNA<-base::c(bp_notNA,bpval[which_notNA[k]]-bpval[which_notNA[k-1]])
      }

      deltabpnotNA<-(bplength_onebin)/1.e6
      which_val<-base::which(which_ones<eval & which_ones>=sval)
      if(base::length(which_notNA)!=0 && length(All_spec1[repeat_val,which_ones[which_val]])>60){
        meanabove<-base::mean(All_spec1[repeat_val,which_ones[which_val]] ,na.rm=TRUE)
        power_mean<-base::c(power_mean,meanabove)
        sumabove<-base::sum(All_spec1[repeat_val,which_ones[which_val]] ,na.rm=TRUE)/deltabpnotNA
        power_sum<-base::c(power_sum,sumabove)
        Nabove<-base::length(All_spec1[repeat_val,which_ones[which_val]] )/deltabpnotNA
        N<-base::c(N,Nabove)
        base::names(power_mean)[base::length(power_mean)]<-sbpval+bplength_onebin/2
      } else{
        meanabove<-NA;power_mean<-base::c(power_mean,NA);power_sum<-base::c(power_sum,NA) ;N<-base::c(N,NA)
        base::names(power_mean)[base::length(power_mean)]<-sbpval+bplength_onebin/2
      }
    } else{
      bplength_onebin<-0
      meanabove<-NA;power_mean<-base::c(power_mean,NA);power_sum<-base::c(power_sum,NA) ;N<-base::c(N,NA)
      base::names(power_mean)[base::length(power_mean)]<-sbpval+bplength_onebin/2
    }
  }
  base::names(power_sum)<-names(N)<-names(power_mean)
  if( !base::is.null(power_sum)){
    if(pflag ){
      graphics::lines(base::as.numeric(names(power_mean)),power_mean,type="l",lwd=3,col="blue",cex=1.5)
    }
    maxspec<-base::max(power_mean,na.rm=TRUE)
    minspec<-base::min(power_mean,na.rm=TRUE)
    mincol<-base::which(minspec==power_mean)[1]
    maxcol<-base::which(maxspec==power_mean)[1]
    min_mean_seqval<-base::as.numeric(names(power_mean[mincol]) )
    max_mean_seqval<-base::as.numeric(names(power_mean[maxcol]) )
    graphics::barplot(power_mean,names.arg = base::names(power_mean),
            main= base::paste("Chromosome",chromnam,"\nMean Power: 1/repeat length (bp) ",repeat_val,"\nmean + ",
                       numstd2,"*sigma= ",base::round(threshspec,digits=1),
                       "\nMid-Points: min at:",min_mean_seqval,"max at:",max_mean_seqval),
            ylim=base::c(threshspec,maxspec),las=2,cex.axis = 0.65,cex.names=0.65)

    maxspec<-base::max(power_sum,na.rm=TRUE)
    minspec<-base::min(power_sum,na.rm=TRUE)
    mincol<-base::which(minspec==power_sum)[1]
    maxcol<-base::which(maxspec==power_sum)[1]       #base::as.numeric
    min_powsum_seqval<-base::as.numeric(names(power_sum[mincol]) )#+deltabp_onebin/2;
    max_powsum_seqval<-base::as.numeric(names(power_sum[maxcol]) )#+deltabp_onebin/2
    graphics::barplot(power_sum,names.arg = base::names(power_sum),
            main= base::paste("Chromosome",chromnam,"\nSum of Power (per Mbp): 1/repeat length (bp) ",repeat_val,"\nmean + ",
                       numstd2,"*sigma= ",base::round(threshspec,digits=1),
                       "\nMid-Points:min at:",min_powsum_seqval,"max at:",max_powsum_seqval),
            ylim=base::c(0,maxspec),las=2,cex.axis = 0.65,cex.names=0.65)

    maxspec<-base::max(N,na.rm=TRUE)
    minspec<-base::min(N,na.rm=TRUE)
    mincol<-base::which(minspec==N)[1]
    maxcol<-base::which(maxspec==N)[1]       #base::as.numeric
    min_N_seqval<-base::as.numeric(names(N[mincol]) )#+deltabp_onebin/2;
    max_N_seqval<-base::as.numeric(names(N[maxcol]) )#+deltabp_onebin/2
    graphics::barplot(N,names.arg = base::names(N),
            main= base::paste("Chromosome",chromnam,"\nNumber above (per Mbp): 1/repeat length (bp) ",repeat_val,"\nmean + ",
                       numstd2,"*sigma= ",base::round(threshspec,digits=1),
                       "\nMid-Points: min at:",min_N_seqval,"max at:",max_N_seqval),
            ylim=base::c(0,maxspec),las=2,cex.axis = 0.65,cex.names=0.65)
  }

  pow_list<-base::list(power_sum=power_sum,power_mean=power_mean,N=N,
                 min_mean_seqval=min_mean_seqval,max_mean_seqval=max_mean_seqval,
                 min_powsum_seqval=min_powsum_seqval,max_powsum_seqval=max_powsum_seqval,
                 min_N_seqval=min_N_seqval,max_N_seqval=max_N_seqval )

  base::return(pow_list)

}

run_sum_bpold<-function(repeat_val,All_spec,chromnam,numstd2, numrange= 1500,pflag=FALSE){
  All_spec1<-All_spec
  All_spec1[base::which(All_spec1==0.0)]<-NA

  meanspec<-base::mean(All_spec1[repeat_val, ],na.rm=TRUE)
  stdspec<-stats::sd(All_spec1[repeat_val, ],na.rm=TRUE)
  threshspec<-meanspec+numstd2*stdspec    #base::length(peak_along_repeat)    numstd<-3

  if(base::length(peak_along_repeat)!=0){
    start_seqval<-base::seq(1, base::ncol(All_spec),by= numrange) #numrange<-100
    which_ones<-base::which(All_spec1[repeat_val, ] >=threshspec)
    lastbp<-sbp<-base::as.numeric(base::colnames(All_spec1)[ base::ncol(All_spec1)])
    power_sum<-NULL;power_mean<-NULL; N<-NULL
    for(sval in start_seqval){   #sval<-41431        #
      eval<-sval+numrange-1
      sbp<-base::as.numeric(base::colnames(All_spec1)[sval])
      if(sval!=start_seqval[base::length(start_seqval)]) {
        ebp<-base::as.numeric(base::colnames(All_spec1)[sval+numrange])
        deltabp_onebin<-ebp-sbp    #find range of single bin
      } else {
        ebp<-lastbp+1
      }
      deltabp<-(ebp-sbp)/1.e6    # normalizes N and power_sum to per Mbp
      which_val<-base::which(which_ones<eval & which_ones>=sval)
      if(base::length(which_val)!=0){
        meanabove<-base::mean(All_spec1[repeat_val,which_ones[which_val]] ,na.rm=TRUE)
        power_mean<-base::c(power_mean,meanabove)
        sumabove<-base::sum(All_spec1[repeat_val,which_ones[which_val]] ,na.rm=TRUE)/deltabp
        power_sum<-base::c(power_sum,sumabove)
        Nabove<-base::length(All_spec1[repeat_val,which_ones[which_val]] )/deltabp
        N<-base::c(N,Nabove)
        base::names(power_mean)[base::length(power_mean)]<-base::colnames(All_spec1)[sval]
      } else{
        meanabove<-NA
      }
    }
    base::names(power_sum)<-names(N)<-names(power_mean)
    if( !base::is.null(power_sum)){
      if(pflag ){
        graphics::lines(base::as.numeric(names(power_mean)),power_mean,type="l",lwd=3,col="blue",cex=1.5)
      }
      maxspec<-base::max(power_mean,na.rm=TRUE)
      minspec<-base::min(power_mean,na.rm=TRUE)
      mincol<-base::which(minspec==power_mean)[1]
      maxcol<-base::which(maxspec==power_mean)[1]      #base::as.numeric
      min_mean_seqval<-base::as.numeric(names(power_mean[mincol]) )+deltabp_onebin/2;
      max_mean_seqval<-base::as.numeric(names(power_mean[maxcol]) )+deltabp_onebin/2
      graphics::barplot(power_mean,names.arg = base::names(power_mean),
              main= base::paste("Chromosome",chromnam,"\nMean Power: 1/repeat length (bp) ",repeat_val,"\nmean + ",
                         numstd2,"*sigma= ",base::round(threshspec,digits=1),
                         "\nmin at:",names(power_mean[mincol]),"max at:",names(power_mean[maxcol])),
              ylim=base::c(threshspec,maxspec),las=2,cex.axis = 0.65,cex.names=0.65)

      maxspec<-base::max(power_sum,na.rm=TRUE)
      minspec<-base::min(power_sum,na.rm=TRUE)
      mincol<-base::which(minspec==power_sum)[1]
      maxcol<-base::which(maxspec==power_sum)[1]       #base::as.numeric
      min_powsum_seqval<-base::as.numeric(names(power_sum[mincol]) )+deltabp_onebin/2;
      max_powsum_seqval<-base::as.numeric(names(power_sum[maxcol]) )+deltabp_onebin/2
      graphics::barplot(power_sum,names.arg = base::names(power_sum),
              main= base::paste("Chromosome",chromnam,"\nSum of Power (per Mbp): 1/repeat length (bp) ",repeat_val,"\nmean + ",
                         numstd2,"*sigma= ",base::round(threshspec,digits=1),
                         "\nmin at:",min_powsum_seqval,"max at:",max_powsum_seqval),
              ylim=base::c(0,maxspec),las=2,cex.axis = 0.65,cex.names=0.65)

      maxspec<-base::max(N,na.rm=TRUE)
      minspec<-base::min(N,na.rm=TRUE)
      mincol<-base::which(minspec==N)[1]
      maxcol<-base::which(maxspec==N)[1]       #base::as.numeric
      min_N_seqval<-base::as.numeric(names(N[mincol]) )+deltabp_onebin/2;
      max_N_seqval<-base::as.numeric(names(N[maxcol]) )+deltabp_onebin/2
      graphics::barplot(N,names.arg = base::names(N),
              main= base::paste("Chromosome",chromnam,"\nNumber above (per Mbp): 1/repeat length (bp) ",repeat_val,"\nmean + ",
                         numstd2,"*sigma= ",base::round(threshspec,digits=1),
                         "\nmin at:",names(N[mincol]),"max at:",names(N[maxcol])),
              ylim=base::c(0,maxspec),las=2,cex.axis = 0.65,cex.names=0.65)
    }
  }
  pow_list<-base::list(power_sum=power_sum,power_mean=power_mean,N=N,
                 min_mean_seqval=min_mean_seqval,max_mean_seqval=max_mean_seqval,
                 min_powsum_seqval=min_powsum_seqval,max_powsum_seqval=max_powsum_seqval,
                 min_N_seqval=min_N_seqval,max_N_seqval=max_N_seqval )

  base::return(pow_list)

}

seq_to_gene<-function(x){
  a1<-x
  a1[base::grep("c",x=x)]<-"g"
  a1[base::grep("g",x=x)]<-"c"
  a1[base::grep("t",x=x)]<-"a"
  a1[base::grep("a",x=x)]<-"t"
  base::return(a1)
}

seq_to_num<-function(x){
  a1<- base::rep(NA,base::length(base::c(x)))
  a1[base::grep("c",x=x)]<-1
  a1[base::grep("g",x=x)]<-2
  a1[base::grep("t",x=x)]<-3
  a1[base::grep("a",x=x)]<-4
  a1[base::grep("n",x=x)]<-0
  base::return(a1)
}

mk_rand<-function(bplength=10,numnuc=4,nuc=base::c("g","c","t","a")){
  ng<-base::sample(1:numnuc, bplength, replace=TRUE)
  graphics::hist(ng)
  gg<-ng
  for (num in 1:numnuc){
    gg[ng==num]<-nuc[num]
  }

  grandom<- base::paste0(gg,collapse="")
  base::print(grandom)
  base::return(gg)
}

mk_palindrome<-function(bplength=10,centre_seglength=0,numnuc=4,nuc=base::c("g","c","t","a"),
                        numnuc_centre=4,nuc_centre=base::c("g","c","t","a")){
  gg<-mk_rand(bplength=bplength,numnuc=numnuc,nuc=nuc)
  if(centre_seglength!=0)centregg<-mk_rand(bplength=centre_seglength,numnuc=numnuc_centre,nuc=nuc_centre) else centregg<-NULL
  ggtruepalindrome<-base::c(gg,gg[base::length(gg):1])
  grandom<- base::paste0(gg,collapse="")
  base::cat("\nmain sequence",grandom,"\n")

  ggtruepalindrome<-base::c(gg,centregg,gg[base::length(gg):1])
  grandomtruepalindrome<- base::paste0(ggtruepalindrome,collapse="")
  base::cat("\ntrue palindrome\n"); base::print(grandomtruepalindrome);base::cat("\n")
  base::return(grandomtruepalindrome)
}

inversion<-function(tdna){
  c_tdna<-complement_strand(tdna,together=TRUE)
  Rc_tdna_sep<-unbase::list(base::strsplit(c_tdna,""))
  Rc_tdna<- base::paste0(Rc_tdna_sep[base::length(Rc_tdna_sep):1],collapse="")
  base::return(Rc_tdna)
}

DNA_Palindrome<-function(tdna,centre_seglength=0,
                         numnuc_centre=4,nuc_centre=base::c("g","c","t","a")){

  Itdna<-inversion(tdna)

  if(centre_seglength!=0)centregg<-mk_rand(bplength=centre_seglength,numnuc=numnuc_centre,nuc=nuc_centre) else centregg<-NULL
  grandom<- base::paste0(tdna,collapse="")
  base::cat("\nmain sequence",grandom,"\n")
  gcentral<- base::paste0(centregg,collapse="")
  if(centre_seglength!=0) base::cat("\ncentral sequence",gcentral,"\n")

  SItdna<-base::c(grandom,gcentral,Itdna)
  base::cat("\nDNA inversion (with central region) palindrome\n"); base::print(SItdna);base::cat("\n")
  base::return(SItdna)
}

true_palindrome<-function(tdna,centre_seglength=0,
                          numnuc_centre=4,nuc_centre=base::c("g","c","t","a"),finversion=FALSE){
  gg<-base::unlist(base::strsplit(tdna,""))
  if(finversion) {
    Itdna<-inversion(tdna)
    Igg<-base::unlist(base::strsplit(Itdna,""))
  } else Igg<-gg

  if(centre_seglength!=0)centregg<-mk_rand(bplength=centre_seglength,numnuc=numnuc_centre,nuc=nuc_centre) else centregg<-NULL
  grandom<- base::paste0(gg,collapse="")
  base::cat("\nmain sequence",grandom,"\n")
  gcentral<- base::paste0(centregg,collapse="")
  if(centre_seglength!=0) base::cat("\ncentral sequence",gcentral,"\n")
  ggtruepalindrome<-base::c(gg,centregg,Igg[base::length(gg):1])
  grandomtruepalindrome<- base::paste0(ggtruepalindrome,collapse="")
  base::cat("\ntrue palindrome\n"); base::print(grandomtruepalindrome);base::cat("\n")
  base::return(grandomtruepalindrome)
}

complement_strand<-function(a,together=FALSE){
  x<-a
  if(together){a1<-base::unlist(base::strsplit(x,""));x<-a1} else a1<-x      #a<- base::paste0(a1,collapse="") x<-c25 together<-TRUE
  a1[base::grep("c",x=x)]<-"g"
  a1[base::grep("g",x=x)]<-"c"
  a1[base::grep("t",x=x)]<-"a"
  a1[base::grep("a",x=x)]<-"t"
  comp<- base::paste0(a1,collapse="")
  base::return(comp)
}

reverse_strand<-function(a,together=TRUE){
  x<-a
  if(together){a1<-base::unlist(base::strsplit(x,""));x<-a1} else a1<-x      #a<- base::paste0(a1,collapse="") x<-c25 together<-TRUE

  a1<-a1[base::length(a1):1]
  comp<- base::paste0(a1,collapse="")
  base::return(comp)
}

imagesub<-function(imagetot,imageave,subval=1,main="",fact=5,flaglog=FALSE){
  imagetot[imagetot==0]<-NA
  sta<-4;
  zmaxtot<-base::max(imagetot,na.rm=TRUE);zmintot<-base::min(imagetot,na.rm=TRUE)
  zlimtot<-base::c(zmintot+(zmaxtot-zmintot)/fact,zmaxtot)
  if(base::length(imagetot)!=0 ){
    if(base::nrow(imagetot)>3*subval){
      freqlim<-1:(base::nrow(imagetot)/subval);

      zmaxtot<-base::max(base::log(imagetot[freqlim,]),na.rm=TRUE);zmintot<-base::mean(base::log(imagetot[freqlim,]),na.rm=TRUE)
      zlimtot<-base::c(zmintot,zmaxtot)

      imagenan(base::log(imagetot[freqlim,]),lasval=2,lnumc=10000,xma=8,xline=4,yma=7,yline=4,
               row_unit="1/repeat length",col_unit="Sequence Window",zlim=zlimtot,
               main= base::paste("LOGARITHM SCALE: all windowssta=",sta,subval,"\n",main))
      rang<-1:(base::length(freqlim)/sta) #top portion
      if(base::length(rang)>3){
        zmaxtot<-base::max(imagetot[rang,],na.rm=TRUE);zmintot<-base::min(imagetot[rang,],na.rm=TRUE)
        zlimtot<-base::c(zmintot+(zmaxtot-zmintot)/fact,zmaxtot)
        imagenan(imagetot[rang,],lasval=2,lnumc=10000,xma=8,xline=4,yma=6,yline=4,outside.below.color='black',
                 row_unit="1/repeat length",col_unit="Sequence Window",zlim=zlimtot,
                 main= base::paste("Top Portion all windows",subval,"fact",fact,"\n",main))
      }
      rang<-(base::length(freqlim)/sta):base::length(freqlim)
      if(base::length(rang)>3){
        zmaxtot<-base::max(imagetot[rang,],na.rm=TRUE);zmintot<-base::min(imagetot[rang,],na.rm=TRUE)
        zlimtot<-base::c(zmintot+(zmaxtot-zmintot)/fact,zmaxtot)
        imagenan(imagetot[rang,],lasval=2,lnumc=10000,xma=8,xline=4,yma=6,yline=4,outside.below.color='black',
                 row_unit="1/repeat length",col_unit="Sequence Window",zlim=zlimtot,
                 main= base::paste("Bottom Portion all windows",subval,"fact",fact,"\n",main))
      }
    }
  }
}

#' largeimagesub
#'
#' Plots fourier spectrum before merging
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
largeimagesub<-function(All_specAve,ofi,rangebp=NULL,rangeseq1=NULL,flaglog=TRUE,
                        main="",pngflag=TRUE,powval=1.2,inmain=""){ # ofi is file without .jpeg attached
  badratio<- 4361/292;
  badratio<-2.5
  badratio1<-1;
  magnif<-magnif1<-2.5

  if(!base::is.null(rangeseq1) & !base::is.numeric(rangeseq1[1]) ){      #All_specAve<-All_spec
    rangeseq<-base::c(base::which(base::colnames(All_specAve)==rangeseq1[1]),base::which(base::colnames(All_specAve)==rangeseq1[2]))
  }
  if(base::is.null(rangeseq1)) rangeseq1<-rangeseq<-base::c(1, base::ncol(All_specAve))

  if(base::is.null(rangeseq)) rangeseq<-base::c(1, base::ncol(All_specAve))
  base::cat("\nrangeseq1,rangeseq",rangeseq1,rangeseq)
  All_specAve[All_specAve==0]<-NA
  repval<-base::unlist(base::strsplit(base::rownames(All_specAve), split="/"))
  repval<-base::as.numeric(repval[base::seq(2,base::length(repval), by=2)])

  whichones<-base::which(repval>=rangebp[1] & repval<rangebp[2])
  newdata<-All_specAve[whichones,rangeseq[1]:rangeseq[2]]
  newrepval<-repval[whichones]
  base::cat("\n rows and columns are:",base::nrow(newdata), base::ncol(newdata))
  minnewrepval<-base::min((newrepval),na.rm=TRUE)
  for (r in 1:base::nrow(newdata)){
    newdata[r,]<-newdata[r,]/((newrepval[r]-2+2 )*1)^powval
  }
  main<- base::paste(main,"\nbp",rangebp[1],"_",rangebp[2],"seq",rangeseq[1],"_",rangeseq[2],flaglog)
  ofil<- base::paste0(ofi,"bp",rangebp[1],"_",rangebp[2],"seq",rangeseq1[1],"_",rangeseq1[2],flaglog,".png")

  if( base::ncol(newdata)<100){
    w<-855;h<-550
    base::cat("\n largeimagesub: normal plot w h",w,h)
    grDevices::png(filename = ofil, width = w, height = h)
    imagenan(newdata,lasval=2,lnumc=10000,topxma=14,xma=16,xline=12,yma=14,yline=12,
             cex.axis = 2,cex.lab = 1.5,cex.main = 1.5,
             row_unit="1/repeat length",col_unit="Sequence Window",
             main= base::paste("NORMALIZED SCALE:",main))
    grDevices::dev.off()
    base::return()
  }

  base::cat("\n writing to",ofil,"\n")
  h<-base::nrow(newdata);w= base::ncol(newdata)
  seqratio<-w/ h
  base::cat("\n badratio seqratio base::nrow(h) base::ncol(w)", badratio,seqratio,h,w,"\n")
  if(w>12000){
    magnif1<-1
    if(w>30000)magnif1<-0.5
  }else {
    magnif1<-magnif
    if(w<1000){
      magnif1<-magnif<-3.5
    }
    if(w<500){
      magnif1<-magnif<-5
      base::cat("\n changing magnif to ",magnif,magnif1)
    }

  }
  if(seqratio>=badratio){
    h<-base::ceiling((w /badratio)*magnif1); w<-base::floor(w*magnif1)

    base::cat("\nseqratio>badratio MAGNIF,new W/H base::nrow(h) base::ncol(w)   ",magnif1,"   ", w/h,h,w,"\n")
  } else {
    if(1/seqratio>badratio1){
      w<-base::ceiling((h /badratio1*magnif1)); h<-base::floor(h*magnif1)
      if(w/h==1 & w>10000){ w<-10000;h<-10000}
      base::cat("\n MAGNIF,new W/H base::nrow(h) base::ncol(w)",magnif, w/h,h,w,"\n")
    } else {
      if((w <500 && h<w) || (h <=500 && w<=h)){
        w<-base::ceiling((w *magnif1)); h<-base::ceiling(h*magnif1); base::cat("\nchanging w h due to image size only\n")
      }

      base::cat("\n1/seqratio<=badratio1 Magnif no change",magnif)
    }

  }

  lnumr<-base::round(base::nrow(newdata)/3) ; lnumc<-base::round( base::ncol(newdata)/5)
  if(lnumr<10)lnumr<-10
  if(lnumc<10)lnumc<-10
  if(flaglog){
    zmaxtot<-base::max(base::log(newdata),na.rm=TRUE);zmintot<-base::mean(base::log(newdata),na.rm=TRUE)
    zlimtot<-base::c(zmintot,zmaxtot)
    base::cat("max min",zlimtot)

    if(magnif!=1 && w>100) {
      if(magnif<3) widths<-base::c(5,5/(5*magnif)) else widths<-base::c(5,5/(1*magnif))
      heights<-base::c(1,0.25)
      base::cat("\nnew widths[1,2] and heights[1,2] are", widths,heights,"\n")
      xma=base::round(10*3*magnif); xline<-base::round(5*2*magnif);yma<-base::round(10*3*magnif);topyma<-base::round(4*2*magnif)
      base::cat("\noutdise png new xma yma,topyma,  xline ",xma, yma,topyma,  xline,"\n")
      cex.axis = (3*magnif); cex.lab = (4*magnif);cex.main=(4*magnif)
    } else{
      widths=base::c(5,1.25); heights=base::c(1,0.5)
      xma=7; xline<-5;yline=5; yma<-7;topyma<-4
      base::cat("\noutside png: xma yma,topyma,  xline ",xma, yma,topyma,  xline,"\n")
      cex.axis = (3); cex.lab = (4);cex.main=(4)

    }
    xma=base::round(10*3*magnif); xline<-base::round(5*2*magnif);yma<-base::round(10*3*magnif);topyma<-base::round(4*2*magnif)
    xma=base::round(10*magnif); xline<-base::round(5*magnif);yma<-base::round(10*magnif);topyma<-base::round(4*magnif)
    cex.axis = (3*magnif); cex.lab = (4*magnif);cex.main=(4*magnif)

    xma=7; xline<-5;yline=5; yma<-7;topyma<-4
    xma=7; xline<-1;yline=1; yma<-7;topyma<-4

    if(pngflag){
      xline<-base::round(30*magnif);yline=base::round(30*magnif)
      widths<-base::c(5,5/(5*magnif))   #image widths and heights
      widths<-base::c(5,5/(7*magnif))   #image widths and heights
      heights<-base::c(1,0.25)
      cex.axis = (2*magnif); cex.lab = (2*magnif);cex.main=(2*magnif)
      if(w> 2000 & w<4000){
        xline<-base::round(2.6*magnif);yline=base::round(2.4*magnif)
        xma=6*magnif;  yma<-6*magnif;topyma<-4*magnif;
        cex.axis = (1*magnif); cex.lab = (1*magnif);cex.main=(1*magnif)
      }
      if(w>= 4000 & w<8000){
        xma=8*magnif;  yma<-8*magnif;topyma<-8*magnif
        cex.axis = (1.5*magnif); cex.lab = (1.5*magnif);cex.main=(1.5*magnif)
      }
      if(w <=2000){

        xline<-base::round(1.8*magnif);yline=base::round(1.8*magnif)
        xma=2.1*magnif;  yma<-2.1*magnif;topyma<-2.1*magnif;
        cex.axis = (0.4*magnif); cex.lab = (0.4*magnif);cex.main=(0.4*magnif)
      }
      if(w <1000){
        xline<-base::round(1*magnif);yline=base::round(1*magnif)
        xma=1.2*magnif;  yma<-1.2*magnif;topyma<-1*magnif;
        cex.axis = (0.18*magnif); cex.lab = (0.18*magnif);cex.main=(0.18*magnif)
      }
      if(w >8000){
        if(w <15000){
          xma=40*magnif;  yma<-40*magnif;topyma<-40*magnif
          cex.axis = (6*magnif); cex.lab = (6*magnif);cex.main=(6*magnif)
          xline<-base::round(30*magnif);yline=base::round(30*magnif)
        }
        if(w >=15000){
          xma=50*magnif;  yma<-50*magnif;topyma<-50*magnif
          cex.axis = (8*magnif); cex.lab = (8*magnif);cex.main=(8*magnif)
          xline<-base::round(40*magnif);yline=base::round(40*magnif)
        }
      }

      base::cat("\npngflag:magnif  cex.axis, cex.lab,cex.main ",magnif,cex.axis, cex.lab,cex.main,"\n")
      base::cat("\npngflag: xma yma,topyma,  xline ",xma, yma,topyma,  xline,"\n")
      base::cat("\npngflag: widths heights",widths, heights, "and w and h for png",w,h,"\n")
      grDevices::png(filename = ofil, width = w, height = h)
    } else{
      widths=base::c(5,1.25); heights=base::c(1,0.5)
      xma=9; xline<-6;yline=7; yma<-9;topyma<-6
      base::cat("\nxma yma,topyma,  xline ",xma, yma,topyma,  xline,"\n")
      cex.axis = (1.2); cex.lab = (1.2);cex.main=(1.2)
    }
    widths=base::c(5,1.25); heights=base::c(1,0.5)
    widths=base::c(5,0.5); heights=base::c(1,0.5)
    base::cat("\nentering imagenan: widths heights",widths, heights, "and w and h for png",w,h,"\n")
    base::cat("\nbefore image: cex.axis, cex.lab,cex.main ",cex.axis, cex.lab,cex.main,"\n")
    imagenan(base::log(newdata),lasval=2,
             xma=xma,xline=xline,yline=yline,yma=yma,topyma=topyma,
             zlim=zlimtot,lnumr=lnumr,lnumc=lnumc, main= base::paste(main,"\n",inmain),
             col_unit="bp Sequence" ,row_unit="1/Repeat length",
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main= cex.main,
             widths=widths, heights=heights)
    if(pngflag) grDevices::dev.off()

  } else {
    zmaxtot<-base::max((newdata),na.rm=TRUE);zmintot<-base::min((newdata),na.rm=TRUE);zsd<-stats::sd((newdata),na.rm=TRUE)
    base::cat("min max sd",zmintot,zmaxtot,zsd)
    zmintot<-zmintot+2*zsd
    zmaxtot<-zmaxtot-7*zsd
    zlimtot<-base::c(zmintot,zmaxtot)
    base::cat("\nrevised min max  sd",zlimtot,zsd)
    imagenan(newdata,lasval=2,
             xma=xma,xline=xline,yma=yma,topyma=topyma,
             zlim=zlimtot,lnumr=lnumr,lnumc=lnumc,main=main,
             xlab="bp Sequence" ,ylab="1/Repeat length",
             cex.axis = 3*magnif, cex.lab = 4*magnif,cex.main=(4*magnif))
    if(pngflag)grDevices::dev.off()
  }

}


#' largeimagesub_NEW
#'
#' Plots fourier spectrum before merging
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
largeimagesub_NEW <-function(All_specAve, fname, chromosome, ofi, rangebp=NULL, rangeseq1=NULL, flaglog=TRUE,
                             main="", pngflag=TRUE,powval=1.2,inmain="", repround=TRUE, part=FALSE, pdf_flag=FALSE){ # ofi is file without .jpeg attached

  # get the range of genome positions to be plotting, if NULL then plot all positions in All_spec_Ave
  if(!base::is.null(rangeseq1) & !base::is.numeric(rangeseq1[1])){  #All_specAve<-All_spec
    rangeseq <- base::c(base::which(base::colnames(All_specAve)==rangeseq1[1]),base::which(base::colnames(All_specAve)==rangeseq1[2]))
  }
  if(base::is.null(rangeseq1)) {
    rangeseq1<-rangeseq<-base::c(1, base::ncol(All_specAve))
  }
  base::cat("\nrangeseq1,rangeseq",rangeseq1,rangeseq)

  # make all zero values NA
  All_specAve[All_specAve==0]<-NA
  All_specAve[,which(colSums(All_specAve, na.rm = FALSE, dims = 1) < 1e-1)] <- NA

  # get the list of repeat lengths
  repval<-as.numeric(stringr::str_split_fixed(base::rownames(All_specAve),"/", 2)[,2])

  # find the repeat lengths in the range we want to plot
  whichones<-base::which(repval>=rangebp[1] & repval<rangebp[2])
  newdata0<-All_specAve[whichones,c(rangeseq[1]:rangeseq[2])]
  newrepval<-as.numeric(repval[whichones])

  # setup plot name and title
  main <- paste0(fname, "_", chromosome)

  ofil<- base::paste0(ofi,"bp",rangebp[1],"_",rangebp[2],"seq",rangeseq1[1],"_",rangeseq1[2],flaglog,".png")

  # scaling of the data to be able to visualize repeats
  for (r in 1:base::nrow(newdata0)){
    newdata0[r,]<-(newdata0[r,]/(newrepval[r])^1.2)
  }

  # for very small sequence files set scale here
  if( base::ncol(newdata0)<100){
    w<-855;h<-550
    grDevices::png(filename = ofil, width = w, height = h)
    # sets the bottom, left, top and right margins
    par(mar=c(15,4.1,4.1,2.1))
    imagenan(log(newdata0),lasval=2,lnumc=10000,topxma=14,xma=16,xline=12,yma=14,yline=12,
             cex.axis = 2,cex.lab = 1.5,cex.main = 1.5,
             row_unit="1/repeat length",col_unit="Sequence Window",
             main= base::paste("NORMALIZED SCALE:",main))
    if(pngflag) grDevices::dev.off()
    base::return()
  }

  # reset names for plotting
  if (rapportools::is.boolean(repround) && repround==TRUE) {rownames(newdata0) <- paste0("1/", base::round(base::as.numeric(stringr::str_split_fixed(base::rownames(newdata0),"/", 2)[,2])))}
  else {rownames(newdata0) <- paste0("1/", base::round(base::as.numeric(stringr::str_split_fixed(base::rownames(newdata0),"/", 2)[,2]), digits=repround))}

  if (rapportools::is.boolean(part) && part==FALSE) {colnames(newdata0) <- round((c(0:(ncol(newdata0)-1))*5000)/1e6)}
  else {colnames(newdata0) <- round((c(0:(ncol(newdata0)-1))*5000)/1e6)+((part-1)*100)}

  h <- base::nrow(newdata0)
  w <- base::ncol(newdata0)

  magnif <- 0.6
  #if(ncol(newdata0)*5000>60e6) magnif <- 0.6
  #if(ncol(newdata0)*5000<60e6) magnif <- 0.9

  lnumr<-base::round(base::ncol(newdata0)/1)
  lnumc<-base::round(base::ncol(newdata0)/200)

  xma=base::round(100*magnif) # x margin
  xline<-base::round(70*magnif) # number of lines from x axis
  yma<-base::round(100*magnif) # y margin
  yline<-base::round(70*magnif) # number of lines away from y axis
  topyma<-base::round(80*magnif) # top margin

  # https://reeddesign.co.uk/test/points-pixels.html

  pixh <- 5000
  pixw <- base::ceiling(w*60e6/(ncol(newdata0)*5000))

  newdata <- newdata0
  # set intensity values
  zmeantot<-base::mean((newdata),na.rm=TRUE)
  zsd<-stats::sd((newdata),na.rm=TRUE)
  zmintot<-zmeantot-(0.8*zsd)+2
  zmaxtot<-zmeantot+(6*zsd)-3

  base::cat("min max sd",zmintot,zmaxtot,zsd)
  zlimtot<-base::c(zmintot,zmaxtot)
  base::cat("\nrevised min max  sd",zlimtot,zsd)

  # plotting
  grDevices::png(filename = ofil, width = pixw, height = pixh)
  imagenan(newdata,lasval=2,
           xma=xma,xline=xline,yline=yline,yma=yma,topyma=topyma,
           zlim=zlimtot,lnumr=lnumr,lnumc=lnumc,main=main,
           col_unit="Genome Position (Mbp)",row_unit="1/Repeat length (1/bp)", zunit="Repeat Abundance",
           cex.axis = (20*magnif), cex.lab =(20*magnif), cex.main=(20*magnif),
           widths=base::c(((pixw-600)/pixw),1-((pixw-600)/pixw)), heights=base::c(1, 0.25))
  grDevices::dev.off()

  if(pdf_flag) {
    pixw_pdf <- pixw/1000
    pixh_pdf <- pixh/1000
    # https://www.rdocumentation.org/packages/grDevices/versions/3.6.2/topics/pdf
    grDevices::pdf(file = base::paste0(ofi,"bp",rangebp[1],"_",rangebp[2],"seq",rangeseq1[1],"_",rangeseq1[2],flaglog,".pdf"), width = pixw_pdf, height = pixh_pdf)
    imagenan(newdata,lasval=2,
             #xma=xma/1000,xline=xline/1000,yline=yline/1000,yma=yma/1000,topyma=topyma/1000,
             zlim=zlimtot,lnumr=lnumr,lnumc=lnumc,main=main,
             col_unit="Genome Position (Mbp)",row_unit="1/Repeat length (1/bp)", zunit="Repeat Abundance")
    grDevices::dev.off()
  }

  #----------------------------
  # try alternative scaling based on row sums - like for Shannon diversity
#
#   All_spec <- as.matrix(All_specAve[whichones,c(rangeseq[1]:rangeseq[2])])
#   All_spec_replen_sum <- apply(All_spec, 1, sum, na.rm=TRUE)
#   All_spec_norm <- All_spec/All_spec_replen_sum
#
#
#   if (rapportools::is.boolean(repround)) {rownames(All_spec_norm) <- paste0("1/", base::round(base::as.numeric(stringr::str_split_fixed(base::rownames(All_spec_norm),"/", 2)[,2])))}
#   else {rownames(All_spec_norm) <- paste0("1/", base::round(base::as.numeric(stringr::str_split_fixed(base::rownames(All_spec_norm),"/", 2)[,2]), digits=repround))}
#
#   if (rapportools::is.boolean(part) && part==FALSE) {colnames(newdata0) <- round((c(0:(ncol(newdata0)-1))*5000)/1e6)}
#   else {colnames(newdata0) <- round((c(0:(ncol(newdata0)-1))*5000)/1e6)+((part-1)*100)}
#
#   ofil1<- base::paste0(ofi,"bp",rangebp[1],"_",rangebp[2],"seq",rangeseq1[1],"_",rangeseq1[2],flaglog,"diff_norm.png")
#
#   # set intensity values
#   # set intensity values
#   zmeantot<-base::mean((All_spec_norm),na.rm=TRUE)
#   zsd<-stats::sd((All_spec_norm),na.rm=TRUE)
#   zmintot<-zmeantot-(0.8*zsd)
#   zmaxtot<-zmeantot+(6*zsd)
#
#   base::cat("min max sd",zmintot,zmaxtot,zsd)
#   zlimtot<-base::c(zmintot,zmaxtot)
#   base::cat("\nrevised min max  sd",zlimtot,zsd)

  # plotting
  # grDevices::png(filename = ofil1, width = pixw, height = pixh)
  # imagenan(All_spec_norm,lasval=2,
  #          xma=xma,xline=xline,yline=yline,yma=yma,topyma=topyma,
  #          zlim=zlimtot,lnumr=lnumr,lnumc=lnumc,main=main,
  #          col_unit="Genome Position (Mbp)",row_unit="1/Repeat length (1/bp)", zunit="Repeat Abundance",
  #          cex.axis = (25*magnif), cex.lab =(25*magnif), cex.main=(25*magnif),
  #          widths=base::c(((pixw-600)/pixw),1-((pixw-600)/pixw)), heights=base::c(1, 0.25))
  # grDevices::dev.off()

}

#' run_chloroplast
#'
#' Creates fourier spectrum of chromosome without plotting
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_chloroplast<-function(nam=nam, fname=fname, outpath=outpath, inpath=inpath,
                          startval=NULL,endval=NULL,fftlength=5000,
                          All_walk=NULL,ptset=NULL,AT_flag=TRUE,
                          sample_every=1,walkflag=FALSE,atflag=TRUE,waxflag=TRUE,
                          printlong=FALSE,printdna=FALSE,pflag=TRUE,plotflag=TRUE,
                          writeflag=TRUE){

  splitname<-base::unlist(stringr::str_split(nam,"_"))
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]

  # create directories
  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)
  base::dir.create(base::file.path(outpath, nam1,chr), showWarnings = FALSE)
  outpathtxt<-base::file.path(outpath, nam1,chr,"spectra.txt")
  outpathpdf<-base::file.path(outpath, nam1,chr,"spectra.pdf")
  outpathwalk<-base::file.path(outpath, nam1,chr,"dnawalk")
  base::dir.create( outpathtxt, showWarnings = FALSE)
  base::dir.create(outpathpdf, showWarnings = FALSE)
  outpathSPECT<-base::file.path(outpath, nam1,chr,"spectra_Table.txt")
  base::dir.create( outpathSPECT, showWarnings = FALSE)
  outpathfractal<-base::file.path(outpath, nam1,chr,"fractalD")
  base::dir.create( outpathfractal, showWarnings = FALSE)
  outpathfractal.bed<-base::file.path(outpath, nam1,chr,"fractalD.bedGraph")
  base::dir.create( outpathfractal.bed, showWarnings = FALSE)
  outpathDNAWALK<-base::file.path(outpath, nam1,chr,"DNAWALK")
  base::dir.create( outpathDNAWALK, showWarnings = FALSE)

  Shortnames<-(base::unlist(base::strsplit(nam," " )))
  if(base::length(Shortnames)>7)Shortnames<-Shortnames[1:7]
  shortname<-NULL; for (sh in Shortnames)shortname<- base::paste0(shortname,"_",sh)
  in_name<- base::paste0(inpath,nam,".fasta")    #
  if(base::is.null(All_walk)){
    dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")

    dn1<-dna_1_to_wax(dna_1)
    if(base::is.null(startval))startval<-1; if(base::is.null(endval))endval<-base::length(dn1)
    base::cat("\nrun_chloroplast : startval endval",startval,endval,"\n")
    dn<-dn1[startval:endval]
    if(printdna){base::cat("\n DNA\n"); base::print( base::paste0(dn,collapse=""))}
  }  else {dn<-NULL}

  base::cat("\nrun_chloroplast : startval endval",startval,endval,"\n")
  dna_1<-NULL; dn1<-NULL; base::gc()
  genname1<-"  AT walk spar 1_"
  main= base::paste0(genname1,"\n",shortname)
  list_plot<-base::c(5,10,20,30,1,2)
  if(printdna){

    if(sample_every!=1){
      sname<- base::file.path(outpathwalk, base::paste0("DNA_WALK_",nam,genname1,startval,"_",endval,"_",sample_every))
    } else{
      sname<-base::file.path(outpathwalk, base::paste0("DNA_WALK_",nam,genname1,startval,"_",endval))
    }
  } else sname<-NULL

  if(sample_every!=1){
    if(plotflag)   grDevices::pdf(base::file.path(outpathpdf, base::paste0(nam,genname1,startval,"_",endval,"_",sample_every,"_",fftlength,".pdf")))
    base::sink(base::file.path(outpathtxt, base::paste0(nam,genname1,startval,"_",endval,"_",sample_every,"_",fftlength,".txt")))
  } else{
    grDevices::pdf(base::file.path(outpathpdf, base::paste0(nam,genname1,startval,"_",endval,"_",fftlength,".pdf")))
    base::sink(base::file.path(outpathtxt, base::paste0(nam,genname1,startval,"_",endval,"_",fftlength,".txt")) )
  }
  # 3.1, 2.2
  main= base::paste(genname1,"\n",shortname," ",sample_every,"_",fftlength)
  spectinfo<-walk_and_plot(tdna=dn,main=main,seqlen=fftlength,sample_every=sample_every,fullname=sname,ptset=ptset,
                           All_walk=All_walk, spar=1,startval=startval,endval=endval,walkflag=walkflag,AT_flag=AT_flag,
                           atflag=atflag,waxflag=waxflag,printlong=printlong,pflag=pflag,plotflag=plotflag,writeflag=writeflag)

  spectimage<-spectinfo[[2]];
  spectmean<-spectinfo[[3]]
  DNAwalk<-spectinfo[[4]]
  spectwalklist<-spectinfo[[5]]
  fractal_list<-spectinfo[[6]]

  if(plotflag){
    for (subval in list_plot){  #subval<-1
      main= base::paste0(startval,"_",endval,"\n",shortname)
      imagesub(spectimage,spectimage,subval=subval,
               main= base::paste("spar=1",main),fact=8)
    }
  }

  if(plotflag) grDevices::dev.off()
  base::sink()

  All_spec<-spectimage
  if(sample_every!=1){
    ofile<-base::file.path(outpathSPECT, base::paste0("concat_ALL_tot_",shortname,"_",
                                         startval,"_",endval,"_",fftlength,"_",sample_every,
                                         "_spar1_Table.txt"))
  } else{
    ofile<-base::file.path(outpathSPECT, base::paste0("concat_ALL_tot_",shortname,"_",
                                         startval,"_",endval,"_",fftlength,
                                         "_spar1_Table.txt"))
  }
  if(writeflag) utils::write.table(x=spectimage, file=ofile, append = FALSE, sep = " ", dec = ".",
                            row.names = TRUE, col.names = TRUE)

  ofile<-base::file.path(outpathDNAWALK, base::paste0("DNA_WALK_",nam,genname1,startval,"_",endval,"_",sample_every,".txt"))

  if(writeflag) utils::write.table(x=DNAwalk, file=ofile, append = FALSE, sep = " ", dec = ".",
                            row.names = TRUE, col.names = TRUE)

  ofile1<-base::file.path(outpathfractal, base::paste0("Fractal_",nam,genname1,startval,"_",endval,"_",sample_every,".txt"))

  if(writeflag)utils::write.table(x=fractal_list, file=ofile1, append = FALSE, sep = " ", dec = ".",
                           row.names = TRUE, col.names = TRUE)
  if(writeflag){
    ofile1<-base::file.path(outpathfractal.bed, base::paste0("Fractal_",nam,genname1,startval,"_",endval,"_",sample_every,
                                                "_",fftlength, ".bedGraph"))
    base::sink(ofile1)
    base::cat("track type=bedGraph")
    for(j in 1:base::length(fractal_list)){
      delta<-base::as.numeric(names(fractal_list)[2])-base::as.numeric(names(fractal_list)[1])  #base::as.numeric(names(fractal)[2])-base::as.numeric(names(fractal)[1])
      base::cat( base::paste0("\n",nam,"\t",base::as.numeric(names(fractal_list)[j]),"\t",
                 base::as.numeric(names(fractal_list)[j])+delta,"\t",
                 fractal_list[j]) )
    }
    base::sink()
  }
  All_list<-base::list(All_spec=All_spec,spectmean=spectmean,DNAwalk=DNAwalk, fractal_list=fractal_list)
  base::return(All_list)
}

#' run_barplots_bedGraph
#'
#' Makes barplots of minimums and maximums for repeats  and bedGraphs of repeat locations
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_barplots_bedGraph <-function(All_spec, nam=nam,fname=fname,inpath=inpath, outpath=outpath,
                                 chromnum=00, atflag=TRUE, numstd=7,  numstd2=3, numrange=200,
                                 repeat_range=NULL,samplesize=1,fftlength=5000){

  xlim<-base::c(base::as.numeric(base::colnames(All_spec)[1]),base::as.numeric(base::colnames(All_spec)[ base::ncol(All_spec)]) )
  chromnam<-nam
  splitname<-base::unlist(stringr::str_split(nam,"_"))   #fname<-"Cannabis"
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]

  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)
  base::dir.create(base::file.path(outpath, nam1,chr), showWarnings = FALSE)
  if((base::as.numeric(base::colnames(All_spec)[2])-base::as.numeric(base::colnames(All_spec)[1]))==5000){
    outpathbarplot<-base::file.path(outpath, nam1,chr,"5000bp_barplots")
    outpathbarplotbedgraph<-base::file.path(outpath, nam1,chr,"5000bp_barplots","bedgraph")
  } else {
    outpathbarplot<-base::file.path(outpath, nam1,chr,"barplots")
    outpathbarplotbedgraph<-base::file.path(outpath, nam1,chr,"barplots","bedgraph")
  }

  base::dir.create( outpathbarplot, showWarnings = FALSE)
  base::dir.create( outpathbarplotbedgraph, showWarnings = FALSE)

  genname1<-"full"
  Shortnames<-(base::unlist(base::strsplit(nam," " )))
  if(base::length(Shortnames)>7)Shortnames<-Shortnames[1:7]
  shortname<-NULL; for (sh in Shortnames)shortname<- base::paste0(shortname,"_",sh)

  in_name<- base::paste0(inpath,"bar",samplesize,"_",fftlength,"_",nam,".txt")


  repeat_list<-base::rownames(All_spec)
  if(!base::is.null(repeat_range)) {
    Nrepeat_list<-base::as.numeric(base::gsub("1/","",repeat_list))
    whichrow<-base::which(Nrepeat_list<=repeat_range[2] & Nrepeat_list>=repeat_range[1])
    repeat_list<-repeat_list[whichrow]
    Nrepeat_list<-Nrepeat_list[whichrow]
  }
  # repeat_val<-"1/416.6667"
  if(atflag){
    chromnam<- base::paste(chromnam,"AT")
    ofile<-base::file.path(outpathbarplot, base::paste0("repeats",
                                           numrange,"_",chromnam,"_s_",numstd2,"std_",numstd,"_",repeat_range[1],"_",repeat_range[2],
                                           "bp.pdf"))
  } else {
    chromnam<- base::paste(chromnam,"CG")
    ofile<-base::file.path(outpathbarplot, base::paste0("repeats_CG",
                                           numrange,"_",chromnam,"_s_",numstd2,"std_",numstd,"_",repeat_range[1],"_",repeat_range[2],
                                           "bp.pdf"))
  }
  base::cat("\n ouput to", ofile)
  grDevices::pdf(file=ofile)
  min_powsum_seqval_list<-NULL; max_powsum_seqval_list<-NULL; goodrepeats<-NULL
  min_mean_seqval_list<-NULL; max_mean_seqval_list<-NULL
  min_N_seqval_list<-NULL; max_N_seqval_list<-NULL
  for( repeat_val in repeat_list){
    base::cat("\n beginning ",repeat_val)
    peak_along_repeat<- run_one_repeat(repeat_val,All_spec,chromnam ,numstd)   #base::colnames(All_spec)
    if(!base::is.null(peak_along_repeat)){
      if(atflag){
        base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_",numrange,"_","_s_",numstd2,"std_",numstd,"_",base::as.numeric(base::gsub("1/",
                                                                                                                                    "",repeat_val)),".bedGraph")))
      } else {
        base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_CG_",numrange,"_","_s_",numstd2,"std_",numstd,"_",base::as.numeric(base::gsub("1/",
                                                                                                                                       "",repeat_val)),".bedGraph")))
      }
      base::cat("track type=bedGraph")
      for(j in 1:base::nrow(peak_along_repeat)){
        base::cat( base::paste0("\n",nam1,"\t",peak_along_repeat[j,"start_peak_bp"],"\t",peak_along_repeat[j,"end_peak_bp"],"\t",
                   peak_along_repeat[j,"Percent of peak"]) )
      }

      pow_list<- run_sum_bp(repeat_val,All_spec,chromnam,numstd2,numrange)

      sum_at_bp<-pow_list$power_sum
      mean_at_bp<-pow_list$power_mean

      goodrepeats<-base::c(goodrepeats,repeat_val)
      min_powsum_seqval_list<-base::c(min_powsum_seqval_list,pow_list$min_powsum_seqval)
      max_powsum_seqval_list<-base::c(max_powsum_seqval_list,pow_list$max_powsum_seqval)

      min_mean_seqval_list<-base::c(min_mean_seqval_list,pow_list$min_mean_seqval)
      max_mean_seqval_list<-base::c(max_mean_seqval_list,pow_list$max_mean_seqval)

      min_N_seqval_list<-base::c(min_N_seqval_list,pow_list$min_N_seqval)
      max_N_seqval_list<-base::c(max_N_seqval_list,pow_list$max_N_seqval)

      N_at_bp<-pow_list$N
      base::sink()
      if(atflag){
        base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_sum_by_bp_std_",numrange,"_s_",numstd2,"std_",
                                                          numstd,"_",base::as.numeric(base::gsub("1/","",repeat_val)),".txt")))
      } else {
        base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_sum_by_bp_std_CG_",numrange,"_s_",numstd2,"std_",
                                                          numstd,"_",base::as.numeric(base::gsub("1/","",repeat_val)),".txt")))
      }
      base::cat("\n1/repeat length ",repeat_val," numrange",numrange,"\nbp sum mean N")
      for(j in 1:base::length(sum_at_bp)) base::cat("\n",names(sum_at_bp)[j],sum_at_bp[j],mean_at_bp[j],N_at_bp[j])
      base::sink()
    }

  }
  if(atflag){
    base::sink(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_POWER_SUM_seqval_",numrange,"_s_",numstd2,"std_",numstd,
                                              "_",repeat_range[1],"_",repeat_range[2],".txt")))
  } else {
    base::sink(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_POWER_SUM_seqval_CG_",numrange,"_s_",numstd2,"std_",numstd,
                                              "_",repeat_range[1],"_",repeat_range[2],".txt")))
  }
  base::cat("\ngoodrepeat min_mean_seqval max_mean_seqval min_powsum_seqval max_powsum_seqval min_N_seqval max_N_seqval\n")
  for(j in 1:base::length(min_powsum_seqval_list)) {
    base::cat(goodrepeats[j],
        min_mean_seqval_list[j],max_mean_seqval_list[j],
        min_powsum_seqval_list[j],max_powsum_seqval_list[j],
        min_N_seqval_list[j],max_N_seqval_list[j],"\n")
  }
  grDevices::dev.off()
  if(atflag){
    grDevices::pdf(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_POWER_SUM_seqval_",numrange,"_s_",numstd2,"std_",numstd,
                                             "_",repeat_range[1],"_",repeat_range[2],".pdf")))
  } else {
    grDevices::pdf(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_POWER_SUM_seqval_CG_",numrange,"_s_",numstd2,"std_",numstd,
                                             "_",repeat_range[1],"_",repeat_range[2],".pdf")))
  }

  base::cat("\n Mean:     Min summary\n");base::print(base::summary(min_mean_seqval_list))
  base::cat("\n Power Sum:Min summary\n");base::print(base::summary(min_powsum_seqval_list))
  base::cat("\n N:        Min summary\n");base::print(base::summary(min_N_seqval_list))
  base::cat("\n Mean:     Max summary\n");base::print(base::summary(max_mean_seqval_list))
  base::cat("\n Power Sum:Max summary\n");base::print(base::summary(max_powsum_seqval_list))
  base::cat("\n N:        Max summary\n");base::print(base::summary(max_N_seqval_list))

  a<-graphics::hist(min_mean_seqval_list,breaks=20, xaxp=base::c(0,base::max(min_mean_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_min_mean_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nMean: Minimum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin Mean:       midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")


  a<-graphics::hist(min_powsum_seqval_list,breaks=20,xaxp=base::c(0,base::max(min_powsum_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_min_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nPower Sum Minimum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin Power Sum:  midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(min_N_seqval_list,breaks=20,xaxp=base::c(0,base::max(min_N_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_min_N_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nNumber: Minimum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin N:          midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<- graphics::hist(max_mean_seqval_list,breaks=20,xaxp=base::c(0,base::max(max_mean_seqval_list),20),xlab="",xlim=xlim,
           main= base::paste0(chromnam,"_max_mean_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                       "\nMean: Maximum values in the Sequence\nbp range ",
                       repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax Mean:      midpointmax counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(max_powsum_seqval_list,breaks=20,xaxp=base::c(0,base::max(max_powsum_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_max_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nPower Sum Maximum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax Power Sum: midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(max_N_seqval_list,breaks=20,xaxp=base::c(0,base::max(max_N_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_max_N_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nNumber: Maximum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax N:         midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")
  base::sink()
  grDevices::dev.off()

}


#' run_barplots
#'
#' Makes barplots of minimums and maximums for repeats  and bedGraphs of repeat locations
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_barplots <-function(All_spec=All_spec, nam=nam,fname=fname,inpath=inpath, outpath=outpath,
                        atflag=TRUE,chromnum=00,numstd=7,  numstd2=3, numrange=200,
                        repeat_range=NULL,samplesize=1,fftlength=5000){

  xlim<-base::c(base::as.numeric(base::colnames(All_spec)[1]),base::as.numeric(base::colnames(All_spec)[ base::ncol(All_spec)]) )
  chromnam<-nam
  splitname<-base::unlist(stringr::str_split(nam,"_"))   #fname<-"Cannabis"
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]
  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)  #creat directory of individual with chromosomes
  base::dir.create(base::file.path(outpath, nam1, chr), showWarnings = FALSE)   # create Chromosome directory

  if((base::as.numeric(base::colnames(All_spec)[2])-base::as.numeric(base::colnames(All_spec)[1]))==5000){
    outpathbarplot<-base::file.path(outpath, nam1,chr,"5000bp_barplots")
    outpathbarplotbedgraph<-base::file.path(outpath, nam1,chr,"5000bp_barplots","bedgraph")
  } else {
    outpathbarplot<-base::file.path(outpath, nam1,chr,"barplots")
    outpathbarplotbedgraph<-base::file.path(outpath, nam1,chr,"barplots","bedgraph")
  }

  base::dir.create( outpathbarplot, showWarnings = FALSE) #create image dir for spectra
  base::dir.create( outpathbarplotbedgraph, showWarnings = FALSE) #create image dir for spectra

  genname1<-"full"                        # 3.1, 2.2
  Shortnames<-(base::unlist(base::strsplit(nam," " )))
  if(base::length(Shortnames)>7)Shortnames<-Shortnames[1:7]
  shortname<-NULL; for (sh in Shortnames)shortname<- base::paste0(shortname,"_",sh)

  in_name<- base::paste0(inpath,"bar",samplesize,"_",fftlength,"_",nam,".txt")    #

  repeat_list<-base::rownames(All_spec)
  if(!base::is.null(repeat_range)) {
    Nrepeat_list<-base::as.numeric(base::gsub("1/","",repeat_list))
    whichrow<-base::which(Nrepeat_list<=repeat_range[2] & Nrepeat_list>=repeat_range[1])
    repeat_list<-repeat_list[whichrow]
    Nrepeat_list<-Nrepeat_list[whichrow]
  }

  if(atflag){
    chromnam<- base::paste(chromnam,"AT")
    ofile<-base::file.path(outpathbarplot, base::paste0("repeats",
                                           numrange,"_",chromnam,"_s_",numstd2,"std_",numstd,"_",repeat_range[1],"_",repeat_range[2],
                                           "bp.pdf"))
  } else {
    chromnam<- base::paste(chromnam,"CG")
    ofile<-base::file.path(outpathbarplot, base::paste0("repeats_CG",
                                           numrange,"_",chromnam,"_s_",numstd2,"std_",numstd,"_",repeat_range[1],"_",repeat_range[2],
                                           "bp.pdf"))
  }
  base::cat("\n ouput to", ofile)
  grDevices::pdf(file=ofile)
  min_powsum_seqval_list<-NULL; max_powsum_seqval_list<-NULL; goodrepeats<-NULL
  min_mean_seqval_list<-NULL; max_mean_seqval_list<-NULL
  min_N_seqval_list<-NULL; max_N_seqval_list<-NULL
  for( repeat_val in repeat_list){
    base::cat("\n beginning ",repeat_val)
    peak_along_repeat<- run_one_repeat(repeat_val,All_spec,chromnam ,numstd)   #base::colnames(All_spec)
    if(!base::is.null(peak_along_repeat)){
      if(atflag){
        base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_",numrange,"_","_s_",numstd2,"std_",numstd,"_",base::as.numeric(base::gsub("1/",
                                                                                                                                    "",repeat_val)),".bedGraph")))
      } else {
        base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_CG_",numrange,"_","_s_",numstd2,"std_",numstd,"_",base::as.numeric(base::gsub("1/",
                                                                                                                                       "",repeat_val)),".bedGraph")))
      }
      base::cat("track type=bedGraph")

      for(j in 1:base::nrow(peak_along_repeat)){
        base::cat( base::paste0("\n",nam1,"\t",peak_along_repeat[j,"start_peak_bp"],"\t",peak_along_repeat[j,"end_peak_bp"],"\t",
                   peak_along_repeat[j,"Percent of peak"]) )
      }
      # perform sum of powers above mean+sigma2 level  (but only for repeats where mean+numstd *sigma holds)
      pow_list<- run_sum_bp(repeat_val,All_spec,chromnam,numstd2,numrange)   #base::colnames(All_spec)

      sum_at_bp<-pow_list$power_sum
      mean_at_bp<-pow_list$power_mean

      goodrepeats<-base::c(goodrepeats,repeat_val)
      min_powsum_seqval_list<-base::c(min_powsum_seqval_list,pow_list$min_powsum_seqval)
      max_powsum_seqval_list<-base::c(max_powsum_seqval_list,pow_list$max_powsum_seqval)

      min_mean_seqval_list<-base::c(min_mean_seqval_list,pow_list$min_mean_seqval)
      max_mean_seqval_list<-base::c(max_mean_seqval_list,pow_list$max_mean_seqval)

      min_N_seqval_list<-base::c(min_N_seqval_list,pow_list$min_N_seqval)
      max_N_seqval_list<-base::c(max_N_seqval_list,pow_list$max_N_seqval)

      N_at_bp<-pow_list$N
      base::sink()
      if(atflag){
        base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_sum_by_bp_std_",numrange,"_s_",numstd2,"std_",
                                                          numstd,"_",base::as.numeric(base::gsub("1/","",repeat_val)),".txt")))
      } else {
        base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_sum_by_bp_std_CG_",numrange,"_s_",numstd2,"std_",
                                                          numstd,"_",base::as.numeric(base::gsub("1/","",repeat_val)),".txt")))
      }
      base::cat("\n1/repeat length ",repeat_val," numrange",numrange,"\nbp sum mean N")
      for(j in 1:base::length(sum_at_bp)) base::cat("\n",names(sum_at_bp)[j],sum_at_bp[j],mean_at_bp[j],N_at_bp[j])
      base::sink()
    }

  }
  if(atflag){
    base::sink(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_POWER_SUM_seqval_",numrange,"_s_",numstd2,"std_",numstd,
                                              "_",repeat_range[1],"_",repeat_range[2],".txt")))
  } else {
    base::sink(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_POWER_SUM_seqval_CG_",numrange,"_s_",numstd2,"std_",numstd,
                                              "_",repeat_range[1],"_",repeat_range[2],".txt")))
  }
  base::cat("\ngoodrepeat min_mean_seqval max_mean_seqval min_powsum_seqval max_powsum_seqval min_N_seqval max_N_seqval\n")
  for(j in 1:base::length(min_powsum_seqval_list)) {
    base::cat(goodrepeats[j],
        min_mean_seqval_list[j],max_mean_seqval_list[j],
        min_powsum_seqval_list[j],max_powsum_seqval_list[j],
        min_N_seqval_list[j],max_N_seqval_list[j],"\n")
  }
  grDevices::dev.off()
  if(atflag){
    grDevices::pdf(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_POWER_SUM_seqval_",numrange,"_s_",numstd2,"std_",numstd,
                                             "_",repeat_range[1],"_",repeat_range[2],".pdf")))
  } else {
    grDevices::pdf(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_POWER_SUM_seqval_CG_",numrange,"_s_",numstd2,"std_",numstd,
                                             "_",repeat_range[1],"_",repeat_range[2],".pdf")))
  }

  base::cat("\n Mean:     Min summary\n");base::print(base::summary(min_mean_seqval_list))
  base::cat("\n Power Sum:Min summary\n");base::print(base::summary(min_powsum_seqval_list))
  base::cat("\n N:        Min summary\n");base::print(base::summary(min_N_seqval_list))
  base::cat("\n Mean:     Max summary\n");base::print(base::summary(max_mean_seqval_list))
  base::cat("\n Power Sum:Max summary\n");base::print(base::summary(max_powsum_seqval_list))
  base::cat("\n N:        Max summary\n");base::print(base::summary(max_N_seqval_list))


  a<-graphics::hist(min_powsum_seqval_list,breaks=20,xaxp=base::c(0,base::max(min_powsum_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_min_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nPower Sum Minimum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin Power Sum:  midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(min_mean_seqval_list,breaks=20, xaxp=base::c(0,base::max(min_mean_seqval_list),20),xlab="",xlim=xlim,
                    main= base::paste0(chromnam,"_min_mean_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                                       "\nMean: Minimum values in the Sequence\nbp range ",
                                       repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin Mean:       midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(min_N_seqval_list,breaks=20,xaxp=base::c(0,base::max(min_N_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_min_N_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nNumber: Minimum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin N:          midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<- graphics::hist(max_mean_seqval_list,breaks=20,xaxp=base::c(0,base::max(max_mean_seqval_list),20),xlab="",xlim=xlim,
           main= base::paste0(chromnam,"_max_mean_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                       "\nMean: Maximum values in the Sequence\nbp range ",
                       repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax Mean:      midpointmax counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(max_powsum_seqval_list,breaks=20,xaxp=base::c(0,base::max(max_powsum_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_max_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nPower Sum Maximum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax Power Sum: midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(max_N_seqval_list,breaks=20,xaxp=base::c(0,base::max(max_N_seqval_list),20),xlab="",xlim=xlim,
          main= base::paste0(chromnam,"_max_N_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                      "\nNumber: Maximum values in the Sequence\nbp range ",
                      repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax N:         midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")
  base::sink()
  grDevices::dev.off()

}

#' run_barplots_chromosome
#'
#' Makes barplots of minimums and maximums for repeats  and bedGraphs of repeat locations
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_barplots_chromosome <-function(All_spec=All_spec, chromosome=chromosome,fname=fname,inpath=inpath, outpath=outpath,
                                   full_length=full_length, atflag=TRUE,chromnum=00,numstd=3,  numstd2=1, numrange=200,
                                   repeat_range=NULL,samplesize=1,fftlength=5000, binnum=20){

  xlim<-base::c(base::as.numeric(base::colnames(All_spec)[1]),base::as.numeric(base::colnames(All_spec)[ base::ncol(All_spec)]) )
  nam <- paste0(fname, "_", chromosome)
  chromnam<-nam
  splitname<-base::unlist(stringr::str_split(nam,"_"))   #fname<-"Cannabis"
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]
  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)  #creat directory of individual with chromosomes
  base::dir.create(base::file.path(outpath, nam1, chr), showWarnings = FALSE)   # create Chromosome directory

  if((base::as.numeric(base::colnames(All_spec)[2])-base::as.numeric(base::colnames(All_spec)[1]))==5000){
    outpathbarplot<-base::file.path(outpath, nam1,chr,"histograms")
    outpathbarplotbedgraph<-base::file.path(outpath, nam1,chr,"histograms","bedgraph")
  } else {
    outpathbarplot<-base::file.path(outpath, nam1,chr,"histograms")
    outpathbarplotbedgraph<-base::file.path(outpath, nam1,chr,"histograms","bedgraph")
  }

  base::dir.create( outpathbarplot, showWarnings = FALSE) #create image dir for spectra
  base::dir.create( outpathbarplotbedgraph, showWarnings = FALSE) #create image dir for spectra

  genname1<-"full"                        # 3.1, 2.2
  Shortnames<-(base::unlist(base::strsplit(nam," " )))
  if(base::length(Shortnames)>7)Shortnames<-Shortnames[1:7]
  shortname<-NULL; for (sh in Shortnames)shortname<- base::paste0(shortname,"_",sh)

  in_name<- base::paste0(inpath,"bar",samplesize,"_",fftlength,"_",nam,".txt")    #

  repeat_list<-base::rownames(All_spec)
  if(!base::is.null(repeat_range)) {
    Nrepeat_list<-base::as.numeric(base::gsub("1/","",repeat_list))
    whichrow<-base::which(Nrepeat_list<=repeat_range[2] & Nrepeat_list>=repeat_range[1])
    repeat_list<-repeat_list[whichrow]
    Nrepeat_list<-Nrepeat_list[whichrow]
  }

  if(atflag){
    #chromnam<- base::paste(chromnam,"AT")
    ofile<-base::file.path(outpathbarplot, base::paste0("repeats",
                                                        numrange,"_",chromnam,"_s_",numstd2,"std_",numstd,"_",repeat_range[1],"_",repeat_range[2],
                                                        "bp.pdf"))
  } else {
    #chromnam<- base::paste(chromnam,"CG")
    ofile<-base::file.path(outpathbarplot, base::paste0("repeats_CG",
                                                        numrange,"_",chromnam,"_s_",numstd2,"std_",numstd,"_",repeat_range[1],"_",repeat_range[2],
                                                        "bp.pdf"))
  }
  base::cat("\n ouput to", ofile)
  grDevices::pdf(file=ofile)

  min_powsum_seqval_list<-NULL; max_powsum_seqval_list<-NULL; goodrepeats<-NULL
  min_mean_seqval_list<-NULL; max_mean_seqval_list<-NULL
  min_N_seqval_list<-NULL; max_N_seqval_list<-NULL
  for( repeat_val in repeat_list){
    base::cat("\n beginning ",repeat_val)
    # peak_along_repeat<- run_one_repeat(repeat_val,All_spec,chromnam ,numstd)   #base::colnames(All_spec)
    # if(!base::is.null(peak_along_repeat)){
    # if(atflag){
    # base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_AT_",numrange,"_","_s_",numstd2,"std_",numstd,"_",base::as.numeric(base::gsub("1/",
    # "",repeat_val)),".bedGraph")))
    # } else {
    # base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_CG_",numrange,"_","_s_",numstd2,"std_",numstd,"_",base::as.numeric(base::gsub("1/",
    # "",repeat_val)),".bedGraph")))
    # }
    # base::cat("track type=bedGraph")

    # for(j in 1:base::nrow(peak_along_repeat)){
    #   base::cat( base::paste0("\n",nam1,"\t",peak_along_repeat[j,"start_peak_bp"],"\t",peak_along_repeat[j,"end_peak_bp"],"\t",
    #                           peak_along_repeat[j,"Percent of peak"]) )
    # }
    # perform sum of powers above mean+sigma2 level  (but only for repeats where mean+numstd *sigma holds)
    pow_list<- run_sum_bp(repeat_val,All_spec,chromnam,numstd2,numrange)   #base::colnames(All_spec)

    sum_at_bp<-pow_list$power_sum
    mean_at_bp<-pow_list$power_mean

    goodrepeats<-base::c(goodrepeats,repeat_val)
    min_powsum_seqval_list<-base::c(min_powsum_seqval_list,pow_list$min_powsum_seqval)
    max_powsum_seqval_list<-base::c(max_powsum_seqval_list,pow_list$max_powsum_seqval)

    min_mean_seqval_list<-base::c(min_mean_seqval_list,pow_list$min_mean_seqval)
    max_mean_seqval_list<-base::c(max_mean_seqval_list,pow_list$max_mean_seqval)

    min_N_seqval_list<-base::c(min_N_seqval_list,pow_list$min_N_seqval)
    max_N_seqval_list<-base::c(max_N_seqval_list,pow_list$max_N_seqval)

    N_at_bp<-pow_list$N
    base::sink()

    if(atflag){
      base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_sum_by_bp_std_",numrange,"_s_",numstd2,"std_",
                                                                           numstd,"_",base::as.numeric(base::gsub("1/","",repeat_val)),".txt")))
    } else {
      base::sink(file=base::file.path(outpathbarplotbedgraph, base::paste0(chromnam,"_sum_by_bp_std_CG_",numrange,"_s_",numstd2,"std_",
                                                                           numstd,"_",base::as.numeric(base::gsub("1/","",repeat_val)),".txt")))
    }
    base::cat("\n1/repeat length ",repeat_val," numrange",numrange,"\nbp sum mean N")
    for(j in 1:base::length(sum_at_bp)) base::cat("\n",names(sum_at_bp)[j],sum_at_bp[j],mean_at_bp[j],N_at_bp[j])
    base::sink()
  }
  grDevices::dev.off()

  if(atflag){
    base::sink(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_Histogram_input_",numrange,"_s_",numstd2,"std_",numstd,
                                                                 "_",repeat_range[1],"_",repeat_range[2],".txt")))
  } else {
    base::sink(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_Histogram_input_CG_",numrange,"_s_",numstd2,"std_",numstd,
                                                                 "_",repeat_range[1],"_",repeat_range[2],".txt")))
  }

  # add single value at end of chromosome to fix plotting - remove any values greater than or less than 1Mbp from ends
  min_powsum_seqval_list <- c(0, min_powsum_seqval_list, full_length)
  min_mean_seqval_list <- c(0, min_mean_seqval_list, full_length)
  min_N_seqval_list <- c(0, min_N_seqval_list, full_length)
  max_mean_seqval_list <- c(0, max_mean_seqval_list, full_length)
  max_powsum_seqval_list <- c(0, max_powsum_seqval_list, full_length)
  max_N_seqval_list <- c(0, max_N_seqval_list, full_length)
  goodrepeats <- c(0, goodrepeats, full_length)

  base::cat("goodrepeat min_mean_seqval max_mean_seqval min_powsum_seqval max_powsum_seqval min_N_seqval max_N_seqval\n")
  for(j in 1:base::length(min_powsum_seqval_list)) {
    base::cat(goodrepeats[j],
              min_mean_seqval_list[j],max_mean_seqval_list[j],
              min_powsum_seqval_list[j],max_powsum_seqval_list[j],
              min_N_seqval_list[j],max_N_seqval_list[j],"\n")
  }
  base::sink()

  # write all histograms out to pdf file
  if(atflag){
    grDevices::pdf(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_Histograms_",numrange,"_s_",numstd2,"std_",numstd,
                                                                     "_",repeat_range[1],"_",repeat_range[2],".pdf")))
  } else {
    grDevices::pdf(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_Histogram_CG_",numrange,"_s_",numstd2,"std_",numstd,
                                                                     "_",repeat_range[1],"_",repeat_range[2],".pdf")))
  }

  # make summary file with mean, min, max of the histograms
  if(atflag){
    base::sink(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_Histogram_Summary_",numrange,"_s_",numstd2,"std_",numstd,
                                                                 "_",repeat_range[1],"_",repeat_range[2],".txt")))
  } else {
    base::sink(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_Histogram_Summary_CG_",numrange,"_s_",numstd2,"std_",numstd,
                                                                 "_",repeat_range[1],"_",repeat_range[2],".txt")))
  }

  base::cat("\n Mean:     Min summary\n");base::print(base::summary(min_mean_seqval_list))
  base::cat("\n Power Sum:Min summary\n");base::print(base::summary(min_powsum_seqval_list))
  base::cat("\n N:        Min summary\n");base::print(base::summary(min_N_seqval_list))
  base::cat("\n Mean:     Max summary\n");base::print(base::summary(max_mean_seqval_list))
  base::cat("\n Power Sum:Max summary\n");base::print(base::summary(max_powsum_seqval_list))
  base::cat("\n N:        Max summary\n");base::print(base::summary(max_N_seqval_list))

  aips<-graphics::hist(min_powsum_seqval_list,breaks=binnum,xaxp=base::c(0,full_length,binnum),xlab="",xlim=xlim,
                    main= base::paste0(chromnam,"_min_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                                       "\nPower Sum Minimum values in the Sequence\nbp range ",
                                       repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin Power Sum:  midpoint max counts",aips$mids[base::which(aips$counts==base::max(aips$counts))],"\n")

  a<-graphics::hist(min_mean_seqval_list,breaks=binnum, xaxp=base::c(0,full_length,binnum),xlab="",xlim=xlim,
                    main= base::paste0(chromnam,"_min_mean_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                                       "\nMean: Minimum values in the Sequence\nbp range ",
                                       repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin Mean:       midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(min_N_seqval_list,breaks=binnum,xaxp=base::c(0,full_length,binnum),xlab="",xlim=xlim,
                    main= base::paste0(chromnam,"_min_N_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                                       "\nNumber: Minimum values in the Sequence\nbp range ",
                                       repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMin N:          midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<- graphics::hist(max_mean_seqval_list,breaks=binnum,xaxp=base::c(0,full_length,binnum),xlab="",xlim=xlim,
                     main= base::paste0(chromnam,"_max_mean_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                                        "\nMean: Maximum values in the Sequence\nbp range ",
                                        repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax Mean:      midpointmax counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(max_powsum_seqval_list,breaks=binnum,xaxp=base::c(0,full_length,binnum),xlab="",xlim=xlim,
                    main= base::paste0(chromnam,"_max_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                                       "\nPower Sum Maximum values in the Sequence\nbp range ",
                                       repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax Power Sum: midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  a<-graphics::hist(max_N_seqval_list,breaks=binnum,xaxp=base::c(0,full_length,binnum),xlab="",xlim=xlim,
                    main= base::paste0(chromnam,"_max_N_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                                       "\nNumber: Maximum values in the Sequence\nbp range ",
                                       repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  base::cat("\nMax N:         midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")

  base::sink()
  grDevices::dev.off() # histograms pdf

  # write out the barplot centromere value
  maxcount_min_powsum <- aips$mids[base::which(aips$counts==base::max(aips$counts))]
  barplot_cent_data <- t(as.matrix(c(fname, chromosome, maxcount_min_powsum, full_length)))
  utils::write.table(barplot_cent_data, file=file.path(outpathbarplot, base::paste0("Centromere_",chromnam,"_MIN_POWER_SUM_",numrange,"_s_",numstd2,"std_",numstd,
                                                                                  "_",repeat_range[1],"_",repeat_range[2],".txt")),
                     append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # plot power sum min in png image
  grDevices::png(file=base::file.path(outpathbarplot, base::paste0(chromnam,"_histogram_POWER_SUM_seqval_",numrange,"_s_",numstd2,"std_",numstd,
                                                                   "_",repeat_range[1],"_",repeat_range[2],".png")), width = 1500, height = 500)
  graphics::hist(min_powsum_seqval_list,breaks=40,xaxp=base::c(0,full_length,40),xlab="",xlim=xlim,
                       main= base::paste0(chromnam,"_min_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
                                          "\nPower Sum Minimum values in the Sequence\nbp range ",
                                          repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
  grDevices::dev.off()
}


#' run_chromosome
#'
#' Creates fourier spectrum of chromosome with plotting
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_chromosome<-function(nam=nam, fname=fname, inpath=inpath, outpath=outpath,
                         startgroup=1, atflag=TRUE, AT_flag=TRUE,
                         majorlen=5000000,
                         length_majorgroup=1000000, name_majorgroup="1Mbp",
                         length_submajor=500000, name_submajorgroup="500Kbp",
                         length_fftgroup=5000,   name_fftgroup="5Kbp",
                         length_minor= 25000, submajor_nam="25Kb",pflag=FALSE,
                         plotflag=TRUE,writeflag=TRUE,
                         samplesize=50){

  splitname<-base::unlist(stringr::str_split(nam,"_"))
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]
  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)
  base::dir.create(base::file.path(outpath, nam1,chr), showWarnings = FALSE)
  outpathtxt<-base::file.path(outpath, nam1,chr,"5000bp_spectra.txt")
  base::dir.create(outpathtxt, showWarnings = FALSE)
  outpathpdf<-base::file.path(outpath, nam1,chr,"5000bp_spectra.pdf")
  base::dir.create(outpathpdf, showWarnings = FALSE)
  outpathfractal<-base::file.path(outpath, nam1,chr,"5000bp_fractal.pdf")
  base::dir.create(outpathfractal, showWarnings = FALSE)
  outpathspectra<-base::file.path(outpath, nam1,chr,"spectra_Table.txt")
  base::dir.create(outpathspectra, showWarnings = FALSE)
  outpathfractal.bed<-base::file.path(outpath, nam1,chr,"5000bpfractalD.bedGraph")
  base::dir.create( outpathfractal.bed, showWarnings = FALSE)
  outpathDNAWALK<-base::file.path(outpath, nam1,chr,"DNAWALK")
  base::dir.create( outpathDNAWALK, showWarnings = FALSE)
  genname1<-"full"                        # 3.1, 2.2
  Shortnames<-(base::unlist(base::strsplit(nam," " )))
  if(base::length(Shortnames)>7)Shortnames<-Shortnames[1:7]
  shortname<-NULL; for (sh in Shortnames)shortname<- base::paste0(shortname,"_",sh)
  in_name<- base::paste0(inpath,nam,".fasta")    #
  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")
  main= base::paste0(genname1,"\n",shortname)
  dn1<-dna_1_to_wax(dna_1)
  list_plot<-base::c(5,10,20,30,1,2)
  full_length<-base::length(dn1)
  start_seq<-base::seq(startgroup,full_length, by=majorlen)
  DNAwalk_FULL<-NULL
  All_list<-NULL;All_listAve<-NULL; fullwalklist<-NULL
  fractal_dimlist<-NULL
  for(startval in start_seq){
    DNAwalk_extended<-NULL
    endval<-startval+majorlen -1
    dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")
    dn1<-dna_1_to_wax(dna_1)
    dn<-dn1[startval:endval]
    dna_1<-dn1<-NULL
    base::gc()
    #put below in 20Mb  loop
    majorlen<- base::length(dn)
    #method to make list of data for major groups etc
    wholelist<-base::seq(1,base::length(dn), by=length_majorgroup )
    #run through whole chromosome
    major_concatenationtotlist<-NULL
    major_concatenationavelist<-NULL
    concat_ALL_tot<-NULL
    concat_ALL_ave<-NULL
    for (startmajor in wholelist){
      endmajor<-startmajor+length_majorgroup-1
      majorlist<- base::seq(startmajor,startmajor+(length_majorgroup-length_submajor), by=length_submajor)
      concat_majorlist_tot<-NULL
      concat_majorlist_ave<-NULL
      for(val1 in majorlist){
        val2<-val1+length_submajor-1
        if(val2<=endval){
          tdna<-dn[val1:val2]
          submajor_num<-base::floor((length_submajor)/length_minor)
          submajorlist<- base::seq(1,(length_submajor), by=length_minor)
          stval<-val1
          sname<-base::file.path(outpathDNAWALK, base::paste0("DNA_WALK_",fname,shortname,"_",startval+val1-1,"_",startval+val1+length_submajor-2,"_",
                                                 submajor_nam,"_",length_fftgroup,"_spar1"))
          if(plotflag){
            if(atflag){
              grDevices::pdf(base::file.path(outpathpdf, base::paste0(fname,shortname,"_",
                                              startval+val1-1,"_",startval+val1+length_submajor-2,"_",submajor_nam,"_",length_fftgroup,
                                              "_spar1.pdf")))

              ofi<-base::file.path(outpathtxt, base::paste0(fname,shortname,"_",
                                               startval+val1-1,"_",startval+val1+length_submajor-2,"_",submajor_nam,"_",length_fftgroup,
                                               "_spar1.txt"))
            } else{
              grDevices::pdf(base::file.path(outpathpdf, base::paste0(fname,shortname,"_",
                                              startval+val1-1,"_",startval+val1+length_submajor-2,"_",submajor_nam,"_",length_fftgroup,
                                              "_spar1_CG.pdf")))

              ofi<-base::file.path(outpathtxt, base::paste0(fname,shortname,"_",
                                               startval+val1-1,"_",startval+val1+length_submajor-2,"_",submajor_nam,"_",length_fftgroup,
                                               "_spar1_CG.txt"))
            }
            base::cat("\nrun_chromosome: output going to",ofi,"\n")
            base::sink(file=ofi)
          }
          imagetot<-NULL;imageave<-NULL
          for (j in submajorlist){
            inc<-length_minor
            fullname<- base::paste0(sname,"_",(startval+stval+j-2),"_",(startval+stval+j+inc-3) )
            if((j+inc-1)<=(length_submajor)){
              genname1<- base::paste0("rng=base::c(",j,"_" ,j+inc-1,")"," ",(startval+stval+j-2),"_" ,(startval+stval+j+inc-3) ,")")
              main= base::paste(genname1,"\n",shortname)
              spectinfo<-walk_and_plot(tdna=tdna,main=main,rng=base::c(j,j+inc-1),seqlen=length_fftgroup,
                                       fullname=fullname, spar=1,startval=(startval+stval+j-2),walkflag=FALSE,waxflag=TRUE,
                                       pflag=FALSE,plotflag=plotflag,writeflag=writeflag,atflag=atflag, AT_flag=AT_flag)
            } else {
              genname1<- base::paste0("rng=base::c(",j,"_" ,length_submajor,")"," ",(startval+stval+j-2),"_" ,(startval+stval+length_submajor) ,")")
              main= base::paste(genname1,"\n",shortname)
              spectinfo<-walk_and_plot(tdna,main=main,rng=base::c(j,j+inc-1),seqlen=5000,
                                       fullname=fullname, spar=1,startval=(startval+stval+j-1),walkflag=FALSE,waxflag=TRUE,
                                       plotflag=plotflag,writeflag=writeflag, atflag=atflag, AT_flag=AT_flag)

            }
            peaks<-spectinfo[[1]];
            spectimage<-spectinfo[[2]];
            spectmean<-spectinfo[[3]]
            DNAwalk<-spectinfo[[4]]
            spectwalklist<-spectinfo[[5]]
            fractalD<-spectinfo[[6]]
            fractal_dimlist<-base::c(fractal_dimlist,fractalD)
            if(writeflag){
              base::cat("\n\n\n ", base::paste0("Final peaks\nrng=base::c(",j,"_" ,j+inc-1,")"," ",
                                   (startval+stval+j),"_" ,(startval+stval+j+inc-1) ,")"),"\n")
              if(base::nrow(peaks)>=20) base::print(peaks[1:20,]) else base::print(peaks[1:base::nrow(peaks),])

            }
            imagetot<-base::cbind(imagetot,spectimage) #imagenan(imagetot)     base::colnames(spectinfo[[2]])

            if(writeflag){
              base::cat("\nstart dwalk calc: DNAwalk\n")
              base::print(DNAwalk[1:10,])
            }
            if(!base::is.null(DNAwalk_extended)){
              if(atflag){
                DNAwalk_extended<- base::rbind(DNAwalk_extended,
                                        base::cbind(DNAwalk_extended[base::nrow(DNAwalk_extended),"AT"]+DNAwalk[,"AT"],
                                              DNAwalk_extended[base::nrow(DNAwalk_extended),"CG"]+DNAwalk[,"CG"]))
              } else {
                DNAwalk_extended<- base::rbind(DNAwalk_extended,
                                        base::cbind(DNAwalk_extended[base::nrow(DNAwalk_extended),"CG"]+DNAwalk[,"CG"],
                                              DNAwalk_extended[base::nrow(DNAwalk_extended),"AT"]+DNAwalk[,"AT"]))
              }
            } else {
              DNAwalk_extended<-DNAwalk
            }
          }

          #plot images for 0.5 Mbp intervals
          if(plotflag){
            if(base::max(imagetot,na.rm=TRUE)>base::min(imagetot,na.rm=TRUE)){
              for (subval in list_plot){  #subval<-1
                main= base::paste0(startval+val1,"_",startval+val1+length_submajor,"\n",shortname)
                imagesub(imagetot,imageave,subval=subval,
                         main= base::paste("spar=1",main),fact=8)
              }
            }
          }

          if(plotflag){
            base::sink()
            grDevices::dev.off()
          }
          #concatenate the 500Kbp images/data
          if(base::max(imagetot,na.rm=TRUE)>base::min(imagetot,na.rm=TRUE)){
            concat_majorlist_tot<-base::cbind(concat_majorlist_tot,imagetot)
            #concat_majorlist_ave<-base::cbind(concat_majorlist_ave,imageave)
          }
        }
      } #end of length_submajor loop  (0.5Mbp)

      if(plotflag){
        if(!base::is.null(base::nrow(concat_majorlist_tot))){
          if(atflag){
            grDevices::pdf(base::file.path(outpathpdf, base::paste0("Concatenate_",fname,shortname,"_",
                                            startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                            "_spar1.pdf")))
            base::sink(file=base::file.path(outpathtxt, base::paste0("Concatenate_",fname,shortname,"_",
                                                  startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                                  "_spar1.txt")))    #total number of bp is 295045
            #plot concatenated data for major group (1.5Mbp)
            main<- base::paste0("Concatenate ",shortname,"_",
                         startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"\n",
                         name_majorgroup,submajor_nam,"_",length_fftgroup)
          } else {
            grDevices::pdf(base::file.path(outpathpdf, base::paste0("Concatenate_CG_",fname,shortname,"_",
                                            startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                            "_spar1.pdf")))
            base::sink(file=base::file.path(outpathtxt, base::paste0("Concatenate_CG_",fname,shortname,"_",
                                                  startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                                  "_spar1.txt")))    #total number of bp is 295045
            #plot concatenated data for major group (1.5Mbp)
            main<- base::paste0("Concatenate CG",shortname,"_",
                         startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"\n",
                         name_majorgroup,submajor_nam,"_",length_fftgroup)
          }

          base::cat("\nrun_chromosome: entering imagesub loop for Concatenate spectra","\n")
          for (subval in list_plot){    # base::nrow(concat_majorlist_tot)
            imagesub(concat_majorlist_tot,imageave=NULL,subval=subval,
                     main= base::paste("spar=1",main),fact=8)                         #? ,sval=startval+startmajor
          }
          base::cat("\nrun_chromosome: exiting imagesub loop for Concatenate spectra","\n")
          base::sink()
          grDevices::dev.off()
          if(atflag){
            ofile<-base::file.path(outpathfractal, base::paste0("FracCon_",fname,shortname,"_",
                                                   startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                                   "_spar1.txt"))

            #fractal_dimlist for1.5Mbp
            s1<-base::which(base::as.numeric(names(fractal_dimlist))>=startval+startmajor-1)[1]
            base::cat("\nrun_chromosome: Writing to",ofile,"\n Starting at s1",s1,"\n")
            utils::write.table(x=fractal_dimlist[s1:base::length(fractal_dimlist)], file=ofile, append = FALSE, sep = " ", dec = ".",
                        row.names = TRUE, col.names = TRUE)


            ofile<-base::file.path(outpathfractal, base::paste0("FraCon_",fname,shortname,"_",
                                                   startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                                   "_spar1.pdf"))
            base::cat("\nrun_chromosome: Writing to",ofile)
            grDevices::pdf(ofile)
            minf<-base::as.numeric(names(fractal_dimlist)[s1])
            maxf<-base::as.numeric(names(fractal_dimlist)[base::length(fractal_dimlist)]);
            numrange= base::round((maxf-minf)/(15*length_fftgroup))  # length_fftgroup
            base::cat("\nnumrange minf maxf set to ",numrange,minf,maxf,"Enter run_sum_Fractal\n")
            run_sum_Fractal(fracdim_list=fractal_dimlist[s1:base::length(fractal_dimlist)], numrange= numrange,pflag=TRUE,main=main)
            grDevices::dev.off()
          }
        }   #end test for null concat_majorlist_tot
      }
      list2<-base::list(concat_majorlist_tot)
      major_concatenationtotlist<- base::c(major_concatenationtotlist, list2=list2)
      base::names(major_concatenationtotlist)[base::length(major_concatenationtotlist)]<- base::paste0(fname,
                                                                                    shortname,"_",
                                                                                    startmajor,"_",startmajor+length_majorgroup-1,"_",submajor_nam,"_",length_fftgroup, "_spar1")

      #plot all together
      concat_ALL_tot<-base::cbind(concat_ALL_tot,concat_majorlist_tot)
      #concat_ALL_ave<-base::cbind(concat_ALL_ave,concat_majorlist_ave)
      base::cat("\n run_chromosome: Ending :",startmajor,"_",startmajor+length_majorgroup-1,
          startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2, "\n")

    } # end length_majorgroup loop   1.5Mbp each run
    # for each 20Mbp loop write out all small windows 5Kbp and then averages over the 25Kbp set

    if(plotflag){
      if(atflag){
        grDevices::pdf(base::file.path(outpathpdf, base::paste0("ConcatALL_",fname,shortname,"_",
                                        startval,"_",endval,"_",submajor_nam,"_",length_fftgroup,
                                        "_spar1.pdf")))

        base::sink(file=base::file.path(outpathtxt, base::paste0("ConcatALL_",fname,shortname,"_",
                                              startval,"_",endval,"_",submajor_nam,"_",length_fftgroup,
                                              "_spar1.txt")))    #total number of bp is 295045
        #plot concatenated data for major group (1.5Mbp)
        main<- base::paste0("Concatenate ",shortname,"_",
                     startval,"_",endval,"\n",
                     name_majorgroup,submajor_nam,"_",length_fftgroup)
      } else {
        grDevices::pdf(base::file.path(outpathpdf, base::paste0("ConcatALL_CG_",fname,shortname,"_",
                                        startval,"_",endval,"_",submajor_nam,"_",length_fftgroup,
                                        "_spar1.pdf")))

        base::sink(file=base::file.path(outpathtxt, base::paste0("ConcatALL_CG_",fname,shortname,"_",
                                              startval,"_",endval,"_",submajor_nam,"_",length_fftgroup,
                                              "_spar1.txt")))    #total number of bp is 295045
        #plot concatenated data for major group (1.5Mbp)
        main<- base::paste0("Concatenate CG",shortname,"_",
                     startval,"_",endval,"\n",
                     name_majorgroup,submajor_nam,"_",length_fftgroup)
      }

      for (subval in list_plot){     # base::nrow(concat_ALL_tot)
        imagesub(concat_ALL_tot,imageave=NULL,subval=subval,
                 main= base::paste("spar=1",main),fact=8)
      }

      base::sink()
      grDevices::dev.off()
    }
    if(atflag){
      ofile<-base::file.path(outpathspectra, base::paste0("concat_ALL_tot_",shortname,"_",
                                             startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1_Table.txt"))
      base::cat("\n writing to ", ofile,"\n")

    } else {
      ofile<-base::file.path(outpathspectra, base::paste0("concat_ALL_tot_CG_",shortname,"_",
                                             startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1_Table.txt"))
      base::cat("\n writing to ", ofile,"\n")

    }
    utils::write.table(x=concat_ALL_tot, file=ofile, append = FALSE, sep = " ", dec = ".",
                row.names = TRUE, col.names = TRUE)
    All_list<-base::c(All_list,ofile)

    base::print(names(major_concatenationtotlist))
    if(atflag){
      ofile<-base::file.path(outpathfractal, base::paste0("Frac_ALL_",shortname,"_",
                                             startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1_Table.txt"))

    } else {
      ofile<-base::file.path(outpathfractal, base::paste0("Frac_ALL_CG",shortname,"_",
                                             startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1_Table.txt"))
    }
    s1<-base::which(base::as.numeric(names(fractal_dimlist))>=startval)[1]
    utils::write.table(x=fractal_dimlist[s1:base::length(fractal_dimlist)], file=ofile, append = FALSE, sep = " ", dec = ".",
                row.names = TRUE, col.names = TRUE)
    if(atflag){
      ofile<-base::file.path(outpathfractal, base::paste0("Frac_ALL_",shortname,"_",
                                             startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1.pdf"))
      main<- base::paste0("Concatenate ",shortname,"_",
                   startval,"_",endval,"\n",
                   name_majorgroup,submajor_nam,"_",length_fftgroup)
    } else {
      ofile<-base::file.path(outpathfractal, base::paste0("Frac_ALL_CG",shortname,"_",
                                             startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1.pdf"))
      main<- base::paste0("Concatenate_CG ",shortname,"_",
                   startval,"_",endval,"\n",
                   name_majorgroup,submajor_nam,"_",length_fftgroup)
    }

    grDevices::pdf(ofile)


    #plot 4.5Mbp DNA walk(code from fracD
    #xlim<-base::c(fCGlim[1],fCGlim[1]+rangeboth-1);ylim<-base::c(fATlim[1],fATlim[1]+rangeboth-1)
    shortseq<-base::seq(1,base::nrow(DNAwalk_extended),20)
    base::plot(DNAwalk_extended[shortseq,"CG"], DNAwalk_extended[shortseq,"AT"],
         main= base::paste(main), type="l",cex.main=0.6,las=2)
    graphics::points(DNAwalk_extended[1,"CG"],DNAwalk_extended[1,"AT"],type="p",pch=83,col="red",cex=2)   #S
    graphics::points(DNAwalk_extended[base::nrow(DNAwalk_extended),"CG"],DNAwalk_extended[base::nrow(DNAwalk_extended),"AT"],
           type="p",pch=69,col="red",cex=2)  #E

    minf<-base::as.numeric(names(fractal_dimlist)[s1])
    maxf<-base::as.numeric(names(fractal_dimlist)[base::length(fractal_dimlist)]);
    numrange= base::round((maxf-minf)/(15*length_fftgroup))  #length_fftgroup
    run_sum_Fractal(fracdim_list=fractal_dimlist[s1:base::length(fractal_dimlist)], numrange= numrange,pflag=TRUE,main=main)

    #plot 4.5Mbp spectra

    for (subval in list_plot){     # base::nrow(concat_ALL_tot)
      imagesub(concat_ALL_tot,imageave=NULL,subval=subval,
               main= base::paste("spar=1",main),fact=8)
    }
    grDevices::dev.off()

    if(atflag){
      DNAwalk_name<-base::file.path(outpathDNAWALK, base::paste0("FULLDNA_WALK_",shortname,"_",
                                                    startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                                    "_spar1.txt"))

    } else {
      DNAwalk_name<-base::file.path(outpathDNAWALK, base::paste0("FULLDNA_WALK_CG_",shortname,"_",
                                                    startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                                    "_spar1.txt"))

    }
    utils::write.table(x=DNAwalk_extended, file=DNAwalk_name, append = FALSE, sep = " ", dec = ".",
                row.names = TRUE, col.names = TRUE)
    fullwalklist<-base::c(fullwalklist,DNAwalk_name)


    #DNAwalk_FULL will be full walk every 50bp
    shortseq<-base::seq(1,base::nrow(DNAwalk_extended),samplesize)  #changed Feb 6 2023
    if(!base::is.null(DNAwalk_FULL)){
      DNAwalk_FULL<- base::rbind(DNAwalk_FULL,
                          base::cbind(DNAwalk_FULL[base::nrow(DNAwalk_FULL),"AT"]+DNAwalk_extended[shortseq,"AT"],
                                DNAwalk_FULL[base::nrow(DNAwalk_FULL),"CG"]+DNAwalk_extended[shortseq,"CG"]))
    } else {
      DNAwalk_FULL<-DNAwalk_extended[shortseq,]
    }

  }#end majorlen 4.5Mbp loop


  if(atflag){
    main<- base::paste0("Concatenate ",shortname,"_",
                 base::names(fractal_dimlist)[1],"_",names(fractal_dimlist)[base::length(fractal_dimlist)],"\n",
                 name_majorgroup,submajor_nam,"_",length_fftgroup)

    ofile<-base::file.path(outpathfractal, base::paste0("Fractal_ALL_",shortname,"_",
                                           submajor_nam,"_",length_fftgroup,
                                           "_spar1_Table.txt"))
  } else {
    main<- base::paste0("Concatenate_CG ",shortname,"_",
                 base::names(fractal_dimlist)[1],"_",names(fractal_dimlist)[base::length(fractal_dimlist)],"\n",
                 name_majorgroup,submajor_nam,"_",length_fftgroup)

    ofile<-base::file.path(outpathfractal, base::paste0("Fractal_ALL_CG_",shortname,"_",
                                           submajor_nam,"_",length_fftgroup,
                                           "_spar1_Table.txt"))
  }
  base::cat("\n runchromosome: writing ofile Fractal txt file", ofile)
  utils::write.table(x=fractal_dimlist[1:base::length(fractal_dimlist)], file=ofile, append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  if(atflag){
    ofile<-base::file.path(outpathfractal, base::paste0("Fractal_ALL_",shortname,"_",
                                           submajor_nam,"_",length_fftgroup,
                                           "_spar1.pdf"))
  } else {
    ofile<-base::file.path(outpathfractal, base::paste0("Fractal_ALL_CG",shortname,"_",
                                           submajor_nam,"_",length_fftgroup,
                                           "_spar1.pdf"))
  }

  base::cat("\n runchromosome: writing ofile Fractal pdf file", ofile)

  grDevices::pdf(ofile)


  base::plot(DNAwalk_FULL[,"CG"], DNAwalk_FULL[,"AT"],
       main= base::paste(main,samplesize), type="l",cex.main=0.6,las=2)
  graphics::points(DNAwalk_FULL[1,"CG"],DNAwalk_FULL[1,"AT"],type="p",pch=83,col="red",cex=2)   #S
  graphics::points(DNAwalk_FULL[base::nrow(DNAwalk_FULL),"CG"],DNAwalk_FULL[base::nrow(DNAwalk_FULL),"AT"],
         type="p",pch=69,col="red",cex=2)  #E



  minf<-base::as.numeric(names(fractal_dimlist)[1])
  maxf<-base::as.numeric(names(fractal_dimlist)[base::length(fractal_dimlist)]);
  numrange= base::round((maxf-minf)/(15*length_fftgroup))  #length_fftgroup
  run_sum_Fractal(fracdim_list=fractal_dimlist, numrange= numrange,pflag=TRUE,main=main)
  grDevices::dev.off()

  base::cat("\n runchromosome: writing ofile Fractal bed file", ofile)
  if(atflag){
    ofile1<-base::file.path(outpathfractal.bed, base::paste0("Fractal_ALL_",shortname,"_",
                                                submajor_nam,"_",length_fftgroup,
                                                ".bedGraph"))
  } else {
    ofile1<-base::file.path(outpathfractal.bed, base::paste0("Fractal_ALL_CG",shortname,"_",
                                                submajor_nam,"_",length_fftgroup,
                                                ".bedGraph"))
  }
  base::sink(ofile1)
  base::cat("track type=bedGraph")
  for(j in 1:base::length(fractal_dimlist)){
    delta<-base::as.numeric(names(fractal_dimlist)[2])-base::as.numeric(names(fractal_dimlist)[1])  #base::as.numeric(names(fractal)[2])-base::as.numeric(names(fractal)[1])
    base::cat( base::paste0("\n",nam,"\t",base::as.numeric(names(fractal_dimlist)[j]),"\t",
               base::as.numeric(names(fractal_dimlist)[j])+delta,"\t",
               fractal_dimlist[j]) )
  }
  base::sink()
  base::cat("\n runchromosome: returning twolists","\n")
  twolists<-base::list(All_list=All_list,All_listAve=All_listAve,fullwalklist=fullwalklist,fractal_dimlist,fractal_dimlist)
  base::return(twolists)
}



run_chromosomelist<-function(nam=nam, fname=fname, inpath=inpath, outpath=outpath,
                             startgroup=1,endgroup=NULL,
                             majorlen=5000000,
                             length_majorgroup=1000000, name_majorgroup="1Mbp",
                             length_submajor=500000, name_submajorgroup="500Kbp",
                             length_fftgroup=5000,   name_fftgroup="5Kbp",
                             length_minor= 25000, submajor_nam="25Kb",pflag=FALSE,
                             plotflag=TRUE, writeflag=TRUE){

  genname1<-"full"                        # 3.1, 2.2
  Shortnames<-(base::unlist(base::strsplit(nam," " )))
  if(base::length(Shortnames)>7)Shortnames<-Shortnames[1:7]
  shortname<-NULL; for (sh in Shortnames)shortname<- base::paste0(shortname,"_",sh)
  in_name<- base::paste0(inpath,nam,".fasta")    #
  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")

  main= base::paste0(genname1,"\n",shortname)
  dn1<-dna_1_to_wax(dna_1)
  if(base::is.null(endgroup))full_length<-base::length(dn1) else full_length<-endgroup

  #run segments of 7.5Mbp each    1.5Mbp groupings   500 Kbp per group 25Kbp ave 5000 bp window
  start_seq<-base::seq(startgroup,full_length, by=majorlen) #starttgroup<-45000001
  All_list<-NULL;All_listAve<-NULL
  for(startval in start_seq){
    endval<-startval+majorlen -1
    dn<-startval:endval
    dna_1<-dn1<-NULL
    base::gc()
    #put below in 20Mb  loop
    majorlen<- base::length(dn)
    wholelist<-base::seq(1,base::length(dn), by=length_majorgroup )
    #run through whole chromosome
    for (startmajor in wholelist){
      #startmajor<-15000001
      endmajor<-startmajor+length_majorgroup-1

      majorlist<- base::seq(startmajor,startmajor+(length_majorgroup-length_submajor), by=length_submajor)
      concat_majorlist_tot<-NULL
      concat_majorlist_ave<-NULL

    }
    ofile<- base::paste0(outpath,"concat_ALL_tot_",fname,shortname,"_",
                  startval,"_",startval+majorlen-1,"_",submajor_nam,"_",length_fftgroup,
                  "_spar1_Table.txt")
    All_list<-base::c(All_list,ofile)

    ofile<- base::paste0(outpath,"concat_ALL_ave_",fname,shortname,"_",
                  startval,"_",startval+majorlen-1,"_",submajor_nam,"_",length_fftgroup,
                  "_spar1_Table.txt")
    All_listAve<-base::c(All_listAve,ofile)

  }#end 7.5Mbp loop

  twolists<-base::list(All_list=All_list,All_listAve=All_listAve)
  base::return(twolists)
}

#' run_chromosomelistNEW
#'
#' rebuilds twolists if already run
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_chromosomelistNEW<-function(nam=nam, fname=fname, inpath=inpath, outpath=outpath,
                                startgroup=1,endgroup=NULL,
                                majorlen=5000000,
                                length_majorgroup=1000000, name_majorgroup="1500Kbp",
                                length_submajor=500000, name_submajorgroup="500Kbp",
                                length_fftgroup=5000,   name_fftgroup="5Kbp",
                                length_minor= 25000, submajor_nam="25Kb",full_length=NULL){

  chromnam<-nam
  splitname<-base::unlist(stringr::str_split(nam,"_"))
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]
  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)
  base::dir.create(base::file.path(outpath, nam1,chr), showWarnings = FALSE)
  outpathspectra<-base::file.path(outpath, nam1,chr,"spectra_Table.txt")
  base::dir.create(outpathspectra, showWarnings = FALSE)
  outpathwalk<-base::file.path(outpath, nam1,chr,"DNAWALK")
  base::dir.create(outpathwalk, showWarnings = FALSE)
  genname1<-"full"
  Shortnames<-(base::unlist(base::strsplit(nam," " )))
  if(base::length(Shortnames)>7)Shortnames<-Shortnames[1:7]
  shortname<-NULL; for (sh in Shortnames)shortname<- base::paste0(shortname,"_",sh)
  in_name<- base::paste0(inpath,nam,".fasta")    #
  if(base::is.null(full_length)){
    dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")
    main= base::paste0(genname1,"\n",shortname)
    dn1<-dna_1_to_wax(dna_1)
    lendna<-base::length(dn1)
  } else lendna<-full_length
  if(base::is.null(endgroup))full_length<-lendna else full_length<-endgroup
  #run segments of 7.5Mbp each    1.5Mbp groupings   500 Kbp per group 25Kbp ave 5000 bp window
  start_seq<-base::seq(startgroup,(full_length), by=majorlen) #starttgroup<-45000001
  All_list<-NULL;All_listAve<-NULL; fullwalklist<-NULL
  for(startval in start_seq){
    endval<-startval+majorlen -1
    dn<-startval:endval
    dna_1<-dn1<-NULL
    base::gc()
    #put below in 20Mb  loop
    majorlen1<- base::length(dn)
    wholelist<-base::seq(1,base::length(dn), by=length_majorgroup )
    #run through whole chromosome
    for (startmajor in wholelist){
      #startmajor<-15000001
      endmajor<-startmajor+length_majorgroup-1

      majorlist<- base::seq(startmajor,startmajor+(length_majorgroup-length_submajor), by=length_submajor)
      concat_majorlist_tot<-NULL
      concat_majorlist_ave<-NULL
      #
      for(val1 in majorlist){     # startval<-1      endval<-startval+majorlen -1
        #lengthsubmajor is incremented by 500 kbp
        val2<-val1+length_submajor-1
        if((startval+val1-1)<=endval){
          submajor_num<-base::floor((length_submajor)/length_minor)
          submajorlist<- base::seq(1,(length_submajor), by=length_minor)    #increments of 25Kbp each
          stval<-val1
        }
      } # end major loop   1.5Mbp each run
      # for each 20Mbp loop write out all small windows 5Kbp and then averages over the 25Kbp set

    }
    # each 4.5Mbp
    # dna walk list for every 25Kbp with sub headings of 500 kbp
    sname<-base::file.path(outpathwalk, base::paste0("FULLDNA_WALK_",shortname,"_",
                                        startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                        "_spar1.txt"))

    fullwalklist<-base::c(fullwalklist,sname)
    ofile<-base::file.path(outpathspectra, base::paste0("concat_ALL_tot_",shortname,"_",
                                           startval,"_",startval+majorlen-1,"_",submajor_nam,"_",length_fftgroup,
                                           "_spar1_Table.txt"))
    base::cat("\n",base::file.path(outpathspectra, base::paste0("concat_ALL_tot_",shortname,"_",
                                             startval,"_",startval+majorlen-1,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1_Table.txt")))

    base::cat("\n",startval+majorlen-1,majorlen,"\n")

    All_list<-base::c(All_list,ofile)

    ofile<-base::file.path(outpathspectra, base::paste0("concat_ALL_ave_",shortname,"_",
                                           startval,"_",startval+majorlen-1,"_",submajor_nam,"_",length_fftgroup,
                                           "_spar1_Table.txt"))
    All_listAve<-base::c(All_listAve,ofile)

  }#end 4.5Mbp loop

  twolists<-base::list(All_list=All_list,All_listAve=All_listAve,fullwalklist=fullwalklist)
  base::return(twolists)
}

#' read_walk_list
#'
#' # Reads list of DNAwalk files
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
read_walk_list<-function(fullwalklist){
  #read the files in and concatenate them togather
  All_walk<-NULL;oldwalk<-base::matrix(base::c(0,0),nrow=1,ncol=2)
  for(ofile in fullwalklist){
    base::cat("\nread_walk_list:  opening ",ofile,"\n")
    walk<-utils::read.table( file=ofile,  sep = " ", dec = ".",header = TRUE)
    walk[,1]<-walk[,1]+oldwalk[,1]
    walk[,2]<-walk[,2]+oldwalk[,2]
    walkbp<- base::as.numeric(base::colnames(walk))

    All_walk<- base::rbind(All_walk,walk)
    oldwalk<-walk[base::nrow(walk),]

  }
  base::return(All_walk)
}

#' read_All_list
#'
#' Reads All list to make All Spec
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
read_All_list<-function(All_list,maxrow=100000){
  #read the files in and concatenate them togather
  All_spec<-NULL
  #maxrow<-721 #maxrow<-688
  for(ofile in All_list){
    base::cat("\n opening ",ofile,"\n")
    #for(ofile in All_list[1:2]){
    concat_ALL_tot<-utils::read.table( file=ofile,  sep = " ", dec = ".",header = TRUE, nrows=maxrow)
    #concat_ALL_tot<-utils::read.table( file=ofile,  sep = " ", dec = ".",header = TRUE, check.names=TRUE)
    #, col.names = 1   row.names = 1
    base::colnames(concat_ALL_tot)<-base::gsub("X", "", base::colnames(concat_ALL_tot))
    #base::nrow(concat_ALL_tot) ; base::ncol(concat_ALL_tot)
    All_spec<-base::cbind(All_spec,base::as.matrix(concat_ALL_tot) )
    #imagenan(base::log(concat_ALL_tot))
    #base::plot(1: base::ncol(concat_ALL_tot),concat_ALL_tot["1/189.87",],ylim=base::c(0,base::max(concat_ALL_tot["1/189.87",])),type="b")

  }

  base::nrow(All_spec); base::ncol(All_spec)
  #jpeg(filename = "Rplot%03d.jpeg", width = base::ncol(Al
  base::return(All_spec)
}

#' run_largeimage
#'
#' Reads All Spec into the full chromosome fourier spectrum
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_largeimage<-function(All_spec,nam=nam, fname=fname, inpath=inpath, outpath=outpath,atflag=TRUE,
                         rangeseq=NULL,pngflag=TRUE,samplesize=1,
                         r0=base::c(35,2000),r1=base::c(15,35),r2=base::c(2,8),#r3=base::c(5,35),r4=base::c(5,2000),
                         inmain=""){

  splitname<-base::unlist(stringr::str_split(nam,"_"))   #fname<-"Cannabis"
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]
  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)  #creat directory of individual with chromosomes
  base::dir.create(base::file.path(outpath, nam1,chr), showWarnings = FALSE)   # create Chromosome directory
  outpathimage<-base::file.path(outpath, nam1,chr,"largeimages.png")

  base::dir.create( outpathimage, showWarnings = FALSE) #create image dir for spectra
  if(base::is.null(rangeseq)){
    s1<-base::colnames(All_spec)[1]; e1<-base::colnames(All_spec)[ base::ncol(All_spec)]
    rangeseq<-base::c(s1,e1)
    base::cat("\ns1 e1 :",s1,e1,"\n")
    rangeseq<-base::c(s1,e1)

  } else {
    s1<-rangeseq[1]; e1<-rangeseq[2]
    base::cat("\nnot NULL rangeseq: s1 e1 :",s1,e1,"\n")
  }
  if(atflag){
    ofi<-base::file.path(outpathimage, base::paste0("All_spec",samplesize,"_",nam,"_"))
    inmain<- base::paste(inmain,"AT")
  } else {
    ofi<-base::file.path(outpathimage, base::paste0("All_spec_CG",samplesize,"_",nam,"_"))
    inmain<- base::paste(inmain,"CG")
  }
  base::cat("\nofi:",ofi,"\n")
  if(!base::is.null(r0))largeimagesub_NEW(All_spec,fname, nam, ofi, rangebp=r0,rangeseq1=rangeseq,flaglog=TRUE,main=nam,pngflag = pngflag)
  if(!base::is.null(r1))largeimagesub_NEW(All_spec,fname, nam, ofi, rangebp=r1,rangeseq1=rangeseq,flaglog=TRUE,main=nam,pngflag = pngflag)
  if(!base::is.null(r2))largeimagesub_NEW(All_spec,fname, nam, ofi, rangebp=r2,rangeseq1=rangeseq,flaglog=TRUE,main=nam,pngflag = pngflag)

}


find_lines_in_large_genome<-function(tdna,rng=NULL){
  if(base::is.null(rng)) rng<- base::c(1:nchar(tdna)) else {
    l1<- base::floor(rng[1]/(81))
    l2<- base::ceiling(rng[2]/(81))
  }
  base::return(base::c(l1,l2))
}

Nested_Palindrome<-function(tdna,iterate=1,N=1, finversion=TRUE){
  if(iterate<=N){
    tdna<-true_palindrome(tdna,finversion=finversion) #inversion
    Itdna<-DNA_Palindrome(tdna) #inversion
    tdna<- base::paste0(tdna,Itdna,collapse="")
    iterate<-iterate+1

    tdna<-Nested_Palindrome(tdna,iterate=iterate,N=N)
  }

  base::return(tdna)
}


#' run_plot_chromosome
#'
#' Creates fourier spectrum of chromosome with plotting
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_plot_chromosome <-function(nam=nam, fname=fname, inpath=inpath, outpath=outpath,
                               atflag=TRUE, AT_flag=TRUE, pflag=FALSE, plotflag=FALSE, writeflag=FALSE){

  base::dir.create(outpath, showWarnings = FALSE)
  in_name<- base::paste0(inpath,nam,".fasta")    #
  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")
  dn1<-dna_1_to_wax(dna_1)
  full_length<-base::length(dn1)

  #this runs CG content
  winbox<-1*1e6
  walklist<- CG_AT_content(dn1,nam,winbox=winbox)

  dn1<-NULL; dna_1<-NULL;base::gc()

  #requires fft window to be multiple of 10 for dnawalk in run_chloroplast
  max_size<-2.0*1e6
  samplesize<-base::round(full_length/max_size)     #samplesize for long spectrum fft windows 150/2=75 100/2=50  50/2=25 20/2=10
  samplesize<-base::floor(samplesize/100)*100
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /50)*50   #full_length<-194*1e6 full_length<-210*1e6
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /20)*20
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /10)*10   #full_length<-6*1e6 full_length<-3*1e6
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /5)*5
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /2)*2     #full_length<-1.5*1e6 full_length<-3*1e6
  if(samplesize==0) samplesize<-1
  base::cat("\nSamplesize is", samplesize,"\n")

  # builds the files for spectrum - 0.5, 1.5, 4.5Mbp at a time
  twolists<-run_chromosome(nam=nam, fname=fname, inpath=inpath, outpath=outpath, atflag=atflag, AT_flag=AT_flag, majorlen=5000000, length_majorgroup=1000000, name_majorgroup="1Mbp",plotflag = plotflag, writeflag=writeflag, samplesize=samplesize)
  All_list<-twolists$All_list
  fullwalk_list<-twolists$fullwalklist

  # sample every 5th
  All_spec<-All_spec[,base::seq(1,full_length, by =5)]

  whival<-(base::which(base::as.numeric(base::colnames(All_spec))>full_length+2500))
  if(base::length(whival)>0){
    All_spec[1:10,(whival)]
    All_spec<-All_spec[,-(whival)]
  }
  base::gc()

  # create large images
  run_largeimage(All_spec,nam=nam, fname=fname, inpath=inpath, outpath=outpath, atflag=atflag, inmain="Individual 5Kbp window")
  base::gc()

  diff<-0
  run_largeimage(All_spec, nam=nam, fname=fname, inpath=inpath, outpath=outpath, atflag=atflag, rangeseq=NULL,pngflag=TRUE,samplesize=1,
                 r0=base::c(5,200),r1=NULL,r2=NULL,inmain="Individual 5Kbp window")
  base::gc()

}

#' run_chromosome_parallel loop
#'
#' Loops over 5Mbp in parallel
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
#' @export
run_chromosome_loop <- function(x, nam, fname, atflag, AT_flag, startgroup, majorlen, length_majorgroup, name_majorgroup,
                                length_submajor, name_submajorgroup,
                                length_fftgroup,   name_fftgroup,
                                length_minor, submajor_nam, pflag,
                                splitname, chr, outpathtxt,outpathpdf,
                                outpathfractal, outpathspectra,outpathfractal.bed,
                                outpathDNAWALK, genname1,Shortnames,shortname,
                                in_name, dna_1,main,dn1, list_plot, full_length,start_seq,
                                All_list, All_listAve, fullwalklist, fractal_dimlist, DNAwalk_FULL,
                                plotflag, writeflag, samplesize){

  #  for(startval in start_seq){
  # https://www.r-bloggers.com/2017/10/running-r-code-in-parallel/
  library(RepeatOBserverV1)

  startval <<- start_seq[x]
  DNAwalk_extended<-NULL
  endval<-startval+majorlen -1

  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")
  dn1<-dna_1_to_wax(dna_1)

  dn<-dn1[startval:endval]
  dna_1<-dn1<-NULL
  base::gc()
  #put below in 20Mb  loop
  majorlen<- base::length(dn)

  #method to make list of data for major groups etc
  wholelist<-base::seq(1,base::length(dn), by=length_majorgroup )

  #run through whole chromosome
  major_concatenationtotlist<-NULL
  major_concatenationavelist<-NULL

  concat_ALL_tot<-NULL
  concat_ALL_ave<-NULL
  for (startmajor in wholelist){
    #startmajor<-15000001
    endmajor<-startmajor+length_majorgroup-1

    majorlist<- base::seq(startmajor,startmajor+(length_majorgroup-length_submajor), by=length_submajor)
    concat_majorlist_tot<-NULL
    concat_majorlist_ave<-NULL
    for(val1 in majorlist){     # startval<-1      endval<-startval+majorlen -1  length_submajor=500000

      val2<-val1+length_submajor-1
      if(val2<=endval){
        tdna<-dn[val1:val2]

        submajor_num<-base::floor((length_submajor)/length_minor)

        submajorlist<- base::seq(1,(length_submajor), by=length_minor)
        stval<-val1
        sname<-base::file.path(outpathDNAWALK, base::paste0("DNA_WALK_",fname,shortname,"_",startval+val1-1,"_",startval+val1+length_submajor-2,"_",
                                               submajor_nam,"_",length_fftgroup,"_spar1"))
        if(plotflag){
          if(atflag){
            grDevices::pdf(base::file.path(outpathpdf, base::paste0(fname,shortname,"_",
                                            startval+val1-1,"_",startval+val1+length_submajor-2,"_",submajor_nam,"_",length_fftgroup,
                                            "_spar1.pdf")))

            ofi<-base::file.path(outpathtxt, base::paste0(fname,shortname,"_",
                                             startval+val1-1,"_",startval+val1+length_submajor-2,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1.txt"))
          } else{
            grDevices::pdf(base::file.path(outpathpdf, base::paste0(fname,shortname,"_",
                                            startval+val1-1,"_",startval+val1+length_submajor-2,"_",submajor_nam,"_",length_fftgroup,
                                            "_spar1_CG.pdf")))
            ofi<-base::file.path(outpathtxt, base::paste0(fname,shortname,"_",
                                             startval+val1-1,"_",startval+val1+length_submajor-2,"_",submajor_nam,"_",length_fftgroup,
                                             "_spar1_CG.txt"))
          }
          base::cat("\nrun_chromosome: output going to",ofi,"\n")
          base::sink(file=ofi)    #total number of bp is 295045
        }
        imagetot<-NULL;imageave<-NULL
        for (j in submajorlist){       #length_minor= 25000
          #j<-470289;   inc<-49503   jlist<-base::c(1,25001)  submajorlist<-base::c(1);stval<-0
          inc<-length_minor
          fullname<- base::paste0(sname,"_",(startval+stval+j-2),"_",(startval+stval+j+inc-3) )
          if((j+inc-1)<=(length_submajor)){
            # base::cat("\nrun_chromosome:running ",
            #      base::paste0("rng=base::c(",j,"_" ,j+inc-1,")"," ",(startval+stval+j-2),"_" ,(startval+stval+j+inc-3) ,")"))
            genname1<- base::paste0("rng=base::c(",j,"_" ,j+inc-1,")"," ",(startval+stval+j-2),"_" ,(startval+stval+j+inc-3) ,")")
            main= base::paste(genname1,"\n",shortname)
            spectinfo<-walk_and_plot(tdna=tdna,main=main,rng=base::c(j,j+inc-1),seqlen=length_fftgroup,
                                     fullname=fullname, spar=1,startval=(startval+stval+j-2),walkflag=FALSE,waxflag=TRUE,
                                     pflag=FALSE,plotflag=plotflag,writeflag=writeflag,atflag=atflag, AT_flag=AT_flag)
          } else {
            # end of region values
            # base::cat("\nrun_chromosome: running ",
            #      base::paste0("rng=base::c(",j,"_" ,length_submajor,")"," ",(startval+stval+j-2),"_" ,(startval+stval+length_submajor) ,")"))
            genname1<- base::paste0("rng=base::c(",j,"_" ,length_submajor,")"," ",(startval+stval+j-2),"_" ,(startval+stval+length_submajor) ,")")
            main= base::paste(genname1,"\n",shortname)
            spectinfo<-walk_and_plot(tdna,main=main,rng=base::c(j,j+inc-1),seqlen=5000,
                                     fullname=fullname, spar=1,startval=(startval+stval+j-1),walkflag=FALSE,waxflag=TRUE,
                                     plotflag=plotflag,writeflag=writeflag, atflag=atflag, AT_flag=AT_flag)

          }
          peaks<-spectinfo[[1]];
          spectimage<-spectinfo[[2]];
          spectmean<-spectinfo[[3]]
          DNAwalk<-spectinfo[[4]]
          spectwalklist<-spectinfo[[5]]
          fractalD<-spectinfo[[6]]
          fractal_dimlist<-base::c(fractal_dimlist,fractalD)

          if(writeflag){
            base::cat("\n\n\n ", base::paste0("Final peaks\nrng=base::c(",j,"_" ,j+inc-1,")"," ",
                                 (startval+stval+j),"_" ,(startval+stval+j+inc-1) ,")"),"\n")
            if(base::nrow(peaks)>=20) base::print(peaks[1:20,]) else base::print(peaks[1:base::nrow(peaks),])
          }
          imagetot<-base::cbind(imagetot,spectimage)

          if(writeflag){
            base::cat("\nstart dwalk calc: DNAwalk\n")
            base::print(DNAwalk[1:10,])
          }
          if(!base::is.null(DNAwalk_extended)){
            if(atflag){
              DNAwalk_extended<- base::rbind(DNAwalk_extended,
                                      base::cbind(DNAwalk_extended[base::nrow(DNAwalk_extended),"AT"]+DNAwalk[,"AT"],
                                            DNAwalk_extended[base::nrow(DNAwalk_extended),"CG"]+DNAwalk[,"CG"]))
            } else {
              DNAwalk_extended<- base::rbind(DNAwalk_extended,
                                      base::cbind(DNAwalk_extended[base::nrow(DNAwalk_extended),"CG"]+DNAwalk[,"CG"],
                                            DNAwalk_extended[base::nrow(DNAwalk_extended),"AT"]+DNAwalk[,"AT"]))
            }
          } else {
            DNAwalk_extended<-DNAwalk
          }
        } #end of length_minor loop  (25Kbp)

        #plot images for 0.5 Mbp intervals
        if(plotflag){
          if(base::max(imagetot,na.rm=TRUE)>base::min(imagetot,na.rm=TRUE)){
            for (subval in list_plot){  #subval<-1
              main= base::paste0(startval+val1,"_",startval+val1+length_submajor,"\n",shortname)
              imagesub(imagetot,imageave,subval=subval,
                       main= base::paste("spar=1",main),fact=8)
            }
          }
        }

        if(plotflag){
          base::sink()
          grDevices::dev.off()
        }
        #concatenate the 500Kbp images/data
        if(base::max(imagetot,na.rm=TRUE)>base::min(imagetot,na.rm=TRUE)){
          concat_majorlist_tot<-base::cbind(concat_majorlist_tot,imagetot)
        }
      }

    } #end of length_submajor loop  (0.5Mbp)
    if(plotflag){
      if(!base::is.null(base::nrow(concat_majorlist_tot))){
        if(atflag){
          grDevices::pdf(base::file.path(outpathpdf, base::paste0("Concatenate_",fname,shortname,"_",
                                          startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                          "_spar1.pdf")))
          base::sink(file=base::file.path(outpathtxt, base::paste0("Concatenate_",fname,shortname,"_",
                                                startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                                "_spar1.txt")))    #total number of bp is 295045
          #plot concatenated data for major group (1.5Mbp)
          main<- base::paste0("Concatenate ",shortname,"_",
                       startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"\n",
                       name_majorgroup,submajor_nam,"_",length_fftgroup)
        } else {
          grDevices::pdf(base::file.path(outpathpdf, base::paste0("Concatenate_CG_",fname,shortname,"_",
                                          startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                          "_spar1.pdf")))
          base::sink(file=base::file.path(outpathtxt, base::paste0("Concatenate_CG_",fname,shortname,"_",
                                                startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                                "_spar1.txt")))    #total number of bp is 295045
          #plot concatenated data for major group (1.5Mbp)
          main<- base::paste0("Concatenate CG",shortname,"_",
                       startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"\n",
                       name_majorgroup,submajor_nam,"_",length_fftgroup)
        }

        base::cat("\nrun_chromosome: entering imagesub loop for Concatenate spectra","\n")
        for (subval in list_plot){    # base::nrow(concat_majorlist_tot)
          imagesub(concat_majorlist_tot,imageave=NULL,subval=subval,
                   main= base::paste("spar=1",main),fact=8)                         #? ,sval=startval+startmajor
        }
        base::cat("\nrun_chromosome: exiting imagesub loop for Concatenate spectra","\n")
        base::sink()
        grDevices::dev.off()
        if(atflag){
          ofile<-base::file.path(outpathfractal, base::paste0("FracCon_",fname,shortname,"_",
                                                 startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                                 "_spar1.txt"))

          #fractal_dimlist for1.5Mbp
          s1<-base::which(base::as.numeric(names(fractal_dimlist))>=startval+startmajor-1)[1]
          base::cat("\nrun_chromosome: Writing to",ofile,"\n Starting at s1",s1,"\n")
          utils::write.table(x=fractal_dimlist[s1:base::length(fractal_dimlist)], file=ofile, append = FALSE, sep = " ", dec = ".",
                      row.names = TRUE, col.names = TRUE)


          ofile<-base::file.path(outpathfractal, base::paste0("FraCon_",fname,shortname,"_",
                                                 startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2,"_",submajor_nam,"_",length_fftgroup,
                                                 "_spar1.pdf"))
          base::cat("\nrun_chromosome: Writing to",ofile)
          grDevices::pdf(ofile)
          minf<-base::as.numeric(names(fractal_dimlist)[s1])
          maxf<-base::as.numeric(names(fractal_dimlist)[base::length(fractal_dimlist)]);
          numrange= base::round((maxf-minf)/(15*length_fftgroup))  # length_fftgroup
          base::cat("\nnumrange minf maxf set to ",numrange,minf,maxf,"Enter run_sum_Fractal\n")
          run_sum_Fractal(fracdim_list=fractal_dimlist[s1:base::length(fractal_dimlist)], numrange= numrange,pflag=TRUE,main=main)
          grDevices::dev.off()
        }

      }   #end test for null concat_majorlist_tot
    }
    list2<-base::list(concat_majorlist_tot)
    major_concatenationtotlist<- base::c(major_concatenationtotlist, list2=list2)
    base::names(major_concatenationtotlist)[base::length(major_concatenationtotlist)]<- base::paste0(fname,
                                                                                  shortname,"_",
                                                                                  startmajor,"_",startmajor+length_majorgroup-1,"_",submajor_nam,"_",length_fftgroup, "_spar1")

    #plot all together
    concat_ALL_tot<-base::cbind(concat_ALL_tot,concat_majorlist_tot)
    #concat_ALL_ave<-base::cbind(concat_ALL_ave,concat_majorlist_ave)
    base::cat("\n run_chromosome: Ending :",startmajor,"_",startmajor+length_majorgroup-1,
        startval+startmajor-1,"_",startval+startmajor+length_majorgroup-2, "\n")

  } # end length_majorgroup loop   1.5Mbp each run
  # for each 20Mbp loop write out all small windows 5Kbp and then averages over the 25Kbp set

  if(plotflag){
    if(atflag){
      grDevices::pdf(base::file.path(outpathpdf, base::paste0("ConcatALL_",fname,shortname,"_",
                                      startval,"_",endval,"_",submajor_nam,"_",length_fftgroup,
                                      "_spar1.pdf")))

      base::sink(file=base::file.path(outpathtxt, base::paste0("ConcatALL_",fname,shortname,"_",
                                            startval,"_",endval,"_",submajor_nam,"_",length_fftgroup,
                                            "_spar1.txt")))    #total number of bp is 295045
      #plot concatenated data for major group (1.5Mbp)
      main<- base::paste0("Concatenate ",shortname,"_",
                   startval,"_",endval,"\n",
                   name_majorgroup,submajor_nam,"_",length_fftgroup)
    } else {
      grDevices::pdf(base::file.path(outpathpdf, base::paste0("ConcatALL_CG_",fname,shortname,"_",
                                      startval,"_",endval,"_",submajor_nam,"_",length_fftgroup,
                                      "_spar1.pdf")))

      base::sink(file=base::file.path(outpathtxt, base::paste0("ConcatALL_CG_",fname,shortname,"_",
                                            startval,"_",endval,"_",submajor_nam,"_",length_fftgroup,
                                            "_spar1.txt")))    #total number of bp is 295045
      #plot concatenated data for major group (1.5Mbp)
      main<- base::paste0("Concatenate CG",shortname,"_",
                   startval,"_",endval,"\n",
                   name_majorgroup,submajor_nam,"_",length_fftgroup)
    }

    for (subval in list_plot){     # base::nrow(concat_ALL_tot)
      imagesub(concat_ALL_tot,imageave=NULL,subval=subval,
               main= base::paste("spar=1",main),fact=8)
    }

    base::sink()
    grDevices::dev.off()
  }
  if(atflag){
    ofile<-base::file.path(outpathspectra, base::paste0("concat_ALL_tot_",shortname,"_",
                                           startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                           "_spar1_Table.txt"))
    base::cat("\n writing to ", ofile,"\n")

  } else {
    ofile<-base::file.path(outpathspectra, base::paste0("concat_ALL_tot_CG_",shortname,"_",
                                           startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                           "_spar1_Table.txt"))
    base::cat("\n writing to ", ofile,"\n")

  }
  utils::write.table(x=concat_ALL_tot, file=ofile, append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  All_list<-base::c(All_list,ofile)

  base::print(names(major_concatenationtotlist))
  if(atflag){
    ofile<-base::file.path(outpathfractal, base::paste0("Frac_ALL_",shortname,"_",
                                           startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                           "_spar1_Table.txt"))

  } else {
    ofile<-base::file.path(outpathfractal, base::paste0("Frac_ALL_CG",shortname,"_",
                                           startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                           "_spar1_Table.txt"))
  }
  s1<-base::which(base::as.numeric(names(fractal_dimlist))>=startval)[1]
  utils::write.table(x=fractal_dimlist[s1:base::length(fractal_dimlist)], file=ofile, append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  if(atflag){
    ofile<-base::file.path(outpathfractal, base::paste0("Frac_ALL_",shortname,"_",
                                           startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                           "_spar1.pdf"))
    main<- base::paste0("Concatenate ",shortname,"_",
                 startval,"_",endval,"\n",
                 name_majorgroup,submajor_nam,"_",length_fftgroup)
  } else {
    ofile<-base::file.path(outpathfractal, base::paste0("Frac_ALL_CG",shortname,"_",
                                           startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                           "_spar1.pdf"))
    main<- base::paste0("Concatenate_CG ",shortname,"_",
                 startval,"_",endval,"\n",
                 name_majorgroup,submajor_nam,"_",length_fftgroup)
  }

  grDevices::pdf(ofile)

  shortseq<-base::seq(1,base::nrow(DNAwalk_extended),20)
  base::plot(DNAwalk_extended[shortseq,"CG"], DNAwalk_extended[shortseq,"AT"],
       main= base::paste(main), type="l",cex.main=0.6,las=2)
  graphics::points(DNAwalk_extended[1,"CG"],DNAwalk_extended[1,"AT"],type="p",pch=83,col="red",cex=2)   #S
  graphics::points(DNAwalk_extended[base::nrow(DNAwalk_extended),"CG"],DNAwalk_extended[base::nrow(DNAwalk_extended),"AT"],
         type="p",pch=69,col="red",cex=2)  #E

  minf<-base::as.numeric(names(fractal_dimlist)[s1])
  maxf<-base::as.numeric(names(fractal_dimlist)[base::length(fractal_dimlist)]);
  numrange= base::round((maxf-minf)/(15*length_fftgroup))  #length_fftgroup
  run_sum_Fractal(fracdim_list=fractal_dimlist[s1:base::length(fractal_dimlist)], numrange= numrange,pflag=TRUE,main=main)

  #plot 4.5Mbp spectra


  for (subval in list_plot){     # base::nrow(concat_ALL_tot)
    imagesub(concat_ALL_tot,imageave=NULL,subval=subval,
             main= base::paste("spar=1",main),fact=8)
  }

  grDevices::dev.off()

  if(atflag){
    DNAwalk_name<-base::file.path(outpathDNAWALK, base::paste0("FULLDNA_WALK_",shortname,"_",
                                                  startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                                  "_spar1.txt"))

  } else {
    DNAwalk_name<-base::file.path(outpathDNAWALK, base::paste0("FULLDNA_WALK_CG_",shortname,"_",
                                                  startval,"_",startval+endmajor-1,"_",submajor_nam,"_",length_fftgroup,
                                                  "_spar1.txt"))

  }
  utils::write.table(x=DNAwalk_extended, file=DNAwalk_name, append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
  fullwalklist<-base::c(fullwalklist,DNAwalk_name)


}#end majorlen 4.5Mbp loop


#' run_chromosome_parallel
#'
#' Creates fourier spectrum of chromosome with plotting
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_chromosome_parallel<-function(nam=nam, fname=fname, inpath=inpath, outpath=outpath,
                                  startgroup=1, atflag=TRUE, AT_flag=TRUE, majorlen=5000000,
                                  length_majorgroup=1000000, name_majorgroup="1Mbp",
                                  length_submajor=500000, name_submajorgroup="500Kbp",
                                  length_fftgroup=5000,   name_fftgroup="5Kbp",
                                  length_minor= 25000, submajor_nam="25Kb",pflag=FALSE,
                                  plotflag=TRUE, writeflag=TRUE, x_cpu=x_cpu, samplesize=50){

  splitname<-base::unlist(stringr::str_split(nam,"_"))
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]
  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)
  base::dir.create(base::file.path(outpath, nam1,chr), showWarnings = FALSE)
  outpathtxt<-base::file.path(outpath, nam1,chr,"5000bp_spectra.txt")
  base::dir.create(outpathtxt, showWarnings = FALSE)
  outpathpdf<-base::file.path(outpath, nam1,chr,"5000bp_spectra.pdf")
  base::dir.create(outpathpdf, showWarnings = FALSE)
  outpathfractal<-base::file.path(outpath, nam1,chr,"5000bp_fractal.pdf")
  base::dir.create(outpathfractal, showWarnings = FALSE)
  outpathspectra<-base::file.path(outpath, nam1,chr,"spectra_Table.txt")
  base::dir.create(outpathspectra, showWarnings = FALSE)   # c
  outpathfractal.bed<-base::file.path(outpath, nam1,chr,"5000bpfractalD.bedGraph")
  base::dir.create( outpathfractal.bed, showWarnings = FALSE)
  outpathDNAWALK<-base::file.path(outpath, nam1,chr,"DNAWALK")
  base::dir.create( outpathDNAWALK, showWarnings = FALSE)
  genname1<-"full"                        # 3.1, 2.2
  Shortnames<-(base::unlist(base::strsplit(nam," " )))
  if(base::length(Shortnames)>7)Shortnames<-Shortnames[1:7]
  shortname<-NULL; for (sh in Shortnames)shortname<- base::paste0(shortname,"_",sh)
  in_name<- base::paste0(inpath,nam,".fasta")    #
  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")
  main= base::paste0(genname1,"\n",shortname)
  dn1<-dna_1_to_wax(dna_1)

  list_plot<-base::c(5,10,20,30,1,2)
  full_length<-base::length(dn1)
  #run segments of 7.5Mbp each    1.5Mbp groupings   500 Kbp per group 25Kbp ave 5000 bp window
  start_seq<-base::seq(startgroup,full_length, by=majorlen) #starttgroup<-45000001   majorlen=4500000
  DNAwalk_FULL<-NULL
  All_list<-NULL;All_listAve<-NULL; fullwalklist<-NULL
  fractal_dimlist<-NULL
  ncpu=full_length/5000000
  cl <- parallel::makeCluster(x_cpu)
  results <- parallel::parSapply(cl, base::seq_along(start_seq), run_chromosome_loop, nam=nam,
                       fname=fname, atflag=atflag, AT_flag=AT_flag, startgroup=startgroup, majorlen=majorlen,
                       length_majorgroup=length_majorgroup, name_majorgroup=name_majorgroup,
                       length_submajor=length_submajor, name_submajorgroup=name_submajorgroup,
                       length_fftgroup=length_fftgroup,   name_fftgroup=name_fftgroup,
                       length_minor=length_minor, submajor_nam=submajor_nam, pflag=pflag,
                       splitname=splitname, chr=chr, outpathtxt=outpathtxt,outpathpdf=outpathpdf,
                       outpathfractal=outpathfractal, outpathspectra=outpathspectra,
                       outpathfractal.bed=outpathfractal.bed,
                       outpathDNAWALK=outpathDNAWALK, genname1=genname1,Shortnames=Shortnames,shortname=shortname,
                       in_name=in_name, dna_1=dna_1,main=main,dn1=dn1, list_plot=list_plot,
                       full_length=full_length,start_seq=start_seq,All_list=All_list, All_listAve=All_listAve,
                       fullwalklist=fullwalklist, fractal_dimlist=fractal_dimlist, DNAwalk_FULL=DNAwalk_FULL,
                       plotflag=plotflag, writeflag=writeflag, samplesize=samplesize)

  base::cat("\n runchromosome: returning twolists","\n")
  twolists<-base::list(All_list=All_list,All_listAve=All_listAve,fullwalklist=fullwalklist,fractal_dimlist,fractal_dimlist)
  base::return(twolists)
}


#' run_plot_chromosome_parallel
#'
#' Creates fourier spectrum of chromosome
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_plot_chromosome_parallel <-function(nam=nam, fname=fname, inpath=inpath, AT_flag=TRUE,
                                        outpath=outpath, atflag=TRUE,  pflag=FALSE, plotflag=FALSE,
                                        writeflag=FALSE, x_cpu=x_cpu){
  splitname<-base::unlist(stringr::str_split(nam,"_"))
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]
  if (base::dir.exists(base::paste0(outpath,"/", fname,"/",chr ))){
    base::cat("\nDirectory already exists for",nam, "\n")
  }else{
    base::dir.create(outpath, showWarnings = FALSE)
    in_name<- base::paste0(inpath,nam,".fasta")

    dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")   #nchar(tdna) base::nchar(waxy1)= 295045
    dn1<-dna_1_to_wax(dna_1)   #base::length(dn) base::nchar(dn)
    full_length<-base::length(dn1)

    #this runs CG content
    winbox<-1*1e6
    #walklist<- CG_AT_content(dn1,nam,winbox=winbox)


    dn1<-NULL; dna_1<-NULL;base::gc()
    max_size<-2.0*1e6
    samplesize<-base::round(full_length/max_size)     #samplesize for long spectrum fft windows 150/2=75 100/2=50  50/2=25 20/2=10
    samplesize<-base::floor(samplesize/100)*100
    if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /50)*50   #full_length<-194*1e6 full_length<-210*1e6
    if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /20)*20
    if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /10)*10   #full_length<-6*1e6 full_length<-3*1e6
    if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /5)*5
    if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /2)*2     #full_length<-1.5*1e6 full_length<-3*1e6
    if(samplesize==0) samplesize<-1
    base::cat("\nSamplesize is", samplesize,"\n")

    # builds the files for spectrum - 0.5, 1.5, 4.5Mbp at a time
    twolists<-run_chromosome_parallel(nam=nam, fname=fname, inpath=inpath, outpath=outpath, atflag=atflag, AT_flag=AT_flag,
                                      majorlen=5000000,
                                      length_majorgroup=1000000, name_majorgroup="1Mbp",plotflag = plotflag,
                                      writeflag=writeflag, samplesize=samplesize, x_cpu=x_cpu)

    twolists<-run_chromosomelistNEW(nam=nam, fname=fname, inpath=inpath, outpath=outpath)

    All_list<-twolists$All_list

    fullwalk_list<-twolists$fullwalklist

    All_spec<-read_All_list(All_list)
    #utils::write.table(All_spec, paste0(outpath,"/", fname,"/",nam1, "/",fname,"_", nam,"_All_spec.txt"))

    whival<-(base::which(base::as.numeric(base::colnames(All_spec))>full_length+2500))
    if(base::length(whival)>0){
      All_spec[1:10,(whival)]
      All_spec<-All_spec[,-(whival)]
    }
    # #imagenan(All_spec)     tail(base::colnames(All_spec))
    base::gc()

    # create large images
    run_largeimage(All_spec,nam=nam, fname=fname, inpath=inpath, outpath=outpath, atflag=atflag, inmain="Individual 5Kbp window")
    base::gc()

    diff<-0
    run_largeimage(All_spec, nam=nam, fname=fname, inpath=inpath, outpath=outpath,atflag=atflag, rangeseq=NULL,pngflag=TRUE,samplesize=1,
                   r0=base::c(5,200),r1=NULL,r2=NULL,inmain="Individual 5Kbp window")
    base::gc()
  }
}



#' write_All_spec_DNAwalk
#'
#' Creates barplots of the spectra and DNAwalks. Needs to be run after the spectra have all been built.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
write_All_spec_DNAwalk <- function(nam=nam, fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath,
                                   pflag=FALSE, plotflag=FALSE, writeflag=FALSE){

  in_name<- base::paste0(inpath,nam,".fasta")    #
  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")   #nchar(tdna) base::nchar(waxy1)= 295045
  dn1<-dna_1_to_wax(dna_1)   #base::length(dn) base::nchar(dn)
  full_length<-base::length(dn1)

  splitname<-base::unlist(stringr::str_split(nam,"_"))   #fname<-fname
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]

  # write All_spec
  out_chromosome <- paste0(outpath, "/", fname,"/", chromosome)
  All_spec_file <- paste0(out_chromosome,"/", nam,"_All_spec.txt")

  if(base::file.exists(All_spec_file)){
    base::cat("\nfile exists: ", All_spec_file)
  } else {
    twolists<-run_chromosomelistNEW(nam, fname,inpath, outpath)
    All_list<-twolists$All_list
    fullwalk_list<-twolists$fullwalklist
    All_spec<-read_All_list(All_list)
    base::dir.create(out_chromosome, showWarnings = FALSE)
    utils::write.table(All_spec, All_spec_file)
    All_spec <- NULL
    gc()
  }

  # write total DNA walk for each part
  bigincre<-50000000;
  samplesizesmall=10
  fftlength<-fftlength1<-50000
  fftlengthlarge<-100000
  fftlengthlargest<-200000
  All_spec_long<-NULL
  DNAwalk_long<-NULL
  fractal_long<-NULL
  fractal_shorter<-NULL
  bprun<-base::seq(1,full_length,bigincre)

  for(sbp in bprun){
    base::gc()
    ebp<-sbp+bigincre-1  #change Nov 25
    if(ebp>full_length) ebp<-full_length

    walk_txtfile<-paste0(out_chromosome,"/dnawalk_",nam,"_", sbp,"_Table.txt")

    if(base::file.exists(walk_txtfile)){
      base::cat("\nfile exists: ", walk_txtfile)
    } else {
      base::cat("\nfile does not exist: construct and write ",walk_txtfile)
      twolists1<-run_chromosomelistNEW(nam=nam, fname=fname, inpath=inpath, outpath=outpath, majorlen=5000000,
                                       length_majorgroup=1000000, name_majorgroup="1Mbp",
                                       startgroup=sbp,endgroup=ebp,full_length=full_length)
      walk_listChr<-twolists1$fullwalklist
      All_walk<-read_walk_list(walk_listChr)
      utils::write.table(All_walk, walk_txtfile)
      gc()
    }
  }
}

#' merge_spectra
#'
#' Merges powers of rounded repeat lengths (currently sums).  Needs to be run after the spectra have all been built and merged for each part.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
merge_spectra <- function(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath){

  nam_list0 <- base::list.files(paste0(outpath,"/", fname,"/", chromosome,"/"))
  nam_list1 <- nam_list0[base::grep("All_spec.txt", nam_list0)]
  nam_list2 <- nam_list1[base::grep("merged_", nam_list1)]
  if (length(nam_list2)>0) {
    base::cat("\nfile exists")
  } else {
    for (i in 1:base::length(nam_list1)){
      file_nam <- nam_list1[i]
      All_spec_merge_file <- paste0(outpath, "/", fname,"/", chromosome,"/", "merged_", file_nam)
      base::cat(file_nam)
      All_spec <- base::as.matrix(utils::read.table(paste0(outpath, "/", fname,"/", chromosome,"/", file_nam), check.names = FALSE))

      All_spec_merged <-NULL
      All_spec_merged1 <-NULL
      All_spec_Freq <-NULL

      for (col in 1:base::ncol(All_spec)){
        base::cat(col, " \n")
        Freq_spec <- NULL
        Freq_spec <- base::cbind(Freq_spec, All_spec[,c(1,col)])
        Freq_spec[,1] <- base::row.names(Freq_spec)
        Freq_spec <- base::as.data.frame(Freq_spec)
        Freq_spec[,3] <-  base::sub("1/", "", Freq_spec[,1] , fixed = TRUE)
        Freq_spec[,4] <- base::as.factor(round(as.numeric(Freq_spec[,3])))
        colnames(Freq_spec) <- base::c("Freq", "Power", "1/Freq", "Freq_round")
        Freq_spec_by_Freq_round <- dplyr::group_by(Freq_spec, by=Freq_round)
        Freq_spec_by_Freq_round$Power <- as.numeric(Freq_spec_by_Freq_round$Power)
        Freq_spec_merged <- dplyr::summarise(Freq_spec_by_Freq_round, Power_merged = mean(Power, na.rm=TRUE))
        Freq_spec_merged <- Freq_spec_merged[base::order(Freq_spec_merged$by, decreasing = TRUE),]
        All_spec_merged <- base::cbind(All_spec_merged, Freq_spec_merged$Power_merged)
        All_spec_Freq <- base::as.character(Freq_spec_merged$by)
      }

      All_spec_merged1 <- base::cbind(All_spec_Freq, All_spec_merged)

      # save the merged version of All_spec
      utils::write.table(All_spec_merged1, paste0(outpath, "/", fname,"/", chromosome,"/", "merged_", file_nam))
    }
  }
}


#' join_chromosome_parts
#'
#' Creates barplots of the spectra and DNAwalks. Needs to be run after the spectra have all been built.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
join_chromosome_parts <- function(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath){

  # get list of part files that make up one chromosome
  nam_list0 <- base::list.files(inpath)
  nam_list1 <- tools::file_path_sans_ext(nam_list0)
  nam_list2 <- stringr::str_split(nam_list1, "_", simplify =TRUE)
  nam_list3 <- nam_list2[base::grep("part", nam_list2[,3]),]
  nam_list4 <- nam_list3[base::grep(paste0(chromosome,"p"), nam_list3[,3]),]

  if(nrow(nam_list3)==0) {
    base::cat("already a chromosome")
    nam_list5 <- nam_list1[base::grep(paste0(chromosome), nam_list1)]
  } else {
    if (is.null(nrow(nam_list4))){
      nam_list5 <- NULL
      nam_list5 <- base::paste0(nam_list4[1], "_", nam_list4[2], "_", nam_list4[3])
    } else {
      nam_list5 <- NULL
      for (i in 1:base::nrow(nam_list4)){
        nam_list5[i]<- base::paste0(nam_list4[i,1], "_", nam_list4[i,2], "_", nam_list4[i,3])
      }
    }
  }

  full_length_total <- NULL
  All_spec_Total <- NULL
  DNAwalk_Total <- NULL

  # get full length of chromosome
  fulllength_file <- paste0(outpath,"/", fname,"/", chromosome,"/",fname, "_", chromosome,"full_length.txt")
  if(base::file.exists(fulllength_file)){
    base::cat("\nfile exists: ", fulllength_file)
  } else {
    for (i in 1:base::length(nam_list5)){
      nam <- nam_list5[i]
      # get total chromosome length
      in_name<- base::paste0(inpath,nam,".fasta")
      dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")
      dn1<-dna_1_to_wax(dna_1)   #base::length(dn) base::nchar(dn)
      full_length_total[i] <- base::length(dn1)
    }
    utils::write.table(full_length_total, paste0(outpath,"/", fname,"/", chromosome,"/",fname, "_", chromosome,"full_length.txt"))
    gc()
  }
  # build chromosome scale spectra
  All_spec_total_file <- paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt")
  full_length <-  sum (full_length_total)
  if(base::file.exists(All_spec_total_file)){
    base::cat("\nfile exists: ", All_spec_total_file)
  } else {
    for (i in 1:base::length(nam_list5)){
      nam <- nam_list5[i]
      cat("Allspec ", nam)
      # join All_spec files for a given chromosome
      All_spec_tmp0 <- base::as.matrix(utils::read.table(paste0(outpath, "/", fname,"/", chromosome,"/merged_", nam,"_All_spec.txt"), header = TRUE, check.names = FALSE))
      All_spec_tmp1 <- All_spec_tmp0[,-1]
      All_spec_Total <-base::cbind(All_spec_Total, All_spec_tmp1)
      row.names(All_spec_Total) <- as.character(All_spec_tmp0[,1])
    }
    colnames(All_spec_Total) <-  as.character(seq.int(2501, ((ncol(All_spec_Total)-1)*5000)+2501, 5000))
    utils::write.table(All_spec_Total, All_spec_total_file)
    #All_spec_Total <- NULL
    gc()
  }

  # build total chromosome DNAwalk
  DNAwalk_total_file <- paste0(outpath, "/", fname,"/", chromosome,"/Total_dnawalk_every50_",fname, "_", chromosome,".txt")
  DNAwalk_endval <- data.frame(0,0)
  if(base::file.exists(DNAwalk_total_file)){
    base::cat("\nfile exists: ", DNAwalk_total_file)
  } else {
    for (i in 1:base::length(nam_list5)){
      #i=1
      nam <- nam_list5[i]
      cat("DNAwalk ",nam, " ")

      # join DNAwalks for all parts of chromosome
      DNAwalk1 <- paste0(outpath, "/", fname,"/", chromosome, "/dnawalk_",nam,"_1_Table.txt")
      if(base::file.exists(DNAwalk1)){
        DNAwalk_tmp1 <- utils::read.table(DNAwalk1, check.names = FALSE)
        DNAwalk_tmp1sub <- DNAwalk_tmp1[c(rep(FALSE,49),TRUE), ]
        DNAwalk_tmp2sub <- DNAwalk_tmp1sub
        DNAwalk_tmp2sub[,1] <- DNAwalk_tmp1sub[,1] + DNAwalk_endval[1,1]
        DNAwalk_tmp2sub[,2] <- DNAwalk_tmp1sub[,2] + DNAwalk_endval[1,2]

        DNAwalk_endval <- DNAwalk_tmp2sub[nrow(DNAwalk_tmp2sub),]
        DNAwalk_Total <- base::rbind(DNAwalk_Total, DNAwalk_tmp2sub)
      }

      DNAwalk2 <- paste0(outpath, "/", fname,"/", chromosome, "/dnawalk_",nam,"_50000001_Table.txt")
      if(base::file.exists(DNAwalk2)){
        DNAwalk_tmp3 <- utils::read.table(DNAwalk2, check.names = FALSE)
        DNAwalk_tmp3sub <- DNAwalk_tmp3[c(rep(FALSE,49),TRUE), ]
        DNAwalk_tmp4sub <- DNAwalk_tmp3sub
        DNAwalk_tmp4sub[,1] <- DNAwalk_tmp3sub[,1] + DNAwalk_endval[1,1]
        DNAwalk_tmp4sub[,2] <- DNAwalk_tmp3sub[,2] + DNAwalk_endval[1,2]

        DNAwalk_endval <- DNAwalk_tmp4sub[nrow(DNAwalk_tmp4sub),]
        DNAwalk_Total <- base::rbind(DNAwalk_Total, DNAwalk_tmp4sub)
        gc()
      }
    }
    utils::write.table(DNAwalk_Total, DNAwalk_total_file)
    gc()
  }
}

#' run_summary_hist
#'
#' Creates barplots of the spectra and DNAwalks. Needs to be run after the spectra have all been built.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_summary_hist <- function(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath, atflag=TRUE){

  # read in full length total
  full_length_total <- base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome,"full_length.txt"), header = TRUE, check.names = FALSE))

  # read in total All_spec_merged
  All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

  # remove blank/zero columns in All_spec
  col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

  # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
  All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
  All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

  # cut telomeres - need to cut exactly one bin though
  col_length1 <- as.numeric(colnames(All_spec1)[ncol(All_spec1)])
  if (ncol(All_spec1)>500) {
    All_spec2 <- All_spec1[,-c(1:400,(ncol(All_spec1)-400):ncol(All_spec1))]

    # cut last bin (Xbp from end ) from histograms ****

    col_length2 <- as.numeric(ncol(All_spec2))*5000
    numbins <-  floor(col_length2/2e+6)
    total_col <- (numbins*2e+6)/5000
    All_spec <- All_spec2#[,c(1:total_col)]

    full_length <-  as.numeric(ncol(All_spec))*5000
    if (full_length>5*1e6) {

      # setup bins
      #if(full_length>(50*1e6))numbins<-round(full_length/3e+6) else numbins<-round(full_length/3e+6)
      binsize<-base::round(full_length/numbins)
      fftlength<-5000
      numrange<-base::round(binsize/fftlength)
      if(full_length>(50*1e6)){numstd_bins<-0.5;numstdred<-1 } else { numstd_bins<-0.5;numstdred<-1 }

      run_barplots_chromosome(All_spec=All_spec,chromosome=chromosome,fname=fname,inpath=inpath, outpath=outpath, full_length=full_length, atflag=atflag, numstd=numstdred,  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(35,2000), binnum=numbins)
    }
  }
#   #-----------------------------
#   repeat_list<-base::rownames(All_spec)
#   if(!base::is.null(repeat_range)) {
#     Nrepeat_list<-base::as.numeric(base::gsub("1/","",repeat_list))
#     whichrow<-base::which(Nrepeat_list<=repeat_range[2] & Nrepeat_list>=repeat_range[1])
#     repeat_list<-repeat_list[whichrow]
#     Nrepeat_list<-Nrepeat_list[whichrow]
#   }
#
#   chromnam <- chromosome
#
#   min_powsum_seqval_list<-NULL
#   max_powsum_seqval_list<-NULL
#   goodrepeats<-NULL
#   min_mean_seqval_list<-NULL
#   max_mean_seqval_list<-NULL
#   min_N_seqval_list<-NULL
#   max_N_seqval_list<-NULL
#
#   for(repeat_val in repeat_list){
#     #repeat_val="5000"
#     base::cat("\n beginning ", repeat_val)
#     sortstr=FALSE
#     All_spec1<-All_spec
#     All_spec1[base::which(All_spec1==0.0)]<-NA
#     maxspec<-base::max(All_spec1[repeat_val, ],na.rm=TRUE)
#     meanspec<-base::mean(All_spec1[repeat_val, ],na.rm=TRUE)
#     stdspec<-stats::sd(All_spec1[repeat_val, ],na.rm=TRUE)
#
#     threshspec<-meanspec+numstd*stdspec    #base::length(peak_along_repeat)
#
#     peak_along_repeat<-pracma::findpeaks(x=All_spec[repeat_val, ], minpeakheight = threshspec,
#                                          minpeakdistance = 1, threshold = 0, npeaks = 0, sortstr = sortstr)  #change to false Feb 5 2023
#
#     if(base::length(peak_along_repeat)!=0){  #change June 17 2022
#       peak_along_repeat<-base::cbind(base::as.numeric(peak_along_repeat[,1])/maxspec,peak_along_repeat)
#       peak_along_repeat<-base::cbind(peak_along_repeat,base::as.numeric(base::colnames(All_spec)[peak_along_repeat[,3]]),
#                                      base::as.numeric(base::colnames(All_spec)[peak_along_repeat[,4]]),
#                                      base::as.numeric(base::colnames(All_spec)[peak_along_repeat[,5]]) )
#
#       base::colnames(peak_along_repeat)<-base::c("Percent of peak","Power","peak index","sindex","eindex","peak_bp","start_peak_bp", "end_peak_bp")
#     } #end of check to ensure peak value found
#
#     if(!base::is.null(peak_along_repeat)){
#       # perform sum of powers above mean+sigma2 level  (but only for repeats where mean+numstd *sigma holds)
#       # pow_list<- run_sum_bp(repeat_val,All_spec,chromnam,numstd2,numrange,pflag=FALSE)   #base::colnames(All_spec)
#
#       #All_spec1[base::which(All_spec1<=1e-5)]<-NA
#       threshspec<-meanspec+numstd2*stdspec
#       bpval<-base::as.numeric(base::colnames(All_spec))
#       deltabp<-(bpval[2]-bpval[1])
#       deltabp_onebin<-deltabp*numrange
#       start_bpseqval<-base::seq(bpval[1],bpval[base::length(bpval)],by= deltabp_onebin)
#       end_bpseqval<-start_bpseqval+deltabp_onebin
#       end_bpseqval[base::length(end_bpseqval)]<-bpval[base::length(bpval)]
#       which_ones<-base::which(All_spec1[repeat_val, ] >=threshspec)
#       power_sum<-NULL;power_mean<-NULL; N<-NULL
#       for(j in 1:base::length(start_bpseqval)){
#         sbpval<-start_bpseqval[j]
#         ebpval<-end_bpseqval[j]
#         whichval<-base::which(bpval>=sbpval & bpval<=ebpval)
#         if(!base::is.na(whichval[1])){
#           sval<-whichval[1]
#           eval<-whichval[base::length(whichval)]
#           which_notNA<-base::which(!base::is.na(All_spec1[repeat_val,sval: eval]))
#           bplength_onebin<-base::length(which_notNA)*deltabp
#           bp_notNA<-NULL
#           for (k in 2: base::length(which_notNA)){
#             bp_notNA<-base::c(bp_notNA,bpval[which_notNA[k]]-bpval[which_notNA[k-1]])
#           }
#
#           deltabpnotNA<-(bplength_onebin)/1.e6
#           which_val<-base::which(which_ones<eval & which_ones>=sval)
#           if(base::length(which_notNA)!=0){
#             meanabove<-base::mean(All_spec1[repeat_val,which_ones[which_val]] ,na.rm=TRUE)
#             power_mean<-base::c(power_mean,meanabove)
#             sumabove<-base::sum(All_spec1[repeat_val,which_ones[which_val]] ,na.rm=TRUE)/deltabpnotNA
#             power_sum<-base::c(power_sum,sumabove)
#             Nabove<-base::length(All_spec1[repeat_val,which_ones[which_val]] )/deltabpnotNA
#             N<-base::c(N,Nabove)
#             base::names(power_mean)[base::length(power_mean)]<-sbpval+bplength_onebin/2
#           } else{
#             meanabove<-NA;power_mean<-base::c(power_mean,NA);power_sum<-base::c(power_sum,NA) ;N<-base::c(N,NA)
#             base::names(power_mean)[base::length(power_mean)]<-sbpval+bplength_onebin/2
#           }
#         } else{
#           bplength_onebin<-0
#           meanabove<-NA;power_mean<-base::c(power_mean,NA);power_sum<-base::c(power_sum,NA) ;N<-base::c(N,NA)
#           base::names(power_mean)[base::length(power_mean)]<-sbpval+bplength_onebin/2
#         }
#       }
#       base::names(power_sum)<-names(N)<-names(power_mean)
#       if( !base::is.null(power_sum)){
# #
#         maxspec<-base::max(power_mean,na.rm=TRUE)
#         minspec<-base::min(power_mean,na.rm=TRUE)
#         mincol<-base::which(minspec==power_mean)[1]
#         maxcol<-base::which(maxspec==power_mean)[1]
#         min_mean_seqval<-base::as.numeric(names(power_mean[mincol]) )
#         max_mean_seqval<-base::as.numeric(names(power_mean[maxcol]) )
#
#         maxspec<-base::max(power_sum,na.rm=TRUE)
#         minspec<-base::min(power_sum,na.rm=TRUE)
#         mincol<-base::which(minspec==power_sum)[1]
#         maxcol<-base::which(maxspec==power_sum)[1]       #base::as.numeric
#         min_powsum_seqval<-base::as.numeric(names(power_sum[mincol]) )#+deltabp_onebin/2;
#         max_powsum_seqval<-base::as.numeric(names(power_sum[maxcol]) )#+deltabp_onebin/2
#
#         maxspec<-base::max(N,na.rm=TRUE)
#         minspec<-base::min(N,na.rm=TRUE)
#         mincol<-base::which(minspec==N)[1]
#         maxcol<-base::which(maxspec==N)[1]       #base::as.numeric
#         min_N_seqval<-base::as.numeric(names(N[mincol]) )#+deltabp_onebin/2;
#         max_N_seqval<-base::as.numeric(names(N[maxcol]) )#+deltabp_onebin/2
#       }
#
#       pow_list<-base::list(power_sum=power_sum,power_mean=power_mean,N=N,
#                            min_mean_seqval=min_mean_seqval,max_mean_seqval=max_mean_seqval,
#                            min_powsum_seqval=min_powsum_seqval,max_powsum_seqval=max_powsum_seqval,
#                            min_N_seqval=min_N_seqval,max_N_seqval=max_N_seqval )
#
#
#       sum_at_bp<-pow_list$power_sum
#       mean_at_bp<-pow_list$power_mean
#
#       goodrepeats<-base::c(goodrepeats,repeat_val)
#
#       min_powsum_seqval_list<-base::c(min_powsum_seqval_list,pow_list$min_powsum_seqval)
#       max_powsum_seqval_list<-base::c(max_powsum_seqval_list,pow_list$max_powsum_seqval)
#
#       min_mean_seqval_list<-base::c(min_mean_seqval_list,pow_list$min_mean_seqval)
#       max_mean_seqval_list<-base::c(max_mean_seqval_list,pow_list$max_mean_seqval)
#
#       min_N_seqval_list<-base::c(min_N_seqval_list,pow_list$min_N_seqval)
#       max_N_seqval_list<-base::c(max_N_seqval_list,pow_list$max_N_seqval)
#
#       N_at_bp<-pow_list$N
#     }
#   }
#
#   ofile<-base::paste0(outpath,"/",fname,"/", chromosome, "/",fname, "_", "histogram", numrange,"_",chromnam,"_s_",numstd2,"std_",numstd,".pdf")
#
#   base::cat("\n ouput to", ofile)
#   grDevices::pdf(file=ofile)
#
#   base::cat("\n Mean:     Min summary\n");base::print(base::summary(min_mean_seqval_list))
#   base::cat("\n Power Sum:Min summary\n");base::print(base::summary(min_powsum_seqval_list))
#   base::cat("\n N:        Min summary\n");base::print(base::summary(min_N_seqval_list))
#   base::cat("\n Mean:     Max summary\n");base::print(base::summary(max_mean_seqval_list))
#   base::cat("\n Power Sum:Max summary\n");base::print(base::summary(max_powsum_seqval_list))
#   base::cat("\n N:        Max summary\n");base::print(base::summary(max_N_seqval_list))
#
#   # https://www.datamentor.io/r-programming/histogram
#
#   # add single value at end of chromosome to fix plotting
#   min_powsum_seqval_list <- c(min_powsum_seqval_list, full_length)
#   min_mean_seqval_list <- c(min_mean_seqval_list, full_length)
#   min_N_seqval_list <- c(min_N_seqval_list, full_length)
#   max_mean_seqval_list <- c(max_mean_seqval_list, full_length)
#   max_powsum_seqval_list <- c(max_powsum_seqval_list, full_length)
#   max_N_seqval_list <- c(max_N_seqval_list, full_length)
#
#   #---------------
#   a<- graphics::hist(min_powsum_seqval_list,breaks=numbins,xaxp=base::c(0,full_length,numbins),xlab="",xlim=c(0,full_length),
#                      main= base::paste0(chromnam,"_min_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,"\nPower Sum Minimum values in the Sequence\nbp range ",
#                                         repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
#   base::cat("\nMin Power Sum:  midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")
#
#   a<-graphics::hist(min_mean_seqval_list, breaks=numbins, xaxp=base::c(0,full_length,numbins),xlab="",xlim=c(0,full_length),
#                     main= base::paste0(chromnam,"_min_mean_seqval_",numrange,"\ns_",numstd2,"std_",numstd,"\nMean: Minimum values in the Sequence\nbp range ",
#                                        repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
#   base::cat("\nMin Mean:       midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")
#
#
#   a<-graphics::hist(min_N_seqval_list,breaks=numbins,xaxp=base::c(0,full_length,numbins),xlab="",xlim=c(0,full_length),
#                     main= base::paste0(chromnam,"_min_N_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
#                                        "\nNumber: Minimum values in the Sequence\nbp range ",
#                                        repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
#   base::cat("\nMin N:          midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")
#
#   a<-graphics::hist(max_mean_seqval_list,breaks=numbins,xaxp=base::c(0,full_length,numbins),xlab="",xlim=c(0,full_length),
#                     main= base::paste0(chromnam,"_max_mean_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
#                                        "\nMean: Maximum values in the Sequence\nbp range ",
#                                        repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
#   base::cat("\nMax Mean:      midpointmax counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")
#
#   a<-graphics::hist(max_powsum_seqval_list,breaks=numbins,xaxp=base::c(0,full_length,numbins),xlab="",xlim=c(0,full_length),
#                     main= base::paste0(chromnam,"_max_pow_sum_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
#                                        "\nPower Sum Maximum values in the Sequence\nbp range ",
#                                        repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
#   base::cat("\nMax Power Sum: midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")
#
#   a<-graphics::hist(max_N_seqval_list,breaks=numbins,xaxp=base::c(0,full_length,numbins),xlab="",xlim=c(0,full_length),
#                     main= base::paste0(chromnam,"_max_N_seqval_",numrange,"\ns_",numstd2,"std_",numstd,
#                                        "\nNumber: Maximum values in the Sequence\nbp range ",
#                                        repeat_range[1],"_",repeat_range[2]),las=2,cex.axis=0.7)
#   base::cat("\nMax N:         midpoint max counts",a$mids[base::which(a$counts==base::max(a$counts))],"\n")
#
#   grDevices::dev.off()

}

#' run_summary_plots
#'
#' Creates spectra and DNAwalks plots for the whole chromosome. Needs to be run after the spectra and DNAwalk have all been built.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_summary_plots <- function(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath){

  # read in total All_spec
  All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"), check.names = FALSE))

  # read in full length total
  full_length_total <- base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome,"full_length.txt"), check.names = FALSE))
  full_length <-  sum(full_length_total)

  # remove blank/zero columns in All_spec
  colnames(All_spec0)[ncol(All_spec0)]
  sum (full_length_total)

  # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
  All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
  All_spec <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

  # https://r-graph-gallery.com/heatmap
  #stats::heatmap(All_spec[,c(1:100)], Rowv=FALSE, Colv=FALSE)

  pngflag=TRUE

  ofi <- paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_All_spec_")

  for (i in 1:length(row.names(All_spec))){
    row.names(All_spec)[i] <- paste0("1/", row.names(All_spec)[i])
  }

  r0=base::c(2,35)
  r1=base::c(35,2000)

  base::cat("\nofi:",ofi,"\n")
  #if(!base::is.null(r0))largeimagesub_NEW(All_spec,fname=fname, chromosome=chromosome,ofi,rangebp=r0,pngflag = pngflag)
  if(!base::is.null(r1))largeimagesub_NEW(All_spec,fname=fname, chromosome=chromosome,ofi,rangebp=r1,pngflag = pngflag, pdf_flag=TRUE)

  #------------------------------------------
  # plot total DNAwalk for whole chromosome

  #from walk_and_plot
  # read DNA walks for all parts of a chromosome
  # read in total_DNAwalk
  DNAwalk_long <-base::as.matrix(utils::read.table(paste0(outpath, "/", fname,"/", chromosome,"/Total_dnawalk_every50_",fname, "_", chromosome,".txt"), check.names = FALSE))

  # write to png/pdf
  ofile2 <-  paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk2D_total.png")
  base::cat("\n ouput to", ofile2)
  grDevices::png(file=ofile2)
  dnawalk <- as.data.frame(DNAwalk_long)
  walk_colours <- NULL
  for (col in grDevices::rainbow(n=round(nrow(dnawalk)/100000))){
    coltmp <- rep(col, 100000)
    walk_colours <- c(walk_colours, coltmp)
  }
  plot(dnawalk$AT, dnawalk$CG, col = walk_colours)
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk2D_total.pdf"))
  dnawalk <- as.data.frame(DNAwalk_long)
  walk_colours <- NULL
  for (col in grDevices::rainbow(n=round(nrow(dnawalk)/100000))){
    coltmp <- rep(col, 100000)
    walk_colours <- c(walk_colours, coltmp)
  }
  plot(dnawalk$AT, dnawalk$CG, col = walk_colours)
  grDevices::dev.off()

  #----------
  ofile3 <-  paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_AT_total.png")
  base::cat("\n ouput to", ofile3)
  grDevices::png(file=ofile3)
  # plot 1D DNAwalk
  dnawalk$rownames <- 50*c(1:nrow(dnawalk))
  plot(dnawalk$rownames, dnawalk$AT, col = walk_colours)
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_AT_total.pdf"))
  # plot 1D DNAwalk
  dnawalk$rownames <- 50*c(1:nrow(dnawalk))
  plot(dnawalk$rownames, dnawalk$AT, col = walk_colours)
  grDevices::dev.off()

  #-------------
  ofile4 <-  paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_CG_total.png")
  base::cat("\n ouput to", ofile4)
  grDevices::png(file=ofile4)
  # plot 1D DNAwalk
  plot(dnawalk$rownames, dnawalk$CG, col = walk_colours)
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_CG_total.pdf"))
  # plot 1D DNAwalk
  plot(dnawalk$rownames, dnawalk$CG, col = walk_colours)
  grDevices::dev.off()

}

#' run_diversity_plots
#'
#' Creates plots of diversity indices (Shannon, Simpson, Fisher's alpha and Pielou's J).
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_diversity_plots <- function(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath){
  # read in total All_spec_merged

  if(base::file.exists(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"))){

    # read in file
    base::cat("\nfile exists: ", paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"))
    All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

    # remove blank/zero columns in All_spec
    col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

    # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
    All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
    All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

    # cut telomeres - need to cut exactly one bin though
    #col_length1 <- as.numeric(colnames(All_spec1)[ncol(All_spec1)])
    All_spec2 <- All_spec1[,-c(1:100,(ncol(All_spec1)-100):ncol(All_spec1))]
    # cut last bin (Xbp from end ) from histograms ****

    col_length2 <- as.numeric(ncol(All_spec2))*5000
    numbins <-  floor(col_length2/2e+6)
    total_col <- (numbins*2e+6)/5000
    All_spec <- All_spec2#[,c(1:total_col)]

    full_length <-  as.numeric(ncol(All_spec))*5000

    full_length_value <-  as.numeric(ncol(All_spec0))*5000

    #-----------------------------
    # calculate species diversity
    #install.packages("vegan")
    #library(vegan)

    genome_pos <- colnames(All_spec)
    All_spec_repeatlengths <- rownames(All_spec)

    # try normalize the rows, for each repeat length, make a fraction of total
    # then run Shannon
    All_spec <- as.matrix(All_spec)
    All_spec_replen_sum <- apply(All_spec, 1, sum)
    All_spec_norm <- All_spec/All_spec_replen_sum
    All_spec_norm_t <- t(as.matrix(All_spec_norm))

    # diversity calculation
    # https://rdrr.io/cran/vegan/man/diversity.html

    Shannon_div <- vegan::diversity(All_spec_norm_t, index = "shannon")#this is the Shannon-Wiener index
    Pielou_div <- vegan::diversity(All_spec_norm_t, index = "shannon")/log(vegan::specnumber(All_spec_norm_t))
    Simpson_div <- vegan::diversity(All_spec_norm_t, index = "simpson")#this is the Simpson index

    Shannon_div <- as.numeric(Shannon_div)
    Pielou_div <- as.numeric(Pielou_div)
    Simpson_div <- as.numeric(Simpson_div)

    Shannon_div[which(Shannon_div<=2)] <- NA
    Pielou_div[which(Pielou_div<=0.3)] <- NA
    Simpson_div[which(Simpson_div<=0.3)] <- NA

    #---------------------
    # make windows and run Shannon on those windows

    rownames(All_spec_norm_t)
    colnames(All_spec_norm_t)
    nrow(All_spec_norm_t)

    wind_factor <- round(full_length/150000)
    Group <- vctrs::vec_rep_each(c(1:floor(nrow(All_spec_norm_t)/wind_factor)), wind_factor)
    Group_plus <- rep(max(Group)+1, nrow(All_spec_norm_t) - length(Group))
    Group_total <- c(Group, Group_plus)
    Shannon_div_wind <- vegan::diversity(All_spec_norm_t, groups=Group_total , index = "shannon")
    Shannon_div_wind <- as.numeric(Shannon_div_wind)
    Shannon_div_wind[which(Shannon_div<=2)] <- NA

    #---------------------
    # plots

    genome_pos_wind <- unique(Group_total)*wind_factor*5000
    Shannon_div_wind_pos <- cbind(genome_pos_wind, Shannon_div_wind)
    cent_wind <- Shannon_div_wind_pos[which(Shannon_div_wind_pos[,2] == min(Shannon_div_wind_pos[,2])),1]

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window_", wind_factor,".png"), width = 3000, height = 1000)
    plot(genome_pos_wind, Shannon_div_wind, type="l")
    grDevices::dev.off()

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm.png"), width = 1500, height = 500)
    plot(genome_pos, Shannon_div)
    grDevices::dev.off()

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Pielou_div_norm.png"), width = 1500, height = 500)
    plot(genome_pos, Pielou_div)
    grDevices::dev.off()

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Simpson_div_norm.png"), width = 1500, height = 500)
    plot(genome_pos, Simpson_div)
    grDevices::dev.off()

    #pdfs
    grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window_", wind_factor,".pdf"))
    plot(genome_pos_wind, Shannon_div_wind, type="l")
    grDevices::dev.off()

    grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm.pdf"))
    plot(genome_pos, Shannon_div)
    grDevices::dev.off()

    grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Pielou_div_norm.pdf"))
    plot(genome_pos, Pielou_div)
    grDevices::dev.off()

    grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Simpson_div_norm.pdf"))
    plot(genome_pos, Simpson_div)
    grDevices::dev.off()

    #-----------------------
    # moving average across Shannon that shows where min is
    # https://www.storybench.org/how-to-calculate-a-rolling-average-in-r/

    #---------------
    if (full_length > 0.5e6){
      bin_size=25
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent25 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_25.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_25.pdf"))
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()
    }else{cent25=0}
    #----------
    if (full_length > 0.5e6){
      bin_size=100
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent100 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_100.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_100.pdf"))
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

    }else{cent100=0}
    #----------
    if (full_length > 1.25e6){
      bin_size=250
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent250 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_250.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_250.pdf"))
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

    }else{cent250=0}
    #----------
    if (full_length > 2.5e6){
      bin_size=500
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent500 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_500.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_500.pdf"))
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

    }else{cent500=0}

    #----------
    if (full_length > 5e6){
      bin_size=1000
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent1000 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_1000.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_1000.pdf"))
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

    }else{cent1000=0}

    #-------------------
    # write out the centromere position
    cent_data0 <- rbind(cent25, cent100, cent250, cent500, cent1000)
    cent_data <- cent_data0*5000+7501
    cent_data[6] <- cent_wind
    cent_data <- as.data.frame(cent_data)

    rownames(cent_data) <- c(paste0(fname, "_", chromosome, "_cent25"),paste0(fname, "_", chromosome, "_cent100"), paste0(fname, "_", chromosome, "_cent250"), paste0(fname, "_", chromosome, "_cent500"), paste0(fname, "_", chromosome, "_cent1000"), paste0(fname, "_", chromosome,"_", wind_factor, "_centwind"))

    cent_data$full_length <- rep(as.character(full_length_value), nrow(cent_data), )
    cent_data$fname <- rep(fname, nrow(cent_data))
    cent_data$chromosome <- rep(chromosome, nrow(cent_data))

    utils::write.table(cent_data, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Centromere_MIN_Shannon.txt"),
                       append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)

    #--------------------
    # write out Shannon diversity data
    Shannon_div1 <- cbind(genome_pos, Shannon_div)

    utils::write.table(Shannon_div1, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Shannon_div.txt"),
                       append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

  } else {
    for (i in 1:sum(grepl(chromosome, nam_list1))){
      Spp_chr_part <- nam_list1[grepl(chromosome, nam_list1)][i]
      chr_part <- stringr::str_split(Spp_chr_part, "_", simplify =TRUE)[,3]

      out_chromosome <- paste0(outpath, "/", fname,"/", chr_part)
      All_spec_file <- paste0(out_chromosome,"/", Spp_chr_part,"_All_spec.txt")

      if(base::file.exists(All_spec_file)){
        base::cat("\nfile exists: ", All_spec_file)
        All_spec0<-base::as.matrix(utils::read.table(All_spec_file, header = TRUE, check.names = FALSE))
      } else {

        nam_list0a <- base::list.files(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/"))
        nam_list1a <- nam_list0a[base::grep("_spar1_Table.txt", nam_list0a)]
        if (length(base::grep("merged", nam_list1a))!=0){
          nam_list1a <- nam_list1a[-base::grep("merged", nam_list1a)]
        }
        for (i in 1:base::length(nam_list1a)){
          file_nam <- nam_list1a[i]
          All_spec_merge_file <- paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/merged_", file_nam)
          if(base::file.exists(All_spec_merge_file)){
            base::cat("\nfile exists: ", All_spec_merge_file)
          } else {
            base::cat(file_nam)
            All_spec <- base::as.matrix(utils::read.table(paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/", file_nam), check.names = FALSE))

            All_spec_merged <-NULL
            All_spec_merged1 <-NULL
            All_spec_Freq <-NULL

            for (col in 1:base::ncol(All_spec)){
              base::cat(col, " \n")
              Freq_spec <- NULL
              Freq_spec <- base::cbind(Freq_spec, All_spec[,c(1,col)])
              Freq_spec[,1] <- base::row.names(Freq_spec)
              Freq_spec <- base::as.data.frame(Freq_spec)
              Freq_spec[,3] <-  base::sub("1/", "", Freq_spec[,1] , fixed = TRUE)
              Freq_spec[,4] <- base::as.factor(round(as.numeric(Freq_spec[,3])))
              colnames(Freq_spec) <- base::c("Freq", "Power", "1/Freq", "Freq_round")
              Freq_spec_by_Freq_round <- dplyr::group_by(Freq_spec, by=Freq_round)
              Freq_spec_by_Freq_round$Power <- as.numeric(Freq_spec_by_Freq_round$Power)
              Freq_spec_merged <- dplyr::summarise(Freq_spec_by_Freq_round, Power_merged = mean(Power, na.rm=TRUE))
              Freq_spec_merged <- Freq_spec_merged[base::order(Freq_spec_merged$by, decreasing = TRUE),]
              All_spec_merged <- base::cbind(All_spec_merged, Freq_spec_merged$Power_merged)
              All_spec_Freq <- base::as.character(Freq_spec_merged$by)
            }

            All_spec_merged1 <- base::cbind(All_spec_Freq, All_spec_merged)

            # save the merged version of All_spec
            utils::write.table(All_spec_merged1, paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/merged_", file_nam))
          }
        }

        # build chromosome scale spectra
        All_spec_total_file <- paste0(outpath,"/", fname,"/",chr_part,"/Total_",fname, "_", chr_part,"_All_spec_merged.txt")
        All_spec_Total <- NULL
        if(base::file.exists(All_spec_total_file)){
          base::cat("\nfile exists: ", All_spec_total_file)
        } else {
          nam_list0a <- base::list.files(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/"))
          nam_list5 <- nam_list0a[base::grep("merged", nam_list0a)]
          for (i in 1:base::length(nam_list5)){
            nam <- nam_list5[i]
            cat("Allspec ", nam)
            # join All_spec files for a given chromosome
            All_spec_tmp0 <- base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/", nam), header = TRUE, check.names = FALSE))
            All_spec_tmp1 <- All_spec_tmp0[,-1]
            All_spec_Total <-base::cbind(All_spec_Total, All_spec_tmp1)
            row.names(All_spec_Total) <- as.character(All_spec_tmp0[,1])
          }
          colnames(All_spec_Total) <-  as.character(seq.int(2501, ((ncol(All_spec_Total)-1)*5000)+2501, 5000))
          utils::write.table(All_spec_Total, All_spec_total_file)

          gc()
        }
      }

      All_spec0 <-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chr_part,"/Total_",fname, "_", chr_part,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

      # remove blank/zero columns in All_spec
      col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

      # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
      All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
      All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

      # cut telomeres - need to cut exactly one bin though
      #col_length1 <- as.numeric(colnames(All_spec1)[ncol(All_spec1)])
      All_spec2 <- All_spec1[,-c(1:100,(ncol(All_spec1)-100):ncol(All_spec1))]
      # cut last bin (Xbp from end ) from histograms ****

      col_length2 <- as.numeric(ncol(All_spec2))*5000
      numbins <-  floor(col_length2/2e+6)
      total_col <- (numbins*2e+6)/5000
      All_spec <- All_spec2#[,c(1:total_col)]

      full_length <-  as.numeric(ncol(All_spec))*5000

      full_length_value <-  as.numeric(ncol(All_spec0))*5000

      #-----------------------------
      # calculate species diversity
      #install.packages("vegan")
      #library(vegan)
      genome_pos <- colnames(All_spec)
      All_spec_repeatlengths <- rownames(All_spec)

      # try normalize the rows, for each repeat length, make a fraction of total
      # then run Shannon

      All_spec <- as.matrix(All_spec)

      All_spec_replen_sum <- apply(All_spec, 1, sum)


      All_spec_norm <- All_spec/All_spec_replen_sum

      All_spec_norm_t <- t(as.matrix(All_spec_norm))

      # diversity calculation

      Shannon_div <- vegan::diversity(All_spec_norm_t, index = "shannon")#this is the Shannon-Wiener index
      Pielou_div <- vegan::diversity(All_spec_norm_t, index = "shannon")/log(vegan::specnumber(All_spec_norm_t))
      Simpson_div <- vegan::diversity(All_spec_norm_t, index = "simpson")#this is the Simpson index

      Shannon_div <- as.numeric(Shannon_div)
      Pielou_div <- as.numeric(Pielou_div)
      Simpson_div <- as.numeric(Simpson_div)

      Shannon_div[which(Shannon_div<=2)] <- NA
      Pielou_div[which(Pielou_div<=0.3)] <- NA
      Simpson_div[which(Simpson_div<=0.3)] <- NA

      #---------------------
      # make windows and run Shannon on those windows

      rownames(All_spec_norm_t)
      colnames(All_spec_norm_t)
      nrow(All_spec_norm_t)

      wind_factor <- round(full_length/150000)
      Group <- vctrs::vec_rep_each(c(1:floor(nrow(All_spec_norm_t)/wind_factor)), wind_factor)
      Group_plus <- rep(max(Group)+1, nrow(All_spec_norm_t) - length(Group))
      Group_total <- c(Group, Group_plus)
      Shannon_div_wind <- vegan::diversity(All_spec_norm_t, groups=Group_total , index = "shannon")
      Shannon_div_wind <- as.numeric(Shannon_div_wind)
      Shannon_div_wind[which(Shannon_div<=2)] <- NA

      #---------------------
      # plots

      genome_pos_wind <- unique(Group_total)*wind_factor*5000
      Shannon_div_wind_pos <- cbind(genome_pos_wind, Shannon_div_wind)
      cent_wind <- Shannon_div_wind_pos[which(Shannon_div_wind_pos[,2] == min(Shannon_div_wind_pos[,2])),1]

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window", wind_factor,".png"), width = 3000, height = 1000)
      plot(genome_pos_wind, Shannon_div_wind, type="l")
      grDevices::dev.off()

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm.png"), width = 1500, height = 500)
      plot(genome_pos, Shannon_div)
      grDevices::dev.off()

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Pielou_div_norm.png"), width = 1500, height = 500)
      plot(genome_pos, Pielou_div)
      grDevices::dev.off()

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Simpson_div_norm.png"),width = 1500, height = 500)
      plot(genome_pos, Simpson_div)
      grDevices::dev.off()

      #pdfs
      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window", wind_factor,".pdf"))
      plot(genome_pos_wind, Shannon_div_wind, type="l")
      grDevices::dev.off()

      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm.pdf"))
      plot(genome_pos, Shannon_div)
      grDevices::dev.off()

      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Pielou_div_norm.pdf"))
      plot(genome_pos, Pielou_div)
      grDevices::dev.off()

      grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Simpson_div_norm.pdf"))
      plot(genome_pos, Simpson_div)
      grDevices::dev.off()

      #-----------------------
      # moving average across Shannon that shows where min is
      # https://www.storybench.org/how-to-calculate-a-rolling-average-in-r/

      #---------------
      if (full_length > 0.5e6){
        bin_size=25
        roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

        cent25 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

        grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_25.png"), width = 1500, height = 500)
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

        grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_25.pdf"))
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

      }else{cent25=0}
      #----------
      if (full_length > 0.5e6){
        bin_size=100
        roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

        cent100 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

        grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_100.png"), width = 1500, height = 500)
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

        grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_100.pdf"))
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

      }else{cent100=0}

      #----------
      if (full_length > 1.25e6){
        bin_size=250
        roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

        cent250 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

        grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_250.png"), width = 1500, height = 500)
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

        grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_250.pdf"))
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

      }else{cent250=0}

      #----------
      if (full_length > 2.5e6){
        bin_size=500
        roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

        cent500 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

        grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_500.png"), width = 1500, height = 500)
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

        grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_500.pdf"))
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

      }else{cent500=0}

      #----------
      if (full_length > 5e6){
        bin_size=1000
        roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

        cent1000 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

        grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_1000.png"), width = 1500, height = 500)
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

        grDevices::pdf(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_1000.png"))
        base::plot(genome_pos, roll_mean_Shannon)
        grDevices::dev.off()

      }else{cent1000=0}

      #-------------------
      # write out the centromere position
      cent_data0 <- rbind(cent25, cent100, cent250, cent500, cent1000)
      cent_data <- cent_data0*5000+7501
      cent_data[6] <- cent_wind
      cent_data <- as.data.frame(cent_data)

      rownames(cent_data) <- c(paste0(fname, "_", chromosome, "_cent25"),paste0(fname, "_", chromosome, "_cent100"), paste0(fname, "_", chromosome, "_cent250"), paste0(fname, "_", chromosome, "_cent500"), paste0(fname, "_", chromosome, "_cent1000"), paste0(fname, "_", chromosome,"_", wind_factor, "_centwind"))

      cent_data$full_length <- rep(as.character(full_length_value), nrow(cent_data), )
      cent_data$fname <- rep(fname, nrow(cent_data))
      cent_data$chromosome <- rep(chromosome, nrow(cent_data))

      utils::write.table(cent_data, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Centromere_MIN_Shannon.txt"),
                         append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)
      #--------------------
      # write out Shannon diversity data
      Shannon_div1 <- cbind(genome_pos, Shannon_div)

      utils::write.table(Shannon_div1, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Shannon_div.txt"),
                         append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

    }
  }
}

#' run_diversity_plots_35_2000
#'
#' Creates plots of diversity indices (Shannon, Simpson, Fisher's alpha and Pielou's J).
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_diversity_plots_35_2000 <- function(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath){
  # read in total All_spec_merged

  if(base::file.exists(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"))){

    # read in file
    base::cat("\nfile exists: ", paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"))
    All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

    # remove blank/zero columns in All_spec
    col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

    # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
    All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
    All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

    # cut telomeres - need to cut exactly one bin though
    #col_length1 <- as.numeric(colnames(All_spec1)[ncol(All_spec1)])
    All_spec2 <- All_spec1[,-c(1:100,(ncol(All_spec1)-100):ncol(All_spec1))]
    # cut last bin (Xbp from end ) from histograms ****

    col_length2 <- as.numeric(ncol(All_spec2))*5000
    numbins <-  floor(col_length2/2e+6)
    total_col <- (numbins*2e+6)/5000
    All_spec <- All_spec2#[,c(1:total_col)]

    full_length <-  as.numeric(ncol(All_spec))*5000

    full_length_value <-  as.numeric(ncol(All_spec0))*5000

    #-----------------------------
    # calculate species diversity
    #install.packages("vegan")
    #library(vegan)

    genome_pos <- colnames(All_spec)
    All_spec_repeatlengths <- rownames(All_spec)

    # try normalize the rows, for each repeat length, make a fraction of total
    # then run Shannon
    All_spec <- as.matrix(All_spec)
    All_spec_replen_sum <- apply(All_spec, 1, sum)
    All_spec_norm <- All_spec/All_spec_replen_sum
    All_spec_norm_t0 <- t(as.matrix(All_spec_norm))

    # subset to only repeats >35bp
    All_spec_norm_t <- All_spec_norm_t0[,-which(as.numeric(colnames(All_spec_norm_t0))<35)]

    # diversity calculation
    # https://rdrr.io/cran/vegan/man/diversity.html

    Shannon_div <- vegan::diversity(All_spec_norm_t, index = "shannon")#this is the Shannon-Wiener index
    Pielou_div <- vegan::diversity(All_spec_norm_t, index = "shannon")/log(vegan::specnumber(All_spec_norm_t))
    Simpson_div <- vegan::diversity(All_spec_norm_t, index = "simpson")#this is the Simpson index

    Shannon_div <- as.numeric(Shannon_div)
    Pielou_div <- as.numeric(Pielou_div)
    Simpson_div <- as.numeric(Simpson_div)

    Shannon_div[which(Shannon_div<=2)] <- NA
    Pielou_div[which(Pielou_div<=0.3)] <- NA
    Simpson_div[which(Simpson_div<=0.3)] <- NA

    #---------------------
    # make windows and run Shannon on those windows

    rownames(All_spec_norm_t)
    colnames(All_spec_norm_t)
    nrow(All_spec_norm_t)

    wind_factor <- round(full_length/150e3)
    Group <- vctrs::vec_rep_each(c(1:floor(nrow(All_spec_norm_t)/wind_factor)), wind_factor)
    Group_plus <- rep(max(Group)+1, nrow(All_spec_norm_t) - length(Group))
    Group_total <- c(Group, Group_plus)
    Shannon_div_wind <- vegan::diversity(All_spec_norm_t, groups=Group_total , index = "shannon")
    Shannon_div_wind <- as.numeric(Shannon_div_wind)
    Shannon_div_wind[which(Shannon_div<=2)] <- NA

    #---------------------
    # plots

    genome_pos_wind <- unique(Group_total)*wind_factor*5000
    Shannon_div_wind_pos <- cbind(genome_pos_wind, Shannon_div_wind)
    cent_wind <- Shannon_div_wind_pos[which(Shannon_div_wind_pos[,2] == min(Shannon_div_wind_pos[,2])),1]

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window_35_", wind_factor,".png"), width = 3000, height = 1000)
    plot(genome_pos_wind, Shannon_div_wind, type="l")
    grDevices::dev.off()

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm_35.png"), width = 3000, height = 1000)
    plot(genome_pos, Shannon_div)
    grDevices::dev.off()

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Pielou_div_norm_35.png"), width = 3000, height = 1000)
    plot(genome_pos, Pielou_div)
    grDevices::dev.off()

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Simpson_div_norm_35.png"), width = 3000, height = 1000)
    plot(genome_pos, Simpson_div)
    grDevices::dev.off()

    #-----------------------
    # moving average across Shannon that shows where min is
    # https://www.storybench.org/how-to-calculate-a-rolling-average-in-r/

    #---------------
    bin_size=25
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent25 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_25_35.png"), width = 1500, height = 500)
    base::plot(genome_pos, roll_mean_Shannon)
    grDevices::dev.off()

    #----------
    bin_size=100
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent100 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_100_35.png"), width = 1500, height = 500)
    base::plot(genome_pos, roll_mean_Shannon)
    grDevices::dev.off()

    #----------
    bin_size=250
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent250 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_250_35.png"), width = 1500, height = 500)
    base::plot(genome_pos, roll_mean_Shannon)
    grDevices::dev.off()

    #----------
    bin_size=500
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent500 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_500_35.png"), width = 1500, height = 500)
    base::plot(genome_pos, roll_mean_Shannon)
    grDevices::dev.off()


    #----------
    bin_size=1000
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent1000 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_1000_35.png"), width = 1500, height = 500)
    base::plot(genome_pos, roll_mean_Shannon)
    grDevices::dev.off()

    #-------------------
    # write out the centromere position
    cent_data0 <- rbind(cent25, cent100, cent250, cent500, cent1000)
    cent_data <- cent_data0*5000+7501
    cent_data[6] <- cent_wind
    cent_data <- as.data.frame(cent_data)

    rownames(cent_data) <- c(paste0(fname, "_", chromosome, "_cent25"),paste0(fname, "_", chromosome, "_cent100"), paste0(fname, "_", chromosome, "_cent250"), paste0(fname, "_", chromosome, "_cent500"), paste0(fname, "_", chromosome, "_cent1000"), paste0(fname, "_", chromosome,"_", wind_factor, "_centwind"))

    cent_data$full_length <- rep(as.character(full_length_value), nrow(cent_data), )
    cent_data$fname <- rep(fname, nrow(cent_data))
    cent_data$chromosome <- rep(chromosome, nrow(cent_data))

    utils::write.table(cent_data, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Centromere_MIN_Shannon_35.txt"),
                       append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)

    #--------------------
    # write out Shannon diversity data
    Shannon_div1 <- cbind(genome_pos, Shannon_div)

    utils::write.table(Shannon_div1, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Shannon_div_35.txt"),
                       append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

  } else {
    for (i in 1:sum(grepl(chromosome, nam_list1))){
      Spp_chr_part <- nam_list1[grepl(chromosome, nam_list1)][i]
      chr_part <- stringr::str_split(Spp_chr_part, "_", simplify =TRUE)[,3]

      out_chromosome <- paste0(outpath, "/", fname,"/", chr_part)
      All_spec_file <- paste0(out_chromosome,"/", Spp_chr_part,"_All_spec.txt")

      if(base::file.exists(All_spec_file)){
        base::cat("\nfile exists: ", All_spec_file)
        All_spec0<-base::as.matrix(utils::read.table(All_spec_file, header = TRUE, check.names = FALSE))
      } else {

        nam_list0a <- base::list.files(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/"))
        nam_list1a <- nam_list0a[base::grep("_spar1_Table.txt", nam_list0a)]
        if (length(base::grep("merged", nam_list1a))!=0){
          nam_list1a <- nam_list1a[-base::grep("merged", nam_list1a)]
        }
        for (i in 1:base::length(nam_list1a)){
          file_nam <- nam_list1a[i]
          All_spec_merge_file <- paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/merged_", file_nam)
          if(base::file.exists(All_spec_merge_file)){
            base::cat("\nfile exists: ", All_spec_merge_file)
          } else {
            base::cat(file_nam)
            All_spec <- base::as.matrix(utils::read.table(paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/", file_nam), check.names = FALSE))

            All_spec_merged <-NULL
            All_spec_merged1 <-NULL
            All_spec_Freq <-NULL

            for (col in 1:base::ncol(All_spec)){
              base::cat(col, " \n")
              Freq_spec <- NULL
              Freq_spec <- base::cbind(Freq_spec, All_spec[,c(1,col)])
              Freq_spec[,1] <- base::row.names(Freq_spec)
              Freq_spec <- base::as.data.frame(Freq_spec)
              Freq_spec[,3] <-  base::sub("1/", "", Freq_spec[,1] , fixed = TRUE)
              Freq_spec[,4] <- base::as.factor(round(as.numeric(Freq_spec[,3])))
              colnames(Freq_spec) <- base::c("Freq", "Power", "1/Freq", "Freq_round")
              Freq_spec_by_Freq_round <- dplyr::group_by(Freq_spec, by=Freq_round)
              Freq_spec_by_Freq_round$Power <- as.numeric(Freq_spec_by_Freq_round$Power)
              Freq_spec_merged <- dplyr::summarise(Freq_spec_by_Freq_round, Power_merged = mean(Power, na.rm=TRUE))
              Freq_spec_merged <- Freq_spec_merged[base::order(Freq_spec_merged$by, decreasing = TRUE),]
              All_spec_merged <- base::cbind(All_spec_merged, Freq_spec_merged$Power_merged)
              All_spec_Freq <- base::as.character(Freq_spec_merged$by)
            }

            All_spec_merged1 <- base::cbind(All_spec_Freq, All_spec_merged)

            # save the merged version of All_spec
            utils::write.table(All_spec_merged1, paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/merged_", file_nam))
          }
        }

        # build chromosome scale spectra
        All_spec_total_file <- paste0(outpath,"/", fname,"/",chr_part,"/Total_",fname, "_", chr_part,"_All_spec_merged.txt")
        All_spec_Total <- NULL
        if(base::file.exists(All_spec_total_file)){
          base::cat("\nfile exists: ", All_spec_total_file)
        } else {
          nam_list0a <- base::list.files(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/"))
          nam_list5 <- nam_list0a[base::grep("merged", nam_list0a)]
          for (i in 1:base::length(nam_list5)){
            nam <- nam_list5[i]
            cat("Allspec ", nam)
            # join All_spec files for a given chromosome
            All_spec_tmp0 <- base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/", nam), header = TRUE, check.names = FALSE))
            All_spec_tmp1 <- All_spec_tmp0[,-1]
            All_spec_Total <-base::cbind(All_spec_Total, All_spec_tmp1)
            row.names(All_spec_Total) <- as.character(All_spec_tmp0[,1])
          }
          colnames(All_spec_Total) <-  as.character(seq.int(2501, ((ncol(All_spec_Total)-1)*5000)+2501, 5000))
          utils::write.table(All_spec_Total, All_spec_total_file)

          gc()
        }
      }

      All_spec0 <-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chr_part,"/Total_",fname, "_", chr_part,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

      # remove blank/zero columns in All_spec
      col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

      # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
      All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
      All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

      # cut telomeres - need to cut exactly one bin though
      #col_length1 <- as.numeric(colnames(All_spec1)[ncol(All_spec1)])
      All_spec2 <- All_spec1[,-c(1:100,(ncol(All_spec1)-100):ncol(All_spec1))]
      # cut last bin (Xbp from end ) from histograms ****

      col_length2 <- as.numeric(ncol(All_spec2))*5000
      numbins <-  floor(col_length2/2e+6)
      total_col <- (numbins*2e+6)/5000
      All_spec <- All_spec2#[,c(1:total_col)]

      full_length <-  as.numeric(ncol(All_spec))*5000

      full_length_value <-  as.numeric(ncol(All_spec0))*5000

      #-----------------------------
      # calculate species diversity
      #install.packages("vegan")
      #library(vegan)
      genome_pos <- colnames(All_spec)
      All_spec_repeatlengths <- rownames(All_spec)

      # try normalize the rows, for each repeat length, make a fraction of total
      # then run Shannon

      All_spec <- as.matrix(All_spec)
      All_spec_replen_sum <- apply(All_spec, 1, sum)
      All_spec_norm <- All_spec/All_spec_replen_sum
      All_spec_norm_t0 <- t(as.matrix(All_spec_norm))

      # subset to only repeats >35bp
      All_spec_norm_t <- All_spec_norm_t0[,-which(as.numeric(colnames(All_spec_norm_t0))<35)]

      # diversity calculation

      Shannon_div <- vegan::diversity(All_spec_norm_t, index = "shannon")#this is the Shannon-Wiener index
      Pielou_div <- vegan::diversity(All_spec_norm_t, index = "shannon")/log(vegan::specnumber(All_spec_norm_t))
      Simpson_div <- vegan::diversity(All_spec_norm_t, index = "simpson")#this is the Simpson index

      Shannon_div <- as.numeric(Shannon_div)
      Pielou_div <- as.numeric(Pielou_div)
      Simpson_div <- as.numeric(Simpson_div)

      Shannon_div[which(Shannon_div<=2)] <- NA
      Pielou_div[which(Pielou_div<=0.3)] <- NA
      Simpson_div[which(Simpson_div<=0.3)] <- NA

      #---------------------
      # make windows and run Shannon on those windows

      rownames(All_spec_norm_t)
      colnames(All_spec_norm_t)
      nrow(All_spec_norm_t)

      wind_factor <- round(full_length/150000)
      Group <- vctrs::vec_rep_each(c(1:floor(nrow(All_spec_norm_t)/wind_factor)), wind_factor)
      Group_plus <- rep(max(Group)+1, nrow(All_spec_norm_t) - length(Group))
      Group_total <- c(Group, Group_plus)
      Shannon_div_wind <- vegan::diversity(All_spec_norm_t, groups=Group_total , index = "shannon")
      Shannon_div_wind <- as.numeric(Shannon_div_wind)
      Shannon_div_wind[which(Shannon_div<=2)] <- NA

      #---------------------
      # plots

      genome_pos_wind <- unique(Group_total)*wind_factor*5000
      Shannon_div_wind_pos <- cbind(genome_pos_wind, Shannon_div_wind)
      cent_wind <- Shannon_div_wind_pos[which(Shannon_div_wind_pos[,2] == min(Shannon_div_wind_pos[,2])),1]

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window_35_", wind_factor,".png"), width = 3000, height = 1000)
      plot(genome_pos_wind, Shannon_div_wind, type="l")
      grDevices::dev.off()


      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm_35.png"), width = 1500, height = 500)
      plot(genome_pos, Shannon_div)
      grDevices::dev.off()

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Pielou_div_norm_35.png"), width = 1500, height = 500)
      plot(genome_pos, Pielou_div)
      grDevices::dev.off()

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Simpson_div_norm_35.png"),width = 1500, height = 500)
      plot(genome_pos, Simpson_div)
      grDevices::dev.off()


      #-----------------------
      # moving average across Shannon that shows where min is
      # https://www.storybench.org/how-to-calculate-a-rolling-average-in-r/

      #---------------
      bin_size=25
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent25 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_25_35.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      #----------
      bin_size=100
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent100 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_100_35.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      #----------
      bin_size=250
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent250 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_250_35.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      #----------
      bin_size=500
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent500 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_500_35.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()


      #----------
      bin_size=1000
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent1000 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_roll_mean_Shannon_1000_35.png"), width = 1500, height = 500)
      base::plot(genome_pos, roll_mean_Shannon)
      grDevices::dev.off()

      #-------------------
      # write out the centromere position
      cent_data0 <- rbind(cent25, cent100, cent250, cent500, cent1000)
      cent_data <- cent_data0*5000+7501
      cent_data[6] <- cent_wind
      cent_data <- as.data.frame(cent_data)

      rownames(cent_data) <- c(paste0(fname, "_", chromosome, "_cent25"),paste0(fname, "_", chromosome, "_cent100"), paste0(fname, "_", chromosome, "_cent250"), paste0(fname, "_", chromosome, "_cent500"), paste0(fname, "_", chromosome, "_cent1000"), paste0(fname, "_", chromosome,"_", wind_factor, "_centwind"))

      cent_data$full_length <- rep(as.character(full_length_value), nrow(cent_data), )
      cent_data$fname <- rep(fname, nrow(cent_data))
      cent_data$chromosome <- rep(chromosome, nrow(cent_data))

      utils::write.table(cent_data, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Centromere_MIN_Shannon_35.txt"),
                         append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)
      #--------------------
      # write out Shannon diversity data
      Shannon_div1 <- cbind(genome_pos, Shannon_div)

      utils::write.table(Shannon_div1, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Shannon_div_35.txt"),
                         append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

    }
  }
}

#' run_diversity_plots_no_telomere
#'
#' Creates plots of diversity indices (Shannon, Simpson, Fisher's alpha and Pielou's J).
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_diversity_plots_no_telomere <- function(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath){
  # read in total All_spec_merged

  if(base::file.exists(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"))){

    # read in file
    base::cat("\nfile exists: ", paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"))
    All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

     # remove blank/zero columns in All_spec
    col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

    # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
    All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
    All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

    # cut telomeres - need to cut exactly one bin though
    telo_cut <- round(ncol(All_spec0)*0.04)
    All_spec2 <- All_spec1[,-c(1:telo_cut,(ncol(All_spec1)-telo_cut):ncol(All_spec1))]
    # cut last bin (Xbp from end ) from histograms ****

    col_length2 <- as.numeric(ncol(All_spec2))*5000
    numbins <-  floor(col_length2/2e+6)
    total_col <- (numbins*2e+6)/5000
    All_spec <- All_spec2#[,c(1:total_col)]

    full_length <-  as.numeric(ncol(All_spec))*5000
    full_length_value <-  as.numeric(ncol(All_spec0))*5000
    #-----------------------------
    # calculate species diversity
    #install.packages("vegan")
    #library(vegan)

    genome_pos <- colnames(All_spec)
    All_spec_repeatlengths <- rownames(All_spec)

    # try normalize the rows, for each repeat length, make a fraction of total
    # then run Shannon
    All_spec <- as.matrix(All_spec)
    All_spec_replen_sum <- apply(All_spec, 1, sum)
    All_spec_norm <- All_spec/All_spec_replen_sum
    All_spec_norm_t <- t(as.matrix(All_spec_norm))

    # diversity calculation
    # https://rdrr.io/cran/vegan/man/diversity.html

    Shannon_div <- vegan::diversity(All_spec_norm_t, index = "shannon")#this is the Shannon-Wiener index
    Pielou_div <- vegan::diversity(All_spec_norm_t, index = "shannon")/log(vegan::specnumber(All_spec_norm_t))
    Simpson_div <- vegan::diversity(All_spec_norm_t, index = "simpson")#this is the Simpson index

    Shannon_div <- as.numeric(Shannon_div)
    Pielou_div <- as.numeric(Pielou_div)
    Simpson_div <- as.numeric(Simpson_div)

    Shannon_div[which(Shannon_div<=2)] <- NA
    Pielou_div[which(Pielou_div<=0.3)] <- NA
    Simpson_div[which(Simpson_div<=0.3)] <- NA

    #---------------------
    # make windows and run Shannon on those windows

    rownames(All_spec_norm_t)
    colnames(All_spec_norm_t)
    nrow(All_spec_norm_t)

    wind_factor <- round(full_length/150000)
    Group <- vctrs::vec_rep_each(c(1:floor(nrow(All_spec_norm_t)/wind_factor)), wind_factor)
    Group_plus <- rep(max(Group)+1, nrow(All_spec_norm_t) - length(Group))
    Group_total <- c(Group, Group_plus)
    Shannon_div_wind <- vegan::diversity(All_spec_norm_t, groups=Group_total , index = "shannon")
    Shannon_div_wind <- as.numeric(Shannon_div_wind)
    Shannon_div_wind[which(Shannon_div<=2)] <- NA

    #---------------------
    # plots

    genome_pos_wind <- unique(Group_total)*wind_factor*5000
    Shannon_div_wind_pos <- cbind(genome_pos_wind, Shannon_div_wind)
    cent_wind <- Shannon_div_wind_pos[which(Shannon_div_wind_pos[,2] == min(Shannon_div_wind_pos[,2])),1]

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window_no_telo_", wind_factor,".png"), width = 3000, height = 1000)
    plot(genome_pos_wind, Shannon_div_wind, type="l")
    grDevices::dev.off()

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm_no_telo.png"), width = 3000, height = 1000)
    plot(genome_pos, Shannon_div)
    grDevices::dev.off()


    #-----------------------
    # moving average across Shannon that shows where min is
    # https://www.storybench.org/how-to-calculate-a-rolling-average-in-r/

    #---------------
    bin_size=25
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent25 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #----------
    bin_size=100
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent100 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #----------
    bin_size=250
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent250 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #----------
    bin_size=500
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent500 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #----------
    bin_size=1000
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent1000 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #-------------------
    # write out the centromere position
    cent_data0 <- rbind(cent25, cent100, cent250, cent500, cent1000)
    cent_data <- cent_data0*5000+7501
    cent_data[6] <- cent_wind
    cent_data <- as.data.frame(cent_data)

    rownames(cent_data) <- c(paste0(fname, "_", chromosome, "_cent25"),paste0(fname, "_", chromosome, "_cent100"), paste0(fname, "_", chromosome, "_cent250"), paste0(fname, "_", chromosome, "_cent500"), paste0(fname, "_", chromosome, "_cent1000"), paste0(fname, "_", chromosome,"_", wind_factor, "_centwind"))

    cent_data$full_length <- rep(as.character(full_length_value), nrow(cent_data), )
    cent_data$fname <- rep(fname, nrow(cent_data))
    cent_data$chromosome <- rep(chromosome, nrow(cent_data))

    utils::write.table(cent_data, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Centromere_MIN_Shannon_no_telo.txt"),
                       append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)
    #--------------------
    # write out Shannon diversity data
    Shannon_div1 <- cbind(genome_pos, Shannon_div)

    utils::write.table(Shannon_div1, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Shannon_div_no_telo.txt"),
                       append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

    } else {
    for (i in 1:sum(grepl(chromosome, nam_list1))){
      Spp_chr_part <- nam_list1[grepl(chromosome, nam_list1)][i]
      chr_part <- stringr::str_split(Spp_chr_part, "_", simplify =TRUE)[,3]

      out_chromosome <- paste0(outpath, "/", fname,"/", chr_part)
      All_spec_file <- paste0(out_chromosome,"/", Spp_chr_part,"_All_spec.txt")

      if(base::file.exists(All_spec_file)){
        base::cat("\nfile exists: ", All_spec_file)
        All_spec0<-base::as.matrix(utils::read.table(All_spec_file, header = TRUE, check.names = FALSE))
      } else {

        nam_list0a <- base::list.files(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/"))
        nam_list1a <- nam_list0a[base::grep("_spar1_Table.txt", nam_list0a)]
        if (length(base::grep("merged", nam_list1a))!=0){
        nam_list1a <- nam_list1a[-base::grep("merged", nam_list1a)]
        }
        for (i in 1:base::length(nam_list1a)){
          file_nam <- nam_list1a[i]
          All_spec_merge_file <- paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/merged_", file_nam)
          if(base::file.exists(All_spec_merge_file)){
            base::cat("\nfile exists: ", All_spec_merge_file)
          } else {
            base::cat(file_nam)
            All_spec <- base::as.matrix(utils::read.table(paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/", file_nam), check.names = FALSE))

            All_spec_merged <-NULL
            All_spec_merged1 <-NULL
            All_spec_Freq <-NULL

            for (col in 1:base::ncol(All_spec)){
              base::cat(col, " \n")
              Freq_spec <- NULL
              Freq_spec <- base::cbind(Freq_spec, All_spec[,c(1,col)])
              Freq_spec[,1] <- base::row.names(Freq_spec)
              Freq_spec <- base::as.data.frame(Freq_spec)
              Freq_spec[,3] <-  base::sub("1/", "", Freq_spec[,1] , fixed = TRUE)
              Freq_spec[,4] <- base::as.factor(round(as.numeric(Freq_spec[,3])))
              colnames(Freq_spec) <- base::c("Freq", "Power", "1/Freq", "Freq_round")
              Freq_spec_by_Freq_round <- dplyr::group_by(Freq_spec, by=Freq_round)
              Freq_spec_by_Freq_round$Power <- as.numeric(Freq_spec_by_Freq_round$Power)
              Freq_spec_merged <- dplyr::summarise(Freq_spec_by_Freq_round, Power_merged = mean(Power, na.rm=TRUE))
              Freq_spec_merged <- Freq_spec_merged[base::order(Freq_spec_merged$by, decreasing = TRUE),]
              All_spec_merged <- base::cbind(All_spec_merged, Freq_spec_merged$Power_merged)
              All_spec_Freq <- base::as.character(Freq_spec_merged$by)
            }

            All_spec_merged1 <- base::cbind(All_spec_Freq, All_spec_merged)

            # save the merged version of All_spec
            utils::write.table(All_spec_merged1, paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/merged_", file_nam))
          }
        }

        # build chromosome scale spectra
        All_spec_total_file <- paste0(outpath,"/", fname,"/",chr_part,"/Total_",fname, "_", chr_part,"_All_spec_merged.txt")
        All_spec_Total <- NULL
        if(base::file.exists(All_spec_total_file)){
          base::cat("\nfile exists: ", All_spec_total_file)
        } else {
          nam_list0a <- base::list.files(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/"))
          nam_list5 <- nam_list0a[base::grep("merged", nam_list0a)]
          for (i in 1:base::length(nam_list5)){
            nam <- nam_list5[i]
            cat("Allspec ", nam)
            # join All_spec files for a given chromosome
            All_spec_tmp0 <- base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/", nam), header = TRUE, check.names = FALSE))
            All_spec_tmp1 <- All_spec_tmp0[,-1]
            All_spec_Total <-base::cbind(All_spec_Total, All_spec_tmp1)
            row.names(All_spec_Total) <- as.character(All_spec_tmp0[,1])
          }
          colnames(All_spec_Total) <-  as.character(seq.int(2501, ((ncol(All_spec_Total)-1)*5000)+2501, 5000))
          utils::write.table(All_spec_Total, All_spec_total_file)

          gc()
        }
      }

      All_spec0 <-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chr_part,"/Total_",fname, "_", chr_part,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

      # remove blank/zero columns in All_spec
      col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

      # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
      All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
      All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

      # cut telomeres - need to cut exactly one bin though
      telo_cut <- round(ncol(All_spec0)*0.04)
      All_spec2 <- All_spec1[,-c(1:telo_cut,(ncol(All_spec1)-telo_cut):ncol(All_spec1))]
      # cut last bin (Xbp from end ) from histograms ****

      col_length2 <- as.numeric(ncol(All_spec2))*5000
      numbins <-  floor(col_length2/2e+6)
      total_col <- (numbins*2e+6)/5000
      All_spec <- All_spec2#[,c(1:total_col)]

      full_length <-  as.numeric(ncol(All_spec))*5000

      full_length_value <-  as.numeric(ncol(All_spec0))*5000

      #-----------------------------
      # calculate species diversity
      #install.packages("vegan")
      #library(vegan)
      genome_pos <- colnames(All_spec)
      All_spec_repeatlengths <- rownames(All_spec)

      # try normalize the rows, for each repeat length, make a fraction of total
      # then run Shannon

      All_spec <- as.matrix(All_spec)

      All_spec_replen_sum <- apply(All_spec, 1, sum)


      All_spec_norm <- All_spec/All_spec_replen_sum

      All_spec_norm_t <- t(as.matrix(All_spec_norm))

      # diversity calculation

      Shannon_div <- vegan::diversity(All_spec_norm_t, index = "shannon")#this is the Shannon-Wiener index
      Pielou_div <- vegan::diversity(All_spec_norm_t, index = "shannon")/log(vegan::specnumber(All_spec_norm_t))
      Simpson_div <- vegan::diversity(All_spec_norm_t, index = "simpson")#this is the Simpson index

      Shannon_div <- as.numeric(Shannon_div)
      Pielou_div <- as.numeric(Pielou_div)
      Simpson_div <- as.numeric(Simpson_div)

      Shannon_div[which(Shannon_div<=2)] <- NA
      Pielou_div[which(Pielou_div<=0.3)] <- NA
      Simpson_div[which(Simpson_div<=0.3)] <- NA

      #---------------------
      # make windows and run Shannon on those windows

      rownames(All_spec_norm_t)
      colnames(All_spec_norm_t)
      nrow(All_spec_norm_t)

      wind_factor <- round(full_length/150000)
      Group <- vctrs::vec_rep_each(c(1:floor(nrow(All_spec_norm_t)/wind_factor)), wind_factor)
      Group_plus <- rep(max(Group)+1, nrow(All_spec_norm_t) - length(Group))
      Group_total <- c(Group, Group_plus)
      Shannon_div_wind <- vegan::diversity(All_spec_norm_t, groups=Group_total , index = "shannon")
      Shannon_div_wind <- as.numeric(Shannon_div_wind)
      Shannon_div_wind[which(Shannon_div<=2)] <- NA

      #---------------------
      # plots

      genome_pos_wind <- unique(Group_total)*wind_factor*5000
      Shannon_div_wind_pos <- cbind(genome_pos_wind, Shannon_div_wind)
      cent_wind <- Shannon_div_wind_pos[which(Shannon_div_wind_pos[,2] == min(Shannon_div_wind_pos[,2])),1]

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window_no_telo_", wind_factor,".png"), width = 3000, height = 1000)
      plot(genome_pos_wind, Shannon_div_wind, type="l")
      grDevices::dev.off()


      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm_no_telo.png"), width = 1500, height = 500)
      plot(genome_pos, Shannon_div)
      grDevices::dev.off()

      #-----------------------
      # moving average across Shannon that shows where min is
      # https://www.storybench.org/how-to-calculate-a-rolling-average-in-r/

      #---------------
      bin_size=25
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent25 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))


      #----------
      bin_size=100
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent100 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))


      #----------
      bin_size=250
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent250 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))


      #----------
      bin_size=500
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent500 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      #----------
      bin_size=1000
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent1000 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))


      #-------------------
      # write out the centromere position
      cent_data0 <- rbind(cent25, cent100, cent250, cent500, cent1000)
      cent_data <- cent_data0*5000+7501
      cent_data[6] <- cent_wind
      cent_data <- as.data.frame(cent_data)

      rownames(cent_data) <- c(paste0(fname, "_", chromosome, "_cent25"),paste0(fname, "_", chromosome, "_cent100"), paste0(fname, "_", chromosome, "_cent250"), paste0(fname, "_", chromosome, "_cent500"), paste0(fname, "_", chromosome, "_cent1000"), paste0(fname, "_", chromosome,"_", wind_factor, "_centwind"))

      cent_data$full_length <- rep(as.character(full_length_value), nrow(cent_data), )
      cent_data$fname <- rep(fname, nrow(cent_data))
      cent_data$chromosome <- rep(chromosome, nrow(cent_data))

      utils::write.table(cent_data, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Centromere_MIN_Shannon_no_telo.txt"),
                         append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)
      #--------------------
      # write out Shannon diversity data
      Shannon_div1 <- cbind(genome_pos, Shannon_div)

      utils::write.table(Shannon_div1, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Shannon_div_no_telo.txt"),
                         append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

    }
  }
}

#' run_diversity_plots_35_2000_no_telomere
#'
#' Creates plots of diversity indices (Shannon, Simpson, Fisher's alpha and Pielou's J).
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_diversity_plots_35_2000_no_telomere <- function(chromosome=chromosome, fname=fname, inpath=inpath, outpath=outpath){
  # read in total All_spec_merged

  if(base::file.exists(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"))){

    # read in file
    base::cat("\nfile exists: ", paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"))
    All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

    # remove blank/zero columns in All_spec
    col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

    # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
    All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
    All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

    # cut telomeres - need to cut exactly one bin though
    telo_cut <- round(ncol(All_spec0)*0.04)
    All_spec2 <- All_spec1[,-c(1:telo_cut,(ncol(All_spec1)-telo_cut):ncol(All_spec1))]
    # cut last bin (Xbp from end ) from histograms ****

    col_length2 <- as.numeric(ncol(All_spec2))*5000
    numbins <-  floor(col_length2/2e+6)
    total_col <- (numbins*2e+6)/5000
    All_spec <- All_spec2#[,c(1:total_col)]

    full_length <-  as.numeric(ncol(All_spec))*5000

    full_length_value <-  as.numeric(ncol(All_spec0))*5000

    #-----------------------------
    # calculate species diversity
    #install.packages("vegan")
    #library(vegan)

    genome_pos <- colnames(All_spec)
    All_spec_repeatlengths <- rownames(All_spec)

    # try normalize the rows, for each repeat length, make a fraction of total
    # then run Shannon
    All_spec <- as.matrix(All_spec)
    All_spec_replen_sum <- apply(All_spec, 1, sum)
    All_spec_norm <- All_spec/All_spec_replen_sum
    All_spec_norm_t0 <- t(as.matrix(All_spec_norm))

    # subset to only repeats >35bp
    All_spec_norm_t <- All_spec_norm_t0[,-which(as.numeric(colnames(All_spec_norm_t0))<35)]

    # diversity calculation
    # https://rdrr.io/cran/vegan/man/diversity.html

    Shannon_div <- vegan::diversity(All_spec_norm_t, index = "shannon")#this is the Shannon-Wiener index
    Pielou_div <- vegan::diversity(All_spec_norm_t, index = "shannon")/log(vegan::specnumber(All_spec_norm_t))
    Simpson_div <- vegan::diversity(All_spec_norm_t, index = "simpson")#this is the Simpson index

    Shannon_div <- as.numeric(Shannon_div)
    Pielou_div <- as.numeric(Pielou_div)
    Simpson_div <- as.numeric(Simpson_div)

    Shannon_div[which(Shannon_div<=2)] <- NA
    Pielou_div[which(Pielou_div<=0.3)] <- NA
    Simpson_div[which(Simpson_div<=0.3)] <- NA

    #---------------------
    # make windows and run Shannon on those windows

    rownames(All_spec_norm_t)
    colnames(All_spec_norm_t)
    nrow(All_spec_norm_t)

    wind_factor <- round(full_length/150e3)
    Group <- vctrs::vec_rep_each(c(1:floor(nrow(All_spec_norm_t)/wind_factor)), wind_factor)
    Group_plus <- rep(max(Group)+1, nrow(All_spec_norm_t) - length(Group))
    Group_total <- c(Group, Group_plus)
    Shannon_div_wind <- vegan::diversity(All_spec_norm_t, groups=Group_total , index = "shannon")
    Shannon_div_wind <- as.numeric(Shannon_div_wind)
    Shannon_div_wind[which(Shannon_div<=2)] <- NA

    #---------------------
    # plots

    genome_pos_wind <- unique(Group_total)*wind_factor*5000
    Shannon_div_wind_pos <- cbind(genome_pos_wind, Shannon_div_wind)
    cent_wind <- Shannon_div_wind_pos[which(Shannon_div_wind_pos[,2] == min(Shannon_div_wind_pos[,2])),1]

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window_35_no_telo_", wind_factor,".png"), width = 3000, height = 1000)
    plot(genome_pos_wind, Shannon_div_wind, type="l")
    grDevices::dev.off()

    grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm_35_no_telo.png"), width = 3000, height = 1000)
    plot(genome_pos, Shannon_div)
    grDevices::dev.off()


    #-----------------------
    # moving average across Shannon that shows where min is
    # https://www.storybench.org/how-to-calculate-a-rolling-average-in-r/

    #---------------
    bin_size=25
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent25 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #----------
    bin_size=100
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent100 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))
    #----------
    bin_size=250
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent250 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #----------
    bin_size=500
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent500 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #----------
    bin_size=1000
    roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

    cent1000 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

    #-------------------
    # write out the centromere position
    cent_data0 <- rbind(cent25, cent100, cent250, cent500, cent1000)
    cent_data <- cent_data0*5000+7501
    cent_data[6] <- cent_wind
    cent_data <- as.data.frame(cent_data)

    rownames(cent_data) <- c(paste0(fname, "_", chromosome, "_cent25"),paste0(fname, "_", chromosome, "_cent100"), paste0(fname, "_", chromosome, "_cent250"), paste0(fname, "_", chromosome, "_cent500"), paste0(fname, "_", chromosome, "_cent1000"), paste0(fname, "_", chromosome,"_", wind_factor, "_centwind"))

    cent_data$full_length <- rep(as.character(full_length_value), nrow(cent_data), )
    cent_data$fname <- rep(fname, nrow(cent_data))
    cent_data$chromosome <- rep(chromosome, nrow(cent_data))

    utils::write.table(cent_data, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Centromere_MIN_Shannon_35_no_telo.txt"),
                       append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)

    #--------------------
    # write out Shannon diversity data
    Shannon_div1 <- cbind(genome_pos, Shannon_div)

    utils::write.table(Shannon_div1, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Shannon_div_35_no_telo.txt"),
                       append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

  } else {
    for (i in 1:sum(grepl(chromosome, nam_list1))){
      Spp_chr_part <- nam_list1[grepl(chromosome, nam_list1)][i]
      chr_part <- stringr::str_split(Spp_chr_part, "_", simplify =TRUE)[,3]

      out_chromosome <- paste0(outpath, "/", fname,"/", chr_part)
      All_spec_file <- paste0(out_chromosome,"/", Spp_chr_part,"_All_spec.txt")

      if(base::file.exists(All_spec_file)){
        base::cat("\nfile exists: ", All_spec_file)
        All_spec0<-base::as.matrix(utils::read.table(All_spec_file, header = TRUE, check.names = FALSE))
      } else {

        nam_list0a <- base::list.files(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/"))
        nam_list1a <- nam_list0a[base::grep("_spar1_Table.txt", nam_list0a)]
        if (length(base::grep("merged", nam_list1a))!=0){
          nam_list1a <- nam_list1a[-base::grep("merged", nam_list1a)]
        }
        for (i in 1:base::length(nam_list1a)){
          file_nam <- nam_list1a[i]
          All_spec_merge_file <- paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/merged_", file_nam)
          if(base::file.exists(All_spec_merge_file)){
            base::cat("\nfile exists: ", All_spec_merge_file)
          } else {
            base::cat(file_nam)
            All_spec <- base::as.matrix(utils::read.table(paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/", file_nam), check.names = FALSE))

            All_spec_merged <-NULL
            All_spec_merged1 <-NULL
            All_spec_Freq <-NULL

            for (col in 1:base::ncol(All_spec)){
              base::cat(col, " \n")
              Freq_spec <- NULL
              Freq_spec <- base::cbind(Freq_spec, All_spec[,c(1,col)])
              Freq_spec[,1] <- base::row.names(Freq_spec)
              Freq_spec <- base::as.data.frame(Freq_spec)
              Freq_spec[,3] <-  base::sub("1/", "", Freq_spec[,1] , fixed = TRUE)
              Freq_spec[,4] <- base::as.factor(round(as.numeric(Freq_spec[,3])))
              colnames(Freq_spec) <- base::c("Freq", "Power", "1/Freq", "Freq_round")
              Freq_spec_by_Freq_round <- dplyr::group_by(Freq_spec, by=Freq_round)
              Freq_spec_by_Freq_round$Power <- as.numeric(Freq_spec_by_Freq_round$Power)
              Freq_spec_merged <- dplyr::summarise(Freq_spec_by_Freq_round, Power_merged = mean(Power, na.rm=TRUE))
              Freq_spec_merged <- Freq_spec_merged[base::order(Freq_spec_merged$by, decreasing = TRUE),]
              All_spec_merged <- base::cbind(All_spec_merged, Freq_spec_merged$Power_merged)
              All_spec_Freq <- base::as.character(Freq_spec_merged$by)
            }

            All_spec_merged1 <- base::cbind(All_spec_Freq, All_spec_merged)

            # save the merged version of All_spec
            utils::write.table(All_spec_merged1, paste0(outpath, "/", fname,"/", chr_part,"/spectra_Table.txt/merged_", file_nam))
          }
        }

        # build chromosome scale spectra
        All_spec_total_file <- paste0(outpath,"/", fname,"/",chr_part,"/Total_",fname, "_", chr_part,"_All_spec_merged.txt")
        All_spec_Total <- NULL
        if(base::file.exists(All_spec_total_file)){
          base::cat("\nfile exists: ", All_spec_total_file)
        } else {
          nam_list0a <- base::list.files(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/"))
          nam_list5 <- nam_list0a[base::grep("merged", nam_list0a)]
          for (i in 1:base::length(nam_list5)){
            nam <- nam_list5[i]
            cat("Allspec ", nam)
            # join All_spec files for a given chromosome
            All_spec_tmp0 <- base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/", chr_part,"/", "spectra_Table.txt/", nam), header = TRUE, check.names = FALSE))
            All_spec_tmp1 <- All_spec_tmp0[,-1]
            All_spec_Total <-base::cbind(All_spec_Total, All_spec_tmp1)
            row.names(All_spec_Total) <- as.character(All_spec_tmp0[,1])
          }
          colnames(All_spec_Total) <-  as.character(seq.int(2501, ((ncol(All_spec_Total)-1)*5000)+2501, 5000))
          utils::write.table(All_spec_Total, All_spec_total_file)

          gc()
        }
      }

      All_spec0 <-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chr_part,"/Total_",fname, "_", chr_part,"_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

      # remove blank/zero columns in All_spec
      col_length0 <- as.numeric(colnames(All_spec0)[ncol(All_spec0)])

      # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
      All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
      All_spec1 <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

      # cut telomeres - need to cut exactly one bin though
      telo_cut <- round(ncol(All_spec0)*0.04)
      All_spec2 <- All_spec1[,-c(1:telo_cut,(ncol(All_spec1)-telo_cut):ncol(All_spec1))]
      # cut last bin (Xbp from end ) from histograms ****

      col_length2 <- as.numeric(ncol(All_spec2))*5000
      numbins <-  floor(col_length2/2e+6)
      total_col <- (numbins*2e+6)/5000
      All_spec <- All_spec2#[,c(1:total_col)]

      full_length <-  as.numeric(ncol(All_spec))*5000

      full_length_value <-  as.numeric(ncol(All_spec0))*5000

      #-----------------------------
      # calculate species diversity
      #install.packages("vegan")
      #library(vegan)
      genome_pos <- colnames(All_spec)
      All_spec_repeatlengths <- rownames(All_spec)

      # try normalize the rows, for each repeat length, make a fraction of total
      # then run Shannon

      All_spec <- as.matrix(All_spec)
      All_spec_replen_sum <- apply(All_spec, 1, sum)
      All_spec_norm <- All_spec/All_spec_replen_sum
      All_spec_norm_t0 <- t(as.matrix(All_spec_norm))

      # subset to only repeats >35bp
      All_spec_norm_t <- All_spec_norm_t0[,-which(as.numeric(colnames(All_spec_norm_t0))<35)]

      # diversity calculation

      Shannon_div <- vegan::diversity(All_spec_norm_t, index = "shannon")#this is the Shannon-Wiener index
      Pielou_div <- vegan::diversity(All_spec_norm_t, index = "shannon")/log(vegan::specnumber(All_spec_norm_t))
      Simpson_div <- vegan::diversity(All_spec_norm_t, index = "simpson")#this is the Simpson index

      Shannon_div <- as.numeric(Shannon_div)
      Pielou_div <- as.numeric(Pielou_div)
      Simpson_div <- as.numeric(Simpson_div)

      Shannon_div[which(Shannon_div<=2)] <- NA
      Pielou_div[which(Pielou_div<=0.3)] <- NA
      Simpson_div[which(Simpson_div<=0.3)] <- NA

      #---------------------
      # make windows and run Shannon on those windows

      rownames(All_spec_norm_t)
      colnames(All_spec_norm_t)
      nrow(All_spec_norm_t)

      wind_factor <- round(full_length/150000)
      Group <- vctrs::vec_rep_each(c(1:floor(nrow(All_spec_norm_t)/wind_factor)), wind_factor)
      Group_plus <- rep(max(Group)+1, nrow(All_spec_norm_t) - length(Group))
      Group_total <- c(Group, Group_plus)
      Shannon_div_wind <- vegan::diversity(All_spec_norm_t, groups=Group_total , index = "shannon")
      Shannon_div_wind <- as.numeric(Shannon_div_wind)
      Shannon_div_wind[which(Shannon_div<=2)] <- NA

      #---------------------
      # plots

      genome_pos_wind <- unique(Group_total)*wind_factor*5000
      Shannon_div_wind_pos <- cbind(genome_pos_wind, Shannon_div_wind)
      cent_wind <- Shannon_div_wind_pos[which(Shannon_div_wind_pos[,2] == min(Shannon_div_wind_pos[,2])),1]

      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_div_window_35_no_telo_", wind_factor,".png"), width = 3000, height = 1000)
      plot(genome_pos_wind, Shannon_div_wind, type="l")
      grDevices::dev.off()


      grDevices::png(file=paste0(outpath,"/", fname,"/",chromosome,"/",fname, "_", chromosome, "_Shannon_plot_norm_35_no_telo.png"), width = 1500, height = 500)
      plot(genome_pos, Shannon_div)
      grDevices::dev.off()


      #-----------------------
      # moving average across Shannon that shows where min is
      # https://www.storybench.org/how-to-calculate-a-rolling-average-in-r/

      #---------------
      bin_size=25
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent25 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      #----------
      bin_size=100
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent100 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      #----------
      bin_size=250
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent250 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      #----------
      bin_size=500
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent500 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      #----------
      bin_size=1000
      roll_mean_Shannon <- zoo::rollapply(Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

      cent1000 <- base::which(roll_mean_Shannon == base::min(roll_mean_Shannon, na.rm=TRUE))

      #-------------------
      # write out the centromere position
      cent_data0 <- rbind(cent25, cent100, cent250, cent500, cent1000)
      cent_data <- cent_data0*5000+7501
      cent_data[6] <- cent_wind
      cent_data <- as.data.frame(cent_data)

      rownames(cent_data) <- c(paste0(fname, "_", chromosome, "_cent25"),paste0(fname, "_", chromosome, "_cent100"), paste0(fname, "_", chromosome, "_cent250"), paste0(fname, "_", chromosome, "_cent500"), paste0(fname, "_", chromosome, "_cent1000"), paste0(fname, "_", chromosome,"_", wind_factor, "_centwind"))

      cent_data$full_length <- rep(as.character(full_length_value), nrow(cent_data), )
      cent_data$fname <- rep(fname, nrow(cent_data))
      cent_data$chromosome <- rep(chromosome, nrow(cent_data))

      utils::write.table(cent_data, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Centromere_MIN_Shannon_35_no_telo.txt"),
                         append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = TRUE, col.names = FALSE)
      #--------------------
      # write out Shannon diversity data
      Shannon_div1 <- cbind(genome_pos, Shannon_div)

      utils::write.table(Shannon_div1, file=paste0(outpath, "/", fname, "/",chromosome, "/", fname, "_", chromosome, "_Shannon_div_35_no_telo.txt"),
                         append = FALSE, sep = " ", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

    }
  }
}



#' run_summary_plots_range
#'
#' Creates spectra and DNAwalks plots for a specific range. Needs to be run after the spectra and DNAwalk have all been built.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_summary_plots_range <- function(fname=fname, chromosome=chromosome, inpath=inpath, outpath=outpath, startbp=startbp, endbp=endbp){

  # in part1
  if (startbp<1e8 && endbp<1e8) {
    part="part01"
    # reset start and end positions for this part
    startbp = startbp
    endbp = endbp
  }
  # in part2
  if (startbp>1e8 && endbp>1e8 && startbp<2e8 && endbp<2e8) {
    part="part02"
    # reset start and end positions for this part
    startbp = startbp - 1e8
    endbp = endbp - 1e8
  }
  # in part3
  if (startbp>2e8 && endbp>2e8 && startbp<3e8 && endbp<3e8) {
    part="part03"
    # reset start and end positions for this part
    startbp = startbp - 2e8
    endbp = endbp - 2e8
  }
  # in part4
  if (startbp>3e8 && endbp>3e8 && startbp<4e8 && endbp<4e8) {
    part="part04"
    # reset start and end positions for this part
    startbp = startbp - 3e8
    endbp = endbp - 3e8
  }
  #in part5
  if (startbp>4e8 && endbp>4e8 && startbp<5e8 && endbp<5e8) {
    part="part05"
    # reset start and end positions for this part
    startbp = startbp - 4e8
    endbp = endbp - 4e8
  }
  #in part6
  if (startbp>5e8 && endbp>5e8 && startbp<6e8 && endbp<6e8) {
    part="part06"
    # reset start and end positions for this part
    startbp = startbp - 5e8
    endbp = endbp - 5e8
  }
  #in part7
  if (startbp>6e8 && endbp>6e8 && startbp<7e8 && endbp<7e8) {
    part="part07"
    # reset start and end positions for this part
    startbp = startbp - 6e8
    endbp = endbp - 6e8
  }
  #in part8
  if (startbp>7e8 && endbp>7e8 && startbp<8e8 && endbp<8e8) {
    part="part08"
    # reset start and end positions for this part
    startbp = startbp - 7e8
    endbp = endbp - 7e8
  }
  #in part9
  if (startbp>8e8 && endbp>8e8 && startbp<9e8 && endbp<9e8) {
    part="part09"
    # reset start and end positions for this part
    startbp = startbp - 8e8
    endbp = endbp - 8e8
  }
  #---------------------
  if (startbp<=5e+06 && endbp<=5e+06){
      startval="1"
      endval="5e+06"
  }
  if (startbp>=5000001 && endbp>=5000001 && startbp<=1e+07 && endbp<=1e+07){
    startval="5000001"
    endval="1e+07"
  }
  if (startbp>10000001 && endbp>10000001 && startbp<1.5e+07 && endbp<1.5e+07){
    startval="10000001"
    endval="1.5e+07"
  }
  if (startbp>15000001 && endbp>15000001 && startbp<2e+07 && endbp<2e+07){
    startval="15000001"
    endval="2e+07"
  }
  if (startbp>20000001 && endbp>20000001 && startbp<2.5e+07 && endbp<2.5e+07){
    startval="20000001"
    endval="2.5e+07"
  }

  if (startbp>15000001 && endbp>15000001 && startbp<2e+07 && endbp<2e+07){
    startval="15000001"
    endval="2e+07"
  }
  if (startbp>20000001 && endbp>20000001 && startbp<2.5e+07 && endbp<2.5e+07){
    startval="20000001"
    endval="2.5e+07"
  }

  if (startbp>25000001 && endbp>25000001 && startbp<3e+07 && endbp<3e+07){
    startval="25000001"
    endval="3e+07"
  }
  if (startbp>30000001 && endbp>30000001 && startbp<3.5e+07 && endbp<3.5e+07){
    startval="30000001"
    endval="3.5e+07"
  }

  if (startbp>35000001 && endbp>35000001 && startbp<4e+07 && endbp<4e+07){
    startval="35000001"
    endval="4e+07"
  }
  if (startbp>40000001 && endbp>40000001 && startbp<5.5e+07 && endbp<5.5e+07){
    startval="40000001"
    endval="5.5e+07"
  }

  if (startbp>45000001 && endbp>45000001 && startbp<5e+07 && endbp<5e+07){
    startval="45000001"
    endval="5e+07"
  }
  if (startbp>50000001 && endbp>50000001 && startbp<5.5e+07 && endbp<5.5e+07){
    startval="50000001"
    endval="5.5e+07"
  }

  if (startbp>55000001 && endbp>55000001 && startbp<6e+07 && endbp<6e+07){
    startval="55000001"
    endval="6e+07"
  }
  if (startbp>60000001 && endbp>60000001 && startbp<6.5e+07 && endbp<6.5e+07){
    startval="60000001"
    endval="6.5e+07"
  }

  if (startbp>65000001 && endbp>65000001 && startbp<7e+07 && endbp<7e+07){
    startval="65000001"
    endval="7e+07"
  }
  if (startbp>70000001 && endbp>70000001 && startbp<7.5e+07 && endbp<7.5e+07){
    startval="70000001"
    endval="7.5e+07"
  }

  if (startbp>75000001 && endbp>75000001 && startbp<8e+07 && endbp<8e+07){
    startval="75000001"
    endval="8e+07"
  }
  if (startbp>80000001 && endbp>80000001 && startbp<8.5e+07 && endbp<8.5e+07){
    startval="80000001"
    endval="8.5e+07"
  }

  if (startbp>85000001 && endbp>85000001 && startbp<9e+07 && endbp<9e+07){
    startval="85000001"
    endval="9e+07"
  }
  if (startbp>90000001 && endbp>90000001 && startbp<9.5e+07 && endbp<9.5e+07){
    startval="90000001"
    endval="9.5e+07"
  }

  if (startbp>95000001 && endbp>95000001 && startbp<1e+08 && endbp<1e+08){
    startval="95000001"
    endval="1e+08"
  }

  #----------------------------
  # read in the relevant All_spec file
  All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome, part,"/spectra_Table.txt/concat_ALL_tot__", fname, "_", chromosome, part, "_", startval,"_", endval, "_25Kb_5000_spar1_Table.txt"), check.names = FALSE))
  # concat_ALL_tot__Ha412_H0_Chr05part01_95000001_1e+08_25Kb_5000_spar1_Table.txt

  # read in the relevant DNAwalk
  DNAwalk_long <-base::as.matrix(utils::read.table(paste0(outpath, "/", fname,"/", chromosome, part,"/DNAWALK/FULLDNA_WALK__",fname, "_", chromosome, part,  "_", startval,"_", endval,"_25Kb_5000_spar1.txt"), check.names = FALSE))
  # FULLDNA_WALK__Ha412_H0_Chr05part01_10000001_1.5e+07_25Kb_5000_spar1.txt

  #--------------------------
  # remove blank/zero columns in All_spec
  colnames(All_spec0)[ncol(All_spec0)]

  # https://stackoverflow.com/questions/21530168/remove-columns-with-zero-values-from-a-dataframe
  All_spec0[,which(colSums(All_spec0, na.rm = FALSE, dims = 1) < 1e-5)] <- NA
  All_spec <- All_spec0[, !is.na(colSums(All_spec0 != 0)) & colSums(All_spec0 != 0) > 0]

  startpos <- round((startbp-as.numeric(startval)+2500)/5000+1)
  endpos <- round((endbp-as.numeric(startval)-2500)/5000)
  All_spec_subset <- All_spec[,c(startpos:endpos)]

  pngflag=TRUE

  ofi <- paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"All_spec_subset_", startbp,"-",endbp)

  #r0=base::c(2,35)
  r1=base::c(35,2000)

  base::cat("\nofi:",ofi,"\n")
  #if(!base::is.null(r0))largeimagesub_NEW(All_spec_subset,fname,chromosome,ofi,rangebp=r0,pngflag = pngflag)
  if(!base::is.null(r1))largeimagesub_NEW(All_spec_subset,fname,chromosome,ofi,rangebp=r1,pngflag = pngflag, pdf_flag=TRUE)

  #------------------------------------------
  # plot total DNAwalk for subset chromosome

  startpos1 <- round(startbp-as.numeric(startval))
  endpos1 <- round(endbp-as.numeric(startval))

  dnawalk <- as.data.frame(DNAwalk_long[c(startpos1:endpos1),])

  walk_colours <- NULL
  for (col in grDevices::rainbow(n=round(nrow(dnawalk)/5000))){
    coltmp <- rep(col, 5000)
    walk_colours <- c(walk_colours, coltmp)
  }

  #------------------
  # write to image
  ofile2 <-  paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome, "_DNAwalk2D_subset_", startbp,"-",endbp,".png")
  base::cat("\n ouput to", ofile2)
  grDevices::png(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome, "_DNAwalk2D_subset_", startbp,"-",endbp,".png"))
  plot(dnawalk$AT, dnawalk$CG, col = walk_colours)
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome, "_DNAwalk2D_subset_", startbp,"-",endbp,".pdf"))
  plot(dnawalk$AT, dnawalk$CG, col = walk_colours)
  grDevices::dev.off()

  #----------
  ofile3 <-  paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_AT_subset_", startbp,"-",endbp,".png")
  base::cat("\n ouput to", ofile3)
  grDevices::png(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_AT_subset_", startbp,"-",endbp,".png"))
  # plot 1D DNAwalk
  dnawalk$rownames <- c(1:nrow(dnawalk))
  plot(dnawalk$rownames, dnawalk$AT, col = walk_colours)
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_AT_subset_", startbp,"-",endbp,".pdf"))
  # plot 1D DNAwalk
  dnawalk$rownames <- c(1:nrow(dnawalk))
  plot(dnawalk$rownames, dnawalk$AT, col = walk_colours)
  grDevices::dev.off()

  #---------------------
  ofile4 <-  paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_CG_subset_", startbp,"-",endbp,".png")
  base::cat("\n ouput to", ofile4)
  grDevices::png(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_CG_subset_", startbp,"-",endbp,".png"))
  # plot 1D DNAwalk
  plot(dnawalk$rownames, dnawalk$CG, col = walk_colours)
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/", chromosome,"/", fname, "_", chromosome,"_DNAwalk1D_CG_subset_", startbp,"-",endbp,".pdf"))
  # plot 1D DNAwalk
  plot(dnawalk$rownames, dnawalk$CG, col = walk_colours)
  grDevices::dev.off()
}


#' run_20Mbpimag
#'
#' Creates 20Mbp images of the spectra. Needs to be run after the spectra have all been built.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_20Mbpimag <- function(nam=nam, fname=fname, inpath=inpath,
                          outpath=outpath, atflag=TRUE, pflag=FALSE, plotflag=FALSE,
                          writeflag=FALSE){

  in_name<- base::paste0(inpath,nam,".fasta")
  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")
  dn1<-dna_1_to_wax(dna_1)
  full_length<-base::length(dn1)

  #this runs CG contant
  winbox<-1*1e6
  #walklist<- CG_AT_content(dn1,nam,winbox=winbox)

  dn1<-NULL; dna_1<-NULL;base::gc()
  max_size<-2.0*1e6
  samplesize<-base::round(full_length/max_size)     #samplesize for long spectrum fft windows 150/2=75 100/2=50  50/2=25 20/2=10
  samplesize<-base::floor(samplesize/100)*100
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /50)*50   #full_length<-194*1e6 full_length<-210*1e6
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /20)*20
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /10)*10   #full_length<-6*1e6 full_length<-3*1e6
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /5)*5
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /2)*2     #full_length<-1.5*1e6 full_length<-3*1e6
  if(samplesize==0) samplesize<-1
  base::cat("\nSamplesize is", samplesize,"\n")

  twolists<-run_chromosomelistNEW(nam, fname,inpath, outpath)
  All_list<-twolists$All_list
  fullwalk_list<-twolists$fullwalklist
  #All_spec<-read_All_list(All_list)
  All_spec<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"_", nam,"_All_spec.txt"), check.names = FALSE))

  whival<-(base::which(base::as.numeric(base::colnames(All_spec))>full_length+2500))
  if(base::length(whival)>0){
    All_spec[1:10,(whival)]
    All_spec<-All_spec[,-(whival)]
  }
  base::gc()

  # 20 Mbp images of spectra
  increm<-20*1e+6
  snum<-2501
  for (N in base::seq(snum,full_length,by=increm)){ #
    enum<-N+increm
    if(enum>full_length)enum<- base::colnames(All_spec)[ base::ncol(All_spec)]
    run_largeimage(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath, atflag=atflag,
                   rangeseq=base::c(base::as.character(N),base::as.character(enum)),inmain="Individual 5Kbp window")
    base::gc()
    run_largeimage(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath, atflag=atflag,
                   r0=base::c(5,200),r1=NULL,r2=NULL,
                   rangeseq=base::c(base::as.character(N),base::as.character(enum)),inmain="Individual 5Kbp window")
  }
  base::gc()
}


#' runs_run_barplots
#'
#' Creates barplots of the spectra and DNAwalks. Needs to be run after the spectra have all been built.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
runs_run_barplots <- function(nam=nam, fname=fname, inpath=inpath, outpath=outpath, atflag=TRUE,
                              pflag=FALSE, plotflag=FALSE, writeflag=FALSE){

  in_name<- base::paste0(inpath,nam,".fasta")    #

  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")   #nchar(tdna) base::nchar(waxy1)= 295045
  dn1<-dna_1_to_wax(dna_1)   #base::length(dn) base::nchar(dn)
  full_length<-base::length(dn1)

  #this runs CG contant
  winbox<-1*1e6

  dn1<-NULL; dna_1<-NULL;base::gc()

  #requires fft window to be multiple of 10 for dnawalk in run_chloroplast
  max_size<-2.0*1e6
  samplesize<-base::round(full_length/max_size)
  samplesize<-base::floor(samplesize/100)*100
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /50)*50
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /20)*20
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /10)*10
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /5)*5
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /2)*2
  if(samplesize==0) samplesize<-1
  base::cat("\nSamplesize is", samplesize,"\n")

  twolists<-run_chromosomelistNEW(nam, fname,inpath, outpath)
  All_list<-twolists$All_list
  fullwalk_list<-twolists$fullwalklist
  All_spec<-read_All_list(All_list)
  #utils::write.table(All_spec, paste0(outpath,"/", fname,"_", nam,"_All_spec.txt"))
  #All_spec<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname,"/",chromosome,"/Total_",fname, "_", chromosome,"_All_spec.txt"), check.names = FALSE))

  whival<-(base::which(base::as.numeric(base::colnames(All_spec))>full_length+2500))
  if(base::length(whival)>0){
    All_spec[1:10,(whival)]
    All_spec<-All_spec[,-(whival)]
  }
  base::gc()

  # Barplots
  if(full_length>(50*1e6))numbins<-25 else numbins<-20
  binsize<-base::round(full_length/(numbins*1e+6))
  fftlength<-5000
  numrange<-base::round(1e6*binsize/fftlength)
  if(full_length>(50*1e6)){numstd_bins<-1;numstdred<-3 } else { numstd_bins<-1;numstdred<-2 }
  #run_barplots(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath,numstd=1,  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(35,2000))
  run_barplots(All_spec=All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath,atflag=atflag, numstd=numstdred,  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(35,2000))
  #run_barplots(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath,numstd=5,  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(35,2000))
  #run_barplots(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath,numstd=7,  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(35,2000))
  #run_barplots(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath,numstd=(numstdred*2),  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(15,35))
  #run_barplots(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath,numstd=(numstdred*2),  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(8,15))
  #run_barplots(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath,numstd=(numstdred*2),  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(3.5,8))
  #run_barplots(All_spec,nam=nam,fname=fname,inpath=inpath, outpath=outpath,numstd=(numstdred*2),  numstd2=numstd_bins, numrange=numrange,repeat_range=base::c(2,3.5))
  All_spec<-NULL
  base::gc()

}


#' run_long_repeats_NEW
#'
#' # Runs large fft windows to look for longer repeats
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
run_long_repeats_NEW <- function(nam=nam, fname=fname, inpath=inpath,
                                 outpath=outpath, AT_flag=TRUE, atflag=TRUE, pflag=FALSE, plotflag=FALSE,
                                 writeflag=FALSE){

  in_name<- base::paste0(inpath,nam,".fasta")    #
  dna_1<-utils::read.csv(in_name, header = FALSE, sep = ",", dec = ".")

  # checks chromosome and gets full length
  dn1<-dna_1_to_wax(dna_1)   #base::length(dn) base::nchar(dn)
  full_length<-base::length(dn1)
  dn1<-NULL; dna_1<-NULL;base::gc()

  splitname<-base::unlist(stringr::str_split(nam,"_"))   #fname<-fname
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]

  #requires fft window to be multiple of 10 for dnawalk in run_chloroplast
  max_size<-2.0*1e6
  samplesize<-base::round(full_length/max_size)
  samplesize<-base::floor(samplesize/100)*100
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /50)*50
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /20)*20
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /10)*10
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /5)*5
  if(samplesize==0) samplesize<-base::floor(base::round(full_length/max_size) /2)*2
  if(samplesize==0) samplesize<-1
  base::cat("\nSamplesize is", samplesize,"\n")

  bigincre<-100000000;
  samplesizesmall=10
  fftlength<-fftlength1<-50000 ; fftlengthlarge<-100000; fftlengthlargest<-200000
  All_spec_long<-NULL
  DNAwalk_long<-NULL   #    DNAwalk_long<-DNAwalk[5:6,]   tail(base::rownames(All_walk))
  fractal_long<-NULL;fractal_shorter<-NULL
  bprun<-base::seq(1,full_length,bigincre)

  for(sbp in bprun){
    base::gc()
    ebp<-sbp+bigincre-1  #change Nov 25
    if(ebp>full_length) ebp<-full_length

    outpathwalk<-base::file.path(outpath, nam1,chr, base::paste("DNAWALK"))
    walk_txtfile<-base::file.path(outpathwalk, base::paste0("dnawalk_",sbp,"_",ebp,"_",nam,"Table.txt"))

    if(base::file.exists(walk_txtfile)){
      base::cat("\nfile exists: reading ",walk_txtfile)
      All_walk<-utils::read.table( walk_txtfile ,header=TRUE)
      vals<-base::which(base::as.numeric(base::rownames(All_walk))>ebp);if(base::length(vals)>0) All_walk<-All_walk[-vals,]
    } else {
      base::cat("\nfile does not exist: construct and write ",walk_txtfile)  #sbp<-1; ebp<-20000000
      twolists1<-run_chromosomelistNEW(nam=nam, fname=fname, inpath=inpath, outpath=outpath, majorlen=5000000,
                                       length_majorgroup=1000000, name_majorgroup="1Mbp",
                                       startgroup=sbp,endgroup=ebp,full_length=full_length)
      walk_listChr<-twolists1$fullwalklist   # twolists$fullwalklist[base::length(twolists$fullwalklist)]
      All_walk<-read_walk_list(walk_listChr)   #twolists$fullwalklist
      vals<-base::which(base::as.numeric(base::rownames(All_walk))>ebp);if(base::length(vals)>0) All_walk<-All_walk[-vals,]
      #All_walk[base::nrow(All_walk),]
      utils::write.table(All_walk, walk_txtfile)    #tail(walk_listChr)
    }

    bp<-base::c(base::as.numeric(base::rownames(All_walk)[1]),base::as.numeric(base::rownames(All_walk)[base::nrow(All_walk)])) #base::as.character(se)[1:10])
    se<-base::seq( bp[1],bp[2] ,by=50)   #different inversion in hetero third strand duplicates first part of palindrome

    if((bp[2]-bp[1])>=fftlengthlarge)fftlength10<-fftlengthlarge else fftlength10<-(bp[2]-bp[1])
    All_list<-run_chloroplast(nam=nam,fname=fname, inpath=inpath, outpath=outpath,
                              startval= bp[1] ,endval=bp[2], fftlength=fftlength10,
                              All_walk=All_walk,sample_every=samplesize,
                              walkflag=FALSE,atflag=atflag,AT_flag=AT_flag, waxflag=FALSE,printlong=FALSE,
                              printdna=FALSE)

    All_spec_sub<-All_list[[1]]
    spectwalklist<-All_list[[2]]
    DNAwalk<-All_list[[3]]  #DNAwalk[1:10,]   attributes(All_list) tail(DNAwalk)
    fractal<-All_list[[4]]    # fractal[1:10]

    run_largeimage(All_spec_sub,fname=fname,nam=nam, inpath=inpath, outpath=outpath,samplesize=samplesize,
                   r0=base::c(samplesize*10,50000),
                   r1=base::c(samplesize*4,samplesize*4+800),
                   r2=base::c(samplesize*3,samplesize*3+400),
                   #r3=base::c(samplesize*2,samplesize*2+200),
                   #r4=base::c(samplesize*2,5000),
                   inmain= base::paste("Individual ", fftlengthlarge/1000,"window sample",samplesize))

    #run subintervals now
    startbp<-base::as.numeric(base::rownames(All_walk)[1])
    small_incre<-2000000 ;  #fftlength<-100000   #base::rownames(All_walk_3)[1:10]
    endbp<-base::as.numeric(base::rownames(All_walk)[base::nrow(All_walk)])-1
    bprun1<-base::seq(startbp,endbp,by=small_incre)
    for(sbp1 in bprun1){ # sbp1<-startbp
      ebp1<-sbp1+small_incre-1   #change Nov 25
      if(ebp1>endbp)ebp1<-endbp
      bp2<-base::c(sbp1,ebp1)  # base::rownames(All_walk)[9*500000]
      base::cat("\nSubinterval runs: sbp1 ebp1",sbp1,ebp1,"\n")

      if((bp2[2]-bp2[1])>=fftlength1){
        All_listsmall<-run_chloroplast(nam=nam, fname=fname,inpath=inpath, outpath=outpath,
                                       startval= bp2[1] ,endval=	 bp2[2], fftlength=fftlength1,
                                       All_walk=All_walk,sample_every=samplesizesmall,
                                       walkflag=FALSE,atflag=atflag, AT_flag=AT_flag, waxflag=FALSE,printlong=FALSE,
                                       printdna=FALSE,pflag=FALSE,plotflag=plotflag,writeflag=writeflag)

        All_spec_small<-All_listsmall[[1]]
        spectwalklistsmall<-All_listsmall[[2]]
        DNAwalksmall<-All_listsmall[[3]]  #DNAwalk[1:10,]   attributes(All_list)
        fractalsmall<-All_listsmall[[4]]    #fractalsmall[1:10]   fracdim_list<- fractalsmall

        if(plotflag) run_largeimage(All_spec=All_spec_small,nam=nam, fname=fname,inpath=inpath, outpath=outpath, samplesize=samplesizesmall,
                                    r0=base::c(700,10000),r1=base::c(200,7000),r2=base::c(1000,20000),#r3=base::c(1000,50000),r4=base::c(100,5000),
                                    inmain= base::paste0("Individual ",fftlength/1000,"Kbp window"))
        fractal_shorter<-base::c(fractal_shorter,fractalsmall)
      } else base::cat("\nfftlength1 larger than range of sequence",fftlength1,bp2[1],bp2[2],"\n")
      All_spec_small<-NULL; DNAwalksmall<-NULL;fractalsmall<-NULL;base::gc()
    }

    if(!base::is.null(DNAwalk_long)){
      DNAwalk[,"CG"]<- DNAwalk_long[base::nrow(DNAwalk_long),"CG"]+DNAwalk[,"CG"]
      DNAwalk[,"AT"]<- DNAwalk_long[base::nrow(DNAwalk_long),"AT"]+DNAwalk[,"AT"]
      DNAwalk_long<- base::rbind(DNAwalk_long,DNAwalk) # base::nrow(DNAwalk_long)*50  head(base::rownames(DNAwalk_long)) full_length
      fractal_long<-base::c(fractal_long,fractal)
    } else {
      DNAwalk_long<-DNAwalk
      fractal_long<-fractal
    }

    All_list<-run_chloroplast(nam=nam,fname=fname, inpath=inpath, outpath=outpath,
                              startval= base::as.numeric(base::rownames(DNAwalk_long)[1]) ,
                              endval=	 base::as.numeric(base::rownames(DNAwalk_long)[base::nrow(DNAwalk_long)]),
                              fftlength=fftlengthlargest,
                              All_walk=DNAwalk_long,sample_every=(-samplesize),
                              walkflag=FALSE,atflag=atflag, AT_flag=AT_flag, waxflag=FALSE,printlong=FALSE,
                              printdna=FALSE,plotflag=plotflag,writeflag=writeflag)

    fractallarge<-All_list[[4]]

    if((bp[2]-bp[1])>=fftlengthlarge)All_spec_long<-base::cbind( All_spec_long , All_spec_sub)

  } #end of 20Mbp loop sbp<-6e+07 bprun<-base::seq(20000001,full_length,bigincre)

  many_list<-run_chloroplast(nam=nam,fname=fname, inpath=inpath, outpath=outpath,
                             startval= base::as.numeric(base::rownames(DNAwalk_long)[1]) ,
                             endval=	 base::as.numeric(base::rownames(DNAwalk_long)[base::nrow(DNAwalk_long)]), fftlength=fftlengthlargest,
                             All_walk=DNAwalk_long,sample_every=(-samplesize),
                             walkflag=FALSE,atflag=atflag, AT_flag=AT_flag,waxflag=FALSE,printlong=FALSE,
                             printdna=FALSE)
  #above uses longer fft window to get complexity

  fftlengthsmallest<-10000

  fftlen_list<-base::c(10, 25 ,100,300,500,1000,5000)  #fftlen_list<-base::c(250,500) fftlen_list<-base::c(5000)
  startval<- base::as.numeric(base::rownames(DNAwalk_long)[1])
  endval<-	 base::as.numeric(base::rownames(DNAwalk_long)[base::nrow(DNAwalk_long)])
  for(fftlen in fftlen_list){
    fftlen1<-fftlen*1e3
    if(fftlen1<(endval-startval)){
      many_list<-run_chloroplast(nam=nam, fname=fname, inpath=inpath, outpath=outpath,
                                 startval= base::as.numeric(base::rownames(DNAwalk_long)[1]) ,
                                 endval=	 base::as.numeric(base::rownames(DNAwalk_long)[base::nrow(DNAwalk_long)]), fftlength=fftlen1,
                                 All_walk=DNAwalk_long,sample_every=(-samplesize),
                                 walkflag=FALSE,atflag=atflag, AT_flag=AT_flag,waxflag=FALSE,printlong=FALSE,
                                 printdna=FALSE)
    }
  }
  All_spec_long5Mbp<-many_list[[1]]

  fftlengthsmallest<-10000

  #creat fractal runs for whole chromosome using different sample
  splitname<-base::unlist(stringr::str_split(nam,"_"))   #fname<-"Cannabis"
  chr<-splitname[base::length(splitname)]
  nam1<-base::unlist(stringr::str_split(nam, base::paste0("_",chr)))[1]
  base::dir.create(base::file.path(outpath, nam1), showWarnings = FALSE)  #creat directory of individual with chromosomes
  base::dir.create(base::file.path(outpath, nam1,chr), showWarnings = FALSE)   # create Chromosome directory
  outpathfractal<-base::file.path(outpath, nam1,chr,"fractalD.other")
  base::dir.create(outpathfractal, showWarnings = FALSE)
  #use 100Kbp fft to gt complexity


  minf<-base::as.numeric(names(fractal_long)[1])      # base::as.numeric(names(fractal_long)[1:10])
  maxf<-base::as.numeric(names(fractal_long)[base::length(fractal_long)]);
  seqlen<-fftlengthlarge;sample_every<-50

  numrange= base::round(((maxf-minf)/sample_every)/(15*seqlen/sample_every))

  ofi<-(base::file.path(outpathfractal, base::paste0("FRACTAL",fname,nam,"_range_",numrange,"_sample_",sample_every,"_fft_",
                                                     seqlen,"_spar1_Table.txt")))
  utils::write.table(fractal_long, ofi)
  ofile1<-(base::file.path(outpathfractal, base::paste0("FRACTAL",fname,nam,"_range_",numrange,"_sample_",sample_every,"_fft_",
                                                        seqlen,"_spar1.bedGraph")))
  base::sink(ofile1)
  for(j in 1:base::length(fractal_long)){
    delta<-base::as.numeric(names(fractal_long)[2])-base::as.numeric(names(fractal_long)[1])  #base::as.numeric(names(fractal)[2])-base::as.numeric(names(fractal)[1])
    base::cat( base::paste0("\n",nam,"\t",base::as.numeric(names(fractal_long)[j]),"\t",
                            base::as.numeric(names(fractal_long)[j])+delta,"\t",
                            fractal_long[j]) )
  }
  base::sink()

  ofi<-(base::file.path(outpathfractal, base::paste0("FRACTAL",fname,nam,"_range_",numrange,"_sample_",sample_every,"_fft_",
                                                     seqlen,"_spar1.pdf")))
  grDevices::pdf(ofi)
  ofi<-(base::file.path(outpathfractal, base::paste0("FRACTAL",fname,nam,"_range_",numrange,"_sample_",sample_every,"_fft_",
                                                     seqlen,"_spar1_info.txt")))
  base::sink(ofi)
  main<- base::paste(nam,"fft window",seqlen)    #fractal_long[1:10]
  run_sum_Fractal(fracdim_list=fractal_long, numrange= numrange,pflag=TRUE,main=main)
  grDevices::dev.off()
  base::sink()

  #use smaller ffts to get complexity
  minf<-base::as.numeric(names(fractal_shorter)[1])
  maxf<-base::as.numeric(names(fractal_shorter)[base::length(fractal_shorter)]);
  seqlen<-fftlength;sample_every<-samplesizesmall
  numrange= base::round(((maxf-minf)/sample_every)/(15*seqlen/sample_every))
  ofi<-(base::file.path(outpathfractal, base::paste0("FRACTAL",fname,nam,"_range_",numrange,"_sample_",sample_every,"_fft_",
                                                     seqlen,"_spar1_Table.txt")))
  utils::write.table(fractal_shorter, ofi)
  ofi<-(base::file.path(outpathfractal, base::paste0("FRACTAL",fname,nam,"_range_",numrange,"_sample_",sample_every,"_fft_",
                                                     seqlen,"_spar1.pdf")))
  ofile1<-(base::file.path(outpathfractal, base::paste0("FRACTAL",fname,nam,"_range_",numrange,"_sample_",sample_every,"_fft_",
                                                        seqlen,"_spar1.bedGraph")))
  base::sink(ofile1)

  for(j in 1:base::length(fractal_shorter)){
    delta<-base::as.numeric(names(fractal_shorter)[2])-base::as.numeric(names(fractal_shorter)[1])  #base::as.numeric(names(fractal)[2])-base::as.numeric(names(fractal)[1])
    base::cat( base::paste0("\n",nam,"\t",base::as.numeric(names(fractal_shorter)[j]),"\t",
                            base::as.numeric(names(fractal_shorter)[j])+delta,"\t",
                            fractal_shorter[j]) )
  }
  base::sink()
  grDevices::pdf(ofi)
  ofi<-(base::file.path(outpathfractal, base::paste0("FRACTAL",fname,nam,"_range_",numrange,"_sample_",sample_every,"_fft_",
                                                     seqlen,"_spar1_info.txt")))
  base::sink(ofi)

  main<- base::paste(nam,"fft window",seqlen)
  run_sum_Fractal(fracdim_list=fractal_shorter, numrange= base::round((maxf-minf)/(15*seqlen)),pflag=TRUE,main=main)
  grDevices::dev.off()
  base::sink()

  run_largeimage(All_spec_long,fname=fname, inpath=inpath, outpath=outpath,nam=nam,samplesize=samplesize,
                 r0=base::c(8000,28000),r1=base::c(1000,8000),r2=base::c(400,10000),#r3=base::c(8000,50000),r4=base::c(1000,20000),
                 inmain= base::paste("Individual", fftlengthlarge/1000,"window"))

}

#' plot_all_chromosomes
#'
#' plot Shannon and barplots for all chromosomes
#'
#' @param fname
#'
#' @return plots of all chromosomes
#'
#' @examples
#' function()
#' @export
plot_all_chromosomes <- function(fname=fname, inpath=inpath, outpath=outpath){

  summary_path <- paste0(outpath,"/", fname, "/", "Summary_output/output_data")
  file_list <- list.files(summary_path, full.names=TRUE)
  Histogram_list <- file_list[grep("Histogram", file_list)]
  Shannon_list <- file_list[grep("Shannon", file_list)]

  #-----------------------------------
  # join Shannon div data
  lsd <- lapply(Shannon_list, read.table)
  sd_chr_list0 <- basename(Shannon_list)
  sd_chr_list1 <- stringr::str_split(sd_chr_list0, "_", simplify =TRUE)
  sd_chr_list2 <- sd_chr_list1[,3]
  names(lsd) <- sd_chr_list2
  Shannon_div_total <- dplyr::bind_rows(lsd, .id = 'chromosome')
  colnames(Shannon_div_total) <- c("Chromosome", "Genome_position", "Shannon_div")

  Shannon_div_total$Chrnum <- as.numeric(stringr::str_split(Shannon_div_total$Chromosome, "r", simplify =TRUE)[,2])

  grDevices::png(file=paste0(outpath,"/", fname,"/Summary_output/",fname, "_Shannon_div.png"), width = 1000, height = 700)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=Shannon_div_total, ggplot2::aes(x=Genome_position, y=Shannon_div))+
      ggplot2::geom_point(ggplot2::aes(x=Genome_position, y=Shannon_div))+
      ggplot2::facet_wrap(~Chrnum, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/Summary_output/",fname, "_Shannon_div.pdf"))
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=Shannon_div_total, ggplot2::aes(x=Genome_position, y=Shannon_div))+
      ggplot2::geom_point(ggplot2::aes(x=Genome_position, y=Shannon_div))+
      ggplot2::facet_wrap(~Chrnum, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  # roll mean for  Shannon div data
  Shannon_div_total$Shannon_roll_mean <- zoo::rollapply(Shannon_div_total$Shannon_div, width = 100, FUN=mean, fill = NA, partial=(100/2))

  grDevices::png(file=paste0(outpath,"/", fname,"/Summary_output/",fname, "_rolling_mean_500Kbp_Shannon_div.png"), width = 1000, height = 700)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=Shannon_div_total, ggplot2::aes(x=Genome_position, y=Shannon_roll_mean))+
      ggplot2::geom_point(ggplot2::aes(x=Genome_position, y=Shannon_roll_mean))+
      ggplot2::facet_wrap(~Chrnum, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/Summary_output/",fname, "_rolling_mean_500Kbp_Shannon_div.pdf"))
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=Shannon_div_total, ggplot2::aes(x=Genome_position, y=Shannon_roll_mean))+
      ggplot2::geom_point(ggplot2::aes(x=Genome_position, y=Shannon_roll_mean))+
      ggplot2::facet_wrap(~Chrnum, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()
  #-----------------------
  # join histogram data
  ldf <- lapply(Histogram_list, read.table)
  chr_list0 <- basename(Histogram_list)
  chr_list1 <- stringr::str_split(chr_list0, "_", simplify =TRUE)
  chr_list2 <- chr_list1[,3]
  names(ldf) <- chr_list2
  Histograms_total <- dplyr::bind_rows(ldf, .id = 'chromosome')
  colnames(Histograms_total) <- c("Chromosome", "Repeat_length", "min_mean_seqval", "max_mean_seqval", "min_powsum_seqval", "max_powsum_seqval", "min_N_seqval", "max_N_seqval")
  Histograms_total <- Histograms_total[-which(Histograms_total[,8]=="max_N_seqval"),]

  Histograms_total$min_powsum_seqval <- as.numeric(Histograms_total$min_powsum_seqval)
  Histograms_total$Chrnum <- as.numeric(stringr::str_split(Histograms_total$Chromosome, "r", simplify =TRUE)[,2])

  # http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
  # https://www3.nd.edu/~steve/computing_with_data/13_Facets/facets.html
  grDevices::png(file=paste0(outpath,"/", fname,"/Summary_output/",fname, "_Histograms.png"), width = 1000, height = 700)
  print(
    ggplot2::ggplot(Histograms_total, ggplot2::aes(x=min_powsum_seqval)) +
      ggplot2::geom_histogram(stat="count", bins=10, position = "stack", colour="black", fill="white")+
      ggplot2::facet_wrap(~Chrnum, scales = "free")+
      ggplot2::theme_classic()
  )
  #ggplot2::geom_density(alpha=.2, fill="#FF6666")
  grDevices::dev.off()

  grDevices::pdf(file=paste0(outpath,"/", fname,"/Summary_output/",fname, "_Histograms.pdf"))
  print(
    ggplot2::ggplot(Histograms_total, ggplot2::aes(x=min_powsum_seqval)) +
      ggplot2::geom_histogram(stat="count", bins=10, position = "stack", colour="black", fill="white")+
      ggplot2::facet_wrap(~Chrnum, scales = "free")+
      ggplot2::theme_classic()
  )
  #ggplot2::geom_density(alpha=.2, fill="#FF6666")
  grDevices::dev.off()
}

#' DNAwalks_with_genes
#'
#' Creates 2D and 1D  DNA walks of segements of chromosomes showing the genes plotted in blue (forward), red (reverse), black (intergenic) and yellow (overlapping forward and reverse).
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
DNAwalks_with_genes <- function(chromosome=chromosome, fname=fname, inpath=inpath,
                                outpath=outpath, gff3_path=gff3_path, cytchr=cytchr,
                                start=start, end=end){

  # import the gff3 annotation file
  cytobands0 <-utils::read.table(gff3_path, sep = "\t", header=FALSE, check.names = FALSE)

  colnames(cytobands0) <- c("chrnum", "Program", "genetype", "Start", "End", "dot1", "Strand", "dot2", "ID")

  cytobands <- cytobands0[,-c(2,6,8)]

  cytobands<- as.data.frame(cytobands)
  unique(cytobands$genetype)
  unique(cytobands$chrnum)
  unique(cytobands$Strand)

  # extract only the CDS
  CDS_cytobands <- cytobands[which(cytobands$chrnum==cytchr),]
  chr_cytobands <- CDS_cytobands[which(CDS_cytobands$genetype=="mRNA"),]

  #nrow(cytobands)
  #nrow(CDS_cytobands)
  #nrow(chr_cytobands)
  # 4133

  # select only genes in the start to end range
  chr_cytobands <- chr_cytobands[which(chr_cytobands$Start>=start),]
  chr_cytobands <- chr_cytobands[which(chr_cytobands$Start<=end),]

  # how many genes?
  nrow(chr_cytobands)
  # 4000
  #--------------------------------
  # set these regions as colours
  #walk_colours <- c("#00FFFF", "#FF00FF", "black", "#CCFF00")
  walk_colours <- c("red", "blue", "black", "#CCFF00")
  cytoband_types <- c("+", "-", "not_gene", "both")
  cyto_type_colours <- as.data.frame(cbind(cytoband_types, walk_colours))

  chr_cytobands
  #chr_cytobands <- Pos_chr_cytobands
  # order by starting position
  chr_cytobands$Start <- as.numeric(chr_cytobands$Start)
  chr_cytobands$End <- as.numeric(chr_cytobands$End)
  chr_cytobands <- chr_cytobands[order(chr_cytobands$Start),]

  walk_bands <- NULL
  skip=FALSE

  # add in the starting values
  chr_cytobands$Start0 <- chr_cytobands$Start
  chr_cytobands$Start <- chr_cytobands$Start-start
  chr_cytobands$End0 <- chr_cytobands$End
  chr_cytobands$End <- chr_cytobands$End-start
  START0 = chr_cytobands$Start[1]
  reptimes <- round(START0)
  colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == "not_gene")]
  coltmp <- rep(colour, reptimes)
  walk_bands <- c(coltmp, walk_bands)

  for (i in 1:nrow(chr_cytobands)) {
    START = chr_cytobands$Start[i]
    END = chr_cytobands$End[i]
    START2 = chr_cytobands$Start[i+1]
    END2 = chr_cytobands$End[i+1]
    print(skip)
    print(i)
    if (!is.na(START2)){
      if (skip==FALSE){
        reptimes <- END-START
        print(paste0("Gene = ",reptimes))
        colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == chr_cytobands$Strand[i])]
        coltmp <- rep(colour, reptimes)
        walk_bands <- c(walk_bands, coltmp)
        print(START)
        print(END)
        print(length(walk_bands))
        skip=FALSE
        if (END < START2){  # genes have space between them
          reptimes2 <- START2-END
          print(paste0("Space = ",reptimes2))
          colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == "not_gene")]
          coltmp <- rep(colour, reptimes2)
          walk_bands <- c(walk_bands, coltmp)
          print(START)
          print(END)
          print(length(walk_bands))
          skip=FALSE
        }
        if (END == START2){ # no overlap but no space, nothing to add but don't skip
          skip=FALSE
        }
        if (END > START2){
          if (END >= END2){	  # genes overlap completely
            skip=TRUE
          }
          if (END < END2){ # genes overlap partially
            reptimes3 <- START2-END
            print(paste0("Overlap = ",reptimes3))
            colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == chr_cytobands$Strand[i])]
            coltmp <- rep(colour, abs(reptimes3))
            walk_bands <- walk_bands[-c((length(walk_bands)- abs(reptimes3)+1):length(walk_bands))]
            walk_bands <- c(walk_bands, coltmp)
            print(START)
            print(END)
            print(length(walk_bands))
            skip=FALSE
          }
        }
      } else{
        START = chr_cytobands$Start[i]
        END = chr_cytobands$End[i]
        START2 = chr_cytobands$Start[i+1]
        END2 = chr_cytobands$End[i+1]
        if (END2 > length(walk_bands)) {
          skip=FALSE
          print("Next Partial overlap")
          if (START2 >= length(walk_bands)) {
            reptimes4 <- START2-length(walk_bands)
            print(paste0("Space = ",reptimes4))
            colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == "not_gene")]
            coltmp <- rep(colour, reptimes4)
            walk_bands <- c(walk_bands, coltmp)
            print(START)
            print(END)
            print(length(walk_bands))
            skip=FALSE
          }
          if (START2 < length(walk_bands)) {
            reptimes5 <- START2-length(walk_bands)
            print(paste0("Overlap = ",reptimes5))
            colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == chr_cytobands$Strand[i])]
            coltmp <- rep(colour, abs(reptimes5))
            walk_bands <- walk_bands[-c((length(walk_bands)- abs(reptimes5)+1):length(walk_bands))]
            walk_bands <- c(walk_bands, coltmp)
            print(START)
            print(END)
            print(length(walk_bands))
            skip=FALSE
          }
        }
        else {
          skip=TRUE
          print("Another complete overlap")
          print(START)
          print(END)
          print(length(walk_bands))
        }
      }
    }else{ # START2 is NA and loop is done
      reptimes6 <- END-length(walk_bands)
      if (reptimes6 >0){
        print(paste0("Space = ",reptimes6))
        colour <- cyto_type_colours$walk_colours[which(cyto_type_colours == "not_gene")]
        coltmp <- rep(colour, reptimes6)
        walk_bands <- c(walk_bands, coltmp)
        print(START)
        print(END)
        print(length(walk_bands))
        print("Loop complete")
      }else{
        print(START)
        print(END)
        print(length(walk_bands))
        print("Loop complete")
      }

    }
  }

  #utils::write.table(walk_bands, file=paste0(outpath, "/", fname, "/", chromosome,"/", fname, "_", chromosome,"_DNAwalkcolours.txt"), append = FALSE, sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)

  walk_bands <- as.data.frame(walk_bands)

  #----------------------
  # load in DNA sequence
  dna1 <- utils::read.table(paste0(inpath, fname, "_", chromosome,".fasta"), sep = "\t", header=FALSE, check.names = FALSE)

  # remove header
  dna1 <- dna1[-1,]

  # concatenate
  waxy <- paste0(dna1, collapse="\n")
  waxy<-stringr::str_replace_all(waxy, "([\n])", "")
  waxy<-(base::unlist(base::strsplit(waxy,split="")))

  # fold upper case
  x <- base::casefold(waxy, upper=F)

  # crop to only keep start-end values
  x <- x[c(start:end)]
  # test with random generated sequence
  #dna1 <- utils::read.table(paste0(inpath,"Random_sequence.csv"), sep = ".", header=FALSE, check.names = FALSE)
  #x <- c(dna1$V1)
  #fname="TEST"

  # build the DNA walk

  ATval<- base::rep(0,base::length(base::c(x)))
  CGval<- base::rep(0,base::length(base::c(x)))
  ACval<- base::rep(0,base::length(base::c(x)))
  TGval<- base::rep(0,base::length(base::c(x)))
  gc()
  MAGval<- base::rep(0,base::length(base::c(x)))
  MATval<- base::rep(0,base::length(base::c(x)))

  Cindices<-base::grep("c",x=x)
  Gindices<-base::grep("g",x=x)
  Tindices<-base::grep("t",x=x)
  Aindices<-base::grep("a",x=x)

  gc()

  ATval[Cindices]<-  0
  ATval[Gindices]<-  0
  ATval[Aindices]<-  1;  MAGval[Aindices]<-  1
  ATval[Tindices]<- (-1);  MAGval[Tindices]<- -1

  CGval[Cindices]<-  1;  MAGval[Cindices]<-  1
  CGval[Gindices]<- (-1);  MAGval[Gindices]<- -1
  CGval[Aindices]<-  0
  CGval[Tindices]<-  0

  gc()
  ACval[Cindices]<-  (-1);  MATval[Cindices]<- -1
  ACval[Gindices]<-  0
  ACval[Aindices]<-  1;  MATval[Aindices]<-  1
  ACval[Tindices]<-  0

  TGval[Cindices]<-  0
  TGval[Gindices]<- (-1);  MATval[Gindices]<- -1
  TGval[Aindices]<-  0
  TGval[Tindices]<-  1;  MATval[Tindices]<-  1

  #Walklist<-base::list(atwalk=ATval,cgwalk=CGval,dnawalk=MAGval)
  Walklist<-base::list(atwalk=ATval,cgwalk=CGval, acwalk=ACval,tgwalk=TGval)

  # convert to cumulative walk

  atwalk<-base::cumsum(Walklist$atwalk)   #cumulative sum of AT single walk values
  cgwalk<-base::cumsum(Walklist$cgwalk)    #cumulative sum of CG single walk values
  acwalk<-base::cumsum(Walklist$acwalk)   #cumulative sum of AT single walk values
  tgwalk<-base::cumsum(Walklist$tgwalk)    #cumulative sum of CG single walk values

  #DNAwalk_long <- cbind(atwalk, cgwalk) #, acwalk, tgwalk)

  #---------------------------------------
  #dnawalk <- as.data.frame(DNAwalk_long)

  gc()
  # set region to plot
  length(atwalk)
  #chr_cytobands$Start
  walk_bands$x <- walk_bands$walk_bands
  length(walk_bands$x)

  # write 2D walk to png
  ofile <-  paste0(outpath, "/", fname, "/", chromosome,"/", "Genes_coloured_", fname, "_", chromosome,"_", start, "_", end)
  # Dryas specific
  # ofile <-  paste0(outpath, "/",  "Genes_coloured_", fname, "_", chromosome,"_", start, "_", end)

  ofile_AT_CG <- paste0(ofile, "_2D_DNAwalk_AT_CG.png")
  base::cat("\n ouput to", ofile_AT_CG)
  grDevices::png(file=ofile_AT_CG, width = 700, height = 500)
  plot(atwalk, cgwalk, cex=0.05, xlab="AT DNA Walk",  ylab="CG DNA Walk", col = walk_bands$x[c(1:length(atwalk))])
  grDevices::dev.off()
  gc()

  genomepos <- c(1:length(atwalk))+start

  # plot 1D DNAwalks
  ofile_AT <-  paste0(ofile, "_1D_DNAwalk_AT.png")
  base::cat("\n ouput to", ofile_AT)
  grDevices::png(file=ofile_AT, width = 700, height = 500)
  plot(genomepos, atwalk, cex=0.05, xlab="Genome Position",  ylab="AT DNA Walk", col =  walk_bands$x[c(1:length(atwalk))])
  grDevices::dev.off()

  ofile_CG <-  paste0(ofile, "_1D_DNAwalk_CG.png")
  base::cat("\n ouput to", ofile_CG)
  grDevices::png(file=ofile_CG, width = 700, height = 500)
  plot(genomepos, cgwalk, cex=0.05, xlab="Genome Position",  ylab="CG DNA Walk", col =  walk_bands$x[c(1:length(atwalk))])
  grDevices::dev.off()

  #------------------
  # plot the smoothed spline for AT
  spar=0.6
  lowpass1.spline <- stats::smooth.spline(genomepos, atwalk, spar = spar) ## Control spar for amount of smoothing
  # plot 1D DNAwalks
  ofile_AT <-  paste0(ofile, "_1D_DNAwalk_AT_with_spline.png")
  base::cat("\n ouput to", ofile_AT)
  grDevices::png(file=ofile_AT, width = 700, height = 500)
  plot(genomepos, atwalk,  cex=0.05, xlab="Genome Position",  ylab="AT DNA Walk", col =  walk_bands$x[c(1:length(atwalk))])
  graphics::lines(stats::predict(lowpass1.spline, genomepos), col = "red", lwd = 1.5)
  grDevices::dev.off()

  ofile_AT <-  paste0(ofile, "_1D_DNAwalk_AT_smoothed.png")
  base::cat("\n ouput to", ofile_AT)
  grDevices::png(file=ofile_AT, width = 700, height = 500)
  highpass1 <- atwalk - stats::predict(lowpass1.spline, genomepos)$y
  #graphics::lines(genomepos, highpass1, lwd =  2)
  base::plot(genomepos, highpass1, type="l", pch=15, lwd =  3, xlab="Genome Position", ylab="Smoothed AT DNA Walk")
  grDevices::dev.off()

  #------------------
  # other

  # ofile_AC_TG <- paste0(ofile, "_2D_DNAwalk_AC_TG.png")
  # base::cat("\n ouput to", ofile_AC_TG)
  # grDevices::png(file=ofile_AC_TG, width = 3000, height = 3000)
  # plot(acwalk, tgwalk, col = walk_bands$x[c(1:length(atwalk))]) #"black") # walk_bands$x[c(start:end)])
  # grDevices::dev.off()
  # gc()
  #
  # ofile_AC <-  paste0(ofile, "_1D_DNAwalk_AC.png")
  # base::cat("\n ouput to", ofile_AC)
  # grDevices::png(file=ofile_AC, width = 3000, height = 1000)
  # plot(genomepos, acwalk, col = walk_bands$x[c(1:length(atwalk))]) #"black") # walk_bands$x[c(start:end)])
  # grDevices::dev.off()
  #
  # ofile_TG <-  paste0(ofile, "_1D_DNAwalk_TG.png")
  # base::cat("\n ouput to", ofile_TG)
  # grDevices::png(file=ofile_TG, width = 3000, height = 1000)
  # plot(genomepos, tgwalk, col = walk_bands$x[c(1:length(atwalk))]) #"black") # walk_bands$x[c(start:end)])
  # grDevices::dev.off()

  # should also try
  # AG and TC

}

#' gene_slopes
#'
#' Plots graphs of the slopes of each DNA walk for each gene
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
gene_slopes <- function(chromosome=chromosome, fname=fname, inpath=inpath,
                                outpath=outpath, gff3_path=gff3_path, cytchr=cytchr,
                                start=start, end=end){


  # import the gff3 annotation file
  cytobands0 <-utils::read.table(gff3_path, sep = "\t", header=FALSE, check.names = FALSE)

  colnames(cytobands0) <- c("chrnum", "Program", "genetype", "Start", "End", "dot1", "Strand", "dot2", "ID")

  cytobands <- cytobands0[,-c(2,6,8)]

  cytobands<- as.data.frame(cytobands)
  unique(cytobands$genetype)
  unique(cytobands$chrnum)
  unique(cytobands$Strand)

  # extract only the CDS
  CDS_cytobands <- cytobands[which(cytobands$chrnum==cytchr),]
  chr_cytobands <- CDS_cytobands[which(CDS_cytobands$genetype=="mRNA"),]

  #nrow(cytobands)
  #nrow(CDS_cytobands)
  #nrow(chr_cytobands)
  # 4133

  # select only genes in the start to end range
  chr_cytobands <- chr_cytobands[which(chr_cytobands$Start>=start),]
  chr_cytobands <- chr_cytobands[which(chr_cytobands$Start<=end),]

  # how many genes?
  nrow(chr_cytobands)
  # 4000
  #--------------------------------
  # set these regions as colours
  #walk_colours <- c("#00FFFF", "#FF00FF", "black", "#CCFF00")
  walk_colours <- c("red", "blue", "black", "#CCFF00")
  cytoband_types <- c("+", "-", "not_gene", "both")
  cyto_type_colours <- as.data.frame(cbind(cytoband_types, walk_colours))

  chr_cytobands
  #chr_cytobands <- Pos_chr_cytobands
  # order by starting position
  chr_cytobands$Start <- as.numeric(chr_cytobands$Start)
  chr_cytobands$End <- as.numeric(chr_cytobands$End)
  chr_cytobands <- chr_cytobands[order(chr_cytobands$Start),]

  walk_bands <- NULL
  skip=FALSE

  #--------------------------------------------------
  # load in DNA sequence
  dna1 <- utils::read.table(paste0(inpath, fname, "_", chromosome,".fasta"), sep = "\t", header=FALSE, check.names = FALSE)

  # remove header
  dna1 <- dna1[-1,]

  # concatenate
  waxy <- paste0(dna1, collapse="\n")
  waxy<-stringr::str_replace_all(waxy, "([\n])", "")
  waxy<-(base::unlist(base::strsplit(waxy,split="")))

  # fold upper case
  x <- base::casefold(waxy, upper=F)

  # crop to only keep start-end values
  x <- x[c(start:end)]
  # test with random generated sequence
  #dna1 <- utils::read.table(paste0(inpath,"Random_sequence.csv"), sep = ".", header=FALSE, check.names = FALSE)
  #x <- c(dna1$V1)
  #fname="TEST"

  # build the DNA walk

  ATval<- base::rep(0,base::length(base::c(x)))
  CGval<- base::rep(0,base::length(base::c(x)))
  ACval<- base::rep(0,base::length(base::c(x)))
  TGval<- base::rep(0,base::length(base::c(x)))
  AGval<- base::rep(0,base::length(base::c(x)))
  TCval<- base::rep(0,base::length(base::c(x)))

  gc()
  MAGval<- base::rep(0,base::length(base::c(x)))
  MATval<- base::rep(0,base::length(base::c(x)))
  MACval<- base::rep(0,base::length(base::c(x)))

  Cindices<-base::grep("c",x=x)
  Gindices<-base::grep("g",x=x)
  Tindices<-base::grep("t",x=x)
  Aindices<-base::grep("a",x=x)

  gc()

  ATval[Cindices]<-  0
  ATval[Gindices]<-  0
  ATval[Aindices]<-  1;  MAGval[Aindices]<-  1
  ATval[Tindices]<- (-1);  MAGval[Tindices]<- -1

  CGval[Cindices]<-  1;  MAGval[Cindices]<-  1
  CGval[Gindices]<- (-1);  MAGval[Gindices]<- -1
  CGval[Aindices]<-  0
  CGval[Tindices]<-  0

  gc()
  ACval[Cindices]<-  (-1);  MATval[Cindices]<- -1
  ACval[Gindices]<-  0
  ACval[Aindices]<-  1;  MATval[Aindices]<-  1
  ACval[Tindices]<-  0

  TGval[Cindices]<-  0
  TGval[Gindices]<- (-1);  MATval[Gindices]<- -1
  TGval[Aindices]<-  0
  TGval[Tindices]<-  1;  MATval[Tindices]<-  1

  gc()
  AGval[Cindices]<-  0
  AGval[Gindices]<-  (-1); MACval[Gindices]<- -1
  AGval[Aindices]<-  1;  MACval[Aindices]<-  1
  AGval[Tindices]<-  0

  TCval[Cindices]<-  (-1);  MACval[Cindices]<- -1
  TCval[Gindices]<-  0
  TCval[Aindices]<-  0
  TCval[Tindices]<-  1;  MACval[Tindices]<-  1

  #Walklist<-base::list(atwalk=ATval,cgwalk=CGval,dnawalk=MAGval)
  Walklist<-base::list(atwalk=ATval,cgwalk=CGval, acwalk=ACval,tgwalk=TGval, agwalk=AGval,tcwalk=TCval)

  # convert to cumulative walk

  atwalk<-base::cumsum(Walklist$atwalk)   #cumulative sum of AT single walk values
  cgwalk<-base::cumsum(Walklist$cgwalk)    #cumulative sum of CG single walk values

  acwalk<-base::cumsum(Walklist$acwalk)   #cumulative sum of AC single walk values
  tgwalk<-base::cumsum(Walklist$tgwalk)    #cumulative sum of TG single walk values

  agwalk<-base::cumsum(Walklist$agwalk)   #cumulative sum of AG single walk values
  tcwalk<-base::cumsum(Walklist$tcwalk)    #cumulative sum of TC single walk values

  #DNAwalk_long <- cbind(atwalk, cgwalk) #, acwalk, tgwalk)

  #------------------------------------------
  # calculate slopes at and cg

  genomepos <- c(1:length(atwalk))+start
  DNAwalk_long<- as.data.frame(cbind(genomepos, atwalk, cgwalk))

  Genes_long <- as.data.frame(chr_cytobands[,c(1,3,4,5)])

  colnames(DNAwalk_long) <- c("Start", "atwalk", "cgwalk")
  Gene_starts <- dplyr::left_join(Genes_long, DNAwalk_long, by="Start")

  colnames(DNAwalk_long) <- c("End", "atwalk", "cgwalk")
  Gene_ends <- dplyr::left_join(Genes_long, DNAwalk_long, by="End")

  colnames(Gene_starts) <- c("chrnum", "Start", "End", "Strand", "Start_atpos", "Start_cgpos")
  colnames(Gene_ends) <- c("chrnum", "Start", "End", "Strand", "End_atpos", "End_cgpos")

  Gene_walks <- as.data.frame(merge(Gene_starts, Gene_ends))

  Gene_walks <- dplyr::distinct(Gene_walks)

  Gene_walks$Length <- (Gene_walks$End-Gene_walks$Start)


  Gene_walks$AT_slope <- (Gene_walks$End_atpos - Gene_walks$Start_atpos)/(Gene_walks$End-Gene_walks$Start)
  Gene_walks$CG_slope <- (Gene_walks$End_cgpos - Gene_walks$Start_cgpos)/(Gene_walks$End-Gene_walks$Start)
  Gene_walks$AT_CG_slope <- (Gene_walks$End_atpos - Gene_walks$Start_atpos)/(Gene_walks$End_cgpos-Gene_walks$Start_cgpos)

  Slope_by_strand <- dplyr::group_by(Gene_walks, by=Strand)
  Slopes_per_strand_AT <- dplyr::summarise(Slope_by_strand, ATmax = max(AT_slope, na.rm = T),
                                           ATmin = min(AT_slope, na.rm = T),
                                           ATmean = mean(AT_slope , na.rm = T))

  Slopes_per_strand_CG <- dplyr::summarise(Slope_by_strand, CGmax = max(CG_slope, na.rm = T),
                                           CGmin = min(CG_slope, na.rm = T),
                                           CGmean = mean(CG_slope , na.rm = T))
  #------------------------------------
  # think about slope of best fit line through walk - might fix the short genes that get affected by a few values
  colnames(DNAwalk_long) <- c("genomepos", "atwalk", "cgwalk")
  Gene_walks$best_fit_AT_slope <- Gene_walks$AT_slope
  Gene_walks$best_fit_CG_slope <- Gene_walks$AT_slope
  for (i in 1:nrow(Gene_walks)){
    if (Gene_walks$Start[i]<Gene_walks$End[i]){
    DNAwalk_subset <- DNAwalk_long[which(DNAwalk_long$genomepos >= Gene_walks$Start[i] & DNAwalk_long$genomepos <=Gene_walks$End[i]),]
    Gene_walks$best_fit_AT_slope[i] <- coef(lm(DNAwalk_subset$atwalk~DNAwalk_subset$genomepos))[2]
    Gene_walks$best_fit_CG_slope[i] <- coef(lm(DNAwalk_subset$cgwalk~DNAwalk_subset$genomepos))[2]
    } else {
      DNAwalk_subset <- DNAwalk_long[which(DNAwalk_long$genomepos <= Gene_walks$Start[i] & DNAwalk_long$genomepos >=Gene_walks$End[i]),]
      Gene_walks$best_fit_AT_slope[i] <- coef(lm(DNAwalk_subset$atwalk~DNAwalk_subset$genomepos))[2]
      Gene_walks$best_fit_CG_slope[i] <- coef(lm(DNAwalk_subset$cgwalk~DNAwalk_subset$genomepos))[2]
    }
  }

  #-----------------------------------
  # plots
  png(paste0(outpath, "/",fname, "_gene_fitslopes_at.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_AT_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="AT Slope")
  stripchart(best_fit_AT_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_fitslopes_cg.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_CG_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="CG Slope")
  stripchart(best_fit_CG_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_fitslopes_at_cg.png"), width = 1000, height = 1000)
  plot(best_fit_CG_slope ~ best_fit_AT_slope, data=Gene_walks, pch=16)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_fitslopes_AT_length.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_AT_slope ~ Length, data=Gene_walks, las=2, xlab="Gene Length", ylab="AT Slope")
  stripchart(best_fit_AT_slope ~ Length, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_fitslopes_AT_length_strand.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_AT_slope ~ Length+Strand, data=Gene_walks, las=2, xlab="Length/Strand", ylab="AT Slope")
  stripchart(best_fit_AT_slope ~ Length+Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()


  png(paste0(outpath, "/",fname, "test.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(best_fit_AT_slope ~ AT_slope, data=Gene_walks, las=2, xlab="Strand", ylab="AT Slope")
  stripchart(best_fit_AT_slope ~ AT_slope, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()
  #---------
  png(paste0(outpath, "/",fname, "_Gene_start_pos.png"), width = 300, height = 15000)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  #boxplot(Start ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="Genome Position")
  stripchart(Start ~ Strand, data=Gene_walks, vertical = TRUE, pch = 16,las = 1, col = "black", cex = 1.5)#, add=TRUE)
  dev.off()

  #-------------------------------------------
  # plots
  png(paste0(outpath, "/",fname, "_gene_slopes_at.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(AT_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="AT Slope")
  stripchart(AT_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  # this was wrong
  png(paste0(outpath, "/",fname, "_gene_slopes_cg.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(CG_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="CG Slope")
  stripchart(CG_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_slopes_at_cg.png"), width = 1000, height = 1000)
  plot(CG_slope ~ AT_slope, data=Gene_walks, pch=16)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_slopes_AT_length.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(AT_slope ~ Length, data=Gene_walks, las=2, xlab="Gene Length", ylab="AT Slope")
  stripchart(AT_slope ~ Length, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()

  png(paste0(outpath, "/",fname, "_gene_slopes_AT_length_strand.png"), width = 1000, height = 600)
  # sets the bottom, left, top and right margins
  par(mar=c(15,4.1,4.1,2.1))
  boxplot(AT_slope ~ Length+Strand, data=Gene_walks, las=2, xlab="Gene Length/Strand", ylab="AT Slope")
  stripchart(AT_slope ~ Length+Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  dev.off()


  # #------------------------------------------
  # # calculate slopes ag and tc
  #
  # genomepos <- c(1:length(agwalk))+start
  # DNAwalk_long<- as.data.frame(cbind(genomepos, agwalk, tcwalk))
  #
  # Genes_long <- as.data.frame(chr_cytobands[,c(1,3,4,5)])
  #
  # colnames(DNAwalk_long) <- c("Start", "agwalk", "tcwalk")
  # Gene_starts <- dplyr::left_join(Genes_long, DNAwalk_long, by="Start")
  #
  # colnames(DNAwalk_long) <- c("End", "agwalk", "tcwalk")
  # Gene_ends <- dplyr::left_join(Genes_long, DNAwalk_long, by="End")
  #
  # colnames(Gene_starts) <- c("chrnum", "Start", "End", "Strand", "Start_agpos", "Start_tcpos")
  # colnames(Gene_ends) <- c("chrnum", "Start", "End", "Strand", "End_agpos", "End_tcpos")
  #
  # Gene_walks <- as.data.frame(merge(Gene_starts, Gene_ends))
  #
  # Gene_walks <- dplyr::distinct(Gene_walks)
  #
  #
  # Gene_walks$ag_slope <- (Gene_walks$End_agpos - Gene_walks$Start_agpos)/(Gene_walks$End-Gene_walks$Start)
  # Gene_walks$tc_slope <- (Gene_walks$End_tcpos - Gene_walks$Start_tcpos)/(Gene_walks$End-Gene_walks$Start)
  # Gene_walks$ag_tc_slope <- (Gene_walks$End_agpos - Gene_walks$Start_agpos)/(Gene_walks$End_tcpos-Gene_walks$Start_tcpos)
  #
  # Gene_walks$Length <- (Gene_walks$End-Gene_walks$Start)
  #
  # Slope_by_strand <- dplyr::group_by(Gene_walks, by=Strand)
  # Slopes_per_strand_ag <- dplyr::summarise(Slope_by_strand, agmax = max(ag_slope, na.rm = T),
  #                                          agmin = min(ag_slope, na.rm = T),
  #                                          agmean = mean(ag_slope , na.rm = T))
  #
  # Slopes_per_strand_tc <- dplyr::summarise(Slope_by_strand, tcmax = max(tc_slope, na.rm = T),
  #                                          tcmin = min(tc_slope, na.rm = T),
  #                                          tcmean = mean(tc_slope , na.rm = T))
  #
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ag.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(ag_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="ag Slope")
  # stripchart(ag_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_tc.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(tc_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="tc Slope")
  # stripchart(ag_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ag_tc.png"), width = 1000, height = 1000)
  # plot(tc_slope ~ ag_slope, data=Gene_walks, pch=16)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ag_length.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(ag_slope ~ Length, data=Gene_walks, las=2, xlab="Strand", ylab="ag Slope")
  # stripchart(ag_slope ~ Length, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # #----------------------
  # LinearAG <- aov(ag_slope ~ Strand, data = Gene_walks)
  # summary(LinearAG)
  #
  # LinearTC <- aov(tc_slope ~ Strand, data = Gene_walks)
  # summary(LinearTC)
  #
  # #------------------------------------------
  # # calculate slopes ac and tg
  #
  # genomepos <- c(1:length(acwalk))+start
  # DNAwalk_long<- as.data.frame(cbind(genomepos, acwalk, tgwalk))
  #
  # Genes_long <- as.data.frame(chr_cytobands[,c(1,3,4,5)])
  #
  # colnames(DNAwalk_long) <- c("Start", "acwalk", "tgwalk")
  # Gene_starts <- dplyr::left_join(Genes_long, DNAwalk_long, by="Start")
  #
  # colnames(DNAwalk_long) <- c("End", "acwalk", "tgwalk")
  # Gene_ends <- dplyr::left_join(Genes_long, DNAwalk_long, by="End")
  #
  # colnames(Gene_starts) <- c("chrnum", "Start", "End", "Strand", "Start_acpos", "Start_tgpos")
  # colnames(Gene_ends) <- c("chrnum", "Start", "End", "Strand", "End_acpos", "End_tgpos")
  #
  # Gene_walks <- as.data.frame(merge(Gene_starts, Gene_ends))
  #
  # Gene_walks <- dplyr::distinct(Gene_walks)
  #
  #
  # Gene_walks$ac_slope <- (Gene_walks$End_acpos - Gene_walks$Start_acpos)/(Gene_walks$End-Gene_walks$Start)
  # Gene_walks$tg_slope <- (Gene_walks$End_tgpos - Gene_walks$Start_tgpos)/(Gene_walks$End-Gene_walks$Start)
  # Gene_walks$ac_tg_slope <- (Gene_walks$End_acpos - Gene_walks$Start_acpos)/(Gene_walks$End_tgpos-Gene_walks$Start_tgpos)
  #
  # Slope_by_strand <- dplyr::group_by(Gene_walks, by=Strand)
  # Slopes_per_strand_ac <- dplyr::summarise(Slope_by_strand, acmax = max(ac_slope, na.rm = T),
  #                                          acmin = min(ac_slope, na.rm = T),
  #                                          acmean = mean(ac_slope , na.rm = T))
  #
  # Slopes_per_strand_tg <- dplyr::summarise(Slope_by_strand, tgmax = max(tg_slope, na.rm = T),
  #                                          tgmin = min(tg_slope, na.rm = T),
  #                                          tgmean = mean(tg_slope , na.rm = T))
  #
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ac.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(ac_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="ac Slope")
  # stripchart(ac_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_tg.png"), width = 1000, height = 600)
  # # sets the bottom, left, top and right margins
  # par(mar=c(15,4.1,4.1,2.1))
  # boxplot(tg_slope ~ Strand, data=Gene_walks, las=2, xlab="Strand", ylab="tg Slope")
  # stripchart(ac_slope ~ Strand, data=Gene_walks, vertical = TRUE, method = "jitter", pch = 16,las = 1, col = "black", cex = 1.5, add=TRUE)
  # dev.off()
  #
  # png(paste0(outpath, "/",fname, "_gene_slopes_ac_tg.png"), width = 1000, height = 1000)
  # plot(tg_slope ~ ac_slope, data=Gene_walks, pch=16)
  # dev.off()
  #
}


#' PCA_repeats
#'
#' Plots a PCA with repeat lengths as variables and genome position as individuals, then forms kmeans groups and colours Shannon div plots based on groups
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
PCA_repeats <- function(haplotype_path=haplotype_path, Spp_chromosome=Spp_chromosome){

  # requires the folowing libraries
  library(RepeatOBserverV1)
  library(corrr)
  library(ggcorrplot)
  library(FactoMineR)
  library(factoextra)
  library(ggbiplot)
  library(ggfortify)
  library(ggplot2)
  library(dplyr)
  library(factoextra)
  library(cluster)
  library(tidyverse)


  #install.packages("corrr")
  #install.packages("ggcorrplot")
  #install.packages("FactoMineR")
  #install.packages("factoextra")
  #install.packages("tidyr")
  #install.packages('ggfortify')
  #install.packages("devtools")
  #library(devtools)
  #install_github("vqv/ggbiplot")
  #install.packages("tidyverse")

  file_list <- list.files(haplotype_path, full.names=TRUE)
  Shannon_list <- file_list[grep("Shannon", file_list)]
  Fourier_list <- file_list[grep("Total", file_list)]

  #-----------------------------------
  # join Fourier data
  lf <- lapply(Fourier_list, read.table)
  lft <- lapply(lf, t)
  lft <- lapply(lft, as.data.frame)
  f_chr_list0 <- basename(Fourier_list)
  f_chr_list1 <- stringr::str_split(f_chr_list0, "_", simplify =TRUE)
  f_chr_list2 <- f_chr_list1[,2]
  names(lft) <- f_chr_list2
  Fourier_total <- dplyr::bind_rows(lft, .id = 'Haplotype')

  head(Fourier_total)

  ################################
  # try to transpose and run on genome position not repeat length
  # to see what positions in the genome cluster

  # add Haplotype to rowname and remove Haplotype column
  rownames(Fourier_total) <- paste0((stringr::str_split(rownames(Fourier_total), '501', simplify =TRUE)[,1]),"-", c(1:nrow(Fourier_total)),"_", Fourier_total[,1])
  Fourier_total <- Fourier_total[,c(2:ncol(Fourier_total))]

  # remove NA values for where they don't line up
  Fourier_total_rmna <- na.omit(Fourier_total)

  # remove all zero columns
  Fourier_total_rm0na<- Fourier_total_rmna[which(rowSums(Fourier_total_rmna) != 0),]

  # data normalization
  data_normalized <- scale(Fourier_total_rm0na)
  head(data_normalized)

  #################################
  # pca option - this one works as expected
  # https://github.com/sinhrks/ggfortify

  Fourier.pca2 <- prcomp(data_normalized, center = TRUE, scale. = TRUE)
  summary(Fourier.pca2)

  Fourier_total_rm0na$Haplotype <- stringr::str_split(rownames(Fourier_total_rm0na), "_", simplify =TRUE)[,2]

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_PCA2.png"), width = 2000, height = 2000)
  print(pca.plot <- autoplot(Fourier.pca2,
                             data = Fourier_total_rm0na,
                             colour = 'Haplotype',
                             loadings = TRUE, loadings.colour = 'black',
                             loadings.label = TRUE, loadings.label.size = 6)+
          theme_classic())
  grDevices::dev.off()

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_PCA2_labels.png"), width = 2000, height = 2000)
  print(pca.plot <- autoplot(Fourier.pca2,
                             data = Fourier_total_rm0na,
                             colour = 'Haplotype', shape = FALSE,
                             label = TRUE, label.size = 3,
                             loadings = TRUE, loadings.colour = 'black',
                             loadings.label = TRUE, loadings.label.size = 6)+
          theme_classic())
  grDevices::dev.off()

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_PCA_biplot.png"), width = 2000, height = 2000)
  print( biplot(Fourier.pca2), cex=c(0.01, 0.01))
  grDevices::dev.off()


  #######################################
  # Clustering in a PCA
  # https://medium.com/@zullinira23/implementation-of-principal-component-analysis-pca-on-k-means-clustering-in-r-794f03ec15f

  pca2_transform = as.data.frame(-Fourier.pca2$x[,1:2])

  # grDevices::png(file=paste0(haplotype_path, "/Arab_Chr4_PCA_cluster_kchoice.png"), width = 700, height = 700)
  # print(fviz_nbclust(pca2_transform, kmeans, method = 'wss'))
  # grDevices::dev.off()

  kmeans_calc <- function(pca2_transform, k){
    kmeans_pca = kmeans(pca2_transform, centers = k, nstart = 50)

    grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_PCA_clustering_k", k,".png"), width = 1000, height = 1000)
    print(fviz_cluster(kmeans_pca, data = pca2_transform))
    grDevices::dev.off()

    # split and add postion names
    Genome_pos_Haplo <- stringr::str_split(names(kmeans_pca$cluster), "_", simplify =TRUE)
    Genome_pos <- stringr::str_split(Genome_pos_Haplo[,1], "-", simplify =TRUE)
    XX <- cbind(Genome_pos_Haplo, Genome_pos, kmeans_pca$cluster)
    colnames(XX) <- c("Genome_pos", "Haplotype", "Postion", "ID", "Kgroup")

    # make a data frame
    XX.df <- as.data.frame(XX)
    head(XX.df)

    # remove the rows that contain X or X1
    XX.df <- XX.df[-which(XX.df$Postion=="X1"),]
    XX.df <- XX.df[-which(XX.df$Postion=="X"),]

    # keep only the relevant columns and make wide
    XX.df <- XX.df[, c(2,3,5)]
    XXX <- tidyr::pivot_wider(XX.df, names_from = 'Haplotype', values_from = 'Kgroup')

    # remove X from the position column
    XX.df$Postion <- as.numeric(gsub("X","",as.character(XX.df$Postion)))*1000+501
    XX.df$Kgroup <- as.numeric(XX.df$Kgroup)

    # plot data
    grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_PCA_Kmeans_k", k, "_along_genome.png"), width = 1500, height = 700)
    print(
      ggplot2::ggplot(data=XX.df, ggplot2::aes(x=Postion, y=Kgroup))+
        ggplot2::geom_point(ggplot2::aes(x=Postion, y=Kgroup))+
        ggplot2::facet_wrap(~Haplotype, scales = "free")+
        ggplot2::theme_classic()
    )
    grDevices::dev.off()

    return(XX.df)
  }

  kmeans_3 <- kmeans_calc(pca2_transform, 3)
  kmeans_6 <- kmeans_calc(pca2_transform, 6)
  kmeans_10 <- kmeans_calc(pca2_transform, 10)
  kmeans_15 <- kmeans_calc(pca2_transform, 15)
  kmeans_25 <- kmeans_calc(pca2_transform, 25)

  ##################################
  #---------------------------
  # join kmeans values

  alldata <- kmeans_3 %>%
    left_join(kmeans_6, by=c('Haplotype', 'Postion')) %>%
    left_join(kmeans_10, by=c('Haplotype', 'Postion')) %>%
    left_join(kmeans_15, by=c('Haplotype', 'Postion')) %>%
    left_join(kmeans_25, by=c('Haplotype', 'Postion'))

  colnames(alldata) <-  c("Haplotype",  "Postion", "Kgroup.3", "Kgroup.6", "Kgroup.10", "Kgroup.15", "Kgroup.25")

  # write out the data
  utils::write.table(alldata, file=paste0(haplotype_path, "/", Spp_chromosome,"_Kmeans_info.txt"), append = FALSE, sep = ",",
                     dec = ".", row.names = FALSE, col.names = TRUE)


  #############################
  # plot Shannon diversity, coloured by PCA group, with line for Shannon cent and barplot cent

  # load Shannon data
  lsd <- lapply(Shannon_list, read.table)
  sd_chr_list0 <- basename(Shannon_list)
  sd_chr_list1 <- stringr::str_split(sd_chr_list0, "_", simplify =TRUE)
  sd_chr_list2 <- sd_chr_list1[,1]
  names(lsd) <- sd_chr_list2
  Shannon_div_total <- dplyr::bind_rows(lsd, .id = 'Haplotype')
  colnames(Shannon_div_total) <- c("Haplotype", "Postion", "Shannon_div")

  #----------------------------
  # join kmeans with Shannon data

  alldata_Shannon <- alldata %>%
    left_join(Shannon_div_total, by=c('Haplotype', 'Postion'))


  #------------------------------
  # plot the data

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_Kmeans6_colours_Shannon.png"), width = 1000, height = 700)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=alldata_Shannon, ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::geom_point(ggplot2::aes(x=Postion, y=Shannon_div), col=alldata_Shannon$Kgroup.6)+
      ggplot2::facet_wrap(~Haplotype, ncol=1, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  #------------------------------
  # plot the data

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_Kmeans10_colours_Shannon.png"), width = 1000, height = 700)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=alldata_Shannon, ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::geom_point(ggplot2::aes(x=Postion, y=Shannon_div), col=alldata_Shannon$Kgroup.10)+
      ggplot2::facet_wrap(~Haplotype, ncol=1, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  #------------------------------
  # plot the data

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_Kmeans15_colours_Shannon.png"), width = 1000, height = 700)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=alldata_Shannon, ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::geom_point(ggplot2::aes(x=Postion, y=Shannon_div), col=alldata_Shannon$Kgroup.15)+
      ggplot2::facet_wrap(~Haplotype, ncol=1, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  #--------------------------
  # plot unique plot for each group

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_K6meansplit_Shannon.png"), width = 2000, height = 1000)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=alldata_Shannon, ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::geom_point(ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::facet_wrap(~Haplotype+Kgroup.6, ncol=6)+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  #--------------------------
  # plot unique plot for each group

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_K10meansplit_Shannon.png"), width = 2000, height = 1000)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=alldata_Shannon, ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::geom_point(ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::facet_wrap(~Haplotype+Kgroup.10, ncol=10)+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  #--------------------------
  # plot unique plot for each group

  grDevices::png(file=paste0(haplotype_path, "/", Spp_chromosome,"_K15meansplit_Shannon.png"), width = 3000, height = 1000)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=alldata_Shannon, ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::geom_point(ggplot2::aes(x=Postion, y=Shannon_div))+
      ggplot2::facet_wrap(~Haplotype+Kgroup.15, ncol=15)+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

}


#' roll_sum_histogram
#'
#' Calculates and plots the roll sum abundance from the summary directory files. This code can only be run once the rest of the progam is complete.
#'
#' @param nam input dataset
#'
#' @return output dataset
#'
#' @examples
#' function()
#' @export
roll_sum_histogram <- function(fname=fname, outpath=outpath){

  # get chromosome count
  summary_path <- paste0(outpath,"/", fname, "/", "Summary_output/output_data")
  file_list <- list.files(summary_path, full.names=TRUE)
  Histogram_list <- file_list[grep("Histogram", file_list)]
  Shannon_list <- file_list[grep("Shannon", file_list)]

  #-------------------------------
  # plot for 4Mbp overlapping windows
  for (i in 1:length(Histogram_list)) {
    chromosome=i
    print(fname)
    print(chromosome)
    # get Fourier data
    All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname, "/Summary_output/output_data/Total_", fname, "_Chr", chromosome, "_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

    # sum columns
    Fourier_sums<- colSums(All_spec0)
    # define window
    wind_size=500
    # run rolling sum
    roll_sum_Fourier_sums <-  zoo::rollsum(Fourier_sums, wind_size, align = "center", fill = NA)

    # join sums with genome positions
    Genome_position <- c(1:length(roll_sum_Fourier_sums))*5000
    Repeat_abund <- cbind(Genome_position, roll_sum_Fourier_sums)

    #plot
    grDevices::png(filename = paste0(outpath,"/", fname, "/Summary_output/histograms/", fname,"_Chr", chromosome, "_Repeat_abundance_sum.png"), width = 1400, height = 700)
    plot(Genome_position, roll_sum_Fourier_sums)
    grDevices::dev.off()

    # write out total file
    utils::write.table(x=Repeat_abund, file=paste0(outpath,"/", fname, "/Summary_output/output_data/", fname,"_Chr", chromosome, "_Repeat_abundance_sum.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)

  }

  #----------------------------------
  # # make non-overlapping 500kbp windows
  # for (i in 1:length(Histogram_list)){
  # chromosome=i
  # All_spec0<-base::as.matrix(utils::read.table(paste0(outpath,"/", fname, "/Summary_output/output_data/Total_", fname, "_Chr", chromosome, "_All_spec_merged.txt"), header = TRUE, check.names = FALSE))

  # Fourier_sums <- NULL
  # for (j in 0:(ncol(All_spec0)/100)-1){
  # x=(j*100)
  # y=x+99
  # Fourier_sums[j] <- sum(colSums(All_spec0[,c(x:y)]))
  # }

  # wind_size=4e6
  # roll_sum_Fourier_sums <-  zoo::rollsum(Fourier_sums, wind_size, align = "center", fill = NA)
  # Genome_position <- c(1:length(roll_sum_Fourier_sums))*500000
  # Repeat_abund <- cbind(Genome_position, roll_sum_Fourier_sums)
  # grDevices::png(filename = paste0(outpath,"/", fname, "/Summary_output/output_data/", fname,"_Chr", chromosome, "_Repeat_abundance_sum_500kbp_", wind_size, ".png"), width = 1400, height = 700)
  # plot(Genome_position, roll_sum_Fourier_sums)
  # grDevices::dev.off()
  # utils::write.table(x=Repeat_abund, file=paste0(outpath,"/", fname, "/Summary_output/output_data/", fname,"_Chr", chromosome, "_Repeat_abundance_sum_500kbp_", wind_size, ".txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = TRUE)
  # }

  #-----------------------
  # Read in data for all chromosomes and plot
  summary_path <- paste0(outpath,"/", fname, "/", "Summary_output/output_data")
  file_list <- list.files(summary_path, full.names=TRUE)
  rollsumhist_list <- file_list[grep("Repeat_abundance_sum", file_list)]

  lsd <- lapply(rollsumhist_list, read.table)
  sd_chr_list0 <- basename(rollsumhist_list)
  sd_chr_list1 <- stringr::str_split(sd_chr_list0, "_", simplify =TRUE)
  sd_chr_list2 <- sd_chr_list1[,3]
  names(lsd) <- sd_chr_list2
  RepeatAbundance_total <- dplyr::bind_rows(lsd, .id = 'chromosome')
  colnames(RepeatAbundance_total) <- c("Chromosome", "Genome_position", "RepeatAbundance")

  RepeatAbundance_total$Chrnum <- as.numeric(stringr::str_split(RepeatAbundance_total$Chromosome, "r", simplify =TRUE)[,2])

  RepeatAbundance_total$Genome_position <- as.numeric(as.character(RepeatAbundance_total$Genome_position))
  RepeatAbundance_total$RepeatAbundance <- as.numeric(as.character(RepeatAbundance_total$RepeatAbundance))

  print("plotting all chromosomes")
  grDevices::png(file=paste0(outpath,"/", fname,"/Summary_output/",fname, "_Repeat_Sum_Abundance.png"), width = 1000, height = 700)
  # plot Shannon on one plot
  # https://www.geeksforgeeks.org/add-vertical-and-horizontal-lines-to-ggplot2-plot-in-r/
  print(
    ggplot2::ggplot(data=RepeatAbundance_total, ggplot2::aes(x=Genome_position, y=RepeatAbundance))+
      ggplot2::geom_point(ggplot2::aes(x=Genome_position, y=RepeatAbundance))+
      ggplot2::facet_wrap(~Chrnum, scales = "free")+
      ggplot2::theme_classic()
  )
  grDevices::dev.off()

  # find start and end of highly repeating regions based 1SD from min
  RepeatAbund_cent <- NULL
  RepeatAbund_min <- NULL
  RepeatAbund_length <- NULL
  for (chromosome in 1:length(rollsumhist_list)){
    RepeatAbundance_chr <- RepeatAbundance_total[which(RepeatAbundance_total$Chrnum == chromosome),]

    # find min position
    cent_min <- RepeatAbundance_chr$Genome_position[which(RepeatAbundance_chr$RepeatAbundance == min(RepeatAbundance_chr$RepeatAbundance, na.rm=TRUE))]
    #15 885 000

    # find SD of data
    SD_repeatAbund <- sd(RepeatAbundance_chr$RepeatAbundance, na.rm=TRUE)
    #29 911 051

    thres_upper = mean(RepeatAbundance_chr$RepeatAbundance, na.rm=TRUE) + (1* SD_repeatAbund)
	thres_lower = mean(RepeatAbundance_chr$RepeatAbundance, na.rm=TRUE) - (1* SD_repeatAbund)
    # 138 796 651

    # find positions of + two SD from min
    cent_range_wind <- RepeatAbundance_chr$Genome_position[which(RepeatAbundance_chr$RepeatAbundance >= thres_upper | RepeatAbundance_chr$RepeatAbundance <= thres_lower )]/5000

    # find range of these values
    # ChemoSpecUtils
    # https://rdrr.io/cran/ChemoSpecUtils/man/check4Gaps.html

    #library(ChemoSpecUtils)

    cent_range <- ChemoSpecUtils::check4Gaps(cent_range_wind)
    cent_range[nrow(cent_range)+1,] <- c(0,0,0,0,0)

    cent_range_pos_start <- cent_range[,1]*5000
    cent_range_pos_end <- cent_range[,2]*5000

    SPP_l <- rep(fname, length(cent_range_pos_start))
    Chr_l <- rep(chromosome, length(cent_range_pos_start))
    typel <- rep("RepeatAbund", length(cent_range_pos_start))

    RepeatAbund_cent_chr <- cbind(typel, SPP_l, Chr_l, cent_range_pos_start, cent_range_pos_end)
    RepeatAbund_cent <- rbind(RepeatAbund_cent, RepeatAbund_cent_chr)

    RepeatAbund_min_chr <- c("RepeatAbund", fname, chromosome, cent_min)
    RepeatAbund_min <- rbind(RepeatAbund_min, RepeatAbund_min_chr)

    chr_length <- max(RepeatAbundance_chr$Genome_position)
    RepeatAbund_length_chr <- c("Length", fname, chromosome, chr_length)
    RepeatAbund_length <- rbind(RepeatAbund_length, RepeatAbund_length_chr)

  }

  # remove zeros
  RepeatAbund_cent <- RepeatAbund_cent[-which(RepeatAbund_cent[,4] == 0),]

  # output final file
  utils::write.table(x=RepeatAbund_cent, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname, "_RepeatAbund_centromere_range.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=RepeatAbund_min, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname, "_RepeatAbund_centromere_prediction_min.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=RepeatAbund_length, file=paste0(outpath,"/", fname,"/Summary_output/histograms/", fname, "_RepeatAbund_centromere_prediction_length.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)

  #-----------------------------------------------

  # get Shannon div data
  lsd <- lapply(Shannon_list, read.table)
  sd_chr_list0 <- basename(Shannon_list)
  sd_chr_list1 <- stringr::str_split(sd_chr_list0, "_", simplify =TRUE)
  sd_chr_list2 <- sd_chr_list1[,3]
  names(lsd) <- sd_chr_list2
  Shannon_div_total <- dplyr::bind_rows(lsd, .id = 'chromosome')
  colnames(Shannon_div_total) <- c("Chromosome", "Genome_position", "Shannon_div")

  Shannon_div_total$Chrnum <- as.numeric(stringr::str_split(Shannon_div_total$Chromosome, "r", simplify =TRUE)[,2])

  Shannon_div_total$Genome_position <- as.numeric(as.character(Shannon_div_total$Genome_position))
  Shannon_div_total$Shannon_div <- as.numeric(as.character(Shannon_div_total$Shannon_div))

  # define window
  bin_size=500
  # run rolling mean
  Shannon_div_total$roll_mean_Shannon <- zoo::rollapply(Shannon_div_total$Shannon_div, width = bin_size, FUN=mean, fill = NA, partial=(bin_size/2))

  # find start and end of highly repeating regions based 1SD from min
  Shannon_cent <- NULL
  Shannon_min <- NULL
  Shannon_length <- NULL
  for (chromosome in 1:length(Shannon_list)){

    Shannon_div_chr <- Shannon_div_total[which(Shannon_div_total$Chrnum == chromosome),]

    # find min position
    cent_min <- Shannon_div_chr$Genome_position[which(Shannon_div_chr$roll_mean_Shannon == min(Shannon_div_chr$roll_mean_Shannon, na.rm=TRUE))]

    # find SD of data
    SD_repeatAbund <- sd(Shannon_div_chr$roll_mean_Shannon, na.rm=TRUE)

    thres = min(Shannon_div_chr$roll_mean_Shannon, na.rm=TRUE) + (1 * SD_repeatAbund)

    # find positions of + two SD from min
    cent_range_wind <- Shannon_div_chr$Genome_position[which(Shannon_div_chr$roll_mean_Shannon <= thres)]/5000

    # find range of these values
    # ChemoSpecUtils
    # https://rdrr.io/cran/ChemoSpecUtils/man/check4Gaps.html

    #library(ChemoSpecUtils)

    cent_range <- ChemoSpecUtils::check4Gaps(cent_range_wind)
    cent_range[nrow(cent_range)+1,] <- c(0,0,0,0,0)

    cent_range_pos_start <- cent_range[,1]*5000
    cent_range_pos_end <- cent_range[,2]*5000

    SPP_l <- rep(fname, length(cent_range_pos_start))
    Chr_l <- rep(chromosome, length(cent_range_pos_start))
    typel <- rep("Shannon", length(cent_range_pos_start))

    Shannon_cent_chr <- cbind(typel, SPP_l, Chr_l, cent_range_pos_start, cent_range_pos_end)
    Shannon_cent <- rbind(Shannon_cent, Shannon_cent_chr)

    Shannon_min_chr <- c("Shannon", fname, chromosome, cent_min)
    Shannon_min <- rbind(Shannon_min, Shannon_min_chr)

    chr_length <- max(Shannon_div_chr$Genome_position)
    Shannon_length_chr <- c("Length", fname, chromosome, chr_length)
    Shannon_length <- rbind(Shannon_length, Shannon_length_chr)

  }

  # remove zeros
  Shannon_cent <- Shannon_cent[-which(Shannon_cent[,4] == 0),]

  # output final file
  utils::write.table(x=Shannon_cent, file=paste0(outpath,"/", fname,"/Summary_output/Shannon_div/", fname, "_Shannon_centromere_range.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=Shannon_min, file=paste0(outpath,"/", fname,"/Summary_output/Shannon_div/", fname, "_Shannon_centromere_prediction_min.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)
  utils::write.table(x=Shannon_length, file=paste0(outpath,"/", fname,"/Summary_output/Shannon_div/", fname, "_Shannon_centromere_prediction_length.txt"), sep = "\t", dec = ".",row.names = FALSE, col.names = FALSE)

}

#----------------------------------------
# Create documentations for functions above
# devtools::document()


