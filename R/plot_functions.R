# Plots out effect intervals for diallel data
#' @export
caterpillar.plot <- function(mcmc.object, 
                             name=NULL, 
                             include.effect.types=c("mu", "female", "add", "mat", "inbred", "epi", "var"),
                             col=c("black", "black", "dodgerblue2", "forestgreen", "darkorange2", "darkorchid1", "black"),
                             inbred.penalty.col="firebrick2",
                             manual.limits=NULL){
  # Processing the plotted variables
  reorder.names <- NULL
  include.col <- rep(FALSE, ncol(mcmc.object))
  use.color <- NULL
  for (i in 1:length(include.effect.types)){
    if (include.effect.types[i] == "var") {
      temp.include.col <- grepl(x=colnames(mcmc.object), pattern="tau|sigma")
    }
    else {
      temp.include.col <- grepl(x=colnames(mcmc.object), pattern=include.effect.types[i]) & !grepl(x=colnames(mcmc.object), pattern="tau|sigma")
    }
    if (include.effect.types[i] == "inbred") {
      use.color <- c(use.color, c(inbred.penalty.col, rep(col[i], sum(temp.include.col) - 1)))
    }
    else {
      use.color <- c(use.color, rep(col[i], sum(temp.include.col)))
    }
    reorder.names <- c(reorder.names, colnames(mcmc.object)[temp.include.col])
  }
  
  mcmc.object <- mcmc.object[,reorder.names]
  num.var <- ncol(mcmc.object)
  ci95.data <- coda::HPDinterval(mcmc.object, prob=0.95)
  ci50.data <- coda::HPDinterval(mcmc.object, prob=0.50)
  par.names <- rownames(ci95.data)
  means.data <- as.vector(apply(mcmc.object, 2, function(x) mean(x)))
  median.data <- as.vector(apply(mcmc.object, 2, function(x) median(x)))
  
  titlename = strsplit(deparse(substitute(mcmc.object)), ".", fixed=T)[1]
  if (!is.null(name)) {
    titlename=name
  }
  
  # Window limits
  if (!is.null(manual.limits)) {
    this.x.lim <- manual.limits
  }
  else {
    this.x.lim <- c(min(ci95.data, na.rm=TRUE), max(ci95.data, na.rm=TRUE))
  }

  plot(ci95.data[1, 1:2], c(1,1), panel.first = abline(h = 1, lty = 3, col = "gray88"), 
       type = "l", ylim = c(0, num.var+1), main = paste(titlename, "parameters"), 
       xlim = this.x.lim, xlab = "HPD intervals of strain effects and model parameters", 
       ylab = "", yaxt = "n", lty=1, lwd=1, col = use.color[1], frame.plot=F)
  if (length(ci95.data[1,]) > 2){
    for (i in seq(3, length(ci95.data[1,])-1, by=2)) {
      lines(ci95.data[1, c(i,i+1)], c(1,1), lty=1, lwd=1, col=use.color[1], lend=2)
    }
  }
  lines(ci50.data[1, 1:2], c(1,1), lty=1, lwd=5, col=use.color[1])
  if (length(ci50.data[1,]) > 2){
    for (i in seq(3, length(ci50.data[1,])-1, by=2)) {
      lines(ci50.data[1, c(i,i+1)], c(1,1), lty=1, lwd=3, col=use.color[1], lend=2)
    }
  }
  for (j in 2:num.var) {
    abline(h = j, lty = 3, col = "gray88")
    for (k in seq(1, length(ci95.data[j,])-1, by=2)) {
      lines(ci95.data[j, c(k,k+1)], c(j,j), lty=1, lwd=1, col=use.color[j], lend=2)
    }
    for (k in seq(1, length(ci50.data[j,])-1, by=2)) {
      lines(ci50.data[j, c(k,k+1)], c(j,j), lty=1, lwd=3, col=use.color[j], lend=2)
    }
  }
  points(median.data, 1:num.var, pch="l", col="white")
  points(means.data, 1:num.var, pch="l", col=use.color)
  abline(v=0, lty=2, col="gray")
  axis(2, c(1:num.var), par.names, c(1:num.var), las = 2, tck = -0.005, cex.axis = 0.5)
}

######################## Plot for all crosses
#' @export
diallelPlotter <- function(results, 
                           cross.type=c("f2", "bc", "rbc1", "rbc2"), 
                           pheno.name="", 
                           col.spectrum=c("blue2red", "gray", "green2red", "blue2green"),
                           path=NULL,
                           height=12,
                           width=12,
                           strains=c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB"),
                           off.x=TRUE,
                           include.biparent.labels=TRUE,
                           include.density=TRUE,
                           include.widgets=TRUE,
                           density.col="white",
                           border.col="black",
                           median.line.col="black",
                           absolute.density.scale=TRUE,
                           ...){
  cross.type <- cross.type[1]
  eu.list <- results$eu
  pheno.list <- results$pheno
  var.list <- results$var
  qtl.num <- results$qtl.num
  n <- results$n
  col.spectrum <- col.spectrum[1]
  
  absolute.max <- NULL
  if (absolute.density.scale) {
    for (i in which(grepl(x=names(eu.list), pattern=cross.type))) {
      absolute.max <- max(absolute.max,
                          max(hist(eu.list[[i]], plot=FALSE, breaks=seq(0, qtl.num, length.out=20))$density))
    }
  }
  
  ## Setting color spectrum
  spectrum <- make.spectrum(col.spectrum=col.spectrum, n=1000)
  
  # Placing labels
  label.indices <- 1:8
  if(pheno.name==""){
    pheno.name <- strsplit(deparse(substitute(eu.list)), split=".", fixed=T)[[1]][1]
  }
  if (!is.null(path)) {
    pdf(paste(path,"/Diallel_", pheno.name, "_", cross.type, "_", qtl.num, "qtl", sep = ""),
        width=12, height = 12)    
  }
  par(mfrow=c(8,8), 
      cex = 0.5, 
      oma=c(1,1,1,0), 
      mar=c(1,1,1,1))
  ## CHOOSE SUBSET OF DATA
  for (i in 1:length(strains)) {
    for (j in 1:length(strains)) {
      if (cross.type == "f2"){
        if(i == j){
          if(i == 1){
            infoPlotter(trait=pheno.name, 
                          experiment=cross.type,
                          n=n,
                          spectrum=spectrum, 
                          qtl.num=qtl.num,
                          include.widgets=include.widgets)
          }
          else if(i == 2){
            if (include.widgets) {
              legendPlotter()
            }
            else {
              emptyPlotter(off.x=off.x)
            }
          }
          else{
            emptyPlotter(off.x=off.x)
          }
        }
        else if(i > j){
          emptyPlotter(off.x=off.x)
        }      	
        else{
          oneParamPlotter(cross.u=eu.list[[paste(strains[i], "x", strains[j], "-f2_eu", sep="")]],
                          cross.type="f2", qtl.perc=median(var.list[[paste(strains[i], "x", strains[j], "-f2_perc", sep="")]]),
                          qtl.num=qtl.num, 
                          homo1.vec=pheno.list[[paste(strains[i], "x", strains[j], "-f2-hom1", sep="")]],
                          homo2.vec=pheno.list[[paste(strains[i], "x", strains[j], "-f2-hom2", sep="")]],
                          hetero.vec=pheno.list[[paste(strains[i], "x", strains[j], "-f2-het", sep="")]],
                          spectrum=spectrum,
                          include.density=include.density,
                          include.widgets=include.widgets,
                          absolute.max=absolute.max,
                          density.col=density.col,
                          border.col=border.col,
                          median.line.col=median.line.col,
                          ...)
        }
      }	
      else if (cross.type == "bc") {
        if (i == j) {
          if (i == 1) {
            infoPlotter(trait=pheno.name, 
                          experiment=cross.type,
                          n=n,
                          spectrum=spectrum, 
                          qtl.num=qtl.num,
                          include.widgets=include.widgets)
          }
          else if (i == 2) {
            if (include.widgets) {
              legendPlotter()
            }
            else {
              emptyPlotter(off.x=off.x)
            }
          }
          else{
            emptyPlotter(off.x=off.x)
          }
        }
        else{
          if(i < j){
            oneParamPlotter(eu.list[[paste(strains[i], "x", strains[j], "-bc1_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[i], "x", strains[j], "-bc1_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[i], "x", strains[j], "-bc1-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[i], "x", strains[j], "-bc1-het", sep="")]],
                            qtl.num=qtl.num, back.allele="B",
                            spectrum=spectrum,
                            absolute.max=absolute.max,
                            include.density=include.density,
                            include.widgets=include.widgets,
                            density.col=density.col,
                            border.col=border.col,
                            median.line.col=median.line.col,
                            ...)          				
          }
          else{
            oneParamPlotter(eu.list[[paste(strains[j], "x", strains[i], "-bc2_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[j], "x", strains[i], "-bc2_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[j], "x", strains[i], "-bc2-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[j], "x", strains[i], "-bc2-het", sep="")]],
                            qtl.num=qtl.num, back.allele="A",
                            spectrum=spectrum,
                            absolute.max=absolute.max,
                            include.density=include.density,
                            include.widgets=include.widgets,
                            density.col=density.col,
                            border.col=border.col,
                            median.line.col=median.line.col,
                            ...)   				
          }
        }
      }
      else if (cross.type == "rbc1") {
        if (i == j){
          if (i == 1) {
            infoPlotter(trait=pheno.name, 
                        experiment=cross.type,
                        n=n,
                        spectrum=spectrum, 
                        qtl.num=qtl.num,
                        include.widgets=include.widgets)
          }
          else if(i == 2) {
            if (include.widgets) {
              legendPlotter()
            }
            else {
              emptyRBCPlotter()
            }
          }
          
          else{
            emptyRBCPlotter()
          }  
        }
        else {
          if (i < j) {
            oneParamPlotter(eu.list[[paste(strains[i], "x", strains[j], "-rbc1_1_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[i], "x", strains[j], "-rbc1_1_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[i], "x", strains[j], "-rbc1_1-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[i], "x", strains[j], "-rbc1_1-het", sep="")]],
                            qtl.num=qtl.num, back.allele="A",
                            spectrum=spectrum,
                            absolute.max=absolute.max,
                            include.density=include.density,
                            include.widgets=include.widgets,
                            density.col=density.col,
                            border.col=border.col,
                            median.line.col=median.line.col,
                            ...)            			
          }
          else {
            oneParamPlotter(eu.list[[paste(strains[j], "x", strains[i], "-rbc1_2_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[j], "x", strains[i], "-rbc1_2_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[j], "x", strains[i], "-rbc1_2-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[j], "x", strains[i], "-rbc1_2-het", sep="")]],
                            qtl.num=qtl.num, back.allele="A",
                            spectrum=spectrum,
                            absolute.max=absolute.max,
                            include.density=include.density,
                            include.widgets=include.widgets,
                            density.col=density.col,
                            border.col=border.col,
                            median.line.col=median.line.col,
                            ...)   				
          }
        }
      }
      else if (cross.type == "rbc2") {
        if (i == j) {
          if (i == 1) {
            infoPlotter(trait=pheno.name, 
                        experiment=cross.type, 
                        n=n, 
                        spectrum=spectrum, 
                        qtl.num=qtl.num,
                        include.widgets=include.widgets)
          }
          else if (i == 2) {
            if (include.widgets) {
              legendPlotter()
            }
            else{
              emptyRBCPlotter()
            }  
          }
          else {
            emptyRBCPlotter()
          }
        }
        else {
          if (i < j) {
            oneParamPlotter(eu.list[[paste(strains[i], "x", strains[j], "-rbc2_1_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[i], "x", strains[j], "-rbc2_1_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[i], "x", strains[j], "-rbc2_1-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[i], "x", strains[j], "-rbc2_1-het", sep="")]],
                            qtl.num=qtl.num, back.allele="B",
                            spectrum=spectrum,
                            absolute.max=absolute.max,
                            include.density=include.density,
                            include.widgets=include.widgets,
                            density.col=density.col,
                            border.col=border.col,
                            median.line.col=median.line.col,
                            ...)              		
          }
          else {
            oneParamPlotter(eu.list[[paste(strains[j], "x", strains[i], "-rbc2_2_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[j], "x", strains[i], "-rbc2_2_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[j], "x", strains[i], "-rbc2_2-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[j], "x", strains[i], "-rbc2_2-het", sep="")]],
                            qtl.num=qtl.num, back.allele="B",
                            spectrum=spectrum,
                            absolute.max=absolute.max,
                            include.density=include.density,
                            include.widgets=include.widgets,
                            density.col=density.col,
                            border.col=border.col,
                            median.line.col=median.line.col,
                            ...)   				
          }
        }
      }
      if (j %in% label.indices & i == 1) {
        if (include.biparent.labels) {
          mtext(paste(strains[j], "(B)"), side=3, cex=1.1)
        }
        else {
          mtext(strains[j], side=3, cex=1.1)
        }
      }
      if(i %in% label.indices & j == 1){
        if (include.biparent.labels) {
          mtext(paste(strains[i], "(A)"), side=2, cex=1.1)
        }
        else {
          mtext(strains[j], side=3, cex=1.1)
        }
      }
    }
  } 
  if (!is.null(path)) {
    dev.off()
  }
}


################## Component plots of Moonrise plot
make.spectrum <- function(col.spectrum, n=1000) {
  if (col.spectrum == "gray") {
    spectrum <- gray(level=n:1/n)
  }
  else {
    spectrum <- do.call(what=eval(parse(text=paste0("colorRamps::", col.spectrum))), args=list(n=n))
  }
  return(spectrum)
}

emptyPlotter <- function(off.x=FALSE, ...){
  if (off.x) {
    this.col <- "white"
  }
  else {
    this.col <- "black"
  }
  plot(x=1, y=1, cex=50, pch=4, xlab="", ylab="",
       frame.plot=FALSE, xaxt="n", yaxt="n", col=this.col)
}

emptyF2Plotter <- function(background){
  plot(x=1, y=1, cex=1, pch=4, xlab="", ylab="", xlim=c(-5, 5), ylim=c(-5,5),
       frame.plot=FALSE, xaxt="n", yaxt="n", col="white")
  text(0, 3, "F2", cex=4)
  
  text(0, -2, labels=paste(background, " (A)", sep=""), cex=2)
}

emptyBCPlotter <- function(background){
  plot(x=1,y=1, cex = 1, pch ="", xlab = "", ylab = "", xlim=c(-10, 10), ylim=c(-10,10),
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")
  arrows(5, 0, -5, 0, length=0.1)
  text(0, -3, labels=paste("Background: ", background, " (A)", sep=""), cex=1.5)
}

emptyRBCPlotter <- function(){
  plot(x=1,y=1, cex = 1, pch ="", xlab = "", ylab = "", xlim=c(-10, 10), ylim=c(-10,10),
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")

  arrows(3, 3, -3 , 3, length=0.1)
  text(0, -1, "Maternal Strain", cex=1.5)
}

legendPlotter <- function(){
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", frame=FALSE, xaxt="n", yaxt="n")
  
  legend("center", legend=c("QTL", "Noise"), 
         fill=c("gray50", "white"), border="black", bty="n", 
         title="Variation type", cex=1.3, pt.cex=1.3) 
}

infoPlotter <- function(trait,
                        experiment,
                        n, 
                        spectrum, 
                        qtl.num=1,
                        include.widgets=TRUE){
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", frame=FALSE, xaxt="n", yaxt="n")
  
  if (experiment == "f2") { this.experiment <- "F2" }
  if (experiment == "bc") { this.experiment <- "BC" }
  if (experiment == "rbc1") { this.experiment <- "RBC (A)" }
  if (experiment == "rbc2") { this.experiment <- "RBC (B)" }
  
  text(x=0.1, y=0.9, labels=paste("Cross type:", this.experiment), adj=0, cex=1.2)
  text(x=0.1, y=0.8, labels=paste("Trait:", trait), adj=0, cex=1.2)
  text(x=0.1, y=0.7, labels=paste("QTL number:", qtl.num), adj=0, cex=1.2)
  text(x=0.1, y=0.6, labels=paste("Number of mice:", n), adj=0, cex=1.2)
  
  barplot(height=rep(0.1, length(spectrum)), width=1/length(spectrum), density=1000,
          angle=90, col=spectrum, 
          border=FALSE, space=FALSE, axes=FALSE, add=TRUE)
  axis(1, at=0:qtl.num, labels=0:qtl.num, tick=TRUE, cex.axis=1.2)
  text(x=0.5, 0.2, labels="Posterior mean utility", cex=1.3)
  if (include.widgets) {
       
  }
}

f2boxPlotter <- function(homo1.vec, 
                         homo2.vec, 
                         hetero.vec, 
                         y.max, 
                         x.max,
                         border.col){
  max.box.y <- max(homo1.vec, homo2.vec, hetero.vec, na.rm=TRUE)
  min.box.y <- min(homo1.vec, homo2.vec, hetero.vec, na.rm=TRUE)
  box.y.range <- max.box.y - min.box.y
  # Scaling and rescaling
  homo1.vec.sc <- (homo1.vec - min.box.y)*(y.max*(3/5)/(max.box.y-min.box.y)) + y.max*(5/4)
  homo2.vec.sc <- (homo2.vec - min.box.y)*(y.max*(3/5)/(max.box.y-min.box.y)) + y.max*(5/4)
  hetero.vec.sc <- (hetero.vec - min.box.y)*(y.max*(3/5)/(max.box.y-min.box.y)) + y.max*(5/4)
  mid.x <- (2/5)*(3/8)*x.max
  shift.x <- (3/4)*mid.x
  boxplot(homo1.vec.sc, 
          hetero.vec.sc, 
          homo2.vec.sc,
          width=rep(1, 3), 
          boxwex=(1/15)*x.max,
          border=border.col,
          outline=FALSE, 
          col="white",
          add=TRUE, 
          at=c(mid.x-shift.x, mid.x, mid.x+shift.x), names=NA, yaxt="n", xaxt="n")
  text(x=mid.x-shift.x, y=y.max*(5/4), labels="A/A", cex=0.8, col=border.col)
  text(x=mid.x, y=y.max*(5/4), labels="A/B", cex=0.8, col=border.col)
  text(x=mid.x+shift.x, y=y.max*(5/4), labels="B/B", cex=0.8, col=border.col)
  text(x=mid.x+(1/2)*shift.x, y=2*y.max-(1/15)*y.max, labels="Phenotypes", cex=1.1, col=border.col)
}

bcboxPlotter <- function(homo.vec, 
                         hetero.vec, 
                         y.max, 
                         x.max, 
                         back.allele="A",
                         border.col){
  max.box.y <- max(homo.vec, hetero.vec, na.rm=TRUE)
  min.box.y <- min(homo.vec, hetero.vec, na.rm=TRUE)
  box.y.range <- max.box.y - min.box.y
  # Scaling and rescaling
  homo.vec.sc <- (homo.vec - min.box.y)*(y.max*(4/7)/(max.box.y-min.box.y)) + y.max*(4/3)
  hetero.vec.sc <- (hetero.vec - min.box.y)*(y.max*(4/7)/(max.box.y-min.box.y)) + y.max*(4/3)
  mid.x <- (2/5)*(3/8)*x.max
  shift.x <- (3/4)*mid.x
  boxplot(homo.vec.sc, 
          hetero.vec.sc,
          width=rep(1, 2), 
          boxwex=(1/15)*x.max, 
          outline=FALSE, 
          col="white",
          add=TRUE, 
          border=border.col,
          at=c(mid.x-shift.x, mid.x+shift.x), names=NA, yaxt="n", xaxt="n")
  if (back.allele=="A") {
    text(x=mid.x-shift.x, y=y.max*(5/4), labels="A/A", cex=0.8, col=border.col)
  }
  else {
    text(x=mid.x-shift.x, y=y.max*(5/4), labels="B/B", cex=0.8, col=border.col)
  }
  text(x=mid.x+shift.x, y=y.max*(5/4), labels="A/B", cex=0.8, col=border.col)
  text(x=mid.x+(1/2)*shift.x, y=2*y.max-(1/15)*y.max, labels="Phenotypes", cex=1.1, col=border.col)
}

oneParamPlotter <- function(cross.u, 
                            cross.type, 
                            qtl.perc, 
                            qtl.num=1, 
                            cross.label1="", 
                            cross.label2="",
                            homo1.vec, 
                            homo2.vec=NULL, 
                            hetero.vec, 
                            back.allele=NULL,
                            include.x.axis=FALSE,
                            spectrum,
                            absolute.max=NULL,
                            include.density,
                            include.widgets,
                            density.col="gray",
                            border.col="black",
                            median.line.col="black",
                            ...){
  x.high <- qtl.num
  post.mean <- mean(cross.u)
  post.median <- median(cross.u)
  bgcolor <- spectrum[round((post.mean/x.high)*length(spectrum))]
  n <- length(cross.u)
  if (is.null(absolute.max)) {
    max.y <- max(hist(cross.u, plot=FALSE, breaks=seq(0, x.high, length.out=20))$density)
  }
  else {
    max.y <- absolute.max
  }
  
  min.y <- 0
  
  ## Adjusting plot window for widgets
  if (include.widgets) {
    y.max <- 1*(max.y - min.y) + max.y
  }
  else {
    y.max <- max.y
  }

  plot(1, 1, xlim=c(0, x.high), 
       ylim=c(0, y.max), 
       lwd=3, xlab="", main="", frame.plot=FALSE, ylab="", yaxt="n", xaxt="n", col="white")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
  y.range <- max.y - min.y
  qtl <- round(qtl.perc/qtl.num, 1)
  total.qtl <- round(qtl.perc, 1)
  total <- round(100 - total.qtl, 1)

  if (include.widgets) {
    if (qtl.num == 1) {
      mapplots::add.pie(z=c(qtl/100, total/100), labels=NA, x=0.5*x.high, y=0.75*y.max, radius=(3/10)*(9/10)*y.range, 
                        col=(c("gray50", "white")), cex=1.3, label.dist=1.2, border=border.col)
    }
    else {
      excess.qtl <- round(qtl*(qtl.num-1), 1)
      mapplots::add.pie(z=c(qtl/100, excess.qtl/100, total/100), labels=NA, x=0.5*x.high, y=0.75*y.max, radius=(3/10)*(9/10)*y.range, 
                        col=(c("gray50", "gray50", "white")), cex=1.3, label.dist=1.2, border=border.col)  
    }
    legend("topright", 
           legend=c(qtl, total), 
           fill=c("gray50", "white"), 
           title="% Variance", 
           border=border.col, 
           bty="n", 
           text.col=border.col, 
           cex=1.2)
    
    if (cross.type == "f2") {
      f2boxPlotter(homo1.vec=homo1.vec, 
                   homo2.vec=homo2.vec, 
                   hetero.vec=hetero.vec, 
                   y.max=max.y, 
                   x.max=x.high,
                   border.col=border.col)
    }
    else if (cross.type == "bc") {
      bcboxPlotter(homo.vec=homo1.vec, 
                   hetero.vec=hetero.vec, 
                   y.max=max.y, 
                   x.max=x.high, 
                   back.allele=back.allele,
                   border.col=border.col)
    }
  }
  if (include.density) {
    hist(cross.u, col=density.col, breaks=seq(0, x.high, by=0.05), ylim=c(0, y.max), xlim=c(0, x.high), add=TRUE, freq=FALSE)
    lines(x=c(post.median, post.median), y=c(0, max.y), lty=5, lwd=2, col=median.line.col)
  }
  if (cross.label1 != "" & cross.label2 != "") {
    text(1, median(dens$y), paste(cross.label1, 'x', cross.label2, sep=""))
  }
  if (include.x.axis) {
    axis(1, at=0:x.high, labels=TRUE, tick=TRUE, cex.axis=1.2)
  }
}

