# Plots out effect intervals for diallel data
#' @export
caterpillar.plot <- function(mcmc.object, name=NULL, full=FALSE)
{
  my.colors <- c(rep("black", 3), rep("blue", 8), rep("deepskyblue", 8), rep("green", 8), rep("cyan", 28), rep("black", 5))
  # Plotting only the random effects
  if(full == F){
    reduced.col <- -c(1, 2, 3, 56:ncol(mcmc.object))
    mcmc.object <- mcmc.object[,reduced.col]
    my.colors <- my.colors[reduced.col]
  }
  no.var <- ncol(mcmc.object)
  ci95.data <- HPDinterval(mcmc.object, prob=0.95)
  ci50.data <- HPDinterval(mcmc.object, prob=0.50)
  par.names <- rownames(ci95.data)
  means.data <- as.vector(apply(mcmc.object, 2, function(x) mean(x)))
  median.data <- as.vector(apply(mcmc.object, 2, function(x) median(x)))
  
  titlename = strsplit(deparse(substitute(mcmc.object)), ".", fixed=T)[1]
  if(!is.null(name))
  {
    titlename=name
  }
  plot(ci95.data[1, 1:2], c(1,1), panel.first = abline(h = 1, lty = 3, col = "gray88"), 
       type = "l", ylim = c(0, no.var+1), main = c("Caterpillar plot of parameters", paste("for ", titlename)), 
       xlim = c(min(ci95.data), max(ci95.data, na.rm=T)), xlab = "HPD intervals of strain effects and parameters", 
       ylab = "", yaxt = "n", lty=1, lwd=1, col = my.colors[1], frame.plot=F)
  if(length(ci95.data[1,]) > 2){
    for(i in seq(3, length(ci95.data[1,])-1, by=2)){
      lines(ci95.data[1, c(i,i+1)], c(1,1), lty=1, lwd=1, col=my.colors[1])
    }
  }
  lines(ci50.data[1, 1:2], c(1,1), lty=1, lwd=5, col=my.colors[1])
  if(length(ci50.data[1,]) > 2){
    for(i in seq(3, length(ci50.data[1,])-1, by=2)){
      lines(ci50.data[1, c(i,i+1)], c(1,1), lty=1, lwd=5, col=my.colors[1])
    }
  }
  for(j in 2:no.var){
    abline(h = j, lty = 3, col = "gray88")
    for(k in seq(1, length(ci95.data[j,])-1, by=2)){
      lines(ci95.data[j, c(k,k+1)], c(j,j), lty=1, lwd=1, col=my.colors[j])
    }
    for(k in seq(1, length(ci50.data[j,])-1, by=2)){
      lines(ci50.data[j, c(k,k+1)], c(j,j), lty=1, lwd=5, col=my.colors[j])
    }
  }
  points(median.data, 1:no.var, pch="l", col="white")
  points(means.data, 1:no.var, pch="l", col=my.colors)
  abline(v=0, lty=2, col="red")
  axis(2, c(1:no.var), par.names, c(1:no.var), las = 2, tck = -0.005, cex.axis = 0.5)
}

################## Component plots of Moonrise plot
emptyPlotter <- function(...){
  plot(x=1,y=1, cex = 50, pch = 4, xlab = "", ylab = "",
       frame.plot = FALSE, xaxt = "n", yaxt = "n")
}

emptyF2Plotter <- function(background){
  plot(x=1,y=1, cex = 1, pch = 4, xlab = "", ylab = "", xlim=c(-5, 5), ylim=c(-5,5),
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")
  text(0, 3, "F2", cex=4)
  
  text(0, -2, labels=paste(background, " (A)", sep=""), cex=2)
}

bclegendPlotter <- function(background, trait, n, n.color, qtl.num=1){
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", frame=F, xaxt="n", yaxt="n")
  
  text(x=0, y=1, labels=paste("Trait: ", trait, sep=""), adj=0, cex=1.2)
  text(x=0, y=0.9, labels=paste("QTL number: ", qtl.num, sep=""), adj=0, cex=1.2)
  text(x=0, y=0.8, labels=paste("Number of mice: ", n, sep=""), adj=0, cex=1.2)
  
  arrows(0.4, 0.55, 0.1, 0.55, length=0.05)
  text(0.25, 0.45, labels=paste(background, " (A)", sep=""), cex=1.1)
  
  palette(gplots::rich.colors(n + 300))
  x <- seq(0, 1, length = n.color)
  grad.col <- gplots::rich.colors(n.color + 300)[-c(1:300)]
  barplot(height=rep(0.1, n.color), width=1/700, density=1000, angle=90, col = grad.col, border = F, space = FALSE, axes = F, add=T)
  axis(1, at=0:qtl.num, labels=0:qtl.num, tick=T, cex.axis=1.2)
  text(x= 0.5, 0.2, labels="Posterior mean utility", cex=1.3)
  legend("right", legend=c("QTL", "Noise"), 
         fill=c("gray50", "white"), border="black", bty="n", 
         title="Variation type", cex=1.3, pt.cex=1.3)
}

emptyBCPlotter <- function(background){
  plot(x=1,y=1, cex = 1, pch ="", xlab = "", ylab = "", xlim=c(-10, 10), ylim=c(-10,10),
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")
  arrows(5, 0, -5, 0, length=0.1)
  text(0, -3, labels=paste("Background: ",background, " (A)", sep=""), cex=1.5)
}

emptyRBC1Plotter <- function(background){
  plot(x=1,y=1, cex = 1, pch ="", xlab = "", ylab = "", xlim=c(-10, 10), ylim=c(-10,10),
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")
  arrows(-3, 7.5, 3, 7.5, length=0.1)
  arrows(-3, 7.5, -3, 1.5, length=0.1)
  text(0, 9.5, paste("Background: ",background, " (A)", sep=""), cex=1.5)
  
  arrows(3, -6, -3 , -6, length=0.1)
  text(0, -8, "Maternal Strain", cex=1.5)
}

emptyRBC2Plotter <- function(background){
  plot(x=1,y=1, cex = 1, pch ="", xlab = "", ylab = "", xlim=c(-10, 10), ylim=c(-10,10),
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")
  arrows(3, 4, -3, 4, length=0.1)
  arrows(3, 4, 3, 10, length=0.1)
  text(0, 2, paste("Background: ",background, " (B)", sep=""), cex=1.5)
  
  arrows(3, -6, -3 , -6, length=0.1)
  text(0, -8, "Maternal Strain", cex=1.5)
}

legendPlotter <- function(trait, n, n.color, qtl.num=1){
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", frame=F, xaxt="n", yaxt="n")
  
  text(x=0, y=1, labels=paste("Trait: ", trait, sep=""), adj=0, cex=1.2)
  text(x=0, y=0.9, labels=paste("QTL number: ", qtl.num, sep=""), adj=0, cex=1.2)
  text(x=0, y=0.8, labels=paste("Number of mice: ", n, sep=""), adj=0, cex=1.2)
  
  palette(gplots::rich.colors(n.color + 300))
  x <- seq(0, 1, length = n.color)
  grad.col <- gplots::rich.colors(n.color + 300)[-c(1:300)]
  barplot(height=rep(0.1, n.color), width=1/700, density=1000, angle=90, col = grad.col, border = F, space = FALSE, axes = F, add=T)
  axis(1, at=0:qtl.num, labels=0:qtl.num, tick=T, cex.axis=1.2)
  text(x= 0.5, 0.2, labels="Posterior mean utility", cex=1.3)
  legend("center", legend=c("QTL", "Noise"), 
         fill=c("gray50", "white"), border="black", bty="n", 
         title="Variation type", cex=1.3, pt.cex=1.3)
}

f2boxPlotter <- function(homo1.vec, homo2.vec, hetero.vec, y.max, x.max)
{
  max.box.y <- max(homo1.vec, homo2.vec, hetero.vec)
  min.box.y <- min(homo1.vec, homo2.vec, hetero.vec)
  box.y.range <- max.box.y - min.box.y
  # Scaling and rescaling
  homo1.vec.sc <- (homo1.vec - min.box.y)*(y.max*(3/5)/(max.box.y-min.box.y)) + y.max*(5/4)
  homo2.vec.sc <- (homo2.vec - min.box.y)*(y.max*(3/5)/(max.box.y-min.box.y)) + y.max*(5/4)
  hetero.vec.sc <- (hetero.vec - min.box.y)*(y.max*(3/5)/(max.box.y-min.box.y)) + y.max*(5/4)
  mid.x <- (2/5)*(3/8)*x.max
  shift.x <- (3/4)*mid.x
  boxplot(homo1.vec.sc, hetero.vec.sc, homo2.vec.sc,
          width=rep(1, 3), boxwex=(1/15)*x.max, outline=F, col="white",
          add=T, at=c(mid.x-shift.x, mid.x, mid.x+shift.x), names=NA, yaxt="n", xaxt="n")
  text(x=mid.x-shift.x, y=y.max*(5/4), labels="A/A", cex=0.8)
  text(x=mid.x, y=y.max*(5/4), labels="A/B", cex=0.8)
  text(x=mid.x+shift.x, y=y.max*(5/4), labels="B/B", cex=0.8)
  text(x=mid.x+(1/2)*shift.x, y=2*y.max-(1/15)*y.max, labels="Phenotypes", cex=1.1)
}

bcboxPlotter <- function(homo.vec, hetero.vec, y.max, x.max, back.allele="A"){
  max.box.y <- max(homo.vec, hetero.vec)
  min.box.y <- min(homo.vec, hetero.vec)
  box.y.range <- max.box.y - min.box.y
  # Scaling and rescaling
  homo.vec.sc <- (homo.vec - min.box.y)*(y.max*(4/7)/(max.box.y-min.box.y)) + y.max*(4/3)
  hetero.vec.sc <- (hetero.vec - min.box.y)*(y.max*(4/7)/(max.box.y-min.box.y)) + y.max*(4/3)
  mid.x <- (2/5)*(3/8)*x.max
  shift.x <- (3/4)*mid.x
  boxplot(homo.vec.sc, hetero.vec.sc,
          width=rep(1, 2), boxwex=(1/15)*x.max, outline=F, col="white",
          add=T, at=c(mid.x-shift.x, mid.x+shift.x), names=NA, yaxt="n", xaxt="n")
  if(back.allele=="A"){
    text(x=mid.x-shift.x, y=y.max*(5/4), labels="A/A", cex=0.8)
  }
  else{
    text(x=mid.x-shift.x, y=y.max*(5/4), labels="B/B", cex=0.8)
  }
  text(x=mid.x+shift.x, y=y.max*(5/4), labels="A/B", cex=0.8)
  text(x=mid.x+(1/2)*shift.x, y=2*y.max-(1/15)*y.max, labels="Phenotypes", cex=1.1)
}

oneParamPlotter <- function(cross.u, cross.type, qtl.perc, qtl.num=1, cross.label1="", cross.label2="",
                            homo1.vec, homo2.vec=NULL, hetero.vec, back.allele=NULL){
  x.high <- qtl.num
  palette(gplots::rich.colors(1000))
  post.mean <- mean(cross.u)
  post.median <- median(cross.u)
  bgcolor <- round(post.mean*(700/x.high) + 300)
  n <- length(cross.u)
  dens <- density(cross.u, from=0, to=x.high)
  max.y <- max(dens$y)
  min.y <- 0
  y.max <- 1*(max.y - min.y) + max.y
  min.x <- which.min(dens$x)
  max.x <- which.max(dens$x)
  plot(dens$x, dens$y, xlim=c(0, x.high), ylim=c(0, y.max), lwd=3, xlab="", main="", frame.plot=F, ylab="", yaxt="n", xaxt="n", col="gray")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=bgcolor)
  y.range <- max.y - min.y
  qtl <- round(qtl.perc/qtl.num, 1)
  total.qtl <- round(qtl.perc, 1)
  total <- round(100 - total.qtl, 1)

  if(qtl.num == 1){
    add.pie(z=c(qtl/100, total/100), labels=NA, x=0.5*x.high, y=0.75*y.max, radius=(3/10)*(9/10)*y.range, 
            col=(c("gray50", "white")), cex=1.3, label.dist=1.2, border="black")
  }
  else{
    excess.qtl <- round(qtl*(qtl.num-1), 1)
    add.pie(z=c(qtl/100, excess.qtl/100, total/100), labels=NA, x=0.5*x.high, y=0.75*y.max, radius=(3/10)*(9/10)*y.range, 
            col=(c("gray50", "gray50", "white")), cex=1.3, label.dist=1.2, border="black")  
  }
  legend("topright", legend=c(qtl, total), fill=c("gray50", "white"), title="% Variance", border="black", bty="n", text.col="black", cex=1.2)
  with(dens, polygon(x=c(x[c(min.x, min.x:max.x, max.x)]), y=c(0, y[min.x:max.x], 0), col="gray", border="black"))
  lines(x=c(post.median, post.median), y=c(0, max.y), lty=5, lwd=2, col="white")
  if(cross.label1 != "" & cross.label2 != ""){
    text(1, median(dens$y), paste(cross.label1, 'x', cross.label2, sep=""))
  }
  if(cross.type=="f2"){
    f2boxPlotter(homo1.vec=homo1.vec, homo2.vec=homo2.vec, hetero.vec=hetero.vec, y.max=max.y, x.max=x.high)
  }
  else if(cross.type=="bc"){
    bcboxPlotter(homo.vec=homo1.vec, hetero.vec=hetero.vec, y.max=max.y, x.max=x.high, back.allele=back.allele)
  }
  axis(1, at=0:x.high, labels=T, tick=T, cex.axis=1.2)
}

######################## Plot for all crosses
#' @export
diallelPlotter <- function(results, cross.type, pheno.name="", path=getwd(), qtl.num=1,
                           strains <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")){
  eu.list <- results$eu
  pheno.list <- results$pheno
  var.list <- results$var
  qtl.num <- results$qtl.num
  n <- results$n
  # Placing labels
  label.indices <- 1:8
  if(pheno.name==""){
    pheno.name <- strsplit(deparse(substitute(eu.list)), split=".", fixed=T)[[1]][1]
  }
  if(qtl.num == 1){
    pdf(paste(path,"/Diallel_", pheno.name, "_", cross.type, ".pdf",sep = ""),
        width=12, height = 12)
  }
  else{
    pdf(paste(path,"/Diallel_", pheno.name, "_", cross.type, "_", qtl.num, "qtl", sep = ""),
        width=12, height = 12)    
  }
  par(mfrow=c(8,8), cex = 0.5, oma=c(1,1,1,0), mar=c(1,1,1,1))
  ## CHOOSE SUBSET OF DATA
  for(i in 1:length(strains)){
    for(j in 1:length(strains)){
      if(cross.type == "f2"){
        if(i == j){
          if(i == length(strains)){
            legendPlotter(trait=pheno.name, n=700, n.color=700, qtl.num=qtl.num)
          }
          else{
            emptyF2Plotter(background=strains[i])
          }
        }
        else if(i > j){
          emptyPlotter()
        }      	
        else{
          oneParamPlotter(cross.u=eu.list[[paste(strains[i], "x", strains[j], "-f2_eu", sep="")]],
                          cross.type="f2", qtl.perc=median(var.list[[paste(strains[i], "x", strains[j], "-f2_perc", sep="")]]),
                          qtl.num=qtl.num, 
                          homo1.vec=pheno.list[[paste(strains[i], "x", strains[j], "-f2-hom1", sep="")]],
                          homo2.vec=pheno.list[[paste(strains[i], "x", strains[j], "-f2-hom2", sep="")]],
                          hetero.vec=pheno.list[[paste(strains[i], "x", strains[j], "-f2-het", sep="")]])
        }
      }	
      else if(cross.type == "bc"){
        if(i == j){
          if(i == 1){
            bclegendPlotter(background=strains[i], pheno.name, n=n, n.color=700, qtl.num=qtl.num)
          }
          else{
            emptyBCPlotter(background=strains[i])
          }
        }
        else{
          if(i < j){
            oneParamPlotter(eu.list[[paste(strains[i], "x", strains[j], "-bc1_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[i], "x", strains[j], "-bc1_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[i], "x", strains[j], "-bc1-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[i], "x", strains[j], "-bc1-het", sep="")]],
                            qtl.num=qtl.num, back.allele="B")          				
          }
          else{
            oneParamPlotter(eu.list[[paste(strains[j], "x", strains[i], "-bc2_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[j], "x", strains[i], "-bc2_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[j], "x", strains[i], "-bc2-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[j], "x", strains[i], "-bc2-het", sep="")]],
                            qtl.num=qtl.num, back.allele="A")   				
          }
        }
      }
      else if(cross.type == "rbc1"){
        if(i == j){
          if(i != length(strains)){
            emptyRBC1Plotter(background=strains[i])
          }
          else{
            legendPlotter(trait=pheno.name, n=n, n.color=700, qtl.num=qtl.num)
          }
        }
        else{
          if(i < j){
            oneParamPlotter(eu.list[[paste(strains[i], "x", strains[j], "-rbc1_1_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[i], "x", strains[j], "-rbc1_1_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[i], "x", strains[j], "-rbc1_1-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[i], "x", strains[j], "-rbc1_1-het", sep="")]],
                            qtl.num=qtl.num, back.allele="A")            			
          }
          else{
            oneParamPlotter(eu.list[[paste(strains[j], "x", strains[i], "-rbc1_2_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[j], "x", strains[i], "-rbc1_2_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[j], "x", strains[i], "-rbc1_2-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[j], "x", strains[i], "-rbc1_2-het", sep="")]],
                            qtl.num=qtl.num, back.allele="A")   				
          }
        }
      }
      else if(cross.type == "rbc2"){
        if(i == j){
          if(i == 1){
            legendPlotter(trait=pheno.name, n=n, n.color=700, qtl.num=qtl.num)
          }
          else{
            emptyRBC2Plotter(background=strains[i])
          }
        }
        else{
          if(i < j){
            oneParamPlotter(eu.list[[paste(strains[i], "x", strains[j], "-rbc2_1_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[i], "x", strains[j], "-rbc2_1_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[i], "x", strains[j], "-rbc2_1-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[i], "x", strains[j], "-rbc2_1-het", sep="")]],
                            qtl.num=qtl.num, back.allele="B")              		
          }
          else{
            oneParamPlotter(eu.list[[paste(strains[j], "x", strains[i], "-rbc2_2_eu", sep="")]],
                            cross.type="bc", qtl.perc=median(var.list[[paste(strains[j], "x", strains[i], "-rbc2_2_perc", sep="")]]),
                            homo1.vec=pheno.list[[paste(strains[j], "x", strains[i], "-rbc2_2-hom", sep="")]],
                            hetero.vec=pheno.list[[paste(strains[j], "x", strains[i], "-rbc2_2-het", sep="")]],
                            qtl.num=qtl.num, back.allele="B")   				
          }
        }
      }
      if(j %in% label.indices & i == 1){
        mtext(strains[j], side=3, cex=1.25)
      }
      if(i %in% label.indices & j == 1){
        mtext(strains[i], side=2, cex=1.25)
      }
    }
  } 
  dev.off()
}

