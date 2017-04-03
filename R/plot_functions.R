# Plots out effect intervals for diallel data
#' @export
caterpillar.plot <- function(mcmc.object, name=NULL, full=FALSE){
  my.colors <- c(rep("black", 3), rep("blue", 8), rep("deepskyblue", 8), rep("cyan", 8), rep("green", 28), rep("black", 5))
  # Plotting only the random effects
  if(full == FALSE){
    mcmc.object <- mcmc.object[,-c(1, 2, 3, 56, 57, 58, 59, 60)]
    my.colors <- my.colors[-c(1, 2, 3, 56, 57, 58, 59, 60)]
  }
  no.var <- ncol(mcmc.object)
  ci95.data <- HPDinterval(mcmc.object, prob=0.95)
  ci50.data <- HPDinterval(mcmc.object, prob=0.50)
  par.names <- rownames(ci95.data)
  means.data <- as.vector(apply(mcmc.object, 2, function(x) mean(x)))
  median.data <- as.vector(apply(mcmc.object, 2, function(x) median(x)))
  
  titlename = strsplit(deparse(substitute(mcmc.object)), ".", fixed=T)[1]
  if(!is.null(name)){
    titlename=name
  }
  plot(ci95.data[1, 1:2], c(1,1), panel.first = abline(h = 1, lty = 3, col = "gray88"), type = "l", ylim = c(0, no.var+1), main = c("Caterpillar plot of parameters", paste("for ", titlename)), xlim = c(min(ci95.data), max(ci95.data, na.rm=T)), xlab = "HPD intervals of strain effects and parameters", ylab = "", yaxt = "n", lty=1, lwd=1, col = my.colors[1])
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

emptyF2Plotter <- function(){
  plot(x=1,y=1, cex = 1, pch = 4, xlab = "", ylab = "",
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")
  text(1, 1, "F2", cex=5)
}

infoPlotter <- function(trait, n, add.qtl.num=1, dom.qtl.num=1){
  plot(x=1,y=1, xlim=c(-1, 1), ylim=c(-1, 1), cex = 50, pch = 4, xlab = "", ylab = "",
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")
  text(x=-0.95, y=0.7, labels=paste("Trait: ", trait, sep=""), adj=0, cex=1.2)
  text(x=-0.95, y=0.2, labels=paste("Additive QTL number: ", add.qtl.num, sep=""), adj=0, cex=1.2)
  text(x=-0.95, y=-0.3, labels=paste("Inbred QTL number: ", dom.qtl.num, sep=""), adj=0, cex=1.2)
  text(x=-0.95, y=-0.8, labels=paste("Number of mice: ", n, sep=""), adj=0, cex=1.2)
}

emptyBCPlotter <- function(){
  plot(x=1,y=1, cex = 1, pch ="", xlab = "", ylab = "",
       frame.plot = FALSE, xaxt = "n", yaxt = "n", col="white")
  text(1, 0.9, "BC Parent", cex=2)
  arrows(1.25, 1.1, 0.75, 1.1, length=0.1)
}

legendPlotter <- function(n, add.num=1, inbred.num=1){
  plot(NA, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", frame=F, xaxt="n", yaxt="n")
  palette(gplots::rich.colors(n + 300))
  x <- seq(0, 1, length = n)
  grad.col <- gplots::rich.colors(n + 300)[-c(1:300)]
  barplot(height=rep(0.25, n), width=1/700, density=1000, angle=90, col = grad.col, border = F, space = FALSE, axes = F, add=T)
  axis(1, at=0:(add.num+inbred.num)/(add.num+inbred.num), labels=0:(add.num+inbred.num), tick=T, cex.axis=1.2)
  text(x= 0.5, 0.35, labels="Posterior mean utility", cex=1.5)
  legend("top", legend=c("Additive", "Inbred", "Noise"), 
         fill=c("gray50", "gray20", "white"), border="black", bty="n", 
         title="QTL type", cex=1.5, pt.cex=1.5)
}

effectPlotter <- function(add.effect, inbred.effect, x.range, y.range)
{
  if(length(inbred.effect) == 2){
    inbred.effect <- inbred.effect[1] + inbred.effect[2]
  }
  effects <- c(add.effect, inbred.effect)
  max.abs.ef <- max(abs(effects))
  scale.factor <- max.abs.ef/(0.15*x.range)
  scaled.effects <- effects/scale.factor
  
  # Add line
  add.shift <- 0.155*x.range + scaled.effects[1]
  lines(x=c(0.155*x.range, add.shift), y=c(0.974*y.range, 0.974*y.range), lend=1, ljoin=2, lwd=5, col="gray50")
  #text(x=add.shift + sign(scaled.effects[1])*0.05*x.range, y=c(0.925*y.range, 0.925*y.range), labels=round(add.effect, 1))
  # Dom line
  lines(x=c(0.155*x.range, 0.155*x.range + scaled.effects[2]), y=c(0.924*y.range, 0.924*y.range), lend=1, ljoin=2, lwd=5, col="gray20")

  # Vertical line
  lines(x=c(0.155*x.range, 0.155*x.range), y=c(0.899*y.range, 0.999*y.range), lwd=1, lend=1, col="black")
}

oneParamPlotter <- function(cross.u, add.effect, inbred.effect, add.perc, inbred.perc, add.num=1, inbred.num=1, cross.label1="", cross.label2="")
{
  x.high <- add.num + inbred.num
  palette(gplots::rich.colors(1000))
  post.mean <- mean(cross.u)
  post.median <- median(cross.u)
  bgcolor <- round(post.mean*(700/x.high) + 300)
  
  dens <- density(cross.u)
  max.y <- max(dens$y)
  min.y <- min(dens$y)
  y.max <- 1*(max.y - min.y) + max.y
  min.x <- which.min(dens$x)
  max.x <- which.max(dens$x)
  plot(dens, xlim=c(0, x.high), ylim=c(0, y.max), lwd=3, xlab="", main="", frame.plot=F, ylab="", yaxt="n", xaxt="n", col="gray")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=bgcolor)
  y.range <- max.y - min.y
  add <- round(add.perc/add.num, 1)
  inbred <- round(inbred.perc/inbred.num, 1)
  total.add <- round(add.perc, 1)
  total.inbred <- round(inbred.perc, 1)
  total <- round(100 - total.add - total.inbred, 1)

  if(add.num == 1 & inbred.num == 1){
    add.pie(z=c(add/100, inbred/100, total/100), labels=NA, x=0.5*x.high, y=0.75*y.max, radius=(3/10)*(9/10)*y.range, 
            col=(c("gray50", "gray20", "white")), cex=1.3, label.dist=1.2, border="black")
  }
  else{
    excess.add <- round(add*(add.num-1), 1)
    excess.inbred <- round(inbred*(inbred.num-1), 1)
    add.pie(z=c(add/100, excess.add/100, inbred/100, excess.inbred/100, total/100), labels=NA, x=0.5*x.high, y=0.75*y.max, radius=(3/10)*(9/10)*y.range, 
            col=(c("gray50", "gray50", "gray20", "gray20", "white")), cex=1.3, label.dist=1.2, border="black")  
  }
  legend("topright", legend=c(add, inbred, total), fill=c("gray50", "gray20", "white"), title="% Variance", border="black", bty="n", text.col="black", cex=1.2)
  legend(x=-0.15*x.high/2, y=0.925*y.max, legend=c(round(add.effect, 2), round(inbred.effect, 2)), fill=c("gray50", "gray20"), title="Effect", border="black", bty="n", text.col="black", cex=1.2)
  with(dens, polygon(x=c(x[c(min.x, min.x:max.x, max.x)]), y=c(0, y[min.x:max.x], 0), col="gray", border="black"))
  lines(x=c(post.median, post.median), y=c(0, max.y), lty=5, lwd=2, col="white")
  if(cross.label1 != "" & cross.label2 != ""){
    text(1, median(dens$y), paste(cross.label1, 'x', cross.label2, sep=""))
  }
  axis(1, at=0:x.high, labels=T, tick=T, cex.axis=1.2)
  # Barplot
  #add.scatter(func=barPlotter(add.effect=add.effect, dom.effect=dom.effect, back.col=bgcolor), posi="topleft", inset=0, ratio=0.1)
  effectPlotter(add.effect=add.effect, inbred.effect=inbred.effect, x.range=x.high, y.range=y.max)
}

######################## Plot for all crosses
#' @export
diallelPlotter <- function(results, cross.type, pheno.name="", path=getwd(), add.num=1, inbred.num=1)
{
  pheno.list <- results[[1]]
  var.list <- results[[2]]
  ef.list <- results[[3]]
  qtl.num <- results[[4]]
  n <- results[[5]]
  # Placing labels
  label.indices <- 1:8
  if(pheno.name==""){
    pheno.name <- strsplit(deparse(substitute(pheno.list)), split=".", fixed=T)[[1]][1]
  }
  if(add.num == 1 & inbred.num == 1){
    pdf(paste(path,"/Diallel_", pheno.name, "_", cross.type, ".pdf",sep = ""),
        width=12, height = 12)
  }
  else{
    pdf(paste(path,"/Diallel_", pheno.name, "_", cross.type, "_", add.num, "add", inbred.num, "inbred.pdf",sep = ""),
        width=12, height = 12)    
  }
  par(mfrow=c(8,8), cex = 0.5, oma=c(1,1,1,0), mar=c(1,1,1,1))
  ## CHOOSE SUBSET OF DATA
  strains <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
  for(i in 1:length(strains)){
    for(j in 1:length(strains)){
      if(cross.type == "f2"){
        if(i == j){
          if(i == 1){
            infoPlotter(pheno.name, n=n, qtl.num[1], qtl.num[2])
          }
          else if(i == 2){
            legendPlotter(n=700, add.num=qtl.num[1], inbred.num=qtl.num[2])
          }
          else{
            emptyF2Plotter()
          }
        }
        else if(i > j){
          emptyPlotter()
        }      	
        else{
          oneParamPlotter(cross.u=pheno.list[[paste(strains[i], "x", strains[j], "-f2_eu", sep="")]],
                          add.effect=median(ef.list[[paste(strains[i], "x", strains[j], "-f2_add_ef", sep="")]]),
                          inbred.effect=median(ef.list[[paste(strains[i], "x", strains[j], "-f2_inbred_ef", sep="")]]),
                          add.perc=median(var.list[[paste(strains[i], "x", strains[j], "-f2_add_perc", sep="")]]),
                          inbred.perc=median(var.list[[paste(strains[i], "x", strains[j], "-f2_inbred_perc", sep="")]]),
                          add.num=add.num, inbred.num=inbred.num)
        }
      }	
      else if(cross.type == "bc"){
        if(i == j){
          if(i == 1){
            infoPlotter(pheno.name, n=n, qtl.num[1], qtl.num[2])
          }
          else if(i == 2){
            legendPlotter(n=700, add.num=qtl.num[1], inbred.num=qtl.num[2])
          }
          else{
            emptyBCPlotter()
          }
        }
        else{
          if(i < j){
            oneParamPlotter(pheno.list[[paste(strains[i], "x", strains[j], "-bc_m_eu", sep="")]],
                            add.effect=median(ef.list[[paste(strains[i], "x", strains[j], "-bc_m_add_ef", sep="")]]),
                            inbred.effect=median(ef.list[[paste(strains[i], "x", strains[j], "-bc_m_inbred_ef", sep="")]]),
                            add.perc=median(var.list[[paste(strains[i], "x", strains[j], "-bc_m_add_perc", sep="")]]),
                            inbred.perc=median(var.list[[paste(strains[i], "x", strains[j], "-bc_m_inbred_perc", sep="")]]),
                            add.num=add.num, inbred.num=inbred.num)          				
          }
          else{
            oneParamPlotter(pheno.list[[paste(strains[j], "x", strains[i], "-bc_p_eu", sep="")]],
                            add.effect=median(ef.list[[paste(strains[j], "x", strains[i], "-bc_p_add_ef", sep="")]]),
                            inbred.effect=median(ef.list[[paste(strains[j], "x", strains[i], "-bc_p_inbred_ef", sep="")]]),
                            add.perc=median(var.list[[paste(strains[j], "x", strains[i], "-bc_p_add_perc", sep="")]]),
                            inbred.perc=median(var.list[[paste(strains[j], "x", strains[i], "-bc_p_inbred_perc", sep="")]]),
                            add.num=add.num, inbred.num=inbred.num)   				
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

############################################## Extra unused plots
# Plotting a single variance explained
oneVarPlotter <- function(cross.var, cross.label1="", cross.label2="")
{
  palette(rev(gplots::rich.colors(1000)))
  post.mean <- mean(cross.var)
  post.median <- median(cross.var)
  bgcolor <- round(post.mean*10)
  
  plot(density(cross.var), xlim=c(0,100), lwd=3, xlab="", main="", frame.plot=F, ylab="", yaxt="n", xaxt="n", col="grey")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=bgcolor)
  lines(density(cross.var), lwd=3, col="grey")
  abline(v=post.median, lty=5, col="white")
  if(cross.label1 != "" & cross.label2 != ""){
    text(1, median(density(cross.var)$y), paste(cross.label1, 'x', cross.label2, sep=""))
  }
  axis(1, at=c(0, 25, 50, 75, 100), labels=T, tick=T, cex.axis=1.2)
}
######################## Plot of all variance explained
varPercPlotter <- function(var.list, cross.type, effect.type, pheno.name="", path=getwd())
{
  # Placing labels
  label.indices <- 1:8
  if(pheno.name==""){
    pheno.name <- strsplit(deparse(substitute(var.list)), split=".", fixed=T)[[1]][1]
  }
  pdf(paste(path,"/VarExp", pheno.name, "_", cross.type, effect.type, ".pdf",sep = ""),
      width=12, height = 12)
  par(mfrow=c(8,8), cex = 0.5, oma=c(1,1,1,0), mar=c(1,1,1,1))
  ## CHOOSE SUBSET OF DATA
  strains <- c("AJ", "B6", "129", "NOD", "NZO", "CAST", "PWK", "WSB")
  for(i in 1:length(strains))
  {
    for(j in 1:length(strains)){
      if(cross.type == "f2"){
        if(i == j){
          emptyF2Plotter()
        }
        else if(i > j){
          emptyPlotter()
        }      	
        else{
          if(effect.type == "add"){
            oneVarPlotter(var.list[[paste(strains[i], "x", strains[j], "-f2_add", sep="")]])
          }
          else if(effect.type == "dom"){
            oneVarPlotter(var.list[[paste(strains[i], "x", strains[j], "-f2_dom", sep="")]])
          }
          else{
            oneVarPlotter(var.list[[paste(strains[i], "x", strains[j], "-f2_add", sep="")]] + var.list[[paste(strains[i], "x", strains[j], "-f2_dom", sep="")]])
          }
        }
      }
      else if(cross.type == "bc"){
        if(i == j){
          emptyBCPlotter()
        }
        else{
          if(i < j){
            if(effect.type == "add"){
              oneVarPlotter(var.list[[paste(strains[i], "x", strains[j], "-bc_m_add", sep="")]])
            }
            else if(effect.type == "dom"){
              oneVarPlotter(var.list[[paste(strains[i], "x", strains[j], "-bc_m_dom", sep="")]])
            }
            else{
              oneVarPlotter(var.list[[paste(strains[i], "x", strains[j], "-bc_m_add", sep="")]] + var.list[[paste(strains[i], "x", strains[j], "-bc_m_dom", sep="")]])
            }           				
          }
          else{
            if(effect.type == "add"){
              oneVarPlotter(var.list[[paste(strains[j], "x", strains[i], "-bc_p_add", sep="")]])
            }
            else if(effect.type == "dom"){
              oneVarPlotter(var.list[[paste(strains[j], "x", strains[i], "-bc_p_dom", sep="")]])
            }
            else{
              oneVarPlotter(var.list[[paste(strains[j], "x", strains[i], "-bc_p_add", sep="")]] + var.list[[paste(strains[j], "x", strains[i], "-bc_p_dom", sep="")]])
            }  				
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

# Allows a different color than black for the text labels - currently unused
my.add.pie <- function (z, x = 0, y = 0, labels = names(z), radius = 1, edges = 200, 
                        clockwise = TRUE, init.angle = 90, density = NULL, angle = 45, 
                        col = NULL, border = NULL, lty = NULL, label.dist = 1.1, 
                        ...) 
{
  if (!is.numeric(z) || any(is.na(z) | z < 0)) 
    stop("'z' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(z))
  else labels <- as.graphicsAnnot(labels)
  z <- c(0, cumsum(z)/sum(z))
  dz <- diff(z)
  nz <- length(dz)
  asp <- get.asp()
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("#737373", "#F15A60", "#7BC36A", "#599BD3", "#F9A75B", 
        "#9E67AB", "#CE7058", "#D77FB4")
  else par("fg")
  col <- rep(col, length.out = nz)
  border <- rep(border, length.out = nz)
  lty <- rep(lty, length.out = nz)
  angle <- rep(angle, length.out = nz)
  density <- rep(density, length.out = nz)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = asp * radius * cos(t2p) + x, y = radius * sin(t2p) + 
           y)
  }
  for (i in 1L:nz) {
    n <- max(2, floor(edges * dz[i]))
    P <- t2xy(seq.int(z[i], z[i + 1], length.out = n))
    polygon(c(P$x, 0 + x), c(P$y, 0 + y), density = density[i], 
            angle = angle[i], border = border[i], col = col[i], 
            lty = lty[i])
    P <- t2xy(mean(z[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      text(label.dist * (P$x - x) + x, label.dist * (P$y - 
                                                       y) + y, labels[i], xpd = TRUE, adj = ifelse(P$x - 
                                                                                                     x < 0, 1, 0), col="gray", ...)
    }
  }
}

barPlotter <- function(add.effect, dom.effect, back.col)
{
  x.lb <- min(0, min(add.effect, dom.effect))
  x.ub <- max(0, max(add.effect, dom.effect))
  plot(NA, xlim=c(x.lb, x.ub), ylim=c(0,2), xaxt="n", yaxt="n", xlab="", ylab="", frame=F)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border=back.col, col=back.col, lwd=20)
  barplot(height=c(dom.effect, add.effect), width=c(1, 1), col=c("gray20", "gray50"), border="black", space=0, horiz=T, axes=F, add=T)
}
