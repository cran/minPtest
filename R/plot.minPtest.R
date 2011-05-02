plot.minPtest <-
function(x, level=0.05, lambda=1, gene.name=FALSE, ...){
  psnp <- x$psnp
  minp <- x$minp
  if(length(which(minp[,"minP"]==0))>0){
    minp[,"minP"][which(minp[,"minP"]==0)] <- 0.0000000001
  }
  SNPtoGene <- x$SNPtoGene
  p.adj.snp <- x$p.adj.psnp
  p.adj.minp <- x$p.adj.minp
  psnp.padj <- cbind(psnp,p.adj.snp)
  genes <- rownames(minp)
  gensnp <- lapply(seq_along(genes), function(i){
    snpnames <- SNPtoGene[which(SNPtoGene[,2]==genes[i]),1]
    psnp.padj.temp <- psnp.padj[snpnames,]
    if(length(snpnames)==1){
      psnp.padj.temp <- matrix(psnp.padj.temp,nrow=1,ncol=2)
      rownames(psnp.padj.temp) <- snpnames
      colnames(psnp.padj.temp) <- c("p_value","p.adjust")
    }
    psnp.padj.temp
  })
  names(gensnp) <- genes
  gensnp.sort <- do.call(rbind,gensnp)
  log.psnp <- -log(gensnp.sort[,"p_value"],10)
  par(mar=c(5,4,4,5))
  plot(x=c(0,x$nrsnp),y=c(0,max(log.psnp)), xlab="SNP", ylab=expression(-log[10](psnp)),xaxt="n",col="white", ...)
  for(i in 1:x$nrsnp){
    if(gensnp.sort[i,"p.adjust"]<=level){
      points(i,log.psnp[i], pch=20, col="red", ...)
    }else{
      points(i,log.psnp[i], pch=20, ...)
    }
  }
  minp.padj <- cbind(minp,p.adj.minp)
  log.minp <- -log(minp[,"minP"],10)
  gene.size <- unlist(lapply(gensnp, function(i){
    l <- nrow(i)
  }))
  x1 <- 1
  x2 <- gene.size[1]
  middle <- median(x1:x2)
  for(i in 1:(length(gene.size))){
    par(new=TRUE)
    plot(x=c(0,x$nrsnp),y=c(0,max(log.minp)), xlab="", ylab="",xaxt="n",col="white", yaxt="n", ...)
    if(minp.padj[i,"p.adjust"]<=level){
      lines(c(x1,x2),c((log.minp[i])*lambda,(log.minp[i])*lambda),col="red", ...)
    }else{
      lines(c(x1,x2),c((log.minp[i])*lambda,(log.minp[i])*lambda), ...)
    }
    x1 <- x2+1
    x2 <- x2+gene.size[i+1]
    if(i!=length(gene.size)){
      middle <- c(middle,median(x1:x2))
    }
  }
  axis(4, at=seq(0,max(log.minp)), labels=seq(0,max(log.minp)), ...)
  if(lambda==1){
    mtext(text=expression(-log[10](minp)), side=4, line=3, ...)
  }else{
    mtext(text=bquote(.(-lambda) * (log[10](minp))), side=4, line=3, ...)
  }
  if(gene.name){
    axis(1, at=middle,labels=genes, ...)
  }
box(...)
}


