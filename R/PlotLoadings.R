PlotLoadings = function( L_pj, whichfactor, Nspecies, addtitle=TRUE ){
  plot(1, type="n", xlim=c(1,nrow(L_pj)+0.5), ylim=range(L_pj)+diff(range(L_pj))*c(0,0.1), xlab="", ylab="", xaxt="n")
  if(addtitle==TRUE) mtext( text=paste("Factor",whichfactor), side=3, line=0.1, adj=0)
  abline(h=0)
  for(p in 1:Nspecies){
    lines(y=c(0,L_pj[p,whichfactor]), x=rep(p,2), lwd=5)
    text(x=(p+0.2), y=0+(0.02*max(L_pj[p,whichfactor])), label=rownames(L_pj)[p], srt=90, pos=4, cex=1.2)
  }
  legend( "top", legend=paste0("Proportion of explained variance= ",round(100*sum(L_pj[,whichfactor]^2)/sum(L_pj^2),1),"%"), bty="n")
}
