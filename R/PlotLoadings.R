#' @export
PlotLoadings = function( L_pj, whichfactor, Nspecies, addtitle=TRUE, LabelPosition="Right", Buffer=c(0,0.1), Labels=rownames(L_pj), Cex=1.2 ){
  plot(1, type="n", xlim=c(0.5,nrow(L_pj)+0.5), ylim=range(L_pj)+diff(range(L_pj))*Buffer, xlab="", ylab="", xaxt="n", xaxs="i" )
  if(LabelPosition=="Xaxis") axis( side=1, at=1:nrow(L_pj), labels=TRUE, las=3)
  if(addtitle==TRUE) mtext( text=paste("Factor",whichfactor), side=3, line=0.1, adj=0)
  abline(h=0)
  for(p in 1:Nspecies){
    lines(y=c(0,L_pj[p,whichfactor]), x=rep(p,2), lwd=5)
    if(LabelPosition=="Right") text(x=(p+0.2), y=0+(0.02*max(L_pj[p,whichfactor])), labels=Labels[p], srt=90, pos=4, cex=Cex)
    if(LabelPosition=="Above") text(x=p, y=max(L_pj), labels=Labels[p], srt=90, pos=3, cex=Cex)
    if(LabelPosition=="Below") text(x=p, y=min(L_pj), labels=Labels[p], srt=90, pos=1, cex=Cex)
  }
  legend( "top", legend=paste0("Proportion of explained variance= ",round(100*sum(L_pj[,whichfactor]^2)/sum(L_pj^2),1),"%"), bty="n")
}
