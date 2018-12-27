#convert TCSAM02 model output to an RData file
require(rTCSAM02)
best<-rTCSAM02::getResLst('./best',verbose=TRUE)
if (!is.null(best)) {
  save(best,file="Results.Model.RData");
  mdfr<-rTCSAM02::getMDFR.OFCs.DataComponents(best,verbose=FALSE)
  write.csv(mdfr,file="OFCs.DataComponents.csv",row.names=FALSE)
  mdfrp<-rTCSAM02::getMDFR.OFCs.NonDataComponents(best,verbose=FALSE)
  write.csv(mdfrp,file="OFCs.NonDataComponents.csv",row.names=FALSE)
}
