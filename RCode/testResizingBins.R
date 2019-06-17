
# cds.min<-32.5;
# cds.bin<-10;
# cds.n  <-5;
# cds<-cds.min+(0:cds.n)*cds.bin;
#
# cms.min<-25;
# cms.bin<-5;
# cms.n  <-10;
# cms<-cms.min+(0:cms.n)*cms.bin;

# cds.min<-25;
# cds.bin<-5;
# cds.n  <-10;
# cds<-cds.min+(0:cds.n)*cds.bin;
#
# cms.min<-25;
# cms.bin<-5;
# cms.n  <-12;
# cms<-cms.min+(0:cms.n)*cms.bin;

cds.min<-25;
cds.bin<-10;
cds.n  <-5;
cds<-cds.min+(0:cds.n)*cds.bin;

cms.min<-25;
cms.bin<-5;
cms.n  <-10;
cms<-cms.min+(0:cms.n)*cms.bin;

cat("cds =",cds,"\n");
cat("cms =",cms,"\n");
p_ij<-matrix(data=0,nrow=cds.n,ncol=cms.n);
for (j in 1:cms.n){
  cat("\n\n");
  cat("cms[",j,"] =",cms[j:(j+1)],"\n")
  for (i in 1:cds.n){
    cat("cds[",i,"] =",cds[i:(i+1)],"\n")
    if (cds[i+1]<=cms[j]) {
      p_ij[i,j]<-0;
    } else if (cms[j+1]<cds[i]){
      p_ij[i,j]<-0;
    } else {
      del<-cms[j+1]-cms[j];
      if ((cds[i]<cms[j])&(cms[j+1]<cds[i+1])) {
        cat("cds[",i,"]=",cds[i]," < cms[",j,"]=",cms[j],"\n",sep="");
        cat("cms[",j+1,"]=",cms[j+1]," <= cds[",i+1,"]=",cds[i+1],"\n",sep="");
        p_ij[i,j]<-1.0;
        cat("\t",p_ij[i,j],"\n")
        #break;
      } else if (cms[j]<=cds[i]) {
        cat("cms[",j,"]=",cms[j]," <= cds[",i,"]=",cds[i],"\n",sep="");
        delp<-cms[j+1]-max(cds[i],cms[j]);
        cat("\t",cms[j+1],max(cds[i],cms[j]),delp,"\n");
        p_ij[i,j]<-delp/del;
        cat("\t",p_ij[i,j],"\n")
        #break;
      } else if (cds[i+1]<=cms[j+1]) {
        cat("cds[",i+1,"]=",cds[i+1]," < cms[",j+1,"]=",cms[j+1],"\n",sep="");
        delp<-min(cds[i+1],cms[j+1])-cms[j];
        cat("\t",min(cds[i+1],cms[j+1]),cms[j],delp,"\n");
        p_ij[i,j]<-delp/del;
        cat("\t",p_ij[i,j],"\n")
        #break;
      }
    }

  }#--j
}#--i
View(p_ij)

sm<-colSums(p_ij);
names(sm)<-paste0(as.character(cms[1:(cms.n)]),"-",as.character(cms[2:(cms.n+1)]));
cat("cds = \n",cds,"\n")
print(sm);


