snglmodels<-list();

load("../AssessmentModelRuns/SingleRuns/0a.17AM.B2b/Results.Model.RData");     snglmodels$`17AM`  <-best; rm(best);
load("../AssessmentModelRuns/18AM17/Results.Model.RData");                     snglmodels$`18AM17`<-best; rm(best);
load("../AssessmentModelRuns/SingleRuns/0b.17AMu/Results.Model.RData");        snglmodels$`17AMu` <-best; rm(best);
load("../AssessmentModelRuns/SingleRuns/1.2018B0.18A/Results.Model.RData");    snglmodels$`18A`   <-best; rm(best);
load("../AssessmentModelRuns/SingleRuns/2.2018B1.18B/Results.Model.RData");    snglmodels$`18B`   <-best; rm(best);

fullmodels<-list();
load("../AssessmentModelRuns/FullRuns/1.2018B0.18A.Full/Results.Model.RData");    fullmodels$`18A`   <-best; rm(best);
load("../AssessmentModelRuns/FullRuns/2.2018B1.18B.Full/Results.Model.RData");    fullmodels$`18B`   <-best; rm(best);
load("../AssessmentModelRuns/FullRuns/3.2018B2.18C0.Full/Results.Model.RData");   fullmodels$`18C0`  <-best; rm(best);
load("../AssessmentModelRuns/FullRuns/3.2018B2.18C0a.Full/Results.Model.RData");  fullmodels$`18C0a` <-best; rm(best);
load("../AssessmentModelRuns/FullRuns/5.2018B4.18C1.Full/Results.Model.RData");   fullmodels$`18C1`  <-best; rm(best);
load("../AssessmentModelRuns/FullRuns/5.2018B4.18C1a.Full/Results.Model.RData");  fullmodels$`18C1a` <-best; rm(best);
load("../AssessmentModelRuns/FullRuns/7.2018B6.18C2a.Full/Results.Model.RData");  fullmodels$`18C2a` <-best; rm(best);
load("../AssessmentModelRuns/FullRuns/8.2018B7.18C3a.Full/Results.Model.RData");  fullmodels$`18C3a` <-best; rm(best);
load("../AssessmentModelRuns/FullRuns/4.2018B3.18D0.Full/Results.Model.RData");   fullmodels$`18D0`  <-best; rm(best);
#load("../AssessmentModelRuns/FullRuns/6.2018B5.18D1.Full/Results.Model.RData");   fullmodels$`18D1`  <-best; rm(best);

save(snglmodels,fullmodels,file="./AllModelRunResults.RData")

