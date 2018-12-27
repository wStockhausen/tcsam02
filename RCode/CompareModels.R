
load("../AllModelRunResults.RData")
models<-fullmodels; rm(fullmodels,snglmodels);
compare<-names(models) %in% c("18A","18B","18C0","18C0a","18C1","18C1a","18C2a","18C3a","18D0");

dfr.PVs   <-rCompTCMs::extractMDFR.Results.ParameterValues(models[compare]);
write.csv(dfr.PVs,"ParamValues.csv",row.names=FALSE)
dfr.PsAtBs<-rCompTCMs::extractMDFR.Results.ParametersAtBounds(models[compare]);
write.csv(dfr.PsAtBs,"ParamsAtBounds.csv",row.names=FALSE)

dfr.OFCs<-rTCSAM02::getMDFR.OFCs.DataComponents(models[compare]);
write.csv(dfr.OFCs,"OFCs.DataComponents.csv",row.names=FALSE)
dfr.OFCs<-rTCSAM02::getMDFR.OFCs.NonDataComponents(models[compare]);
write.csv(dfr.OFCs,"OFCs.NonDataComponents.csv",row.names=FALSE)

dfr.OFLs<-rTCSAM02::getMDFR.OFLResults(models[compare],verbose=TRUE);
dfr.OFLs$avgRec<-dfr.OFLs$avgRecM + dfr.OFLs$avgRecF;
dfr.OFLs<-wtsUtilities::deleteCols(dfr.OFLs,cols=c("avgRecM","avgRecF"),debug=TRUE);
write.csv(dfr.OFLs,"oflResults.csv",row.names=FALSE)

# rCompTCMs::compareResults.Pop.NaturalMortality(models[compare],dodge=0)
# rCompTCMs::compareResults.Pop.MeanGrowth(models[compare],dodge=0)
# rCompTCMs::compareResults.Pop.GrowthMatrices.LinePlots(models[compare],dodge=0)
# rCompTCMs::compareResults.Pop.PrM2M(models[compare],dodge=0)
# rCompTCMs::compareFits.MaturityData(models[compare],dodge=0,types="fits",plot1stObs=TRUE)
# rCompTCMs::compareResults.Pop.Recruitment(models[compare],dodge=0,numRecent=30)
# rCompTCMs::compareResults.Pop.MatureBiomass(models[compare],dodge=0,numRecent=30)
#
# rCompTCMs::compareFits.EffectiveNs(models[compare],
#                                    fleet.type="survey",
#                                    facet_grid="x+m+s~.",
#                                    dodge=0)
#
# rCompTCMs::compareFits.EffectiveNs(models[compare],
#                                    fleet.type="fishery",
#                                    category="retained",
#                                    facet_grid="x+m+s~.",
#                                    dodge=0)
#
# rCompTCMs::compareFits.EffectiveNs(models[compare],
#                                    fleet.type="fishery",
#                                    category="total",
#                                    facet_grid="x+m+s~.",
#                                    dodge=0)
#
# rCompTCMs::compareResults.Surveys.Catchability(models[compare],dodge=0)
# rCompTCMs::compareResults.Surveys.SelFcns(models[compare],dodge=0,years=c(1981,2018))
# rCompTCMs::compareResults.Surveys.CaptureProbs(models[compare],dodge=0,years=c(1981,2018))
# rCompTCMs::compareFits.BiomassData(models[compare],fleet.type="survey",catch.type="index",
#                                    plot1stObs=TRUE,numRecent=30)
#
#create pdf output
rCompTCMs::modelComparisons.ModelFits.OtherData(models[compare],
                                               plot1stObs=TRUE,
                                               output_format="pdf_document",
                                               clean=TRUE);
rCompTCMs::modelComparisons.ModelFits.ACD(models[compare],
                                         plot1stObs=TRUE,
                                         output_format="pdf_document",
                                         clean=TRUE);
rCompTCMs::modelComparisons.ModelFits.ZCsByYear(models[compare],
                                         plot1stObs=TRUE,
                                         output_format="pdf_document",
                                         clean=TRUE);
rCompTCMs::modelComparisons.ModelFits.ZCs(models[compare],
                                          type="Surveys",
                                          plot1stObs=TRUE,
                                          output_format="pdf_document",
                                          clean=TRUE);
rCompTCMs::modelComparisons.ModelFits.ZCs(models[compare],
                                          type="Fisheries",
                                          plot1stObs=TRUE,
                                          output_format="pdf_document",
                                          clean=TRUE);
rCompTCMs::modelComparisons.PopProcesses(models[compare],
                                         output_format="pdf_document",
                                         clean=TRUE);
rCompTCMs::modelComparisons.PopQuantities(models[compare],
                                         output_format="pdf_document",
                                         clean=TRUE);
rCompTCMs::modelComparisons.Characteristics.Surveys(models[compare],
                                                    output_format="pdf_document",
                                                    clean=TRUE);
rCompTCMs::modelComparisons.Characteristics.Fisheries.(models[compare],
                                                       output_format="pdf_document",
                                                       clean=TRUE);
# rCompTCMs::modelComparisons.ParameterTables(models[compare],
#                                            output_format="pdf_document");

# rCompTCMs::modelComparisons(models[compare],
#                             plot1stObs=FALSE,
#                             output_format="pdf_document");

