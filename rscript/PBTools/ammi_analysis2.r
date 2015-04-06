ammi.analysis2 <- function (ENV, GEN, REP, Y, SIGMA2, number = TRUE, graphtype = "biplot", respVar) {
# ammi.analysis2 <- function (ENV, GEN, REP, Y, SIGMA2, number = TRUE, graphtype = "none", graphfor = "fixed", respVar) {
  
  #dir.create("plots")

  numrepEnv <- tapply(REP, ENV, mean)
  sigma2Env <- tapply(SIGMA2, ENV, mean)
  MSE <- sum(numrepEnv*sigma2Env)/sum(numrepEnv)
  NREP <- 1/mean(1/numrepEnv)

#   ammi2 <- ammi.analysis(ENV, GEN, REP, Y, MSE, number, graph = graphtype)
#   if (graphtype != "none") {
#     png(file = paste(getwd(),"/plots/",graphtype,"Mea2S_",respVar,".png",sep = ""))
#     ammiOut <- ammi.analysis(ENV, GEN, NREP, Y, MSE, graph = graphtype)
#     ammi2$biplot
#     dev.off()
#   } else {
#     ammi2 <- ammi.analysis(ENV, GEN, NREP, Y, MSE, graph = graphtype, gname = graphfor, yVar = respVar)
    ammi12Out <- ammi.analysis(ENV, GEN, NREP, Y, MSE, number, graph = graphtype, yVar = respVar)
#     ammi1 <- ammi1.plot(ammi12Out$biplot, Y, respVar) 
#   }
}
# bplot = ammiOut$biplot
