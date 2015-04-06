#--------------------------------------------------------------#
#     PRINCIPAL COMPONENT ANALYSIS(PCA) FUNCTION 09.25.2012    #
#     PCA(pc, graph)				                     #
#--------------------------------------------------------------#

PCAgraph <- function(pc, graph=c("scree", "biplot", "pcaplot")){

	if(graph == "scree"){
	plot(pc , main = "Scree Plot")
	}
	if(graph == "biplot"){
	biplot(pc, choice=1:2, cex=0.6, expand = 1, main = "Biplot")
	}
	if(graph == "pcaplot"){
	plot(pc$x, pch =20, main = "PCA Plot")
	}
}#--- end statement---#


