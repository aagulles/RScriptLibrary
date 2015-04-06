ANOVA <-
function (Y, data, cont=1, DESIGN, TRT, Block, Env, Envname)  {
	
	if (DESIGN == "CRD")    {model <- aov(Y ~ Env+TRT+Env:TRT, data=data) }
	if (DESIGN == "RCB")    {model <- aov(Y ~ Env+Block:Env+TRT+Env:TRT, data=data) }
	
	out <- anova(model)
	SSenv <- out[rownames(out)==Envname, "Sum Sq"]/cont
	SSg <- out[rownames(out)=="TRT", "Sum Sq"]/cont
	envtrt<-paste(Envname,":TRT", sep="")
	SSgxenv <- out[rownames(out)==envtrt, "Sum Sq"]/cont
	SSe <- out[rownames(out)=="Residuals", "Sum Sq"]/cont
	
	DFenv <- out[rownames(out)==envtrt, "Df"]
	DFg <- out[rownames(out)=="TRT", "Df"]
	DFgxenv <- out[rownames(out)==envtrt, "Df"]
	DFe <- out[rownames(out)=="Residuals", "Df"]
	
	MSenv <- SSenv/DFenv
	MSg <- SSg/DFg
	MSgxenv <- SSgxenv/DFgxenv
	MSe <- SSe/DFe
	
	return (list(SSenv=SSenv, SSg=SSg, SSgxenv=SSgxenv, SSe=SSe, 
					DFenv=DFenv, DFg=DFg, DFgxenv=DFgxenv, DFe=DFe, 
					MSenv=MSenv, MSg=MSg, MSgxenv=MSgxenv, MSe=MSe))
}

