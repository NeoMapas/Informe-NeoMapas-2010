mi.ss <- NM.m1$Especieid!=1379
mi.s2 <- NM.m2$Especieid!=1379
for (transid in unique(NM.m1$IDTransecta)) {
  mi.ss <- NM.m1$Especieid!=1379 & NM.m1$IDTransecta %in% transid
  dtp <- with(NM.m1[mi.ss,],tapply(Ntotal,list(idpunto,Especieid),sum,na.rm=T))
  dtp[is.na(dtp)] <- 0
  write.csv(file=sprintf("AvesNM2010_%s.csv",transid),dtp)
}

mi.jmp <- jmp.NM[(jmp.NM$fase %in% 2) &jmp.NM$fam=="Pieridae" & !is.na(jmp.NM$genero) & jmp.NM$especie != "sp.",]
mi.jmp$esp <- paste(mi.jmp$genero,mi.jmp$especie)

for (transid in unique(mi.jmp$NM)) {
  mi.ss <- mi.jmp$NM %in% transid
  mi.tvn <-  tvn.NM[tvn.NM$fase %in% 2 & tvn.NM$NM %in% transid,]
  dtp <- with(mi.jmp[mi.ss,],tapply(jmp,list(factor(vst,levels=unique(mi.tvn$vst)),esp),function(x){length(unique(x))}))
  dtp[is.na(dtp)] <- 0
  rownames(dtp) <- paste("Sample",1:nrow(dtp),sep="")
  colnames(dtp) <- 1:ncol(dtp)
  write.csv(file=sprintf("PieridaeNM2006_NM%s.csv",transid),dtp)
}

##mi.jmp$esp[mi.jmp$esp %in% c("Eurema daira","Eurema elathea")] <- "Eurema sp. [daira o elathea]"

