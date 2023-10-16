###################################################
### chunk number 1: 
###################################################
#line 58 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
load(file="~/NeoMapas/Rdata/Anexos.rda")


###################################################
### chunk number 2: 
###################################################
#line 62 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
nms <- unique(c(scrb.rsm$NM,trans.info$NM,unique(tvn.NM$NM),unique(avs.NM$NM)))
nms <- as.numeric(nms[!is.na(nms)])
mi.info <- info.NM[info.NM$NM %in% nms,c("NM","Nombre","ADM1","Bioregión","CNEB")]
mi.info <- merge(mi.info, CNEB.nm[,c("cdg","x","y","PC1","PC2","PC3")], by.x="CNEB", by.y="cdg")

mi.info$NM <- sprintf("%02s",mi.info$NM)

mi.info <- mi.info[order(mi.info$Bioregión,mi.info$NM),]

mi.info$lon <- sprintf("%02iº%02i' W",floor(abs(mi.info$x)),round((1-(mi.info$x %% 1))*60))
mi.info$lat <- sprintf("%02iº%02i' N",floor(abs(mi.info$y)),round(((mi.info$y %% 1))*60))

#Bird method calibration (2001-2003)
mi.info$act <- ""
slc <- unique(avs.NM$NM)

mi.info[mi.info$NM %in% slc,"act"] <- paste(mi.info[mi.info$NM %in% slc,"act"],"AV01", sep=", ")

#Butterfly method calibration (2003-2005)
slc <- tvn.NM$NM[tvn.NM$fase %in% 1]
mi.info[mi.info$NM %in% slc,"act"] <- paste(mi.info[mi.info$NM %in% slc,"act"],"PP03", sep=", ")

##Dung Beetles method calibration (2005-2008)
slc <- scrb.rsm$NM[scrb.rsm$yr %in% c(2005,2008)]
mi.info[mi.info$NM %in% slc,"act"] <- paste(mi.info[mi.info$NM %in% slc,"act"],"SC05", sep=", ")

##First National Butterfly Survey (2006)
slc <- tvn.NM$NM[tvn.NM$fase %in% 2]
mi.info[mi.info$NM %in% slc,"act"] <- paste(mi.info[mi.info$NM %in% slc,"act"],"PP06", sep=", ")

##Dung Beetles pilot study (2006) 
slc <- scrb.rsm$NM[scrb.rsm$yr %in% c(2006)]
mi.info[mi.info$NM %in% slc,"act"] <- paste(mi.info[mi.info$NM %in% slc,"act"],"SC06", sep=", ")

##Second National Butterfly Survey (2009-2010)
slc <- tvn.NM$NM[tvn.NM$fase %in% 3]
mi.info[mi.info$NM %in% slc,"act"] <- paste(mi.info[mi.info$NM %in% slc,"act"],"PP09", sep=", ")

##First National Dung Beetle Survey (2009-2010) 
slc <- scrb.rsm$NM[scrb.rsm$yr %in% c(2009,2010)]
mi.info[mi.info$NM %in% slc,"act"] <- paste(mi.info[mi.info$NM %in% slc,"act"],"SC09", sep=", ")

##First National Bird Survey (2010)
slc <- trans.info$NM
mi.info[mi.info$NM %in% slc,"act"] <- paste(mi.info[mi.info$NM %in% slc,"act"],"AV10", sep=", ")
mi.info$act <- trim(sub("^, ","",mi.info$act))


regs <- cumsum(table(mi.info$Bioregión)[-1])


###################################################
### chunk number 3: tabla general
###################################################
#line 118 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
print(xtable(mi.info[,c("NM","Nombre","ADM1","CNEB","PC1","PC2","PC3","act")], caption=paste("Transects included in NeoMaps' activities from 2001 to 2010, grouped by bioregions. NM: NeoMaps transect code and name, State: Venezuelan major administrative units, VBG: code of rows (letter) and column (number) of the corresponding cell in the Venezuelan Biodiversity Grid (VBG, see Fig.~1 in main document), PC1-PC3: scores of the corresponding VBG cell in the principal components analysis used for sampling design (see Methods in main document). Activities: See text for explanations."), label="TAB:RSMINFO",  align="c|c|lp{2.6cm}c|ccc|p{3cm}|", digits=c(0,0,0,0,0,3,3,3,0)), caption.placement="top", include.rownames=F, include.colnames=F, add.to.row=list(pos=as.list(c(-1,-1,0,regs[-5])), command=c("\\hline\n NM & Name & State & VBG & PC1 & PC2 & PC3 & Activities\\\\", "\\endfirsthead 
\\hline\n NM & Name & State & VBG & PC1 & PC2 & PC3 & Activities\\\\
\\hline \\multicolumn{8}{l}{Continued from previous page}\\\\\\hline\\endhead\n \\endlastfoot\n\\hline\\multicolumn{8}{l}{Continues on next page}\\\\ \\hline\\endfoot\n","\\hline \\multicolumn{8}{l}{Occident}\\\\ \\hline","\\hline \\multicolumn{8}{l}{Andean mountains}\\\\","\\hline \\multicolumn{8}{l}{Coastal mountains}\\\\","\\hline \\multicolumn{8}{l}{Orinoco floodplain}\\\\","\\hline \\multicolumn{8}{l}{Guayana shield}\\\\")), hline.after=regs[-5], floating=FALSE, tabular.environment="longtable", table.placement="H", size="small")


###################################################
### chunk number 4: 
###################################################
#line 128 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
sd.mrps <- apply(log1p(mrps.NM),1,sd)
abnd.mrps <- apply(log1p(mrps.NM),1,mean)
H0.mrps <- diversity(mrps.NM, index = "shannon")
H1.mrps <- renyi(mrps.NM,scales=1)
ACE.mrps <- estimateR(mrps.NM)["S.ACE",]
## realmente, según el código de estimateR es standard deviation
se.mrps <- estimateR(mrps.NM)["se.ACE",]


tt <- aggregate(tvn.NM$sfrz,list(NM=tvn.NM$NM,fase=tvn.NM$fase),sum,na.rm=T)
tt$x <- tt$x/60

t2 <- aggregate(jmp.NM[,c("jmp","spp")],list(NM=jmp.NM$NM,fase=jmp.NM$fase),luq)

mrps.info <- merge(tt,t2,by=c("NM","fase"),all.x=T)

mrps.tab <- cbind(rbind(mrps.info[mrps.info$fase==1,],matrix(NA,ncol=5, nrow=sum(mrps.info$fase==3)-sum(mrps.info$fase==1), dimnames=list(1:(sum(mrps.info$fase==3)-sum(mrps.info$fase==1)),colnames(mrps.info)))),
rbind(mrps.info[mrps.info$fase==2,],matrix(NA,ncol=5, nrow=sum(mrps.info$fase==3)-sum(mrps.info$fase==2), dimnames=list(1:(sum(mrps.info$fase==3)-sum(mrps.info$fase==2)),colnames(mrps.info)))),
mrps.info[mrps.info$fase==3,])


###################################################
### chunk number 5: 
###################################################
#line 152 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
t3 <- aggregate(jmp.NM[jmp.NM$fam %in% "Pieridae" & jmp.NM$fase %in% 2, c("jmp","spp")],list(NM=jmp.NM$NM[jmp.NM$fam %in% "Pieridae" & jmp.NM$fase %in% 2],fase=jmp.NM$fase[jmp.NM$fam %in% "Pieridae" & jmp.NM$fase %in% 2]),luq)

mrps.rsm <- merge(tt,t3,by=c("NM","fase"),all.x=F)

colnames(mrps.rsm) <- c("NM","fase","sfrz","jmp","spp")



mrps.tmp <- data.frame(mrps.rsm[,c("NM","sfrz","jmp","spp")])
mrps.rslt <- data.frame(log.abnd=abnd.mrps, sd.log.abnd=sd.mrps, H0=H0.mrps, ACE=ACE.mrps, se.ACE=se.mrps, NM=sub("NM","",sapply(rownames(mrps.NM),function(x) strsplit(x," :: ")[[1]][1], simplify=T)))

mrps.rslt <- merge(mrps.tmp, mrps.rslt, by.x=c("NM"), by.y=c("NM"), all.x=T)

mrps.rslt <- mrps.rslt[order(mrps.rslt$NM),]



###################################################
### chunk number 6: 
###################################################
#line 178 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
print(xtable(mrps.tab[,c(1,3:6,8:11,13:15)],
caption=paste("Transects included in NeoMaps butterfly surveys from 2003 to 2010, $E$: sampling effort in hours*collector (h*p), $N_{total}$: total number of individuals of sampled, $S_{obs}$: number of species of six families of butterflies identified so far.\\vspace{1cm}\\mbox{}"), label="TAB:MRPSINFO", align="c|cccc|cccc|cccc|", digits=c(0,0,1,0,0,0,1,0,0,0,1,0,0)), caption.placement="top", include.rownames=F, include.colnames=F, add.to.row=list(pos=list(-1), command=c("\n\n\\multicolumn{4}{c}{Method calibration 2003-2005} & \\multicolumn{4}{c}{First survey 2006}& \\multicolumn{4}{c}{Second survey 2009-2010}\\\\\n
\\hline\n
 NM & $E$ & $N_{total}$ & $S_{obs}$ & 
 NM & $E$ & $N_{total}$ & $S_{obs}$ & 
 NM & $E$ & $N_{total}$ & $S_{obs}$ \\\\\n")), hline.after=c(0,nrow(mrps.tab)),  table.placement="H", size="normalsize")


###################################################
### chunk number 7: 
###################################################
#line 187 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
print(xtable(mrps.rslt, caption="Summary results for butterfly transects from the first national survey (2006) for the family Pieridae. NM: NeoMaps transect code (see table \\ref{TAB:RSMINFO} for details), $E$: sampling effort in hours*collector (h*p), $N_{total}$: total number of individuals of Pieridae sampled, $S_{obs}$: number of species detected, Estimated mean and standard deviation of dung beetle species abundance in log-scale  ($\\log N (\\mu \\pm \\sigma)$), Shannon-Weaver diversity index, ($H$), and abundance coverage estimator of dung beetle species richness (ACE) for NeoMaps transects.\\vspace{1cm}", label="TAB:MRPSEST", align="|l|cccc|rlcrl|", digits=c(0,0,1,0,0,3,3,3,2,2)), caption.placement="top", include.rownames=F, include.colnames=F, add.to.row=list(pos=list(-1), command=c("\\hline\n NM & $E$ (h*t) & $N_{total}$ & $S_{obs}$& \\multicolumn{2}{|c}{$\\log N (\\mu \\pm \\sigma)$} & $H$ & \\multicolumn{2}{c|}{ACE $\\pm$ s.e.} \\\\")), table.placement="H", hline.after=c(-1,nrow(mrps.rslt)), size="normalsize")


###################################################
### chunk number 8: MTRZmrps
###################################################
#line 198 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"

mx.scrb <- as.matrix(d1.mrps)
rownames(mx.scrb) <- substr(rownames(mx.scrb),3,4)
colnames(mx.scrb) <- substr(colnames(mx.scrb),3,4)
mx.scrb <- mx.scrb[order(info.NM[match(as.numeric(rownames(mx.scrb))
,info.NM$NM), "Bioregión"]), order(info.NM[match(as.numeric(colnames(mx.scrb))
,info.NM$NM), "Bioregión"])]

mx.scrb[lower.tri(mx.scrb)] <- NA

mis.cols <- brewer.pal(11,"BrBG")
mis.cols <- grey(seq(14,25,length=11)/30)

par(mar=c(3,3,3,0),oma=c(0,0,0,0),xpd=T)
image(1:nrow(mx.scrb),1:ncol(mx.scrb),mx.scrb, axes=F, col=c(NA,rev(mis.cols)), xlab="", ylab="")

polygon(cumsum(table(info.NM[match(as.numeric(rownames(mx.scrb)),info.NM$NM),"Bioregión"]))[c(1,2,2,3,3,4,4,5,5,6,6,5,5,4,4,3,3,2,2,1,1)]+.5,cumsum(table(info.NM[match(as.numeric(rownames(mx.scrb)),info.NM$NM),"Bioregión"]))[c(1,1,2,2,3,3,4,4,5,5,6,6,5,5,4,4,3,3,2,2,1)]+.5,lwd=2)

atta <- aggregate(1:ncol(mx.scrb),by=list(info.NM[match(as.numeric(rownames(mx.scrb)),info.NM$NM),"Bioregión"]),mean)$x
text(atta+c(0,0.5,0,0,0),c(-.5,6.5,9.5,13.5,21.5),c("Occident","\nAndean\nmountains","Coastal\nmountains","Orinoco\nfloodplain","Guayana\nshield"),cex=1)

axis(3,las=2,1:ncol(mx.scrb),rownames(mx.scrb),cex.axis=1)
axis(2,las=2,1:ncol(mx.scrb),rownames(mx.scrb),cex.axis=1)

text(rep(1:nrow(mx.scrb),ncol(mx.scrb)), rep(1:nrow(mx.scrb),rep(ncol(mx.scrb),ncol(mx.scrb))) ,round(as.vector(mx.scrb),2),cex=.5)



###################################################
### chunk number 9: 
###################################################
#line 237 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"

spps.NM <- colnames(mrps.NM)
genero <- unname(sapply(spps.NM,function(x) {strsplit(x," ")[[1]][1]}))
epiteto <- unname(sapply(spps.NM,function(x) {paste(strsplit(x," ")[[1]][-1],collapse=" ")}))

epiteto[grep("sp.",spps.NM)] <- "elathea"
epiteto[grep("nise",spps.NM)] <- "venusta"
epiteto[grep("daira",spps.NM)] <- "daira"

NM.pird <- Mariposas[match(paste(genero,epiteto),paste(Mariposas$gen,Mariposas$esp)),]

NM.pird <- NM.pird[order(NM.pird$ecdg),]

NM.pird$aut <- sub("�f","äf",sub("�b","üb",NM.pird$aut))



###################################################
### chunk number 10: 
###################################################
#line 255 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
cat(paste("\\begin{itemize}"))

for (j in unique(NM.pird$sfam)) {
	cat(sprintf("\\item Subfamily %s",j))
	strb <- unique(NM.pird$trb[NM.pird$sfam %in% j])
	cat(sprintf("\\begin{itemize}"))
	for (k in strb) {
		if (!is.na(k) & k!="NA") {
			cat(sprintf("\\item Tribe %s",k))
		}	
		cat(paste("\\begin{enumerate}"))
		for (i in seq(along=NM.pird$trb)[NM.pird$sfam %in% j & NM.pird$trb %in% k]) {
			cat(paste("\n\\item", sub("<i>", "\\\\textit{",sub("</i>","}", sub("<A>",sub("&","\\&",NM.pird$aut[i]), sub("<S>", "", sub("<F>",NM.pird$fch[i], sub("<E>",NM.pird$esp[i], sub("<G>",NM.pird$gen[i],NM.pird$stx[i])))))))))
		}
		cat(paste("\n\\end{enumerate}\n"))
	}
	cat(sprintf("\\end{itemize}"))
}
	cat(sprintf("\\end{itemize}"))



###################################################
### chunk number 11: 
###################################################
#line 282 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
sd.scrb <- apply(log1p(scrb.NM),1,sd)
abnd.scrb <- apply(log1p(scrb.NM),1,mean)
H0.scrb <- diversity(scrb.NM, index = "shannon")
H1.scrb <- renyi(scrb.NM,scales=1)
ACE.scrb <- estimateR(scrb.NM)["S.ACE",]
## realmente, según el código de estimateR es standard deviation
se.scrb <- estimateR(scrb.NM)["se.ACE",]

 scrb.rsm$jmp[is.na(scrb.rsm$jmp)] <- scrb.rsm$njmp[is.na(scrb.rsm$jmp)]


scrb.info <- data.frame(scrb.rsm[,c("NM","yr","sfrz","jmp","spp")])
scrb.rslt <- data.frame(log.abnd=abnd.scrb, sd.log.abnd=sd.scrb, H0=H0.scrb, ACE=ACE.scrb, se.ACE=se.scrb, t(sapply(rownames(scrb.NM),function(x) strsplit(x,"::")[[1]], simplify=T)))

scrb.rslt <- merge(scrb.info, scrb.rslt, by.x=c("NM","yr"), by.y=c("X1","X2"), all.x=T)

NMNs <- paste(scrb.rslt$NM[is.na(scrb.rslt$spp)],collapse=", ")
regs <- cumsum(table(scrb.rslt$bioreg))[-1]
NM2s <- paste(scrb.rslt$NM[duplicated(scrb.rslt$NM)],collapse=", ")


scrb.rslt <- scrb.rslt[order(scrb.rslt$NM),]



###################################################
### chunk number 12: 
###################################################
#line 316 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"

 print(xtable(scrb.rslt, caption="Transects included in NeoMaps' Dung beetle surveys from 2005 to 2010. NM: NeoMaps transect code (see table \\ref{TAB:RSMINFO} for details), Year of sampling, $E$: sampling effort in hours*trap (h*t), $N_{total}$: total number of individuals of Scarabaeinae sampled, $S_{obs}$: number of species detected. Summary results are presented for a subset of transects: estimated mean and standard deviation of dung beetle species abundance in log-scale  ($\\log N (\\mu \\pm \\sigma)$), Shannon-Weaver diversity index, ($H$), and abundance coverage estimator of dung beetle species richness (ACE) for NeoMaps transects.", label="TAB:SCRBEST", align="|l|ccccc|rlcrl|", digits=c(0,0,0,1,0,0,3,3,3,2,2)), caption.placement="top", include.rownames=F, include.colnames=F, add.to.row=list(pos=list(-1), command=c("\\hline\n NM & Year & $E$ (h*t) & $N_{total}$ & $S_{obs}$& \\multicolumn{2}{|c}{$\\log N (\\mu \\pm \\sigma)$} & $H$ & \\multicolumn{2}{c|}{ACE $\\pm$ s.e.} \\\\")), table.placement="H", hline.after=c(-1,nrow(scrb.rslt)), size="small")


###################################################
### chunk number 13: MTRZscrb
###################################################
#line 328 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"

mx.scrb <- as.matrix(d1.scrb)
rownames(mx.scrb) <- substr(rownames(mx.scrb),0,2)
colnames(mx.scrb) <- substr(colnames(mx.scrb),0,2)
mx.scrb <- mx.scrb[order(info.NM[match(as.numeric(rownames(mx.scrb))
,info.NM$NM), "Bioregión"]), order(info.NM[match(as.numeric(colnames(mx.scrb))
,info.NM$NM), "Bioregión"])]

mx.scrb[lower.tri(mx.scrb)] <- NA


par(mar=c(3,3,3,0),oma=c(0,0,0,0),xpd=T)
image(1:nrow(mx.scrb),1:ncol(mx.scrb),mx.scrb, axes=F, col=c(NA,rev(mis.cols)), xlab="", ylab="")

polygon(cumsum(table(info.NM[match(as.numeric(rownames(mx.scrb)),info.NM$NM),"Bioregión"]))[c(1,2,2,3,3,4,4,5,5,6,6,5,5,4,4,3,3,2,2,1,1)]+.5,cumsum(table(info.NM[match(as.numeric(rownames(mx.scrb)),info.NM$NM),"Bioregión"]))[c(1,1,2,2,3,3,4,4,5,5,6,6,5,5,4,4,3,3,2,2,1)]+.5,lwd=2)

atta <- aggregate(1:ncol(mx.scrb),by=list(info.NM[match(as.numeric(rownames(mx.scrb)),info.NM$NM),"Bioregión"]),mean)$x
text(atta,c(-.5,3.5,7.5,11.5,15.5),c("Occident","\nAndean\nmountains","Coastal\nmountains","Orinoco\nfloodplain","Guayana\nshield"),cex=1)

axis(3,las=2,1:ncol(mx.scrb),rownames(mx.scrb),cex.axis=1)
axis(2,las=2,1:ncol(mx.scrb),rownames(mx.scrb),cex.axis=1)

text(rep(1:nrow(mx.scrb),ncol(mx.scrb)), rep(1:nrow(mx.scrb),rep(ncol(mx.scrb),ncol(mx.scrb))) ,round(as.vector(mx.scrb),2),cex=.5)



###################################################
### chunk number 14: 
###################################################
#line 366 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
colnames(scrb.NM)[grep("bordoni",colnames(scrb.NM))] <- "Genieridium.bordoni"
colnames(scrb.NM)[grep("atrosericus",colnames(scrb.NM))] <-"Eurysternus.atrosericus"
colnames(scrb.NM)[grep("convexa",colnames(scrb.NM))] <- "Tetramereia.convexus"
colnames(scrb.NM)[grep("viridanum",colnames(scrb.NM))] <- "Oxysternon.festivum.ssp.viridanum"

spps.NM <- colnames(scrb.NM)
genero <- unname(sapply(spps.NM,function(x) {strsplit(x,"\\.")[[1]][1]}))
epiteto <- unname(sapply(spps.NM,function(x) {paste(strsplit(x,"\\.")[[1]][-1],collapse=" ")}))

NM.scarab <- Scarab[match(spps.NM,paste(Scarab$genero,Scarab$epiteto,sep=".")),]
NM.scarab$status <- "A"
tmp <- Scarab[match(genero[is.na(NM.scarab$cdg_genero)],Scarab$genero),]
tmp$cdg_especie <- NA
tmp$epiteto <- epiteto[is.na(NM.scarab$cdg_genero)]
tmp$status <- "B"
tmp$sintax <- "<i><G></i> sp. [code <E>]"
NM.scarab[is.na(NM.scarab$cdg_genero),] <- tmp

NM.scarab <- NM.scarab[order(NM.scarab$tribu,NM.scarab$subtribu,NM.scarab$genero,NM.scarab$status,NM.scarab$epiteto),]

NM.scarab$sintax[grep("guildingii",NM.scarab$epiteto)] <- "<i><G> <E></i> [different from <i>D. icarus</i>]"
NM.scarab$sintax[grep("viridanum",NM.scarab$epiteto)] <- "<i><G> <E></i>"


###################################################
### chunk number 15: 
###################################################
#line 391 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
cat(paste("\\begin{itemize}"))

for (j in unique(NM.scarab$trb)) {
	cat(sprintf("\\item Tribu %s",j))
	strb <- unique(NM.scarab$strb[NM.scarab$trb %in% j])
	cat(sprintf("\\begin{itemize}"))
	for (k in strb) {
		if (!is.na(k)) {
			cat(sprintf("\\item Subribu %s",k))
		}	
		cat(paste("\\begin{enumerate}"))
		for (i in seq(along=NM.scarab$trb)[NM.scarab$trb %in% j & NM.scarab$strb %in% k]) {
			cat(paste("\n\\item", sub("<i>", "\\\\textit{",sub("</i>","}", sub("<A>",sub("&","et",NM.scarab$autor[i]), sub("<F>",NM.scarab$fecha[i], sub("<E>",NM.scarab$epiteto[i], sub("<G>",NM.scarab$genero[i],NM.scarab$sintax[i]))))))))
		}
		cat(paste("\n\\end{enumerate}\n"))
	}
	cat(sprintf("\\end{itemize}"))
}
cat(paste("\n\\end{itemize}\n"))



###################################################
### chunk number 16: 
###################################################
#line 418 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
tt <- aggregate(avs.NM$dur,list(NM=avs.NM$NM),sum)
tt$sfrz <- tt$x/60

t0 <- merge(merge(merge(aggregate(NM.m1[,c("ID", "Especieid")], by=list(NM=NM.m1$IDTransecta), luq), aggregate(NM.m1[,c("Ntotal")], by=list(NM=NM.m1$IDTransecta),sum), by="NM"),
aggregate(NM.m2[,c("ID","Especieid")],by=list(NM=NM.m2$IDTransecta),luq),  by="NM"),
aggregate(c(NM.m2$Especieid,NM.m1$Especieid), by=list(NM=c(as.character(NM.m2$IDTransecta),as.character(NM.m1$IDTransecta))),luq), by="NM")

colnames(t0) <- c("NM","nreg.m1","nspp.m1","nind.m1","nreg.m2","nspp.m2","nspp")
t0$NM <- trans.info$NM[match(t0$NM,trans.info$IDTransecta)]

t0$sfrz.m1 <- 150/60
t0$sfrz.m2 <- 90/60

t0 <- t0[order(t0$NM),c("NM","sfrz.m1","nreg.m1","nind.m1", "nspp.m1", "sfrz.m2","nreg.m2", "nspp.m2", "nspp")]



###################################################
### chunk number 17: 
###################################################
#line 440 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
 print(xtable(t0, caption="Transects included in NeoMaps' Bird survey from 2010. NM: NeoMaps transect code (see table \\ref{TAB:RSMINFO} for details), $E$: sampling effort in hours*observer (h*p), $C$: number of unique records, $N$: total number of individuals counted (first day sample), $S_{obs}$: number of species detected, $S_{total}$: total number of species detected in both days. Summary results are presented for a subset of transects: estimated mean and standard deviation of dung beetle species abundance in log-scale  ($\\log N (\\mu \\pm \\sigma)$), Shannon-Weaver diversity index, ($H$), and abundance coverage estimator of dung beetle species richness (ACE) for NeoMaps transects.", label="TAB:AVESRSM", align="|l|c|cccc|ccc|c|", digits=c(0,0,1,0,0,0,1,0,0,0)), caption.placement="top", include.rownames=F, include.colnames=F, add.to.row=list(pos=list(-1), command=c("
\\multicolumn{1}{c}{} & \\multicolumn{4}{c}{First day count} &  \\multicolumn{3}{c}{Second day count} & \\multicolumn{1}{c}{} \\\\\\hline\n
 NM & $E$ (h*p) & $C$ & $N$ & $S_{obs}$ &  $E$ (h*p) & $C$ & $S_{obs}$& $S_{total}$ \\\\")), table.placement="H", hline.after=c(-1,nrow(t0)), size="normalsize")



###################################################
### chunk number 18: MTRZaves
###################################################
#line 454 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"

d1.aves <- vegdist(aves.NM,"chao")
d2.aves <- dist.binary(aves.LIT.um,method=5)## Sørensen

mx.aves <- as.matrix(d1.aves)

mx.aves <- mx.aves[order(trans.info[match(rownames(aves.NM),trans.info$Nombretransecta),"Bioregión"]),order(trans.info[match(rownames(aves.NM),trans.info$Nombretransecta),"Bioregión"])]

rownames(mx.aves) <- trans.info[match(rownames(mx.aves),trans.info$Nombretransecta),"NM"]
colnames(mx.aves) <- trans.info[match(colnames(mx.aves),trans.info$Nombretransecta),"NM"]

tmp <- as.matrix(d2.aves)
mx.aves[upper.tri(mx.aves)] <- NA #tmp[upper.tri(tmp)]


par(mar=c(3,3,3,0),oma=c(0,0,0,0),xpd=T)
image(1:nrow(mx.aves),1:ncol(mx.aves),t(mx.aves), axes=F, col=c(NA,rev(mis.cols)), xlab="", ylab="")

polygon(cumsum(table(trans.info[match(rownames(mx.aves),trans.info$NM),"Bioregión"]))[c(1,2,2,3,3,4,4,5,5,6,6,5,5,4,4,3,3,2,2,1,1)]+.5,cumsum(table(trans.info[match(rownames(mx.aves),trans.info$NM),"Bioregión"]))[c(1,1,2,2,3,3,4,4,5,5,6,6,5,5,4,4,3,3,2,2,1)]+.5,lwd=2)

atta <- aggregate(1:ncol(mx.aves),by=list(trans.info[match(rownames(mx.aves),trans.info$NM),"Bioregión"]),mean)$x
text(atta,c(-.5,4.5,8.5,13.5,19.5),c("Occident","Andean\nmountains","Coastal\nmountains","Orinoco\nfloodplain","Guayana\nshield"),cex=1)

axis(3,las=2,1:ncol(mx.aves),rownames(mx.aves),cex.axis=1)
axis(2,las=2,1:ncol(mx.aves),rownames(mx.aves),cex.axis=1)

text(rep(1:nrow(mx.aves),ncol(mx.aves)), rep(1:nrow(mx.aves),rep(ncol(mx.aves),ncol(mx.aves))) ,round(as.vector(t(mx.aves)),2),cex=.5)



###################################################
### chunk number 19:  eval=FALSE
###################################################
## #line 497 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
## require(ROpenOffice)
## load("~/NeoMapas/Rdata/RevisionLiteratura.rda")
## if(!exists("lamas")) {
## 	lamas <- read.ods("~/NeoMapas/etc/1_Coleccion/130_Taxonomia/Lamas2004_2011_06_22_Todo.ods",TRUE,stringsAsFactors=F)
## 	}
## 	tbl.mrps <- data.frame()
## 	cols <- colnames(lamas[[2]])
## 	for (i in names(lamas)[-1]) {
## 		tt <- lamas[[i]]
## 		colnames(tt) <- cols
## 		tbl.mrps <- rbind(tbl.mrps, tt)
## 	}
## 
## 
## spps <- tbl.mrps[tbl.mrps$scdg %in% "a" & tbl.mrps$ecdg %in% mrps.cneb$cdg[mrps.cneb$validacion %in% c("exacto", "ortografico", "soundex epiteto",  "soundex genero", "manual") & mrps.cneb$adm0 %in% "VEN" & grepl("^097-",mrps.cneb$cdg)],]
## 
## spps <- spps[order(spps$ecdg),]


###################################################
### chunk number 20:  eval=FALSE
###################################################
## #line 522 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
## cat(sprintf("\\nocite{%s}\n\n",paste(unique(mrps.cneb$cdg_ref[mrps.cneb$validacion %in% c("exacto", "ortografico", "soundex epiteto",  "soundex genero", "manual") & mrps.cneb$adm0 %in% "VEN" & grepl("^097-",mrps.cneb$cdg)]),collapse=",")))


###################################################
### chunk number 21:  eval=FALSE
###################################################
## #line 526 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
## cat(paste("\\begin{itemize}"))
## 
## for (j in unique(spps$sfam)) {
## 	cat(sprintf("\\item Subfamilia %s",j))
## 	strb <- unique(spps$trb[spps$sfam %in% j])
## 	cat(sprintf("\\begin{itemize}"))
## 	for (k in strb) {
## 		if (!is.na(k) & k != "NA") {
## 			cat(sprintf("\\item Tribu %s",k))
## 		}	
## 		cat(paste("\\begin{enumerate}"))
## 		for (i in seq(along=spps$trb)[spps$sfam %in% j & spps$trb %in% k]) {
## 			cat(paste("\n\\item", sub("<i>", "\\\\textit{",sub("</i>","}", sub("<A>",sub("&","\\&",spps$aut[i]), sub("<F>",spps$fch[i], sub("<E>",spps$esp[i], sub(" <S>","",  sub("<G>",spps$gen[i],spps$stx[i])))))))))
## 
## 			tmp <- unique(scrb.cneb[scrb.cneb$ecdg %in% spps$cdg_especie[i],c("cdg_ref", "especie", "ecdg", "validacion", "cneb")])
## 			mi.ss <- tmp$validacion %in% c("vld-auto :: concordancia exacta", "heredada","manual")
## 
## 			if (sum(mi.ss)>0) {
## 				cns <- unique(tmp$cneb)
## 				cns <- cns[!is.na(cns)]
## 				tmp2 <- aggregate(tmp$especie[mi.ss], by=list(tmp$cdg_ref[mi.ss]), paste, collapse="; ")
## 				if (length(cns)==0) {
## 					cat(paste("\n\\newline {\\small Records without detailled geographic information."))
## 				} else {
## 					cat(paste("\n\\newline {\\small Distribution: ",paste(cns,collapse=", "),".",sep=""))
## 				}
## 
## 			cat(paste("\n\\newline References: "))
## 	
## 				cat(paste("\\citet{",tmp2[,1],"}", sep="",collapse="; "))
## 				cat(paste("}"))
## 				
## 			}
## 
## 			if (sum(!mi.ss)>0) {
## 				tmp3 <- aggregate(tmp$especie[!mi.ss], by=list(tmp$cdg_ref[!mi.ss]), paste, collapse="; ")
## 				cat(paste("\n{\\textcolor{red}{Nombres por verificar: "))
## 				cat(paste("\\textbf{",tmp3[,1],"}: ", tmp3[,2], sep="",collapse="; "))
## 				cat(paste("}"))
## 			}
## 		}
## 		cat(paste("\n\\end{enumerate}\n"))
## 	}
## 	cat(sprintf("\\end{itemize}"))
## }
## cat(paste("\n\\end{itemize}\n"))
## 


###################################################
### chunk number 22:  eval=FALSE
###################################################
## #line 582 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
## require(ROpenOffice)
## if(!exists("tbl.scrb")) {
## 	tbl.scrb <- read.ods(file="~/NeoMapas/etc/1_Coleccion/130_Taxonomia/20110701_Scarabaeinae.ods",TRUE,stringsAsFactors=F)
## }
## 	load("~/NeoMapas/Rdata/RevisionLiteratura.rda")
## 
## Scarab <- tbl.scrb[["Taxonomia"]]
## 
## spps <- scrb.vzla <- Scarab[Scarab$cdg_especie %in% scrb.cneb$ecdg[scrb.cneb$validacion %in% c("heredada","manual","vld-auto :: concordancia exacta") & scrb.cneb$adm0 %in% "VEN"],]
## 
## spps <- spps[order(spps$tribu,spps$subtribu,spps$genero,spps$epiteto),]


###################################################
### chunk number 23:  eval=FALSE
###################################################
## #line 600 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
## cat(sprintf("\\nocite{%s}\n\n",paste(unique(scrb.cneb$cdg_ref[scrb.cneb$validacion %in% c("heredada","manual","vld-auto :: concordancia exacta") & scrb.cneb$adm0 %in% "VEN"]),collapse=",")))


###################################################
### chunk number 24:  eval=FALSE
###################################################
## #line 604 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento3_AnexoArticuloNeoMapas.Rnw"
## cat(paste("\\begin{itemize}"))
## 
## for (j in unique(spps$trb)) {
## 	cat(sprintf("\\item Tribu %s",j))
## 	strb <- unique(spps$strb[spps$trb %in% j])
## 	cat(sprintf("\\begin{itemize}"))
## 	for (k in strb) {
## 		if (!is.na(k)) {
## 			cat(sprintf("\\item Subribu %s",k))
## 		}	
## 		cat(paste("\\begin{enumerate}"))
## 		for (i in seq(along=spps$trb)[spps$trb %in% j & spps$strb %in% k]) {
## 			cat(paste("\n\\item", sub("<i>", "\\\\textit{",sub("</i>","}", sub("<A>",sub("&","et",spps$autor[i]), sub("<F>",spps$fecha[i], sub("<E>",spps$epiteto[i], sub("<G>",spps$genero[i],spps$sintax[i]))))))))
## 
## 			tmp <- unique(scrb.cneb[scrb.cneb$ecdg %in% spps$cdg_especie[i],c("cdg_ref", "especie", "ecdg", "validacion", "cneb")])
## 			mi.ss <- tmp$validacion %in% c("vld-auto :: concordancia exacta", "heredada","manual")
## 
## 			if (sum(mi.ss)>0) {
## 				cns <- unique(tmp$cneb)
## 				cns <- cns[!is.na(cns)]
## 				tmp2 <- aggregate(tmp$especie[mi.ss], by=list(tmp$cdg_ref[mi.ss]), paste, collapse="; ")
## 				if (length(cns)==0) {
## 					cat(paste("\n\\newline {\\small Records without detailled geographic information."))
## 				} else {
## 					cat(paste("\n\\newline {\\small Distribution: ",paste(cns,collapse=", "),".",sep=""))
## 				}
## 
## 			cat(paste("\n\\newline References: "))
## 	
## 				cat(paste("\\citet{",tmp2[,1],"}", sep="",collapse="; "))
## 				cat(paste("}"))
## 				
## 			}
## 
## 			if (sum(!mi.ss)>0) {
## 				tmp3 <- aggregate(tmp$especie[!mi.ss], by=list(tmp$cdg_ref[!mi.ss]), paste, collapse="; ")
## 				cat(paste("\n{\\textcolor{red}{Nombres por verificar: "))
## 				cat(paste("\\textbf{",tmp3[,1],"}: ", tmp3[,2], sep="",collapse="; "))
## 				cat(paste("}"))
## 			}
## 		}
## 		cat(paste("\n\\end{enumerate}\n"))
## 	}
## 	cat(sprintf("\\end{itemize}"))
## }
## cat(paste("\n\\end{itemize}\n"))
## 


