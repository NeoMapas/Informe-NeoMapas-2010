###################################################
### chunk number 1: 
###################################################
#line 61 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
require(foreign)
require(vegan)
require(vegan)
require(cluster)
require(xtable)
require(labdsv)

load(file="~/NeoMapas/Rdata/SIG.rda")
##VBG <- read.dbf("/var/local/gis/mi.gis/Venezuela/NeoMapas/dbf/VBG.dbf")
VBG <- read.dbf("/Users/jferrer/mi.gis/Venezuela/NeoMapas/dbf/VBG.dbf")
CNEB.nm <- merge(CNEB@data,VBG, by=1:5)



###################################################
### chunk number 2: 
###################################################
#line 78 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

rngs <- CNEB.nm[,grep("max",colnames(CNEB.nm),value=T)]-CNEB.nm[,grep("min",colnames(CNEB.nm),value=T)]

colnames(rngs) <- sub("max","rng",colnames(rngs))
CNEB.nm <- cbind(CNEB.nm,rngs)

vars.nm <- c("lon", "lat", "altr_avrg", "altr_rng", "tmpm_avrg", "tmpm_rng", "prcp_avrg", "prcp_rng", "mscs_avrg", "mscs_rng", "bosq_avrg", "bosq_rng", "decd_avrg", "decd_rng")
vars.nm <- c("altr_avrg", "prcp_rng","prcp_avrg", "tmpm_rng", "mscs_avrg", "mscs_rng", "bosq_avrg", "bosq_rng", "decd_avrg", "decd_rng")

mi.cneb <- CNEB.nm[CNEB.nm$UM==1,vars.nm]
rownames(mi.cneb) <- CNEB.nm[CNEB.nm$UM==1,"cdg"]
VIFs <- c()
for (mv in vars.nm) {
	VIFs <- c(VIFs,(sqrt(1/(1-summary(step(lm(formula(paste(mv,"~",".")),mi.cneb)))$adj.r.squared))))
}
VIFs

env.pc <- rda(mi.cneb,scale=T)
##env.dudi <- dudi.pca(mi.cneb,scannf=F,nf=3)
CNEB.nm$PC1 <- CNEB.nm$PC2 <- CNEB.nm$PC3 <- numeric(nrow(CNEB.nm))
CNEB.nm[CNEB.nm$UM==1,c("PC1","PC2","PC3")] <- scores(env.pc, choices=1:3,display="sites")



###################################################
### chunk number 3: 
###################################################
#line 106 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

## PC1 altura precipitacion y temperatura
## PC2 cobertura boscosa
## PC3 meses secos
cor(scores(env.pc,choice=1:3,display="wa"),mi.cneb)



###################################################
### chunk number 4: 
###################################################
#line 115 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
require(randomForest)
mi.rf <- randomForest(y=as.factor(CNEB.nm$est[CNEB.nm$UM==1 & CNEB.nm$est>0]),x=CNEB.nm[CNEB.nm$UM==1  & CNEB.nm$est>0,15:47])
CNEB.nm$est.pr <- predict(mi.rf,CNEB.nm)

CNEB.nm$est.pr[CNEB.nm$UM==1 & CNEB.nm$est>0] <- CNEB.nm$est[CNEB.nm$UM==1 & CNEB.nm$est>0]



###################################################
### chunk number 5: 
###################################################
#line 126 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

##CNEB.nm$bioreg <- factor(CNEB.nm$bioreg)
CNEB.nm$bioreg <- factor(CNEB.nm$bioreg,labels=c("Otro","Occident","Andean mountains","Coastal mountains","Llanos","Guayana"),levels=c(0,1,5,2,3,4))

##CNEB.nm$bioreg <- factor(CNEB.nm$bioreg,labels=c("Otro","Occidente","Centro y Costa","Orinoco floodplain","Guayana shield","Andes"))
CNEB.nm$est.pr <- factor(CNEB.nm$est.pr)


## combinaciones de estratos y bioregiones 
dim(unique(CNEB.nm[CNEB.nm$UM==1,c("bioreg","est.pr")]))

table(CNEB.nm[CNEB.nm$UM==1,c("bioreg","est.pr")])


###################################################
### chunk number 6: 
###################################################
#line 142 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
 xtable(rbind(env.pc$CA$v[,1:3],env.pc$CA$eig[1:3]))


###################################################
### chunk number 7: 
###################################################
#line 151 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
## datos de las bases de datos de NeoMapas
load(file="~/NeoMapas/Rdata/NM.rda")

## datos de identificaciones de escarabajos de Solís
tmp <- read.csv("~/NeoMapas/data/AngelSolisDatos/DatosScarabSetiembre2010.csv", sep="\t", header=T, as.is=T,nrow=5)
scrb.solis <- read.csv("~/NeoMapas/data/AngelSolisDatos/DatosScarabSetiembre2010.csv", sep="\t", header=F, as.is=T, skip=2)
colnames(scrb.solis) <- c("NM", "Prd", "Fecha", "Caja", "ADN", colnames(tmp)[6:ncol(tmp)])
##head(tmp)
scrb.solis$yr <- sapply(scrb.solis$Fecha,function(x){strsplit(x,"/")[[1]][3]},simplify=T)
scrb.solis$ttl <- rowSums(scrb.solis[,colnames(tmp)[6:(ncol(tmp)-1)]])
scrb.solis <- scrb.solis[!is.na(scrb.solis$yr) & !(scrb.solis$NM %in% c("1000IVIC","41")),]

## datos de observación de mariposas de Gustavo Rodríguez
NM.m1 <- read.csv("~/NeoMapas/data/NMAves2010/Muestreo1.total.csv")
NM.m2 <- read.csv("~/NeoMapas/data/NMAves2010/Muestreo2.total.csv")
spp.aves <- read.csv("~/NeoMapas/data/NMAves2010/Venezuela.total.csv",as.is=T)
trans.info <- read.csv("~/NeoMapas/data/NMAves2010/Transectas.total.csv",as.is=T)

mi.ss <- NM.m1$Especieid!=1379
mi.s2 <- NM.m2$Especieid!=1379
##aves.NM <- table(NM.m1$IDTransecta[mi.ss],NM.m1$Especieid[mi.ss])

##aves.NM <- tapply(c(NM.m1$Ntotal[mi.ss],NM.m2$Ntotal[mi.s2]),list(c(NM.m1$IDTransecta[mi.ss],NM.m2$IDTransecta[mi.s2]),c(NM.m1$Especieid[mi.ss],NM.m2$Especieid[mi.s2])),sum,na.rm=T)
aves.NM <- tapply(c(NM.m1$Ntotal[mi.ss]),list(NM.m1$IDTransecta[mi.ss],NM.m1$Especieid[mi.ss]),sum,na.rm=T)
aves.NM[is.na(aves.NM)] <- 0


###################################################
### chunk number 8: 
###################################################
#line 181 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

luq <- function(x){length(unique(x))}
tapply(tvn.NM$NM,tvn.NM$yr,luq)
tapply(tvn.NM$CN,tvn.NM$yr,luq)
tapply(tvn.NM$sfrz,tvn.NM$yr,sum,na.rm=T)
tapply(jmp.NM$jmp,jmp.NM$yr,luq)
tapply(obs.NM$jmp,substr(obs.NM$fch,1,4),length)
tapply(obs.NM$NM,substr(obs.NM$fch,1,4),luq)
tapply(trmp.NM$NM,trmp.NM$yr,luq)
tapply(trmp.NM$sfrz,trmp.NM$yr,sum)

sfrz <- jmps <- nms <- sfrze <- nme <- jmpe <- sfrza <- nma <- jmpa <- c()
sfrz["2005"] <- sum(tvn.NM$sfrz[tvn.NM$yr %in% 2003:2005],na.rm=T)/60
sfrz["2006"] <- sum(tvn.NM$sfrz[tvn.NM$yr %in% 2006],na.rm=T)/60
sfrz["2009"] <- sum(tvn.NM$sfrz[tvn.NM$yr %in% 2009:2010],na.rm=T)/60

sfrze["2005"] <- sum(trmp.NM$sfrz[trmp.NM$yr %in% 2003:2005],na.rm=T)
sfrze["2006"] <- sum(trmp.NM$sfrz[trmp.NM$yr %in% 2006],na.rm=T)
sfrze["2009"] <- sum(trmp.NM$sfrz[trmp.NM$yr %in% 2009:2010],na.rm=T)

nms["2005"] <- luq(tvn.NM$NM[tvn.NM$yr %in% 2003:2005])
nms["2006"] <- luq(tvn.NM$NM[tvn.NM$yr %in% 2006])
nms["2009"] <- luq(tvn.NM$NM[tvn.NM$yr %in% 2009:2010])

nme["2005"] <- luq(trmp.NM$NM[trmp.NM$yr %in% 2003:2005])
nme["2006"] <- luq(trmp.NM$NM[trmp.NM$yr %in% 2006])
nme["2009"] <- luq(trmp.NM$NM[trmp.NM$yr %in% 2009:2010])

jmps["2005"] <- length(jmp.NM$jmp[jmp.NM$tvn %in% tvn.NM$tvn[tvn.NM$yr %in% 2003:2005]])
jmps["2006"] <- length(jmp.NM$jmp[jmp.NM$tvn %in% tvn.NM$tvn[tvn.NM$yr %in% 2006]])
jmps["2009"] <- length(jmp.NM$jmp[jmp.NM$tvn %in% tvn.NM$tvn[tvn.NM$yr %in% 2009:2010]])


jmpe["2005"] <- sum(scrb.solis[scrb.solis$yr %in% "05",6:(ncol(scrb.solis)-3)],na.rm=T)+nrow(scr.NM[scr.NM$yr =="2005" & scr.NM$NM %in% c("09","24","26"),])

jmpe["2006"] <- sum(scrb.solis[scrb.solis$yr %in% "06",6:(ncol(scrb.solis)-3)],na.rm=T)
jmpe["2009"] <- sum(scrb.solis[scrb.solis$yr %in% "09",6:(ncol(scrb.solis)-3)],na.rm=T)


sfrza["2002"] <- sum(avs.NM$sfrz[substr(avs.NM$fch,1,4) %in% c("2001","2002")])/60

sfrza["2010"] <- ((nrow(unique(NM.m1[,c("IDTransecta","Punto")]))*3)+(nrow(unique(NM.m2[,c("IDTransecta","Punto")]))*9))/60

nma["2002"] <- luq(avs.NM$NM[substr(avs.NM$fch,1,4) %in% c("2001","2002")])
nma["2010"] <- luq(NM.m1$IDTransecta)

jmpa["2010"] <- nrow(unique(NM.m1))+nrow(unique(NM.m2))

jmpa["2002"] <- nrow(obs.NM)

sfrze <- round(sfrze,1)
sfrz <- round(sfrz,1)
sfrza <- round(sfrza,1)
##esfuerzo aves 2010
((nrow(unique(NM.m1[,c("IDTransecta","Punto")]))*3)+(nrow(unique(NM.m2[,c("IDTransecta","Punto")]))*9))/60
nrow(NM.m1)+nrow(NM.m2)
nrow(NM.m1)
dim(aves.NM)



###################################################
### chunk number 9: 
###################################################
#line 260 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

mi.jmp <- jmp.NM[(jmp.NM$tvn %in% tvn.NM$tvn[tvn.NM$yr %in% "2006"]),]
mi.jmp <- mi.jmp[mi.jmp$NM!="92" & !is.na(mi.jmp$NM),]
mi.jmp$CN <- info.NM$CNEB[match(as.numeric(mi.jmp$NM),info.NM$NM)]

mi.jmp$bioreg <- CNEB.nm[match(mi.jmp$CN,CNEB.nm$cdg),"bioreg"]

mi.tt <- table(mi.jmp$bioreg,mi.jmp$familia,useNA="ifany")
colnames(mi.tt) <- c("Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae",  "Riodinidae", "No identificado")

fams <- c("Hesperiidae", "Papilionidae", "Pieridae", "Lycaenidae", "Riodinidae", "Nymphalidae","No identificado")

mi.tt <- mi.tt[,fams]
colnames(mi.tt) <- abbreviate(fams)
mi.tt <- cbind(mi.tt,total=rowSums(mi.tt))
mi.tt <- rbind(mi.tt,total=colSums(mi.tt))
mi.tt <- mi.tt[c(2,6,3,4,5,7),]


##tapply(mi.jmp$especie,mi.jmp$familia,function(x){length(unique(x))})
xtable(mi.tt,digits=0,caption=paste(abbreviate(unique(fams)),": ",unique(fams),sep="", collapse="; "))


###################################################
### chunk number 10: 
###################################################
#line 285 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
mi.tvn <- tvn.NM[tvn.NM$yr %in% "2006",]
mi.tvn <- mi.tvn[mi.tvn$NM!="92",]
mi.tvn$CN <- info.NM$CNEB[match(as.numeric(mi.tvn$NM),info.NM$NM)]

mi.tvn$bioreg <- CNEB.nm[match(mi.tvn$CN,CNEB.nm$cdg),"bioreg"]

tapply(mi.tvn$vst,mi.tvn$bioreg, function(x){length(unique(x))})
mrps.sfrz <- tapply(mi.tvn$sfrz,mi.tvn$NM,sum,na.rm=T)/60



###################################################
### chunk number 11: 
###################################################
#line 297 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
mi.jmp <- jmp.NM[(jmp.NM$tvn %in% tvn.NM$tvn[tvn.NM$yr %in% "2006"])&jmp.NM$fam=="Pieridae" & !is.na(jmp.NM$genero),]
mi.jmp$esp <- paste(mi.jmp$genero,mi.jmp$especie)
mi.jmp <- mi.jmp[mi.jmp$NM!="92",]

mrps.NM <- tapply(mi.jmp$jmp,list(mi.jmp$NM,mi.jmp$esp),function(x){length(unique(x))})
mrps.NM[is.na(mrps.NM)] <- 0
mrps.NM <- mrps.NM[,colSums(mrps.NM)>0]
mrps.CN <- mrps.NM
rownames(mrps.NM) <- info.NM$Nombre[match(as.numeric(rownames(mrps.NM)),info.NM$NM)]
rownames(mrps.CN) <- info.NM$CNEB[match(as.numeric(rownames(mrps.CN)),info.NM$NM)]


###################################################
### chunk number 12: 
###################################################
#line 313 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"



scrb.06 <- rowsum(scrb.solis[scrb.solis$yr %in% "06",6:(ncol(scrb.solis)-4)],
                  group=scrb.solis$NM[scrb.solis$yr %in% "06"])
scrb.09 <- rowsum(scrb.solis[scrb.solis$yr %in% c("09","10"),6:(ncol(scrb.solis)-4)],
                  group=scrb.solis$NM[scrb.solis$yr %in%  c("09","10")])
scrb.NM <- rbind(scrb.09,
                 scrb.06[rownames(scrb.06) %in% c(13,16,18,2,27,5),])

tmp001 <- tapply(trmp.NM$sfrz,list(trmp.NM$NM,trmp.NM$yr),sum)
scrb.sfrz <- c(tmp001[sub(" ","0",sprintf("%02s",rownames(scrb.09))),"2009"],tmp001[c(13,16,18,"02",27,"05"),"2006"])

scrb.sfrz["15"] <- sum(trmp.NM[trmp.NM$NM %in% "15" & trmp.NM$fch %in% "2009-08-12","sfrz"])
scrb.sfrz["66"] <- tmp001[c("66"),"2010"]


scrb.NM.c <- (scrb.NM/scrb.sfrz)



##rownames(scrb.06) <- info.NM$Nombre[match(rownames(scrb.06),info.NM$NM)
scrb.CN <- info.NM$CNEB[match(rownames(scrb.NM),info.NM$NM)]

rownames(scrb.NM) <- info.NM$Nombre[match(rownames(scrb.NM),info.NM$NM)]

dim(scrb.NM[,colSums(scrb.NM)>0])



###################################################
### chunk number 13: 
###################################################
#line 345 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
load("~/NeoMapas/Rdata/LIT.rda")
load("~/NeoMapas/Rdata/MSM.rda")

scrb.cneb <- CNEB.nm[match(scrb.CN,CNEB.nm$cdg),]
table(scrb.cneb[,c("bioreg","est.pr")])

fams <- c()
for (i in colnames(scrb.NM)) {
  fams <- c(fams,unique(Escarabajos$trb[match(strsplit(i,"\\.")[[1]][1],Escarabajos$gen)]))
}
fams[grep("Genieridium",colnames(scrb.NM))] <- "Dichotomiini"
fams[grep("Onoreidium",colnames(scrb.NM))] <- "Dichotomiini"
fams[grep("Trichillidium",colnames(scrb.NM))] <- "Dichotomiini"
fams[grep("Oxysternum",colnames(scrb.NM))] <- "Phanaeini"
fams[grep("Tetramerteia",colnames(scrb.NM))] <- "Phanaeini"
fams[grep("Euristernus",colnames(scrb.NM))] <- "Eurysternini"
mi.tt <- data.frame()
for (i in unique(fams)) {
  mi.tt <- rbind(mi.tt,i=rowSums( rowsum(scrb.NM[,fams %in% i],scrb.cneb$bioreg)))
}
mi.tt <- t(mi.tt)
colnames(mi.tt) <- abbreviate(unique(fams))
rownames(mi.tt) <- levels(scrb.cneb$bioreg)[-1]
mi.tt <- cbind(mi.tt,total=rowSums(mi.tt))
mi.tt <- rbind(mi.tt,total=colSums(mi.tt))
mi.tt <- mi.tt[c(1,5,2,3,4,6),]
xtable(mi.tt,digits=0,caption=paste(abbreviate(unique(fams)),": ",unique(fams),sep="", collapse="; "))


###################################################
### chunk number 14: 
###################################################
#line 376 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
sfrz <- round(tapply(trmp.NM$sfrz,list(trmp.NM$NM,trmp.NM$yr),sum))
sfrz[is.na(sfrz)] <- 0
tmp <- info.NM[match(as.numeric(rownames(sfrz)),info.NM$NM),"CNEB"]
tmp2 <- CNEB.nm$bioreg[match(tmp,CNEB.nm$cdg)]

xtable(rowsum(sfrz,tmp2),digits=0)

sfrz <- cbind(CNEB=tmp,bioreg=as.character(tmp2),sfrz)

xtable(sfrz[order(tmp2),],digits=0)

sfrz <- round(tapply(trmp.NM$vst,list(trmp.NM$NM,trmp.NM$yr),function(x){length(unique(x))}))
sfrz[is.na(sfrz)] <- 0
tmp <- info.NM[match(as.numeric(rownames(sfrz)),info.NM$NM),"CNEB"]
tmp2 <- CNEB.nm$bioreg[match(tmp,CNEB.nm$cdg)]

xtable(rowsum(sfrz,tmp2),digits=0)

sfrz <- cbind(CNEB=tmp,bioreg=as.character(tmp2),sfrz)
xtable(sfrz[order(tmp2),],digits=0)




###################################################
### chunk number 15: 
###################################################
#line 403 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

mi.ss <- NM.m1$Especieid!=1379
mi.s2 <- NM.m2$Especieid!=1379
##aves.NM <- table(NM.m1$IDTransecta[mi.ss],NM.m1$Especieid[mi.ss])

##aves.NM <- tapply(c(NM.m1$Ntotal[mi.ss],NM.m2$Ntotal[mi.s2]),list(c(NM.m1$IDTransecta[mi.ss],NM.m2$IDTransecta[mi.s2]),c(NM.m1$Especieid[mi.ss],NM.m2$Especieid[mi.s2])),sum,na.rm=T)
aves.NM <- tapply(c(NM.m1$Ntotal[mi.ss]),list(NM.m1$IDTransecta[mi.ss],NM.m1$Especieid[mi.ss]),sum,na.rm=T)
aves.NM[is.na(aves.NM)] <- 0

colnames(aves.NM) <- spp.aves$Latname[match(colnames(aves.NM),spp.aves$N_ave)]
aves.CN <- aves.NM
rownames(aves.NM) <- trans.info[match(rownames(aves.NM),trans.info$IDTransecta),"Nombretransecta"]

rownames(aves.CN) <- info.NM$CNEB[match(rownames(aves.NM),info.NM$Nombre)]




mi.ss <- NM.m2$Especieid!=1379
m2 <- table(NM.m2$IDTransecta[mi.ss],NM.m2$Especieid[mi.ss])


###################################################
### chunk number 16: 
###################################################
#line 428 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
load(file="~/NeoMapas/Rdata/MSM.rda")
fams <- ordo <- c()
for (i in colnames(aves.NM)) {
  fams <- c(fams,unique(Aves$fam[match(strsplit(i,"\\ ")[[1]][1],Aves$gen)]))
  ordo <- c(ordo,unique(Aves$ordo[match(strsplit(i,"\\ ")[[1]][1],Aves$gen)]))
}


colnames(aves.NM)[is.na(fams)]
fams[grep("Columba",colnames(aves.NM))] <- "Columbidae"
fams[grep("Ajaia",colnames(aves.NM))] <- "Threskiornithidae"
fams[grep("Megaceryle",colnames(aves.NM))] <- "Cerylidae"
fams[grep("Pelecanus",colnames(aves.NM))] <- "Pelecanidae"
fams[fams=="dontophoridae"] <- "Odontophoridae"
fams[grep("Pachyramphus",colnames(aves.NM))] <- "Cotingidae"
fams[grep("Laniocera",colnames(aves.NM))] <- "Cotingidae"
fams[grep("Piprites",colnames(aves.NM))] <- "Cotingidae"


ordo[fams %in% c("Pelecanidae", "Cerylidae")] <- "Pelecaniformes"

for (i in unique(c(fams[ordo==""],fams[is.na(ordo)]))) {
  ordo[fams==i] <- unique(Aves$ordo[match(i,Aves$fam)])
}
unique(fams[ordo==""])

aves.cneb <- CNEB.nm[match(rownames(aves.CN),CNEB.nm$cdg),]

ordo[ordo %in% c('Tinamiformes','Pelecaniformes','Anseriformes','Craciformes','Galliformes','Gruiformes','Charadriiformes','Cuculiformes','Strigiformes','Trogoniformes','Coraciiformes','Galbuliformes','Piciformes')] <- "Otros"
mi.tt <- data.frame()
for (i in unique(ordo)) {
  mi.tt <- rbind(mi.tt,i=rowSums( rowsum(aves.NM[,ordo %in% i],aves.cneb$bioreg)))
}
mi.tt <- t(mi.tt)
colnames(mi.tt) <- abbreviate(unique(ordo))
##colnames(mi.tt) <-  unique(ordo)
rownames(mi.tt) <- levels(aves.cneb$bioreg)[-1]
mi.tt <- cbind(mi.tt,total=rowSums(mi.tt))
mi.tt <- rbind(mi.tt,total=colSums(mi.tt))
mi.tt <- mi.tt[c(1,5,2,3,4,6),]
xtable(mi.tt,digits=0,caption=paste(paste(abbreviate(unique(ordo)),": ",unique(ordo),sep="", collapse="; "),". Otros=Tinamiformes, Pelecaniformes, Anseriformes, Craciformes, Galliformes, Gruiformes, Charadriiformes, Cuculiformes, Strigiformes, Trogoniformes, Coraciiformes, Galbuliformes, Piciformes",sep=""))


###################################################
### chunk number 17: 
###################################################
#line 474 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
tapply(scrb.solis$ttl,scrb.solis$yr,sum)



## escarabajos
sum(scrb.06[rownames(scrb.06) %in% c(13,16,18,2,27,5),])
dim(scrb.09)
 sum(scrb.09)

##mariposas
mi.tt <- table(jmp.NM$status_id %in% c("definitiva","preliminar"),jmp.NM$yr)
mi.tt/colSums(mi.tt)
round(t(mi.tt)*100/colSums(mi.tt),2)
round(100*rowSums(mi.tt[,1:3])/sum(rowSums(mi.tt[,1:3])),2)
round(100*rowSums(mi.tt[,6:7])/sum(rowSums(mi.tt[,6:7])),2)



###################################################
### chunk number 18: mapaUM
###################################################
#line 494 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

slc.avs <- info.NM$CNEB[match(rownames(aves.NM),info.NM$Nombre)]
slc.scrb <- info.NM$CNEB[match(rownames(scrb.NM),info.NM$Nombre)]
slc.mrps <- info.NM$CNEB[match(rownames(mrps.NM),info.NM$Nombre)]
um <- unique(c(as.character(CNEB.nm$cdg[CNEB.nm$UM==1]),c(slc.scrb,slc.avs,slc.mrps)))

par(mar=c(1,0,2,0))

##Cuadrícula con los códigos en los bordes
plot(vzla,col=NA,border=NA,ylim=c(0,13))
plot(CNEB[CNEB@data$cdg %in% VBG$cdg[VBG$vzla>0],],border="grey77",add=T)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)

## títulos
title(main="Universo Muestral",line=-1)
mi.ds <- c(0,30,40,20,35,50)
mi.ag <- c(0,45,135,135,45,135)
## Datos del universo muestral celdas sombreadas de la cuadrícula
##plot(CNEB[CNEB@data$cdg %in% um,],add=T,col="grey83",border="grey77")
br <- CNEB.nm$bioreg[order(CNEB.nm$row,CNEB.nm$col)]
br[!(CNEB.nm$cdg[order(CNEB.nm$row,CNEB.nm$col)] %in% um)] <- "Otro"
plot(CNEB,density=mi.ds[br],
     border=c(0,"grey77","grey77","grey77","grey77","grey77")[br],
     angle=mi.ag[br],col="grey77",add=T)

## División política
plot(vzla,border="maroon",add=T)
plot(vzla[vzla@data$ID==36,],add=T,border="maroon",col="grey87")

xs <- seq(-72,-60,by=2)
ys <- rep(0.5,length(xs))
points(xs,ys,pch=3,cex=.7)
text(xs,ys-.3,paste(abs(xs),"º00' W",sep=""),cex=.7)

ys <- seq(1,11,by=2)
xs <- rep(-58,length(ys))
points(xs,ys,pch=3,cex=.7)
text(xs+.3,ys,paste(abs(ys),"º00' N",sep=""),cex=.7,srt=90)

legend(-73,5,levels(br)[-1],density=mi.ds[-1],
     angle=mi.ag[-1],bty="n",cex=1.3)



###################################################
### chunk number 19: mapaBR1 eval=FALSE
###################################################
## #line 540 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
## par(mar=c(1,0,2,0))
## 
## ##Cuadrícula con los códigos en los bordes
## br <- CNEB.nm$bioreg[order(CNEB.nm$row,CNEB.nm$col)]
## plot(CNEB,density=c(0,10,15,20,25,30)[br],
##      border=c(0,"grey77","grey77","grey77","grey77","grey77")[br],
##      angle=c(0,45,30,135,60,45)[br],col="grey77")
## ## División política
## plot(vzla[vzla@data$ID!=36,],border=1,add=T)
## 


###################################################
### chunk number 20: mapaBR2
###################################################
#line 554 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
par(mar=c(1,0,2,0))

##Cuadrícula con los códigos en los bordes
br <- CNEB.nm$bioreg[order(CNEB.nm$row,CNEB.nm$col)]
mi.col <- c(0,"thistle","palegoldenrod","peachpuff3","whitesmoke","wheat")
plot(CNEB,col=mi.col[br],
     border=mi.col[br])
## División política
plot(vzla[vzla@data$ID!=36,],border="grey57",add=T)

tds <- data.frame(scrb=CNEB@data$cdg %in% slc.scrb,mrps=CNEB@data$cdg %in% slc.mrps,aves=CNEB@data$cdg %in% slc.avs)

stars(tds[rowSums(tds)>0,],draw.segments=T,locations=coordinates(CNEB[rowSums(tds)>0,]),labels="",add=T,col.segments=NA,len=.25,lwd=2,key.labels=colnames(tds),key.loc=c(-72,1))

legend(-73,5,levels(br)[-1],fill=mi.col[-1],bty="n")

##symbols(coordinates(CNEB[rowSums(tds)>0,]),stars=as.matrix(tds+0)[rowSums(tds)>0,],inches=.12,add=T)

##x.mrps <- coordinates(CNEB[tds$mrps,])
##x.mrps[,2] <- x.mrps[,2] + .1
##x.scrb <- coordinates(CNEB[tds$scrb,])
##x.scrb[,2] <- x.scrb[,2] - .1
##x.scrb[,1] <- x.scrb[,1] + .1
##x.aves <- coordinates(CNEB[tds$aves,])
##x.aves[,2] <- x.aves[,2] - .1
##x.aves[,1] <- x.aves[,1] - .1

##invisible(text(x.aves, labels=rep("A",nrow(x.aves)), cex=0.7,col=1))
##invisible(text(x.mrps, labels=rep("P",nrow(x.mrps)), cex=0.7,col=1))
##invisible(text(x.scrb, labels=rep("S",nrow(x.scrb)), cex=0.7,col=1))



###################################################
### chunk number 21: 
###################################################
#line 589 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
require(sp)

lit.scrb <- rownames(scrb.LIT)
tds.scrb <- info.NM$CNEB[match(as.numeric(unique(trmp.NM$NM[!is.na(trmp.NM$NM)])),info.NM$NM)]
tds.mrps <- info.NM$CNEB[match(as.numeric(unique(tvn.NM$NM[!is.na(tvn.NM$NM)])),info.NM$NM)]
msm.mrps <- rownames(mrps.MSM)

##pdf(file=paste("~/NeoMapas/img/Figs/",mi.dir,"/MAPA_ColectasAves.pdf",sep=""), width=10, height=8)
par(mar=c(0,0,3,0))
##plot(CNEB,border="grey77")
##title(main="Aves 2010")
##plot(CNEB[CNEB@data$cdg %in% slc.avs,],add=T,col="slateblue3")
##plot(vzla,border="maroon",add=T)
##text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
##text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)

##dev.off()

### abrimos un pdf para hacer el mapa
##pdf(file=paste("~/NeoMapas/img/Figs/",mi.dir,"/MAPA_ColectasMariposas.pdf",sep=""), width=10, height=8)
## cerramos el archivo
##dev.off()

##pdf(file=paste("~/NeoMapas/img/Figs/",mi.dir,"/MAPA_ColectasEscarabajos.pdf",sep=""), width=10, height=8)
##dev.off()



###################################################
### chunk number 22: mapas
###################################################
#line 620 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
layout(matrix(1:4,ncol=2,byrow=T))

## colocamos los márgenes
par(mar=c(1,0,2,0))

## Mariposas 
##Cuadrícula con los códigos en los bordes
plot(CNEB,border=NA)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)

## títulos
title(main="Mariposas",sub="NeoMapas 2005 a 2010 y Museos",line=-1)


## Datos de los museos como celdas sombreadas de la cuadrícula
plot(CNEB[CNEB@data$cdg %in% msm.mrps,],add=T,col="grey83",border="grey77")

## División política
plot(vzla,border="maroon",add=T)

## mostramos las celdas muestreadas y la selección de celdas utilizadas
## en este caso

##plot(CNEB[CNEB@data$cdg %in% tds.mrps,],add=T,col="grey77")
##plot(CNEB[CNEB@data$cdg %in% slc.mrps,],add=T,border="slateblue3")
symbols(coordinates(CNEB[CNEB@data$cdg %in% tds.mrps,]),circle=rep(1,length(unique(tds.mrps))),inches=.06,add=T)
##invisible(text(coordinates(CNEB[CNEB@data$cdg %in% tds.mrps,]), labels=CNEB@data[CNEB@data$cdg %in% tds.mrps,"cdg"], cex=0.7,col=1))

symbols(coordinates(CNEB[CNEB@data$cdg %in% slc.mrps,]),circle=rep(1,length(slc.mrps)),inches=.06,add=T,fg=1,bg=1)
##invisible(text(coordinates(CNEB[CNEB@data$cdg %in% slc.mrps,]), labels=CNEB@data[CNEB@data$cdg %in% slc.mrps,"cdg"], cex=0.7,col="white"))


## Escarabajos

##Cuadrícula con los códigos en los bordes
plot(CNEB,border=NA)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)
title(main="Escarabajos",sub="NeoMapas 2005 a 2010 y Literatura",line=-1)

## Datos de los museos como celdas sombreadas de la cuadrícula
plot(CNEB[CNEB@data$cdg %in% lit.scrb,],add=T,col="grey83",border="grey77")

## División política
plot(vzla,border="maroon",add=T)

## mostramos las celdas muestreadas y la selección de celdas utilizadas
## en este caso

symbols(coordinates(CNEB[CNEB@data$cdg %in% tds.scrb,]),circle=rep(1,length(unique(tds.scrb))),inches=.06,add=T)
##invisible(text(coordinates(CNEB[CNEB@data$cdg %in% tds.scrb,]), labels=CNEB@data[CNEB@data$cdg %in% tds.scrb,"cdg"], cex=0.7,col=1))

symbols(coordinates(CNEB[CNEB@data$cdg %in% slc.scrb,]),circle=rep(1,length(unique(slc.scrb))),inches=.06,add=T,fg=1,bg=1)
##invisible(text(coordinates(CNEB[CNEB@data$cdg %in% slc.scrb,]), labels=CNEB@data[CNEB@data$cdg %in% slc.scrb,"cdg"], cex=0.7,col="white"))

##Aves
##Cuadrícula con los códigos en los bordes
plot(CNEB,border=NA)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)

## títulos
title(main="Aves",sub="NeoMapas 2010",line=-1)


## División política
plot(vzla,border="maroon",add=T)

## mostramos las celdas muestreadas 
symbols(coordinates(CNEB[CNEB@data$cdg %in% slc.avs,]),circle=rep(1,length(slc.avs)),inches=.06,add=T,fg=1,bg=1)
##invisible(text(coordinates(CNEB[CNEB@data$cdg %in% slc.avs,]), labels=CNEB@data[CNEB@data$cdg %in% slc.avs,"cdg"], cex=0.7,col="white"))





###################################################
### chunk number 23: 
###################################################
#line 705 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
load("~/NeoMapas/Rdata/LIT.rda")
load("~/NeoMapas/Rdata/MSM.rda")


aves.LIT <- t(spp.aves[,grep("NA[1-9]",colnames(spp.aves))])
colnames(aves.LIT) <- spp.aves$Latname

aves.LIT.um <- aves.LIT[!(rownames(aves.LIT) %in% "NA5a"),]


scrb.LIT.um <- scrb.LIT[rownames(scrb.LIT) %in% CNEB.nm$cdg[CNEB.nm$UM==1],]
mrps.MSM.um <- mrps.MSM[rownames(mrps.MSM) %in% CNEB.nm$cdg[CNEB.nm$UM==1],]

luq(scrb.BDV$cdg_ref)
dim(scrb.BDV)
dim(scrb.BDV0)
luq(scrb.BDV$cdg_taxon)
luq(scrb.BDV0$cdg_taxon)
sum(scrb.LIT.um)
sum(mrps.MSM.um)
sum(colSums(scrb.LIT.um)>0)
sum(colSums(scrb.LIT.um)>0)
sum(colSums(mrps.MSM.um)>0)
sum(colSums(mrps.MSM)>0)



###################################################
### chunk number 24: 
###################################################
#line 740 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
## solo las especies con nombre
scrb.tau <- specaccum(scrb.NM[,-grep(".sp",colnames(scrb.NM))],method="exact",conditioned=FALSE,gamma="chao")
scrb.tds.tau <- specaccum(scrb.NM,method="exact",conditioned=FALSE,gamma="chao")
scrb.cor.tau <- specaccum(scrb.NM.c,method="exact",conditioned=FALSE,gamma="chao")
##scrb.lit.tau <- specaccum(scrb.LIT,method="exact",conditioned=FALSE,gamma="chao")
scrb.lit.tau <- specaccum(scrb.LIT.um,method="exact",conditioned=FALSE,gamma="chao")

mrps.tau <- specaccum(mrps.NM,method="exact",conditioned=FALSE,gamma="chao")
mrps.msm.tau <- specaccum(mrps.MSM.um,method="exact",conditioned=FALSE,gamma="chao")

aves.tau <- specaccum(aves.NM,method="exact",conditioned=FALSE,gamma="chao")
aves.lit.tau <- specaccum(aves.LIT,method="exact",conditioned=FALSE,gamma="chao")
  
max.scrb <- length(unique(scrb.BDV$nombre[scrb.BDV$cdg_ref=="scarabnet"]))
max.mrps <- 106
max.aves <- 1382



###################################################
### chunk number 25: 
###################################################
#line 760 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

mrps.NM.alpha <- rowSums(mrps.NM>0)
scrb.NM.alpha <- rowSums(scrb.NM>0)
aves.NM.alpha <- rowSums(aves.NM>0)

mrps.MSM.alpha <- rowSums(mrps.MSM.um>0)
scrb.LIT.alpha <- rowSums(scrb.LIT.um>0)
aves.LIT.alpha <- rowSums(aves.LIT>0)[!(rownames(aves.LIT) %in% "NA5a")]

mrps.cneb <- CNEB.nm[match(rownames(mrps.CN),CNEB.nm$cdg),]
mrps.msm.cneb <- CNEB.nm[match(rownames(mrps.MSM.um),CNEB.nm$cdg),]
scrb.cneb <- CNEB.nm[match(scrb.CN,CNEB.nm$cdg),]
scrb.lit.cneb <- CNEB.nm[match(rownames(scrb.LIT.um),CNEB.nm$cdg),]
aves.cneb <- aves.lit.cneb <- CNEB.nm[match(rownames(aves.CN),CNEB.nm$cdg),]



###################################################
### chunk number 26: 
###################################################
#line 778 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
cols <- c("grey42","grey77")
## realizamos una figura con las curvas para cada grupo
layout(matrix(1:4,ncol=2,byrow=T))
par(mar=c(3,3,2,1))
plot(scrb.lit.tau, main="Scarabaeinae",
     xlab="", ylab="", 
     ylim=c(0,max.scrb*1.10),
     ci.type="polygon", col = cols[2], ci.col = cols[2],border=cols[2])
plot(scrb.lit.tau,add=T,ci=0)
##plot(scrb.tds.tau,add=T,ci=0,lty=2)
plot(scrb.tds.tau,add=T,ci.type="polygon", col = 1, ci.col = cols[1],border=cols[1])
plot(scrb.tau,add=T,ci.type="line", col = 1, ci.col = cols[1],ci.lty=2)
abline(h=max.scrb,lty=3,col=2)
text(60,max.scrb+3,paste("ScarabNet:",max.scrb,"spp."),col=2,cex=.8)

plot(mrps.msm.tau,
     xlab="", ylab="",  
     ylim=c(0,max.mrps*1.10),
     ci.type="polygon", col = cols[2], ci.col = cols[2],border=cols[2],main="Pieridae")
plot(mrps.msm.tau,add=T,ci=0)
plot(mrps.tau,add=T,ci.type="polygon", col = 1, ci.col = cols[1],border=cols[1])
abline(h=max.mrps,lty=3,col=2)
text(60,max.mrps+3,paste("Viloria (1990):",max.mrps,"spp."),col=2,cex=.8)

plot(aves.lit.tau,
     xlab="", ylab="", 
     ylim=c(0,max.aves*1.10),
     ci.type="polygon", col = cols[2], ci.col = cols[2],border=cols[2],
     main="Aves")
plot(aves.lit.tau,add=T,ci=0)
plot(aves.tau,add=T,ci.type="polygon", col = 1, ci.col = cols[1],border=cols[1])

abline(h=max.aves,lty=3,col=2)
text(20,max.aves+30,paste("Hilty (2003):",max.aves,"spp."),col=2,cex=.8)
title(xlab="Número de celdas", ylab="Número de especies", outer=T)

symbols(c(1,2),c(1,1),boxplots=matrix(c(.5,2,0,0,.5,.5,2,0,0,.5),ncol=5,byrow=T),bg=cols,fg=cols,lwd=2,xlim=c(0.7,3),axes=F)
text(c(1.5,2.5),c(1,1),c("NeoMapas","Other\nsources"))


###################################################
### chunk number 27: 
###################################################
#line 819 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
cols <- c("grey42","grey77")
## realizamos una figura con las curvas para cada grupo
layout(matrix(1:4,ncol=2,byrow=T))
par(mar=c(3,3,2,1))


boxplot(mrps.NM.alpha~mrps.cneb$bioreg,at=c(1:6)-.25,width=rep(.1,6),xlim=c(1.6,6.47),varwidth=T,boxwex=.3,ylim=c(0,39),notch=F,axes=F,col = cols[1])
axis(2)
axis(1,1:6,levels(mrps.cneb$bioreg))
##mrps.cneb$bioreg[mrps.NM.alpha>40]
##boxplot(rowSums(mrps.MSM>0)~m2.cneb$bioreg)
boxplot(mrps.MSM.alpha~mrps.msm.cneb$bioreg,at=c(1:6)+.25,add=T,varwidth=T,boxwex=.3,col="grey74",notch=F,axes=F)
box()

text(as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>39])+.25,rep(39,sum(mrps.MSM.alpha>39)),mrps.MSM.alpha[mrps.MSM.alpha>39])
arrows(as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>39])+.25, rep(37,sum(mrps.MSM.alpha>39)), as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>39])+.25, rep(38.6,sum(mrps.MSM.alpha>39)),length = 0.09)


boxplot(scrb.NM.alpha~scrb.cneb$bioreg,at=c(1:6)-.25,width=rep(.1,6),xlim=c(1.6,6.47),varwidth=T,boxwex=.3,ylim=c(0,39),notch=F,axes=F,col = cols[1])
axis(2)
axis(1,1:6,levels(scrb.cneb$bioreg))
##scrb.cneb$bioreg[scrb.NM.alpha>40]
##boxplot(rowSums(mrps.MSM>0)~m2.cneb$bioreg)
boxplot(scrb.LIT.alpha~scrb.lit.cneb$bioreg,at=c(1:6)+.25,add=T,varwidth=T,boxwex=.3,col="grey74",notch=F,axes=F)
box()

text(as.numeric(scrb.lit.cneb$bioreg[scrb.LIT.alpha>39])+.25,rep(39,sum(scrb.LIT.alpha>39)),scrb.LIT.alpha[scrb.LIT.alpha>39])
arrows(as.numeric(scrb.lit.cneb$bioreg[scrb.LIT.alpha>39])+.25, rep(37,sum(scrb.LIT.alpha>39)), as.numeric(scrb.lit.cneb$bioreg[scrb.LIT.alpha>39])+.25, rep(38.6,sum(scrb.LIT.alpha>39)),length = 0.09)

boxplot(aves.NM.alpha~aves.cneb$bioreg,at=c(1:6)-.25,width=rep(.1,6),xlim=c(1.6,6.47),varwidth=T,boxwex=.3,ylim=c(0,540),notch=F,axes=F,col = cols[1])
axis(2)
axis(1,1:6,levels(aves.cneb$bioreg))
boxplot(aves.LIT.alpha~aves.lit.cneb$bioreg,at=c(1:6)+.25,add=T,varwidth=T,boxwex=.3,col=cols[2],notch=F,axes=F)
box()


symbols(c(1,2),c(1,1),boxplots=matrix(c(.5,2,1,1,.5,.5,2,1,1,.5),ncol=5,byrow=T),bg=cols,fg=1,lwd=1,xlim=c(0.7,2.3),ylim=c(.25,1.95),axes=F,lty=1)

points(c(1,2),c(1.85,1.85))
text(c(1.5,1.5,1.5,1.5,1.5,1.5),c(0.35,0.75,1,1.25,1.65, 1.85),c("mínimum within\n1.5 interquantile\nrange","25%","50%","75%","maximum within\n1.5 interquantile\nrange","outlying observation"),cex=.75)

text(c(1,2),c(.29,.29),c("NeoMapas","Other\nsources"))


###################################################
### chunk number 28: 
###################################################
#line 865 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
dts.scrb <- dts.mrps <- dts.aves <- data.frame()
for (br in levels(CNEB.nm$bioreg)[-1]) {
 
  spp.NM <- colnames(scrb.NM)[colSums(scrb.NM[scrb.cneb$bioreg %in% br,])>0]
  tmp <- Escarabajos[Escarabajos$cdg %in% colnames(scrb.LIT.um)[colSums(scrb.LIT.um[scrb.lit.cneb$bioreg %in% br,])>0],c("gen","esp")]
  spp.LIT <- paste(tmp$gen,tmp$esp,sep=".")
 spp.chao <- specpool(scrb.NM.c[scrb.cneb$bioreg %in% br,])$chao
  se.chao <- specpool(scrb.NM.c[scrb.cneb$bioreg %in% br,])$chao.se
 spp2.chao <- specpool(scrb.LIT.um[scrb.lit.cneb$bioreg %in% br,])$chao
  se2.chao <- specpool(scrb.LIT.um[scrb.lit.cneb$bioreg %in% br,])$chao.se
  dts.scrb <- rbind(dts.scrb,
                    data.frame(bioreg=br,
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                               spp.LIT=length(unique(spp.LIT)),
                               chao.LIT=spp2.chao,
                               se.LIT=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))

    spp.NM <- colnames(mrps.NM)[colSums(mrps.NM[mrps.cneb$bioreg %in% br,])>0]
  tmp <- Mariposas[Mariposas$cdg %in% colnames(mrps.MSM.um)[colSums(mrps.MSM.um[mrps.msm.cneb$bioreg %in% br,])>0],c("gen","esp")]
  spp.LIT <- paste(tmp$gen,tmp$esp,sep=" ")
 spp.chao <- specpool(mrps.NM[mrps.cneb$bioreg %in% br,])$chao
  se.chao <- specpool(mrps.NM[mrps.cneb$bioreg %in% br,])$chao.se
 spp2.chao <- specpool(mrps.MSM.um[mrps.msm.cneb$bioreg %in% br,])$chao
  se2.chao <- specpool(mrps.MSM.um[mrps.msm.cneb$bioreg %in% br,])$chao.se

  dts.mrps <- rbind(dts.mrps,
                    data.frame(bioreg=br,
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                               spp.LIT=length(unique(spp.LIT)),
                               chao.LIT=spp2.chao,
                               se.LIT=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))


  spp.NM <- colnames(aves.NM)[colSums(aves.NM[aves.cneb$bioreg %in% br,])>0]
  spp.LIT <- colnames(aves.LIT.um)[colSums(aves.LIT.um[aves.lit.cneb$bioreg %in% br,])>0]
  spp.chao <- specpool(aves.NM[aves.cneb$bioreg %in% br,])$chao
  se.chao <- specpool(aves.NM[aves.cneb$bioreg %in% br,])$chao.se
  spp2.chao <- specpool(aves.LIT.um[aves.cneb$bioreg %in% br,])$chao
  se2.chao <- specpool(aves.LIT.um[aves.cneb$bioreg %in% br,])$chao.se

  dts.aves <- rbind(dts.aves,
                    data.frame(bioreg=br,
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                               spp.LIT=length(unique(spp.LIT)),
                               chao.LIT=spp2.chao,
                               se.LIT=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))

}
spp.NM <- colnames(scrb.NM)
tmp <- Escarabajos[Escarabajos$cdg %in% colnames(scrb.LIT.um),]
spp.LIT <- paste(tmp$gen,tmp$esp,sep=".")
spp.chao <- specpool(scrb.NM)$chao
se.chao <- specpool(scrb.NM)$chao.se
spp2.chao <- specpool(scrb.LIT.um)$chao
se2.chao <- specpool(scrb.LIT.um)$chao.se

  dts.scrb <- rbind(dts.scrb,
                    data.frame(bioreg="total",
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                               spp.LIT=length(unique(spp.LIT)),
                               chao.LIT=spp2.chao,
                               se.LIT=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))

spp.NM <- colnames(mrps.NM)
tmp <- Mariposas[Mariposas$cdg %in% colnames(mrps.MSM.um),]
spp.LIT <- paste(tmp$gen,tmp$esp,sep=" ")
spp.chao <- specpool(mrps.NM)$chao
se.chao <- specpool(mrps.NM)$chao.se
spp2.chao <- specpool(mrps.MSM.um)$chao
se2.chao <- specpool(mrps.MSM.um)$chao.se
dts.mrps <- rbind(dts.mrps,
                  data.frame(bioreg="total",
                             spp.NM=length(unique(spp.NM)),
                             chao.NM=spp.chao,
                             se.NM=se.chao,
                             spp.LIT=length(unique(spp.LIT)),
                             chao.LIT=spp2.chao,
                             se.LIT=se2.chao,
                             total=length(unique(c(spp.NM,spp.LIT)))))

  spp.NM <- colnames(aves.NM)[colSums(aves.NM)>0]
  spp.LIT <- colnames(aves.LIT.um)[colSums(aves.LIT.um)>0]
   spp.chao <- specpool(aves.NM)$chao
  se.chao <- specpool(aves.NM)$chao.se
   spp2.chao <- specpool(aves.LIT.um)$chao
  se2.chao <- specpool(aves.LIT.um)$chao.se
 dts.aves <- rbind(dts.aves,
                    data.frame(bioreg="total",
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                               spp.LIT=length(unique(spp.LIT)),
                               chao.LIT=spp2.chao,
                               se.LIT=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))

xtable(dts.scrb)
xtable(dts.mrps)
xtable(dts.aves)


###################################################
### chunk number 29: 
###################################################
#line 992 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
dts <- c()
mtds <- c("mountford","jaccard","chao")
for (mtd in mtds) {
  mi.d0 <- vegdist(mrps.NM,mtd)
  mi.d1 <- vegdist(scrb.NM,mtd)
  mi.d2 <- vegdist(aves.NM,mtd)
  mi.h0 <- agnes(mi.d0)
  mi.h1 <- agnes(mi.d1)
  mi.h2 <- agnes(mi.d2)
  dts <- c(dts,mi.h0$ac,mi.h1$ac,mi.h2$ac)
}

xtable(matrix(dts,ncol=3,
       dimnames=list(c("Pieridae","Scarabaeinae","Aves"),
         mtds)),
       caption="Coeficiente de aglomeración según un algoritmo de aglomeración jerárquica (agnes) para tres medidas de disimilitud biótica y los tres grupos de estudio de NeoMapas.")



###################################################
### chunk number 30: 
###################################################
#line 1014 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
plot(mi.h0,which.plots=2)


###################################################
### chunk number 31: 
###################################################
#line 1019 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
plot(mi.h1,which.plots=2)


###################################################
### chunk number 32: 
###################################################
#line 1024 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
plot(mi.h2,which.plots=2)


###################################################
### chunk number 33: 
###################################################
#line 1032 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
mrps.cneb <- CNEB.nm[match(rownames(mrps.CN),CNEB.nm$cdg),]
table(mrps.cneb[,c("bioreg","est.pr")])
##adonis(mi.d1~bioreg+est.pr,scrb.cneb)
adonis(mi.d0~bioreg+(PC1+PC2+PC3),mrps.cneb)


###################################################
### chunk number 34: 
###################################################
#line 1038 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
##adonis(mi.d1~bioreg+est.pr,scrb.cneb)
adonis(mi.d1~bioreg+(PC1+PC2+PC3),scrb.cneb)


###################################################
### chunk number 35: 
###################################################
#line 1043 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
table(aves.cneb[,c("bioreg","est.pr")])
##adonis(mi.d1~bioreg+est.pr,scrb.cneb)
adonis(mi.d2~bioreg+(PC1+PC2+PC3),aves.cneb)


###################################################
### chunk number 36: 
###################################################
#line 1051 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

x <- capscale(scrb.NM~PC1+PC2+PC3+bioreg,data=scrb.cneb,
              sqrt.dist = TRUE,
              distance = "chao", dfun=vegdist)
##x <- rda(decostand(scrb.NM,method="log"))
y <- plot(x,scaling=-1,col=2,pch=4,type="none")
text(x,display="cn",scaling=-1,col="blue",cex=.7)
text(x, display = "wa",cex=.5,col="grey56",scaling=-1)
points(x,scaling=-1,display="sp",col=2,pch=3,cex=.6,select=colnames(scrb.NM)[rowSums(indpower(scrb.NM),na.rm=T)<=43])
##identify(y,what="species")
text(x,scaling=-1,display="sp",select=colnames(scrb.NM)[rowSums(indpower(scrb.NM),na.rm=T)>43],cex=.5,col=2)




###################################################
### chunk number 37: 
###################################################
#line 1072 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
mi.tt <- aggregate(jmp.Museos,list(jmp.Museos$especie),function(x){length(unique(x))})
mi.tt


###################################################
### chunk number 38: 
###################################################
#line 1078 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
plot(rowSums(mrps.NM),rowSums(mrps.NM>0),xlab="Nr. ejemplares",ylab="Nr. Especies")


###################################################
### chunk number 39: 
###################################################
#line 1082 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
##plot(rowSums(scrb.NM),rowSums(scrb.NM>0))
plot(log(rowSums(scrb.NM)),log(rowSums(scrb.NM>0)),xlab="Nr. ejemplares",ylab="log(H)",ylim=c(0,5))
points(log(rowSums(scrb.NM)),renyi(scrb.NM,scale=2),col=2)
points(log(rowSums(scrb.NM)),renyi(scrb.NM,scale=4),col=4)


###################################################
### chunk number 40: 
###################################################
#line 1089 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

##plot(scrb.sfrz,log(rowSums(scrb.NM>0)))
plot(rowSums(scrb.NM)/scrb.sfrz,rowSums(scrb.NM>0),col=1+(scrb.sfrz<4500 | scrb.sfrz>5500))



###################################################
### chunk number 41: 
###################################################
#line 1095 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"

##plot(mrps.sfrz,log(rowSums(mrps.NM>0)))
plot(rowSums(mrps.NM)/mrps.sfrz,rowSums(mrps.NM>0),col=1+(mrps.sfrz<35 | mrps.sfrz>65))



###################################################
### chunk number 42: 
###################################################
#line 1102 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
plot(rowSums(aves.NM),rowSums(aves.NM>0))


###################################################
### chunk number 43: 
###################################################
#line 1108 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
abnd.scrb <- apply(log1p(scrb.NM),1,sd)
abnd.aves <- apply(log1p(aves.NM),1,sd)
abnd.mrps <- apply(log1p(mrps.NM),1,sd)

##abnd.scrb <- apply(log1p(scrb.NM),1,mean)
##abnd.aves <- apply(log1p(aves.NM),1,mean)
##abnd.mrps <- apply(log1p(mrps.NM),1,mean)


H0.scrb <- diversity(scrb.NM, index = "shannon")
H0.aves <- diversity(aves.NM, index = "shannon")
H0.mrps <- diversity(mrps.NM, index = "shannon")

H1.scrb <- renyi(scrb.NM,scales=1)
H1.aves <- renyi(aves.NM,scales=1)
H1.mrps <- renyi(mrps.NM,scales=1)

ACE.scrb <- estimateR(scrb.NM)["S.ACE",]
ACE.aves <- estimateR(aves.NM)["S.ACE",]
ACE.mrps <- estimateR(mrps.NM)["S.ACE",]

se.scrb <- 1/estimateR(scrb.NM)["se.ACE",]
se.aves <- 1/estimateR(aves.NM)["se.ACE",]
se.mrps <- 1/estimateR(mrps.NM)["se.ACE",]



###################################################
### chunk number 44: 
###################################################
#line 1136 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
lm.abnd.aves <- lm(abnd.aves~bioreg+PC1+PC2+PC3,data=aves.cneb)
lm.abnd.scrb <- lm(abnd.scrb~bioreg+PC1+PC2+PC3,data=scrb.cneb,weights=scrb.sfrz)
lm.abnd.mrps <- lm(abnd.mrps~bioreg+PC1+PC2+PC3,data=mrps.cneb,weights=mrps.sfrz)

lm.H0.aves <- lm(H0.aves~bioreg+PC1+PC2+PC3,data=aves.cneb)
lm.H0.scrb <- lm(H0.scrb~bioreg+PC1+PC2+PC3,data=scrb.cneb,weights=scrb.sfrz)
lm.H0.mrps <- lm(H0.mrps~bioreg+PC1+PC2+PC3,data=mrps.cneb,weights=mrps.sfrz)

##lm.H1.aves <- lm(H1.aves~bioreg+PC1+PC2+PC3,data=aves.cneb)
##lm.H1.scrb <- lm(H1.scrb~bioreg+PC1+PC2+PC3,data=scrb.cneb)
##lm.H1.mrps <- lm(H1.mrps~bioreg+PC1+PC2+PC3,data=mrps.cneb)

lm.ACE.aves <- lm(ACE.aves~bioreg+PC1+PC2+PC3,data=aves.cneb,weights=se.aves)
lm.ACE.scrb <- lm(ACE.scrb~bioreg+PC1+PC2+PC3,data=scrb.cneb,weights=se.scrb)
lm.ACE.mrps <- lm(ACE.mrps~bioreg+PC1+PC2+PC3,data=mrps.cneb,weights=se.mrps)


###################################################
### chunk number 45: 
###################################################
#line 1154 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
step.abnd.mrps <- step(lm.abnd.mrps)
step.abnd.scrb <- step(lm.abnd.scrb)
step.abnd.aves <- step(lm.abnd.aves)

step.H0.mrps <- step(lm.H0.mrps)
step.H0.scrb <- step(lm.H0.scrb)
step.H0.aves <- step(lm.H0.aves)

##step.H1.mrps <- step(lm.H1.mrps)
##step.H1.scrb <- step(lm.H1.scrb)
##step.H1.aves <- step(lm.H1.aves)

step.ACE.mrps <- step(lm.ACE.mrps)
step.ACE.scrb <- step(lm.ACE.scrb)
step.ACE.aves <- step(lm.ACE.aves)



###################################################
### chunk number 46:  eval=FALSE
###################################################
## #line 1173 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
## summary(lm.abnd.mrps)
## summary(lm.abnd.aves)
## summary(lm.abnd.scrb)
## 
## ##summary(lm.H1.mrps)
## ##summary(lm.H1.aves)
## ##summary(lm.H1.scrb)
## 
## summary(lm.H0.mrps)
## summary(lm.H0.aves)
## summary(lm.H0.scrb)
## 
## summary(lm.ACE.mrps)
## summary(lm.ACE.aves)
## summary(lm.ACE.scrb)
## 


###################################################
### chunk number 47: 
###################################################
#line 1194 "~/NeoMapas/doc/200_InformeNeoMapas2010/Documento1_InformeActividadesNeoMapas.Rnw"
summary(step.abnd.mrps)
summary(step.abnd.aves)
summary(step.abnd.scrb)

##summary(step.H1.mrps)
##summary(step.H1.aves)
##summary(step.H1.scrb)

summary(step.H0.mrps)
summary(step.H0.aves)
summary(step.H0.scrb)

summary(step.ACE.mrps)
summary(step.ACE.aves)
summary(step.ACE.scrb)


