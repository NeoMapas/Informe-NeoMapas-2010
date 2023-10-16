###################################################
### chunk number 1: 
###################################################
#line 59 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
require(foreign)
require(ade4)
require(vegan)
require(cluster)
require(xtable)
require(labdsv)
require(randomForest)
require(gtools)
require(sp)
require(gdata)
require(RColorBrewer)
require(spatstat)
require(ROpenOffice)

paquetes <- (.packages())
paquetes <- paquetes[!(paquetes %in% c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base", "deldir", "DBI", "RMySQL"))]

luq <- function(x,contar.NA=FALSE) {
	if (contar.NA==F) {
	x <- x[!is.na(x)]
	}
 length(unique(x))
 }



###################################################
### chunk number 2: citas paquetes
###################################################
#line 87 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
cat(paste("\\emph{",paquetes,"} \\citep{pqt::",paquetes,"}",sep="",collapse="; "))


###################################################
### chunk number 3: Datos de otras fuentes
###################################################
#line 94 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

##load("~/NeoMapas/Rdata/LIT.rda") ## sustituimos esto por el siguiente: 
load("~/NeoMapas/Rdata/RevisionLiteratura.rda")

load("~/NeoMapas/Rdata/MSM.rda")
##load("~/NeoMapas/Rdata/GBIF.rda")
load("~/NeoMapas/Rdata/GBIFcsv.rda")
## para no ocupar tanto espacio de memoria
rm(list=grep("Plantae",ls(),value=T))
rm(list=grep("Mammalia",ls(),value=T))


if(!exists("lamas")) {
	lamas <- read.ods("~/NeoMapas/etc/1_Coleccion/130_Taxonomia/Lamas2004_2011_06_22_Todo.ods",TRUE,stringsAsFactors=F)
}
if(!exists("tbl.scrb")) {
	tbl.scrb <- read.ods(file="~/NeoMapas/etc/1_Coleccion/130_Taxonomia/20110701_Scarabaeinae.ods",TRUE,stringsAsFactors=F)
}

	Mariposas <- data.frame()
	cols <- colnames(lamas[[2]])
	for (i in names(lamas)[-1]) {
		tt <- lamas[[i]]
		colnames(tt) <- cols
		Mariposas <- rbind(Mariposas, tt)
	}

Escarabajos <- tbl.scrb[["Taxonomia"]]


mrps.BDV <- mrps.cneb[mrps.cneb$adm0 %in% "VEN" & mrps.cneb$validacion %in% c("exacto", "ortografico", "soundex epiteto",  "soundex genero","manual"),]

mrps.BDV$cdg_ref[mrps.BDV$cdg_ref == "Nei 92"] <- "Nei92"
mrps.BDV$cdg_ref[mrps.BDV$cdg_ref == "miza"] <- "MIZA"
mrps.BDV$cdg_ref[mrps.BDV$cdg_ref == "RAC91"] <- "Rac91"
mrps.BDV$cdg_ref[mrps.BDV$cdg_ref == "VIL94a"] <- "Vil94a"
mrps.BDV$cdg_ref[mrps.BDV$cdg_ref == "ViL98"] <- "Vil98"

mrps.BDV0 <- mrps.BDV[!is.na(mrps.BDV$lat) & !is.na(mrps.BDV$lon),]

scrb.BDV <- scrb.cneb[scrb.cneb$adm0 %in% "VEN" & scrb.cneb$validacion %in% c("vld-auto :: concordancia exacta", "heredada", "manual"),]
scrb.BDV0 <- scrb.BDV[!is.na(scrb.BDV$lat) & !is.na(scrb.BDV$lon),]

pird.BDV0 <- mrps.BDV0[grepl("^097",mrps.BDV0$cdg),c("cdg","lon","lat","cneb")]
colnames(pird.BDV0) <- c("especie","lon","lat","cneb")

scrb.otrs <- scrb.BDV0[, c("ecdg","lon","lat","cneb")]
colnames(scrb.otrs) <- c("especie","lon","lat","cneb")

mrps.otrs <- rbind(pird.BDV0,
	jmp.Museos[,c("especie","lon","lat","cneb")])

mrps.otrs$spp <- paste(Mariposas$gen[match(mrps.otrs$especie,Mariposas$ecdg)], Mariposas$esp[match(mrps.otrs$especie,Mariposas$ecdg)], sep=" ")

##mrps.otrs$spp[mrps.otrs$especie %in% c("097-00103","097-00101")] <- "Eurema sp. [daira o elathea]"
mrps.otrs$spp[mrps.otrs$especie %in% c("097-00095","097-00091")] <- "Pyrisitia sp. [nise/venusta]"

mrps.OTR <- table(mrps.otrs$cneb,mrps.otrs$spp)

scrb.otrs$spp <- paste(Escarabajos$genero[match(scrb.otrs$especie,Escarabajos$cdg_especie)], Escarabajos$epiteto[match(scrb.otrs$especie,Escarabajos$cdg_especie)], sep=" ")

scrb.OTR <- table(scrb.otrs$cneb,scrb.otrs$spp)


###################################################
### chunk number 4: Cargar Datos NM
###################################################
#line 160 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## datos de las bases de datos de NeoMapas
load(file="~/NeoMapas/Rdata/NM.rda")

##Definimos tres fases del trabajo de campo
tvn.NM$fase <- NA
tvn.NM$fase[tvn.NM$fecha < "2006-06-30"] <- 1
tvn.NM$fase[tvn.NM$yr %in% c(2009, 2010) & tvn.NM$NM!="97"] <- 3
tvn.NM$fase[tvn.NM$fecha > "2006-06-30" & tvn.NM$yr==2006] <- 2
jmp.NM$fase <- NA
jmp.NM$fase[jmp.NM$tvn %in% tvn.NM$tvn[tvn.NM$fase %in% 1]] <- 1
jmp.NM$fase[jmp.NM$tvn %in% tvn.NM$tvn[tvn.NM$fase %in% 2]] <- 2
jmp.NM$fase[jmp.NM$tvn %in% tvn.NM$tvn[tvn.NM$fase %in% 3]] <- 3
##table(tvn.NM$fase,useNA="always")

mi.tvn <- tvn.NM[tvn.NM$fase %in% 2,]

mi.jmp <- jmp.NM[(jmp.NM$fase %in% 2) &jmp.NM$fam=="Pieridae" & !is.na(jmp.NM$genero) & jmp.NM$especie != "sp.",]

mi.jmp$esp <- paste(mi.jmp$genero,mi.jmp$especie)

##mi.jmp$esp[mi.jmp$esp %in% c("Eurema daira","Eurema elathea")] <- "Eurema sp. [daira o elathea]"

mrps.NM <- tapply(mi.jmp$jmp,list(mi.jmp$NM,mi.jmp$esp),function(x){length(unique(x))})
mrps.NM[is.na(mrps.NM)] <- 0
mrps.NM <- mrps.NM[,colSums(mrps.NM)>0]
mrps.sfrz <- aggregate(mi.tvn$sfrz,by=list(mi.tvn$NM),sum,na.rm=T)
mrps.sfrz <- mrps.sfrz$x[match(rownames(mrps.NM),mrps.sfrz$Group.1)]

mrps.CN <- mrps.NM
rownames(mrps.NM) <- paste("NM",rownames(mrps.NM)," :: ",info.NM$Nombre[match(as.numeric(rownames(mrps.NM)),info.NM$NM)],sep="")
rownames(mrps.CN) <- info.NM$CNEB[match(as.numeric(rownames(mrps.CN)),info.NM$NM)]


###################################################
### chunk number 5: Cargar SIG
###################################################
#line 194 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
load(file="~/NeoMapas/Rdata/SIG.rda")
##VBG <- read.dbf("/var/local/gis/mi.gis/Venezuela/NeoMapas/dbf/VBG.dbf")
VBG <- read.dbf("/Users/jferrer/mi.gis/Venezuela/NeoMapas/dbf/VBG.dbf")
CNEB.nm <- merge(CNEB@data,VBG, by=1:5)
CNEB.nm <- CNEB.nm[match(CNEB@data$cdg,CNEB.nm$cdg),]

##CNEB.nm$bioreg <- factor(CNEB.nm$bioreg)
##CNEB.nm$bioreg <- factor(CNEB.nm$bioreg,labels=c("Otro","Occident","Andean mountains","Coastal mountains","Llanos","Guayana"),levels=c(0,1,5,2,3,4))

CNEB.nm$bioreg <- factor(CNEB.nm$bioreg,labels=c("Otro","Occidente","Andes", "Centro y Costa","Llanos","Guayana"),levels=c(0,1,5,2,3,4))
info.NM[,"Bioregión"] <- CNEB.nm[match(info.NM$CNEB,CNEB.nm$cdg),"bioreg"]




###################################################
### chunk number 6: Cargar Datos Solis
###################################################
#line 211 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

## datos de identificaciones de escarabajos de Solís
tmp <- read.csv("~/NeoMapas/data/AngelSolisDatos/DatosScarabSetiembre2010.csv", sep="\t", header=T, as.is=T,nrow=5)
scrb.solis <- read.csv("~/NeoMapas/data/AngelSolisDatos/DatosScarabSetiembre2010.csv", sep="\t", header=F, as.is=T, skip=2)
colnames(scrb.solis) <- c("NM", "Prd", "Fecha", "Caja", "ADN", colnames(tmp)[6:ncol(tmp)])
##head(tmp)
scrb.solis$yr <- sapply(scrb.solis$Fecha,function(x){strsplit(x,"/")[[1]][3]},simplify=T)
scrb.solis$ttl <- rowSums(scrb.solis[,colnames(tmp)[6:(ncol(tmp)-1)]])
scrb.solis <- scrb.solis[!is.na(scrb.solis$yr) ,]
scrb.solis$NM[scrb.solis$NM=="1000IVIC"] <- "97"
scrb.solis$NM[scrb.solis$NM=="10"] <- "92"
scrb.solis$NM <- sprintf("%02s",scrb.solis$NM)


###################################################
### chunk number 7: Cargar Datos NMAves
###################################################
#line 227 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## datos de observación de aves de Gustavo Rodríguez
NM.m1 <- read.csv("~/NeoMapas/data/NMAves2010/Muestreo1.total.csv")
NM.m2 <- read.csv("~/NeoMapas/data/NMAves2010/Muestreo2.total.csv")
NM.lista <-  read.csv("~/NeoMapas/data/NMAves2010/AvesVistas.total.csv",as.is=T)

spp.aves <- read.csv("~/NeoMapas/data/NMAves2010/Venezuela.total.csv",as.is=T)
trans.info <- read.csv("~/NeoMapas/data/NMAves2010/Transectas.total.csv",as.is=T)
trans.info[,c("NM","CNEB","Bioregión")] <- 
info.NM[match(trans.info$Nombretransecta,info.NM$Nombre),c("NM", "CNEB", "Bioregión")]
trans.info$NM <- sprintf("%02s",trans.info$NM)


mi.ss <- NM.m1$Especieid!=1379
mi.s2 <- NM.m2$Especieid!=1379
##aves.NM <- table(NM.m1$IDTransecta[mi.ss],NM.m1$Especieid[mi.ss])

##aves.NM <- tapply(c(NM.m1$Ntotal[mi.ss],NM.m2$Ntotal[mi.s2]),list(c(NM.m1$IDTransecta[mi.ss],NM.m2$IDTransecta[mi.s2]),c(NM.m1$Especieid[mi.ss],NM.m2$Especieid[mi.s2])),sum,na.rm=T)
aves.NM <- tapply(c(NM.m1$Ntotal[mi.ss]),list(NM.m1$IDTransecta[mi.ss],NM.m1$Especieid[mi.ss]),sum,na.rm=T)
aves.NM[is.na(aves.NM)] <- 0


###################################################
### chunk number 8: Cargar Datos GPS NMAves
###################################################
#line 249 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
pps <- data.frame()
for (ff in dir(pattern=".gdb$","~/NeoMapas/lib/gps/Mapsource_NMAves2010/")) {
system(paste("/Applications/gpsbabel -i gdb -f ~/NeoMapas/lib/gps/Mapsource_NMAves2010/",gsub(" ","\\\\\\ ",ff)," -o csv,prefer_shortnames -F tmp.csv",sep=""))
pps <- rbind(pps,read.csv("tmp.csv",header=F, as.is=T))
}
pps$V3 <- trim(tolower(pps$V3))

muestreos <- pps[grep("^[1-7][a-g]-[0-9]",pps$V3),c("V2","V1","V3")]
muestreos$NM <- substr(muestreos$V3,0,2)
muestreos <- unique(muestreos)

NM.m1$idpunto <- sprintf("%s-%03d",NM.m1$IDTransecta,NM.m1$Punto)
NM.m2$idpunto <- sprintf("%s-%03d",NM.m2$IDTransecta,NM.m2$Punto)

NM.m1[,c("lat","lon")] <-  pps[match(NM.m1$idpunto,pps$V3),c("V1","V2")]
NM.m1[is.na(NM.m1$lat),c("lat","lon")] <-  pps[match(sub("-0","-",NM.m1$idpunto[is.na(NM.m1$lat)]),pps$V3),c("V1","V2")]

NM.m2[,c("lat","lon")] <-  pps[match(NM.m2$idpunto,pps$V3),c("V1","V2")]
NM.m2[is.na(NM.m2$lat),c("lat","lon")] <-  pps[match(sub("-0","-",NM.m2$idpunto[is.na(NM.m2$lat)]),pps$V3),c("V1","V2")]

mi.ss <- NM.m1$Especieid!=1379
mi.s2 <- NM.m2$Especieid!=1379
##aves.NM <- table(NM.m1$IDTransecta[mi.ss],NM.m1$Especieid[mi.ss])

aves.NM <- tapply(c(NM.m1$Ntotal[mi.ss]),list(NM.m1$IDTransecta[mi.ss],NM.m1$Especieid[mi.ss]),sum,na.rm=T)
aves.NM[is.na(aves.NM)] <- 0

colnames(aves.NM) <- spp.aves$Latname[match(colnames(aves.NM),spp.aves$N_ave)]
aves.CN <- aves.NM
rownames(aves.NM) <- trans.info[match(rownames(aves.NM),trans.info$IDTransecta),"Nombretransecta"]

rownames(aves.CN) <- info.NM$CNEB[match(rownames(aves.NM),info.NM$Nombre)]

mi.ss <- NM.m2$Especieid!=1379
m2 <- table(NM.m2$IDTransecta[mi.ss],NM.m2$Especieid[mi.ss])


###################################################
### chunk number 9: taxonomia aves
###################################################
#line 288 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
load(file="~/CEBA/Rdata/taxonomia/aves.rda")

spps <- spp.aves[,c("Latname","Name")]

spps$status <- "no determinado"
spps$altname <- spp.aves$Latname

mi.ss <- iocn$Species..Scientific. != ""

spps$status[tolower(spps$Name) %in%  tolower(iconv(iocn$Species..English.[mi.ss],"latin1","utf8"))] <- "nombre en ingles"

spps$status[sub("us$","a",spps$Latname) %in%  iconv(iocn$spp[mi.ss],"latin1","utf8")] <- "correccion ortografica"
spps$status[sub("a$","us",spps$Latname) %in%  iconv(iocn$spp[mi.ss],"latin1","utf8")] <- "correccion ortografica"

spps$status[spps$Latname %in%  iconv(iocn$spp[mi.ss],"latin1","utf8")] <- "valido"


spps$altname[spps$status=="nombre en ingles"] <- iocn[match(spps$Name[spps$status=="nombre en ingles"], iocn$Species..English.),"spp"]

spps$altname[spps$status=="correccion ortografica"] <- iocn[match(sub("a$","us",spps$Latname)[spps$status=="correccion ortografica"], iocn$spp),"spp"]

spps$altname[spps$status=="correccion ortografica" & is.na(spps$altname)] <- iocn[match(sub("us$","a",spps$Latname)[spps$status=="correccion ortografica" & is.na(spps$altname)], iocn$spp),"spp"]


table(spps$Latname==spps$altname,spps$status)

for (j in spps[spps$status=="no determinado","Latname"]) {
guess <- iocn[grep(strsplit(j," ")[[1]][2],iocn$spp),]
if (nrow(guess)==1) {
	cat(paste(j,"==>",guess$spp,"\n"))
	spps[spps$Latname==j,c("status","altname")] <- c("auto: match epiteto",guess$spp)
	} else {
	cat(paste(j," tiene ",nrow(guess), " ", paste(guess$spp, collapse="; ")," matches\n"))
	}
}



spps[spps$Latname=="Caprimulgus carolinensis",c("status","altname")] <- c("cambio manual","Antrostomus carolinensis")
spps[spps$Latname=="Icterus chrysocephalus",c("status","altname")] <- c("cambio manual","Icterus cayanensis")
 
spps[spps$Latname=="Otus roraimae",c("status","altname")] <- c("cambio manual","Megascops roraimae")
spps[spps$Latname=="Otus ingens",c("status","altname")] <- c("cambio manual","Megascops ingens")
spps[spps$Latname=="Otus albogularis",c("status","altname")] <- c("cambio manual","Megascops albogularis")
spps[spps$Latname=="Lophotriccus pilaris",c("status","altname")] <- c("cambio manual","Atalotriccus pilaris")

spps[spps$Latname=="Coeligena eos",c("status","altname")] <- c("cambio manual","Coeligena bonapartei")
spps[spps$Latname=="Amazilia cupreicauda",c("status","altname")] <- c("cambio manual","Amazilia viridigaster")

spps[spps$Latname=="Heliangelus clarisse",c("status","altname")] <- c("cambio manual","Heliangelus amethysticollis")

spps[spps$Latname=="Buarremon torquatus",c("status","altname")] <- c("cambio manual","Arremon torquatus")
spps[spps$Latname=="Buarremon brunneinuchus",c("status","altname")] <- c("cambio manual","Arremon brunneinucha")
spps[spps$Latname=="Ramphotrigon megacephala",c("status","altname")] <- c("cambio manual","Ramphotrigon megacephalum")
spps[spps$Latname=="Pyrrhura caeruleiceps",c("status","altname")] <- c("cambio manual","Pyrrhura picta")
spps[spps$Latname=="Dysithamnus tucuyensis",c("status","altname")] <- c("cambio manual","Dysithamnus leucostictus")
spps[spps$Latname=="Basileuterus roraimae",c("status","altname")] <- c("cambio manual","Basileuterus bivittatus")

table(spps$Latname==spps$altname,spps$status)

table(iocn[iocn$spp %in% spps[match(spp.aves$Latname,spps$Latname),"altname"],"Order"])
table(iocn[iocn$spp %in% spps[match(colnames(aves.NM),spps$Latname),"altname"],"Order"])

Aves.iocn <- merge(spps,iocn,by.x="altname",by.y="spp")

colnames(aves.NM) <- tolower(spps[match(colnames(aves.NM),spps$Latname),"altname"])



###################################################
### chunk number 10: 
###################################################
#line 401 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"


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


###################################################
### chunk number 11: 
###################################################
#line 472 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

ttl.loc <- nrow(unique(round(rbind(tvn.NM[,c("lon","lat")], trmp.NM[,c("lon","lat")], avs.NM[,c("lon","lat")], NM.m1[,c("lon","lat")], NM.m2[,c("lon","lat")]),3)))

ttl.spp <- luq(NM.m1$Especieid)+luq(jmp.NM$especie)+luq(paste(scr.NM$genero,scr.NM$especie))



###################################################
### chunk number 12: mapaResumen
###################################################
#line 483 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
par(mar=c(1,0,2,0))

slc.avs <- c(info.NM$CNEB[match(rownames(aves.NM),info.NM$Nombre)],info.NM$CNEB[match(as.numeric(unique(avs.NM$NM)),info.NM$NM)])

slc.scrb <- info.NM$CNEB[match(as.numeric(unique(trmp.NM$NM)),info.NM$NM)]
slc.scrb <- slc.scrb[!is.na(slc.scrb)]

slc.mrps <- info.NM$CNEB[match(as.numeric(unique(tvn.NM$NM)),info.NM$NM)]
slc.mrps <- slc.mrps[!is.na(slc.mrps)]

br <- CNEB.nm$bioreg[order(CNEB.nm$row,CNEB.nm$col)]
mi.col <- c(0,"thistle","palegoldenrod","peachpuff3","whitesmoke","wheat")
plot(CNEB,col=mi.col[br],
     border=mi.col[br])

plot(vzla[!(vzla@data$ID %in% 36),],border="grey57",add=T)
plot(vzla[vzla@data$ID %in% 36,],add=T,border="grey57",col="grey87")

tds <- data.frame(scrb=CNEB@data$cdg %in% slc.scrb,mrps=CNEB@data$cdg %in% slc.mrps,aves=CNEB@data$cdg %in% slc.avs)

stars(tds[rowSums(tds)>0,],draw.segments=T,locations=coordinates(CNEB[rowSums(tds)>0,]),labels="",add=T,col.segments=NA,len=.25,lwd=2,key.labels=colnames(tds),key.loc=c(-72,1))
text(coordinates(vzla),sub("Dependencias Federales","",iconv(as.character(vzla@data$ESTADO),"latin1","utf8")),cex=.5)

legend(-73,5,levels(br)[-1],fill=mi.col[-1],bty="n")


###################################################
### chunk number 13: mapaUM1
###################################################
#line 557 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
um <- as.character(CNEB.nm$cdg[CNEB.nm$UM==1])
par(mar=c(1,0,1,0))

##Cuadrícula con los códigos en los bordes
plot(vzla,col=NA,border=NA,ylim=c(0,13))

text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)
## División política
plot(CNEB[CNEB@data$cdg %in% um,],border="grey77",col="grey77",add=T)
plot(vzla,border="maroon",add=T)
plot(vzla[vzla@data$ID %in% 36,],add=T,border="maroon",col="grey87")

plot(CNEB[CNEB@data$cdg %in% VBG$cdg[VBG$vzla>0],],border="grey47",add=T)
text(coordinates(vzla),sub("NA","",sub("Dependencias Federales","",iconv(as.character(vzla@data$ESTADO),"latin1","utf8"))),cex=.5)

xs <- seq(-72,-60,by=2)
ys <- rep(0.5,length(xs))
points(xs,ys,pch=3,cex=.7)
text(xs,ys-.3,paste(abs(xs),"º00' W",sep=""),cex=.7)

ys <- seq(1,11,by=2)
xs <- rep(-58,length(ys))
points(xs,ys,pch=3,cex=.7)
text(xs+.3,ys,paste(abs(ys),"º00' N",sep=""),cex=.7,srt=90)

## agregar mapa de ubicación continental [bajar resolución...]
##op <- par(fig=c(0.03,.3,0.15,.45),new=TRUE,xpd=F,mar=c(0.5,2,0.5,2))
##plot(amgnsn,border="grey73", xlim=c(-4500004,4000004),ylim=c(-6568465,8264709))
##abline(h=0,lty=3,col=2)
##plot(amgnsn[amgnsn@data$COUNTRY %in% "VENEZUELA",],col=1,add=T)
##box()
##par(op)



###################################################
### chunk number 14: tauUM
###################################################
#line 605 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tmp1 <- CNEB.mtz.Aves[rownames(CNEB.mtz.Aves) %in% CNEB.nm$cdg[CNEB.nm$UM==1],]
tmp1 <- tmp1[,colSums(tmp1)>0]
tmp2 <- CNEB.mtz.Aves[rownames(CNEB.mtz.Aves) %in% CNEB.nm$cdg[CNEB.nm$VEN==1],]
tmp2 <- tmp2[,colSums(tmp2)>0]

tau.avs.um <- specaccum(tmp1,method="exact",conditioned=FALSE,gamma="chao")
tau.avs.ve <- specaccum(tmp2,method="exact",conditioned=FALSE,gamma="chao")

tmp1 <- scrb.OTR[rownames(scrb.OTR) %in% CNEB.nm$cdg[CNEB.nm$UM==1],]
tmp1 <- tmp1[,colSums(tmp1)>0]
tmp2 <- scrb.OTR[rownames(scrb.OTR) %in% CNEB.nm$cdg[CNEB.nm$VEN==1],]
tmp2 <- tmp2[,colSums(tmp2)>0]

tau.scrb.um <- specaccum(tmp1,method="exact",conditioned=FALSE,gamma="chao")
tau.scrb.ve <- specaccum(tmp2,method="exact",conditioned=FALSE,gamma="chao")

tmp1 <- mrps.OTR[rownames(mrps.OTR) %in% CNEB.nm$cdg[CNEB.nm$UM==1],]
tmp1 <- tmp1[,colSums(tmp1)>0]
tmp2 <- mrps.OTR[rownames(mrps.OTR) %in% CNEB.nm$cdg[CNEB.nm$VEN==1],]
tmp2 <- tmp2[,colSums(tmp2)>0]


tau.mrps.um <- specaccum(tmp1,method="exact",conditioned=FALSE,gamma="chao")
tau.mrps.ve <- specaccum(tmp2,method="exact",conditioned=FALSE,gamma="chao")

layout(matrix(1:4,ncol=2))
par(oma=c(0,0,0,0))
 plot(tau.avs.ve, main="Aves")
plot(tau.avs.um,add=T,col=2)
 plot(tau.scrb.ve,main="Scarabaeinae")
plot(tau.scrb.um,add=T,col=2)
 plot(tau.mrps.ve,main="Pieridae")
plot(tau.mrps.um,add=T,col=2)


###################################################
### chunk number 15: 
###################################################
#line 651 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## variables que faltaban
rngs <- CNEB.nm[,grep("max",colnames(CNEB.nm),value=T)]-CNEB.nm[,grep("min",colnames(CNEB.nm),value=T)]

colnames(rngs) <- sub("max","rng",colnames(rngs))
CNEB.nm <- cbind(CNEB.nm,rngs)

vars.nm <- c("lon", "lat", "altr_avrg", "altr_rng", "tmpm_avrg", "tmpm_rng", "prcp_avrg", "prcp_rng", "mscs_avrg", "mscs_rng", "bosq_avrg", "bosq_rng", "decd_avrg", "decd_rng")

mi.cneb <- CNEB.nm[CNEB.nm$UM==1,vars.nm]
rownames(mi.cneb) <- CNEB.nm[CNEB.nm$UM==1,"cdg"]
VIFs.a <- c()
for (mv in vars.nm) {
	VIFs.a <- c(VIFs.a,(sqrt(1/(1-summary(step(lm(formula(paste(mv,"~",".")),mi.cneb)))$adj.r.squared))))
}
names(VIFs.a) <- vars.nm

vars.nm <- c("altr_avrg", "prcp_rng","prcp_avrg", "tmpm_rng", "mscs_avrg", "mscs_rng", "bosq_avrg", "bosq_rng", "decd_avrg", "decd_rng")

mi.cneb <- CNEB.nm[CNEB.nm$UM==1,vars.nm]
rownames(mi.cneb) <- CNEB.nm[CNEB.nm$UM==1,"cdg"]
VIFs <- c()
for (mv in vars.nm) {
	VIFs <- c(VIFs,(sqrt(1/(1-summary(step(lm(formula(paste(mv,"~",".")),mi.cneb)))$adj.r.squared))))
}
names(VIFs) <- vars.nm



###################################################
### chunk number 16: VIFs
###################################################
#line 680 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tmp <- merge(VIFs.a, VIFs, by="row.names", all=T)
tmp$var.nombre <- gsub("lat", "Latitud", gsub("lon", "Longitud", gsub("prcp", "Precipitación total anual", gsub("tmpm", "Temperatura media anual", gsub("altr", "Altitud", gsub("decd", "Cobertura bosque deciduos", gsub("mscs", "Número de meses secos", gsub("bosq", "Cobertura boscosa total", gsub("_", " ", gsub("rng", "(intervalo)", gsub("avrg", "(promedio)", tmp[,1])))))))))))

colnames(tmp) <- c("Codigo.variable","VIF.previo","VIF.seleccion","Variable")
tab <- xtable(tmp[c(7,8,1:6,9:14),c(4,1:3)],caption="Factor de inflación de la varianza (VIF) para las variables consideradas en el diseño muestral",label="TAB:VIF")

print(tab,include.rownames=F, caption.placement="top")


###################################################
### chunk number 17: 
###################################################
#line 690 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
env.pc <- rda(mi.cneb,scale=T)
env.pc2 <- dudi.pca(mi.cneb,nf=10,scale=T, scannf=F)

env.sco <- scores(env.pc,choice=1:3,display="wa")

## hay que hacer esto para incluir celdas que inicialmente no fueron consideradas en el universo muestral
env.sco.full <- predict(env.pc,CNEB.nm,type="wa",scaling=2)

##env.sco.full[CNEB.nm$UM==1,1]

env.str <- env.sco 
env.str[,1] <- env.sco[,1]>median(env.sco[,1])
env.str[,2] <- env.sco[,2]>median(env.sco[,2])
env.str[,3] <- env.sco[,3]>median(env.sco[,3])

##env.dudi <- dudi.pca(mi.cneb,scannf=F,nf=3)
##CNEB.nm$PC1 <- CNEB.nm$PC2 <- CNEB.nm$PC3 <- numeric(nrow(CNEB.nm))
##CNEB.nm[CNEB.nm$UM==1,c("PC1","PC2","PC3")] <- scores(env.pc, choices=1:3,display="sites")

for (i in 1:10) {
CNEB.nm[,paste("PC",i,sep="")] <- numeric(nrow(CNEB.nm))
}
##CNEB.nm[CNEB.nm$UM==1,paste("PC",1:10,sep="")] <- scores(env.pc, choices=1:10,display="sites")

CNEB.nm[,paste("PC",1:10,sep="")] <- env.sco.full



###################################################
### chunk number 18: 
###################################################
#line 722 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

## PC1 altura precipitacion y temperatura
## PC2 cobertura boscosa
## PC3 meses secos
print(xtable(t(cor(scores(env.pc,choice=1:3,display="wa"),mi.cneb)),caption="Correlación entre las variables originales y los scores de los componentes principales",label="TAB:COVPC"), caption.placement="top")



###################################################
### chunk number 19: pc2
###################################################
#line 736 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
d.cp1 <- density(env.sco[,2])
q.cp1 <- quantile(env.sco[,2],probs=seq(0,1,1/2))
d.cp1.1 <- density(env.sco[,2],from=q.cp1[1],to=q.cp1[2])
d.cp1.2 <- density(env.sco[,2],from=q.cp1[2],to=q.cp1[3])

plot(d.cp1,type="n",main="",xlab="Valor en el primer componente principal",ylab="Densidad de observaciones")

polygon(c(q.cp1[1],d.cp1.1$x,q.cp1[2]),c(0,d.cp1.1$y,0),col="lightskyblue",border="lightskyblue")
polygon(c(q.cp1[2],d.cp1.2$x,q.cp1[3]),c(0,d.cp1.2$y,0),col="lightblue",border="lightblue")
polygon(d.cp1)
abline(v=q.cp1,lty=3)
rug(env.sco[,2])


###################################################
### chunk number 20: 
###################################################
#line 759 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## Utilizamos \emph{random forests} para llenar la cuadrícula con una predicción de los estratos (así removemos ceros).
mi.rf <- randomForest(y=as.factor(CNEB.nm$est[CNEB.nm$UM==1 & CNEB.nm$est>0]),x=CNEB.nm[CNEB.nm$UM==1  & CNEB.nm$est>0,15:47])
CNEB.nm$est.pr <- predict(mi.rf,CNEB.nm)
CNEB.nm$est.pr[CNEB.nm$UM==1 & CNEB.nm$est>0] <- CNEB.nm$est[CNEB.nm$UM==1 & CNEB.nm$est>0]
CNEB.nm$est.pr <- factor(CNEB.nm$est.pr)

info.NM[,"Estrato"] <- CNEB.nm[match(info.NM$CNEB,CNEB.nm$cdg),"est.pr"]



###################################################
### chunk number 21: mapaEST
###################################################
#line 774 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

layout(matrix(1:2,ncol=1))
par(mar=c(.5,.5,.5,.5))

um <- as.character(CNEB.nm$cdg[CNEB.nm$UM==1])

par(mar=c(1,0,1,0))
##Cuadrícula con los códigos en los bordes
plot(vzla,col=NA,border=NA,ylim=c(0,13))
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)

cols <- brewer.pal(8,"Accent")
plot(CNEB[CNEB@data$cdg %in% um,],border=cols[CNEB.nm$est.pr[CNEB.nm$cdg %in% um]], col=cols[CNEB.nm$est.pr[CNEB.nm$cdg %in% um]],add=T)


##división política
title(main="(a)",line=-1)

plot(vzla,border="maroon",add=T)
plot(vzla[vzla@data$ID %in% 36,],add=T,border="maroon",col="grey87")

plot(CNEB[CNEB@data$cdg %in% VBG$cdg[VBG$vzla>0],],border="grey47",add=T)
legend(-73,5, paste("Estrato",1:8),fill=cols,cex=.7)


slc.avs <- info.NM$CNEB[match(rownames(aves.NM),info.NM$Nombre)]
slc.scrb <- info.NM$CNEB[match(as.numeric(unique(trmp.NM$NM[!is.na(trmp.NM$NM)])),info.NM$NM)]
slc.mrps <- info.NM$CNEB[match(rownames(mrps.NM),info.NM$Nombre)]
um <- unique(c(as.character(CNEB.nm$cdg[CNEB.nm$UM==1]),c(slc.scrb,slc.avs,slc.mrps)))

par(mar=c(1,0,2,0))
mi.ds <- c(0,30,40,20,35,70)
mi.ag <- c(0,45,135,135,45,135)
mi.cl <- c(0,"grey36","grey56","grey77","grey56","grey77")

plot(vzla,col=NA,border=NA,ylim=c(0,13))
plot(CNEB[CNEB@data$cdg %in% VBG$cdg[VBG$vzla>0],],border="grey77",add=T)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)

title(main="(b)",line=-1)

br <- CNEB.nm$bioreg[order(CNEB.nm$row,CNEB.nm$col)]
br[!(CNEB.nm$cdg[order(CNEB.nm$row,CNEB.nm$col)] %in% um)] <- "Otro"
plot(CNEB,density=mi.ds[br],
     border=mi.cl[br],
     angle=mi.ag[br],col=mi.cl[br],add=T)

## División política
plot(vzla,border="maroon",add=T)
plot(vzla[vzla@data$ID %in% 36,],add=T,border="maroon",col="grey87")

legend(-73,5,levels(br)[-1],density=mi.ds[-1],
     angle=mi.ag[-1],bty="n",cex=.7,col=mi.cl[-1],border=mi.cl[-1])



###################################################
### chunk number 22: 
###################################################
#line 840 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
 xtable(rbind(env.pc$CA$v[,1:3],env.pc$CA$eig[1:3]))


###################################################
### chunk number 23: UMenPCs
###################################################
#line 847 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
plot(env.pc, type="n", xlab="PC1",ylab="PC2")
##points(env.sco[,1],env.sco[,2], col=c("blue3","red3")[env.str[,3]+1])

cols <- brewer.pal(8,"Accent")
symbols(env.sco[,1],env.sco[,2],circle=env.sco[,3]+abs(min(env.sco[,3]))+.01,add=T,inches=.1, bg=cols[CNEB.nm$est.pr[match(rownames(env.sco),CNEB.nm$cdg)]], fg=cols[CNEB.nm$est.pr[match(rownames(env.sco),CNEB.nm$cdg)]])
mi.ss <- CNEB.nm$M2006[CNEB.nm$UM==1]==1
symbols(env.sco[mi.ss,1],env.sco[mi.ss,2],circle=env.sco[mi.ss,3]+abs(min(env.sco[,3]))+.01,add=T,inches=.1, fg=1, bg=cols[CNEB.nm$est.pr[match(rownames(env.sco[mi.ss,]),CNEB.nm$cdg)]],lwd=2)


###################################################
### chunk number 24: 
###################################################
#line 872 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## combinaciones de estratos y bioregiones 
##dim(unique(CNEB.nm[CNEB.nm$UM==1,c("bioreg","est.pr")]))
table(CNEB.nm[CNEB.nm$UM==1,c("bioreg","est.pr")])


###################################################
### chunk number 25: 
###################################################
#line 884 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
slc.mrps <- as.numeric(unique(tvn.NM$NM[tvn.NM$yr %in% 2006]))
slc.mrps <- slc.mrps[!is.na(slc.mrps) & slc.mrps!=92]
slc.mrps <- info.NM$CNEB[match(slc.mrps,info.NM$NM)]

table(CNEB.nm[CNEB.nm$cdg %in% slc.mrps,c("bioreg","est.pr")])


###################################################
### chunk number 26: 
###################################################
#line 896 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
slc.scrb <- info.NM$CNEB[match(as.numeric(unique(trmp.NM$NM[trmp.NM$yr %in% c(2006, 2009, 2010)])),info.NM$NM)]
slc.scrb <- slc.scrb[!is.na(slc.scrb)]

table(CNEB.nm[CNEB.nm$cdg %in% slc.scrb,c("bioreg","est.pr")])


###################################################
### chunk number 27: 
###################################################
#line 905 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
slc.avs <- info.NM$CNEB[match(rownames(aves.NM),info.NM$Nombre)]
table(CNEB.nm[CNEB.nm$cdg %in% slc.avs,c("bioreg","est.pr")])


###################################################
### chunk number 28: tabla URA
###################################################
#line 916 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tbl <- info.NM[1:97,c("NM","CNEB","Nombre","ADM1", "Bioregión", "Estrato")]

tbl$NM <- sprintf("NM%02d",tbl$NM)
tbl[89,3] <- "UNET"
colnames(tbl) <- c("NM","CNEB","Nombre","Estado","Bioregión", "Estrato")

print(xtable(tbl, label="TAB:URAs", caption="Lista de transecciones de NeoMapas planificadas. Para una lista de las transecciones activas, y más detalles de cada una, consulte \\url{http://www.neomapas.org}."), size="footnotesize", caption.placement="top", include.rownames=FALSE, floating=FALSE, tabular.environment="longtable", add.to.row=list(pos=list(-1), command="\\hline\\endhead\n\\hline\\endlastfoot\n\\hline\\multicolumn{4}{l}{Continúa en la página siguiente}\\\\hline\\endfoot\n"))


###################################################
### chunk number 29: 
###################################################
#line 927 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
table(CNEB.nm[CNEB.nm$cdg %in% tbl$CNEB,c("bioreg","est.pr")])


###################################################
### chunk number 30: 
###################################################
#line 939 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
xys <- rbind(tvn.NM[,c("lon","lat")], trmp.NM[,c("lon","lat")], avs.NM[,c("lon","lat")], NM.m1[,c("lon","lat")], NM.m2[,c("lon","lat")])

xys <- unique(xys[!is.na(xys$lon) & !is.na(xys$lat) & xys$lon < 0 & xys$lon > -74 & xys$lat >2,]) 


###################################################
### chunk number 31: 
###################################################
#line 951 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.ppp <- ppp(xys$lon,xys$lat,window=convexhull.xy(xys))
Z <- density(mi.ppp, sigma=.75)
mi.rnd <- rpoint(n=nrow(xys), f=median(Z), win=convexhull.xy(xys))
Z0 <- density(mi.rnd, sigma=.75)

ZD <- Z
ZD$v <- Z$v-Z0$v

par(mar=c(0,0,0,1.6))
plot(ZD,main="")
plot(vzla,add=T,border="maroon")
contour(ZD,add=T,col="grey17")


###################################################
### chunk number 32: SesgoPC3
###################################################
#line 974 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tmp <- crossdist(xys$lon,xys$lat,CNEB@data$x,CNEB@data$y)
xys$cneb <-  as.character(CNEB@data$cdg[apply(tmp,1,which.min)])
tmp2 <- table(xys$cneb)

layout(matrix(1:6,ncol=2, byrow=T))
par(mar=c(3,3,2,0))
d1 <- density(CNEB.nm$PC1[CNEB.nm$UM==1])
d2 <- density(CNEB.nm$PC1[match(names(tmp2),CNEB.nm$cdg)], weights=tmp2/sum(tmp2))

plot(d1,log="", main="PC1",ylim=c(0,1))
lines(d2,col=2)

par(mar=c(0,0,0,0))
CNEB.nm$bias.PC1 <- 0
qq <- table(cc <- cumsum((d2$y > d1$y)+0))
for (i in names(qq)[qq>1]) {
	CNEB.nm$bias.PC1[CNEB.nm$PC1 > min(d2$x[cc == i]) & CNEB.nm$PC1 < max(d2$x[cc == i])] <- 1
	}
plot(CNEB,col=c("white","slateblue4")[CNEB.nm$bias.PC1*CNEB.nm$UM+1], border=c(NA,"slateblue4")[CNEB.nm$UM+1])
	plot(vzla,add=T)

d1 <- density(CNEB.nm$PC2[CNEB.nm$UM==1])
d2 <- density(CNEB.nm$PC2[match(names(tmp2),CNEB.nm$cdg)], weights=tmp2/sum(tmp2))

par(mar=c(3,3,2,0))
plot(d1,log="", main="PC2",ylim=c(0,1))
lines(d2,col=2)

par(mar=c(0,0,0,0))
CNEB.nm$bias.PC2 <- 0
qq <- table(cc <- cumsum((d2$y > d1$y)+0))
for (i in names(qq)[qq>1]) {
	CNEB.nm$bias.PC2[CNEB.nm$PC2 > min(d2$x[cc == i]) & CNEB.nm$PC2 < max(d2$x[cc == i])] <- 1
	}
plot(CNEB,col=c("white","slateblue4")[CNEB.nm$bias.PC2*CNEB.nm$UM+1], border=c(NA,"slateblue4")[CNEB.nm$UM+1])
	plot(vzla,add=T)


d1 <- density(CNEB.nm$PC3[CNEB.nm$UM==1])
d2 <- density(CNEB.nm$PC3[match(names(tmp2),CNEB.nm$cdg)], weights=tmp2/sum(tmp2))

par(mar=c(3,3,2,0))
plot(d1,log="", main="PC3",ylim=c(0,1))
lines(d2,col=2)

par(mar=c(0,0,0,0))
CNEB.nm$bias.PC3 <- 0
qq <- table(cc <- cumsum((d2$y > d1$y)+0))
for (i in names(qq)[qq>1]) {
	CNEB.nm$bias.PC3[CNEB.nm$PC3 > min(d2$x[cc == i]) & CNEB.nm$PC3 < max(d2$x[cc == i])] <- 1
	}
plot(CNEB,col=c("white","slateblue4")[CNEB.nm$bias.PC3*CNEB.nm$UM+1], border=c(NA,"slateblue4")[CNEB.nm$UM+1])
	plot(vzla,add=T)



###################################################
### chunk number 33: SesgoPCs
###################################################
#line 1038 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
CNEB.nm$bias.PCs <- 0
for (j in paste("PC",1:10,sep="")) {
d1 <- density(CNEB.nm[CNEB.nm$UM==1,j],)
d2 <- density(CNEB.nm[match(names(tmp2),CNEB.nm$cdg),j], weights=tmp2/sum(tmp2))

qq <- table(cc <- cumsum((d2$y > d1$y)+0))
for (i in names(qq)[qq>1]) {
	CNEB.nm$bias.PCs[CNEB.nm[,j] > min(d2$x[cc == i]) & CNEB.nm[,j] < max(d2$x[cc == i])] <- 
	CNEB.nm$bias.PCs[CNEB.nm[,j] > min(d2$x[cc == i]) & CNEB.nm[,j] < max(d2$x[cc == i])]+1
	}
	}
cols <- c(NA,NA,NA,NA,NA,brewer.pal(max(CNEB.nm$bias.PCs)-4,"YlOrRd"))
plot(CNEB,col=cols[CNEB.nm$bias.PCs*CNEB.nm$UM+1], border=c(NA,"grey54")[CNEB.nm$UM+1])
		plot(vzla,add=T)
##legend(-73,5,...)



###################################################
### chunk number 34: fase 1 mrps
###################################################
#line 1093 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.tvn <- tvn.NM[tvn.NM$fase %in% 1,]
mi.jmp <- jmp.NM[jmp.NM$fase %in% 1,]
n.vst <- luq(mi.tvn$vst)
sfrz <- sum(mi.tvn$sfrz,na.rm=T)
n.jmp <- nrow(mi.jmp)
p.id <- sum(mi.jmp$status_id %in% c("definitiva","preliminar","morfo-especie"))
n.CN <- luq(mi.tvn$CN)
n.NM <- luq(mi.tvn$NM)


###################################################
### chunk number 35: EsfuerzoFase
###################################################
#line 1111 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tt <- tapply(tvn.NM$sfrz,list(tvn.NM$CN,tvn.NM$fase),sum,na.rm=T)
tt <- tt[rowSums(tt,na.rm=T)>100,]
par(mar=c(0,0,0,0))
plot(CNEB[CNEB@data$cdg %in% rownames(tt),],border=NA)
plot(vzla,add=T,border="maroon")

stars((tt>0)+0, draw.segment=T, scale=F, col.segments=c("grey77","pink2","palegreen"), add=T, len=.23, 
key.labels=c("2003-2005","2006","2009-2010"), key.loc=c(-72,4.8),
locations=coordinates(CNEB)[match(rownames(tt),CNEB@data$cdg),], cex=1,labels=NULL,lwd=2)

stars((tt>1000)+0, draw.segment=T, scale=F, col.segments=1:3,add=T, len=.23, key.labels=c("2003-2005","2006","2009-2010"), key.loc=c(-72,5.8), locations=coordinates(CNEB)[match(rownames(tt),CNEB@data$cdg),], cex=1,labels=NULL,lwd=2)


###################################################
### chunk number 36: fase 2 mrps
###################################################
#line 1135 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.tvn <- tvn.NM[tvn.NM$fase %in% 2,]
mi.jmp <- jmp.NM[jmp.NM$fase %in% 2,]
n.vst <- luq(mi.tvn$vst)
sfrz <- sum(mi.tvn$sfrz,na.rm=T)
n.jmp <- nrow(mi.jmp)
p.id <- sum(mi.jmp$status_id %in% c("definitiva","preliminar","morfo-especie"))
n.CN <- luq(mi.tvn$CN)
n.NM <- luq(mi.tvn$NM)


###################################################
### chunk number 37: fase 3 mrps
###################################################
#line 1155 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.tvn <- tvn.NM[tvn.NM$fase %in% 3,]
mi.jmp <- jmp.NM[jmp.NM$fase %in% 3,]
n.vst <- luq(mi.tvn$vst)
sfrz <- sum(mi.tvn$sfrz,na.rm=T)
n.jmp <- nrow(mi.jmp)
p.id <- sum(mi.jmp$status_id %in% c("definitiva","preliminar","morfo-especie"))
n.CN <- luq(mi.tvn$CN)
n.NM <- luq(mi.tvn$NM)

##table(substr(mi.tvn$fecha,1,7),mi.tvn$NM)
##tapply(mi.tvn$vst,mi.tvn$NM,luq)



###################################################
### chunk number 38: EjemplaresEsfuerzo
###################################################
#line 1175 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

t2 <- tapply(jmp.NM$jmp,list(jmp.NM$CN,jmp.NM$fase),luq)
t2 <- t2[rownames(t2) %in% rownames(tt),]
##t2[is.na(t2)] <- 0

plot(t2~tt,log="y", xlab="Esfuerzo de muestreo", ylab="Ejemplares capturados", col=c(rep(1,nrow(t2)),rep(2,nrow(t2)),rep(3,nrow(t2))),cex=1.9)
text(tt,t2,rownames(t2),col=c(rep(1,nrow(t2)),rep(2,nrow(t2)),rep(3,nrow(t2))),cex=.38)

lines(supsmu(tt[!is.na(tt)],t2[!is.na(tt)],bass=1), lwd=2, col="slateblue4")



###################################################
### chunk number 39: 
###################################################
#line 1198 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tt <- tapply(tvn.NM$sfrz,list(tvn.NM$CN,tvn.NM$fase),sum,na.rm=T)
tt <- tt[rowSums(tt,na.rm=T)>100,]
ttl.CN <- nrow(tt)
cmpl.CN <- sum(rowSums(tt>0,na.rm=T)==3)
## esto muestra algo similar a la figura {MAPA:MRPSFASE}
##tt[rowSums(tt>0,na.rm=T)==3,]


###################################################
### chunk number 40: tablaFamiliaFase
###################################################
#line 1213 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
jmp.NM$bioreg <- CNEB.nm[match(jmp.NM$CN,CNEB.nm$cdg),"bioreg"]

mi.tt <- table(jmp.NM$fase,jmp.NM$familia,useNA="ifany")
fams <- c("Hesperiidae", "Papilionidae", "Pieridae", "Lycaenidae", "Riodinidae", "Nymphalidae")
colnames(mi.tt)[is.na(colnames(mi.tt)) | colnames(mi.tt) %in% "por asignar"] <- "Sin Id"
mi.tt <- mi.tt[1:3,c(fams,"Sin Id")]
colnames(mi.tt) <- c(abbreviate(fams),"Sin Id")
rownames(mi.tt) <- paste("Fase",rownames(mi.tt))
## mosaicplot(mi.tt)

mi.tt <- cbind(mi.tt,total=rowSums(mi.tt))
mi.tt <- rbind(mi.tt,total=colSums(mi.tt))
xtab <- xtable(mi.tt,digits=0,caption=paste("Ejemplares capturados por familia en cada una de las fases del trabajo de campo. ", paste(abbreviate(unique(fams)),": ",unique(fams),sep="", collapse="; "),"; Sin Id: ejemplares no revisados o con identificación dudosa.",sep=""), label="TAB:FAMMRPS")

print(xtab, caption.placement="top")
##mi.tt <- table(jmp.NM$bioreg,jmp.NM$familia,useNA="ifany")
##colnames(mi.tt) <- c("No identificado", "Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae",  "Riodinidae")
##mi.tt <- mi.tt[c(2,6,3,4,5,7),]


###################################################
### chunk number 41: MrpsID
###################################################
#line 1239 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
jmp.NM$familia[is.na(jmp.NM$familia)] <- "por asignar"
jmp.NM$familia <- factor(jmp.NM$familia,levels=c("por asignar", "Hesperiidae", "Papilionidae", "Pieridae", "Lycaenidae", "Riodinidae", "Nymphalidae"))
jmp.NM$spp <- paste(jmp.NM$genero,jmp.NM$especie)
jmp.NM$spp[jmp.NM$spp %in% c("Eurema daira","Eurema elathea", "Eurema sp. [daira / elathea]")] <- "Eurema sp. [daira o elathea]"

jmp.NM$valido <- "SI"
jmp.NM$valido[is.na(jmp.NM$genero)] <- "NO"
jmp.NM$valido[is.na(jmp.NM$especie)] <- "NO"
jmp.NM$valido[grep("NA",jmp.NM$spp)] <- "NO"

mi.ss <- jmp.NM$valido == "SI"

tmp <- data.frame(id=tapply(mi.ss,jmp.NM$CN,sum), ttl=tapply(mi.ss,jmp.NM$CN,length))
tmp$p.id <- tmp$id*100/tmp$ttl
xys <- coordinates(CNEB)[match(rownames(tmp),CNEB@data$cdg),]
grps <-cut(tmp$p.id,breaks=c(-1,20,40,60,80,101))

par(mar=c(0,0,0,0))
plot(CNEB[CNEB@data$vzla>0,],border="grey86")
plot(vzla[vzla@data$ESTADO %in% "Zona en Reclamaci\xf3n",],density=29,add=T,col="maroon",border="maroon")
plot(vzla,border="maroon",add=T)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)

symbols(xys[,1],xys[,2],circle=sqrt(tmp$ttl),bg=c("red","orange","yellow","green","blue")[grps],fg=c("red","orange","yellow","green","blue")[grps],inches=.12,add=T)
legend(-73.4, 5, c("0-20%","20-40%","40-60%","60-80%","80-100%"), fill=c("red","orange","yellow","green","blue"),title="Porcentaje de\n identificación",cex=.7,bty="n")
 
 text(-69.5,5,"Número de\nejemplares",font=2,cex=.7)
 symbols(rep(-70,5),seq(2.6,4.4,length=5),circle=sqrt(quantile(tmp$ttl)),inches=.12,add=T)
 text(rep(-69,5),seq(2.6,4.4,length=5),round(quantile(tmp$ttl)),cex=.7)


###################################################
### chunk number 42: FamIDMrps
###################################################
#line 1277 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
##mariposas
grp <- paste(jmp.NM$familia, "::", jmp.NM$subfamilia)
grp[is.na(jmp.NM$subfamilia)] <- paste(jmp.NM$familia[is.na(jmp.NM$subfamilia)], " por asignar a subfamilia")

grp[is.na(jmp.NM$familia)] <- paste(" por asignar a familia")

mi.tt <- table(grp,jmp.NM$status_id %in% c("definitiva","preliminar"))
colnames(mi.tt) <- c("id. pendiente","id. confiable")
xtab <- xtable(mi.tt, caption="Número de ejemplares con identificación confiable e identificación pendiente por familia y subfamilia.", label="TAB:FAMID")
print(xtab,caption.placement="top", size="footnotesize")


###################################################
### chunk number 43:  eval=FALSE
###################################################
## #line 1294 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## tvn.NM$bioreg <- CNEB.nm[match(tvn.NM$CN,CNEB.nm$cdg),"bioreg"]
## tapply(tvn.NM$vst,list(tvn.NM$bioreg,tvn.NM$fase),luq)
## tapply(tvn.NM$sfrz,list(tvn.NM$bioreg,tvn.NM$fase),sum,na.rm=T)


###################################################
### chunk number 44: RADmrps
###################################################
#line 1303 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tt <- table(jmp.NM$fase[jmp.NM$valido=="SI"],jmp.NM$spp[jmp.NM$valido=="SI"])
tt <- tt[,colSums(tt)>0]
tt <- tt[,rev(order(colSums(tt)))]
layout(matrix(1:2,ncol=1))
par(mar=c(4,4,1,2))
plot(1:ncol(tt),colSums(tt,na.rm=T),type="l",log="y",xlab="Especies ordenadas por conteos totales", ylab="Nr. de ejemplares capturados")
matpoints(1:ncol(tt),t(tt),log="y",type="p",pch=1,cex=.5)

par(mar=c(6,4,1,2))
plot(1:ncol(tt),colSums(tt,na.rm=T),type="l",log="y",xlim=c(1,30),axes=F, xlab="", ylab="Nr. de ejemplares capturados")
matpoints(1:ncol(tt),t(tt),log="y",type="p",pch=1,cex=1)
axis(2)
axis(1,1:30,colnames(tt)[1:30],las=2,cex.axis=.5)
box()


###################################################
### chunk number 45: 
###################################################
#line 1331 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## hacer una comparación por fase/celda, de la composición de especies o buscar especies con potencial indicador para los estratos ambientales o las bioregiones...
mi.ss <- jmp.NM$valido=="SI"
jmp.NM$URA <- paste(jmp.NM$CN,jmp.NM$fase,sep=".")
tt <- tapply(jmp.NM$jmp[mi.ss],list(jmp.NM$URA[mi.ss],jmp.NM$spp[mi.ss]),luq)
tt <- data.frame(tt)
tt[is.na(tt)] <- 0

mi.iv <- indval(tt, as.numeric(CNEB.nm[match(substr(rownames(tt),1,3), CNEB@data$cdg), "bioreg"]))

summary(mi.iv, p=0.01)


###################################################
### chunk number 46: MrpsCRE
###################################################
#line 1358 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

mrps.MSM.um <- mrps.OTR[rownames(mrps.OTR) %in% CNEB.nm$cdg[CNEB.nm$UM==1],]
mrps.MSM.um <- mrps.MSM.um[,colSums(mrps.MSM.um)>0]


max.mrps <- 106
mrps.NM.alpha <- rowSums(mrps.NM>0)
mrps.MSM.alpha <- rowSums(mrps.MSM.um>0)
mrps.cneb <- CNEB.nm[match(rownames(mrps.CN),CNEB.nm$cdg),]
mrps.msm.cneb <- CNEB.nm[match(rownames(mrps.MSM.um),CNEB.nm$cdg),]

mrps.tau <- specaccum(mrps.NM,method="exact",conditioned=FALSE,gamma="chao")
mrps.msm.tau <- specaccum(mrps.MSM.um,method="exact",conditioned=FALSE,gamma="chao")

cols <- c("grey42","grey77","pink3")
cols <- c("slateblue4","lightseagreen","pink3")
## realizamos una figura con las curvas para cada grupo
par(mar=c(4,4,2,1))

plot(mrps.msm.tau,
     xlab="", ylab="",  
     ylim=c(0,max.mrps*1.10),
     ci.type="polygon", col = cols[2], ci.col = cols[2],border=cols[2],main="")
plot(mrps.msm.tau,add=T,ci=0)
plot(mrps.tau,add=T,ci.type="polygon", col = 1, ci.col = cols[1],border=cols[1])
abline(h=max.mrps,lty=3,col=2)
text(60,max.mrps+3,paste("Viloria (1990):",max.mrps,"spp."),col=2,cex=.8)
title(xlab="Celdas de la CNEB", ylab="Nr. de especies detectadas")



###################################################
### chunk number 47: mrpsBoxPlot
###################################################
#line 1397 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
cols <- c("grey42","grey77")
cols <- c("slateblue4","lightseagreen","pink3")
## realizamos una figura con las curvas para cada grupo
par(mar=c(4,4,1,1),oma=rep(0,4))

max.y <- 43
boxplot(mrps.NM.alpha~mrps.cneb$bioreg,at=c(1:6)-.25,width=rep(.1,6),xlim=c(1.6,6.47),varwidth=T,boxwex=.3,ylim=c(0,max.y),notch=F,axes=F,col = cols[1])
title(xlab="Bioregión", ylab="Número de especies por celda", line=3)
axis(2)
axis(1,1:6,rep("",6))
axis(1,seq(1,6,by=2),levels(mrps.cneb$bioreg)[seq(1,6,by=2)],line=-.4,lty=0)
axis(1,seq(2,6,by=2),levels(mrps.cneb$bioreg)[seq(2,6,by=2)],line=.4,lty=0)
##mrps.cneb$bioreg[mrps.NM.alpha>40]
##boxplot(rowSums(mrps.OTR>0)~m2.cneb$bioreg)
boxplot(mrps.MSM.alpha~mrps.msm.cneb$bioreg,at=c(1:6)+.25,add=T,varwidth=T,boxwex=.3,col=cols[2],notch=F,axes=F)
box()

text(as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>max.y])+.25,rep(max.y,sum(mrps.MSM.alpha>max.y)),mrps.MSM.alpha[mrps.MSM.alpha>max.y], cex=.7)
arrows(as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>max.y])+.25, rep(max.y-4,sum(mrps.MSM.alpha>max.y)), as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>max.y])+.25, rep(max.y-1.2,sum(mrps.MSM.alpha>max.y)),length = 0.03)




###################################################
### chunk number 48: TablaSppPieridae
###################################################
#line 1427 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
dts.mrps <- data.frame()
for (br in levels(CNEB.nm$bioreg)[-1]) {

    spp.NM <- colnames(mrps.NM)[colSums(mrps.NM[mrps.cneb$bioreg %in% br,])>0]
   spp.LIT <- colnames(mrps.MSM.um)[colSums(mrps.MSM.um[mrps.msm.cneb$bioreg %in% br,])>0]
 spp.chao <- specpool(mrps.NM[mrps.cneb$bioreg %in% br,])$chao
  se.chao <- specpool(mrps.NM[mrps.cneb$bioreg %in% br,])$chao.se
 spp2.chao <- specpool(mrps.MSM.um[mrps.msm.cneb$bioreg %in% br,])$chao
  se2.chao <- specpool(mrps.MSM.um[mrps.msm.cneb$bioreg %in% br,])$chao.se

  dts.mrps <- rbind(dts.mrps,
                    data.frame(bioreg=br,                    			n.NM=sum(mrps.cneb$bioreg %in% br),
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                               
                    n.LIT=sum(mrps.msm.cneb$bioreg %in% br),
                               spp.LIT=length(unique(spp.LIT)),
                               chao.LIT=spp2.chao,
                               se.LIT=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))
}

spp.NM <- colnames(mrps.NM)
spp.LIT <- colnames(mrps.MSM.um)
spp.chao <- specpool(mrps.NM)$chao
se.chao <- specpool(mrps.NM)$chao.se
spp2.chao <- specpool(mrps.MSM.um)$chao
se2.chao <- specpool(mrps.MSM.um)$chao.se
dts.mrps <- rbind(dts.mrps,
                  data.frame(bioreg="total",
                            	n.NM=nrow(mrps.cneb),
                             spp.NM=length(unique(spp.NM)),
                             chao.NM=spp.chao,
                             se.NM=se.chao,
                             n.LIT=nrow(mrps.MSM.um),
                             spp.LIT=length(unique(spp.LIT)),
                             chao.LIT=spp2.chao,
                             se.LIT=se2.chao,
                             total=length(unique(c(spp.NM,spp.LIT)))))

print(xtable(dts.mrps, caption="Total de especies por bioregión ($S_{obs}$) y estimado de riqueza basados en incidencia de especies ($S_{chao}$) para los muestreos de mariposas Pieridae de NeoMapas, y datos de literatura y colecciones. N: Número de transecciones o celdas de la CNEB. $S_{total}$: total de especies detectadas entre ambas fuentes de datos.",label="TAB:MRPSBRG", align="l|l|ccrl|ccrl|c|"), caption.placement="top", include.rownames=F, include.colnames=F, add.to.row=list(pos=list(-1,5), command=c("\\hline\nBioregión & \\multicolumn{4}{c|}{NeoMapas} & \\multicolumn{4}{c|}{Literatura + Colecciones} & $S_{total}$ \\\\\n & N & $S_{obs}$ & \\multicolumn{2}{c|}{$S_{chao} \\pm SE$} & N & $S_{obs}$ & \\multicolumn{2}{c|}{$S_{chao} \\pm SE$} & \\\\\n ","\\hline\n")), size="small")



###################################################
### chunk number 49: especiesRaras
###################################################
#line 1477 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tt <- colSums(mrps.OTR>0)

t2 <- colSums(mrps.NM>0)
t3 <- colSums(mrps.NM)

df <- merge(data.frame(tt),data.frame(t2),by="row.names",all=T)
df <- merge(df,data.frame(t3),by.x="Row.names",by.y="row.names",all=T)

##layout(matrix(1:2,ncol=1))
##truehist(df[!is.na(df$t2),"tt"],prob=F)
##truehist(df[is.na(df$t2),"tt"],prob=F)
amplia <- paste("\\\\emph{",df[is.na(df$t2) & df$tt>=8,1],"}", sep="", collapse=", ")
reducida <- paste("\\\\textit{",df[is.na(df$t2) & df$tt<8 & df$tt>3,1],"}", sep="", collapse=", ")
restringida <- paste("\\\\textit{",df[is.na(df$t2) & df$tt<4,1],"}", sep="", collapse=", ")

ss <- colSums(mrps.OTR[rownames(mrps.OTR) %in% rownames(mrps.CN),!(colnames(mrps.OTR) %in% colnames(mrps.CN))])



###################################################
### chunk number 50: VariablesDependientesMariposas
###################################################
#line 1504 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
sd.mrps <- apply(log1p(mrps.NM),1,sd)
abnd.mrps <- apply(log1p(mrps.NM),1,mean)
H0.mrps <- diversity(mrps.NM, index = "shannon")
H1.mrps <- renyi(mrps.NM,scales=1)

ACE.mrps <- estimateR(mrps.NM)["S.ACE",]
## realmente, según el código de estimateR es standard deviation
se.mrps <- estimateR(mrps.NM)["se.ACE",]

 print(xtable(data.frame(sfrz=mrps.sfrz, N.jmp=rowSums(mrps.NM), S.obs=rowSums(mrps.NM>0), log.abnd=abnd.mrps,sd.log.abnd=sd.mrps, H0=H0.mrps, ACE=ACE.mrps, se.ACE=se.mrps), caption="Estimados de abundancia, riqueza de especies y diversidad  para los datos de Pieridae de 2006 de NeoMapas. $E$: esfuerzo de muestreo en minutos*colector (m*c), $N_{total}$: número total de ejemplares de Pieridae, $S_{obs}$: número de especies detectadas.", label="TAB:MRPSEST", align="|l|cccrlcrl|", digits=c(0,1,0,0,2,2,2,2,2)), caption.placement="top", include.rownames=T, include.colnames=F, add.to.row=list(pos=list(-1), command=c("\\hline\n & $E$ (m*c) & $N_{total}$ & $S_{obs}$& \\multicolumn{2}{c}{$\\log N (\\mu \\pm \\sigma)$} & $H$ & \\multicolumn{2}{c|}{ACE $\\pm$ s.e.} \\\\")), hline.after=c(-1,length(abnd.mrps)), size="footnotesize")


###################################################
### chunk number 51: DDRGmrps
###################################################
#line 1528 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.mtz <- mrps.NM
d1.mrps <- vegdist(mi.mtz,"chao")
mi.h1 <- as.dendrogram(as.hclust(agnes(d1.mrps)), hang=.2)
mi.mrpp <- meandist(d1.mrps,mrps.cneb$bioreg)
clrs <- brewer.pal(6,"Set1")
names(clrs) <- levels(mrps.cneb$bioreg)
pchs <- c(20,21,22,23,24,25)
names(pchs) <- levels(mrps.cneb$bioreg)

mis.clrs <- mi.col <- clrs[mrps.cneb$bioreg]
##c(0,"thistle","palegoldenrod","peachpuff3","whitesmoke","wheat")
names(mis.clrs) <- rownames(mi.mtz)


mis.pchs <- pchs[mrps.cneb$bioreg]
names(mis.pchs) <- rownames(mi.mtz)

 colLab <- function(n) {
           if(is.leaf(n)) {
             a <- attributes(n)
             attr(n, "nodePar") <- c(a$nodePar,
                                     list(lab.font=1,
                                          col=1,bg=mis.clrs[a$label],
                                          lab.col=1,
                                          pch=mis.pchs[a$label],cex=1.15,
                                          lab.cex=.68))

           }
           n
     }
     h0e <- dendrapply(mi.h1, colLab)

par(mar=c(7,4,1,1))
plot(h0e,ylab="Índice de disimilitud (Chao)")

op <- par(fig=c(0.3,.6,0.55,.85),new=TRUE,xpd=NA,mar=c(1,1,1,1),cex=.8)
cl <- hclust(as.dist(mi.mrpp), method="average")
cl <- as.dendrogram(cl, hang = 0)
w <- diag(mi.mrpp)[labels(cl)]
tr <- unlist(dendrapply(cl, function(n) attr(n, "height")))
root <- attr(cl, "height")
ylim <- range(c(w, tr, root), na.rm = TRUE)
plot(cl, ylim = ylim, leaflab = "none", axes = T, main="Disimilitud media por bioregión")
box(fg="grey66")
for (i in 1:length(w)) segments(i, tr[i], i, w[i])
pos <- ifelse(w < tr, 1, 3)
pos[is.na(pos)] <- 1
w[is.na(w)] <- tr[is.na(w)]
text(1:length(w), w, labels = labels(cl), pos = pos, 
            srt = 0,col=clrs[labels(cl)], font=2)
points(1:length(w), w,bg=clrs[labels(cl)],pch=pchs[labels(cl)])

##legend(4.3,.65,levels(mrps.cneb$bioreg)[-1],pt.bg=brewer.pal(6,"Set1")[-1],pch=21:25,cex=.85,title="Bioregiones")



###################################################
### chunk number 52: tablaAICmrps
###################################################
#line 1593 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## según sugerencia de Kate
	aics.ttl <- data.frame()
	mis.vars <- c("bioreg","PC1","PC2","PC3")
	##mis.vars <- c("bioreg","bosq_avrg","prcp_avrg","tmpm_avrg")
	k <- length(mis.vars)
	tds.vars <- "bioreg+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
	##tds.vars <- "bioreg+PC1+PC2+PC3"
aves.sfrz <- NULL
for (grp in c("mrps")) {
	mi.cneb <- get(paste(grp,"cneb",sep="."))
for (mi.vd in c("abnd","sd","H0","ACE")) {
	
	switch(mi.vd, ACE= {## usamos 1/sigma^2
		mi.wg <- 1/(get(paste("se",grp,sep="."))^2)
		},
		H0={
		mi.wg <- (get(paste(grp,"sfrz",sep=".")))
		},
		sd={
		mi.wg <- (get(paste(grp,"sfrz",sep=".")))
		},
		abnd={
		mi.wg <- (get(paste(grp,"sfrz",sep=".")))
		})
		
	aics <- data.frame()
	mi.lm <- lm(formula(paste(paste(mi.vd, grp, sep="."), "~", 1)),data=mi.cneb, weights=mi.wg)
	aics <- data.frame(grp=grp, vd=mi.vd, formula="~1", R2.adj=summary(mi.lm)$adj.r.squared, logLik=logLik(mi.lm)[[1]], k=mi.lm$rank+1, n=nrow(mi.lm$model), AIC=AIC(mi.lm), deltaAIC=NA, wi=NA, aic.weights=NA)			
	
	for (i in 1:k) {
		kiis <- combinations(k,i)
		for (j in 1:nrow(kiis)) {
			mi.frm <- paste(paste(mi.vd, grp, sep="."), "~", paste(mis.vars[kiis[j,]], collapse="+"))
			mi.lm <- lm(formula(mi.frm),data=mi.cneb, weights=mi.wg)
			aics <- rbind(aics, data.frame(grp=grp,vd=mi.vd, formula=paste(mis.vars[kiis[j,]], collapse="+"), R2.adj=summary(mi.lm)$adj.r.squared, logLik=logLik(mi.lm)[[1]], k=mi.lm$rank+1,n=nrow(mi.lm$model), AIC=AIC(mi.lm), deltaAIC=NA,wi=NA, aic.weights=NA))			
		}
	}
	

	aics$AICc <- aics$AIC + (2*aics$k*(aics$k+1))/(aics$n-aics$k-1)
	
	aics$deltaAIC <- aics$AICc - min(aics$AICc)
	aics$wi<-exp(1)^(-0.5*aics$deltaAIC)
	aics$aic.weights<-aics$wi/sum(aics$wi)

	mi.lm <- lm(formula(paste(paste(mi.vd,grp,sep="."),"~",
			tds.vars)),data=mi.cneb, weights=mi.wg)
	mi.st <- step(mi.lm, trace=FALSE)			
	aics <- rbind(aics, data.frame(grp=grp,vd=mi.vd, formula=paste(colnames(mi.st$model)[2:(ncol(mi.st$model))]
, collapse=" + "), R2.adj=summary(mi.st)$adj.r.squared, logLik=logLik(mi.st)[[1]],k=mi.st$rank+1, n=nrow(mi.st$model),  	AIC=AIC(mi.st),deltaAIC=NA, wi=NA, aic.weights=NA, AICc=NA))

	aics <- rbind(aics, data.frame(grp=grp, vd=mi.vd, formula=tds.vars, R2.adj=summary(mi.lm)$adj.r.squared, logLik=logLik(mi.lm)[[1]],k=mi.lm$rank+1, n=nrow(mi.lm$model), AIC=AIC(mi.lm), deltaAIC=NA, wi=NA,aic.weights=NA, AICc=NA))

	aics.ttl <- rbind(aics.ttl,aics)
}
}

aics.mrps <- aics.ttl[!is.na(aics.ttl$aic.weights), ]

print(xtable(aics.mrps[aics.mrps$aic.weights>0.1, c("vd","formula","k","n","logLik","AICc","aic.weights","R2.adj")], caption="Pesos de Akaike para todos los modelos ajustados a la riqueza, diversidad y abundancia de mariposas Pieridae en los muestreos de 2006.", label="TAB:AICCMRPS", digits=c(0,0,0,0,0,2,2,3,2)), caption.placement="top", size="footnotesize", include.rownames=F)



###################################################
### chunk number 53: ACEPCs
###################################################
#line 1663 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
##plot(env.pc, display = "sites", choices = c(1, 3), scaling = 2, col="grey77",pch=19)
plot(PC3~PC1,CNEB.nm,col="grey77", xlim=c(-2,1))
abline(h=0, v=0, col="aliceblue")
abline(h=0, v=0, lty=3, col="slateblue4")
arrows(x0=0,y0=0,x1=scores(env.pc, display = "species", choices = c(1), scaling = 2),y1=scores(env.pc, display = "species", choices = c(3), scaling = 2), col="maroon", length=.08, angle=25)
abline(h=0,v=0,col="aliceblue")
symbols(mrps.cneb$PC1,mrps.cneb$PC3,circles=ACE.mrps,inches=.1,bg=3,add=T)
text(env.pc, display = "species", choices = c(1, 3), scaling = 2, col=1,cex=.7,font=2)



###################################################
### chunk number 54: 
###################################################
#line 1682 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mrps.cneb <- CNEB.nm[match(rownames(mrps.CN),CNEB.nm$cdg),]
##table(mrps.cneb[,c("bioreg","est.pr")])
##adonis(mi.d1~bioreg+est.pr,scrb.cneb)
tmp <- adonis(d1.mrps~bioreg+(PC1+PC2+PC3),mrps.cneb)
print(tmp)


###################################################
### chunk number 55: ResumenScrb
###################################################
#line 1702 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
trmp.NM$cb <- "sin cebo"
trmp.NM$cb[trmp.NM$cebo %in% c("heces humanas", "bosta de caballo", "heces de perro", "Heces Humanas", "HH", "HV", "Estiercol", "Heces Cerdo", "bosta de vaca",  "heces perro")] <- "heces"

trmp.NM$cb[trmp.NM$cebo %in% c("higado", "pulmon", "Viceras", "insecto", "hongos", "fruta")] <- "otros"


trmp.NM$cb[trmp.NM$cebo %in% c("pescado", "pollo", "rata", "Carne Podrida", "Polo", "Pollo",  "carne podrida", "Pollo Podrido", "carne", "pescado ", "CP", "PP", "Carne de Cerdo", "carne molida",   "Pescado", "POLLO", "PARECE CARNE EN ", "Carne molida")] <- "carne"

##trmp.NM$cebo %in% c("sin cebo", "robado", "??", "", "NA")

tt <- aggregate(scrb.solis[,6:(ncol(scrb.solis)-3)],by=list(scrb.solis$NM,scrb.solis$yr),sum,na.rm=T)
colnames(tt)[1:2] <- c("NM","yr")

tmp1 <- aggregate(trmp.NM[,c("sfrz","njmp")],  by=list(trmp.NM$yr, trmp.NM$NM), sum)
colnames(tmp1) <- c("yr","NM","sfrz","njmp")

tmp2 <- aggregate(trmp.NM[,c("Prd","vst","fch","cb")], by=list(trmp.NM$yr,trmp.NM$NM),luq)
colnames(tmp2)[1:2] <- c("yr","NM")

tmp3 <- merge(tmp1,tmp2)

tmp4 <- data.frame(NM=tt[,1],
			yr=sprintf("20%s",tt[,2]),
			solis.jmp=rowSums(tt[,-(1:2)]),
			solis.spp=rowSums(tt[,-(1:2)]>0))

tmp5 <- aggregate(scr.NM[,c("jmp","spp")],by=list(scr.NM$NM,scr.NM$yr),luq)
colnames(tmp5)[1:2] <- c("NM","yr")

tmp6 <- aggregate(scr.NM[,c("vst")],by=list(scr.NM$NM,scr.NM$yr),luq)
colnames(tmp5)[1:2] <- c("NM","yr")


scrb.rsm <- merge(tmp5,merge(tmp3,tmp4,by=c("NM","yr"),all=T), by=c("NM","yr"), all=T)
##scrb.rsm <- merge(tmp3,tmp4,by=c("NM","yr"),all=T)
scrb.rsm <- scrb.rsm[order(scrb.rsm$yr,scrb.rsm$NM),]


scrb.rsm$jmp[!is.na(scrb.rsm$solis.jmp)] <- scrb.rsm$solis.jmp[!is.na(scrb.rsm$solis.jmp)]

scrb.rsm$spp[!is.na(scrb.rsm$solis.spp)] <- scrb.rsm$solis.spp[!is.na(scrb.rsm$solis.spp)]

scrb.rsm <- scrb.rsm[!is.na(scrb.rsm$sfrz),]
 scrb.rsm$p.id <- scrb.rsm$jmp/scrb.rsm$njmp
scrb.rsm[is.na(scrb.rsm$p.id),"p.id"] <- 0 
scrb.rsm[scrb.rsm$NM %in% c("15","41"),"p.id"] <- 0.01 

scrb.rsm$CNEB <- info.NM$CNEB[match(as.numeric(scrb.rsm$NM),info.NM$NM)]

scrb.rsm$cb <- NA
for (j in 1:nrow(scrb.rsm)) {
	scrb.rsm$cb[j] <- luq(trmp.NM$vst[trmp.NM$NM %in% scrb.rsm$NM[j] & trmp.NM$yr %in% scrb.rsm$yr[j] & trmp.NM$cb=="heces"])*100/luq(trmp.NM$vst[trmp.NM$NM %in% scrb.rsm$NM[j] & trmp.NM$yr %in% scrb.rsm$yr[j]])
}

scr.NM$esp <- paste(scr.NM$genero,scr.NM$especie)
tt <- aggregate(scrb.solis[,6:(ncol(scrb.solis)-3)],by=list(scrb.solis$NM,scrb.solis$yr),sum,na.rm=T)
colnames(tt)[1:2] <- c("NM","yr")
mi.ss <- scr.NM$yr=="2005" & grepl("Scarabaeinae",scr.NM$familia)
t2 <- data.frame(tapply(scr.NM[mi.ss,c("jmp")],list(scr.NM$NM[mi.ss],scr.NM$esp[mi.ss]),luq))
t2[is.na(t2)] <- 0

t2$NM <- rownames(t2)
t2$yr <- "05"

scrb.NM <- merge(tt,t2[c("09","24","26"),],all=T)
scrb.NM[is.na(scrb.NM)] <- 0
scrb.cneb <- CNEB.nm[match(info.NM$CNEB[match(as.numeric(scrb.NM$NM),info.NM$NM)],CNEB.nm$cdg),]

lit.scrb <- rownames(scrb.OTR)
tds.scrb <- info.NM$CNEB[match(as.numeric(unique(scrb.rsm$NM)),info.NM$NM)]
slc.scrb <- info.NM$CNEB[match(as.numeric(unique(scrb.rsm$NM[scrb.rsm$p.id>.1])),info.NM$NM)]

tds.mrps <- info.NM$CNEB[match(as.numeric(unique(tvn.NM$NM[!is.na(tvn.NM$NM)])),info.NM$NM)]
msm.mrps <- rownames(mrps.OTR)


###################################################
### chunk number 56: TablaResumenScrb
###################################################
#line 1779 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tt <- scrb.rsm[, c("NM", "CNEB", "yr", "vst","cb","sfrz","jmp", "spp")]
tt$notas <- ""
tt$notas[scrb.rsm$p.id <.1] <- "**"
tt$notas[scrb.rsm$p.id==0] <- "*"

tab <- xtable(tt,caption="Transecciones incluídas en los muestreos de escarabajos desde 2005 hasta 2010. NM: código de la transección, yr: año de muestreo, vst: número de trampas colocadas, cb: porcentaje de trampas con cebo de heces, sfrz: esfuerzo de muestreo en horas trampa, jmp: número de ejemplares capturados, spp: número de especies identificadas hasta la fecha, notas: * Muestras no identificadas, ** Muestras identificadas parcialmente", label="TAB:SCRB", digits=c(0,0,0,0,0,1,1,0,0,0))

print(tab,include.rownames=F, caption.placement="top", size="footnotesize")



###################################################
### chunk number 57: mapaScrb
###################################################
#line 1795 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
scrb.OTR.um <- scrb.OTR[rownames(scrb.OTR) %in% CNEB.nm$cdg[CNEB.nm$UM==1],]
scrb.OTR.um <- scrb.OTR.um[,colSums(scrb.OTR.um)>0]

## colocamos los márgenes
par(mar=c(1,0,2,0))
##Cuadrícula con los códigos en los bordes
plot(CNEB,border=NA)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)
##title(main="Escarabajos",sub="NeoMapas 2005 a 2010 y Literatura",line=-1)

## Datos de los museos como celdas sombreadas de la cuadrícula
plot(CNEB[CNEB@data$cdg %in% lit.scrb,],add=T,col="grey83",border="grey77")

## División política
plot(vzla,border="maroon",add=T)

## mostramos las celdas muestreadas y la selección de celdas utilizadas
## en este caso
tmp <- table(scrb.rsm$CNEB)

symbols(coordinates(CNEB[CNEB@data$cdg %in% names(tmp),]),circle=rep(1,length(tmp)),inches=.06,add=T)

symbols(coordinates(CNEB[CNEB@data$cdg %in% slc.scrb,]),circle=rep(1,length(unique(slc.scrb))),inches=.06,add=T,fg=1,bg=1)

symbols(coordinates(CNEB[CNEB@data$cdg %in% names(tmp)[tmp>1],]),circle=rep(1,sum(tmp>1)),inches=.06,bg=2,fg=2,add=T)
##invisible(text(coordinates(CNEB[CNEB@data$cdg %in% tds.scrb,]), labels=CNEB@data[CNEB@data$cdg %in% tds.scrb,"cdg"], cex=0.7,col=1))



###################################################
### chunk number 58: 
###################################################
#line 1840 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
fams <- c()
for (i in colnames(scrb.NM)) {
  fams <- c(fams,unique(Escarabajos$trb[match(strsplit(i,"\\.")[[1]][1],Escarabajos$gen)]))
}
fams[grep("Oxysternum",colnames(scrb.NM))] <- "Phanaeini"
fams[grep("Tetramerteia",colnames(scrb.NM))] <- "Phanaeini"
fams[grep("Euristernus",colnames(scrb.NM))] <- "Oniticellini"

mi.tt <- data.frame()
for (i in unique(fams[!is.na(fams)])) {
  mi.tt <- rbind(mi.tt,i=rowSums( rowsum(scrb.NM[,fams %in% i],scrb.cneb$bioreg)))
}
mi.tt <- t(mi.tt)
colnames(mi.tt) <- abbreviate(unique(fams[!is.na(fams)]))
rownames(mi.tt) <- levels(scrb.cneb$bioreg)[-1]
mi.tt <- cbind(mi.tt,total=rowSums(mi.tt))
mi.tt <- rbind(mi.tt,total=colSums(mi.tt))
mi.tt <- mi.tt[c(1,5,2,3,4,6),]
print(xtable(mi.tt,digits=0,caption=paste("Ejemplares capturados por tribu y bioregión.",paste(abbreviate(unique(fams[!is.na(fams)])),": ",unique(fams[!is.na(fams)]),sep="", collapse="; ")), label="TAB:SCRBFAM"),caption.placement="top", size="footnotesize")


###################################################
### chunk number 59: scrbfase1
###################################################
#line 1872 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tmp <- scrb.rsm[scrb.rsm$yr < 2006,]
mtz <- scrb.NM[scrb.NM$yr=="05",-(1:2)]


###################################################
### chunk number 60: scrbfase2
###################################################
#line 1886 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tmp <- scrb.rsm[scrb.rsm$yr == 2006,]
mtz <- scrb.NM[scrb.NM$yr=="06",-(1:2)]


###################################################
### chunk number 61: EjemplaresEsfuerzoScarabaeinae
###################################################
#line 1898 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
cols <- c("2005"=1,"2006"=2,"2008"=1,"2009"=3,"2010"=3)
plot(jmp~sfrz, scrb.rsm,log="y", xlab="Esfuerzo de muestreo", ylab="Ejemplares capturados", col=cols[scrb.rsm$yr],cex=1.9)
text(scrb.rsm$sfrz,scrb.rsm$jmp,scrb.rsm$NM, col=cols[scrb.rsm$yr], cex=.38)
lines(supsmu(scrb.rsm$sfrz,scrb.rsm$jmp,bass=1), lwd=2, col="slateblue4")


###################################################
### chunk number 62: scrbfase3
###################################################
#line 1916 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tmp <- scrb.rsm[scrb.rsm$yr > 2008,]
mtz <- scrb.NM[scrb.NM$yr %in% c("09","10"),-(1:2)]


###################################################
### chunk number 63:  eval=FALSE
###################################################
## #line 1932 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## aggregate(colSums(scrb.OTR),list(Escarabajos[Escarabajos$cdg_especie %in% colnames(scrb.OTR),"tribu"]),sum)
## table(Escarabajos[Escarabajos$cdg_especie %in% unique(scrb.BDV$cdg_taxon[!(scrb.BDV$cdg_taxon %in% scrb.BDV0$cdg_taxon)]),"trb"])


###################################################
### chunk number 64: Selección datos escarabajos
###################################################
#line 1946 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
slc.09 <- scrb.rsm[scrb.rsm$sfrz>3000 & scrb.rsm$p.id>.1 & scrb.rsm$yr %in% c("2009","2010"),"NM"]
##slc.09 <- slc.09[slc.09!="15"]
slc.06 <- scrb.rsm[scrb.rsm$yr %in% "2006" & !(scrb.rsm$NM %in% slc.09),"NM"]

tmp <- rbind(
scrb.NM[scrb.NM$yr %in% c("06") & scrb.NM$NM %in% slc.06,],
scrb.NM[scrb.NM$yr %in% c("09","10") & scrb.NM$NM %in% slc.09,])

rownames(tmp) <- paste(tmp$NM,"::20",tmp$yr,sep="")
tmp <- tmp[,-(1:2)]

scrb.sfrz <- scrb.rsm[match(rownames(tmp), paste(scrb.rsm$NM,scrb.rsm$yr,sep="::")),"sfrz"]
scrb.CN <- scrb.rsm[match(rownames(tmp), paste(scrb.rsm$NM,scrb.rsm$yr,sep="::")),"CNEB"]

scrb.NM <- tmp[,colSums(tmp)>0]
scrb.NM.c <- scrb.NM/(scrb.sfrz/24)

scrb.NM.alpha <- rowSums(scrb.NM>0)
scrb.OTR.alpha <- rowSums(scrb.OTR.um>0)
scrb.cneb <- CNEB.nm[match(scrb.CN,CNEB.nm$cdg),]
scrb.OTR.cneb <- CNEB.nm[match(rownames(scrb.OTR.um),CNEB.nm$cdg),]

##rownames(scrb.NM) <- info.NM$Nombre[match(as.numeric(rownames(scrb.NM)),info.NM$NM)]


###################################################
### chunk number 65: SCRBCRE
###################################################
#line 1976 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
scrb.tau <- specaccum(scrb.NM[,-grep(".sp",colnames(scrb.NM))],method="exact",conditioned=FALSE,gamma="chao")
scrb.tds.tau <- specaccum(scrb.NM,method="exact",conditioned=FALSE,gamma="chao")
scrb.cor.tau <- specaccum(scrb.NM.c,method="exact",conditioned=FALSE,gamma="chao")
##scrb.OTR.tau <- specaccum(scrb.OTR,method="exact",conditioned=FALSE,gamma="chao")
scrb.OTR.tau <- specaccum(scrb.OTR.um,method="exact",conditioned=FALSE,gamma="chao")

max.scrb <- 120
##length(unique(scrb.BDV$especie[scrb.BDV$cdg_ref=="scarabnet"]))

cols <- c("grey42","grey77","pink3")
cols <- c("slateblue4","lightseagreen","pink3")
## realizamos una figura con las curvas para cada grupo
par(mar=c(4,4,2,1))

plot(scrb.OTR.tau, main="Scarabaeinae",
     xlab="", ylab="", 
     ylim=c(0,max.scrb*1.15),
     ci.type="polygon", col = cols[2], ci.col = cols[2],border=cols[2])
plot(scrb.OTR.tau,add=T,ci=0)
##plot(scrb.tds.tau,add=T,ci=0,lty=2)
plot(scrb.tds.tau,add=T,ci.type="polygon", col = 1, ci.col = cols[1], border=cols[1])
plot(scrb.tau,add=T,ci.type="line", col = 1, ci.col = cols[1],ci.lty=2)
abline(h=max.scrb,lty=3,col=2)
text(60,max.scrb+3,paste("ScarabNet:",max.scrb,"spp."),col=2,cex=.8)



###################################################
### chunk number 66: scrbBoxPlot
###################################################
#line 2013 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
cols <- c("grey42","grey77")
cols <- c("slateblue4","lightseagreen","pink3")
## realizamos una figura con las curvas para cada grupo
par(mar=c(4,4,1,1),oma=rep(0,4))

max.y <- 40
boxplot(scrb.NM.alpha~scrb.cneb$bioreg,at=c(1:6)-.25,width=rep(.1,6),xlim=c(1.6,6.47),varwidth=T,boxwex=.3,ylim=c(0,max.y),notch=F,axes=F,col = cols[1])
title(xlab="Bioregión", ylab="Número de especies por celda", line=3)
axis(2)
axis(1,1:6,rep("",6))
axis(1,seq(1,6,by=2),levels(scrb.cneb$bioreg)[seq(1,6,by=2)],line=-.4,lty=0)
axis(1,seq(2,6,by=2),levels(scrb.cneb$bioreg)[seq(2,6,by=2)],line=.4,lty=0)
boxplot(scrb.OTR.alpha~scrb.OTR.cneb$bioreg,at=c(1:6)+.25,add=T,varwidth=T,boxwex=.3,col=cols[2],notch=F,axes=F)
box()

if (any(scrb.OTR.alpha>max.y)) {
	text(as.numeric(scrb.OTR.cneb$bioreg[scrb.OTR.alpha>max.y])+.25,rep(max.y,sum(scrb.OTR.alpha>max.y)),scrb.OTR.alpha[scrb.OTR.alpha>max.y], cex=.7)
	arrows(as.numeric(scrb.OTR.cneb$bioreg[scrb.OTR.alpha>max.y])+.25, rep(max.y-4,sum(scrb.OTR.alpha>max.y)), as.numeric(scrb.OTR.cneb$bioreg[scrb.OTR.alpha>max.y])+.25, rep(max.y-1.2,sum(scrb.OTR.alpha>max.y)),length = 0.03)
}



###################################################
### chunk number 67: TablaSppScrb
###################################################
#line 2044 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
dts.scrb <- data.frame()
for (br in levels(CNEB.nm$bioreg)[-1]) {
 
  spp.NM <- colnames(scrb.NM)[colSums(scrb.NM[scrb.cneb$bioreg %in% br,])>0]
  tmp <- Escarabajos[paste(Escarabajos$genero,Escarabajos$epiteto)  %in% colnames(scrb.OTR.um)[colSums(scrb.OTR.um[scrb.OTR.cneb$bioreg %in% br,])>0],c("genero","epiteto")]
  spp.LIT <- paste(tmp$genero,tmp$epiteto,sep=".")
 spp.chao <- specpool(scrb.NM.c[scrb.cneb$bioreg %in% br,])$chao
  se.chao <- specpool(scrb.NM.c[scrb.cneb$bioreg %in% br,])$chao.se
 spp2.chao <- specpool(scrb.OTR.um[scrb.OTR.cneb$bioreg %in% br,])$chao
  se2.chao <- specpool(scrb.OTR.um[scrb.OTR.cneb$bioreg %in% br,])$chao.se
  dts.scrb <- rbind(dts.scrb,
                    data.frame(bioreg=br,
                    			n.NM=sum(scrb.cneb$bioreg %in% br),
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                               n.LIT=sum(scrb.OTR.cneb$bioreg %in% br),
                    			spp.LIT=length(unique(spp.LIT)),
                               chao.LIT=spp2.chao,
                               se.LIT=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))
}
spp.NM <- colnames(scrb.NM)
tmp <- Escarabajos[paste(Escarabajos$genero,Escarabajos$epiteto) %in% colnames(scrb.OTR.um),]
spp.LIT <- paste(tmp$genero,tmp$epiteto,sep=".")
spp.chao <- specpool(scrb.NM)$chao
se.chao <- specpool(scrb.NM)$chao.se
spp2.chao <- specpool(scrb.OTR.um)$chao
se2.chao <- specpool(scrb.OTR.um)$chao.se

  dts.scrb <- rbind(dts.scrb,
                    data.frame(bioreg="total",
                    			n.NM=nrow(scrb.cneb),
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                               n.LIT=nrow(scrb.OTR.um),
                               spp.LIT=length(unique(spp.LIT)),
                               chao.LIT=spp2.chao,
                               se.LIT=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))


print(xtable(dts.scrb, caption="Total de especies por bioregión ($S_{obs}$) y estimado de riqueza basados en incidencia de especies ($S_{chao}$) para los muestreos de escarabajos coprófagos (Scarabaeinae) de NeoMapas, y datos de literatura. N: Número de transecciones o celdas de la CNEB. $S_{total}$: total de especies detectadas entre ambas fuentes de datos.",label="TAB:SCRBBRG", align="l|l|ccrl|ccrl|c|"), caption.placement="top", include.rownames=F, include.colnames=F, add.to.row=list(pos=list(-1,5), command=c("\\hline\nBioregión & \\multicolumn{4}{c|}{NeoMapas} & \\multicolumn{4}{c|}{Literatura} & $S_{total}$ \\\\\n & N & $S_{obs}$ & \\multicolumn{2}{c|}{$S_{chao} \\pm SE$} & N & $S_{obs}$ & \\multicolumn{2}{c|}{$S_{chao} \\pm SE$} & \\\\\n ","\\hline\n")), size="small")


###################################################
### chunk number 68: VariablesDependientesScarabaeinae
###################################################
#line 2093 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
sd.scrb <- apply(log1p(scrb.NM),1,sd)
abnd.scrb <- apply(log1p(scrb.NM),1,mean)
H0.scrb <- diversity(scrb.NM, index = "shannon")
H1.scrb <- renyi(scrb.NM,scales=1)
ACE.scrb <- estimateR(scrb.NM)["S.ACE",]
## realmente, según el código de estimateR es standard deviation
se.scrb <- estimateR(scrb.NM)["se.ACE",]

scrb.rslt <- xtable(data.frame(sfrz=scrb.sfrz, N.jmp=rowSums(scrb.NM), S.obs=rowSums(scrb.NM>0),log.abnd=abnd.scrb,sd.log.abnd=sd.scrb, H0=H0.scrb, ACE=ACE.scrb, se.ACE=se.scrb), caption="Estimados de abundancia, riqueza de especies y diversidad  para los muestreos de escarabajos coprófago de NeoMapas. $E$: esfuerzo de muestreo en horas*trampa (h*t), $N_{total}$: número total de ejemplares de Scarabaeinae, $S_{obs}$: número de especies detectadas.", label="TAB:SCRBEST", align="|l|cccrlcrl|", digits=c(0,1,0,0,2,2,2,2,2))

 print(scrb.rslt, caption.placement="top", include.rownames=T, include.colnames=F, add.to.row=list(pos=list(-1), command=c("\\hline\n & $E$ (h*t) & $N_{total}$ & $S_{obs}$& \\multicolumn{2}{c}{$\\log N (\\mu \\pm \\sigma)$} & $H$ & \\multicolumn{2}{c|}{ACE $\\pm$ s.e.} \\\\")), hline.after=c(-1,length(abnd.scrb)), size="footnotesize")


###################################################
### chunk number 69: tablaAICscrb
###################################################
#line 2112 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

## según sugerencia de Kate
	aics.ttl <- data.frame()
	mis.vars <- c("bioreg","PC1","PC2","PC3")
	k <- length(mis.vars)
	tds.vars <- "bioreg+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
aves.sfrz <- NULL
for (grp in c("scrb")) {
	mi.cneb <- get(paste(grp,"cneb",sep="."))
for (mi.vd in c("abnd","sd","H0","ACE")) {
	
	switch(mi.vd, ACE= {## usamos 1/sigma^2
		mi.wg <- 1/(get(paste("se",grp,sep="."))^2)
		},
		H0={
		mi.wg <- (get(paste(grp,"sfrz",sep=".")))
		},
		sd={
		mi.wg <- (get(paste(grp,"sfrz",sep=".")))
		},
		abnd={
		mi.wg <- (get(paste(grp,"sfrz",sep=".")))
		})
		
	aics <- data.frame()
	mi.lm <- lm(formula(paste(paste(mi.vd, grp, sep="."), "~", 1)),data=mi.cneb, weights=mi.wg)
	aics <- data.frame(grp=grp, vd=mi.vd, formula="~1", R2.adj=summary(mi.lm)$adj.r.squared, logLik=logLik(mi.lm)[[1]], k=mi.lm$rank+1, n=nrow(mi.lm$model), AIC=AIC(mi.lm), deltaAIC=NA, wi=NA, aic.weights=NA)			
	
	for (i in 1:k) {
		kiis <- combinations(k,i)
		for (j in 1:nrow(kiis)) {
			mi.frm <- paste(paste(mi.vd, grp, sep="."), "~", paste(mis.vars[kiis[j,]], collapse="+"))
			mi.lm <- lm(formula(mi.frm),data=mi.cneb, weights=mi.wg)
			aics <- rbind(aics, data.frame(grp=grp,vd=mi.vd, formula=paste(mis.vars[kiis[j,]], collapse="+"), R2.adj=summary(mi.lm)$adj.r.squared, logLik=logLik(mi.lm)[[1]], k=mi.lm$rank+1,n=nrow(mi.lm$model), AIC=AIC(mi.lm), deltaAIC=NA,wi=NA, aic.weights=NA))			
		}
	}
	

	aics$AICc <- aics$AIC + (2*aics$k*(aics$k+1))/(aics$n-aics$k-1)
	
	aics$deltaAIC <- aics$AICc - min(aics$AICc)
	aics$wi<-exp(1)^(-0.5*aics$deltaAIC)
	aics$aic.weights<-aics$wi/sum(aics$wi)

	aics.ttl <- rbind(aics.ttl,aics)
}
}

aics.scrb <- aics.ttl[ !is.na(aics.ttl$aic.weights), ]

print(xtable(aics.scrb[aics.scrb$aic.weights>0.1, c("vd","formula","k","n","logLik","AICc","aic.weights","R2.adj")], caption="Pesos de Akaike para todos los modelos ajustados", label="TAB:AICCSCRB", digits=c(0,0,0,0,0,2,2,3,2)), caption.placement="top", size="footnotesize", include.rownames=F)



###################################################
### chunk number 70: 
###################################################
#line 2172 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
d1.scrb <- vegdist(scrb.NM,"chao")
##table(scrb.cneb[,c("bioreg","est.pr")])
##adonis(mi.d1~bioreg+est.pr,scrb.cneb)
(tmp <- adonis(d1.scrb~bioreg+(PC1+PC2+PC3),scrb.cneb))


###################################################
### chunk number 71: Fuentes alternativas Aves
###################################################
#line 2236 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
aves.LIT <- t(spp.aves[,grep("NA[1-9]",colnames(spp.aves))])
colnames(aves.LIT) <- tolower(Aves.iocn[match(spp.aves$Latname,Aves.iocn$Latname),"altname"])

aves.LIT.um <- aves.LIT[!(rownames(aves.LIT) %in% "NA5a"),]

rownames(aves.LIT.um) <- trans.info[match(rownames(aves.LIT.um),paste("NA",trans.info$IDTransecta,sep="")),"CNEB"]

aves.GBIF <- CNEB.mtz.Aves[rownames(CNEB.mtz.Aves) %in% CNEB.nm$cdg[CNEB.nm$VEN==1],]
aves.GBIF <- aves.GBIF[,colSums(aves.GBIF)>0]
aves.GBIF.completo <- aves.GBIF
aves.GBIF <- aves.GBIF[,tolower(colnames(aves.GBIF)) %in% tolower(Aves.iocn$altname)]

aves.GBIF.um <- CNEB.mtz.Aves[rownames(CNEB.mtz.Aves) %in% CNEB.nm$cdg[CNEB.nm$UM==1],]
aves.GBIF.um <- aves.GBIF.um[,colSums(aves.GBIF.um)>0]

aves.GBIF.um <- aves.GBIF.um[,tolower(colnames(aves.GBIF.um)) %in% tolower(Aves.iocn$altname)]


colnames(aves.GBIF.um) <- tolower(Aves.txn[match(colnames(aves.GBIF.um),Aves.txn$spp),"altname"])

spps <- unique(c(colnames(aves.LIT.um),colnames(aves.NM), colnames(aves.GBIF.um)))

table(spps %in% colnames(aves.NM),
spps %in% colnames(aves.GBIF.um))



###################################################
### chunk number 72: mapas
###################################################
#line 2269 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## colocamos los márgenes
par(mar=c(1,0,2,0))

##Aves
##Cuadrícula con los códigos en los bordes
plot(CNEB,border=NA)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)

## Datos de GBIF como celdas sombreadas de la cuadrícula
plot(CNEB[CNEB@data$cdg %in% rownames(aves.GBIF),],add=T,col="grey83",border="grey77")

## División política
plot(vzla,border="maroon",add=T)

## mostramos las celdas muestreadas 
symbols(coordinates(CNEB[CNEB@data$cdg %in% slc.avs,]),circle=rep(1,length(slc.avs)),inches=.06,add=T,fg=1,bg=1)
##invisible(text(coordinates(CNEB[CNEB@data$cdg %in% slc.avs,]), labels=CNEB@data[CNEB@data$cdg %in% slc.avs,"cdg"], cex=0.7,col="white"))


###################################################
### chunk number 73: MTRZaves
###################################################
#line 2300 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

d1.aves <- vegdist(aves.NM,"chao")
d2.aves <- dist.binary(aves.LIT.um,method=5)## Sørensen

mx.aves <- as.matrix(d1.aves)
tmp <- as.matrix(d2.aves)
mx.aves[upper.tri(mx.aves)] <- tmp[upper.tri(tmp)]

mx.aves <- mx.aves[order(trans.info[match(rownames(aves.NM),trans.info$Nombretransecta),"Bioregión"]),order(trans.info[match(rownames(aves.NM),trans.info$Nombretransecta),"Bioregión"])]

par(mar=c(7,3,1,7),oma=c(0,0,0,0),xpd=F)
image(1:nrow(mx.aves),1:ncol(mx.aves),mx.aves, axes=F, col=c(NA,rev(brewer.pal(11,"BrBG"))), xlab="", ylab="")

abline(v=cumsum(table(trans.info[match(rownames(mx.aves),trans.info$Nombretransecta),"Bioregión"]))+.5, lty=3)
abline(h=cumsum(table(trans.info[match(rownames(mx.aves),trans.info$Nombretransecta),"Bioregión"]))+.5, lty=3)

polygon(cumsum(table(trans.info[match(rownames(mx.aves),trans.info$Nombretransecta),"Bioregión"]))[c(1,2,2,3,3,4,4,5,5,6,6,5,5,4,4,3,3,2,2,1,1)]+.5,cumsum(table(trans.info[match(rownames(mx.aves),trans.info$Nombretransecta),"Bioregión"]))[c(1,1,2,2,3,3,4,4,5,5,6,6,5,5,4,4,3,3,2,2,1)]+.5,lwd=2)

atta <- aggregate(1:ncol(mx.aves),by=list(trans.info[match(rownames(mx.aves),trans.info$Nombretransecta),"Bioregión"]),mean)$x

axis(2,atta[c(1,3,5)],levels(trans.info$Bioregión)[c(2,4,6)],line=.01,lty=0)
axis(2,atta[c(2,4)],levels(trans.info$Bioregión)[c(3,5)],line=.7,lty=0)
axis(2,atta,rep("",5))
axis(1,las=2,1:ncol(mx.aves),rownames(mx.aves),cex.axis=.5)
axis(4,las=2,1:ncol(mx.aves),rownames(mx.aves),cex.axis=.5)
text(rep(1:nrow(mx.aves),ncol(mx.aves)), rep(1:nrow(mx.aves),rep(ncol(mx.aves),ncol(mx.aves))) ,round(as.vector(mx.aves),2),cex=.5)



###################################################
### chunk number 74: Tabla ordenes aves
###################################################
#line 2339 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
  
aves.cneb <- CNEB.nm[match(rownames(aves.CN),CNEB.nm$cdg),]

##table(Aves.iocn$Order[match(colnames(aves.NM),Aves.iocn$Latname)])

ordo <- Aves.iocn$Order[match(colnames(aves.NM),tolower(Aves.iocn$altname))]
ordo <- paste(toupper(substring(ordo,1,1)),tolower(substring(ordo,2)),sep="")
ordo[ordo %in% 
(otros <- c('Anseriformes','Ciconiiformes','Coraciiformes', 'Galliformes',         'Gruiformes','Opisthocomiformes', 'Cuculiformes', 'Falconiformes',
 'Podicipediformes', 'Strigiformes', 'Suliformes','Tinamiformes',     'Trogoniformes'))] <- "Otros"

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
xtab <- xtable(mi.tt,digits=0,caption=paste("Número de ejemplares observados/escuchados por orden y bioregión. ",paste(abbreviate(unique(ordo)),": ",unique(ordo), " (",table(ordo)[unique(ordo)], " spp.)",
sep="", collapse="; "),". Otros=",paste(otros, collapse=", "),sep=""))

print(xtab, caption.placement="top")


###################################################
### chunk number 75: EspecieEjemplares
###################################################
#line 2374 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.ss <- NM.m1$Especieid != 1379
njmp <- aggregate(NM.m1[mi.ss,"Ntotal"],by=list(NM.m1$IDTransecta[mi.ss]),sum)
nspp <- aggregate(NM.m1[mi.ss,"Especieid"],by=list(NM.m1$IDTransecta[mi.ss]), luq)

mi.dt <- data.frame(NM=trans.info[match(njmp$Group.1,trans.info$IDTransecta),"NM"],
njmp=njmp$x,nspp=nspp$x)

plot(nspp~njmp, mi.dt, cex=1.9,log="y", xlab="Individuos registrados", ylab="Especies detectadas")
text(mi.dt$njmp,mi.dt$nspp,mi.dt$NM,cex=.39)
lines(supsmu(mi.dt$njmp,mi.dt$nspp,bass=1), lwd=2, col="slateblue4")


###################################################
### chunk number 76: 
###################################################
#line 2393 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.ss <- NM.m1$Especieid != 1379
rsm <- aggregate(NM.m1[mi.ss, c("Nvistos.100", "Nvistos.100.1", "Noidos", "Ntotal")], by=list(NM.m1$IDTransecta[mi.ss]),sum)

rsm <- merge(rsm, nspp)
rsm$NM <- trans.info[match(njmp$Group.1,trans.info$IDTransecta),"NM"]
rsm$CNEB <- trans.info[match(njmp$Group.1,trans.info$IDTransecta),"CNEB"]

colnames(rsm) <- c("Id", "cerca", "lejos", "oidos", "total", "especies", "NM", "CNEB")
print(xtable(rsm[,c("NM","CNEB","cerca","lejos","oidos","total","especies")], caption="Número de individuos observados (cerca: $< 100m$, lejos: $> 100m$) y escuchados, y número de especies detectadas en las transecciones de NeoMapas.", label="TAB:OBSAVS"),caption.placement="top",include.rownames=F)


###################################################
### chunk number 77:  eval=FALSE
###################################################
## #line 2405 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## mi.ss <- NM.m2$Especieid != 1379
## rsm <- aggregate(NM.m2[mi.ss, c("Especieid")], by=list(NM.m2$IDTransecta[mi.ss],NM.m2$Lapso[mi.ss]),luq)
## 


###################################################
### chunk number 78: RqzAlfa
###################################################
#line 2411 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

aves.NM.alpha <- rowSums(aves.NM>0)
aves.LIT.alpha <- rowSums(aves.LIT>0)[!(rownames(aves.LIT) %in% "NA5a")]
aves.GBIF.alpha <- rowSums(aves.GBIF.um>0)

aves.cneb <- aves.lit.cneb <- CNEB.nm[match(rownames(aves.CN),CNEB.nm$cdg),]
aves.gbif.cneb <- CNEB.nm[match(rownames(aves.GBIF.um),CNEB.nm$cdg),]


###################################################
### chunk number 79: TablaSppAves
###################################################
#line 2424 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
dts.aves <- data.frame()
for (br in levels(CNEB.nm$bioreg)[-1]) {

  spp.NM <- colnames(aves.NM)[colSums(aves.NM[aves.cneb$bioreg %in% br,])>0]
  spp.LIT <- colnames(aves.LIT.um)[colSums(aves.LIT.um[aves.lit.cneb$bioreg %in% br,])>0]
  spp.GBIF <- colnames(aves.GBIF.um)[colSums(aves.GBIF.um[aves.gbif.cneb$bioreg %in% br,])>0]
  spp.chao <- specpool(aves.NM[aves.cneb$bioreg %in% br,])$chao
  se.chao <- specpool(aves.NM[aves.cneb$bioreg %in% br,])$chao.se
  ##spp2.chao <- specpool(aves.LIT.um[aves.cneb$bioreg %in% br,])$chao
  ##se2.chao <- specpool(aves.LIT.um[aves.cneb$bioreg %in% br,])$chao.se
  spp2.chao <- specpool(aves.GBIF.um[aves.gbif.cneb$bioreg %in% br,])$chao
  se2.chao <- specpool(aves.GBIF.um[aves.gbif.cneb$bioreg %in% br,])$chao.se

  dts.aves <- rbind(dts.aves,
                    data.frame(bioreg=br,                    			n.NM=sum(aves.cneb$bioreg %in% br),
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                    n.GBIF=sum(aves.gbif.cneb$bioreg %in% br),
                               spp.GBIF=length(unique(spp.GBIF)),
                               chao.GBIF=spp2.chao,
                               se.GBIF=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))

}

  spp.NM <- colnames(aves.NM)[colSums(aves.NM)>0]
  spp.LIT <- colnames(aves.LIT.um)[colSums(aves.LIT.um)>0]
  spp.GBIF <- colnames(aves.GBIF.um)[colSums(aves.GBIF.um)>0]
   spp.chao <- specpool(aves.NM)$chao
  se.chao <- specpool(aves.NM)$chao.se
##   spp2.chao <- specpool(aves.LIT.um)$chao
##  se2.chao <- specpool(aves.LIT.um)$chao.se
   spp2.chao <- specpool(aves.GBIF.um)$chao
  se2.chao <- specpool(aves.GBIF.um)$chao.se
 dts.aves <- rbind(dts.aves,
                    data.frame(bioreg="total",
                    	n.NM=nrow(aves.cneb),
                               spp.NM=length(unique(spp.NM)),
                               chao.NM=spp.chao,
                               se.NM=se.chao,
                        n.GBIF=nrow(aves.GBIF.um),
                               spp.GBIF=length(unique(spp.GBIF)),
                               chao.GBIF=spp2.chao,
                               se.GBIF=se2.chao,
                               total=length(unique(c(spp.NM,spp.LIT)))))


print(xtable(dts.aves, caption="Total de especies por bioregión ($S_{obs}$) y estimado de riqueza basados en incidencia de especies ($S_{chao}$) para los muestreos de aves de NeoMapas, y datos de GBIF. N: Número de transecciones o celdas de la CNEB. $S_{esp}$: total de especies esperadas según la revisión de literatura previa a los muestreos.",label="TAB:AVSBRG", align="l|l|ccrl|ccrl|c|"), caption.placement="top", include.rownames=F, include.colnames=F, add.to.row=list(pos=list(-1,5), command=c("\\hline\nBioregión & \\multicolumn{4}{c|}{NeoMapas} & \\multicolumn{4}{c|}{GBIF} & $S_{esp}$ \\\\\n & N & $S_{obs}$ & \\multicolumn{2}{c|}{$S_{chao} \\pm SE$} & N & $S_{obs}$ & \\multicolumn{2}{c|}{$S_{chao} \\pm SE$} & \\\\\n ","\\hline\n")), size="small")


###################################################
### chunk number 80: 
###################################################
#line 2479 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
tmp <- rowSums(aves.GBIF[rownames(aves.GBIF) %in% rownames(aves.LIT.um),]>0)
tmp2 <- rowSums(aves.LIT.um[match(names(tmp),rownames(aves.LIT.um)),]>0)


###################################################
### chunk number 81: VariablesDependientes
###################################################
#line 2487 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
sd.aves <- apply(log1p(aves.NM),1,sd)
abnd.aves <- apply(log1p(aves.NM),1,mean)
H0.aves <- diversity(aves.NM, index = "shannon")
H1.aves <- renyi(aves.NM,scales=1)
ACE.aves <- estimateR(aves.NM)["S.ACE",]
## realmente, según el código de estimateR es standard deviation
se.aves <- estimateR(aves.NM)["se.ACE",]

print(xtable(data.frame(N.jmp=rowSums(aves.NM), S.obs=rowSums(aves.NM>0),log.abnd=abnd.aves,sd.log.abnd=sd.aves, H0=H0.aves, ACE=ACE.aves, se.ACE=se.aves), caption="Estimados de abundancia, riqueza de especies y diversidad  para los muestreos de aves de NeoMapas. $N_{total}$: número total de individuos de aves detectados durante el muestreo 1, $S_{obs}$: número de especies detectadas. El esfuerzo de muestreo fue de 150 min. de observación para todas las transecciones.", label="TAB:AVSEST", align="|l|ccrlcrl|", digits=c(0,0,0,2,2,2,2,2)), caption.placement="top", include.rownames=T, include.colnames=F, add.to.row=list(pos=list(-1), command=c("\\hline\n & $N_{total}$ & $S_{obs}$& \\multicolumn{2}{c}{$\\log N (\\mu \\pm \\sigma)$} & $H$ & \\multicolumn{2}{c|}{ACE $\\pm$ s.e.} \\\\")), hline.after=c(-1,length(abnd.aves)), size="footnotesize")


###################################################
### chunk number 82: tablaAICaves
###################################################
#line 2510 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

## según sugerencia de Kate
	aics.ttl <- data.frame()
	mis.vars <- c("bioreg","PC1","PC2","PC3")
	k <- length(mis.vars)
	tds.vars <- "bioreg+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
aves.sfrz <- NULL
for (grp in c("aves")) {
	mi.cneb <- get(paste(grp,"cneb",sep="."))
for (mi.vd in c("sd","H0","ACE")) {
	
	switch(mi.vd, ACE= {## usamos 1/sigma^2
		mi.wg <- 1/(get(paste("se",grp,sep="."))^2)
		},
		H0={
		mi.wg <- NULL
		},
		sd={
		mi.wg <- NULL
		})
		
	aics <- data.frame()
	mi.lm <- lm(formula(paste(paste(mi.vd, grp, sep="."), "~", 1)),data=mi.cneb, weights=mi.wg)
	aics <- data.frame(grp=grp, vd=mi.vd, formula="~1", R2.adj=summary(mi.lm)$adj.r.squared, logLik=logLik(mi.lm)[[1]], k=mi.lm$rank+1, n=nrow(mi.lm$model), AIC=AIC(mi.lm), deltaAIC=NA, wi=NA, aic.weights=NA)			
	
	for (i in 1:k) {
		kiis <- combinations(k,i)
		for (j in 1:nrow(kiis)) {
			mi.frm <- paste(paste(mi.vd, grp, sep="."), "~", paste(mis.vars[kiis[j,]], collapse="+"))
			mi.lm <- lm(formula(mi.frm),data=mi.cneb, weights=mi.wg)
			aics <- rbind(aics, data.frame(grp=grp,vd=mi.vd, formula=paste(mis.vars[kiis[j,]], collapse="+"), R2.adj=summary(mi.lm)$adj.r.squared, logLik=logLik(mi.lm)[[1]], k=mi.lm$rank+1,n=nrow(mi.lm$model), AIC=AIC(mi.lm), deltaAIC=NA,wi=NA, aic.weights=NA))			
		}
	}
	

	aics$AICc <- aics$AIC + (2*aics$k*(aics$k+1))/(aics$n-aics$k-1)
	
	aics$deltaAIC <- aics$AICc - min(aics$AICc)
	aics$wi<-exp(1)^(-0.5*aics$deltaAIC)
	aics$aic.weights<-aics$wi/sum(aics$wi)

	aics.ttl <- rbind(aics.ttl,aics)
}
}

aics.aves <- aics.ttl[ !is.na(aics.ttl$aic.weights), ]
print(xtable(aics.aves[aics.aves$aic.weights>0.1 , c("vd","formula","k","n","logLik","AICc","aic.weights","R2.adj")], caption="Pesos de Akaike para todos los modelos ajustados", label="TAB:AICCAVS", digits=c(0,0,0,0,0,2,2,3,2)), caption.placement="top", size="footnotesize", include.rownames=F)



###################################################
### chunk number 83: 
###################################################
#line 2566 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
d1.aves <- vegdist(aves.NM,"chao")
(tmp <- adonis(d1.aves~bioreg+(PC1+PC2+PC3),aves.cneb))


###################################################
### chunk number 84: 
###################################################
#line 2592 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## solo las especies con nombre

aves.tau <- specaccum(aves.NM,method="exact",conditioned=FALSE,gamma="chao")
aves.lit.tau <- specaccum(aves.LIT.um,method="exact",conditioned=FALSE,gamma="chao")
aves.gbif <- specaccum(aves.GBIF.um,method="exact",conditioned=FALSE,gamma="chao")
  
max.aves <- 1383



###################################################
### chunk number 85: 
###################################################
#line 2640 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
bib1 <- mrps.BDV[substr(mrps.BDV$cdg,0,3) %in% "097",]
bib0 <- mrps.BDV0[substr(mrps.BDV0$cdg,0,3) %in% "097",]


###################################################
### chunk number 86: Bibliografía Mariposas
###################################################
#line 2662 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
kk <-1
tmp <- Mariposas[Mariposas$scdg %in% "a" & Mariposas$ecdg %in% unique(bib1$cdg),]
for (j in 1:nrow(tmp)){
	nns <- paste("\\textit{",tmp$gen[j]," ",tmp$esp[j],"}",sep="")	
	ccs <- unique(bib0$cneb[bib0$cdg %in% tmp$ecdg[j]])
	ccs <- ccs[!is.na(ccs)]
	ccs <- paste(ccs[order(ccs)],collapse="; ")
	tts <- paste("\\citet{",paste(unique(bib1$cdg_ref[bib1$cdg %in% tmp$ecdg[j]]), collapse=","),"}",sep="")
	cat(paste("{\\footnotesize ", kk,"} & {\\small ",nns," & {\\footnotesize ",tts,"} & {\\small ",ccs,"} \\\\" ,sep=""))
	kk <- kk+1
}



###################################################
### chunk number 87: Bibliografía Escarabajos
###################################################
#line 2702 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
kk <-1
tmp <- Escarabajos[Escarabajos$cdg_especie %in% unique(scrb.BDV$ecdg),]
for (j in 1:nrow(tmp)){
	nns <- paste("\\textit{",tmp$gen[j]," ",tmp$epiteto[j],"}",sep="")	
	ccs <- unique(scrb.BDV0$cneb[scrb.BDV0$ecdg %in% tmp$cdg_especie[j]])
	ccs <- ccs[!is.na(ccs)]
	ccs <- paste(ccs[order(ccs)],collapse="; ")
	tts <- paste("\\citet{",paste(unique(scrb.BDV$cdg_ref[scrb.BDV$ecdg %in% tmp$cdg_especie[j]]), collapse=","),"}",sep="")
	cat(paste("{\\footnotesize ", kk,"} & {\\small ",nns," & {\\footnotesize ",tts,"} & {\\small ",ccs,"} \\\\" ,sep=""))
	kk <- kk+1
}


###################################################
### chunk number 88: Composicion
###################################################
#line 2743 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
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
### chunk number 89: ManuscritoMapa1
###################################################
#line 2783 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
##plot(vegnsn)
##polygon(c(750000, 950000,950000,750000),c(390000,390000,400000,400000))
##polygon(c(750000, 850000,850000,750000),c(390000,390000,400000,400000),col=1)
##text(c(750000,850000,950000),c(420000,420000,420000),c("0","100","200"),cex=.5)
##text(c(850000),c(370000),c("km"),cex=.5)

xy1 <- data.frame(x=c(750000, 950000,950000,750000),y=c(290000,290000,300000,300000))
coordinates(xy1) <- c("x","y")
proj4string(xy1)<- vegnsn@proj4string

xy2 <- data.frame(x=c(750000, 850000,850000,750000),y=c(290000,290000,300000,300000))
coordinates(xy2) <- c("x","y")
proj4string(xy2)<- vegnsn@proj4string

xy3 <- data.frame(x=c(750000,850000,950000,850000),y=c(320000,320000,320000,270000))
coordinates(xy3) <- c("x","y")
proj4string(xy3)<- vegnsn@proj4string

ll1 <- rgdal::spTransform(xy1,vzla@proj4string)
ll2 <- rgdal::spTransform(xy2,vzla@proj4string)
ll3 <- rgdal::spTransform(xy3,vzla@proj4string)

op <- par(fig=c(0.20,1,0,1),xpd=F,mar=c(0.1,0.1,0.1,0.1))

mi.ds <- c(0,30,70,20,35,40)
mi.ag <- c(0,45,135,135,45,135)
mi.cl <- c(0,"grey36","grey56","grey46","grey56","grey46")
plot(vzla,col=NA,border=NA,ylim=c(0,13))
plot(CNEB[CNEB@data$cdg %in% VBG$cdg[VBG$vzla>0],],border="grey77",add=T)
text(rep(-73.5,23),unique(CNEB@data$y),unique(CNEB@data$fila),cex=.5)
text(unique(CNEB@data$x),rep(12.4,27),unique(CNEB@data$col),cex=.5)
br <- CNEB.nm$bioreg[order(CNEB.nm$row,CNEB.nm$col)]
br[!(CNEB.nm$cdg[order(CNEB.nm$row,CNEB.nm$col)] %in% um)] <- "Otro"
plot(CNEB,density=mi.ds[br],
     border=mi.cl[br],
     angle=mi.ag[br],col=mi.cl[br],add=T)
## División política
plot(vzla,border=1,add=T)
plot(vzla[vzla@data$ID %in% 36,],add=T,border=1,col=NA)

polygon(coordinates(ll1))
polygon(coordinates(ll2),col=1)
text(coordinates(ll3),c("0","100","200","km"),cex=.9)
text(coordinates(ll3)[2,],"km",cex=.8)

xs <- seq(-72,-60,by=2)
ys <- rep(0.5,length(xs))
points(xs,ys,pch=3,cex=.7)
text(xs,ys-.3,paste(abs(xs),"º00' W",sep=""),cex=.7)

ys <- seq(1,11,by=2)
xs <- rep(-58,length(ys))
points(xs,ys,pch=3,cex=.7)
text(xs+.3,ys,paste(abs(ys),"º00' N",sep=""),cex=.7,srt=90)


op <- par(fig=c(0,.19,0,.5),new=TRUE,xpd=T,mar=c(0,0,0,0))
plot(c(0,1),c(0,1),axes=F,xlab=NA,ylab=NA,pch=NA)
legend(0,1,c("Occident","Andean mountains","Coastal mountains","Orinoco floodplain","Guayana shield"),density=mi.ds[-1], angle=mi.ag[-1], bty="n", cex=.9, col=mi.cl[-1],border=mi.cl[-1])



## agregar mapa de ubicación continental [bajar resolución...]
op <- par(fig=c(0,.19,.5,1),new=TRUE,xpd=F,mar=c(0,0,0,0))
plot(amgnsn,border="grey73", xlim=c(-4500004,4000004),ylim=c(-6568465,8264709))
abline(h=0,lty=3,col=1)
plot(amgnsn[amgnsn@data$COUNTRY %in% "VENEZUELA",],col=1,add=T)
box()
par(op)



###################################################
### chunk number 90: ManuscritoMapas
###################################################
#line 2859 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
cols <- c("grey42","grey77","pink3")
## realizamos una figura con los mapas de cada grupo
layout(matrix(1:4,ncol=2,byrow=T))
par(mar=c(1,0,2,0))
    


##Pieridae
plot(CNEB,border=NA)
plot(CNEB[CNEB@data$cdg %in% unique(mrps.otrs$cneb) & CNEB@data$vzla>0,], add=T, col="grey83",border="grey77")
plot(vzla,border="grey33",add=T)
tmp <- info.NM$CNEB[match(as.numeric(unique(tvn.NM$NM)),info.NM$NM)]
mrps.nms <- unique(tmp[!is.na(tmp)])
symbols(coordinates(CNEB[CNEB@data$cdg %in% mrps.nms,]),circle=rep(1,length(mrps.nms)),inches=.06,add=T,fg=1)
symbols(coordinates(CNEB[CNEB@data$cdg %in% rownames(mrps.CN),]),circle=rep(1,nrow(mrps.CN)),inches=.06,add=T,fg=1,bg=1)
title(main="(a)",adj=0)

polygon(coordinates(ll1),lwd=5)
text(coordinates(ll3)[c(1,3),1],coordinates(ll3)[c(1,3),2]+.3,c("0","200km"),cex=1.1)

##Escarabajos
plot(CNEB,border=NA)
plot(CNEB[CNEB@data$cdg %in% lit.scrb & CNEB@data$vzla>0,], add=T, col="grey83", border="grey77")
plot(vzla,border="grey33",add=T)
scrb.nms <- unique(scrb.rsm$CNEB[scrb.rsm$NM %in% c(slc.06,slc.09)])
symbols(coordinates(CNEB[CNEB@data$cdg %in% unique(scrb.rsm$CNEB),]),circle=rep(1,luq(scrb.rsm$CNEB)),inches=.06,add=T,fg=1)
symbols(coordinates(CNEB[CNEB@data$cdg %in% scrb.nms,]),circle=rep(1,length(scrb.nms)),inches=.06,add=T,fg=1,bg=1)
title(main="(b)",adj=0)


##Aves
plot(CNEB,border=NA)
plot(CNEB[CNEB@data$cdg %in% rownames(aves.GBIF),],add=T,col="grey83",border="grey77")
plot(vzla,border="grey33",add=T)
symbols(coordinates(CNEB[CNEB@data$cdg %in% slc.avs,]),circle=rep(1,length(slc.avs)),inches=.06,add=T,fg=1,bg=1)
title(main="(c)",adj=0)


###################################################
### chunk number 91: ManuscritoCRE
###################################################
#line 2905 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
cols <- c("grey42","grey77","pink3")
## realizamos una figura con las curvas para cada grupo
layout(matrix(1:4,ncol=2,byrow=T))
par(mar=c(3,3,2,1), oma=c(3,3,0,0))

plot(mrps.msm.tau,
     xlab="", ylab="",  
     ylim=c(0,max.mrps*1.10),
     ci.type="polygon", col = cols[2], ci.col = cols[2],border=cols[2])
     title(main="(a)",adj=0)

plot(mrps.msm.tau,add=T,ci=0)
plot(mrps.tau,add=T,ci.type="polygon", col = 1, ci.col = cols[1],border=cols[1])
abline(h=max.mrps,lty=3,col=1)
text(25,max.mrps+3,paste("Viloria (1990):",max.mrps,"spp."),col=1,cex=.8)


plot(scrb.OTR.tau, 
     xlab="", ylab="", 
     ylim=c(0,max.scrb*1.10),
     ci.type="polygon", col = cols[2], ci.col = cols[2],border=cols[2])
     title(main="(b)",adj=0)

plot(scrb.OTR.tau,add=T,ci=0)
##plot(scrb.tds.tau,add=T,ci=0,lty=2)
plot(scrb.tds.tau,add=T,ci.type="polygon", col = 1, ci.col = cols[1],border=cols[1])
plot(scrb.tau,add=T,ci.type="line", col = 1, ci.col = cols[1],ci.lty=2)
abline(h=max.scrb,lty=3,col=1)
text(60,max.scrb+3,paste("ScarabNet:",max.scrb,"spp."),col=1,cex=.8)

plot(aves.gbif,
     xlab="", ylab="", 
     ylim=c(0,max.aves*1.10),
     ci.type="polygon", col = cols[2], ci.col = cols[2],border=cols[2])
     title(main="(c)",adj=0)

plot(aves.gbif,add=T,ci=0)
plot(aves.tau,add=T,ci.type="polygon", col = 1, ci.col = cols[1],border=cols[1])
##plot(aves.lit.tau,add=T,ci.type="polygon", col = 2, ci.col = cols[3],border=cols[3])

abline(h=max.aves,lty=3,col=1)
text(40,max.aves+30,paste("Hilty (2003):",max.aves,"spp."),col=1,cex=.8)
title(xlab="VBG cells", ylab="Accumulated species richness", outer=T,line=1)

symbols(c(1,2),c(1,1),boxplots=matrix(c(.5,2,0,0,.5,.5,2,0,0,.5),ncol=5,byrow=T),bg=cols,fg=cols,lwd=2,xlim=c(0.7,3),axes=F)
text(c(1.5,2.5),c(1,1),c("NeoMapas","Other\nsources"))



###################################################
### chunk number 92: ManuscritoBoxplot
###################################################
#line 2963 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
cols <- c("grey42","grey77")
labs <- c("otro","Occident","Andean mountains","Coastal mountains","Llanos","Guayana")

## realizamos una figura con las curvas para cada grupo
layout(matrix(1:4,ncol=2,byrow=T))
par(mar=c(3,3,2,1), oma=c(3,3,0,0))


max.y <- 60
boxplot(mrps.NM.alpha~mrps.cneb$bioreg,at=c(1:6)-.25,width=rep(.1,6),xlim=c(1.6,6.47),varwidth=T,boxwex=.3,ylim=c(0,max.y),notch=F,axes=F,col = cols[1])
title(main="(a)",adj=0)

axis(2)
axis(1,1:6,rep("",6))
##mrps.cneb$bioreg[mrps.NM.alpha>40]
##boxplot(rowSums(mrps.OTR>0)~m2.cneb$bioreg)
boxplot(mrps.MSM.alpha~mrps.msm.cneb$bioreg,at=c(1:6)+.25,add=T,varwidth=T,boxwex=.3,col="grey74",notch=F,axes=F)
box()

text(as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>max.y])+.25,rep(max.y,sum(mrps.MSM.alpha>max.y)),mrps.MSM.alpha[mrps.MSM.alpha>max.y], cex=.7)
arrows(as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>max.y])+.25, rep(max.y-4,sum(mrps.MSM.alpha>max.y)), as.numeric(mrps.msm.cneb$bioreg[mrps.MSM.alpha>max.y])+.25, rep(max.y-1.2,sum(mrps.MSM.alpha>max.y)),length = 0.03)

max.y <- 40
boxplot(scrb.NM.alpha~scrb.cneb$bioreg,at=c(1:6)-.25,width=rep(.1,6),xlim=c(1.6,6.47),varwidth=T,boxwex=.3,ylim=c(0,max.y),notch=F,axes=F,col = cols[1])
title(main="(b)",adj=0)
axis(2)
axis(1,1:6,rep("",6))
##scrb.cneb$bioreg[scrb.NM.alpha>40]
##boxplot(rowSums(mrps.OTR>0)~m2.cneb$bioreg)
boxplot(scrb.OTR.alpha~scrb.OTR.cneb$bioreg,at=c(1:6)+.25,add=T,varwidth=T,boxwex=.3,col="grey74",notch=F,axes=F)
box()

text(as.numeric(scrb.OTR.cneb$bioreg[scrb.OTR.alpha>max.y])+.25,rep(max.y,sum(scrb.OTR.alpha>max.y)),scrb.OTR.alpha[scrb.OTR.alpha>max.y], cex=.7)
arrows(as.numeric(scrb.OTR.cneb$bioreg[scrb.OTR.alpha>max.y])+.25, rep(max.y-4,sum(scrb.OTR.alpha>max.y)), as.numeric(scrb.OTR.cneb$bioreg[scrb.OTR.alpha>max.y])+.25, rep(max.y-1.2,sum(scrb.OTR.alpha>max.y)),length = 0.03)

max.y <- 350
boxplot(aves.NM.alpha~aves.cneb$bioreg,at=c(1:6)-.25,width=rep(.1,6),xlim=c(1.6,6.47),varwidth=T,boxwex=.3,ylim=c(0,max.y),notch=F,axes=F,col = cols[1])
title(main="(c)",adj=0)
axis(2)
axis(1,seq(1,6,by=2),labs[seq(1,6,by=2)])
axis(1,seq(2,6,by=2),labs[seq(2,6,by=2)],line=1.1,lty=0)
boxplot(aves.GBIF.alpha~aves.gbif.cneb$bioreg,at=c(1:6)+.25,add=T,varwidth=T,boxwex=.3,col=cols[2],notch=F,axes=F)
##boxplot(aves.LIT.alpha~aves.lit.cneb$bioreg,at=c(1:6)+.5,add=T,varwidth=T,boxwex=.3,col=cols[3],border=2,notch=F,axes=F)
box()

text(as.numeric(aves.gbif.cneb$bioreg[aves.GBIF.alpha>max.y])+.25,rep(max.y,sum(aves.GBIF.alpha>max.y)),aves.GBIF.alpha[aves.GBIF.alpha>max.y], cex=.7)
arrows(as.numeric(aves.gbif.cneb$bioreg[aves.GBIF.alpha>max.y])+.25, rep(max.y-24,sum(aves.GBIF.alpha>max.y)), as.numeric(aves.gbif.cneb$bioreg[aves.GBIF.alpha>max.y])+.25, rep(max.y-10.2,sum(aves.GBIF.alpha>max.y)),length = 0.03)

par(xpd=T)
symbols(c(1,2),c(1,1),boxplots=matrix(c(.5,2,1,1,.5,.5,2,1,1,.5),ncol=5,byrow=T),bg=cols,fg=1,lwd=1,xlim=c(0.7,2.3),ylim=c(.25,1.95),axes=F,lty=1,xlab="",ylab="")

points(c(1,2),c(1.85,1.85))
arrows(2,1.95,2,2.05,length = 0.03)

text(c(1.5,1.5,1.5,1.5,1.5,1.5),c(0.35,0.75,1,1.25,1.65, 1.90),c("minimum within\n1.5 interquantile\nrange","25%","50%","75%","maximum within\n1.5 interquantile\nrange","outlying observations"),cex=.75)

text(c(1,2),c(.19,.19),c("NeoMapas","Other\nsources"))

title(xlab=NA, ylab="Number of species per cell", outer=T,line=1)


###################################################
### chunk number 93: PrediccionEspacialConteos
###################################################
#line 3033 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.ss <- jmp.NM$valido=="SI"
tt <- aggregate(jmp.NM$jmp[mi.ss],by=list(fase=jmp.NM$fase[mi.ss], CN=jmp.NM$CN[mi.ss], spp=jmp.NM$spp[mi.ss]),luq)
colnames(tt)[4] <- "njmp"

t2 <- aggregate(tvn.NM$sfrz,by=list(fase=tvn.NM$fase, CN=tvn.NM$CN),sum,na.rm=T)
colnames(t2)[3] <- "sfrz"

t3 <- merge(tt,t2,by=c("fase","CN"),all.x=T)
t3$tjmp <- t3$njmp*60/t3$sfrz
t3[,c("x","y")] <- coordinates(CNEB)[match(t3$CN,CNEB@data$cdg),]
t2[,c("x","y")] <- coordinates(CNEB)[match(t2$CN,CNEB@data$cdg),]

spp="Anartia jatrophae"
spp="Anartia amathea"
spp <- "Euptoieta hegesia"
spp <- "Dryadula phaetusa"
spp <- "Phoebis sennae"

spp <- "Eurema phiale"
spp="Anartia amathea"

layout(matrix(1:2,ncol=1))
par(mar=c(0,0,0,0))
##for (spp in spps) {
	mi.ss <- t3$spp==spp
	plot(vzla,border="maroon",ylim=c(0,13))
##	title(main=spp,line=-3)
	symbols(t3$x[mi.ss],t3$y[mi.ss], circle=log1p(t3$tjmp[mi.ss]), inches=.1,add=T, bg=t3$fase[mi.ss])
##}

pp <- t3[mi.ss,]
aa <- t2[!(paste(t2$fase,t2$CN) %in% paste(t3$fase[mi.ss],t3$CN[mi.ss])),]
aa$njmp <- 0

mi.dt <- rbind(pp[,colnames(aa)],aa)
mi.dt$fase <- as.factor(mi.dt$fase)

nw.dt <- data.frame(coordinates(CNEB),fase=factor(2,levels=levels(mi.dt$fase)))
colnames(nw.dt) <- c("x","y","fase")

mi.gam <- gam(njmp~s(x,y)+fase, weights=sfrz/1000, data=mi.dt, subset=mi.dt$sfrz>100, family=poisson(log))

CNEB@data$abnd <- round(predict(mi.gam, nw.dt,type="response"),3)

grps <-cut(CNEB@data$abnd,breaks=c(-1,quantile(CNEB@data$abnd[CNEB@data$abnd>0], p=(1:9)/9)))
cols <- c(NA,brewer.pal(9, "YlOrRd"))
plot(CNEB,col=cols[as.numeric(grps)],border=cols[as.numeric(grps)])
plot(vzla,add=T)

symbols(t3$x[mi.ss],t3$y[mi.ss], squares=log1p(t3$tjmp[mi.ss]), inches=.2,add=T)
##symbols(aa$x,aa$y, circles=log1p(aa$sfrz), inches=.07,add=T)


###################################################
### chunk number 94: 
###################################################
#line 3100 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

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
### chunk number 95: 
###################################################
#line 3124 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
mi.mtz <- scrb.NM
mi.d1 <- vegdist(mi.mtz,"chao")
mi.h1 <- as.dendrogram(as.hclust(agnes(mi.d1)), hang=.2)

mis.clrs <- c(1,"brown")[1+(scrb.rsm$p.id<.1)]
mis.pchs <- c(24,21)[1+(scrb.rsm$cb>50)]
names(mis.clrs) <- paste("NM",scrb.rsm$NM, " :: ",scrb.rsm$yr,sep="")
names(mis.pchs) <- paste("NM",scrb.rsm$NM, " :: ",scrb.rsm$yr,sep="")

 colLab <- function(n) {
           if(is.leaf(n)) {
             a <- attributes(n)
             attr(n, "nodePar") <- c(a$nodePar,
                                     list(lab.font=1,
                                          col=1,bg=mis.clrs[a$label],
                                          lab.col=mis.clrs[a$label],
                                          pch=mis.pchs[a$label],cex=1.35,
                                          lab.cex=.88))

           }
           n
     }
     h0e <- dendrapply(mi.h1, colLab)

par(mar=c(5,4,1,1))
plot(h0e,ylab="Índice de disimilitud (Chao)")


###################################################
### chunk number 96: Todos los AICS
###################################################
#line 3159 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
aics.mrps$grp <- "Pieridae"
aics.scrb$grp <- "Scarabaeinae"
aics.aves$grp <- "Aves"

tabla.aics <- rbind(aics.mrps,aics.scrb,aics.aves)
tabla.aics$R2.adj[tabla.aics$R2.adj<0] <- 0
tabla.aics <- tabla.aics[tabla.aics$vd!="abnd",c("grp", "vd", "formula", "k", "n", "logLik", "AICc", "aic.weights", "R2.adj")]
tabla.aics <- tabla.aics[order(tabla.aics$grp,tabla.aics$vd,tabla.aics$AICc),]

print(xtable(tabla.aics[tabla.aics$aic.weights>.1,],digits=c(0,0,0,0,0,0,2,2,3,2),caption="Comparación de pesos de Akaike para todos los modelos ajustados a los tres grupos taxonómicos"), include.rownames=F, caption.placement="top", size="footnotesize")




###################################################
### chunk number 97: 
###################################################
#line 3174 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
aggregate(tabla.aics$R2.adj[tabla.aics$aic.weights>.1],list(tabla.aics$grp[tabla.aics$aic.weights>.1],tabla.aics$vd[tabla.aics$aic.weights>.1]),range)


###################################################
### chunk number 98: 
###################################################
#line 3179 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
save(file="~/NeoMapas/Rdata/Anexos.rda", scrb.NM, scrb.rsm, scrb.cneb, d1.scrb, mrps.NM, mrps.cneb, d1.mrps, aves.NM, aves.cneb, d1.aves)


###################################################
### chunk number 99: GenerarKML eval=FALSE
###################################################
## #line 3185 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"
## cat("<?xml version='1.0' encoding='UTF-8'?>
##      <kml xmlns='http://earth.google.com/kml/2.2'>
##      <Document>
## <name>Cuadrícula Nacional para el Estudio de la Biodiversidad en Venezuela (CNEB)</name>
##   <open>1</open>
## <description>Cuadrícula de referencia para la iniciativa NeoMapas 
## http://www.neomapas.org</description>
##   <Style id='transPurpleLineGreenPoly'>
##       <LineStyle>
##         <color>7fff00ff</color>
##         <width>4</width>
##       </LineStyle>
##       <PolyStyle>
##         <color>7f00ff00</color>
##       </PolyStyle>
##     </Style>
## ",      file="CNEB.kml",append=F)
## 
## slc <- unique(c(unique(tvn.NM$NM), unique(scrb.rsm$NM), trans.info$NM))
## slc <- slc[!is.na(slc)]
## 
## cns <- info.NM[info.NM$NM %in% as.numeric(slc),]
## 
## mis.cc <- CNEB.nm[CNEB.nm$cdg %in% cns$CNEB,c("cat","cdg","lon","lat","bioreg","est.pr")]
## 
## for (j in 1:nrow(mis.cc)) {
## dd <- expand.grid(x=mis.cc[j,c("lon")]+c(0.25,-.25),
##                   y=mis.cc[j,c("lat")]+c(-0.25,.25),
##                   alt=1000)
## 
## 	 cat(paste("<Folder>
##         <name>",mis.cc$cdg[j],"</name>
##         <description>
##         Celda ",mis.cc$cdg[j]," de la CNEB, en la bioregión ",mis.cc$bioreg[j]," y el estrato ",mis.cc$est.pr[j],"
##        </description>
##        <Placemark>
##         <name>",mis.cc$cdg[j],"</name>
##         <description>
##         Celda ",mis.cc$cdg[j],"
##         </description>
##         <styleUrl>#transPurpleLineGreenPoly</styleUrl>
##        <Polygon>
##             <extrude>1</extrude>
##             <altitudeMode>relativeToGround</altitudeMode>
##             <outerBoundaryIs>
##             <LinearRing>
##             <coordinates>\n",
##        paste(apply(dd[c(1,2,4,3,1),],1,paste,collapse=","),collapse="\n")     
##             ,"\n</coordinates>
##             </LinearRing>
##             </outerBoundaryIs>
##             <innerBoundaryIs>
##             <LinearRing>
##             <coordinates>\n",
##        paste(apply(dd[c(1,2,4,3,1),],1,paste,collapse=","),collapse="\n")     
##             ,"\n</coordinates>
##             </LinearRing>
##             </innerBoundaryIs>
##         </Polygon>
##       </Placemark>",sep=""),
##       file="CNEB.kml",append=T)
## 
## 
## 
## 	for (nm in cns$NM[cns$CNEB ==mis.cc$cdg[j]]) {
## 		nmbr <- cns[cns$NM==nm,"Nombre"]
## 		adm1 <- cns[cns$NM==nm,"ADM1"]
## 
## 	 cat(paste("<Folder>
##         <name>NM",nm," :: ",nmbr,"</name>
##         <description>
##         Transecta NM",nm," en el estado ",adm1,"
##        </description>",sep=""),
##         file="CNEB.kml",append=T)
## 
## 		mi.tvn <- tvn.NM[tvn.NM$NM %in% nm,]
## 		mi.trmp <- trmp.NM[trmp.NM$NM %in% nm,]
## 		mi.avs <- trans.info[trans.info$NM==nm,]
## 
## 		if (nrow(mi.tvn)>0) {
## 			t1 <- aggregate(mi.tvn$sfrz,by=list(mi.tvn$yr),sum,na.rm=T)
## 			t2 <- aggregate(mi.tvn$vst,by=list(mi.tvn$yr),luq)
## 
## 			for (k in t1$Group.1) {
## 
## 			dscr <- paste("Año ",k,": ",t2[t2$Group.1==k,2]," muestreos para un esfuezo total de ",round(t1[t1$Group.1==k,2]/60,2)," horas*colector<br/>",sep="")
## 			
## 			t3 <- aggregate(mi.tvn[mi.tvn$yr==k,c("lon","lat")],by=list(mi.tvn$vst[mi.tvn$yr==k]),median,na.rm=T)
## 
## 		   cat(paste("<Folder>
##     	    <name> NM",nm," año ",k," Muestreos de mariposas</name>
##         	<description><![CDATA[
## 	        ",dscr,"
##     	    ]]> </description>",
## 			paste("<Placemark><name>",t3$Group.1,"</name><Point><coordinates>", t3$lon, ",", t3$lat, "\n</coordinates></Point></Placemark>",sep="",collapse="\n"),
##       		"",sep=""),
##         	file="CNEB.kml",append=T)
## 			cat(paste("</Folder>",sep=""),
##     	    file="CNEB.kml",append=T)
##     			
## 			}
## 		}
## 
## 		if (nrow(mi.avs)>0) {
## 			t3 <- pps[grep(mi.avs$IDTransecta,pps$V3),]
## 			colnames(t3) <- c("lat","lon","Group.1")
## 
## 			cat(paste("<Folder>
##     	    <name> NM",nm," año 2010 Muestreos de aves</name>
##         	<description><![CDATA[
## 	        protocolo de observación de aves de NeoMapas aplicado en 2010
##     	    ]]> </description>",
## 			paste("<Placemark><name>",t3$Group.1,"</name><Point><coordinates>", t3$lon, ",", t3$lat, "\n</coordinates></Point></Placemark>",sep="",collapse="\n"),
##       		"",sep=""),
##         	file="CNEB.kml",append=T)
## 			cat(paste("</Folder>",sep=""),
##     	    file="CNEB.kml",append=T)
## 		}		
## 
## 		if (nrow(mi.trmp)>0) {
## 			t1 <- aggregate(mi.trmp$sfrz,by=list(mi.trmp$yr),sum,na.rm=T)
## 			t2 <- aggregate(mi.trmp$vst,by=list(mi.trmp$yr),luq)
## 
## 			for (k in t1$Group.1) {
## 
## 			dscr <- paste("Año ",k,": ",t2[t2$Group.1==k,2]," trampas colocadas para un esfuezo total de ",round(t1[t1$Group.1==k,2]/24,2)," dias*trampa<br/>",sep="")
## 			
## 			t3 <- aggregate(mi.trmp[mi.trmp$yr==k,c("lon","lat")],by=list(mi.trmp$vst[mi.trmp$yr==k]),median,na.rm=T)
## 
## 		   cat(paste("<Folder>
##     	    <name> NM",nm," año ",k," Muestreos de escarabajos</name>
##         	<description><![CDATA[
## 	        ",dscr,"
##     	    ]]> </description>",
## 			paste("<Placemark><name>",t3$Group.1,"</name><Point><coordinates>", t3$lon, ",", t3$lat, "\n</coordinates></Point></Placemark>",sep="",collapse="\n"),
##       		"",sep=""),
##         	file="CNEB.kml",append=T)
## 			cat(paste("</Folder>",sep=""),
##     	    file="CNEB.kml",append=T)
## 			}
## 		}
## 
##      cat(paste("</Folder>",sep=""),
##         file="CNEB.kml",append=T)
## 
##   }
##      cat(paste("</Folder>",sep=""),
##         file="CNEB.kml",append=T)
##  
## }
## cat("</Document></kml>",
##       file="CNEB.kml",append=T)
## 
## 
## 


###################################################
### chunk number 100: GenerarKML
###################################################
#line 3343 "~/Dropbox/NeoMapas/doc/200_InformeNeoMapas2010/Documento2_InformeResultadosNeoMapas.Rnw"

arch.kml <- "20110626_NeoMapas_version0.3_preliminar.kml"

cat("<?xml version='1.0' encoding='UTF-8'?>
     <kml xmlns='http://earth.google.com/kml/2.2'>
     <Document>
<name>Muestreos de NeoMapas 2005 al 2010</name>
  <open>1</open>
<description>Ubicación de las localidades de muestreo de la
<b>Iniciativa para el Mapeo de la Biodiversidad Neotropical</b>
http://www.neomapas.org
Versión 0.3 (preliminar, puede contener errores en las coordenadas, por favor no difundir)
</description>
",      file=arch.kml,append=F)

slc <- unique(c(unique(tvn.NM$NM), unique(scrb.rsm$NM), trans.info$NM))
slc <- slc[!is.na(slc)]

cns <- info.NM[info.NM$NM %in% as.numeric(slc),]

	for (nm in cns$NM) {
		nmbr <- cns[cns$NM==nm,"Nombre"]
		adm1 <- cns[cns$NM==nm,"ADM1"]
		bioreg <- cns[cns$NM==nm,"Bioregión"]
		estrato <- cns[cns$NM==nm,"Estrato"]

	 cat(paste("<Folder>
        <name>NM",nm," :: ",nmbr,"</name>
        <description>
        Transecta NM",nm," en el estado ",adm1,"<br/>
        Ubicada en la bioregión ",bioreg," y el estrato ",estrato,"<br/>
       </description>",sep=""),
        file=arch.kml,append=T)

		mi.tvn <- tvn.NM[as.numeric(tvn.NM$NM) %in% nm,]
		mi.trmp <- trmp.NM[as.numeric(trmp.NM$NM) %in% nm,]
		mi.avs <- trans.info[as.numeric(trans.info$NM)==nm,]

		if (nrow(mi.tvn)>0) {
			t1 <- aggregate(mi.tvn$sfrz,by=list(mi.tvn$yr),sum,na.rm=T)
			t2 <- aggregate(mi.tvn$vst,by=list(mi.tvn$yr),luq)

			for (k in t1$Group.1) {

			dscr <- paste("Año ",k,": ",t2[t2$Group.1==k,2]," muestreos para un esfuezo total de ",round(t1[t1$Group.1==k,2]/60,2)," horas*colector<br/>",sep="")
			
			t3 <- aggregate(mi.tvn[mi.tvn$yr==k,c("lon","lat")],by=list(mi.tvn$vst[mi.tvn$yr==k]),median,na.rm=T)

		   cat(paste("<Folder>
    	    <name> NM",nm," año ",k," Muestreos de mariposas</name>
        	<description><![CDATA[
	        ",dscr,"
    	    ]]> </description>",
			paste("<Placemark><name>",t3$Group.1,"</name><Point><coordinates>", t3$lon, ",", t3$lat, "\n</coordinates></Point></Placemark>",sep="",collapse="\n"),
      		"",sep=""),
        	file=arch.kml,append=T)
			cat(paste("</Folder>",sep=""),
    	    file=arch.kml,append=T)
    			
			}
		}

		if (nrow(mi.avs)>0) {
			t3 <- pps[grep(mi.avs$IDTransecta,pps$V3),]
			colnames(t3) <- c("lat","lon","Group.1")

			cat(paste("<Folder>
    	    <name> NM",nm," año 2010 Muestreos de aves</name>
        	<description><![CDATA[
	        protocolo de observación de aves de NeoMapas aplicado en 2010
    	    ]]> </description>",
			paste("<Placemark><name>",t3$Group.1,"</name><Point><coordinates>", t3$lon, ",", t3$lat, "\n</coordinates></Point></Placemark>",sep="",collapse="\n"),
      		"",sep=""),
        	file=arch.kml,append=T)
			cat(paste("</Folder>",sep=""),
    	    file=arch.kml,append=T)
		}		

		if (nrow(mi.trmp)>0) {
			t1 <- aggregate(mi.trmp$sfrz,by=list(mi.trmp$yr),sum,na.rm=T)
			t2 <- aggregate(mi.trmp$vst,by=list(mi.trmp$yr),luq)

			for (k in t1$Group.1) {

			dscr <- paste("Año ",k,": ",t2[t2$Group.1==k,2]," trampas colocadas para un esfuezo total de ",round(t1[t1$Group.1==k,2]/24,2)," dias*trampa<br/>",sep="")
			
			t3 <- aggregate(mi.trmp[mi.trmp$yr==k,c("lon","lat")],by=list(mi.trmp$vst[mi.trmp$yr==k]),median,na.rm=T)

		   cat(paste("<Folder>
    	    <name> NM",nm," año ",k," Muestreos de escarabajos</name>
        	<description><![CDATA[
	        ",dscr,"
    	    ]]> </description>",
			paste("<Placemark><name>",t3$Group.1,"</name><Point><coordinates>", t3$lon, ",", t3$lat, "\n</coordinates></Point></Placemark>",sep="",collapse="\n"),
      		"",sep=""),
        	file=arch.kml,append=T)
			cat(paste("</Folder>",sep=""),
    	    file=arch.kml,append=T)
			}
		}

     cat(paste("</Folder>",sep=""),
        file=arch.kml,append=T)

}
cat("</Document></kml>",
      file=arch.kml,append=T)


