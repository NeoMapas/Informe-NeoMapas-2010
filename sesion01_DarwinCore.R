### Generar Darwin Core


#############
## EML file
############
## documentado en:
##https://knb.ecoinformatics.org/#external//emlparser/docs/eml-2.1.1/index.html#N100DE
## interfase en R
##https://github.com/ropensci/EML

##geographic content
 UM <- CNEB[CNEB.nm$UM>0,]
UM@bbox
## temporal range, format: http://www.w3.org/TR/NOTE-datetime
range(trans.info$Fecha.Llegada)
range(trans.info$Fecha.Partida)


#############
## Simple Darwin Core
############
## documentado en http://rs.tdwg.org/dwc/terms/simple/index.htm
## http://rs.tdwg.org/dwc/terms/#dcterms:type
## http://dublincore.org/documents/dcmi-terms/#terms-type
## https://code.google.com/p/darwincore/wiki/Examples
## https://code.google.com/p/darwincore/wiki/Occurrence
##http://rs.tdwg.org/dwc/terms/#eventTime

prs <- read.csv("~/NeoMapas/data/NMAves2010/Observadores.total.csv",as.is=T)
prs$Observador[prs$Codigobservador=="JRF"] <- "Jose R. Ferrer Paris"

DwC.data <- data.frame(oID = sprintf("NMAVS2010:M1:OBS%05d", NM.m1$ID),
                       eID = sprintf("NMAVS2010:M1:VST%s", NM.m1$idpunto),
                       lID = sprintf("NMAVSLOC:%s", NM.m1$idpunto),
                       iC = NM.m1$Ntotal,
                       eD = sprintf("20%s-%s-%s",substr(trans.info$Fecha.Muestreo.1[match(NM.m1$IDTransecta,trans.info$IDTransecta)],7,8),substr(trans.info$Fecha.Muestreo.1[match(NM.m1$IDTransecta,trans.info$IDTransecta)],0,2),substr(trans.info$Fecha.Muestreo.1[match(NM.m1$IDTransecta,trans.info$IDTransecta)],4,5)),
                       eT = sprintf("%s-0430",substr(NM.m1$Hora,10,14)),
                       lat=NM.m1$lat,
                       lon=NM.m1$lon,
                       nmbr=spp.aves$Latname[match(NM.m1$Especieid,spp.aves$N_ave)],
                       epiteto=sub("^[A-za-z]+ ","",spp.aves$Latname[match(NM.m1$Especieid,spp.aves$N_ave)]),
                       genero=sub(" [A-za-z]+$","",spp.aves$Latname[match(NM.m1$Especieid,spp.aves$N_ave)]),
                       rB=prs$Observador[match(NM.m1$Codigobservador,prs$Codigobservador)],
                       oR=paste(ifelse(NM.m1$Noidos==0,"Visual ID",
                         ifelse(NM.m1$Noidos==NM.m1$Ntotal,
                                "Auditory ID","Both visual and auditory ID")),
                         ifelse(NM.m1$Volando==0,"",", flying overhead"),
                         ifelse(NM.m1$Observaciones=="",".",
                                sprintf(". Aditional notes (spanish): %s",
                                        NM.m1$Observaciones)),sep=""),
                       iB=prs$Observador[match(as.numeric(substr(NM.m1$IDTransecta,0,1)),prs$Region)],
                       iVS=ifelse(NM.m1$Seguro==1,"Certain",ifelse(NM.m1$Seguridad>90,"Almost certain (>90% confidence)",ifelse(NM.m1$Seguridad>50,"Probably  (>50% confidence)","Doubtful (<50% confidence)"))))


DwC.dt <- subset(DwC.data,!NM.m1$Especieid %in% 1379)

vcblry <- c("dcterms","dcterms","dwc","dwc","dcterms",rep("dwc",28))
vcblry <- c("http://purl.org/dc/terms/","http://purl.org/dc/terms/","http://rs.tdwg.org/dwc/terms/","http://rs.tdwg.org/dwc/terms/","http://purl.org/dc/terms/",rep("http://rs.tdwg.org/dwc/terms/",28))

  
  
DwC <- with(DwC.dt,data.frame(SDR=oID,language="en",
                              type="Event",
                              basisOfRecord="HumanObservation",
                              occurrenceID=oID,
                              modified="2014-06-09T12:00:00Z",
                              institutionCode="IVIC",
                              collectionCode="NMaves",
                              individualCount=iC,
                              eventID=eID,
                              samplingProtocol="point count",
                              samplingEffort="3 minute * observer",
                              eventDate=eD,
                              eventTime=eT,
                              locationID=lID,
                              country="Venezuela",
                              countryCode="VE",
                              decimalLatitude=lat,
                              decimalLongitude=lon,
                              geodeticDatum="WGS84",
                              georeferenceProtocol="GPS reading",
                              coordinateUncertaintyInMeters="15",
                              scientificName=nmbr,
                              class="Aves",
                              genus=genero,
                              specificEpithet=epiteto,
                              identifiedBy=iB,
                              dateIdentified=eD,
                              nameAccordingTo="Hilty, S. L. 2003. Birds of Venezuela. Second Edition. Princeton University Press, Princeton, New Jersey, USA,",
                              recordedBy=rB,
                              occurrenceStatus="present",
                              occurrenceRemarks=oR,
                              identificationReferences="Hilty, S. L. 2003. Birds of Venezuela. Second Edition. Princeton University Press, Princeton, New Jersey, USA,",
                              identificationVerificationStatus=iVS))

 cat(sprintf("<field  index='%s' term='%s%s'/>\n",1:33,vcblry,colnames(DwC)[-1]))


write.table(file="~/NeoMapas/DwC/NMavesM1/NMaves2010M1.tab",
            sep = "\t",
            DwC,row.names=F)


######################
###
#####################


##excluir estos
table(NM.m1$Especieid %in% 1379)

#############
## Core
###############


#############
## extensions
############
Event

eventID | samplingProtocol | samplingEffort | eventDate | eventTime | startDayOfYear | endDayOfYear | year | month | day | verbatimEventDate | habitat | fieldNumber | fieldNotes | eventRemarks





c("GUID", "Title", "Continent/Ocean", "Verbatim Coordinate System", "Coordinate uncertainty in meters", "Country (ISO alpha-2)", "County", "Geodetic datum", "Georeference protocol", "Georeference remarks", "Island group", "Island", "Locality", "Map", "Maximum depth", "Maximum elevation", "Minimum depth", "Minimum elevation", "State/Province")


require(RMySQL)

drv <- dbDriver("MySQL")
cn1 <- dbConnect(drv,user="NMerre",password="datos",dbname="NM_Escarabajos")
tmp <- dbGetQuery(cn1,"SET NAMES UTF8")
tmp <- dbGetQuery(cn1,"SET CHARACTER SET UTF8")


## localidades tienen un formato de nombre, y las trampas otro... revisar

rslt <- dbGetQuery(cn1,"select * from LOCALIDADES where NM='02'")
jmps <- dbGetQuery(cn1,"select * from EJEMPLARES where cdg_trampa like '02-%'")
table(jmps$cdg_trampa)

##select cdg_localidad, CONCATENATE("Paraguaná (NM",NM,", VBG ",CN,") ",Prd), "Southern America", "decimal degrees", epe, "VE", municipio, "WGS84", "GPS reading", "","","", CONCATENATE(nombre,, ", Paraguaná, estado ",estado), CONCATENATE("POINT:(",lat,",",lon,")","",alt,"",alt,estado) from LOCALIDADES where NM="02"


write.csv(file="PruebaDCloc.csv",
          unique(with(rslt,
                      data.frame(GUID=cdg_localidad,
                                 Title=sprintf("Paraguaná (NM%s, VBG%s) %s",NM,CN,Prd),
                                 Continent="Southern America",
                                 Coords="decimal degrees",
                                 Coord.un=epe,
                                 Country="VE",
                                 County=municipio,
                                 Datum="WGS84",
                                 Protocol="GPS reading",
                                 remarks="",
                                 Igroup="",
                                 Island="",
                                 Locality=sprintf("%s, Paraguaná",nombre),
                                 Map=sprintf("POINT:(%s,%s)",lat,lon),
                                 mxdepth="",
                                 mxalt=alt,
                                 mndepth="",
                                 mnalt=alt,
                                 state="Falcón"))))
