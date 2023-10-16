##R --vanilla
setwd("~/NeoMapas/tmp")
hoy <- format(Sys.time(), "%Y%m%d")
options(width=70)
mi.dir <- "200_InformeNeoMapas2010"
mi.arch <- "Documento1_InformeActividadesNeoMapas"
titulo <- "InformeNeoMapas2010"
 options(device.ask.default=FALSE)

system(paste("mkdir ~/NeoMapas/img/Figs/",mi.dir,sep=""))


Sweave(file=paste("~/NeoMapas/",mi.dir,"/",mi.arch,".Rnw",sep=""),eps=F)
Stangle(file=paste("~/NeoMapas/",mi.dir,"/",mi.arch,".Rnw",sep=""))
tools::texi2dvi(paste(mi.arch,".tex",sep=""), pdf=TRUE)

system(paste("mv ",mi.arch,".pdf ~/NeoMapas/",mi.dir,"/",hoy,"_",titulo,".pdf",sep=""))
system(paste("mv ",mi.arch,".R ~/NeoMapas/",mi.dir,"/",hoy,"_",titulo,".R",sep=""))



