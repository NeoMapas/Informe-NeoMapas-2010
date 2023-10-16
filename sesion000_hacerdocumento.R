##export TEXINPUTS="/Library/Frameworks/R.framework/Versions/Current/ Resources/share/texmf:"
##R --vanilla
setwd("~/NeoMapas/tmp")
hoy <- format(Sys.time(), "%Y%m%d")
options(width=70)
mi.dir <- "200_InformeNeoMapas2010"
mi.arch <- "Documento1_InformeActividadesNeoMapas"
titulo <- "InformeNeoMapas2010"
 options(device.ask.default=FALSE)

system(paste("mkdir ~/NeoMapas/img/Figs/",mi.dir,sep=""))


Sweave(file=paste("~/NeoMapas/doc/",mi.dir,"/",mi.arch,".Rnw",sep=""),eps=F)
Stangle(file=paste("~/NeoMapas/doc/",mi.dir,"/",mi.arch,".Rnw",sep=""))
tools::texi2dvi(paste(mi.arch,".tex",sep=""), pdf=TRUE)

system(paste("mv ",mi.arch,".pdf ~/NeoMapas/doc/",mi.dir,"/",hoy,"_",titulo,".pdf",sep=""))
system(paste("mv ",mi.arch,".R ~/NeoMapas/doc/",mi.dir,"/",hoy,"_",titulo,".R",sep=""))


system(paste("tar -cjvf ",hoy,"_",titulo,"_PDFS.tar.bz2 ",mi.arch,"*pdf",sep=""))
##system(paste("zip ",hoy,"_",titulo,"_JPGS.zip ",mi.arch,"*jpg",sep=""))

system(paste("mv ",hoy,"_",titulo,"_PDFS.tar.bz2 ~/NeoMapas/doc/",mi.dir,"/",sep=""))


system("rm *jpg tmp.ppm")
for (jj in dir(pattern=paste(mi.arch,"-[0-9A-Za-z]+.pdf",sep=""))) {
  system(paste("pdftoppm ",jj," > tmp.ppm",sep=""))
  system(paste("ppmtojpeg --quality=70 tmp.ppm > ",sub(".pdf",".jpg",jj),sep=""))
## system(paste("ppmtojpeg --quality=70 --greyscale tmp.ppm > ",sub(".pdf","gris.jpg",jj),sep=""))
  
}

system(paste("tar -cjvf ",hoy,"_",titulo,"_JPGS.tar.bz2 ",mi.arch,"*jpg",sep=""))
system(paste("mv ",hoy,"_",titulo,"_JPGS.tar.bz2 ~/NeoMapas/",mi.dir,"/",sep=""))
