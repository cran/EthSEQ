get.Model <- function(model.available,model.folder)
{
  if(!model.available%in%c("SS2","SS2.Light","SS4","HALO","NimblegenV3","Exonic"))
  {
    message.Date(paste("Model ",model.available," not available",sep=""))
    return(NA)
  }
  if(!file.exists(file.path(model.folder,paste(model.available,".Model.gds",sep=""))))
  {
    download.file(paste("https://github.com/aromanel/EthSEQ_Data/raw/master/EthSEQ_Models/",model.available,".Model.gds",sep=""),
                  file.path(model.folder,paste(model.available,".Model.gds",sep="")))
  }
  return(file.path(model.folder,paste(model.available,".Model.gds",sep="")))
}