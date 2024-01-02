#===============================================================================

#Revising Kruschke code for the preattentiveness project --- Wei
  # 1. Test for one pilot study data

#===============================================================================

#-------------------------------------------------------------------------------
# MACRO VARIABLES


#-------------------------------------------------------------------------------





# Example for Jags-Ymet-Xnom2fac-MnormalHom.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
#Load The data file 

fileNameRoot = "SalaryNormalHom--" 
graphFileType = "eps" 
myDataFrame = read.csv("pilot_5.csv")

# filter the null data out. Those data come from the participants' wrong answer
#myDataFrame = subset(originDdata, ReactionTime != "null")
# Specify the column names in the data file relevant to the analysis:
yName="ReactionTime" 
# x1 should be factor with fewer ElementNums, to plot in single pane:
x1Name="Feature" 
x2Name="ElemNum" 
# Specify desired contrasts.
# Each main-effect contrast is a list of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (Kruschke sets (-1000,1000), 
# while I use NULL):
x1contrasts = list( 
  list( c("Hue") , c("Length") , compVal=0.0 , ROPE=NULL ) ,
  list( c("Length") , c("Shape") , compVal=0.0 , ROPE=NULL ) 
)
x2contrasts = list( 
  list( c("12") , c("48") , compVal=0.0 , ROPE=NULL ) ,
  list( c("12") , c("192") , compVal=0.0 , ROPE=NULL ) ,
  list( c("768") , c("192","12","48") , compVal=0.0 , ROPE=NULL ) 
)
# Each interaction contrast is a list of 2 lists of 2 vectors of level names, 
# a comparison value (typically 0.0), and a ROPE (which could be NULL)::
x1x2contrasts = list( 
  list( list( c("Hue") , c("Shape") ) ,
        list( c("12") , c("48") ) ,
        compVal=0.0 , ROPE=NULL ) ,
  list( list( c("Hue") , c("Shape") ) ,
        list( c("12") , c("192") ) ,
        compVal=0.0 , ROPE=NULL ) ,
  list( list( c("Hue") , c("Length","Shape") ) ,
        list( c("768") , c("192","12","48") ) , 
        compVal=0.0 , ROPE=NULL )
) 

y = myDataFrame[,yName]
print(y)
print(max(y))
print(min(y))

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-Xnom2fac-MnormalHom.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , 
                    yName=yName , x1Name=x1Name , x2Name=x2Name ,
                    numSavedSteps=15000 , thinSteps=5 , 
                    saveName=fileNameRoot )
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) 
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("b0","b1[1]","b2[1]","b1b2[1,1]","ySigma") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        datFrm=myDataFrame , x1Name=x1Name , x2Name=x2Name ,
                        x1contrasts=x1contrasts , 
                        x2contrasts=x2contrasts , 
                        x1x2contrasts=x1x2contrasts ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , x1Name=x1Name , x2Name=x2Name ,
          x1contrasts=x1contrasts , 
          x2contrasts=x2contrasts , 
          x1x2contrasts=x1x2contrasts ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
# Other specific comparisons of cells:
if ( fileNameRoot == "SalaryNormalHom-" ) {
  # THIS x1level minus THAT x1level at AT x2level:
  THISx1 = "Hue"
  THATx1 = "Shape"
  ATx2 = "12"
  THISidx = which(ElementNums(myDataFrame[,x1Name])==THISx1)
  THATidx = which(ElementNums(myDataFrame[,x1Name])==THATx1)
  ATidx   = which(ElementNums(myDataFrame[,x2Name])==ATx2)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] , 
    main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=NULL )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx1,"-",THATx1,"At",ATx2,sep=""),
            type=graphFileType)
  # THIS x1level minus THAT x1level at AT x2level:
  THISx1 = "Hue"
  THATx1 = "Shape"
  ATx2 = "192"
  THISidx = which(ElementNums(myDataFrame[,x1Name])==THISx1)
  THATidx = which(ElementNums(myDataFrame[,x1Name])==THATx1)
  ATidx   = which(ElementNums(myDataFrame[,x2Name])==ATx2)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",THISidx,",",ATidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",THATidx,",",ATidx,"]",sep="")] , 
    main=paste(THISx1,"-",THATx1,"@",ATx2) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=NULL )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx1,"-",THATx1,"At",ATx2,sep=""),
            type=graphFileType)
  # THIS x2level minus THAT x2level at AT x1level:
  THISx2 = "192"
  THATx2 = "48"
  ATx1 = "Hue"
  THISidx = which(ElementNums(myDataFrame[,x2Name])==THISx2)
  THATidx = which(ElementNums(myDataFrame[,x2Name])==THATx2)
  ATidx   = which(ElementNums(myDataFrame[,x1Name])==ATx1)
  openGraph(height=4,width=4)
  compInfo = plotPost( 
    as.matrix(mcmcCoda)[,paste("m[",ATidx,",",THISidx,"]",sep="")] -
      as.matrix(mcmcCoda)[,paste("m[",ATidx,",",THATidx,"]",sep="")] , 
    main=paste(THISx2,"-",THATx2,"@",ATx1) , 
    xlab=paste("Difference in",yName) , 
    compVal=0 ,ROPE=NULL )
  show(compInfo)
  saveGraph(file=paste(fileNameRoot,THISx2,"-",THATx2,"At",ATx1,sep=""),
            type=graphFileType)
}
#------------------------------------------------------------------------------- 
