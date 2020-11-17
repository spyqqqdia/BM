source("R/BM.R")
source("R/Tests.R")

doRunSimulation <- TRUE
doRunTests <- FALSE


# source("R/run.R")



if(doRunTests){
   
   test_H_function()
}



if(doRunSimulation){

   nPaths <-2^22    # about 4 million
   runSimulation(nPaths)   
}
