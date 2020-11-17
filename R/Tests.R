cat("Reading Tests.R\n")

# Compare the two computations of the function 
# 
#       H = H(j) := P( Y(t_{j+1})>=d_{j+1} | Y(t_j)>=d_j )
#         = fcn_H(j) 
#         = fcn_H1(j)
#
# from docs/BrownianMotionProblem.pdf, p1, eq(2),
# for all j=1,...,nSteps-1.
# Report any discrepancies.
#
test_H_function <- function(){
   
   outFile <- "results/H_function_Test.txt"
   msg <-paste("\n\nChecking if fcnH(j)=fcn_H1(j) for all j in 1:",nSteps-1,"\n",sep="")
   cat(msg); cat(msg,file=outFile)
   OK <- TRUE
   deviation <- double(0)
   fH <- double(0)
   fH1 <- double(0)
   step <- double(0)
   for(j in 1:(nSteps-1)){
      
      f_H <- fcn_H(j)
      f_H1 <- fcn_H1(j)
      if(abs(f_H-f_H1)/f_H1 > 1e-8){
         
         OK <- FALSE
         step <- c(step,j)
         fH <- c(fH,f_H)
         fH1 <- c(fH1,f_H1)
         deviation <- c(deviation,100*abs(f_H-f_H1)/f_H1)
      }
   }

   if(OK){
      
     cat("No deviation between fcn_H and fcn_H1 detected.\n", file=outFile, append=TRUE) 
      
   } else {
      
     result <- data.frame(
        j = step,
        fcn_H_at_j = fH,
        fcn_H1_at_j = fH,
        deviationPcnt =deviation
     ) 
     cat("Deviations between fcn_H and fcn_H1 detected:\n\n", file=outFile, append=TRUE)
     capture.output(print(result), file=outFile, append=TRUE)
   }
   cat("Finished, results in",outFile,"\n")
}