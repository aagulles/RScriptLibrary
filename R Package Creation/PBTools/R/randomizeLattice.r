#-------------------------------------------------
# Author: Alaine A. Gulles
#-------------------------------------------------

randomizeLattice <- function(k, r) UseMethod("randomizeLattice")

randomizeLattice.default <- function(k, r) {
     
     # -- checking input by the user:
     if (!is.numeric(k) || length(k) != 1) { stop("The object 'k' should be numeric of length 1.") }
     if (!is.numeric(r) || length(r) != 1) { stop("The object 'r' should be numeric of length 1.") }
     if (k > 12 || k < 3) { stop("This function can only construct lattice design for k between 3 and 12.") }
     if (k %% 2 != 0 || k == 4) {
          if (r < 1 || r > k + 1) { stop("The object 'r' should be between 2 and k + 1.") }
          repMax <- k + 1
     } else { 
          if (k == 8) { repMax <- 6 } else { repMax <- 3 }
          if (r < 1 || r > repMax) { stop(paste("The object 'r' should be between 2 and ", repMax,".", sep = "")) }
     } 
     
     
     # k is the block size (number of treatment per block) 
     ntrmt <- k**2   	# treatment levels
     b <- ntrmt/k	# number of blocks per rep
     #repMax <- k + 1
     
     
     #if (ntrmt %% k != 0) { stop("The size of the block is not appropriate. The number of treatments must be multiple of k (block size).") }
     if (r %% repMax == 0) { designType <- "Balanced"
     } else { designType <- "Partially Balanced" }
     
     tempLayout <- list()
     
     trmtRandomize <- sample(paste("T", 1:ntrmt, sep = ""), ntrmt, replace = FALSE)
     #trmtRandomize <- paste("T", 1:ntrmt, sep = "")
     tempLayout[[1]] <- matrix(trmtRandomize, nrow = b, ncol = k, byrow = TRUE)
     tempLayout[[2]] <- t(tempLayout[[1]])
     
     if (r > 2) {
          # CASE: K is odd
          if (k %% 2 != 0) {
               temp <- tempLayout[[2]]
               for (i in (3:r)) {
                    temp <- cbind(temp, temp[,1:(k-1)])
                    tempLayout[[i]] <- matrix(0, nrow = b, ncol = k)
                    for (j in (1:k)) { tempLayout[[i]][j,] <- diag(temp[, j:(k+j-1)]) }
                    temp <- tempLayout[[i]]
               }
          } else {
               if (k == 4) {
                    temp1 <- matrix(c(1,2,2,1),nrow = 2, ncol = 2)
                    temp2 <- temp1 + 2
                    temp <- rbind(cbind(temp1, temp2), cbind(temp2, temp1))
                    for (i in (3:r)) {
                         tempLayout[[i]] <- matrix(tempLayout[[i-1]][order(temp)], nrow = b, ncol = k, byrow = TRUE)
                         templine <- temp[2,]
                         for (j in (2:nrow(temp))) {
                              if (j < nrow(temp)) { temp[j,] <- temp[j+1,]
                              } else { temp[j,] <- templine }
                         }
                    }
               } # end stmt -- if (k == 4)
               
               if (k == 8) {
                    temp <- matrix(c(1,2,3,4,5,6,7,8,2,1,7,6,8,4,3,5,3,7,1,5,4,8,2,6,4,6,5,1,3,2,8,7,5,8,4,3,1,7,6,2,6,4,8,2,7,1,5,3,7,3,2,8,6,5,1,4,8,5,6,7,2,3,4,1), nrow = b, ncol = k, byrow = TRUE)
                    tempLayout[[3]] <- matrix(tempLayout[[1]][order(temp)], nrow = b, ncol = k, byrow = TRUE)
                    if (r > 4) {
                         for (i in (4:r)) {
                              if (i == 4) { index <- c(1,7,5,2,8,3,6,4) } else { index <- c(1,3,4,6,2,8,5,7) }
                              if (i != 6) { 
                                   temp <- temp[, index] 
                                   tempLayout[[i]] <- matrix(tempLayout[[i - 1]][order(temp)], nrow = b, ncol = k, byrow = TRUE)
                              } else { tempLayout[[i]] <- matrix(tempLayout[[i - 1]][order(t(temp))], nrow = b, ncol = k, byrow = TRUE) }
                         }
                    }
               } # end stmt -- if (k == 8)
               
               if (k == 6 || k == 10 || k == 12) {
                    tempLayout[[3]] <- matrix(0, nrow = b, ncol = k) 
                    temp <- rbind(tempLayout[[2]], tempLayout[[2]])
                    for (j in (1:b)) { tempLayout[[3]][j,] <- diag(temp[j:(b+j-1),]) }
                    if (k == 12 && r >= 4) {
                         tempLayout[[4]] <- matrix(0, nrow = b, ncol = k) 
                         temp <- rbind(tempLayout[[3]], tempLayout[[3]])
                         for (j in (1:b)) { tempLayout[[4]][j,] <- diag(temp[j:(b+j-1),]) }
                    }
               } # end stmt -- if (k == 6 || k == 10 || k == 12)
          } # end stmt -- else if (k %% 2 != 0)
     }
     
     names(tempLayout) <- paste("Rep", 1:length(tempLayout), sep = "")
     book <- NULL
     for (i in (1:length(tempLayout))){
          tempLayout[[i]] <- tempLayout[[i]][sample(1:b, b, replace = FALSE),]
          rownames(tempLayout[[i]]) <- paste("Block", ((b*i)-(b-1)):(b * i), sep = "")	
          colnames(tempLayout[[i]]) <- paste(1:ncol(tempLayout[[i]]))
          book <- rbind(book, data.frame(Rep = i, as.data.frame.table(tempLayout[[i]])))
     }
     book <- book[,c("Rep","Var1","Freq")]
     book[,"Var1"] <- as.numeric(book[,"Var1"])
     names(book) <- c("Rep", "Block", "Treatment")
     book[,"Rep"] <- factor(book[,"Rep"])
     book[,"Block"] <- factor(book[,"Block"])
     book[,"Treatment"] <- factor(book[,"Treatment"])
     book <- book[order(book$Rep, book$Block),]
     
     return(list(design = paste(k, "x", k, designType, "Lattice Design"), 
                 plan = tempLayout, 
                 randomization = book))
}