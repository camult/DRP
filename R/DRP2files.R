#' @title Deregressing estimated breeding values - Two long format file
#' 
#' @description This package is easy to use and can be helpful to calculate deregressed proofs, 
#'              and their reliabilities and weights.
#' 
#' @param animalData It is animal data file
#' @param parentData It is parents data file
#' @param animalCol It is the name of the animal's column
#' @param sireCol It is the name of the animal's dam column
#' @param damCol It is the name of the animal's dam column
#' @param parentCol It the name of the parents' column in the parents data file
#' @param ebvName It is the name of the EBV column
#' @param r2Name It is the name of the accuracy column
#' @param c It is the fraction of genetic variance not explained by markers
#' @param h2 It the heritability of the trait
#' 
#' @return A data frame with deregressed proofs, reliability and weights. 
#' 
#' @references Garrick, D. J., J. F. Taylor, and R. L. Fernando. 2009. 
#'             Deregressing estimated breeding values and weighting information 
#'             for genomic regression analyses. Genet. Sel. Evol. 41:55.
#'
#' @examples
#' ## Not to run ##
#' 
#' ## Example from Garrick et al., (2009)
#' 
#' animalData=data.frame(ID="Animal", sire="Sire", dam="Dam", EBV=15, r2=0.68)
#' parentData=data.frame(ID=c("Sire", "Dam"), EBV=c(10, 2), r2=c(0.97, 0.36))
#'
#' DRP2files(animalData=animalData,
#'           parentData=parentData,
#'           animalCol = "ID",
#'           sireCol   = "sire",
#'           damCol    = "dam",
#'           parentCol = "ID",
#'           ebvName   = "EBV",
#'           r2Name   = "r2",
#'           c         = 0.5,
#'           h2        = 0.25)
#'
#' ## End(Not run)
#' 
#' @export DRP2files
DRP2files <- function(animalData, parentData, animalCol, sireCol, damCol, 
                 parentCol, ebvName, r2Name, c=0.5, h2){
  cat("\n")
  centerText <- function() {
    width <- getOption("width")
    A <- ("                                                       ._____.    \n")
    B <- ("    _/////_            Fernando Brito Lopes           _|_____|_   \n")
    C <- ("   (' o o ')     Animal Scientist (Zootechnician)     (' o o ')   \n")
    D <- ("__ooO_(_)_Ooo_____ Animal Breeding and Genetics _____ooO_(_)_Ooo__\n")
    E <- ("                    e-mail: <camult@gmail.com>                    \n")
    ws <- rep(" ", floor((width - nchar(A))/2))
    cat(ws,A,ws,B,ws,C,ws,D,ws,E,ws,sep = "")
  }
  #-------------------------------------------------------------------------------------------------#
  # Garrick et al., (2009) - De-regressed EPDs
  #-------------------------------------------------------------------------------------------------#
  # Personal comunication (Garrick and Taylor, 2016)
  #  If the weight is negative, it means the accuracy of the individual EBV is less than the accuracy
  #  of its parent average EBV. In theory this cannot occur!
  #  In practice, this can occur when there is something wrong with the way in which the
  #  reliabilities were calculated.
  #  Most software approximates the reliabilities, sometimes not very well, and the storage and
  #  reporting of reliabilities may also involve rounding procedures.
  #  We usually discard any samples that dont have a weight sufficiently large to indicate the DEBV
  #  was computed with some real information.
  #  If you have good-sized contemporary groups, you might want the weight to be equivalent to say
  #  half what would be achieved if the animal had a single record on itself.
  #  This will mean you dont use sire records that correspond to them having just one or two progeny
  #  Some breeding programs have a habit of reporting accuracies as 0.05 for animals that have
  #  absolutely no individual data, and for which the parents have no data.
  #  So you have to trap these and chuck them out as you go.
  #-------------------------------------------------------------------------------------------------#
  centerText()
  cat("\n\n")
  animalData[c(damCol,sireCol)][animalData[c(damCol,sireCol)]=="0"] <- NA
  animalData <- na.omit(animalData)
  nani <- nrow(animalData)
  DRP  <- matrix(NA, nrow = nani, ncol=1)
  DRP_PA <- matrix(NA, nrow = nani, ncol=1)
  r2   <- matrix(NA, nrow = nani, ncol=1)
  w    <- matrix(NA, nrow = nani, ncol=1)
  ERC  <- matrix(NA, nrow = nani, ncol=1)
  Dam  <- matrix(NA, nrow = nani, ncol=2)
  Sire <- matrix(NA, nrow = nani, ncol=2)
  for(i in 1:nani){
    animalName <- as.character(animalData[i, animalCol])
    sireName   <- as.character(animalData[i, sireCol])
    damName    <- as.character(animalData[i, damCol])
    if((sireName%in%parentData[,parentCol])&(damName%in%parentData[,parentCol])){
      EBV_paSire<- parentData[parentData[,parentCol]==sireName, ebvName]
      EBV_paDam <- parentData[parentData[,parentCol]==damName, ebvName]
      r2_paSire <- parentData[parentData[,parentCol]==sireName, r2Name]
      r2_paDam <- parentData[parentData[,parentCol]==damName, r2Name]
      EBVpa <- (EBV_paSire+EBV_paDam)/2
      r2pa <- (r2_paSire+r2_paDam)/4
      EBV_anim <- animalData[animalData[,animalCol]==animalName,ebvName]
      r2_anim <- animalData[animalData[,animalCol]==animalName,r2Name]
      delta <- (0.5-r2pa)/(1-r2_anim)
      alpha <- 1/(0.5-r2pa)
      lambda<- (1-h2)/h2
      ZpaZpa<- lambda*((0.5*alpha)-4) + (0.5*lambda*sqrt((alpha^2)+(16/delta)))
      ZiZi  <- (delta*ZpaZpa) + 2*lambda*((2*delta)-1)
      yi <- ((-2*lambda*EBVpa)+(ZiZi+(2*lambda))*EBV_anim)
      di <- yi/ZiZi
      diPA <- di+EBVpa
      r2i<- 1-lambda/(ZiZi+lambda)
      wi <- (1-h2)/((c+(1-r2i)/r2i)*h2)
      ERCi <- ((1-h2)/h2) * (r2i/(1-r2i))
      # Deregressed effective record contribution (ERC)
      ERC[i,1] <- ERCi
      DRP[i,1]<- di
      DRP_PA[i,1]<- diPA
      r2[i,1] <- r2i
      w[i,1]  <- wi
      Sire[i,]  <- c(EBV_paSire, r2_paSire)
      Dam[i,]  <- c(EBV_paDam, r2_paDam)
    }
  }
  parents <- data.frame(Sire, Dam)
  colnames(parents) <- c("sireEBV","sire_r2","damEBV","dam_r2")
  colnames(DRP) <- paste0("DRP_", ebvName)
  colnames(DRP_PA) <- paste0("DRP_PA_", ebvName)
  colnames(r2) <- paste0("DRP_", ebvName, "_r2")
  colnames(w) <- paste0("DRP_", ebvName, "_w")
  colnames(ERC) <- paste0("DRP_", ebvName, "_ERC")
  dEBV <- data.frame(r2, w)
  dEBV <- data.frame(animalData, parents, DRP, dEBV, ERC, DRP_PA)
  #-------------------------------------------------------------------------------------------------#
  if(any(wi<0, na.rm=TRUE)){
    centerText2 <- function(){
      cat("\n")
      width <- getOption("width")
      A <- ("Some weights are negative, it means the accuracy of the individual \n")
      B <- ("individual EBV is less than the accuracy of its parent average EBV.\n")
      C <- ("In theory this cannot occur! In practice, this can occur when there\n")
      D <- ("  is something wrong with the way in which the reliabilities were  \n")
      E <- ("calculated. Most software approximates the reliabilities, sometimes\n")
      G <- (" not very well, and the storage and reporting of reliabilities may \n")
      H <- ("                  also involve rounding procedures                 \n")
      ws <- rep(" ", floor((width - nchar(A))/2))
      message(ws,A,ws,B,ws,C,ws,D,ws,E,ws,G,ws,H,ws,sep = "")
    }
    centerText2()
  }
  #-------------------------------------------------------------------------------------------------#
  return(dEBV)
}
