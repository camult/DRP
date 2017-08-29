#' @title Deregressing estimated breeding values - One wide format file
#' 
#' @description This package is easy to use and can be helpful to calculate deregressed proofs, 
#'              and their reliabilities and weights.
#' 
#' @param Data It is the name of the data file
#' @param animalId It is the name of the animal's column
#' @param sireId It is the name of the sire's column
#' @param damId It is the name of the dam's column
#' @param animalEBV It is the name of the animal's EBV column
#' @param sireEBV It is the name of the sire's EBV column
#' @param damEBV It is the name of the dam's EBV column
#' @param animalr2 It is the name of the animal's accuracy column
#' @param sirer2 It is the name of the sire's accuracy column
#' @param damr2 It is the name of the dam's accuracy column
#' @param traitName It the name of the trait
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
#' Dataset=data.frame(animal="A1000", sire="S10", dam="D100", ebv_anim=15, ebv_sire=10, ebv_dam=2,
#'                    r2_anim=0.68, r2_sire=0.97, r2_dam=0.36, trait="Trait", c=0.5, h2=0.25)
#'
#' wideDRP(Data     =  Dataset,
#'         animalId  = "animal",
#'         sireId    = "sire",
#'         damId     = "dam",
#'         animalEBV = "ebv_anim",
#'         sireEBV   = "ebv_sire",
#'         damEBV    = "ebv_dam",
#'         animalr2  = "r2_anim",
#'         sirer2    = "r2_sire",
#'         damr2     = "r2_dam",
#'         traitName = "trait",
#'         c         =  0.5,
#'         h2        =  0.25)
#'
#' ## End(Not run)
#' 
#' @export wideDRP
#' @import stats
wideDRP <- function(Data, animalId, sireId, damId, c=0.5, h2, traitName=NULL,
                 animalEBV, sireEBV, damEBV, animalr2, sirer2, damr2){
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
  if(is.null(traitName)) traitName <- "Trait"
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
  #-------------------------------------------------------------------------------------------------#
  r2_gm=as.matrix((Data[,sirer2] + Data[,damr2])/4)
  colnames(r2_gm)<-c("r2gm")
  alfa=1/(0.5 - r2_gm)
  colnames(alfa)<-c("alfa")
  delta=(0.5 - r2_gm)/(1-Data[,animalr2])
  colnames(delta)<-c("delta")
  alfa_delta=(alfa^2)+(16/delta)
  lambda_star=(1-h2)/h2
  Zlgm_Zgm=lambda_star*(0.5*alfa - 4) + 0.5*lambda_star*sqrt(alfa_delta)
  colnames(Zlgm_Zgm)<-c("Zlgm_Zgm")
  Zli_Zi=delta*Zlgm_Zgm+2*lambda_star*(2*delta-1)
  colnames(Zli_Zi)<-c("Zli_Zi")
  r2i<- 1-lambda_star/(Zli_Zi+lambda_star)
  colnames(r2i) <- paste0("DRP_", traitName, "_r2")
  gm=as.matrix((Data[,sireEBV] + Data[,damEBV])/2)
  yi=-2*lambda_star*gm + (Zli_Zi + 2*lambda_star)*as.matrix(Data[,animalEBV])
  DRP <- yi/Zli_Zi
  DRP_PA <- DRP + gm
  colnames(DRP) <- paste0("DRP_", traitName)
  colnames(DRP_PA) <- paste0("DRP_PA_", traitName)
  wi <- (1-h2)/((c+(1-r2i)/r2i)*h2)
  colnames(wi) <- paste0("DRP_", traitName, "_w")
  # Deregressed effective record contribution (ERC)
  ERC <- ((1-h2)/h2) * (r2i/(1-r2i))
  colnames(ERC) <- paste0("DRP_", traitName, "_ERC")
  # Computing the de-regressed values without removing parent average
  Anim_lambda=(1-h2)/h2
  Anim_Zli_Zi=Anim_lambda/(1-Data[,animalr2])
  Anim_yi= (Anim_Zli_Zi + Anim_lambda)*Data[,animalEBV]
  Anim_DRP <- as.matrix(Anim_yi/Anim_Zli_Zi)
  colnames(Anim_DRP) <- paste0("Anim_DRP_", traitName)
  Anim_r2i <- as.matrix(1-Anim_lambda/(Anim_Zli_Zi+Anim_lambda))
  colnames(Anim_r2i) <- paste0("Anim_DRP_", traitName, "_r2")
  dEBV <- cbind(Data[, c(animalId, sireId, damId, animalEBV, sireEBV, damEBV,
                      animalr2, sirer2, damr2)], DRP, r2i, wi, ERC, DRP_PA, Anim_DRP, Anim_r2i)
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
