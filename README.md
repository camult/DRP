# DRP

R package to calculate deregressed proofs, and its reliabilities and weights.


# How to Install

To install this package, use devtools:

devtools::install_github("camult/DRP")


# `wideDRP`: Deregressing estimated breeding values - One wide format file

## Description


 This package is easy to use and can be helpful to calculate deregressed proofs,
 and their reliabilities and weights.


## Usage

```r
wideDRP(Data, animalId, sireId, damId, c = 0.5, h2, traitName = NULL,
  animalEBV, sireEBV, damEBV, animalr2, sirer2, damr2)
```


## Arguments

Argument      |Description
------------- |----------------
```Data```     |     It is the name of the data file
```animalId```     |     It is the name of the animal's column
```sireId```     |     It is the name of the sire's column
```damId```     |     It is the name of the dam's column
```c```     |     It is the fraction of genetic variance not explained by markers
```h2```     |     It the heritability of the trait
```traitName```     |     It the name of the trait
```animalEBV```     |     It is the name of the animal's EBV column
```sireEBV```     |     It is the name of the sire's EBV column
```damEBV```     |     It is the name of the dam's EBV column
```animalr2```     |     It is the name of the animal's accuracy column
```sirer2```     |     It is the name of the sire's accuracy column
```damr2```     |     It is the name of the dam's accuracy column

## Value


 A data frame with deregressed proofs, reliability and weights.


## References


 Garrick, D. J., J. F. Taylor, and R. L. Fernando. 2009.
 Deregressing estimated breeding values and weighting information
 for genomic regression analyses. Genet. Sel. Evol. 41:55.


## Examples

```r 
 ## Not to run ##
 
 ## Example from Garrick et al., (2009)
 
 Dataset=data.frame(animal="A1000", sire="S10", dam="D100", ebv_anim=15, ebv_sire=10, ebv_dam=2,
 r2_anim=0.68, r2_sire=0.97, r2_dam=0.36, trait="Trait", c=0.5, h2=0.25)
 
 wideDRP(Data     =  Dataset,
 animalId  = "animal",
 sireId    = "sire",
 damId     = "dam",
 animalEBV = "ebv_anim",
 sireEBV   = "ebv_sire",
 damEBV    = "ebv_dam",
 animalr2  = "r2_anim",
 sirer2    = "r2_sire",
 damr2     = "r2_dam",
 traitName = "trait",
 c         =  0.5,
 h2        =  0.25)
 
 ## End(Not run)
 
 ``` 

# `DRP2files`: Deregressing estimated breeding values - Two long format file

## Description


 This package is easy to use and can be helpful to calculate deregressed proofs,
 and their reliabilities and weights.


## Usage

```r
DRP2files(animalData, parentData, animalCol, sireCol, damCol, parentCol,
  ebvName, r2Name, c = 0.5, h2)
```


## Arguments

Argument      |Description
------------- |----------------
```animalData```     |     It is animal data file
```parentData```     |     It is parents data file
```animalCol```     |     It is the name of the animal's column
```sireCol```     |     It is the name of the animal's dam column
```damCol```     |     It is the name of the animal's dam column
```parentCol```     |     It the name of the parents' column in the parents data file
```ebvName```     |     It is the name of the EBV column
```r2Name```     |     It is the name of the accuracy column
```c```     |     It is the fraction of genetic variance not explained by markers
```h2```     |     It the heritability of the trait

## Value


 A data frame with deregressed proofs, reliability and weights.


## References


 Garrick, D. J., J. F. Taylor, and R. L. Fernando. 2009.
 Deregressing estimated breeding values and weighting information
 for genomic regression analyses. Genet. Sel. Evol. 41:55.


## Examples

```r 
 ## Not to run ##
 
 ## Example from Garrick et al., (2009)
 
 animalData=data.frame(ID="Animal", sire="Sire", dam="Dam", EBV=15, r2=0.68)
 parentData=data.frame(ID=c("Sire", "Dam"), EBV=c(10, 2), r2=c(0.97, 0.36))
 
 DRP2files(animalData=animalData,
 parentData=parentData,
 animalCol = "ID",
 sireCol   = "sire",
 damCol    = "dam",
 parentCol = "ID",
 ebvName   = "EBV",
 r2Name   = "r2",
 c         = 0.5,
 h2        = 0.25)
 
 ## End(Not run)
 
 ``` 

