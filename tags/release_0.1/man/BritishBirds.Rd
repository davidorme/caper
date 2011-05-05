\name{BritishBirds}
\alias{BritishBirds}
\alias{BritishBirds.tree}
\alias{BritishBirds.data}
\docType{data}
\title{
Conservation status of British birds (Thomas 2008)
}
\description{
The dataset contains a molecular phylogeny of 249 species of British birds and a data frame containing information on the conservation status of 181 of those species. The dataset is taken from Thomas (2008)
}
\usage{data(BritishBirds)}
\details{
The data frame contains 26 variables:

binomial    
common_name 

Red_list    
Amber_list  
Green_list  
Red_amber_list
IUCN

Red_HD
Amber_HD 
Red_amber_HD  

Red_list_BDp
amb_BDp 
Red_amb_BDp 

Red_list_BDr
Red_list_WDp
     
amb_BDr     
amb_WDMp    
amb_spec    
amb_BR      
amb_BL      
amb_WL      
amb_BI      
amb_WI      



Red_amb_BDr 
pop_size    
range_size  

}
\references{
Thomas, G. H. (2008). Phylogenetic distributions of british birds of conservation concern. Proceedings of the Royal Society B-Biological Sciences, 275(1647):2077-2083.
}
\examples{
data(BritishBirds)
BritishBirds <- comparative.data(BritishBirds.tree, BritishBirds.data, binomial)
}
\keyword{datasets}
