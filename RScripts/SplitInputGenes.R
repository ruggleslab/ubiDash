splitInputGenes = function(inputGenes){
  inputGenes = toupper(inputGenes)
  inputGenes = gsub(",", " ", inputGenes)
  inputGenes = gsub("  ", " ", inputGenes)
  inputGenes = unlist(strsplit(inputGenes, " "))
  return(inputGenes)
}