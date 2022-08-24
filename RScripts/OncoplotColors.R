get_mutation_colors = function(mutations){
  mutationColors = brewer.pal(n = 8, name = "Dark2")
  mutationColors[9] = "#FF0000"
  mutationColors = setNames(mutationColors, mutations)
  return(mutationColors)
}

get_mutation_alter_fun = function(mutation_colors){
  return(list(
    background = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#DEDEDE", col = NA)),
    
    Missense_Mutation = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["Missense_Mutation"], col = NA)),
    Frame_Shift_Del = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["Frame_Shift_Del"], col = NA)),
    Frame_Shift_Ins = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["Frame_Shift_Ins"], col = NA)),
    Nonsense_Mutation = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["Nonsense_Mutation"], col = NA)),
    Nonstop_Mutation = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["Nonstop_Mutation"], col = NA)),
    Splice_Site = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["Splice_Site"], col = NA)),
    In_Frame_Ins = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["In_Frame_Ins"], col = NA)),
    In_Frame_Del = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["In_Frame_Del"], col = NA)),
    Translation_Start_Site = function(x, y, w, h) 
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["Translation_Start_Site"], col = NA))
    
  ))
}

get_cohort_colors = function() {
  return(c(BRCA = "royalblue4",
              CCRCC = "maroon4",
              COAD = "goldenrod1",
              GBM = "aquamarine",
              HCC = "pink3",
              HNSCC = "coral4",
              LSCC = "slategray3",
              LUAD = "mediumpurple",
              OV = "seagreen",
              PDA = "lightgoldenrod2",
              UCEC = "orangered3"))
  }


get_cohort_alter_fun =function(cohort_colors){
  return(list(
    background = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#DEDEDE", col = NA)),
    
    BRCA = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["BRCA"], col = NA)),
    CCRCC = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["CCRCC"], col = NA)),
    COAD = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["COAD"], col = NA)),
    GBM = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["GBM"], col = NA)),
    HNSCC = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["HNSCC"], col = NA)),
    LSCC = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["LSCC"], col = NA)),
    LUAD = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["LUAD"], col = NA)),
    OV = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["OV"], col = NA)),
    PDA = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["PDA"], col = NA)),
    UCEC = function(x, y, w, h)
      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = cohort_colors["UCEC"], col = NA))
  ))
}
