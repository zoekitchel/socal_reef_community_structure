#add dunn test to plot
generate_boxplot_nonpar_dunn_letters <- function(dt, var, plot, divide_var_1000 = F, richness = F){
  dt[, mainland_island := factor(ifelse(type_wAR %in% c("Artificial_reef", "Natural_mainland"), "Mainland","Island"))]
  
  #Run Kruskal-Wallis with Groups (not testing for interactions)
  #Create combined variable
  dt[,type_location := factor(paste(DepthZone,mainland_island,sep = "."), levels = 
                                c("Inner.Island","Inner.Mainland","Middle.Island","Middle.Mainland","Outer.Island","Outer.Mainland", "Deep.Island","Deep.Mainland",
                                  "AR_PVR.Mainland","AR_SM.Mainland"))]
  
  kruskal_test <- kruskal.test(get(var)~ type_location, data = dt)
  
  #Identify response var
  response_var <- dt[,get(var)] 
  
  #Pairwise comparisons with Dunn's test
  dunn_test <- dunn.test(response_var, dt$type_location, method = "Holm") #Using Holm method
  
  #Assign letters
  #Extract p-values from Dunn's test result
  p_values <- dunn_test$P.adjusted < 0.05
  
  #Extract names without spaces (multicomp can't deal with spaces)
  Names <- gsub(" ","", dunn_test$comparisons)
  
  #Name p_values
  names(p_values) <- Names
  
  #Perform letter assignment based on p-values
  letters <- multcompView::multcompLetters(p_values)
  
  #range of points on plot
  range <- max(dt[[var]], na.rm = TRUE) - min(dt[[var]], na.rm = TRUE)
  
  #Add positional label
  dt[,true_label := -(range*0.04)]
  
  #Divid var by 1000?
  if(divide_var_1000 ==T){
    dt[,true_label := true_label/1000]
  }
  #Is this for richenss plot? adjust height 
  if(richness ==T){
    dt[,true_label := -true_label]
  }
  
  #unique label positions
  boxplot.dt <- unique(dt[complete.cases(dt$true_label),.(DepthZone, mainland_island, true_label)])
  
  #Make plot.labels into data table too
  plot.labels.dt <- data.table(full_label = names(letters$Letters), letter = letters$Letters)
  plot.labels.dt[, c("DepthZone", "mainland_island") := tstrsplit(full_label, ".", fixed = TRUE)][,full_label := NULL]
  
  #Link back
  plot.labels.dt <- boxplot.dt[plot.labels.dt, on = c("DepthZone","mainland_island")]
  
  #Set order
  plot.labels.dt[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))][,mainland_island := factor(mainland_island, levels = c("Island","Mainland"))]
  setkey(plot.labels.dt, DepthZone, mainland_island)
  
  plot.labels.dt[, x_location := c(0.82,
                                   1.18,
                                   1.82,
                                   2.18,
                                   2.8,
                                   3.18,
                                   3.81,
                                   4.18,
                                   5.01,
                                   6.01
  )]
  
  full_plot <- plot + geom_text(data = plot.labels.dt, aes(x = x_location, y = true_label, label = letter), size = 3)
  
  return(full_plot)
}