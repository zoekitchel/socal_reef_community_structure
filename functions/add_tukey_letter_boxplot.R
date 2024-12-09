#add tukey test to plot
generate_boxplot_tukey_labels <- function(dt, var, plot, divide_var_1000 = F, richness = F, non_parametric = F){
  dt[, mainland_island := ifelse(type_wAR %in% c("Artificial_reef", "Natural_mainland"), "Mainland","Island")]
  #run anova
  anova <- aov(data = dt, get(var)~DepthZone*mainland_island)
  
      #Alt non_parametric test
      if(non_parametric == T){
      #Add interaction term to data table
        dt$interaction_depthzone_type <- interaction(dt$mainland_island,dt$DepthZone,sep = ":")
      
       kruskal <- kruskal.test(data = dt, get(var)~interaction_depthzone_type) #p = 0.00000008
      }
  #run tukey
  tukey <- TukeyHSD(anova)
  
  #run non-parametric tukey
    if(non_parametric == T){
      #extract response variable
      var_extract <- dt[[var]]
      
      #Run wilcoxin text (multiple pairwise comparisons; exact = F deals with ties)
      wilcox <- pairwise.wilcox.test(x = var_extract, g = dt$interaction_depthzone_type, p.adjust.method = "bonferroni", exact = F)
      
      #Create full matrix from lower triangle
      p_vals_matrix <-wilcox$p.value
      
      #Convert from triangular to square matrix
      tri.to.squ<-function(x)
      {
        rn<-row.names(x)
        cn<-colnames(x)
        an<-unique(c(cn,rn))
        myval<-x[!is.na(x)]
        mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
        for(ext in 1:length(cn))
        {
          for(int in 1:length(rn))
          {
            if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
            mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
            mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
          }
          
        }
        return(mymat)
      }
      
      mymat <- tri.to.squ(p_vals_matrix)
      
      #Generate comparison letters
      multcompLetters(mymat)
      
      ggboxplot(dt, x = "interaction_depthzone_type", y = var, color = "type_wAR") +theme(axis.text.x = element_text(angle = 90, hjust = 1)) + stat_compare_means()
      View(compare_means(total_abundance_depthzone_site~interaction_depthzone_type, data = dt))

      }
  
  # Extract labels and factor levels from Tukey post-hoc 
  tukey$`DepthZone:mainland_island`[!complete.cases(tukey$`DepthZone:mainland_island`),] <- 0
  Tukey.labels <- multcompLetters4(anova, tukey)
  plot.labels <- Tukey.labels$`DepthZone:mainland_island`
  
  #range of points on plot
  range <- max(dt[[var]], na.rm = TRUE) - min(dt[[var]], na.rm = TRUE)
  
  #DIDN'T END UP USING BUT GOOD TO HAVE
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  #calculate 0.75 and 0.25 quantiles
 # dt[,twentyfive := quantile(get(var),0.25),.(DepthZone, mainland_island)][,seventyfive := quantile(get(var),0.75),.(DepthZone, mainland_island)][,IQR := seventyfive-twentyfive]
 # dt[, label_location := seventyfive + 1.5*IQR]
 # #identify max value less than IQR
 # dt[, IQR_dif := label_location-get(var)]
 # dt[, IQR_dif := ifelse(IQR_dif<0,NA,IQR_dif)]
 # dt[,max_value_lessthan_IQR := ifelse(IQR_dif == min(IQR_dif, na.rm = T),T,F),.(DepthZone,mainland_island)]
 # dt[,true_label := ifelse(max_value_lessthan_IQR == T,get(var)+0.8,NA)]
  dt[,true_label := -(range*0.04)]
  if(divide_var_1000 ==T){
  dt[,true_label := true_label/1000]
  }
  if(richness ==T){
    dt[,true_label := -true_label]
  }
  
  boxplot.dt <- unique(dt[complete.cases(dt$true_label),.(DepthZone, mainland_island, true_label)])
  
  #make plot.labels into data table too
  plot.labels.dt <- data.table(`DepthZone:mainland_island` = names(plot.labels$monospacedLetters), letter = plot.labels$Letters)
  plot.labels.dt[,DepthZone := sub(":.*","",`DepthZone:mainland_island`)][,mainland_island := sub(".*:", "", `DepthZone:mainland_island`)]
  
  #link
  plot.labels.dt <- boxplot.dt[plot.labels.dt, on = c("DepthZone","mainland_island")]
  
  plot.labels.dt <- plot.labels.dt[!(`DepthZone:mainland_island`%in% c("AR_PVR:Island","AR_SM:Island")),]
  
  #set order
  plot.labels.dt[, `DepthZone:mainland_island`:= factor(`DepthZone:mainland_island`, levels = c("Inner:Island","Inner:Mainland", "Middle:Island","Middle:Mainland","Outer:Island", "Outer:Mainland", "Deep:Island","Deep:Mainland","AR_PVR:Mainland","AR_SM:Mainland"))]
  setkey(plot.labels.dt, `DepthZone:mainland_island`)
  
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