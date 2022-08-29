library(gplots)
library(dplyr) #required for modification of data
library(pryr) #required to Betave basic rplots as object
library(ggplot2)
library(gt)
library(ggpubr)
library(gridExtra)
library(reshape2)



filenames <- list.files(pattern="group_")
filenames
#Strain Name
AlphaVOC <- "Alpha: B.1.1.7"
BetaVOC <- "beta: B.1.351"
GammaVOC <- "Gamma: P.1"
DeltaVOC <-  "Delta: B.1.617.2"
OMVOC <- "Omicron:B.1.1.529"

Report_Type <- c("Alpha", "gamma")
#Report_Type <- c("Beta", "Gamma")
#Report_Type <- c("Delta", "OM")

VOC <- c(AlphaVOC, GammaVOC)
#VOC <- c(BetaVOC, GammaVOC)
#VOC <- c(BZLVOC, DeltaVOC)
#VOC <- c(DeltaVOC, OMVOC)

#Provide variant positions
Alphapos <- factor(c(8782,14408,14676,15279,23063,23403,28881,28882))
Betapos <- factor(c(1059,8782,14408,23403,25563,23063,23012))
Gammapos <- factor(c(733, 2749, 3828, 5648, 12778,13860))
Deltapos <- factor(c( 4181, 5184, 9891, 11201,11332,  11418, 11514,19220, 27874, 22227, 22917, 22995, 27638,28881, 29402))
OMpos <- factor(c(15240, 23202, 23525, 23599, 24424, 27259))



#OMpos <- factor(c(241, 2832, 3037, 5386, 5924, 6513, 6514, 6515, 8393, 10029, 10449, 11285, 11286, 11287, 11288, 11289, 11290, 11291, 
#                11292, 11293, 11537, 13195, 14408, 15240, 18163, 21762, 21765, 21766, 21767, 21768, 21769, 21770, 21846, 21987, 21988, 
#                21989, 21990, 21991, 21992, 21993, 21994, 21995, 22194, 22195, 22196, 22578, 22673, 22674, 22679, 22813, 22882, 22898, 
#                23202, 23403, 23525, 23599, 23604, 23854, 23948, 24130, 24424, 24469, 24503, 25000, 25584, 26270, 26530, 26577, 26709, 
#                27259, 27527, 27807, 28271, 28311, 28362, 28363, 28364, 28365, 28366, 28367, 28368, 28369, 28370, 28881, 28882, 28883))



make_report <- function(sample_name1,VOC){

  #sample_name1 <- filenames[1]
  if(VOC== AlphaVOC){
    my_positions <- Alphapos
    Num_VOC <- as.numeric(length(Alphapos))
    report_type <- c("Alpha")
  } else if(VOC== BetaVOC){
    my_positions <- Betapos
    Num_VOC <- as.numeric(length(Betapos))
    report_type <- c("Beta")
  } else if(VOC== GammaVOC){
    my_positions <- Gammapos
    Num_VOC <- as.numeric(length(Gammapos))
    report_type <- c("Gamma")
  } else if(VOC== DeltaVOC){
    my_positions <- Deltapos
    Num_VOC <- as.numeric(length(Deltapos))
    report_type <- c("Delta")
  } else if(VOC== OMVOC){
    my_positions <- OMpos
    Num_VOC <- as.numeric(length(OMpos))
    report_type <- c("Omicron")
  }
  
  report_type1 <- as.character(paste0("_",report_type,".pdf"))
  report_name <- toString(sub(".txt", report_type1, sample_name1))
  
  print(report_type)
  print(my_positions)
  print(Betapos)
  print(report_name)
  
  #Read data having variant coverage information
  mydata1 <- read.table(sample_name1, sep="\t", header=TRUE)
  #mydata1 <- read.table("sample_1_18-11-17-2021_S8_L001_1__pair_combined.pnBase.txt", sep="\t", header=TRUE)
  mydata1
  #remove insertion type column
  mydata1 <- mydata1[,-10]
  #Compute coverage
  Coverage <- rowSums(mydata1[,4:9])
  
  #Add coverage column to the data
  mydata1 <- cbind(mydata1,Coverage)
  
  #Extract selected position columns
  mydata2 <- mydata1[mydata1$Pos %in% my_positions,-1]
  
  #Remove extra column and make final data with vital info
  cleanCoverage <- mydata2[,-9]
  cleanCoverage
  
  #Melt data
  mydataMelt <- melt(cleanCoverage, id=c("Pos","Nucl"))
  mydataMelt
  
  #Remove coverage rows
  mydataMelt1 <- subset(mydataMelt, !(variable == "Coverage"))
  
  # FDeltavariant information
  TopVariants <- mydataMelt1 %>% 
    group_by(Pos) %>% 
    mutate(max_score = max(value)) %>% 
    ungroup() %>% 
    filter(value == max_score)
  
  #Print top variants
  TopVariants
  
  #Extract coverage
  mydataMelt_cov <- mydataMelt %>% filter(variable == "Coverage")
  
  #Complete data with coverage
  cov <- mydataMelt_cov %>% 
    group_by(Pos) %>% 
    mutate(max_score = max(value)) %>% 
    ungroup() %>% 
    filter(value == max_score)
  
  #Print mat
  cov
  
  #Create Betame order based on pos for variant and coverage 
  TopVariants1 <- TopVariants[order(TopVariants$Pos),]
  cov1 <- cov[order(cov$Pos),]
  
  #Make final table with key information of variant and coverage
  TopVariants_new <- cbind(TopVariants1[1:4], cov1$max_score)
  TopVariants_new 
  
  #Add Column names
  colnames(TopVariants_new ) <- c("Pos", "Nucl", "Variant", "Var_reads", "Coverage")
  
  #Compute %age for coverage
  Fraction <- round((TopVariants_new$Var_reads/ TopVariants_new$Coverage) * 100, 1)
  #Fraction <- paste(Fraction,"%")
  
  #Combine variant data with %age column
  TopVariants_new1 <-cbind(TopVariants_new, Fraction)
  
  #Compute whether it is a variant
  TopVariants_new1$VarYN <- ifelse(as.character(TopVariants_new1$Variant) == as.character(TopVariants_new1$Nucl),"N","Y")


  #compute proportion of mutations
  Var_num<- as.data.frame(subset(TopVariants_new1, VarYN=="Y"))
  #Var_num
  num <- as.numeric(count(Var_num))
  
 
  #pro <- round(num/Num_VOC, 2)
  #Coverage level for variant positions
  high_pos <- as.data.frame(subset(Var_num, Fraction > 20 ))
  high <- as.numeric(count(as.data.frame(subset(Var_num, Fraction > 20 ))))
  high 
  if (high == 0){ high_pos1 <- "None"   }
 
 else { high_pos1 <- high_pos$Pos }
  
  low_pos <- as.data.frame(subset(Var_num, Fraction <= 20 ))
  low <- as.numeric(count(as.data.frame(subset(Var_num, Fraction <= 20))))
  low
  if (low == 0){  low_pos1 <-  "None"  }
  else {  low_pos1 <- low_pos$Pos   }

  
  
  #Add Column names
  colnames(TopVariants_new1) <- c("Pos", "Nucl", "Variant", "Var_reads", "Coverage", "Fraction(%)", "VarYN")
  
  
  #make a table of top variants for report
  TopVariants_new1 %>% gt()
  
  #make a plotting theme
  mytheme <- theme(
    text=element_text(family="Helvetica", size=10),
    axis.text.x = element_text(angle = 90,size = 5),
    axis.text.y = element_text(size = 5),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, color = "grey75"),
    axis.ticks = element_line(color = "grey85"),
    panel.grid.major = element_line(color = "grey95", size = 0.2),
    panel.grid.minor = element_line(color = "grey95", size = 0.2),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
  
  #summary table
  variant_summary <- ggtexttable(TopVariants_new1, theme=ttheme("light"))
  variant_summary <- variant_summary %>%
    tab_add_title(text = paste("Variant Summary of ", sample_name1, "for", VOC, "\n", "Proportion of Variant mutations in sample:", num, "/", Num_VOC,
                               "\n", "Number of Variant Position with high Coverage (>20%):", high,
                               "\n", "Number of Variant Position with Low Coverage ( <20%):", low), 
                               face = "plain", size = 10)
  
  #plot coverage
  p1 <-ggplot(data=mydata1, aes(x=Pos, y=Coverage, aspect_ratio = 2.5)) +
    geom_area(color="lightgray",fill="lightgray",alpha=0.5) +
    geom_point(aes(color = Coverage)) +
    scale_color_gradientn(colours=c('red','orange',"cyan","blue",'blueviolet')) +
    scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 50)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggtitle(paste("Coverage for ",sample_name1)) +
    mytheme
  
  #plot coverage per variant
  mydata3 <- mydata2[,-10]
  mydata4 <- cbind(mydata3,Var = TopVariants1$variable)
  mydata4 <- melt(mydata4, id=c("Pos","Nucl", "Var"))
  
  p2 <-ggplot(mydata4, aes(x=as.factor(Pos)), y=value) +
    geom_col(aes(y=value, fill=variable, group=Nucl), width = 0.7) +
    ggtitle(paste("Variant details for ",VOC)) +
    geom_text(aes(y = -100, label = Nucl), color = "black") +
    geom_text(aes(y = max(value)+200, label = Var), color = "red") +
    mytheme
  
 #grid <- grid.arrange(p1,arrangeGrob(p2,variant_summary,ncol=2), nrow=3)
 grid <- ggarrange(p1, p2, variant_summary, heights = c(10,10,10), ncol=1)
 ggsave(grid, file=report_name, height = 30, width = 10)
  
  
  #grid <- ggarrange(variant_summary, heights = 10, ncol=1)
  #ggsave(grid, file=report_name, height = 10, width = 10)
  
  
  #sample_name2 <- as.character(paste0(sample_name1,VOC))
  #grid <- grid.arrange(variant_summary, nrow=1)
  #ggsave(grid, file=paste0(sample_name2,".pdf"), height = 15 , width = 12, device = pdf,  dpi = 300) ## size of pdf file can be changed depending upon the no. of variants
  
  
  #grid <- ggarrange(variant_summary, heights = c(10), ncol=1)
  #ggBetave(grid, file=report_name, height = 10, width = 10)
  
  
   
  #grid <- ggarrange(p1, p2, variant_summary, heights = c(10,10,35), ncol=1)
  #ggBetave(grid, file=report_name, height = 45, width = 12)
  
}

for (i in 1:length(filenames)){
  make_report(sample_name1=filenames[i],VOC=AlphaVOC)
  make_report(sample_name1=filenames[i],VOC=BetaVOC)
  make_report(sample_name1=filenames[i],VOC=GammaVOC)
  make_report(sample_name1=filenames[i],VOC=DeltaVOC)
  make_report(sample_name1=filenames[i],VOC=OMVOC)
}

dev.off()

