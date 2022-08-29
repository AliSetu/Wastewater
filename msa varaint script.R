################# MSA for variants ############
#load libraries
library(gplots)
library(dplyr) #required for modification of data
library(pryr) #required to save basic rplots as object
library(ggplot2)
library(gt)
library(ggpubr)


seq = "Spike Protein"
startnt <- 21563
endnt <- 25384

report_name <- toString(sub(".txt", ".pdf", seq))


#load and prepare data
msa <- read.table("input_data.txt", sep="\t",  header=TRUE, check.names = F)


#prepare a numeric matrix for msa1 data frame
msa1m <- data.matrix(msa) 

#remove the "POS" column
msa1 <- msa[,-1] 

seq = "Spike Protein"
startnt <- 21563
endnt <- 25384

#select from the letter data frame
msa1sel <- msa1[startnt:endnt,]

#select from the numeric matrix
msa1msel <- msa1m[startnt:endnt,]

#for numeric matrix, make sure to add the position to row names
row.names(msa1msel) <- msa1msel[,1] 

#now, remove the POS column from the input data
msa1msel <- msa1msel[,-1]

#select only data where any sample has a variant to create a report
new_df <- msa1sel[!(msa1sel[2] == msa1sel[3] & msa1sel[2] == msa1sel[4] & msa1sel[2] == msa1sel[5] & msa1sel[2] == msa1sel[6] & msa1sel[2] == msa1sel[7] & msa1sel[2] == msa1sel[8] & msa1sel[2] == msa1sel[9]),] 
new_dfT <- t(new_df)

library(reshape2)
melted_mat <- melt(new_dfT)

head(melted_mat)


#draw MSA the plot
plot1 <- ggplot(melted_mat) + #add data and fill cells
  geom_tile(data = melted_mat, aes(x=factor(Var2), y=Var1, fill=value), colour = "black") + #add black border around cells
  geom_text(data = melted_mat, aes(x=factor(Var2), y=Var1,label = value), size=4, family = "Avenir Next") + 
  scale_fill_manual(values = c("-" = "lightgray", "A" = "lightgreen", "C" = "pink", "G" = "lightblue", "T" = "yellow", "V"="red", "*"="white")) +
  #coord_equal() +
  xlab("Nucleotide Position") +
  ylab("Samples") +
  labs(fill = "NT") +
  labs(title= paste("No. of Variants found in", seq, ":" , ncol(new_dfT)  )) +
  theme(axis.text = element_text(size=9, family = "Avenir Next"),
        legend.title = element_text(color = "black", size = 8),
        legend.text = element_text(color = "black"),
        axis.title = element_text(size=6, vjust = 2, face="italic"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )

#make the first plot
plot(plot1)

######################## Variants table #####################

x <- ncol(new_df)
x
y<- nrow(new_df)
y
#create an empty data frame to store the new data
df_new1 <- data.frame(matrix(0, ncol = x+1, nrow = 11))
colnames(df_new1) <- colnames(new_df)

#define i value
i=2

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  df_new1[,i] <- ifelse(new_df[,1]== new_df[,i],"-", paste0(new_df[,1],">",new_df[,i]))
}


#Remove reference genome column and last column (NA column) from the dataframe 
df_new2 <- select(df_new1,-1, -ncol(df_new1))


Variant_pos <- rownames(new_df[1])
#new_pos <- as.data.frame(t(row.names(new_df)))
Variant_pos
#Add position column
final_df <- cbind(Variant_pos, df_new2)
final_df
#Check whether output is correct
head(final_df)

write.table(new_dfT, file= paste(seq, "_variant_info.txt"), row.names = F, sep ="\t")

final_df %>% gt()


#summary table

main.title <- paste("Summary of mutations in", seq)
mutation_summary <- ggtexttable(final_df, rows = NULL, theme = ttheme("light", padding = unit(c(2, 2), "mm"), tbody.style = tbody_style(fill = "white", size = 7 ), colnames.style = colnames_style(fill = "white", size = 7, color = "Black")))
mutation_summary1 <- mutation_summary %>% tab_add_title(text = main.title, face = "bold", padding = unit(0.1, "line"))

mutation_summary1



########### mutations barplot #############

#create an empty data frame to store the new data
df_mut <- data.frame(matrix(0, ncol = x+1, nrow = y))
colnames(df_mut) <- colnames(new_df)


#define i value
i=2

#Create a for loop, to assign variant (i.e. nucleotide mutated to which new variant/nucleotide). Here, we assigning "-" if there is no mutation, else assigning to mutation

for (i in 2:x) {
  
  df_mut[,i] <- ifelse(new_df[,1]== new_df[,i],0,1)
}



#Load library
library(dplyr) #required for modification of data

#Remove reference genome column and last column (NA column) from the dataframe 
df_mut <- select(df_mut, -ncol(df_mut))
#df_mut%>% gt()

#Compute sum of mutations per samples
summ1 <- as.data.frame(colSums(df_mut))

#Samples
samples<- as.data.frame(row.names(summ1))

#Combine samples and sum
summ2 <- cbind(samples, summ1)

#Add column names
colnames(summ2) <- c("Sample", "Mutations")

#make a barplot
bb <- ggplot(summ2, aes(x=Sample, y=Mutations)) + geom_bar(stat="identity", fill="steelblue") +   
  labs(title= paste("No. of mutations per Sample in", seq )) +
  theme(axis.text = element_text(size=10, family = "Avenir Next"),
        legend.title = element_text(color = "black", size = 8),
        legend.text = element_text(color = "black"),
        axis.title = element_text(size=8, vjust = 2, face="italic"),
        axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, size = 8),
        axis.text.y = element_text(color = "black",size = 8, angle = 90,vjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

bb                                                                             

###Grid all samples

grid <- grid.arrange(plot1,arrangeGrob(bb,mutation_summary1,ncol=2), nrow=2)
ggsave(grid, file=report_name, height = 8 , width = 10)

dev.off()