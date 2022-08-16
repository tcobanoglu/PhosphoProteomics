inka <- read.delim("inka_urea.xlsx",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# shows you the sample names in the table
library(readxl)
inka_urea <- read_excel("C:/Users/Tugce Su/OneDrive - Vrije Universiteit Amsterdam/Desktop/INKA_code/inka-urea/InKA_20220330/inka_urea.xlsx")

colnames(inka_urea)

m1 <- (match("ATR",inka_urea$Kinase))
m2 <- c("inka_save_Urea_A" ,"inka_save_Urea_B", "inka_save_Urea_E" ,"inka_save_Urea_F",
                  "inka_save_Urea_I", "inka_save_Urea_J")
m3 <- (inka_urea[(m1) ,(m2)])

#pdf(file="CDK2.pdf")
barplot(as.matrix(m3), main="ATR", xlab="Samples", ylab="INKA Score")


