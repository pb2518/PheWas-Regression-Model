######## Function: Merge Genotype and Phenotype File Tables #########

#check latest versions of R on WinSCP/HPC PuTTy
module avail R

#load R module
module load R/3.4.0
R

#set wd
setwd("/rdsgpfs/general")

#load fread (read in) file associated packages
library(data.table)
library(bit64)

#### Step 1: Reformat data type for each data class ####

#load and select the correct column list of current data types for each column of the UKBB phenotype data
type<-read.table("/rdsgpfs/general/user/pb2518/home/data_type.txt",header=TRUE,sep="\t") #read in data type summary file
type1<-as.character(type[,4]) #select out data type column information

#convert data types into R data classes
type1[type1 == "Sequence"]<-"character"
type1[type1 == "Categorical (single)"]<-"factor"
type1[type1 == "Categorical (multiple)"]<-"factor"
type1[type1 == "Text"]<-"character"
type1[type1 == "Integer"]<-"integer" 
type1[type1 == "Date"]<-"character"
type1[type1 == "Continuous"]<-"numeric"
type1[type1 == "Curve"]<-"character"

#### Step 2: Load phenotype and genotype files ####

#read in phenotype with new type1 R data classes 
#read in genotype file
	#QC: Identify and convert 'fake' NA values into real NA values that can be detected by R using na.strings
d <- fread("/rdsgpfs/general/project/lms-ware-raw/live/UKBB/r/ukb29171.tab", header=TRUE, sep="\t", na.strings=c("Na","NA"), colClasses=type1) #3 mins
c <- read.table("/rdsgpfs/general/user/pb2518/home/genotype_dom_3cols.txt", header = TRUE, sep = " ", na.strings=c("Na","NA"))

#### Step 3: Merge phenotype file to genotype file ####

#merge the two files by f.eid
e <- merge(c, d, by = "f.eid", all = TRUE, sort = FALSE)

	#QC: Check structure format
	dim(e) ##502,536 rows and 4712 columns
	head(e[ ,c(1:3)], n=10)
	head(e[ ,c(1:5)], n=10)
	head(e[ ,c("f.eid", "f.31.0.0", "BAG3", "FHOD3")], n=10)

	#QC: check class at 'c'=original genotype, 'd'=original phenotype, 'e'=after merge
	class(c$BAG3) #was: integer
	class(e$BAG3) #now: integer 
	class(c$FHOD3)
	class(e$FHOD3) 
	class(d$f.31.0.0) #should be catergorical single therefore was: factor
	class(e$f.31.0.0) #now: factor
	class(d$f.36.0.0) #should be text therefore was: character
	class(e$f.36.0.0) #now: character
	class(d$f.48.0.0) #should be numeric therefore was: numeric
	class(e$f.48.0.0) #now: numeric
	class(e$f.53.0.0) #now: character
	
#### Step 4: Output the combined table ####

#output a table for the merged file 'e'
write.table(e, "/rdsgpfs/general/user/pb2518/home/combined_geno_pheno.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Note: file should not be read in as data classes may change, therefore merged variable 'e' was used directly
