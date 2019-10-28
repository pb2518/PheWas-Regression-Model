######## Function: Apply a Manual Filter to select out specific phenotypes #########

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

#### Step 4: Apply Manual filter for relavent selected phenotypes ####

#provide a list of selected data to be included
f_list<-read.table("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/f_list.txt", header = FALSE,colClasses="character")
e_f<-e[,f_list[,1]] #filter out those select data from the entire 'e' phenotype data table