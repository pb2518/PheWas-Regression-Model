######### Final Script for PheWAS Regression Model #########

#In R studio

#set wd
setwd("/rdsgpfs/general")

#load fread (read in) file associated packages
library(data.table)
library(bit64)
#library(PheWAS)

#### Step 1: Reformat data type for each data class ####

#load and select the correct column list of current data types for each column of the UKBB phenotype data
type<-as.data.frame(fread("/rdsgpfs/general/user/pb2518/home/data_type.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)) #read in data type summary file
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
c <- read.table("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/genotype/genotype_add.txt", header = TRUE, sep = " ", na.strings=c("Na","NA"))

	#QC: Automatic filter to remove background non-caucasian population by ethnicity
	#get genetic ethnic grouping information from UKBB phenotype file (f.22006.0.0), 1 is Caucasian. NA are non-Caucasians
	d<-d[!is.na(d$f.22006.0.0), ] #select out unique values i.e exclude rows with NA values

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

#### Step 4: Apply Manual filter for relavent selected phenotypes ####

#provide a list of selected data to be included
f_list<-read.table("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/f_list.txt", header = FALSE,colClasses="character")
e_f<-e[,f_list[,1]] #filter out those select data from the entire 'e' phenotype data table

#### Step 5: Regression analysis between BAG3 genotype data and the filtered (104 samples of) phenotype data ####

out_all<-matrix(NA,ncol(e_f),7) #define a matrix comprising 7 columns and rows matching the numer of rows in e_f

i=2 #for specified BAG3 genotype column

#for(i in 2:ncol(e_f[,2:4]))
#{
	for(j in 3:ncol(e_f)) #and for every other filtered phenotype column
	{
		print(paste("column combination: ",i," and ",j,sep="")) #print column combination numbers
		print(colnames(e_f)[j]) #print selected phenotype column name
		t<-class(e_f[,j]) #print selected phenotype column class
		out_all[j,1]<-colnames(e_f)[j] #for this phenotype into the first column of table add the phenotype column name
		out_all[j,2]<-t #for this phenotype into the second column of table add the phenotype class
		out_all[j,3]<-length(unique(e_f[,j])) #for this phenotype into the third column of table add the number of unique values for this phenotype variable

		fmla <- as.formula(paste0(colnames(e_f)[j]," ~ ", colnames(e_f)[i],"+f.31.0.0")) #formula vector with sex as a covariate
		#print(fmla)
		
		tryCatch({ #function allows handling of unusual errors and warnings
		
		if (t=="integer" | t=="numeric") #if the class of this phenotype is integer or numeric
		{
		
		#print(paste0("data type is: ",t))
		#print(length(unique(e_f[,j])))
		out_all[j,4]<-"linear" #for this phenotype into the fourth column of table add the 'linear' relating to the regression type
		out<-lm(formula = fmla, data=e_f,family=binomial(), na.action=na.exclude) #linear regression #QC: remove rows with NA values
		out1<-summary(out) 
		#print(out1) #print output of regression summary statistics
		
		} else { #otherwise if class of this phenotype is any other class
	
		out_all[j,4]<-"logistic" #for this phenotype into the fourth column of table add the 'logistic' relating to the regression type
		out<-glm(formula = fmla, data=e_f,family=binomial(), na.action=na.exclude) #logistic regression #QC: remove rows with NA values
		out1<-summary(out) 
		#print(out1) #print output of regression summary statistics
		
		}
		
		#p<-min(out1$coefficients[-1,4])
		#beta<-out1$coefficients[which(out1$coefficients[,4]==p),1]
		#name<-rownames(out1$coefficients)[which(out1$coefficients[,4]==p)]
		
		p<-out1$coefficients["BAG3",4] #select the p value from the summary statistics
		beta<-out1$coefficients["BAG3",1] #select the beta value from the summary statistics
		name<-rownames(out1$coefficients)[which(out1$coefficients[,4]==p)] #select the name of the type of regression according to the descriptive statistics
		
		out_all[j,5]<-name #for this phenotype and genotype add name to the fifth column 
		out_all[j,6]<-beta #for this regression add the beta value into the sixth column
		out_all[j,7]<-p #for this regression add the p value into the seventh column
		
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #in case of error, skip error and continue loop
	}
#}

out_all<-out_all[-c(1,2),] #display output minus the first 2 

type[,2] <- paste0("f.",gsub('-', '.', type[,2]))
out_all[,8]<-type[match(out_all[,1], type[,2], nomatch = NA_integer_, incomparables = NULL),5] #add phenotype description to column 8

colnames(out_all)<-c("phenotype","type","unique_value","regression","genotype","beta","p","description")
write.table(out_all,"/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/out_phewas_FHOD3.txt",col.names=FALSE, row.names=FALSE,quote=FALSE,sep="\t")

#out_plot<-as.data.frame(out_all[,c(1,7)])
#out_plot[,1]<-as.character(out_plot[,1])
#out_plot[,2]<-log10(as.numeric(as.character(out_plot[,2]))

#colnames(out_plot)<-c("phenotype","value")

#png("/rdsgpfs/general/project/lms-ware-analysis/live/xiao/M.sc/phewas_result_plot.png", res=150, width = 1000, height = 1000)
#phenotypePlot(out_plot,max.y=9,suggestive.line=-log10(1e-05), significant.line=-log10(5e-08), title="BAG3",use.color=F,x.group.labels=F)
#dev.off()