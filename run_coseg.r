library(pedtools)
library(segregatr)

PEN_FILE <- "penetrances/penetrances_Belman_1998_2002.tsv"

COLNAMES <- c("FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age", "Yob", "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn", "BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1", "ER:PR:HER2:CK14:CK56")

DATA <- read.table("example_data/example1.txt", col.names=COLNAMES)


#TODO define cancer entities to consider
BrCa <- T
OvCa <- T
PanCa <- T
ProCa <- T
#TODO Gene/Frequency parameter handling
GENE <- "BRCA1"
FREQ <- 0.001

#TODO check for unique individual names
#TODO check for !=1 probands: proband=x$ID[which(DATA$Target==1)]

### write ped data.frame

id <- DATA["Name"][,1] # ids have to be unique (and will be used for fid/mid columns
#TODO check for doubled ids
id_ind = 1
for (i in seq(length(id)))
	if (id[i] == "0"){
		id[i] <- id_ind
		id_ind <- id_ind+1
	}


fid <- DATA["FathID"][,1]
mid <- DATA["MothID"][,1]
# fill fid/mid with ids
for (i in seq(length(id))){
	if (fid[i] != "0"){
		#cat(i, fid[i], id[as.integer(fid[i])], class(id[as.integer(fid[i])]), "\n")
		fid[i] <- id[as.integer(fid[i])]
	}
	if (mid[i] != "0"){
		#cat(i, mid[i], id[mid[i]], class(id[mid[i]]), "\n")
		mid[i] <- id[as.integer(mid[i])]
	}
}

sex <- DATA["Sex"] # M/F due to CanRisk, has to be converted into M=1, F=2
sex[which(sex == "M"),] <- 1
sex[which(sex == "F"),] <- 2

PED <- data.frame(id, fid, mid, sex)
#TODO define file name & flag for saving file
write.table(PED, 'tmp.ped', append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE)

x <- readPed("tmp.ped") #TODO delete

#carriers: A character vector (or coercible to such), containing the ID
#          labels of pedigree members known to carry the variant in
#          question.

#noncarriers: A character vector (or coercible to such), containing the
#          ID labels of pedigree members known _not_ to carry the
#          variant in question.

#    freq: A single number strictly between 0 and 1: the population
#          frequency of the observed allele.

#affected: The affected pedigree members.

# unknown: Pedigree members with unknown affection status.

# proband: The ID label of the proband. This person must also be in both
#          ‘carriers’ and ‘affected’.

#penetrances: Either a numeric vector of length 3, corresponding to (f0,
#          f1, f2) or a matrix or data frame with 3 columns. Each row
#          contains the penetrance values of a liability class.

#liability: A vector of length ‘pedsize(x)’, containing for each
#          pedigree member the row number of ‘penetrances’ which should
#          be used for that individual. (If ‘penetrances’ is just a
#          vector, it will be used for all classes.) If ‘liability’ is
#          NULL (the default), it is set to ‘1’ for all individuals.


PENETRANCES <- read.table(PEN_FILE, header=T)

CARRIERS <- x$ID[c(which(DATA[GENE] == "T:P"), which(DATA[GENE] == "S:P"))] #DONE
NONCARRIERS <- x$ID[c(which(DATA[GENE] == "T:N"), which(DATA[GENE] == "S:N"))]
AFFECTED <- c()
UNKNOWN <- c() #refers to affection status

# Penetrance matrix will have one entry per individual, hence liability vector is pre-defined
LIAB <- seq(pedsize(x)) #DONE
PENMAT <- matrix(ncol=3, nrow=0)
colnames(PENMAT) <- c("f0", "f1", "f2")

### TEST
#for (i in seq(pedsize(x))){
#	PENMAT <- rbind(PENMAT, c(f0=0, f1=0.5, f2=0.5))
#}
for (i in seq(pedsize(x))){

	if (DATA[i,]$Age == 0){
		UNKNOWN <- c(UNKNOWN, x$ID[i])
		mat <- PENETRANCES[PENETRANCES$GENE == GENE,]
		if (x$SEX[i] == 1){
			# male
			mat <- mat[mat$SEX == "male",]
		}
		else {
			mat <- mat[mat$SEX == "female",]
		}
		mat <- mat[mat$CANCER == "unaff",]
		mat <- mat[(mat$AGE_LOW <= 20) & (20 <= mat$AGE_HIGH),]
		#PENMAT[i,1:3] = mat[1,6:8]
		PENMAT <- rbind(PENMAT, mat[1,6:8])
	}
	else{
		agevec = c(DATA[i,]$Age, DATA[i,]$BC1, DATA[i,]$OC, DATA[i,]$PRO, DATA[i,]$PAN)
		agevec[which(agevec==0)] <- Inf
		
		entity <- c("BrCa", "OvCa", "ProCa", "PanCa")
		age <- c(DATA[i,]$BC1, DATA[i,]$OC, DATA[i,]$PRO, DATA[i,]$PAN)
		flag <- c(BrCa, OvCa, ProCa, PanCa)
		df <- data.frame(entity, age, flag)
		aff_df <- df[(df$age>0) & df$flag,]
		
		# UNAFFECTED
		if (length(aff_df[,1]) == 0){
			cat(x$ID[i], "unaffected", '\n')
		
			#_mat <- PENETRANCES[PENETRANCES$GENE == GENE & PENETRANCES$SEX == "male" & PENETRANCES$CANCER == "unaff",]
			mat <- PENETRANCES[PENETRANCES$GENE == GENE,]
			rownames(mat) <- c()
			if (x$SEX[i] == 1){
				# male
				mat <- mat[mat$SEX == "male",]
				cat('... and male\n') #DEBUG
			}
			else {
				mat <- mat[mat$SEX == "female",]
				cat('... and female\n') #DEBUG
			}
			mat <- mat[mat$CANCER == "unaff",]
			mat <- mat[(mat$AGE_LOW <= agevec[1]) & (agevec[1] <= mat$AGE_HIGH),]
			#PENMAT[i,1:3] = mat[1,6:8]
			PENMAT <- rbind(PENMAT, mat[1,6:8])
		}
		
		# AFFECTED
		else{
			# get min age
			aff_df <- aff_df[which(aff_df$age == min(aff_df$age)),]
			print(aff_df)
			if (length(aff_df[,1]) == 1){
				ent <- aff_df$entity[1]
				cat(x$ID[i], "affected with", ent, '\n')
				mat <- PENETRANCES[PENETRANCES$GENE == GENE,]
				if (x$SEX[i] == 1){
					# male
					mat <- mat[mat$SEX == "male",]
					cat('... and male\n') #DEBUG
				}
				else {
					mat <- mat[mat$SEX == "female",]
					cat('... and female\n') #DEBUG
				}
				mat <- mat[mat$CANCER == ent,]
				mat <- mat[(mat$AGE_LOW <= aff_df$age[1]) & (aff_df$age[1] <= mat$AGE_HIGH),]
				PENMAT <- rbind(PENMAT, mat[1,6:8])
				
				AFFECTED <- c(AFFECTED, x$ID[i])
				
			}
			else{
				#TODO
				cat(x$ID[i], "affected with several cancer entities at same age", '\n')
				print(aff_df)
				stop()
			}
		
		}
		
	}
}


#FLB(x, carriers=CARRIERS, noncarriers=NONCARRIERS, freq=FREQ, affected=AFFECTED, unknown=UNKNOWN, proband=x$ID[which(DATA$Target==1)], penetrances = PENMAT, liability=LIAB)
