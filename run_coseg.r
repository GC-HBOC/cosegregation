library(pedtools)

COLNAMES <- c("FamID", "Name", "Target", "IndivID", "FathID", "MothID", "Sex", "MZtwin", "Dead", "Age", "Yob", "BC1", "BC2", "OC", "PRO", "PAN", "Ashkn", "BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "BARD1", "RAD51D", "RAD51C", "BRIP1", "ER:PR:HER2:CK14:CK56")

DATA <- read.table("example_data/09-00015-005.txt", col.names=COLNAMES)
#TODO check for unique individual names


### write ped data.frame

id <- DATA["Name"][,1] # ids have to be unique (and will be used for fid/mid columns
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

