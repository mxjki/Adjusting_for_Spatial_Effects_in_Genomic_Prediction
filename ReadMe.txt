#Date: 2019-07-20

Description: this contains the implementations of the proposed method and other existing methods in the empirical study for Data2.

########################################
Data (The wheat data which is public available at https://www.g3journal.org/content/suppl/2013/09/30/g3.113.007807.DC1.)
########################################
GRF_wheat_breeding_foreach_1.R:      		This R function implement the proposed GRF model with respect to the phenotype DH under full irrigated.

GRF_no_spatial_wheat_breeding_foreach_1.R:      This R function implement the proposed GRF model without the consideration of spatial kernel in the 
						covariance structure with respect to the phenotype DH under full irrigated.

########################################
Auxiliary Dataset
########################################
Kin_VanRaden_1.csv:       This is the Kinship matrix by adopting the formula suggested by VanRaden (2008).

PC_10.CSV:  		  This contains the first ten principal components of the SNP genotypes.


########################################
Existing Method
########################################

GAPIT_wheat_breeding.R:      		This R function implement the existing CMLM model with respect to all eight phenotypes.

IB_wheat_breeding_1.R:      		This R function implement the existing IB model with respect to phenotype DH under full irrigated.