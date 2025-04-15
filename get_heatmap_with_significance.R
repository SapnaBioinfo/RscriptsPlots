## Name : Sapna Sharma
## Email : sapna.sharma@helmholtz-munich.de
## This file produces heatmap of matrix estimates and put significance information using *, **, ***
## Step : 1 read in regression result files which has "Protein", "Estimate", "se", "p_value", "p_adjusted"
## Step : 2 Collect all the proteins names which shows significance association
## Step : 3 For signicicant proteins collect Estimates and P-adjust values.
## Step 4 : Make a column with gene names from Protein ID
## Step 5: Make an additional column of indicator
## Step 6: Using ggplot function make heatmap


## Step 1 reading regression model result files

ANTI = read.table("model_d_U3TANTIHY.txt",header=T)
MHYP = read.table("model_d_U3TMHYPOL.txt",header=T)
MIBZ = read.table("model_d_U3TMIBU_R.txt",header=T)
MPOI = read.table("model_d_U3TMOPIO.txt",header=T)

## makeing protein a rowmane (it will help)
row.names(ANTI) = ANTI$Protein
row.names(MHYP) = MHYP$Protein
row.names(MIBZ) = MIBZ$Protein
row.names(MPOI) = MPOI$Protein

## Step 2: Collecting all the sifnificant proteins 
p_ANTI = ANTI$Protein[which(ANTI$p_adjusted < 0.05)]
p_MHYP = MHYP$Protein[which(MHYP$p_adjusted < 0.05)]
p_MIBZ = MIBZ$Protein[which(MIBZ$p_adjusted < 0.05)]
p_MPOI = MPOI$Protein[which(MPOI$p_adjusted < 0.05)]


col_p = c(p_ANTI,p_MHYP,p_MIBZ,p_MPOI)
col_pu = unique(col_p)

### Step 3: Collect model information for significant proteins
pu_ANTI = NULL
pu_MHYP = NULL
pu_MIBZ = NULL
pu_MPOI =  NULL
  
pu_ANTI = ANTI[col_pu,c("Protein","Estimate","p_value","p_adjusted")]
pu_MHYP = MHYP[col_pu,c("Protein","Estimate","p_value","p_adjusted")]
pu_MIBZ = MIBZ[col_pu,c("Protein","Estimate","p_value","p_adjusted")]
pu_MPOI = MPOI[col_pu,c("Protein","Estimate","p_value","p_adjusted")]

pu_ANTI$id = "ANTI"
pu_MHYP$id = "MHYP"
pu_MIBZ$id = "MIBZ"
pu_MPOI$id = "MPOI"

getData = NULL
getData = rbind(pu_ANTI,pu_MHYP,pu_MIBZ,pu_MPOI)

getData = data.frame(getData)

## Step 4 : Make a column with gene names
coll_pn = NULL
for (p in getData$Protein){
  a= strsplit(p,"_")[[1]]
  b = a[length(a)-1] 
  coll_pn = rbind(coll_pn,b)
}

## Add the column
nData = getData[,c("Protein","id","Estimate","p_adjusted")]
nData$Protein = coll_pn


## Step 5: Make an additional column of indicator
nData$Significance <- ifelse(nData$p_adjusted < 0.001, "***",
                            ifelse(nData$p_adjusted < 0.01, "**",
                                   ifelse(nData$p_adjusted < 0.05, "*", "")))

## Step 6: using ggplot geom_tile function make a heatmap

g1 = ggplot(nData, aes(x = Protein, y = id, fill = Estimate)) +
  geom_tile() +
  geom_text(aes(label = Significance), color = "white", size = 4, angle = 90, hjust = 1, vjust = 1) + # Add significance asterisks
  #geom_text(aes(label = sprintf("%.3f", p_adjusted)), color = "white", size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Heatmap of Estimates with p-values", fill = "Estimate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=4))
g1

## Save the figure as a pdf
ggsave(filename = "significant_drug_proteine_association.pdf",g1,dpi=600,width=20,height=10,units="cm")
