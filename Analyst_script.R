#Cory Williams
 #BF528
 #Project1
 #Role: Analyst
 
 #all scripts were ran in Rstudio on the SCC cluseter
 
 #install limma (microarray anaylysis)
BiocManager::install("limma")
 
#all packages that need to be loaded
library(biobase)
library(limma)
library(tidy)
library(broom)
librarry(affy)  
#import normalized sample CSV file to datafram in R 
 exprs_rma.csv > exprs_rma
 
 #First filter, produce dataframe
 firstfilter <- filter(exprs_rma, GSM971959_JS_05_U133_2.CEL.gz > log2(15))
	> write.csv(firstfilter,file = file.choose(new = T))


#convert first filter datafram to matrix
matrixfirstfilter <- as.matrix(firstfilter)	

#covert matrix to expressionset
firstfilterexprs <- ExpressionSet(matrixfirstfilter)

#obtain pvalue using firstfilter expressionset
pvaliuefirstfilter <- rowttests(, fac, tstatOnly = FALSE) 

#at this point is where the issue started to occur due to the row names (genes) bot appearing after the pvalue was run


#heatmap creation
pheatmap::pheatmap(firstfilter)
	
#test script for alternative way to acheive p value
test%>% 
select(rmadata, GSM971959_JS_05_U133_2.CEL.gz, X1) %>% gather(key = variable, value = value, -X1 %>% group_by(X1, variable) %>% 
summarise(value = list(value)) %>% 
spread(X1, value) %>% 
group_by(variable) %>% 
mutate(p_value = t.test(unlist(1007_s_at), unlist(1053_at))$p.value,
t_value = t.test(unlist(1007_s_at), unlist(1053_at))$statistic)

