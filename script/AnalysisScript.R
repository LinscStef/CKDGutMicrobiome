###This is for RG Kidney function Gut project###

# set data -------------------------------------------------------------
##out
outPath <- "./v3/"
for (i in c("baseChar","BacAbun/Maaslin2","diet","SV","Path")) {if (!dir.exists(file.path(outPath,i))) {dir.create(file.path(outPath,i),recursive = T)}}
##pkg
pkg <- c("dplyr","ggplot2","taxize","pheatmap","RColorBrewer","ggbiplot","vegan","plyr","parallel","doParallel")
for (i in pkg) library(i,character.only = T)

# fun -------------------------------------------------------------
options(stringsAsFactors = F,max.print=10000)
source("./myfun.R")

# set -------------------------------------------------------------
load("D:/Projects/RData/RG2021_Sampleinfo.RData")

testVarRG <- c("Creatinine","Cystatin_C","Impaired_Kidney_Function","eGFR")
eGFR <- "eGFR"

fdr_cut <- 0.1
p_cut <- 0.05
cor_m <- "spearman"

targetBacS <- "Prevotella_copri"


# load-------------------------------------------------------------
pd_RG_bac <- openxlsx::read.xlsx("./pd_RG_bac.xlsx")
Species_abun_RG <- openxlsx::read.xlsx("./Species_abun_RG.xlsx")
metacyc_abun_RG <- openxlsx::read.xlsx("./metacyc_abun_RG.xlsx")
ec_abun_RG <- openxlsx::read.xlsx("./ec_abun_RG.xlsx")
asnA_RG <- openxlsx::read.xlsx("./asnA_RG.xlsx")
RG_SV <- openxlsx::read.xlsx("./RG_SV.xlsx")
SV_tax <- openxlsx::read.xlsx("./SV_tax.xlsx")

# P.copri detected
P.copri.clades.RG <- colnames(Species_abun_RG) %>% str_subset(.,"copri.+[A-Z]+")

# AST-transform ------------------------------------------------------------
RG_Species_ast <- apply(t(Species_abun_RG/100), 2, AST) %>% t %>% as.data.frame()
RG_Path_ast <- apply(t(metacyc_abun_RG/100), 2, AST) %>% t %>% as.data.frame()
RG_EC_ast <- apply(t(ec_abun_RG/100), 2, AST) %>% t %>% as.data.frame()
RG_asnA_ast <- apply(asnA_RG/100, 2, AST) %>% as.data.frame()

###############################analysis###############################


# gut species and kdiney function in RG -------------------------------------------------------------
abun_test_RG <- NULL
for (i in testVarRG) {
  cat("\r", i)
  tmp <- Maaslin2::Maaslin2(input_data=Species_abun_RG/100,
                            input_metadata=pd_RG_bac,
                            output=file.path(outPath,"BacAbun/Maaslin2",paste0(i,"_species")),
                            normalization = "NONE",
                            transform = "AST",
                            fixed_effects = c(i,"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL"),
                            min_abundance = 0,
                            min_prevalence = 0,
                            standardize = T,
                            max_significance = fdr_cut,
                            plot_scatter = F,
                            cores = 14)
  abun_test_RG <- tmp$results %>% subset(.,metadata==i) %>% transform(.,FDR = p.adjust(.$pval,method = "BH")) %>% transform(.,X=feature, Y=i, beta=coef, se=stderr, P=pval, sample_size=N) %>% rbind(abun_test_RG,.)
};my.dev.off()

## format sig species
RG_Bac_sig <- abun_test_RG[abun_test_RG$FDR < fdr_cut,]
sig.kidney.bac.RG <- RG_Bac_sig[c("X","Y","beta")] %>% spread(.,"Y","beta",fill = NA) %>% data.frame(.,row.names = "X")
idx <- colnames(sig.kidney.bac.RG) == eGFR
sig.kidney.bac.RG$Class <- apply(sig.kidney.bac.RG, 1, function(x){if(all(is.na(x))){return(NA)};if(all(na.omit(x[idx]) > 0) & all(na.omit(x[!idx]) < 0)){return("Protect")}else if(all(na.omit(x[idx]) < 0) & all(na.omit(x[!idx]) > 0)){return("Damage")}else{NA}}) %>% factor(.,levels = c("Protect","Damage"))
sig.kidney.bac.RG.rank <- abun_test_RG %>% subset(.,Y==eGFR) %>% data.frame(.,row.names = .$X) %>% .[order(abs(.$beta),decreasing = T),] %>% subset(.,X %in% rownames(sig.kidney.bac.RG)) %>% {sig.kidney.bac.RG[rownames(.),]}
sig.kidney.bac.RG.rank <- sig.kidney.bac.RG.rank[!rownames(sig.kidney.bac.RG.rank) %in% P.copri.clades.RG,] #<<---
eGFR.species <- rownames(sig.kidney.bac.RG.rank)[!is.na(sig.kidney.bac.RG.rank$eGFR)]

## plot by rank Figure 2c
all_p <- RG_Bac_sig %>% .[c("X","Y","beta")] %>% spread(.,"Y","beta",fill = NA) %>% data.frame(.,row.names = "X") %>% .[testVarRG] %>% t
all_p <- all_p[,rownames(sig.kidney.bac.RG.rank)[!is.na(sig.kidney.bac.RG.rank$eGFR)]]
all_p <- all_p[c("Impaired_Kidney_Function","Cystatin_C","Creatinine","eGFR"),]
all_p <- apply(all_p, 1, FUN = function(x){all(is.na(x))}) %>% `!` %>% all_p[.,]
rownames(all_p) <- str_replace_all(rownames(all_p),"_"," ")
colnames(all_p) <- str_replace_all(colnames(all_p),"_"," ")

pdf(sprintf("%s/BacAbun/Figure 2c.sig.RG.eGFR.species.pdf",outPath), width = 10, height = 4) #plot new
idx <- c(-0.03,-0.01, 0, 0.01, 0.03)
ComplexHeatmap::Heatmap(as.matrix(all_p),
                        col = circlize::colorRamp2(idx, c("darkblue","blue", "white", "red","darkred")),
                        cluster_rows = F,
                        cluster_columns = F,
                        na_col = NULL,
                        rect_gp = gpar(col="white",lty=1,lwd=2),
                        
                        show_column_names = T,
                        column_names_gp  = gpar(fontsize = 9,fontface="italic",col="black"),
                        
                        show_row_names = T,
                        row_names_gp = gpar(fontsize = 9,col="black"),
                        
                        heatmap_legend_param = list(title = "Maaslin Coef",at=seq(-0.03,0.03,by=0.015),title_gp = gpar(fontsize = 6, fontface = "bold"),labels_gp  = gpar(fontsize = 6),legend_height=unit(16.5, "mm"),grid_width = unit(3.5, "mm")),
                        
                        width = unit(20, "cm"), height = unit(2.5, "cm")
)
my.dev.off()

# gut pathways and kdiney function in RG -------------------------------------------------------------
path_test_RG <- NULL
for (i in testVarRG) {
  cat("\r", i)
  tmp <- Maaslin2::Maaslin2(input_data=metacyc_abun_RG/100,
                            input_metadata=pd_RG_bac,
                            output=file.path(outPath,"BacAbun/Maaslin2",paste0(i,"_pathways")),
                            normalization = "NONE",
                            transform = "AST",
                            fixed_effects = c(i,"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL"),
                            min_abundance = 0,
                            min_prevalence = 0,
                            standardize = T,
                            max_significance = fdr_cut,
                            plot_scatter = F,
                            cores = 14)
  path_test_RG <- tmp$results %>% subset(.,metadata==i) %>% transform(.,FDR = p.adjust(.$pval,method = "BH")) %>% transform(.,X=feature, Y=i, beta=coef, se=stderr, P=pval, sample_size=N) %>% rbind(path_test_RG,.)
};my.dev.off()
path_test_RG$X <- colnames(metacyc_abun_RG)[match(path_test_RG$X,make.names(colnames(metacyc_abun_RG)))]
Path.eGFR.Sig <- path_test_RG %>% subset(.,FDR < fdr_cut & Y==eGFR) %>% data.frame(.,row.names = .$X)

# gut EC and kdiney function in RG -------------------------------------------------------------
ec_test_RG <- NULL
for (i in testVarRG) {
  cat("\r", i)
  tmp <- Maaslin2::Maaslin2(input_data=ec_abun_RG/100,
                            input_metadata=pd_RG_bac,
                            output=file.path(outPath,"BacAbun/Maaslin2",paste0(i,"_ec")),
                            normalization = "NONE",
                            transform = "AST",
                            fixed_effects = c(i,"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","LDL"),
                            min_abundance = 0,
                            min_prevalence = 0,
                            standardize = T,
                            max_significance = fdr_cut,
                            plot_heatmap = F,
                            plot_scatter = F,
                            cores = 14)
  ec_test_RG <- tmp$results %>% subset(.,metadata==i) %>% transform(.,FDR = p.adjust(.$pval,method = "BH")) %>% transform(.,X=feature, Y=i, beta=coef, se=stderr, P=pval, sample_size=N) %>% rbind(ec_test_RG,.)
};my.dev.off()
ec_test_RG$X <- colnames(ec_abun_RG)[match(ec_test_RG$X,make.names(colnames(ec_abun_RG)))]
for (i in testVarRG) {ec_test_RG[ec_test_RG$Y==i,"FDR"] <- p.adjust(ec_test_RG$P[ec_test_RG$Y==i],method = "BH")}
EC.eGFR.Sig <- ec_test_RG %>% subset(., Y==eGFR & FDR < fdr_cut) %>% data.frame(.,row.names = .$X)
EC.eGFR <- ec_test_RG %>% subset(., Y==eGFR) %>% data.frame(.,row.names = .$X)

# asnA and kdiney function in RG -------------------------------------------------------------
asnA_test_RG <- NULL
for (i in testVarRG) {
  tmp <-Maaslin2::Maaslin2(input_data=asnA_RG/100,
                           input_metadata=pd_RG_bac,
                           output=file.path(outPath,"BacAbun/Maaslin2/tmp"),
                           normalization = "NONE",
                           transform = "AST",
                           fixed_effects = c(i,"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL"),
                           min_abundance = 0,
                           min_prevalence = 0,
                           standardize = T,
                           plot_heatmap = F,
                           plot_scatter = F,
                           max_significance = fdr_cut,
                           cores = 1)
  asnA_test_RG <- tmp$results %>% subset(.,metadata==i) %>% transform(.,FDR = p.adjust(.$pval,method = "BH")) %>% transform(.,X=feature, Y=i, beta=coef, se=stderr, P=pval, sample_size=N) %>% rbind(asnA_test_RG,.)
}

# ammo Path and EC -------------------------------------------------------------
ammoPathEC <- list(
  "PWY-5030: L-histidine degradation III"=c("4.3.1.3","4.2.1.49","3.5.2.7","2.1.2.5","4.3.1.4","3.5.4.9"), 
  "HISDEG-PWY: L-histidine degradation I"=c("4.3.1.3","4.2.1.49","3.5.2.7","3.5.3.8"), 
  "PWY0-845: superpathway of pyridoxal 5'-phosphate biosynthesis and salvage"=c("1.4.3.5"),
  "CITRULBIO-PWY: L-citrulline biosynthesis"=c("3.5.1.2","3.5.1.38"),
  "PWY-6703: preQ0 biosynthesis"=c("4.3.99.3"),
  "PWY-7323: superpathway of GDP-mannose-derived O-antigen building blocks biosynthesis"=c("4.2.1.168"),
  "PWY-5497: purine nucleobases degradation II (anaerobic)"=c("3.5.4.3","3.5.4.8","4.3.1.4","4.3.1.17"),
  "PWY-5022: 4-aminobutanoate degradation V"=c("1.4.1.2"),
  "PWY-6906: chitin derivatives degradation"=c("3.5.99.6")
)
ammoPathEC <- ammoPathEC[names(ammoPathEC) %in% Path.eGFR.Sig$X]
ammoPathEC <- sapply(ammoPathEC,function(x){any(x%in%EC.eGFR.Sig$X)}) %>% ammoPathEC[.]
ammoEC <- unlist(ammoPathEC,use.names = T) %>% unique() %>%  setdiff(.,c("4.2.1.49","3.5.2.7","2.1.2.5","3.5.4.9","3.5.3.8")) %>% .[.%in%EC.eGFR.Sig$X]
ammoPathEC <- lapply(ammoPathEC, function(x){x[x%in%c(ammoEC,c("4.2.1.49","3.5.2.7","2.1.2.5","3.5.4.9","3.5.3.8"))]})
ammoEC <- c(ammoEC,"6.3.1.1") ##target EC path #ammoPathEC ammoEC
ammoPath <- names(ammoPathEC)
## re-names
ammoPathName <- ammoPath %>% sub(".+\\:\\s+",x = .,"") %>% sub("\\s+\\(.+\\)$",x = .,"") %>% sub("-phosphate",x = .,"\r\nphosphate") %>% setNames(.,ammoPath)
ammoECName <- ammoEC %>% paste0("EC ",.) %>% setNames(.,ammoEC) %>% c(.,asnA="asnA")

## ammo Path and EC and eGFR -------------------------------------------------------------
## format
stopifnot(all(Reduce(intersect,list(rownames(Species_abun_RG),rownames(metacyc_abun_RG),rownames(ec_abun_RG),rownames(asnA_RG)))==rownames(pd_RG_bac)))
eGFR.features.abun.RG <- cbind(Species_abun_RG[,eGFR.species],metacyc_abun_RG[,ammoPath],ec_abun_RG[,ammoEC],asnA_RG)
eGFR.features.ast.RG <- cbind(RG_Species_ast[,eGFR.species],RG_Path_ast[,ammoPath],RG_EC_ast[,ammoEC],RG_asnA_ast)
## Figure 3a : eGFR ~ ammo path
targetMeta <- pd_RG_bac %>% transform(.,ID=rownames(.))  %>% .[,c("ID","eGFR")] %>% cbind(.,eGFR.features.ast.RG)
xx <- "PWY0-845: superpathway of pyridoxal 5'-phosphate biosynthesis and salvage" %>%  c(setdiff(ammoPath,.),.) %>% .[. %in% ammoPath]
tmp <- targetMeta %>% gather(.,"Var","Value",-c(eGFR,"ID"))  %>% subset(.,Var %in% c(xx)) %>% transform(.,Path=factor(Var,levels=xx,labels=gsub("(.+\\:\\s*)|(\\(.+\\))","",xx)))
p1 <- ggplot(data = tmp,aes(x=eGFR,y=Value,group=Path,color=Path,fill=Path)) + geom_smooth(method = MASS::rlm,se=T,alpha=.1) + xlab("eGFR, ml/min/1.73 m2") + ylab("AST transformed abundance of\r\npathways related to ammonia production") + GTplus2 + scale_color_manual(values = setNames(pool_col2[1:length(levels(tmp$Path))],levels(tmp$Path))) + scale_fill_manual(values = setNames(pool_col2[1:length(levels(tmp$Path))],levels(tmp$Path))) + guides(color=guide_legend(nrow = 4,byrow = T),fill=guide_legend(nrow = 4,byrow = T)) + theme(plot.margin=unit(rep(.5,4),"cm"),legend.position = "top") + GTs
ggsave(plot = p1,filename = sprintf("%s/Path/Figure 3a.%s and %s.pdf",outPath,"ammo.path",eGFR), width = 5.5, height = 5.3);my.dev.off()
## Figure 3d  : p.c ~ ammo path
tmp <- targetMeta %>% gather(.,"Var","Value",-c(targetBacS,"ID"))  %>% subset(.,Var %in% c(xx)) %>% transform(.,Path=factor(Var,levels=xx,labels=gsub("(.+\\:\\s*)|(\\(.+\\))","",xx)))
p1 <- ggplot(data = tmp,aes(x=get(targetBacS),y=Value,group=Path,color=Path,fill=Path)) + geom_smooth(method = MASS::rlm,se=T,alpha=.1) + xlab(paste("AST transformed abundance of",str_replace(targetBacS,"_"," "))) + ylab("AST transformed abundance of\r\npathways related to ammonia production") + GTplus2 + scale_color_manual(values = setNames(pool_col2[1:length(levels(tmp$Path))],levels(tmp$Path))) + scale_fill_manual(values = setNames(pool_col2[1:length(levels(tmp$Path))],levels(tmp$Path))) + guides(color=guide_legend(nrow = 4,byrow = T),fill=guide_legend(nrow = 4,byrow = T)) + theme(plot.margin=unit(rep(.5,4),"cm"),legend.position = "top") + GTs
ggsave(plot = p1,filename = sprintf("%s/Path/Figure 3d.%s and %s.pdf",outPath,"ammo.path",str_replace_all(string = targetBacS,pattern = "(\\:)|(\\/)",replacement = " ")), width = 5.5, height = 5.3);my.dev.off()
## Figure 3b : p.c ~ ammo ec
for (i in unique(unlist(ammoPathEC,use.names = F))) {
  tmp <- EC.eGFR %>% subset(.,X == i);if(nrow(tmp)!=1) return(next)
  xx <- data.frame(X=pd_RG_bac[,eGFR],Y=RG_EC_ast[,tmp$X])
  p1 <- ggplot(data = xx, aes(x=X,y=Y)) + geom_smooth(method = MASS::rlm,se = T,alpha=.1,color=pool_col0[2],fill=pool_col0[2]) + GTplus3 + xlab("eGFR, ml/min/1.73 m2") + GTplus2 + ylab(paste("EC",i)) + geom_text(data = data.frame(X=quantile(xx$X,0.95,na.rm = T),Y=quantile(xx$Y,0.7,na.rm = T),L=sub("= <","<",sprintf("P = %s",ifelse(tmp$P < 0.001,"< 0.001",ifelse(tmp$P < 0.01,round(tmp$P,3),round(tmp$P,2)))))),inherit.aes = F,aes(x=X,y=Y,label=L,size=10),show.legend = F) + GTs + theme(plot.margin=unit(rep(0.5,4),"cm"))
  ggsave(plot = p1,filename = sprintf("%s/Path/Figure 3b.eGFR and EC.%s.pdf",outPath,str_replace_all(string = i,pattern = "(\\:)|(\\/)",replacement = " ")), width = 3.6, height = 3.3)
}
## Figure 3f : p.c ~ asnA
xx <- with(eGFR.features.ast.RG,cor.test(get(targetBacS),asnA,method=cor_m))
p1 <- ggplot(eGFR.features.ast.RG,aes(x=get(targetBacS),y=asnA)) + geom_point(alpha = .7,size =1.3,color="black",stroke=.01,shape=21,fill="#B6C59E") + geom_smooth(method = "lm",se = T) + GTplus3 + ggtitle(sprintf("Spearman's r = %s",round(xx$estimate,2))) + xlab(paste("AST transformed abundance of",str_replace(targetBacS,"_"," "))) + ylab(paste0("AST transformed abundance of asnA")) +  GTs + theme(plot.margin=unit(rep(.5,4),"cm")) + guides(color=guide_legend(title = element_blank()))  + scale_y_continuous(breaks = seq(0,0.012,by=0.003),labels = seq(0,0.012,by=0.003),limits = c(0,0.012)) + theme(axis.title = element_text(size = 15,color="black"),axis.text = element_text(size=15,color="black"))
ggsave(plot = p1,filename = sprintf("%s/Path/Figure 3f.%s and %s.pdf",outPath,"asnA",str_replace_all(string = targetBacS,pattern = "(\\:)|(\\/)",replacement = " ")), width = 5.5, height = 5.1);my.dev.off()
## Figure 3g : p.c ~ asnA
p1 <- ggplot(targetMeta,aes(x=eGFR,y=asnA)) + geom_point(alpha = .7,size =1.3,color="black",stroke=.01,shape=21,fill="#B6C59E") + geom_smooth(method = "lm",se = T) + GTplus3 + ggtitle(sprintf("P = %s",round(asnA_test_RG[asnA_test_RG$Y==eGFR,"P"],3))) + xlab("eGFR, ml/min/1.73 m2") + ylab(paste0("AST transformed abundance of asnA")) +  GTs + theme(plot.margin=unit(rep(.5,4),"cm")) + guides(color=guide_legend(title = element_blank())) + scale_y_continuous(breaks = seq(0,0.012,by=0.003),labels = seq(0,0.012,by=0.003),limits = c(0,0.012)) + theme(axis.title = element_text(size = 15,color="black"),axis.text = element_text(size=15,color="black"))
ggsave(plot = p1,filename = sprintf("%s/Path/Figure 3g.%s and %s.pdf",outPath,"asnA",eGFR), width = 5.5, height = 5.1);my.dev.off()

# enterotypes in RG -------------------------------------------------------------
testData <- Species_abun_RG[,setdiff(colnames(Species_abun_RG),P.copri.clades.RG)]
data <- t(testData) #test data
data.dist <- dist.JSD(data)
## detected
nclusters <- NULL
for (k in 1:5) { 
  if (k==1) { nclusters[k]=NA } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=clusterSim::index.G1(t(data),data.cluster_temp,  d = data.dist,centrotypes = "medoids")}}
pdf(sprintf("%s/BacAbun/enterotypes.cluster.pdf",outPath), width = 4.5, height = 4.6);plot(nclusters, type="h", xlab="k clusters", ylab="CH index");my.dev.off()
rm(list = c("data","data.dist"))
## cluster DMN base
nclusters <- 2
best <- DirichletMultinomial::dmn(k = nclusters,count=as.matrix(testData),verbose=T)
RG_enterotype <- DirichletMultinomial::mixture(best, assign=T) %>% as.data.frame() %>% setNames(.,"ET") %>% transform(.,ET=paste0("ET",ET))
## format
RG_enterotype$ET <- RG_enterotype$ET %>% str_replace(.,"1","3") %>% str_replace(.,"2","1")  %>% str_replace(.,"3","2")
pd_RG_bac$enterotype <- RG_enterotype$ET[match(rownames(pd_RG_bac),rownames(RG_enterotype))] %>% factor()
## enterotypes contri
RG_enterotype_contri <- ddply(data.frame(enterotype=RG_enterotype$ET,testData,check.rows = F,check.names = F),"enterotype",function(x){colMeans(x[,-1,drop=F])}) %>% data.frame(.,row.names = "enterotype",check.rows = F,check.names = F) %>% t %>% as.data.frame(.,check.rows = F,check.names = F)
## enterotypes - species test
tmp <- Maaslin2::Maaslin2(input_data=testData/100,
                          input_metadata=pd_RG_bac,
                          output=file.path(outPath,"BacAbun/Maaslin2",paste0("enterotype","_species")),
                          normalization = "NONE",
                          transform = "AST",
                          fixed_effects = c("enterotype","Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL"),
                          min_abundance = 0,
                          min_prevalence = 0,
                          standardize = T,
                          max_significance = fdr_cut,
                          plot_scatter = F,
                          cores = 14);my.dev.off()
ET_bac_test <-  tmp$results %>% subset(.,metadata=="enterotype") %>% transform(.,FDR = p.adjust(.$pval,method = "BH")) %>% transform(.,X=feature, Y="enterotype", beta=coef, se=stderr, P=pval, sample_size=N)

# PERMANOVA & kidney function -------------------------------------------------------------
## distance - species
testData <- Species_abun_RG[,setdiff(colnames(Species_abun_RG),P.copri.clades.RG)]/100
RG_species_dist <- as.matrix(vegdist(testData,method = "bray"))
## adonis2 - species
adonis_res.abun <- list()
testDataDist <- RG_species_dist
for (i in testVarRG) {
  cat(i,"\r\n")
  pd <- pd_RG_bac[,c(i,"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL","C_reactive_protein")]
  set.seed(1234);xx <- vegan::adonis2(testDataDist ~ .,data = pd, parallel = 14)
  adonis_res.abun[[i]] <- as.data.frame(xx,check.names=F)[,c(3,5)] %>% setNames(.,c("R2","P")) %>% subset(.,!is.na(P))
}
## plot
for (i in testVarRG) {
  xx <- adonis_res.abun[[i]] %>% .[order(.$R2,decreasing = T),] %>% mutate(.,R2=R2*100,term=factor(rownames(.),levels = rownames(.),labels=str_replace_all(rownames(.),"_"," "))) %>% transform(.,term=str_replace(term,"L$","L-c")) %>% transform(.,term=str_replace(term,"C r","C-r")) %>% mutate(.,term=factor(term,levels = as.character(term)))
  p1 <- ggplot(data = xx,aes(x=term,y=R2)) + geom_bar(stat = "identity",fill = ifelse(xx$P<p_cut,"red","black"),width = .8)  +  xlab("") + ylab("Explained variance proportions (R2,%)") + GTplus + theme(axis.title = element_text(size = 10,color="black"),axis.text.y = element_text(size=10,color="black"),axis.text.x = element_text(size=10,color="black",angle = 35)) + theme(plot.margin=unit(rep(1,4),"cm")) + scale_y_continuous(limits = c(0,0.36),breaks=seq(0,0.3,by=0.1),labels = seq(0,0.3,by=0.1))
  ggsave(plot = p1,filename = sprintf("%s/BacAbun/Figure %s.PERMANOVA.species.%s.pdf",outPath,ifelse(i=="eGFR","2b","S"),i),width = 4.5,height = 3.8) 
}

## distance - pathway
testData <- metacyc_abun_RG/100
RG_pathway_dist <- as.matrix(vegdist(testData,method = "bray"))
## adonis2 - pathway
adonis_res.path <- list()
testDataDist <- RG_pathway_dist
for (i in testVarRG) {
  cat(i,"\r\n")
  pd <- pd_RG_bac[,c(i,"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL","C_reactive_protein")]
  set.seed(1234);xx <- vegan::adonis2(testDataDist ~ .,data = pd, parallel = 14)
  adonis_res.path[[i]] <- as.data.frame(xx,check.names=F)[,c(3,5)] %>% setNames(.,c("R2","P")) %>% subset(.,!is.na(P))
}
## plot
for (i in testVarRG) {
  xx <- adonis_res.path[[i]] %>% .[order(.$R2,decreasing = T),] %>% mutate(.,R2=R2*100,term=factor(rownames(.),levels = rownames(.),labels=str_replace_all(rownames(.),"_"," "))) %>% transform(.,term=str_replace(term,"L$","L-c")) %>% transform(.,term=str_replace(term,"C r","C-r")) %>% mutate(.,term=factor(term,levels = as.character(term)))
  p1 <- ggplot(data = xx,aes(x=term,y=R2)) + geom_bar(stat = "identity",fill = ifelse(xx$P<p_cut,"red","black"),width = .8)  +  xlab("") + ylab("Explained variance proportions (R2,%)") + GTplus + theme(axis.title = element_text(size = 10,color="black"),axis.text.y = element_text(size=10,color="black"),axis.text.x = element_text(size=10,color="black",angle = 35)) + theme(plot.margin=unit(rep(1,4),"cm")) + scale_y_continuous(limits = c(0,0.52),breaks=seq(0,0.5,by=0.1),labels = seq(0,0.5,by=0.1))
  ggsave(plot = p1,filename = sprintf("%s/Path/Figure S.PERMANOVA.pathway.%s.pdf",outPath,i),width = 4.5,height = 3.8) 
}

# PCoA & kidney function -------------------------------------------------------------
## species - Impaired_Kidney_Function
PCoA_abun <- cmdscale(RG_species_dist, k=2, eig = T)
PCoA_abun.point <- data.frame(PCoA_abun$points) %>% setNames(.,paste0("abun_PCo",1:2))
pd_RG_bac <- PCoA_abun.point[match(rownames(pd_RG_bac),rownames(PCoA_abun.point)),] %>% cbind(pd_RG_bac,.)
## plot 
testData <- PCoA_abun 
testDataDist <- RG_species_dist
for (idx in c("enterotype","Impaired_Kidney_Function")) {
  pd <- pd_RG_bac %>% .[,paste0("abun_PCo",1:2)] %>% setNames(.,c("PCo1","PCo2")) %>% cbind(.,pd_RG_bac[,idx,drop=F])
  if (idx=="Impaired_Kidney_Function") {pd[,idx] <- factor(pd[,idx],labels = c("No","Yes"))}
  xx <- vegan::adonis2(testDataDist ~ .,data = pd_RG_bac[,idx,drop=F], parallel = 14)
  p1 <- ggplot(data = pd,aes(PCo1,PCo2, color = get(idx))) + geom_point(size = 1,alpha = 0.8) + ggtitle(sprintf("P = %s",round(as.data.frame(xx)[1,5],3))) + xlab(paste("PCo1 = ",round((testData$eig[1]/sum(testData$eig))*100,digits = 2),"%",sep = "")) + ylab(paste("PCo2 = ",round((testData$eig[2]/sum(testData$eig))*100,digits = 2),"%",sep = "")) + GTplus3 + GTs + theme(plot.margin=unit(rep(.3,4),"cm")) +
    stat_ellipse(aes(group = pd[,idx], fill = pd[,idx], color = pd[,idx]) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F) + scale_color_manual(name=str_replace_all(idx,"_"," "), values = if(idx=="enterotype"){setNames(pool_colET,levels(pd[,idx]))}else{setNames(pool_colBR,levels(pd[,idx]))})
  ggsave(plot = p1,filename = sprintf("%s/BacAbun/Figure %s.PCoA.species.%s.pdf",outPath,ifelse(idx=="Impaired_Kidney_Function","2a","2e"),idx),width = 4.3,height = 4.4) 
}

## pathway - Impaired_Kidney_Function
PCoA_pathway <- cmdscale(RG_pathway_dist, k=2, eig = T)
PCoA_pathway.point <- data.frame(PCoA_pathway$points) %>% setNames(.,paste0("pathway_PCo",1:2))
pd_RG_bac <- PCoA_pathway.point[match(rownames(pd_RG_bac),rownames(PCoA_pathway.point)),] %>% cbind(pd_RG_bac,.)
## plot 
testData <- PCoA_pathway 
testDataDist <- RG_pathway_dist
for (idx in c("enterotype","Impaired_Kidney_Function")) {
  pd <- pd_RG_bac %>% .[,paste0("pathway_PCo",1:2)] %>% setNames(.,c("PCo1","PCo2")) %>% cbind(.,pd_RG_bac[,idx,drop=F])
  if (idx=="Impaired_Kidney_Function") {pd[,idx] <- factor(pd[,idx],labels = c("No","Yes"))}
  xx <- vegan::adonis2(testDataDist ~ .,data = pd_RG_bac[,idx,drop=F], parallel = 14)
  p1 <- ggplot(data = pd,aes(PCo1,PCo2, color = get(idx))) + geom_point(size = 1,alpha = 0.8) + ggtitle(sprintf("P = %s",round(as.data.frame(xx)[1,5],3))) + xlab(paste("PCo1 = ",round((testData$eig[1]/sum(testData$eig))*100,digits = 2),"%",sep = "")) + ylab(paste("PCo2 = ",round((testData$eig[2]/sum(testData$eig))*100,digits = 2),"%",sep = "")) + GTplus3 + GTs + theme(plot.margin=unit(rep(.3,4),"cm")) +
    stat_ellipse(aes(group = pd[,idx], fill = pd[,idx], color = pd[,idx]) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F) + scale_color_manual(name=str_replace_all(idx,"_"," "), values = if(idx=="enterotype"){setNames(pool_colET,levels(pd[,idx]))}else{setNames(pool_colBR,levels(pd[,idx]))})
  ggsave(plot = p1,filename = sprintf("%s/Path/Figure S.PCoA.pathway.%s.pdf",outPath,idx),width = 4.3,height = 4.4) 
}

# alpha & kidney function -------------------------------------------------------------
## species
testData <- Species_abun_RG[,setdiff(colnames(Species_abun_RG),P.copri.clades.RG)]
alpha_species_RG <- calc_alpha(testData) %>% as.data.frame()
## test
alpha_species_RG_p <- data.frame() #p-value adjust
pd <- pd_RG_bac
pd$Impaired_Kidney_Function <- pd$Impaired_Kidney_Function %>% factor(.,levels=c("0","1"),labels = c("No","Yes"))
xcVar <- c("Impaired_Kidney_Function", "eGFR_Grade")
for (i in colnames(alpha_species_RG)) {
  cat(i,"\r\n")
  tmp <- pd[rownames(alpha_species_RG),xcVar]
  pd_f <- colnames(tmp)[sapply(as.list(tmp),is.factor)] %>% tmp[,.,drop=F] %>% cbind(alpha_species_RG[,i,drop=F],.)
  for (j in 2:ncol(pd_f)) {
    #p-value adjust  alpha ~ vars
    tmp <- combn(levels(pd_f[,j]),2)
    for (z in 1:ncol(tmp)) {
      tmp_p <- pd_f[pd_f[,j] %in% tmp[,z],c(1,j)] %>% cbind(.,pd[rownames(.),"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL"]) %>% lm(data = .,sprintf("%s ~ %s + %s",names(.)[1],names(.)[2],paste0("Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL",collapse = "+"))) %>% summary() %>% {.$coefficients[2,4]}
      alpha_species_RG_p <- data.frame(Y=i,X=paste0(tmp[,z],collapse = " vs "),P=tmp_p) %>% rbind(alpha_species_RG_p,.)
    }
    #plot - adjust P
    p1 <- ggplot(data = pd_f[,c(1,j)],aes(x=get(names(pd_f)[j]),y=get(i),fill=get(names(pd_f)[j]))) + geom_boxplot(outlier.size = .1,width=.5,show.legend = F) + ylab(i) + xlab(str_replace_all(str_replace(colnames(pd_f)[j],"eGFR_Grade","eGFR, ml/min/1.73 m2"),"_"," ")) + GTplus2 + theme(plot.margin=unit(rep(.2,4),"cm")) + GTs + scale_fill_manual(values = if(length(levels(pd_f[,j]))==3){setNames(pool_colRG,levels(pd_f[,j]))}else{setNames(pool_colBR,levels(pd_f[,j]))})
    p1 <- p1 + ggsignif::geom_signif( comparisons = as.list(as.data.frame(tmp)),map_signif_level = F,step_increase = 0.12,textsize = 3,tip_length=0.015,size=.5,
                                      annotations = data.frame(subset(alpha_species_RG_p,Y==names(pd_f)[1])) %>% data.frame(.,row.names = .$X) %>% .[apply(tmp, 2, paste0,collapse = " vs "),"P"] %>% {ifelse(. < 0.001, "P < 0.001", ifelse(. < 0.01, sprintf("P = %s",round(.,3)), ifelse(. < 0.05, sprintf("P = %s",round(.,2)), "ns")))})
    ggsave(plot = p1,filename = sprintf("%s/BacAbun/Figure S.alpha_species_%s_%s.pdf",outPath,i,names(pd_f)[j]),width = 3.5,height = 3.3)
  }
}

## pathway
testData <- metacyc_abun_RG
alpha_pathway_RG <- calc_alpha(testData) %>% as.data.frame()
## test
alpha_pathway_RG_p <- data.frame() #p-value adjust
for (i in colnames(alpha_pathway_RG)) {
  cat(i,"\r\n")
  tmp <- pd[rownames(alpha_pathway_RG),xcVar]
  pd_f <- colnames(tmp)[sapply(as.list(tmp),is.factor)] %>% tmp[,.,drop=F] %>% cbind(alpha_pathway_RG[,i,drop=F],.)
  for (j in 2:ncol(pd_f)) {
    #p-value adjust  alpha ~ vars
    tmp <- combn(levels(pd_f[,j]),2)
    for (z in 1:ncol(tmp)) {
      tmp_p <- pd_f[pd_f[,j] %in% tmp[,z],c(1,j)] %>% cbind(.,pd[rownames(.),"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL"]) %>% lm(data = .,sprintf("%s ~ %s + %s",names(.)[1],names(.)[2],paste0("Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL",collapse = "+"))) %>% summary() %>% {.$coefficients[2,4]}
      alpha_pathway_RG_p <- data.frame(Y=i,X=paste0(tmp[,z],collapse = " vs "),P=tmp_p) %>% rbind(alpha_pathway_RG_p,.)
    }
    #plot
    p1 <- ggplot(data = pd_f[,c(1,j)],aes(x=get(names(pd_f)[j]),y=get(i),fill=get(names(pd_f)[j]))) + geom_boxplot(outlier.size = .1,width=.5,show.legend = F) + ylab(i) + xlab(str_replace_all(str_replace(colnames(pd_f)[j],"eGFR_Grade","eGFR, ml/min/1.73 m2"),"_"," ")) + GTplus2 + theme(plot.margin=unit(rep(.2,4),"cm")) + GTs + scale_fill_manual(values = if(length(levels(pd_f[,j]))==3){setNames(pool_colRG,levels(pd_f[,j]))}else{setNames(pool_colBR,levels(pd_f[,j]))})
    p1 <- p1 + ggsignif::geom_signif( comparisons = as.list(as.data.frame(tmp)),map_signif_level = F,step_increase = 0.12,textsize = 3,tip_length=0.015,size=.5,
                                      annotations = data.frame(subset(alpha_pathway_RG_p,Y==names(pd_f)[1])) %>% data.frame(.,row.names = .$X) %>% .[apply(tmp, 2, paste0,collapse = " vs "),"P"] %>% {ifelse(. < 0.001, "P < 0.001", ifelse(. < 0.01, sprintf("P = %s",round(.,3)), ifelse(. < 0.05, sprintf("P = %s",round(.,2)), "ns")))})
    ggsave(plot = p1,filename = sprintf("%s/Path/Figure S.alpha_pathway_%s_%s.pdf",outPath,i,names(pd_f)[j]),width = 3.5,height = 3.3)
  }
}

# SV & kidney function -------------------------------------------------------------
## SV - species distance matrix -------------------------------------------------------------
all_msv_dist <-NULL
for (i in c(1:nrow(SV_tax))){
  cat(i,"\r")
  tmp <- RG_SV[,getSVname(SV_tax$taxid[i],colnames(RG_SV))]
  stopifnot(all(unique(sub(x=colnames(tmp),pattern = "(\\d+)\\..*","\\1"))==SV_tax$taxid[i]))
  tmp <- apply(tmp, 2, function(x){length(unique(na.omit(x)))>1}) %>% tmp[,.]
  tmp <- as.matrix(vegdist(as.data.frame(tmp),method = "canberra",na.rm = T))
  tmp <- dist_rmna(tmp)
  all_msv_dist[[SV_tax$taxid[i]]] <- tmp
}

## SV - species ~ kidney function -------------------------------------------------------------
SV_bac_test <- NULL
for (i in 1:nrow(SV_tax)) {
  cat("\r",i)
  tmp <- SV_tax$spec_abun_name[i]
  if(!is.na(tmp)){testcoVar<- data.frame(pd_RG_bac,Species_abun_RG[,tmp,drop=F])}else{testcoVar<- data.frame(pd_RG_bac)}
  for (j in testVarRG) {
    tmp <- my_adonis_terms(inDist = all_msv_dist[SV_tax$taxid[i]],inDf = pd_RG_bac[,j,drop=F],covDf = testcoVar)
    SV_bac_test <- rbind(SV_bac_test,tmp)
  }
}
SV_bac_test$bac_name <- SV_tax$spec_abun_SV_name[match(SV_bac_test$OutBacSV,SV_tax$taxid)]
for (i in testVarRG) SV_bac_test$FDR[SV_bac_test$testVar==i] <- p.adjust(SV_bac_test$P[SV_bac_test$testVar==i],method = "BH")

## PCoA SV -------------------------------------------------------------
RG_SV_dist <- as.matrix(vegdist(as.data.frame(RG_SV),method = "canberra",na.rm = T)) %>% dist_rmna()
PCoA_SV <- cmdscale(RG_SV_dist, k=2, eig = T)
PCoA_SV.point <- data.frame(PCoA_SV$points) %>% setNames(.,paste0("SV_PCo",1:2))
pd_RG_bac <- PCoA_SV.point[match(rownames(pd_RG_bac),rownames(PCoA_SV.point)),] %>% cbind(pd_RG_bac,.)

## PCoA SV & P.copri & kidney function -------------------------------------------------------------
testDataDist <- SV_tax$taxid[SV_tax$spec_abun_SV_name==targetBacS] %>% all_msv_dist[[.]]
testData <- cmdscale(testDataDist, k=2, eig = T)
for (idx in c("Impaired_Kidney_Function",eGFR)) {
  pd <- pd_RG_bac %>% .[rownames(testDataDist),idx,drop=F] %>% cbind(setNames(as.data.frame(testData$points[rownames(.),]),c("PCo1","PCo2")),.)
  if (idx=="Impaired_Kidney_Function") {pd[,idx] <- factor(pd[,idx],labels = c("No","Yes"))}
  # xx <- vegan::adonis2(testDataDist ~ .,data = pd[,idx,drop=F], parallel = 14) %>%  as.data.frame(xx)[1,5]
  xx <- SV_bac_test[(SV_bac_test$OutBacSV == SV_tax$taxid[SV_tax$spec_abun_SV_name==targetBacS]) & SV_bac_test$testVar==idx,"P"]
  if (is.factor(pd[,idx])) {
    p1 <- ggplot(data = pd,aes(PCo1,PCo2, color = get(idx))) + geom_point(size = 1,alpha = 0.8) + ggtitle(sprintf("P = %s",round(xx,3))) + xlab(paste("PCo1 = ",round((testData$eig[1]/sum(testData$eig))*100,digits = 2),"%",sep = "")) + ylab(paste("PCo2 = ",round((testData$eig[2]/sum(testData$eig))*100,digits = 2),"%",sep = "")) + GTplus3 + GTs + theme(plot.margin=unit(rep(.3,4),"cm")) +
      stat_ellipse(aes(group = pd[,idx], fill = pd[,idx], color = pd[,idx]) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F) + scale_color_manual(name=str_replace_all(idx,"_"," "), values = setNames(pool_colBR,levels(pd[,idx])))
  }else{
    p1 <- ggplot(data = pd,aes(PCo1,PCo2, color = get(idx))) + geom_point(size = 1,alpha = 0.8) + ggtitle(sprintf("P = %s",round(xx,3))) + xlab(paste("PCo1 = ",round((testData$eig[1]/sum(testData$eig))*100,digits = 2),"%",sep = "")) + ylab(paste("PCo2 = ",round((testData$eig[2]/sum(testData$eig))*100,digits = 2),"%",sep = "")) + GTplus3 + GTs + theme(plot.margin=unit(rep(.3,4),"cm")) +
      scale_color_gradientn(name=str_replace_all(idx,"_"," "), colours = (colorRampPalette(c("dodgerblue","blue", "red","indianred1"))(nrow(pd))))
  }
  ggsave(plot = p1,filename = sprintf("%s/SV/Figure S.PCoA.SV.%s.%s.pdf",outPath,targetBacS,idx),width = 4.3,height = 4.4) 
}

# R2 & kidney function -------------------------------------------------------------
lifeStyleFactor <- c("smoking","drinking","light_PA")
pd <- pd_RG_bac
testpd <- list()
testpd$Covariates <- pd[,"Hypertension","Sex","Diabetes","BMI","Cancer","Medication_use","HDL","Age","LDL"]
testpd$Life_style_features <-  pd[,lifeStyleFactor]
testpd$Microbial_features <- str_subset(colnames(pd),"PCo") %>% str_subset(.,paste0(1:2,collapse = "|")) %>% sort %>% pd[,.]  %>% setNames(.,str_replace_all(colnames(.),"abun","species")) # %>% cbind(alpha_div,.)
kidney_R2 <- data.frame()
for (i in testVarRG) {
  # tmp <- lm(data = pd,sprintf("%s ~ %s",i,paste0(unlist(sapply(testpd, colnames),use.names = F),collapse = "+"))) %>% summary()
  for (j in names(testpd)) {
    tmp_pd <- testpd[[j]]
    for (x in colnames(tmp_pd)) {
      tmp <- data.frame(X=tmp_pd[,x],Y=pd[rownames(tmp_pd),i],row.names = rownames(tmp_pd)) %>% na.omit()
      if (is.factor(pd[,i])) {
        x1 <- glm(data = tmp,formula = "Y~X",family = binomial()) %>% summary()
        tmp <- rms::lrm(data = tmp, formula = Y~X)
        kidney_R2 <- data.frame(Type=j,Y=i,X=x,p=x1$coefficients[2,4],R2=tmp$stats["R2"]) %>% rbind.fill(kidney_R2,.)
      }else{
        tmp <- lm(data = tmp,formula = "Y~X") %>% summary()
        kidney_R2 <- data.frame(Type=j,Y=i,X=x,p=tmp$coefficients[2,4],R2=tmp$r.squared) %>% rbind.fill(kidney_R2,.)
      }
    }
  }
}
kidney_R2 <- kidney_R2[order(kidney_R2$p,-kidney_R2$R2,kidney_R2$Type,kidney_R2$Y,decreasing = F),] %>% mutate(.,X=factor(X,levels = (unlist(sapply(testpd, colnames),use.names = F))),Y=factor(Y,levels = c("eGFR","Impaired_Kidney_Function","Creatinine" ,"Cystatin_C"),labels = c("eGFR","Impaired kidney function","Creatinine" ,"Cystatin C")),type=ifelse(p<p_cut,Type,"non.Sig")) %>% transform(.,type=factor(type,levels = c(c(names(testpd),"non.Sig")))) %>% transform(.,R2=R2*100)
levels(kidney_R2$X) <- levels(kidney_R2$X) %>% str_replace_all("_"," ") %>% str_replace_all("DL","DL-c")  %>% Hmisc::capitalize() %>% str_replace_all(" PA"," Physical Activity")
for (i in levels(kidney_R2$Y)) {
  p1 <- ggplot(data = subset(kidney_R2,Y==i),aes(x=X,y=R2,fill=type,color=type)) + geom_bar(stat = "identity",width = .65) + GT + ylab("") + scale_fill_manual(values =setNames(c(pool_col1[1:(length(levels(kidney_R2$type))-1)],"grey"),levels(kidney_R2$type))) + scale_color_manual(values =setNames(c(pool_col1[1:(length(levels(kidney_R2$type))-1)],"grey"),levels(kidney_R2$type))) + GTs + theme(axis.text.x = element_text(angle = 35, hjust = 1,size=7)) + ylab("Explained variance proportions (R2,%)") + xlab("") + ggtitle(str_replace(str_replace(i,"eGFR","eGFR, ml/min/1.73 m2"),"_"," ")) + theme(plot.title = element_text(hjust = .5),axis.title.y = element_text(size=10),axis.text.y =  element_text(size=8))
  p1 <- gg.gap::gg.gap(plot = p1,segments = c(0.05,0.12)*100,ylim = c(0,.16)*100,rel_heights = c(0.05,0,0.011)*100,tick_width = c(0.01,0.02)*100)
  ggsave(plot = p1,filename = sprintf("%s/baseChar/Figure S.R2 for %s.pdf",outPath,i),width = 5.5,height = 4.5)
}



###############################end###############################
gc()


















