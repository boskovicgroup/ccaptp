library(ggplot2)
library(reshape)

fda = read.delim('drugs_pred.txt', sep = '\t')
burck10micro = read.delim('Burckhalter-10microMolar.txt', sep ='\t')
burck100nano = read.delim('Burckhalter-0,1microMolar.txt', sep = '\t')
cmld = read.delim('cmld_smiles.smi_out_predictions_20200408-220935.txt', sep = '\t')
gefitinibs = read.delim('gefitinib-series.smi_out_predictions_20200409-024832.txt', sep = '\t')

colnames(gefitinibs)[17:30] <- c("Gefitinib", "Analog 1", "Analog 2", "Analog 3", "Analog 4", 
                              "Analog 5", "Analog 6", "Analog 7", "Analog 8", "Analog 9", 
                              "Analog 10", "Analog 11", "Analog 12", "Analog 13")
plot(hclust(dist(cor(gefitinibs[,17:30], use = 'pairwise.complete.obs'))), 
     main = 'Gefitinib analogs', xlab = 'Distance', ylab = 'Height')


cmld.compounds = colnames(cmld)[17:9720]
fda.compounds = colnames(fda)[17:1405]
burck.compounds = colnames(burck10micro)[17:1891]

hits = list()
nfda.predicted= list()
for (i in seq(0,1,0.05)){
  hits = list()
  for (comp in fda.compounds){
    hits.temp = list()
    hits.temp = as.data.frame(fda[which(fda[comp]>i),'Name'])
      if(nrow(hits.temp)>0){
        hits[[comp]] = hits.temp
      }
      else {
        next
      }
  }
 nfda.predicted = append(nfda.predicted, length(hits))
}

nfda.predicted.frac = lapply(nfda.predicted, function(row_val) row_val/length(fda.compounds))

hits = list()
cmld.predicted= list()
for (i in seq(0,1,0.05)){
  hits = list()
  for (comp in cmld.compounds){
    hits.temp = list()
    hits.temp = as.data.frame(cmld[which(cmld[comp]>i),'Name'])
    if(nrow(hits.temp)>0){
      hits[[comp]] = hits.temp
    }
    else {
      next
    }
  }
  cmld.predicted = append(cmld.predicted, length(hits))
}

cmld.predicted.frac = lapply(cmld.predicted, function(row_val) row_val/length(cmld.compounds))


hits = list()
burck.predicted= list()
for (i in seq(0,1,0.05)){
  hits = list()
  for (comp in burck.compounds){
    hits.temp = list()
    hits.temp = as.data.frame(burck10micro[which(burck10micro[comp]>i),'Name'])
    if(nrow(hits.temp)>0){
      hits[[comp]] = hits.temp
    }
    else {
      next
    }
  }
  burck.predicted = append(burck.predicted, length(hits))
}

burck.predicted.frac = lapply(burck.predicted, function(row_val) row_val/length(burck.compounds))



hits.fda = list()
for (comp in fda.compounds){
  hits.temp = list()
  hits.temp = as.data.frame(fda[which(fda[comp]>0.95),c('Name', 'Protein_Classification')])
  if(nrow(hits.temp)>0){
    hits.fda = rbind(hits.fda, cbind(hits.temp, as.character(comp)))
  }
  else {
    next
  }
}

hits.burck = list()
for (comp in burck.compounds){
  hits.temp = list()
  hits.temp = as.data.frame(burck10micro[which(burck10micro[comp]>0.95),c('Name', 'Protein_Classification')])
  if(nrow(hits.temp)>0){
    hits.burck = rbind(hits.burck, cbind(hits.temp, as.character(comp)))
  }
  else {
    next
  }
}

hits.cmld = list()
for (comp in cmld.compounds){
  hits.temp = list()
  hits.temp = as.data.frame(cmld[which(cmld[comp]>0.95),c('Name', 'Protein_Classification')])
  if(nrow(hits.temp)>0){
    hits.cmld = rbind(hits.cmld, cbind(hits.temp, as.character(comp)))
  }
  else {
    next
  }
}

confidence = seq(0,1,0.05)
df = cbind(confidence, as.numeric(nfda.predicted.frac), as.numeric(cmld.predicted.frac), as.numeric(burck.predicted.frac))
df = as.data.frame(df)
colnames(df) = c('confidence', 'FDA', 'CMLD', 'BHC')
df.melt <- melt(df, id.vars = 'confidence', variable_name = 'Collection')

ggplot(data = df.melt, aes(color = Collection, x = confidence, y = value)) +                    # basic graphical object
  geom_point( size = 2) +
  labs(title = "Target Prediction Confidence", x = 'Prediction confidence', 
       y = 'Fraction of molecules with predicted target') +
  theme(plot.title = element_text(hjust = 0.5))

