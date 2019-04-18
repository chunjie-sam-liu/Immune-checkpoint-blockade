library(coXpress)

# Step 1: reading your data
read.table("/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/pre_PD1_filtered_symbol_expr.txt",header = T,as.is = TRUE) ->all_expr

readxl::read_excel("/data/liull/immune-checkpoint-blockade/all_metadata_available.xlsx",col_names = TRUE,sheet="SRA") %>%
  dplyr::filter(Library_strategy=="RNA-Seq") %>%
  dplyr::filter(Cancer=="melanoma") %>%
  dplyr::filter(Anti_target=="anti-PD1") %>%
  dplyr::filter(Biopsy_Time=="pre-treatment")%>%
  dplyr::select(SRA_Study,Run,Response,Biopsy_Time) ->metadata
dplyr::filter(metadata,Response %in% c("CR","PR","R"))-> response
dplyr::filter(metadata,Response %in% c("SD","PD","NR")) -> non_response
dplyr::select(all_expr,response$Run,non_response$Run)->ordered_all_expr

# Step 2: Cluster data based on a subset of experiments.
hc.gene <- cluster.gene(ordered_all_expr[,1:nrow(response)], "pearson", "average")

# Step 3: cut the tree at a height of 0.4 
g <- cutree(hc.gene, h=0.4)

# Step 4: examine the difference between R and NR
cox <- coXpress(ordered_all_expr, g, 1:nrow(response), (nrow(response)+1):(nrow(response)+nrow(non_response)), times=1000)
write.table(cox,"/data/liull/immune-checkpoint-blockade/New_batch_effect_pipeline/melanoma_PD1/survival/modules/coXpress/all_groups.txt",quote = FALSE,col.names = TRUE)

# Step 5: the results are a data.frame, with one row for each group of genes.
cox[cox$pr.g1 <= .05 & cox$pr.g2 >= 0.05 & cox$N>=7,]## look for groups of genes with >= 7 members that are non-random in group 1 and random in group 2

# Step 6: examine groups of interest graphically     
# look at group 239
plot.compare.group(ordered_all_expr,g,239,1:nrow(response), (nrow(response)+1):(nrow(response)+nrow(non_response)),
                   scale.center=TRUE,scale.scale=TRUE,
                   ylim=c(-5,5))
inspect.group(ordered_all_expr,g,239,1:nrow(response), (nrow(response)+1):(nrow(response)+nrow(non_response)))#group 239's paied members' correlation in the R and NR
