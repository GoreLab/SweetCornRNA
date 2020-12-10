# calculate BLUP via lme4

# source
library(lme4)

# file I/O
dir.in.rlog <- "RAWDATA/Seetcorn_TagSeq"
file.in.rlog <- "htseq_count_matrix_sweetcorn_B73_RLOG_all_info_v1.txt"
dir.in.key <- "RAWDATA"
file.in.key <- "master_key.csv"
dir.in.expinfo <- "RAWDATA/MetaData"
file.in.expinfo <- "rnaseq_trial_2019_upload.csv"
dir.save <- "RESULT/2.2-BLUP_lmer"

# objects
i <- 5332 # i-th gene to be used
name.tech.ctrl <- c("Control_NA", "fill_c1", "fill_c2", "fill_c3", "fill_c4")
name.check <- c("CHECK1", "CHECK2", "CHECK3", "CHECK4")

# make dir to save result
dir.create(dir.save, recursive = TRUE)

# --------------------------------------------------------------------------- #
# ----- 1. Load data & make a few objects
# --------------------------------------------------------------------------- #
# load rlog data
rlog.df <- read.delim(file = paste0(dir.in.rlog, "/", file.in.rlog)) # a few min

# rlog data matrix
gene.info.df <- rlog.df[, 1:5]
rlog.mat <- as.matrix(rlog.df[, 6:ncol(rlog.df)])
rownames(rlog.mat) <- gene.info.df$gene_id

# add sample name
h <- readLines(con = paste0(dir.in.rlog, "/", file.in.rlog), n = 1) # get header
name.header <- strsplit(h, "\t")[[1]]
colnames(rlog.mat) <- name.header[6:ncol(rlog.df)]

# load master data
key.df <- read.csv(file = paste0(dir.in.key, "/", file.in.key))

# check: list of names are identical? -- yes!
setequal(key.df$rlog_ids, colnames(rlog.mat))

# load experiment info data
expinfo.df <- read.csv(file = paste0(dir.in.expinfo, "/", file.in.expinfo))

# check: info is available for all samples?
setdiff(key.df$sample_name_metadata, expinfo.df$plot_name) # OK -- "Pos" is positive control

# --------------------------------------------------------------------------- #
# ----- 2. make a data frame to fit BLUE/BLUP model
# --------------------------------------------------------------------------- #
# add rlog value of the i-th gene to the key file
df.model <- data.frame("rlog.id" = colnames(rlog.mat),
                       "Expr.rlog" = rlog.mat[i, ],
                       row.names = NULL)
df.model <- merge.data.frame(x = df.model, y = key.df, by.x = "rlog.id", by.y = "rlog_ids")

# remove technical control in RNAseq
df.model <- df.model[!(df.model$accession_name %in% name.tech.ctrl), ]

# merge metadata
length(expinfo.df$plot_name) == length(unique(expinfo.df$plot_name)) # uniqueness: OK
length(df.model$sample_name_metadata) == length(unique(df.model$sample_name_metadata)) # uniqueness: OK
setdiff(df.model$sample_name_metadata, expinfo.df$plot_name) # check inclusion: OK
df.model <- merge.data.frame(x = df.model, y = expinfo.df,
                             by.x = "sample_name_metadata", by.y = "plot_name",
                             all.x = TRUE, all.y = FALSE)

# rename columns for simplicity
colnames(df.model)[colnames(df.model) == "accession_name.y"] <- "accession_name_meta"
colnames(df.model)[colnames(df.model) == "accession_name.x"] <- "accession_name"

# Harvest.Date
days.harvest <- as.integer(as.Date(df.model$tissue_harvest_date) - min(as.Date(df.model$tissue_harvest_date)))

# days.rna.extract should be two-level factor
rna.extract <- rep("spring", nrow(df.model))
rna.extract["2020-05-01" < as.Date(df.model$nucleic_acid_extraction_date)] <- "summer"
rna.extract <- as.factor(rna.extract)

# final data
df.fit <- data.frame("y" = df.model$Expr.rlog,
                     "Switch" = df.model$is_a_control,
                     "Genotype" = df.model$accession_name,
                     "Range" = as.factor(df.model$range_number),
                     "Block" = as.factor(df.model$block_number),
                     "Column" = as.factor(df.model$col_number),
                     "Row" = as.factor(df.model$row_number),
                     "Plate" = as.factor(df.model$plate_number),
                     "Harvest.Date" = days.harvest,
                     "RnaExtract" = as.factor(rna.extract))

# --------------------------------------------------------------------------- #
# ----- 3. fit the model & get results
# --------------------------------------------------------------------------- #
# fit the model (BLUP)
res <- lmer(y ~ dummy(Genotype, name.check) + Harvest.Date
            + (0+dummy(Switch, "0")|Genotype)
            + (1|Range)
            + (1|Range:Block)
            + (1|Range:Block:Column)
            + (1|Row)
            + (1|Plate)
            + (1|RnaExtract), data = df.fit)

# estimated effects & variances
fe <- fixef(res)
re <- ranef(res)
vc <- VarCorr(res)

# variance comp
geno.var <- attr(vc$Genotype, "stddev") ^ 2
range.var <- attr(vc$Range, "stddev") ^ 2
block.var <- attr(vc$`Range:Block`, "stddev") ^ 2
column.var <- attr(vc$`Range:Block:Column`, "stddev") ^ 2
row.var <- attr(vc$Row, "stddev") ^ 2
plate.var <- attr(vc$Plate, "stddev") ^ 2
extract.var <- attr(vc$RnaExtract, "stddev") ^ 2
resid.var <- attr(vc, "sc") ^ 2
vec.var <- c(geno.var, range.var, block.var, column.var, row.var,
             plate.var, extract.var, resid.var) # add plate
names(vec.var) <- c("var.geno", "var.range", "var.block", "var.column", "var.row",
                    "var.plate", "var.extract", "var.resid") # add plate

# grand mean
vec.int <- fe[1]
names(vec.int) <- "Intercept"

# genotypic values of check
fix.ef.check <- fe[2:(length(fe)-1)]
names(fix.ef.check) <- substr(names(fix.ef.check), 28, 999)

# slope for RNA extract date
fix.ef.extract <- tail(fe, 1)

# gentypic values
ran.ef.geno.tmp <- re$Genotype; colnames(ran.ef.geno.tmp) <- NULL
ran.ef.geno <- unlist(ran.ef.geno.tmp)
names(ran.ef.geno) <- rownames(ran.ef.geno.tmp)

# range effects
ran.ef.range.tmp <- re$Range; colnames(ran.ef.range.tmp) <- NULL
ran.ef.range <- unlist(ran.ef.range.tmp)
names(ran.ef.range) <- paste0("range", rownames(ran.ef.range.tmp))

# block effects
ran.ef.block.tmp <- re$`Range:Block`; colnames(ran.ef.block.tmp) <- NULL
ran.ef.block <- unlist(ran.ef.block.tmp)
num.block.vec <- sapply(X = rownames(ran.ef.block.tmp), FUN = function(x){strsplit(x, ":")[[1]][2]} )
names(ran.ef.block) <- paste0("block", num.block.vec)

# column effects
ran.ef.col.tmp <- re$`Range:Block:Column`; colnames(ran.ef.col.tmp) <- NULL
ran.ef.col <- unlist(ran.ef.col.tmp)
num.col.vec <- sapply(X = rownames(ran.ef.col.tmp), FUN = function(x){strsplit(x, ":")[[1]][3]} )
names(ran.ef.col) <- paste0("column", num.col.vec)

# row effects
ran.ef.row.tmp <- re$`Row`; colnames(ran.ef.row.tmp) <- NULL
ran.ef.row <- unlist(ran.ef.row.tmp)
names(ran.ef.row) <- paste0("row", rownames(ran.ef.row.tmp))

# plate effects
ran.ef.plate.tmp <- re$Plate; colnames(ran.ef.plate.tmp) <- NULL
ran.ef.plate <- unlist(ran.ef.plate.tmp)
names(ran.ef.plate) <- paste0("plate", rownames(ran.ef.plate.tmp))

# RnaExtract date effect
ran.ef.extract.tmp <- re$RnaExtract; colnames(ran.ef.extract.tmp) <- NULL
ran.ef.extract <- unlist(ran.ef.extract.tmp)
names(ran.ef.extract) <- rownames(ran.ef.extract.tmp)

# make a long vector for output
vec.out <- c(vec.var, vec.int, fix.ef.check, ran.ef.geno, 
             ran.ef.range, ran.ef.block, ran.ef.col, ran.ef.row,
             ran.ef.plate, ran.ef.extract)

# write the result as a text file
df.save <- data.frame("Model.Term" = names(vec.out), "Est.effect" = vec.out, row.names = NULL)
write.csv(df.save, 
          file = paste0(dir.save, "/BLUP_", formatC(i, digits = 5, flag = "0"), ".csv"),
          row.names = FALSE)
