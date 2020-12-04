cd /workdir/rt475/TWAS_2020-07
mkdir LOGFILE

# Run 1.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.1-MakeDataset.R --args v1.1 B73 > LOGFILE/log_1.1_v1.1_B73.txt &
nohup R --vanilla --slave < 1.1-MakeDataset.R --args v1.1 PH207 > LOGFILE/log_1.1_v1.1_PH207.txt &
nohup R --vanilla --slave < 1.1-MakeDataset.R --args v1 B73 > LOGFILE/log_1.1_v1_B73.txt &
nohup R --vanilla --slave < 1.1-MakeDataset.R --args v1 PH207 > LOGFILE/log_1.1_v1_PH207.txt &

# Run 1.2
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.2-AttachGDD.R --args v1.1 B73 > LOGFILE/log_1.2_v1.1_B73.txt &
nohup R --vanilla --slave < 1.2-AttachGDD.R --args v1.1 PH207 > LOGFILE/log_1.2_v1.1_PH207.txt &
nohup R --vanilla --slave < 1.2-AttachGDD.R --args v1 B73 > LOGFILE/log_1.2_v1_B73.txt &
nohup R --vanilla --slave < 1.2-AttachGDD.R --args v1 PH207 > LOGFILE/log_1.2_v1_PH207.txt &

# Run 1.3.1 -- this part is done by other shell scripts
cd /workdir/rt475/TWAS_2020-07
nohup sh 1.3-RunMultipleR_v1.1_B73.sh > LOGFILE/log_1.3_v1.1_B73.txt &
nohup sh 1.3-RunMultipleR_v1_B73.sh > LOGFILE/log_1.3_v1_B73.txt &
nohup sh 1.3-RunMultipleR_v1.1_PH207.sh > LOGFILE/log_1.3_v1.1_PH207.txt &
nohup sh 1.3-RunMultipleR_v1_PH207.sh > LOGFILE/log_1.3_v1_PH207.txt &

# Run 1.3.2
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.3.2-BLUE_merge.R --args v1.1 B73 > LOGFILE/log_1.3.2_v1.1_B73.txt &
nohup R --vanilla --slave < 1.3.2-BLUE_merge.R --args v1 B73 > LOGFILE/log_1.3.2_v1_B73.txt &
nohup R --vanilla --slave < 1.3.2-BLUE_merge.R --args v1.1 PH207 > LOGFILE/log_1.3.2_v1.1_PH207.txt &
nohup R --vanilla --slave < 1.3.2-BLUE_merge.R --args v1 PH207 > LOGFILE/log_1.3.2_v1_PH207.txt &

# Run 1.3.3
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.3.3-BLUE_RmErr.R --args v1.1 B73 > LOGFILE/log_1.3.3_v1.1_B73.txt &
nohup R --vanilla --slave < 1.3.3-BLUE_RmErr.R --args v1 B73 > LOGFILE/log_1.3.3_v1_B73.txt &
nohup R --vanilla --slave < 1.3.3-BLUE_RmErr.R --args v1.1 PH207 > LOGFILE/log_1.3.3_v1.1_PH207.txt &
nohup R --vanilla --slave < 1.3.3-BLUE_RmErr.R --args v1 PH207 > LOGFILE/log_1.3.3_v1_PH207.txt &

# Run 1.3.4
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 1.3.4-BLUE_SampleFiltering.R --args v1.1 B73 > LOGFILE/log_1.3.4_v1.1_B73.txt &
nohup R --vanilla --slave < 1.3.4-BLUE_SampleFiltering.R --args v1 B73 > LOGFILE/log_1.3.4_v1_B73.txt &
nohup R --vanilla --slave < 1.3.4-BLUE_SampleFiltering.R --args v1.1 PH207 > LOGFILE/log_1.3.4_v1.1_PH207.txt &
nohup R --vanilla --slave < 1.3.4-BLUE_SampleFiltering.R --args v1 PH207 > LOGFILE/log_1.3.4_v1_PH207.txt &

# Run 2.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 2.1-Peer_Use25Fact.R --args v1.1 B73 > LOGFILE/log_2.1_v1.1_B73.txt &
nohup R --vanilla --slave < 2.1-Peer_Use25Fact.R --args v1 B73 > LOGFILE/log_2.1_v1_B73.txt &
nohup R --vanilla --slave < 2.1-Peer_Use25Fact.R --args v1.1 PH207 > LOGFILE/log_2.1_v1.1_PH207.txt &
nohup R --vanilla --slave < 2.1-Peer_Use25Fact.R --args v1 PH207 > LOGFILE/log_2.1_v1_PH207.txt &

# Run 2.2
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 2.2-MakeScreePlot.R --args v1.1 B73 > LOGFILE/log_2.2_v1.1_B73.txt &
nohup R --vanilla --slave < 2.2-MakeScreePlot.R --args v1 B73 > LOGFILE/log_2.2_v1_B73.txt &
nohup R --vanilla --slave < 2.2-MakeScreePlot.R --args v1.1 PH207 > LOGFILE/log_2.2_v1.1_PH207.txt &
nohup R --vanilla --slave < 2.2-MakeScreePlot.R --args v1 PH207 > LOGFILE/log_2.2_v1_PH207.txt &

# Run 2.3
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 2.3-Peer_UseOptFact.R --args v1.1 B73 11 > LOGFILE/log_2.3_v1.1_B73.txt &
nohup R --vanilla --slave < 2.3-Peer_UseOptFact.R --args v1 B73 5 > LOGFILE/log_2.3_v1_B73.txt &
nohup R --vanilla --slave < 2.3-Peer_UseOptFact.R --args v1.1 PH207 12 > LOGFILE/log_2.3_v1.1_PH207.txt &
nohup R --vanilla --slave < 2.3-Peer_UseOptFact.R --args v1 PH207 6 > LOGFILE/log_2.3_v1_PH207.txt &

# Run 3.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 3.1-OutlierRemoval.R --args v1.1 B73 11 > LOGFILE/log_3.1_v1.1_B73.txt &
nohup R --vanilla --slave < 3.1-OutlierRemoval.R --args v1 B73 5 > LOGFILE/log_3.1_v1_B73.txt &
nohup R --vanilla --slave < 3.1-OutlierRemoval.R --args v1.1 PH207 12 > LOGFILE/log_3.1_v1.1_PH207.txt &
nohup R --vanilla --slave < 3.1-OutlierRemoval.R --args v1 PH207 6 > LOGFILE/log_3.1_v1_PH207.txt &

# Run 3.2 (summary of outlier detection)
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 3.2-OutlierRemoval_MakeSummary.R --args v1.1 B73 > LOGFILE/log_3.2_v1.1_B73.txt &
nohup R --vanilla --slave < 3.2-OutlierRemoval_MakeSummary.R --args v1 B73 > LOGFILE/log_3.2_v1_B73.txt &
nohup R --vanilla --slave < 3.2-OutlierRemoval_MakeSummary.R --args v1.1 PH207 > LOGFILE/log_3.2_v1.1_PH207.txt &
nohup R --vanilla --slave < 3.2-OutlierRemoval_MakeSummary.R --args v1 PH207 > LOGFILE/log_3.2_v1_PH207.txt &

# Run 3.3
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 3.3-OutlierRemoval_RemoveExtreme.R --args v1 B73 > LOGFILE/log_3.3_v1_B73.txt &
nohup R --vanilla --slave < 3.3-OutlierRemoval_RemoveExtreme.R --args v1 PH207 > LOGFILE/log_3.3_v1_PH207.txt &

# Merge overlapped genes: do this before FC test
cd /workdir/rt475/TWAS_2020-07/RAWDATA/Annotation
nohup R --vanilla --slave < 1-MergeGenes_B73.R > log_1.txt &
nohup R --vanilla --slave < 2-MergeGenes_PH207.R > log_2.txt &
nohup R --vanilla --slave < 3-CompareGff_PH207.R > log_3.txt &

# Run 4.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 4.1-TWAS.R --args v1.1 B73 > LOGFILE/log_4.1_v1.1_B73.txt &
nohup R --vanilla --slave < 4.1-TWAS.R --args v1 B73 > LOGFILE/log_4.1_v1_B73.txt &
nohup R --vanilla --slave < 4.1-TWAS.R --args v1.1 PH207 > LOGFILE/log_4.1_v1.1_PH207.txt &
nohup R --vanilla --slave < 4.1-TWAS.R --args v1 PH207 > LOGFILE/log_4.1_v1_PH207.txt &

# Run 5.1
cd /workdir/rt475/TWAS_2020-07
nohup R --vanilla --slave < 5.1-CombinedTest.R --args v1.1 B73 > LOGFILE/log_5.1_v1.1_B73.txt &
nohup R --vanilla --slave < 5.1-CombinedTest.R --args v1 B73 > LOGFILE/log_5.1_v1_B73.txt &
nohup R --vanilla --slave < 5.1-CombinedTest.R --args v1.1 PH207 > LOGFILE/log_5.1_v1.1_PH207.txt &
nohup R --vanilla --slave < 5.1-CombinedTest.R --args v1 PH207 > LOGFILE/log_5.1_v1_PH207.txt &

########## NEED RERUN ##########




########## NOT DONE ##########
