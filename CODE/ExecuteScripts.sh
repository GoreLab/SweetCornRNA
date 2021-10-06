cd /workdir/rt475/Sweetcorn
mkdir LOGFILE

# Run 1.1
nohup R --vanilla --slave < 1.1-LookTagSeq_MakeData.R --args B73 > LOGFILE/log_1.1_B73.txt &

# Run 1.2
nohup R --vanilla --slave < 1.2-GeneFiltering.R --args B73 0.5 > LOGFILE/log_1.2_B73.txt &

# Run 1.3
nohup R --vanilla --slave < 1.3-OutlierRemoval_and_Imputation.R --args B73 100 > LOGFILE/log_1.3_B73.txt &

# Run 3.4 -- this part is done by other shell scripts
nohup sh 3.4-RunMultipleR_B73_NoHarvDate.sh > LOGFILE/log_3.4_B73.txt &

# Run 3.5
nohup R --vanilla --slave < 3.5-BLUE_merge_NoHarvDate.R --args B73 > LOGFILE/log_3.5_B73.txt &

# Run 3.6
nohup R --vanilla --slave < 3.6-BLUE_RmErr_NoHarvDate.R --args B73 > LOGFILE/log_3.6_B73.txt &

# Run 4.2
nohup R --vanilla --slave < 4.2-MakeBlueDatasets_NoHarvDate.R --args B73 > LOGFILE/log_4.2_B73.txt &

# Run 5.5
nohup R --vanilla --slave < 5.5-Peer_Use25Fact_NoHarvDate.R --args B73 toco > LOGFILE/log_5.5_B73.txt &
nohup R --vanilla --slave < 5.5-Peer_Use25Fact_NoHarvDate.R --args B73 carot > LOGFILE/log_5.5_B73.txt &

# Run 5.6 (this script make plot for all)
nohup R --vanilla --slave < 5.6-MakeScreePlot_NoHarvDate.R > LOGFILE/log_5.6.txt &

# Run 5.7
nohup R --vanilla --slave < 5.7-Peer_UseOptFact_NoHarvDate.R --args B73 toco 17 > LOGFILE/log_5.7_B73.txt &
nohup R --vanilla --slave < 5.7-Peer_UseOptFact_NoHarvDate.R --args B73 carot 16 > LOGFILE/log_5.7_B73.txt &

# Run 5.8
nohup R --vanilla --slave < 5.8-OutlierRemoval_NoHarvDate.R --args B73 toco > LOGFILE/log_5.8_B73.txt &
nohup R --vanilla --slave < 5.8-OutlierRemoval_NoHarvDate.R --args B73 carot > LOGFILE/log_5.8_B73.txt &

# Run 6.1 & 6.2 -- I DID NOT RUN ON THE LINUX SERVER. I RAN IT IN LOCAL ENV OF MY R STUDIO (R version 4.0.2; GAPIT3)
# if execute it, run like this: nohup R --vanilla --slave < 6.1-CalcKinship_GAPIT.R &
# same for 6.2

# Run 7.3
nohup R --vanilla --slave < 7.3-TWAS_NoHarvDate.R --args B73 toco > LOGFILE/log_7.3_B73_toco.txt &
nohup R --vanilla --slave < 7.3-TWAS_NoHarvDate.R --args B73 carot > LOGFILE/log_7.3_B73_carot.txt &

# Run 9.6 (kernel mutant type, use BLUE/PEER without using harvest date, without the double mutant)
nohup R --vanilla --slave < 9.6-ModelComparison_BIC_mutant_type_ReAnalysis.R --args toco > LOGFILE/log_9.6_toco.txt &
nohup R --vanilla --slave < 9.6-ModelComparison_BIC_mutant_type_ReAnalysis.R --args carot > LOGFILE/log_9.6_carot.txt &

# Run 9.7
nohup R --vanilla --slave < 9.7-TWAS_Kernel_mutant_type_ReAnalysis.R --args B73 toco > LOGFILE/log_9.7_B73_toco.txt &
nohup R --vanilla --slave < 9.7-TWAS_Kernel_mutant_type_ReAnalysis.R --args B73 carot > LOGFILE/log_9.7_B73_carot.txt &

################################################################################
# for the new RNAseq data
nohup R --vanilla --slave < 1.1-LookTagSeq_MakeData_Ia453.R > LOGFILE/log_1.1_Ia453.txt &
nohup R --vanilla --slave < 1.2-GeneFiltering.R --args Ia453 0.5 > LOGFILE/log_1.2_Ia453.txt &
nohup R --vanilla --slave < 1.3-OutlierRemoval_and_Imputation.R --args Ia453 100 > LOGFILE/log_1.3_Ia453.txt &
nohup sh 3.4-RunMultipleR_Ia453_NoHarvDate.sh > LOGFILE/log_3.4_Ia453.txt &
nohup R --vanilla --slave < 3.5-BLUE_merge_NoHarvDate.R --args Ia453 > LOGFILE/log_3.5_Ia453.txt &
nohup R --vanilla --slave < 3.6-BLUE_RmErr_NoHarvDate.R --args Ia453 > LOGFILE/log_3.6_Ia453.txt &
nohup R --vanilla --slave < 4.2-MakeBlueDatasets_NoHarvDate.R --args Ia453 > LOGFILE/log_4.2_Ia453.txt &
nohup R --vanilla --slave < 5.5-Peer_Use25Fact_NoHarvDate.R --args Ia453 toco > LOGFILE/log_5.5_Ia453.txt &
nohup R --vanilla --slave < 5.5-Peer_Use25Fact_NoHarvDate.R --args Ia453 carot > LOGFILE/log_5.5_Ia453.txt &
nohup R --vanilla --slave < 5.6-MakeScreePlot_NoHarvDate.R > LOGFILE/log_5.6.txt &
nohup R --vanilla --slave < 5.7-Peer_UseOptFact_NoHarvDate.R --args Ia453 toco 15 > LOGFILE/log_5.7_Ia453.txt &
nohup R --vanilla --slave < 5.7-Peer_UseOptFact_NoHarvDate.R --args Ia453 carot 16 > LOGFILE/log_5.7_Ia453.txt &
nohup R --vanilla --slave < 5.8-OutlierRemoval_NoHarvDate.R --args Ia453 toco > LOGFILE/log_5.8_Ia453.txt &
nohup R --vanilla --slave < 5.8-OutlierRemoval_NoHarvDate.R --args Ia453 carot > LOGFILE/log_5.8_Ia453.txt &
nohup R --vanilla --slave < 7.3-TWAS_NoHarvDate.R --args Ia453 toco > LOGFILE/log_7.3_Ia453_toco.txt &
nohup R --vanilla --slave < 7.3-TWAS_NoHarvDate.R --args Ia453 carot > LOGFILE/log_7.3_Ia453_carot.txt &
nohup R --vanilla --slave < 9.7-TWAS_Kernel_mutant_type_ReAnalysis.R --args Ia453 toco > LOGFILE/log_9.7_Ia453_toco.txt &
nohup R --vanilla --slave < 9.7-TWAS_Kernel_mutant_type_ReAnalysis.R --args Ia453 carot > LOGFILE/log_9.7_Ia453_carot.txt &
nohup R --vanilla --slave < 7.3-TWAS_NoHarvDate_Ia453_aT.R > LOGFILE/log_7.3_aT.txt &
nohup R --vanilla --slave < 9.8-MakeManhattan_kernel_mutant_type_ReAnalysis.R > LOGFILE/log_9.8.txt &

################################################################################
# stratified sampling
# Run 13.1
module load R/3.5.0
nohup R --vanilla --slave < 13.1-Prediction_stratified.R --args toco UseMu > LOGFILE/log_13.1_toco_UseMu.txt &
nohup R --vanilla --slave < 13.1-Prediction_stratified.R --args carot UseMu > LOGFILE/log_13.1_carot_UseMu.txt &
nohup R --vanilla --slave < 13.1-Prediction_stratified.R --args toco NoMu > LOGFILE/log_13.1_toco_NoMu.txt &
nohup R --vanilla --slave < 13.1-Prediction_stratified.R --args carot NoMu > LOGFILE/log_13.1_carot_NoMu.txt &

# Run 13.2
module load R/3.5.0
nohup R --vanilla --slave < 13.2-Prediction_stratified_useCand.R --args toco UseMu > LOGFILE/log_13.2_toco_UseMu.txt &
nohup R --vanilla --slave < 13.2-Prediction_stratified_useCand.R --args carot UseMu > LOGFILE/log_13.2_carot_UseMu.txt &
nohup R --vanilla --slave < 13.2-Prediction_stratified_useCand.R --args toco NoMu > LOGFILE/log_13.2_toco_NoMu.txt &
nohup R --vanilla --slave < 13.2-Prediction_stratified_useCand.R --args carot NoMu > LOGFILE/log_13.2_carot_NoMu.txt &
