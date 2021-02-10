#!/bin/bash
cd /workdir/rt475/Sweetcorn
nohup R --vanilla --slave < 4.3-Peer_UseOptFact.R --args BLUE lmer 12 > LOGFILE/out_43.txt &
