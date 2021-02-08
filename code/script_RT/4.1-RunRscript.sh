#!/bin/bash
cd /workdir/rt475/Sweetcorn
nohup R --vanilla --slave < 4.1-Peer_Use25Fact.R --args BLUE lmer > LOGFILE/out_41.txt &
