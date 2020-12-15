#!/bin/bash
nohup R --vanilla --slave < 4.1-Peer_Use25Fact.R --args BLUE lmer > LOGFILE/out_411.txt &
nohup R --vanilla --slave < 4.1-Peer_Use25Fact.R --args BLUP lmer > LOGFILE/out_412.txt &