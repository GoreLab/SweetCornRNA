#!/bin/bash
nohup R --vanilla --slave < 4.3-Peer_UseOptFact.R --args BLUE lmer 13 > LOGFILE/out_431.txt &
nohup R --vanilla --slave < 4.3-Peer_UseOptFact.R --args BLUP lmer 15 > LOGFILE/out_432.txt &