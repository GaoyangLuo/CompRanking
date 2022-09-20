#!/usr/bin/env python
# title             :combined.py
# description       :process hmm output for virulence prediction
# author            :Gaoyang Luo
# date              :202209017
# version           :1.0
# usage             :import combined
# required packages :re, pandas, numpy 
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import re
import glob
import os
import cpr_run as cpr

class Combined():
    #preset feature
    def __init__(self, rgiInputDir, dvfInputDir,plasflowInputDir,seekerInputDir, VFInputDir):
        self.rgi_input_dir=rgiInputDir
        self.dvf_input_dir=dvfInputDir
        self.plasflow_input_dir=plasflowInputDir
        self.seeker_input_dir=seekerInputDir
        self.virulence_dir=VFInputDir
        
   
        
        
        