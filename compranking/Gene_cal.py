#!/usr/bin/env python
# title             :Gene_cal.py -> Gene_cal.ipynb
# description       :calculate relative abundance of genes
# author            :Gaoyang Luo
# date              :202201101
# version           :1.0
# usage             :import Gene_cal
# required packages :Bio, pandas, os
# notification: enjoy yourself
#==============================================================================
#import modules
import pandas as pd
import re
import glob
import os
import path
from Bio import SeqIO
