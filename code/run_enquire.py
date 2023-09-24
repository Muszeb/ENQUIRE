


import os
import sys
import numpy as np 
import pandas as pd
import subprocess
import glob 
import re
import bz2

# ===================================

# import ENQUIRE
from ENQUIRE import run

# ------

# output_=ENQUIRE.run('SplicingFactors_Neoplasms_Antigens', 'pmid-ferroptosis_ImmuneSystem.txt',      1,     "/Users/surya/Documents/GITHUB-REPOSITORIES/ENQUIRE/",     4,2,3,"all",6)
output_=run('SplicingFactors_Neoplasms_Antigens', 'pmid-ferroptosis_ImmuneSystem.txt',      1,     "/Users/surya/Documents/GITHUB-REPOSITORIES/ENQUIRE/",     4,2,3,"all",6)

# ===================================

print(output_)