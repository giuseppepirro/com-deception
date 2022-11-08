#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 20:08:44 2022

@author: saifaldeen
"""

import pandas as pd

df = pd.read_csv("./freeman.txt", delimiter=" ")

print(df)

df.to_csv("./freeman.csv", index=False)
