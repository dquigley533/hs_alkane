#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 13:39:16 2021

@author: phseal
"""

import hs_alkane.alkane as alk
import time


#########################
#  Testing timer module #
#########################
alk.timer_init() # initialise timer
time.sleep(1)

alk.cvar.timer_qtime = 10    # 10 seconds of allowed execution time
alk.cvar.timer_closetime = 1 #  5 seconds needed to stop safely

safe = alk.timer_check_continuation()
if (safe == 1):
    print("Safe to continue")


##########################
#  Testing random module #
##########################
seed = 234289
print("Seeding Fortran RNG with seed = ",seed)
alk.random_set_random_seed(seed)
xi = alk.random_uniform_random()
print("Random number xi = ",xi)


