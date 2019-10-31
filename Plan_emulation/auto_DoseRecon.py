##! /usr/bin/python
# auto_DoseRecon.py
"""reconstruct RT plan on a surrogate CT and calculate the DVH metrics."""
# Copyright (c) 2018- Ziyuan
# Copyright (c) 2018- Marco Virgolin

import os
import sys
import shutil
import threading
import time

import pydicom
import numpy as np
from math import log10, floor
import glob
from dvhcalc import * 
import logging
from rtplanchanger import *
from auto_DoseCalc import *
from Dose_accumulation import *
 
d_ref = '../data/ref_plans/'
d_plan = '../data/autoPlans/'
d_meas = '../data/auto_measures/'
d_patients = '../data/CT_data/'
d_RS = '../data/data_RS/'
d_DVH = '../data/dvh/'
d_dose = '../data/auto_Dose/'


def printpage():
    print "---------------Help info--------------"
    print "Script for auto dose reconstruction!"
    print "python auto_DoseRecon.py refplan surID"
    print "..."
    print "default: dir_ref: ../ref_plans/; dir_plan: ../autoPlans/; dir_patients: ../2D3D_data/; dir_DVH = ../dvh/"
    print "optional (set options at the end): -dir_ref PATH (end with /) -dir_plan PATH -dir_patients PATH -dir_RS PATH -dir_DVH PATH"
    quit()

if not os.path.exists(d_dose):
	os.makedirs(d_dose)

if len(sys.argv) == 1: #if no parameter
    printpage()
else:
    if ( len(sys.argv) == 2 and sys.argv[1] == '-help' ): #if user ask for -help
        printpage()


# set custom options if present
for i in range(4, len(sys.argv)):
    if sys.argv[i] == '-dir_plan':
        d_plan = sys.argv[i+1]
    elif sys.argv[i] == '-dir_ref':
        d_ref = sys.argv[i+1]
    elif sys.argv[i] == '-dir_patients':
        d_patients = sys.argv[i+1]
    elif sys.argv[i] == '-dir_RS':
        d_RS = sys.argv[i+1]
    elif sys.argv[i] == '-dir_DVH':
        d_DVH = sys.argv[i+1]


def DoseRecon(refplan,surID):
    if "block" in refplan:
        refID = refplan[0:refplan.find('.')]
        block = refplan[refplan.find('.')+1:len(refplan)]+'.'
    elif "sampled" in refplan:
        refID = refplan[0:refplan.find('_sampled')]    
    else:
        refID = refplan
        block = ''
    refRP = d_ref + 'RP.ref' + refplan +'.'+refID+ '.dcm'
    refRS = d_RS + 'RS.' + refID + '.dcm'
    surRS = d_RS + 'RS.' + surID + '.dcm'
    measure_ref = d_meas + 'measurement_' + refID + '.txt'
    measure_sur = d_meas + 'measurement_' + surID + '.txt'
    IDpair = refplan+'.'+surID
    autoRP = 'RP.ref' + IDpair + '.dcm'
    surPatient = d_patients+surID+'/'
    RD_pt = "RD.ref"+IDpair+".*" 
    rtdosefile = d_dose + "RD.00.ref" + IDpair + ".dcm"
   
    try:
        rtplangenerator(refRP, refRS, measure_ref, surRS, measure_sur, d_plan+autoRP)
        DoseCalc(surPatient, IDpair,d_plan) # this function needs the Dose Engine to be stalled in your computer !!!
        print "Dose calculation finished!"
    except ValueError as e:
        print e
        shutil.copyfile(d_plan+autoRP,'../autoPlans_error/'+autoRP)
    try:
        print "Start dose accumulation ..."
        Dose_accumulation(refplan,surID,d_plan,d_dose) 
    except ValueError as e:
        print e

    outputfile = d_DVH + 'DVHmetrics_plan' + refplan + '.txt' 
    calDVHmetrics(surRS,rtdosefile,outputfile)
    
    print 'save dose dvh for ', refplan, 'on ', surID, ' at', outputfile 
    return 0


# reconstruct the specific reference plan (ID) and surrogate patient (ID)
refplan = sys.argv[1]
surID = sys.argv[2]
try:
	DoseRecon(refplan,surID)
except ValueError as e:
	print e



