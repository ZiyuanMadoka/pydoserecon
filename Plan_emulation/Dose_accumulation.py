import os
import pydicom
import sys
import numpy as np
from math import log10, floor
import glob
from rtplanchanger import *
# import matplotlib.pyplot as plt
# import datetime

dose_dir = "../data/auto_Dose/"
d_plan = '../data/autoPlans/'
def printpage():
    print "Script for accumalating dose files saved at ../auto_Dose"
    print "./Dose_accumulation.py refplan surID"
    print "output dose file name is RD.00.ref+(refplan).(surID).dcm"
    quit()

def round_sig(x,sig=2):
	return round(x, sig-int(floor(log10(abs(x))))-1)

def Dose_accumulation(refplan,surID,d_plan,dose_dir=dose_dir):
    RD_pt="RD.ref"+str(refplan)+"."+str(surID)+".*" # dose file name pattern
    dose_array=0
    RDfiles = glob.glob(dose_dir+RD_pt)
    #print "dose find!"
    if len(RDfiles)!=0:
        for file in RDfiles:
            dose = pydicom.dcmread(file)
            # Dose Grid Scaling: when multiplied by the dose grid data, yields grid doses specified by dose units (3004, 0002), in our case Gy 
            dose_array = dose.pixel_array * dose.DoseGridScaling + dose_array
    else:
        print("cannot find dose file "+RD_pt + "*. dcm")
        print("Please check the input patient IDs: (institute) + (_)+ (digit) ")
        raise ValueError('DoseCalc failed')

    new_dose = dose
    
    print "Accumulated maximum dose based on original MU is", round_sig(dose_array.max(),4), "Gy"
    isoc = dose_ratio_max_isoc(refplan, surID,d_plan)
    dose_pos = np.squeeze(np.array(new_dose.ImagePositionPatient)).astype(float)
    slice_thickness = dose.GridFrameOffsetVector[1]-dose.GridFrameOffsetVector[0]
    pixel_spacing = dose.PixelSpacing
    #print slice_thickness, pixel_spacing
    dose_idx = np.array([(isoc[0] - dose_pos[0])/pixel_spacing[0], (isoc[1]-dose_pos[1])/pixel_spacing[1], (isoc[2]-dose_pos[2])/slice_thickness]).astype(int)
    dose_isoc = dose_array[dose_idx[2],dose_idx[1],dose_idx[0]]
    dose_slice = dose_array[dose_idx[2],:]
    # print dose_slice.shape, dose_idx[1], dose_idx[0]
    # plt.imshow(dose_slice)
    # plt.show()
    print 'dose isoc: ', dose_isoc
    #change the max dose to keep isoc dose as 14.4 Gy
    dose_array = dose_array/dose_isoc*14.4 
    max_dose = dose_array.max()

    print "scale dose grid to deliver 14.4 Gy to isocenter...new maximum dose", round_sig(max_dose,4)
    if max_dose > 0:
        new_DoseGridScaling = max_dose / (2 **16-1) # to make sure after the rounded value , it is still less than 2 ** 16 -1
        new_DoseGridScaling = round_sig(new_DoseGridScaling, 7) # rounded to keep 7 significant decimal
        new_dose.DoseGridScaling = str(new_DoseGridScaling)
        dose_array = dose_array/new_DoseGridScaling
    
    new_pixel_array = np.array(dose_array,dtype=np.uint16)
    
    new_dose.DoseSummationType = "PLAN"
    new_dose.PixelData = new_pixel_array.tobytes()
    #new_private= round_sig(new_dose.DoseGridScaling/766.8,7)  #what is 766.8?
    #print(new_private)
    #dose_new[0x3005,0x1023] = str(new_private)
    #test = list(dose_new[0x3005,0x1023])
    #print test
    new_dose.SOPInstanceUID = new_dose.SOPInstanceUID[:-5]+'12345' 
    del new_dose[0x300c,0x0002][0][0x300c,0x0020]  # delete referenced Fraction Group Sequence, because the summation type changed.:
    outfile = dose_dir + "RD.00.ref"+refplan+'.'+surID+".dcm"
    new_dose.save_as(outfile)
    print("Finish dose acculation. Saved as "+outfile)

def dose_ratio_max_isoc(refplan, surID,d_plan):
    file_plan = d_plan + "RP.ref"+str(refplan) + '.'+str(surID) + ".dcm"
    features_txt = planfeatures(file_plan)
    line0 = ''
    planf = {}
    #print features_txt
    for line in features_txt.split('\n'):
        if 'isocenter' in line:
            line = line.replace(']','')
            line = line.replace('[','')
            line = line.split()
            isocstr = line[1:]
                   
    isoc = np.squeeze(np.array(isocstr)).astype(float)
    return isoc


def main():
    if len(sys.argv) == 1: # if no input parameter 
        printpage()
    elif ( len(sys.argv) == 2 and sys.argv[1] == '-help' ): #if user ask for -help
        printpage()
    refplan = sys.argv[1]
    surID = sys.argv[2]
    try:
	Dose_accumulation(refplan,surID,d_plan)
    except ValueError as e:
	print e
        
if __name__ == "__main__":
    main()
