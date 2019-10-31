##! /usr/bin/python
# plot_RTplan.py
"""plot RT plan on the corresponding DRR with multiple options"""
# Copyright (c) 2018- Ziyuan
# Copyright (c) 2018- Marco Virgolin
import os
import sys
import shutil
import threading
import time
import pydicom
import numpy as np
import math
import glob
import logging
from rtplanchanger import *
#from medpy.io import load
import SimpleITK as sitk
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import subprocess
import shlex
from shapely.geometry import Point, LineString, Polygon
import shapely.ops

def plotmhd(file_DRR_AP, title=None, margin = 0.0, dpi = 300):
    #plot the mhd DRR and return the ax and plt
    itkimage = sitk.ReadImage(file_DRR_AP)    
    nda = sitk.GetArrayFromImage(itkimage)
    spacing = itkimage.GetSpacing()
    xsize = nda.shape[1]
    ysize = nda.shape[0]
    figsize = (1+margin) * xsize / dpi, (1+margin) * ysize /dpi
    fig = plt.figure(figsize=figsize,dpi=dpi)
    #print "image pixel numbers", xsize, ysize dpi:dots per inch
    #print "figure size", figsize
    #ax = fig.add_subplot(111)
    ax = plt.Axes(fig,[0,0,1,1], frameon=False)
    fig.add_axes(ax)
 
    
    extent= (0,xsize * spacing[0], ysize *spacing[1],0) # the same orientation as DRR
    #print "extent", extent
    t = ax.imshow(nda,extent=extent,interpolation=None)
    if nda.ndim ==2:
        t.set_cmap("gray")
    
    if(title):
        plt.tile(title)
    #plt.show()
    return ax, fig
    
def rotate(x,y,theta_deg):
    #apply rotation matrix to vector that rotates by theta contour-lockwise
    theta = theta_deg *math.pi/180
    x_new = math.cos(theta) * x - math.sin(theta) * y
    y_new = math.sin(theta) * x + math.cos(theta) * y
    return x_new, y_new

def append_point(x,y,point):
    # append the point's x-y coordinates to the x,y list
    x_new = np.append(x,point[0])
    y_new = np.append(y,point[1])
    return x_new, y_new

def insert_point(x,y,point):
    # insert the point's x-y coordinates to first of the x,y list
    x_new = np.insert(x,0,point[0])
    y_new = np.insert(y,0,point[1])
    return x_new, y_new

def add_border(x, y, point1, point2):
    #if x y is empty, add point1/point2 as borders
    if len(x):
        return x, y
    else:
        x_new, y_new = append_point(x,y,point1)
        x_new, y_new = insert_point(x_new,y_new,point2)
    return x_new, y_new

def CT_coordinate(isoc,point):
    # translate to CT coordinate system
    x = isoc[0]+point[0]
    y = isoc[1]-point[1]
    return x, y

def limit_in_field(borders, positions, Fbounds):
    #print Fbounds, borders[0], borders[-1]
    #print 'before', positions

    if borders[0] != Fbounds[0]:
        positions[0] = positions[1]+(Fbounds[0]-borders[1])*(positions[2]-positions[1])/10
        #positions[0] = positions[0]+(Fbounds[0]-borders[0])*(positions[1]-positions[0])/10
        borders[0] = Fbounds[0]
        
    if borders[-1] != Fbounds[1]:
        positions[-1] = positions[-2] + (Fbounds[1]-borders[-2])*(positions[-2] - positions[-3])/10
        borders[-1] = Fbounds[1]
    
    #print 'after', positions, borders
    return borders, positions

def planplot(file_plan,surID,d_fig, d_DRR, mode="mask",showaxis = "off", updateDRR = "N"):
    if not file_plan.endswith('.dcm'):
        print('What? File',file_plan,'is not a .dcm file!')
        return

    # plot plan field on top of DRR
    features_txt = planfeatures(file_plan)
    line0 = ''
    planf = {}
    #print features_txt
    for line in features_txt.split('\n'):
        if ']' in line:
            line = line0+line
            line = line.replace(']','')
            line = line.replace('[','')
            line = line.split()
            planf[line[0]] = line[1:]
            line0 = ''
        else:
            line0 = line0+line  
    isoc = np.squeeze(np.array(planf["isocenter"])).astype(float)
    if mode == "shape":
        set_DRR_string = ' -isoc '
    else:
        set_DRR_string = ' -bone_enh 1 -window 3000 -level 1300 -isoc '
    #command_DRR = '/U_Disk/Algemeen/ZWG/Pat_pydicom/RTplan_changer/gen_DRRs.sh FILE ' + surID +' -crop N -RS N -dir_DRR '+d_DRR+set_DRR_string +str(isoc[0]) +' '+ str(isoc[1]) +' '+ str(isoc[2])
    # the commented line requires the DRR generation software to be intalled in your pc
    
    file_DRR_AP = d_DRR+'DRR' + surID + '.mhd'
    #generate new DRR based on plan isocenter when updateDRR flag is "Y"
    if updateDRR == "Y":
        try:
            subprocess.check_call(command_DRR,shell=True,env=my_env)
        except subprocess.CalledProcessError as err:
            print('ERROR:', err)
    
    assert os.path.isfile(file_DRR_AP), "mhd file does not exist"
    itkimage = sitk.ReadImage(file_DRR_AP)
  
    origin = itkimage.GetOrigin()
    spacing = itkimage.GetSpacing()
    ax, fig = plotmhd(file_DRR_AP)
    theta = np.squeeze(planf["CollimatorAngle"]).astype(float)
    FX = np.squeeze(np.array(planf["field_borders_x"])).astype(float)
    FZ = np.squeeze(np.array(planf["field_borders_z"])).astype(float)
    LU_point = rotate(FX[0], FZ[1], theta) # vector of field left upper point relative to isocenter in CT system
    LU_point0 = LU_point
    LB_point = rotate(FX[0], FZ[0], theta) # vector of field left bottom point relative to isocenter in CT system
    RU_point = rotate(FX[1], FZ[1], theta) # vector of field right upper point relative to isocenter in CT system
    RB_point = rotate(FX[1], FZ[0], theta) # vector of field right bottom point relative to isocetner in CT system 
    isoc_loc = np.array([(isoc[0]-origin[0]),(origin[1]-isoc[2])])
    #print "isoc_loc in DRR", isoc_loc, LU_point0, LB_point, RU_point
    Leftborder_x = []
    Rightborder_x = []
    Leftborder_y = []
    Rightborder_y = []
    Upperborder_x = []
    Bottomborder_x = []
    Upperborder_y = []
    Bottomborder_y = []
    if "field_MLC_x" in planf:
        MLC_x = np.squeeze(np.array(planf["field_MLC_x"])).astype(float)
        ys_L = np.array(10*(np.arange(40)-20)+5).astype(float)
        ys = np.append(ys_L,ys_L)
         # rotates MLC vectors relative to isocenter back to  CT system (same orientation)
        idxL = range(0,40)
        idxR = range(40,80)
        if mode == "mask" or mode == "shape":
            idxL = [i for i in idxL if ys[i] >= FZ[0]-5 and ys[i] <= FZ[1]+5 and MLC_x[i] >= FX[0] and MLC_x[i] <= FX[1]]
            idxR = [i for i in idxR if ys[i] >= FZ[0]-5 and ys[i] <= FZ[1]+5 and MLC_x[i] >= FX[0] and MLC_x[i] <= FX[1]]
        ys[idxL], MLC_x[idxL] = limit_in_field(ys[idxL],MLC_x[idxL],FZ)
        #print MLC_x[idxR]
        ys[idxR], MLC_x[idxR] = limit_in_field(ys[idxR],MLC_x[idxR],FZ)
        #print MLC_x[idxR]
        MLC_xnew,ys_new = rotate(MLC_x,ys,theta)
        
        if idxL:
            Leftborder_x = np.array(isoc_loc[0]+MLC_xnew[idxL])
            Leftborder_y = np.array(isoc_loc[1]-ys_new[idxL])
            LU_point = [MLC_xnew[idxL[-1]], LU_point[1]]
            LB_point = [MLC_xnew[idxL[0]], LB_point[1]]

        
        if idxR:
            Rightborder_x = np.array(isoc_loc[0]+MLC_xnew[idxR])
            Rightborder_y = np.array(isoc_loc[1]-ys_new[idxR])
            RU_point = [MLC_xnew[idxR[-1]], RU_point[1]]
            RB_point = [MLC_xnew[idxR[0]], RB_point[1]]

    #if no MLCx, add the field borders to the MLCx border arrays
    Leftborder_x, Leftborder_y = add_border(Leftborder_x, Leftborder_y, CT_coordinate(isoc_loc,LU_point),CT_coordinate(isoc_loc,LB_point))
    Rightborder_x, Rightborder_y = add_border(Rightborder_x, Rightborder_y, CT_coordinate(isoc_loc,RU_point),CT_coordinate(isoc_loc,RB_point))

    if "field_MLC_y" in planf:
        MLC_y = np.squeeze(np.array(planf["field_MLC_y"])).astype(float)
        xs_U = np.array(10*(np.arange(40)-20)+5).astype(float)
        xs = np.append(xs_U,xs_U)  
        idyU = range(0,40)
        idyB = range(40,80)
        
        border_y = np.array(isoc_loc[0]+xs_new[idxL])
        if mode == "mask" or mode == "shape":
            idyU = [i for i in idyL if xs[i] >= FX[0]-5 and xs[i] <= FX[1]+5 and MLC_y[i] >= FY[0] and MLC_y[i] <= FY[1]]
            idyB = [i for i in idyR if xs[i] >= FX[0]-5 and xs[i] <= FX[1]+5 and MLC_y[i] >= FY[0] and MLC_y[i] <= FY[1]]   
            
        xs[idxU], MLC_y[idyU] = limit_in_field(xs[idxU], MLC_y[idyU],FX)
        xs[idxB], MLC_y[idyB] = limit_in_field(xs[idxB], MLC_y[idyB],FX)

        xs_new, MLC_ynew = rotate(xs,MLC_y,theta)  # rotates MLC vectors relative to isocenter back to  CT system (same orientation)
        if idyU:
            Upperborder_x = np.array(isoc_loc[0]+xs_new[idyU])
            Upperborder_y = np.array(isoc_loc[1]-MLC_ynew[idyU])
            LU_point = [LU_point[0],MLC_ynew[idyU[0]]]
            RU_point = [RU_point[0], MLC_ynew[idyU[-1]]]
            
        if idyB:
            Bottomborder_x = np.array(isoc_loc[0]+xs_new[idyB])
            Bottomborder_y = np.array(isoc_loc[1]-MLC_ynew[idyB])
            LB_point = [LB_point[0],MLC_ynew[idyB[0]]]
            RB_point = [RB_point[0], MLC_ynew[idyB[-1]]]
    
    x = np.concatenate((Leftborder_x,np.flipud(Upperborder_x),np.flipud(Rightborder_x),Bottomborder_x),axis=0) # upperborder was from right to left, need to be flipped

    y = np.concatenate((Leftborder_y,np.flipud(Upperborder_y), np.flipud(Rightborder_y),Bottomborder_y),axis=0)

    ls = LineString(np.c_[x,y])
    # make line rings of the borders
    lr = LineString(ls.coords[:] + ls.coords[0:1])

    
    mls = shapely.ops.unary_union(lr)
    
    max_area = float("-inf")
    # if the polygon is self-intersected, select the one with largest area:
    for polygon in shapely.ops.polygonize(mls):
        if polygon.area > max_area:
            max_area = polygon.area
            ef_plg = polygon

    x_test, y_test = ef_plg.exterior.coords.xy
      

    ax.plot(x_test,y_test, linewidth=0.5) #, color = color )
        #ax.plot(isoc_loc[0]+MLC_xnew,isoc_loc[1]-ys_new,linewidth=0.5) # to check the original MLCs including outside the field
    rect = patches.Rectangle((isoc_loc[0]+LU_point0[0], isoc_loc[1]-LU_point0[1]),np.sum(np.absolute(FX)),np.sum(np.absolute(FZ)),angle=-theta,linewidth=0.5,edgecolor='y',facecolor='none')
        #  we use -theta here because the DRR coordinate system is up-donw flipped compared to the CT system, so rotation appears contour-clock
        #  the xy value in rectangle function is the left upper,point of the field 

    if showaxis == "off":
        plt.axis('off')

    ax.plot(isoc_loc[0],isoc_loc[1],'ro',markersize=2)
    
    if mode == "shape":
        print("effective field shape!")
        ax.add_patch(rect)
    #plt.show()
    #print file_plan
    plantxt =file_plan.split('RP.')[1]
    figtitle = d_fig+'RP.'+ plantxt.split('.dcm')[0] + '.png'
    # figtitle = d_fig+'DRR.'+ plantxt.split('.dcm')[0] + '.png' # temp
    print 'save figure to ', figtitle
    fig.savefig(figtitle)
    plt.close()


def printpage():
    print "---------------Help info--------------"
    print "Script for plot effective field shape on DRR! Mode: FILE or DIR. Two types of usage:"
    print "python plot_RT.py FILE refID surID (optional)"
    print "python plot_RT.py DIR (optional)"
    print "..."
    print "default: type: mask dir_plan: ../autoPlans/; dir_DRRs: ../auto_DRRs/; dir_fig: ../plotted_field"
    print "optional: -type shape -dir_plan PATH(end with /); -dir_DRR PATH; -dir_fig PATH -updateDRR Y"
    quit()

def main():
    d_plan = '../data/autoPlans/'
    d_DRR = '../data/auto_DRRs/'
    d_fig = '../data/plotted_field/'
    mode = "shape"
    updateDRR = "N"

    if len(sys.argv) == 1: #if no parameter
        printpage()
    else:
        if ( len(sys.argv) == 2 and sys.argv[1] == '-help' ): #if user ask for -help
            printpage()
    # set custom options if present
    for i in range(2, len(sys.argv)):
        if sys.argv[i] == '-dir_plan':
            d_plan = sys.argv[i+1]
        elif sys.argv[i] == '-type':
            mode = sys.argv[i+1]
        elif sys.argv[i] == '-dir_DRR':
            d_DRR = sys.argv[i+1]
        elif sys.argv[i] == '-updateDRR':
            updateDRR = sys.argv[i+1]
        elif sys.argv[i]=='-dir_fig':
            d_fig = sys.argv[i+1]

    if not os.path.exists(d_fig):
        os.makedirs(d_fig)
    # FILE mode
    # reconstruct the specific reference plan (ID) and surrogate patient (ID)
    if sys.argv[1] == 'FILE':
        refplan = sys.argv[2]
        surID = sys.argv[3]
        autoRP = 'RP.ref' + refplan +'.'+surID+ '.dcm'
        planplot(os.path.join(d_plan+autoRP),surID,d_fig = d_fig, d_DRR = d_DRR, mode = mode,updateDRR=updateDRR)
      
    # DIR mode
    # calculate dose for all the auto plans saved dir_plan on correspondent surrogate patients 
    elif sys.argv[1] == 'DIR':
        all_plans = os.listdir(d_plan)
        # figs = os.listdir(d_fig)
        # paths = [os.path.join(d_fig, basename) for basename in figs]
        # print len(figs)
        # print all_plans.index('RP.refAMC_2_sampled_52.UMCU_35.dcm')
        # print max(paths, key=os.path.getctime)
        # quit()
        for autoRP in all_plans:
            tempIDs = autoRP[ autoRP.rfind('ref')+3 : autoRP.rfind('.dcm') ]
            refplan = tempIDs[0: tempIDs.rfind('.')]
            surID = tempIDs[tempIDs.rfind('.')+1: len(tempIDs)]
            # skip non-.dcm files                                                                           
            try:
                planplot(d_plan+autoRP,surID,d_fig = d_fig, d_DRR = d_DRR, mode = mode,updateDRR=updateDRR)
            except ValueError as e:
                print e
                continue
            except AssertionError as e:
                print e 
                continue
            
        print 'Finished plotting! '

    else:
        print "Wrong mode, please specify between FILE and DIR: ..."

        printpage()


if __name__ == "__main__":
    main()