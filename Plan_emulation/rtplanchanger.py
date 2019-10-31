#!/usr/bin/env python
# -*- coding: utf-8 -*-
# rtplanchanger.py
"""Emulate a plan on a surrogate anatomy based on measurements from DRR."""
# Copyright (c) 2018- Ziyuan
# Copyright (c) 2018- Marco Virgolin

import pydicom
import sys
import numpy as np
import math
from ParameterHandler import ParameterHandler as PR
from TransformationParamMaker import TransformationParamMaker as TPM
from IsocAPRatios import IsocAPRatiosFinder
from shapely.geometry import box, LineString, Polygon
from shapely.ops import split
import pdb
#file_plan_ref = sys.argv[1]
#file_struct_ref = sys.argv[2]
#file_landmarks_ref = sys.argv[3]

#file_struct_sur = sys.argv[4]
#file_landmarks_sur = sys.argv[5]
#outfile_plan_sur = sys.argv[6]


#---------------test end---------------------
    
def rolled_position(roll, roll_neighbor, pos_roll_neighbor):
	pos_roll=[0 for x in range(len(roll))]
	for ii in range(0,len(roll)):
		pos_roll[ii] = pos_roll_neighbor[0] + (roll[ii]-roll_neighbor)* (pos_roll_neighbor[1]-pos_roll_neighbor[0])
	return pos_roll
def reshapemlc(old_mlc, old_bds, new_bds, old_vbds, new_vbds,tpm,xoryMLC,collim,gantry):
	num_pp = int(round(len(old_mlc)/4))
	vs_L = np.array(10*(np.arange(40)-20)+5).astype(float)
	vs = np.append(vs_L,vs_L)
	for l in range(0,4):
		temp = old_mlc[np.add(range(0,num_pp),l*num_pp)]
		temp_v = vs[np.add(range(0,num_pp),l*num_pp)]
		if l% 2 == 0:
			roll_nn = int(math.ceil((new_vbds[0]-old_vbds[0])/10))
		else:
			roll_nn = int(math.ceil((new_vbds[1]-old_vbds[1])/10))
		new_temp = np.roll(temp,roll_nn)
		idxblock = [i for i in range(0,len(new_temp)) if new_temp[i] >old_bds[0] and new_temp[i] < old_bds[1] and temp_v[i] >= new_vbds[0]-5 and temp_v[i] <= new_vbds[1]+5]
		if roll_nn<0:
			new_temp[range(roll_nn,0)] = rolled_position(range(roll_nn,0), 2*roll_nn-1, new_temp[range(2*roll_nn-1,2*roll_nn+1)])
			if len(idxblock) >2*abs(roll_nn) and l%2 == 1:
				
				blocks = new_temp[idxblock]
				blocks[range(roll_nn,0)] = rolled_position(range(roll_nn,0), 2*roll_nn-1, blocks[range(2*roll_nn-1,2*roll_nn+1)])
				new_temp[idxblock] = blocks
		elif roll_nn>0:
			new_temp[range(0,roll_nn)] = rolled_position(range(0,roll_nn), 2*roll_nn-1, new_temp[range(2*roll_nn-1,2*roll_nn+1)])
			if len(idxblock) > 2*abs(roll_nn) and l%2 == 0:
				blocks = new_temp[idxblock]
				blocks[range(0,roll_nn)] = rolled_position(range(0,roll_nn), 2*roll_nn-1, blocks[range(2*roll_nn-1,2*roll_nn+1)])
				new_temp[idxblock] = blocks
		old_mlc[np.add(range(0,num_pp),l*num_pp)] = new_temp
	
	new_mlc = np.squeeze(tpm.getNewFieldBorders(old_mlc,xoryMLC,collim,gantry))
	for mm in range(0, len(new_mlc)):
		new_mlc[mm] = str(round(new_mlc[mm],2))
	str_mlc_file = new_mlc.tolist()
	return str_mlc_file

def decompose(idxblock):
	block1 = []
	ii = 0
	block2 = []
	for ii in range(1,len(idxblock)):
		if (idxblock[ii]-idxblock[ii-1]) == 1:
			block1.append(idxblock[ii-1])
			if len(idxblock)-1 == ii:
				block1.append(idxblock[ii])
		else:
			block1.append(idxblock[ii-1])
			break

	if ii < len(idxblock)-1:
		for jj in range(ii+1,len(idxblock)):
			if (idxblock[jj]-idxblock[jj-1]) == 1:
				block2.append(idxblock[jj-1])
				if len(idxblock)-1==jj:
					block2.append(idxblock[jj])
			else:
				block2.append(idxblock[jj-1])
				break
				
	return block1, block2
def leafindex(x):
	if x!=0:
		return int(abs(x)/x * math.ceil(abs(x-5)/10)+20)
	else:
		return 0

def modelmlc(old_mlc, old_bds, new_bds, old_vbds, new_vbds,tpm,xoryMLC,collim,gantry):
	num_pp = int(round(len(old_mlc)/2))
	vs_half = np.array(10*(np.arange(40)-20)+5).astype(float)
	vs = np.append(vs_half,vs_half)
	new_mlc = np.squeeze(tpm.getNewFieldBorders(old_mlc,xoryMLC,collim,gantry)) # first scaling the MLC borders 
	for l in range(0,2):
		temp = old_mlc[np.add(range(0,num_pp),l*num_pp)]
		temp_v = vs_half
		idxblock = [i for i in range(0,len(temp)) if temp[i] >old_bds[0] and temp[i] < old_bds[1] and temp_v[i] > old_vbds[0]-5 and temp_v[i] < old_vbds[1]+5]
		block1, block2 = decompose(idxblock)
		# print block1, block2
		# print temp
		idx_in_field = [i for i in range(0,len(temp)) if temp_v[i] >= new_vbds[0]-5 and temp_v[i] <= new_vbds[1]+5]
		new_mlc[np.add(idx_in_field,l*num_pp)] = new_bds[l] # first make a rectangular field
		for block in block1, block2:
			if len(block) <2:
				continue
			elif len(block)> 2:
				scaling = tpm.getScalingSurToRef()
				if temp_v[block[0]-1]<old_vbds[0] and temp_v[block[-1]+1] < old_vbds[1]:
					block_temp = range(block[0], block[-1]+2)
				elif temp_v[block[0]-1]<old_vbds[0] and temp_v[block[-1]+1] >= old_vbds[1]:
					block_temp = range(block[0], block[-1]+1)
				elif temp_v[block[-1]+1]>old_vbds[1]:
					block_temp = range(block[0]-1,block[-1]+1)
				else:
					block_temp = range(block[0]-1,block[-1]+2)
				params= np.polyfit(temp_v[block_temp],temp[block_temp],2)
				#print "this is a big block!", params, block, block_temp, temp[block_temp],temp_v[block_temp],temp
				# fit a quadrant curve of y-x; (MLC borders is a function of leaf vertical position)
				new_params = [ scaling[0]/(scaling[1]**2)*params[0], scaling[0]/scaling[1]*params[1],scaling[0]*params[2]]
				#new_params = [ scaling[0]/(scaling[1]**3)*params[0], scaling[0]/(scaling[1]**2)*params[1], scaling[0]/scaling[1]*params[2],scaling[0]*params[3]]
				p = np.poly1d(new_params) # the quadrant curve after scaling in both directions
				#p0 = np.poly1d(params)
				new_block = range(leafindex(temp_v[block_temp[0]]*scaling[1]),leafindex(temp_v[block_temp[-1]]*scaling[1])+1)
			else:
				scaling = tpm.getScalingSurToRef()
				block_temp = range(block[0]-1,block[-1]+2)
				params = np.polyfit(temp_v[block_temp],temp[block_temp],1) #
				new_params = [scaling[0]/scaling[1]*params[0], scaling[0]*params[1]]
				p = np.poly1d(new_params)
				new_block = range(leafindex(temp_v[block_temp[0]]*scaling[1]),leafindex(temp_v[block_temp[-1]]*scaling[1])+1)
			#print new_block, p(temp_v[new_block]),new_bds
			new_block_f = [i for i in new_block if p(temp_v[i])>new_bds[0] and p(temp_v[i])<new_bds[1]] # only update blocks within field
			new_block_f1, new_block_f2 = decompose(new_block_f)
			if len(new_block_f2)>len(new_block_f1):
				new_block_f = new_block_f2
			else:
				new_block_f = new_block_f1
			
			new_mlc[np.add(new_block_f,l*num_pp)] = np.round(p(temp_v[new_block_f]),2)
			#print "new_blocks:",temp_v[new_block_f], np.round(p(temp_v[new_block_f]),2)

	for mm in range(0, len(new_mlc)):
		new_mlc[mm] = str(round(new_mlc[mm],2))
	str_mlc_file = new_mlc.tolist()
	return str_mlc_file

def rtplangenerator(file_plan_ref, file_struct_ref, file_landmarks_ref, file_struct_sur,
	file_landmarks_sur,
	outfile_plan_sur):
	refID = file_struct_ref[file_struct_ref .rfind('RS.') + 3 : file_struct_ref.rfind('.dcm') ]
	surID = file_struct_sur[file_struct_sur .rfind('RS.') + 3 : file_struct_sur.rfind('.dcm') ]
	pr = PR()
	landmark_params_ref = pr.readFile(file_landmarks_ref)
	landmark_params_sur = pr.readFile(file_landmarks_sur)

	plan = pydicom.dcmread(file_plan_ref)
	struct_ref = pydicom.read_file(file_struct_ref)
	struct_sur  = pydicom.dcmread(file_struct_sur)
	# change plan name
	old_PlanName = plan.RTPlanName
	old_PlanName = old_PlanName.split('_')[0]
	plan.RTPlanName = 'auto_'+old_PlanName
	plan.RTPlanLabel= 'auto_'+old_PlanName

	# change studyID of the plan according to the surrogate patient
	plan.StudyID = struct_sur.StudyID
	plan.StudyInstanceUID = struct_sur.StudyInstanceUID
	plan.StudyTime = struct_sur.StudyTime

	# change patientID, name, Date-of_birth, sex to those surrogate patient (even it's anonymized ID/name)
	plan.PatientID = struct_sur.PatientID
	plan.PatientName = struct_sur.PatientName
	plan.PatientBirthDate = struct_sur.PatientBirthDate
	plan.PatientSex = struct_sur.PatientSex


	plan.FrameOfReferenceUID = struct_sur.ReferencedFrameOfReferenceSequence[0].FrameOfReferenceUID
	plan.ReferencedStructureSetSequence[0].ReferencedSOPInstanceUID = struct_sur.SOPInstanceUID
	if len(plan.DoseReferenceSequence) > 0:
		del plan.DoseReferenceSequence[0] # delete the original dose reference sequence

	# read some parameters to change the plan
	#isoc = np.array( plan.BeamSequence[0].ControlPointSequence[0].IsocenterPosition ).astype(float)
	field_borders = np.empty([2,2])

	nbeams = len(plan.BeamSequence)

	# find new isoc beforehand, to avoid repeating computations
	isoc = np.squeeze( np.array( plan.BeamSequence[0].ControlPointSequence[0].IsocenterPosition ) ) .astype(float)
	isoc_sur = np.array([0.0, 0.0, 0.0])
	tpm = TPM(landmark_params_ref, landmark_params_sur, struct_ref, struct_sur,isoc)
	for landmark in range(0,4):
		isoc_sur += tpm.getNewIsocenter( landmark )
		#print tpm.getNewIsocenter( landmark )

	isoc_sur = np.round(isoc_sur/4,2)
	tpm.isoc_sur = isoc_sur # update the isoc_sur in tpm

	# shift the isocenter UP
	for i in range(0,nbeams):
		#plan.BeamSequence[i][0x300b,0x100a].value = 'AMC-U3-v2'  
		cps = plan.BeamSequence[i].ControlPointSequence[0]
		# Remove pre-computed values
		cps.SourceToSurfaceDistance=""
		nseq = len(plan.BeamSequence[i].ControlPointSequence)
		# remove the referenced dose tag in all the control point sequence
		for j in range(0,nseq):
			cps_j = plan.BeamSequence[i].ControlPointSequence[j]
			if hasattr(cps, 'ReferencedDoseReferenceSequence'):
				del cps_j.ReferencedDoseReferenceSequence

		# Change isocenter position
		cps.IsocenterPosition = [isoc_sur[0],isoc_sur[1],isoc_sur[2]]
	        
		# Change collimator angle
		collim = cps.BeamLimitingDeviceAngle # the initial collimator angle
		#print collim
		diff_collim = round(tpm.getNewCollimatorAngle() * 180 / np.pi)
		rot_list = [diff_collim % 360, -diff_collim % 360]
		min_ang = min(rot_list) #  was attenuated by /2
		rot_sig = 1
		if diff_collim !=0:
			rot_sig = -2*rot_list.index(min(rot_list))+1
		
		if cps.GantryAngle==0:
			collim = (rot_sig*min_ang+collim) % 360
		else:
			collim = (-rot_sig*min_ang+collim) % 360  # take a modulo operation 
		#print 'new collim: ', collim
		cps.BeamLimitingDeviceAngle = collim
		
		#Change Plan Borders
		lds = cps.BeamLimitingDevicePositionSequence
		for k in range(0,len(lds)):
			if lds[k].RTBeamLimitingDeviceType =='ASYMX':
				old_xbds = np.squeeze(np.array(lds[k].LeafJawPositions)).astype(float)
				new_xbds = tpm.getNewFieldBorders( old_xbds, 0 ,collim,i)
				lds[k].LeafJawPositions = [new_xbds[0], new_xbds[1]]
				#print old_xbds,new_xbds
			elif lds[k].RTBeamLimitingDeviceType =='ASYMY':
				old_ybds = np.squeeze(np.array(lds[k].LeafJawPositions)).astype(float)
				new_ybds = tpm.getNewFieldBorders( old_ybds, 1 , collim,i)
				lds[k].LeafJawPositions = [new_ybds[0], new_ybds[1]]
				#print old_ybds,new_ybds
		# change Plan MLCs (dependent on new border's information)
		for k in range(0,len(lds)):
			if lds[k].RTBeamLimitingDeviceType =='MLCX':
				old_mlcx = np.array(np.array(lds[k].LeafJawPositions))
				#print "method old:", new_mlcx_file
				if collim > 345 or collim <15:
					new_mlcx_file = modelmlc(old_mlcx, old_xbds, new_xbds, old_ybds, new_ybds,tpm,0,collim,i)
				else:
					new_mlcx_file = reshapemlc(old_mlcx, old_xbds, new_xbds, old_ybds, new_ybds,tpm,0,collim,i)
				#print "method new:", new_mlcx_file	
				lds[k].LeafJawPositions = new_mlcx_file
				#print "...MLC in x updated", old_mlcx, str_mlcx_file
			elif lds[k].RTBeamLimitingDeviceType =='MLCY':
				print "MLCY!"
				old_mlcy = np.array(np.array(lds[k].LeafJawPositions))
				new_mlcy_file = modelmlc(old_mlcy, old_ybds, new_ybds,old_xbds, new_xbds,tpm,1,collim,i)
				lds[k].LeafJawPositions = str_mlcy
				#print "... MLC in y updated"
	print "plan emulated!"
	plan.save_as(outfile_plan_sur)
	return True

def featuremlc(old_mlc, old_bds, old_vbds):
	num_pp = int(round(len(old_mlc)/2))
	vs_half = np.array(10*(np.arange(40)-20)+5).astype(float)
	vs = np.append(vs_half,vs_half)
	for l in range(0,2):
		temp = old_mlc[np.add(range(0,num_pp),l*num_pp)]
		temp_v = vs_half
		idxblock = [i for i in range(0,len(temp)) if temp[i] >old_bds[0] and temp[i] < old_bds[1] and temp_v[i] > old_vbds[0]-5 and temp_v[i] < old_vbds[1]+5]
		block1, block2 = decompose(idxblock)
		block = block1
		if len(block2) > len(block1):
			block = block2                    # only consider the larger block (the upper one)
		if len(block) <2:
			continue
		elif len(block)> 1:
			if temp_v[block[0]-1]<old_vbds[0] and temp_v[block[-1]+1] < old_vbds[1]:
				block_temp = range(block[0], block[-1]+2)
			elif temp_v[block[0]-1]<old_vbds[0] and temp_v[block[-1]+1] >= old_vbds[1]:
				block_temp = range(block[0], block[-1]+1)
			elif temp_v[block[-1]+1]>old_vbds[1]:
				block_temp = range(block[0]-1,block[-1]+1)
			else:
				block_temp = range(block[0]-1,block[-1]+2)
			params= np.polyfit(temp_v[block_temp],temp[block_temp],1)
			p = np.poly1d(params)
			x_pos = p(temp_v)
			xy = [(x_pos[i],temp_v[i]) for i in range(0,len(x_pos))]
			mlc_path = LineString(xy)
			field_shape = box(old_bds[0], old_vbds[0],old_bds[1],old_vbds[1]) # a rectangular polygon
			# fig = plt.figure(1,figsize=(15,15),dpi=90)
			# ax = fig.add_subplot(111)
			# ax.plot(x_pos,temp_v,color='#6699cc', alpha=1, linewidth=3)
			# x, y =field_shape.exterior.coords.xy
			# ax.plot(x,y,color='red')
			# plt.savefig('temp.png')
			#pdb.set_trace()
			parts = split(field_shape,mlc_path)
			remain_area = 0
			for temp in parts:
				if temp.area > remain_area:
					remain_area = temp.area
			block_ratio = np.round(1 - remain_area/field_shape.area,2)
			#print block_ratio,np.round(abs(1/params[0]),2)
			return block_ratio, np.round(abs(1/params[0]),2)
	return 0,0
			#print "new_blocks:",temp_v[new_block_f], np.round(p(temp_v[new_block_f]),2)

def planfeatures(file_plan):
	plan = pydicom.dcmread(file_plan, force=True)
	isoc = np.squeeze( np.array( plan.BeamSequence[0].ControlPointSequence[0].IsocenterPosition ) ) .astype(float)
	nbeams = len(plan.BeamSequence)
	i = 0 # only consider AP
	cps = plan.BeamSequence[i].ControlPointSequence[0]
	collimangle = cps.BeamLimitingDeviceAngle 
	lds = cps.BeamLimitingDevicePositionSequence
	num_px =0
	num_py =0
	xbds = None
	for k in range(0,len(lds)):
		if lds[k].RTBeamLimitingDeviceType =='ASYMX':
			xbds = np.squeeze(np.array(lds[k].LeafJawPositions)).astype(float)	
		elif lds[k].RTBeamLimitingDeviceType =='ASYMY':
			ybds = np.squeeze(np.array(lds[k].LeafJawPositions)).astype(float)
		elif lds[k].RTBeamLimitingDeviceType =='MLCX':
			mlcx = np.array(np.array(lds[k].LeafJawPositions))
			num_px = int(round(len(mlcx)/4))
		else: 
			mlcy = np.array(np.array(lds[k].LeafJawPositions))
			num_py = int(round(len(mlcy)/4))
	
	file_APdiam = file_plan+".txt"
	f_newline = ''                                                                              
	f_newline = f_newline + "isocenter"+ '\t' + str(isoc) +'\n'
	f_newline = f_newline + "CollimatorAngle"+ '\t[' + str(collimangle) +']\n'
	f_newline = f_newline + "field_borders_x"+ '\t' + str(xbds) +'\n'
	f_newline = f_newline + "field_borders_z"+ '\t' + str(ybds) +'\n'
	if num_px >0:
		f_newline = f_newline + "field_MLC_x"+ '\t' + str(mlcx) +'\n'
		block_ratio, block_slope = featuremlc(mlcx, xbds, ybds)
		f_newline = f_newline + "block_ratio" + '\t' + str(block_ratio)+'\n'
		f_newline = f_newline + "block_slope" + '\t' + str(block_slope) + '\n'  
	if num_py >0:
		f_newline = f_newline + "field_MLC_z"+ '\t' + str(mlcy) +'\n'
		block_ratio, block_slope = featuremlc(mlcy, ybds, xbds)
		f_newline = f_newline + "block_ratio" + '\t' + str(block_ratio)+'\n'
		f_newline = f_newline + "block_slope" + '\t' + str(block_slope) + '\n' 

	return(f_newline)

def printpage():
    print "---------------Help info--------------"
    print "Plan emulation function! "
    print "parameters: file_plan_ref,file_struct_ref,file_landmarks_ref,file_struct_sur,file_landmarks_sur,outfile_plan_sur"
    quit()
    
def main():
	if len(sys.argv) == 1:
		printpage()
	else:
		if (len(sys.argv) == 2 and sys.argv[1] == '-help'):
			printpage()

	file_plan_ref = sys.argv[1]
	file_struct_ref = sys.argv[2]
	file_landmarks_ref = sys.argv[3]

	file_struct_sur = sys.argv[4]
	file_landmarks_sur = sys.argv[5]
	outfile_plan_sur = sys.argv[6]
	rtplangenerator(file_plan_ref,file_struct_ref,file_landmarks_ref,file_struct_sur,file_landmarks_sur,outfile_plan_sur)

if __name__ == "__main__":
    main()
