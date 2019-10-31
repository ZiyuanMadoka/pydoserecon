import os
import sys
import numpy as np
from IsocAPRatios import IsocAPRatiosFinder
import math

class TransformationParamMaker(object):

	def __init__(self, params_ref, params_sur, rt_struct_ref, rt_struct_sur,isoc_ref):
		self.params_ref = params_ref
		self.params_sur = params_sur
		self.rt_struct_ref = rt_struct_ref
		self.rt_struct_sur = rt_struct_sur
		self.isoc_ref = isoc_ref
		self.isoc_sur = None
		if isoc_ref[0] - params_ref["Th10_bottom_xyz"][0]>0:
			self.side = "left"
			self.land_set = ["L2_right_cor_xyz", "Rib_left_xyz","Th12_right_cor_xyz","L4_bottom_xyz", "Rib_right_xyz", "L2_left_cor_xyz"]
		else:
			self.side = "right"
			self.land_set= ["L2_left_cor_xyz","Rib_right_xyz","Th12_left_cor_xyz","L4_bottom_xyz","Rib_left_xyz","L2_right_cor_xyz"]
		#self.sizes_ref = np.array([params_ref["CT_size"][0], params_ref["CT_size"][2]]);
		#self.sizes_sur = np.array([params_sur["CT_size"][0], params_sur["CT_size"][2]]);
		#CT_origin is not in measurement txt files. select reference point as L2_right_cor_xy, can be changed later
	

	def rotateCoordinates(self, coords, theta):
	 # theta is contourclock-wise --updated Jan31 2019
		xr = coords[0] * np.cos (theta) - coords[1] * np.sin (theta)
		yr = coords[1] * np.cos (theta) + coords[0] * np.sin (theta)

		return np.array([xr,yr]);

	def getNewCollimatorAngle(self):
		# the collimator here is contour-clockwise
		collimator_angle_ref = self.params_ref['Collimator_angle'] * np.pi / 180
		collimator_angle_sur = self.params_sur['Collimator_angle'] * np.pi / 180
		diff_angle = collimator_angle_sur - collimator_angle_ref
		return diff_angle
	def getScalingSurToRef(self):
		if self.side == "left":
			x_scaling =  np.divide( self.params_sur['Rib_width_left'],  self.params_ref['Rib_width_left'] )
		else:
			x_scaling =  np.divide( self.params_sur['Rib_width_right'],  self.params_ref['Rib_width_right'] )
        # updated the scaling reference length as T11L4, can be changed later
		ref_length = "Length_T11L4"
		y_scaling = np.divide( self.params_sur[ref_length],  self.params_ref[ref_length] )
		return np.array( [x_scaling, y_scaling] );

	def getScalingVerteSurToRef(self):
		if self.side == "left": # the part over spine
			x_scaling =  np.divide( self.params_sur['Rib_width_right'] ,  self.params_ref['Rib_width_right'])
			#x_scaling =  np.divide( self.params_sur['Rib_width_left'] +self.params_sur['Verte_width']/2,  self.params_ref['Rib_width_left']+self.params_ref['Verte_width']/2)
		else:
			x_scaling =  np.divide( self.params_sur['Rib_width_left'],  self.params_ref['Rib_width_left'])
			#x_scaling =  np.divide( self.params_sur['Rib_width_right']+self.params_sur['Verte_width']/2,  self.params_ref['Rib_width_right']+self.params_ref['Verte_width']/2)
		#x_scaling =  np.divide( self.params_sur['Rib_width'],  self.params_ref['Rib_width'] )
		return x_scaling


	def getNewIsocenter(self, landmark):
		collimator_angle_ref = self.params_ref['Collimator_angle'] * np.pi / 180
		collimator_angle_sur = self.params_sur['Collimator_angle'] * np.pi / 180

		scaling = self.getScalingSurToRef()
		isoc_ref = self.isoc_ref
		# Find ratios in AP direction for the isocenter position
		is_ap_rf = IsocAPRatiosFinder(isoc_ref, self.rt_struct_ref) # initialize is_ap_rf with reference contour !!
		
		isoc_ap_ratios_ref = is_ap_rf.getIsocRatios()
		isoc_ap = isoc_ref[1]
		
		# Use only coords LR and CC for what follows
		isoc_ref_xz = np.array( [isoc_ref[0], isoc_ref[2]] )
		land_ref_para_set = self.land_set # get the set of landmarks for the reference CT
		
		land_ref_para = land_ref_para_set[landmark]
		land_ref_xz = np.array([self.params_ref[land_ref_para][0], self.params_ref[land_ref_para][2]])
		
		isoc_vector_CT =  isoc_ref_xz - land_ref_xz # the vector from landmark to isocenter in CT coordintate system 
		isoc_vector_sur = self.rotateCoordinates(isoc_vector_CT, -collimator_angle_ref + collimator_angle_sur) # rotates the vector from reference CT system to surrogate CT system  [CT rotates -theata relative to collimator system] 
		land_sur_xz = np.array([self.params_sur[land_ref_para][0], self.params_sur[land_ref_para][2]])

		# Delta of the surrogate, obtained by scaling the delfta of ref
		isoc_sur_delta = isoc_vector_sur * scaling
		#calculate isocenter position in surrogate CT system
		isoc_sur_xz = land_sur_xz + isoc_sur_delta
		# Compute AP isoc surrogate
		isoc_ap_sur = is_ap_rf.getIsocAPFromRatios(isoc_ap_ratios_ref, self.rt_struct_sur, isoc_sur_xz)
		# checking ap value usidng surrogate contour!!

		# Set isoc sur back to 3D, adding the AP
		isoc_sur = np.round(np.array([ isoc_sur_xz[0], isoc_ap_sur, isoc_sur_xz[1] ]),2) 
		
		return isoc_sur
	
	def getNewSpineSideBorder(self, b_ref,scaling,collim, gantry):
		flip = -2 * gantry +1
		if self.side =="left":
			landmark ="L2_right_cor_xyz"
		else:
			landmark = "L2_left_cor_xyz"
		
		dis_ref = (flip*self.isoc_ref[0]  -flip*self.params_ref[landmark][0])*math.cos(math.radians(collim)) + b_ref
		dis_sur = dis_ref * scaling
		b_sur = dis_sur + (flip*self.params_sur[landmark][0]-flip*self.isoc_sur[0])*math.cos(math.radians(collim))
		#print collim,math.cos(math.radians(collim)), scaling
		#print b_ref, b_sur
		return b_sur

	def getNewFieldBorders(self, fb_ref, coord,collim, gantry):
		if collim < 15 or collim >345:
			coord_id = coord
		else:
			coord_id = abs(coord-1) 
		#print coord, coord_id

		# field border is defined in collimator system, relative to the isocenter.
		scaling = self.getScalingSurToRef()
		#print scaling
		
		fb_sur = fb_ref * scaling[coord_id]
		scaling2 = self.getScalingVerteSurToRef(); # the scaling factor for border along the spine side
		if coord_id ==0 and len(fb_sur)==2: # for borders along RL
			if self.side =="left":
				spine_side = gantry # for gantry angle 0(AP), the first index, for gantry 180 (PA), the second index
			else:
				spine_side = abs(gantry-1)
			#print gantry, fb_ref,fb_sur
			fb_sur[spine_side] = self.getNewSpineSideBorder(fb_ref[spine_side], scaling2, collim,gantry) # the field border at the spine side is scaled indipendently
			#print "new: ", scaling[0],scaling2,fb_sur
		fb_sur = np.round(fb_sur,1) 
		return fb_sur

	def getAPdiameter(self,isoc):
		is_ap_rf = IsocAPRatiosFinder(isoc, self.rt_struct_sur)
		b = is_ap_rf.getBodyContourBoundsAtLRCCIsoc() 
		#ap_top: b[1], ap_bottom: b[3] the ap direction goes from anterior to posterior, so the top has lower values
		diameter = np.round(b[3]-b[1],1)
		return diameter
