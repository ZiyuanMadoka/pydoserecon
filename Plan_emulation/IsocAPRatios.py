### Reads RT plan and RT struct, determines AP distance between isoc and top of body, and  AP distance between isoc and bottom of body'''
import numpy as np
import pydicom
import sys
from shapely.geometry import Point, LineString, Polygon, LinearRing



class IsocAPRatiosFinder(object):

	def __init__(self, isoc=None, rt_struct=None):
		self.isoc = isoc
		self.rt_struct = rt_struct  # whicht rt_struct?? 

	def group(self, lst, n):
	    return zip(*[lst[i::n] for i in range(n)]) 

	def getContour(self, contour_name, rt_struct) :
		#contour_name can be Body or other OAR
		rt_contour_idx = -1
		found_contours = []
		for i in range(0, len(rt_struct.StructureSetROISequence)):
			cname = rt_struct.StructureSetROISequence[i].ROIName
			found_contours.append(cname)
			cname = cname.replace('_',' ')
			if cname == contour_name:
				rt_contour_idx = i
				break

		if rt_contour_idx == -1:
			print contour_name, "not found in RT DICOM. "
		#	for cname in found_contours:
		#		print cname
			return None # if not found, return nan

		# sort along z, increasing
		contour_slices = sorted(rt_struct.ROIContourSequence[rt_contour_idx].ContourSequence, key=lambda x: float( x.ContourData[-1]))
		# in sorted should not be self.rt_struct, but should be the rt_struct passed to this function !!!
		# convert to np format 
		contour_slices_np = []

		for cs in contour_slices:
			r = list( self.group(cs.ContourData, 3) )
			for i in range(0,len(r)):
				r[i] = [float(j) for j in r[i]]
			contour_slices_np.append( r )


		contour_slices_np = np.squeeze( np.array(contour_slices_np) )

		return contour_slices_np


	def getClosestContourSliceCC(self, cc, contour):
		best = None
		best_dist = 999999
		""" give an original cz, len, to check if there are multiple closed contours at the cz value
		"""
		slic_test = Polygon (list (tuple (map (tuple, contour[0] ))))
		area_p = slic_test.area    # an initialization of the previous contour area and previous cz

		cz_p = contour[0][0][2]
		for slic in contour:
			cz = slic[0][2] # z coordinate (CC)
			d = np.abs(cc - cz)
			if d < best_dist:
				best_dist = d
				best = slic
				cz_p = cz
				slic_test = Polygon (list (tuple (map (tuple, slic ))))
				area_p = slic_test.area
			elif d == best_dist:  # check if one cc corresponds to multiple closed contours, find the body
				slic_test = Polygon (list (tuple (map (tuple, slic ))))
				if slic_test.area>area_p:   # assume closed contour with the largest area is the body
					best = slic
					cz_p = cz
					area_p = slic_test.area
			else :
				break	# because they are sorted according to z
		return best
    


	#def getBodyContourBoundsAtLRCCIsoc(self, rt_struct=self.rt_struct, isoc=self.isoc): give an error, self is not defined:
	def getBodyContourBoundsAtLRCCIsoc(self, rt_struct=None, isoc=None):
		"""
		compute the body contour bounds along AP (anteriro-posterior) direction at the isocenter plane z (CC) - x (LR) 
		"""
		if rt_struct is None:
			rt_struct = self.rt_struct
		
		if isoc is None:
			isoc = self.isoc
		
		body_contour = self.getContour('Body', rt_struct )
		if body_contour is None:
			exit()
		# get the slice points at given isocenter plane (CC, isoc[2])
		slic_points = self.getClosestContourSliceCC( isoc[2], body_contour )
		slic = Polygon (list (tuple (map (tuple, slic_points ))))
		# get the bounds of this slice (RL & AP)
		b = slic.bounds
		line_ap = LineString([ (isoc[0], b[1]), (isoc[0], b[3]) ]) # line between min AP and max AP, at fixed LR
		#intersection_line = LineString (list(slic.intersection(line_ap).coords) )
		ring = LinearRing(list(slic.exterior.coords)) # convert the polygon to linearRing: see gis.stackexchange.com/questions/128154/multi-part-geometrical-intersection-error
		intersection_line = LineString(list(line_ap.intersection(ring))) # intersection of line_ap on the slice and the contour on the slice
		b = intersection_line.bounds
		return b


	#def getAPPABodyBoundsBodyContourAtIsoc(rt_struct, isoc):
	#	b = self.getBodyContourBoundsAtLRCCIsoc(rt_struct, isoc) # is a member function of self object
		
	#	return b[1], b[3]

	def getIsocRatios(self):
		b = self.getBodyContourBoundsAtLRCCIsoc() # is a member function of self object
		ratio_bottom = (b[3] - self.isoc[1])/(b[3] - b[1])
		ratio_top = (self.isoc[1] - b[1])/(b[3] - b[1])
		return ratio_top, ratio_bottom


	# given the AP ratios, the struct file, and the LR and CC coords of the isoc, gives the AP coord of the isoc
	def getIsocAPFromRatios(self, ap_ratios, rt_struct, lr_cc_isoc):
		isoc_w_AP_placeholder = [  lr_cc_isoc[0], 0.0,  lr_cc_isoc[1] ] 
		b = self.getBodyContourBoundsAtLRCCIsoc(rt_struct, isoc_w_AP_placeholder )
        #ap_top: b[1], ap_bottom: b[3] the ap direction goes from anterior to posterior, so the top has lower values
		isoc_ap = b[1] + (b[3] - b[1]) * ap_ratios[0]
		#isoc_ap2 = b[3] - (b[3] - b[1]) * ap_ratios[1]  # it's exactly the same. since ap_ratios[0]+ap_ratios[1] = 1
	    
		return isoc_ap
