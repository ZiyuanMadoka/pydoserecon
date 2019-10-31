import os
import sys
import numpy as np

## changed the class from ParameterReader to ParameterHandler

class ParameterHandler(object):
	
	""" Reads in the information from 2D X-ray (or DRR) """
	def readFile(self, fname):
		params = {}
		with open(fname) as f:
		    content = f.readlines()
		# you may also want to remove whitespace characters like `\n` at the end of each line
		content = [x.strip() for x in content] 
		for line in content:
			if len(line) < 1 or line.startswith('%'):
				continue
			else:
				x = line.split(' ')
				v = []
				for i in range(1,len(x)):
					if (len(x[i]) > 0):
						v.append( float(x[i]) )
				params[x[0]] = np.squeeze( np.array(v) )  # 0-D array could be directly used as a scalar value, 
 				# can be indexed as a.item() or a[()]
		return params
        """ writes in the information to a txt file """
	def writeFile(self, params,fname):
		file_ob = open(fname,'w') 
	       	for paraname in params:
		    line = paraname + ' ' 
		    paravalue = params[paraname]
		    if paravalue.ndim==0:
		    	word = paravalue.item()
		    	line+=str(word)
		    else:
		    	for word in params[paraname]:
		    		line+=str(word)
		    		line+=' '
		    line+= '\n'
		    file_ob.write(line)
		file_ob.close
		return (1)
