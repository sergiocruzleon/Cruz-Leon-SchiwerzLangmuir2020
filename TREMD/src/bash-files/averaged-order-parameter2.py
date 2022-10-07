#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import os


#Load the data from the user file:



#path = '../LIO/dist-demux/'
#retunr an numpy array that contains the distances between the ion and the position in the RNA.
def load_data(filename):
	return  np.loadtxt(filename, comments=('#','@', '&'))




#This method returns an array with the order (bound or unbound calculated from the distnaces) 
def order_parameter(data1, data2, cut, unbound):
    counter =0
    index=0
    OP_av=0
    OP=0
    OP_out=[]
    # Exponential decay parameters
    lambd=1.05
    beta=5.	
    #counter
    counter=0
    for i in range(0, len(data1[:,1])):
        if data1[i,1] <= cut or data2[i,1] <= cut:
            OP=1.
        else:
			if ( data1[i,1] > unbound and data2[i,1] > unbound): 
				OP=0.
			else:
				OP=1/(1+np.exp(beta*(0.5*(data1[i,1]+data2[i,1])-lambd*cut)))
        #Average the OP
        if counter<2:#20):
            OP_av+=OP
	    counter+=1
        if(counter==2):#20):
            for n in range(0, 2):#20):
                OP_out.append([index, OP_av/2])
		index+=4 # Trajectories are written every 4ps
            counter=0
            OP_av=0
    return np.array(OP_out)
            
# Save the Order-parameter in a file
# @name: name of the file to store the data
# @data: is the array of the order parameter 
def save_OP(name, data):
	#print data
	np.savetxt(name_file, np.c_[data[:,0],data[:,1]])



###################################
## Main program
###################################
#To run this program is required the files, cut_on (nm), cut_off(nm), and output_name

# Open the file file for reading
f = open(sys.argv[1])
#Filename or Path?
filename=sys.argv[1]
filename2=sys.argv[2]
#path=sys.argv[1]

#parameters
c_on=float(sys.argv[3])
c_off=float(sys.argv[4])

# Load the data from the file
parameters=load_data(filename)
parameters2=load_data(filename2)

data=order_parameter(parameters, parameters2, c_on, c_off)
#Save the file
name_file=sys.argv[5]

save_OP(name_file, data)

