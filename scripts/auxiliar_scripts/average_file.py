#!/bin/python3

#This short script reads a colon of numbers from a file and returns their average


import numpy as np
import argparse as ap

#Get name of the file
parser=ap.ArgumentParser()
parser.add_argument("-name",type=str,help='Name of the data file as column of numbers')
args=parser.parse_args()
filename=args.name

#Load data
data=np.loadtxt(filename)

#Calculate its average
average=np.mean(data)

#Get standard dev
std_dev=np.std(data)

#Show average result
print(average,std_dev)


