#! /usr/bin/env python3
# ver. 20210906


import os
import sys
from sys import stdin
import starseparator
import auto_main

#######################################################################################
def printHeader():
    """ Print the Header """
    print(" ")
    print("*******************************************************************************")
    print("                         autophot: ver. 20210906")
    print(" ")

def printUsage():
    """ Print the usage info """
    print("Usage:")
    print("autophot.py ObsFitsImage1")
    print(" ")
    print("Parameters:")
    print("  ObsFitsImage1: FITs-image.")
    print(" ")
    print("*******************************************************************************")
    print(" ")

##########################################################################


if __name__ == "__main__":
	printHeader()
	fits_file2 = None
	if (len(sys.argv) < 2):
		printUsage()
	if (len(sys.argv) > 2):
		fits_file2 = sys.argv[2]
	if (len(sys.argv) > 1):
		fits_file1 = sys.argv[1]
	
	else:
		print(" ")
		sys.stdout.write("Enter the filename of the FITs-image: ")
		sys.stdout.flush()
		output_file_path = stdin.readline()
		fits_file1 = output_file_path.rstrip('\n')
	if not os.path.exists(fits_file1):
		print("The file", fits_file1, "doesn't exist. Stop executing!")
		exit()
	if fits_file2 != None and not os.path.exists(fits_file2):
		print("The file", fits_file2, "doesn't exist. Stop executing!")
		exit()

	if len(sys.argv) > 2:
		auto_main.MainProject_dual(fits_file1, fits_file2)
	elif len(sys.argv) > 1:
		auto_main.MainProject(fits_file1)
	else:
		print('Something went wrong')


