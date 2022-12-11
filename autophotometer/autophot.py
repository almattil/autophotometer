#! /usr/bin/env python3
LocalVersion='20221211'

import os
import sys
from sys import stdin
from pkg_resources import get_distribution

try:
    import autophotometer.auto_main as auto_main
    isModuleInstalled=1
except ModuleNotFoundError:
    isModuleInstalled=0
    import auto_main


#######################################################################################
def printHeader(version):
    """ Print the Header """
    print("\n")
    print("*******************************************************************************")
    print("                         autophotometer: ver.",version)
    print(" ")

def printUsage():
    """ Print the usage info """
    print("Usage:")
    print("autophot.py ObsFitsImage1 [ObsFitsImage2]\n")
    print("Parameters:")
    print("  ObsFitsImage1: FITS-image.")
    print("  ObsFitsImage2: optional second FITS-image.")
    print("\n*******************************************************************************\n")
#    print(" ")

##########################################################################


def cli_main():
	if isModuleInstalled:
		version = get_distribution('autophotometer').version
	else:
		version = LocalVersion+' local\n'
	printHeader(version)
	fits_file2 = None
	if (len(sys.argv) < 2):
		printUsage()
	if (len(sys.argv) > 2):
		fits_file2 = sys.argv[2]
	if (len(sys.argv) > 1):
		fits_file1 = sys.argv[1]
	
	else:
		print(" ")
		sys.stdout.write("Enter the filename of the FITS-image: ")
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


if __name__ == "__main__":
	cli_main()
