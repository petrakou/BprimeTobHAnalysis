#!/usr/bin/env python

import sys, os, re, subprocess
from optparse import OptionParser


def make_filelist(input_dir):

    proc = subprocess.Popen( [ '/afs/cern.ch/project/eos/installation/cms/bin/eos.select', 'ls -l', input_dir ], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
    output = proc.communicate()[0]
    if proc.returncode != 0:
        print output
        sys.exit(1)

    return output.splitlines()


def process_input_dir(input_dir, sizelist):

    input_dir = input_dir.rstrip('/')+'/'
    filelist = make_filelist(input_dir)

    for line in filelist:
        sizelist.append(float(line.split()[4]))

    return


# usage description
usage = """Usage: ./getDatasetSize.py [options]\n
Example: ./getDatasetSize.py -d datasetList.txt\n
For more help: ./getDatasetSize.py --help
"""

def main():
  # input parameters
  parser = OptionParser(usage=usage)

  parser.add_option("-d", "--dataset_list", dest="dataset_list", action='store', help="Text file containing a list of datasets to be processed", metavar="DATASET_LIST")

  (options, args) = parser.parse_args()

  # make sure all necessary input parameters are provided
  if not options.dataset_list:
    print usage
    sys.exit()

  dataset_list = options.dataset_list

  # open and read the dataset_list file
  dataset_list_file = open(dataset_list,"r")
  dataset_list_lines = dataset_list_file.readlines()

  # loop over datasets
  for line in dataset_list_lines:
    line_elements = line.split()
    if (len(line_elements)==0 or line_elements[0][0]=='#'): continue

    print 'Dataset ' + line_elements[0]

    sizelist = []
    process_input_dir(line_elements[2], sizelist)

    size = 0.
    for s in sizelist:
      size = size + s

    print "  Size: %.1f kB, %.1f MB, %.1f GB"%(size/1024.,size/(1024.**2),size/(1024.**3))


  # close all open files
  dataset_list_file.close()


if __name__ == "__main__":
  main()
