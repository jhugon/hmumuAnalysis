#!/usr/bin/env python

import os
import sys
import re

if len(sys.argv)<=2:
  print "Error: At least 2 filename inputs required, exiting"
  print "Usage ./printSystematic.py nominalFile systematicFile1 [systematicFile2] ..."
  print "where files are Anna's text efficiency files"
  print "The maximum relative difference from nominal in the systematic files will be the systematic output"
  sys.exit(1)

effFileName = sys.argv[1]
sysFileNames = sys.argv[2:]

print "Using %s as nominal efficiency" % effFileName
print "Using %s for systematics" % sysFileNames

effFile = open(effFileName, 'r')
sysFiles = [open(i, 'r') for i in sysFileNames]

matchString = r"([\w]+)[\s]+([-0-9.eE]+)[\s]+([-0-9.eE]+)[\s]*"

while True:
  try:
    nomLine = effFile.next()
    sysLines = [i.next() for i in sysFiles]
    if not "Jet" in nomLine: # to only do Baseline++ cats
      continue
    nomMatch = re.match(matchString,nomLine)
    if not nomMatch:
      print("Error: Nom line didn't match: %s" % nomLine)
      continue
    nomName = nomMatch.group(1)
    nomEff = float(nomMatch.group(2))
    statErr = float(nomMatch.group(3))
    relStatErr = statErr/nomEff
    errEffs = []
    correlation = None
    for i in sysLines:
      errMatch = re.match(matchString,i)
      if not errMatch:
        print("Error: Err line didn't match: %s" % i)
        errEffs.append(0.)
        continue
      tmpErrName = errMatch.group(1)
      tmpErrEff = float(errMatch.group(2))
      if tmpErrName != tmpErrName:
        print("Error: Err cat name didn't match nom: %s != %s" % (tmpErrName, nomName))
        errEff.append(0.)
        continue
      errEffs.append(abs(tmpErrEff-nomEff))
      if correlation == None:
        if tmpErrEff>=nomEff:
          correlation = ""
        else:
          correlation = "-"
        
    relErr = max(errEffs)/nomEff
    print "%s  %s%.4f" % (nomName, correlation,relErr)
    #print "%s  %s%.4f  %.4f" % (nomName, correlation, relErr,relStatErr)
  except StopIteration:
    break

# cleanup
for i in sysFiles:
  i.close()
effFile.close()
