#!/usr/bin/python

import sys
import os
import re
import argparse
import boto

bucketName = "cloud-hep-testing-1"
folderNames = [
"WHmumu",
"ZHmumu",
"vbfHmumu3",
"Zmumujets",
"ZmumujetsMgt100",
"ggHmumu"
]

bucketName = args.bucket
s3 = boto.connect_s3()
try:
  bucket = s3.get_bucket(bucketName)
except:
  print("Error: Bucket "+bucketName+"not found, exiting.")
  sys.exit(1)

for dirname in folderNames:
  while dirname[0]=="/":
    dirname = dirname[1:]
  
  needToMakeDir = False
  
  if not os.path.exists(dirname):
    needToMakeDir = True
  elif not os.path.isdir(dirname):
    print("Error: local path: "+dirname+" is not a directory")
  
  keyFound = False
  for key in bucket.get_all_keys():
    if re.match(r"^"+dirname+"/.*",key.name):
      keyFound = True
      break
  
  if not keyFound:
    print("Error: S3 Directory Not Found")
    sys.exit(1)
  
  if needToMakeDir:
    os.makedirs(dirname)
  
  for key in bucket.get_all_keys():
    if re.match(r"^"+dirname+"/.*",key.name):
      print("Getting file: "+key.name)
      key.get_contents_to_filename(key.name)
