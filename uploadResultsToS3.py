#!/usr/bin/python

# Usage: ./uploadResultsToS3.py outbucket outname outfiles...
print("Starting uploadResultsToS3.py")
import sys
import boto

s3 = boto.connect_s3()
bucket = s3.get_bucket(sys.argv[1])
outDir = sys.argv[2] + "/"
for i in sys.argv[3:]:
  key = bucket.new_key(outDir+i+".root")
  key.set_contents_from_filename(i+".root")
  key.set_acl("public-read")
  key.change_storage_class("REDUCED_REDUNDANCY")
  print("Created file {0} on bucket {1}".format(outDir+i+".root",sys.argv[1]))
  
print("Done uploadResultsToS3.py")
