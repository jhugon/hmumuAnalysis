
import sys
import ROOT as root

sys.argv
if len(sys.argv) != 3:
  print("Error: Must have command line arguments:")
  print("getEventNums.py <file1.root> <file2.root>")

print("File1: "+sys.argv[1])
print("File2: "+sys.argv[2])

f1 = root.TFile(sys.argv[1])
f2 = root.TFile(sys.argv[2])

t1 = f1.Get("outtree")
t2 = f2.Get("outtree")

print("nEvents1: "+str(t1.GetEntries()))
print("nEvents2: "+str(t2.GetEntries()))

s1 = set()
s2 = set()
m1 = dict()
m2 = dict()

for t,s,m in zip([t1,t2],[s1,s2],[m1,m2]):
  nEntries = t.GetEntries()
  for i in range(nEntries):
    t.GetEntry(i)
    tmp = str(t.eventInfo_run)
    tmp += ":"
    tmp += str(t.eventInfo_event)
    s.add(tmp)
    m[tmp] = i
    if i > 10000:
      break

missingIn1 = s2.difference(s1)
missingIn2 = s1.difference(s2)
print missingIn1
print missingIn2
missingIndex1 = [m2[i] for i in missingIn1]
missingIndex2 = [m1[i] for i in missingIn2]
print missingIndex1
print missingIndex2

first = True
for mi,t in zip([missingIndex2,missingIndex1],[t1,t2]):
  nEntries = t.GetEntries()
  if first:
    print("Missing Events in Tree 2 Present in Tree 1:")
    first = False
  else:
    print("Missing Events in Tree 1 Present in Tree 2:")
  for i in mi:
    t.GetEntry(i)
    #print(str(i)+" "+str(t.eventInfo_run)+":"+str(t.eventInfo_event))
    print "Event: %i:%i" % (t.eventInfo_run,t.eventInfo_event)
    for b in t.GetListOfBranches():
      bName = b.GetName()
      if "dimuon" not in bName and bName != "muonSub_pt" and bName != "muonLead_pt":
      #if "dijet" not in bName:
        continue
      print "  %25s %10.4g" % (bName,getattr(t,bName))

print("Missing Events in Tree 1 Present in Tree 2:")
for i in missingIn1:
  print i
print("Missing Events in Tree 2 Present in Tree 1:")
for i in missingIn2:
  print i

