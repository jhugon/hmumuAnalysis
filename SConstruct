
# SPODS--Simple Physics Object based Detector Simulation
# Copyright (C) 2011, 2012 Justin Hugon

# Main SCONS Build Configuration File

import os
import glob

env = Environment(ENV = {'PATH':os.environ['PATH']} )

includes = []
libpath = []
libs = []
env.MergeFlags('-fPIC -O2 -lm')

#root
env.ParseConfig('root-config --cflags')
env.ParseConfig('root-config --libs')
libs.append("EG")
libs.append("TMVA")
libs.append("Minuit2")
libs.append("Minuit")
libs.append("MathMore")
libs.append("XMLIO")
libs.append("MLP")
libs.append("TreePlayer")

env.Append(CPPPATH=includes)
env.Append(LIBPATH=libpath)
env.Append(LIBS=libs)

#Checking for things to exist (~autoConf)
if not env.GetOption("clean"):
  conf = Configure(env)
  if not conf.CheckCXX():
    print("Error: C++ Compiler Broken")
    Exit(1)
  if not conf.CheckSHCXX():
    print("Error: C++ Compiler Broken")
    Exit(1)
  if not conf.CheckLibWithHeader("Hist","TH1F.h","c++","TH1F h;"):
    print("Error: ROOT libs must be installed!")
    Exit(1)
  if not conf.CheckLibWithHeader("EG","TGenerator.h","c++","TGenerator g;"):
    print("Error: ROOT lib libEG.a must be installed!")
    Exit(1)
  if not conf.CheckLibWithHeader("TMVA",["TFile.h","TMVA/Factory.h"],"c++",'TMVA::Factory f("",&(TFile("")),"");'):
    print("Error: ROOT lib libTMVA.a must be installed!")
    Exit(1)
  env = conf.Finish()
 
env.Program(target="analyzer", source=["analyzer.cc"])
env.Program(target="trainingTreeMaker", source=["trainingTreeMaker.cc"])
