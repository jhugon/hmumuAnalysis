
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
 
env.Program(target="analyzer", source=["analyzer.cc"])
env.Program(target="trainingTreeMaker", source=["trainingTreeMaker.cc"])
