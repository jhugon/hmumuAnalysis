
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

env.Append(CPPPATH=includes)
env.Append(LIBPATH=libpath)
env.Append(LIBS=libs)
 
env.Program(target="analyzer", source=["analyzer.cc"])
