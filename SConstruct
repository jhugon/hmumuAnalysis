import os
import sys
import glob

env = Environment(ENV = {'PATH':os.environ['PATH']} )

includes = []
libpath = []
libs = []
env.MergeFlags('-fPIC -O2 -lm')

includes.append("src/")

#boost
#includes.append("/usr/include/boost141")
libs.append("boost_program_options")
libs.append("boost_regex")
#libpath.append("/usr/lib/boost141")

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
  if not conf.CheckLibWithHeader("TMVA",["TMVA/Tools.h"],"c++",'TMVA::Tools::Instance();'):
    print("Error: ROOT lib libTMVA.a must be installed!")
    Exit(1)
  if not conf.CheckCXXHeader("boost/program_options.hpp"):
    print("Error: boost/program_options.hpp header not installed!")
    Exit(1)
  if not conf.CheckLibWithHeader("boost_program_options","boost/program_options.hpp","c++",'boost::program_options::options_description optionDesc("options");'):
    print("Error: boost/program_options.hpp header and lib must be installed!")
  if not conf.CheckLibWithHeader("boost_regex","boost/regex.hpp","c++",'boost::regex re("aregex");'):
    print("Error: boost/regex.hpp header and lib must be installed!")
    Exit(1)
  env = conf.Finish()
 
env.Library(targer="src/mva",source=["src/mva.cc"])
env.Library(targer="src/helpers",source=["src/helpers.cc"])
env.Program(target="analyzer", source=["analyzer.cc","src/libmva.a","src/libhelpers.a"])
env.Program(target="mvaTrain", source=["mvaTrain.cc"])
