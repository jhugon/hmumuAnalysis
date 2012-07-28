#ifndef SPODS_PhysicsObjectFunctions
#define SPODS_PhysicsObjectFunctions

// SPODS--Simple Physics Object based Detector Simulation
// Copyright (C) 2011, 2012 Justin Hugon

#include <TParticle.h>
#include <iostream>
#include <vector>
#include <algorithm>

inline
void setGenScalingFactor(TParticle* particle, double scaleFactor)
	{particle->SetCalcMass(scaleFactor);};

inline
double getGenScalingFactor(TParticle* particle)
	{return particle->GetCalcMass();};

inline
void setIso(TParticle* particle, double iso)
{
	particle->SetWeight(iso);
};

inline
double getIso(TParticle* particle)
{
	return particle->GetWeight();
};

inline
void setGenB(TParticle* particle, bool genB)
{
	if (genB)
	  particle->SetFirstDaughter(1);
	else
	  particle->SetFirstDaughter(0);
};

inline
bool getGenB(TParticle* particle)
{
        int i = particle->GetFirstDaughter();
	if(i==1)
		return true;
	else if(i==0)
		return false;
	else
		std::cout << "Warning: genGenB got unknown value!!\n";
		return false;
};

inline
void setBTag(TParticle* particle, bool bTag)
{
	if (bTag)
	  particle->SetLastDaughter(1);
	else
	  particle->SetLastDaughter(0);
};

inline
bool getBTag(TParticle* particle)
{
        int i = particle->GetLastDaughter();
	if(i==1)
		return true;
	else if(i==0)
		return false;
	else
		std::cout << "Warning: genBTag got unknown value!!\n";
		return false;
};
#endif
