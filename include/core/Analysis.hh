// ============================================================================
// Analysis.hh
// Analysis manager configuration for NNBAR detector simulation
// ============================================================================

#ifndef Analysis_h
#define Analysis_h 1

// Select analysis manager based on Geant4 version
#include "G4Version.hh"

#if G4VERSION_NUMBER >= 1100
// Geant4 11.x uses unified analysis manager
#include "G4AnalysisManager.hh"
#else
// Geant4 10.x uses ROOT-specific analysis manager
#include "G4RootAnalysisManager.hh"
using G4AnalysisManager = G4RootAnalysisManager;
#endif

#endif // Analysis_h
