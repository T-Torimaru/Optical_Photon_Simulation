//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B5EventAction.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file B5EventAction.hh
/// \brief Definition of the B5EventAction class

#ifndef OpNoviceEventAction_h
#define OpNoviceEventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"
#include "OpNoviceScintHit.hh"
#include "OpNoviceTriggerHit.hh"

/// Event action

class OpNoviceRunAction;
class OpNoviceEventAction : public G4UserEventAction
{
public:
    OpNoviceEventAction();
    virtual ~OpNoviceEventAction();
    
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  void AddDep(G4double de);
  void AddDep_trigger(G4double de_tri);
    
private:
  // methods                                                                                    
  OpNoviceTriggerHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;

    G4int fmppcID;
    G4int fHHC1ID;
    G4int fHHC2ID;
    G4int fDHC1ID;
    G4int fDHC2ID;
    G4int fECHCID;
    G4int fHCHCID;
    G4int ftriggerHCID;
  G4int fAbsHCID;
  
  G4int fMPPCCollID;
  G4int fTRIGGERCollID;
  G4double fEnergyDep;
  G4double fEnergyDep_trigger;

};  
  inline void OpNoviceEventAction::AddDep(G4double de){
    fEnergyDep += de;
}
  inline void OpNoviceEventAction::AddDep_trigger(G4double de_tri){
    fEnergyDep_trigger += de_tri;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
