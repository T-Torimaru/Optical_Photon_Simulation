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
// $Id: B5HodoscopeHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file B5HodoscopeHit.hh
/// \brief Definition of the B5HodoscopeHit class

#ifndef OpNoviceTriggerHit_h
#define OpNoviceTriggerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "tls.hh"

class OpNoviceTriggerHit : public G4VHit
{
public:
  OpNoviceTriggerHit();
  OpNoviceTriggerHit(const OpNoviceTriggerHit&);
  virtual ~OpNoviceTriggerHit();

  // operators                                                                                                    
  const OpNoviceTriggerHit& operator=(const OpNoviceTriggerHit&);
  G4int operator==(const OpNoviceTriggerHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  // methods from base class                                                                                      
  virtual void Draw() {}
  virtual void Print();

  // methods to handle data                                                                                       
  void Add(G4double de, G4double dl);

  // get methods                                                                                                  
  G4double GetEdep() const;
  G4double GetTrackLength() const;

private:
  G4double fEdep;        ///< Energy deposit in the sensitive volume                                              
  G4double fTrackLength; ///< Track length in the  sensitive volume                                               
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                    

typedef G4THitsCollection<OpNoviceTriggerHit> OpNoviceTriggerHitsCollection;

extern G4ThreadLocal G4Allocator<OpNoviceTriggerHit>* OpNoviceTriggerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                    

inline void* OpNoviceTriggerHit::operator new(size_t)
{
  if(!OpNoviceTriggerHitAllocator)
    OpNoviceTriggerHitAllocator = new G4Allocator<OpNoviceTriggerHit>;
  void *hit;
  hit = (void *) OpNoviceTriggerHitAllocator->MallocSingle();

  return hit;
}

inline void OpNoviceTriggerHit::operator delete(void *hit)
{
  if(!OpNoviceTriggerHitAllocator)
    OpNoviceTriggerHitAllocator = new G4Allocator<OpNoviceTriggerHit>;
  OpNoviceTriggerHitAllocator->FreeSingle((OpNoviceTriggerHit*) hit);
}

inline void OpNoviceTriggerHit::Add(G4double de, G4double dl) {
  fEdep += de;
  fTrackLength += dl;
}

inline G4double OpNoviceTriggerHit::GetEdep() const {
  return fEdep;
}

inline G4double OpNoviceTriggerHit::GetTrackLength() const {
  return fTrackLength;
}

#endif
