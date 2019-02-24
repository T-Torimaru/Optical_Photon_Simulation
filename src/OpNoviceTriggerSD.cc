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
// $Id: B5HodoscopeSD.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file B5HodoscopeSD.cc
/// \brief Implementation of the B5HodoscopeSD class

#include "OpNoviceTriggerSD.hh"
#include "OpNoviceTriggerHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                            

OpNoviceTriggerSD::OpNoviceTriggerSD(
				   const G4String& name,
				   const G4String& hitsCollectionName,
				   G4int nofCells)
  : G4VSensitiveDetector(name),
    fHitsCollection(0),
    fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                            

OpNoviceTriggerSD::~OpNoviceTriggerSD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                            

void OpNoviceTriggerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection                                                                                                 
  fHitsCollection
    = new OpNoviceTriggerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce                                                                                             
  G4int hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );

  // Create hits                                                                                                            
  // fNofCells for cells + one more for total sums                                                                          
  for (G4int i=0; i<fNofCells+1; i++ ) {
    fHitsCollection->insert(new OpNoviceTriggerHit());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                            

G4bool OpNoviceTriggerSD::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{
  // energy deposit                                                                                                         
  G4double edep = step->GetTotalEnergyDeposit();

  // step length                                                                                                            
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }

  if ( edep==0. && stepLength == 0. ) return false;

  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());

  // Get calorimeter cell id                                                                                                
  G4int layerNumber = touchable->GetReplicaNumber(1);

  // Get hit accounting data for this cell                                                               
  OpNoviceTriggerHit* hit = (*fHitsCollection)[layerNumber];
  if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << layerNumber;
    G4Exception("OpNoviceTriggerSD::ProcessHits()",
		"MyCode0004", FatalException, msg);
  }

  // Get hit for total accounting                                                                                           
  OpNoviceTriggerHit* hitTotal
    = (*fHitsCollection)[fHitsCollection->entries()-1];

  // Add values                                                                                                             
  hit->Add(edep, stepLength);
  hitTotal->Add(edep, stepLength);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                            

void OpNoviceTriggerSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) {
    G4int nofHits = fHitsCollection->entries();
    G4cout << "\n-------->Hits Collection: in this event they are " << nofHits
	   << " hits in the tracker chambers: " << G4endl;
    for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}                   
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// OpNoviceTriggerSD::OpNoviceTriggerSD(G4String name)
// : G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
// {
//   //    G4String HCname = "hodoscopeColl";
//     G4String HCname = "triggerColl";
//     collectionName.insert(HCname);
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// OpNoviceTriggerSD::~OpNoviceTriggerSD()
// {}

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void OpNoviceTriggerSD::Initialize(G4HCofThisEvent* hce)
// {
//   //    fHitsCollection = new B5HodoscopeHitsCollection
//     fHitsCollection = new OpNoviceTriggerHitsCollection
//     (SensitiveDetectorName,collectionName[0]);
//     if (fHCID<0)
//     { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
//     hce->AddHitsCollection(fHCID,fHitsCollection);
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4bool OpNoviceTriggerSD::ProcessHits(G4Step* step, G4TouchableHistory*)
// {
//   if (step == NULL) return false;
//   // G4Track* theTrack = step-> GetTrack();
  
//   // if(theTrack->GetDefinition()
//   //    != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
  
//   G4double edep = step->GetTotalEnergyDeposit();
//   if (edep==0.) return true;
  
//   // G4StepPoint* thePostPoint = step->GetPostStepPoint();
  
//   // G4TouchableHistory* touchable
//   //   = (G4TouchableHistory*)(thePostPoint->GetTouchable());
//   // G4int copyNo = touchable->GetVolume()->GetCopyNo();
//   // G4double hitTime = thePostPoint->GetGlobalTime();
  
//   // // check if this finger already has a hit
//   // G4int ix = -1;
//   // for (G4int i=0;i<fHitsCollection->entries();i++)
//   //   {
//   //     if ((*fHitsCollection)[i]->GetID()==copyNo)
//   //       {
//   // 	  ix = i;
//   // 	  break;
//   //       }
//   //   }
  
//   // if (ix>=0)
//   //   // if it has, then take the earlier time
//   //   {
//   //     if ((*fHitsCollection)[ix]->GetTime()>hitTime)
//   //       { (*fHitsCollection)[ix]->SetTime(hitTime); }
//   //   }
//   // //else
//   //   // if not, create a new hit and set it to the collection
//   //    //{
//   //     OpNoviceTriggerHit* hit = new OpNoviceTriggerHit(copyNo,hitTime);
//   //     G4VPhysicalVolume* physical = touchable->GetVolume();
//   //     hit->SetLogV(physical->GetLogicalVolume());
//   //     G4AffineTransform transform 
//   // 	= touchable->GetHistory()->GetTopTransform();
//   //     transform.Invert();
//   //     hit->SetRot(transform.NetRotation());
//   //     hit->SetPos(transform.NetTranslation());
//   //     fHitsCollection->insert(hit);
//   //     // }    
      
//   G4StepPoint* preStepPoint = step->GetPreStepPoint();
  
//   G4TouchableHistory* touchable
//     = (G4TouchableHistory*)(preStepPoint->GetTouchable());
//   G4int copyNo = touchable->GetVolume()->GetCopyNo();
//   G4double hitTime = preStepPoint->GetGlobalTime();
  
//   // check if this finger already has a hit
//   // G4int ix = -1;
//   // for (G4int i=0;i<fHitsCollection->entries();i++)
//   //   {
//   //     if ((*fHitsCollection)[i]->GetID()==copyNo)
//   //       {
//   // 	  ix = i;
//   // 	  break;
//   //       }
//   //   }
  
//   // if (ix>=0)
//   //   // if it has, then take the earlier time
//   //   {
//   //     if ((*fHitsCollection)[ix]->GetTime()>hitTime)
//   //       { (*fHitsCollection)[ix]->SetTime(hitTime); }
//   //   }
//   // else
//     // if not, create a new hit and set it to the collection
//     //  {
//       OpNoviceTriggerHit* hit = new OpNoviceTriggerHit(copyNo,hitTime);
//       G4VPhysicalVolume* physical = touchable->GetVolume();
//       hit->SetLogV(physical->GetLogicalVolume());
//       G4AffineTransform transform 
//   	= touchable->GetHistory()->GetTopTransform();
//       transform.Invert();
//       hit->Add(edep);
//       hit->SetRot(transform.NetRotation());
//       hit->SetPos(transform.NetTranslation());
//       fHitsCollection->insert(hit);
//       //    }    
//   return true;
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
