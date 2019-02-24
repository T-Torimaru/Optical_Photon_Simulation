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
// $Id: OpNoviceSteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file OpNoviceSteppingAction.cc
/// \brief Implementation of the OpNoviceSteppingAction class

#include "OpNoviceSteppingAction.hh"
#include "OpNoviceEventAction.hh"
#include "OpNoviceRunAction.hh"
#include "OpNoviceDetectorConstruction.hh"

#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ThreeVector.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4Run.hh"
#include "G4StepStatus.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "OpNoviceAnalysis.hh"
#include <sstream>
//#include "G4LogicalVolume.hh"
//#include "G4SystemOfUnits.hh"

#include "G4OpticalPhoton.hh"
#include "G4UserLimits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceSteppingAction::OpNoviceSteppingAction(
					       const OpNoviceDetectorConstruction* detectorConstruction,
					       OpNoviceEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
    //,fScoringVolume(0)
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
  //  prevID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceSteppingAction::~OpNoviceSteppingAction()
{ ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceSteppingAction::UserSteppingAction(const G4Step* step)
{
  
  G4int eventNumber = G4RunManager::GetRunManager()->
    GetCurrentEvent()->GetEventID();
  
  G4StepPoint* PrePoint = step->GetPreStepPoint();
  G4StepPoint* PostPoint = step->GetPostStepPoint();
  
  G4Track* track = step->GetTrack();
  G4int trackID = track->GetTrackID();
  G4ThreeVector pos_pre = PrePoint->GetPosition();
  G4double x_pre = pos_pre.x();
  G4double y_pre = pos_pre.y();
  G4double z_pre = pos_pre.z();
  G4ThreeVector pos_post = PostPoint->GetPosition();
  G4double x_post = pos_post.x();
  G4double y_post = pos_post.y();
  G4double z_post = pos_post.z();
  G4ThreeVector mom_pre = PrePoint->GetMomentum();
  G4double mom_preX = mom_pre.x();
  G4double mom_preY = mom_pre.y();
  G4double mom_preZ = mom_pre.z();
  G4ThreeVector mom_post = PostPoint->GetMomentum();
  G4double mom_postX = mom_post.x();
  G4double mom_postY = mom_post.y();
  G4double mom_postZ = mom_post.z();
  
  G4double depth;
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // get volume of the current step
  G4VPhysicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  G4double edep = step->GetTotalEnergyDeposit();
  
  if ( volume == fDetConstruction->GetScintPV() ){
    fEventAction->AddDep(edep);
  }
  if ( volume == fDetConstruction->GetTriggerPV() ){
    fEventAction->AddDep_trigger(edep);
  }
  // if ( step->GetTrack()->GetDefinition()->GetPDGCharge()!=0.){
  // // G4cout << "test_pos : " << -z_pre << G4endl;
  // // G4cout << "energy : " << edep << G4endl;
  //   analysisManager->FillH2(0,z_pre,edep);
  
  // if ( -z_pre>-1.5 && -z_pre<-1.2 ){
  //   depth = -1.35;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=-1.2 && -z_pre<-0.9 ){
  //   depth = -1.05;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=-0.9 && -z_pre<-0.6 ){
  //   depth = -0.75;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=-0.6 && -z_pre<-0.3 ){
  //   depth = -0.45;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=-0.3 && -z_pre<0. ){
  //   depth = -0.15;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=0. && -z_pre<0.3 ){
  //   depth = 0.15;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=0.3 && -z_pre<0.6 ){
  //   depth = 0.45;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=0.6 && -z_pre<0.9 ){
  //   depth = 0.75;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=0.9 && -z_pre<1.2 ){
  //   depth = 1.05;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // else if ( -z_pre>=1.2 && -z_pre<1.5 ){
  //   depth = 1.35;
  //   analysisManager->FillH1(0,depth);
  //   //      analysisManager->FillNtupleDColumn(6,depth);
  // }
  // }
  
  if (eventNumber != fEventNumber) {
    // G4cout << " Number of Scintillation Photons in previous event: "
    //        << fScintillationCounter << G4endl;
    fEventNumber = eventNumber;
    //  G4cout << "test : " << eventNumber << " " << fEventNumber << G4endl;
    fScintillationCounter = 0;
  }
  
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge()!=0.){
    if (z_pre == 1.5 ){
      //      G4cout << "test_pos : " << x_pre << " " << y_pre << " " << z_pre << G4endl;
      analysisManager->FillNtupleDColumn(2,x_pre);
      analysisManager->FillNtupleDColumn(3,y_pre);
      analysisManager->FillNtupleDColumn(4,z_pre);
    }
  }  
  
  //  if ( trackID==1 ){
  //      if ( PrePoint->GetStepStatus()==fGeomBoundary ){
  // if (volume == fDetConstruction->GetScintPV() ){
  //        G4cout << "track : " << trackID << G4endl;
  
  // G4cout << "test_mom : " << mom_preX << " " << mom_preY << " " << mom_preZ << G4endl;
  // G4cout << "test_pos2 : " << x_post << " " << y_post << " " << z_post << G4endl;
  // G4cout << "test_mom2 : " << mom_postX << " " << mom_postY << " " << mom_postZ << G4endl;
  //}
  
  
  G4String ParticleName = track->GetDynamicParticle()->
    GetParticleDefinition()->GetParticleName();
  
  
  // G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
  // G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
  
  
  // G4int parentID = track->GetTrackID();
  
  
  // if (trackID==1 && parentID==0){
  //   analysisManager->FillNtupleDColumn(3,x_pre);
  //   analysisManager->FillNtupleDColumn(4,y_pre);
  //   analysisManager->FillNtupleDColumn(5,z_pre);
  // }
  
  if (ParticleName == "opticalphoton") return;
  
  const std::vector<const G4Track*>* secondaries =
    step->GetSecondaryInCurrentStep();
  
  if (secondaries->size()>0) {
    for(unsigned int i=0; i<secondaries->size(); ++i) {
      if (secondaries->at(i)->GetParentID()>0) {
	if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
	   == G4OpticalPhoton::OpticalPhotonDefinition()){
	  if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
	      == "Scintillation")fScintillationCounter++;
	}
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
