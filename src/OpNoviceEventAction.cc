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
// $Id: B5EventAction.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file B5EventAction.cc
/// \brief Implementation of the B5EventAction class

#include "OpNoviceEventAction.hh"
#include "OpNoviceScintHit.hh"
#include "OpNoviceTriggerHit.hh"
#include "OpNoviceScintSD.hh"
#include "OpNoviceTriggerSD.hh"

#include "OpNoviceAnalysis.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceEventAction::OpNoviceEventAction()
  : G4UserEventAction(), 
    fmppcID(-1),
    ftriggerHCID(-1),
    fAbsHCID(-1),
    fEnergyDep(0.),
    fEnergyDep_trigger(0.)
{
  fMPPCCollID = 0;
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceEventAction::~OpNoviceEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceTriggerHitsCollection*
OpNoviceEventAction::GetHitsCollection(G4int hcID,
				       const G4Event* event) const
{
  OpNoviceTriggerHitsCollection* hitsCollection
    = static_cast<OpNoviceTriggerHitsCollection*>(
						  event->GetHCofThisEvent()->GetHC(hcID));
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("OpNoviceEventAction::GetHitsCollection()",
		"MyCode0003", FatalException, msg);
  }
  
  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceEventAction::BeginOfEventAction(const G4Event*)
{
  // initialisation per event 
  fEnergyDep = 0.;
  fEnergyDep_trigger = 0.;
  
}     

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceEventAction::EndOfEventAction(const G4Event* event)
{
  
  if (fmppcID==-1) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    G4String colName2 = "MPPCDS/mppcColl";
    fMPPCCollID = sdManager->GetCollectionID(colName2);
  }	
  
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  if (!hce) 
    {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found.\n"; 
      G4Exception("OpNoviceEventAction::EndOfEventAction()",
		  "OpNoviceCode001", JustWarning, msg);
      return;
    }   
  
  
  // Get the hit collections
  
  OpNoviceScintHitsCollection* mppcHC = (OpNoviceScintHitsCollection*)hce->GetHC(fMPPCCollID);
  if(!mppcHC) return;
  
  
  // Get hit information about photons that reached the detector in this event
  G4int n_hit1 = mppcHC->entries();
  
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  //Fill ntuple
  analysisManager->FillNtupleDColumn(0, n_hit1);
  analysisManager->FillNtupleDColumn(1, fEnergyDep);
  analysisManager->FillNtupleDColumn(6, fEnergyDep_trigger);
  analysisManager->AddNtupleRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
