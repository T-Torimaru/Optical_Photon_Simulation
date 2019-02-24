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
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make this appear first!
#include "G4Timer.hh"

#include "OpNoviceRunAction.hh"
#include "G4Run.hh"
#include "OpNoviceAnalysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceRunAction::OpNoviceRunAction()
  : G4UserRunAction()
{
  //  fTimer = new G4Timer;
  //  fRunMessenger = new WLSRunActionMessenger(this);
   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // //    analysisManager->SetNtupleMerging(true);
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFileName("test5");

    analysisManager->CreateNtuple("opTree","Optical_Photon");
    analysisManager->CreateNtupleDColumn("detPhoton");
    analysisManager->CreateNtupleDColumn("ScintEDep");
    analysisManager->CreateNtupleDColumn("EinitX");
    analysisManager->CreateNtupleDColumn("EinitY");
    analysisManager->CreateNtupleDColumn("EinitZ");
    analysisManager->CreateNtupleDColumn("excitedPhoton");
    analysisManager->CreateNtupleDColumn("TriggerEDep");
    analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceRunAction::~OpNoviceRunAction()
{
  delete G4AnalysisManager::Instance();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceRunAction::BeginOfRunAction(const G4Run* )
{
  // G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  // fTimer->Start();
   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   G4String fileName = "OpNovice";
   analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNoviceRunAction::EndOfRunAction(const G4Run* )
{
  // fTimer->Stop();
  // G4cout << "number of event = " << aRun->GetNumberOfEvent()
  //        << " " << *fTimer << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}
