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

#ifndef OpNoviceDetectorConstruction_h
#define OpNoviceDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include <vector>
#include "G4Material.hh"

class G4VSensitiveDetector;
class G4VisAttributes;
class G4VPhysicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class OpNoviceDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    OpNoviceDetectorConstruction();
    virtual ~OpNoviceDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    const G4VPhysicalVolume* GetScintPV() const;
    const G4VPhysicalVolume* GetTriggerPV() const;

  private:
    std::vector<G4VisAttributes*> fVisAttributes;

    G4VPhysicalVolume* Scinti_phys;
    G4VPhysicalVolume* Trigger_phys;

    G4LogicalVolume* Apd_log;
    G4LogicalVolume* Trigger_log;

    G4double fExpHall_x;
    G4double fExpHall_y;
    G4double fExpHall_z;

    G4double fScinti_x;
    G4double fScinti_y;
    G4double fScinti_z;

    G4double Rmax;
    G4double Rmin;
    G4double startPhi;
    G4double endPhi;
    G4double startTheta;
    G4double endTheta;
    G4double center;
    G4double Ax,Ay,Az;
    G4double Epo_x;
    G4double Epo_y;
    G4double Epo_z;
    G4double z_epo;
    G4double apd_z;

    G4double trigger_x;
    G4double trigger_y;
    G4double trigger_z;

    G4Material* Sci;

};

inline const G4VPhysicalVolume* OpNoviceDetectorConstruction::GetScintPV() const {
  return Scinti_phys;
}
inline const G4VPhysicalVolume* OpNoviceDetectorConstruction::GetTriggerPV() const {
  return Trigger_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*OpNoviceDetectorConstruction_h*/
