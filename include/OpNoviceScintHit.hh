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

#ifndef OpNoviceScintHit_h
#define OpNoviceScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class OpNoviceScintHit : public G4VHit
{
public:
    OpNoviceScintHit(G4int i,G4double t);
    OpNoviceScintHit(const OpNoviceScintHit &right);
    virtual ~OpNoviceScintHit();

    const OpNoviceScintHit& operator=(const OpNoviceScintHit &right);
    int operator==(const OpNoviceScintHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void*aHit);
    
    void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    void Print();
    
    G4int GetID() const { return fId; }

    void SetTime(G4double val) { fTime = val; }
    G4double GetTime() const { return fTime; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }
    
private:
    G4int fId;
    G4double fTime;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
    const G4LogicalVolume* fPLogV;
};

typedef G4THitsCollection<OpNoviceScintHit> OpNoviceScintHitsCollection;

extern G4ThreadLocal G4Allocator<OpNoviceScintHit>* OpNoviceScintHitAllocator;

inline void* OpNoviceScintHit::operator new(size_t)
{
    if (!OpNoviceScintHitAllocator)
        OpNoviceScintHitAllocator = new G4Allocator<OpNoviceScintHit>;
    return (void*)OpNoviceScintHitAllocator->MallocSingle();
}

inline void OpNoviceScintHit::operator delete(void*aHit)
{
    OpNoviceScintHitAllocator->FreeSingle((OpNoviceScintHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
