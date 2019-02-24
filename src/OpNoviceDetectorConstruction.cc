
#include "OpNoviceDetectorConstruction.hh"

#include "G4ios.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "OpNoviceScintSD.hh"
#include "OpNoviceTriggerSD.hh"
#include "G4VSolid.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstruction::OpNoviceDetectorConstruction()
  : G4VUserDetectorConstruction(),fVisAttributes(),Apd_log(0)
{
  fExpHall_x = fExpHall_y = 30.5*mm;
  fExpHall_z =  1*cm;
  fScinti_x  = fScinti_y               =  30.0*mm;
  fScinti_z                            =  1.50*mm;
  Rmax       = 4.5*mm;
  Rmin       = 0.0*mm;
  startPhi   = 0.*deg;
  endPhi     = 360*deg;
  startTheta = 0.*deg;
  endTheta   = 180*deg;
  center    = 4.8*mm;
  Ax = Ay = 0.50*mm;
  Az = 0.085*mm;
  Epo_x = 1.3125*mm;
  Epo_y = 1.225*mm;
  Epo_z = 0.25*mm;
  z_epo = -3.55*mm;
  apd_z = 0.165*mm;
  trigger_x = trigger_y = trigger_z = 2.5*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNoviceDetectorConstruction::~OpNoviceDetectorConstruction(){
  for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) 
    {
      delete fVisAttributes[i];
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* OpNoviceDetectorConstruction::Construct()
{
  
  // ------------- Materials -------------
  
  G4double a, z, density;
  G4int nelements, natoms;
  
  G4Element* N  = new G4Element("Nitrogen" , "N", z=7 , a=14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   , "O", z=8 , a=16.00*g/mole);
  G4Element* C  = new G4Element("Carbon"   , "C", z=6 , a=12.01*g/mole);
  G4Element* H  = new G4Element("Hydrogen" , "H", z=1 , a=1.01*g/mole);
  G4Element* Si = new G4Element("Silicon"   ,"Si", z=14, a=28.0855*g/mole);
  
  // Air
  G4Material* air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  air->AddElement(N, 70.*perCent);
  air->AddElement(O, 30.*perCent);
  
  // Scintillator
  G4Material* Sci = new G4Material("Scintillator", density= 1.032*g/cm3, nelements=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);
  // Sci=G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  // Epoxy
  G4Material* epoxy = new G4Material("Epoxy", density= 1.85*g/cm3, nelements=3);
  epoxy->AddElement(C, natoms=18);
  epoxy->AddElement(H, natoms=20);
  epoxy->AddElement(O, natoms=2);
  
  // Si
  G4Material* apd = new G4Material("Si", density= 2.3290*g/cm3, nelements=1);
  apd->AddElement(Si, natoms=1);
  
  // Mylar
  //  Mylar=G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR");
  
  //
  // ------------ Generate & Add Material Properties Table ------------
  //
  const G4int nEntries = 32;
  
  G4double photonwavelength[nEntries]={400.,403.75,407.5,411.25,
				       415.,418.75,422.5,426.25,
				       430.,433.75,437.5,441.25,
				       445.,448.75,452.5,456.25,
				       460.,463.75,467.5,471.25,
				       475.,478.75,482.5,486.25,
				       490.,493.75,497.5,501.25,
				       505.,508.75,512.5,516.25};
  G4double Amplitude_fast[nEntries]={0.07,0.09,0.17,0.22,
  				     0.52,0.81,0.95,1.,
  				     0.8,0.72,0.71,0.7,
  				     0.68,0.63,0.5,0.42,
  				     0.37,0.33,0.3,0.28,
  				     0.23,0.21,0.19,0.18,
  				     0.16,0.15,0.13,0.11,
  				     0.1,0.09,0.08,0.07};
  G4double Amplitude_slow[nEntries]={1.,1.,1.,1.,1.,1.,1.,1.,
  				     1.,1.,1.,1.,1.,1.,1.,1.,
  			             1.,1.,1.,1.,1.,1.,1.,1.,
  				     1.,1.,1.,1.,1.,1.,1.,1.};
  G4double pde[nEntries]={0.230,0.231,0.232,0.233,
			  0.234,0.235,0.236,0.237,
			  0.238,0.239,0.240,0.241,
			  0.242,0.243,0.244,0.245,
			  0.246,0.247,0.248,0.249,
			  0.250,0.250,0.249,0.248,
			  0.247,0.246,0.245,0.244,
			  0.243,0.242,0.241,0.240};
  G4double photonenergy[nEntries];
  G4double Reflectivity_Si[nEntries];
  //  G4double Efficiency_Si[nEntries];
  G4double Epoxy_Refractive[nEntries];
  G4double Epoxy_Absorp_leng[nEntries];
  // G4double Mylar_Absorp_leng[nEntries];
  // G4double Mylar_Refractive[nEntries];
  G4double Air_Refractive[nEntries];
  
  G4double h_planck=6.62606e-34*joule*s;//4.136*eV*s;
  G4double c=3e8*m/s;
  G4double e=1.6021765e-19*joule/eV;
  
  for(G4int j=0;j<nEntries;j++){
    photonenergy[j]=(((h_planck*c/(photonwavelength[j]*1e-9))/e)*1e3)*eV;
  }
  
  //EJ-212 refractive index: fixed at 1.58 (no wavelength dependent)
  //from the Home Page of the ELJEN
  G4double Scinti_Refractive[nEntries];
  
  //EJ-212 absorption length: fixed at 80cm (no wavelength dependent)
  //from the thesis of Ms. Ieki
  G4double Scinti_Absorp_leng[nEntries];
  
  //if FAST/SLOW COMPONENTS are not defined, gammas are not generated! (to be understood)
  //from a paper in Ohio
  // const  G4int FastSlow = 4;
  // G4double FastSlowE[FastSlow]={2.*eV,2.87*eV,2.9*eV,3.47*eV};
  // G4double ScintiFast[FastSlow] = {1.,1.,1.,1.};
  // G4double ScintiSlow[FastSlow]= {0.,0.,1.,1.};
  
  for(G4int i=0;i<nEntries;i++){
    
    Reflectivity_Si[i] = 0.;
    //    Efficiency_Si[i] = 1.;
    Scinti_Refractive[i] = 1.58; //from ELJEN
    Epoxy_Refractive[i] = 1.55; //from hamamatsu
    Air_Refractive[i] = 1.;
    Scinti_Absorp_leng[i] = 80*cm;
    //    Epoxy_Absorp_leng[i] = 70*cm; //arbitrary
    Epoxy_Absorp_leng[i] = 220*cm; //Angela
    //    Mylar_Refractive[i] = 220*cm; //Angela
    //    Mylar_Absorp_leng[i] = 0.21*mm; //Angela
  }
  
  //Properties Table for EJ212
  
  G4MaterialPropertiesTable* EJ212_table = new G4MaterialPropertiesTable();
  
  EJ212_table->AddProperty("RINDEX", photonenergy, Scinti_Refractive, nEntries);
  EJ212_table->AddProperty("ABSLENGTH", photonenergy, Scinti_Absorp_leng, nEntries);
  //  EJ212_table->AddProperty("FASTCOMPONENT", FastSlowE, ScintiFast, FastSlow);
  //  EJ212_table->AddProperty("SLOWCOMPONENT", FastSlowE, ScintiSlow, FastSlow);
  EJ212_table->AddProperty("FASTCOMPONENT", photonenergy, Amplitude_fast, nEntries);
  EJ212_table->AddProperty("SLOWCOMPONENT", photonenergy, Amplitude_slow, nEntries);
  
  G4double LY = 10000./MeV;
  
  EJ212_table->AddConstProperty("SCINTILLATIONYIELD",LY); //from ELJEN
  
  EJ212_table->AddConstProperty("RESOLUTIONSCALE",1.);
  EJ212_table->AddConstProperty("FASTTIMECONSTANT", 2.7*ns);
  EJ212_table->AddConstProperty("SLOWTIMECONSTANT",2.7*ns);
  EJ212_table->AddConstProperty("YIELDRATIO",1.);
  Sci->GetIonisation()->SetBirksConstant(0.126*mm/MeV); //constant
  
  Sci->SetMaterialPropertiesTable(EJ212_table);
  
  //Properties Table for Air
  
  G4MaterialPropertiesTable* air_table = new G4MaterialPropertiesTable();
  
  air_table->AddProperty("RINDEX", photonenergy, Air_Refractive, nEntries);
  
  air->SetMaterialPropertiesTable(air_table);
  
  //Properties Table for Epoxy
  
  G4MaterialPropertiesTable* Epo_table = new G4MaterialPropertiesTable();
  
  Epo_table->AddProperty("RINDEX", photonenergy, Epoxy_Refractive, nEntries);
  Epo_table->AddProperty("ABSLENGTH", photonenergy, Epoxy_Absorp_leng, nEntries);   
  epoxy->SetMaterialPropertiesTable(Epo_table);
  
  //Properties Table for Silicon
  
  G4MaterialPropertiesTable* Silicon_table = new G4MaterialPropertiesTable();
  
  Silicon_table->AddProperty("REFLECTIVITY", photonenergy, Reflectivity_Si, nEntries);
  Silicon_table->AddProperty("EFFICIENCY", photonenergy, pde, nEntries);
  
  apd->SetMaterialPropertiesTable(Silicon_table);
  
  
  //
  // ------------- Volumes --------------
  
  // The experimental Hall
  //
  G4VSolid* expHall_box = new G4Box("World",fExpHall_x,fExpHall_y,fExpHall_z);
  
  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,air,"World_logical");
  
  G4VPhysicalVolume* expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);
  
  // The Scintillator
  //
  G4VSolid* Scinti_box = new G4Box("Tile",fScinti_x,fScinti_y,fScinti_z);
  
  G4LogicalVolume* Scinti_log
    = new G4LogicalVolume(Scinti_box,Sci,"Tile_logical");
  
  G4double maxStep = 0.1*mm;
  G4UserLimits *fstepLimit = new G4UserLimits(maxStep);
  Scinti_log->SetUserLimits(fstepLimit);
  
  Scinti_phys
    = new G4PVPlacement(0,G4ThreeVector(),Scinti_log,"Tile",
			expHall_log,true,0);
  
  // Trigger Counter
  //
  G4VSolid* Trigger_box = new G4Box("Trigger",trigger_x,trigger_y,trigger_z);
  
  Trigger_log
    = new G4LogicalVolume(Trigger_box,Sci,"Trigger_logical");
  
  Trigger_phys
    = new G4PVPlacement(0,G4ThreeVector(12.5*mm,12.5*mm,-6.*mm),Trigger_log,
  			"Trigger",expHall_log,false,0);
  
  // The Sphere
  //
  G4Sphere* sphere = new G4Sphere("Sphere",Rmin,Rmax,
				  startPhi,endPhi,startTheta,endTheta);
  
  G4LogicalVolume* sphere_log1
    = new G4LogicalVolume(sphere,air,"Sphere_logical1");
  
  G4VPhysicalVolume* sphere_phys1
    = new G4PVPlacement(0,G4ThreeVector(-15,15,center),
			sphere_log1,"Sphere",expHall_log,true,0);
  
  // G4LogicalVolume* sphere_log2
  //   = new G4LogicalVolume(sphere,air,"Sphere_logical2");
  
  // G4VPhysicalVolume* sphere_phys2
  //   = new G4PVPlacement(0,G4ThreeVector(0,0,center),
  // 			sphere_log2,"Sphere",expHall_log,true,0);
  
  // G4LogicalVolume* sphere_log3
  //   = new G4LogicalVolume(sphere,air,"Sphere_logical3");

  // G4VPhysicalVolume* sphere_phys3
  //   = new G4PVPlacement(0,G4ThreeVector(15,-15,center),
  // 			sphere_log3,"Sphere",expHall_log,true,0);

  
  // The Tile Design
  //
  G4SubtractionSolid* subtraction 
    = new G4SubtractionSolid("Scintillator-sphere",Scinti_box,sphere);
  
  //Epoxy
  //
  G4VSolid* epoxy_box = new G4Box("Epo",Epo_x,Epo_y,Epo_z);
  
  G4LogicalVolume* Epoxy_log1
    = new G4LogicalVolume(epoxy_box,epoxy,"Epo_logical1");
  
  // G4LogicalVolume* Epoxy_log2
  //   = new G4LogicalVolume(epoxy_box,epoxy,"Epo_logical2");
  
  G4VPhysicalVolume* Epoxy_phys1
    = new G4PVPlacement(0,G4ThreeVector(0,0,z_epo),Epoxy_log1,"Epo",
  			sphere_log1,false,0);
  
  // G4VPhysicalVolume* Epoxy_phys2
  //   = new G4PVPlacement(0,G4ThreeVector(0,0,z_epo),Epoxy_log2,"Epo",
  // 			sphere_log3,false,0);
  
  // G4VPhysicalVolume* Epoxy_phys1
  //   = new G4PVPlacement(0,G4ThreeVector(-0.09,0,z_epo),Epoxy_log1,"Epo",
  // 			sphere_log1,false,0);
  
  
  // APD
  //
  G4VSolid* Apd_box = new G4Box("Apd",Ax,Ay,Az);
  
  Apd_log
    = new G4LogicalVolume(Apd_box,apd,"apd_logical");
  
  // G4VPhysicalVolume* Apd_phys1
  //   = new G4PVPlacement(0,G4ThreeVector(0,0,apd_z),Apd_log,"Apd",
  // 			Epoxy_log1,false,0);
  
  // G4VPhysicalVolume* Apd_phys2
  //   = new G4PVPlacement(0,G4ThreeVector(0,0,apd_z),Apd_log,"Apd",
  // 			Epoxy_log2,false,0);
  
  G4VPhysicalVolume* Apd_phys1
    = new G4PVPlacement(0,G4ThreeVector(0.69,0.6,apd_z),Apd_log,"Apd",
  			Epoxy_log1,false,0);
  
  G4VPhysicalVolume* Apd_phys2
    = new G4PVPlacement(0,G4ThreeVector(0.69,-0.6,apd_z),Apd_log,"Apd",
  			Epoxy_log1,false,0);
  
  G4VPhysicalVolume* Apd_phys3
    = new G4PVPlacement(0,G4ThreeVector(-0.51,0.6,apd_z),Apd_log,"Apd",
  			Epoxy_log1,false,0);
  
  G4VPhysicalVolume* Apd_phys4
    = new G4PVPlacement(0,G4ThreeVector(-0.51,-0.6,apd_z),Apd_log,"Apd",
  			Epoxy_log1,false,0);
  
  
  // ------------- Surfaces --------------
  
  const G4int num = 2;
  const G4int num2 = 7;
  G4double ephoton[num] = {2.034*eV, 4.136*eV};
  G4double ephoton2[num2] = {2.48*eV, 2.59*eV, 2.70*eV, 2.82*eV, 2.95*eV, 3.11*eV, 3.27*eV};
  G4double refractiveIndex[num] = {1., 1.};
  G4double specularLobe[num]    = {0.9, 0.9};
  G4double specularSpike[num]   = {0., 0.};
  G4double backScatter[num]     = {0., 0.};
  G4double diffuseLobe[num]     = {0.1, 0.1};
  G4double reflectivity1[num]    = {0.95, 0.95};
  G4double reflectivity2[num]    = {1.,1.};  
  G4double reflectivity3[num]    = {0.,0.};  
  G4double efficiency[num]      = {0., 0.};
  //  G4double efficiency2[num2] = {0.454, 0.473, 0.485, 0.473, 0.454, 0.419, 0.377};
  G4double efficiency2[num2] = {1., 1., 1., 1., 1., 1., 1.};
  G4double sigma_alpha = 0.1;
  
  //Scintillator -> Air (-> Reflector)
  //
  G4OpticalSurface* ScintToAir = new G4OpticalSurface("SciSurface1");
  G4LogicalBorderSurface* sciSurface1 =
    new G4LogicalBorderSurface("SciSurface1",Scinti_phys,expHall_phys,ScintToAir);
  
  ScintToAir->SetType(dielectric_dielectric);
  ScintToAir->SetFinish(polishedbackpainted);
  ScintToAir->SetModel(unified);
  //  ScintToAir->SetSigmaAlpha(sigma_alpha);
  
  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  
  myST1->AddProperty("RINDEX",                ephoton, refractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);
  myST1->AddProperty("DIFFUSELOBE",           ephoton, diffuseLobe,     num);
  myST1->AddProperty("REFLECTIVITY",          ephoton, reflectivity1,   num);
  //  myST1->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);
  
  ScintToAir->SetMaterialPropertiesTable(myST1);
  
  //Scintillator -> Sphere (Air)
  //
  G4OpticalSurface* ScintToSphere1 = new G4OpticalSurface("SciSurface2");
  G4LogicalBorderSurface* sciSurface2 =
    new G4LogicalBorderSurface("SciSurface2", Scinti_phys, sphere_phys1, ScintToSphere1);
  
  ScintToSphere1->SetType(dielectric_dielectric);
  ScintToSphere1->SetFinish(polished);
  // ScintToSphere->SetModel(unified);
  // ScintToSphere->SetSigmaAlpha(sigma_alpha);
  
  // G4MaterialPropertiesTable* myST3 = new G4MaterialPropertiesTable();
  
  // myST3->AddProperty("REFLECTIVITY",          ephoton, reflectivity2,   num);
  // myST3->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);
  
  // ScintToSphere->SetMaterialPropertiesTable(myST3);
  
  // //Scintillator -> Sphere (Air)
  // //
  // G4OpticalSurface* ScintToSphere2 = new G4OpticalSurface("SciSurface3");
  // G4LogicalBorderSurface* sciSurface3 =
  //   new G4LogicalBorderSurface("SciSurface3", Scinti_phys, sphere_phys2, ScintToSphere2);
  
  // ScintToSphere2->SetType(dielectric_dielectric);
  // ScintToSphere2->SetFinish(polished);
  // // ScintToSphere->SetModel(unified);
  // // ScintToSphere->SetSigmaAlpha(sigma_alpha);

  // // G4MaterialPropertiesTable* myST3 = new G4MaterialPropertiesTable();

  // // myST3->AddProperty("REFLECTIVITY",          ephoton, reflectivity2,   num);
  // // myST3->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);

  // // // ScintToSphere->SetMaterialPropertiesTable(myST3);

  // //Scintillator -> Sphere (Air)
  // //
  // G4OpticalSurface* ScintToSphere3 = new G4OpticalSurface("SciSurface4");
  // G4LogicalBorderSurface* sciSurface4 =
  //   new G4LogicalBorderSurface("SciSurface4", Scinti_phys, sphere_phys3, ScintToSphere3);
  
  // ScintToSphere3->SetType(dielectric_dielectric);
  // ScintToSphere3->SetFinish(polished);
  // // ScintToSphere->SetModel(unified);
  // // ScintToSphere->SetSigmaAlpha(sigma_alpha);

  // // G4MaterialPropertiesTable* myST3 = new G4MaterialPropertiesTable();

  // // myST3->AddProperty("REFLECTIVITY",          ephoton, reflectivity2,   num);
  // // myST3->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);

  // // ScintToSphere->SetMaterialPropertiesTable(myST3);

  //Epoxy
  //
  G4OpticalSurface* surface_Epoxy = new G4OpticalSurface("EpoxySurface");
  G4LogicalSkinSurface* epoxySurface =
    new G4LogicalSkinSurface("EpoxySurface", Epoxy_log1, surface_Epoxy);
  
  surface_Epoxy->SetType(dielectric_dielectric);
  surface_Epoxy->SetFinish(polished);
  surface_Epoxy->SetModel(unified);
  surface_Epoxy->SetSigmaAlpha(sigma_alpha);
  //  surface_Epoxy->SetPolish(1);
  
  // G4MaterialPropertiesTable* myST9 = new G4MaterialPropertiesTable();

  // myST9->AddProperty("REFLECTIVITY",          ephoton, reflectivity2,   num);
  // myST9->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);

  // // surface_Epoxy->SetMaterialPropertiesTable(myST9);

  // G4OpticalSurface* surface_Epoxy2 = new G4OpticalSurface("EpoxySurface2");
  // G4LogicalSkinSurface* epoxySurface2 =
  //   new G4LogicalSkinSurface("EpoxySurface2", Epoxy_log2, surface_Epoxy2);
  
  // surface_Epoxy2->SetType(dielectric_dielectric);
  // surface_Epoxy2->SetFinish(polished);
  // surface_Epoxy2->SetModel(unified);
  // surface_Epoxy2->SetSigmaAlpha(sigma_alpha);
  // //  surface_Epoxy->SetPolish(1);

  // // G4MaterialPropertiesTable* myST9 = new G4MaterialPropertiesTable();

  // // myST9->AddProperty("REFLECTIVITY",          ephoton, reflectivity2,   num);
  // // myST9->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);

  // surface_Epoxy->SetMaterialPropertiesTable(myST9);

  // APD
  //
  G4OpticalSurface* APDSur = new G4OpticalSurface("APDSurface");
  G4LogicalSkinSurface* APDSurface =
    new G4LogicalSkinSurface("APDSurface",  Apd_log, APDSur);
  
  APDSur->SetType(dielectric_metal);
  // APDSur->SetFinish(polishedlumirrorair);
  // APDSur->SetModel(LUT);
  
  G4MaterialPropertiesTable* myST10 = new G4MaterialPropertiesTable();
  
  myST10->AddProperty("REFLECTIVITY",          ephoton, reflectivity3,   num);
  myST10->AddProperty("EFFICIENCY",           photonenergy,   pde,   nEntries);
  
  APDSur->SetMaterialPropertiesTable(myST10);
  
  // //Sphere (Air) -> Scintillator
  // //
  // G4OpticalSurface* SphereToScint = new G4OpticalSurface("SciSurface4");
  // G4LogicalBorderSurface* sciSurface4 =
  //   new G4LogicalBorderSurface("SciSurface4", sphere_phys, Scinti_phys, SphereToScint);
  
  // SphereToScint->SetType(dielectric_dielectric);
  // SphereToScint->SetFinish(polished);
  // SphereToScint->SetModel(unified);

  // G4MaterialPropertiesTable* myST4 = new G4MaterialPropertiesTable();
  
  // myST4->AddProperty("REFLECTIVITY",          ephoton, reflectivity2,   num);
  // myST4->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);
  
  // ScintToSphere->SetMaterialPropertiesTable(myST4);

  // //Sphere (Air) -> Epoxy
  // //
  // G4OpticalSurface* SphereToEpoxy = new G4OpticalSurface("EpoSurface");
  // G4LogicalBorderSurface* epoSurface =
  //   new G4LogicalBorderSurface("EpoSurface", sphere_phys, Epoxy_phys, SphereToEpoxy);
  
  // SphereToEpoxy->SetType(dielectric_dielectric);
  // SphereToEpoxy->SetFinish(polished);
  // SphereToEpoxy->SetModel(glisur);
  // SphereToEpoxy->SetPolish(1);

  // G4MaterialPropertiesTable* myST5 = new G4MaterialPropertiesTable();

  // myST5->AddProperty("REFLECTIVITY",          ephoton, reflectivity2,   num);
  // myST5->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);

  // SphereToEpoxy->SetMaterialPropertiesTable(myST5);

  // //Epoxy -> Sphere (Air)
  // //
  // G4OpticalSurface* EpoxyToSphere = new G4OpticalSurface("SphereSurface");
  // G4LogicalBorderSurface* sphereSurface =
  //   new G4LogicalBorderSurface("SphereSurface", Epoxy_phys, Sphere_phys, EpoxyToSphere);
  
  // EpoxyToSphere->SetType(dielectric_dielectric);
  // EpoxyToSphere->SetFinish(polished);
  // EpoxyToSphere->SetModel(glisur);
  // EpoxyToSphere->SetPolish(1);

  // G4MaterialPropertiesTable* myST6 = new G4MaterialPropertiesTable();

  // myST6->AddProperty("REFLECTIVITY",          ephoton, reflectivity2,   num);
  // myST6->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);

  // EpoxyToSphere->SetMaterialPropertiesTable(myST6);



  // //Epoxy -> APD
  // //
  // G4OpticalSurface* EpoxyToApd = new G4OpticalSurface("ApdSurface");
  // G4LogicalBorderSurface* apdSurface =
  //   new G4LogicalBorderSurface("ApdSurface", Epoxy_phys, Apd_phys, EpoxyToApd);
 
  // EpoxyToApd->SetType(dielectric_metal);
  // EpoxyToApd->SetFinish(polishedlumirrorair);
  // EpoxyToApd->SetModel(LUT);

  // G4MaterialPropertiesTable* myST7 = new G4MaterialPropertiesTable();

  // myST7->AddProperty("REFLECTIVITY",          ephoton, reflectivity3,   num);
  // myST7->AddProperty("EFFICIENCY",           ephoton2,   efficiency2,   num);

  // EpoxyToApd->SetMaterialPropertiesTable(myST7);

  // //Sphere -> Air
  // //
  // G4OpticalSurface* SphereToAir = new G4OpticalSurface("SpheSurface");
  // G4LogicalBorderSurface* SpheSurface =
  //   new G4LogicalBorderSurface("SpheSurface", sphere_phys, expHall_phys, SphereToAir);
 
  // SphereToAir->SetType(dielectric_dielectric);
  // //  SphereToAir->SetFinish(polishedlumirrorair);
  // SphereToAir->SetModel(glisur);
  // SphereToAir->SetPolish(1);

  // G4MaterialPropertiesTable* myST8 = new G4MaterialPropertiesTable();

  // myST8->AddProperty("REFLECTIVITY",          ephoton, reflectivity3,   num);
  // myST8->AddProperty("EFFICIENCY",           ephoton,   efficiency,   num2);

  // SphereToAir->SetMaterialPropertiesTable(myST8);
  
  //Surface of trigger counter
  //
  G4OpticalSurface* SurfaceTrigger = new G4OpticalSurface("TriggerSurface");
  G4LogicalSkinSurface* triggersurface =
    new G4LogicalSkinSurface("TriggerSurface",Trigger_log,SurfaceTrigger);
  
  SurfaceTrigger->SetType(dielectric_dielectric);
  SurfaceTrigger->SetFinish(polishedbackpainted);
  SurfaceTrigger->SetModel(unified);
  //  ScintToAir->SetSigmaAlpha(sigma_alpha);
  
  G4MaterialPropertiesTable* myST20 = new G4MaterialPropertiesTable();
  
  myST20->AddProperty("RINDEX",                ephoton, refractiveIndex, num);
  myST20->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
  myST20->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
  myST20->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);
  myST20->AddProperty("DIFFUSELOBE",           ephoton, diffuseLobe,     num);
  myST20->AddProperty("REFLECTIVITY",          ephoton, reflectivity1,   num);
  //  myST1->AddProperty("EFFICIENCY",            ephoton, efficiency,      num);
  
  SurfaceTrigger->SetMaterialPropertiesTable(myST20);
  
  //---------- visualization attributes ---------------------------------
  
  G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1,0,0));
  visAttributes->SetVisibility(false);
  expHall_log->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0,1,1));
  Scinti_log->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0,0,1));
  Trigger_log->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0,0,0));
  sphere_log1->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  // visAttributes = new G4VisAttributes(G4Colour(0,0,0));
  // sphere_log2->SetVisAttributes(visAttributes);
  // fVisAttributes.push_back(visAttributes);
  
  // visAttributes = new G4VisAttributes(G4Colour(0,0,0));
  // sphere_log3->SetVisAttributes(visAttributes);
  // fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0,1,0));
  Epoxy_log1->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  // visAttributes = new G4VisAttributes(G4Colour(0,1,0));
  // Epoxy_log2->SetVisAttributes(visAttributes);
  // fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(1,1,1));
  Apd_log->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  //always return the physical World
  return expHall_phys;
}

void OpNoviceDetectorConstruction::ConstructSDandField()
{
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname1;
  G4String SDname2;
  
  G4VSensitiveDetector* MPPCDS
    = new OpNoviceScintSD(SDname1="/MPPCDS");
  SDman->AddNewDetector(MPPCDS);
  Apd_log->SetSensitiveDetector(MPPCDS);
  
  G4VSensitiveDetector* TRIGGER
    = new OpNoviceScintSD(SDname2="/TRIGGER");
  SDman->AddNewDetector(TRIGGER);
  Trigger_log->SetSensitiveDetector(TRIGGER);
  
}
