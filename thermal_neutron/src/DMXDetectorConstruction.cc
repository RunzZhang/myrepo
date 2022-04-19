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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo
//                    (27th November 2001)
//
// DetectorConstruction program
// --------------------------------------------------------------

#include "DMXDetectorConstruction.hh"
#include "DMXDetectorMessenger.hh"

#include "DMXScintSD.hh"
#include "DMXPmtSD.hh"


#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include <math.h>

#define PI 3.14159265
using namespace std;

//Ryan added package
#include "G4OpticalSurface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::DMXDetectorConstruction()
{


  theUserLimitsForRoom     = 0;
  theUserLimitsForDetector = 0;
  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theRoomTimeCut      = 1000. * nanosecond;
  theMinEkine         = 250.0*eV; // minimum kinetic energy required in volume
  theRoomMinEkine     = 250.0*eV; // minimum kinetic energy required in volume

  //Zero the G4Cache objects to contain logical volumes
  LXeSD.Put(0);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::~DMXDetectorConstruction()
{
  delete theUserLimitsForRoom;
  delete theUserLimitsForDetector;
  delete detectorMessenger;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXDetectorConstruction::DefineMaterials()
{

 #include "DMXDetectorMaterial.icc"

 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DMXDetectorConstruction::Construct() {

  DefineMaterials();

  // DefineField();

  // make colours
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.85, .85, .85) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ;
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  orange  (.75, .55, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  G4Colour  lgreen  (0.0, .75, 0.0) ;
  G4Colour  green   (0.0, 1.0, 0.0) ;
  G4Colour  brown   (0.7, 0.4, 0.1) ;


  //  un-used colours:
  //  G4Colour  black   (0.0, 0.0, 0.0) ;



  // Universe


  G4double worldWidth  = 2*m ;
  G4double worldLength = 2*m ;
  G4double worldHeight = 1*m ;

  G4Box* world_box = new G4Box
     ("world_box", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight );
  world_log  = new G4LogicalVolume(world_box, vacuum_mat, "world_log");
  world_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),


     "world_phys", world_log, NULL, false,0);

G4double ArcontainerWidth = 1*m ;
G4double ArcontainerLength = 1*m ;
G4double ArcontainerHeight = 3*cm ;

G4Box* Ar_box = new G4Box("Ar_box",0.5*ArcontainerWidth, 0.5*ArcontainerLength, 0.5*ArcontainerHeight);
Arbox_log = new G4LogicalVolume(Ar_box, LAr_mat,"Arbox_log");
Arbox_phys =  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Arbox_log, "Arbox_phys", world_log, false, 0);
G4VisAttributes* argon_blue = new G4VisAttributes(blue);
argon_blue->SetVisibility(true);
Arbox_log -> SetVisAttributes(argon_blue);

//add opticalphoton cut
//G4double minEkin = 10*MeV;
//Arbox_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX,
 //                                         minEkin));


return world_phys;



}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXDetectorConstruction::ConstructSDandField()
{
    // ......................................................................
    // sensitive detectors ..................................................
    // ......................................................................

   if (LXeSD.Get() == 0)                                           // Aquí detecto los eventos del argón
      {
        G4String name="/DMXDet/LXeSD";
        DMXScintSD* aSD = new DMXScintSD(name);
        LXeSD.Put(aSD);
      }
    G4SDManager::GetSDMpointer()->AddNewDetector(LXeSD.Get());
 if(Arbox_log)
      SetSensitiveDetector(Arbox_log,LXeSD.Get());


    /*if (pmtSD.Get() == 0)                                        //Aquí detecto los eventos en el SiPM
    {
      G4String name="/DMXDet/pmtSD";
      DMXPmtSD* aSD = new DMXPmtSD(name);
      pmtSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD.Get());
  if (SIPM1_Si_log)
    SetSensitiveDetector(SIPM1_Si_log,pmtSD.Get());
  if (SIPM2_Si_log)
    SetSensitiveDetector(SIPM2_Si_log,pmtSD.Get());
  if (SIPM3_Si_log)
    SetSensitiveDetector(SIPM3_Si_log,pmtSD.Get());
  if (SIPM4_Si_log)
    SetSensitiveDetector(SIPM4_Si_log,pmtSD.Get());
  if (SIPM5_Si_log)
    SetSensitiveDetector(SIPM5_Si_log,pmtSD.Get());


*/


    return;

}
