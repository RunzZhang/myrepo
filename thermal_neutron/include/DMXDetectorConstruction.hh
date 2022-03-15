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
// DetectorConstruction header
// --------------------------------------------------------------

#ifndef DMXDetectorConstruction_h
#define DMXDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class G4UserLimits;

class DMXScintSD;
class DMXPmtSD;

class DMXDetectorMessenger;

class DMXDetectorConstruction : public G4VUserDetectorConstruction 
{
public:

  DMXDetectorConstruction();
  ~DMXDetectorConstruction();

public:

  G4VPhysicalVolume* Construct();
  void ConstructSDandField();

  void SetRoomEnergyCut(G4double);
  void SetEnergyCut(G4double);
  void SetTimeCut(G4double);
  void SetRoomTimeCut(G4double);
 
private:

  void DefineMaterials();

  G4UserLimits*    theUserLimitsForRoom; 
  G4UserLimits*    theUserLimitsForDetector; 
  //  G4UserLimits*    theUserLimitsForXenon; 

  G4double         theMaxTimeCuts;
  G4double         theMaxStepSize;
  G4double         theDetectorStepSize;
  G4double         theMinEkine;
  G4double         theRoomMinEkine;
  
  G4double         theRoomTimeCut;


#include "DMXDetectorMaterial.ihh"  // materials used

  G4double sourceZ;

  G4LogicalVolume*   Arbox_log;        // pointers
  G4VPhysicalVolume* Arbox_phys;  

#include "DMXDetectorRoom.ihh"

  G4LogicalVolume*   world_log;        // pointers
  G4VPhysicalVolume* world_phys; 
	

  /*G4LogicalVolume*   Vacuum_vessel_log;
  G4VPhysicalVolume* Vacuum_vessel_phys;
  G4LogicalVolume*   Camera_port_log;
  G4VPhysicalVolume* Camera_port_phys;
  G4LogicalVolume*   Camera_port2_log;
  G4VPhysicalVolume* Camera_port2_phys;
  G4LogicalVolume*   Camera_port3_log;
  G4VPhysicalVolume* Camera_port3_phys;
  G4LogicalVolume*   Camera_System_log;
  G4VPhysicalVolume* Camera_System_phys;
  G4LogicalVolume*   Sapphire_log;
  G4VPhysicalVolume* Sapphire_phys;
  G4LogicalVolume*   Sapphire_Ssteal_log;
  G4VPhysicalVolume* Sapphire_Ssteal_phys;
  G4LogicalVolume*   Lens_Holder_log;
  G4VPhysicalVolume* Lens_Holder_phys;
  G4LogicalVolume*   Lens_log;
  G4VPhysicalVolume* Lens_phys;  
  G4LogicalVolume*   Alignment_Ring_log;
  G4VPhysicalVolume* Alignment_Ring_phys;    
  G4LogicalVolume*   Front_Nanoguide_Holder_log;
  G4VPhysicalVolume* Front_Nanoguide_Holder_phys; 
  G4LogicalVolume*   Peak_Rod1_log;
  G4VPhysicalVolume* Peak_Rod1_phys; 
  G4LogicalVolume*   Peak_Rod2_log;
  G4VPhysicalVolume* Peak_Rod2_phys; 
  G4LogicalVolume*   Rear_Nanoguide_Holder_log;
  G4VPhysicalVolume* Rear_Nanoguide_Holder_phys;   
  G4LogicalVolume*   Nanoguide_log;
  G4VPhysicalVolume* Nanoguide_phys;  
  G4LogicalVolume*   Sensor_Plate_log;
  G4VPhysicalVolume* Sensor_Plate_phys;    
  G4LogicalVolume*   Sensor_Cover_log;
  G4VPhysicalVolume* Sensor_Cover_phys;  
  G4LogicalVolume*   Screw1_log;
  G4VPhysicalVolume* Screw1_phys;   
  G4LogicalVolume*   Screw2_log;
  G4VPhysicalVolume* Screw2_phys;  
  G4LogicalVolume*   Screw3_log;
  G4VPhysicalVolume* Screw3_phys; 
  G4LogicalVolume*   Screw4_log;
  G4VPhysicalVolume* Screw4_phys;     
  G4LogicalVolume*   Camera_spring1_log;
  G4VPhysicalVolume* Camera_spring1_phys;   
  G4LogicalVolume*   Camera_spring2_log;
  G4VPhysicalVolume* Camera_spring2_phys;   
  G4LogicalVolume*   Camera_log;
  G4VPhysicalVolume* Camera_phys;   
  G4LogicalVolume*   Camera_PCB_log;
  G4VPhysicalVolume* Camera_PCB_phys;  
  G4LogicalVolume*   looking_log;
  G4VPhysicalVolume* looking_phys;  
  
  G4LogicalVolume*   Camera_System1_log;
  G4VPhysicalVolume* Camera_System1_phys;
  G4LogicalVolume*   Sapphire1_log;
  G4VPhysicalVolume* Sapphire1_phys; 
  G4LogicalVolume*   Sapphire1_Ssteal_log;
  G4VPhysicalVolume* Sapphire1_Ssteal_phys;
  G4LogicalVolume*   Plastic_Flange_log;
  G4VPhysicalVolume* Plastic_Flange_phys;
  G4LogicalVolume*   Lens_Holder2_log;
  G4VPhysicalVolume* Lens_Holder2_phys; 
  G4LogicalVolume*   Lens1_log;
  G4VPhysicalVolume* Lens1_phys;  
  G4LogicalVolume*   Adjustment_log;
  G4VPhysicalVolume* Adjustment_phys; 
  G4LogicalVolume*   Adjustment1_log;
  G4VPhysicalVolume* Adjustment1_phys; 
  G4LogicalVolume*   Top_Plastic_log;
  G4VPhysicalVolume* Top_Plastic_phys;     
  G4LogicalVolume*   Top_Plastic1_log;
  G4VPhysicalVolume* Top_Plastic1_phys;
  G4LogicalVolume*   Aspheric_lens_log;
  G4VPhysicalVolume* Aspheric_lens_phys;    
  G4LogicalVolume*   Aspheric_lens1_log;
  G4VPhysicalVolume* Aspheric_lens1_phys; 
  G4LogicalVolume*   Lens_Holder3_log;
  G4VPhysicalVolume* Lens_Holder3_phys;   
  G4LogicalVolume*   Lens_Holder4_log;
  G4VPhysicalVolume* Lens_Holder4_phys;
  G4LogicalVolume*    Iris_Holder_log;
  G4VPhysicalVolume*  Iris_Holder_phys; 
  G4LogicalVolume*    Iris_Holder1_log;
  G4VPhysicalVolume*  Iris_Holder1_phys; 
  G4LogicalVolume*    Iris_log;
  G4VPhysicalVolume*  Iris_phys;  
  G4LogicalVolume*    Iris_brass_log;
  G4VPhysicalVolume*  Iris_brass_phys; 
  G4LogicalVolume*    Iris_blades_log;
  G4VPhysicalVolume*  Iris_blades_phys; 
  G4LogicalVolume*    Sensor_Holder_log;
  G4VPhysicalVolume*  Sensor_Holder_phys;  
  G4LogicalVolume*     Copper_Rod1_log;
  G4VPhysicalVolume*   Copper_Rod1_phys;  
  G4LogicalVolume*     Copper_Rod2_log;
  G4VPhysicalVolume*   Copper_Rod2_phys;   
  G4LogicalVolume*     SS_Rod1_log;
  G4VPhysicalVolume*   SS_Rod1_phys;  
  G4LogicalVolume*     SS_Rod2_log;
  G4VPhysicalVolume*   SS_Rod2_phys;  
  G4LogicalVolume*   Camera1_log;
  G4VPhysicalVolume* Camera1_phys;   
  G4LogicalVolume*   Camera_PCB1_log;
  G4VPhysicalVolume* Camera_PCB1_phys; 

  G4LogicalVolume*   Holder1_VV_log;
  G4VPhysicalVolume* Holder1_VV_phys;
  G4LogicalVolume*   Holder2_VV_log;
  G4VPhysicalVolume* Holder2_VV_phys;
  G4LogicalVolume*   Holder3_VV_log;
  G4VPhysicalVolume* Holder3_VV_phys;
  G4LogicalVolume*   Inside_vacuum_vessel_log;
  G4VPhysicalVolume* Inside_vacuum_vessel_phys;
  G4LogicalVolume*   pressure_vessel_log;
  G4VPhysicalVolume* pressure_vessel_phys;
  G4LogicalVolume*   PV_spool_log;
  G4VPhysicalVolume* PV_spool_phys;
  G4LogicalVolume*   hydraulic_fluid_log;
  G4VPhysicalVolume* hydraulic_fluid_phys;
  G4LogicalVolume*   contenedor_emanacion_log;
  G4VPhysicalVolume* contenedor_emanacion_phys;
  G4LogicalVolume*   contenedorp_log;
  G4VPhysicalVolume* contenedorp_phys;
  G4LogicalVolume*   outer_jar_log;
  G4VPhysicalVolume* outer_jar_phys;
  G4LogicalVolume*   LAr_log;
  G4VPhysicalVolume* LAr_phys;
  G4LogicalVolume*   inner_jar_log;
  G4VPhysicalVolume* inner_jar_phys;
  G4LogicalVolume*   top_flange_log;
  G4VPhysicalVolume* top_flange_phys;
  G4LogicalVolume*   Top_sf_log;
  G4VPhysicalVolume* Top_sf_phys;
  G4LogicalVolume*   Guide_rod_flange_log;
  G4VPhysicalVolume* Guide_rod_flange_phys;
  G4LogicalVolume*   OJ_Spacer_log;
  G4VPhysicalVolume* OJ_Spacer_phys;
  G4LogicalVolume*   Bellows_weldment_log;
  G4VPhysicalVolume* Bellows_weldment_phys;
  G4LogicalVolume*   CF4_inside_log;
  G4VPhysicalVolume* CF4_inside_phys;
  G4LogicalVolume*   CF4_outside_log;
  G4VPhysicalVolume* CF4_outside_phys;
  G4LogicalVolume*   HDPE_inner_jar_log;
  G4VPhysicalVolume* HDPE_inner_jar_phys;
  G4LogicalVolume*   Bottom_Spacer_log;
  G4VPhysicalVolume* Bottom_Spacer_phys;
  G4LogicalVolume*   Bottom_Flange_log;
  G4VPhysicalVolume* Bottom_Flange_phys;  
  G4LogicalVolume*   Base_SF_log;
  G4VPhysicalVolume* Base_SF_phys; 
  G4LogicalVolume*   Side_support_log ;
  G4VPhysicalVolume* Side_support_phys; 
  G4LogicalVolume*   Guide_Rod_log;
  G4VPhysicalVolume* Guide_Rod_phys;  
  G4LogicalVolume*   Hyspan_bellow_log;
  G4VPhysicalVolume* Hyspan_bellow_phys;   
  G4LogicalVolume*   HDPE_pressure_vessel_log;
  G4VPhysicalVolume* HDPE_pressure_vessel_phys;
  G4LogicalVolume*   HDPE_CF4_log;
  G4VPhysicalVolume* HDPE_CF4_phys;
  G4LogicalVolume*   reflector_PTFE_log;
  G4VPhysicalVolume* reflector_PTFE_phys;
  G4LogicalVolume*   reflector_Ny_log;
  G4VPhysicalVolume* reflector_Ny_phys;
  G4LogicalVolume*   reflector_Cu_log;
  G4VPhysicalVolume* reflector_Cu_phys;
  G4LogicalVolume*   reflector_top_Cu_log;
  G4VPhysicalVolume* reflector_top_Cu_phys;
  G4LogicalVolume*   reflector_top_PTFE_log;
  G4VPhysicalVolume* reflector_top_PTFE_phys;
  G4LogicalVolume*   reflector_top_Ny_log;
  G4VPhysicalVolume* reflector_top_Ny_phys;
  G4LogicalVolume*   SIPM1_Si_log;
  G4VPhysicalVolume* SIPM1_Si_phys;
  G4LogicalVolume*   SIPM2_Si_log;
  G4VPhysicalVolume* SIPM2_Si_phys;
  G4LogicalVolume*   SIPM3_Si_log;
  G4VPhysicalVolume* SIPM3_Si_phys;
  G4LogicalVolume*   SIPM4_Si_log;
  G4VPhysicalVolume* SIPM4_Si_phys;
  G4LogicalVolume*   SIPM5_Si_log;
  G4VPhysicalVolume* SIPM5_Si_phys;
  G4LogicalVolume*   SiPM1_log;
  G4VPhysicalVolume* SiPM1_phys;
  G4LogicalVolume*   SiPM2_log;
  G4VPhysicalVolume* SiPM2_phys;
  G4LogicalVolume*   SiPM3_log;
  G4VPhysicalVolume* SiPM3_phys;
  G4LogicalVolume*   SiPM4_log;
  G4VPhysicalVolume* SiPM4_phys;
  G4LogicalVolume*   SiPM5_log;
  G4VPhysicalVolume* SiPM5_phys;
  G4LogicalVolume*   SiPM1_Inside_log;
  G4VPhysicalVolume* SiPM1_Inside_phys;
  G4LogicalVolume*   SiPM2_Inside_log;
  G4VPhysicalVolume* SiPM2_Inside_phys;
  G4LogicalVolume*   SiPM3_Inside_log;
  G4VPhysicalVolume* SiPM3_Inside_phys;
  G4LogicalVolume*   SiPM4_Inside_log;
  G4VPhysicalVolume* SiPM4_Inside_phys;
  G4LogicalVolume*   SiPM5_Inside_log;
  G4VPhysicalVolume* SiPM5_Inside_phys;
  G4LogicalVolume*   SiPM_Holder1_log;
  G4VPhysicalVolume* SiPM_Holder1_phys;
  G4LogicalVolume*   SiPM_Holder2_log;
  G4VPhysicalVolume* SiPM_Holder2_phys;
  G4LogicalVolume*   SiPM_Holder3_log;
  G4VPhysicalVolume* SiPM_Holder3_phys;
  G4LogicalVolume*   SiPM_Holder4_log;
  G4VPhysicalVolume* SiPM_Holder4_phys;
  G4LogicalVolume*   SiPM_Holder5_log;
  G4VPhysicalVolume* SiPM_Holder5_phys;
  G4LogicalVolume*   Quartz1_log;
  G4VPhysicalVolume* Quartz1_phys;
  G4LogicalVolume*   Quartz2_log;
  G4VPhysicalVolume* Quartz2_phys;
  G4LogicalVolume*   Quartz3_log;
  G4VPhysicalVolume* Quartz3_phys;
  G4LogicalVolume*   Quartz4_log;
  G4VPhysicalVolume* Quartz4_phys;
  G4LogicalVolume*   Quartz5_log;
  G4VPhysicalVolume* Quartz5_phys;  
  G4LogicalVolume*   SiPM_Plastic1_log;
  G4VPhysicalVolume* SiPM_Plastic1_phys;
  G4LogicalVolume*   SiPM_PCB1_log;
  G4VPhysicalVolume* SiPM_PCB1_phys;
  G4LogicalVolume*   SiPM_Plastic2_log;
  G4VPhysicalVolume* SiPM_Plastic2_phys;
  G4LogicalVolume*   SiPM_PCB2_log;
  G4VPhysicalVolume* SiPM_PCB2_phys;
  G4LogicalVolume*   SiPM_Plastic3_log;
  G4VPhysicalVolume* SiPM_Plastic3_phys;
  G4LogicalVolume*   SiPM_PCB3_log;
  G4VPhysicalVolume* SiPM_PCB3_phys;
  G4LogicalVolume*   SiPM_Plastic4_log;
  G4VPhysicalVolume* SiPM_Plastic4_phys;
  G4LogicalVolume*   SiPM_PCB4_log;
  G4VPhysicalVolume* SiPM_PCB4_phys;  
  G4LogicalVolume*   SiPM_Plastic5_log;
  G4VPhysicalVolume* SiPM_Plastic5_phys;
  G4LogicalVolume*   SiPM_PCB5_log;
  G4VPhysicalVolume* SiPM_PCB5_phys;  

  G4LogicalVolume*   Emanation_cable_log;
  G4VPhysicalVolume* Emanation_cable_phys;
  G4LogicalVolume*   Emanation_cable2_log;
  G4VPhysicalVolume* Emanation_cable2_phys;
  G4LogicalVolume*   Emanation_cable3_log;
  G4VPhysicalVolume* Emanation_cable3_phys;

  G4LogicalVolume*   Emanation_Holder_log;
  G4VPhysicalVolume* Emanation_Holder_phys;

  G4LogicalVolume*   cover_cable_log;
  G4VPhysicalVolume* cover_cable_phys;
  G4LogicalVolume*   cover1_cable_log;
  G4VPhysicalVolume* cover1_cable_phys;
  G4LogicalVolume*   cover2_cable_log;
  G4VPhysicalVolume* cover2_cable_phys;
  G4LogicalVolume*   shield_cable_log;
  G4VPhysicalVolume* shield_cable_phys;
  G4LogicalVolume*   shield1_cable_log;
  G4VPhysicalVolume* shield1_cable_phys;
  G4LogicalVolume*   shield2_cable_log;
  G4VPhysicalVolume* shield2_cable_phys;
  G4LogicalVolume*   inside_cable_log;
  G4VPhysicalVolume* inside_cable_phys;
  G4LogicalVolume*   inside1_cable_log;
  G4VPhysicalVolume* inside1_cable_phys;
  G4LogicalVolume*   inside2_cable_log;
  G4VPhysicalVolume* inside2_cable_phys;
  G4LogicalVolume*   core_cable_log;
  G4VPhysicalVolume* core_cable_phys;
  G4LogicalVolume*   core1_cable_log;
  G4VPhysicalVolume* core1_cable_phys;
  G4LogicalVolume*   core2_cable_log;
  G4VPhysicalVolume* core2_cable_phys;

  G4LogicalVolume*   LED_log;
  G4VPhysicalVolume* LED_phys;
 


  G4LogicalVolume*   Cu_log;
  G4VPhysicalVolume* Piezo_Cu_phys;
  G4LogicalVolume*   Vacio_log;
  G4VPhysicalVolume* Vacio_phys;
  G4LogicalVolume*   PBC_log;
  G4VPhysicalVolume* PBC_phys;
  G4LogicalVolume*   Piezo_log;
  G4VPhysicalVolume* Piezo_phys;
  G4LogicalVolume*   Cable_log;
  G4VPhysicalVolume* Cable_phys;
  G4LogicalVolume*   Cobre_log;
  G4VPhysicalVolume* Cobre_phys;
  G4LogicalVolume*   S_Epoxy_log;
  G4VPhysicalVolume* S_Epoxy_phys;

  G4LogicalVolume*   calibration_port_log;
  G4VPhysicalVolume* calibration_port_phys;
  G4LogicalVolume*   calibration_air_log;
  G4VPhysicalVolume* calibration_air_phys;    
  G4LogicalVolume*   calibration_Be_log;
  G4VPhysicalVolume* calibration_Be_phys;

  G4LogicalVolume*   Container_camara_log;
  G4VPhysicalVolume* Container_camara_phys;
  G4LogicalVolume*   Lente_log;
  G4VPhysicalVolume* Lente_phys;
  G4LogicalVolume*   Camara_log;
  G4VPhysicalVolume* Camara_phys;
  G4LogicalVolume*   Guide_log;
  G4VPhysicalVolume* Guide_phys;

  G4LogicalVolume*   Al_box_log;
  G4VPhysicalVolume* Al_box_phys;*/



  G4Cache<DMXScintSD*> LXeSD; //pointer to sensitive detectors
  G4Cache<DMXPmtSD*> pmtSD;


  // pointer to the Detector Messenger:
  DMXDetectorMessenger*  detectorMessenger;

};

#endif

