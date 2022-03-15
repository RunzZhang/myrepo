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
// PmtSD (sensitive PMT) program
// --------------------------------------------------------------

#include "DMXPmtSD.hh"


#include "DMXDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

//added added some package by Ryan////

#include "DMXPmtHit.hh"
#include "G4OpticalPhoton.hh"

#include <iostream>
#include <fstream>
#include <cstring>
#include "G4UnitsTable.hh"

#include "G4HCofThisEvent.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"

/////and one parameter////////
int m=1;
//////////added finished//////////////

////////////////////////////////////////////////////////////////////////////


DMXPmtSD::DMXPmtSD(G4String name) 
  :G4VSensitiveDetector(name) {

  G4String HCname="pmtCollection";
  collectionName.insert(HCname);
//////added by Ryan///////////////////

  PmtInfo.open("PmtInformation.txt");
  PmtInfo1.open("PmtGeneral.txt");
  PmtInfo << "Energia_Kinetic   " << "Energia_Deposit  " << "Particula  " << "Pos_Interaccion  "<<"Pre Momentum    "<<"Delta Momentum        "<<"Physics Volume"<< std::endl;
  PmtInfo1 << "Evento   " << "Hits_Generados  " << std::endl;

}

DMXPmtSD::~DMXPmtSD() {PmtInfo.close();  PmtInfo1.close();}

////////////finish editing//////////////////


////////////////////////////////////////////////////////////////////////////
void DMXPmtSD::Initialize(G4HCofThisEvent*) {

  pmtCollection = new DMXPmtHitsCollection
    (SensitiveDetectorName,collectionName[0]); 

  HitID = -1;


}



////////////////////////////////////////////////////////////////////////////
G4bool DMXPmtSD::ProcessHits
  (G4Step* aStep, G4TouchableHistory*){

  ////added by Ryan///////////////////

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double ek = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ParticleDefinition* particleType = aStep->GetTrack()->GetDefinition();
  G4String particleName = particleType->GetParticleName();

  G4double posx = aStep->GetPostStepPoint()->GetPosition().x();
  G4double posy = aStep->GetPostStepPoint()->GetPosition().y();
  G4double posz = aStep->GetPostStepPoint()->GetPosition().z();

  G4double momx = aStep->GetPreStepPoint()->GetMomentum().x();
  G4double momy = aStep->GetPreStepPoint()->GetMomentum().y();
  G4double momz = aStep->GetPreStepPoint()->GetMomentum().z();

  G4double dmomx = aStep->GetDeltaMomentum().x();
  G4double dmomy = aStep->GetDeltaMomentum().y();
  G4double dmomz = aStep->GetDeltaMomentum().z();


  
  G4String Volume = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();

  G4double stepl = 0.;
  if (particleType->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();

  ////////////finish editing//////////////

  // make known hit position
  DMXPmtHit* aPmtHit = new DMXPmtHit();
  aPmtHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  aPmtHit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
  aPmtHit->SetParticle(particleName);
  aPmtHit->SetParticleEnergy(aStep->GetPreStepPoint()->GetKineticEnergy() );
  HitID = pmtCollection->insert(aPmtHit);

  ////added by Ryan////////////////

  //PmtInfo << m <<"  "<<G4BestUnit(ek,"Energy")<<"  "<<G4BestUnit(edep,"Energy") <<"  " <<particleName << "  "<<G4BestUnit(posx,"Length") << G4BestUnit(posy,"Length")<< G4BestUnit(posz,"Length") << std::endl;

PmtInfo << m <<"  "<<G4BestUnit(ek,"Energy")<<"  "<<G4BestUnit(edep,"Energy") <<"  " <<particleName << "  "<<G4BestUnit(posx,"Length") << G4BestUnit(posy,"Length")<< G4BestUnit(posz,"Length") <<G4BestUnit(momx,"Energy")<<"    "<<G4BestUnit(momy,"Energy")<<"    "<<G4BestUnit(momz,"Energy")<<"    "<<G4BestUnit(dmomx,"Energy")<<"    "<<G4BestUnit(dmomy,"Energy")<<"    "<<G4BestUnit(dmomz,"Energy")<<"    "<<Volume<< std::endl;



  /////finish editing////////////////

  return true;
 
}



////////////////////////////////////////////////////////////////////////////
void DMXPmtSD::EndOfEvent(G4HCofThisEvent* HCE) {

  G4String HCname = collectionName[0];

  static G4int HCID = -1;
  if(HCID<0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
  HCE->AddHitsCollection(HCID,pmtCollection);
  
  G4int nHits = pmtCollection->entries();
  if (verboseLevel>=1) {
    G4cout << "     PMT collection: " << nHits << " hits" << G4endl;

    ///added by Ryan//////////////
    PmtInfo1 << m<<"  "<<nHits << std::endl;
    m++;
    //////finish editing////////////
    if (verboseLevel>=2)
      pmtCollection->PrintAllHits();
  }


}


////////////////////////////////////////////////////////////////////////////
void DMXPmtSD::clear()    {;}


void DMXPmtSD::DrawAll()  {;}


void DMXPmtSD::PrintAll() {;}




