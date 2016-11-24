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
// $Id: ExN01RunAction.cc,v 1.19 2008/01/17 17:31:32 maire Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN01RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::ExN01RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01RunAction::~ExN01RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
  //initialize cumulative quantities
  //
  sumElAr_phys = sum2ElAr_phys =sumEPCGe_Ge_1_top_phys = sum2EPCGe_Ge_1_top_phys = 0.;
  sumLlAr_phys = sum2LlAr_phys =sumLPCGe_Ge_1_top_phys = sum2LPCGe_Ge_1_top_phys = 0.; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::fillPerEvent(G4double ElAr_phys, G4double EPCGe_Ge_1_top_phys,
                                  G4double LlAr_phys, G4double LPCGe_Ge_1_top_phys)
{
  //accumulate statistic
  //
  sumElAr_phys += ElAr_phys;  sum2ElAr_phys += ElAr_phys*ElAr_phys;
  sumEPCGe_Ge_1_top_phys += EPCGe_Ge_1_top_phys;  sum2EPCGe_Ge_1_top_phys += EPCGe_Ge_1_top_phys*EPCGe_Ge_1_top_phys;
  
  sumLlAr_phys += LlAr_phys;  sum2LlAr_phys += LlAr_phys*LlAr_phys;
  sumLPCGe_Ge_1_top_phys += LPCGe_Ge_1_top_phys;  sum2LPCGe_Ge_1_top_phys += LPCGe_Ge_1_top_phys*LPCGe_Ge_1_top_phys;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //compute statistics: mean and rms
  //
  sumElAr_phys /= NbOfEvents; sum2ElAr_phys /= NbOfEvents;
  G4double rmsElAr_phys = sum2ElAr_phys - sumElAr_phys*sumElAr_phys;
  if (rmsElAr_phys >0.) rmsElAr_phys = std::sqrt(rmsElAr_phys); else rmsElAr_phys = 0.;
  
  sumEPCGe_Ge_1_top_phys /= NbOfEvents; sum2EPCGe_Ge_1_top_phys /= NbOfEvents;
  G4double rmsEPCGe_Ge_1_top_phys = sum2EPCGe_Ge_1_top_phys - sumEPCGe_Ge_1_top_phys*sumEPCGe_Ge_1_top_phys;
  if (rmsEPCGe_Ge_1_top_phys >0.) rmsEPCGe_Ge_1_top_phys = std::sqrt(rmsEPCGe_Ge_1_top_phys); else rmsEPCGe_Ge_1_top_phys = 0.;
  
  sumLlAr_phys /= NbOfEvents; sum2LlAr_phys /= NbOfEvents;
  G4double rmsLlAr_phys = sum2LlAr_phys - sumLlAr_phys*sumLlAr_phys;
  if (rmsLlAr_phys >0.) rmsLlAr_phys = std::sqrt(rmsLlAr_phys); else rmsLlAr_phys = 0.;
  
  sumLPCGe_Ge_1_top_phys /= NbOfEvents; sum2LPCGe_Ge_1_top_phys /= NbOfEvents;
  G4double rmsLPCGe_Ge_1_top_phys = sum2LPCGe_Ge_1_top_phys - sumLPCGe_Ge_1_top_phys*sumLPCGe_Ge_1_top_phys;
  if (rmsLPCGe_Ge_1_top_phys >0.) rmsLPCGe_Ge_1_top_phys = std::sqrt(rmsLPCGe_Ge_1_top_phys); else rmsLPCGe_Ge_1_top_phys = 0.;
  
  //print
  //
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << "\n mean Energy in lAr_phys :  " << G4BestUnit(sumElAr_phys,"Energy")
     << " +- "                          << G4BestUnit(rmsElAr_phys,"Energy")  
     << "\n mean Energy in PCGe_Ge_1_top_phys : " << G4BestUnit(sumEPCGe_Ge_1_top_phys,"Energy")
     << " +- "                          << G4BestUnit(rmsEPCGe_Ge_1_top_phys,"Energy")
     << G4endl;
     
  G4cout
     << "\n mean trackLength in lAr_phys :  " << G4BestUnit(sumLlAr_phys,"Length")
     << " +- "                               << G4BestUnit(rmsLlAr_phys,"Length")  
     << "\n mean trackLength in PCGe_Ge_1_top_phys : " << G4BestUnit(sumLPCGe_Ge_1_top_phys,"Length")
     << " +- "                               << G4BestUnit(rmsLPCGe_Ge_1_top_phys,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
