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
// $Id: ExN01EventAction.hh,v 1.12 2007/07/02 13:22:08 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN01EventAction_h
#define ExN01EventAction_h 1

#include "G4UserEventAction.hh"
#include "ExN01PrimaryGeneratorAction.hh"

class ExN01PrimaryGeneratorAction;
class ExN01RunAction;
class ExN01MultichannelTally;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN01EventAction : public G4UserEventAction
{
public:
  ExN01EventAction(ExN01RunAction*,ExN01PrimaryGeneratorAction*,ExN01MultichannelTally*,G4int);
  virtual ~ExN01EventAction();

  void  BeginOfEventAction(const G4Event*);
  void    EndOfEventAction(const G4Event*);
    
  void addEnergyAndTrackInDetector(G4int detectorNumber,G4double de, G4double dl) 
	  {
		  EnergyDepositionInDetector[detectorNumber] += de;
		  TrackLengthInDetector[detectorNumber] += dl;
	  };

  void SetPrintModulo(G4int val)  {printModulo = val;};
    
	void Settmp_outputname(G4String tmp){tmp_outputname=tmp;}

  G4int getNumberOfSourceNuclidesDecay() {return NumberOfSourceNuclidesDecay;};

private:
   
   ExN01EventAction* eventaction;
   ExN01RunAction*  runAct;
   ExN01PrimaryGeneratorAction* primaryGenerator;
   ExN01MultichannelTally* multiChannel;
                             
   G4int  printModulo,saveModulo;
   G4int  NumberOfSourceNuclidesDecay; //这是用于归一化到每个源核素贡献的系数，与event数目不同
   G4String tmp_outputname;
   G4double  EnergyDepositionInDetector[30];
   G4double  TrackLengthInDetector[30];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

    
