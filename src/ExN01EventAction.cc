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
// $Id: ExN01EventAction.cc,v 1.29 2008/01/17 17:31:32 maire Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN01EventAction.hh"

#include "ExN01RunAction.hh"

#include "ExN01MultichannelTally.hh"

#include "G4Event.hh"
#include "G4VVisManager.hh"
#include "G4SystemOfUnits.hh"

G4double TriggerSignalThreshold[30]={500*eV,
			500*eV,500*eV,500*eV,500*eV,500*eV,500*eV,500*eV,
			500*eV,500*eV,500*eV,500*eV,500*eV,500*eV,500*eV,
			500*eV,500*eV,500*eV,500*eV,500*eV,500*eV,500*eV,
			0,0,0,0,0,0,0,0};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01EventAction::ExN01EventAction
	(ExN01RunAction* run,ExN01PrimaryGeneratorAction* pri,ExN01MultichannelTally* mul,G4int sav)
:runAct(run),primaryGenerator(pri),multiChannel(mul),printModulo(100000),saveModulo(sav)
{
	NumberOfSourceNuclidesDecay = 0;
	memset(EnergyDepositionInDetector,0,30*sizeof(G4double));
	memset(TrackLengthInDetector,0,30*sizeof(G4double));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01EventAction::~ExN01EventAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  //if (evtNb>0 && evtNb%1000000 == 0)//保存阶段性结果
	//{ 
		//G4cout<<"region of particle is "<<primaryGenerator->getR_Min()-50<<" - "<<primaryGenerator->getR_Max()-50<<G4endl;
		//multiChannel->OutputFile("f",false);
		//multiChannel->ClearAllData();
	//}
  if (evtNb>0 && evtNb%saveModulo == 0)//保存阶段性结果
  { 
	  //multiChannel->SetMessage(inputset); //向文件中加入一个说明信息
	  //multiChannel->SetNPS(evtNb);
	  multiChannel->SetNPS(NumberOfSourceNuclidesDecay);
	  multiChannel->SetFactor(1.); //将数据除以一个因子，比如面积、体积之类。
	  multiChannel->OutputFile(tmp_outputname,true); //将数据保存到指定的文件中去。


	  G4cout << "<--- Phased result of " << evtNb << " events saved.\n" << G4endl;
  }
  if (evtNb%printModulo == 0)
  { 
	  G4cout << "---> Begin of event: " << evtNb << G4endl;
	  //CLHEP::HepRandom::showEngineStatus();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01EventAction::EndOfEventAction(const G4Event* evt)
{
  //accumulates statistic
  G4int evtNb = evt->GetEventID();
  //
  //runAct->fillPerEvent(EnergyDepositionInDetector[0], EnergyDepositionInDetector[1],
		//			   TrackLengthInDetector[0], TrackLengthInDetector[1]);
  //if (evtNb>0 && evtNb%saveModulo == 0)//保存阶段性结果
  if(primaryGenerator->getEnergyLevel()==0)
	  NumberOfSourceNuclidesDecay++;//衰变次数加1

  if(primaryGenerator->getIfEndThisCount())
  {
	  G4double positionweight=primaryGenerator->getPositionWeight();
	  G4double directionweight=primaryGenerator->getDirectionWeight();
	  G4double activityfactor=primaryGenerator->getActivityFactor();
	  G4double WeightToAddData = positionweight * directionweight * activityfactor;

	  //if(primaryGenerator->getIfAlpha())
	  //{
		//  G4double energyadd = primaryGenerator->getEnergy();
		//  multiChannel->AddData(0,0,energyadd,WeightToAddData);
	  //}
	  //else
	  //{
		//  multiChannel->AddData(0,0,EnergyDepositionInDetector[0],WeightToAddData);
	  //}
	  //G4bool coutalpha = true;
	  if(primaryGenerator->getIfAlpha()==1)
	  {
		  G4double alphaene = primaryGenerator->getAlphaEnergy();
		  G4double alphaweight = primaryGenerator->getAlphaFactor();
		  alphaweight = 0.01430*alphaene+0.23751084;
		  addEnergyAndTrackInDetector(0,alphaene*alphaweight,0);
		  primaryGenerator->SetIfAlpha(false);
		  //countalpha = false;
	  }
	  if (EnergyDepositionInDetector[0]>0 && EnergyDepositionInDetector[1]==0)
	  {
		  multiChannel->AddData(static_cast<int>(evtNb / saveModulo), EnergyDepositionInDetector[0], 0);
	  }
// 	  for (G4int i=0;i<=1;i++)
// 	  {
// 		multiChannel->AddData(static_cast<int>(evtNb/saveModulo),EnergyDepositionInDetector[i],i);
// 	  }
// 	  if (1000-static_cast<int>(1000-EnergyDepositionInDetector[0]/keV)==662)
// 	  {
// 		  multiChannel->SetNumOfParticle(static_cast<int>(evtNb / saveModulo), 0);
// 	  }
	  if (EnergyDepositionInDetector[0] > 0 && EnergyDepositionInDetector[1] == 0)
	  {
		  multiChannel->SetNumOfParticle(static_cast<int>(evtNb / saveModulo), 0);
	  }
	  if (EnergyDepositionInDetector[0] > 0 && EnergyDepositionInDetector[1] > 0)
	  {
		  multiChannel->SetNumOfParticle(static_cast<int>(evtNb / saveModulo), 1);
	  }
// 	  if (EnergyDepositionInDetector[0] == 0 && EnergyDepositionInDetector[1] > 0)
// 	  {
// 		  multiChannel->SetNumOfParticle(static_cast<int>(evtNb / saveModulo), 3);
// 	  }
// 	  if (1000 - static_cast<int>(1000 - EnergyDepositionInDetector[1] / keV) == 662)
// 	  {
// 		  multiChannel->SetNumOfParticle(static_cast<int>(evtNb / saveModulo), 4);
// 	  }


	  //multiChannel->UpdateThisEvent();
	   
	  // 为下一个计数初始化
	  memset(EnergyDepositionInDetector,0,30*sizeof(G4double));
	  memset(TrackLengthInDetector,0,30*sizeof(G4double));

  }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
