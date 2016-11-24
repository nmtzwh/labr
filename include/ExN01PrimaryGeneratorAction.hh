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
// $Id: ExN01PrimaryGeneratorAction.hh,v 1.5 2006/06/29 17:47:17 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#ifndef ExN01PrimaryGeneratorAction_h
#define ExN01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

class ExN01DetectorConstruction;
class ExN01EventAction;
class G4ParticleGun;
class G4Event;
class ExN01MultichannelTally;

class ExN01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
	ExN01PrimaryGeneratorAction(ExN01MultichannelTally* mul,
		ExN01DetectorConstruction* det,G4int srcregion,G4int srcparticle,
		G4double dettime,G4int,G4double,G4double);
	~ExN01PrimaryGeneratorAction();

	G4double getDirectionWeight() {return DirectionWeight;};
	G4double getPositionWeight() {return PositionWeight;};
	G4double getActivityFactor() {return ActivityFactor;};
	G4int getEnergyLevel() {return EnergyLevel;};
	G4int getIfEndThisCount() {return IfEndThisCount;};
	G4int getIfAlpha() {return IfAlpha;};
	G4double getAlphaEnergy() {return alphaEnergy;};
	G4double getAlphaFactor(){return alphafactor;};
	G4double getR_Min(){return R_Min;};
	G4double getR_Max(){return R_Max;};
	//G4int getparticle(){return particletype;};

	//void SetNPS(G4int nps){NPS=nps;}
	void SetIfAlpha(G4bool alpha) {IfAlpha = alpha;};
	void GeneratePrimaries(G4Event* anEvent);

  private:
	G4ParticleGun* particleGun;
	ExN01DetectorConstruction* detector;
	ExN01MultichannelTally* multiChannel;

	G4int sourceregion; //要抽样的源项区域；
	G4double source_rBiasPowerExponent;//径向偏倚指数；
	G4double source_zBiasPowerExponent;//竖向偏倚指数；
	G4int particlesort,saveModulo; //要抽样的源项粒子类型；
	G4double alphaEnergy;
	G4double alphafactor;

	G4double DirectionWeight; //粒子的方向偏倚权重；
	G4double PositionWeight; //粒子的位置偏倚权重；
	G4double ActivityFactor; //归一化到每活度源项的系数；

	G4int EnergyLevel; //核素跃迁前的能级；
	G4int EnergyLevel2; //核素跃迁后的能级；
	G4ThreeVector Position; //源粒子位置坐标；
	G4ThreeVector Direction; //源粒子发射方向坐标；
	G4double DetectorResponseTime; //探测器计数响应时间；
	G4bool IfEndThisCount; //探测器对本次计数是否终结；
	G4bool IfAlpha;//产生α粒子
	G4bool samplealpha;//α衰变分支
	G4bool samplebeta;//β衰变分支
	G4bool sampleec;//ec电子分支
	G4int numofparticle;
	G4double R_Min,R_Max;
	//G4int particletype;//粒子种类，用于输出文件命名
	G4double tube_rmin;
	G4double tube_rmax;                        
	G4double tube_zmin;
	G4double tube_zmax;
	const G4double* a;
	G4double sizeofcrystal,Distance;
};

#endif


