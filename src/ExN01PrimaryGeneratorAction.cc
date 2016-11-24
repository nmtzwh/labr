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
// $Id: ExN01PrimaryGeneratorAction.cc,v 1.6 2006/06/29 17:47:23 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01DetectorConstruction.hh"
#include "ExN01MultichannelTally.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "Randomize.hh"

#include "G4TransportationManager.hh" 
#include "G4Navigator.hh" 

#include <fstream>
#include <iomanip>
using namespace std;

  //--------------- 源项抽样需调用的函数定义
  //--------------- source sample function definition

  //对一段区间内的均匀分布抽样出值：
  G4double sampleHomogeneous(G4double interval_min,G4double interval_max) 
	{
		return interval_min+(interval_max-interval_min)*G4UniformRand();
	}

  //对一段区间内的梯形分布抽样出值：（参照《蒙特卡罗方法在实验核物理中的应用》第40页）
  G4double sampleTrapezoidal(G4double x0,G4double x1,G4double y0,G4double y1) 
	{
		G4double ksi = G4UniformRand();
		if( ksi <= std::min(y0,y1)*2/(y0+y1) ) return x0+(x1-x0)*G4UniformRand();
		else if(y0<=y1) return x0+(x1-x0)*std::max(G4UniformRand(),G4UniformRand());
			 else		return x0+(x1-x0)*std::min(G4UniformRand(),G4UniformRand());
	}

  //对[0,+∞)区间内以λ为参数（数学期望为1/λ）的指数分布抽样出值：（直接抽样方法）
  G4double sampleExponential(G4double lamda) 
	{
		return -std::log(G4UniformRand())/lamda;
	}

  //对在0~2π均匀分布的角，抽样出一对cos值和sin值：（参照课件3《由已知分布的随机抽样》第37张幻灯片）
  G4ThreeVector sampleCosAndSin()
	{
		G4double ksi1,ksi2;
		G4double sample_cos,sample_sin;
		do
		{
			ksi1 = G4UniformRand();
			ksi2 = G4UniformRand();
			if (3*ksi1>2-ksi2) {ksi1=1-ksi1; ksi2=ksi2-1;}
		}while(3*ksi1*ksi1+ksi2*ksi2>1);
		sample_cos=(3*ksi1*ksi1-ksi2*ksi2)/(3*ksi1*ksi1+ksi2*ksi2);
		sample_sin=(3.4641016*ksi1*ksi2)/(3*ksi1*ksi1+ksi2*ksi2);
		return G4ThreeVector(sample_cos,sample_sin,0);
	}

  //对平面圆环内的均匀分布，抽样出半径：（参照《蒙特卡罗方法在实验核物理中的应用》第41页）
  G4double sampleCircleR(G4double circle_rmin,G4double circle_rmax)
	{
		G4double r_;
		G4double P1=2*circle_rmin/(circle_rmin+circle_rmax);
		if (G4UniformRand()<=P1) r_=G4UniformRand();
		else r_=std::max(G4UniformRand(),G4UniformRand());
		return circle_rmin+(circle_rmax-circle_rmin)*r_;
	}

  //对平面圆环内的均匀分布，抽样出坐标(x,y)：
  G4ThreeVector sampleCircleXY(G4double circle_rmin,G4double circle_rmax)
	{
		G4double sample_r=sampleCircleR(circle_rmin,circle_rmax);
		G4ThreeVector sample_cos_sin=sampleCosAndSin();
		G4double sample_x = sample_r*sample_cos_sin.getX();
		G4double sample_y = sample_r*sample_cos_sin.getY();

		return G4ThreeVector(sample_x,sample_y,0);
	}

  //抽样出长方体内的均匀位置分布坐标(x,y,z)：
  G4ThreeVector samplePositionInBox
		(G4double box_xmin,G4double box_xmax,G4double box_ymin,G4double box_ymax,G4double box_zmin,G4double box_zmax)
	{
		G4double sample_x = sampleHomogeneous(box_xmin,box_xmax);
		G4double sample_y = sampleHomogeneous(box_ymin,box_ymax);
		G4double sample_z = sampleHomogeneous(box_zmin,box_zmax);

		return G4ThreeVector(sample_x,sample_y,sample_z);
	}

  //抽样出圆管体内的均匀位置分布坐标(x,y,z)：
  G4ThreeVector samplePositionInTube
		(G4double tube_rmin,G4double tube_rmax,G4double tube_zmin,G4double tube_zmax)
	{
		G4ThreeVector sample_cos_sin=sampleCosAndSin();

		G4double sample_r = sampleCircleR(tube_rmin,tube_rmax);
		G4double sample_x = sample_r*sample_cos_sin.getX();
		G4double sample_y = sample_r*sample_cos_sin.getY();
		G4double sample_z = sampleHomogeneous(tube_zmin,tube_zmax);

		return G4ThreeVector(sample_x,sample_y,sample_z);
	}

  //对球壳内的均匀分布，抽样出半径：（参照《蒙特卡罗方法在实验核物理中的应用》第42页）
  G4double sampleSphereR(G4double sphere_rmin,G4double sphere_rmax) 
	{
		G4double r_;
		G4double lamda=sphere_rmin*sphere_rmin+sphere_rmin*sphere_rmax+sphere_rmax*sphere_rmax;
		G4double ksi1=G4UniformRand();
		if (ksi1<=3*sphere_rmin*sphere_rmin/lamda) r_=G4UniformRand();
		else{if (ksi1<=3*sphere_rmin*sphere_rmax/lamda) r_=std::max(G4UniformRand(),G4UniformRand());
			 else r_=std::max(std::max(G4UniformRand(),G4UniformRand()),G4UniformRand());
			}
		return sphere_rmin+(sphere_rmax-sphere_rmin)*r_;
	}
  
  //抽样出球体(半径为r_max-r_min的球壳,球心位于(x,y,z))内的均匀位置分布坐标(x,y,z)：
  G4ThreeVector samplePositionInSphere
		(G4double r_min,G4double r_max,G4double x,G4double y,G4double z)
	{
		G4double sample_r = sampleSphereR(r_min,r_max);//半径
		G4ThreeVector sample_cos_sin=sampleCosAndSin();
		G4double costheta = sampleHomogeneous(-1.0,1.0);
		G4double sintheta = sqrt(1-costheta*costheta);

		G4double sample_x = sample_r*sample_cos_sin.getX()*sintheta;
		G4double sample_y = sample_r*sample_cos_sin.getY()*sintheta;
		G4double sample_z = sample_r*costheta;

		return G4ThreeVector(sample_x,sample_y,sample_z);
	}

  //抽样源粒子发射方向,θ范围在theta_min~theta_max之间,φ角2π均匀分布
  G4ThreeVector sampleDirection(G4double theta_max, G4double theta_min)
  {
	  G4double costheta = (theta_max - theta_min)*G4UniformRand() + theta_min;
	  G4double sintheta = sqrt(1 - costheta*costheta);
	  G4ThreeVector phi = sampleCosAndSin();
	  return G4ThreeVector(sintheta*phi.getY(), sintheta*phi.getX(), costheta);
  }

  //从用区间上界表示的累积概率表中抽样出值：
  G4double sampleFromStaircaseDistributionTable(G4double* binUpperBound, G4double* cumProb) 
	{
		G4int i=0;
		G4double ksi=G4UniformRand();
		while (*(cumProb+i)<ksi) i++;
		return sampleHomogeneous(*(binUpperBound+i-1),*(binUpperBound+i));
	}

  //从描点表示的折线形概率密度函数中抽样出值：
  G4double sampleFromLineDistributionTable(G4double* tracingPoint, G4double* probDen, G4double* cumProb) 
	{
		G4int i=0;
		G4double ksi=G4UniformRand();
		while (*(cumProb+i)<ksi) i++;
		return sampleTrapezoidal(*(tracingPoint+i-1),*(tracingPoint+i),*(probDen+i-1),*(probDen+i));
	}

  //从离散分布累积概率表中抽样出值：
  G4double sampleFromDiscreteDistributionTable(G4double* discreteValue, G4double* cumProb) 
	{
		G4int i=0;
		G4double ksi=G4UniformRand();
		while (*(cumProb+i)<ksi) i++;
		return *(discreteValue+i);
	}

  //从离散分布累积概率表中抽样出序号：
  G4int sampleRank(G4double* cumProb) 
	{
		G4int i=0;
		G4double ksi=G4UniformRand();
		while (*(cumProb+i)<ksi) i++;
		return i;
	}

  //从离散分布累积概率表中抽样出序号和该序号对应的概率值：
  G4ThreeVector sampleRankAndProbability(G4double* cumProb) 
	{
		G4int i=0;
		G4double ksi=G4UniformRand();
		while (*(cumProb+i)<ksi) i++;
		if(i==0)
			return G4ThreeVector(i,*(cumProb+i),0);
		else
			return G4ThreeVector(i,*(cumProb+i)-*(cumProb+i-1),0);
	}

  //抽样出各向同性的方向向量：（参照《蒙特卡罗方法在实验核物理中的应用》第61页）
  G4ThreeVector sampleIsotropicVector()
	{
		G4double sample_costhita,sample_sinthita; //(cosθ)在[-1,1]上均匀分布
		sample_costhita = 2*G4UniformRand()-1;
		sample_sinthita = sqrt(1-sample_costhita*sample_costhita);

		G4double sample_cosphi,sample_sinphi;	  //φ在[0,2π]上均匀分布
		G4ThreeVector sample_cosphi_sinphi=sampleCosAndSin();
		sample_cosphi = sample_cosphi_sinphi.getX();
		sample_sinphi = sample_cosphi_sinphi.getY();

		G4double vector_u,vector_v,vector_w;
		vector_u = sample_sinthita*sample_cosphi; //u=sinθ·cosφ
		vector_v = sample_sinthita*sample_sinphi; //v=sinθ·sinφ
		vector_w = sample_costhita;				  //w=cosθ

		return G4ThreeVector(vector_u,vector_v,vector_w);
	}

  //从用区间上界表示的累积概率表、偏倚概率表中抽样出值和偏倚权重：
  G4ThreeVector sampleBiasFromStaircaseDistributionTable
		(G4double* binUpperBound, G4double* cumProb,  G4double* biasCumProb) 
	{
		G4int i=0;
		G4double samplevalue,biasweight=1.;
		G4double ksi=G4UniformRand();
		while (*(biasCumProb+i)<ksi) i++;
		samplevalue = sampleHomogeneous(*(binUpperBound+i-1),*(binUpperBound+i));
		biasweight = ( *(cumProb+i)-*(cumProb+i-1) ) /
				 ( *(biasCumProb+i)-*(biasCumProb+i-1) );
		return G4ThreeVector(samplevalue,biasweight,0);
	}


  G4bool sampleIfEndThisCount(G4double halflife, G4double dettime)
	{
		G4double LevelTransitionTime = sampleExponential(0.69314718/halflife);
		if(LevelTransitionTime > dettime) return true;
		else return false;
	}

  //--------------- 本段只在程序开始时运行一次，用来写不需抽样的源项参数（Geant4固有段）

  ExN01PrimaryGeneratorAction::ExN01PrimaryGeneratorAction
	  (ExN01MultichannelTally* mul,ExN01DetectorConstruction* det,G4int srcregion,G4int srcparticle,G4double dettime,G4int sav,G4double size,G4double dis)
	  :multiChannel(mul),detector(det),sourceregion(srcregion),particlesort(srcparticle),DetectorResponseTime(dettime),saveModulo(sav),sizeofcrystal(size),Distance(dis)
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  if(sourceregion==-1)
  {
	  G4cout<<"Please input source region (4):"<<G4endl;
	  G4cin>>sourceregion;
  }

  if(particlesort==-1)
  {
	  G4cout<<"Please input particle sort (78):"<<G4endl;
	  G4cin>>particlesort;
	  //G4cin>>particletype;
  }

  EnergyLevel=0; //能级初始化
  EnergyLevel2=0;
  samplealpha=false;
  samplebeta=false;
  sampleec=false;
  IfAlpha=false;
  numofparticle = 0;
  R_Min = 0;
  R_Max = 0;
}

  ExN01PrimaryGeneratorAction::~ExN01PrimaryGeneratorAction()
{
  delete particleGun;
}

  //--------------- 本段在每次产生新源粒子时运行一次，用来写需要抽样的源项参数（Geant4固有段）

  void ExN01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

 //--------------- lAr source Position sampling start ---------------
  PositionWeight = 1.;

  switch(sourceregion)
  {
  case 1:
	  {
		  a = detector -> GetSample_Info();      
		  tube_rmin = a[0];
		  tube_rmax = a[1];                                 
		  tube_zmin = a[5]-a[2];
		  tube_zmax = a[5]+a[2];
		  G4VPhysicalVolume* LocateVolume;

		  Position = samplePositionInTube(tube_rmin,tube_rmax,tube_zmin,tube_zmax);
		  Position = G4ThreeVector((Position.getX())+a[3],(Position.getY())+a[4],(Position.getZ()));
		  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
		  particleGun->SetParticlePosition(Position);


		  LocateVolume = theNavigator->LocateGlobalPointAndSetup(Position);
		  if(LocateVolume != detector->GetLaBr3Crystal_phys()){
			  G4int aa = 3;
		  }
		  Direction = sampleIsotropicVector();
		  DirectionWeight = 1.;
	  }break;

  case 2:
	  {
		  G4double Labr_rmin = 0.*mm,Labr_rmax = 76.2/2*mm;
		  G4double Labr_zmin = -76.2/2*mm,Labr_zmax = 76.2/2*mm;
		  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
		  G4VPhysicalVolume* LocateVolume;

		  if(EnergyLevel==0)
		  {
			  //do
			  //{
				  Position = samplePositionInTube(Labr_rmin,Labr_rmax,Labr_zmin,Labr_zmax);
				  LocateVolume = theNavigator->LocateGlobalPointAndSetup(Position);
			  //}
			  //while(LocateVolume != detector->GetDetector_phys());
		  }
		  PositionWeight = 1.;
		  Direction = sampleIsotropicVector();
		  DirectionWeight = 1.;
	  }break;

  case 3:
	  {
		  G4double Labr_rmin = 0.*mm,Labr_rmax = 38.1/2*mm;
		  G4double Labr_zmin = -38.1/2*mm,Labr_zmax = 38.1/2*mm;
		  G4double Position_x,Position_y,Position_z;

		  if(EnergyLevel==0)
		  {
			  G4ThreeVector pos = samplePositionInTube(Labr_rmin,Labr_rmax,Labr_zmin,Labr_zmax);
			  Position_x = pos.getX();
			Position_y = pos.getY();
			Position_z = pos.getZ() + 79.5*mm;
			Position = G4ThreeVector(Position_x,Position_y,Position_z);
		  }
			  //Position = samplePositionInTube(Labr_rmin,Labr_rmax,Labr_zmin,Labr_zmax);

		  PositionWeight = 1.;
		  Direction = sampleIsotropicVector();
		  DirectionWeight = 1.;
	  }break;

  case 4:
  {
	  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	  G4VPhysicalVolume* LocateVolume;
	  R_Min = 50.*(static_cast<int>(numofparticle / saveModulo))*mm;
	  R_Max = 50.*(static_cast<int>(numofparticle / saveModulo) + 1)*mm;
	  if (EnergyLevel == 0)
	  {
		  do
		  {
			  Position = samplePositionInSphere(R_Min, R_Max, 0.*mm, 0.*mm, 0.*mm);//10cm~20cm抽样
			  LocateVolume = theNavigator->LocateGlobalPointAndSetup(Position);
		  } while ((LocateVolume != detector->GetWater_phys()));
	  }
	  PositionWeight = 1.;
	  numofparticle++;
	  Direction = sampleIsotropicVector();
	  DirectionWeight = 1.;
  }break;

  case 5:
  {
	  G4double inch = 25.4*mm;
	  Position = G4ThreeVector(0, 0, -Distance*mm);
	  PositionWeight = 1.;
	  Direction = sampleDirection(1, Distance / (sqrt(Distance*Distance + sizeofcrystal*25.4 / 2 * sizeofcrystal*25.4 / 2)));
	  DirectionWeight = 1.;
  }break;

	  

  }//end of switch
  particleGun->SetParticlePosition(Position);
 //--------------- lAr source Position sampling end ---------------

 //--------------- source direction sampling start ---------------
  DirectionWeight = 1.;

  particleGun->SetParticleMomentumDirection(Direction);
 //--------------- source direction sampling end ---------------

 //--------------- source energy sampling start ---------------
  //ActivityFactor = 1.;
  G4double Energy=0;
  if(EnergyLevel==0) IfEndThisCount=false;  //计数终结标志初始化
  particleGun->SetParticleDefinition(particleTable->
	  FindParticle(particleName="gamma"));
  ActivityFactor = 1.0;		
  EnergyLevel2=0;

  switch(particlesort)
  {
  case 1:  Energy = 657.750*keV;  break;//Ag110  
  case 2:  Energy = 884.670*keV;  break;//Ag110  
  case 3:  Energy = 810.766*keV;  break;//Co58  
  case 4:  Energy = 1173.240*keV; break;//Co60  
  case 5:  Energy = 1332.508*keV; break;//Co60  
  case 6:  Energy = 569.331*keV;  break;//Cs134  
  case 7:  Energy = 604.722*keV;  break;//Cs134  
  case 8:  Energy = 795.868*keV;  break;//Cs134  
  case 9:  Energy = 661.659*keV;  break;//Cs137  
  case 10: Energy = 364.490*keV;  break;//I131  
  case 11: Energy = 834.855*keV;  break;//Mn54  
  case 12: Energy = 602.728*keV;  break;//Sb124  
  case 13: Energy = 1690.984*keV; break;//Sb124
  case 14: Energy = 63.30*keV;    break;//Th-234 3.75%
  case 15: Energy = 92.38*keV;	  break;//Th-234 2.18%
  case 16: Energy = 92.80*keV;	  break;//Th-234 2.15%
  case 17: Energy = 1460.822*keV;	  break;//K-40 10.55%

  case 18: Energy = 50.*keV;	  break;
  case 19: Energy = 100.*keV;	  break;
  case 20: Energy = 150.*keV;	  break;
  case 21: Energy = 200.*keV;	  break;
  case 22: Energy = 250.*keV;	  break;
  case 23: Energy = 300.*keV;	  break;
  case 24: Energy = 350.*keV;	  break;
  case 25: Energy = 400.*keV;     break;
  case 26: Energy = 500.*keV;     break;
  case 27: Energy = 600.*keV;     break;
  case 28: Energy = 700.*keV;     break;
  case 29: Energy = 800.*keV;     break;
  case 30: Energy = 900.*keV;     break;
  case 31: Energy = 1000.*keV;     break;
  case 32: Energy = 1100.*keV;     break;
  case 33: Energy = 1200.*keV;     break;
  case 34: Energy = 1300.*keV;     break;
  case 35: Energy = 1400.*keV;     break;
  case 36: Energy = 1500.*keV;     break;
  case 37: Energy = 1600.*keV;     break;
  }//end of switch

  particleGun->SetParticleEnergy(Energy);
  multiChannel->SetEnergyOfParticle(Energy);
  //particleEnergy = Energy;
  EnergyLevel = EnergyLevel2;
  if(EnergyLevel==0) IfEndThisCount = true;

 //--------------- source energy sampling end ---------------

  particleGun->GeneratePrimaryVertex(anEvent);
}
