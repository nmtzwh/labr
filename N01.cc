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
// $Id: exampleN01.cc,v 1.6 2006/06/29 17:47:10 gunter Exp $
// GEANT4 tag $Name: geant4-09-02-patch-01 $
//
//
// --------------------------------------------------------------
//      GEANT 4 - exampleN01
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4SystemOfUnits.hh"

#include "ExN01DetectorConstruction.hh"
#include "ExN01PhysicsList.hh"
#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01RunAction.hh"
#include "ExN01EventAction.hh"
#include "ExN01SteppingAction.hh"
#include "ExN01MultichannelTally.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#endif

int main(int argc,char* argv[])
{
  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // set mandatory user action class
  G4int numberOfEvent,numberOfSave=0,src_region=-1,src_particle=-1,numberOfChannel,myseed=0;
  G4double det_time=30*microsecond,size_LaBr3,size_NaI,height_NaI;
  //G4double src_rBias=1,src_zBias=0;

  //����������ʽ��N01 ������ �м䱣�������� Դ������ Դ�������� �廯�羧���ߴ� ���������� ̽������Ӧʱ��(��λ��s)
  if(argc>1)
  {
	  numberOfEvent=atoi(argv[1]);
	  numberOfSave=atoi(argv[2]);
	  src_region=atoi(argv[3]);
	  src_particle=atoi(argv[4]);
	  size_LaBr3 = atof(argv[5]);
	  size_NaI = atof(argv[6]);
    height_NaI = atof(argv[7]); // weihe added 2016.11.16, height of surrounding NaI
	  if (argc==8)
	  {
		  argv[8]="0";
		  argv[9]="30";
	  }
	  myseed=atoi(argv[8]);
	  if(argc==9)argv[9]="30"; //��������̽������Ӧʱ��,Ĭ��Ϊ30��s,argc==6��Ϊ����������5��
	  det_time=atof(argv[9]);
  }
  else
  {
	G4cout<<"Please input nps:"<<G4endl;
	G4cin>>numberOfEvent;
	if(numberOfEvent!=-1)
	{
		G4cout<<"Please input save period:"<<G4endl;
		G4cin>>numberOfSave;
	}
	G4cout << "Please input source region (4):" << G4endl;
	G4cin >> src_region;
	G4cout << "Please input particle sort (1~17):" << G4endl;
	G4cin >> src_particle;
	G4cout << "Please input size of LaBr3:" << G4endl;
	G4cin >> size_LaBr3;
	G4cout << "Please input size of NaI:" << G4endl;
	G4cin >> size_NaI;
  G4cout << "Please input height of NaI:" << G4endl;
  G4cin >> height_NaI;
  }
  if(numberOfSave==0) numberOfSave=INT_MAX;

  // set mandatory initialization classes
  //
  ExN01DetectorConstruction* detector = new ExN01DetectorConstruction(size_LaBr3,size_NaI,height_NaI);
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new ExN01PhysicsList;
  runManager->SetUserInitialization(physics);

  numberOfChannel=2000;
  G4double Energy1 = 0;
  G4double Energy2 = numberOfChannel*1.*keV;
  G4int numberOfDetector = 12;

  ExN01MultichannelTally* multiChannel = new ExN01MultichannelTally
	  (Energy1,Energy2,numberOfChannel,numberOfEvent,numberOfDetector,numberOfSave);

  //ExN01MultichannelTally* multiChannel2 = new ExN01MultichannelTally
	  //(Energy1,Energy2,numberOfChannel,numberOfEvent,numberOfDetector,numberOfVetoSort);

  ExN01PrimaryGeneratorAction* gen_action =
					new ExN01PrimaryGeneratorAction(multiChannel,detector,
						src_region,src_particle,det_time*microsecond,numberOfSave,size_LaBr3,200.0);
  runManager->SetUserAction(gen_action);


  ExN01RunAction* run_action = new ExN01RunAction;
  runManager->SetUserAction(run_action);
  //
  ExN01EventAction* event_action =
					new ExN01EventAction(run_action,gen_action,multiChannel,numberOfSave);
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action =
					new ExN01SteppingAction(detector, event_action);
  runManager->SetUserAction(stepping_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  //
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/run/verbose 1");

  if(numberOfEvent<=100&&numberOfEvent>0)
  {
	  UI->ApplyCommand("/event/verbose 1");
	  UI->ApplyCommand("/tracking/verbose 1");
  }
  // Start a run
  //
  if(numberOfEvent==-1)
  {
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);
	#ifdef G4VIS_USE
		  G4VisManager* visManager = new G4VisExecutive;
		  visManager->Initialize();
	#endif
	#ifdef G4VIS_USE
		  UI->ApplyCommand("/control/execute vis.mac");
	#endif
		  ui->SessionStart();
		  delete ui;
	#ifdef G4VIS_USE
		  delete visManager;
	#endif
  }
  else
  {
	  char inputset[80]={" "};
	  if(argc>1)
	  {
		  for(int i=1;i<=5;i++)
		  {
			strcat(inputset,argv[i]);
			strcat(inputset," ");
		  }
	  }

	  //char outputname[80]={"CDEX-1 "};
	  char outputname[80]={""};
	  char LaBr3Size[10], NaISize[10], NaIHeight[10];
	  sprintf(LaBr3Size, "%.2f", size_LaBr3);
	  sprintf(NaISize, "%.2f", size_NaI);
    sprintf(NaIHeight, "%.2f", height_NaI);
	  strcat(outputname, LaBr3Size);
	  strcat(outputname, "LaBr3+");
	  strcat(outputname, NaISize);
    strcat(outputname, "_h");
    strcat(outputname, NaIHeight);
	  strcat(outputname, "NaI_");
	  //src_particle = gen_action->getparticle();
	  switch(src_particle)
		  {
			case 1: strcat(outputname,"110Ag1"); break;
			case 2: strcat(outputname,"110Ag2"); break;
			case 3: strcat(outputname,"58Co"); break;
			case 4: strcat(outputname,"60Co1"); break;
			case 5: strcat(outputname,"60Co2"); break;
			case 6: strcat(outputname,"134Cs1"); break;
			case 7: strcat(outputname,"134Cs2"); break;
			case 8: strcat(outputname,"134Cs3"); break;
			case 9: strcat(outputname,"137Cs"); break;
			case 10: strcat(outputname,"131I"); break;
			case 11: strcat(outputname,"54Mn"); break;
			case 12: strcat(outputname,"124Sb1"); break;
			case 13: strcat(outputname,"124Sb2"); break;
			case 14: strcat(outputname,"234Th1"); break;
			case 15: strcat(outputname,"234Th2"); break;
			case 16: strcat(outputname, "234Th3"); break;
			case 17: strcat(outputname, "40K"); break;

			case 18: strcat(outputname, "50"); break;
			case 19: strcat(outputname, "100"); break;
			case 20: strcat(outputname, "150"); break;
			case 21: strcat(outputname, "200"); break;
			case 22: strcat(outputname, "250"); break;
			case 23: strcat(outputname, "300"); break;
			case 24: strcat(outputname, "350"); break;
			case 25: strcat(outputname, "400"); break;
			case 26: strcat(outputname, "500"); break;
			case 27: strcat(outputname, "600"); break;
			case 28: strcat(outputname, "700"); break;
			case 29: strcat(outputname, "800"); break;
			case 30: strcat(outputname, "900"); break;
			case 31: strcat(outputname, "1000"); break;
			case 32: strcat(outputname, "1100"); break;
			case 33: strcat(outputname, "1200"); break;
			case 34: strcat(outputname, "1300"); break;
			case 35: strcat(outputname, "1400"); break;
			case 36: strcat(outputname, "1500"); break;
			case 37: strcat(outputname, "1600"); break;
		  }

	  switch(src_region)
		  {
			case 4: strcat(outputname,"_total"); break;
			case 5: strcat(outputname, "_200mm"); break;
		  }

	  //strcat(outputname, "  50eV.xls");
	  strcat(outputname, ".xls");

	  char tmp_outputname[80]={"tmp_"};
	  strcat(tmp_outputname,outputname);
	  event_action->Settmp_outputname(tmp_outputname);

	  runManager->BeamOn(numberOfEvent);

	  multiChannel->SetMessage(inputset); //���ļ��м���һ��˵����Ϣ
	  //multiChannel->SetNPS(numberOfEvent);
	  multiChannel->SetNPS(event_action->getNumberOfSourceNuclidesDecay());
	  G4cout<<event_action->getNumberOfSourceNuclidesDecay()<<G4endl;
	  multiChannel->SetFactor(1.); //�����ݳ���һ�����ӣ���������������֮�ࡣ
	  multiChannel->OutputFile(outputname,false); //�����ݱ��浽ָ�����ļ���ȥ��

	  remove(tmp_outputname);
  }

  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //

  delete runManager;
  delete multiChannel;
  //G4cin>>numberOfEvent;
  return 0;
}
