#include "ExN01MultichannelTally.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iomanip>
using namespace std;

//Single-channel constructor.
//ExN01MultichannelTally::ExN01MultichannelTally(G4int numberofevent)
//{
//	E1 = E2 = dE = 0;
//	NCHN = 1;
//	NPS = numberofevent;
//	message = "";
//	factor = 1.;
//}

//Multi-channel constructor.
ExN01MultichannelTally::ExN01MultichannelTally
		(G4double e1,G4double e2, G4int numofchannel, G4int numofevent,
		 G4int numofdetector,G4int sav)
:E1(e1),E2(e2),NCHN(numofchannel+1),NPS(numofevent),detectorNumberTotal(numofdetector),message(""),factor(1.),saveModulo(sav)
{
	dE=(E2-E1)/(NCHN-1);
}

//Deconstructor.
ExN01MultichannelTally::~ExN01MultichannelTally()
{;}

//Clear all tally data(multi-channel).
void ExN01MultichannelTally::ClearAllData()
{
	for(G4int k=0;k<=detectorNumberTotal-1;k++)
	{
		for (G4int i = 0; i < NCHN; i++)
		{
			if (LaBr3Data[k][i]>0)
				LaBr3Data[k][i] = 0;
			if (NaIData[k][i] > 0)
				NaIData[k][i] = 0;
		}
	}
}

//Record tally data(multi-channel).
void ExN01MultichannelTally::AddData(G4int detectorNumber,G4double Energy,G4int flag)
{
	if(Energy>0)
	{
		G4int channel = static_cast<int> ((E2-Energy)/dE);
		channel = NCHN-1-channel;
		//if(channel<0) channel = 0;
		if (channel > 0)
		{
			if (flag == 0)
				LaBr3Data[detectorNumber][channel] += 1;
			else if(flag == 1)
				NaIData[detectorNumber][channel] += 1;
		}
		//TData[channel] += Data;
	}
}


void ExN01MultichannelTally::OutputFile(G4String FileName,G4bool ifCoverFile)
{
   while(true)
   {
	   ifstream f;//f(FileName);
	   f.open(FileName);
	   if (f && ifCoverFile==false)
	   {
		   //G4cout<< "File \"" << FileName << "\" already exists." <<G4endl;
		   FileName = "z_"+FileName;
		   f.close();
	   }
	   else
		   break;
   }
   fstream dataFile(FileName, ios::out);// | ios::trunc);
   dataFile<<"Spectrum of different areas(cps/Bq)"<<"\t\t";
   dataFile << "Number of source particles: \t" << NPS << "\t\t";
   dataFile << "Number of source particles in every region: \t" << saveModulo << "\t\t" << G4endl;
//    dataFile << "EnergyInLaBr3>660keV\t";
//    for (G4int k = 0; k < detectorNumberTotal; k++)
// 	   dataFile << NumofParticle[k][0] << "\t\t";
//    dataFile << G4endl;
   dataFile << "EnergyOnlyInLaBr3\t";
   for (G4int k = 0; k<detectorNumberTotal; k++)
	   dataFile << NumofParticle[k][0] << "\t";
   dataFile << G4endl;
   dataFile << "EnergyInLaBr3AndNaI\t";
   for (G4int k = 0; k<detectorNumberTotal; k++)
	   dataFile << NumofParticle[k][1] << "\t";
   dataFile << G4endl;
//    dataFile << "EnergyOnlyInNaI\t";
//    for (G4int k = 0; k<detectorNumberTotal; k++)
// 	   dataFile << NumofParticle[k][3] << "\t\t";
//    dataFile << G4endl;
//    dataFile << "EnergyInNaI>660keV\t";
//    for (G4int k = 0; k<detectorNumberTotal; k++)
// 	   dataFile << NumofParticle[k][4] << "\t\t";
//    dataFile << G4endl;
   dataFile << "R_max(cm)\t";
   for (G4int k = 0; k<detectorNumberTotal; k++)
   {
	   dataFile << 50 * (k + 1) << "\t";
   }
   dataFile << endl;
   dataFile << "Region\t";
   for (G4int k = 0; k<detectorNumberTotal; k++)
   {
	   dataFile << "Region- " << k << "\t";
   }
   dataFile << G4endl;
   for(G4int i=0;i<NCHN;i++)
   {
	   dataFile << std::fixed << (E1+i*dE)/keV << " \t ";
	   for (G4int j=0;j<detectorNumberTotal;j++)
	   {
		   dataFile << std::fixed << LaBr3Data[j][i]<< " \t ";
	   }
	   dataFile<<endl;
   }
   dataFile.close();
   G4cout<< "File \"" << FileName << "\" created." << G4endl;
}
