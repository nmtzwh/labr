#ifndef ExN01MultichannelTally_hh
#define ExN01MultichannelTally_hh
 
#include "globals.hh"
 
class ExN01MultichannelTally
{
	public:
		ExN01MultichannelTally(G4double e1,G4double e2, G4int numofchannel,
							   G4int numofevent,G4int numofdetector,G4int);
		//ExN01MultichannelTally(G4int numberofevent);
		~ExN01MultichannelTally();
		void ClearAllData();
		void AddData(G4int detectorNumber,G4double Energy,G4int flag);
		//void AddData2(G4int detectorNumber,G4double Energy,G4double Data);
		//void AddData(G4double Data);
		void OutputFile(G4String FileName,G4bool ifCoverFile);
		void SetMessage(G4String msg){message=msg;}
		void SetFactor(G4double f){factor=f;}
		void SetNPS(G4int nps){NPS=nps;}
		//void UpdateThisEvent();
		void SetEnergyOfParticle(G4double energy) { EnergyOfParticle = energy; }
		void SetNumOfParticle(G4int position, G4int flag) { NumofParticle[position][flag] += 1; }

	private:
		//void CalErr();
		G4double LaBr3Data[20][10000], LaBr3TData[20][10000];
		G4double NaIData[20][10000], NaITData[20][10000];
		G4double E1,E2,dE;
		G4int NCHN,NPS,saveModulo;
		G4int detectorNumberTotal;
		G4String message;  //additional message
		G4double factor;     //constant factor like area or volume
		G4double EnergyOfParticle;
		G4int NumofParticle[20][5];//不同行代表不同位置的粒子；列依次代表：在LaBr3中沉积全部能量,仅在LaBr3中沉积能量,在LaBr3和NaI中沉积能量,仅在NaI中沉积能量,在NaI中沉积全部能量
 
};
#endif
