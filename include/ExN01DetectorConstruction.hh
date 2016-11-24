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
// $Id: ExN01DetectorConstruction.hh,v 1.6 2006/06/29 17:47:13 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#ifndef ExN01DetectorConstruction_H
#define ExN01DetectorConstruction_H 1

class G4Box;
class G4Tubs;
class G4Cons;
class G4Sphere;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UnionSolid;

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"

class ExN01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

	ExN01DetectorConstruction(G4double,G4double,G4double);
	~ExN01DetectorConstruction();

	G4VPhysicalVolume* Construct();
	const G4VPhysicalVolume* GetLaBr3Crystal_phys() { return LaBr3crystal_phys; };
	const G4VPhysicalVolume* GetNaICrystal_phys() { return NaIcrystal_phys; };
	const G4VPhysicalVolume* GetWater_phys() {return Water_phys;};
	const G4double* GetSample_Info()
	{
		b[0]=innerRadiusOfTheSample;
		b[1]=outerRadiusOfTheSample;
		b[2]=hightOfTheSample;
// 		b[3]=location_Sample.getX();
// 		b[4]=location_Sample.getY();
// 		b[5]=location_Sample.getZ();
		b[3] = 0;
		b[4] = 0;
		b[5] = 0;
		return b;
	};
  private:

	G4double b[6];
	G4double sizeofLaBr3, sizeofNaI, heightofNaI;

	G4ThreeVector re_location,location_Sample;
	G4double innerRadiusOfTheSample,outerRadiusOfTheSample,hightOfTheSample;

	G4LogicalVolume* Tank_log;
	G4LogicalVolume* Water_log;
	G4LogicalVolume* TankAir_log;
	G4LogicalVolume* detectorPVC_log;
	G4LogicalVolume* detectorAir_log;
	// G4LogicalVolume* detectorAl_log;
	// G4LogicalVolume* crystalouter_log;
	G4LogicalVolume* LaBr3crystal_log;
	G4LogicalVolume* NaIcrystal_log;
	//G4LogicalVolume* detectorFe_log;

	G4VPhysicalVolume* Tank_phys;
	G4VPhysicalVolume* Water_phys;
	G4VPhysicalVolume* TankAir_phys;
	G4VPhysicalVolume* detectorPVC_phys;
	G4VPhysicalVolume* detectorAir_phys;
	// G4VPhysicalVolume* detectorAl_phys;
	// G4VPhysicalVolume* crystalouter_phys;
	G4VPhysicalVolume* LaBr3crystal_phys;
	G4VPhysicalVolume* NaIcrystal_phys;
	//G4VPhysicalVolume* detectorFe_phys;
};

#endif
