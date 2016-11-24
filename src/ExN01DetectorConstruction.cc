// $Id: ExN01DetectorConstruction.cc,v 1.9 2006/06/29 17:47:19 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#include "ExN01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4AssemblyVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
ExN01DetectorConstruction::ExN01DetectorConstruction
	(G4double size1,G4double size2, G4double size3)
	:sizeofLaBr3(size1),sizeofNaI(size2),heightofNaI(size3)
{;}

ExN01DetectorConstruction::~ExN01DetectorConstruction()
{;}

G4VPhysicalVolume* ExN01DetectorConstruction::Construct()
{

   //------------------------------------------------------ materials

   G4String symbol;             //a=mass of a mole;
   G4double a, z, density;      //z=mean number of protons;

   G4int ncomponents, natoms;
   G4double fractionmass;
   //
   // define Elements
   //

   G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
   G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
   G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
   G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
   G4Element* F  = new G4Element("Fluorine",symbol="F" , z= 9., a= 19.00*g/mole);
   G4Element* Al = new G4Element("Aluminum",symbol="Al", z= 13.,a= 26.98*g/mole);
   G4Element* Cl = new G4Element("Chlorine",symbol="Cl", z= 17.,a= 35.45*g/mole);
   G4Element* Cu = new G4Element("Copper"  ,symbol="Cu", z= 29.,a= 63.55*g/mole);
   G4Element* Zn = new G4Element("Zinc"    ,symbol="Zn", z= 30.,a= 65.39*g/mole);
   G4Element* La = new G4Element("Lanthanum",symbol="La",z= 57.,a= 138.91*g/mole);
   G4Element* Br = new G4Element("Bromine",symbol="Br",z= 35.,a= 79.90*g/mole);
   G4Element* Mg = new G4Element("Magnesium",symbol="Mg",z= 12.,a= 24.31*g/mole);
   G4Element* Na = new G4Element("Sodium",symbol="Na",z= 11.,a= 22.99*g/mole);
   G4Element* S = new G4Element("Sulphur",symbol="S",z= 16.,a= 32.06*g/mole);
   G4Element* Ca = new G4Element("Calcium",symbol="Ca",z= 20.,a= 40.08*g/mole);
   G4Element* K = new G4Element("Potassium",symbol="K",z= 19.,a= 39.10*g/mole);
   G4Element* Sr = new G4Element("Strontium",symbol="Sr",z= 38.,a= 87.62*g/mole);
   G4Element* B = new G4Element("Boron",symbol="B",z= 5.,a= 10.81*g/mole);
   G4Element* Si = new G4Element("Silicon",symbol="Si",z= 14.,a= 28.09*g/mole);
   G4Element* Fe = new G4Element("Iron",symbol="Fe",z= 26.,a= 55.85*g/mole);
   G4Element* Mn = new G4Element("Manganese",symbol="Mn",z= 25.,a= 54.94*g/mole);
   G4Element* P = new G4Element("Phosphorus",symbol="P",z= 15.,a= 30.97*g/mole);
   G4Element* Cr = new G4Element("Chrome",symbol="Cr",z= 24.,a= 52.00*g/mole);
   G4Element* Ni = new G4Element("Nickel", symbol = "Ni", z = 28., a = 58.69*g / mole);
   G4Element* Mo = new G4Element("Molybdenum", symbol = "Mo", z = 42., a = 95.94*g / mole);
   G4Element* I = new G4Element("Iodine", symbol = "I", z = 53., a = 126.90447*g / mole);

   G4Material* Vacuum =//����
	new G4Material("Vacuum", z=1., a=1.01*g/mole,density= universe_mean_density,
		kStateGas, 2.73*kelvin, 3.e-18*pascal);

   G4Material* LaBr3 = //LaBr3
	   new G4Material("LaBr3", density = 5.08*g / cm3, ncomponents = 2);
   LaBr3->AddElement(La, natoms = 1);
   LaBr3->AddElement(Br, natoms = 3);
//    new G4Material("LaBr3", density = 3.67*g / cm3, ncomponents = 2);
//    LaBr3->AddElement(Na, natoms = 1);
//    LaBr3->AddElement(I, natoms = 1);


   G4Material* LY12 = //���Ͻ�
	new G4Material("LY12" , density= 2.78*g/cm3, ncomponents=3);
   LY12->AddElement(Al, fractionmass=0.94);
   LY12->AddElement(Cu, fractionmass=0.045);
   LY12->AddElement(Mg, fractionmass=0.015);

   G4Material* Seawater = //��ˮ
	   new G4Material("Seawater",density = 1.025*g/cm3,ncomponents=15);
   Seawater->AddElement(H, fractionmass=0.108323963);
   Seawater->AddElement(O, fractionmass=0.859713989);
   Seawater->AddElement(Cl,fractionmass=0.018419512);
   Seawater->AddElement(Na,fractionmass=0.010507317);
   Seawater->AddElement(Mg,fractionmass=0.001258537);
   Seawater->AddElement(S, fractionmass=0.000862439);
   Seawater->AddElement(Ca,fractionmass=0.000402049);
   Seawater->AddElement(K, fractionmass=0.000389268);
   Seawater->AddElement(Br,fractionmass=0.000065659);
   Seawater->AddElement(C, fractionmass=0.000027317);
   Seawater->AddElement(N, fractionmass=0.000014634);
   Seawater->AddElement(Sr,fractionmass=0.000007707);
   Seawater->AddElement(B, fractionmass=0.000004390);
   Seawater->AddElement(Si,fractionmass=0.000001951);
   Seawater->AddElement(F, fractionmass=0.000001268);

   G4Material* Steel = //������
	   new G4Material("Steel",density = 7.89*g/cm3,ncomponents=9);
   Steel->AddElement(C, fractionmass=0.0008);
   Steel->AddElement(Si, fractionmass=0.01);
   Steel->AddElement(Mn,fractionmass=0.02);
   Steel->AddElement(P,fractionmass=0.00035);
   Steel->AddElement(S,fractionmass=0.0003);
   Steel->AddElement(Ni, fractionmass=0.12);
   Steel->AddElement(Cr,fractionmass=0.17);
   Steel->AddElement(Mo, fractionmass=0.015);
   Steel->AddElement(Fe, fractionmass=0.66355);

   G4Material* Air = //����
   new G4Material("Air" , density= 1.205*mg/cm3, ncomponents=2);
   Air->AddElement(N, fractionmass=0.7);
   Air->AddElement(O, fractionmass=0.3);

   G4Material* FeZn =
   new G4Material("FeZn" , density= 7.87*g/cm3, ncomponents=2);
   FeZn->AddElement(Fe, fractionmass=0.99);
   FeZn->AddElement(Zn, fractionmass=0.01);

   G4Material* PVC = //������ϩ
   new G4Material("PVC" , density= 1.38*g/cm3, ncomponents=3);
   PVC->AddElement(C, natoms=2);
   PVC->AddElement(H, natoms=3);
   PVC->AddElement(Cl, natoms=1);;

   G4Material* PE = //����ϩ
   new G4Material("PE" , density= 0.95*g/cm3, ncomponents=2);
   PE->AddElement(C, natoms=1);
   PE->AddElement(H, natoms=2);

   G4Material* NaI =
	   new G4Material("NaI", density = 3.67*g / cm3, ncomponents = 2);
   NaI->AddElement(Na, natoms = 1);
   NaI->AddElement(I, natoms = 1);

	 // Option to switch on/off checking of volumes overlaps
   //
   G4bool checkOverlaps = true;

   //
   // World
   //
   G4double world_sizeXY = 3000.*mm;
   G4double world_sizeZ  = 3000.*mm;
  //  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

   G4Box* solidWorld =
     new G4Box("World",                       //its name
        0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

   G4LogicalVolume* logicWorld =
     new G4LogicalVolume(solidWorld,          //its solid
                         Vacuum,           //its material
                         "World");            //its name

   G4VPhysicalVolume* physWorld =
     new G4PVPlacement(0,                     //no rotation
                       G4ThreeVector(),       //at (0,0,0)
                       logicWorld,            //its logical volume
                       "World",               //its name
                       0,                     //its mother  volume
                       false,                 //no boolean operation
                       0,                     //copy number
                       checkOverlaps);        //overlaps checking


   //����ˮ

   //eff = 1.;
   G4double thickofTank = 10.*mm;
   //-------ˮ��----a Tank tube
   G4double innerRadiusOfTheTube = 0./2*mm;
   G4double outerRadiusOfTheTube = 2000./2*mm;
   G4double halfHeightOfTheTube = 2300./2*mm;
   G4double startAngleOfTheTube = 0.*deg;
   G4double spanningAngleOfTheTube = 360.*deg;
   G4Tubs* Tank_tube = new G4Tubs("Tank_tube",innerRadiusOfTheTube,
	   outerRadiusOfTheTube,halfHeightOfTheTube,startAngleOfTheTube,spanningAngleOfTheTube);
   //G4Sphere* Tank_tube = new G4Sphere("Tank_tube",0.*mm,1000.*mm,0.*deg,360.*deg,0.*deg,360.*deg);
   Tank_log = new G4LogicalVolume(Tank_tube,PE,"Tank_log",0,0,0);
   Tank_phys = new G4PVPlacement(0,G4ThreeVector(),Tank_log,"Tank",logicWorld,false,0, checkOverlaps);

   //-------��ˮ----a Water tube
   innerRadiusOfTheTube = 0./2*mm;
   outerRadiusOfTheTube = outerRadiusOfTheTube - thickofTank;
   halfHeightOfTheTube = halfHeightOfTheTube - thickofTank;
   startAngleOfTheTube = 0.*deg;
   //G4Sphere* Water_tube = new G4Sphere("Tank_tube",0.*mm,990.*mm,0.*deg,360.*deg,0.*deg,360.*deg);
   G4Tubs* Water_tube = new G4Tubs("Water_tube",innerRadiusOfTheTube,
	   outerRadiusOfTheTube,halfHeightOfTheTube,startAngleOfTheTube,spanningAngleOfTheTube);
   Water_log = new G4LogicalVolume(Water_tube,Seawater,"Water_log",0,0,0);
   Water_phys = new G4PVPlacement(0,G4ThreeVector(), Water_log,"Water",Tank_log,false,0, checkOverlaps);

   G4double thickofAl = 0.*mm,thickofPVC = 4.0*mm,thickofAir=0.5*mm;
	 G4double heightOfPMT  = 300 *mm;
   G4double inch = 25.4*mm;
   G4double eff = sizeofLaBr3 / 3;

   //------����----a detectorPVC
 //  G4double spanningAngleOfTheTube = 180.*deg;
   outerRadiusOfTheTube = (sizeofLaBr3 + sizeofNaI)*inch / 2 + thickofAl + thickofPVC + thickofAir;
   halfHeightOfTheTube =  heightofNaI*inch / 2. + thickofAl + thickofPVC + thickofAir + heightOfPMT/2.;
   G4Tubs* PVCtube1 = new G4Tubs("PVCtube1",0.*mm,outerRadiusOfTheTube,halfHeightOfTheTube,0.*deg,spanningAngleOfTheTube);
  //  G4Tubs* PVCtube2 = new G4Tubs("PVCtube2",0.*mm,outerRadiusOfTheTube,360.*eff/2*mm,0.*deg,spanningAngleOfTheTube);
  //  G4Tubs* PVCtube3 = new G4Tubs("PVCtube3",0.*mm,160.*eff/2*mm,75.*eff/2*mm,0.*deg,spanningAngleOfTheTube);
  //  G4UnionSolid* detectorPVC1 = new G4UnionSolid("detectorPVC1",PVCtube1, PVCtube2,0,G4ThreeVector(0.*mm,0.*mm, halfHeightOfTheTube+360.*eff/2*mm));
  //  G4UnionSolid* detectorPVC = new G4UnionSolid("detectorPVC", detectorPVC1, PVCtube3, 0, G4ThreeVector(0.*mm, 0.*mm, halfHeightOfTheTube + 360.*eff*mm + 75. *eff / 2 * mm));
	 G4double shiftToLaBrCenter = (-sizeofLaBr3 + heightofNaI) * inch / 2 + heightOfPMT/2.;
   G4ThreeVector location_PVCOuter =  G4ThreeVector(0,0,-shiftToLaBrCenter);
  //  re_location = location_PVCOuter - G4ThreeVector(0,0,0);

   detectorPVC_log = new G4LogicalVolume(PVCtube1,PVC,"detectorPVC_log",0,0,0);
   detectorPVC_phys = new G4PVPlacement(0, location_PVCOuter,detectorPVC_log,"detectorPVC",Water_log,false,0, checkOverlaps);


   //------�����ڿ���----detectorAir
   outerRadiusOfTheTube = (sizeofLaBr3 + sizeofNaI)*inch / 2 + thickofAl + thickofAir;
   halfHeightOfTheTube = heightofNaI*inch / 2. + thickofAl + thickofAir + heightOfPMT/2.;
   G4Tubs* Airtube1 = new G4Tubs("Airtube1",0.*mm, outerRadiusOfTheTube, halfHeightOfTheTube,0.*deg,spanningAngleOfTheTube);
  //  G4Tubs* Airtube2 = new G4Tubs("Airtube2",0.*mm,96.*eff/2*mm,(360.*mm+60.*mm+thickofPVC)*eff/2,0.*deg,spanningAngleOfTheTube);
  //  G4UnionSolid* detectorAir = new G4UnionSolid("detectorAir", Airtube1, Airtube2,0,G4ThreeVector(0.*mm,0.*mm,halfHeightOfTheTube+ (360.*mm + 60.*mm + thickofPVC)*eff / 2));

   G4ThreeVector location_PVCInner =  G4ThreeVector(0,0,0);
  //  re_location = location_PVCInner + location_PVCOuter;
   detectorAir_log = new G4LogicalVolume(Airtube1,Air,"detectorAir_log",0,0,0);
   detectorAir_phys = new G4PVPlacement(0,location_PVCInner, detectorAir_log,"detectorAir",detectorPVC_log,false,0, checkOverlaps);

  //  G4double pSPhi = 0.*deg;
  //  G4double pDPhi = 360.*deg;

  //  //����(������)
  //  outerRadiusOfTheTube = (sizeofLaBr3 + sizeofNaI)*inch / 2 + thickofAl;
  //  halfHeightOfTheTube = sizeofLaBr3*inch / 2 + thickofAl;
  //  G4Tubs* Altube1 = new G4Tubs("Altube1", 0.*mm, outerRadiusOfTheTube, halfHeightOfTheTube, 0.*deg, spanningAngleOfTheTube);
  //  G4Tubs* Altube2 = new G4Tubs("Altube2", 0.*mm, 81.*eff / 2 * mm, 180.*mm*eff / 2, 0.*deg, spanningAngleOfTheTube);
  //  G4UnionSolid* detectorAl = new G4UnionSolid("detectorAir", Altube1, Altube2, 0, G4ThreeVector(0.*mm, 0.*mm, halfHeightOfTheTube + 180.*eff / 2 * mm));
	 //
  //  G4ThreeVector location_Al = G4ThreeVector(0, 0, 0);
  //  re_location = location_Al - location_PVCInner;
  //  detectorAl_log = new G4LogicalVolume(detectorAl, LY12, "detectorAl_log", 0, 0, 0);
  //  detectorAl_phys = new G4PVPlacement(0, re_location, detectorAl_log, "detectorAl", detectorAir_log, false, 0);


  //  //�����ⲿ����
  //  outerRadiusOfTheTube = (sizeofLaBr3 + sizeofNaI)*inch / 2;
  //  halfHeightOfTheTube = sizeofLaBr3*inch / 2;
  //  G4Tubs* Airtube3 = new G4Tubs("Altube3", 0.*mm, outerRadiusOfTheTube, halfHeightOfTheTube, 0.*deg, spanningAngleOfTheTube);
  //  G4Tubs* Airtube4 = new G4Tubs("Altube4", 0.*mm, (81.*mm/2*eff)-thickofAl, 180.*mm*eff / 2, 0.*deg, spanningAngleOfTheTube);
  //  G4UnionSolid* crystalouter = new G4UnionSolid("detectorAir", Airtube3, Airtube4, 0, G4ThreeVector(0.*mm, 0.*mm, halfHeightOfTheTube + 180.*eff / 2 * mm));
	 //
  //  G4ThreeVector location_Air = G4ThreeVector(0, 0, 0);
  //  re_location = location_Air - location_Al;
  //  crystalouter_log = new G4LogicalVolume(crystalouter, Air, "crystalouter_log", 0, 0, 0);
  //  crystalouter_phys = new G4PVPlacement(0, re_location, crystalouter_log, "crystalouter", detectorAl_log, false, 0);

   //NaI����
   outerRadiusOfTheTube = (sizeofLaBr3 + sizeofNaI)*inch / 2;
   halfHeightOfTheTube = heightofNaI*inch / 2;
   G4Tubs* NaItube = new G4Tubs("NaItube", sizeofLaBr3*inch / 2, outerRadiusOfTheTube, halfHeightOfTheTube, 0.*deg, spanningAngleOfTheTube);

   G4ThreeVector location_NaI = G4ThreeVector(0, 0, (sizeofLaBr3 - heightofNaI)*inch / 2 + shiftToLaBrCenter);
  //  re_location = location_NaI - location_Air;
   NaIcrystal_log = new G4LogicalVolume(NaItube, NaI, "NaIcrystal_log", 0, 0, 0);
   NaIcrystal_phys = new G4PVPlacement(0, location_NaI, NaIcrystal_log, "NaIcrystal", detectorAir_log, false, 0, checkOverlaps);

   //LabR3����
   outerRadiusOfTheTube = sizeofLaBr3*inch / 2;
   halfHeightOfTheTube = sizeofLaBr3*inch / 2;
   G4Tubs* LaBr3tube = new G4Tubs("LaBr3tube", 0.*mm, outerRadiusOfTheTube, halfHeightOfTheTube, 0.*deg, spanningAngleOfTheTube);

   G4ThreeVector location_LaBr3 = G4ThreeVector(0, 0, shiftToLaBrCenter);
  //  re_location = location_LaBr3 - location_NaI;
   LaBr3crystal_log = new G4LogicalVolume(LaBr3tube, LaBr3, "LaBr3crystal_log", 0, 0, 0);
   LaBr3crystal_phys = new G4PVPlacement(0, location_LaBr3, LaBr3crystal_log, "LaBr3crystal", detectorAir_log, false, 0, checkOverlaps);


	G4Colour white(1.0,1.0,1.0);	//��
	G4Colour gray(0.5,0.5,0.5);		//��
	G4Colour black(0.0,0.0,0.0);	//��
	G4Colour red(1.0,0.0,0.0);		//��
	G4Colour green(0.0,1.0,0.0);	//��
	G4Colour blue(0.0,0.0,1.0);		//��
	G4Colour cyan(0.0,1.0,1.0);		//��
	G4Colour magenta(1.0,0.0,1.0);	//��
	G4Colour yellow(1.0,1.0,0.0);	//��
	G4Colour Cu_colour(1.0,0.588,0.314);	//ͭɫ
	G4Colour Vespel_colour(0.353,0.843,0.);	//Vespelɫ
	G4Colour PTFCE_colour(0.5,0.,0.882);	//PTFCEɫ
	G4Colour Al_colour(0.863,0.875,0.89);//��ɫ
	G4Colour Tank_colour(0.918,0.918,0.918);//ˮ��
	G4Colour Sea_colour(0.573,0.859,0.961);//��ˮ
	G4Colour PVC_colour(0.204,0.435,0.714);//����
	G4Colour Steel_colour(0.502,0.502,0.494);//������
	G4Colour color(0.5,1.0,1.0);//�ʺ���ˮ
	G4Colour color1(0.25,0.5,0.75);//�ʺ���ˮ
	G4Colour color2(0.75,0.75,0.75);//�ʺ���ˮ

	G4VisAttributes* TankVisAtt = new G4VisAttributes(gray);
	Tank_log ->SetVisAttributes(TankVisAtt);

	G4VisAttributes* WaterVisAtt = new G4VisAttributes(color);
	Water_log ->SetVisAttributes(WaterVisAtt);

	G4VisAttributes* PVCVisAtt = new G4VisAttributes(color2);
	detectorPVC_log ->SetVisAttributes(PVCVisAtt);

	G4VisAttributes* AirVisAtt = new G4VisAttributes(white);
	detectorAir_log ->SetVisAttributes(AirVisAtt);
	// crystalouter_log->SetVisAttributes(AirVisAtt);

	G4VisAttributes* AlVisAtt = new G4VisAttributes(yellow);
	// detectorAl_log ->SetVisAttributes(AlVisAtt);

	G4VisAttributes* NaIVisAtt = new G4VisAttributes(cyan);
	NaIcrystal_log ->SetVisAttributes(NaIVisAtt);

	G4VisAttributes* Labr3VisAtt = new G4VisAttributes(blue);
	LaBr3crystal_log ->SetVisAttributes(Labr3VisAtt);

	return Tank_phys;
}
