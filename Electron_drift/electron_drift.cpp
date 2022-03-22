#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "TRandom3.h"
#include "TBufferJSON.h"
#include "TNtupleD.h"
#include "TFile.h"

#include <vector>
#include <map>
#include <fstream>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ElectronsCoordinates;

std::vector<ElectronsCoordinates> GetInitialPositions(size_t num = 10, double avalancheSize = 20 /*will be removed once coordinates are given*/) {
  TRandom3 rand;
  std::vector<ElectronsCoordinates> electronsPositions(num, ElectronsCoordinates(0, 0, 0, 0));
  // So far random
  size_t counter = 0;
  // For a signal in the top bin of central pitch
  for (auto& electron : electronsPositions)
  {
      // electron = ElectronsCoordinates((1./4 + 1./4 * rand.Uniform())*avalancheSize, (1./4 + 1./4 * rand.Uniform())*avalancheSize, 0.003, 0); // 3nm in
      electron = ElectronsCoordinates((1./4 + 1./8)*avalancheSize, (1./4 + 1./8)*avalancheSize, 0.003, 0); // 3nm in
  }
  // for the central pitch:
  // for (auto& electron : electronsPositions)
  // {
  //     electron = ElectronsCoordinates(rand.Uniform()*pitchSize - pitchSize/2., rand.Uniform()*pitchSize - pitchSize/2., 0.003, 0); // 3nm in
  // }

  return electronsPositions;
}

void electron_drift() {

  // Importning the geometry
  string line;
  ifstream input_file("../Data/Geometry/Test_data.json");
  getline(input_file, line);

  std::map<string, double> *sipmGeomMap = nullptr;
  TBufferJSON::FromJSON(sipmGeomMap, line.data());

  // for (auto& s:*sipmGeomMap) {
  //   // s.second = s.second * 1e3; // From mm to mkm
  //   std::cout << s.first << std::endl;
  // }

  // Physcis params
  double DiffConst=3.9;//um^2/ns D = kT/q mu - Einstein relation
  double Dr=0.05;//um - Step size
  double Dt=Dr*Dr/2./DiffConst; //
  std::cout << "Dt " << Dt << std::endl;
  double tau=100;// e- lifetime in ns
  TRandom3 rand;
  TNtupleD* ntp = new TNtupleD("ntp","ntp","initialX:initialY:initialZ:initialT:finalX:finalY:finalZ:finalT");


  // Defining a start position
  std::vector<ElectronsCoordinates> electronsInitialPositions = GetInitialPositions(100000, (*sipmGeomMap)["av_width"]);
  std::vector<ElectronsCoordinates> electronsFinalPositions;
  electronsFinalPositions.reserve( electronsInitialPositions.size() ); // Actually we can fill ntp right away but just in case store them so far
  size_t counter = 0;
  for (auto electronPos : electronsInitialPositions){
    double x = electronPos.X();
    double y = electronPos.Y();
    double z = electronPos.Z();
    double t = electronPos.t();
    while(rand.Uniform()>Dt/tau && z < (*sipmGeomMap)["electrons_drift_depth"]){
      double Dx,Dy,Dz;
      rand.Sphere(Dx,Dy,Dz,Dr);
      x+=Dx; 
      y+=Dy;
      z+=Dz;
      // Maybe change to the condition that Dz < 0 -> Dz *= -1?
      if(z<0.) z*=-1.; //reflect at surface - push back by surface Efield
      t+=Dt;
    }
    electronsFinalPositions.push_back(ElectronsCoordinates(x, y, z, t));
    counter++;
    if (counter % 5000 == 0) std::cout << "Event number " << counter << std::endl;
  }

  for (size_t i = 0; i < electronsInitialPositions.size(); ++i) {
    double toTuple[8] = {electronsInitialPositions[i].X(), electronsInitialPositions[i].Y(), electronsInitialPositions[i].Z(), electronsInitialPositions[i].T(),
                         electronsFinalPositions[i].X(), electronsFinalPositions[i].Y(), electronsFinalPositions[i].Z(), electronsFinalPositions[i].T()};
    ntp->Fill(toTuple);
  }
  TFile* fOutFile = new TFile("../Data/ElectronDriftOutput/100k_point.root","RECREATE");
  // TFile* fOutFile = new TFile("test_fabrice.root","RECREATE");
  // std::cout << "Output filename " << fOutFileName << std::endl;
  ntp->Write();
  TObject DtWrite;
  DtWrite.SetUniqueID(Dt);
  DtWrite.Write("Dt");
  fOutFile->Close();

  // for (auto electron:electronsFinalPositions) {
  //   std::cout << electron.Z() << std::endl;
  // }
}

// void electron_drift(){

// // Physcis params
// double DiffConst=3.9;//um^2/ns D = kT/q mu - Einstein relation
// double Dr=0.05;//um - Step size
// double Dt=Dr*Dr/2./DiffConst; //
// cout << "Dt " << Dt << endl;
// double tau=100;// e- lifetime in ns
// int nMC=100;

// char ntpName[20];
// sprintf(ntpName,"ntp%i",0);
// TNtuple* ntp = (TNtuple*) gROOT->FindObjectAny(ntpName);
// if(ntp) ntp->Delete();
// ntp = new TNtuple(ntpName,ntpName,"pitch:th:ne:tf:nspad");
// TRandom3 rand;
// // start at 0,0,0. End either if e- dies or reach high field region

// //TH2D* HNSPAD = (TH2D*) gROOT->FindObjectAny("HNSPAD");
// //if(HNSPAD) HTTS->Delete();
// //HNSPAD = new TH2D("HNSPAD","HNSPAD",6,5.,55.,10,10.,210.);

// // for(double pitch=10.; pitch<=40.; pitch+=5.){
// //   cout << "Pitch " << pitch << endl;
// //   for(double th=1.; th<=6.; th+=0.5){
// //     cout << "Thickness " << th << endl;
// //     for(int nEl=1; nEl<=100; nEl+=99){
// //       cout << "# electrons " << nEl << endl;

//       double pitch = 40;
//       double th = 0.1;
//       double nEl = 10;
//       // for(int iMC=0; iMC<nMC; iMC++){
//         double tFirst=1e6;

//         // Making spad zeros
//         int iSPAD[20][20];
//         for(int ii=0; ii<20; ii++){
//           for(int jj=0; jj<20; jj++){
//             iSPAD[ii][jj]=0;
//           }
//         }

//         // Defining a start position
//         double xStart=rand.Uniform()*pitch;
//         double yStart=rand.Uniform()*pitch;
//         for(int iEl=0; iEl<nEl; iEl++){
//         double x = xStart;
//         double y = yStart;
//         double z=0.003; //3nm in
//         double t=0.;
//         while(rand.Uniform()>Dt/tau && z<th){
//           double Dx,Dy,Dz;
//           rand.Sphere(Dx,Dy,Dz,Dr);
//           x+=Dx;
//           y+=Dy;
//           z+=Dz;
//           // Maybe change to the condition that Dz < 0 -> Dz *= -1?
//           if(z<0.) z*=-1.; //reflect at surface - push back by surface Efield
//           t+=Dt;
//         }
//         cout << "Z last " << z << std::endl;  
//         int ix = (int)x/pitch+10;
//         cout << "ix " << ix << std::endl;  

//         if(ix<0 || ix>=20) cout << "x out of bound" << ix << endl;
//         else{
//          int iy = (int)y/pitch+10;
//           cout << "iy " << iy << std::endl;  
//          if(iy<0 || iy>=20) cout << "y out of bound" << iy << endl;
//          else iSPAD[ix][iy]++;
//         }
//         if(t<tFirst) tFirst=t;
//         }
//         int nSPAD=0;
//         for(int ii=0; ii<20; ii++){
//           for(int jj=0; jj<20; jj++){
//             if(iSPAD[ii][jj]) nSPAD++;
//           }
//         }
//         ntp->Fill(pitch,th,nEl,tFirst,nSPAD);

// TFile* fOutFile = new TFile("test.root","RECREATE");
// // std::cout << "Output filename " << fOutFileName << std::endl;
// ntp->Write();
// fOutFile->Close();

// }
