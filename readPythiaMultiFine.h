#include "TH1D.h"
#include "TH3D.h"
#include "THnSparse.h"
#define PI 3.14159

// hists for all particles
// TH1I *pdgDist = new TH1I("pdgDist", "PDG Distribution of All Particles", 2751, -5500, 5500);

// Single particle distributions (axes are pt, phi, eta)
TH3D* allParticleDist = new TH3D("allParticleDist", "All Particle Distribution", 100, 0, 10, 100, 0, 2*PI, 100, -10, 10); 
TH3D* chargedHadronDist = new TH3D("chargedHadronDist", "Charged Hadron Distribution", 100, 0, 10, 100, 0, 2*PI, 100, -10, 10); 
TH3D* lambdaDist = new TH3D("lambdaDist", "Lambda Distribution", 100, 0, 10, 100, 0, 2*PI, 100, -10, 10);
TH3D* triggeredLambdaDist = new TH3D("triggeredLambdaDist", "Lambda Distribution in events with a trigger", 100, 0, 10, 100, 0, 2*PI, 100, -10, 10);


// Correlation distributions (axes are trigger eta, associated eta, trigger pt, associated pt, delta phi, delta eta)
// doubled number of delta phi bins from readPythiaMulti
int correlation_num_bins[6] = {20, 20, 20, 20, 32, 40};
double correlation_mins[6] = {-2, -2 , 0, 0, -PI/2, -4};
double correlation_maxes[6] = {2, 2 , 10, 10, 3*PI/2, 4};

THnSparseD *hLambdaDist = new THnSparseD("hLambdaDist", "h-#Lambda Correlation Distribution", 6, correlation_num_bins, correlation_mins, correlation_maxes);
THnSparseD *hhDist = new THnSparseD("hhDist", "h-h Correlation Distribution", 6, correlation_num_bins, correlation_mins, correlation_maxes);