#include "TH1D.h"
#include <map>
#define PI 3.14159

// hists for all particles
TH1I *pdgDist = new TH1I("pdgDist", "PDG Distribution of All Particles", 2751, -5500, 5500);
TH1D *ptDist = new TH1D("ptDist", "Pt Distribution of All Particles", 30, 0, 8);
TH1D *phiDist = new TH1D("phiDist", "Phi Distribution of All Particles", 1000, 0, 2*PI);
TH1D *etaDist = new TH1D("etaDist", "Eta Distribution of All Particles", 1000, -10, 10);


// hists with no cuts
TH1I *pdgDistHadLam = new TH1I("pdgDistHadLam", "H-L PDG Distribution", 2751, -5500, 5500);
TH1D *hadPtDist = new TH1D("hadPtDist", "Hadron PT Distribution", 100, -10, 30);
TH1D *lamPtDist = new TH1D("lamPtDist", "Hadron PT Distribution", 100, -10, 30);
TH1D *phiDistHadLam = new TH1D("phiDistHadLam", "H-L Phi Distribution", 1000, 0, 2*PI);
TH1D *etaDistHadLam = new TH1D("etaDistHadLam", "H-L Eta Distribution", 1000, -10, 10);

TH1D *dPhiHadHad = new TH1D("dPhiHadHad", "H-H Delta Phi Distribution", 128, -PI/2, 3*PI/2);
TH1D *dPhiLamLam = new TH1D("dPhiLamLam", "L-L Delta Phi Distribution", 128, -PI/2, 3*PI/2);
TH1D *dPhiHadLam = new TH1D("dPhiHadLam", "H-L Delta Phi Distribution", 128, -PI/2, 3*PI/2);
TH1D *dPhiLamDivHad = new TH1D("dPhiLamDivHad", "L/H Delta Phi distribution", 128, -PI/2, 3*PI/2);

TH2D *dPhiDEtaHadHad = new TH2D("dPhiDEtaHadHad", "H-H Delta Phi Delta Eta", 128, -PI/2, 3*PI/2, 1000, -10, 10);
TH2D *dPhiDEtaLamLam = new TH2D("dPhiDEtaLamLam", "L-L Delta Phi Delta Eta", 128, -PI/2, 3*PI/2, 1000, -10, 10);
TH2D *dPhiDEtaHadLam = new TH2D("dPhiDEtaHadLam", "H-L Delta Phi Delta Eta", 128, -PI/2, 3*PI/2, 1000, -10, 10);
TH2D *dPhiDEtaLamDivHad = new TH2D("dPhiDEtaLamDivHad", "L/H Delta Phi Delta Eta", 128, -PI/2, 3*PI/2, 1000, -10, 10);

TH1D *lamWithAntiLamCounter = new TH1D("lamWithAntiLamCounter", "Number of L/Lbar Trigger w/ Lbar/L Assoc", 12, -0.5, 2);


// hists with eta cuts
TH1I *pdgDistCutHadLam = new TH1I("pdgDistCuHadLamt", "H-L PDG Distribution for Cut", 2751, -5500, 5500);
TH1D *hadPtDistCut = new TH1D("hadPtDistCut", "H pT Distribution for Cut", 100, -10, 30);
TH1D *lamPtDistCut = new TH1D("lamPtDistCut", "L pT Distribution for Cut", 100, -10, 30);
TH1D *phiDistCutHadLam = new TH1D("phiDistCutHadLam", "H-L Phi Distribution for Cut", 1000, 0, 2*PI);
TH1D *etaDistCutHadLam = new TH1D("etaDisCutHadLam", "H-L Eta Distribution for Cut", 1000, -10, 10);

TH1D *dPhiCutHadHad = new TH1D("dPhiCutHadHad", "H-H Delta Phi for Cut", 128, -PI/2, 3*PI/2);
TH1D *dPhiCutLamLam = new TH1D("dPhiCutLamLam", "L-L Delta Phi for Cut", 128, -PI/2, 3*PI/2);
TH1D *dPhiCutHadLam = new TH1D("dPhiCutHadLam", "H-L Delta Phi for Cut", 128, -PI/2, 3*PI/2);
TH1D *dPhiCutLamDivHad = new TH1D();

TH2D *dPhiDEtaCutHadHad = new TH2D("dPhiDEtaCutHadHad", "H-H Delta Phi Delta Eta for Cut", 128, -PI/2, 3*PI/2, 1000, -10, 10);
TH2D *dPhiDEtaCutLamLam = new TH2D("dPhiDEtaCutLamLam", "L-L Delta Phi Delta Eta for Cut", 128, -PI/2, 3*PI/2, 1000, -10, 10);
TH2D *dPhiDEtaCutHadLam = new TH2D("dPhiDEtaCutHadLam", "H-L Delta Phi Delta Eta for Cut", 128, -PI/2, 3*PI/2, 1000, -10, 10);
TH2D *dPhiDEtaCutLamDivHad = new TH2D();

TH1D *lamWithAntiLamCutCounter = new TH1D("lamWithAntiLamCutCounter", "Number of L/Lbar Trigger w/ Lbar/L Assoc for Cut", 12, -0.5, 2);

// hists for trigger & assoc eta dists
TH1D *etaDistHad = new TH1D("etaDistHad", "Eta Distribution of All Hadrons", 1000, -10, 10);
TH1D *etaDistLam = new TH1D("etaDistLam", "Eta Distribution of All Lambdas", 1000, -10, 10);

TH1D *etaDistTrigHad = new TH1D("etaDistTrigHad", "Eta Distribution of Trigger Hadrons", 1000, -10, 10);
TH1D *etaDistTrigLam = new TH1D("etaDistTrigLam", "Eta Distribution of Trigger Lambdas", 1000, -10, 10);
TH1D *etaDistAssocHad = new TH1D("etaDistAssocHad", "Eta Distribution of Associate Hadrons", 1000, -10, 10);
TH1D *etaDistAssocLam = new TH1D("etaDistAssocLam", "Eta Distribution of Associate Lambdas", 1000, -10, 10);

// particle dictionary, made by Amanda
std::map<int, string> PARTICLEDICT = { {1,"d"}, {2,"u"}, {3,"s"}, {4,"c"}, {5,"b"}, {6,"t"}, {7,"b'"}, {8,"t'"}, {11,"e-"}, {12,"nu_e"}, {13,"mu-"}, {14,"nu_mu"}, {15,"tau-"}, {16,"nu_tau"}, {17,"tau'"}, {18,"nu'_tau"}, {21,"g"}, {22,"photon"}, {23,"Z0"}, {24,"W+"}, {25,"h0"}, {81,"spectflav"}, {82,"rndmflav"}, {83,"phasespa"}, {84,"c-hadron"}, {85,"b-hadron"}, {88,"junction"}, {90,"system"}, {91,"cluster"}, {92,"string"}, {93,"indep"}, {94,"CMshower"}, {95,"SPHEaxis"}, {96,"THURaxis"}, {97,"CLUSjet"}, {98,"CELLjet"}, {99,"table"}, {111,"pi0"}, {113,"rho0"}, {130,"K_L0"}, {211,"pi+-"}, {213,"rho+-"}, {221,"eta"}, {223,"omega"}, {310,"K_S0"}, {311,"K0"}, {313,"K*0"}, {321,"K+-"}, {323,"K*+-"}, {331,"eta'"}, {333,"phi"}, {411,"D+-"}, {413,"D*+-"}, {421,"D0"}, {423,"D*0"}, {431,"D_s+-"}, {433,"D*_s+-"}, {441,"eta_c"}, {443,"J/psi"}};
