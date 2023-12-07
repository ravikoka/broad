#include "TH1D.h"
#include <map>
#define PI 3.14159

// hists for trigger & assoc eta dists
TH1D *etaDistHad = new TH1D("etaDistHad", "Eta Distribution of All Hadrons", 1000, -10, 10);
TH1D *etaDistLam = new TH1D("etaDistLam", "Eta Distribution of All Lambdas", 1000, -10, 10);

TH1D *etaDistTrigHad = new TH1D("etaDistTrigHad", "Eta Distribution of Trigger Hadrons", 1000, -10, 10);
TH1D *etaDistTrigLam = new TH1D("etaDistTrigLam", "Eta Distribution of Trigger Lambdas", 1000, -10, 10);
TH1D *etaDistAssocHad = new TH1D("etaDistAssocHad", "Eta Distribution of Associate Hadrons", 1000, -10, 10);
TH1D *etaDistAssocLam = new TH1D("etaDistAssocLam", "Eta Distribution of Associate Lambdas", 1000, -10, 10);

// particle dictionary, made by Amanda
std::map<int, string> PARTICLEDICT = { {1,"d"}, {2,"u"}, {3,"s"}, {4,"c"}, {5,"b"}, {6,"t"}, {7,"b'"}, {8,"t'"}, {11,"e-"}, {12,"nu_e"}, {13,"mu-"}, {14,"nu_mu"}, {15,"tau-"}, {16,"nu_tau"}, {17,"tau'"}, {18,"nu'_tau"}, {21,"g"}, {22,"photon"}, {23,"Z0"}, {24,"W+"}, {25,"h0"}, {81,"spectflav"}, {82,"rndmflav"}, {83,"phasespa"}, {84,"c-hadron"}, {85,"b-hadron"}, {88,"junction"}, {90,"system"}, {91,"cluster"}, {92,"string"}, {93,"indep"}, {94,"CMshower"}, {95,"SPHEaxis"}, {96,"THURaxis"}, {97,"CLUSjet"}, {98,"CELLjet"}, {99,"table"}, {111,"pi0"}, {113,"rho0"}, {130,"K_L0"}, {211,"pi+-"}, {213,"rho+-"}, {221,"eta"}, {223,"omega"}, {310,"K_S0"}, {311,"K0"}, {313,"K*0"}, {321,"K+-"}, {323,"K*+-"}, {331,"eta'"}, {333,"phi"}, {411,"D+-"}, {413,"D*+-"}, {421,"D0"}, {423,"D*0"}, {431,"D_s+-"}, {433,"D*_s+-"}, {441,"eta_c"}, {443,"J/psi"}};
