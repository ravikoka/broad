#include "test_readPythia.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "AliStack.h"

// Functions:
void writeOut() 
{
    TFile *outFile = new TFile("test_analysis.root", "RECREATE");


    pdgDist->Write("pdgDist");

    outFile->Close();
} 


// corrects for tangent range in angular correlations
double triggerphi,assophi,deltaphi;
double caldeltaphi(double triggerphi, double assophi) 
{
	double deltaphi = triggerphi - assophi;
	if (deltaphi<-PI/2) {deltaphi += 2*PI;}
	if (deltaphi>=3*PI/2) {deltaphi -= 2*PI;}
	return deltaphi;
}

// calculates delta eta
double trigeta,assoeta,deltata;
double caldeltaeta(double trigeta, double assoeta) {
	double deltaeta = trigeta - assoeta;
	return deltaeta;
}


void test_readPythia()
{
    TStopwatch timer;
    timer.Start();

    // load ROOT file generated with minbiasgen_fixed.c
    string file = "/home/ravikkoka/alice/broad/data/outloop_io_test_1mil.root";
    AliRunLoader* inputRun = AliRunLoader::Open(file.c_str()); //.c_str()
    int numEvents = inputRun->GetNumberOfEvents(); 
    inputRun->LoadKinematics();
    inputRun->LoadHeader();


    // checkpoint: print the number of events
    std::cout << numEvents << std::endl;

    // loop through all events in our run
    for (int event=0; event < numEvents; event++)
        {
        
        // checkpoint: print statement for every 10,000 events
        if ((event % 10000) == 0){std::cout << event << endl;}

        // for each event, create list of particles, do correlations, add these to histogram
        //std::vector<TParticle*> hadron_list;
        //std::vector<TParticle*> lambda_list;

        //std::vector<TParticle*> hadron_trigger_list;
        //std::vector<TParticle*> lambda_trigger_list;

        //std::vector<TParticle*> hadron_assoc_list;
        //std::vector<TParticle*> lambda_assoc_list;

        // for each event, load the stack (list of all particles in an event)
        inputRun->GetEvent(event);
        AliStack *theStack = inputRun->Stack();

        // loop through the stack.
        // GetNTracks() returns the length of our "list" (stack)
        for (int i = 0; i < theStack->GetNtrack(); i++)
            {

            TParticle *particle = theStack->Particle(i);
            
            int pdg = particle->GetPdgCode();
            //double pT = particle->Pt();
            //double phi = particle->Phi();
            //double eta = particle->Eta();

            // fill histograms for all particles
            pdgDist->Fill(pdg);
            //ptDist->Fill(pT);
            //phiDist->Fill(phi);
            //etaDist->Fill(eta);
                
            }                  
            
        
        }

    
    writeOut();

    timer.Stop();
    timer.Print();

}


