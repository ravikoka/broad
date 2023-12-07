#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TParticle.h"
#include "TStopwatch.h"
#include "AliRunLoader.h"
#include "AliStack.h"

#include "readPythiaMulti.h"

// Functions:
void writeOut(TString outFileName)
{
    TFile *outFile = new TFile(outFileName, "RECREATE");
    outFile->cd();

    // write out dists
    allParticleDist->Write();
    chargedHadronDist->Write();
    lambdaDist->Write();
    triggeredLambdaDist->Write();

    hLambdaDist->Write();
    hhDist->Write();

    outFile->Close();
} 

// corrects for tangent range in angular correlations
double calcDeltaPhi(double triggerphi, double assophi) 
{
	double deltaphi = triggerphi - assophi;
	if (deltaphi<-PI/2) {deltaphi += 2*PI;}
	if (deltaphi>=3*PI/2) {deltaphi -= 2*PI;}
	return deltaphi;
}

// calculates delta eta
double calcDeltaEta(double trigeta, double assoeta) 
{
	double deltaeta = trigeta - assoeta;
	return deltaeta;
}

void fillSingleParticleDist(std::vector<TParticle*> particle_list, TH3D* hist)
{
    for (int i = 0; i < particle_list.size(); i++)
    {
        TParticle* particle = particle_list[i];
        // Not really sure what Phi() returns (either -pi to pi or 0 to 2pi), forcing to 0-2pi
        double shifted_phi;
        if(particle->Phi() < 0) shifted_phi = particle->Phi() + 2*PI;
        else shifted_phi = particle->Phi();

        hist->Fill(particle->Pt(), shifted_phi, particle->Eta());
    }
}

void fillCorrelationDist(std::vector<TParticle*> trigger_list, std::vector<TParticle*> associated_list, THnSparseD* hist)
{
    double correlation_array[6];

    for (int i = 0; i < trigger_list.size(); i++)
    {
        TParticle* trigger = trigger_list[i];
        for (int j = 0; j < associated_list.size(); j++)
        {
            // if(trigger[i]->GetUniqueID() == associated[j]-GetUniqueID()) continue; 
            if(i == j) continue; 

            TParticle* associated = associated_list[j];
            double delta_phi = calcDeltaPhi(trigger->Phi(), associated->Phi());
            double delta_eta = calcDeltaEta(trigger->Eta(), associated->Eta());

            correlation_array[0] = trigger->Eta();
            correlation_array[1] = associated->Eta();
            correlation_array[2] = trigger->Pt();
            correlation_array[3] = associated->Pt();
            correlation_array[4] = delta_phi;
            correlation_array[5] = delta_eta;

            hist->Fill(correlation_array);
        }
    }
}


void readPythiaMulti()
{
    TStopwatch timer;
    timer.Start();

    // load ROOT file generated with minbiasgen_fixed.c
    AliRunLoader* inputRun = AliRunLoader::Open("outloop_io_test_1mil.root");
    int numEvents = inputRun->GetNumberOfEvents();
    inputRun->LoadKinematics();
    inputRun->LoadHeader();

    hLambdaDist->Sumw2();
    hhDist->Sumw2();


    // checkpoint: print the number of events
    std::cout << "The number of events: " << numEvents << std::endl;

    // loop through all events in our run
    for (int event=0; event < numEvents; event++)
        {
        
        // checkpoint: print statement for every 10,000 events
        if ((event % 10000) == 0){std::cout << "Processing event: " << event << std::endl;}

        // for each event, create list of particles, fill list of particles, run functions on said lists
        std::vector<TParticle*> all_particle_list;
        std::vector<TParticle*> hadron_list;
        std::vector<TParticle*> lambda_list;

        // bool to keep track of whether or not there was a trigger in the event
        bool found_trigger_hadron = false;

        // for each event, load the stack (list of all particles in an event)
        inputRun->GetEvent(event);
        AliStack *theStack = inputRun->Stack();

        // loop through the stack.
        // GetNTracks() returns the length of our "list" (stack)
        for (int i = 0; i < theStack->GetNtrack(); i++)
        {
            TParticle *particle = theStack->Particle(i);

            all_particle_list.push_back(particle);

            int absPdg = TMath::Abs( particle->GetPdgCode() );

            // select for protons, charged kaons, charged pions, electrons, and their corresponding antiparticles
            // append to list
            if ( (absPdg==2212) || (absPdg==321) || (absPdg==211) || (absPdg==11) ) 
            {
                hadron_list.push_back(particle);
                if(particle->Pt() > 4.0 && particle->Pt() < 8.0) found_trigger_hadron = true;
            }

            // select for lambdas, append to list
            else if ( absPdg==3122 )
            {
                lambda_list.push_back(particle);
            }
        }

        // fill our single particle distributions
        fillSingleParticleDist(all_particle_list, allParticleDist);
        fillSingleParticleDist(hadron_list, chargedHadronDist);
        fillSingleParticleDist(lambda_list, lambdaDist);
        if(found_trigger_hadron) fillSingleParticleDist(lambda_list, triggeredLambdaDist);

        // fill our correlation distributions
        fillCorrelationDist(hadron_list, hadron_list, hhDist);
        fillCorrelationDist(hadron_list, lambda_list, hLambdaDist);
    }

    writeOut("analysisMulti.root");

    timer.Stop();
    timer.Print();

}