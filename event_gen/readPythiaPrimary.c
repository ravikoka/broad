#include "readPythia.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "AliStack.h"
#include "AliMCParticle.h"

// Functions:
void writeOut() 
{
    TFile *outFile = new TFile("analysisprimary.root", "RECREATE");


    pdgDist->Write("pdgDist");
    ptDist->Write("ptDistTotal");
    phiDist->Write("phiDistTotal");
    etaDist->Write("etaDistTotal");

    pdgDistHadLam->Write("pdgDistHadLam");
    hadPtDist->Write("hadronPtDist");
    lamPtDist->Write("lambdaPtDist");
    phiDistHadLam->Write("phiDistHadLam");
    etaDistHadLam->Write("etaDistHadLam");
    dPhiHadHad->Write("dPhiHadHad");
    dPhiLamLam->Write("dPhiLamLam");
    dPhiHadLam->Write("dPhiHadLam");
    dPhiLamDivHad->Write("dPhiLamDivHad");
    dPhiDEtaHadHad->Write("dPhiDEtaHadHad");
    dPhiDEtaLamLam->Write("dPhiDEtaLamLam");
    dPhiDEtaHadLam->Write("dPhiDEtaHadLam");
    dPhiDEtaLamDivHad->Write("dPhiDEtaLamDivHad");
    lamWithAntiLamCounter->Write("lamWithAntiLamCounter");

    pdgDistCutHadLam->Write("pdgDistCutHadLam");
    hadPtDistCut->Write("hadPtDistCut");
    lamPtDistCut->Write("lambdaPtDistCut");
    phiDistCutHadLam->Write("phiDistCutHadLam");
    etaDistCutHadLam->Write("etaDistCutHadLam");
    dPhiCutHadHad->Write("dPhiCutHadHad");
    dPhiCutLamLam->Write("dPhiCutLamLam");
    dPhiCutHadLam->Write("dPhiCutHadLam");
    dPhiCutLamDivHad->Write("dPhiCutLamDivHad");
    dPhiDEtaCutHadHad->Write("dPhiDEtaCutHadHad");
    dPhiDEtaCutLamLam->Write("dPhiDEtaCutLamLam");
    dPhiDEtaCutHadLam->Write("dPhiDEtaCutHadLam");
    dPhiDEtaCutLamDivHad->Write("dPhiDEtaCutLamDivHad");
    lamWithAntiLamCutCounter->Write("lamWithAntiLamCutCounter");

    //outFile->Close();
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


void readPythiaPrimary()
{
    TStopwatch timer;
    timer.Start();

    // load ROOT file generated with minbiasgen_fixed.c
    AliRunLoader* inputRun = AliRunLoader::Open("outloop_io_test_1mil.root");
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
        std::vector<TParticle*> hadron_list;
        std::vector<TParticle*> lambda_list;

        std::vector<TParticle*> hadron_trigger_list;
        std::vector<TParticle*> lambda_trigger_list;

        std::vector<TParticle*> hadron_assoc_list;
        std::vector<TParticle*> lambda_assoc_list;

        // for each event, load the stack (list of all particles in an event)
        inputRun->GetEvent(event);
        AliStack *theStack = inputRun->Stack();

        // loop through the stack.
        // GetNTracks() returns the length of our "list" (stack)
        for (int i = 0; i < theStack->GetNtrack(); i++)
            {

            TParticle *particle = theStack->Particle(i);
            //AliMCParticle *particle = theStack->Particle(i);
            
            bool primary = theStack->IsPhysicalPrimary(i);                

            int pdg = particle->GetPdgCode();
            double pT = particle->Pt();
            double phi = particle->Phi();
            double eta = particle->Eta();
            

            // fill histograms for all particles
            pdgDist->Fill(pdg);
            ptDist->Fill(pT);
            phiDist->Fill(phi);
            etaDist->Fill(eta);

            // select for PRIMARY HADRONS
            // protons, charged kaons, charged pions, electrons, and their corresponding antiparticles
            // append to list
            if ( (primary) && ( (pdg==2212) || (pdg==-2212) || (pdg==321) || (pdg==-321) || (pdg==211) || (pdg==-211) || (pdg==11) || (pdg==-11) ) ) 
                {
                hadron_list.push_back(particle);

                pdgDistHadLam->Fill(pdg);
                hadPtDist->Fill(pT);
                phiDistHadLam->Fill(phi);
                etaDistHadLam->Fill(eta);

                if ((pT <= 8) && (pT >= 4)){
                    hadron_trigger_list.push_back(particle);
                }
                else if ((pT <= 4) && (pT >= 2)){
                    hadron_assoc_list.push_back(particle);
                }

                }

            // select for PRIMARY lambdas, append to list
            else if ( (primary) && ( (pdg==3122) || (pdg==-3122) ) )
                {
                    lambda_list.push_back(particle);

                    pdgDistHadLam->Fill(pdg);
                    lamPtDist->Fill(pT);
                    phiDistHadLam->Fill(phi);
                    etaDistHadLam->Fill(eta);

                if ((pT <= 8) && (pT >= 4)){
                    lambda_trigger_list.push_back(particle);
                }
                else if ((pT <= 4) && (pT >= 2)){
                    lambda_assoc_list.push_back(particle);
                }
                } 
                
                
            }


        int numHadrons = hadron_list.size(); 
        int numLambdas = lambda_list.size();
        int numHadTriggers = hadron_trigger_list.size();
        int numLamTriggers = lambda_trigger_list.size();
        int numHadAssoc = hadron_assoc_list.size();
        int numLamAssoc = lambda_assoc_list.size();

        // Had-Had Correlations
        for (int i=0; i < numHadTriggers; i++)
            {
            TParticle *trigger = hadron_trigger_list[i];

            double trig_pt = trigger->Pt();
            double trig_phi = trigger->Phi();
            double trig_eta = trigger->Eta();

            if ((TMath::Abs(trig_eta) < 0.8))
            {
            pdgDistHadLam->Fill(trigger->GetPdgCode());  
            }

            
            for (int j=0; j < numHadAssoc; j++)
                {

                TParticle *assoc = hadron_assoc_list[j];

                double assoc_pt = assoc->Pt();
                double assoc_phi = assoc->Phi();
                double assoc_eta = assoc->Eta();

                double delta_phi = caldeltaphi(trig_phi, assoc_phi);
                double delta_eta = caldeltaeta(trig_eta, assoc_eta);

                dPhiHadHad->Fill(delta_phi);
                dPhiDEtaHadHad->Fill(delta_phi, delta_eta);

                if ((TMath::Abs(trig_eta) < 0.8) && (TMath::Abs(assoc_eta) < 0.8))
                {
                    etaDistCutHadLam->Fill(trig_eta);
                    phiDistCutHadLam->Fill(trig_phi);
                    pdgDistCutHadLam->Fill(trigger->GetPdgCode());
                    
                    dPhiCutHadHad->Fill(delta_phi);
                    dPhiDEtaCutHadHad->Fill(delta_phi, delta_eta);
                }
                }
            }

            
        
        // Lam-Lam Correlations
        for (int i=0; i < numLamTriggers; i++)
            {
            TParticle *trigger = lambda_trigger_list[i];

            int trig_pdg = trigger->GetPdgCode();
            double trig_pt = trigger->Pt();
            double trig_phi = trigger->Phi();
            double trig_eta = trigger->Eta();

            if ((TMath::Abs(trig_eta) < 0.8))
            {
            pdgDistHadLam->Fill(trigger->GetPdgCode());  
            }

            
            for (int j=0; j < numLamAssoc; j++)
                {

                TParticle *assoc = lambda_assoc_list[j];

                int assoc_pdg = assoc->GetPdgCode();
                int pdg_ratio = trig_pdg / assoc_pdg;
                double assoc_pt = assoc->Pt();
                double assoc_phi = assoc->Phi();
                double assoc_eta = assoc->Eta();

                double delta_phi = caldeltaphi(trig_phi, assoc_phi);
                double delta_eta = caldeltaeta(trig_eta, assoc_eta);

                dPhiLamLam->Fill(delta_phi);
                dPhiDEtaLamLam->Fill(delta_phi, delta_eta);

                if (pdg_ratio == -1)
                {
                lamWithAntiLamCounter->Fill(1);
                }
                else 
                {
                lamWithAntiLamCounter->Fill(0);
                }

                if ((TMath::Abs(trig_eta) < 0.8) && (TMath::Abs(assoc_eta) < 0.8))
                {
                    etaDistCutHadLam->Fill(trig_eta);
                    phiDistCutHadLam->Fill(trig_phi);
                    pdgDistCutHadLam->Fill(trigger->GetPdgCode());
                    
                    dPhiCutLamLam->Fill(delta_phi);
                    dPhiDEtaCutLamLam->Fill(delta_phi, delta_eta);

                    if (pdg_ratio == -1)
                    {
                    lamWithAntiLamCutCounter->Fill(1);
                    }
                    else 
                    {
                    lamWithAntiLamCutCounter->Fill(0);
                    }
                }
                }
            } 

        // Had-Lam Correlations    
        for (int i=0; i < numHadTriggers; i++)
            {
            TParticle *trigger = hadron_trigger_list[i];

            double trig_pt = trigger->Pt();
            double trig_phi = trigger->Phi();
            double trig_eta = trigger->Eta();

            for (int j=0; j < numLamAssoc; j++)
                {

                TParticle *assoc = lambda_assoc_list[j];

                double assoc_pt = assoc->Pt();
                double assoc_phi = assoc->Phi();
                double assoc_eta = assoc->Eta();

                double delta_phi = caldeltaphi(trig_phi, assoc_phi);
                double delta_eta = caldeltaeta(trig_eta, assoc_eta);

                dPhiHadLam->Fill(delta_phi);
                dPhiDEtaHadLam->Fill(delta_phi, delta_eta);

                if ((TMath::Abs(trig_eta) < 0.8) && (TMath::Abs(assoc_eta) < 0.8))
                {   
                    dPhiCutHadLam->Fill(delta_phi);
                    dPhiDEtaCutHadLam->Fill(delta_phi, delta_eta);
                }
                }
            }                   
            
        
        }

    
    *dPhiLamDivHad = (*dPhiLamLam) / (*dPhiHadHad);
    *dPhiDEtaLamDivHad = (*dPhiDEtaLamLam) / (*dPhiDEtaHadHad);

    *dPhiCutLamDivHad = (*dPhiCutLamLam) / (*dPhiCutHadHad);
    *dPhiDEtaCutLamDivHad = (*dPhiDEtaCutLamLam) / (*dPhiDEtaCutHadHad);

    /*
    float bin1 = dPhiCutHadHad->FindBin(3.14159265/2);
    float bin2 = dPhiCutHadHad->FindBin(3*3.14159265/2);

    std::cout << dPhiCutHadHad->Integral(bin1, bin2)/dPhiCut
    */

    writeOut();

    timer.Stop();
    timer.Print();

}


