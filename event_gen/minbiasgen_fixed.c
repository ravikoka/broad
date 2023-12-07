R__LOAD_LIBRARY(libpythia6.so)
#include "AliHeader.h"
#include "AliGenPythia.h"
#define LOG_NO_INFO

void minbiasgen_fixed(Int_t nev = 10000000, char* filename = (char*) "outloop_io_test_1mil.root", char* folderName = (char*) "test"){

    TStopwatch timer;
    timer.Start();

    gRandom->SetSeed(0);//put 0 to use system time
    AliRunLoader* rl = AliRunLoader::Open(filename, folderName, "recreate");

    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(10000);
    rl->LoadKinematics("RECREATE");
    rl->MakeTree("E");
    gAlice->SetRunLoader(rl);

    //Create stack
    rl->MakeStack();
    AliStack* stack = rl->Stack();

    //Header
    AliHeader* header = rl->GetHeader();

    AliGenPythia * gener = new AliGenPythia(-1);   //something like event scaling
    gener->SetEnergyCMS(7000.);
    gener->SetPtHard(4,10000);
    gener->SetOrigin(0.,0.,0.); 
    gener->SetSigma(0.,0.,5.3);
    gener->SetVertexSmear(kPerEvent);
    gener->SetTriggerChargedMultiplicity(2, 0.9, 0.150); //only events with 2 charged particles, in eta range <0.9, >150 MeV
    gener->SetTriggerParticle(0, 1., 5., 1000.);
    gener->SetChildPtRange(5, 1000);
    gener->SetYRange(-2, 2);

    //Initialize generator
    gener->Init();
    gener->SetStack(stack);

    //Event Loop
    for (int iev = 0; iev < nev; iev++) {
        if (iev % 500 == 0)  Printf("%s", Form(" Event number %d", iev));
        //Initialize event
        header->Reset(0,iev);
        rl->SetEventNumber(iev);
        stack->Reset();
        rl->MakeTree("K");
        rl->MakeStack();
        stack->ConnectTree(rl->TreeK());

        //Generate event
        gener->Generate();

        //Analysis
        Int_t npart = stack->GetNprimary();

        //Finish event
        header->SetNprimary(stack->GetNprimary());
        header->SetNtrack(stack->GetNtrack());

        //I/O
        stack->FinishEvent();
        header->SetStack(stack);
        rl->TreeE()->Fill();
        rl->WriteKinematics("OVERWRITE");

    }   //end event loop

    //Termination

    //Generator
    gener->FinishRun();

    //Write file
    rl->WriteHeader("OVERWRITE");
    gener->Write();
    rl->Write();

    rl->RemoveEventFolder();
    rl->UnloadAll();
  
    timer.Stop();
    timer.Print();
}
