#include <iostream>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

using namespace std;

void dat2root(const char* dir, const char* filename, unsigned int fNcycles=100000000) {
 
    bool invert0 = false;
    unsigned int nsamples = 1024; 
    unsigned int offset = 830; // FIXME
    unsigned int adc=0;
    unsigned int ch_n=0;
    unsigned int evt_n=0;
    unsigned long t_stmp=0;
    double Vpp=10.;
    double adc2v=Vpp/TMath::Power(2.,14);
    double sample2mu=0.01;

    double voltage0;
    unsigned long n(0);

    unsigned int ncycles = fNcycles;

    FILE* fp0 = fopen(Form("raw/%s/%s",dir,filename),"rb");
    if(fp0==NULL){
        cout << "Error open dat file:" << Form("raw/%s/%s",dir,filename) << endl;
        exit(1);
    }
    fseek(fp0, 0, SEEK_END);
    int filesize = ftell(fp0);
    cout << "file size " << filesize/1e6 << " MB" << endl;
    fseek(fp0, 0, SEEK_SET);
    
    gSystem->Exec(Form("mkdir -p root/%s", dir));
    TFile* out = new TFile(Form("root/%s/%s.root",dir,filename),"RECREATE");

    double mus=0;
    unsigned int sample=0; 

    TTree* t0 = new TTree("conditions","conditions record");
    t0->Branch("invert0",&invert0);
    t0->Branch("offset",&offset);
    t0->Branch("sample2mu",&sample2mu);
    t0->Branch("nsamples",&nsamples);
    t0->Branch("ncycles",&ncycles);
    t0->Branch("adc2v",&adc2v);
    t0->Fill();

    TTree* t = new TTree("tree","tree");
    t->Branch("sample",&sample);
    t->Branch("voltage",&voltage0);
    t->Branch("mus",&mus);
    t->Branch("evt_n",&evt_n);

    bool isEOF = false;

    while(sample<nsamples*(ncycles)) {
        for(unsigned int i=0;i<nsamples;++i) {
            fread((unsigned char*)&adc,2,1,fp0);
            if(feof(fp0)) {
              isEOF=true;
              break;
            };
            if(sample!=0 && sample%nsamples==0) {
                //cout<<"Reading evt_n="<<evt_n<<endl;
                evt_n++;
            }
            voltage0 = (adc*1.-offset)*adc2v;
            mus = sample2mu*sample - nsamples*evt_n*sample2mu;
            t->Fill();
            sample++;
        }
        if((sample/nsamples)%(filesize/nsamples/2/10)==0){
            printf("Converting to root: %.0f %% done.\r", sample*2*100./filesize);
            fflush(stdout);
        }
        if(isEOF) break;
    }

    fclose(fp0);

    t0->Write();
    t->Write();
    out->Close();
    return;

}

int main(int argc, char** argv){

    if(argc!=2 && argc!=3){
        cout << "syntax: " << argv[0] << " <input dir> [input filename]" << endl;
        exit(1);
    }

    const char* infile = argv[1];
    const char* outfile = "wave0.dat";
    if(argc==3) outfile = argv[2];

    dat2root(infile, outfile);

}
