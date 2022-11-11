#include <algorithm> // for std::sort

void peak_finding(){

    const char* filename = "root/20180209/vuv4_2_47V_source/charge.root";
    TFile* f = new TFile(filename, "read");
    TH1F* h = (TH1F*)f->Get("h");
    double gain, ped;
    get_peaks(h);
    return;

}

void get_peaks(TH1* h, const char* outfilename="./tmp_peak_pos.txt"){
    // take a histogram as input
    // find peak positions with TSpectrum 
    // write found peaks to outfile

    TSpectrum* sp = new TSpectrum();
    int nfound = sp->Search(h, 2, "", 0.005);
    //sp->Print(); //print found peaks
    //h->Draw();
    cout << "nfound " << nfound << endl;

    int npeaks = sp->GetNPeaks();

    // only float works for this version of ROOT
    float *xx = sp->GetPositionX();
    //float *yy = sp->GetPositionY();

    //// sort array
    std::sort(xx, xx+npeaks);

    //// I don't know why this does not give the right fitted positions?!
    //for(int i=0;i<npeaks;i++){
    //    cout << i << " " << xx[i] << " " << yy[i] << endl;
    //}
    //
    FILE* pfile = fopen(outfilename, "w");

    TGraph *gr = new TGraph();
    for(int i=0;i<npeaks;i++){
        //cout << i << " " << x[i] << " " << y[i] << endl;
        gr->SetPoint(i, i, xx[i]);
        fprintf(pfile, "%f\n", xx[i]);
    }
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);

    fclose(pfile);

    //delete gr;
    
    return;

}

void get_gain_calib(const char* filename, double &gain, double &ped){

    TGraph* gr = new TGraph();

    ifstream infile(filename);
    if(!infile.is_open()) {
        cout << filename << " is not opened successfully, exit now ..." << endl;
        exit(1);
        //return;
    }
    TString name;
    int npt = 0;    
    while(true){
        string line;
        if(!infile.good()) break; // one more empty line will be loaded
        getline(infile, line);
        if(line.empty()) continue;
        
        istringstream sline(line);
        if(sline.peek()=='#'){
            //cout << "skip line " << line.data() << endl;
            continue;
        }
        TString tmp_line(line.data());
        strs = tmp_line.Tokenize(" ");
        int size = strs->GetEntries();
        double x, y;

        //cout << npt << " " << line.data() << endl;

        if(size==1){
            // only one value, the first value is by default 0 pe
            sline >> y;
            gr->SetPoint(npt, npt, y);
        }
        if(size==2){
            // two values format, first col is number of pe, sec col
            // is measured value
            sline >> x >> y;
            gr->SetPoint(npt, x, y);
        }
        npt++;
    }

    // fit gain
    new TCanvas();
    gr->Draw("AP");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.0);
    TF1* func = new TF1("func", "[0]+[1]*x", 0, 10);
    func->SetParameters(0, 1);
    gr->Fit(func, "Q");
    ped = func->GetParameter(0);
    gain = func->GetParameter(1);
    cout << "ped= " << ped << " gain=" << gain << endl;

}
