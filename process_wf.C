void process_wf(TString filename,
                TString outfilename){

    TFile* file = new TFile(filename, "read");
    TTree* tree = (TTree*)file->Get("tree");
    int nEntries = tree->GetEntries();
    cout << "Total entries: " << nEntries << endl;

    int pos = outfilename.Last('/');
    TString outdir = filename(0, pos);
    //cout << outdir.Data() << endl;

    unsigned int evt_n;
    double v;
    tree->SetBranchAddress("evt_n",     &evt_n);
    tree->SetBranchAddress("voltage",   &v);

    TFile* outfile = new TFile(outfilename, "recreate");
    TTree* out_tree = new TTree("wf_info", "waveform info");
    double info[10];
    out_tree->Branch("charge",  &info[0],   "charge/D");
    out_tree->Branch("ped",     &info[1],   "ped/D");
    out_tree->Branch("post_ped",&info[2],   "post_ped/D");
    out_tree->Branch("height",  &info[3],   "height/D");
    out_tree->Branch("ped_mean",&info[4],   "ped_mean/D");
    out_tree->Branch("ped_rms", &info[5],   "ped_rms/D");
    out_tree->Branch("ped_err", &info[6],   "ped_err/D");
    out_tree->Branch("sig_mean",&info[7],   "sig_mean/D");
    out_tree->Branch("sig_rms", &info[8],   "sig_rms/D");
    out_tree->Branch("peak_pos",&info[9],   "peak_pos/D");

    const int nsamples = 1024;
    double vs[nsamples];
    double xs[nsamples];

    int last_evt_n = 0;
    int sample_idx = 0;

    // charge
    TH1F* h = new TH1F("h", "charge", 2000, -5, 15);

    TCanvas* c1 = new TCanvas();
    c1->Print(outdir+"/waveform.pdf[");

    for(int i=0;i<nEntries;i++){
        if(i%(nEntries/10)==0){
            printf("Progress: %.0f %% finished.\r", i*100.0/nEntries);
            fflush(stdout);
        }
        tree->GetEntry(i);
        
        if(evt_n!=last_evt_n){
            // next waveform
            TGraph* gr = new TGraph(nsamples, xs, vs);
            if(evt_n<100){
                // draw only the first 20 waveforms
                gr->Draw("Al");
                c1->Print(outdir+"/waveform.pdf");
            }
            //gr->SetMarkerStyle(20);
            //gr->SetMarkerSize(0.5);
            //gr->Draw("AL");
            analyze_wf(gr, info);
            out_tree->Fill();
            delete gr;
            h->Fill(info[0]);
            if(sample_idx!=nsamples){
                cout << "Error: more than " << nsamples 
                     << " samples in one waveform found!" << endl;
            }
            sample_idx = 0;
        }
        vs[sample_idx] = v;
        xs[sample_idx] = sample_idx;
        sample_idx++;
        last_evt_n = evt_n;
    }

    c1->Print(outdir+"/waveform.pdf]");

    //h->Draw();
    outfile->cd();
    out_tree->Write();
    h->Write();
    outfile->Close();
    file->Close();

}

void get_stats(double* v, 
               int start_idx, 
               int end_idx, 
               double &mean,
               double &rms){

    int n = end_idx - start_idx;
    mean = 0.;
    for(int i=start_idx;i<end_idx;i++){
        mean += v[i];
    }
    mean /= n;

    rms = 0.;
    for(int i=start_idx;i<end_idx;i++){
        rms += (v[i]-mean)*(v[i]-mean);
    }
    rms = TMath::Sqrt(rms/n);
    return;


}

void analyze_wf(TGraph* gr, double* info){

    // flat line fit to baseline
    int npts = gr->GetN();
    double* y = gr->GetY();
    //cout << npts << endl;
    
    //const int ped_start = 0;
    //const int ped_end = 250;

    //const int sig_start = 250;
    //const int sig_end = 1000;

    //const int post_ped_start = 1000;
    //const int post_ped_end = 1000;

    //const int ped_start = 0;
    //const int ped_end = 600;

    //const int sig_start = 650;
    //const int sig_end = 1000;

    //const int post_ped_start = 1000;
    //const int post_ped_end = 1000;
    //
    const int ped_start = 0;
    const int ped_end = 300;

    const int sig_start = 350;
    const int sig_end = 450;

    const int post_ped_start = 500;
    const int post_ped_end = 1000;

    double x, v;
    double ped_mean, ped_rms;
    double sig_mean, sig_rms;
    get_stats(y, ped_start, ped_end, ped_mean, ped_rms);
    get_stats(y, sig_start, sig_end, sig_mean, sig_rms);

    TF1* func1 = new TF1("func1","[0]", 0, npts);
    gr->Fit(func1, "Q", "", ped_start, ped_end);
    double ped = func1->GetParameter(0);
    double ped_err = func1->GetParError(0);
    gr->Fit(func1, "Q", "", post_ped_start, post_ped_end);
    double post_ped = func1->GetParameter(0);
    double charge = 0;
    for(int i=sig_start;i<sig_end;i++){
        gr->GetPoint(i, x, v);
        charge += (v-ped);
    }
    double ymax = 0;
    int xmax = -1;
    for(int i=0;i<npts;i++){
        if(y[i]>ymax){
            ymax = y[i];
            xmax = i;
        }
    }
    double height = ymax - ped;
    
    //cout << charge << endl;
    info[0] = charge;
    info[1] = ped;
    info[2] = post_ped;
    info[3] = height;
    info[4] = ped_mean;
    info[5] = ped_rms;
    info[6] = ped_err;
    info[7] = sig_mean;
    info[8] = sig_rms;
    info[9] = xmax;

    delete func1;

    return;

}
