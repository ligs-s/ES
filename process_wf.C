bool ENABLE_PEAK_FINDING = true;
bool JUST_DRAW = false;

void process_wf(TString filename,
                TString outfilename,
                int nsamples=1024,
                int sig_start=350,
                int sig_end=500){
    
    if(ENABLE_PEAK_FINDING){
        cout << "******* PEAK FINDING IS ON *******" << endl;
    }
    if(JUST_DRAW){
        cout << "******* PLOTTING WAVEFORM ONLY *******" << endl;
    }

    bool isPMT = false;
    if ((filename.Contains("PMT") && filename.Contains("shaped") && !filename.Contains("inverted")) || filename.Contains("pmt")){
        cout << "PMT shaped signal" << endl;
        isPMT = true;
    }
        

    TFile* file = new TFile(filename, "read");
    TTree* tree = (TTree*)file->Get("tree");
    long int nEntries = tree->GetEntries();
    cout << "Total entries: " << nEntries << endl;

    int pos = outfilename.Last('/');
    TString outdir = filename(0, pos);
    //cout << outdir.Data() << endl;

    unsigned int evt_n;
    double v;
    tree->SetBranchAddress("evt_n",     &evt_n);
    tree->SetBranchAddress("voltage",   &v);

    //int ped_end = TMath::Min(int(nsamples*0.2), sig_start-20);
    //int post_ped_start = TMath::Max(int(nsamples*0.8), sig_end+20);
    int ped_end = sig_start-20;
    int post_ped_start = sig_end+20;
    int pre_peak = 10;
    int post_peak = 10;

    TFile* outfile = new TFile(outfilename, "recreate");
    TTree* cond_tree = new TTree("condition", "condition");
    cond_tree->Branch("nsamples", &nsamples, "nsamples/I");
    cond_tree->Branch("sig_start", &sig_start, "sig_start/I");
    cond_tree->Branch("sig_end", &sig_end, "sig_end/I");
    cond_tree->Branch("ped_end", &ped_end, "ped_end/I");
    cond_tree->Branch("post_ped_start", &post_ped_start, "post_ped_start/I");
    cond_tree->Branch("pre_peak", &pre_peak, "pre_peak/I");
    cond_tree->Branch("post_peak", &post_peak, "post_peak/I");
    cond_tree->Fill();


    TTree* out_tree = new TTree("wf_info", "waveform info");
    double info[18];
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
    out_tree->Branch("charge_near_peak",    &info[10],   "charge_near_peak/D"); // charge int around peak pos
    out_tree->Branch("peak_pos_local",      &info[11],   "peak_pos_local/D"); // charge int around peak pos
    out_tree->Branch("height_local",        &info[12],   "height_local/D"); // charge int around peak pos
    out_tree->Branch("min_pos",             &info[13],   "min_pos/D"); 
    out_tree->Branch("min_height",          &info[14],   "min_height/D"); 
    out_tree->Branch("charge_near_peak_local",    &info[15],   "charge_near_peak_local/D"); // charge int around peak pos
    out_tree->Branch("height_left_local",        &info[16],   "height_left_local/D");
    out_tree->Branch("height_right_local",        &info[17],   "height_right_local/D");

    // peak_finding
    const int MAXPEAKS = 30;
    int pf_npeaks;
    float pf_peak_pos[MAXPEAKS];
    float pf_height[MAXPEAKS];
    if(ENABLE_PEAK_FINDING){
        out_tree->Branch("pf_npeaks",       &pf_npeaks,     "pf_npeaks/I");
        out_tree->Branch("pf_peak_pos",     pf_peak_pos,    "pf_peak_pos[pf_npeaks]/F");
        out_tree->Branch("pf_height",       pf_height,      "pf_height[pf_npeaks]/F");
    }

    double* vs = new double[nsamples];
    double* xs = new double[nsamples];

    int last_evt_n = 0;
    int sample_idx = 0;

    // charge
    TH1F* h = new TH1F("h", "charge", 2000, -5, 15);

    TCanvas* c1 = new TCanvas();
    const int n_plotted_wfs = 100;
    bool pdf_saved = false;
    c1->Print(outdir+"/waveform.pdf[");

    for(int i=0;i<nEntries;i++){
        if(i%(nEntries/10)==0){
            printf("Progress: %.0f %% finished.\r", i*100.0/nEntries);
            fflush(stdout);
        }
        tree->GetEntry(i);
        
        if(evt_n!=last_evt_n){
            //cout << "Event #" << evt_n << endl;
            // next waveform
            TGraph* gr = new TGraph(nsamples, xs, vs);
            if(JUST_DRAW){
                if(last_evt_n==250){
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerSize(0.5);
                    gr->Draw("ALp");
                    return;
                }
            }
            else{

                if(evt_n<n_plotted_wfs){
                    // draw only the first 20 waveforms
                    gr->Draw("Al");
                    c1->Print(outdir+"/waveform.pdf");
                }

                // process waveform and extract information
                analyze_wf_2(gr, info, sig_start, sig_end, ped_end, post_ped_start, pre_peak, post_peak, isPMT);
                
                if(ENABLE_PEAK_FINDING){
                    TH1F* h_wf = new TH1F("h_wf", "h_wf", nsamples, -0.5, nsamples-0.5);
                    for(int ii=0;ii<nsamples;ii++) h_wf->SetBinContent(ii+1, vs[ii]);
                    peak_finding(h_wf, info, MAXPEAKS, pf_npeaks, pf_peak_pos, pf_height); // peak_finding
                    //for(int ii=0;ii<pf_npeaks;ii++){
                    //    cout << ii << " " << pf_peak_pos[ii] << " " << pf_height[ii] << endl;
                    //
                    
                    //if(evt_n<n_plotted_wfs){
                    //    // draw only the first 20 waveforms
                    //    h_wf->Draw();
                    //    c1->Print(outdir+"/waveform.pdf");
                    //}

                    delete h_wf;
                }

                out_tree->Fill();

                h->Fill(info[0]);
                if(sample_idx!=nsamples){
                    cout << "Error: more than " << nsamples 
                         << " samples in one waveform found!" << endl;
                }
            }
            delete gr;
            sample_idx = 0;
        }
        vs[sample_idx] = v;
        xs[sample_idx] = sample_idx;
        sample_idx++;
        last_evt_n = evt_n;
        if(evt_n==n_plotted_wfs && !pdf_saved){
            c1->Print(outdir+"/waveform.pdf]");
            pdf_saved = true;
        }
    }


    //h->Draw();
    outfile->cd();
    cond_tree->Write();
    out_tree->Write();
    h->Write();
    outfile->Close();
    file->Close();

    delete vs;
    delete xs;

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

void analyze_wf(TGraph* gr, double* info, bool isPMT=false){

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
    //
    // latet set up
    const int ped_start = 0;
    const int ped_end = 300;

    const int sig_start = 350;
    const int sig_end = 500;

    const int post_ped_start = 500;
    const int post_ped_end = 1000;

    //// Ako's 201606 setup
    //const int ped_start = 0;
    //const int ped_end = 250;

    //const int sig_start = 300;
    //const int sig_end = 450;

    //const int post_ped_start = 460;
    //const int post_ped_end = 500;

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

    double ypeak;
    int xpeak;
    if(isPMT){
        ypeak = TMath::MinElement(npts, y);
        xpeak = TMath::LocMin(npts, y);
    }
    else{
        ypeak = TMath::MaxElement(npts, y);
        xpeak = TMath::LocMax(npts, y);
    }

    double height = ypeak - ped;

    double charge_near_peak = 0.;
    for(int i=xpeak-65;i<xpeak+85;i++){
        gr->GetPoint(i, x, v);
        charge_near_peak += (v-ped);
    }
    if(isPMT){
        charge *= -1;
        height *= -1;
        charge_near_peak *= -1;
    }
    
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
    info[9] = xpeak;
    info[10] = charge_near_peak;

    delete func1;

    return;

}

void analyze_wf_2(TGraph* gr, 
                  double* info, 
                  int sig_start,
                  int sig_end,
                  int ped_end,
                  int post_ped_start,
                  int pre_peak, // window defined for charge integration
                  int post_peak, // window defined for charge integration
                  bool isPMT=false){

    // improved reconstruction
    //
    int npts = gr->GetN();
    double* y = gr->GetY();

    

    // 256 sample
    const int ped_start = 0;
    //const int ped_end = TMath::Min(npts*0.2, sig_start-50);

    //const int sig_start = 120;
    //const int sig_end = 150;

    //const int post_ped_start = TMath::Max(npts*0.8, sig_end+50);
    const int post_ped_end = npts-1;

    //// latet set up
    //const int ped_start = 0;
    //const int ped_end = 300;

    //const int sig_start = 350;
    //const int sig_end = 500;

    //const int post_ped_start = 500;
    //const int post_ped_end = 1000;

    //const int pre_peak = 65;
    //const int post_peak = 85;
    
    //const int pre_peak = 10;
    //const int post_peak = 10;

    TF1* func1 = new TF1("func1","[0]", 0, npts);
    gr->Fit(func1, "Q", "", ped_start, ped_end);
    double ped = func1->GetParameter(0); // ped
    double ped_err = func1->GetParError(0);
    gr->Fit(func1, "Q", "", post_ped_start, post_ped_end);
    double post_ped = func1->GetParameter(0);

    for(int i=0;i<npts;i++){
        // remove baseline 
        y[i] -= ped;
        if(isPMT){
            // and invert waveform if it is PMT signal
            y[i] *= -1;
        }
    }

    // calculate useful info based on baseline subtracted waveform
    double x, v;
    double ped_mean, ped_rms;
    double sig_mean, sig_rms;
    get_stats(y, ped_start, ped_end, ped_mean, ped_rms);
    get_stats(y, sig_start, sig_end, sig_mean, sig_rms);

    double charge = 0;
    for(int i=sig_start;i<sig_end;i++){
        gr->GetPoint(i, x, v);
        charge += v;
    }

    double ypeak = TMath::MaxElement(npts, y);
    double xpeak = TMath::LocMax(npts, y);

    double ymin = TMath::MinElement(npts, y);
    double xmin = TMath::LocMin(npts, y);

    int sig_window = sig_end - sig_start;

    double ypeak_local = TMath::MaxElement(sig_window, &y[sig_start]);
    double xpeak_local = TMath::LocMax(sig_window, &y[sig_start]) + sig_start;

    // calculate local maximum at side band
    double ypeak_left_local = TMath::MaxElement(sig_window, &y[sig_start-20-sig_window]);
    double ypeak_right_local = TMath::MaxElement(sig_window, &y[sig_end+20]);

    double height = ypeak;
    double height_local = ypeak_local;

    double charge_near_peak = 0.;
    for(int i=xpeak-pre_peak;i<xpeak+post_peak;i++){
        if (i<0 || i>=npts) continue;
        gr->GetPoint(i, x, v);
        charge_near_peak += v;
    }

    double charge_near_peak_local = 0.;
    for(int i=xpeak_local-pre_peak;i<xpeak_local+post_peak;i++){
        if (i<0 || i>=npts) continue;
        gr->GetPoint(i, x, v);
        charge_near_peak_local += v;
    }

    info[0] = charge;
    info[1] = ped;
    info[2] = post_ped;
    info[3] = height;
    info[4] = ped_mean;
    info[5] = ped_rms;
    info[6] = ped_err;
    info[7] = sig_mean;
    info[8] = sig_rms;
    info[9] = xpeak;
    info[10] = charge_near_peak;
    info[11] = xpeak_local;
    info[12] = height_local;
    info[13] = xmin;
    info[14] = ymin;
    info[15] = charge_near_peak_local;
    info[16] = ypeak_left_local;
    info[17] = ypeak_right_local;

    delete func1;

}

void peak_finding(TH1* h, double* info, int max_peaks, int &npeaks, float* peak_pos, float* height){

    // peak finding with TSpectrum
    //TH1F* h_wf = new TH1F("h_wf", "h_wf", n, 0, n);
    //for(int i=0;i<n;i++){
    //    h_wf->SetBinContent(i+1, y[i]);
    //}
    //
    double ymax = info[3];
    double rms = info[5];
    double baseline = info[1];
    double thresh = 5*rms/ymax;
    //cout << "threshold: ymax=" << ymax << " rms=" << rms << " thresh" << thresh << endl;

    TSpectrum* sp = new TSpectrum();
    int nfound = sp->Search(h, 2, "", thresh);
    //sp->Print(); //print found peaks
    //h->Draw();
    //cout << "nfound " << nfound << endl;

    npeaks = sp->GetNPeaks();
    if (npeaks>max_peaks) npeaks = max_peaks;
    // only float works for this version of ROOT
    float* xx = sp->GetPositionX();
    float* yy = sp->GetPositionY();

    for(int i=0;i<npeaks;i++){
        peak_pos[i] = xx[i];
        height[i] = yy[i];
        //cout << i << " 1. " << xx[i] << " " << yy[i] << endl;
        //cout << i << " 2. " << info[9] << " " << info[3] << endl;
    }

    delete sp;
    
    // ========


}
