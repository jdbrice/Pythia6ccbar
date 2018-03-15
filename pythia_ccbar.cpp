
#include "TPythia6.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"


// STL
#include <iostream>
#include <string> 



using namespace std;


TPythia6 		*pythia = NULL;
TFile 			*tFile = nullptr;

int 			nParticles = 0;

#define PID_c 4
#define PID_mu 13
#define PID_string 92

#define MASS_mu 0.1056583745 // GeV/c^2

#define K_ID 2
#define K_PARENT 3



TH1 * hPassCuts;
TH2 * hFullPS, *hRapCut, *hDEtaCut, *hDEtaPtCut;
TH2 * hPt1Pt2, *hEta1Eta2, *hPhi1Phi2;

void makeHistos(){
	hPassCuts = new TH1D( "hPassCuts", "", 10, 0, 10 );
	hFullPS = new TH2D( "hFullPS", "Full PS", 500, 0, 5, 150, 0, 15 );
	hRapCut = new TH2D( "hRapCut", "Pair y cut", 500, 0, 5, 150, 0, 15 );
	hDEtaCut = new TH2D( "hDEtaCut", "Daughter eta cut", 500, 0, 5, 150, 0, 15 );
	hDEtaPtCut = new TH2D( "hDEtaPtCut", "Daughter eta cut", 500, 0, 5, 150, 0, 15 );

	hPt1Pt2 = new TH2D( "hPt1Pt2", "pt1 vs pt2", 150, 0, 15, 150, 0, 15 );
	hEta1Eta2 = new TH2D( "hEta1Eta2", "eta1 vs et2", 120, -0.6, 0.6, 120, -0.6, 0.6 );
	hPhi1Phi2 = new TH2D( "hPhi1Phi2", "phi1 vs. phi2", 300, -3.2, 3.2, 300, -3.2, 3.2 );
}

void setupPythia( int trigger, long int iseed = 0 ){
	pythia = new TPythia6();

	const int MSEL_MINBIAS = 1;
	const int MSEL_CCBAR_TRIG = 4;
	if ( MSEL_MINBIAS != trigger && MSEL_CCBAR_TRIG != trigger ){
		cout << "WARNING: trigger must be (1=mb) or (4=ccbar), got " << trigger << endl;
	}

	pythia->SetMSEL(trigger);

	// Tune from JOEY BUTTERWORTH
	pythia->SetPARP(91,1.0); //<kT>
	pythia->SetPARP(67,1.0);  //mstp*4


	// // TURN ON relavent ccbar decays
	for(int i=673; i<=683; i++) pythia->SetMDME(i,1,0); //D+ ->e
	for(int i=684; i<=694; i++) pythia->SetMDME(i,1,1); //D+ -> mu
	for(int i=695; i<=735; i++) pythia->SetMDME(i,1,0); //D+ -> other

	// 	// D0
	for(int i=747; i<=754; i++) pythia->SetMDME(i,1,0); //D0 -> e
	for(int i=755; i<=762; i++) pythia->SetMDME(i,1,1); //D0 -> mu
	for(int i=763; i<=807; i++) pythia->SetMDME(i,1,0); //D0 -> other

	for(int i=818; i<=823; i++) pythia->SetMDME(i,1,0); //Ds -> e
	for(int i=824; i<=828; i++) pythia->SetMDME(i,1,1); //Ds -> mu
	for(int i=829; i<=850; i++) pythia->SetMDME(i,1,0); //Ds -> other

	for(int i=857; i<=862; i++) pythia->SetMDME(i,1,0); //eta_c,J/psi,chi_2c
	for(int i=1090; i<=1096; i++) pythia->SetMDME(i,1,0); //Lambda_c -> e
	for(int i=1097; i<=1103; i++) pythia->SetMDME(i,1,1); //Lambda_c -> mu
	for(int i=1104; i<=1165; i++) pythia->SetMDME(i,1,0); //Lambda_c -> other
	

	pythia->SetMRPY(1, iseed); 
	pythia->Initialize("CMS","p","p",200);

	pythia->Pylist(12);
	pythia->Pystat(1);
	pythia->Pystat(4);
	pythia->Pystat(5);


	string name = "./pythia_ccbar_TRIG_" + to_string( trigger ) + "_seed_" + to_string( iseed ) +".root" ;
	tFile = new TFile( name.c_str(), "RECREATE" );
	tFile->cd();
}

void printPlc( int i ){
	cout << "K(1)=" << pythia->GetK(i, 1) << ", K(2)=" << pythia->GetK(i, 2) << ", K(3)=" << pythia->GetK(i, 3) << ", K(4)=" << pythia->GetK(i, 4) << ", K(5)=" << pythia->GetK(i, 5) << endl;
}

int state( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 1 );
	return -1;
}
int plcId( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 2 );
	return -1;
}
int parentIndex( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 3 );
	return -1;
}
int posX( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 4 );
	return -1;
}
int posY( int i ){
	if ( i>0 && i<= nParticles )
		return pythia->GetK( i, 5 );
	return -1;
}
TLorentzVector lvec( int i ){
	TLorentzVector lv;
	if ( i>0 && i<= nParticles ){
		lv.SetPxPyPzE(  pythia->GetP(i,1),
						pythia->GetP(i,2),
						pythia->GetP(i,3),
						pythia->GetP(i,4));
	}
	return lv;
}


/* findStrings
 * counts the number of strings in the event
 * for ccbar we want exactly 2 strings to be consistent with dielectron papers
 * for leading order correlated pair
 */
int findStrings(){
	int nStrings = 0;
	for(Int_t i = 0; i < nParticles; i++){
		int id = abs(pythia->GetK(i+1, K_ID ) );
		if(id == PID_string ){  
			int parentIdIndex = pythia->GetK( i+1, K_PARENT);
			
			if(abs(pythia->GetK(parentIdIndex,2)) == PID_c){
				nStrings++;
			}//charm

		}//string
	}
	return nStrings;
}


/* isMuon
 * checks particle i to see if it is decayed muon from some charmed meson
 */
bool isMuon( int i ){

	if ( abs( plcId( i ) ) == PID_mu ){

		if ( (posX(i) != 0 || posY(i) != 0) && state( i ) != 1 ){
			cout << "ERROR state pos mismatch" << endl;
		}

		int pIndex = parentIndex( i );
		if ( pIndex >= 1 ){
			// printPlc( pIndex );

			int parentId = abs(plcId( pIndex ));
			if(parentId==411||parentId==421||parentId==431||parentId==4122){
				// printPlc( i );
				// cout << "\tPARENT: "; printPlc( pIndex );
				return true;
			}

		} else {
			cout << "\t reject no parent" << endl;
		}
	}
	return false;
}

void findMuons(){
	
	TLorentzVector plv, nlv, lv;
	int pParentId = -1, nParentId = -1;
	bool foundPos = false; 
	bool foundNeg = false;

	for(Int_t i = 0; i < nParticles; i++){
		int iPlc = i+1;
		bool isMu = isMuon( iPlc );
		int pId = plcId( iPlc );
		if (  isMu && pId == PID_mu ){
			// set pos
			nlv = lvec( iPlc );
			nParentId = plcId( parentIndex( iPlc ) );
			foundNeg = true;
		} else if ( isMu && pId == -PID_mu ){
			plv = lvec( iPlc );
			pParentId = plcId( parentIndex( iPlc ) );
			foundPos = true;
		}

		if ( foundPos && foundNeg ) break;

	}

	hPassCuts->Fill( 0.5 );

	if ( foundPos && foundNeg ){
		lv = plv + nlv;

		// full phase space
		hPassCuts->Fill( 1.5 );
		hFullPS->Fill( lv.M(), lv.Pt() );

		if ( abs(lv.Rapidity()) < 0.5 ){
			hPassCuts->Fill( 2.5 );
			hRapCut->Fill( lv.M(), lv.Pt() );

			hPt1Pt2->Fill( plv.Pt(), nlv.Pt() );
			hEta1Eta2->Fill( plv.PseudoRapidity(), nlv.PseudoRapidity() );
			hPhi1Phi2->Fill( plv.Phi(), nlv.Phi() );

			if ( abs(plv.PseudoRapidity() ) < 0.5 && abs(nlv.PseudoRapidity() ) < 0.5 ){
				hPassCuts->Fill( 3.5 );
				hDEtaCut->Fill( lv.M(), lv.Pt() );

				if ( plv.Pt() > 1.0 && nlv.Pt() > 1.0 ){
					hPassCuts->Fill( 4.5 );
					hDEtaPtCut->Fill( lv.M(), lv.Pt() );
				} // pt > 1.0
			} // |eta| < 0.5
		} // pair |y| < 0.5
	}
}


void genEvents( ULong_t _nEvents ){

	ULong_t iEvent = 0;
	while ( iEvent < _nEvents ){

		pythia->GenerateEvent();
		nParticles = pythia->GetNumberOfParticles();

		int nStrings = findStrings();
		if ( 2 == nStrings )
			findMuons( );

		if ( iEvent % 1000 == 0 )
			cout << "." << std::flush;

		iEvent++;
	}
	cout << endl;

}


// NOTE: at first i thought it looked wrong to assume only 1 ccbar -> mu mu but it is so rare that this is fine

Int_t main( Int_t argc, Char_t **argv){
	cout << "argc = " << argc << endl;
	if ( argc < 4 ) {
		cout  << "USAGE\n GENERATOR trig(1=mb,4=ccbar) nEvents rndSeed" << endl;
		return 1;
	}

	int trigger = atoi( argv[1] );
	int nEvents = atoi( argv[2] );
	long int seed = atol( argv[3] );
	setupPythia( trigger, seed );
	makeHistos();
	cout << "Generating " << nEvents << " events" << endl;
	genEvents( nEvents );

	cout<<"::::::::Printing out stats:::::::"<<endl;
	pythia->Pystat(1);
	pythia->Pystat(4);
	pythia->Pystat(5);

	cout << "SAVE THEN END" << endl;
	
	tFile->Write();
	tFile->Close();
	delete tFile;
	cout << "END" << endl;



	return 0;
}
