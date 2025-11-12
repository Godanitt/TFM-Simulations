#include <iostream>
#include <fstream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <time.h>
#include <TH2F.h>
#include <TRandom3.h>

// #include "Garfield/ViewField.hh"
// #include "Garfield/ViewCell.hh"
// #include "Garfield/ViewDrift.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"

#include <time.h>

using namespace Garfield;

using namespace std;



//Define a public tree that can be access in all functions and the global variables to be shared
double xExc,yExc,zExc,tExc;
int ie, eventNumber;

TTree *dataExc = new TTree("dataExc", "dataExc"); //THE BRANCHS ARE DEFINED LATTER


void userHandle(double x, double y, double z, double t, int type, int level,Medium * m)
{
	xExc=x;
	yExc=y;
	zExc=z;
	tExc=t;
	dataExc->Fill();
}


// the program must be called like:
// ./ fatGem "ansysfilesdir"
// in my case ./fatGem "/home/arouca/Desktop/ansys/cutansys_results/"
// ./fatGem ../fieldMaps2bar teste.root 1 5 10 xe 99. cf4 1.
int main(int argc, char *argv[]){

	if (argc<9){
		cout<<"Wrong number of arguments."<<endl;
		cout<<"./fatGem fieldMapFolder rootFileName.root pressure(bar) pitch(mm) npe(#) gas1() mixture(%) gas2() mixture(%)"<<endl;
		exit(EXIT_FAILURE);
	}

	time_t tstart, tend;

	double dif;

	time (&tstart);

	clock_t start, end = 0;

	double cpuTime = 0;

	start = clock ();

	// fundamental_constants fundks;


	//get simulation parameters:


	char gasstr[5];

	std::strcpy (gasstr, argv[6]);

	// Pasamos la presión a float
	TString pressureInputString=argv[3];
	double pressureBar = pressureInputString.Atof();
	double pressure= double(pressureBar)* 750.061683; //bar to torr
	
	// Pasamos la npe a int
	TString npeInputString=argv[5];
	int npe = npeInputString.Atoi();
	// Pasamos la mixture1 a float
	TString mixture1InputString=argv[7];
	// Pasamos gas1 a TString
	std::string gas1_str = argv[6];
	double mixture1 = mixture1InputString.Atof();
	// Pasamos gas1 a TString
	std::string gas2_str = argv[8];
	// Pasamos la mixture2 a float
	TString mixture2InputString=argv[9];
	double mixture2 = mixture2InputString.Atof();
    std::cout << "SetComposition: mezcla " << mixture2 << "\n";


	double temp = 293;

	char rootfile[200];


	std::strcpy (rootfile, argv[2]);//segundo argumento é o nome do ficheiro ROOT

	char elist_file[500];

	char nlist_file[500];

	char mplist_file[500];

	char prnsol_file[500];

	std::strcpy (elist_file, argv[1]);

	std::strcat (elist_file, "/ELIST.lis");

	std::strcpy (nlist_file, argv[1]);

	std::strcat (nlist_file, "/NLIST.lis");

	std::strcpy (mplist_file, argv[1]);

	std::strcat (mplist_file, "/MPLIST.lis");

	std::strcpy (prnsol_file, argv[1]);

	std::strcat (prnsol_file, "/PRNSOL.lis");

	ofstream outfile;

	// Make a gas medium
	MediumMagboltz * gas = new MediumMagboltz ();

	//const double lambdaPenning = 0.0;
	//gas->EnablePenningTransfer(70/100., lambdaPenning, "ar");
	//gas->LoadIonMobility("/home/pmendes/garf++/Data/IonMobility_Ne+_Ne.txt");
	//gas->EnableDeexcitation();
	// Set the temperature [K] and pressure [Torr]
	//MediumMagboltz gas("ar", 80., "ch4", 20.); 
	if (mixture2 <= 0.0) {
		gas->SetComposition(gas1_str, mixture1);  // solo un gas
		std::cout << "SetComposition: " << argv[6]
				<< " " << mixture1 << "% (gas único)\n";
	} else {
		gas->SetComposition(gas1_str, mixture1, gas2_str, mixture2);
		std::cout << "SetComposition: mezcla " << argv[6] << "/" << argv[8]
				<< " = " << mixture1 << "/" << mixture2 << "\n";
	}

	gas->EnableDebugging ();

	gas->Initialise ();

	gas->DisableDebugging ();

	gas->SetTemperature (temp);

	gas->SetPressure (pressure);

	//information about collison lelevs:
	int nLevels = gas->GetNumberOfLevels ();

	//std::cout<<"0"<<endl;  
	printf("# of levels: %i\n",nLevels);
	int ngas = 0;

	int type = 0;

	std::string description;

	double en = 0;

	//std::cout << "Information about collision levels: \n";
	for (int il = 0; il < nLevels; il++)
	{

		gas->GetLevel (il, ngas, type, description, en);

	} 
	//std::cout<<"1"<<endl;  
	// Make a component from ansys file
	ComponentAnsys123 * fm = new ComponentAnsys123 ();

	fm->Initialise (elist_file, nlist_file, mplist_file, prnsol_file, "mm");

	fm->EnableMirrorPeriodicityX ();

	fm->EnableMirrorPeriodicityY ();

	fm->PrintRange ();


	const int nMaterials = fm->GetNumberOfMaterials ();

	for (int i = 0; i < nMaterials; ++i)
	{

		const double eps = fm->GetPermittivity (i);

		cout << "permitivity: " << eps << "\n";

		if (eps == 1.)
			fm->SetMedium (i, gas);

	}

	fm->PrintMaterials ();


	TString pitchInputString=argv[4];//input im mm
	Int_t pitch_mm=pitchInputString.Atoi();
	const double pitch = double(pitch_mm)*0.1; //convert to mm

	cout<<"\n PITCH "<<pitch<<endl;

	double acrylic = 0.46;

	// Make a sensor
	Sensor * sensor = new Sensor ();

	sensor->AddComponent (fm);

	sensor->SetArea (-4 * pitch, -4 * pitch, -0.7317, 4 * pitch, 4 * pitch,1.7317);


	cout<<"\n Sensor is set \n"<<endl;

	//create a Tree file tree4.root

	TFile f (rootfile, "RECREATE");


	int nBinsGain = 10;

	int gmin = 0;

	int gmax=10;

	TH1F * hElectrons =
		new TH1F ("hElectrons", "Number of electrons", 
				nBinsGain, gmin, gmax);


	int nBinsChrg = 400;
	

	TH1F * hEstartz =
		new TH1F ("hEstart_position_z", "hEstart_z", 200, 0,1.7317);


	TH1F * hEstarty =
		new TH1F ("hEstart_position_y", "hEstart_y", 200,
				-3 * sqrt (3) * pitch / 2, 3 * sqrt (3) * pitch / 2);


	TH1F * hEstartx =
		new TH1F ("hEstart_position_x", "hEstart_x", 200, -3 * pitch / 2,
				3 * pitch / 2);


	TH1F * nend_elec_hist =
		new TH1F ("hChrgE_endz", "Charges_elec_endz", nBinsChrg * 100, -0.1,
				0.15);


	TRandom3 *rndistx = new TRandom3(0);
	TRandom3 *rndisty = new TRandom3(0);


	int ne, ni;

	double x0, y0, z0, t0, e0;

	double x1, y1, z1, t1, e1;

	int status;

	unsigned int nElastic, nIonising, nAttachment, nInelastic, nExcitation, nSuperelastic;

	// Create a ROOT Tree
	//TTree tree ("tree", "A Tree with Events");
	TTree *tree= new TTree("DataCharge","DataCharge");
	tree->Branch("eventNumber", &eventNumber,"eventNumber/I");  
	tree->Branch("PrimaryElectrons", &ne,"ne/I");  
	tree->Branch("ie", &ie,"ie/I");//Numero do electrao gerado 
	tree->Branch("nExc", &nExcitation, "nExc/I");
	tree->Branch("nIon", &nIonising, "nIon/I");
	tree->Branch("x0", &x0,"x0/D");  
	tree->Branch("y0", &y0,"y0/D");
	tree->Branch("z0", &z0,"z0/D");
	tree->Branch("t0", &t0,"t0/D");
	tree->Branch("e0", &e0,"e0/D");
	tree->Branch("x1", &x1,"x1/D");  
	tree->Branch("y1", &y1,"y1/D");
	tree->Branch("z1", &z1,"z1/D");
	tree->Branch("t1", &t1,"t1/D");
	tree->Branch("e1", &e1,"e1/D");
	tree->Branch("status", &status,"status/I");
	tree->Branch("pressure", &pressure,"pressure/D");
	tree->Branch("temp", &temp,"temp/D");

	dataExc->Branch("eventNumber", &eventNumber,"eventNumber/I"); 
	dataExc->Branch("ie", &ie,"ie/I");//Numero do electrao gerado 
	dataExc->Branch("xExc",&xExc,"xExc/D");
	dataExc->Branch("yExc",&yExc,"yExc/D");
	dataExc->Branch("zExc",&zExc,"zExc/D");
	dataExc->Branch("tExc",&tExc,"tExc/D");



	// Make a microscopic tracking class for electron transport
	AvalancheMicroscopic * aval = new AvalancheMicroscopic ();

	aval->SetSensor (sensor);

	aval->EnableAvalancheSizeLimit (1);	//limite no ganho: The idea gor this project is to not have gain


	////// Set User Handle to follow the electron collicion-by-collision
	aval->SetUserHandleInelastic(userHandle);

	//energy histogram:
	TH1F * histen = new TH1F ("hen", "energy distribution", 1000, 0.0, 100.0);


	aval->EnableElectronEnergyHistogramming (histen);



	Int_t nelec_total = 0;

	// Int_t ncollected = 0;

	// Int_t n_out = 0;


	// Calculate a few avalanches

	// int nAbort = 0;

	//double maxCPUTime = 86000;

	for (eventNumber = 0; eventNumber < npe; eventNumber++)
	{

		gas->ResetCollisionCounters();

		//if (cpuTime >= maxCPUTime)break;

		//distancias em cm
		x0 = -pitch+2*pitch*rndistx->Rndm();//to have a random distribution in the interval [-pitch,pitch]
		y0 = -pitch+2*pitch*rndisty->Rndm();
		//NOTE: Confirm if Z0 is inside sensor. otherwhise garfield will claim that ther is no valid eletric field at the launchin point  
		z0 = 1.5;	//larga os electroes a 8 mm do topo da fatGem. (1+0.2, ja que o referencial esta no meio do buraco)

		t0 = 0.0;

		if (eventNumber == 0){
			e0=0.1;}
		else{
				e0 = histen->GetRandom(); // generate energy randomly accordingly to previous collisions
			}				

	

		std::cout << eventNumber + 1 << " electrões num total de " << npe << endl;

		aval->AvalancheElectron (x0, y0, z0, t0, e0, 0., 0., 0.);


		std::cout << eventNumber + 1 << " finished " << npe << endl;

		// Get the number of electrons and ions in the avalanche.

		aval->GetAvalancheSize (ne, ni);

		//ep.ne = ne;

		nelec_total = nelec_total + ne;

		hElectrons->Fill (ne);

		for (ie = 0; ie < ne; ie++)
		{

			float timeused = (float) clock () / CLOCKS_PER_SEC;
			
			aval->GetElectronEndpoint (ie, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
			gas->GetNumberOfElectronCollisions(nElastic, nIonising, nAttachment, nInelastic, nExcitation, nSuperelastic);
			
			// Fill the tree
			tree->Fill ();

			cout << "Electron " << ie << " drifting time was " << (t1 - t0) <<" ns, with the status " << status << endl;

			//std::cout<<"electron "<<ie<<" xf "<< x1 <<" yf "<<y1 << " zf "<< z1<<" status= "<<status<<endl;
			//std::cout << "electron " << ie << " rf " << sqrt (pow (x1, 2) +  pow (y1,2)) << " zf " << z1 << " status= " << status << endl;

			hEstarty->Fill (y0);

			hEstartz->Fill (z0);

			hEstartx->Fill (x0);

			nend_elec_hist->Fill (z1);




			std::cout << "\nInformation about electron " << ie <<" from avalanche " << eventNumber << " :\n";

			std::cout << "(x0,y0,z0)= (" << x0 << "," << y0 << "," << z0 <<"), t0= " << t0 << ", e0= " << e0 << "\n";

			std::cout << "(x1,y1,z1)= (" << x1 << "," << y1 << "," << z1 <<	"), t1= " << t1 << ", e1= " << e1 << "\n";

			std::cout << "status = " << status << "\n\n";

			
			
		}
	}

	


	tree->Write ();

	dataExc->Write();

	histen->Write ();

	hElectrons->Write ();

	nend_elec_hist->Write ();

	hEstartz->Write ();

	hEstarty->Write ();

	hEstartx->Write ();


	printf ("average # of electrons produced: %g\n",
			((double) nelec_total) / ((double) npe));

	// //printf("average # of electrons on the acrylic: %g\n",((double)nacrylic)/((double) npe));
	// printf ("average # of electrons on the electrodes: %g\n",
	// 		((double) ncollected) / ((double) npe));

	// printf
	// 	("average # of primary electrons out of the zone of the central hole: %g\n",
	// 	 ((double) n_out) / ((double) npe));

	time (&tend);

	dif = difftime (tend, tstart);

	std::cout << "It took you " << dif << " seconds to finish.\n";


	end = clock ();

	cpuTime = (end - start) / (double) (CLOCKS_PER_SEC);

	printf ("CPU time = %gs\n", cpuTime);

	// printf ("Number of aborted avalanches: %d \n", nAbort);

	//return 0;

	delete hEstartx;

	delete hEstarty;

	delete hEstartz;

	delete hElectrons;


	delete nend_elec_hist;


	delete gas;

	delete fm;

	delete sensor;

	delete aval;

	delete histen;

	f.Close ();f.Close ();

	std::cout << "Finalizando correctamente fatGem()." << std::endl;
	return 0;


} 

