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

#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"

#include <time.h>

using namespace Garfield;
using namespace std;

//Define a public tree that can be access in all functions and the global variables to be shared
double xExc,yExc,zExc,tExc,levelExc;
int ie, eventNumber;

TTree *dataExc = new TTree("dataExc", "dataExc"); //THE BRANCHS ARE DEFINED LATTER


void userHandle(double x, double y, double z, double t, int type, int level,Medium * m)
{
	xExc=x;
	yExc=y;
	zExc=z;
	tExc=t;
	levelExc = level;
	dataExc->Fill();
}


int main(int argc, char *argv[]){


	///////////////////////////////////////////////////////////////////////////////////
	// Le pedimos valores de entrada al usuario

	if (argc<9){
		cout<<"Wrong number of arguments."<<endl;
		cout<<"./fatGem rootFileName.root fieldE(V/cm) pitch(mm) pressure(bar) npe(#) gas1() mixture(%) gas2() mixture(%)"<<endl;
		exit(EXIT_FAILURE);
	}

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Ahora recuperamos los datos de entrada y los hacemos datos reales de interés

	// Definimos temperatura: 
	double temp = 293.15;
	
	// Archivo de Root en el que queremos imprimir la información
	char rootfile[200];
	std::strcpy (rootfile, argv[1]);//segundo argumento é o nome do ficheiro ROOT


	// Pasamos el campo eléctrico a float
	TString uniformEInputStuniformEring=argv[2];
	double uniformE = uniformEInputStuniformEring.Atof();


	// Pasamos el pitch a float (pitch == gap, tamaño del cámpo eléctrico)
	TString pitchInputString=argv[3];//input im mm
	Int_t pitch_mm=pitchInputString.Atoi();
	const double pitch = double(pitch_mm)*0.1; //convert to mm


	// Pasamos la presión a float
	TString pressureInputString=argv[4];
	double pressureBar = pressureInputString.Atof();
	double pressure= double(pressureBar)* 750.061683; //bar to torr

	// Pasamos la npe a int
	TString npeInputString=argv[5];
	int npe = npeInputString.Atoi();

	// Pasamos gas1 a TString
	std::string gas1_str = argv[6];

	// Pasamos la mixture1 a float
	TString mixture1InputString=argv[7];
	double mixture1 = mixture1InputString.Atof();

	// Pasamos gas1 a TString
	std::string gas2_str = argv[8];

	// Pasamos la mixture2 a float
	TString mixture2InputString=argv[9];
	double mixture2 = mixture2InputString.Atof();

	
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Inicializamos el tiempo

	time_t tstart, tend;

	double dif;

	time (&tstart);

	clock_t start, end = 0;

	double cpuTime = 0;

	start = clock ();

	
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Inicializamos gas

	MediumMagboltz * gas = new MediumMagboltz ();

	// Detecta si tenemos un solo gas o mezcla de dos gases: 
	if (mixture2 <= 0.0) {
		gas->SetComposition(gas1_str, mixture1);  // solo un gas
		std::cout << "SetComposition: " << argv[6]
				<< " " << mixture1 << "% (gas único)\n";
	} else {
		gas->SetComposition(gas1_str, mixture1, gas2_str, mixture2);
		std::cout << "SetComposition: mezcla " << argv[6] << "/" << argv[8]
				<< " = " << mixture1 << "/" << mixture2 << "\n";
	}

	// Definimos temperatura
	gas->SetTemperature (temp);

	// Definimos presión
	gas->SetPressure (pressure);

	// Generamos las secciones eficaces de Magboltz (Steve) para generar el archivo
	gas->EnableDebugging ();
	gas->Initialise ();
	gas->DisableDebugging ();
	
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Leemos el número de niveles

	// Información sobre el número de niveles:
	int nLevels = gas->GetNumberOfLevels ();
	printf("# of levels: %i\n",nLevels);

	// Iniciamos variables que serán devueltas por gas->GetLevel
	int ngas = 0;
	int type = 0;
	std::string description;
	double en = 0;

	// Ahora leemos la información de cada nivel, recorriendo todos los posibles niveles
	for (int il = 0; il < nLevels; il++)
	{
		gas->GetLevel (il, ngas, type, description, en);
	} 
	

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Creamos el campo eléctrico

	ComponentUser* componentUniformElectricField = new ComponentUser();

	//componentUniformElectricField->EnableMirrorPeriodicityX ();
	//componentUniformElectricField->EnableMirrorPeriodicityY ();


	auto efield = [uniformE](const double x, const double y, const double z, double& ex, double& ey, double& ez) {
		ex = ey = 0;
		ez = uniformE;
	};

	componentUniformElectricField->SetElectricField(efield);
	componentUniformElectricField->SetArea(-4 * pitch, -4 * pitch, 0.0, 4 * pitch, 4 * pitch, pitch);
	componentUniformElectricField->SetMedium(gas);

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Creamos el Sensor
	// La clase Sensor es básicamente la compensación de componentes del detector

	Sensor * sensor = new Sensor ();

	// Añadimos el campo eléctrico
	sensor->AddComponent (componentUniformElectricField);

	// Añadimos el mismo area
	sensor->SetArea(-4* pitch, -4 * pitch, 0, 4 * pitch, 4 * pitch, pitch);
	
	
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Inicializamos variables de interés
	int ne, ni;

	double x0, y0, z0, t0, e0;

	double x1, y1, z1, t1, e1;

	int status;

	unsigned int nElastic, nIonising, nAttachment, nInelastic, nExcitation, nSuperelastic;

	int eventNumber; 

	
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Creamos el archivo Root y los histogramas que vamos a necesitar

		
	TFile f (rootfile, "RECREATE");


	int nBinsGain = 10;

	int gmin = 0;

	int gmax=10;

	TH1F * hElectrons =
		new TH1F ("hElectrons", "Number of electrons", 
				nBinsGain, gmin, gmax);


	int nBinsChrg = 400;
	

	TH1F * hEstartz =
		new TH1F ("hEstart_position_z", "hEstart_z", 200, 0,pitch);


	TH1F * hEstarty =
		new TH1F ("hEstart_position_y", "hEstart_y", 200,
				-3 * sqrt (3) * pitch / 2, 3 * sqrt (3) * pitch / 2);


	TH1F * hEstartx =
		new TH1F ("hEstart_position_x", "hEstart_x", 200, -3 * pitch / 2,
				3 * pitch / 2);


	TH1F * nend_elec_hist =
		new TH1F ("hChrgE_endz", "Charges_elec_endz", nBinsChrg * 100, -0.1,
				0.15);


	// Semilla para generar variables aleatorias en las posiciones x,y
	TRandom3 *rndistx = new TRandom3(0);
	TRandom3 *rndisty = new TRandom3(0);

	// Histograma de la energía de los electrones
	TH1F * histen = new TH1F ("hen", "energy distribution", 1000, 0.0, 100.0);

	
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Creamos las ramas que vamos a necesitar (tree)
	

	TTree *tree= new TTree("DataCharge","DataCharge");
	tree->Branch("eventNumber", &eventNumber,"eventNumber/I");  
	tree->Branch("PrimaryElectrons", &ne,"ne/I");  
	tree->Branch("ie", &ie,"ie/I");								// Numero de electrones generados
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
	dataExc->Branch("levelExc", &levelExc, "levelExc/D");

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Hacemos la simulación, para el número de electrones deseados

	// Usamos AvalancheMicroscopic para simular: 
	AvalancheMicroscopic * aval = new AvalancheMicroscopic ();
	aval->SetSensor (sensor);
	aval->EnableElectronEnergyHistogramming(histen);
	Int_t nelec_total = 0;

	// Set User Handle to follow the electron collicion-by-collision
	aval->SetUserHandleInelastic(userHandle);


	for (eventNumber = 0; eventNumber < npe; eventNumber++)
	{

		gas->ResetCollisionCounters();

		
		// Distancias en cm; aleatorio para no tener siempre los mismos datos
		x0 = -pitch+2*pitch*rndistx->Rndm();
		y0 = -pitch+2*pitch*rndisty->Rndm();

		// Emitimos el electrón en la parte de arriba
		z0 =  pitch*0.99;	

		t0 = 0.0;

		if (eventNumber == 0){
			e0=0.1;}

		else{
				e0 = histen->GetRandom(); // generate energy randomly accordingly to previous collisions
			}				

	

		std::cout << eventNumber + 1 << " electrones en un total de " << npe << endl;

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

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Guardamos las ramas e histogramas en el archivo ROOT

	

	tree->Write ();

	dataExc->Write();

	histen->Write ();

	hElectrons->Write ();

	nend_elec_hist->Write ();

	hEstartz->Write ();

	hEstarty->Write ();

	hEstartx->Write ();


	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// Borramos todos los elementos necesarios (importante) e imprimimos el tiempo de simulación


	printf ("average # of electrons produced: %g\n",
			((double) nelec_total) / ((double) npe));

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

	delete componentUniformElectricField;

	delete sensor;

	delete aval;

	delete histen;

	f.Close ();f.Close ();

	std::cout << "Finalizando correctamente uniformE()." << std::endl;
	return 0;
}
