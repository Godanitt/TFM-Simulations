#include <iostream>
#include <cmath> 
#include <cstring>
#include <fstream>
#include <TCanvas.h>
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>

#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ComponentElmer.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/AvalancheMicroscopic.hh"

using namespace Garfield;

// the program must be called like: ./ViewGeometry

int main(int argc, char * argv[]) {

    TApplication app("app", &argc, argv);

    // Define the medium.
    MediumMagboltz* gas     = new MediumMagboltz();
    double pressure_bar     = 0.1;
    double cf_bartoTorr     = 750.062;
    double pressure_Torr    = pressure_bar*cf_bartoTorr;
    double temperature_Cdeg = 20;
    double cf_CtoK          = 273.15;
    double temperature_Kdeg = temperature_Cdeg+cf_CtoK;
    // gas->EnableDrift();                           // Allow for drifting in this medium
    gas->SetTemperature(temperature_Cdeg);
    gas->SetPressure(pressure_Torr);
    gas->SetComposition("Xe", 100.);              // Specify the gas mixture (pure Xe)

    // gas->Initialise();

    // Load the field map from ansys
    char elist_file[500]="ELIST.lis";
    char nlist_file[500]="NLIST.lis";
    char mplist_file[500]= "MPLIST.lis";
    char prnsol_file[500]="PRNSOL.lis";
    ComponentAnsys123 * fm = new ComponentAnsys123 ();
    
    fm->Initialise (elist_file, nlist_file, mplist_file, prnsol_file, "mm");

    const int nMaterials = fm->GetNumberOfMaterials();
    for (int i = 0; i < nMaterials; ++i) {
        const double eps = fm->GetPermittivity(i);
        if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);}
    // Print a list of the field map materials (for information).
    fm->PrintMaterials();
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
	fm->PrintRange ();
    // fm->SetGas(gas);

    // Set up a sensor object
    const double axis_x = 0.5;  // X-width will cover between +/- axis_x
    const double axis_y = 0.5;  // Y-width will cover between +/- axis_y
    const double axis_z = 1.5;
    Sensor* sensor = new Sensor();
    sensor->AddComponent(fm);
    sensor->SetArea(-axis_x,-axis_y,-0.8,axis_x,axis_y,axis_z);

    AvalancheMicroscopic *aval=new AvalancheMicroscopic();
    aval->SetSensor(sensor);
    // aval->SetCollisionSteps(10000);

    // Set the initial position [cm] and starting time [ns].
    double X0 = 0., Y0 = 0., Z0 = 1.0, T0 = 0.; //nota x0 é um array definido acima. Assim defino como X=
// Set the initial energy [eV].
    double E0 = 0.1;
// Set the initial direction (x, y, z).
// In case of a null vector, the direction is randomized.
    double dx0 = 0., dy0 = 0., dz0 = 0.;
// Calculate an electron avalanche.

    ViewDrift *viewDrift=new ViewDrift();
    viewDrift->SetArea(-axis_x, -axis_y, -0.8, axis_x, axis_y, axis_z);
    viewDrift->SetColourExcitations(kGreen);
    aval->EnablePlotting(viewDrift);
   
   //NOTA: Só posso lançar o electrão depois da linha acima, senão não faz o plot, pois não guardou a informação
    aval->AvalancheElectron(X0, Y0, Z0, T0, E0, dx0, dy0, dz0);


    TCanvas * cGeom = new TCanvas("geom","Geometry/Avalanche/Fields");
    //// Set up the object for FE mesh visualization.
    ViewFEMesh * vFE = new ViewFEMesh();
    vFE->SetCanvas(cGeom);
    vFE->SetComponent(fm);
    vFE->SetPlane(0,-1,0,0,0,0);//->Front View
    //vFE->SetPlane(0,1,0,0,0,0);//->Front View
    // vFE->SetPlane(0,0,-1,0,0,0);//->top View
    vFE->SetFillMesh(true);
    vFE->SetColor(0, kGray);//Metal
    vFE->SetColor(1, 3);//FR4
    vFE->SetColor(2, kRed+1);
    //void SetFillColor (int matID, int colorID);

    vFE->SetArea(-1*axis_x,-1*axis_y,-0.8,axis_x,axis_y,axis_z);      //       z-axis is irrelevant for 2D projections
    vFE->EnableAxes();             // comment this to disable creation of independent axes when contours are plotted
    vFE->SetXaxisTitle("y (cm)");
    vFE->SetYaxisTitle("z (cm)");
    vFE->SetViewDrift(viewDrift);
    vFE->Plot();

    //TESTES viewDrift

  



//     //Calculate one avalanche


// TCanvas *canvasDrift =new TCanvas("canvasDrift","canvasDrift");
viewDrift->Plot();


//DEBUG REGION
int ne, ni;
// Get the number of electrons and ions in the avalanche.
aval->GetAvalancheSize(ne, ni);
std::cout<<ne<<std::endl;

int np = aval->GetNumberOfElectronEndpoints();
// Initial position and time
double x1, y1, z1, t1;
// Final position and time
double x2, y2, z2, t2;
// Initial and final energy
double e1, e2;
// Flag indicating why the tracking of an electron was stopped.
int status;
// Loop over the avalanche electrons.
for (int i = 0; i < np; ++i) {
  aval->GetElectronEndpoint(i, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, status);
} 

std::cout<<x1<<" "<<y1<<" "<<z1<<" "<<std::endl;
std::cout<<x2<<" "<<y2<<" "<<z2<<" "<<std::endl;
//END DEBUG



	std::cout<<"END"<<std::endl;
	
    app.Run(kTRUE);

    return 0;
}
