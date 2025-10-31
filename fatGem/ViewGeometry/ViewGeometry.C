
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
    double pressure_bar     = 1.0;
    double cf_bartoTorr     = 750.062;
    double pressure_Torr    = pressure_bar*cf_bartoTorr;
    double temperature_Cdeg = 20;
    double cf_CtoK          = 273.15;
    double temperature_Kdeg = temperature_Cdeg+cf_CtoK;
    gas->EnableDrift();                           // Allow for drifting in this medium
     gas->SetTemperature(temperature_Cdeg);
    gas->SetPressure(pressure_Torr);
    gas->SetComposition("Xe", 100.);              // Specify the gas mixture (pure Xe)


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


    // Mega importante: permite repetir el campo eléctrico y geometría a lo largo de los ejes y y z
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
	fm->PrintRange ();


    // Create several canvases for the plots.
   
    TCanvas * cField = new TCanvas("Field","Field");
    


    // Set up a sensor object
    const double axis_x = 1.0;  // X-width will cover between +/- axis_x
    const double axis_y = 1.0;  // Y-width will cover between +/- axis_y
    const double axis_z = 1.7;
    Sensor* sensor = new Sensor();
    sensor->AddComponent(fm);
    sensor->SetArea(-2*axis_x,-2*axis_y,-1*axis_z,2*axis_x,2*axis_y,axis_z);

    // Set up the object for field visualization.
    ViewField * vf = new ViewField();
    vf->SetComponent(fm);
    vf->SetCanvas(cField);
    vf->SetArea(-1*axis_x,-1.0,axis_x,2.);
    vf->SetNumberOfContours(100);
    // vf->SetNumberOfSamples2d(200,200);
    vf->SetPlane(0,-1,0,0,0,0);
    //vf->SetVoltageRange(-15000.,0.);
	//vf-> SetElectricFieldRange(0.,2000.0);
    vf->PlotContour("v");//"e"->electric field
    //vf->PlotContour("e");//"e"->electric field

 
    TCanvas * cFieldLines = new TCanvas("fieldLines","fieldLines");

    ViewField *vfLines = new ViewField();
    vfLines->SetComponent(fm);
    vfLines->SetSensor(sensor);
    vfLines->SetCanvas(cFieldLines);
    vfLines->SetArea(-1*axis_x,-1,1*axis_x,axis_z);
    //vfLines->SetNumberOfContours(100);
    vfLines->SetPlane(0,-1,0,0,0,0);

    const int nLines=100;
    std::vector<double> x0 (nLines);
    std::vector<double> y0 (nLines);
    std::vector<double> z0 (nLines);

    std::fill(z0.begin(),z0.end(),1.5);


//NOTA: HA QQ COISA QUE HA MEDIDA QUE AUMENTO O NUMERO D ELINHAS TENHO DE DIMINUIR O FACTOR DE MULTRIPLICACAO ABAIXO
//TALVEZ TENHA DE TER UM SENSOR MAIOR

    for (int i = 0; i < x0.size(); ++i){
        x0[i] = -0.8*axis_x+0.8*2*i*axis_x/nLines;
        y0[i] = -0.8*axis_y+0.8*2*i*axis_y/nLines;
    }
    vfLines->PlotFieldLines(x0,y0,z0,true,true);


// gPad->ls();


    TCanvas * cGeom = new TCanvas("geom","Geometry/Avalanche/Fields");
    //// Set up the object for FE mesh visualization.
    ViewFEMesh * vFE = new ViewFEMesh();
    vFE->SetCanvas(cGeom);
    vFE->SetComponent(fm);
    // vFE->SetPlane(0,-1,0,0,0,0);//->Front View
    //vFE->SetPlane(0,1,0,0,0,0);//->Front View
    vFE->SetPlane(0,0,-1,0,0,0);//->top View
    vFE->SetFillMesh(true);
    vFE->SetColor(0, kGray);//Metal
    vFE->SetColor(1, 3);//FR4
    vFE->SetColor(2, kRed+1);
    //void SetFillColor (int matID, int colorID);

    vFE->SetArea(-1*axis_x,-1*axis_y,-1.,axis_x,axis_y,2.);      //       z-axis is irrelevant for 2D projections
    vFE->EnableAxes();             // comment this to disable creation of independent axes when contours are plotted
    vFE->SetXaxisTitle("y (cm)");
    vFE->SetYaxisTitle("z (cm)");
    vFE->Plot();


	

    TCanvas *cFieldProjection = new TCanvas("fieldProjection","fieldProjection");
    vf->SetCanvas(cFieldProjection);
    vf->PlotProfile(0,0,-axis_z,0,0,axis_z,"e");//plot eletric field profile along z, at position x=0 and y=0

    


	std::cout<<"END"<<std::endl;
	
    app.Run(kTRUE);

    return 0;
}
