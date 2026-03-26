#include <fstream>
#include <string>
#include <vector>
#include "Garfield/MediumMagboltz.hh"


using namespace Garfield;
using namespace std;

std::string CollisionTypeToString(const int type) {
  switch (type) {
    case 0: return "elastic";
    case 1: return "ionisation";
    case 2: return "attachment";
    case 3: return "inelastic";
    case 4: return "excitation";
    case 5: return "superelastic";
    default: return "unknown";
  }
}

void ExportAllLevelsCsv(Garfield::MediumMagboltz* gas,
                        const std::string& outCsv,
                        const std::vector<std::string>& gasNames) {
  std::ofstream fout(outCsv);
  fout << "level,gas,state_name,type,energy_eV,n_collisions\n";

  const int nLevels = gas->GetNumberOfLevels();
  for (int i = 0; i < nLevels; ++i) {
    int ngas = -1, type = -1;
    double e = 0.;
    std::string descr;

    if (!gas->GetLevel(i, ngas, type, descr, e)) continue;

    std::string gasName = "unknown";
    if (ngas >= 0 && ngas < (int)gasNames.size()) {
      gasName = gasNames[ngas];
    }

    const int nColl = gas->GetNumberOfElectronCollisions(i);

    fout << i << ","
         << "\"" << gasName << "\","
         << "\"" << descr << "\","
         << "\"" << CollisionTypeToString(type) << "\","
         << e << ","
         << nColl << "\n";
  }

  fout.close();
}

int main() {

  MediumMagboltz * gas = new MediumMagboltz ();
  gas->SetComposition("ar", 90., "cf4", 10.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->Initialise(true);

  const std::vector<std::string> gasNames = {"Ar", "CF4"};

  ExportAllLevelsCsv(gas, "level_gas_state.csv", gasNames);

  return 0;
}