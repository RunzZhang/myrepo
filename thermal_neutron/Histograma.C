#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TBox.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TProfile.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring> 
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include <stdint.h>
#include "TMath.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include <TRandom.h>


using namespace std;


void Histograma()
{


int i=0;
char name[50]="LED_U238_Up.txt";
char t[50];
float t1;
char eV [10]="eV";
char MeV [10]="MeV";
char keV [10]="keV";
char mm [10]="mm";
char Ar36 [10]="Ar36";
char Ar38 [10]="Ar38";
char Ar40 [10]="Ar40";
char part [12];
char comp [10];
char comp1 [10];
float edep, ek ;



ifstream file;
TH1* h1 = new TH1F("Datos", "Energia Cinetica de los Neutrones", 500, 0.0, 40000);             //Declaro el 1er histograma
TH1* h2 = new TH1F("Datos1", "Energia Depositada por los Neutrones", 500, 0.0, 6000);             //Declaro el 1er histograma

file.open(name);   //Abriendo el archivo txt con los datos
file>>t>>t>>t>>t;        //Paso por los nombres de las columnas
file>>ek>>comp>>edep>>comp1>>part>>t1>>t>>t1>>t>>t1>>t;

while (!file.eof())    //Se detiene cuando el fichero termina
{

 if(strcmp(part,Ar36) == 0 || strcmp(part,Ar38) == 0 || strcmp(part,Ar40) == 0)
  {    cout<<part<<endl;
       i++;

    if (strcmp(comp,MeV) == 0)    //Llevo todas las unidades a eV
    {
        ek=ek*1000000;
       // cout<<"jajaj"<<endl;
    }

    if (strcmp(comp,keV) == 0)    //Llevo todas las unidades a eV
    {
        ek=ek*1000;
       // cout<<"jajaj"<<endl;
    }

    if (strcmp(comp1,MeV) == 0)    //Llevo todas las unidades a eV
    {
        edep=edep*1000000;
       // cout<<"jajaj"<<endl;
    }

    if (strcmp(comp1,keV) == 0)    //Llevo todas las unidades a eV
    {
        edep=edep*1000;
       // cout<<"jajaj"<<endl;
    }

  //cout<<edep<<endl;
  h1->Fill(ek);
  h2->Fill(edep);
 }

 file>>ek>>comp>>edep>>comp1>>part>>t1>>t>>t1>>t>>t1>>t;


}

cout<<"Interactuaron "<< i << " neutrones"<<endl;

TCanvas *c1=new TCanvas("c1","Histograma",700,500);

h1->GetXaxis()->SetTitle("Energia_eV");
h1->GetYaxis()->SetTitle("Frecuencia");
h1->Draw("colz");
c1->Print("Histo_Energia_Cinetica_n.pdf");                //Para guardar el gráfico

TCanvas *c2=new TCanvas("c2","Histograma2",700,500);

h2->GetXaxis()->SetTitle("Energia_eV");
h2->GetYaxis()->SetTitle("Frecuencia");
h2->Draw("colz");
c2->Print("Histo_Energia_Depositada_n.pdf");                //Para guardar el gráfico
  
}
