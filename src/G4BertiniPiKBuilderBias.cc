#include "G4BertiniPiKBuilderBias.hh"
#include <G4SystemOfUnits.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4PionPlus.hh>
#include <G4PionMinus.hh>
#include <G4KaonPlus.hh>
#include <G4ProcessManager.hh>
#include <G4BGGPionInelasticXS.hh>
#include <G4ComponentGGHadronNucleusXsc.hh>
#include <G4CrossSectionInelastic.hh>
#include <G4HadronicParameters.hh>


G4BertiniPiKBuilderBias::
G4BertiniPiKBuilderBias(double pion_plus_bias, double kaon_plus_bias)
  : fPionPlusBias(pion_plus_bias), fKaonPlusBias(kaon_plus_bias) {
   kaonxs = new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
   theMin = 0.0;
   theMax = G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade();
   theModel = new G4CascadeInterface;
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax); 
 }

void G4BertiniPiKBuilderBias::
Build(G4HadronInelasticProcess * aP)
 {
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
   // Geant4 11 has one generic Build; identify the particle by process name.
   const G4String& name = aP->GetProcessName();
   auto is = [&name](G4ParticleDefinition* p) {
       return name == p->GetParticleName() + "Inelastic";
   };
   if (is(G4PionPlus::Definition())) {
       aP->AddDataSet( new G4BGGPionInelasticXS( G4PionPlus::Definition() ) );
       std::cout << "Biasing PiPlus: " << fPionPlusBias << std::endl;
       aP->MultiplyCrossSectionBy(fPionPlusBias);
   }
   else if (is(G4PionMinus::Definition())) {
       aP->AddDataSet( new G4BGGPionInelasticXS( G4PionMinus::Definition() ) );
   }
   else {
       aP->AddDataSet(kaonxs);
       if (is(G4KaonPlus::Definition())) {
           std::cout << "Multiplying kaon by " << fKaonPlusBias << std::endl;
           aP->MultiplyCrossSectionBy(fKaonPlusBias);
       }
   }
   aP->RegisterMe(theModel);
 }
