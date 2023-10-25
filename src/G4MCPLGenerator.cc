#include "G4MCPLGenerator.hh"
#include "G4ParticleGun.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ios.hh"
#include <cassert>
#include <G4ParticleDefinition.hh>
#include "G4Threading.hh"
#include <boost/lexical_cast.hpp>
#include "RunAction.hh"
#include "G4ParticleTable.hh"

#include <arrow/io/file.h>
#include <parquet/stream_writer.h>
#include "parquet_writer.h"

using namespace std;

extern G4double event_number_global;
extern G4int run_number;

mcpl_file_t m_mcplfile;
const mcpl_particle_t* m_p;

extern parquet::StreamWriter Particle_os;
// extern std::ofstream Particle_outfile;

using boost::lexical_cast;
#include "G4AutoLock.hh"

G4ThreadLocal G4double event_number;
G4ThreadLocal G4int local_event_number_MCPL;
G4ThreadLocal int Event_found = 0; 

namespace { G4Mutex MCPLMutex = G4MUTEX_INITIALIZER;}


G4MCPLGenerator::G4MCPLGenerator(const G4String& inputFile)
  : G4VUserPrimaryGeneratorAction(),
    m_currentPDG(0),
    m_currentPartDef(0),
    m_nUnknownPDG(0),
    m_inputFile(inputFile)
{
  G4AutoLock lock(&MCPLMutex);
  m_mcplfile.internal = 0;
  m_mcplfile = mcpl_open_file(m_inputFile.c_str());
  m_gun = new G4ParticleGun(1);
  // m_p is a null pointer, run FindNext to get to the first entry in the mcpl list
  FindNext();
}

G4MCPLGenerator::~G4MCPLGenerator()
{
  if (m_nUnknownPDG) {
    std::ostringstream cmt;
    cmt << "Ignored a total of " << m_nUnknownPDG << " particles in input due to untranslatable pdg codes";
    G4Exception("G4MCPLGenerator::~G4MCPLGenerator()", "G4MCPLGenerator07",JustWarning, cmt.str().c_str());
  }
  if (m_mcplfile.internal && G4Threading::G4GetThreadId() == 0){mcpl_close_file(m_mcplfile);} // only one file, can only be closed by one thread! 
  delete m_gun;
}

bool G4MCPLGenerator::UseParticle(const mcpl_particle_t*) const{return true;}

void G4MCPLGenerator::ModifyParticle(G4ThreeVector&, G4ThreeVector&, G4ThreeVector&, G4double&, G4double&) const{}

void G4MCPLGenerator::GeneratePrimaries(G4Event* evt)
{
  
  auto particleTable = G4ParticleTable::GetParticleTable();

  //Open the mcpl file specified
  G4AutoLock lock(&MCPLMutex);
  local_event_number_MCPL = event_number_global; 
  evt->SetEventID(event_number_global);

  std::cout << "Thread ID " << G4Threading::G4GetThreadId() <<  " ######################### Event ID: " << evt->GetEventID() << " mcpl flag " << m_p->userflags << std::endl;
  
  event_number_global++;
  
  Event_found = 0; 

  // If the particle it has is exactly what we are looking for: shoot it, and check the next
  if (m_p->userflags == local_event_number_MCPL) {
    
      Event_found = 1;
        
        auto PartDef = particleTable->FindParticle(m_p->pdgcode);
        if ( !PartDef && ((m_p->pdgcode)/100000000 == 10)) {
          //Not in ParticleTable and pdgcode is of form 10xxxxxxxx, so look for ion:
          PartDef = G4IonTable::GetIonTable()->GetIon(m_p->pdgcode);
        }

        m_gun->SetParticleDefinition(PartDef);
        G4ThreeVector pos(m_p->position[0],m_p->position[1],m_p->position[2]); // the output from mcpl file is mm
        pos*=CLHEP::cm;
        G4ThreeVector dir(m_p->direction[0],m_p->direction[1],m_p->direction[2]);
        G4ThreeVector pol(m_p->polarisation[0],m_p->polarisation[1],m_p->polarisation[2]);
        G4double KE = m_p->ekin;
        G4double time = 0.;//m_p->time*CLHEP::millisecond;
        G4double weight = m_p->weight;
        ModifyParticle(pos,dir,pol,time,weight);

        m_gun->SetParticleMomentumDirection(dir);
        m_gun->SetParticlePosition(pos);
        m_gun->SetParticleEnergy(KE);//already in MeV and CLHEP::MeV=1
        m_gun->SetParticleTime(time); //time
        m_gun->SetParticlePolarization(pol);
      
          Particle_os << evt->GetEventID()<< m_p->pdgcode << particleTable -> FindParticle(m_p->pdgcode) -> GetPDGMass() 
          << particleTable -> FindParticle(m_p->pdgcode)->GetParticleName()<< particleTable -> FindParticle(m_p->pdgcode) -> GetPDGCharge() 
          << KE << 0.0  
          << pos[0]/CLHEP::cm << pos[1]/CLHEP::cm << pos[2]/CLHEP::cm << time/CLHEP::ms 
          << dir[0] << dir[1] << dir[2] << weight 
          << parquet::EndRow;                    
     
          m_gun->GeneratePrimaryVertex(evt);
        
        std::cout << "Thread ID " << G4Threading::G4GetThreadId() << 
                  " :: Event ID from software" << evt->GetEventID()<<
                  " Particle:" << m_p->pdgcode << ", KE: " << m_p->ekin << "pos : " <<  pos/CLHEP::cm << " time: " << time << std::endl;
    
    // then start searching the next entry in the mcpl list
      while ((m_p = mcpl_read(m_mcplfile))){
        //std::cout<< "Thread "<< G4Threading::G4GetThreadId() << " Looking for event " << local_event_number_MCPL << " :: MCPL reading -- 1" << std::endl;
        // if the next entry event number is not the same as what we are searching -> the end of our event
        if (m_p->userflags != local_event_number_MCPL && Event_found==1){break;}
        if (m_p->userflags == local_event_number_MCPL){
              auto PartDef = particleTable->FindParticle(m_p->pdgcode);
              if ( !PartDef && ((m_p->pdgcode)/100000000 == 10)) {
                //Not in ParticleTable and pdgcode is of form 10xxxxxxxx, so look for ion:
                PartDef = G4IonTable::GetIonTable()->GetIon(m_p->pdgcode);
              }
              m_gun->SetParticleDefinition(PartDef);
              G4ThreeVector pos(m_p->position[0],m_p->position[1],m_p->position[2]);
              pos*=CLHEP::cm;
              G4ThreeVector dir(m_p->direction[0],m_p->direction[1],m_p->direction[2]);
              G4ThreeVector pol(m_p->polarisation[0],m_p->polarisation[1],m_p->polarisation[2]);
              G4double KE = m_p->ekin;
              G4double time = 0.; //m_p->time*CLHEP::millisecond;
              G4double weight = m_p->weight;
              ModifyParticle(pos,dir,pol,time,weight);

              m_gun->SetParticleMomentumDirection(dir);
              m_gun->SetParticlePosition(pos);
              m_gun->SetParticleEnergy(KE);//already in MeV and CLHEP::MeV=1
              m_gun->SetParticleTime(time); //time
              m_gun->SetParticlePolarization(pol);
              
                  Particle_os << evt->GetEventID()<< m_p->pdgcode << particleTable -> FindParticle(m_p->pdgcode) -> GetPDGMass() 
                  << particleTable -> FindParticle(m_p->pdgcode)->GetParticleName() << particleTable -> FindParticle(m_p->pdgcode) -> GetPDGCharge() 
                  << KE << 0.0
                  << pos[0]/CLHEP::cm << pos[1]/CLHEP::cm << pos[2]/CLHEP::cm << time/ms 
                  << dir[0] << dir[1] << dir[2] << weight 
                  << parquet::EndRow;                    
                   

                  m_gun->GeneratePrimaryVertex(evt);

                  std::cout << "Thread ID " << G4Threading::G4GetThreadId() << 
                  " :: Event ID from software" << evt->GetEventID()<<
                  " Particle:" << m_p->pdgcode << ", KE: " << m_p->ekin << "pos : " <<  pos/CLHEP::cm << " time: " << time << std::endl;

        }
      }
  }
  
  // If the particle is not what we want -> We keep looking for the next entry
  else{
    while ((m_p = mcpl_read(m_mcplfile))){
      //std::cout << "Thread " << G4Threading::G4GetThreadId() << " Looking for event " << local_event_number_MCPL << " :: MCPL reading -- 2 ::" << Event_found<< std::endl;
      if (m_p->userflags != local_event_number_MCPL && Event_found==1){break;}
      if (m_p->userflags == local_event_number_MCPL){
        Event_found = 1;
          auto PartDef = particleTable->FindParticle(m_p->pdgcode);
          if ( !PartDef && ((m_p->pdgcode)/100000000 == 10)) {
            //Not in ParticleTable and pdgcode is of form 10xxxxxxxx, so look for ion:
            PartDef = G4IonTable::GetIonTable()->GetIon(m_p->pdgcode);
          }
          m_gun->SetParticleDefinition(PartDef);
          G4ThreeVector pos(m_p->position[0],m_p->position[1],m_p->position[2]);
          pos*=CLHEP::cm;
          G4ThreeVector dir(m_p->direction[0],m_p->direction[1],m_p->direction[2]);
          G4ThreeVector pol(m_p->polarisation[0],m_p->polarisation[1],m_p->polarisation[2]);
          G4double KE = m_p->ekin;
          G4double time = 0.; //m_p->time*CLHEP::millisecond;
          G4double weight = m_p->weight;
          ModifyParticle(pos,dir,pol,time,weight);

          m_gun->SetParticleMomentumDirection(dir);
          m_gun->SetParticlePosition(pos);
          m_gun->SetParticleEnergy(KE);//already in MeV and CLHEP::MeV=1
          m_gun->SetParticleTime(time); //time
          m_gun->SetParticlePolarization(pol);
          
            Particle_os << evt->GetEventID()<< m_p->pdgcode << particleTable -> FindParticle(m_p->pdgcode) -> GetPDGMass() 
            << particleTable -> FindParticle(m_p->pdgcode)->GetParticleName() << particleTable -> FindParticle(m_p->pdgcode) -> GetPDGCharge() 
            << KE << 0.0
            << pos[0]/CLHEP::cm << pos[1]/CLHEP::cm << pos[2]/CLHEP::cm << time/ms 
            << dir[0] << dir[1] << dir[2] << weight 
            << parquet::EndRow;                    
            
          m_gun->GeneratePrimaryVertex(evt);
              
          std::cout << "Thread ID " << G4Threading::G4GetThreadId() << 
                       " :: Event ID from software" << evt->GetEventID()<<
                       " Particle:" << m_p->pdgcode << ", KE: " << m_p->ekin << "pos : " <<  pos/CLHEP::cm << " time: " << time << std::endl;
      }
    }
  }

  if (Event_found==0){
    std::cout<< "MCPL list running out -- run must be aborted" << std::endl;
    G4RunManager::GetRunManager()->AbortRun(true);//hard abort
    return;
  }
}

void G4MCPLGenerator::FindNext()
{

  while( ( m_p = mcpl_read(m_mcplfile))) {

    if (!UseParticle(m_p)) {continue;}
    if (!(m_p->weight>0.0)) {continue;}
    if (m_p->pdgcode==0) {continue;}
    if (!LookupPDG(m_p->pdgcode)) {
      ++m_nUnknownPDG;
      if (m_nUnknownPDG<=100) {
        std::ostringstream cmt;
        cmt << "Ignoring particle in input with untranslatable pdg code ("<< m_p->pdgcode <<")";
        G4Exception("G4MCPLGenerator::GeneratePrimaries()", "G4MCPLGenerator05",JustWarning, cmt.str().c_str());
        if (m_nUnknownPDG==100)
          G4Exception("G4MCPLGenerator::GeneratePrimaries()", "G4MCPLGenerator06",
                      JustWarning, "Limit reached. Suppressing further warnings"
                      " regarding untranslatable pdg codes");
      }
      continue;
    }
    break;
  }
}

G4ParticleDefinition* G4MCPLGenerator::LookupPDG(G4int pdgcode)
{
  if (m_currentPDG == pdgcode)
    return m_currentPartDef;
  m_currentPDG = pdgcode;
  std::map<G4int,G4ParticleDefinition*>::const_iterator it = m_pdg2pdef.find(pdgcode);
  if (it!=m_pdg2pdef.end()) {
    m_currentPartDef = it->second;
  } else {
    m_currentPartDef = G4ParticleTable::GetParticleTable()->FindParticle(pdgcode);
    if ( !m_currentPartDef && (pdgcode/100000000 == 10)) {
      //Not in ParticleTable and pdgcode is of form 10xxxxxxxx, so look for ion:
      m_currentPartDef = G4IonTable::GetIonTable()->GetIon(pdgcode);
    }
    m_pdg2pdef[pdgcode] = m_currentPartDef;
  }
  return m_currentPartDef;
}
