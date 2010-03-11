
// -*- C++ -*-
//
// Package:    DiJetAnalysis
// Class:      VertexAnalyzer
// 
/**\class VertexAnalyzer VertexAnalyzer.cc JSturdy/DiJetAnalysis/src/VertexAnalyzer.cc

Description: Collects variables related to vertices, performs a primary vertex check, 
             If successful, it stores the variables and returns the value of the check

*/
//
// Original Author:  Jared Sturdy
//         Created:  Fri Jan 29 16:10:31 PDT 2010
// $Id: VertexAnalyzer.cpp,v 1.1 2010/01/29 16:10:31 sturdy Exp $
//
//

#include "JSturdy/DiJetAnalysis/interface/VertexAnalyzer.h"
#include <TMath.h>
#include <sstream>

using namespace std;
using namespace reco;
using namespace edm;

//________________________________________________________________________________________
VertexAnalyzer::VertexAnalyzer(const edm::ParameterSet& pset)
{ 
  vertexParams = pset;
  //defaults
  _minNVtx    = 1;    //minimum number of vertices
  _minVtxTrks = 3;    //minimum number of tracks that contribute to the vertex
  _minVtxNdof = 4;    //minimum number of degrees of freedom for the vertex
  _maxVtxChi2  = 2.4;  //max chi2 of vertex is 2.4
  _maxVtxZ     = 15.0; //max z of vertex is 15 cm

  if (vertexParams.exists("minNVtx"))
    _minNVtx = vertexParams.getParameter<int>("minNVtx");
  if (vertexParams.exists("minVtxTrks"))
    _minVtxTrks = vertexParams.getParameter<int>("minVtxTrks");
  if (vertexParams.exists("minVtxNdof"))
    _minVtxNdof = vertexParams.getParameter<int>("minVtxNdof");
  if (vertexParams.exists("maxVtxChi2"))
    _maxVtxChi2    = vertexParams.getParameter<double>("maxVtxChi2");
  if (vertexParams.exists("maxVtxZ"))
    _maxVtxZ    = vertexParams.getParameter<double>("maxVtxZ");
    
  vtxTag_    = vertexParams.getParameter<edm::InputTag>("vtxTag"); 

  initTuple();
}


//________________________________________________________________________________________
VertexAnalyzer::~VertexAnalyzer() {}


//________________________________________________________________________________________
// Method called to for each event
bool
VertexAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::LogVerbatim("VertexAnalyzer") << " Start  " << std::endl;

  std::ostringstream dbg;

  // get the Vertex collection

  LogDebug("VertexAnalyzer") << "Vertex results for InputTag" << vtxTag_;
  Handle<VertexCollection> vertices;
  iEvent.getByLabel(vtxTag_, vertices);
  if ( !vertices.isValid() ) {
    LogDebug("VertexAnalyzer") << "No Vertex results for InputTag" << vtxTag_;
    return false;
  } 

  m_nVtx = (*vertices).size();

  const reco::Vertex* pVertex = &(*vertices)[0];
  m_VtxNormalizedChi2[0] = pVertex->normalizedChi2();
  m_VtxIsValid[0]        = pVertex->isValid();
  m_VtxChi2[0]           = pVertex->chi2();
  m_VtxNdof[0]           = pVertex->ndof();
  m_VtxX[0]              = pVertex->x();
  m_VtxY[0]              = pVertex->y();
  m_VtxZ[0]              = pVertex->z();
  m_VtxdX[0]             = pVertex->xError();
  m_VtxdY[0]             = pVertex->yError();
  m_VtxdZ[0]             = pVertex->zError();
  
  bool vertexDecision = false;
  if(m_nVtx>=_minNVtx)
    if(pVertex->isValid())
      //if(pVertex->Tracks()>=_minVtxTrks)
      if(pVertex->ndof()>=_minVtxNdof)
	//if(pVertex->chi2()<=_maxVtxChi2)
	if(pVertex->z()<=_maxVtxZ)
	  vertexDecision = true;
  
  
  if ((*vertices).size()>1) {
    for (int i=1; i< (int) (*vertices).size(); i++){  
      const reco::Vertex* pVertex = &(*vertices)[i];
      m_VtxNormalizedChi2[i] = pVertex->normalizedChi2();
      m_VtxIsValid[i]        = pVertex->isValid();
      m_VtxChi2[i]           = pVertex->chi2();
      m_VtxNdof[i]           = pVertex->ndof();
      m_VtxX[i]              = pVertex->x();
      m_VtxY[i]              = pVertex->y();
      m_VtxZ[i]              = pVertex->z();
      m_VtxdX[i]             = pVertex->xError();
      m_VtxdY[i]             = pVertex->yError();
      m_VtxdZ[i]             = pVertex->zError();
    } 
  }
  
  //mVertexData->Fill();
  return vertexDecision;
}

//________________________________________________________________________________________
void 
VertexAnalyzer::beginJob(const edm::EventSetup&) {}

//________________________________________________________________________________________
void 
VertexAnalyzer::endJob() {

}

//________________________________________________________________________________________
void
VertexAnalyzer::initTuple() {

  std::ostringstream variables; // Container for all variables
  
  // 1. Event variables
  variables << "weight:process";
  
  // Register this ntuple
  edm::Service<TFileService> fs;

  mVertexData = fs->make<TTree>( "VertexData", "data after cuts" );
  mVertexData->SetAutoSave(10);

  //other information

  mVertexData->Branch("nVtx",                &m_nVtx,            "nVtx/int");
  mVertexData->Branch("VertexChi2",          m_VtxChi2,          "VertexChi2[nVtx]/double");
  mVertexData->Branch("VertexNdof",          m_VtxNdof,          "VertexNdof[nVtx]/double");
  mVertexData->Branch("VertexIsValid",       m_VtxIsValid,       "VertexIsValid[nVtx]/double");
  mVertexData->Branch("VertexNormalizedChi2",m_VtxNormalizedChi2,"VertexNormalizedChi2[nVtx]/double");

  mVertexData->Branch("VertexX", m_VtxX, "VertexX[nVtx]/double");
  mVertexData->Branch("VertexY", m_VtxY, "VertexY[nVtx]/double");
  mVertexData->Branch("VertexZ", m_VtxZ, "VertexZ[nVtx]/double");
  mVertexData->Branch("VertexdX",m_VtxdX,"VertexdX[nVtx]/double");
  mVertexData->Branch("VertexdY",m_VtxdY,"VertexdY[nVtx]/double");
  mVertexData->Branch("VertexdZ",m_VtxdZ,"VertexdZ[nVtx]/double");
    
  edm::LogInfo("VertexAnalyzer") << "Ntuple variables " << variables.str();
  
}

//_______________________________________________________________________________________
// Define this as a plug-in
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(VertexAnalyzer);
