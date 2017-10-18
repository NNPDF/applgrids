// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"

#include "mcgrid/mcgrid.hh"

namespace Rivet {


  class D0_2014_I1333394 : public Analysis {
  public:

    /// Constructor
    D0_2014_I1333394()
      : Analysis("D0_2014_I1333394")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {

      const double minMass = 50.0*GeV;
      const double maxMass = 130.0*GeV;
      const double missingET = 25.0*GeV;

      FinalState fs;
      WFinder wfinder(fs, Cuts::abseta < 3.2 && Cuts::pT > 25*GeV, PID::ELECTRON, minMass, maxMass, missingET, 0, WFinder::NOCLUSTER, WFinder::NOTRACK, WFinder::TRANSMASS);
      addProjection(wfinder, "WFinder");
      
      // A bit of a hack - but the binning is the same
      _tmp_h_plus = bookHisto1D(1, 1, 2);
      _tmp_h_minus = bookHisto1D(2, 1, 2);

      _h_asym = bookHisto1D(1, 1, 1);

      // MCgrid
      const std::string configname = "D0_2014_I1333394.config";
      MCgrid::bookPDF(configname, histoDir(), MCgrid::BEAM_PROTON, MCgrid::BEAM_ANTIPROTON);

      MCgrid::gridArch arch(50,1,5,0);
      const double MW2 = 6463.838404;
      _a_plus = MCgrid::bookGrid(_tmp_h_plus, histoDir(), configname, 0, 1E-5, 1, MW2, MW2, arch);
      _a_minus = MCgrid::bookGrid(_tmp_h_minus, histoDir(), configname, 0, 1E-5, 1, MW2, MW2, arch);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      MCgrid::PDFHandler::HandleEvent( event, histoDir());

      const WFinder& wfinder = applyProjection<WFinder>(event, "WFinder");

      if (wfinder.bosons().size() != 1) 
        vetoEvent;
 
      // Get lepton constituent
      Particle l=wfinder.constituentLeptons()[0];

      const bool isMinus = l.pid() > 0;
      const bool isPosEta = l.eta() > 0;

      // Includes -A(-eta) = A(eta);
      Histo1DPtr      htmp = ( isMinus == isPosEta ) ?  _tmp_h_minus:_tmp_h_plus;
      MCgrid::gridPtr atmp = ( isMinus == isPosEta ) ?  _a_minus:_a_plus;
      htmp->fill(fabs(l.eta()), event.weight());
      atmp->fill(fabs(l.eta()), event);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double xsec_norm = crossSection()/sumOfWeights();
      scale ( _tmp_h_plus ,  xsec_norm ) ; 
      scale ( _tmp_h_minus , xsec_norm ) ; 
      _a_plus->scale  ( xsec_norm ) ;
      _a_minus->scale ( xsec_norm ) ;

      // Use height rather than sumW here to get properly normalised distributions
      for (size_t i = 0; i < _tmp_h_plus->numBins(); ++i) {
        const double num   = _tmp_h_plus->bin(i).height() - _tmp_h_minus->bin(i).height();
        const double denom = _tmp_h_plus->bin(i).height() + _tmp_h_minus->bin(i).height();
        const double asym = (num != 0 && denom != 0) ? num / denom : 0;
        const double norm = 100.0*_tmp_h_plus->bin(i).width(); // Data scaled by 100
        _h_asym->fill(_tmp_h_plus->bin(i).midpoint(), norm*asym);
      }

      // Normalise and export APPLgrids
      _a_plus->exportgrid();
      _a_minus->exportgrid();

      MCgrid::PDFHandler::CheckOutAnalysis(histoDir());
    }



  private:

    Histo1DPtr _h_asym;
    Histo1DPtr  _tmp_h_plus, _tmp_h_minus;

    MCgrid::gridPtr _a_plus;
    MCgrid::gridPtr _a_minus;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2014_I1333394);


}
