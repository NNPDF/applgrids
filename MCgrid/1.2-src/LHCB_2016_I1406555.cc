// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  class LHCB_2016_I1406555 : public Analysis {
  public:

    /// Constructor
    LHCB_2016_I1406555()
      : Analysis("LHCB_2016_I1406555")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {

      _h_Wplus = bookHisto1D(1, 1, 1);
      _h_Wminus = bookHisto1D(1, 1, 3);
      _h_Z = bookHisto1D(2, 1, 1);

      FinalState fs;
      Cut cuts = Cuts::eta > 2.0 && Cuts::eta < 4.5 && Cuts::pT > 20*GeV;

      ZFinder zfinder_mu(fs, cuts, PID::MUON, 60*GeV, 120*GeV, 0);
      addProjection(zfinder_mu, "ZFinder_mu");

      WFinder wfinder_mu(fs,cuts, PID::MUON, 0*GeV, 1000*GeV, 25*GeV, 0);
      addProjection(wfinder_mu, "WFinder_mu");
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {
      const double weight = e.weight();

      const ZFinder& zfinder_mu = applyProjection<ZFinder>(e, "ZFinder_mu");
      const WFinder& wfinder_mu = applyProjection<WFinder>(e, "WFinder_mu");

      if (wfinder_mu.bosons().size()==1)
      {
        const Particle& W = wfinder_mu.boson();
        if (W.charge() > 0)
          _h_Wplus->fill(W.eta(), weight);
        else
          _h_Wminus->fill(W.eta(), weight);
      }

      if (zfinder_mu.bosons().size()==1)
        _h_Z->fill(zfinder_mu.boson().eta(), weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {


      scale(_h_Wplus, crossSection()/sumOfWeights()); // norm to cross section
      scale(_h_Wminus, crossSection()/sumOfWeights()); // norm to cross section
      scale(_h_Z, crossSection()/sumOfWeights()); // norm to cross section
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Wplus;
    Histo1DPtr _h_Wminus;
    Histo1DPtr _h_Z;

    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2016_I1406555);


}
