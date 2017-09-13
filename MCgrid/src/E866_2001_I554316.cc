// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "mcgrid/mcgrid.hh"

namespace Rivet {

  static const double MAXRAPIDITY = 10000000;
  static const double sqrts = 38.8;
  class E866_2001_I554316 : public Analysis {
  public:

    /// Constructor
    E866_2001_I554316()
      : Analysis("E866_2001_I554316")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Set up projections
      const ChargedFinalState cfs(-MAXRAPIDITY, MAXRAPIDITY, 0);
      const ChargedFinalState clfs(-MAXRAPIDITY, MAXRAPIDITY, 0);
      addProjection(cfs, "FS");
      addProjection(ChargedLeptons(clfs), "CL");

      // Booking histograms
      _h_average = bookHisto1D(4, 1, 1);

      const MCgrid::gridArch arch(50,30,3,3);
      const std::string configname = "E866_2001_I554316.config";
      MCgrid::bookPDF(configname, histoDir(), MCgrid::BEAM_PROTON, MCgrid::BEAM_PROTON);
      _a_average = MCgrid::bookGrid(_h_average, histoDir(), configname, 0, 1E-2, 1, 4*4, sqrts*sqrts, arch);

    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {
      MCgrid::PDFHandler::HandleEvent( e, histoDir());
      const double weight = e.weight();
      const Particles& leptons = applyProjection<ChargedLeptons>(e, "CL").chargedLeptons();
      if (leptons.size() != 2 || leptons[0].pid() != -leptons[1].pid() ) vetoEvent;

      const FourMomentum dilepton = leptons[0].momentum() + leptons[1].momentum();
      const double M  = dilepton.mass();
      const double M2 = dilepton.mass2();
      const double pZ = dilepton.pz();

      if (M/GeV < 4.0 || ( M/GeV > 8.7 && M/GeV < 10.85) ) // Cut out J/psi and upsilon resonances
        vetoEvent;

      const double x_2 = -(1.0/sqrts)*(pZ - sqrt(pZ*pZ + M2));
      _h_average->fill(x_2, weight);
      _a_average->fill(x_2, e);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_average, crossSection()/sumOfWeights()); // norm to cross section
      _a_average->scale(crossSection()/sumOfWeights());
      _a_average->exportgrid();

      MCgrid::PDFHandler::CheckOutAnalysis(histoDir());
    }

    //@}


  private:
    // Average over mass-bins
    Histo1DPtr      _h_average; 
    MCgrid::gridPtr _a_average;
  };



  DECLARE_RIVET_PLUGIN(E866_2001_I554316);


}
