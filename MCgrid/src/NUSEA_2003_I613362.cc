// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

#include "mcgrid/mcgrid.hh"
#include "mcgrid/mcgrid_binned.hh"
#include "HepMC/GenEvent.h"
#include <algorithm>

namespace Rivet {

  static const double MAXRAPIDITY = 10000000;
  static const double sqrts = 38.8;
  class NUSEA_2003_I613362 : public Analysis {
  public:

    /// Constructor
    NUSEA_2003_I613362()
      : Analysis("NUSEA_2003_I613362")
    {    }


    void init() {
      const ChargedFinalState cfs(-MAXRAPIDITY, MAXRAPIDITY, 0);
      const ChargedFinalState clfs(-MAXRAPIDITY, MAXRAPIDITY, 0);
      addProjection(cfs, "FS");
      addProjection(ChargedLeptons(clfs), "CL");

      // Initialise the subprocess PDF and grid architecture
      const string PDFname("NUSEA_2003_I613362.config");
      MCgrid::bookPDF(PDFname, histoDir(), MCgrid::BEAM_PROTON, MCgrid::BEAM_PROTON);
      MCgrid::gridArch arch(50,20,3,3);
      const double xmin = 1E-2; const double Q2min = 4.2*4.2;
      const double xmax = 1;    const double Q2max = 16.85*16.85;

      for (int i=0; i<bin_edge_xF.size()-2; i++)
      {
        Histo1DPtr      histogram = bookHisto1D(i+1, 1, 1);
        MCgrid::gridPtr applgrid  = MCgrid::bookGrid(histogram, histoDir(), PDFname, 0, xmin, xmax, Q2min, Q2max, arch);
        _hist_sigma.addHistogram(bin_edge_xF[i], bin_edge_xF[i+1], histogram);
        _appl_sigma.addGrid(bin_edge_xF[i], bin_edge_xF[i+1], applgrid);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {
      MCgrid::PDFHandler::HandleEvent(e, histoDir());

      const Particles& leptons = applyProjection<ChargedLeptons>(e, "CL").chargedLeptons();
      if (leptons.size() != 2 || leptons[0].pid() != -leptons[1].pid() ) vetoEvent;

      const FourMomentum dilepton = leptons[0].momentum() + leptons[1].momentum();
      const double M  = dilepton.mass();
      const double pZ = dilepton.pz();
      const double xF = 2.0*pZ / sqrts;

      if ( M > 8.7 && M < 10.85 ) vetoEvent;

      // Scale weights by M**3
      HepMC::GenEvent event_copy(*e.genEvent());
      for (int i=0; i<event_copy.weights().size(); i++)
        if (std::find(protect_wgts.begin(), protect_wgts.end(), i) == protect_wgts.end())
          event_copy.weights()[i] *= pow(M,3.0)/nanobarn;
      Event scaled_event(event_copy);

      _hist_sigma.fill(xF, M, scaled_event.weight());
      _appl_sigma.fill(xF, M, scaled_event);
    }


    void finalize() {
      _hist_sigma.scale(crossSection()/sumOfWeights(), this);
      _appl_sigma.scale(crossSection()/sumOfWeights());
      _appl_sigma.exportgrids();
      MCgrid::PDFHandler::CheckOutAnalysis(histoDir());
    }


  private:
    const std::array<int, 4> protect_wgts = {4,6,7,8};
    const std::array<double, 17>  bin_edge_xF = {-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8};
    BinnedHistogram<double> _hist_sigma;
    MCgrid::BinnedGrid<double> _appl_sigma;
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(NUSEA_2003_I613362);


}
