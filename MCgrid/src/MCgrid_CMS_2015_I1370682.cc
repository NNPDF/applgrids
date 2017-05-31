// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "mcgrid/mcgrid.hh"

namespace Rivet {


  class MCgrid_CMS_2015_I1370682 : public Analysis {
  public:

    /// Constructor
    MCgrid_CMS_2015_I1370682()
      : Analysis("MCgrid_CMS_2015_I1370682")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Total xsec
      _h_xsec = bookHisto1D(50,1,1);

      _h_m_ttb  = bookHisto1D(23,1,1);
      _h_y_ttb  = bookHisto1D(22,1,1);
      _h_pt_ttb = bookHisto1D(21,1,1);
      _h_y_t  = bookHisto1D(17,1,1);
      _h_pt_t = bookHisto1D(15,1,1);

      IdentifiedFinalState topFinder; topFinder.acceptIdPair(6);
      addProjection(topFinder, "tops");

      // MCgrid
      const std::string configname = "MCgrid_CMS_2015_I1370682.config";
      MCgrid::bookPDF(configname, histoDir(), MCgrid::BEAM_PROTON, MCgrid::BEAM_PROTON);

      const double Q2min = 7500; // (0.25*(2*m_t))^2
      const double Q2max = 1.7E7; // 16007508.2225 GeV^2

      const double xmin = 1E-4;
      const double xmax = 1;

      MCgrid::gridArch arch(40,10,5,5);
      _a_xsec  = MCgrid::bookGrid(_h_xsec, histoDir(), configname, 2, xmin, xmax, Q2min, Q2max, arch);
      _a_m_ttb  = MCgrid::bookGrid(_h_m_ttb, histoDir(), configname, 2, xmin, xmax, Q2min, Q2max, arch);
      _a_y_ttb  = MCgrid::bookGrid(_h_y_ttb, histoDir(), configname, 2, xmin, xmax, Q2min, Q2max, arch);
      _a_pt_ttb = MCgrid::bookGrid(_h_pt_ttb, histoDir(), configname, 2, xmin, xmax, Q2min, Q2max, arch);
      _a_y_t  = MCgrid::bookGrid(_h_y_t, histoDir(), configname, 2, xmin, xmax, Q2min, Q2max, arch);
      _a_pt_t = MCgrid::bookGrid(_h_pt_t, histoDir(), configname, 2, xmin, xmax, Q2min, Q2max, arch);

      evt_Qmax = 0; evt_Qmin = 9000;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      MCgrid::PDFHandler::HandleEvent( event, histoDir());
      const double weight = event.weight();

      const FinalState& topFinder = applyProjection<IdentifiedFinalState>(event, "tops");
        if (topFinder.particles().size() != 2)
          vetoEvent;

      const FourMomentum top1 = topFinder.particles()[0].momentum();
      const FourMomentum top2 = topFinder.particles()[1].momentum();
      const FourMomentum ttbar = top1+top2;

      // Total inclusive
      _h_xsec->fill(0.5, weight);
      _a_xsec->fill(0.5, event);

      // ttbar system histograms
      _h_m_ttb->fill(ttbar.mass(),    weight);
      _h_y_ttb->fill(ttbar.rapidity(), weight);
      _h_pt_ttb->fill(ttbar.pt(),     weight);

      // ttbar system applgrids
      _a_m_ttb->fill(ttbar.mass(),    event);
      _a_y_ttb->fill(ttbar.rapidity(), event);
      _a_pt_ttb->fill(ttbar.pt(),     event);

      // single-top histograms
      _h_y_t->fill(top1.rapidity(),   weight);
      _h_pt_t->fill(top1.pt(),        weight);
      _h_y_t->fill(top2.rapidity(),   weight);
      _h_pt_t->fill(top2.pt(),        weight);

      // single-top applgrids
      _a_y_t->fill(top1.rapidity(),   event);
      _a_pt_t->fill(top1.pt(),        event);
      _a_y_t->fill(top2.rapidity(),   event);
      _a_pt_t->fill(top2.pt(),        event);

      evt_Qmin = min(evt_Qmin, top1.pt() + top2.pt());
      evt_Qmax = max(evt_Qmax, top1.pt() + top2.pt());
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      std::cout << "Q2min: " << pow(evt_Qmin,2)*0.25 <<  "Q2max: "<< pow(evt_Qmax,2)*0.25 <<endl;

      const double xsec_norm = crossSection()/sumOfWeights();
      const double xsec_plot = 1.0/sumOfWeights(); // Data is normalised by xsec

      // Scale histograms
      scale( _h_xsec, xsec_norm );

      scale ( _h_m_ttb,  xsec_plot ) ; 
      scale ( _h_y_ttb,  xsec_plot ) ; 
      scale ( _h_pt_ttb, xsec_plot ) ; 

      scale ( _h_y_t, xsec_plot/2.0 ) ; 
      scale ( _h_pt_t, xsec_plot/2.0 ) ; 

      // Scale APPLgrids
      _a_xsec->scale(xsec_norm);

      _a_m_ttb->scale(xsec_norm);
      _a_y_ttb->scale(xsec_norm);
      _a_pt_ttb->scale(xsec_norm);

      _a_y_t->scale(xsec_norm/2.0);
      _a_pt_t->scale(xsec_norm/2.0);

      // Export APPLgrids
      _a_xsec->exportgrid();
      _a_m_ttb->exportgrid();
      _a_y_ttb->exportgrid();
      _a_pt_ttb->exportgrid();
      _a_y_t->exportgrid();
      _a_pt_t->exportgrid();

      MCgrid::PDFHandler::CheckOutAnalysis(histoDir());
    }

    //@}


  private:

    double evt_Qmin;
    double evt_Qmax;

    /// @name Histograms
    //@{
    Histo1DPtr _h_m_ttb;
    Histo1DPtr _h_y_ttb;
    Histo1DPtr _h_pt_ttb;

    Histo1DPtr _h_y_t;
    Histo1DPtr _h_pt_t;

    Histo1DPtr _h_xsec;
    //@}

    /// @name APPLgrids
    //@{
    MCgrid::gridPtr _a_m_ttb;
    MCgrid::gridPtr _a_y_ttb;
    MCgrid::gridPtr _a_pt_ttb;

    MCgrid::gridPtr _a_y_t;
    MCgrid::gridPtr _a_pt_t;

    MCgrid::gridPtr _a_xsec;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MCgrid_CMS_2015_I1370682);


}
