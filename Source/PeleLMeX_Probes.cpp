#include <PeleLMeX.H>
#include <PeleLMeX_BPatch.H>

void
PeleLM::probesSetup()
{
  std::string pele_prefix = "peleLM";
  amrex::ParmParse pp(pele_prefix);

  int n_probes = 0;
  n_probes = pp.countval("probes");
  if (n_probes > 0) {
    m_have_probes = 1;
    m_probes_coord.resize(n_probes);
    m_probes_name.resize(n_probes);
    m_probes_cell.resize(n_probes);
    m_probes_weight.resize(n_probes);

    std::string probe_prefix = pele_prefix + ".probes";
    amrex::ParmParse ppp(probe_prefix);
    for (int n = 0; n < n_probes; ++n) {
      pp.get("probes", m_probes_name[n], n);
      amrex::Vector<amrex::Real> coord(AMREX_SPACEDIM);
      ppp.getarr(m_probes_name[n], coord, 0, AMREX_SPACEDIM);
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_probes_coord[n][idim] = coord[idim];
      }
    }

    ppp.query("int", m_probes_int);

    int n_probes_var = 0;
    n_probes_var = ppp.countval("variables");
    if (n_probes_var > 0) {
      m_probes_vars.resize(n_probes_var);
      for (int ivar = 0; ivar < n_probes_var; ++ivar) {
        ppp.get("variables", m_probes_vars[ivar], ivar);
        // Need to check the availability of the variables, but the derived
        // arent initialized yet
      }
    } else {
        amrex::Abort("Using probes explicitly requires to probe of list in peleLM.probes.variables");
    }
  }
}

void
PeleLM::updateProbesCellsAndWeights()
{
  if (m_have_probes == 0) {
    return;
  }
}

void
PeleLM::writeProbes()
{
  tmpProbesFile << m_nstep << "," << m_cur_time << "," << m_dt // Time
               << "," << "Test"
               << "\n";
  tmpProbesFile.flush();
}

void
PeleLM::openProbesFile()
{
  if (m_have_probes == 0) {
    return;
  }

  // Create the temporal directory
  amrex::UtilCreateDirectory("temporals", 0755);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::string probesFileName = "temporals/probes";
    tmpProbesFile.open(
      probesFileName.c_str(),
      std::ios::out | std::ios::app | std::ios_base::binary);
    tmpProbesFile.precision(12);
    tmpProbesFile << "Tests \n";
    tmpProbesFile.flush();
  }
}

void
PeleLM::closeProbesFile()
{
  if (m_have_probes == 0) {
    return;
  }

  if (amrex::ParallelDescriptor::IOProcessor()) {
    tmpProbesFile.flush();
    tmpProbesFile.close();
  }
}
