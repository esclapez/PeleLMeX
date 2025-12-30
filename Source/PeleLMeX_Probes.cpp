#include "PeleLMeX_Index.H"
#include <PeleLMeX.H>

void
PeleLM::probesSetup()
{
  std::string pele_prefix = "peleLM";
  amrex::ParmParse pp(pele_prefix);

  m_nprobes = pp.countval("probes");
  if (m_nprobes <= 0) {
      return;
  }

  if (m_nprobes > 0) {
    // Input data: probes coordinates and labels
    m_probes_coord.resize(m_nprobes);
    m_probes_name.resize(m_nprobes);

    // Internal data, every rank has full
    // knowledge of where all the probes are
    m_probes_cell.resize(m_nprobes);
    m_probes_weight.resize(m_nprobes);
    m_tmpProbesFile.resize(m_nprobes);

    // Each rank keep a list of the probes it owns
    // thus up to a full m_nprobes if all the probes on a given rank
    m_myprobes.resize(m_nprobes);
    for (auto &id: m_myprobes) {
        id = -1;
    }

    std::string probe_prefix = pele_prefix + ".probes";
    amrex::ParmParse ppp(probe_prefix);
    for (int n = 0; n < m_nprobes; ++n) {
      pp.get("probes", m_probes_name[n], n);
      amrex::Vector<amrex::Real> coord(AMREX_SPACEDIM);
      ppp.getarr(m_probes_name[n], coord, 0, AMREX_SPACEDIM);
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_probes_coord[n][idim] = coord[idim];
      }
    }

    ppp.query("int", m_probes_int);

    // TODO: For now only dump the state variables
    int n_probes_var = NVAR;
    //n_probes_var = ppp.countval("variables");
    if (n_probes_var > 0) {
      m_probes_vars.resize(n_probes_var);
      //for (int ivar = 0; ivar < n_probes_var; ++ivar) {
      //  ppp.get("variables", m_probes_vars[ivar], ivar);
      //  // Need to check the availability of the variables, but the derived
      //  // arent initialized yet
      //}
    } else {
        amrex::Abort("Using probes explicitly requires to probe of list in peleLM.probes.variables");
    }
  }
}

void
PeleLM::updateProbesCellsAndWeights()
{
  if (m_nprobes < 0) {
    return;
  }

  // Have rank close the ofstream to the probes they previously owned
  closeProbesFile();

  // Reset the list of probes owned
  for (auto &id: m_myprobes) {
      id = -1;
  }
  m_mynprobes = 0;

  // Set all the probes levels at -1
  for (auto& probe_arr : m_probes_cell) {
      probe_arr[0] = -1;
  }

  // Find the cell, level and proc the probe is on
  for (int lev = finest_level; lev >= 0; --lev) {

    for (int i = 0; i < m_nprobes; ++i) {

      // Probe already located on finer level, skip
      if (m_probes_cell[i][0] >= 0) {
          continue;
      }

      // Find cell index
      const auto dxinv = geom[lev].InvCellSizeArray();
      amrex::IntVect iv;
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          iv[idim] = static_cast<int>(m_probes_coord[i][idim] * dxinv[idim]);
      }

      // Probe is not on this level
      if (!grids[lev].contains(iv)) {
          continue;
      }

      // Locate box
      int box_id = -1;
      for (int k = 0; k < grids[lev].size(); ++k) {
          if (grids[lev][k].contains(iv)) {
              box_id = k;
              break;
          }
      }

      // Update the data
      m_probes_cell[i][0] = lev;
      m_probes_cell[i][1] = box_id;
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_probes_cell[i][idim+2] = iv[idim];
      }

      // Set the weight to 1.0 on the center and 0.0
      // on surrounding cells for now
      // We might interpolate later
      m_probes_weight[i].fill(0.0);
      m_probes_weight[i][0] = 1.0;

      // Rank
      int owner_rank = dmap[lev][box_id];
      if (owner_rank == amrex::ParallelDescriptor::MyProc()) {
          m_myprobes[m_mynprobes] = i;
          m_mynprobes += 1; 
      }
    }
  }

  // Each rank re-open the ofstreams it now owns
  openProbesFile();
}

void
PeleLM::writeProbes()
{
  if (m_mynprobes < 1) {
      return;
  }

  for (int k=0; k < m_mynprobes; ++k) {
    // Probe index
    int pidx = m_myprobes[k];
    // Get the data to write
    amrex::Vector<amrex::Real> pdata(NVAR, 0.0);

    // TODO: GPU
    auto* ldata_p = getLevelDataPtr(m_probes_cell[pidx][0], AmrNewTime);
    int ip = m_probes_cell[pidx][2];
    int jp = m_probes_cell[pidx][3];
    int kp = (AMREX_SPACEDIM == 3) ? m_probes_cell[pidx][4] : 0;
    auto state_celldata = ldata_p->state.array(m_probes_cell[pidx][1]).cellData(ip,jp,kp);
    for (int i = 0; i < NVAR; ++i) {
        pdata[i] = state_celldata[i];
    }

    m_tmpProbesFile[pidx] << m_nstep << "," << m_cur_time << "," << m_dt; // Time
    for (const auto &v: pdata) {
        m_tmpProbesFile[pidx] << "," << v;
    }
    m_tmpProbesFile[pidx] << "\n";
    m_tmpProbesFile[pidx].flush();
  }
}

void
PeleLM::openProbesFile()
{
  if (m_nprobes < 0) {
    return;
  }

  // Create the temporal directory
  amrex::UtilCreateDirectory("temporals", 0755);

  for (int k=0; k < m_mynprobes; ++k) {
    std::string probesFileName = "temporals/probes_" + m_probes_name[m_myprobes[k]];
    m_tmpProbesFile[m_myprobes[k]].open(probesFileName.c_str(),
            std::ios::out | std::ios::app | std::ios_base::binary);
    m_tmpProbesFile[m_myprobes[k]].precision(12);
  }
}

void
PeleLM::closeProbesFile()
{
  if (m_nprobes < 0) {
    return;
  }

  for (int k=0; k < m_mynprobes; ++k) {
    m_tmpProbesFile[m_myprobes[k]].flush();
    m_tmpProbesFile[m_myprobes[k]].close();
  }
}
