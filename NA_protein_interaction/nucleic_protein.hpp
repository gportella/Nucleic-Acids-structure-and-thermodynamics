#include "nucleic_protein_db.hpp"
#include <gromacs/trajectoryanalysis.h>
#include <gromacs/utility/futil.h>
#include <gromacs/utility/smalloc.h>
#include <iostream>
#include <map>
#include <stddef.h>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <vector>

// YOU COULD CONSIDER USING THE FOLLOWING
// Inits the IDs for use with SelectionPosition::mappedId() for grouping.
// (SEARCH FOR MAPPEDID, READ WHAT IT SAYS ABOUT USING type == INDEX_RES
// THIS MIGHT BE THE SOLUTION...
//
// You compute everything atomwise, based on the index, store it somewhere
// in some array, and finally postprocess everything to put it in the right
// place
//

namespace gpc {

struct Tag_nuc {};
struct Tag_aa {};

template <typename Tstr, typename Tstrvect>
bool in_array(Tstr const &value, Tstrvect const &array) {
  return std::find(array.begin(), array.end(), value) != array.end();
}

class AtomInfo {

public:
  unsigned ind;
  real q;
  real c6;
  real c12;
};

//////////////////////////////////////////////////////////////////////////////
// Class  ResInfo (could actually just be a typdef struct...)
//////////////////////////////////////////////////////////////////////////////
class ResInfo {
  // See if you should declare it private and friend InterActor
  //  friend class InterActor;
public:
  std::vector<AtomInfo> a_inf;
  std::string resname; // For convenience
  unsigned resid;
};

//////////////////////////////////////////////////////////////////////////////
// Class InterActor
//
// -> Populates data structures containing q and LJ parameters
// -> Computes Interactions
//
// Dependencies:
// Inherits the analysis settings from TrajectoryAnalysisSettings
// Inherits from ResInfo... This might not be brilliant
//
//////////////////////////////////////////////////////////////////////////////

class InterActor : TrajectoryAnalysisSettings, ResInfo {

public:
  InterActor();
  ~InterActor();

  template <typename Ttop, typename Tref, typename Tsel>
  void populate(Ttop const &top, Tref const &refsel, Tsel const &sel);
  void interact(AnalysisNeighborhoodSearch nbseach, Selection ref_,
                SelectionList sel_);

  std::vector<ResInfo> ref_nuc;             // Store the reference nucleic acid
  std::vector<std::vector<ResInfo>> sel_aa; // Store each AA for each selection
  std::map<unsigned, std::tuple<unsigned, unsigned>> mmap_ref;
  std::vector<std::map<unsigned, std::tuple<unsigned, unsigned>>> mmap_sel;
  typedef struct {
    float lj;
    float coul;
  } energy;
  energy ***mat_inter;
  unsigned ref_num;
  std::vector<unsigned> sel_num;

private:
  template <typename Ttop, typename Tref>
  unsigned fill_pop_r_(Ttop const &top, Tref const &ref, Tag_nuc const &);
  template <typename Ttop, typename Tsel>
  std::vector<unsigned> fill_pop_s_(Ttop const &top, Tsel const &sel,
                                    Tag_aa const &);

  real e_coul(real inv_d, real q1, real q2);
  real e_lj(real inv_d, real sig, real eps);
};

//////////////////////////////////////////////////////////////////////////////
InterActor::InterActor() {
  mat_inter = NULL;
  ref_num = 0;
}
// To be determined
InterActor::~InterActor() {
  for (unsigned i = 0; i < sel_num.size(); ++i) {
    for (unsigned j = 0; j < ref_num; ++j) {
      sfree(mat_inter[i][j]);
    }
    sfree(mat_inter[i]);
  }
  sfree(mat_inter);
}
//////////////////////////////////////////////////////////////////////////////
// Call this function to populate InterActor data
// We allow for one reference group and N selections, though it is
// capped at two selections by the gmx parsing function right now
//////////////////////////////////////////////////////////////////////////////

template <typename Ttop, typename Tref, typename Tsel>
void InterActor::populate(Ttop const &top, Tref const &refsel,
                          Tsel const &sel) {

  // actually you could fill the ref_nuc and sel_aa inside the function
  // just like you do with maps, and no need to return anything
  ref_num = fill_pop_r_(top, refsel, Tag_nuc());
  sel_num = fill_pop_s_(top, sel, Tag_aa());
  /*
  for (unsigned i = 0; i < sel_num.size(); ++i) {
    std::cout << "Found SELNUM- " << i << " : " << sel_num[i] << " residues"
              << std::endl;
  }
  */
  snew(mat_inter, sel_num.size());
  for (unsigned i = 0; i < sel_num.size(); ++i) {
    snew(mat_inter[i], ref_num);
    for (unsigned j = 0; j < ref_num; ++j) {
      snew(mat_inter[i][j], sel_num[i]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
// Do the population job, one template per nucleotides (reference) and the
// other for aminoacids (aa, selections)
// Keep a map (array of tuples) such that given an atom index we can know
// in with ResInfo element we have the data
//////////////////////////////////////////////////////////////////////////////

template <typename Ttop, typename Tref>
unsigned InterActor::fill_pop_r_(Ttop const &top, Tref const &ref,
                                 Tag_nuc const &) {
  auto t_atoms = top.topology()->atoms;
  auto resinfo = t_atoms.resinfo;
  auto iparams = top.topology()->idef.iparams;
  auto ntype = top.topology()->idef.atnr;

  unsigned ct_at = 0;
  unsigned ct_res = 0;
  unsigned rsnr_b4 = 0;
  ResInfo r;
  AtomInfo ainf;

  std::string resname_b4 = "";
  for (unsigned i = 0; i < ref.atomCount(); ++i) {
    unsigned atomindex = ref.atomIndices()[i];
    std::string resname = *(resinfo[t_atoms.atom[atomindex].resind].name);
    unsigned rsnr = resinfo[t_atoms.atom[atomindex].resind].nr;
    if (i == 0) {
      rsnr_b4 = rsnr;
      resname_b4 = resname;
    }
    // We fill up residuewise
    r.resname = resname_b4;
    r.resid = rsnr_b4;
    if (in_array(resname, nuc_resnm)) {
      // If there has been a residue change, attach atoms
      if (rsnr != rsnr_b4) {
        ref_nuc.push_back(r);
        // delete r...
        r = ResInfo();
        ct_at = 0;
        ct_res++;
      }
      // if not, then record or keep recording
      ainf.ind = atomindex;
      ainf.q = t_atoms.atom[atomindex].q;
      auto itype = t_atoms.atom[atomindex].type;
      ainf.c6 = iparams[itype * (ntype + 1)].lj.c6;
      ainf.c12 = iparams[itype * (ntype + 1)].lj.c12;
      r.a_inf.push_back(ainf);
      mmap_ref[atomindex] =
          std::make_tuple(ct_res, ct_at); // fill up tuple to map
      ct_at++;
      if (i == ref.atomCount() - 1) {
        ref_nuc.push_back(r);
        // delete r
        r = ResInfo();
        ct_res++;
      }
    }
    rsnr_b4 = rsnr;
    resname_b4 = resname;
  }
  return ct_res;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename Ttop, typename Tsel>
std::vector<unsigned> InterActor::fill_pop_s_(Ttop const &top, Tsel const &sel,
                                              Tag_aa const &) {

  auto t_atoms = top.topology()->atoms;
  auto resinfo = t_atoms.resinfo;
  auto iparams = top.topology()->idef.iparams;
  auto ntype = top.topology()->idef.atnr;

  unsigned ct_at = 0;
  unsigned ct_res = 0;
  unsigned rsnr_b4 = 0;

  std::vector<unsigned> sel_res_n;
  std::vector<ResInfo> ll_sel_aa;
  std::map<unsigned, std::tuple<unsigned, unsigned>> l_mmap_sel;
  std::string resname_b4 = "";

  ResInfo r;
  for (unsigned g = 0; g < sel.size(); ++g) {
    for (unsigned i = 0; i < sel[g].atomCount(); ++i) {
      unsigned atomindex = sel[g].atomIndices()[i];
      std::string resname = *(resinfo[t_atoms.atom[atomindex].resind].name);
      unsigned rsnr = resinfo[t_atoms.atom[atomindex].resind].nr;
      if (i == 0) {
        rsnr_b4 = rsnr;
        resname_b4 = resname;
      }
      // We fill up residuewise
      r.resname = resname_b4;
      r.resid = rsnr_b4;
      if (in_array(resname, aa_resnm)) {
        if (rsnr != rsnr_b4) {
          ll_sel_aa.push_back(r);
          // delete r...
          r = ResInfo();
          ct_at = 0;
          ct_res++;
        }
        // if not, then record or keep recording
        AtomInfo ainf;
        ainf.ind = atomindex;
        ainf.q = t_atoms.atom[atomindex].q;
        auto itype = t_atoms.atom[atomindex].type;
        ainf.c6 = iparams[itype * (ntype + 1)].lj.c6;
        ainf.c12 = iparams[itype * (ntype + 1)].lj.c12;
        r.a_inf.push_back(ainf);
        l_mmap_sel[atomindex] =
            std::make_tuple(ct_res, ct_at); // fill up tuple to map
        ct_at++;
        // if it is the last element record it, don't wait
        if (i == sel[g].atomCount() - 1) {
          ll_sel_aa.push_back(r);
          // delete r...
          r = ResInfo();
        }
      }
      rsnr_b4 = rsnr;
      resname_b4 = resname;
    }
    sel_aa.push_back(ll_sel_aa);
    mmap_sel.push_back(l_mmap_sel);
    sel_res_n.push_back(ct_res);
    ct_at = 0;
    ct_res = 0;
  }
  return sel_res_n;
}

//////////////////////////////////////////////////////////////////////////////
// This will compute the interactions using a neighbour list
// Requires a neighbour search object already initialized
//////////////////////////////////////////////////////////////////////////////

void InterActor::interact(AnalysisNeighborhoodSearch nbsearch,
                          Selection refsel_, SelectionList sel_) {

  for (unsigned g = 0; g < sel_.size(); ++g) {

    AnalysisNeighborhoodPairSearch pairSearch =
        nbsearch.startPairSearch(sel_[g]);
    AnalysisNeighborhoodPair pair;

    while (pairSearch.findNextPair(&pair)) {

      auto ref_index = pair.refIndex();
      auto pair_index = pair.testIndex();
      auto dd = refsel_.atomIndices()[ref_index];
      auto rr = sel_[g].atomIndices()[pair_index];

      if (mmap_sel[g].find(rr) == mmap_sel[g].end()) {
        std::cout << "Woops!, key sel " << rr << " not found " << std::endl;
      } else if (mmap_ref.find(dd) == mmap_ref.end()) {
        std::cout << "Woops!, key ref" << dd << " not found " << std::endl;
      } else {
        std::tuple<unsigned, unsigned> t_r = mmap_ref[dd];
        std::tuple<unsigned, unsigned> t_s = mmap_sel[g][rr];

        real inv_d = gmx::invsqrt(pair.distance2());

        // recall that in map tuple the first entry is the residue, and
        // the second is the atom
        unsigned r_r = std::get<0>(t_r);
        unsigned a_r = std::get<1>(t_r);
        unsigned r_s = std::get<0>(t_s);
        unsigned a_s = std::get<1>(t_s);

        // Coulomb
        real q1 = ref_nuc[r_r].a_inf[a_r].q;
        real q2 = sel_aa[g][r_s].a_inf[a_s].q;
        real e_pot = e_coul(inv_d, q1, q2);
        mat_inter[g][r_r][r_s].coul += e_pot;

        // recover c6/c12, convert to sigma/epsilon
        // combine with LB rule and compute LennarJones
        real c6_1 = ref_nuc[r_r].a_inf[a_r].c6;
        real c12_1 = ref_nuc[r_r].a_inf[a_r].c12;
        real c6_2 = sel_aa[g][r_s].a_inf[a_s].c6;
        real c12_2 = sel_aa[g][r_s].a_inf[a_s].c12;

        if (c6_1 != 0 && c12_1 != 0 && c6_2 != 0 && c12_2 != 0) {
          real sig_1 = gmx::sixthroot(c12_1 / c6_1);
          real eps_1 = 0.25 * c6_1 / gmx::power6(sig_1);
          real sig_2 = gmx::sixthroot(c12_2 / c6_2);
          real eps_2 = 0.25 * c6_2 / gmx::power6(sig_2);
          // mixing rules LB are used in amber/charmm
          // compute LJ
          real sig_12 = 0.5 * (sig_1 + sig_2);
          real eps_12 = std::sqrt(eps_1 * eps_2);
          real lj = e_lj(inv_d, sig_12, eps_12);
          mat_inter[g][r_r][r_s].lj += lj;
        }
      }
    }
  }
}

/*! |brief
 * Compute Lennard Jones given the right sigma and eps
 * Sigma and epsilon should have been properly combined
 * according to the mixing/combination rules
 * For this function to work ok remember that
 * sig_12 and eps_12 should be already based on combination rules
 * sig_ij = 1/2 (sig_i + sig_j)
 * eps_ij = sqrt(eps_i * eps_j)
 */
real InterActor::e_lj(real inv_d, real sig, real eps) {
  real sig_over_r_6 = sig * inv_d * sig * inv_d * sig * inv_d;
  sig_over_r_6 *= sig * inv_d * sig * inv_d * sig * inv_d;
  real sig_over_r_12 = sig_over_r_6 * sig_over_r_6;

  return (4 * eps * ((sig_over_r_12) - (sig_over_r_6)));
}
/*! |brief
 * Compute electrostatic potential
 */
real InterActor::e_coul(real inv_d, real q1, real q2) {
  return (inv_d * (q1 * q2));
}

//////////////////////////////////////////////////////////////////////////////
// End of class InterActor
//////////////////////////////////////////////////////////////////////////////

} // end of namespace gpc

////////////////////////////////Tricky trick////////////////////////////////////
