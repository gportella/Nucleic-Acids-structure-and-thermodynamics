#include "nucleic_protein.hpp"
#include <gromacs/trajectoryanalysis.h>
#include <gromacs/utility/futil.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <vector>

using namespace gmx;

class AnalysisTemplate : public TrajectoryAnalysisModule {
public:
  AnalysisTemplate();

  virtual void initOptions(IOptionsContainer *options,
                           TrajectoryAnalysisSettings *settings);

  virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                            const TopologyInformation &top);

  virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                            TrajectoryAnalysisModuleData *pdata);

  virtual void finishAnalysis(int nframes);
  virtual void writeOutput();

private:
  class ModuleData;

  std::string fnInter_;
  std::string fnInter_lj_;
  std::string fnInter_c_;
  std::string fnInter_p_;
  std::string fnInter_c_p_;
  std::string fnInter_lj_p_;
  double cutoff_;
  Selection refsel_;
  SelectionList sel_;
  AnalysisNeighborhood nb_;
  TrajectoryAnalysisModuleData *pdata;
  AnalysisData data_;
  gpc::InterActor inter_;
};

AnalysisTemplate::AnalysisTemplate() : cutoff_(0.0) {}

void AnalysisTemplate::initOptions(IOptionsContainer *options,
                                   TrajectoryAnalysisSettings *settings) {
  static const char *const desc[] = {
      "Computes the interaction energies between a reference (nucleic acids)",
      "and a selection (protein).",
      "Multiple selections could be easily implemented (the code is ready),",
      "it would only require adapting the output bit. "};

  settings->setHelpText(desc);

  options->addOption(FileNameOption("o")
                         .filetype(eftGenericData)
                         .outputFile()
                         .store(&fnInter_)
                         .required()
                         .defaultBasename("interaction")
                         .description("Interaction matrix"));

  options->addOption(FileNameOption("olj")
                         .filetype(eftGenericData)
                         .outputFile()
                         .store(&fnInter_lj_)
                         .required()
                         .defaultBasename("inter_lj")
                         .description("Interaction matrix with only LJ"));
  options->addOption(FileNameOption("oc")
                         .filetype(eftGenericData)
                         .outputFile()
                         .store(&fnInter_c_)
                         .required()
                         .defaultBasename("inter_coul")
                         .description("Interaction matrix with only Coulomb"));
  options->addOption(
      FileNameOption("o_p")
          .filetype(eftGenericData)
          .outputFile()
          .store(&fnInter_p_)
          .required()
          .defaultBasename("inter_proj")
          .description("Interaction profile for each NA residue"));

  options->addOption(
      FileNameOption("oc_p")
          .filetype(eftGenericData)
          .outputFile()
          .store(&fnInter_c_p_)
          .required()
          .defaultBasename("inter_coul_proj")
          .description("Interaction profile Coulomb for each NA residue"));

  options->addOption(
      FileNameOption("olj_p")
          .filetype(eftGenericData)
          .outputFile()
          .store(&fnInter_lj_p_)
          .required()
          .defaultBasename("inter_lj_proj")
          .description("Interaction profile LJ for each NA residue"));

  options->addOption(
      SelectionOption("reference")
          .store(&refsel_)
          .required()
          .description(
              "Reference nucleic acid group to calculate distances from."));

  options->addOption(
      SelectionOption("select")
          .storeVector(&sel_)
          .required()
          //.valueCount(1)
          .multiValue()
          .description("Protein group to calculate distances to."));

  options->addOption(DoubleOption("cutoff").store(&cutoff_).description(
      "Cutoff for neighbour calculation (0 = no cutoff)"));

  settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void AnalysisTemplate::initAnalysis(const TrajectoryAnalysisSettings &settings,
                                    const TopologyInformation &top /*top*/)

{
  nb_.setCutoff(cutoff_);
  inter_.populate(top, refsel_, sel_);
}

void AnalysisTemplate::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                    TrajectoryAnalysisModuleData *pdata) {
  const Selection &refsel = pdata->parallelSelection(refsel_);
  AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, refsel);

  /*! |brief
    * I removed the parallelSelection bit (see template.cpp example from
    * in the share/gromacs/template installation) because it
    * would imply changing InterActor to just care for one selection
    * and create an array of those, which could be done, but the parallelism
    * is actually not yet implemented according to the documentation anyhow
  */
  // TODO --> perhaps it would be best to have this as constructor?
  inter_.interact(nbsearch, refsel_, sel_);
}

void AnalysisTemplate::finishAnalysis(int frnr) {

  for (unsigned i = 0; i < sel_.size(); ++i) {
    for (unsigned j = 0; j < inter_.ref_num; ++j) {
      for (unsigned k = 0; k < inter_.sel_num[i]; ++k) {
        inter_.mat_inter[i][j][k].lj /= frnr;
        inter_.mat_inter[i][j][k].coul /= frnr;
      }
    }
  }
}

void AnalysisTemplate::writeOutput() {
  FILE *fp = NULL;
  FILE *fp_lj = NULL;
  FILE *fp_c = NULL;
  FILE *fp_proj = NULL;
  FILE *fp_c_proj = NULL;
  FILE *fp_lj_proj = NULL;
  for (unsigned g = 0; g < sel_.size(); ++g) {
    std::string prefix = "sel_" + std::to_string(g + 1) + "_";
    std::string filen = prefix + fnInter_;
    fp = fopen(filen.c_str(), "w");
    filen = prefix + fnInter_lj_;
    fp_lj = fopen(filen.c_str(), "w");
    filen = prefix + fnInter_c_;
    fp_c = fopen(filen.c_str(), "w");
    filen = prefix + fnInter_p_;
    fp_proj = fopen(filen.c_str(), "w");
    filen = prefix + fnInter_lj_p_;
    fp_c_proj = fopen(filen.c_str(), "w");
    filen = prefix + fnInter_c_p_;
    fp_lj_proj = fopen(filen.c_str(), "w");
    real sum_lj = 0;
    real sum_c = 0;
    for (unsigned j = 0; j < inter_.ref_num; ++j) {
      auto ref_resid = inter_.ref_nuc[j].resid;
      auto ref_resname = inter_.ref_nuc[j].resname;
      for (unsigned k = 0; k < inter_.sel_num[g]; ++k) {
        auto sel_resid = inter_.sel_aa[g][k].resid;
        auto sel_resname = inter_.sel_aa[g][k].resname;
        fprintf(fp, "%s %d %s %d %7.4f\n", ref_resname.c_str(), ref_resid,
                sel_resname.c_str(), sel_resid,
                inter_.mat_inter[g][j][k].coul + inter_.mat_inter[g][j][k].lj);
        fprintf(fp_lj, "%s %d %s %d %7.4f\n", ref_resname.c_str(), ref_resid,
                sel_resname.c_str(), sel_resid, inter_.mat_inter[g][j][k].lj);
        fprintf(fp_c, "%s %d %s %d %7.4f\n", ref_resname.c_str(), ref_resid,
                sel_resname.c_str(), sel_resid, inter_.mat_inter[g][j][k].coul);
        sum_c += inter_.mat_inter[g][j][k].coul;
        sum_lj += inter_.mat_inter[g][j][k].lj;
      }
      fprintf(fp_proj, "%s %d %7.4f\n", ref_resname.c_str(), ref_resid,
              sum_c + sum_lj);
      fprintf(fp_c_proj, "%s %d %7.4f\n", ref_resname.c_str(), ref_resid,
              sum_c);
      fprintf(fp_lj_proj, "%s %d %7.4f\n", ref_resname.c_str(), ref_resid,
              sum_lj);
      fprintf(fp, "\n");
      fprintf(fp_lj, "\n");
      fprintf(fp_c, "\n");
      sum_c = 0;
      sum_lj = 0;
    }
    fclose(fp);
    fclose(fp_lj);
    fclose(fp_proj);
    fclose(fp_c_proj);
    fclose(fp_lj_proj);
    fclose(fp_c);
  }
}

/*! \brief
 * The main function for the analysis template.
 */
int main(int argc, char *argv[]) {
  return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<AnalysisTemplate>(
      argc, argv);
}
