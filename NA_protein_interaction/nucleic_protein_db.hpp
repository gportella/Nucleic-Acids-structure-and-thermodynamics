#include <gromacs/trajectoryanalysis.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

#define PI 3.14159265358979323846    /* PI = 4*atan(1)       */
#define CAL2JOULE (4.184)            /* id           */
#define ONE_4PI_EPS0 (FACEL * 0.1)   /* 1/(4*pi*e0)*/
#define ONE_4PI_EPS0_KCAL (33.20636) /* 1/(4*pi*e0)*/

using namespace gmx;

namespace gpc {

// Probably incomplete, might be better to read aminoacids.dat
// and residue type eventually, or make up another filename, such that
// it can be done by the user. Right now I will do it like that until
// I get the code ready
// This is how to find the filename
// #include <gromacs/utility/futil.h>
//  std::cout << "Aminoacids " << gmxlibfn("aminoacids.dat") << std::endl;
// You can open it with libopen

// Vector with all names allowed as nucleic acids
std::vector<const std::string> nuc_resnm = {
    "DA",  "DA3", "DA5", "DAB", "DALA", "DAN", "DC",  "DC3", "DC5", "DCN",
    "DG",  "DG3", "DG5", "DGN", "DT",   "DT3", "DT5", "DTN", "RA",  "RA3",
    "RA5", "RAN", "RC",  "RC3", "RC5",  "RCN", "RG",  "RG3", "RG5", "RGN",
    "RU",  "RU3", "RU5", "RUN", "NA5",  "NA3", "NA",  "NG5", "NG3", "NG",
    "NC5", "NC3", "NC",  "NU5", "NU3",  "NU",  "NT5", "NT3", "NT",  "FA5",
    "FA3", "FA",  "FG5", "FG3", "FG",   "FC5", "FC3", "FC",  "FU5", "FU3",
    "FU",  "FT5", "FT3", "FT",  "A",    "T",   "U",   "C",   "G",   "AMO",
    "GMO", "TMO", "CMO", "PMO", "DP5",  "DP",  "DP3", "URA", "ADE", "CYT",
    "GUA", "THY", "UR2", "AD2"};

// Vector with all names allowed as aminoacids
std::vector<const std::string> aa_resnm = {
    "ABU",  "ACE",  "ACE",  "AIB",  "ALA",  "ARG",  "ARGN",  "ASH",  "ASN",
    "ASN1", "ASP",  "ASP1", "ASPH", "CALA", "CARG", "CASN",  "CASP", "CCYN",
    "CCYX", "CGLN", "CGLU", "CGLY", "CHID", "CHIE", "CHIP",  "CILE", "CLEU",
    "CLYP", "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP",  "CTYR", "CVAL",
    "CYM",  "CYN",  "CYS",  "CYS1", "CYS2", "CYSH", "CYX",   "DA",   "DA3",
    "DA5",  "DAB",  "DALA", "DAN",  "GLH",  "GLN",  "GLU",   "GLUH", "GLY",
    "HID",  "HIE",  "HIP",  "HIS",  "HIS1", "HISA", "HISB",  "HISH", "HYP",
    "ILE",  "LEU",  "LYN",  "LYP",  "LYS",  "LYSH", "MELEU", "MET",  "MEVAL",
    "NAC",  "NALA", "NARG", "NASN", "NASP", "NCYN", "NCYX",  "NGLN", "NGLU",
    "NGLY", "NH2",  "NHE",  "NHID", "NHIE", "NHIP", "NILE",  "NLEU", "NLYP",
    "NME",  "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP",  "NTYR", "NVAL",
    "ORN",  "PHE",  "PHEH", "PHEU", "PHL",  "PRO",  "SER",   "THR",  "TRP",
    "TRPH", "TRPU", "TYR",  "TYRH", "TYRU", "VAL"};

std::vector<const std::string> Ade_sc_at = {"N1", "C2", "N3", "C4", "C5",
                                            "C6", "N7", "C8", "N9"};
std::vector<const std::string> Gua_sc_at = {"N1", "C2", "N3", "C4", "C5",
                                            "C6", "N7", "C8", "N9"};
std::vector<const std::string> Thy_sc_at = {"N1", "C2", "N3", "C4", "C5", "C6"};
std::vector<const std::string> Cyt_sc_at = {"N1", "C2", "N3", "C4", "C5", "C6"};

std::vector<const std::string> nucleic_at = {
    "C2",  "C4",  "C5", "C6", "C7",  "C8",  "H1",  "H2",  "H21", "H22", "H3",
    "H41", "H42", "H5", "H6", "H61", "H62", "H71", "H72", "H73", "H8",  "N1",
    "N2",  "N3",  "N4", "N6", "N7",  "N9",  "O2",  "O4",  "O6"};

std::vector<const std::string> nucleic_5p_resnm = {
    "DG5", "DA5", "DT5", "DC5", "U5", "A5", "C5", "G5", "RT5", "DU5"};
std::vector<const std::string> nucleic_3p_resnm = {
    "DG3", "DA3", "DT3", "DC3", "U3", "A3", "C3", "G3", "RT3", "DU3"};

std::vector<const std::string> AA_bb_at = {"C", "O", "N", "HN", "CA", "HA"};
}
