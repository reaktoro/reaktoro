//--------------------------------------------------------------------
// $Id: node_format.cpp 879 2013-10-10 14:47:33Z dmitrieva $
//
/// \file node_format.cpp
/// Interface for writing/reading DBR and DCH I/O files of GEMS3K
/// Works  with DATACH and DATABR structures
//
// Copyright (c) 2006-2012 S.Dmytriyeva, D.Kulik
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------

#include <iomanip>
#include  <iostream>
#include "io_arrays.h"
#include "node.h"
#include "gdatastream.h"

#ifdef NODEARRAYLEVEL
#include "nodearray.h"
#endif

//extern bool _comment;
extern const char* _GEMIPM_version_stamp;

//===============================================================
// in the arrays below, the first field of each structure contains a string
// which is put into <> to comprise a data object tag, e.g. <IterDone>, in
// free text input files. The second field (0 or 1) denotes whether the data
// object can be skipped from the file (0) and default value(s) can be used,
// or (1) the data object must be always present in the file. The third
// field is used internally and must be set to 0 here. The fourth field contains
// the text of the comment for this data object, optionally written into the
// text-format output DBR or DCH file.
//
outField DataBR_fields[f_lga+1/*60*/] =  {
  { "NodeHandle",  0, 0, 1, "# NodeHandle: Node identification handle"},
  { "NodeTypeHY",  0, 0, 1, "# NodeTypeHY:  Node type code (hydraulic), not used on TNode level; see typedef NODETYPE" },
  { "NodeTypeMT",  0, 0, 1, "# NodeTypeMT:  Node type (mass transport), not used on TNode level; see typedef NODETYPE" },
  { "NodeStatusFMT",  1, 0, 1, "# NodeStatusFMT:  Node status code in FMT part, not used on TNode level; see typedef NODECODEFMT" },
  { "NodeStatusCH",  1, 0, 1, "# NodeStatusCH: Node status code and control in GEM input and output; see typedef NODECODECH"},
  { "IterDone",  0, 0, 1, "# IterDone:  Number of iterations performed by GEM IPM in the last run (GEM output)" },
  { "TK",   1, 0, 1, "# TK: Node temperature T, Kelvin. This value must always be provided (GEM input)"  },
  { "P",    1, 0, 1, "# P:  Node Pressure P, Pa. This value must always be provided (GEM input)" },
  { "Vs",   0, 0, 1, "# Vs: Volume V of reactive subsystem, m3 (GEM output)" },
  { "Vi",   0, 0, 1, "# Vi: Volume of inert subsystem, m3 (mass transport)" },
  { "Ms",   0, 0, 1, "# Ms: Mass of reactive subsystem, kg (GEM output)" },
  { "Mi",   0, 0, 1, "# Mi: Mass of inert subsystem, kg (mass transport)" },
  { "Hs",   0, 0, 1, "# Hs: Total enthalpy of reactive subsystem, J (reserved)" },
  { "Hi",   0, 0, 1, "# Hi: Total enthalpy of inert subsystem, J (reserved, mass transport) " },
  { "Gs",   0, 0, 1, "# Gs: Total Gibbs energy of the reactive subsystem, J/(RT) (GEM output)" },
  { "IS",   0, 0, 1, "# IS: Effective aqueous ionic strength, molal (GEM output)" },
  { "pH",   0, 0, 1, "# pH: pH of aqueous solution in molal activity scale (GEM output)" },
  { "pe",   0, 0, 1, "# pe: pe of aqueous solution in molal activity scale (GEM output)" },
  { "Eh",   0, 0, 1, "# Eh: Eh of aqueous solution, V (GEM output)" },
  { "Tm",   0, 0, 1, "# Tm: Actual total simulation time, s (kinetics, metastability, transport)" },
  { "dt",   0, 0, 1, "# dt: Actual time step, s (kinetics, metastability, transport)" },
//#ifdef NODEARRAYLEVEL
// Scalar parameters below are only used at TNodeArray level
  { "Dif",  0, 0, 1, "# Dif: General diffusivity of disolved matter, m2/s (mass transport)" },
  { "Vt",   0, 0, 1, "# Vt: Total volume of the node, m3 (mass transport)" },
  { "vp",   0, 0, 1, "# vp: Advection velocity in pores, m/s (mass transport)" },
  { "eps",  0, 0, 1, "# eps: Effective (actual) porosity normalized to 1 (mass transport)" },
  { "Km",   0, 0, 1, "# Km: Actual permeability, m2 (mass transport)" },
  { "Kf",   0, 0, 1, "# Kf: Actual Darcy`s constant, (m2/s (mass transport)" },
  { "S",    0, 0, 1, "# S: Specific storage coefficient, dimensionless (mass transport)" },
  { "Tr",   0, 0, 1, "# Tr: Transmissivity, m2/s (mass transport)" },
  { "h",    0, 0, 1, "# h: Actual hydraulic head (hydraulic potential), m (mass transport)" },
  { "rho",  0, 0, 1, "# rho: Actual carrier density for density-driven flow, kg/m3 (mass transport)" },
  { "al",   0, 0, 1, "# al: Specific longitudinal dispersivity of porous media, m (mass transport)" },
  { "at",   0, 0, 1, "# at: Specific transversal dispersivity of porous media, m (mass transport)" },
  { "av",   0, 0, 1, "# av: Specific vertical dispersivity of porous media, m (mass transport)" },
  { "hDl",  0, 0, 1, "# hDl: Hydraulic longitudinal dispersivity, m2/s (mass transport)" },
  { "hDt",  0, 0, 1, "# hDt: Hydraulic transversal dispersivity, m2/s (mass transport)" },
  { "hDv",  0, 0, 1, "# hDv: Hydraulic vertical dispersivity, m2/s (mass transport)" },
  { "nto",  0, 0, 1, "# nto: Tortuosity factor, dimensionless (mass transport)" },
//#endif
 // dynamic arrays (52-38=14)
  { "bIC",  1, 0, nICbi, "# bIC: Bulk composition of reactive subsystem (main GEM input), moles of ICs [nICb]" },
  { "rMB",  0, 0, nICbi, "\n# rMB: Mass balance residuals, moles (GEM output) [nICb]" },
  { "uIC",  0, 0, nICbi, "\n# uIC: Chemical potentials of ICs in equilibrium (dual solution), J/(RT) (GEM output) [nICb]" },
  { "xDC",  0, 0, nDCbi, "# xDC: Speciation - amounts of DCs in equilibrium (primal solution), moles (GEM output/input) [nDCb]" },
  { "gam",  0, 0, nDCbi, "\n# gam: Activity coefficients of DCs in their respective phases (GEM output/input) [nDCb]" },
  { "dll",  0, 0, nDCbi, "\n# dll: Lower metastability restrictions on amounts of DCs, moles (GEM input) [nDCb]" },
  { "dul",  0, 0, nDCbi, "\n# dul: Upper metastability constraints on amounts of DCs, moles (GEM input) [nDCb]" },
  { "aPH",  0, 0, nPHbi, "# aPH: Specific surface areas of phases, m2/kg (GEM input) [nPHb]" },
  { "xPH",  0, 0, nPHbi, "\n# xPH: Amounts of phases in equilibrium state, moles (GEM output) [nPHb]" },
  { "vPS",  0, 0, nPSbi, "\n# vPS: Volumes of multicomponent phases, m3 (GEM output) [nPSb]" },
  { "mPS",  0, 0, nPSbi, "\n# mPS: Masses of multicomponent phases, kg (GEM output) [nPSb]" },
  { "bPS",  0, 0, nPSbnICbi, "\n\n# bPS: Bulk elemental compositions of multicomponent phases, moles (GEM output) [nPSb*nICb]"},
  { "xPA",  0, 0, nPSbi, "\n# xPA: Amount of carrier (sorbent or solvent) in multicomponent phases, moles (GEM output) [nPSb]" },
  { "bSP",  0, 0, nICbi, "\n# bSP: Output bulk composition of the equilibrium solid part of the system, moles " },
    { "amru",  0, 0, nPSbi, "\n# amru: Upper AMRs on masses of DCs (kg) [nPSb]  " },
    { "amrl",  0, 0, nPSbi, "\n# amrl: Lower AMRs on masses of DCs (kg) [nPSb]" },

// only for VTK format output
    { "mPH",  0, 0, nPHbi, "# mPH: Masses of phases in equilibrium, kg [nPHb]" },
    { "vPH",  0, 0, nPHbi, "# vPH: Volumes of phases in equilibrium, m3 [nPHb]" },
    { "m_t",  0, 0, nICbi, "# m_t: Total dissolved molality of independent components, m [nICb]" },
    { "con",  0, 0, nDCbi, "# con: DC concentrations in phases (molal, mole fraction) [nDCb]" },
    { "mju",  0, 0, nDCbi, "# mju: DC chemical potentials in equilibrium, J/mol [nDCb]" },
    { "lga",  0, 0, nDCbi, "# lga: DC activities in equilibrium, in log10 scale [nDCb]" }
};

outField DataCH_static_fields[14] =  {
  { "nIC",   1, 0, 0, "# nIC: Number of Independent Components (usually chemical elements and charge)" },
  { "nDC",   1, 0, 0, "# nDC: Number of Dependent Components (chemical species made of Independent Components)" },
  { "nPH",   1, 0, 0, "# nPH: Number of phases (into which Dependent Components are grouped)" },
  { "nPS",   1, 0, 0, "# nPS: Number of phases-solutions (multicomponent phases) <= nPH" },
  { "nDCs",  1, 0, 0, "# nDCs: Number of Dependent Components in phases-solutions <= nDC" },
  { "nICb",  1, 0, 0, "# nICb: Number of ICs kept in the DBR file and DATABR memory structure (<= nIC)" },
  { "nDCb",  1, 0, 0, "# nDCb: Number of DCs kept in the DBR file and DATABR memory structure (<=nDC)"  },
  { "nPHb",  1, 0, 0, "# nPHb: Number of phases kept in the DBR file and DATABR structure (<=nPH)" },
  { "nPSb",  1, 0, 0, "# nPSb: Number of phases-solutions kept in the DBR file and DATABR structure (<=nPS)" },
  { "nTp",   1, 0, 0, "# nTp: Number of temperature grid points in lookup arrays for data interpolation, >=1" },
  { "nPp",   1, 0, 0, "# nPp: Number of pressure grid points in lookup arrays for data interpolation, >=1" },
  { "iGrd",  1, 0, 0, "# iGrd: Flag for allocation of array of diffusition coefficients in DATACH structure (DCH file)" },
  { "fAalp", 1, 0, 0, "# fAalp: Flag for keeping specific surface areas of phases in DATABR structure (1) or ignoring them (0)" },
  { "mLook", 1, 0, 0, "# mLook: Lookup mode: 0 interpolation over nTp*nPp grid; 1 data for T,P pairs, no interpolation"}
};

outField DataCH_dynamic_fields[30] =  { //+4
   { "xic",   1, 0, 0, "# xIC: DATACH access index list for ICs kept in the DATABR structure and in DBR files [nICb]" },
   { "xdc",   1, 0, 0, "# xDC: DATACH access index list of DCs kept in the DATABR  structure and in DBR files [nDCb]" },
   { "xph",   1, 0, 0, "# xPH: DATACH access index list for Phases kept in the DATABR structure and in DBR files [nPHb]" },
   { "ICNL",  1, 0, 0, "\n# ICNL: Name list of ICs (Independent Components, up to 5 characters per name) [nIC]" },
   { "ccIC",  1, 0, 0, "# ccIC: Class codes of ICs (Independent Components) [nIC]" },
   { "ICmm",  1, 0, 0, "# ICmm: Atomic (molar) masses of ICs,  kg/mol [nIC]" },
   { "DCNL",  1, 0, 0, "\n# DCNL: Name list of DCs (Dependent Components, up to 16 characters per name) [nDC]" },
   { "ccDC",  1, 0, 0, "# ccDC: Class codes of DCs (Dependent Components) [nDC]" },
   { "DCmm",  0, 0, 0, "\n# DCmm: Molar masses of DCs, kg/mol [nDC]" },
   { "PHNL",  1, 0, 0, "# PHNL: List of Phase names (up to 16 characters per name) [nPH]" },
   { "ccPH",  1, 0, 0, "# ccPH: Codes of phase aggregate state [nPH]" },
   { "nDCinPH",  1, 0, 0, "# nDCinPH: Number of DCs included in each phase [nPH]" },
   { "A",     1, 0, 0, "# A: Stoichiometry matrix A (expanded formulae) for DCs [nDC*nIC]"},
   { "Ttol",  0, 0, 0, "# Ttol: Tolerance for the temperature interpolation, K" },
   { "TKval", 1, 0, 0, "# TKval: Temperature values, K for lookup arrays of thermodynamic data [nTp]" },
   { "Ptol",  0, 0, 0, "# Ptol: Tolerance for the pressure interpolation, Pa" },
   { "Pval",  1, 0, 0, "# Pval: Pressure values, Pa for lookup arrays of thermodynamic data [nPp]" },
   { "denW",  1, 0, 0, "\n# denW: Look-up array for the density of water-solvent, kg/m3, and its derivatives [5*nPp*nTp]" },
   { "denWg", 1, 0, 0, "\n# denWg: Look-up array for the density of water vapour, kg/m3, and its derivatives [5*nPp*nTp]" },
   { "epsW",  1, 0, 0, "\n# epsW: Look-up array for the dielectric constant of water-solvent and its derivatives [5*nPp*nTp]" },
   { "epsWg", 1, 0, 0, "\n# epsWg: Look-up array for the dielectric constant of water vapour and its derivatives [5*nPp*nTp]" },
//   { "visW",  1, 0, 0 },
   { "V0",    1, 0, 0, "\n# V0: Look-up array for DC (standard) molar volumes, J/Pa [nDC*nPp*nTp]" },
   { "G0",    1, 0, 0, "\n# G0: Look-up array for DC molar Gibbs energy function g(T,P), J/mol [nDC*nPp*nTp]" },
   { "H0",    0, 0, 0, "\n# H0: Look-up array for DC molar enthalpy h(T,P), J/mol [nDC*nPp*nTp]" },
   { "S0",    0, 0, 0, "\n# S0: Look-up array for DC absolute entropy S(T,P), J/K/mol [nDC*nPp*nTp] " },
   { "Cp0",   0, 0, 0, "\n# Cp0: Look-up array for DC heat capacity Cp(T,P), J/K/mol [nDC*nPp*nTp]" },
   { "A0",    0, 0, 0, "\n# A0: reserved: Look-up array for DC Helmholtz energy function, J/mol [nDC*nPp*nTp]" },
   { "U0",    0, 0, 0, "\n# U0: reserved: Look-up array for DC internal energy function, J/mol [nDC*nPp*nTp]" },
   { "DD",    0, 0, 0, "\n# DD: reserved: Look-up array for DC diffusion coefficients [nDC*nPp*nTp]" },
   { "Psat",  0, 0, 0, "# Psat: Pressure Pa at saturated H2O vapour at given temperature [nTp]" }
};

//===============================================================

void TNode::databr_to_text_file( fstream& ff, bool with_comments, bool brief_mode, const char* path )
{
// fstream ff("DataBR.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");
  bool _comment = with_comments;

  TPrintArrays  prar(f_amrl+1/*54*/, DataBR_fields, ff);

   if( _comment )
   {
      ff << "# " << _GEMIPM_version_stamp << endl << "# File: " << path << endl;
      ff << "# Comments can be marked with # $ ; as the first character in the line" << endl;
      ff << "# DBR text input file for node system recipe and speciation data" << endl;
      ff << "# (should be read only after the DCH and the IPM files)" << endl << endl;
      ff << "# (1): Flags controlling GEM IPM-3 operation and data exchange";
   }

#ifndef NODEARRAYLEVEL
   CNode->NodeStatusFMT = No_nodearray;
#endif
   prar.writeField(f_NodeHandle, CNode->NodeHandle, _comment, brief_mode  );
   prar.writeField(f_NodeTypeHY, CNode->NodeTypeHY, _comment, brief_mode  );
   prar.writeField(f_NodeTypeMT, CNode->NodeTypeMT, _comment, brief_mode  );
   prar.writeField(f_NodeStatusFMT, CNode->NodeStatusFMT, _comment, brief_mode  );
   prar.writeField(f_NodeStatusCH, CNode->NodeStatusCH, _comment, brief_mode  );
   prar.writeField(f_IterDone, CNode->IterDone, _comment, brief_mode  );

  if( _comment )
      ff << "\n\n## (2) Chemical scalar properies of the node system";

  prar.writeField(f_TK, CNode->TK, _comment, brief_mode  );
  prar.writeField(f_P, CNode->P, _comment, brief_mode  );

  prar.writeField(f_Vs, CNode->Vs, _comment, brief_mode  );
  prar.writeField(f_Vi, CNode->Vi, _comment, brief_mode  );

  prar.writeField(f_Ms, CNode->Ms, _comment, brief_mode  );
  prar.writeField(f_Mi, CNode->Mi, _comment, brief_mode  );


  prar.writeField(f_Hs, CNode->Hs, _comment, brief_mode  );
  prar.writeField(f_Hi, CNode->Hi, _comment, brief_mode  );

  prar.writeField(f_Gs, CNode->Gs, _comment, brief_mode  );

 if( CSD->ccPH[0] == PH_AQUEL )
 {
     prar.writeField(f_IS, CNode->IC, _comment, brief_mode  );
     prar.writeField(f_pH, CNode->pH, _comment, brief_mode  );
     prar.writeField(f_pe, CNode->pe, _comment, brief_mode  );
     prar.writeField(f_Eh, CNode->Eh, _comment, brief_mode  );
 }

  prar.writeField(f_Tm, CNode->Tm, _comment, brief_mode  );
  prar.writeField(f_dt, CNode->dt, _comment, brief_mode  );

#ifdef NODEARRAYLEVEL
  if( CNode->NodeStatusFMT != No_nodearray /*TNodeArray::na->nNodes() > 1*/ )
  {
   if( _comment )
      ff << "\n\n## (3) Scalar mass-trasport properties (used only at NodeArray level)";
   prar.writeField(f_Dif, CNode->Dif, _comment, brief_mode  );
   prar.writeField(f_Vt, CNode->Vt, _comment, brief_mode  );
   prar.writeField(f_vp, CNode->vp, _comment, brief_mode  );
   prar.writeField(f_eps, CNode->eps, _comment, brief_mode  );
   prar.writeField(f_Km, CNode->Km, _comment, brief_mode  );
   prar.writeField(f_Kf, CNode->Kf, _comment, brief_mode  );
   prar.writeField(f_S, CNode->S, _comment, brief_mode  );
   prar.writeField(f_Tr, CNode->Tr, _comment, brief_mode  );
   prar.writeField(f_h, CNode->h, _comment, brief_mode  );
   prar.writeField(f_rho, CNode->rho, _comment, brief_mode  );
   prar.writeField(f_al, CNode->al, _comment, brief_mode  );
   prar.writeField(f_at, CNode->at, _comment, brief_mode  );
   prar.writeField(f_av, CNode->av, _comment, brief_mode  );
   prar.writeField(f_hDl, CNode->hDl, _comment, brief_mode  );
   prar.writeField(f_hDt, CNode->hDt, _comment, brief_mode  );
   prar.writeField(f_hDv, CNode->hDv, _comment, brief_mode  );
   prar.writeField(f_nto, CNode->nto, _comment, brief_mode  );
  }
#endif

  if( _comment )
   {   ff << "\n\n### Arrays: for dimensions and index lists, see Section (2) of DCH file" << endl << endl;
       ff << "## (4) Data for Independent Components";
       prar.writeArray(  NULL, CSD->ICNL[0], CSD->nIC, MaxICN );
       //ff << endl;
   }

  prar.writeArray(  f_bIC,  CNode->bIC, CSD->nICb, -1L,_comment, brief_mode );
  prar.writeArray(  f_rMB,  CNode->rMB, CSD->nICb, -1L,_comment, brief_mode );
  prar.writeArray(  f_uIC,  CNode->uIC, CSD->nICb, -1L,_comment, brief_mode );
  prar.writeArray(  f_bSP,  CNode->bSP, CSD->nICb, -1L,_comment, brief_mode );

  if( _comment )
  {    ff << "\n\n## (5) Data for Dependent Components";
       prar.writeArray(  NULL, CSD->DCNL[0], CSD->nDC, MaxDCN );
       //ff << endl;
  }

  prar.writeArray(  f_xDC,  CNode->xDC, CSD->nDCb, -1L,_comment, brief_mode  );
  prar.writeArray(  f_gam,  CNode->gam, CSD->nDCb, -1L,_comment, brief_mode  );
  prar.writeArray(  f_dll,  CNode->dll, CSD->nDCb, -1L,_comment, brief_mode  );
  prar.writeArray(  f_dul,  CNode->dul, CSD->nDCb, -1L,_comment, brief_mode  );

  if( _comment )
  {    ff << "\n\n## (6) Data for Phases";
        prar.writeArray(  NULL, CSD->PHNL[0], CSD->nPH, MaxPHN );
       // ff << endl;
  }

  prar.writeArray(  f_aPH,  CNode->aPH, CSD->nPHb, -1L,_comment, brief_mode );
  prar.writeArray(  f_xPH,  CNode->xPH, CSD->nPHb, -1L,_comment, brief_mode );
  prar.writeArray(  f_vPS,  CNode->vPS, CSD->nPSb, -1L,_comment, brief_mode );
  prar.writeArray(  f_mPS,  CNode->mPS, CSD->nPSb, -1L,_comment, brief_mode );
  prar.writeArray(  f_xPA,  CNode->xPA, CSD->nPSb, -1L,_comment, brief_mode );
  prar.writeArray(  f_amru,  CNode->amru, CSD->nPSb, -1L,_comment, brief_mode );
  prar.writeArray(  f_amrl,  CNode->amrl, CSD->nPSb, -1L,_comment, brief_mode );

  if(!brief_mode || prar.getAlws( f_bPS ))
  {  if( _comment )
     {
          ff << DataBR_fields[f_bPS].comment.c_str();
	  prar.writeArray(  NULL, CSD->ICNL[0], CSD->nIC, MaxICN );
      }
     prar.writeArray(  f_bPS,  CNode->bPS, CSD->nPSb*CSD->nICb, CSD->nICb,false, brief_mode );
  }

  ff << endl;
  if( _comment )
   {   //  ff << "\n# reserved" << endl;
         ff << "\n# End of file"<< endl;
   }
}

// Reading work dataBR structure from text file
void TNode::databr_from_text_file( fstream& ff )
{
#ifndef NODEARRAYLEVEL
   double tmpVal;
#endif
   // fstream ff("DataBR.out", ios::out );
   // ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

 // mem_set( &CNode->Tm, 0, 19*sizeof(double));
 databr_reset( CNode );
 TReadArrays  rdar(f_amrl+1/*54*/, DataBR_fields, ff);
 long int nfild = rdar.findNext();
 while( nfild >=0 )
 {
   switch( nfild )
   {
    case f_NodeHandle: rdar.readArray( "NodeHandle",  &CNode->NodeHandle, 1);
            break;
    case f_NodeTypeHY: rdar.readArray( "NodeTypeHY",  &CNode->NodeTypeHY, 1);
            break;
    case f_NodeTypeMT: rdar.readArray( "NodeTypeMT",  &CNode->NodeTypeMT, 1);
            break;
    case f_NodeStatusFMT: rdar.readArray( "NodeStatusFMT",  &CNode->NodeStatusFMT, 1);
            break;
    case f_NodeStatusCH: rdar.readArray( "NodeStatusCH",  &CNode->NodeStatusCH, 1);
            break;
    case f_IterDone: rdar.readArray( "IterDone",  &CNode->IterDone, 1);
            break;
    case f_TK: rdar.readArray( "TK",  &CNode->TK, 1);
            break;
    case f_P: rdar.readArray( "P",  &CNode->P, 1);
            break;
    case f_Vs: rdar.readArray( "Vs", &CNode->Vs, 1);
            break;
    case f_Vi: rdar.readArray( "Vi",  &CNode->Vi, 1);
            break;
    case f_Ms: rdar.readArray( "Ms",  &CNode->Ms, 1);
            break;
    case f_Mi: rdar.readArray( "Mi",  &CNode->Mi, 1);
            break;
    case f_Hs: rdar.readArray( "Hs",  &CNode->Hs, 1);
            break;
    case f_Hi: rdar.readArray( "Hi",  &CNode->Hi, 1);
            break;
    case f_Gs: rdar.readArray( "Gs",  &CNode->Gs, 1);
             break;
    case f_IS: rdar.readArray( "IS",  &CNode->IC, 1);
            break;
    case f_pH: rdar.readArray( "pH",  &CNode->pH, 1);
            break;
    case f_pe: rdar.readArray( "pe",  &CNode->pe, 1);
            break;
    case f_Eh: rdar.readArray( "Eh",  &CNode->Eh, 1);
            break;
    case f_Tm: rdar.readArray( "Tm",  &CNode->Tm, 1);
            break;
    case f_dt: rdar.readArray( "dt",  &CNode->dt, 1);
            break;
#ifdef NODEARRAYLEVEL
   case f_Dif: rdar.readArray( "Dif",  &CNode->Dif, 1);
            break;
    case f_Vt: rdar.readArray( "Vt",  &CNode->Vt, 1);
            break;
    case f_vp: rdar.readArray( "vp",  &CNode->vp, 1);
            break;
    case f_eps: rdar.readArray( "eps",  &CNode->eps, 1);
            break;
    case f_Km: rdar.readArray( "Km",  &CNode->Km, 1);
            break;
    case f_Kf: rdar.readArray( "Kf",  &CNode->Kf, 1);
            break;
    case f_S: rdar.readArray( "S",  &CNode->S, 1);
            break;
    case f_Tr: rdar.readArray( "Tr",  &CNode->Tr, 1);
            break;
    case f_h: rdar.readArray( "h",  &CNode->h, 1);
            break;
    case f_rho: rdar.readArray( "rho",  &CNode->rho, 1);
            break;
    case f_al: rdar.readArray( "al",  &CNode->al, 1);
            break;
    case f_at: rdar.readArray( "at",  &CNode->at, 1);
            break;
    case f_av: rdar.readArray( "av",  &CNode->av, 1);
            break;
    case f_hDl: rdar.readArray( "hDl",  &CNode->hDl, 1);
            break;
    case f_hDt: rdar.readArray( "hDt",  &CNode->hDt, 1);
            break;
    case f_hDv: rdar.readArray( "hDv",  &CNode->hDv, 1);
            break;
    case f_nto: rdar.readArray( "nto",  &CNode->nto, 1);
            break;
#else
   case f_Dif: rdar.readArray( "Dif",  &tmpVal, 1);
            break;
    case f_Vt: rdar.readArray( "Vt",  &tmpVal, 1);
            break;
    case f_vp: rdar.readArray( "vp",  &tmpVal, 1);
            break;
    case f_eps: rdar.readArray( "eps",  &tmpVal, 1);
            break;
    case f_Km: rdar.readArray( "Km",  &tmpVal, 1);
            break;
    case f_Kf: rdar.readArray( "Kf",  &tmpVal, 1);
            break;
    case f_S: rdar.readArray( "S",  &tmpVal, 1);
            break;
    case f_Tr: rdar.readArray( "Tr",  &tmpVal, 1);
            break;
    case f_h: rdar.readArray( "h",  &tmpVal, 1);
            break;
    case f_rho: rdar.readArray( "rho",  &tmpVal, 1);
            break;
    case f_al: rdar.readArray( "al",  &tmpVal, 1);
            break;
    case f_at: rdar.readArray( "at",  &tmpVal, 1);
            break;
    case f_av: rdar.readArray( "av",  &tmpVal, 1);
            break;
    case f_hDl: rdar.readArray( "hDl",  &tmpVal, 1);
            break;
    case f_hDt: rdar.readArray( "hDt",  &tmpVal, 1);
            break;
    case f_hDv: rdar.readArray( "hDv",  &tmpVal, 1);
            break;
    case f_nto: rdar.readArray( "nto",  &tmpVal, 1);
            break;
#endif
    case f_bIC: rdar.readArray( "bIC",  CNode->bIC, CSD->nICb );
            break;
    case f_rMB: rdar.readArray( "rMB",  CNode->rMB, CSD->nICb );
            break;
    case f_uIC: rdar.readArray( "uIC",  CNode->uIC, CSD->nICb );
            break;
    case f_xDC: rdar.readArray( "xDC",  CNode->xDC, CSD->nDCb );
            break;
    case f_gam: rdar.readArray( "gam",  CNode->gam, CSD->nDCb );
            break;
    case f_dll: rdar.readArray( "dll",  CNode->dll, CSD->nDCb );
            break;
    case f_dul: rdar.readArray( "dul",  CNode->dul, CSD->nDCb );
            break;
    case f_aPH: rdar.readArray( "aPH",  CNode->aPH, CSD->nPHb );
            break;
    case f_xPH: rdar.readArray( "xPH",  CNode->xPH, CSD->nPHb );
            break;
    case f_vPS: rdar.readArray( "vPS",  CNode->vPS, CSD->nPSb );
            break;
    case f_mPS: rdar.readArray( "mPS",  CNode->mPS, CSD->nPSb );
            break;
    case f_bPS: rdar.readArray( "bPS",  CNode->bPS, CSD->nPSb*CSD->nICb );
            break;
    case f_xPA: rdar.readArray( "xPA",  CNode->xPA, CSD->nPSb );
            break;
    case f_bSP: rdar.readArray( "bSP",  CNode->bSP, CSD->nICb );
           break;
   case f_amru: rdar.readArray( "amru",  CNode->amru, CSD->nPSb );
           break;
   case f_amrl: rdar.readArray( "arml",  CNode->amrl, CSD->nPSb );
           break;
   }
   nfild = rdar.findNext();
 }

 // testing read
 gstring ret = rdar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from DataBR structure";
    Error( "Error", ret);
  }
}

//==============================================================================

void TNode::datach_to_text_file( fstream& ff, bool with_comments, bool brief_mode, const char* path )
{
// fstream ff("DataCH.out", ios::out );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");
  bool _comment = with_comments;
  TPrintArrays  prar1(14, DataCH_static_fields, ff);
  TPrintArrays  prar(30, DataCH_dynamic_fields, ff);
  if( CSD->nIC == CSD->nICb )
          prar.setNoAlws( f_xic);
  if(CSD->nDC == CSD->nDCb )
          prar.setNoAlws( f_xdc);
  if(CSD->nPH == CSD->nPHb )
          prar.setNoAlws( f_xph );

  if( _comment )
  {
     ff << "# " << _GEMIPM_version_stamp << endl << "# File: " << path << endl;
     ff << "# Comments can be marked with # $ ; as the first character in the line" << endl;
     ff << "# DCH text input file (should be read before IPM and DBR files)" << endl << endl;
     ff << "## (1) Dimensions for memory allocation";
  }
  prar1.writeField(f_nIC, CSD->nIC, _comment, brief_mode  );
  prar1.writeField(f_nDC, CSD->nDC, _comment, brief_mode  );
  prar1.writeField(f_nPH, CSD->nPH, _comment, brief_mode  );
  prar1.writeField(f_nPS, CSD->nPS, _comment, brief_mode  );
  prar1.writeField(f_nDCs, CSD->nDCs, _comment, brief_mode  );

  if( _comment )
    ff << endl << "\n## (2) Dimensions for DBR node recipe (memory allocation)";
  prar1.writeField(f_nICb, CSD->nICb, _comment, brief_mode  );
  prar1.writeField(f_nDCb, CSD->nDCb, _comment, brief_mode  );
  prar1.writeField(f_nPHb, CSD->nPHb, _comment, brief_mode  );
  prar1.writeField(f_nPSb, CSD->nPSb, _comment, brief_mode  );

  if( _comment )
    ff << endl << "\n## (3) Dimensions for thermodynamic data arrays";
  prar1.writeField(f_nTp, CSD->nTp, _comment, brief_mode  );
  prar1.writeField(f_nPp, CSD->nPp, _comment, brief_mode  );
  prar1.writeField(f_iGrd, CSD->iGrd, _comment, brief_mode  );
  prar1.writeField(f_fAalp, CSD->nAalp, _comment, brief_mode  );
  prar1.writeField(f_mLook, CSD->mLook, _comment, brief_mode  );

  ff << endl << "\n<END_DIM>" << endl;

// dynamic arrays - must follow static data
  if( _comment )
     ff << "\n## (4) DBR node recipe connection index lists";
  prar.writeArray(  f_xic, CSD->xic, CSD->nICb, -1L,_comment, brief_mode);
  prar.writeArray(  f_xdc, CSD->xdc, CSD->nDCb, -1L,_comment, brief_mode);
  prar.writeArray(  f_xph, CSD->xph, CSD->nPHb, -1L,_comment, brief_mode);

  if( _comment )
     ff << "\n\n## (5) Independent Components and their properties";
  if(!brief_mode || prar.getAlws( f_ICNL ))
  {
     if( _comment )
         ff << "\n# ICNL: List of Independent Component names (<=4 characters per name) [nIC]";
      prar.writeArray(  "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
  }
  prar.writeArrayF(  f_ccIC, CSD->ccIC, CSD->nIC, 1L,_comment, brief_mode );
  prar.writeArray(  f_ICmm, CSD->ICmm, CSD->nIC, -1L,_comment, brief_mode);

  if( _comment )
    ff << "\n\n## (6) Dependent Components and their codes";
  if(!brief_mode || prar.getAlws( f_DCNL ))
  {	  if( _comment )
       ff << "\n# DCNL: Name list of Dependent Components (<=16 characters per name) [nDC]";
     prar.writeArray(  "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
  }
  prar.writeArrayF(  f_ccDC, CSD->ccDC, CSD->nDC, 1L,_comment, brief_mode );
  prar.writeArray(  f_DCmm, CSD->DCmm, CSD->nDC, -1L,_comment, brief_mode);

  if( _comment )
    ff << "\n\n## (7) Phases and their codes" << endl;
  if(!brief_mode || prar.getAlws( f_PHNL ))
  { if( _comment )
      ff << "# PHNL: List of Phase names (<=16 characters per name) [nPH]";
    prar.writeArray(  "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
  }
  prar.writeArrayF(  f_ccPH, CSD->ccPH, CSD->nPH, 1L,_comment, brief_mode );
  prar.writeArray(  f_nDCinPH, CSD->nDCinPH, CSD->nPH, -1L,_comment, brief_mode);

  if( _comment )
    ff << "\n\n# (8) Data for Dependent Components";
  prar.writeArray(  f_A, CSD->A, CSD->nDC*CSD->nIC, CSD->nIC, _comment, brief_mode );
  ff << endl;

  if( _comment )
    ff << "\n## (9) Thermodynamic data for Dependent Components";
  prar.writeField(  f_Ttol, CSD->Ttol, _comment, brief_mode  );
  prar.writeArray(  f_TKval, CSD->TKval, CSD->nTp, -1L,_comment, brief_mode );
  prar.writeArray(  f_Psat, CSD->Psat, CSD->nTp, -1L,_comment, brief_mode );
  ff << endl;

  prar.writeField(  f_Ptol, CSD->Ptol, _comment, brief_mode  );
  prar.writeArray(  f_Pval, CSD->Pval, CSD->nPp,  -1L,_comment, brief_mode );

  if( CSD->ccPH[0] == PH_AQUEL )
  {
    prar.writeArray(  f_denW, CSD->denW, 5*(gridTP()), gridTP(), _comment, brief_mode );
    prar.writeArray(  f_denWg, CSD->denWg, 5*(gridTP()), gridTP(), _comment, brief_mode );
    prar.writeArray(  f_epsW, CSD->epsW,  5*(gridTP()), gridTP(), _comment, brief_mode );
    prar.writeArray(  f_epsWg, CSD->epsWg, 5*(gridTP()),  gridTP(), _comment, brief_mode );
  }
  prar.writeArray(  f_V0, CSD->V0,  CSD->nDC*gridTP(), gridTP(), _comment, brief_mode );
  prar.writeArray(  f_G0, CSD->G0, CSD->nDC*gridTP(), gridTP(), _comment, brief_mode );
  prar.writeArray(  f_H0, CSD->H0,  CSD->nDC*gridTP(),gridTP(), _comment, brief_mode );
  prar.writeArray(  f_S0, CSD->S0,CSD->nDC*gridTP(),  gridTP(), _comment, brief_mode  );
  prar.writeArray(  f_Cp0, CSD->Cp0,CSD->nDC*gridTP(), gridTP(), _comment, brief_mode  );
  prar.writeArray(  f_A0, CSD->A0, CSD->nDC*gridTP(), gridTP(), _comment, brief_mode  );
  prar.writeArray(  f_U0, CSD->U0, CSD->nDC*gridTP(), gridTP(), _comment, brief_mode  );

  if( CSD->iGrd  )
    prar.writeArray(  f_DD, CSD->DD, CSD->nDCs*gridTP(),  gridTP(),  _comment, brief_mode);

  ff << endl;
  if( _comment )
      ff << "\n# End of file";
}

// Reading dataCH structure from text file
void TNode::datach_from_text_file(fstream& ff)
{
  long int ii;
// fstream ff("DataCH.out", ios::in );
// ErrorIf( !ff.good() , "DataCH.out", "Fileopen error");

// static arrays
 TReadArrays  rdar( 14, DataCH_static_fields, ff);
 long int nfild = rdar.findNext();
 while( nfild >=0 )
 {
   switch( nfild )
   {
    case f_nIC: rdar.readArray( "nIC", &CSD->nIC, 1);
            break;
    case f_nDC: rdar.readArray( "nDC", &CSD->nDC, 1);
            break;
    case f_nPH: rdar.readArray( "nPH", &CSD->nPH, 1);
            break;
    case f_nPS: rdar.readArray( "nPS", &CSD->nPS, 1);
            break;
    case f_nDCs: rdar.readArray( "nDCs", &CSD->nDCs, 1);
            break;
    case f_nICb: rdar.readArray( "nICb", &CSD->nICb, 1);
            break;
    case f_nDCb: rdar.readArray( "nDCb", &CSD->nDCb, 1);
            break;
    case f_nPHb: rdar.readArray( "nPHb", &CSD->nPHb, 1);
            break;
    case f_nPSb: rdar.readArray( "nPSb", &CSD->nPSb, 1);
            break;
    case f_nTp: rdar.readArray( "nTp", &CSD->nTp, 1);
            break;
    case f_nPp: rdar.readArray( "nPp", &CSD->nPp, 1);
            break;
    case f_iGrd: rdar.readArray( "iGrd", &CSD->iGrd, 1);
            break;
    case f_fAalp: rdar.readArray( "fAalp", &CSD->nAalp, 1);
             break;
    case f_mLook: rdar.readArray( "mLook", &CSD->mLook, 1);
             break;
  }
   nfild = rdar.findNext();
 }

 // testing read
 gstring ret = rdar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from DataCH structure";
    Error( "Error", ret);
  }

  datach_realloc();
  databr_realloc();

//dynamic data
 TReadArrays  rddar( 30, DataCH_dynamic_fields, ff);

   if( CSD->iGrd  )
      rddar.setAlws( f_DD /*28 "DD"*/);

  // default set up
  for( ii=0; ii< CSD->nDCs*gridTP(); ii++ )
  {
    if( CSD->DD ) CSD->DD[ii] = 0.;
    if( CSD->Cp0) CSD->Cp0[ii] = 0.;
    if( CSD->H0 ) CSD->H0[ii] = 0.;
    if( CSD->S0 ) CSD->S0[ii] = 0.;
    if( CSD->A0 ) CSD->A0[ii] = 0.;
    if( CSD->U0 ) CSD->U0[ii] = 0.;
  }
  CSD->Ttol = 0.1;
  CSD->Ptol = 10000;
  if( CSD->nIC == CSD->nICb )
  {
    rddar.setNoAlws( f_xic /*"xic"*/ );
    for( ii=0; ii< CSD->nICb; ii++ )
      CSD->xic[ii] = ii;
  }
  if(CSD->nDC == CSD->nDCb )
  {
    rddar.setNoAlws( f_xdc /*"xdc"*/);
    for( ii=0; ii< CSD->nDCb; ii++ )
      CSD->xdc[ii] = ii;
  }
  if(CSD->nPH == CSD->nPHb )
  {
    rddar.setNoAlws( f_xph /*"xph"*/);
    for( ii=0; ii< CSD->nPHb; ii++ )
      CSD->xph[ii] = ii;
  }

  nfild = rddar.findNext();
  while( nfild >=0 )
  {
   switch( nfild )
   {
    case f_xic: rddar.readArray( "xic", CSD->xic, CSD->nICb);
            break;
    case f_xdc: rddar.readArray( "xdc", CSD->xdc, CSD->nDCb);
            break;
    case f_xph: rddar.readArray( "xph", CSD->xph, CSD->nPHb);
            break;
    case f_ICNL: rddar.readArray( "ICNL", CSD->ICNL[0], CSD->nIC, MaxICN );
            break;
    case f_ccIC: rddar.readArray( "ccIC", CSD->ccIC, CSD->nIC, 1 );
            break;
    case f_ICmm: rddar.readArray( "ICmm", CSD->ICmm, CSD->nIC);
            break;
    case f_DCNL: rddar.readArray( "DCNL", CSD->DCNL[0], CSD->nDC, MaxDCN );
            break;
    case f_ccDC: rddar.readArray( "ccDC", CSD->ccDC, CSD->nDC, 1 );
            break;
    case f_DCmm: rddar.readArray( "DCmm", CSD->DCmm, CSD->nDC);
            break;
    case f_PHNL: rddar.readArray( "PHNL", CSD->PHNL[0], CSD->nPH, MaxPHN );
            break;
    case f_ccPH: rddar.readArray( "ccPH", CSD->ccPH, CSD->nPH, 1 );
            break;
    case f_nDCinPH: rddar.readArray( "nDCinPH", CSD->nDCinPH, CSD->nPH);
            break;
    case f_A: rddar.readArray( "A", CSD->A, CSD->nDC*CSD->nIC );
            break;
    case f_Ttol: rddar.readArray( "Ttol", &CSD->Ttol, 1);
            break;
    case f_TKval: rddar.readArray( "TKval", CSD->TKval, CSD->nTp );
            break;
    case f_Ptol: rddar.readArray( "Ptol", &CSD->Ptol, 1);
            break;
    case f_Pval: rddar.readArray( "Pval", CSD->Pval, CSD->nPp );
              break;
    case f_denW: if( !CSD->denW )
                   Error( "Error", "Array denW is not allocated in DCH!");
             rddar.readArray( "denW", CSD->denW, 5*gridTP() );
              break;
    case f_denWg: if( !CSD->denWg )
                   Error( "Error", "Array denWg is not allocated in DCH!");
             rddar.readArray( "denWg", CSD->denWg, 5*gridTP() );
              break;
    case f_epsW: if( !CSD->epsW )
                   Error( "Error", "Array epsW is not allocated in DCH!");
             rddar.readArray( "epsW", CSD->epsW,  5*gridTP() );
            break;
    case f_epsWg: if( !CSD->epsWg )
                   Error( "Error", "Array epsWg is not allocated in DCH!");
             rddar.readArray( "epsWg", CSD->epsWg,  5*gridTP() );
            break;
    case f_V0: rddar.readArray( "V0", CSD->V0,  CSD->nDC*gridTP() );
            break;
    case f_G0: rddar.readArray( "G0", CSD->G0, CSD->nDC*gridTP() );
              break;
    case f_H0: if( !CSD->H0 )
                   Error( "Error", "Array HO is not allocated in DCH!");
            rddar.readArray( "H0", CSD->H0,  CSD->nDC*gridTP());
            break;
    case f_S0: if( !CSD->S0 )
                   Error( "Error", "Array S0 is not allocated in DCH!");
            rddar.readArray( "S0", CSD->S0,CSD->nDC*gridTP());
            break;
    case f_Cp0: if( !CSD->Cp0 )
                   Error( "Error", "Array CpO is not allocated in DCH!");
            rddar.readArray( "Cp0", CSD->Cp0,CSD->nDC*gridTP() );
            break;
    case f_A0: if( !CSD->A0 )
                   Error( "Error", "Array AO is not allocated in DCH!");
            rddar.readArray( "A0", CSD->A0, CSD->nDC*gridTP() );
            break;
    case f_U0: if( !CSD->U0 )
                   Error( "Error", "Array UO is not allocated in DCH!");
            rddar.readArray( "U0", CSD->U0, CSD->nDC*gridTP() );
            break;
    case f_DD: if( !CSD->DD )
                    Error( "Error", "Array DD is not allocated in DCH!");
            rddar.readArray( "DD", CSD->DD, CSD->nDCs*gridTP());
           break;
    case f_Psat: rddar.readArray( "Psat", CSD->Psat, CSD->nTp );
           break;
  }
     nfild = rddar.findNext();
 }

  // Set up flags
  if( CSD->ccPH[0] != PH_AQUEL )
  {
        rddar.setNoAlws( f_denW );
        rddar.setNoAlws( f_denWg );
        rddar.setNoAlws( f_epsW );
        rddar.setNoAlws( f_epsWg );
  }

 // testing read
 ret = rddar.testRead();
 if( !ret.empty() )
  { ret += " - fields must be read from DataCH structure";
    Error( "Error", ret);
  }

}

//---------------------------------------------------------------
// new i/o structures

// Writing DataCH to binary file
void TNode::datach_to_file( GemDataStream& ff )
{
// const data
   ff.writeArray( &CSD->nIC, 14 );
   ff.writeArray( &CSD->Ttol, 4 );

//dynamic data
   ff.writeArray( CSD->nDCinPH, CSD->nPH );
//   if( CSD->nICb >0 )
   ff.writeArray( CSD->xic, CSD->nICb );
   ff.writeArray( CSD->xdc, CSD->nDCb );
   ff.writeArray( CSD->xph, CSD->nPHb );

   ff.writeArray( CSD->A, CSD->nIC*CSD->nDC );
   ff.writeArray( CSD->ICmm, CSD->nIC );
   ff.writeArray( CSD->DCmm, CSD->nDC );

   ff.writeArray( CSD->TKval,  CSD->nTp );
   ff.writeArray( CSD->Psat,  CSD->nTp );
   ff.writeArray( CSD->Pval,  CSD->nPp );

   ff.writeArray( CSD->ccIC, CSD->nIC );
   ff.writeArray( CSD->ccDC, CSD->nDC );
   ff.writeArray( CSD->ccPH, CSD->nPH );

   if( CSD->ccPH[0] == PH_AQUEL )
   { ff.writeArray( CSD->denW,  5*gridTP() );
     ff.writeArray( CSD->denWg,  5*gridTP() );
     ff.writeArray( CSD->epsW, 5*gridTP() );
     ff.writeArray( CSD->epsWg, 5*gridTP() );
   }
   ff.writeArray( CSD->G0,  CSD->nDC*gridTP() );
   ff.writeArray( CSD->V0,  CSD->nDC*gridTP() );
   ff.writeArray( CSD->H0,  CSD->nDC*gridTP() );
   ff.writeArray( CSD->S0, CSD->nDC*gridTP() );
   ff.writeArray( CSD->Cp0, CSD->nDC*gridTP() );
   ff.writeArray( CSD->A0, CSD->nDC*gridTP() );
   ff.writeArray( CSD->U0, CSD->nDC*gridTP() );
   if(  CSD->iGrd  )
      ff.writeArray( CSD->DD, CSD->nDCs*gridTP() );

   ff.writeArray( (char *)CSD->ICNL, MaxICN*CSD->nIC );
   ff.writeArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC );
   ff.writeArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH );
}

// Reading DataCH structure from binary file
void TNode::datach_from_file( GemDataStream& ff )
{
// const data
   ff.readArray( &CSD->nIC, 14 );
   ff.readArray( &CSD->Ttol, 4 );

  datach_realloc();
  databr_realloc();

//dynamic data
   ff.readArray( CSD->nDCinPH, CSD->nPH );
//   if( CSD->nICb >0 )
   ff.readArray( CSD->xic, CSD->nICb );
   ff.readArray( CSD->xdc, CSD->nDCb );
   ff.readArray( CSD->xph, CSD->nPHb );

   ff.readArray( CSD->A, CSD->nIC*CSD->nDC );
   ff.readArray( CSD->ICmm, CSD->nIC );
   ff.readArray( CSD->DCmm, CSD->nDC );

   ff.readArray( CSD->TKval,  CSD->nTp );
   ff.readArray( CSD->Psat,  CSD->nTp );
   ff.readArray( CSD->Pval,  CSD->nPp );

   ff.readArray( CSD->ccIC, CSD->nIC );
   ff.readArray( CSD->ccDC, CSD->nDC );
   ff.readArray( CSD->ccPH, CSD->nPH );

   if( CSD->ccPH[0] == PH_AQUEL )
   {
         ff.readArray( CSD->denW,  5*gridTP() );
         ff.readArray( CSD->denWg,  5*gridTP() );
     ff.readArray( CSD->epsW, 5*gridTP() );
     ff.readArray( CSD->epsWg, 5*gridTP() );
   }
   ff.readArray( CSD->G0,  CSD->nDC*gridTP() );
   ff.readArray( CSD->V0,  CSD->nDC*gridTP() );
     ff.readArray( CSD->H0,  CSD->nDC*gridTP() );
     ff.readArray( CSD->S0, CSD->nDC*gridTP() );
     ff.readArray( CSD->Cp0, CSD->nDC*gridTP() );
     ff.readArray( CSD->A0, CSD->nDC*gridTP() );
     ff.readArray( CSD->U0, CSD->nDC*gridTP() );
   if(  CSD->iGrd  )
     ff.readArray( CSD->DD, CSD->nDCs*gridTP() );

   ff.readArray( (char *)CSD->ICNL, MaxICN*CSD->nIC );
   ff.readArray( (char *)CSD->DCNL, MaxDCN*CSD->nDC );
   ff.readArray( (char *)CSD->PHNL, MaxPHN*CSD->nPH );


}

// allocating DataCH structure
void TNode::datach_realloc()
{
  if( CSD->mLook == 1 &&  (CSD->nPp != CSD->nTp) )
     Error( "No-interpolation mode",
           "Different number of points for temperature and pressure ");

 CSD->nDCinPH = new long int[CSD->nPH];

 if( CSD->nICb >0 )
   CSD->xic = new long int[CSD->nICb];
 else  CSD->xic = 0;
 if( CSD->nDCb >0 )
   CSD->xdc = new long int[CSD->nDCb];
 else  CSD->xdc = 0;
 if( CSD->nPHb >0 )
   CSD->xph = new long int[CSD->nPHb];
 else  CSD->xph = 0;

  CSD->A = new double[CSD->nIC*CSD->nDC];
  CSD->ICmm = new double[CSD->nIC];
  CSD->DCmm = new double[CSD->nDC];
CSD->DCmm[0] = 0.0;   // Added by DK on 03.03.2007

  CSD->TKval = new double[CSD->nTp];
  CSD->Psat = new double[CSD->nTp];
  CSD->Pval = new double[CSD->nPp];

  CSD->denW = new double[ 5*gridTP()];
  CSD->denWg = new double[ 5*gridTP()];
  CSD->epsW = new double[ 5*gridTP()];
  CSD->epsWg = new double[ 5*gridTP()];

  CSD->G0 = new double[CSD->nDC*gridTP()];
  CSD->V0 = new double[CSD->nDC*gridTP()];
  CSD->H0 = new double[CSD->nDC*gridTP()];
  CSD->S0 = new double[CSD->nDC*gridTP()];
  CSD->Cp0 = new double[CSD->nDC*gridTP()];
  CSD->A0 = new double[CSD->nDC*gridTP()];
  CSD->U0 = new double[CSD->nDC*gridTP()];

  if(  CSD->iGrd  )
       CSD->DD = new double[CSD->nDCs*gridTP()];
  else
       CSD->DD = 0;
  CSD->ICNL = new char[CSD->nIC][MaxICN];
  CSD->DCNL = new char[CSD->nDC][MaxDCN];
  CSD->PHNL = new char[CSD->nPH][MaxPHN];

  CSD->ccIC = new char[CSD->nIC];
  CSD->ccDC = new char[CSD->nDC];
  CSD->ccPH = new char[CSD->nPH];
}

// free dynamic memory
void TNode::datach_free()
{
 if( CSD->nDCinPH )
  { delete[] CSD->nDCinPH;
    CSD->nDCinPH = 0;
  }
 if( CSD->xic )
  { delete[] CSD->xic;
    CSD->xic = 0;
  }
 if( CSD->xdc )
  { delete[] CSD->xdc;
    CSD->xdc = 0;
  }
 if( CSD->xph )
  { delete[] CSD->xph;
    CSD->xph = 0;
  }
 if( CSD->A )
  { delete[] CSD->A;
    CSD->A = 0;
  }
 if( CSD->ICmm )
  { delete[] CSD->ICmm;
    CSD->ICmm = 0;
  }
 if( CSD->DCmm )
  { delete[] CSD->DCmm;
    CSD->DCmm = 0;
  }

 if( CSD->TKval )
  { delete[] CSD->TKval;
    CSD->TKval = 0;
  }
 if( CSD->Psat )
  { delete[] CSD->Psat;
    CSD->Psat = 0;
  }
 if( CSD->Pval )
  { delete[] CSD->Pval;
    CSD->Pval = 0;
  }

 if( CSD->denW )
  { delete[] CSD->denW;
    CSD->denW = 0;
  }
 if( CSD->denWg )
  { delete[] CSD->denWg;
    CSD->denWg = 0;
  }
 if( CSD->epsW )
  { delete[] CSD->epsW;
    CSD->epsW = 0;
  }
 if( CSD->epsWg )
  { delete[] CSD->epsWg;
    CSD->epsWg = 0;
  }
 if( CSD->G0 )
  { delete[] CSD->G0;
    CSD->G0 = 0;
  }
 if( CSD->V0 )
  { delete[] CSD->V0;
    CSD->V0 = 0;
  }
 if( CSD->H0 )
  { delete[] CSD->H0;
    CSD->H0 = 0;
  }
 if( CSD->Cp0 )
  { delete[] CSD->Cp0;
    CSD->Cp0 = 0;
  }
  if( CSD->S0 )
  { delete[] CSD->S0;
     CSD->S0 = 0;
  }
  if( CSD->A0 )
  { delete[] CSD->A0;
     CSD->A0 = 0;
  }
  if( CSD->U0 )
  { delete[] CSD->U0;
     CSD->U0 = 0;
  }
  if( CSD->DD )
  { delete[] CSD->DD;
     CSD->DD = 0;
  }

 if( CSD->ICNL )
  { delete[] CSD->ICNL;
    CSD->ICNL = 0;
  }
 if( CSD->DCNL )
  { delete[] CSD->DCNL;
    CSD->DCNL = 0;
  }
 if( CSD->PHNL )
  { delete[] CSD->PHNL;
    CSD->PHNL = 0;
  }

 if( CSD->ccIC )
  { delete[] CSD->ccIC;
    CSD->ccIC = 0;
  }
 if( CSD->ccDC )
  { delete[] CSD->ccDC;
    CSD->ccDC = 0;
  }
 if( CSD->ccPH )
  { delete[] CSD->ccPH;
    CSD->ccPH = 0;
  }
 // delete[] CSD;
}

// writing DataBR to binary file
void TNode::databr_to_file( GemDataStream& ff )
{

#ifndef NODEARRAYLEVEL
   CNode->NodeStatusFMT = No_nodearray;
#endif
// const data
   ff.writeArray( &CNode->NodeHandle, 6 );

#ifdef NODEARRAYLEVEL
   if( CNode->NodeStatusFMT != No_nodearray )
       ff.writeArray( &CNode->TK, 32 );
   else
      ff.writeArray( &CNode->TK, 15 );
#else
      ff.writeArray( &CNode->TK, 15 );
#endif

//dynamic data
   ff.writeArray( CNode->bIC, CSD->nICb );
   ff.writeArray( CNode->rMB, CSD->nICb );
   ff.writeArray( CNode->uIC, CSD->nICb );
   ff.writeArray( CNode->bSP, CSD->nICb );

   ff.writeArray( CNode->xDC, CSD->nDCb );
   ff.writeArray( CNode->gam, CSD->nDCb );
   ff.writeArray( CNode->dul, CSD->nDCb );
   ff.writeArray( CNode->dll, CSD->nDCb );

   if( CSD->nAalp >0 )
        ff.writeArray( CNode->aPH, CSD->nPHb );
   ff.writeArray( CNode->xPH, CSD->nPHb );
   ff.writeArray( CNode->vPS, CSD->nPSb );
   ff.writeArray( CNode->mPS, CSD->nPSb );
   ff.writeArray( CNode->bPS, CSD->nPSb*CSD->nICb );
   ff.writeArray( CNode->xPA, CSD->nPSb );
   ff.writeArray( CNode->amru, CSD->nPSb );
   ff.writeArray( CNode->amrl, CSD->nPSb );
//   datach_to_text_file();
//   databr_to_text_file();
}

// Reading work dataBR structure from binary file
void TNode::databr_from_file( GemDataStream& ff )
{
// const data
   ff.readArray( &CNode->NodeHandle, 6 );

#ifdef NODEARRAYLEVEL
   if( CNode->NodeStatusFMT != No_nodearray )
       ff.readArray( &CNode->TK, 32 );
   else
      ff.readArray( &CNode->TK, 15 );
#else
   ErrorIf(CNode->NodeStatusFMT != No_nodearray, ff.GetPath(),
     "Error reading work dataBR structure from binary file (No_nodearray)");
   ff.readArray( &CNode->TK, 15 );
#endif
//dynamic data
   ff.readArray( CNode->bIC, CSD->nICb );
   ff.readArray( CNode->rMB, CSD->nICb );
   ff.readArray( CNode->uIC, CSD->nICb );
   ff.readArray( CNode->bSP, CSD->nICb );

   ff.readArray( CNode->xDC, CSD->nDCb );
   ff.readArray( CNode->gam, CSD->nDCb );
   ff.readArray( CNode->dul, CSD->nDCb );
   ff.readArray( CNode->dll, CSD->nDCb );

   if( CSD->nAalp >0 )
        ff.readArray( CNode->aPH, CSD->nPHb );
   ff.readArray( CNode->xPH, CSD->nPHb );
   ff.readArray( CNode->vPS, CSD->nPSb );
   ff.readArray( CNode->mPS, CSD->nPSb );
   ff.readArray( CNode->bPS, CSD->nPSb*CSD->nICb );
   ff.readArray( CNode->xPA, CSD->nPSb );
   ff.readArray( CNode->amru, CSD->nPSb );
   ff.readArray( CNode->amrl, CSD->nPSb );
}

// Allocates DataBR structure
void TNode::databr_realloc()
{
  long int j,k;
  CNode->bIC = new double[CSD->nICb];
  CNode->rMB = new double[CSD->nICb];
  CNode->uIC = new double[CSD->nICb];
  CNode->bSP = new double[CSD->nICb];

  for(  j=0; j<CSD->nICb; j++ )
  {
      CNode->rMB[j] = 0.;
      CNode->uIC[j] = 0.;
      CNode->bSP[j] = 0.;
   }

  CNode->xDC = new double[CSD->nDCb];
  CNode->gam = new double[CSD->nDCb];

  for(  j=0; j<CSD->nDCb; j++ )
  {
    CNode->xDC[j] = 0.;
    CNode->gam[j] = 1.;
  }

  //  default assignment
 CNode->dul = new double[CSD->nDCb];
 for(  j=0; j<CSD->nDCb; j++ )
   CNode->dul[j] = 1.0e6;            // default assignment
 CNode->dll = new double[CSD->nDCb];
 for(  j=0; j<CSD->nDCb; j++ )
   CNode->dll[j] = 0.0;              // default assignment

 if( CSD->nAalp >0 )
 {
    CNode->aPH = new double[CSD->nPHb];
    for(  k=0; k<CSD->nPHb; k++ )
      CNode->aPH[k] = 0.0;       // default assignment
 }
 else
    CNode->aPH = 0;

 CNode->xPH = new double[CSD->nPHb];

 for(  k=0; k<CSD->nPHb; k++ )
   CNode->xPH[k] = 0.0;       // default assignment

 CNode->vPS = new double[CSD->nPSb];
 CNode->mPS = new double[CSD->nPSb];
 CNode->bPS = new double[CSD->nPSb*CSD->nICb];
 CNode->xPA = new double[CSD->nPSb];
 CNode->amru = new double[CSD->nPSb];
 CNode->amrl = new double[CSD->nPSb];

 for(  k=0; k<CSD->nPSb; k++ )
 {
     CNode->vPS[k] = 0.0;
     CNode->mPS[k] = 0.0;
     CNode->xPA[k] = 0.0;
     for(  j=0; j<CSD->nICb; j++ )
        CNode->bPS[k*CSD->nICb+j] = 0.0;
     CNode->amru[k] = 0.0;
     CNode->amrl[k] = 0.0;
 }
}

// free dynamic memory
DATABR * TNode::databr_free( DATABR *CNode_ )
{
  if( CNode_ == 0)
    CNode_ = CNode;

 if( CNode_->bIC )
 { delete[] CNode_->bIC;
   CNode_->bIC = 0;
 }
 if( CNode_->rMB )
 { delete[] CNode_->rMB;
   CNode_->rMB = 0;
 }
 if( CNode_->uIC )
 { delete[] CNode_->uIC;
   CNode_->uIC = 0;
 }

 if( CNode_->xDC )
  { delete[] CNode_->xDC;
    CNode_->xDC = 0;
  }
 if( CNode_->gam )
  { delete[] CNode_->gam;
    CNode_->gam = 0;
  }
 if( CNode_->dul )
   { delete[] CNode_->dul;
     CNode_->dul = 0;
   }
 if( CNode_->dll )
   { delete[] CNode_->dll;
     CNode_->dll = 0;
   }

 if( CNode_->aPH )
 { delete[] CNode_->aPH;
   CNode_->aPH = 0;
 }
 if( CNode_->xPH )
  { delete[] CNode_->xPH;
    CNode_->xPH = 0;
  }
 if( CNode_->vPS )
  { delete[] CNode_->vPS;
    CNode_->vPS = 0;
  }
 if( CNode_->mPS )
  { delete[] CNode_->mPS;
    CNode_->mPS = 0;
  }
 if( CNode_->bPS )
  { delete[] CNode_->bPS;
    CNode_->bPS = 0;
  }
 if( CNode_->xPA )
  { delete[] CNode_->xPA;
    CNode_->xPA = 0;
  }

 if( CNode_->bSP )
 { delete[] CNode_->bSP;
   CNode_->bSP = 0;
 }
 if( CNode_->amru )
  { delete[] CNode_->amru;
    CNode_->amru = 0;
  }
 if( CNode_->amrl )
  { delete[] CNode_->amrl;
    CNode_->amrl = 0;
  }
 delete CNode_;
 return NULL;
}

// set default values(zeros) for DATABR structure
void TNode::databr_reset( DATABR *CNode, long int level )
{
	//  FMT variables (units or dimensionsless) - to be used for storing them
	//  at the nodearray level = 0.; normally not used in the single-node FMT-GEM coupling
		CNode->Tm = 0.;
		CNode->dt = 0.;
#ifdef NODEARRAYLEVEL
		CNode->Dif = 0.;
		CNode->Vt = 0.;
		CNode->vp = 0.;
		CNode->eps = 0.;
		CNode->Km = 0.;
		CNode->Kf = 0.;
		CNode->S = 0.;
		CNode->Tr = 0.;
		CNode->h = 0.;
		CNode->rho = 0.;
		CNode->al = 0.;
		CNode->at = 0.;
		CNode->av = 0.;
		CNode->hDl = 0.;
		CNode->hDt = 0.;
		CNode->hDv = 0.;
        CNode->nto = 0.; //19
#endif
		if(level <1 )
          return;

   CNode->NodeHandle = 0;
   CNode->NodeTypeHY = normal;
   CNode->NodeTypeMT = normal;
   CNode->NodeStatusFMT = Initial_RUN;
   CNode->NodeStatusCH = NEED_GEM_AIA;
   CNode->IterDone = 0;      //6

// Chemical scalar variables
	CNode->TK = 0.;
	CNode->P = 0.;
	CNode->Vs = 0.;
	CNode->Vi = 0.;
	CNode->Ms = 0.;
	CNode->Mi = 0.;
	CNode->Gs = 0.;
	CNode->Hs = 0.;
	CNode->Hi = 0.;
	CNode->IC = 0.;
	CNode->pH = 0.;
	CNode->pe = 0.;
	CNode->Eh = 0.; //13

	if( level < 2 )
       return;

// Data arrays - dimensions nICb, nDCb, nPHb, nPSb see in the DATACH structure
	CNode->bIC = 0;
	CNode->rMB = 0;
	CNode->uIC = 0;
	CNode->xDC = 0;
	CNode->gam = 0;
   CNode->dul = 0;
   CNode->dll = 0;
   CNode->aPH = 0;
   CNode->xPH = 0;
   CNode->vPS = 0;
   CNode->mPS = 0;
   CNode->bPS = 0;
   CNode->xPA = 0;
   CNode->bSP = 0;
   CNode->amru = 0;
   CNode->amrl = 0;
}

// set default values(zeros) for DATACH structure
void TNode::datach_reset()
{
	CSD->nIC = 0;
	CSD->nDC = 0;
	CSD->nPH = 0;
	CSD->nPS = 0;
	CSD->nDCs = 0;
	CSD->nTp = 0;
	CSD->nPp = 0;
	CSD->iGrd = 0;
	CSD->nAalp = 0;
	CSD->nICb = 0;
	CSD->nDCb = 0;
	CSD->nPHb = 0;
	CSD->nPSb = 0;
    CSD->mLook = 0;
// Lists = 0; vectors and matrices
	CSD->nDCinPH = 0;
	CSD->xic = 0;
	CSD->xdc = 0;
	CSD->xph = 0;  //18

	CSD->TKval = 0;
    CSD->Psat = 0;
    CSD->Pval = 0;
	CSD->A = 0;
	CSD->Ttol = 0.;
	CSD->Ptol = 0.;
	CSD->dRes1 = 0.;
	CSD->dRes2 = 0.;
    CSD->ICmm = 0;
    CSD->DCmm = 0;
    CSD->DD = 0;
    CSD->denW = 0;
    CSD->epsW = 0;
    CSD->denWg = 0;
    CSD->epsWg = 0;
    CSD->G0 = 0;
    CSD->V0 = 0;
    CSD->S0 = 0;
    CSD->H0 = 0;
    CSD->Cp0 = 0;
    CSD->A0 = 0;
    CSD->U0 = 0;
    CSD->ICNL = 0;
    CSD->DCNL = 0;
    CSD->PHNL = 0;
    CSD->ccIC = 0;
    CSD->ccDC = 0;
    CSD->ccPH = 0;
}

//===============================================================

void TNode::databr_element_to_vtk( fstream& ff, DATABR *CNode_, long int nfild, long int ndx )
{

  TPrintArrays  prar(f_lga+1/*58*/, DataBR_fields, ff);

  switch( nfild )
  {
   case f_NodeHandle: prar.writeValue( CNode_->NodeHandle);
           break;
   case f_NodeTypeHY: prar.writeValue( CNode_->NodeTypeHY);
           break;
   case f_NodeTypeMT: prar.writeValue( CNode_->NodeTypeMT);
           break;
   case f_NodeStatusFMT: prar.writeValue( CNode_->NodeStatusFMT);
           break;
   case f_NodeStatusCH: prar.writeValue( CNode_->NodeStatusCH);
           break;
   case f_IterDone: prar.writeValue( CNode_->IterDone);
           break;
   case f_TK: prar.writeValue( CNode_->TK);
           break;
   case f_P: prar.writeValue( CNode_->P);
           break;
   case f_Vs: prar.writeValue( CNode_->Vs);
           break;
   case f_Vi: prar.writeValue( CNode_->Vi);
           break;
   case f_Ms: prar.writeValue( CNode_->Ms);
           break;
   case f_Mi: prar.writeValue( CNode_->Mi);
           break;
   case f_Hs: prar.writeValue( CNode_->Hs);
           break;
   case f_Hi: prar.writeValue( CNode_->Hi);
           break;
   case f_Gs: prar.writeValue( CNode_->Gs);
            break;
   case f_IS: prar.writeValue( CNode_->IC);
           break;
   case f_pH:prar.writeValue( CNode_->pH);
           break;
   case f_pe: prar.writeValue( CNode_->pe);
           break;
   case f_Eh: prar.writeValue( CNode_->Eh);
           break;
   case f_Tm: prar.writeValue( CNode_->Tm);
           break;
   case f_dt:prar.writeValue( CNode_->dt);
           break;
 #ifdef NODEARRAYLEVEL
   case f_Dif: prar.writeValue( CNode_->Dif);
           break;
   case f_Vt: prar.writeValue( CNode_->Vt);
           break;
   case f_vp: prar.writeValue( CNode_->vp);
           break;
   case f_eps: prar.writeValue( CNode_->eps);
           break;
   case f_Km:  prar.writeValue( CNode_->Km);
           break;
   case f_Kf:  prar.writeValue( CNode_->Kf);
           break;
   case f_S:  prar.writeValue( CNode_->S);
           break;
   case f_Tr:  prar.writeValue( CNode_->Tr);
           break;
   case f_h:  prar.writeValue( CNode_->h);
           break;
   case f_rho:  prar.writeValue( CNode_->rho);
           break;
   case f_al:  prar.writeValue( CNode_->al);
           break;
   case f_at:  prar.writeValue( CNode_->at);
           break;
   case f_av:  prar.writeValue( CNode_->av);
           break;
   case f_hDl:  prar.writeValue( CNode_->hDl);
           break;
   case f_hDt:  prar.writeValue( CNode_->hDt);
           break;
   case f_hDv: prar.writeValue( CNode_->hDv);
           break;
   case f_nto: prar.writeValue( CNode_->nto);
           break;
#endif
   case f_bIC: prar.writeValue(   CNode_->bIC[ndx] );
           break;
   case f_rMB: prar.writeValue(   CNode_->rMB[ndx] );
           break;
   case f_uIC: prar.writeValue(   CNode_->uIC[ndx] );
           break;
   case f_xDC: prar.writeValue(   CNode_->xDC[ndx] );
           break;
   case f_gam: prar.writeValue(   CNode_->gam[ndx] );
           break;
   case f_dll: prar.writeValue(   CNode_->dll[ndx] );
           break;
   case f_dul: prar.writeValue(   CNode_->dul[ndx] );
           break;
   case f_aPH: prar.writeValue(   CNode_->aPH[ndx] );
           break;
   case f_xPH: prar.writeValue(   CNode_->xPH[ndx] );
           break;
   case f_vPS: prar.writeValue(   CNode_->vPS[ndx] );
           break;
   case f_mPS: prar.writeValue(   CNode_->mPS[ndx] );
           break;
   case f_bPS: prar.writeValue(   CNode_->bPS[ndx] );
           break;
   case f_xPA: prar.writeValue(   CNode_->xPA[ndx] );
           break;
  case f_bSP: prar.writeValue(   CNode_->bSP[ndx] );
         break;
  case f_amru: prar.writeValue(   CNode_->amru[ndx] );
         break;
  case f_amrl: prar.writeValue(   CNode_->amrl[ndx] );
         break;
   // CNode_ must be pointer to a work node data bridge structure CNode
  case f_mPH: prar.writeValue(   Ph_Mass(ndx) );
         break;
  case f_vPH: prar.writeValue(   Ph_Volume(ndx) );
         break;
  case f_m_t: prar.writeValue(   Get_mIC( ndx ) );
         break;
  case f_con: prar.writeValue(   Get_cDC(ndx) );
         break;
  case f_mju: prar.writeValue(   Get_muDC(ndx, false ) );
         break;
   case f_lga: prar.writeValue(   Get_aDC(ndx, false ) );
         break;
   default: break;
  }
  ff << endl;
}


void TNode::databr_name_to_vtk( fstream& ff, long int nfild, long int ndx, long int nx2 )
{
  gstring str="", str2="";
  // full name of data fiels (add index name)
  ff << "SCALARS " <<  DataBR_fields[nfild].name.c_str();

  switch( DataBR_fields[nfild].indexation )
  {
    case 1: break;
    case nICbi: str = gstring( CSD->ICNL[ IC_xDB_to_xCH( ndx ) ], 0,MaxICN );
                break;
    case nDCbi: str = gstring( CSD->DCNL[ DC_xDB_to_xCH( ndx ) ], 0,MaxDCN );
              break;
    case nPHbi:
    case nPSbi: str = gstring(  CSD->PHNL[ Ph_xDB_to_xCH( ndx ) ], 0,MaxPHN );
            break;
    case nPSbnICbi:
                str = gstring(  CSD->PHNL[ Ph_xDB_to_xCH( ndx/nx2 ) ], 0,MaxPHN );
                str2 = gstring( CSD->ICNL[ IC_xDB_to_xCH( ndx%nx2 ) ], 0, MaxICN );
          break;
    default: str = gstring( "UNDEFINED");
  }

#ifdef IPMGEMPLUGIN
  strip(str);
  strip(str2);
#else
  str.strip();
  str2.strip();
#endif

  if( !str.empty() )
  { ff << "_" << str.c_str();
    if( !str2.empty() )
      ff << "_" << str2.c_str();
    ff << "_";
  }

  if( nfild < 6)
     ff << " long";
  else
     ff << " double";
  ff << " 1" << endl;

  ff << "LOOKUP_TABLE default" << endl;
}

void TNode::databr_size_to_vtk(  long int nfild, long int& nel, long int& nel2 )
{

    nel = 1;
    nel2 = 1;
    switch( DataBR_fields[nfild].indexation )
    {
      case 1: break;
      case nICbi: nel = CSD->nICb;
                  break;
      case nDCbi: nel = CSD->nDCb;
                  break;
      case nPHbi: nel = CSD->nPHb;
                  break;
      case nPSbi: nel = CSD->nPSb;
                  break;
      case nPSbnICbi: nel = CSD->nPSb; nel2 = CSD->nICb;
        break;
   }

}

void TNode::databr_head_to_vtk( fstream& ff, const char*name, double time, long cycle,
                               long nx, long ny, long nz )
{
 ff << "# vtk DataFile Version 3.0" <<  endl;
 ff << "GEM2MT " << name <<  endl;
 ff << "ASCII" <<  endl;
 ff << "DATASET STRUCTURED_POINTS" <<  endl;
 ff << "DIMENSIONS " <<  nx << " " << ny << " " << nz << endl;
 ff << "ORIGIN 0.0 0.0 0.0" <<  endl;
 ff << "SPACING 1.0 1.0 1.0" <<  endl;
 //????
 ff << "FIELD TimesAndCycles 2" << endl;
 ff << "TIME 1 1 double" <<  endl << time << endl;
 ff << "CYCLE 1 1 long" <<  endl << cycle << endl;
 ff << "POINT_DATA " << nx*ny*nz << endl;
}

void TNode::databr_to_vtk( fstream& ff, const char*name, double time, long int  cycle,
                          long int  nFilds, long int  (*Flds)[2])
{
   bool all = false;
   long int kk, ii, nf, nel, nel2;

   // write header of file
   databr_head_to_vtk( ff, name, time, cycle, 1, 1, 1 );

   if( nFilds < 1 || !Flds )
   {
      all = true;
      nFilds = 51;
   }

   for( kk=0; kk<nFilds; kk++)
   {
       if( all )
         nf = kk;
       else nf= Flds[kk][0];

       databr_size_to_vtk(  nf, nel, nel2 );

       if( all )
         { ii=0; }
       else
         { ii = Flds[kk][1];
           nel = ii+1;
         }

       for( ; ii<nel; ii++ )
       {
        databr_name_to_vtk( ff, nf, ii, nel2 );
        // cycle for TNode array
        databr_element_to_vtk( ff, CNode, nf, ii );
       }
   }
}

//-----------------------End of node_format.cpp--------------------------
