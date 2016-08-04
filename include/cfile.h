//========================================================================
// FILE: cread.h
//
// This an utility program to manipulate and genarate structure files
//
// Copyright 2011-2016 by Edmanuel Torres
// email:   eetorres@gmail.com
//
//========================================================================
//  This file is part of xmolview                                       //
//                                                                      //
//  xmolview is free software: you can redistribute it and/or modify    //
//  it under the terms of the GNU General Public License as published by//
//  the Free Software Foundation, either version 3 of the License, or   //
//  (at your option) any later version.                                 //
//                                                                      //
//  xmolview is distributed in the hope that it will be useful,         //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of      //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       //
//  GNU General Public License for more details.                        //
//                                                                      //
//  You should have received a copy of the GNU General Public License   //
//  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.     //
//======================================================================//

#ifndef _CFILE_H_
#define _CFILE_H_

#include <config_debug.h>
#include <global.h>
#include <ctopmol.h>
#include <cfragment.h>
#include <cpotcar.h>
#include <cposcar.h>
#include <cxyz.h>
#include <cgau.h>
#include <czmat.h>
#include <cpdb.h>
#include <cdlp.h>
#include <czmat.h>

class CFile : public CTopmol{

public:

  CFile();
  ~CFile(){};
  // NOTE: public members while debuging
  TVector<CFragment> v_fragments;
  //
  void clear(void);
  void clear_xyz(void);
  void clear_gau(void);
  void clear_pdb(void);
  void clear_dlp(void);
  void clear_poscar(void);
  void clear_potcar(void);
  //
  bool read_input_file(void);
  bool read_poscar(void);
  bool read_xyz(void);
  bool read_gau(void);
  bool read_pdb(void);
  bool read_dlp(void);
  bool read_zmt(void);
  bool read_potcar(void);
  //
  void save_input_file(void);
  void save_xyz(void);
  void save_as_file(std::string);
  void save_as_file(std::string,std::string);
  void save_poscar_as(std::string,std::string,uint);
  void save_xyz_as(std::string,std::string,uint);
  void save_gau_as(std::string,std::string,uint);
  void save_pdb_as(std::string,std::string,uint);
  void save_dlp_as(std::string,std::string,uint);
  void save_zmt_as(std::string,std::string,uint);
  //
  void init_fragments(void);
  void clear_fragments(void);
  void cast_fragments(void);
  void eval_connections(const TMatrix<uint>&,uint);
  void eval_fragments(void);
  bool eval_new_fragment(const TVector<uint>&);
  bool eval_scaled_fragment(const uint,bool,real);
  bool eval_scaled_fragment(const uint,real);
  void eval_scaled_fragments(real);
  void eval_vdw_fragments(void);
  void comp_all_cartesian(void);
  void comp_cartesian(uint);
  void comp_all_direct(void);
  void comp_direct(uint);
  void update_cell_table(void);
  void update_cartesian(void);
  void update_direct(void);
  //
  void set_input_file(std::string);
  void set_output_file_name(std::string);
  void set_potcar_file(std::string);
  void set_input_dir(std::string);
  std::string get_input_dir(void);
  //
  void set_input_type(uint);
  void set_input_units(uint);
  void set_input_format(uint);
  //
  void set_output_file_type(uint);
  void set_output_file_format(uint);
  void set_export_format(uint);
  //
  void set_active_fragment(uint u){ u_active_fragment=u;};
  uint get_active_fragment(void){ return u_active_fragment;};
  //
  void set_bounding_box(bool);
  void set_labels(bool);
  void set_numbers(bool);
  void set_fragments(bool);
  void set_fragment_table(uint u){v_fragment_table.resize(u);};
  void set_fragment_table(uint u, uint v){v_fragment_table[u]=v;};
  void set_fragment_table(TVector<uint> v, uint u);
  void set_modified(bool b);
  void set_is_periodic(bool b){ b_periodic=b;}
  //
  void set_xyz(TMatrix<real> m);
  void set_xyz_cells(int x, int y, int z, int t);

  uint get_total_atoms(void){ return u_total_atoms;}
  uint get_atomic_species(void){ return u_atomic_species;}
  uint get_input_type(void){ return u_input_file_type;}
  uint get_input_units(void){ return u_input_file_units;}
  uint get_input_format(void){ return u_input_format;}
  uint get_output_file_format(void){ return u_output_file_format;}
  uint get_output_file_type(void){ return u_output_file_type;}
  uint get_number_of_fragments(void){ return u_number_of_fragments;}
  //
  TVector<uint> get_atomic_composition_table(void){ return v_atomic_composition_table;}
  TVector<uint> get_atomic_numbers(void){ return v_atomic_numbers;}
  //
  TVector<uint> get_fragment_table(void){ return v_fragment_table;}
  uint get_fragment_table(uint u){ return v_fragment_table[u];}

  TVector<uint> get_atom_table(void){ return v_atom_table;}
  uint get_atom_table(uint u){ return v_atom_table[u];}

  TVector<uint> get_atomic_number_table(void){ return v_atomic_number_table;}
  uint get_atomic_number_table(uint u){ return v_atomic_number_table[u];}
  //
  TVector<std::string> get_atomic_labels(void){ return v_atomic_labels;}
  strg get_atomic_label(uint u){ return v_atomic_labels[u];}
  //
  TVector<std::string> get_atomic_symbols(void){ return v_atomic_symbols;}
  strg get_atomic_symbol(uint u){ return v_atomic_symbols[u];}
  //
  TVector<std::string> get_atomic_symbol_table(void){ return v_atomic_symbol_table;}
  strg get_atomic_symbol_table(uint u){ return v_atomic_symbol_table[u];}
  //  
  uint get_fragment_size(uint u){ return v_fragments[u].size();}
  uint get_atomic_composition(uint u){ return v_atomic_composition_table[u];}
  //
  uint get_atomic_number(uint u){ return v_atomic_numbers[u];}
  uint get_atom_cell_table(uint u){ return v_atom_cell_table[u];}
  // fragment tools
  void eval_initial_position(void){
    v_fragments[u_active_fragment].eval_initial_position();
  }
  void eval_initial_orientation(void){
    v_fragments[u_active_fragment].eval_initial_orientation();
  }
  void compute_origin_cartesian(void){
    v_fragments[u_active_fragment].compute_origin_cartesian();
  }
  void is_initialized(bool b){
    v_fragments[u_active_fragment].is_initialized(b);
  }
  void is_pbc(bool b){
    v_fragments[u_active_fragment].is_pbc(b);
  }
  void set_axis_index(const TVector<uint>& _v){
    v_fragments[u_active_fragment].set_axis_index(_v);
  }
  void set_axis_index(const uint u1, const uint u2){
    v_fragments[u_active_fragment].set_axis_index(u1,get_atom_cell_table(u2));
  }
  void set_axis_twist(real r){
    v_fragments[u_active_fragment].set_axis_twist(r);
  }
  void set_axis_precession(real r){
    v_fragments[u_active_fragment].set_axis_precession(r);
  }
  void set_axis_tilt(real r){
    v_fragments[u_active_fragment].set_axis_tilt(r);
  }
  void set_origin_u(real r){
    v_fragments[u_active_fragment].set_origin_u(r);
  }
  void set_origin_v(real r){
    v_fragments[u_active_fragment].set_origin_v(r);
  }
  void set_origin_w(real r){
    v_fragments[u_active_fragment].set_origin_w(r);
  }
  real get_backbone_tilt(void){
    return v_fragments[u_active_fragment].get_backbone_tilt();
  }
  real get_backbone_precession(void){
    return v_fragments[u_active_fragment].get_backbone_precession();
  }
  real get_axis_tilt(void){
    return v_fragments[u_active_fragment].get_axis_tilt();
  }
  real get_axis_precession(void){
    return v_fragments[u_active_fragment].get_axis_precession();
  }
  bool is_initialized(void){
    return v_fragments[u_active_fragment].is_initialized();
  }
  TVector<real> get_axis_angles(){
    return v_fragments[u_active_fragment].get_axis_angles();
  }
  TVector<real> get_basis_direct(){
    return v_fragments[u_active_fragment].get_basis_direct();
  }
  TVector<real> get_origin_direct(){
    return v_fragments[u_active_fragment].get_origin_direct();
  }
  TVector<real> get_origin_uvw(){
    return v_fragments[u_active_fragment].get_origin_uvw();
  }
  TVector<real> get_origin_cartesian(){
    return v_fragments[u_active_fragment].get_origin_cartesian();
  }
  TVector<real> get_centered_origin_cartesian(void){
    return v_fragments[u_active_fragment].get_centered_origin_cartesian();
  }
  TVector<real> get_fragment_centered_cartesian(uint i,uint j){
    return v_fragments[i].get_centered_cartesian(j);
  }
  TVector<real> get_fragment_direct(uint i,uint j){
    return v_fragments[i].get_direct(j);
  }
  TVector<real> get_fragment_cartesian(uint i,uint j){
    return v_fragments[i].get_cartesian(j);
  }
  //
  void is_direct(bool b){ b_direct=b;}
  bool is_direct(void){ return b_direct;}
  void is_periodic(bool b){ b_periodic=b;}
  bool is_periodic(void){ return b_periodic;}
  bool get_is_direct(void){ return b_direct;}
  bool get_is_periodic(void){ return b_periodic;}
  //
  TMatrix<real> get_xyz(void){ return m_xyz;}
  TVector<real> get_xyz(uint u){ return m_xyz[u];}
  TMatrix<real> get_uvw(void){ return m_uvw;}
  TVector<real> get_uvw(uint u){ return m_uvw[u];}
  TMatrix<real> get_uvw_to_xyz_u(void){ return m_uvw_to_xyz_u;}
  TMatrix<real> get_uvw_to_xyz(void){ return m_uvw_to_xyz;}
  TVector<real> get_uvw_to_xyz(uint u){ return m_uvw_to_xyz[u];}

private:

  bool b_direct;
  bool b_periodic;
  bool is_potcar;
  bool is_bounding_box;
  bool is_labels;
  bool is_numbers;
  bool is_fragments;
  bool is_charges;
  bool is_modified;
  //
  int x_cells;
  int y_cells;
  int z_cells;
  int total_cells;
  //
  uint u_input_format;
  uint u_input_file_type;
  uint u_input_file_units;
  uint u_output_file_type;
  uint u_output_file_format;
  uint u_export_format;
  uint u_total_atoms;
  uint u_active_fragment;
  uint u_total_fragments;
  uint u_atomic_species;
  uint u_number_of_fragments;
  //
  std::string inputfile;
  std::string output_filename;
  std::string potcarfile;
  std::string s_actual_dir;
  //
  CPotcar file_potcar;
  CPoscar file_poscar;
  CXyz    file_xyz;
  CGau    file_gau;
  CZmat   file_zmt;
  CPdb    file_pdb;
  CDlp    file_dlp;
  //
  // TVector<CFragment> v_fragments;
  //
  TVector<std::string> v_atomic_symbol_table;
  TVector<uint> v_atomic_composition_table;
  TVector<uint> v_atomic_number_table;
  TVector<uint> v_fragment_table;
  TVector<uint> v_atom_table;
  TVector<uint> v_atomic_numbers;
  TVector<uint> v_atom_cell_table;
  TVector<real> v_atomic_charges;
  TVector<int>  v_connections;
  //
  TVector<std::string> v_atomic_symbols;
  TVector<std::string> v_atomic_labels;
  //
  TMatrix<real> m_file_input;
  TMatrix<real> m_xyz;
  TMatrix<real> m_uvw;
  TMatrix<real> m_uvw_to_xyz;
  TMatrix<real> m_uvw_to_xyz_u;
};

#endif
