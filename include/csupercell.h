//========================================================================
// FILE: csupercell.h
//
// This an utility program to manipulate and genarate structure files
//
// Copyright 2006-2015 by Edmanuel Torres
// email:   eetorres@gmail.com
//
// Lastest update: 13/07/2015
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

#ifndef _CSUPERCELL_H_
#define _CSUPERCELL_H_

#include<config_debug.h>
#include<cfragment.h>
#include<cfile.h>

// Neighbouring cells
const int neighbor_cells[27][3] = {
{  0, 0, 0}, // 0
{  1, 0, 0}, // 1
{  1, 1, 0}, // 2
{  0, 1, 0}, // 3
{ -1, 1, 0}, // 4
{  0, 0, 1}, // 5
{  1, 0, 1}, // 6
{ -1, 0, 1}, // 7
{  1, 1, 1}, // 8
{  0, 1, 1}, // 9
{ -1, 1, 1}, // 10
{ -1,-1, 1}, // 11
{  0,-1, 1}, // 12
{  1,-1, 1}, // 13
{  0, 0,-1}, // 14
{  1, 0,-1}, // 15
{ -1, 0,-1}, // 16
{  1, 1,-1}, // 17
{  0, 1,-1}, // 18
{ -1, 1,-1}, // 19
{ -1,-1,-1}, // 20
{  0,-1,-1}, // 21
{  1,-1,-1}, // 22
{ -1, 0, 0}, // 23
{ -1,-1, 0}, // 24
{  0,-1, 0}, // 25
{  1,-1, 0}  // 26
};

class CSupercell {

public:

  CSupercell();
  ~CSupercell();
  //
  CFile gsf;
  //
  void clear(void);
  bool read_input_file(void);                  // read the input file
  void initialize_fragments(void);             //
  //
  /////////////////////////////////////////////////////////////////////////////////
  // Save functions
  /////////////////////////////////////////////////////////////////////////////////
  void save_input_file(void);                  // save the file
  void save_as_file(std::string,bool,bool,TMatrix<real> m,int,int,int,int,bool);
  void save_as_file(std::string,std::string,bool,bool,TMatrix<real> m,int,int,int,int,bool);
  void save_topmol_file(std::string,std::string);

  /////////////////////////////////////////////////////////////////////////////////
  // Compute functions
  /////////////////////////////////////////////////////////////////////////////////
  //void compute_position_direct(void);
  void compute_position_cartesian(void);
  //
  /////////////////////////////////////////////////////////////////////////////////
  // Evaluation functions
  /////////////////////////////////////////////////////////////////////////////////
  void eval_initial_position(void);
  void eval_initial_orientation(void);
  //void eval_cell_table(void);
  void eval_connections(uint);
  bool eval_new_fragment(const TVector<uint>&);               // Create all posible scaled framgent
  bool eval_scaled_fragment(uint,bool,real);      // Create a framgent using scaled distance
  bool eval_radial_fragment(uint,bool,real);      // Create a framgent using scaled distance
  bool eval_scaled_fragment(uint,real);           // Create single scaled distance fragment
  void eval_scaled_fragments(real);               // Create all posible scaled framgent
  // Begin deprecated
  void eval_vdw_fragments(void);               // Create all vdW framgents
  void eval_atom_fragments(void);               // Create all atom framgents
  // End deprecated
  bool eval_merge_fragment(uint,bool);         // Merge two fragments
  //
  bool delete_atom(uint u);                    // Delete an atom
  //
  void update_fragmol_cartesian(void);
  void update_fragmol_direct(void);
  //
  void is_fragmol_initialized(bool);
  void apply_pbc(bool);

  /////////////////////////////////////////////////////////////////////////////////
  // Set functions
  /////////////////////////////////////////////////////////////////////////////////
  void set_cartesian(void);
  //
  void set_gsf_modified(const bool);
  //
  void set_active_fragment(const uint);
  void set_input_file_type(const uint);
  void set_input_file_units(const uint);
  void set_output_file_type(const uint);
  void set_output_file_format(const uint);
  void set_export_format(const uint);
  //
  void set_fragment_twist(const real);
  void set_fragmol_fragment_precession(const real);
  void set_fragmol_fragment_tilt(const real);
  //
  void set_fragmol_fragment_position_u(const real);
  void set_fragmol_fragment_position_v(const real);
  void set_fragmol_fragment_position_w(const real);
  //
  void set_input_file(std::string);
  void set_dir(std::string);
  void set_topmol_directory(std::string);
  //
  void set_fragment_axis(const TVector<uint>&);
  void set_fragment_axis(const uint, const uint);
  void set_fragment_direct(const uint,const TMatrix<real>&);
  void set_fragment_cartesian(const uint,const TMatrix<real>&);
  //
  /////////////////////////////////////////////////////////////////////////////////
  // Is  functions
  /////////////////////////////////////////////////////////////////////////////////
  bool is_direct(void);
  bool is_periodic(void);
  void is_periodic(bool);
  bool is_potmol(void);
  bool is_fragmol_initialized();

  /////////////////////////////////////////////////////////////////////////////////
  //  Get  functions
  /////////////////////////////////////////////////////////////////////////////////
  std::string get_dir(void);
  std::string get_fragment_atomic_label(uint,uint);
  std::string get_fragment_atomic_symbol(uint,uint);
  //
  uint get_number_of_fragments(void);
  uint get_fragmol_active_fragment(void);
  uint get_fragment_size(uint);
  uint get_fragment_atomic_number(uint,uint);
  //
  uint get_fragmol_total_atoms(void);
  uint get_atomic_species(void);
  uint get_atomic_composition(uint);
  uint get_total_atoms(void);
  uint get_input_file_format(void);
  uint get_output_file_type(void);
  uint get_output_file_format(void);
  //
  real get_fragmol_axis_tilt(void);
  real get_fragmol_axis_precession(void);
  real get_fragmol_backbone_tilt(void);
  real get_fragmol_backbone_precession(void);
  //
  TVector<std::string> get_atomic_symbol_table(void);
  std::string get_atomic_symbol_table(uint);
  //
  TVector<uint> get_fragmol_atomic_composition_table(void);
  TVector<uint> get_atomic_number_table(void);
  uint get_atomic_number_table(uint);
  TVector<uint> get_atom_table(void);
  uint get_atom_table(uint);
  TVector<uint> get_fragment_table(void);
  void eval_atomic_bonds(void);
  uint get_fragment_table(uint);
  //
  TVector<real> get_fragmol_axis_angles(void);
  TVector<real> get_fragmol_basis_direct(void);
  TVector<real> get_position_direct(void);
  TVector<real> get_fragmol_position_uvw(void);
  TVector<real> get_position_cartesian(void);
  TVector<real> get_fragment_centered_position_cartesian(void);
  //
  TVector<real> get_cartesian(uint);
  TVector<real> get_direct(uint);
  TVector<real> get_vector_diff(uint,uint);
  //
  real get_distance(uint,uint);
  real get_angle(uint,uint,uint);
  real get_dihedral(uint,uint,uint,uint);
  //
  TVector<real> get_fragment_direct(uint,uint);
  TVector<real> get_fragment_cartesian(uint,uint);
  TVector<real> get_fragment_centered_cartesian(uint,uint);
  TMatrix<real> get_cartesian(void);
  TMatrix<real> get_direct(void);
  TMatrix<real> get_fragment_direct(uint);
  TMatrix<real> get_fragment_cartesian(uint);
  TMatrix<real> get_fragment_centered_cartesian(uint);
  TMatrix<real> get_uvw_to_xyz(void);
  TVector<real> get_uvw_to_xyz(uint);
  TMatrix<real> get_unit_uvw_to_xyz(void);
  TMatrix<real> get_bounding_box(void);

  // MD special functions
  // cell list
  int  i_neighbor_cells;
  bool b_linked_cell;
  uint u_cell_number;
  bool is_linked_cell(void){ return b_linked_cell;};
  void is_linked_cell(bool b){ b_linked_cell=b;};
  void set_cells(void);
  void set_inverse_cell(void);
  void set_cell_list(void);
  void eval_linked_list(void);
  int  get_cell_list(int i){ return v_cell_list[i];};
  int  get_cell_head(int i){ return v_cell_head[i];};
  int  get_neighbor_cells(void){ return i_neighbor_cells;};
  uint get_cell_number(void){ return u_cell_number;};
  void set_bbox(const TVector<real>& v){ v_bbox = v;};
  real get_bbox(uint u){ return v_bbox[u];};
  TVector<real> get_bbox(void){ return v_bbox;};
  TVector<real> get_radius_color(uint u){ return m_radius_color[u];};
  real get_radius_color(uint u1, uint u2){ return m_radius_color[u1][u2];};
  void set_box_size(const TVector<real>& v){ v_box_size = 2.0*v;};
  TVector<real> get_box_size(void){ return v_box_size;};
  TVector<real> get_cell_frac(void){ return v_cell_frac;};
  TVector<int>  get_neighbor_cells_xyz(uint u){ return neighbor_cells_xyz[u];};
  TMatrix<real> get_inv_bbox(void){ return m_inv_bbox;};

  void set_inv_bbox(void);
  void set_radius_color(void);
  TMatrix<real> get_radius_color(void);
  // Temporal public variables
  TVector<real> v_bbox;
  TVector<int>  v_cell_side, v_cell_list, v_cell_head;
  TMatrix<real> m_inv_bbox;
  TMatrix<int>  neighbor_cells_xyz;
  // bond functions
  uint get_number_of_bonds(void){ return i_number_of_bonds;};
  uint get_number_of_bonds_pbc(void){ return i_number_of_bonds_pbc;};
  uint get_bond_number(uint u){ return v_bond_number[u];};
  uint get_bond_number_pbc(uint u){ return v_bond_number_pbc[u];};
  uint get_bond_index(const real f){ return (uint)(f/f_atom_bond_delta);};
  int  get_bond_boundary_pbc(uint u1, uint u2){ return m_bond_boundary_pbc[u1][u2];};
  uint get_bond_indices(uint u1, uint u2){ return m_bond_indices[u1][u2];};
  uint get_bond_indices_pbc(uint u1, uint u2){ return m_bond_indices_pbc[u1][u2];};
  uint get_bond_table(uint u){ return v_bond_table[u];};
  uint get_bond_types(void){ return v_bond_table.size();};
  inline uint check_bond(uint);
  // bond variables
  real f_atom_bond_delta;
  TVector<uint> v_bond_table;
  uint i_number_of_bonds;
  uint i_number_of_bonds_pbc;
  TVector<uint> v_bond_number;
  TVector<uint> v_bond_number_pbc;
  TMatrix<int>  m_bond_boundary_pbc;
  TMatrix<real> m_radius_color;
  TMatrix<uint> m_bond_indices;
  TMatrix<uint> m_bond_indices_pbc;
  // Temporal functions
  void set_cut_radius(real r){ r_cut_radius=r; r_cut_radius_2=r*r;};

private:

  bool __is_direct;
  bool __is_periodic;
  bool __is_potmol;
  bool __is_input_fragment;
  //
  uint __input_format;
  uint __output_format;
  uint __total_atoms;
  uint __atomic_species;
  uint __active_fragment;
  //
  real r_cut_radius, r_cut_radius_2;
  //
  std::ifstream iposmol;
  std::ifstream itopmol;
  std::ofstream oposmol;
  std::ofstream otopmol;
  ///
  TVector<real> v_cell_frac, v_box_size, v_box_middle;
  //
  TMatrix<real> m_xyz;
  TMatrix<real> m_uvw;
  //
  std::string inputfile, output_filename, potmolfile, __sdir;
  //
};

#endif

