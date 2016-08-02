//======================================================================//
// FILE: csupercell.cxx -> csupercell                                   //
//                                                                      //
// Abstraction layer of the supercell.                                  //
// The manipulation of individual fragments is hidden.                  //
//                                                                      //
// Copyright 2011-2016 by Edmanuel Torres                               //
// email:   eetorres@gmail.com                                          //
//                                                                      //
// Comment:                                                             //
//                                                                      //
//======================================================================//
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

#include <config_debug.h>
#include <csupercell.h>

CSupercell::CSupercell(){
  __is_potmol=false;
  f_atom_bond_delta          = 0.1;
  clear();
}

CSupercell::~CSupercell(void){
}

void CSupercell::clear(void){
  __is_potmol=false;
  __active_fragment=0;
  is_linked_cell(true);
}

bool CSupercell::read_input_file(void){
  bool res=false;
  clear();
  // check if a POTCAR is available
  gsf.read_potcar();
  // Read the input structure file
  res = gsf.read_input_file();
  if(res){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: total atoms = "<<gsf.get_total_atoms()<<std::endl;
    std::cout<<" FRAGMOL: total species = "<<gsf.get_atomic_species()<<std::endl;
    std::cout<<" FRAGMOL: is direct? "<<gsf.get_is_direct()<<std::endl;
    std::cout<<" FRAGMOL: is periodic? "<<gsf.get_is_periodic()<<std::endl;
    std::cout<<" FRAGMOL: input format = "<<gsf.get_input_format()<<std::endl;
#endif
    //
    //v_atomic_composition_table=gsf.get_atomic_composition_table();
    //v_atomic_number_table=gsf.get_atomic_number_table();
    //v_atomic_symbol_table=gsf.get_atomic_symbol_table();
    //v_atom_type_table=gsf.get_atom_table();
    //v_atomic_labels=gsf.get_atomic_labels();
    //v_atomic_symbols=gsf.get_atomic_symbols();
    //v_atomic_numbers=gsf.get_atomic_numbers();
    //v_atom_cell_table.resize(gsf.get_total_atoms());
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    //std::cout<<" FRAGMOL: Composition ="<<v_atomic_composition_table;
    //std::cout<<" FRAGMOL: Z-Number Table ="<<v_atomic_number_table;
    //std::cout<<" FRAGMOL: Symbol Table ="<<v_atomic_symbol_table;
    //std::cout<<" FRAGMOL: Type ="<<v_atom_type_table;
    //std::cout<<" FRAGMOL: Labels ="<<v_atomic_labels;
    //std::cout<<" FRAGMOL: Symbols ="<<v_atomic_symbols;
    //std::cout<<" FRAGMOL: Z-Numbers ="<<v_atomic_numbers;
#endif
    //
    m_xyz=gsf.get_xyz();
    m_uvw=gsf.get_uvw();
#ifdef _FRAGMOL_DATA_MESSAGES_
    //std::cout<<" FRAGMOL: XYZ = "<<m_xyz;
    //std::cout<<" FRAGMOL: UVW = "<<m_uvw;
    //std::cout<<" FRAGMOL: uvwTxyz U = "<<m_uvw_to_xyz_u;
    //std::cout<<" FRAGMOL: uvwTxyz  = "<<m_uvw_to_xyz;
#endif
    return res;
  }
  return res;
}

bool CSupercell::delete_atom(uint u){
  bool res=false;
  return res;
}

// the code below should be migrated to a new cread.h file API
void CSupercell::save_input_file(void){
  gsf.save_input_file();
}

void CSupercell::save_as_file(std::string _f, bool _l, bool _n,TMatrix<real> m, int xc, int yc, int zc, int tc, bool _bb){
  //if(__output_format == OUTPUT_FORMAT_ATM_FRG || __output_format == OUTPUT_FORMAT_NAT_FRG )
    //gsf.set_fragment_table(v_fragment_table,__number_of_fragments);
  gsf.set_bounding_box(_bb);
  gsf.set_labels(_l);
  gsf.set_numbers(_n);
  gsf.set_xyz_cells(xc,yc,zc,tc);
  gsf.set_xyz(m);
  gsf.save_as_file(_f);
}

void CSupercell::save_as_file(std::string _p, std::string _f, bool _l, bool _n, TMatrix<real> m, int xc, int yc, int zc, int tc, bool _bb){
  //if(__output_format == OUTPUT_FORMAT_ATM_FRG || __output_format == OUTPUT_FORMAT_NAT_FRG )
    //gsf.set_fragment_table(v_fragment_table,__number_of_fragments);
  gsf.set_bounding_box(_bb);
  gsf.set_labels(_l);
  gsf.set_numbers(_n);
  gsf.set_xyz_cells(xc,yc,zc,tc);
  gsf.set_xyz(m);
  gsf.save_as_file(_p,_f);
}

void CSupercell::initialize_fragments(void){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGCAR: initialize fragments"<<std::endl;
#endif
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGCAR: !!! No topology definition to use !!!"<<std::endl;
  std::cout<<" FRAGCAR: Build a single fragment system"<<std::endl;
#endif
  gsf.init_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGCAR:  ["<<get_number_of_fragments()<<"] fragments loaded"<<std::endl;
#endif
}

/*
void CSupercell::eval_cell_table(void){
  TVector<uint> v_l;
  uint _s;
  for(uint i=0;i<gsf.get_number_of_fragments();i++){
    _s=gsf.get_fragment_size(i);
    v_l = gsf.get_topology_atoms(i);
    for(uint j=0;j<_s;j++){
      // store the position of the current atom inside the fragment as a table
      v_atom_cell_table[v_l[j]]=j;
    }
  }
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: updated atom cell table: "<<v_atom_cell_table;
#endif
}
*/

void CSupercell::eval_connections(uint u){
  gsf.eval_connections(m_bond_indices,u);
}

void CSupercell::eval_initial_position(void){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: eval_initial_position"<<std::endl;
#endif
  gsf.v_fragments[__active_fragment].eval_initial_position();
}

void CSupercell::eval_initial_orientation(void){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: eval_initial_orientation"<<std::endl;
#endif
  gsf.v_fragments[__active_fragment].eval_initial_orientation();
}

void CSupercell::set_cartesian(void){
  TVector<real> _v;
  TVector<unsigned int> v_l;
  unsigned int _n, _s;
  _n = gsf.get_number_of_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<"FRAGMOL: TOTAL FRAGMENTS: "<<_n<<std::endl;
#endif
  for(unsigned int i=0; i<_n; i++){
    _s = get_fragment_size(i);
    v_l = gsf.get_topology_atoms(i);
    for(unsigned int j=0;j<_s;j++){
      _v = get_fragment_cartesian(i,j);
      // put each atomic coordinate in the right place inside the matrix
      m_xyz[v_l[j]]=_v;
    }
  }
#ifdef _FRAGMOL_FORMAT_MESSAGES_
  std::cout<<" FRAGMOL: Cartesian coordinates ready"<<std::endl;
#endif
#ifdef _FRAGMOL_DATA_MESSAGES_
  std::cout<<" FRAGMOL: m_cartesian: "<<m_xyz;
#endif
}

void CSupercell::update_fragmol_cartesian(void){
  TVector<real> _v;
  TVector<unsigned int> v_l;
  unsigned int _n, _s;
  _n = gsf.get_number_of_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<"FRAGMOL: TOTAL FRAGMENTS: "<<_n<<std::endl;
#endif
  for(unsigned int i=0; i<_n; i++){
    _s = get_fragment_size(i);
    v_l = gsf.get_topology_atoms(i);
    for(unsigned int j=0;j<_s;j++){
      _v = get_fragment_centered_cartesian(i,j);
      // put each atomic coordinate in the right place inside the matrix
      m_xyz[v_l[j]]=_v;
    }
  }
#ifdef _FRAGMOL_FORMAT_MESSAGES_
  std::cout<<" FRAGMOL: Cartesian coordinates ready"<<std::endl;
#endif
#ifdef _FRAGMOL_DATA_MESSAGES_
  std::cout<<" FRAGMOL: centered cartesian: "<<m_xyz;
#endif
}

void CSupercell::update_fragmol_direct(void){
  TVector<real> _v;
  TVector<unsigned int> v_l;
  unsigned int _n, _s;
  _n = gsf.get_number_of_fragments();
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<"FRAGMOL: TOTAL FRAGMENTS: "<<_n<<std::endl;
#endif
  for(unsigned int i=0; i<_n; i++){
    _s = get_fragment_size(i);
    v_l = gsf.get_topology_atoms(i);
    for(unsigned int j=0;j<_s;j++){
      _v = get_fragment_direct(i,j);
      // put each atomic coordinate in the right place inside the matrix
      m_uvw[v_l[j]]=_v;
    }
  }
#ifdef _FRAGMOL_FORMAT_MESSAGES_
  std::cout<<" FRAGMOL: Direct coordinates ready"<<std::endl;
#endif
#ifdef _FRAGMOL_DATA_MESSAGES_
  std::cout<<" FRAGMOL: m_direct: "<<m_uvw;
#endif
}

void CSupercell::compute_position_cartesian(void){
  gsf.v_fragments[__active_fragment].compute_origin_cartesian();
}

void CSupercell::is_fragmol_initialized(bool b){
  gsf.v_fragments[__active_fragment].is_initialized(b);
}

void CSupercell::apply_pbc(bool b){
  gsf.v_fragments[__active_fragment].is_pbc(b);
}

void CSupercell::set_radius_color(void){
  int i_z=0;
  TVector<real> _rcolor(4);
  m_radius_color.resize(get_total_atoms(),4);
  for(int i=0; i<get_total_atoms(); i++){
      i_z=gsf.get_atomic_number(i);
      _rcolor[0] = atom_rrgb[i_z][0];
      _rcolor[1] = atom_rrgb[i_z][1];
      _rcolor[2] = atom_rrgb[i_z][2];
      _rcolor[3] = atom_rrgb[i_z][3];
      m_radius_color[i]=_rcolor;
  }
  //m_radius_color=m;
}

TMatrix<real> CSupercell::get_radius_color(void){
  return m_radius_color;
}

void CSupercell::set_gsf_modified(const bool b){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: the gsf was modified"<<std::endl;
#endif
  gsf.set_modified(b);
}

void CSupercell::set_active_fragment(const uint i){
  __active_fragment=i;
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: active fragment= "<<__active_fragment<<std::endl;
#endif
}

void CSupercell::set_input_file_type(const uint u){
  //__input_type=u;
  gsf.set_input_type(u);
}

void CSupercell::set_input_file_units(const uint u){
  //__input_type=u;
  gsf.set_input_units(u);
}

void CSupercell::set_output_file_type(const uint u){
  //__output_file_type=u;
  gsf.set_output_file_type(u);
}

void CSupercell::set_output_file_format(const uint u){
  __output_format=u;
  gsf.set_output_file_format(u);
}

void CSupercell::set_export_format(const uint u){
  //__export_format=u;
//#ifdef _FRAGMOL_DEBUG_MESSAGES_
//  std::cout<<" FRAGMOL: export format = "<<__export_format<<std::endl;
//#endif
  gsf.set_export_format(u);
}

void CSupercell::set_fragment_axis(const TVector<uint>& _v){
  is_fragmol_initialized(false);
  gsf.v_fragments[__active_fragment].set_axis_index(_v);
}

void CSupercell::set_fragment_axis(const uint idx, const uint val){
  is_fragmol_initialized(false);
  gsf.v_fragments[__active_fragment].set_axis_index(idx,gsf.get_atom_table(val));
}

void CSupercell::set_fragment_twist(const real r){
  gsf.v_fragments[__active_fragment].set_axis_twist(r);
}

void CSupercell::set_fragmol_fragment_precession(const real r){
  gsf.v_fragments[__active_fragment].set_axis_precession(r);
}

void CSupercell::set_fragmol_fragment_tilt(const real r){
  gsf.v_fragments[__active_fragment].set_axis_tilt(r);
}

void CSupercell::set_fragmol_fragment_position_u(const real r){
  gsf.v_fragments[__active_fragment].set_origin_u(r);
}

void CSupercell::set_fragmol_fragment_position_v(const real r){
  gsf.v_fragments[__active_fragment].set_origin_v(r);
}

void CSupercell::set_fragmol_fragment_position_w(const real r){
  gsf.v_fragments[__active_fragment].set_origin_w(r);
}

void CSupercell::set_input_file(std::string s){
  gsf.set_input_file(s);
  gsf.topmol_filename(s);
}

void CSupercell::set_dir(std::string s){
  // use__sdir in case you would like to change the default directory
  __sdir = s;
  gsf.set_input_dir(s);
}


void CSupercell::set_topmol_directory(std::string s){
  gsf.topmol_dir(s);
}

void CSupercell::save_topmol_file(std::string _d,std::string _f){
    gsf.save_topmol(_d,_f);
}

void CSupercell::set_fragment_direct(const uint i, const TMatrix<real>& _m){
  gsf.v_fragments[i].set_position_direct(_m);
}
void CSupercell::set_fragment_cartesian(const uint i, const TMatrix<real>& _m){
  gsf.v_fragments[i].set_position_cartesian(_m);
}

real CSupercell::get_fragmol_backbone_tilt(void){
  return gsf.v_fragments[__active_fragment].get_backbone_tilt();
}

real CSupercell::get_fragmol_backbone_precession(void){
  return gsf.v_fragments[__active_fragment].get_backbone_precession();
}

real CSupercell::get_fragmol_axis_tilt(void){
  return gsf.v_fragments[__active_fragment].get_axis_tilt();
}

real CSupercell::get_fragmol_axis_precession(void){
  return gsf.v_fragments[__active_fragment].get_axis_precession();
}

bool CSupercell::is_direct(void){
  return gsf.get_is_direct();
}

bool CSupercell::is_periodic(void){
  return gsf.get_is_periodic();
}

void CSupercell::is_periodic(bool b){
  gsf.set_is_periodic(b);
}

bool CSupercell::is_potmol(void){
  return __is_potmol;
}

bool CSupercell::is_fragmol_initialized(void){
  return gsf.v_fragments[__active_fragment].is_initialized();
}

std::string CSupercell::get_dir(void){
  return __sdir;
}

std::string CSupercell::get_fragment_atomic_label(uint i, uint j){
  return gsf.v_fragments[i].get_atomic_label(j);
}

std::string CSupercell::get_fragment_atomic_symbol(uint i, uint j){
  return gsf.v_fragments[i].get_atomic_symbol(j);
}

uint CSupercell::get_number_of_fragments(void){
  return gsf.get_number_of_fragments();
}

uint CSupercell::get_fragmol_active_fragment(void){
  return __active_fragment;
}

uint CSupercell::get_fragment_size(uint i){
  return gsf.v_fragments[i].size();
}

uint CSupercell::get_fragment_atomic_number(uint i,uint j){
  return gsf.v_fragments[i].get_atomic_number(j);
}

TVector<real> CSupercell::get_fragmol_axis_angles(void){
  return gsf.v_fragments[__active_fragment].get_axis_angles();
}

TVector<real> CSupercell::get_fragmol_basis_direct(void){
  return gsf.v_fragments[__active_fragment].get_basis_direct();
}

TVector<real> CSupercell::get_position_direct(void){
  return gsf.v_fragments[__active_fragment].get_origin_direct();
}

TVector<real> CSupercell::get_fragmol_position_uvw(void){
  return gsf.v_fragments[__active_fragment].get_origin_uvw();
}

TVector<real> CSupercell::get_position_cartesian(void){
  return gsf.v_fragments[__active_fragment].get_origin_cartesian();
}

TVector<real> CSupercell::get_fragment_centered_position_cartesian(void){
  return gsf.v_fragments[__active_fragment].get_centered_origin_cartesian();
}

TVector<real> CSupercell::get_fragment_direct(uint i,uint j){
  return gsf.v_fragments[i].get_direct(j);
}

TVector<real> CSupercell::get_fragment_cartesian(uint i,uint j){
  return gsf.v_fragments[i].get_cartesian(j);
  //return v_fragments[i].get_centered_cartesian(j);
}

TVector<real> CSupercell::get_fragment_centered_cartesian(uint i,uint j){
  return gsf.v_fragments[i].get_centered_cartesian(j);
}

TMatrix<real> CSupercell::get_fragment_direct(uint i){
  return gsf.v_fragments[i].get_direct();
}

TMatrix<real> CSupercell::get_fragment_cartesian(uint i){
  return gsf.v_fragments[i].get_cartesian();
}

TMatrix<real> CSupercell::get_fragment_centered_cartesian(uint i){
  return gsf.v_fragments[i].get_centered_cartesian();
}

// STRUCTURE

uint CSupercell::get_atomic_species(void){
  return gsf.get_atomic_species();
}

uint CSupercell::get_fragmol_total_atoms(void){
  return gsf.get_total_atoms();
}

uint CSupercell::get_total_atoms(void){
  return gsf.get_total_atoms();
}

uint CSupercell::get_input_file_format(void){
  return gsf.get_input_format();
}

uint CSupercell::get_output_file_type(void){
  return gsf.get_output_file_type();
}

uint CSupercell::get_output_file_format(void){
  return gsf.get_output_file_format();
}

uint CSupercell::get_atomic_composition(uint i){
  return gsf.get_atomic_composition(i);
  //v_atomic_composition_table[i];
}

TVector<real> CSupercell::get_cartesian(uint i){
  return m_xyz[i];
}

TVector<real> CSupercell::get_direct(uint i){
  return m_uvw[i];
}

real CSupercell::get_distance(uint i, uint j){
  //TVector<real> v1, v2, v3;
  TVector<real> v3;
  //v1 = get_cartesian(i);
  //v2 = get_cartesian(j);
  // v3 = v2-v1;
  v3 = get_vector_diff(i,j);
  return v3.magnitude();
}

TVector<real> CSupercell::get_vector_diff(uint i, uint j){
  TVector<real> v1, v2, v3;
  v1 = get_cartesian(i);
  v2 = get_cartesian(j);
  v3 = v2-v1;
  return v3;
}

real CSupercell::get_angle(uint i, uint j, uint k){
  TVector<real> v1, v2, v3;
  v1 = get_cartesian(i);
  v2 = get_cartesian(j);
  v3 = get_cartesian(k);
  v1 = (v2-v1);
  v3 = (v2-v3);
  real res = (v1*v3)/(v1.magnitude()*v3.magnitude());
  res = acos(res)*RAD_DEG;
  return res;
}

real CSupercell::get_dihedral(uint i, uint j, uint k, uint l){
  TVector<real> v1, v2, v3, v4, vd;
  real d12, d13, d14, d23, d24, d34;
  v1 = get_cartesian(i);
  v2 = get_cartesian(j);
  v3 = get_cartesian(k);
  v4 = get_cartesian(l);
  vd = v2-v1;
  d12 = vd.magnitude();
  vd = v3-v1;
  d13 = vd.magnitude();
  vd = v4-v1;
  d14 = vd.magnitude();
  vd = v3-v2;
  d23 = vd.magnitude();
  vd = v4-v2;
  d24 = vd.magnitude();
  vd = v4-v3;
  d34 = vd.magnitude();
  real P = SQ(d12) * ( SQ(d23)+SQ(d34)-SQ(d24)) +
            SQ(d23) * (-SQ(d23)+SQ(d34)+SQ(d24)) +
            SQ(d13) * ( SQ(d23)-SQ(d34)+SQ(d24)) -
            2 * SQ(d23) * SQ(d14);
  real Q = (d12 + d23 + d13) * ( d12 + d23 - d13) *
            (d12 - d23 + d13) * (-d12 + d23 + d13 ) *
            (d23 + d34 + d24) * ( d23 + d34 - d24 ) *
            (d23 - d34 + d24) * (-d23 + d34 + d24 );
  //v1 = (v1-v2);
  //v3 = (v4-v3);
  real res = P/sqrt(Q);
  //real res = (v1*v3)/(v1.magnitude()*v3.magnitude());
  res = acos(res)*RAD_DEG;
  return res;
}

//////
TVector<std::string> CSupercell::get_atomic_symbol_table(void){
  return gsf.get_atomic_symbol_table();
}

std::string CSupercell::get_atomic_symbol_table(uint u){
  return gsf.get_atomic_symbol_table(u);
}

TVector<uint> CSupercell::get_fragmol_atomic_composition_table(void){
  return gsf.get_atomic_composition_table();
}

TVector<uint> CSupercell::get_atomic_number_table(void){
  return gsf.get_atomic_number_table();
}

uint CSupercell::get_atomic_number_table(uint u){
  return gsf.get_atomic_number_table(u);
}

TVector<uint> CSupercell::get_atom_table(void){
  return gsf.get_atom_table();
}

uint CSupercell::get_atom_table(uint u){
  return gsf.get_atom_table(u);
}

//////

TVector<uint> CSupercell::get_fragment_table(void){
  return gsf.get_fragment_table();
}

uint CSupercell::get_fragment_table(uint u){
  return gsf.get_fragment_table(u);
}

TMatrix<real> CSupercell::get_unit_uvw_to_xyz(void){
  return gsf.get_uvw_to_xyz_u();
}

TMatrix<real> CSupercell::get_uvw_to_xyz(void){
  return gsf.get_uvw_to_xyz();
}

TVector<real> CSupercell::get_uvw_to_xyz(uint u){
  return gsf.get_uvw_to_xyz(u);
}

// deprecated
TMatrix<real> CSupercell::get_bounding_box(void){
  return gsf.get_uvw_to_xyz();
}

void CSupercell::set_inv_bbox(void){
  TMatrix<real> u_bbox;
  u_bbox=get_unit_uvw_to_xyz();
  m_inv_bbox=u_bbox.inverse();
}

TMatrix<real> CSupercell::get_cartesian(void){
  return m_xyz;
}

TMatrix<real> CSupercell::get_direct(void){
  return m_uvw;
}

// Special MD functions

void CSupercell::set_cut_radius(real r){
  r_cut_radius=r; r_cut_radius_2=r*r;
}

void CSupercell::set_cut_radius(void){
  //v_atomic_number_table_gl = v;
  //real r_cut_radius=0;
  //std::cout<<" v_atomic_number_table_gl = "<<v_atomic_number_table_gl;
  for(uint i=0; i<get_atomic_species(); i++){
    r_cut_radius=maxi(r_cut_radius,atom_rrgb[get_atomic_number_table(i)][0]);
  }
  r_cut_radius*=2.0;
  //supercell.set_cut_radius(r_cut_radius);
  r_cut_radius_2 = (r_cut_radius*r_cut_radius);
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" Cut Radius = "<<r_cut_radius<<std::endl;
  std::cout<<" Cut Radius2 = "<<r_cut_radius_2<<std::endl;
#endif
  //r_cut_radius=r; r_cut_radius_2=r*r;
};

// Linked and shell cell configuration functions
void CSupercell::set_cells(void){
  TVector<int> v1(3);
  int xpcb=-2, ypcb=-2, zpcb=-2;
  v_cell_side = iVScale(v_box_size, (1.0/r_cut_radius));
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" Cut Radius = "<<r_cut_radius<<std::endl;
  std::cout<<" LINKED: box = "<<v_box_size;
  std::cout<<" LINKED: Cell = "<<v_cell_side;
#endif
  // set the necesary neighbor cells
  if(v_cell_side[0]==2) xpcb=-1;
  if(v_cell_side[1]==2) ypcb=-1;
  if(v_cell_side[2]==2) zpcb=-1;
  neighbor_cells_xyz.resize(0,3);
  i_neighbor_cells=0;
  for(int k=0; k<27; k++){
	v1[0] = neighbor_cells[k][0];
	v1[1] = neighbor_cells[k][1];
	v1[2] = neighbor_cells[k][2];
	if((v1[0] > xpcb) && (v1[1] > ypcb) && (v1[2] > zpcb)){
	  neighbor_cells_xyz.add_row(v1);
	  i_neighbor_cells++;
	}
  }
  b_linked_cell=true;
  for(int i=0; i<3; i++){
    if(v_cell_side[1]<=2){
      b_linked_cell=false;
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" LINKED CELL too small [turned off]"<<std::endl;
#endif
    }
  }
}

void CSupercell::set_inverse_cell(void){
  v_cell_frac = fVDiv(v_cell_side,v_box_size);
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: Cell Frac = "<<v_cell_frac;
#endif
}

void CSupercell::set_cell_list(void){
  u_cell_number = (uint)vVol(v_cell_side);
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: Number of used cells = "<<u_cell_number<<std::endl;
#endif
  v_cell_head.resize(u_cell_number);
  v_cell_list.resize(get_total_atoms());
}

void CSupercell::eval_linked_list(void){
  set_cells();
  set_inverse_cell();
  set_cell_list();
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: BEGIN: eval_cell_list"<<std::endl;
#endif
  uint _n, u_icell;
  TVector<real> v_positive_r;
  TVector<int>  v_integer_r;
  for(_n=0; _n<u_cell_number; _n++)
    v_cell_head[_n] = -1;
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" LINKED: build the linked cell"<<std::endl;
#endif
  for(_n=0; _n<(uint)get_total_atoms(); _n++){
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" atom coordinates["<<_n<<"] = "<<get_cartesian(_n); //m_atom_coordinates[_n];
#endif
    v_positive_r=get_cartesian(_n); //m_atom_coordinates[_n];
    // uvw coordinates
    v_positive_r =  (v_positive_r*m_inv_bbox);
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" direct r = "<<v_positive_r;
#endif
    // use half of the box here
    v_positive_r = fVAdd(v_positive_r,v_bbox);
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" positive r = "<<v_positive_r;
#endif
    // apply PBC to place all the atoms inside the box
    v_integer_r = iVMul(v_positive_r,v_cell_frac);
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" v_integer_r = "<<v_integer_r;
#endif
    u_icell = iVLinear(v_integer_r,v_cell_side);
#ifdef _SHOW_DEBUG_LINKED_EVAL_
    std::cout<<" u_icell = "<<u_icell<<std::endl;
#endif
    v_cell_list[_n] = v_cell_head[u_icell];
    v_cell_head[u_icell] = _n;
  }
#ifdef _SHOW_DEBUG_LINKED_
  std::cout<<" v_cell_list="<<v_cell_list<<std::endl;
  std::cout<<" END: eval_cell_list"<<std::endl;
#endif
}

uint CSupercell::check_bond(uint u){
  uint size = v_bond_table.size();
  for(uint i=0; i<size; i++){
    if(v_bond_table[i]==u){
      return i;
    }
  }
  v_bond_table.push_back(u);
  return size;
}

// Fri Jan 13 16:55:51 MST 2012
// beta version
// find bonds between atoms closer than the sum of their van der Waals radius
void CSupercell::eval_atomic_bonds(void){
#ifdef _SUPERCELL_DEBUG_MESSAGES_
  std::cout<<" SUPERCELL: eval_atomic_bonds "<<std::endl;
#endif
#ifdef _SHOW_TIME_
  gl_atom_clock.start();
#endif
  int j;
  bool use_pbc;
  real rl, r, r2, rr, ri, rj, rlz, rlxy;
  TVector<uint> vidx(2);
  TVector<real> vi, vj, vij, vang(2);
  TVector<real> vi_uvw, vj_uvw, vij_uvw;
  TVector<uint> v_ft;
  v_ft = get_fragment_table();
  //
  int u_icell;
  uint k;
  TVector<real> v_positive_r;
  TVector<int>  v_integer_r;
  TVector<int>  v_neighbor_cell;
  TVector<int>  v_pbc;
  //
  uint max_bonds = 15*get_total_atoms();
#ifdef _SUPERCELL_DEBUG_MESSAGES_
  std::cout<<" SUPERCELL: estimated max bonds "<<max_bonds<<std::endl;
#endif
  v_bond_table.resize(0);
  i_number_of_bonds = 0;
  i_number_of_bonds_pbc = 0;
  ////////////////////////////////////////
  v_bond_number.resize(max_bonds);
  m_bond_indices.resize(max_bonds,2);
  v_bond_number_pbc.resize(uint(max_bonds/2));
  m_bond_indices_pbc.resize(uint(max_bonds/2),2);
  m_bond_boundary_pbc.resize(uint(max_bonds/2),2);
  ////////////////////////////////////////
  v_pbc.resize(3);
//#ifdef _SUPERCELL_DEBUG_BONDS_
  //std::cout<<" FL_GL_SUPERCELL: LINKED: Box = "<<v_box_size;
//#endif
  if(is_linked_cell()){
#ifdef _SUPERCELL_DEBUG_BONDS_
    std::cout<<" FL_GL_SUPERCELL: LINKED CELL USED"<<std::endl;
#endif
    eval_linked_list();
#ifdef _SUPERCELL_DEBUG_BONDS_
    std::cout<<" FL_GL_SUPERCELL: LINKED CELL READY"<<std::endl;
#endif
  }
#ifdef _SUPERCELL_DEBUG_BONDS_
  else{
    std::cout<<" FL_GL_SUPERCELL: LINKED CELL NOT USED"<<std::endl;
  }
#endif
  for(int i=0; i<get_total_atoms()-1; i++){
#ifdef _SUPERCELL_DEBUG_BONDS_
    std::cout<<" SUPERCELL: i="<<i<<std::endl;
#endif
    if(strcmp(get_atomic_symbol_table(get_atom_table(i)).c_str(),"X")){
      vi = get_cartesian(i);//m_atom_coordinates[i];
      vi_uvw = (vi*get_inv_bbox());
      ri = m_radius_color[i][0];
#ifdef _SUPERCELL_DEBUG_BONDS_
    std::cout<<" SUPERCELL: ii="<<vi<<std::endl;
    std::cout<<" SUPERCELL: vi_uvw="<<vi_uvw<<std::endl;
    std::cout<<" SUPERCELL: ii="<<vi<<std::endl;
#endif
      if(is_linked_cell()){
        v_positive_r = fVAdd(vi_uvw,get_bbox());
        v_integer_r  = iVMul(v_positive_r,get_cell_frac());
#ifdef _SUPERCELL_DEBUG_BONDS_
        std::cout<<" SUPERCELL: v_positive_r="<<v_positive_r<<std::endl;
        std::cout<<" SUPERCELL: v_integer_r="<<v_integer_r<<std::endl;
#endif
        // using the neighbour list
        // searching inside the neighbour cells
        for (int _m=0; _m<get_neighbor_cells(); _m++){
          v_neighbor_cell = iVAdd(v_integer_r,get_neighbor_cells_xyz(_m));
#ifdef _SUPERCELL_DEBUG_BONDS_
          std::cout<<" SUPERCELL: v_neighbor_cell="<<v_neighbor_cell<<std::endl;
#endif
          // Used to apply PCB to cells in each dimension
          /////////////////////////////////////////////////////////////////////////////////
          for (uint coord=0; coord<3; coord++){
            if(v_neighbor_cell[coord] >= v_cell_side[coord]){ // check if  PBC is necessary
              v_neighbor_cell[coord]= 0;                      // apply PBC to each  cell
              //if(is_pbc) use_pbc = false;
            }else if(v_neighbor_cell[coord] < 0){             // check if  PBC is necessary
              v_neighbor_cell[coord]= v_cell_side[coord]-1;   // apply PBC to each cell
              //if(is_pbc) use_pbc = false;
            }
          }
          /////////////////////////////////////////////////////////////////////////////////
          u_icell = iVLinear(v_neighbor_cell,v_cell_side);    // head atom index in the the cell
#ifdef _SUPERCELL_DEBUG_BONDS_
          std::cout<<" SUPERCELL: v_neighbor_cell="<<v_neighbor_cell<<std::endl;
          std::cout<<" SUPERCELL: u_icell ="<<u_icell<<std::endl;
#endif
          if(u_icell>=0 && u_icell < (int)get_cell_number()){          // inside of a cells
            j = get_cell_head(u_icell);                         // head atom in the actual cell
          }else{                                              // out of the box
            j = -1;                                           // outside of a cells
          }
          while(1){                                           // over all the particles in the cell
            if(j<0) break;                                    // stop searching in the cell
            //if(_m!=0 || j>i){                               // avoid self-interaction
            if((j>i) && (v_ft[j] == v_ft[i])){                // avoid self-interaction and double bond
            //if(j>i){                                        // avoid self-interaction and double bond
#ifdef _SUPERCELL_DEBUG_BONDS_
              std::cout<<" FL_GL_SUPERCELL: test bond for [i-j]="<<i<<","<<j<<std::endl;
#endif
              r2 = 0;                                         // set distance to cero
              vj = get_cartesian(j);  //m_atom_coordinates[j];
              vj_uvw = vj*get_inv_bbox();
              rj = m_radius_color[j][0];
              r = (ri+rj);
              rr = 1.2*(r*r);
              use_pbc = false;
              vij_uvw = (vj_uvw-vi_uvw);
              for(uint coord=0; coord<3; coord++){
                v_pbc[coord] = 0;
                if(vij_uvw[coord] <= -get_bbox(coord)){
                  vj += get_uvw_to_xyz(coord); //2.0*m_bbox[coord];       // PBC
                  use_pbc = true;
                  v_pbc[coord] = 1;
                }else if(vij_uvw[coord] > get_bbox(coord)){
                  vj -= get_uvw_to_xyz(coord); //2.0*m_bbox[coord];       // PBC
                  use_pbc = true;
                  v_pbc[coord] = -1;
                }
              }
              vij = (vj-vi);
              r2 = vij.norm();
              if(r2 < rr){                                    // atoms inside de cut radius
#ifdef _SUPERCELL_DEBUG_BONDS_
                std::cout<<" FL_GL_SUPERCELL: Bond found [i-j]="<<i<<","<<j<<std::endl;
#endif
                vidx[0]=i;
                vidx[1]=j;
                r = sqrt(r2);
                vij = (vj-vi);
                rlz = vij[2];
                vij[2] = 0;
                rlxy = vij.magnitude();
                vang[0] = RAD_DEG*atan2(vij[1],vij[0]);       // bond precession
                vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz));      // bond tilt
                k = get_bond_index(r);
                if(use_pbc){
                  v_bond_number_pbc[i_number_of_bonds_pbc]=check_bond(k);
                  m_bond_indices_pbc[i_number_of_bonds_pbc]=vidx;
                  m_bond_boundary_pbc[i_number_of_bonds_pbc]=v_pbc;
                  i_number_of_bonds_pbc++;
                }else{
                  v_bond_number[i_number_of_bonds]=check_bond(k);
                  m_bond_indices[i_number_of_bonds]=vidx;
                  i_number_of_bonds++;
                }
              }
            }
            j = get_cell_list(j);
          }
        }
      }else{
#ifdef _SUPERCELL_DEBUG_BONDS_
        std::cout<<" FL_GL_SUPERCELL: NOT LINKED CELL USED"<<std::endl;
#endif
        // the code below can be used for small number of atoms.
        // searching inside the neighbour cells
        for (int _m=0; _m<27; _m++){
          for(int j=i+1; j<get_total_atoms(); j++){
            if(v_ft[j] == v_ft[i]){  // avoid bonds between fragments
#ifdef _SUPERCELL_DEBUG_BONDS_
            std::cout<<" SUPERCELL: j="<<j<<" i="<<i<<" m="<<_m<<std::endl;
#endif
            if(strcmp(get_atomic_symbol_table(get_atom_table(j)).c_str(),"X")){
              vj = get_cartesian(j); //m_atom_coordinates[j];
              for(uint coord=0; coord<3; coord++){
                vj += neighbor_cells[_m][coord]*get_uvw_to_xyz(coord); //  2.0*m_bbox[coord];       // PBC
                v_pbc[coord] = neighbor_cells[_m][coord];
              }
              rj = m_radius_color[j][0];
              vj_uvw = vj*get_inv_bbox();
              vij = (vj-vi);
              r = (ri+rj);
              rr = 1.2*(r*r);
              rl = vij.norm();
              if( (rl <= rr) &&  (rl > 0.1) ){
#ifdef _SUPERCELL_DEBUG_BONDS_
                std::cout<<" SUPERCELL: j="<<j<<" i="<<i<<" bonded"<<std::endl;
#endif
                vidx[0]=i;
                vidx[1]=j;
                rlz = vij[2];
                vij[2] = 0;
                rlxy = vij.magnitude();
                vang[0] = RAD_DEG*atan2(vij[1],vij[0]);  // precession
                vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz)); // tilt
                r = sqrt(rl);
                uint k = get_bond_index(r);
                if(_m>0){
                  v_bond_number_pbc[i_number_of_bonds_pbc]=check_bond(k);
                  m_bond_indices_pbc[i_number_of_bonds_pbc]=vidx;
                  m_bond_boundary_pbc[i_number_of_bonds_pbc]=v_pbc;
                  i_number_of_bonds_pbc++;
#ifdef _SUPERCELL_DEBUG_BONDS_
                  std::cout<<" SUPERCELL: bond pcb type: "<<check_bond(k)<<std::endl;
#endif
                }else{
                  v_bond_number[i_number_of_bonds]=check_bond(k);
                  m_bond_indices[i_number_of_bonds]=vidx;
                  i_number_of_bonds++;
#ifdef _SHOW_DEBUG_NOPBC_BONDS_
                  std::cout<<" SUPERCELL: bond type: "<<check_bond(k)<<std::endl;
#endif
                }
              }
            }
            }
          }
        }
      }
    }
  }
  //m_bond_rcolor_0.resize(i_number_of_bonds,4);
  //m_bond_rcolor_1.resize(i_number_of_bonds,4);
  //m_bond_rcolor_pbc_0.resize(i_number_of_bonds_pbc,4);
  //m_bond_rcolor_pbc_1.resize(i_number_of_bonds_pbc,4);
  //
  //m_bond_angles.resize(i_number_of_bonds,2);
  //m_bond_angles_pbc.resize(i_number_of_bonds_pbc,2);
  //m_bond_position.resize(i_number_of_bonds,3);
  //m_bond_position_pbc.resize(i_number_of_bonds_pbc,3);
#ifdef _SUPERCELL_DEBUG_BONDS_
  std::cout<<" SUPERCELL: i_number_of_bonds="<<i_number_of_bonds<<std::endl;
  std::cout<<" SUPERCELL: i_number_of_bonds_pbc="<<i_number_of_bonds_pbc<<std::endl;
  std::cout<<" SUPERCELL: v_bond_number="<<v_bond_number;
  std::cout<<" SUPERCELL: v_bond_number_pbc="<<v_bond_number_pbc;
  std::cout<<" SUPERCELL: v_bond_table="<<v_bond_table;
  //std::cout<<" SUPERCELL: m_bond_position_pbc="<<m_bond_position_pbc;
#endif
  //u_bond_types=v_bond_table.size();
#ifdef _SUPERCELL_DEBUG_MESSAGES_
  std::cout<<" SUPERCELL: max number of bonds "<<max_bonds<<std::endl;
  std::cout<<" SUPERCELL: number of bonds "<<i_number_of_bonds<<std::endl;
  std::cout<<" SUPERCELL: number of PBC bonds "<<i_number_of_bonds_pbc<<std::endl;
  std::cout<<" SUPERCELL: eval_atomic_bonds "<<std::endl;
#endif
  //update_atomic_bonds();
#ifdef _SHOW_TIME_
  gl_atom_clock.stop();
  gl_atom_clock.show();
#endif
  // send update color flag
  //update_bonds_color=true;
}

// SUPERCELL END
