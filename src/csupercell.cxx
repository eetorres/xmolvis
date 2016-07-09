//======================================================================//
// FILE: csupercell.cxx -> csupercell                                   //
//                                                                      //
// Abstraction layer of the supercell and fragments.                    //
// The manipulation of individual fragments is hidden.                  //
//                                                                      //
// Copyright 2011-2015 by Edmanuel Torres                               //
// email:   eetorres@gmail.com                                          //
//                                                                      //
// Comment: change this class/file names to CSupercell/csupercell       //
//                                                                      //
// Lastest update: Tue Jul 14 16:22:34 EDT 2015                         //
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
}

CSupercell::~CSupercell(void){
}

void CSupercell::clear(void){
  __is_potmol=false;
  __active_fragment=0;
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
    v_atomic_composition_table=gsf.get_atomic_composition_table();
    v_atomic_number_table=gsf.get_atomic_number_table();
    v_atomic_symbol_table=gsf.get_atomic_symbol_table();
    v_atom_type_table=gsf.get_atom_table();
    v_atomic_labels=gsf.get_atomic_labels();
    v_atomic_symbols=gsf.get_atomic_symbols();
    v_atomic_numbers=gsf.get_atomic_numbers();
    v_atom_cell_table.resize(gsf.get_total_atoms());
#ifdef _FRAGMOL_DEBUG_MESSAGES_
    std::cout<<" FRAGMOL: Composition ="<<v_atomic_composition_table;
    std::cout<<" FRAGMOL: Z-Number Table ="<<v_atomic_number_table;
    std::cout<<" FRAGMOL: Symbol Table ="<<v_atomic_symbol_table;
    std::cout<<" FRAGMOL: Type ="<<v_atom_type_table;
    std::cout<<" FRAGMOL: Labels ="<<v_atomic_labels;
    std::cout<<" FRAGMOL: Symbols ="<<v_atomic_symbols;
    std::cout<<" FRAGMOL: Z-Numbers ="<<v_atomic_numbers;
#endif
    //
    m_xyz=gsf.get_xyz();
    m_uvw=gsf.get_uvw();
#ifdef _FRAGMOL_DATA_MESSAGES_
    std::cout<<" FRAGMOL: XYZ = "<<m_xyz;
    std::cout<<" FRAGMOL: UVW = "<<m_uvw;
    std::cout<<" FRAGMOL: uvwTxyz U = "<<m_uvw_to_xyz_u;
    std::cout<<" FRAGMOL: uvwTxyz  = "<<m_uvw_to_xyz;
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
  std::cout<<" FRAGCAR:  ["<<get_fragmol_number_of_fragments()<<"] fragments loaded"<<std::endl;
#endif
}

void CSupercell::eval_cell_table(void){
  TVector<uint> v_l;
  uint _s;
  for(uint i=0;i<gsf.get_number_of_fragments();i++){
    _s=v_fragments[i].size();
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

void CSupercell::eval_connections(const TMatrix<uint>& _m, uint u){
  gsf.eval_connections(_m,u);
}

void CSupercell::eval_fragmol_initial_position(void){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: eval_fragmol_initial_position"<<std::endl;
#endif
  gsf.v_fragments[__active_fragment].eval_initial_position();
}

void CSupercell::eval_fragmol_initial_orientation(void){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: eval_fragmol_initial_orientation"<<std::endl;
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

void CSupercell::compute_fragmol_position_cartesian(void){
  gsf.v_fragments[__active_fragment].compute_origin_cartesian();
}

void CSupercell::is_fragmol_initialized(bool b){
  gsf.v_fragments[__active_fragment].is_initialized(b);
}

void CSupercell::apply_pbc(bool b){
  gsf.v_fragments[__active_fragment].is_pbc(b);
}

void CSupercell::set_gsf_modified(const bool b){
#ifdef _FRAGMOL_DEBUG_MESSAGES_
  std::cout<<" FRAGMOL: the gsf was modified"<<std::endl;
#endif
  gsf.set_modified(b);
}

void CSupercell::set_fragmol_active_fragment(const uint i){
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
  gsf.v_fragments[__active_fragment].set_axis_index(idx,v_atom_cell_table[val]);
}

void CSupercell::set_fragmol_fragment_twist(const real r){
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

uint CSupercell::get_fragmol_number_of_fragments(void){
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

uint CSupercell::get_fragmol_total_atoms(void){
  return gsf.get_total_atoms();
}

TVector<real> CSupercell::get_fragmol_axis_angles(void){
  return gsf.v_fragments[__active_fragment].get_axis_angles();
}

TVector<real> CSupercell::get_fragmol_basis_direct(void){
  return gsf.v_fragments[__active_fragment].get_basis_direct();
}

TVector<real> CSupercell::get_fragmol_position_direct(void){
  return gsf.v_fragments[__active_fragment].get_origin_direct();
}

TVector<real> CSupercell::get_fragmol_position_uvw(void){
  return gsf.v_fragments[__active_fragment].get_origin_uvw();
}

TVector<real> CSupercell::get_fragmol_position_cartesian(void){
  return gsf.v_fragments[__active_fragment].get_origin_cartesian();
}

TVector<real> CSupercell::get_fragmol_centered_position_cartesian(void){
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

uint CSupercell::get_species(void){
  return gsf.get_atomic_species();
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
  return v_atomic_composition_table[i];
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
TVector<std::string> CSupercell::get_fragmol_atomic_symbol_table(void){
  return gsf.get_atomic_symbol_table();
}

TVector<uint> CSupercell::get_fragmol_atomic_composition_table(void){
  return gsf.get_atomic_composition_table();
}

TVector<uint> CSupercell::get_fragmol_atomic_number_table(void){
  return gsf.get_atomic_number_table();
}

TVector<uint> CSupercell::get_fragmol_atom_table(void){
  return gsf.get_atom_table();
}
//////

TVector<uint> CSupercell::get_fragmol_fragment_table(void){
  return gsf.get_fragment_table();
}

TMatrix<real> CSupercell::get_unit_uvw_to_xyz(void){
  return gsf.get_uvw_to_xyz_u();
}

TMatrix<real> CSupercell::get_uvw_to_xyz(void){
  return gsf.get_uvw_to_xyz();
}

TMatrix<real> CSupercell::get_bounding_box(void){
  return gsf.get_uvw_to_xyz();
}

TMatrix<real> CSupercell::get_cartesian(void){
  return m_xyz;
}

TMatrix<real> CSupercell::get_direct(void){
  return m_uvw;
}

// SUPERCELL END
