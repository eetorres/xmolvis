//========================================================================
// FILE - fl_gl_atom.cxx                                              //
// For the Fast Light Tool Kit (FLTK) - www.fltk.org                    //
//========================================================================
//                                                                      //
// OpenGL 3D visualization widget to display atomic structures.         //
//                                                                      //
// Copyright 2002-2015 by Edmanuel Torres                               //
// email: eetorres@gmail.com                                            //
//                                                                      //
// Lastest update: 13/07/2015
//======================================================================//

#include <config_debug.h>
#include <fl_gl_atom.h>

#if HAVE_GL
Fl_Gl_Atom::Fl_Gl_Atom(int x,int y,int w,int h,const char *l) : Fl_Gl_Window(x,y,w,h,l)
#else
Fl_Gl_Atom::Fl_Gl_Atom(int x,int y,int w,int h,const char *l) : Fl_Box(x,y,w,h,l)
#endif // HAVE_GL
{
  cell.pos_x_cells = 0;
  cell.pos_y_cells = 0;
  cell.pos_z_cells = 0;
  cell.neg_x_cells = 0;
  cell.neg_y_cells = 0;
  cell.neg_z_cells = 0;
  param.u_sphere_resolution = 0;
  param.u_cylinder_resolution = 10;
  is_eval_sphere=true;
  cell.total_cells = 1;
  cell.x_cells     = 1;
  cell.y_cells     = 1;
  cell.z_cells     = 1;
  is_update_bonds  = false;
#if !HAVE_GL
  label("OpenGL is required for this demo to operate.");
  align(FL_ALIGN_WRAP | FL_ALIGN_INSIDE);
#endif // !HAVE_GL
}

int Fl_Gl_Atom::handle(int event){
  int ret = Fl_Gl_Window::handle(event);
  /*switch(event){
  case FL_FOCUS:
    //std::cout<<" got the focus"<<std::endl;
    //if(Fl::focus() != this)
    //    Fl::focus(this);
    //return 1;
  case FL_UNFOCUS:
    //if(Fl::focus() != this)
    //    Fl::focus(this);
    //... Return 1 if you want keyboard events, 0 otherwise
    //redraw();
    return 1;
  }*/
  return ret;
}

void Fl_Gl_Atom::eval_initial_properties(void){
  real _radius;
  TVector<real> _atom_xyz;
  TVector<real> _rcolor(4);
  index_palette.set_color(0);
  index_palette.initialize(0,get_total_atoms()+MENU_RESERVED_IDS,get_total_atoms()+MENU_RESERVED_IDS);
  index_palette.update_palette_index();
  int i_z=0;
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: Number of atoms: "<<get_total_atoms()<<std::endl;
#endif
  if(is_eval_sphere){
    set_sphere_resolution(2);
    is_eval_sphere=false;
  }
  if(is_draw_atoms_){
    supercell.set_radius_color();
  }
  m_atom_rcolor = supercell.get_radius_color();
}

//void Fl_Gl_Atom::set_axis_position(const TVector<real>& v){
//  fragment.v_axis_position=v;
//}

void Fl_Gl_Atom::update_atomic_coordinates(void){
  is_update_bonds = true;
  is_update_atomic_properties = true;
  is_update_mask_rcolor = true;
  set_xyz_cells();
}

void Fl_Gl_Atom::initialize_atomic_coordinates(void){
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: Initialize"<<std::endl;
#endif
  update_atomic_coordinates();
  if(get_total_atoms()<100){
    supercell.is_linked_cell(false);
  }
  if(is_first_structure_){
#ifdef _ATOM_DEBUG_MESSAGES_
    std::cout<<" ATOM: First structures"<<std::endl;
#endif
    if(is_initialize_rot)
      initialize_rotation_matrix();
    initialize_transform_matrix();
    is_draw_atoms_=true;
    is_first_structure_=false;
  }
  is_update_atomic_properties = true;
  if(is_draw_bonds_){
    is_update_bonds=true;
  }
  is_update_radius=true;
  eval_initial_properties(); // <---------------------!!!@
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: Initialized"<<std::endl;
#endif
}

// this the function in charge of the scene actualization
void Fl_Gl_Atom::update_data(void){
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: begin update_data"<<std::endl;
#endif
  update_view();
  if(is_update_coordinates){
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: update_atomic_coordinates "<<std::endl;
#endif
    update_atomic_coordinates();
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: set update "<<std::endl;
#endif
    fragment.set_axis_position(get_fragment_centered_position_cartesian());
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: set update 1"<<std::endl;
#endif
    fragment.set_axis_precession(get_axis_precession());
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: set update 2"<<std::endl;
#endif
    fragment.set_axis_tilt(get_axis_tilt());
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: set update 3"<<std::endl;
#endif
    fragment.set_backbone_precession(get_backbone_precession());
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: set update 4"<<std::endl;
#endif
    fragment.set_backbone_tilt(get_backbone_tilt());
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: set update 5"<<std::endl;
#endif
    is_update_coordinates=false;
  }
  Fl::wait(0.1);
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: end update_data"<<std::endl;
#endif
}

void Fl_Gl_Atom::update_atomic_bonds(void){
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: start - update_atomic_bonds(void)"<<std::endl;
#endif
  uint i, j;
  real rlz, rlxy;
  TVector<real> vi, vj, vij, vang(2);
  TVector<real> vi_uvw, vj_uvw, vij_uvw;
  TVector<int>  v_pbc;
  // IN CELL BONDS
  for(uint n=0; n<supercell.get_number_of_bonds(); n++){
    i=supercell.get_bond_indices(n,0); //m_bond_indices[n][0];
    j=supercell.get_bond_indices(n,1);  //m_bond_indices[n][1];
    vi = get_cartesian(i);
    vj = get_cartesian(j);
    vij = (vj-vi);
    rlz = vij[2];
    vij[2] = 0;
    rlxy = vij.magnitude();
    vang[0] = RAD_DEG*atan2(vij[1],vij[0]);       // bond precession
    vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz));      // bond tilt
    vij = 0.5*(vj+vi);
    m_bond_angles[n]=vang;
    m_bond_position[n]=vij;
  }
  // PBC BONDS
  for(uint n=0; n<supercell.get_number_of_bonds_pbc(); n++){
    i=supercell.get_bond_indices_pbc(n,0);  //m_bond_indices_pbc[n][0];
    j=supercell.get_bond_indices_pbc(n,1);  //m_bond_indices_pbc[n][1];
    vi = get_cartesian(i);
    vi_uvw = (vi*supercell.get_inv_bbox());
    vj = get_cartesian(j);
    for(uint coord=0; coord<3; coord++){
      if ( supercell.get_bond_boundary_pbc(n,coord) != 0 )
        vj += supercell.get_bond_boundary_pbc(n,coord)*supercell.get_uvw_to_xyz(coord); //2.0*m_bbox[coord];       // PBC
    }
    vj_uvw = (vj*supercell.get_inv_bbox());
    vij = (vj-vi);
    rlz = vij[2];
    vij[2] = 0;
    rlxy = vij.magnitude();
    vang[0] = RAD_DEG*atan2(vij[1],vij[0]);       // bond precession
    vang[1] = RAD_DEG*fabs(atan2(rlxy,rlz));      // bond tilt
    vij = 0.5*(vj+vi);
    m_bond_angles_pbc[n]=vang;
    m_bond_position_pbc[n]=vij;
  }
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: end - update_atomic_bonds(void)"<<std::endl;
#endif
}

void Fl_Gl_Atom::set_x_cells(int i){
  cell.set_x_cells(i);
  set_xyz_cells();
  supercell.set_gsf_modified(true);
}

void Fl_Gl_Atom::set_y_cells(int i){
  cell.set_y_cells(i);
  set_xyz_cells();
  supercell.set_gsf_modified(true);
}

void Fl_Gl_Atom::set_z_cells(int i){
  cell.set_z_cells(i);
  set_xyz_cells();
  supercell.set_gsf_modified(true);
}



void Fl_Gl_Atom::set_xyz_cells(void){
  TVector<real> _x(3);
  TVector<real> _xyz;
  TVector<real> e(3),p(3);
  uint cont = 0;
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" total cells = "<<cell.get_total_cells()<<std::endl;
#endif
  m_atom_position.resize(cell.get_total_cells()*get_total_atoms(),3);
  if(is_draw_atoms_){
    for(int x=cell.neg_x_cells; x<cell.pos_x_cells+1; x++){ // repetition in x
      for(int y=cell.neg_y_cells; y<cell.pos_y_cells+1; y++){ // repetition in y
        for(int z=cell.neg_z_cells; z<cell.pos_z_cells+1; z++){ // repetition in z
          for(int i=0; i<get_total_atoms(); i++){
            _xyz=get_cartesian(i); //m_atom_coordinates[i];
            _xyz=_xyz+2.0*(x*_vu+y*_vv+z*_vw);
            m_atom_position[cont]=_xyz;
            cont++;
          }
        }
      }
    }
  }
}



void Fl_Gl_Atom::save_wysiwyg_as(std::string _p, std::string _f){
  // in case the bonds have not been evaluated
  if(((supercell.get_output_file_format() == OUTPUT_FORMAT_ATM_NFR) || (supercell.get_output_file_format() == OUTPUT_FORMAT_ATM_FRG)) && (supercell.get_output_file_type() == OUTPUT_FILE_TYPE_DLP)){
    //std::cout<<"bonds have not been evaluated"<<std::endl;
    if(is_eval_bonds){
      eval_atomic_bonds();  // <----------
      is_eval_bonds=false;
    }
    supercell.eval_connections(supercell.get_number_of_bonds());
  }
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: do save_wysiwyg_as (1)"<<std::endl;
  std::cout<<" ATOM: m_atom_position = "<<m_atom_position;
#endif
  supercell.save_as_file(_p,_f,is_draw_symbols_,is_draw_numbers_,m_atom_position,cell.x_cells,cell.y_cells,cell.z_cells,cell.total_cells,is_draw_bbox_);
}

void Fl_Gl_Atom::save_wysiwyg_as(std::string _f){
  //if((get_output_file_format() == OUTPUT_FORMAT_ATM_NFR) || (get_output_file_format() == OUTPUT_FORMAT_ATM_FRG))
    //std::cout<<"check";
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: do save_wysiwyg_as (2)"<<std::endl;
#endif
  supercell.save_as_file(_f,is_draw_symbols_,is_draw_numbers_,m_atom_position,cell.x_cells,cell.y_cells,cell.z_cells,cell.total_cells,is_draw_bbox_);
}

void Fl_Gl_Atom::save_wysiwyg_extension(std::string _f){
  std::string _fext;
  switch(supercell.get_output_file_type()){
    case OUTPUT_FILE_TYPE_VSP:
      _fext = _f + ".vsp";
    break;
    case OUTPUT_FILE_TYPE_XYZ:
      _fext = _f + ".xyz";
    break;
    case OUTPUT_FILE_TYPE_GAU:
      _fext = _f + ".gau";
    break;
    case OUTPUT_FILE_TYPE_PDB:  // reserved for xmatrix
      _fext = _f + ".pdb";
    break;
    case OUTPUT_FILE_TYPE_DLP:
      _fext = _f + ".dlp";
    break;
  }
  supercell.save_as_file(_fext,is_draw_symbols_,is_draw_numbers_,m_atom_position,cell.x_cells,cell.y_cells,cell.z_cells,cell.total_cells,is_draw_bbox_);
}

real Fl_Gl_Atom::set_bounding_box(real r){
  real r_view=r;
  TVector<real> v_bb;
  v_bb.resize(3);
  _vu = 0.5*supercell.get_uvw_to_xyz(0); //m_bbox[0];
  _vv = 0.5*supercell.get_uvw_to_xyz(1); //m_bbox[1];
  _vw = 0.5*supercell.get_uvw_to_xyz(2); //m_bbox[2];
  // set up the half box sides
  v_bb[0]=_vu.magnitude();
  v_bb[1]=_vv.magnitude();
  v_bb[2]=_vw.magnitude();
  supercell.set_bbox(v_bb);
  // the full box is used for the bond search method
  supercell.set_box_size(v_bb);
  supercell.set_inv_bbox();
  // set view size
  r_view = maxi(_vu.magnitude(), _vv.magnitude());
  r_view = maxi(r_view, _vw.magnitude());
  //r_axes_position = 0.9*r_view;
  r_view *= 1.1;
  return r_view;
}

void Fl_Gl_Atom::set_sphere_resolution(uint u){
  if(param.u_sphere_resolution != u){
    param.u_sphere_resolution = u;
    eval_sphere(param.u_sphere_resolution);
    param.u_cylinder_resolution=(5*param.u_sphere_resolution+10);
    eval_cylinder(param.u_cylinder_resolution);
  }
}

void Fl_Gl_Atom::initialize_sphere(real r){
    int s;
    real scale = r;
    // iterate over the 20 sides of the icosahedron
    int cont = 0;
    for(s = 0; s < 20; s++){
        for(int i = 0; i < param.u_sphere_rows; i++){
            // strip the ith trapezoid block
            glBegin(GL_TRIANGLE_STRIP);
            for(int j=0; j<(2*i+3); j++){
                // add one vertice at a time
                GLV(m_sphere[cont]);
                cont++;
            }
            glEnd();
        }
    }
}

void Fl_Gl_Atom::initialize_cylinder(real r,real d){
  TVector<real> t(3);
  glBegin(GL_QUAD_STRIP);
  real scale = 0.2*d;
  for (int i=0;i<param.u_cylinder_strip_length;i++){
    t = m_cylinder[i];
    t[2]=0; // -r
    glNormal3f(t[0],t[1],t[2]);
    glVertex3f(scale*t[0],scale*t[1],t[2]);
    t[2]=r;
    glNormal3f(t[0],t[1],t[2]);
    glVertex3f(scale*t[0],scale*t[1],t[2]);
  }
  glEnd();
}

void Fl_Gl_Atom::add_stick(const TVector<real>& c,real l, real r, real a1, real a2){
  TVector<real> e(3), p(3), t(3);
#ifdef _SHOW_DEBUG_ADD_STICK_
  std::cout<<" ATOM: add_stick"<<std::endl;
#endif
  glBegin(GL_QUAD_STRIP);
  for (int i=0;i<=param.u_cylinder_resolution;i++) {
    e = m_cylinder_e1[i];
    p[0] = r * e[0];
    p[1] = r * e[1];
    p[2] = 0;
    t = vrot_y(p,a2);
    t = vrot_z(t,a1);
    t = c+t;
    e = vrot_y(e,a2);
    e = vrot_z(e,a1);
    glNormal3f(e[0],e[1],e[2]);
    glVertex3f(t[0],t[1],t[2]);
    p[2] = l;
    t = vrot_y(p,a2);
    t = vrot_z(t,a1);
    t = c+t;
    glNormal3f(e[0],e[1],e[2]);
    glVertex3f(t[0],t[1],t[2]);
  }
#ifdef _SHOW_DEBUG_ADD_STICK_
  std::cout<<" ATOM: add_stick end"<<std::endl;
#endif
  glEnd();
}

void Fl_Gl_Atom::add_axis(const TVector<real>& c,real l, real r, real a1, real a2){
#ifdef _SHOW_DEBUG_ADD_AXIS_
  std::cout<<" ATOM: add_axis"<<std::endl;
#endif
  TVector<real> e(3),p(3), t1(2), t2(2);
#ifdef _SHOW_DEBUG_ADD_AXIS_
  std::cout<<" ATOM: check r"<<std::endl;
#endif
  if(r < 0) r = -r;
  if(r <= 0) {
    glBegin(GL_POINTS);
    glVertex3f(c[0],c[1],c[2]);
    glEnd();
    return;
  }
#ifdef _SHOW_DEBUG_ADD_AXIS_
  std::cout<<" ATOM: init"<<std::endl;
#endif
  add_stick(c,l,r,a1,a2);
#ifdef _SHOW_DEBUG_ADD_AXIS_
  std::cout<<" ATOM: strip"<<std::endl;
#endif
  glBegin(GL_QUAD_STRIP);
  real tip = l/3.0;
  for (int i=0;i<=param.u_cylinder_resolution;i++) {
    e = m_cylinder_e1[i];
    p[0] = 2.5*r* e[0];
    p[1] = 2.5*r* e[1];
    p[2] = l;
    p = vrot_y(p,a2);
    p = vrot_z(p,a1);
    p = c+p;
    glNormal3f(e[0],e[1],0);
    glTexCoord2f(m_cylinder_texture1[i][0],m_cylinder_texture1[i][1]);
    glVertex3f(p[0],p[1],p[2]);
    p[0] = 0;
    p[1] = 0;
    p[2] = l+tip;
    p = vrot_y(p,a2);
    p = vrot_z(p,a1);
    p = c+p;
    glNormal3f(e[0],e[1],0);
    glTexCoord2f(m_cylinder_texture1[i][0],m_cylinder_texture1[i][1]);
    glVertex3f(p[0],p[1],p[2]);
  }
  glEnd();
}

void Fl_Gl_Atom::eval_atomic_bonds(void){
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: eval_atomic_bonds "<<std::endl;
#endif
#ifdef _SHOW_TIME_
  gl_atom_clock.start();
#endif
  TVector<uint> v_ft;
  v_ft = supercell.get_fragment_table();
  uint max_bonds = 15*get_total_atoms();
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: estimated max bonds "<<max_bonds<<std::endl;
#endif
  ////////////////////////////////////////
  supercell.eval_atomic_bonds();
  //
  m_bond_rcolor_0.resize(supercell.get_number_of_bonds(),4);
  m_bond_rcolor_1.resize(supercell.get_number_of_bonds(),4);
  m_bond_rcolor_pbc_0.resize(supercell.get_number_of_bonds_pbc(),4);
  m_bond_rcolor_pbc_1.resize(supercell.get_number_of_bonds_pbc(),4);
  //
  m_bond_angles.resize(supercell.get_number_of_bonds(),2);
  m_bond_angles_pbc.resize(supercell.get_number_of_bonds_pbc(),2);
  m_bond_position.resize(supercell.get_number_of_bonds(),3);
  m_bond_position_pbc.resize(supercell.get_number_of_bonds_pbc(),3);
#ifdef _ATOM_DEBUG_BONDS_
  std::cout<<" ATOM: i_number_of_bonds="<<supercell.get_number_of_bonds()<<std::endl;
  std::cout<<" ATOM: i_number_of_bonds_pbc="<<supercell.get_number_of_bonds_pbc()<<std::endl;
  std::cout<<" ATOM: m_bond_position_pbc="<<m_bond_position_pbc;
#endif
  //u_bond_types=supercell.get_bond_types();
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" ATOM: max number of bonds "<<max_bonds<<std::endl;
  std::cout<<" ATOM: number of bonds "<<supercell.get_number_of_bonds()<<std::endl;
  std::cout<<" ATOM: number of PBC bonds "<<supercell.get_number_of_bonds_pbc()<<std::endl;
  std::cout<<" ATOM: eval_atomic_bonds "<<std::endl;
#endif
  update_atomic_bonds();
#ifdef _SHOW_TIME_
  gl_atom_clock.stop();
  gl_atom_clock.show();
#endif
  // send update color flag
  //update_bonds_color=true;
}

void Fl_Gl_Atom::set_palette(uint u){
  palette.set(u);
  palette.set_color(4);
  palette.initialize(0,u,u);
  palette.update_palette_real();
}

void Fl_Gl_Atom::update_fragments(uint _u, bool _sw){
#ifdef _ATOM_FRAGMENT_MESSAGES_
  std::cout<<" ATOM: begin update_fragments ="<<supercell.get_number_of_fragments()<<std::endl;
#endif
  set_palette(supercell.get_number_of_fragments());
  if(_sw)
    fragment.set_active(supercell.get_fragment_table(_u));
  else
    fragment.set_active(_u);
#ifdef _ATOM_FRAGMENT_MESSAGES_
  std::cout<<" total fragments: "<<supercell.get_number_of_fragments()<<std::endl;
#endif
  set_update_coordinates(true);
  // fragments are counted from 0
#ifdef _ATOM_FRAGMENT_MESSAGES_
  std::cout<<" set active fragment: "<<fragment.get_active()-1<<std::endl;
#endif
  set_active_fragment(0);
  is_eval_bonds=true;
  is_update_bonds=true;
  //update_bonds_color=true;
  update_data();
#ifdef _ATOM_FRAGMENT_MESSAGES_
  std::cout<<" ATOM: end update_fragments ="<<supercell.get_number_of_fragments()<<std::endl;
#endif
}

void Fl_Gl_Atom::compute_atom_fragment(const uint _u){
  // Use fragment table
  update_fragments(_u,true);
}

void Fl_Gl_Atom::compute_radial_fragment(const uint _u, const real _r){
  // Use fragment table
  update_fragments(_u,true);
}

void Fl_Gl_Atom::compute_atom_fragments(void){
  // Use fragment number
  update_fragments(1,false);
}

void Fl_Gl_Atom::compute_vdw_fragment(const uint u){
  supercell.eval_vdw_fragment(u);
  update_fragments(1,false);
}

void Fl_Gl_Atom::compute_vdw_fragments(void){
  // Use fragment number
  supercell.eval_vdw_fragments();
  update_fragments(1,false);
}

void Fl_Gl_Atom::compute_merge_fragments(const uint _u){
  // Use fragment table
  update_fragments(_u,true);
}

void Fl_Gl_Atom::set_active_fragment(const uint u){
  // fragments are counted from 0
  fragment.set_active(u);
  supercell.set_active_fragment(u);
  set_update_coordinates(true);
  update_data(); //<-------------------------
}

void Fl_Gl_Atom::set_atom_active_fragment(const uint u){
  uint _u;
  // fragments are counted from 0
  _u= supercell.get_fragment_table(u);
  fragment.set_active(_u);
  supercell.set_active_fragment(_u);
  set_update_coordinates(true);
  update_data(); //<-------------------------
}

void Fl_Gl_Atom::eval_sphere(uint maxlevel){
    param.u_sphere_rows = 1 << maxlevel;
    int s, cont = 0;
    param.u_sphere_strip_length=20*(param.u_sphere_rows*(param.u_sphere_rows-1)+(param.u_sphere_rows*3));
    m_sphere.resize(param.u_sphere_strip_length,3);
    TVector<real> vt(3);
    /* iterate over the 20 sides of the icosahedron */
    for(s = 0; s < 20; s++) {
        int i;
        triangle *t = (triangle *)&icosahedron[s];
        for(i = 0; i < param.u_sphere_rows; i++) {
            /* create a tstrip for each row */
            /* number of triangles in this row is number in previous +2 */
            /* strip the ith trapezoid block */
            point v0, v1, v2, v3, va, vb;
            int j;
            linearly_interpolate(&t->pt[1], &t->pt[0], (float)(i+1)/param.u_sphere_rows, &v0);
            linearly_interpolate(&t->pt[1], &t->pt[0], (float)i/param.u_sphere_rows, &v1);
            linearly_interpolate(&t->pt[1], &t->pt[2], (float)(i+1)/param.u_sphere_rows, &v2);
            linearly_interpolate(&t->pt[1], &t->pt[2], (float)i/param.u_sphere_rows, &v3);
            NORMV(v0,cont);
            cont++;
            NORMV(v1,cont);
            cont++;
            for(j = 0; j < i; j++) {
                /* calculate 2 more vertices at a time */
                linearly_interpolate(&v0, &v2, (float)(j+1)/(i+1), &va);
                linearly_interpolate(&v1, &v3, (float)(j+1)/i, &vb);
                NORMV(va,cont);
                cont++;
                NORMV(vb,cont);
                cont++;
            }
            NORMV(v2,cont);
            cont++;
        }
    }
#ifdef _ATOM_DEBUG_MESSAGES_
  std::cout<<" FL_GL_ATOM: sphere rows = "<<param.u_sphere_rows<<std::endl;
  std::cout<<" FL_GL_ATOM: sphere_strip_length = "<<param.u_sphere_strip_length<<std::endl;
#endif
}

// cylinder with unitary length
void Fl_Gl_Atom::eval_cylinder(uint n){
  //real theta1,theta2;
  real theta3;
  TVector<real> e(3),t1(2);
  m_cylinder_e1.resize(0,3);
  m_cylinder_texture1.resize(0,2);
  m_cylinder.resize(0,3);
  param.u_cylinder_strip_length=0;
  // check resolution
  if(n < 4){
      n = 4;
  }
  param.u_cylinder_strip_length=n;
  for (uint i=0;i<=n;i++){
    theta3 = i * C_2PI / n;
    e[0] = cos(theta3);
    e[1] = sin(theta3);
    e[2] = 0;
    m_cylinder_e1.add_row(e); // still used
    m_cylinder.add_row(e);
    t1[0] = i/(real)n;
    t1[1] = 1;
    m_cylinder_texture1.add_row(t1);
  }
  m_cylinder.add_row(m_cylinder[n-1]);
  param.u_cylinder_strip_length++;
}

void Fl_Gl_Atom::initialize_transform_matrix(void){
  tb.tb_angle = 0.0;
  glMatrixMode(GL_PROJECTION);
  // put the identity in the trackball transform
  glPushMatrix();
  glLoadIdentity();
  glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat*)tb.tb_transform);
  glPopMatrix();
}

void Fl_Gl_Atom::initialize_rotation_matrix(void){
  tb.tb_angle = 0.0;
  glMatrixMode(GL_MODELVIEW);
  /* put the identity in the trackball transform */
  glPushMatrix();
  glLoadIdentity();
  glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat*)tb.rot_matrix);
  glPopMatrix();
  is_initialize_rot=false;
}

// linearly interpolate between a & b, by fraction f
void Fl_Gl_Atom::linearly_interpolate(point *a, point *b, float f, point *r) {
    r->x=a->x+f*(b->x-a->x);
    r->y=a->y+f*(b->y-a->y);
    r->z=a->z+f*(b->z-a->z);
}

// normalize_point point r
void Fl_Gl_Atom::normalize_point(point *r) {
    float mag;
    mag = r->x * r->x + r->y * r->y + r->z * r->z;
    if (mag != 0.0f) {
        mag = 1.0f / sqrt(mag);
        r->x *= mag;
        r->y *= mag;
        r->z *= mag;
    }
}

// END
