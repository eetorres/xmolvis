//========================================================================
// FILE - Fl_Atom_view.cxx                                              //
// For the Fast Light Tool Kit (FLTK) - www.fltk.org                    //
//========================================================================
//                                                                      //
// OpenGL 3D visualization widget to display atomic structures.         //
//                                                                      //
// Copyright 2002-2016 by Edmanuel Torres                               //
// email: eetorres@gmail.com                                            //
//                                                                      //
//======================================================================//
//  This file is part of xmolvis                                        //
//                                                                      //
//  xmolvis is free software: you can redistribute it and/or modify     //
//  it under the terms of the GNU General Public License as published by//
//  the Free Software Foundation, either version 3 of the License, or   //
//  (at your option) any later version.                                 //
//                                                                      //
//  xmolvis is distributed in the hope that it will be useful,          //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of      //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       //
//  GNU General Public License for more details.                        //
//                                                                      //
//  You should have received a copy of the GNU General Public License   //
//  along with xmolvis.  If not, see <http://www.gnu.org/licenses/>.    //
//======================================================================//

#ifndef _FL_GLMOL_VIEW_H_
#define _FL_GLMOL_VIEW_H_

#include <stdlib.h>
#include <stdio.h>

#include <atom_color.h>
#include <atom_symbol.h>
#include <atom_name.h>
#include <msmvtl/const.h>
#include <msmvtl/linalg.h>
#include <assert.h>

#include <FL/Fl.H>


#define HAVE_GL 1
#if HAVE_GL
    #include <fl_gl_atom.h>
//    #include <stroke.h>
#else
    #include <FL/Fl_Box.H>
#endif // HAVE_GL

#define MODE_RENDER    1
#define MODE_SELECT    2
#define MODE_MENU      3

#define BUFSIZE     1024

// Main menues
#define MAIN_MENU    0
#define CONTROL_MENU 1
#define ATOM_MENU    2
#define NOT_MENU     3

const uint NUMBER_OF_SLIDERS  = 4;
const uint NUMBER_OF_RADIOS   = 9;

// Default menu
const std::string lblank[6] = { "blank", "blank", "blank", "blank",   "close", "blank" }; // 0
//                                 0        1        2        3          4        5
// Main menu
const std::string l0[6] = { "mode", "view",  "show",  "tools", "close", "colors" }; // 0
//                             0       1        2        3        4        5
// Controls menu
const std::string l1[6] = { "menu", "color", "zoom",  "move", "close", "copy" };   // 1
//                            0        1        2        3      4        5
// Atom menu
const std::string l2[6] = { "edit", "axis", "select", "frag", "close", "fix" };   // 2
//                            0      1        2        3        4       5
//
// Main menu buttons
#define MAIN_MODE_SUBMENU     0
#define MAIN_VIEW_SUBMENU     1
#define MAIN_SHOW_SUBMENU     2
#define MAIN_TOOLS_SUBMENU    3
#define MAIN_COLOR_SUBMENU    5
// Atom menu buttons
#define ATOM_ELEMENT_SUBMENU  0
#define ATOM_AXIS_SUBMENU     1
#define ATOM_ACTIVE_SUBMENU   2
#define ATOM_FRAGMENT_SUBMENU 3
#define ATOM_DELETE_SUBMENU   5
//
// Submenues                      0       1         2         3        4         5
// Main menu
// Mode submenu                     PBC     symetry
const std::string sl_mode[6] = { "--",    "--",    "frags",   "atoms", "cancel", "--" };
// View submenu
const std::string sl_view[6] = { "XY (z)", "YZ (x)", "ZX (y)", "YX (Z)", "ZY (X)", "XZ (Y)" };
// Show submenu
const std::string sl_show[6] = { "labels", "box",  "symbols", "bonds", "cancel", "numbers" };
// Tools submenu
const std::string sl_tools[6] = { "params", "box",  "axes", "split", "cancel", "atoms" };


// Submenues                       0       1         2         3        4         5
// Atom menu
// Fragment submenu
const std::string l_frag[6]  = { "vdW",   "atom", "delete", "merge", "cancel", "radial" };
// Axis submenu
const std::string l_axis[6]  = { "show",   "head", "tail", "plane", "cancel", "zero" };
// Edit submenu
const std::string l_edit[6]  = { "delete",   "copy", "symbol", "color", "cancel", "size" };

#define CLOSE_MENU             4
#define CLOSE_SUBMENU          4

// Mode submenu buttons
#define MODE_PBC_BUTTON        0
#define MODE_SYMMETRY_BUTTON   1
#define MODE_FRAGMENTS_BUTTON  2
#define MODE_ATOMS_BUTTON      3
#define MODE_BLANK_BUTTON      5

// View submenu buttons
#define VIEW_XY_BUTTON         0
#define VIEW_YZ_BUTTON         1
#define VIEW_ZX_BUTTON         2
#define VIEW_YX_BUTTON         3
#define VIEW_ZY_BUTTON         4
#define VIEW_XZ_BUTTON         5

// Show submenu buttons
#define SHOW_LABELS_BUTTON     0
#define SHOW_BOX_BUTTON        1
#define SHOW_SYMBOLS_BUTTON    2
#define SHOW_BONDS_BUTTON      3
#define SHOW_NUMBERS_BUTTON    5

// Tools submenu buttons
#define TOOLS_MESURE_BUTTON    0
#define TOOLS_BOX_BUTTON       1
#define TOOLS_AXES_BUTTON      2
#define TOOLS_FRAGMENT_BUTTON  3
#define TOOLS_ATOMS_BUTTON     5

// Fragment submenu buttons
#define ATOM_FRAGMENT_VDW_BUTTON    0
#define ATOM_FRAGMENT_ATOM_BUTTON   1
#define ATOM_FRAGMENT_DELETE_BUTTON 2
#define ATOM_FRAGMENT_MERGE_BUTTON  3
#define ATOM_FRAGMENT_RADIAL_BUTTON 5

// Axis submenu buttons
#define AXIS_SHOW_BUTTON       0
#define AXIS_HEAD_BUTTON       1
#define AXIS_TAIL_BUTTON       2
#define AXIS_PLANE_BUTTON      3
#define AXIS_ZERO_BUTTON       5

//#define FRAGMENT_MENU 3
#define AUTO_SUBMENU           11

class CSelection {

public:

  void clear(void){
    __highlight_atom        = 0;
    __last_highlight_atom   = 0;
    __highlight_atom_a      = 0;
    __highlight_atom_b      = 0;
    __select_begin          = -1;
    __select_end            = -1;
  }

  int  __highlight_atom_a;
  int  __highlight_atom_b;
  int  __select_begin;
  int  __select_end;
  uint __highlight_atom;
  uint __last_highlight_atom;
};

class CTools {

public:

  CTools(void){
    v_distance1.resize(3);
    v_distance2.resize(3);
    v_distance3.resize(3);
  }

  void clear(void){
    u_selected_index = 0;
    r_distance1      = 0.0;
    r_distance2      = 0.0;
    r_distance3      = 0.0;
    r_angle1         = 0.0;
    r_angle2         = 0.0;
    r_dihedral       = 0.0;
    v_distance1.zero();
    v_distance2.zero();
    v_distance3.zero();
  }

//private:
  uint u_selected_index;

  real r_distance1;
  real r_distance2;
  real r_distance3;
  real r_angle1;
  real r_angle2;
  real r_dihedral;

  TVector<real> v_distance1;
  TVector<real> v_distance2;
  TVector<real> v_distance3;
};


class CViewsetup {

public:

  CViewsetup(void){
    // Foreground colors
    fgred   = 1.0; //0.5;
    fggreen = 1.0;
    fgblue  = 1.0; //0.5;
    // Background colors
    bgred   = 0.0;
    bggreen = 0.0;
    bgblue  = 0.1;
    // Appareance
    f_bond_brightness          = 0.8;
    f_background_brightness    = 0.0;
    f_highlight_brightness     = 1.0;
    f_select_brightness        = 1.0;
    f_highlight_brightness_max = 1.5;
    f_select_brightness_max    = 1.5;
    f_atom_brightness          = 0.8;
    f_atom_brightness_max      = 1.5;

  }

  void set_background(real r1, real r2, real r3){
    bgred=r1; bggreen=r2; bgblue=r3;
  }

  void set_forekground(real r1, real r2, real r3){
    fgred=r1; fggreen=r2; fgblue=r3;
  }

  real fgred, fggreen, fgblue;
  real bgred, bggreen, bgblue;
  real f_atom_brightness;
  real f_atom_brightness_max;
  real f_bond_brightness;
  real f_select_brightness;
  real f_select_brightness_max;
  real f_highlight_brightness;
  real f_highlight_brightness_max;
  real f_background_brightness;
};


class CGLWindow {

  public:

  CGLWindow(void){
    base_view   =  10.0;
    view_left   =  base_view;
    view_right  =  base_view;
    view_bottom =  base_view;
    view_top    =  base_view;
    view_near   = -1000;
    view_far    =  1000;
    view_axis_x = -base_view+8.0;
    view_axis_y = -base_view+2.0;
    zoom        = 1.0;
    shift_factor = 0.05;
    //
    // View direction
    __x_ang = 0.0;
    __y_ang = 0.0;
    __z_ang = 0.0;
    // View position and scale
    //size    = 10.0;
    x_shift = 0.0;
    y_shift = 0.0;
    z_shift = 0.0;
    y_off   = 1.0;
  }

  // Gl window parameters
  real base_view;
  real view_left;
  real view_right;
  real view_bottom;
  real view_top;
  real view_near;
  real view_far;
  real x_factor, y_factor;
  real view_axis_x, view_axis_y;
  //
  real x_shift, y_shift, z_shift;
  real __x_ang, __y_ang, __z_ang;
  real y_off; //, _scl, 
  real zoom_step;
  real zoom;
  real shift_factor;
};


class Fl_Gl_Mol_View : public Fl_Gl_Atom{

public:

  //real size;
  //gm_rgb oldrgb, rgb;
  Fl_Gl_Mol_View();
  Fl_Gl_Mol_View(int,int,int,int,const char* l=0);
  ~Fl_Gl_Mol_View();
  //void initData(void);
  void clear(void);
  bool initialize(void);
  void graph_cb(void);
  //  void redraw(void){draw();};
  void set_background_color(real,real,real);
  void set_foreground_color(real,real,real);
  // Scale function
  void set_zoom(real f){ if(f>0.0001) glview.zoom = f;};
  void set_zoom_step(real f){ if(f>0.0001) glview.zoom_step = f;};
  // the rotation about the vertical (x) axis
  void x_angle(real f){ glview.__x_ang = f;redraw();};
  real x_angle(){return glview.__x_ang;};
  // the rotation about the horizontal (y) axis
  void y_angle(real f){ glview.__y_ang = f;redraw();};
  real y_angle(){return glview.__y_ang ;};
  // the rotation about the (z) axis
  void z_angle(real f){ glview.__z_ang = f;redraw();};
  real z_angle(){return glview.__z_ang ;};
  // Begin GUI controls
  void set_position_x(real f){glview.x_shift=f;};
  void set_position_y(real f){glview.y_shift=f;};
  void set_position_z(real f){glview.z_shift=f;};
  //
  void set_view_xy_front(void);
  void set_view_yz_front(void);
  void set_view_zx_front(void);
  void set_view_xy_back(void);
  void set_view_yz_back(void);
  void set_view_zx_back(void);
  void set_view(real,real,real);
  void set_highlight_atom(int);
  void set_highlight_atom_a(int);
  void set_highlight_atom_b(int);
  void set_select_begin(int);
  void set_select_end(int);
  //
  void set_selected_atom(uint);
  void set_active_slider(uint);
  void set_active_radio(uint,bool);
  //
  void set_update_active_fragment(void);
  //
  void set_atom_brightness(real);
  void set_bond_brightness(real);
  void set_background_brightness(real);
  void set_highlihght_brightness(real);
  void set_select_brightness(real);

  // Set radius scaling factor
  void set_atom_radius_scale(real);
  void set_bond_radius_scale(real);
  //
  //void set_basis_vectors(const TVector<real>&,const TVector<real>&,const TVector<real>&);
  // End GUI controls
  void set_atomic_cut_radius(void);
  //
  uint get_highlight_atom(void);
  //uint get_action(void){ return last_action;};
  //
  void is_graphics(bool);
  void is_highlight_atom(bool);
  void is_mask_atoms(bool);
  void is_draw_bbox(bool);
  void is_draw_world_axes(bool);
  void is_draw_molecular_axes(bool);
  void is_draw_molecular_axis(bool);
  void is_draw_bonds(bool);
  void is_draw_labels(bool);
  void is_draw_symbols(bool);
  void is_draw_tools(bool);
  void is_draw_numbers(bool);
  void set_lock_controls(bool);
  void is_highlight_fragment(bool);
  void set_pcb(bool);
  //inline int sign_of(int i){ return (i==0)?0:(i<0?-1:1); }
  int WindowDump(void);
  int handle(int);
  void handle_main_menu(void);
  void handle_atom_menu(void);
  void handle_controls_menu(void);
  //
  // void set_update_coordinates(bool b){ update_coordinates=b;};
  virtual void view_redraw(void){ redraw();};
  //

#if HAVE_GL
#else
#endif // HAVE_GL
#if HAVE_GL
    void draw();
#endif // HAVE_GL

private:

  // GL parameters
  CGLWindow glview;
  CViewsetup setup;
  // Tools
  CTools tools;
  // Selections
  CSelection marker;
  //
  int font_size_symbol;
  int font_size_panel_label;
  int font_size_slider_label;
  int font_size_pie_label;
  //
  int u_slider_index;
  int u_radio_index;
  //
  uint u_menu_index;
  uint u_submenu_index;
  uint u_unselected_atom;
  uint u_active_menu;
  //
  real f_atom_bond_delta;
  real f_atom_radius_scale;
  real f_ball_radius1_scale;
  real f_ball_radius2_scale;
  real f_bond_radius_scale;
  //
  real r_cut_radius, r_cut_radius_2;
  //
  TVector<real> v_axes_position;
  //
  GLfloat  scaled_light_position[4];
  GLdouble menu_pos_cx, menu_pos_cy;
  GLdouble menu_pos_x, menu_pos_y, menu_pos_z;
  GLdouble click_pos_x, click_pos_y;
  GLdouble side_pos_x, side_pos_y;
  GLdouble submenu_pos_cx, submenu_pos_cy;
  GLdouble submenu_pos_x, submenu_pos_y, submenu_pos_z;
  GLdouble label_atom_pos_x, label_atom_y, label_atom_z;
  GLdouble label_symbol_pos_x, label_symbol_y, label_symbol_z;
  GLdouble label_menu_pos_x, label_menu_y, label_menu_z;
  //
  uint v_selected_atoms[4];
  unsigned char pixel[3];
  //
  std::string label;
  std::string legends[6];
  std::string sub_label;
  std::string sub_legends[6];
  // atomic spheres
  GLuint    sphere_dl;
  // atomic bonds
  GLuint    cylinder_dl;
  // what is displayed
  bool is_graphics_on, is_draw_world_axes_;
  bool is_draw_molecular_axes_, is_draw_molecular_axis_;
  bool is_normal_color_, is_highlight_atom_;
  bool is_dark_mask_, is_highlight_fragment_;
  bool is_pbc;
  bool is_draw_processing;
  bool is_draw_menu;
  bool is_draw_controls;
  bool is_lock_controls;
  bool is_draw_tools_;
  bool is_draw_pie_menu;
  bool is_draw_pie_submenu;
  bool is_draw_line;
  bool is_draw_point;
  bool is_atom_picked;
  bool is_background_picked;
  bool is_menu_picked;
  bool is_menu_pie_picked;
  bool is_submenu_pie_picked;
  bool is_menu_position;
  bool is_lock_dragging;
  bool is_mode_atom;
  bool is_mode_fragment;
  bool is_left_click;
  bool is_right_click;
  //
  bool is_update_position;
  bool is_update_bonds_color;
  bool is_update_normal_color;
  bool is_update_dark_mask;
  bool is_update_highlight_fragment;
  bool is_update_highlight_atom;
  bool is_update_selected_atoms;
  //
  bool is_unselected_atom;
  bool is_handle_atom_;
  bool is_handle_main_;
  bool is_control_left_on;
  bool is_slider_active[NUMBER_OF_SLIDERS];
  bool is_radio_active[NUMBER_OF_RADIOS];
  //
  void eval_system_properties(void);
  void eval_mask_rcolor(void);
  void eval_tool_parameters(void);
  // drawing functions
  void clear_all(void);
  void clear_scene(void);
  void draw_scene(void);
  void draw_atoms(void);
  void draw_axes(void);
  void draw_bonds(void);
  void draw_symbols(void);
  void draw_selected_numbers(void);
  void draw_box(void);
  // GL-GUI /////////////////////////////////
  void resize(int X,int Y,int W,int H);
  void draw_pie_menu(GLfloat cx, GLfloat cy, GLfloat z, GLfloat r, GLint num_segments);
  void draw_pie_submenu( GLfloat z, GLfloat r, GLint num_segments);
  void set_pie_labels(const std::string l[], std::string);
  void draw_controls(GLfloat z);
  void draw_message(GLfloat z);
  void draw_tools(GLfloat z);
  void draw_settings(GLfloat z);
  void draw_information(GLfloat z);
  void draw_slider(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, GLfloat val, bool active, char* l);
  void draw_radio_button(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, GLfloat val, bool active, char* l);
  void draw_switch_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, bool val, char* l);
  void widget_float_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, float val, char* l);
  void widget_float_output_xy(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, int x, int y, float val, char* l);
  void widget_text_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, char* text, char* l);
  void widget_int_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, int val, char* l);
  void widget_vector_output(GLfloat x1, GLfloat y1, GLfloat y2, GLfloat z, TVector<real> v, char* l);
  //void draw_pie_submenu(GLfloat cx, GLfloat cy, GLfloat z, GLfloat r, GLint num_segments);
  //void draw_pie_sub_menu(GLfloat,GLfloat,GLfloat,GLint);
  //void draw_pie(GLfloat x1, GLfloat y1, GLfloat z, std::string l[], GLint nl, GLfloat r, GLint n);
  //void draw_sub_pie(GLfloat x1, GLfloat y1, GLfloat z, std::string l[], GLint nl, GLfloat r, GLint n);
  //void draw_pie_disk(GLfloat x, GLfloat y, GLfloat z, GLfloat r, GLint n);
  //void draw_pie_labels(GLfloat cx, GLfloat cy, GLfloat z, GLfloat r, std::string l[], std::string m, GLint nl);
  //
  // GUI utils
  inline void view_reshape(int width, int height);
  void set_font_size(void);
  void initialize_opengl(void);
  void process_picking(unsigned char pc[3]);
  //
  void create_sphere_dl(void);
  void delete_sphere_dl(void);
  void create_cylinder_dl(void);
  //
  inline int sign_of(int i){
    return (i<0 ?-1:1);
  };
  //
  inline void set_mouse_motion(int,int);
  inline void point_to_vector(int, int, int, int,GLfloat v[3]);
  // font functions
  //  drawLetter() interprets the instructions from the array
  //  for that letter and renders the letter with line segments.
  //void drawLetter(CPfont *l);
  // Create a display list for each of 6 characters
  //void init_font (void);
  //void printStrokedString(char *s);
};

#endif //

