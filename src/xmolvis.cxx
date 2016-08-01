//========================================================================
// XMOLview - www.molecular-explorer.com
//
// This program can manipulate and genarate structure files
//
//
// Features:
//   + visualization and edition
//   + Topology definition of fragments
//   + Configuration structure scanning
//
//
// Copyrigth 2002-2015 by Edmanuel Torres
//
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

#include<fl_xmol_view.h>
#include<limits.h>
#include<unistd.h>

std::string getexepath()
{
  char result[ PATH_MAX ];
  char* bf = getcwd(result, PATH_MAX);
  //ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  //return std::string(result, (count > 0) ? count : 0 );
  return result;
}

int main(int argc, char* argv[]){
#if defined (BUILD_FOR_LINUX) || defined (BUILD_FOR_MACOS)
  std::string dir_file = getexepath();
#endif
  fl_xmol_view * mv = new fl_xmol_view();
#if defined (BUILD_FOR_LINUX) || defined (BUILD_FOR_MACOS)
  //fl_filename_absolute(full_dir_file, sizeof(full_dir_file), dir_file.c_str());
  //fl_filename_relative(full_dir_file, sizeof(full_dir_file), dir_file.c_str());
  std::cout<<" directory :"<<dir_file<<std::endl;
  mv->open_directory(dir_file.c_str());
  //if ( argc > 1){
    //file_name = argv[1];
    //std::cout<<" file :"<<file_name<<std::endl;
    //mv->open_file(file_name.c_str(),full_dir_file);
  //}
#endif
  mv->show();
  //mv->show_view();
  Fl::run();
  return 0;
}

