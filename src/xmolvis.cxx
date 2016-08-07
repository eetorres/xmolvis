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
//  along with xmolvis.  If not, see <http://www.gnu.org/licenses/>.     //
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
  std::string str_file;
  std::string str_dir;
  char chr_file[PATH_MAX];
  char chr_dir_file[PATH_MAX];
#endif
  fl_xmol_view * mv = new fl_xmol_view();
#if defined (BUILD_FOR_LINUX) || defined (BUILD_FOR_MACOS)
  str_dir = getexepath();
  //mv->open_directory(str_dir.c_str());
  if ( argc > 1){
    str_file = argv[1];
    fl_filename_absolute(chr_dir_file,sizeof(chr_dir_file),str_file.c_str());
    std::cout<<" file :"<<str_file<<std::endl;
    std::cout<<" directory :"<<chr_dir_file<<std::endl;
    mv->open_directory(chr_dir_file);
    //fl_filename_relative(full_dir_file, sizeof(full_dir_file), dir_file.c_str());
    //mv->open_file(file_name.c_str(),full_dir_file);
  }
#endif
  mv->show();
  //mv->show_view();
  Fl::run();
  return 0;
}

