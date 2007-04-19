/*
 *  This file is part of BGX, the Bayesian Gene eXpression program.
 *  Copyright (c) 2003-2004  Graeme Ambler <graeme@ambler.me.uk>
 *
 *  BGX is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 2, as
 *  published by the Free Software Foundation.
 *
 *  BGX is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "rundir.hh"

// ASCII itoa.
void int_to_string(int i,std::string& str)
{
  int digit=i%10;
  i/=10;
  std::string reverse;
  reverse=digit+48;
  while(i!=0)
    {
      digit=i%10;
      i/=10;
      reverse+=digit+48;
    }
  std::string::iterator r_begin=reverse.begin();
  std::string::iterator r_end=reverse.end();
  str=*--r_end;
  while(r_begin!=r_end) str+=*--r_end;
}

std::string rundir(const char* str)
{
  std::string basename;
  basename=str;
  basename+='.';
  std::string dirname;
  std::string num;
  for(int i=1; true; ++i){
    int_to_string(i,num);
    dirname=basename+num;
    DIR* out_dir=opendir(dirname.c_str());
    if(out_dir == NULL){
#ifdef WIN32
      mkdir(dirname.c_str());
#else
      mkdir(dirname.c_str(),S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif
      break;
    }
    else
      closedir(out_dir);
  }
  return dirname;
}
