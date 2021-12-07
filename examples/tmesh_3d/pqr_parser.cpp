#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "nanoshaper.h"
using namespace NS;

std::basic_istream<char>& operator>>(std::basic_istream<char>& inputfile, Atom &a) {

  int Atom_number;
  std::string Field_name;

  inputfile >> Field_name
            >> Atom_number
            >> a.ai.name
            >> a.ai.resName
            >> a.ai.resNum
            >> a.pos[0]
            >> a.pos[1]
            >> a.pos[2]
            >> a.charge
            >> a.radius;

  a.radius2 = a.radius*a.radius;
  return inputfile;
}

void
read_atoms_from_pqr (std::basic_istream<char> &inputfile,
                     std::vector<Atom> &atoms) {

  static Atom a;
  atoms.clear ();
  while (inputfile >> a)          
    atoms.push_back (a);
  
}

void
write_atoms_to_pqr (std::basic_ostream<char> &outputfile,
                    const std::vector<Atom> &atoms) {

  int Atom_number = 1;

  outputfile << std::setw(10) << std::left << "fieldname" << std::setw(12) 
             << std::left <<"Atom_number" << std::setw(12) << std::left << "Atom_name" << std::setw(16) << std::left
             << "Residue_name" << std::setw(16) << std::left << "Residue_number" << std::setw(10) << std::left << "X" 
             << std::setw(10) << std::left << "Y" << std::setw(10) << std::left
             << "Z" << std::setw(10) << std::left << "Charge" << std::setw(10) << std::left << "Radius" << std::endl;

  for (auto & ii : atoms) 
    outputfile << std::setw(10) << std::left << "ATOM" << std::setw(12) 
               << std::left << Atom_number++ << std::setw(12) << std::left << ii.ai.name << std::setw(16) << std::left
               << ii.ai.resName << std::setw(16) << std::left << ii.ai.resNum << std::setw(10) << std::left << ii.pos[0] 
               << std::setw(10) << std::left << ii.pos[1] << std::setw(10) << std::left
               << ii.pos[2] << std::setw(10) << std::left << ii.charge << std::setw(10) << std::left << ii.radius << std::endl;
}
