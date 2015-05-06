/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

const char *cite =
  "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers\n"
  "Guillaume Marcais; Carl Kingsford\n"
  "Bioinformatics (2011) 27(6): 764-770 first published online January 7, 2011 doi:10.1093/bioinformatics/btr011\n";

const char *url =
  "http://www.cbcb.umd.edu/software/jellyfish\n"
  "http://bioinformatics.oxfordjournals.org/content/early/2011/01/07/bioinformatics.btr011";

const char *bibtex =
  "@article{Jellyfish2010,\n"
  "         author = {Mar\\c{c}ais, Guillaume and Kingsford, Carl},\n"
  "         title = {A fast, lock-free approach for efficient parallel counting of occurrences of k-mers},\n"
  "         volume = {27},\n"
  "         number = {6},\n"
  "         pages = {764-770},\n"
  "         year = {2011},\n"
  "         doi = {10.1093/bioinformatics/btr011},\n"
  "         URL = {http://bioinformatics.oxfordjournals.org/content/27/6/764.abstract},\n"
  "         eprint = {http://bioinformatics.oxfordjournals.org/content/27/6/764.full.pdf+html},\n"
  "         journal = {Bioinformatics}\n"
  "}";

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/fstream_default.hpp>
#include <sub_commands/cite_main_cmdline.hpp>

int cite_main(int argc, char *argv[])
{
  cite_main_cmdline args(argc, argv);

  ofstream_default out(args.output_given ? args.output_arg : 0, std::cout);
  if(!out.good())
    die << "Can't open output file '" << args.output_arg << "'" << jellyfish::err::no;

  if(args.bibtex_flag) {
    out << bibtex << std::endl;
  } else {
    out << "This software has been published. If you use it for your research, cite:\n\n"
        << cite << "\n\n" << url << std::endl;
  }
  out.close();

  return 0;
}
