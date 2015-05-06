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

#include <config.h>
#include <iostream>
#include <string>
#include <string.h>
#include <jellyfish/misc.hpp>

typedef int (main_func_t)(int argc, char *argv[]);

main_func_t count_main;
main_func_t bc_main;
main_func_t info_main;
main_func_t stats_main;
main_func_t merge_main;
main_func_t histo_main;
main_func_t query_main;
main_func_t dump_main;
main_func_t cite_main;
// main_func_t dump_fastq_main;
// main_func_t histo_fastq_main;
// main_func_t hash_fastq_merge_main;
main_func_t sos;
main_func_t version;
main_func_t jf_main;

struct cmd_func {
  const char  *cmd;
  //  std::string  cmd;
  main_func_t *func;
};
cmd_func cmd_list[] = {
  {"count",             &count_main},
  {"bc",                &bc_main},
  {"info",              &info_main},
  {"stats",             &stats_main},
  {"histo",             &histo_main},
  {"dump",              &dump_main},
  {"merge",             &merge_main},
  {"query",             &query_main},
  {"cite",              &cite_main},
  // {"qhisto",            &histo_fastq_main},
  // {"qdump",             &dump_fastq_main},
  // {"qmerge",            &hash_fastq_merge_main},
  {"jf",                &jf_main},

  /* help in all its form. Must be first non-command */
  {"help",              &sos},
  {"-h",                &sos},
  {"-help",             &sos},
  {"--help",            &sos},
  {"-?",                &sos},
  {"--version",         &version},
  {"-V",                &version},
  {"",                  0}
};



void __sos(std::ostream *os)
{
  *os << "Usage: jellyfish <cmd> [options] arg..."  << std::endl <<
    "Where <cmd> is one of: ";
  bool comma = false;
  for(cmd_func *ccmd = cmd_list; ccmd->func != sos; ccmd++) {
    *os << (comma ? ", " : "") << ccmd->cmd;
    comma = true;
  }
  *os << "." << std::endl;
  *os << "Options:" << std::endl <<
    "  --version        Display version" << std::endl <<
    "  --help           Display this message" << std::endl;
}

int jf_main(int argc, char* argv[]) {
  const char* aa =
    "                   .......\n"
    "          ..........      .....\n"
    "       ....                   ....\n"
    "      ..     /-+       +---\\     ...\n"
    "      .     /--|       +----\\      ...\n"
    "     ..                              ...\n"
    "     .                                 .\n"
    "     ..      +----------------+         .\n"
    "      .      |. AAGATGGAGCGC .|         ..\n"
    "      .      |---.        .--/           .\n"
    "     ..          \\--------/     .        .\n"
    "     .     .            ..     ..        .\n"
    "     .    ... .....   .....    ..        ..\n"
    "     .   .. . .   .  ..   .   ....        .\n"
    "     .  ..  . ..   . .    ..  .  .         .\n"
    "     . ..   .  .   ...     . ..  ..        .\n"
    "    ....    . ..   ..      ...    ..       .\n"
    "   .. .     ...     .      ..      ..      .\n"
    "   . ..      .      .       .       ...    ..\n"
    "   ...       .      .      ..         ...   .\n"
    "   .         ..     .      ..           .....\n"
    "  ____  ____  ._    __   _  _  ____  ____  ___  _   _\n"
    " (_  _)( ___)(  )  (  ) ( \\/ )( ___)(_  _)/ __)( )_( )\n"
    ".-_)(   )__)  )(__  )(__ \\  /  )__)  _)(_ \\__ \\ ) _ ( \n"
    "\\____) (____)(____)(____)(__) (__)  (____)(___/(_) (_)\n";
  std::cout << aa;
  return 0;
}

int sos(int argc, char *argv[])
{
  __sos(&std::cout);
  return 0;
}

int version(int argc, char *argv[])
{
  std::cout << PACKAGE_STRING << std::endl;
  return 0;
}

int main(int argc, char *argv[])
{
  std::string error;

  if(argc < 2) {
    error = "Too few arguments";
  } else {
    for(cmd_func *ccmd = cmd_list; ccmd->func != 0; ccmd++) {
      if(!strcmp(ccmd->cmd, argv[1]))
      //      if(!ccmd->cmd.compare(argv[1]))
        return ccmd->func(argc - 1, argv + 1);
    }
    error = "Unknown command '";
    error += argv[1];
    error += "'\n";
  }

  std::cerr << error << std::endl;
  __sos(&std::cerr);
  return 1;
}
