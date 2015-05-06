#include <iostream>
#include <fstream>
#include <string>

#include <jellyfish/err.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/misc.hpp>
#include <sub_commands/info_main_cmdline.hpp>

static info_main_cmdline args;

std::string get_command(const jellyfish::generic_file_header& h) {
  std::string cmd(h["exe_path"]);
  std::vector<std::string> cmdline = h.cmdline();
  for(auto it = cmdline.cbegin(); it != cmdline.cend(); ++it)
    (cmd += " ") += jellyfish::quote_arg(*it);

  return cmd;
}

std::string get_where(const jellyfish::generic_file_header& h) {
  std::string res(jellyfish::quote_arg(h["hostname"]));
  if(!res.empty())
    res += ":";
  res += jellyfish::quote_arg(h["pwd"]);
  return res;
}

int info_main(int argc, char *argv[]) {
  args.parse(argc, argv);

  std::ifstream file(args.file_arg);
  if(!file.good())
    die << "Can't open '" << args.file_arg << "'" << jellyfish::err::no;

  jellyfish::file_header header;
  header.read(file);

  if(args.skip_flag)
    std::cout << file.rdbuf();
  else if(args.json_flag)
    std::cout << header;
  else if(args.cmd_flag)
    std::cout << get_command(header) << "\n";
  else
    std::cout << "command: " << get_command(header) << "\n"
              << "where: " << get_where(header) << "\n"
              << "when: " << header["time"] << "\n"
              << "canonical: " << (header.canonical() ? "yes" : "no") << "\n";

  return 0;
}
