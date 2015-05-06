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

#include <signal.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

#include <jellyfish/err.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/generator_manager.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/file_header.hpp>
#include <sub_commands/bc_main_cmdline.hpp>

using std::chrono::system_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

template<typename DtnType>
inline double as_seconds(DtnType dtn) { return duration_cast<duration<double>>(dtn).count(); }

static bc_main_cmdline args; // Command line switches and arguments
typedef std::vector<const char*> file_vector;
using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<file_vector::const_iterator> > sequence_parser;
typedef jellyfish::mer_iterator<sequence_parser, jellyfish::mer_dna> mer_iterator;

template<typename PathIterator>
class mer_bloom_counter : public jellyfish::thread_exec {
  int                                     nb_threads_;
  mer_dna_bloom_counter&                  filter_;
  jellyfish::stream_manager<PathIterator> streams_;
  sequence_parser                         parser_;

public:
  mer_bloom_counter(int nb_threads, mer_dna_bloom_counter& filter,
                    PathIterator file_begin, PathIterator file_end,
                    PathIterator pipe_begin, PathIterator pipe_end,
                    uint32_t concurent_files) :
    filter_(filter),
    streams_(file_begin, file_end, pipe_begin, pipe_end, concurent_files),
    parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
  { }

  virtual void start(int thid) {
    for(mer_iterator mers(parser_, args.canonical_flag) ; mers; ++mers) {
      filter_.insert(*mers);
    }
  }
};

// If get a termination signal, kill the manager and then kill myself.
static pid_t manager_pid = 0;
static void signal_handler(int sig) {
  if(manager_pid)
    kill(manager_pid, SIGTERM);
  signal(sig, SIG_DFL);
  kill(getpid(), sig);
  _exit(EXIT_FAILURE); // Should not be reached
}

int bc_main(int argc, char *argv[])
{
  auto start_time = system_clock::now();

  jellyfish::file_header header;
  header.fill_standard();
  header.set_cmdline(argc, argv);

  args.parse(argc, argv);
  mer_dna::k(args.mer_len_arg);

  std::unique_ptr<jellyfish::generator_manager> generator_manager;
  if(args.generator_given) {
    auto gm =
      new jellyfish::generator_manager(args.generator_arg, args.Generators_arg,
                                       args.shell_given ? args.shell_arg : (const char*)0);
    generator_manager.reset(gm);
    generator_manager->start();
    manager_pid = generator_manager->pid();
    struct sigaction act;
    memset(&act, '\0', sizeof(act));
    act.sa_handler = signal_handler;
    assert(sigaction(SIGTERM, &act, 0) == 0);
  }

  header.canonical(args.canonical_flag);
  std::ofstream output(args.output_arg);
  if(!output.good())
    die << "Can't open output file '" << args.output_arg << "'";

  header.format("bloomcounter");
  header.key_len(args.mer_len_arg * 2);
  jellyfish::hash_pair<mer_dna> hash_fns;
  header.matrix(hash_fns.m1, 1);
  header.matrix(hash_fns.m2, 2);

  mer_dna_bloom_counter filter(args.fpr_arg, args.size_arg, hash_fns);
  header.size(filter.m());
  header.nb_hashes(filter.k());
  header.write(output);

  auto after_init_time = system_clock::now();

  // Iterators to the multi pipe paths. If no generator manager,
  // generate an empty range.
  auto pipes_begin = generator_manager.get() ? generator_manager->pipes().begin() : args.file_arg.end();
  auto pipes_end = (bool)generator_manager ? generator_manager->pipes().end() : args.file_arg.end();

  mer_bloom_counter<file_vector::const_iterator> counter(args.threads_arg, filter,
                                                         args.file_arg.begin(), args.file_arg.end(),
                                                         pipes_begin, pipes_end, args.Files_arg);
  counter.exec_join(args.threads_arg);

  // If we have a manager, wait for it
  if(generator_manager) {
    signal(SIGTERM, SIG_DFL);
    manager_pid = 0;
    if(!generator_manager->wait())
      die << "Some generator commands failed";
    generator_manager.reset();
  }

  auto after_count_time = system_clock::now();

  filter.write_bits(output);
  output.close();

  auto after_dump_time = system_clock::now();

  if(args.timing_given) {
    std::ofstream timing_file(args.timing_arg);
    timing_file << "Init     " << as_seconds(after_init_time - start_time) << "\n"
                << "Counting " << as_seconds(after_count_time - after_init_time) << "\n"
                << "Writing  " << as_seconds(after_dump_time - after_count_time) << "\n";
  }

  return 0;
}
