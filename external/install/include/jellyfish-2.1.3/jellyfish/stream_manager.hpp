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

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <vector>
#include <list>
#include <set>

#include <jellyfish/locks_pthread.hpp>
#include <jellyfish/err.hpp>

namespace jellyfish {
template<typename PathIterator>
class stream_manager {
  /// A wrapper around an ifstream for a standard file. Standard in
  /// opposition to a pipe_stream below, but the file may be a regular
  /// file or a pipe. The file is opened once and notifies the manager
  /// that it is closed upon destruction.
  class file_stream : public std::ifstream {
    stream_manager& manager_;
  public:
    file_stream(const char* path, stream_manager& manager) :
      std::ifstream(path),
      manager_(manager)
    {
      manager_.take_file();
    }
    virtual ~file_stream() { manager_.release_file(); }
  };
  friend class file_stream;

  /// A wrapper around an ifstream for a "multi pipe". The multi pipe
  /// are connected to generators (external commands generating
  /// sequence). They are opened repeatedly, until they are unlinked
  /// from the file system.
  class pipe_stream : public std::ifstream {
    const char*     path_;
    stream_manager& manager_;
  public:
    pipe_stream(const char* path, stream_manager& manager) :
      std::ifstream(path),
      path_(path),
      manager_(manager)
    { }
    virtual ~pipe_stream() { manager_.release_pipe(path_); }
  };
  friend class pipe_stream;

  typedef std::unique_ptr<std::istream> stream_type;

  PathIterator           paths_cur_, paths_end_;
  int                    files_open_;
  const int              concurrent_files_;
  std::list<const char*> free_pipes_;
  std::set<const char*>  busy_pipes_;
  locks::pthread::mutex_recursive  mutex_;

public:
  define_error_class(Error);

  stream_manager(PathIterator paths_begin, PathIterator paths_end, int concurrent_files = 1) :
    paths_cur_(paths_begin), paths_end_(paths_end),
    files_open_(0),
    concurrent_files_(concurrent_files)
  { }

  stream_manager(PathIterator paths_begin, PathIterator paths_end,
                 PathIterator pipe_begin, PathIterator pipe_end,
                 int concurrent_files = 1) :
    paths_cur_(paths_begin), paths_end_(paths_end),
    files_open_(0),
    concurrent_files_(concurrent_files),
    free_pipes_(pipe_begin, pipe_end)
  { }

  stream_type next() {
    locks::pthread::mutex_lock lock(mutex_);
    stream_type res;
    open_next_file(res);
    if(!res)
      open_next_pipe(res);
    return res;
  }

  int concurrent_files() const { return concurrent_files_; }
  // Number of pipes available. Not thread safe
  int concurrent_pipes() const { return free_pipes_.size() + busy_pipes_.size(); }
  // Number of streams available. Not thread safe
  int nb_streams() const { return concurrent_files() + concurrent_pipes(); }

protected:
  void open_next_file(stream_type& res) {
    if(files_open_ >= concurrent_files_)
      return;
    while(paths_cur_ != paths_end_) {
      std::string path = *paths_cur_;
      ++paths_cur_;
      res.reset(new file_stream(path.c_str(), *this));
      if(res->good())
        return;
      res.reset();
      eraise(std::runtime_error) << "Can't open file '" << path << "'";
    }
  }

  void open_next_pipe(stream_type& res) {
    while(!free_pipes_.empty()) {
      const char* path = free_pipes_.front();
      free_pipes_.pop_front();
      res.reset(new pipe_stream(path, *this));
      if(res->good()) {
        busy_pipes_.insert(path);
        return;
      }
      // The pipe failed to open, so it is not marked as busy. This
      // reset will make us forget about this path.
      res.reset();
    }
  }

  void take_file() {
    locks::pthread::mutex_lock lock(mutex_);
    ++files_open_;
  }

  void release_file() {
    locks::pthread::mutex_lock lock(mutex_);
    --files_open_;
  }

  // void take_pipe(const char* path) {
  //   locks::pthread::mutex_lock lock(mutex_);
  // }
  void release_pipe(const char* path) {
    locks::pthread::mutex_lock lock(mutex_);
    if(busy_pipes_.erase(path) == 0)
      return; // Nothing erased. We forget about that path
    free_pipes_.push_back(path);
  }
};
} // namespace jellyfish
