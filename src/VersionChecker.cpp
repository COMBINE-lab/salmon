#include "VersionChecker.hpp"
#include "SalmonConfig.hpp"
#include "httplib.hpp"

std::string getVersionMessage() {
  std::stringstream ss;
  try {
    // NOTE: getaddrinfo / freeaddrinfo will cause a "memory leak"
    // once per address, program pair.  This is a known issue 
    // https://lists.debian.org/debian-glibc/2016/03/msg00243.html.
    // If valgrind leads you here, best not to worry about it.
    httplib::Client cli("combine-lab.github.io");
    std::string path{"/salmon/version_info/"};
    path += salmon::version;
    cli.set_timeout_sec(2); // timeouts in 2 seconds
    auto res = cli.Get(path.c_str());
    if (res) { // non-null response
      if (res->status == 200) { // response OK
        ss << "Version Info: " << res->body;
      } else { // response something else
        ss << "Version Server Response: " << httplib::detail::status_message(res->status) << "\n";
      }
    } else { // null response
      ss << "Version Info Exception: server did not respond before timeout\n";
    }
  } catch (std::exception& e) {
    ss << "Version Info Exception: " << e.what() << "\n";
  }

  return ss.str();
}
