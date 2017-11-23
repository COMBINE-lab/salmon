/*
Copyright (C) 2011,2012 Remik Ziemlinski. See MIT-LICENSE.

CHANGELOG

v0.0.0 20110502 rsz Created.
V1.0.0 20110522 rsz Extended to show eta with growing bar.
v2.0.0 20110525 rsz Added time elapsed.
v2.0.1 20111006 rsz Added default constructor value.
v2.0.2 20130123 rob Switched over to C++11 timer facilities
*/

#ifndef EZ_ETAPROGRESSBAR_H
#define EZ_ETAPROGRESSBAR_H

#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/ioctl.h>
#include <sys/time.h>
#endif

namespace ez {
// One-line refreshing progress bar inspired by wget that shows ETA (time
// remaining). 90% [========================================>     ] ETA 12d 23h
// 56s
class ezETAProgressBar {
private:
  using system_clock = std::chrono::system_clock;
  using duration = std::chrono::duration<size_t, system_clock::period>;
  using partialDuration = std::chrono::duration<double, system_clock::period>;

public:
  ezETAProgressBar(unsigned int _n = 0) : n(_n), cur(0), pct(0), width(80) {}
  void reset(uint64_t _n) {
    n = _n;
    pct = 0;
    cur = 0;
  }
  void start() {
    startTime = system_clock::now();
    lastCheck = startTime;
    setPct(0);
  }

  void operator++() {
    if (cur >= n)
      return;
    ++cur;
    endTime = system_clock::now();
    if ((endTime - lastCheck) >= std::chrono::seconds(1) or (cur == n)) {
      setPct(static_cast<double>(cur) / n);
      lastCheck = endTime;
    }
  };

  void operator+=(const unsigned int d) {
    if (cur >= n)
      return;
    cur += d;
    endTime = system_clock::now();
    if ((endTime - lastCheck) >= std::chrono::seconds(1) or (cur == n)) {
      setPct(static_cast<double>(cur) / n);
      lastCheck = endTime;
    }
  };

  bool isDone() { return cur == n; }
  void done() {
    cur = n;
    setPct(1.0);
  }

  std::string durationString(duration t) {
    using std::chrono::duration_cast;
    using days = std::chrono::duration<size_t, std::ratio<86400>>;
    using std::chrono::hours;
    using std::chrono::minutes;
    using std::chrono::seconds;

    std::stringstream out(std::stringstream::out);
    // std::string out;

    if (t >= days(1)) {
      auto numDays = duration_cast<days>(t);
      out << numDays.count() << "d ";
      t -= numDays;
    }

    if (t >= hours(1)) {
      auto numHours = duration_cast<hours>(t);
      out << numHours.count() << "h ";
      t -= numHours;
    }

    if (t >= minutes(1)) {
      auto numMins = duration_cast<minutes>(t);
      out << numMins.count() << "m ";
      t -= numMins;
    }

    if (t >= seconds(1)) {
      auto numSecs = duration_cast<seconds>(t);
      out << numSecs.count() << "s";
    }

    std::string tstring = out.str();
    if (tstring.empty()) {
      tstring = "0s";
    }
    return tstring;
  }

  // Set 0.0-1.0, where 1.0 equals 100%.
  void setPct(double Pct) {
    using std::chrono::duration_cast;
    using std::chrono::seconds;
    using weeks = std::chrono::duration<size_t, std::ratio<604800>>;

    endTime = system_clock::now();
    char pctstr[5];
    sprintf(pctstr, "%3d%%", (int)(100 * Pct));

    // Compute how many tics we can display.
    int nticsMax = (width - 27);
    int ntics = std::max(1, static_cast<int>(nticsMax * Pct));
    std::string out(pctstr);
    out.append(" [");

#ifdef HAVE_ANSI_TERM
    // Green!
    out.append("\e[0;32m");
#endif // HAVE_ANSI_TERM

    out.append(std::max(0, ntics - 1), '=');
    out.append(Pct == 1.0 ? "=" : ">");
    out.append(nticsMax - ntics, ' ');

#ifdef HAVE_ANSI_TERM
    // Not-Green :(
    out.append("\e[0m");
#endif // HAVE_ANSI_TERM

    out.append("] ");
    out.append((Pct < 1.0) ? "ETA " : "in ");

    // Time since we started the progress bar (or reset)
    auto dt = endTime - startTime;

    std::string tstr;
    if (Pct >= 1.0) {
      // Print overall time and newline.
      tstr = durationString(dt);
      out.append(tstr);
      size_t effLen = out.length();
#ifdef HAVE_ANSI_TERM
      effLen -= 11;
#endif // HAVE_ANSI_TERM
      if (effLen < width) {
        out.append(width - effLen, ' ');
      }

      out.append("\r\n");
      std::cerr << out;
      return;
    } else {
      duration eta = std::chrono::seconds::max();

      if (Pct > 0.0) {
        duration esecs = duration_cast<seconds>(dt);
        eta = duration_cast<duration>(((esecs * n) / cur) - esecs);
      }

      if (eta > weeks(1)) {
        out.append("> 1 week");
      } else {
        tstr = durationString(eta);
        out.append(tstr);
      }
    }

    size_t effLen = out.length();
#ifdef HAVE_ANSI_TERM
    effLen -= 11;
#endif // HAVE_ANSI_TERM

    // Pad end with spaces to overwrite previous string that may have been
    // longer.
    if (effLen < width) {
      out.append(width - effLen, ' ');
    }

    out.append("\r");
    std::cerr << out;
    std::cerr.flush();
  }

private:
  uint64_t n;
  uint64_t cur;
  unsigned short pct;  // Stored as 0-1000, so 2.5% is encoded as 25.
  unsigned char width; // How many chars the entire line can be.
  std::chrono::system_clock::time_point startTime, endTime, lastCheck;
};
} // namespace ez
#endif // EZ_ETAPROGRESSBAR_H
