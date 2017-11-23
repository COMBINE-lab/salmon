// based off of
// async_client.cpp
// ~~~~~~~~~~~~~~~~
//
// Copyright (c) 2003-2012 Christopher M. Kohlhoff (chris at kohlhoff dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef VERSION_CHECKER_HPP
#define VERSION_CHECKER_HPP

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>

using boost::asio::ip::tcp;

class VersionChecker {
public:
  VersionChecker(boost::asio::io_service& io_service, const std::string& server,
                 const std::string& path);
  std::string message();

private:
  void cancel_upgrade_check(const boost::system::error_code& err);
  void handle_resolve(const boost::system::error_code& err,
                      tcp::resolver::iterator endpoint_iterator);
  void handle_connect(const boost::system::error_code& err);
  void handle_write_request(const boost::system::error_code& err);
  void handle_read_status_line(const boost::system::error_code& err);
  void handle_read_headers(const boost::system::error_code& err);
  void handle_read_content(const boost::system::error_code& err);

  tcp::resolver resolver_;
  tcp::socket socket_;
  boost::asio::streambuf request_;
  boost::asio::streambuf response_;
  boost::asio::deadline_timer deadline_;
  std::stringstream messageStream_;
};

std::string getVersionMessage();

#endif // VERSION_CHECKER_HPP