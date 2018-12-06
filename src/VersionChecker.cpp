// based off of
// async_client.cpp
// ~~~~~~~~~~~~~~~~
//
// Copyright (c) 2003-2012 Christopher M. Kohlhoff (chris at kohlhoff dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include "VersionChecker.hpp"
#include "SalmonConfig.hpp"

using boost::asio::ip::tcp;

VersionChecker::VersionChecker(boost::asio::io_service& io_service,
                               const std::string& server,
                               const std::string& path)
    : resolver_(io_service), socket_(io_service), deadline_(io_service) {
  // Form the request. We specify the "Connection: close" header so that the
  // server will close the socket after transmitting the response. This will
  // allow us to treat all data up until the EOF as the content.
  std::ostream request_stream(&request_);
  request_stream << "GET " << path << " HTTP/1.0\r\n";
  request_stream << "Host: " << server << "\r\n";
  request_stream << "Accept: */*\r\n";
  request_stream << "Connection: close\r\n\r\n";

  deadline_.async_wait(boost::bind(&VersionChecker::cancel_upgrade_check, this,
                                   boost::asio::placeholders::error));
  deadline_.expires_from_now(boost::posix_time::seconds(1));

  // Start an asynchronous resolve to translate the server and service names
  // into a list of endpoints.
  tcp::resolver::query query(server, "http");
  resolver_.async_resolve(query,
                          boost::bind(&VersionChecker::handle_resolve, this,
                                      boost::asio::placeholders::error,
                                      boost::asio::placeholders::iterator));
}

std::string VersionChecker::message() { return messageStream_.str(); }

void VersionChecker::cancel_upgrade_check(
    const boost::system::error_code& err) {
  if (err != boost::asio::error::operation_aborted) {
    deadline_.cancel();
    messageStream_
        << "Could not resolve upgrade information in the alotted time.\n";
    messageStream_ << "Check for upgrades manually at "
                      "https://combine-lab.github.io/salmon\n";
    socket_.close();
  }
}

void VersionChecker::handle_resolve(const boost::system::error_code& err,
                                    tcp::resolver::iterator endpoint_iterator) {
  if (!err) {
    // Attempt a connection to each endpoint in the list until we
    // successfully establish a connection.
    boost::asio::async_connect(socket_, endpoint_iterator,
                               boost::bind(&VersionChecker::handle_connect,
                                           this,
                                           boost::asio::placeholders::error));
  } else {
    deadline_.cancel();
    cancel_upgrade_check(err);
  }
}

void VersionChecker::handle_connect(const boost::system::error_code& err) {
  if (!err) {
    // The connection was successful. Send the request.
    boost::asio::async_write(socket_, request_,
                             boost::bind(&VersionChecker::handle_write_request,
                                         this,
                                         boost::asio::placeholders::error));
  } else {
    deadline_.cancel();
    cancel_upgrade_check(err);
  }
}

void VersionChecker::handle_write_request(
    const boost::system::error_code& err) {
  if (!err) {
    // Read the response status line. The response_ streambuf will
    // automatically grow to accommodate the entire line. The growth may be
    // limited by passing a maximum size to the streambuf constructor.
    boost::asio::async_read_until(
        socket_, response_, "\r\n",
        boost::bind(&VersionChecker::handle_read_status_line, this,
                    boost::asio::placeholders::error));
  } else {
    deadline_.cancel();
    cancel_upgrade_check(err);
  }
}

void VersionChecker::handle_read_status_line(
    const boost::system::error_code& err) {
  if (!err) {
    // Check that response is OK.
    std::istream response_stream(&response_);
    std::string http_version;
    response_stream >> http_version;
    unsigned int status_code;
    response_stream >> status_code;
    std::string status_message;
    std::getline(response_stream, status_message);
    if (!response_stream || http_version.substr(0, 5) != "HTTP/") {
      deadline_.cancel();
      cancel_upgrade_check(err);
      return;
    }
    if (status_code != 200) {
      deadline_.cancel();
      cancel_upgrade_check(err);
      return;
    }

    // Read the response headers, which are terminated by a blank line.
    boost::asio::async_read_until(
        socket_, response_, "\r\n\r\n",
        boost::bind(&VersionChecker::handle_read_headers, this,
                    boost::asio::placeholders::error));
  } else {
    deadline_.cancel();
    cancel_upgrade_check(err);
  }
}

void VersionChecker::handle_read_headers(const boost::system::error_code& err) {
  if (!err) {
    deadline_.cancel();
    // Process the response headers.
    std::istream response_stream(&response_);
    std::string header;
    while (std::getline(response_stream, header) && header != "\r") {
    }

    // Write whatever content we already have to output.
    if (response_.size() > 0)
      messageStream_ << &response_;

    // Start reading remaining data until EOF.
    boost::asio::async_read(
        socket_, response_, boost::asio::transfer_at_least(1),
        boost::bind(&VersionChecker::handle_read_content, this,
                    boost::asio::placeholders::error));
  } else {
    deadline_.cancel();
    cancel_upgrade_check(err);
  }
}

void VersionChecker::handle_read_content(const boost::system::error_code& err) {
  if (!err) {
    // Write all of the data that has been read so far.
    messageStream_ << &response_;

    // Continue reading remaining data until EOF.
    boost::asio::async_read(
        socket_, response_, boost::asio::transfer_at_least(1),
        boost::bind(&VersionChecker::handle_read_content, this,
                    boost::asio::placeholders::error));
  } else if (err != boost::asio::error::eof) {
    deadline_.cancel();
    cancel_upgrade_check(err);
  }
}

std::string getVersionMessage() {
  std::string baseSite{"combine-lab.github.io"};
  std::string path{"/salmon/version_info/"};
  path += salmon::version;

  std::stringstream ss;
  try {
    boost::asio::io_service io_service;
    VersionChecker c(io_service, baseSite, path);
    io_service.run();
    ss << "Version Info: " << c.message();
  } catch (std::exception& e) {
    ss << "Version Info Exception: " << e.what() << "\n";
  }

  return ss.str();
}
