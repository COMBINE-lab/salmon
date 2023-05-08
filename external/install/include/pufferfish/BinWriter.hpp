#ifndef __BIN_WRITER_H__
#define __BIN_WRITER_H__

#include <vector>
#include <ostream>
#include <mutex>
#include <iterator>

#pragma once

#include "spdlog/details/null_mutex.h"
#include "spdlog/sinks/base_sink.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/sinks/ansicolor_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"
#include "string_view.hpp"

class BinWriter
{
    private:
        std::vector<char> _bin_data;
        
    public:
        BinWriter() {};
        BinWriter(size_t reserve_size) {_bin_data.reserve(reserve_size);}
        //copy bin data to this record
        BinWriter(const std::vector<char>& bin_data) :_bin_data(bin_data) {};
        //or just move if possible
        BinWriter(std::vector<char>&& bin_data) :_bin_data(move(bin_data)) {};

        void clear() {
            _bin_data.clear();
        }

        BinWriter& operator<<(const bool &inval) {
            char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
            std::copy(inCharPtr, inCharPtr+sizeof(inval),
                    std::back_inserter(_bin_data));
            return *this;
        }
        BinWriter& operator<<(const uint8_t &inval) {
            char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
            std::copy(inCharPtr, inCharPtr+sizeof(inval),
                    std::back_inserter(_bin_data));
            return *this;
        }
        BinWriter& operator<<(const uint16_t &inval) {
            char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
            std::copy(inCharPtr, inCharPtr+sizeof(inval),
                    std::back_inserter(_bin_data));
            return *this;
        }
        BinWriter& operator<<(const uint32_t &inval) {
            char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
            std::copy(inCharPtr, inCharPtr+sizeof(inval),
                    std::back_inserter(_bin_data));
            return *this;
        }
        BinWriter& operator<<(const uint64_t &inval) {
            char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
            std::copy(inCharPtr, inCharPtr+sizeof(inval),
                    std::back_inserter(_bin_data));
            return *this;
        }
    BinWriter& operator<<(const int32_t &inval) {
        char* inCharPtr = const_cast<char*>(reinterpret_cast<const char*>(&inval));
        std::copy(inCharPtr, inCharPtr+sizeof(inval),
                  std::back_inserter(_bin_data));
        return *this;
    }

    BinWriter &operator<<(const double &inval) {
        char *inCharPtr = const_cast<char *>(reinterpret_cast<const char *>(&inval));
        std::copy(inCharPtr, inCharPtr + sizeof(inval),
                  std::back_inserter(_bin_data));
        return *this;
    }
        BinWriter& operator<<(const std::string &inval) {
            if (inval.size() >= 0x100) {
                std::cerr << "ERROR!! DOESN'T SUPPORT STRING LENGTH LONGER THAN 255. String length: " 
                          << inval.size() << "\n";
                std::exit(1);
            }
            char tmp = static_cast<uint8_t>(inval.size());
            _bin_data.push_back(tmp);
            //(*this) << inval.size();
            //std::cout << inval.size() << " " << inval.c_str()  << " " << inval << "\n";
            char* inCharPtr = const_cast<char*>(inval.c_str());
            std::copy(inCharPtr, inCharPtr+inval.size(),
                    std::back_inserter(_bin_data));
            return *this;
        }
        BinWriter& operator<<(const stx::string_view &inval) {
            if (inval.size() >= 0x100) {
                std::cerr << "ERROR!! DOESN'T SUPPORT STRING LENGTH LONGER THAN 255. String length: " 
                          << inval.size() << "\n";
                std::exit(1);
            }
            char tmp = static_cast<uint8_t>(inval.size());
            _bin_data.push_back(tmp);
            //(*this) << inval.size();
            //std::cout << inval.size() << " " << inval.data()  << " " << inval << "\n";
            char* inCharPtr = const_cast<char*>(inval.data());
            std::copy(inCharPtr, inCharPtr+inval.size(),
                    std::back_inserter(_bin_data));
            return *this;
        }

    uint64_t getBytes() { return _bin_data.size(); }
        //support for logging directly from spdlog
        template<typename OStream>
        friend OStream& operator<<(OStream& os, const BinWriter &bin_record)
        {        
            std::ostream_iterator<char> out_iter(os);		
            std::copy(bin_record._bin_data.begin(), bin_record._bin_data.end(), out_iter);		
            return os;
        }
};



template<class Mutex>
class ostream_bin_sink: public spdlog::sinks::base_sink<std::mutex> {
public:
    explicit ostream_bin_sink(std::ostream& os, bool force_flush=false) :_ostream(os), _force_flush(force_flush) {}
    ostream_bin_sink(const ostream_bin_sink&) = delete;
    ostream_bin_sink& operator=(const ostream_bin_sink&) = delete;
    virtual ~ostream_bin_sink() = default;

protected:
    void _sink_it(const spdlog::details::log_msg& msg) override
    {
        //const char* bin_data = msg.raw.data();
        //size_t bin_size = msg.raw.size();
        //std::cerr << "msg.size " << msg.raw.size() << " ";
        _ostream.write(msg.raw.data(), msg.raw.size());
        //_ostream.write(msg.formatted.data(), msg.formatted.size());
        if (_force_flush)
            _ostream.flush();
    }

    void _flush() override
    {
        _ostream.flush();
    }

    std::ostream& _ostream;
    bool _force_flush;
};

typedef ostream_bin_sink<std::mutex> ostream_bin_sink_mt;
typedef ostream_bin_sink<spdlog::details::null_mutex> ostream_bin_sink_st;

#endif
