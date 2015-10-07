#ifndef __BOOTSTRAP_WRITER_HPP__
#define __BOOTSTRAP_WRITER_HPP__

#include <vector>
#include <mutex>

class BootstrapWriter {
    public:
        virtual ~BootstrapWriter() {}
        virtual bool writeHeader(std::string& comments, std::vector<Transcript>& transcripts) = 0;
        virtual bool writeBootstrap(std::vector<double>& abund) = 0;
};

#endif // __BOOTSTRAP_WRITER_HPP__

