/**
 *    @file  kseq++.hpp
 *   @brief  C++ implementation of kseq library.
 *
 *  This is a header-only library re-implementing the original kseq library.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Jul 15, 2018  19:15
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  KSEQPP_KSEQPP_HPP__
#define  KSEQPP_KSEQPP_HPP__

#include <cassert>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <ios>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>

//#include "config.hpp"

namespace klibpp {
  template< typename TFile,
            typename TFunc,
            typename TSpec >
              class KStream;

  class KStreamBase_ {
    protected:
      /* Typedefs */
      using size_type = long int;
      using char_type = char;
  };

  struct KSeq {  // kseq_t
    std::string name;
    std::string comment;
    std::string seq;
    std::string qual;
    inline void clear( ) {
      name.clear();
      comment.clear();
      seq.clear();
      qual.clear();
    }
  };

  namespace mode {
    struct In_ { };
    struct Out_ { };

    constexpr In_ in;
    constexpr Out_ out;
  }  /* -----  end of namespace mode  ----- */

  namespace format {
    enum Format { mix, fasta, fastq };
  }

  struct KEnd_ {};
  constexpr KEnd_ kend;

  template< typename TFile,
            typename TFunc >
    class KStream< TFile, TFunc, mode::Out_ > : public KStreamBase_ {
      public:
        /* Typedefs */
        using base_type = KStreamBase_;
        using spec_type = mode::Out_;
        using size_type = base_type::size_type;
        using char_type = base_type::char_type;
        using close_type = int(*)( TFile );
      protected:
        /* Consts */
        constexpr static std::make_unsigned_t< size_type > DEFAULT_BUFSIZE = 131072;
        constexpr static unsigned int DEFAULT_WRAPLEN = 60;
        constexpr static format::Format DEFAULT_FORMAT = format::mix;
        /* Data members */
        char_type* m_buf;                               /**< @brief character buffer */
        char_type* w_buf;                               /**< @brief second character buffer */
        size_type bufsize;                              /**< @brief buffer size */
        std::thread worker;                             /**< @brief worker thread */
        std::unique_ptr< std::mutex > bufslock;         /**< @brief buffers mutex */
        std::unique_ptr< std::condition_variable > cv;  /**< @brief consumer/producer condition variable */
        bool terminate;                                 /**< @brief thread terminate flag XXX: set before notify */
        bool produced;                                  /**< @brief produced flag. XXX: SHARED (data race) */
        size_type m_begin;                              /**< @brief begin buffer index */
        size_type m_end;                                /**< @brief end buffer index or error flag if -1 */
        size_type w_end;                                /**< @brief end second buffer index or error flag if -1 */
        unsigned int wraplen;                           /**< @brief line wrap length */
        unsigned long int counter;                      /**< @brief number of records written so far */
        format::Format fmt;                             /**< @brief format of the output records */
        TFile f;                                        /**< @brief file handler */
        TFunc func;                                     /**< @brief write function */
        close_type close;                               /**< @brief close function */
      public:
        KStream( TFile f_,
            TFunc func_,
            spec_type=mode::out,
            format::Format fmt_=DEFAULT_FORMAT,
            std::make_unsigned_t< size_type > bs_=DEFAULT_BUFSIZE,
            close_type cfunc_=nullptr )
          : m_buf( new char_type[ bs_ ] ), w_buf( new char_type[ bs_ ] ),
          bufsize( bs_ ), bufslock( new std::mutex ), cv( new std::condition_variable ),
          wraplen( DEFAULT_WRAPLEN ), fmt( fmt_ ), f( std::move( f_ ) ),
          func( std::move(  func_  ) ), close( cfunc_ )
        {
          this->m_begin = 0;
          this->m_end = 0;
          this->terminate = false;
          this->produced = false;
          this->w_end = 0;
          this->counter = 0;
          this->worker_start();
        }

        KStream( TFile f_,
            TFunc func_,
            format::Format fmt_,
            std::make_unsigned_t< size_type > bs_=DEFAULT_BUFSIZE,
            close_type cfunc_=nullptr )
          : KStream( std::move( f_ ), std::move( func_ ), mode::out, fmt_, bs_, cfunc_ )
        { }

        KStream( TFile f_,
            TFunc func_,
            spec_type,
            format::Format fmt_,
            close_type cfunc_ )
          : KStream( std::move( f_ ), std::move( func_ ), mode::out, fmt_, DEFAULT_BUFSIZE, cfunc_ )
        { }

        KStream( TFile f_,
            TFunc func_,
            format::Format fmt_,
            close_type cfunc_ )
          : KStream( std::move( f_ ), std::move( func_ ), mode::out, fmt_, DEFAULT_BUFSIZE, cfunc_ )
        { }

        KStream( TFile f_,
            TFunc func_,
            spec_type,
            std::make_unsigned_t< size_type > bs_,
            close_type cfunc_=nullptr )
          : KStream( std::move( f_ ), std::move( func_ ), mode::out, DEFAULT_FORMAT, bs_, cfunc_ )
        { }

        KStream( TFile f_,
            TFunc func_,
            std::make_unsigned_t< size_type > bs_,
            close_type cfunc_=nullptr )
          : KStream( std::move( f_ ), std::move( func_ ), mode::out, DEFAULT_FORMAT, bs_, cfunc_ )
        { }

        KStream( TFile f_,
            TFunc func_,
            spec_type,
            close_type cfunc_ )
          : KStream( std::move( f_ ), std::move( func_ ), mode::out, DEFAULT_FORMAT, DEFAULT_BUFSIZE, cfunc_ )
        { }

        KStream( TFile f_,
            TFunc func_,
            close_type cfunc_ )
          : KStream( std::move( f_ ), std::move( func_ ), mode::out, DEFAULT_FORMAT, DEFAULT_BUFSIZE, cfunc_ )
        { }

        KStream( KStream const& ) = delete;
        KStream& operator=( KStream const& ) = delete;

        KStream( KStream&& other ) noexcept
        {
          other.worker_join();
          this->m_buf = other.m_buf;
          this->w_buf = other.w_buf;
          other.m_buf = nullptr;
          other.w_buf = nullptr;
          this->bufsize = other.bufsize;
          this->bufslock = std::move( other.bufslock );
          this->cv = std::move( other.cv );
          this->terminate = false;
          this->produced = other.produced;
          this->m_begin = other.m_begin;
          this->m_end = other.m_end;
          this->w_end = other.w_end;
          this->wraplen = other.wraplen;
          this->counter = other.counter;
          this->fmt = other.fmt;
          this->f = std::move( other.f );
          this->func = std::move( other.func );
          this->close = other.close;
          this->worker_start();
        }

        KStream& operator=( KStream&& other ) noexcept
        {
          if ( this == &other ) return *this;
          other.worker_join();
          delete[] this->m_buf;
          delete[] this->w_buf;
          this->m_buf = other.m_buf;
          this->w_buf = other.w_buf;
          other.m_buf = nullptr;
          other.w_buf = nullptr;
          this->bufsize = other.bufsize;
          this->bufslock = std::move( other.bufslock );
          this->cv = std::move( other.cv );
          this->terminate = false;
          this->produced = other.produced;
          this->m_begin = other.m_begin;
          this->m_end = other.m_end;
          this->w_end = other.w_end;
          this->wraplen = other.wraplen;
          this->counter = other.counter;
          this->fmt = other.fmt;
          this->f = std::move( other.f );
          this->func = std::move( other.func );
          this->close = other.close;
          this->worker_start();
          return *this;
        }

        ~KStream( ) noexcept
        {
          this->worker_join();
          delete[] this->m_buf;
          delete[] this->w_buf;
          if ( this->close != nullptr ) this->close( this->f );
        }
        /* Accessors */
          inline unsigned long int
        counts( ) const
        {
          return this->counter;
        }

          inline format::Format
        get_format( ) const
        {
          return this->fmt;
        }
        /* Mutators */
          inline void
        set_wraplen( unsigned int len )
        {
          this->wraplen = len;
        }

          inline void
        set_format( format::Format fmt_ )
        {
          this->fmt = fmt_;
        }
        /* Methods */
          inline bool
        fail( ) const
        {
          return this->m_end == -1;
        }

          inline KStream&
        operator<<( const KSeq& rec )
        {
          if ( ( this->fmt == format::mix && rec.qual.empty() ) ||  // FASTA record
               ( this->fmt == format::fasta ) ) this->puts( '>' );      // Forced FASTA
          else {
            if ( rec.qual.size() != rec.seq.size() ) {
              throw std::runtime_error( "the sequence length doesn't match with"
                                        " the length of its quality string.");
            }
            this->puts( '@' );  // FASTQ record
          }
          this->puts( rec.name );
          if ( !rec.comment.empty() ) {
            this->puts( ' ' );
            this->puts( rec.comment );
          }
          this->puts( '\n' );
          this->puts( rec.seq, true );
          if ( ( this->fmt == format::mix && !rec.qual.empty() ) ||  // FASTQ record
               ( this->fmt == format::fastq ) ) {                        // Forced FASTQ
            this->puts( '\n' );
            this->puts( '+' );
            this->puts( '\n' );
            this->puts( rec.qual, true );
          }
          this->puts( '\n' );
          if ( *this ) this->counter++;
          return *this;
        }

          inline KStream&
        operator<<( format::Format fmt_ )
        {
          this->fmt = fmt_;
          return *this;
        }

          inline KStream&
        operator<<( KEnd_ )
        {
          this->flush();
          return *this;
        }

        operator bool( ) const
        {
          return !this->fail();
        }
        /* Low-level methods */
          inline bool
        puts( std::string const& s, bool wrap=false ) noexcept
        {
          if ( this->fail() ) return false;

          std::string::size_type cursor = 0;
          std::string::size_type len = 0;
          while ( cursor != s.size() ) {
            assert( cursor < s.size() );
            if ( this->m_begin >= this->bufsize ) this->async_write();
            if ( this->fail() ) break;
            if ( wrap && cursor != 0 && cursor % this->wraplen == 0 ) {
              this->m_buf[ this->m_begin++ ] = '\n';
            }
            len = std::min( s.size() - cursor,
                static_cast< std::string::size_type >( this->bufsize - this->m_begin ) );
            if ( wrap )
              len = std::min( len, this->wraplen - cursor % this->wraplen );
            std::copy( &s[ cursor ], &s[ cursor ] + len, this->m_buf + this->m_begin );
            this->m_begin += len;
            cursor += len;
          }
          return !this->fail();
        }

          inline bool
        puts( char_type c ) noexcept
        {
          if ( this->fail() ) return false;
          if ( this->m_begin >= this->bufsize ) this->async_write();
          this->m_buf[ this->m_begin++ ] = c;
          return !this->fail();
        }

          inline void
        flush( ) noexcept
        {
          this->async_write( );
          {
            // wait until it is actually written to the file.
            std::unique_lock< std::mutex > lock( *this->bufslock );
            this->cv->wait( lock, [this]{ return !this->produced; } );
          }
        }
      private:
        /* Methods */
          inline void
        async_write( bool term=false ) noexcept
        {
          if ( this->fail() || this->terminate ) return;

          {
            std::unique_lock< std::mutex > lock( *this->bufslock );
            this->cv->wait( lock, [this]{ return !this->produced; } );
            this->m_end = this->w_end;
            if ( !this->fail() ) {
              this->w_end = this->m_begin;
              std::copy( this->m_buf, this->m_buf + this->m_begin, this->w_buf );
              this->produced = true;
              if ( term ) this->terminate = true;  /**< XXX: only set here! */
            }
            this->m_begin = 0;
          }
          this->cv->notify_one();
        }

          inline void
        writer( ) noexcept
        {
          bool term = false;
          do {
            {
              std::unique_lock< std::mutex > lock( *this->bufslock );
              this->cv->wait( lock, [this]{ return this->produced; } );
              if ( !this->func( this->f, this->w_buf, this->w_end ) && this->w_end ) {
                this->w_end = -1;
              }
              this->produced = false;
              if ( this->terminate || this->w_end < 0 ) term = true;
            }
            this->cv->notify_one();
          } while ( !term );
        }

          inline void
        worker_join( )
        {
          this->async_write( true );
          if ( this->worker.joinable() ) this->worker.join();
        }

          inline void
        worker_start( )
        {
          this->worker = std::thread( &KStream::writer, this );
        }
    };

  template< typename TFile,
            typename TFunc >
    class KStream< TFile, TFunc, mode::In_ > : public KStreamBase_ {  // kstream_t
      public:
        /* Typedefs */
        using base_type = KStreamBase_;
        using spec_type = mode::In_;
        using size_type = base_type::size_type;
        using char_type = base_type::char_type;
        using close_type = int(*)( TFile );
      protected:
        /* Separators */
        constexpr static char_type SEP_SPACE = 0;  // isspace(): \t, \n, \v, \f, \r
        constexpr static char_type SEP_TAB = 1;    // isspace() && !' '
        constexpr static char_type SEP_LINE = 2;   // line separator: "\n" (Unix) or "\r\n" (Windows)
        constexpr static char_type SEP_MAX = 2;
        /* Consts */
        constexpr static std::make_unsigned_t< size_type > DEFAULT_BUFSIZE = 16384;
        /* Data members */
        char_type* buf;                      /**< @brief character buffer */
        size_type bufsize;                   /**< @brief buffer size */
        size_type begin;                     /**< @brief begin buffer index */
        size_type end;                       /**< @brief end buffer index or error flag if -1 */
        bool is_eof;                         /**< @brief eof flag */
        bool is_tqs;                         /**< @brief truncated quality string flag */
        bool is_ready;                       /**< @brief next record ready flag */
        bool last;                           /**< @brief last read was successful */
        unsigned long int counter;           /**< @brief number of parsed records so far */
        TFile f;                             /**< @brief file handler */
        TFunc func;                          /**< @brief read function */
        close_type close;                    /**< @brief close function */
      public:
        KStream( TFile f_,
            TFunc func_,
            spec_type=mode::in,
            std::make_unsigned_t< size_type > bs_=DEFAULT_BUFSIZE,
            close_type cfunc_=nullptr )  // ks_init
          : buf( new char_type[ bs_ ] ), bufsize( bs_ ),
          f( std::move( f_ ) ), func( std::move(  func_  ) ), close( cfunc_ )
        {
          this->begin = 0;
          this->end = 0;
          this->is_eof = false;
          this->is_tqs = false;
          this->is_ready = false;
          this->last = false;
          this->counter = 0;
        }

        KStream( TFile f_,
            TFunc func_,
            std::make_unsigned_t< size_type > bs_,
            close_type cfunc_=nullptr )
          : KStream( std::move( f_ ), std::move( func_ ), mode::in, bs_, cfunc_ )
        { }

        KStream( TFile f_,
            TFunc func_,
            spec_type,
            close_type cfunc_ )
          : KStream( std::move( f_ ), std::move( func_ ), mode::in, DEFAULT_BUFSIZE, cfunc_ )
        { }

        KStream( TFile f_,
            TFunc func_,
            close_type cfunc_ )
          : KStream( std::move( f_ ), std::move( func_ ), mode::in, DEFAULT_BUFSIZE, cfunc_ )
        { }

        KStream( KStream const& ) = delete;
        KStream& operator=( KStream const& ) = delete;

        KStream( KStream&& other ) noexcept
        {
          this->buf = other.buf;
          other.buf = nullptr;
          this->bufsize = other.bufsize;
          this->begin = other.begin;
          this->end = other.end;
          this->is_eof = other.is_eof;
          this->is_tqs = other.is_tqs;
          this->is_ready = other.is_ready;
          this->last = other.last;
          this->counter = other.counter;
          this->f = std::move( other.f );
          this->func = std::move( other.func );
          this->close = other.close;
        }

        KStream& operator=( KStream&& other ) noexcept
        {
          if ( this == &other ) return *this;
          delete[] this->buf;
          this->buf = other.buf;
          other.buf = nullptr;
          this->bufsize = other.bufsize;
          this->begin = other.begin;
          this->end = other.end;
          this->is_eof = other.is_eof;
          this->is_tqs = other.is_tqs;
          this->is_ready = other.is_ready;
          this->last = other.last;
          this->counter = other.counter;
          this->f = std::move( other.f );
          this->func = std::move( other.func );
          this->close = other.close;
          return *this;
        }

        ~KStream( ) noexcept
        {
          delete[] this->buf;
          if ( this->close != nullptr ) this->close( this->f );
        }
        /* Accessors */
          inline unsigned long int
        counts( ) const
        {
          return this->counter;
        }
        /* Methods */
          inline bool
        err( ) const  // ks_err
        {
          return this->end == -1;
        }

          inline bool
        eof( ) const  // ks_eof
        {
          return this->is_eof && this->begin >= this->end;
        }

          inline bool
        tqs( ) const
        {
          return this->is_tqs;
        }

          inline bool
        fail( ) const
        {
          return this->err() || this->tqs() || ( this->eof() && !this->last );
        }

          inline KStream&
        operator>>( KSeq& rec )  // kseq_read
        {
          char_type c;
          this->last = false;
          if ( !this->is_ready ) {  // then jump to the next header line
            while ( ( c = this->getc( ) ) && c != '>' && c != '@' );
            if ( this->fail() ) return *this;
            this->is_ready = true;
          }  // else: the first header char has been read in the previous call
          rec.clear();  // reset all members
          if ( !this->getuntil( KStream::SEP_SPACE, rec.name, &c ) ) return *this;
          if ( c != '\n' ) {  // read FASTA/Q comment
            this->getuntil( KStream::SEP_LINE, rec.comment, nullptr );
          }
          while ( ( c = this->getc( ) ) && c != '>' && c != '@' && c != '+' ) {
            if ( c == '\n' ) continue;  // skip empty lines
            rec.seq += c;
            this->getuntil( KStream::SEP_LINE, rec.seq, nullptr, true ); // read the rest of the line
          }
          this->last = true;
          ++this->counter;
          if ( c == '>' || c == '@' ) this->is_ready = true;  // the first header char has been read
          if ( c != '+' ) return *this;  // FASTA
          while ( ( c = this->getc( ) ) && c != '\n' );  // skip the rest of '+' line
          if ( this->eof() ) {  // error: no quality string
            this->is_tqs = true;
            return *this;
          }
          while ( this->getuntil( KStream::SEP_LINE, rec.qual, nullptr, true ) &&
              rec.qual.size() < rec.seq.size() );
          if ( this->err() ) return *this;
          this->is_ready = false;  // we have not come to the next header line
          if ( rec.seq.size() != rec.qual.size() ) {  // error: qual string is of a different length
            this->is_tqs = true;  // should return here
          }

          return *this;
        }

        operator bool( ) const
        {
          return !this->fail();
        }

          inline std::vector< KSeq >
        read( std::vector< KSeq >::size_type const size )
        {
          std::vector< KSeq > ret;
          ret.reserve( size );
          for ( std::vector< KSeq >::size_type i = 0; i < size; ++i ) {
            ret.emplace_back();
            *this >> ret.back();
            if ( !( *this ) ) {
              ret.pop_back();
              break;
            }
          }
          return ret;
        }

          inline std::vector< KSeq >
        read( )
        {
          std::vector< KSeq > ret;
          while ( ( ret.emplace_back(), true ) && *this >> ret.back() );
          ret.pop_back();
          return ret;
        }
        /* Low-level methods */
          inline char_type
        getc( ) noexcept  // ks_getc
        {
          // error
          if ( this->err() || this->eof() ) return 0;
          // fetch
          if ( this->begin >= this->end ) {
            this->begin = 0;
            this->end = this->func( this->f, this->buf, this->bufsize );
            if ( this->end <= 0 ) {  // err if end == -1 and eof if 0
              this->is_eof = true;
              return 0;
            }
          }
          // ready
          return this->buf[ this->begin++ ];
        }

          inline bool
        getuntil( char_type delimiter, std::string& str, char_type *dret, bool append=false )  // ks_getuntil
          noexcept
        {
          char_type c;
          bool gotany = false;
          if ( dret ) *dret = 0;
          if ( !append ) str.clear();
          size_type i = -1;
          do {
            if ( !( c = this->getc( ) ) ) break;
            --this->begin;
            if ( delimiter == KStream::SEP_LINE ) {
              for ( i = this->begin; i < this->end; ++i ) {
                if ( this->buf[ i ] == '\n' ) break;
              }
            }
            else if ( delimiter > KStream::SEP_MAX ) {
              for ( i = this->begin; i < this->end; ++i ) {
                if ( this->buf[ i ] == delimiter ) break;
              }
            }
            else if ( delimiter == KStream::SEP_SPACE ) {
              for ( i = this->begin; i < this->end; ++i ) {
                if ( std::isspace( this->buf[ i ] ) ) break;
              }
            }
            else if ( delimiter == KStream::SEP_TAB ) {
              for ( i = this->begin; i < this->end; ++i ) {
                if ( std::isspace( this->buf[ i ] ) && this->buf[ i ] != ' ' ) break;
              }
            }
            else {
              assert( false );  // it should not reach here
              return false;  // when assert is replaced by NOOP
            }

            gotany = true;
            str.append( this->buf + this->begin, i - this->begin );
            this->begin = i + 1;
          } while ( i >= this->end );

          if ( this->err() || ( this->eof() && !gotany ) ) return false;

          assert( i != -1 );
          if ( !this->eof() && dret ) *dret = this->buf[ i ];
          if ( delimiter == KStream::SEP_LINE && !str.empty() && str.back() == '\r' ) {
            str.pop_back();
          }
          return true;
        }
    };

  template< typename TFile, typename TFunc >
    using KStreamIn = KStream< TFile, TFunc, mode::In_ >;

  template< typename TFile, typename TFunc >
    using KStreamOut = KStream< TFile, TFunc, mode::Out_ >;

  template< typename TFile, typename TFunc, typename TSpec, typename... Args >
      inline KStream< std::decay_t< TFile >, std::decay_t< TFunc >, TSpec >
    make_kstream( TFile&& file, TFunc&& func, TSpec, Args&&... args )
    {
      return KStream< std::decay_t< TFile >, std::decay_t< TFunc >, TSpec >(
          std::forward< TFile >( file ), std::forward< TFunc >( func ), TSpec(),
          std::forward< Args >( args )... );
    }

  template< typename TFile, typename TFunc, typename... Args >
      inline KStream< std::decay_t< TFile >, std::decay_t< TFunc >, mode::In_ >
    make_ikstream( TFile&& file, TFunc&& func, Args&&... args )
    {
      return KStream< std::decay_t< TFile >, std::decay_t< TFunc >, mode::In_ >(
          std::forward< TFile >( file ), std::forward< TFunc >( func ), mode::in,
          std::forward< Args >( args )... );
    }

  template< typename TFile, typename TFunc, typename... Args >
      inline KStream< std::decay_t< TFile >, std::decay_t< TFunc >, mode::Out_ >
    make_okstream( TFile&& file, TFunc&& func, Args&&... args )
    {
      return KStream< std::decay_t< TFile >, std::decay_t< TFunc >, mode::Out_ >(
          std::forward< TFile >( file ), std::forward< TFunc >( func ), mode::out,
          std::forward< Args >( args )... );
    }
}  /* -----  end of namespace klibpp  ----- */
#endif  /* ----- #ifndef KSEQPP_KSEQPP_HPP__  ----- */
