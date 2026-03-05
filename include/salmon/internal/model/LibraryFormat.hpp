#ifndef LIBRARY_FORMAT_HPP
#define LIBRARY_FORMAT_HPP

#include <cstdint>
#include <ostream>

enum class ReadType : std::uint8_t { SINGLE_END = 0, PAIRED_END = 1 };
enum class ReadOrientation : std::uint8_t {
  SAME = 0,
  AWAY = 1,
  TOWARD = 2,
  NONE = 3
};
enum class ReadStrandedness : std::uint8_t {
  SA = 0,
  AS = 1,
  S = 2,
  A = 3,
  U = 4
};

class LibraryFormat {
public:
  LibraryFormat(ReadType type_in, ReadOrientation orientation_in,
                ReadStrandedness strandedness_in);

  LibraryFormat(const LibraryFormat& o) = default;
  LibraryFormat(LibraryFormat&& o) = default;
  LibraryFormat& operator=(LibraryFormat& o) = default;
  LibraryFormat& operator=(LibraryFormat&& o) = default;

  ReadType type;
  ReadOrientation orientation;
  ReadStrandedness strandedness;

  /**
   * Returns true if the specified library format is OK, false otherwise.
   */
  bool check();
  friend std::ostream& operator<<(std::ostream& os, const LibraryFormat& lf);
  // The maximum index that could be returned for a library type
  constexpr static uint8_t maxLibTypeID() {
    return 0 | ((static_cast<uint8_t>(ReadType::PAIRED_END)) |
                (static_cast<uint8_t>(ReadOrientation::NONE)) << 1 |
                (static_cast<uint8_t>(ReadStrandedness::U)) << 3);
  }

  inline static LibraryFormat formatFromID(uint8_t id) {
    ReadType rt = ReadType::SINGLE_END;
    ReadOrientation ro = ReadOrientation::NONE;
    ReadStrandedness rs = ReadStrandedness::U;

    switch (id & 0x01) {
    case 0:
      rt = ReadType::SINGLE_END;
      break;
    case 1:
      rt = ReadType::PAIRED_END;
      break;
    }

    switch ((id >> 1) & 0x3) {
    case static_cast<uint8_t>(ReadOrientation::SAME):
      ro = ReadOrientation::SAME;
      break;
    case static_cast<uint8_t>(ReadOrientation::AWAY):
      ro = ReadOrientation::AWAY;
      break;
    case static_cast<uint8_t>(ReadOrientation::TOWARD):
      ro = ReadOrientation::TOWARD;
      break;
    case static_cast<uint8_t>(ReadOrientation::NONE):
      ro = ReadOrientation::NONE;
      break;
    }

    switch ((id >> 3) & 0x7) {
    case static_cast<uint8_t>(ReadStrandedness::SA):
      rs = ReadStrandedness::SA;
      break;
    case static_cast<uint8_t>(ReadStrandedness::AS):
      rs = ReadStrandedness::AS;
      break;
    case static_cast<uint8_t>(ReadStrandedness::S):
      rs = ReadStrandedness::S;
      break;
    case static_cast<uint8_t>(ReadStrandedness::A):
      rs = ReadStrandedness::A;
      break;
    case static_cast<uint8_t>(ReadStrandedness::U):
      rs = ReadStrandedness::U;
      break;
    }

    return LibraryFormat(rt, ro, rs);
  }

  std::string toString() const {
    std::string desc;
    if (type == ReadType::SINGLE_END) {
      if (orientation != ReadOrientation::NONE) {
        return desc;
      } else if (strandedness == ReadStrandedness::SA or
                 strandedness == ReadStrandedness::AS) {
        return desc;
      }
    }
    if (type == ReadType::PAIRED_END) {
      if (orientation == ReadOrientation::NONE) {
        return desc;
      } else if (strandedness == ReadStrandedness::S or
                 strandedness == ReadStrandedness::A) {
        return desc;
      }
    }

    /*
    switch (type) {
    case ReadType::SINGLE_END:
        desc += "single:";
        break;
    case ReadType::PAIRED_END:
        desc += "paired:";
        break;
    }
    */

    switch (orientation) {
    case ReadOrientation::TOWARD:
      desc += "I";
      break;
    case ReadOrientation::AWAY:
      desc += "O";
      break;
    case ReadOrientation::SAME:
      desc += "M";
      break;
    case ReadOrientation::NONE:
      break;
    }

    switch (strandedness) {
    case ReadStrandedness::SA:
      desc += "SF";
      break;
    case ReadStrandedness::AS:
      desc += "SR";
      break;
    case ReadStrandedness::S:
      desc += "SF";
      break;
    case ReadStrandedness::A:
      desc += "SR";
      break;
    case ReadStrandedness::U:
      desc += "U";
      break;
    }
    return desc;
  }

  // Assigns a unique ID to each potential library
  // type.  The IDs are such that 0 <= formatID(lib) < num possible formats
  inline uint8_t formatID() const {
    uint8_t id = 0;
    // mask: 00000001
    id |= (0x01 & static_cast<uint8_t>(type));
    // mask: 00000110
    id |= (0x3 & static_cast<uint8_t>(orientation)) << 1;
    // mask: 00111000
    id |= (0x7 & static_cast<uint8_t>(strandedness)) << 3;
    return id;
  }
};

inline bool operator==(const LibraryFormat& lhs, const LibraryFormat& rhs) {
  return ((lhs.type == rhs.type) and (lhs.orientation == rhs.orientation) and
          (lhs.strandedness == rhs.strandedness));
}

#endif // LIBRARY_FORMAT_HPP
