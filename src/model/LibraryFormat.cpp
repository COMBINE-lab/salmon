#include "LibraryFormat.hpp"

LibraryFormat::LibraryFormat(ReadType type_in, ReadOrientation orientation_in,
                             ReadStrandedness strandedness_in)
    : type(type_in), orientation(orientation_in),
      strandedness(strandedness_in) {}

bool LibraryFormat::check() {
  bool valid = true;
  switch (type) {

  case ReadType::SINGLE_END:
    valid = valid and (orientation == ReadOrientation::NONE);
    switch (strandedness) {
    // If we're dealing with single end reads, then separate orientations
    // are nonsensical.
    case ReadStrandedness::SA:
    case ReadStrandedness::AS:
      valid = false;
      break;
    default:
      break;
    }
    break;

  case ReadType::PAIRED_END:
    switch (orientation) {
    // For paired-end reads, an orientation must be defined
    case ReadOrientation::NONE:
      valid = false;
      break;
    // If the reads are oriented in the same direction, then they
    // have to come from the same strand.
    case ReadOrientation::SAME:
      valid = valid and ((strandedness == ReadStrandedness::S) or
                         (strandedness == ReadStrandedness::A) or
                         (strandedness == ReadStrandedness::U));

      break;
    // If the reads are oriented away from each other or toward each other,
    // then they must come from different strands.
    case ReadOrientation::AWAY:
    case ReadOrientation::TOWARD:
      valid = valid and ((strandedness == ReadStrandedness::SA) or
                         (strandedness == ReadStrandedness::AS) or
                         (strandedness == ReadStrandedness::U));
      break;
    }
    break;
  }
  return valid;
}

std::ostream& operator<<(std::ostream& os, const LibraryFormat& lf) {
  os << "Library format { type:";
  switch (lf.type) {
  case ReadType::SINGLE_END:
    os << "single end";
    break;
  case ReadType::PAIRED_END:
    os << "paired end";
    break;
  }

  os << ", relative orientation:";
  switch (lf.orientation) {
  case ReadOrientation::TOWARD:
    os << "inward";
    break;
  case ReadOrientation::AWAY:
    os << "outward";
    break;
  case ReadOrientation::SAME:
    os << "matching";
    break;
  case ReadOrientation::NONE:
    os << "none";
    break;
  }

  os << ", strandedness:";
  switch (lf.strandedness) {
  case ReadStrandedness::SA:
    os << "(sense, antisense)";
    break;
  case ReadStrandedness::AS:
    os << "(antisense, sense)";
    break;
  case ReadStrandedness::S:
    os << "sense";
    break;
  case ReadStrandedness::A:
    os << "antisense";
    break;
  case ReadStrandedness::U:
    os << "unstranded";
    break;
  }
  os << " }";
  return os;
}
