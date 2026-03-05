#include "StadenUtils.hpp"

#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

char* dupCString(const char* str) {
  if (str == nullptr) {
    return nullptr;
  }

  const auto len = std::strlen(str);
  auto* out = static_cast<char*>(std::malloc(len + 1));
  if (out != nullptr) {
    std::memcpy(out, str, len + 1);
  }
  return out;
}

char* dupCString(const std::string& str) { return dupCString(str.c_str()); }

void appendTag(salmon::io::SAM_hdr_ref& ref, const std::string& token) {
  auto* tag = new SAM_hdr_tag;
  tag->str = dupCString(token);
  tag->len = static_cast<int>(token.size());
  tag->next = ref.tag;
  ref.tag = tag;
}

void populateReferenceTags(SAM_hdr* header) {
  if (header == nullptr || header->raw == nullptr) {
    return;
  }

  const char* text = sam_hdr_str(header->raw);
  if (text == nullptr) {
    return;
  }

  std::unordered_map<std::string, int32_t> refNameToIndex;
  for (int32_t i = 0; i < header->nref; ++i) {
    refNameToIndex.emplace(header->ref[i].name, i);
  }

  std::string headerText{text};
  size_t start = 0;
  while (start < headerText.size()) {
    size_t end = headerText.find('\n', start);
    if (end == std::string::npos) {
      end = headerText.size();
    }

    std::string line = headerText.substr(start, end - start);
    start = end + 1;

    if (line.rfind("@SQ\t", 0) != 0) {
      continue;
    }

    std::string refName;
    std::vector<std::string> tags;
    size_t tokenStart = 4;
    while (tokenStart <= line.size()) {
      size_t tokenEnd = line.find('\t', tokenStart);
      if (tokenEnd == std::string::npos) {
        tokenEnd = line.size();
      }
      auto token = line.substr(tokenStart, tokenEnd - tokenStart);
      if (token.rfind("SN:", 0) == 0) {
        refName = token.substr(3);
      } else if (!token.empty()) {
        tags.push_back(token);
      }
      tokenStart = tokenEnd + 1;
    }

    auto refIt = refNameToIndex.find(refName);
    if (refIt == refNameToIndex.end()) {
      continue;
    }

    auto& ref = header->ref[refIt->second];
    for (auto tagIt = tags.rbegin(); tagIt != tags.rend(); ++tagIt) {
      appendTag(ref, *tagIt);
    }
  }
}

} // namespace

namespace salmon::io {

SAM_hdr* wrap_header(sam_hdr_t* raw) {
  if (raw == nullptr) {
    return nullptr;
  }

  auto* header = new SAM_hdr;
  header->raw = raw;
  header->nref = raw->n_targets;
  if (header->nref > 0) {
    header->ref = new SAM_hdr_ref[header->nref];
    for (int32_t i = 0; i < header->nref; ++i) {
      header->ref[i].name = dupCString(raw->target_name[i]);
      header->ref[i].len = raw->target_len[i];
      header->ref[i].tag = nullptr;
    }
  }
  populateReferenceTags(header);
  return header;
}

void destroy_header(SAM_hdr* header) {
  if (header == nullptr) {
    return;
  }

  if (header->ref != nullptr) {
    for (int32_t i = 0; i < header->nref; ++i) {
      std::free(header->ref[i].name);
      auto* tag = header->ref[i].tag;
      while (tag != nullptr) {
        auto* next = tag->next;
        std::free(tag->str);
        delete tag;
        tag = next;
      }
    }
    delete[] header->ref;
  }

  if (header->raw != nullptr) {
    sam_hdr_destroy(header->raw);
  }
  delete header;
}

} // namespace salmon::io

namespace staden {
namespace utils {

bam_seq_t* bam_init() {
  return bam_init1();
}

void bam_destroy(bam_seq_t* b) {
  if (b == 0) {
    return;
  }
  bam_destroy1(b);
}

} // namespace utils
} // namespace staden

scram_fd* scram_open(const char* path, const char* mode) {
  auto* raw = sam_open(path, mode);
  if (raw == nullptr) {
    return nullptr;
  }

  auto* file = new scram_fd;
  file->raw = raw;

  if (mode != nullptr && mode[0] == 'r') {
    auto* rawHeader = sam_hdr_read(raw);
    if (rawHeader == nullptr) {
      sam_close(raw);
      delete file;
      return nullptr;
    }
    file->header = salmon::io::wrap_header(rawHeader);
  }

  return file;
}

int scram_close(scram_fd* file) {
  if (file == nullptr) {
    return 0;
  }

  int rc = 0;
  if (file->raw != nullptr) {
    rc = sam_close(file->raw);
  }
  delete file;
  return rc;
}

SAM_hdr* scram_get_header(scram_fd* file) {
  return (file == nullptr) ? nullptr : file->header;
}

int scram_set_option(scram_fd* file, int option, int value) {
  if (file == nullptr || file->raw == nullptr) {
    return -1;
  }
  if (option == CRAM_OPT_NTHREADS) {
    return hts_set_threads(file->raw, value);
  }
  return 0;
}

int scram_get_seq(scram_fd* file, bam_seq_t** read) {
  if (file == nullptr || file->raw == nullptr || file->header == nullptr) {
    return -1;
  }
  if (*read == nullptr) {
    *read = bam_init1();
  }
  return sam_read1(file->raw, file->header->raw, *read);
}

int scram_put_seq(scram_fd* file, bam_seq_t* read) {
  if (file == nullptr || file->raw == nullptr || file->header == nullptr) {
    return -1;
  }
  return (sam_write1(file->raw, file->header->raw, read) >= 0) ? 0 : -1;
}

void scram_set_header(scram_fd* file, SAM_hdr* header) {
  if (file != nullptr) {
    file->header = header;
  }
}

int scram_write_header(scram_fd* file) {
  if (file == nullptr || file->raw == nullptr || file->header == nullptr) {
    return -1;
  }
  return sam_hdr_write(file->raw, file->header->raw);
}
