#include "salmon/internal/io/AlignmentIO.hpp"

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

void appendTag(salmon::io::AlignmentHeaderRef& ref, const std::string& token) {
  auto* tag = new AlignmentHeaderTag;
  tag->str = dupCString(token);
  tag->len = static_cast<int>(token.size());
  tag->next = ref.tag;
  ref.tag = tag;
}

void populateReferenceTags(AlignmentHeader* header) {
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

AlignmentHeader* wrapHeader(sam_hdr_t* raw) {
  if (raw == nullptr) {
    return nullptr;
  }

  auto* header = new AlignmentHeader;
  header->raw = raw;
  header->nref = raw->n_targets;
  if (header->nref > 0) {
    header->ref = new AlignmentHeaderRef[header->nref];
    for (int32_t i = 0; i < header->nref; ++i) {
      header->ref[i].name = dupCString(raw->target_name[i]);
      header->ref[i].len = raw->target_len[i];
      header->ref[i].tag = nullptr;
    }
  }
  populateReferenceTags(header);
  return header;
}

void destroyHeader(AlignmentHeader* header) {
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

AlignmentFileHandle* openAlignmentFile(const char* path, const char* mode) {
  auto* raw = sam_open(path, mode);
  if (raw == nullptr) {
    return nullptr;
  }

  auto* file = new AlignmentFileHandle;
  file->raw = raw;

  if (mode != nullptr && mode[0] == 'r') {
    auto* rawHeader = sam_hdr_read(raw);
    if (rawHeader == nullptr) {
      sam_close(raw);
      delete file;
      return nullptr;
    }
    file->header = salmon::io::wrapHeader(rawHeader);
  }

  return file;
}

int closeAlignmentFile(AlignmentFileHandle* file) {
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
