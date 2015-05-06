/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __JELLYFISH_MERGE_FILES_HPP__
#define __JELLYFISH_MERGE_FILES_HPP__

#include <vector>
#include <jellyfish/err.hpp>
#include <jellyfish/file_header.hpp>

define_error_class(MergeError);

/// Merge files. Throw a MergeError in case of error.
void merge_files(std::vector<const char*> input_files, const char* out_file,
                 jellyfish::file_header& h, uint64_t min, uint64_t max);

#endif /* __JELLYFISH_MERGE_FILES_HPP__ */
