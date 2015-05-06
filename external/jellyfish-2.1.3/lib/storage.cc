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

#include "config.h"
#include <jellyfish/storage.hpp>

namespace jellyfish {
const size_t _quadratic_reprobes[257] = {
      1,
      1,     3,     6,    10,    15,    21,    28,    36,    45,    55,
     66,    78,    91,   105,   120,   136,   153,   171,   190,   210,
    231,   253,   276,   300,   325,   351,   378,   406,   435,   465,
    496,   528,   561,   595,   630,   666,   703,   741,   780,   820,
    861,   903,   946,   990,  1035,  1081,  1128,  1176,  1225,  1275,
   1326,  1378,  1431,  1485,  1540,  1596,  1653,  1711,  1770,  1830,
   1891,  1953,  2016,  2080,  2145,  2211,  2278,  2346,  2415,  2485,
   2556,  2628,  2701,  2775,  2850,  2926,  3003,  3081,  3160,  3240,
   3321,  3403,  3486,  3570,  3655,  3741,  3828,  3916,  4005,  4095,
   4186,  4278,  4371,  4465,  4560,  4656,  4753,  4851,  4950,  5050,
   5151,  5253,  5356,  5460,  5565,  5671,  5778,  5886,  5995,  6105,
   6216,  6328,  6441,  6555,  6670,  6786,  6903,  7021,  7140,  7260,
   7381,  7503,  7626,  7750,  7875,  8001,  8128,  8256,  8385,  8515,
   8646,  8778,  8911,  9045,  9180,  9316,  9453,  9591,  9730,  9870,
  10011, 10153, 10296, 10440, 10585, 10731, 10878, 11026, 11175, 11325,
  11476, 11628, 11781, 11935, 12090, 12246, 12403, 12561, 12720, 12880,
  13041, 13203, 13366, 13530, 13695, 13861, 14028, 14196, 14365, 14535,
  14706, 14878, 15051, 15225, 15400, 15576, 15753, 15931, 16110, 16290,
  16471, 16653, 16836, 17020, 17205, 17391, 17578, 17766, 17955, 18145,
  18336, 18528, 18721, 18915, 19110, 19306, 19503, 19701, 19900, 20100,
  20301, 20503, 20706, 20910, 21115, 21321, 21528, 21736, 21945, 22155,
  22366, 22578, 22791, 23005, 23220, 23436, 23653, 23871, 24090, 24310,
  24531, 24753, 24976, 25200, 25425, 25651, 25878, 26106, 26335, 26565,
  26796, 27028, 27261, 27495, 27730, 27966, 28203, 28441, 28680, 28920,
  29161, 29403, 29646, 29890, 30135, 30381, 30628, 30876, 31125, 31375,
  31626, 31878, 32131, 32385, 32640, 32896
};
const size_t *quadratic_reprobes = _quadratic_reprobes;
}
