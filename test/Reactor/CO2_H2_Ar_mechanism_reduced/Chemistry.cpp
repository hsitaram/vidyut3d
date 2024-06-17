#include "Chemistry.H"
const int rmap[NUM_REACTIONS] = {
    12,  13,  16,  19,  24,  41,  43,  44,  50,  51,  52,  56,  84,  86,  87,
    138, 140, 143, 148, 151, 154, 158, 159, 164, 167, 173, 179, 184, 192, 195,
    196, 199, 202, 203, 261, 263, 264, 270, 274, 275, 276, 277, 278, 279, 280,
    281, 283, 284, 285, 286, 287, 288, 289, 291, 292, 293, 294, 295, 297, 300,
    301, 308, 324, 334, 374, 376, 378, 0,   1,   2,   3,   4,   5,   6,   7,
    8,   9,   10,  11,  14,  15,  17,  18,  20,  21,  22,  23,  25,  26,  27,
    28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  42,  45,
    46,  47,  48,  49,  53,  54,  55,  57,  58,  59,  60,  61,  62,  63,  64,
    65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
    80,  81,  82,  83,  85,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,
    98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112,
    113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 139, 141, 142, 144, 145,
    146, 147, 149, 150, 152, 153, 155, 156, 157, 160, 161, 162, 163, 165, 166,
    168, 169, 170, 171, 172, 174, 175, 176, 177, 178, 180, 181, 182, 183, 185,
    186, 187, 188, 189, 190, 191, 193, 194, 197, 198, 200, 201, 204, 205, 206,
    207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221,
    222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236,
    237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251,
    252, 253, 254, 255, 256, 257, 258, 259, 260, 262, 265, 266, 267, 268, 269,
    271, 272, 273, 282, 290, 296, 298, 299, 302, 303, 304, 305, 306, 307, 309,
    310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 325,
    326, 327, 328, 329, 330, 331, 332, 333, 335, 336, 337, 338, 339, 340, 341,
    342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356,
    357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371,
    372, 373, 375, 377, 379, 380, 381, 382, 383, 384};

// Returns 0-based map of reaction order
void GET_RMAP(int *_rmap) {
  for (int j = 0; j < NUM_REACTIONS; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(const int i, int &nspec, int ki[], int nu[]) {
  const int ns[NUM_GAS_REACTIONS] = {
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 2,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 4, 4, 4, 4, 4,
      3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 5, 4, 4, 4, 4, 4, 4,
      5, 4, 4, 4, 4, 4, 5, 5, 4, 3, 3, 2, 3, 5, 5, 4, 5, 5, 4, 4, 2, 2, 3, 4, 4,
      4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      5, 5, 5, 5, 5, 4, 4, 4, 3, 4, 4, 5, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4, 4, 4, 5,
      4, 5, 4, 4, 4, 4, 5, 4, 5, 4, 5, 5, 5, 5, 4, 4, 3, 4, 4, 4, 4, 3, 4, 4, 4,
      4, 5, 4, 5, 4, 5, 4, 5, 4, 5, 5, 5, 6, 6, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 2,
      2, 2, 2, 3, 3, 3, 3, 4, 2, 2, 2, 2, 2, 2, 3, 4, 3, 3, 3, 3, 3, 4, 2, 3, 4,
      3, 3, 4, 3, 3, 3, 4, 4, 3, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4,
      3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 3, 3,
      3, 3, 4, 3, 4, 5, 5, 4, 4, 5};
  const int kiv[NUM_GAS_REACTIONS * 6] = {
      12, 0,  12, 0,  0,  0, 12, 0,  13, 0,  0,  0, 12, 0,  14, 0,  0,  0,
      12, 0,  14, 0,  0,  0, 12, 0,  15, 0,  0,  0, 12, 0,  16, 0,  0,  0,
      12, 0,  16, 0,  0,  0, 12, 0,  16, 0,  0,  0, 12, 0,  16, 0,  0,  0,
      12, 0,  16, 0,  0,  0, 12, 0,  16, 0,  0,  0, 12, 0,  16, 0,  0,  0,
      12, 11, 7,  0,  0,  0, 12, 17, 0,  0,  0,  0, 12, 0,  11, 8,  0,  0,
      11, 0,  11, 0,  0,  0, 11, 21, 7,  0,  0,  0, 38, 0,  38, 0,  0,  0,
      38, 0,  39, 0,  0,  0, 38, 40, 0,  0,  0,  0, 0,  3,  0,  3,  0,  0,
      0,  3,  0,  4,  0,  0, 0,  3,  0,  4,  0,  0, 0,  3,  0,  4,  0,  0,
      3,  1,  0,  0,  0,  0, 0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0,
      0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0,
      0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0,
      0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0,
      0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0,
      0,  3,  0,  3,  0,  0, 0,  3,  0,  3,  0,  0, 3,  0,  5,  0,  0,  0,
      0,  3,  1,  2,  0,  0, 0,  9,  10, 0,  0,  0, 0,  7,  8,  0,  0,  0,
      7,  8,  0,  9,  0,  0, 1,  8,  0,  22, 0,  0, 3,  8,  0,  24, 0,  0,
      21, 8,  11, 0,  0,  0, 11, 8,  12, 0,  0,  0, 10, 0,  9,  0,  0,  0,
      10, 0,  9,  0,  0,  0, 10, 0,  9,  0,  0,  0, 1,  2,  0,  3,  0,  0,
      2,  7,  0,  22, 0,  0, 2,  22, 0,  24, 0,  0, 2,  0,  1,  0,  0,  0,
      31, 2,  33, 0,  0,  0, 29, 2,  31, 0,  0,  0, 11, 20, 12, 0,  0,  0,
      17, 0,  11, 7,  0,  0, 17, 0,  21, 9,  0,  0, 18, 0,  11, 12, 0,  0,
      19, 0,  12, 0,  0,  0, 0,  5,  1,  0,  0,  0, 0,  6,  1,  0,  0,  0,
      0,  6,  1,  3,  0,  0, 32, 0,  29, 1,  0,  0, 32, 0,  21, 1,  3,  0,
      30, 0,  21, 3,  0,  0, 30, 0,  21, 1,  0,  0, 28, 0,  21, 1,  0,  0,
      0,  26, 1,  24, 0,  0, 0,  26, 3,  22, 0,  0, 0,  26, 1,  22, 0,  0,
      0,  26, 1,  3,  7,  0, 0,  25, 1,  22, 0,  0, 0,  25, 3,  7,  0,  0,
      0,  25, 1,  7,  0,  0, 0,  23, 1,  7,  0,  0, 41, 0,  38, 1,  0,  0,
      42, 39, 38, 40, 0,  0, 42, 39, 38, 43, 0,  0, 42, 38, 43, 0,  0,  0,
      38, 39, 42, 0,  0,  0, 38, 39, 38, 0,  0,  0, 42, 38, 0,  0,  0,  0,
      42, 38, 39, 0,  0,  0, 39, 12, 38, 11, 7,  0, 39, 11, 38, 21, 7,  0,
      39, 3,  38, 1,  0,  0, 39, 22, 38, 1,  7,  0, 39, 24, 38, 1,  22, 0,
      39, 9,  38, 7,  0,  0, 42, 9,  38, 7,  0,  0, 39, 38, 0,  0,  0,  0,
      42, 38, 0,  0,  0,  0, 12, 13, 12, 0,  0,  0, 11, 13, 11, 12, 0,  0,
      13, 9,  12, 9,  0,  0, 38, 13, 38, 12, 0,  0, 13, 3,  12, 3,  0,  0,
      12, 14, 12, 0,  0,  0, 11, 14, 11, 12, 0,  0, 14, 9,  12, 9,  0,  0,
      38, 14, 38, 12, 0,  0, 14, 3,  12, 3,  0,  0, 12, 14, 12, 13, 0,  0,
      11, 14, 11, 13, 0,  0, 14, 9,  13, 9,  0,  0, 38, 14, 38, 13, 0,  0,
      14, 3,  13, 3,  0,  0, 12, 15, 12, 14, 0,  0, 11, 15, 11, 14, 0,  0,
      15, 9,  14, 9,  0,  0, 38, 15, 38, 14, 0,  0, 15, 3,  14, 3,  0,  0,
      12, 15, 12, 16, 0,  0, 11, 15, 11, 16, 0,  0, 15, 9,  16, 9,  0,  0,
      38, 15, 38, 16, 0,  0, 15, 3,  16, 3,  0,  0, 12, 15, 13, 14, 0,  0,
      12, 15, 12, 13, 0,  0, 11, 15, 11, 13, 0,  0, 15, 9,  13, 9,  0,  0,
      38, 15, 38, 13, 0,  0, 15, 3,  13, 3,  0,  0, 12, 16, 12, 14, 0,  0,
      11, 16, 11, 14, 0,  0, 16, 9,  14, 9,  0,  0, 38, 16, 38, 14, 0,  0,
      16, 3,  14, 3,  0,  0, 12, 16, 12, 13, 0,  0, 11, 16, 11, 13, 0,  0,
      16, 9,  13, 9,  0,  0, 38, 16, 38, 13, 0,  0, 16, 3,  13, 3,  0,  0,
      4,  3,  0,  0,  0,  0, 40, 8,  38, 7,  0,  0, 40, 8,  38, 7,  0,  0,
      40, 10, 38, 9,  0,  0, 40, 10, 38, 7,  0,  0, 40, 10, 38, 9,  0,  0,
      43, 8,  38, 7,  0,  0, 43, 10, 38, 9,  0,  0, 43, 10, 38, 7,  0,  0,
      40, 2,  38, 1,  0,  0, 40, 2,  38, 1,  0,  0, 43, 2,  38, 1,  0,  0,
      41, 8,  38, 1,  7,  0, 41, 8,  38, 1,  7,  0, 41, 10, 38, 1,  9,  0,
      41, 10, 38, 1,  7,  0, 41, 10, 38, 1,  9,  0, 41, 2,  38, 3,  0,  0,
      41, 2,  38, 1,  0,  0, 5,  8,  1,  7,  0,  0, 5,  8,  24, 0,  0,  0,
      5,  8,  3,  7,  0,  0, 5,  8,  3,  7,  0,  0, 6,  8,  1,  3,  7,  0,
      8,  23, 7,  22, 0,  0, 8,  23, 1,  7,  0,  0, 8,  23, 7,  22, 0,  0,
      25, 8,  24, 7,  0,  0, 25, 8,  1,  7,  22, 0, 25, 8,  24, 7,  0,  0,
      26, 8,  1,  24, 7,  0, 5,  10, 3,  9,  0,  0, 5,  10, 1,  9,  0,  0,
      5,  10, 3,  7,  0,  0, 5,  10, 1,  7,  0,  0, 5,  10, 3,  9,  0,  0,
      6,  10, 1,  3,  9,  0, 10, 23, 9,  22, 0,  0, 10, 23, 1,  7,  9,  0,
      10, 23, 7,  22, 0,  0, 10, 23, 1,  7,  0,  0, 10, 23, 9,  22, 0,  0,
      25, 10, 24, 9,  0,  0, 25, 10, 1,  9,  22, 0, 25, 10, 24, 7,  0,  0,
      25, 10, 1,  7,  22, 0, 25, 10, 24, 9,  0,  0, 26, 10, 1,  24, 9,  0,
      26, 10, 1,  24, 7,  0, 17, 10, 11, 7,  9,  0, 18, 10, 11, 12, 9,  0,
      19, 10, 12, 9,  0,  0, 2,  5,  1,  3,  0,  0, 2,  5,  1,  0,  0,  0,
      2,  5,  1,  3,  0,  0, 2,  6,  1,  3,  0,  0, 2,  23, 1,  7,  0,  0,
      2,  23, 1,  22, 0,  0, 2,  23, 24, 0,  0,  0, 2,  25, 1,  24, 0,  0,
      2,  25, 1,  22, 0,  0, 2,  25, 1,  24, 0,  0, 2,  26, 3,  24, 0,  0,
      2,  26, 1,  3,  22, 0, 2,  26, 3,  24, 0,  0, 2,  26, 1,  3,  22, 0,
      17, 2,  12, 1,  0,  0, 17, 2,  11, 1,  7,  0, 17, 20, 12, 7,  0,  0,
      18, 20, 11, 12, 7,  0, 19, 20, 12, 7,  0,  0, 20, 23, 12, 7,  22, 0,
      20, 23, 12, 1,  7,  0, 20, 25, 12, 24, 7,  0, 20, 25, 12, 1,  7,  22,
      20, 26, 12, 1,  24, 7, 40, 3,  41, 1,  0,  0, 40, 3,  38, 5,  0,  0,
      40, 24, 38, 25, 0,  0, 40, 24, 41, 22, 0,  0, 40, 12, 38, 17, 0,  0,
      43, 24, 38, 25, 0,  0, 43, 24, 38, 41, 22, 0, 43, 12, 38, 17, 0,  0,
      41, 1,  38, 5,  0,  0, 41, 3,  38, 6,  0,  0, 41, 24, 38, 26, 0,  0,
      38, 5,  41, 1,  0,  0, 38, 5,  40, 3,  0,  0, 38, 6,  41, 3,  0,  0,
      3,  5,  1,  6,  0,  0, 5,  7,  1,  23, 0,  0, 5,  22, 3,  23, 0,  0,
      5,  22, 1,  25, 0,  0, 5,  24, 3,  25, 0,  0, 5,  24, 1,  26, 0,  0,
      33, 5,  32, 1,  3,  0, 29, 5,  32, 1,  0,  0, 29, 5,  30, 3,  0,  0,
      21, 5,  28, 1,  0,  0, 6,  7,  1,  25, 0,  0, 6,  7,  3,  23, 0,  0,
      6,  22, 3,  25, 0,  0, 24, 6,  3,  26, 0,  0, 29, 6,  32, 3,  0,  0,
      21, 6,  28, 3,  0,  0, 3,  23, 1,  25, 0,  0, 22, 23, 25, 7,  0,  0,
      24, 23, 25, 22, 0,  0, 24, 23, 26, 7,  0,  0, 33, 23, 29, 26, 0,  0,
      29, 23, 32, 7,  0,  0, 29, 23, 30, 22, 0,  0, 21, 23, 28, 7,  0,  0,
      3,  25, 1,  26, 0,  0, 25, 22, 26, 7,  0,  0, 33, 25, 31, 26, 0,  0,
      29, 25, 32, 22, 0,  0, 29, 25, 30, 24, 0,  0, 21, 25, 28, 22, 0,  0,
      3,  26, 24, 6,  0,  0, 29, 26, 32, 24, 0,  0, 17, 24, 12, 25, 0,  0,
      12, 17, 19, 0,  0,  0, 19, 11, 18, 12, 0,  0, 19, 11, 18, 12, 0,  0,
      19, 12, 17, 0,  0,  0, 28, 3,  30, 1,  0,  0, 28, 24, 21, 26, 0,  0,
      30, 33, 31, 32, 0,  0, 30, 3,  32, 1,  0,  0, 8,  9,  7,  10, 0,  0,
      12, 8,  20, 0,  0,  0, 7,  10, 8,  9,  0,  0, 1,  10, 2,  9,  0,  0,
      20, 7,  12, 10, 0,  0, 7,  9,  0,  0,  0,  0, 7,  9,  0,  0,  0,  0,
      7,  9,  0,  0,  0,  0, 7,  9,  0,  0,  0,  0, 1,  7,  22, 0,  0,  0,
      1,  7,  22, 0,  0,  0, 1,  7,  22, 0,  0,  0, 1,  7,  22, 0,  0,  0,
      1,  9,  7,  22, 0,  0, 9,  7,  0,  0,  0,  0, 1,  3,  0,  0,  0,  0,
      1,  3,  0,  0,  0,  0, 1,  3,  0,  0,  0,  0, 1,  3,  0,  0,  0,  0,
      3,  1,  0,  0,  0,  0, 22, 1,  7,  0,  0,  0, 1,  22, 3,  7,  0,  0,
      1,  22, 24, 0,  0,  0, 1,  22, 24, 0,  0,  0, 1,  22, 24, 0,  0,  0,
      1,  22, 24, 0,  0,  0, 24, 1,  22, 0,  0,  0, 1,  24, 3,  22, 0,  0,
      9,  7,  0,  0,  0,  0, 3,  9,  22, 0,  0,  0, 3,  22, 1,  24, 0,  0,
      22, 1,  7,  0,  0,  0, 24, 1,  22, 0,  0,  0, 7,  22, 1,  9,  0,  0,
      22, 24, 7,  0,  0,  0, 24, 7,  22, 0,  0,  0, 29, 33, 31, 0,  0,  0,
      33, 1,  31, 3,  0,  0, 31, 1,  29, 3,  0,  0, 31, 1,  33, 0,  0,  0,
      33, 7,  31, 22, 0,  0, 31, 7,  34, 1,  0,  0, 31, 7,  11, 1,  3,  0,
      29, 7,  11, 3,  0,  0, 29, 7,  11, 1,  0,  0, 29, 9,  12, 3,  0,  0,
      29, 9,  11, 24, 0,  0, 29, 9,  34, 7,  0,  0, 21, 9,  11, 7,  0,  0,
      33, 22, 31, 24, 0,  0, 33, 27, 34, 31, 0,  0, 35, 33, 31, 36, 0,  0,
      31, 22, 29, 24, 0,  0, 31, 22, 37, 1,  0,  0, 31, 22, 35, 1,  0,  0,
      31, 22, 36, 0,  0,  0, 31, 27, 33, 11, 0,  0, 31, 35, 34, 33, 0,  0,
      29, 12, 34, 11, 0,  0, 29, 22, 34, 1,  0,  0, 29, 34, 31, 27, 0,  0,
      29, 27, 31, 11, 0,  0, 29, 35, 34, 31, 0,  0, 3,  27, 34, 1,  0,  0,
      12, 1,  11, 22, 0,  0, 11, 1,  27, 0,  0,  0, 1,  27, 11, 3,  0,  0,
      35, 1,  34, 3,  0,  0, 34, 7,  27, 22, 0,  0, 27, 7,  11, 22, 0,  0,
      27, 7,  12, 1,  0,  0, 35, 7,  31, 9,  0,  0, 35, 7,  34, 22, 0,  0,
      35, 11, 31, 12, 0,  0, 24, 27, 34, 22, 0,  0, 34, 22, 24, 27, 0,  0,
      36, 22, 35, 24, 0,  0, 27, 22, 11, 24, 0,  0, 35, 22, 34, 24, 0,  0,
      27, 34, 11, 0,  0,  0, 35, 27, 36, 11, 0,  0, 35, 34, 36, 0,  0,  0,
      37, 33, 31, 36, 0,  0, 37, 31, 34, 33, 0,  0, 29, 36, 31, 35, 0,  0,
      29, 36, 37, 31, 0,  0, 29, 37, 34, 31, 0,  0, 37, 3,  36, 1,  0,  0,
      36, 1,  35, 3,  0,  0, 37, 1,  34, 3,  0,  0, 37, 1,  35, 1,  0,  0,
      35, 1,  36, 0,  0,  0, 36, 7,  37, 22, 0,  0, 36, 7,  35, 22, 0,  0,
      37, 7,  34, 22, 0,  0, 36, 22, 37, 24, 0,  0, 37, 22, 34, 24, 0,  0,
      34, 37, 36, 27, 0,  0, 37, 27, 34, 0,  0,  0, 37, 27, 36, 11, 0,  0,
      35, 36, 37, 36, 0,  0, 37, 35, 34, 36, 0,  0, 37, 34, 36, 0,  0,  0,
      29, 7,  1,  27, 0,  0, 31, 7,  35, 0,  0,  0, 12, 11, 7,  0,  0,  0,
      21, 12, 11, 0,  0,  0, 11, 7,  12, 0,  0,  0, 11, 9,  12, 7,  0,  0,
      21, 7,  11, 0,  0,  0, 29, 9,  12, 1,  0,  0, 29, 9,  11, 1,  22, 0,
      36, 22, 34, 1,  24, 0, 31, 9,  34, 22, 0,  0, 31, 9,  24, 27, 0,  0,
      34, 7,  11, 1,  22, 0};
  const int nuv[NUM_GAS_REACTIONS * 6] = {
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, 1,  1, 0, 0, 0, -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, 2,  0, 0, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, 1,  1, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 0, 0, 0, -1, -1, 1, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, 1,  1, 0, 0, 0,
      -1, 1,  1, 0, 0, 0, -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, 1,  1, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 2, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 2, 0, 0, 0, -1, -1, 2, 0, 0, 0, -1, -1, 3, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 2, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 2, 1, 0, 0,
      -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 2, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 2, 1, 1, 0, -1, -1, 1, 1, 1, 0, -2, 2,  1, 1, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, -1, 2, 0, 0, 0, -1, 2,  0, 0, 0, 0,
      -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 2, 0, 0, -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 2, 0, 0, -1, -1, 2, 2, 0, 0, -1, 1,  0, 0, 0, 0,
      -1, 2,  0, 0, 0, 0, -1, -1, 2, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 2, 0, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, 1,  0, 0, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 2, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 2, 1, 0, 0, -1, -1, 2, 1, 0, 0, -1, -1, 2, 2, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 2, 1, 0, 0,
      -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 1, 2, 0, -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 2, 0, 0, -1, -1, 2, 1, 0, 0, -1, -1, 1, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 2, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 0, -1, -1, 2, 1, 0, 0,
      -1, -1, 1, 2, 0, 0, -1, -1, 2, 2, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 2, 1, 0, 0, -1, -1, 1, 3, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0, -1, -1, 1, 2, 0, 0,
      -1, -1, 1, 2, 1, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 1, 2, 0, -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 2, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 3, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 2, 1, 0, 0, -1, -1, 2, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 2, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0, -1, -1, 2, 1, 0, 0,
      -1, -1, 1, 2, 1, 0, -1, -1, 3, 1, 0, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 1, 2, 0, -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 1, 1,
      -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 2, 1, 0, 0, -1, -1, 1, 1, 1, 0, -1, -1, 2, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -2, 1,  0, 0, 0, 0, -2, 1,  0, 0, 0, 0,
      -2, 1,  0, 0, 0, 0, -2, 1,  0, 0, 0, 0, -1, -1, 1, 0, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, -1, 1, 0, 0, 0, -1, -1, 1, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, 2,  0, 0, 0, 0, -2, 1,  0, 0, 0, 0,
      -2, 1,  0, 0, 0, 0, -2, 1,  0, 0, 0, 0, -2, 1,  0, 0, 0, 0,
      -1, 2,  0, 0, 0, 0, -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, -1, 1, 0, 0, 0, -1, -1, 1, 0, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, 2,  0, 0, 0, 0, -1, -1, 2, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, 1,  1, 0, 0, 0, -1, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -2, 1,  1, 0, 0, 0, -1, -1, 2, 0, 0, 0, -1, -1, 2, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 2, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -2, 1,  1, 0, 0, 0, -1, -1, 1, 1, 0, 0, -2, 1,  1, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 2, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0, -2, 1,  1, 0, 0, 0,
      -1, -1, 1, 1, 0, 0, -1, -1, 1, 0, 0, 0, -1, 1,  1, 0, 0, 0,
      -1, -1, 2, 0, 0, 0, -1, -1, 1, 0, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 0, 0, 0, -1, -1, 1, 2, 0, 0, -1, -1, 1, 1, 1, 0,
      -1, -1, 1, 1, 1, 0, -1, -1, 1, 1, 0, 0, -1, -1, 1, 1, 0, 0,
      -1, -1, 1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 6;
  } else {
    if (i > NUM_GAS_REACTIONS) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 6 + j] + 1;
        nu[j] = nuv[(i - 1) * 6 + j];
      }
    }
  }
}

// save atomic weights into array
void atomicWeight(amrex::Real *awt) {
  awt[0] = 0.000549;  // E
  awt[1] = 15.999000; // O
  awt[2] = 12.011000; // C
  awt[3] = 39.950000; // Ar
  awt[4] = 1.008000;  // H
}

// get atomic weight for all elements
void CKAWT(amrex::Real *awt) { atomicWeight(awt); }

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void CKNCF(int *ncf) {
  int kd = 5;
  // Zero ncf
  for (int id = 0; id < kd * 44; ++id) {
    ncf[id] = 0;
  }

  // E
  ncf[0 * kd + 0] = 1; // E

  // H
  ncf[1 * kd + 4] = 1; // H

  // H-
  ncf[2 * kd + 0] = 1; // E
  ncf[2 * kd + 4] = 1; // H

  // H2
  ncf[3 * kd + 4] = 2; // H

  // H2v
  ncf[4 * kd + 4] = 2; // H

  // H2+
  ncf[5 * kd + 0] = -1; // E
  ncf[5 * kd + 4] = 2;  // H

  // H3+
  ncf[6 * kd + 0] = -1; // E
  ncf[6 * kd + 4] = 3;  // H

  // O
  ncf[7 * kd + 1] = 1; // O

  // O-
  ncf[8 * kd + 0] = 1; // E
  ncf[8 * kd + 1] = 1; // O

  // O2
  ncf[9 * kd + 1] = 2; // O

  // O2-
  ncf[10 * kd + 0] = 1; // E
  ncf[10 * kd + 1] = 2; // O

  // CO
  ncf[11 * kd + 2] = 1; // C
  ncf[11 * kd + 1] = 1; // O

  // CO2
  ncf[12 * kd + 2] = 1; // C
  ncf[12 * kd + 1] = 2; // O

  // CO2v1
  ncf[13 * kd + 2] = 1; // C
  ncf[13 * kd + 1] = 2; // O

  // CO2v2
  ncf[14 * kd + 2] = 1; // C
  ncf[14 * kd + 1] = 2; // O

  // CO2v3
  ncf[15 * kd + 2] = 1; // C
  ncf[15 * kd + 1] = 2; // O

  // CO2v4
  ncf[16 * kd + 2] = 1; // C
  ncf[16 * kd + 1] = 2; // O

  // CO2+
  ncf[17 * kd + 2] = 1;  // C
  ncf[17 * kd + 0] = -1; // E
  ncf[17 * kd + 1] = 2;  // O

  // C2O3+
  ncf[18 * kd + 2] = 2;  // C
  ncf[18 * kd + 0] = -1; // E
  ncf[18 * kd + 1] = 3;  // O

  // C2O4+
  ncf[19 * kd + 2] = 2;  // C
  ncf[19 * kd + 0] = -1; // E
  ncf[19 * kd + 1] = 4;  // O

  // CO3-
  ncf[20 * kd + 2] = 1; // C
  ncf[20 * kd + 0] = 1; // E
  ncf[20 * kd + 1] = 3; // O

  // C
  ncf[21 * kd + 2] = 1; // C

  // OH
  ncf[22 * kd + 4] = 1; // H
  ncf[22 * kd + 1] = 1; // O

  // OH+
  ncf[23 * kd + 0] = -1; // E
  ncf[23 * kd + 4] = 1;  // H
  ncf[23 * kd + 1] = 1;  // O

  // H2O
  ncf[24 * kd + 4] = 2; // H
  ncf[24 * kd + 1] = 1; // O

  // H2O+
  ncf[25 * kd + 0] = -1; // E
  ncf[25 * kd + 4] = 2;  // H
  ncf[25 * kd + 1] = 1;  // O

  // H3O+
  ncf[26 * kd + 0] = -1; // E
  ncf[26 * kd + 4] = 3;  // H
  ncf[26 * kd + 1] = 1;  // O

  // HCO
  ncf[27 * kd + 2] = 1; // C
  ncf[27 * kd + 4] = 1; // H
  ncf[27 * kd + 1] = 1; // O

  // CH+
  ncf[28 * kd + 2] = 1;  // C
  ncf[28 * kd + 0] = -1; // E
  ncf[28 * kd + 4] = 1;  // H

  // CH2
  ncf[29 * kd + 2] = 1; // C
  ncf[29 * kd + 4] = 2; // H

  // CH2+
  ncf[30 * kd + 2] = 1;  // C
  ncf[30 * kd + 0] = -1; // E
  ncf[30 * kd + 4] = 2;  // H

  // CH3
  ncf[31 * kd + 2] = 1; // C
  ncf[31 * kd + 4] = 3; // H

  // CH3+
  ncf[32 * kd + 2] = 1;  // C
  ncf[32 * kd + 0] = -1; // E
  ncf[32 * kd + 4] = 3;  // H

  // CH4
  ncf[33 * kd + 2] = 1; // C
  ncf[33 * kd + 4] = 4; // H

  // CH2O
  ncf[34 * kd + 2] = 1; // C
  ncf[34 * kd + 4] = 2; // H
  ncf[34 * kd + 1] = 1; // O

  // CH3O
  ncf[35 * kd + 2] = 1; // C
  ncf[35 * kd + 4] = 3; // H
  ncf[35 * kd + 1] = 1; // O

  // CH3OH
  ncf[36 * kd + 2] = 1; // C
  ncf[36 * kd + 4] = 4; // H
  ncf[36 * kd + 1] = 1; // O

  // CH2OH
  ncf[37 * kd + 2] = 1; // C
  ncf[37 * kd + 4] = 3; // H
  ncf[37 * kd + 1] = 1; // O

  // AR
  ncf[38 * kd + 3] = 1; // Ar

  // ARe
  ncf[39 * kd + 3] = 1; // Ar

  // AR+
  ncf[40 * kd + 3] = 1;  // Ar
  ncf[40 * kd + 0] = -1; // E

  // ARH+
  ncf[41 * kd + 3] = 1;  // Ar
  ncf[41 * kd + 0] = -1; // E
  ncf[41 * kd + 4] = 1;  // H

  // AR2e
  ncf[42 * kd + 3] = 2; // Ar

  // AR2+
  ncf[43 * kd + 3] = 2;  // Ar
  ncf[43 * kd + 0] = -1; // E
}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string> &ename) {
  ename.resize(5);
  ename[0] = "E";
  ename[1] = "O";
  ename[2] = "C";
  ename[3] = "Ar";
  ename[4] = "H";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string> &kname) {
  kname.resize(44);
  kname[0] = "E";
  kname[1] = "H";
  kname[2] = "H-";
  kname[3] = "H2";
  kname[4] = "H2v";
  kname[5] = "H2+";
  kname[6] = "H3+";
  kname[7] = "O";
  kname[8] = "O-";
  kname[9] = "O2";
  kname[10] = "O2-";
  kname[11] = "CO";
  kname[12] = "CO2";
  kname[13] = "CO2v1";
  kname[14] = "CO2v2";
  kname[15] = "CO2v3";
  kname[16] = "CO2v4";
  kname[17] = "CO2+";
  kname[18] = "C2O3+";
  kname[19] = "C2O4+";
  kname[20] = "CO3-";
  kname[21] = "C";
  kname[22] = "OH";
  kname[23] = "OH+";
  kname[24] = "H2O";
  kname[25] = "H2O+";
  kname[26] = "H3O+";
  kname[27] = "HCO";
  kname[28] = "CH+";
  kname[29] = "CH2";
  kname[30] = "CH2+";
  kname[31] = "CH3";
  kname[32] = "CH3+";
  kname[33] = "CH4";
  kname[34] = "CH2O";
  kname[35] = "CH3O";
  kname[36] = "CH3OH";
  kname[37] = "CH2OH";
  kname[38] = "AR";
  kname[39] = "ARe";
  kname[40] = "AR+";
  kname[41] = "ARH+";
  kname[42] = "AR2e";
  kname[43] = "AR2+";
}