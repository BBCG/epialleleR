// Global definitions


// FNV-1a hash
// http://www.isthe.com/chongo/tech/comp/fnv/
#define FNV1a_OFFSET_BASIS 14695981039346656037u
#define FNV1a_PRIME 1099511628211u
#define fnv_add(hash, pointer, size) {     /* hash starts with offset_basis */ \
  for (unsigned int idx=0; idx<size; idx++) {        /* cycle through bytes */ \
    hash ^= *(pointer+idx);                       /* hash xor octet_of_data */ \
    hash *= FNV1a_PRIME;                                /* hash * FNV_prime */ \
  }                                                                            \
}                                                                              

// Cytosine context to index (ctx_to_idx) conversion:
// ctx  bin       +2        >>2&15  idx
// +    00101011  00101101  1011    11
// -    00101101  00101111  1011    11
// .    00101110  00110000  1100    12
// H    01001000  01001010  0010    2
// U    01010101  01010111  0101    5
// X    01011000  01011010  0110    6
// Z    01011010  01011100  0111    7
// h    01101000  01101010  1010    10
// u    01110101  01110111  1101    13
// x    01111000  01111010  1110    14
// z    01111010  01111100  1111    15
#define ctx_to_idx(c) ((c+2)>>2) & 15




