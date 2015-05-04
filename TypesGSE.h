#ifndef TYPESGSE_H
#define TYPESGSE_H

//#define USE_NEW_K

#ifdef USE_NEW_K
#define AK (97.0/120.0)
#define BK (12.0/120.0)
#define CK (-1.0/240.0)
#endif

#ifdef USE_NEW_KK
#define AK (7033.0/3780.0)
#define BK (-1741.0/2520.0)
#define CK (-223.0/1260.0)
#define DK (11.0/1512.0)
#endif

template <typename T>
struct xyzq_t {
  T x, y, z, q;
};

const double pi_dbl = 3.14159265358979323846;

#endif // TYPESGSE_H
