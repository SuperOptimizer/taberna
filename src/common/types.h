/* types.h — shared fixed-width scalar aliases used across taberna's first-party
 * sources. Matches the convention in src/segmentation/snic.h (which predates this
 * header and defines the same set inline; C11 permits the identical typedefs to
 * coexist if both are included in one translation unit). */
#ifndef TABERNA_TYPES_H
#define TABERNA_TYPES_H

#include <stdint.h>

typedef uint8_t  u8;
typedef int8_t   s8;
typedef uint16_t u16;
typedef int16_t  s16;
typedef uint32_t u32;
typedef int32_t  s32;
typedef uint64_t u64;
typedef int64_t  s64;
typedef float    f32;
typedef double   f64;

#endif // TABERNA_TYPES_H
