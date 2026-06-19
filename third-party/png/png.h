/* png.h — minimal, dependency-free PNG reader/writer.
 *
 * Scope: read and write baseline non-interlaced PNG (the format every viewer
 * and library emits by default), with a self-contained DEFLATE/INFLATE so there
 * is no zlib dependency. Mirrors the style of the sibling `tiff` module.
 *
 * Pixel model: row-major, chunky (interleaved) samples. Supported:
 *
 *     bitdepth 8  or 16   (16-bit samples held NATIVE-endian in memory)
 *     channels 1 2 3 4    (gray, gray+alpha, RGB, RGBA -> PNG color types 0/4/2/6)
 *
 * The writer emits one IDAT (zlib: fixed-Huffman DEFLATE over LZ77), adaptive
 * per-row filtering (None/Sub/Up/Average/Paeth, min-sum heuristic), CRC32 per
 * chunk. The reader parses signature + IHDR/IDAT/IEND, runs full INFLATE
 * (stored/fixed/dynamic Huffman), and reverses the filters. It rejects what it
 * cannot safely return: palette (color type 3), interlaced (Adam7), or bit
 * depths other than 8/16.
 *
 * 16-bit note: PNG stores samples big-endian on disk; this API exchanges them
 * as native-endian uint16_t in img->data (write converts out, read converts in).
 *
 * All I/O is via stdio. png_read allocates img->data; release with png_free.
 */
#ifndef PNG_H
#define PNG_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Return codes: 0 on success, negative on failure. */
enum {
    PNG_OK           =  0,
    PNG_EOPEN        = -1,   /* could not open/create the file           */
    PNG_EIO          = -2,   /* short read/write or stdio error          */
    PNG_EFORMAT      = -3,   /* not a PNG / malformed / bad CRC          */
    PNG_EUNSUPPORTED = -4,   /* valid PNG but a layout we don't handle   */
    PNG_EMEM         = -5,   /* allocation failed                        */
    PNG_EINVAL       = -6,   /* bad argument to the API                  */
};

typedef struct {
    uint32_t width;      /* pixels per row                                */
    uint32_t height;     /* number of rows                                */
    uint16_t channels;   /* 1=gray 2=gray+alpha 3=RGB 4=RGBA              */
    uint16_t bitdepth;   /* 8 or 16                                       */
    void    *data;       /* width*height*channels samples, row-major chunky;
                          * 8-bit => uint8_t, 16-bit => native uint16_t   */
} png_image;

/* Human-readable name for a return code. */
const char *png_strerror(int code);

/* Write `img` as a baseline PNG. img->data must hold
 * width*height*channels*(bitdepth/8) bytes. Returns PNG_OK or < 0. */
int png_write(const char *path, const png_image *img);

/* Read a baseline PNG. On success fills *img and allocates img->data
 * (free with png_free). On failure *img is zeroed and returns < 0. */
int png_read(const char *path, png_image *img);

/* Free img->data and zero the struct. Safe on a zeroed/already-freed image. */
void png_free(png_image *img);

#ifdef __cplusplus
}
#endif

#endif /* PNG_H */
