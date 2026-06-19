/* jpg.h — minimal, dependency-free baseline JPEG reader/writer.
 *
 * Scope: write and read baseline (sequential DCT, Huffman) JFIF JPEG — the
 * universal interchange profile. No progressive, arithmetic, or 12-bit support.
 * Self-contained: own DCT, own standard Huffman/quant tables, no libjpeg.
 * Mirrors the style of the sibling `tiff` / `png` modules.
 *
 * Pixel model: row-major, chunky, 8-bit samples. Supported:
 *     channels 1  (grayscale)          -> 1-component JPEG
 *     channels 3  (RGB)                -> YCbCr 4:4:4 (writer); reader also
 *                                         handles 4:2:2 / 4:2:0 / 4:4:4
 *
 * The writer emits SOI/APP0/DQT/SOF0/DHT/SOS/EOI with the standard Annex-K
 * Huffman tables and quality-scaled quant tables, no chroma subsampling (4:4:4).
 * The reader parses baseline files from any encoder: 1 or 3 components, the
 * common sampling factors, restart intervals, standard byte-stuffing.
 *
 * JPEG is lossy: a written-then-read image will differ from the original within
 * the quantization error set by `quality` (1..100; ~90 is visually clean).
 *
 * All I/O is via stdio. jpg_read allocates img->data; release with jpg_free.
 */
#ifndef JPG_H
#define JPG_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Return codes: 0 on success, negative on failure. */
enum {
    JPG_OK           =  0,
    JPG_EOPEN        = -1,   /* could not open/create the file           */
    JPG_EIO          = -2,   /* short read/write or stdio error          */
    JPG_EFORMAT      = -3,   /* not a JPEG / malformed                   */
    JPG_EUNSUPPORTED = -4,   /* valid JPEG but a profile we don't handle */
    JPG_EMEM         = -5,   /* allocation failed                        */
    JPG_EINVAL       = -6,   /* bad argument to the API                  */
};

typedef struct {
    uint32_t width;      /* pixels per row                                */
    uint32_t height;     /* number of rows                                */
    uint16_t channels;   /* 1 = grayscale, 3 = RGB                        */
    void    *data;       /* width*height*channels uint8_t, row-major chunky */
} jpg_image;

/* Human-readable name for a return code. */
const char *jpg_strerror(int code);

/* Write `img` as baseline JPEG at `quality` (1..100). img->data must hold
 * width*height*channels bytes. Returns JPG_OK or < 0. */
int jpg_write(const char *path, const jpg_image *img, int quality);

/* Read a baseline JPEG. On success fills *img and allocates img->data
 * (free with jpg_free). On failure *img is zeroed and returns < 0. */
int jpg_read(const char *path, jpg_image *img);

/* Free img->data and zero the struct. Safe on a zeroed/already-freed image. */
void jpg_free(jpg_image *img);

#ifdef __cplusplus
}
#endif

#endif /* JPG_H */
