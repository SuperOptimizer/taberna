/* jpg.c — see jpg.h. Baseline JFIF JPEG, self-contained (own DCT/Huffman). */
#include "jpg.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const char *jpg_strerror(int code) {
    switch (code) {
        case JPG_OK:           return "ok";
        case JPG_EOPEN:        return "could not open file";
        case JPG_EIO:          return "I/O error";
        case JPG_EFORMAT:      return "not a JPEG / malformed";
        case JPG_EUNSUPPORTED: return "unsupported JPEG profile";
        case JPG_EMEM:         return "out of memory";
        case JPG_EINVAL:       return "invalid argument";
    }
    return "unknown error";
}

/* ===================================================== shared tables ===== */

static const int ZIGZAG[64] = {
     0, 1, 8,16, 9, 2, 3,10, 17,24,32,25,18,11, 4, 5,
    12,19,26,33,40,48,41,34, 27,20,13, 6, 7,14,21,28,
    35,42,49,56,57,50,43,36, 29,22,15,23,30,37,44,51,
    58,59,52,45,38,31,39,46, 53,60,61,54,47,55,62,63 };

static const uint8_t STD_LUM_Q[64] = {
    16,11,10,16,24,40,51,61, 12,12,14,19,26,58,60,55,
    14,13,16,24,40,57,69,56, 14,17,22,29,51,87,80,62,
    18,22,37,56,68,109,103,77, 24,35,55,64,81,104,113,92,
    49,64,78,87,103,121,120,101, 72,92,95,98,112,100,103,99 };
static const uint8_t STD_CHR_Q[64] = {
    17,18,24,47,99,99,99,99, 18,21,26,66,99,99,99,99,
    24,26,56,99,99,99,99,99, 47,66,99,99,99,99,99,99,
    99,99,99,99,99,99,99,99, 99,99,99,99,99,99,99,99,
    99,99,99,99,99,99,99,99, 99,99,99,99,99,99,99,99 };

/* standard Annex-K Huffman specs: BITS[1..16] then HUFFVAL */
static const uint8_t DC_LUM_BITS[17] = {0,0,1,5,1,1,1,1,1,1,0,0,0,0,0,0,0};
static const uint8_t DC_LUM_VAL[12]  = {0,1,2,3,4,5,6,7,8,9,10,11};
static const uint8_t DC_CHR_BITS[17] = {0,0,3,1,1,1,1,1,1,1,1,1,0,0,0,0,0};
static const uint8_t DC_CHR_VAL[12]  = {0,1,2,3,4,5,6,7,8,9,10,11};
static const uint8_t AC_LUM_BITS[17] = {0,0,2,1,3,3,2,4,3,5,5,4,4,0,0,1,125};
static const uint8_t AC_LUM_VAL[162] = {
    0x01,0x02,0x03,0x00,0x04,0x11,0x05,0x12,0x21,0x31,0x41,0x06,0x13,0x51,0x61,0x07,
    0x22,0x71,0x14,0x32,0x81,0x91,0xa1,0x08,0x23,0x42,0xb1,0xc1,0x15,0x52,0xd1,0xf0,
    0x24,0x33,0x62,0x72,0x82,0x09,0x0a,0x16,0x17,0x18,0x19,0x1a,0x25,0x26,0x27,0x28,
    0x29,0x2a,0x34,0x35,0x36,0x37,0x38,0x39,0x3a,0x43,0x44,0x45,0x46,0x47,0x48,0x49,
    0x4a,0x53,0x54,0x55,0x56,0x57,0x58,0x59,0x5a,0x63,0x64,0x65,0x66,0x67,0x68,0x69,
    0x6a,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7a,0x83,0x84,0x85,0x86,0x87,0x88,0x89,
    0x8a,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9a,0xa2,0xa3,0xa4,0xa5,0xa6,0xa7,
    0xa8,0xa9,0xaa,0xb2,0xb3,0xb4,0xb5,0xb6,0xb7,0xb8,0xb9,0xba,0xc2,0xc3,0xc4,0xc5,
    0xc6,0xc7,0xc8,0xc9,0xca,0xd2,0xd3,0xd4,0xd5,0xd6,0xd7,0xd8,0xd9,0xda,0xe1,0xe2,
    0xe3,0xe4,0xe5,0xe6,0xe7,0xe8,0xe9,0xea,0xf1,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,
    0xf9,0xfa };
static const uint8_t AC_CHR_BITS[17] = {0,0,2,1,2,4,4,3,4,7,5,4,4,0,1,2,119};
static const uint8_t AC_CHR_VAL[162] = {
    0x00,0x01,0x02,0x03,0x11,0x04,0x05,0x21,0x31,0x06,0x12,0x41,0x51,0x07,0x61,0x71,
    0x13,0x22,0x32,0x81,0x08,0x14,0x42,0x91,0xa1,0xb1,0xc1,0x09,0x23,0x33,0x52,0xf0,
    0x15,0x62,0x72,0xd1,0x0a,0x16,0x24,0x34,0xe1,0x25,0xf1,0x17,0x18,0x19,0x1a,0x26,
    0x27,0x28,0x29,0x2a,0x35,0x36,0x37,0x38,0x39,0x3a,0x43,0x44,0x45,0x46,0x47,0x48,
    0x49,0x4a,0x53,0x54,0x55,0x56,0x57,0x58,0x59,0x5a,0x63,0x64,0x65,0x66,0x67,0x68,
    0x69,0x6a,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7a,0x82,0x83,0x84,0x85,0x86,0x87,
    0x88,0x89,0x8a,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9a,0xa2,0xa3,0xa4,0xa5,
    0xa6,0xa7,0xa8,0xa9,0xaa,0xb2,0xb3,0xb4,0xb5,0xb6,0xb7,0xb8,0xb9,0xba,0xc2,0xc3,
    0xc4,0xc5,0xc6,0xc7,0xc8,0xc9,0xca,0xd2,0xd3,0xd4,0xd5,0xd6,0xd7,0xd8,0xd9,0xda,
    0xe2,0xe3,0xe4,0xe5,0xe6,0xe7,0xe8,0xe9,0xea,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,
    0xf9,0xfa };

/* ============================================================ DCT ======= */

static double COST[8][8];   /* COST[x][u] = cos((2x+1)u pi/16) */
static int cost_ready = 0;
static void cost_init(void) {
    for (int x = 0; x < 8; x++)
        for (int u = 0; u < 8; u++)
            COST[x][u] = cos((2 * x + 1) * u * M_PI / 16.0);
    cost_ready = 1;
}
static void fdct8(const double in[64], double out[64]) {
    double tmp[64];
    for (int y = 0; y < 8; y++)               /* rows */
        for (int u = 0; u < 8; u++) {
            double s = 0; for (int x = 0; x < 8; x++) s += in[y * 8 + x] * COST[x][u];
            tmp[y * 8 + u] = s * (u == 0 ? 0.70710678118654752 : 1.0) * 0.5;
        }
    for (int u = 0; u < 8; u++)               /* cols */
        for (int v = 0; v < 8; v++) {
            double s = 0; for (int y = 0; y < 8; y++) s += tmp[y * 8 + u] * COST[y][v];
            out[v * 8 + u] = s * (v == 0 ? 0.70710678118654752 : 1.0) * 0.5;
        }
}
static void idct8(const double in[64], double out[64]) {
    /* in[v*8+u] = F(u,v); pair u<->x and v<->y, consistent with fdct8 */
    double tmp[64];
    for (int v = 0; v < 8; v++)               /* sum over u -> spatial x */
        for (int x = 0; x < 8; x++) {
            double s = 0; for (int u = 0; u < 8; u++)
                s += (u == 0 ? 0.70710678118654752 : 1.0) * in[v * 8 + u] * COST[x][u];
            tmp[v * 8 + x] = s * 0.5;
        }
    for (int x = 0; x < 8; x++)               /* sum over v -> spatial y */
        for (int y = 0; y < 8; y++) {
            double s = 0; for (int v = 0; v < 8; v++)
                s += (v == 0 ? 0.70710678118654752 : 1.0) * tmp[v * 8 + x] * COST[y][v];
            out[y * 8 + x] = s * 0.5;
        }
}

static int iclamp(int v) { return v < 0 ? 0 : v > 255 ? 255 : v; }

/* ===================================================== quant scaling ==== */

static void scale_qtable(const uint8_t base[64], int quality, uint8_t out[64]) {
    if (quality < 1) quality = 1;
    if (quality > 100) quality = 100;
    int scale = quality < 50 ? 5000 / quality : 200 - 2 * quality;
    for (int i = 0; i < 64; i++) {
        int q = (base[i] * scale + 50) / 100;
        out[i] = (uint8_t)(q < 1 ? 1 : q > 255 ? 255 : q);
    }
}

/* =================================================== Huffman (encode) === */

typedef struct { uint16_t code[256]; uint8_t size[256]; } ehuff;

static void build_ehuff(const uint8_t bits[17], const uint8_t *val, ehuff *e) {
    memset(e, 0, sizeof *e);
    uint8_t huffsize[257]; int k = 0;
    for (int l = 1; l <= 16; l++) for (int i = 0; i < bits[l]; i++) huffsize[k++] = (uint8_t)l;
    huffsize[k] = 0; int total = k;
    uint16_t code = 0; int si = huffsize[0]; k = 0;
    uint16_t huffcode[257];
    while (huffsize[k]) {
        while (huffsize[k] == si) { huffcode[k] = code; code++; k++; }
        code <<= 1; si++;
    }
    for (int i = 0; i < total; i++) { e->code[val[i]] = huffcode[i]; e->size[val[i]] = huffsize[i]; }
}

/* =================================================== Huffman (decode) === */

typedef struct {
    int mincode[17], maxcode[17], valptr[17];
    const uint8_t *val;
} dhuff;

static void build_dhuff(const uint8_t bits[17], const uint8_t *val, dhuff *d) {
    int huffsize[257]; int k = 0;
    for (int l = 1; l <= 16; l++) for (int i = 0; i < bits[l]; i++) huffsize[k++] = l;
    huffsize[k] = 0;
    int huffcode[257]; int code = 0, si = huffsize[0]; k = 0;
    while (huffsize[k]) { while (huffsize[k] == si) { huffcode[k++] = code++; } code <<= 1; si++; }
    int p = 0;
    for (int l = 1; l <= 16; l++) {
        if (bits[l]) { d->valptr[l] = p; d->mincode[l] = huffcode[p]; p += bits[l]; d->maxcode[l] = huffcode[p - 1]; }
        else d->maxcode[l] = -1;
    }
    d->val = val;
}

/* ===================================================== bit writer ======= */

typedef struct { FILE *f; uint32_t acc; int nbits; int err; } jbw;

static void jbw_putbits(jbw *w, int code, int size) {     /* MSB-first, 0xFF stuffing */
    if (size == 0) return;
    w->acc = (w->acc << size) | (code & ((1 << size) - 1));
    w->nbits += size;
    while (w->nbits >= 8) {
        uint8_t b = (uint8_t)((w->acc >> (w->nbits - 8)) & 0xff);
        if (fputc(b, w->f) == EOF) w->err = 1;
        if (b == 0xff) { if (fputc(0, w->f) == EOF) w->err = 1; }
        w->nbits -= 8;
    }
}
static void jbw_flush(jbw *w) {                            /* pad with 1-bits */
    if (w->nbits > 0) {
        int pad = 8 - w->nbits;
        jbw_putbits(w, (1 << pad) - 1, pad);
    }
}

/* magnitude category and the bits to emit for a signed coefficient */
static int magcat(int v) { int a = v < 0 ? -v : v, n = 0; while (a) { a >>= 1; n++; } return n; }
static int magbits(int v, int size) { return v >= 0 ? v : v + (1 << size) - 1; }

/* ===================================================== file helpers ===== */

static void w16(FILE *f, int v) { fputc((v >> 8) & 0xff, f); fputc(v & 0xff, f); }

static void write_marker_seg(FILE *f, int marker, const uint8_t *payload, int n) {
    fputc(0xff, f); fputc(marker, f); w16(f, n + 2); if (n) fwrite(payload, 1, n, f);
}

/* ======================================================== encoder ======= */

static void encode_block(jbw *w, const double pix[64], const uint8_t qz[64] /*natural*/,
                         const ehuff *dc, const ehuff *ac, int *pred) {
    double F[64]; fdct8(pix, F);
    int q[64];
    for (int i = 0; i < 64; i++) {
        double v = F[i] / qz[i];
        q[i] = (int)lrint(v);
    }
    /* DC */
    int diff = q[0] - *pred; *pred = q[0];
    int s = magcat(diff);
    jbw_putbits(w, dc->code[s], dc->size[s]);
    if (s) jbw_putbits(w, magbits(diff, s), s);
    /* AC in zigzag order */
    int run = 0;
    for (int k = 1; k < 64; k++) {
        int v = q[ZIGZAG[k]];
        if (v == 0) { run++; continue; }
        while (run > 15) { jbw_putbits(w, ac->code[0xf0], ac->size[0xf0]); run -= 16; }
        int sz = magcat(v);
        int sym = (run << 4) | sz;
        jbw_putbits(w, ac->code[sym], ac->size[sym]);
        jbw_putbits(w, magbits(v, sz), sz);
        run = 0;
    }
    if (run) jbw_putbits(w, ac->code[0x00], ac->size[0x00]);   /* EOB */
}

int jpg_write(const char *path, const jpg_image *img, int quality) {
    if (!path || !img || !img->data) return JPG_EINVAL;
    if (img->channels != 1 && img->channels != 3) return JPG_EINVAL;
    if (!img->width || !img->height) return JPG_EINVAL;
    if (!cost_ready) cost_init();
    int ncomp = img->channels, W = img->width, H = img->height;

    uint8_t lq[64], cq[64];
    scale_qtable(STD_LUM_Q, quality, lq);
    scale_qtable(STD_CHR_Q, quality, cq);

    ehuff edl, eal, edc, eac;
    build_ehuff(DC_LUM_BITS, DC_LUM_VAL, &edl);
    build_ehuff(AC_LUM_BITS, AC_LUM_VAL, &eal);
    build_ehuff(DC_CHR_BITS, DC_CHR_VAL, &edc);
    build_ehuff(AC_CHR_BITS, AC_CHR_VAL, &eac);

    FILE *f = fopen(path, "wb");
    if (!f) return JPG_EOPEN;

    fputc(0xff, f); fputc(0xd8, f);                           /* SOI */
    { uint8_t app0[14] = {'J','F','I','F',0, 1,1, 0, 0,1,0,1, 0,0};
      write_marker_seg(f, 0xe0, app0, 14); }

    /* DQT (zigzag order, precision 8) */
    { uint8_t dqt[65]; dqt[0] = 0x00; for (int i = 0; i < 64; i++) dqt[1 + i] = lq[ZIGZAG[i]];
      write_marker_seg(f, 0xdb, dqt, 65); }
    if (ncomp == 3) { uint8_t dqt[65]; dqt[0] = 0x01; for (int i = 0; i < 64; i++) dqt[1 + i] = cq[ZIGZAG[i]];
      write_marker_seg(f, 0xdb, dqt, 65); }

    /* SOF0 baseline, 4:4:4 */
    { uint8_t sof[6 + 3 * 3]; int p = 0;
      sof[p++] = 8; sof[p++] = (H >> 8) & 0xff; sof[p++] = H & 0xff;
      sof[p++] = (W >> 8) & 0xff; sof[p++] = W & 0xff; sof[p++] = (uint8_t)ncomp;
      for (int c = 0; c < ncomp; c++) { sof[p++] = (uint8_t)(c + 1); sof[p++] = 0x11; sof[p++] = (uint8_t)(c == 0 ? 0 : 1); }
      write_marker_seg(f, 0xc0, sof, p); }

    /* DHT */
    #define PUT_DHT(cls_id, bits, val, nval) do { \
        uint8_t hdr[1 + 16 + 256]; int p = 0; hdr[p++] = (cls_id); \
        for (int i = 1; i <= 16; i++) hdr[p++] = (bits)[i]; \
        for (int i = 0; i < (nval); i++) hdr[p++] = (val)[i]; \
        write_marker_seg(f, 0xc4, hdr, p); } while (0)
    PUT_DHT(0x00, DC_LUM_BITS, DC_LUM_VAL, 12);
    PUT_DHT(0x10, AC_LUM_BITS, AC_LUM_VAL, 162);
    if (ncomp == 3) { PUT_DHT(0x01, DC_CHR_BITS, DC_CHR_VAL, 12);
                      PUT_DHT(0x11, AC_CHR_BITS, AC_CHR_VAL, 162); }
    #undef PUT_DHT

    /* SOS */
    { uint8_t sos[1 + 2 * 3 + 3]; int p = 0; sos[p++] = (uint8_t)ncomp;
      for (int c = 0; c < ncomp; c++) { sos[p++] = (uint8_t)(c + 1); sos[p++] = (uint8_t)(c == 0 ? 0x00 : 0x11); }
      sos[p++] = 0; sos[p++] = 63; sos[p++] = 0;
      write_marker_seg(f, 0xda, sos, p); }

    /* entropy-coded data, MCU = 8x8 per component (4:4:4) */
    jbw w = { f, 0, 0, 0 };
    int predY = 0, predCb = 0, predCr = 0;
    const uint8_t *src = img->data;
    int bx = (W + 7) / 8, by = (H + 7) / 8;
    double bY[64], bCb[64], bCr[64];
    for (int byi = 0; byi < by; byi++) {
        for (int bxi = 0; bxi < bx; bxi++) {
            for (int yy = 0; yy < 8; yy++) {
                int sy = byi * 8 + yy; if (sy >= H) sy = H - 1;
                for (int xx = 0; xx < 8; xx++) {
                    int sx = bxi * 8 + xx; if (sx >= W) sx = W - 1;
                    const uint8_t *px = src + ((size_t)sy * W + sx) * ncomp;
                    if (ncomp == 1) {
                        bY[yy * 8 + xx] = px[0] - 128.0;
                    } else {
                        double r = px[0], g = px[1], b = px[2];
                        bY[yy * 8 + xx]  = (0.299 * r + 0.587 * g + 0.114 * b) - 128.0;
                        bCb[yy * 8 + xx] = (-0.168736 * r - 0.331264 * g + 0.5 * b);
                        bCr[yy * 8 + xx] = (0.5 * r - 0.418688 * g - 0.081312 * b);
                    }
                }
            }
            encode_block(&w, bY, lq, &edl, &eal, &predY);
            if (ncomp == 3) {
                encode_block(&w, bCb, cq, &edc, &eac, &predCb);
                encode_block(&w, bCr, cq, &edc, &eac, &predCr);
            }
        }
    }
    jbw_flush(&w);
    fputc(0xff, f); fputc(0xd9, f);                           /* EOI */
    int err = w.err || ferror(f);
    if (fclose(f) != 0) err = 1;
    return err ? JPG_EIO : JPG_OK;
}

/* ======================================================== decoder ======= */

typedef struct {
    const uint8_t *p; size_t len, pos;
    uint32_t acc; int nbits; int marker; int err;
} jbr;

static int jbr_bit(jbr *r) {
    if (r->nbits == 0) {
        if (r->pos >= r->len) { r->err = 1; return 0; }
        uint8_t b = r->p[r->pos++];
        if (b == 0xff) {
            uint8_t b2 = r->pos < r->len ? r->p[r->pos] : 0xd9;
            if (b2 == 0x00) { r->pos++; }                     /* stuffed FF */
            else { r->marker = b2; r->acc = 0; r->nbits = 0; return 0; }  /* hit marker */
        }
        r->acc = b; r->nbits = 8;
    }
    r->nbits--;
    return (r->acc >> r->nbits) & 1;
}
static int jbr_bits(jbr *r, int n) { int v = 0; while (n--) v = (v << 1) | jbr_bit(r); return v; }
static int jbr_recv_ext(jbr *r, int s) {                      /* signed magnitude decode */
    if (s == 0) return 0;
    int v = jbr_bits(r, s);
    if (v < (1 << (s - 1))) v += (-(1 << s)) + 1;
    return v;
}
static int jbr_decode(jbr *r, const dhuff *h) {
    int code = 0;
    for (int l = 1; l <= 16; l++) {
        code = (code << 1) | jbr_bit(r);
        if (h->maxcode[l] >= 0 && code <= h->maxcode[l])
            return h->val[h->valptr[l] + code - h->mincode[l]];
    }
    r->err = 1; return 0;
}

typedef struct { int id, h, v, tq, td, ta, pred; } jcomp;

int jpg_read(const char *path, jpg_image *img) {
    if (img) memset(img, 0, sizeof *img);
    if (!path || !img) return JPG_EINVAL;
    if (!cost_ready) cost_init();

    FILE *f = fopen(path, "rb");
    if (!f) return JPG_EOPEN;
    fseek(f, 0, SEEK_END); long fsz = ftell(f); fseek(f, 0, SEEK_SET);
    if (fsz < 4) { fclose(f); return JPG_EFORMAT; }
    uint8_t *buf = malloc(fsz);
    if (!buf) { fclose(f); return JPG_EMEM; }
    if (fread(buf, 1, fsz, f) != (size_t)fsz) { free(buf); fclose(f); return JPG_EIO; }
    fclose(f);

    int rc = JPG_EFORMAT;
    uint8_t qt[4][64]; int have_qt[4] = {0};
    dhuff dc_h[4], ac_h[4]; int have_dc[4] = {0}, have_ac[4] = {0};
    /* stash DHT byte specs so dhuff.val pointers stay valid */
    uint8_t hbits[8][17]; uint8_t hval[8][256]; int hb_used = 0;
    jcomp comp[3]; int ncomp = 0, W = 0, H = 0;
    int restart = 0;
    uint8_t *planes[3] = {0}; int pw[3] = {0}, ph[3] = {0};
    int Hmax = 1, Vmax = 1;
    uint8_t *out = NULL;

    if (buf[0] != 0xff || buf[1] != 0xd8) { free(buf); return JPG_EFORMAT; }
    size_t p = 2;
    while (p + 1 < (size_t)fsz) {
        if (buf[p] != 0xff) { p++; continue; }
        int m = buf[p + 1]; p += 2;
        if (m == 0xd9) break;                                 /* EOI */
        if (m == 0x01 || (m >= 0xd0 && m <= 0xd7)) continue;  /* standalone */
        if (p + 2 > (size_t)fsz) { rc = JPG_EFORMAT; goto fail; }
        int seg = (buf[p] << 8) | buf[p + 1];
        const uint8_t *d = buf + p + 2; int dn = seg - 2;
        size_t next = p + seg;
        if (next > (size_t)fsz) { rc = JPG_EFORMAT; goto fail; }

        if (m == 0xdb) {                                      /* DQT */
            int o = 0;
            while (o < dn) {
                int pq = d[o] >> 4, tq = d[o] & 15; o++;
                if (pq != 0 || tq > 3) { rc = JPG_EUNSUPPORTED; goto fail; }
                for (int i = 0; i < 64; i++) qt[tq][ZIGZAG[i]] = d[o + i];
                o += 64; have_qt[tq] = 1;
            }
        } else if (m == 0xc0) {                               /* SOF0 baseline */
            if (d[0] != 8) { rc = JPG_EUNSUPPORTED; goto fail; }
            H = (d[1] << 8) | d[2]; W = (d[3] << 8) | d[4]; ncomp = d[5];
            if ((ncomp != 1 && ncomp != 3) || !W || !H) { rc = JPG_EUNSUPPORTED; goto fail; }
            for (int c = 0; c < ncomp; c++) {
                comp[c].id = d[6 + c * 3];
                comp[c].h = d[7 + c * 3] >> 4; comp[c].v = d[7 + c * 3] & 15;
                comp[c].tq = d[8 + c * 3]; comp[c].pred = 0;
                if (comp[c].h < 1 || comp[c].v < 1 || comp[c].h > 2 || comp[c].v > 2) { rc = JPG_EUNSUPPORTED; goto fail; }
            }
        } else if (m == 0xc1 || m == 0xc2 || m == 0xc3 || (m >= 0xc5 && m <= 0xcf && m != 0xc8)) {
            rc = JPG_EUNSUPPORTED; goto fail;                 /* progressive/arith/etc. */
        } else if (m == 0xc4) {                               /* DHT */
            int o = 0;
            while (o < dn) {
                int tc = d[o] >> 4, th = d[o] & 15; o++;
                if (th > 3 || hb_used >= 8) { rc = JPG_EUNSUPPORTED; goto fail; }
                uint8_t *bits = hbits[hb_used]; uint8_t *val = hval[hb_used];
                bits[0] = 0; int tot = 0;
                for (int i = 1; i <= 16; i++) { bits[i] = d[o + i - 1]; tot += bits[i]; }
                o += 16;
                for (int i = 0; i < tot; i++) val[i] = d[o + i];
                o += tot;
                if (tc == 0) { build_dhuff(bits, val, &dc_h[th]); have_dc[th] = 1; }
                else         { build_dhuff(bits, val, &ac_h[th]); have_ac[th] = 1; }
                hb_used++;
            }
        } else if (m == 0xdd) {                               /* DRI */
            restart = (d[0] << 8) | d[1];
        } else if (m == 0xda) {                               /* SOS -> entropy scan */
            int ns = d[0];
            for (int i = 0; i < ns; i++) {
                int cs = d[1 + i * 2], tdta = d[2 + i * 2];
                for (int c = 0; c < ncomp; c++) if (comp[c].id == cs) { comp[c].td = tdta >> 4; comp[c].ta = tdta & 15; }
            }
            /* set up planes */
            for (int c = 0; c < ncomp; c++) { if (comp[c].h > Hmax) Hmax = comp[c].h; if (comp[c].v > Vmax) Vmax = comp[c].v; }
            int mcux = (W + 8 * Hmax - 1) / (8 * Hmax), mcuy = (H + 8 * Vmax - 1) / (8 * Vmax);
            for (int c = 0; c < ncomp; c++) {
                pw[c] = mcux * comp[c].h * 8; ph[c] = mcuy * comp[c].v * 8;
                planes[c] = calloc((size_t)pw[c] * ph[c], 1);
                if (!planes[c]) { rc = JPG_EMEM; goto fail; }
            }
            for (int c = 0; c < ncomp; c++)
                if (!have_qt[comp[c].tq] || !have_dc[comp[c].td] || !have_ac[comp[c].ta]) { rc = JPG_EFORMAT; goto fail; }

            jbr r = { buf, (size_t)fsz, p + seg, 0, 0, 0, 0 };
            int rst_cnt = restart;
            for (int my = 0; my < mcuy; my++) {
                for (int mx = 0; mx < mcux; mx++) {
                    for (int c = 0; c < ncomp; c++) {
                        for (int vy = 0; vy < comp[c].v; vy++) {
                            for (int vx = 0; vx < comp[c].h; vx++) {
                                double F[64]; for (int i = 0; i < 64; i++) F[i] = 0;
                                int t = jbr_decode(&r, &dc_h[comp[c].td]);
                                int diff = jbr_recv_ext(&r, t);
                                comp[c].pred += diff;
                                const uint8_t *Q = qt[comp[c].tq];
                                F[0] = (double)comp[c].pred * Q[0];
                                for (int k = 1; k < 64; ) {
                                    int rs = jbr_decode(&r, &ac_h[comp[c].ta]);
                                    int run = rs >> 4, sz = rs & 15;
                                    if (sz == 0) { if (run != 15) break; k += 16; continue; }
                                    k += run; if (k >= 64) break;
                                    int val = jbr_recv_ext(&r, sz);
                                    int zz = ZIGZAG[k];
                                    F[zz] = (double)val * Q[zz];
                                    k++;
                                }
                                double pix[64]; idct8(F, pix);
                                int ox = (mx * comp[c].h + vx) * 8, oy = (my * comp[c].v + vy) * 8;
                                for (int yy = 0; yy < 8; yy++)
                                    for (int xx = 0; xx < 8; xx++)
                                        planes[c][(size_t)(oy + yy) * pw[c] + (ox + xx)] = (uint8_t)iclamp((int)lrint(pix[yy * 8 + xx]) + 128);
                                if (r.err) { rc = JPG_EFORMAT; goto fail; }
                            }
                        }
                    }
                    if (restart && --rst_cnt == 0 && !(my == mcuy - 1 && mx == mcux - 1)) {
                        /* align to the RSTn marker and reset predictors */
                        r.nbits = 0; r.acc = 0;
                        while (r.pos + 1 < r.len && !(r.p[r.pos] == 0xff && r.p[r.pos + 1] >= 0xd0 && r.p[r.pos + 1] <= 0xd7)) r.pos++;
                        if (r.pos + 1 < r.len) r.pos += 2;
                        for (int c = 0; c < ncomp; c++) comp[c].pred = 0;
                        r.marker = 0; rst_cnt = restart;
                    }
                }
            }
            rc = JPG_OK;
            break;                                            /* one scan (baseline) */
        }
        p = next;
    }
    if (rc != JPG_OK) goto fail;

    /* assemble output: upsample chroma, color-convert */
    out = malloc((size_t)W * H * ncomp);
    if (!out) { rc = JPG_EMEM; goto fail; }
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            if (ncomp == 1) {
                out[(size_t)y * W + x] = planes[0][(size_t)y * pw[0] + x];
            } else {
                int yv = planes[0][(size_t)(y) * pw[0] + x];
                int cbx = x * comp[1].h / Hmax, cby = y * comp[1].v / Vmax;
                int crx = x * comp[2].h / Hmax, cry = y * comp[2].v / Vmax;
                int cb = planes[1][(size_t)cby * pw[1] + cbx] - 128;
                int cr = planes[2][(size_t)cry * pw[2] + crx] - 128;
                uint8_t *o = out + ((size_t)y * W + x) * 3;
                o[0] = (uint8_t)iclamp((int)lrint(yv + 1.402 * cr));
                o[1] = (uint8_t)iclamp((int)lrint(yv - 0.344136 * cb - 0.714136 * cr));
                o[2] = (uint8_t)iclamp((int)lrint(yv + 1.772 * cb));
            }
        }
    }
    img->width = W; img->height = H; img->channels = (uint16_t)ncomp; img->data = out;
    for (int c = 0; c < ncomp; c++) free(planes[c]);
    free(buf);
    return JPG_OK;

fail:
    for (int c = 0; c < 3; c++) free(planes[c]);
    free(out); free(buf);
    return rc;
}

void jpg_free(jpg_image *img) {
    if (!img) return;
    free(img->data);
    memset(img, 0, sizeof *img);
}
