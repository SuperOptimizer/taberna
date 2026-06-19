/* png.c — see png.h. Baseline PNG with a self-contained DEFLATE/INFLATE. */
#include "png.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *png_strerror(int code) {
    switch (code) {
        case PNG_OK:           return "ok";
        case PNG_EOPEN:        return "could not open file";
        case PNG_EIO:          return "I/O error";
        case PNG_EFORMAT:      return "not a PNG / malformed";
        case PNG_EUNSUPPORTED: return "unsupported PNG layout";
        case PNG_EMEM:         return "out of memory";
        case PNG_EINVAL:       return "invalid argument";
    }
    return "unknown error";
}

/* ============================================================ checksums == */

static uint32_t crc_table[256];
static int crc_ready = 0;
static void crc_init(void) {
    for (uint32_t n = 0; n < 256; n++) {
        uint32_t c = n;
        for (int k = 0; k < 8; k++) c = (c & 1) ? 0xedb88320u ^ (c >> 1) : c >> 1;
        crc_table[n] = c;
    }
    crc_ready = 1;
}
static uint32_t crc32_buf(const uint8_t *p, size_t n) {
    if (!crc_ready) crc_init();
    uint32_t c = 0xffffffffu;
    for (size_t i = 0; i < n; i++) c = crc_table[(c ^ p[i]) & 0xff] ^ (c >> 8);
    return c ^ 0xffffffffu;
}
static uint32_t adler32_buf(const uint8_t *p, size_t n) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < n; i++) { a = (a + p[i]) % 65521; b = (b + a) % 65521; }
    return (b << 16) | a;
}

/* ===================================================== DEFLATE tables ==== */

static const uint16_t LEN_BASE[29] = {
    3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258 };
static const uint8_t  LEN_EXTRA[29] = {
    0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0 };
static const uint16_t DIST_BASE[30] = {
    1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,
    1025,1537,2049,3073,4097,6145,8193,12289,16385,24577 };
static const uint8_t  DIST_EXTRA[30] = {
    0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13 };

static uint32_t bitrev(uint32_t v, int n) {
    uint32_t r = 0;
    for (int i = 0; i < n; i++) { r = (r << 1) | (v & 1); v >>= 1; }
    return r;
}

/* ====================================================== bit writer ====== */

typedef struct { uint8_t *buf; size_t len, cap; uint32_t acc; int nbits; int err; } bitw;

static void bw_byte(bitw *w, uint8_t b) {
    if (w->len + 1 > w->cap) {
        size_t nc = w->cap ? w->cap * 2 : 4096;
        uint8_t *nb = realloc(w->buf, nc);
        if (!nb) { w->err = 1; return; }
        w->buf = nb; w->cap = nc;
    }
    w->buf[w->len++] = b;
}
static void bw_bits(bitw *w, uint32_t v, int n) {           /* LSB-first */
    w->acc |= (v & ((1u << n) - 1)) << w->nbits;
    w->nbits += n;
    while (w->nbits >= 8) { bw_byte(w, (uint8_t)(w->acc & 0xff)); w->acc >>= 8; w->nbits -= 8; }
}
static void bw_flush(bitw *w) { if (w->nbits > 0) { bw_byte(w, (uint8_t)(w->acc & 0xff)); w->acc = 0; w->nbits = 0; } }

/* canonical Huffman codes (MSB-first) from code lengths */
static void build_codes(const uint8_t *len, int n, uint16_t *code) {
    int bl_count[16] = {0};
    for (int i = 0; i < n; i++) bl_count[len[i]]++;
    bl_count[0] = 0;
    uint16_t next[16] = {0}; uint16_t c = 0;
    for (int b = 1; b < 16; b++) { c = (c + bl_count[b - 1]) << 1; next[b] = c; }
    for (int i = 0; i < n; i++) if (len[i]) code[i] = next[len[i]]++;
}

/* fixed-Huffman lit/len code lengths */
static void fixed_litlen(uint8_t *len) {
    int i = 0;
    for (; i < 144; i++) len[i] = 8;
    for (; i < 256; i++) len[i] = 9;
    for (; i < 280; i++) len[i] = 7;
    for (; i < 288; i++) len[i] = 8;
}

/* greedy LZ77 + fixed-Huffman DEFLATE (raw stream, no zlib wrapper) */
#define HBITS 15
#define HSIZE (1 << HBITS)
#define WSIZE 32768
#define MINM  3
#define MAXM  258
#define MAXCHAIN 256

static uint32_t hash3(const uint8_t *p) {
    return (((uint32_t)p[0] << 12) ^ ((uint32_t)p[1] << 6) ^ p[2]) & (HSIZE - 1);
}

static uint8_t *deflate_fixed(const uint8_t *src, size_t n, size_t *outlen) {
    bitw w = {0};
    uint8_t llen[288], dlen[30]; uint16_t lcode[288], dcode[30];
    fixed_litlen(llen);
    for (int i = 0; i < 30; i++) dlen[i] = 5;
    build_codes(llen, 288, lcode);
    build_codes(dlen, 30, dcode);

    int *head = malloc(sizeof(int) * HSIZE);
    int *prev = (n ? malloc(sizeof(int) * n) : NULL);
    if (!head || (n && !prev)) { free(head); free(prev); return NULL; }
    for (int i = 0; i < HSIZE; i++) head[i] = -1;

    bw_bits(&w, 1, 1);   /* BFINAL */
    bw_bits(&w, 1, 2);   /* BTYPE = 01 fixed */

    #define EMIT_SYM(s) bw_bits(&w, bitrev(lcode[s], llen[s]), llen[s])
    size_t i = 0;
    while (i < n) {
        int best = 0, bestdist = 0;
        if (i + MINM <= n) {
            uint32_t h = hash3(src + i);
            int cand = head[h], chain = MAXCHAIN;
            long limit = (long)i - WSIZE; if (limit < 0) limit = 0;
            size_t maxl = n - i; if (maxl > MAXM) maxl = MAXM;
            while (cand >= 0 && cand >= limit && chain--) {
                size_t l = 0;
                while (l < maxl && src[cand + l] == src[i + l]) l++;
                if ((int)l > best) { best = (int)l; bestdist = (int)(i - cand); if (l >= maxl) break; }
                cand = prev[cand];
            }
        }
        if (best >= MINM) {
            int li = 28; while (LEN_BASE[li] > best) li--;
            EMIT_SYM(257 + li);
            bw_bits(&w, best - LEN_BASE[li], LEN_EXTRA[li]);
            int di = 29; while (DIST_BASE[di] > bestdist) di--;
            bw_bits(&w, bitrev(dcode[di], dlen[di]), dlen[di]);
            bw_bits(&w, bestdist - DIST_BASE[di], DIST_EXTRA[di]);
            for (int k = 0; k < best; k++) {
                size_t p = i + k;
                if (p + MINM <= n) { uint32_t hh = hash3(src + p); prev[p] = head[hh]; head[hh] = (int)p; }
            }
            i += best;
        } else {
            EMIT_SYM(src[i]);
            if (i + MINM <= n) { uint32_t hh = hash3(src + i); prev[i] = head[hh]; head[hh] = (int)i; }
            i++;
        }
    }
    EMIT_SYM(256);   /* end of block */
    #undef EMIT_SYM
    bw_flush(&w);
    free(head); free(prev);
    if (w.err) { free(w.buf); return NULL; }
    *outlen = w.len;
    return w.buf;
}

/* zlib wrapper around the raw deflate stream */
static uint8_t *zlib_compress(const uint8_t *src, size_t n, size_t *outlen) {
    size_t dl; uint8_t *d = deflate_fixed(src, n, &dl);
    if (!d) return NULL;
    uint8_t *out = malloc(2 + dl + 4);
    if (!out) { free(d); return NULL; }
    out[0] = 0x78; out[1] = 0x9c;                /* CMF, FLG (valid, %31==0) */
    memcpy(out + 2, d, dl);
    uint32_t a = adler32_buf(src, n);
    out[2 + dl + 0] = (a >> 24) & 0xff; out[2 + dl + 1] = (a >> 16) & 0xff;
    out[2 + dl + 2] = (a >> 8) & 0xff;  out[2 + dl + 3] = a & 0xff;
    free(d);
    *outlen = 2 + dl + 4;
    return out;
}

/* ======================================================= INFLATE ======== */

typedef struct { const uint8_t *p; size_t len, pos; uint32_t acc; int nbits; int err; } bitr;

static uint32_t br_bits(bitr *r, int n) {
    while (r->nbits < n) {
        uint8_t b = 0;
        if (r->pos < r->len) b = r->p[r->pos++]; else r->err = 1;
        r->acc |= (uint32_t)b << r->nbits; r->nbits += 8;
    }
    uint32_t v = r->acc & ((1u << n) - 1);
    r->acc >>= n; r->nbits -= n;
    return v;
}

typedef struct { short count[16]; short symbol[288]; } huff;

static void huff_build(huff *h, const uint8_t *len, int n) {
    for (int i = 0; i < 16; i++) h->count[i] = 0;
    for (int i = 0; i < n; i++) h->count[len[i]]++;
    h->count[0] = 0;
    short offs[16]; offs[0] = 0; offs[1] = 0;
    for (int b = 1; b < 15; b++) offs[b + 1] = offs[b] + h->count[b];
    for (int i = 0; i < n; i++) if (len[i]) h->symbol[offs[len[i]]++] = (short)i;
}

static int huff_decode(bitr *r, const huff *h) {
    int code = 0, first = 0, index = 0;
    for (int len = 1; len < 16; len++) {
        code |= (int)br_bits(r, 1);
        int count = h->count[len];
        if (code - first < count) return h->symbol[index + (code - first)];
        index += count; first += count; first <<= 1; code <<= 1;
        if (r->err) return -1;
    }
    return -1;
}

/* INFLATE into a buffer of exactly `outcap` bytes (PNG knows the size). */
static int inflate_raw(const uint8_t *in, size_t inlen, uint8_t *out, size_t outcap) {
    bitr r = { in, inlen, 0, 0, 0, 0 };
    size_t o = 0;
    huff lit, dist;
    int final;
    do {
        final = (int)br_bits(&r, 1);
        int type = (int)br_bits(&r, 2);
        if (r.err) return PNG_EFORMAT;
        if (type == 0) {                                  /* stored */
            r.acc = 0; r.nbits = 0;                        /* align to byte */
            if (r.pos + 4 > r.len) return PNG_EFORMAT;
            uint32_t len = r.p[r.pos] | (r.p[r.pos + 1] << 8); r.pos += 4;
            if (r.pos + len > r.len || o + len > outcap) return PNG_EFORMAT;
            memcpy(out + o, r.p + r.pos, len); r.pos += len; o += len;
            continue;
        } else if (type == 1) {                           /* fixed */
            uint8_t ll[288], dl[30];
            fixed_litlen(ll); for (int i = 0; i < 30; i++) dl[i] = 5;
            huff_build(&lit, ll, 288); huff_build(&dist, dl, 30);
        } else if (type == 2) {                           /* dynamic */
            int hlit = (int)br_bits(&r, 5) + 257;
            int hdist = (int)br_bits(&r, 5) + 1;
            int hclen = (int)br_bits(&r, 4) + 4;
            static const int ord[19] = {16,17,18,0,8,7,9,6,10,5,11,4,12,3,13,2,14,1,15};
            uint8_t cl[19] = {0};
            for (int i = 0; i < hclen; i++) cl[ord[i]] = (uint8_t)br_bits(&r, 3);
            huff clh; huff_build(&clh, cl, 19);
            uint8_t lens[288 + 32] = {0};
            int n = 0, total = hlit + hdist;
            while (n < total) {
                int s = huff_decode(&r, &clh);
                if (s < 0) return PNG_EFORMAT;
                if (s < 16) lens[n++] = (uint8_t)s;
                else if (s == 16) { if (n == 0) return PNG_EFORMAT;
                    int rep = 3 + (int)br_bits(&r, 2); uint8_t pv = lens[n - 1];
                    while (rep-- && n < total) lens[n++] = pv; }
                else if (s == 17) { int rep = 3 + (int)br_bits(&r, 3); while (rep-- && n < total) lens[n++] = 0; }
                else { int rep = 11 + (int)br_bits(&r, 7); while (rep-- && n < total) lens[n++] = 0; }
                if (r.err) return PNG_EFORMAT;
            }
            huff_build(&lit, lens, hlit);
            huff_build(&dist, lens + hlit, hdist);
        } else return PNG_EFORMAT;

        for (;;) {                                        /* decode symbols */
            int s = huff_decode(&r, &lit);
            if (s < 0) return PNG_EFORMAT;
            if (s < 256) { if (o >= outcap) return PNG_EFORMAT; out[o++] = (uint8_t)s; }
            else if (s == 256) break;
            else {
                s -= 257; if (s >= 29) return PNG_EFORMAT;
                int len = LEN_BASE[s] + (int)br_bits(&r, LEN_EXTRA[s]);
                int ds = huff_decode(&r, &dist);
                if (ds < 0 || ds >= 30) return PNG_EFORMAT;
                int d = DIST_BASE[ds] + (int)br_bits(&r, DIST_EXTRA[ds]);
                if ((size_t)d > o || o + len > outcap) return PNG_EFORMAT;
                for (int k = 0; k < len; k++) { out[o] = out[o - d]; o++; }
            }
            if (r.err) return PNG_EFORMAT;
        }
    } while (!final);
    return o == outcap ? PNG_OK : PNG_EFORMAT;
}

static int zlib_decompress(const uint8_t *in, size_t inlen, uint8_t *out, size_t outcap) {
    if (inlen < 6) return PNG_EFORMAT;
    if ((in[0] & 0x0f) != 8) return PNG_EUNSUPPORTED;     /* CM=deflate */
    return inflate_raw(in + 2, inlen - 6, out, outcap);   /* skip 2 hdr, 4 adler */
}

/* ======================================================== filters ======= */

static int paeth(int a, int b, int c) {
    int p = a + b - c, pa = abs(p - a), pb = abs(p - b), pc = abs(p - c);
    if (pa <= pb && pa <= pc) return a;
    return pb <= pc ? b : c;
}

/* ====================================================== chunk I/O ======== */

static int put_chunk(FILE *f, const char *type, const uint8_t *data, size_t n) {
    uint8_t be[4] = { (uint8_t)(n >> 24), (uint8_t)(n >> 16), (uint8_t)(n >> 8), (uint8_t)n };
    if (fwrite(be, 1, 4, f) != 4) return PNG_EIO;
    if (fwrite(type, 1, 4, f) != 4) return PNG_EIO;
    if (n && fwrite(data, 1, n, f) != n) return PNG_EIO;
    uint8_t *tmp = malloc(4 + n);
    if (!tmp) return PNG_EMEM;
    memcpy(tmp, type, 4); if (n) memcpy(tmp + 4, data, n);
    uint32_t crc = crc32_buf(tmp, 4 + n); free(tmp);
    uint8_t cb[4] = { (uint8_t)(crc >> 24), (uint8_t)(crc >> 16), (uint8_t)(crc >> 8), (uint8_t)crc };
    return fwrite(cb, 1, 4, f) == 4 ? PNG_OK : PNG_EIO;
}

static int color_type(int ch) {
    switch (ch) { case 1: return 0; case 2: return 4; case 3: return 2; case 4: return 6; }
    return -1;
}
static int channels_of(int ct) {
    switch (ct) { case 0: return 1; case 4: return 2; case 2: return 3; case 6: return 4; }
    return -1;
}

/* ======================================================== writer ======== */

int png_write(const char *path, const png_image *img) {
    if (!path || !img || !img->data) return PNG_EINVAL;
    if (img->bitdepth != 8 && img->bitdepth != 16) return PNG_EINVAL;
    int ct = color_type(img->channels);
    if (ct < 0) return PNG_EINVAL;
    if (!img->width || !img->height) return PNG_EINVAL;

    size_t W = img->width, H = img->height, ch = img->channels;
    size_t bpp = ch * (img->bitdepth / 8);                /* bytes per pixel */
    size_t rb = W * bpp;                                  /* row bytes */

    /* big-endian raw scanlines */
    uint8_t *raw = malloc(H * rb);
    if (!raw) return PNG_EMEM;
    if (img->bitdepth == 8) {
        memcpy(raw, img->data, H * rb);
    } else {
        const uint16_t *s = img->data;
        for (size_t i = 0; i < H * W * ch; i++) { raw[2 * i] = (uint8_t)(s[i] >> 8); raw[2 * i + 1] = (uint8_t)s[i]; }
    }

    /* adaptive filtering -> filtered stream (1 filter byte + rb per row) */
    uint8_t *filt = malloc(H * (1 + rb));
    uint8_t *line = malloc(rb), *best = malloc(rb);
    if (!filt || !line || !best) { free(raw); free(filt); free(line); free(best); return PNG_EMEM; }
    for (size_t y = 0; y < H; y++) {
        const uint8_t *cur = raw + y * rb;
        const uint8_t *prev = y ? raw + (y - 1) * rb : NULL;
        int bestf = 0; unsigned long bestsum = ~0UL;
        for (int f = 0; f < 5; f++) {
            unsigned long sum = 0;
            for (size_t x = 0; x < rb; x++) {
                int a = x >= bpp ? cur[x - bpp] : 0;
                int b = prev ? prev[x] : 0;
                int c = (prev && x >= bpp) ? prev[x - bpp] : 0;
                int v;
                switch (f) {
                    case 0: v = cur[x]; break;
                    case 1: v = cur[x] - a; break;
                    case 2: v = cur[x] - b; break;
                    case 3: v = cur[x] - ((a + b) >> 1); break;
                    default: v = cur[x] - paeth(a, b, c); break;
                }
                line[x] = (uint8_t)v;
                int sv = (signed char)v; sum += sv < 0 ? -sv : sv;
            }
            if (sum < bestsum) { bestsum = sum; bestf = f; memcpy(best, line, rb); }
        }
        filt[y * (1 + rb)] = (uint8_t)bestf;
        memcpy(filt + y * (1 + rb) + 1, best, rb);
    }
    free(raw); free(line); free(best);

    size_t zl; uint8_t *z = zlib_compress(filt, H * (1 + rb), &zl);
    free(filt);
    if (!z) return PNG_EMEM;

    FILE *f = fopen(path, "wb");
    if (!f) { free(z); return PNG_EOPEN; }
    static const uint8_t sig[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };
    int rc = PNG_OK;
    if (fwrite(sig, 1, 8, f) != 8) { rc = PNG_EIO; goto done; }
    {
        uint8_t ihdr[13];
        ihdr[0] = (uint8_t)(W >> 24); ihdr[1] = (uint8_t)(W >> 16); ihdr[2] = (uint8_t)(W >> 8); ihdr[3] = (uint8_t)W;
        ihdr[4] = (uint8_t)(H >> 24); ihdr[5] = (uint8_t)(H >> 16); ihdr[6] = (uint8_t)(H >> 8); ihdr[7] = (uint8_t)H;
        ihdr[8] = (uint8_t)img->bitdepth; ihdr[9] = (uint8_t)ct;
        ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;         /* deflate, adaptive filter, no interlace */
        if ((rc = put_chunk(f, "IHDR", ihdr, 13)) != PNG_OK) goto done;
    }
    if ((rc = put_chunk(f, "IDAT", z, zl)) != PNG_OK) goto done;
    rc = put_chunk(f, "IEND", NULL, 0);
done:
    free(z);
    if (fclose(f) != 0 && rc == PNG_OK) rc = PNG_EIO;
    return rc;
}

/* ======================================================== reader ======== */

static uint32_t be32(const uint8_t *p) { return ((uint32_t)p[0] << 24) | (p[1] << 16) | (p[2] << 8) | p[3]; }

int png_read(const char *path, png_image *img) {
    if (img) memset(img, 0, sizeof *img);
    if (!path || !img) return PNG_EINVAL;
    FILE *f = fopen(path, "rb");
    if (!f) return PNG_EOPEN;
    fseek(f, 0, SEEK_END); long fsz = ftell(f); fseek(f, 0, SEEK_SET);
    if (fsz < 8) { fclose(f); return PNG_EFORMAT; }
    uint8_t *file = malloc(fsz);
    if (!file) { fclose(f); return PNG_EMEM; }
    if (fread(file, 1, fsz, f) != (size_t)fsz) { free(file); fclose(f); return PNG_EIO; }
    fclose(f);

    static const uint8_t sig[8] = { 137, 80, 78, 71, 13, 10, 26, 10 };
    if (memcmp(file, sig, 8) != 0) { free(file); return PNG_EFORMAT; }

    int rc = PNG_EFORMAT, got_ihdr = 0;
    uint32_t W = 0, H = 0; int bitd = 0, ct = 0, ch = 0;
    uint8_t *idat = NULL; size_t idat_len = 0, idat_cap = 0;

    size_t p = 8;
    while (p + 8 <= (size_t)fsz) {
        uint32_t len = be32(file + p);
        const uint8_t *type = file + p + 4;
        if (p + 12 + len > (size_t)fsz) { rc = PNG_EFORMAT; break; }
        const uint8_t *data = file + p + 8;
        uint32_t crc = be32(file + p + 8 + len);
        if (crc32_buf(file + p + 4, 4 + len) != crc) { rc = PNG_EFORMAT; break; }

        if (!memcmp(type, "IHDR", 4)) {
            if (len != 13) { rc = PNG_EFORMAT; break; }
            W = be32(data); H = be32(data + 4);
            bitd = data[8]; ct = data[9];
            int comp = data[10], filt = data[11], interlace = data[12];
            if (comp != 0 || filt != 0) { rc = PNG_EFORMAT; break; }
            if (interlace != 0) { rc = PNG_EUNSUPPORTED; break; }   /* Adam7 */
            if (bitd != 8 && bitd != 16) { rc = PNG_EUNSUPPORTED; break; }
            ch = channels_of(ct);
            if (ch < 0) { rc = PNG_EUNSUPPORTED; break; }           /* palette etc. */
            if (!W || !H) { rc = PNG_EFORMAT; break; }
            got_ihdr = 1;
        } else if (!memcmp(type, "IDAT", 4)) {
            if (!got_ihdr) { rc = PNG_EFORMAT; break; }
            if (idat_len + len > idat_cap) {
                size_t nc = (idat_len + len) * 2 + 64;
                uint8_t *nb = realloc(idat, nc);
                if (!nb) { rc = PNG_EMEM; break; }
                idat = nb; idat_cap = nc;
            }
            memcpy(idat + idat_len, data, len); idat_len += len;
        } else if (!memcmp(type, "IEND", 4)) {
            rc = got_ihdr ? PNG_OK : PNG_EFORMAT;
            break;
        }
        p += 12 + len;
    }
    free(file);
    if (rc != PNG_OK) { free(idat); return rc; }
    if (!idat_len) { free(idat); return PNG_EFORMAT; }

    size_t bpp = (size_t)ch * (bitd / 8), rb = (size_t)W * bpp;
    uint8_t *filt = malloc(H * (1 + rb));
    if (!filt) { free(idat); return PNG_EMEM; }
    rc = zlib_decompress(idat, idat_len, filt, H * (1 + rb));
    free(idat);
    if (rc != PNG_OK) { free(filt); return rc; }

    /* unfilter in place into a raw scanline buffer */
    uint8_t *raw = malloc(H * rb);
    if (!raw) { free(filt); return PNG_EMEM; }
    for (size_t y = 0; y < H; y++) {
        int ftype = filt[y * (1 + rb)];
        const uint8_t *fl = filt + y * (1 + rb) + 1;
        uint8_t *cur = raw + y * rb;
        const uint8_t *prev = y ? raw + (y - 1) * rb : NULL;
        for (size_t x = 0; x < rb; x++) {
            int a = x >= bpp ? cur[x - bpp] : 0;
            int b = prev ? prev[x] : 0;
            int c = (prev && x >= bpp) ? prev[x - bpp] : 0;
            int v = fl[x];
            switch (ftype) {
                case 0: break;
                case 1: v += a; break;
                case 2: v += b; break;
                case 3: v += (a + b) >> 1; break;
                case 4: v += paeth(a, b, c); break;
                default: free(filt); free(raw); return PNG_EFORMAT;
            }
            cur[x] = (uint8_t)v;
        }
    }
    free(filt);

    /* hand back; convert 16-bit big-endian -> native */
    if (bitd == 8) {
        img->data = raw;
    } else {
        uint16_t *u = malloc(H * (size_t)W * ch * 2);
        if (!u) { free(raw); return PNG_EMEM; }
        for (size_t i = 0; i < H * (size_t)W * ch; i++) u[i] = (uint16_t)((raw[2 * i] << 8) | raw[2 * i + 1]);
        free(raw); img->data = u;
    }
    img->width = W; img->height = H; img->channels = (uint16_t)ch; img->bitdepth = (uint16_t)bitd;
    return PNG_OK;
}

void png_free(png_image *img) {
    if (!img) return;
    free(img->data);
    memset(img, 0, sizeof *img);
}
