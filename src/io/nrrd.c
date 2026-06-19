/* nrrd.c — see nrrd.h. */
#define _POSIX_C_SOURCE 200809L  // strtok_r
#include "io/nrrd.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

// ---- small helpers ----------------------------------------------------------

static char *read_whole_file(const char *path, size_t *len_out) {
  FILE *f = fopen(path, "rb");
  if (!f) return NULL;
  fseek(f, 0, SEEK_END);
  long sz = ftell(f);
  if (sz < 0) { fclose(f); return NULL; }
  fseek(f, 0, SEEK_SET);
  char *buf = (char *)malloc((size_t)sz);
  if (!buf) { fclose(f); return NULL; }
  size_t got = fread(buf, 1, (size_t)sz, f);
  fclose(f);
  if (got != (size_t)sz) { free(buf); return NULL; }
  *len_out = (size_t)sz;
  return buf;
}

// Find the byte offset just past the header-terminating blank line ("\n\n").
// Returns the offset of the first payload byte, or 0 if not found.
static size_t find_payload_start(const char *buf, size_t len) {
  for (size_t i = 0; i + 1 < len; i++)
    if (buf[i] == '\n' && buf[i + 1] == '\n') return i + 2;
  return 0;
}

static void rstrip(char *s) {
  size_t n = strlen(s);
  while (n && (s[n - 1] == ' ' || s[n - 1] == '\t' || s[n - 1] == '\r' || s[n - 1] == '\n'))
    s[--n] = '\0';
}
static char *lstrip(char *s) {
  while (*s == ' ' || *s == '\t') s++;
  return s;
}

static int parse_type(const char *v, nrrd_type *t, int *sz) {
  if (!strcmp(v, "unsigned char") || !strcmp(v, "uchar") || !strcmp(v, "uint8") ||
      !strcmp(v, "uint8_t")) { *t = NRRD_U8;  *sz = 1; return 0; }
  if (!strcmp(v, "unsigned short") || !strcmp(v, "ushort") || !strcmp(v, "uint16") ||
      !strcmp(v, "uint16_t")) { *t = NRRD_U16; *sz = 2; return 0; }
  if (!strcmp(v, "float")) { *t = NRRD_F32; *sz = 4; return 0; }
  return -1;
}

// gzip/zlib auto-detecting inflate into a pre-sized destination.
static int gunzip(const u8 *src, size_t srclen, u8 *dst, size_t dstlen) {
  z_stream s;
  memset(&s, 0, sizeof s);
  s.next_in = (Bytef *)src;
  s.avail_in = (uInt)srclen;
  s.next_out = dst;
  s.avail_out = (uInt)dstlen;
  if (inflateInit2(&s, 15 + 32) != Z_OK) return -1;  // +32: auto gzip/zlib header
  int r = inflate(&s, Z_FINISH);
  inflateEnd(&s);
  return (r == Z_STREAM_END) ? 0 : -1;
}

// ---- read -------------------------------------------------------------------

int nrrd_read(const char *path, nrrd *out) {
  memset(out, 0, sizeof *out);
  size_t len = 0;
  char *buf = read_whole_file(path, &len);
  if (!buf) return -1;

  size_t payload = find_payload_start(buf, len);
  if (!payload) { free(buf); return -2; }

  // Parse the ASCII header (everything before the blank line).
  int have_type = 0, have_dim = 0, have_sizes = 0;
  char encoding[32] = "raw";
  int dim = 0;
  char *hdr = buf;
  hdr[payload - 1] = '\0';  // terminate at the blank line for safe line walking

  char *save = NULL;
  for (char *line = strtok_r(hdr, "\n", &save); line; line = strtok_r(NULL, "\n", &save)) {
    if (line[0] == '\0' || line[0] == '#') continue;
    if (!strncmp(line, "NRRD", 4)) continue;            // magic
    if (strstr(line, ":=")) continue;                   // key/value pair, ignore
    char *colon = strchr(line, ':');
    if (!colon) continue;
    *colon = '\0';
    char *key = line;
    char *val = lstrip(colon + 1);
    rstrip(key);
    rstrip(val);

    if (!strcmp(key, "type")) {
      if (parse_type(val, &out->type, &out->type_size)) { free(buf); return -3; }
      have_type = 1;
    } else if (!strcmp(key, "dimension")) {
      dim = atoi(val);
      have_dim = 1;
    } else if (!strcmp(key, "sizes")) {
      int i = 0;
      char *s2 = NULL;
      for (char *tok = strtok_r(val, " ", &s2); tok && i < 4; tok = strtok_r(NULL, " ", &s2))
        out->sizes[i++] = atoi(tok);
      out->ndim = i;
      have_sizes = 1;
    } else if (!strcmp(key, "encoding")) {
      strncpy(encoding, val, sizeof encoding - 1);
    } else if (!strcmp(key, "endian")) {
      if (strcmp(val, "little") != 0 && out->type_size > 1) {
        // We only handle little-endian multi-byte data.
        free(buf);
        return -4;
      }
    }
  }

  if (!have_type || !have_dim || !have_sizes || dim != out->ndim || out->ndim < 1) {
    free(buf);
    return -5;
  }

  out->count = 1;
  for (int i = 0; i < out->ndim; i++) out->count *= (size_t)out->sizes[i];
  size_t nbytes = out->count * (size_t)out->type_size;

  out->data = malloc(nbytes);
  if (!out->data) { free(buf); return -6; }

  const u8 *pl = (const u8 *)buf + payload;
  size_t pl_len = len - payload;

  if (!strcmp(encoding, "raw")) {
    if (pl_len < nbytes) { nrrd_free(out); free(buf); return -7; }
    memcpy(out->data, pl, nbytes);
  } else if (!strcmp(encoding, "gzip") || !strcmp(encoding, "gz")) {
    if (gunzip(pl, pl_len, (u8 *)out->data, nbytes)) { nrrd_free(out); free(buf); return -8; }
  } else {
    nrrd_free(out);
    free(buf);
    return -9;  // unsupported encoding (ascii/bzip2/...)
  }

  free(buf);
  return 0;
}

void nrrd_free(nrrd *n) {
  if (n && n->data) { free(n->data); n->data = NULL; }
}

f32 *nrrd_read_f32(const char *path, int *nz, int *ny, int *nx) {
  nrrd n;
  if (nrrd_read(path, &n)) return NULL;
  if (n.ndim != 3) { nrrd_free(&n); return NULL; }

  *nx = n.sizes[0];
  *ny = n.sizes[1];
  *nz = n.sizes[2];

  f32 *vol = (f32 *)malloc(n.count * sizeof(f32));
  if (!vol) { nrrd_free(&n); return NULL; }

  switch (n.type) {
    case NRRD_U8: {
      const u8 *s = (const u8 *)n.data;
      for (size_t i = 0; i < n.count; i++) vol[i] = (f32)s[i];
      break;
    }
    case NRRD_U16: {
      const u16 *s = (const u16 *)n.data;
      for (size_t i = 0; i < n.count; i++) vol[i] = (f32)s[i];
      break;
    }
    case NRRD_F32: {
      memcpy(vol, n.data, n.count * sizeof(f32));
      break;
    }
  }
  nrrd_free(&n);
  return vol;
}

// ---- write (raw) ------------------------------------------------------------

static int nrrd_write_raw(const char *path, const char *type_str, int type_size,
                          const void *data, int nz, int ny, int nx) {
  FILE *f = fopen(path, "wb");
  if (!f) return -1;
  fprintf(f, "NRRD0004\n");
  fprintf(f, "type: %s\n", type_str);
  fprintf(f, "dimension: 3\n");
  fprintf(f, "sizes: %d %d %d\n", nx, ny, nz);  // fastest-first: x y z
  fprintf(f, "encoding: raw\n");
  if (type_size > 1) fprintf(f, "endian: little\n");
  fprintf(f, "\n");  // blank line terminates header
  size_t count = (size_t)nx * (size_t)ny * (size_t)nz;
  size_t wrote = fwrite(data, type_size, count, f);
  fclose(f);
  return wrote == count ? 0 : -2;
}

int nrrd_write_f32(const char *path, const f32 *vol, int nz, int ny, int nx) {
  return nrrd_write_raw(path, "float", 4, vol, nz, ny, nx);
}
int nrrd_write_u8(const char *path, const u8 *vol, int nz, int ny, int nx) {
  return nrrd_write_raw(path, "unsigned char", 1, vol, nz, ny, nx);
}
