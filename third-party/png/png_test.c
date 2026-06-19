/* png_test.c — round-trip gray/GA/RGB/RGBA at 8 and 16 bit. */
#include "png.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int run(int ch, int bitd) {
    const uint32_t W = 67, H = 41;
    size_t ns = (size_t)W * H * ch;
    size_t bs = bitd / 8;
    void *src = malloc(ns * bs);
    if (bitd == 8) { uint8_t *s = src; for (size_t i = 0; i < ns; i++) s[i] = (uint8_t)((i * 137 + 11) & 0xff); }
    else { uint16_t *s = src; for (size_t i = 0; i < ns; i++) s[i] = (uint16_t)((i * 40503 + 12345) & 0xffff); }

    png_image w = { .width = W, .height = H, .channels = (uint16_t)ch, .bitdepth = (uint16_t)bitd, .data = src };
    const char *path = "/tmp/png_rt.png";
    int rc = png_write(path, &w);
    if (rc != PNG_OK) { printf("  WRITE ch%d %db: %s\n", ch, bitd, png_strerror(rc)); free(src); return 1; }

    png_image r;
    rc = png_read(path, &r);
    if (rc != PNG_OK) { printf("  READ  ch%d %db: %s\n", ch, bitd, png_strerror(rc)); free(src); return 1; }

    int fail = 0;
    if (r.width != W || r.height != H || r.channels != ch || r.bitdepth != bitd) {
        printf("  META  ch%d %db mismatch: %ux%u c%u d%u\n", ch, bitd, r.width, r.height, r.channels, r.bitdepth);
        fail = 1;
    } else if (memcmp(src, r.data, ns * bs) != 0) {
        printf("  DATA  ch%d %db mismatch\n", ch, bitd);
        fail = 1;
    }
    png_free(&r); free(src);
    if (!fail) printf("  ok    ch%d %db\n", ch, bitd);
    return fail;
}

int main(void) {
    int fails = 0;
    printf("png round-trip:\n");
    for (int bitd = 8; bitd <= 16; bitd += 8)
        for (int ch = 1; ch <= 4; ch++)
            fails += run(ch, bitd);
    printf(fails ? "FAIL (%d)\n" : "PASS\n", fails);
    return fails ? 1 : 0;
}
