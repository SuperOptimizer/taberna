/* jpg_test.c — lossy round-trip: error must be small at high quality. */
#include "jpg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int run(int ch, int q, double tol) {
    const uint32_t W = 96, H = 64;
    size_t ns = (size_t)W * H * ch;
    uint8_t *src = malloc(ns);
    for (uint32_t y = 0; y < H; y++)               /* smooth gradients + gentle texture (no wrap) */
        for (uint32_t x = 0; x < W; x++)
            for (int c = 0; c < ch; c++) {
                double v = 128 + 90.0 * sin(x * 0.09 + c) * cos(y * 0.07)
                               + 0.4 * (double)x - 0.3 * (double)y;
                int iv = (int)(v + 0.5); iv = iv < 0 ? 0 : iv > 255 ? 255 : iv;
                src[((size_t)y * W + x) * ch + c] = (uint8_t)iv;
            }

    jpg_image w = { .width = W, .height = H, .channels = (uint16_t)ch, .data = src };
    const char *path = "/tmp/jpg_rt.jpg";
    int rc = jpg_write(path, &w, q);
    if (rc != JPG_OK) { printf("  WRITE ch%d q%d: %s\n", ch, q, jpg_strerror(rc)); free(src); return 1; }

    jpg_image r;
    rc = jpg_read(path, &r);
    if (rc != JPG_OK) { printf("  READ  ch%d q%d: %s\n", ch, q, jpg_strerror(rc)); free(src); return 1; }

    int fail = 0;
    if (r.width != W || r.height != H || r.channels != ch) {
        printf("  META  ch%d q%d mismatch: %ux%u c%u\n", ch, q, r.width, r.height, r.channels); fail = 1;
    } else {
        double se = 0, mae = 0; int mx = 0;
        for (size_t i = 0; i < ns; i++) { int d = (int)((uint8_t*)r.data)[i] - src[i]; se += (double)d*d; mae += abs(d); if (abs(d) > mx) mx = abs(d); }
        mae /= ns; double rmse = sqrt(se / ns);
        printf("  ch%d q%d: RMSE=%.2f MAE=%.2f max=%d  %s\n", ch, q, rmse, mae, mx, mae <= tol ? "ok" : "TOO LOSSY");
        if (mae > tol) fail = 1;
    }
    jpg_free(&r); free(src);
    return fail;
}

int main(void) {
    int fails = 0;
    printf("jpg lossy round-trip:\n");
    fails += run(1, 90, 3.0);
    fails += run(3, 90, 4.0);
    fails += run(1, 75, 6.0);
    fails += run(3, 75, 7.0);
    printf(fails ? "FAIL (%d)\n" : "PASS\n", fails);
    return fails ? 1 : 0;
}
