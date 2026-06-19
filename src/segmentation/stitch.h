/* stitch.h — cross-brick label stitching for the block-wise pipeline
 * (docs/unwrapping-plan.md §2 stage 5). Two adjacent bricks are processed
 * independently (their local cluster ids are unrelated); in the shared halo
 * overlap the same physical voxel carries a label from each brick. We match those
 * labels (majority co-occurrence) and union them into a global labeling.
 *
 * This is the delayed-merge primitive: stitch resolves identity at brick
 * boundaries so the global result is independent of the chunking. v0 does
 * pairwise majority matching over the overlap; the full Lu/Zlateski/Seung
 * recursive-doubling schedule composes these pairwise stitches.
 */
#ifndef TABERNA_STITCH_H
#define TABERNA_STITCH_H

#include <stddef.h>
#include "common/types.h"

// Given two label arrays over the SAME overlap region (same length, voxel-aligned),
// produce a remap so that brick-B labels are merged with the brick-A labels they
// coincide with. `remap_b` (size nb_labels) maps each B label to a global id; A
// labels keep their ids and B labels that don't match any A label get fresh ids
// starting at na_labels. Returns the total number of global labels.
//   a_over, b_over : label per overlap voxel (length `n`)
//   na_labels, nb_labels : (max label + 1) in A and B respectively
int stitch_overlap(const u32 *a_over, const u32 *b_over, size_t n,
                   int na_labels, int nb_labels, u32 *remap_b);

#endif // TABERNA_STITCH_H
