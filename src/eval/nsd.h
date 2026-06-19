/* nsd.h — Normalized Surface Dice at tolerance, a native port of Google
 * DeepMind's surface_distance (the exact SurfaceDice the official metric uses).
 *
 * Each voxel-corner gets an 8-bit code from its 2x2x2 neighborhood; corners with
 * a mixed code (not all-0, not all-1) are surface points, each weighted by a
 * marching-cubes surfel AREA (256-entry table for unit spacing). NSD@tol =
 * (gt surface area within tol of pred + pred area within tol of gt) / total area,
 * with distances from an exact Euclidean distance transform. Matches
 * sd.compute_surface_dice_at_tolerance for spacing (1,1,1).
 *
 * Masks are u8 (nonzero = foreground), z-major / x-fastest.
 */
#ifndef TABERNA_NSD_H
#define TABERNA_NSD_H

#include "common/types.h"

double surface_dice_nsd(const u8 *gt, const u8 *pred, int nz, int ny, int nx,
                        double tol);

#endif // TABERNA_NSD_H
