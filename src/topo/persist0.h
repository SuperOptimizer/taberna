/* persist0.h — fast dimension-0 cubical persistence (union-find merge tree) and
 * persistence-based topological simplification of a scalar field.
 *
 * The generic boundary-matrix engine (cubical.c) is O(cells) and only viable on
 * small/cropped volumes; this is the near-linear path for dim-0 on full 320^3
 * fields. Connectivity matches the oracle (6-conn, sublevel). dim-0 features are
 * the components/basins of the sublevel sets — the merge tree.
 *
 * Application: simplify_dim0 removes shallow basins (spurious minima whose
 * persistence is below tau) by raising them to their merge saddle — principled,
 * scale-aware denoising of the continuous sheetness field that removes
 * topological noise at the source instead of as binary post-hoc cleanup.
 *
 * Volumes z-major, x-fastest.
 */
#ifndef TABERNA_PERSIST0_H
#define TABERNA_PERSIST0_H

#include "common/types.h"
#include "topo/cubical.h"   // pers_pair

/* Dimension-0 persistence pairs (sublevel, 6-conn) via union-find. Off-diagonal
 * pairs match cubical_persistence dim 0; the one essential class (global min) has
 * death = TOPO_INF. Returns malloc'd array; count in *npairs. */
pers_pair *dim0_persistence(const f32 *field, int nz, int ny, int nx, int *npairs);

/* In-place dim-0 persistence simplification: every basin whose persistence
 * (saddle - min) is < tau is raised to its merge saddle, removing it. The result
 * has all dim-0 features of persistence >= tau (the global min survives). */
void simplify_dim0(f32 *field, int nz, int ny, int nx, f32 tau);

#endif // TABERNA_PERSIST0_H
