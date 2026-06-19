/* partition.h — Mutex Watershed (Wolf et al. 2018) on a signed graph.
 *
 * Parameter-free signed-graph partitioning: process edges by descending |weight|;
 * an attractive edge (w>0) merges its two clusters unless a mutex (cannot-link)
 * forbids it; a repulsive edge (w<0) plants a mutex between its clusters. The
 * repulsion is exactly what keeps touching wraps in separate segments — the
 * failure mode (whole scroll collapsing into one label) that unsigned
 * agglomeration cannot prevent. Near-linearithmic via union-find + mutex lists.
 */
#ifndef TABERNA_PARTITION_H
#define TABERNA_PARTITION_H

#include "common/types.h"
#include "segmentation/affinity.h"

// Partition g into clusters. `labels` (caller-allocated, size g->nnodes) receives
// a dense cluster id per node. Returns the number of clusters, or -1 on error.
int mws_partition(const sgraph *g, u32 *labels);

#endif // TABERNA_PARTITION_H
