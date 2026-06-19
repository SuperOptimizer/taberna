/* mimalloc_tune.c — when taberna is built against mimalloc (TABERNA_MIMALLOC=ON),
 * the library is linked into every tool via taberna_io and its malloc/free override
 * the system allocator automatically. This constructor runs before main() and sets
 * mimalloc options programmatically, so behavior doesn't depend on MIMALLOC_* env.
 *
 * Why these: the default (purge_delay ~10ms) is fast but holds reserved segments,
 * inflating RSS between pipeline phases. purge_delay=0 returns memory eagerly but
 * pays a syscall per free and measured ~50% slower. 50ms is the sweet spot: freed
 * segments go back to the OS promptly (bounded RSS) without per-free churn.
 */
#if TABERNA_HAVE_MIMALLOC
#include <mimalloc.h>

__attribute__((constructor)) static void taberna_mimalloc_tune(void) {
  mi_option_set(mi_option_purge_delay, 50);   // ms before returning freed mem to OS
}
#endif
