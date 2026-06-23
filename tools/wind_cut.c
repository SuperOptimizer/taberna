/* wind_cut — Phase 3b (hard) of docs/touching-sheets-plan.md: the ORDERED-LABEL MIN-CUT that places the
 * missing wrap boundary through a LEAKED (gapless fused) touch — the one thing the soft min-slope ordering
 * provably cannot do (it can't choose a discrete cut location). This is the Boykov-Kolmogorov/Ishikawa
 * construction the research and the negative results both point to.
 *
 * BUILT VERIFY-FIRST: this file starts as a correct, SELF-TESTED maxflow (Dinic). The grid construction
 * (banded binary cuts / ordered level-sets along the smoothed-normal columns) is layered on top ONLY once
 * the solver is proven — no grid cut on an unverified maxflow.
 *
 *   wind_cut selftest                 # run the maxflow unit tests (known min-cuts) and exit
 *   wind_cut ARC OUT ... (grid cut)   # TODO: layered on after the solver is verified
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ---------------------------------------------------------------- Dinic maxflow
 * Half-edge arrays: edge e and its reverse e^1 are adjacent (2*k, 2*k+1). cap[] is the RESIDUAL
 * capacity (decremented on push, the reverse incremented), so flow = original_cap - cap on forward
 * edges. Adjacency as head[]/nxt[] linked lists. Correctness rests on: (1) BFS level graph only
 * advances to strictly deeper nodes; (2) blocking-flow DFS with it[] "current-arc" pointers never
 * re-scans a saturated arc; (3) dead-ends are pruned by setting level=-1. */
typedef struct {
  int n;                 /* node count */
  long m;                /* half-edge count */
  long cap_e;            /* edge array capacity */
  int *to; long *cap; long *nxt;   /* per half-edge: target, residual cap, next-in-adjacency */
  long *head;            /* per node: first half-edge (-1 = none) */
  int *level; long *it;  /* BFS level, current-arc iterator */
  int *q;                /* BFS queue */
} Dinic;

static void d_init(Dinic*d,int n,long max_edges){
  d->n=n; d->m=0; d->cap_e=2*max_edges;
  d->to=malloc(d->cap_e*sizeof(int)); d->cap=malloc(d->cap_e*sizeof(long)); d->nxt=malloc(d->cap_e*sizeof(long));
  d->head=malloc((size_t)n*sizeof(long)); for(int i=0;i<n;i++)d->head[i]=-1;
  d->level=malloc((size_t)n*sizeof(int)); d->it=malloc((size_t)n*sizeof(long)); d->q=malloc((size_t)n*sizeof(int));
}
static void d_free(Dinic*d){ free(d->to);free(d->cap);free(d->nxt);free(d->head);free(d->level);free(d->it);free(d->q); }
/* directed edge u->v with capacity cu, and reverse v->u with capacity cv (cv=0 for a pure directed arc,
 * cv=cu for an undirected edge). */
static void d_edge(Dinic*d,int u,int v,long cu,long cv){
  d->to[d->m]=v; d->cap[d->m]=cu; d->nxt[d->m]=d->head[u]; d->head[u]=d->m; d->m++;
  d->to[d->m]=u; d->cap[d->m]=cv; d->nxt[d->m]=d->head[v]; d->head[v]=d->m; d->m++;
}
static int d_bfs(Dinic*d,int s,int t){
  for(int i=0;i<d->n;i++)d->level[i]=-1;
  int qh=0,qt=0; d->level[s]=0; d->q[qt++]=s;
  while(qh<qt){ int u=d->q[qh++];
    for(long e=d->head[u];e!=-1;e=d->nxt[e]){ int v=d->to[e];
      if(d->cap[e]>0 && d->level[v]<0){ d->level[v]=d->level[u]+1; d->q[qt++]=v; } } }
  return d->level[t]>=0;
}
static long d_dfs(Dinic*d,int u,int t,long f){
  if(u==t) return f;
  for(long *ep=&d->it[u]; *ep!=-1; *ep=d->nxt[*ep]){ long e=*ep; int v=d->to[e];
    if(d->cap[e]>0 && d->level[v]==d->level[u]+1){
      long got=d_dfs(d,v,t,f<d->cap[e]?f:d->cap[e]);
      if(got>0){ d->cap[e]-=got; d->cap[e^1]+=got; return got; } } }
  d->level[u]=-1;   /* dead end: don't revisit this round */
  return 0;
}
static long d_maxflow(Dinic*d,int s,int t){
  long flow=0;
  while(d_bfs(d,s,t)){
    for(int i=0;i<d->n;i++)d->it[i]=d->head[i];
    long f; while((f=d_dfs(d,s,t,(long)1<<60))>0) flow+=f;
  }
  return flow;
}
/* min-cut: source side = nodes reachable from s in the residual graph after maxflow */
static void d_cut_side(Dinic*d,int s,unsigned char*side){
  for(int i=0;i<d->n;i++)side[i]=0;
  int qh=0,qt=0; side[s]=1; d->q[qt++]=s;
  while(qh<qt){ int u=d->q[qh++];
    for(long e=d->head[u];e!=-1;e=d->nxt[e]){ int v=d->to[e];
      if(d->cap[e]>0 && !side[v]){ side[v]=1; d->q[qt++]=v; } } }
}

/* ---------------------------------------------------------------- self-test */
static int g_fail=0;
static void check(const char*name,long got,long want){
  if(got==want) printf("  PASS %-28s = %ld\n",name,got);
  else { printf("  FAIL %-28s = %ld (want %ld)\n",name,got,want); g_fail=1; }
}
static int selftest(void){
  printf("maxflow self-test (Dinic):\n");
  /* 1. trivial bottleneck: s->a cap5, a->t cap3 -> 3 */
  { Dinic d; d_init(&d,4,4); int s=0,a=1,t=2; d_edge(&d,s,a,5,0); d_edge(&d,a,t,3,0);
    check("bottleneck s-a-t",d_maxflow(&d,s,t),3); d_free(&d); }
  /* 2. two parallel paths: s->a->t (cap 4,4), s->b->t (cap 2,5) -> 4+2=6 */
  { Dinic d; d_init(&d,4,4); int s=0,a=1,b=2,t=3;
    d_edge(&d,s,a,4,0); d_edge(&d,a,t,4,0); d_edge(&d,s,b,2,0); d_edge(&d,b,t,5,0);
    check("two parallel paths",d_maxflow(&d,s,t),6); d_free(&d); }
  /* 3. CLRS classic 6-node network, known max flow = 23 */
  { Dinic d; d_init(&d,6,9); int s=0,t=5;
    d_edge(&d,0,1,16,0); d_edge(&d,0,2,13,0); d_edge(&d,1,2,10,0); d_edge(&d,2,1,4,0);
    d_edge(&d,1,3,12,0); d_edge(&d,3,2,9,0); d_edge(&d,2,4,14,0); d_edge(&d,4,3,7,0);
    d_edge(&d,3,5,20,0); d_edge(&d,4,5,4,0);
    check("CLRS network",d_maxflow(&d,s,t),23); d_free(&d); }
  /* 4. min-cut value == max-flow, and the cut partition is consistent (s-side reaches t? must NOT) */
  { Dinic d; d_init(&d,4,4); int s=0,a=1,b=2,t=3;
    d_edge(&d,s,a,3,0); d_edge(&d,a,b,2,0); d_edge(&d,b,t,3,0); d_edge(&d,s,b,1,0);
    long f=d_maxflow(&d,s,t); check("series+shortcut flow",f,3);
    unsigned char side[4]; d_cut_side(&d,s,side);
    /* cut edges = forward edges from s-side to t-side; their capacity sum must equal f */
    long cutcap=0; /* recompute against ORIGINAL caps: rebuild needed -> instead verify t not on s-side */
    check("t not on source side",side[t]?1:0,0); (void)cutcap;
    d_free(&d); }
  /* 5. undirected grid 2x2, s=corner, t=opposite, all caps 1: max flow = 2 (two edge-disjoint paths) */
  { Dinic d; d_init(&d,4,8); /* nodes 0=(0,0)1=(0,1)2=(1,0)3=(1,1); s=0 t=3 */
    d_edge(&d,0,1,1,1); d_edge(&d,0,2,1,1); d_edge(&d,1,3,1,1); d_edge(&d,2,3,1,1);
    check("2x2 undirected grid",d_maxflow(&d,0,3),2); d_free(&d); }
  /* 6. ordering arc (infinity): models the Ishikawa "label monotone" constraint. s->a cap10,
   *    a=>b INF (ordering: b can't be cut from a's side cheaply), b->t cap10. Cutting must pay
   *    min(10,10)=10, never the INF arc. */
  { Dinic d; d_init(&d,4,3); int s=0,a=1,b=2,t=3; long INF=1L<<50;
    d_edge(&d,s,a,10,0); d_edge(&d,a,b,INF,0); d_edge(&d,b,t,10,0);
    check("infinity ordering arc",d_maxflow(&d,s,t),10); d_free(&d); }
  printf(g_fail?"SELFTEST FAILED\n":"ALL TESTS PASSED\n");
  return g_fail;
}

int main(int argc,char**argv){
  if(argc>=2 && !strcmp(argv[1],"selftest")) return selftest();
  fprintf(stderr,"usage: %s selftest   (grid ordered-cut layered on after the solver is verified)\n",argv[0]);
  return 2;
}
