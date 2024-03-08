#include <string.h>
#include <zlib.h>

#include "gfa.h"
#include "mgpriv.h"
#include "minigraph.h"

int main(int argc, char *argv[]) {
  char *gfa_path = argv[1];

  // k and w for minimizers
  uint k = 17, w = 11;

  gfa_t *g = gfa_read(gfa_path);
  // gfa_print(g, stdout, 0);

  mg128_v minimizers = {0, 0, 0};
  char mmer[w];
  int i, j;
  uint8_t *c_tag;
  for (i = 0; i < g->n_seg; ++i) {
    gfa_seg_t *s = &g->seg[i];
    if (s->len < k)
      continue;
    for (j = 0; j < s->len; ++j)
      if (s->seq[j] >= 'a' && s->seq[j] <= 'z')
        s->seq[j] -= 32;
    printf("%d %s\n", i, s->name);
    printf("%s\n", s->seq);

    c_tag = gfa_aux_get(s->aux.l_aux, s->aux.aux, "GL");
    // if (c_tag)
    printf("GL: %d\n", *(int32_t *)(c_tag + 1));

    minimizers.n = 0;
    mg_sketch(0, s->seq, s->len, w, k, i,
              &minimizers); // TODO: this can be parallelized

    /*
     * p->a[i].x = kMer << 8 | kmerSpan
     * p->a[i].y = rid << 32 | lastPos << 1 | strand
     *
     * where lastPos is the position of the last base of the i-th minimizer, and
     * strand indicates whether the minimizer comes from the top or the bottom
     * strand.
     */

    printf("#minimizers: %d\n", minimizers.n);
    for (j = 0; j < minimizers.n; ++j) {
      mg128_t minimizer = minimizers.a[j];
      if (!(minimizer.y & 1))
        continue; // FIXME
      printf("minimizer: %d\n", minimizer.x >> 8);
      printf("kmerSpan: %d\n", minimizer.x & 0xFF);

      printf("rid: %d\n", minimizer.y >> 32);
      printf("lastPos: %d\n", minimizer.y >> 1 & 0x7FFFFFFF);
      printf("Strand: %d\n", minimizer.y & 1);

      strncpy(mmer, s->seq + (minimizer.y >> 1 & 0x7FFFFFFF) - w, w);
      printf("mmer: %s\n", mmer);
    }
    printf("\n");
  }

  free(minimizers.a);

  return 0;
}
