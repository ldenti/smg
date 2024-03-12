#include <string.h>
#include <zlib.h>

#include "gfa.h"
#include "kseq.h"
#include "mgpriv.h"
#include "minigraph.h"

KSEQ_INIT(gzFile, gzread)

int main_gfa(int argc, char *argv[]) {
  char *gfa_path = argv[1];

  /* k and w for minimizers */
  uint k = 17, w = 11;

  gfa_t *g = gfa_read(gfa_path);
  /* gfa_print(g, stdout, 0); */

  mg128_v minimizers = {0, 0, 0};
  char mmer[w];
  int i, j;
  uint8_t *c_tag;
  for (i = 0; i < g->n_seg; ++i) {
    gfa_seg_t *s = &g->seg[i];
    if (s->len < k)
      continue; // FIXME skip or not?
    for (j = 0; j < s->len; ++j)
      if (s->seq[j] >= 'a' && s->seq[j] <= 'z')
        s->seq[j] -= 32;
    /* printf("%d %s\n", i, s->name); */
    /* printf("%s\n", s->seq); */

    c_tag = gfa_aux_get(s->aux.l_aux, s->aux.aux, "GL");
    // if (c_tag)
    // printf("GL: %d\n", *(int32_t *)(c_tag + 1));

    minimizers.n = 0;
    mg_sketch(0, s->seq, s->len, w, k, i,
              &minimizers); // TODO this can be parallelized

    /*
     * p->a[i].x = kMer << 8 | kmerSpan
     * p->a[i].y = rid << 32 | lastPos << 1 | strand
     *
     * where lastPos is the position of the last base of the i-th minimizer, and
     * strand indicates whether the minimizer comes from the top or the bottom
     * strand.
     */

    /* printf("#minimizers: %d\n", minimizers.n); */
    for (j = 0; j < minimizers.n; ++j) {
      mg128_t minimizer = minimizers.a[j];
      /* if (!(minimizer.y & 1)) */
      /*   continue; // FIXME */
      printf("%ld %d\n", minimizer.x >> 8, *(int32_t *)(c_tag + 1));
      /* printf("minimizer: %d\n", minimizer.x >> 8); */
      /* printf("kmerSpan: %d\n", minimizer.x & 0xFF); */
      /* printf("rid: %d\n", minimizer.y >> 32); */
      /* printf("lastPos: %d\n", minimizer.y >> 1 & 0x7FFFFFFF); */
      /* printf("Strand: %d\n", minimizer.y & 1); */
      /* // conversion to char */
      /* strncpy(mmer, s->seq + (minimizer.y >> 1 & 0x7FFFFFFF) - w, w); */
      /* printf("mmer: %s\n", mmer); */
    }
  }

  free(minimizers.a);
  /* kh_destroy(hm64, hm); */

  return 0;
}

int main_fx(int argc, char *argv[]) {
  char *fx_fn = argv[1];

  /* k and w for minimizers */
  uint k = 17, w = 11;

  gzFile fp = gzopen(fx_fn, "rb");
  kseq_t *ks = kseq_init(fp);
  int l, i, j;
  mg128_v minimizers = {0, 0, 0};
  while ((l = kseq_read(ks)) >= 0) {
    for (j = 0; j < ks->seq.l; ++j)
      if (ks->seq.s[j] >= 'a' && ks->seq.s[j] <= 'z')
        ks->seq.s[j] -= 32;

    minimizers.n = 0;
    mg_sketch(0, ks->seq.s, ks->seq.l, w, k, i, &minimizers);

    printf("%s", ks->name.s);
    for (j = 0; j < minimizers.n; ++j) {
      mg128_t minimizer = minimizers.a[j];
      /* if (!(minimizer.y & 1)) */
      /*   continue; // FIXME */
      printf(" %d", minimizer.x >> 8);
    }
    printf("\n");
  }

  free(minimizers.a);
  kseq_destroy(ks);
  gzclose(fp);

  return 0;
}

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "gfa") == 0)
    return main_gfa(argc - 1, argv + 1);
  else if (strcmp(argv[1], "fx") == 0)
    return main_fx(argc - 1, argv + 1);
  else
    return 1;
}
