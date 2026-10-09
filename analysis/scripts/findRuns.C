#include "clas12reader.h"
#include <cstdio>
using namespace clas12;

void findRuns(const char* fname, int expected = 18540) {
  clas12reader c12(fname);
  long n = 0, nforeign = 0;
  while (c12.next()) {
    int rn = c12.runconfig()->getRun();
    if (rn != expected) {
      printf("entry %ld : run %d, event %d\n", n, rn, c12.runconfig()->getEvent());
      nforeign++;
    }
    n++;
    if (n % 10000000 == 0) fprintf(stderr, "read %ld events\n", n);
  }
  printf("total events read: %ld, foreign events: %ld\n", n, nforeign);
}
