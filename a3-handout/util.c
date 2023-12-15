#include "util.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

double distance(int d, const double *x, const double *y) {
  double sum = 0.0;

  for (int i = 0; i < d; i++) {
    sum += (x[i] - y[i]) * (x[i] - y[i]);
  }

  return sqrt(sum);
}

int insert_if_closer(int k, int d,
                     const double *points, int *closest, const double *query,
                     int candidate) {
  double candidate_distance = distance(d, &points[candidate * d], query);
  int farthest_index = 0;
  double farthest_distance = 0.0;

  for (int i = 0; i < k; i++) {
    if (closest[i] != -1) {
      double current_distance = distance(d, &points[i * d], query);
      if (current_distance > farthest_distance) {
        farthest_index = i;
        farthest_distance = current_distance;
      }
    }
    else {
      farthest_index = i;
      break;
    }
  }

  if (candidate_distance < farthest_distance) {
    closest[farthest_index] = candidate;
    return 1;
  }

  return 0;
}