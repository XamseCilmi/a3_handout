#include "kdtree.h"
#include "sort.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

struct node {
  int point_index;
  int axis;
  struct node *left;
  struct node *right;
};

struct kdtree {
  int d;
  const double *points;
  struct node* root;
};

static int double_compare(const void *x, const void *y) {
  int g = (*(double*)x > *(double*)y);
  int l = (*(double*)x < *(double*)y);
  return g-l;
}

struct node* kdtree_create_node(int d, const double *points,
                                int depth, int n, int *indexes) {
  struct node m = (struct node){.axis = depth%d, .left=NULL, .right=NULL};
  if(n < 2){
    m.point_index = indexes[0];
  } else {
    double* points_d = malloc(sizeof(double)*n);
    for (int i = 0; i < n; i++){
      points_d[i]=points[indexes[i]*d+m.axis];
    }
    qsort(points_d, n, sizeof(double), double_compare);
    int n_med = n/2+((n/2)<(n-n/2));
    double med_val = points_d[n_med-1];
    free(points_d);
    int* left_indexes = malloc(sizeof(int)*(n_med-1));
    int n_left = 0;
    int* right_indexes = malloc(sizeof(int)*(n-n_med));
    int n_right = 0;
    int median_assigned = 0;
    int *medians = malloc(n*sizeof(int));
    int like_median = 0;
    for (int i = 0; i < n; i++){
      if(points[indexes[i]*d+m.axis]<med_val){
        assert(n_med>1);
        left_indexes[n_left] = indexes[i];
        n_left ++;
      } else if(points[indexes[i]*d+m.axis]>med_val){
        assert((n-n_med)>0);
        right_indexes[n_right] = indexes[i];
        n_right++;
      } else if (!median_assigned) {
        m.point_index = indexes[i];
        median_assigned++;
      } else {
        medians[like_median] = indexes[i];
        like_median++;
      }
    }
    if (like_median){
      while(n_left<(n_med-1)){
        left_indexes[n_left] = medians[like_median];
        like_median--;
        n_left++;
      }
      while(n_right<(n-n_med)){
        right_indexes[n_right] = medians[like_median];
        like_median--;
        n_right++;
      }
    }
    free(medians);
    if (n_left>0){
      m.left = kdtree_create_node(d, points, depth+1, n_left, left_indexes);
    }
    free(left_indexes);
    if (n_right>0){
      m.right = kdtree_create_node(d, points, depth+1, n_right, right_indexes);
    }
    free(right_indexes);
  }
  struct node *mm = malloc(sizeof(struct node));
  mm[0] = m;
  return mm;
}

struct kdtree *kdtree_create(int d, int n, const double *points) {
  struct kdtree *tree = malloc(sizeof(struct kdtree));
  tree->d = d;
  tree->points = points;

  int *indexes = malloc(sizeof(int) * n);

  for (int i = 0; i < n; i++) {
    indexes[i] = i;
  }

  tree->root = kdtree_create_node(d, points, 0, n, indexes);
  free(indexes);

  return tree;
}

void kdtree_free_node(struct node *node) {
  if((node->left) != NULL){
    kdtree_free_node(node->left);
  }
  if((node->right) != NULL){
    kdtree_free_node(node->right);
  }
  free(node);
}

void kdtree_free(struct kdtree *tree) {
  kdtree_free_node(tree->root);
  free(tree);
}

void kdtree_knn_node(const struct kdtree *tree, int k, const double* query,
                     int *closest, double *radius,
                     const struct node *node) {
  if(node==NULL){
    return;
  }
  int updated = insert_if_closer(k, tree->d, tree->points, closest, query, node->point_index);
  double diff = tree->points[(node->point_index)*(tree->d)+(node->axis)] -
          query[(node->axis)];
  if (updated & (closest[k-1] != -1)){
    *radius = distance(tree->d,
                    &(tree->points)[(tree->d)*closest[k-1]], query);
  }
  if ((diff>=0) | (*radius>fabs(diff))){
    kdtree_knn_node(tree, k, query, closest, radius, node->left);
  }
  if ((diff<=0) | (*radius>fabs(diff))){
    kdtree_knn_node(tree, k, query, closest, radius, node->right);
  }
}

int* kdtree_knn(const struct kdtree *tree, int k, const double* query) {
  int* closest = malloc(k * sizeof(int));
  double radius = INFINITY;

  for (int i = 0; i < k; i++) {
    closest[i] = -1;
  }

  kdtree_knn_node(tree, k, query, closest, &radius, tree->root);

  return closest;
}

void kdtree_svg_node(double scale, FILE *f, const struct kdtree *tree,
                     double x1, double y1, double x2, double y2,
                     const struct node *node) {
  if (node == NULL) {
    return;
  }

  double coord = tree->points[node->point_index*2+node->axis];
  if (node->axis == 0) {
    fprintf(f, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"black\" />\n",
            coord*scale, y1*scale, coord*scale, y2*scale);
    kdtree_svg_node(scale, f, tree,
                    x1, y1, coord, y2,
                    node->left);
    kdtree_svg_node(scale, f, tree,
                    coord, y1, x2, y2,
                    node->right);
  } else {
    fprintf(f, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke-width=\"1\" stroke=\"black\" />\n",
            x1*scale, coord*scale, x2*scale, coord*scale);
    kdtree_svg_node(scale, f, tree,
                    x1, y1, x2, coord,
                    node->left);
    kdtree_svg_node(scale, f, tree,
                    x1, coord, x2, y2,
                    node->right);
  }
}

void kdtree_svg(double scale, FILE* f, const struct kdtree *tree) {
  assert(tree->d == 2);
  kdtree_svg_node(scale, f, tree, 0, 0, 1, 1, tree->root);
}
