#include "csc.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "cjson/cJSON.h"

int csc_Nonzeros(const SparseMatrixCSC *A) { return A->colptr[A->n]; }

void csc_FreeSparseMatrixCSC(SparseMatrixCSC *A) {
  if (A->colptr) {
    free(A->colptr);
  }
  if (A->rowval) {
    free(A->rowval);
  }
  if (A->nzval) {
    free(A->nzval);
  }
}

void kkt_FreeKKTSystem(KKTSystem* kkt) {
  csc_FreeSparseMatrixCSC(&kkt->A);
  free(kkt->b);
  free(kkt->x);
}

/**
 * @brief Read the contents of a file into a heap-allocated `char` array.
 *
 * It is the user's responsibility to call `free` on `out` after the string data
 * is allocated.
 *
 * @param filename Name of the file to read.
 * @param out      Pointer to the array (pointer) to the heap-allocated string
 * data.
 * @param len      Length of the string data.
 * @return int     0 if successful, -1 otherwise.
 */
int ReadFile(const char *filename, char **out, int *len) {
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Couldn't open file\n");
    return -1;
  }

  fseek(fp, 0L, SEEK_END);
  int size = ftell(fp);
  rewind(fp);

  char *buf = (char *)malloc(size + 1);
  if (!buf) {
    fclose(fp);
    fprintf(stderr, "Couldn't allocate memory for the file contents.");
    return -1;
  }

  if (1 != fread(buf, size, 1, fp)) {
    fprintf(stderr, "Failed to read the entire file.");
    free(buf);
    fclose(fp);
    return -1;
  }

  // Add null-termination
  buf[size] = '\0';

  *out = buf;
  *len = size;
  fclose(fp);
  return 0;
}

/**
 * @brief Read a JSON array into a vector of doubles.
 * 
 * @param json JSON object containing the array
 * @param name Name of the JSON array
 * @param buf  Storage location for the data
 * @param len  Expected length of the array
 * @return int 0 if successful, -1 otherwise.
 */
int ReadJSONVectord(cJSON* json, const char* name, double* buf, int len) {
  cJSON* array = cJSON_GetObjectItemCaseSensitive(json, name);
  cJSON* number;
  int i = 0;
  if (cJSON_IsArray(array)) {
    cJSON_ArrayForEach(number, array) {
      if (cJSON_IsNumber(number) && i < len) {
        buf[i] = number->valuedouble;
      }
      ++i;
    }
    if (i == len) {
      return 0;
    } else if (i > len) {
      fprintf(stderr, "JSON array was longer than expected (%d instead of %d).\n", i, len);
    } else {
      fprintf(stderr, "JSON array was shorter than expected (%d instead of %d).\n", i, len);
    }
  } else {
    fprintf(stderr, "Couldn't find an array of name %s.\n", name);
  }
  return -1;
}

int ReadJSONVectori(cJSON* json, const char* name, csc_int* buf, int len) {
  cJSON* array = cJSON_GetObjectItemCaseSensitive(json, name);
  cJSON* number;
  int i = 0;
  if (cJSON_IsArray(array)) {
    cJSON_ArrayForEach(number, array) {
      if (cJSON_IsNumber(number) && i < len) {
        buf[i] = number->valueint;
      }
      ++i;
    }
    if (i == len) {
      return 0;
    } else if (i > len) {
      fprintf(stderr, "JSON array was longer than expected (%d instead of %d).\n", i, len);
    } else {
      fprintf(stderr, "JSON array was shorter than expected (%d instead of %d).\n", i, len);
    }
  } else {
    fprintf(stderr, "Couldn't find an array of name %s.\n", name);
  }
  return -1;
}


KKTSystem kkt_ReadFromFile(const char *filename) {
  SparseMatrixCSC Anull = {.n = 0, .colptr = NULL, .rowval = NULL, .nzval = NULL};
  KKTSystem kkt = {
    .nprimals = 0,
    .nduals = 0,
    .A = Anull,
    .b = NULL,
    .x = NULL,
  };

  // Read the file as a string
  char *jsondata = NULL;
  int len;
  if (0 != ReadFile(filename, &jsondata, &len)) {
    fprintf(stderr, "ERROR: Reading CSC file failed.\n");
    return kkt;
  }

  cJSON *json = cJSON_Parse(jsondata);
  free(jsondata);

  if (json == NULL) {
    const char *error_ptr = cJSON_GetErrorPtr();
    if (error_ptr != NULL) {
      fprintf(stderr, "ERROR: Error parsing JSON file: %s\n", error_ptr);
    }
    return kkt;
  }

  // Read Integer Data
  cJSON* item;
  // int nstates = 0;
  // item = cJSON_GetObjectItemCaseSensitive(json, "nx");
  // if (cJSON_IsNumber(item)) {
  //   nstates = item->valueint;
  // }

  // int ninputs = 0;
  // item = cJSON_GetObjectItemCaseSensitive(json, "nu");
  // if (cJSON_IsNumber(item)) {
  //   ninputs = item->valueint;
  // }

  int nprimals = 0;
  item = cJSON_GetObjectItemCaseSensitive(json, "nprimals");
  if (cJSON_IsNumber(item)) {
    nprimals = item->valueint;
  }

  int nduals = 0;
  item = cJSON_GetObjectItemCaseSensitive(json, "nduals");
  if (cJSON_IsNumber(item)) {
    nduals = item->valueint;
  }

  int nnz = 0;
  item = cJSON_GetObjectItemCaseSensitive(json, "nnz");
  if (cJSON_IsNumber(item)) {
    nnz = item->valueint;
  }
  // printf("Read the following data:\n");
  // printf("  nx = %d\n", nstates);
  // printf("  nu = %d\n", ninputs);
  // printf("  nprimals = %d\n", nprimals);
  // printf("  nduals = %d\n", nduals);
  // printf("  nnz = %d\n", nnz);

  // Read arrays
  int n = nprimals + nduals;
  int status = 0;
  csc_int* colptr = (csc_int*) malloc((n + 1) * sizeof(csc_int));
  status = ReadJSONVectori(json, "colptr", colptr, n + 1);
  if (status != 0) {
    fprintf(stderr, "ERROR: Failed to read the colptr data.");
    free(colptr);
    cJSON_Delete(json);
    return kkt;
  }
  if (colptr[n] != nnz) {
    fprintf(stderr, "ERROR: bad terminal colptr data.");
    free(colptr);
    cJSON_Delete(json);
    return kkt;
  }

  csc_int* rowval = (csc_int*) malloc(nnz * sizeof(csc_int));
  status = ReadJSONVectori(json, "rowval", rowval, nnz);
  if (status != 0) {
    fprintf(stderr, "ERROR: Failed to read the rowval data.");
    free(colptr);
    free(rowval);
    return kkt;
  }

  double* nzval = (double*) malloc(nnz * sizeof(double));
  status = ReadJSONVectord(json, "nzval", nzval, nnz);
  if (status != 0) {
    fprintf(stderr, "ERROR: Failed to read the nzval data.");
    free(colptr);
    free(rowval);
    free(nzval);
    cJSON_Delete(json);
    return kkt;
  }

  double* b = (double*) malloc(n * sizeof(double));
  status = ReadJSONVectord(json, "b", b, n);
  if (status != 0) {
    fprintf(stderr, "ERROR: Failed to read the b data.");
    free(colptr);
    free(rowval);
    free(nzval);
    free(b);
    cJSON_Delete(json);
    return kkt;
  }

  double* x = (double*) malloc(n * sizeof(double));
  status = ReadJSONVectord(json, "x", x, n);
  if (status != 0) {
    fprintf(stderr, "ERROR: Failed to read the x data.");
    free(colptr);
    free(rowval);
    free(nzval);
    free(b);
    free(x);
    cJSON_Delete(json);
    return kkt;
  }
  cJSON_Delete(json);

  SparseMatrixCSC A = {
    .n = n,
    .colptr = colptr,
    .rowval = rowval, 
    .nzval = nzval,
  };

  kkt.nprimals = nprimals;
  kkt.nduals = nduals;
  kkt.A = A;
  kkt.b = b;
  kkt.x = x;
  return kkt;
}
