
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "mmio.h"

#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "Maxfiles.h"
#include <MaxSLiCInterface.h>

using namespace std;

typedef float value_t;
typedef uint32_t index_t;
typedef double value_rom_t;

const int numPipes      = BasicCSRPSpMV_numPipes;
const int numCols       = BasicCSRPSpMV_numCols;
const int m_streamLength= 4;

int read_matrix(char *f_name, int *m, int *n, int *nnz, 
    index_t **rows, index_t **cols, value_t **vals) 
{
    printf("reading matrix %s ...\n", f_name);
    int r;
    FILE *f;
    MM_typecode matcode;
    if ((f = fopen(f_name, "r")) == NULL) {
        fprintf(stderr, "error opening file.\n");
        exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0) {
        fprintf(stderr, "error processing matrix banner.\n");
        exit(1);
    }

    if ((r = mm_read_mtx_crd_size(f, m, n, nnz)) != 0) {
        fprintf(stderr, "error reading matrix size info.\n");
        return r;
    }

    *rows = (index_t *) malloc(sizeof(index_t) * (*nnz));
    *cols = (index_t *) malloc(sizeof(index_t) * (*nnz));
    *vals = (value_t *) malloc(sizeof(value_t) * (*nnz));

    for (int i = 0; i < *nnz; i++) {
        fscanf(f, "%d %d %f\n", (*rows)+i, (*cols)+i, (*vals)+i);
        (*vals)[i] = 1.0;
        (*rows)[i] --;
        (*cols)[i] --;
    }
    return 0;
}

// will transform the original COO format into our format
// based on numPipes value.
void transform(int m, int n, int nnz, index_t *rows, index_t *cols, value_t *vals,
    int numPipes, int *_maxColWidth, int *dataSize, index_t **index, value_t **value)
{
    int ceil_m = (int)ceil((double)m/numPipes) * numPipes;
    printf("ceil m: %d\n", ceil_m);
    vector< queue<value_t> > valueQueue(ceil_m);
    vector< queue<index_t> > indexQueue(ceil_m);
    int maxColWidth = 0;
    for (int i = 0; i < nnz; i++) {
        valueQueue[rows[i]].push(vals[i]);
        indexQueue[rows[i]].push(cols[i]);
    }
    for (int i = 0; i < m; i++)
        maxColWidth = max((int)valueQueue[i].size(), maxColWidth);
    
    *dataSize = maxColWidth * ceil_m;
    int streamLength = (int)ceil((double)(*dataSize) / m_streamLength) * m_streamLength;
    *dataSize = streamLength;
    *index = (index_t *) malloc(sizeof(index_t) * streamLength);
    *value = (value_t *) malloc(sizeof(value_t) * streamLength);
   
    // assume n % numPipes == 0
    for (int i = 0; i < ceil_m/numPipes; i ++) {
        for (int j = 0; j < maxColWidth; j++) {
            for (int k = 0; k < numPipes; k++) {
                int idx = i * (numPipes * maxColWidth) + j * numPipes + k;
                int rowIdx = i * numPipes + k;
                (*index)[idx] = (indexQueue[rowIdx].empty()) ? 0 : indexQueue[rowIdx].front();
                (*value)[idx] = (valueQueue[rowIdx].empty()) ? 0.0 : valueQueue[rowIdx].front();

                if (!indexQueue[rowIdx].empty()) {
                    indexQueue[rowIdx].pop();
                    valueQueue[rowIdx].pop();
                }
            }
        }
    }
    *_maxColWidth = maxColWidth;
    printf("finished");
}

void generate_vector(int n, value_rom_t **vector) 
{
    *vector = (value_rom_t *) malloc(sizeof(value_rom_t) * n);
    for (int i = 0; i < n; i++)
        (*vector)[i] = (value_rom_t) i;
}

void basic_spmv(int m, int n, int nnz, index_t *rows, index_t *cols, 
    value_t *vals, value_rom_t *vector, value_t **result)
{
    (*result) = (value_t *) malloc(sizeof(value_t) * m);
    memset(*result, 0, sizeof(value_t) * m);
    for (int i = 0; i < nnz; i++)
        (*result)[rows[i]] += vector[cols[i]] * vals[i];
}

void csrp_spmv(int m, int n, int dataSize, int numPipes, index_t *index, value_t *value, 
    value_rom_t *vector, value_t **result)
{
    (*result) = (value_t *) malloc(sizeof(value_t) * m);
    memset(*result, 0, sizeof(value_t) * m);
    
    int maxColWidth = dataSize / m;
    int startRowIdx = 0;
    for (int i = 0; i < dataSize; i++) {
        if (i != 0 && i % (maxColWidth * numPipes) == 0) {
            startRowIdx += numPipes;
        }
        int biasInPipes = i % numPipes;
        int rowIdx = startRowIdx + biasInPipes;
        (*result)[rowIdx] += vector[index[i]] * value[i];
    }
}

int main(int argc, char *argv[]) 
{
    //if (argc <= 1) {
    //    fprintf(stderr, "usage: %s <matrix name>\n", argv[0]);
    //    exit(1);
    //}

    char *f_name = "../matrix/examples/gemat11/gemat11.mtx";
    int r;

    int m, n, nnz;
    value_t *vals;
    index_t *rows, *cols;
    if ((r = read_matrix(f_name, &m, &n, &nnz, &rows, &cols, &vals)) < 0) {
        fprintf(stderr, "error reading matrix %s, ret %d.\n", f_name, r);
        exit(1);
    }

    printf("\nMatrix m %d n %d nnz %d\n", m, n, nnz);
    
    value_rom_t *vector;
    value_t *expected;
    generate_vector(n, &vector);
    basic_spmv(m, n, nnz, rows, cols, vals, vector, &expected);

    printf("transforming ...\n");

    int numPipes = BasicCSRPSpMV_numPipes;
    int dataSize;
    int maxColWidth;
    index_t *index;
    value_t *value;
    transform(m, n, nnz, rows, cols, vals, numPipes, &maxColWidth, &dataSize, &index, &value);

    printf("transform result:\n");
    printf("max column width: %d\n", maxColWidth);
    printf("dataSize: %d\n", dataSize);

    assert(BasicCSRPSpMV_numCols >= n);
    //for (int i = 0; i < dataSize; i++) {
    //    printf("(%d, %f)\n", index[i], value[i]);
    //}

    value_t *result;
    csrp_spmv(m, n, dataSize, numPipes, index, value, vector, &result);

    printf("Comparing CSRP software version with primitive version ...\n");
    for (int i = 0; i < m; i++) {
        if (abs(expected[i]-result[i]) > 1e-6)
            printf("[%4d] %f %f\n", i, expected[i], result[i]);
    }
    printf("Initializing MaxJ DFE running environment...\n");
    value_t *output_dfe = (value_t *) malloc(sizeof(value_t) * dataSize);
    value_t *result_dfe = (value_t *) malloc(sizeof(value_t) * m);
    
    BasicCSRPSpMV_actions_t actions;
    actions.param_length = dataSize;
    actions.param_maxColWidth = maxColWidth;
    actions.instream_index = index;
    actions.instream_input = value;
    actions.outstream_output = output_dfe;
    
    value_rom_t **valueRom = (value_rom_t **) &(actions.inmem_BasicCSRPSpMVKernel_vectorRom0000);
    int maxRomBytes = sizeof(value_rom_t) * BasicCSRPSpMV_numCols;
    for (int i = 0; i < numPipes; i++) {
        valueRom[i] = (value_rom_t *) malloc(maxRomBytes);
        memcpy(valueRom[i], vector, maxRomBytes);
        //for (int j = 0; j < n; j++)
        //    printf("valueRom[%d][%3d] = %f\n", i, j, valueRom[i][j]);
    }

    max_file_t *mf =  BasicCSRPSpMV_init();
    max_engine_t *me = max_load(mf, "*");
    
    printf("Running DFE.\n");
    struct timeval t0, t1;
    gettimeofday(&t0, 0);
    BasicCSRPSpMV_run(me, &actions);
    gettimeofday(&t1, 0);
    double duration = (double)(t1.tv_sec-t0.tv_sec)+(double)(t1.tv_usec-t0.tv_usec)/1e6;
    printf("Total time %lf s per element time %lf ms %lf GFlops\n", 
        duration,
        duration * 1e6 / nnz,
        (double)(2 * nnz) / duration / 1e9);

    max_unload(me);

    //for (int i = 0; i < dataSize; i++) {
    //    printf("output_dfe[%4d] = %f\n", i, output_dfe[i]);
    //}

    printf("transform output_dfe to result_dfe...\n");
    for (int i = 0; i < m; i++) {
        int blockIdx = i / numPipes;
        int blockBias = i % numPipes;
        int outputIdx = blockIdx * maxColWidth * numPipes + (maxColWidth-1) * numPipes + blockBias;
        result_dfe[i] = output_dfe[outputIdx];
    }

    for (int i = 0; i < m; i++) {
        if (abs(expected[i]-result_dfe[i])/expected[i] >= 1e-5) {
            printf("\x1B[31mERROR\x1B[0m\n");
            printf("expected[%4d] = %f %f\n", i, expected[i], result_dfe[i]);
            exit(1);
        }
    }
    printf("\x1B[32mEverything is fine.\x1B[0m\n\n");

    return 0;
}