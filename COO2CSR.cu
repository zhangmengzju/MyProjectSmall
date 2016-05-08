#include <cuda_runtime.h>
#include <iostream>
#include "cusparse.h"
using namespace std;
   
int SpasrseCOO2CSR(
		const int *h_cooRowInd,
		int nnz, 
		int m, 
		int *h_csrRowPtr
		)
{
	////start1////
	cusparseHandle_t handle = 0;
	cusparseStatus_t cusparseStatus = cusparseCreate(&handle);
	////end1////

	////start2////
	int *d_cooRowInd;
	cudaMalloc((void**)&d_cooRowInd, sizeof(int)*(nnz));
	cudaMemcpy(d_cooRowInd, h_cooRowInd, (nnz)*sizeof(int), cudaMemcpyHostToDevice);

	int *d_csrRowPtr;
	cudaMalloc((void**)&d_csrRowPtr, sizeof(int)*(m+1));
	cudaMemcpy(d_csrRowPtr, h_csrRowPtr, (m+1)*sizeof(int), cudaMemcpyHostToDevice);
	////end2////
		
	cusparseXcoo2csr(handle, 	//cusparseHandle_t
					d_cooRowInd,
					nnz, 
					m, 
					d_csrRowPtr, 
					CUSPARSE_INDEX_BASE_ZERO);	//cusparseIndexBase_t			
					 
	////start3//// 
	cudaMemcpy(h_csrRowPtr, d_csrRowPtr, (m+1)*sizeof(int), cudaMemcpyDeviceToHost); 
	cusparseDestroy(handle);	 
	cudaFree(d_cooRowInd); 	 
	cudaFree(d_csrRowPtr); 
	cudaDeviceReset();
	////end3////
	
	return 1;
}

int main()
{
	float* cooVal = new float[9];
    int* cooCol = new int[9]; 	
    int* cooRowInd = new int[9];

    float* pCooVal = cooVal;
    int* pCooCol = cooCol;	
    int* pCooRowInd = cooRowInd;

    *pCooVal = 1.0F; pCooVal++;
    *pCooVal = 4.0F; pCooVal++;
    *pCooVal = 2.0F; pCooVal++;
	*pCooVal = 3.0F; pCooVal++;
	*pCooVal = 5.0F; pCooVal++;
	*pCooVal = 7.0F; pCooVal++;
	*pCooVal = 8.0F; pCooVal++;
	*pCooVal = 9.0F; pCooVal++;
	*pCooVal = 6.0F; pCooVal++;

	*pCooRowInd = 0; pCooRowInd++;
	*pCooRowInd = 0; pCooRowInd++;
	*pCooRowInd = 1; pCooRowInd++;
	*pCooRowInd = 1; pCooRowInd++;
	*pCooRowInd = 2; pCooRowInd++;
	*pCooRowInd = 2; pCooRowInd++;
	*pCooRowInd = 2; pCooRowInd++;
	*pCooRowInd = 3; pCooRowInd++;
	*pCooRowInd = 3; pCooRowInd++;

	*pCooCol = 0; pCooCol++;
	*pCooCol = 1; pCooCol++;
	*pCooCol = 1; pCooCol++;
	*pCooCol = 2; pCooCol++;
	*pCooCol = 0; pCooCol++;
	*pCooCol = 3; pCooCol++;
	*pCooCol = 4; pCooCol++;
	*pCooCol = 2; pCooCol++;
	*pCooCol = 4; pCooCol++;
	//以上代码在给矩阵A赋值，使用COO格式
	
	int *csrRowPtr = new int[5];//m+1	

	int c = SpasrseCOO2CSR(
		cooRowInd,//const int *h_cooRowInd
		9,//int nnz,
		4,//int m,		
		csrRowPtr//int *h_csrRowPtr
	   );
	   
	for(int i=0; i<5; i++)
		cout<< csrRowPtr[i]<<endl;	   
	cout<< "COO2CSR Result" << c << endl;  
	   
	return 0;
} 
