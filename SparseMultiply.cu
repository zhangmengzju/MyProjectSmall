#include <cuda_runtime.h>
#include <iostream>
#include "cusparse.h"
using namespace std;

int SparseMultiply(int m,int n, int k,int nnzA,int nnzB,float *h_A,int *h_RowA,int *h_ColA,float *h_B,int *h_RowB,int *h_ColB,float *h_C,int *h_RowC,int *h_ColC)
{
	int baseC,nnzC;
	// nnzTotalDevHostPtr points to host memory
	int *nnzTotalDevHostPtr = &nnzC;

	cusparseHandle_t handle=0;
	cusparseStatus_t cusparseStatus;

	cusparseMatDescr_t descrA=0;
	cusparseMatDescr_t descrB=0;
	cusparseMatDescr_t descrC=0;

	float *d_A;
	int *d_RowA;
	int *d_ColA;

	float *d_B;
	int *d_RowB;
	int *d_ColB;

	float *d_C;
	int *d_RowC;
	int *d_ColC;


	cusparseStatus = cusparseCreate(&handle);

	cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);

	cudaMalloc((void**)&d_A, sizeof(float)*(nnzA));
	cudaMalloc((void**)&d_RowA, sizeof(int)*(m+1));
	cudaMalloc((void**)&d_ColA, sizeof(int)*(nnzA));

	cudaMalloc((void**)&d_B, sizeof(float)*(nnzB));
	cudaMalloc((void**)&d_RowB, sizeof(int)*(k+1));
	cudaMalloc((void**)&d_ColB, sizeof(int)*(nnzB));


	cudaMemcpy(d_A, h_A, nnzA*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_RowA, h_RowA, (m+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ColA, h_ColA, nnzA*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(d_B, h_B, nnzB*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_RowB, h_RowB, (k+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ColB, h_ColB, nnzB*sizeof(int), cudaMemcpyHostToDevice);

	cusparseStatus = cusparseCreateMatDescr(&descrA);
	cusparseSetMatType(descrA,CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(descrA,CUSPARSE_INDEX_BASE_ZERO);

	cusparseStatus = cusparseCreateMatDescr(&descrB);
	cusparseSetMatType(descrB,CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(descrB,CUSPARSE_INDEX_BASE_ZERO);

	cusparseStatus = cusparseCreateMatDescr(&descrC);
	cusparseSetMatType(descrC,CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(descrC,CUSPARSE_INDEX_BASE_ZERO);
	//////////////////////////////////////////////////////////////////////////
	cudaMalloc((void**)&d_RowC, sizeof(int)*(m+1));

	cusparseXcsrgemmNnz(
		handle,
		CUSPARSE_OPERATION_NON_TRANSPOSE, 
		CUSPARSE_OPERATION_NON_TRANSPOSE, 
		m,
		n,
		k,
		descrA,
		nnzA,
		d_RowA,
		d_ColA,
		descrB,
		nnzB,
		d_RowB,
		d_ColB,
		descrC,
		d_RowC,
		nnzTotalDevHostPtr);

	if (NULL != nnzTotalDevHostPtr)
	{
		nnzC = *nnzTotalDevHostPtr;
	} else {
		cudaMemcpy(&nnzC, d_RowC+m, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&baseC, d_RowC, sizeof(int), cudaMemcpyDeviceToHost);
		nnzC -= baseC;
	}

	
	cudaMalloc((void**)&d_ColC, sizeof(int)*nnzC);
	cudaMalloc((void**)&d_C, sizeof(float)*nnzC);

	
	cusparseScsrgemm(
		handle, 
		CUSPARSE_OPERATION_NON_TRANSPOSE, 
		CUSPARSE_OPERATION_NON_TRANSPOSE, 
		m, 
		n, 
		k,
		descrA, 
		nnzA,
		d_A, 
		d_RowA, 
		d_ColA,
		descrB,
		nnzB,
		d_B, 
		d_RowB, 
		d_ColB,
		descrC,
		d_C, 
		d_RowC, 
		d_ColC);


	cudaMemcpy(h_C, d_C, nnzC*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_RowC, d_RowC, (m+1)*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_ColC, d_ColC, nnzC*sizeof(int), cudaMemcpyDeviceToHost);

	cusparseDestroy(handle);

	cudaFree(d_A);
	cudaFree(d_RowA);
	cudaFree(d_ColA);

	cudaFree(d_B);
	cudaFree(d_RowB);
	cudaFree(d_ColB);

	cudaFree(d_C);
	cudaFree(d_RowC);
	cudaFree(d_ColC);

	cudaDeviceReset();
	return 1;
}

int main()
{
	float* ValA = new float[9];
    int* RowA = new int[5];
    int* ColA = new int[9]; 

    float* pValA = ValA;
    int* pRowA = RowA;
    int* pColA = ColA;

    *pValA = 1.0F; pValA++;
    *pValA = 4.0F; pValA++;
    *pValA = 2.0F; pValA++;
	*pValA = 3.0F; pValA++;
	*pValA = 5.0F; pValA++;
	*pValA = 7.0F; pValA++;
	*pValA = 8.0F; pValA++;
	*pValA = 9.0F; pValA++;
	*pValA = 6.0F; pValA++;

	*pRowA = 0; pRowA++;
	*pRowA = 2; pRowA++;
	*pRowA = 4; pRowA++;
	*pRowA = 7; pRowA++;
	*pRowA = 9; pRowA++;

	*pColA = 0; pColA++;
	*pColA = 1; pColA++;
	*pColA = 1; pColA++;
	*pColA = 2; pColA++;
	*pColA = 0; pColA++;
	*pColA = 3; pColA++;
	*pColA = 4; pColA++;
	*pColA = 2; pColA++;
	*pColA = 4; pColA++;
	//以上代码在给矩阵A赋值，使用CSR格式
	
	float* ValB = new float[5];
	int* RowB = new int[6];
	int* ColB = new int[5]; 

	float* pValB = ValB;
	int* pRowB = RowB;
	int* pColB = ColB;

	*pValB = 1.0F; pValB++;
	*pValB = 2.0F; pValB++;
	*pValB = 3.0F; pValB++;
	*pValB = 4.0F; pValB++;
	*pValB = 5.0F; pValB++;

	*pRowB = 0; pRowB++;
	*pRowB = 1; pRowB++;
	*pRowB = 2; pRowB++;
	*pRowB = 3; pRowB++;
	*pRowB = 4; pRowB++;
	*pRowB = 5; pRowB++;
	//以上代码在给矩阵B赋值，使用CSR格式
	
	*pColB = 0; pColB++;
	*pColB = 1; pColB++;
	*pColB = 2; pColB++;
	*pColB = 0; pColB++;
	*pColB = 1; pColB++;

	float* ValC = new float[20];
	int* RowC = new int[10];
	int* ColC = new int[20]; 

	int c = SparseMultiply(
		4,// int m,
		3,//int n,
		5,//int k,
		9,//int nnzA,
		5,//int nnzB,
		ValA,//float *h_A,
		RowA,// int *h_RowA,
		ColA,//int *h_ColA,
		ValB,// float *h_B,
		RowB,//  int *h_RowB,
		ColB,// int *h_ColB,
		ValC,//float *h_C,
		RowC,//int *h_RowC,
		ColC//int *h_ColC)
	   );
	   
	cout<< "SparseMultiply Result" << c << endl;  
	   
	return 0;
} 
