#include <cuda_runtime.h>
#include <iostream>
#include "cusparse.h"
#include <mpi.h>
#include <vector>
#include <fstream>
#include <sstream>//在string后连接int等类型
#include <string> 
#include <stdlib.h>//使能在g++中编译通过atoi
//#include <io.h>//windows下得到目录下的文件信息

#include <dirent.h>//linux下得到目录下的文件信息opendir
#include <sys/stat.h>//linux下得到目录下的文件信息
#include <sys/types.h>//linux创建目录mkdir

#include <algorithm>//sort
#include <iomanip>//cout setprecision
#include <map>
#include <time.h>
#include <errno.h>

using namespace std;//g++中要对string、vector、ifstream、stringstream、endl、cout等显式添加std:: 

/////////////////////class COO begin///////////////////////////////////////////////
class COO {
	public:
		vector<int> coo_rows;//避免像使用数组时，需要先确定大小
		vector<int> coo_cols;
		vector<float> coo_vals;
		int coo_rows_max ;//m=coo_rows_max+1 coo_rows下标的最大值，若下标最大为3，则实际有4行
		int coo_cols_max ;//n=coo_cols_max+1 coo_cols下标的最大值，若下标最大为3，则实际有4列
		
		COO(int default_max) {
			coo_rows_max = default_max;
			coo_cols_max = default_max;
		}
};
/////////////////////class COO end///////////////////////////////////////////////


/////////////////////class CSR begin///////////////////////////////////////////////
class CSR {
	public:
		vector<int> csr_row_ptrs;//避免像使用数组时，需要先确定大小
		vector<int> csr_cols;
		vector<float> csr_vals;
		int csr_rows_max ;//m=csr_rows_max+1 csr_rows下标的最大值，若下标最大为3，则实际有4行
		int csr_cols_max ;//n=csr_cols_max+1 csr_cols下标的最大值，若下标最大为3，则实际有4列
		
		CSR( int default_max) {
			csr_rows_max = default_max;
			csr_cols_max = default_max;
		}


	//	!!!!!!!!!!!!!!!!!!!!!!
    /////////////////LongID与IntID的映射表
	//	map<long,int> Long2IntMap;
    //	vector<long> Int2LongVector;
    //	!!!!!!!!!!!!!!!!!!!!!!

};
/////////////////////class CSR end///////////////////////////////////////////////



/////////////////////class CSC begin///////////////////////////////////////////////
class CSC {
        public:
                vector<int> csc_col_ptrs;//避免像使用数组时，需要先确定大小
                vector<int> csc_rows;
                vector<float> csc_vals;
                int csc_rows_max ;//m=csc_rows_max+1 csc_rows下标的最大值，若下标最大为3，则实际有4行
                int csc_cols_max ;//n=csc_cols_max+1 csc_cols下标的最大值，若下标最大为3，则实际有4列

                CSC( int default_max) {
                        csc_rows_max = default_max;
                        csc_cols_max = default_max;
                }
		
		CSR use_CSC_Create_CSR() {
			CSR csr(-1);
                        csr.csr_rows_max = csc_rows_max;
                        csr.csr_cols_max = csc_cols_max;

			csr.csr_row_ptrs = csc_col_ptrs;
			csr.csr_cols = csc_rows;
			csr.csr_vals = csc_vals;

			return csr;
		}
};
/////////////////////class CSC end///////////////////////////////////////////////


/*////////////////////////////////////////////////////////////////////////////////
//////////////////saveCSRAsFiles begin////////////////////////////////////////////////
void saveCSRAsFiles(string fileNameDir, CSR csr) {

        //目录不为空则删除，然后新建目录
        if(NULL!=opendir(fileNameDir.c_str())){
                string cmd = "rm -rf ";
                cmd += fileNameDir.c_str();
                cout << "[cmd]" <<cmd.c_str() <<endl;
                system( cmd.c_str() );
        }
        mkdir(fileNameDir.c_str(), 0775);


	/////////save file_csr_row_ptrs_and_rows_max_and_cols_max
        ostringstream oss_1;//在string后连接int等类型
        oss_1 << fileNameDir << "//row_max_and_col_max_and_row_ptrs.txt";
        string fileNameTmp_1 = oss_1.str();

        ofstream fout_1(fileNameTmp_1.c_str());

	//save csr_rows_max
	fout_1 << csr.csr_rows_max <<endl;

 	//save csr_cols_max	
	fout_1 << csr.csr_rows_max <<endl; 

 	//save csr_row_ptrs
        for (int k= 0; k < csr.csr_row_ptrs.size(); k++){
		
		////////////////////////////////////
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		要保存原始的LongID，而非映射后的IntID
		///////////////////////////////////

		cout << csr.csr_row_ptrs.get(k) << endl;
                fout_1 << csr.csr_row_ptrs.get(k) << endl;
        }
        fout_1.close();




 	/////////save file_csr_cols_and_csr_vals
        ostringstream oss_2;//在string后连接int等类型
        oss_2 << fileNameDir << "//cols_and_vals";
        string fileNameTmp_2 = oss_2.str();

        ofstream fout_2(fileNameTmp_1.c_str());

        for (int k= 0; k < csr.csr_cols.size(); k++){
                cout << csr.csr_cols.get(k) << "," << csr.csr_vals.get(k)<<endl;
                fout_2 << csr.csr_cols.get(k) << "," << csr.csr_vals.get(k)<<endl;

        }
        fout_2.close();


	!!!!!!!!!!!!!!!!!!!!!!
	/////////////////save LongID与IntID的映射表
	!!!!!!!!!!!!!!!!!!!!!!

}
//////////////////saveCSRAsFiles end//////////////////////////////////////////////////



//////////////////CSRFileReader begin///////////////////////////////////////////
CSR CSRFileReader(string file_csr_row_ptrs_and_rows_max_and_cols_max, 
		  string file_csr_cols_and_csr_vals,		   
		  const char regex)
{
        CSR csr(-1);

        
        map<long,int> Long2IntMap;
        vector<long> Int2LongVector;




	//这里从文件中加载csr_row_ptrs和csr_rows_max和csr_cols_max
        ifstream fin_1(file_csr_cols_and_csr_vals.c_str());//在g++中使用字符串str时，要用str.c_str()
        string line_1;
        int ind = 0;
        while (getline(fin_1, line_1))
        {
                //cout << "Read from file: " << line << endl; 
                long row = atol(line_1.c_str());//string->long
                if(ind==0)
                        csr.csr_rows_max = row;
                else if(ind==1)
                        csr.csr_cols_max = row;
                else{
                        int new_row = LongID_To_IntID(Long2IntMap, Int2LongVector, row);
                        cout<< "[" << row << "][" << new_row << "]" << endl;
                        csr.csr_rows.push_back(new_row);
                }
                ind++;
        }
        fin_1.close();



	//这里从文件中加载csr_cols和csr_vals两个vector
	ifstream fin_2(file_csr_cols_and_csr_vals.c_str());//在g++中使用字符串str时，要用str.c_str()
        string line_2;
        while (getline(fin_2, line_2))
        {
                vector <string> fields = split(line_2, regex);

                if (fields.size() == 2 ) {
         		int col = atoi(fields[0].c_str());//string->int
                        float val = atof(fields[1].c_str());//string->float

                        csr.csr_cols.push_back(col);
                        csr.csr_vals.push_back(val);
                }
        }
	fin_2.close();


	/////
        int* csr_rows_arr = &csr.csr_rows[0];//vector转成array
        int* csr_cols_arr = &csr.csr_cols[0];
        float* csr_vals_arr = &csr.csr_vals[0];

        //for (int i = 0; i < coo.coo_rows.size(); i++)
        //      cout << coo_rows_arr[i] << endl;
        //for (int i = 0; i < coo.coo_cols.size(); i++)
        //      cout << coo_cols_arr[i] << endl;
        //for (int i = 0; i < coo.coo_vals.size(); i++)
        //      cout << coo_vals_arr[i] << endl; 

        cout << "csr_rows_max " << csr.csr_rows_max << endl;
        cout << "csr_cols_max " << csr.csr_cols_max << endl;

        return csr;
}


CSR useCSRFileReader(string file_csr_row_ptrs_and_rows_max_and_cols_max,
                     string file_csr_cols_and_csr_vals, 
	             const char regex)
{
        return CSRFileReader(file_csr_row_ptrs_and_rows_max_and_cols_max,
                             file_csr_cols_and_csr_vals, 
                             regex);
}


//////////////////CSRFileReader end/////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////*/






//////////////////COOFileReader begin//////////////////////////////////////////////////////////////////////////////
vector<string> &split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

vector<string> split(const string &s, char delim) {
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}


int LongID_To_IntID(	map<long,int>& Long2IntMap,
			vector<long>& Int2LongVector , 
			long oriID)
{
	int nextLoc = Long2IntMap.size();
	if(Long2IntMap.find(oriID) == Long2IntMap.end())//将以前没有出现过的long型的原始ID，加入到HashMap中
	{
		Long2IntMap[oriID] = nextLoc ;
		Int2LongVector.push_back(oriID);
	}
	
	return Long2IntMap[oriID];
}


COO COOFileReader(string filename, const char regex)
{
	COO coo(-1);

	ifstream fin(filename.c_str());//在g++中使用字符串str时，要用str.c_str()


        map<long,int> Long2IntMap;
        vector<long> Int2LongVector;

	
	string line;
	while (getline(fin, line))
	{ 
		//cout << "Read from file: " << line << endl; 
		vector <string> fields = split(line, regex);

		if (fields.size() == 3 ) {
			long row = atol(fields[0].c_str());//string->long
			int col = atoi(fields[1].c_str());//string->int
			float val = atof(fields[2].c_str());//string->float
			
			int new_row = LongID_To_IntID(Long2IntMap, Int2LongVector, row);
			//cout<< "[" << row << "][" << new_row << "][" << col << "][" << setprecision(20) << val << "]" << endl;	


			coo.coo_rows.push_back(new_row);
			coo.coo_cols.push_back(col);
			coo.coo_vals.push_back(val);

			coo.coo_rows_max = new_row > coo.coo_rows_max ? new_row : coo.coo_rows_max;
			coo.coo_cols_max = col > coo.coo_cols_max ? col : coo.coo_cols_max; 		
		}
	}

	int* coo_rows_arr = &coo.coo_rows[0];//vector转成array
	int* coo_cols_arr = &coo.coo_cols[0];
	float* coo_vals_arr = &coo.coo_vals[0];

	//for (int i = 0; i < coo.coo_rows.size(); i++)
	//	cout << coo_rows_arr[i] << endl;
	//for (int i = 0; i < coo.coo_cols.size(); i++)
	//	cout << coo_cols_arr[i] << endl;
	//for (int i = 0; i < coo.coo_vals.size(); i++)
	//	cout << coo_vals_arr[i] << endl; 
		
	cout << "coo_rows_max " << coo.coo_rows_max << endl;
	cout << "coo_cols_max " << coo.coo_cols_max << endl;	
	
	fin.close();////////////
	return coo;
}

COO useCOOFileReader(string filename, const char regex)
{
	return COOFileReader(filename, regex);	
}
//////////////////COOFileReader end//////////////////////////////////////////////////////////////////////////////


//////////////////SpasrseCOO2CSR begin//////////////////////////////////////////////////////////////////////////////
int * SpasrseCOO2CSR(
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
	
	return h_csrRowPtr;
}

CSR useSparseCOO2CSR(COO coo){

	CSR csr(-1);
	csr.csr_cols = coo.coo_cols;
	csr.csr_vals = coo.coo_vals;	
	csr.csr_rows_max = coo.coo_rows_max;//coo_rows_max为下标的最大值，若下标最大为3，则实际有4行m=coo_rows_max+1
	csr.csr_cols_max = coo.coo_cols_max;	
	
	//csr_row_ptrs的vector中共m+1=coo_rows_max+2个元素
	for(int i=0; i<=coo.coo_rows_max+1; i++)//这里先创建出一个连续存放了coo_rows_max+2个元素的vector,然后将vector转成array
		csr.csr_row_ptrs.push_back(-1);
	
	int * csrRowPtr = SpasrseCOO2CSR(&coo.coo_rows[0], csr.csr_vals.size(), csr.csr_rows_max+1, &csr.csr_row_ptrs[0]);
	//for(int i=0 ;i<=coo.coo_rows_max+1; i++ ) {
	//	cout << csr.csr_row_ptrs[i] << endl;		
	//}
	cout<<"csrRowPtr.size" << csr.csr_row_ptrs.size() <<endl;
	return csr;
}
//////////////////SpasrseCOO2CSR end//////////////////////////////////////////////////////////////////////////////



//////////////////SpasrseCSR2CSC begin//////////////////////////////////////////////////////////////////////////////
CSC SpasrseCSR2CSC( CSR csr )
{
  	int m = csr.csr_rows_max + 1;
        int n = csr.csr_cols_max + 1;
	int nnz = csr.csr_vals.size();

	//////CSC init start///////////////////////
	CSC csc(-1);
        csc.csc_rows_max = csr.csr_cols_max;//coo_rows_max为下标的最大值，若下标最大为3，则实际有4行m=coo_rows_max+1
        csc.csc_cols_max = csr.csr_rows_max;
	 
        for(int i=0; i<= n; i++)//csc_col_ptrs的vector中共n+1个元素
                csc.csc_col_ptrs.push_back(-1);//这里先创建出一个连续存放了n+1个元素的vector,然后将vector转成array
	for(int i=0; i < nnz; i++)
		csc.csc_rows.push_back(-1);
	for(int i=0; i < nnz; i++)
                csc.csc_vals.push_back(-1);
	//////CSC init end/////////////////////////


        ////start1////
        cusparseHandle_t handle = 0;
        cusparseStatus_t cusparseStatus = cusparseCreate(&handle);
        ////end1////


        ////start2////
        /////csr
	int *h_csrRowPtrs = &csr.csr_row_ptrs[0];
	int *d_csrRowPtrs;
        cudaMalloc((void**)&d_csrRowPtrs, sizeof(int)*(m+1));
        cudaMemcpy(d_csrRowPtrs, h_csrRowPtrs, (m+1)*sizeof(int), cudaMemcpyHostToDevice);

	int *h_csrCols = &csr.csr_cols[0];
        int *d_csrCols;
        cudaMalloc((void**)&d_csrCols, sizeof(int)*(nnz));
        cudaMemcpy(d_csrCols, h_csrCols, (nnz)*sizeof(int), cudaMemcpyHostToDevice);
        
	float *h_csrVals = &csr.csr_vals[0];
	float *d_csrVals;
        cudaMalloc((void**)&d_csrVals, sizeof(float)*(nnz));
        cudaMemcpy(d_csrVals, h_csrVals, (nnz)*sizeof(float), cudaMemcpyHostToDevice);

	/////csc
	int *h_cscColPtrs = &csc.csc_col_ptrs[0];
	int *d_cscColPtrs;
        cudaMalloc((void**)&d_cscColPtrs, sizeof(int)*(n+1));
        cudaMemcpy(d_cscColPtrs, h_cscColPtrs, (n+1)*sizeof(int), cudaMemcpyHostToDevice);

	int *h_cscRows = &csc.csc_rows[0]; 
	int *d_cscRows;
        cudaMalloc((void**)&d_cscRows, sizeof(int)*(nnz));
        cudaMemcpy(d_cscRows, h_cscRows, (nnz)*sizeof(int), cudaMemcpyHostToDevice);

	float *h_cscVals = &csc.csc_vals[0];
	float *d_cscVals;
        cudaMalloc((void**)&d_cscVals, sizeof(float)*(nnz));
        cudaMemcpy(d_cscVals, h_cscVals, (nnz)*sizeof(float), cudaMemcpyHostToDevice);
	////end2////


        cusparseScsr2csc(
			handle,//handle
			m,//m
			n,//n
			nnz,//nnz

			d_csrVals,//csrval
			d_csrRowPtrs,//csrRowPtr
			d_csrCols,//csrColInd

			d_cscVals,//cscVal
			d_cscRows,//cscRowInd
			d_cscColPtrs,//cscColPtr

			CUSPARSE_ACTION_NUMERIC,//copyValues, CUSPARSE_ACTION_NUMERIC: operation performed on data and indices
 			CUSPARSE_INDEX_BASE_ZERO//idxBase
                        );             

        ////start3//// 
        cudaMemcpy(h_cscColPtrs, d_cscColPtrs, (n+1)*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_cscRows, d_cscRows, (nnz)*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_cscVals, d_cscVals, (nnz)*sizeof(float), cudaMemcpyDeviceToHost);

        cusparseDestroy(handle);

	cudaFree(d_csrRowPtrs);
        cudaFree(d_csrCols);
 	cudaFree(d_csrVals);

        cudaFree(d_cscColPtrs);        
	cudaFree(d_cscRows);
        cudaFree(d_cscVals);

        cudaDeviceReset();
        ////end3////

        return csc;
}

CSC useSparseCSR2CSC(CSR csr){

       	cout << "+++++++++++++++++++++++" << endl;
  	cout << "csr_row_ptrs.size "<< csr.csr_row_ptrs.size() << endl;
        //for(int i=0; i<csr.csr_row_ptrs.size(); i++ ) 
        //     cout << csr.csr_row_ptrs[i] << endl;       

        cout << "csr_cols.size " << csr.csr_cols.size() << endl;
        //for(int i=0; i<csr.csr_cols.size(); i++ )    
        //      cout << csr.csr_cols[i] << endl;    

        cout << "csr_vals.size " << csr.csr_vals.size() << endl;
        //for(int i=0; i<csr.csr_vals.size(); i++ )    
        //      cout << csr.csr_vals[i] << endl; 
	cout << "-----------------------" << endl;


        CSC csc = SpasrseCSR2CSC(csr);
        

	cout << "+++++++++++++++++++++++" << endl;
	cout << "csc_col_ptrs.size " << csc.csc_col_ptrs.size() << endl;
	//for(int i=0; i<csc.csc_col_ptrs.size(); i++ ) 
        //      cout << csc.csc_col_ptrs[i] << endl;          

	cout << "csc_rows.size " << csc.csc_rows.size() << endl;
        //for(int i=0; i<csc.csc_rows.size(); i++ ) 
        //      cout << csc.csc_rows[i] << endl;       

	cout << "csc_vals.size " << csc.csc_vals.size() << endl;
        //for(int i=0; i<csc.csc_vals.size(); i++ ) 
        //      cout << csc.csc_vals[i] << endl;       
	cout << "----------------------" << endl;

        return csc;
}
//////////////////SpasrseCSR2CSC end//////////////////////////////////////////////////////////////////////////////




//////////////////SparseMultiply begin//////////////////////////////////////////////////////////////////////////////
CSR SparseMultiply(int m,int n, int k,int nnzA,int nnzB,
					float *h_A,int *h_RowA,int *h_ColA,
					float *h_B,int *h_RowB,int *h_ColB)
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
	cudaMalloc((void**)&d_RowB, sizeof(int)*(n+1));
	cudaMalloc((void**)&d_ColB, sizeof(int)*(nnzB));


	cudaMemcpy(d_A, h_A, nnzA*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_RowA, h_RowA, (m+1)*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ColA, h_ColA, nnzA*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(d_B, h_B, nnzB*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_RowB, h_RowB, (n+1)*sizeof(int), cudaMemcpyHostToDevice);
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
		cout << "[null != nnzTotalDevHostPtr]" << endl;
		nnzC = *nnzTotalDevHostPtr;
	} else {
		cudaMemcpy(&nnzC, d_RowC+m, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&baseC, d_RowC, sizeof(int), cudaMemcpyDeviceToHost);

		cout << "[null == nnzTotalDevHostPtr]" << endl;
		cout << "[null == nnzTotalDevHostPtr nnzC ]" << nnzC << endl;
		cout << "[null == nnzTotalDevHostPtr baseC ]" << baseC << endl;
		nnzC -= baseC;
	}

	
	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	cout << "[nnzC First]" << nnzC << endl;

	vector<int> c_csr_row_ptrs;
	for(int i=0; i<=m; i++)//这里先创建出一个连续存放了m+1个元素的vector,然后将vector转成array
		c_csr_row_ptrs.push_back(-1);
		
	vector<int> c_cols;	
	vector<float> c_vals;
	for(int i=0; i<nnzC; i++) { //这里先创建出一个连续存放了nnzC个元素的vector,然后将vector转成array
		c_cols.push_back(-1);
		c_vals.push_back(-1);
	}		
	///host上的矩阵相乘后的结果
	float *h_C = &c_vals[0];
	int *h_RowC = &c_csr_row_ptrs[0];
	int *h_ColC = &c_cols[0];	
	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	
	
	
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
	
	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////	
	CSR csr_C(-1);
	csr_C.csr_rows_max = m-1;//csr_rows下标的最大值，若下标最大为3，则实际有4行
	csr_C.csr_cols_max = k-1;//csr_cols下标的最大值，若下标最大为3，则实际有4列
		
	cout << "[nnzC]" << nnzC << endl;	
	for(int i=0; i<=m; i++) {
	//	cout << "[h_RowC]"<< h_RowC[i] << endl; 
		csr_C.csr_row_ptrs.push_back(h_RowC[i]);	
	}
	for(int i=0; i<nnzC; i++) {
	//	cout << "[h_ColC]"<< h_ColC[i] << endl;
		csr_C.csr_cols.push_back(h_ColC[i]);	
	}
	for(int i=0; i<nnzC; i++) {
	//	cout << "[h_C]"<< h_C[i] << endl;
		csr_C.csr_vals.push_back(h_C[i]);
	}
	//////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	
	return csr_C;
}

CSR useSparseMultiply(CSR csr_A, CSR csr_B){

	float* ValA = &csr_A.csr_vals[0];
    	int* RowA = &csr_A.csr_row_ptrs[0];
    	int* ColA = &csr_A.csr_cols[0];     
	//以上代码在给矩阵A赋值，使用CSR格式
	
	float* ValB = &csr_B.csr_vals[0];
	int* RowB = &csr_B.csr_row_ptrs[0];
	int* ColB = &csr_B.csr_cols[0]; 
	//以上代码在给矩阵B赋值，使用CSR格式 
	

	int n_new = (csr_A.csr_cols_max+1) >  (csr_B.csr_rows_max+1) ?  (csr_A.csr_cols_max+1) : (csr_B.csr_rows_max+1) ;

	CSR csr_C = SparseMultiply(
		csr_A.csr_rows_max+1,//int m,
		n_new,		     //int n,
		csr_B.csr_cols_max+1,//int k,
		csr_A.csr_vals.size(),//int nnzA,
		csr_B.csr_vals.size(),//int nnzB,
		
		ValA,//float *h_A,
		RowA,//int *h_RowA,
		ColA,//int *h_ColA,
		
		ValB,//float *h_B,
		RowB,//int *h_RowB,
		ColB //int *h_ColB
	   );
	   
	   
	   	
	cout << "+++1++++" << endl;
	//for(int i=0; i<= csr_C.csr_rows_max+1; i++) //行数m = 下标最大值csr_rows_max+1
	//	cout << csr_C.csr_row_ptrs[i] << endl;
 	cout << "multiply csr_c.csr_row_ptrs.size " << csr_C.csr_row_ptrs.size() << endl;	

	cout << "+++2++++" << endl;
	//for(int i=0; i< csr_C.csr_cols.size(); i++) 
	//	cout << csr_C.csr_cols[i] << endl;
 	cout <<"mulitiply csr_c.csr_cols.size " << csr_C.csr_cols.size() << endl;	

	cout << "+++3++++" << endl;
	//for(int i=0; i< csr_C.csr_vals.size(); i++) 
	//	cout << csr_C.csr_vals[i] << endl;    
 	cout<< "multiply csr_c.csr_vals.size " << csr_C.csr_vals.size() << endl;
	return csr_C;
}
//////////////////SparseMultiply end//////////////////////////////////////////////////////////////////////////////




//////////////////SparseAddition start///////////////////////////////////////////////////////////////////////////
CSR SparseAddition(int m,int n, int nnzA,int nnzB,
                                        float *h_A,int *h_RowA,int *h_ColA,
                                        float *h_B,int *h_RowB,int *h_ColB)
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
        cudaMalloc((void**)&d_RowB, sizeof(int)*(m+1));
        cudaMalloc((void**)&d_ColB, sizeof(int)*(nnzB));


        cudaMemcpy(d_A, h_A, nnzA*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_RowA, h_RowA, (m+1)*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_ColA, h_ColA, nnzA*sizeof(int), cudaMemcpyHostToDevice);

        cudaMemcpy(d_B, h_B, nnzB*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_RowB, h_RowB, (m+1)*sizeof(int), cudaMemcpyHostToDevice);
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

        cusparseXcsrgeamNnz(
                handle,
                m,
                n,
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


        //////////////////////////////////////////////////////
        //////////////////////////////////////////////////////
        cout << "[nnzC First]" << nnzC << endl;

        vector<int> c_csr_row_ptrs;
        for(int i=0; i<=m; i++)//这里先创建出一个连续存放了m+1个元素的vector,然后将vector转成array
                c_csr_row_ptrs.push_back(-1);

        vector<int> c_cols;
        vector<float> c_vals;
        for(int i=0; i<nnzC; i++) { //这里先创建出一个连续存放了nnzC个元素的vector,然后将vector转成array
                c_cols.push_back(-1);
                c_vals.push_back(-1);
        }
        ///host上的矩阵相乘后的结果
        float *h_C = &c_vals[0];
        int *h_RowC = &c_csr_row_ptrs[0];
        int *h_ColC = &c_cols[0];
        //////////////////////////////////////////////////////
        //////////////////////////////////////////////////////



        cudaMalloc((void**)&d_ColC, sizeof(int)*nnzC);
        cudaMalloc((void**)&d_C, sizeof(float)*nnzC);


	float alpha = 0.0;
	float beta = 0.0;

        cusparseScsrgeam(
                handle,
                m,
                n,
                &alpha,////////
                descrA,
                nnzA,
                d_A,
                d_RowA,
                d_ColA,
		&beta,///////
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

        //////////////////////////////////////////////////////
        //////////////////////////////////////////////////////  
        CSR csr_C(-1);
        csr_C.csr_rows_max = m-1;//csr_rows下标的最大值，若下标最大为3，则实际有4行
        //csr_C.csr_cols_max = k-1;//csr_cols下标的最大值，若下标最大为3，则实际有4列

        cout << "[nnzC]" << nnzC << endl;
        for(int i=0; i<=m; i++) {
        //      cout << "[h_RowC]"<< h_RowC[i] << endl; 
                csr_C.csr_row_ptrs.push_back(h_RowC[i]);
        }
        for(int i=0; i<nnzC; i++) {
        //      cout << "[h_ColC]"<< h_ColC[i] << endl;
                csr_C.csr_cols.push_back(h_ColC[i]);
        }
        for(int i=0; i<nnzC; i++) {
        //      cout << "[h_C]"<< h_C[i] << endl;
                csr_C.csr_vals.push_back(h_C[i]);
        }
        //////////////////////////////////////////////////////
        //////////////////////////////////////////////////////

        return csr_C;
}


CSR useSparseAddition(CSR csr_A, CSR csr_B){


        float* ValA = &csr_A.csr_vals[0];
        int* RowA = &csr_A.csr_row_ptrs[0];
        int* ColA = &csr_A.csr_cols[0];
        //以上代码在给矩阵A赋值，使用CSR格式

        float* ValB = &csr_B.csr_vals[0];
        int* RowB = &csr_B.csr_row_ptrs[0];
        int* ColB = &csr_B.csr_cols[0];
        //以上代码在给矩阵B赋值，使用CSR格式 


        CSR csr_C = SparseAddition(
                csr_A.csr_rows_max+1,//int m,
                csr_A.csr_cols_max+1,//int n,
                csr_A.csr_vals.size(),//int nnzA,
                csr_B.csr_vals.size(),//int nnzB,

                ValA,//float *h_A,
                RowA,//int *h_RowA,
                ColA,//int *h_ColA,

                ValB,//float *h_B,
                RowB,//int *h_RowB,
                ColB //int *h_ColB
           );



        cout << "+++1++++" << endl;
        //for(int i=0; i<= csr_C.csr_rows_max+1; i++) //行数m = 下标最大值csr_rows_max+1
        //      cout << csr_C.csr_row_ptrs[i] << endl;
        cout << "addition csr_c.csr_row_ptrs.size " << csr_C.csr_row_ptrs.size() << endl;

        cout << "+++2++++" << endl;
        //for(int i=0; i< csr_C.csr_cols.size(); i++) 
        //      cout << csr_C.csr_cols[i] << endl;
        cout <<"addition csr_c.csr_cols.size " << csr_C.csr_cols.size() << endl;

        cout << "+++3++++" << endl;
        //for(int i=0; i< csr_C.csr_vals.size(); i++) 
        //      cout << csr_C.csr_vals[i] << endl;    
        cout<< "addition csr_c.csr_vals.size " << csr_C.csr_vals.size() << endl;
        return csr_C;
}
////////////////SparseAddition end/////////////////////////////////////////////////////////////////////////////////////

//////////////////MPI_GPUTaskAssigner begin//////////////////////////////////////////////////////////////////////////////
void getFileInfosInCurrentDir(const char* dir, vector< pair< string, unsigned int > >& fileInfos)
{
	/*/文件句柄  
	long  hFile = 0;
	//文件信息  
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(dir).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			if (!(fileinfo.attrib &  _A_SUBDIR))
			{
				string filePath = p.assign(dir).append("\\").append(fileinfo.name);
				unsigned long fileSize = fileinfo.size;
				fileInfos.push_back( make_pair(filePath, fileSize) );
			} 
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}*/

	cout << "CurDir"<< dir << endl;	
	
    	struct stat file_stats;
    	DIR *dirp;
    	struct dirent* dent;
    	dirp=opendir(dir); // specify directory here: "." is the "current directory"
   	do {
        	dent = readdir(dirp);
     
   		if (dent) {  
			if(strcmp(dent->d_name,".")==0 || strcmp(dent->d_name,"..")==0)
				continue;

			//extern int errno;
			errno = 0; 

			string cooFile ="/"; 
			cooFile = dir + cooFile + dent->d_name;
	
			stat(cooFile.c_str(), &file_stats);
			if(errno !=0 )
				printf("%s\n",strerror(errno));


			if (!stat(cooFile.c_str(), &file_stats))
			{
				//string fileName = dent->d_name;
				unsigned int fileSize = ( unsigned int )file_stats.st_size;
				//fileInfos.push_back( make_pair(fileName, fileSize) );
 				fileInfos.push_back( make_pair(cooFile, fileSize) );
				cout << "cooFile " << cooFile << " bytes " << fileSize << endl;				
			}
			else
			{
				printf("(stat() failed for this file)\n");
			}
        	}
    	} while (dent);
    	closedir(dirp);	
}


bool Less(const pair<string, unsigned int>& p1, const pair<string, unsigned int>& p2){
	return p1.second < p2.second;
}



void useGetFileInfosInCurrentDir(const char* dir_str, vector< pair< string, unsigned int > >& fileInfos) {
	/*char * dir = "D:\\SparkData";
	vector< pair<string, unsigned long> > fileInfos;
	getFileInfosInCurrentDir(dir, fileInfos);
	for (int i = 0; i<fileInfos.size(); i++)
	{
		cout << "filePath " << fileInfos[i].first << " fileSize " << fileInfos[i].second << endl;
	}*/
	//char* dir = new char[dir_str.length() + 1];
	getFileInfosInCurrentDir(dir_str, fileInfos);

	sort(fileInfos.begin(), fileInfos.end(), Less);
}



///////
void saveAsTaskFile(string fileNameDir, int devCount, vector< pair<string, unsigned int> > fileInfos) {
	
	//目录不为空则删除，然后新建目录
	if(NULL!=opendir(fileNameDir.c_str())){ 
		string cmd = "rm -rf ";
		cmd += fileNameDir.c_str();
		cout << "[cmd]" <<cmd.c_str() <<endl;
		system( cmd.c_str() );
	}
	mkdir(fileNameDir.c_str(), 0775);
	
	
	//int eachDevTaskSize = fileInfos.size() / devCount;

	cout <<" devCount " << devCount <<endl;
	for (int i = 0; i< devCount ; i++) {

		ostringstream oss;//在string后连接int等类型
		oss << fileNameDir << "//" << i << ".txt";
		string fileNameTmp = oss.str();		 

		ofstream fout(fileNameTmp.c_str());

		for (int j = i; j < fileInfos.size(); j = j + devCount){
			cout << fileInfos[j].first << endl;
			fout << fileInfos[j].first << endl;
		}
		fout.close();
	}
}




///////
vector<string> getLocalCOOFiles(string taskTmpDataDir, int dev){
	vector<string> localCOOFiles;

 	ostringstream oss;//在string后连接int等类型
        oss << taskTmpDataDir << "/" << dev << ".txt";
        string taskTmpDataFile = oss.str();


	ifstream fin(taskTmpDataFile.c_str());//在g++中使用字符串str时，要用str.c_str()

        string line;
        while (getline(fin, line))
        {
                cout << "Read from file: " << line << endl;
                localCOOFiles.push_back(line.c_str());
        }


	return localCOOFiles;
}




template < class T >
void ClearVector( vector< T >& vt ) 
{	////swap()是交换函数，使vector离开其自身的作用域，从而强制释放vector所占的内存空间
    	vector< T > vecTemp; 
    	vecTemp.swap( vt );
}


void MPI_GPUTaskAssigner(int argc, char **argv){
	int ierr, num_procs, my_id;


	///////////////////////////////////////////////////////////////////////////////////                             
        int devCount;
        cudaGetDeviceCount(&devCount);
        printf("There are %d CUDA devices.\n", devCount);
        ///////////////////////////////////////////////////////////////////////////////////     


	ierr = MPI_Init(&argc, &argv);

	/* find out MY process ID, and how many processes were started. */

	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	cout<< "num_procs " << num_procs <<endl;
	
	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name(processor_name, &namelen); 	
	printf("out Spawning from %s %d\n", processor_name,my_id);


	///////////////////////////////////////////////////////////////////////////////////
	int dev = (my_id-1)%devCount ;//my_id==0是主进程不干活
	cudaSetDevice(dev);	
	//////////////////////////////////////////////////////////////
        string taskTmpDataDir = "/gruntdata/app_data/zhangmeng.zm/TuTongKuang_01265/TmpDataNew";	
	int tag = 0;
	int number;
 	
	if(my_id == 0) {	
	
		//把处理的文件，按大小排序后发送给2个GPU进行均衡分配
		//MPI_Send(&fileNames[0], fileNames.size(), MPI_INT, 1, 0 , MPI_COMM_WORLD);
	

 		MPI_Get_processor_name(processor_name, &namelen);
    		printf("my_id == 0  Spawning from %s \n", processor_name);


		vector< pair<string, unsigned int> > fileInfos;
		useGetFileInfosInCurrentDir("/gruntdata/app_data/zhangmeng.zm/TuTongKuang_01265/COOFilesNew", fileInfos);
	
		for (int i = 0; i<fileInfos.size(); i++)
            		cout << "sortedFilePath " << fileInfos[i].first << " fileSize " << fileInfos[i].second << endl;
        	

		saveAsTaskFile(taskTmpDataDir.c_str(), devCount, fileInfos); 

		////////
        	//把处理的文件，按大小排序后发送给2个GPU进行均衡分配
		number = -1 ; 
		for(int i=1; i<num_procs; i++)
        		MPI_Send(&number, 1, MPI_INT, i, tag , MPI_COMM_WORLD);
	
	} else
	{
		////////
		MPI_Recv(&number, 1, MPI_INT, 0, tag , MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Get_processor_name(processor_name, &namelen);
                printf("my_id == %d  Spawning from %s \n", my_id, processor_name);

		///////
		vector<string> localCOOFiles = getLocalCOOFiles(taskTmpDataDir.c_str() , dev);
		printf("my_id == %d  localCOOFiles.size() : %d\n", my_id, localCOOFiles.size() );
		for(int i=0; i<localCOOFiles.size(); i++) {

			///////////////////////////////////////////////////////////////////////////////////		
			cout << "### " << localCOOFiles[i].c_str() << endl;
			clock_t t1 = clock();

			COO coo = useCOOFileReader(localCOOFiles[i].c_str(), ','); //COO文件中分隔符为英文逗号
			clock_t t2 = clock();
			cout << "[useCOOFileReader cost: " << (t2-t1)/CLOCKS_PER_SEC << " seconds]" << endl;
				

			CSR csr = useSparseCOO2CSR(coo);
                        clock_t t3 = clock();
                        cout << "[useSparseCOO2CSR cost: " << (t3-t2)/CLOCKS_PER_SEC << " seconds]" << endl;


			CSC csc = useSparseCSR2CSC(csr);
                        clock_t t4 = clock();
                        cout << "[useSparseCSR2CSC cost: " << (t4-t3)/CLOCKS_PER_SEC << " seconds]" << endl;


			CSR csr_C = useSparseMultiply(csr, csc.use_CSC_Create_CSR());
			//CSR csr_C = useSparseAddition(csr, csr);
			clock_t t5 = clock();
                        cout << "[useSparseMultiply cost: " << (t5-t4)/CLOCKS_PER_SEC << " seconds]" << endl;
			//cout << "[useSparseAddition cost: " << (t5-t4)/CLOCKS_PER_SEC << " seconds]" << endl;
			///////////////////////////////////////////////////////////////////////////////////	
		}

		printf("Hello world! I'm process %i out of %i processes , I am on Dev %i \n", 
				my_id, num_procs, dev);		
	}
	
	ierr = MPI_Finalize();
}
//////////////////MPI_GPUTaskAssigner end//////////////////////////////////////////////////////////////////////////////

//////////////////div_by_row_and_multiply begin////////////////////////////////////////////////////////////////////////
void MPI_GPUTaskAssigner_New(int argc, char **argv){
        int ierr, num_procs, my_id;


        ///////////////////////////////////////////////////////////////////////////////////                             
        int devCount;
        cudaGetDeviceCount(&devCount);
        printf("There are %d CUDA devices.\n", devCount);
        ///////////////////////////////////////////////////////////////////////////////////     


        ierr = MPI_Init(&argc, &argv);

        /* find out MY process ID, and how many processes were started. */

        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
        ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        cout<< "num_procs " << num_procs <<endl;

        int namelen;
        char processor_name[MPI_MAX_PROCESSOR_NAME];
        MPI_Get_processor_name(processor_name, &namelen);
        printf("out Spawning from %s %d\n", processor_name,my_id);


        ///////////////////////////////////////////////////////////////////////////////////
        int dev = (my_id-1)%devCount ;//my_id==0是主进程不干活
        cudaSetDevice(dev);
        //////////////////////////////////////////////////////////////
        string taskTmpDataDir = "/gruntdata/app_data/zhangmeng.zm/TuTongKuang_01265/TmpDataNew";
        int tag = 0;
        int number;



	////所有的COO文件及其大小
 	vector< pair<string, unsigned int> > fileInfos;
        useGetFileInfosInCurrentDir("/gruntdata/app_data/zhangmeng.zm/TuTongKuang_01265/COOFilesNew", fileInfos);




        if(my_id == 0) {
                //把处理的文件，按大小排序后，进行均衡分配并写入两个info文件中，然后通知2个GPU去读各自的info文件
                MPI_Get_processor_name(processor_name, &namelen);
                printf("my_id == 0  Spawning from %s \n", processor_name);



                for (int i = 0; i<fileInfos.size(); i++)
                        cout << "sortedFilePath " << fileInfos[i].first << " fileSize " << fileInfos[i].second << endl;


                saveAsTaskFile(taskTmpDataDir.c_str(), devCount, fileInfos);

                ////////
                //把处理的文件，按大小排序后发送给2个GPU进行均衡分配
                number = -1 ;
                for(int i=1; i<num_procs; i++)
                        MPI_Send(&number, 1, MPI_INT, i, tag , MPI_COMM_WORLD);/////这里只是发一个同步操作，并不传实际数据

        } else {
                /////这里只是进行一个同步操作，具体数据在my_id==0的机器的文件系统中
                MPI_Recv(&number, 1, MPI_INT, 0, tag , MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                MPI_Get_processor_name(processor_name, &namelen);
                printf("my_id == %d  Spawning from %s \n", my_id, processor_name);

                ///////vector中元素数量不多时，没必要使用set、hash_set
        	vector<string> localCOOFiles = getLocalCOOFiles(taskTmpDataDir.c_str() , dev); 

		printf("my_id == %d  localCOOFiles.size() : %d\n", my_id, localCOOFiles.size() );
                printf("my_id == %d  fileInfos.size() : %d\n", my_id, fileInfos.size() );

		vector<CSC> csc_all ;
		vector<CSR> csr_part ;
		for(int i=0; i < fileInfos.size(); i++) {
		
                        ///////////////////////////////////////////////////////////////////////////////////             
                        string tmp_all_file = fileInfos[i].first.c_str();
			cout << "[### "<< i << " ] " << tmp_all_file << endl;
                        clock_t t1 = clock();

                        COO coo = useCOOFileReader(tmp_all_file, ','); //COO文件中分隔符为英文逗号
                        clock_t t2 = clock();
                        cout << "[useCOOFileReader cost: " << (t2-t1)/CLOCKS_PER_SEC << " seconds]" << endl;
                                

                        CSR csr = useSparseCOO2CSR(coo);
 			ClearVector(coo.coo_rows);//将COO对象占用的内存空间释放掉
                        ClearVector(coo.coo_cols);
                        ClearVector(coo.coo_vals);
                        clock_t t3 = clock();
			cout << "[useSparseCOO2CSR cost: " << (t3-t2)/CLOCKS_PER_SEC << " seconds]" << endl;


                        if(std::find(localCOOFiles.begin(), localCOOFiles.end(), tmp_all_file) != localCOOFiles.end())
			{	//只将本组需要的csr载入内存
				csr_part.push_back(csr);
				cout << "[csr_part.size]" << csr_part.size() << "[csr_part add]" << tmp_all_file  << endl;
			}

			
                        CSC csc = useSparseCSR2CSC(csr);
                        ClearVector(csr.csr_row_ptrs);//将CSR对象占用的内存空间释放掉
                        ClearVector(csr.csr_cols);
                        ClearVector(csr.csr_vals);
			clock_t t4 = clock();
 			cout << "[useSparseCSR2CSC cost: " << (t4-t3)/CLOCKS_PER_SEC << " seconds]" << endl;

			csc_all.push_back(csc);//需要将所有的csc载入内存
			cout << "[csc_all.size]" << csc_all.size() << endl;
                        ///////////////////////////////////////////////////////////////////////////////////     
                }

		for(int k=0; k<csc_all.size(); k++)
		{
			cout << "[csc_all " << k << "][csc_col_prts.size]" << csc_all[k].csc_col_ptrs.size() << endl;
			cout << "[csc_all " << k << "][csc_rows.size]" << csc_all[k].csc_rows.size() << endl;
			cout << "[csc_all " << k << "][csc_vals.size]" << csc_all[k].csc_vals.size() << endl;
			
			CSR csc_all_csr =  csc_all[k].use_CSC_Create_CSR();
			cout << "[csc_all_csr " << k << "][csr_row_prts.size]" << csc_all_csr.csr_row_ptrs.size() << endl;
			cout << "[csc_all_csr " << k << "][csr_cols.size]" << csc_all_csr.csr_cols.size() << endl;
			cout << "[csc_all_csr " << k << "][csr_vals.size]" << csc_all_csr.csr_vals.size() << endl;
		}

                //////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////
		for(int i=0; i< csr_part.size(); i++){
			for(int j=0; j<csc_all.size(); j++){
		
				///////multiply csr_a csc_b
				clock_t t6 = clock();
                		CSR csr_c = useSparseMultiply(csr_part[i], csc_all[j].use_CSC_Create_CSR());
				clock_t t7 = clock();
				cout << "[csr " << i <<"][left_cols_max"<<csr_part[i].csr_cols_max << "]"<< endl;
				cout << "[csc " << j <<"][right_rows_max"<<csc_all[j].use_CSC_Create_CSR().csr_rows_max<<"]"<<endl;
                		cout << "[useSparseMultiply cost: " << (t7-t6)/CLOCKS_PER_SEC << " seconds]" << endl;
				cout << "==================================================================" << endl;
					
				ClearVector(csr_c.csr_row_ptrs);//将CSR对象占用的内存空间释放掉
                        	ClearVector(csr_c.csr_cols);
                        	ClearVector(csr_c.csr_vals);
			}
			
			ClearVector(csr_part[i].csr_row_ptrs);//将CSR对象占用的内存空间释放掉
                        ClearVector(csr_part[i].csr_cols);
                        ClearVector(csr_part[i].csr_vals);
		}






		/*
                string file_a = "/gruntdata/app_data/zhangmeng.zm/TuTongKuang_01265/COOFilesNew/div_by_row_suffix_mod_11_3";
                string file_b = "/gruntdata/app_data/zhangmeng.zm/TuTongKuang_01265/COOFilesNew/div_by_row_suffix_mod_11_9";

                cout << "### file_a " << file_a << endl;
                cout << "### file_b " << file_b << endl;
                clock_t t1 = clock();


                //////left csr_a
                COO coo_a = useCOOFileReader(file_a, ','); //COO文件中分隔符为英文逗号
                clock_t t2 = clock();
                cout << "[useCOOFileReader coo_a cost: " << (t2-t1)/CLOCKS_PER_SEC << " seconds]" << endl;


                CSR csr_a = useSparseCOO2CSR(coo_a);
                clock_t t3 = clock();
                cout << "[useSparseCOO2CSR csr_a cost: " << (t3-t2)/CLOCKS_PER_SEC << " seconds]" << endl;
                ClearVector(coo_a.coo_rows);//将COO对象占用的内存空间释放掉
                ClearVector(coo_a.coo_cols);
                ClearVector(coo_a.coo_vals);


                ///////right csc_b
                COO coo_b = useCOOFileReader(file_b, ','); //COO文件中分隔符为英文逗号
                clock_t t4 = clock();
                cout << "[useCOOFileReader coo_b cost: " << (t4-t3)/CLOCKS_PER_SEC << " seconds]" << endl;


                CSR csr_b = useSparseCOO2CSR(coo_b);
                clock_t t5 = clock();
                cout << "[useSparseCOO2CSR csr_b cost: " << (t5-t4)/CLOCKS_PER_SEC << " seconds]" << endl;
                ClearVector(coo_b.coo_rows);//将COO对象占用的内存空间释放掉
                ClearVector(coo_b.coo_cols);
                ClearVector(coo_b.coo_vals);


                CSC csc_b = useSparseCSR2CSC(csr_b);
                clock_t t6 = clock();
                cout << "[useSparseCSR2CSC csc_b cost: " << (t6-t5)/CLOCKS_PER_SEC << " seconds]" << endl;


                ///////multiply csr_a csc_b
                CSR csr_C = useSparseMultiply(csr_a, csc_b.use_CSC_Create_CSR());
                clock_t t7 = clock();
                cout << "[useSparseMultiply cost: " << (t7-t6)/CLOCKS_PER_SEC << " seconds]" << endl;*/
                //////////////////////////////////////////////////////////////////////////////////////
                //////////////////////////////////////////////////////////////////////////////////////


                printf("Hello world! I'm process %i out of %i processes , I am on Dev %i \n",
                                my_id, num_procs, dev);
        }

        ierr = MPI_Finalize();
}
//////////////////div_by_row_and_multiply end////////////////////////////////////////////////////////////////////////



int main(int argc, char **argv)
{
	MPI_GPUTaskAssigner_New(argc, argv);
	//MPI_GPUTaskAssigner(argc, argv); 
	return 0;
} 
