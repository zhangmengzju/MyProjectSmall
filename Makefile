NVCC = nvcc
MPI_PATH = /gruntdata/app_data/zhisong.fzs/openmpi-1.8.5rc3-build
CAFFE_PATH = /gruntdata/app_data/zhisong.fzs/caffe
CAFFE_BUILD = /gruntdata/app_data/zhisong.fzs/caffe/build
NVCC_OPTS = -O3 -lineinfo -Xcompiler -fopenmp -Xcompiler -g -Xcompiler -rdynamic\
	-I. -I/usr/local/include -I$(MPI_PATH)/include -I/usr/include \

LD_OPTIONS = -L/usr/lib/ -L/usr/local/lib -L/gruntdata/app_data/zhisong.fzs/caffe/build/lib -L/grunt/anaconda/pkgs/hdf5-1.8.9-1/lib \
	-L/grunt/anaconda/pkgs/opencv-2.4.6-np18py27_0/lib \
	-L/grunt/anaconda/pkgs/python-2.7.6-1/lib \
	-L/opt/intel/mkl/lib/intel64 \
	-L/grunt/anaconda/pkgs/jpeg-8d-0/lib \
	-L/grunt/anaconda/pkgs/libpng-1.5.13-1/lib \
	-L/gruntdata/app_data/zhisong.fzs/tfs/restful/lib \
	-L../preprocess \
	-L.

LD_LIBS = -lcusparse -lmpi -lstdc++ -lgomp \

GEN_SM35 = -gencode=arch=compute_35,code=\"sm_35,compute_35\" 
GEN_SM30 = -gencode=arch=compute_30,code=\"sm_30,compute_30\" 
GEN_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\" 
GEN_SM13 = -gencode=arch=compute_13,code=\"sm_13,compute_13\" 
GEN_SM10 = -gencode=arch=compute_10,code=\"sm_10,compute_10\" 
SM_TARGETS = $(GEN_SM20) $(GEN_SM30) $(GEN_SM35) 

all: Style 

MultiSparseMultiply.o : MultiSparseMultiply.cu Makefile
	$(NVCC) --compiler-bindir $(MPI_PATH)/bin/mpicc -c -o $@ $< $(NVCC_OPTS) $(SM_TARGETS)

Style: MultiSparseMultiply.o
	$(NVCC) --compiler-bindir $(MPI_PATH)/bin/mpicc -o $@ $^ $(LD_OPTIONS) $(LD_LIBS)

clean:
	rm -f Style *.o
