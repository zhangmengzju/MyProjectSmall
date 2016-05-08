#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string> 
#include <stdlib.h>//使能在g++中编译通过atoi

//using namespace std;//g++中要对string、vector、ifstream、stringstream、endl、cout等显式添加std:: 

std::vector<int> coo_rows;//避免像使用数组时，需要先确定大小
std::vector<int> coo_cols;
std::vector<int> coo_vals;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}


int main()
{
	std::string filename = "myCOO.txt";
	std::ifstream fin(filename.c_str());//在g++中使用字符串str时，要用str.c_str()
	
	char regex = ' ';//COO文件中分隔符

	std::string line;
	while (getline(fin, line))
	{ 
		std::cout << "Read from file: " << line << std::endl; 
		std::vector <std::string> fields = split(line, regex);

		if (fields.size() == 3){
			coo_rows.push_back(atoi(fields[0].c_str()));
			coo_cols.push_back(atoi(fields[1].c_str()));
			coo_vals.push_back(atoi(fields[2].c_str()));
		}
	}


	int* coo_rows_arr = &coo_rows[0];//vector转成array
	int* coo_cols_arr = &coo_cols[0];
	int* coo_vals_arr = &coo_vals[0];


	for (int i = 0; i < coo_rows.size(); i++)
		std::cout << coo_rows_arr[i] << std::endl;
	for (int i = 0; i < coo_cols.size(); i++)
		std::cout << coo_cols_arr[i] << std::endl;
	for (int i = 0; i < coo_vals.size(); i++)
		std::cout << coo_vals_arr[i] << std::endl; 
	 
	//system("pause");
	return 0;
}
