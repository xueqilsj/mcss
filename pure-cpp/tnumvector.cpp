/*
 ============================================================================
 Name        : NVTPotts.cpp
 Description : Parallel tempering in MPI C++
 ============================================================================
 */

#include "nvector.h"
#include <iostream>

int main() {
	int myints[] = { 16, 2, 77, 220 };
	/*
	 NV_FromArray(a,myints,int);
	 myints[2] = 12;
	 NV_FromArray(a2,myints,int);
	 NV_DumpFile(a, a2, "a2.txt");
	 std::vector<int> b;
	 NV_LoadFileCol("TEST.DAT", 1, b);
	 //NV_LoadFile("TEST.DAT", b);
	 NV_Print<int> (b, ",");
	 std::cout << NV_Mean<int> (b) << '\t' << NV_StdDev<int> (b) << std::endl;
	 NV_ToArray(b,c,int);
	 */
	nvector<int> a(myints, 4);
	nvector<char> d(2);
	a.dump("a.txt");
	d.load("WLPotts12x12Q10_Conf.txt", 0, true, ';',";");
	std::cout << d.len << d;
	nvector<double> vLnGe;
	nvector<unsigned long> vGeN;
	vLnGe.load("WLPotts12x12Q10_LnGe_GeN.txt", 0);
	vGeN.load("WLPotts12x12Q10_LnGe_GeN.txt", 1);

	std::cout << vLnGe.len<<'\t'<<vGeN.len;

printf("unsigned long%d,int%d,unsigned long long%d,uint32_t%d,uint64_t%d",sizeof(unsigned long),sizeof(int),sizeof(unsigned long long),sizeof(uint32_t),sizeof(uint64_t));
	return 0;

}

