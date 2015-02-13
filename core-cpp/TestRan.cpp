#include "Random.h"
#include <iostream>
#include "nvector.h"
#include "Statis.h"
#include <iomanip>
#include <stdlib.h>

using namespace std;
int main(int argc, char *argv[]) {
	clock_t __start = clock();
	char name[256];
	double a = 256568985.564;
	double b = 23330223223897546302;
	sprintf(name, "gcewfp%.15g_%.15lg_%.18gConf.txt", a, a, b);
	cout << name << endl;
	return 1;
}
