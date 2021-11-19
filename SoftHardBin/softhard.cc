#include "CompareMBandHB.h"
#include "aux_func.h"


int main(int argc, char* argv[]){
	int argnum = argc;
	vector<string> args(argv, argv + argc);
	hard(argnum, args);
	soft(argnum, args);
	CompareMBandHB(args);
	return 0;
}