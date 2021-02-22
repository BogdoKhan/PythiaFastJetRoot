#include "aux_func.h"

Double_t RecalcAngle(const double& _angle){
	Double_t angle = (Double_t) _angle;
	if (angle >= 0 && angle < 2*pi){
		return angle;
	}
	while (angle >= 2*pi) {
		angle -= 2*pi;
	}
	while (angle < 0) {
		angle += 2*pi;
	}
	return angle;
}

void MakeDir(const string& dirpath, const string& dirname){
	string path = dirpath + dirname;
	const char* folder = path.c_str();
	//folder = "./Results";
	struct stat sb;

	if (stat(folder, &sb) != 0 && !S_ISDIR(sb.st_mode)) {
		mkdir(folder, 0755);
	} 
}