#pragma once
#include <sys/stat.h>
#include <string>
#include <ctime>
#include <vector>

#include "TMath.h"
#include "RtypesCore.h"
using namespace std;

const double pi = TMath::Pi();

Double_t RecalcAngle(const double& _angle);
void MakeDir(const string& dirpath, const string& dirname);
int soft(int argc, vector<string> argv);
int hard(int argc, vector<string> argv);