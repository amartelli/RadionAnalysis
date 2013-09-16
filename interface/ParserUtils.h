#ifndef ParserUtils_h
#define ParserUtils_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#include <algorithm>
#include <iomanip>

/*#include "TFile.h"
#include "TH1F.h"
#include "TChain.h"
#include "TVector3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Boost.h"*/
#include "ConfigParser.h"

/** get the parameters from a congiguration file */
int parseConfigFile (const TString& config) ;

#endif
