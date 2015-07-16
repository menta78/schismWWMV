/* ------------------------------ */
#include <ctype.h>
#include <malloc.h>
#include <unistd.h>
#include <getopt.h>
#include <chrono>
#include <ctime>
/* ------------------------------ */


/* ------------------------------ */
#include <math.h>
/* ------------------------------ */


/* ------------------------------ */
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
/* ------------------------------ */


/* ------------------------------ */
#include <exception>
#include <vector>
#include <list>
#include <set>
#include <map>
/* ------------------------------ */


/* ------------------------------ */
#include <fstream>
#include <iostream>
#include <sstream>
/* ------------------------------ */


/* ------------------------------ */
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include "gmpxx.h"
/* ------------------------------ */


/* ------------------------------ */
#include "tbb/tbb.h"
/* ------------------------------ */


/* ------------------------------ */
typedef boost::dynamic_bitset<> Face;
/* ------------------------------ */


/* ------------------------------ */
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
/* ------------------------------ */


/* ------------------------------ */
//{
//  boost::filesystem::create_directories(eDir.c_str() );
/* ------------------------------ */


/* ------------------------------ */
#include <netcdf>
/* ------------------------------ */


/* ------------------------------ */
#include "grib_api.h"
/* ------------------------------ */


/* ------------------------------ */
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
template <typename T> using MyVector=Eigen::Matrix<T,Eigen::Dynamic,1>;
template <typename T> using MyMatrix=Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template <typename T> using MySparseMatrix=Eigen::SparseMatrix<T,Eigen::ColMajor>;
/* ------------------------------ */


/* ------------------------------ */
#include <thread>
#include <mutex>
#include <condition_variable>
/* ------------------------------ */


/* ------------------------------ */
void PrintMMA_FCT(MyMatrix<double> const& F, MyMatrix<int> const& MSK, std::string const& VarName)
{
  int eta=F.rows();
  int xi=F.cols();
  double minval, maxval;
  double sum=0;
  int nb=0;
  bool IsFirst=true;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++)
      if (MSK(i,j) == 1) {
	double eVal=F(i,j);
	sum+=eVal;
	nb++;
	if (IsFirst) {
	  IsFirst=false;
	  maxval=eVal;
	  minval=eVal;
	}
	else {
	  if (eVal > maxval)
	    maxval=eVal;
	  if (eVal < minval)
	    minval=eVal;
	}
      }
  double eMean=sum/double(nb);
  std::cerr << "  " << VarName << " min=" << minval << " max=" << maxval << " avg=" << eMean << "\n";
}
/* ------------------------------ */


/* ------------------------------ */
struct SingleBlock {
  std::map<std::string, int> ListIntValues;
  std::map<std::string, bool> ListBoolValues;
  std::map<std::string, double> ListDoubleValues;
  std::map<std::string, std::vector<double> > ListListDoubleValues;
  std::map<std::string, std::vector<int> > ListListIntValues;
  std::map<std::string, std::string> ListStringValues;
  std::map<std::string, std::vector<std::string> > ListListStringValues;
  std::string BlockName;
};
/* ------------------------------ */


/* ------------------------------ */
struct FullNamelist {
  std::map<std::string,SingleBlock> ListBlock;
  std::string FileName;
};
/* ------------------------------ */


/* ------------------------------ */
std::string NAMELIST_FindPositionVariableInBlock(std::string const& eVarName,
						 SingleBlock &eSingleBlock)
{
  auto search1=eSingleBlock.ListIntValues.find(eVarName);
  if (search1 != eSingleBlock.ListIntValues.end())
    return "int";
  auto search2=eSingleBlock.ListBoolValues.find(eVarName);
  if (search2 != eSingleBlock.ListBoolValues.end())
    return "bool";
  auto search3=eSingleBlock.ListDoubleValues.find(eVarName);
  if (search3 != eSingleBlock.ListDoubleValues.end())
    return "double";
  auto search4=eSingleBlock.ListListDoubleValues.find(eVarName);
  if (search4 != eSingleBlock.ListListDoubleValues.end())
    return "listdouble";
  auto search4b=eSingleBlock.ListListIntValues.find(eVarName);
  if (search4b != eSingleBlock.ListListIntValues.end())
    return "listint";
  auto search5=eSingleBlock.ListStringValues.find(eVarName);
  if (search5 != eSingleBlock.ListStringValues.end())
    return "string";
  auto search6=eSingleBlock.ListListStringValues.find(eVarName);
  if (search6 != eSingleBlock.ListListStringValues.end())
    return "liststring";
  return "not found";
}
/* ------------------------------ */


/* ------------------------------ */
std::string NAMELIST_RemoveAfterCommentChar(std::string const&eStr, std::string &eChar)
{
  int WeFound=0;
  std::string RetStr;
  int len=eStr.size();
  for (int i=0; i<len; i++) {
    std::string fChar=eStr.substr(i,1);
    if (fChar == eChar) {
      WeFound=1;
    }
    if (WeFound == 0) {
      RetStr=RetStr + eStr.at(i);
    }
  }
  return RetStr;
}
/* ------------------------------ */


/* ------------------------------ */
std::string NAMELIST_RemoveAfterLastChar(std::string const& eStr, std::string const& eLastChar)
{
  int iPos=-1;
  int len=eStr.size();
  for (int i=0; i<len; i++) {
    int j=len-1-i;
    if (iPos == -1) {
      std::string eChar=eStr.substr(j,1);
      //      std::cerr << "j=" << j << " eChar=" << eChar << "\n";
      if (eChar == eLastChar)
	iPos=j;
    }
  }
  if (iPos == -1)
    return eStr;
  //  std::cerr << "iPos=" << iPos << "\n";
  std::string RetStr;
  for (int i=0; i<iPos; i++)
    RetStr=RetStr + eStr.at(i);
  return RetStr;
}
/* ------------------------------ */


/* ------------------------------ */
bool NAMELIST_ReadBoolValue(std::string const& eVarValue)
{
  if (eVarValue == ".F.")
    return false;
  if (eVarValue == ".T.")
    return true;
  if (eVarValue == "F")
    return false;
  if (eVarValue == "T")
    return true;
  std::cerr << "Boolean value has not been found\n";
  std::cerr << "eVarValue = " << eVarValue << "\n";
  std::cerr << "Allowed: T / F / .T. / .F.\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
std::string NAMELIST_ConvertFortranStringToCppString(std::string const& eStr)
{
  int len=eStr.length();
  std::string eFirstChar=eStr.substr(0, 1);
  std::string eLastChar=eStr.substr(len-1, 1);
  int RemovableEnding=0;
  if (eFirstChar == "'") {
    RemovableEnding=1;
    if (eFirstChar != eLastChar) {
      std::cerr << "Inconsistency in the string\n";
      exit(1);
    }
  }
  if (eFirstChar == "\"") {
    RemovableEnding=1;
    if (eFirstChar != eLastChar) {
      std::cerr << "Inconsistency in the string\n";
      exit(1);
    }
  }
  if (RemovableEnding == 1) {
    std::string RetStr=eStr.substr(1,len-2);
    return RetStr;
  }
  return eStr;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<std::string> NAMELIST_ConvertFortranListStringToCppListString(std::string const& eStr)
{
  int len=eStr.length();
  std::string eFirstChar=eStr.substr(0, 1);
  std::string eLastChar=eStr.substr(len-1,1);
  if (eFirstChar != "'" && eFirstChar != "\"") {
    std::cerr << "For list of strings, one should use string \"  \"   or '    '   \n";
    exit(1);
  }
  if (eLastChar != "'" && eLastChar != "\"") {
    std::cerr << "For list of strings, one should use string \"  \"   or '    '   \n";
    exit(1);
  }
  if (eFirstChar != eLastChar) {
    std::cerr << "No coherency in endings\n";
    exit(1);
  }
  std::string eSepChar=eFirstChar;
  int IsInString=0;
  std::string eFound="";
  std::vector<std::string> eListStr;
  for (int i=0; i<len; i++) {
    std::string eChar=eStr.substr(i,1);
    if (eChar == eSepChar) {
      eFound=eFound + eChar;
      if (IsInString == 1) {
	IsInString=0;
	std::string eCppStr=NAMELIST_ConvertFortranStringToCppString(eFound);
	eListStr.push_back(eCppStr);
	eFound="";
      }
      else {
	IsInString=1;
      }
    }
    else {
      if (IsInString == 1)
	eFound=eFound + eChar;
    }
  }
  return eListStr;
}
/* ------------------------------ */


/* ------------------------------ */
void NAMELIST_WriteNamelistFile(std::ostream &os, FullNamelist const& eFullNamelist)
{
  int iBlock=0;
  for (std::map<std::string,SingleBlock>::const_iterator itB=eFullNamelist.ListBlock.begin(); itB!=eFullNamelist.ListBlock.end(); ++itB) {
    std::string eBlockName=itB->first;
    SingleBlock eBlock=itB->second;
    if (iBlock > 0) {
      os << "\n";
      os << "\n";
    }
    os << "&" << eBlockName << "\n";
    for (std::map<std::string,int>::const_iterator it=eBlock.ListIntValues.begin(); it!=eBlock.ListIntValues.end(); ++it)
      os << "  " << it->first << " = " << it->second << "\n";
    for (std::map<std::string,bool>::const_iterator it=eBlock.ListBoolValues.begin(); it!=eBlock.ListBoolValues.end(); ++it)
      os << "  " << it->first << " = " << it->second << "\n";
    for (std::map<std::string,double>::const_iterator it=eBlock.ListDoubleValues.begin(); it!=eBlock.ListDoubleValues.end(); ++it)
      os << "  " << it->first << " = " << it->second << "\n";
    for (std::map<std::string,std::vector<double> >::const_iterator it=eBlock.ListListDoubleValues.begin(); it!=eBlock.ListListDoubleValues.end(); ++it) {
      os << "  " << it->first << " = ";
      std::vector<double> eListDoubl=it->second;
      int nbDoubl=eListDoubl.size();
      for (int iDoubl=0; iDoubl<nbDoubl; iDoubl++) {
	if (iDoubl > 0)
	  os << ",";
	os << eListDoubl[iDoubl];
      }
      os << "\n";
    }
    for (std::map<std::string,std::vector<int> >::const_iterator it=eBlock.ListListIntValues.begin(); it!=eBlock.ListListIntValues.end(); ++it) {
      os << "  " << it->first << " = ";
      std::vector<int> eListInt=it->second;
      int nbInt=eListInt.size();
      for (int iInt=0; iInt<nbInt; iInt++) {
	if (iInt > 0)
	  os << ",";
	os << eListInt[iInt];
      }
      os << "\n";
    }
    for (std::map<std::string,std::string>::const_iterator it=eBlock.ListStringValues.begin(); it!=eBlock.ListStringValues.end(); ++it)
      os << "  " << it->first << " = " << it->second << "\n";
    for (std::map<std::string,std::vector<std::string> >::const_iterator it=eBlock.ListListStringValues.begin(); it!=eBlock.ListListStringValues.end(); ++it) {
      os << "  " << it->first << " = ";
      std::vector<std::string> eListStr=it->second;
      int nbString=eListStr.size();
      for (int iString=0; iString<nbString; iString++) {
	if (iString > 0)
	  os << ",";
	os << eListStr[iString];
      }
      os << "\n";
    }
    os << "&END\n";
    iBlock++;
  }
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<std::string> GetAllNamesOfSatelliteAltimeter()
{
  return {"ERS1", "ERS2", "ENVISAT", "TOPEX", "POSEIDON", "JASON1", "GFO", "JASON2", "CRYOSAT", "SARAL"};
}
/* ------------------------------ */


/* ------------------------------ */
void CHECK_Model_Allowedness(std::string const& eModelName)
{
  std::vector<std::string> vec={"COSMO", "WAM", "ROMS", "WWM", "WW3", "GRIB_DWD", "GRIB_ECMWF", "GRIB_GFS", "GRIB_COSMO"};
  bool isPresent = (std::find(vec.begin(), vec.end(), eModelName) != vec.end());
  if (isPresent == false) {
    std::cerr << "We did not find the MODEL NAME\n";
    std::cerr << "MODELNAME = " << eModelName << "\n";
    std::cerr << "List of allowed models\n";
    for (int iModel=0; iModel<int(vec.size()); iModel++) {
      std::cerr << "iModel=" << iModel << " eModel=" << vec[iModel] << "\n";
    }
    exit(1);
  }
}
/* ------------------------------ */


/* ------------------------------ */
double TheSignFct(double const& eVal)
{
  if (eVal > 0)
    return 1;
  if (eVal < 0)
    return -1;
  return 0;
}
/* ------------------------------ */


/* ------------------------------ */
void DifferenceLonRenormalize(double & Lon)
{
  if (Lon > 180)
    Lon=Lon - 360;
  if (Lon < -180)
    Lon=Lon + 360;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> COMPUTE_NORM(MyMatrix<double> const& U, MyMatrix<double> const& V)
{
  int eta=U.rows();
  int xi=U.cols();
  MyMatrix<double> WINDMAG(eta, xi);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double eU=U(i,j);
      double eV=V(i,j);
      WINDMAG(i,j)=sqrt(eU*eU + eV*eV);
    }
  return WINDMAG;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> FreqPeriodChange(MyMatrix<double> const& F)
{
  int eta=F.rows();
  int xi=F.cols();
  MyMatrix<double> Fret(eta, xi);
  double pi=3.1415926535;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double eVal=F(i,j);
      double NewVal=2*pi/eVal;
      Fret(i,j)=NewVal;
    }
  return Fret;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<std::string> GetAllPossibleVariables()
{
  std::vector<std::string> ListVarOut={"WIND10", "SurfCurr", "Hwave", "WINDMAG", "TempSurf", "SaltSurf", "AIRT2", "Rh2", "ZetaOcean", "MwaveFreq", "PwaveFreq", "AIRD", "CdWave", "AlphaWave", "rain", "swrad", "lwrad", "latent", "sensible", "shflux", "ssflux", "evaporation", "MwavePer", "PwavePer", "SurfPres", "SstOcean"};
  return ListVarOut;
}
/* ------------------------------ */


/* ------------------------------ */
std::string NCL_bool(bool eBool)
{
  if (eBool)
    return "True";
  return "False";
}
/* ------------------------------ */


/* ------------------------------ */
struct AnnotationRec {
  bool DrawAnnotation = false;
  double AnnotationLon;
  double AnnotationLat;
  std::string AnnotationText;
};
/* ------------------------------ */


/* ------------------------------ */
struct PairLL {
  double eLon;
  double eLat;
};
/* ------------------------------ */


/* ------------------------------ */
struct SeqLineSegment {
  std::vector<PairLL> ListPairLL;
  bool IsClosed;
};
/* ------------------------------ */


/* ------------------------------ */
void ADD_ANNOTATION_TEXT(std::ofstream & os, AnnotationRec const& TheAnnot)
{
  if (TheAnnot.DrawAnnotation) {
    os << "  label=\"" << TheAnnot.AnnotationText << "\"\n";
    os << "  Xpos=" << TheAnnot.AnnotationLon << "\n";
    os << "  Ypos=" << TheAnnot.AnnotationLat << "\n";
    os << "  txres             = True                         ; Text resources desired\n";
    os << "  txres@txFont        = \"helvetica\"\n";
    os << "  txres@txFontHeightF=0.02\n";
    os << "  text = gsn_add_text(wks,plot,label, Xpos, Ypos, txres)\n";
  }
}
/* ------------------------------ */


/* ------------------------------ */
void GET_DOUBLE(mpq_class const& eQ, double & eD)
{
  eD=eQ.get_d();
}
/* ------------------------------ */


/* ------------------------------ */
void GET_DOUBLE(double const& eQ, double & eD)
{
  eD=eQ;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T1, typename T2>
void TYPE_CONVERSION(T1 const& a1, T2 & a2)
{
  a2=a1;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T2>
void TYPE_CONVERSION(int const& a1, T2 & a2)
{
  a2=a1;
}
/* ------------------------------ */


/* ------------------------------ */
template<>
void TYPE_CONVERSION(mpq_class const& a1, int & a2)
{
  mpz_class a1_z=a1.get_num();
  long a1_long=a1_z.get_si();
  a2=a1_long;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
T T_abs(T const& eVal)
{
  if (eVal > 0)
    return eVal;
  T fVal= - eVal;
  return fVal;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
T T_min(T const& eVal1, T const& eVal2)
{
  if (eVal1 > eVal2)
    return eVal2;
  return eVal1;
}
/* ------------------------------ */


/* ------------------------------ */
std::string random_string( size_t length )
{
  srand ( time(NULL) );
  auto randchar = []() -> char
    {
      const char charset[] =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";
      const size_t max_index = (sizeof(charset) - 1);
      return charset[ rand() % max_index ];
    };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
struct CollectedResult {
  std::vector<T> LVal;
  std::vector<int> LMult;
};
/* ------------------------------ */


/* ------------------------------ */
int NextIdx(int const& len,int const& i)
{
  if (i == len-1)
    return 0;
  return i+1;
}
/* ------------------------------ */


/* ------------------------------ */
int PrevIdx(int const& len,int const& i)
{
  if (i == 0)
    return len-1;
  return i-1;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<std::string> FILE_GetDirectoryListFile(std::string const& eDir)
{
  std::string ePath=eDir + ".";
  DIR* dirp=opendir(ePath.c_str());
  struct dirent *dp;
  std::vector<std::string> ListFile;
  while ((dp = readdir(dirp)) != NULL) {
    std::string eName=dp->d_name;
    //    free(dp); // not sure this is portable
    if (eName != ".." && eName != ".")
      ListFile.push_back(eName);
  }
  int err=closedir(dirp);
  if (err != 0) {
    std::cerr << "err=" << err << "\n";
    printf("Oh dear, something went wrong with ls! %s\n", strerror(errno));
    exit(1);
  }
  return ListFile;
}
/* ------------------------------ */


/* ------------------------------ */
bool FILE_IsRegularFile(std::string const& eFile)
{
  int status;
  struct stat st_buf;  
  status = stat(eFile.c_str(), &st_buf);
  if (status != 0) {
    std::cerr << "Problem in FILE_IsRegularFile\n";
    std::cerr << "Error, errno = " << errno << "\n";
    std::cerr << "eFile=" << eFile << "\n";
    exit(1);
  }
  if (S_ISREG (st_buf.st_mode)) {
    return true;
  }
  return false;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<std::string> FILE_GetDirectoryFilesRecursively(std::string const& eDir)
{
  std::cerr << "Beginning of FILE_GetDirectoryFilesRecursively\n";
  std::vector<std::string> ListDir={eDir};
  std::vector<std::string> ListFile;
  while(1) {
    std::vector<std::string> NewListDir;
    for (auto & eDir : ListDir) {
      std::vector<std::string> LocalListFile=FILE_GetDirectoryListFile(eDir);
      for (auto & eFile : LocalListFile) {
	std::string NewEnt=eDir + eFile;
	if (FILE_IsRegularFile(NewEnt) == true) {
	  ListFile.push_back(NewEnt);
	}
	else {
	  std::string NewDir=NewEnt + "/";
	  NewListDir.push_back(NewDir);
	}
      }
    }
    if (NewListDir.size() == 0)
      break;
    ListDir=NewListDir;
  }
  std::cerr << "Ending of FILE_GetDirectoryFilesRecursively\n";
  return ListFile;
}
/* ------------------------------ */


/* ------------------------------ */
bool IsExistingDirectory(std::string const& ThePrefix)
{
  if (0 != access(ThePrefix.c_str(), F_OK)) {
    if (ENOENT == errno) {
      // does not exist
      return false;
    }
    if (ENOTDIR == errno) {
      return false;
      // not a directory
    }
    std::cerr << "Should not happen a priori\n";
    exit(1);
  }
  return true;
}
/* ------------------------------ */


/* ------------------------------ */
bool IsExistingFile(std::string const& eFile)
{
  std::ifstream f(eFile.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }   
}
/* ------------------------------ */


/* ------------------------------ */
void RemoveFile(std::string const& eFile)
{
  std::remove(eFile.c_str());
}
/* ------------------------------ */


/* ------------------------------ */
void RemoveFileIfExist(std::string const& eFile)
{
  if (IsExistingFile(eFile))
    RemoveFile(eFile);
}
/* ------------------------------ */


/* ------------------------------ */
void RemoveFileSpecificExtension(std::string const& ThePrefix, std::string const& TheExtension)
{
  bool test=IsExistingDirectory(ThePrefix);
  if (test == false)
    return;
  std::vector<std::string> ListFile=FILE_GetDirectoryListFile(ThePrefix);
  int nbCharEnd=TheExtension.size();
  for (auto & eFile : ListFile) {
    int len=eFile.size();
    if (len > nbCharEnd) {
      std::string TheEnd=eFile.substr(len-nbCharEnd, nbCharEnd);
      if (TheEnd == TheExtension) {
	std::string eFileTot=ThePrefix + eFile;
	//	std::cerr << "Removing eFileTot=" << eFileTot << "\n";
	RemoveFile(eFileTot);
      }
    }
  }
}
/* ------------------------------ */


/* ------------------------------ */
#ifndef WINDOWS
std::string GetCurrentDirectory()
{
  int size = pathconf(".", _PC_PATH_MAX);
  char *buf;
  char *ptr;
  if ((buf = (char *)malloc((size_t)size)) != NULL) {
    ptr=getcwd(buf, (size_t)size);
    //    std::cerr << "ptr=" << (void*)ptr << "\n";
    //    std::cerr << "buf=" << (void*)buf << "\n";

    if (ptr == NULL && errno != ERANGE) {
      std::cerr << "Error while trying to use getcwd\n";
      exit(1);
    }
    //    std::cerr << "ptr=" << ptr << "\n";
    std::string eRet = buf;
    eRet=eRet + "/";
    free(buf);
    if (ptr != NULL) {
      if (ptr != buf) {
	std::cerr << "Before ptr freeing\n";
	free(ptr);
	std::cerr << "After ptr freeing\n";
      }
    }
    return eRet;
  }
  else {
    std::cerr << "Not enough memory\n";
    exit(1);
  }
}
#endif
/* ------------------------------ */


/* ------------------------------ */
#ifndef WINDOWS
std::string FILE_GetAbsoluteDirectory(std::string const& ePrefix)
{
  std::string FirstChar=ePrefix.substr(0, 1);
  if (FirstChar == "/") {
    return ePrefix;
  }
  else {
    std::string ePWD=GetCurrentDirectory();
    return ePWD + ePrefix;
  }
}
#endif
/* ------------------------------ */


/* ------------------------------ */
void CreateDirectory(std::string const& eDir)
{
  std::cerr << "eDir=" << eDir << "\n";
  const char *dir=eDir.c_str();
  char tmp[256];
  char *p = NULL;
  size_t len;
  snprintf(tmp, sizeof(tmp),"%s",dir);
  len = strlen(tmp);
  if(tmp[len - 1] == '/')
    tmp[len - 1] = 0;
  for(p = tmp + 1; *p; p++)
    if(*p == '/') {
      *p = 0;
      //      std::cerr << "Before mkdir, tmp=" << tmp << "\n";
      mkdir(tmp, S_IRWXU);
      *p = '/';
    }
  //  std::cerr << "Before mkdir, tmp=" << tmp << "\n";
  mkdir(tmp, S_IRWXU);
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<std::string> ls_operation(std::string const& ThePrefix)
{
  std::string TmpFile="/tmp/file" + random_string(20);
  std::string eOrder="ls " + ThePrefix + " > " + TmpFile;
  int iret=system(eOrder.c_str() );
  if (iret == -1) {
    std::cerr << "Error in ls_operation\n";
    std::cerr << "ThePrefix=" << ThePrefix << "\n";
    std::cerr << "unable to run the process\n";
    exit(1);
  }
  //
  std::ifstream os;
  os.open(TmpFile);
  std::vector<std::string> ListFile;
  while(1) {
    if (os.eof()) {
      os.close();
      RemoveFile(TmpFile);
      return ListFile;
    }
    std::string eFile;
    os >> eFile;
    ListFile.push_back(eFile);
  }
}
/* ------------------------------ */


/* ------------------------------ */
bool STRING_IsStringReduceToSpace(std::string const& eStr)
{
  int len=eStr.length();
  std::string eChar=" ";
  for (int i=0; i<len; i++) {
    std::string eSubChar=eStr.substr(i,1);
    if (eSubChar != eChar)
      return false;
  }
  return true;
}
/* ------------------------------ */


/* ------------------------------ */
int STRING_GetCharPositionInString(std::string const& eStr, std::string const& eChar)
{
  int len=eStr.length();
  for (int i=0; i<len; i++) {
    std::string eSubChar=eStr.substr(i,1);
    if (eSubChar == eChar)
      return i;
  }
  return -1;
}
/* ------------------------------ */


/* ------------------------------ */
std::string DoubleTo4dot2f(double const& x)
{
  char buffer[50];
  int n=sprintf(buffer, "%4.2f", x);
  if (n == 0) {
    std::cerr << "Clear error in DoubleTo4dot2f\n";
    exit(1);
  }
  return std::string(buffer);
}
/* ------------------------------ */


/* ------------------------------ */
std::string IntToString(int const & x)
{
  std::stringstream s;
  s << x;
  std::string converted(s.str());
  return converted;
}
/* ------------------------------ */


/* ------------------------------ */
std::string LongToString(long const & x)
{
  std::stringstream s;
  s << x;
  std::string converted(s.str());
  return converted;
}
/* ------------------------------ */


/* ------------------------------ */
std::string StringNumber(int const& nb, int const& nbDigit)
{
  if (nb > pow(10, nbDigit) ) {
    std::stringstream s;
    s << "Critical error in StringNumber\n";
    s << "nb=" << nb << "\n";
    s << "nbDigit=" << nbDigit << "\n";
    std::string eStr(s.str());
    throw eStr;
  }
  int idx=1;
  while(1) {
    if (nb < pow(10, idx) ) {
      std::string TheStr="";
      for (int i=0; i<nbDigit-idx; i++) {
	TheStr=TheStr + "0";
      }
      TheStr=TheStr + IntToString(nb);
      return TheStr;
    }
    idx++;
  }
}
/* ------------------------------ */


/* ------------------------------ */
std::string STRING_RemoveSpacesBeginningEnd(std::string const& eStr)
{
  int len=eStr.size();
  int ListIsSpace[len];
  std::string eSpace=" ";
  for (int i=0; i<len; i++) {
    std::string eChar=eStr.substr(i, 1);
    if (eChar == eSpace)
      ListIsSpace[i]=1;
    else
      ListIsSpace[i]=0;
  }
  int PosLow=-1;
  for (int i=0; i<len; i++)
    if (PosLow == -1)
      if (ListIsSpace[i] == 0)
	PosLow=i;
  int PosUpp=-1;
  for (int i=0; i<len; i++) {
    int j=len-1-i;
    if (PosUpp == -1)
      if (ListIsSpace[j] == 0)
	PosUpp=j;
  }
  std::string RetStr;
  if (PosLow == -1) {
    return RetStr;
  }
  for (int iPos=PosLow; iPos<PosUpp+1; iPos++) {
    RetStr=RetStr + eStr.at(iPos);
  }
  return RetStr;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<std::string> STRING_Split(std::string const& eStrA, std::string const& eStrB)
{
  int lenA=eStrA.length();
  int lenB=eStrB.length();
  int ListStatus[lenA];
  //  std::cerr << "eStrA=" << eStrA << "\n";
  //  std::cerr << "eStrB=" << eStrB << "\n";
  for (int iA=0; iA<lenA; iA++)
    ListStatus[iA]=1;
  for (int iA=0; iA<lenA - lenB; iA++)
    if (ListStatus[iA] == 1) {
      int IsMatch=1;
      for (int iB=0; iB<lenB; iB++) {
	std::string eCharA=eStrA.substr(iA+iB,1);
	std::string eCharB=eStrB.substr(iB,1);
	if (eCharA != eCharB) {
	  IsMatch=0;
	}
      }
      if (IsMatch == 1)
	for (int iB=0; iB<lenB; iB++)
	  ListStatus[iA + iB]=0;
    }
  std::vector<std::string> RetList;
  std::string eFound;
  for (int iA=0; iA<lenA; iA++) {
    std::string eChar=eStrA.substr(iA, 1);
    if (ListStatus[iA] == 1) {
      eFound=eFound + eChar;
    }
    if (ListStatus[iA] == 0) {
      int siz=eFound.length();
      if (siz > 0) {
	RetList.push_back(eFound);
      }
      eFound="";
    }
  }
  int siz=eFound.size();
  if (siz > 0) {
    RetList.push_back(eFound);
  }
  //  std::cerr << "Leaving STRING_Split\n";
  //  int nbStr=RetList.size();
  //  for (int iStr=0; iStr<nbStr; iStr++)
  //    std::cerr << "iStr=" << iStr << " eStr=" << RetList[iStr] << "\n";
  return RetList;
}
/* ------------------------------ */


/* ------------------------------ */
std::string FILE_GetExtension(std::string const& eFile)
{
  std::vector<std::string> LStr=STRING_Split(eFile, "/");
  std::string eFinal=LStr[LStr.size()-1];
  std::vector<std::string> LBlck=STRING_Split(eFile, ".");
  return LBlck[LBlck.size()-1];
}
/* ------------------------------ */


/* ------------------------------ */
bool NC_IsVar(std::string const& eFile, std::string const& eVar)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Error in NC_IsVar\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile=" << eFile << "\n";
    exit(1);
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if(data.isNull()) {
    return false;
  }
  return true;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> NC_Read2Dvariable(std::string const& eFile, std::string const& eVar)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Error in NC_Read2Dvariable\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile=" << eFile << "\n";
    exit(1);
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  netCDF::NcType eType=data.getType();
  if(data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "eVar=" << eVar << "\n";
    exit(1);
  }
  int nbDim=data.getDimCount();
  if (nbDim != 2) {
    std::cerr << "The number of dimensions is not correct\n";
    exit(1);
  }
  netCDF::NcDim eDim=data.getDim(0);
  int eta=eDim.getSize();
  netCDF::NcDim fDim=data.getDim(1);
  int xi=fDim.getSize();
  MyMatrix<double> eArr(eta, xi);
  int IsMatch=false;
  if (eType == netCDF::NcType::nc_DOUBLE) {
    double *eVal;
    eVal=new double[eta*xi];
    data.getVar(eVal);
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	eArr(i,j)=eVal[idx];
	idx++;
      }
    delete [] eVal;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    float *eValFLOAT;
    eValFLOAT=new float[eta*xi];
    data.getVar(eValFLOAT);
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	float eValF=eValFLOAT[idx];
	double eValD=double(eValF);
	eArr(i,j)=eValD;
	idx++;
      }
    delete [] eValFLOAT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_INT) {
    int *eValINT;
    eValINT=new int[eta*xi];
    data.getVar(eValINT);
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	int eValI=eValINT[idx];
	double eValD=double(eValI);
	eArr(i,j)=eValD;
	idx++;
      }
    delete [] eValINT;
    IsMatch=true;
  }
  if (IsMatch == false) {
    std::cerr << "Did not find the right number type\n";
    exit(1);
  }
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<int> NC_Read2Dvariable_int(std::string const& eFile, std::string const& eVar)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Error in NC_Read2Dvariable_int\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile=" << eFile << "\n";
    exit(1);
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  netCDF::NcType eType=data.getType();
  if(data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "eVar=" << eVar << "\n";
    exit(1);
  }
  int nbDim=data.getDimCount();
  if (nbDim != 2) {
    std::cerr << "The number of dimensions is not correct\n";
    exit(1);
  }
  netCDF::NcDim eDim=data.getDim(0);
  int eta=eDim.getSize();
  netCDF::NcDim fDim=data.getDim(1);
  int xi=fDim.getSize();
  MyMatrix<int> eArr(eta, xi);
  if (eType == netCDF::NcType::nc_INT) {
    int *eValINT;
    eValINT=new int[eta*xi];
    data.getVar(eValINT);
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	int eValI=eValINT[idx];
	double eValD=double(eValI);
	eArr(i,j)=eValD;
	idx++;
      }
    delete [] eValINT;
  }
  else {
    std::cerr << "Error in the call\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "eVar=" << eVar << "\n";
    exit(1);
  }
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
MyVector<double> NC_Read1Dvariable(std::string const& eFile, std::string const& eVar)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Error in NC_Read1Dvariable\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile=" << eFile << "\n";
    exit(1);
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  netCDF::NcType eType=data.getType();
  if(data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "eVar=" << eVar << "\n";
    exit(1);
  }
  int nbDim=data.getDimCount();
  if (nbDim != 1) {
    std::cerr << "The number of dimensions is not correct\n";
    exit(1);
  }
  netCDF::NcDim eDim=data.getDim(0);
  int dim=eDim.getSize();
  MyVector<double> eArr(dim);
  // Fetching the variable from the file
  bool IsMatch=false;
  if (eType == netCDF::NcType::nc_DOUBLE) {
    double *eVal;
    eVal=new double[dim];
    data.getVar(eVal);
    for (int i=0; i<dim; i++)
      eArr(i)=eVal[i];
    delete [] eVal;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    float *eValFLOAT;
    eValFLOAT=new float[dim];
    data.getVar(eValFLOAT);
    for (int i=0; i<dim; i++) {
      float eValF=eValFLOAT[i];
      double eValD=double(eValF);
      eArr(i)=eValD;
    }
    delete [] eValFLOAT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_INT) {
    int *eValINT;
    eValINT=new int[dim];
    data.getVar(eValINT);
    for (int i=0; i<dim; i++) {
      int eValI=eValINT[i];
      double eValD=double(eValI);
      eArr(i)=eValD;
    }
    delete [] eValINT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_INT) {
    int *eValINT;
    eValINT=new int[dim];
    data.getVar(eValINT);
    for (int i=0; i<dim; i++) {
      int eValI=eValINT[i];
      double eValD=double(eValI);
      eArr(i)=eValD;
    }
    delete [] eValINT;
    IsMatch=true;
  }
  if (eType == netCDF::NcType::nc_SHORT) {
    signed short int *eValINT;
    eValINT=new signed short int[dim];
    data.getVar(eValINT);
    for (int i=0; i<dim; i++) {
      double eValD=double(eValINT[i]);
      eArr(i)=eValD;
    }
    delete [] eValINT;
    IsMatch=true;
  }
  if (IsMatch == false) {
    std::cerr << "Did not find any matching number type\n";
    exit(1);
  }
  // Now reading the offset and scaling_factor
  double eScal, eOff;
  try {
    //    std::cerr << "Before reading scale_factor\n";
    netCDF::NcVarAtt eScalAtt=data.getAtt("scale_factor");
    //    std::cerr << "After reading scale_factor\n";
    if (eScalAtt.isNull()) {
      eScal=1;
    }
    else {
      eScalAtt.getValues(&eScal);
    }
  }
  catch (...) {
    eScal=1;
  }
  try {
    netCDF::NcVarAtt eOffAtt=data.getAtt("add_offset");
    if (eOffAtt.isNull()) {
      eOff=0;
    }
    else {
      eOffAtt.getValues(&eOff);
    }
  }
  catch (...) {
    eOff=0;
  }
  for (int i=0; i<dim; i++) {
    eArr(i)=eOff + eScal*eArr(i);
  }
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
MyVector<int> NC_Read1Dvariable_int(std::string const& eFile, std::string const& eVar)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Error in NC_Read1Dvariable_int\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile=" << eFile << "\n";
    exit(1);
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  netCDF::NcType eType=data.getType();
  if(data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "eVar=" << eVar << "\n";
    exit(1);
  }
  int nbDim=data.getDimCount();
  if (nbDim != 1) {
    std::cerr << "The number of dimensions is not correct\n";
    exit(1);
  }
  netCDF::NcDim eDim=data.getDim(0);
  int dim=eDim.getSize();
  //  std::cerr << "dim=" << dim << "\n";
  MyVector<int> eArr(dim);
  int *eValINT;
  eValINT=new int[dim];
  data.getVar(eValINT);
  for (int i=0; i<dim; i++) {
    int eValI=eValINT[i];
    //    std::cerr << "i=" << i << " eValI=" << eValI << "\n";
    eArr(i)=eValI;
    //  std::cerr << "After the write\n";
  }
  delete [] eValINT;
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> My_u2rho(MyMatrix<double> const& eVar_u, MyMatrix<int> const& MSK_rho)
{
  int eta_rho=MSK_rho.rows();
  int xi_rho=MSK_rho.cols();
  int eta_u=eVar_u.rows();
  int xi_u=eVar_u.cols();
  if (eta_u != eta_rho || xi_u != xi_rho-1) {
    std::cerr << "Dimension error in My_u2rho\n";
    exit(1);
  }
  MyMatrix<double> eVar_rho(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eSumMsk=0;
      double eSumVal=0;
      if (MSK_rho(i,j) == 1) {
	eSumMsk++;
	eSumVal=eSumVal + eVar_u(i,j);
      }
      if (j < xi_u) {
	if (MSK_rho(i,j+1) == 1) {
	  eSumMsk++;
	  eSumVal=eSumVal + eVar_u(i,j+1);
	}
      }
      if (eSumMsk == 0) {
	eVar_rho(i,j)=0;
      }
      else {
	double eVal=eSumVal/double(eSumMsk);
	eVar_rho(i,j)=eVal;
      }
    }
  return eVar_rho;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> My_v2rho(MyMatrix<double> const& eVar_v, MyMatrix<int> const& MSK_rho)
{
  int eta_rho=MSK_rho.rows();
  int xi_rho=MSK_rho.cols();
  int eta_v=eVar_v.rows();
  int xi_v=eVar_v.cols();
  if (eta_v != eta_rho-1 || xi_v != xi_rho) {
    std::cerr << "Dimension error in My_v2rho\n";
    exit(1);
  }
  MyMatrix<double> eVar_rho(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      int eSumMsk=0;
      double eSumVal=0;
      if (MSK_rho(i,j) == 1) {
	eSumMsk++;
	eSumVal=eSumVal + eVar_v(i,j);
      }
      if (i < eta_v) {
	if (MSK_rho(i+1,j) == 1) {
	  eSumMsk++;
	  eSumVal=eSumVal + eVar_v(i+1,j);
	}
      }
      if (eSumMsk == 0) {
	eVar_rho(i,j)=0;
      }
      else {
	double eVal=eSumVal/double(eSumMsk);
	eVar_rho(i,j)=eVal;
      }
    }
  return eVar_rho;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<std::string> GRIB_GetAllFilesInDirectory(std::string const& ePrefix)
{
  std::vector<std::string> ListFile=FILE_GetDirectoryFilesRecursively(ePrefix);
  std::vector<std::string> RetListFile;
  for (auto & eFile : ListFile) {
    int len=eFile.size();
    if (len > 3) {
      std::string eEnd=eFile.substr(len-3,3);
      if (eEnd == "grb")
	RetListFile.push_back(eFile);
    }
  }
  return RetListFile;
}
/* ------------------------------ */


/* ------------------------------ */
struct CosmoGridInfo {
  double latitudeOfSouthernPoleInDegrees;
  double longitudeOfSouthernPoleInDegrees;
  double angleOfRotationInDegrees;
  //
  double longitudeOfFirstGridPointInDegrees;
  double latitudeOfFirstGridPointInDegrees;
  double longitudeOfLastGridPointInDegrees;
  double latitudeOfLastGridPointInDegrees;
  //
  double iDirectionIncrementInDegrees;
  double jDirectionIncrementInDegrees;
};
/* ------------------------------ */


/* ------------------------------ */
double phirot2phi(double const& phirot, double const& rlarot, double const& polphi, double const& pollam, double const& polgam)
{
  double pi=3.1415926535;
  double eMult=pi/double(180);
  double eMultInv=double(180)/pi;
  double zsinpol = sin(eMult * polphi);
  double zcospol = cos(eMult * polphi);
  double zphis  = eMult * phirot;
  double zrlas;
  if (rlarot > double(180)) {
    zrlas = rlarot - double(360);
  }
  else {
    zrlas = rlarot;
  }
  zrlas = eMult * zrlas;
  double zarg;
  if (fabs(polgam) > 0) {
    double zgam  = eMult * polgam;
    zarg = zsinpol*sin(zphis) + zcospol*cos(zphis) * ( cos(zrlas)*cos(zgam) - sin(zgam)*sin(zrlas));
  }
  else {
    zarg = zcospol * cos(zphis) * cos(zrlas) + zsinpol * sin(zphis);
  }
  double phirot2phi  = eMultInv * asin(zarg);
  return phirot2phi;
}
/* ------------------------------ */


/* ------------------------------ */
double rlarot2rla(double const& phirot, double const& rlarot, double const& polphi, double const& pollam, double const& polgam)
{
  double pi=3.1415926535;
  double eMult=pi/double(180);
  double eMultInv=double(180)/pi;
  double zsinpol = sin(eMult * polphi);
  double zcospol = cos(eMult * polphi);
  double zphis  = eMult * phirot;
  double zrlas;
  if (rlarot > double(180)) {
    zrlas = rlarot - double(360);
  }
  else {
    zrlas = rlarot;
  }
  zrlas = eMult * zrlas;
  double zlampol = eMult * pollam;
  double zarg1, zarg2;
  if (fabs(polgam) > 0) {
    double zgam    = eMult * polgam;
    zarg1   = sin (zlampol) *                                            
      (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam)) 
       + zcospol * sin(zphis))                                              
      - cos (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam));
    zarg2   = cos (zlampol) *                                               
      (- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))  
       + zcospol * sin(zphis))                                               
      + sin (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam));
  }
  else {
    zarg1   = sin (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis)  +    
			       zcospol *              sin(zphis)) -    
      cos (zlampol) *             sin(zrlas) * cos(zphis);
    zarg2   = cos (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis)  +    
			       zcospol *              sin(zphis)) +   
      sin (zlampol) *             sin(zrlas) * cos(zphis);
  }
  if (zarg2 == 0) zarg2=1.0e-20;
  double rlarot2rla = eMultInv * atan2(zarg1,zarg2);
  return rlarot2rla;
}
/* ------------------------------ */


/* ------------------------------ */
void Apply_COSMO_Transformation(MyMatrix<double> & LON, MyMatrix<double> & LAT, CosmoGridInfo const& eCosmoGrid)
{
  double pollat_sp=eCosmoGrid.latitudeOfSouthernPoleInDegrees;
  double pollon_sp=eCosmoGrid.longitudeOfSouthernPoleInDegrees;
  double polgam=eCosmoGrid.angleOfRotationInDegrees;
  double zstartlon_tot=eCosmoGrid.longitudeOfFirstGridPointInDegrees;
  double zstartlat_tot=eCosmoGrid.latitudeOfFirstGridPointInDegrees;
  double zendlon_tot=eCosmoGrid.longitudeOfLastGridPointInDegrees;
  double zendlat_tot=eCosmoGrid.latitudeOfLastGridPointInDegrees;
  double dlon=eCosmoGrid.iDirectionIncrementInDegrees;
  double dlat=eCosmoGrid.jDirectionIncrementInDegrees;
  if (zendlon_tot == zstartlon_tot || zendlat_tot == zstartlat_tot) {
    std::cerr << "Error of consistency in zstartlat / zendlat\n";
    exit(1);
  }
  //
  int eta_rho=LON.rows();
  int xi_rho=LON.cols();
  double pollat= - pollat_sp;
  double pollon= pollon_sp - double(180);
  double startlon_tot=zstartlon_tot; // This part is quite unsure. Maybe there is a shift
  double startlat_tot=zstartlat_tot; // same remark
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      double eLonR=startlon_tot + double(i-1)*dlon;
      double eLatR=startlat_tot + double(j-1)*dlat;
      double eLat=phirot2phi(eLatR, eLonR, pollat, pollon, polgam);
      double eLon=rlarot2rla(eLatR, eLonR, pollat, pollon, polgam);
      LON(i,j)=eLon;
      LAT(i,j)=eLat;
    }
  std::cerr << "After LON min=" << LON.minCoeff() << " max=" << LON.maxCoeff() << "\n";
  std::cerr << "After LAT min=" << LAT.minCoeff() << " max=" << LAT.maxCoeff() << "\n";
}
/* ------------------------------ */


/* ------------------------------ */
struct GRIB_MessageInfo {
  std::string shortName;
  std::string name;
  int idx;
  double time;
  double timeStart;
  int stepRange;
};
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> GRIB_Read2DVariable(std::vector<std::string> const& ListFileName, std::string const& VarName)
{
  for (auto & eFileName : ListFileName) {
    grib_handle *h = NULL;
    int err;
    FILE* in = NULL;
    in = fopen(eFileName.c_str(),"r");
    unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;
    while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) { 
      //    std::cerr << "err=" << err << "\n";
      if (err != GRIB_SUCCESS)
	GRIB_CHECK(err,0);
      //
      char name_space[3]="ls";
      grib_keys_iterator* kiter=NULL;
      kiter=grib_keys_iterator_new(h,key_iterator_filter_flags,name_space);
      if (!kiter) {
	printf("ERROR: Unable to create keys iterator\n");
	exit(1);
      }
      std::string ShortNameValue;
      while(grib_keys_iterator_next(kiter)) {
	int MAX_VAL_LEN=1024;
	char value[MAX_VAL_LEN];
	size_t vlen=MAX_VAL_LEN;
	const char* name = grib_keys_iterator_get_name(kiter);
	vlen=MAX_VAL_LEN;
	bzero(value,vlen);
	GRIB_CHECK(grib_get_string(h,name,value,&vlen),name);
	std::string nameStr=name;
	std::string valueStr=value;
	if (nameStr == "shortName") {
	  ShortNameValue=valueStr;
	}
	//	std::cerr << "name=" << nameStr << " value=" << valueStr << "\n";
	//printf("%s = %s\n",name,value);
      }
      grib_keys_iterator_delete(kiter);
      //      std::cerr << "ShortNameValue=" << ShortNameValue << "\n";
      //
      if (ShortNameValue == VarName) {
	long Ni, Nj, numberOfDataPoints;
	GRIB_CHECK(grib_get_long(h,"Ni",&Ni),0);
	//	std::cerr << "Ni=" << Ni << "\n";
	GRIB_CHECK(grib_get_long(h,"Nj",&Nj),0);
	//	std::cerr << "Nj=" << Nj << "\n";
	GRIB_CHECK(grib_get_long(h,"numberOfDataPoints",&numberOfDataPoints),0);
	//	std::cerr << "NumberOfDataPoints=" << numberOfDataPoints << "\n";
	double *lats, *lons, *values;
	size_t size=numberOfDataPoints;
	size_t* sizePtr=&size;
	lats=(double*)malloc(size*sizeof(double));
	lons=(double*)malloc(size*sizeof(double));
	values=(double*)malloc(size*sizeof(double));
	err=grib_get_data(h, lats, lons, values, sizePtr);
	grib_handle_delete(h);
	if (err != GRIB_SUCCESS)
	  GRIB_CHECK(err,0);
	int eta_rho=Ni;
	int xi_rho=Nj;
	MyMatrix<double> VAL(eta_rho, xi_rho);
	//	std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
	int idx=0;
	for (int j=0; j<xi_rho; j++)
	  for (int i=0; i<eta_rho; i++) {
	    VAL(i,j)=values[idx];
	    idx++;
	  }
	free(lats);
	free(lons);
	free(values);
	fclose(in);
	return VAL;
      }
      grib_handle_delete(h);
    }
    fclose(in);
  }
  std::cerr << "Error in GRIB_Read2DVariable\n";
  std::cerr << "Failed to find the variable =" << VarName << "\n";
  std::cerr << "Exiting\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
double GeodesicDistance(double const& LonDeg1, double const& LatDeg1, double const& LonDeg2, double const& LatDeg2)
{
  double pi=3.141592653589792;
  double lon1=pi*LonDeg1/double(180);
  double lat1=pi*LatDeg1/double(180);
  double x1=cos(lon1)*cos(lat1);
  double y1=sin(lon1)*cos(lat1);
  double z1=sin(lat1);

  double lon2=pi*LonDeg2/double(180);
  double lat2=pi*LatDeg2/double(180);
  double x2=cos(lon2)*cos(lat2);
  double y2=sin(lon2)*cos(lat2);
  double z2=sin(lat2);
  double scalprod=x1*x2+y1*y2+z1*z2;
  if (scalprod > 1)
    return 0;
  else
    return acos(scalprod);
}
/* ------------------------------ */


/* ------------------------------ */
double GeodesicDistanceKM(double const& LonDeg1, double const& LatDeg1, double const& LonDeg2, double const& LatDeg2)
{
  double EarthRadius=6370;
  return EarthRadius*GeodesicDistance(LonDeg1, LatDeg1, LonDeg2, LatDeg2);
}
/* ------------------------------ */


/* ------------------------------ */
struct Point {
  double eLon;
  double eLat;
};
/* ------------------------------ */


/* ------------------------------ */
struct DiscInfo {
  Point eSample;
  Point avgPoint;
  double SpreadLon;
  double SpreadLat;
  double MaxSpread;
};
/* ------------------------------ */


/* ------------------------------ */
struct KTreeElt {
  std::vector<Point> ListPt;
  DiscInfo eDisc;
  //
  bool IsSplit;
  int iSub1;
  int iSub2;
};
/* ------------------------------ */


/* ------------------------------ */
DiscInfo KTree_ComputeDisc(std::vector<Point> const& ListPt)
{
  double MinLon=ListPt[0].eLon;
  double MaxLon=ListPt[0].eLon;
  double MinLat=ListPt[0].eLat;
  double MaxLat=ListPt[0].eLat;
  int siz=ListPt.size();
  double SumLon=0;
  double SumLat=0;
  for (int i=0; i<siz; i++) {
    double eLon=ListPt[i].eLon;
    double eLat=ListPt[i].eLat;
    SumLon += eLon;
    SumLat += eLat;
    if (eLon < MinLon)
      MinLon=eLon;
    if (eLon > MaxLon)
      MaxLon=eLon;
    if (eLat < MinLat)
      MinLat=eLat;
    if (eLat > MaxLat)
      MaxLat=eLat;
  }
  double avgLon=SumLon/double(siz);
  double avgLat=SumLat/double(siz);
  double dist12=GeodesicDistanceKM(MinLon, MinLat, MinLon, MaxLat);
  double dist23=GeodesicDistanceKM(MinLon, MaxLat, MaxLon, MaxLat);
  double dist34=GeodesicDistanceKM(MaxLon, MaxLat, MaxLon, MinLat);
  double dist41=GeodesicDistanceKM(MaxLon, MinLat, MinLon, MinLat);
  //
  double SpreadLat=(dist12 + dist34)/double(2);
  double SpreadLon=(dist23 + dist41)/double(2);
  double MinDistKM=20000;
  int idxMin=-1;
  for (int i=0; i<siz; i++) {
    double eLon=ListPt[i].eLon;
    double eLat=ListPt[i].eLat;
    double dist=GeodesicDistanceKM(eLon, eLat, avgLon, avgLat);
    if (dist < MinDistKM) {
      MinDistKM=dist;
      idxMin=i;
    }
  }
  Point eSample=ListPt[idxMin];
  double eSampleLon=eSample.eLon;
  double eSampleLat=eSample.eLat;
  double MaxSpread=0;
  for (int i=0; i<siz; i++) {
    double eLon=ListPt[i].eLon;
    double eLat=ListPt[i].eLat;
    double dist=GeodesicDistanceKM(eLon, eLat, eSampleLon, eSampleLat);
    if (dist > MaxSpread)
      MaxSpread=dist;
  }
  Point avgPoint{avgLon, avgLat};
  return {eSample, avgPoint, SpreadLon, SpreadLat, MaxSpread};
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<KTreeElt> SplitKTreeElt(KTreeElt const& eKD)
{
  std::vector<Point> NewListPt1;
  std::vector<Point> NewListPt2;
  if (eKD.eDisc.SpreadLon > eKD.eDisc.SpreadLat) {
    for (auto & ePt : eKD.ListPt) {
      if (ePt.eLon < eKD.eDisc.avgPoint.eLon)
	NewListPt1.push_back(ePt);
      else
	NewListPt2.push_back(ePt);
    }
  }
  else {
    for (auto & ePt : eKD.ListPt) {
      if (ePt.eLat < eKD.eDisc.avgPoint.eLat)
	NewListPt1.push_back(ePt);
      else
	NewListPt2.push_back(ePt);
    }
  }
  //
  KTreeElt eKD1{{}, eKD.eDisc, true, -2, -2};
  KTreeElt eKD2{NewListPt1, KTree_ComputeDisc(NewListPt1), false, -1, -1};
  KTreeElt eKD3{NewListPt2, KTree_ComputeDisc(NewListPt2), false, -1, -1};
  return {eKD1, eKD2, eKD3};
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<KTreeElt> KTree_GetDecomposition(std::vector<Point> const& ListPtCoast)
{
  std::vector<KTreeElt> TList;
  std::vector<int> IsTreated;
  int MaxNumberPerCell=100;
  std::function<void(int)> SplitComponent=[&](int const& iComp) -> void {
    int len=TList.size();
    std::vector<KTreeElt> NList=SplitKTreeElt(TList[iComp]);
    TList[iComp]=NList[0];
    TList[iComp].iSub1=len;
    TList[iComp].iSub2=len+1;
    TList.push_back(NList[1]);
    TList.push_back(NList[2]);
    IsTreated.push_back(0);
    IsTreated.push_back(0);
  };
  KTreeElt eElt{ListPtCoast, KTree_ComputeDisc(ListPtCoast), false, -1, -1};
  TList.push_back(eElt);
  IsTreated.push_back(0);
  while(true) {
    int siz=TList.size();
    //    std::cerr << "siz=" << siz << "\n";
    bool IsFinished=true;
    for (int i=0; i<siz; i++)
      if (IsTreated[i] == 0) {
	IsTreated[i]=1;
	IsFinished=false;
	int len=TList[i].ListPt.size();
	if (len > MaxNumberPerCell)
	  SplitComponent(i);
      }
    if (IsFinished)
      break;
  }
  //  for (int i=0; i<siz; i++) {
  //    std::cerr << "i=" << i << " split=" << TList[i].IsSplit << " |ListPt|=" << TList[i].ListPt.size() << " MaxSpread=" << TList[i].eDisc.MaxSpread << "\n";
  //  }
  return TList;
}
/* ------------------------------ */


/* ------------------------------ */
double ShortestDistance(std::vector<KTreeElt> const& ListKT, Point const& ePt, double & UpperEstimate)
{
  double TheDist=UpperEstimate;
  double eLon=ePt.eLon;
  double eLat=ePt.eLat;
  //  std::cerr << "|ListKT|=" << ListKT.size() << "\n";
  std::vector<int> ListIdx{0};
  while(true) {
    if (ListIdx.size() == 0)
      break;
    //    std::cerr << "Before creating NewListIdx\n";
    std::vector<int> NewListIdx;
    for (int & eVal : ListIdx) {
      if (ListKT[eVal].IsSplit == false) {
	int len=ListKT[eVal].ListPt.size();
	for (int i=0; i<len; i++) {
	  double nLon=ListKT[eVal].ListPt[i].eLon;
	  double nLat=ListKT[eVal].ListPt[i].eLat;
	  double nDist=GeodesicDistanceKM(nLon, nLat, eLon, eLat);
	  if (nDist < TheDist)
	    TheDist=nDist;
	}
      }
      else {
	double SampleLon=ListKT[eVal].eDisc.eSample.eLon;
	double SampleLat=ListKT[eVal].eDisc.eSample.eLat;
	double nDist=GeodesicDistanceKM(SampleLon, SampleLat, eLon, eLat);
	double LowerBound=nDist - ListKT[eVal].eDisc.MaxSpread;
	if (LowerBound < TheDist) {
	  NewListIdx.push_back(ListKT[eVal].iSub1);
	  NewListIdx.push_back(ListKT[eVal].iSub2);
	}
      }
    }
    ListIdx=NewListIdx;
  }
  return TheDist;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<double> GetUpperEstimateMinDist(std::vector<Point> const& ListPt1, std::vector<Point> const& ListPt2)
{
  int nbPt2=ListPt2.size();
  int nbPt1=ListPt1.size();
  if (nbPt1 == 0) {
    std::cerr << "The list ListPt1 should not be empty\n";
    std::cerr << "nbPt1=" << nbPt1 << "\n";
    exit(1);
  }
  //  std::cerr << "nbPt1=" << nbPt1 << " nbPt2=" << nbPt2 << "\n";
  std::function<double(int,int)> fDist=[&](int const& i1, int const& i2) -> double {
    double eLon1=ListPt1[i1].eLon;
    double eLat1=ListPt1[i1].eLat;
    double eLon2=ListPt2[i2].eLon;
    double eLat2=ListPt2[i2].eLat;
    return GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
  };
  int i1=0;
  std::vector<double> ListUpperEst(nbPt2);
  for (int i2=0; i2<nbPt2; i2++) {
    double eDist=fDist(i1, i2);
    //    std::cerr << "i2=" << i2 << "\n";
    for (int iter=0; iter<4; iter++) {
      int i1New =rand() % nbPt1;
      double nDist=fDist(i1New, i2);
      if (nDist < eDist) {
	eDist=nDist;
	i1=i1New;
	//	std::cerr << "  i1New=" << i1 << " eDist=" << eDist << "\n";
      }
    }
    while(true) {
      bool DoSomething=false;
      //      std::cerr << "  i1=" << i1 << " eDist=" << eDist << "\n";
      if (i1>0) {
	double nDist=fDist(i1-1,i2);
	if (nDist < eDist) {
	  eDist=nDist;
	  i1=i1 - 1;
	  DoSomething=true;
	}
      }
      if (i1 < nbPt1-1) {
	double nDist=fDist(i1+1,i2);
	if (nDist < eDist) {
	  eDist=nDist;
	  i1=i1 + 1;
	  DoSomething=true;
	}
      }
      if (DoSomething == false)
	break;
    }
    ListUpperEst[i2]=eDist;
  }
  //  std::cerr << "Now leaving returning ListUpperEst\n";
  return ListUpperEst;
}
/* ------------------------------ */


/* ------------------------------ */
void TwoPiNormalization(double & TheAng)
{
  double ThePi=3.141592653589792;
  if (TheAng < -ThePi) {
    TheAng += double(2)*ThePi;
  }
  if (TheAng > ThePi) {
    TheAng -= double(2)*ThePi;
  }
}
/* ------------------------------ */


/* ------------------------------ */
int MySign(double & TheVal)
{
  if (TheVal > 0)
    return 1;
  if (TheVal < 0)
    return -1;
  return 0;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> CreateAngleMatrix(MyMatrix<double> const& LON_rho, MyMatrix<double> const& LAT_rho)
{
  int eta_rho=LON_rho.rows();
  int xi_rho=LON_rho.cols();
  int eta_v=eta_rho-1;
  int xi_v=xi_rho;
  MyMatrix<double> LONrad_v(eta_v,xi_v);
  MyMatrix<double> LATrad_v(eta_v,xi_v);
  MyMatrix<double> azim(eta_v-1,xi_v);
  double ThePi=3.141592653589792;
  double DegTwoRad=ThePi/double(180);
  for (int iEta=0; iEta<eta_v; iEta++)
    for (int iXi=0; iXi<xi_v; iXi++) {
      double eLon=(LON_rho(iEta,iXi)+LON_rho(iEta+1,iXi))/double(2);
      double eLat=(LAT_rho(iEta,iXi)+LAT_rho(iEta+1,iXi))/double(2);
      LONrad_v(iEta,iXi)=eLon*DegTwoRad;
      LATrad_v(iEta,iXi)=eLat*DegTwoRad;
    }
  for (int iEta=1; iEta<eta_v-1; iEta++)
    for (int iXi=0; iXi<xi_v; iXi++) {
      double phi1=LATrad_v(iEta,iXi);
      double xlam1=LONrad_v(iEta,iXi);
      double phi2=LATrad_v(iEta+1,iXi);
      double xlam2=LONrad_v(iEta+1,iXi);
      double TPSI2=tan(phi2);
      double dlam=xlam2-xlam1;
      TwoPiNormalization(dlam);
      double cta12=(cos(phi1)*TPSI2 - sin(phi1)*cos(dlam))/sin(dlam);
      double eAzim=atan(double(1)/cta12);
      int signAzim=MySign(eAzim);
      int signDlam=MySign(dlam);
      int eFact2;
      if (signDlam != signAzim) {
	eFact2=1;
      }
      else {
	eFact2=0;
      }
      int eFact1=-signAzim;
      double fAzim=eAzim+ThePi*eFact1*eFact2;
      azim(iEta,iXi)=fAzim;
    }
  MyMatrix<double> ANG_rho(eta_rho, xi_rho);
  for (int iEta=1; iEta<eta_v; iEta++)
    for (int iXi=0; iXi<xi_v; iXi++)
      ANG_rho(iEta,iXi)=ThePi/double(2) - azim(iEta-1,iXi);
  for (int iXi=0; iXi<xi_v; iXi++) {
    ANG_rho(0,iXi)=ANG_rho(1,iXi);
    ANG_rho(eta_rho-1,iXi)=ANG_rho(eta_rho-2,iXi);
  }
  return ANG_rho;
}     
/* ------------------------------ */


/* ------------------------------ */
struct T_stat {
  int nbMeas;
  double MaxMeas;
  double MinMeas;
  double MaxModel;
  double MinModel;
  double MeanMeas;
  double MeanModel;
  double MeanError;
  double AbsoluteError;
  double RMSE;
  double CenteredRMSE;
  double Correlation;
  double ScatterIndex;
  double CenteredScatterIndex;
  double Slope;
  std::string strMaxMeas;
  std::string strMinMeas;
  std::string strMaxModel;
  std::string strMinModel;
  std::string strMeanMeas;
  std::string strMeanModel;
  std::string strMeanError;
  std::string strAbsoluteError;
  std::string strRMSE;
  std::string strCenteredRMSE;
  std::string strCorrelation;
  std::string strScatterIndex;
  std::string strCenteredScatterIndex;
  std::string strSlope;
  std::string strNature="ME    AE    RMSE CRMSE  CORR   SCI   CSCI";
  std::string str;
};
/* ------------------------------ */


/* ------------------------------ */
struct PairMM {
  double Meas;
  double Model;
};
/* ------------------------------ */


/* ------------------------------ */
T_stat ComputeStatistics_Pair(std::vector<PairMM> const& eVect)
{
  T_stat eStat;
  int nbMeas=0;
  double SumAbs=0;
  double SumSqr=0;
  double eSum1=0;
  double eSum2=0;
  double eSum11=0;
  double eSum12=0;
  double eSum22=0;
  double MaxMeas=-10^(31);
  double MaxModel=-10^(31);
  double MinMeas=10^(31);
  double MinModel=10^(31);
  for (auto& ePair : eVect) {
    nbMeas++;
    double eMeas=ePair.Meas;
    double eModel=ePair.Model;
    MaxMeas=std::max(MaxMeas, eMeas);
    MaxModel=std::max(MaxModel, eModel);
    MinMeas=std::min(MinMeas, eMeas);
    MinModel=std::min(MinModel, eModel);
    eSum1 =eSum1  + eMeas;
    eSum2 =eSum2  + eModel;
    eSum11=eSum11 + eMeas*eMeas;
    eSum12=eSum12 + eMeas*eModel;
    eSum22=eSum22 + eModel*eModel;
    SumAbs=SumAbs + fabs(eMeas - eModel);
    double eDiff=eMeas-eModel;
    SumSqr=SumSqr + eDiff*eDiff;
  }
  double eME=(eSum2 - eSum1)/double(nbMeas);
  double eRMSE=sqrt(SumSqr / double(nbMeas));
  double eCentRMSE=sqrt(eRMSE*eRMSE - eME*eME);
  double eAE=SumAbs/double(nbMeas);
  double avgSum1=eSum1/double(nbMeas);
  double avgSum2=eSum2/double(nbMeas);
  double avgSum11=eSum11/double(nbMeas);
  double avgSum12=eSum12/double(nbMeas);
  double avgSum22=eSum22/double(nbMeas);
  double eProd11=avgSum11 - avgSum1*avgSum1;
  double eProd12=avgSum12 - avgSum1*avgSum2;
  double eProd22=avgSum22 - avgSum2*avgSum2;
  double TheCorr=eProd12/sqrt(eProd11*eProd22);
  double eScat=eRMSE/avgSum1;
  double eCentScat=eCentRMSE/avgSum1;
  double eSlope=eSum12/eSum11;
  //
  eStat.nbMeas=nbMeas;
  eStat.MaxMeas=MaxMeas;
  eStat.MinMeas=MinMeas;
  eStat.MaxModel=MaxModel;
  eStat.MinModel=MinModel;
  eStat.MeanMeas=avgSum1;
  eStat.MeanModel=avgSum2;
  eStat.MeanError=eME;
  eStat.AbsoluteError=eAE;
  eStat.RMSE=eRMSE;
  eStat.CenteredRMSE=eCentRMSE;
  eStat.Correlation=TheCorr;
  eStat.ScatterIndex=eScat;
  eStat.CenteredScatterIndex=eCentScat;
  eStat.Slope=eSlope;
  //
  eStat.strMaxMeas=DoubleTo4dot2f(MaxMeas);
  eStat.strMinMeas=DoubleTo4dot2f(MinMeas);
  eStat.strMaxModel=DoubleTo4dot2f(MaxModel);
  eStat.strMinModel=DoubleTo4dot2f(MinModel);
  eStat.strMeanMeas=DoubleTo4dot2f(avgSum1);
  eStat.strMeanModel=DoubleTo4dot2f(avgSum2);
  eStat.strMeanError=DoubleTo4dot2f(eME);
  eStat.strAbsoluteError=DoubleTo4dot2f(eAE);
  eStat.strRMSE=DoubleTo4dot2f(eRMSE);
  eStat.strCenteredRMSE=DoubleTo4dot2f(eCentRMSE);
  eStat.strCorrelation=DoubleTo4dot2f(TheCorr);
  eStat.strScatterIndex=DoubleTo4dot2f(eScat);
  eStat.strCenteredScatterIndex=DoubleTo4dot2f(eCentScat);
  eStat.strSlope=DoubleTo4dot2f(eSlope);
  //
  eStat.str=eStat.strMeanError + " " + eStat.strAbsoluteError + " " + eStat.strRMSE + " " + eStat.strCenteredRMSE + " " + eStat.strCorrelation + " " + eStat.strScatterIndex + " " + eStat.strCenteredScatterIndex;
  return eStat;
}
/* ------------------------------ */


/* ------------------------------ */
T_stat ComputeStatistics_MyVector(MyVector<double> const& ListMeas, MyVector<double> const& ListModel)
{
  if (ListMeas.size() != ListModel.size()) {
    std::cerr << "Error in ComputeStatistics_MyVector\n";
    std::cerr << "Discrepancy in number of measurements\n";
    std::cerr << "Please solve the problem\n";
    exit(1);
  }
  std::vector<PairMM> ListPair;
  int nbEnt=ListMeas.size();
  for (int iEnt=0; iEnt<nbEnt; iEnt++) {
    ListPair.push_back({ListMeas(iEnt), ListModel(iEnt)});
  }
  return ComputeStatistics_Pair(ListPair);
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
struct ThreeDimArray {
public:

  ThreeDimArray()
  {
    n1=0;
    n2=0;
    n3=0;
    ListElt=NULL;
  }
  ThreeDimArray(size_t p1, size_t p2, size_t p3)
  {
    if (p1 < 0 || p2 < 0 || p3 < 0)
      {
	std::cerr << "p123=" << p1 << "," << p2 << "," << p3 << "\n";
	std::cerr <<"(zero dimension allowed for things like nullspace)\n";
	exit(1);
      }
    n1=p1;
    n2=p2;
    n3=p3;
    ListElt=new T[n1*n2*n3];
  }
  ThreeDimArray(ThreeDimArray<T> const& eArr)
  {
    std::vector<size_t> LDim=eArr.GetDims();
    n1=LDim[0];
    n2=LDim[1];
    n3=LDim[2];
    size_t len=n1*n2*n3;
    ListElt=new T[n1*n2*n3];
    for (size_t ielt=0; ielt<len; ielt++) {
      T eVal=eArr.getelt(ielt);
      ListElt[ielt]=eVal;
    }
  }
  T getelt(size_t idx) const
  {
    return ListElt[idx];
  }
  ThreeDimArray<T> operator=(ThreeDimArray<T> const& eArr)
  {
    std::vector<size_t> LDim=eArr.GetDims();
    n1=LDim[0];
    n2=LDim[1];
    n3=LDim[2];
    size_t len=n1*n2*n3;
    ListElt=new T[n1*n2*n3];
    for (size_t ielt=0; ielt<len; ielt++) {
      T eVal=eArr.getelt(ielt);
      ListElt[ielt]=eVal;
    }
    return *this;
  }

  ~ThreeDimArray()
  {
    delete [] ListElt;
  }
  // easier parts of the job
  std::vector<size_t> GetDims(void) const
  {
    std::vector<size_t> eRet={n1, n2, n3};
    return eRet;
  }
  void assign(size_t i1, size_t i2, size_t i3, T const& eVal)
  {
    size_t idx=i1 + n1*i2 + n1*n2*i3;
    ListElt[idx]=eVal;
  }
  T get(size_t i1, size_t i2, size_t i3) const
  {
    size_t idx=i1 + n1*i2 + n1*n2*i3;
    if (i1 > n1-1 || i1 < 0 ) {
      std::cerr << "i1=" << i1 << " n1=" << n1 << "\n";
      exit(1);
    }
    if (i2 > n2-1 || i2 < 0 ) {
      std::cerr << "i2=" << i2 << " n2=" << n2 << "\n";
      exit(1);
    }
    if (i3 > n3-1 || i3 < 0 ) {
      std::cerr << "i3=" << i3 << " n3=" << n3 << "\n";
      exit(1);
    }
    return ListElt[idx];
  }
  MyMatrix<T> DimensionExtraction(size_t const& iDim, size_t const& eDim) const
  {
    if (iDim == 0) {
      MyMatrix<T> eMat(n2, n3);
      for (size_t i2=0; i2<n2; i2++)
	for (size_t i3=0; i3<n3; i3++) {
	  size_t idx=eDim + n1*i2 + n1*n2*i3;
	  eMat(i2, i3)=ListElt[idx];
	}
      return eMat;
    }
    if (iDim == 1) {
      MyMatrix<T> eMat(n1, n3);
      for (size_t i1=0; i1<n1; i1++)
	for (size_t i3=0; i3<n3; i3++) {
	  size_t idx=i1 + n1*eDim + n1*n2*i3;
	  eMat(i1, i3)=ListElt[idx];
	}
      return eMat;
    }
    if (iDim == 2) {
      MyMatrix<T> eMat(n1, n2);
      for (size_t i1=0; i1<n1; i1++)
	for (size_t i2=0; i2<n2; i2++) {
	  size_t idx=i1 + n1*i2 + n1*n2*eDim;
	  eMat(i1, i1)=ListElt[idx];
	}
      return eMat;
    }
    std::cerr << "Wrong input in ThreeDimArray\n";
    std::cerr << "iDim=" << iDim << "\n";
    std::cerr << "Allowed values: 0, 1, 2\n";
    exit(1);
  }
private:
  size_t n1, n2, n3;
  T *ListElt;
};
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
T ScalarProduct(MyVector<T> const& V1, MyVector<T> const & V2)
{
  if (V1.size() != V2.size()) {
    std::cerr << "Vectors of wrong sizes\n";
    exit(1);
  }
  size_t siz=V1.size();
  T eSum=0;
  for (size_t i=0; i<siz; i++) {
    T eVal1=V1(i);
    T eVal2=V2(i);
    eSum += eVal1*eVal2;
  }
  return eSum;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyMatrix<T> ZeroMatrix(int const& nbRow, int const& nbCol)
{
  MyMatrix<T> retMat(nbRow, nbCol);
  T eZero;
  eZero=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      retMat(iRow, iCol)=eZero;
  return retMat;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyVector<T> ZeroVector(int const& nbRow)
{
  MyVector<T> retVect(nbRow);
  T eZero;
  eZero=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    retVect(iRow)=eZero;
  return retVect;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
void TMat_Copy(MyMatrix<T> const&eMatI, MyMatrix<T> &eMatO)
{
  int nbRowI, nbColI, nbRowO, nbColO;
  nbRowI=eMatI.rows();
  nbRowO=eMatO.rows();
  nbColI=eMatI.cols();
  nbColO=eMatO.cols();
  if (nbRowI != nbRowO || nbColI != nbColO) {
    std::cerr << "Error in the input\n";
    exit(1);
  }
  for (int iRow=0; iRow<nbRowI; iRow++)
    for (int iCol=0; iCol<nbColI; iCol++) {
      T eVal=eMatI(iRow, iCol);
      eMatO(iRow, iCol)=eVal;
    }
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
void ZeroAssignation(MyMatrix<T> &TheMat)
{
  int nbRow=TheMat.rows();
  int nbCol=TheMat.cols();
  T eVal;
  eVal=0;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      TheMat(iRow, iCol)=eVal;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyMatrix<T> TransposedMat(MyMatrix<T> const&TheMat)
{
  int nbCol=TheMat.cols();
  int nbRow=TheMat.rows();
  MyMatrix<T> TheTrans(nbCol, nbRow);
  for (int iCol=0; iCol<nbCol; iCol++)
    for (int iRow=0; iRow<nbRow; iRow++) {
      T eVal=TheMat(iRow, iCol);
      TheTrans(iCol, iRow)=eVal;
    }
  return TheTrans;
}
/* ------------------------------ */


/* ------------------------------ */
/* We allocate M3 and compute the product M3=M1 M2 */
template<typename T>
MyMatrix<T> TMat_Product(MyMatrix<T> const&M1, MyMatrix<T> const&M2)
{
  int nbCol1=M1.cols();
  int nbRow1=M1.rows();
  int nbCol2=M2.cols();
  int nbRow2=M2.rows();
  if (nbCol1 != nbRow2) {
    std::cerr << "Error in matrix sizes\n";
    exit(1);
  }
  MyMatrix<T> TheProd(nbRow1, nbCol2);
  for (int iCol=0; iCol<nbCol2; iCol++)
    for (int iRow=0; iRow<nbRow1; iRow++) {
      T eSum=0;
      for (int i=0; i<nbCol1; i++) {
	T eVal1=M1(iRow, i);
	T eVal2=M2(i, iCol);
	eSum += eVal1*eVal2;
      }
      TheProd(iRow, iCol)=eSum;
    }
  return TheProd;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyVector<T> VectorMatrix(MyVector<T> const& eVect, MyMatrix<T> const& eMat)
{
  int nbCol=eMat.cols();
  int nbRow=eMat.rows();
  int n=eVect.size();
  if (n != nbRow) {
    std::cerr << "n should be equal to nbRow\n";
    exit(1);
  }
  MyVector<T> rVect(nbCol);
  for (int iCol=0; iCol<nbCol; iCol++) {
    T eSum=0;
    for (int iRow=0; iRow<nbRow; iRow++) {
      T eVal=eMat(iRow, iCol);
      T fVal=eVect(iRow);
      eSum += eVal*fVal;
    }
    rVect(iCol)=eSum;
  }
  return rVect;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
void SwapValues(T& val1, T& val2)
{
  T prov;
  prov=val1;
  val1=val2;
  val2=prov;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
struct Inverse_exception {
  std::string errmsg;
  T pivot;
};
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
void TMat_Inverse_destroy(MyMatrix<T> &Input, MyMatrix<T> &Output)
{
  int nbRow, nbCol;
  int iCol, iRow, WeFound;
  int iRowFound;
  int iColB;
  nbRow=Input.rows();
  nbCol=Input.cols();
  T prov1, prov2, eVal;
  if (nbRow != nbCol) {
    std::cerr << "Error on nbRow, nbCol in TMat_Inverse_destroy";
    exit(1);
  }
  for (iRow=0; iRow<nbRow; iRow++)
    for (iCol=0; iCol<nbRow; iCol++)
      {
	if (iRow == iCol)
	  prov1=1;
	else
	  prov1=0;
	Output(iRow,iCol)=prov1;
      }
  iRowFound=-400;
  for (iCol=0; iCol<nbCol; iCol++)
    {
      WeFound=0;
      for (iRow=iCol; iRow<nbRow; iRow++)
	if (WeFound == 0)
	  {
	    eVal=Input(iRow,iCol);
	    if (eVal != 0)
	      {
		WeFound=1;
		iRowFound=iRow;
		prov1=1/eVal;
	      }
	  }
      if (WeFound == 0) {
	Inverse_exception<T> eExcept;
	eExcept.errmsg="Error in matrix inversion";
	eExcept.pivot=0;
	throw eExcept;
      }
      for (iColB=0; iColB<nbCol; iColB++)
	{
	  eVal=prov1*Input(iRowFound,iColB);
	  Input(iRowFound,iColB)=eVal;
	  eVal=prov1*Output(iRowFound,iColB);
	  Output(iRowFound,iColB)=eVal;
	}
      for (iRow=0; iRow<nbRow; iRow++)
	if (iRow != iRowFound) {
	  prov2=Input(iRow, iCol);
	  for (iColB=0; iColB<nbCol; iColB++) {
	    prov1=prov2*Input(iRowFound,iColB);
	    eVal=Input(iRow,iColB) - prov1;
	    Input(iRow, iColB)=eVal;
	    //
	    prov1=prov2*Output(iRowFound,iColB);
	    eVal=Output(iRow,iColB) - prov1;
	    Output(iRow,iColB)=eVal;
	  }
	}
      if (iRowFound != iCol) {
	for (iColB=0; iColB<nbCol; iColB++) {
	  prov1=Input(iRowFound, iColB);
	  prov2=Input(iCol, iColB);
	  SwapValues(prov1, prov2);
	  Input(iRowFound, iColB)=prov1;
	  Input(iCol     , iColB)=prov2;
	  //
	  prov1=Output(iRowFound, iColB);
	  prov2=Output(iCol, iColB);
	  SwapValues(prov1, prov2);
	  Output(iRowFound, iColB)=prov1;
	  Output(iCol     , iColB)=prov2;
	}
      }
    }
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyMatrix<T> Inverse(MyMatrix<T> const&Input)
{
  int nbRow, nbCol;
  nbRow=Input.rows();
  nbCol=Input.cols();
  MyMatrix<T> provMat(nbRow, nbCol);
  TMat_Copy(Input, provMat);
  MyMatrix<T> Output(nbRow, nbRow);
  TMat_Inverse_destroy(provMat, Output);
  return Output;
}
/* ------------------------------ */


/* ------------------------------ */
/* This function is for rank calculation.
   Of course, it can be used for many other purpose:
   1> Selecting specific sets of rows and columns for reduction 
   2> Computing a set of generators of the zero set.
   Initial matrix is of the form
   Input (nbRow, nbCol)
   The zero matrix is of the form
   NSP (dimKer, nbCol)
*/
template<typename T>
struct SelectionRowCol {
  int TheRank;
  MyMatrix<T> NSP;
  std::vector<int> ListColSelect;
  std::vector<int> ListRowSelect;
};
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
SelectionRowCol<T> TMat_SelectRowCol(MyMatrix<T> const&Input)
{
  int nbRow, nbCol;
  T eVal, eVal1, eVal2, eVal3;
  int iRank;
  int sizMat, nbVect, iRow, iCol;
  int eCol, IsFinished;
  int nbVectZero, maxRank, eRank, FirstNonZeroCol;
  nbRow=Input.rows();
  nbCol=Input.cols();
  maxRank=nbRow;
  if (nbCol < maxRank)
    maxRank=nbCol;
  sizMat=maxRank+1;
  MyMatrix<T> provMat(sizMat, nbCol);
  std::vector<int> ListColSelect;
  std::vector<int> ListRowSelect;
  std::vector<int> ListColSelect01(nbCol);
  for (iCol=0; iCol<nbCol; iCol++)
    ListColSelect01[iCol]=0;
  eRank=0;
  for (iRow=0; iRow<nbRow; iRow++) {
    for (iCol=0; iCol<nbCol; iCol++) {
      eVal=Input(iRow, iCol);
      provMat(eRank, iCol)=eVal;
    }
    for (iRank=0; iRank<eRank; iRank++) {
      eCol=ListColSelect[iRank];
      eVal1=provMat(eRank, eCol);
      for (iCol=0; iCol<nbCol; iCol++) {
	eVal2=provMat(iRank, iCol);
	eVal3=provMat(eRank, iCol);
	eVal=eVal1*eVal2;
	eVal3=eVal3 - eVal;
	  provMat(eRank, iCol)=eVal3;
      }
    }
    IsFinished=1;
    FirstNonZeroCol=-1;
    for (iCol=0; iCol<nbCol; iCol++)
      if (IsFinished == 1) {
	eVal=provMat(eRank, iCol);
	if (eVal != 0) {
	  FirstNonZeroCol=iCol;
	  IsFinished=0;
	}
      }
    if (IsFinished == 0) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListRowSelect.push_back(iRow);
      ListColSelect01[FirstNonZeroCol]=1;
      eVal=provMat(eRank, FirstNonZeroCol);
      eVal2=1/eVal;
      for (iCol=0; iCol<nbCol; iCol++) {
	eVal=provMat(eRank, iCol);
	eVal=eVal*eVal2;
	provMat(eRank, iCol)=eVal;
      }
      for (iRank=0; iRank<eRank; iRank++) {
	eVal1=provMat(iRank, FirstNonZeroCol);
	for (iCol=0; iCol<nbCol; iCol++) {
	  eVal2=provMat(iRank, iCol);
	  eVal3=provMat(eRank, iCol);
	  eVal=eVal1*eVal3;
	  eVal2=eVal2 - eVal;
	  provMat(iRank, iCol)=eVal2;
	}
      }
      eRank++;
    }
  }
  nbVectZero=nbCol-eRank;
  MyMatrix<T> NSP(nbVectZero, nbCol);
  //  std::cerr << "nbVectZero=" << nbVectZero << " nbCol=" << nbCol << "\n";
  ZeroAssignation(NSP);
  nbVect=0;
  for (iCol=0; iCol<nbCol; iCol++)
    if (ListColSelect01[iCol] == 0)
      {
	eVal=1;
	NSP(nbVect, iCol)=eVal;
	for (iRank=0; iRank<eRank; iRank++)
	  {
	    eCol=ListColSelect[iRank];
	    eVal=provMat(iRank, iCol);
	    eVal=-eVal;
	    NSP(nbVect, eCol)=eVal;
	  }
	nbVect++;
      }
  SelectionRowCol<T> retSelect{eRank, NSP, ListColSelect, ListRowSelect};
  return retSelect;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyVector<T> TVec_AplusBC(MyVector<T> const& V1, MyVector<T> const& V2, T const& alpha)
{
  int n=V1.size();
  MyVector<T> VRet=MyVector<T>(n);
  for (int i=0; i<n; i++)
    {
      T eVal=V1(i) + alpha*V2(i);
      VRet(i)=eVal;
    }
  return VRet;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
void TVec_AddBCtoA(MyVector<T> &V1, MyVector<T> const& V2, T const& alpha)
{
  T eVal;
  int n=V1.size();
  for (int i=0; i<n; i++)
    {
      eVal=V1(i) + alpha*V2(i);
      V1(i)=eVal;
    }
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyMatrix<T> TMat_SelectRow(MyMatrix<T> const&TheMat, std::vector<int> const& eList)
{
  int nbRowRed=eList.size();
  int nbCol=TheMat.cols();
  MyMatrix<T> TheProv(nbRowRed, nbCol);
  for (int iRow=0; iRow<nbRowRed; iRow++) {
    int jRow=eList[iRow];
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=TheMat(jRow, iCol);
      TheProv(iRow, iCol)=eVal;
    }
  }
  return TheProv;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyMatrix<T> TMat_SelectRow(MyMatrix<T> const&TheMat, Face const& eList)
{
  int nbRowRed=eList.count();
  int nbCol=TheMat.cols();
  MyMatrix<T> TheProv(nbRowRed, nbCol);
  int jRow=eList.find_first();
  for (int iRow=0; iRow<nbRowRed; iRow++) {
    for (int iCol=0; iCol<nbCol; iCol++) {
      T eVal=TheMat(jRow, iCol);
      TheProv(iRow, iCol)=eVal;
    }
    jRow=eList.find_next(jRow);
  }
  return TheProv;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyMatrix<T> TMat_SelectCol(MyMatrix<T> const& TheMat, std::vector<int> const& eList)
{
  int nbRow=TheMat.rows();
  int nbColRed=eList.size();
  MyMatrix<T> TheProv(nbRow, nbColRed);
  for (int iCol=0; iCol<nbColRed; iCol++) {
    int jCol=eList[iCol];
    for (int iRow=0; iRow<nbRow; iRow++)
      TheProv(iRow, iCol)=TheMat(iRow, jCol);
  }
  return TheProv;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyVector<T> TVec_SelectCol(MyVector<T> const& TheV, std::vector<int> const& eList)
{
  int nbColRed=eList.size();
  MyVector<T> TheProv(nbColRed);
  for (int iCol=0; iCol<nbColRed; iCol++) {
    int jCol=eList[iCol];
    TheProv(iCol)=TheV(jCol);
  }
  return TheProv;
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
struct SolMatResult {
  bool result;
  MyVector<T> eSol;
};
/* ------------------------------ */


/* ------------------------------ */
template<typename T> 
struct resultFT {
  bool test;
  MyMatrix<T> eMat;
};
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
resultFT<T> FindTransformationGeneral(MyMatrix<T> const& EXT1, MyMatrix<T> const& EXT2, std::vector<int> const& ePerm)
{
  if (EXT1.cols() != EXT2.cols() )
    return {false, {}};
  if (EXT1.rows() != EXT2.rows() )
    return {false, {}};
  int nbCol=EXT1.cols();
  int nbRow=EXT1.rows();
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(EXT1);
  int eRank=eSelect.TheRank;
  if (eRank != nbCol)
    return {false, {}};
  MyMatrix<T> eMat1(nbCol, nbCol);
  MyMatrix<T> eMat2(nbCol, nbCol);
  for (int iRow=0; iRow<nbCol; iRow++) {
    int iRow1=eSelect.ListRowSelect[iRow];
    int iRow2=ePerm[iRow1];
    eMat1.row(iRow)=EXT1.row(iRow1);
    eMat2.row(iRow)=EXT2.row(iRow2);
  }
  MyMatrix<T> eMat1inv=Inverse(eMat1);
  MyMatrix<T> RetMat=TMat_Product(eMat1inv, eMat2);
  MyMatrix<T> CheckMat=TMat_Product(EXT1, RetMat);
  for (int iRow1=0; iRow1<nbRow; iRow1++) {
    int iRow2=ePerm[iRow1];
    for (int iCol=0; iCol<nbCol; iCol++)
      if (CheckMat(iRow1,iCol) != EXT2(iRow2,iCol))
	return {false, {}};
  }
  return {true, RetMat};
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
void WriteStdVector(std::ostream &os, std::vector<T> const & TheVec)
{
  T eVal;
  int i, n;
  n=TheVec.size();
  for (i=0; i<n; i++) {
    eVal=TheVec[i];
    os << " " << eVal;
  }
  os << "\n";
}
/* ------------------------------ */


/* ------------------------------ */
template<typename T>
MyVector<T> CanonicalizeVector(MyVector<T> const& V)
{
  int n=V.size();
  T TheMin=0;
  int iSelect=-1;
  for (int i=0; i<n; i++) {
    T eVal=V(i);
    if (eVal != 0) {
      if (iSelect == -1) {
	TheMin=T_abs(eVal);
	iSelect=i;
      }
      else {
	T eAbs=T_abs(eVal);
	if (eAbs < TheMin) {
	  TheMin=eAbs;
	  iSelect=i;
	}
      }
    }
  }
  if (iSelect == -1)
    return V;
  MyVector<T> Vret(n);
  T eQuot=1/V(iSelect);
  for (int i=0; i<n; i++)
    Vret(i)=V(i)*eQuot;
  return Vret;
}
/* ------------------------------ */


/* ------------------------------ */
/* return true if V1 < V2 according to lexicographic order */
template<typename T>
bool IsLower(MyVector<T> const& V1, MyVector<T> const& V2)
{
  int n=V1.size();
  for (int i=0; i<n; i++) {
    if (V1(i) != V2(i)) {
      if (V1(i) < V2(i)) {
	return true;
      }
      else {
	return false;
      }
    }
  }
  return false;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<int> GetEdgeSet(MyMatrix<int> const& INE, int nbNode)
{
  int ListDegree[nbNode];
  for (int iNode=0; iNode<nbNode; iNode++) {
    ListDegree[iNode]=0;
  }
  int mne=INE.rows();
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int ip=INE(ie,i);
      ListDegree[ip]+=2;
    }
  int MaxDeg=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    int eDeg=ListDegree[iNode];
    if (eDeg > MaxDeg)
      MaxDeg=eDeg;
  }
  for (int iNode=0; iNode<nbNode; iNode++) {
    ListDegree[iNode]=0;
  }
  int ListAdjacency[nbNode][MaxDeg];
  for (int ie=0; ie<mne; ie++) {
    int i1=INE(ie,0);
    int i2=INE(ie,1);
    int i3=INE(ie,2);
    int eDeg1=ListDegree[i1];
    int eDeg2=ListDegree[i2];
    int eDeg3=ListDegree[i3];
    ListAdjacency[i1][eDeg1]=i2;
    ListAdjacency[i1][eDeg1+1]=i3;
    ListAdjacency[i2][eDeg2]=i1;
    ListAdjacency[i2][eDeg2+1]=i3;
    ListAdjacency[i3][eDeg3]=i1;
    ListAdjacency[i3][eDeg3+1]=i2;
    ListDegree[i1] = eDeg1 + 2;
    ListDegree[i2] = eDeg2 + 2;
    ListDegree[i3] = eDeg3 + 2;
  }
  int nbEdge=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    std::set<int> eSet;
    int eDeg=ListDegree[iNode];
    for (int iAdj=0; iAdj<eDeg; iAdj++) {
      int eAdj=ListAdjacency[iNode][iAdj];
      if (eAdj > iNode)
	eSet.insert(eAdj);
    }
    nbEdge+=eSet.size();
  }
  MyMatrix<int> ListEdges(nbEdge,2);
  int iEdge=0;
  for (int iNode=0; iNode<nbNode; iNode++) {
    std::set<int> eSet;
    int eDeg=ListDegree[iNode];
    for (int iAdj=0; iAdj<eDeg; iAdj++) {
      int eAdj=ListAdjacency[iNode][iAdj];
      if (eAdj > iNode)
	eSet.insert(eAdj);
    }
    for (auto eAdj : eSet) {
      ListEdges(iEdge,0)=iNode;
      ListEdges(iEdge,1)=eAdj;
      iEdge++;
    }
  }
  return ListEdges;
}
/* ------------------------------ */


/* ------------------------------ */
double DATE2JD(std::vector<int> const& Date)
{
  double eJDbase, eFracDay, eJD;
  int year, month, day, hour, min, sec;
  year=Date[0];
  month=Date[1];
  day=Date[2];
  hour=Date[3];
  min=Date[4];
  sec=Date[5];
  int a, y, m;
  a = floor((double(14) - double(month))/double(12));
  y = year + 4800 - a;
  m = month + 12*a - 3;
  // For a date in the Gregorian calendar:
  eJDbase = double(day) 
    + double(floor((double(153)*double(m) + double(2))/double(5)))
    + double(y)*double(365)                                       
    + double(floor(double(y)/double(4)))                           
    - double(floor(double(y)/double(100)))                          
    + double(floor(double(y)/double(400))) - double(32045);
  eFracDay=(double(sec) +                                         
	    double(60)*double(min) +                               
	    double(3600)*(double(hour) - double(12))               
	    )/double(86400);
  eJD=eJDbase + eFracDay;
  return eJD;
}
/* ------------------------------ */


/* ------------------------------ */
double DATE_ConvertSix2mjd(std::vector<int> const& eDate)
{
  double eJD1=DATE2JD(eDate);
  double eJD2=DATE2JD({1858, 11, 17, 0, 0, 0});
  double eMJD=eJD1-eJD2;
  return eMJD;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<int> DATE_ConvertString2six(std::string const& eTimeStr)
{
  std::string eYear, eMonth, eDay, eHour, eMin, eSec;
  eYear=eTimeStr.substr(0,4);
  eMonth=eTimeStr.substr(4,2);
  eDay=eTimeStr.substr(6,2);
  eHour=eTimeStr.substr(9,2);
  eMin=eTimeStr.substr(11,2);
  eSec=eTimeStr.substr(13,2);

  int year, month, day, hour, min, sec;
  std::istringstream(eYear) >> year;
  std::istringstream(eMonth) >> month;
  std::istringstream(eDay) >> day;
  std::istringstream(eHour) >> hour;
  std::istringstream(eMin) >> min;
  std::istringstream(eSec) >> sec;
  std::vector<int> Date={year, month, day, hour, min, sec};
  return Date;
}
/* ------------------------------ */


/* ------------------------------ */
std::string DATE_ConvertSix2string(std::vector<int> const& Date)
{
  int year, month, day, hour, min, sec;
  year=Date[0];
  month=Date[1];
  day=Date[2];
  hour=Date[3];
  min=Date[4];
  sec=Date[5];
  std::string eTimeStr=StringNumber(year, 4) + 
    StringNumber(month, 2) + 
    StringNumber(day, 2) + "." + 
    StringNumber(hour, 2) + 
    StringNumber(min, 2) + 
    StringNumber(sec, 2);
  return eTimeStr;
}
/* ------------------------------ */


/* ------------------------------ */
std::string DATE_ConvertSix2mystringPres(std::vector<int> const& Date)
{
  try {
    int year, month, day, hour, min, sec;
    year=Date[0];
    month=Date[1];
    day=Date[2];
    hour=Date[3];
    min=Date[4];
    sec=Date[5];
    std::string eTimeStr=StringNumber(year, 4) + "-" + 
      StringNumber(month, 2) + "-" +
      StringNumber(day, 2) + " " + 
      StringNumber(hour, 2) + ":" + StringNumber(min, 2) + ":" + StringNumber(sec, 2);
    return eTimeStr;
  }
  catch (std::string & eStr) {
    std::stringstream s;
    s << "Error in DATE_ConvertSix2mystringFile\n";
    s << "Date.size()=" << Date.size() << "\n";
    s << "Date=";
    WriteStdVector(s, Date);
    s << "-----------------------------------------\n";
    s << "exception eStr=\n";
    s << eStr;
    std::string vStr(s.str());
    std::cerr << vStr;
    exit(1);
  }
}
/* ------------------------------ */


/* ------------------------------ */
std::string DATE_ConvertSix2mystringFile(std::vector<int> const& Date)
{
  try {
    int year, month, day, hour, min, sec;
    year=Date[0];
    month=Date[1];
    day=Date[2];
    hour=Date[3];
    min=Date[4];
    sec=Date[5];
    std::string eTimeStr=StringNumber(year, 4) + 
      StringNumber(month, 2) + StringNumber(day, 2) + "_" +
      StringNumber(hour, 2) + StringNumber(min, 2) + StringNumber(sec, 2);
    return eTimeStr;
  }
  catch (std::string & eStr) {
    std::stringstream s;
    s << "Error in DATE_ConvertSix2mystringFile\n";
    s << "Date.size()=" << Date.size() << "\n";
    s << "Date=";
    WriteStdVector(s, Date);
    s << "-----------------------------------------\n";
    s << "exception eStr=\n";
    s << eStr;
    std::string vStr(s.str());
    std::cerr << vStr;
    exit(1);
  }
}
/* ------------------------------ */


/* ------------------------------ */
int MONTH_LEN(int const& year, int const& month)
{
  if (month == 1 || month == 3 || month == 5 || month == 7 || month == 8 || month == 10 || month == 12)
    return 31;
  if (month == 4 || month == 6 || month == 9 || month == 11)
    return 30;
  if (month == 2) {
    int res4=year % 4;
    int res100=year % 100;
    int res400=year % 400;
    if (res4 != 0) {
      return 28;
    }
    else {
      if (res100 != 0) {
	return 29;
      }
      else {
	if (res400 != 0) {
	  return 28;
	}
	else {
	  return 29;
	}
      }
    }
  }
  std::cerr << "Error happened in LEN_MONTH\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<int> JD2DATE(double const& eJD)
{
  int year, month, day, hour, min, sec;
  int ijd, a, b, c, d, e, m;
  double fjd, second;
  ijd = floor(eJD + 0.5);
  //
  a = ijd + 32044;
  b = floor((double(4)*double(a) + double(3)) / double(146097));
  c = a - floor((double(b) * double(146097)) / double(4));
  //
  d = floor((double(4)*double(c) + double(3)) / double(1461));
  e = c - floor((double(1461)*double(d)) / double(4));
  m = floor((double(5) * double(e) + double(2)) / double(153));
  //
  day   = e - floor((double(153) * double(m) + double(2)) / double(5)) + 1;
  month = m + 3 - 12 * floor(double(m) / double(10));
  year  = b * 100 + d - 4800 + floor(double(m) / double(10));
  //
  fjd    = eJD - double(ijd) + 0.5;
  second = double(86400) * fjd;
  hour   = floor(second/double(3600));
  second = second - double(3600)*double(hour);
  min    = floor(second/double(60));
  sec    = floor(second - double(60)*min);
  // Now renormalizing
  int secNear=round(second - double(60)*min);
  if (secNear == 60) {
    sec=0;
    min=min+1;
  }
  if (min == 60) {
    min=0;
    hour=hour+1;
  }
  if (hour == 24) {
    hour=0;
    day=day+1;
  }
  int lenmonth=MONTH_LEN(year, month);
  if (day == lenmonth+1) {
    day=1;
    month=month+1;
  }
  if (month == 13) {
    month=1;
    year=year+1;
  }
  std::vector<int> Date={year, month, day, hour, min, sec};
  return Date;
}
/* ------------------------------ */


/* ------------------------------ */
double CT2MJD(std::string const& STIME)
{
  double XMJD;
  std::vector<int> eDate=DATE_ConvertString2six(STIME);
  XMJD=DATE_ConvertSix2mjd(eDate);
  return XMJD;
}
/* ------------------------------ */


/* ------------------------------ */
std::string MJD2CT(double const& XMJD)
{
  std::string STIME;
  double XMJD_1858, eMJD;
  XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  STIME=DATE_ConvertSix2string(eDate);
  return STIME;
}
/* ------------------------------ */


/* ------------------------------ */
std::string DATE_ConvertMjd2mystringPres(double const& XMJD)
{
  std::string STIME;
  double XMJD_1858, eMJD;
  XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  STIME=DATE_ConvertSix2mystringPres(eDate);
  return STIME;
}
/* ------------------------------ */


/* ------------------------------ */
std::string DATE_ConvertMjd2mystringFile(double const& XMJD)
{
  std::string STIME;
  double XMJD_1858, eMJD;
  XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
  eMJD = XMJD + XMJD_1858;
  std::vector<int> eDate=JD2DATE(eMJD);
  STIME=DATE_ConvertSix2mystringFile(eDate);
  return STIME;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<double> GetInterval(double const& FirstTime, double const& LastTime, double const& DeltaInterval)
{
  double eTime=FirstTime;
  std::vector<double> ListTime;
  double tolDay=double(1)/double(10000);
  //  std::cerr.width(15);
  //  std::cerr << "DeltaInterval=" << DeltaInterval << "\n";
  while(1) {
    ListTime.push_back(eTime);
    eTime=eTime + DeltaInterval;
    if (eTime > LastTime + tolDay)
      return ListTime;
  }
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<double> GetInterval(std::string const& BEGTC, std::string const& ENDTC, double const& eInterval, std::string const& UNITC)
{
  int IsDone=0;
  double eMult;
  if (UNITC == "DAY") {
    IsDone=1;
    eMult=double(1);
  }
  if (UNITC == "HOUR") {
    IsDone=1;
    eMult=double(1)/double(24);
  }
  if (UNITC == "MIN") {
    IsDone=1;
    eMult=double(1)/double(1440);
  }
  if (UNITC == "SEC") {
    IsDone=1;
    eMult=double(1)/double(86400);
  }
  if (IsDone == 0) {
    std::cerr << "UNITC has not been found\n";
    std::cerr << "Allowed: DAY, HOUR, MIN, SEC\n";
    std::cerr << "UNITC=" << UNITC << "\n";
    exit(1);
  }
  double DeltaInterval=eInterval*eMult;
  double FirstTime=CT2MJD(BEGTC);
  double LastTime=CT2MJD(ENDTC);
  return GetInterval(FirstTime, LastTime, DeltaInterval);
}
/* ------------------------------ */


/* ------------------------------ */
struct InterpInfo {
  int iTimeLow, iTimeUpp;
  double alphaLow, alphaUpp;
  bool UseSingleEntry;
};
/* ------------------------------ */


/* ------------------------------ */
InterpInfo GetTimeInterpolationInfo(std::vector<double> const& LTime, double const& eTimeDay)
{
  InterpInfo eInterpInfo;
  double tolDay=double(1)/double(100000);
  int nbTime=LTime.size();
  for (int iTime=0; iTime<nbTime; iTime++) {
    double eDist=fabs(LTime[iTime] - eTimeDay);
    if (eDist < tolDay) {
      //      std::cerr << "eDist=" << eDist << "\n";
      //      std::cerr << "eTimeDay=" << eTimeDay << "\n";
      //      std::cerr << "LTime[iTime]=" << LTime[iTime] << "\n";
      //      std::cerr << "iTime=" << iTime << "\n";
      eInterpInfo.UseSingleEntry=true;
      eInterpInfo.iTimeLow=iTime;
      return eInterpInfo;
    }
  }
  eInterpInfo.UseSingleEntry=false;
  if (eTimeDay < LTime[0] - tolDay) {
    std::cerr << "The asked entry is before the first time\n";
    std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
    std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
    std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
    exit(1);
  }
  if (eTimeDay > LTime[nbTime-1] + tolDay) {
    std::cerr << "The asked entry is after the last time\n";
    std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
    std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
    std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
    exit(1);
  }
  for (int iTimeUpp=1; iTimeUpp<nbTime; iTimeUpp++) {
    int iTimeLow=iTimeUpp-1;
    double eTimeUpp=LTime[iTimeUpp];
    double eTimeLow=LTime[iTimeLow];
    double alphaLow=(eTimeDay - eTimeUpp)/(eTimeLow - eTimeUpp);
    double alphaUpp=(eTimeLow - eTimeDay)/(eTimeLow - eTimeUpp);
    if (alphaLow >= 0 && alphaUpp >= 0) {
      eInterpInfo.iTimeLow=iTimeLow;
      eInterpInfo.iTimeUpp=iTimeUpp;
      eInterpInfo.alphaLow=alphaLow;
      eInterpInfo.alphaUpp=alphaUpp;
      return eInterpInfo;
    }
  }
  std::cerr << "Failed to find matching record\n";
  std::cerr << "Please debug\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<int> GetIntervalListITime(std::vector<double> const& LTime, double const& eTimeDay, double const& TimeFrameDay)
{
  double epsilon=0.0001;
  std::vector<int> ListRelITime;
  double eTimeLow=eTimeDay - epsilon;
  double eTimeUpp=eTimeDay + TimeFrameDay - epsilon;
  int nbTime=LTime.size();
  for (int iTime=0; iTime<nbTime; iTime++) {
    double eTime=LTime[iTime];
    if (eTime > eTimeLow && eTime < eTimeUpp)
      ListRelITime.push_back(iTime);
  }
  return ListRelITime;
}
/* ------------------------------ */


/* ------------------------------ */
struct PermanentInfoDrawing {
  std::string Extension;
  std::string PicPrefix;
  std::string PrefixTemp;
  bool KeepNC_NCL;
  std::string eDir;
  int NPROC;
};
/* ------------------------------ */


/* ------------------------------ */
std::string PLOT_CreatePrefixTemp(FullNamelist const& eFull)
{
  std::map<std::string, SingleBlock> ListBlock=eFull.ListBlock;
  SingleBlock eBlPROC=ListBlock.at("PROC");
  bool KeepNC_NCL=eBlPROC.ListBoolValues.at("KeepNC_NCL");
  std::string eRand=random_string(20);
  std::string Nature=eBlPROC.ListStringValues.at("__NaturePlot");
  std::string PrefixTemp="/tmp/PLOT_" + Nature + "_" + eRand + "/";
  if (KeepNC_NCL)
    std::cerr << "PrefixTemp = " << PrefixTemp << "\n";
  CreateDirectory(PrefixTemp);
  return PrefixTemp;
}
/* ------------------------------ */


/* ------------------------------ */
PermanentInfoDrawing GET_PERMANENT_INFO(FullNamelist const& eFull)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::string PicPrefix=eBlPROC.ListStringValues.at("PicPrefix");
  std::string Extension=eBlPROC.ListStringValues.at("Extension");
  bool KeepNC_NCL=eBlPROC.ListBoolValues.at("KeepNC_NCL");
  int NPROC=eBlPROC.ListIntValues.at("NPROC");
  std::string PrefixTemp=PLOT_CreatePrefixTemp(eFull);
  if (eBlPROC.ListBoolValues.at("FirstCleanDirectory")) {
    RemoveFileSpecificExtension(PicPrefix, Extension);
    /*
    std::string order="rm -f " + PicPrefix + "*" + Extension;
    int iret=system(order.c_str() );
    if (iret == -1) {
      std::cerr << "order = " << order << "\n";
      std::cerr << "Unable to launch the process\n";
      exit(1);
      }*/
  }
  std::string eDir=FILE_GetAbsoluteDirectory(PicPrefix);
  CreateDirectory(eDir);
  //
  PermanentInfoDrawing ePerm;
  ePerm.PrefixTemp=PrefixTemp;
  ePerm.eDir=eDir;
  ePerm.Extension=Extension;
  ePerm.PicPrefix=PicPrefix;
  ePerm.KeepNC_NCL=KeepNC_NCL;
  ePerm.NPROC=NPROC;
  return ePerm;
}
/* ------------------------------ */


/* ------------------------------ */
struct RecVar {
  double eTimeDay;
  int iTime;
  std::string strPres;
  std::string strFile;
  std::string strAll;
  //
  std::string VarName1;
  std::string VarName2;
  double minval;
  double maxval;
  double mindiff;
  double maxdiff;
  std::string Unit;
  MyMatrix<double> U;
  MyMatrix<double> V;
  MyMatrix<double> F;
  std::string VarNature;
  std::string nameU, nameV;
};
/* ------------------------------ */


/* ------------------------------ */
struct CoordGridArrayFD {
  int eta, xi;
  MyMatrix<int> MSK;
  MyMatrix<double> LON, LAT, DEP, ANG;
  int nbWet;
  bool HaveDEP;
  std::vector<int> Idx, Jdx;
};
/* ------------------------------ */


/* ------------------------------ */
struct GridArray {
  std::string ModelName;
  int IsFE;
  int IsLonLat;
  CoordGridArrayFD GrdArrRho, GrdArrU, GrdArrV, GrdArrPsi;
  MyMatrix<int> INE;
  bool L_IndexSelect;
  std::vector<int> I_IndexSelect;
};
/* ------------------------------ */


/* ------------------------------ */
struct ArrayHistory {
  std::string KindArchive;
  int nbFile, nbTime;
  double FirstTime, LastTime;
  std::string FirstTimeStr, LastTimeStr;
  std::vector<std::string> ListFileNames;
  std::vector<std::vector<std::string> > ListListFileNames;
  std::vector<int> ListIFile;
  std::vector<int> ListIRec;
  std::vector<double> ListTime;
  bool AppendVarName;
};
/* ------------------------------ */


/* ------------------------------ */
struct TotalArrGetData {
  GridArray GrdArr;
  ArrayHistory eArr;
};
/* ------------------------------ */


/* ------------------------------ */
struct VarQuery {
  double eTimeDay;
  int iTime;
  std::string NatureQuery; // Can be "instant", "average", "swathMax", "swathMin"
  std::string NaturePlot; // Can be "single" or "diff"
  double TimeFrameDay;
};
/* ------------------------------ */


/* ------------------------------ */
struct PlotBound {
  std::vector<std::string> BoundSingle_var;
  std::vector<double> BoundSingle_min;
  std::vector<double> BoundSingle_max;
  std::vector<std::string> BoundDiff_var;
  std::vector<double> BoundDiff_min;
  std::vector<double> BoundDiff_max;
};
/* ------------------------------ */


/* ------------------------------ */
struct PairRecVar {
  RecVar RecVar1;
  RecVar RecVar2;
};
/* ------------------------------ */


/* ------------------------------ */
struct GraphSparseImmutable {
public:
  GraphSparseImmutable() = delete;
GraphSparseImmutable(int const& _nbVert, std::vector<int> const& _ListStart, std::vector<int> const& _ListListAdj) : nbVert(_nbVert), ListStart(_ListStart), ListListAdj(_ListListAdj)
  {
    HasVertexColor=false;
  }
  ~GraphSparseImmutable()
  {
  }
  GraphSparseImmutable(GraphSparseImmutable const& eG)
  {
    nbVert=eG.GetNbVert();
    ListStart=eG.GetListStart();
    ListListAdj=eG.GetListListAdj();
    HasVertexColor=eG.GetHasVertexColor();
    ListVertexColor=eG.GetListVertexColor();
  }
  GraphSparseImmutable operator=(GraphSparseImmutable const& eG)
  {
    nbVert=eG.GetNbVert();
    ListStart=eG.GetListStart();
    ListListAdj=eG.GetListListAdj();
    HasVertexColor=eG.GetHasVertexColor();
    ListVertexColor=eG.GetListVertexColor();
    return *this;
  }
  // lighter stuff
  int GetNbVert() const
  {
    return nbVert;
  }
  std::vector<int> GetListListAdj() const
  {
    return ListListAdj;
  }
  std::vector<int> GetListStart() const
  {
    return ListStart;
  }
  bool GetHasVertexColor() const
  {
    return HasVertexColor;
  }
  std::vector<int> GetListVertexColor() const
  {
    return ListVertexColor;
  }
  //
  void SetHasColor(bool const& TheVal)
  {
    if (TheVal == HasVertexColor) {
      return;
    }
    HasVertexColor=TheVal;
    if (TheVal == true) 
      ListVertexColor=std::vector<int>(nbVert);
    if (TheVal == false)
      ListVertexColor.clear();
  }
  void SetColor(int const& iVert, int const& eColor)
  {
    ListVertexColor[iVert]=eColor;
  }
  std::vector<int> Adjacency(int const& iVert) const
  {
    int eStart=ListStart[iVert];
    int eEnd=ListStart[iVert+1];
    std::vector<int> TheRet;
    for (int i=eStart; i<eEnd; i++)
      TheRet.push_back(ListListAdj[i]);
    return TheRet;
  }
  bool IsAdjacent(int const& iVert, int const& jVert) const
  {
    int eStart=ListStart[iVert];
    int eEnd=ListStart[iVert+1];
    for (int i=eStart; i<eEnd; i++)
      if (ListListAdj[i] == jVert)
	return true;
    return false;
  }
  int GetColor(int const& iVert) const
  {
    if (HasVertexColor == false) {
      std::cerr << "Call to GetColor while HasVertexColor=false\n";
      exit(1);
    }
    return ListVertexColor[iVert];
  }
private:
  int nbVert;
  std::vector<int> ListStart;
  std::vector<int> ListListAdj;
  bool HasVertexColor;
  std::vector<int> ListVertexColor;
};
/* ------------------------------ */


/* ------------------------------ */
PlotBound ReadPlotBound(FullNamelist const& eFull)
{
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  std::vector<std::string> BoundSingle_var=eBlPLOT.ListListStringValues.at("BoundSingle_var");
  std::vector<double> BoundSingle_min=eBlPLOT.ListListDoubleValues.at("BoundSingle_min");
  std::vector<double> BoundSingle_max=eBlPLOT.ListListDoubleValues.at("BoundSingle_max");
  std::vector<std::string> BoundDiff_var=eBlPLOT.ListListStringValues.at("BoundDiff_var");
  std::vector<double> BoundDiff_min=eBlPLOT.ListListDoubleValues.at("BoundDiff_min");
  std::vector<double> BoundDiff_max=eBlPLOT.ListListDoubleValues.at("BoundDiff_max");
  return {BoundSingle_var, BoundSingle_min, BoundSingle_max, BoundDiff_var, BoundDiff_min, BoundDiff_max};
}
/* ------------------------------ */


/* ------------------------------ */
std::string NAMELIST_ClearEndOfLine(std::string const& eStr)
{
  std::string eCharCommentB="!";
  std::string eStr3=NAMELIST_RemoveAfterLastChar(eStr, eCharCommentB);
  //
  int iPos=-1;
  int len=eStr3.size();
  std::string eLastChar=",";
  for (int i=0; i<len; i++) {
    int j=len-1-i;
    if (iPos == -1) {
      std::string eChar=eStr3.substr(j,1);
      if (eChar == eLastChar)
	iPos=j;
    }
  }
  if (iPos == -1)
    return eStr3;
  std::string eStrPrior=eStr3.substr(0, iPos);
  std::string eStrPosterior=eStr3.substr(iPos+1, len-iPos-1);
  bool test=STRING_IsStringReduceToSpace(eStrPosterior);
  if (test) {
    return eStrPrior;
  }
  else {
    return eStr3;
  }
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<double> NAMELIST_ConvertFortranStringListDoubleToCppVectorDouble(std::string const& eVarValue)
{
  //  std::cerr << "eVarValue=" << eVarValue << "\n";
  std::string eSepChar=",";
  std::vector<std::string> ListStr=STRING_Split(eVarValue, eSepChar);
  std::vector<double> eListRetDouble;
  int siz=ListStr.size();
  for (int i=0; i<siz; i++) {
    std::string eStr1=ListStr[i];
    std::string eStr2=STRING_RemoveSpacesBeginningEnd(eStr1);
    double eVal;
    std::istringstream(eStr2) >> eVal;
    eListRetDouble.push_back(eVal);
  }
  //  std::cerr << "eListRetDouble=";
  //  for (auto& eVal : eListRetDouble)
  //    std::cerr << " " << eVal;
  //  std::cerr << "\n";
  return eListRetDouble;
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<int> NAMELIST_ConvertFortranStringListIntToCppVectorInt(std::string const& eVarValue)
{
  //  std::cerr << "eVarValue=" << eVarValue << "\n";
  std::string eSepChar=",";
  std::vector<std::string> ListStr=STRING_Split(eVarValue, eSepChar);
  std::vector<int> eListRetInt;
  int siz=ListStr.size();
  for (int i=0; i<siz; i++) {
    std::string eStr1=ListStr[i];
    std::string eStr2=STRING_RemoveSpacesBeginningEnd(eStr1);
    int eVal;
    std::istringstream(eStr2) >> eVal;
    eListRetInt.push_back(eVal);
  }
  //  std::cerr << "eListRetDouble=";
  //  for (auto& eVal : eListRetDouble)
  //    std::cerr << " " << eVal;
  //  std::cerr << "\n";
  return eListRetInt;
}
/* ------------------------------ */


/* ------------------------------ */
void NAMELIST_ReadNamelistFile(std::string const& eFileName, FullNamelist &eFullNamelist)
{
  if (IsExistingFile(eFileName) == false) {
    std::cerr << "The following namelist file is missing\n";
    std::cerr << "eFileName = " << eFileName << "\n";
    exit(1);
  }
  std::ifstream INfs;
  INfs.open(eFileName);
  int InBlock=0;
  std::string eBlockName;
  int iLine=0;
  while(!INfs.eof()) {
    std::string Ampersand="&";
    std::string strTab="\t";
    int siz=1024;
    char eLine[siz];
    INfs.getline(eLine, siz);
    std::string PreStr=eLine;
    //    std::cerr << "iLine=" << iLine << "\n";
    //    std::cerr << "PreStr=" << PreStr << "\n";
    std::string eCharComment="!";
    std::string PreStrB=NAMELIST_RemoveAfterCommentChar(PreStr, eCharComment);
    //    std::cerr << "PreStrB=" << PreStrB << "\n";
    std::string eStr=STRING_RemoveSpacesBeginningEnd(PreStrB);
    //    std::cerr << "eStr=" << eStr << "\n";
    int len=eStr.length();
    if (eStr.find(strTab) != std::string::npos) {
      std::cerr << "Tabs are not allowed\n";
      std::cerr << "LINE=" << eStr << "\n";
      exit(1);
    }
    if (len> 0) {
      if (eStr.find(Ampersand) != std::string::npos) {
	std::string eFirstChar=eStr.substr(0,1);
	if (eFirstChar != "&") {
	  std::cerr << "Error while reading namelist file = " << eFileName << "\n";
	  std::cerr << "Error, Ampersand (&) should be only in the first character\n";
	  std::cerr << "LINE=" << eStr << "\n";
	  exit(1);
	}
	std::string strRed=eStr.substr(1, len-1);
	if (InBlock == 0) {
	  eBlockName=strRed;
          auto search=eFullNamelist.ListBlock.find(eBlockName);
          if (search == eFullNamelist.ListBlock.end() ) {
	    std::cerr << "Find BlockName = " << eBlockName << "\n";
	    std::cerr << "which is not in the authorized list\n";
	    std::cerr << "LINE=" << eStr << "\n";
	    std::cerr << "List of authorized block names:\n";
            for (auto & eBlock : eFullNamelist.ListBlock) {
	      std::cerr << "Block name=" << eBlock.first << "\n";
	    }
	    exit(1);
	  }
	  InBlock=1;
	}
	else {
	  if (strRed != "END") {
	    std::cerr << "Ampersand detected. We should leave with a END\n";
	    std::cerr << "LINE=" << eStr << "\n";
	    exit(1);
	  }
	  InBlock=0;
	}
      }
      else {
	if (eStr != "/") {
	  std::string eStr3=NAMELIST_ClearEndOfLine(eStr);
	  //	  std::cerr << "eStr3=" << eStr3 << "\n";
	  std::string strEqual="=";
	  int posEqual=STRING_GetCharPositionInString(eStr3, strEqual);
	  if (posEqual != -1) {
	    int len3=eStr3.length();
	    //	    std::cerr << "eStr=" << eStr << "\n";
	    //	    std::cerr << "eStr3=" << eStr3 << "\n";
	    //	    std::cerr << "posEqual=" << posEqual << " len3=" << len3 << "\n";
	    std::string eStrPrior=eStr3.substr(0, posEqual);
	    std::string eStrPosterior=eStr3.substr(posEqual+1, len3-posEqual-1);
	    //	    std::cerr << "eStrPrior=" << eStrPrior << "\n";
	    //	    std::cerr << "eStrPosterior=" << eStrPosterior << "\n";
	    std::string eVarName=STRING_RemoveSpacesBeginningEnd(eStrPrior);
	    std::string eVarValue=STRING_RemoveSpacesBeginningEnd(eStrPosterior);
	    //	    std::cerr << "eVarName=" << eVarName << "\n";
	    //	    std::cerr << "eVarValue=" << eVarValue << "\n";
	    std::string eVarNature=NAMELIST_FindPositionVariableInBlock(
		      eVarName, eFullNamelist.ListBlock[eBlockName]);
	    if (eVarNature == "not found") {
	      NAMELIST_WriteNamelistFile(std::cerr, eFullNamelist);
	      std::cerr << "Error in reading the NAMELIST file. See above the allowed\n";
	      std::cerr << "The variable " << eVarName << "\n";
	      std::cerr << "is in block " << eBlockName << "\n";
	      std::cerr << "of the file " << eFileName << "\n";
	      std::cerr << "but it is not allowed for the chosen application\n";

	      exit(1);
	    }
	    if (eVarNature == "int") {
	      int eVal;
	      std::istringstream(eVarValue) >> eVal;
	      eFullNamelist.ListBlock[eBlockName].ListIntValues[eVarName]=eVal;
	    }
	    if (eVarNature == "bool") {
	      bool eVal=NAMELIST_ReadBoolValue(eVarValue);
	      eFullNamelist.ListBlock[eBlockName].ListBoolValues[eVarName]=eVal;
	    }
	    if (eVarNature == "double") {
	      double eVal;
	      std::istringstream(eVarValue) >> eVal;
	      eFullNamelist.ListBlock[eBlockName].ListDoubleValues[eVarName]=eVal;
	    }
	    if (eVarNature == "string") {
	      std::string eVal=NAMELIST_ConvertFortranStringToCppString(eVarValue);
	      eFullNamelist.ListBlock[eBlockName].ListStringValues[eVarName]=eVal;
	    }
	    if (eVarNature == "listdouble") {
	      std::vector<double> eVal=NAMELIST_ConvertFortranStringListDoubleToCppVectorDouble(eVarValue);
	      eFullNamelist.ListBlock[eBlockName].ListListDoubleValues[eVarName]=eVal;
	    }
	    if (eVarNature == "listint") {
	      std::vector<int> eVal=NAMELIST_ConvertFortranStringListIntToCppVectorInt(eVarValue);
	      eFullNamelist.ListBlock[eBlockName].ListListIntValues[eVarName]=eVal;
	    }
	    if (eVarNature == "liststring") {
	      std::vector<std::string> eVal=NAMELIST_ConvertFortranListStringToCppListString(eVarValue);
	      eFullNamelist.ListBlock[eBlockName].ListListStringValues[eVarName]=eVal;
	    }
	  }
	}
	else {
	  InBlock=0;	  
	}
      }
    }
    iLine++;
  }
  INfs.close();
}
/* ------------------------------ */


/* ------------------------------ */
FullNamelist NAMELIST_GetStandard_PlotRoutine_common()
{
  std::map<std::string, SingleBlock> ListBlock;
  // PROC
  std::string BlockName1="PROC";
  std::map<std::string, int> ListIntValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, double> ListDoubleValues1;
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, std::vector<std::string> > ListListStringValues1;
  ListStringValues1["MODELNAME"]="COSMO or WAM or ROMS";
  ListStringValues1["BEGTC"]="20110915.000000";
  ListStringValues1["ENDTC"]="20110925.000000";
  ListDoubleValues1["DELTC"]=600;
  ListStringValues1["UNITC"]="SEC";
  ListStringValues1["GridFile"]="roms_grid.nc";
  ListBoolValues1["CutWorldMap"]=false;
  ListBoolValues1["HigherLatitudeCut"]=false;
  ListDoubleValues1["MinLatCut"]=-80;
  ListDoubleValues1["MaxLatCut"]=80;
  ListStringValues1["PicPrefix"]="Pictures/DIR_plot/";
  ListStringValues1["Extension"]="png";
  ListListStringValues1["ListNatureQuery"]={"instant"}; // By default instantaneous values
  ListDoubleValues1["TimeFrameDay"]=1;
  ListBoolValues1["FirstCleanDirectory"]=true;
  ListBoolValues1["KeepNC_NCL"]=false;
  ListIntValues1["NPROC"]=1;
  SingleBlock BlockPROC;
  BlockPROC.ListIntValues=ListIntValues1;
  BlockPROC.ListBoolValues=ListBoolValues1;
  BlockPROC.ListDoubleValues=ListDoubleValues1;
  BlockPROC.ListStringValues=ListStringValues1;
  BlockPROC.ListListStringValues=ListListStringValues1;
  BlockPROC.BlockName=BlockName1;
  ListBlock["PROC"]=BlockPROC;
  // PLOT
  std::string BlockName2="PLOT";
  std::map<std::string, int> ListIntValues2;
  std::map<std::string, bool> ListBoolValues2;
  std::map<std::string, double> ListDoubleValues2;
  std::map<std::string, std::string> ListStringValues2;
  std::map<std::string, std::vector<double> > ListListDoubleValues2;
  std::map<std::string, std::vector<std::string> > ListListStringValues2;
  ListStringValues2["ColorMap"]="BlAqGrYeOrReVi200";
  ListStringValues2["ColorMapDiff"]="BlWhRe";
  ListStringValues2["cnFillMode"]="RasterFill";
  ListBoolValues2["cnSmoothingOn"]=true;
  ListIntValues2["nbLevelSpa"]=50;
  ListIntValues2["nbLabelStride"]=10;
  ListBoolValues2["UseNativeGrid"]=true;
  ListBoolValues2["DoTitle"]=true;
  ListStringValues2["GridResolution"]="HighRes";
  ListBoolValues2["DrawRiver"]=false;
  ListBoolValues2["PrintMMA"]=false;
  ListBoolValues2["DrawContourBathy"]=false;
  ListBoolValues2["DrawAnnotation"]=false;
  ListDoubleValues2["AnnotationLon"]=0;
  ListDoubleValues2["AnnotationLat"]=0;
  ListStringValues2["AnnotationText"]="something to write";
  ListListStringValues2["BoundSingle_var"]={};
  ListListDoubleValues2["BoundSingle_min"]={};
  ListListDoubleValues2["BoundSingle_max"]={};
  ListListStringValues2["BoundDiff_var"]={};
  ListListDoubleValues2["BoundDiff_min"]={};
  ListListDoubleValues2["BoundDiff_max"]={};
  ListBoolValues2["FillLand"]=true;
  SingleBlock BlockPLOT;
  BlockPLOT.ListIntValues=ListIntValues2;
  BlockPLOT.ListBoolValues=ListBoolValues2;
  BlockPLOT.ListDoubleValues=ListDoubleValues2;
  BlockPLOT.ListStringValues=ListStringValues2;
  BlockPLOT.ListListStringValues=ListListStringValues2;
  BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
  BlockPLOT.BlockName=BlockName2;
  ListBlock["PLOT"]=BlockPLOT;
  // VARS
  std::string BlockName3="VARS";
  std::map<std::string, int> ListIntValues3;
  std::map<std::string, bool> ListBoolValues3;
  std::map<std::string, double> ListDoubleValues3;
  std::map<std::string, std::string> ListStringValues3;
  std::map<std::string, std::vector<std::string> > ListListStringValues3;
  std::vector<std::string> ListVarOut=GetAllPossibleVariables();
  for (auto& eVal : ListVarOut)
    ListBoolValues3[eVal]=false;
  SingleBlock BlockVARS;
  BlockVARS.ListIntValues=ListIntValues3;
  BlockVARS.ListBoolValues=ListBoolValues3;
  BlockVARS.ListDoubleValues=ListDoubleValues3;
  BlockVARS.ListStringValues=ListStringValues3;
  BlockVARS.ListListStringValues=ListListStringValues3;
  BlockVARS.BlockName=BlockName3;
  ListBlock["VARS"]=BlockVARS;
  // Final part
  FullNamelist eFullNamelist;
  eFullNamelist.ListBlock=ListBlock;
  eFullNamelist.FileName="undefined";
  return eFullNamelist;
}
/* ------------------------------ */


/* ------------------------------ */
FullNamelist NAMELIST_GetStandard_PlotRoutine_single()
{
  FullNamelist eFull=NAMELIST_GetStandard_PlotRoutine_common();
  eFull.ListBlock["PROC"].ListStringValues["__NaturePlot"]="SINGLE";
  eFull.ListBlock["PROC"].ListStringValues["HisPrefix"]="COSMO_output_ or WAM_output_ or ROMS_output_";
  return eFull;
}
/* ------------------------------ */


/* ------------------------------ */
void InitializeIdxJdxWet(CoordGridArrayFD & eCoordGrdArr)
{
  int eta=eCoordGrdArr.eta;
  int xi=eCoordGrdArr.xi;
  int nbWet=0;
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++)
      if (eCoordGrdArr.MSK(i, j) == 1) {
        nbWet++;
        eCoordGrdArr.Idx.push_back(i);
        eCoordGrdArr.Jdx.push_back(j);
      }
  eCoordGrdArr.nbWet=nbWet;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray NC_ReadRomsGridFile(std::string const& eFile)
{
  GridArray eRomsGridArray;
  eRomsGridArray.ModelName="ROMS";
  eRomsGridArray.IsFE=0;
  eRomsGridArray.IsLonLat=1;
  // Rho part of the arrays
  eRomsGridArray.GrdArrRho.LON=NC_Read2Dvariable(eFile, "lon_rho");
  eRomsGridArray.GrdArrRho.LAT=NC_Read2Dvariable(eFile, "lat_rho");
  eRomsGridArray.GrdArrRho.DEP=NC_Read2Dvariable(eFile, "h");
  eRomsGridArray.GrdArrRho.HaveDEP=true;
  eRomsGridArray.GrdArrRho.ANG=NC_Read2Dvariable(eFile, "angle");
  int eta_rho=eRomsGridArray.GrdArrRho.LON.rows();
  int xi_rho=eRomsGridArray.GrdArrRho.LON.cols();
  eRomsGridArray.GrdArrRho.eta=eta_rho;
  eRomsGridArray.GrdArrRho.xi =xi_rho;
  MyMatrix<double> eMSK_rho_double=NC_Read2Dvariable(eFile, "mask_rho");
  MyMatrix<int> eMSK_rho(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++)
      eMSK_rho(i, j)=int(eMSK_rho_double(i, j));
  eRomsGridArray.GrdArrRho.MSK=eMSK_rho;
  InitializeIdxJdxWet(eRomsGridArray.GrdArrRho);
  // U
  eRomsGridArray.GrdArrU.LON=NC_Read2Dvariable(eFile, "lon_u");
  eRomsGridArray.GrdArrU.LAT=NC_Read2Dvariable(eFile, "lat_u");
  int eta_u=eRomsGridArray.GrdArrU.LON.rows();
  int xi_u=eRomsGridArray.GrdArrU.LON.cols();
  eRomsGridArray.GrdArrU.eta=eta_u;
  eRomsGridArray.GrdArrU.xi =xi_u;
  MyMatrix<int> MSKu(eta_u, xi_u);
  MyMatrix<double> DEPu(eta_u, xi_u);
  MyMatrix<double> ANGu(eta_u, xi_u);
  for (int i=0; i<eta_u; i++)
    for (int j=0; j<xi_u; j++) {
      DEPu(i, j)=
	(eRomsGridArray.GrdArrRho.DEP(i, j)+
	 eRomsGridArray.GrdArrRho.DEP(i, j+1) )/double(2);
      ANGu(i, j)=
	(eRomsGridArray.GrdArrRho.ANG(i, j)+
	 eRomsGridArray.GrdArrRho.ANG(i, j+1) )/double(2);
      MSKu(i, j)=
	eRomsGridArray.GrdArrRho.MSK(i, j)*
	eRomsGridArray.GrdArrRho.MSK(i, j+1);
    }
  eRomsGridArray.GrdArrU.MSK=MSKu;
  eRomsGridArray.GrdArrU.DEP=DEPu;
  eRomsGridArray.GrdArrU.HaveDEP=true;
  eRomsGridArray.GrdArrU.ANG=ANGu;
  InitializeIdxJdxWet(eRomsGridArray.GrdArrU);
  // V
  eRomsGridArray.GrdArrV.LON=NC_Read2Dvariable(eFile, "lon_v");
  eRomsGridArray.GrdArrV.LAT=NC_Read2Dvariable(eFile, "lat_v");
  int eta_v=eRomsGridArray.GrdArrV.LON.rows();
  int xi_v=eRomsGridArray.GrdArrV.LON.cols();
  eRomsGridArray.GrdArrV.eta=eta_v;
  eRomsGridArray.GrdArrV.xi =xi_v;
  MyMatrix<int> MSKv(eta_v, xi_v);
  MyMatrix<double> DEPv(eta_v, xi_v);
  MyMatrix<double> ANGv(eta_v, xi_v);
  for (int i=0; i<eta_v; i++)
    for (int j=0; j<xi_v; j++) {
      DEPv(i, j)=
	(eRomsGridArray.GrdArrRho.DEP(i, j)+
	 eRomsGridArray.GrdArrRho.DEP(i+1,j) )/double(2);
      ANGv(i, j)=
	(eRomsGridArray.GrdArrRho.ANG(i, j)+
	 eRomsGridArray.GrdArrRho.ANG(i+1, j) )/double(2);
      MSKv(i, j)=
	eRomsGridArray.GrdArrRho.MSK(i, j)*
	eRomsGridArray.GrdArrRho.MSK(i+1, j);
    }
  eRomsGridArray.GrdArrV.MSK=MSKv;
  eRomsGridArray.GrdArrV.DEP=DEPv;
  eRomsGridArray.GrdArrV.HaveDEP=true;
  eRomsGridArray.GrdArrV.ANG=ANGv;
  InitializeIdxJdxWet(eRomsGridArray.GrdArrV);
  // PSI
  eRomsGridArray.GrdArrPsi.LON=NC_Read2Dvariable(eFile, "lon_psi");
  eRomsGridArray.GrdArrPsi.LAT=NC_Read2Dvariable(eFile, "lat_psi");
  int eta_psi=eRomsGridArray.GrdArrPsi.LON.rows();
  int xi_psi=eRomsGridArray.GrdArrPsi.LON.cols();
  eRomsGridArray.GrdArrPsi.eta=eta_psi;
  eRomsGridArray.GrdArrPsi.xi =xi_psi;
  MyMatrix<int> MSKp(eta_psi, xi_psi);
  MyMatrix<double> DEPp(eta_psi, xi_psi);
  MyMatrix<double> ANGp(eta_psi, xi_psi);
  for (int i=0; i<eta_psi; i++)
    for (int j=0; j<xi_psi; j++) {
      DEPp(i, j)=
	(eRomsGridArray.GrdArrRho.DEP(i  ,j+1)+
	 eRomsGridArray.GrdArrRho.DEP(i+1,j+1)+
	 eRomsGridArray.GrdArrRho.DEP(i  ,j  )+
	 eRomsGridArray.GrdArrRho.DEP(i+1,j  ))/double(4);
      ANGp(i, j)=
	(eRomsGridArray.GrdArrRho.ANG(i  ,j+1)+
	 eRomsGridArray.GrdArrRho.ANG(i+1,j+1)+
	 eRomsGridArray.GrdArrRho.ANG(i  ,j  )+
	 eRomsGridArray.GrdArrRho.ANG(i+1,j  ))/double(4);
      MSKp(i, j)=
	 eRomsGridArray.GrdArrRho.MSK(i  ,j+1)*
	 eRomsGridArray.GrdArrRho.MSK(i+1,j+1)*
	 eRomsGridArray.GrdArrRho.MSK(i  ,j  )*
	 eRomsGridArray.GrdArrRho.MSK(i+1,j  );
    }
  eRomsGridArray.GrdArrPsi.MSK=MSKp;
  eRomsGridArray.GrdArrPsi.DEP=DEPp;
  eRomsGridArray.GrdArrPsi.HaveDEP=true;
  eRomsGridArray.GrdArrPsi.ANG=ANGp;
  InitializeIdxJdxWet(eRomsGridArray.GrdArrPsi);
  std::cerr << "The ROMS grid has been read\n";
  return eRomsGridArray;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray NC_ReadCosmoWamStructGridFile(std::string const& eFile, std::string const& postfix)
{
  GridArray eCosmoWamGridArray;
  if (postfix == "atm") {
    eCosmoWamGridArray.ModelName="COSMO";
  }
  else {
    eCosmoWamGridArray.ModelName="WAM";
  }
  eCosmoWamGridArray.IsFE=0;
  eCosmoWamGridArray.IsLonLat=1;
  // Rho part of the arrays
  std::string LONstr="LON_" + postfix;
  //  std::cerr << "Before reading " << LONstr << "\n";
  eCosmoWamGridArray.GrdArrRho.LON=NC_Read2Dvariable(eFile, LONstr);
  std::string LATstr="LAT_" + postfix;
  //  std::cerr << "Before reading " << LATstr << "\n";
  eCosmoWamGridArray.GrdArrRho.LAT=NC_Read2Dvariable(eFile, LATstr);
  int eta_rho=eCosmoWamGridArray.GrdArrRho.LON.rows();
  int xi_rho=eCosmoWamGridArray.GrdArrRho.LON.cols();
  //  std::cerr << "eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
  // The bathymetry if available
  std::string DEPstr="DEP_" + postfix;
  if (NC_IsVar(eFile, DEPstr) ) {
    eCosmoWamGridArray.GrdArrRho.DEP=NC_Read2Dvariable(eFile, DEPstr);
    eCosmoWamGridArray.GrdArrRho.HaveDEP=true;
  }
  else {
    eCosmoWamGridArray.GrdArrRho.DEP=ZeroMatrix<double>(eta_rho, xi_rho);
    eCosmoWamGridArray.GrdArrRho.HaveDEP=false;
  }
  // The mask if available
  std::string MSKstr="MSK_" + postfix;
  MyMatrix<double> MSK_double(eta_rho, xi_rho);
  if (NC_IsVar(eFile, DEPstr) ) {
    MSK_double=NC_Read2Dvariable(eFile, MSKstr);
  }
  else {
    for (int i=0; i<eta_rho; i++)
      for (int j=0; j<xi_rho; j++) {
	MSK_double(i,j)=double(1);
      }
  }
  MyMatrix<int> MSK_int(eta_rho, xi_rho);
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      MSK_int(i,j)=int(MSK_double(i,j));
    }
  eCosmoWamGridArray.GrdArrRho.MSK=MSK_int;
  // The angle if available
  //  std::cerr << "Before reading angle\n";
  if (postfix == "atm") {
    eCosmoWamGridArray.GrdArrRho.ANG=NC_Read2Dvariable(eFile, "ANG_atm");
  }
  else {
    eCosmoWamGridArray.GrdArrRho.ANG=ZeroMatrix<double>(eta_rho, xi_rho);
  }
  eCosmoWamGridArray.GrdArrRho.eta=eta_rho;
  eCosmoWamGridArray.GrdArrRho.xi =xi_rho;
  // U / V / PSI we do not need a priori
  return eCosmoWamGridArray;
}
/* ------------------------------ */


/* ------------------------------ */
void CHECK_UnstructuredGrid(GridArray const& GrdArr)
{
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int IP=GrdArr.INE(ie,i);
      if (IP < 0 || IP >= mnp) {
	std::cerr << "Error in the unstructured grid\n";
	std::cerr << "mnp=" << mnp << "  mne=" << mne << "\n";
	std::cerr << "ie=" << ie << "\n";
	std::cerr << "INE=[" << GrdArr.INE(ie,0) << " , " << GrdArr.INE(ie,1) << " , " << GrdArr.INE(ie,2) << "]\n";
	exit(1);
      }
    }
}
/* ------------------------------ */


/* ------------------------------ */
void CHECK_CombinatorialGrid(GridArray const& GrdArr)
{
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  int CCON[mnp];
  int POS_TRICK[3][2];
  POS_TRICK[0][0] = 1;
  POS_TRICK[1][0] = 2;
  POS_TRICK[2][0] = 0;
  POS_TRICK[0][1] = 2;
  POS_TRICK[1][1] = 0;
  POS_TRICK[2][1] = 1;
  for (int ip=0; ip<mnp; ip++) {
    CCON[ip]=0;
  }
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int ip=GrdArr.INE(ie,i);
      CCON[ip]++;
    }
  int MAXMNECON=0;
  for (int ip=0; ip<mnp; ip++) {
    int eCon=CCON[ip];
    if (eCon > MAXMNECON)
      MAXMNECON=eCon;
  }
  int CHILF[mnp];
  for (int ip=0; ip<mnp; ip++) {
    CHILF[ip]=0;
  }
  int CELLVERTEX[mnp][MAXMNECON][2];
  for (int ie=0; ie<mne; ie++)
    for (int j=0; j<3; j++) {
      int i=GrdArr.INE(ie,j);
      CELLVERTEX[i][ CHILF[i] ][0] = ie;
      CELLVERTEX[i][ CHILF[i] ][1] = j;
      CHILF[i]++;
    }
  int COUNT_MAX=0;
  for (int ip=0; ip<mnp; ip++)
    COUNT_MAX += CCON[ip];
  //  int IE_CELL[COUNT_MAX];
  //  int POS_CELL[COUNT_MAX];
  int IE_CELL2[mnp][MAXMNECON];
  int POS_CELL2[mnp][MAXMNECON];
  int j=0;
  for (int ip=0; ip<mnp; ip++)
    for (int i=0; i<CCON[ip]; i++) {
      //      IE_CELL[j] = CELLVERTEX[ip][i][0];
      //      POS_CELL[j] = CELLVERTEX[ip][i][1];
      IE_CELL2[ip][i] = CELLVERTEX[ip][i][0];
      POS_CELL2[ip][i] = CELLVERTEX[ip][i][1];
      j++;
    }
  for (int ie=0; ie<mne; ie++)
    for (int i=0; i<3; i++) {
      int INEXT=POS_TRICK[i][0];
      int ip=GrdArr.INE(ie,i);
      int IP_NEXT=GrdArr.INE(ie, INEXT);
      int nbMatch=0;
      std::vector<int> Lmatch;
      for (int icon=0; icon<CCON[ip]; icon++) {
	int ie2=IE_CELL2[ip][icon];
	if (ie != ie2) {
	  int POS=POS_CELL2[ip][icon];
	  int POS_NEXT=POS_TRICK[POS][0];
	  int IP_ADJ_NEXT=GrdArr.INE(ie2, POS_NEXT);
	  if (IP_ADJ_NEXT == IP_NEXT) {
	    std::cerr << "Combinatorial orientability problem\n";
	    std::cerr << "IE=" << ie << " IE2=" << ie2 << "\n";
	    std::cerr << "IP=" << ip << " IP_NEXT=" << IP_NEXT << "\n";
	    exit(1);
	  }
	  int POS_PREV=POS_TRICK[POS][1];
	  int IP_ADJ_PREV=GrdArr.INE(ie2, POS_PREV);
	  if (IP_ADJ_PREV == IP_NEXT) {
	    nbMatch++;
	    Lmatch.push_back(ie2);
	  }
	}
      }
      if (nbMatch > 1) {
	std::cerr << "nbMatch is too large.\n";
	std::cerr << "Should be 0 for boundary edge\n";
	std::cerr << "Should be 1 for interior edges\n";
	std::cerr << "ie=" << ie << " i=" << i << "\n";
	std::cerr << " ip=" << ip << " IP_NEXT=" << IP_NEXT << "\n";
	std::cerr << "nbMatch=" << nbMatch << "\n";
	for (int iMatch=0; iMatch<nbMatch; iMatch++) {
	  int iem=Lmatch[iMatch];
	  std::cerr << "  iMatch=" << iMatch << " ie=" << iem << "\n";
	  std::cerr << "     ine=[" << GrdArr.INE(iem,0) << "," << GrdArr.INE(iem,1) << "," << GrdArr.INE(iem,2) << "]\n";
	}
	exit(1);
      }
    }
  std::cerr << "Now leaving the combinatorial check\n";
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<int> NC_ReadElements(std::string const& eFile, std::string const& eStr)
{
  MyMatrix<int> INE=NC_Read2Dvariable_int(eFile, eStr);
  int nbRow=INE.rows();
  int nbCol=INE.cols();
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++) {
      int IP=INE(iRow, iCol);
      INE(iRow, iCol)=IP-1;
    }
  return INE;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray NC_ReadWamGridFile(std::string const& eFile)
{
  MyVector<int> eVectLLUNSTR=NC_Read1Dvariable_int(eFile, "LLUNSTR");
  int LLUNSTR=eVectLLUNSTR(0);
  //  std::cerr << "LLUNSTR=" << LLUNSTR << "\n";
  if (LLUNSTR == 0)
    return NC_ReadCosmoWamStructGridFile(eFile, "wav");
  GridArray eGridArray;
  eGridArray.IsFE=1;
  eGridArray.IsLonLat=1;
  eGridArray.L_IndexSelect=false;
  //
  eGridArray.INE=NC_ReadElements(eFile, "ele");
  //  std::cerr << "NC_ReadWam, step 1\n";
  MyVector<double> LON=NC_Read1Dvariable(eFile, "LON_wav");
  //  std::cerr << "NC_ReadWam, step 2\n";
  MyVector<double> LAT=NC_Read1Dvariable(eFile, "LAT_wav");
  //  std::cerr << "NC_ReadWam, step 3\n";
  MyVector<double> DEP=NC_Read1Dvariable(eFile, "DEP_wav");
  //  std::cerr << "NC_ReadWam, step 4\n";
  int nbPoint=LON.size();
  MyMatrix<double> LONarr(nbPoint,1);
  MyMatrix<double> LATarr(nbPoint,1);
  MyMatrix<double> DEParr(nbPoint,1);
  MyMatrix<double> ANGarr(nbPoint,1);
  MyMatrix<int> MSKarr(nbPoint,1);
  //  std::cerr << "NC_ReadWam, step 5\n";
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    LONarr(iPoint,0)=LON(iPoint);
    LATarr(iPoint,0)=LAT(iPoint);
    DEParr(iPoint,0)=DEP(iPoint);
    ANGarr(iPoint,0)=0;
    MSKarr(iPoint,0)=1;
  }
  eGridArray.GrdArrRho.LON=LONarr;
  eGridArray.GrdArrRho.LAT=LATarr;
  eGridArray.GrdArrRho.DEP=DEParr;
  eGridArray.GrdArrRho.ANG=ANGarr;
  eGridArray.GrdArrRho.MSK=MSKarr;
  //  std::cerr << "NC_ReadWam, step 6\n";
  CHECK_UnstructuredGrid(eGridArray);
  CHECK_CombinatorialGrid(eGridArray);
  //  std::cerr << "NC_ReadWam, step 7\n";
  return eGridArray;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray WWM_ReadGridFile_netcdf(std::string const& GridFile)
{
  GridArray eGridArray;
  eGridArray.ModelName="WWM";
  eGridArray.IsFE=1;
  eGridArray.IsLonLat=1;
  eGridArray.L_IndexSelect=false;
  //
  if (IsExistingFile(GridFile) == false) {
    std::cerr << "Error in WWM_ReadGridFile_netcdf\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    exit(1);
  }
  //  std::cerr << "WWM_ReadGridFile_netcdf, step 1\n";
  eGridArray.INE=NC_ReadElements(GridFile, "ele");
  //  std::cerr << "WWM_ReadGridFile_netcdf, step 2\n";
  MyVector<int> LType=NC_Read1Dvariable_int(GridFile, "LSPHE");
  //  std::cerr << "WWM_ReadGridFile_netcdf, step 3\n";
  int LSPHE=LType(0);
  std::string Xname, Yname;
  if (LSPHE == 1) {
    Xname="lon";
    Yname="lat";
  }
  else {
    Xname="x";
    Yname="y";
  }
  //  std::cerr << "WWM_ReadGridFile_netcdf, step 4\n";
  MyVector<double> LON=NC_Read1Dvariable(GridFile, Xname);
  //  std::cerr << "WWM_ReadGridFile_netcdf, step 5\n";
  MyVector<double> LAT=NC_Read1Dvariable(GridFile, Yname);
  //  std::cerr << "WWM_ReadGridFile_netcdf, step 6\n";
  MyVector<double> DEP=NC_Read1Dvariable(GridFile, "depth");
  //  std::cerr << "WWM_ReadGridFile_netcdf, step 7\n";
  int nbPoint=LON.size();
  MyMatrix<double> LONarr(nbPoint,1);
  MyMatrix<double> LATarr(nbPoint,1);
  MyMatrix<double> DEParr(nbPoint,1);
  MyMatrix<double> ANGarr(nbPoint,1);
  MyMatrix<int> MSKarr(nbPoint,1);
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    LONarr(iPoint,0)=LON(iPoint);
    LATarr(iPoint,0)=LAT(iPoint);
    DEParr(iPoint,0)=DEP(iPoint);
    ANGarr(iPoint,0)=0;
    MSKarr(iPoint,0)=1;
  }
  eGridArray.GrdArrRho.LON=LONarr;
  eGridArray.GrdArrRho.LAT=LATarr;
  eGridArray.GrdArrRho.DEP=DEParr;
  eGridArray.GrdArrRho.ANG=ANGarr;
  eGridArray.GrdArrRho.MSK=MSKarr;
  CHECK_UnstructuredGrid(eGridArray);
  return eGridArray;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray WWM_ReadGridFile_gr3(std::string const& GridFile)
{
  GridArray eGridArray;
  eGridArray.ModelName="WWM";
  eGridArray.IsFE=1;
  eGridArray.IsLonLat=1;
  eGridArray.L_IndexSelect=false;
  //
  if (IsExistingFile(GridFile) == false) {
    std::cerr << "Error in WWM_ReadGridFile_gr3\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    exit(1);
  }
  std::ifstream IN(GridFile);
  // read first line
  std::string line;
  std::getline(IN, line);
  std::cerr << "line=" << line << "\n";
  //
  int mne, mnp;
  IN >> mne;
  IN >> mnp;
  std::cerr << "mne=" << mne << " mnp=" << mnp << "\n";
  eGridArray.INE=MyMatrix<int>(mne,3);
  eGridArray.GrdArrRho.LON=MyMatrix<double>(mnp,1);
  eGridArray.GrdArrRho.LAT=MyMatrix<double>(mnp,1);
  eGridArray.GrdArrRho.DEP=MyMatrix<double>(mnp,1);
  eGridArray.GrdArrRho.ANG=MyMatrix<double>(mnp,1);
  eGridArray.GrdArrRho.MSK=MyMatrix<int>(mnp,1);
  for (int iP=0; iP<mnp; iP++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    //    std::cerr << "iP=" << iP << " XYZ=" << XPDTMP << " " << YPDTMP << " " << ZPDTMP << "\n";
    eGridArray.GrdArrRho.LON(iP,0)=XPDTMP;
    eGridArray.GrdArrRho.LAT(iP,0)=YPDTMP;
    eGridArray.GrdArrRho.DEP(iP,0)=ZPDTMP;
    eGridArray.GrdArrRho.ANG(iP,0)=0;
    eGridArray.GrdArrRho.MSK(iP,0)=1;
  }
  for (int iE=0; iE<mne; iE++) {
    int KTMP, LTMP, ip1, ip2, ip3;
    IN >> KTMP >> LTMP >> ip1 >> ip2 >> ip3;
    //    std::cerr << "iE=" << iE << " IP123=" << ip1 << " " << ip2 << " " << ip3 << "\n";
    eGridArray.INE(iE,0)=ip1 - 1;
    eGridArray.INE(iE,1)=ip2 - 1;
    eGridArray.INE(iE,2)=ip3 - 1;
  }
  IN.close();
  CHECK_UnstructuredGrid(eGridArray);
  return eGridArray;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray WWM_ReadGridFile_xfn(std::string const& GridFile)
{
  GridArray eGridArray;
  eGridArray.ModelName="WWM";
  eGridArray.IsFE=1;
  eGridArray.IsLonLat=1;
  eGridArray.L_IndexSelect=false;
  //
  if (IsExistingFile(GridFile) == false) {
    std::cerr << "Error in WWM_ReadGridFile_xfn\n";
    std::cerr << "GridFile = " << GridFile << "\n";
    std::cerr << "is missing\n";
    exit(1);
  }
  std::string line;
  int ITMP, JTMP;
  std::ifstream IN(GridFile.c_str());
  for (int i=0; i<2; i++)
    std::getline(IN, line);
  IN >> ITMP;
  std::getline(IN, line);
  IN >> JTMP;
  int mnp=ITMP + JTMP;
  for (int i=0; i<7; i++)
    std::getline(IN, line);
  eGridArray.GrdArrRho.LON=MyMatrix<double>(mnp,1);
  eGridArray.GrdArrRho.LAT=MyMatrix<double>(mnp,1);
  eGridArray.GrdArrRho.DEP=MyMatrix<double>(mnp,1);
  eGridArray.GrdArrRho.ANG=MyMatrix<double>(mnp,1);
  eGridArray.GrdArrRho.MSK=MyMatrix<int>(mnp,1);
  for (int iP=0; iP<mnp; iP++) {
    int KTMP;
    double XPDTMP, YPDTMP, ZPDTMP;
    IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
    eGridArray.GrdArrRho.LON(iP,0)=XPDTMP;
    eGridArray.GrdArrRho.LAT(iP,0)=YPDTMP;
    eGridArray.GrdArrRho.DEP(iP,0)=ZPDTMP;
    eGridArray.GrdArrRho.ANG(iP,0)=0;
    eGridArray.GrdArrRho.MSK(iP,0)=1;
  }
  for (int i=0; i<2; i++)
    std::getline(IN, line);
  int mne;
  IN >> mne;
  for (int i=0; i<3; i++)
    std::getline(IN, line);
  for (int iE=0; iE<mne; iE++) {
    int KTMP, LTMP, ip1, ip2, ip3;
    IN >> KTMP >> LTMP >> ip1 >> ip2 >> ip3;
    eGridArray.INE(iE,0)=ip1 - 1;
    eGridArray.INE(iE,1)=ip2 - 1;
    eGridArray.INE(iE,2)=ip3 - 1;
  }
  IN.close();
  CHECK_UnstructuredGrid(eGridArray);
  return eGridArray;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray NC_ReadWW3_GridFile(std::string const& eFile)
{
  GridArray eGridArray;
  eGridArray.ModelName="WWM";
  eGridArray.IsFE=1;
  eGridArray.IsLonLat=1;
  eGridArray.L_IndexSelect=false;
  //
  eGridArray.INE=NC_ReadElements(eFile, "tri");
  MyVector<double> LON=NC_Read1Dvariable(eFile, "longitude");
  MyVector<double> LAT=NC_Read1Dvariable(eFile, "latitude");
  int nbPoint=LON.size();
  MyMatrix<double> LONarr(nbPoint,1);
  MyMatrix<double> LATarr(nbPoint,1);
  MyMatrix<double> DEParr(nbPoint,1);
  MyMatrix<double> ANGarr(nbPoint,1);
  MyMatrix<int> MSKarr(nbPoint,1);
  for (int iPoint=0; iPoint<nbPoint; iPoint++) {
    LONarr(iPoint,0)=LON(iPoint);
    LATarr(iPoint,0)=LAT(iPoint);
    DEParr(iPoint,0)=0;
    ANGarr(iPoint,0)=0;
    MSKarr(iPoint,0)=1;
  }
  eGridArray.GrdArrRho.LON=LONarr;
  eGridArray.GrdArrRho.LAT=LATarr;
  eGridArray.GrdArrRho.DEP=DEParr;
  eGridArray.GrdArrRho.ANG=ANGarr;
  eGridArray.GrdArrRho.MSK=MSKarr;
  CHECK_UnstructuredGrid(eGridArray);
  return eGridArray;
}
/* ------------------------------ */


/* ------------------------------ */
void CutWorldMap(GridArray & GrdArr)
{
  // We cut at -180 - eps
  double eps=1e-8;
  int nbPoint=GrdArr.GrdArrRho.LON.rows();
  double LonSplit=0;
  while(1) {
    double MinDist=2400;
    for (int iPoint=0; iPoint<nbPoint; iPoint++) {
      double eLon=GrdArr.GrdArrRho.LON(iPoint,0);
      if (eLon > 0)
	eLon=eLon - 360;
      LonSplit= - 180 - eps;
      double dist=fabs(eLon - LonSplit);
      if (dist < MinDist)
	MinDist=dist;
    }
    std::cerr << "eps=" << eps << " MinDist=" << MinDist << "\n";
    if (MinDist > eps)
      break;
    eps=eps*2;
  }
  int nbTrig=GrdArr.INE.rows();
  int ListStatus[nbTrig];
  int SumStatus=0;
  for (int iTrig=0; iTrig<nbTrig; iTrig++) {
    int i1=GrdArr.INE(iTrig,0);
    int i2=GrdArr.INE(iTrig,1);
    int i3=GrdArr.INE(iTrig,2);
    double eLon1=GrdArr.GrdArrRho.LON(i1,0);
    double eLon2=GrdArr.GrdArrRho.LON(i2,0);
    double eLon3=GrdArr.GrdArrRho.LON(i3,0);
    eLon1=eLon1 - LonSplit;
    eLon2=eLon2 - LonSplit;
    eLon3=eLon3 - LonSplit;
    DifferenceLonRenormalize(eLon1);
    DifferenceLonRenormalize(eLon2);
    DifferenceLonRenormalize(eLon3);
    int eStatus=1;
    double UpperLimit=90;
    if (fabs(eLon1) < UpperLimit &&
	fabs(eLon2) < UpperLimit &&
	fabs(eLon3) < UpperLimit) {
      double eProd12=eLon1*eLon2;
      double eProd23=eLon2*eLon3;
      double eProd31=eLon3*eLon1;
      if (eProd12 < 0 || eProd23 < 0 || eProd31 < 0) {
	eStatus=0;
      }
    }
    ListStatus[iTrig]=eStatus;
    SumStatus+=eStatus;
  }
  std::cerr << "SumStatus=" << SumStatus << " nbTrig=" << nbTrig << "\n";
  int nbTrigNew=0;
  for (int iTrig=0; iTrig<nbTrig; iTrig++)
    if (ListStatus[iTrig] == 1)
      nbTrigNew++;
  MyMatrix<int> INEnew(nbTrigNew,3);
  int iTrigNew=0;
  for (int iTrig=0; iTrig<nbTrig; iTrig++)
    if (ListStatus[iTrig] == 1) {
      int i1=GrdArr.INE(iTrig,0);
      int i2=GrdArr.INE(iTrig,1);
      int i3=GrdArr.INE(iTrig,2);
      INEnew(iTrigNew,0)=i1;
      INEnew(iTrigNew,1)=i2;
      INEnew(iTrigNew,2)=i3;
      iTrigNew++;
    }
  GrdArr.INE=INEnew;
}
/* ------------------------------ */


/* ------------------------------ */
void CUT_HigherLatitude(GridArray & GrdArr, double MinLatCut, double MaxLatCut)
{
  int mnp=GrdArr.GrdArrRho.LON.rows();
  int mne=GrdArr.INE.rows();
  int ListStatus[mnp];
  int Index[mnp];
  int RevIndex[mnp];
  int iNodeNew=0;
  std::cerr << "CUT_HigherLatitude, step 1\n";
  std::vector<int> I_IndexSelectOld;
  if (GrdArr.L_IndexSelect) {
    I_IndexSelectOld=GrdArr.I_IndexSelect;
  }
  else {
    for (int i=0; i<mnp; i++)
      I_IndexSelectOld.push_back(i);
  }
  std::cerr << "CUT_HigherLatitude, step 2\n";
  std::vector<int> I_IndexSelect;
  for (int iNode=0; iNode<mnp; iNode++) {
    double eLat=GrdArr.GrdArrRho.LAT(iNode,0);
    if (MinLatCut < eLat && eLat < MaxLatCut) {
      ListStatus[iNode]=1;
      Index[iNode]=iNodeNew;
      RevIndex[iNodeNew]=iNode;
      int iNodeMain=I_IndexSelectOld[iNode];
      I_IndexSelect.push_back(iNodeMain);
      iNodeNew++;
    }
    else {
      ListStatus[iNode]=0;
    }
  }
  std::cerr << "CUT_HigherLatitude, step 3\n";
  int nbNodeNew=iNodeNew;
  MyMatrix<double> LONnew(nbNodeNew,1);
  MyMatrix<double> LATnew(nbNodeNew,1);
  MyMatrix<double> DEPnew(nbNodeNew,1);
  for (int iNodeNew=0; iNodeNew<nbNodeNew; iNodeNew++) {
    int iNode=RevIndex[iNodeNew];
    double eLon=GrdArr.GrdArrRho.LON(iNode,0);
    double eLat=GrdArr.GrdArrRho.LAT(iNode,0);
    double eDep=GrdArr.GrdArrRho.DEP(iNode,0);
    LONnew(iNodeNew,0)=eLon;
    LATnew(iNodeNew,0)=eLat;
    DEPnew(iNodeNew,0)=eDep;
  }
  std::cerr << "CUT_HigherLatitude, step 4\n";
  GrdArr.GrdArrRho.LON=LONnew;
  GrdArr.GrdArrRho.LAT=LATnew;
  GrdArr.GrdArrRho.DEP=DEPnew;
  int nbTrigNew=0;
  for (int iTrig=0; iTrig<mne; iTrig++) {
    int i1=GrdArr.INE(iTrig,0);
    int i2=GrdArr.INE(iTrig,1);
    int i3=GrdArr.INE(iTrig,2);
    if (ListStatus[i1] == 1 && ListStatus[i2] == 1 && ListStatus[i3] == 1)
      nbTrigNew++;
  }
  std::cerr << "CUT_HigherLatitude, step 5\n";
  MyMatrix<int> INEnew(nbTrigNew, 3);
  int iTrigNew=0;
  for (int iTrig=0; iTrig<mne; iTrig++) {
    int i1=GrdArr.INE(iTrig,0);
    int i2=GrdArr.INE(iTrig,1);
    int i3=GrdArr.INE(iTrig,2);
    if (ListStatus[i1] == 1 && ListStatus[i2] == 1 && ListStatus[i3] == 1) {
      INEnew(iTrigNew,0)=Index[i1];
      INEnew(iTrigNew,1)=Index[i2];
      INEnew(iTrigNew,2)=Index[i3];
      iTrigNew++;
    }
  }
  std::cerr << "CUT_HigherLatitude, step 6\n";
  GrdArr.INE=INEnew;
  GrdArr.L_IndexSelect = true;
  GrdArr.I_IndexSelect = I_IndexSelect;
}
/* ------------------------------ */


/* ------------------------------ */
void CHECK_COORDINATE_ORIENTATION(GridArray const& GrdArr)
{
  int mne=GrdArr.INE.rows();
  int nbPlus=0;
  int nbMinus=0;
  for (int ie=0; ie<mne; ie++) {
    int i1=GrdArr.INE(ie, 0);
    int i2=GrdArr.INE(ie, 1);
    int i3=GrdArr.INE(ie, 2);
    double eLon1=GrdArr.GrdArrRho.LON(i1,0);
    double eLon2=GrdArr.GrdArrRho.LON(i2,0);
    double eLon3=GrdArr.GrdArrRho.LON(i3,0);
    double eLat1=GrdArr.GrdArrRho.LAT(i1,0);
    double eLat2=GrdArr.GrdArrRho.LAT(i2,0);
    double eLat3=GrdArr.GrdArrRho.LAT(i3,0);
    double deltaLON12=eLon2 - eLon1;
    double deltaLAT12=eLat2 - eLat1;
    double deltaLON13=eLon3 - eLon1;
    double deltaLAT13=eLat3 - eLat1;
    double eArea=deltaLON13 * deltaLAT12 - deltaLON12 * deltaLAT13;
    if (eArea > 0) {
      nbPlus++;
    }
    else {
      nbMinus++;
    }
  }
  if (nbPlus > 0 && nbMinus > 0) {
    std::cerr << "Orientation error\n";
    std::cerr << "nbPlus =" << nbPlus << "\n";
    std::cerr << "nbMinus=" << nbMinus << "\n";
    exit(1);
  }
  std::cerr << "nbPlus = " << nbPlus << "  nbMinus = " << nbMinus << "\n";
}
/* ------------------------------ */


/* ------------------------------ */
struct TripleModelDesc {
  std::string ModelName;
  std::string GridFile;
  std::string HisPrefix;
  bool CutWorldMap;
  bool HigherLatitudeCut;
  double MinLatCut;
  double MaxLatCut;
};
/* ------------------------------ */


/* ------------------------------ */
std::string GET_GRID_FILE(TripleModelDesc const& eTriple)
{
  std::string eModelName=eTriple.ModelName;
  std::string HisPrefix=eTriple.HisPrefix;
  if (eModelName == "COSMO")
    return HisPrefix + "0001.nc";
  if (eModelName == "WAM")
    return HisPrefix + "0001.nc";
  if (eModelName == "ROMS")
    return eTriple.GridFile;
  if (eModelName == "WWM")
    return eTriple.GridFile;
  if (eModelName == "WW3") {
    std::string ThePrefix=HisPrefix + "*";
    std::vector<std::string> ListFile=ls_operation(ThePrefix);
    return ListFile[0];
  }
  if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO") {
    std::vector<std::string> ListFile=GRIB_GetAllFilesInDirectory(HisPrefix);
    if (ListFile.size() == 0) {
      std::cerr << "The list of files is empty\n";
      std::cerr << "Error happened in GRIB_GetAllFilesInDirectory\n";
      exit(1);
    }
    return ListFile[0];
  }
  std::cerr << "Error in GET_GRID_FILE\n";
  std::cerr << "Did not find the matching model for the grid\n";
  std::cerr << "Please correct\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
std::string GetStrAllOfPlot(VarQuery const& eQuery)
{
  int iTime=eQuery.iTime;
  std::string strAll;
  std::string strFile=DATE_ConvertMjd2mystringFile(eQuery.eTimeDay);
  if (iTime == -1) {
    strAll=strFile;
  }
  else {
    if (iTime > 10000) {
      std::cerr << "Error in the code\n";
      std::cerr << "iTime is too large\n";
      exit(1);
    }
    strAll=StringNumber(iTime, 4) + "_" + strFile;
  }
  if (eQuery.NatureQuery != "instant")
    strAll=strAll + "_" + eQuery.NatureQuery;
  return strAll;
}
/* ------------------------------ */


/* ------------------------------ */
struct DrawArr {
  PermanentInfoDrawing ePerm;
  bool DoTitle;
  std::string TitleStr;
  bool DrawRiver = false;
  bool DrawContourBathy = false;
  // colormaps
  std::string ColorMap = "BlWhRe";
  std::string cnFillMode = "RasterFill";
  bool cnSmoothingOn = true;
  int nbLevelSpa = 40;
  int nbLabelStride = 10;
  // annotations
  AnnotationRec TheAnnot;
  bool FillLand = true;
  bool UseNativeGrid = true;
  std::string GridResolution = "HighRes";
  std::vector<SeqLineSegment> ListLineSegment;
};
/* ------------------------------ */


/* ------------------------------ */
void CALL_NCL(bool const& KeepNC_NCL, 
	      std::string const& TargetFile, 
	      std::string const& eFileNC, 
	      std::string const& eFileNCL)
{
  //
  std::string eComm="ncl";
  std::string eCommand=eComm + " " + eFileNCL + " > /dev/null";
  int iret=system(eCommand.c_str());
  //  std::cerr << "iret=" << iret << "\n";
  if (iret == -1) {
    printf("Oh dear, something went wrong with ncl! %s\n", strerror(errno));
    exit(1);
  }
  /*  if (iret == -1) {
    std::cerr << "Failed to execute the NCL program\n";
    exit(1);
    }*/
  if (IsExistingFile(TargetFile) == false) {
    std::cerr << "The following TargetFile was not created\n";
    std::cerr << "TargetFile=" << TargetFile << "\n";
    std::cerr << "eFileNC = " << eFileNC << "\n";
    std::cerr << "eFileNCL = " << eFileNCL << "\n";
    std::cerr << "Please debug\n";
    exit(1);
  }
  if (KeepNC_NCL == false) {
    RemoveFileIfExist(eFileNC);
    RemoveFileIfExist(eFileNCL);
  }
}
/* ------------------------------ */


/* ------------------------------ */
void ADD_RIVER(std::ofstream & os, DrawArr const& eDrawArr)
{
  if (eDrawArr.DrawRiver) {
    os << "  riv_data = asciiread(\"JadranRivers_extractor.dat\",(/10583,2/),\"float\")\n";
    os << "  lon=riv_data(:,0)\n";
    os << "  lat=riv_data(:,1)\n";
    os << "  segments=ind(lon.eq.-999)\n";
    os << "  ns=dimsizes(segments)\n";
    os << "  resP = True\n";
    os << "  resP@gsLineThicknessF = 1.5\n";
    os << "  resP@gsLineColor  = \"dodgerblue1\"\n";
    os << "  resP@tfPolyDrawOrder = \"PostDraw\"\n";
    os << "  lines = new(ns(0)-1,graphic)   ; array to hold polylines\n";
    os << "  do i=0,ns-2\n";
    os << "    xp=lon(segments(i)+2 : segments(i+1)-2)\n";
    os << "    yp=lat(segments(i)+2 : segments(i+1)-2)\n";
    os << "    lines(i)=gsn_add_polyline(wks,plot,xp,yp,resP)\n";
    os << "    delete(xp)\n";
    os << "    delete(yp)\n";
    os << "  end do\n";
    os << "  delete(segments)\n";
  }
}
/* ------------------------------ */


/* ------------------------------ */
struct DrawScatterArr {
  std::string VarNameAB_file;
  PermanentInfoDrawing ePerm;
  // Drawing itself
  AnnotationRec TheAnnot;
  bool DoTitle;
  bool AddStatMeasModel;
  std::string NameA_plot;
  std::string NameB_plot;
  std::vector<double> data_rangeA;
  std::vector<double> data_rangeB;
  MyVector<double> eVectA;
  MyVector<double> eVectB;
  int aSize;
  int bSize;
};
/* ------------------------------ */


/* ------------------------------ */
void DEFINE_SCATTER_NC(std::string const& eFileNC, 
		       DrawScatterArr const& eDrawScatter)
{
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  int nb=eDrawScatter.eVectA.size();
  int aSize=eDrawScatter.aSize;
  int bSize=eDrawScatter.bSize;
  netCDF::NcDim eDimOne=dataFile.addDim("one", 1);
  netCDF::NcDim eDimTwo=dataFile.addDim("two", 2);
  netCDF::NcDim eDimNb=dataFile.addDim("nb", nb);
  netCDF::NcDim eDimASize=dataFile.addDim("aSize", aSize);
  netCDF::NcDim eDimBSize=dataFile.addDim("bSize", bSize);

  std::vector<std::string> ListDimOne={"one"};
  std::vector<std::string> ListDimTwo={"two"};
  std::vector<std::string> ListDimNb={"nb"};
  std::vector<std::string> ListDimABsize={"aSize", "bSize"};

  netCDF::NcVar eVarData_rangeA=dataFile.addVar("data_rangeA", "double", ListDimTwo);
  netCDF::NcVar eVarData_rangeB=dataFile.addVar("data_rangeB", "double", ListDimTwo);
  netCDF::NcVar eVarListX=dataFile.addVar("ListX", "double", ListDimTwo);
  netCDF::NcVar eVarListY=dataFile.addVar("ListY", "double", ListDimTwo);
  netCDF::NcVar eVarX=dataFile.addVar("X", "double", ListDimNb);
  netCDF::NcVar eVarY=dataFile.addVar("Y", "double", ListDimNb);
  netCDF::NcVar eVarX2D=dataFile.addVar("X2D", "double", ListDimABsize);
  netCDF::NcVar eVarY2D=dataFile.addVar("Y2D", "double", ListDimABsize);
  netCDF::NcVar eVarCanvas=dataFile.addVar("canvas", "double", ListDimABsize);
  //
  double *eFieldX, *eFieldY;
  eFieldX=new double[nb];
  eFieldY=new double[nb];
  double TheMaxX=0;
  double TheMaxY=0;
  for (int i=0; i<nb; i++) {
    double eVal=eDrawScatter.eVectA(i);
    if (eVal > TheMaxX)
      TheMaxX=eVal;
    eFieldX[i]=eVal;
  }
  for (int i=0; i<nb; i++) {
    double eVal=eDrawScatter.eVectB(i);
    if (eVal > TheMaxY)
      TheMaxY=eVal;
    eFieldY[i]=eVal;
  }
  eVarX.putVar(eFieldX);
  eVarY.putVar(eFieldY);
  delete [] eFieldX;
  delete [] eFieldY;
  //
  double ePair[2];
  ePair[0]=0;
  ePair[1]=TheMaxX;
  eVarListX.putVar(ePair);
  ePair[0]=0;
  ePair[1]=TheMaxY;
  eVarListY.putVar(ePair);
  //
  std::vector<double> data_rangeA=eDrawScatter.data_rangeA;
  std::vector<double> data_rangeB=eDrawScatter.data_rangeB;
  double eFrangeA[2];
  for (int i=0; i<2; i++)
    eFrangeA[i]=data_rangeA[i];
  eVarData_rangeA.putVar(eFrangeA);
  double eFrangeB[2];
  for (int i=0; i<2; i++)
    eFrangeB[i]=data_rangeB[i];
  eVarData_rangeB.putVar(eFrangeB);
  //
  double *X2D, *Y2D;
  X2D=new double[aSize*bSize];
  Y2D=new double[aSize*bSize];
  double deltaA=(data_rangeA[1] - data_rangeA[0])/double(aSize - 1);
  double deltaB=(data_rangeB[1] - data_rangeB[0])/double(bSize - 1);
  int idx=0;
  for (int iA=0; iA<aSize; iA++)
    for (int iB=0; iB<bSize; iB++) {
      double eA=data_rangeA[0] + iA*deltaA;
      double eB=data_rangeB[0] + iB*deltaB;
      X2D[idx]=eA;
      Y2D[idx]=eB;
      idx++;
    }
  eVarX2D.putVar(X2D);
  eVarY2D.putVar(Y2D);
  delete [] X2D;
  delete [] Y2D;
  //
  int canvasInt[aSize][bSize];
  for (int iA=0; iA<aSize; iA++)
    for (int iB=0; iB<bSize; iB++)
      canvasInt[iA][iB]=0;
  for (int iEnt=0; iEnt<nb; iEnt++) {
    double eA=eDrawScatter.eVectA(iEnt);
    double eB=eDrawScatter.eVectB(iEnt);
    double iA_d=(eA - data_rangeA[0])/deltaA;
    double iB_d=(eB - data_rangeB[0])/deltaB;
    int iA=floor(iA_d);
    int iB=floor(iB_d);
    if (iA >=0 && iA<aSize && iB>=0 && iB<bSize)
      canvasInt[iA][iB]++;
  }
  double *canvas;
  canvas=new double[aSize*bSize];
  double MissVal=0;
  idx=0;
  for (int iA=0; iA<aSize; iA++)
    for (int iB=0; iB<bSize; iB++) {
      double eVal;
      int eValI=canvasInt[iA][iB];
      if (eValI == 0)
	eVal=MissVal;
      else
	eVal=log10(double(eValI));
      //      std::cerr << "iA=" << iA << " iB=" << iB << " eVal=" << eVal << "\n";
      canvas[idx]=eVal;
      idx++;
    }
  eVarCanvas.putVar(canvas);
  delete [] canvas;
  //
}
/* ------------------------------ */


/* ------------------------------ */
void DEFINE_QUIVER_NC(std::string const& eFileNC, 
		      GridArray const& GrdArr, 
		      MyMatrix<double> const& U_rho, MyMatrix<double> const& V_rho)
{
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  int eta=GrdArr.GrdArrRho.LON.rows();
  int xi =GrdArr.GrdArrRho.LON.cols();
  int eta_lat=GrdArr.GrdArrRho.LAT.rows();
  int xi_lat =GrdArr.GrdArrRho.LAT.cols();
  int eta_U=U_rho.rows();
  int xi_U =U_rho.cols();
  int eta_V=V_rho.rows();
  int xi_V =V_rho.cols();
  if (eta != eta_lat || eta != eta_U || eta != eta_V) {
    std::cerr << "All dimensions should be the same.\n";
    std::cerr << "Now we have following:\n";
    std::cerr << "eta_lon=" << eta << "\n";
    std::cerr << "eta_lat=" << eta_lat << "\n";
    std::cerr << "eta_U=" << eta_U << "\n";
    std::cerr << "eta_V=" << eta_V << "\n";
    exit(1);
  }
  if (xi != xi_lat || xi != xi_U || xi != xi_V) {
    std::cerr << "All dimensions should be the same.\n";
    std::cerr << "Now we have following:\n";
    std::cerr << "xi_lon=" << xi << "\n";
    std::cerr << "xi_lat=" << xi_lat << "\n";
    std::cerr << "xi_U=" << xi_U << "\n";
    std::cerr << "xi_V=" << xi_V << "\n";
    exit(1);
  }
  std::string typeName="double";
  std::string typeNameInt="int";
  std::string eEta="eta_rho";
  std::string eXi ="xi_rho";
  std::string Lon="lon";
  std::string Lat="lat";
  std::string Uvar="u";
  std::string Vvar="v";
  if (xi > 1) {
    netCDF::NcDim eDimEta=dataFile.addDim(eEta, eta);
    netCDF::NcDim eDimXi =dataFile.addDim(eXi, xi);
    std::vector<std::string> ListDim={eEta, eXi};
    netCDF::NcVar eVarLON=dataFile.addVar(Lon, typeName, ListDim);
    netCDF::NcVar eVarLAT=dataFile.addVar(Lat, typeName, ListDim);
    netCDF::NcVar eVarU=dataFile.addVar(Uvar, typeName, ListDim);
    netCDF::NcVar eVarV=dataFile.addVar(Vvar, typeName, ListDim);
    double *valLON, *valLAT, *valU, *valV;
    valLON=new double[eta*xi];
    valLAT=new double[eta*xi];
    valU  =new double[eta*xi];
    valV  =new double[eta*xi];
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	valLON[idx]=GrdArr.GrdArrRho.LON(i, j);
	valLAT[idx]=GrdArr.GrdArrRho.LAT(i, j);
	valU[idx]=U_rho(i,j);
	valV[idx]=V_rho(i,j);
	idx++;
      }
    eVarLON.putVar(valLON);
    eVarLAT.putVar(valLAT);
    eVarU.putVar(valU);
    eVarV.putVar(valV);
    delete [] valLON;
    delete [] valLAT;
    delete [] valU;
    delete [] valV;
  }
  else {
    int mnp=eta;
    int mne=GrdArr.INE.rows();
    std::string eMnp="mnp";
    std::string eMne="mne";
    std::string eThree="three";
    netCDF::NcDim eDimMnp=dataFile.addDim(eMnp, mnp);
    netCDF::NcDim eDimThree=dataFile.addDim(eThree, 3);
    netCDF::NcDim eDimMne=dataFile.addDim(eMne, mne);
    std::vector<std::string> ListDim={eMnp};
    netCDF::NcVar eVarLON=dataFile.addVar(Lon, typeName, ListDim);
    netCDF::NcVar eVarLAT=dataFile.addVar(Lat, typeName, ListDim);
    netCDF::NcVar eVarU=dataFile.addVar(Uvar, typeName, ListDim);
    netCDF::NcVar eVarV=dataFile.addVar(Vvar, typeName, ListDim);
    std::string Fine="ele";
    std::vector<std::string> ListDimINE={eMne, eThree};
    netCDF::NcVar eVarINE=dataFile.addVar(Fine, typeNameInt, ListDimINE);
    double *valLON, *valLAT, *valU, *valV;
    valLON=new double[mnp];
    valLAT=new double[mnp];
    valU  =new double[mnp];
    valV  =new double[mnp];
    for (int i=0; i<eta; i++) {
      valLON[i]=GrdArr.GrdArrRho.LON(i,0);
      valLAT[i]=GrdArr.GrdArrRho.LAT(i,0);
      valU[i]=U_rho(i,0);
      valV[i]=V_rho(i,0);
    }
    eVarLON.putVar(valLON);
    eVarLAT.putVar(valLAT);
    eVarU.putVar(valU);
    eVarV.putVar(valV);
    delete [] valLON;
    delete [] valLAT;
    delete [] valU;
    delete [] valV;
    //
    int *valI;
    valI=new int[3*mne];
    int idx=0;
    for (int ie=0; ie<mne; ie++)
      for (int i=0; i<3; i++) {
	int eConn=GrdArr.INE(ie,i);
	valI[idx]=eConn;
	idx++;
      }
    eVarINE.putVar(valI);
    delete [] valI;
  }
}
/* ------------------------------ */


/* ------------------------------ */
void PLOT_QUIVER(std::string const& FileName, 
		 GridArray const& GrdArr, 
		 DrawArr const& eDrawArr, 
		 RecVar const& eRecVar)
{
  std::string eFileNC=eDrawArr.ePerm.PrefixTemp + "DataQuiver_" + eRecVar.VarName1 + "_" + eRecVar.strAll + ".nc";
  std::string eFileNCL=eDrawArr.ePerm.PrefixTemp + "ScriptQuiver_" + eRecVar.VarName1 + "_" + eRecVar.strAll + ".ncl";
  std::string TargetFile=FileName + "." + eDrawArr.ePerm.Extension;
  DEFINE_QUIVER_NC(eFileNC, GrdArr, eRecVar.U, eRecVar.V);
  //
  std::ofstream OUTncl;
  OUTncl.open(eFileNCL);
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
  OUTncl << "begin\n";
  OUTncl << "  ; \n";
  OUTncl << "  ; Loading fundamental data\n";
  OUTncl << "  ; \n";
  OUTncl << "  f = addfile(\"" << eFileNC << "\", \"r\")\n";
  OUTncl << "  lat2d  = f->lat\n";
  OUTncl << "  lon2d  = f->lon\n";
  OUTncl << "  eVarU  = f->u\n";
  OUTncl << "  eVarV  = f->v\n";
  OUTncl << "  speed=sqrt(eVarU^2 + eVarV^2)\n";
  OUTncl << "  wks  = gsn_open_wks (\"" << eDrawArr.ePerm.Extension << "\",\"" << FileName << "\")\n";
  OUTncl << "  vres1=True\n";
  OUTncl << "  vres1@gsnDraw = False\n";
  OUTncl << "  vres1@gsnFrame = False\n";
  OUTncl << "  vres1@gsnPaperOrientation  = \"Portrait\"\n";
  OUTncl << "  vres1@gsnScalarContour     = True\n";
  OUTncl << "  vres1@gsnSpreadColors      = True\n";
  OUTncl << "  vres1@gsnSpreadColorEnd    = -3\n";
  OUTncl << "  ;  vres1@gsnSpreadColorStart  = 1\n";
  OUTncl << "  vres1@mpLandFillColor      = \"grey\"\n";
  OUTncl << "  vres1@cnFillDrawOrder      = \"PreDraw\"\n";
  OUTncl << "  vres1@cnFillOn             = True\n";
  OUTncl << "  vres1@cnLinesOn            = False\n";
  OUTncl << "  vres1@cnLineLabelsOn       = False\n";
  OUTncl << "  vres1@cnLevelSelectionMode = \"ManualLevels\"\n";
  OUTncl << "  vres1@cnMinLevelValF       = " << eRecVar.minval << "\n";
  OUTncl << "  vres1@cnMaxLevelValF       = " << eRecVar.maxval << "\n";
  OUTncl << "  vres1@lbLabelStride        = 8\n";
  int nbLevelSpa=eDrawArr.nbLevelSpa;
  double TheLevelSpa=(eRecVar.maxval - eRecVar.minval)/double(nbLevelSpa);
  OUTncl << "  vres1@cnLevelSpacingF      = " << TheLevelSpa << "\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; The vertical label bar on the right\n";
  OUTncl << "  ;\n";
  OUTncl << "  vres1@lbOrientation        = \"Vertical\"     ; Vertical label bar\n";
  OUTncl << "  vres1@pmLabelBarOrthogonalPosF = 0.025          ; move label bar closer\n";
  OUTncl << "  vres1@pmLabelBarDisplayMode = \"Always\"          ; Turn on a label bar.\n";
  OUTncl << "  vres1@lbPerimOn             = False             ; no box around it\n";
  OUTncl << "  vres1@lbBoxLinesOn         = False               ; Yes/No labelbar box lines.\n";
  OUTncl << "  vres1@lbTitleString    = \"" << eRecVar.Unit << "\"\n";
  OUTncl << "  vres1@lbTitleFont      = \"Helvetica\"\n";
  OUTncl << "  vres1@lbTitleFontHeightF = 0.015\n";
  OUTncl << "  vres1@lbTitleDirection     = \"Across\" \n";
  OUTncl << "  vres1@lbTitlePosition = \"Right\"\n";
  OUTncl << "  vres1@lbTitleAngleF = 90\n";
  OUTncl << "  vres1@mpProjection = \"Mercator\"\n";
  OUTncl << "  vres1@mpLimitMode         = \"Corners\"             ; choose range of map\n";
  OUTncl << "  vres1@mpLeftCornerLatF    = min(lat2d)\n";
  OUTncl << "  vres1@mpLeftCornerLonF    = min(lon2d)\n";
  OUTncl << "  vres1@mpRightCornerLatF   = max(lat2d)\n";
  OUTncl << "  vres1@mpRightCornerLonF   = max(lon2d)\n";
  OUTncl << "  vres1@mpDataBaseVersion      = \"" << eDrawArr.GridResolution << "\"          ; use high resolution coast\n";
  OUTncl << "  vres1@pmTickMarkDisplayMode  = \"Always\"           ; turn on tickmarks\n";
  OUTncl << "  vres1@mpLandFillColor      = \"LightGrey\"\n";
  OUTncl << "  vres1@vcMonoLineArrowColor  = True             ; colored by their mag\n";
  OUTncl << "  vres1@pmLabelBarWidthF = 0.02\n";
  if (eDrawArr.DoTitle) {
    OUTncl << "  vres1@tiMainString    = \"" << eDrawArr.TitleStr << "\"\n";
    OUTncl << "  vres1@tiMainFont      = \"Helvetica\"\n";
    OUTncl << "  vres1@tiMainFontHeightF=0.015\n";
  }
  OUTncl << "  ; vres1@gsnRightString    = \"Sea surface elevation\"\n";
  OUTncl << "  ; vres1@gsnLeftString     = \"Difference\"\n";
  OUTncl << "  ; vres1@gsnRightString    = \"\"\n";
  OUTncl << "  ; vres1@gsnLeftString     = \"\"\n";
  OUTncl << "  ;  vres1@vcRefMagnitudeF   = 10.\n";
  OUTncl << "  vres1@vcGlyphStyle      = \"CurlyVector\"\n";
  OUTncl << "  ;  vres1@vcGlyphStyle      = \"LineArrow\"\n";
  OUTncl << "  vres1@vcMinDistanceF    = 0.014\n";
  OUTncl << "  vres1@vcRefLengthF      = 0.02\n";
  OUTncl << "  vres1@vcRefAnnoOn       = False\n";
  OUTncl << "  vres1@vcLineArrowHeadMaxSizeF = 0.005\n";
  OUTncl << "  ;  vres1@vcRefAnnoOrthogonalPosF = -0.8\n";
  OUTncl << "  ;  vres1@vcRefAnnoParallelPosF = 1.0\n";
  OUTncl << "  ;  vres1@vcRefMagnitudeF         = 10\n";
  OUTncl << "  ;  vres1@vcRefAnnoString1 = \"10 m/s\"\n";
  OUTncl << "  ;  gsn_define_colormap (wks,\"testcmap\")\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; Specifying the colormap\n";
  OUTncl << "  ;\n";
  OUTncl << "  gsn_define_colormap (wks,\"BlAqGrYeOrReVi200\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"BlGrYeOrReVi200\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"hotres\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"rainbow\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"ViBlGrWhYeOrRe\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"BlWhRe\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"GrayWhiteGray\")\n";
  OUTncl << "  ;  gsn_define_colormap(wks,\"BlWhRe\")\n";
  OUTncl << "  i = NhlNewColor(wks,0.8,0.8,0.8)      ; add gray to colormap\n";
  OUTncl << "  i = NhlNewColor(wks,0.9,0.9,0.9)      ; add gray to colormap\n";
  OUTncl << "  eVarU@lat2d=lat2d\n";
  OUTncl << "  eVarU@lon2d=lon2d\n";
  OUTncl << "  eVarV@lat2d=lat2d\n";
  OUTncl << "  eVarV@lon2d=lon2d\n";
  OUTncl << "  speed@lat2d=lat2d\n";
  OUTncl << "  speed@lon2d=lon2d\n";
  bool UseFvar=false;
  if (UseFvar) {
    OUTncl << "  fVarU=eVarU/speed\n";
    OUTncl << "  fVarV=eVarV/speed\n";
    OUTncl << "  fVarU@lat2d=lat2d\n";
    OUTncl << "  fVarU@lon2d=lon2d\n";
    OUTncl << "  fVarV@lat2d=lat2d\n";
    OUTncl << "  fVarV@lon2d=lon2d\n";
    OUTncl << "  plot = gsn_csm_vector_scalar_map(wks,fVarU,fVarV,speed,vres1)\n";
  }
  else {
    OUTncl << "  plot = gsn_csm_vector_scalar_map(wks,eVarU,eVarV,speed,vres1)\n";
  }
  ADD_RIVER(OUTncl, eDrawArr);
  ADD_ANNOTATION_TEXT(OUTncl, eDrawArr.TheAnnot);
  OUTncl << "  draw(plot)\n";
  OUTncl << "  frame(wks)\n";
  OUTncl << "end\n";
  OUTncl.close();
  //
  CALL_NCL(eDrawArr.ePerm.KeepNC_NCL, 
	   TargetFile, eFileNC, eFileNCL);
}
/* ------------------------------ */


/* ------------------------------ */
void DEFINE_PCOLOR_NC_NCL(std::string const& eFileNC, 
			  GridArray const& GrdArr, 
			  MyMatrix<double> const& F_rho, 
			  bool const& WriteDEP,
			  std::vector<SeqLineSegment> const& ListLineSegment)
{
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  if (GrdArr.IsFE == 0) {
    int eta=GrdArr.GrdArrRho.LON.rows();
    int xi =GrdArr.GrdArrRho.LON.cols();
    netCDF::NcDim eDimEta=dataFile.addDim("eta_rho", eta);
    netCDF::NcDim eDimXi =dataFile.addDim("xi_rho", xi);
    std::vector<std::string> ListDim={"eta_rho", "xi_rho"};
    netCDF::NcVar eVarLON=dataFile.addVar("lon", "double", ListDim);
    netCDF::NcVar eVarLAT=dataFile.addVar("lat", "double", ListDim);
    netCDF::NcVar eVarF=dataFile.addVar("field", "double", ListDim);
    double *valLON, *valLAT, *valF;
    valLON=new double[eta*xi];
    valLAT=new double[eta*xi];
    valF  =new double[eta*xi];
    int idx=0;
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	valLON[idx]=GrdArr.GrdArrRho.LON(i, j);
	valLAT[idx]=GrdArr.GrdArrRho.LAT(i, j);
	double eValF=F_rho(i,j);
	valF[idx]=eValF;
	idx++;
      }
    eVarLON.putVar(valLON);
    eVarLAT.putVar(valLAT);
    eVarF.putVar(valF);
    delete [] valLON;
    delete [] valLAT;
    delete [] valF;
    if (WriteDEP) {
      double *valD;
      valD=new double[eta*xi];
      netCDF::NcVar eVarDEP=dataFile.addVar("dep", "double", ListDim);
      int idx=0;
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++) {
	  double eValD=GrdArr.GrdArrRho.DEP(i, j);
	  valD[idx]=eValD;
	  idx++;
	}
      eVarDEP.putVar(valD);
      delete [] valD;
    }
  }
  else {
    int mnp=GrdArr.GrdArrRho.LON.rows();
    int mne=GrdArr.INE.rows();
    netCDF::NcDim eDimMnp=dataFile.addDim("mnp", mnp);
    netCDF::NcDim eDimThree=dataFile.addDim("three", 3);
    netCDF::NcDim eDimMne=dataFile.addDim("mne", mne);
    std::vector<std::string> ListDim={"mnp"};
    netCDF::NcVar eVarLON=dataFile.addVar("lon", "double", ListDim);
    netCDF::NcVar eVarLAT=dataFile.addVar("lat", "double", ListDim);
    netCDF::NcVar eVarF=dataFile.addVar("field", "double", ListDim);
    std::vector<std::string> ListDimINE={"mne", "three"};
    netCDF::NcVar eVarINE=dataFile.addVar("ele", "int", ListDimINE);
    double *valLON, *valLAT, *valF;
    valLON=new double[mnp];
    valLAT=new double[mnp];
    valF  =new double[mnp];
    for (int i=0; i<mnp; i++) {
      valLON[i]=GrdArr.GrdArrRho.LON(i,0);
      valLAT[i]=GrdArr.GrdArrRho.LAT(i,0);
      double eValF=F_rho(i,0);
      valF[i]=eValF;
    }
    eVarLON.putVar(valLON);
    eVarLAT.putVar(valLAT);
    eVarF.putVar(valF);
    delete [] valLON;
    delete [] valLAT;
    delete [] valF;
    //
    int *valI;
    valI=new int[3*mne];
    int idx=0;
    for (int ie=0; ie<mne; ie++)
      for (int i=0; i<3; i++) {
	int eConn=GrdArr.INE(ie,i);
	valI[idx]=eConn;
	idx++;
      }
    eVarINE.putVar(valI);
    delete [] valI;
    //
    if (WriteDEP) {
      netCDF::NcVar eVarDEP=dataFile.addVar("dep", "double", ListDim);
      double *valD;
      valD=new double[mnp];
      for (int i=0; i<mnp; i++) {
	double eValD=GrdArr.GrdArrRho.DEP(i,0);
	valD[i]=eValD;
      }
      eVarDEP.putVar(valD);
      delete [] valD;
    }
  }
  int nbLineSeq=ListLineSegment.size();
  if (nbLineSeq > 0) {
    int TotalLen=0;
    for (auto& eSeqLineSegment : ListLineSegment) {
      int len=eSeqLineSegment.ListPairLL.size();
      if (eSeqLineSegment.IsClosed == false) {
	len--;
      }
      TotalLen += len;
    }
    int RelSiz=2*TotalLen;
    double *ListLon, *ListLat;
    ListLon=new double[RelSiz];
    ListLat=new double[RelSiz];
    int idx=0;
    for (auto& eSeqLineSegment : ListLineSegment) {
      int len=eSeqLineSegment.ListPairLL.size();
      for (int i=0; i<len-1; i++) {
	double eLon1=eSeqLineSegment.ListPairLL[i].eLon;
	double eLat1=eSeqLineSegment.ListPairLL[i].eLat;
	double eLon2=eSeqLineSegment.ListPairLL[i+1].eLon;
	double eLat2=eSeqLineSegment.ListPairLL[i+1].eLat;
	ListLon[idx]=eLon1;
	ListLat[idx]=eLat1;
	idx++;
	ListLon[idx]=eLon2;
	ListLat[idx]=eLat2;
	idx++;
      }
      if (eSeqLineSegment.IsClosed) {
	double eLon1=eSeqLineSegment.ListPairLL[0].eLon;
	double eLat1=eSeqLineSegment.ListPairLL[0].eLat;
	double eLon2=eSeqLineSegment.ListPairLL[len-1].eLon;
	double eLat2=eSeqLineSegment.ListPairLL[len-1].eLat;
	ListLon[idx]=eLon1;
	ListLat[idx]=eLat1;
	idx++;
	ListLon[idx]=eLon2;
	ListLat[idx]=eLat2;
	idx++;
      }
    }
    netCDF::NcDim eDimMnp=dataFile.addDim("RelSiz", RelSiz);
    std::vector<std::string> ListDim={"RelSiz"};
    netCDF::NcVar eVarListLon=dataFile.addVar("ListLon", "double", ListDim);
    netCDF::NcVar eVarListLat=dataFile.addVar("ListLat", "double", ListDim);
    eVarListLon.putVar(ListLon);
    eVarListLat.putVar(ListLat);
    delete [] ListLon;
    delete [] ListLat;
  }
}
/* ------------------------------ */


/* ------------------------------ */
void DEFINE_MESH_NC(std::string const& eFileNC, 
		    GridArray const& GrdArr)
{
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  if (GrdArr.IsFE == 1) {
    int mnp=GrdArr.GrdArrRho.LON.rows();
    int mne=GrdArr.INE.rows();
    MyMatrix<int> ListEdges=GetEdgeSet(GrdArr.INE, mnp);
    int nbEdge=ListEdges.rows();
    netCDF::NcDim eDimMnp=dataFile.addDim("mnp", mnp);
    netCDF::NcDim eDimThree=dataFile.addDim("three", 3);
    netCDF::NcDim eDimTwo=dataFile.addDim("two", 2);
    netCDF::NcDim eDimMne=dataFile.addDim("mne", mne);
    netCDF::NcDim eDimEdge=dataFile.addDim("nbEdges", nbEdge);
    netCDF::NcVar eVarLON=dataFile.addVar("lon", "double", std::vector<std::string>({"mnp"}));
    netCDF::NcVar eVarLAT=dataFile.addVar("lat", "double", std::vector<std::string>({"mnp"}));
    netCDF::NcVar eVarINE=dataFile.addVar("ele", "int", std::vector<std::string>({"mne", "three"}));
    netCDF::NcVar eVarEDGE=dataFile.addVar("edges", "int", std::vector<std::string>({"nbEdges", "two"}));
    double *valLON, *valLAT;
    valLON=new double[mnp];
    valLAT=new double[mnp];
    for (int i=0; i<mnp; i++) {
      valLON[i]=GrdArr.GrdArrRho.LON(i,0);
      valLAT[i]=GrdArr.GrdArrRho.LAT(i,0);
    }
    eVarLON.putVar(valLON);
    eVarLAT.putVar(valLAT);
    delete [] valLON;
    delete [] valLAT;
    //
    int *valI;
    valI=new int[mne*3];
    int idx=0;
    for (int ie=0; ie<mne; ie++)
      for (int i=0; i<3; i++) {
	int eConn=GrdArr.INE(ie,i);
	valI[idx]=eConn;
	idx++;
      }
    eVarINE.putVar(valI);
    delete [] valI;
    //
    int *valEDGE;
    valEDGE=new int[2*nbEdge];
    idx=0;
    for (int iedge=0; iedge<nbEdge; iedge++)
      for (int i=0; i<2; i++) {
	int eConn=ListEdges(iedge,i);
	valEDGE[idx]=eConn;
	idx++;
      }
    eVarEDGE.putVar(valEDGE);
    delete [] valEDGE;
  }
  else {
    std::cerr << "The corresponding code for finite differences need to be written\n";
    exit(1);
  }
}
/* ------------------------------ */


/* ------------------------------ */
void DEFINE_MESH_SVG(std::string const& eFileSVG, 
		     GridArray const& GrdArr)
{
  if (GrdArr.IsFE == 1) {
    int mnp=GrdArr.GrdArrRho.LON.rows();
    MyMatrix<int> ListEdges=GetEdgeSet(GrdArr.INE, mnp);
    int nbEdge=ListEdges.rows();
    //
    std::ofstream OUTsvg;
    OUTsvg.open(eFileSVG);
    OUTsvg << "<svg height=\"180\" width=\"360\">\n";
    for (int iEdge=0; iEdge<nbEdge; iEdge++) {
      int i1=ListEdges(iEdge,0);
      int i2=ListEdges(iEdge,1);
      double x1=180 + GrdArr.GrdArrRho.LON(i1);
      double x2=180 + GrdArr.GrdArrRho.LON(i2);
      double y1=90 - GrdArr.GrdArrRho.LAT(i1);
      double y2=90 - GrdArr.GrdArrRho.LAT(i2);
      OUTsvg << "  <line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" style=\"stroke:rgb(255,0,0);stroke-width:2\" />\n";
    }
    OUTsvg << "</svg>\n";
    OUTsvg.close();
    //
  }
  else {
    std::cerr << "The corresponding code for finite differences need to be written\n";
    exit(1);
  }
}
/* ------------------------------ */


/* ------------------------------ */
void PLOT_FD_PCOLOR(std::string const& FileName, 
		    GridArray const& GrdArr, 
		    DrawArr const& eDrawArr, 
		    RecVar const& eRecVar)
{
  std::string eFileNC=eDrawArr.ePerm.PrefixTemp + "DataPcolor_" + eRecVar.VarName1 + "_" + eRecVar.strAll + ".nc";
  std::string eFileNCL=eDrawArr.ePerm.PrefixTemp + "ScriptPcolor_" + eRecVar.VarName1 + "_" + eRecVar.strAll +".ncl";
  DEFINE_PCOLOR_NC_NCL(eFileNC, GrdArr, eRecVar.F, 
		       eDrawArr.DrawContourBathy,
		       eDrawArr.ListLineSegment);
  int IsFE=GrdArr.IsFE;
  //
  std::ofstream OUTncl;
  OUTncl.open(eFileNCL);
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
  OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
  OUTncl << "begin\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; Data reading\n";
  OUTncl << "  ;\n";
  OUTncl << "  f = addfile(\"" << eFileNC << "\", \"r\")\n";
  OUTncl << "  lat  = f->lat\n";
  OUTncl << "  lon  = f->lon\n";
  OUTncl << "  eVarF  = f->field\n";
  OUTncl << "  wks  = gsn_open_wks (\"" << eDrawArr.ePerm.Extension << "\",\"" << FileName << "\")\n";
  if (eDrawArr.DrawContourBathy) {
    OUTncl << "  DEP  = f->dep\n";
    OUTncl << "  res1 = True\n";
    OUTncl << "  res1@gsnTickMarksOn   = False; no tickmarks\n";
    OUTncl << "  res1@gsnDraw          = False; don't draw\n";
    OUTncl << "  res1@gsnFrame         = False; don't advance frame\n";
    OUTncl << "  res1@gsnLeftString    = \"\"; no titles\n";
    OUTncl << "  res1@gsnRightString   = \"\"\n";
    OUTncl << "  res1@tiXAxisString    = \"\"\n";
    OUTncl << "  res1@tiYAxisString    = \"\"\n";
    OUTncl << "  res1@cnLineThicknessF = 1.5; thicker contours\n";
    OUTncl << "  res1@cnLineLabelsOn   = False; no line labels\n";
    OUTncl << "  plot2 = gsn_csm_contour(wks,DEP,res1)\n";
  }
  OUTncl << "  res2 = True               ; plot mods desired\n";
  OUTncl << "  res2@gsnDraw   = False\n";
  OUTncl << "  res2@gsnFrame  = False\n";
  OUTncl << "  res2@gsnMaximize     = True    ; Maximize plot in frame\n";
  //  OUTncl << "  res2@gsnPaperOrientation  = \"Portrait\"\n";
  OUTncl << "  res2@gsnPaperOrientation = \"Landscape\"\n";
  OUTncl << "  ;\n";
  OUTncl << "  ; General frame information\n";
  OUTncl << "  ;\n";
  OUTncl << "  res2@mpProjection = \"Mercator\"\n";
  OUTncl << "  res2@mpLimitMode         = \"Corners\"             ; choose range of map\n";
  OUTncl << "  res2@mpLeftCornerLatF    = min(lat)\n";
  OUTncl << "  res2@mpLeftCornerLonF    = min(lon)\n";
  OUTncl << "  res2@mpRightCornerLatF   = max(lat)\n";
  OUTncl << "  res2@mpRightCornerLonF   = max(lon)\n";
  OUTncl << "  res2@pmTickMarkDisplayMode  = \"Always\"           ; turn on tickmarks\n";
  if (eDrawArr.FillLand) {
    OUTncl << "  res2@mpFillOn      = True\n";
    OUTncl << "  res2@mpDataBaseVersion      = \"" << eDrawArr.GridResolution << "\"          ; use high resolution coast\n";
    OUTncl << "  res2@mpLandFillColor      = \"LightGrey\"\n";
    OUTncl << "  res2@mpLandFillColor       = \"gray\"            ; set land to be gray\n";
  }
  else {
    OUTncl << "  res2@mpFillOn      = False\n";
  }
  OUTncl << "  ;\n";
  OUTncl << "  ; Contour map information\n";
  OUTncl << "  ;\n";
  OUTncl << "  res2@cnFillDrawOrder        = \"PreDraw\"\n";
  OUTncl << "  res2@cnFillOn             = True               ; turn on color for contours\n";
  OUTncl << "  res2@cnLinesOn            = False              ; turn off contour lines\n";
  OUTncl << "  res2@cnLineLabelsOn       = False              ; turn off contour line labels\n";
  OUTncl << "  res2@cnFillMode           = \"" << eDrawArr.cnFillMode << "\"\n";
  OUTncl << "            ; AreaFill : slow and buggy but maybe more beautiful\n";
  OUTncl << "            ; RasterFill : fast and efficient\n";
  OUTncl << "            ; CellFill : similar tp RaterFill\n";
  OUTncl << "  ;  res2@cnRasterSmoothingOn  = True\n";
  OUTncl << "  res2@cnSmoothingOn = " << NCL_bool(eDrawArr.cnSmoothingOn) << "\n";
  OUTncl << "  ;  res2@cnSmoothingDistanceF  = 0.05\n";
  OUTncl << "  res2@cnSmoothingTensionF  = -1\n";
  OUTncl << "  res2@cnLevelSelectionMode = \"ManualLevels\"     ; set manual contour levels\n";
  OUTncl << "  res2@cnMinLevelValF       = " << eRecVar.minval << "  ; set min contour level\n";
  OUTncl << "  res2@cnMaxLevelValF       = " << eRecVar.maxval << "  ; set max contour level\n";
  int nbLevelSpa=eDrawArr.nbLevelSpa;
  double TheLevelSpa=(eRecVar.maxval - eRecVar.minval)/double(nbLevelSpa);
  OUTncl << "  res2@cnLevelSpacingF      = " << TheLevelSpa << "     ; set contour spacing\n";
  int nbLabelStride=eDrawArr.nbLabelStride;
  OUTncl << "  res2@lbLabelStride            = " << nbLabelStride << "\n";
  OUTncl << "  ;  res2@gsnScalarContour     = False               ; contours desired\n";
  if (eDrawArr.DoTitle) {
    OUTncl << "  res2@tiMainString    = \"" << eDrawArr.TitleStr << "\"\n";
    OUTncl << "  res2@tiMainFont      = \"Helvetica\"\n";
    OUTncl << "  res2@tiMainFontHeightF = 0.015\n";
    OUTncl << "  ;  res2@cnTitlePosition  = \"Top\"\n";
  }
  OUTncl << "  res2@gsnSpreadColors      = True               ; use full color map\n";
  OUTncl << "  res2@gsnSpreadColorEnd     = -3\n";
  // LABEL BAR
  OUTncl << "  ;\n";
  OUTncl << "  ; Label bar plotting\n";
  OUTncl << "  ;\n";
  OUTncl << "  res2@lbLabelBarOn = True\n";
  OUTncl << "  res2@lbTitleString    = \"" << eRecVar.Unit << "\"\n";
  OUTncl << "  res2@lbTitleFont      = \"Helvetica\"\n";
  OUTncl << "  res2@lbTitleFontHeightF = 0.015\n";
  OUTncl << "  res2@lbTitleDirection     = \"Across\"\n";
  OUTncl << "  res2@lbTitlePosition = \"Right\"\n";
  OUTncl << "  res2@lbTitleAngleF = 90\n";
  OUTncl << "  res2@lbOrientation        = \"Vertical\"     ; Vertical label bar\n";
  OUTncl << "  res2@pmLabelBarOrthogonalPosF = 0.025          ; move label bar closer\n";
  OUTncl << "  ;  res2@lbHeightF               = 0.7          ; move label bar closer\n";
  OUTncl << "  res2@pmLabelBarDisplayMode = \"Always\"          ; Turn on a label bar.\n";
  OUTncl << "  res2@lbPerimOn             = False             ; no box around it\n";
  OUTncl << "  res2@lbBoxLinesOn         = False               ; Yes/No labelbar box lines\n";
  OUTncl << "  res2@pmLabelBarWidthF = 0.03\n";
  OUTncl << "  ; res2@gsnRightString  = \"Sea surface elevation\"\n";
  OUTncl << "  ; res2@gsnLeftString    = \"Difference\"\n";
  OUTncl << "  ; res2@gsnRightString  = \"\"\n";
  OUTncl << "  ; res2@gsnLeftString    = \"\"\n";
  // COLORMAP
  OUTncl << "  ;\n";
  OUTncl << "  ; Colormap assignation\n";
  OUTncl << "  ;\n";
  OUTncl << "  gsn_define_colormap (wks,\"" << eDrawArr.ColorMap << "\")\n";
  OUTncl << "  ;     other possibilities: hotres, rainbow, ViBlGrWhYeOrRe, BlWhRe, GrayWhiteGray, BlGrYeOrReVi200\n";
  OUTncl << "  i = NhlNewColor(wks,0.8,0.8,0.8)      ; add gray to colormap\n";
  OUTncl << "  i = NhlNewColor(wks,0.9,0.9,0.9)      ; add gray to colormap\n";
  // THE PLOTTING ITSELF
  OUTncl << "  ;\n";
  OUTncl << "  ; Pcolor kind of plot\n";
  OUTncl << "  ;\n";
  if (IsFE == 1) {
    OUTncl << "  res2@sfXArray            = lon\n";
    OUTncl << "  res2@sfYArray            = lat\n";
    if (eDrawArr.UseNativeGrid) {
      OUTncl << "  res2@sfElementNodes      = f->ele\n";
      OUTncl << "  res2@sfFirstNodeIndex    = 0\n";
    }
  }
  else {
    OUTncl << "  eVarF@lat2d=lat\n";
    OUTncl << "  eVarF@lon2d=lon\n";
  }
  OUTncl << "  plot = gsn_csm_contour_map(wks,eVarF,res2)\n";
  int nbLineSeq=eDrawArr.ListLineSegment.size();
  if (nbLineSeq > 0) {
    OUTncl << "  ListLon = f->ListLon\n";
    OUTncl << "  ListLat = f->ListLat\n";
    OUTncl << "  ns = dimsizes(ListLon)\n";
    OUTncl << "  nbLine=ns(0)/2\n";
    OUTncl << "  resLine = True\n";
    OUTncl << "  resLine@gsLineThicknessF = 3.0\n";
    OUTncl << "  resLine@gsLineColor  = \"dodgerblue1\"\n";
    OUTncl << "  resLine@tfPolyDrawOrder = \"PostDraw\"\n";
    OUTncl << "  linesRect = new(nbLine,graphic)\n";
    OUTncl << "  do iLine=0,nbLine-1\n";
    OUTncl << "    xp=ListLon(2*iLine : 2*iLine+1)\n";
    OUTncl << "    yp=ListLat(2*iLine : 2*iLine+1)\n";
    OUTncl << "    linesRect(iLine)=gsn_add_polyline(wks,plot,xp,yp,resLine)\n";
    OUTncl << "    delete(xp)\n";
    OUTncl << "    delete(yp)\n";
    OUTncl << "  end do\n";
  }
  ADD_RIVER(OUTncl, eDrawArr);
  ADD_ANNOTATION_TEXT(OUTncl, eDrawArr.TheAnnot);
  OUTncl << "  draw(plot)\n";
  OUTncl << "  frame(wks)\n";
  OUTncl << "end\n";
  OUTncl.close();
  //
  std::string TargetFile=FileName + "." + eDrawArr.ePerm.Extension;
  CALL_NCL(eDrawArr.ePerm.KeepNC_NCL, 
	   TargetFile, eFileNC, eFileNCL);
}
/* ------------------------------ */


/* ------------------------------ */
struct DrawLinesArr {
  PermanentInfoDrawing ePerm;
  bool DoTitle;
  std::string TitleStr;
  std::string VarName;
  AnnotationRec TheAnnot;
  // Drawing itself
  bool DoTitleLines;
  double TheMax;
  double TheMin;
  bool IsTimeSeries;
  bool PairComparison;
  std::string YAxisString;
  std::vector<std::string> ListName_plot;
  MyVector<double> ListX;
  std::vector<MyVector<double> > ListListVect;
};
/* ------------------------------ */


/* ------------------------------ */
void LINES_DEFINE_NC(std::string const& eFileNC, 
		     DrawLinesArr const& eDrawArr)
{
  netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  int nbArr=eDrawArr.ListListVect.size();
  int nbEntry=eDrawArr.ListListVect[0].size();
  netCDF::NcDim eDimArr=dataFile.addDim("nbArr", nbArr);
  netCDF::NcDim eDimEntry=dataFile.addDim("nbEntry", nbEntry);
  std::vector<std::string> ListDim={"nbArr", "nbEntry"};
  std::vector<std::string> ListDimX={"nbEntry"};
  netCDF::NcVar eVar=dataFile.addVar("ListListVect", "double", ListDim);
  netCDF::NcVar eVarX=dataFile.addVar("ListX", "double", ListDimX);
  //
  double *val;
  val=new double[nbArr*nbEntry];
  int idx=0;
  for (int iArr=0; iArr<nbArr; iArr++) {
    MyVector<double> eVect=eDrawArr.ListListVect[iArr];
    for (int i=0; i<nbEntry; i++) {
      val[idx]=eVect(i);
      idx++;
    }
  }
  eVar.putVar(val);
  delete [] val;
  //
  double *valX;
  valX=new double[nbEntry];
  for (int i=0; i<nbEntry; i++) {
    double eX=eDrawArr.ListX(i);
    valX[i]=eX;
  }
  eVarX.putVar(valX);
  delete [] valX;
}
/* ------------------------------ */


/* ------------------------------ */
void AngleRhoRot(MyMatrix<double> & U_rho, MyMatrix<double> & V_rho, 
		 GridArray const& GrdArr)
{
  int eta_rho=U_rho.rows();
  int xi_rho=U_rho.cols();
  for (int i=0; i<eta_rho; i++)
    for (int j=0; j<xi_rho; j++) {
      double eU=U_rho(i,j);
      double eV=V_rho(i,j);
      double eAng=GrdArr.GrdArrRho.ANG(i,j);
      double NewU=eU*cos(eAng) - eV*sin(eAng);
      double NewV=eU*sin(eAng) + eV*cos(eAng);
      U_rho(i,j)=NewU;
      V_rho(i,j)=NewV;
    }
}
/* ------------------------------ */


/* ------------------------------ */
void CF_EXTRACT_TIME(std::string const& eStrUnitTime, double & ConvertToDay, double & eTimeStart)
{
  std::string YnameYear, YnameMonth, YnameDay;
  std::string YnameHour, YnameMin, YnameSec;
  std::string YnameD, YnameE;
  std::string YnameDate, YnameTime, YnameTimeP;
  std::string eStrTime;
  std::string strSpace=" ";
  int posBlank=STRING_GetCharPositionInString(eStrUnitTime, strSpace);
  std::string Xname=eStrUnitTime.substr(0, posBlank);
  int IsDone=0;
  if (Xname == "days") {
    IsDone=1;
    ConvertToDay=double(1);
  }
  if (Xname == "hours") {
    IsDone=1;
    ConvertToDay=double(1)/double(24);
  }
  if (Xname == "seconds") {
    IsDone=1;
    ConvertToDay=double(1)/double(86400);
  }
  if (IsDone == 0) {
    std::cerr << "We did not find a match for the time unit\n";
    std::cerr << "eStrUnitTime=" << eStrUnitTime << "\n";
    std::cerr << "Xname=" << Xname << "\n";
    std::cerr << "allowed Xname=days/hours/seconds\n";
    exit(1);
  }
  //
  int alen=eStrUnitTime.length();
  std::string Yname=eStrUnitTime.substr(posBlank+1, alen - 1 - posBlank);
  int alenB=Yname.length();
  int posBlankB=STRING_GetCharPositionInString(Yname, strSpace);
  std::string YnameB=Yname.substr(posBlankB+1, alenB - 1 - posBlankB);
  // YnameB should be 1990-01-01 0:0:0 or something like that

  std::string strT="T";
  std::string strZ="Z";
  std::vector<std::string> LStrDateT=STRING_Split(YnameB, strT);
  int sizStrDateT=LStrDateT.size();
  if (sizStrDateT > 1) {
    // case of WW3 that has dates with the file format
    // "days since 1990-01-01T00:00:00Z"
    YnameDate=LStrDateT[0]; // should be 1990-01-01
    std::string eStrB=LStrDateT[1]; // 00:00:00Z
    int alen=eStrUnitTime.length();
    YnameTime=eStrB.substr(0,alen-2);
    std::cerr << "Case of WW3\n";
    std::cerr << "YnameDate=" << YnameDate << "\n";
    std::cerr << "YnameTime=" << YnameTime << "\n";
  }
  else {
    std::vector<std::string> LStrDate=STRING_Split(YnameB, strSpace);
    int sizStrDate=LStrDate.size();
    if (sizStrDate > 1) {
      YnameDate=LStrDate[0]; // should be 1990-01-01
      YnameTime=LStrDate[1]; // should be 0:0:0
      if (sizStrDate > 2) {
	std::string StrShift=LStrDate[2];
	if (StrShift != "GMT") {
	  std::cerr << "Yname=" << Yname << "\n";
	  std::cerr << "YnameB=" << YnameB << "\n";
	  std::cerr << "YnameDate=" << YnameDate << "\n";
	  std::cerr << "YnameTimeP=" << YnameTimeP << "\n";
	  std::cerr << "Need to program that case\n";
	  std::cerr << "Basically time can be of the form 0:0:0 GMT\n";
	  std::cerr << "or other stuff like that\n";
	  exit(1);
	}
      }
    }
    else {
      YnameDate=LStrDate[0];
      YnameTime="00:00:00";
    }
  }
  //
  std::vector<std::string> eVectDate=STRING_Split(YnameDate,"-");
  std::string eStrYear, eStrMonth, eStrDay;
  eStrYear=eVectDate[0];
  eStrMonth=eVectDate[1];
  eStrDay=eVectDate[2];
  int year, month, day;
  std::istringstream(eStrYear) >> year;
  std::istringstream(eStrMonth) >> month;
  std::istringstream(eStrDay) >> day;
  // 
  std::vector<std::string> eVectTime=STRING_Split(YnameTime,":");
  std::string eStrHour, eStrMin, eStrSec;
  eStrHour=eVectTime[0];
  eStrMin=eVectTime[1];
  eStrSec=eVectTime[2];
  int hour, min, sec;
  std::istringstream(eStrHour) >> hour;
  std::istringstream(eStrMin) >> min;
  std::istringstream(eStrSec) >> sec;
  // Now collating
  eTimeStart=DATE_ConvertSix2mjd({year, month, day, hour, min, sec});
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<double> NC_ReadTimeFromFile(std::string const& eFile, std::string const& StringTime)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Error in NC_ReadTimeFromFile\n";
    std::cerr << "Trying to open non-existing file\n";
    std::cerr << "eFile=" << eFile << "\n";
    exit(1);
  }
  std::cerr << "eFile=" << eFile << "\n";
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(StringTime);
  if(data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "StringTime=" << StringTime << "\n";
    exit(1);
  }
  int nbDim=data.getDimCount();
  if (nbDim != 1) {
    std::cerr << "The number of dimensions is not correct\n";
    exit(1);
  }
  netCDF::NcDim eDim=data.getDim(0);
  int siz=eDim.getSize();
  //  std::cerr << "siz=" << siz << "\n";
  double *eVal;
  eVal=new double[siz];
  data.getVar(eVal);
  netCDF::NcVarAtt eTimeAtt=data.getAtt("units");
  char eString[1024]="";
  //  fprintf(stderr, "Before call to getValues\n");
  eTimeAtt.getValues(eString);
  //  fprintf(stderr, "After call to getValues\n");
  //  std::cerr << "eString=" << eString << "\n";
  std::string eStrUnitTime=eString;
  //  std::cerr << "eStrUnitTime=" << eStrUnitTime << "\n";
  //  fprintf(stderr, "After conversion to std::string\n");
  double ConvertToDay, eTimeStart;
  CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart);
  //  fprintf(stderr, "After extraction of key variables from attribute\n");
  std::vector<double> LTime;
  double minTime, maxTime;
  for (int i=0; i<siz; i++) {
    //    std::cerr << "i=" << i << " eVal=" << eVal[i] << "\n";
    //    std::cerr << "ConvertToDay=" << ConvertToDay << "\n";
    double eTimeDay = eVal[i]*ConvertToDay + eTimeStart;
    if (i == 0) {
      minTime=eTimeDay;
      maxTime=eTimeDay;
    }
    else {
      if (eTimeDay > maxTime)
	maxTime=eTimeDay;
      if (eTimeDay < minTime)
	minTime=eTimeDay;
    }
    LTime.push_back(eTimeDay);
  }
  delete [] eVal;
  bool ShowMinMax=false;
  if (ShowMinMax) {
    double minTime, maxTime;
    for (int i=0; i<siz; i++) {
      double eTimeDay = LTime[i];
      if (i == 0) {
	minTime=eTimeDay;
	maxTime=eTimeDay;
      }
      else {
	if (eTimeDay > maxTime)
	  maxTime=eTimeDay;
	if (eTimeDay < minTime)
	  minTime=eTimeDay;
      }
      LTime.push_back(eTimeDay);
    }
    std::cerr << "minTime=" << minTime << " maxTime=" << maxTime << "\n";
    std::string strPresMin=DATE_ConvertMjd2mystringPres(minTime);
    std::string strPresMax=DATE_ConvertMjd2mystringPres(maxTime);
    std::cerr << "strPresMin=" << strPresMin << "\n";
    std::cerr << "strPresMax=" << strPresMax << "\n";
  }
  return LTime;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> Get2DvariableSpecEntry_FD(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Get2DvariableSpecEntry_FD\n";
    std::cerr << "The file eFile=" << eFile << "\n";
    std::cerr << "does not exist\n";
    exit(1);
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if(data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "eVar=" << eVar << "\n";
    std::cerr << "iRec=" << iRec << "\n";
    exit(1);
  }
  int nbDim=data.getDimCount();
  if (nbDim == 3) {
    netCDF::NcDim eDim;
    eDim=data.getDim(0);
    int nbRec=eDim.getSize();
    if (iRec >= nbRec) {
      std::cerr << "Error, iRec is too large\n";
      std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
      exit(1);
    }
    eDim=data.getDim(1);
    size_t eta=eDim.getSize();
    eDim=data.getDim(2);
    size_t xi=eDim.getSize();
    std::vector<size_t> start{size_t(iRec), 0, 0};
    std::vector<size_t> count{1, eta, xi};
    netCDF::NcType eType=data.getType();
    MyMatrix<double> eArr(eta, xi);
    bool IsDone=false;
    if (eType == netCDF::NcType::nc_DOUBLE) {
      IsDone=true;
      double *eVal;
      eVal=new double[eta*xi];
      data.getVar(start, count, eVal);
      int idx=0;
      for (size_t i=0; i<eta; i++)
	for (size_t j=0; j<xi; j++) {
	  eArr(i, j)=eVal[idx];
	  idx++;
	}
      delete [] eVal;
    }
    if (eType == netCDF::NcType::nc_FLOAT) {
      IsDone=true;
      float *eVal;
      eVal=new float[eta*xi];
      data.getVar(start, count, eVal);
      int idx=0;
      for (size_t i=0; i<eta; i++)
	for (size_t j=0; j<xi; j++) {
	  float eValF=eVal[idx];
	  double eValD=double(eValF);
	  eArr(i, j)=eValD;
	  idx++;
	}
      delete [] eVal;
    }
    if (IsDone == false) {
      std::cerr << "no good type founds\n";
      exit(1);
    }
    return eArr;  
  }
  // The file format uses only WET points.
  netCDF::NcDim eDim;
  eDim=data.getDim(0);
  int nbRec=eDim.getSize();
  if (iRec >= nbRec) {
    std::cerr << "Error, iRec is too large\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    exit(1);
  }
  eDim=data.getDim(1);
  int nbWet=eDim.getSize();
  std::vector<double> eVal(nbWet);
  std::vector<size_t> start{size_t(iRec), 0};
  std::vector<size_t> count{1, size_t(nbWet)};
  netCDF::NcType eType=data.getType();
  bool IsDone=false;
  if (eType == netCDF::NcType::nc_DOUBLE) {
    IsDone=true;
    double *eValD;
    eValD=new double[nbWet];
    data.getVar(start, count, eValD);
    for (int iWet=0; iWet<nbWet; iWet++)
      eVal[iWet]=eValD[iWet];
    delete [] eValD;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    IsDone=true;
    float *eValF;
    eValF=new float[nbWet];
    data.getVar(start, count, eValF);
    for (int iWet=0; iWet<nbWet; iWet++)
      eVal[iWet]=double(eValF[iWet]);
    delete [] eValF;
  }
  if (IsDone == false) {
    std::cerr << "no good type founds\n";
    exit(1);
  }
  if (nbWet == GrdArr.GrdArrRho.nbWet) {
    int eta=GrdArr.GrdArrRho.eta;
    int xi=GrdArr.GrdArrRho.xi;
    MyMatrix<double> eArr(eta, xi);
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++)
	eArr(i, j)=0;
    for (int iWet=0; iWet<nbWet; iWet++) {
      int i=GrdArr.GrdArrRho.Idx[iWet];
      int j=GrdArr.GrdArrRho.Jdx[iWet];
      eArr(i, j)=eVal[iWet];
    }
    return eArr;
  }
  if (nbWet == GrdArr.GrdArrU.nbWet) {
    int eta=GrdArr.GrdArrU.eta;
    int xi=GrdArr.GrdArrU.xi;
    MyMatrix<double> eArr(eta, xi);
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++)
	eArr(i, j)=0;
    for (int iWet=0; iWet<nbWet; iWet++) {
      int i=GrdArr.GrdArrU.Idx[iWet];
      int j=GrdArr.GrdArrU.Jdx[iWet];
      eArr(i, j)=eVal[iWet];
    }
    return eArr;
  }
  if (nbWet == GrdArr.GrdArrV.nbWet) {
    int eta=GrdArr.GrdArrV.eta;
    int xi=GrdArr.GrdArrV.xi;
    MyMatrix<double> eArr(eta, xi);
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++)
	eArr(i, j)=0;
    for (int iWet=0; iWet<nbWet; iWet++) {
      int i=GrdArr.GrdArrV.Idx[iWet];
      int j=GrdArr.GrdArrV.Jdx[iWet];
      eArr(i, j)=eVal[iWet];
    }
    return eArr;
  }
  
  std::cerr << "Routine is Get2DvariableSpecEntry\n";
  std::cerr << "eVar=" << eVar << "\n";
  std::cerr << "We did not find the size\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> Get2DvariableSpecEntry_FE(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Get2DvariableSpecEntry_FE\n";
    std::cerr << "The file eFile=" << eFile << "\n";
    std::cerr << "does not exist\n";
    exit(1);
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar data=dataFile.getVar(eVar);
  if(data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "eVar=" << eVar << "\n";
    std::cerr << "iRec=" << iRec << "\n";
    exit(1);
  }
  int nbDim=data.getDimCount();
  if (nbDim != 2) {
    std::cerr << "This command will certainly not work\n";
    std::cerr << "Dimensions are not correct\n";
    exit(1);
  }
  netCDF::NcDim eDim;
  eDim=data.getDim(0);
  int nbRec=eDim.getSize();
  if (iRec >= nbRec) {
    std::cerr << "Error, iRec is too large\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    exit(1);
  }
  eDim=data.getDim(1);
  size_t mnp=eDim.getSize();
  std::vector<size_t> start{size_t(iRec), 0};
  std::vector<size_t> count{1, mnp};
  netCDF::NcType eType=data.getType();
  MyMatrix<double> eArr(mnp, 1);
  bool IsDone=false;
  if (eType == netCDF::NcType::nc_DOUBLE) {
    double *eVal;
    eVal=new double[mnp];
    data.getVar(start, count, eVal);
    for (size_t i=0; i<mnp; i++)
      eArr(i,0)=eVal[i];
    delete [] eVal;
    IsDone=true;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    float *eVal;
    eVal=new float[mnp];
    data.getVar(start, count, eVal);
    for (size_t i=0; i<mnp; i++) {
      float eValF=eVal[i];
      double eValD=double(eValF);
      eArr(i,0)=eValD;
    }
    delete [] eVal;
    IsDone=true;
  }
  if (IsDone == false) {
    std::cerr << "Data reading failed for 2D finite element\n";
    exit(1);
  }
  if (GrdArr.L_IndexSelect) {
    int siz=GrdArr.I_IndexSelect.size();
    MyMatrix<double> eArrRet(siz, 1);
    for (int i=0; i<siz; i++) {
      int iGlob=GrdArr.I_IndexSelect[i];
      eArrRet(i,0)=eArr(iGlob,0);
    }
    return eArrRet;
  }
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> Get2DvariableSpecEntry(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (GrdArr.IsFE == 1) {
    return Get2DvariableSpecEntry_FE(eFile, GrdArr, eVar, iRec);
  }
  return Get2DvariableSpecEntry_FD(eFile, GrdArr, eVar, iRec);
}
/* ------------------------------ */


/* ------------------------------ */
bool NC_TOTALARR_IsVar(TotalArrGetData const& TotalArr, std::string const& eVar)
{
  std::string HisFile;
  int iFile=0;
  ArrayHistory eArr=TotalArr.eArr;
  if (eArr.AppendVarName) {
    HisFile=eArr.ListFileNames[iFile] + eVar + ".nc";;
  }
  else {
    HisFile=eArr.ListFileNames[iFile];
  }
  return NC_IsVar(HisFile, eVar);
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> NETCDF_Get2DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& eVar, double const& eTimeDay)
{
  ArrayHistory eArr=TotalArr.eArr;
  GridArray GrdArr=TotalArr.GrdArr;
  InterpInfo eInterpInfo=GetTimeInterpolationInfo(eArr.ListTime, eTimeDay);
  if (eInterpInfo.UseSingleEntry) {
    int iTime=eInterpInfo.iTimeLow;
    int iFile=eArr.ListIFile[iTime];
    int iRec=eArr.ListIRec[iTime];
    std::string HisFile;
    if (eArr.AppendVarName) {
      HisFile=eArr.ListFileNames[iFile] + eVar + ".nc";
    }
    else {
      HisFile=eArr.ListFileNames[iFile];
    }
    return Get2DvariableSpecEntry(HisFile, GrdArr, eVar, iRec);
  }
  double alphaLow=eInterpInfo.alphaLow;
  int iTimeLow=eInterpInfo.iTimeLow;
  int iFileLow=eArr.ListIFile[iTimeLow];
  int iRecLow=eArr.ListIRec[iTimeLow];
  double alphaUpp=eInterpInfo.alphaUpp;
  int iTimeUpp=eInterpInfo.iTimeUpp;
  int iFileUpp=eArr.ListIFile[iTimeUpp];
  int iRecUpp=eArr.ListIRec[iTimeUpp];
  std::string HisFileLow, HisFileUpp;
  if (eArr.AppendVarName) {
    HisFileLow=eArr.ListFileNames[iFileLow] + eVar + ".nc";
  }
  else {
    HisFileLow=eArr.ListFileNames[iFileLow];
  }
  if (eArr.AppendVarName) {
    HisFileUpp=eArr.ListFileNames[iFileUpp] + eVar + ".nc";;
  }
  else {
    HisFileUpp=eArr.ListFileNames[iFileUpp];
  }
  MyMatrix<double> eVarLow=Get2DvariableSpecEntry(HisFileLow, GrdArr, eVar, iRecLow);
  MyMatrix<double> eVarUpp=Get2DvariableSpecEntry(HisFileUpp, GrdArr, eVar, iRecUpp);
  int eta=eVarLow.rows();
  int xi=eVarLow.cols();
  MyMatrix<double> RetVar(eta, xi);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double eVal=alphaLow*eVarLow(i, j) + alphaUpp*eVarUpp(i, j);
      RetVar(i, j)=eVal;
    }
  return RetVar;
}
/* ------------------------------ */


/* ------------------------------ */
ThreeDimArray<double> Get3DvariableSpecEntry_FD(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (IsExistingFile(eFile) == false) {
    std::cerr << "Get3DvariableSpecEntry_FD\n";
    std::cerr << "The file eFile=" << eFile << "\n";
    std::cerr << "does not exist\n";
    exit(1);
  }
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcDim eDim;
  eDim=dataFile.getDim("s_rho");
  int s_rho=eDim.getSize();
  eDim=dataFile.getDim("s_w");
  int s_w=eDim.getSize();
  netCDF::NcVar data=dataFile.getVar(eVar);
  if(data.isNull()) {
    std::cerr << "Error in accessing to the file\n";
    std::cerr << "eFile=" << eFile << "\n";
    std::cerr << "eVar=" << eVar << "\n";
    exit(1);
  }
  int nbDim=data.getDimCount();
  if (nbDim == 4) {
    eDim=data.getDim(0);
    int nbRec=eDim.getSize();
    if (iRec >= nbRec) {
      std::cerr << "Error, iRec is too large\n";
      std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
      exit(1);
    }
    eDim=data.getDim(1);
    int s_vert=eDim.getSize();
    eDim=data.getDim(2);
    int eta=eDim.getSize();
    eDim=data.getDim(3);
    int xi=eDim.getSize();
    std::vector<size_t> start{size_t(iRec), 0, 0, 0};
    std::vector<size_t> count{1, size_t(s_vert), size_t(eta), size_t(xi)};
    std::vector<double> eVal(s_vert*eta*xi);
    ThreeDimArray<double> eArr(s_vert, eta, xi);
    netCDF::NcType eType=data.getType();
    bool IsDone=false;
    if (eType == netCDF::NcType::nc_DOUBLE) {
      IsDone=true;
      double *eValD;
      eValD=new double[s_vert*eta*xi];
      data.getVar(start, count, eValD);
      for (int i=0; i<s_vert*eta*xi; i++)
	eVal[i]=eValD[i];
      delete [] eValD;
    }
    if (eType == netCDF::NcType::nc_FLOAT) {
      IsDone=true;
      float *eValF;
      eValF=new float[s_vert*eta*xi];
      data.getVar(start, count, eValF);
      for (int i=0; i<s_vert*eta*xi; i++)
	eVal[i]=double(eValF[i]);
      delete [] eValF;
    }
    if (IsDone == false) {
      std::cerr << "no good type founds\n";
      exit(1);
    }
    int idx=0;
    for (int k=0; k<s_vert; k++)
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++) {
	  eArr.assign(k, i, j, eVal[idx]);
	  idx++;
	}
    return eArr;  
  }
  // The file format uses only WET points.
  eDim=data.getDim(0);
  int nbRec=eDim.getSize();
  if (iRec >= nbRec) {
    std::cerr << "Error, iRec is too large\n";
    std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
    exit(1);
  }
  eDim=data.getDim(1);
  int nbWet=eDim.getSize();
  std::vector<double> eVal(nbWet);
  std::vector<size_t> start{size_t(iRec), 0};
  std::vector<size_t> count{1, size_t(nbWet)};
  netCDF::NcType eType=data.getType();
  bool IsDone=false;
  if (eType == netCDF::NcType::nc_DOUBLE) {
    IsDone=true;
    double *eValD;
    eValD=new double[nbWet];
    data.getVar(start, count, eValD);
    for (int i=0; i<nbWet; i++)
      eVal[i]=eValD[i];
    delete [] eValD;
  }
  if (eType == netCDF::NcType::nc_FLOAT) {
    IsDone=true;
    float *eValF;
    eValF=new float[nbWet];
    data.getVar(start, count, eValF);
    for (int iWet=0; iWet<nbWet; iWet++)
      eVal[iWet]=double(eValF[iWet]);
    delete [] eValF;
  }
  if (IsDone == false) {
    std::cerr << "no good type founds\n";
    exit(1);
  }
  if (nbWet == s_rho*GrdArr.GrdArrRho.nbWet) {
    int eta=GrdArr.GrdArrRho.eta;
    int xi=GrdArr.GrdArrRho.xi;
    ThreeDimArray<double> eArr(s_rho, eta, xi);
    for (int k=0; k<s_rho; k++)
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++)
	  eArr.assign(k, i, j, 0);
    int idx=0;
    for (int k=0; k<s_rho; k++) {
      for (int iWet=0; iWet<GrdArr.GrdArrRho.nbWet; iWet++) {
	int i=GrdArr.GrdArrRho.Idx[iWet];
	int j=GrdArr.GrdArrRho.Jdx[iWet];
	eArr.assign(k, i, j, eVal[idx]);
	idx++;
      }
    }
    return eArr;
  }
  if (nbWet == s_rho*GrdArr.GrdArrU.nbWet) {
    int eta=GrdArr.GrdArrU.eta;
    int xi=GrdArr.GrdArrU.xi;
    ThreeDimArray<double> eArr(s_rho, eta, xi);
    for (int k=0; k<s_rho; k++)
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++)
	  eArr.assign(k, i, j, 0);
    int idx=0;
    for (int k=0; k<s_rho; k++) {
      for (int iWet=0; iWet<GrdArr.GrdArrU.nbWet; iWet++) {
	int i=GrdArr.GrdArrU.Idx[iWet];
	int j=GrdArr.GrdArrU.Jdx[iWet];
	eArr.assign(k, i, j, eVal[idx]);
	idx++;
      }
    }
    return eArr;
  }
  if (nbWet == s_rho*GrdArr.GrdArrV.nbWet) {
    int eta=GrdArr.GrdArrV.eta;
    int xi=GrdArr.GrdArrV.xi;
    ThreeDimArray<double> eArr(s_rho, eta, xi);
    for (int k=0; k<s_rho; k++)
      for (int i=0; i<eta; i++)
	for (int j=0; j<xi; j++)
	  eArr.assign(k, i, j, 0);
    int idx=0;
    for (int k=0; k<s_rho; k++) {
      for (int iWet=0; iWet<GrdArr.GrdArrV.nbWet; iWet++) {
	int i=GrdArr.GrdArrV.Idx[iWet];
	int j=GrdArr.GrdArrV.Jdx[iWet];
	eArr.assign(k, i, j, eVal[idx]);
	idx++;
      }
    }
    return eArr;
  }
  std::cerr << "Routine is Get2DvariableSpecEntry\n";
  std::cerr << "eVar=" << eVar << "\n";
  std::cerr << "s_rho=" << s_rho << " s_w=" << s_w << "\n";
  std::cerr << "nbWet=" << nbWet << "\n";
  std::cerr << "nbWetRho=" << GrdArr.GrdArrRho.nbWet << "\n";
  std::cerr << "  nbWetU=" << GrdArr.GrdArrU.nbWet << "\n";
  std::cerr << "  nbWetV=" << GrdArr.GrdArrV.nbWet << "\n";
  std::cerr << "We did not find the size\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
ThreeDimArray<double> Get3DvariableSpecEntry(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
  if (GrdArr.IsFE == 1) {
    std::cerr << "You need to program this part of the program\n";
    exit(1);
  }
  return Get3DvariableSpecEntry_FD(eFile, GrdArr, eVar, iRec);
}
/* ------------------------------ */


/* ------------------------------ */
ThreeDimArray<double> Get3DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& eVar, double const& eTimeDay)
{
  ArrayHistory eArr=TotalArr.eArr;
  GridArray GrdArr=TotalArr.GrdArr;
  InterpInfo eInterpInfo=GetTimeInterpolationInfo(eArr.ListTime, eTimeDay);
  //  std::cerr << "eTimeDay=" << eTimeDay << "\n";
  //  std::cerr << "UseSingleEntry=" << eInterpInfo.UseSingleEntry << "\n";
  if (eInterpInfo.UseSingleEntry) {
    double iTime=eInterpInfo.iTimeLow;
    int iFile=eArr.ListIFile[iTime];
    int iRec=eArr.ListIRec[iTime];
    std::string HisFile=eArr.ListFileNames[iFile];
    //    std::cerr << "HisFile=" << HisFile << "\n";
    //    std::cerr << "iTime=" << iTime << "\n";
    //    std::cerr << "iFile=" << iFile << " iRec=" << iRec << "\n";
    return Get3DvariableSpecEntry(HisFile, GrdArr, eVar, iRec);
  }
  double alphaLow=eInterpInfo.alphaLow;
  int iTimeLow=eInterpInfo.iTimeLow;
  int iFileLow=eArr.ListIFile[iTimeLow];
  int iRecLow=eArr.ListIRec[iTimeLow];
  double alphaUpp=eInterpInfo.alphaUpp;
  int iTimeUpp=eInterpInfo.iTimeUpp;
  int iFileUpp=eArr.ListIFile[iTimeUpp];
  int iRecUpp=eArr.ListIRec[iTimeUpp];
  std::string HisFileLow=eArr.ListFileNames[iFileLow];
  std::string HisFileUpp=eArr.ListFileNames[iFileUpp];
  ThreeDimArray<double> eVarLow=Get3DvariableSpecEntry(HisFileLow, GrdArr, eVar, iRecLow);
  ThreeDimArray<double> eVarUpp=Get3DvariableSpecEntry(HisFileUpp, GrdArr, eVar, iRecUpp);
  std::vector<size_t> LDim=eVarLow.GetDims();
  int s_vert=LDim[0];
  int eta=LDim[1];
  int xi=LDim[2];
  ThreeDimArray<double> RetVar(s_vert, eta, xi);
  for (int k=0; k<s_vert; k++)
    for (int i=0; i<eta; i++)
      for (int j=0; j<xi; j++) {
	double eVal=alphaLow*eVarLow.get(k, i, j) + alphaUpp*eVarUpp.get(k, i, j);
	RetVar.assign(k, i, j, eVal);
      }
  return RetVar;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray GRIB_ReadGridArray(std::string const& FileName, std::string const& eModelName)
{
  //  grib_context *c;
  grib_handle *h = NULL;
  int err;
  FILE* in = NULL;
  in = fopen(FileName.c_str(),"r");
  //  unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;

  int nbMessage=0;
  while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) { 
    nbMessage++;
    std::cerr << "nbMessage=" << nbMessage << "\n";
    //    std::cerr << "err=" << err << "\n";
    if (err != GRIB_SUCCESS)
      GRIB_CHECK(err,0);
    //
    /*
    char* name_space="ls";
    grib_keys_iterator* kiter=NULL;
    kiter=grib_keys_iterator_new(h,key_iterator_filter_flags,name_space);
    if (!kiter) {
      printf("ERROR: Unable to create keys iterator\n");
      exit(1);
    }
    std::string ShortNameValue;
    while(grib_keys_iterator_next(kiter)) {
      int MAX_VAL_LEN=1024;
      char value[MAX_VAL_LEN];
      size_t vlen=MAX_VAL_LEN;
      const char* name = grib_keys_iterator_get_name(kiter);
      vlen=MAX_VAL_LEN;
      bzero(value,vlen);
      GRIB_CHECK(grib_get_string(h,name,value,&vlen),name);
      std::string nameStr=name;
      std::string valueStr=value;
      if (nameStr == "shortName") {
	ShortNameValue=valueStr;
      }
      //printf("%s = %s\n",name,value);
    }
    grib_keys_iterator_delete(kiter); */
    //
    long Ni, Nj, numberOfDataPoints;
    GRIB_CHECK(grib_get_long(h,"Ni",&Ni),0);
    //    std::cerr << "Ni=" << Ni << "\n";
    GRIB_CHECK(grib_get_long(h,"Nj",&Nj),0);
    //    std::cerr << "Nj=" << Nj << "\n";
    GRIB_CHECK(grib_get_long(h,"numberOfDataPoints",&numberOfDataPoints),0);
    CosmoGridInfo eCosmoGrid;
    if (eModelName == "GRIB_COSMO") {
      double latitudeOfSouthernPoleInDegrees, longitudeOfSouthernPoleInDegrees, angleOfRotationInDegrees;
      double latitudeOfFirstGridPointInDegrees, longitudeOfFirstGridPointInDegrees, latitudeOfLastGridPointInDegrees, longitudeOfLastGridPointInDegrees;
      double iDirectionIncrementInDegrees, jDirectionIncrementInDegrees;
      //
      // The southern pole coordinates
      //
      GRIB_CHECK(grib_get_double(h,"latitudeOfSouthernPoleInDegrees",&latitudeOfSouthernPoleInDegrees),0);
      std::cerr << "latitudeOfSouthernPoleInDegrees=" << latitudeOfSouthernPoleInDegrees << "\n";
      eCosmoGrid.latitudeOfSouthernPoleInDegrees=latitudeOfSouthernPoleInDegrees;
      //
      GRIB_CHECK(grib_get_double(h,"longitudeOfSouthernPoleInDegrees",&longitudeOfSouthernPoleInDegrees),0);
      std::cerr << "longitudeOfSouthernPoleInDegrees=" << longitudeOfSouthernPoleInDegrees << "\n";
      eCosmoGrid.longitudeOfSouthernPoleInDegrees=longitudeOfSouthernPoleInDegrees;
      //
      GRIB_CHECK(grib_get_double(h,"angleOfRotationInDegrees",&angleOfRotationInDegrees),0);
      std::cerr << "angleOfRotationInDegrees=" << angleOfRotationInDegrees << "\n";
      eCosmoGrid.angleOfRotationInDegrees=angleOfRotationInDegrees;
      //
      // The first and last longitudes
      //
      GRIB_CHECK(grib_get_double(h,"latitudeOfFirstGridPointInDegrees",&latitudeOfFirstGridPointInDegrees),0);
      std::cerr << "latitudeOfFirstGridPointInDegrees=" << latitudeOfFirstGridPointInDegrees << "\n";
      eCosmoGrid.latitudeOfFirstGridPointInDegrees=latitudeOfFirstGridPointInDegrees;
      //
      GRIB_CHECK(grib_get_double(h,"longitudeOfFirstGridPointInDegrees",&longitudeOfFirstGridPointInDegrees),0);
      std::cerr << "longitudeOfFirstGridPointInDegrees=" << longitudeOfFirstGridPointInDegrees << "\n";
      eCosmoGrid.longitudeOfFirstGridPointInDegrees=longitudeOfFirstGridPointInDegrees;
      //
      GRIB_CHECK(grib_get_double(h,"latitudeOfLastGridPointInDegrees",&latitudeOfLastGridPointInDegrees),0);
      std::cerr << "latitudeOfLastGridPointInDegrees=" << latitudeOfLastGridPointInDegrees << "\n";
      eCosmoGrid.latitudeOfLastGridPointInDegrees=latitudeOfLastGridPointInDegrees;
      //
      GRIB_CHECK(grib_get_double(h,"longitudeOfLastGridPointInDegrees",&longitudeOfLastGridPointInDegrees),0);
      std::cerr << "longitudeOfLastGridPointInDegrees=" << longitudeOfLastGridPointInDegrees << "\n";
      eCosmoGrid.longitudeOfLastGridPointInDegrees=longitudeOfLastGridPointInDegrees;
      //
      // the increments
      //
      GRIB_CHECK(grib_get_double(h,"iDirectionIncrementInDegrees",&iDirectionIncrementInDegrees),0);
      std::cerr << "iDirectionIncrementInDegrees=" << iDirectionIncrementInDegrees << "\n";
      eCosmoGrid.iDirectionIncrementInDegrees=iDirectionIncrementInDegrees;
      //
      GRIB_CHECK(grib_get_double(h,"jDirectionIncrementInDegrees",&jDirectionIncrementInDegrees),0);
      std::cerr << "jDirectionIncrementInDegrees=" << jDirectionIncrementInDegrees << "\n";
      eCosmoGrid.jDirectionIncrementInDegrees=jDirectionIncrementInDegrees;
      //
    }
    //    std::cerr << "NumberOfDataPoints=" << numberOfDataPoints << "\n";
    double *lats, *lons, *values;
    size_t size=numberOfDataPoints;
    size_t* sizePtr=&size;
    lats=(double*)malloc(size*sizeof(double));
    lons=(double*)malloc(size*sizeof(double));
    values=(double*)malloc(size*sizeof(double));
    err=grib_get_data(h, lats, lons, values, sizePtr);
    grib_handle_delete(h);
    if (err != GRIB_SUCCESS)
      GRIB_CHECK(err,0);
    int eta_rho=Ni;
    int xi_rho=Nj;
    MyMatrix<double> LON(eta_rho, xi_rho);
    MyMatrix<double> LAT(eta_rho, xi_rho);
    int idx=0;
    for (int j=0; j<xi_rho; j++)
      for (int i=0; i<eta_rho; i++) {
	LON(i,j)=lons[idx];
	LAT(i,j)=lats[idx];
	idx++;
      }
    if (eModelName == "GRIB_COSMO") {
      Apply_COSMO_Transformation(LON, LAT, eCosmoGrid);
    }
    double MinLON=LON.minCoeff();
    double MaxLON=LON.maxCoeff();
    double MinLAT=LAT.minCoeff();
    double MaxLAT=LAT.maxCoeff();
    std::cerr << "LON min=" << MinLON << " max=" << MaxLON << "\n";
    std::cerr << "LAT min=" << MinLAT << " max=" << MaxLAT << "\n";
    std::cerr << "[0,0]            lon=" << LON(0,0) << " lat=" << LAT(0,0) << "\n";
    std::cerr << "[eta_rho,0]      lon=" << LON(eta_rho-1,0) << " lat=" << LAT(eta_rho-1,0) << "\n";
    std::cerr << "[eta_rho,xi_rho] lon=" << LON(eta_rho-1,xi_rho-1) << " lat=" << LAT(eta_rho-1,xi_rho-1) << "\n";
    std::cerr << "[0,xi_rho]       lon=" << LON(0,xi_rho-1) << " lat=" << LAT(0,xi_rho-1) << "\n";
    // need to assign ANG
    MyMatrix<double> ANG=CreateAngleMatrix(LON, LAT);
    GridArray eGridArrGRIB;
    eGridArrGRIB.IsFE=0;
    eGridArrGRIB.GrdArrRho.LON=LON;
    eGridArrGRIB.GrdArrRho.LAT=LAT;
    eGridArrGRIB.GrdArrRho.ANG=ANG;
    eGridArrGRIB.GrdArrRho.HaveDEP=false;
    free(lats);
    free(lons);
    free(values);
    return eGridArrGRIB;
  }
  std::cerr << "Failed to find the variable\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
std::vector<GRIB_MessageInfo> GRIB_GetAllListPairTime(std::string const& FileName)
{
  grib_handle *h = NULL;
  int err;
  FILE* in = NULL;
  in = fopen(FileName.c_str(),"r");
  unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;
  std::vector<GRIB_MessageInfo>  ListInfo;
  int idx=0;
  while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) { 
    //    std::cerr << "idx=" << idx << "\n";
    if (err != GRIB_SUCCESS)
      GRIB_CHECK(err,0);
    //
    //    std::string eStr="ls";
    char name_space[3]="ls";
    grib_keys_iterator* kiter=NULL;
    kiter=grib_keys_iterator_new(h,key_iterator_filter_flags,name_space);
    if (!kiter) {
      printf("ERROR: Unable to create keys iterator\n");
      exit(1);
    }
    std::string ShortNameValue;
    std::string NameValue;
    std::string DataDateValue = "unset";
    std::string DataTimeValue = "unset";
    std::string StepRangeValue;
    /*
    std::string latitudeOfSouthernPoleInDegreesValue = "unset";
    std::string longitudeOfSouthernPoleInDegreesValue = "unset";
    std::string angleOfRotationInDegreesValue = "unset";*/
    while(grib_keys_iterator_next(kiter)) {
      int MAX_VAL_LEN=1024;
      char value[MAX_VAL_LEN];
      size_t vlen=MAX_VAL_LEN;
      const char* name = grib_keys_iterator_get_name(kiter);
      vlen=MAX_VAL_LEN;
      bzero(value,vlen);
      GRIB_CHECK(grib_get_string(h,name,value,&vlen),name);
      std::string nameStr=name;
      std::string valueStr=value;
      if (nameStr == "shortName") {
	ShortNameValue=valueStr;
      }
      if (nameStr == "name") {
	NameValue=valueStr;
      }
      if (nameStr == "dataDate") {
	DataDateValue=valueStr;
      }
      if (nameStr == "dataTime") {
	DataTimeValue=valueStr;
      }
      if (nameStr == "stepRange") {
	StepRangeValue=valueStr;
      }
      if (nameStr == "stepRange") {
	StepRangeValue=valueStr;
      }
      /*      std::cerr << "nameStr=" << nameStr << "\n";
      if (nameStr == "latitudeOfSouthernPoleInDegrees") {
	latitudeOfSouthernPoleInDegreesValue=valueStr;
      }
      if (nameStr == "longitudeOfSouthernPoleInDegrees") {
	longitudeOfSouthernPoleInDegreesValue=valueStr;
      }
      if (nameStr == "angleOfRotationInDegrees") {
	angleOfRotationInDegreesValue=valueStr;
	}*/
      //      std::cerr << "name=" << nameStr << " value=" << valueStr << "\n";
      //printf("%s = %s\n",name,value);
    }
    grib_keys_iterator_delete(kiter);
    //    std::cerr << "ShortNameValue=" << ShortNameValue << "\n";
    //
    //    std::cerr << "FileName=" << FileName << "\n";
    //    std::cerr << "StepRangeValue=" << StepRangeValue << "\n";
    /*    std::cerr << "latitudeOfSouthernPoleInDegrees=" << latitudeOfSouthernPoleInDegreesValue << "\n";
    std::cerr << "longitudeOfSouthernPoleInDegrees=" << longitudeOfSouthernPoleInDegreesValue << "\n";
    double latitudeOfSouthernPoleInDegrees;
    double longitudeOfSouthernPoleInDegrees;
    double angleOfRotationInDegrees;
    if (latitudeOfSouthernPoleInDegreesValue == "unset") {
      latitudeOfSouthernPoleInDegrees = double(-400);
    }
    else {
      latitudeOfSouthernPoleInDegrees = stod(latitudeOfSouthernPoleInDegreesValue);
    }
    if (longitudeOfSouthernPoleInDegreesValue == "unset") {
      longitudeOfSouthernPoleInDegrees = double(-400);
    }
    else {
      longitudeOfSouthernPoleInDegrees = stod(longitudeOfSouthernPoleInDegreesValue);
    }
    if (angleOfRotationInDegreesValue == "unset") {
      angleOfRotationInDegrees = double(-400);
    }
    else {
      angleOfRotationInDegrees = stod(angleOfRotationInDegreesValue);
      }*/
    int stepRange;
    if (StepRangeValue == "0") {
      stepRange=0;
    }
    else {
      std::vector<std::string> LStr=STRING_Split(StepRangeValue, "-");
      int siz=LStr.size();
      if (siz == 1) {
	stepRange=stoi(StepRangeValue);
      }
      else {
	if (siz != 2) {
	  std::cerr << "Inconsistency in our assumptions\n";
	  exit(1);
	}
	//	std::cerr << "LStr[1]=" << LStr[1] << "\n";
	stepRange=stoi(LStr[1]);
      }
    }
    //    std::cerr << "StepRangeValue=" << StepRangeValue << " stepRange=" << stepRange << "\n";
    //
    //    std::cerr << "DataDateValue=" << DataDateValue << "\n";
    int year, month, day;
    if (DataDateValue == "unset") {
      long dataDate;
      GRIB_CHECK(grib_get_long(h,"dataDate",&dataDate),0);
      //      std::cerr << "dataDate=" << dataDate << "\n";
      day=dataDate % 100;
      int res1=(dataDate - day)/100;
      month=res1 % 100;
      int res2=(res1 - month)/100;
      year=res2;
      //      std::cerr << "year=" << year << " month=" << month << " day=" << day << "\n";
      //      std::cerr << "We should think of something\n";
      //      exit(1);
    }
    else {
      std::string yearStr=DataDateValue.substr(0,4);
      //    std::cerr << "yearStr=" << yearStr << "\n";
      std::string monthStr=DataDateValue.substr(4,2);
      //    std::cerr << "monthStr=" << monthStr << "\n";
      std::string dayStr=DataDateValue.substr(6,2);
      //    std::cerr << "dayStr=" << dayStr << "\n";
      year=stoi(yearStr);
      month=stoi(monthStr);
      day=stoi(dayStr);
    }
    //
    if (DataTimeValue == "unset") {
      long dataTime;
      GRIB_CHECK(grib_get_long(h,"dataTime",&dataTime),0);
      //      std::cerr << "dataTime=" << dataTime << "\n";
      if (dataTime == 0) {
	DataTimeValue="0000";
      }
      else {
	DataTimeValue=LongToString(dataTime);
      }
    }
    grib_handle_delete(h);
    //    std::cerr << "DataTimeValue=" << DataTimeValue << "\n";
    std::string HourStr=DataTimeValue.substr(0,2);
    //    std::cerr << "HourStr=" << HourStr << "\n";
    std::string MinStr=DataTimeValue.substr(2,2);
    //    std::cerr << "MinStr=" << MinStr << "\n";
    int hour=stoi(HourStr);
    int min=stoi(MinStr);
    //
    double eTimeStart=DATE_ConvertSix2mjd({year, month, day, hour, min, 0});
    double eTime=eTimeStart + double(stepRange)/double(24);
    //
    GRIB_MessageInfo eInfo{ShortNameValue, NameValue, idx, eTime, eTimeStart, stepRange};
    ListInfo.push_back(eInfo);
    idx++;
  }
  fclose(in);
  return ListInfo;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> GRIB_Get2DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& eVar, double const& eTimeDay)
{
  ArrayHistory eArr=TotalArr.eArr;
  GridArray GrdArr=TotalArr.GrdArr;
  InterpInfo eInterpInfo=GetTimeInterpolationInfo(eArr.ListTime, eTimeDay);
  if (eInterpInfo.UseSingleEntry) {
    int iTime=eInterpInfo.iTimeLow;
    return GRIB_Read2DVariable(eArr.ListListFileNames[iTime], eVar);
  }
  double alphaLow=eInterpInfo.alphaLow;
  int iTimeLow=eInterpInfo.iTimeLow;
  double alphaUpp=eInterpInfo.alphaUpp;
  int iTimeUpp=eInterpInfo.iTimeUpp;
  MyMatrix<double> eVarLow=GRIB_Read2DVariable(eArr.ListListFileNames[iTimeLow], eVar);
  MyMatrix<double> eVarUpp=GRIB_Read2DVariable(eArr.ListListFileNames[iTimeUpp], eVar);
  int eta=eVarLow.rows();
  int xi=eVarLow.cols();
  MyMatrix<double> RetVar(eta, xi);
  for (int i=0; i<eta; i++)
    for (int j=0; j<xi; j++) {
      double eVal=alphaLow*eVarLow(i, j) + alphaUpp*eVarUpp(i, j);
      RetVar(i, j)=eVal;
    }
  return RetVar;
}
/* ------------------------------ */


/* ------------------------------ */
DrawArr CommonAssignation_DrawArr(FullNamelist const& eFull, GridArray const& GrdArr, PermanentInfoDrawing const& ePerm)
{
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  DrawArr eDrawArr;
  eDrawArr.ePerm=ePerm;
  int nbLevelSpa=eBlPLOT.ListIntValues.at("nbLevelSpa");
  std::string strNature=eBlPROC.ListStringValues["__NaturePlot"];
  /*
  if (strNature == "PAIR") {
    int res=nbLevelSpa % 2;
    if (res == 0)
      nbLevelSpa++;
      }*/
  eDrawArr.nbLevelSpa=nbLevelSpa;
  int nbLabelStride=eBlPLOT.ListIntValues.at("nbLabelStride");
  eDrawArr.nbLabelStride=nbLabelStride;
  eDrawArr.DrawRiver=eBlPLOT.ListBoolValues.at("DrawRiver");
  eDrawArr.TheAnnot.DrawAnnotation=eBlPLOT.ListBoolValues.at("DrawAnnotation");
  eDrawArr.TheAnnot.AnnotationLon=eBlPLOT.ListDoubleValues.at("AnnotationLon");
  eDrawArr.TheAnnot.AnnotationLat=eBlPLOT.ListDoubleValues.at("AnnotationLat");
  eDrawArr.TheAnnot.AnnotationText=eBlPLOT.ListStringValues.at("AnnotationText");
  eDrawArr.FillLand=eBlPLOT.ListBoolValues.at("FillLand");
  eDrawArr.GridResolution=eBlPLOT.ListStringValues.at("GridResolution");
  eDrawArr.UseNativeGrid=eBlPLOT.ListBoolValues.at("UseNativeGrid");
  eDrawArr.DoTitle=eBlPLOT.ListBoolValues.at("DoTitle");
  eDrawArr.cnFillMode=eBlPLOT.ListStringValues.at("cnFillMode");
  eDrawArr.cnSmoothingOn=eBlPLOT.ListBoolValues.at("cnSmoothingOn");
  eDrawArr.ColorMap=eBlPLOT.ListStringValues.at("ColorMap");
  if (eBlPLOT.ListBoolValues.at("DrawContourBathy")) {
    eDrawArr.DrawContourBathy=true;
    if (GrdArr.GrdArrRho.HaveDEP == false) {
      std::cerr << "You asked for DrawContourBathy\n";
      std::cerr << "However, the bathymetry is not available\n";
      exit(1);
    }
  }
  else {
    eDrawArr.DrawContourBathy=false;
  }
  return eDrawArr;
}
/* ------------------------------ */


/* ------------------------------ */
void PLOT_DIFF_FD_RHO_PCOLOR(GridArray const& GrdArr, 
			     RecVar const& eRecVar1, 
			     RecVar const& eRecVar2, 
			     PermanentInfoDrawing const& ePerm,
			     FullNamelist const& eFull)
{
  DrawArr eDrawArr=CommonAssignation_DrawArr(eFull, GrdArr, ePerm);
  std::string FileName=ePerm.eDir + eRecVar1.VarName1 + "_" + eRecVar1.strAll;
  std::string Name1=eFull.ListBlock.at("PROC").ListStringValues.at("Name1");
  std::string Name2=eFull.ListBlock.at("PROC").ListStringValues.at("Name2");
  MyMatrix<double> eDiff = eRecVar1.F - eRecVar2.F;
  bool PrintMMA=eFull.ListBlock.at("PLOT").ListBoolValues.at("PrintMMA");
  if (PrintMMA) {
    PrintMMA_FCT(eDiff, GrdArr.GrdArrRho.MSK, eRecVar1.VarName1);
  }
  std::string TitleStr="diff. " + eRecVar1.VarName2 + " (" + Name1 + "-" + Name2 + ")";
  TitleStr=TitleStr + " at " + eRecVar1.strPres;
  eDrawArr.TitleStr=TitleStr;
  SingleBlock eBlPLOT=eFull.ListBlock.at("PLOT");
  std::string strColorMapDiff=eBlPLOT.ListStringValues.at("ColorMapDiff");
  eDrawArr.ColorMap=strColorMapDiff;
  RecVar eRecVar=eRecVar1;
  eRecVar.F=eDiff;
  eRecVar.minval=eRecVar1.mindiff;
  eRecVar.maxval=eRecVar1.maxdiff;
  PLOT_FD_PCOLOR(FileName, 
		 GrdArr, 
		 eDrawArr,
		 eRecVar);
}
/* ------------------------------ */


/* ------------------------------ */
void SINGLE_PLOT_QUIVER(GridArray const& GrdArr, 
			RecVar const& eRecVar, 
			PermanentInfoDrawing const& ePerm, 
			FullNamelist const& eFull)
{
  DrawArr eDrawArr=CommonAssignation_DrawArr(eFull, GrdArr, ePerm);
  std::string FileName=ePerm.eDir + eRecVar.VarName1 + "_" + eRecVar.strAll;
  std::string TitleStr=eRecVar.VarName2 + " at " + eRecVar.strPres;
  eDrawArr.TitleStr=TitleStr;
  PLOT_QUIVER(FileName,
	      GrdArr,
	      eDrawArr,
	      eRecVar);
}
/* ------------------------------ */


/* ------------------------------ */
void SINGLE_PLOT_PCOLOR(GridArray const& GrdArr, 
			RecVar const& eRecVar, 
			PermanentInfoDrawing const& ePerm, 
			FullNamelist const& eFull)
{
  DrawArr eDrawArr=CommonAssignation_DrawArr(eFull, GrdArr, ePerm);
  std::string FileName=ePerm.eDir + eRecVar.VarName1 + "_" + eRecVar.strAll;
  std::string TitleStr=eRecVar.VarName2 + " at " + eRecVar.strPres;
  eDrawArr.TitleStr=TitleStr;
  bool PrintMMA=eFull.ListBlock.at("PLOT").ListBoolValues.at("PrintMMA");
  if (PrintMMA) {
    PrintMMA_FCT(eRecVar.F, GrdArr.GrdArrRho.MSK, eRecVar.VarName1);
  }
  PLOT_FD_PCOLOR(FileName, 
		 GrdArr, 
		 eDrawArr,
		 eRecVar);
}
/* ------------------------------ */


/* ------------------------------ */
void GENERAL_PLOT_SINGLE(GridArray const& GrdArr,
			 RecVar const& eRecVar, 
			 PermanentInfoDrawing const& ePerm, 
			 FullNamelist const& eFull)
{
  if (eRecVar.VarNature == "uv") {
    SINGLE_PLOT_QUIVER(GrdArr, eRecVar, ePerm, eFull);
    return;
  }
  if (eRecVar.VarNature == "rho") {
    SINGLE_PLOT_PCOLOR(GrdArr, eRecVar, ePerm, eFull);
    return;
  }
  std::cerr << "No matching VarNature\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
void GRID_PLOTTING(TotalArrGetData const& TotalArr, std::string const& GridFile, PermanentInfoDrawing const& ePerm, FullNamelist const& eFull)
{
  DrawArr eDrawArr=CommonAssignation_DrawArr(eFull, TotalArr.GrdArr, ePerm);
  if (TotalArr.eArr.KindArchive == "NETCDF") {
    if (NC_IsVar(GridFile, "MSK_att_ocn")) {
      MyMatrix<double> MSK_att_ocn=NC_Read2Dvariable(GridFile, "MSK_att_ocn");
      RecVar eRecVar;
      eRecVar.strAll="unset";
      eRecVar.VarName1="MSK_att_ocn";
      eRecVar.VarName2="Mask attainment oceanic model";
      eRecVar.minval=0;
      eRecVar.maxval=1;
      eRecVar.Unit="nondim.";
      eRecVar.F=MSK_att_ocn;
      std::string FileName=ePerm.eDir + "MSK_att_ocn";
      eDrawArr.TitleStr=eRecVar.VarName2;
      PLOT_FD_PCOLOR(FileName, 
		     TotalArr.GrdArr, 
		     eDrawArr,
		     eRecVar);
    }
    if (NC_IsVar(GridFile, "MSK_att_wav")) {
      MyMatrix<double> MSK_att_wav=NC_Read2Dvariable(GridFile, "MSK_att_wav");
      RecVar eRecVar;
      eRecVar.strAll="unset";
      eRecVar.VarName1="MSK_att_wav";
      eRecVar.VarName2="Mask attainment wav model";
      eRecVar.minval=0;
      eRecVar.maxval=1;
      eRecVar.Unit="nondim.";
      eRecVar.F=MSK_att_wav;
      std::string FileName=ePerm.eDir + "MSK_att_wav";
      eDrawArr.TitleStr=eRecVar.VarName2;
      PLOT_FD_PCOLOR(FileName, 
		     TotalArr.GrdArr, 
		     eDrawArr,
		     eRecVar);
    }
  }
}
/* ------------------------------ */


/* ------------------------------ */
void GENERAL_PLOT_PAIR(GridArray const& GrdArr, DrawArr const& eDrawArr, RecVar const& eRecVar1, RecVar const& eRecVar2, PermanentInfoDrawing const& ePerm, FullNamelist const& eFull)
{
  if (eRecVar1.VarNature == "uv") {
    std::cerr << "For the plotting of pair and the field of the type UV\n";
    std::cerr << "we do not yet have relevant code\n";
    std::cerr << "VarName1=" << eRecVar1.VarName1 << "\n";
    return;
  }
  if (eRecVar1.VarNature == "rho") {
    PLOT_DIFF_FD_RHO_PCOLOR(GrdArr, eRecVar1, eRecVar2, ePerm, eFull);
    return;
  }
  std::cerr << "No matching VarNature\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
ArrayHistory NC_ReadArrayHistory_Kernel(std::string const& HisPrefix, std::string const& StringTime)
{
  double FirstTime, LastTime;
  std::vector<std::string> ListFileNames;
  std::vector<int> ListIFile;
  std::vector<int> ListIRec;
  std::vector<double> ListTime;
  if (IsExistingFile(HisPrefix)) {
    ArrayHistory eArr;
    std::vector<double> LTime=NC_ReadTimeFromFile(HisPrefix, StringTime);
    ListFileNames.push_back(HisPrefix);
    int siz=LTime.size();
    for (int i=0; i<siz; i++) {
      ListIFile.push_back(0);
      ListIRec.push_back(i);
      ListTime.push_back(LTime[i]);
    }
    FirstTime=ListTime[0];
    LastTime=ListTime[siz-1];
  }
  else {
    int iFileBegin=0;
    while(1) {
      iFileBegin++;
      std::string TheHisFile=HisPrefix + StringNumber(iFileBegin, 4) + ".nc";
      if (IsExistingFile(TheHisFile) == 1) 
	break;
      if (iFileBegin == 9999) {
	std::cerr << "maybe you specified wrong HisPrefix\n";
	std::cerr << "HisPrefix = " << HisPrefix << "\n";
	std::cerr << " there is no files  HisPrefix????.nc\n";
	std::cerr << "Please correct\n";
	exit(1);
      }
    }
    //    std::cerr << "iFileBegin=" << iFileBegin << "\n";
    int iFileEnd=iFileBegin;
    while(1) {
      std::string TheHisFile=HisPrefix + StringNumber(iFileEnd+1, 4) + ".nc";
      if (IsExistingFile(TheHisFile) == 0)
	break;
      iFileEnd++;
    }
    //    std::cerr << "iFileEnd=" << iFileEnd << "\n";
    std::string TheHisFileBegin=HisPrefix + StringNumber(iFileBegin, 4) + ".nc";
    std::vector<double> LTimeBegin=NC_ReadTimeFromFile(TheHisFileBegin, StringTime);
    int nbRecBegin=LTimeBegin.size();
    double DeltaTime;
    if (nbRecBegin > 1) {
      DeltaTime=LTimeBegin[1]-LTimeBegin[0];
    }
    else {
      if (iFileEnd > iFileBegin) {
	std::string TheHisFile=HisPrefix + StringNumber(iFileBegin+1, 4) + ".nc";
	std::vector<double> LTimeBP1=NC_ReadTimeFromFile(TheHisFile, StringTime);
	DeltaTime=LTimeBP1[0] - LTimeBegin[0];
      }
      else {
	DeltaTime=0;
      }
    }
    // Determination of number of arrays
    int iFileMiddle, nbRecMiddle;
    if (iFileEnd != iFileBegin) {
      iFileMiddle=iFileBegin+1;
      std::string TheHisFile=HisPrefix + StringNumber(iFileMiddle, 4) + ".nc";
      std::vector<double> LTimeMiddle=NC_ReadTimeFromFile(TheHisFile, StringTime);
      nbRecMiddle=LTimeMiddle.size();
    }
    else {
      iFileMiddle=iFileBegin;
      nbRecMiddle=nbRecBegin;
    }
    //    std::cerr << "iFileMiddle=" << iFileMiddle << "\n";
    int nbRecEnd;
    if (iFileEnd != iFileMiddle) {
      std::string TheHisFile=HisPrefix + StringNumber(iFileEnd, 4) + ".nc";
      std::vector<double> LTimeEnd=NC_ReadTimeFromFile(TheHisFile, StringTime);
      nbRecEnd=LTimeEnd.size();
    }
    else {
      nbRecEnd=nbRecMiddle;
    }
    int nbFile=1 + iFileEnd - iFileBegin;
    int NbPerArray[nbFile];
    for (int iFile=0; iFile<nbFile; iFile++)
      NbPerArray[iFile]=nbRecMiddle;
    NbPerArray[0]=nbRecBegin;
    NbPerArray[nbFile-1]=nbRecEnd;
    /*
    for (int iFile=0; iFile<nbFile; iFile++) {
      int nbPer=NbPerArray[iFile];
            std::cerr << "iFile=" << iFile << " nbPer=" << nbPer << "\n";
    }
    */
    double currentTime=LTimeBegin[0];
    FirstTime=currentTime;
    for (int iFile=0; iFile<nbFile; iFile++) {
      int iFileTot=iFile+iFileBegin;
      std::string TheHisFile=HisPrefix + StringNumber(iFileTot, 4) + ".nc";
      ListFileNames.push_back(TheHisFile);
      int siz=NbPerArray[iFile];
      for (int i=0; i<siz; i++) {
	ListIFile.push_back(iFile);
	ListIRec.push_back(i);
	ListTime.push_back(currentTime);
	currentTime=currentTime + DeltaTime;
      }
    }
    LastTime=currentTime;
  }
  ArrayHistory eArr;
  int nbFile=ListFileNames.size();
  int nbTime=ListTime.size();
  eArr.nbFile=nbFile;
  eArr.nbTime=nbTime;
  eArr.FirstTime=FirstTime;
  eArr.LastTime=LastTime;
  eArr.ListFileNames=ListFileNames;
  eArr.ListIFile=ListIFile;
  eArr.ListIRec=ListIRec;
  eArr.ListTime=ListTime;
  eArr.AppendVarName=false;
  eArr.KindArchive="NETCDF";
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
ArrayHistory WW3_ReadArrayHistory(std::string const& HisFile, std::string const& HisPrefix)
{
  std::vector<double> LTime=NC_ReadTimeFromFile(HisFile, "time");
  int nbTime=LTime.size();
  std::vector<int> ListIRec, ListIFile;
  for (int iTime=0; iTime<nbTime; iTime++) {
    ListIFile.push_back(0);
    ListIRec.push_back(iTime);
  }
  double FirstTime=LTime[0];
  double LastTime=LTime[nbTime-1];
  std::vector<std::string> ListFileNames;
  std::cerr << "NC_ReadArrayHistory, HisPrefix=" << HisPrefix << "\n";
  ListFileNames.push_back(HisPrefix);
  ArrayHistory eArr;
  eArr.nbFile=1;
  eArr.nbTime=nbTime;
  eArr.FirstTime=FirstTime;
  eArr.LastTime=LastTime;
  eArr.ListFileNames=ListFileNames;
  std::cerr << "|ListFileNames|=" << ListFileNames.size() << "\n";
  eArr.ListIFile=ListIFile;
  eArr.ListIRec=ListIRec;
  eArr.ListTime=LTime;
  eArr.AppendVarName=true;
  eArr.KindArchive="NETCDF";
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
GridArray PRE_RETRIEVE_GRID_ARRAY(TripleModelDesc const& eTriple)
{
  std::cerr << "PRE_RETRIEVE_GRID_ARRAY, step 1\n";
  std::string eModelName=eTriple.ModelName;
  std::cerr << "PRE_RETRIEVE_GRID_ARRAY, step 2\n";
  std::string GridFile=GET_GRID_FILE(eTriple);
  std::cerr << "eModelName=" << eModelName << "\n";
  if (eModelName == "COSMO") {
    return NC_ReadCosmoWamStructGridFile(GridFile, "atm");
  }
  if (eModelName == "WAM") {
    return NC_ReadWamGridFile(GridFile);
  }
  if (eModelName == "ROMS") {
    return NC_ReadRomsGridFile(GridFile);
  }
  if (eModelName == "WWM") {
    std::string eExtension=FILE_GetExtension(GridFile);
    std::cerr << "eExtension=" << eExtension << "\n";
    if (eExtension == "gr3")
      return WWM_ReadGridFile_gr3(GridFile);
    if (eExtension == "dat")
      return WWM_ReadGridFile_xfn(GridFile);
    if (eExtension == "nc")
      return WWM_ReadGridFile_netcdf(GridFile);
    std::cerr << "Error in reading grid for WWM\n";
    std::cerr << "We did not find the right kind\n";
    exit(1);
  }
  if (eModelName == "WW3") {
    return NC_ReadWW3_GridFile(GridFile);
  }
  if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO") {
    std::string HisPrefix=eTriple.HisPrefix;
    std::vector<std::string> ListFile=GRIB_GetAllFilesInDirectory(HisPrefix);
    if (ListFile.size() == 0) {
      std::cerr << "The list of files is empty\n";
      std::cerr << "Error happened in GRIB_GetAllFilesInDirectory\n";
      exit(1);
    }
    std::string eFileName=ListFile[0];
    GridArray GrdArr=GRIB_ReadGridArray(eFileName, eModelName);
    GrdArr.ModelName=eModelName;
    return GrdArr;
  }
  std::cerr << "Error in PRE_RETRIEVE_GRID_ARRAY\n";
  std::cerr << "Did not find the matching model for the grid\n";
  std::cerr << "Please correct\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
GridArray RETRIEVE_GRID_ARRAY(TripleModelDesc const& eTriple)
{
  std::cerr << "Before PRE_RETRIEVE_GRID_ARRAY\n";
  GridArray GrdArr=PRE_RETRIEVE_GRID_ARRAY(eTriple);
  std::cerr << "After PRE_RETRIEVE_GRID_ARRAY\n";
  if (GrdArr.IsFE == 0)
    return GrdArr;
  if (eTriple.CutWorldMap)
    CutWorldMap(GrdArr);
  std::cerr << "After CutWorldMap\n";
  if (eTriple.HigherLatitudeCut) {
    double MinLatCut=eTriple.MinLatCut;
    double MaxLatCut=eTriple.MaxLatCut;
    CUT_HigherLatitude(GrdArr, MinLatCut, MaxLatCut);
  }
  std::cerr << "After CUT_HigherLatitude\n";
  CHECK_COORDINATE_ORIENTATION(GrdArr);
  return GrdArr;
}
/* ------------------------------ */


/* ------------------------------ */
ArrayHistory NC_ReadArrayHistory(TripleModelDesc const& eTriple)
{
  std::string StringTime="ocean_time";
  std::string eModelName=eTriple.ModelName;
  std::string HisPrefix=eTriple.HisPrefix;
  if (eModelName != "WW3")
    return NC_ReadArrayHistory_Kernel(HisPrefix, StringTime);
  std::string HisFile=GET_GRID_FILE(eTriple);
  return WW3_ReadArrayHistory(HisFile, HisPrefix);
}
/* ------------------------------ */


/* ------------------------------ */
ArrayHistory GRIB_ReadArrayHistory(std::string const& HisPrefix)
{
  std::cerr << "Beginning of GRIB_ReadArrayHistory\n";
  std::vector<std::string> ListFile=GRIB_GetAllFilesInDirectory(HisPrefix);
  int nbFile=ListFile.size();
  std::cerr << "nbFile=" << nbFile << "\n";
  std::vector<GRIB_MessageInfo> TotalListMessage;
  struct PairTime {
    double time;
    double timeStart;
  };
  std::vector<PairTime> ListPairTime;
  double MaxErrorTime=0.01;
  std::function<void(PairTime)> fInsert=[&](PairTime const& ePair) -> void {
    int len=ListPairTime.size();
    //    std::cerr << "fInsert time=" << ePair.time << " timeStart=" << ePair.timeStart << "\n";
    for (int i=0; i<len; i++) {
      //      std::cerr << "i=" << i << " ListPairTime[i],time=" << ListPairTime[i].time << " timeStart=" << ListPairTime[i].timeStart << "\n";
      if (fabs(ePair.time - ListPairTime[i].time) < MaxErrorTime) {
	//	std::cerr << "Pass the test\n";
	double deltaTime=ListPairTime[i].time - ListPairTime[i].timeStart;
	double deltaTimeNew=ePair.time - ePair.timeStart;
	//	std::cerr << "  deltaTime=" << deltaTime << " deltaTimeNew=" << deltaTimeNew << "\n";
	if (deltaTimeNew < deltaTime) {
	  ListPairTime[i]=ePair;
	  //	  std::cerr << "Update of ListPairTime\n";
	}
	return;
      }
    }
    //    std::cerr << "One insertion\n";
    ListPairTime.push_back(ePair);
  };
  std::function<int(PairTime)> fPosition=[&](PairTime const& ePair) -> int {
    int len=ListPairTime.size();
    for (int i=0; i<len; i++)
      if (fabs(ePair.time - ListPairTime[i].time) < MaxErrorTime && fabs(ePair.timeStart - ListPairTime[i].timeStart) < MaxErrorTime)
	return i;
    return -1;
  };
  std::vector<GRIB_MessageInfo> ListSingleEntry;
  for (auto & eFile : ListFile) {
    //    std::cerr << "Before getting ListMessage eFile=" << eFile << "\n";
    std::vector<GRIB_MessageInfo> ListMessage=GRIB_GetAllListPairTime(eFile);
    //    std::cerr << "After getting ListMessage\n";
    int nbMessage=ListMessage.size();
    if (nbMessage == 0) {
      std::cerr << "nbMessage=0\n";
      std::cerr << "No messages in the file\n";
      std::cerr << "eFile = " << eFile << "\n";
      exit(1);
    }
    ListSingleEntry.push_back(ListMessage[0]);
    bool IsFirst=true;
    double time=0;
    double timeStart=0;
    double InconsTime=0;
    double InconsTimeStart=0;
    for (auto & eMesg : ListMessage) {
      if (IsFirst == false) {
	InconsTime += fabs(time - eMesg.time);
	InconsTimeStart += fabs(timeStart - eMesg.timeStart);
      }
      else {
	time=eMesg.time;
	timeStart=eMesg.timeStart;
      }
      IsFirst=false;
    }
    if (InconsTime > MaxErrorTime || InconsTimeStart > MaxErrorTime) {
      std::cerr << "Error in time for the GRIB.\n";
      std::cerr << "We are working under the assumption that the grib files\n";
      std::cerr << "contain a single time and that hypothesis is not satisfied\n";
      std::cerr << "InconsTime=" << InconsTime << "\n";
      std::cerr << "InconsTimeStart=" << InconsTimeStart << "\n";
      exit(1);
    }
    fInsert({time, timeStart});
  }
  int nbTime=ListPairTime.size();
  /*
  std::cerr << "Before sorting\n";
  for (int iTime=0; iTime<nbTime; iTime++) {
    std::cerr << "iTime=" << iTime << "/" << nbTime << " time=" << ListPairTime[iTime].time << " timestart=" << ListPairTime[iTime].timeStart << "\n";
    }*/
  sort(ListPairTime.begin(), ListPairTime.end(), 
       [&](PairTime const& a, PairTime const& b) -> bool {
	 if (a.time < b.time)
	   return true;
	 return false;
       });
  /*
  std::cerr << "After sorting\n";
  for (int iTime=0; iTime<nbTime; iTime++) {
    std::cerr << "iTime=" << iTime << "/" << nbTime << " time=" << ListPairTime[iTime].time << " timestart=" << ListPairTime[iTime].timeStart << "\n";
    }*/


  /*
  for (int i=1; i<nbTime; i++) {
    double eDiff=ListPairTime[i].time - ListPairTime[i-1].time;
    std::cerr << "i=" << i << " diff=" << eDiff << "\n";
  }
  */
  std::vector<double> ListTime(nbTime);
  for (int i=0; i<nbTime; i++)
    ListTime[i]=ListPairTime[i].time;
  std::vector<std::vector<std::string> > ListListFileNames(nbTime);
  for (int iFile=0; iFile<nbFile; iFile++) {
    std::string eFile=ListFile[iFile];
    GRIB_MessageInfo eMesg=ListSingleEntry[iFile];
    int pos=fPosition({eMesg.time, eMesg.timeStart});
    if (pos != -1)
      ListListFileNames[pos].push_back(eFile);
  }
  ArrayHistory eArr;
  eArr.nbFile=nbFile;
  eArr.nbTime=nbTime;
  //  eArr.FirstTime=FirstTime;
  //  eArr.LastTime=LastTime;
  eArr.ListListFileNames=ListListFileNames;
  eArr.ListTime=ListTime;
  eArr.KindArchive="GRIB";
  std::cerr << "Ending of GRIB_ReadArrayHistory\n";
  //  exit(1);
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
ArrayHistory ReadArrayHistory(TripleModelDesc const& eTriple)
{
  ArrayHistory eArr;
  std::string HisPrefix=eTriple.HisPrefix;
  std::string eModelName=eTriple.ModelName;
  CHECK_Model_Allowedness(eModelName);
  if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_COSMO" || eModelName == "GRIB_ECMWF") {
    eArr=GRIB_ReadArrayHistory(HisPrefix);
  }
  else {
    eArr=NC_ReadArrayHistory(eTriple);
  }
  return eArr;
}
/* ------------------------------ */


/* ------------------------------ */
MyMatrix<double> Get2DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& VarName, double const& eTimeDay)
{
  if (TotalArr.eArr.KindArchive == "NETCDF") {
    return NETCDF_Get2DvariableSpecTime(TotalArr, VarName, eTimeDay);
  }
  if (TotalArr.eArr.KindArchive == "GRIB") {
    return GRIB_Get2DvariableSpecTime(TotalArr, VarName, eTimeDay);
  }
  std::cerr << "The KindArchive does not allow to find the nature\n";
  std::cerr << "KindArchive=" << TotalArr.eArr.KindArchive << "\n";
  exit(1);
}
/* ------------------------------ */


/* ------------------------------ */
RecVar ModelSpecificVarSpecificTime(TotalArrGetData const& TotalArr, std::string const& eVarName, double const& eTimeDay)
{
  std::string eModelName=TotalArr.GrdArr.ModelName;
  RecVar eRecVar;
  std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDay);
  std::string strFile=DATE_ConvertMjd2mystringFile(eTimeDay);
  eRecVar.eTimeDay=eTimeDay;
  eRecVar.strPres=strPres;
  eRecVar.strFile=strFile;
  eRecVar.VarNature="rho";
  MyMatrix<double> F;
  MyMatrix<double> Fret;
  MyMatrix<double> U;
  MyMatrix<double> V;
  eRecVar.VarName1=eVarName;
  eRecVar.VarName2="unset";
  if (eVarName == "WIND10") {
    if (eModelName == "ROMS" || eModelName == "WWM") {
      U=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      U=Get2DvariableSpecTime(TotalArr, "U_10", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "V_10", eTimeDay);
    }
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO") {
      U=Get2DvariableSpecTime(TotalArr, "10u", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "10v", eTimeDay);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr);
    eRecVar.VarName2="10m wind";
    eRecVar.minval=0;
    eRecVar.maxval=13;
    eRecVar.mindiff=-2;
    eRecVar.maxdiff=2;
    eRecVar.Unit="m/s";
    eRecVar.VarNature="uv";
    eRecVar.nameU="Uwind";
    eRecVar.nameV="Vwind";
  }
  if (eVarName == "SurfCurr") {
    if (eModelName == "ROMS") {
      ThreeDimArray<double> Utot=Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
      ThreeDimArray<double> Vtot=Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
      int s_rho=Utot.GetDims()[0];
      MyMatrix<double> Usurf=Utot.DimensionExtraction(0, s_rho-1);
      MyMatrix<double> Vsurf=Vtot.DimensionExtraction(0, s_rho-1);
      U=My_u2rho(Usurf, TotalArr.GrdArr.GrdArrRho.MSK);
      V=My_v2rho(Vsurf, TotalArr.GrdArr.GrdArrRho.MSK);
    }
    if (eModelName == "WWM") {
      U=Get2DvariableSpecTime(TotalArr, "CURTX", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "CURTY", eTimeDay);
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      U=Get2DvariableSpecTime(TotalArr, "ucurr", eTimeDay);
      V=Get2DvariableSpecTime(TotalArr, "vcurr", eTimeDay);
    }
    AngleRhoRot(U, V, TotalArr.GrdArr);
    eRecVar.VarName2="surface current";
    eRecVar.minval=0;
    eRecVar.maxval=0.5;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="m/s";
    eRecVar.VarNature="uv";
    eRecVar.nameU="UsurfCurr";
    eRecVar.nameV="VsurfCurr";
  }
  if (eVarName == "Hwave") {
    if (eModelName == "WWM")
      F=Get2DvariableSpecTime(TotalArr, "HS", eTimeDay);
    if (eModelName == "WW3")
      F=Get2DvariableSpecTime(TotalArr, "hs", eTimeDay);
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "Hwave", eTimeDay);
    eRecVar.VarName2="Significant wave height";
    eRecVar.minval=0;
    eRecVar.maxval=4.5;
    eRecVar.mindiff=-0.5;
    eRecVar.maxdiff=0.5;
    eRecVar.Unit="m";
  }
  if (eVarName == "WINDMAG") {
    if (eModelName == "ROMS" || eModelName == "WWM") {
      if (NC_TOTALARR_IsVar(TotalArr, "Uwind") == true) {
	MyMatrix<double> U=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
	MyMatrix<double> V=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
	F=COMPUTE_NORM(U, V);
      }
      else {
	F=Get2DvariableSpecTime(TotalArr, "WNDMAG", eTimeDay);
      }
    }
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> U=Get2DvariableSpecTime(TotalArr, "U_10", eTimeDay);
      MyMatrix<double> V=Get2DvariableSpecTime(TotalArr, "V_10", eTimeDay);
      F=COMPUTE_NORM(U, V);
    }
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO") {
      MyMatrix<double> U=Get2DvariableSpecTime(TotalArr, "10u", eTimeDay);
      MyMatrix<double> V=Get2DvariableSpecTime(TotalArr, "10v", eTimeDay);
      F=COMPUTE_NORM(U, V);
    }
    eRecVar.VarName2="10m wind speed";
    eRecVar.minval=0;
    eRecVar.maxval=13;
    eRecVar.mindiff=-2;
    eRecVar.maxdiff=2;
    eRecVar.Unit="m/s";
  }
  if (eVarName == "TempSurf") {
    if (eModelName == "ROMS") {
      ThreeDimArray<double> TheTemp=Get3DvariableSpecTime(TotalArr, "temp", eTimeDay);
      int s_rho=TheTemp.GetDims()[0];
      F=TheTemp.DimensionExtraction(0, s_rho-1);
    }
    if (eModelName == "COSMO") {
      F=Get2DvariableSpecTime(TotalArr, "t_s", eTimeDay);
      int siz=F.size();
      for (int i=0; i<siz; i++)
	F(i) -= double(273.15);
    }
    eRecVar.VarName2="sea surface temperature";
    eRecVar.minval=10;
    eRecVar.maxval=20;
    eRecVar.mindiff=-2;
    eRecVar.maxdiff=2;
    eRecVar.Unit="deg";
  }
  if (eVarName == "SaltSurf") {
    if (eModelName == "ROMS") {
      ThreeDimArray<double> TheSalt=Get3DvariableSpecTime(TotalArr, "salt", eTimeDay);
      int s_rho=TheSalt.GetDims()[0];
      F=TheSalt.DimensionExtraction(0, s_rho-1);
    }
    eRecVar.VarName2="sea surface salinity";
    eRecVar.minval=30;
    eRecVar.maxval=40;
    eRecVar.mindiff=-2;
    eRecVar.maxdiff=2;
    eRecVar.Unit="PSU";
  }
  if (eVarName == "AIRT2") {
    if (eModelName == "COSMO")
      F=Get2DvariableSpecTime(TotalArr, "t_2m", eTimeDay);
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_GFS" || eModelName == "GRIB_COSMO")
      F=Get2DvariableSpecTime(TotalArr, "2t", eTimeDay);
    eRecVar.VarName2="2m air temperature";
    eRecVar.minval=10;
    eRecVar.maxval=20;
    eRecVar.mindiff=-2;
    eRecVar.maxdiff=2;
    eRecVar.Unit="deg";
  }
  if (eVarName == "Rh2") {
    if (eModelName == "COSMO")
      F=Get2DvariableSpecTime(TotalArr, "rh_2m", eTimeDay);
    if (eModelName == "GRIB_DWD")
      F=Get2DvariableSpecTime(TotalArr, "RELHUM_2M", eTimeDay);
    eRecVar.VarName2="2m relative humidity";
    eRecVar.minval=0;
    eRecVar.maxval=100;
    eRecVar.mindiff=-20;
    eRecVar.maxdiff=20;
    eRecVar.Unit="nondimensional";
  }
  if (eVarName == "ZetaOcean") {
    if (eModelName == "COSMO")
      F=Get2DvariableSpecTime(TotalArr, "ZetaOcean", eTimeDay);
    eRecVar.VarName2="free surface elevation";
    eRecVar.minval=-0.2;
    eRecVar.maxval=0.2;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="m";
  }
  if (eVarName == "MwaveFreq") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "MwaveFreq", eTimeDay);
    eRecVar.VarName2="mean wave frequency";
    eRecVar.minval=0;
    eRecVar.maxval=0.9;
    eRecVar.mindiff=-0.2;
    eRecVar.maxdiff=0.2;
    eRecVar.Unit="Hz";
  }
  if (eVarName == "PwaveFreq") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "PwaveFreq", eTimeDay);
    eRecVar.VarName2="peak wave frequency";
    eRecVar.minval=0;
    eRecVar.maxval=0.9;
    eRecVar.mindiff=-0.2;
    eRecVar.maxdiff=0.2;
    eRecVar.Unit="Hz";
  }
  if (eVarName == "MwavePer") {
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "MwaveFreq", eTimeDay);
      F=FreqPeriodChange(Fin);
    }
    eRecVar.VarName2="mean wave period";
    eRecVar.minval=2;
    eRecVar.maxval=10;
    eRecVar.mindiff=-1;
    eRecVar.maxdiff=1;
    eRecVar.Unit="s";
  }
  if (eVarName == "PwavePer") {
    if (eModelName == "COSMO" || eModelName == "WAM") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "PwaveFreq", eTimeDay);
      F=FreqPeriodChange(Fin);
    }
    eRecVar.VarName2="peak wave period";
    eRecVar.minval=2;
    eRecVar.maxval=10;
    eRecVar.mindiff=-1;
    eRecVar.maxdiff=1;
    eRecVar.Unit="s";
  }
  if (eVarName == "AIRD") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "AIRD", eTimeDay);
    eRecVar.VarName2="air density";
    eRecVar.minval=1.12;
    eRecVar.maxval=1.20;
    eRecVar.mindiff=-0.02;
    eRecVar.maxdiff=0.02;
    eRecVar.Unit="kg/m3";
  }
  if (eVarName == "CdWave") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "CdWave", eTimeDay);
    eRecVar.VarName2="drag coefficient from the wave model";
    eRecVar.minval=0.10;
    eRecVar.maxval=0.20;
    eRecVar.mindiff=-0.05;
    eRecVar.maxdiff=0.05;
    eRecVar.Unit="nondimensional";
  }
  if (eVarName == "AlphaWave") {
    if (eModelName == "COSMO" || eModelName == "WAM")
      F=Get2DvariableSpecTime(TotalArr, "AlphaWave", eTimeDay);
    eRecVar.VarName2="Charnock coefficient from the wave model";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="nondimensional";
  }
  if (eVarName == "rain") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "rain", eTimeDay);
    eRecVar.VarName2="rain of all types";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="kg/m2s";
  }
  if (eVarName == "swrad") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "swrad", eTimeDay);
    eRecVar.VarName2="Shortwave flux";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="W/m2";
  }
  if (eVarName == "lwrad") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "lwrad", eTimeDay);
    eRecVar.VarName2="Longwave flux";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="W/m2";
  }
  if (eVarName == "latent") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "latent", eTimeDay);
    eRecVar.VarName2="Latent flux";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="W/m2";
  }
  if (eVarName == "SurfPres") {
    if (eModelName == "ROMS") {
      MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "Pair", eTimeDay);
      F=100*Fin;
    }
    if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS")
      F=Get2DvariableSpecTime(TotalArr, "prmsl", eTimeDay);
    if (eModelName == "GRIB_ECMWF")
      F=Get2DvariableSpecTime(TotalArr, "msl", eTimeDay);
    if (eModelName == "GRIB_COSMO")
      F=Get2DvariableSpecTime(TotalArr, "pmsl", eTimeDay);
    eRecVar.VarName2="mean sea level pressure";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="W/m2";
  }
  if (eVarName == "sensible") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "sensible", eTimeDay);
    eRecVar.VarName2="Sensible heat flux";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="W/m2";
  }
  if (eVarName == "shflux") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "shflux", eTimeDay);
    eRecVar.VarName2="Surface heat flux";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="W/m2";
  }
  if (eVarName == "ssflux") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "ssflux", eTimeDay);
    eRecVar.VarName2="Surface salinity flux";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="PSU/m2s";
  }
  if (eVarName == "evaporation") {
    if (eModelName == "ROMS")
      F=Get2DvariableSpecTime(TotalArr, "evaporation", eTimeDay);
    eRecVar.VarName2="Evaporation rate";
    eRecVar.minval=0;
    eRecVar.maxval=0.033;
    eRecVar.mindiff=-0.1;
    eRecVar.maxdiff=0.1;
    eRecVar.Unit="kg/m2s";
  }
  //
  // Now error parsing and assignations
  // 
  if (eRecVar.VarName2 == "unset") {
    std::cerr << "We did not find the variable\n";
    std::cerr << "eVarName = " << eVarName << "\n";
    std::cerr << "in the list of allowed ones\n";
    std::cerr << "possibly missspelling or lack of relevant code\n";
    exit(1);
  }
  if (eModelName != "TRIVIAL") {
    if (eRecVar.VarNature == "rho") {
      if (F.size() == 0) {
	std::cerr << "Variable eVarName = " << eVarName << "\n";
	std::cerr << "is recognized by the program\n";
	std::cerr << "But it has not been assigned.\n";
	std::cerr << "Possibly because of missing facility for\n";
	std::cerr << "eModelName = " << eModelName << "\n";
	exit(1);
      }
      eRecVar.F=F;
    }
    if (eRecVar.VarNature == "uv") {
      if (U.size() == 0 || V.size() == 0) {
	std::cerr << "Variable eVarName = " << eVarName << "\n";
	std::cerr << "is recognized by the program\n";
	std::cerr << "But it has not been assigned.\n";
	std::cerr << "Possibly because of missing facility for\n";
	std::cerr << "eModelName = " << eModelName << "\n";
	exit(1);
      }
      eRecVar.U=U;
      eRecVar.V=V;
    }
  }
  return eRecVar;
}
/* ------------------------------ */


/* ------------------------------ */
RecVar ModelSpecificVarSpecificTime(TotalArrGetData const& TotalArr, std::string const& eVarName, double const& eTimeDay, PlotBound const& ePlotBound)
{
  RecVar eRecVar=ModelSpecificVarSpecificTime(TotalArr, eVarName, eTimeDay);
  //
  // Setting up bounds for the plots.
  //
  int nbSingle=ePlotBound.BoundSingle_var.size();
  int nbSingleMin=ePlotBound.BoundSingle_min.size();
  int nbSingleMax=ePlotBound.BoundSingle_max.size();
  //  std::cerr << "nbSingle=" << nbSingle << "\n";
  if (nbSingle != nbSingleMin || nbSingle != nbSingleMax) {
    std::cerr << "Number of entries in BoundSingle_var, BoundSingle_min, BoundSingle_max\n";
    std::cerr << "Should all be the same. Now,\n";
    std::cerr << "nbSingle    = " << nbSingle << "\n";
    std::cerr << "nbSingleMin = " << nbSingleMin << "\n";
    std::cerr << "nbSingleMax = " << nbSingleMax << "\n";
  }
  for (int iS=0; iS<nbSingle; iS++)
    if (ePlotBound.BoundSingle_var[iS] == eVarName) {
      eRecVar.minval=ePlotBound.BoundSingle_min[iS];
      eRecVar.maxval=ePlotBound.BoundSingle_max[iS];
    }
  int nbDiff=ePlotBound.BoundDiff_var.size();
  int nbDiffMin=ePlotBound.BoundDiff_min.size();
  int nbDiffMax=ePlotBound.BoundDiff_max.size();
  //  std::cerr << "nbDiff=" << nbDiff << "\n";
  if (nbDiff != nbDiffMin || nbDiff != nbDiffMax) {
    std::cerr << "Number of entries in BoundDiff_var, BoundDiff_min, BoundDiff_max\n";
    std::cerr << "Should all be the same. Now,\n";
    std::cerr << "nbDiff    = " << nbDiff << "\n";
    std::cerr << "nbDiffMin = " << nbDiffMin << "\n";
    std::cerr << "nbDiffMax = " << nbDiffMax << "\n";
  }
  for (int iD=0; iD<nbDiff; iD++)
    if (ePlotBound.BoundDiff_var[iD] == eVarName) {
      eRecVar.mindiff=ePlotBound.BoundDiff_min[iD];
      eRecVar.maxdiff=ePlotBound.BoundDiff_max[iD];
    }
  return eRecVar;
}
/* ------------------------------ */


/* ------------------------------ */
RecVar RetrieveTrivialRecVar(std::string const& eVarName)
{
  TotalArrGetData TotalArrTrivial;
  TotalArrTrivial.GrdArr.ModelName="TRIVIAL";
  double eTimeDayTrivial=0;
  return ModelSpecificVarSpecificTime(TotalArrTrivial, eVarName, eTimeDayTrivial);
}
/* ------------------------------ */


/* ------------------------------ */
RecVar ModelSpecificVarSpecificTimeGeneral(TotalArrGetData const& TotalArr, std::string const& eVarName, VarQuery const& eQuery, PlotBound const& ePlotBound)
{
  //
  // Check correctness 
  //
  std::vector<std::string> ListAllow{"instant", "average", "swathMax", "swathMin"};
  if (std::find(ListAllow.begin(), ListAllow.end(), eQuery.NatureQuery) == ListAllow.end()) {
    std::cerr << "We failed to find NatureQuery=" << eQuery.NatureQuery << "\n";
    std::cerr << "List of allowed queries:\n";
    for (auto & eStr : ListAllow)
      std::cerr << "  eStr=" << eStr << "\n";
    exit(1);
  }
  std::string strPres=DATE_ConvertMjd2mystringPres(eQuery.eTimeDay);
  std::cerr << "Query ModelSpecificVarSpecificTimeGeneral NatureQuery=" << eQuery.NatureQuery << " date=" << strPres << " VarName=" << eVarName << "\n";
  //
  // Reading array and doing operations "instant", "average", etc.
  //
  RecVar eRecVar;
  if (eQuery.NatureQuery == "instant") {
    eRecVar=ModelSpecificVarSpecificTime(TotalArr, eVarName, eQuery.eTimeDay, ePlotBound);
  }
  else {
    double eTimeDay=eQuery.eTimeDay;
    double TimeFrameDay=eQuery.TimeFrameDay;
    std::vector<int> ListRelITime=GetIntervalListITime(TotalArr.eArr.ListTime, eTimeDay, TimeFrameDay);
    int nbTimeRel=ListRelITime.size();
    RecVar RecVarTrivial=RetrieveTrivialRecVar(eVarName);
    MyMatrix<double> F, U, V;
    for (int iTimeRel=0; iTimeRel<nbTimeRel; iTimeRel++) {
      int iTime=ListRelITime[iTimeRel];
      double eTimeDay=TotalArr.eArr.ListTime[iTime];
      eRecVar=ModelSpecificVarSpecificTime(TotalArr, eVarName, eTimeDay, ePlotBound);
      if (iTimeRel == 0) {
	if (RecVarTrivial.VarNature == "rho") {
	  F=eRecVar.F;
	}
	else {
	  U=eRecVar.U;
	  V=eRecVar.V;
	}
      }
      else {
	if (eQuery.NatureQuery == "average") {
	  if (RecVarTrivial.VarNature == "rho") {
	    //	  std::cerr << "Before += sum(F)=" << F.sum() << "\n";
	    F += eRecVar.F;
	    //	  std::cerr << " After += sum(F)=" << F.sum() << "\n";
	  }
	  if (RecVarTrivial.VarNature == "uv") {
	    U += eRecVar.U;
	    V += eRecVar.V;
	  }
	}
	if (eQuery.NatureQuery == "swathMax") {
	  if (RecVarTrivial.VarNature == "rho") {
	    F=F.cwiseMax(eRecVar.F);
	  }
	  if (RecVarTrivial.VarNature == "uv") {
	    std::cerr << "swathMax for uv does not have any sense\n";
	    exit(1);
	  }
	}
	if (eQuery.NatureQuery == "swathMin") {
	  if (RecVarTrivial.VarNature == "rho") {
	    F=F.cwiseMin(eRecVar.F);
	  }
	  if (RecVarTrivial.VarNature == "uv") {
	    std::cerr << "swathMin for uv does not have any sense\n";
	    exit(1);
	  }
	}
      }
    }
    if (eQuery.NatureQuery == "average") {
      if (RecVarTrivial.VarNature == "rho") {
	F /= double(nbTimeRel);
      }
      if (RecVarTrivial.VarNature == "uv") {
	U /= double(nbTimeRel);
	V /= double(nbTimeRel);
      }
    }
    if (RecVarTrivial.VarNature == "rho") {
      eRecVar.F=F;
    }
    if (RecVarTrivial.VarNature == "uv") {
      eRecVar.U=U;
      eRecVar.V=V;
    }
  }
  std::string strAll=GetStrAllOfPlot(eQuery);
  eRecVar.strAll=strAll;
  return eRecVar;
}
/* ------------------------------ */


/* ------------------------------ */
PairRecVar ModelPairSpecificVarSpecificTimeGeneral(TotalArrGetData const& TotalArr1, TotalArrGetData const& TotalArr2, std::string const& eVarName, VarQuery const& eQuery, PlotBound const& ePlotBound)
{
  //  std::cerr << "NatureQuery=" << eQuery.NatureQuery << "\n";
  std::vector<std::string> ListAllow{"instant", "average", "swathMax", "swathMin"};
  RecVar eRecVar1, eRecVar2;
  if (std::find(ListAllow.begin(), ListAllow.end(), eQuery.NatureQuery) != ListAllow.end()) {
    eRecVar1=ModelSpecificVarSpecificTimeGeneral(TotalArr1, eVarName, eQuery, ePlotBound);
    eRecVar2=ModelSpecificVarSpecificTimeGeneral(TotalArr2, eVarName, eQuery, ePlotBound);
  }
  else {
    std::vector<std::string> ListAllowB{"MaxDiff", "MinDiff"};
    if (std::find(ListAllowB.begin(), ListAllowB.end(), eQuery.NatureQuery) == ListAllowB.end()) {
      std::cerr << "Only two difference operator are allowed: MaxDiff and MinDiff\n";
      std::cerr << "NatureQuery=" << eQuery.NatureQuery << "\n";
    }
    RecVar RecVarTrivial=RetrieveTrivialRecVar(eVarName);
    //
    double eTimeDay=eQuery.eTimeDay;
    double TimeFrameDay=eQuery.TimeFrameDay;
    std::vector<int> ListRelITime1=GetIntervalListITime(TotalArr1.eArr.ListTime, eTimeDay, TimeFrameDay);
    std::vector<int> ListRelITime2=GetIntervalListITime(TotalArr2.eArr.ListTime, eTimeDay, TimeFrameDay);
    int nbTimeRel1=ListRelITime1.size();
    int nbTimeRel2=ListRelITime2.size();
    double DeltaTime1=(TotalArr1.eArr.ListTime[nbTimeRel1-1] - TotalArr1.eArr.ListTime[0])/double(nbTimeRel1-1);
    double DeltaTime2=(TotalArr2.eArr.ListTime[nbTimeRel2-1] - TotalArr2.eArr.ListTime[0])/double(nbTimeRel2-1);
    double DeltaTime=T_min(DeltaTime1, DeltaTime2);
    double FirstTime=eTimeDay;
    double LastTime=eTimeDay + TimeFrameDay;
    std::vector<double> ListRelTimeDay=GetInterval(FirstTime, LastTime, DeltaTime);
    int nbTimeRel=ListRelTimeDay.size();
    MyMatrix<double> diffF;
    if (RecVarTrivial.VarNature == "uv") {
      std::cerr << "Cannot do MaxDiff or MinDiff for variables of type uv\n";
      exit(1);
    }
    for (int iTimeRel=0; iTimeRel<nbTimeRel; iTimeRel++) {
      double eTimeDay=ListRelTimeDay[iTimeRel];
      eRecVar1=ModelSpecificVarSpecificTime(TotalArr1, eVarName, eTimeDay, ePlotBound);
      eRecVar2=ModelSpecificVarSpecificTime(TotalArr2, eVarName, eTimeDay, ePlotBound);
      if (iTimeRel == 0) {
	if (RecVarTrivial.VarNature == "rho") {
	  diffF = eRecVar1.F - eRecVar2.F;
	}
      }
      else {
	if (eQuery.NatureQuery == "MaxDiff") {
	  if (RecVarTrivial.VarNature == "rho")
	    diffF=diffF.cwiseMax(eRecVar1.F - eRecVar2.F);
	}
	if (eQuery.NatureQuery == "MinDiff") {
	  if (RecVarTrivial.VarNature == "rho")
	    diffF=diffF.cwiseMin(eRecVar1.F - eRecVar2.F);
	}
      }
    }
    if (RecVarTrivial.VarNature == "rho") {
      eRecVar1.F=diffF;
      eRecVar2.F.fill(double(0));
    }
  }
  std::string strAll=GetStrAllOfPlot(eQuery);
  eRecVar1.strAll=strAll;
  eRecVar2.strAll=strAll;
  return {eRecVar1, eRecVar2};
}
/* ------------------------------ */


/* ------------------------------ */
void SINGLE_Plotting_Function(FullNamelist const& eFull)
{
  //
  // Creating the triple for data reading (grid and history)
  // 
  std::cerr << "SINGLE_Plotting_Function, step 0\n";
  std::map<std::string, SingleBlock> ListBlock=eFull.ListBlock;
  SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
  std::string eModelName=eBlPROC.ListStringValues.at("MODELNAME");
  std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
  std::string HisPrefix=eBlPROC.ListStringValues.at("HisPrefix");
  bool CutWorldMap=eBlPROC.ListBoolValues.at("CutWorldMap");
  bool HigherLatitudeCut=eBlPROC.ListBoolValues.at("HigherLatitudeCut");
  double MinLatCut=eBlPROC.ListDoubleValues.at("MinLatCut");
  double MaxLatCut=eBlPROC.ListDoubleValues.at("MaxLatCut");
  TripleModelDesc eTriple{eModelName, GridFile, HisPrefix, CutWorldMap, HigherLatitudeCut, MinLatCut, MaxLatCut};
  //
  // Retrieving the grid array
  //
  GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
  //
  // Setting up the timings.
  //
  ArrayHistory eArr=ReadArrayHistory(eTriple);
  std::vector<double> ListTime=GetInterval(eBlPROC.ListStringValues.at("BEGTC"),
					   eBlPROC.ListStringValues.at("ENDTC"), 
					   eBlPROC.ListDoubleValues.at("DELTC"), 
					   eBlPROC.ListStringValues.at("UNITC"));
  int nbTime=ListTime.size();
  std::cerr << "nbTime=" << nbTime << "\n";
  //
  TotalArrGetData TotalArr{GrdArr, eArr};
  //
  // Setting up DrawArr
  // It contains the permanent feature of plots that will not change from variable
  // to variable and time to time
  //
  PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
  //
  // Preliminary drawings 
  //
  GRID_PLOTTING(TotalArr, GridFile, ePerm, eFull);
  //
  // Now the major time loop
  //
  SingleBlock eBlockVAR=eFull.ListBlock.at("VARS");
  std::vector<std::string> ListVarOut=GetAllPossibleVariables();
  VarQuery eQuery;
  std::vector<std::string> ListNatureQuery=eBlPROC.ListListStringValues.at("ListNatureQuery");
  eQuery.TimeFrameDay=eBlPROC.ListDoubleValues.at("TimeFrameDay");
  eQuery.NaturePlot="single";
  PlotBound ePlotBound=ReadPlotBound(eFull);
  for (int iTime=0; iTime<nbTime; iTime++) {
    double eTimeDay=ListTime[iTime];
    eQuery.iTime=iTime;
    eQuery.eTimeDay=eTimeDay;
    std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDay);
    std::cerr << "iTime=" << iTime << "/" << nbTime << " date=" << strPres << "\n";
    for (auto & eVarName : ListVarOut) {
      if (eBlockVAR.ListBoolValues.at(eVarName)) {
	for (auto & eNatureQuery : ListNatureQuery) {
	  eQuery.NatureQuery=eNatureQuery;
	  RecVar eRecVar=ModelSpecificVarSpecificTimeGeneral(TotalArr, eVarName, eQuery, ePlotBound);
	  GENERAL_PLOT_SINGLE(GrdArr, eRecVar, ePerm, eFull);
	}
      }
    }
  }
}
/* ------------------------------ */


/* ------------------------------ */
int main(int argc, char *argv[])
{
  FullNamelist eFull=NAMELIST_GetStandard_PlotRoutine_single();
  if (argc != 2) {
    fprintf(stderr, "PLOT_results is used as\n");
    fprintf(stderr, "PLOT_results [file.nml]\n");
    fprintf(stderr, "with file.nml the file describing the plotting routines\n");
    return -1;
  }
  std::string eFileName=argv[1];
  NAMELIST_ReadNamelistFile(eFileName, eFull);
  SINGLE_Plotting_Function(eFull);
}
/* ------------------------------ */


