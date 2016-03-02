# 1 "/home/mathieu/GIT/wwmIII/CppOcean/AltimeterComparison.cpp"
#include <netcdf>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <unsupported/Eigen/CXX11/Tensor>
template <typename T> using MyVector = Eigen::Matrix<T,Eigen::Dynamic,1>;
template <typename T> using MyMatrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
template <typename T> using MySparseMatrix = Eigen::SparseMatrix<T,Eigen::ColMajor>;
#include <ctype.h>
#include <malloc.h>
#include <unistd.h>
#include <getopt.h>
#include <chrono>
#include <ctime>
#include <math.h>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <functional>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
typedef unsigned long ulong;
typedef unsigned int uint;
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include "grib_api.h"
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
struct SingleEntryAltimeter {
double Time, Lon, Lat;
double WindSpeed, WindSpeed_cor, WindSpeed_used;
double Sigma0, Sigma0_cal, Sigma0std, Sigma0second, Sigma0secondstd;
double Swh, SwhStd, Swh_cor, Swh_used;
double Swh_model, WindSpeed_model;
int Satellite;
double DistToCoast;
};
struct SingleSearchEntry {
size_t iEntry;
int eEta;
int eXi;
double TimeCoeff;
double SpatialCoeff;
};
struct SmoothArr {
int LenTotal;
std::vector<int> ListShift;
std::vector<double> ListWeight;
};
SmoothArr GetSmoothingArray(double const& avgDistKM_model, double const& avgDistKM_track)
{
double TheSizeReal=avgDistKM_model / avgDistKM_track;
double w=(TheSizeReal-1)/2;
int wlow=int(floor(w));
int LenTotal=1+2*wlow + 2;
std::vector<int> ListShift;
std::vector<double> ListWeight;
double eOne=1.0;
ListShift.push_back(0);
ListWeight.push_back(eOne);
for (int widx=1; widx<=wlow; widx++) {
ListShift.push_back(widx);
ListWeight.push_back(eOne);
ListShift.push_back(-widx);
ListWeight.push_back(eOne);
}
double wres=(TheSizeReal-eOne-double(2*wlow))/double(2);
ListShift.push_back(wlow+1);
ListWeight.push_back(wres);
ListShift.push_back(-1-wlow);
ListWeight.push_back(wres);
return {LenTotal, ListShift, ListWeight};
}
std::vector<SingleEntryAltimeter> SpatialAveragingTrack(std::vector<SingleEntryAltimeter> const& eVectEnt, SmoothArr const& eSmoothArr)
{
int len=eVectEnt.size();
std::vector<SingleEntryAltimeter> eVectRet(len);
int LenTotal=eSmoothArr.LenTotal;
for (int i=0; i<len; i++) {
SingleEntryAltimeter eRefEnt=eVectEnt[i];
double Time=eRefEnt.Time;
double Lon=eRefEnt.Lon;
double Lat=eRefEnt.Lat;
double WindSpeed_model=eRefEnt.WindSpeed_model;
double Swh_model=eRefEnt.Swh_model;
double SumSwh=0;
double SumSwh_cor=0;
double SumSwh_used=0;
double SumWindSpeed=0;
double SumWindSpeed_cor=0;
double SumWindSpeed_used=0;
double SumWeight=0;
for (int iShift=0; iShift<LenTotal; iShift++) {
int iNew=i + eSmoothArr.ListShift[iShift];
double eWeight=eSmoothArr.ListWeight[iShift];
if (iNew >=0 && iNew < len) {
auto eRef=eVectEnt[iNew];
SumWeight += eWeight;
SumSwh += eRef.Swh;
SumSwh_cor += eRef.Swh_cor;
SumSwh_used += eRef.Swh_used;
SumWindSpeed += eRef.WindSpeed;
SumWindSpeed_cor += eRef.WindSpeed_cor;
SumWindSpeed_used += eRef.WindSpeed_used;
}
}
SingleEntryAltimeter eNewEnt;
eNewEnt.Time=Time;
eNewEnt.Lon=Lon;
eNewEnt.Lat=Lat;
eNewEnt.WindSpeed_model=WindSpeed_model;
eNewEnt.Swh_model=Swh_model;
eNewEnt.Swh=SumSwh/SumWeight;
eNewEnt.Swh_cor=SumSwh_cor/SumWeight;
eNewEnt.Swh_used=SumSwh_used/SumWeight;
eNewEnt.WindSpeed=SumWindSpeed/SumWeight;
eNewEnt.WindSpeed_cor=SumWindSpeed_cor/SumWeight;
eNewEnt.WindSpeed_used=SumWindSpeed_used/SumWeight;
eVectRet[i]=eNewEnt;
}
return eVectRet;
}
std::vector<SingleEntryAltimeter> SelectByStatus(std::vector<std::vector<int> > const& ListListStatus, std::vector<SingleEntryAltimeter> const& eVectEnt)
{
int siz=eVectEnt.size();
std::vector<SingleEntryAltimeter> RetVect;
int nbStatus=ListListStatus.size();
for (int i=0; i<siz; i++) {
int eStatus=1;
for (int iStatus=0; iStatus<nbStatus; iStatus++)
eStatus=eStatus*ListListStatus[iStatus][i];
if (eStatus == 1)
RetVect.push_back(eVectEnt[i]);
}
return RetVect;
}
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
struct PairMM {
double Meas;
double Model;
};
void Print_Down_Statistics(std::ostream & os, std::string const& eName, T_stat const& eStat)
{
os << "        eName=" << eName << "\n";
os << "       nbMeas=" << eStat.nbMeas << "\n";
int nbMeas=eStat.nbMeas;
if (nbMeas >0) {
os << "      MaxMeas=" << eStat.MaxMeas << "\n";
os << "      MinMeas=" << eStat.MinMeas << "\n";
os << "     MaxModel=" << eStat.MaxModel << "\n";
os << "     MinModel=" << eStat.MinModel << "\n";
os << "     MeanMeas=" << eStat.MeanMeas << "\n";
os << "    MeanModel=" << eStat.MeanModel << "\n";
os << "    MeanError=" << eStat.MeanError << "\n";
os << "AbsoluteError=" << eStat.AbsoluteError << "\n";
os << "         RMSE=" << eStat.RMSE << "\n";
os << " CenteredRMSE=" << eStat.CenteredRMSE << "\n";
os << "  Correlation=" << eStat.Correlation << "\n";
os << " ScatterIndex=" << eStat.ScatterIndex << "\n";
os << "  CenteredSci=" << eStat.CenteredScatterIndex << "\n";
}
}
bool IsPointInside(double const& testx, double const& testy, std::vector<double> const& vertx, std::vector<double> const& verty)
{
int nvert=vertx.size();
int i, j;
bool c=false;
for (i = 0, j = nvert-1; i < nvert; j = i++) {
if ( ((verty[i]>testy) != (verty[j]>testy)) &&
(testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i] ) + vertx[i]) )
c = !c;
}
return c;
}
struct PairCoord {
int i;
int j;
};
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
double GeodesicDistanceKM(double const& LonDeg1, double const& LatDeg1, double const& LonDeg2, double const& LatDeg2)
{
double EarthRadius=6370;
return EarthRadius*GeodesicDistance(LonDeg1, LatDeg1, LonDeg2, LatDeg2);
}
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
int MySign(double & TheVal)
{
if (TheVal > 0)
return 1;
if (TheVal < 0)
return -1;
return 0;
}
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
ANG_rho(iEta,iXi) = ThePi/double(2) - azim(iEta-1,iXi);
for (int iXi=0; iXi<xi_v; iXi++) {
ANG_rho(0,iXi) = ANG_rho(1,iXi);
ANG_rho(eta_rho-1,iXi) = ANG_rho(eta_rho-2,iXi);
}
return ANG_rho;
}
MyMatrix<double> My_u2rho(MyMatrix<double> const& eVar_u, MyMatrix<int> const& MSK_u)
{
int eta_u=eVar_u.rows();
int xi_u=eVar_u.cols();
int eta_rho = eta_u;
int xi_rho = xi_u + 1;
MyMatrix<double> eVar_rho(eta_rho, xi_rho);
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
int eSumMsk=0;
double eSumVal=0;
if (j<xi_u) {
if (MSK_u(i,j) == 1) {
eSumMsk++;
eSumVal += eVar_u(i,j);
}
}
if (j > 0) {
if (MSK_u(i,j-1) == 1) {
eSumMsk++;
eSumVal += eVar_u(i,j-1);
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
MyMatrix<double> My_v2rho(MyMatrix<double> const& eVar_v, MyMatrix<int> const& MSK_v)
{
int eta_v=eVar_v.rows();
int xi_v=eVar_v.cols();
int xi_rho = xi_v;
int eta_rho = eta_v + 1;
MyMatrix<double> eVar_rho(eta_rho, xi_rho);
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
int eSumMsk=0;
double eSumVal=0;
if (i < eta_v) {
if (MSK_v(i,j) == 1) {
eSumMsk++;
eSumVal += eVar_v(i,j);
}
}
if (i > 0) {
if (MSK_v(i-1,j) == 1) {
eSumMsk++;
eSumVal += eVar_v(i-1,j);
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
Eigen::Tensor<double,3> My_u2rho_3D(Eigen::Tensor<double,3> const& eVar_u, MyMatrix<int> const& MSK_u)
{
auto LDim=eVar_u.dimensions();
int s_vert=LDim[0];
int eta_u=LDim[1];
int xi_u=LDim[2];
int eta_rho = eta_u;
int xi_rho = xi_u + 1;
std::vector<double> VertColumn(s_vert);
Eigen::Tensor<double,3> eVar_rho(s_vert, eta_rho, xi_rho);
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
int eSumMsk=0;
for (int k=0; k<s_vert; k++)
VertColumn[k]=0;
if (j<xi_u) {
if (MSK_u(i,j) == 1) {
eSumMsk++;
for (int k=0; k<s_vert; k++)
VertColumn[k] += eVar_u(k,i,j);
}
}
if (j > 0) {
if (MSK_u(i,j-1) == 1) {
eSumMsk++;
for (int k=0; k<s_vert; k++)
VertColumn[k] += eVar_u(k,i,j-1);
}
}
if (eSumMsk == 0) {
for (int k=0; k<s_vert; k++)
eVar_rho(k,i,j)=0;
}
else {
for (int k=0; k<s_vert; k++) {
double eVal=VertColumn[k]/double(eSumMsk);
eVar_rho(k,i,j)=eVal;
}
}
}
return eVar_rho;
}
Eigen::Tensor<double,3> My_v2rho_3D(Eigen::Tensor<double,3> const& eVar_v, MyMatrix<int> const& MSK_v)
{
auto LDim=eVar_v.dimensions();
int s_vert=LDim[0];
int eta_v=LDim[1];
int xi_v=LDim[2];
int xi_rho = xi_v;
int eta_rho = eta_v + 1;
std::vector<double> VertColumn;
Eigen::Tensor<double,3> eVar_rho(s_vert,eta_rho, xi_rho);
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
int eSumMsk=0;
for (int k=0; k<s_vert; k++)
VertColumn[k]=0;
if (i < eta_v) {
if (MSK_v(i,j) == 1) {
eSumMsk++;
for (int k=0; k<s_vert; k++)
VertColumn[k] += eVar_v(k,i,j);
}
}
if (i > 0) {
if (MSK_v(i-1,j) == 1) {
eSumMsk++;
for (int k=0; k<s_vert; k++)
VertColumn[k] += eVar_v(k,i-1,j);
}
}
if (eSumMsk == 0) {
for (int k=0; k<s_vert; k++)
eVar_rho(k,i,j)=0;
}
else {
for (int k=0; k<s_vert; k++) {
double eVal=VertColumn[k]/double(eSumMsk);
eVar_rho(k,i,j)=eVal;
}
}
}
return eVar_rho;
}
std::vector<std::string> GetAllPossibleModels()
{
std::vector<std::string> vec{"COSMO", "WAM", "ROMS", "ROMS_IVICA", "WWM", "WWM_DAILY", "WW3", "GRIB_DWD", "GRIB_ECMWF", "GRIB_GFS", "GRIB_COSMO", "GRIB_WAM_FORT30", "SCHISM_SFLUX"};
return vec;
}
int TheSignFct(double const& eVal)
{
if (eVal > 0)
return 1;
if (eVal < 0)
return -1;
return 0;
}
void DifferenceLonRenormalize(double & Lon)
{
if (Lon > 180)
Lon=Lon - 360;
if (Lon < -180)
Lon=Lon + 360;
}
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
std::vector<std::string> GetAllPossibleVariables()
{
std::vector<std::string> ListVarOut{
"IOBPWW3", "MAPSTA", "FieldOut1", "CFL1", "CFL2", "CFL3", "ThreeDfield1", "NbIterSolv",
"WIND10", "Uwind", "Vwind","WINDMAG",
"SurfCurr", "UsurfCurr", "VsurfCurr", "SurfCurrMag",
"Curr", "CurrMag",
"TempSurf", "SaltSurf", "AIRT2", "AIRT2K", "Rh2", "Rh2frac", "AIRD", "SurfPres",
"ZetaOcean", "ZetaOceanDerivative", "DynBathy", "ZetaSetup",
"CdWave", "AlphaWave", "AirZ0", "AirFricVel",
"shflux", "ssflux", "evaporation",
"Hwave", "BreakingFraction",
"rain", "swrad", "lwrad", "latent", "sensible",
"MeanWaveFreq", "PeakWaveFreq", "TM02",
"MeanWavePer", "PeakWavePer",
"MeanWaveDirSpread", "PeakWaveDirSpread",
"MeanWaveDir", "PeakWaveDir", "MeanWaveDirVect", "PeakWaveDirVect",
"DiscPeakWaveDir",
"MeanWaveLength", "PeakWaveLength", "MeanWaveNumber", "PeakWaveNumber",
"TotSurfStr", "WaveSurfStr", "SurfStrHF"};
return ListVarOut;
}
struct PairMinMax {
double TheMin;
double TheMax;
};
void PrintMyScriptSubtitle(std::ostream &os)
{
os << "procedure subtitles(wks:graphic,plot:graphic,lstr:string,cstr:string, \\\n";
os << "                    rstr:string,tres)\n";
os << "local txres, font_height, amres\n";
os << "begin\n";
os << "  if(tres) then\n";
os << "    txres = tres     ; Copy resources\n";
os << "  else\n";
os << "    txres = True\n";
os << "  end if\n";
os << ";\n";
os << "; Retrieve font height of left axis string and use to\n";
os << "; calculate size of subtitles.\n";
os << ";\n";
os << "  if(.not.isatt(txres,\"txFontHeightF\")) then\n";
os << "    getvalues plot\n";
os << "      \"tiXAxisFontHeightF\" : font_height\n";
os << "    end getvalues\n";
os << "    txres@txFontHeightF = font_height*0.9\n";
os << "  end if\n";
os << ";\n";
os << "; Set some some annotation resources.\n";
os << ";\n";
os << "  amres                  = True\n";
os << "  if(.not.isatt(txres,\"amOrthogonalPosF\")) then\n";
os << "    amres@amOrthogonalPosF = -0.53   ; Top of plot plus a little extra\n";
os << "                                     ; to stay out of the tickmarks.\n";
os << "  else\n";
os << "    amres@amOrthogonalPosF = txres@amOrthogonalPosF\n";
os << "  end if\n";
os << ";\n";
os << "; Create three strings to put at the top, using a slightly\n";
os << "; smaller font height than the axis titles.\n";
os << ";\n";
os << "  if(lstr.ne.\"\") then\n";
os << "    txidl = gsn_create_text(wks, lstr, txres)\n";
os << "    amres@amJust           = \"BottomLeft\"\n";
os << "    amres@amParallelPosF   = -0.5   ; Left-justified\n";
os << "    annoidl = gsn_add_annotation(plot, txidl, amres)\n";
os << "  end if\n";
os << "  if(cstr.ne.\"\") then\n";
os << "    txidc = gsn_create_text(wks, cstr, txres)\n";
os << "    amres@amJust           = \"BottomCenter\"\n";
os << "    amres@amParallelPosF   = 0.0   ; Centered\n";
os << "    annoidc = gsn_add_annotation(plot, txidc, amres)\n";
os << "  end if\n";
os << "  if(rstr.ne.\"\") then\n";
os << "    txidr = gsn_create_text(wks, rstr, txres)\n";
os << "    amres@amJust           = \"BottomRight\"\n";
os << "    amres@amParallelPosF   = 0.5   ; Right-justifed\n";
os << "    annoidr = gsn_add_annotation(plot, txidr, amres)\n";
os << "  end if\n";
os << "end\n";
}
std::vector<std::string> GetAllNamesOfSatelliteAltimeter()
{
return {"ERS1", "ERS2", "ENVISAT", "TOPEX", "POSEIDON", "JASON1", "GFO", "JASON2", "CRYOSAT", "SARAL"};
}
std::string GetNameOfSatelliteAltimeter(int iSat)
{
return GetAllNamesOfSatelliteAltimeter()[iSat-1];
}
struct SinglePartInterp {
int eEta, eXi;
double eCoeff;
};
struct SingleRecInterp {
bool status;
std::vector<SinglePartInterp> LPart;
};
struct QuadCoordinate {
double MinLon;
double MaxLon;
double MinLat;
double MaxLat;
};
std::vector<double> DetermineCoefficient(std::vector<double> const& X, std::vector<double> const& Y, double const& eX, double const& eY)
{
MyMatrix<double> A(3,3);
for (int i=0; i<3; i++) {
A(0,i)=1;
A(1,i)=X[i];
A(2,i)=Y[i];
}
MyVector<double> V(3);
V(0)=1;
V(1)=eX;
V(2)=eY;
MyMatrix<double> eInv=A.inverse().eval();
MyVector<double> eProduct=eInv*V;
std::vector<double> LCoeff(3);
for (int i=0; i<3; i++) {
LCoeff[i]=eProduct(i);
}
# 818 "/home/mathieu/GIT/wwmIII/CppOcean/AltimeterComparison.cpp"
return LCoeff;
}
bool TestFeasibilityByQuad(QuadCoordinate const& eQuad, double const& eLon, double const& eLat)
{
if (eLon > eQuad.MinLon && eLon < eQuad.MaxLon && eLat > eQuad.MinLat && eLat < eQuad.MaxLat)
return true;
return false;
}
MyMatrix<int> GetDirection()
{
MyMatrix<int> eMat(2,4);
eMat(0,0)=0;
eMat(1,0)=0;
eMat(0,1)=1;
eMat(1,1)=0;
eMat(0,2)=0;
eMat(1,2)=1;
eMat(0,3)=1;
eMat(1,3)=1;
return eMat;
}
template<typename T>
T VectorSum(MyVector<T> const& eVect)
{
T eSum=0;
int siz=eVect.size();
for (int i=0; i<siz; i++)
eSum += eVect(i);
return eSum;
}
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
template<typename T>
void SwapValues(T& val1, T& val2)
{
T prov;
prov=val1;
val1=val2;
val2=prov;
}
template<typename T>
struct Inverse_exception {
std::string errmsg;
T pivot;
};
# 958 "/home/mathieu/GIT/wwmIII/CppOcean/AltimeterComparison.cpp"
template<typename T>
struct SelectionRowCol {
int TheRank;
MyMatrix<T> NSP;
std::vector<int> ListColSelect;
std::vector<int> ListRowSelect;
};
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
std::vector<int> ListColSelect01(nbCol,0);
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
ZeroAssignation(NSP);
nbVect=0;
for (iCol=0; iCol<nbCol; iCol++)
if (ListColSelect01[iCol] == 0) {
eVal=1;
NSP(nbVect, iCol)=eVal;
for (iRank=0; iRank<eRank; iRank++) {
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
template<typename T>
MyMatrix<T> SelectRow(MyMatrix<T> const&TheMat, std::vector<int> const& eList)
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
template<typename T>
MyMatrix<T> SelectColumn(MyMatrix<T> const& TheMat, std::vector<int> const& eList)
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
template<typename T>
MyVector<T> SelectColumnVector(MyVector<T> const& TheV, std::vector<int> const& eList)
{
int nbColRed=eList.size();
MyVector<T> TheProv(nbColRed);
for (int iCol=0; iCol<nbColRed; iCol++) {
int jCol=eList[iCol];
TheProv(iCol)=TheV(jCol);
}
return TheProv;
}
template<typename T>
struct SolMatResult {
bool result;
MyVector<T> eSol;
};
template<typename T1, typename T2>
MyMatrix<T2> ConvertMatrix(MyMatrix<T1> const& M, std::function<T2(const T1&)> const& f)
{
int eta_rho=M.rows();
int xi_rho=M.cols();
MyMatrix<T2> eRet(eta_rho, xi_rho);
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
T1 eVal1=M(i,j);
T2 eVal2=f(eVal1);
eRet(i,j)=eVal2;
}
return eRet;
}
template<typename T>
void WriteMatrixGAP(std::ostream &os, MyMatrix<T> const&TheMat)
{
long nbRow=TheMat.rows();
long nbCol=TheMat.cols();
os << "[ ";
for (long iRow=0; iRow<nbRow; iRow++) {
if (iRow > 0)
os << ",\n";
os << "[";
for (long iCol=0; iCol<nbCol; iCol++) {
T eVal=TheMat(iRow, iCol);
if (iCol > 0)
os << ",";
os << " " << eVal;
}
os << " ]";
}
os << " ]";
}
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
std::string DoubleToString(double const& x)
{
std::stringstream s;
s << x;
std::string converted(s.str());
return converted;
}
std::string IntToString(int const & x)
{
std::stringstream s;
s << x;
std::string converted(s.str());
return converted;
}
std::string LongToString(long const & x)
{
std::stringstream s;
s << x;
std::string converted(s.str());
return converted;
}
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
TheStr += "0";
}
TheStr += IntToString(nb);
return TheStr;
}
idx++;
}
}
std::string STRING_RemoveSpacesBeginningEnd(std::string const& eStr)
{
int len=eStr.size();
std::vector<int> ListIsSpace(len,0);
std::string eSpace=" ";
for (int i=0; i<len; i++) {
std::string eChar=eStr.substr(i, 1);
if (eChar == eSpace)
ListIsSpace[i]=1;
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
RetStr += eStr.at(iPos);
}
return RetStr;
}
std::vector<std::string> STRING_Split(std::string const& eStrA, std::string const& eStrB)
{
int lenA=eStrA.length();
int lenB=eStrB.length();
std::vector<int> ListStatus(lenA,1);
for (int iA=0; iA<lenA - lenB; iA++)
if (ListStatus[iA] == 1) {
bool IsMatch=true;
for (int iB=0; iB<lenB; iB++) {
std::string eCharA=eStrA.substr(iA+iB,1);
std::string eCharB=eStrB.substr(iB,1);
if (eCharA != eCharB)
IsMatch=false;
}
if (IsMatch)
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
if (siz > 0)
RetList.push_back(eFound);
eFound="";
}
}
int siz=eFound.size();
if (siz > 0)
RetList.push_back(eFound);
return RetList;
}
std::string FILE_GetExtension(std::string const& eFile)
{
std::vector<std::string> LStr=STRING_Split(eFile, "/");
std::string eFinal=LStr[LStr.size()-1];
std::vector<std::string> LBlck=STRING_Split(eFile, ".");
return LBlck[LBlck.size()-1];
}
struct TerminalException {
int eVal;
};
template <typename T>
struct is_ring_field {
static const bool value = false;
};
template <>
struct is_ring_field<short> {
static const bool value = false;
};
template <>
struct is_ring_field<long> {
static const bool value = false;
};
template <>
struct is_ring_field<int> {
static const bool value = false;
};
template <>
struct is_ring_field<long long> {
static const bool value = false;
};
template <>
struct is_ring_field<double> {
static const bool value = true;
};
template <>
struct is_ring_field<float> {
static const bool value = true;
};
template<typename T>
T VectorSum(std::vector<T> const& eVect)
{
T eSum=0;
for (T eVal : eVect)
eSum += eVal;
return eSum;
}
template<typename T>
T VectorMin(std::vector<T> const& eVect)
{
T eMin=eVect[0];
for (T eVal : eVect)
if (eVal < eMin)
eMin=eVal;
return eMin;
}
template<typename T>
T VectorMax(std::vector<T> const& eVect)
{
T eMax=eVect[0];
for (T eVal : eVect)
if (eVal > eMax)
eMax=eVal;
return eMax;
}
template<typename T>
T T_abs(T const& eVal)
{
if (eVal > 0)
return eVal;
T fVal= - eVal;
return fVal;
}
template<typename T>
T T_min(T const& eVal1, T const& eVal2)
{
if (eVal1 > eVal2)
return eVal2;
return eVal1;
}
std::string random_string( size_t length )
{
srand ( time(NULL) );
auto randchar = []() -> char {
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
template<typename T>
void WriteStdVector(std::ostream& os, std::vector<T> const& V)
{
for (auto & eVal : V)
os << " " << eVal;
}
template<typename T>
void WriteStdVectorGAP(std::ostream& os, std::vector<T> const& V)
{
os << "[";
bool IsFirst=true;
for (auto & eVal : V) {
if (IsFirst == false)
os << ",";
IsFirst=false;
os << eVal;
}
os << "]";
}
template<typename T>
struct CollectedResult {
std::vector<T> LVal;
std::vector<int> LMult;
};
template<typename T>
CollectedResult<T> Collected(std::vector<T> const& eVect)
{
std::set<T> SetVal;
for (auto & eVal : eVect)
SetVal.insert(eVal);
std::vector<T> LVal;
for (auto & eVal : SetVal)
LVal.push_back(eVal);
int eSize=LVal.size();
std::vector<int> LMult(eSize,0);
auto UpPosition=[&](T const& eVal) -> void {
for (int i=0; i<eSize; i++)
if (LVal[i] == eVal) {
LMult[i] += 1;
return;
}
std::cerr << "Should never reach that stage\n";
throw TerminalException{1};
};
for (auto & eVal : eVect)
UpPosition(eVal);
return {LVal, LMult};
}
int NextIdx(int const& len,int const& i)
{
if (i == len-1)
return 0;
return i+1;
}
int PrevIdx(int const& len,int const& i)
{
if (i == 0)
return len-1;
return i-1;
}
struct PairLL {
double eLon;
double eLat;
};
struct RecSymbolic {
double eTimeDay;
int iTime;
std::string strPres;
std::string strFile;
std::string strAll;
std::string VarName1;
std::string VarName2;
double minval;
double maxval;
double mindiff;
double maxdiff;
std::string Unit;
std::string VarNature;
std::string nameU, nameV;
};
struct RecVar {
RecSymbolic RecS;
MyMatrix<double> U;
MyMatrix<double> V;
MyMatrix<double> F;
Eigen::Tensor<double,3> Uthree;
Eigen::Tensor<double,3> Vthree;
Eigen::Tensor<double,3> Tens3;
};
struct CoordGridArrayFD {
int eta, xi;
MyMatrix<int> MSK;
MyMatrix<double> LON, LAT, DEP, ANG;
int nbWet;
bool HaveDEP;
std::vector<int> Idx, Jdx;
};
struct GridArray {
std::string ModelName;
int IsFE;
bool IsSpherical;
CoordGridArrayFD GrdArrRho, GrdArrU, GrdArrV, GrdArrPsi;
MyMatrix<int> INE;
MyVector<int> IOBP;
bool L_IndexSelect;
std::vector<int> I_IndexSelect;
};
struct GRIB_MessageInfo {
std::string shortName;
std::string name;
int idx;
double time;
double timeStart;
int stepRange;
std::string FileName;
};
struct ArrayHistory {
std::string KindArchive;
int nbFile, nbTime;
double FirstTime, LastTime;
std::string FirstTimeStr, LastTimeStr;
std::vector<std::string> ListFileNames;
std::vector<std::vector<GRIB_MessageInfo> > ListListMessages;
std::vector<GRIB_MessageInfo> ListAllMessage;
std::vector<std::string> RawVarNames;
std::vector<int> ListITime;
std::vector<double> ListStartTime;
std::vector<int> ListIStartTime;
std::vector<int> ListIFile;
std::vector<int> ListIRec;
std::vector<double> ListTime;
std::string TimeSteppingInfo;
std::string HisPrefix;
double SeparationTime;
int nbDigit;
int nbRecBegin;
int nbRecMiddle;
bool AppendVarName;
};
struct TotalArrGetData {
GridArray GrdArr;
ArrayHistory eArr;
};
struct VarQuery {
double eTimeDay;
int iTime;
std::string NatureQuery;
std::string NaturePlot;
double TimeFrameDay;
};
struct PlotBound {
bool VariableRange;
std::vector<std::string> BoundSingle_var;
std::vector<double> BoundSingle_min;
std::vector<double> BoundSingle_max;
std::vector<std::string> BoundDiff_var;
std::vector<double> BoundDiff_min;
std::vector<double> BoundDiff_max;
};
struct QuadArray {
double MinLon;
double MaxLon;
double MinLat;
double MaxLat;
};
std::vector<std::string> FILE_GetDirectoryListFile(std::string const& eDir)
{
std::string ePath=eDir + ".";
DIR* dirp=opendir(ePath.c_str());
if (dirp == NULL) {
std::cerr << "Error in routine FILE_GetDirectoryListFile\n";
std::cerr << "Error in call to opendir\n";
throw TerminalException{1};
}
struct dirent *dp;
std::vector<std::string> ListFile;
while ((dp = readdir(dirp)) != NULL) {
std::string eName=dp->d_name;
if (eName != ".." && eName != ".")
ListFile.push_back(eName);
}
int err=closedir(dirp);
if (err != 0) {
std::cerr << "err=" << err << "\n";
printf("Oh dear, something went wrong with ls! %s\n", strerror(errno));
throw TerminalException{1};
}
return ListFile;
}
bool FILE_IsDirectoryEmpty(std::string const& eDir)
{
std::vector<std::string> TheList = FILE_GetDirectoryListFile(eDir);
if (TheList.size() == 0)
return true;
return false;
}
bool FILE_IsRegularFile(std::string const& eFile)
{
int status;
struct stat st_buf;
status = stat(eFile.c_str(), &st_buf);
if (status != 0) {
std::cerr << "Problem in FILE_IsRegularFile\n";
std::cerr << "Error, errno = " << errno << "\n";
std::cerr << "eFile=" << eFile << "\n";
throw TerminalException{1};
}
if (S_ISREG (st_buf.st_mode)) {
return true;
}
return false;
}
std::vector<std::string> FILE_GetDirectoryFilesRecursively(std::string const& eDir)
{
std::vector<std::string> ListDir={eDir};
std::vector<std::string> ListFile;
while(1) {
std::vector<std::string> NewListDir;
for (auto & fDir : ListDir) {
std::vector<std::string> LocalListFile=FILE_GetDirectoryListFile(fDir);
for (auto & eFile : LocalListFile) {
std::string NewEnt=fDir + eFile;
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
return ListFile;
}
bool IsExistingDirectory(std::string const& ThePrefix)
{
if (0 != access(ThePrefix.c_str(), F_OK)) {
if (ENOENT == errno) {
return false;
}
if (ENOTDIR == errno) {
return false;
}
std::cerr << "Should not happen a priori\n";
throw TerminalException{1};
}
return true;
}
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
void RemoveFile(std::string const& eFile)
{
std::remove(eFile.c_str());
}
void RemoveFileIfExist(std::string const& eFile)
{
if (IsExistingFile(eFile))
RemoveFile(eFile);
}
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
RemoveFile(eFileTot);
}
}
}
}
int FILE_GetNumberLine(std::string const& eFile)
{
int number_of_lines = 0;
std::string line;
std::ifstream myfile(eFile);
while (std::getline(myfile, line))
++number_of_lines;
return number_of_lines;
}
#ifndef WINDOWS
std::string GetCurrentDirectory()
{
int size = pathconf(".", _PC_PATH_MAX);
char *buf;
char *ptr;
if ((buf = (char *)malloc((size_t)size)) != NULL) {
ptr=getcwd(buf, (size_t)size);
if (ptr == NULL && errno != ERANGE) {
std::cerr << "Error while trying to use getcwd\n";
throw TerminalException{1};
}
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
throw TerminalException{1};
}
}
#endif
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
void CreateDirectory(std::string const& eDir)
{
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
mkdir(tmp, S_IRWXU);
*p = '/';
}
mkdir(tmp, S_IRWXU);
}
std::vector<std::string> ls_operation(std::string const& ThePrefix)
{
std::string TmpFile="/tmp/file" + random_string(20);
std::string eOrder="ls " + ThePrefix + " > " + TmpFile;
int iret=system(eOrder.c_str() );
if (iret == -1) {
std::cerr << "Error in ls_operation\n";
std::cerr << "ThePrefix=" << ThePrefix << "\n";
std::cerr << "unable to run the process\n";
throw TerminalException{1};
}
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
a = int(floor((double(14) - double(month))/double(12)));
y = year + 4800 - a;
m = month + 12*a - 3;
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
double DATE_ConvertSix2mjd(std::vector<int> const& eDate)
{
double eJD1=DATE2JD(eDate);
double eJD2=DATE2JD({1858, 11, 17, 0, 0, 0});
double eMJD=eJD1-eJD2;
return eMJD;
}
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
throw TerminalException{1};
}
}
std::string DATE_ConvertSix2mystringPresReduced(std::vector<int> const& Date)
{
if (Date[3] != 0 || Date[4] != 0 || Date[5] != 0) {
return DATE_ConvertSix2mystringPres(Date);
}
try {
int year, month, day;
year=Date[0];
month=Date[1];
day=Date[2];
std::string eTimeStr=StringNumber(year, 4) + "-" +
StringNumber(month, 2) + "-" +
StringNumber(day, 2);
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
throw TerminalException{1};
}
}
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
throw TerminalException{1};
}
}
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
throw TerminalException{1};
}
std::vector<int> JD2DATE(double const& eJD)
{
int year, month, day, hour, min, sec;
int ijd, a, b, c, d, e, m;
double fjd, second;
ijd = int(floor(eJD + 0.5));
a = ijd + 32044;
b = int(floor((double(4)*double(a) + double(3)) / double(146097)));
c = a - int(floor((double(b) * double(146097)) / double(4)));
d = int(floor((double(4)*double(c) + double(3)) / double(1461)));
e = c - int(floor((double(1461)*double(d)) / double(4)));
m = int(floor((double(5) * double(e) + double(2)) / double(153)));
day = e - int(floor((double(153) * double(m) + double(2)) / double(5))) + 1;
month = m + 3 - 12 * int(floor(double(m) / double(10)));
year = b * 100 + d - 4800 + int(floor(double(m) / double(10)));
fjd = eJD - double(ijd) + 0.5;
second = double(86400) * fjd;
hour = int(floor(second/double(3600)));
second = second - double(3600)*double(hour);
min = int(floor(second/double(60)));
sec = int(floor(second - double(60)*min));
int secNear=int(round(second - double(60)*min));
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
double CT2MJD(std::string const& STIME)
{
double XMJD;
std::vector<int> eDate=DATE_ConvertString2six(STIME);
XMJD=DATE_ConvertSix2mjd(eDate);
return XMJD;
}
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
std::string DATE_ConvertMjd2mystringPresReduced(double const& XMJD)
{
std::string STIME;
double XMJD_1858, eMJD;
XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
eMJD = XMJD + XMJD_1858;
std::vector<int> eDate=JD2DATE(eMJD);
STIME=DATE_ConvertSix2mystringPresReduced(eDate);
return STIME;
}
std::vector<int> DATE_ConvertMjd2six(double const& XMJD)
{
double XMJD_1858, eMJD;
XMJD_1858=DATE2JD({1858, 11, 17, 0, 0, 0});
eMJD = XMJD + XMJD_1858;
return JD2DATE(eMJD);
}
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
std::vector<double> GetInterval(double const& FirstTime, double const& LastTime, double const& DeltaInterval)
{
double eTime=FirstTime;
std::vector<double> ListTime;
double tolDay=double(1)/double(10000);
while(1) {
ListTime.push_back(eTime);
eTime=eTime + DeltaInterval;
if (eTime > LastTime + tolDay)
return ListTime;
}
}
std::vector<double> GetInterval(std::string const& BEGTC, std::string const& ENDTC, double const& eInterval, std::string const& UNITC)
{
int IsDone=0;
double eMult = 0;
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
throw TerminalException{1};
}
double DeltaInterval=eInterval*eMult;
double FirstTime=CT2MJD(BEGTC);
double LastTime=CT2MJD(ENDTC);
double tolDay= double(1) / double(10000);
if (LastTime < FirstTime - tolDay) {
std::cerr << "We should have ENDTC >= BEGTC. But instead we have:\n";
std::cerr << "BEGTC = " << BEGTC << "\n";
std::cerr << "ENDTC = " << ENDTC << "\n";
std::cerr << "Please correct\n";
throw TerminalException{1};
}
return GetInterval(FirstTime, LastTime, DeltaInterval);
}
struct InterpInfo {
int iTimeLow, iTimeUpp;
double alphaLow, alphaUpp;
bool UseSingleEntry;
};
InterpInfo GetTimeDifferentiationInfo(std::vector<double> const& LTime, double const& eTimeDay)
{
InterpInfo eInterpInfo;
double tolDay=double(1)/double(1000000);
int nbTime=LTime.size();
eInterpInfo.UseSingleEntry=false;
if (eTimeDay < LTime[0] - tolDay) {
std::cerr << "The asked entry is before the first time\n";
std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
throw TerminalException{1};
}
if (eTimeDay > LTime[nbTime-1] + tolDay) {
std::cerr << "The asked entry is after the last time\n";
std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
throw TerminalException{1};
}
if (nbTime <= 1) {
std::cerr << "We need at least two entries in order to do the time differential\n";
std::cerr << "nbTime=" << nbTime << "\n";
throw TerminalException{1};
}
for (int iTimeUpp=1; iTimeUpp<nbTime; iTimeUpp++) {
int iTimeLow=iTimeUpp-1;
double eTimeUpp=LTime[iTimeUpp];
double eTimeLow=LTime[iTimeLow];
if (eTimeLow - tolDay < eTimeDay && eTimeDay < eTimeUpp + tolDay) {
double alphaLow=(eTimeDay - eTimeUpp)/(eTimeLow - eTimeUpp);
double alphaUpp=(eTimeLow - eTimeDay)/(eTimeLow - eTimeUpp);
eInterpInfo.iTimeLow=iTimeLow;
eInterpInfo.iTimeUpp=iTimeUpp;
eInterpInfo.alphaLow=alphaLow;
eInterpInfo.alphaUpp=alphaUpp;
return eInterpInfo;
}
}
std::cerr << "Failed to find matching record\n";
std::cerr << "Please debug\n";
throw TerminalException{1};
}
InterpInfo GetTimeInterpolationInfo(std::vector<double> const& LTime, double const& eTimeDay)
{
InterpInfo eInterpInfo;
double tolDay=double(1)/double(1000000);
int nbTime=LTime.size();
for (int iTime=0; iTime<nbTime; iTime++) {
if (fabs(LTime[iTime] - eTimeDay) < tolDay) {
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
throw TerminalException{1};
}
if (eTimeDay > LTime[nbTime-1] + tolDay) {
std::cerr << "The asked entry is after the last time\n";
std::cerr << "AskedTime=" << DATE_ConvertMjd2mystringPres(eTimeDay) << "\n";
std::cerr << "FirstTime=" << DATE_ConvertMjd2mystringPres(LTime[0]) << "\n";
std::cerr << " LastTime=" << DATE_ConvertMjd2mystringPres(LTime[nbTime-1]) << "\n";
throw TerminalException{1};
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
throw TerminalException{1};
}
InterpInfo GetTimeInterpolationInfo_infinite(double const& FirstTime, double const& TheSep, double const& eTimeDay)
{
double tolDay=double(1)/double(1000000);
if (eTimeDay < FirstTime - tolDay) {
std::cerr << "Error in GetTimeInterpolationInfo_infinite\n";
std::cerr << "We have FirstTime = " << FirstTime << "\n";
std::cerr << "     and eTimeDay = " << eTimeDay << "\n";
std::cerr << "i.e. eTimeDay < FirstTime\n";
throw TerminalException{1};
}
if (TheSep < 0) {
std::cerr << "We need TheSep > 0\n";
std::cerr << "But we have TheSep = " << TheSep << "\n";
throw TerminalException{1};
}
InterpInfo eInterpInfo;
int iTime=1;
while(1) {
int iTimeUpp = iTime;
int iTimeLow = iTime-1;
double eTimeLow = FirstTime + double(iTimeLow)*TheSep;
double eTimeUpp = FirstTime + double(iTimeUpp)*TheSep;
if (fabs(eTimeLow - eTimeDay) < tolDay) {
eInterpInfo.UseSingleEntry=true;
eInterpInfo.iTimeLow=iTimeLow;
return eInterpInfo;
}
if (fabs(eTimeUpp - eTimeDay) < tolDay) {
eInterpInfo.UseSingleEntry=true;
eInterpInfo.iTimeLow=iTimeUpp;
return eInterpInfo;
}
double alphaLow=(eTimeDay - eTimeUpp)/(eTimeLow - eTimeUpp);
double alphaUpp=(eTimeLow - eTimeDay)/(eTimeLow - eTimeUpp);
if (alphaLow >= 0 && alphaUpp >= 0) {
eInterpInfo.iTimeLow = iTimeLow;
eInterpInfo.iTimeUpp = iTimeUpp;
eInterpInfo.alphaLow = alphaLow;
eInterpInfo.alphaUpp = alphaUpp;
eInterpInfo.UseSingleEntry=false;
return eInterpInfo;
}
iTime++;
if (iTime > 1000000) {
std::cerr << "FirstTime = " << FirstTime << "\n";
std::cerr << "eTimeFay  = " << eTimeDay << "\n";
std::cerr << "iTime     = " << iTime << "\n";
std::cerr << "TheSep    = " << TheSep << "\n";
std::cerr << "Probably a bug in the infinite loop\n";
throw TerminalException{1};
}
}
}
std::vector<int> GetIFileIRec(int const& nbRecBegin, int const& nbRecMiddle, int const& iTime)
{
int iFile=0;
int iRec;
while(1) {
if (iTime < nbRecBegin + iFile*nbRecMiddle) {
if (iFile == 0) {
iRec=iTime;
}
else {
iRec = iTime - nbRecBegin - (iFile-1)*nbRecMiddle;
}
return {iFile, iRec};
}
iFile++;
}
}
std::vector<int> GetIntervalListITime(std::vector<double> const& LTime, double const& eTimeDay, double const& TimeFrameDay)
{
double epsilon=0.0000001;
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
double GetListTimeSeparation(std::vector<double> const& ListTime)
{
int nbTime=ListTime.size();
std::vector<double> ListVal;
std::vector<int> ListNb;
double tolDay = double(1) / double(100000);
auto InsertDiff=[&](double const& eVal) -> void {
int len=ListVal.size();
for (int i=0; i<len; i++) {
if (fabs(eVal - ListVal[i]) < tolDay) {
ListNb[i]++;
return;
}
}
ListVal.push_back(eVal);
ListNb.push_back(1);
};
for (int iTime=1; iTime<nbTime; iTime++) {
double eDiff=ListTime[iTime] - ListTime[iTime-1];
InsertDiff(eDiff);
}
int siz=ListVal.size();
int eNb=0;
double eVal = -1;
for (int i=0; i<siz; i++) {
if (ListNb[i] > eNb) {
eNb = ListNb[i];
eVal = ListVal[i];
}
}
return eVal;
# 2694 "/home/mathieu/GIT/wwmIII/CppOcean/AltimeterComparison.cpp"
}
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
struct CosmoGridInfo {
double latitudeOfSouthernPoleInDegrees;
double longitudeOfSouthernPoleInDegrees;
double angleOfRotationInDegrees;
double longitudeOfFirstGridPointInDegrees;
double latitudeOfFirstGridPointInDegrees;
double longitudeOfLastGridPointInDegrees;
double latitudeOfLastGridPointInDegrees;
double iDirectionIncrementInDegrees;
double jDirectionIncrementInDegrees;
};
double phirot2phi(double const& phirot, double const& rlarot, double const& polphi, double const& pollam, double const& polgam)
{
double pi=3.1415926535;
double eMult=pi/double(180);
double eMultInv=double(180)/pi;
double zsinpol = sin(eMult * polphi);
double zcospol = cos(eMult * polphi);
double zphis = eMult * phirot;
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
double zgam = eMult * polgam;
zarg = zsinpol*sin(zphis) + zcospol*cos(zphis) * ( cos(zrlas)*cos(zgam) - sin(zgam)*sin(zrlas));
}
else {
zarg = zcospol * cos(zphis) * cos(zrlas) + zsinpol * sin(zphis);
}
double phirot2phi = eMultInv * asin(zarg);
return phirot2phi;
}
double rlarot2rla(double const& phirot, double const& rlarot, double const& polphi, double const& pollam, double const& polgam)
{
double pi=3.1415926535;
double eMult=pi/double(180);
double eMultInv=double(180)/pi;
double zsinpol = sin(eMult * polphi);
double zcospol = cos(eMult * polphi);
double zphis = eMult * phirot;
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
double zgam = eMult * polgam;
zarg1 = sin (zlampol) *
(- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))
+ zcospol * sin(zphis))
- cos (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam));
zarg2 = cos (zlampol) *
(- zsinpol*cos(zphis) * (cos(zrlas)*cos(zgam) - sin(zrlas)*sin(zgam))
+ zcospol * sin(zphis))
+ sin (zlampol)*cos(zphis) * (sin(zrlas)*cos(zgam) + cos(zrlas)*sin(zgam));
}
else {
zarg1 = sin (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis) +
zcospol * sin(zphis)) -
cos (zlampol) * sin(zrlas) * cos(zphis);
zarg2 = cos (zlampol) * (-zsinpol * cos(zrlas) * cos(zphis) +
zcospol * sin(zphis)) +
sin (zlampol) * sin(zrlas) * cos(zphis);
}
double rlarot2rla = eMultInv * atan2(zarg1,zarg2);
return rlarot2rla;
}
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
double tolLL = double(1) / double(100000);
if ( fabs(zendlon_tot - zstartlon_tot) < tolLL || fabs(zendlat_tot - zstartlat_tot) < tolLL) {
std::cerr << "Error of consistency in zstartlat / zendlat\n";
throw TerminalException{1};
}
int eta_rho=LON.rows();
int xi_rho=LON.cols();
double pollat= - pollat_sp;
double pollon= pollon_sp - double(180);
double startlon_tot=zstartlon_tot;
double startlat_tot=zstartlat_tot;
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
double eLonR=startlon_tot + double(i)*dlon;
double eLatR=startlat_tot + double(j)*dlat;
double eLat=phirot2phi(eLatR, eLonR, pollat, pollon, polgam);
double eLon=rlarot2rla(eLatR, eLonR, pollat, pollon, polgam);
LON(i,j)=eLon;
LAT(i,j)=eLat;
}
}
GridArray GRIB_ReadGridArray(std::string const& FileName, std::string const& eModelName)
{
grib_handle *h = NULL;
int err;
FILE* in = NULL;
in = fopen(FileName.c_str(),"r");
int nbMessage=0;
while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) {
nbMessage++;
if (err != GRIB_SUCCESS)
GRIB_CHECK(err,0);
long Ni, Nj, numberOfDataPoints;
GRIB_CHECK(grib_get_long(h,"Ni",&Ni),0);
GRIB_CHECK(grib_get_long(h,"Nj",&Nj),0);
GRIB_CHECK(grib_get_long(h,"numberOfDataPoints",&numberOfDataPoints),0);
CosmoGridInfo eCosmoGrid;
if (eModelName == "GRIB_COSMO") {
double latitudeOfSouthernPoleInDegrees, longitudeOfSouthernPoleInDegrees, angleOfRotationInDegrees;
double latitudeOfFirstGridPointInDegrees, longitudeOfFirstGridPointInDegrees, latitudeOfLastGridPointInDegrees, longitudeOfLastGridPointInDegrees;
double iDirectionIncrementInDegrees, jDirectionIncrementInDegrees;
GRIB_CHECK(grib_get_double(h,"latitudeOfSouthernPoleInDegrees",&latitudeOfSouthernPoleInDegrees),0);
eCosmoGrid.latitudeOfSouthernPoleInDegrees=latitudeOfSouthernPoleInDegrees;
GRIB_CHECK(grib_get_double(h,"longitudeOfSouthernPoleInDegrees",&longitudeOfSouthernPoleInDegrees),0);
eCosmoGrid.longitudeOfSouthernPoleInDegrees=longitudeOfSouthernPoleInDegrees;
GRIB_CHECK(grib_get_double(h,"angleOfRotationInDegrees",&angleOfRotationInDegrees),0);
eCosmoGrid.angleOfRotationInDegrees=angleOfRotationInDegrees;
GRIB_CHECK(grib_get_double(h,"latitudeOfFirstGridPointInDegrees",&latitudeOfFirstGridPointInDegrees),0);
eCosmoGrid.latitudeOfFirstGridPointInDegrees=latitudeOfFirstGridPointInDegrees;
GRIB_CHECK(grib_get_double(h,"longitudeOfFirstGridPointInDegrees",&longitudeOfFirstGridPointInDegrees),0);
eCosmoGrid.longitudeOfFirstGridPointInDegrees=longitudeOfFirstGridPointInDegrees;
GRIB_CHECK(grib_get_double(h,"latitudeOfLastGridPointInDegrees",&latitudeOfLastGridPointInDegrees),0);
eCosmoGrid.latitudeOfLastGridPointInDegrees=latitudeOfLastGridPointInDegrees;
GRIB_CHECK(grib_get_double(h,"longitudeOfLastGridPointInDegrees",&longitudeOfLastGridPointInDegrees),0);
eCosmoGrid.longitudeOfLastGridPointInDegrees=longitudeOfLastGridPointInDegrees;
GRIB_CHECK(grib_get_double(h,"iDirectionIncrementInDegrees",&iDirectionIncrementInDegrees),0);
eCosmoGrid.iDirectionIncrementInDegrees=iDirectionIncrementInDegrees;
GRIB_CHECK(grib_get_double(h,"jDirectionIncrementInDegrees",&jDirectionIncrementInDegrees),0);
eCosmoGrid.jDirectionIncrementInDegrees=jDirectionIncrementInDegrees;
}
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
MyMatrix<int> MSK(eta_rho, xi_rho);
int idx=0;
for (int j=0; j<xi_rho; j++)
for (int i=0; i<eta_rho; i++) {
LON(i,j)=lons[idx];
LAT(i,j)=lats[idx];
MSK(i,j)=1;
idx++;
}
if (eModelName == "GRIB_COSMO") {
Apply_COSMO_Transformation(LON, LAT, eCosmoGrid);
}
std::cerr << "GRIB: [0,0]            lon=" << LON(0,0) << " lat=" << LAT(0,0) << "\n";
std::cerr << "      [eta_rho,0]      lon=" << LON(eta_rho-1,0) << " lat=" << LAT(eta_rho-1,0) << "\n";
std::cerr << "      [eta_rho,xi_rho] lon=" << LON(eta_rho-1,xi_rho-1) << " lat=" << LAT(eta_rho-1,xi_rho-1) << "\n";
std::cerr << "      [0,xi_rho]       lon=" << LON(0,xi_rho-1) << " lat=" << LAT(0,xi_rho-1) << "\n";
MyMatrix<double> ANG=CreateAngleMatrix(LON, LAT);
GridArray GrdArr;
GrdArr.ModelName=eModelName;
GrdArr.IsFE=0;
GrdArr.IsSpherical=true;
GrdArr.GrdArrRho.LON=LON;
GrdArr.GrdArrRho.LAT=LAT;
GrdArr.GrdArrRho.MSK=MSK;
GrdArr.GrdArrRho.ANG=ANG;
GrdArr.GrdArrRho.HaveDEP=false;
free(lats);
free(lons);
free(values);
return GrdArr;
}
std::cerr << "Failed to find the variable\n";
throw TerminalException{1};
}
std::vector<GRIB_MessageInfo> GRIB_GetAllListPairTime(std::string const& FileName)
{
grib_handle *h = NULL;
int err;
FILE* in = NULL;
in = fopen(FileName.c_str(),"r");
unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;
std::vector<GRIB_MessageInfo> ListInfo;
int idx=0;
while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) {
if (err != GRIB_SUCCESS)
GRIB_CHECK(err,0);
char name_space[3]="ls";
grib_keys_iterator* kiter=NULL;
kiter=grib_keys_iterator_new(h,key_iterator_filter_flags,name_space);
if (!kiter) {
printf("ERROR: Unable to create keys iterator\n");
throw TerminalException{1};
}
std::string ShortNameValue;
std::string NameValue;
std::string DataDateValue = "unset";
std::string DataTimeValue = "unset";
std::string StepRangeValue;
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
}
grib_keys_iterator_delete(kiter);
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
throw TerminalException{1};
}
stepRange=stoi(LStr[1]);
}
}
int year, month, day;
if (DataDateValue == "unset") {
long dataDate;
GRIB_CHECK(grib_get_long(h,"dataDate",&dataDate),0);
day=dataDate % 100;
int res1=(dataDate - day)/100;
month=res1 % 100;
int res2=(res1 - month)/100;
year=res2;
}
else {
std::string yearStr=DataDateValue.substr(0,4);
std::string monthStr=DataDateValue.substr(4,2);
std::string dayStr=DataDateValue.substr(6,2);
year=stoi(yearStr);
month=stoi(monthStr);
day=stoi(dayStr);
}
if (DataTimeValue == "unset") {
long dataTime;
GRIB_CHECK(grib_get_long(h,"dataTime",&dataTime),0);
if (dataTime == 0) {
DataTimeValue="0000";
}
else {
DataTimeValue=LongToString(dataTime);
}
}
grib_handle_delete(h);
std::string HourStr=DataTimeValue.substr(0,2);
std::string MinStr=DataTimeValue.substr(2,2);
int hour=stoi(HourStr);
int min=stoi(MinStr);
double eTimeStart=DATE_ConvertSix2mjd({year, month, day, hour, min, 0});
double eTime=eTimeStart + double(stepRange)/double(24);
GRIB_MessageInfo eInfo{ShortNameValue, NameValue, idx, eTime, eTimeStart, stepRange, FileName};
ListInfo.push_back(eInfo);
idx++;
}
fclose(in);
return ListInfo;
}
MyMatrix<double> GRIB_ReadFromMessageInfo(GRIB_MessageInfo const& eMesg)
{
grib_handle *h = NULL;
int err;
FILE* in = NULL;
in = fopen(eMesg.FileName.c_str(),"r");
unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;
int idx=0;
while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) {
if (err != GRIB_SUCCESS)
GRIB_CHECK(err,0);
char name_space[3]="ls";
grib_keys_iterator* kiter=NULL;
kiter=grib_keys_iterator_new(h,key_iterator_filter_flags,name_space);
if (!kiter) {
std::cerr << "ERROR: Unable to create keys iterator\n";
throw TerminalException{1};
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
if (nameStr == "shortName")
ShortNameValue=valueStr;
}
grib_keys_iterator_delete(kiter);
if (idx == eMesg.idx) {
long Ni, Nj, numberOfDataPoints;
GRIB_CHECK(grib_get_long(h,"Ni",&Ni),0);
GRIB_CHECK(grib_get_long(h,"Nj",&Nj),0);
GRIB_CHECK(grib_get_long(h,"numberOfDataPoints",&numberOfDataPoints),0);
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
int idxPos=0;
for (int j=0; j<xi_rho; j++)
for (int i=0; i<eta_rho; i++) {
VAL(i,j)=values[idxPos];
idxPos++;
}
free(lats);
free(lons);
free(values);
fclose(in);
return VAL;
}
grib_handle_delete(h);
idx++;
}
fclose(in);
std::cerr << "Failed to find the matching GRIB_MessageInfo\n";
throw TerminalException{1};
}
MyMatrix<double> GRIB_Read2DVariable(std::vector<GRIB_MessageInfo> const& ListMessage, std::string const& VarName)
{
for (auto & eMesg : ListMessage)
if (eMesg.shortName == VarName)
return GRIB_ReadFromMessageInfo(eMesg);
std::cerr << "|ListMessage|=" << ListMessage.size() << "\n";
for (auto & eMesg : ListMessage) {
std::cerr << "Messsage:\n";
std::cerr << "  eMesg.shortName = " << eMesg.shortName << "\n";
std::cerr << "  eMesg.FileName = " << eMesg.FileName << "\n";
std::cerr << "  eMesg.idx = " << eMesg.idx << "\n";
}
std::cerr << "Error in GRIB_Read2DVariable\n";
std::cerr << "Failed to find the variable =" << VarName << "\n";
std::cerr << "Exiting\n";
throw TerminalException{1};
}
MyMatrix<double> GRID_Get2DVariableTimeDifferentiate(TotalArrGetData const& TotalArr, std::string const& eVar, double const& eTimeDay)
{
InterpInfo eInterpInfo=GetTimeDifferentiationInfo(TotalArr.eArr.ListTime, eTimeDay);
double eTimeUpp=TotalArr.eArr.ListTime[eInterpInfo.iTimeUpp];
double eTimeLow=TotalArr.eArr.ListTime[eInterpInfo.iTimeLow];
double DeltaTimeSec = (eTimeUpp - eTimeLow)*double(86400);
int nbTimeStart=TotalArr.eArr.ListStartTime.size();
int nbMesgTotal=TotalArr.eArr.ListAllMessage.size();
auto GetStatusVector=[&](int const& iTimeSearch) -> std::vector<int> {
std::vector<int> ListStatus(nbTimeStart,-1);
for (int iMesg=0; iMesg<nbMesgTotal; iMesg++) {
if (TotalArr.eArr.ListITime[iMesg] == iTimeSearch && TotalArr.eArr.ListAllMessage[iMesg].shortName == eVar) {
int iTimeStart=TotalArr.eArr.ListIStartTime[iMesg];
ListStatus[iTimeStart]=iMesg;
}
}
return ListStatus;
};
std::vector<int> ListStatusLow=GetStatusVector(eInterpInfo.iTimeLow);
std::vector<int> ListStatusUpp=GetStatusVector(eInterpInfo.iTimeUpp);
auto GetITimeStart=[&]() -> int {
for (int iTimeStart=0; iTimeStart<nbTimeStart; iTimeStart++)
if (ListStatusLow[iTimeStart] != -1 && ListStatusUpp[iTimeStart] != -1)
return iTimeStart;
std::cerr << "We failed to find a correct iTimeStart\n";
throw TerminalException{1};
};
int iTimeStart=GetITimeStart();
int iMesgLow=ListStatusLow[iTimeStart];
int iMesgUpp=ListStatusUpp[iTimeStart];
MyMatrix<double> Flow=GRIB_ReadFromMessageInfo(TotalArr.eArr.ListAllMessage[iMesgLow]);
MyMatrix<double> Fupp=GRIB_ReadFromMessageInfo(TotalArr.eArr.ListAllMessage[iMesgUpp]);
int eta_rho=Flow.rows();
int xi_rho=Fupp.cols();
MyMatrix<double> Fret(eta_rho, xi_rho);
for (int iEta=0; iEta<eta_rho; iEta++)
for (int iXi=0; iXi<xi_rho; iXi++) {
double eValLow=Flow(iEta, iXi);
double eValUpp=Fupp(iEta, iXi);
double eDiff=(eValUpp - eValLow) / DeltaTimeSec;
Fret(iEta, iXi) = eDiff;
}
return Fret;
}
MyMatrix<double> GRIB_Get2DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& eVar, double const& eTimeDay)
{
InterpInfo eInterpInfo=GetTimeInterpolationInfo(TotalArr.eArr.ListTime, eTimeDay);
if (eInterpInfo.UseSingleEntry) {
int iTime=eInterpInfo.iTimeLow;
return GRIB_Read2DVariable(TotalArr.eArr.ListListMessages[iTime], eVar);
}
double alphaLow=eInterpInfo.alphaLow;
int iTimeLow=eInterpInfo.iTimeLow;
double alphaUpp=eInterpInfo.alphaUpp;
int iTimeUpp=eInterpInfo.iTimeUpp;
MyMatrix<double> eVarLow=GRIB_Read2DVariable(TotalArr.eArr.ListListMessages[iTimeLow], eVar);
MyMatrix<double> eVarUpp=GRIB_Read2DVariable(TotalArr.eArr.ListListMessages[iTimeUpp], eVar);
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
MyVector<int> GetBoundaryStatus(MyMatrix<int> const& INE, int nbNode)
{
MyVector<int> Status(nbNode);
for (int i=0; i<nbNode; i++)
Status[i]=0;
std::vector<int> PrevVert(nbNode), NextVert(nbNode), Collected(nbNode);
int mne=INE.rows();
for (int ie=0; ie<mne; ie++) {
for (int i=0; i<3; i++) {
int iPrev=PrevIdx(3, i);
int iNext=NextIdx(3, i);
int ip=INE(ie,i);
int ipnext=INE(ie,iNext);
int ipprev=INE(ie,iPrev);
if (Status[ip] == 0) {
Status[ip]=1;
PrevVert[ip]=ipprev;
NextVert[ip]=ipnext;
}
}
}
for (int i=0; i<nbNode; i++)
Status[i]=0;
while(1) {
for (int i=0; i<nbNode; i++)
Collected[i]=0;
for (int ie=0; ie<mne; ie++) {
for (int i=0; i<3; i++) {
int iPrev=PrevIdx(3, i);
int iNext=NextIdx(3, i);
int ip=INE(ie,i);
int ipnext=INE(ie,iNext);
int ipprev=INE(ie,iPrev);
if (Status[ip] == 0) {
int zNext=NextVert[ip];
if (zNext == ipprev) {
Collected[ip]=1;
NextVert[ip]=ipnext;
if (NextVert[ip] == PrevVert[ip])
Status[ip]=1;
}
}
}
}
int IsFinished=1;
for (int i=0; i<nbNode; i++) {
if (Collected[i] == 0 && Status[i] == 0)
Status[i]=-1;
if (Status[i] == 0)
IsFinished=0;
}
if (IsFinished == 1)
break;
}
return Status;
}
MyMatrix<int> GetEdgeSet(MyMatrix<int> const& INE, int nbNode)
{
std::vector<int> ListDegree(nbNode, 0);
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
for (int iNode=0; iNode<nbNode; iNode++)
ListDegree[iNode]=0;
MyMatrix<int> ListAdjacency(nbNode, MaxDeg);
for (int ie=0; ie<mne; ie++) {
int i1=INE(ie,0);
int i2=INE(ie,1);
int i3=INE(ie,2);
int eDeg1=ListDegree[i1];
int eDeg2=ListDegree[i2];
int eDeg3=ListDegree[i3];
ListAdjacency(i1, eDeg1 )=i2;
ListAdjacency(i1, eDeg1+1)=i3;
ListAdjacency(i2, eDeg2 )=i1;
ListAdjacency(i2, eDeg2+1)=i3;
ListAdjacency(i3, eDeg3 )=i1;
ListAdjacency(i3, eDeg3+1)=i2;
ListDegree[i1] = eDeg1 + 2;
ListDegree[i2] = eDeg2 + 2;
ListDegree[i3] = eDeg3 + 2;
}
int nbEdge=0;
for (int iNode=0; iNode<nbNode; iNode++) {
std::set<int> eSet;
int eDeg=ListDegree[iNode];
for (int iAdj=0; iAdj<eDeg; iAdj++) {
int eAdj=ListAdjacency(iNode, iAdj);
if (eAdj > iNode)
eSet.insert(eAdj);
}
nbEdge += eSet.size();
}
MyMatrix<int> ListEdges(nbEdge,2);
int iEdge=0;
for (int iNode=0; iNode<nbNode; iNode++) {
std::set<int> eSet;
int eDeg=ListDegree[iNode];
for (int iAdj=0; iAdj<eDeg; iAdj++) {
int eAdj=ListAdjacency(iNode, iAdj);
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
std::vector<int> GetUnstructuredTriangleAdjInfo_vectint(MyMatrix<int> const& INE, int nbNode)
{
int nbEle=INE.rows();
std::vector<int> ListDegree(nbNode, 0);
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
for (int iNode=0; iNode<nbNode; iNode++)
ListDegree[iNode]=0;
MyMatrix<int> ListAdjacency(nbNode, MaxDeg);
for (int ie=0; ie<mne; ie++) {
for (int i=0; i<3; i++) {
int iNext=NextIdx(3,i);
int iPrev=PrevIdx(3,i);
int i1=INE(ie,i);
int i2=INE(ie,iNext);
int i3=INE(ie,iPrev);
int eDeg=ListDegree[i1];
ListAdjacency(i1, eDeg )=i2;
ListAdjacency(i1, eDeg+1)=i3;
ListDegree[i1] = eDeg + 2;
}
}
int nbEdge=0;
for (int iNode=0; iNode<nbNode; iNode++) {
std::set<int> eSet;
int eDeg=ListDegree[iNode];
for (int iAdj=0; iAdj<eDeg; iAdj++) {
int eAdj=ListAdjacency(iNode, iAdj);
if (eAdj > iNode)
eSet.insert(eAdj);
}
nbEdge += eSet.size();
}
MyMatrix<int> ListEdges(nbEdge,2);
int posEdge=0;
std::vector<int> IndexStart(nbNode);
std::vector<int> IndexEnd (nbNode);
for (int iNode=0; iNode<nbNode; iNode++) {
std::set<int> eSet;
int eDeg=ListDegree[iNode];
for (int iAdj=0; iAdj<eDeg; iAdj++) {
int eAdj=ListAdjacency(iNode, iAdj);
if (eAdj > iNode)
eSet.insert(eAdj);
}
IndexStart[iNode]=posEdge;
for (auto eAdj : eSet) {
ListEdges(posEdge,0)=iNode;
ListEdges(posEdge,1)=eAdj;
posEdge++;
}
IndexEnd [iNode]=posEdge;
}
auto GetIEdge=[&](int const& eVert1, int const& eVert2) -> int {
int idxStart=IndexStart[eVert1];
int idxEnd=IndexEnd[eVert1];
for (int iEdge=idxStart; iEdge<idxEnd; iEdge++) {
if (ListEdges(iEdge,0) != eVert1) {
std::cerr << "Clear inconsistency in code\n";
throw TerminalException{1};
}
if (ListEdges(iEdge,1) == eVert2)
return iEdge;
}
std::cerr << "Failed to find the correct indexes iEdgeIter\n";
throw TerminalException{1};
};
MyMatrix<int> LEdge=GetEdgeSet(INE, nbNode);
std::vector<int> NumberMatch(nbEdge, 0);
MyMatrix<int> IncidenceTrigEdge(nbEdge,2);
for (int iEle=0; iEle<nbEle; iEle++) {
for (int i=0; i<3; i++) {
int iNext=NextIdx(3,i);
int eVert=INE(iEle, i);
int eNext=INE(iEle, iNext);
int eVert1, eVert2;
if (eVert < eNext) {
eVert1=eVert;
eVert2=eNext;
}
else {
eVert1=eNext;
eVert2=eVert;
}
int iEdge=GetIEdge(eVert1, eVert2);
int pos=NumberMatch[iEdge];
IncidenceTrigEdge(iEdge,pos)=iEle;
NumberMatch[iEdge]=pos+1;
}
}
std::vector<int> ListAdj(3*nbEle,-1);
std::vector<int> DegreeTriangle(nbEle,0);
for (int iEdge=0; iEdge<nbEdge; iEdge++)
if (NumberMatch[iEdge] == 2) {
int iEle1=IncidenceTrigEdge(iEdge,0);
int iEle2=IncidenceTrigEdge(iEdge,1);
int pos1=DegreeTriangle[iEle1];
int pos2=DegreeTriangle[iEle2];
ListAdj[3*iEle1 + pos1] = iEle2;
ListAdj[3*iEle2 + pos2] = iEle1;
DegreeTriangle[iEle1]=pos1+1;
DegreeTriangle[iEle2]=pos2+1;
}
return ListAdj;
}
void CHECK_UnstructuredGrid(GridArray const& GrdArr)
{
int mnp=GrdArr.GrdArrRho.LON.rows();
int mne=GrdArr.INE.rows();
std::cerr << "mne=" << mne << "\n";
int nbPlus=0;
int nbMinus=0;
for (int ie=0; ie<mne; ie++) {
int i1=GrdArr.INE(ie,0);
int i2=GrdArr.INE(ie,1);
int i3=GrdArr.INE(ie,2);
if (i1 == i2 || i1 == i3 || i2 == i3) {
std::cerr << "For ie=" << ie << "\n";
std::cerr << "We have i123=" << i1 << "," << i2 << "," << i3 << "\n";
}
double xi = GrdArr.GrdArrRho.LON(i1);
double yi = GrdArr.GrdArrRho.LAT(i1);
double xj = GrdArr.GrdArrRho.LON(i2);
double yj = GrdArr.GrdArrRho.LAT(i2);
double xk = GrdArr.GrdArrRho.LON(i3);
double yk = GrdArr.GrdArrRho.LAT(i3);
double area=xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj);
if (area > 0)
nbPlus++;
if (area < 0)
nbMinus++;
for (int i=0; i<3; i++) {
int IP=GrdArr.INE(ie,i);
if (IP < 0 || IP >= mnp) {
std::cerr << "Error in the unstructured grid\n";
std::cerr << "mnp=" << mnp << "  mne=" << mne << "\n";
std::cerr << "ie=" << ie << "\n";
std::cerr << "INE=[" << GrdArr.INE(ie,0) << " , " << GrdArr.INE(ie,1) << " , " << GrdArr.INE(ie,2) << "]\n";
throw TerminalException{1};
}
}
}
if (nbPlus > 0 && nbMinus > 0) {
std::cerr << "The grid is incorrectly oriented\n";
std::cerr << "mne=" << mne << " : nbPlus=" << nbPlus << "  nbMinus=" << nbMinus << "\n";
throw TerminalException{1};
}
MyVector<int> Status=GetBoundaryStatus(GrdArr.INE, mnp);
int nbStatusNormal=0;
int nbStatusBound=0;
for (int i=0; i<mnp; i++) {
if (Status[i] == 1)
nbStatusNormal++;
if (Status[i] == -1)
nbStatusBound++;
}
std::cerr << "nbStatusNormal=" << nbStatusNormal << " nbStatusBound=" << nbStatusBound << "\n";
}
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
std::vector<int> CHILF(mnp,0);
Eigen::Tensor<int,3> CELLVERTEX(mnp,MAXMNECON,2);
for (int ie=0; ie<mne; ie++)
for (int j=0; j<3; j++) {
int i=GrdArr.INE(ie,j);
CELLVERTEX(i, CHILF[i] ,0) = ie;
CELLVERTEX(i, CHILF[i] ,1) = j;
CHILF[i]++;
}
int COUNT_MAX=0;
for (int ip=0; ip<mnp; ip++)
COUNT_MAX += CCON[ip];
MyMatrix<int> IE_CELL2 (mnp,MAXMNECON);
MyMatrix<int> POS_CELL2(mnp,MAXMNECON);
int j=0;
for (int ip=0; ip<mnp; ip++)
for (int i=0; i<CCON[ip]; i++) {
IE_CELL2 (ip,i) = CELLVERTEX(ip,i,0);
POS_CELL2(ip,i) = CELLVERTEX(ip,i,1);
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
int ie2=IE_CELL2(ip,icon);
if (ie != ie2) {
int POS=POS_CELL2(ip,icon);
int POS_NEXT=POS_TRICK[POS][0];
int IP_ADJ_NEXT=GrdArr.INE(ie2, POS_NEXT);
if (IP_ADJ_NEXT == IP_NEXT) {
std::cerr << "Combinatorial orientability problem\n";
std::cerr << "IE=" << ie << " IE2=" << ie2 << "\n";
std::cerr << "IP=" << ip << " IP_NEXT=" << IP_NEXT << "\n";
throw TerminalException{1};
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
throw TerminalException{1};
}
}
std::cerr << "Now leaving the combinatorial check\n";
}
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
throw TerminalException{1};
}
std::cerr << "nbPlus = " << nbPlus << "  nbMinus = " << nbMinus << "\n";
}
std::string NCL_bool(bool eBool)
{
if (eBool)
return "True";
return "False";
}
struct AnnotationRec {
bool DrawAnnotation = false;
double AnnotationLon;
double AnnotationLat;
std::string AnnotationText;
};
struct SeqLineSegment {
std::vector<PairLL> ListPairLL;
bool IsClosed;
};
struct DrawArr {
bool DoTitle;
std::string TitleStr;
std::string VarNameUF;
bool DrawRiver;
bool DrawContourBathy;
bool PrintMMA;
bool DoColorBar;
std::string ColorMap;
std::string cnFillMode;
bool cnSmoothingOn;
int nbLevelSpa;
int nbLabelStride;
QuadArray eQuadFrame;
AnnotationRec TheAnnot;
bool FillLand;
bool UseNativeGrid;
std::string GridResolution;
std::vector<SeqLineSegment> ListLineSegment;
};
struct QuadDrawInfo {
std::string eFrameName;
int iFrame;
QuadArray eQuad;
};
struct InterpolToUVpoints {
GridArray GrdArr;
MySparseMatrix<double> InterpMat;
};
void CALL_NCL(bool const& KeepNC_NCL,
std::string const& TargetFile,
std::string const& eFileNC,
std::string const& eFileNCL)
{
std::string eComm="ncl";
std::string eCommand=eComm + " " + eFileNCL + " > /dev/null";
if (KeepNC_NCL == true) {
std::cerr << "eCommand = " << eCommand << "\n";
}
int iret=system(eCommand.c_str());
if (iret == -1) {
printf("Oh dear, something went wrong with ncl! %s\n", strerror(errno));
throw TerminalException{1};
}
if (KeepNC_NCL == true) {
std::cerr << "eFileNC  = " << eFileNC << "\n";
std::cerr << "eFileNCL = " << eFileNCL << "\n";
}
if (IsExistingFile(TargetFile) == false) {
std::cerr << "The following TargetFile was not created\n";
std::cerr << "TargetFile = " << TargetFile << "\n";
std::cerr << "eFileNC    = " << eFileNC << "\n";
std::cerr << "eFileNCL   = " << eFileNCL << "\n";
std::cerr << "Please debug\n";
throw TerminalException{1};
}
if (KeepNC_NCL == false) {
RemoveFileIfExist(eFileNC);
RemoveFileIfExist(eFileNCL);
}
}
struct TripleNCL {
std::string TargetFile;
std::string eFileNC;
std::string eFileNCL;
};
struct NCLcaller {
NCLcaller() = delete;
NCLcaller(bool const& KeepNC_NCL, int const& nproc) : eKeep(KeepNC_NCL), NPROC(nproc)
{
if (eKeep && nproc > 1) {
std::cerr << "Cannot have KeepNC_NCL = T and NPROC > 1\n";
throw TerminalException{1};
}
ListExch=new int[NPROC];
ListTerm=new int[NPROC];
ListTripleNCL=new TripleNCL[NPROC];
ListCond=std::vector<std::condition_variable>(nproc);
NbRunningJob=0;
InWhile=0;
for (int iProc=0; iProc<NPROC; iProc++) {
ListExch[iProc]=0;
ListTerm[iProc]=0;
}
auto IterationLoop=[&](int iproc, int *Exch, int* Term, TripleNCL *eCall) {
int IsFirst=1;
while(1) {
if (IsFirst == 1) {
InWhile++;
IsFirst=0;
}
if (*Exch == 1) {
CALL_NCL(eKeep, eCall->TargetFile, eCall->eFileNC, eCall->eFileNCL);
*Exch=0;
NbRunningJob--;
sub_cv.notify_one();
}
if (*Term == -1) {
break;
}
if (*Exch == 0) {
std::unique_lock<std::mutex> lk(inst_mtx);
ListCond[iproc].wait(lk, [&]{return *Exch == 1 || *Term == -1;});
}
}
InWhile--;
fin_cv.notify_one();
};
for (int iProc=0; iProc<NPROC; iProc++) {
ListThr.push_back(std::thread(IterationLoop, iProc, &(ListExch[iProc]), &(ListTerm[iProc]), &(ListTripleNCL[iProc])));
}
for (int iProc=0; iProc<NPROC; iProc++)
ListThr[iProc].detach();
}
~NCLcaller()
{
for (int iProc=0; iProc<NPROC; iProc++) {
ListTerm[iProc]=-1;
ListCond[iProc].notify_one();
}
std::unique_lock<std::mutex> lk(fin_mtx);
fin_cv.wait(lk, [&]{return InWhile == 0;});
delete [] ListExch;
delete [] ListTerm;
delete [] ListTripleNCL;
}
void SubmitJob(std::string const& TargetFile, std::string const& eFileNC, std::string const& eFileNCL)
{
std::unique_lock<std::mutex> lk(sub_mtx);
sub_cv.wait(lk, [&]{return NbRunningJob < NPROC;});
int iProcFound=-1;
for (int iProc=0; iProc<NPROC; iProc++)
if (ListExch[iProc] == 0)
iProcFound=iProc;
if (iProcFound == -1) {
std::cerr << "Failed to find the processor.\n";
std::cerr << "Bug is in NCLcaller\n";
std::cerr << "NbRunningJob=" << NbRunningJob << " NPROC=" << NPROC << "\n";
for (int iProc=0; iProc<NPROC; iProc++)
std::cerr << "iProc=" << iProc << " Exch=" << ListExch[iProc] << " Term=" << ListTerm[iProc] << "\n";
throw TerminalException{1};
}
ListTripleNCL[iProcFound].TargetFile=TargetFile;
ListTripleNCL[iProcFound].eFileNC=eFileNC;
ListTripleNCL[iProcFound].eFileNCL=eFileNCL;
ListExch[iProcFound]=1;
NbRunningJob++;
ListCond[iProcFound].notify_one();
}
private:
std::mutex fin_mtx;
std::condition_variable fin_cv;
std::atomic<int> NbRunningJob;
std::atomic<int> InWhile;
std::mutex sub_mtx;
std::condition_variable sub_cv;
std::mutex inst_mtx;
bool eKeep;
int NPROC;
std::vector<std::thread> ListThr;
int* ListExch;
int* ListTerm;
TripleNCL* ListTripleNCL;
std::vector<std::condition_variable> ListCond;
};
struct coor {
double x;
double y;
};
coor operator+(coor const& c1, coor const& c2)
{
return {c1.x + c2.x, c1.y + c2.y};
}
coor MultScal(coor const& c, double const& eScal)
{
return {eScal*c.x, eScal*c.y};
}
struct SVGqualInfo {
std::vector<int> color;
int Size;
std::string MarkerEnd;
};
struct SVGbezier {
coor pointM;
coor pointC;
coor point1;
coor point2;
SVGqualInfo eQual;
};
struct SVGpolyline {
std::vector<coor> ListCoor;
SVGqualInfo eQual;
};
struct SVGline {
coor ePt;
coor fPt;
SVGqualInfo eQual;
};
struct SVGplotDescription {
std::vector<SVGpolyline> ListPolyline;
std::vector<SVGbezier> ListBezier;
std::vector<SVGline> ListLine;
double height;
double width;
double scale_factor;
double add_offsetX;
double add_offsetY;
int FrameOption;
int RoundMethod;
};
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
struct FullNamelist {
std::map<std::string,SingleBlock> ListBlock;
std::string FileName;
};
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
std::string NAMELIST_RemoveAfterCommentChar(std::string const&eStr, std::string &eChar)
{
bool WeFound=false;
std::string RetStr;
int len=eStr.size();
for (int i=0; i<len; i++) {
std::string fChar=eStr.substr(i,1);
if (fChar == eChar)
WeFound=true;
if (WeFound == false)
RetStr += eStr.at(i);
}
return RetStr;
}
std::string NAMELIST_RemoveAfterLastChar(std::string const& eStr, std::string const& eLastChar)
{
int iPos=-1;
int len=eStr.size();
for (int i=0; i<len; i++) {
int j=len-1-i;
if (iPos == -1) {
std::string eChar=eStr.substr(j,1);
if (eChar == eLastChar)
iPos=j;
}
}
if (iPos == -1)
return eStr;
std::string RetStr;
for (int i=0; i<iPos; i++)
RetStr=RetStr + eStr.at(i);
return RetStr;
}
std::string NAMELIST_ClearEndOfLine(std::string const& eStr)
{
std::string eCharCommentB="!";
std::string eStr3=NAMELIST_RemoveAfterLastChar(eStr, eCharCommentB);
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
if (test)
return eStrPrior;
return eStr3;
}
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
throw TerminalException{1};
}
std::string NAMELIST_ConvertFortranStringToCppString(std::string const& eStr)
{
int len=eStr.length();
std::string eFirstChar=eStr.substr(0, 1);
std::string eLastChar=eStr.substr(len-1, 1);
int RemovableEnding=0;
if (eFirstChar == "'" || eFirstChar == "\"") {
RemovableEnding=1;
if (eFirstChar != eLastChar) {
std::cerr << "eFirstChar = " << eFirstChar << "\n";
std::cerr << " eLastChar = " << eLastChar << "\n";
std::cerr << "The character used for noting beginning and end of string should be identical\n";
throw TerminalException{1};
}
}
if (RemovableEnding == 1)
return eStr.substr(1,len-2);
return eStr;
}
std::vector<std::string> NAMELIST_ConvertFortranListStringToCppListString(std::string const& eStr)
{
int len=eStr.length();
std::string eFirstChar=eStr.substr(0, 1);
std::string eLastChar=eStr.substr(len-1,1);
if (eFirstChar != "'" && eFirstChar != "\"") {
std::cerr << "eStr=" << eStr << "\n";
std::cerr << "For list of strings, one should use string \"  \"   or '    '   \n";
throw TerminalException{1};
}
if (eLastChar != "'" && eLastChar != "\"") {
std::cerr << "eStr=" << eStr << "\n";
std::cerr << "For list of strings, one should use string \"  \"   or '    '   \n";
throw TerminalException{1};
}
if (eFirstChar != eLastChar) {
std::cerr << "eStr=" << eStr << "\n";
std::cerr << "eFirstChar=" << eFirstChar << "\n";
std::cerr << "eLastChar=" << eLastChar << "\n";
std::cerr << "No coherency in endings\n";
throw TerminalException{1};
}
std::string eSepChar=eFirstChar;
int IsInString=0;
std::string eFound="";
std::vector<std::string> eListStr;
for (int i=0; i<len; i++) {
std::string eChar=eStr.substr(i,1);
if (eChar == eSepChar) {
eFound += eChar;
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
eFound += eChar;
}
}
return eListStr;
}
std::vector<double> NAMELIST_ConvertFortranStringListDoubleToCppVectorDouble(std::string const& eVarValue)
{
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
return eListRetDouble;
}
std::vector<int> NAMELIST_ConvertFortranStringListIntToCppVectorInt(std::string const& eVarValue)
{
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
return eListRetInt;
}
void NAMELIST_WriteBlock(std::ostream &os, std::string const& eBlockName, SingleBlock const& eBlock)
{
os << "&" << eBlockName << "\n";
for (std::map<std::string,int>::const_iterator it=eBlock.ListIntValues.begin(); it!=eBlock.ListIntValues.end(); ++it)
os << "  " << it->first << " = " << it->second << "\n";
for (std::map<std::string,bool>::const_iterator it=eBlock.ListBoolValues.begin(); it!=eBlock.ListBoolValues.end(); ++it) {
bool eVal=it->second;
std::string eValStr;
if (eVal == false) {
eValStr="F";
}
else {
eValStr="T";
}
os << "  " << it->first << " = " << eValStr << "\n";
}
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
os << "  " << it->first << " = \"" << it->second << "\"\n";
for (std::map<std::string,std::vector<std::string> >::const_iterator it=eBlock.ListListStringValues.begin(); it!=eBlock.ListListStringValues.end(); ++it) {
os << "  " << it->first << " = ";
std::vector<std::string> eListStr=it->second;
int nbString=eListStr.size();
for (int iString=0; iString<nbString; iString++) {
if (iString > 0)
os << ",";
os << "\"" << eListStr[iString] << "\"";
}
os << "\n";
}
os << "&END\n";
}
void NAMELIST_WriteNamelistFile(std::ostream &os, FullNamelist const& eFullNamelist)
{
int iBlock=0;
for (std::map<std::string,SingleBlock>::const_iterator itB=eFullNamelist.ListBlock.begin(); itB!=eFullNamelist.ListBlock.end(); ++itB) {
std::string eBlockName=itB->first;
SingleBlock eBlock=itB->second;
if (iBlock > 0)
os << "\n\n";
NAMELIST_WriteBlock(os, eBlockName, eBlock);
iBlock++;
}
}
void NAMELIST_ReadNamelistFile(std::string const& eFileName, FullNamelist &eFullNamelist)
{
std::set<std::pair<std::string, std::string>> ListInsertValues;
if (IsExistingFile(eFileName) == false) {
std::cerr << "The following namelist file is missing\n";
std::cerr << "eFileName = " << eFileName << "\n";
throw TerminalException{1};
}
std::ifstream INfs(eFileName);
bool InBlock=false;
std::string eBlockName;
while(!INfs.eof()) {
std::string Ampersand="&";
std::string strTab="\t";
int siz=1024;
char eLine[siz];
INfs.getline(eLine, siz);
std::string PreStr=eLine;
std::string eCharComment="!";
std::string PreStrB=NAMELIST_RemoveAfterCommentChar(PreStr, eCharComment);
std::string eStr=STRING_RemoveSpacesBeginningEnd(PreStrB);
int len=eStr.length();
if (eStr.find(strTab) != std::string::npos) {
std::cerr << "Tabs are not allowed\n";
std::cerr << "LINE=" << eStr << "\n";
throw TerminalException{1};
}
if (len> 0) {
if (eStr.find(Ampersand) != std::string::npos) {
std::string eFirstChar=eStr.substr(0,1);
if (eFirstChar != "&") {
std::cerr << "Error while reading namelist file = " << eFileName << "\n";
std::cerr << "Error, Ampersand (&) should be only in the first character\n";
std::cerr << "LINE=" << eStr << "\n";
throw TerminalException{1};
}
std::string strRed=eStr.substr(1, len-1);
if (InBlock == false) {
eBlockName=strRed;
auto search=eFullNamelist.ListBlock.find(eBlockName);
if (search == eFullNamelist.ListBlock.end() ) {
std::cerr << "Find BlockName = " << eBlockName << "\n";
std::cerr << "which is not in the authorized list\n";
std::cerr << "LINE=" << eStr << "\n";
std::cerr << "List of authorized block names:\n";
for (auto & eBlock : eFullNamelist.ListBlock)
std::cerr << "Block name=" << eBlock.first << "\n";
throw TerminalException{1};
}
InBlock=true;
}
else {
if (strRed != "END") {
std::cerr << "Ampersand detected. We should leave with a END\n";
std::cerr << "LINE=" << eStr << "\n";
throw TerminalException{1};
}
InBlock=false;
}
}
else {
if (eStr != "/") {
std::string eStr3=NAMELIST_ClearEndOfLine(eStr);
std::string strEqual="=";
int posEqual=STRING_GetCharPositionInString(eStr3, strEqual);
if (posEqual != -1) {
int len3=eStr3.length();
std::string eStrPrior=eStr3.substr(0, posEqual);
std::string eStrPosterior=eStr3.substr(posEqual+1, len3-posEqual-1);
std::string eVarName=STRING_RemoveSpacesBeginningEnd(eStrPrior);
std::string eVarValue=STRING_RemoveSpacesBeginningEnd(eStrPosterior);
std::string eVarNature=NAMELIST_FindPositionVariableInBlock(
eVarName, eFullNamelist.ListBlock[eBlockName]);
std::pair<std::string, std::string> ePair{eBlockName, eVarName};
auto searchB=ListInsertValues.find(ePair);
if (searchB != ListInsertValues.end()) {
std::cerr << "In the block " << eBlockName << "\n";
std::cerr << "the entry " << eVarName << "\n";
std::cerr << "is defined two times\n";
throw TerminalException{1};
}
ListInsertValues.insert(ePair);
if (eVarNature == "not found") {
NAMELIST_WriteBlock(std::cerr, eBlockName, eFullNamelist.ListBlock[eBlockName]);
std::cerr << "Error in reading the NAMELIST file. See above allowed entries\n";
std::cerr << "The variable " << eVarName << "\n";
std::cerr << "is in block " << eBlockName << "\n";
std::cerr << "of the file " << eFileName << "\n";
std::cerr << "but it is not allowed for the chosen application\n";
throw TerminalException{1};
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
else {
int nbChar=eStr3.size();
if (nbChar != 0) {
std::cerr << "If lines has no = sign then it should be empty\n";
throw TerminalException{1};
}
}
}
else {
InBlock=false;
}
}
}
}
if (InBlock == true) {
std::cerr << "Error. When leaving namelist reading, we should be out of block\n";
throw TerminalException{1};
}
}
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
throw TerminalException{1};
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
void COMPUTE_STATUS_LON_LAT(std::vector<int> &ListStatus, MyVector<double> const& ListLon, MyVector<double> const& ListLat, FullNamelist const& eFull)
{
int eStatus;
std::map<std::string, SingleBlock> ListBlock=eFull.ListBlock;
int nbMeas=ListStatus.size();
int GEOSELECTION=ListBlock.at("SELECT").ListIntValues.at("GEOSELECTION");
double MinLON=ListBlock.at("SELECT").ListDoubleValues.at("MinLON");
double MaxLON=ListBlock.at("SELECT").ListDoubleValues.at("MaxLON");
double MinLAT=ListBlock.at("SELECT").ListDoubleValues.at("MinLAT");
double MaxLAT=ListBlock.at("SELECT").ListDoubleValues.at("MaxLAT");
std::vector<double> LonPoly=ListBlock.at("SELECT").ListListDoubleValues.at("LONPOLY");
std::vector<double> LatPoly=ListBlock.at("SELECT").ListListDoubleValues.at("LATPOLY");
for (int iMeas=0; iMeas<nbMeas; iMeas++)
if (ListStatus[iMeas] == 1) {
double eLon=ListLon(iMeas);
double eLat=ListLat(iMeas);
if (GEOSELECTION == 1) {
if (eLon >= MinLON && eLon <= MaxLON && eLat >= MinLAT && eLat <= MaxLAT) {
eStatus=1;
}
else {
eStatus=0;
}
ListStatus[iMeas]=eStatus;
}
if (GEOSELECTION == 2) {
bool IsInside=IsPointInside(eLon, eLat, LonPoly, LatPoly);
eStatus=0;
if (IsInside)
eStatus=1;
ListStatus[iMeas]=eStatus;
}
}
}
void COMPUTE_STATUS_TIME(std::vector<int> &ListStatus, std::vector<double> const& ListTime, FullNamelist const& eFull)
{
std::map<std::string, SingleBlock> ListBlock=eFull.ListBlock;
SingleBlock eBlSELECT=ListBlock.at("SELECT");
double BeginTime=CT2MJD(eBlSELECT.ListStringValues.at("BEGTC"));
double EndTime =CT2MJD(eBlSELECT.ListStringValues.at("ENDTC"));
int nbMeas=ListStatus.size();
for (int iMeas=0; iMeas<nbMeas; iMeas++)
if (ListStatus[iMeas] == 1) {
double eTime=ListTime[iMeas];
int eStatus=0;
if (eTime >= BeginTime && eTime <= EndTime) {
eStatus=1;
}
ListStatus[iMeas]=eStatus;
}
}
void CheckSatelliteList(FullNamelist const& eFull)
{
std::vector<std::string> AllSatNames=GetAllNamesOfSatelliteAltimeter();
std::vector<std::string> AllowedSatNames=eFull.ListBlock.at("SELECT").ListListStringValues.at("AllowedSatellites");
for (auto& eSatName : AllowedSatNames) {
std::cerr << "eSatName=" << eSatName << "\n";
bool WeMatch=false;
for (auto& fSatName : AllSatNames)
if (fSatName == eSatName)
WeMatch=true;
if (WeMatch == false) {
std::cerr << "Error in input.\n";
std::cerr << "AllowedSatellites was assigned to:\n";
WriteStdVector(std::cerr, AllowedSatNames);
std::cerr << "On the other the only satellites that you can put are following:\n";
WriteStdVector(std::cerr, AllSatNames);
throw TerminalException{1};
}
}
}
std::vector<int> GetListSatelliteId_vect(FullNamelist const& eFull)
{
std::set<int> SatelliteId;
CheckSatelliteList(eFull);
std::vector<std::string> AllSatNames=GetAllNamesOfSatelliteAltimeter();
int nbSatellite=AllSatNames.size();
std::vector<int> ListStatusSatellite(nbSatellite,0);
std::vector<std::string> AllowedSatNames=eFull.ListBlock.at("SELECT").ListListStringValues.at("AllowedSatellites");
for (auto& eSatName : AllowedSatNames)
for (int iSat=0; iSat<nbSatellite; iSat++)
if (AllSatNames[iSat] == eSatName)
ListStatusSatellite[iSat]=1;
return ListStatusSatellite;
}
struct PairListWindWave {
int eSat;
std::vector<PairMM> ListPairWind;
std::vector<PairMM> ListPairWave;
};
struct SatelliteListTrack {
int eSat;
std::vector<std::vector<SingleEntryAltimeter> > ListListEntAltimeter;
double avgDistKM;
};
std::vector<int> GetListStatusTrackLength(std::vector<SingleEntryAltimeter> const& eVectEnt, FullNamelist const& eFull)
{
int siz=eVectEnt.size();
std::vector<int> ListStatus(siz,1);
double MaxDistTrackPointKM=eFull.ListBlock.at("SELECT").ListDoubleValues.at("MaxDistTrackPointKM");
int MinimalTrackSize=eFull.ListBlock.at("SELECT").ListIntValues.at("MinimalTrackSize");
std::set<int> SatelliteId;
for (auto& eEnt : eVectEnt)
SatelliteId.insert(eEnt.Satellite);
for (auto& eSat : SatelliteId) {
std::vector<SingleEntryAltimeter> ListEnt;
std::vector<int> ListIL;
int nbL=eVectEnt.size();
for (int iL=0; iL<nbL; iL++) {
SingleEntryAltimeter eEnt=eVectEnt[iL];
if (eEnt.Satellite == eSat) {
ListEnt.push_back(eEnt);
ListIL.push_back(iL);
}
}
int nbEnt=ListEnt.size();
for (int iEnt=0; iEnt<nbEnt-1; iEnt++) {
double eDiff=ListEnt[iEnt+1].Time - ListEnt[iEnt].Time;
if (eDiff < 0) {
std::cerr << "We have a decrease in time when it should increase\n";
throw TerminalException{1};
}
}
std::vector<int> ListSep;
for (int iEnt=0; iEnt<nbEnt-1; iEnt++) {
double eLon1=ListEnt[iEnt].Lon;
double eLat1=ListEnt[iEnt].Lat;
double eLon2=ListEnt[iEnt+1].Lon;
double eLat2=ListEnt[iEnt+1].Lat;
double eDistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
if (eDistKM >= MaxDistTrackPointKM)
ListSep.push_back(iEnt);
}
int nbScene=ListSep.size() + 1;
for (int iScene=0; iScene<nbScene; iScene++) {
int iBegin, iEnd;
if (nbScene == 0) {
iBegin=0;
iEnd=nbEnt-1;
}
else {
if (iScene == 0) {
iBegin=0;
iEnd=ListSep[0];
}
else {
if (iScene == nbScene-1) {
iBegin=ListSep[nbScene-2]+1;
iEnd=nbEnt-1;
}
else {
iBegin=ListSep[iScene-1]+1;
iEnd=ListSep[iScene];
}
}
}
int len=1 + iEnd - iBegin;
if (len < MinimalTrackSize) {
for (int i=iBegin; i<=iEnd; i++) {
int iL=ListIL[i];
ListStatus[iL]=0;
}
}
}
}
return ListStatus;
}
std::vector<SatelliteListTrack> GetListTrackAltimeter(std::vector<SingleEntryAltimeter> const& eVectEnt, double const& avgDistKM_target, FullNamelist const& eFull)
{
std::vector<SatelliteListTrack> RetList;
double MaxDistTrackPointKM=eFull.ListBlock.at("SELECT").ListDoubleValues.at("MaxDistTrackPointKM");
bool DoTrackSmoothing=eFull.ListBlock.at("SELECT").ListBoolValues.at("DoTrackSmoothing");
std::set<int> SatelliteId;
for (auto& eEnt : eVectEnt)
SatelliteId.insert(eEnt.Satellite);
for (auto& eSat : SatelliteId) {
std::vector<SingleEntryAltimeter> ListEnt;
for (auto & eEnt : eVectEnt)
if (eEnt.Satellite == eSat)
ListEnt.push_back(eEnt);
int nbEnt=ListEnt.size();
for (int iEnt=0; iEnt<nbEnt-1; iEnt++) {
double eDiff=ListEnt[iEnt+1].Time - ListEnt[iEnt].Time;
if (eDiff < 0) {
std::cerr << "We have a decrease in time when it should increase\n";
throw TerminalException{1};
}
}
std::vector<int> ListSep;
double SumDistKM=0;
int nbPair=0;
for (int iEnt=0; iEnt<nbEnt-1; iEnt++) {
double eLon1=ListEnt[iEnt].Lon;
double eLat1=ListEnt[iEnt].Lat;
double eLon2=ListEnt[iEnt+1].Lon;
double eLat2=ListEnt[iEnt+1].Lat;
double eDistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
if (eDistKM < MaxDistTrackPointKM) {
SumDistKM += eDistKM;
nbPair++;
}
else {
ListSep.push_back(iEnt);
}
}
double avgDistKM_track=SumDistKM/double(nbPair);
SmoothArr eSmoothArr=GetSmoothingArray(avgDistKM_target, avgDistKM_track);
int nbScene=ListSep.size() + 1;
std::vector<std::vector<SingleEntryAltimeter> > ListListEntAltimeter;
for (int iScene=0; iScene<nbScene; iScene++) {
int iBegin, iEnd;
if (nbScene == 0) {
iBegin=0;
iEnd=nbEnt-1;
}
else {
if (iScene == 0) {
iBegin=0;
iEnd=ListSep[0];
}
else {
if (iScene == nbScene-1) {
iBegin=ListSep[nbScene-2]+1;
iEnd=nbEnt-1;
}
else {
iBegin=ListSep[iScene-1]+1;
iEnd=ListSep[iScene];
}
}
}
int len=1 + iEnd - iBegin;
std::vector<SingleEntryAltimeter> ListEntAltimeter(len);
for (int i=iBegin; i<=iEnd; i++)
ListEntAltimeter[i-iBegin]=ListEnt[i];
if (DoTrackSmoothing) {
ListListEntAltimeter.push_back(SpatialAveragingTrack(ListEntAltimeter, eSmoothArr));
}
else {
ListListEntAltimeter.push_back(ListEntAltimeter);
}
}
RetList.push_back({eSat, ListListEntAltimeter, avgDistKM_track});
}
return RetList;
}
std::set<int> GetListSatelliteId_set(std::vector<SingleEntryAltimeter> const& eVectEnt, FullNamelist const& eFull)
{
std::set<int> PreSatelliteId;
for (auto& eEnt : eVectEnt)
PreSatelliteId.insert(eEnt.Satellite);
std::set<int> SatelliteId;
CheckSatelliteList(eFull);
std::vector<std::string> AllSatNames=GetAllNamesOfSatelliteAltimeter();
std::vector<std::string> AllowedSatNames=eFull.ListBlock.at("SELECT").ListListStringValues.at("AllowedSatellites");
for (auto& eEnt : PreSatelliteId) {
std::string eSatName=AllSatNames[eEnt-1];
bool WeMatch=false;
for (auto& fSatName : AllowedSatNames)
if (eSatName == fSatName)
WeMatch=true;
if (WeMatch)
SatelliteId.insert(eEnt);
}
return SatelliteId;
}
std::vector<PairListWindWave> MergeTracksForRawStatistics(std::vector<SatelliteListTrack> const& LTrack, FullNamelist const& eFull)
{
SingleBlock eBlSEL=eFull.ListBlock.at("SELECT");
double MinWind=eBlSEL.ListDoubleValues.at("MinWIND");
double MaxWind=eBlSEL.ListDoubleValues.at("MaxWIND");
double MinHs=eBlSEL.ListDoubleValues.at("MinHS");
double MaxHs=eBlSEL.ListDoubleValues.at("MaxHS");
std::vector<PairListWindWave> RetList;
for (auto& eRecTrack : LTrack) {
std::vector<PairMM> ListPairWind;
std::vector<PairMM> ListPairWave;
int nbMatch=0;
int nbCorrWind=0;
int nbCorrWave=0;
for (auto & eListEnt : eRecTrack.ListListEntAltimeter)
for (auto & eEnt : eListEnt) {
nbMatch++;
if (eEnt.Swh_used < MaxHs && eEnt.Swh_used > MinHs) {
nbCorrWave++;
PairMM ePairHs{eEnt.Swh_used, eEnt.Swh_model};
ListPairWave.push_back(ePairHs);
}
if (eEnt.WindSpeed_used < MaxWind && eEnt.WindSpeed_used > MinWind) {
nbCorrWind++;
PairMM ePairWind{eEnt.WindSpeed_used, eEnt.WindSpeed_model};
ListPairWind.push_back(ePairWind);
}
}
std::cerr << "nbMatch=" << nbMatch << "  nbCorr(wind/wave)=" << nbCorrWind << "/" << nbCorrWave << "\n";
PairListWindWave eSS{eRecTrack.eSat, ListPairWind, ListPairWave};
RetList.push_back(eSS);
}
return RetList;
}
std::vector<PairLL> ReadLonLatDiscFile(std::string const& eFile)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "Error in ReadLonLatDiscFile\n";
std::cerr << "The file eFile=" << eFile << "\n";
std::cerr << "is missing\n";
throw TerminalException{1};
}
int nbLine=FILE_GetNumberLine(eFile);
std::ifstream is;
is.open(eFile);
std::vector<PairLL> ListPt(nbLine);
for (int iLine=0; iLine<nbLine; iLine++) {
double eLon, eLat;
is >> eLon;
is >> eLat;
PairLL ePt{eLon, eLat};
ListPt[iLine]=ePt;
}
return ListPt;
}
bool IsPointInside_Point(PairLL const& ePt, std::vector<PairLL> const& ListPt)
{
int nvert=ListPt.size();
int i, j;
bool c=false;
for (i = 0, j = nvert-1; i < nvert; j = i++) {
if ( ((ListPt[i].eLat > ePt.eLat) != (ListPt[j].eLat > ePt.eLat)) &&
(ePt.eLon < (ListPt[j].eLon-ListPt[i].eLon) * (ePt.eLat-ListPt[i].eLat) / (ListPt[j].eLat-ListPt[i].eLat ) + ListPt[i].eLon) )
c = !c;
}
return c;
}
std::vector<PairLL> GetGridBoundary(MyMatrix<double> const& LON, MyMatrix<double> const& LAT, int const& iStart, int const& iEnd, int const& jStart, int const& jEnd)
{
int len=2*(iEnd - iStart) + 2*(jEnd - jStart);
std::vector<PairLL> ListPt(len);
int idx=0;
for (int iEta=iStart; iEta<=iEnd-1; iEta++) {
int iXi=jStart;
PairLL ePt{LON(iEta, iXi), LAT(iEta, iXi)};
ListPt[idx]=ePt;
idx++;
}
for (int iXi=jStart; iXi<=jEnd-1; iXi++) {
int iEta=iEnd;
PairLL ePt{LON(iEta, iXi), LAT(iEta, iXi)};
ListPt[idx]=ePt;
idx++;
}
for (int iEta=iEnd; iEta>=iStart+1; iEta--) {
int iXi=jEnd;
PairLL ePt{LON(iEta, iXi), LAT(iEta, iXi)};
ListPt[idx]=ePt;
idx++;
}
for (int iXi=jEnd; iXi>=jStart+1; iXi--) {
int iEta=iStart;
PairLL ePt{LON(iEta, iXi), LAT(iEta, iXi)};
ListPt[idx]=ePt;
idx++;
}
return ListPt;
}
PairCoord FindContaining(PairLL const& ePt, MyMatrix<double> const& LON, MyMatrix<double> const& LAT)
{
int eta_rho=LON.rows();
int xi_rho=LON.cols();
int iStart, iEnd, jStart, jEnd;
iStart=0;
iEnd=eta_rho-1;
jStart=0;
jEnd=xi_rho-1;
while(1) {
std::vector<PairLL> ListPt1=GetGridBoundary(LON, LAT, iStart, iEnd, jStart, jEnd);
bool test1=IsPointInside_Point(ePt, ListPt1);
if (test1 == false)
return {-1, -1};
int iDiff=iEnd - iStart;
int jDiff=jEnd - jStart;
if (iDiff == 1 && jDiff == 1)
break;
if (iDiff > 1) {
int iMid=int(roundl((float(iStart) + float(iEnd))/double(2)));
std::vector<PairLL> ListPt2=GetGridBoundary(LON, LAT, iStart, iMid, jStart, jEnd);
bool test2=IsPointInside_Point(ePt, ListPt2);
if (test2 == true) {
iEnd=iMid;
}
else {
iStart=iMid;
}
}
if (jDiff > 1) {
int jMid=int(roundl((float(jStart) + float(jEnd))/double(2)));
std::vector<PairLL> ListPt3=GetGridBoundary(LON, LAT, iStart, iEnd, jStart, jMid);
bool test3=IsPointInside_Point(ePt, ListPt3);
if (test3 == true) {
jEnd=jMid;
}
else {
jStart=jMid;
}
}
}
return {iStart, jStart};
}
struct DiscInfo {
PairLL eSample;
PairLL avgPoint;
double SpreadLon;
double SpreadLat;
double MaxSpread;
};
struct KTreeElt {
std::vector<PairLL> ListPt;
DiscInfo eDisc;
bool IsSplit;
int iSub1;
int iSub2;
};
DiscInfo KTree_ComputeDisc(std::vector<PairLL> const& ListPt)
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
PairLL eSample=ListPt[idxMin];
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
PairLL avgPoint{avgLon, avgLat};
return {eSample, avgPoint, SpreadLon, SpreadLat, MaxSpread};
}
std::vector<KTreeElt> SplitKTreeElt(KTreeElt const& eKD)
{
std::vector<PairLL> NewListPt1;
std::vector<PairLL> NewListPt2;
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
KTreeElt eKD1{{}, eKD.eDisc, true, -2, -2};
KTreeElt eKD2{NewListPt1, KTree_ComputeDisc(NewListPt1), false, -1, -1};
KTreeElt eKD3{NewListPt2, KTree_ComputeDisc(NewListPt2), false, -1, -1};
return {eKD1, eKD2, eKD3};
}
std::vector<KTreeElt> KTree_GetDecomposition(std::vector<PairLL> const& ListPtCoast)
{
std::vector<KTreeElt> TList;
std::vector<int> IsTreated;
int MaxNumberPerCell=100;
auto SplitComponent=[&](int const& iComp) -> void {
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
return TList;
}
double ShortestDistance(std::vector<KTreeElt> const& ListKT, PairLL const& ePt, double & UpperEstimate)
{
double TheDist=UpperEstimate;
double eLon=ePt.eLon;
double eLat=ePt.eLat;
std::vector<int> ListIdx{0};
while(true) {
if (ListIdx.size() == 0)
break;
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
std::vector<double> GetUpperEstimateMinDist(std::vector<PairLL> const& ListPt1, std::vector<PairLL> const& ListPt2)
{
int nbPt2=ListPt2.size();
int nbPt1=ListPt1.size();
if (nbPt1 == 0) {
std::cerr << "The list ListPt1 should not be empty\n";
std::cerr << "nbPt1=" << nbPt1 << "\n";
throw TerminalException{1};
}
auto fDist=[&](int const& i1, int const& i2) -> double {
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
for (int iter=0; iter<4; iter++) {
int i1New =rand() % nbPt1;
double nDist=fDist(i1New, i2);
if (nDist < eDist) {
eDist=nDist;
i1=i1New;
}
}
while(true) {
bool DoSomething=false;
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
return ListUpperEst;
}
std::vector<double> GetListMinimalDistances(std::vector<PairLL> const& ListPtCoast, std::vector<PairLL> const& ListPt)
{
int nbPt=ListPt.size();
std::vector<double> ListUpperEst=GetUpperEstimateMinDist(ListPtCoast, ListPt);
std::vector<KTreeElt> ListKT=KTree_GetDecomposition(ListPtCoast);
std::vector<double> ListShortest(nbPt);
double TotalDefect=0;
for (int iPt=0; iPt<nbPt; iPt++) {
double eEst=ListUpperEst[iPt];
PairLL ePt=ListPt[iPt];
double eMinDist=ShortestDistance(ListKT, ePt, eEst);
ListShortest[iPt]=eMinDist;
double eDefect=eEst - eMinDist;
TotalDefect += eDefect;
}
return ListShortest;
}
bool NC_IsVar(std::string const& eFile, std::string const& eVar)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "Error in NC_IsVar\n";
std::cerr << "Trying to open non-existing file\n";
std::cerr << "eFile = " << eFile << "\n";
throw TerminalException{1};
}
try {
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(eVar);
if(data.isNull()) {
return false;
}
return true;
}
catch (...) {
return false;
}
}
MyMatrix<double> NC_Read2Dvariable(std::string const& eFile, std::string const& eVar)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "Error in NC_Read2Dvariable\n";
std::cerr << "Trying to open non-existing file\n";
std::cerr << "eFile = " << eFile << "\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(eVar);
netCDF::NcType eType=data.getType();
if (data.isNull()) {
std::cerr << "Error in accessing to the file (Case 1)\n";
std::cerr << "eFile = " << eFile << "\n";
std::cerr << "eVar  = " << eVar << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim != 2) {
std::cerr << "The number of dimensions is not correct\n";
throw TerminalException{1};
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
throw TerminalException{1};
}
return eArr;
}
MyMatrix<int> NC_Read2Dvariable_int(std::string const& eFile, std::string const& eVar)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "Error in NC_Read2Dvariable_int\n";
std::cerr << "Trying to open non-existing file\n";
std::cerr << "eFile = " << eFile << "\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(eVar);
netCDF::NcType eType=data.getType();
if(data.isNull()) {
std::cerr << "Error in accessing to the file (Case 2)\n";
std::cerr << "eFile = " << eFile << "\n";
std::cerr << "eVar  = " << eVar << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim != 2) {
std::cerr << "The number of dimensions is not correct\n";
throw TerminalException{1};
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
throw TerminalException{1};
}
return eArr;
}
MyVector<double> NC_Read1Dvariable(std::string const& eFile, std::string const& eVar)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "Error in NC_Read1Dvariable\n";
std::cerr << "Trying to open non-existing file\n";
std::cerr << "eFile = " << eFile << "\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(eVar);
netCDF::NcType eType=data.getType();
if(data.isNull()) {
std::cerr << "Error in accessing to the file (Case 3)\n";
std::cerr << "eFile = " << eFile << "\n";
std::cerr << "eVar  = " << eVar << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim != 1) {
std::cerr << "The number of dimensions is not correct\n";
throw TerminalException{1};
}
netCDF::NcDim eDim=data.getDim(0);
int dim=eDim.getSize();
MyVector<double> eArr(dim);
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
throw TerminalException{1};
}
double eScal, eOff;
try {
netCDF::NcVarAtt eScalAtt=data.getAtt("scale_factor");
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
MyVector<int> NC_Read1Dvariable_int(std::string const& eFile, std::string const& eVar)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "Error in NC_Read1Dvariable_int\n";
std::cerr << "Trying to open non-existing file\n";
std::cerr << "eFile = " << eFile << "\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(eVar);
netCDF::NcType eType=data.getType();
if(data.isNull()) {
std::cerr << "Error in accessing to the file (Case 4)\n";
std::cerr << "eFile = " << eFile << "\n";
std::cerr << "eVar  = " << eVar << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim != 1) {
std::cerr << "The number of dimensions is not correct\n";
throw TerminalException{1};
}
netCDF::NcDim eDim=data.getDim(0);
int dim=eDim.getSize();
MyVector<int> eArr(dim);
int *eValINT;
eValINT=new int[dim];
data.getVar(eValINT);
for (int i=0; i<dim; i++) {
int eValI=eValINT[i];
eArr(i)=eValI;
}
delete [] eValINT;
return eArr;
}
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
throw TerminalException{1};
}
int alen=eStrUnitTime.length();
std::string Yname=eStrUnitTime.substr(posBlank+1, alen - 1 - posBlank);
int alenB=Yname.length();
int posBlankB=STRING_GetCharPositionInString(Yname, strSpace);
std::string YnameB=Yname.substr(posBlankB+1, alenB - 1 - posBlankB);
std::string strT="T";
std::string strZ="Z";
std::vector<std::string> LStrDateT=STRING_Split(YnameB, strT);
int sizStrDateT=LStrDateT.size();
if (sizStrDateT > 1) {
YnameDate=LStrDateT[0];
std::string eStrB=LStrDateT[1];
int alenC=eStrUnitTime.length();
YnameTime=eStrB.substr(0,alenC-2);
std::cerr << "Case of WW3\n";
std::cerr << "YnameDate=" << YnameDate << "\n";
std::cerr << "YnameTime=" << YnameTime << "\n";
}
else {
std::vector<std::string> LStrDate=STRING_Split(YnameB, strSpace);
int sizStrDate=LStrDate.size();
if (sizStrDate > 1) {
YnameDate=LStrDate[0];
YnameTime=LStrDate[1];
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
throw TerminalException{1};
}
}
}
else {
YnameDate=LStrDate[0];
YnameTime="00:00:00";
}
}
std::vector<std::string> eVectDate=STRING_Split(YnameDate,"-");
std::string eStrYear, eStrMonth, eStrDay;
eStrYear=eVectDate[0];
eStrMonth=eVectDate[1];
eStrDay=eVectDate[2];
int year, month, day;
std::istringstream(eStrYear) >> year;
std::istringstream(eStrMonth) >> month;
std::istringstream(eStrDay) >> day;
std::vector<std::string> eVectTime=STRING_Split(YnameTime,":");
std::string eStrHour, eStrMin, eStrSec;
eStrHour=eVectTime[0];
eStrMin=eVectTime[1];
eStrSec=eVectTime[2];
int hour, min, sec;
std::istringstream(eStrHour) >> hour;
std::istringstream(eStrMin) >> min;
std::istringstream(eStrSec) >> sec;
eTimeStart=DATE_ConvertSix2mjd({year, month, day, hour, min, sec});
}
std::vector<double> NC_ReadTimeFromFile(std::string const& eFile, std::string const& StringTime)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "Error in NC_ReadTimeFromFile\n";
std::cerr << "Trying to open non-existing file\n";
std::cerr << "eFile = " << eFile << "\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(StringTime);
if(data.isNull()) {
std::cerr << "Error in accessing to the file (Case 5)\n";
std::cerr << "eFile = " << eFile << "\n";
std::cerr << "StringTime = " << StringTime << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim != 1) {
std::cerr << "The number of dimensions is not correct\n";
throw TerminalException{1};
}
netCDF::NcDim eDim=data.getDim(0);
int siz=eDim.getSize();
double *eVal;
eVal=new double[siz];
data.getVar(eVal);
netCDF::NcVarAtt eTimeAtt=data.getAtt("units");
char eString[1024]="";
eTimeAtt.getValues(eString);
std::string eStrUnitTime=eString;
double ConvertToDay, eTimeStart;
CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart);
std::vector<double> LTime(siz);
if (siz == 0) {
std::cerr << "We found siz=0\n";
std::cerr << "This means that we find zero times in the file\n";
std::cerr << "which is not what we expected\n";
throw TerminalException{1};
}
double minTime=0, maxTime=0;
for (int i=0; i<siz; i++) {
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
LTime[i]=eTimeDay;
}
delete [] eVal;
bool ShowMinMax=false;
if (ShowMinMax) {
std::cerr << "minTime=" << minTime << " maxTime=" << maxTime << "\n";
std::string strPresMin=DATE_ConvertMjd2mystringPres(minTime);
std::string strPresMax=DATE_ConvertMjd2mystringPres(maxTime);
std::cerr << "strPresMin=" << strPresMin << "\n";
std::cerr << "strPresMax=" << strPresMax << "\n";
}
return LTime;
}
MyMatrix<double> NETCDF_Get2DvariableSpecEntry_FD(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "NETCDF_Get2DvariableSpecEntry_FD\n";
std::cerr << "The file eFile=" << eFile << "\n";
std::cerr << "does not exist\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(eVar);
if(data.isNull()) {
std::cerr << "Error in accessing to the file (Case 6)\n";
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "eVar=" << eVar << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim == 3) {
netCDF::NcDim eDim;
eDim=data.getDim(0);
int nbRec=eDim.getSize();
if (iRec >= nbRec) {
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "Error, iRec is too large (Case 1)\n";
std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
std::cerr << "We need C-convention iRec < nbRec\n";
throw TerminalException{1};
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
throw TerminalException{1};
}
return eArr;
}
netCDF::NcDim eDim;
eDim=data.getDim(0);
int nbRec=eDim.getSize();
if (iRec >= nbRec) {
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "Error, iRec is too large (Case 2)\n";
std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
std::cerr << "We need C-convention iRec < nbRec\n";
throw TerminalException{1};
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
throw TerminalException{1};
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
std::cerr << "Routine is NETCDF_Get2DvariableSpecEntry_FD\n";
std::cerr << "eVar=" << eVar << "\n";
std::cerr << "We did not find the size\n";
throw TerminalException{1};
}
MyMatrix<double> NETCDF_Get2DvariableSpecEntry_FE(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "NETCDF_Get2DvariableSpecEntry_FE\n";
std::cerr << "eVar=" << eVar << " iRec=" << iRec << "\n";
std::cerr << "The file eFile=" << eFile << "\n";
std::cerr << "does not exist\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(eVar);
if(data.isNull()) {
std::cerr << "Error in accessing to the file (Case 7)\n";
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "eVar=" << eVar << "\n";
std::cerr << "iRec=" << iRec << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim != 2) {
std::cerr << "This command will certainly not work\n";
std::cerr << "Dimensions are not correct\n";
throw TerminalException{1};
}
netCDF::NcDim eDim;
eDim=data.getDim(0);
int nbRec=eDim.getSize();
if (iRec >= nbRec) {
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "Error, iRec is too large (Case 3)\n";
std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
std::cerr << "We need C-convention iRec < nbRec\n";
throw TerminalException{1};
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
if (eType == netCDF::NcType::nc_INT) {
int *eVal;
eVal=new int[mnp];
data.getVar(start, count, eVal);
for (size_t i=0; i<mnp; i++) {
int eValF=eVal[i];
double eValD=double(eValF);
eArr(i,0)=eValD;
}
delete [] eVal;
IsDone=true;
}
if (eType == netCDF::NcType::nc_UINT) {
unsigned int *eVal;
eVal=new unsigned int[mnp];
data.getVar(start, count, eVal);
for (size_t i=0; i<mnp; i++) {
unsigned int eValF=eVal[i];
double eValD=double(eValF);
eArr(i,0)=eValD;
}
delete [] eVal;
IsDone=true;
}
if (IsDone == false) {
std::cerr << "Data reading failed for 2D finite element\n";
throw TerminalException{1};
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
Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecEntry_FE(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "NETCDF_Get2DvariableSpecEntry_FE\n";
std::cerr << "The file eFile=" << eFile << "\n";
std::cerr << "does not exist\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcVar data=dataFile.getVar(eVar);
if(data.isNull()) {
std::cerr << "Error in accessing to the file (Case 8)\n";
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "eVar=" << eVar << "\n";
std::cerr << "iRec=" << iRec << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim != 3) {
std::cerr << "This command will certainly not work\n";
std::cerr << "Dimensions are not correct\n";
throw TerminalException{1};
}
netCDF::NcDim eDim0=data.getDim(0);
int nbRec=eDim0.getSize();
if (iRec >= nbRec) {
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "Error, iRec is too large (Case 4)\n";
std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
std::cerr << "We need C-convention iRec < nbRec\n";
throw TerminalException{1};
}
netCDF::NcDim eDim1=data.getDim(1);
size_t mnp=eDim1.getSize();
netCDF::NcDim eDim2=data.getDim(2);
size_t NTR=eDim2.getSize();
std::vector<size_t> start{size_t(iRec), 0, 0};
std::vector<size_t> count{1, mnp, NTR};
netCDF::NcType eType=data.getType();
Eigen::Tensor<double,3> eArr(int(NTR), int(mnp), 1);
bool IsDone=false;
if (eType == netCDF::NcType::nc_DOUBLE) {
double *eVal;
eVal=new double[mnp*NTR];
data.getVar(start, count, eVal);
int idx=0;
for (size_t i=0; i<mnp; i++) {
for (size_t iTr=0; iTr<NTR; iTr++) {
eArr(iTr,i,0)=eVal[idx];
idx++;
}
}
delete [] eVal;
IsDone=true;
}
if (eType == netCDF::NcType::nc_FLOAT) {
float *eVal;
eVal=new float[mnp*NTR];
data.getVar(start, count, eVal);
int idx=0;
for (size_t i=0; i<mnp; i++) {
for (size_t iTr=0; iTr<NTR; iTr++) {
float eValF=eVal[idx];
double eValD=double(eValF);
eArr(iTr,i,0)=eValD;
idx++;
}
}
delete [] eVal;
IsDone=true;
}
if (eType == netCDF::NcType::nc_INT) {
int *eVal;
eVal=new int[mnp*NTR];
data.getVar(start, count, eVal);
int idx=0;
for (size_t i=0; i<mnp; i++) {
for (size_t iTr=0; iTr<NTR; iTr++) {
int eValF=eVal[idx];
double eValD=double(eValF);
eArr(iTr,i,0)=eValD;
idx++;
}
}
delete [] eVal;
IsDone=true;
}
if (eType == netCDF::NcType::nc_UINT) {
unsigned int *eVal;
eVal=new unsigned int[mnp];
data.getVar(start, count, eVal);
int idx=0;
for (size_t i=0; i<mnp; i++) {
for (size_t iTr=0; iTr<NTR; iTr++) {
unsigned int eValF=eVal[idx];
double eValD=double(eValF);
eArr(iTr,i,0)=eValD;
idx++;
}
}
delete [] eVal;
IsDone=true;
}
if (IsDone == false) {
std::cerr << "Data reading failed for 2D finite element\n";
throw TerminalException{1};
}
if (GrdArr.L_IndexSelect) {
int siz=GrdArr.I_IndexSelect.size();
Eigen::Tensor<double,3> eArrRet(siz, int(NTR), 1);
for (int i=0; i<siz; i++) {
int iGlob=GrdArr.I_IndexSelect[i];
for (size_t iTr=0; iTr<NTR; iTr++)
eArrRet(iTr,i,0)=eArr(iTr,iGlob,0);
}
return eArrRet;
}
return eArr;
}
MyMatrix<double> NETCDF_Get2DvariableSpecEntry(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
if (GrdArr.IsFE == 1) {
return NETCDF_Get2DvariableSpecEntry_FE(eFile, GrdArr, eVar, iRec);
}
return NETCDF_Get2DvariableSpecEntry_FD(eFile, GrdArr, eVar, iRec);
}
MyMatrix<double> NETCDF_Get2DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& eVar, double const& eTimeDay)
{
ArrayHistory eArr=TotalArr.eArr;
GridArray GrdArr=TotalArr.GrdArr;
auto GetHisFileName=[&](int const& iFile) -> std::string {
if (eArr.TimeSteppingInfo == "multiplenetcdf") {
return eArr.HisPrefix + StringNumber(iFile+1,eArr.nbDigit) + ".nc";
}
int len=eArr.ListFileNames.size();
if (iFile >= len) {
std::cerr << "iFile=" << iFile << " len=" << len << "\n";
std::cerr << "We need iFile < len\n";
std::cerr << "Error. trying to get eArr.ListFileNames\n";
std::cerr << "After the last values\n";
throw TerminalException{1};
}
if (eArr.AppendVarName) {
return eArr.ListFileNames[iFile] + eVar + ".nc";
}
else {
return eArr.ListFileNames[iFile];
}
};
std::string HisFileLow, HisFileUpp;
int iRecLow, iRecUpp;
double alphaLow=0, alphaUpp=0;
bool IsDone=false;
if (eArr.TimeSteppingInfo == "classic") {
InterpInfo eInterpInfo=GetTimeInterpolationInfo(eArr.ListTime, eTimeDay);
if (eInterpInfo.UseSingleEntry) {
int iTime=eInterpInfo.iTimeLow;
int iFile=eArr.ListIFile[iTime];
int iRec=eArr.ListIRec[iTime];
std::string HisFile=GetHisFileName(iFile);
return NETCDF_Get2DvariableSpecEntry(HisFile, GrdArr, eVar, iRec);
}
alphaLow=eInterpInfo.alphaLow;
alphaUpp=eInterpInfo.alphaUpp;
int iTimeLow=eInterpInfo.iTimeLow;
int iTimeUpp=eInterpInfo.iTimeUpp;
int iFileLow=eArr.ListIFile[iTimeLow];
int iFileUpp=eArr.ListIFile[iTimeUpp];
iRecLow=eArr.ListIRec[iTimeLow];
iRecUpp=eArr.ListIRec[iTimeUpp];
HisFileLow=GetHisFileName(iFileLow);
HisFileUpp=GetHisFileName(iFileUpp);
IsDone=true;
}
if (eArr.TimeSteppingInfo == "singlefile") {
InterpInfo eInterpInfo=GetTimeInterpolationInfo_infinite(eArr.FirstTime, eArr.SeparationTime, eTimeDay);
if (eInterpInfo.UseSingleEntry) {
int iTime=eInterpInfo.iTimeLow;
int iFile=0;
int iRec=iTime;
std::string HisFile=GetHisFileName(iFile);
return NETCDF_Get2DvariableSpecEntry(HisFile, GrdArr, eVar, iRec);
}
alphaLow=eInterpInfo.alphaLow;
alphaUpp=eInterpInfo.alphaUpp;
int iTimeLow=eInterpInfo.iTimeLow;
int iTimeUpp=eInterpInfo.iTimeUpp;
int iFileLow=0;
int iFileUpp=0;
iRecLow=iTimeLow;
iRecUpp=iTimeUpp;
HisFileLow=GetHisFileName(iFileLow);
HisFileUpp=GetHisFileName(iFileUpp);
IsDone=true;
}
if (eArr.TimeSteppingInfo == "multiplenetcdf") {
InterpInfo eInterpInfo=GetTimeInterpolationInfo_infinite(eArr.FirstTime, eArr.SeparationTime, eTimeDay);
if (eInterpInfo.UseSingleEntry) {
int iTime=eInterpInfo.iTimeLow;
std::vector<int> eRec=GetIFileIRec(eArr.nbRecBegin, eArr.nbRecMiddle, iTime);
int iFile=eRec[0];
int iRec=eRec[1];
std::string HisFile=GetHisFileName(iFile);
return NETCDF_Get2DvariableSpecEntry(HisFile, GrdArr, eVar, iRec);
}
alphaLow=eInterpInfo.alphaLow;
alphaUpp=eInterpInfo.alphaUpp;
int iTimeLow=eInterpInfo.iTimeLow;
int iTimeUpp=eInterpInfo.iTimeUpp;
std::vector<int> eRecLow=GetIFileIRec(eArr.nbRecBegin, eArr.nbRecMiddle, iTimeLow);
std::vector<int> eRecUpp=GetIFileIRec(eArr.nbRecBegin, eArr.nbRecMiddle, iTimeUpp);
int iFileLow=eRecLow[0];
int iFileUpp=eRecUpp[0];
iRecLow=eRecLow[1];
iRecUpp=eRecUpp[1];
HisFileLow=GetHisFileName(iFileLow);
HisFileUpp=GetHisFileName(iFileUpp);
IsDone=true;
}
if (IsDone == false) {
std::cerr << "Failed to find matching entry for TimeSteppingInfo = " << eArr.TimeSteppingInfo << "\n";
throw TerminalException{1};
}
MyMatrix<double> eVarLow=NETCDF_Get2DvariableSpecEntry(HisFileLow, GrdArr, eVar, iRecLow);
MyMatrix<double> eVarUpp=NETCDF_Get2DvariableSpecEntry(HisFileUpp, GrdArr, eVar, iRecUpp);
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
Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecEntry_FD(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
if (IsExistingFile(eFile) == false) {
std::cerr << "NETCDF_Get3DvariableSpecEntry_FD\n";
std::cerr << "The file eFile=" << eFile << "\n";
std::cerr << "does not exist\n";
throw TerminalException{1};
}
netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
netCDF::NcDim eDim;
eDim=dataFile.getDim("s_rho");
int s_rho=eDim.getSize();
eDim=dataFile.getDim("s_w");
int s_w=eDim.getSize();
netCDF::NcVar data=dataFile.getVar(eVar);
if(data.isNull()) {
std::cerr << "Error in accessing to the file (Case 10)\n";
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "eVar=" << eVar << "\n";
throw TerminalException{1};
}
int nbDim=data.getDimCount();
if (nbDim == 4) {
eDim=data.getDim(0);
int nbRec=eDim.getSize();
if (iRec >= nbRec) {
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "Error, iRec is too large (Case 6)\n";
std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
std::cerr << "We need C-convention iRec < nbRec\n";
throw TerminalException{1};
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
Eigen::Tensor<double,3> eArr(s_vert, eta, xi);
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
throw TerminalException{1};
}
int idx=0;
for (int k=0; k<s_vert; k++)
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++) {
eArr(k, i, j)=eVal[idx];
idx++;
}
return eArr;
}
eDim=data.getDim(0);
int nbRec=eDim.getSize();
if (iRec >= nbRec) {
std::cerr << "eFile=" << eFile << "\n";
std::cerr << "Error, iRec is too large (Case 7)\n";
std::cerr << "iRec=" << iRec << " nbRec=" << nbRec << "\n";
std::cerr << "We need C-convention iRec < nbRec\n";
throw TerminalException{1};
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
throw TerminalException{1};
}
if (nbWet == s_rho*GrdArr.GrdArrRho.nbWet) {
int eta=GrdArr.GrdArrRho.eta;
int xi=GrdArr.GrdArrRho.xi;
Eigen::Tensor<double,3> eArr(s_rho, eta, xi);
for (int k=0; k<s_rho; k++)
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++)
eArr(k, i, j)=0;
int idx=0;
for (int k=0; k<s_rho; k++) {
for (int iWet=0; iWet<GrdArr.GrdArrRho.nbWet; iWet++) {
int i=GrdArr.GrdArrRho.Idx[iWet];
int j=GrdArr.GrdArrRho.Jdx[iWet];
eArr(k, i, j)=eVal[idx];
idx++;
}
}
return eArr;
}
if (nbWet == s_rho*GrdArr.GrdArrU.nbWet) {
int eta=GrdArr.GrdArrU.eta;
int xi=GrdArr.GrdArrU.xi;
Eigen::Tensor<double,3> eArr(s_rho, eta, xi);
for (int k=0; k<s_rho; k++)
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++)
eArr(k, i, j)=0;
int idx=0;
for (int k=0; k<s_rho; k++) {
for (int iWet=0; iWet<GrdArr.GrdArrU.nbWet; iWet++) {
int i=GrdArr.GrdArrU.Idx[iWet];
int j=GrdArr.GrdArrU.Jdx[iWet];
eArr(k, i, j)=eVal[idx];
idx++;
}
}
return eArr;
}
if (nbWet == s_rho*GrdArr.GrdArrV.nbWet) {
int eta=GrdArr.GrdArrV.eta;
int xi=GrdArr.GrdArrV.xi;
Eigen::Tensor<double,3> eArr(s_rho, eta, xi);
for (int k=0; k<s_rho; k++)
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++)
eArr(k, i, j)=0;
int idx=0;
for (int k=0; k<s_rho; k++) {
for (int iWet=0; iWet<GrdArr.GrdArrV.nbWet; iWet++) {
int i=GrdArr.GrdArrV.Idx[iWet];
int j=GrdArr.GrdArrV.Jdx[iWet];
eArr(k, i, j)=eVal[idx];
idx++;
}
}
return eArr;
}
std::cerr << "Routine is NETCDF_Get2DvariableSpecEntry\n";
std::cerr << "eVar=" << eVar << "\n";
std::cerr << "s_rho=" << s_rho << " s_w=" << s_w << "\n";
std::cerr << "nbWet=" << nbWet << "\n";
std::cerr << "nbWetRho=" << GrdArr.GrdArrRho.nbWet << "\n";
std::cerr << "  nbWetU=" << GrdArr.GrdArrU.nbWet << "\n";
std::cerr << "  nbWetV=" << GrdArr.GrdArrV.nbWet << "\n";
std::cerr << "We did not find the size\n";
throw TerminalException{1};
}
Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecEntry(std::string const& eFile, GridArray const& GrdArr, std::string const& eVar, int const& iRec)
{
if (GrdArr.IsFE == 1) {
return NETCDF_Get3DvariableSpecEntry_FE(eFile, GrdArr, eVar, iRec);
}
return NETCDF_Get3DvariableSpecEntry_FD(eFile, GrdArr, eVar, iRec);
}
Eigen::Tensor<double,3> NETCDF_Get3DvariableSpecTime(TotalArrGetData const& TotalArr, std::string const& eVar, double const& eTimeDay)
{
ArrayHistory eArr=TotalArr.eArr;
GridArray GrdArr=TotalArr.GrdArr;
InterpInfo eInterpInfo=GetTimeInterpolationInfo(eArr.ListTime, eTimeDay);
if (eInterpInfo.UseSingleEntry) {
double iTime=eInterpInfo.iTimeLow;
int iFile=eArr.ListIFile[iTime];
int iRec=eArr.ListIRec[iTime];
std::string HisFile=eArr.ListFileNames[iFile];
return NETCDF_Get3DvariableSpecEntry(HisFile, GrdArr, eVar, iRec);
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
Eigen::Tensor<double,3> eVarLow=NETCDF_Get3DvariableSpecEntry(HisFileLow, GrdArr, eVar, iRecLow);
Eigen::Tensor<double,3> eVarUpp=NETCDF_Get3DvariableSpecEntry(HisFileUpp, GrdArr, eVar, iRecUpp);
auto LDim=eVarLow.dimensions();
int s_vert=LDim[0];
int eta=LDim[1];
int xi=LDim[2];
Eigen::Tensor<double,3> RetVar(s_vert, eta, xi);
for (int k=0; k<s_vert; k++)
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++)
RetVar(k, i, j)=alphaLow*eVarLow(k, i, j) + alphaUpp*eVarUpp(k, i, j);
return RetVar;
}
QuadArray GetQuadArray(GridArray const& GrdArr)
{
double MinLon=0, MaxLon=0, MinLat=0, MaxLat=0;
if (GrdArr.IsFE == 1) {
MinLon=GrdArr.GrdArrRho.LON.minCoeff();
MaxLon=GrdArr.GrdArrRho.LON.maxCoeff();
MinLat=GrdArr.GrdArrRho.LAT.minCoeff();
MaxLat=GrdArr.GrdArrRho.LAT.maxCoeff();
}
else {
bool IsFirst=true;
int eta_rho=GrdArr.GrdArrRho.LON.rows();
int xi_rho =GrdArr.GrdArrRho.LON.cols();
int eta_rho_msk=GrdArr.GrdArrRho.MSK.rows();
int xi_rho_msk =GrdArr.GrdArrRho.MSK.cols();
if (eta_rho_msk != eta_rho || xi_rho_msk != xi_rho) {
std::cerr << "Dimension error in the arrays\n";
throw TerminalException{1};
}
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++)
if (GrdArr.GrdArrRho.MSK(i,j) == 1) {
double eLon=GrdArr.GrdArrRho.LON(i,j);
double eLat=GrdArr.GrdArrRho.LAT(i,j);
if (IsFirst == true) {
MinLon=eLon;
MaxLon=eLon;
MinLat=eLat;
MaxLat=eLat;
IsFirst=false;
}
else {
if (eLon < MinLon)
MinLon=eLon;
if (eLon > MaxLon)
MaxLon=eLon;
if (eLat < MinLat)
MinLat=eLat;
if (eLat > MaxLat)
MaxLat=eLat;
}
}
}
return {MinLon, MaxLon, MinLat, MaxLat};
}
void CHECK_Model_Allowedness(std::string const& eModelName)
{
std::vector<std::string> vec=GetAllPossibleModels();
bool isPresent = (std::find(vec.begin(), vec.end(), eModelName) != vec.end());
if (isPresent == false) {
std::cerr << "We did not find the MODEL NAME\n";
std::cerr << "MODELNAME = " << eModelName << "\n";
std::cerr << "List of allowed models\n";
for (int iModel=0; iModel<int(vec.size()); iModel++) {
std::cerr << "iModel=" << iModel << " eModel=" << vec[iModel] << "\n";
}
throw TerminalException{1};
}
}
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
GridArray NC_ReadRomsGridFile(std::string const& eFile)
{
std::function<int(double const&)> fConv=[](double const& x) -> int {
return int(x);
};
GridArray eRomsGridArray;
eRomsGridArray.ModelName="ROMS";
eRomsGridArray.IsFE=0;
eRomsGridArray.IsSpherical=true;
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
eRomsGridArray.GrdArrRho.MSK=ConvertMatrix(eMSK_rho_double, fConv);
InitializeIdxJdxWet(eRomsGridArray.GrdArrRho);
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
(eRomsGridArray.GrdArrRho.DEP(i ,j+1)+
eRomsGridArray.GrdArrRho.DEP(i+1,j+1)+
eRomsGridArray.GrdArrRho.DEP(i ,j )+
eRomsGridArray.GrdArrRho.DEP(i+1,j ))/double(4);
ANGp(i, j)=
(eRomsGridArray.GrdArrRho.ANG(i ,j+1)+
eRomsGridArray.GrdArrRho.ANG(i+1,j+1)+
eRomsGridArray.GrdArrRho.ANG(i ,j )+
eRomsGridArray.GrdArrRho.ANG(i+1,j ))/double(4);
MSKp(i, j)=
eRomsGridArray.GrdArrRho.MSK(i ,j+1)*
eRomsGridArray.GrdArrRho.MSK(i+1,j+1)*
eRomsGridArray.GrdArrRho.MSK(i ,j )*
eRomsGridArray.GrdArrRho.MSK(i+1,j );
}
eRomsGridArray.GrdArrPsi.MSK=MSKp;
eRomsGridArray.GrdArrPsi.DEP=DEPp;
eRomsGridArray.GrdArrPsi.HaveDEP=true;
eRomsGridArray.GrdArrPsi.ANG=ANGp;
InitializeIdxJdxWet(eRomsGridArray.GrdArrPsi);
std::cerr << "The ROMS grid has been read\n";
return eRomsGridArray;
}
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
eCosmoWamGridArray.IsSpherical=true;
std::string LONstr="LON_" + postfix;
eCosmoWamGridArray.GrdArrRho.LON=NC_Read2Dvariable(eFile, LONstr);
std::string LATstr="LAT_" + postfix;
eCosmoWamGridArray.GrdArrRho.LAT=NC_Read2Dvariable(eFile, LATstr);
int eta_rho=eCosmoWamGridArray.GrdArrRho.LON.rows();
int xi_rho=eCosmoWamGridArray.GrdArrRho.LON.cols();
std::string DEPstr="DEP_" + postfix;
if (NC_IsVar(eFile, DEPstr) ) {
eCosmoWamGridArray.GrdArrRho.DEP=NC_Read2Dvariable(eFile, DEPstr);
eCosmoWamGridArray.GrdArrRho.HaveDEP=true;
}
else {
eCosmoWamGridArray.GrdArrRho.DEP=ZeroMatrix<double>(eta_rho, xi_rho);
eCosmoWamGridArray.GrdArrRho.HaveDEP=false;
}
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
if (postfix == "atm") {
eCosmoWamGridArray.GrdArrRho.ANG=NC_Read2Dvariable(eFile, "ANG_atm");
}
else {
eCosmoWamGridArray.GrdArrRho.ANG=CreateAngleMatrix(eCosmoWamGridArray.GrdArrRho.LON, eCosmoWamGridArray.GrdArrRho.LAT);
}
eCosmoWamGridArray.GrdArrRho.eta=eta_rho;
eCosmoWamGridArray.GrdArrRho.xi =xi_rho;
return eCosmoWamGridArray;
}
GridArray NC_ReadSCHISM_sflux_grid(std::string const& eFile)
{
GridArray GrdArr;
GrdArr.ModelName="SCHISM_SFLUX";
GrdArr.IsFE=0;
GrdArr.IsSpherical=true;
std::string LONstr="lon";
GrdArr.GrdArrRho.LON=NC_Read2Dvariable(eFile, LONstr);
std::string LATstr="lat";
GrdArr.GrdArrRho.LAT=NC_Read2Dvariable(eFile, LATstr);
int eta_rho=GrdArr.GrdArrRho.LON.rows();
int xi_rho=GrdArr.GrdArrRho.LON.cols();
GrdArr.GrdArrRho.HaveDEP=false;
MyMatrix<int> MSK_int(eta_rho, xi_rho);
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++)
MSK_int(i,j)=1;
GrdArr.GrdArrRho.MSK=MSK_int;
GrdArr.GrdArrRho.ANG=CreateAngleMatrix(GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT);
GrdArr.GrdArrRho.eta=eta_rho;
GrdArr.GrdArrRho.xi =xi_rho;
return GrdArr;
}
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
GridArray NC_ReadWamGridFile(std::string const& eFile)
{
MyVector<int> eVectLLUNSTR=NC_Read1Dvariable_int(eFile, "LLUNSTR");
int LLUNSTR=eVectLLUNSTR(0);
if (LLUNSTR == 0)
return NC_ReadCosmoWamStructGridFile(eFile, "wav");
GridArray GrdArr;
GrdArr.ModelName="WAM";
GrdArr.IsFE=1;
GrdArr.IsSpherical=true;
GrdArr.L_IndexSelect=false;
GrdArr.INE=NC_ReadElements(eFile, "ele");
MyVector<double> LON=NC_Read1Dvariable(eFile, "LON_wav");
MyVector<double> LAT=NC_Read1Dvariable(eFile, "LAT_wav");
MyVector<double> DEP=NC_Read1Dvariable(eFile, "DEP_wav");
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
GrdArr.GrdArrRho.LON=LONarr;
GrdArr.GrdArrRho.LAT=LATarr;
GrdArr.GrdArrRho.DEP=DEParr;
GrdArr.GrdArrRho.ANG=ANGarr;
GrdArr.GrdArrRho.MSK=MSKarr;
return GrdArr;
}
GridArray WWM_ReadGridFile_netcdf(std::string const& GridFile)
{
GridArray GrdArr;
GrdArr.ModelName="WWM";
GrdArr.IsFE=1;
GrdArr.L_IndexSelect=false;
if (IsExistingFile(GridFile) == false) {
std::cerr << "Error in WWM_ReadGridFile_netcdf\n";
std::cerr << "GridFile = " << GridFile << "\n";
std::cerr << "is missing\n";
throw TerminalException{1};
}
GrdArr.INE=NC_ReadElements(GridFile, "ele");
MyVector<int> LType=NC_Read1Dvariable_int(GridFile, "LSPHE");
int LSPHE=LType(0);
std::string Xname, Yname;
if (LSPHE == 1) {
Xname="lon";
Yname="lat";
GrdArr.IsSpherical=true;
}
else {
Xname="x";
Yname="y";
GrdArr.IsSpherical=false;
}
MyVector<double> LON=NC_Read1Dvariable(GridFile, Xname);
MyVector<double> LAT=NC_Read1Dvariable(GridFile, Yname);
MyVector<double> DEP=NC_Read1Dvariable(GridFile, "depth");
MyVector<double> IOBP=NC_Read1Dvariable(GridFile, "IOBP");
int nbPoint=LON.size();
MyMatrix<double> LONarr(nbPoint,1);
MyMatrix<double> LATarr(nbPoint,1);
MyMatrix<double> DEParr(nbPoint,1);
MyMatrix<double> ANGarr(nbPoint,1);
MyMatrix<int> MSKarr(nbPoint,1);
MyVector<int> IOBParr(nbPoint);
for (int iPoint=0; iPoint<nbPoint; iPoint++) {
LONarr(iPoint,0)=LON(iPoint);
LATarr(iPoint,0)=LAT(iPoint);
DEParr(iPoint,0)=DEP(iPoint);
ANGarr(iPoint,0)=0;
MSKarr(iPoint,0)=1;
IOBParr(iPoint)=int(IOBP(iPoint));
}
GrdArr.GrdArrRho.LON=LONarr;
GrdArr.GrdArrRho.LAT=LATarr;
GrdArr.GrdArrRho.DEP=DEParr;
GrdArr.GrdArrRho.ANG=ANGarr;
GrdArr.GrdArrRho.MSK=MSKarr;
GrdArr.IOBP=IOBParr;
return GrdArr;
}
MyVector<int> WWM_ReadBoundFile_gr3(std::string const& BoundFile)
{
if (IsExistingFile(BoundFile) == false) {
std::cerr << "Error in WWM_ReadBoundFile_gr3\n";
std::cerr << "Missing BoundFile=" << BoundFile << "\n";
throw TerminalException{1};
}
std::ifstream IN(BoundFile);
std::string line;
std::getline(IN, line);
int mne, mnp;
IN >> mne;
IN >> mnp;
MyVector<int> eVect(mnp);
for (int i=0; i<mnp; i++) {
int KTMP;
double XPDTMP, YPDTMP, ZPDTMP;
IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
if (KTMP != i+1) {
std::cerr << "Inconsistency at this level\n";
throw TerminalException{1};
}
int eIOBP=int(ZPDTMP);
eVect(i)=eIOBP;
}
return eVect;
}
GridArray WWM_ReadGridFile_gr3(std::string const& GridFile)
{
GridArray GrdArr;
GrdArr.ModelName="WWM";
GrdArr.IsFE=1;
GrdArr.IsSpherical=true;
GrdArr.L_IndexSelect=false;
if (IsExistingFile(GridFile) == false) {
std::cerr << "Error in WWM_ReadGridFile_gr3\n";
std::cerr << "GridFile = " << GridFile << "\n";
std::cerr << "is missing\n";
throw TerminalException{1};
}
std::ifstream IN(GridFile);
std::string line;
std::getline(IN, line);
std::cerr << "line=" << line << "\n";
int mne, mnp;
IN >> mne;
IN >> mnp;
std::cerr << "mne=" << mne << " mnp=" << mnp << "\n";
GrdArr.INE=MyMatrix<int>(mne,3);
GrdArr.GrdArrRho.LON=MyMatrix<double>(mnp,1);
GrdArr.GrdArrRho.LAT=MyMatrix<double>(mnp,1);
GrdArr.GrdArrRho.DEP=MyMatrix<double>(mnp,1);
GrdArr.GrdArrRho.ANG=MyMatrix<double>(mnp,1);
GrdArr.GrdArrRho.MSK=MyMatrix<int>(mnp,1);
for (int iP=0; iP<mnp; iP++) {
int KTMP;
double XPDTMP, YPDTMP, ZPDTMP;
IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
GrdArr.GrdArrRho.LON(iP,0)=XPDTMP;
GrdArr.GrdArrRho.LAT(iP,0)=YPDTMP;
GrdArr.GrdArrRho.DEP(iP,0)=ZPDTMP;
GrdArr.GrdArrRho.ANG(iP,0)=0;
GrdArr.GrdArrRho.MSK(iP,0)=1;
}
for (int iE=0; iE<mne; iE++) {
int KTMP, LTMP, ip1, ip2, ip3;
IN >> KTMP >> LTMP >> ip1 >> ip2 >> ip3;
GrdArr.INE(iE,0)=ip1 - 1;
GrdArr.INE(iE,1)=ip2 - 1;
GrdArr.INE(iE,2)=ip3 - 1;
}
return GrdArr;
}
MyVector<int> WWM_ReadBoundFile_xfn(std::string const& BoundFile)
{
if (IsExistingFile(BoundFile) == false) {
std::cerr << "Error in WWM_ReadBoundFile_xfn\n";
std::cerr << "Missing BoundFile=" << BoundFile << "\n";
throw TerminalException{1};
}
std::ifstream IN(BoundFile);
std::string line;
for (int i=0; i<2; i++)
std::getline(IN, line);
int ITMP, JTMP;
IN >> ITMP;
std::getline(IN, line);
IN >> JTMP;
int mnp=ITMP + JTMP;
for (int i=0; i<7; i++)
std::getline(IN, line);
MyVector<int> eVect(mnp);
for (int i=0; i<mnp; i++) {
int KTMP;
double XPDTMP, YPDTMP, ZPDTMP;
IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
if (KTMP != i+1) {
std::cerr << "Inconsistency error\n";
throw TerminalException{1};
}
int eIOBP=int(ZPDTMP);
eVect(i)=eIOBP;
}
return eVect;
}
GridArray WWM_ReadGridFile_xfn(std::string const& GridFile)
{
GridArray GrdArr;
GrdArr.ModelName="WWM";
GrdArr.IsFE=1;
GrdArr.IsSpherical=true;
GrdArr.L_IndexSelect=false;
if (IsExistingFile(GridFile) == false) {
std::cerr << "Error in WWM_ReadGridFile_xfn\n";
std::cerr << "GridFile = " << GridFile << "\n";
std::cerr << "is missing\n";
throw TerminalException{1};
}
std::string line;
int ITMP, JTMP;
std::ifstream IN(GridFile);
for (int i=0; i<2; i++)
std::getline(IN, line);
IN >> ITMP;
std::getline(IN, line);
IN >> JTMP;
int mnp=ITMP + JTMP;
for (int i=0; i<7; i++)
std::getline(IN, line);
GrdArr.GrdArrRho.LON=MyMatrix<double>(mnp,1);
GrdArr.GrdArrRho.LAT=MyMatrix<double>(mnp,1);
GrdArr.GrdArrRho.DEP=MyMatrix<double>(mnp,1);
GrdArr.GrdArrRho.ANG=MyMatrix<double>(mnp,1);
GrdArr.GrdArrRho.MSK=MyMatrix<int>(mnp,1);
for (int iP=0; iP<mnp; iP++) {
int KTMP;
double XPDTMP, YPDTMP, ZPDTMP;
IN >> KTMP >> XPDTMP >> YPDTMP >> ZPDTMP;
GrdArr.GrdArrRho.LON(iP,0)=XPDTMP;
GrdArr.GrdArrRho.LAT(iP,0)=YPDTMP;
GrdArr.GrdArrRho.DEP(iP,0)=ZPDTMP;
GrdArr.GrdArrRho.ANG(iP,0)=0;
GrdArr.GrdArrRho.MSK(iP,0)=1;
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
GrdArr.INE(iE,0)=ip1 - 1;
GrdArr.INE(iE,1)=ip2 - 1;
GrdArr.INE(iE,2)=ip3 - 1;
}
return GrdArr;
}
GridArray NC_ReadWW3_GridFile(std::string const& eFile)
{
GridArray GrdArr;
GrdArr.ModelName="WWM";
GrdArr.IsFE=1;
GrdArr.IsSpherical=true;
GrdArr.L_IndexSelect=false;
GrdArr.INE=NC_ReadElements(eFile, "tri");
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
GrdArr.GrdArrRho.LON=LONarr;
GrdArr.GrdArrRho.LAT=LATarr;
GrdArr.GrdArrRho.DEP=DEParr;
GrdArr.GrdArrRho.ANG=ANGarr;
GrdArr.GrdArrRho.MSK=MSKarr;
GrdArr.ModelName="WW3";
return GrdArr;
}
GridArray TRIVIAL_GRID_ARRAY(QuadArray const& eQuad, int const& nbSplitLon, int const& nbSplitLat)
{
double MinLon = eQuad.MinLon;
double MinLat = eQuad.MinLat;
double MaxLon = eQuad.MaxLon;
double MaxLat = eQuad.MaxLat;
double deltaLon=(MaxLon - MinLon)/double(nbSplitLon-1);
double deltaLat=(MaxLat - MinLat)/double(nbSplitLat-1);
MyMatrix<double> LON(nbSplitLon, nbSplitLat);
MyMatrix<double> LAT(nbSplitLon, nbSplitLat);
MyMatrix<int> MSK(nbSplitLon, nbSplitLat);
MyMatrix<double> DEP(nbSplitLon, nbSplitLat);
MyMatrix<double> ANG(nbSplitLon, nbSplitLat);
for (int iLon=0; iLon<nbSplitLon; iLon++)
for (int iLat=0; iLat<nbSplitLat; iLat++) {
LON(iLon, iLat)=MinLon + deltaLon*iLon;
LAT(iLon, iLat)=MinLat + deltaLat*iLat;
MSK(iLon, iLat)=1;
DEP(iLon, iLat)=0;
ANG(iLon, iLat)=0;
}
CoordGridArrayFD GrdArrRho;
GrdArrRho.eta=nbSplitLon;
GrdArrRho.xi=nbSplitLat;
GrdArrRho.LON=LON;
GrdArrRho.LAT=LAT;
GrdArrRho.MSK=MSK;
GrdArrRho.DEP=DEP;
GrdArrRho.ANG=ANG;
GridArray GrdArr;
GrdArr.ModelName="trivial";
GrdArr.IsFE=0;
GrdArr.IsSpherical=true;
GrdArr.GrdArrRho=GrdArrRho;
return GrdArr;
}
void CutWorldMap(GridArray & GrdArr)
{
double eps=1e-8;
int nbPoint=GrdArr.GrdArrRho.LON.rows();
double LonSplit=0;
while(1) {
double MinDist=2400;
for (int iPoint=0; iPoint<nbPoint; iPoint++) {
double eLon=GrdArr.GrdArrRho.LON(iPoint,0);
if (eLon > 0)
eLon -= 360;
LonSplit= - 180 - eps;
double dist=fabs(eLon - LonSplit);
if (dist < MinDist)
MinDist=dist;
}
std::cerr << "eps=" << eps << " MinDist=" << MinDist << "\n";
if (MinDist > eps/2)
break;
eps *= 2;
}
int nbTrig=GrdArr.INE.rows();
std::vector<int> ListStatus(nbTrig);
int SumStatus=0;
for (int iTrig=0; iTrig<nbTrig; iTrig++) {
int i1=GrdArr.INE(iTrig,0);
int i2=GrdArr.INE(iTrig,1);
int i3=GrdArr.INE(iTrig,2);
double eLon1=GrdArr.GrdArrRho.LON(i1,0);
double eLon2=GrdArr.GrdArrRho.LON(i2,0);
double eLon3=GrdArr.GrdArrRho.LON(i3,0);
eLon1 -= LonSplit;
eLon2 -= LonSplit;
eLon3 -= LonSplit;
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
if (eProd12 < 0 || eProd23 < 0 || eProd31 < 0)
eStatus=0;
}
ListStatus[iTrig]=eStatus;
SumStatus += eStatus;
}
std::cerr << "SumStatus = " << SumStatus << "   nbTrig = " << nbTrig << "\n";
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
void CUT_HigherLatitude(GridArray & GrdArr, double MinLatCut, double MaxLatCut)
{
int mnp=GrdArr.GrdArrRho.LON.rows();
int mne=GrdArr.INE.rows();
std::vector<int> ListStatus(mnp);
std::vector<int> Index(mnp);
std::vector<int> RevIndex(mnp);
int iNodeNew=0;
std::vector<int> I_IndexSelectOld;
if (GrdArr.L_IndexSelect) {
I_IndexSelectOld=GrdArr.I_IndexSelect;
}
else {
for (int i=0; i<mnp; i++)
I_IndexSelectOld.push_back(i);
}
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
int nbNodeNew=iNodeNew;
MyMatrix<double> LONnew(nbNodeNew,1);
MyMatrix<double> LATnew(nbNodeNew,1);
MyMatrix<double> DEPnew(nbNodeNew,1);
MyMatrix<double> ANGnew(nbNodeNew,1);
MyMatrix<int> MSKnew(nbNodeNew,1);
for (int iNodeNewB=0; iNodeNewB<nbNodeNew; iNodeNewB++) {
int iNode=RevIndex[iNodeNewB];
double eLon=GrdArr.GrdArrRho.LON(iNode,0);
double eLat=GrdArr.GrdArrRho.LAT(iNode,0);
double eDep=GrdArr.GrdArrRho.DEP(iNode,0);
double eAng=GrdArr.GrdArrRho.ANG(iNode,0);
int eMsk=GrdArr.GrdArrRho.MSK(iNode,0);
LONnew(iNodeNewB,0)=eLon;
LATnew(iNodeNewB,0)=eLat;
DEPnew(iNodeNewB,0)=eDep;
ANGnew(iNodeNewB,0)=eAng;
MSKnew(iNodeNewB,0)=eMsk;
}
GrdArr.GrdArrRho.LON=LONnew;
GrdArr.GrdArrRho.LAT=LATnew;
GrdArr.GrdArrRho.DEP=DEPnew;
GrdArr.GrdArrRho.ANG=ANGnew;
GrdArr.GrdArrRho.MSK=MSKnew;
int nbTrigNew=0;
for (int iTrig=0; iTrig<mne; iTrig++) {
int i1=GrdArr.INE(iTrig,0);
int i2=GrdArr.INE(iTrig,1);
int i3=GrdArr.INE(iTrig,2);
if (ListStatus[i1] == 1 && ListStatus[i2] == 1 && ListStatus[i3] == 1)
nbTrigNew++;
}
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
GrdArr.INE=INEnew;
GrdArr.L_IndexSelect = true;
GrdArr.I_IndexSelect = I_IndexSelect;
}
double GetGridSpacing(GridArray const& GrdArr)
{
int IsFE=GrdArr.IsFE;
double SumDistKM=0;
int SumNb=0;
if (IsFE == 1) {
int nbEle=GrdArr.INE.rows();
for (int iEle=0; iEle<nbEle; iEle++)
for (int i=0; i<3; i++) {
int j=NextIdx(3,i);
int iNode1=GrdArr.INE(iEle,i);
int iNode2=GrdArr.INE(iEle,j);
double eLon1=GrdArr.GrdArrRho.LON(iNode1,0);
double eLat1=GrdArr.GrdArrRho.LAT(iNode1,0);
double eLon2=GrdArr.GrdArrRho.LON(iNode2,0);
double eLat2=GrdArr.GrdArrRho.LAT(iNode2,0);
double DistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
SumDistKM += DistKM;
SumNb += 1;
}
}
else {
int nbRow=GrdArr.GrdArrRho.LON.rows();
int nbCol=GrdArr.GrdArrRho.LON.cols();
for (int iRow=0; iRow<nbRow; iRow++) {
for (int iCol=1; iCol<nbCol; iCol++) {
double eLon1=GrdArr.GrdArrRho.LON(iRow, iCol);
double eLat1=GrdArr.GrdArrRho.LAT(iRow, iCol);
double eLon2=GrdArr.GrdArrRho.LON(iRow, iCol-1);
double eLat2=GrdArr.GrdArrRho.LAT(iRow, iCol-1);
double DistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
SumDistKM += DistKM;
SumNb += 1;
}
}
for (int iRow=1; iRow<nbRow; iRow++) {
for (int iCol=0; iCol<nbCol; iCol++) {
double eLon1=GrdArr.GrdArrRho.LON(iRow, iCol);
double eLat1=GrdArr.GrdArrRho.LAT(iRow, iCol);
double eLon2=GrdArr.GrdArrRho.LON(iRow-1, iCol);
double eLat2=GrdArr.GrdArrRho.LAT(iRow-1, iCol);
double DistKM=GeodesicDistanceKM(eLon1, eLat1, eLon2, eLat2);
SumDistKM += DistKM;
SumNb++;
}
}
}
double avgDistKM=SumDistKM / double(SumNb);
return avgDistKM;
}
ArrayHistory Sequential_ReadArrayHistory(std::string const& HisPrefix)
{
int len=HisPrefix.length();
std::vector<int> ListPos;
for (int iChar=0; iChar<len; iChar++) {
std::string eChar=HisPrefix.substr(iChar,1);
if (eChar == "/")
ListPos.push_back(iChar);
}
if (ListPos.size() == 0) {
std::cerr << "We should use / in HisPrefix for ROMS_IVICA\n";
throw TerminalException{1};
}
int iCharLast=ListPos[ListPos.size() - 1];
std::string eDir=HisPrefix.substr(0,iCharLast+1);
std::string RawPrefix=HisPrefix.substr(iCharLast+1,len-iCharLast-1);
std::cerr << "HisPrefix=" << HisPrefix << "\n";
std::cerr << "eDir=" << eDir << "\n";
std::cerr << "RawPrefix=" << RawPrefix << "\n";
std::vector<std::string> PreListFile=FILE_GetDirectoryListFile(eDir);
std::cerr << "|PreListFile|=" << PreListFile.size() << "\n";
std::vector<std::string> ListFileNames;
for (auto & eFile : PreListFile) {
std::string eFileTot=eDir + eFile;
std::vector<std::string> LStr=STRING_Split(eFileTot, RawPrefix);
int nbBlock=LStr.size();
if (nbBlock == 2)
ListFileNames.push_back(eFileTot);
}
std::cerr << "|ListFileNames|=" << ListFileNames.size() << "\n";
struct FullEntry {
int iFile;
int iTime;
double eTime;
};
std::vector<FullEntry> ListFull;
int nbFile=ListFileNames.size();
std::cerr << "nbFile=" << nbFile << "\n";
for (int iFile=0; iFile<nbFile; iFile++) {
std::string eFile=ListFileNames[iFile];
std::vector<double> LTime=NC_ReadTimeFromFile(eFile, "ocean_time");
int nbTime=LTime.size();
for (int iTime=0; iTime<nbTime; iTime++) {
double eTime=LTime[iTime];
FullEntry eFull{iFile, iTime, eTime};
ListFull.push_back(eFull);
}
}
std::cerr << "Now |ListFull|=" << ListFull.size() << "\n";
sort(ListFull.begin(), ListFull.end(),
[&](FullEntry const& a, FullEntry const& b) -> bool {
if (a.eTime < b.eTime)
return true;
return false;
});
int nbFull=ListFull.size();
std::cerr << "nbFull=" << nbFull << "\n";
std::vector<double> ListTime(nbFull);
std::vector<int> ListIFile(nbFull);
std::vector<int> ListIRec(nbFull);
for (int i=0; i<nbFull; i++) {
ListTime[i]=ListFull[i].eTime;
ListIFile[i]=ListFull[i].iFile;
ListIRec[i]=ListFull[i].iTime;
}
ArrayHistory eArr;
eArr.nbFile=nbFile;
eArr.nbTime=nbFull;
eArr.ListIFile=ListIFile;
eArr.ListIRec=ListIRec;
eArr.ListFileNames=ListFileNames;
eArr.ListTime=ListTime;
eArr.AppendVarName=false;
eArr.KindArchive="NETCDF";
eArr.TimeSteppingInfo="classic";
std::cerr << "Sequential array has been completed. Leaving\n";
return eArr;
}
ArrayHistory NC_ReadArrayHistory_Kernel(std::string const& HisPrefix, std::string const& StringTime, int const&nbDigit)
{
double FirstTime, LastTime;
std::vector<std::string> ListFileNames;
std::vector<int> ListIFile;
std::vector<int> ListIRec;
std::vector<double> ListTime;
ArrayHistory eArr;
if (IsExistingFile(HisPrefix)) {
std::cerr << "StringTime=" << StringTime << "\n";
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
double TheSep=GetListTimeSeparation(ListTime);
std::cerr << "TheSep=" << TheSep << "\n";
if (TheSep < 0) {
eArr.TimeSteppingInfo="classic";
}
else {
eArr.TimeSteppingInfo="singlefile";
eArr.SeparationTime=TheSep;
}
std::cerr << "eArr.SeparationTime = " << eArr.SeparationTime << "\n";
}
else {
int iFileBegin=0;
while(1) {
iFileBegin++;
std::string TheHisFile=HisPrefix + StringNumber(iFileBegin, nbDigit) + ".nc";
if (IsExistingFile(TheHisFile) == 1)
break;
if (iFileBegin == 9999) {
std::cerr << "maybe you specified wrong HisPrefix\n";
std::cerr << "HisPrefix = " << HisPrefix << "\n";
std::cerr << " there is no files  HisPrefix????.nc\n";
std::cerr << "Please correct\n";
throw TerminalException{1};
}
}
std::cerr << "iFileBegin=" << iFileBegin << "\n";
int iFileEnd=iFileBegin;
while(1) {
std::string TheHisFile=HisPrefix + StringNumber(iFileEnd+1, nbDigit) + ".nc";
if (IsExistingFile(TheHisFile) == 0)
break;
iFileEnd++;
}
std::cerr << "iFileEnd=" << iFileEnd << "\n";
std::string TheHisFileBegin=HisPrefix + StringNumber(iFileBegin, nbDigit) + ".nc";
std::vector<double> LTimeBegin=NC_ReadTimeFromFile(TheHisFileBegin, StringTime);
int nbRecBegin=LTimeBegin.size();
double DeltaTime;
if (nbRecBegin > 1) {
DeltaTime=LTimeBegin[1]-LTimeBegin[0];
}
else {
if (iFileEnd > iFileBegin) {
std::string TheHisFile=HisPrefix + StringNumber(iFileBegin+1, nbDigit) + ".nc";
std::vector<double> LTimeBP1=NC_ReadTimeFromFile(TheHisFile, StringTime);
DeltaTime=LTimeBP1[0] - LTimeBegin[0];
}
else {
DeltaTime=0;
}
}
int iFileMiddle, nbRecMiddle;
if (iFileEnd != iFileBegin) {
iFileMiddle=iFileBegin+1;
std::string TheHisFile=HisPrefix + StringNumber(iFileMiddle, nbDigit) + ".nc";
std::vector<double> LTimeMiddle=NC_ReadTimeFromFile(TheHisFile, StringTime);
nbRecMiddle=LTimeMiddle.size();
}
else {
iFileMiddle=iFileBegin;
nbRecMiddle=nbRecBegin;
}
std::cerr << "iFileMiddle=" << iFileMiddle << "\n";
int nbRecEnd;
if (iFileEnd != iFileMiddle) {
std::string TheHisFile=HisPrefix + StringNumber(iFileEnd, nbDigit) + ".nc";
std::vector<double> LTimeEnd=NC_ReadTimeFromFile(TheHisFile, StringTime);
nbRecEnd=LTimeEnd.size();
}
else {
nbRecEnd=nbRecMiddle;
}
int nbFile=1 + iFileEnd - iFileBegin;
std::vector<int> NbPerArray(nbFile);
for (int iFile=0; iFile<nbFile; iFile++)
NbPerArray[iFile]=nbRecMiddle;
NbPerArray[0]=nbRecBegin;
NbPerArray[nbFile-1]=nbRecEnd;
double currentTime=LTimeBegin[0];
FirstTime=currentTime;
for (int iFile=0; iFile<nbFile; iFile++) {
int iFileTot=iFile+iFileBegin;
std::string TheHisFile=HisPrefix + StringNumber(iFileTot, nbDigit) + ".nc";
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
eArr.TimeSteppingInfo="multiplenetcdf";
eArr.SeparationTime=DeltaTime;
eArr.nbRecBegin=nbRecBegin;
eArr.nbRecMiddle=nbRecMiddle;
eArr.HisPrefix=HisPrefix;
}
int nbFile=ListFileNames.size();
std::cerr << "nbFile=" << nbFile << "\n";
eArr.nbTime=ListTime.size();
eArr.FirstTime=FirstTime;
eArr.LastTime=LastTime;
if (eArr.TimeSteppingInfo != "multiplenetcdf") {
eArr.nbFile=nbFile;
eArr.ListFileNames=ListFileNames;
eArr.ListIFile=ListIFile;
eArr.ListIRec=ListIRec;
eArr.ListTime=ListTime;
}
eArr.nbDigit=nbDigit;
eArr.AppendVarName=false;
eArr.KindArchive="NETCDF";
std::cerr << "Leaving NC_ReadArrayHistory_Kernel\n";
return eArr;
}
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
eArr.TimeSteppingInfo="classic";
return eArr;
}
struct TripleModelDesc {
std::string ModelName;
std::string GridFile;
std::string BoundFile;
std::string HisPrefix;
bool CutWorldMap;
bool HigherLatitudeCut;
double MinLatCut;
double MaxLatCut;
};
std::string GET_GRID_FILE(TripleModelDesc const& eTriple)
{
std::string eModelName=eTriple.ModelName;
std::string HisPrefix=eTriple.HisPrefix;
if (eModelName == "COSMO")
return HisPrefix + "0001.nc";
if (eModelName == "WAM")
return HisPrefix + "0001.nc";
if (eModelName == "ROMS" || eModelName == "ROMS_IVICA")
return eTriple.GridFile;
if (eModelName == "WWM" || eModelName == "WWM_DAILY")
return eTriple.GridFile;
if (eModelName == "SCHISM_SFLUX")
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
throw TerminalException{1};
}
return ListFile[0];
}
if (eModelName == "GRIB_WAM_FORT30") {
std::string eFile=HisPrefix;
if (IsExistingFile(eFile) == false) {
std::cerr << "The file eFile = " << eFile << " is missing\n";
std::cerr << "It serves as grid and should be put in HisPrefix\n";
throw TerminalException{1};
}
return eFile;
}
std::cerr << "Error in GET_GRID_FILE\n";
std::cerr << "Did not find the matching model for the grid\n";
std::cerr << "Please correct\n";
throw TerminalException{1};
}
GridArray ReadUnstructuredGrid(std::string const& GridFile, std::string const& BoundFile)
{
std::string eExtension=FILE_GetExtension(GridFile);
if (eExtension == "gr3" || eExtension == "ll") {
GridArray GrdArr=WWM_ReadGridFile_gr3(GridFile);
if (BoundFile != "unset") {
std::cerr << "BoundFile=" << BoundFile << "\n";
MyVector<int> eVect=WWM_ReadBoundFile_gr3(BoundFile);
if (eVect.size() != GrdArr.GrdArrRho.LON.size()) {
std::cerr << "not same number of vertices between grid file and boundary file\n";
std::cerr << "nbVert(grid)=" << GrdArr.GrdArrRho.LON.size() << "\n";
std::cerr << "nbVert(bound)=" << eVect.size() << "\n";
GrdArr.IOBP=eVect;
}
}
return GrdArr;
}
if (eExtension == "dat") {
GridArray GrdArr=WWM_ReadGridFile_xfn(GridFile);
if (BoundFile != "unset") {
MyVector<int> eVect=WWM_ReadBoundFile_xfn(BoundFile);
if (eVect.size() != GrdArr.GrdArrRho.LON.size()) {
std::cerr << "not same number of vertices between grid file and boundary file\n";
std::cerr << "nbVert(grid)=" << GrdArr.GrdArrRho.LON.size() << "\n";
std::cerr << "nbVert(bound)=" << eVect.size() << "\n";
GrdArr.IOBP=eVect;
}
}
return GrdArr;
}
if (eExtension == "nc")
return WWM_ReadGridFile_netcdf(GridFile);
std::cerr << "Error in reading grid for WWM\n";
std::cerr << "We did not find the right kind\n";
throw TerminalException{1};
}
GridArray PRE_RETRIEVE_GRID_ARRAY(TripleModelDesc const& eTriple)
{
std::string eModelName=eTriple.ModelName;
CHECK_Model_Allowedness(eModelName);
std::string GridFile=GET_GRID_FILE(eTriple);
if (eModelName == "COSMO") {
return NC_ReadCosmoWamStructGridFile(GridFile, "atm");
}
if (eModelName == "WAM") {
return NC_ReadWamGridFile(GridFile);
}
if (eModelName == "ROMS" || eModelName == "ROMS_IVICA") {
return NC_ReadRomsGridFile(GridFile);
}
if (eModelName == "WWM" || eModelName == "WWM_DAILY") {
std::string BoundFile=eTriple.BoundFile;
return ReadUnstructuredGrid(GridFile, BoundFile);
}
if (eModelName == "WW3") {
return NC_ReadWW3_GridFile(GridFile);
}
if (eModelName == "SCHISM_SFLUX") {
return NC_ReadSCHISM_sflux_grid(GridFile);
}
if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO") {
std::string HisPrefix=eTriple.HisPrefix;
std::vector<std::string> ListFile=GRIB_GetAllFilesInDirectory(HisPrefix);
if (ListFile.size() == 0) {
std::cerr << "The list of files is empty\n";
std::cerr << "Error happened in GRIB_GetAllFilesInDirectory\n";
throw TerminalException{1};
}
std::string eFileName=ListFile[0];
return GRIB_ReadGridArray(eFileName, eModelName);
}
if (eModelName == "GRIB_WAM_FORT30") {
std::string eFileName=eTriple.HisPrefix;
if (IsExistingFile(eFileName) == false) {
std::cerr << "The file eFileName = " << eFileName << " is missing\n";
std::cerr << "This is set by HisPRefix and serves for the data storage\n";
throw TerminalException{1};
}
return GRIB_ReadGridArray(eFileName, eModelName);
}
std::cerr << "Error in PRE_RETRIEVE_GRID_ARRAY\n";
std::cerr << "Did not find the matching model for the grid\n";
std::cerr << "Please correct\n";
throw TerminalException{1};
}
GridArray RETRIEVE_GRID_ARRAY(TripleModelDesc const& eTriple)
{
GridArray GrdArr=PRE_RETRIEVE_GRID_ARRAY(eTriple);
if (GrdArr.IsFE == 0)
return GrdArr;
if (eTriple.CutWorldMap)
CutWorldMap(GrdArr);
if (eTriple.HigherLatitudeCut) {
double MinLatCut=eTriple.MinLatCut;
double MaxLatCut=eTriple.MaxLatCut;
CUT_HigherLatitude(GrdArr, MinLatCut, MaxLatCut);
}
std::cerr << "After CUT_HigherLatitude\n";
CHECK_UnstructuredGrid(GrdArr);
CHECK_CombinatorialGrid(GrdArr);
CHECK_COORDINATE_ORIENTATION(GrdArr);
std::cerr << "Before returning GrdArr in RETRIEVE_GRID_ARRAY\n";
return GrdArr;
}
ArrayHistory NC_ReadArrayHistory(TripleModelDesc const& eTriple)
{
std::string StringTime="ocean_time";
std::string eModelName=eTriple.ModelName;
std::string HisPrefix=eTriple.HisPrefix;
if (eModelName == "WW3") {
std::string HisFile=GET_GRID_FILE(eTriple);
return WW3_ReadArrayHistory(HisFile, HisPrefix);
}
if (eModelName == "ROMS_IVICA" || eModelName == "WWM_DAILY")
return Sequential_ReadArrayHistory(HisPrefix);
if (eModelName == "SCHISM_SFLUX")
return NC_ReadArrayHistory_Kernel(HisPrefix, "time", 3);
return NC_ReadArrayHistory_Kernel(HisPrefix, StringTime, 4);
}
ArrayHistory GRIB_ReadArrayHistory(std::string const& HisPrefix)
{
std::vector<std::string> ListFile;
if (IsExistingFile(HisPrefix) == true && FILE_IsRegularFile(HisPrefix) == true) {
ListFile = {HisPrefix};
}
else {
ListFile = GRIB_GetAllFilesInDirectory(HisPrefix);
}
int nbFile=ListFile.size();
double MaxErrorTime=0.01;
std::vector<GRIB_MessageInfo> ListAllMessage;
for (auto & eFile : ListFile) {
std::vector<GRIB_MessageInfo> ListMessage=GRIB_GetAllListPairTime(eFile);
int nbMessage=ListMessage.size();
if (nbMessage == 0) {
std::cerr << "Remark: eFile = " << eFile << " has zero messages\n";
}
for (auto & eMesg : ListMessage)
ListAllMessage.push_back(eMesg);
}
int TotalNbMessage=ListAllMessage.size();
if (TotalNbMessage == 0) {
std::cerr << "TotalNbMessage=" << TotalNbMessage << "\n";
std::cerr << "|ListFile|=" << ListFile.size() << "\n";
std::cerr << "We have zero messages. No work can be done\n";
throw TerminalException{1};
}
sort(ListAllMessage.begin(), ListAllMessage.end(),
[&](GRIB_MessageInfo const& a, GRIB_MessageInfo const& b) -> bool {
if (a.time < b.time)
return true;
return false;
});
double TimePrev=ListAllMessage[0].time;
std::vector<double> ListTime{TimePrev};
std::vector<int> ListITime(TotalNbMessage);
int posTime=0;
for (int iMesg=0; iMesg<TotalNbMessage; iMesg++) {
GRIB_MessageInfo eMesg=ListAllMessage[iMesg];
double eTime=eMesg.time;
double TimeDiff = fabs(eTime - TimePrev);
if (TimeDiff > MaxErrorTime) {
ListTime.push_back(eTime);
TimePrev=eTime;
posTime++;
}
ListITime[iMesg]=posTime;
}
std::vector<double> ListStartTime;
std::vector<int> ListIStartTime(TotalNbMessage);
double tolDay=double(1)/double(10000);
auto GetStartTime=[&](double const& timeStart) -> int {
int nbTimeStart=ListStartTime.size();
for (int iTimeStart=0; iTimeStart<nbTimeStart; iTimeStart++)
if (fabs(timeStart - ListStartTime[iTimeStart]) < tolDay)
return iTimeStart;
ListStartTime.push_back(timeStart);
return nbTimeStart;
};
for (int iMesg=0; iMesg<TotalNbMessage; iMesg++) {
double eStartTime=ListAllMessage[iMesg].timeStart;
int iTimeStart=GetStartTime(eStartTime);
ListIStartTime[iMesg]=iTimeStart;
}
int nbTime=ListTime.size();
std::vector<std::vector<GRIB_MessageInfo> > ListListMessages(nbTime);
for (int iMesg=0; iMesg<TotalNbMessage; iMesg++) {
GRIB_MessageInfo eMesg=ListAllMessage[iMesg];
int iTime=ListITime[iMesg];
if (iTime < 0 || iTime >= nbTime) {
std::cerr << "iTime=" << iTime << " but nbTime=" << nbTime << "\n";
throw TerminalException{1};
}
ListListMessages[iTime].push_back(eMesg);
}
std::set<std::string> SetRawNames;
for (int iTime=0; iTime<nbTime; iTime++) {
std::vector<GRIB_MessageInfo> ListMessages=ListListMessages[iTime];
std::set<std::string> ListShortName;
for (auto & eMesg : ListMessages) {
ListShortName.insert(eMesg.shortName);
SetRawNames.insert(eMesg.shortName);
}
std::vector<GRIB_MessageInfo> NewListMessages;
for (auto & eShortName : ListShortName) {
bool IsFirst=true;
GRIB_MessageInfo NewMesg;
for (auto & eMesg : ListMessages) {
if (eMesg.shortName == eShortName) {
if (IsFirst) {
NewMesg = eMesg;
IsFirst=false;
}
else {
if (eMesg.timeStart < NewMesg.timeStart)
NewMesg = eMesg;
}
}
}
NewListMessages.push_back(NewMesg);
}
ListListMessages[iTime] = NewListMessages;
}
std::vector<std::string> RawVarNames;
for (auto& eName : SetRawNames)
RawVarNames.push_back(eName);
ArrayHistory eArr;
eArr.nbFile=nbFile;
eArr.nbTime=nbTime;
eArr.ListListMessages=ListListMessages;
eArr.ListAllMessage=ListAllMessage;
eArr.ListStartTime=ListStartTime;
eArr.ListIStartTime=ListIStartTime;
eArr.RawVarNames=RawVarNames;
eArr.ListITime=ListITime;
eArr.ListTime=ListTime;
eArr.KindArchive="GRIB";
eArr.TimeSteppingInfo="classic";
return eArr;
}
ArrayHistory ReadArrayHistory(TripleModelDesc const& eTriple)
{
ArrayHistory eArr;
std::string HisPrefix=eTriple.HisPrefix;
std::string eModelName=eTriple.ModelName;
CHECK_Model_Allowedness(eModelName);
if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_COSMO" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_WAM_FORT30") {
eArr=GRIB_ReadArrayHistory(HisPrefix);
}
else {
eArr=NC_ReadArrayHistory(eTriple);
}
return eArr;
}
bool TOTALARR_IsVar(TotalArrGetData const& TotalArr, std::string const& eVar)
{
if (TotalArr.eArr.KindArchive == "NETCDF") {
std::string HisFile;
int iFile=0;
ArrayHistory eArr=TotalArr.eArr;
if (eArr.TimeSteppingInfo == "multiplenetcdf") {
HisFile=eArr.HisPrefix + StringNumber(iFile+1,4) + ".nc";
}
else {
if (eArr.AppendVarName) {
HisFile=eArr.ListFileNames[iFile] + eVar + ".nc";;
}
else {
HisFile=eArr.ListFileNames[iFile];
}
}
return NC_IsVar(HisFile, eVar);
}
if (TotalArr.eArr.KindArchive == "GRIB") {
for (auto & eVarName : TotalArr.eArr.RawVarNames)
if (eVarName == eVar)
return true;
return false;
}
std::cerr << "Error in TOTALARR_IsVar\n";
std::cerr << "The KindArchive does not allow to find the nature\n";
std::cerr << "KindArchive=" << TotalArr.eArr.KindArchive << "\n";
throw TerminalException{1};
}
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
throw TerminalException{1};
}
PairMinMax ComputeMinMax(GridArray const& GrdArr, MyMatrix<double> const& F)
{
bool IsFirst=true;
int eta_rho=F.rows();
int xi_rho=F.cols();
int eta_rho_msk=GrdArr.GrdArrRho.MSK.rows();
int xi_rho_msk=GrdArr.GrdArrRho.MSK.cols();
if (eta_rho != eta_rho_msk || xi_rho != xi_rho_msk) {
std::cerr << "ComputeMinMax error : Inconsistency in dimension\n";
std::cerr << "  F: eta_rho=" << eta_rho << " xi_rho=" << xi_rho << "\n";
std::cerr << "MSK: eta_rho=" << eta_rho_msk << " xi_rho=" << xi_rho_msk << "\n";
throw TerminalException{1};
}
double TheMin=0;
double TheMax=0;
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
int eMsk;
if (GrdArr.IsFE == 1) {
eMsk=1;
}
else {
eMsk=GrdArr.GrdArrRho.MSK(i,j);
}
if (eMsk == 1) {
double eVal=F(i,j);
if (IsFirst == true) {
TheMin=eVal;
TheMax=eVal;
}
else {
if (eVal < TheMin)
TheMin=eVal;
if (eVal > TheMax)
TheMax=eVal;
}
IsFirst=false;
}
}
return {TheMin, TheMax};
}
void ApplyPlotBound(TotalArrGetData const& TotalArr, RecVar & eRecVar, std::string const& eVarName, PlotBound const& ePlotBound)
{
int nbSingle=ePlotBound.BoundSingle_var.size();
int nbSingleMin=ePlotBound.BoundSingle_min.size();
int nbSingleMax=ePlotBound.BoundSingle_max.size();
if (nbSingle != nbSingleMin || nbSingle != nbSingleMax) {
std::cerr << "Number of entries in BoundSingle_var, BoundSingle_min, BoundSingle_max\n";
std::cerr << "Should all be the same. Now,\n";
std::cerr << "nbSingle    = " << nbSingle << "\n";
std::cerr << "nbSingleMin = " << nbSingleMin << "\n";
std::cerr << "nbSingleMax = " << nbSingleMax << "\n";
}
for (int iS=0; iS<nbSingle; iS++)
if (ePlotBound.BoundSingle_var[iS] == eVarName) {
eRecVar.RecS.minval=ePlotBound.BoundSingle_min[iS];
eRecVar.RecS.maxval=ePlotBound.BoundSingle_max[iS];
}
int nbDiff=ePlotBound.BoundDiff_var.size();
int nbDiffMin=ePlotBound.BoundDiff_min.size();
int nbDiffMax=ePlotBound.BoundDiff_max.size();
if (nbDiff != nbDiffMin || nbDiff != nbDiffMax) {
std::cerr << "Number of entries in BoundDiff_var, BoundDiff_min, BoundDiff_max\n";
std::cerr << "Should all be the same. Now,\n";
std::cerr << "nbDiff    = " << nbDiff << "\n";
std::cerr << "nbDiffMin = " << nbDiffMin << "\n";
std::cerr << "nbDiffMax = " << nbDiffMax << "\n";
}
for (int iD=0; iD<nbDiff; iD++)
if (ePlotBound.BoundDiff_var[iD] == eVarName) {
eRecVar.RecS.mindiff=ePlotBound.BoundDiff_min[iD];
eRecVar.RecS.maxdiff=ePlotBound.BoundDiff_max[iD];
}
int eSize=eRecVar.F.size();
if (ePlotBound.VariableRange == true && eSize > 0) {
PairMinMax ePair=ComputeMinMax(TotalArr.GrdArr, eRecVar.F);
eRecVar.RecS.mindiff=ePair.TheMin;
eRecVar.RecS.maxdiff=ePair.TheMax;
eRecVar.RecS.minval=ePair.TheMin;
eRecVar.RecS.maxval=ePair.TheMax;
}
}
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
throw TerminalException{1};
}
strAll=StringNumber(iTime, 4) + "_" + strFile;
}
if (eQuery.NatureQuery != "instant")
strAll += "_" + eQuery.NatureQuery;
return strAll;
}
std::string GetStrPresOfPlot(VarQuery const& eQuery)
{
std::string strPres1=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay);
if (eQuery.NatureQuery == "instant") {
return "at " + strPres1;
}
if (eQuery.NatureQuery == "swathMax") {
double TimeFrameDay=eQuery.TimeFrameDay;
std::string strPres2=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay + TimeFrameDay);
return "max from " + strPres1 + " to " + strPres2;
}
if (eQuery.NatureQuery == "swathMin") {
double TimeFrameDay=eQuery.TimeFrameDay;
std::string strPres2=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay + TimeFrameDay);
return "min from " + strPres1 + " to " + strPres2;
}
if (eQuery.NatureQuery == "average") {
double TimeFrameDay=eQuery.TimeFrameDay;
std::string strPres2=DATE_ConvertMjd2mystringPresReduced(eQuery.eTimeDay + TimeFrameDay);
return "avg. from " + strPres1 + " to " + strPres2;
}
std::cerr << "Failed to find NatureQuery in list of available options\n";
std::cerr << "eQuery.NatureQuery=" << eQuery.NatureQuery << "\n";
throw TerminalException{1};
}
RecVar GetTrivialArrayPlot(GridArray const& GrdArr)
{
int eta_rho=GrdArr.GrdArrRho.LON.rows();
int xi_rho =GrdArr.GrdArrRho.LON.cols();
MyMatrix<double> F(eta_rho, xi_rho);
double MaxLon=GrdArr.GrdArrRho.LON.maxCoeff();
double MinLon=GrdArr.GrdArrRho.LON.minCoeff();
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
double eLon=GrdArr.GrdArrRho.LON(i,j);
double eVal=(eLon - MinLon)/(MaxLon - MinLon);
F(i, j)=eVal;
}
RecVar eRecVar;
eRecVar.RecS.strAll="unset";
eRecVar.RecS.VarName1="Track";
eRecVar.RecS.VarName2="Track";
eRecVar.RecS.minval=2;
eRecVar.RecS.maxval=3;
eRecVar.RecS.Unit="nondim.";
eRecVar.F=F;
return eRecVar;
}
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
struct DrawScatterArr {
std::string VarNameAB_file;
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
double ePair[2];
ePair[0]=0;
ePair[1]=TheMaxX;
eVarListX.putVar(ePair);
ePair[0]=0;
ePair[1]=TheMaxY;
eVarListY.putVar(ePair);
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
MyMatrix<int> canvasInt(aSize,bSize);
for (int iA=0; iA<aSize; iA++)
for (int iB=0; iB<bSize; iB++)
canvasInt(iA,iB)=0;
for (int iEnt=0; iEnt<nb; iEnt++) {
double eA=eDrawScatter.eVectA(iEnt);
double eB=eDrawScatter.eVectB(iEnt);
double iA_d=(eA - data_rangeA[0])/deltaA;
double iB_d=(eB - data_rangeB[0])/deltaB;
int iA=int(floor(iA_d));
int iB=int(floor(iB_d));
if (iA >=0 && iA<aSize && iB>=0 && iB<bSize)
canvasInt(iA,iB)++;
}
double *canvas;
canvas=new double[aSize*bSize];
double MissVal=0;
idx=0;
for (int iA=0; iA<aSize; iA++)
for (int iB=0; iB<bSize; iB++) {
double eVal;
int eValI=canvasInt(iA,iB);
if (eValI == 0)
eVal=MissVal;
else
eVal=log10(double(eValI));
canvas[idx]=eVal;
idx++;
}
eVarCanvas.putVar(canvas);
delete [] canvas;
}
void DEFINE_QUIVER_NC(std::string const& eFileNC,
GridArray const& GrdArr,
MyMatrix<double> const& U_rho, MyMatrix<double> const& V_rho,
MyMatrix<double> const& F_rho)
{
int idx;
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
throw TerminalException{1};
}
if (xi != xi_lat || xi != xi_U || xi != xi_V) {
std::cerr << "All dimensions should be the same.\n";
std::cerr << "Now we have following:\n";
std::cerr << "xi_lon=" << xi << "\n";
std::cerr << "xi_lat=" << xi_lat << "\n";
std::cerr << "xi_U=" << xi_U << "\n";
std::cerr << "xi_V=" << xi_V << "\n";
throw TerminalException{1};
}
std::string typeName="double";
std::string typeNameInt="int";
std::string eEta="eta_rho";
std::string eXi ="xi_rho";
std::string Lon="lon";
std::string Lat="lat";
std::string Uvar="u";
std::string Vvar="v";
std::string Fvar="F";
if (xi > 1) {
bool ApplyCritValue=true;
double eCritValue = -10000;
double dataMiss[1];
dataMiss[0]=eCritValue;
std::string MissVal="_FillValue";
netCDF::NcDim eDimEta=dataFile.addDim(eEta, eta);
netCDF::NcDim eDimXi =dataFile.addDim(eXi, xi);
std::vector<std::string> ListDim={eEta, eXi};
netCDF::NcVar eVarLON=dataFile.addVar(Lon, typeName, ListDim);
netCDF::NcVar eVarLAT=dataFile.addVar(Lat, typeName, ListDim);
netCDF::NcVar eVarU=dataFile.addVar(Uvar, typeName, ListDim);
eVarU.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
netCDF::NcVar eVarV=dataFile.addVar(Vvar, typeName, ListDim);
eVarV.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
netCDF::NcVar eVarF=dataFile.addVar(Fvar, typeName, ListDim);
eVarF.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
double *valLON, *valLAT, *valU, *valV, *valF;
valLON=new double[eta*xi];
valLAT=new double[eta*xi];
valU =new double[eta*xi];
valV =new double[eta*xi];
valF =new double[eta*xi];
idx=0;
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++) {
valLON[idx]=GrdArr.GrdArrRho.LON(i, j);
valLAT[idx]=GrdArr.GrdArrRho.LAT(i, j);
if (GrdArr.GrdArrRho.MSK(i,j) == 1 || ApplyCritValue == false) {
valU[idx]=U_rho(i,j);
valV[idx]=V_rho(i,j);
valF[idx]=F_rho(i,j);
}
else {
valU[idx]=eCritValue;
valV[idx]=eCritValue;
valF[idx]=eCritValue;
}
idx++;
}
eVarLON.putVar(valLON);
eVarLAT.putVar(valLAT);
eVarU.putVar(valU);
eVarV.putVar(valV);
eVarF.putVar(valF);
delete [] valLON;
delete [] valLAT;
delete [] valU;
delete [] valV;
delete [] valF;
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
netCDF::NcVar eVarF=dataFile.addVar(Fvar, typeName, ListDim);
std::string Fine="ele";
std::vector<std::string> ListDimINE={eMne, eThree};
netCDF::NcVar eVarINE=dataFile.addVar(Fine, typeNameInt, ListDimINE);
double *valLON, *valLAT, *valU, *valV, *valF;
valLON=new double[mnp];
valLAT=new double[mnp];
valU =new double[mnp];
valV =new double[mnp];
valF =new double[mnp];
for (int i=0; i<eta; i++) {
valLON[i]=GrdArr.GrdArrRho.LON(i,0);
valLAT[i]=GrdArr.GrdArrRho.LAT(i,0);
valU[i]=U_rho(i,0);
valV[i]=V_rho(i,0);
valF[i]=F_rho(i,0);
}
eVarLON.putVar(valLON);
eVarLAT.putVar(valLAT);
eVarU.putVar(valU);
eVarV.putVar(valV);
eVarF.putVar(valF);
delete [] valLON;
delete [] valLAT;
delete [] valU;
delete [] valV;
delete [] valF;
int *valI;
valI=new int[3*mne];
idx=0;
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
void DEFINE_PCOLOR_NC_NCL(std::string const& eFileNC,
GridArray const& GrdArr,
MyMatrix<double> const& F_rho,
bool const& WriteDEP,
std::vector<SeqLineSegment> const& ListLineSegment)
{
int idx;
netCDF::NcFile dataFile(eFileNC, netCDF::NcFile::replace, netCDF::NcFile::nc4);
int eta=GrdArr.GrdArrRho.LON.rows();
int xi =GrdArr.GrdArrRho.LON.cols();
if (GrdArr.GrdArrRho.LAT.rows() != eta || GrdArr.GrdArrRho.LAT.cols() != xi) {
std::cerr << "Dimension errors\n";
std::cerr << "dim(LON)=" << eta << "/" << xi << "\n";
std::cerr << "dim(LAT)=" << GrdArr.GrdArrRho.LAT.rows() << "/" << GrdArr.GrdArrRho.LAT.cols() << "\n";
throw TerminalException{1};
}
if (F_rho.rows() != eta || F_rho.cols() != xi) {
std::cerr << "Dimension errors\n";
std::cerr << "dim(LON)=" << eta << "/" << xi << "\n";
std::cerr << "dim(F  )=" << F_rho.rows() << "/" << F_rho.cols() << "\n";
throw TerminalException{1};
}
if (GrdArr.IsFE == 0) {
bool ApplyCritValue=true;
double eCritValue = -10000;
double dataMiss[1];
dataMiss[0]=eCritValue;
std::string MissVal="_FillValue";
netCDF::NcDim eDimEta=dataFile.addDim("eta_rho", eta);
netCDF::NcDim eDimXi =dataFile.addDim("xi_rho", xi);
std::vector<std::string> ListDim={"eta_rho", "xi_rho"};
netCDF::NcVar eVarLON=dataFile.addVar("lon", "double", ListDim);
netCDF::NcVar eVarLAT=dataFile.addVar("lat", "double", ListDim);
netCDF::NcVar eVarF=dataFile.addVar("field", "double", ListDim);
eVarF.putAtt(MissVal, netCDF::NcType::nc_DOUBLE, 1, dataMiss);
double *valLON, *valLAT, *valF;
valLON=new double[eta*xi];
valLAT=new double[eta*xi];
valF =new double[eta*xi];
idx=0;
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++) {
valLON[idx]=GrdArr.GrdArrRho.LON(i, j);
valLAT[idx]=GrdArr.GrdArrRho.LAT(i, j);
double eValF;
if (GrdArr.GrdArrRho.MSK(i,j) == 1 || ApplyCritValue == false) {
eValF = F_rho(i,j);
}
else {
eValF = eCritValue;
}
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
idx=0;
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
int mnp=eta;
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
valF =new double[mnp];
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
int *valI;
valI=new int[3*mne];
idx=0;
for (int ie=0; ie<mne; ie++)
for (int i=0; i<3; i++) {
int eConn=GrdArr.INE(ie,i);
valI[idx]=eConn;
idx++;
}
eVarINE.putVar(valI);
delete [] valI;
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
idx=0;
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
throw TerminalException{1};
}
}
struct DrawLinesArr {
bool DoTitle;
std::string TitleStr;
std::string VarName;
AnnotationRec TheAnnot;
bool DoTitleLines;
double TheMax;
double TheMin;
bool IsTimeSeries;
bool PairComparison;
std::string XAxisString;
std::string YAxisString;
std::vector<std::string> ListName_plot;
MyVector<double> ListX;
std::vector<MyVector<double> > ListListVect;
};
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
double *valX;
valX=new double[nbEntry];
for (int i=0; i<nbEntry; i++) {
double eX=eDrawArr.ListX(i);
valX[i]=eX;
}
eVarX.putVar(valX);
delete [] valX;
}
FullNamelist NAMELIST_GetStandardALTIMETRY_COMPARISON()
{
std::map<std::string, SingleBlock> ListBlock;
std::string BlockName1="PROC";
std::map<std::string, int> ListIntValues1;
std::map<std::string, bool> ListBoolValues1;
std::map<std::string, double> ListDoubleValues1;
std::map<std::string, std::vector<double> > ListListDoubleValues1;
std::map<std::string, std::string> ListStringValues1;
std::map<std::string, std::vector<std::string> > ListListStringValues1;
ListStringValues1["MODELNAME"]="unset MODELNAME";
ListStringValues1["GridFile"]="unset GridFile";
ListStringValues1["BoundFile"]="unset";
ListStringValues1["HisPrefix"]="unset HisPrefix";
ListStringValues1["PicPrefix"]="unset PicPrefix";
ListStringValues1["Extension"]="png";
ListStringValues1["__NaturePlot"]="ALTIMETRY";
ListBoolValues1["FirstCleanDirectory"]=true;
ListBoolValues1["KeepNC_NCL"]=false;
ListIntValues1["NPROC"]=1;
SingleBlock BlockPROC;
BlockPROC.ListIntValues=ListIntValues1;
BlockPROC.ListBoolValues=ListBoolValues1;
BlockPROC.ListDoubleValues=ListDoubleValues1;
BlockPROC.ListListDoubleValues=ListListDoubleValues1;
BlockPROC.ListStringValues=ListStringValues1;
BlockPROC.ListListStringValues=ListListStringValues1;
BlockPROC.BlockName=BlockName1;
ListBlock["PROC"]=BlockPROC;
std::string BlockName2="SELECT";
std::map<std::string, int> ListIntValues2;
std::map<std::string, bool> ListBoolValues2;
std::map<std::string, double> ListDoubleValues2;
std::map<std::string, std::vector<double> > ListListDoubleValues2;
std::map<std::string, std::string> ListStringValues2;
std::map<std::string, std::vector<std::string> > ListListStringValues2;
ListIntValues2["GEOSELECTION"]=1;
ListDoubleValues2["MinLON"]=-7;
ListDoubleValues2["MaxLON"]=37;
ListDoubleValues2["MinLAT"]=30;
ListDoubleValues2["MaxLAT"]=46;
ListListDoubleValues2["LONPOLY"]={10, 10, 10};
ListListDoubleValues2["LATPOLY"]={10, 10, 10};
ListDoubleValues2["MaxDistTrackPointKM"]=20;
ListDoubleValues2["MinWIND"]=0;
ListDoubleValues2["MaxWIND"]=998;
ListDoubleValues2["MinHS"]=0;
ListDoubleValues2["MaxHS"]=998;
ListStringValues2["BEGTC"]="20110915.000000";
ListStringValues2["ENDTC"]="20110925.000000";
ListListStringValues2["AllowedSatellites"]=GetAllNamesOfSatelliteAltimeter();
ListBoolValues2["DoTrackSmoothing"]=false;
ListBoolValues2["EliminationShortTrack"]=false;
ListIntValues2["MinimalTrackSize"]=30;
ListBoolValues2["DoMinDistCoast"]=false;
ListDoubleValues2["MinDistCoastKM"]=60;
ListStringValues2["LonLatDiscFile"]="LonLatDisc.txt";
SingleBlock BlockSELECT;
BlockSELECT.ListIntValues=ListIntValues2;
BlockSELECT.ListBoolValues=ListBoolValues2;
BlockSELECT.ListDoubleValues=ListDoubleValues2;
BlockSELECT.ListListDoubleValues=ListListDoubleValues2;
BlockSELECT.ListStringValues=ListStringValues2;
BlockSELECT.ListListStringValues=ListListStringValues2;
BlockSELECT.BlockName=BlockName2;
ListBlock["SELECT"]=BlockSELECT;
std::string BlockName3="PROCESS";
std::map<std::string, int> ListIntValues3;
std::map<std::string, bool> ListBoolValues3;
std::map<std::string, double> ListDoubleValues3;
std::map<std::string, std::string> ListStringValues3;
std::map<std::string, std::vector<std::string> > ListListStringValues3;
ListBoolValues3["DO_WNDMAG"]=false;
ListBoolValues3["DO_HS"]=false;
ListBoolValues3["DO_STAT"]=false;
ListBoolValues3["DO_TXTRAW"]=false;
ListBoolValues3["DO_NCOUT"]=false;
ListBoolValues3["PLOT_ALL_TRACKS"]=false;
ListBoolValues3["PLOT_TRACKS"]=false;
ListBoolValues3["USE_CORRECTED"]=true;
ListBoolValues3["DO_SCATTERPLOT"]=true;
ListBoolValues3["SPATIALAVER"]=false;
ListIntValues3["MinEntryTrackPlot"]=30;
ListBoolValues3["DO_SAVE_TXT"]=false;
ListStringValues3["FILE_SAVE_TXT"]="alldata.txt";
ListBoolValues3["DO_SAVE_NC"]=false;
ListStringValues3["FILE_SAVE_NC"]="NeededFileName.nc";
SingleBlock BlockPROCESS;
BlockPROCESS.ListIntValues=ListIntValues3;
BlockPROCESS.ListBoolValues=ListBoolValues3;
BlockPROCESS.ListDoubleValues=ListDoubleValues3;
BlockPROCESS.ListStringValues=ListStringValues3;
BlockPROCESS.ListListStringValues=ListListStringValues3;
BlockPROCESS.BlockName=BlockName3;
ListBlock["PROCESS"]=BlockPROCESS;
FullNamelist eFullNamelist;
eFullNamelist.ListBlock=ListBlock;
eFullNamelist.FileName="undefined";
return eFullNamelist;
}
std::vector<SingleRecInterp> TRIG_FIND_ELE(MyMatrix<int> const& INE, MyMatrix<double> const& X, MyMatrix<double> const& Y, QuadCoordinate const& eQuad, MyMatrix<double> const& ListXY)
{
double THR=1e-10;
int mnp=X.rows();
auto IsCorrect=[&](int const& ie, double const& Xp, double const& Yp) -> bool {
int ki = INE(ie,0);
int kj = INE(ie,1);
int kk = INE(ie,2);
double xi = X(ki);
double yi = Y(ki);
double xj = X(kj);
double yj = Y(kj);
double xk = X(kk);
double yk = Y(kk);
# 9308 "/home/mathieu/GIT/wwmIII/CppOcean/AltimeterComparison.cpp"
double f1, f2, f3;
f1 = xi*(yj-Yp) + xj*(Yp-yi) + Xp*(yi-yj);
f2 = xj*(yk-Yp) + xk*(Yp-yj) + Xp*(yj-yk);
f3 = xk*(yi-Yp) + xi*(Yp-yk) + Xp*(yk-yi);
if (f1 > -THR && f2 > -THR && f3 > -THR)
return true;
return false;
};
double dx=0;
double dy=0;
int nbEle=INE.rows();
for (int ie=0; ie<nbEle; ie++) {
int ki = INE(ie,0);
int kj = INE(ie,1);
int kk = INE(ie,2);
double xi = X(ki);
double yi = Y(ki);
double xj = X(kj);
double yj = Y(kj);
double xk = X(kk);
double yk = Y(kk);
dx=std::max(dx, std::abs(xi - xj));
dx=std::max(dx, std::abs(xi - xk));
dx=std::max(dx, std::abs(xj - xk));
dy=std::max(dy, std::abs(yi - yj));
dy=std::max(dy, std::abs(yi - yk));
dy=std::max(dy, std::abs(yj - yk));
}
std::vector<int> ListAdjTrig=GetUnstructuredTriangleAdjInfo_vectint(INE, mnp);
auto SearchElement_V1=[&](double const& eX, double const& eY) -> int {
for (int iele=0; iele<nbEle; iele++)
if (IsCorrect(iele, eX, eY))
return iele;
return -1;
};
auto DistCentTriangle=[&](double const& eX, double const& eY, int const& iEle) -> double {
int i1=INE(iEle,0);
int i2=INE(iEle,1);
int i3=INE(iEle,2);
double eXcent=(X(i1) + X(i2) + X(i3))/double(3);
double eYcent=(X(i1) + X(i2) + X(i3))/double(3);
double delX=eXcent - eX;
double delY=eYcent - eY;
return sqrt(delX*delX + delY*delY);
};
auto SearchElement=[&](double const& eX, double const& eY, int const& iEltStart) -> int {
if (IsCorrect(iEltStart, eX, eY))
return iEltStart;
int iEleWork=iEltStart;
double distCurr=DistCentTriangle(eX, eY, iEleWork);
int nbIter=0;
while(1) {
bool DoSomething=false;
nbIter++;
for (int i=0; i<3; i++) {
int iEleAdj=ListAdjTrig[3*iEleWork + i];
if (iEleAdj != -1 && DoSomething == false) {
double eDist=DistCentTriangle(eX, eY, iEleAdj);
if (eDist < distCurr) {
iEleWork=iEleAdj;
distCurr=eDist;
DoSomething=true;
if (IsCorrect(iEleWork, eX, eY)) {
return iEleWork;
}
}
}
}
if (DoSomething == false)
break;
}
return SearchElement_V1(eX, eY);
};
int nbPoint=ListXY.cols();
int ielePrev=0;
std::vector<SingleRecInterp> LRec(nbPoint);
for (int iPoint=0; iPoint<nbPoint; iPoint++) {
double Xp=ListXY(0,iPoint);
double Yp=ListXY(1,iPoint);
int eElt=SearchElement(Xp, Yp, ielePrev);
if (eElt >= 0)
ielePrev=eElt;
std::cerr << "iPoint=" << iPoint << " eElt=" << eElt << "\n";
SingleRecInterp eRec;
if (eElt == -1) {
eRec={false, {}};
}
else {
std::vector<int> LEta(3);
std::vector<double> Xcall(3), Ycall(3);
for (int i=0; i<3; i++) {
int IP=INE(eElt,i);
LEta[i]=IP;
Xcall[i]=X(IP);
Ycall[i]=Y(IP);
}
std::vector<double> LCoeff=DetermineCoefficient(Xcall, Ycall, Xp, Yp);
std::vector<SinglePartInterp> LPart(3);
for (int i=0; i<3; i++) {
LPart[i]={LEta[i], 0, LCoeff[i]};
}
eRec={true, LPart};
};
LRec[iPoint]=eRec;
}
return LRec;
}
std::vector<SingleRecInterp> FD_FIND_ELE(CoordGridArrayFD const& CoordGridArr, QuadCoordinate const& eQuad, MyMatrix<double> const& ListXY)
{
double THR=1e-10;
MyMatrix<int> MatDir=GetDirection();
auto FindCoefficient=[&](int const& eEta, int const& eXi, double const& eX, double const& eY) -> std::vector<double> {
std::vector<double> X, Y, LCoeff;
double TheMin;
X={CoordGridArr.LON(eEta, eXi), CoordGridArr.LON(eEta+1, eXi), CoordGridArr.LON(eEta, eXi+1)};
Y={CoordGridArr.LAT(eEta, eXi), CoordGridArr.LAT(eEta+1, eXi), CoordGridArr.LAT(eEta, eXi+1)};
LCoeff=DetermineCoefficient(X, Y, eX, eY);
TheMin=VectorMin(LCoeff);
if (TheMin > -THR) {
double lambda1=LCoeff[1];
double lambda2=LCoeff[2];
double eCoeff00=(1-lambda1)*(1-lambda2);
double eCoeff10=lambda1*(1-lambda2);
double eCoeff01=(1-lambda1)*lambda2;
double eCoeff11=lambda1*lambda2;
std::vector<double> LCoeffRet{eCoeff00, eCoeff10, eCoeff01, eCoeff11};
return LCoeffRet;
}
X={CoordGridArr.LON(eEta+1, eXi+1), CoordGridArr.LON(eEta+1, eXi), CoordGridArr.LON(eEta, eXi+1)};
Y={CoordGridArr.LAT(eEta+1, eXi+1), CoordGridArr.LAT(eEta+1, eXi), CoordGridArr.LAT(eEta, eXi+1)};
LCoeff=DetermineCoefficient(X, Y, eX, eY);
TheMin=VectorMin(LCoeff);
if (TheMin > -THR) {
double lambda1=LCoeff[2];
double lambda2=LCoeff[1];
double eCoeff00=lambda1*lambda2;
double eCoeff10=(1-lambda1)*lambda2;
double eCoeff01=lambda1*(1-lambda2);
double eCoeff11=(1-lambda1)*(1-lambda2);
return {eCoeff00, eCoeff10, eCoeff01, eCoeff11};
}
return {-1,-1,-1,-1};
};
auto FindRecordArray=[&](int const& eEta, int const& eXi, double const& eX, double const& eY, SingleRecInterp& eRec) -> bool {
std::vector<double> LCoeff=FindCoefficient(eEta, eXi, eX, eY);
if (LCoeff[0] != -1) {
std::vector<SinglePartInterp> LPart(4);
double deltaX=eX;
double deltaY=eY;
for (int i=0; i<4; i++) {
int nEta=eEta + MatDir(0,i);
int nXi=eXi + MatDir(1,i);
SinglePartInterp ePart={nEta, nXi, LCoeff[i]};
deltaX=deltaX - LCoeff[i]*CoordGridArr.LON(nEta, nXi);
deltaY=deltaY - LCoeff[i]*CoordGridArr.LAT(nEta, nXi);
LPart[i]=ePart;
}
eRec={true, LPart};
return true;
}
return false;
};
int nbEta=CoordGridArr.LON.rows();
int nbXi =CoordGridArr.LON.cols();
auto AdmissibleEtaXi=[&](int const& eEta, int const& eXi) -> bool {
if (eEta >= 0 && eEta < nbEta-1 && eXi >= 0 && eXi < nbXi-1)
return true;
return false;
};
auto FindRecord=[&](int const& eEta, int const& eXi, double const& eX, double const& eY) -> SingleRecInterp {
bool testQuad=TestFeasibilityByQuad(eQuad, eX, eY);
if (testQuad == false) {
return {false, {}};
}
auto fEvaluateCorrectness=[&](int const& fEta, int const& fXi, bool& DoSomething, SingleRecInterp & eRec) -> bool {
if (AdmissibleEtaXi(fEta, fXi)) {
DoSomething=true;
bool test=FindRecordArray(fEta, fXi, eX, eY, eRec);
if (test == true) {
return true;
}
}
return false;
};
SingleRecInterp eRec;
int siz=0;
int sizExp=3;
for (siz=0; siz<sizExp; siz++) {
std::cerr << "siz=" << siz << "\n";
std::vector<int> ePair;
bool DoSomething=false;
if (siz == 0) {
if (fEvaluateCorrectness(eEta, eXi, DoSomething, eRec))
return eRec;
}
for (int i=-siz; i<siz; i++) {
if (fEvaluateCorrectness(eEta-siz, eXi+i , DoSomething, eRec))
return eRec;
if (fEvaluateCorrectness(eEta+i , eXi+siz, DoSomething, eRec))
return eRec;
if (fEvaluateCorrectness(eEta+siz , eXi-i, DoSomething, eRec))
return eRec;
if (fEvaluateCorrectness(eEta-i , eXi-siz, DoSomething, eRec))
return eRec;
}
if (DoSomething == false)
return {false, {}};
}
PairLL ePt{eX, eY};
PairCoord ePair=FindContaining(ePt, CoordGridArr.LON, CoordGridArr.LAT);
if (ePair.i == -1)
return {false, {}};
bool test=FindRecordArray(ePair.i, ePair.j, eX, eY, eRec);
if (test == false) {
std::cerr << "Inconsistency in the computation\n";
throw TerminalException{1};
}
return eRec;
};
int nbPoint=ListXY.cols();
int iEtaPrev=0;
int iXiPrev=0;
std::vector<SingleRecInterp> LRec;
for (int iPoint=0; iPoint<nbPoint; iPoint++) {
double Xp=ListXY(0,iPoint);
double Yp=ListXY(1,iPoint);
SingleRecInterp eEnt=FindRecord(iEtaPrev, iXiPrev, Xp, Yp);
LRec.push_back(eEnt);
std::cerr << "iPoint=" << iPoint << " / " << nbPoint << " x/y=" << Xp << "/" << Yp << " status=" << eEnt.status << "\n";
if (eEnt.status == true) {
iEtaPrev=eEnt.LPart[0].eEta;
iXiPrev=eEnt.LPart[0].eXi;
}
}
return LRec;
}
void Print_InterpolationError(std::vector<SingleRecInterp> const& LRec, GridArray const& GrdArr, MyMatrix<double> const& ListXY)
{
int nbPoint=LRec.size();
double TotalErr=0;
int nbPointInside=0;
for (int iPoint=0; iPoint<nbPoint; iPoint++) {
double deltaX=ListXY(0, iPoint);
double deltaY=ListXY(1, iPoint);
SingleRecInterp eSing=LRec[iPoint];
if (eSing.status == true) {
nbPointInside++;
for (auto& ePart : eSing.LPart) {
int eEta=ePart.eEta;
int eXi=ePart.eXi;
double eCoeff=ePart.eCoeff;
deltaX = deltaX - eCoeff*GrdArr.GrdArrRho.LON(eEta, eXi);
deltaY = deltaY - eCoeff*GrdArr.GrdArrRho.LAT(eEta, eXi);
}
double eErr=std::abs(deltaX) + std::abs(deltaY);
TotalErr=TotalErr + eErr;
}
}
std::cerr << "Total Interpolation error = " << TotalErr << "\n";
std::cerr << "nbPoint=" << nbPoint << " nbPointInside=" << nbPointInside << "\n";
}
QuadCoordinate Get_QuadCoordinate(GridArray const& GrdArr)
{
double MinLon=0;
double MaxLon=0;
double MinLat=0;
double MaxLat=0;
bool IsFirst=true;
int eta_rho=GrdArr.GrdArrRho.LON.rows();
int xi_rho =GrdArr.GrdArrRho.LON.cols();
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++) {
double eLon=GrdArr.GrdArrRho.LON(i,j);
double eLat=GrdArr.GrdArrRho.LAT(i,j);
if (IsFirst == true) {
MinLon=eLon;
MaxLon=eLon;
MinLat=eLat;
MaxLat=eLat;
}
else {
if (eLon > MaxLon)
MaxLon=eLon;
if (eLon < MinLon)
MinLon=eLon;
if (eLat > MaxLat)
MaxLat=eLat;
if (eLat < MinLat)
MinLat=eLat;
}
IsFirst=false;
}
return {MinLon, MaxLon, MinLat, MaxLat};
}
std::vector<SingleRecInterp> General_FindInterpolationWeight(GridArray const& GrdArr, MyMatrix<double> const& ListXY)
{
std::vector<SingleRecInterp> LRec;
QuadCoordinate eQuad=Get_QuadCoordinate(GrdArr);
if (GrdArr.IsFE == 0) {
std::cerr << "Before FD_FIND_ELE\n";
LRec=FD_FIND_ELE(GrdArr.GrdArrRho, eQuad, ListXY);
}
else {
std::cerr << "Before TRIG_FIND_ELE\n";
LRec=TRIG_FIND_ELE(GrdArr.INE, GrdArr.GrdArrRho.LON, GrdArr.GrdArrRho.LAT, eQuad, ListXY);
}
Print_InterpolationError(LRec, GrdArr, ListXY);
return LRec;
}
template<typename T>
MyMatrix<T> DimensionExtraction(Eigen::Tensor<T, 3> const& eT, size_t const& iDim, int const& eDim)
{
int n1=eT.dimension(0);
int n2=eT.dimension(1);
int n3=eT.dimension(2);
if (iDim == 0) {
MyMatrix<T> eMat(n2, n3);
for (int i2=0; i2<n2; i2++)
for (int i3=0; i3<n3; i3++)
eMat(i2,i3)=eT(eDim,i2,i3);
return eMat;
}
if (iDim == 1) {
MyMatrix<T> eMat(n1, n3);
for (int i1=0; i1<n1; i1++)
for (int i3=0; i3<n3; i3++)
eMat(i1,i3)=eT(i1,eDim,i3);
return eMat;
}
if (iDim == 2) {
MyMatrix<T> eMat(n1, n2);
for (int i1=0; i1<n1; i1++)
for (int i2=0; i2<n2; i2++)
eMat(i1,i2)=eT(i1,i2,eDim);
return eMat;
}
std::cerr << "Wrong input in ThreeDimArray\n";
std::cerr << "iDim=" << iDim << "\n";
std::cerr << "Allowed values: 0, 1, 2\n";
throw TerminalException{1};
}
template<typename T>
T ScalarProduct(MyVector<T> const& V1, MyVector<T> const & V2)
{
if (V1.size() != V2.size()) {
std::cerr << "Vectors of wrong sizes\n";
throw TerminalException{1};
}
size_t siz=V1.size();
T eSum=0;
for (size_t i=0; i<siz; i++)
eSum += V1(i)*V2(i);
return eSum;
}
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
throw TerminalException{1};
}
for (int iRow=0; iRow<nbRowI; iRow++)
for (int iCol=0; iCol<nbColI; iCol++) {
T eVal=eMatI(iRow, iCol);
eMatO(iRow, iCol)=eVal;
}
}
template<typename T>
MyVector<T> ProductVectorMatrix(MyVector<T> const& X, MyMatrix<T> const& M)
{
int nbCol=M.cols();
int nbRow=M.rows();
if (X.size() != nbRow) {
std::cerr << "Error in the product X A\n";
throw TerminalException{1};
}
MyVector<T> Vret(nbCol);
for (int iCol=0; iCol<nbCol; iCol++) {
T sum=0;
for (int iRow=0; iRow<nbRow; iRow++)
sum += M(iRow,iCol)*X(iRow);
Vret(iCol)=sum;
}
return Vret;
}
template<typename T>
MyVector<T> VectorMatrix(MyVector<T> const& eVect, MyMatrix<T> const& eMat)
{
int nbCol=eMat.cols();
int nbRow=eMat.rows();
int n=eVect.size();
if (n != nbRow) {
std::cerr << "n should be equal to nbRow\n";
throw TerminalException{1};
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
throw TerminalException{1};
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
Input(iCol , iColB)=prov2;
prov1=Output(iRowFound, iColB);
prov2=Output(iCol, iColB);
SwapValues(prov1, prov2);
Output(iRowFound, iColB)=prov1;
Output(iCol , iColB)=prov2;
}
}
}
}
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
template<typename T>
MyVector<T> CanonicalizeVector(MyVector<T> const& V)
{
int n=V.size();
T TheMin=0;
int iSelect=-1;
for (int i=0; i<n; i++) {
T eVal=V(i);
if (eVal != 0) {
T eAbs=T_abs(eVal);
if (iSelect == -1) {
TheMin=eAbs;
iSelect=i;
}
else {
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
T eQuot=1/T_abs(V(iSelect));
for (int i=0; i<n; i++)
Vret(i)=V(i)*eQuot;
return Vret;
}
std::string STRING_GETENV(std::string const& eStr)
{
char *ePre=std::getenv(eStr.c_str());
if (ePre == NULL) {
std::cerr << "Error in reading the environment variable : " << eStr << "\n";
throw TerminalException{1};
}
std::string eStrRet=ePre;
return eStrRet;
}
std::string DoubleTo4dot2f(double const& x)
{
char buffer[50];
int n=sprintf(buffer, "%4.2f", x);
if (n == 0) {
std::cerr << "Clear error in DoubleTo4dot2f\n";
throw TerminalException{1};
}
return std::string(buffer);
}
struct TempDirectory {
private:
std::string DirName;
public:
TempDirectory()
{
DirName="unset_and_irrelevant";
}
TempDirectory(std::string const& eDir)
{
DirName=eDir;
CreateDirectory(DirName);
}
TempDirectory& operator=(TempDirectory&& eTemp)
{
DirName=eTemp.str();
eTemp.DirName="unset_and_irrelevant";
return *this;
}
TempDirectory(TempDirectory && eTemp) : DirName(eTemp.str())
{
}
~TempDirectory()
{
if (DirName != "unset_and_irrelevant") {
if (IsExistingDirectory(DirName) == true) {
if (FILE_IsDirectoryEmpty(DirName) == false) {
std::cerr << "Keeping " << DirName << " since it is not empty\n";
}
else {
RemoveFile(DirName);
}
}
}
}
bool IsExisting() const
{
return IsExistingDirectory(DirName);
}
std::string str() const
{
return DirName;
}
};
struct PermanentInfoDrawing {
TempDirectory PrefixTemp;
std::string eDir;
std::string Extension;
std::string PicPrefix;
FullNamelist eFull;
bool KeepNC_NCL;
int NPROC;
DrawArr eDrawArr;
std::vector<QuadDrawInfo> ListQuadInfo;
std::vector<InterpolToUVpoints> ListInterpol;
};
TempDirectory PLOT_CreatePrefixTemp(FullNamelist const& eFull)
{
std::map<std::string, SingleBlock> ListBlock=eFull.ListBlock;
SingleBlock eBlPROC=ListBlock.at("PROC");
bool KeepNC_NCL=eBlPROC.ListBoolValues.at("KeepNC_NCL");
std::string eRand=random_string(20);
std::string Nature=eBlPROC.ListStringValues.at("__NaturePlot");
std::string PrefixTemp="/tmp/PLOT_" + Nature + "_" + eRand + "/";
if (KeepNC_NCL)
std::cerr << "PrefixTemp = " << PrefixTemp << "\n";
return TempDirectory(PrefixTemp);
}
PermanentInfoDrawing GET_PERMANENT_INFO(FullNamelist const& eFull)
{
SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
std::string PicPrefix=eBlPROC.ListStringValues.at("PicPrefix");
std::string Extension=eBlPROC.ListStringValues.at("Extension");
bool KeepNC_NCL=eBlPROC.ListBoolValues.at("KeepNC_NCL");
int NPROC=eBlPROC.ListIntValues.at("NPROC");
TempDirectory PrefixTemp=PLOT_CreatePrefixTemp(eFull);
if (eBlPROC.ListBoolValues.at("FirstCleanDirectory")) {
RemoveFileSpecificExtension(PicPrefix, Extension);
}
std::string eDir=FILE_GetAbsoluteDirectory(PicPrefix);
CreateDirectory(eDir);
DrawArr eDrawArr;
std::vector<QuadDrawInfo> ListQuadInfo;
std::vector<InterpolToUVpoints> ListInterpol;
PermanentInfoDrawing ePerm;
ePerm.PrefixTemp=std::move(PrefixTemp);
ePerm.eDir=eDir;
ePerm.Extension=Extension;
ePerm.PicPrefix=PicPrefix;
ePerm.eFull=eFull;
ePerm.KeepNC_NCL=KeepNC_NCL;
ePerm.NPROC=NPROC;
return ePerm;
}
void GeneralWriteSVGfile(std::string const& eFile, SVGplotDescription const& eSVGplot)
{
double MinX=0, MaxX=0, MinY=0, MaxY=0;
auto UpdateMinMaxXY=[&](coor const& pt) -> void {
if (pt.x > MaxX)
MaxX=pt.x;
if (pt.x < MinX)
MinX=pt.x;
if (pt.y > MaxY)
MaxY=pt.y;
if (pt.y < MinY)
MinY=pt.y;
};
for (auto & eLine : eSVGplot.ListLine) {
UpdateMinMaxXY(eLine.ePt);
UpdateMinMaxXY(eLine.fPt);
}
for (auto & ePolyline : eSVGplot.ListPolyline)
for (auto & eCoor : ePolyline.ListCoor)
UpdateMinMaxXY(eCoor);
for (auto& eBez : eSVGplot.ListBezier) {
UpdateMinMaxXY(eBez.pointM);
UpdateMinMaxXY(eBez.point2);
}
std::cerr << "SVG: X(min/max)=" << MinX << " / " << MaxX << "\n";
std::cerr << "SVG: Y(min/max)=" << MinY << " / " << MaxY << "\n";
double scale_factor, add_offsetX, add_offsetY;
double height=0, width=0;
bool FrameInit=false;
if (eSVGplot.FrameOption == 0) {
height=eSVGplot.height;
width=eSVGplot.width;
scale_factor=eSVGplot.scale_factor;
add_offsetX=eSVGplot.add_offsetX;
add_offsetY=eSVGplot.add_offsetY;
FrameInit=true;
}
if (eSVGplot.FrameOption == 1) {
height=eSVGplot.height;
width=eSVGplot.width;
double FrameX=eSVGplot.width;
double FrameY=eSVGplot.height;
double scale_factorX=FrameY / (MaxX - MinX);
double scale_factorY=FrameX / (MaxY - MinY);
double MidX=(MaxX + MinX) / 2;
double MidY=(MaxY + MinY) / 2;
scale_factor=T_min(scale_factorX, scale_factorY);
add_offsetX=FrameX/2 - scale_factor*MidX;
add_offsetY=FrameY/2 - scale_factor*MidY;
FrameInit=true;
}
if (FrameInit == false) {
std::cerr << "FrameOption has not been used\n";
std::cerr << "Please correct this\n";
throw TerminalException{1};
}
auto GetStringValue=[&](double const& eVal, double const& add_offset) -> std::string {
double eValM=add_offset + eVal*scale_factor;
if (eSVGplot.RoundMethod == 1)
return DoubleTo4dot2f(eValM);
if (eSVGplot.RoundMethod == 2)
return DoubleToString(eValM);
if (eSVGplot.RoundMethod == 3)
return DoubleToString(eValM);
std::cerr << "Failed to find relevant function\n";
throw TerminalException{1};
};
auto GetStringValueX=[&](double const& eVal) -> std::string {
return GetStringValue(eVal, add_offsetX);
};
auto GetStringValueY=[&](double const& eVal) -> std::string {
return GetStringValue(eVal, add_offsetY);
};
auto GetStringPair=[&](coor const& pt) -> std::string {
return GetStringValueX(pt.x) + " " + GetStringValueY(pt.y);
};
auto StringColor=[&](std::vector<int> const& eV) -> std::string {
return "rgb(" + IntToString(eV[0]) + "," + IntToString(eV[1]) + "," + IntToString(eV[2]) + ")";
};
auto GetQualityString=[&](SVGqualInfo const& eQual) -> std::string {
std::string eRet="style=\"stroke:" + StringColor(eQual.color) + ";stroke-width:" + IntToString(eQual.Size) + "\"";
if (eQual.MarkerEnd != "") {
eRet += " marker-end=\"url(#" + eQual.MarkerEnd + ")\"";
}
return eRet;
};
std::ofstream os(eFile);
os << "<svg height=\"" << height << "\" width=\"" << width << "\">\n";
std::cerr << "|ListLine|=" << eSVGplot.ListLine.size() << "\n";
for (auto & eLine : eSVGplot.ListLine) {
coor ePt=eLine.ePt;
coor fPt=eLine.fPt;
os << "  <line x1=\"" << GetStringValueX(ePt.x) << "\" y1=\"" << GetStringValueY(ePt.y) << "\" x2=\"" << GetStringValueX(fPt.x) << "\" y2=\"" << GetStringValueY(fPt.y) << "\" " << GetQualityString(eLine.eQual) << " />\n";
}
std::cerr << "|ListPolyline|=" << eSVGplot.ListPolyline.size() << "\n";
for (auto & ePolyline : eSVGplot.ListPolyline) {
os << "<polyline points=\"";
bool IsFirst=true;
for (auto & ePt : ePolyline.ListCoor) {
if (IsFirst == false)
os << " ";
os << GetStringValueX(ePt.x) << "," << GetStringValueY(ePt.y);
}
os << "\" " << GetQualityString(ePolyline.eQual) << " />\n";
}
std::cerr << "|ListBezier|=" << eSVGplot.ListBezier.size() << "\n";
for (auto& eBez : eSVGplot.ListBezier) {
os << "  <path d=\"M" << GetStringPair(eBez.pointM) << " C " << GetStringPair(eBez.pointC) << ", " << GetStringPair(eBez.point1) << ", " << GetStringPair(eBez.point2) << "\" fill=\"none\" " << GetQualityString(eBez.eQual) << " />\n";
}
os << "</svg>\n";
}
std::vector<SingleEntryAltimeter> READ_ALTI_FILE(std::string const& eFileAlti, std::vector<int> const& ListStatusSatellite, FullNamelist const& eFull)
{
double MinWind=eFull.ListBlock.at("SELECT").ListDoubleValues.at("MinWIND");
double MaxWind=eFull.ListBlock.at("SELECT").ListDoubleValues.at("MaxWIND");
double MinHs=eFull.ListBlock.at("SELECT").ListDoubleValues.at("MinHS");
double MaxHs=eFull.ListBlock.at("SELECT").ListDoubleValues.at("MaxHS");
std::cerr << "Read Altimetri file : " << eFileAlti << "\n";
std::vector<double> ListTime=NC_ReadTimeFromFile(eFileAlti, "time");
MyVector<double> ListLon=NC_Read1Dvariable(eFileAlti, "lon");
MyVector<double> ListLat=NC_Read1Dvariable(eFileAlti, "lat");
MyVector<double> ListWindSpeed=NC_Read1Dvariable(eFileAlti, "wind_speed");
MyVector<double> ListWindSpeed_cor=NC_Read1Dvariable(eFileAlti, "wind_speed_cor");
MyVector<double> ListSigma0=NC_Read1Dvariable(eFileAlti, "sigma0");
MyVector<double> ListSigma0cal=NC_Read1Dvariable(eFileAlti, "sigma0_cal");
MyVector<double> ListSigma0std=NC_Read1Dvariable(eFileAlti, "sigma0std");
MyVector<double> ListSigma0second=NC_Read1Dvariable(eFileAlti, "sigma0second");
MyVector<double> ListSigma0secondstd=NC_Read1Dvariable(eFileAlti, "sigma0secondstd");
MyVector<double> ListSwh=NC_Read1Dvariable(eFileAlti, "swh");
MyVector<double> ListSwhstd=NC_Read1Dvariable(eFileAlti, "swhstd");
MyVector<double> ListSwhcor=NC_Read1Dvariable(eFileAlti, "swhcor");
MyVector<int> ListSatellite=NC_Read1Dvariable_int(eFileAlti, "satellite");
int nbMeas=ListTime.size();
std::vector<int> ListStatus(nbMeas,1);
std::cerr << "sum(ListStatus:Before check) 1: " << VectorSum(ListStatus) << "\n";
COMPUTE_STATUS_LON_LAT(ListStatus, ListLon, ListLat, eFull);
std::cerr << "sum(ListStatus:After geography) 2: " << VectorSum(ListStatus) << "\n";
COMPUTE_STATUS_TIME(ListStatus, ListTime, eFull);
std::cerr << "sum(ListStatus:After time) 3: " << VectorSum(ListStatus) << "\n";
for (int iMeas=0; iMeas<nbMeas; iMeas++)
if (ListStatus[iMeas] == 1) {
if (ListStatusSatellite[ListSatellite[iMeas]-1] == 0)
ListStatus[iMeas]=0;
}
std::cerr << "sum(ListStatus:After chosen satellite) 4: " << VectorSum(ListStatus) << "\n";
bool USE_CORRECTED=eFull.ListBlock.at("PROCESS").ListBoolValues.at("USE_CORRECTED");
double LargeValue=1.0e31;
std::vector<SingleEntryAltimeter> ListEnt;
for (int iMeas=0; iMeas<nbMeas; iMeas++)
if (ListStatus[iMeas] == 1) {
SingleEntryAltimeter eEnt;
eEnt.Satellite=ListSatellite[iMeas];
eEnt.Time=ListTime[iMeas];
eEnt.Lon=ListLon(iMeas);
eEnt.Lat=ListLat(iMeas);
eEnt.WindSpeed=ListWindSpeed(iMeas);
eEnt.WindSpeed_cor=ListWindSpeed_cor(iMeas);
double WindSpeed_used;
if (USE_CORRECTED) {
WindSpeed_used=ListWindSpeed_cor(iMeas);
}
else {
WindSpeed_used=ListWindSpeed(iMeas);
}
if (WindSpeed_used > MaxWind || WindSpeed_used < MinWind) {
WindSpeed_used=LargeValue;
}
eEnt.WindSpeed_used=WindSpeed_used;
eEnt.Sigma0=ListSigma0(iMeas);
eEnt.Sigma0_cal=ListSigma0cal(iMeas);
eEnt.Sigma0std=ListSigma0std(iMeas);
eEnt.Sigma0second=ListSigma0second(iMeas);
eEnt.Sigma0secondstd=ListSigma0secondstd(iMeas);
eEnt.Swh=ListSwh(iMeas);
eEnt.SwhStd=ListSwhstd(iMeas);
eEnt.Swh_cor=ListSwhcor(iMeas);
double Swh_used;
if (USE_CORRECTED) {
Swh_used=ListSwhcor(iMeas);
}
else {
Swh_used=ListSwh(iMeas);
}
if (Swh_used > MaxHs || Swh_used < MinHs) {
Swh_used=LargeValue;
}
eEnt.Swh_used=Swh_used;
ListEnt.push_back(eEnt);
}
return ListEnt;
}
std::vector<int> FilterByMinDistCoast(std::vector<SingleEntryAltimeter> const& eList, FullNamelist const& eFull)
{
int nbEntry=eList.size();
std::vector<int> ListStatus(nbEntry,1);
std::string eFileCoast=eFull.ListBlock.at("SELECT").ListStringValues.at("LonLatDiscFile");
double MinDistCoastKM=eFull.ListBlock.at("SELECT").ListDoubleValues.at("MinDistCoastKM");
std::cerr << "eFileCoast=" << eFileCoast << "\n";
std::vector<PairLL> ListPtCoast=ReadLonLatDiscFile(eFileCoast);
std::cerr << "|ListPtCoast|=" << ListPtCoast.size() << "\n";
std::vector<PairLL> ListPt(nbEntry);
for (int iEntry=0; iEntry<nbEntry; iEntry++) {
double eLon=eList[iEntry].Lon;
double eLat=eList[iEntry].Lat;
PairLL ePt{eLon, eLat};
ListPt[iEntry]=ePt;
}
std::vector<double> ListMinDist=GetListMinimalDistances(ListPtCoast, ListPt);
std::vector<SingleEntryAltimeter> RetList;
for (int iEntry=0; iEntry<nbEntry; iEntry++)
if (ListMinDist[iEntry] < MinDistCoastKM)
ListStatus[iEntry]=0;
return ListStatus;
}
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
eSum1 =eSum1 + eMeas;
eSum2 =eSum2 + eModel;
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
eStat.str=eStat.strMeanError + " " + eStat.strAbsoluteError + " " + eStat.strRMSE + " " + eStat.strCenteredRMSE + " " + eStat.strCorrelation + " " + eStat.strScatterIndex + " " + eStat.strCenteredScatterIndex;
return eStat;
}
T_stat ComputeStatistics_vector(std::vector<double> const& ListMeas, std::vector<double> const& ListModel)
{
if (ListMeas.size() != ListModel.size()) {
std::cerr << "Error in ComputeStatistics_vector\n";
std::cerr << "Discrepancy in number of measurements\n";
std::cerr << "Please solve the problem\n";
throw TerminalException{1};
}
std::vector<PairMM> ListPair;
int nbEnt=ListMeas.size();
for (int iEnt=0; iEnt<nbEnt; iEnt++) {
ListPair.push_back({ListMeas[iEnt], ListModel[iEnt]});
}
return ComputeStatistics_Pair(ListPair);
}
T_stat ComputeStatistics_MyVector(MyVector<double> const& ListMeas, MyVector<double> const& ListModel)
{
if (ListMeas.size() != ListModel.size()) {
std::cerr << "Error in ComputeStatistics_MyVector\n";
std::cerr << "Discrepancy in number of measurements\n";
std::cerr << "Please solve the problem\n";
throw TerminalException{1};
}
std::vector<PairMM> ListPair;
int nbEnt=ListMeas.size();
for (int iEnt=0; iEnt<nbEnt; iEnt++) {
ListPair.push_back({ListMeas(iEnt), ListModel(iEnt)});
}
return ComputeStatistics_Pair(ListPair);
}
RecVar ModelSpecificVarSpecificTime_Kernel(TotalArrGetData const& TotalArr, std::string const& eVarName, double const& eTimeDay)
{
std::string eModelName=TotalArr.GrdArr.ModelName;
if (eModelName == "ROMS_IVICA")
eModelName = "ROMS";
if (eModelName == "WWM_DAILY")
eModelName = "WWM";
int eta_rho=TotalArr.GrdArr.GrdArrRho.LON.rows();
int xi_rho=TotalArr.GrdArr.GrdArrRho.LON.cols();
std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDay);
std::string strFile=DATE_ConvertMjd2mystringFile(eTimeDay);
RecSymbolic RecS;
RecS.eTimeDay=eTimeDay;
RecS.strPres="at " + strPres;
RecS.strFile=strFile;
RecS.VarNature="rho";
RecS.VarName1=eVarName;
RecS.VarName2="unset";
MyMatrix<double> F;
MyMatrix<double> U;
MyMatrix<double> V;
Eigen::Tensor<double,3> Tens3;
Eigen::Tensor<double,3> Uthree;
Eigen::Tensor<double,3> Vthree;
if (eVarName == "NbIterSolv") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "NB_ITER_SOLV", eTimeDay);
RecS.VarName2="nb Iteration Solver";
RecS.minval=0;
RecS.maxval=50;
RecS.mindiff=-5;
RecS.maxdiff=5;
RecS.Unit="nondim.";
}
if (eVarName == "CFL1") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "CFL1", eTimeDay);
RecS.VarName2="CFL1";
RecS.minval=0;
RecS.maxval=5;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="nondim.";
}
if (eVarName == "CFL2") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "CFL2", eTimeDay);
RecS.VarName2="CFL2";
RecS.minval=0;
RecS.maxval=5;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="nondim.";
}
if (eVarName == "CFL3") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "CFL3", eTimeDay);
RecS.VarName2="CFL3";
RecS.minval=0;
RecS.maxval=5;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="nondim.";
}
if (eVarName == "ThreeDfield1") {
if (eModelName == "WWM" || eModelName == "ROMS")
Tens3 = NETCDF_Get3DvariableSpecTime(TotalArr, "ThreeDfield1", eTimeDay);
RecS.VarName2="Generic three dim. field 1";
RecS.minval=0;
RecS.maxval=1;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.VarNature="3Drho";
RecS.Unit="unspecified";
}
if (eVarName == "FieldOut1") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "FieldOut1", eTimeDay);
RecS.VarName2="Generic Field Out 1";
RecS.minval=0;
RecS.maxval=1;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="unspecified";
}
if (eVarName == "IOBPWW3") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "IOBP_WW3", eTimeDay);
RecS.VarName2="IOBP of wavewatchIII";
RecS.minval=0;
RecS.maxval=1;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="nondim.";
}
if (eVarName == "MAPSTA") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "MAPSTA", eTimeDay);
RecS.VarName2="MAPSTA of wavewatchIII";
RecS.minval=-2;
RecS.maxval=2;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="nondim.";
}
if (eVarName == "Uwind") {
RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "WIND10", eTimeDay);
F=RecVarWork.U;
RecS.VarName2="Eastward wind";
RecS.minval=-10;
RecS.maxval=10;
RecS.mindiff=-2;
RecS.maxdiff=2;
RecS.Unit="m/s";
}
if (eVarName == "Vwind") {
RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "WIND10", eTimeDay);
F=RecVarWork.V;
RecS.VarName2="Northward wind";
RecS.minval=-10;
RecS.maxval=10;
RecS.mindiff=-2;
RecS.maxdiff=2;
RecS.Unit="m/s";
}
if (eVarName == "WIND10") {
if (eModelName == "SCHISM_SFLUX") {
U=Get2DvariableSpecTime(TotalArr, "uwind", eTimeDay);
V=Get2DvariableSpecTime(TotalArr, "vwind", eTimeDay);
}
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
RecS.VarName2="10m wind";
RecS.minval=0;
RecS.maxval=13;
RecS.mindiff=-2;
RecS.maxdiff=2;
RecS.Unit="m/s";
RecS.VarNature="uv";
RecS.nameU="Uwind";
RecS.nameV="Vwind";
}
if (eVarName == "WINDMAG") {
if (eModelName == "ROMS" || eModelName == "WWM") {
if (TOTALARR_IsVar(TotalArr, "Uwind") && TOTALARR_IsVar(TotalArr, "Vwind") ) {
MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "Uwind", eTimeDay);
MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "Vwind", eTimeDay);
F=COMPUTE_NORM(Us, Vs);
}
else {
if (eModelName == "WWM") {
F=Get2DvariableSpecTime(TotalArr, "WINDMAG", eTimeDay);
}
else {
F=Get2DvariableSpecTime(TotalArr, "WNDMAG", eTimeDay);
}
}
}
if (eModelName == "SCHISM_SFLUX") {
MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "uwind", eTimeDay);
MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "vwind", eTimeDay);
F=COMPUTE_NORM(Us, Vs);
}
if (eModelName == "COSMO" || eModelName == "WAM") {
MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "U_10", eTimeDay);
MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "V_10", eTimeDay);
F=COMPUTE_NORM(Us, Vs);
}
if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_COSMO") {
MyMatrix<double> Us=Get2DvariableSpecTime(TotalArr, "10u", eTimeDay);
MyMatrix<double> Vs=Get2DvariableSpecTime(TotalArr, "10v", eTimeDay);
F=COMPUTE_NORM(Us, Vs);
}
if (eModelName == "GRIB_WAM_FORT30") {
F=Get2DvariableSpecTime(TotalArr, "wind", eTimeDay);
}
RecS.VarName2="10m wind speed";
RecS.minval=0;
RecS.maxval=13;
RecS.mindiff=-2;
RecS.maxdiff=2;
RecS.Unit="m/s";
}
if (eVarName == "AIRD") {
if (eModelName == "COSMO" || eModelName == "WAM")
F=Get2DvariableSpecTime(TotalArr, "AIRD", eTimeDay);
RecS.VarName2="air density";
RecS.minval=1.12;
RecS.maxval=1.20;
RecS.mindiff=-0.02;
RecS.maxdiff=0.02;
RecS.Unit="kg/m3";
}
if (eVarName == "rain") {
if (eModelName == "SCHISM_SFLUX")
F=Get2DvariableSpecTime(TotalArr, "prate", eTimeDay);
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "rain", eTimeDay);
if (eModelName == "GRIB_COSMO" || eModelName == "GRIB_DWD" || eModelName == "GRIB_ECMWF")
F=GRID_Get2DVariableTimeDifferentiate(TotalArr, "tp", eTimeDay);
RecS.VarName2="rainfall rate";
RecS.minval=0;
RecS.maxval=0.001;
RecS.mindiff=-0.001;
RecS.maxdiff=0.001;
RecS.Unit="kg/m^2/s";
}
if (eVarName == "swrad") {
if (eModelName == "SCHISM_SFLUX")
F=Get2DvariableSpecTime(TotalArr, "dswrf", eTimeDay);
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "swrad", eTimeDay);
if (eModelName == "GRIB_COSMO")
F=Get2DvariableSpecTime(TotalArr, "sobs_rad", eTimeDay);
if (eModelName == "GRIB_ECMWF")
F=GRID_Get2DVariableTimeDifferentiate(TotalArr, "ssrd", eTimeDay);
RecS.VarName2="Shortwave flux";
RecS.minval=100;
RecS.maxval=1000;
RecS.mindiff=-100;
RecS.maxdiff=100;
RecS.Unit="W/m2";
}
if (eVarName == "lwrad") {
if (eModelName == "SCHISM_SFLUX")
F=Get2DvariableSpecTime(TotalArr, "dlwrf", eTimeDay);
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "lwrad", eTimeDay);
if (eModelName == "GRIB_COSMO")
F=Get2DvariableSpecTime(TotalArr, "thbs_rad", eTimeDay);
if (eModelName == "GRIB_ECMWF")
F=GRID_Get2DVariableTimeDifferentiate(TotalArr, "strd", eTimeDay);
RecS.VarName2="Longwave flux";
RecS.minval=200;
RecS.maxval=500;
RecS.mindiff=-50;
RecS.maxdiff= 50;
RecS.Unit="W/m2";
}
if (eVarName == "latent") {
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "latent", eTimeDay);
RecS.VarName2="Latent flux";
RecS.minval=0;
RecS.maxval=0.033;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="W/m2";
}
if (eVarName == "SurfPres") {
if (eModelName == "ROMS") {
MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "Pair", eTimeDay);
F=100*Fin;
}
if (eModelName == "SCHISM_SFLUX") {
F=Get2DvariableSpecTime(TotalArr, "prmsl", eTimeDay);
}
if (eModelName == "GRIB_DWD" || eModelName == "GRIB_GFS")
F=Get2DvariableSpecTime(TotalArr, "prmsl", eTimeDay);
if (eModelName == "GRIB_ECMWF")
F=Get2DvariableSpecTime(TotalArr, "msl", eTimeDay);
if (eModelName == "GRIB_COSMO")
F=Get2DvariableSpecTime(TotalArr, "pmsl", eTimeDay);
RecS.VarName2="mean sea level pressure";
RecS.minval=100000;
RecS.maxval=103000;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="Pa";
}
if (eVarName == "sensible") {
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "sensible", eTimeDay);
RecS.VarName2="Sensible heat flux";
RecS.minval=0;
RecS.maxval=0.033;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="W/m2";
}
if (eVarName == "shflux") {
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "shflux", eTimeDay);
RecS.VarName2="Surface heat flux";
RecS.minval=0;
RecS.maxval=0.033;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="W/m2";
}
if (eVarName == "ssflux") {
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "ssflux", eTimeDay);
RecS.VarName2="Surface salinity flux";
RecS.minval=0;
RecS.maxval=0.033;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="PSU/m2s";
}
if (eVarName == "evaporation") {
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "evaporation", eTimeDay);
RecS.VarName2="Evaporation rate";
RecS.minval=0;
RecS.maxval=0.033;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="kg/m2s";
}
if (eVarName == "AIRT2K") {
RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "AIRT2", eTimeDay);
F=RecVarWork.F;
int siz=F.size();
for (int i=0; i<siz; i++)
F(i) += double(273.15);
RecS.VarName2="2m air temperature (K)";
RecS.minval=273.15 + 10;
RecS.maxval=273.15 + 20;
RecS.mindiff=-2;
RecS.maxdiff=2;
RecS.Unit="deg";
}
if (eVarName == "AIRT2") {
if (eModelName == "SCHISM_SFLUX") {
F=Get2DvariableSpecTime(TotalArr, "stmp", eTimeDay);
int siz=F.size();
for (int i=0; i<siz; i++)
F(i) -= double(273.15);
}
if (eModelName == "COSMO") {
F=Get2DvariableSpecTime(TotalArr, "t_2m", eTimeDay);
int siz=F.size();
for (int i=0; i<siz; i++)
F(i) -= double(273.15);
}
if (eModelName == "GRIB_DWD" || eModelName == "GRIB_ECMWF" || eModelName == "GRIB_GFS" || eModelName == "GRIB_COSMO") {
F=Get2DvariableSpecTime(TotalArr, "2t", eTimeDay);
int siz=F.size();
for (int i=0; i<siz; i++)
F(i) -= double(273.15);
}
RecS.VarName2="2m air temperature";
RecS.minval=10;
RecS.maxval=20;
RecS.mindiff=-2;
RecS.maxdiff=2;
RecS.Unit="deg";
}
if (eVarName == "Rh2frac") {
RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "Rh2", eTimeDay);
F=RecVarWork.F / double(100);
RecS.VarName2="2m relative humidity";
RecS.minval=0;
RecS.maxval=1;
RecS.mindiff=-0.2;
RecS.maxdiff=0.2;
RecS.Unit="nondim.";
}
if (eVarName == "Rh2") {
if (eModelName == "COSMO")
F=Get2DvariableSpecTime(TotalArr, "rh_2m", eTimeDay);
if (eModelName == "GRIB_DWD")
F=Get2DvariableSpecTime(TotalArr, "RELHUM_2M", eTimeDay);
if (eModelName == "SCHISM_SFLUX") {
MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "spfh", eTimeDay);
F=100 * Fin;
}
if (eModelName == "GRIB_ECMWF") {
if (TOTALARR_IsVar(TotalArr, "2r")) {
F=Get2DvariableSpecTime(TotalArr, "2r", eTimeDay);
}
else {
int MethodRH = 2;
if (MethodRH == 1) {
double airDens=1.225;
double waterDens=0.804;
double quot=waterDens / (airDens - waterDens);
MyMatrix<double> Fspecific=Get2DvariableSpecTime(TotalArr, "q", eTimeDay);
double TheMult=100 / ( 0.622 * quot);
F = Fspecific * TheMult;
}
if (MethodRH == 2) {
MyMatrix<double> F_q=Get2DvariableSpecTime(TotalArr, "q", eTimeDay);
MyMatrix<double> F_p;
if (TOTALARR_IsVar(TotalArr, "msl")) {
F_p=Get2DvariableSpecTime(TotalArr, "msl", eTimeDay);
}
else {
int eta=F_q.rows();
int xi=F_q.cols();
F_p=MyMatrix<double>(eta,xi);
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++)
F_p(i,j)=103000;
}
MyMatrix<double> F_TK=Get2DvariableSpecTime(TotalArr, "2t", eTimeDay);
int eta=F_q.rows();
int xi=F_q.cols();
F=MyMatrix<double>(eta,xi);
for (int i=0; i<eta; i++)
for (int j=0; j<xi; j++) {
double eT=F_TK(i,j);
double eQ=F_q(i,j);
double eP=F_p(i,j);
double eT0=double(273.15);
double TheQuot=double(17.67) * (eT - eT0)/(eT - double(29.65));
double eRH=0.263 * eP *eQ /(exp(TheQuot));
F(i,j)=T_min(eRH, double(100));
}
}
}
}
RecS.VarName2="2m relative humidity";
RecS.minval=0;
RecS.maxval=100;
RecS.mindiff=-20;
RecS.maxdiff=20;
RecS.Unit="nondim.";
}
if (eVarName == "UsurfCurr") {
RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "SurfCurr", eTimeDay);
F=RecVarWork.U;
RecS.VarName2="Eastward current";
RecS.minval=-0.3;
RecS.maxval=0.3;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="m/s";
}
if (eVarName == "VsurfCurr") {
RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "SurfCurr", eTimeDay);
F=RecVarWork.V;
RecS.VarName2="Northward current";
RecS.minval=-0.3;
RecS.maxval=0.3;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="m/s";
}
if (eVarName == "Curr") {
if (eModelName == "ROMS") {
Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
Uthree=My_u2rho_3D(Utot, TotalArr.GrdArr.GrdArrU.MSK);
Vthree=My_v2rho_3D(Vtot, TotalArr.GrdArr.GrdArrV.MSK);
}
if (eModelName == "WWM") {
Uthree=NETCDF_Get3DvariableSpecTime(TotalArr, "Ucurr", eTimeDay);
Vthree=NETCDF_Get3DvariableSpecTime(TotalArr, "Vcurr", eTimeDay);
}
RecS.VarName2="baroclinic current";
RecS.minval=0;
RecS.maxval=0.2;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.VarNature="3Duv";
RecS.Unit="m/s";
}
if (eVarName == "CurrMag") {
RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "Curr", eTimeDay);
auto LDim=RecVarWork.Uthree.dimensions();
int dim0=LDim[0];
int dim1=LDim[1];
int dim2=LDim[2];
Eigen::Tensor<double,3> Tens3(dim0, dim1, dim2);
for (int i0=0; i0<dim0; i0++)
for (int i1=0; i1<dim1; i1++)
for (int i2=0; i2<dim2; i2++) {
double eU=RecVarWork.Uthree(i0, i1, i2);
double eV=RecVarWork.Vthree(i0, i1, i2);
double eNorm=sqrt(eU*eU + eV*eV);
Tens3(i0, i1, i2) = eNorm;
}
RecS.VarName2="baroclinic current magnitude";
RecS.minval=0;
RecS.maxval=0.2;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.VarNature="3Drho";
RecS.Unit="m/s";
}
if (eVarName == "SurfCurr") {
if (eModelName == "ROMS") {
Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "u", eTimeDay);
Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "v", eTimeDay);
int s_rho=Utot.dimension(0);
MyMatrix<double> Usurf=DimensionExtraction(Utot, 0, s_rho-1);
MyMatrix<double> Vsurf=DimensionExtraction(Vtot, 0, s_rho-1);
U=My_u2rho(Usurf, TotalArr.GrdArr.GrdArrU.MSK);
V=My_v2rho(Vsurf, TotalArr.GrdArr.GrdArrV.MSK);
}
if (eModelName == "WWM") {
if (TOTALARR_IsVar(TotalArr, "CURTX") && TOTALARR_IsVar(TotalArr, "CURTX") ) {
U=Get2DvariableSpecTime(TotalArr, "CURTX", eTimeDay);
V=Get2DvariableSpecTime(TotalArr, "CURTY", eTimeDay);
}
else {
if (TOTALARR_IsVar(TotalArr, "UsurfCurr") && TOTALARR_IsVar(TotalArr, "VsurfCurr") ) {
U=Get2DvariableSpecTime(TotalArr, "UsurfCurr", eTimeDay);
V=Get2DvariableSpecTime(TotalArr, "VsurfCurr", eTimeDay);
}
else {
Eigen::Tensor<double,3> Utot=NETCDF_Get3DvariableSpecTime(TotalArr, "Ucurr", eTimeDay);
Eigen::Tensor<double,3> Vtot=NETCDF_Get3DvariableSpecTime(TotalArr, "Vcurr", eTimeDay);
int s_rho=Utot.dimension(0);
U=DimensionExtraction(Utot, 0, s_rho-1);
V=DimensionExtraction(Vtot, 0, s_rho-1);
}
}
}
if (eModelName == "COSMO" || eModelName == "WAM") {
U=Get2DvariableSpecTime(TotalArr, "ucurr", eTimeDay);
V=Get2DvariableSpecTime(TotalArr, "vcurr", eTimeDay);
}
AngleRhoRot(U, V, TotalArr.GrdArr);
RecS.VarName2="surface current";
RecS.minval=0;
RecS.maxval=0.5;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="m/s";
RecS.VarNature="uv";
RecS.nameU="UsurfCurr";
RecS.nameV="VsurfCurr";
}
if (eVarName == "SurfCurrMag") {
RecVar RecVarWork=ModelSpecificVarSpecificTime_Kernel(TotalArr, "SurfCurr", eTimeDay);
F=COMPUTE_NORM(RecVarWork.U, RecVarWork.V);
RecS.VarName2="surface current magnitude";
RecS.minval=0;
RecS.maxval=0.5;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="m/s";
}
if (eVarName == "TempSurf") {
if (eModelName == "ROMS") {
Eigen::Tensor<double,3> TheTemp=NETCDF_Get3DvariableSpecTime(TotalArr, "temp", eTimeDay);
int s_rho=TheTemp.dimension(0);
F=DimensionExtraction(TheTemp, 0, s_rho-1);
}
if (eModelName == "COSMO") {
F=Get2DvariableSpecTime(TotalArr, "t_s", eTimeDay);
int siz=F.size();
for (int i=0; i<siz; i++)
F(i) -= double(273.15);
}
RecS.VarName2="sea surface temperature";
RecS.minval=10;
RecS.maxval=20;
RecS.mindiff=-2;
RecS.maxdiff=2;
RecS.Unit="deg";
}
if (eVarName == "SaltSurf") {
if (eModelName == "ROMS") {
Eigen::Tensor<double,3> TheSalt=NETCDF_Get3DvariableSpecTime(TotalArr, "salt", eTimeDay);
int s_rho=TheSalt.dimension(0);
F=DimensionExtraction(TheSalt, 0, s_rho-1);
}
RecS.VarName2="sea surface salinity";
RecS.minval=30;
RecS.maxval=40;
RecS.mindiff=-2;
RecS.maxdiff=2;
RecS.Unit="PSU";
}
if (eVarName == "ZetaOcean") {
if (eModelName == "COSMO")
F=Get2DvariableSpecTime(TotalArr, "ZetaOcean", eTimeDay);
if (eModelName == "ROMS")
F=Get2DvariableSpecTime(TotalArr, "zeta", eTimeDay);
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
RecS.VarName2="free surface elevation";
RecS.minval=-0.2;
RecS.maxval=0.2;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="m";
}
if (eVarName == "ZetaOceanDerivative") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "DEPDT", eTimeDay);
RecS.VarName2="free surface elevation derivative";
RecS.minval=-0.01;
RecS.maxval=0.01;
RecS.mindiff=-0.001;
RecS.maxdiff=0.001;
RecS.Unit="m/s";
}
if (eVarName == "MeanWaveLength") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "WLM", eTimeDay);
RecS.VarName2="mean wave length";
RecS.minval=2;
RecS.maxval=30;
RecS.mindiff=-5;
RecS.maxdiff=5;
RecS.Unit="m";
}
if (eVarName == "PeakWaveLength") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "LPP", eTimeDay);
RecS.VarName2="peak wave length";
RecS.minval=2;
RecS.maxval=30;
RecS.mindiff=-5;
RecS.maxdiff=5;
RecS.Unit="m";
}
if (eVarName == "MeanWaveNumber") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "KLM", eTimeDay);
RecS.VarName2="mean wave number";
RecS.minval=0;
RecS.maxval=1;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="m-1";
}
if (eVarName == "PeakWaveNumber") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "KPP", eTimeDay);
RecS.VarName2="peak wave number";
RecS.minval=0;
RecS.maxval=1;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="m-1";
}
if (eVarName == "MeanWaveDir") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "DM", eTimeDay);
RecS.VarName2="mean wave direction";
RecS.minval=0;
RecS.maxval=360;
RecS.mindiff=-30;
RecS.maxdiff=30;
RecS.Unit="deg";
}
if (eVarName == "PeakWaveDir") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "PEAKD", eTimeDay);
RecS.VarName2="peak wave direction";
RecS.minval=0;
RecS.maxval=360;
RecS.mindiff=-30;
RecS.maxdiff=30;
RecS.Unit="deg";
}
if (eVarName == "MeanWaveDirVect") {
if (eModelName == "WWM") {
F=Get2DvariableSpecTime(TotalArr, "DM", eTimeDay);
int nbRow=F.rows();
int nbCol=F.cols();
double deg2rad=3.1415926535 / double(180);
U=MyMatrix<double>(nbRow,nbCol);
V=MyMatrix<double>(nbRow,nbCol);
for (int iRow=0; iRow<nbRow; iRow++)
for (int iCol=0; iCol<nbCol; iCol++) {
double eAngRad=deg2rad*F(iRow,iCol);
U(iRow,iCol)=cos(eAngRad);
V(iRow,iCol)=sin(eAngRad);
}
}
RecS.VarName2="mean wave direction";
RecS.minval=0;
RecS.maxval=360;
RecS.mindiff=-30;
RecS.maxdiff=30;
RecS.VarNature="uv";
RecS.Unit="deg";
}
if (eVarName == "PeakWaveDirVect") {
if (eModelName == "WWM") {
F=Get2DvariableSpecTime(TotalArr, "PEAKD", eTimeDay);
int nbRow=F.rows();
int nbCol=F.cols();
double deg2rad=3.1415926535 / double(180);
U=MyMatrix<double>(nbRow,nbCol);
V=MyMatrix<double>(nbRow,nbCol);
for (int iRow=0; iRow<nbRow; iRow++)
for (int iCol=0; iCol<nbCol; iCol++) {
double eAngRad=deg2rad*F(iRow,iCol);
U(iRow,iCol)=cos(eAngRad);
V(iRow,iCol)=sin(eAngRad);
}
}
RecS.VarName2="peak wave direction";
RecS.minval=0;
RecS.maxval=360;
RecS.mindiff=-30;
RecS.maxdiff=30;
RecS.VarNature="uv";
RecS.Unit="deg";
}
if (eVarName == "DiscPeakWaveDir") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "DPEAK", eTimeDay);
RecS.VarName2="discrete peak wave direction";
RecS.minval=0;
RecS.maxval=360;
RecS.mindiff=-30;
RecS.maxdiff=30;
RecS.Unit="deg";
}
if (eVarName == "ZetaSetup") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "ZETA_SETUP", eTimeDay);
RecS.VarName2="free surface setup";
RecS.minval=0;
RecS.maxval=0.76;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="m";
}
if (eVarName == "BreakingFraction") {
MyMatrix<double> Fhs, Fzeta;
if (eModelName == "WWM")
Fhs=Get2DvariableSpecTime(TotalArr, "HS", eTimeDay);
if (eModelName == "WWM")
Fzeta=Get2DvariableSpecTime(TotalArr, "WATLEV", eTimeDay);
F=MyMatrix<double>(eta_rho, xi_rho);
for (int i=0; i<eta_rho; i++)
for (int j=0; j<xi_rho; j++)
F(i,j)=Fhs(i,j) / (Fzeta(i,j) + TotalArr.GrdArr.GrdArrRho.DEP(i,j));
RecS.VarName2="Breaking fraction";
RecS.minval=0;
RecS.maxval=0.76;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="nondim.";
}
if (eVarName == "Hwave") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "HS", eTimeDay);
if (eModelName == "WW3")
F=Get2DvariableSpecTime(TotalArr, "hs", eTimeDay);
if (eModelName == "COSMO" || eModelName == "WAM")
F=Get2DvariableSpecTime(TotalArr, "Hwave", eTimeDay);
if (eModelName == "GRIB_WAM_FORT30")
F=Get2DvariableSpecTime(TotalArr, "swh", eTimeDay);
RecS.VarName2="Significant wave height";
RecS.minval=0;
RecS.maxval=4.5;
RecS.mindiff=-0.5;
RecS.maxdiff=0.5;
RecS.Unit="m";
}
if (eVarName == "MeanWaveFreq") {
if (eModelName == "COSMO" || eModelName == "WAM")
F=Get2DvariableSpecTime(TotalArr, "MwaveFreq", eTimeDay);
if (eModelName == "WWM") {
MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "TM01", eTimeDay);
F=FreqPeriodChange(Fin);
}
RecS.VarName2="mean wave frequency";
RecS.minval=0;
RecS.maxval=0.9;
RecS.mindiff=-0.2;
RecS.maxdiff=0.2;
RecS.Unit="Hz";
}
if (eVarName == "PeakWaveFreq") {
if (eModelName == "COSMO" || eModelName == "WAM")
F=Get2DvariableSpecTime(TotalArr, "PwaveFreq", eTimeDay);
if (eModelName == "WWM") {
MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "TPP", eTimeDay);
F=FreqPeriodChange(Fin);
}
RecS.VarName2="peak wave frequency";
RecS.minval=0;
RecS.maxval=0.9;
RecS.mindiff=-0.2;
RecS.maxdiff=0.2;
RecS.Unit="Hz";
}
if (eVarName == "MeanWavePer") {
if (eModelName == "COSMO" || eModelName == "WAM") {
MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "MwaveFreq", eTimeDay);
F=FreqPeriodChange(Fin);
}
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "TM01", eTimeDay);
RecS.VarName2="mean wave period";
RecS.minval=2;
RecS.maxval=10;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="s";
}
if (eVarName == "PeakWavePer") {
if (eModelName == "COSMO" || eModelName == "WAM") {
MyMatrix<double> Fin=Get2DvariableSpecTime(TotalArr, "PwaveFreq", eTimeDay);
F=FreqPeriodChange(Fin);
}
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "TPP", eTimeDay);
RecS.VarName2="peak wave period";
RecS.minval=2;
RecS.maxval=10;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="s";
}
if (eVarName == "TM02") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "TM02", eTimeDay);
RecS.VarName2="zero crossing wave period";
RecS.minval=2;
RecS.maxval=10;
RecS.mindiff=-1;
RecS.maxdiff=1;
RecS.Unit="s";
}
if (eVarName == "DynBathy") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "DW", eTimeDay);
RecS.VarName2="dynamic bathymetry";
RecS.minval=0;
RecS.maxval=30;
RecS.mindiff=-5;
RecS.maxdiff=5;
RecS.Unit="deg";
}
if (eVarName == "MeanWaveDirSpread") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "DSPR", eTimeDay);
RecS.VarName2="directional spreading";
RecS.minval=0;
RecS.maxval=30;
RecS.mindiff=-5;
RecS.maxdiff=5;
RecS.Unit="deg";
}
if (eVarName == "PeakWaveDirSpread") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "PEAKDSPR", eTimeDay);
RecS.VarName2="peak directional spreading";
RecS.minval=0;
RecS.maxval=30;
RecS.mindiff=-5;
RecS.maxdiff=5;
RecS.Unit="deg";
}
if (eVarName == "AirZ0") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "Z0", eTimeDay);
RecS.VarName2="air roughness length";
RecS.minval=0;
RecS.maxval= 0.0002;
RecS.mindiff=-0.00005;
RecS.maxdiff= 0.00005;
RecS.Unit="m";
}
if (eVarName == "AirFricVel") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "UFRIC", eTimeDay);
RecS.VarName2="air roughness length";
RecS.minval=0;
RecS.maxval=0.3;
RecS.mindiff=-0.05;
RecS.maxdiff= 0.05;
RecS.Unit="m";
}
if (eVarName == "CdWave") {
if (eModelName == "COSMO" || eModelName == "WAM")
F=Get2DvariableSpecTime(TotalArr, "CdWave", eTimeDay);
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "CD", eTimeDay);
RecS.VarName2="drag coefficient from the wave model";
RecS.minval=0;
RecS.maxval=0.20;
RecS.mindiff=-0.05;
RecS.maxdiff=0.05;
RecS.Unit="nondim.";
}
if (eVarName == "AlphaWave") {
if (eModelName == "COSMO" || eModelName == "WAM")
F=Get2DvariableSpecTime(TotalArr, "AlphaWave", eTimeDay);
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "ALPHA_CH", eTimeDay);
RecS.VarName2="Charnock coefficient from the wave model";
RecS.minval=0;
RecS.maxval=0.033;
RecS.mindiff=-0.1;
RecS.maxdiff=0.1;
RecS.Unit="nondim.";
}
if (eVarName == "TotSurfStr") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "TAUTOT", eTimeDay);
RecS.VarName2="Total Surface stress";
RecS.minval=0;
RecS.maxval=0.06;
RecS.mindiff=-0.01;
RecS.maxdiff= 0.01;
RecS.Unit="unknown";
}
if (eVarName == "WaveSurfStr") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "TAUW", eTimeDay);
RecS.VarName2="wave supported Surface stress";
RecS.minval=0;
RecS.maxval=0.06;
RecS.mindiff=-0.01;
RecS.maxdiff= 0.01;
RecS.Unit="unknown";
}
if (eVarName == "SurfStrHF") {
if (eModelName == "WWM")
F=Get2DvariableSpecTime(TotalArr, "TAUHF", eTimeDay);
RecS.VarName2="high frequency Surface stress";
RecS.minval=0;
RecS.maxval=0.06;
RecS.mindiff=-0.01;
RecS.maxdiff= 0.01;
RecS.Unit="unknown";
}
RecVar eRecVar;
eRecVar.RecS=RecS;
if (RecS.VarName2 == "unset") {
std::cerr << "We did not find the variable\n";
std::cerr << "eVarName = " << eVarName << "\n";
std::cerr << "in the list of allowed ones\n";
std::cerr << "possibly missspelling or lack of relevant code\n";
throw TerminalException{1};
}
if (eModelName != "TRIVIAL") {
if (RecS.VarNature == "rho") {
if (F.size() == 0) {
std::cerr << "VarNature = " << RecS.VarNature << "\n";
std::cerr << "Variable eVarName = " << eVarName << "\n";
std::cerr << "is recognized by the program as authorized variable\n";
std::cerr << "But it has not been assigned.\n";
std::cerr << "Possibly because of missing facility for\n";
std::cerr << "eModelName = " << eModelName << "\n";
throw TerminalException{1};
}
eRecVar.F=F;
}
if (RecS.VarNature == "uv") {
if (U.size() == 0 || V.size() == 0) {
std::cerr << "VarNature = " << RecS.VarNature << "\n";
std::cerr << "Variable eVarName = " << eVarName << "\n";
std::cerr << "is recognized by the program\n";
std::cerr << "But it has not been assigned.\n";
std::cerr << "Possibly because of missing facility for\n";
std::cerr << "eModelName = " << eModelName << "\n";
throw TerminalException{1};
}
eRecVar.U=U;
eRecVar.V=V;
int nbRow=U.rows();
int nbCol=U.cols();
MyMatrix<double> Fwr(nbRow, nbCol);
for (int iRow=0; iRow<nbRow; iRow++)
for (int iCol=0; iCol<nbCol; iCol++) {
double eU=U(iRow,iCol);
double eV=V(iRow,iCol);
double eNorm=sqrt(eU*eU + eV*eV);
Fwr(iRow,iCol) = eNorm;
}
eRecVar.F=Fwr;
}
if (RecS.VarNature == "3Drho") {
auto LDim=Tens3.dimensions();
if (LDim[0] == 0) {
std::cerr << "VarNature = " << RecS.VarNature << "\n";
std::cerr << "Variable eVarName = " << eVarName << "\n";
std::cerr << "is recognized by the program\n";
std::cerr << "But it has not been assigned.\n";
std::cerr << "Possibly because of missing facility for\n";
std::cerr << "eModelName = " << eModelName << "\n";
throw TerminalException{1};
}
eRecVar.Tens3=Tens3;
}
if (RecS.VarNature == "3Duv") {
auto LDim=Uthree.dimensions();
if (LDim[0] == 0) {
std::cerr << "VarNature = " << RecS.VarNature << "\n";
std::cerr << "Variable eVarName = " << eVarName << "\n";
std::cerr << "is recognized by the program\n";
std::cerr << "But it has not been assigned.\n";
std::cerr << "Possibly because of missing facility for\n";
std::cerr << "eModelName = " << eModelName << "\n";
throw TerminalException{1};
}
int dim0=LDim[0];
int dim1=LDim[1];
int dim2=LDim[2];
Eigen::Tensor<double,3> Fwr(dim0, dim1, dim2);
for (int i0=0; i0<dim0; i0++)
for (int i1=0; i1<dim1; i1++)
for (int i2=0; i2<dim2; i2++) {
double eU=Uthree(i0, i1, i2);
double eV=Vthree(i0, i1, i2);
double eNorm=sqrt(eU*eU + eV*eV);
Fwr(i0, i1, i2) = eNorm;
}
eRecVar.Tens3=Fwr;
}
}
return eRecVar;
}
RecVar ModelSpecificVarSpecificTime(TotalArrGetData const& TotalArr, std::string const& eVarName, double const& eTimeDay)
{
std::string eSep="_";
std::vector<std::string> ListStr=STRING_Split(eVarName, eSep);
int len=ListStr.size();
if (len == 1)
return ModelSpecificVarSpecificTime_Kernel(TotalArr, eVarName, eTimeDay);
std::string eVar_rho=ListStr[0];
std::string eVar_uv=ListStr[1];
RecVar RecVar_rho=ModelSpecificVarSpecificTime_Kernel(TotalArr, eVar_rho, eTimeDay);
RecVar RecVar_uv =ModelSpecificVarSpecificTime_Kernel(TotalArr, eVar_uv , eTimeDay);
std::string VarNat_rho=RecVar_rho.RecS.VarNature;
std::string VarNat_uv =RecVar_uv.RecS.VarNature;
if (VarNat_rho != "rho" && VarNat_rho != "3Drho") {
std::cerr << "The RecVar_rho is not a rho type variable. Error!\n";
std::cerr << "Correct way to call is Var_rho _ Var_uv\n";
std::cerr << "for Example Hwave_SurfCurr for Hwave as rho variable and SurfCurr as uv variable\n";
std::cerr << "The call was with eVarName=" << eVarName << "\n";
throw TerminalException{1};
}
if (VarNat_uv != "uv" && VarNat_uv != "3Duv") {
std::cerr << "The RecVar_uv is not a uv type variable. Error!\n";
std::cerr << "Correct way to call is Var_rho _ Var_uv\n";
std::cerr << "for Example Hwave_SurfCurr for Hwave as rho variable and SurfCurr as uv variable\n";
std::cerr << "The call was with eVarName=" << eVarName << "\n";
throw TerminalException{1};
}
if ((VarNat_rho == "3Drho" && VarNat_uv == "uv") || (VarNat_rho == "rho" && VarNat_uv == "3Duv") ) {
std::cerr << "Error. variables do not have the same dimensionality\n";
std::cerr << "It should be both 3D or both 2D\n";
std::cerr << "Right now, we have VarNat_rho=" << VarNat_rho << "\n";
std::cerr << "Right now, we have VarNat_uv =" << VarNat_uv << "\n";
throw TerminalException{1};
}
if (VarNat_rho == "rho") {
RecVar_rho.U = RecVar_uv.U;
RecVar_rho.V = RecVar_uv.V;
}
if (VarNat_rho == "3Drho") {
RecVar_rho.Uthree = RecVar_uv.Uthree;
RecVar_rho.Vthree = RecVar_uv.Vthree;
}
RecVar_rho.RecS.VarName1 += "_" + RecVar_uv.RecS.VarName1;
RecVar_rho.RecS.VarName2 += " + " + RecVar_uv.RecS.VarName2;
RecVar_rho.RecS.VarNature="uv";
return RecVar_rho;
}
RecVar RetrieveTrivialRecVar(std::string const& eVarName)
{
TotalArrGetData TotalArrTrivial;
TotalArrTrivial.GrdArr.ModelName="TRIVIAL";
double eTimeDayTrivial=0;
return ModelSpecificVarSpecificTime(TotalArrTrivial, eVarName, eTimeDayTrivial);
}
std::vector<std::string> GetAllPossibleVariables_with_pairs()
{
std::vector<std::string> ListVar=GetAllPossibleVariables();
std::vector<std::string> ListVar_rho;
std::vector<std::string> ListVar_uv;
std::vector<std::string> ListVar_3Drho;
std::vector<std::string> ListVar_3Duv;
std::vector<std::string> ListVar_Ret;
for (auto & eVar : ListVar) {
RecVar eRec=RetrieveTrivialRecVar(eVar);
if (eRec.RecS.VarNature == "rho")
ListVar_rho.push_back(eVar);
if (eRec.RecS.VarNature == "uv")
ListVar_uv.push_back(eVar);
if (eRec.RecS.VarNature == "3Drho")
ListVar_3Drho.push_back(eVar);
if (eRec.RecS.VarNature == "3Duv")
ListVar_3Duv.push_back(eVar);
ListVar_Ret.push_back(eVar);
}
for (auto & eVar_uv : ListVar_uv) {
for (auto & eVar_rho : ListVar_rho) {
std::string eVarTot = eVar_rho + "_" + eVar_uv;
ListVar_Ret.push_back(eVarTot);
}
}
for (auto & eVar_uv : ListVar_3Duv) {
for (auto & eVar_rho : ListVar_3Drho) {
std::string eVarTot = eVar_rho + "_" + eVar_uv;
ListVar_Ret.push_back(eVarTot);
}
}
return ListVar_Ret;
}
RecVar ModelSpecificVarSpecificTimeBound(TotalArrGetData const& TotalArr, std::string const& eVarName, double const& eTimeDay, PlotBound const& ePlotBound)
{
RecVar eRecVar=ModelSpecificVarSpecificTime(TotalArr, eVarName, eTimeDay);
ApplyPlotBound(TotalArr, eRecVar, eVarName, ePlotBound);
return eRecVar;
}
RecVar ModelSpecificVarSpecificTimeGeneral(TotalArrGetData const& TotalArr, std::string const& eVarName, VarQuery const& eQuery, PlotBound const& ePlotBound)
{
std::vector<std::string> ListAllow{"instant", "average", "swathMax", "swathMin"};
if (std::find(ListAllow.begin(), ListAllow.end(), eQuery.NatureQuery) == ListAllow.end()) {
std::cerr << "We failed to find NatureQuery=" << eQuery.NatureQuery << "\n";
std::cerr << "List of allowed queries:\n";
for (auto & eStr : ListAllow)
std::cerr << "  eStr=" << eStr << "\n";
throw TerminalException{1};
}
std::string strPres=DATE_ConvertMjd2mystringPres(eQuery.eTimeDay);
std::cerr << "Query ModelSpecificVarSpecificTimeGeneral NatureQuery=" << eQuery.NatureQuery << " date=" << strPres << " VarName=" << eVarName << "\n";
RecVar eRecVar;
if (eQuery.NatureQuery == "instant") {
eRecVar=ModelSpecificVarSpecificTimeBound(TotalArr, eVarName, eQuery.eTimeDay, ePlotBound);
}
else {
double eTimeDay=eQuery.eTimeDay;
double TimeFrameDay=eQuery.TimeFrameDay;
std::vector<int> ListRelITime=GetIntervalListITime(TotalArr.eArr.ListTime, eTimeDay, TimeFrameDay);
int nbTimeRel=ListRelITime.size();
RecVar RecVarTrivial=RetrieveTrivialRecVar(eVarName);
MyMatrix<double> F, U, V;
Eigen::Tensor<double,3> Tens3;
Eigen::Tensor<double,3> Uthree;
Eigen::Tensor<double,3> Vthree;
for (int iTimeRel=0; iTimeRel<nbTimeRel; iTimeRel++) {
int iTime=ListRelITime[iTimeRel];
double eTimeDayB=TotalArr.eArr.ListTime[iTime];
eRecVar=ModelSpecificVarSpecificTimeBound(TotalArr, eVarName, eTimeDayB, ePlotBound);
if (iTimeRel == 0) {
if (RecVarTrivial.RecS.VarNature == "rho") {
F=eRecVar.F;
}
if (RecVarTrivial.RecS.VarNature == "uv") {
U=eRecVar.U;
V=eRecVar.V;
F=eRecVar.F;
}
if (RecVarTrivial.RecS.VarNature == "3Drho") {
Tens3=eRecVar.Tens3;
}
if (RecVarTrivial.RecS.VarNature == "3Duv") {
Uthree=eRecVar.Uthree;
Vthree=eRecVar.Vthree;
Tens3=eRecVar.Tens3;
}
}
else {
if (eQuery.NatureQuery == "average") {
if (RecVarTrivial.RecS.VarNature == "rho") {
F += eRecVar.F;
}
if (RecVarTrivial.RecS.VarNature == "uv") {
U += eRecVar.U;
V += eRecVar.V;
F += eRecVar.F;
}
if (RecVarTrivial.RecS.VarNature == "3Drho") {
Tens3 += eRecVar.Tens3;
}
if (RecVarTrivial.RecS.VarNature == "3Drho") {
Uthree += eRecVar.Uthree;
Vthree += eRecVar.Vthree;
Tens3 += eRecVar.Tens3;
}
}
if (eQuery.NatureQuery == "swathMax") {
if (RecVarTrivial.RecS.VarNature == "rho") {
F=F.cwiseMax(eRecVar.F);
}
if (RecVarTrivial.RecS.VarNature == "3Drho") {
Tens3=Tens3.cwiseMax(eRecVar.Tens3);
}
if (RecVarTrivial.RecS.VarNature == "uv" || RecVarTrivial.RecS.VarNature == "3Duv") {
std::cerr << "swathMax for uv does not have any sense\n";
throw TerminalException{1};
}
}
if (eQuery.NatureQuery == "swathMin") {
if (RecVarTrivial.RecS.VarNature == "rho") {
F=F.cwiseMin(eRecVar.F);
}
if (RecVarTrivial.RecS.VarNature == "3Drho") {
Tens3=Tens3.cwiseMin(eRecVar.Tens3);
}
if (RecVarTrivial.RecS.VarNature == "uv" || RecVarTrivial.RecS.VarNature == "3Duv") {
std::cerr << "swathMin for uv does not have any sense\n";
throw TerminalException{1};
}
}
}
}
if (eQuery.NatureQuery == "average") {
if (RecVarTrivial.RecS.VarNature == "rho") {
F /= double(nbTimeRel);
}
if (RecVarTrivial.RecS.VarNature == "uv") {
U /= double(nbTimeRel);
V /= double(nbTimeRel);
F /= double(nbTimeRel);
}
if (RecVarTrivial.RecS.VarNature == "3Drho") {
auto LDim=Tens3.dimensions();
int dim0=LDim[0];
int dim1=LDim[1];
int dim2=LDim[2];
for (int i0=0; i0<dim0; i0++)
for (int i1=0; i1<dim1; i1++)
for (int i2=0; i2<dim2; i2++)
Tens3(i0, i1, i2) /= double(nbTimeRel);
}
if (RecVarTrivial.RecS.VarNature == "3Drho") {
auto LDim=Tens3.dimensions();
int dim0=LDim[0];
int dim1=LDim[1];
int dim2=LDim[2];
for (int i0=0; i0<dim0; i0++)
for (int i1=0; i1<dim1; i1++)
for (int i2=0; i2<dim2; i2++) {
Uthree(i0, i1, i2) /= double(nbTimeRel);
Vthree(i0, i1, i2) /= double(nbTimeRel);
Tens3 (i0, i1, i2) /= double(nbTimeRel);
}
}
}
if (RecVarTrivial.RecS.VarNature == "rho") {
eRecVar.F=F;
}
if (RecVarTrivial.RecS.VarNature == "3Drho") {
eRecVar.Tens3=Tens3;
}
if (RecVarTrivial.RecS.VarNature == "uv") {
eRecVar.U=U;
eRecVar.V=V;
eRecVar.F=F;
}
if (RecVarTrivial.RecS.VarNature == "3Duv") {
eRecVar.Uthree=Uthree;
eRecVar.Vthree=Vthree;
eRecVar.Tens3=Tens3;
}
}
std::string strAll=GetStrAllOfPlot(eQuery);
std::string strPresPlot=GetStrPresOfPlot(eQuery);
eRecVar.RecS.strAll=strAll;
eRecVar.RecS.strPres=strPresPlot;
ApplyPlotBound(TotalArr, eRecVar, eVarName, ePlotBound);
return eRecVar;
}
void PLOT_SCATTER(DrawScatterArr const& eDrawScatter,
NCLcaller & eCall,
PermanentInfoDrawing const& ePerm)
{
std::string TitleStr;
std::string VarNameAB_file=eDrawScatter.VarNameAB_file;
std::string eFileNC=ePerm.PrefixTemp.str() + "DataScatter_" + VarNameAB_file + ".nc";
std::string eFileNCL=ePerm.PrefixTemp.str() + "ScriptScatter_" + VarNameAB_file + ".ncl";
std::string FileName=ePerm.eDir + VarNameAB_file;
std::string TargetFile=FileName + "." + ePerm.Extension;
DEFINE_SCATTER_NC(eFileNC, eDrawScatter);
std::ofstream OUTncl;
OUTncl.open(eFileNCL);
OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
PrintMyScriptSubtitle(OUTncl);
OUTncl << "begin\n";
OUTncl << "  f = addfile(\"" << eFileNC << "\", \"r\")\n";
OUTncl << "  X = f->X2D\n";
OUTncl << "  Y = f->Y2D\n";
OUTncl << "  Vtot = f->canvas\n";
OUTncl << "  stride=20\n";
OUTncl << "  spacing=0.025\n";
OUTncl << "  minValue=min(Vtot)\n";
OUTncl << "  maxValue=2\n";
OUTncl << "  wks  = gsn_open_wks (\"" << ePerm.Extension << "\",\"" << FileName << "\")\n";
OUTncl << "  gsn_define_colormap(wks,\"BlAqGrYeOrRevi200\")\n";
OUTncl << "  vres1 = True               ; plot mods desired\n";
OUTncl << "  vres1@gsnDraw   = False\n";
OUTncl << "  vres1@gsnFrame  = False\n";
OUTncl << "  vres1@gsnMaximize     = True    ; Maximize plot in frame\n";
OUTncl << "  vres1@gsnPaperOrientation = \"Portrait\"\n";
OUTncl << "  vres1@cnFillDrawOrder        = \"PreDraw\"\n";
OUTncl << "  vres1@cnFillOn             = True               ; turn on color for contours\n";
OUTncl << "  vres1@cnLinesOn            = False              ; turn off contour lines\n";
OUTncl << "  vres1@cnLineLabelsOn       = False              ; turn off contour line labels\n";
OUTncl << "  vres1@cnFillMode           = \"AreaFill\"\n";
OUTncl << "  vres1@gsnSpreadColors      = True               ; use full color map\n";
OUTncl << "  vres1@gsnSpreadColorEnd     = -2\n";
OUTncl << "  LLabel=new(81,string)\n";
OUTncl << "  LLabel(0)=\"1\"\n";
OUTncl << "  LLabel(20)=\"3\"\n";
OUTncl << "  LLabel(40)=\"10\"\n";
OUTncl << "  LLabel(60)=\"30\"\n";
OUTncl << "  LLabel(80)=\"100\"\n";
OUTncl << "  vres1@lbLabelBarOn =  True\n";
OUTncl << "  vres1@lbLabelStride            = stride\n";
OUTncl << "  vres1@lbOrientation        = \"Vertical\"     ; Vertical label bar\n";
OUTncl << "  vres1@lbLabelStrings = LLabel\n";
OUTncl << "  vres1@cnLevelSelectionMode = \"ManualLevels\"     ; set manual contour levels\n";
OUTncl << "  vres1@cnMinLevelValF       = minValue                ; set min contour level\n";
OUTncl << "  vres1@cnMaxLevelValF       = maxValue              ; set max contour level\n";
OUTncl << "  vres1@cnLevelSpacingF      = spacing                 ; set contour spacing\n";
OUTncl << "  vres1@pmLabelBarOrthogonalPosF = -0.01          ; move label bar closer\n";
OUTncl << "  ;  vres1@pmLabelBarDisplayMode = \"Always\"          ; Turn on a label bar.\n";
OUTncl << "  vres1@lbPerimOn             = False             ; no box around it\n";
OUTncl << "  vres1@lbBoxLinesOn         = False               ; Yes/No labelbar box lines.\n";
OUTncl << "  vres1@tiXAxisString  = \"" << eDrawScatter.NameA_plot << "\"\n";
OUTncl << "  vres1@tiYAxisString    = \"" << eDrawScatter.NameB_plot << "\"\n";
OUTncl << "  vres1@tiXAxisOffsetYF = 0.0\n";
OUTncl << "  ; vres1@tmYLPrecision = 0\n";
OUTncl << "  ; First part Hwave\n";
OUTncl << "  vres1@sfXArray            = X\n";
OUTncl << "  vres1@sfYArray            = Y\n";
OUTncl << "  vres1@trGridType          = \"TriangularMesh\"\n";
OUTncl << "  ;  vres1@trGridType         = \"curvilinear\"\n";
OUTncl << "  plot = gsn_csm_contour(wks,Vtot,vres1)\n";
OUTncl << "  resP = True\n";
OUTncl << "  resP@gsLineThicknessF = 1.5\n";
OUTncl << "  resP@gsLineColor  = \"black\"\n";
OUTncl << "  resP@tfPolyDrawOrder = \"PostDraw\"\n";
if (eDrawScatter.AddStatMeasModel) {
OUTncl << "  xp=(/-180, 180/)\n";
OUTncl << "  yp=(/-180, 180/)\n";
OUTncl << "  line0=gsn_add_polyline(wks,plot,xp,yp,resP)\n";
}
OUTncl << "  data_rangeA = f->data_rangeA\n";
OUTncl << "  data_rangeB = f->data_rangeB\n";
OUTncl << "  line1=gsn_add_polyline(wks,plot,data_rangeA,data_rangeB,resP)\n";
ADD_ANNOTATION_TEXT(OUTncl, eDrawScatter.TheAnnot);
if (eDrawScatter.AddStatMeasModel) {
T_stat eStat=ComputeStatistics_MyVector(eDrawScatter.eVectA, eDrawScatter.eVectB);
std::string strWrite="m=" + eStat.strSlope + " c=" + eStat.strCorrelation + " s=" + eStat.strScatterIndex;
OUTncl << "  txresB             = True\n";
OUTncl << "  txresB@txFontHeightF = 0.02\n";
OUTncl << "  txresB@txFontColor = \"black\"\n";
OUTncl << "  strLeft=\"\"\n";
OUTncl << "  strMid=\"" << strWrite << "\"\n";
OUTncl << "  strRight=\"\"\n";
OUTncl << "  subtitles(wks, plot, strLeft, strMid, strRight, txresB)\n";
}
OUTncl << "  draw(plot)\n";
OUTncl << "  frame(wks)\n";
OUTncl << "end\n";
OUTncl.close();
eCall.SubmitJob(TargetFile, eFileNC, eFileNCL);
}
void PLOT_PCOLOR(std::string const& FileName,
GridArray const& GrdArr,
DrawArr const& eDrawArr,
RecVar const& eRecVar,
NCLcaller & eCall,
PermanentInfoDrawing const& ePerm)
{
RecSymbolic RecS=eRecVar.RecS;
std::string eFileNC=ePerm.PrefixTemp.str() + "DataPcolor_" + eDrawArr.VarNameUF + "_" + RecS.strAll + ".nc";
std::string eFileNCL=ePerm.PrefixTemp.str() + "ScriptPcolor_" + eDrawArr.VarNameUF + "_" + RecS.strAll + ".ncl";
DEFINE_PCOLOR_NC_NCL(eFileNC, GrdArr, eRecVar.F,
eDrawArr.DrawContourBathy,
eDrawArr.ListLineSegment);
int IsFE=GrdArr.IsFE;
bool IsSpherical=GrdArr.IsSpherical;
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
OUTncl << "  wks  = gsn_open_wks (\"" << ePerm.Extension << "\",\"" << FileName << "\")\n";
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
OUTncl << "  res2@gsnPaperOrientation = \"Landscape\"\n";
OUTncl << "  ;\n";
OUTncl << "  ; General frame information\n";
OUTncl << "  ;\n";
if (IsSpherical) {
OUTncl << "  res2@mpProjection = \"Mercator\"\n";
OUTncl << "  res2@mpLimitMode         = \"Corners\"             ; choose range of map\n";
OUTncl << "  res2@mpLeftCornerLatF    = " << eDrawArr.eQuadFrame.MinLat << "\n";
OUTncl << "  res2@mpLeftCornerLonF    = " << eDrawArr.eQuadFrame.MinLon << "\n";
OUTncl << "  res2@mpRightCornerLatF   = " << eDrawArr.eQuadFrame.MaxLat << "\n";
OUTncl << "  res2@mpRightCornerLonF   = " << eDrawArr.eQuadFrame.MaxLon << "\n";
}
OUTncl << "  res2@pmTickMarkDisplayMode  = \"Always\"           ; turn on tickmarks\n";
if (IsSpherical) {
if (eDrawArr.FillLand) {
OUTncl << "  res2@mpFillOn      = True\n";
OUTncl << "  res2@mpDataBaseVersion      = \"" << eDrawArr.GridResolution << "\"          ; use high resolution coast\n";
OUTncl << "  res2@mpLandFillColor      = \"LightGrey\"\n";
OUTncl << "  res2@mpLandFillColor       = \"gray\"            ; set land to be gray\n";
}
else {
OUTncl << "  res2@mpFillOn      = False\n";
}
}
OUTncl << "  ;\n";
OUTncl << "  ; Contour map information\n";
OUTncl << "  ;\n";
OUTncl << "  res2@cnFillDrawOrder        = \"PreDraw\"\n";
OUTncl << "  res2@cnFillOn             = True               ; turn on color for contours\n";
OUTncl << "  res2@cnLinesOn            = False              ; turn off contour lines\n";
OUTncl << "  res2@cnLineLabelsOn       = False              ; turn off contour line labels\n";
if (IsFE == 1 && eDrawArr.cnFillMode == "CellFill") {
std::cerr << "The \"CellFill\" option can only be used in finite difference grids\n";
throw TerminalException{1};
}
if (IsFE == 1 && eDrawArr.cnFillMode == "AreaFill") {
OUTncl << "; This is an option to emulate \"AreaFill\" method\n";
OUTncl << "  res2@cnFillMode           = \"RasterFill\"\n";
OUTncl << "  res2@cnRasterSmoothingOn = True\n";
}
else {
OUTncl << "  res2@cnFillMode           = \"" << eDrawArr.cnFillMode << "\"\n";
}
OUTncl << "            ; AreaFill : slow and buggy but maybe more beautiful\n";
OUTncl << "            ; RasterFill : fast and efficient\n";
OUTncl << "            ; CellFill : similar to RasterFill but only for finite difference\n";
OUTncl << "  ;  res2@cnRasterSmoothingOn  = True\n";
OUTncl << "  res2@cnSmoothingOn = " << NCL_bool(eDrawArr.cnSmoothingOn) << "\n";
OUTncl << "  ;  res2@cnSmoothingDistanceF  = 0.05\n";
OUTncl << "  res2@cnSmoothingTensionF  = -1\n";
OUTncl << "  res2@cnLevelSelectionMode = \"ManualLevels\"     ; set manual contour levels\n";
OUTncl << "  res2@cnMinLevelValF       = " << RecS.minval << "  ; set min contour level\n";
OUTncl << "  res2@cnMaxLevelValF       = " << RecS.maxval << "  ; set max contour level\n";
int nbLevelSpa=eDrawArr.nbLevelSpa;
double TheLevelSpa=(RecS.maxval - RecS.minval)/double(nbLevelSpa);
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
OUTncl << "  ;\n";
OUTncl << "  ; Label bar plotting\n";
OUTncl << "  ;\n";
if (eDrawArr.DoColorBar == true && IsSpherical == true) {
OUTncl << "  res2@lbLabelBarOn = True\n";
}
else {
OUTncl << "  res2@lbLabelBarOn = False\n";
}
OUTncl << "  res2@lbTitleString    = \"" << RecS.VarName1 << " [" << RecS.Unit << "]\"\n";
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
if (IsSpherical == true) {
OUTncl << "  res2@pmLabelBarWidthF = 0.03\n";
}
OUTncl << "  ; res2@gsnRightString  = \"Sea surface elevation\"\n";
OUTncl << "  ; res2@gsnLeftString    = \"Difference\"\n";
OUTncl << "  ; res2@gsnRightString  = \"\"\n";
OUTncl << "  ; res2@gsnLeftString    = \"\"\n";
OUTncl << "  ;\n";
OUTncl << "  ; Colormap assignation\n";
OUTncl << "  ;\n";
OUTncl << "  gsn_define_colormap (wks,\"" << eDrawArr.ColorMap << "\")\n";
OUTncl << "  ;     other possibilities: hotres, rainbow, ViBlGrWhYeOrRe, BlWhRe, GrayWhiteGray, BlGrYeOrReVi200\n";
OUTncl << "  i = NhlNewColor(wks,0.8,0.8,0.8)      ; add gray to colormap\n";
OUTncl << "  i = NhlNewColor(wks,0.9,0.9,0.9)      ; add gray to colormap\n";
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
if (IsSpherical) {
OUTncl << "  plot = gsn_csm_contour_map(wks,eVarF,res2)\n";
}
else {
OUTncl << "  plot = gsn_csm_contour(wks,eVarF,res2)\n";
}
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
std::string TargetFile=FileName + "." + ePerm.Extension;
eCall.SubmitJob(TargetFile, eFileNC, eFileNCL);
}
void LINES_PLOT(std::string const& FileName,
DrawLinesArr const& eDrawArr,
NCLcaller & eCall,
PermanentInfoDrawing const& ePerm)
{
std::string eFileNC=ePerm.PrefixTemp.str() + "DataLines_" + eDrawArr.VarName + ".nc";
std::string eFileNCL=ePerm.PrefixTemp.str() + "ScriptLines_" + eDrawArr.VarName + ".ncl";
LINES_DEFINE_NC(eFileNC, eDrawArr);
int nbArr=eDrawArr.ListListVect.size();
std::ofstream OUTncl;
OUTncl.open(eFileNCL);
OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl\"\n";
OUTncl << "load \"$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl\"\n";
PrintMyScriptSubtitle(OUTncl);
OUTncl << ";************************************************\n";
OUTncl << "begin\n";
OUTncl << "  ;\n";
OUTncl << "  ; Data reading\n";
OUTncl << "  ;\n";
OUTncl << "  f = addfile(\"" << eFileNC << "\", \"r\")\n";
OUTncl << "  ListX=f->ListX\n";
OUTncl << "  Data=f->ListListVect\n";
OUTncl << "  TheMax=max(Data)\n";
OUTncl << "  TheMin=0\n";
OUTncl << "  wks  = gsn_open_wks (\"" << ePerm.Extension << "\",\"" << FileName << "\")\n";
OUTncl << "  res                  = True   ; plot mods desired\n";
OUTncl << "  res@gsnFrame          = False\n";
OUTncl << "  res@xyMonoDashPattern = False\n";
OUTncl << "  res@xyDashPatterns = (/0/)\n";
OUTncl << "  res@tiMainFontHeightF  = 0.015\n";
int nbCharX=eDrawArr.XAxisString.size();
if (nbCharX > 0) {
OUTncl << "  res@tiXAxisString   = \"x (m)\"\n";
OUTncl << "  res@tiXAxisFontHeightF = 0.020\n";
}
int nbCharY=eDrawArr.YAxisString.size();
if (nbCharY > 0) {
OUTncl << "  res@tiYAxisString   = \"" << eDrawArr.YAxisString << "\"\n";
OUTncl << "  res@tiYAxisFontHeightF = 0.015\n";
}
std::vector<std::string> ListColors={"black", "red", "blue", "purple", "green"};
int nbColor=ListColors.size();
if (nbArr > nbColor) {
std::cerr << "Number of colors is insufficient\n";
throw TerminalException{1};
}
std::string TheTotal;
for (int iArr=0; iArr<nbArr; iArr++) {
if (iArr > 0)
TheTotal=TheTotal + ",";
TheTotal=TheTotal + "\"" + ListColors[iArr] + "\"";
}
OUTncl << "  res@xyLineColors = (/" << TheTotal << "/)\n";
OUTncl << "  res@trYMaxF=" << eDrawArr.TheMax << "\n";
OUTncl << "  res@trYMinF=" << eDrawArr.TheMin << "\n";
OUTncl << "  res@trXMaxF=max(ListX)\n";
OUTncl << "  res@trXMinF=min(ListX)\n";
bool DoExplicit=false;
if (DoExplicit) {
OUTncl << "  res@tmXBMode=\"Explicit\" \n";
OUTncl << "  res@tmXBValues=(/0,10,20,max(ListTimeDayStart)/)\n";
OUTncl << "  res@tmXBLabels=(/\"15 nov\", \"25 nov\", \"5 dec\", \"15 dec\"/)\n";
}
OUTncl << "  plot=gsn_csm_xy(wks,ListX,Data,res)\n";
if (eDrawArr.DoTitle) {
OUTncl << "  txresB             = True\n";
OUTncl << "  txresB@txFontHeightF = 0.02\n";
OUTncl << "  txresB@txFontColor = \"black\"\n";
OUTncl << "  strLeft=\"\"\n";
OUTncl << "  strMid=\"" << eDrawArr.TitleStr << "\"\n";
OUTncl << "  strRight=\"\"\n";
OUTncl << "  subtitles(wks, plot, strLeft, strMid, strRight, txresB)\n";
}
bool DoBackground=true;
if (DoBackground) {
OUTncl << "  lgres                    = True\n";
OUTncl << "  lgres@lgLineColors     = (/" << TheTotal << "/)\n";
OUTncl << "  lgres@lgItemType         = \"Lines\"\n";
OUTncl << "  lgres@lgDashIndexes = (/";
for (int iArr=0; iArr<nbArr; iArr++) {
if (iArr > 0)
OUTncl << ",";
OUTncl << "0";
}
OUTncl << "/)\n";
OUTncl << "  lgres@lgLabelFontHeightF = .07\n";
OUTncl << "  lgres@vpWidthF           = 0.11\n";
OUTncl << "  lgres@vpHeightF          = 0.12\n";
OUTncl << "  lgres@lgPerimOn = False\n";
OUTncl << "  ;   lgres@lgPerimColor       = \"orange\"\n";
OUTncl << "  lgres@lgPerimThicknessF  = 5.0\n";
OUTncl << "  ListLabels= (/";
for (int iArr=0; iArr<nbArr; iArr++) {
if (iArr > 0)
OUTncl << ",";
OUTncl << "\"" << eDrawArr.ListName_plot[iArr] << "\"";
}
OUTncl << "/)\n";
OUTncl << "  lbid = gsn_create_legend(wks," << nbArr << ",ListLabels,lgres)\n";
OUTncl << "  amres = True\n";
OUTncl << "  amres@amParallelPosF   = 0.35\n";
OUTncl << "  amres@amOrthogonalPosF = 0.30\n";
OUTncl << "  annoid1 = gsn_add_annotation(plot,lbid,amres)\n";
}
ADD_ANNOTATION_TEXT(OUTncl, eDrawArr.TheAnnot);
OUTncl << "  draw(plot)\n";
OUTncl << "  frame(wks)\n";
OUTncl << "end\n";
OUTncl.close();
std::string TargetFile=FileName + "." + ePerm.Extension;
eCall.SubmitJob(TargetFile, eFileNC, eFileNCL);
}
FullNamelist NAMELIST_GetStandard_PlotRoutine_common()
{
std::map<std::string, SingleBlock> ListBlock;
std::string BlockName1="PROC";
std::map<std::string, int> ListIntValues1;
std::map<std::string, bool> ListBoolValues1;
std::map<std::string, double> ListDoubleValues1;
std::map<std::string, std::string> ListStringValues1;
std::map<std::string, std::vector<std::string> > ListListStringValues1;
std::string LPoss="Possibilities:";
bool IsFirst=true;
for (auto & eStr : GetAllPossibleModels()) {
if (IsFirst == false)
LPoss += ",";
LPoss += " " + eStr;
}
ListStringValues1["MODELNAME"]=LPoss;
ListStringValues1["BEGTC"]="20110915.000000";
ListStringValues1["ENDTC"]="20110925.000000";
ListDoubleValues1["DELTC"]=600;
ListStringValues1["UNITC"]="SEC";
ListStringValues1["GridFile"]="unset GridFile";
ListStringValues1["BoundFile"]="unset";
ListBoolValues1["CutWorldMap"]=false;
ListBoolValues1["HigherLatitudeCut"]=false;
ListBoolValues1["SplittingAt180"]=false;
ListDoubleValues1["MinLatCut"]=-80;
ListDoubleValues1["MaxLatCut"]=80;
ListStringValues1["PicPrefix"]="Pictures/DIR_plot/";
ListStringValues1["Extension"]="png";
ListListStringValues1["ListNatureQuery"]={"instant"};
ListDoubleValues1["TimeFrameDay"]=1;
ListBoolValues1["FirstCleanDirectory"]=true;
ListBoolValues1["KeepNC_NCL"]=false;
ListBoolValues1["OverwritePrevious"]=false;
ListBoolValues1["WriteITimeInFileName"]=true;
ListIntValues1["NPROC"]=1;
ListDoubleValues1["ThresholdApplyAverage"]=10000;
ListBoolValues1["ApplyThresholdAveraging"]=false;
SingleBlock BlockPROC;
BlockPROC.ListIntValues=ListIntValues1;
BlockPROC.ListBoolValues=ListBoolValues1;
BlockPROC.ListDoubleValues=ListDoubleValues1;
BlockPROC.ListStringValues=ListStringValues1;
BlockPROC.ListListStringValues=ListListStringValues1;
BlockPROC.BlockName=BlockName1;
ListBlock["PROC"]=BlockPROC;
std::string BlockName2="PLOT";
std::map<std::string, int> ListIntValues2;
std::map<std::string, bool> ListBoolValues2;
std::map<std::string, double> ListDoubleValues2;
std::map<std::string, std::string> ListStringValues2;
std::map<std::string, std::vector<double> > ListListDoubleValues2;
std::map<std::string, std::vector<int> > ListListIntValues2;
std::map<std::string, std::vector<std::string> > ListListStringValues2;
ListStringValues2["ColorMap"]="BlAqGrYeOrReVi200";
ListStringValues2["ColorMapDiff"]="BlWhRe";
ListStringValues2["cnFillMode"]="RasterFill";
ListBoolValues2["DoColorBar"]=true;
ListBoolValues2["cnSmoothingOn"]=true;
ListIntValues2["nbLevelSpa"]=50;
ListIntValues2["nbLabelStride"]=10;
ListBoolValues2["UseNativeGrid"]=true;
ListBoolValues2["DoTitle"]=true;
ListStringValues2["GridResolution"]="HighRes";
ListBoolValues2["DrawRiver"]=false;
ListBoolValues2["PrintMMA"]=false;
ListBoolValues2["LocateMM"]=false;
ListBoolValues2["DoMain"]=true;
ListBoolValues2["PlotDepth"]=true;
ListBoolValues2["PlotMesh"]=false;
ListBoolValues2["DrawContourBathy"]=false;
ListBoolValues2["DrawAnnotation"]=false;
ListBoolValues2["ExcludeLargeValues"]=false;
ListDoubleValues2["ThresholdExclusionPlot"]=100000;
ListDoubleValues2["MultiplierResolutionFE_FD"]=1;
ListListIntValues2["Tens3ListLevel"]={};
ListDoubleValues2["AnnotationLon"]=0;
ListDoubleValues2["AnnotationLat"]=0;
ListStringValues2["AnnotationText"]="something to write";
ListListStringValues2["RenameVariable_VarName1"]={};
ListListStringValues2["RenameVariable_VarName2"]={};
ListListStringValues2["BoundSingle_var"]={};
ListListDoubleValues2["BoundSingle_min"]={};
ListListDoubleValues2["BoundSingle_max"]={};
ListListStringValues2["BoundDiff_var"]={};
ListListDoubleValues2["BoundDiff_min"]={};
ListListDoubleValues2["BoundDiff_max"]={};
ListBoolValues2["VariableRange"]=false;
ListBoolValues2["FillLand"]=true;
ListListDoubleValues2["ListFrameMinLon"]={};
ListListDoubleValues2["ListFrameMinLat"]={};
ListListDoubleValues2["ListFrameMaxLon"]={};
ListListDoubleValues2["ListFrameMaxLat"]={};
ListBoolValues2["UseFDgrid"]=false;
ListBoolValues2["DoMain"]=true;
SingleBlock BlockPLOT;
BlockPLOT.ListIntValues=ListIntValues2;
BlockPLOT.ListBoolValues=ListBoolValues2;
BlockPLOT.ListDoubleValues=ListDoubleValues2;
BlockPLOT.ListStringValues=ListStringValues2;
BlockPLOT.ListListStringValues=ListListStringValues2;
BlockPLOT.ListListDoubleValues=ListListDoubleValues2;
BlockPLOT.ListListIntValues=ListListIntValues2;
BlockPLOT.BlockName=BlockName2;
ListBlock["PLOT"]=BlockPLOT;
std::string BlockName3="VARS";
std::map<std::string, int> ListIntValues3;
std::map<std::string, bool> ListBoolValues3;
std::map<std::string, double> ListDoubleValues3;
std::map<std::string, std::string> ListStringValues3;
std::map<std::string, std::vector<std::string> > ListListStringValues3;
std::vector<std::string> ListVarOut=GetAllPossibleVariables_with_pairs();
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
FullNamelist eFullNamelist;
eFullNamelist.ListBlock=ListBlock;
eFullNamelist.FileName="undefined";
return eFullNamelist;
}
void InterpolateAltimeterData(std::vector<SingleEntryAltimeter> & ListEntry, GridArray const& GrdArr, FullNamelist const& eFull)
{
SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
std::string HisPrefix=eBlPROC.ListStringValues.at("HisPrefix");
std::string eModelName=eBlPROC.ListStringValues.at("MODELNAME");
std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, false, false, 0, 0};
ArrayHistory eArr=ReadArrayHistory(eTriple);
TotalArrGetData TotalArr{GrdArr, eArr};
size_t nbEntry=ListEntry.size();
MyMatrix<double> ListXY(2, nbEntry);
std::vector<double> ListTime(nbEntry);
std::vector<int> ListStatusTime(nbEntry,0);
for (size_t iEntry=0; iEntry<nbEntry; iEntry++) {
ListXY(0,iEntry)=ListEntry[iEntry].Lon;
ListXY(1,iEntry)=ListEntry[iEntry].Lat;
}
int nbTime=eArr.ListTime.size();
std::vector<SingleRecInterp> LRec=General_FindInterpolationWeight(GrdArr, ListXY);
std::vector<std::vector<SingleSearchEntry> > ListListCases(nbTime);
for (size_t iEntry=0; iEntry<nbEntry; iEntry++) {
double eTime=ListEntry[iEntry].Time;
ListTime[iEntry]=eTime;
InterpInfo eInterpInfo=GetTimeInterpolationInfo(eArr.ListTime, eTime);
SingleRecInterp eRec=LRec[iEntry];
int iTimeLow=eInterpInfo.iTimeLow;
int iTimeUpp=eInterpInfo.iTimeUpp;
double alphaLow=eInterpInfo.alphaLow;
double alphaUpp=eInterpInfo.alphaUpp;
double eValAssign;
if (eRec.status == 1) {
eValAssign=0;
for (auto & ePart : eRec.LPart) {
int eEta=ePart.eEta;
int eXi=ePart.eXi;
double eSpatialCoeff=ePart.eCoeff;
SingleSearchEntry eEnt={iEntry, eEta, eXi, alphaLow, eSpatialCoeff};
SingleSearchEntry fEnt={iEntry, eEta, eXi, alphaUpp, eSpatialCoeff};
ListListCases[iTimeLow].push_back(eEnt);
ListListCases[iTimeUpp].push_back(fEnt);
}
}
else {
eValAssign=-7e2;
}
ListEntry[iEntry].Swh_model=eValAssign;
ListEntry[iEntry].WindSpeed_model=eValAssign;
}
for (int iTime=0; iTime<nbTime; iTime++) {
int nbCase=ListListCases[iTime].size();
double eTimeDay=eArr.ListTime[iTime];
if (nbCase > 0) {
if (eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_WNDMAG")) {
RecVar eRecVar=ModelSpecificVarSpecificTime(TotalArr, "WINDMAG", eTimeDay);
for (auto& eSingEntry : ListListCases[iTime]) {
int iEntry=eSingEntry.iEntry;
int eEta=eSingEntry.eEta;
int eXi=eSingEntry.eXi;
double TimeCoeff=eSingEntry.TimeCoeff;
double SpatialCoeff=eSingEntry.SpatialCoeff;
ListEntry[iEntry].WindSpeed_model=ListEntry[iEntry].WindSpeed_model + TimeCoeff*SpatialCoeff*eRecVar.F(eEta,eXi);
}
}
if (eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_HS")) {
RecVar eRecVar=ModelSpecificVarSpecificTime(TotalArr, "Hwave", eTimeDay);
for (auto& eSingEntry : ListListCases[iTime]) {
int iEntry=eSingEntry.iEntry;
int eEta=eSingEntry.eEta;
int eXi=eSingEntry.eXi;
double TimeCoeff=eSingEntry.TimeCoeff;
double SpatialCoeff=eSingEntry.SpatialCoeff;
ListEntry[iEntry].Swh_model=ListEntry[iEntry].Swh_model + TimeCoeff*SpatialCoeff*eRecVar.F(eEta,eXi);
}
}
}
}
}
std::vector<SingleEntryAltimeter> RETRIEVE_RELEVANT_ALTI_DATA(GridArray const& GrdArr, FullNamelist const& eFull)
{
std::map<std::string, SingleBlock> ListBlock=eFull.ListBlock;
SingleBlock eBlSELECT=eFull.ListBlock.at("SELECT");
double BeginTime=CT2MJD(eBlSELECT.ListStringValues.at("BEGTC"));
double EndTime =CT2MJD(eBlSELECT.ListStringValues.at("ENDTC"));
double LowestDay=floor(BeginTime);
double HighestDay=ceil(EndTime)-1;
std::string eStrEnvVar="ALTIMETER_DIRECTORY";
std::string ePrefixAlti=STRING_GETENV(eStrEnvVar);
std::cerr << "ALTIMETRY_DIRECTORY=" << ePrefixAlti << "\n";
std::vector<SingleEntryAltimeter> eRetList;
std::vector<int> ListStatusSatellite=GetListSatelliteId_vect(eFull);
for (int eDay=int(LowestDay); eDay<=int(HighestDay); eDay++) {
double eDayDoubl=double(eDay);
std::vector<int> eDate=DATE_ConvertMjd2six(eDayDoubl);
int year=eDate[0];
int month=eDate[1];
int day=eDate[2];
std::string eFileAlti=ePrefixAlti + StringNumber(year,4) + "/" + StringNumber(month,2) + "/wm_" + StringNumber(year,4) + StringNumber(month,2) + StringNumber(day,2) + ".nc";
if (IsExistingFile(eFileAlti) == false) {
std::cerr << "Missing file " + eFileAlti + "\n";
std::cerr << "Environment variable ALTIMETER_DIRECTORY=" << ePrefixAlti << "\n";
std::cerr << "Please use perl script DownloadAltimeterIfremer\n";
std::cerr << "for downloading the data\n";
throw TerminalException{1};
}
std::vector<SingleEntryAltimeter> eV=READ_ALTI_FILE(eFileAlti, ListStatusSatellite, eFull);
eRetList.insert(eRetList.end(), eV.begin(), eV.end());
}
InterpolateAltimeterData(eRetList, GrdArr, eFull);
return eRetList;
}
void RAW_STATISTICS_ALTIMETER(std::vector<PairListWindWave> const & eSS,
std::vector<SatelliteListTrack> const& LTrack,
PermanentInfoDrawing const& ePerm)
{
std::cerr << "Running RAW_STATISTICS_ALTIMETER\n";
std::string FileName=ePerm.eDir + "Statistics.txt";
std::string FileNameShort=ePerm.eDir + "StatisticsShort.txt";
std::ofstream os;
os.open(FileName);
std::ofstream osShort;
osShort.open(FileNameShort);
os << "Comparison of model results with altimeter\n";
SingleBlock eBlSELECT=ePerm.eFull.ListBlock.at("SELECT");
double BeginTime=CT2MJD(eBlSELECT.ListStringValues.at("BEGTC"));
double EndTime =CT2MJD(eBlSELECT.ListStringValues.at("ENDTC"));
std::string strBegin=DATE_ConvertMjd2mystringPres(BeginTime);
std::string strEnd=DATE_ConvertMjd2mystringPres(EndTime);
os << "Beginning time = " << strBegin << "\n";
os << "Ending    time = " << strEnd << "\n";
osShort << "Beginning time = " << strBegin << "\n";
osShort << "Ending    time = " << strEnd << "\n";
SingleBlock eBlPROC=ePerm.eFull.ListBlock.at("PROC");
os << "HisPrefix = " << eBlPROC.ListStringValues.at("HisPrefix") << "\n";
for (auto & eRec : LTrack) {
std::string SatName=GetNameOfSatelliteAltimeter(eRec.eSat);
os << "Satellite = " << SatName << " TrackDistKM=" << eRec.avgDistKM << "\n";
osShort << "Satellite = " << SatName << " TrackDistKM=" << eRec.avgDistKM << "\n";
}
os << "Do wind statistic = " << ePerm.eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_WNDMAG") << "\n";
os << "Do wave statistic = " << ePerm.eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_HS") << "\n";
if (ePerm.eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_WNDMAG")) {
std::vector<PairMM> TotalListPairWind;
for (auto & eRec : eSS) {
TotalListPairWind.insert(TotalListPairWind.end(), eRec.ListPairWind.begin(), eRec.ListPairWind.end());
}
std::string eNameWind="All satellites wind statistics";
T_stat eStatWind=ComputeStatistics_Pair(TotalListPairWind);
Print_Down_Statistics(os, eNameWind, eStatWind);
for (auto & eRec : eSS) {
std::string SatName=GetNameOfSatelliteAltimeter(eRec.eSat);
std::string eNameWindB=SatName + " wind statistics";
T_stat eStatWindB=ComputeStatistics_Pair(eRec.ListPairWind);
Print_Down_Statistics(os, eNameWindB, eStatWindB);
osShort << eNameWindB << "\n";
osShort << eStatWindB.strNature << "\n";
osShort << eStatWindB.str << "\n";
}
}
if (ePerm.eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_HS")) {
std::vector<PairMM> TotalListPairWave;
for (auto & eRec : eSS) {
TotalListPairWave.insert(TotalListPairWave.end(), eRec.ListPairWave.begin(), eRec.ListPairWave.end());
}
std::string eNameWave="All satellites wave statistics";
T_stat eStatWave=ComputeStatistics_Pair(TotalListPairWave);
Print_Down_Statistics(os, eNameWave, eStatWave);
for (auto & eRec : eSS) {
std::string SatName=GetNameOfSatelliteAltimeter(eRec.eSat);
std::string eNameWaveB=SatName + " wave statistics";
T_stat eStatWaveB=ComputeStatistics_Pair(eRec.ListPairWave);
Print_Down_Statistics(os, eNameWaveB, eStatWaveB);
osShort << eNameWaveB << "\n";
osShort << eStatWaveB.strNature << "\n";
osShort << eStatWaveB.str << "\n";
}
}
os.close();
osShort.close();
}
void RAW_SCATTER_ALTIMETER(std::ostream & os, std::vector<PairListWindWave> const & eSS,
NCLcaller & eCall,
PermanentInfoDrawing const& ePerm)
{
std::cerr << "Running RAW_SCATTER_ALTIMETER\n";
SingleBlock eBlSELECT=ePerm.eFull.ListBlock.at("SELECT");
double BeginTime=CT2MJD(eBlSELECT.ListStringValues.at("BEGTC"));
double EndTime =CT2MJD(eBlSELECT.ListStringValues.at("ENDTC"));
std::string strBegin=DATE_ConvertMjd2mystringPres(BeginTime);
std::string strEnd=DATE_ConvertMjd2mystringPres(EndTime);
int idWind=1;
int idWave=2;
auto fPlot=[&](std::vector<PairMM> const& LPair, int const& idWindWave, int const& idSat) -> void {
DrawScatterArr eDrw;
int siz=LPair.size();
MyVector<double> eVectA(siz);
MyVector<double> eVectB(siz);
for (int i=0; i<siz; i++) {
PairMM ePair=LPair[i];
eVectA[i]=ePair.Meas;
eVectB[i]=ePair.Model;
}
std::vector<double> data_rangeA(2);
std::vector<double> data_rangeB(2);
std::string SatName, SatNameFile;
if (idSat == -1) {
SatName="all satellites";
SatNameFile="allsat";
}
else {
SatName=GetNameOfSatelliteAltimeter(idSat);
SatNameFile=SatName;
}
std::string VarType, eUnit;
double TheMax;
if (idWindWave == idWind) {
VarType="Wind";
eUnit="m/s";
TheMax=20;
}
else {
VarType="Wave";
eUnit="m";
TheMax=4;
}
data_rangeA[0]=0;
data_rangeA[1]=TheMax;
data_rangeB[0]=0;
data_rangeB[1]=TheMax;
eDrw.VarNameAB_file="Scatter_" + SatNameFile + "_" + VarType;
eDrw.DoTitle=false;
eDrw.AddStatMeasModel=true;
eDrw.NameA_plot="Data (" + eUnit + ")";
eDrw.NameB_plot="Model (" + eUnit + ")";
eDrw.data_rangeA=data_rangeA;
eDrw.data_rangeB=data_rangeB;
eDrw.eVectA=eVectA;
eDrw.eVectB=eVectB;
eDrw.aSize=100;
eDrw.bSize=100;
PLOT_SCATTER(eDrw, eCall, ePerm);
};
if (ePerm.eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_WNDMAG")) {
std::vector<PairMM> TotalListPairWind;
for (auto & eRec : eSS)
TotalListPairWind.insert(TotalListPairWind.end(), eRec.ListPairWind.begin(), eRec.ListPairWind.end());
fPlot(TotalListPairWind, idWind, -1);
for (auto & eRec : eSS)
fPlot(eRec.ListPairWind, idWind, eRec.eSat);
}
if (ePerm.eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_HS")) {
std::vector<PairMM> TotalListPairWave;
for (auto & eRec : eSS)
TotalListPairWave.insert(TotalListPairWave.end(), eRec.ListPairWave.begin(), eRec.ListPairWave.end());
fPlot(TotalListPairWave, idWave, -1);
for (auto & eRec : eSS)
fPlot(eRec.ListPairWave, idWave, eRec.eSat);
}
}
void RAW_PLOT_ALTIMETER_TRACKS(std::ostream & os, std::vector<SatelliteListTrack> const& LTrack,
NCLcaller & eCall,
PermanentInfoDrawing const& ePerm)
{
std::cerr << "Running RAW_PLOT_ALTIMETER_TRACKS\n";
double MinLon, MinLat, MaxLon, MaxLat;
SingleBlock eBlSEL=ePerm.eFull.ListBlock.at("SELECT");
if (eBlSEL.ListIntValues.at("GEOSELECTION") == 1) {
MinLon=eBlSEL.ListDoubleValues.at("MinLON");
MaxLon=eBlSEL.ListDoubleValues.at("MaxLON");
MinLat=eBlSEL.ListDoubleValues.at("MinLAT");
MaxLat=eBlSEL.ListDoubleValues.at("MaxLAT");
}
else {
std::vector<double> ListLon=eBlSEL.ListListDoubleValues.at("LONPOLY");
std::vector<double> ListLat=eBlSEL.ListListDoubleValues.at("LATPOLY");
int siz=ListLon.size();
MinLon=ListLon[0];
MaxLon=ListLon[0];
MinLat=ListLat[0];
MaxLat=ListLat[0];
for (int i=1; i<siz; i++) {
double eLon=ListLon[i];
double eLat=ListLat[i];
if (MinLon > eLon)
MinLon=eLon;
if (MaxLon < eLon)
MaxLon=eLon;
if (MinLat > eLat)
MinLat=eLat;
if (MaxLat < eLat)
MaxLat=eLat;
}
}
QuadArray eQuad{MinLon, MaxLon, MinLat, MaxLat};
int nbSplitLon=100;
int nbSplitLat=100;
GridArray GrdArr=TRIVIAL_GRID_ARRAY(eQuad, nbSplitLon, nbSplitLat);
RecVar eRecVar=GetTrivialArrayPlot(GrdArr);
auto fPlot=[&](std::vector<SeqLineSegment> const& TheList, int const& idSat) -> void {
std::string SatName, SatNameFile;
if (idSat == -1) {
SatName="all satellites";
SatNameFile="allsat";
}
else {
SatName=GetNameOfSatelliteAltimeter(idSat);
SatNameFile=SatName;
}
std::string TitleStr="Tracks for " + SatName;
std::string FileName=ePerm.eDir + "Tracks_for_" + SatNameFile;
DrawArr eDrw;
eDrw.DoColorBar=false;
eDrw.eQuadFrame=eQuad;
eDrw.DoTitle=true;
eDrw.TitleStr=TitleStr;
eDrw.ColorMap="WhBlGrYeRe";
eDrw.DoTitle=false;
eDrw.ListLineSegment=TheList;
PLOT_PCOLOR(FileName, GrdArr, eDrw, eRecVar, eCall, ePerm);
};
struct SatelliteListTrack {
int eSat;
std::vector<SeqLineSegment> ListLineSegment;
};
std::vector<SatelliteListTrack> ListListLineSegment;
for (auto& eRecTrack : LTrack) {
std::vector<SeqLineSegment> eListLineSegment;
for (auto& eListEnt : eRecTrack.ListListEntAltimeter) {
std::vector<PairLL> ListPairLL;
for (auto & eEnt : eListEnt) {
double eLon=eEnt.Lon;
double eLat=eEnt.Lat;
PairLL ePairLL{eLon, eLat};
ListPairLL.push_back(ePairLL);
}
eListLineSegment.push_back({ListPairLL, false});
}
ListListLineSegment.push_back({eRecTrack.eSat, eListLineSegment});
}
std::vector<SeqLineSegment> ListTotal;
for (auto& eRec : ListListLineSegment)
ListTotal.insert(ListTotal.end(), eRec.ListLineSegment.begin(), eRec.ListLineSegment.end());
fPlot(ListTotal, -1);
for (auto& eRec : ListListLineSegment)
fPlot(eRec.ListLineSegment, eRec.eSat);
}
void RAW_PLOT_VALUE_TRACKS(std::ostream & os, std::vector<SatelliteListTrack> const& LTrack,
NCLcaller & eCall,
PermanentInfoDrawing const& ePerm)
{
std::cerr << "Running RAW_PLOT_VALUE_TRACKS\n";
int idWind=1;
int idWave=2;
int MinEntryTrackPlot=ePerm.eFull.ListBlock.at("PROCESS").ListIntValues.at("MinEntryTrackPlot");
auto fPlot=[&](MyVector<double> const& ListLat, MyVector<double> const& ListMeas, MyVector<double> const& ListModel, int const& idWindWave, int const& eSat, int const& iTrack, double const& eTimeDay) -> void {
std::string eVarName;
double TheMin, TheMax;
if (idWindWave == idWind) {
eVarName="Wind";
TheMin=0;
TheMax=20;
}
else {
eVarName="Wave";
TheMin=0;
TheMax=4;
}
std::string SatName=GetNameOfSatelliteAltimeter(eSat);
std::string SatNameFile=SatName;
std::string strPres=DATE_ConvertMjd2mystringPres(eTimeDay);
std::string strFile=DATE_ConvertMjd2mystringFile(eTimeDay);
DrawLinesArr eDrawArr;
eDrawArr.DoTitle=true;
eDrawArr.TitleStr="Track Nr" + IntToString(iTrack) + " of " + eVarName + " for " + SatName + " at " + strPres;
std::string fVarName=SatName + "_" + eVarName + "_Track" + StringNumber(iTrack, 4);
eDrawArr.IsTimeSeries=false;
eDrawArr.PairComparison=true;
eDrawArr.VarName=fVarName;
eDrawArr.TheMax=TheMax;
eDrawArr.TheMin=TheMin;
eDrawArr.ListX=ListLat;
eDrawArr.ListListVect={ListMeas, ListModel};
eDrawArr.ListName_plot={"meas.", "model"};
std::string FileName=ePerm.eDir + SatNameFile + "_" + eVarName + "_Track" + StringNumber(iTrack,4) + "_at_" + strFile;
LINES_PLOT(FileName, eDrawArr, eCall, ePerm);
};
for (auto& eRec : LTrack) {
int eSat=eRec.eSat;
int iTrack=0;
for (auto & SingBlock : eRec.ListListEntAltimeter) {
int TheSize=SingBlock.size();
if (TheSize > MinEntryTrackPlot) {
MyVector<double> ListLat(TheSize);
MyVector<double> ListMeasWind(TheSize);
MyVector<double> ListModelWind(TheSize);
MyVector<double> ListMeasWave(TheSize);
MyVector<double> ListModelWave(TheSize);
int idx=0;
double SumTimeDay=0;
for (auto& eEnt : SingBlock) {
ListLat[idx]=eEnt.Lat;
ListMeasWind[idx]=eEnt.WindSpeed_used;
ListModelWind[idx]=eEnt.WindSpeed_model;
ListMeasWave[idx]=eEnt.Swh_used;
ListModelWave[idx]=eEnt.Swh_model;
SumTimeDay += eEnt.Time;
idx++;
}
double eTimeDay=SumTimeDay/double(idx);
if (ePerm.eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_WNDMAG"))
fPlot(ListLat, ListMeasWind, ListModelWind, idWind, eSat, iTrack, eTimeDay);
if (ePerm.eFull.ListBlock.at("PROCESS").ListBoolValues.at("DO_HS"))
fPlot(ListLat, ListMeasWave, ListModelWave, idWave, eSat, iTrack, eTimeDay);
iTrack++;
}
}
}
}
void SAVE_TO_TEXT_FILE(std::vector<SatelliteListTrack> const& LTrack,
PermanentInfoDrawing const& ePerm)
{
std::cerr << "Running SAVE_TO_TEXT_FILE\n";
std::string eName=ePerm.eFull.ListBlock.at("PROCESS").ListStringValues.at("FILE_SAVE_TXT");
std::string eFileName=ePerm.eDir + eName;
std::ofstream os;
os.open(eFileName);
double LargeValue=30000;
for (auto& eRec : LTrack) {
int eSat=eRec.eSat;
int iTrack=0;
std::string SatName=GetNameOfSatelliteAltimeter(eSat);
int nbTrack=eRec.ListListEntAltimeter.size();
os << "Satellite=" << SatName << " nbTrack=" << nbTrack << "\n";
for (auto & SingBlock : eRec.ListListEntAltimeter) {
std::vector<double> ListMeasWind;
std::vector<double> ListModelWind;
std::vector<double> ListMeasWave;
std::vector<double> ListModelWave;
std::string strMeasWave, strModelWave, strMeasWind, strModelWind;
int idx=0;
double SumTimeDay=0;
double MinTime=0, MaxTime=0;
std::vector<std::string> ListStrMeasWind;
std::vector<std::string> ListStrModelWind;
std::vector<std::string> ListStrMeasWave;
std::vector<std::string> ListStrModelWave;
for (auto& eEnt : SingBlock) {
if (eEnt.WindSpeed_used < LargeValue && eEnt.WindSpeed_model < LargeValue) {
ListMeasWind.push_back(eEnt.WindSpeed_used);
ListModelWind.push_back(eEnt.WindSpeed_model);
strMeasWind=DoubleTo4dot2f(eEnt.WindSpeed_used);
strModelWind=DoubleTo4dot2f(eEnt.WindSpeed_model);
}
else {
strMeasWind="NA";
strModelWind="NA";
}
ListStrMeasWind.push_back(strMeasWind);
ListStrModelWind.push_back(strModelWind);
if (eEnt.Swh_used < LargeValue && eEnt.Swh_model < LargeValue) {
ListMeasWave.push_back(eEnt.Swh_used);
ListModelWave.push_back(eEnt.Swh_model);
strMeasWave=DoubleTo4dot2f(eEnt.Swh_used);
strModelWave=DoubleTo4dot2f(eEnt.Swh_model);
}
else {
strMeasWave="NA";
strModelWave="NA";
}
ListStrMeasWave.push_back(strMeasWave);
ListStrModelWave.push_back(strModelWave);
SumTimeDay += eEnt.Time;
if (idx == 0) {
MinTime=eEnt.Time;
MaxTime=eEnt.Time;
}
else {
if (eEnt.Time > MaxTime)
MaxTime=eEnt.Time;
if (eEnt.Time < MinTime)
MinTime=eEnt.Time;
}
idx++;
}
T_stat eStatWave=ComputeStatistics_vector(ListMeasWave, ListModelWave);
T_stat eStatWind=ComputeStatistics_vector(ListMeasWind, ListModelWind);
double eTimeMid=SumTimeDay/double(idx);
std::string strPresMid=DATE_ConvertMjd2mystringPres(eTimeMid);
std::string strPresMin=DATE_ConvertMjd2mystringPres(MinTime);
std::string strPresMax=DATE_ConvertMjd2mystringPres(MaxTime);
os << "--------------------------------------------------------------\n";
os << "Satellite=" << SatName << "\n";
os << "iTrack=" << iTrack << "/ " << nbTrack << " TimeMid=" << strPresMid << " len=" << idx << "\n";
os << "Time min=" << strPresMin << " max=" << strPresMax << "\n";
os << "Wind: " << eStatWind.strNature << "\n";
os << "Wind: " << eStatWind.str << "\n";
os << "Wave: " << eStatWave.strNature << "\n";
os << "Wave: " << eStatWave.str << "\n";
os << "All data len=" << idx << "\n";
os << "Wave(Model/Meas) Wind(Model/Meas)\n";
for (int i=0; i<idx; i++) {
os << ListStrModelWave[i] << " / " << ListStrMeasWave[i] << "     " << ListStrModelWind[i] << " / " << ListStrMeasWind[i] << "\n";
}
iTrack++;
}
}
os.close();
}
void Process_Altimetry_Comparison_Request(FullNamelist const& eFull)
{
PermanentInfoDrawing ePerm=GET_PERMANENT_INFO(eFull);
NCLcaller eCall(ePerm.KeepNC_NCL, ePerm.NPROC);
SingleBlock eBlPROC=eFull.ListBlock.at("PROC");
std::string eModelName=eBlPROC.ListStringValues.at("MODELNAME");
std::string GridFile=eBlPROC.ListStringValues.at("GridFile");
std::string HisPrefix=eBlPROC.ListStringValues.at("HisPrefix");
TripleModelDesc eTriple{eModelName, GridFile, "unset", HisPrefix, false, false, 0, 0};
GridArray GrdArr=RETRIEVE_GRID_ARRAY(eTriple);
std::vector<std::vector<int> > ListListStatus;
std::vector<SingleEntryAltimeter> PreListSingleEntry=RETRIEVE_RELEVANT_ALTI_DATA(GrdArr, eFull);
std::set<int> SatelliteId=GetListSatelliteId_set(PreListSingleEntry, eFull);
if (eFull.ListBlock.at("SELECT").ListBoolValues.at("DoMinDistCoast") == true)
ListListStatus.push_back(FilterByMinDistCoast(PreListSingleEntry, eFull));
if (eFull.ListBlock.at("SELECT").ListBoolValues.at("EliminationShortTrack") == true)
ListListStatus.push_back(GetListStatusTrackLength(PreListSingleEntry, eFull));
std::vector<SingleEntryAltimeter> ListSingleEntry=SelectByStatus(ListListStatus, PreListSingleEntry);
double avgDistKM_target=GetGridSpacing(GrdArr);
std::vector<SatelliteListTrack> ListTrackInfo=GetListTrackAltimeter(ListSingleEntry, avgDistKM_target, eFull);
std::vector<PairListWindWave> eSS=MergeTracksForRawStatistics(ListTrackInfo, eFull);
if (eBlPROC.ListBoolValues.at("DO_STAT")) {
RAW_STATISTICS_ALTIMETER(eSS, ListTrackInfo, ePerm);
}
if (eBlPROC.ListBoolValues.at("DO_SCATTERPLOT")) {
RAW_SCATTER_ALTIMETER(std::cout, eSS, eCall, ePerm);
}
if (eBlPROC.ListBoolValues.at("DO_SAVE_TXT")) {
SAVE_TO_TEXT_FILE(ListTrackInfo, ePerm);
}
if (eBlPROC.ListBoolValues.at("PLOT_ALL_TRACKS")) {
RAW_PLOT_ALTIMETER_TRACKS(std::cout, ListTrackInfo, eCall, ePerm);
}
if (eBlPROC.ListBoolValues.at("PLOT_TRACKS")) {
RAW_PLOT_VALUE_TRACKS(std::cout, ListTrackInfo, eCall, ePerm);
}
}
int main(int argc, char *argv[])
{
try {
FullNamelist eFull=NAMELIST_GetStandardALTIMETRY_COMPARISON();
if (argc != 2) {
std::cerr << "AltimeterComparison is used as\n";
std::cerr << "AltimeterComparison [alti.nml]\n";
std::cerr << "with alti.nml the file describing the chosen options\n";
return -1;
}
std::string eFileName=argv[1];
NAMELIST_ReadNamelistFile(eFileName, eFull);
Process_Altimetry_Comparison_Request(eFull);
}
catch (TerminalException const& e) {
exit(e.eVal);
}
}
