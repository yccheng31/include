#pragma once

#ifndef DEF_H
#define DEF_H

#include <stdio.h>
//#include <direct.h>

#include <complex>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <limits>
#include <cassert>
using namespace std;

#ifdef _MSC_VER
#include <io.h>
#include <direct.h>
#include <hash_map>
using namespace stdext;
#define unordered_map hash_map
#elif defined(__GNUC__)
#include <tr1/unordered_map>
using namespace std::tr1;
#endif

#include "connected.h"
#include "matrix.h"
using namespace math;


#include <sys/stat.h>
#include <sys/types.h> 
#include <float.h>
#include <stdlib.h>
#include <stdio.h>


#include "cv.h"
#include "cxcore.h"
#include "cvaux.h"
#include "highgui.h"
#include "ml.h"
using namespace cv;

//#ifdef _MSC_VER
//#include <crtdbg.h>  
////#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
//#ifdef _DEBUG
////#define CRTDBG_MAP_ALLOC
////#define new DEBUG_NEW
////#undef THIS_FILE
////static char THIS_FILE[] = __FILE__;
//#endif
//#endif

#define MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN 15
#define PI 3.14159265358979323846264338327950288
#define MY_MAX_DOUBLE 1.7976931348623158e+308

// From HTK
#define LZERO  (-1.0E10)   /* ~log(0) */
#define LSMALL (-0.5E10)   /* log values < LSMALL are set to LZERO */
#define DIAGONAL
#define LOGHALF -0.69314718055994530941723212145818
#define TWOPI 6.28318530717958647692528676655900576
#define LOGTWOPI 1.8378770664093454835606594728112
#define LOGLOGTWOPI 0.60861113726210809872881467715665
static double minLogExp = -log(-LZERO);

#ifndef DWORD
typedef unsigned long       DWORD;
#endif

#ifndef HRESULT
typedef long       HRESULT;
#define S_OK                             0L
#define S_FALSE                          1L
#define ERROR_OUT_OF_PAPER               28L
#define ERROR_OPEN_FAILED                110L
#define ERROR_INVALID_ADDRESS            487L
#define ERROR_INVALID_INDEX              1413L
#define ERROR_INVALID_DATA		 0x8007000DL
#define ERROR_EMPTY                      4306L
#define E_POINTER                        0x80004003L
#define E_NOTIMPL                        0x80004001L // somehow visual studio 2010 cannot find it, so we define it, 2014/04/16
#endif

#ifdef _DEBUG
   #define DEBUG_CLIENTBLOCK   new( _CLIENT_BLOCK, __FILE__, __LINE__)
#else
   #define DEBUG_CLIENTBLOCK
#endif


#ifndef isnan
inline bool Myisnan(double x) {
	return x != x;
}
#endif

#ifndef isinf
inline bool Myisinf(double x) {	
	return numeric_limits<double>::has_infinity && x == numeric_limits<double>::infinity();
}
#endif

#ifndef _MAX_FNAME
#define _MAX_FNAME 4096
#endif

inline bool file_exist (const std::string& name) {
  struct stat buf;
  char lastChar = name.at( name.length() - 1 );
  string name_copy = name;
  if(lastChar == '\\' || lastChar == '/')
  {
	  name_copy = name.substr(0, name.length()-1);
  }
  int nresult = stat (name_copy.c_str(), &buf);
  return (nresult == 0); 
}

#ifndef FENode
struct FENode
{
	vector<IplImage*> vpImgs;
	string szOutFeatureFile;
	string szOutFeatureFile_uhalf;
	string szOutFeatureFile_lhalf;
	string szErrLog;
	string szMLFName;
	string szLabelName;
	string szKeyName;
	bool IsCreated;
};
#endif

#ifndef DoublePoint
struct DoublePoint
{
	double x;
	double y;
}typedef DoublePoint;
#endif
#ifndef VideoNode
struct VideoNode
{
	int nStart;
	int nStop;
	string szName;
}typedef VideoNode;
#endif

// For sorting
struct LessDereference {
	template <class T>
	bool operator()(const T * lhs, const T * rhs) const {
	return *lhs > *rhs;
	}
};
template<class T> struct index_cmp {
	index_cmp(const T arr) : arr(arr) {}
	bool operator()(const size_t a, const size_t b) const
	{ return arr[a] > arr[b]; }
	const T arr;
};
template<class T> struct index_cmp_ascend {
	index_cmp_ascend(const T arr) : arr(arr) {}
	bool operator()(const size_t a, const size_t b) const
	{ return arr[a] < arr[b]; }
	const T arr;
};


// For HMM stuff
class Mixture
{
public:
	Mixture()
	{
		m_dbWeight = 0.0;
		m_dbgconst = 0.0;
	}
	double m_dbWeight;
	vector<double> m_vMean;
	vector<double> m_vDiagVar;
	double m_dbgconst;
};

class State
{
public:
	State()
	{
		m_nNumMixs = 0;
	}
	int m_nNumMixs;
	vector<Mixture> m_Mixtures;
};

class Model
{
public:
	Model()
	{
		m_nNumStates = 0;
	}
	string m_szName;
	int m_nNumStates;
	vector<double> m_vPis;
	vector< vector<double> > m_vvTransp;
	vector<State> m_States;
};
class GaussianNode
{

public:

	GaussianNode()
	{
		m_nPrevClosestMix = 0;
	}
	void SetDimMix(int nDim, int nMixs)
	{
		m_vdbWeight.resize(nMixs, 1.0/(double)nMixs);
		m_vvdbMean.resize(nMixs);
		m_vvdbDiagvar.resize(nMixs);
		for(int n = 0; n < nMixs; n++)
		{
			m_vvdbMean[n].resize(nDim, 0.0);
			m_vvdbDiagvar[n].resize(nDim, 0.0);
		}
	}
	void Mixup(int nNexMixs)
	{
		int nCurrMixSize = (int)m_vdbWeight.size();
		m_vdbWeight.resize(nNexMixs);
		m_vvdbMean.resize(nNexMixs);
		m_vvdbDiagvar.resize(nNexMixs);
		for(int n = nCurrMixSize; n < nNexMixs; n++)
		{
			m_vvdbMean[n].resize(m_vvdbMean[0].size(), 0.0);
			m_vvdbDiagvar[n].resize(m_vvdbMean[0].size(), 0.0);
		}
	}
	vector<double> m_vdbWeight;
	vector< vector<double> > m_vvdbMean;
	vector< vector<double> > m_vvdbDiagvar;
	int m_nPrevClosestMix;
};

#endif
