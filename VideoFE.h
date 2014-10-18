#pragma once

#include "LoadJPGLabelInfo.h"
#include "GMM.h"
#include "GMMAdaptTool.h"
#include "LocalFeatuers.h"

#ifdef _MSC_VER
#include "windows.h"  // needed for Visual studio 2010, 2014/04/16
#include "Mmsystem.h"
#pragma comment( lib, "Winmm.lib")
#endif

//#include <stdio.h>
//#include <direct.h>
//#include <io.h>
//#include "matrix.h"
//using namespace math;
//#include <cctype>
//#include <vector>
//#include <list>
//#include <fstream>
//#include <iostream>
//#include <hash_map>
//
//#include "cv.h"
//#include "cxcore.h"
//#include "cvaux.h"
//#include "highgui.h"
//#include "ml.h"
//
////#include "include/cv.h"
////#include "include/cxcore.h"
////#include "include/cvaux.h"
////#include "include/highgui.h"
////#include "include/ml.h"
//
//using namespace std;
//using namespace stdext;
//
//#define MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN 15
//#define PI 3.14159265358979323846264338327950288
//#define MY_MAX_DOUBLE 1.7976931348623158e+308
//
//// From HTK
//#define LZERO  (-1.0E10)   /* ~log(0) */
//#define LSMALL (-0.5E10)   /* log values < LSMALL are set to LZERO */
//#define DIAGONAL
//#define LOGHALF -0.69314718055994530941723212145818
//#define TWOPI 6.28318530717958647692528676655900576
//#define LOGTWOPI 1.8378770664093454835606594728112
//#define LOGLOGTWOPI 0.60861113726210809872881467715665
//static double minLogExp = -log(-LZERO);
//
//
//#ifndef FENode
//struct FENode
//{
//	vector<IplImage*> vpImgs;
//	string szOutFeatureFile;
//	string szOutFeatureFile_uhalf;
//	string szOutFeatureFile_lhalf;
//	string szErrLog;
//	string szMLFName;
//	string szLabelName;
//	string szKeyName;
//};
//#endif
//
//#ifndef DoublePoint
//struct DoublePoint
//{
//	double x;
//	double y;
//}typedef DoublePoint;
//#endif
//#ifndef VideoNode
//struct VideoNode
//{
//	int nStart;
//	int nStop;
//	string szName;
//}typedef VideoNode;
//#endif
//
//// For sorting
//struct LessDereference {
//	template <class T>
//	bool operator()(const T * lhs, const T * rhs) const {
//	return *lhs > *rhs;
//	}
//};
//
//// For HMM stuff
//class Mixture
//{
//public:
//	Mixture()
//	{
//		m_dbWeight = 0.0;
//		m_dbgconst = 0.0;
//	}
//	double m_dbWeight;
//	vector<double> m_vMean;
//	vector<double> m_vDiagVar;
//	double m_dbgconst;
//};
//
//class State
//{
//public:
//	State()
//	{
//		m_nNumMixs = 0;
//	}
//	int m_nNumMixs;
//	vector<Mixture> m_Mixtures;
//};
//
//class Model
//{
//public:
//	Model()
//	{
//		m_nNumStates = 0;
//	}
//	string m_szName;
//	int m_nNumStates;
//	vector<double> m_vPis;
//	vector< vector<double> > m_vvTransp;
//	vector<State> m_States;
//};
//class GaussianNode
//{
//
//public:
//
//	GaussianNode()
//	{
//		m_nPrevClosestMix = 0;
//	}
//	void SetDimMix(int nDim, int nMixs)
//	{
//		m_vdbWeight.resize(nMixs, 1.0/(double)nMixs);
//		m_vvdbMean.resize(nMixs);
//		m_vvdbDiagvar.resize(nMixs);
//		for(int n = 0; n < nMixs; n++)
//		{
//			m_vvdbMean[n].resize(nDim, 0.0);
//			m_vvdbDiagvar[n].resize(nDim, 0.0);
//		}
//	}
//	void Mixup(int nNexMixs)
//	{
//		int nCurrMixSize = (int)m_vdbWeight.size();
//		m_vdbWeight.resize(nNexMixs);
//		m_vvdbMean.resize(nNexMixs);
//		m_vvdbDiagvar.resize(nNexMixs);
//		for(int n = nCurrMixSize; n < nNexMixs; n++)
//		{
//			m_vvdbMean[n].resize(m_vvdbMean[0].size(), 0.0);
//			m_vvdbDiagvar[n].resize(m_vvdbMean[0].size(), 0.0);
//		}
//	}
//	vector<double> m_vdbWeight;
//	vector< vector<double> > m_vvdbMean;
//	vector< vector<double> > m_vvdbDiagvar;
//	int m_nPrevClosestMix;
//};
//
//// Skin color GMM class
//class GMM
//{
//public:
//	GMM()
//	{
//		m_nVecSize = 0;
//		m_nMixes = 0;
//		m_szCovType = "diagonal";
//	}
//	GMM(int nVecSize, int nMixes = 1, string szCovType = "diagonal")
//	{
//		m_nVecSize = nVecSize;
//		m_nMixes = nMixes;
//		m_vmdMean.resize(m_nMixes);
//		m_szCovType = szCovType;
//		m_vMixWeights = vector<double>(nMixes, (1.0/(double)nMixes));
//		for(int n = 0; n < m_nMixes; n++)
//		{
//			m_vmdMean[n].SetSize(m_nVecSize, 1);
//			m_vmdMean[n].Null();
//			m_vmdCovar[n].SetSize(m_nVecSize, m_nVecSize);
//			m_vmdInvCovar[n].SetSize(m_nVecSize, m_nVecSize);
//			m_vmdCovar[n].Null();
//			m_vmdInvCovar[n].Null();
//			for(int d = 0; d < m_nVecSize; d++)
//			{
//				m_vmdCovar[n](d,d) = 1.0;
//				m_vmdInvCovar[n](d,d) = 1.0;
//			}
//		}
//	}
//	GMM(const GMM& rhs)
//	{
//		m_nVecSize = rhs.m_nVecSize;
//		m_nMixes = rhs.m_nMixes;
//		m_szCovType = rhs.m_szCovType;
//		m_szFeatureType = rhs.m_szFeatureType;
//		m_vdgconst = rhs.m_vdgconst;
//		m_vMixWeights = rhs.m_vMixWeights;
//		m_vmdMean = rhs.m_vmdMean;
//		m_vmdCovar = rhs.m_vmdCovar;
//		m_vmdInvCovar = rhs.m_vmdInvCovar;
//		m_vlogMixWeights = rhs.m_vlogMixWeights;
//		m_vvddiaggconst = rhs.m_vvddiaggconst;
//
//	}
//	double LogAdd(double x, double y)
//	{
//	   double temp,diff,z;
//	   if (x<y) {
//		  temp = x; x = y; y = temp;
//	   }
//	   diff = y-x;
//	   if (diff<minLogExp) 
//		  return  (x<LSMALL)?LZERO:x;
//	   else {
//		  z = exp(diff);
//		  return x+log(1.0+z);
//	   }
//	}
//	void InitGMMMeans(int nVecSize, int nMixes = 1)
//	{
//		m_nVecSize = nVecSize;
//		m_nMixes = nMixes;
//		m_vmdMean.resize(m_nMixes);
//		double oneovern = 1.0/(double)nMixes;
//		double log1n = log(oneovern);
//		m_vMixWeights = vector<double>(nMixes, oneovern);
//		m_vlogMixWeights = vector<double>(nMixes, log1n);
//		for(int n = 0; n < m_nMixes; n++)
//		{
//			m_vmdMean[n].SetSize(m_nVecSize, 1);
//			m_vmdMean[n].Null();
//		}
//	}
//	void InitGMMVars(int nVecSize, int nMixIndex, string szCovType = "diagonal")
//	{
//		m_nVecSize = nVecSize;
//		m_szCovType = szCovType;
//		if(nMixIndex >= m_vmdCovar.size())
//		{
//			m_vmdCovar.resize(nMixIndex + 1);
//			m_vmdInvCovar.resize(nMixIndex + 1);
//		}
//		for(int n = 0; n < m_nMixes; n++)
//		{
//			m_vmdCovar[nMixIndex].SetSize(m_nVecSize, m_nVecSize);
//			m_vmdInvCovar[nMixIndex].SetSize(m_nVecSize, m_nVecSize);
//			m_vmdCovar[nMixIndex].Null();
//			m_vmdInvCovar[nMixIndex].Null();
//			for(int d = 0; d < m_nVecSize; d++)
//			{
//				m_vmdCovar[nMixIndex](d,d) = 1.0;
//				m_vmdInvCovar[nMixIndex](d,d) = 1.0;
//			}
//		}
//	}
//	void CompGconsts()
//	{
//		m_vdgconst.resize(m_nMixes, 0.0);
//		m_vvddiaggconst.resize(m_nMixes, vector<double>(m_nVecSize, 0.0));
//		for(int n = 0; n < m_nMixes; n++)
//		{
//			m_vdgconst[n] = (double)m_nVecSize * log(TWOPI) + log(m_vmdCovar[n].Det());
//			for(int d = 0; d < m_nVecSize; d++)
//			{
//				m_vvddiaggconst[n][d] = log(TWOPI) + log(m_vmdCovar[n](d,d));
//			}
//		}
//	}
//	// Actually, return log likelihood instead of likelihood
//	double CompLikelihood(const vector<double>& vec)
//	{
//		matrix<double> mvec;
//		size_t nsize = vec.size();
//		mvec.SetSize(nsize, 1);
//		for(size_t n = 0; n < nsize; n++)
//		{
//			mvec(n, 0) = vec[n];
//		}
//		vector<double> vdTemp(m_nMixes, 0.0);
//		double dTemp = LZERO;
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			matrix<double> x = (mvec - m_vmdMean[n]);
//			if(m_szCovType == "diagonal")
//			{
//#ifndef DIAGONAL			
//				for(size_t j = 0; j < m_nVecSize; j++)
//				{
//					for(size_t i = j+1; i < m_nVecSize; i++)
//					{
//						vdTemp[n] += x(i, 0) * x(j, 0) * m_vmdInvCovar[n](i,j);
//					}
//				}
//				vdTemp[n] *= 2;
//				for(size_t i = 0; i < m_nVecSize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
//				}
//				vdTemp[n] += m_vdgconst[n];
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//#else
//				vdTemp[n] += m_vdgconst[n];
//				for(size_t i = 0; i < m_nVecSize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
//				}
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//#endif
//			}
//			else if(m_szCovType == "full")
//			{
//				vdTemp[n] += m_vdgconst[n];
//				vdTemp[n] += (~x * m_vmdInvCovar[n] * x)(0, 0);
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//			}
//			dTemp = LogAdd(dTemp, vdTemp[n]);
//		}
//		return dTemp; //exp(dTemp)
//	}
//	double CompLikelihood(const matrix<double>& mvec)
//	{
//		vector<double> vdTemp(m_nMixes, 0.0);
//		size_t nsize = mvec.ColNo();
//		double dTemp = LZERO;
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			matrix<double> x = (mvec - m_vmdMean[n]);
//			if(m_szCovType == "diagonal")
//			{
//#ifndef DIAGONAL			
//				for(size_t j = 0; j < m_nVecSize; j++)
//				{
//					for(size_t i = j+1; i < m_nVecSize; i++)
//					{
//						vdTemp[n] += x(i, 0) * x(j, 0) * m_vmdInvCovar[n](i,j);
//					}
//				}
//				vdTemp[n] *= 2;
//				for(size_t i = 0; i < m_nVecSize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
//				}
//				vdTemp[n] += m_vdgconst[n];
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//#else
//				vdTemp[n] += m_vdgconst[n];
//				for(size_t i = 0; i < m_nVecSize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
//				}
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//#endif
//			}
//			else if(m_szCovType == "full")
//			{
//				vdTemp[n] += m_vdgconst[n];
//				vdTemp[n] += (~x * m_vmdInvCovar[n] * x)(0, 0);
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//			}
//			dTemp = LogAdd(dTemp, vdTemp[n]);
//		}
//		return dTemp; //exp(dTemp)
//	}
//	// the posterior probability ( not log )
//	void CompPosteriori(const vector<double>& vec, vector<double>& vPost)
//	{
//		matrix<double> mvec;
//		size_t nsize = vec.size();
//		mvec.SetSize(nsize, 1);
//		for(size_t n = 0; n < nsize; n++)
//		{
//			mvec(n, 0) = vec[n];
//		}
//		vector<double> vdTemp(m_nMixes, 0.0);
//		double dTemp = LZERO;
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			matrix<double> x = (mvec - m_vmdMean[n]);
//			if(m_szCovType == "diagonal")
//			{
//#ifndef DIAGONAL			
//				for(size_t j = 0; j < m_nVecSize; j++)
//				{
//					for(size_t i = j+1; i < m_nVecSize; i++)
//					{
//						vdTemp[n] += x(i, 0) * x(j, 0) * m_vmdInvCovar[n](i,j);
//					}
//				}
//				vdTemp[n] *= 2;
//				for(size_t i = 0; i < m_nVecSize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
//				}
//				vdTemp[n] += m_vdgconst[n];
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//#else
//				vdTemp[n] += m_vdgconst[n];
//				for(size_t i = 0; i < m_nVecSize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
//				}
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//#endif
//			}
//			else if(m_szCovType == "full")
//			{
//				vdTemp[n] += m_vdgconst[n];
//				vdTemp[n] += (~x * m_vmdInvCovar[n] * x)(0, 0);
//				vdTemp[n] *= -0.5;
//				vdTemp[n] += m_vlogMixWeights[n];
//			}
//			dTemp = LogAdd(dTemp, vdTemp[n]);
//		}
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			vdTemp[n] -= dTemp;
//			vdTemp[n] = exp(vdTemp[n]);
//		}
//		vPost = vdTemp;
//	}
//	// Get mean vector for subspace, the space chosen is indicated by vDimConsidered
//	void GetPartialMeans(vector< matrix<double> >& PartialMeans, vector<int>& vDimConsidered)
//	{
//		PartialMeans.resize(m_nMixes);
//		for(size_t k = 0; k < m_vmdMean.size(); k++)
//		{
//			PartialMeans[k].SetSize(vDimConsidered.size(), 1);
//			for(size_t n = 0; n < vDimConsidered.size(); n++)
//			{
//				PartialMeans[k](n, 0) = m_vmdMean[k](vDimConsidered[n], 0);
//			}
//		}
//	}
//	// Get covariance matrix for subspace, the space chosen is indicated by vDimConsidered
//	void GetPartialVars(vector< matrix<double> >& PartialVars, vector< matrix<double> >& PartialInvVars, vector<double>& vPartialGconst, vector<int>& vDimConsidered)
//	{
//		PartialVars.resize(m_nMixes);
//		PartialInvVars.resize(m_nMixes);
//		for(size_t k = 0; k < m_vmdMean.size(); k++)
//		{
//			PartialVars[k].SetSize(vDimConsidered.size(), vDimConsidered.size());
//			PartialInvVars[k].SetSize(vDimConsidered.size(), vDimConsidered.size());
//			for(size_t n = 0; n < vDimConsidered.size(); n++)
//			{
//				for(size_t m = 0; m < vDimConsidered.size(); m++)
//				{
//					PartialVars[k](n, m) = m_vmdCovar[k](vDimConsidered[n], vDimConsidered[m]);
//					PartialInvVars[k](n, m) = m_vmdInvCovar[k](vDimConsidered[n], vDimConsidered[m]);
//				}
//			}
//		}
//		vPartialGconst.resize(m_nMixes, 0.0);
//		for(int n = 0; n < m_nMixes; n++)
//		{
//			vPartialGconst[n] = (double)vDimConsidered.size() * log(TWOPI)+log(PartialVars[n].Det());
//		}
//	}	
//    // obsolete function, don't use it
//	void CompPosterioriPartOld(const vector<double>& vec, vector<double>& vPost, vector<int>& vDimConsidered)
//	{
//		matrix<double> mvec;
//		vector< matrix<double> > PartialMeans;
//		vector< matrix<double> > PartialVars;
//		vector< matrix<double> > PartialInvVars;
//		vector<double> vPartialGconst;
//
//		size_t nsize = vDimConsidered.size();
//		mvec.SetSize(nsize, 1);
//		for(size_t n = 0; n < nsize; n++)
//		{
//			mvec(n, 0) = vec[vDimConsidered[n]];
//		}
//		
//		
//		GetPartialMeans(PartialMeans, vDimConsidered);
//		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
//		
//		vector<double> vdTemp(m_nMixes, 0.0);
//		double dTemp = 0.0;
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			matrix<double> x = (mvec - PartialMeans[n]);
//			matrix<double> xt = ~x;
//			matrix<double> mval = 0.5 * xt * PartialInvVars[n] * x;
//			vdTemp[n] = -0.5 * vPartialGconst[n] - mval(0,0);
//			vdTemp[n] = m_vMixWeights[n] * exp(vdTemp[n]);
//			dTemp += vdTemp[n];
//		}
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			vdTemp[n] /= dTemp;
//		}
//		vPost = vdTemp;
//	}
//
//	// Compute posterior probability from subspace
//	void CompPosterioriPart(const vector<double>& vec, vector<double>& vPost, vector<int>& vDimConsidered)
//	{
//		matrix<double> mvec;
//		vector< matrix<double> > PartialMeans;
//		vector< matrix<double> > PartialVars;
//		vector< matrix<double> > PartialInvVars;
//		vector<double> vPartialGconst;
//
//		size_t nsize = vDimConsidered.size();
//
//		//mvec.SetSize(nsize, 1);
//		//for(size_t n = 0; n < nsize; n++)
//		//{
//		//	mvec(n, 0) = vec[vDimConsidered[n]];
//		//}
//		vector<double> vdTemp(m_nMixes, 0.0);
//		double dTemp = LZERO;
//#ifndef DIAGONAL		
//		//mvec.SetSize(nsize, 1);
//		//for(size_t n = 0; n < nsize; n++)
//		//{
//		//	mvec(n, 0) = mvec0(vDimConsidered[n], 0);
//		//}
//		//dw2 = timeGetTime();
//		//cout<<"CompLikelihoodPart: Prepare vector takes "<<dw2-dw1<<" ms"<<endl;
//		//dw1 = timeGetTime();
//		//GetPartialMeans(PartialMeans, vDimConsidered);
//		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
//#endif
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			matrix<double> x(nsize, 1);
//			for(size_t d = 0; d < nsize; d++)
//				x(d,0) = (vec[vDimConsidered[d]] - m_vmdMean[n](vDimConsidered[d], 0));
//#ifndef DIAGONAL			
//			for(size_t j = 0; j < nsize; j++)
//			{
//				for(size_t i = j+1; i < nsize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(j, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[j]);
//				}
//			}
//			vdTemp[n] *= 2;
//			for(size_t i = 0; i < nsize; i++)
//			{
//				vdTemp[n] += x(i, 0) * x(i, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[i]);
//			}
//			vdTemp[n] += vPartialGconst[n];
//			vdTemp[n] *= -0.5;
//			vdTemp[n] += m_vlogMixWeights[n];
//#else
//			for(size_t d = 0; d < nsize; d++)
//				vdTemp[n] += m_vvddiaggconst[n][vDimConsidered[d]];
//			for(size_t i = 0; i < nsize; i++)
//			{
//				vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](vDimConsidered[i],vDimConsidered[i]);
//			}
//			vdTemp[n] *= -0.5;
//			vdTemp[n] += m_vlogMixWeights[n];
//#endif
//			dTemp = LogAdd(dTemp, vdTemp[n]);
//		}
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			vdTemp[n] -= dTemp;
//			vdTemp[n] = exp(vdTemp[n]);
//		}
//		vPost = vdTemp;
//	}
//	// Obsolete function, don't use it
//	double CompLikelihoodPartOld(const vector<double>& vec, vector<int>& vDimConsidered)
//	{
//		DWORD dw1 = 0, dw2 = 0;
//		//dw1 = timeGetTime();
//		matrix<double> mvec;
//		vector< matrix<double> > PartialMeans;
//		vector< matrix<double> > PartialVars;
//		vector< matrix<double> > PartialInvVars;
//		vector<double> vPartialGconst;
//		
//		size_t nsize = vDimConsidered.size();
//		mvec.SetSize(nsize, 1);
//		for(size_t n = 0; n < nsize; n++)
//		{
//			mvec(n, 0) = vec[vDimConsidered[n]];
//		}
//		
//		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
//		GetPartialMeans(PartialMeans, vDimConsidered);
//		vector<double> vdTemp(m_nMixes, 0.0);
//		double dTemp = 0;
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			matrix<double> x = mvec - PartialMeans[n];
//			matrix<double> xt = ~x;
//			matrix<double> mval = 0.5 * xt * PartialInvVars[n] * x;
//			vdTemp[n] = -0.5 * vPartialGconst[n] - mval(0,0);
//			vdTemp[n] = m_vMixWeights[n] * exp(vdTemp[n]);
//			dTemp += vdTemp[n];
//		}
//		return log(dTemp); //exp(dTemp);
//	}
//	// Compuete log likelihood from subspace
//	double CompLikelihoodPart(const vector<double>& vec, vector<int>& vDimConsidered)
//	{
//		DWORD dw1 = 0, dw2 = 0;
//		//dw1 = timeGetTime();
//		matrix<double> mvec;
//		vector< matrix<double> > PartialMeans;
//		vector< matrix<double> > PartialVars;
//		vector< matrix<double> > PartialInvVars;
//		vector<double> vPartialGconst;
//		
//		size_t nsize = vDimConsidered.size();
//#ifndef DIAGONAL		
//		//mvec.SetSize(nsize, 1);
//		//for(size_t n = 0; n < nsize; n++)
//		//{
//		//	mvec(n, 0) = vec[vDimConsidered[n]];
//		//}
//		//dw2 = timeGetTime();
//		//cout<<"CompLikelihoodPart: Prepare vector takes "<<dw2-dw1<<" ms"<<endl;
//		//dw1 = timeGetTime();
//		//GetPartialMeans(PartialMeans, vDimConsidered);
//		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
//#endif		
//		vector<double> vdTemp(m_nMixes, 0.0);
//		double dTemp = LZERO;
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			matrix<double> x(nsize, 1);
//			for(size_t d = 0; d < nsize; d++)
//				x(d,0) = (vec[vDimConsidered[d]] - m_vmdMean[n](vDimConsidered[d], 0));
//#ifndef DIAGONAL			
//			for(size_t j = 0; j < nsize; j++)
//			{
//				for(size_t i = j+1; i < nsize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(j, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[j]);
//				}
//			}
//			vdTemp[n] *= 2;
//			for(size_t i = 0; i < nsize; i++)
//			{
//				vdTemp[n] += x(i, 0) * x(i, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[i]);
//			}
//			vdTemp[n] += vPartialGconst[n];
//			vdTemp[n] *= -0.5;
//			vdTemp[n] += m_vlogMixWeights[n];
//#else
//			for(size_t d = 0; d < nsize; d++)
//				vdTemp[n] += m_vvddiaggconst[n][vDimConsidered[d]];
//			for(size_t i = 0; i < nsize; i++)
//			{
//				vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](vDimConsidered[i],vDimConsidered[i]);
//			}
//			vdTemp[n] *= -0.5;
//			vdTemp[n] += m_vlogMixWeights[n];
//#endif
//			dTemp = LogAdd(dTemp, vdTemp[n]);
//		}
//		return dTemp; //exp(dTemp);
//	}
//	// Compute log likelihood from subspace
//	double CompLikelihoodPart(const matrix<double>& mvec0, vector<int>& vDimConsidered)
//	{
//		matrix<double> mvec;
//		vector< matrix<double> > PartialMeans;
//		vector< matrix<double> > PartialVars;
//		vector< matrix<double> > PartialInvVars;
//		vector<double> vPartialGconst;
//
//		size_t nsize = vDimConsidered.size();
//		
//#ifndef DIAGONAL		
//		//mvec.SetSize(nsize, 1);
//		//for(size_t n = 0; n < nsize; n++)
//		//{
//		//	mvec(n, 0) = mvec0(vDimConsidered[n], 0);
//		//}
//		//dw2 = timeGetTime();
//		//cout<<"CompLikelihoodPart: Prepare vector takes "<<dw2-dw1<<" ms"<<endl;
//		//dw1 = timeGetTime();
//		//GetPartialMeans(PartialMeans, vDimConsidered);
//		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
//#endif
//		vector<double> vdTemp(m_nMixes, 0.0);
//		double dTemp = LZERO;
//		for(size_t n = 0; n < m_nMixes; n++)
//		{
//			matrix<double> x(nsize, 1);
//			for(size_t d = 0; d < nsize; d++)
//				x(d,0) = (mvec0(vDimConsidered[d], 0) - m_vmdMean[n](vDimConsidered[d], 0));
//#ifndef DIAGONAL			
//			for(size_t j = 0; j < nsize; j++)
//			{
//				for(size_t i = j+1; i < nsize; i++)
//				{
//					vdTemp[n] += x(i, 0) * x(j, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[j]);
//				}
//			}
//			vdTemp[n] *= 2;
//			for(size_t i = 0; i < nsize; i++)
//			{
//				vdTemp[n] += x(i, 0) * x(i, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[i]);
//			}
//			vdTemp[n] += vPartialGconst[n];
//			vdTemp[n] *= -0.5;
//			vdTemp[n] += m_vlogMixWeights[n];
//#else
//			for(size_t d = 0; d < nsize; d++)
//				vdTemp[n] += m_vvddiaggconst[n][vDimConsidered[d]];
//			for(size_t i = 0; i < nsize; i++)
//			{
//				vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](vDimConsidered[i],vDimConsidered[i]);
//			}
//			vdTemp[n] *= -0.5;
//			vdTemp[n] += m_vlogMixWeights[n];
//#endif
//			dTemp = LogAdd(dTemp, vdTemp[n]);
//		}
//		return dTemp;
//	}
//	int m_nVecSize;
//	int m_nMixes;
//	string m_szCovType;
//	string m_szFeatureType;
//	vector<double> m_vMixWeights;
//	vector< matrix<double> > m_vmdMean;
//	vector< matrix<double> > m_vmdCovar;
//	vector< matrix<double> > m_vmdInvCovar;
//	vector<double> m_vdgconst;
//	vector<double> m_vlogMixWeights;
//	vector< vector<double> > m_vvddiaggconst;
//};
//

class CVideoFE
{
public:
	// CVideoFE(void);
	// ~CVideoFE(void);
	// HRESULT GMMAdaptTool::AccBackground(const CvArr* src, CvArr*& dst);
	// void AggMerge(vector<CvRect>& vRect, vector<CvRect>& vAggRect, int& nROI);
	// HRESULT AvgBackground(const CvArr* src, CvArr*& dst, unsigned int nnum);
	// void AutoLabelRGBOnly(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio);
	// void AutoLabelPart(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio);
	// void AutoLabelPartSM(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio);
	// void CalcDist(vector<CvRect>& vRect, list< list<double> >& llRectDist, double& dMin, int& argminx, int& argminy);
	// bool CheckRectInt9r2(CvRect r1, CvRect r2);
	// HRESULT CloseImage(const CvArr* src, CvArr*& dst);
	// double CompRectDistNN(CvRect r1, CvRect r2);
	// double ComputeCvVectorCosine(const CvPoint& vec1, const CvPoint& vec2);
	// HRESULT GetKCurvaturePoint(CvPoint* pArr, int nSize, unsigned int i, int k, CvPoint& Pik, int MaxX, int MaxY);
	// HRESULT EstCOMandFingerLoc1(IplImage*& src, vector<double>& vec);
	// HRESULT EstCOMandFingerLoc2(IplImage*& src, vector<double>& vec, const vector<double>& lastvec);
	// int GMMBiasAdapt(vector< vector<double> >& vvData, matrix<double>& newmub, vector<GMM>& oldGMMs, vector<GMM>& newGMMs);
	// int GMMAdaptRcogAFrame(vector<GMM>& oldGMMs, vector< vector<double> >& vvData, vector<int>& vrecogclasses);
	// int LoadGMM(string szGMM, hash_map<string, GMM>& hGMMs);
	// HRESULT Makein01(const vector< vector<double> >& vvData, vector< vector<double> >& vvNewData);
	// HRESULT GetvvRangeAtnD(const vector< vector<double> >& vvData, unsigned int n, double& dbMin, double& dbMax);
	// int Test2ClassesGMMs(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses);
	// int Test2ClassesGMMsFull(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses);
	// bool UpdateDistRect(vector<CvRect>& vRect, list< list<double> >& llRectDist, double& dMin, int& argminx, int& argminy);
	// HRESULT GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vector<IplImage*> vpImgs, int& nDetStopBkGnd, int& nStop);
	// HRESULT OfflineFEImgSeq2_2_3_part(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs);
	// HRESULT OfflineFEImgSeq2_2_3_part_cheat(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs, LoadJPGLabelInfo& lab);
	// HRESULT WriteFeatureAsHTKMFCC(const vector< vector<double> >& vvData, string szOutFileName);
	// HRESULT WriteFeatureAsHTKMFCCAndItsHalves(const vector< vector<double> >& vvData, string szOutFileName, string szOutFileName_uhalf, string szOutFileName_lhalf);
	// int MkdirHier(string sztmp);
	// vector<int> vDimConsidered;
	// int m_nHandModelIndex;
	// int m_nBkGndModelIndex;
	// int m_nK;
	// double m_duth;
	// //==================================================================================================================================================================================================//
	// // New functions
	// void MakeDebugOutputName(string szPathName0, string szSubFolder, string szKeyName, string szExt, string& szOutFileName);
	// void SavecvImageForDebug(string szPathName0, string szSubFolder, string szKeyName, int n, int nLength, string szExt, IplImage* pImg);
	// // Calculate the mean vector of a vector< vector<double> > data matrix by trating each vector<double> as a data point
	// HRESULT CalcvvMean(const vector< vector<double> >& vvData, vector<double>& vMean);
	// // Get the minimum and maximum value of matrix stored in vector< vector<double> > form at n-th dimension
	// //HRESULT GetvvRangeAtnD(const vector< vector<double> >& vvData, unsigned int n, double& dbMin, double& dbMax);
	// // Make vector< vector<double> > data zero mean and range from [0, 1]
	// HRESULT MakeZeroMeanAndUnitRange(const vector< vector<double> >& vvData, vector< vector<double> >& vvNewData);
	// HRESULT WriteFeatureAsDEBUGTRAJ(const vector< vector<double> >& vvData, string szOutFileName1, string szOutFileName2);
	// CVideoFE(void);
	// ~CVideoFE(void);
	// HRESULT GMMAdaptTool::AccBackground(const CvArr* src, CvArr*& dst);
	// void AggMerge(vector<CvRect>& vRect, vector<CvRect>& vAggRect, int& nROI);
	// HRESULT AvgBackground(const CvArr* src, CvArr*& dst, unsigned int nnum);
	// void AutoLabelRGBOnly(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio);
	// void AutoLabelPart(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio);
	// void AutoLabelPartSM(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio);
	// void CalcDist(vector<CvRect>& vRect, list< list<double> >& llRectDist, double& dMin, int& argminx, int& argminy);
	// bool CheckRectInt9r2(CvRect r1, CvRect r2);
	// HRESULT CloseImage(const CvArr* src, CvArr*& dst);
	// double CompRectDistNN(CvRect r1, CvRect r2);
	// double ComputeCvVectorCosine(const CvPoint& vec1, const CvPoint& vec2);
	// HRESULT GetKCurvaturePoint(CvPoint* pArr, int nSize, unsigned int i, int k, CvPoint& Pik, int MaxX, int MaxY);
	// HRESULT EstCOMandFingerLoc1(IplImage*& src, vector<double>& vec);
	// HRESULT EstCOMandFingerLoc2(IplImage*& src, vector<double>& vec, const vector<double>& lastvec);
	// int GMMBiasAdapt(vector< vector<double> >& vvData, matrix<double>& newmub, vector<GMM>& oldGMMs, vector<GMM>& newGMMs);
	// int GMMAdaptRcogAFrame(vector<GMM>& oldGMMs, vector< vector<double> >& vvData, vector<int>& vrecogclasses);
	// int LoadGMM(string szGMM, hash_map<string, GMM>& hGMMs);
	// HRESULT Makein01(const vector< vector<double> >& vvData, vector< vector<double> >& vvNewData);
	// HRESULT GetvvRangeAtnD(const vector< vector<double> >& vvData, unsigned int n, double& dbMin, double& dbMax);
	// int Test2ClassesGMMs(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses);
	// int Test2ClassesGMMsFull(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses);
	// bool UpdateDistRect(vector<CvRect>& vRect, list< list<double> >& llRectDist, double& dMin, int& argminx, int& argminy);
	// HRESULT GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vector<IplImage*> vpImgs, int& nDetStopBkGnd, int& nStop);
	// HRESULT OfflineFEImgSeq2_2_3_part(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs);
	// HRESULT OfflineFEImgSeq2_2_3_part_cheat(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs, LoadJPGLabelInfo& lab);
	// HRESULT WriteFeatureAsHTKMFCC(const vector< vector<double> >& vvData, string szOutFileName);
	// HRESULT WriteFeatureAsHTKMFCCAndItsHalves(const vector< vector<double> >& vvData, string szOutFileName, string szOutFileName_uhalf, string szOutFileName_lhalf);
	// int MkdirHier(string sztmp);
	// //==================================================================================================================================================================================================//
	// // New functions
	// void MakeDebugOutputName(string szPathName0, string szSubFolder, string szKeyName, string szExt, string& szOutFileName);
	// void SavecvImageForDebug(string szPathName0, string szSubFolder, string szKeyName, int n, int nLength, string szExt, IplImage* pImg);
	// // Calculate the mean vector of a vector< vector<double> > data matrix by trating each vector<double> as a data point
	// HRESULT CalcvvMean(const vector< vector<double> >& vvData, vector<double>& vMean);
	// // Get the minimum and maximum value of matrix stored in vector< vector<double> > form at n-th dimension
	// //HRESULT GetvvRangeAtnD(const vector< vector<double> >& vvData, unsigned int n, double& dbMin, double& dbMax);
	// // Make vector< vector<double> > data zero mean and range from [0, 1]
	// HRESULT MakeZeroMeanAndUnitRange(const vector< vector<double> >& vvData, vector< vector<double> >& vvNewData);
	// HRESULT WriteFeatureAsDEBUGTRAJ(const vector< vector<double> >& vvData, string szOutFileName1, string szOutFileName2);

	// int LoadGMM(string szGMM, hash_map<string, GMM>& hGMMs)
	// {
		// ifstream fin;
		// fin.open(szGMM.c_str());
		// string szBuf;
		// if(fin.is_open())
		// {
			// int nVecSize = 0;
			// while(getline(fin, szBuf))
			// {
				// string strName = "";
				// string szFeatureType = "";
				// int nfindtilt = -1;
				// int nMixes = 0;

				// if(szBuf.npos != (nfindtilt = szBuf.find("~o")))
				// {
					
					// while(1)
					// {
						// getline(fin, szBuf);
						// int nfoundpos = -1;
						// if(szBuf.npos != (nfoundpos = szBuf.find("<VECSIZE>")))
						// {
							// int npos01 = szBuf.find_last_of("<");
							// int npos02 = szBuf.find_last_of(">");
							// string szVecSize = szBuf.substr(nfoundpos + 9, npos01 - nfoundpos - 9);
							// nVecSize = atoi(szVecSize.c_str());
							// szFeatureType = szBuf.substr(npos01 + 1, npos02 - npos01 - 1);
							// break;
						// }
					// }
				// }
				// else if(szBuf.npos != (nfindtilt = szBuf.find("~h")))
				// {
					// strName = szBuf.substr(nfindtilt+2, szBuf.length() - nfindtilt - 2);
					// int npos1 = strName.find_first_of('\"');
					// int npos2 = strName.find_last_of('\"');
					// strName = strName.substr(npos1+1, npos2-npos1-1);
					// getline(fin, szBuf);
					// while(1)
					// {
						// int nfoundpos = -1;
						// if(szBuf.npos != szBuf.find('<') && szBuf.npos != szBuf.find('>'))
						// {
							// std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
							// if(szBuf.npos != szBuf.find("<BEGINHMM>"))
							// {
							// }
							// else if(szBuf.npos != (nfoundpos = szBuf.find("<VECSIZE>")))
							// {
								// int npos01 = szBuf.find_last_of("<");
								// int npos02 = szBuf.find_last_of(">");
								// string szVecSize = szBuf.substr(nfoundpos + 9, npos01 - nfoundpos - 9);
								// nVecSize = atoi(szVecSize.c_str());
								// szFeatureType = szBuf.substr(npos01 + 1, npos02 - npos01 - 1);
							// }
							// else if(szBuf.npos != (nfoundpos = szBuf.find("<NUMMIXES>")))
							// {
								// int npos01 = szBuf.length() - 1;
								// string szMixes = szBuf.substr(nfoundpos + 10, npos01 - nfoundpos - 9);
								// nMixes = atoi(szMixes.c_str());
								// hGMMs[strName].InitGMMMeans(nVecSize, nMixes);
								// int nMixIndex = 0;
								// int nDims = 0;
								// double dMixWeight;
								// for(int kk = 0; kk < nMixes; kk++)
								// {
									// fin>>szBuf>>nMixIndex>>dMixWeight;
									// nMixIndex--;
									// fin>>szBuf>>nDims;
									// hGMMs[strName].m_vMixWeights[nMixIndex] = dMixWeight;
									// hGMMs[strName].m_vlogMixWeights[nMixIndex] = log(dMixWeight);
									// std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
									// if(szBuf.npos != szBuf.find("<MEAN>"))
									// {
										// //hGMMs[strName].InitGMMMeans(nVecSize, nMixes);
										// for(int d = 0; d < nDims; d++)
										// {
											// double db1 = 0.0;
											// fin>>db1;
											// hGMMs[strName].m_vmdMean[nMixIndex](d, 0) = db1;
										// }
										// fin>>szBuf>>nDims;
										// std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
										// if(szBuf.npos != szBuf.find("<VARIANCE>"))
										// {
											// hGMMs[strName].InitGMMVars(nVecSize, nMixIndex);
											// for(int d = 0; d < nDims; d++)
											// {
												// fin>>hGMMs[strName].m_vmdCovar[nMixIndex](d,d);
												// hGMMs[strName].m_vmdInvCovar[nMixIndex](d,d) = 1.0 / hGMMs[strName].m_vmdCovar[nMixIndex](d,d);
											// }
										// }
										// else if(szBuf.npos != szBuf.find("<INVCOVAR>"))
										// {
											// hGMMs[strName].InitGMMVars(nVecSize, nMixIndex, "full");
// #ifdef DIAGONAL
// #undef DIAGONAL
// #endif
											// for(int d = 0; d < nDims; d++)
											// {
												// for(int dd = d; dd < nDims; dd++)
												// {
													// fin>>hGMMs[strName].m_vmdInvCovar[nMixIndex](d,dd);
													// hGMMs[strName].m_vmdInvCovar[nMixIndex](dd, d) = hGMMs[strName].m_vmdInvCovar[nMixIndex](d,dd);
												// }
											// }
											// hGMMs[strName].m_vmdCovar[nMixIndex] = !hGMMs[strName].m_vmdInvCovar[nMixIndex];
										// }
										// double dbgconst = 0.0;
										// fin>>szBuf>>dbgconst;
									// }
									// else if(szBuf.npos != szBuf.find("<VARIANCE>"))
									// {
										// hGMMs[strName].InitGMMVars(nVecSize, nMixIndex);
										// for(int d = 0; d < nDims; d++)
										// {
											// fin>>hGMMs[strName].m_vmdCovar[nMixIndex](d,d);
											// hGMMs[strName].m_vmdInvCovar[nMixIndex](d,d) = 1.0 / hGMMs[strName].m_vmdCovar[nMixIndex](d,d);
										// }
										// fin>>szBuf>>nDims;
										// std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
										// if(szBuf.npos != szBuf.find("<MEAN>"))
										// {
											// for(int d = 0; d < nDims; d++)
											// {
												// fin>>hGMMs[strName].m_vmdMean[nMixIndex](d, 0);
											// }
										// }
										// double dbgconst = 0.0;
										// fin>>szBuf>>dbgconst;
									// }
									// else if(szBuf.npos != szBuf.find("<INVCOVAR>"))
									// {
										// hGMMs[strName].InitGMMVars(nVecSize, nMixIndex, "full");
										// for(int d = 0; d < nDims; d++)
										// {
											// for(int dd = d; dd < nDims; dd++)
											// {
												// fin>>hGMMs[strName].m_vmdInvCovar[nMixIndex](d,dd);
												// hGMMs[strName].m_vmdInvCovar[nMixIndex](dd, d) = hGMMs[strName].m_vmdInvCovar[nMixIndex](d,dd);
											// }
										// }									
										// hGMMs[strName].m_vmdCovar[nMixIndex] = !hGMMs[strName].m_vmdInvCovar[nMixIndex];

										// fin>>szBuf>>nDims;
										// std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
										// if(szBuf.npos != szBuf.find("<MEAN>"))
										// {
											// for(int d = 0; d < nDims; d++)
											// {
												// fin>>hGMMs[strName].m_vmdMean[nMixIndex](d, 0);
											// }
										// }
										// double dbgconst = 0.0;
										// fin>>szBuf>>dbgconst;
									// }

								// }
							// }
						// }
						// getline(fin, szBuf);
						// if(szBuf.npos != szBuf.find("<ENDHMM>"))
						// {
							// break;
						// }
					// }
				// }
			// }
		// }
		// fin.close();
		// hash_map<string, GMM>::iterator itr;
		// for(itr = hGMMs.begin(); itr != hGMMs.end(); itr++)
		// {
			// string szKey = itr->first;
			// hGMMs[szKey].CompGconsts();
		// }
		// return 1;

	// }
	
	CVideoFE(void)
	{
		vDimConsidered.resize(7);
		for(int n = 0; n < 7; n++)
			vDimConsidered[n] = n+1;
		m_nHandModelIndex = 0;
		m_nBkGndModelIndex = 0;
		m_nK = 20;
		m_duth = 0.75;
		m_szGMMPath = "";
	}

	~CVideoFE(void)
	{
	}

	// maximum likelihood decision for 2 classes
	int Test2ClassesGMMs(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses)
	{
		return (HRESULT) GMMAdaptTool::Test2ClassesGMMsPart(vvData, GMMs, recogClasses);
		// int N = (int)vvData.size();
		// recogClasses.resize(N, 0);
		// int I = (int)GMMs.size();
		// for(int n = 0; n < N; n++)
		// {
			// double maxni = GMMs[0].CompLikelihoodPart(vvData[n], vDimConsidered);
			// int argmaxni = 0;
			// for(int i = 1; i < I ;i++)
			// {
				// double dbi = GMMs[i].CompLikelihoodPart(vvData[n], vDimConsidered);
				// if(dbi > maxni)
				// {
					// maxni = dbi;
					// argmaxni = i;
				// }
			// }
			// recogClasses[n] = argmaxni;
		// }
		// return 1;
	}

	// maximum likelihood decision for 2 classes, without considering partial dimension
	int Test2ClassesGMMsFull(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses)
	{
		return GMMAdaptTool::Test2ClassesGMMsFull(vvData, GMMs, recogClasses);
	}

	// Stochastic matching on feature space
	int GMMAdaptRcogAFrame(vector<GMM>& oldGMMs, vector< vector<double> >& vvData, vector<int>& vrecogclasses)
	{
		return GMMAdaptTool::GMMAdaptRcogAFrame(oldGMMs, vvData, vrecogclasses);
	}

	void AutoLabelRGBOnly(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio)
	{
		return GMMAdaptTool::AutoLabelRGBOnly(szPath, vGMMs, mask, shrinkratio, m_nHandModelIndex);
	}
	void AutoLabel_BS0(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio, cv::Mat& ref)
	{
		return GMMAdaptTool::AutoLabel_BS0(szPath, vGMMs, mask, shrinkratio, m_nHandModelIndex, ref);
	}
	void Approx_AutoLabel_BS0(string szPath, bool*** vvvDecisionTable, cv::Mat& mask, int shrinkratio, cv::Mat& ref)
	{
		return GMMAdaptTool::Approx_AutoLabel_BS0(szPath, vvvDecisionTable, mask, shrinkratio, m_nHandModelIndex, ref);
	}
	// Get skin/nonskin color label, only consider RGB mean
	void AutoLabelPart(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio)
	{
		return GMMAdaptTool::AutoLabelPart(szPath, vGMMs, mask, shrinkratio, m_nHandModelIndex);
	}
	// Get skin/nonskin color label with stochastic feature matching
	void AutoLabelPartSM(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio)
	{
		return GMMAdaptTool::AutoLabelPartSM(szPath, vGMMs, mask, shrinkratio, m_nHandModelIndex);
	}
	double CompRectDistNN(CvRect r1, CvRect r2)
	{
		// compute the nearest distance between 2 rectangles
		double x1 = (double)r1.x;
		double y1 = (double)r1.y;
		double x2 = (double)(r1.x + r1.width);
		double y2 = (double)(r1.y + r1.height);
		double x3 = (double)r2.x;
		double y3 = (double)r2.y;
		double x4 = (double)(r2.x + r2.width);
		double y4 = (double)(r2.y + r2.height);
		double minx = MAX(x1,x3);
		double miny = MAX(y1,y3);
		double maxx = MIN(x2,x4);
		double maxy = MIN(y2,y4);
		double compx = ((minx > maxx)?(minx - maxx):0.0);
		double compy = ((miny > maxy)?(miny - maxy):0.0);
		double comp = ((compx > 0 && compy > 0)?sqrt(compx*compx+compy*compy):MAX(compx,compy));
		double area1 = (double)(r1.width * r1.height);
		double area2 =  (double)(r2.width * r2.height);
		comp /= sqrt(MIN(area1, area2));
		return comp;
	}

	void CalcDist(vector<CvRect>& vRect, list< list<double> >& llRectDist, double& dMin, int& argminx, int& argminy)
	{
		int N = (int)vRect.size();
		dMin = 1.7976931348623158e+308;
		llRectDist.resize(N-1);
		list< list<double> >::iterator llitr = llRectDist.begin();
		for(int i = 0; i < N-1; i++)
		{	
			for(int j = i+1; j < N; j++)
			{
				double dij = CompRectDistNN(vRect[i], vRect[j]);
				if(dij < dMin)
				{
					dMin = dij;
					argminx = j;
					argminy = i;
				}
				(*llitr).push_back(dij);
			}
			llitr++;
		}
	}

	bool CheckRectInt9r2(CvRect r1, CvRect r2)
	{
		// check if 4x rect2 intersect r1
		double x1 = (double)r1.x;
		double y1 = (double)r1.y;
		double x2 = (double)(r1.x + r1.width);
		double y2 = (double)(r1.y + r1.height);

		double x3 = (double)r2.x - (double)r2.width;
		double y3 = (double)r2.y - (double)r2.height;
		double x4 = x3 + 3.0 * (double)r2.width;
		double y4 = y3 + 3.0 * (double)r2.height;

		double minx = MAX(x1,x3);
		double miny = MAX(y1,y3);
		double maxx = MIN(x2,x4);
		double maxy = MIN(y2,y4);
		if((minx>maxx)|(miny>maxy))
		{
			return false;
		}
		return true;
	}

	bool UpdateDistRect(vector<CvRect>& vRect, list< list<double> >& llRectDist, double& dMin, int& argminx, int& argminy)
	{
		int N = (int)vRect.size();
		if(N == 1)
		{
			return false;
		}
		dMin = DBL_MAX; //1.7976931348623158e+308;
		//llRectDist.resize(N-1);
		list< list<double> >::iterator llitr = llRectDist.begin();
		int i = MIN(argminx, argminy);
		int j = MAX(argminx, argminy);
		double Areai = (double)(vRect[i].width * vRect[i].height);
		double Areaj = (double)(vRect[j].width * vRect[j].height);
		CvRect LRect = vRect[i];
		CvRect SRect = vRect[j];
		if(Areai < Areaj)
		{
			LRect = vRect[j];
			SRect = vRect[i];
		}
		bool bCanMerge = CheckRectInt9r2(LRect, SRect);
		if(bCanMerge)
		{
			CvRect Rectij = ::cvMaxRect(&LRect, &SRect);
			// Update the ith column of 0 to i-1 th rows,
			for(int m = 0; m < i; m++)
			{
				list<double>::iterator itrmn = (*llitr).begin();
				int mi = 0;
				while(mi < i - m - 1)
				{
					itrmn++;
					mi++;
				}
				(*itrmn) = CompRectDistNN(Rectij, vRect[m]);
				llitr++;
			}
			// Update ith row
			list<double>::iterator itrmn = (*llitr).begin();
			for(int n = i+1; n < N; n++)
			{
				if(n != j)
				{
					(*itrmn) = CompRectDistNN(Rectij, vRect[n]);
				}
				itrmn++;
			}
			// Delete jth row
			for(int k = i; k < j; k++)
			{
				llitr++;
			}
			int N1 = N - 1;
			if(llitr != llRectDist.end())
			{
				llRectDist.erase(llitr);
			}
			// Delete jth column for 0 to j-1 th row
			llitr = llRectDist.begin();
			for(int k = 0; k < j; k++)
			{
				itrmn = (*llitr).begin();
				int mi = 0;
				while(mi < j - k - 1)
				{
					itrmn++;
					mi++;
				}
				(*llitr).erase(itrmn);
				if(0 == (*llitr).size())
				{
					llRectDist.erase(llitr);
				}
				else
				{
					llitr++;
				}
			}
			// Find the new minimum distance
			llitr = llRectDist.begin();
			for(int n = 0; n < N1-1; n++)
			{	
				itrmn = (*llitr).begin();
				for(int m = n+1; m < N1; m++)
				{
					double dnm = (*itrmn);
					if(dnm < dMin)
					{
						dMin = dnm;
						argminx = m;
						argminy = n;
					}
					itrmn++;
				}
				llitr++;
			}
			// Update vRect
			vRect[i] = Rectij;
			for(int k = j; k < N1; k++)
			{
				vRect[k] = vRect[k+1];
			}
			vRect.pop_back();
		}
		return bCanMerge;
	}

	void AggMerge(vector<CvRect>& vRect, vector<CvRect>& vAggRect, int& nROI)
	{
		int N = (int) vRect.size();
		list< list<double> > llRectDist;
		double dMin = 1.7976931348623158e+308;
		int argminx = 0;
		int argminy = 0;
		for(int n = 0; n < N; n++)
		{
			vAggRect.push_back(vRect[n]);
		}
		CalcDist(vAggRect, llRectDist, dMin, argminx, argminy);
		while(UpdateDistRect(vAggRect, llRectDist, dMin, argminx, argminy))
		{
		}
		double dArea = 0.0;
		nROI = 0;
		for(int n = 0; n < (int)vAggRect.size(); n++)
		{
			double dArean = (double)(vAggRect[n].width * vAggRect[n].height);
			if(dArean > dArea)
			{
				dArea = dArean;
				nROI = n;
			}
		}

	}

	HRESULT GetKCurvaturePoint(CvPoint* pArr, int nSize, unsigned int i, int k, CvPoint& Pik, int MaxX, int MaxY)
	{
		if(pArr == 0)
		{
			return E_POINTER;
		}
		if((unsigned int)nSize <= i)
		{
			return ERROR_INVALID_ADDRESS;
		}
		CvPoint Pi = pArr[i];

		if(Pi.x >= MaxX - 2 || Pi.y >= MaxY - 2)
		{
			// Boundary point, set its cosine to -1 outside this function
			return ERROR_INVALID_INDEX;
		}
		if(k == 0)
		{
			// k must be integer other than 0
			return ERROR_INVALID_INDEX;
		}
			int nnk = 1;
		int nk = (k > 0)?i+1:i-1;
		nk = nk % nSize;
		if(nk < 0)
			nk += nSize;
		int ak = abs(k);
		bool bTouchBoundary = false;
		while(nnk <= ak)
		{
			CvPoint Pnk = pArr[nk];
			if(Pnk.x >= MaxX - 2 || Pnk.y >= MaxY - 2)
			{
				Pik.x = ak / nnk * (Pnk.x - Pi.x) + Pi.x;
				Pik.y = ak / nnk * (Pnk.y - Pi.y) + Pi.y;
				bTouchBoundary = true;
				break;
			}
			nnk++;
			nk = (k > 0)?nk+1:nk-1;
			nk = nk % nSize;
			if(nk < 0)
				nk += nSize;		
		}
		if(bTouchBoundary == false)
		{
			int nIndex = i + k;
			nIndex %= nSize;
			while(nIndex < 0)
			{
				nIndex += nSize;
			}
			Pik.x = pArr[nIndex].x;
			Pik.y = pArr[nIndex].y;
		}
		return S_OK;
	}

	// return the cosine value between 2 vectors
	double ComputeCvVectorCosine(const CvPoint& vec1, const CvPoint& vec2)
	{
		double norm1 = sqrt((double)(vec1.x * vec1.x + vec1.y * vec1.y));
		double norm2 = sqrt((double)(vec2.x * vec2.x + vec2.y * vec2.y));
		double dot = (double)(vec1.x * vec2.x + vec1.y * vec2.y);
		double dbcos = dot/(norm1 * norm2);
		return dbcos;
	}	
	HRESULT COMPoints(Mat& src, vector<double>& vec)
	{
		cv::Mat srccpy;
		src.copyTo(srccpy);
		vector< vector<Point> > contours;
		vector<double> vcontour_area;
		cv::findContours(srccpy, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);
		vcontour_area.resize(contours.size(),0);
		if(0 == contours.size())
		{
			vec.resize(2,-1.0);
			return S_FALSE;
		}
		double dbmaxarea = 0;
		for(size_t nn = 0; nn < contours.size(); ++nn)
		{
			vcontour_area[nn] = cv::contourArea(cv::Mat(contours[nn]));
			if(vcontour_area[nn] > dbmaxarea)
				dbmaxarea = vcontour_area[nn];
		}
		double area_tr = dbmaxarea*0.1;
		area_tr = (area_tr>10)?area_tr:10;
		vector<size_t> valid_idx;
		for(size_t nn = 0; nn < contours.size(); ++nn)
		{
			if(vcontour_area[nn]>area_tr)
			{
				valid_idx.push_back(nn);
			}
		}
		vec.resize(2,0.0);
		double dbcount = 0.0;
		for(size_t nn = 0; nn < src.cols; ++nn)
		{
			for(size_t mm = 0; mm < src.rows; ++mm)
			{
				CvMat tmpmat = src.operator CvMat();
				double dbVal = cvGetReal2D(&tmpmat, mm,nn);
				if(dbVal != 0)
				{
					Point2f pt((float)nn, (float)mm);
					for(size_t nn1 = 0; nn1 < valid_idx.size(); ++nn1)
					{
						double resnn1 = pointPolygonTest(cv::Mat(contours[valid_idx[nn1]]), pt, false);
						if(resnn1 >=0)
						{
							vec[0] += (double)nn;
							vec[1] += (double)mm;
							dbcount += 1.0;
							break;
						}
					}
				}
			}
		}
		vec[0]/=dbcount;
		vec[1]/=dbcount;
		return S_OK;
	}
	// Estimate center of mass (vec[0], vec[1]) and finger location (vec[2], vec[3])
	HRESULT EstCOMandFingerLoc1(IplImage*& src, vector<double>& vec)
	{
		CvMemStorage* storage = cvCreateMemStorage(0);
		CvMemStorage* chstorage = cvCreateMemStorage(0);
		CvSeq* first_contour = 0;
		double dbAngleThreshold = 0.1;
		int M = src->height;
		int N = src->width;
		int mode = CV_RETR_EXTERNAL;
		IplImage* srcpy = cvCreateImage(cvSize(src->width, src->height), src->depth, src->nChannels);
		cvCopy(src, srcpy);
		cvFindContours( srcpy, storage, &first_contour, sizeof(CvContour), 
			  mode, CV_CHAIN_APPROX_NONE/*CV_CHAIN_APPROX_SIMPLE*/, cvPoint(0,0));

		if(first_contour == 0) // No contour found
		{
			cvReleaseMemStorage(&storage);
			cvReleaseImage(&srcpy);
			return ERROR_OUT_OF_PAPER;
		}
		double dMaxArea = fabs(cvContourArea(first_contour));
		CvSeq* Maxcontour = first_contour;
		CvSeq* contour = first_contour->h_next;
		
		// For debug only
		/*
		if(m_pImg4DrawContour != 0)
		{
			cvDrawContours(m_pImg4DrawContour, 
				first_contour, 
				CV_RGB(0, 255, 255), 
				CV_RGB(255, 0, 255), 0, 2 
				);
		}
		*/

		while(contour != NULL)
		{
			double dArea = fabs(cvContourArea(contour));
			// For debug only
			/*
			if(m_pImg4DrawContour != 0)
			{
				cvDrawContours(m_pImg4DrawContour, 
					contour, 
					CV_RGB(0, 255, 255), 
					CV_RGB(255, 0, 255), 0, 2 
					);
			}
			*/
			if(dArea > dMaxArea)
			{
				Maxcontour = contour;
				dMaxArea = dArea;
			}
			contour = contour->h_next;
		}
		CvPoint* Arr = new CvPoint[Maxcontour->total];
		cvCvtSeqToArray(Maxcontour, Arr);
		int* pnCH = new int [Maxcontour->total];
		int nhullsize = 0;
		cvConvexHull(Arr, Maxcontour->total, NULL, CV_CLOCKWISE, pnCH, &nhullsize);

		CvPoint* Arrch = new CvPoint[nhullsize];
		for(int n = 0; n < nhullsize; n++)
		{
			Arrch[n].x = Arr[pnCH[n]].x;
			Arrch[n].y = Arr[pnCH[n]].y;
		}

		int my = src->height;
		int My = 0;
		int mx = src->width; 
		int Mx = 0;
		for(int n = 0; n < Maxcontour->total; n++)
		{
			if(Arr[n].y < my)
			{
				my = Arr[n].y;
			}
			if(Arr[n].x < mx)
			{
				mx = Arr[n].x;
			}
			if(Arr[n].x > Mx)
			{
				Mx = Arr[n].x;
			}
			if(Arr[n].y > My)
			{
				My = Arr[n].y;
			}
		}
		int tempx = Mx-mx+1;
		int tempy = My-my+1;
		cvSetImageROI(src, cvRect(mx, my, tempx, tempy));
		int tempx0 = 0;
		int tempy0 = 0;	
		int npts = 0;
		for(int m = 0; m < tempx; m++)
		{
			for(int n = 0; n < tempy; n++)
			{
				CvScalar val = cvGet2D(src, n, m);
				if(val.val[0] == 255.0)
				{
					tempx0 += m;
					tempy0 += n;
					npts++;
				}
			}
		}
		double cx = ((double)tempx0/(double)npts);
		double cy = ((double)tempy0/(double)npts);

		cvResetImageROI(src);
		cx += (double)mx;
		cy += (double)my;

		vector<double> vKK(Maxcontour->total);

		// For debug only
		/*
		if(m_pImg4DrawContour != 0)
		{
			cvDrawContours(m_pImg4DrawContour, 
				Maxcontour, 
				CV_RGB(0, 255, 255), 
				CV_RGB(255, 0, 255), 0 
				);
			for(int nch = 0; nch < nhullsize - 1; nch++)
			{
				cvLine(m_pImg4DrawContour, 
					Arrch[nch],
					Arrch[nch+1],
					CV_RGB(255, 0, 255));
			}

			cvLine(m_pImg4DrawContour, 
				Arrch[nhullsize - 1],
				Arrch[0],
				CV_RGB(255, 0, 255));
		}
		*/

		double MaxKK = -1.0;
		int MaxKKPos = 0;
		double MaxKK2 = -1.0;
		int MaxKKPos2 = 0;
		vector<double> vKK2(nhullsize);

		int nNumPts = nhullsize;
		double dbrK = (double)m_nK * (double)nNumPts / (double)Maxcontour->total;
		int nrK = (int)dbrK;
		if(nrK <= 2)
		{
			nrK = 2;
		}
		CvPoint* ArrToTrack = Arrch;
		for(int n = 0; n < nNumPts; n++)
		{
			CvPoint Pi = ArrToTrack[n];

			CvPoint Pip1;
			CvPoint Pim1;
			CvPoint Pip2;
			CvPoint Pim2;
			HRESULT hrp1 = GetKCurvaturePoint(ArrToTrack, nNumPts, n, 1, Pip1, N, M);
			HRESULT hrm1 = GetKCurvaturePoint(ArrToTrack, nNumPts, n, -1, Pim1, N, M);
			HRESULT hrp2 = GetKCurvaturePoint(ArrToTrack, nNumPts, n, 2, Pip2, N, M);
			HRESULT hrm2 = GetKCurvaturePoint(ArrToTrack, nNumPts, n, -2, Pim2, N, M);
			CvPoint Pip1ip2 = cvPoint(Pip2.x - Pip1.x, Pip2.y - Pip1.y);
			CvPoint Pip1i = cvPoint(Pi.x - Pip1.x, Pi.y - Pip1.y);
			CvPoint Pim1im2 = cvPoint(Pim2.x - Pim1.x, Pim2.y - Pim1.y);
			CvPoint Pim1i = cvPoint(Pi.x - Pim1.x, Pi.y - Pim1.y);
			double cosp1 = (double)(Pip1ip2.x * Pip1i.x + Pip1ip2.y * Pip1i.y)/
					( sqrt((double)(Pip1ip2.x * Pip1ip2.x + Pip1ip2.y * Pip1ip2.y)) * 
					sqrt((double)(Pip1i.x * Pip1i.x + Pip1i.y * Pip1i.y)));
			double cosm1 = (double)(Pim1im2.x * Pim1i.x + Pim1im2.y * Pim1i.y)/
					( sqrt((double)(Pim1im2.x * Pim1im2.x + Pim1im2.y * Pim1im2.y)) * 
					sqrt((double)(Pim1i.x * Pim1i.x + Pim1i.y * Pim1i.y)));

			CvPoint Pipk;
			CvPoint Pimk;
			
			HRESULT hrpk = GetKCurvaturePoint(ArrToTrack, nNumPts, n, nrK, Pipk, N, M);
			HRESULT hrmk = GetKCurvaturePoint(ArrToTrack, nNumPts, n, -1 * nrK, Pimk, N, M);
			if(hrpk == ERROR_INVALID_INDEX ||
				hrmk == ERROR_INVALID_INDEX ||
				hrp1 == ERROR_INVALID_INDEX ||
				hrm1 == ERROR_INVALID_INDEX ||
				hrp2 == ERROR_INVALID_INDEX ||
				hrm2 == ERROR_INVALID_INDEX)
			{
				vKK[n] = -1.0;
				vKK2[n] = -1.0;
			}
			else
			{
				CvPoint DPipk = cvPoint(Pipk.x - Pi.x, Pipk.y - Pi.y);
				CvPoint DPimk = cvPoint(Pimk.x - Pi.x, Pimk.y - Pi.y);
				// For debug only
				/*
				if(m_pImg4DrawContour != 0)
				{
					cvCircle(m_pImg4DrawContour,
						Pi,
						2,
						CV_RGB(255, 64, 255));
				}
				*/
				double dyi = (double)(Pi.y - my);
				vKK2[n] = (double)(DPipk.x * DPimk.x + DPipk.y * DPimk.y)/
					( sqrt((double)(DPipk.x * DPipk.x + DPipk.y * DPipk.y)) * 
					sqrt((double)(DPimk.x * DPimk.x + DPimk.y * DPimk.y)));
				if(vKK2[n] > MaxKK2 && Pi.x != Mx)
				{
					MaxKKPos2 = n;
					MaxKK2 = vKK2[n];
				}
				if(dyi >= 0 && dyi < (1 - m_duth) * (double)(M - my))
				{
					vKK[n] = (double)(DPipk.x * DPimk.x + DPipk.y * DPimk.y)/
						( sqrt((double)(DPipk.x * DPipk.x + DPipk.y * DPipk.y)) * 
						sqrt((double)(DPimk.x * DPimk.x + DPimk.y * DPimk.y)));

					if(fabs(cosp1) >= dbAngleThreshold &&
						fabs(cosm1) >= dbAngleThreshold &&
						Pi.x < N - 2
		/*				!(Pip1.x == Pip2.x || Pip1.y == Pip2.y) &&
						!(Pim1.x == Pim2.x || Pim1.y == Pim2.y) && 
						Pim1.x != Pi.x*/)
					{
						if(vKK[n] > MaxKK)
						{
							MaxKKPos = n;
							MaxKK = vKK[n];
						}
					}
				}
				else
				{
					vKK[n] = -1.0;
				}
			}
		}
		
		if(MaxKK < -0.9 /*&& (M - my) < 0.2 * M*/)
		{
			MaxKK = MaxKK2;
			MaxKKPos = MaxKKPos2;
		}

		double dfx = ArrToTrack[MaxKKPos].x;
		double dfy = ArrToTrack[MaxKKPos].y;
		// For debug only
		/*
		if(m_pImg4DrawContour != 0 )
		{
			cvCircle(m_pImg4DrawContour,cvPoint(cx, cy), 2, CV_RGB(255, 0, 0));
			cvCircle(m_pImg4DrawContour,cvPoint(dfx, dfy), 2, CV_RGB(0, 255, 0));
		}
		*/
		vec.resize(4);
		vec[0] = cx;
		vec[1] = cy;
		vec[2] = dfx;
		vec[3] = dfy;
		cvReleaseMemStorage(&storage);
		cvReleaseMemStorage(&chstorage);
		cvReleaseImage(&srcpy);
		delete [] Arr;
		delete [] Arrch;
		return S_OK;
	}

	// Estimate center of mass (vec[0], vec[1]) and finger location (vec[2], vec[3])
	HRESULT EstCOMandFingerLoc2(IplImage*& src, vector<double>& vec, const vector<double>& lastvec)
	{
		CvMemStorage* storage = cvCreateMemStorage(0);
		CvMemStorage* chstorage = cvCreateMemStorage(0);
		CvSeq* first_contour = 0;
		double dbAngleThreshold = 0.1;
		int M = src->height;
		int N = src->width;
		int mode = CV_RETR_EXTERNAL;
		IplImage* srcpy = cvCreateImage(cvSize(src->width, src->height), src->depth, src->nChannels);
		cvCopy(src, srcpy);
		cvFindContours( srcpy, storage, &first_contour, sizeof(CvContour), 
			  mode, CV_CHAIN_APPROX_NONE/*CV_CHAIN_APPROX_SIMPLE*/, cvPoint(0,0));
		if(first_contour == 0) // No contour found
		{
			cvReleaseMemStorage(&storage);
			cvReleaseImage(&srcpy);
			return ERROR_OUT_OF_PAPER;
		}
		double dMaxArea = fabs(cvContourArea(first_contour));
		CvSeq* Maxcontour = first_contour;
		CvSeq* contour = first_contour->h_next;
		// For debug only
		/*
		if(m_pImg4DrawContour != 0)
		{
			cvDrawContours(m_pImg4DrawContour, 
				first_contour, 
				CV_RGB(0, 255, 255), 
				CV_RGB(255, 0, 255), 0, 2 
				);
		}
		*/
		while(contour != NULL)
		{
			double dArea = fabs(cvContourArea(contour));
			// For debug only
			/*
			if(m_pImg4DrawContour != 0)
			{
				cvDrawContours(m_pImg4DrawContour, 
					contour, 
					CV_RGB(0, 255, 255), 
					CV_RGB(255, 0, 255), 0);
			}
			*/
			if(dArea > dMaxArea)
			{
				Maxcontour = contour;
				dMaxArea = dArea;
			}
			contour = contour->h_next;
		}

		// K-curvature Algorithm
		CvPoint* Arr = new CvPoint[Maxcontour->total];
		cvCvtSeqToArray(Maxcontour, Arr);
		int* pnCH = new int [Maxcontour->total];
		int nhullsize = 0;
		cvConvexHull(Arr, Maxcontour->total, NULL, CV_CLOCKWISE, pnCH, &nhullsize);

		CvPoint* Arrch = new CvPoint[nhullsize];
		for(int n = 0; n < nhullsize; n++)
		{
			Arrch[n].x = Arr[pnCH[n]].x;
			Arrch[n].y = Arr[pnCH[n]].y;
		}

		int my = src->height;
		int My = 0;
		int mx = src->width; 
		int Mx = 0;
		for(int n = 0; n < Maxcontour->total; n++)
		{
			if(Arr[n].y < my)
			{
				my = Arr[n].y;
			}
			if(Arr[n].x < mx)
			{
				mx = Arr[n].x;
			}
			if(Arr[n].x > Mx)
			{
				Mx = Arr[n].x;
			}
			if(Arr[n].y > My)
			{
				My = Arr[n].y;
			}
		}
		int tempx = Mx-mx+1;
		int tempy = My-my+1;
		cvSetImageROI(src, cvRect(mx, my, tempx, tempy));
		int tempx0 = 0;
		int tempy0 = 0;	
		int npts = 0;
		for(int m = 0; m < tempx; m++)
		{
			for(int n = 0; n < tempy; n++)
			{
				CvScalar val = cvGet2D(src, n, m);
				if(val.val[0] == 255.0)
				{
					tempx0 += m;
					tempy0 += n;
					npts++;
				}
			}
		}
		double cx = ((double)tempx0/(double)npts);
		double cy = ((double)tempy0/(double)npts);
		cvResetImageROI(src);
		cx += (double)mx;
		cy += (double)my;

		vector<double> vKK(Maxcontour->total);
		// For debug only
		/*
		if(m_pImg4DrawContour != 0)
		{
			cvDrawContours(m_pImg4DrawContour, 
				Maxcontour, 
				CV_RGB(0, 255, 255), 
				CV_RGB(255, 0, 255), 
				0);
			for(int nch = 0; nch < nhullsize - 1; nch++)
			{
				cvLine(m_pImg4DrawContour, 
					Arrch[nch],
					Arrch[nch+1],
					CV_RGB(255, 0, 255));
			}
			cvLine(m_pImg4DrawContour, 
				Arrch[nhullsize - 1],
				Arrch[0],
				CV_RGB(255, 0, 255));
		}
		*/
		double MaxKK = -1.0;
		int MaxKKPos = 0;
		double MaxKK2 = -1.0;
		int MaxKKPos2 = 0;
		vector<double> vKK2(nhullsize);

		int nNumPts = nhullsize;
		double dbrK = (double)m_nK * (double)nNumPts / (double)Maxcontour->total;
		int nrK = (int)dbrK;
		if(nrK <= 2)
		{
			nrK = 2;
		}
		CvPoint* ArrToTrack = Arrch;
		CvPoint mcMove = cvPoint((int)(cx - lastvec[0]), (int)(cy - lastvec[1]));
		CvPoint lastmcfg = cvPoint((int)(lastvec[2] - lastvec[0]), (int)(lastvec[3] - lastvec[1]));
		ofstream outf;

		for(int n = 0; n < nNumPts; n++)
		{
			CvPoint Pi = ArrToTrack[n];		
			CvPoint fgMove = cvPoint((int)(Pi.x - lastvec[2]), (int)(Pi.y - lastvec[3]));
			CvPoint mcfg = cvPoint((int)(Pi.x - cx), (int)(Pi.y - cy));
			double dMovecos = ComputeCvVectorCosine(mcMove, fgMove);
			double dmcfgcos = ComputeCvVectorCosine(lastmcfg, mcfg);

			if(dmcfgcos > 0/*mcMove.x * fgMove.x >= 0 /*&& 
				mcMove.y * fgMove.y >= 0 &&
				lastmcfg.x * mcfg.x >= 0 &&
				lastmcfg.y * mcfg.y >= 0*/)
			{
				CvPoint Pip1;
				CvPoint Pim1;
				CvPoint Pip2;
				CvPoint Pim2;
				HRESULT hrp1 = GetKCurvaturePoint(ArrToTrack, nNumPts, n, 1, Pip1, N, M);
				HRESULT hrm1 = GetKCurvaturePoint(ArrToTrack, nNumPts, n, -1, Pim1, N, M);
				HRESULT hrp2 = GetKCurvaturePoint(ArrToTrack, nNumPts, n, 2, Pip2, N, M);
				HRESULT hrm2 = GetKCurvaturePoint(ArrToTrack, nNumPts, n, -2, Pim2, N, M);
				CvPoint Pip1ip2 = cvPoint(Pip2.x - Pip1.x, Pip2.y - Pip1.y);
				CvPoint Pip1i = cvPoint(Pi.x - Pip1.x, Pi.y - Pip1.y);
				CvPoint Pim1im2 = cvPoint(Pim2.x - Pim1.x, Pim2.y - Pim1.y);
				CvPoint Pim1i = cvPoint(Pi.x - Pim1.x, Pi.y - Pim1.y);
				double cosp1 = (double)(Pip1ip2.x * Pip1i.x + Pip1ip2.y * Pip1i.y)/
						( sqrt((double)(Pip1ip2.x * Pip1ip2.x + Pip1ip2.y * Pip1ip2.y)) * 
						sqrt((double)(Pip1i.x * Pip1i.x + Pip1i.y * Pip1i.y)));
				double cosm1 = (double)(Pim1im2.x * Pim1i.x + Pim1im2.y * Pim1i.y)/
						( sqrt((double)(Pim1im2.x * Pim1im2.x + Pim1im2.y * Pim1im2.y)) * 
						sqrt((double)(Pim1i.x * Pim1i.x + Pim1i.y * Pim1i.y)));

				CvPoint Pipk;
				CvPoint Pimk;
				
				HRESULT hrpk = GetKCurvaturePoint(ArrToTrack, nNumPts, n, nrK, Pipk, N, M);
				HRESULT hrmk = GetKCurvaturePoint(ArrToTrack, nNumPts, n, -1 * nrK, Pimk, N, M);
				if(hrpk == ERROR_INVALID_INDEX ||
					hrmk == ERROR_INVALID_INDEX ||
					hrp1 == ERROR_INVALID_INDEX ||
					hrm1 == ERROR_INVALID_INDEX ||
					hrp2 == ERROR_INVALID_INDEX ||
					hrm2 == ERROR_INVALID_INDEX)
				{
					vKK[n] = -1.0;
					vKK2[n] = -1.0;
				}
				else
				{
					CvPoint DPipk = cvPoint(Pipk.x - Pi.x, Pipk.y - Pi.y);
					CvPoint DPimk = cvPoint(Pimk.x - Pi.x, Pimk.y - Pi.y);
					// For debug only
					/*
					if(m_pImg4DrawContour != 0)
					{
						cvCircle(m_pImg4DrawContour,
							Pi,
							2,
							CV_RGB(255, 64, 255));
					}
					*/
					double dyi = (double)(Pi.y - my);
					vKK2[n] = (double)(DPipk.x * DPimk.x + DPipk.y * DPimk.y)/
						( sqrt((double)(DPipk.x * DPipk.x + DPipk.y * DPipk.y)) * 
						sqrt((double)(DPimk.x * DPimk.x + DPimk.y * DPimk.y)));
					if(vKK2[n] > MaxKK2 && Pi.x != Mx)
					{
						MaxKKPos2 = n;
						MaxKK2 = vKK2[n];
					}
					if(dyi >= 0 && dyi < (1 - m_duth) * (double)(M - my))
					{
						vKK[n] = (double)(DPipk.x * DPimk.x + DPipk.y * DPimk.y)/
							( sqrt((double)(DPipk.x * DPipk.x + DPipk.y * DPipk.y)) * 
							sqrt((double)(DPimk.x * DPimk.x + DPimk.y * DPimk.y)));

						if(1/*fabs(cosp1) >= dbAngleThreshold &&
							fabs(cosm1) >= dbAngleThreshold &&
							Pi.x < N - 2*/ 
							/*
							!(Pip1.x == Pip2.x || Pip1.y == Pip2.y) &&
							!(Pim1.x == Pim2.x || Pim1.y == Pim2.y) && 
							Pim1.x != Pi.x*/)
						{
							if(vKK[n] > MaxKK)
							{
								MaxKKPos = n;
								MaxKK = vKK[n];
							}
						}
					}
					else
					{
						vKK[n] = -1.0;
					}
				}
			}
			else
			{
				vKK[n] = -1.0;
				vKK2[n] = -1.0;
			}
		}
		if(MaxKK < -0.9 /*&& (M - my) < 0.2 * M*/)
		{
			MaxKK = MaxKK2;
			MaxKKPos = MaxKKPos2;
		}
		double dfx = ArrToTrack[MaxKKPos].x;
		double dfy = ArrToTrack[MaxKKPos].y;
		// For debug only
		/*
		if(m_pImg4DrawContour != 0 )
		{
			cvCircle(m_pImg4DrawContour,cvPoint(cx, cy), 2, CV_RGB(255, 0, 0));
			cvCircle(m_pImg4DrawContour,cvPoint(dfx, dfy), 2, CV_RGB(0, 255, 0));
		}
		*/
		vec.resize(4);
		vec[0] = cx;
		vec[1] = cy;
		vec[2] = dfx;
		vec[3] = dfy;
		cvReleaseMemStorage(&storage);
		cvReleaseMemStorage(&chstorage);
		cvReleaseImage(&srcpy);
		delete [] Arr;
		delete [] Arrch;
		return S_OK;
	}

	// Morphological closing operation
	HRESULT CloseImage(const CvArr* src, CvArr*& dst)
	{
		if(dst == NULL)
		{
			if(CV_IS_MAT_HDR(src))
			{
				dst = cvCreateMat(((CvMat*)src)->rows, ((CvMat*)src)->cols, ((CvMat*)src)->type);
			}
			else if(CV_IS_IMAGE_HDR(src))
			{
				dst = cvCreateImage(cvSize(((IplImage*)src)->width, ((IplImage*)src)->height), ((IplImage*)src)->depth, ((IplImage*)src)->nChannels );
			}
		}
		// construct the structured element
		IplConvKernel* pelem = ::cvCreateStructuringElementEx(4, 4, 0, 0, CV_SHAPE_ELLIPSE);
		//Mat elem = cv::getStructuringElement(MORPH_ELLIPSE, Size(6, 6));
		//Mat temp = Mat();
		CvArr* temp = 0 ;
		if(CV_IS_MAT_HDR(src))
		{
			temp = cvCreateMat(((CvMat*)src)->rows, ((CvMat*)src)->cols, ((CvMat*)src)->type);
		}
		else if(CV_IS_IMAGE_HDR(src))
		{
			temp = cvCreateImage(cvSize(((IplImage*)src)->width, ((IplImage*)src)->height), ((IplImage*)src)->depth, ((IplImage*)src)->nChannels );
		}
		cvMorphologyEx(src, dst, temp, pelem, CV_MOP_CLOSE, 2); 
		cvReleaseStructuringElement(&pelem);
		if(CV_IS_MAT_HDR(src))
		{
			cvReleaseMat((CvMat**)&temp);
		}
		else if(CV_IS_IMAGE_HDR(src))
		{
			cvReleaseImage((IplImage**)&temp);
		}	
		return S_OK;
	}


	HRESULT OfflineFEImgSeq2_2_3_part_baseline(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		//for(n = nStart; n <= nStopBkGnd; n++)
		//{
		//	pcurr = vpImgs[n];
		//	if(n == nStart)
		//	{
		//		pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//		pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//		cvSet(pBkGnd, cvScalar(0.0));
		//		cvSet(pBkGndTmp, cvScalar(0.0));
		//	}
		//	pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//	cvConvert(pcurr, pcurr2);
		//	GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
		//	cvReleaseImage(&pcurr);
		//	cvReleaseImage(&pcurr2);
		//}
		//AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		//cvReleaseImage(&pBkGndTmp);
		//=============================================================================//

		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask2, 4);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask0, 32);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask, 16);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask2, 4);

			//AutoLabelPart(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPart(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPart(szImgNamen, vGMMs, refmask2, 4);
			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			IplImage* prefmask2 = &iplrefmask2;
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}
			//=============================================================================//
			vector<double> comvec;
			// ===============================================================================//
			// Output combined mask for durther processing
			// ===============================================================================//
			SaveMATForDebug(szPathName0, "mask_combine", szKeyName, n, 5, ".dat", prefmask2);
			// ===============================================================================//
			COMPoints(refmask2, comvec);
			if(2 <= comvec.size())
				if(comvec[0] >= 0 && comvec[1] >= 0)
					vvFeature.push_back(comvec);
			
			// pcurr = vpImgs[n];

			// pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			// pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			// //pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			// pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			// pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			// //pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			// pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			// cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			//cvConvert(pcurr2b, pcurr2b_tempb);
			// //cvAnd(pcurr, pcurr2b_tempb, pMask1);
			//cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			//IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			//IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			//cvSet(pMask1, CV_RGB(64, 128, 64));
			// //cv::Mat Mask1(pMask1);
			//cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			// //cv::Mat Tmp(pTmp);
			// IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			// IplImage* piplcurrtemp = &iplcurrtemp;
			// cvCopy(piplcurrtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			// SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			// cvReleaseImage(&pTmp);
			// cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			// cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			// CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			// cvConvert(pcurr2c, pImgb);
			// IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			// cvAnd(pcurr, pImgb, pMask2);
			// SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			// cvReleaseImage(&pMask2);
			 cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			 cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			// hrTrackFinger = S_OK;
			// if(vvFeature.size() == 0)
			// {
			//  	hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			// }
			// else
			// {
			// 	hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			// }

/* 			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// // simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// // Means the gesture ends in online FE				
				}
				// // When it comes to the frame before the last frame,
				// // we just assume the point stay at the last observation 
				// // all the time
				// //!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// // Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			} */
			refmask0.release();
			refmask.release();
			refmask2.release();
			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		// cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}

	HRESULT OfflineFEImgSeq2_2_3_part(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		//for(n = nStart; n <= nStopBkGnd; n++)
		//{
		//	pcurr = vpImgs[n];
		//	if(n == nStart)
		//	{
		//		pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//		pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//		cvSet(pBkGnd, cvScalar(0.0));
		//		cvSet(pBkGndTmp, cvScalar(0.0));
		//	}
		//	pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//	cvConvert(pcurr, pcurr2);
		//	GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
		//	cvReleaseImage(&pcurr);
		//	cvReleaseImage(&pcurr2);
		//}
		//AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		//cvReleaseImage(&pBkGndTmp);
		//=============================================================================//

		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask2, 4);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask0, 32);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask, 16);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask2, 4);

			//AutoLabelPart(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPart(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPart(szImgNamen, vGMMs, refmask2, 4);
			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			SavecvImageForDebug(szPathName0, "Clr_32_32", szKeyName, n, 5, ".jpg", &iplrefmask0);
			SavecvImageForDebug(szPathName0, "Clr_16_16", szKeyName, n, 5, ".jpg", &iplrefmask);
			SavecvImageForDebug(szPathName0, "Clr_4_4", szKeyName, n, 5, ".jpg", &iplrefmask2);
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}
			//=============================================================================//
		
			pcurr = vpImgs[n];

			pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			//cv::Mat Mask1(pMask1);
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			//cv::Mat Tmp(pTmp);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			cvAnd(pcurr, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();
			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		//cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}

	HRESULT OfflineFEImgSeq2_2_3_part_BS0(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//

		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask2, 4);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask0, 32, BkGnd);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask, 16, BkGnd);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask2, 4, BkGnd);

//#ifdef _MSC_VER
			//DWORD dw1 = ::timeGetTime();
//#endif
			AutoLabel_BS0(szImgNamen, vGMMs, refmask0, 32, BkGnd);
//#ifdef _MSC_VER
			//DWORD dw2 = ::timeGetTime();
//#endif
			AutoLabel_BS0(szImgNamen, vGMMs, refmask, 16, BkGnd);
//#ifdef _MSC_VER
//			DWORD dw3 = ::timeGetTime();
//#endif
			AutoLabel_BS0(szImgNamen, vGMMs, refmask2, 4, BkGnd);

//#ifdef _MSC_VER
//			DWORD dw4 = ::timeGetTime();
			//cout<<"Scale 32*32 takes "<<dw2-dw1<<" ms"<<endl;
			//cout<<"Scale 16*16 takes "<<dw3-dw2<<" ms"<<endl;
			//cout<<"Scale 4*4 takes "<<dw4-dw3<<" ms"<<endl;
//#endif
			//AutoLabelPart(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPart(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPart(szImgNamen, vGMMs, refmask2, 4);
			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}
			vector<double> comvec;
			//COMPoints(refmask2, comvec);
			//=============================================================================//
		
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask, 16, BkGnd);
			//cv::resize(refmask, refmask2, cv::Size(), 4, 4, INTER_NEAREST);
			//iplrefmask2 = refmask2.operator IplImage();
			pcurr = vpImgs[n];

			pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			//cv::Mat Mask1(pMask1);
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			//cv::Mat Tmp(pTmp);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			cvAnd(pcurr, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();
			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}



	// bFrameBatch: true: framewise, false: batch
	// nAdaptMode: 0: MAPOnly 1: ScoreFusion Only 2: MAP+Score Fusion
	// old version
	//HRESULT OfflineFEImgSeq2_2_3_part_BS0_Fusion(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode)
	//{
	//#ifdef _WIN32
	//	string delim = "\\";
	//#else
	//	string delim = "/";
	//#endif
	//	string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
	//	string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
	//	string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
	//	string szMLFName = vFEN[nID].szMLFName;
	//	string szLabelName = vFEN[nID].szLabelName;
	//	string szKeyName = vFEN[nID].szKeyName; 
	//	string szErrLog = vFEN[nID].szErrLog;
	//	size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
	//	string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

	//	IplImage* pcurr = 0;
	//	IplImage* pBkGndTmp = 0;
	//	IplImage* pBkGnd = 0;
	//	IplImage* pRef1 = 0;
	//	IplImage* pRef2 = 0;
	//	IplImage* pcurr2 = 0;
	//	IplImage* pcurr2c = 0;
	//	IplImage* pcurr2b = 0;
	//	//IplImage* pcurr2bcpy = 0;
	//	IplImage* pcurr2b_tempb = 0;
	//	//IplImage* pcurr2b_tempb2 = 0;
	//	IplImage* pcurr2b_tempb3 = 0;
	//	IplImage* pImgb = 0;
	//	IplImage* pImgbs = 0;
	//	IplImage iplrefmask0;
	//	IplImage iplrefmask;
	//	IplImage iplrefmask2;
	//	CvMemStorage* storage = 0;
	//	CvSeq* first_contour = 0;
	//	CvSeq* contour = 0;
	//	cv::Mat curr2b;
	//	cv::Mat curr2b2;
	//	vector< vector<double> > vvFeature;
	//	vector<double> vTmpTmp(4);
	//	vector<double> vTmp;
	//	vector<int> vHandVanishIndices;
	//	vector<CvRect> vRect;
	//	vector<CvRect> vAggRect;
	//	CvRect rectROI;
	//	CvRect rect1;
	//	CvRect recti;
	//	IplImage* pImgbSWM = 0;
	//	IplImage* pSWMask = 0;

	//	HRESULT hr = S_OK;
	//	HRESULT hrTrackFinger = S_OK;
	//	int nStopBkGnd = 0;
	//	int nStart = 0;
	//	int nStop = (int)vpImgs.size() - 1;
	//	int n = 0;
	//	int kk1 = 0;
	//	int kk2 = 0;
	//	int nn1 = 0;
	//	int nn2 = 0;
	//	int nk = 0;
	//	int nrectnum = 0;
	//	int nROIIndex = 0;
	//	
	//	
	//	double db1 = 0.0;
	//	double cxn_k = 0.0;
	//	double cyn_k = 0.0;
	//	double dxn_k = 0.0;
	//	double dyn_k = 0.0;
	//	double cxn = 0.0;
	//	double cyn = 0.0;;
	//	double dxn = 0.0;
	//	double dyn = 0.0;					
	//	double r = 0.0;
	//	//=============================================================================//
	//	// Begin feature extraction
	//	//=============================================================================//
	//	if(nStop <= 0)
	//	{
	//		return ERROR_EMPTY;
	//	}
	//	//=============================================================================//

	//	//=============================================================================//
	//	// Check where hand is inside the frame and outside 
	//	//=============================================================================//
	//	hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
	//	if(hr == ERROR_EMPTY)
	//	{
	//		return hr;
	//	}
	//	//=============================================================================//


	//	//=============================================================================//
	//	//Accumulate and average background images
	//	//=============================================================================//
	//	for(n = nStart; n <= nStopBkGnd; n++)
	//	{
	//		pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
	//		if(n == nStart)
	//		{
	//			pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
	//			pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
	//			cvSet(pBkGnd, cvScalar(0.0));
	//			cvSet(pBkGndTmp, cvScalar(0.0));
	//		}
	//		pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
	//		cvConvert(pcurr, pcurr2);
	//		GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
	//		cvReleaseImage(&pcurr);
	//		cvReleaseImage(&pcurr2);
	//	}
	//	AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
	//	cvReleaseImage(&pBkGndTmp);
	//	cv::Mat BkGnd(pBkGnd);
	//	//=============================================================================//

	//	//=============================================================================//
	//	// SMAP for frames
	//	//=============================================================================//
	//	vector<string> vszImagePathNames;
	//	vector<cv::Mat> vMasks1;
	//	vector<cv::Mat> vMasks2;
	//	vector<cv::Mat> vMasks3;
	//	vector<cv::Mat> vMasks01;
	//	for(n = 0/*nStopBkGnd + 1*/; n < vpImgs.size()/*= nStop*/; n++)
	//	{
	//		vszImagePathNames.push_back(ImgPathNames[n]);
	//	}
	//	if(m_nHandModelIndex == 0)
	//		m_nHandModelIndex = 1;
	//	//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion(vszImagePathNames, vGMMs, vMasks1, 32, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1);
	//	//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion(vszImagePathNames, vGMMs, vMasks2, 16, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1);
	//	//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion(vszImagePathNames, vGMMs, vMasks3, 4, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1);
	//	vector<cv::Mat> vScore0_1;
	//	vector<cv::Mat> vScore1_1;
	//	vector<cv::Mat> vScore0_2;
	//	vector<cv::Mat> vScore1_2;
	//	vector<cv::Mat> vScore0_3;
	//	vector<cv::Mat> vScore1_3;

	//	GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks1, 32, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_1, vScore1_1);
	//	GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks2, 16, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_2, vScore1_2);
	//	GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks3, 4, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_3, vScore1_3);

	//	//IplImage tmp1 = vMasks1[16].operator IplImage();
	//	//IplImage tmp01 = vMasks01[16].operator IplImage();

	//	//SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, 16, 5, "_o.jpg", &(tmp1));
	//	//SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, 16, 5, "_ws.jpg", &(tmp01));

	//	
	//	
	//	//=============================================================================//
	//	// Do feature extraction for those frames with hands
	//	//=============================================================================//
	//	for(n = nStopBkGnd + 1; n <= nStop; n++)
	//	{
	//		string szImgNamen = ImgPathNames[n];
	//		cv::Mat refmask, refmaskq, refmask2, refmask0;
	//		//=============================================================================//
	//		// AutoLabelPartSM: label image in 3 different block sizes with feature space 
	//		// stochastic matching:
	//		// AutoLabelPart: direct label with current GMM in 3 different block sizes
	//		// around 1.5% performance gain using AutoLabelPartSM
	//		//=============================================================================//
	//		//AutoLabelPartSM(szImgNamen, vGMMs, refmask0, 32);
	//		//AutoLabelPartSM(szImgNamen, vGMMs, refmask, 16);
	//		//AutoLabelPartSM(szImgNamen, vGMMs, refmask2, 4);
	//		//AutoLabel_BS0(szImgNamen, vGMMs, refmask0, 32, BkGnd);
	//		//AutoLabel_BS0(szImgNamen, vGMMs, refmask, 16, BkGnd);
	//		//AutoLabel_BS0(szImgNamen, vGMMs, refmask2, 4, BkGnd);
	//		refmask0 = vMasks1[n];// - nStopBkGnd - 1];
	//		refmask = vMasks2[n];// - nStopBkGnd - 1];
	//		refmask2 = vMasks3[n];// - nStopBkGnd - 1];
	//		IplImage iplrefmask0 = refmask0.operator IplImage();
	//		IplImage iplrefmask1 = refmask.operator IplImage();
	//		IplImage iplrefmask2 = refmask2.operator IplImage();
	//		IplImage* prefmask0 = &iplrefmask0;
	//		IplImage* prefmask1 = &iplrefmask1;
	//		IplImage* prefmask2 = &iplrefmask2;
	//		SaveMATForDebug(szPathName0, "mask_32_32", szKeyName, n, 5, ".dat", prefmask0);
	//		SaveMATForDebug(szPathName0, "mask_16_16", szKeyName, n, 5, ".dat", prefmask1);
	//		SaveMATForDebug(szPathName0, "mask_4_4", szKeyName, n, 5, ".dat", prefmask2);
	//		IplImage iplscore0_1 = vScore0_1[n].operator IplImage();
	//		IplImage iplscore1_1 = vScore1_1[n].operator IplImage();
	//		IplImage iplscore0_2 = vScore0_2[n].operator IplImage();
	//		IplImage iplscore1_2 = vScore1_2[n].operator IplImage();
	//		IplImage iplscore0_3 = vScore0_3[n].operator IplImage();
	//		IplImage iplscore1_3 = vScore1_3[n].operator IplImage();
	//		IplImage* piplscore0_1 = &iplscore0_1;
	//		IplImage* piplscore1_1 = &iplscore1_1;
	//		IplImage* piplscore0_2 = &iplscore0_2;
	//		IplImage* piplscore1_2 = &iplscore1_2;
	//		IplImage* piplscore0_3 = &iplscore0_3;
	//		IplImage* piplscore1_3 = &iplscore1_3;
	//		SaveMATForDebug(szPathName0, "mask_32_32_score0", szKeyName, n, 5, ".dat", piplscore0_1);
	//		SaveMATForDebug(szPathName0, "mask_32_32_score1", szKeyName, n, 5, ".dat", piplscore1_1);
	//		SaveMATForDebug(szPathName0, "mask_16_16_score0", szKeyName, n, 5, ".dat", piplscore0_2);
	//		SaveMATForDebug(szPathName0, "mask_16_16_score1", szKeyName, n, 5, ".dat", piplscore1_2);
	//		SaveMATForDebug(szPathName0, "mask_4_4_score0", szKeyName, n, 5, ".dat", piplscore0_3);
	//		SaveMATForDebug(szPathName0, "mask_4_4_score1", szKeyName, n, 5, ".dat", piplscore1_3);



	//		//AutoLabelPart(szImgNamen, vGMMs, refmask0, 32);
	//		//AutoLabelPart(szImgNamen, vGMMs, refmask, 16);
	//		//AutoLabelPart(szImgNamen, vGMMs, refmask2, 4);
	//		//=============================================================================//

	//		//=============================================================================//
	//		// Then make all 3 labels consistent, that is,
	//		// If 4x4 block some blocks are detected as skin, 
	//		// but corresponding and adjacent 16x16 block 
	//		// didn't have skin label, then discard it.
	//		//=============================================================================//
	//		iplrefmask0 = refmask0.operator IplImage();
	//		iplrefmask = refmask.operator IplImage();
	//		iplrefmask2 = refmask2.operator IplImage();
	//		for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
	//		{
	//			for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
	//			{
	//				db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
	//				nn1 = (int)floor((double)kk1 / 2.0);
	//				nn2 = (int)floor((double)kk2 / 2.0);
	//				cv::Point pt0 = cv::Point(nn1, nn2);
	//				cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
	//				cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
	//				cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
	//				cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

	//				if(cvGetReal2D(&iplrefmask, kk1, kk2))
	//				{
	//					if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
	//						cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
	//						cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
	//						cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
	//						cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
	//					{
	//					}
	//					else
	//					{
	//						cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
	//					}
	//				}
	//			}
	//		}

	//		for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
	//		{
	//			for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
	//			{
	//				db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
	//				nn1 = (int)floor((double)kk1 / 4.0);
	//				nn2 = (int)floor((double)kk2 / 4.0);
	//				cv::Point pt0 = cv::Point(nn1, nn2);
	//				cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
	//				cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
	//				cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
	//				cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

	//				if(cvGetReal2D(&iplrefmask2, kk1, kk2))
	//				{
	//					if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
	//						cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
	//						cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
	//						cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
	//						cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
	//					{
	//					}
	//					else
	//					{
	//						cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
	//					}
	//				}
	//			}
	//		}
	//		//=============================================================================//
	//	
	//		pcurr = vpImgs[n];

	//		pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
	//		pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
	//		//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
	//		pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
	//		pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
	//		//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
	//		pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

	//		cvConvert(&iplrefmask2,pcurr2b);
	//		//curr2b = cv::Mat(cv::abs(pcurr2b));
	//		//curr2b.convertTo(curr2b2, CV_8UC1);

	//		//=============================================================================//
	//		// Find ROI
	//		//=============================================================================//
	//		//storage = cvCreateMemStorage(0);
	//		//cvConvert(pcurr2b, pcurr2bcpy);
	//		//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
	//		//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
	//		//vRect.clear();
	//		//rectROI = CvRect();
	//		////=============================================================================//
	//		//// Compute the bounding rectangle so that ROI can be assured	
	//		////=============================================================================//
	//		//if(first_contour == 0) // No contour found, no hand found
	//		//{
	//		//	cvReleaseMemStorage(&storage);
	//		//	cvReleaseImage(&pcurr2bcpy);
	//		//}
	//		//else
	//		//{
	//		//	rect1 = cvBoundingRect(first_contour);
	//		//	vRect.push_back(rect1);
	//		//	contour = first_contour->h_next;			
	//		//	while(contour != NULL)
	//		//	{
	//		//		recti = cvBoundingRect(contour);
	//		//		vRect.push_back(recti);
	//		//		contour = contour->h_next;
	//		//	}
	//		//	// Now, try to merge these rectangles
	//		//	nrectnum = (int)vRect.size();
	//		//	// Aggolomerate clustering here:
	//		//	vAggRect.clear();
	//		//	nROIIndex = 0;
	//		//	AggMerge(vRect, vAggRect, nROIIndex);
	//		//	// Take he region with largest area
	//		//	rectROI = vAggRect[nROIIndex];
	//		//}
	//		//=============================================================================//

	//;

	//		//=============================================================================//
	//		// Now, combine the mask result with ROI
	//		//=============================================================================//
	//		cvConvert(pcurr2b, pcurr2b_tempb);
	//		//cvAnd(pcurr, pcurr2b_tempb, pMask1);
	//		cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
	//		IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
	//		IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
	//		cvSet(pMask1, CV_RGB(64, 128, 64));
	//		//cv::Mat Mask1(pMask1);
	//		cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
	//		//cv::Mat Tmp(pTmp);
	//		IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
	//		IplImage* piplcurrtemp = &iplcurrtemp;
	//		cvCopy(piplcurrtemp, pMask1, pTmp);
	//		//cvCurrImgTemp.copyTo(Mask1, Tmp);
	//		SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
	//		cvReleaseImage(&pTmp);
	//		cvReleaseImage(&pMask1);
	//		//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
	//		//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
	//		cvCopy(pcurr2b, pcurr2b_tempb3);
	//		//=============================================================================//

	//		//=============================================================================//
	//		// Morphological closing/ hole filling
	//		//=============================================================================//
	//		CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
	//		cvConvert(pcurr2c, pImgb);
	//		IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
	//		cvAnd(pcurr, pImgb, pMask2);
	//		SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
	//		cvReleaseImage(&pMask2);
	//		cvReleaseImage(&pcurr2b_tempb);
	//		//cvReleaseImage(&pcurr2b_tempb2);
	//		cvReleaseImage(&pcurr2b_tempb3);
	//		//=============================================================================//
	//		
	//		//=============================================================================//
	//		// Now, try to detect center of mass of the hand region along with
	//		// fingertip.
	//		// However, currently, finger detection is not good in in car dataset,
	//		// we only use the center of mass, so this part should be simplified in the 
	//		// future if only COM is used
	//		//=============================================================================//
	//		hrTrackFinger = S_OK;
	//		if(vvFeature.size() == 0)
	//		{
	//			hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
	//		}
	//		else
	//		{
	//			hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
	//		}

	//		if( S_OK == hrTrackFinger)
	//		{
	//			size_t nWait = vHandVanishIndices.size();
	//			if(nWait > 0)
	//			{
	//				// simply discard those you cannot save
	//				if(0 < vvFeature.size())
	//				{
	//					cxn_k = vvFeature[vvFeature.size()-1][0];
	//					cyn_k = vvFeature[vvFeature.size()-1][1];
	//					dxn_k = vvFeature[vvFeature.size()-1][2];
	//					dyn_k = vvFeature[vvFeature.size()-1][3];
	//					cxn = vTmp[0];
	//					cyn = vTmp[1];
	//					dxn = vTmp[2];
	//					dyn = vTmp[3];					
	//					for(nk = 1; nk <= (int)nWait; nk++)
	//					{
	//						r = (double)nk / ((double)nWait + 1.0);
	//						vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
	//						vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
	//						vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
	//						vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
	//						vvFeature.push_back(vTmpTmp);
	//					}
	//				}
	//				nWait = 0;
	//				vHandVanishIndices.clear();
	//			}
	//			vvFeature.push_back(vTmp);
	//		}
	//		else
	//		{
	//			vHandVanishIndices.push_back(n);
	//			if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
	//			{	
	//				// Means the gesture ends in online FE				
	//			}
	//			// When it comes to the frame before the last frame,
	//			// we just assume the point stay at the last observation 
	//			// all the time
	//			//!!!! Can be wrong!!!!
	//			else if(n == nStop - 1)
	//			{
	//				for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
	//				{
	//					if(vTmp.size() == 4)
	//					{
	//						vvFeature.push_back(vTmp);
	//					}
	//				}
	//				vHandVanishIndices.clear();
	//			}
	//			else if(n == nStop)
	//			{
	//				// Do nothing, coz nothing is observed 
	//				vHandVanishIndices.clear();
	//			}
	//			cvCurrImgTemp.release();
	//		}
	//		refmask0.release();
	//		refmask.release();
	//		refmask2.release();
	//		//cvReleaseImage(&pcurr);
	//		cvReleaseImage(&pcurr2c);
	//		cvReleaseImage(&pcurr2b);
	//		//cvReleaseImage(&pcurr2bcpy);
	//		cvReleaseImage(&pImgb);
	//		//cvReleaseMemStorage(&storage);
	//		//cvReleaseImage(&pcurr2bcpy);

	//	}
	//	//for(size_t nn = 0; nn < vpImgs.size(); nn++)
	//	//{
	//	//	IplImage* pimg = vpImgs[nn];
	//	//	cvReleaseImage(&pimg);
	//	//	pimg = NULL;
	//	//}
	//	cvReleaseImage(&pBkGnd);
	//	// Write features
	//	if(0 < vvFeature.size())
	//	{
	//		HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
	//		if(hrlast == S_OK)
	//		{
	//			vFEN[nID].IsCreated = true;
	//			string szOutFileName001, szOutFileName002;
	//			MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
	//			MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
	//			WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
	//			//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
	//			//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
	//			size_t nSizet = vvFeature.size() - 1;
	//			nSizet *= 666666;
	//			htVN[szKeyName][0].nStart = 0;
	//			htVN[szKeyName][0].nStop = (int)nSizet;
	//			htVN[szKeyName][0].szName = szLabelName;
	//		}
	//		else if(hrlast == ERROR_INVALID_DATA)
	//		{
	//			vFEN[nID].IsCreated = false;
	//			ofstream feout;
	//			feout.open(szErrLog.c_str(), ios_base::app);
	//			feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
	//			cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
	//			feout.close();
	//		}
	//		else
	//		{
	//			vFEN[nID].IsCreated = false;
	//			ofstream feout;
	//			feout.open(szErrLog.c_str(), ios_base::app);
	//			feout<<szOutFeatureFile<<" is not created for some other reason!\n";
	//			cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
	//			feout.close();
	//		}
	//	}
	//	else
	//	{
	//		vFEN[nID].IsCreated = false;
	//		ofstream feout;
	//		feout.open(szErrLog.c_str(), ios_base::app);
	//		feout<<szOutFeatureFile<<" is not created!\n";
	//		cerr<<szOutFeatureFile<<" is not created!\n";
	//		feout.close();
	//	}
	//	return hr;
	//}

	HRESULT OfflineFEImgSeq2_2_3_part_BS0_Fusion(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		cout<<"Begin Processing Feature file "<<szKeyName<<endl;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				IplImage* pBkGnd_tmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				IplImage* pBkGndTmp_tmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				// ============================================================================================================/
				// For color normalization, brute force!, or pBkGnd = pBkGnd_tmp, pBkGndTmp_tmp = pBkGndTmp
				// ============================================================================================================/
				NormalizeImageRGB(pBkGnd_tmp, pBkGnd);
				NormalizeImageRGB(pBkGndTmp_tmp, pBkGndTmp);
				cvReleaseImage(&pBkGnd_tmp);
				cvReleaseImage(&pBkGndTmp_tmp);
				// ============================================================================================================/
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//

		//=============================================================================//
		// MAP for frames
		//=============================================================================//
		vector<string> vszImagePathNames;
		vector<cv::Mat> vMasks1;
		vector<cv::Mat> vMasks2;
		vector<cv::Mat> vMasks3;
		vector<cv::Mat> vMasks01;
		for(n = 0/*nStopBkGnd + 1*/; n < vpImgs.size()/*= nStop*/; n++)
		{
			vszImagePathNames.push_back(ImgPathNames[n]);
		}
		if(m_nHandModelIndex == 0)
			m_nHandModelIndex = 1;
		//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion(vszImagePathNames, vGMMs, vMasks1, 32, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1);
		//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion(vszImagePathNames, vGMMs, vMasks2, 16, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1);
		//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion(vszImagePathNames, vGMMs, vMasks3, 4, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1);
		vector<cv::Mat> vScore0_1;
		vector<cv::Mat> vScore1_1;
		vector<cv::Mat> vScore0_2;
		vector<cv::Mat> vScore1_2;
		vector<cv::Mat> vScore0_3;
		vector<cv::Mat> vScore1_3;
		
		vector<cv::Mat> vLLR_1;
		vector<cv::Mat> vCorr_1;
		vector<cv::Mat> vKS_1;
		vector<cv::Mat> vLLR_2;
		vector<cv::Mat> vCorr_2;
		vector<cv::Mat> vKS_2;
		vector<cv::Mat> vLLR_3;
		vector<cv::Mat> vCorr_3;
		vector<cv::Mat> vKS_3;

		//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks1, 32, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_1, vScore1_1, vLLR_1, vCorr_1);
		//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks2, 16, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_2, vScore1_2, vLLR_2, vCorr_2);
		//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks3, 4, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_3, vScore1_3, vLLR_3, vCorr_3);
		GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks1, 32, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_1, vScore1_1, vLLR_1, vCorr_1, vKS_1);
		GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks2, 16, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_2, vScore1_2, vLLR_2, vCorr_2, vKS_2);
		GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores(vszImagePathNames, vGMMs, vMasks3, 4, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_3, vScore1_3, vLLR_3, vCorr_3, vKS_3);

		//IplImage tmp1 = vMasks1[16].operator IplImage();
		//IplImage tmp01 = vMasks01[16].operator IplImage();

		//SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, 16, 5, "_o.jpg", &(tmp1));
		//SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, 16, 5, "_ws.jpg", &(tmp01));

		
		
		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			//cout<<"Doing frame "<<n<<endl;
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask2, 4);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask0, 32, BkGnd);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask, 16, BkGnd);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask2, 4, BkGnd);
			refmask0 = vMasks1[n];// - nStopBkGnd - 1];
			refmask = vMasks2[n];// - nStopBkGnd - 1];
			refmask2 = vMasks3[n];// - nStopBkGnd - 1];
			IplImage iplrefmask0 = refmask0.operator IplImage();
			IplImage iplrefmask1 = refmask.operator IplImage();
			IplImage iplrefmask2 = refmask2.operator IplImage();
			IplImage* prefmask0 = &iplrefmask0;
			IplImage* prefmask1 = &iplrefmask1;
			IplImage* prefmask2 = &iplrefmask2;
			SaveMATForDebug(szPathName0, "mask_32_32", szKeyName, n, 5, ".dat", prefmask0);
			SaveMATForDebug(szPathName0, "mask_16_16", szKeyName, n, 5, ".dat", prefmask1);
			SaveMATForDebug(szPathName0, "mask_4_4", szKeyName, n, 5, ".dat", prefmask2);
			// Save scores from LLR and Corr
			IplImage iplsLLR1 = vLLR_1[n].operator IplImage();
			IplImage* piplsLLR1 = &iplsLLR1;
			SaveMATForDebug(szPathName0, "sLLR_32_32", szKeyName, n, 5, ".dat", piplsLLR1);

			IplImage iplsLLR2 = vLLR_2[n].operator IplImage();
			IplImage* piplsLLR2 = &iplsLLR2;
			SaveMATForDebug(szPathName0, "sLLR_16_16", szKeyName, n, 5, ".dat", piplsLLR2);

			IplImage iplsLLR3 = vLLR_3[n].operator IplImage();
			IplImage* piplsLLR3 = &iplsLLR3;
			SaveMATForDebug(szPathName0, "sLLR_4_4", szKeyName, n, 5, ".dat", piplsLLR3);

			IplImage iplsCorr1 = vCorr_1[n].operator IplImage();
			IplImage* piplsCorr1 = &iplsCorr1;
			SaveMATForDebug(szPathName0, "Corr_32_32", szKeyName, n, 5, ".dat", piplsCorr1);

			IplImage iplsCorr2 = vCorr_2[n].operator IplImage();
			IplImage* piplsCorr2 = &iplsCorr2;
			SaveMATForDebug(szPathName0, "Corr_16_16", szKeyName, n, 5, ".dat", piplsCorr2);

			IplImage iplsCorr3 = vCorr_3[n].operator IplImage();
			IplImage* piplsCorr3 = &iplsCorr3;
			SaveMATForDebug(szPathName0, "Corr_4_4", szKeyName, n, 5, ".dat", piplsCorr3);

			IplImage iplsKS1 = vKS_1[n].operator IplImage();
			IplImage* piplsKS1 = &iplsKS1;
			SaveMATForDebug(szPathName0, "KS_32_32", szKeyName, n, 5, ".dat", piplsKS1);

			IplImage iplsKS2 = vKS_2[n].operator IplImage();
			IplImage* piplsKS2 = &iplsKS2;
			SaveMATForDebug(szPathName0, "KS_16_16", szKeyName, n, 5, ".dat", piplsKS2);

			IplImage iplsKS3 = vKS_3[n].operator IplImage();
			IplImage* piplsKS3 = &iplsKS3;
			SaveMATForDebug(szPathName0, "KS_4_4", szKeyName, n, 5, ".dat", piplsKS3);

			IplImage iplscore0_1 = vScore0_1[n].operator IplImage();
			IplImage iplscore1_1 = vScore1_1[n].operator IplImage();
			IplImage iplscore0_2 = vScore0_2[n].operator IplImage();
			IplImage iplscore1_2 = vScore1_2[n].operator IplImage();
			IplImage iplscore0_3 = vScore0_3[n].operator IplImage();
			IplImage iplscore1_3 = vScore1_3[n].operator IplImage();
			IplImage* piplscore0_1 = &iplscore0_1;
			IplImage* piplscore1_1 = &iplscore1_1;
			IplImage* piplscore0_2 = &iplscore0_2;
			IplImage* piplscore1_2 = &iplscore1_2;
			IplImage* piplscore0_3 = &iplscore0_3;
			IplImage* piplscore1_3 = &iplscore1_3;
			SaveMATForDebug(szPathName0, "mask_32_32_score0", szKeyName, n, 5, ".dat", piplscore0_1);
			SaveMATForDebug(szPathName0, "mask_32_32_score1", szKeyName, n, 5, ".dat", piplscore1_1);
			SaveMATForDebug(szPathName0, "mask_16_16_score0", szKeyName, n, 5, ".dat", piplscore0_2);
			SaveMATForDebug(szPathName0, "mask_16_16_score1", szKeyName, n, 5, ".dat", piplscore1_2);
			SaveMATForDebug(szPathName0, "mask_4_4_score0", szKeyName, n, 5, ".dat", piplscore0_3);
			SaveMATForDebug(szPathName0, "mask_4_4_score1", szKeyName, n, 5, ".dat", piplscore1_3);



			//AutoLabelPart(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPart(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPart(szImgNamen, vGMMs, refmask2, 4);
			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}
			vector<double> comvec;
			// ===============================================================================//
			// Output combined mask for durther processing
			// ===============================================================================//
			SaveMATForDebug(szPathName0, "mask_combine", szKeyName, n, 5, ".dat", prefmask2);
			// ===============================================================================//
			COMPoints(refmask2, comvec);
			vvFeature.push_back(comvec);
			//=============================================================================//
		
			//pcurr = vpImgs[n];

			//pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			////pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			////pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			//cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			//cvConvert(pcurr2b, pcurr2b_tempb);
			////cvAnd(pcurr, pcurr2b_tempb, pMask1);
			//cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			//IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			//IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			//cvSet(pMask1, CV_RGB(64, 128, 64));
			////cv::Mat Mask1(pMask1);
			//cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			////cv::Mat Tmp(pTmp);
			//IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			//IplImage* piplcurrtemp = &iplcurrtemp;
			//cvCopy(piplcurrtemp, pMask1, pTmp);
			////cvCurrImgTemp.copyTo(Mask1, Tmp);
			//SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			//cvReleaseImage(&pTmp);
			//cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			//cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			//CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			//cvConvert(pcurr2c, pImgb);
			//IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//cvAnd(pcurr, pImgb, pMask2);
			//SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			//cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			//hrTrackFinger = S_OK;
			//if(vvFeature.size() == 0)
			//{
			//	hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			//}
			//else
			//{
			//	hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			//}

			//if( S_OK == hrTrackFinger)
			//{
			//	size_t nWait = vHandVanishIndices.size();
			//	if(nWait > 0)
			//	{
			//		// simply discard those you cannot save
			//		if(0 < vvFeature.size())
			//		{
			//			cxn_k = vvFeature[vvFeature.size()-1][0];
			//			cyn_k = vvFeature[vvFeature.size()-1][1];
			//			dxn_k = vvFeature[vvFeature.size()-1][2];
			//			dyn_k = vvFeature[vvFeature.size()-1][3];
			//			cxn = vTmp[0];
			//			cyn = vTmp[1];
			//			dxn = vTmp[2];
			//			dyn = vTmp[3];					
			//			for(nk = 1; nk <= (int)nWait; nk++)
			//			{
			//				r = (double)nk / ((double)nWait + 1.0);
			//				vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
			//				vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
			//				vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
			//				vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
			//				vvFeature.push_back(vTmpTmp);
			//			}
			//		}
			//		nWait = 0;
			//		vHandVanishIndices.clear();
			//	}
			//	vvFeature.push_back(vTmp);
			//}
			//else
			//{
			//	vHandVanishIndices.push_back(n);
			//	if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
			//	{	
			//		// Means the gesture ends in online FE				
			//	}
			//	// When it comes to the frame before the last frame,
			//	// we just assume the point stay at the last observation 
			//	// all the time
			//	//!!!! Can be wrong!!!!
			//	else if(n == nStop - 1)
			//	{
			//		for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
			//		{
			//			if(vTmp.size() == 4)
			//			{
			//				vvFeature.push_back(vTmp);
			//			}
			//		}
			//		vHandVanishIndices.clear();
			//	}
			//	else if(n == nStop)
			//	{
			//		// Do nothing, coz nothing is observed 
			//		vHandVanishIndices.clear();
			//	}
			//	cvCurrImgTemp.release();
			//}
			refmask0.release();
			refmask.release();
			refmask2.release();
			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		cvReleaseImage(&pBkGnd);

		if(0 < vvFeature.size())
		{
			vector< vector<double> > vvFeature2 = vvFeature;
			vector<int> lastcom(vvFeature.size(), 0);
			vector<int> nextcom(vvFeature.size(), 0);
			int lastcomtmp = -1;
			int nextcomtmp = -1;
			for(int nn1 = 0; nn1 < (int)vvFeature.size(); nn1++)
			{
				if(vvFeature[nn1][0] > 0 && vvFeature[nn1][1] > 0)
					lastcomtmp = nn1;
				lastcom[nn1] = lastcomtmp;
			}

			for(int nn1 = (int)vvFeature.size() - 1; nn1 >=0 ; nn1--)
			{
				if(vvFeature[nn1][0] > 0 && vvFeature[nn1][1] > 0)
					nextcomtmp = (int)nn1;
				nextcom[nn1] = nextcomtmp;
			}
			// interpolate COM
			for(int nn1 = 0; nn1 < (int)vvFeature.size(); nn1++)
			{
				if(lastcom[nn1] >= 0 && nextcom[nn1] >= 0 && nextcom[nn1] != lastcom[nn1])
				{
					double distcomtmpd = (double)(nextcom[nn1] - lastcom[nn1]);
					double distcomtmpn = (double)(nn1 - lastcom[nn1]);
					vector<double> distcomtmpv(2, 0.0);
					distcomtmpv[0] = vvFeature[nextcom[nn1]][0] - vvFeature[lastcom[nn1]][0];
					distcomtmpv[1] = vvFeature[nextcom[nn1]][1] - vvFeature[lastcom[nn1]][1];
					vvFeature2[nn1][0] = vvFeature[nextcom[nn1]][0] + distcomtmpv[0]  * distcomtmpn / distcomtmpd;	
					vvFeature2[nn1][1] = vvFeature[nextcom[nn1]][1] + distcomtmpv[1]  * distcomtmpn / distcomtmpd;	
				}
				else if(lastcom[nn1] < 0 && nextcom[nn1] >= 0)
				{
					vvFeature2[nn1][0] = vvFeature[nextcom[nn1]][0];
					vvFeature2[nn1][1] = vvFeature[nextcom[nn1]][1];
				}
				else if(lastcom[nn1] >= 0 && nextcom[nn1] < 0)
				{
					vvFeature2[nn1][0] = vvFeature[lastcom[nn1]][0];
					vvFeature2[nn1][1] = vvFeature[lastcom[nn1]][1];
				}
			}
			//cout<<"ready to write HTK feature\n";
			HRESULT hrlast = WriteFeatureAsHTKMFCC(vvFeature2, szOutFeatureFile_uhalf);
			//cout<<"HTK feature written\n";
			if(hrlast == S_OK)
			{
				//cout<<"before vFEN[nID].IsCreated = true;"<<endl;
				vFEN[nID].IsCreated = true;
				//cout<<"before string szOutFileName001, szOutFileName002;"<<endl;
				string szOutFileName001, szOutFileName002;
				//cout<<"before MakeDebugOutputName(szPathName0, \"COM\", szKeyName, \".jpg\", szOutFileName001);"<<endl;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				//cout<<"before writing trajectory to jpg\n";
				WriteFeatureAsDEBUGTRAJ(vvFeature2, szOutFileName001);
				//cout<<"after writing trajectory to jpg\n";
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		//cvReleaseImage(&pBkGnd);
		// Write features


		//if(0 < vvFeature.size())
		//{
		//	HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
		//	if(hrlast == S_OK)
		//	{
		//		vFEN[nID].IsCreated = true;
		//		string szOutFileName001, szOutFileName002;
		//		MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
		//		MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
		//		WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
		//		//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
		//		//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
		//		size_t nSizet = vvFeature.size() - 1;
		//		nSizet *= 666666;
		//		htVN[szKeyName][0].nStart = 0;
		//		htVN[szKeyName][0].nStop = (int)nSizet;
		//		htVN[szKeyName][0].szName = szLabelName;
		//	}
		//	else if(hrlast == ERROR_INVALID_DATA)
		//	{
		//		vFEN[nID].IsCreated = false;
		//		ofstream feout;
		//		feout.open(szErrLog.c_str(), ios_base::app);
		//		feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
		//		cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
		//		feout.close();
		//	}
		//	else
		//	{
		//		vFEN[nID].IsCreated = false;
		//		ofstream feout;
		//		feout.open(szErrLog.c_str(), ios_base::app);
		//		feout<<szOutFeatureFile<<" is not created for some other reason!\n";
		//		cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
		//		feout.close();
		//	}
		//}
		//else
		//{
		//	vFEN[nID].IsCreated = false;
		//	ofstream feout;
		//	feout.open(szErrLog.c_str(), ios_base::app);
		//	feout<<szOutFeatureFile<<" is not created!\n";
		//	cerr<<szOutFeatureFile<<" is not created!\n";
		//	feout.close();
		//}
		//cout<<"Finish properly\n";
		return hr;
	}


	HRESULT OfflineFEImgSeq2_2_3_part_BS0_SMAPPerFrame(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//

		//=============================================================================//
		// SMAP for frames
		//=============================================================================//
		vector<string> vszImagePathNames;
		vector<cv::Mat> vMasks1;
		vector<cv::Mat> vMasks2;
		vector<cv::Mat> vMasks3;
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			vszImagePathNames.push_back(ImgPathNames[n]);
		}
		if(m_nHandModelIndex == 0)
			m_nHandModelIndex = 1;
		GMMAdaptTool::AutoLabelSMAPPerFrame_BS0(vszImagePathNames, vGMMs, vMasks1, 32, m_nHandModelIndex, m_szGMMPath, BkGnd);
		GMMAdaptTool::AutoLabelSMAPPerFrame_BS0(vszImagePathNames, vGMMs, vMasks2, 16, m_nHandModelIndex, m_szGMMPath, BkGnd);
		GMMAdaptTool::AutoLabelSMAPPerFrame_BS0(vszImagePathNames, vGMMs, vMasks3, 4, m_nHandModelIndex, m_szGMMPath, BkGnd);
		
		
		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask2, 4);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask0, 32, BkGnd);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask, 16, BkGnd);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask2, 4, BkGnd);
			refmask0 = vMasks1[n - nStopBkGnd - 1];
			refmask = vMasks2[n - nStopBkGnd - 1];
			refmask2 = vMasks3[n - nStopBkGnd - 1];

			//AutoLabelPart(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPart(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPart(szImgNamen, vGMMs, refmask2, 4);
			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}
			//=============================================================================//
		
			pcurr = vpImgs[n];

			pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			//cv::Mat Mask1(pMask1);
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			//cv::Mat Tmp(pTmp);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			cvAnd(pcurr, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();
			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}
	
	HRESULT OfflineFEImgSeq2_2_3_part_BS0_SMAPASeq(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//

		//=============================================================================//
		// SMAP for frames
		//=============================================================================//
		vector<string> vszImagePathNames;
		vector<cv::Mat> vMasks1;
		vector<cv::Mat> vMasks2;
		vector<cv::Mat> vMasks3;
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			vszImagePathNames.push_back(ImgPathNames[n]);
		}
		if(m_nHandModelIndex == 0)
			m_nHandModelIndex = 1;
		GMMAdaptTool::AutoLabelSMAPASeq_BS0(vszImagePathNames, vGMMs, vMasks1, 32, m_nHandModelIndex, m_szGMMPath, BkGnd);
		GMMAdaptTool::AutoLabelSMAPASeq_BS0(vszImagePathNames, vGMMs, vMasks2, 16, m_nHandModelIndex, m_szGMMPath, BkGnd);
		GMMAdaptTool::AutoLabelSMAPASeq_BS0(vszImagePathNames, vGMMs, vMasks3, 4, m_nHandModelIndex, m_szGMMPath, BkGnd);
		
		
		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask2, 4);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask0, 32, BkGnd);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask, 16, BkGnd);
			//AutoLabel_BS0(szImgNamen, vGMMs, refmask2, 4, BkGnd);
			refmask0 = vMasks1[n - nStopBkGnd - 1];
			refmask = vMasks2[n - nStopBkGnd - 1];
			refmask2 = vMasks3[n - nStopBkGnd - 1];

			//AutoLabelPart(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPart(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPart(szImgNamen, vGMMs, refmask2, 4);
			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}
			//=============================================================================//
		
			pcurr = vpImgs[n];

			pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			//cv::Mat Mask1(pMask1);
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			//cv::Mat Tmp(pTmp);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			cvAnd(pcurr, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();
			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}
	

	HRESULT OfflineFEImgSeq2_2_3_part_BS0_Rough(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector< GMM >& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//
		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			AutoLabel_BS0(szImgNamen, vGMMs, refmask0, 32, BkGnd);
			AutoLabel_BS0(szImgNamen, vGMMs, refmask, 16, BkGnd);
			//AutoLabel_BS0(szImgNamen, vvvDecisionTable, refmask2, 4, BkGnd);

			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			//iplrefmask2 = refmask2.operator IplImage();
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}
			cv::resize(refmask, refmask2, cv::Size(refmask.cols * 4, refmask.rows * 4), 0, 0, INTER_NEAREST); 
			iplrefmask2 = refmask2.operator IplImage();
			//=============================================================================//
		
			pcurr = vpImgs[n];

			pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			//cv::Mat Mask1(pMask1);
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			//cv::Mat Tmp(pTmp);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			cvAnd(pcurr, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();
			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}

	HRESULT OfflineFEImgSeq2_2_3_part_BS0_Table(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, bool*** vvvDecisionTable)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//
		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//



#ifdef _MSC_VER
			DWORD dw1 = ::timeGetTime();
#endif
			Approx_AutoLabel_BS0(szImgNamen, vvvDecisionTable, refmask0, 32, BkGnd);
#ifdef _MSC_VER
			DWORD dw2 = ::timeGetTime();
#endif
			Approx_AutoLabel_BS0(szImgNamen, vvvDecisionTable, refmask, 16, BkGnd);
#ifdef _MSC_VER
			DWORD dw3 = ::timeGetTime();
#endif
			Approx_AutoLabel_BS0(szImgNamen, vvvDecisionTable, refmask2, 4, BkGnd);
#ifdef _MSC_VER
			DWORD dw4 = ::timeGetTime();
			cout<<"Scale 32*32 takes "<<dw2-dw1<<" ms"<<endl;
			cout<<"Scale 16*16 takes "<<dw3-dw2<<" ms"<<endl;
			cout<<"Scale 4*4 takes "<<dw4-dw3<<" ms"<<endl;
#endif
			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}
			//=============================================================================//
		
			pcurr = vpImgs[n];

			pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			//cv::Mat Mask1(pMask1);
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			//cv::Mat Tmp(pTmp);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			cvAnd(pcurr, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();
			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}


	HRESULT OfflineFEImgSeq2_2_3_with_MCD(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		//for(n = nStart; n <= nStopBkGnd; n++)
		//{
		//	pcurr = vpImgs[n];
		//	if(n == nStart)
		//	{
		//		pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//		pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//		cvSet(pBkGnd, cvScalar(0.0));
		//		cvSet(pBkGndTmp, cvScalar(0.0));
		//	}
		//	pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
		//	cvConvert(pcurr, pcurr2);
		//	GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
		//	cvReleaseImage(&pcurr);
		//	cvReleaseImage(&pcurr2);
		//}
		//AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		//cvReleaseImage(&pBkGndTmp);
		//=============================================================================//

		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPartSM(szImgNamen, vGMMs, refmask2, 4);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask0, 32);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask, 16);
			AutoLabelRGBOnly(szImgNamen, vGMMs, refmask2, 4);

			//AutoLabelPart(szImgNamen, vGMMs, refmask0, 32);
			//AutoLabelPart(szImgNamen, vGMMs, refmask, 16);
			//AutoLabelPart(szImgNamen, vGMMs, refmask2, 4);
			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}
			//=============================================================================//
		
			pcurr = vpImgs[n];

			pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);

			cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			//cv::Mat Mask1(pMask1);
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			//cv::Mat Tmp(pTmp);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			cvAnd(pcurr, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();

			//cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//for(size_t nn = 0; nn < vpImgs.size(); nn++)
		//{
		//	IplImage* pimg = vpImgs[nn];
		//	cvReleaseImage(&pimg);
		//	pimg = NULL;
		//}
		//cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			vector< vector<double> > vvAdjustedFeature;
			CheckTrajValidityBYU2D(vvFeature, vvAdjustedFeature);
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvAdjustedFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
MYRETURN:
		return hr;
	}



	HRESULT OfflineFEImgSeq2_2_3_part_cheat(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs, LoadJPGLabelInfo& lab)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		////IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		for(size_t nn = 0; nn < ImgPathNames.size(); nn++)
		{
			if(lab.m_hJPGLavelInfo.find(ImgPathNames[nn]) == lab.m_hJPGLavelInfo.end())
			{
			}
			else
			{
				if(lab.m_hJPGLavelInfo[ImgPathNames[nn]].m_nSkinCount)
				{
					nStopBkGnd = nn - 1;
					break;
				}
			}
		}
		for(size_t nn = ImgPathNames.size() - 1; nn >= 0; nn--)
		{
			if(lab.m_hJPGLavelInfo.find(ImgPathNames[nn]) == lab.m_hJPGLavelInfo.end())
			{
			}
			else
			{
				if(lab.m_hJPGLavelInfo[ImgPathNames[nn]].m_nSkinCount)
				{
					nStop = nn;
					break;
				}
			}
		}
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = vpImgs[n];
			cvReleaseImage(&pcurr);
			vpImgs[n] = NULL;
		}
		for(n = nStop + 1; n < vpImgs.size(); n++)
		{
			pcurr = vpImgs[n];
			cvReleaseImage(&pcurr);
			vpImgs[n] = NULL;
		}
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			// Load label:
			int nBW = lab.m_hJPGLavelInfo[szImgNamen].m_nBlockWidth;
			int nBH = lab.m_hJPGLavelInfo[szImgNamen].m_nBlockHeight;
			cv::Mat refmask00(nBH, nBW, CV_8UC1);
			refmask00.setTo(cv::Scalar(0,0,0,0));
			for(int kk = 0; kk < nBH; kk++)
			{
				for(int ll = 0; ll < nBW; ll++)
				{
					if(lab.m_hJPGLavelInfo[szImgNamen].m_ppnLabel[kk][ll] == 1)
					{
						refmask00.at<uchar>(kk, ll) = 255;
					}
				}
			}

			cv::Mat refmask, refmaskq, refmask2, refmask0;
			cv::resize(refmask00, refmask2, cvSize(nBW*4, nBH*4), 0, 0, cv::INTER_NEAREST);
		
			pcurr = vpImgs[n];

			pImgb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			//pcurr2bcpy = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2c = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			//pcurr2b_tempb2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
			IplImage iplrefmask2 = refmask2.operator IplImage();
			cvConvert(&iplrefmask2,pcurr2b);
			//curr2b = cv::Mat(cv::abs(pcurr2b));
			//curr2b.convertTo(curr2b2, CV_8UC1);

			//=============================================================================//
			// Find ROI
			//=============================================================================//
			//storage = cvCreateMemStorage(0);
			//cvConvert(pcurr2b, pcurr2bcpy);
			//cvFindContours( pcurr2bcpy, storage, &first_contour, sizeof(CvContour), 
			//				CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0));
			//vRect.clear();
			//rectROI = CvRect();
			////=============================================================================//
			//// Compute the bounding rectangle so that ROI can be assured	
			////=============================================================================//
			//if(first_contour == 0) // No contour found, no hand found
			//{
			//	cvReleaseMemStorage(&storage);
			//	cvReleaseImage(&pcurr2bcpy);
			//}
			//else
			//{
			//	rect1 = cvBoundingRect(first_contour);
			//	vRect.push_back(rect1);
			//	contour = first_contour->h_next;			
			//	while(contour != NULL)
			//	{
			//		recti = cvBoundingRect(contour);
			//		vRect.push_back(recti);
			//		contour = contour->h_next;
			//	}
			//	// Now, try to merge these rectangles
			//	nrectnum = (int)vRect.size();
			//	// Aggolomerate clustering here:
			//	vAggRect.clear();
			//	nROIIndex = 0;
			//	AggMerge(vRect, vAggRect, nROIIndex);
			//	// Take he region with largest area
			//	rectROI = vAggRect[nROIIndex];
			//}
			//=============================================================================//

	;

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			//cv::Mat Mask1(pMask1);
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			//cv::Mat Tmp(pTmp);
			IplImage iplcurrimgtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrimgtemp = &iplcurrimgtemp;
			cvCopy(piplcurrimgtemp, pMask1, pTmp);
			//cvCurrImgTemp.copyTo(Mask1, Tmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//cvAnd(pcurr2b_tempb, &iplrefmask2, pcurr2b_tempb2);
			//cvConvert(pcurr2b_tempb2, pcurr2b_tempb3);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			cvAnd(pcurr, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			//cvReleaseImage(&pcurr2b_tempb);
			//cvReleaseImage(&pcurr2b_tempb2);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
			}
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			//cvReleaseImage(&pcurr2bcpy);
			cvReleaseImage(&pImgb);
			//cvReleaseMemStorage(&storage);
			//cvReleaseImage(&pcurr2bcpy);

		}
		//cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(S_OK == hrlast)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}

	int MkdirHier(string sztmp)
	{
		//////////////////////////////////////////////////////////////////////
		// Make directory hierachically
		//////////////////////////////////////////////////////////////////////
		int nFound1 = 0;
		int nFound2 = 0;
		int MakeResult = 0;
		int IfSuccess = 1;
#ifdef _WIN32
		string szdelim = "\\";\
		char delim = '\\';
#else
		string szdelim = "/";
		char delim = '/';
#endif
		nFound1 = (int)sztmp.find(szdelim, 0);
		nFound2 = (int)sztmp.find(szdelim, nFound1 + 1);
		string szdir = sztmp.substr(0, nFound1);
		string szdmdtmp = szdir; //"";//sztmp.substr(0, nFound1);
		vector<string> vsz1;
		while(nFound2 >= 0)
		{
			vsz1.push_back(sztmp.substr(nFound1 + 1, nFound2-nFound1-1));
			nFound1 = nFound2;
			nFound2 = (int)sztmp.find(szdelim, nFound1 + 1);
		}
		
		if(sztmp[sztmp.length() - 1] != delim)
		{
			vsz1.push_back(sztmp.substr(nFound1 + 1, sztmp.size()-nFound1-1));
		}
		string syscom = "";
		for(size_t nn = 0; nn < vsz1.size(); nn++)
		{
			//cout<<"MkdirHier, index="<<nn<<", the string to be appended is ";
			//cout<<vsz1[nn]<<endl;
			szdmdtmp.append(szdelim);
			szdmdtmp.append(vsz1[nn].c_str());
			syscom.append(szdmdtmp);
			syscom.append(szdelim);
			string s1 = szdir;
			s1.append(":\\");
			struct stat st;
#ifdef _WIN32
			string mkdir_str = "mkdir " + syscom;
			if(!file_exist(syscom))
				system(mkdir_str.c_str());
#else
			string mkdir_str = "mkdir -p " + syscom;
			if(!file_exist(syscom))
				system(mkdir_str.c_str());
#endif
/*
			if(stat(syscom.c_str(),&st) != 0)
			{
#ifdef _WIN32
				MakeResult = _mkdir(syscom.c_str());
#else
				MakeResult = mkdir(syscom.c_str(), 0755);

#endif
			}
			if(MakeResult == -1)
			{
				int nerrno = 0;
#ifdef _WIN32
				::_get_errno(&nerrno);
#else
				nerrno = 2;				
#endif
				if(nerrno == 2)
				{
					cerr<<"Cannot create path while processing command "<<syscom<<endl;
					IfSuccess = -1;
				}

			}
*/
			syscom = "";
		}
		//////////////////////////////////////////////////////////////////////
		return IfSuccess;
	}	

	//=======================================================
	// New functions
	//=======================================================
	void MakeDebugOutputName(string szPathName0, string szSubFolder, string szKeyName, string szExt, string& szOutFileName)
	{
		char delim;
		string szdelim;
	#ifdef _WIN32
		delim = '\\';
		szdelim = "\\";
	#else
		delim = '/';
		szdelim = "/";
	#endif
		szOutFileName = szPathName0;
		//cout<<"In MakeDebugOutputName\n";
		//cout<<"szPathName0 is ["<<szPathName0<<"]\n";
		//cout<<"szSubFolder is ["<<szSubFolder<<"]\n";
		//cout<<"szKeyName is ["<<szKeyName<<"]\n";
		//cout<<"szExt is ["<<szExt<<"]\n";
		//cout<<"szOutFileName is ["<<szOutFileName<<"]\n";
		if(szSubFolder[0] != delim && szPathName0[szPathName0.size() - 1] != delim)
		{
			string szTemp = szdelim;
			szSubFolder = szTemp + szSubFolder;
		}
		if(szSubFolder[szSubFolder.size() - 1] != delim)
			szSubFolder.append(szdelim);
		szOutFileName.append(szSubFolder);
		//cout<<"before mkdirhier\n";
		MkdirHier(szOutFileName);
		//cout<<"after mkdirhier\n";
        #ifdef _WIN32
        #else
		//cout<<"doing chmod\n";
                string szchmod = "chmod 755 " + szOutFileName;
		system(szchmod.c_str());
		//cout<<"chmod done\n";
        #endif

		szOutFileName.append(szKeyName);
		szOutFileName.append(szExt);
		//cout<<"after append szKeyName and szExt\n";
			
	}

	void SavecvImageForDebug(string szPathName0, string szSubFolder, string szKeyName, int n, int nLength, string szExt, IplImage* pImg)
	{
		//SavecvImageForDebug(szPathName0, "Debug01", szKeyName, n, 5, ".jpg", pImg);
		char delim;
		string szdelim;
	#ifdef _WIN32
		delim = '\\';
		szdelim = "\\";
	#else
		delim = '/';
		szdelim = "/";
	#endif
		string st111 = szPathName0;
		if(szSubFolder[0] != delim && szPathName0[szPathName0.size() - 1] != delim)
		{
			string szTemp = szdelim;
			szSubFolder = szTemp + szSubFolder;
		}
		if(szSubFolder[szSubFolder.size() - 1] != delim)
			szSubFolder.append(szdelim);
		st111.append(szSubFolder);
		MkdirHier(st111);
        #ifdef _WIN32
        #else
                string szchmod1 = "chmod 755 " + st111;
		system(szchmod1.c_str());
        #endif
		st111.append(szKeyName);
		char* buf01 = new char[nLength + 1];
		memset(buf01, 0, nLength + 1);
		string szTmpformat = "%0";
		char* buf02 = new char[nLength + 1];
		memset(buf02, 0, nLength + 1);
		sprintf(buf02, "%i", nLength);
		//itoa(nLength, buf02, 10);
		szTmpformat.append(buf02);
		szTmpformat.append("d");
		sprintf(buf01, szTmpformat.c_str(), n);
		st111.append(buf01);
		st111.append(szExt);
		if(pImg != NULL)
		{
			cvSaveImage(st111.c_str(), pImg);
		}
		delete [ ] buf02;
		delete [ ] buf01;
	}

	//=======================================================

	void SaveMATForDebug(string szPathName0, string szSubFolder, string szKeyName, int n, int nLength, string szExt, IplImage* pImg)
	{
		//SavecvImageForDebug(szPathName0, "Debug01", szKeyName, n, 5, ".jpg", pImg);
		char delim;
		string szdelim;
	#ifdef _WIN32
		delim = '\\';
		szdelim = "\\";
	#else
		delim = '/';
		szdelim = "/";
	#endif
		string st111 = szPathName0;
		if(szSubFolder[0] != delim && szPathName0[szPathName0.size() - 1] != delim)
		{
			string szTemp = szdelim;
			szSubFolder = szTemp + szSubFolder;
		}
		if(szSubFolder[szSubFolder.size() - 1] != delim)
			szSubFolder.append(szdelim);
		st111.append(szSubFolder);
		MkdirHier(st111);
        #ifdef _WIN32
        #else
                //string szchmod1 = "chmod 777 " + st111;
		//system(szchmod1.c_str());
        #endif
		st111.append(szKeyName);
		char* buf01 = new char[nLength + 1];
		memset(buf01, 0, nLength + 1);
		string szTmpformat = "%0";
		char* buf02 = new char[nLength + 1];
		memset(buf02, 0, nLength + 1);
		sprintf(buf02, "%i", nLength);
		//itoa(nLength, buf02, 10);
		szTmpformat.append(buf02);
		szTmpformat.append("d");
		sprintf(buf01, szTmpformat.c_str(), n);
		st111.append(buf01);
		st111.append(szExt);
		if(pImg != NULL)
		{
			FILE* f1 = NULL;
			f1 = fopen(st111.c_str(), "wb");
			if(!f1)
			{
				cerr<<"Cannot open file "<<st111<<" for writing"<<endl;
				delete [ ] buf02;
				delete [ ] buf01;
				return;
			}
			int height = pImg->height;
			int width = pImg->width;
			fwrite(&height, sizeof(int), 1, f1);
			fwrite(&width, sizeof(int), 1, f1);
			for(int n = 0; n < pImg->height; ++n)
			{
				for(int m = 0; m < pImg->width; ++m)
				{
					double db1 = cvGetReal2D(pImg, n, m);
					fwrite(&db1, sizeof(double), 1, f1);
				}
			}
			fclose(f1);
		}
		delete [ ] buf02;
		delete [ ] buf01;
	}

	//=======================================================


	// Core of Stochastic matching on feature space
	int GMMBiasAdapt(vector< vector<double> >& vvData, matrix<double>& newmub, vector<GMM>& oldGMMs, vector<GMM>& newGMMs)
	{
		return GMMAdaptTool::GMMBiasAdapt(vvData, newmub, oldGMMs, newGMMs);
	}

	// Get the minimum and maximum value of matrix stored in vector< vector<double> > form at n-th dimension
	HRESULT GetvvRangeAtnD(const vector< vector<double> >& vvData, unsigned int n, double& dbMin, double& dbMax)
	{
		return (HRESULT) GMMAdaptTool::GetvvRangeAtnD(vvData, n, dbMin, dbMax);
	}

	// Calculate the mean vector of a vector< vector<double> > data matrix by trating each vector<double> as a data point
	HRESULT CalcvvMean(const vector< vector<double> >& vvData, vector<double>& vMean)
	{
		return (HRESULT) GMMAdaptTool::CalcvvMean(vvData, vMean);
	}
	
	// normalize so that the sequence are in [0,1]
	HRESULT Makein01(const vector< vector<double> >& vvData, vector< vector<double> >& vvNewData)
	{
		return (HRESULT) GMMAdaptTool::Makein01(vvData, vvNewData);
	}

	// Make vector< vector<double> > data zero mean and range from [-1, 1]
	HRESULT MakeZeroMeanAndUnitRange(const vector< vector<double> >& vvData, vector< vector<double> >& vvNewData)
	{
		return (HRESULT) GMMAdaptTool::MakeZeroMeanAndUnitRange(vvData, vvNewData);
	}	
	
	// Write the feature vector as HTK format MFCC feature file
	HRESULT WriteFeatureAsHTKMFCC(const vector< vector<double> >& vvData, string szOutFileName, short parmKind=2054, int sampPeriod=666666)
	{
		return (HRESULT) GMMAdaptTool::WriteFeatureAsHTKMFCC(vvData, szOutFileName, parmKind, sampPeriod);
	}

	// Write the feature vector as HTK format MFCC feature file and its lower and upper halves
	HRESULT WriteFeatureAsHTKMFCCAndItsHalves(const vector< vector<double> >& vvData, string szOutFileName, string szOutFileName_uhalf, string szOutFileName_lhalf)
	{
		return (HRESULT) GMMAdaptTool::WriteFeatureAsHTKMFCCAndItsHalves(vvData, szOutFileName, szOutFileName_uhalf, szOutFileName_lhalf);
	}

	HRESULT WriteFeatureAsDEBUGTRAJ(const vector< vector<double> >& vvData, string szOutFileName1, string szOutFileName2)
	{
		return (HRESULT) GMMAdaptTool::WriteFeatureAsDEBUGTRAJ(vvData, szOutFileName1, szOutFileName2);
	}
	HRESULT WriteFeatureAsDEBUGTRAJ(const vector< vector<double> >& vvData, string szOutFileName1)
	{
		return (HRESULT) GMMAdaptTool::WriteFeatureAsDEBUGTRAJ(vvData, szOutFileName1);
	}

	// Algorithms
	bool CompmatrixMeanVar(matrix<double>& X, matrix<double>& mu, matrix<double>& S, matrix<double>& invS)
	{
		int d = (int)X.ColNo();
		int N = (int)X.RowNo();
		double dbN = (double)N;
		mu.SetSize(1, d);
		mu.Null();
		for(int n1 = 0; n1 < N; n1++)
		{
			for(int n2 = 0; n2 < d; n2++)
				mu(0, n2) += X(n1, n2);
		}
		mu /= dbN;
		matrix<double> dXTmp(N, d);
		dXTmp.Null();
		for(int n1 = 0; n1 < N; n1++)
		{
			for(int n2 = 0; n2 < d; n2++)
				dXTmp(n1, n2) = X(n1, n2) - mu(0, n2);
		}
		S = ~dXTmp * dXTmp / (dbN - 1);
		for(int n2 = 0; n2 < d; n2++)
		{
			S(n2, n2) += 1.0e-6;
		}
		invS = !S;
		return true;
	}
	// 90 % CI test for outlier
	bool MCD(matrix<double>& X, vector<bool>& vbOutlier, double dbSigVal = 4.605, int h = 0)
	{
		int d = (int)X.ColNo();
		int N = (int)X.RowNo();
		if(N <= 0 || d <= 0)
		{
			cerr<<"Invalid data, dimension 0!"<<endl;
			return ERROR_INVALID_DATA;
		}
		//int dd = d*(d+1)/2;
		if(N < d)
		{
			cerr<<"A singular dataset, need more data points"<<endl;
			return false;
		}
		if(h < d)
		{
			int hh = (int)floor(0.75 * (double)N);			
			h = (d < hh)?hh:d;
		}
		double dbh = (double)h;
		vector<int> dummy(N);
		
		for(size_t nn = 0; nn < N ; nn++)
		{
			dummy[nn] = nn;
		}
		random_shuffle(dummy.begin(), dummy.end());
		matrix<double> XTmp(h, d);
		for(size_t nn = 0; nn < h ; nn++)
		{
			for(size_t dd = 0; dd < d; dd++)
			{
				XTmp(nn, dd) = X(dummy[nn], dd);
			}
		}
		
		matrix<double> mu1, s1_prev, s1, is1;
		CompmatrixMeanVar(XTmp, mu1, s1, is1);
		s1_prev = s1;
		bool bStop = false;
		vector<double> vmd(N, 0.0);
		vector<int> dummy2(N);
		vbOutlier.resize(N, false);
		do
		{
			for(size_t nn = 0; nn < N ; nn++)
				dummy2[nn] = nn;		
			
			// Compute sample mean
			// Compute Mahalanobis distance
			for(int n1 = 0; n1 < N; n1++)
			{
				matrix<double> xn1(1, d);
				for(int n2 = 0; n2 < d; n2++)
					xn1(0, n2) = X(n1, n2) - mu1(0, n2);
				vmd[n1] = (xn1 * is1 * ~xn1)(0, 0);
			}
			sort(dummy2.begin(), dummy2.end(), index_cmp_ascend< vector<double>& >(vmd));
			for(size_t nn = 0; nn < h ; nn++)
			{
				for(size_t dd = 0; dd < d; dd++)
				{
					XTmp(nn, dd) = X(dummy2[nn], dd);
				}
			}
			CompmatrixMeanVar(XTmp, mu1, s1, is1);
			double db1 = s1.Det();
			double db0 = s1_prev.Det();
			if(fabs(db1 - db0) < 1e-6)
			{
				bStop = true;
			}
			s1_prev = s1;
//			CompmatrixMeanVar(XTmp, mu1, s1, is1);
		}while(!bStop);
		for(int n = 0; n < N; n++)
		{
			if(dbSigVal < vmd[dummy2[n]])
			{		
				vbOutlier[dummy2[n]] = true;
			}
		}
		return true;
	}

	HRESULT MCD(vector< vector<double> >& vvX, vector<bool>& vbOutlier, double dbSigVal = 4.605, int h = 0)
	{
		size_t N = vvX.size();
		if(N < 1)
		{
			cerr<<"Empty Data for MCD!"<<endl;
			return ERROR_EMPTY;
		}
		size_t d0 = vvX[0].size();
		for(size_t n = 1; n < N; n++)
		{
			size_t d = vvX[n].size();
			if(d != d0)
			{
				cerr<<"Data dimension mismatched in MCD"<<endl;
				return ERROR_INVALID_DATA;
			}
		}
		matrix<double> X(N, d0);
		for(size_t n = 0; n < N; n++)
		{
			for(size_t d = 0; d < d0; d++)
			{
				X(n, d) = vvX[n][d];
			}
		}
		bool b = MCD(X, vbOutlier, dbSigVal, h);
		if(b)
			return S_OK;
		else
			return ERROR_INVALID_DATA;
	}

	HRESULT GetvvDelta(vector< vector<double> >& vvData, vector< vector<double> >& vvDData)
	{
		size_t N = vvData.size();
		if(N < 1)
		{
			cerr<<"Empty Data for GetvvDelta!"<<endl;
			return ERROR_EMPTY;
		}
		size_t d0 = vvData[0].size();
		for(size_t n = 1; n < N; n++)
		{
			size_t d = vvData[n].size();
			if(d != d0)
			{
				cerr<<"Data dimension mismatched in GetvvDelta"<<endl;
				return ERROR_INVALID_DATA;
			}
		}
		vvDData.resize(N, vector<double>(d0, 0.0));
		for(size_t n = 1; n < N; n++)
		{
			for(size_t d = 0; d < d0; d++)
			{
				vvDData[n][d] = vvData[n][d] - vvData[n-1][d];
			}
		}
		return S_OK;
	}


	HRESULT CheckTrajValidity(vector< vector<double> >& vvData, vector< vector<double> >& vvNewData)
	{
		size_t NN = vvData.size();
		if(NN < 16)
		{
			cerr<<"The length of the sequence is already too short in CheckTrajValidity"<<endl;
			vvNewData = vvData;
			return ERROR_INVALID_DATA;
		}
		vector< vector<double> > vvDData;
		HRESULT hr = GetvvDelta(vvData, vvDData);
		if(hr != S_OK)	return hr;
		vector<bool> vbOutlier, vbOutlierD;
		hr = MCD(vvData,vbOutlier);
		if(hr != S_OK)	return hr;
		hr = MCD(vvDData,vbOutlierD);
		if(hr != S_OK)	return hr;
		size_t N =  vbOutlier.size();
		vector<bool> vbNewOutlier(N, false);
		for(size_t n = 0; n < N; n++)
		{
			if(vbOutlier[n])
			{
				size_t leftindex = (n >= 1)?n:1;
				size_t rightindex = (n < N - 1)?(n+1):(N-1);
				vbNewOutlier[n] = vbOutlier[n] & (vbOutlier[leftindex] | vbOutlier[rightindex]);
			}
		}
		size_t nStart = 0;
		size_t nStop = N-1;

		for(size_t n = 0; n < N; n++)
		{
			if(vbNewOutlier[n] == false)
			{
				nStart = n;
				break;
			}
		}
		for(size_t n = nStop; n >= 0; n--)
		{
			if(vbNewOutlier[n] == false)
			{
				nStop = n;
				break;
			}
		}
		if(nStop-nStart < 15)
		{
			cerr<<"The sequence will beome too short after truncation, so aborted"<<endl;
			vvNewData = vvData;
			return ERROR_INVALID_DATA;			
		}
		vvNewData.resize(nStop - nStart + 1);
		for(size_t n = nStart; n < nStop; n++)
		{
			vvNewData[n - nStart] = vvData[n];
		}
		return S_OK;
	}


	HRESULT CheckTrajValidityBYU2D(vector< vector<double> >& vvOriginalData, vector< vector<double> >& vvNewData)
	{
		size_t NN = vvOriginalData.size();
		if(NN < 16)
		{
			cerr<<"The length of the sequence is already too short in CheckTrajValidityBYU2D"<<endl;
			vvNewData = vvOriginalData;
			return ERROR_INVALID_DATA;
		}
		vector< vector<double> > vvData(NN, vector<double>(2, 0.0) );
		for(size_t n = 0; n < NN ; n++)
		{
			vvData[n][0] = vvOriginalData[n][0];
			vvData[n][1] = vvOriginalData[n][1];
		}
		// Check validity only on first 2 dimensions
		vector< vector<double> > vvDData;
		HRESULT hr = GetvvDelta(vvData, vvDData);
		if(hr != S_OK)	return hr;
		vector<bool> vbOutlier, vbOutlierD;
		hr = MCD(vvData,vbOutlier);
		if(hr != S_OK)	return hr;
		hr = MCD(vvDData,vbOutlierD);
		if(hr != S_OK)	return hr;
		size_t N =  vbOutlier.size();
		vector<bool> vbNewOutlier(N, false);
		for(size_t n = 0; n < N; n++)
		{
			if(vbOutlier[n])
			{
				size_t leftindex = (n >= 1)?n:1;
				size_t rightindex = (n < N - 1)?(n+1):(N-1);
				vbNewOutlier[n] = vbOutlier[n] & (vbOutlier[leftindex] | vbOutlier[rightindex]);
			}
		}
		size_t nStart = 0;
		size_t nStop = N-1;

		for(size_t n = 0; n < N; n++)
		{
			if(vbNewOutlier[n] == false)
			{
				nStart = n;
				break;
			}
		}
		for(size_t n = nStop; n >= 0; n--)
		{
			if(vbNewOutlier[n] == false)
			{
				nStop = n;
				break;
			}
		}
		if(nStop-nStart < 15)
		{
			cerr<<"The sequence will beome too short after truncation, so aborted"<<endl;
			vvNewData = vvOriginalData;
			return ERROR_INVALID_DATA;			
		}
		vvNewData.resize(nStop - nStart + 1);
		for(size_t n = nStart; n <= nStop; n++)
		{
			vvNewData[n - nStart] = vvOriginalData[n];
		}
		return S_OK;
	}

 	// Divide the number of accumulated background by its number of frames
	 HRESULT AvgBackground(const CvArr* src, CvArr*& dst, unsigned int nnum)
	 {
		 if(dst == NULL)
		 {
			 if(CV_IS_MAT_HDR(src))
			 {
				 dst = cvCreateMat(((CvMat*)src)->rows, ((CvMat*)src)->cols, ((CvMat*)src)->type);
			 }
			 else if(CV_IS_IMAGE_HDR(src))
			 {
				 dst = cvCreateImage(cvSize(((IplImage*)src)->width, ((IplImage*)src)->height), ((IplImage*)src)->depth, ((IplImage*)src)->nChannels );
			 }
		 }
		 double dscale = 1.0;
		 if(nnum > 0)
		 {
			 dscale = 1.0/((double)(nnum));
		 }
		 cvConvertScale(src, dst, dscale);
		 return S_OK;
	 }

	void PreCompute2ClassesGMMTable(vector<GMM>& vGMMs, bool*** vvvDecisionTable)
	{
		//vvvDecisionTable = vector< vector< vector<bool> > >(256, vector< vector<bool> >(256, vector<bool>(256, false)));
		vector<double> vtmp(3, 0.0);
		for(int i = 0; i < 256; ++i)
		{
			for(int j = 0; j < 256; ++j)
			{
				for(int k = 0; k < 256; ++k)
				{
					vtmp[0] = (double)i;
					vtmp[1] = (double)j;
					vtmp[2] = (double)k;
					double db0 = vGMMs[0].CompLikelihood(vtmp);
					double db1 = vGMMs[1].CompLikelihood(vtmp);
					vvvDecisionTable[i][j][k] = (db1 > db0)?true:false;
				}
			}
		}
	}
	
	// Normalize image color by RGB histogram
	void NormalizeImageRGB(IplImage* pImg, IplImage*& pNewImg, double ptr = 0.01)
	{
		GMMAdaptTool::NormalizeImageRGB(pImg, pNewImg, ptr);
		/*int height = pImg->height;
		int width  = pImg->width;
		int nC = pImg->nChannels;
		int depth = pImg->depth;
		double dPts = (double)(width * height);
		double dptrnum = dPts * ptr;
		if(pNewImg != NULL)
		{
			cvReleaseImage(&pNewImg);
		}
		pNewImg = cvCreateImage(cvSize(width, height), depth, nC);
		vector< vector<double> > vHist(nC, vector<double>(256,0.0));
		vector<double> vtr(nC, 255.0);
		for(int n = 0; n < height; ++n)
		{
			for(int m = 0; m < width; ++m)
			{
				CvScalar s = cvGet2D(pImg, n, m);
				for(int c = 0; c < nC; ++c)
				{
					int nIndex = (int)float(s.val[c]+0.5);
					nIndex = (nIndex <0)?0:nIndex;
					nIndex = (nIndex>255)?255:nIndex;
					vHist[c][nIndex]++;
				}
			}
		}
		
		for(int c = 0; c < nC; ++c)
		{
			double dbtmpcb = 0.0;
			double dbtmptr = 255.0;
			for(int b = 255; b >=0; b--)
			{
				dbtmpcb += vHist[c][b];
				if(dbtmpcb >= dptrnum)
				{
					dbtmptr = (double)b;
					break;
				}
			}
			vtr[c] = dbtmptr;
		}
		
		for(int n = 0; n < height; ++n)
		{
			for(int m = 0; m < width; ++m)
			{
				CvScalar s = cvGet2D(pImg, n, m);
				for(int c = 0; c < nC; ++c)
				{
					if(s.val[c] >= vtr[c])
						s.val[c] = 255.0;
					else
						s.val[c] = s.val[c] / vtr[c] * 255;
				}
				cvSet2D(pNewImg, n, m, s);
			}
		}*/
	}

	void AutoLabelSLICRGBOnly(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int nsp)
	{
		return GMMAdaptTool::AutoLabelSLICRGBOnly(szPath, vGMMs, mask, nsp, m_nHandModelIndex);
	}

	void AutoLabel_SLIC_BS0(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int nsp, cv::Mat& ref)
	{
		return GMMAdaptTool::AutoLabel_SLIC_BS0(szPath, vGMMs, mask, nsp, m_nHandModelIndex, ref);
	}

	HRESULT OfflineFEImgSeq2_2_3_part_SLIC(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//

		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			AutoLabelSLICRGBOnly(szImgNamen, vGMMs, refmask, 200);
			iplrefmask = refmask.operator IplImage();
			SavecvImageForDebug(szPathName0, "Clr_16_16", szKeyName, n, 5, ".jpg", &iplrefmask);
		
			pcurr = vpImgs[n];
			CvSize cvsize = cvSize(pcurr->width*4, pcurr->height*4);
			pImgb = cvCreateImage(cvsize, IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvsize, IPL_DEPTH_32F, 1);
			pcurr2c = cvCreateImage(cvsize, IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvsize, IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvsize, IPL_DEPTH_32F, 1);

			// cvConvert(&iplrefmask2,pcurr2b);
			cvConvert(&iplrefmask,pcurr2b);

			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			cv::Mat cvCurrImgTempGray(cvCurrImgTemp.size(), CV_8UC1);
			cv::cvtColor(cvCurrImgTemp, cvCurrImgTempGray, CV_BGR2GRAY);
			IplImage curr_orig = cvCurrImgTempGray.operator IplImage();
			IplImage* pcurr_orig = &curr_orig;
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			IplImage* pMask2 = cvCreateImage(cvsize, IPL_DEPTH_8U, 1);
			//cvAnd(pcurr, pImgb, pMask2);
			cvAnd(pcurr_orig, pImgb, pMask2);
			
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			cvReleaseImage(&pImgb);
		}

		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				//WriteFeatureAsHTKMFCC(vvFeature, szOutFeatureFile);
				//WriteFeatureAsDEBUGTRAJ(szOutFeatureFile, szOutFileName1, szOutFileName2);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}

	HRESULT OfflineFEImgSeq2_2_3_part_SLIC_BS0(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				pBkGnd = cvCreateImage	(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				pBkGndTmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//

		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			AutoLabel_SLIC_BS0(szImgNamen, vGMMs, refmask, 200, BkGnd);

			//=============================================================================//

			iplrefmask = refmask.operator IplImage();
			//vector<double> comvec;
			//COMPoints(refmask2, comvec);
			//=============================================================================//
		
			pcurr = vpImgs[n];

			CvSize cvsize = cvSize(pcurr->width*4, pcurr->height*4);
			pImgb = cvCreateImage(cvsize, IPL_DEPTH_8U, 1);
			pcurr2b = cvCreateImage(cvsize, IPL_DEPTH_32F, 1);
			pcurr2c = cvCreateImage(cvsize, IPL_DEPTH_32F, 1);
			pcurr2b_tempb = cvCreateImage(cvsize, IPL_DEPTH_8U, 1);
			pcurr2b_tempb3 = cvCreateImage(cvsize, IPL_DEPTH_32F, 1);

			//cvConvert(&iplrefmask2,pcurr2b);
			cvConvert(&iplrefmask,pcurr2b);
			//=============================================================================//
			// Now, combine the mask result with ROI
			//=============================================================================//
			cvConvert(pcurr2b, pcurr2b_tempb);
			//cvAnd(pcurr, pcurr2b_tempb, pMask1);
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			cv::Mat cvCurrImgTempGray(cvCurrImgTemp.size(), CV_8UC1);
			cv::cvtColor(cvCurrImgTemp, cvCurrImgTempGray, CV_BGR2GRAY);
			IplImage curr_orig = cvCurrImgTempGray.operator IplImage();
			IplImage* pcurr_orig = &curr_orig;
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			cvResize(pcurr2b_tempb, pTmp, CV_INTER_NN);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			cvCopy(pcurr2b, pcurr2b_tempb3);
			//=============================================================================//

			//=============================================================================//
			// Morphological closing/ hole filling
			//=============================================================================//
			CloseImage(pcurr2b_tempb3, (CvArr*&)pcurr2c);
			cvConvert(pcurr2c, pImgb);
			//IplImage* pMask2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_8U, 1);
			IplImage* pMask2 = cvCreateImage(cvsize, IPL_DEPTH_8U, 1);
			//cvAnd(pcurr, pImgb, pMask2);
			cvAnd(pcurr_orig, pImgb, pMask2);
			SavecvImageForDebug(szPathName0, "ClosedMask", szKeyName, n, 5, ".jpg", pImgb);
			cvReleaseImage(&pMask2);
			cvReleaseImage(&pcurr2b_tempb);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			//=============================================================================//
			// Now, try to detect center of mass of the hand region along with
			// fingertip.
			// However, currently, finger detection is not good in in car dataset,
			// we only use the center of mass, so this part should be simplified in the 
			// future if only COM is used
			//=============================================================================//
			hrTrackFinger = S_OK;
			if(vvFeature.size() == 0)
			{
				hrTrackFinger = EstCOMandFingerLoc1(pImgb, vTmp);
			}
			else
			{
				hrTrackFinger = EstCOMandFingerLoc2(pImgb, vTmp, vvFeature[vvFeature.size() - 1]);
			}

			if( S_OK == hrTrackFinger)
			{
				size_t nWait = vHandVanishIndices.size();
				if(nWait > 0)
				{
					// simply discard those you cannot save
					if(0 < vvFeature.size())
					{
						cxn_k = vvFeature[vvFeature.size()-1][0];
						cyn_k = vvFeature[vvFeature.size()-1][1];
						dxn_k = vvFeature[vvFeature.size()-1][2];
						dyn_k = vvFeature[vvFeature.size()-1][3];
						cxn = vTmp[0];
						cyn = vTmp[1];
						dxn = vTmp[2];
						dyn = vTmp[3];					
						for(nk = 1; nk <= (int)nWait; nk++)
						{
							r = (double)nk / ((double)nWait + 1.0);
							vTmpTmp[0] = (cxn - cxn_k) * r + cxn_k;
							vTmpTmp[1] = (cyn - cyn_k) * r + cyn_k;
							vTmpTmp[2] = floor((dxn - dxn_k) * r + dxn_k + 0.5);
							vTmpTmp[3] = floor((dyn - dyn_k) * r + dyn_k + 0.5);
							vvFeature.push_back(vTmpTmp);
						}
					}
					nWait = 0;
					vHandVanishIndices.clear();
				}
				vvFeature.push_back(vTmp);
			}
			else
			{
				vHandVanishIndices.push_back(n);
				if(vHandVanishIndices.size() >= MAX_NUM_FRAME_WAIT_FOR_HAND_APPEAR_AGAIN)
				{	
					// Means the gesture ends in online FE				
				}
				// When it comes to the frame before the last frame,
				// we just assume the point stay at the last observation 
				// all the time
				//!!!! Can be wrong!!!!
				else if(n == nStop - 1)
				{
					for(nk = 1; nk <= (int)vHandVanishIndices.size(); nk++)
					{
						if(vTmp.size() == 4)
						{
							vvFeature.push_back(vTmp);
						}
					}
					vHandVanishIndices.clear();
				}
				else if(n == nStop)
				{
					// Do nothing, coz nothing is observed 
					vHandVanishIndices.clear();
				}
				cvCurrImgTemp.release();
			}
			refmask0.release();
			refmask.release();
			refmask2.release();
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			cvReleaseImage(&pImgb);
		}
		cvReleaseImage(&pBkGnd);
		// Write features
		if(0 < vvFeature.size())
		{
			HRESULT hrlast = WriteFeatureAsHTKMFCCAndItsHalves(vvFeature, szOutFeatureFile, szOutFeatureFile_uhalf, szOutFeatureFile_lhalf);
			if(hrlast == S_OK)
			{
				vFEN[nID].IsCreated = true;
				string szOutFileName001, szOutFileName002;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				MakeDebugOutputName(szPathName0, "FNG", szKeyName, ".jpg", szOutFileName002);
				WriteFeatureAsDEBUGTRAJ(vvFeature, szOutFileName001, szOutFileName002);
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}

	HRESULT OfflineFEImgSeq2_2_3_part_SLIC_BS0_Fusion(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		//IplImage* pcurr2bcpy = 0;
		IplImage* pcurr2b_tempb = 0;
		//IplImage* pcurr2b_tempb2 = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		cout<<"Begin Processing Feature file "<<szKeyName<<endl;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				IplImage* pBkGnd_tmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				IplImage* pBkGndTmp_tmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				// ============================================================================================================/
				// For color normalization, brute force!, or pBkGnd = pBkGnd_tmp, pBkGndTmp_tmp = pBkGndTmp
				// ============================================================================================================/
				NormalizeImageRGB(pBkGnd_tmp, pBkGnd);
				NormalizeImageRGB(pBkGndTmp_tmp, pBkGndTmp);
				cvReleaseImage(&pBkGnd_tmp);
				cvReleaseImage(&pBkGndTmp_tmp);
				// ============================================================================================================/
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//

		//=============================================================================//
		// MAP for frames
		//=============================================================================//
		vector<string> vszImagePathNames;
		vector<cv::Mat> vMasks1;
		vector<cv::Mat> vMasks2;
		vector<cv::Mat> vMasks3;
		vector<cv::Mat> vMasks01;
		for(n = 0/*nStopBkGnd + 1*/; n < vpImgs.size()/*= nStop*/; n++)
		{
			vszImagePathNames.push_back(ImgPathNames[n]);
		}
		if(m_nHandModelIndex == 0)
			m_nHandModelIndex = 1;
		vector<cv::Mat> vScore0_1;
		vector<cv::Mat> vScore1_1;
		vector<cv::Mat> vScore0_2;
		vector<cv::Mat> vScore1_2;
		vector<cv::Mat> vScore0_3;
		vector<cv::Mat> vScore1_3;
		
		vector<cv::Mat> vLLR_1;
		vector<cv::Mat> vCorr_1;
		vector<cv::Mat> vKS_1;
		vector<cv::Mat> vLLR_2;
		vector<cv::Mat> vCorr_2;
		vector<cv::Mat> vKS_2;
		vector<cv::Mat> vLLR_3;
		vector<cv::Mat> vCorr_3;
		vector<cv::Mat> vKS_3;

		//GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC
		//	(vszImagePathNames, 
		//	vGMMs, 
		//	vMasks2, 
		//	200, 
		//	m_nHandModelIndex, 
		//	m_szGMMPath, 
		//	BkGnd, 
		//	vFusionGMMs, 
		//	bFrameBatch, 
		//	nAdaptMode, 
		//	nStopBkGnd + 1, 
		//	nStop + 1, 
		//	vScore0_2, 
		//	vScore1_2, 
		//	vLLR_2, 
		//	vCorr_2);
		GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC
			(vszImagePathNames, 
			vGMMs, 
			vMasks2, 
			200, 
			m_nHandModelIndex, 
			m_szGMMPath, 
			BkGnd, 
			vFusionGMMs, 
			bFrameBatch, 
			nAdaptMode, 
			nStopBkGnd + 1, 
			nStop + 1, 
			vScore0_2, 
			vScore1_2, 
			vLLR_2, 
			vCorr_2,
			vKS_2);

		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			//cout<<"Doing frame "<<n<<endl;
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			refmask = vMasks2[n];// - nStopBkGnd - 1];
			IplImage iplrefmask1 = refmask.operator IplImage();
			IplImage* prefmask1 = &iplrefmask1;
			SaveMATForDebug(szPathName0, "mask_200_sp", szKeyName, n, 5, ".dat", prefmask1);
			// Save scores from LLR and Corr

			IplImage iplsLLR2 = vLLR_2[n].operator IplImage();
			IplImage* piplsLLR2 = &iplsLLR2;
			SaveMATForDebug(szPathName0, "sLLR_200_sp", szKeyName, n, 5, ".dat", piplsLLR2);

			IplImage iplsCorr2 = vCorr_2[n].operator IplImage();
			IplImage* piplsCorr2 = &iplsCorr2;
			SaveMATForDebug(szPathName0, "Corr_200_sp", szKeyName, n, 5, ".dat", piplsCorr2);
			

			IplImage iplscore0_2 = vScore0_2[n].operator IplImage();
			IplImage iplscore1_2 = vScore1_2[n].operator IplImage();
			IplImage* piplscore0_2 = &iplscore0_2;
			IplImage* piplscore1_2 = &iplscore1_2;
			SaveMATForDebug(szPathName0, "mask_200_sp_score0", szKeyName, n, 5, ".dat", piplscore0_2);
			SaveMATForDebug(szPathName0, "mask_200_sp_score1", szKeyName, n, 5, ".dat", piplscore1_2);
			//=============================================================================//

			iplrefmask = refmask.operator IplImage();
			
			//=============================================================================//
			// Save masked images
			//=============================================================================//
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask1 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask1, CV_RGB(64, 128, 64));
			IplImage* piplrefmask = &iplrefmask;
			cvResize(piplrefmask, pTmp, CV_INTER_NN);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask1, pTmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask1);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask1);
			//=============================================================================//

			vector<double> comvec;
			// ===============================================================================//
			//COMPoints(refmask2, comvec);
			COMPoints(refmask, comvec);
			vvFeature.push_back(comvec);
			//=============================================================================//
		
			cvReleaseImage(&pcurr2b_tempb);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			refmask.release();
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			cvReleaseImage(&pImgb);
		}
		cvReleaseImage(&pBkGnd);

		if(0 < vvFeature.size())
		{
			vector< vector<double> > vvFeature2 = vvFeature;
			vector<int> lastcom(vvFeature.size(), 0);
			vector<int> nextcom(vvFeature.size(), 0);
			int lastcomtmp = -1;
			int nextcomtmp = -1;
			for(int nn1 = 0; nn1 < (int)vvFeature.size(); nn1++)
			{
				if(vvFeature[nn1][0] > 0 && vvFeature[nn1][1] > 0)
					lastcomtmp = nn1;
				lastcom[nn1] = lastcomtmp;
			}

			for(int nn1 = (int)vvFeature.size() - 1; nn1 >=0 ; nn1--)
			{
				if(vvFeature[nn1][0] > 0 && vvFeature[nn1][1] > 0)
					nextcomtmp = (int)nn1;
				nextcom[nn1] = nextcomtmp;
			}
			// interpolate COM
			for(int nn1 = 0; nn1 < (int)vvFeature.size(); nn1++)
			{
				if(lastcom[nn1] >= 0 && nextcom[nn1] >= 0 && nextcom[nn1] != lastcom[nn1])
				{
					double distcomtmpd = (double)(nextcom[nn1] - lastcom[nn1]);
					double distcomtmpn = (double)(nn1 - lastcom[nn1]);
					vector<double> distcomtmpv(2, 0.0);
					distcomtmpv[0] = vvFeature[nextcom[nn1]][0] - vvFeature[lastcom[nn1]][0];
					distcomtmpv[1] = vvFeature[nextcom[nn1]][1] - vvFeature[lastcom[nn1]][1];
					vvFeature2[nn1][0] = vvFeature[nextcom[nn1]][0] + distcomtmpv[0]  * distcomtmpn / distcomtmpd;	
					vvFeature2[nn1][1] = vvFeature[nextcom[nn1]][1] + distcomtmpv[1]  * distcomtmpn / distcomtmpd;	
				}
				else if(lastcom[nn1] < 0 && nextcom[nn1] >= 0)
				{
					vvFeature2[nn1][0] = vvFeature[nextcom[nn1]][0];
					vvFeature2[nn1][1] = vvFeature[nextcom[nn1]][1];
				}
				else if(lastcom[nn1] >= 0 && nextcom[nn1] < 0)
				{
					vvFeature2[nn1][0] = vvFeature[lastcom[nn1]][0];
					vvFeature2[nn1][1] = vvFeature[lastcom[nn1]][1];
				}
			}
			//cout<<"ready to write HTK feature\n";
			HRESULT hrlast = WriteFeatureAsHTKMFCC(vvFeature2, szOutFeatureFile_uhalf);
			//cout<<"HTK feature written\n";
			if(hrlast == S_OK)
			{
				//cout<<"before vFEN[nID].IsCreated = true;"<<endl;
				vFEN[nID].IsCreated = true;
				//cout<<"before string szOutFileName001, szOutFileName002;"<<endl;
				string szOutFileName001, szOutFileName002;
				//cout<<"before MakeDebugOutputName(szPathName0, \"COM\", szKeyName, \".jpg\", szOutFileName001);"<<endl;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				//cout<<"before writing trajectory to jpg\n";
				WriteFeatureAsDEBUGTRAJ(vvFeature2, szOutFileName001);
				//cout<<"after writing trajectory to jpg\n";
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}
	
	HRESULT OfflineFEImgSeq2_2_3_part_BS0_Fusion_SLIC_MASK(DWORD nID, hash_map<int, FENode>& vFEN, vector<IplImage*>& vpImgs, hash_map<string, vector<VideoNode> >& htVN, vector<string>& ImgPathNames, vector<GMM>& vGMMs, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode)
	{
	#ifdef _WIN32
		string delim = "\\";
	#else
		string delim = "/";
	#endif
		string szOutFeatureFile = vFEN[nID].szOutFeatureFile;
		string szOutFeatureFile_uhalf = vFEN[nID].szOutFeatureFile_uhalf;
		string szOutFeatureFile_lhalf = vFEN[nID].szOutFeatureFile_lhalf;
		string szMLFName = vFEN[nID].szMLFName;
		string szLabelName = vFEN[nID].szLabelName;
		string szKeyName = vFEN[nID].szKeyName; 
		string szErrLog = vFEN[nID].szErrLog;
		size_t noutfeslash = szOutFeatureFile.find_last_of(delim);
		string szPathName0 = szOutFeatureFile.substr(0, noutfeslash);

		IplImage* pcurr = 0;
		IplImage* pBkGndTmp = 0;
		IplImage* pBkGnd = 0;
		IplImage* pRef1 = 0;
		IplImage* pRef2 = 0;
		IplImage* pcurr2 = 0;
		IplImage* pcurr2c = 0;
		IplImage* pcurr2b = 0;
		IplImage* pcurr2b_tempb = 0;
		IplImage* pcurr2b_tempb3 = 0;
		IplImage* pImgb = 0;
		IplImage* pImgbs = 0;
		IplImage iplrefmask0;
		IplImage iplrefmask;
		IplImage iplrefmask2;
		CvMemStorage* storage = 0;
		CvSeq* first_contour = 0;
		CvSeq* contour = 0;
		cv::Mat curr2b;
		cv::Mat curr2b2;
		vector< vector<double> > vvFeature;
		vector<double> vTmpTmp(4);
		vector<double> vTmp;
		vector<int> vHandVanishIndices;
		vector<CvRect> vRect;
		vector<CvRect> vAggRect;
		CvRect rectROI;
		CvRect rect1;
		CvRect recti;
		IplImage* pImgbSWM = 0;
		IplImage* pSWMask = 0;

		HRESULT hr = S_OK;
		HRESULT hrTrackFinger = S_OK;
		int nStopBkGnd = 0;
		int nStart = 0;
		int nStop = (int)vpImgs.size() - 1;
		int n = 0;
		int kk1 = 0;
		int kk2 = 0;
		int nn1 = 0;
		int nn2 = 0;
		int nk = 0;
		int nrectnum = 0;
		int nROIIndex = 0;
		
		
		double db1 = 0.0;
		double cxn_k = 0.0;
		double cyn_k = 0.0;
		double dxn_k = 0.0;
		double dyn_k = 0.0;
		double cxn = 0.0;
		double cyn = 0.0;;
		double dxn = 0.0;
		double dyn = 0.0;					
		double r = 0.0;
		cout<<"Begin Processing Feature file "<<szKeyName<<endl;
		//=============================================================================//
		// Begin feature extraction
		//=============================================================================//
		if(nStop <= 0)
		{
			return ERROR_EMPTY;
		}
		//=============================================================================//

		//=============================================================================//
		// Check where hand is inside the frame and outside 
		//=============================================================================//
		hr = GMMAdaptTool::TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
		if(hr == ERROR_EMPTY)
		{
			return hr;
		}
		//=============================================================================//


		//=============================================================================//
		//Accumulate and average background images
		//=============================================================================//
		for(n = nStart; n <= nStopBkGnd; n++)
		{
			pcurr = cvLoadImage(ImgPathNames[n].c_str()); //vpImgs[n];
			if(n == nStart)
			{
				IplImage* pBkGnd_tmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				IplImage* pBkGndTmp_tmp = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
				// ============================================================================================================/
				// For color normalization, brute force!, or pBkGnd = pBkGnd_tmp, pBkGndTmp_tmp = pBkGndTmp
				// ============================================================================================================/
				NormalizeImageRGB(pBkGnd_tmp, pBkGnd);
				NormalizeImageRGB(pBkGndTmp_tmp, pBkGndTmp);
				cvReleaseImage(&pBkGnd_tmp);
				cvReleaseImage(&pBkGndTmp_tmp);
				// ============================================================================================================/
				cvSet(pBkGnd, cvScalar(0.0));
				cvSet(pBkGndTmp, cvScalar(0.0));
			}
			pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, pcurr->nChannels);
			cvConvert(pcurr, pcurr2);
			GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGndTmp);
			cvReleaseImage(&pcurr);
			cvReleaseImage(&pcurr2);
		}
		AvgBackground(pBkGndTmp, (CvArr*&)pBkGnd, nStopBkGnd - nStart + 1);
		cvReleaseImage(&pBkGndTmp);
		cv::Mat BkGnd(pBkGnd);
		//=============================================================================//

		//=============================================================================//
		// MAP for frames
		//=============================================================================//
		vector<string> vszImagePathNames;
		vector<cv::Mat> vMasks1;
		vector<cv::Mat> vMasks2;
		vector<cv::Mat> vMasks3;
		vector<cv::Mat> vMasks01;
		for(n = 0/*nStopBkGnd + 1*/; n < vpImgs.size()/*= nStop*/; n++)
		{
			vszImagePathNames.push_back(ImgPathNames[n]);
		}
		if(m_nHandModelIndex == 0)
			m_nHandModelIndex = 1;
		vector<cv::Mat> vScore0_1;
		vector<cv::Mat> vScore1_1;
		vector<cv::Mat> vScore0_2;
		vector<cv::Mat> vScore1_2;
		vector<cv::Mat> vScore0_3;
		vector<cv::Mat> vScore1_3;
		
		vector<cv::Mat> vLLR_1;
		vector<cv::Mat> vCorr_1;
		vector<cv::Mat> vLLR_2;
		vector<cv::Mat> vCorr_2;
		vector<cv::Mat> vLLR_3;
		vector<cv::Mat> vCorr_3;
		
		// Find SLIC mask first.
		vector<string> vszLoadedFiles;
		vector< vector< vector<cv::Point> > > vvvSLICLabels;
		vector< vector< vector<double> > > vvvFeatures_SLIC;
		vector< vector< vector<double> > > vvvCorrs_SLIC;
		vector< vector<double> > vvKSs_SLIC;
		vector<cv::Mat> vSLICMasks;
		//GMMAdaptTool::FindSLIC(vszImagePathNames, vszLoadadFiles, vvvSLICLabels, 200);
		//GMMAdaptTool::FindSLICCorrs(vszImagePathNames, vszLoadedFiles, vvvFeatures_SLIC, vvvSLICLabels, vvvCorrs_SLIC, vSLICMasks, BkGnd, nStopBkGnd + 1, nStop + 1, 200);
		GMMAdaptTool::FindSLICCorrs(vszImagePathNames, vszLoadedFiles, vvvFeatures_SLIC, vvvSLICLabels, vvvCorrs_SLIC, vvKSs_SLIC, vSLICMasks, BkGnd, nStopBkGnd + 1, nStop + 1, 200);
					
		GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC_MASK(vszImagePathNames, vGMMs, vMasks1, 32, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_1, vScore1_1, vSLICMasks, vLLR_1, vCorr_1);
		GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC_MASK(vszImagePathNames, vGMMs, vMasks2, 16, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_2, vScore1_2, vSLICMasks, vLLR_2, vCorr_2);
		GMMAdaptTool::AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC_MASK(vszImagePathNames, vGMMs, vMasks3, 4, m_nHandModelIndex, m_szGMMPath, BkGnd, vFusionGMMs, bFrameBatch, nAdaptMode, nStopBkGnd + 1, nStop + 1, vScore0_3, vScore1_3, vSLICMasks, vLLR_3, vCorr_3);

		//=============================================================================//
		// Do feature extraction for those frames with hands
		//=============================================================================//
		for(n = nStopBkGnd + 1; n <= nStop; n++)
		{
			string szImgNamen = ImgPathNames[n];
			cv::Mat refmask, refmaskq, refmask2, refmask0;
			//=============================================================================//
			// AutoLabelPartSM: label image in 3 different block sizes with feature space 
			// stochastic matching:
			// AutoLabelPart: direct label with current GMM in 3 different block sizes
			// around 1.5% performance gain using AutoLabelPartSM
			//=============================================================================//
			refmask0 = vMasks1[n];// - nStopBkGnd - 1];
			refmask = vMasks2[n];// - nStopBkGnd - 1];
			refmask2 = vMasks3[n];// - nStopBkGnd - 1];
			IplImage iplrefmask0 = refmask0.operator IplImage();
			IplImage iplrefmask1 = refmask.operator IplImage();
			IplImage iplrefmask2 = refmask2.operator IplImage();
			IplImage* prefmask0 = &iplrefmask0;
			IplImage* prefmask1 = &iplrefmask1;
			IplImage* prefmask2 = &iplrefmask2;
			SaveMATForDebug(szPathName0, "mask_32_32", szKeyName, n, 5, ".dat", prefmask0);
			SaveMATForDebug(szPathName0, "mask_16_16", szKeyName, n, 5, ".dat", prefmask1);
			SaveMATForDebug(szPathName0, "mask_4_4", szKeyName, n, 5, ".dat", prefmask2);
			// Save scores from LLR and Corr
			IplImage iplsLLR1 = vLLR_1[n].operator IplImage();
			IplImage* piplsLLR1 = &iplsLLR1;
			SaveMATForDebug(szPathName0, "sLLR_32_32", szKeyName, n, 5, ".dat", piplsLLR1);

			IplImage iplsLLR2 = vLLR_2[n].operator IplImage();
			IplImage* piplsLLR2 = &iplsLLR2;
			SaveMATForDebug(szPathName0, "sLLR_16_16", szKeyName, n, 5, ".dat", piplsLLR2);

			IplImage iplsLLR3 = vLLR_3[n].operator IplImage();
			IplImage* piplsLLR3 = &iplsLLR3;
			SaveMATForDebug(szPathName0, "sLLR_4_4", szKeyName, n, 5, ".dat", piplsLLR3);

			IplImage iplsCorr1 = vCorr_1[n].operator IplImage();
			IplImage* piplsCorr1 = &iplsCorr1;
			SaveMATForDebug(szPathName0, "Corr_32_32", szKeyName, n, 5, ".dat", piplsCorr1);

			IplImage iplsCorr2 = vCorr_2[n].operator IplImage();
			IplImage* piplsCorr2 = &iplsCorr2;
			SaveMATForDebug(szPathName0, "Corr_16_16", szKeyName, n, 5, ".dat", piplsCorr2);

			IplImage iplsCorr3 = vCorr_3[n].operator IplImage();
			IplImage* piplsCorr3 = &iplsCorr3;
			SaveMATForDebug(szPathName0, "Corr_4_4", szKeyName, n, 5, ".dat", piplsCorr3);

			

			IplImage iplscore0_1 = vScore0_1[n].operator IplImage();
			IplImage iplscore1_1 = vScore1_1[n].operator IplImage();
			IplImage iplscore0_2 = vScore0_2[n].operator IplImage();
			IplImage iplscore1_2 = vScore1_2[n].operator IplImage();
			IplImage iplscore0_3 = vScore0_3[n].operator IplImage();
			IplImage iplscore1_3 = vScore1_3[n].operator IplImage();
			IplImage* piplscore0_1 = &iplscore0_1;
			IplImage* piplscore1_1 = &iplscore1_1;
			IplImage* piplscore0_2 = &iplscore0_2;
			IplImage* piplscore1_2 = &iplscore1_2;
			IplImage* piplscore0_3 = &iplscore0_3;
			IplImage* piplscore1_3 = &iplscore1_3;
			SaveMATForDebug(szPathName0, "mask_32_32_score0", szKeyName, n, 5, ".dat", piplscore0_1);
			SaveMATForDebug(szPathName0, "mask_32_32_score1", szKeyName, n, 5, ".dat", piplscore1_1);
			SaveMATForDebug(szPathName0, "mask_16_16_score0", szKeyName, n, 5, ".dat", piplscore0_2);
			SaveMATForDebug(szPathName0, "mask_16_16_score1", szKeyName, n, 5, ".dat", piplscore1_2);
			SaveMATForDebug(szPathName0, "mask_4_4_score0", szKeyName, n, 5, ".dat", piplscore0_3);
			SaveMATForDebug(szPathName0, "mask_4_4_score1", szKeyName, n, 5, ".dat", piplscore1_3);


			//=============================================================================//

			//=============================================================================//
			// Then make all 3 labels consistent, that is,
			// If 4x4 block some blocks are detected as skin, 
			// but corresponding and adjacent 16x16 block 
			// didn't have skin label, then discard it.
			//=============================================================================//
			iplrefmask0 = refmask0.operator IplImage();
			iplrefmask = refmask.operator IplImage();
			iplrefmask2 = refmask2.operator IplImage();
			for(kk1 = 0; kk1 < iplrefmask.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 2.0);
					nn2 = (int)floor((double)kk2 / 2.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask0.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask0.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask0, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask0, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask0, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask0, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask0, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask, kk1, kk2, 0.0);
						}
					}
				}
			}

			for(kk1 = 0; kk1 < iplrefmask2.height; kk1++)
			{
				for(kk2 = 0; kk2 < iplrefmask2.width; kk2++)
				{
					db1 = cvGetReal2D(&iplrefmask2, kk1, kk2);
					nn1 = (int)floor((double)kk1 / 4.0);
					nn2 = (int)floor((double)kk2 / 4.0);
					cv::Point pt0 = cv::Point(nn1, nn2);
					cv::Point pt1 = cv::Point(max(0, nn1-1), nn2);
					cv::Point pt2 = cv::Point(nn1, max(0, nn2-1));
					cv::Point pt3 = cv::Point(nn1, min(iplrefmask.width - 1, nn2 + 1));
					cv::Point pt4 = cv::Point(min(iplrefmask.height - 1, nn1 + 1), nn2);

					if(cvGetReal2D(&iplrefmask2, kk1, kk2))
					{
						if(cvGetReal2D(&iplrefmask, pt0.x, pt0.y) + 
							cvGetReal2D(&iplrefmask, pt1.x, pt1.y) + 
							cvGetReal2D(&iplrefmask, pt2.x, pt2.y) + 
							cvGetReal2D(&iplrefmask, pt3.x, pt3.y) + 
							cvGetReal2D(&iplrefmask, pt4.x, pt4.y))
						{
						}
						else
						{
							cvSetReal2D(&iplrefmask2, kk1, kk2, 0.0);
						}
					}
				}
			}

			//=============================================================================//
			// Save masked images
			//=============================================================================//
			cv::Mat cvCurrImgTemp = cv::imread(szImgNamen);
			
			IplImage* pTmp = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 1);
			IplImage* pMask2 = cvCreateImage(cvSize(cvCurrImgTemp.cols, cvCurrImgTemp.rows), IPL_DEPTH_8U, 3);
			cvSet(pMask2, CV_RGB(64, 128, 64));
			IplImage* piplrefmask2 = &iplrefmask2;
			cvResize(piplrefmask2, pTmp, CV_INTER_NN);
			IplImage iplcurrtemp = cvCurrImgTemp.operator IplImage();
			IplImage* piplcurrtemp = &iplcurrtemp;
			cvCopy(piplcurrtemp, pMask2, pTmp);
			SavecvImageForDebug(szPathName0, "ClrMask", szKeyName, n, 5, ".jpg", pMask2);
			cvReleaseImage(&pTmp);
			cvReleaseImage(&pMask2);
			//=============================================================================//


			vector<double> comvec;
			// ===============================================================================//
			// Output combined mask for durther processing
			// ===============================================================================//
			SaveMATForDebug(szPathName0, "mask_combine", szKeyName, n, 5, ".dat", prefmask2);
			// ===============================================================================//
			COMPoints(refmask2, comvec);
			vvFeature.push_back(comvec);
			//=============================================================================//
		
			cvReleaseImage(&pcurr2b_tempb);
			cvReleaseImage(&pcurr2b_tempb3);
			//=============================================================================//
			
			refmask0.release();
			refmask.release();
			refmask2.release();
			cvReleaseImage(&pcurr2c);
			cvReleaseImage(&pcurr2b);
			cvReleaseImage(&pImgb);

		}
		cvReleaseImage(&pBkGnd);

		if(0 < vvFeature.size())
		{
			vector< vector<double> > vvFeature2 = vvFeature;
			vector<int> lastcom(vvFeature.size(), 0);
			vector<int> nextcom(vvFeature.size(), 0);
			int lastcomtmp = -1;
			int nextcomtmp = -1;
			for(int nn1 = 0; nn1 < (int)vvFeature.size(); nn1++)
			{
				if(vvFeature[nn1][0] > 0 && vvFeature[nn1][1] > 0)
					lastcomtmp = nn1;
				lastcom[nn1] = lastcomtmp;
			}

			for(int nn1 = (int)vvFeature.size() - 1; nn1 >=0 ; nn1--)
			{
				if(vvFeature[nn1][0] > 0 && vvFeature[nn1][1] > 0)
					nextcomtmp = (int)nn1;
				nextcom[nn1] = nextcomtmp;
			}
			// interpolate COM
			for(int nn1 = 0; nn1 < (int)vvFeature.size(); nn1++)
			{
				if(lastcom[nn1] >= 0 && nextcom[nn1] >= 0 && nextcom[nn1] != lastcom[nn1])
				{
					double distcomtmpd = (double)(nextcom[nn1] - lastcom[nn1]);
					double distcomtmpn = (double)(nn1 - lastcom[nn1]);
					vector<double> distcomtmpv(2, 0.0);
					distcomtmpv[0] = vvFeature[nextcom[nn1]][0] - vvFeature[lastcom[nn1]][0];
					distcomtmpv[1] = vvFeature[nextcom[nn1]][1] - vvFeature[lastcom[nn1]][1];
					vvFeature2[nn1][0] = vvFeature[nextcom[nn1]][0] + distcomtmpv[0]  * distcomtmpn / distcomtmpd;	
					vvFeature2[nn1][1] = vvFeature[nextcom[nn1]][1] + distcomtmpv[1]  * distcomtmpn / distcomtmpd;	
				}
				else if(lastcom[nn1] < 0 && nextcom[nn1] >= 0)
				{
					vvFeature2[nn1][0] = vvFeature[nextcom[nn1]][0];
					vvFeature2[nn1][1] = vvFeature[nextcom[nn1]][1];
				}
				else if(lastcom[nn1] >= 0 && nextcom[nn1] < 0)
				{
					vvFeature2[nn1][0] = vvFeature[lastcom[nn1]][0];
					vvFeature2[nn1][1] = vvFeature[lastcom[nn1]][1];
				}
			}
			//cout<<"ready to write HTK feature\n";
			HRESULT hrlast = WriteFeatureAsHTKMFCC(vvFeature2, szOutFeatureFile_uhalf);
			//cout<<"HTK feature written\n";
			if(hrlast == S_OK)
			{
				//cout<<"before vFEN[nID].IsCreated = true;"<<endl;
				vFEN[nID].IsCreated = true;
				//cout<<"before string szOutFileName001, szOutFileName002;"<<endl;
				string szOutFileName001, szOutFileName002;
				//cout<<"before MakeDebugOutputName(szPathName0, \"COM\", szKeyName, \".jpg\", szOutFileName001);"<<endl;
				MakeDebugOutputName(szPathName0, "COM", szKeyName, ".jpg", szOutFileName001);
				//cout<<"before writing trajectory to jpg\n";
				WriteFeatureAsDEBUGTRAJ(vvFeature2, szOutFileName001);
				//cout<<"after writing trajectory to jpg\n";
				size_t nSizet = vvFeature.size() - 1;
				nSizet *= 666666;
				htVN[szKeyName][0].nStart = 0;
				htVN[szKeyName][0].nStop = (int)nSizet;
				htVN[szKeyName][0].szName = szLabelName;
			}
			else if(hrlast == ERROR_INVALID_DATA)
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				cerr<<szOutFeatureFile<<" is not created due to all the same on some coordinate!\n";
				feout.close();
			}
			else
			{
				vFEN[nID].IsCreated = false;
				ofstream feout;
				feout.open(szErrLog.c_str(), ios_base::app);
				feout<<szOutFeatureFile<<" is not created for some other reason!\n";
				cerr<<szOutFeatureFile<<" is not created for some other reason!\n";
				feout.close();
			}
		}
		else
		{
			vFEN[nID].IsCreated = false;
			ofstream feout;
			feout.open(szErrLog.c_str(), ios_base::app);
			feout<<szOutFeatureFile<<" is not created!\n";
			cerr<<szOutFeatureFile<<" is not created!\n";
			feout.close();
		}
		return hr;
	}


	vector<int> vDimConsidered;
	int m_nHandModelIndex;
	int m_nBkGndModelIndex;
	int m_nK;
	double m_duth;
	string m_szGMMPath;
};
