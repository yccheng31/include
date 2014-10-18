#pragma once

#include "Defs.h"


// Skin color GMM class
class GMM
{
//friend class GMMHyperParam;
public:
	GMM()
	{
		m_nVecSize = 0;
		m_nMixes = 0;
		m_szCovType = "diagonal";
		//m_ptheta = NULL;
		m_vvTransp = vector< vector<double> >(3, vector<double>(3, 0.0));
	}
	GMM(int nVecSize, int nMixes = 1, string szCovType = "diagonal")
	{
		m_nVecSize = nVecSize;
		m_nMixes = nMixes;
		m_vmdMean.resize(m_nMixes);
		if(0 == m_vmdCovar.size())
			m_vmdCovar.resize(m_nMixes);
		if(0 == m_vmdInvCovar.size())
			m_vmdInvCovar.resize(m_nMixes);
		m_szCovType = szCovType;
		m_vMixWeights = vector<double>(nMixes, (1.0/(double)nMixes));
		m_vlogMixWeights = vector<double>(nMixes, log(1.0/(double)nMixes));
		for(int n = 0; n < m_nMixes; n++)
		{
			m_vmdMean[n].SetSize(m_nVecSize, 1);
			m_vmdMean[n].Null();
			m_vmdCovar[n].SetSize(m_nVecSize, m_nVecSize);
			m_vmdInvCovar[n].SetSize(m_nVecSize, m_nVecSize);
			m_vmdCovar[n].Null();
			m_vmdInvCovar[n].Null();
			for(int d = 0; d < m_nVecSize; d++)
			{
				m_vmdCovar[n](d,d) = 1.0;
				m_vmdInvCovar[n](d,d) = 1.0;
			}
		}
		m_vvTransp = vector< vector<double> >(3, vector<double>(3, 0.0));
		//m_ptheta = NULL;
	}
	GMM(const GMM& rhs)
	{
		m_nVecSize = rhs.m_nVecSize;
		m_nMixes = rhs.m_nMixes;
		m_szCovType = rhs.m_szCovType;
		m_szFeatureType = rhs.m_szFeatureType;
		m_vdgconst = rhs.m_vdgconst;
		m_vMixWeights = rhs.m_vMixWeights;
		m_vmdMean = rhs.m_vmdMean;
		m_vmdCovar = rhs.m_vmdCovar;
		m_vmdInvCovar = rhs.m_vmdInvCovar;
		m_vlogMixWeights = rhs.m_vlogMixWeights;
		m_vvddiaggconst = rhs.m_vvddiaggconst;
		m_vvTransp = rhs.m_vvTransp;
		//m_ptheta = NULL;
	}
	~GMM()
	{
		//if(m_ptheta != NULL)
		//{
		//	delete m_ptheta;
		//}
	}
	double LogAdd(double x, double y)
	{
	   double temp,diff,z;
	   if (x<y) {
		  temp = x; x = y; y = temp;
	   }
	   diff = y-x;
	   if (diff<minLogExp) 
		  return  (x<LSMALL)?LZERO:x;
	   else {
		  z = exp(diff);
		  return x+log(1.0+z);
	   }
	}
	void InitGMMMeans(int nVecSize, int nMixes = 1)
	{
		m_nVecSize = nVecSize;
		m_nMixes = nMixes;
		m_vmdMean.resize(m_nMixes);
		double oneovern = 1.0/(double)nMixes;
		double log1n = log(oneovern);
		m_vMixWeights = vector<double>(nMixes, oneovern);
		m_vlogMixWeights = vector<double>(nMixes, log1n);
		for(int n = 0; n < m_nMixes; n++)
		{
			m_vmdMean[n].SetSize(m_nVecSize, 1);
			m_vmdMean[n].Null();
		}
	}
	void InitGMMVars(int nVecSize, int nMixIndex, string szCovType = "diagonal")
	{
		m_nVecSize = nVecSize;
		m_szCovType = szCovType;
		if(nMixIndex >= m_vmdCovar.size())
		{
			m_vmdCovar.resize(nMixIndex + 1);
			m_vmdInvCovar.resize(nMixIndex + 1);
		}
		for(int n = 0; n < m_nMixes; n++)
		{
			m_vmdCovar[nMixIndex].SetSize(m_nVecSize, m_nVecSize);
			m_vmdInvCovar[nMixIndex].SetSize(m_nVecSize, m_nVecSize);
			m_vmdCovar[nMixIndex].Null();
			m_vmdInvCovar[nMixIndex].Null();
			for(int d = 0; d < m_nVecSize; d++)
			{
				m_vmdCovar[nMixIndex](d,d) = 1.0;
				m_vmdInvCovar[nMixIndex](d,d) = 1.0;
			}
		}
	}
	void CompGconsts()
	{
		m_vdgconst.resize(m_nMixes, 0.0);
		m_vvddiaggconst.resize(m_nMixes, vector<double>(m_nVecSize, 0.0));
		for(int n = 0; n < m_nMixes; n++)
		{
			//cout<<n<<"-th mixture, "<<m_vmdCovar[n].ColNo()<<"x"<<m_vmdCovar[n].RowNo()<<endl;
			if(m_vMixWeights[n] > 0)
			{
				size_t covtype = m_szCovType.compare("diagonal");
				if(covtype)
				{
					m_vdgconst[n] = (double)m_nVecSize * log(TWOPI) + log(m_vmdCovar[n].Det());
				}
				else
				{
					m_vdgconst[n] = (double)m_nVecSize * log(TWOPI);
					for(int d = 0; d < m_nVecSize; d++)
					{
						m_vdgconst[n] += log(m_vmdCovar[n](d,d));
					}
				}
				for(int d = 0; d < m_nVecSize; d++)
				{
					m_vvddiaggconst[n][d] = log(TWOPI) + log(m_vmdCovar[n](d,d));
				}			
			}
			else
			{
				m_vdgconst[n] = LZERO;			     
				for(int d = 0; d < m_nVecSize; d++)
				{
					m_vvddiaggconst[n][d] = LZERO;
				}
			}
		}
	}
	// Actually, return log likelihood instead of likelihood
	double CompLikelihood(const vector<double>& vec)
	{
		matrix<double> mvec;
		size_t nsize = vec.size();
		mvec.SetSize(nsize, 1);
		for(size_t n = 0; n < nsize; n++)
		{
			mvec(n, 0) = vec[n];
		}
		vector<double> vdTemp(m_nMixes, 0.0);
		double dTemp = LZERO;
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				matrix<double> x = (mvec - m_vmdMean[n]);
				if(m_szCovType == "diagonal")
				{
#ifndef DIAGONAL			
					for(size_t j = 0; j < m_nVecSize; j++)
					{
						for(size_t i = j+1; i < m_nVecSize; i++)
						{
							vdTemp[n] += x(i, 0) * x(j, 0) * m_vmdInvCovar[n](i,j);
						}
					}
					vdTemp[n] *= 2;
					for(size_t i = 0; i < m_nVecSize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
					}
					vdTemp[n] += m_vdgconst[n];
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#else
					vdTemp[n] += m_vdgconst[n];
					for(size_t i = 0; i < m_nVecSize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
					}
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#endif
				}
				else if(m_szCovType == "full")
				{
					vdTemp[n] += m_vdgconst[n];
					matrix<double> mvectemp = ~x * m_vmdInvCovar[n] * x;
					vdTemp[n] += mvectemp(0, 0);
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
				}
				dTemp = LogAdd(dTemp, vdTemp[n]);
			}
		}
		return dTemp; //exp(dTemp)
	}
	double CompLikelihood(const matrix<double>& mvec)
	{
		vector<double> vdTemp(m_nMixes, 0.0);
		size_t nsize = mvec.ColNo();
		double dTemp = LZERO;
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				matrix<double> x = (mvec - m_vmdMean[n]);
				if(m_szCovType == "diagonal")
				{
#ifndef DIAGONAL			
					for(size_t j = 0; j < m_nVecSize; j++)
					{
						for(size_t i = j+1; i < m_nVecSize; i++)
						{
							vdTemp[n] += x(i, 0) * x(j, 0) * m_vmdInvCovar[n](i,j);
						}
					}
					vdTemp[n] *= 2;
					for(size_t i = 0; i < m_nVecSize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
					}
					vdTemp[n] += m_vdgconst[n];
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#else
					vdTemp[n] += m_vdgconst[n];
					for(size_t i = 0; i < m_nVecSize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
					}
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#endif
				}
				else if(m_szCovType == "full")
				{
					vdTemp[n] += m_vdgconst[n];
					vdTemp[n] += (~x * m_vmdInvCovar[n] * x)(0, 0);
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
				}
				dTemp = LogAdd(dTemp, vdTemp[n]);
			}		
		}
		return dTemp; //exp(dTemp)
	}
	// the posterior probability ( not log )
	void CompPosteriori(const vector<double>& vec, vector<double>& vPost)
	{
		matrix<double> mvec;
		size_t nsize = vec.size();
		mvec.SetSize(nsize, 1);
		for(size_t n = 0; n < nsize; n++)
		{
			mvec(n, 0) = vec[n];
		}
		vector<double> vdTemp(m_nMixes, 0.0);
		double dTemp = LZERO;
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				matrix<double> x = (mvec - m_vmdMean[n]);
				if(m_szCovType == "diagonal")
				{
#ifndef DIAGONAL			
					for(size_t j = 0; j < m_nVecSize; j++)
					{
						for(size_t i = j+1; i < m_nVecSize; i++)
						{
							vdTemp[n] += x(i, 0) * x(j, 0) * m_vmdInvCovar[n](i,j);
						}
					}
					vdTemp[n] *= 2;
					for(size_t i = 0; i < m_nVecSize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
					}
					vdTemp[n] += m_vdgconst[n];
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#else
					vdTemp[n] += m_vdgconst[n];
					for(size_t i = 0; i < m_nVecSize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](i,i);
					}
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#endif
				}
				else if(m_szCovType == "full")
				{
					vdTemp[n] += m_vdgconst[n];
					vdTemp[n] += (~x * m_vmdInvCovar[n] * x)(0, 0);
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
				}
				dTemp = LogAdd(dTemp, vdTemp[n]);
			}
		}
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				vdTemp[n] -= dTemp;
				vdTemp[n] = exp(vdTemp[n]);
			}
			else
			{
				vdTemp[n] = 0.0;
			}
		}
		vPost = vdTemp;
	}
	// Get mean vector for subspace, the space chosen is indicated by vDimConsidered
	void GetPartialMeans(vector< matrix<double> >& PartialMeans, vector<int>& vDimConsidered)
	{
		PartialMeans.resize(m_nMixes);
		for(size_t k = 0; k < m_vmdMean.size(); k++)
		{
			PartialMeans[k].SetSize(vDimConsidered.size(), 1);
			for(size_t n = 0; n < vDimConsidered.size(); n++)
			{
				PartialMeans[k](n, 0) = m_vmdMean[k](vDimConsidered[n], 0);
			}
		}
	}
	// Get covariance matrix for subspace, the space chosen is indicated by vDimConsidered
	void GetPartialVars(vector< matrix<double> >& PartialVars, vector< matrix<double> >& PartialInvVars, vector<double>& vPartialGconst, vector<int>& vDimConsidered)
	{
		PartialVars.resize(m_nMixes);
		PartialInvVars.resize(m_nMixes);
		for(size_t k = 0; k < m_vmdMean.size(); k++)
		{
			PartialVars[k].SetSize(vDimConsidered.size(), vDimConsidered.size());
			PartialInvVars[k].SetSize(vDimConsidered.size(), vDimConsidered.size());
			for(size_t n = 0; n < vDimConsidered.size(); n++)
			{
				for(size_t m = 0; m < vDimConsidered.size(); m++)
				{
					PartialVars[k](n, m) = m_vmdCovar[k](vDimConsidered[n], vDimConsidered[m]);
					PartialInvVars[k](n, m) = m_vmdInvCovar[k](vDimConsidered[n], vDimConsidered[m]);
				}
			}
		}
		vPartialGconst.resize(m_nMixes, 0.0);
		for(int n = 0; n < m_nMixes; n++)
		{
			vPartialGconst[n] = (double)vDimConsidered.size() * log(TWOPI)+log(PartialVars[n].Det());
		}
	}	
    // obsolete function, don't use it
	void CompPosterioriPartOld(const vector<double>& vec, vector<double>& vPost, vector<int>& vDimConsidered)
	{
		matrix<double> mvec;
		vector< matrix<double> > PartialMeans;
		vector< matrix<double> > PartialVars;
		vector< matrix<double> > PartialInvVars;
		vector<double> vPartialGconst;

		size_t nsize = vDimConsidered.size();
		mvec.SetSize(nsize, 1);
		for(size_t n = 0; n < nsize; n++)
		{
			mvec(n, 0) = vec[vDimConsidered[n]];
		}
		
		
		GetPartialMeans(PartialMeans, vDimConsidered);
		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
		
		vector<double> vdTemp(m_nMixes, 0.0);
		double dTemp = 0.0;
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				matrix<double> x = (mvec - PartialMeans[n]);
				matrix<double> xt = ~x;
				matrix<double> mval = 0.5 * xt * PartialInvVars[n] * x;
				vdTemp[n] = -0.5 * vPartialGconst[n] - mval(0,0);
				vdTemp[n] = m_vMixWeights[n] * exp(vdTemp[n]);
				dTemp += vdTemp[n];
			}
		}
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				vdTemp[n] /= dTemp;
			}
		}
		vPost = vdTemp;
	}

	// Compute posterior probability from subspace
	void CompPosterioriPart(const vector<double>& vec, vector<double>& vPost, vector<int>& vDimConsidered)
	{
		matrix<double> mvec;
		vector< matrix<double> > PartialMeans;
		vector< matrix<double> > PartialVars;
		vector< matrix<double> > PartialInvVars;
		vector<double> vPartialGconst;

		size_t nsize = vDimConsidered.size();

		// mvec.SetSize(nsize, 1);
		// for(size_t n = 0; n < nsize; n++)
		// {
			// mvec(n, 0) = vec[vDimConsidered[n]];
		// }
		vector<double> vdTemp(m_nMixes, 0.0);
		double dTemp = LZERO;
#ifndef DIAGONAL		
		// mvec.SetSize(nsize, 1);
		// for(size_t n = 0; n < nsize; n++)
		// {
			// mvec(n, 0) = mvec0(vDimConsidered[n], 0);
		// }
		//dw2 = timeGetTime();
		//cout<<"CompLikelihoodPart: Prepare vector takes "<<dw2-dw1<<" ms"<<endl;
		//dw1 = timeGetTime();
		//GetPartialMeans(PartialMeans, vDimConsidered);
		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
#endif
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				matrix<double> x(nsize, 1);
				for(size_t d = 0; d < nsize; d++)
					x(d,0) = (vec[vDimConsidered[d]] - m_vmdMean[n](vDimConsidered[d], 0));
				if(m_szCovType == "diagonal")
				{
#ifndef DIAGONAL			
				for(size_t j = 0; j < nsize; j++)
				{
					for(size_t i = j+1; i < nsize; i++)
					{
						vdTemp[n] += x(i, 0) * x(j, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[j]);
					}
				}
				vdTemp[n] *= 2;
				for(size_t i = 0; i < nsize; i++)
				{
					vdTemp[n] += x(i, 0) * x(i, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[i]);
				}
				vdTemp[n] += vPartialGconst[n];
				vdTemp[n] *= -0.5;
				vdTemp[n] += m_vlogMixWeights[n];
#else
				for(size_t d = 0; d < nsize; d++)
					vdTemp[n] += m_vvddiaggconst[n][vDimConsidered[d]];
				for(size_t i = 0; i < nsize; i++)
				{
					vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](vDimConsidered[i],vDimConsidered[i]);
				}
				vdTemp[n] *= -0.5;
				vdTemp[n] += m_vlogMixWeights[n];
#endif
				}
				else if(m_szCovType == "full")
				{
					vdTemp[n] += vPartialGconst[n];
					vdTemp[n] += (~x * PartialInvVars[n] * x)(0, 0);
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];					
				}
				dTemp = LogAdd(dTemp, vdTemp[n]);
			}
		}
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				vdTemp[n] -= dTemp;
				vdTemp[n] = exp(vdTemp[n]);
			}
		}
		vPost = vdTemp;
	}
	// Obsolete function, don't use it
	double CompLikelihoodPartOld(const vector<double>& vec, vector<int>& vDimConsidered)
	{
		DWORD dw1 = 0, dw2 = 0;
		//dw1 = timeGetTime();
		matrix<double> mvec;
		vector< matrix<double> > PartialMeans;
		vector< matrix<double> > PartialVars;
		vector< matrix<double> > PartialInvVars;
		vector<double> vPartialGconst;
		
		size_t nsize = vDimConsidered.size();
		mvec.SetSize(nsize, 1);
		for(size_t n = 0; n < nsize; n++)
		{
			mvec(n, 0) = vec[vDimConsidered[n]];
		}
		
		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
		GetPartialMeans(PartialMeans, vDimConsidered);
		vector<double> vdTemp(m_nMixes, 0.0);
		double dTemp = 0;
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				matrix<double> x = mvec - PartialMeans[n];
				matrix<double> xt = ~x;
				matrix<double> mval = 0.5 * xt * PartialInvVars[n] * x;
				vdTemp[n] = -0.5 * vPartialGconst[n] - mval(0,0);
				vdTemp[n] = m_vMixWeights[n] * exp(vdTemp[n]);
				dTemp += vdTemp[n];
			}
		}
		return log(dTemp); //exp(dTemp);
	}
	// Compuete log likelihood from subspace
	double CompLikelihoodPart(const vector<double>& vec, vector<int>& vDimConsidered)
	{
		DWORD dw1 = 0, dw2 = 0;
		//dw1 = timeGetTime();
		matrix<double> mvec;
		vector< matrix<double> > PartialMeans;
		vector< matrix<double> > PartialVars;
		vector< matrix<double> > PartialInvVars;
		vector<double> vPartialGconst;
		
		size_t nsize = vDimConsidered.size();
#ifndef DIAGONAL		
		// mvec.SetSize(nsize, 1);
		// for(size_t n = 0; n < nsize; n++)
		// {
			// mvec(n, 0) = vec[vDimConsidered[n]];
		// }
		//dw2 = timeGetTime();
		//cout<<"CompLikelihoodPart: Prepare vector takes "<<dw2-dw1<<" ms"<<endl;
		//dw1 = timeGetTime();
		//GetPartialMeans(PartialMeans, vDimConsidered);
		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
#endif		
		vector<double> vdTemp(m_nMixes, 0.0);
		double dTemp = LZERO;
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)	
			{
				matrix<double> x(nsize, 1);
				for(size_t d = 0; d < nsize; d++)
					x(d,0) = (vec[vDimConsidered[d]] - m_vmdMean[n](vDimConsidered[d], 0));
				if(m_szCovType == "diagonal")
				{
#ifndef DIAGONAL			
					for(size_t j = 0; j < nsize; j++)
					{
						for(size_t i = j+1; i < nsize; i++)
						{
							vdTemp[n] += x(i, 0) * x(j, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[j]);
						}
					}
					vdTemp[n] *= 2;
					for(size_t i = 0; i < nsize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[i]);
					}
					vdTemp[n] += vPartialGconst[n];
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#else
					for(size_t d = 0; d < nsize; d++)
						vdTemp[n] += m_vvddiaggconst[n][vDimConsidered[d]];
					for(size_t i = 0; i < nsize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](vDimConsidered[i],vDimConsidered[i]);
					}
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#endif
				}
				else if(m_szCovType == "full")
				{
					vdTemp[n] += vPartialGconst[n];
					vdTemp[n] += (~x * PartialInvVars[n] * x)(0, 0);
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];					
				}
				dTemp = LogAdd(dTemp, vdTemp[n]);
			}
		}
		return dTemp; //exp(dTemp);
	}
	// Compute log likelihood from subspace
	double CompLikelihoodPart(const matrix<double>& mvec0, vector<int>& vDimConsidered)
	{
		matrix<double> mvec;
		vector< matrix<double> > PartialMeans;
		vector< matrix<double> > PartialVars;
		vector< matrix<double> > PartialInvVars;
		vector<double> vPartialGconst;

		size_t nsize = vDimConsidered.size();
		
#ifndef DIAGONAL		
		// mvec.SetSize(nsize, 1);
		// for(size_t n = 0; n < nsize; n++)
		// {
			// mvec(n, 0) = mvec0(vDimConsidered[n], 0);
		// }
		//dw2 = timeGetTime();
		//cout<<"CompLikelihoodPart: Prepare vector takes "<<dw2-dw1<<" ms"<<endl;
		//dw1 = timeGetTime();
		//GetPartialMeans(PartialMeans, vDimConsidered);
		GetPartialVars(PartialVars, PartialInvVars, vPartialGconst, vDimConsidered);
#endif
		vector<double> vdTemp(m_nMixes, 0.0);
		double dTemp = LZERO;
		for(size_t n = 0; n < m_nMixes; n++)
		{
			if(m_vMixWeights[n] > 0)
			{
				matrix<double> x(nsize, 1);
				for(size_t d = 0; d < nsize; d++)
					x(d,0) = (mvec0(vDimConsidered[d], 0) - m_vmdMean[n](vDimConsidered[d], 0));
				if(m_szCovType == "diagonal")
				{
#ifndef DIAGONAL			
					for(size_t j = 0; j < nsize; j++)
					{
						for(size_t i = j+1; i < nsize; i++)
						{
							vdTemp[n] += x(i, 0) * x(j, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[j]);
						}
					}
					vdTemp[n] *= 2;
					for(size_t i = 0; i < nsize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * PartialInvVars(vDimConsidered[i],vDimConsidered[i]);
					}
					vdTemp[n] += vPartialGconst[n];
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#else
					for(size_t d = 0; d < nsize; d++)
						vdTemp[n] += m_vvddiaggconst[n][vDimConsidered[d]];
					for(size_t i = 0; i < nsize; i++)
					{
						vdTemp[n] += x(i, 0) * x(i, 0) * m_vmdInvCovar[n](vDimConsidered[i],vDimConsidered[i]);
					}
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];
#endif
				}
				else if(m_szCovType == "full")
				{
					vdTemp[n] += vPartialGconst[n];
					vdTemp[n] += (~x * PartialInvVars[n] * x)(0, 0);
					vdTemp[n] *= -0.5;
					vdTemp[n] += m_vlogMixWeights[n];					
				}				
				dTemp = LogAdd(dTemp, vdTemp[n]);
			}
		}
		return dTemp;
	}
	//void InitHyperParams()
	//{
	//	m_theta.InitHyperParambyGMM(&this);
	//}
	int m_nVecSize;
	int m_nMixes;
	string m_szCovType;
	string m_szFeatureType;
	vector<double> m_vMixWeights;
	vector< matrix<double> > m_vmdMean;
	vector< matrix<double> > m_vmdCovar;
	vector< matrix<double> > m_vmdInvCovar;
	vector<double> m_vdgconst;
	vector<double> m_vlogMixWeights;
	vector< vector<double> > m_vvddiaggconst;
	vector< vector<double> > m_vvTransp;
	//GMMHyperParam m_hyperparam;
	//GMMHyperParam* m_ptheta;
};

