#pragma once
#include "GMMAdaptTool.h"

#define MYNULLINDEX 0x80000000

class BlockFeatures
{
public:
	static void cvmat2vec(cv::Mat& mat, vector< vector<double> >& vvData)
	{
		int channels = mat.channels();
		int ndepth = mat.depth();
		vector<Mat> planes;
		cv::split(mat, planes);
		vvData.resize(channels, vector<double>(mat.rows * mat.cols, 0.0));
				
		for(int c = 0; c < channels; c++)
		{
			int ncount = 0;
			for(int n = 0; n < mat.rows; n++)
			{
				if(ndepth == CV_8U)
				{
					uchar* pfirst = planes[c].ptr<uchar>(n);
					for(int m = 0; m < mat.cols; m++)
					{
						vvData[c][ncount] = (double)pfirst[m];
						ncount++;
					}
				}
				else if(ndepth == CV_8S)
				{
					char* pfirst = planes[c].ptr<char>(n);
					for(int m = 0; m < mat.cols; m++)
					{
						vvData[c][ncount] = (double)pfirst[m];
						ncount++;
					}
				}
				else if(ndepth == CV_16U)
				{
					ushort* pfirst = planes[c].ptr<ushort>(n);
					for(int m = 0; m < mat.cols; m++)
					{
						vvData[c][ncount] = (double)pfirst[m];
						ncount++;
					}
				}
				else if(ndepth == CV_16S)
				{
					short* pfirst = planes[c].ptr<short>(n);
					for(int m = 0; m < mat.cols; m++)
					{
						vvData[c][ncount] = (double)pfirst[m];
						ncount++;
					}
				}
				else if(ndepth == CV_32S)
				{
					int* pfirst = planes[c].ptr<int>(n);
					for(int m = 0; m < mat.cols; m++)
					{
						vvData[c][ncount] = (double)pfirst[m];
						ncount++;
					}
				}
				else if(ndepth == CV_32F)
				{
					float* pfirst = planes[c].ptr<float>(n);
					for(int m = 0; m < mat.cols; m++)
					{
						vvData[c][ncount] = (double)pfirst[m];
						ncount++;
					}
				}
				else if(ndepth == CV_64F)
				{
					double* pfirst = planes[c].ptr<double>(n);
					for(int m = 0; m < mat.cols; m++)
					{
						vvData[c][ncount] = (double)pfirst[m];
						ncount++;
					}
				}
			}
		}
	}

	void setmask(cv::Mat& mat)
	{
		int channels = mat.channels();
		int ndepth = mat.depth();
		if(channels != 1)
		{
			cerr<<"Number of channels must be 1 for a mask"<<endl;
			return;
		}
		if(ndepth != CV_8U)
		{
			cerr<<"Data type must be CV_8U for a mask"<<endl;
			return;		
		}
		vector<Mat> planes;
		cv::split(mat, planes);
		m_vMask = vector<double>(mat.rows * mat.cols, 0.0);
				
		for(int c = 0; c < channels; c++)
		{
			int ncount = 0;
			for(int n = 0; n < mat.rows; n++)
			{
				uchar* pfirst = planes[c].ptr<uchar>(n);
				for(int m = 0; m < mat.cols; m++)
				{
					uchar myval = mat.at<uchar>(n, m);
					m_vMask[ncount] = (double)myval;
					ncount++;
				}
			}
		}
	}


	static void compquantilesrangemode(vector<double>& vData, double& quad25, double& median, double& quad75, double& range, double& mode)
	{
		int N = (int)vData.size();
		if(N <= 0)
		{
			cerr<<"Empty data, cannot compute quantiles"<<endl;
			return;
		}
		vector<int> dummy(N, 0);
		for(size_t n = 0; n < N; n++)
		{
			dummy[n] = n;
		}
		sort(dummy.begin(), dummy.end(), index_cmp_ascend< vector<double>& >(vData));
		int q25n1 = (int)floor((double)(N + 1) * 0.25);
		int q25n2 = (int)ceil((double)(N + 1) * 0.25);
		int q50n1 = (int)floor((double)(N + 1) * 0.50);
		int q50n2 = (int)ceil((double)(N + 1) * 0.50);
		int q75n1 = (int)floor((double)(N + 1) * 0.75);
		int q75n2 = (int)ceil((double)(N + 1) * 0.75);
		quad25 = (vData[dummy[q25n1]] + vData[dummy[q25n2]])/2.0;
		median = (vData[dummy[q50n1]] + vData[dummy[q50n2]])/2.0;
		quad75 = (vData[dummy[q75n1]] + vData[dummy[q75n2]])/2.0;
		range = vData[dummy[N-1]] - vData[dummy[0]];
		// compute mode
		unordered_map<int, int> hCounts;
		for(size_t n = 0; n < N; ++n)
		{
			hCounts[(int)vData[n]]++;
		}
		int nmaxcount = 0;
		for(unordered_map<int, int>::iterator itr = hCounts.begin(); itr != hCounts.end(); ++itr)
		{
			int ncountnow = itr->second;
			if(nmaxcount < ncountnow)
			{
				nmaxcount = ncountnow;
				mode = itr->first;
			}
		}
	}


	void CompQuantiles_n_Range_Mode()
	{
		size_t C = m_vvData.size();
		m_vQuad25 = vector<double>(3, 0.0);
		m_vMedian = vector<double>(3, 0.0);
		m_vQuad75 = vector<double>(3, 0.0);
		m_vMode = vector<double>(3, 0.0);
		m_vRange = vector<double>(3, 0.0);
		for(size_t c = 0; c < C; c++)
		{
			compquantilesrangemode(m_vvData[c], m_vQuad25[c], m_vMedian[c], m_vQuad75[c], m_vRange[c], m_vMode[c]); 
		}
	}
	void extractBlockMean()
	{
		size_t C = m_vvData.size();
		m_vMeans = vector<double>(C, 0.0);
		if(m_vMask.size() > 0)
		{
			for(size_t c = 0; c < C; ++c)
			{
				size_t Nc = m_vvData[c].size();
				if(m_vMask.size() != Nc)
				{
					cerr<<"Mask size did not match the block, ignored"<<endl;
					for(size_t n = 0; n < Nc; ++n)
					{
						m_vMeans[c] += m_vvData[c][n];
					}
					m_vMeans[c] /= (double)Nc;
				}
				else
				{
					double dbNc = 0;
					for(size_t n = 0; n < Nc; ++n)
					{
						if(m_vMask[n] == 255)
						{
							m_vMeans[c] += m_vvData[c][n];
							dbNc += 1.0;
						}
					}
					m_vMeans[c] /= dbNc;

				}
			}
		}
		else
		{
			for(size_t c = 0; c < C; ++c)
			{
				size_t Nc = m_vvData[c].size();
				for(size_t n = 0; n < Nc; ++n)
				{
					m_vMeans[c] += m_vvData[c][n];
				}
				m_vMeans[c] /= (double)Nc;
			}
		}
	}

	void extractBlockVars()
	{
		size_t C = m_vvData.size();
		m_vVars = vector<double>(C, 0.0);
		if(m_vMask.size() > 0)
		{
			for(size_t c = 0; c < C; ++c)
			{
				size_t Nc = m_vvData[c].size();
				if(m_vMask.size() != Nc)
				{
					cerr<<"Mask size did not match the block, ignored"<<endl;
					for(size_t n = 0; n < Nc; ++n)
					{
						m_vVars[c] += pow((m_vvData[c][n] - m_vMeans[c]), 2.0);
					}
					m_vVars[c] /= (double)((Nc > 1)?(Nc - 1):Nc);
				}
				else
				{
					double dbNc = 0;
					for(size_t n = 0; n < Nc; ++n)
					{
						if(m_vMask[n] == 255)
						{
							m_vVars[c] += pow((m_vvData[c][n] - m_vMeans[c]), 2.0);
							dbNc += 1.0;
						}
					}
					m_vVars[c] /= (double)((dbNc > 1.0)?(dbNc - 1):dbNc);

				}
			}
		}
		else
		{
			for(size_t c = 0; c < C; ++c)
			{
				size_t Nc = m_vvData[c].size();
				for(size_t n = 0; n < Nc; ++n)
				{
					m_vVars[c] += pow((m_vvData[c][n] - m_vMeans[c]), 2.0);
				}
				m_vVars[c] /= (double)((Nc > 1)?(Nc - 1):Nc);
			}
		}
	}
	
	void compCorrelation(BlockFeatures& rhs, vector<double>& vCorr)
	{
		if(rhs.m_vvData.size() != m_vvData.size())
		{
			cerr<<"number of channel mismatched in compCorrelation"<<endl;
			return;
		}

		for(size_t n = 0; n < m_vvData.size(); n++)
		{
			if(rhs.m_vvData[n].size() != m_vvData[n].size())
			{
				cerr<<"number of data mismatched in channel "<<n<<" in compCorrelation"<<endl;
				return;
			}
		}
		bool bWillApplyMask = true;
		for(size_t n = 0; n < m_vvData.size(); ++n)
		{
			if(m_vvData[n].size() != m_vMask.size())
				bWillApplyMask = false;
		}
		for(size_t n = 0; n < rhs.m_vvData.size(); ++n)
		{
			if(rhs.m_vvData[n].size() != rhs.m_vMask.size())
				bWillApplyMask = false;
		}
		if(rhs.m_vMask.size() != m_vMask.size())
		{
			bWillApplyMask = false;
			cerr<<"mask sizes mismatch, ignored\n";
		}
		if(bWillApplyMask)
		{
			for(size_t n = 0; n < m_vMask.size(); ++n)
			{
				if(m_vMask[n] != rhs.m_vMask[n])
				{
					bWillApplyMask = false;
					cerr<<"different masks used, ignored\n";
					break;
				}
			}
		}

		size_t C = m_vvData.size();
		vCorr.reserve(C);
		vCorr.resize(C);
		for(size_t c = 0; c < C; ++c)
		{
			if(bWillApplyMask)
			{
				size_t Nc = m_vvData[c].size();
				double dbVal = 0.0;
				double dbNc = 0.0;
				for(size_t i = 0; i < Nc; ++i)
				{
					if(m_vMask[i] == 255)
					{
						dbVal += (m_vvData[c][i] - m_vMeans[c]) * (rhs.m_vvData[c][i] - rhs.m_vMeans[c]);
						dbNc += 1.0;
					}
				}
				double norm = (double)((dbNc > 1)?(dbNc-1):dbNc);
				double trnorm = 10.0;
				double dbVar1 = sqrt(m_vVars[c] * norm);
				double dbVar2 = sqrt(rhs.m_vVars[c] * norm);
				double dbVar = dbVar1 * dbVar2;
				if(dbVar1 < trnorm && dbVar2 < trnorm) // constant to constant map
				{
					vCorr[c] = 1.0;
				}
				else
				{
					vCorr[c] = fabs((dbVar > 0)?(dbVal / dbVar):0.0);
				}

			}
			else
			{
				size_t Nc = m_vvData[c].size();
				double dbVal = 0.0;
				for(size_t i = 0; i < Nc; ++i)
				{
					dbVal += (m_vvData[c][i] - m_vMeans[c]) * (rhs.m_vvData[c][i] - rhs.m_vMeans[c]);
				}
				double norm = (double)((Nc > 1)?(Nc-1):Nc);
				double trnorm = 10.0;
				double dbVar1 = sqrt(m_vVars[c] * norm);
				double dbVar2 = sqrt(rhs.m_vVars[c] * norm);
				double dbVar = dbVar1 * dbVar2;
				if(dbVar1 < trnorm && dbVar2 < trnorm) // constant to constant map
				{
					vCorr[c] = 1.0;
				}
				else
				{
					vCorr[c] = fabs((dbVar > 0)?(dbVal / dbVar):0.0);
				}
			}
		}
	}

	void compDiffCorrelation(BlockFeatures& rhs)
	{
		compCorrelation(rhs, m_vDiffCorr);
	}

	void compBSCorrelation(BlockFeatures& rhs)
	{
		compCorrelation(rhs, m_vBSCorr);
	}

	//;

	void CompDeltaQuantiles_n_Range_Mode()
	{
		size_t C = m_vvDeltaData.size();
		m_vDeltaQuad25 = vector<double>(3, 0.0);
		m_vDeltaMedian = vector<double>(3, 0.0);
		m_vDeltaQuad75 = vector<double>(3, 0.0);
		m_vDeltaMode = vector<double>(3, 0.0);
		m_vDeltaRange = vector<double>(3, 0.0);
		for(size_t c = 0; c < C; c++)
		{
			compquantilesrangemode(m_vvDeltaData[c], m_vDeltaQuad25[c], m_vDeltaMedian[c], m_vDeltaQuad75[c], m_vDeltaRange[c], m_vDeltaMode[c]); 
		}
	}
	void extractDeltaBlockMean()
	{
		size_t C = m_vvDeltaData.size();
		m_vDeltaMeans = vector<double>(C, 0.0);
		for(size_t c = 0; c < C; ++c)
		{
			size_t Nc = m_vvDeltaData[c].size();
			for(size_t n = 0; n < Nc; ++n)
			{
				m_vDeltaMeans[c] += m_vvDeltaData[c][n];
			}
			m_vDeltaMeans[c] /= (double)Nc;
		}
	}

	void extractDeltaBlockVars()
	{
		size_t C = m_vvDeltaData.size();
		m_vDeltaVars = vector<double>(C, 0.0);
		for(size_t c = 0; c < C; ++c)
		{
			size_t Nc = m_vvDeltaData[c].size();
			for(size_t n = 0; n < Nc; ++n)
			{
				m_vDeltaVars[c] += pow((m_vvDeltaData[c][n] - m_vDeltaMeans[c]), 2.0);
			}
			m_vDeltaVars[c] /= (double)((Nc > 1)?(Nc - 1):Nc);
		}
	}

	void CompBSQuantiles_n_Range_Mode()
	{
		size_t C = m_vvBSData.size();
		m_vBSQuad25 = vector<double>(3, 0.0);
		m_vBSMedian = vector<double>(3, 0.0);
		m_vBSQuad75 = vector<double>(3, 0.0);
		m_vBSMode = vector<double>(3, 0.0);
		m_vBSRange = vector<double>(3, 0.0);
		for(size_t c = 0; c < C; c++)
		{
			compquantilesrangemode(m_vvBSData[c], m_vBSQuad25[c], m_vBSMedian[c], m_vBSQuad75[c], m_vBSRange[c], m_vBSMode[c]); 
		}
	}
	void extractBSBlockMean()
	{
		size_t C = m_vvBSData.size();
		m_vBSMeans = vector<double>(C, 0.0);
		for(size_t c = 0; c < C; ++c)
		{
			size_t Nc = m_vvBSData[c].size();
			for(size_t n = 0; n < Nc; ++n)
			{
				m_vBSMeans[c] += m_vvBSData[c][n];
			}
			m_vBSMeans[c] /= (double)Nc;
		}
	}

	void extractBSBlockVars()
	{
		size_t C = m_vvBSData.size();
		m_vBSVars = vector<double>(C, 0.0);
		for(size_t c = 0; c < C; ++c)
		{
			size_t Nc = m_vvBSData[c].size();
			for(size_t n = 0; n < Nc; ++n)
			{
				m_vBSVars[c] += pow((m_vvBSData[c][n] - m_vBSMeans[c]), 2.0);
			}
			m_vBSVars[c] /= (double)((Nc > 1)?(Nc - 1):Nc);
		}
	}

	void extractBlock_V_CDF()
	{
		size_t C = m_vvData.size();
		m_vVCDF = vector<double>(256, 0.0);
		if(C < 1) {
			cerr<<"Did not have color channel!!!\n";
			return;
		}
		size_t Nc = m_vvData[0].size();
		if(m_vMask.size() > 0)
		{
			if(m_vMask.size() != Nc)
			{
				cerr<<"Mask size did not match the block, ignored"<<endl;
				for(size_t n = 0; n < Nc; ++n)
				{
					double v = 0.0;
					for(size_t c = 0; c < C; ++c) {
						if(m_vvData[c][n] > v)
							v = m_vvData[c][n];
					}
					m_vVCDF[(int)floor(v+0.5)] += 1.0;
				}
				for(size_t i = 0; i < 256; ++i)
				{
					m_vVCDF[i] /= (double)Nc;
					if(i >= 1)
						m_vVCDF[i] += m_vVCDF[i-1];
				}
			}
			else
			{
				double dbNc = 0;
				for(size_t n = 0; n < Nc; ++n)
				{
					if(m_vMask[n] == 255)
					{
						double v = 0.0;
						for(size_t c = 0; c < C; ++c) {
							if(m_vvData[c][n] > v)
								v = m_vvData[c][n];
						}
						m_vVCDF[(int)floor(v+0.5)] += 1.0;
						dbNc += 1.0;
					}
				}
				for(size_t i = 0; i < 256; ++i)
				{
					m_vVCDF[i] /= (double)Nc;
					if(i >= 1)
						m_vVCDF[i] += m_vVCDF[i-1];
				}
			}
		}
		else
		{
			for(size_t n = 0; n < Nc; ++n)
			{
				double v = 0.0;
				for(size_t c = 0; c < C; ++c) {
					if(m_vvData[c][n] > v)
						v = m_vvData[c][n];
				}
				m_vVCDF[(int)floor(v+0.5)] += 1.0;
			}
			for(size_t i = 0; i < 256; ++i)
			{
				m_vVCDF[i] /= (double)Nc;
				if(i >= 1)
					m_vVCDF[i] += m_vVCDF[i-1];
			}
		}
	}


	int m_nx;
	int m_ny;
	int m_nLabel; 
	int m_nIndex; // -1/-2, forward/backword background, others are purely indices
	vector< vector<double> > m_vvData;
	vector<double> m_vMask;
	// variables to store local features
	vector<double> m_vMeans;
	vector<double> m_vVars;
	vector<double> m_vQuad25;
	vector<double> m_vMedian;
	vector<double> m_vQuad75;
	vector<double> m_vMode;
	vector<double> m_vRange;

	vector< vector<double> > m_vvDeltaData;
	// variables to store local features
	vector<double> m_vDeltaMeans;
	vector<double> m_vDeltaVars;
	vector<double> m_vDeltaQuad25;
	vector<double> m_vDeltaMedian;
	vector<double> m_vDeltaQuad75;
	vector<double> m_vDeltaMode;
	vector<double> m_vDeltaRange;

	vector< vector<double> > m_vvBSData;
	// variables to store local features
	vector<double> m_vBSMeans;
	vector<double> m_vBSVars;
	vector<double> m_vBSQuad25;
	vector<double> m_vBSMedian;
	vector<double> m_vBSQuad75;
	vector<double> m_vBSMode;
	vector<double> m_vBSRange;
	
	vector<double> m_vDiffCorr;
	vector<double> m_vBSCorr;
	vector<double> m_vVCDF;
};

class SingleFrameFeatures
{
public:
	static void cvmat2vec(cv::Mat& mat, vector< vector<double> >& vvData)
	{
		BlockFeatures::cvmat2vec(mat, vvData);
	}

	void init(cv::Mat& mat, JPGLabelInfo* pinfo, int nIndex, int BlockSize = 16)
	{
		cv::Mat fmat;
		int C = mat.channels();
		int newtype = CV_MAKE_TYPE(CV_64F, C);
		mat.convertTo(fmat, newtype);
		m_nBlockSize = BlockSize;
		m_nBlockWidth = mat.cols / m_nBlockSize;
		m_nBlockHeight = mat.rows / m_nBlockSize;
		int nsize = m_nBlockWidth * m_nBlockHeight;
		if(nsize <= 0)
		{
			cerr<<"the image is empty, so cannot extract feature!"<<endl;
			return;
		}
		m_nBlocks = nsize;
		m_vBlocks.reserve(nsize);
		m_vBlocks.resize(nsize);
		int ncount = 0;
		for(int n = 0; n < m_nBlockHeight; ++n)
		{
			cv::Mat fmatn = fmat.rowRange(n*m_nBlockSize, (n+1)*m_nBlockSize);
			for(int m = 0; m < m_nBlockWidth; ++m)
			{
				cv::Mat fmatnm = fmatn.colRange(m*m_nBlockSize, (m+1)*m_nBlockSize);
				m_vBlocks[ncount].m_nx = m;
				m_vBlocks[ncount].m_ny = n;
				m_vBlocks[ncount].m_nLabel = pinfo->m_ppnLabel[n][m];
				m_vBlocks[ncount].m_nIndex = m_nIndex = nIndex;
				cvmat2vec(fmatnm, m_vBlocks[ncount].m_vvData);
				ncount++;
			}
		}
	}
	// should only be called after init
	void ExtractSingleFrameStaticFeatures(cv::Mat& mat, int BlockSize = 16)
	{
		// Extract block mean and variance;
		cv::Scalar mean1, stddev1;
		cv::meanStdDev(mat, mean1, stddev1);
		m_vMeans.push_back(mean1.val[0]);
		m_vMeans.push_back(mean1.val[1]);
		m_vMeans.push_back(mean1.val[2]);
		m_vVars.push_back(stddev1.val[0] * stddev1.val[0]);
		m_vVars.push_back(stddev1.val[1] * stddev1.val[1]);
		m_vVars.push_back(stddev1.val[2] * stddev1.val[2]);
		for(int n = 0; n < m_nBlocks; ++n)
		{
			// Call feature extraction function of BlockFeatures
			m_vBlocks[n].extractBlockMean();
			m_vBlocks[n].extractBlockVars();
			m_vBlocks[n].CompQuantiles_n_Range_Mode();
		}
	}

	string m_szSeqName;
	int m_nBlockSize;
	int m_nBlockWidth;
	int m_nBlockHeight;
	int m_nBlocks;
	int m_nIndex;
	int m_nPrevIndex;
	int m_nNextIndex;


	vector<double> m_vMeans;
	vector<double> m_vVars;
	vector< BlockFeatures > m_vBlocks;

	vector<double> m_vDeltaMeans;
	vector<double> m_vDeltaVars;

	vector<double> m_vBSMeans;
	vector<double> m_vBSVars;
};

class LocalFeatures
{
public:
	LocalFeatures()
	{
		unsigned long nan[2]={0xffffffff, 0x7fffffff};
		m_nan = *((double*)nan);
		m_nCurrIndex = 0;
	}
	static void cvmat2vec(cv::Mat& mat, vector< vector<double> >& vvData)
	{
		BlockFeatures::cvmat2vec(mat, vvData);
	}
	void init(int N)
	{
		m_vFeatures.reserve(N);
		m_vFeatures.resize(N);
		m_nCurrIndex = 0;
	}
	void CompFrameDiff_and_BS_Features()
	{
		int nIndex = m_nCurrIndex;
		int nPrevIndex = nIndex - 1;
		int nBkGndIndex = 0;
		if(nPrevIndex < 0)
		{
			size_t Cn = m_vFeatures[nIndex].m_vMeans.size();
			m_vFeatures[nIndex].m_vDeltaMeans = vector<double>(Cn, m_nan);
			m_vFeatures[nIndex].m_vDeltaVars = vector<double>(Cn, m_nan);
			m_vFeatures[nIndex].m_vBSMeans = vector<double>(Cn, m_nan);
			m_vFeatures[nIndex].m_vBSVars = vector<double>(Cn, m_nan);
			for(size_t nn = 0; nn < m_vFeatures[nIndex].m_nBlocks; ++nn)
			{
				size_t Cnn = m_vFeatures[nIndex].m_vBlocks[nn].m_vvData.size();
				size_t nnSize = (size_t)m_vFeatures[nIndex].m_nBlockSize * m_vFeatures[nIndex].m_nBlockSize;
				m_vFeatures[nIndex].m_vBlocks[nn].m_vvDeltaData = vector< vector<double> >(Cnn, vector<double>(nnSize, m_nan));
				m_vFeatures[nIndex].m_vBlocks[nn].m_vDeltaMeans = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vDeltaVars = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vDeltaQuad25 = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vDeltaMedian = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vDeltaQuad75 = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vDeltaMode = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vDeltaRange = vector<double>(Cnn, m_nan);

				m_vFeatures[nIndex].m_vBlocks[nn].m_vvBSData = vector< vector<double> >(Cnn, vector<double>(nnSize, m_nan));
				m_vFeatures[nIndex].m_vBlocks[nn].m_vBSMeans = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vBSVars = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vBSQuad25 = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vBSMedian = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vBSQuad75 = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vBSMode = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vBSRange = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vDiffCorr = vector<double>(Cnn, m_nan);
				m_vFeatures[nIndex].m_vBlocks[nn].m_vBSCorr = vector<double>(Cnn, m_nan);
			}
		}
		else
		{
			size_t Cn = m_vFeatures[nIndex].m_vMeans.size();
			m_vFeatures[nIndex].m_vDeltaMeans = vector<double>(Cn, 0.0);
			m_vFeatures[nIndex].m_vDeltaVars = vector<double>(Cn, 0.0);
			m_vFeatures[nIndex].m_vBSMeans = vector<double>(Cn, 0.0);
			m_vFeatures[nIndex].m_vBSVars = vector<double>(Cn, 0.0);
			for(size_t cn0 = 0; cn0 < Cn; cn0++)
			{
				m_vFeatures[nIndex].m_vDeltaMeans[cn0] = m_vFeatures[nIndex].m_vMeans[cn0] - 
															m_vFeatures[nPrevIndex].m_vMeans[cn0];
				m_vFeatures[nIndex].m_vDeltaVars[cn0] = m_vFeatures[nIndex].m_vVars[cn0] + 
															m_vFeatures[nPrevIndex].m_vVars[cn0];
				m_vFeatures[nIndex].m_vBSMeans[cn0] = m_vFeatures[nIndex].m_vMeans[cn0] - 
															m_vFeatures[nBkGndIndex].m_vMeans[cn0];
				m_vFeatures[nIndex].m_vBSVars[cn0] = m_vFeatures[nIndex].m_vVars[cn0] + 
															m_vFeatures[nBkGndIndex].m_vVars[cn0];
			}

			for(size_t nn = 0; nn < m_vFeatures[nIndex].m_nBlocks; ++nn)
			{
				size_t Cnn = m_vFeatures[nIndex].m_vBlocks[nn].m_vvData.size();
				size_t nnSize = (size_t)m_vFeatures[nIndex].m_nBlockSize * m_vFeatures[nIndex].m_nBlockSize;
				m_vFeatures[nIndex].m_vBlocks[nn].m_vvDeltaData = vector< vector<double> >(Cnn, vector<double>(nnSize, 0.0));
				m_vFeatures[nIndex].m_vBlocks[nn].m_vvBSData = vector< vector<double> >(Cnn, vector<double>(nnSize, 0.0));
				for(size_t cnnn = 0; cnnn < Cnn; ++cnnn)
				{
					for(size_t nnn = 0; nnn < nnSize; ++nnn)
					{
						m_vFeatures[nIndex].m_vBlocks[nn].m_vvDeltaData[cnnn][nnn] = 
							m_vFeatures[nIndex].m_vBlocks[nn].m_vvData[cnnn][nnn] - 
							m_vFeatures[nPrevIndex].m_vBlocks[nn].m_vvData[cnnn][nnn];
						m_vFeatures[nIndex].m_vBlocks[nn].m_vvBSData[cnnn][nnn] = 
							m_vFeatures[nIndex].m_vBlocks[nn].m_vvData[cnnn][nnn] - 
							m_vFeatures[nBkGndIndex].m_vBlocks[nn].m_vvData[cnnn][nnn];
					}
				}
				m_vFeatures[nIndex].m_vBlocks[nn].extractDeltaBlockMean();
				m_vFeatures[nIndex].m_vBlocks[nn].extractDeltaBlockVars();
				m_vFeatures[nIndex].m_vBlocks[nn].CompDeltaQuantiles_n_Range_Mode();
				m_vFeatures[nIndex].m_vBlocks[nn].extractBSBlockMean();
				m_vFeatures[nIndex].m_vBlocks[nn].extractBSBlockVars();
				m_vFeatures[nIndex].m_vBlocks[nn].CompBSQuantiles_n_Range_Mode();
				m_vFeatures[nIndex].m_vBlocks[nn].compDiffCorrelation(m_vFeatures[nPrevIndex].m_vBlocks[nn]);
				m_vFeatures[nIndex].m_vBlocks[nn].compBSCorrelation(m_vFeatures[nBkGndIndex].m_vBlocks[nn]);
			}
		}
		//int nIndex = m_vFeatures[m_nCurrIndex].m_nIndex;
		//int nPrevIndex = m_vFeatures[m_nCurrIndex].m_nPrevIndex;

	}

	void ExtFeaturesForFrame(int nIndex, cv::Mat& mat, JPGLabelInfo* pinfo, int nBlockSize)
	{
		// Check if this is the first data
		m_vFeatures[m_nCurrIndex].init(mat, pinfo, nIndex, nBlockSize);
		m_vFeatures[m_nCurrIndex].ExtractSingleFrameStaticFeatures(mat, nBlockSize);
		if(nIndex == -1)
		{
			m_vFeatures[m_nCurrIndex].m_nPrevIndex = MYNULLINDEX;
		}
		else if(nIndex == -2)
		{
			m_vFeatures[m_nCurrIndex].m_nPrevIndex = m_vFeatures[m_nCurrIndex - 1].m_nIndex;
			m_vFeatures[m_nCurrIndex].m_nNextIndex = MYNULLINDEX;
		}
		else
		{
			m_vFeatures[m_nCurrIndex].m_nPrevIndex = m_vFeatures[m_nCurrIndex - 1].m_nIndex;
			m_vFeatures[m_nCurrIndex - 1].m_nNextIndex = nIndex;
		}
		// Compute frame diff features and bs features
		CompFrameDiff_and_BS_Features();
		m_nCurrIndex++;
	}
	friend ostream& operator<<(ostream& ostrm, const LocalFeatures& fe)
	{
		size_t N = fe.m_vFeatures.size();
		ostrm<<"Sequence_Name\t"
			<<"Index\t"
			<<"Position_X\t"
			<<"Position_Y\t"
			<<"Label\t"
			<<"Frame_B_Mean\t"
			<<"Frame_G_Mean\t"
			<<"Frame_R Mean\t"
			<<"Frame_B_Variance\t"
			<<"Frame_G_Variance\t"
			<<"Frame_R_Variance\t"
			<<"Frame_Diff_B_Mean\t"
			<<"Frame_Diff_G_Mean\t"
			<<"Frame_Diff_R Mean\t"
			<<"Frame_Diff_B_Variance\t"
			<<"Frame_Diff_G_Variance\t"
			<<"Frame_Diff_R_Variance\t"
			<<"Frame_BkGnd_Sub_B_Mean\t"
			<<"Frame_BkGnd_Sub_G_Mean\t"
			<<"Frame_BkGnd_Sub_R Mean\t"
			<<"Frame_BkGnd_Sub_B_Variance\t"
			<<"Frame_BkGnd_Sub_G_Variance\t"
			<<"Frame_BkGnd_Sub_R_Variance\t"
			<<"Block_B_Mean\t"
			<<"Block_G_Mean\t"
			<<"Block_R_Mean\t"
			<<"Block_B_Median\t"
			<<"Block_G_Median\t"
			<<"Block_R_Meddian\t"
			<<"Block_B_25_Quantile\t"
			<<"Block_G_25_Quantile\t"
			<<"Block_R_25_Quantile\t"
			<<"Block_B_75_Quantile\t"
			<<"Block_G_75_Quantile\t"
			<<"Block_R_75_Quantile\t"
			<<"Block_B_Mode\t"
			<<"Block_G_Mode\t"
			<<"Block_R_Mode\t"
			<<"Block_B_Range\t"
			<<"Block_G_Range\t"
			<<"Block_R_Range\t"
			<<"Block_B_Variance\t"
			<<"Block_G_Variance\t"
			<<"Block_R_Variance\t"
			<<"Block_Diff_B_Mean\t"
			<<"Block_Diff_G_Mean\t"
			<<"Block_Diff_R_Mean\t"
			<<"Block_Diff_B_Median\t"
			<<"Block_Diff_G_Median\t"
			<<"Block_Diff_R_Meddian\t"
			<<"Block_Diff_B_25_Quantile\t"
			<<"Block_Diff_G_25_Quantile\t"
			<<"Block_Diff_R_25_Quantile\t"
			<<"Block_Diff_B_75_Quantile\t"
			<<"Block_Diff_G_75_Quantile\t"
			<<"Block_Diff_R_75_Quantile\t"
			<<"Block_Diff_B_Mode\t"
			<<"Block_Diff_G_Mode\t"
			<<"Block_Diff_R_Mode\t"
			<<"Block_Diff_B_Range\t"
			<<"Block_Diff_G_Range\t"
			<<"Block_Diff_R_Range\t"
			<<"Block_Diff_B_Variance\t"
			<<"Block_Diff_G_Variance\t"
			<<"Block_Diff_R_Variance\t"
			<<"Block_Diff_B_Corr\t"
			<<"Block_Diff_G_Corr\t"
			<<"Block_Diff_R_Corr\t"
			<<"Block_BkGnd_Sub_B_Mean\t"
			<<"Block_BkGnd_Sub_G_Mean\t"
			<<"Block_BkGnd_Sub_R_Mean\t"
			<<"Block_BkGnd_Sub_B_Median\t"
			<<"Block_BkGnd_Sub_G_Median\t"
			<<"Block_BkGnd_Sub_R_Meddian\t"
			<<"Block_BkGnd_Sub_B_25_Quantile\t"
			<<"Block_BkGnd_Sub_G_25_Quantile\t"
			<<"Block_BkGnd_Sub_R_25_Quantile\t"
			<<"Block_BkGnd_Sub_B_75_Quantile\t"
			<<"Block_BkGnd_Sub_G_75_Quantile\t"
			<<"Block_BkGnd_Sub_R_75_Quantile\t"
			<<"Block_BkGnd_Sub_B_Mode\t"
			<<"Block_BkGnd_Sub_G_Mode\t"
			<<"Block_BkGnd_Sub_R_Mode\t"
			<<"Block_BkGnd_Sub_B_Range\t"
			<<"Block_BkGnd_Sub_G_Range\t"
			<<"Block_BkGnd_Sub_R_Range\t"
			<<"Block_BkGnd_Sub_B_Variance\t"
			<<"Block_BkGnd_Sub_G_Variance\t"
			<<"Block_BkGnd_Sub_R_Variance\t"
			<<"Block_BkGnd_Sub_B_Corr\t"
			<<"Block_BkGnd_Sub_G_Corr\t"
			<<"Block_BkGnd_Sub_R_Corr\n";
		for(size_t n = 0; n < N; ++n)
		{
			size_t NB = fe.m_vFeatures[n].m_nBlocks;
			for(size_t m = 0; m < NB; m++)
			{
				ostrm<<fe.m_vFeatures[n].m_szSeqName<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_nIndex<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_nx<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_ny<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_nLabel<<"\t"
					<<fe.m_vFeatures[n].m_vMeans[0]<<"\t"
					<<fe.m_vFeatures[n].m_vMeans[1]<<"\t"
					<<fe.m_vFeatures[n].m_vMeans[2]<<"\t"
					<<fe.m_vFeatures[n].m_vVars[0]<<"\t"
					<<fe.m_vFeatures[n].m_vVars[1]<<"\t"
					<<fe.m_vFeatures[n].m_vVars[2]<<"\t"
					<<fe.m_vFeatures[n].m_vDeltaMeans[0]<<"\t"
					<<fe.m_vFeatures[n].m_vDeltaMeans[1]<<"\t"
					<<fe.m_vFeatures[n].m_vDeltaMeans[2]<<"\t"
					<<fe.m_vFeatures[n].m_vDeltaVars[0]<<"\t"
					<<fe.m_vFeatures[n].m_vDeltaVars[1]<<"\t"
					<<fe.m_vFeatures[n].m_vDeltaVars[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBSMeans[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBSMeans[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBSMeans[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBSVars[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBSVars[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBSVars[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMeans[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMeans[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMeans[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMedian[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMedian[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMedian[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vQuad25[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vQuad25[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vQuad25[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vQuad75[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vQuad75[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vQuad75[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMode[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMode[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vMode[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vRange[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vRange[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vRange[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vVars[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vVars[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vVars[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMeans[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMeans[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMeans[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMedian[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMedian[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMedian[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaQuad25[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaQuad25[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaQuad25[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaQuad75[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaQuad75[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaQuad75[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMode[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMode[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaMode[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaRange[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaRange[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaRange[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaVars[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaVars[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDeltaVars[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDiffCorr[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDiffCorr[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vDiffCorr[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMeans[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMeans[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMeans[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMedian[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMedian[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMedian[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSQuad25[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSQuad25[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSQuad25[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSQuad75[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSQuad75[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSQuad75[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMode[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMode[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSMode[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSRange[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSRange[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSRange[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSVars[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSVars[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSVars[2]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSCorr[0]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSCorr[1]<<"\t"
					<<fe.m_vFeatures[n].m_vBlocks[m].m_vBSCorr[2]<<"\n";
			}
		}
		return ostrm;
	}

	//static void DataSubtraction(vector<double> )
	//void cvmat2vec(cv::Mat& mat, vector< vector<double> >& vvData)
	//{
	//	int channels = mat.channels();
	//	int ndepth = mat.depth();
	//	vector<Mat> planes;
	//	cv::split(mat, planes);
	//	vvData.resize(channels, vector<double>(mat.rows * mat.cols, 0.0));
	//			
	//	for(int c = 0; c < channels; c++)
	//	{
	//		int ncount = 0;
	//		for(int n = 0; n < mat.rows; n++)
	//		{
	//			if(ndepth == CV_8U)
	//			{
	//				uchar* pfirst = planes[c].ptr<uchar>(n);
	//				for(int m = 0; m < mat.cols; m++)
	//				{
	//					vvData[c][ncount] = (double)pfirst[m];
	//					ncount++;
	//				}
	//			}
	//			else if(ndepth == CV_8S)
	//			{
	//				char* pfirst = planes[c].ptr<char>(n);
	//				for(int m = 0; m < mat.cols; m++)
	//				{
	//					vvData[c][ncount] = (double)pfirst[m];
	//					ncount++;
	//				}
	//			}
	//			else if(ndepth == CV_16U)
	//			{
	//				ushort* pfirst = planes[c].ptr<ushort>(n);
	//				for(int m = 0; m < mat.cols; m++)
	//				{
	//					vvData[c][ncount] = (double)pfirst[m];
	//					ncount++;
	//				}
	//			}
	//			else if(ndepth == CV_16S)
	//			{
	//				short* pfirst = planes[c].ptr<short>(n);
	//				for(int m = 0; m < mat.cols; m++)
	//				{
	//					vvData[c][ncount] = (double)pfirst[m];
	//					ncount++;
	//				}
	//			}
	//			else if(ndepth == CV_32S)
	//			{
	//				int* pfirst = planes[c].ptr<int>(n);
	//				for(int m = 0; m < mat.cols; m++)
	//				{
	//					vvData[c][ncount] = (double)pfirst[m];
	//					ncount++;
	//				}
	//			}
	//			else if(ndepth == CV_32F)
	//			{
	//				float* pfirst = planes[c].ptr<float>(n);
	//				for(int m = 0; m < mat.cols; m++)
	//				{
	//					vvData[c][ncount] = (double)pfirst[m];
	//					ncount++;
	//				}
	//			}
	//			else if(ndepth == CV_64F)
	//			{
	//				double* pfirst = planes[c].ptr<double>(n);
	//				for(int m = 0; m < mat.cols; m++)
	//				{
	//					vvData[c][ncount] = (double)pfirst[m];
	//					ncount++;
	//				}
	//			}
	//		}
	//	}
	//}
	vector< SingleFrameFeatures > m_vFeatures;
	string m_szSeqName;
	int m_nCurrIndex;
	double m_nan;
};

