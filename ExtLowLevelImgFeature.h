#pragma once
#ifndef EXTLOWLEVELIMGFEATURE_H
#define EXTLOWLEVELIMGFEATURE_H

#include "Defs.h"
#include "LocalFeatuers.h"
#include "My_SLIC_Interface.h"


class ExtLowLevelImgFeature
{
public:
	ExtLowLevelImgFeature(void)
	{
	}

	~ExtLowLevelImgFeature(void)
	{
	}

	static void ParseFileName(string& szFileName, string& szFolder, string& szName, string& szExt, string& szDelim)
	{
		int nPos0001 = szFileName.find_last_of("/");
		int nPob0001 = szFileName.find_last_of("\\");
		int nPos001 = (nPob0001>nPos0001)?nPob0001:nPos0001;
		szDelim = (nPob0001>nPos0001)?"\\":"/";
		char delim = (nPob0001>nPos0001)?'\\':'/';
		while(delim == szFileName[szFileName.length()-1])
		{
			szFileName = szFileName.substr(0, szFileName.length()-1);
		}
		nPos001 = szFileName.find_last_of(szDelim);
		szFolder = (nPos001 <= 0)?"":szFileName.substr(0, nPos001);
		int nPos002 = szFileName.find_last_of(".");
		nPos002 = (nPos002 > 0)?nPos002:(szFileName.length());
		szName = szFileName.substr(nPos001 + 1, nPos002 - nPos001 - 1);
		szExt = szFileName.substr(nPos002);
	}

	// Guess SLIC data path
	void GuessSLICPath(string& szPathName, string& szSLICPathName)
	{
		string szFoldertmp;
		string szNametmp;
		string szDelimtmp;
		string szExttmp;
		string szSLICdirtmp;
		string szSLICdattmp;
		ParseFileName(szPathName, szFoldertmp, szNametmp, szExttmp, szDelimtmp);
		int TrainPos = szFoldertmp.find("Train");
		int OldTrainPos = TrainPos;
		while(0 < (TrainPos = (int)szFoldertmp.find("Train", TrainPos+1)))
		{
			OldTrainPos = TrainPos;
		}
		TrainPos = OldTrainPos;
		int TestPos = szFoldertmp.find("Test");
		int OldTestPos = TestPos;
		while(0 < (TestPos = (int)szFoldertmp.find("Test", TrainPos+1)))
		{
			OldTestPos = TestPos;
		}
		TestPos = OldTestPos;
		if(TrainPos >= 0 && TrainPos > TestPos)
		{
			szSLICdirtmp = szFoldertmp.substr(0, TrainPos) + "SLIC" + szDelimtmp + "Train" + szDelimtmp;
			szSLICdattmp = szSLICdirtmp + szNametmp + ".dat";
		}
		else if(TestPos >= 0 && TrainPos < TestPos)
		{
			szSLICdirtmp = szFoldertmp.substr(0, TestPos) + "SLIC" + szDelimtmp + "Test" + szDelimtmp;
			szSLICdattmp = szSLICdirtmp + szNametmp + ".dat";
		}
		szSLICPathName = szSLICdattmp;
	}


	static int ExtColorFeatures(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat fimg;
		img.convertTo(fimg, CV_32FC3);
		//===========================================================================================//
		// 08/01/2011, prepare for color histogram, You-Chi
		//===========================================================================================//
		bool bColorHist = false;
		const int channels[] = {0, 1, 2};
		const cv::Mat defaultmat = cv::Mat();
		int ndims = 3;
		const int histSize[] = {4, 8, 8};
		float rrange[] = {0,256};
		float grange[] = {0,256};
		float brange[] = {0,256};
		//const int histSize[] = {16, 4, 4};
		const float* ranges[] = {brange, grange, rrange};
		//float hrange[] = {0,180};
		//float srange[] = {0,1};
		//float vrange[] = {0,256};
		//const float* ranges[] = {hrange, srange, vrange};
		//===========================================================================================//
		matrix<double> UT30, mu;
		if(bColorHist)
		{
#ifdef _MSC_VER
			LoadMyBinMat("C:\\CollectedData\\Sample004\\Sample004\\Individual_RGB\\FEList_Covar_EigVecs_30.dat", UT30);
			LoadMyBinMat("C:\\CollectedData\\Sample004\\Sample004\\Individual_RGB\\FEList_Mean_debug.dat", mu);
#else
        	        LoadMyBinMat("/nv/pec1/ycheng31/local/include/AdjustGMMNew/FEList_Covar_EigVecs_30.dat", UT30);
                	LoadMyBinMat("/nv/pec1/ycheng31/local/include/AdjustGMMNew/FEList_Mean_debug.dat", mu);		
#endif
		}
		double dbNormalizedFactor = 256.0 / (double)(BLKSIZE * BLKSIZE);
		for(int x = 0; x < img.cols / BLKSIZE; x++)
		{
			cv::Scalar fmeanxy, stddevxy;
			cv::Scalar fmeanxyhsv, stddevxyhsv;
			cv::Scalar fmeanxyx, stddevxyx;
			cv::Scalar fmeanxyy, stddevxyy;
			cv::Mat colx = fimg/*fres*/.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/); 
			for(int y = 0; y < img.rows / BLKSIZE; y++)
			{	
				cv::Mat sub = colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				cv::Mat subhsv;
				//sub.convertTo(subf, CV_32F);
				//cv::cvtColor(sub, subhsv, CV_BGR2HSV);
				//cv::Scalar ss = cv::mean(subhsv);
				//cv::Scalar ss0 = cv::mean(sub);
				//===========================================================================================//
				// 08/01/2011, compute color histogram instead, You-Chi
				//===========================================================================================//
				if(bColorHist)
				{
					// Let's calculate color histogram instead
					cv::MatND hist;
					cv::calcHist(&sub, 1, channels, defaultmat, hist, ndims, histSize, ranges, true, false);
					vector<double> vdbTemp(4*8*8, 0.0);
					size_t ii = 0;
					double dbCount = 0;
					for(size_t i1 = 0; i1 < 4; ++i1)
					{
						for(size_t i2 = 0; i2 < 8; ++i2)
						{
							for(size_t i3 = 0; i3 < 8; ++i3)
							{
								//SparseMat::
								double val = (double) hist.at<float>(i1, i2, i3);
								dbCount += val;
								vdbTemp[ii] = val;
								++ii;
							}
						}
					}
					matrix<double> mvec;
					mvec.SetSize(vdbTemp.size(), 1);
					for(size_t ii = 0; ii < vdbTemp.size(); ++ii)
					{
						mvec(ii, 0) = vdbTemp[ii] * dbNormalizedFactor - mu(ii, 0);
					}
					matrix<double> mvec30 = UT30 * mvec;
					vector<double> vdbTemp30(mvec30.RowNo(), 0.0);
					for(size_t ii = 0; ii < mvec30.RowNo(); ++ii)
					{
						vdbTemp30[ii] = mvec30(ii, 0);
					}
					vvFeatures.push_back(vdbTemp30);
				}
				else
				{
					// Old version
					cv::Mat sub8bits;
					sub.convertTo(sub8bits, CV_8UC3);
					cv::Mat subhsv;
					cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
					cv::meanStdDev(sub, fmeanxy, stddevxy);
					cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);

					BlockFeatures blk;
					blk.cvmat2vec(sub, blk.m_vvData);
					blk.extractBlockMean();
					blk.extractBlockVars();

					//double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
					//double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
					double dbBMean = fmeanxy.val[0];
					double dbGMean = fmeanxy.val[1];
					double dbRMean = fmeanxy.val[2];
					//double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
					//double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
					//double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
					//double dbVMean = fmeanxyhsv.val[2];
					//double dbSMean = fmeanxyhsv.val[1];
					//double dbHMean = fmeanxyhsv.val[0];
					vector<double> vdbTemp;
					//vdbTemp.push_back(dbAvgStdev);
					//vdbTemp.push_back(dbIMean);
					vdbTemp.push_back(dbBMean);
					vdbTemp.push_back(dbGMean);
					vdbTemp.push_back(dbRMean);
					//double dbFE1 = (-0.5668) * (dbBMean - 93.72874) + (-0.5715) * (dbGMean - 96.96279) + (-0.5934) * (dbRMean - 100.3968);
					//double dbFE2 = (0.5921) * (dbBMean - 93.72874) + (0.2182) * (dbGMean - 96.96279) + (-0.7758) * (dbRMean - 100.3968);
					//double dbFE3 = (0.5729) * (dbBMean - 93.72874) + (-0.7910) * (dbGMean - 96.96279) + (0.2147) * (dbRMean - 100.3968);
					
					//vdbTemp.push_back(dbFE1);
					//vdbTemp.push_back(dbFE2);
					//vdbTemp.push_back(dbFE3);					
					//vdbTemp.push_back(dbExBMean);
					//vdbTemp.push_back(dbExGMean);
					//vdbTemp.push_back(dbExRMean);
					//vdbTemp.push_back(dbVMean);
					//vdbTemp.push_back(dbSMean);
					//vdbTemp.push_back(dbHMean);
					vvFeatures.push_back(vdbTemp);
				}
			}
		}
		int yrange = img.rows / BLKSIZE;
		return yrange;
	}

	static int ExtColorFeaturesDiag(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat fimg;
		img.convertTo(fimg, CV_32FC3);
		//===========================================================================================//
		// 08/01/2011, prepare for color histogram, You-Chi
		//===========================================================================================//
		bool bColorHist = false;
		const int channels[] = {0, 1, 2};
		const cv::Mat defaultmat = cv::Mat();
		int ndims = 3;
		const int histSize[] = {8, 8, 8};
		float rrange[] = {0,256};
		float grange[] = {0,256};
		float brange[] = {0,256};
		const float* ranges[] = {brange, grange, rrange};

		//===========================================================================================//
		for(int x = 0; x < img.cols / BLKSIZE; x++)
		{
			cv::Scalar fmeanxy, stddevxy;
			cv::Scalar fmeanxyhsv, stddevxyhsv;
			cv::Scalar fmeanxyx, stddevxyx;
			cv::Scalar fmeanxyy, stddevxyy;
			cv::Mat colx = fimg/*fres*/.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/); 
			for(int y = 0; y < img.rows / BLKSIZE; y++)
			{	
				cv::Mat sub = colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				//===========================================================================================//
				// 08/01/2011, compute color histogram instead, You-Chi
				//===========================================================================================//
				if(bColorHist)
				{
					// Let's calculate color histogram instead
					cv::MatND hist;
					cv::calcHist(&sub, 1, channels, defaultmat, hist, ndims, histSize, ranges, true, false);
					vector<double> vdbTemp(8*8*8, 0.0);
					size_t ii = 0;
					double dbCount = 0;
					for(size_t i1 = 0; i1 < 8; ++i1)
					{
						for(size_t i2 = 0; i2 < 8; ++i2)
						{
							for(size_t i3 = 0; i3 < 8; ++i3)
							{
								//SparseMat::
								double val = (double) hist.at<float>(i1, i2, i3);
								dbCount += val;
								vdbTemp[ii] = val;
								++ii;
								
							}
						}
					}
					vvFeatures.push_back(vdbTemp);
				}
				else
				{
					// Old version
					cv::Mat sub8bits;
					sub.convertTo(sub8bits, CV_8UC3);
					cv::Mat subhsv;
					cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
					cv::meanStdDev(sub, fmeanxy, stddevxy);
					cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);

					BlockFeatures blk;
					blk.cvmat2vec(sub, blk.m_vvData);
					blk.extractBlockMean();
					blk.extractBlockVars();

					//double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
					//double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
					double dbBMean = fmeanxy.val[0];
					double dbGMean = fmeanxy.val[1];
					double dbRMean = fmeanxy.val[2];
					double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
					double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
					double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
					//double dbVMean = fmeanxyhsv.val[2];
					//double dbSMean = fmeanxyhsv.val[1];
					//double dbHMean = fmeanxyhsv.val[0];
					vector<double> vdbTemp;
					//vdbTemp.push_back(dbAvgStdev);
					//vdbTemp.push_back(dbIMean);
					vdbTemp.push_back(dbBMean);
					vdbTemp.push_back(dbGMean);
					vdbTemp.push_back(dbRMean);
					//vdbTemp.push_back(dbExBMean);
					//vdbTemp.push_back(dbExGMean);
					//vdbTemp.push_back(dbExRMean);
					//vdbTemp.push_back(dbVMean);
					//vdbTemp.push_back(dbSMean);
					//vdbTemp.push_back(dbHMean);
					vvFeatures.push_back(vdbTemp);
				}
			}
		}
		int yrange = img.rows / BLKSIZE;
		return yrange;
	}

	
	//static int ExtColorFeatures_BS0(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pRefImg, vector< vector<double> >& vvCorrs)
	static int ExtColorFeatures_BS0(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pRefImg, vector< vector<double> >& vvCorrs, vector<double>& vKSs)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat refimg(pRefImg);
		cv::Mat fimg;
		img.convertTo(fimg, CV_32FC3);
		//===========================================================================================//
		// 08/01/2011, prepare for color histogram, You-Chi
		//===========================================================================================//
		bool bColorHist = false;
		const int channels[] = {0, 1, 2};
		const cv::Mat defaultmat = cv::Mat();
		int ndims = 3;
		const int histSize[] = {4, 8, 8};
		float rrange[] = {0,256};
		float grange[] = {0,256};
		float brange[] = {0,256};
		const float* ranges[] = {brange, grange, rrange};
		//const int histSize[] = {16, 4, 4};
		//float hrange[] = {0,180};
		//float srange[] = {0,1};
		//float vrange[] = {0,256};
		//const float* ranges[] = {hrange, srange, vrange};
		matrix<double> UT30, mu;
		if(bColorHist)
		{
#ifdef _MSC_VER
			LoadMyBinMat("C:\\CollectedData\\Sample004\\Sample004\\Individual_RGB\\FEList_Covar_EigVecs_30.dat", UT30);
			LoadMyBinMat("C:\\CollectedData\\Sample004\\Sample004\\Individual_RGB\\FEList_Mean_debug.dat", mu);
#else
        	        LoadMyBinMat("/nv/pec1/ycheng31/local/include/AdjustGMMNew/FEList_Covar_EigVecs_30.dat", UT30);
                	LoadMyBinMat("/nv/pec1/ycheng31/local/include/AdjustGMMNew/FEList_Mean_debug.dat", mu);		
#endif
		}
		double dbNormalizedFactor = 256.0 / (double)(BLKSIZE * BLKSIZE);
		//===========================================================================================//
		for(int x = 0; x < img.cols / BLKSIZE; x++)
		{
			cv::Scalar fmeanxy, stddevxy;
			cv::Scalar fmeanxyhsv, stddevxyhsv;
			cv::Scalar fmeanxyx, stddevxyx;
			cv::Scalar fmeanxyy, stddevxyy;
			cv::Mat colx = fimg/*fres*/.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/); 
			cv::Mat refcolx = refimg/*fres*/.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/); 
			for(int y = 0; y < img.rows / BLKSIZE; y++)
			{	
				cv::Mat sub = colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				cv::Mat refsub = refcolx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				//===========================================================================================//
				// 08/01/2011, compute color histogram instead, You-Chi
				//===========================================================================================//
				// Old version
				cv::Mat sub8bits;
				sub.convertTo(sub8bits, CV_8UC3);
				cv::Mat subhsv;
				//cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
				cv::meanStdDev(sub, fmeanxy, stddevxy);
				//cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);

				BlockFeatures blk;
				blk.cvmat2vec(sub, blk.m_vvData);
				blk.extractBlockMean();
				blk.extractBlockVars();

				BlockFeatures refblk;
				refblk.cvmat2vec(refsub, refblk.m_vvData);
				refblk.extractBlockMean();
				refblk.extractBlockVars();
				blk.compBSCorrelation(refblk);
				
				vvCorrs.push_back(blk.m_vBSCorr);

				refblk.extractBlock_V_CDF();
				blk.extractBlock_V_CDF();					
				vKSs.push_back(CompKSStat(refblk.m_vVCDF, blk.m_vVCDF));

				if(bColorHist)
				{
					// Let's calculate color histogram instead
					cv::MatND hist;
					cv::calcHist(&sub, 1, channels, defaultmat, hist, ndims, histSize, ranges, true, false);
					vector<double> vdbTemp(4*8*8, 0.0);
					size_t ii = 0;
					double dbCount = 0;
					for(size_t i1 = 0; i1 < 4; ++i1)
					{
						for(size_t i2 = 0; i2 < 8; ++i2)
						{
							for(size_t i3 = 0; i3 < 8; ++i3)
							{
								//SparseMat::
								double val = (double) hist.at<float>(i1, i2, i3);
								dbCount += val;
								vdbTemp[ii] = val;
								++ii;
								
							}
						}
					}
					matrix<double> mvec;
					mvec.SetSize(vdbTemp.size(), 1);
					for(size_t ii = 0; ii < vdbTemp.size(); ++ii)
					{
						mvec(ii, 0) = vdbTemp[ii] * dbNormalizedFactor - mu(ii, 0);
					}
					matrix<double> mvec30 = UT30 * mvec;
					vector<double> vdbTemp30(mvec30.RowNo(), 0.0);
					for(size_t ii = 0; ii < mvec30.RowNo(); ++ii)
					{
						vdbTemp30[ii] = mvec30(ii, 0);
					}
					vvFeatures.push_back(vdbTemp30);
				}
				else
				{
					//double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
					//double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
					double dbBMean = fmeanxy.val[0];
					double dbGMean = fmeanxy.val[1];
					double dbRMean = fmeanxy.val[2];
					//double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
					//double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
					//double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
					//double dbVMean = fmeanxyhsv.val[2];
					//double dbSMean = fmeanxyhsv.val[1];
					//double dbHMean = fmeanxyhsv.val[0];
					vector<double> vdbTemp;
					//vdbTemp.push_back(dbAvgStdev);
					//vdbTemp.push_back(dbIMean);
					//double dbFE1 = (-0.5668) * (dbBMean - 93.72874) + (-0.5715) * (dbGMean - 96.96279) + (-0.5934) * (dbRMean - 100.3968);
					//double dbFE2 = (0.5921) * (dbBMean - 93.72874) + (0.2182) * (dbGMean - 96.96279) + (-0.7758) * (dbRMean - 100.3968);
					//double dbFE3 = (0.5729) * (dbBMean - 93.72874) + (-0.7910) * (dbGMean - 96.96279) + (0.2147) * (dbRMean - 100.3968);
					//vdbTemp.push_back(dbFE1);
					//vdbTemp.push_back(dbFE2);
					//vdbTemp.push_back(dbFE3);					
					vdbTemp.push_back(dbBMean);
					vdbTemp.push_back(dbGMean);
					vdbTemp.push_back(dbRMean);
					//vdbTemp.push_back(dbExBMean);
					//vdbTemp.push_back(dbExGMean);
					//vdbTemp.push_back(dbExRMean);
					//vdbTemp.push_back(dbVMean);
					//vdbTemp.push_back(dbSMean);
					//vdbTemp.push_back(dbHMean);
					vvFeatures.push_back(vdbTemp);
				}
			}
		}
		int yrange = img.rows / BLKSIZE;
		return yrange;
	}




	static int ExtColorFeatureRGBOnly(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat fimg;
		img.convertTo(fimg, CV_32FC3);

		for(int x = 0; x < img.cols / BLKSIZE; x++)
		{
			cv::Scalar fmeanxy, stddevxy;
			cv::Scalar fmeanxyhsv, stddevxyhsv;
			cv::Scalar fmeanxyx, stddevxyx;
			cv::Scalar fmeanxyy, stddevxyy;
			cv::Mat colx = fimg/*fres*/.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/); 
			for(int y = 0; y < img.rows / BLKSIZE; y++)
			{	
				cv::Mat sub = colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				cv::Mat sub8bits;
				sub.convertTo(sub8bits, CV_8UC3);
				cv::Mat subhsv;
				//cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
				cv::meanStdDev(sub, fmeanxy, stddevxy);
				//cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);
				//double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
				//double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
				double dbBMean = fmeanxy.val[0];
				double dbGMean = fmeanxy.val[1];
				double dbRMean = fmeanxy.val[2];
				//double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
				//double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
				//double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
				//double dbVMean = fmeanxyhsv.val[2];
				//double dbSMean = fmeanxyhsv.val[1];
				//double dbHMean = fmeanxyhsv.val[0];
				vector<double> vdbTemp;
				//vdbTemp.push_back(dbAvgStdev);
				//vdbTemp.push_back(dbIMean);
				vdbTemp.push_back(dbBMean);
				vdbTemp.push_back(dbGMean);
				vdbTemp.push_back(dbRMean);
				//vdbTemp.push_back(dbExBMean);
				//vdbTemp.push_back(dbExGMean);
				//vdbTemp.push_back(dbExRMean);
				//vdbTemp.push_back(dbVMean);
				//vdbTemp.push_back(dbSMean);
				//vdbTemp.push_back(dbHMean);
				vvFeatures.push_back(vdbTemp);
			}
		}
		int yrange = img.rows / BLKSIZE;
		fimg.release();
		return yrange;
	}

	static int ExtColorFeatureEXRGBHSV(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat fimg;
		img.convertTo(fimg, CV_32FC3);

		for(int x = 0; x < img.cols / BLKSIZE; x++)
		{
			cv::Scalar fmeanxy, stddevxy;
			cv::Scalar fmeanxyhsv, stddevxyhsv;
			cv::Scalar fmeanxyx, stddevxyx;
			cv::Scalar fmeanxyy, stddevxyy;
			cv::Mat colx = fimg/*fres*/.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/); 
			for(int y = 0; y < img.rows / BLKSIZE; y++)
			{	
				cv::Mat sub = colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				cv::Mat sub8bits;
				sub.convertTo(sub8bits, CV_8UC3);
				cv::Mat subhsv;
				//cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
				cv::meanStdDev(sub, fmeanxy, stddevxy);
				//cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);
				//double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
				//double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
				//double dbBMean = fmeanxy.val[0];
				//double dbGMean = fmeanxy.val[1];
				//double dbRMean = fmeanxy.val[2];
				double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
				double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
				double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
				double dbVMean = fmeanxyhsv.val[2];
				double dbSMean = fmeanxyhsv.val[1];
				double dbHMean = fmeanxyhsv.val[0];
				vector<double> vdbTemp;
				//vdbTemp.push_back(dbAvgStdev);
				//vdbTemp.push_back(dbIMean);
				//vdbTemp.push_back(dbBMean);
				//vdbTemp.push_back(dbGMean);
				//vdbTemp.push_back(dbRMean);
				vdbTemp.push_back(dbExBMean);
				vdbTemp.push_back(dbExGMean);
				vdbTemp.push_back(dbExRMean);
				vdbTemp.push_back(dbVMean);
				vdbTemp.push_back(dbSMean);
				vdbTemp.push_back(dbHMean);
				vvFeatures.push_back(vdbTemp);
			}
		}
		int yrange = img.rows / BLKSIZE;
		return yrange;
	}




	static int ExtColorFeatures(string szFileName, vector< vector<double> >& vvFeatures, int BLKSIZE)
	{
		IplImage* pimg;
	#ifdef LPF
		CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
		cvSetReal2D(hlpf, 0, 0, 0.1);
		cvSetReal2D(hlpf, 0, 1, 0.1);
		cvSetReal2D(hlpf, 0, 2, 0.1);
		cvSetReal2D(hlpf, 1, 0, 0.1);
		cvSetReal2D(hlpf, 1, 1, 0.8);
		cvSetReal2D(hlpf, 1, 2, 0.1);
		cvSetReal2D(hlpf, 2, 0, 0.1);
		cvSetReal2D(hlpf, 2, 1, 0.1);
		cvSetReal2D(hlpf, 2, 2, 0.1);
		IplImage* pimg0 = cvLoadImage(szFileName.c_str());
		IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
		cvConvert(pimg0, pimg1);
		cvFilter2D(pimg1, pimg2, hlpf);
		cvConvert(pimg2, pimg);
		cvReleaseImage(&pimg0);
		cvReleaseImage(&pimg1);
		cvReleaseImage(&pimg2);
		cvReleaseMat(&hlpf);
	#else
		pimg = cvLoadImage(szFileName.c_str());
	#endif

		int yrange = ExtColorFeatures(pimg, vvFeatures, BLKSIZE);
		cvReleaseImage(&pimg);
		return yrange;
	}

	//static int ExtColorFeatures_BS0(string szFileName,  vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pref, vector< vector<double> >& vvCorrs)
	static int ExtColorFeatures_BS0(string szFileName,  vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pref, vector< vector<double> >& vvCorrs, vector<double>& vKSs)
	{
		IplImage* pimg;
	#ifdef LPF
		CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
		cvSetReal2D(hlpf, 0, 0, 0.1);
		cvSetReal2D(hlpf, 0, 1, 0.1);
		cvSetReal2D(hlpf, 0, 2, 0.1);
		cvSetReal2D(hlpf, 1, 0, 0.1);
		cvSetReal2D(hlpf, 1, 1, 0.8);
		cvSetReal2D(hlpf, 1, 2, 0.1);
		cvSetReal2D(hlpf, 2, 0, 0.1);
		cvSetReal2D(hlpf, 2, 1, 0.1);
		cvSetReal2D(hlpf, 2, 2, 0.1);
		IplImage* pimg0 = cvLoadImage(szFileName.c_str());
		IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
		cvConvert(pimg0, pimg1);
		cvFilter2D(pimg1, pimg2, hlpf);
		cvConvert(pimg2, pimg);
		cvReleaseImage(&pimg0);
		cvReleaseImage(&pimg1);
		cvReleaseImage(&pimg2);
		cvReleaseMat(&hlpf);
	#else
		pimg = cvLoadImage(szFileName.c_str());
	#endif
		//ExtColorFeatures_BS0(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pRefImg, vector< vector<double> >& vvCorrs)
		//int yrange = ExtColorFeatures_BS0(pimg, vvFeatures, BLKSIZE, pref, vvCorrs);
		int yrange = ExtColorFeatures_BS0(pimg, vvFeatures, BLKSIZE, pref, vvCorrs, vKSs);
		cvReleaseImage(&pimg);
		return yrange;
	}

	static int ExtColorFeatureRGBOnly(string szFileName, vector< vector<double> >& vvFeatures, int BLKSIZE)
	{
		IplImage* pimg;
	#ifdef LPF
		CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
		cvSetReal2D(hlpf, 0, 0, 0.1);
		cvSetReal2D(hlpf, 0, 1, 0.1);
		cvSetReal2D(hlpf, 0, 2, 0.1);
		cvSetReal2D(hlpf, 1, 0, 0.1);
		cvSetReal2D(hlpf, 1, 1, 0.8);
		cvSetReal2D(hlpf, 1, 2, 0.1);
		cvSetReal2D(hlpf, 2, 0, 0.1);
		cvSetReal2D(hlpf, 2, 1, 0.1);
		cvSetReal2D(hlpf, 2, 2, 0.1);
		IplImage* pimg0 = cvLoadImage(szFileName.c_str());
		IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
		cvConvert(pimg0, pimg1);
		cvFilter2D(pimg1, pimg2, hlpf);
		cvConvert(pimg2, pimg);
		cvReleaseImage(&pimg0);
		cvReleaseImage(&pimg1);
		cvReleaseImage(&pimg2);
		cvReleaseMat(&hlpf);
	#else
		pimg = cvLoadImage(szFileName.c_str());
	#endif
		int yrange = ExtColorFeatureRGBOnly(pimg, vvFeatures, BLKSIZE);
		cvReleaseImage(&pimg);
		return yrange;
	}

	static int ExtColorFeatureEXRGBHSV(string szFileName, vector< vector<double> >& vvFeatures, int BLKSIZE)
	{
		IplImage* pimg;
	#ifdef LPF
		CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
		cvSetReal2D(hlpf, 0, 0, 0.1);
		cvSetReal2D(hlpf, 0, 1, 0.1);
		cvSetReal2D(hlpf, 0, 2, 0.1);
		cvSetReal2D(hlpf, 1, 0, 0.1);
		cvSetReal2D(hlpf, 1, 1, 0.8);
		cvSetReal2D(hlpf, 1, 2, 0.1);
		cvSetReal2D(hlpf, 2, 0, 0.1);
		cvSetReal2D(hlpf, 2, 1, 0.1);
		cvSetReal2D(hlpf, 2, 2, 0.1);
		IplImage* pimg0 = cvLoadImage(szFileName.c_str());
		IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
		cvConvert(pimg0, pimg1);
		cvFilter2D(pimg1, pimg2, hlpf);
		cvConvert(pimg2, pimg);
		cvReleaseImage(&pimg0);
		cvReleaseImage(&pimg1);
		cvReleaseImage(&pimg2);
		cvReleaseMat(&hlpf);
	#else
		pimg = cvLoadImage(szFileName.c_str());
	#endif
		int yrange = ExtColorFeatureEXRGBHSV(pimg, vvFeatures, BLKSIZE);
		cvReleaseImage(&pimg);
		return yrange;
	}
	static void LoadMyBinMat(string szPathName, matrix<double>& mat)
	{
		FILE* f = NULL;
		int n1 = 0;
		int n2 = 0;
		double dbBuf = 0.0;
		f = fopen(szPathName.c_str(), "rb");
		if(!f)
		{
			cerr<<"Cannot open "<<szPathName<<endl;
			return;
		}
		fread(&n1, sizeof(int), 1, f);
		fread(&n2, sizeof(int), 1, f);
		mat.SetSize(n1, n2);
		for(int i1 = 0; i1 < n1; ++i1)
		{
			for(int i2 = 0; i2 < n2; ++i2)
			{
				fread(&dbBuf, sizeof(double), 1, f);
				mat(i1, i2) = dbBuf;
			}
		}
		fclose(f);
	}
	
	static void WriteMyBinMat(string szPathName, matrix<double>& mat)
	{
		FILE* f = NULL;
		int n1 = (int)mat.RowNo();
		int n2 = (int)mat.ColNo();
		double dbBuf = 0.0;
		f = fopen(szPathName.c_str(), "wb");
		if(!f)
		{
			cerr<<"Cannot open "<<szPathName<<endl;
			return;
		}
		fwrite(&n1, sizeof(int), 1, f);
		fwrite(&n2, sizeof(int), 1, f);
		for(int i1 = 0; i1 < n1; ++i1)
		{
			for(int i2 = 0; i2 < n2; ++i2)
			{
				dbBuf = mat(i1, i2);
				fwrite(&dbBuf, sizeof(double), 1, f);
			}
		}
		fclose(f);
	}

	static int ExtColorFeatureSLICRGBOnly(IplImage*& pimg, vector< vector<double> >& vvFeatures, vector< vector<cv::Point> >& vvSLICLabels, int nsp, cv::Size& imsize, bool SLIC_is_given=false)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat fimg;
		img.convertTo(fimg, CV_32FC3);
		imsize = img.size();
		if(!SLIC_is_given) {
			My_SLIC_Interface slic_interface;
			slic_interface.Eval_SLICO_Single_Frame(img, nsp, vvSLICLabels);
		}
		int nsps = (int)vvSLICLabels.size();
		int n_my_sps = 0;
		for(int n = 0; n < nsps; ++n)
		{
			cv::Scalar fmeanxy;
			fmeanxy.val[0] = 0.0;
			fmeanxy.val[1] = 0.0;
			fmeanxy.val[2] = 0.0;
			int nnps = (int)vvSLICLabels[n].size();
			if(nnps > 0)
			{
				double dbnps = (double)nnps;
				for(int m = 0; m < nnps; ++m)
				{
					Vec3f fval = fimg.at<Vec3f>(vvSLICLabels[n][m].x, vvSLICLabels[n][m].y);
					fmeanxy.val[0] += fval.val[0];
					fmeanxy.val[1] += fval.val[1];
					fmeanxy.val[2] += fval.val[2];
				}
				fmeanxy.val[0] /= dbnps;
				fmeanxy.val[1] /= dbnps;
				fmeanxy.val[2] /= dbnps;
				double dbBMean = fmeanxy.val[0];
				double dbGMean = fmeanxy.val[1];
				double dbRMean = fmeanxy.val[2];
				vector<double> vdbTemp;
				vdbTemp.push_back(dbBMean);
				vdbTemp.push_back(dbGMean);
				vdbTemp.push_back(dbRMean);
				vvFeatures.push_back(vdbTemp);
				n_my_sps++;
			}
		}
		fimg.release();
		return n_my_sps;
	}


	static int ExtColorFeatureSLICRGBOnly(string szFileName, vector< vector<double> >& vvFeatures, vector< vector<cv::Point> >& vvSLICLabels, int nsp, cv::Size& imsize, bool SLIC_is_given=false)
	{
		IplImage* pimg;
	#ifdef LPF
		CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
		cvSetReal2D(hlpf, 0, 0, 0.1);
		cvSetReal2D(hlpf, 0, 1, 0.1);
		cvSetReal2D(hlpf, 0, 2, 0.1);
		cvSetReal2D(hlpf, 1, 0, 0.1);
		cvSetReal2D(hlpf, 1, 1, 0.8);
		cvSetReal2D(hlpf, 1, 2, 0.1);
		cvSetReal2D(hlpf, 2, 0, 0.1);
		cvSetReal2D(hlpf, 2, 1, 0.1);
		cvSetReal2D(hlpf, 2, 2, 0.1);
		IplImage* pimg0 = cvLoadImage(szFileName.c_str());
		IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
		cvConvert(pimg0, pimg1);
		cvFilter2D(pimg1, pimg2, hlpf);
		cvConvert(pimg2, pimg);
		cvReleaseImage(&pimg0);
		cvReleaseImage(&pimg1);
		cvReleaseImage(&pimg2);
		cvReleaseMat(&hlpf);
	#else
		pimg = cvLoadImage(szFileName.c_str());
	#endif
		int nsps = ExtColorFeatureSLICRGBOnly(pimg, vvFeatures, vvSLICLabels, nsp, imsize, SLIC_is_given);
		cvReleaseImage(&pimg);
		return nsps;
	}

	static int ExtColorFeatureSLICRGBOnly_LoadLabel(IplImage*& pimg, vector< vector<double> >& vvFeatures, string& szSLICLabelPathName, vector< vector<cv::Point> >& vvSLICLabels, cv::Size& imsize)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat fimg;
		img.convertTo(fimg, CV_32FC3);
		imsize = img.size();
		
		My_SLIC_Interface::LoadSLICDat(szSLICLabelPathName, imsize.width, imsize.height, vvSLICLabels);
		int nsps = (int)vvSLICLabels.size();
		int n_my_sps = 0;
		for(int n = 0; n < nsps; ++n)
		{
			cv::Scalar fmeanxy;
			fmeanxy.val[0] = 0.0;
			fmeanxy.val[1] = 0.0;
			fmeanxy.val[2] = 0.0;
			int nnps = (int)vvSLICLabels[n].size();
			if(nnps > 0)
			{
				double dbnps = (double)nnps;
				for(int m = 0; m < nnps; ++m)
				{
					Vec3f fval = fimg.at<Vec3f>(vvSLICLabels[n][m].x, vvSLICLabels[n][m].y);
					fmeanxy.val[0] += fval.val[0];
					fmeanxy.val[1] += fval.val[1];
					fmeanxy.val[2] += fval.val[2];
				}
				fmeanxy.val[0] /= dbnps;
				fmeanxy.val[1] /= dbnps;
				fmeanxy.val[2] /= dbnps;
				double dbBMean = fmeanxy.val[0];
				double dbGMean = fmeanxy.val[1];
				double dbRMean = fmeanxy.val[2];
				vector<double> vdbTemp;
				vdbTemp.push_back(dbBMean);
				vdbTemp.push_back(dbGMean);
				vdbTemp.push_back(dbRMean);
				vvFeatures.push_back(vdbTemp);
				n_my_sps++;
			}
		}
		fimg.release();
		return n_my_sps;
	}

	static int ExtColorFeatureSLICRGBOnly_LoadLabel(string szFileName, vector< vector<double> >& vvFeatures, string& szSLICLabelPathName)
	{
		IplImage* pimg;
	#ifdef LPF
		CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
		cvSetReal2D(hlpf, 0, 0, 0.1);
		cvSetReal2D(hlpf, 0, 1, 0.1);
		cvSetReal2D(hlpf, 0, 2, 0.1);
		cvSetReal2D(hlpf, 1, 0, 0.1);
		cvSetReal2D(hlpf, 1, 1, 0.8);
		cvSetReal2D(hlpf, 1, 2, 0.1);
		cvSetReal2D(hlpf, 2, 0, 0.1);
		cvSetReal2D(hlpf, 2, 1, 0.1);
		cvSetReal2D(hlpf, 2, 2, 0.1);
		IplImage* pimg0 = cvLoadImage(szFileName.c_str());
		IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
		cvConvert(pimg0, pimg1);
		cvFilter2D(pimg1, pimg2, hlpf);
		cvConvert(pimg2, pimg);
		cvReleaseImage(&pimg0);
		cvReleaseImage(&pimg1);
		cvReleaseImage(&pimg2);
		cvReleaseMat(&hlpf);
	#else
		pimg = cvLoadImage(szFileName.c_str());
	#endif
		vector< vector<cv::Point> > vvSLICLabels;
		cv::Size imsize;
		int nsps = ExtColorFeatureSLICRGBOnly_LoadLabel(pimg, vvFeatures, szSLICLabelPathName, vvSLICLabels, imsize);
		cvReleaseImage(&pimg);
		return nsps;
	}

	//static int ExtColorFeatures_SLIC_BS0(string szFileName,  vector< vector<double> >& vvFeatures, vector< vector<cv::Point> >& vvSLICLabels, int nsp, IplImage*& pref, vector< vector<double> >& vvCorrs, cv::Size& imsize)
	static int ExtColorFeatures_SLIC_BS0(string szFileName,  vector< vector<double> >& vvFeatures, vector< vector<cv::Point> >& vvSLICLabels, int nsp, IplImage*& pref, vector< vector<double> >& vvCorrs, vector<double>& vKSs, cv::Size& imsize, bool SLIC_is_given = false)
	{
		IplImage* pimg;
	#ifdef LPF
		CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
		cvSetReal2D(hlpf, 0, 0, 0.1);
		cvSetReal2D(hlpf, 0, 1, 0.1);
		cvSetReal2D(hlpf, 0, 2, 0.1);
		cvSetReal2D(hlpf, 1, 0, 0.1);
		cvSetReal2D(hlpf, 1, 1, 0.8);
		cvSetReal2D(hlpf, 1, 2, 0.1);
		cvSetReal2D(hlpf, 2, 0, 0.1);
		cvSetReal2D(hlpf, 2, 1, 0.1);
		cvSetReal2D(hlpf, 2, 2, 0.1);
		IplImage* pimg0 = cvLoadImage(szFileName.c_str());
		IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
		pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
		cvConvert(pimg0, pimg1);
		cvFilter2D(pimg1, pimg2, hlpf);
		cvConvert(pimg2, pimg);
		cvReleaseImage(&pimg0);
		cvReleaseImage(&pimg1);
		cvReleaseImage(&pimg2);
		cvReleaseMat(&hlpf);
	#else
		pimg = cvLoadImage(szFileName.c_str());
	#endif
		//ExtColorFeatures_BS0(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pRefImg, vector< vector<double> >& vvCorrs)
		//int nsps = ExtColorFeatures_SLIC_BS0(pimg, vvFeatures, vvSLICLabels, nsp, pref, vvCorrs, imsize);
		int nsps = ExtColorFeatures_SLIC_BS0(pimg, vvFeatures, vvSLICLabels, nsp, pref, vvCorrs, vKSs, imsize, SLIC_is_given);
		cvReleaseImage(&pimg);
		return nsps;
	}
	static void CompVCDF(cv::Mat& img, vector<double>& vVCDF)
	{
		cv::Mat imghsv, imgv;
		cv::cvtColor(img, imghsv, CV_BGR2HSV); 
		double dbImPix = (double)(imghsv.rows * imghsv.cols);
		const int channels[] = {2};
		const cv::Mat defaultmat = cv::Mat();
		int ndims = 1;
		const int histSize[] = {256};
		float vrange[] = {0,256};
		const float* ranges[] = {vrange};
		//vector<float> hist;
		cv::Mat mathist;
		cv::calcHist(&imghsv, 1, channels, defaultmat, mathist, ndims, histSize, ranges, true, false);
		//cv::normalize(mathist, mathist);
		vVCDF = vector<double>(256, 0);
		for(int n = 0; n < 256; ++n)
		{
			float* elem = mathist.ptr<float>(n, 0);
			vVCDF[n] = (((double)(*elem))/dbImPix);
			if(n >= 1)
				vVCDF[n] += vVCDF[n-1];
		}
		//vVCDF.resize(hist.size());
		//for(size_t n = 0; n < hist.size(); ++n)
		//	vVCDF[n] = (double)hist[n];
	}
	static void CompVCDF(vector<cv::Point>& vSLICLabels, cv::Mat& img, vector<double>& vVCDF)
	{
		vector<double> vhist(256, 0.0);
		vVCDF.reserve(256);
		vVCDF.resize(256);
		size_t nSize = vSLICLabels.size();
		double dbN = (double)nSize;
		for(size_t m = 0; m < nSize; ++m)
		{
			Vec3b bval = img.at<Vec3b>(vSLICLabels[m].x, vSLICLabels[m].y);
			uchar vval = bval.val[0];
			vval = (bval.val[1] > vval)?bval.val[1]:vval;
			vval = (bval.val[2] > vval)?bval.val[2]:vval;
			vhist[vval]++;
		}
		for(size_t i = 0; i < 256; ++i)
		{
			vhist[i] /= dbN;
			if(i > 1)
				vVCDF[i] += vVCDF[i-1];
			vVCDF[i] += vhist[i];
		}
	}
	static double CompKSStat(vector<double>& vVCDFRef, vector<double>& vVCDF)
	{
		size_t nSize = vVCDFRef.size();
		if(nSize != vVCDF.size())
		{
			cerr<<"Histogram size mismatched for computing KS Stat\n";
			return 0.0;
		}
		double dbMax = 0.0;
		for(size_t n = 0; n < nSize; ++n)
		{
			double diff = fabs(vVCDFRef[n] - vVCDF[n]);
			if(diff > dbMax)
				dbMax = diff;
		}
		return dbMax;
	}

	//static int ExtColorFeatures_SLIC_BS0(IplImage*& pimg, vector< vector<double> >& vvFeatures, vector< vector<cv::Point> >& vvSLICLabels, int nsp, IplImage*& pRefImg, vector< vector<double> >& vvCorrs, cv::Size& imsize)
	static int ExtColorFeatures_SLIC_BS0(IplImage*& pimg, vector< vector<double> >& vvFeatures, vector< vector<cv::Point> >& vvSLICLabels, int nsp, IplImage*& pRefImg, vector< vector<double> >& vvCorrs, vector<double>& vKSs, cv::Size& imsize, bool SLIC_is_given=false)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat refimg(pRefImg);
		cv::Mat fimg, frefimg;
		//vector<double> vVCDFAll;
		//CompVCDF(img, vVCDFAll);
		img.convertTo(fimg, CV_32FC3);
		imsize = img.size();
		refimg.convertTo(frefimg, CV_32FC3);
		if(!SLIC_is_given) {
			My_SLIC_Interface slic_interface;
			slic_interface.Eval_SLICO_Single_Frame(img, nsp, vvSLICLabels);
		}
		int nsps = (int)vvSLICLabels.size();
		for(int n = 0; n < nsps; ++n)
		{
			cv::Scalar fmeanxy, fvarxy, frefmeanxy, frefvarxy, fcovxy;
			vector<double> vCorr(3);
			for(int cc = 0; cc < 3; ++cc)
			{
				fmeanxy.val[cc] = 0.0;
				fvarxy.val[cc] = 0.0;
				frefmeanxy.val[cc] = 0.0;
				frefvarxy.val[cc] = 0.0;
				fcovxy.val[cc] = 0.0;			
			}
			int nnps = (int)vvSLICLabels[n].size();
			double dbnps = (double)nnps;
			if(nnps > 0)
			{
				vector<double> vVCDF;
				CompVCDF(vvSLICLabels[n], img, vVCDF);
				vector<double> vVCDFRef;
				CompVCDF(vvSLICLabels[n], refimg, vVCDFRef);
				vKSs.push_back(CompKSStat(vVCDFRef, vVCDF));
				for(int m = 0; m < nnps; ++m)
				{
					Vec3f fval = fimg.at<Vec3f>(vvSLICLabels[n][m].x, vvSLICLabels[n][m].y);
					for(int cc = 0; cc < 3; ++cc)
						fmeanxy.val[cc] += fval.val[cc];
					Vec3f frefval = frefimg.at<Vec3f>(vvSLICLabels[n][m].x, vvSLICLabels[n][m].y);
					for(int cc = 0; cc < 3; ++cc)
						frefmeanxy.val[cc] += frefval.val[cc];
				}
				for(int cc = 0; cc < 3; ++cc)
				{
					fmeanxy.val[cc] /= dbnps;
					frefmeanxy.val[cc] /= dbnps;
				}
				for(int m = 0; m < nnps; ++m)
				{
					Vec3f fval = fimg.at<Vec3f>(vvSLICLabels[n][m].x, vvSLICLabels[n][m].y);
					Vec3f frefval = frefimg.at<Vec3f>(vvSLICLabels[n][m].x, vvSLICLabels[n][m].y);
					for(int cc = 0; cc < 3; ++cc)
					{
						fvarxy.val[cc] += (fval.val[cc] - fmeanxy.val[cc]) * (fval.val[cc] - fmeanxy.val[cc]);
						frefvarxy.val[cc] += (frefval.val[cc] - frefmeanxy.val[cc]) * (frefval.val[cc] - frefmeanxy.val[cc]);
						fcovxy.val[cc] += (fval.val[cc] - fmeanxy.val[cc]) * (frefval.val[cc] - frefmeanxy.val[cc]);						
					}
				}
				double trnorm = 10.0;
				for(int cc = 0; cc < 3; ++cc)
				{
					double dbVarcc = sqrt(fvarxy.val[cc] * frefvarxy.val[cc]);
					if(fvarxy.val[cc] < trnorm && frefvarxy.val[cc] < trnorm) // constant to constant map
					{
						vCorr[cc] = 1.0;
					}
					else
					{
						vCorr[cc] = fabs((dbVarcc > 0)?(fcovxy.val[cc] / dbVarcc):0.0);
					}
				}
				vvCorrs.push_back(vCorr);
				if(dbnps > 1.0)
				{
					for(int cc = 0; cc < 3; ++cc)
					{
						fvarxy.val[cc] /= (dbnps - 1.0);
						frefvarxy.val[cc] /= (dbnps - 1.0);
						fcovxy.val[cc] /= (dbnps - 1.0);
					}
				}
				else
				{
					for(int cc = 0; cc < 3; ++cc)
					{
						fvarxy.val[cc] /= (dbnps);
						frefvarxy.val[cc] /= (dbnps);
						fcovxy.val[cc] /= (dbnps);
					}
				}
				vector<double> vdbTemp;
				for(int cc = 0; cc < 3; ++cc)
					vdbTemp.push_back(fmeanxy.val[cc]);
				vvFeatures.push_back(vdbTemp);
			}
		}
		return nsps;
	}

	//static int ExtColorFeatures_BS0_SLIC_MASK(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pRefImg, vector< vector<double> >& vvCorrs, unordered_map<int, int>& hValids, cv::Mat& SLIC_mask)
	static int ExtColorFeatures_BS0_SLIC_MASK(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pRefImg, vector< vector<double> >& vvCorrs, vector<double>& vKSs, unordered_map<int, int>& hValids, cv::Mat& SLIC_mask)
	{
		if(pimg == NULL)
		{
			cerr<<"No image is loaded for feature extraction\n";
			return 0;
		}
		cv::Mat img(pimg);
		cv::Mat refimg(pRefImg);
		cv::Mat fimg;
		
		img.convertTo(fimg, CV_32FC3);
		//===========================================================================================//
		int nIndex_Tmp = 0;
		for(int x = 0; x < img.cols / BLKSIZE; x++)
		{
			cv::Scalar fmeanxy, stddevxy;
			cv::Scalar fmeanxyhsv, stddevxyhsv;
			cv::Scalar fmeanxyx, stddevxyx;
			cv::Scalar fmeanxyy, stddevxyy;
			cv::Mat colx = fimg/*fres*/.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/); 
			cv::Mat refcolx = refimg/*fres*/.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/);
			cv::Mat slic_mask_colx = SLIC_mask.colRange(x*BLKSIZE, (x+1)*BLKSIZE/*min((x+1)*8, img.cols - 1)*/);
			for(int y = 0; y < img.rows / BLKSIZE; y++)
			{	
				cv::Mat sub = colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				cv::Mat refsub = refcolx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				cv::Mat slic_mask_sub = slic_mask_colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
				//===========================================================================================//
				// 08/01/2011, compute color histogram instead, You-Chi
				//===========================================================================================//
				// Old version
				cv::Mat sub8bits;
				sub.convertTo(sub8bits, CV_8UC3);
				cv::Mat subhsv;
				cv::meanStdDev(sub, fmeanxy, stddevxy, slic_mask_sub);

				BlockFeatures blk;
				blk.cvmat2vec(sub, blk.m_vvData);
				blk.setmask(slic_mask_sub);
				blk.extractBlockMean();
				blk.extractBlockVars();

				BlockFeatures refblk;
				refblk.cvmat2vec(refsub, refblk.m_vvData);
				refblk.setmask(slic_mask_sub);
				refblk.extractBlockMean();
				refblk.extractBlockVars();
				blk.compBSCorrelation(refblk);
				if(!(Myisnan(blk.m_vMeans[0]) || Myisnan(blk.m_vMeans[1]) || Myisnan(blk.m_vMeans[2])))
				{
					refblk.extractBlock_V_CDF();
					blk.extractBlock_V_CDF();					
					vKSs.push_back(CompKSStat(refblk.m_vVCDF, blk.m_vVCDF));
					vvCorrs.push_back(blk.m_vBSCorr);
					double dbBMean = fmeanxy.val[0];
					double dbGMean = fmeanxy.val[1];
					double dbRMean = fmeanxy.val[2];
					vector<double> vdbTemp;
					vdbTemp.push_back(dbBMean);
					vdbTemp.push_back(dbGMean);
					vdbTemp.push_back(dbRMean);
					hValids[nIndex_Tmp] = vvFeatures.size();
					vvFeatures.push_back(vdbTemp);
				}
				nIndex_Tmp++;
				
			}
		}
		int yrange = img.rows / BLKSIZE;
		return yrange;
	}

};
#endif


