#pragma once

#include "Defs.h"
#include "SLIC.h"
#define DEFAULT_NUM_SP 200

#ifndef UINT
typedef unsigned int        UINT;
#endif

class My_SLIC_Interface
{
public:
	My_SLIC_Interface(void);
	~My_SLIC_Interface(void);
	// main interface for doing SLIC
	int main(const std::string& sz_image_list, const std::string& sz_out_dir, const int n_sp_count);
	// store the list of image paths
	std::vector<std::string> m_vsz_picvec;
	// string for output directory
	std::string m_sz_out_dir;
	// number of super pixels
	int m_n_sp_count;
	// compactness factor
	int m_n_compactness;
	// convert cv::Mat data into a space storing image data (new memory)
	void cvMat2UINTbuf(const cv::Mat& im, UINT*& p_img_data);

	// Fix path by adding / or \ at the end if the path string does not end with / or \

	bool FixDirPath(std::string& szDirPath)
	{
		char lastChar = szDirPath.at( szDirPath.length() - 1 );
#ifdef WIN32
		if(lastChar != '\\')
		{
			szDirPath = szDirPath + "\\";
		}
#else
		if(lastChar != '/')
		{
			szDirPath = szDirPath + "/";
		}
#endif
		return false;
	}
	// append UINT buffer data into a new cv::Mat object
	void UINTbuf2cvMat(UINT* p_img_data, cv::Mat& im, const int width, const int height);
	// Save SLIC result
	void SavePicture(const std::string& szOldImgPath, const int width, const int height, UINT* p_img_data, const string& prefix, const string& suffix);
	// Evaluate SLICO label and store the result in a vector of vector, 1st layer: label index, 2nd layer, points
	void Eval_SLICO_Single_Frame(const std::string& szImagePathName, const int n_sp_count , std::vector< std::vector<cv::Point> >& vvLabels);
	// Evaluate SLICO label and store the result in a vector of vector, 1st layer: label index, 2nd layer, points
	void Eval_SLICO_Single_Frame(const cv::Mat& im, const int n_sp_count , std::vector< std::vector<cv::Point> >& vvLabels);
	// Load SLIC label file
	static bool LoadSLICDat(std::string& szPathName, int W, int H, std::vector< std::vector<cv::Point> >& vvLabels);
};

