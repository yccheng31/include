#include "My_SLIC_Interface.h"


My_SLIC_Interface::My_SLIC_Interface(void)
	: m_n_sp_count(0)
	, m_n_compactness(0)
{
}


My_SLIC_Interface::~My_SLIC_Interface(void)
{
}


// main interface for doing SLIC
int My_SLIC_Interface::main(const std::string& sz_image_list, const std::string& sz_out_dir, const int n_sp_count)
{
	std::ifstream fin;
	fin.open(sz_image_list.c_str());
	if(fin.is_open())
	{
		std::string szbuf;
		while(fin>>szbuf)
		{
			m_vsz_picvec.push_back(szbuf);
		}
	}
	else
	{
		std::cerr<<"Cannot open image list file "<<sz_image_list<<std::endl;
		return -1;
	}
	if(!file_exist(sz_out_dir))
	{
		std::string sz_mkdir_sz_out_dir = "mkdir " + sz_out_dir;
		system(sz_mkdir_sz_out_dir.c_str());
	}
	m_n_sp_count = n_sp_count;
	m_sz_out_dir = sz_out_dir;
	FixDirPath(m_sz_out_dir);
	if(n_sp_count < 0)
	{
		std::cerr<<"No negative super-pixel number is allowed, set to default value "<<DEFAULT_NUM_SP<<std::endl;
		m_n_sp_count = DEFAULT_NUM_SP;
	}
	m_n_compactness = 0;
	for(size_t k = 0; k < m_vsz_picvec.size(); ++k)
	{
		cv::Mat im = cv::imread(m_vsz_picvec[k]);
		cv::Size imsize = im.size();
		int width = (int)imsize.width;
		int height = (int)imsize.height;
		int sz = width * height;
		if(m_n_sp_count > sz) 
		{
			cerr<<"Number of superpixels ("<<m_n_sp_count
				<<") exceeds number of pixels ("<<sz<<") in the image"<<endl;
		}
		else
		{
			int* labels = new int[sz];
			int numlabels = 0;
			// Actually do SLIC processing
			SLIC slic;
			UINT* p_img_data = NULL;
			cvMat2UINTbuf(im, p_img_data);
			//ofstream f_tmp_out;
			//f_tmp_out.open("C:\\debug_hand_detection\\debug_SLIC_output\\debug\\p_my_img_data.log");
			//for(UINT nn = 0; nn < sz; ++nn)
			//{
			//	UINT nn_a = p_img_data[nn]/16777216;
			//	UINT nn_r = (p_img_data[nn] - nn_a * 16777216)/65536;
			//	UINT nn_g = (p_img_data[nn] - nn_a * 16777216 - 65536 * nn_r)/256;
			//	UINT nn_b = (p_img_data[nn] - nn_a * 16777216 - 65536 * nn_r - 256 * nn_g);
			//	f_tmp_out<<nn_a<<"\t"<<nn_r<<"\t"<<nn_g<<"\t"<<nn_b<<endl;
			//}
			//f_tmp_out.close();
			//UINT* p_img_data = *pp_img_data;

			slic.PerformSLICO_ForGivenK(p_img_data, width, height, labels, numlabels, m_n_sp_count, m_n_compactness);//for a given number K of superpixels
			//slic.PerformSLICO_ForGivenStepSize(img, width, height, labels, numlabels, m_stepsize, m_compactness);//for a given grid step size
			//slic.DrawContoursAroundSegments(img, labels, width, height, 0);//for black contours around superpixels
			slic.DrawContoursAroundSegmentsTwoColors(p_img_data, labels, width, height);//for black-and-white contours around superpixels
			slic.SaveSuperpixelLabels(labels,width,height,m_vsz_picvec[k],m_sz_out_dir);
			SavePicture(m_vsz_picvec[k], width, height, p_img_data,"","_SLICO");
			// free temporary memory allocations
			if(labels) delete [] labels;
			if(p_img_data) {
				delete [] p_img_data;
				//pp_img_data = NULL;
			}
		}

		im.release();
		std::cout<<"Processed "<<k<<"-th file: \t"<<m_vsz_picvec[k]<<std::endl;
	}
	return 0;
}


// convert cv::Mat data into a space storing image data (new memory)
void My_SLIC_Interface::cvMat2UINTbuf(const cv::Mat& im, UINT*& p_img_data)
{
	if(p_img_data)
	{
		//UINT* p_img_data = *pp_img_data;
		delete [] p_img_data;
		p_img_data = NULL;
	}
	int height = im.size().height;
	int width = im.size().width;
	int sz = width * height;
	if(sz > 0)
	{
		
		//UINT* 
		p_img_data = new UINT [sz];
		//pp_img_data = &p_img_data;
		for(size_t h = 0; h < height; ++h)
		{
			for(size_t w = 0; w < width; ++w)
			{
				
				UINT bval = im.at<cv::Vec3b>(h,w)[0];
				UINT gval = im.at<cv::Vec3b>(h,w)[1];
				UINT rval = im.at<cv::Vec3b>(h,w)[2];
				p_img_data[(h * width + w)] = 255 * 16777216 + rval * 65536 + gval * 256 + bval;
			}
		}
	}
}


// append UINT buffer data into a new cv::Mat object
void My_SLIC_Interface::UINTbuf2cvMat(UINT* p_img_data, cv::Mat& im, const int width, const int height)
{
	//im = cv::Mat(height, width, CV_8UC4);
	im = cv::Mat(height, width, CV_8UC3);
	int sz = width * height;
	for(int nn = 0; nn < sz; ++nn)
	{
		UINT nn_a = p_img_data[nn]/16777216;
		UINT nn_r = (p_img_data[nn] - nn_a * 16777216)/65536;
		UINT nn_g = (p_img_data[nn] - nn_a * 16777216 - 65536 * nn_r)/256;
		UINT nn_b = (p_img_data[nn] - nn_a * 16777216 - 65536 * nn_r - 256 * nn_g);
		int h = nn / width;
		int w = nn  - h * width;
		cv::Vec3b val3b;
		val3b.val[0] = (uchar)nn_b;
		val3b.val[1] = (uchar)nn_g;
		val3b.val[2] = (uchar)nn_r;
		im.at<Vec3b>(h, w) = val3b;
		//cv::Vec4b val4b;
		//val4b.val[0] = (uchar)nn_a;
		//val4b.val[1] = (uchar)nn_r;
		//val4b.val[2] = (uchar)nn_g;
		//val4b.val[3] = (uchar)nn_b;
		//im.at<Vec4b>(h, w) = val4b;
	}
}


// Save SLIC result
void My_SLIC_Interface::SavePicture(const std::string& szOldImgPath, const int width, const int height, UINT* p_img_data, const string& prefix, const string& suffix)
{
	cv::Mat im_out;
	char fname[_MAX_FNAME];
	char extn[_MAX_FNAME];
	_splitpath(szOldImgPath.c_str(), NULL, NULL, fname, extn);
	UINTbuf2cvMat(p_img_data, im_out, width, height);
	std::string szOutImgName = m_sz_out_dir + prefix + fname + suffix + ".jpg";
	cv::imwrite(szOutImgName, im_out);
}


// Evaluate SLICO label and store the result in a vector of vector, 1st layer: label index, 2nd layer, points
void My_SLIC_Interface::Eval_SLICO_Single_Frame(const std::string& szImagePathName, const int n_sp_count , std::vector< std::vector<cv::Point> >& vvLabels)
{
	if(n_sp_count < 0)
	{
		std::cerr<<"No negative super-pixel number is allowed, set to default value "<<DEFAULT_NUM_SP<<std::endl;
		m_n_sp_count = DEFAULT_NUM_SP;
	}
	m_n_compactness = 0;
	cv::Mat im = cv::imread(szImagePathName);
	cv::Size imsize = im.size();
	int width = (int)imsize.width;
	int height = (int)imsize.height;
	int sz = width * height;
	if(m_n_sp_count > sz) 
	{
		cerr<<"Number of superpixels ("<<m_n_sp_count
			<<") exceeds number of pixels ("<<sz<<") in the image"<<endl;
	}
	else
	{
		int* labels = new int[sz];
		int numlabels = 0;
		// Actually do SLIC processing
		SLIC slic;
		UINT* p_img_data = NULL;
		cvMat2UINTbuf(im, p_img_data);
		slic.PerformSLICO_ForGivenK(p_img_data, width, height, labels, numlabels, m_n_sp_count, m_n_compactness);//for a given number K of superpixels
		vvLabels.resize(m_n_sp_count);
		for(UINT nn = 0; nn < sz; ++nn)
		{
			int h = nn / width;
			int w = nn - h * width;
			vvLabels[labels[nn]].push_back(cv::Point(h,w));
		}
		if(labels) delete [] labels;
		if(p_img_data) {
			delete [] p_img_data;
		}
		im.release();
	}
}

// Evaluate SLICO label and store the result in a vector of vector, 1st layer: label index, 2nd layer, points
void My_SLIC_Interface::Eval_SLICO_Single_Frame(const cv::Mat& im, const int n_sp_count , std::vector< std::vector<cv::Point> >& vvLabels)
{
	m_n_sp_count = n_sp_count;
	if(n_sp_count < 0)
	{
		std::cerr<<"No negative super-pixel number is allowed, set to default value "<<DEFAULT_NUM_SP<<std::endl;
		m_n_sp_count = DEFAULT_NUM_SP;
	}
	m_n_compactness = 0;

	cv::Size imsize = im.size();
	int width = (int)imsize.width;
	int height = (int)imsize.height;
	int sz = width * height;
	if(m_n_sp_count > sz) 
	{
		cerr<<"Number of superpixels ("<<m_n_sp_count
			<<") exceeds number of pixels ("<<sz<<") in the image"<<endl;
	}
	else
	{
		int* labels = new int[sz];
		int numlabels = 0;
		// Actually do SLIC processing
		SLIC slic;
		UINT* p_img_data = NULL;
		cvMat2UINTbuf(im, p_img_data);
		slic.PerformSLICO_ForGivenK(p_img_data, width, height, labels, numlabels, m_n_sp_count, m_n_compactness);//for a given number K of superpixels
		vvLabels.resize(m_n_sp_count);
		for(UINT nn = 0; nn < sz; ++nn)
		{
			int h = nn / width;
			int w = nn - h * width;
			vvLabels[labels[nn]].push_back(cv::Point(h,w));
		}
		if(labels) delete [] labels;
		if(p_img_data) {
			delete [] p_img_data;
		}
	}
}

// Load SLIC label file
bool My_SLIC_Interface::LoadSLICDat(std::string& szPathName, int W, int H, std::vector< std::vector<cv::Point> >& vvLabels)
{
	FILE* pFile = NULL;
	pFile = fopen(szPathName.c_str(), "rb");
	if(pFile)
	{
		int nbuf = 0;
		int nmax = 0;
		int count = 0;
		while(fread(&nbuf, 4, 1, pFile))
		{
			count++;
			if(nbuf > nmax) {
				nmax = nbuf;
			}
		}
		if(count != W*H)
		{
			cerr<<"Number of elements is "<<count<<", not "<<(W*H)<<endl;
			fclose(pFile);
			return false;
		}
		vvLabels.resize(count+1);
		count = 0;
		rewind(pFile);
		for(int y = 0; y < H; ++y)
		{
			for(int x = 0; x < W; ++x)
			{
				fread(&nbuf, 4, 1, pFile);
				cv::Point pt = cv::Point(y, x);
				vvLabels[nbuf].push_back(pt);
			}
		}
		fclose(pFile);
		return true;
	}
	return false;
}