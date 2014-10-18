#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;
#ifdef _WIN32
#include <hash_map>
#define unordered_map hash_map 
using namespace stdext;
#else
#include <tr1/unordered_map>
#define hash_map unordered_map
using namespace tr1;
#endif


#define MYBLKSIZE 16
#define NUMFRAMETR 4
class JPGLabelInfo
{
public:
	int m_nWidth;
	int m_nHeight;
	int m_nBlockWidth;
	int m_nBlockHeight;
	int** m_ppnLabel;
	int m_nSkinCount;
	int BLOCKSIZE;
	vector< vector<int> > m_vvOnes;
	string m_szPathName;

	JPGLabelInfo()
	{
		m_nWidth=0;
		m_nHeight=0;
		m_ppnLabel=0;
		BLOCKSIZE = MYBLKSIZE;
		m_nSkinCount = 0;
	}
	~JPGLabelInfo()
	{
		if(m_ppnLabel != NULL)
		{
			for(int n = 0; n < m_nBlockHeight; n++)
			{
				delete [] m_ppnLabel[n];
			}
			delete [] m_ppnLabel;
		}
		m_nWidth = 0;
		m_nHeight = 0;
		m_nBlockWidth = 0;
		m_nBlockHeight = 0;
	}
	void SetWidthHeight(int nWidth, int nHeight)
	{
		m_nWidth = nWidth;
		m_nHeight = nHeight;
		m_nBlockWidth = (int)floor((double)nWidth/(double)BLOCKSIZE);
		m_nBlockHeight = (int)floor((double)nHeight/(double)BLOCKSIZE);
	
		m_ppnLabel = new int*[m_nBlockHeight];
		for(int n = 0; n < m_nBlockHeight; n++)
		{
			m_ppnLabel[n] = new int[m_nBlockWidth];
			for(int m = 0; m < m_nBlockWidth; m++)
			{
				m_ppnLabel[n][m] = 0;
			}
		}
	}
	void SetWidthHeight(int nWidth, int nHeight, int BLKSIZE)
	{
		BLOCKSIZE = BLKSIZE;
		m_nWidth = nWidth;
		m_nHeight = nHeight;
		m_nBlockWidth = (int)floor((double)nWidth/(double)BLOCKSIZE);
		m_nBlockHeight = (int)floor((double)nHeight/(double)BLOCKSIZE);
	
		m_ppnLabel = new int*[m_nBlockHeight];
		for(int n = 0; n < m_nBlockHeight; n++)
		{
			m_ppnLabel[n] = new int[m_nBlockWidth];
			for(int m = 0; m < m_nBlockWidth; m++)
			{
				m_ppnLabel[n][m] = 0;
			}
		}
	}
	void SetAsSkinPixel(int x, int y)
	{
		if(m_ppnLabel && x < m_nBlockWidth && y < m_nBlockHeight)
		{
			m_ppnLabel[y][x] = 1;
			m_nSkinCount++;
		}
	}
	void SetAsNotSurePixel(int x, int y)
	{
		if(m_ppnLabel && x < m_nBlockWidth && y < m_nBlockHeight)
		{
			m_ppnLabel[y][x] = -1;
		}
	}
	void SetAsNonSkinPixel(int x, int y)
	{
		if(m_ppnLabel && x < m_nBlockWidth && y < m_nBlockHeight)
		{
			m_ppnLabel[y][x] = 0;
		}
	}
};

class LoadJPGLabelInfo
{
public:
	LoadJPGLabelInfo()
	{
		BLOCKSIZE = MYBLKSIZE;
		m_nCountX = 0;
		m_nCountO = 0;

	}
	unordered_map<string, int> m_hszJPGFiles;
	unordered_map<string, int> m_hszCurrJPGFiles;
	unordered_map<string, JPGLabelInfo> m_hJPGLavelInfo;
	unordered_map< string, vector<int> > m_hJPGClasses;
	vector<string> m_vszKeyOrder;
	unordered_map<string, int> m_hHyperKeys;
	string m_szLabelPathName;
	unordered_map<string, unordered_map<string, int> > m_hhTrainFiles;
	unordered_map<string, vector<int> > m_hvSizes;

	int BLOCKSIZE;
	int m_nCountX;
	int m_nCountO;
	void LoadJPGInfo()
	{
		//ifstream fin;
		//fin.open(m_szLabelPathName.c_str());
		//string szBuf;
		//while(fin>>szBuf)
		//{
		//	size_t stdot1 = szBuf.find_first_of(',');
		//	string szFileName = szBuf.substr(0, stdot1);
		//	szBuf = szBuf.substr(stdot1+1, szBuf.length() - stdot1 - 1);
		//	size_t stdot2 = szBuf.find_first_of(',');
		//	string szWidth = szBuf.substr(0, stdot2);
		//	int nWidth = atoi(szWidth.c_str());
		//	szBuf = szBuf.substr(stdot2+1, szBuf.length() - stdot2 - 1);
		//	size_t stdot3 = szBuf.find_first_of(',');
		//	string szHeight = szBuf.substr(0, stdot3);
		//	int nHeight = atoi(szHeight.c_str());
		//	szBuf = szBuf.substr(stdot3+1, szBuf.length() - stdot3 - 1);
		//	string szHyperKeyTmp = szFileName.substr(0, szFileName.find_last_of("_") + 1);
		//	m_hHyperKeys[szHyperKeyTmp]++;
		//	if(m_hszJPGFiles.find(szFileName) == m_hszJPGFiles.end())
		//	{
		//		m_vszKeyOrder.push_back(szFileName);
		//	}
		//	m_hszJPGFiles[szFileName] = atoi(szBuf.c_str());
		//	m_hJPGLavelInfo[szFileName].SetWidthHeight(nWidth * BLOCKSIZE, nHeight * BLOCKSIZE);
		//	m_hJPGLavelInfo[szFileName].m_szPathName = szFileName;
		//	fin>>szBuf;
		//	string szFileNameKey = szFileName.substr(szFileName.find_last_of("\\") + 1, szFileName.find_last_of("_")  - szFileName.find_last_of("\\"));
		//	int nCountLabHand = 0;
		//	while(0 != szBuf.compare("."))
		//	{
		//		size_t stdot01 = szBuf.find_first_of(',');
		//		string szx = szBuf.substr(0, stdot01);
		//		szBuf = szBuf.substr(stdot01+1, szBuf.length() - stdot01 - 1);
		//		size_t stdot02 = szBuf.find_first_of(',');
		//		string szy = szBuf.substr(0, stdot02);
		//		szBuf = szBuf.substr(stdot02+1, szBuf.length() - stdot02 - 1);
		//		string szv = szBuf;
		//		if(-1 == atoi(szv.c_str()))
		//		{
		//			m_hJPGLavelInfo[szFileName].SetAsNotSurePixel(atoi(szx.c_str()), atoi(szy.c_str()));
		//			m_nCountX++;
		//		}
		//		else if(1 == atoi(szv.c_str()))
		//		{
		//			m_hJPGLavelInfo[szFileName].SetAsSkinPixel(atoi(szx.c_str()), atoi(szy.c_str()));
		//			m_nCountO++;
		//			nCountLabHand++;
		//		}
		//		fin>>szBuf;
		//	}
		//	if(nCountLabHand > 0)
		//	{
		//		
		//		if(m_hhTrainFiles.find(szFileNameKey) == m_hhTrainFiles.end())
		//		{
		//			m_hhTrainFiles[szFileNameKey][szFileName] = nCountLabHand;
		//		}
		//		else
		//		{
		//			if(NUMFRAMETR > m_hhTrainFiles[szFileNameKey].size())
		//			{
		//				m_hhTrainFiles[szFileNameKey][szFileName] = nCountLabHand;
		//			}
		//		}
		//	}
		//}
		//fin.close();
		ifstream fin;
		fin.open(m_szLabelPathName.c_str());
		string szBuf;
		while(fin>>szBuf)
		{
			if(szBuf[szBuf.length()-1] == '\r')
				szBuf = szBuf.substr(0, szBuf.length()-1);
			size_t stdot1 = szBuf.find_first_of(',');
			string szFileName = szBuf.substr(0, stdot1);
			szBuf = szBuf.substr(stdot1+1, szBuf.length() - stdot1 - 1);
			size_t stdot2 = szBuf.find_first_of(',');
			string szWidth = szBuf.substr(0, stdot2);
			int nWidth = atoi(szWidth.c_str());
			szBuf = szBuf.substr(stdot2+1, szBuf.length() - stdot2 - 1);
			size_t stdot3 = szBuf.find_first_of(',');
			string szHeight = szBuf.substr(0, stdot3);
			int nHeight = atoi(szHeight.c_str());
			szBuf = szBuf.substr(stdot3+1, szBuf.length() - stdot3 - 1);
			if(m_hszJPGFiles.find(szFileName) == m_hszJPGFiles.end())
			{
				m_vszKeyOrder.push_back(szFileName);
			}
			int nnpostmp1 =(int)szFileName.find_last_of("\\");
			int nnpostmp11 =(int)szFileName.find_last_of("/");
			int nnpostmp01 = (nnpostmp1 > nnpostmp11)?nnpostmp1:nnpostmp11;
			int nnpostmp02 = (int)szFileName.find_last_of("_");
			string szTmpKey = szFileName.substr(nnpostmp01 + 1, nnpostmp02 - nnpostmp01);
			m_hvSizes[szTmpKey].push_back(nWidth);
			m_hvSizes[szTmpKey].push_back(nHeight);
			m_hszJPGFiles[szFileName] = atoi(szBuf.c_str());
			m_hJPGLavelInfo[szFileName].SetWidthHeight(nWidth * BLOCKSIZE, nHeight * BLOCKSIZE, BLOCKSIZE);
			m_hJPGLavelInfo[szFileName].m_szPathName = szFileName;
			m_hJPGClasses[szFileName] = vector<int>(nWidth*nHeight, 0);
			fin>>szBuf;
			if(szBuf[szBuf.length()-1] == '\r')
				szBuf = szBuf.substr(0, szBuf.length()-1);
			//string szFileNameKey = szFileName.substr(szFileName.find_last_of("\\") + 1, szFileName.find_last_of("_")  - szFileName.find_last_of("\\"));
			string szFileNameKey = szFileName.substr(szFileName.find_last_of("/") + 1, szFileName.find_last_of("_")  - szFileName.find_last_of("/"));
			int nCountLabHand = 0;
			while(0 != szBuf.compare("."))
			{
				size_t stdot01 = szBuf.find_first_of(',');
				string szx = szBuf.substr(0, stdot01);
				szBuf = szBuf.substr(stdot01+1, szBuf.length() - stdot01 - 1);
				size_t stdot02 = szBuf.find_first_of(',');
				string szy = szBuf.substr(0, stdot02);
				szBuf = szBuf.substr(stdot02+1, szBuf.length() - stdot02 - 1);
				string szv = szBuf;
				int ncurrx = atoi(szx.c_str());
				int ncurry = atoi(szy.c_str());
				if(-1 == atoi(szv.c_str()))
				{
					m_hJPGLavelInfo[szFileName].SetAsNotSurePixel(ncurrx, ncurry);
					m_hJPGClasses[szFileName][ncurrx * nHeight + ncurry] = -1;
					m_nCountX++;
				}
				else if(1 == atoi(szv.c_str()))
				{
					m_hJPGLavelInfo[szFileName].SetAsSkinPixel(ncurrx, ncurry);
					m_hJPGClasses[szFileName][ncurrx * nHeight + ncurry] = 1;
					vector<int> v2(2,0);
					v2[0] = ncurrx;
					v2[1] = ncurry;
					m_hJPGLavelInfo[szFileName].m_vvOnes.push_back(v2);
					m_nCountO++;
					nCountLabHand++;
				}
				fin>>szBuf;
				if(szBuf[szBuf.length()-1] == '\r')
					szBuf = szBuf.substr(0, szBuf.length()-1);
			}
			if(nCountLabHand > 0)
			{
				
				if(m_hhTrainFiles.find(szFileNameKey) == m_hhTrainFiles.end())
				{
					m_hhTrainFiles[szFileNameKey][szFileName] = nCountLabHand;
				}
				else
				{
					if(NUMFRAMETR > m_hhTrainFiles[szFileNameKey].size())
					{
						m_hhTrainFiles[szFileNameKey][szFileName] = nCountLabHand;
					}
				}
			}
		}
		fin.close();

	}
	void BruteForcelyFixKeys(vector<string>& vszKeys)
	{
		if(0 == vszKeys.size())
		{
			cerr<<"Empty vector, aborted!"<<endl;
			return;
		}
		for(size_t n = 0; n < vszKeys.size(); n++)
		{
			if(m_hJPGLavelInfo.find(vszKeys[n]) == m_hJPGLavelInfo.end())
			{
				int nnpostmp1 =(int)vszKeys[n].find_last_of("\\");
				int nnpostmp11 =(int)vszKeys[n].find_last_of("/");
				int nnpostmp01 = (nnpostmp1 > nnpostmp11)?nnpostmp1:nnpostmp11;
				int nnpostmp02 = (int)vszKeys[n].find_last_of("_");
				string szTmpKey = vszKeys[n].substr(nnpostmp01 + 1, nnpostmp02 - nnpostmp01);
				int nSize1 = 0;
				int nSize2 = 0;
				if(m_hvSizes.find(szTmpKey) == m_hvSizes.end())
				{
					cerr<<"Warning, the sequence was not found in original stream, size set to the same as the fist one"<<endl;
					if(2 == m_hvSizes.begin()->second.size())
					{
						nSize1 = m_hvSizes.begin()->second[0];
						nSize2 = m_hvSizes.begin()->second[1];
					}
					else
					{
						cerr<<"OMG, you don't even have the first one, lets give you a fake size"<<endl;
						nSize1 = 16;
						nSize2 = 12;
					}
				}
				else
				{
					nSize1 = m_hvSizes[szTmpKey][0];
					nSize2 = m_hvSizes[szTmpKey][1];
				}
				m_hJPGLavelInfo[vszKeys[n]].SetWidthHeight(nSize1 * BLOCKSIZE, nSize2 * BLOCKSIZE, BLOCKSIZE);
				m_hJPGLavelInfo[vszKeys[n]].m_szPathName = vszKeys[n];
				m_hJPGClasses[vszKeys[n]] = vector<int>(nSize1*nSize2, 0);
			}
			else
			{

			}
			
		}
		m_vszKeyOrder = vszKeys;
	}
	void SaveJPGInfo(string szLabelPathName)
	{
		//string& szLabelPathName = m_szLabelPathName;
		unordered_map<string, JPGLabelInfo>& hJPGLavelInfo = m_hJPGLavelInfo;
		unordered_map<string, int>& hszJPGFiles = m_hszJPGFiles;
		vector<string>& vszJPGFiles = m_vszKeyOrder;
		ofstream fout;
		fout.open(szLabelPathName.c_str());
		size_t stFiles = hJPGLavelInfo.size();
		size_t stCurrFileLists = vszJPGFiles.size();
		vector<string> vszOthers;
		int nIndex = 0;
		for(int nnn = 0; nnn < (int)vszJPGFiles.size(); nnn++)
		{
			if(hJPGLavelInfo.find(vszJPGFiles[nnn]) != hJPGLavelInfo.end())
			{
				nIndex = nnn;
			}
		}
		//for(unordered_map<string, JPGLabelInfo>::iterator itr = m_hJPGLavelInfo.begin();itr!=m_hJPGLavelInfo.end(); itr++)
		for(size_t st = 0; st < (size_t)nIndex; st++)
		{
			//string szKey = itr->first;
			string szKey = vszJPGFiles[st];
			fout<<hJPGLavelInfo[szKey].m_szPathName<<","<<hJPGLavelInfo[szKey].m_nBlockWidth
				<<","<<hJPGLavelInfo[szKey].m_nBlockHeight<<","<<hszJPGFiles[szKey]<<endl;
			for(int mm = 0; mm < hJPGLavelInfo[szKey].m_nBlockWidth; mm++)
			{
				for(int nn = 0; nn < hJPGLavelInfo[szKey].m_nBlockHeight; nn++)
				{
					if(hJPGLavelInfo[szKey].m_ppnLabel[nn][mm] == 1)
					{
						fout<<mm<<","<<nn<<","<<1<<endl;
					}
				}
			}
			for(int mm = 0; mm < hJPGLavelInfo[szKey].m_nBlockWidth; mm++)
			{
				for(int nn = 0; nn < hJPGLavelInfo[szKey].m_nBlockHeight; nn++)
				{
					if(hJPGLavelInfo[szKey].m_ppnLabel[nn][mm] == -1)
					{
						fout<<mm<<","<<nn<<","<<-1<<endl;
					}
				}
			}
			fout<<"."<<endl;
		}
		fout.close();
	}

	void MarkSelectedPointWithoutDrawing(string szCurrFile, int x, int y, int val, 
		int nWidth,
		int nHeight)
	{
		unordered_map<string, JPGLabelInfo>& hJPGLavelInfo = m_hJPGLavelInfo;
		int& nCountO = m_nCountO;
		int& nCountX = m_nCountX;

		int BLKSIZE = BLOCKSIZE;

		if(hJPGLavelInfo.find(szCurrFile) == hJPGLavelInfo.end())
		{
			hJPGLavelInfo[szCurrFile].SetWidthHeight(nWidth*BLKSIZE, nHeight*BLKSIZE, BLKSIZE);
			hJPGLavelInfo[szCurrFile].m_szPathName = szCurrFile;
		}
		if(val == 0)
		{
			hJPGLavelInfo[szCurrFile].SetAsNonSkinPixel(x, y);
		}
		else if(val == 1)
		{
			if(hJPGLavelInfo[szCurrFile].m_ppnLabel[y][x] == 1)
			{
				hJPGLavelInfo[szCurrFile].SetAsNonSkinPixel(x, y);
				nCountO--;
			}
			else
			{
				if(hJPGLavelInfo[szCurrFile].m_ppnLabel[y][x] == -1)
				{
					nCountX--;
				}
				hJPGLavelInfo[szCurrFile].SetAsSkinPixel(x, y);
				nCountO++;
			}
		}
		else if(val == -1)
		{
			if(hJPGLavelInfo[szCurrFile].m_ppnLabel[y][x] == -1)
			{
				hJPGLavelInfo[szCurrFile].SetAsNonSkinPixel(x, y);
				nCountX--;
			}
			else
			{
				if(hJPGLavelInfo[szCurrFile].m_ppnLabel[y][x] == 1)
				{
					nCountO--;
				}
				hJPGLavelInfo[szCurrFile].SetAsNotSurePixel(x, y);
				nCountX++;
			}
		}
	}


};


class LoadJPGSLICLabelInfo
{
public:
	LoadJPGSLICLabelInfo()
	{
	}
	unordered_map<string, string> m_hszJPGFiles;
	unordered_map< string, vector<int> > m_hJPGClasses;
	unordered_map<string, unordered_map<int, int> > m_hJPGSLICLabelInfo;
	unordered_map<string, vector<int> > m_hJPGSLICLabelID;
	vector<string> m_vszKeyOrder;
	string m_szLabelPathName;

	void LoadJPGSLICLabel()
	{
		ifstream fin;
		fin.open(m_szLabelPathName.c_str());
		string szBuf;
		while(fin>>szBuf)
		{
			if(szBuf[szBuf.length()-1] == '\r')
				szBuf = szBuf.substr(0, szBuf.length()-1);
			string szFileName = szBuf;
			if(m_hszJPGFiles.find(szFileName) == m_hszJPGFiles.end())
			{
				m_vszKeyOrder.push_back(szFileName);
			}
			int nnpostmp1 =(int)szFileName.find_last_of("\\");
			int nnpostmp11 =(int)szFileName.find_last_of("/");
			int nnpostmp01 = (nnpostmp1 > nnpostmp11)?nnpostmp1:nnpostmp11;
			int nnpostmp02 = (int)szFileName.find_last_of("_");
			string szTmpKey = szFileName.substr(nnpostmp01 + 1, nnpostmp02 - nnpostmp01);
			m_hszJPGFiles[szFileName] = szTmpKey;
			fin>>szBuf;
			if(szBuf[szBuf.length()-1] == '\r')
				szBuf = szBuf.substr(0, szBuf.length()-1);
			string szFileNameKey = szFileName.substr(szFileName.find_last_of("/") + 1, szFileName.find_last_of("_")  - szFileName.find_last_of("/"));
			unordered_map<int, int> htmp;
			while(0 != szBuf.compare("."))
			{
				int SLIC_ID = atoi(szBuf.c_str());
				htmp[SLIC_ID] = 1;
				m_hJPGSLICLabelID[szFileName].push_back(SLIC_ID);
				fin>>szBuf;
				if(szBuf[szBuf.length()-1] == '\r')
					szBuf = szBuf.substr(0, szBuf.length()-1);
			}
			m_hJPGSLICLabelInfo[szFileName] = htmp;
		}
		fin.close();
	}

	void SaveJPGSLICLabel(string szLabelPathName)
	{
		unordered_map<string, unordered_map<int, int> >& hJPGLabelInfo = m_hJPGSLICLabelInfo;
		unordered_map<string, string>& hszJPGFiles = m_hszJPGFiles;
		vector<string>& vszJPGFiles = m_vszKeyOrder;

		ofstream fout;
		fout.open(szLabelPathName.c_str());
		size_t stFiles = hJPGLabelInfo.size();
		size_t stCurrFileLists = vszJPGFiles.size();
		vector<string> vszOthers;
		int nIndex = 0;
		for(int nnn = 0; nnn < (int)vszJPGFiles.size(); nnn++)
		{
			if(hJPGLabelInfo.find(vszJPGFiles[nnn]) != hJPGLabelInfo.end())
			{
				string szKey = vszJPGFiles[nnn];
				fout<<szKey<<endl;
				for(size_t st = 0; st < m_hJPGSLICLabelID[szKey].size(); ++st)
				{
					fout<<m_hJPGSLICLabelID[szKey][st]<<endl;
				}
				fout<<"."<<endl;
				
			}
		}
		fout.close();
	}
};
