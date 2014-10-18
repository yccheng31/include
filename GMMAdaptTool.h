#pragma once
#include "GMM.h"
#include "connected.h"
#include "LoadJPGLabelInfo.h"
#include "ExtLowLevelImgFeature.h"
#include "VideoFE.h"


//class CVideoFE;
//#include <tr1/unordered_map>
//using namespace std::tr1;
#define NUMITERATION 5
#define WIDTH 16
#define HEIGHT 12
//==============================================================================//
// inserted, 6/6/2011 You-Chi Cheng, and all 4096 are changed to MAXSTRLENGTH
//===============================================================================//
#define MAXSTRLENGTH 4096
#define CCThreshold 0.5
#define GMMITER 5
//===============================================================================//
extern vector<int> vDimConsidered;

//class JPGLabelInfo
//{
//public:
//	int m_nWidth;
//	int m_nHeight;
//	int m_nBlockWidth;
//	int m_nBlockHeight;
//	int** m_ppnLabel;
//	vector< vector<int> > m_vvOnes;
//	string m_szPathName;
//
//	JPGLabelInfo()
//	{
//		m_nWidth=0;
//		m_nHeight=0;
//		m_nBlockWidth = 0;
//		m_nBlockHeight = 0;
//		m_ppnLabel=0;
//	}
//	~JPGLabelInfo()
//	{
//		if(m_ppnLabel)
//		{
//			for(int n = 0; n < m_nBlockHeight; n++)
//			{
//				delete [] m_ppnLabel[n];
//			}
//			delete [] m_ppnLabel;
//		}
//		m_nWidth = 0;
//		m_nHeight = 0;
//		m_nBlockWidth = 0;
//		m_nBlockHeight = 0;
//	}
//	void SetWidthHeight(int nWidth, int nHeight, int BLOCKSIZE)
//	{
//		m_nWidth = nWidth;
//		m_nHeight = nHeight;
//		m_nBlockWidth = (int)floor((double)nWidth/(double)BLOCKSIZE);
//		m_nBlockHeight = (int)floor((double)nHeight/(double)BLOCKSIZE);
//	
//		m_ppnLabel = new int*[m_nBlockHeight];
//		for(int n = 0; n < m_nBlockHeight; n++)
//		{
//			m_ppnLabel[n] = new int[m_nBlockWidth];
//			for(int m = 0; m < m_nBlockWidth; m++)
//			{
//				m_ppnLabel[n][m] = 0;
//			}
//		}
//	}
//	void SetAsSkinPixel(int x, int y)
//	{
//		if(m_ppnLabel && x < m_nBlockWidth && y < m_nBlockHeight)
//		{
//			m_ppnLabel[y][x] = 1;
//			vector<int> vPt(2,0);
//			vPt[0] = x;
//			vPt[1] = y;
//			m_vvOnes.push_back(vPt);
//		}
//	}
//	void SetAsNotSurePixel(int x, int y)
//	{
//		if(m_ppnLabel && x < m_nBlockWidth && y < m_nBlockHeight)
//		{
//			m_ppnLabel[y][x] = -1;
//		}
//	}
//	void SetAsNonSkinPixel(int x, int y)
//	{
//		if(m_ppnLabel && x < m_nBlockWidth && y < m_nBlockHeight)
//		{
//			m_ppnLabel[y][x] = 0;
//		}
//	}
//};
//

class GMMAdaptTool
{
public:
	GMMAdaptTool(void)
	{	
	}

	~GMMAdaptTool(void)
	{
	}

	static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
		std::stringstream ss(s);
		std::string item;
		while(std::getline(ss, item, delim)) {
			elems.push_back(item);
		}
		return elems;
	}
	static std::vector<std::string> split(const std::string &s, char delim) {
		std::vector<std::string> elems;
		return split(s, delim, elems);
	}
	static std::vector<double> vecstr2vecdb(vector<string>& vdata)
	{
		vector<double> results = vector<double>(vdata.size(), 0.0);
		for(size_t n = 0; n < vdata.size(); n++)
		{
			results[n] = atof(vdata[n].c_str());
		}
		return results;
	}
	static void LoadLogData(string szPathName, vector< vector<double> >& vvdata, vector<int>& vclasses)
	{
		ifstream fin;
		fin.open(szPathName.c_str());
		char szbuf[65535];
		memset(szbuf, 0, 65535);
		if(fin.is_open())
		{
			while(fin.getline(szbuf, 65535))
			{
				string strbuf = string(szbuf);
				vector<string> data = split(strbuf, ',');
				vector<double> vdata = vecstr2vecdb(data);
				vector<double> vdatan = vector<double>(vdata.size() - 1, 0.0);
				for(size_t n = 1; n < vdata.size(); n++)
				{
					vdatan[n-1] = vdata[n];
				}
				vvdata.push_back(vdatan);
				vclasses.push_back(atoi(data[0].c_str()));
			}
			fin.close();
		}
		else
		{
			cerr<<"Cannot open "<<szPathName<<endl;
		}
	}

	static int LoadGMM(string szGMM, unordered_map<string, GMM>& hGMMs)
	{
		ifstream fin;
		fin.open(szGMM.c_str());
		string szBuf;
		unordered_map<string, GMM> hMixtureMacro;
		if(fin.is_open())
		{
			int nVecSize = 0;
			while(getline(fin, szBuf))
			{
				string strName = "";
				string szFeatureType = "";
				int nfindtilt = -1;
				int nMixes = 0;

				if(szBuf.npos != (nfindtilt = szBuf.find("~o")))
				{	
					while(1)
					{
						getline(fin, szBuf);
						int nfoundpos = -1;
						if(szBuf.npos != (nfoundpos = szBuf.find("<VECSIZE>")))
						{
							int npos01 = szBuf.find_last_of("<");
							int npos02 = szBuf.find_last_of(">");
							string szVecSize = szBuf.substr(nfoundpos + 9, npos01 - nfoundpos - 9);
							nVecSize = atoi(szVecSize.c_str());
							szFeatureType = szBuf.substr(npos01 + 1, npos02 - npos01 - 1);
							break;
						}
					}
				}
				// mixture macro
				else if(szBuf.npos != (nfindtilt = szBuf.find("~m")))
				{
					string szMacroName = szBuf.substr(nfindtilt+2, szBuf.length() - nfindtilt - 2);
					int npos1 = szMacroName.find_first_of('\"');
					int npos2 = szMacroName.find_last_of('\"');
					szMacroName = szMacroName.substr(npos1, npos2-npos1 + 1);
					int nVecSizeMacro = 0;
					fin>>szBuf>>nVecSizeMacro;
					std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
					if(szBuf.npos != szBuf.find("<MEAN>"))
					{
						hMixtureMacro[szMacroName].InitGMMMeans(nVecSizeMacro, 1);
						hMixtureMacro[szMacroName].m_vmdCovar.resize(1);
						hMixtureMacro[szMacroName].m_vmdInvCovar.resize(1);
						hMixtureMacro[szMacroName].m_vMixWeights[0] = 0;
						hMixtureMacro[szMacroName].m_vlogMixWeights[0] = LZERO;
						hMixtureMacro[szMacroName].m_vmdMean[0].Null();										
						for(int d = 0; d < nVecSizeMacro; d++)
						{
							double db1 = 0.0;
							fin>>db1;
							hMixtureMacro[szMacroName].m_vmdMean[0](d, 0) = db1;
						}
						fin>>szBuf>>nVecSizeMacro;
						std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
						if(szBuf.npos != szBuf.find("<VARIANCE>"))
						{
							hMixtureMacro[szMacroName].InitGMMVars(nVecSizeMacro, 0);
							//cout<<"in hGMMs["<<strName<<"], the size of the var is "<<hGMMs[strName].m_vmdCovar[nMixIndex].RowNo()<<"x"<<hGMMs[strName].m_vmdCovar[nMixIndex].ColNo()<<endl;
							for(int d = 0; d < nVecSizeMacro; d++)
							{
								fin>>hMixtureMacro[szMacroName].m_vmdCovar[0](d,d);
								hMixtureMacro[szMacroName].m_vmdInvCovar[0](d,d) = 1.0 /hMixtureMacro[szMacroName].m_vmdCovar[0](d,d);
							}
						}
						else if(szBuf.npos != szBuf.find("<INVCOVAR>"))
						{
							hMixtureMacro[szMacroName].InitGMMVars(nVecSizeMacro, 0, "full");
							for(int d = 0; d < nVecSizeMacro; d++)
							{
								for(int dd = d; dd < nVecSizeMacro; dd++)
								{
									fin>>hMixtureMacro[szMacroName].m_vmdInvCovar[0](d,dd);
									hMixtureMacro[szMacroName].m_vmdInvCovar[0](dd, d) = hMixtureMacro[szMacroName].m_vmdInvCovar[0](d,dd);
								}
							}
							hMixtureMacro[szMacroName].m_vmdCovar[0] = !hMixtureMacro[szMacroName].m_vmdInvCovar[0];
						}
						double dbgconst = 0.0;
						fin>>szBuf>>dbgconst;
						hMixtureMacro[szMacroName].CompGconsts();
					}
				}
				else if(szBuf.npos != (nfindtilt = szBuf.find("~h")))
				{
					strName = szBuf.substr(nfindtilt+2, szBuf.length() - nfindtilt - 2);
					int npos1 = strName.find_first_of('\"');
					int npos2 = strName.find_last_of('\"');
					strName = strName.substr(npos1+1, npos2-npos1-1);
					getline(fin, szBuf);
					while(1)
					{
						int nfoundpos = -1;
						bool bforcestop = false;
						if(szBuf.npos != szBuf.find('<') && szBuf.npos != szBuf.find('>'))
						{
							std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
							if(szBuf.npos != szBuf.find("<BEGINHMM>"))
							{
							}
							else if(szBuf.npos != (nfoundpos = szBuf.find("<VECSIZE>")))
							{
								int npos01 = szBuf.find_last_of("<");
								int npos02 = szBuf.find_last_of(">");
								string szVecSize = szBuf.substr(nfoundpos + 9, npos01 - nfoundpos - 9);
								nVecSize = atoi(szVecSize.c_str());
								szFeatureType = szBuf.substr(npos01 + 1, npos02 - npos01 - 1);
							}
							else if(szBuf.npos != (nfoundpos = szBuf.find("<NUMMIXES>")))
							{
								int npos01 = szBuf.length() - 1;
								string szMixes = szBuf.substr(nfoundpos + 10, npos01 - nfoundpos - 9);
								nMixes = atoi(szMixes.c_str());
								hGMMs[strName].InitGMMMeans(nVecSize, nMixes);
								int nMixIndex = 0;
								int nDims = 0;
								double dMixWeight;
								int nlastload = -1;
								hGMMs[strName].m_vmdCovar.resize(nMixes);
								hGMMs[strName].m_vmdInvCovar.resize(nMixes);
								
								for(int kk = 0; kk < nMixes; kk++)
								{
									//fin>>szBuf>>nMixIndex>>dMixWeight;
									fin>>szBuf;
									if(szBuf.npos != szBuf.find("<TRANSP>"))
									{
										bforcestop = true;
										for(int kkk = nlastload + 1; kkk < nMixes; kkk++)
										{
											hGMMs[strName].m_vMixWeights[kkk] = 0;
											hGMMs[strName].m_vlogMixWeights[kkk] = LZERO;
											hGMMs[strName].m_vmdMean[kkk].Null();										
											hGMMs[strName].m_vmdCovar[kkk].SetSize(nVecSize, nVecSize);
											hGMMs[strName].m_vmdCovar[kkk].Null();
											hGMMs[strName].m_vmdInvCovar[kkk].SetSize(nVecSize, nVecSize);
											hGMMs[strName].m_vmdInvCovar[kkk].Null();
										}
										break;
									}
									else if(szBuf.npos != szBuf.find("<MIXTURE>"))
									{
										fin>>nMixIndex>>dMixWeight;
										nMixIndex--;
										if(nMixIndex - nlastload > 1)
										{
											int nstartkkk = max(0, nlastload + 1);
											for(int kkk = nstartkkk; kkk < nMixIndex; kkk++)
											{
												hGMMs[strName].m_vMixWeights[kkk] = 0;
												hGMMs[strName].m_vlogMixWeights[kkk] = LZERO;
												hGMMs[strName].m_vmdMean[kkk].Null();
												hGMMs[strName].m_vmdCovar[kkk].SetSize(nVecSize, nVecSize);
												hGMMs[strName].m_vmdCovar[kkk].Null();
												hGMMs[strName].m_vmdInvCovar[kkk].SetSize(nVecSize, nVecSize);
												hGMMs[strName].m_vmdInvCovar[kkk].Null();
											}
											kk = nMixIndex;
										}
										nlastload = kk;										
										hGMMs[strName].m_vMixWeights[nMixIndex] = dMixWeight;
										hGMMs[strName].m_vlogMixWeights[nMixIndex] = log(dMixWeight);
										//fin>>szBuf>>nDims;
										getline(fin, szBuf);
										while(0 == szBuf.length() || 0 == szBuf.compare("\r"))
											getline(fin, szBuf);
										string szBuf02 = szBuf;
										std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
										if(szBuf.npos != szBuf.find("<MEAN>"))
										{
											stringstream ss;
											ss.str(szBuf);
											string sz1;
											ss>>sz1>>nDims;											
											//ss.open()
											//hGMMs[strName].InitGMMMeans(nVecSize, nMixes);
											for(int d = 0; d < nDims; d++)
											{
												double db1 = 0.0;
												fin>>db1;
												hGMMs[strName].m_vmdMean[nMixIndex](d, 0) = db1;
											}
											fin>>szBuf>>nDims;
											std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
											if(szBuf.npos != szBuf.find("<VARIANCE>"))
											{
												hGMMs[strName].InitGMMVars(nVecSize, nMixIndex);
												//cout<<"in hGMMs["<<strName<<"], the size of the var is "<<hGMMs[strName].m_vmdCovar[nMixIndex].RowNo()<<"x"<<hGMMs[strName].m_vmdCovar[nMixIndex].ColNo()<<endl;
												for(int d = 0; d < nDims; d++)
												{
													fin>>hGMMs[strName].m_vmdCovar[nMixIndex](d,d);
													hGMMs[strName].m_vmdInvCovar[nMixIndex](d,d) = 1.0 / hGMMs[strName].m_vmdCovar[nMixIndex](d,d);
												}
											}
											else if(szBuf.npos != szBuf.find("<INVCOVAR>"))
											{
												hGMMs[strName].InitGMMVars(nVecSize, nMixIndex, "full");
												for(int d = 0; d < nDims; d++)
												{
													for(int dd = d; dd < nDims; dd++)
													{
														fin>>hGMMs[strName].m_vmdInvCovar[nMixIndex](d,dd);
														hGMMs[strName].m_vmdInvCovar[nMixIndex](dd, d) = hGMMs[strName].m_vmdInvCovar[nMixIndex](d,dd);
													}
												}
												hGMMs[strName].m_vmdCovar[nMixIndex] = !hGMMs[strName].m_vmdInvCovar[nMixIndex];
											}
											double dbgconst = 0.0;
											fin>>szBuf>>dbgconst;
										}
										else if(szBuf.npos != szBuf.find("<VARIANCE>"))
										{
											stringstream ss;
											ss.str(szBuf);
											string sz1;
											ss>>sz1>>nDims;
											hGMMs[strName].InitGMMVars(nVecSize, nMixIndex);
											for(int d = 0; d < nDims; d++)
											{
												fin>>hGMMs[strName].m_vmdCovar[nMixIndex](d,d);
												hGMMs[strName].m_vmdInvCovar[nMixIndex](d,d) = 1.0 / hGMMs[strName].m_vmdCovar[nMixIndex](d,d);
											}
											fin>>szBuf>>nDims;
											std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
											if(szBuf.npos != szBuf.find("<MEAN>"))
											{
												for(int d = 0; d < nDims; d++)
												{
													fin>>hGMMs[strName].m_vmdMean[nMixIndex](d, 0);
												}
											}
											double dbgconst = 0.0;
											fin>>szBuf>>dbgconst;
										}
										else if(szBuf.npos != szBuf.find("<INVCOVAR>"))
										{
											stringstream ss;
											ss.str(szBuf);
											string sz1;
											ss>>sz1>>nDims;
											hGMMs[strName].InitGMMVars(nVecSize, nMixIndex, "full");
											for(int d = 0; d < nDims; d++)
											{
												for(int dd = d; dd < nDims; dd++)
												{
													fin>>hGMMs[strName].m_vmdInvCovar[nMixIndex](d,dd);
													hGMMs[strName].m_vmdInvCovar[nMixIndex](dd, d) = hGMMs[strName].m_vmdInvCovar[nMixIndex](d,dd);
												}
											}									
											hGMMs[strName].m_vmdCovar[nMixIndex] = !hGMMs[strName].m_vmdInvCovar[nMixIndex];
											fin>>szBuf>>nDims;
											std::transform(szBuf.begin(), szBuf.end(), szBuf.begin(), (int(*)(int)) std::toupper);
											if(szBuf.npos != szBuf.find("<MEAN>"))
											{
												for(int d = 0; d < nDims; d++)
												{
													fin>>hGMMs[strName].m_vmdMean[nMixIndex](d, 0);
												}
											}
											double dbgconst = 0.0;
											fin>>szBuf>>dbgconst;
										
										}
										else if(szBuf.npos != szBuf.find("~M"))
										{
											stringstream ss;
											ss.str(szBuf02);
											string sz1, sz2;
											ss>>sz1>>sz2;
											if(hMixtureMacro.find(sz2) == hMixtureMacro.end())
											{
												cerr<<"The specified macro was not defined!!!"<<endl;
												return -1;
											}
											else
											{
												
												hGMMs[strName].m_vmdMean[nMixIndex] = hMixtureMacro[sz2].m_vmdMean[0];
												hGMMs[strName].m_vmdCovar[nMixIndex] = hMixtureMacro[sz2].m_vmdCovar[0];
												hGMMs[strName].m_vmdInvCovar[nMixIndex] = hMixtureMacro[sz2].m_vmdInvCovar[0];
												if(hGMMs[strName].m_vdgconst.size() != nMixes)
													hGMMs[strName].m_vdgconst.resize(nMixes);
												if(hGMMs[strName].m_vvddiaggconst.size() != nMixes)
													hGMMs[strName].m_vvddiaggconst.resize(nMixes);
												hGMMs[strName].m_vdgconst[nMixIndex] =  hMixtureMacro[sz2].m_vdgconst[0];
												hGMMs[strName].m_vvddiaggconst[nMixIndex] =  hMixtureMacro[sz2].m_vvddiaggconst[0];
											}
										}
									}
								}
							}
						}
						int nNumStates = 0;
						if(bforcestop)
						{
							fin>>nNumStates;
						}
						else
						{
							getline(fin, szBuf);
						}
						if(szBuf.npos != szBuf.find("<TRANSP>"))
						{
							if(!bforcestop)
							{
								int npos00 = szBuf.find_last_of(">");
								string szNumStates = szBuf.substr(npos00 + 1, szBuf.length() - npos00);
								nNumStates = atoi(szNumStates.c_str());
							}
							hGMMs[strName].m_vvTransp.resize(nNumStates);
							for(int nt = 0; nt < nNumStates; nt++)
							{
								hGMMs[strName].m_vvTransp[nt].resize(nNumStates);
								double db1 = 0.0;
								for(int nt1 = 0; nt1 < nNumStates; nt1++)
								{
									fin>>db1;
									hGMMs[strName].m_vvTransp[nt][nt1] = db1;
								}
							}

						}
						if(szBuf.npos != szBuf.find("<ENDHMM>"))
						{
							break;
						}
					}
				}
			}
		}
		fin.close();
		unordered_map<string, GMM>::iterator itr;
		for(itr = hGMMs.begin(); itr != hGMMs.end(); itr++)
		{
			string szKey = itr->first;
			hGMMs[szKey].CompGconsts();
		}
		return 1;
	}

	static int WriteGMM(string szOutGMM, unordered_map<string, GMM>& hGMMs)
	{
			// Current version for Diag covariance only, with dormat tag "USER"
			ofstream fout;
		if(0 == hGMMs.size())
		{
			cerr<<"Cannot write empty model!"<<endl;
			return -1;
		}
			fout.open(szOutGMM.c_str());
		unordered_map<string, GMM>::iterator itr = hGMMs.begin();
			string firstgmm = itr->first;
			int nVecSize = hGMMs[firstgmm].m_nVecSize;
			if(fout.is_open())
			{
					fout<<"~o"<<endl;
					fout<<"<STREAMINFO> 1 "<<nVecSize<<endl;
					fout<<"<VECSIZE> "<<nVecSize<<"<NULLD><USER><DIAGC>"<<endl;
					for(;itr != hGMMs.end();itr++)
					{
							string szKey = itr->first;
							fout<<"~h \""<<szKey<<"\""<<endl;
							fout<<"<BEGINHMM>"<<endl;
							// We are gmms, we only have one state
							fout<<"<NUMSTATES> 3"<<endl;
							fout<<"<STATE> 2"<<endl;
				fout<<"<NUMMIXES> "<<hGMMs[szKey].m_nMixes<<endl;
				for(int n = 0; n < hGMMs[szKey].m_nMixes; n++)
				{
					if(hGMMs[szKey].m_vMixWeights[n] > 0)
					{
						fout<<"<MIXTURE> "<<(n+1)<<" ";
						fout.setf(ios::scientific, ios::floatfield);
						fout<<hGMMs[szKey].m_vMixWeights[n]<<endl;
						fout.unsetf(ios::floatfield);
						fout<<"<MEAN> "<<hGMMs[szKey].m_nVecSize<<endl;
						fout.setf(ios::scientific, ios::floatfield);
						for(int m = 0; m < hGMMs[szKey].m_nVecSize; m++)
						{
							fout<<" "<<hGMMs[szKey].m_vmdMean[n](m, 0);
						}
						fout<<endl;
						fout.unsetf(ios::floatfield);
						fout<<"<VARIANCE> "<<hGMMs[szKey].m_nVecSize<<endl;
						fout.setf(ios::scientific, ios::floatfield);
						for(int m = 0; m < hGMMs[szKey].m_nVecSize; m++)
										{
												fout<<" "<<hGMMs[szKey].m_vmdCovar[n](m, m);
										}
						fout<<endl;
						fout<<"<GCONST> "<<hGMMs[szKey].m_vdgconst[n]<<endl;
						fout.unsetf(ios::floatfield);
					}
				}
				fout<<"<TRANSP> "<<hGMMs[szKey].m_vvTransp.size()<<endl;
				fout.setf(ios::scientific, ios::floatfield);
				for(int n = 0; n < hGMMs[szKey].m_vvTransp.size(); n++)
				{
					for(int n1 = 0; n1 < hGMMs[szKey].m_vvTransp[n].size(); n1++)
					{
						fout<<" "<<hGMMs[szKey].m_vvTransp[n][n1];
					}
					fout<<endl;
				}
				fout<<"<ENDHMM>"<<endl;
					}
			}
			else
			{
					cerr<<"Cannot open "<<szOutGMM<<" for writing GMM"<<endl;
					return -1;
			}

			return 0;
	}

	static int GMMHash2VectBH(unordered_map<string, GMM>& hGMMs, vector<GMM>& vGMMs)
	{
		vGMMs.resize(2);
		vGMMs[0] = GMM(hGMMs["BkGnd"]);
		vGMMs[1] = GMM(hGMMs["Hand"]);
		return 1;
	}
	
	static int RecognizeByGMMs(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses)
	{
		int N = (int)vvData.size();
		if(N == 0)
		{
			cerr<<"Empty data for GMM recognition!"<<endl;
			return 0;
		}
		recogClasses.resize(N, 0);
		int I = (int)GMMs.size();
		for(int n = 0; n < N; n++)
		{
			double maxni = GMMs[0].CompLikelihood(vvData[n]);
			int argmaxni = 0;
			for(int i = 1; i < I ;i++)
			{
				double dbi = GMMs[i].CompLikelihood(vvData[n]);
				if(dbi > maxni)
				{
					maxni = dbi;
					argmaxni = i;
				}
			}
			recogClasses[n] = argmaxni;
		}
		return 1;
	}

	static int GMMsScores(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector< vector<double> >& vvrecogScores)
	{
		int N = (int)vvData.size();
		if(N == 0)
		{
			cerr<<"Empty data for GMM recognition!"<<endl;
			return 0;
		}
		int I = (int)GMMs.size();
		vvrecogScores.resize(N, vector<double>(I, 0.0));
		for(int n = 0; n < N; n++)
		{
			for(int i = 0; i < I ;i++)
			{
				vvrecogScores[n][i] = GMMs[i].CompLikelihood(vvData[n]);				
			}
		}
		return 1;
	}


	static int GMMsScoresRecog(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses, vector< vector<double> >& vvrecogScores)
	{
		int N = (int)vvData.size();
		if(N == 0)
		{
			cerr<<"Empty data for GMM recognition!"<<endl;
			return 0;
		}
		int I = (int)GMMs.size();
		vvrecogScores.resize(N, vector<double>(I, 0.0));
		recogClasses.resize(N, 0);
		for(int n = 0; n < N; n++)
		{
			double dbMaxLL = LZERO;
			for(int i = 0; i < I ;i++)
			{
				vvrecogScores[n][i] = GMMs[i].CompLikelihood(vvData[n]);
				if(vvrecogScores[n][i] > dbMaxLL)
				{
					dbMaxLL = vvrecogScores[n][i];
					recogClasses[n] = i;
				}
			}
		}
		return 1;
	}
	// Must make sure vvvDecisionTable is 256*256*256
	static int ApproxRecognizeByTable(vector< vector<double> >& vvData, bool*** vvvDecisionTable, vector<int>& recogClasses)
	{
		int N = (int)vvData.size();
		if(N == 0)
		{
			cerr<<"Empty data for recognition!"<<endl;
			return 0;
		}
		recogClasses.resize(N, 0);
		for(int n = 0; n < N; n++)
		{
			int nIndex1 = (int)floor(vvData[n][0] + 0.5);
			int nIndex2 = (int)floor(vvData[n][1] + 0.5);
			int nIndex3 = (int)floor(vvData[n][2] + 0.5);
			if(nIndex1 > 255)
				nIndex1 = 255;
			if(nIndex2 > 255)
				nIndex2 = 255;
			if(nIndex3 > 255)
				nIndex3 = 255;
				recogClasses[n] = (vvvDecisionTable[nIndex1][nIndex2][nIndex3])?1:0;
		}
		return 1;
	}	

	static int RecognizeByGMMsPart(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses)
	{
		int N = (int)vvData.size();
		if(N == 0)
		{
			cerr<<"Empty data for GMM recognition!"<<endl;
			return 0;
		}
		recogClasses.resize(N, 0);
		int I = (int)GMMs.size();
		for(int n = 0; n < N; n++)
		{
			double maxni = GMMs[0].CompLikelihoodPart(vvData[n], vDimConsidered);
			int argmaxni = 0;
			for(int i = 1; i < I ;i++)
			{
				double dbi = GMMs[i].CompLikelihoodPart(vvData[n], vDimConsidered);
				if(dbi > maxni)
				{
					maxni = dbi;
					argmaxni = i;
				}
			}
			recogClasses[n] = argmaxni;
		}
		return 1;
	}

	static int Test2ClassesGMMs(vector< vector<double> >& vvData, vector<int>& classes, vector<GMM>& GMMs, vector<int>& recogClasses, double& Acc, vector<int>& confres)
	{
		if(!RecognizeByGMMs(vvData, GMMs, recogClasses))
		{
			return 0;
		}
		int N = (int)vvData.size();

		//recogClasses.resize(N, 0);
		//int I = (int)GMMs.size();
		//int ErrCount = 0;
		//for(int n = 0; n < N; n++)
		//{
		//	//double maxni = GMMs[0].CompLikelihoodPart(vvData[n], vDimConsidered);//GMMs[0].CompLikelihood(vvData[n]);
		//	double maxni = GMMs[0].CompLikelihood(vvData[n]);
		//	int argmaxni = 0;
		//	for(int i = 1; i < I ;i++)
		//	{
		//		//double dbi = GMMs[i].CompLikelihoodPart(vvData[n], vDimConsidered);//GMMs[i].CompLikelihood(vvData[n]);
		//		double dbi = GMMs[i].CompLikelihood(vvData[n]);
		//		if(dbi > maxni)
		//		{
		//			maxni = dbi;
		//			argmaxni = i;
		//		}
		//	}
		//	recogClasses[n] = argmaxni;
		//	//if(argmaxni != 0)
		//	//{
		//	//	if(maxni - GMMs[0].CompLikelihood(vvData[n]) < log(1.0))
		//	//		recogClasses[n] = 0;
		//	//}
		//	
		//	//if(recogClasses[n] != classes[n])
		//	//{
		//	//	ErrCount++;
		//	//}
		//}
		//int nN = N / (WIDTH*HEIGHT);
		//for(int nn = 0; nn < nN; nn++)
		//{
		//	vector<int> recogresn = vector<int>(WIDTH*HEIGHT, 0);
		//	for(int mm = 0; mm < WIDTH*HEIGHT; mm++)
		//	{
		//		recogresn[mm] = recogClasses[nn * WIDTH*HEIGHT + mm];
		//	}
		//	vector<int> recogresn2;
		//	GetLargestConnectedComponents(recogresn, recogresn2);
		//	for(int mm = 0; mm < WIDTH*HEIGHT; mm++)
			//        {
			//                recogClasses[nn * WIDTH*HEIGHT + mm] = recogresn2[mm];
			//        }
		//}
		int nTP = 0;
		int nFP = 0;
		int nFN = 0;
		int nTN = 0;
		int ErrCount = 0;
		for(int n = 0; n < N; n++)
		{
			if(classes[n] == 1 && recogClasses[n] == 1)
			{
				nTP++;
			}
			if(classes[n] != 1 && recogClasses[n] == 1)
			{
				nFP++;
			}
			if(classes[n] == 1 && recogClasses[n] != 1)
			{
				nFN++;
			}
			if(classes[n] != 1 && recogClasses[n] != 1)
			{
				nTN++;
			}
			if(recogClasses[n] != classes[n])
			{
				ErrCount++;
			}	
		}
		confres.resize(4);
		confres[0] = nTP;
		confres[1] = nFP;
		confres[2] = nFN;
		confres[3] = nTN;
		Acc = (double)(N - ErrCount) / (double) N;
		return 1;
	}

	static int Test2ClassesGMMs(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses)
	{
		return RecognizeByGMMs(vvData, GMMs, recogClasses);
	}

	static int Test2ClassesGMMsPart(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses)
	{
		return RecognizeByGMMsPart(vvData, GMMs, recogClasses);
	}
	
	static int Test2ClassesGMMsFull(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses)
	{
		return RecognizeByGMMs(vvData, GMMs, recogClasses);
	}	
	
	static int GMMBiasAdapt(vector< vector<double> >& vvData, matrix<double>& newmub, vector<GMM>& oldGMMs, vector<GMM>& newGMMs)
	{
		//==================================================================================//
		// Get dimension and number of data
		//==================================================================================//
		int N = (int)vvData.size();
		if(N < 1)
		{
			cerr<<"Empty Data vector for adaptation"<<endl;
			return 0;
		}
		int D = (int)vvData[1].size();
		for(int n = 0; n < N; n++)
		{
			int Dn = (int)vvData[n].size();
			if(Dn != D)
			{
				cerr<<"Data components have different dimensions"<<endl;
				return 0;
			}
		}
		int I = (int) oldGMMs.size();
		//==================================================================================//
		// Compute posterior probability and get first recognition result
		//==================================================================================//
		vector< vector< vector<double> > > gammank = vector< vector< vector<double> > >(I);
		for(int i = 0; i < I ; i++)
		{
			gammank[i] = vector< vector<double> >(N, vector<double>(D, 0.0));
			for(int n = 0; n < N; n++)
			{
				oldGMMs[i].CompPosterioriPart(vvData[n], gammank[i][n], vDimConsidered);
			}
		}
		vector<int> modelseq(N, 1);
		for(int n = 0; n < N; n++)
		{
			double dbLik = LZERO;
			for(int i = 0; i < I ; i++)
			{
				double db = oldGMMs[i].CompLikelihoodPart(vvData[n], vDimConsidered);
				if(db > dbLik)
				{
					dbLik = db;
					modelseq[n] = i;
				}
			}
		}
		//==================================================================================//
		// EM core
		//==================================================================================//
		matrix<double> b1 = matrix<double>(D, D);
		matrix<double> b2 = matrix<double>(D, 1);
		b1.Null();
		b2.Null();
		matrix<double> oldb1 = matrix<double>(D, D);
		matrix<double> oldb2 = matrix<double>(D, 1);
		oldb1.Null();
		oldb2.Null();
		for(int n = 0; n < N; n++)
		{
			int modeli = modelseq[n];

			matrix<double> mvec;
			mvec.SetSize(D, 1);
			for(size_t d = 0; d < D; d++)
			{
				mvec(d, 0) = vvData[n][d];
			}

			for(int k = 0; k < oldGMMs[modeli].m_nMixes; k++)
			{
	#ifndef DIAGONAL
				b1 += gammank[modeli][n][k] * oldGMMs[modeli].m_vmdInvCovar[k];
				b2 += gammank[modeli][n][k] * oldGMMs[modeli].m_vmdInvCovar[k] * (mvec - oldGMMs[modeli].m_vmdMean[k]);
	#else
				for(int d1 = 0; d1 < D; d1++)
				{
					b1(d1, d1) += gammank[modeli][n][k] * oldGMMs[modeli].m_vmdInvCovar[k](d1, d1);
					b2(d1, 0) += gammank[modeli][n][k] * oldGMMs[modeli].m_vmdInvCovar[k](d1, d1) * (mvec(d1, 0) - oldGMMs[modeli].m_vmdMean[k](d1, 0));
				}
	#endif
			}
		}
	#ifndef DIAGONAL
		newmub = !b1 * b2;
	#else
		newmub.SetSize(D,1);
		for(int d1 = 0; d1 < D; d1++)
		{
			newmub(d1, 0) = b2(d1, 0)/b1(d1, d1);
		}
	#endif
		//==================================================================================//
		// Update Model, M step
		//==================================================================================//
		newGMMs.resize(I);
		for(int i = 0; i < I; i++)
		{
			newGMMs[i] = GMM(oldGMMs[i]);
			for(int k = 0; k < newGMMs[i].m_nMixes; k++)
			{
				newGMMs[i].m_vmdMean[k] += newmub;
			}
		}
		return 1;
	}

	static int GMMBiasAdapt2(vector< vector<double> >& vvData, matrix<double>& Sigmab, matrix<double>& newmub, matrix<double>& newSigmab, vector<GMM>& oldGMMs, vector<GMM>& newGMMs, int nIter)
	{
		//==================================================================================//
		// Get dimension and number of data
		//==================================================================================//
		int N = (int)vvData.size();
		if(N < 1)
		{
			cerr<<"Empty Data vector for adaptation"<<endl;
			return 0;
		}
		int D = (int)vvData[1].size();
		for(int n = 0; n < N; n++)
		{
			int Dn = (int)vvData[n].size();
			if(Dn != D)
			{
				cerr<<"Data components have different dimensions"<<endl;
				return 0;
			}
		}
		int I = (int) oldGMMs.size();
		//==================================================================================//
		// Compute posterior probability and get first recognition result
		//==================================================================================//
		vector< vector< vector<double> > > gammank = vector< vector< vector<double> > >(I);
		for(int i = 0; i < I ; i++)
		{
			gammank[i] = vector< vector<double> >(N);
			for(int n = 0; n < N; n++)
			{
				 oldGMMs[i].CompPosterioriPart(vvData[n], gammank[i][n], vDimConsidered);//oldGMMs[i].CompPosteriori(vvData[n], gammank[i][n]);
			}
		}
		vector<int> modelseq(N, 0);
		for(int n = 0; n < N; n++)
		{
			double dbLik = 0;
			for(int i = 0; i < I ; i++)
			{
				double db = oldGMMs[i].CompLikelihoodPart(vvData[n], vDimConsidered);//oldGMMs[i].CompLikelihood(vvData[n]);
				if(db > dbLik)
				{
					dbLik = db;
					modelseq[n] = i;
				}
			}
		}
		if(nIter == 0)
		{
			Sigmab.SetSize(D, D);
			Sigmab.Null();
			matrix<double> diagVar;
			CompDiagVar(vvData, diagVar);
			for(int n = 0; n < N; n++)
			{
				for(int k = 0; k < oldGMMs[modelseq[n]].m_nMixes; k++)
				{
					for(int d = 0; d < D; d++)
					{
						double dbnkd = oldGMMs[modelseq[n]].m_vmdCovar[k](d,d);
						double dbnkd1 = diagVar(d,d);
						double dbnkd2 = 0.001 * dbnkd;
						if(dbnkd1 - dbnkd >= 0.001 * dbnkd)
						{
							dbnkd2 = dbnkd1 - dbnkd;
						}
						Sigmab(d,d) += gammank[modelseq[n]][n][k] * dbnkd2 / (double)N;
					}
				}
			}
		}

		//==================================================================================//
		// EM core
		//==================================================================================//
		matrix<double> b1 = matrix<double>(D, D);
		matrix<double> b2 = matrix<double>(D, 1);
		b1.Null();
		b2.Null();
		vector< matrix<double> > vmvec = vector< matrix<double> >(N);
		double N2 = 0.0;
		for(int n = 0; n < N; n++)
		{
			int modeli = modelseq[n];
			vmvec[n].SetSize(D, 1);
			for(size_t d = 0; d < D; d++)
			{
				vmvec[n](d, 0) = vvData[n][d];
			}
			for(int k = 0; k < oldGMMs[modeli].m_nMixes; k++)
			{
				N2 += gammank[modeli][n][k];
				b1 += gammank[modeli][n][k] * oldGMMs[modeli].m_vmdInvCovar[k];
				b2 += gammank[modeli][n][k] * oldGMMs[modeli].m_vmdInvCovar[k] * (vmvec[n] - oldGMMs[modeli].m_vmdMean[k]);
			}
		}
		// Compute posteroir probability
		newmub = !b1 * b2;
		newSigmab.SetSize(D, D);
		newSigmab.Null();
		for(int n = 0; n < N; n++)
		{
			int modeli = modelseq[n];
			for(int k = 0; k < oldGMMs[modeli].m_nMixes; k++)
			{
				matrix<double> Sigman = Sigmab + oldGMMs[modeli].m_vmdCovar[k];
				matrix<double> InvSigman = !Sigman;
				matrix<double> Eb = newmub + Sigmab * InvSigman * (vmvec[n] - newmub - oldGMMs[modeli].m_vmdMean[k]);
				matrix<double> Eb2 = Sigmab *  InvSigman * oldGMMs[modeli].m_vmdCovar[k] + Eb * ~Eb;
				newSigmab += gammank[modeli][n][k] / N2 * Eb2;
			}
			//Eb = mubp + (Sigmab * inv(Sigmab + sk)*(X(n,:) - mubp - gmms{i}.mu(k,:))')';
			//Eb2 = Sigmab * inv(Sigmab + sk) * sk + Eb' * Eb;
			//Sigmab2 = Sigmab2 + gammank{i}(n, k) * Eb2;

		}
		newSigmab -= newmub * ~newmub;
		//==================================================================================//
		// Update Model
		//==================================================================================//
		newGMMs.resize(I);
		for(int i = 0; i < I; i++)
		{
			newGMMs[i] = GMM(oldGMMs[i]);
			for(int k = 0; k < newGMMs[i].m_nMixes; k++)
			{
				newGMMs[i].m_vmdMean[k] += newmub;
				if(newGMMs[i].m_szCovType == "diagonal")
				{
					for(int d = 0; d < D; d++)
					{
						if(newGMMs[i].m_vmdCovar[k](d,d) + newSigmab(d,d) > 0)
						{
							newGMMs[i].m_vmdCovar[k](d,d) += newSigmab(d,d);
						}
						else
						{
							newGMMs[i].m_vmdCovar[k](d,d) += 0;
						}
					}
				}
				else
				{
					newGMMs[i].m_vmdCovar[k] += newSigmab;
				}
			}
		}
		return 1;
	}

	static int GenEpsilon(matrix<double>& mat, int nDims, double eps)
	{
		mat.SetSize(nDims, nDims);
		mat.Null();
		for(int d = 0; d < nDims; d++)
		{
			mat(d, d) = eps;
		}
		return 1;
	}

	static int GMMAdaptRcogbyFrame(vector<GMM>& oldGMMs, vector< vector<double> >& vvData, vector<int>& vclasses, double& Acc, int nFrameSize)
	{
		vector<int> vrecogclasses(vclasses.size(), 0);
		int N = (int)vvData.size();
		int nN = (int) ceil((double)N / (double)nFrameSize );
		int nErrCount = 0;
		for(int n = 0; n < nN; n++)
		{
			int nStart = n * nFrameSize;
			int nStop = min((n + 1) * nFrameSize - 1, N - 1);
			if(nStart > nStop)
			{
				break;
			}
			vector< vector<double> > vvSub = vector< vector<double> >(nStop-nStart+1, vector<double>(vvData[0].size(), 0.0));
			for(int i = nStart; i <= nStop; i++)
			{
				vvSub[i-nStart] = vvData[i];
			}
			//matrix<double> newmub;
			//vector<GMM> newGMMs;
			//GMMAdaptTool::GMMBiasAdapt(vvSub, newmub, oldGMMs, newGMMs);
			vector<int> recogClassesn;
			//GMMAdaptTool::Test2ClassesGMMs(vvSub, newGMMs, recogClassesn);
			vector<int> recogClassesn2;
			GMMAdaptTool::GMMAdaptRcogAFrame(oldGMMs, vvSub, recogClassesn);
			GMMAdaptTool::GetLargestConnectedComponents(recogClassesn, recogClassesn2);
			for(int i = nStart; i <= nStop; i++)
			{
				vrecogclasses[i] = recogClassesn2[i-nStart];
				if(vrecogclasses[i] != vclasses[i])
				{
					nErrCount++;
				}
			}
		}
		Acc = (double)(N - nErrCount)/(double)N;
		return 1;
	}

	static int GMMAdaptRcogbyFrame2(vector<GMM>& oldGMMs, vector< vector<double> >& vvData, vector<int>& vclasses, double& Acc, int nFrameSize)
	{
		vector<int> vrecogclasses(vclasses.size(), 0);
		int N = (int)vvData.size();
		int nN = (int) ceil((double)N / (double)nFrameSize );
		int nErrCount = 0;
		for(int n = 0; n < nN; n++)
		{
			int nStart = n * nFrameSize;
			int nStop = min((n + 1) * nFrameSize - 1, N - 1);
			if(nStart > nStop)
			{
				break;
			}
			vector< vector<double> > vvSub = vector< vector<double> >(nStop-nStart+1, vector<double>(vvData[0].size(), 0.0));
			for(int i = nStart; i <= nStop; i++)
			{
				vvSub[i-nStart] = vvData[i];
			}
			matrix<double> newmub;
			vector< vector<GMM> > vnewGMMs;
			vnewGMMs.push_back(oldGMMs);
			vector< matrix<double> > vnewSigmab;
			vector< matrix<double> > vnewmub;
			matrix<double> Sigmab;
			// init Sigmab here
			GMMAdaptTool::GenEpsilon(Sigmab, vvData[0].size(), 1);
			vnewSigmab.push_back(Sigmab);
			
			for(int itr = 0; itr < NUMITERATION; itr++)
			{
				vector<GMM> newGMMs;
				matrix<double> newmub;
				matrix<double> newSigmab;
				GMMAdaptTool::GMMBiasAdapt2(vvSub, vnewSigmab[itr], newmub, newSigmab, vnewGMMs[itr], newGMMs, itr);
				vnewmub.push_back(newmub);
				vnewSigmab.push_back(newSigmab);

				//GMMBiasAdapt2(vvSub, newmub, vnewGMMs[itr], newGMMs);
				vnewGMMs.push_back(newGMMs);
			}
			vector<int> recogClassesn;
			GMMAdaptTool::Test2ClassesGMMs(vvSub, vnewGMMs[vnewGMMs.size() - 1], recogClassesn);
			vector<int> recogClassesn2;
			GMMAdaptTool::GetLargestConnectedComponents(recogClassesn, recogClassesn2);
			for(int i = nStart; i <= nStop; i++)
			{
				vrecogclasses[i] = recogClassesn2[i-nStart];
				if(vrecogclasses[i] != vclasses[i])
				{
					nErrCount++;
				}
			}
		}
		Acc = (double)(N - nErrCount)/(double)N;
		return 1;
	}

	static int GMMAdaptRcogbyFrame(string GMMPaths, string DataPaths, double& Acc, int nFrameSize)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		vector< vector<double> > vvData;
		vector<int> vClasses;
		GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
		GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
		GMMAdaptTool::LoadLogData(DataPaths, vvData, vClasses);
		return GMMAdaptTool::GMMAdaptRcogbyFrame(GMMs, vvData, vClasses, Acc, nFrameSize);
	}

	static int GMMAdaptRcogbyFrame2(string GMMPaths, string DataPaths, double& Acc, int nFrameSize)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		vector< vector<double> > vvData;
		vector<int> vClasses;
		GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
		GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
		GMMAdaptTool::LoadLogData(DataPaths, vvData, vClasses);
		return GMMAdaptTool::GMMAdaptRcogbyFrame2(GMMs, vvData, vClasses, Acc, nFrameSize);
	}

	// Stochastic matching on feature space
	static int GMMAdaptRcogAFrame(vector<GMM>& oldGMMs, vector< vector<double> >& vvData, vector<int>& vrecogclasses)
	{
		vrecogclasses.resize(vvData.size());
		int nFrameSize = (int)vvData.size();
		int nErrCount = 0;
		int n = 0;
		int nStart = 0;
		int nStop = nFrameSize - 1;
		matrix<double> newmub;
		vector<GMM> newGMMs;
		GMMAdaptTool::GMMBiasAdapt(vvData, newmub, oldGMMs, newGMMs);
		GMMAdaptTool::Test2ClassesGMMs(vvData, newGMMs, vrecogclasses);
		// Doing grouping elimination, pick the largest one
		return 1;
	}
	
	static void AutoLabelRGBOnly(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio, int& nHandModelIndex)
	{
		vector< vector<double> > vvData;
		int yrange = ExtLowLevelImgFeature::ExtColorFeatureRGBOnly(szPath, vvData, shrinkratio);

		//DWORD dw1 = 0, dw2 = 0;
		//IplImage* pimg;
		//pimg = cvLoadImage(szPath.c_str());
		////IplImage* pfimg = cvCreateImage(cvSize(pimg->width, pimg->height), IPL_DEPTH_32F, pimg->nChannels);
		//cv::Mat img = cv::Mat(pimg);
		////cvConvert(pimg, pfimg);
		//cv::Mat fimg; // = cv::Mat(pfimg);
		//img.convertTo(fimg, CV_32FC3);
		//
		//int nratio = (int) shrinkratio;
		//double dbBlkNum = (double)(img.cols * img.rows/ shrinkratio);
		//double dbLLR = 1.0;
		//int xrange = img.cols / shrinkratio;
		//int yrange = img.rows / shrinkratio;
		//mask = cv::Mat(yrange, xrange, CV_8UC1, cv::Scalar(0, 0, 0, 0));
		//DWORD dwcol = 0; 
		//DWORD dwrow = 0; 
		//DWORD dwcvt = 0;
		//DWORD dwmv = 0; 
		//DWORD dwfet = 0; 
		//DWORD dwlik = 0; 
		//vector< vector<double> > vvData;

		//for(int x = 0; x < xrange; x++)
		//{
		//	cv::Scalar fmeanxy, stddevxy;
		//	cv::Scalar fmeanxyhsv, stddevxyhsv;
		//	cv::Scalar fmeanxyx, stddevxyx;
		//	cv::Scalar fmeanxyy, stddevxyy;
		//	cv::Mat colx = fimg.colRange(x*shrinkratio, (x+1)*shrinkratio); 
		//	// Compute color feature, note that 
		//	// we only need some part of the feature
		//	// but the current code is not flexible enough to just compute them
		//	// so this may need to be changed to enhance speed
		//	for(int y = 0; y < yrange; y++)
		//	{	
		//		cv::Mat sub = colx.rowRange(y*shrinkratio, (y+1)*shrinkratio);
		//		cv::Mat sub8bits;
		//		sub.convertTo(sub8bits, CV_8UC3);
		//		cv::Mat subhsv;
		//		cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
		//		cv::meanStdDev(sub, fmeanxy, stddevxy);
		//		cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);
		//		double dbBMean = fmeanxy.val[0];
		//		double dbGMean = fmeanxy.val[1];
		//		double dbRMean = fmeanxy.val[2];
		//		vector<double> vdbTemp;
		//		vdbTemp.push_back(dbBMean);
		//		vdbTemp.push_back(dbGMean);
		//		vdbTemp.push_back(dbRMean);
		//		vvData.push_back(vdbTemp);
		//	}
		//}

		vector<int> vrecogclasses;
		// classifying skin/nonskin for thes blocks
		vrecogclasses.resize(vvData.size());
		Test2ClassesGMMsFull(vvData, vGMMs, vrecogclasses);
		size_t nResultsSize = vrecogclasses.size();
		mask = cv::Mat(yrange, nResultsSize/yrange, CV_8UC1);
		for(size_t nn = 0; nn < nResultsSize; nn++)
		{
			int x = nn / yrange;
			int y = nn - x * yrange;
			uchar* uptr = mask.ptr<uchar>(y);
			if(vrecogclasses[nn] == nHandModelIndex)
			{
				uptr[x] = 255;
			}
			else
			{
				uptr[x] = 0;
			}
		}
		//cvReleaseImage(&pimg);
		//cvReleaseImage(&pfimg);
	}



	static void AutoLabel_BS0(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio, int& nHandModelIndex, cv::Mat& ref)
	{
		vector< vector<double> > vvData;
		//ExtLowLevelImgFeature::ExtColorFeatures_BS0(IplImage*& pimg, vector< vector<double> >& vvFeatures, int BLKSIZE, IplImage*& pRefImg, vector< vector<double> >& vvCorrs);
		IplImage iplref = ref.operator IplImage();
		IplImage* pref = &iplref;
		vector< vector<double> > vvCorrs;
		vector<double> vKSs;
		//int yrange = ExtLowLevelImgFeature::ExtColorFeatures_BS0(szPath, vvData, shrinkratio, pref, vvCorrs);
		int yrange = ExtLowLevelImgFeature::ExtColorFeatures_BS0(szPath, vvData, shrinkratio, pref, vvCorrs, vKSs);

		//DWORD dw1 = 0, dw2 = 0;
		//IplImage* pimg;
		//pimg = cvLoadImage(szPath.c_str());
		////IplImage* pfimg = cvCreateImage(cvSize(pimg->width, pimg->height), IPL_DEPTH_32F, pimg->nChannels);
		//cv::Mat img = cv::Mat(pimg);
		////cvConvert(pimg, pfimg);
		//cv::Mat fimg; // = cv::Mat(pfimg);
		//img.convertTo(fimg, CV_32FC3);
		//
		//int nratio = (int) shrinkratio;
		//double dbBlkNum = (double)(img.cols * img.rows/ shrinkratio);
		//double dbLLR = 1.0;
		//int xrange = img.cols / shrinkratio;
		//int yrange = img.rows / shrinkratio;
		//mask = cv::Mat(yrange, xrange, CV_8UC1, cv::Scalar(0, 0, 0, 0));
		//DWORD dwcol = 0; 
		//DWORD dwrow = 0; 
		//DWORD dwcvt = 0;
		//DWORD dwmv = 0; 
		//DWORD dwfet = 0; 
		//DWORD dwlik = 0; 
		//vector< vector<double> > vvData;

		//for(int x = 0; x < xrange; x++)
		//{
		//	cv::Scalar fmeanxy, stddevxy;
		//	cv::Scalar fmeanxyhsv, stddevxyhsv;
		//	cv::Scalar fmeanxyx, stddevxyx;
		//	cv::Scalar fmeanxyy, stddevxyy;
		//	cv::Mat colx = fimg.colRange(x*shrinkratio, (x+1)*shrinkratio); 
		//	// Compute color feature, note that 
		//	// we only need some part of the feature
		//	// but the current code is not flexible enough to just compute them
		//	// so this may need to be changed to enhance speed
		//	for(int y = 0; y < yrange; y++)
		//	{	
		//		cv::Mat sub = colx.rowRange(y*shrinkratio, (y+1)*shrinkratio);
		//		cv::Mat sub8bits;
		//		sub.convertTo(sub8bits, CV_8UC3);
		//		cv::Mat subhsv;
		//		cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
		//		cv::meanStdDev(sub, fmeanxy, stddevxy);
		//		cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);
		//		double dbBMean = fmeanxy.val[0];
		//		double dbGMean = fmeanxy.val[1];
		//		double dbRMean = fmeanxy.val[2];
		//		vector<double> vdbTemp;
		//		vdbTemp.push_back(dbBMean);
		//		vdbTemp.push_back(dbGMean);
		//		vdbTemp.push_back(dbRMean);
		//		vvData.push_back(vdbTemp);
		//	}
		//}

		vector<int> vrecogclasses;
		// classifying skin/nonskin for thes blocks
		vrecogclasses.resize(vvData.size());
		Test2ClassesGMMsFull(vvData, vGMMs, vrecogclasses);
		size_t nResultsSize = vrecogclasses.size();
		mask = cv::Mat(yrange, nResultsSize/yrange, CV_8UC1);
		for(size_t nn = 0; nn < nResultsSize; nn++)
		{
			int x = nn / yrange;
			int y = nn - x * yrange;
			uchar* uptr = mask.ptr<uchar>(y);
			double Corr = 0.1140 * vvCorrs[nn][0] + 0.5870 * vvCorrs[nn][1] + 0.2989 * vvCorrs[nn][2];
			if(vrecogclasses[nn] == nHandModelIndex && Corr < CCThreshold)
			{
				uptr[x] = 255;
			}
			else
			{
				uptr[x] = 0;
			}
		}
		//cvReleaseImage(&pimg);
		//cvReleaseImage(&pfimg);
	}

	static void AutoLabelHTKMAP_BS0_Fusion(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
					int x = nn / yrange;
					int y = nn - x * yrange;
					uchar* uptr = vMasks[n].ptr<uchar>(y);
					if(vvRecogClasses[n][nn] == nHandModelIndex)
					{
							uptr[x] = 255;
					}
					else
					{
							uptr[x] = 0;
					}
			}
		}

	}

	static void AutoLabelHTKMAP_BS0_Fusion_with_Scores(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq, vector<cv::Mat>& vScore0, vector<cv::Mat>& vScore1, vector<cv::Mat>& vLLR, vector<cv::Mat>& vCorr, vector<cv::Mat>& vKS)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
	#ifdef _MSC_VER
		string mkdirFolder = "mkdir " + szFEFolder;
	#elif defined(__GNUC__)
		string mkdirFolder = "mkdir -p " + szFEFolder;
	#endif
		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		vector< vector<double> > vvScore0;
		vector< vector<double> > vvScore1;
		vector< vector< vector<double> > > vvvScoreForFusions;
		//HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		//HTKGMMMAP_BS0_Fusion_With_Scores(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1);
		HTKGMMMAP_BS0_Fusion_With_Scores(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1, vvvScoreForFusions);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		vScore0.resize(N);
		vScore1.resize(N);
		vLLR.resize(N);
		vCorr.resize(N);
		vKS.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			// Scores0
			vScore0[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Scores1
			vScore1[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Decision
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			vLLR[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			vCorr[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			vKS[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
					int x = nn / yrange;
					int y = nn - x * yrange;
					uchar* uptr = vMasks[n].ptr<uchar>(y);
					double* ptr0 = vScore0[n].ptr<double>(y);
					double* ptr1 = vScore1[n].ptr<double>(y);
					double* ptrLLR = vLLR[n].ptr<double>(y);
					double* ptrCorr = vCorr[n].ptr<double>(y);
					double* ptrKS = vKS[n].ptr<double>(y);
					ptr0[x] = vvScore0[n][nn];
					ptr1[x] = vvScore1[n][nn];
					//cout<<"["<<n<<","<<nn<<"]:"<<vvvScoreForFusions[n][nn][0]<<","<<vvvScoreForFusions[n][nn][1]<<endl;
					ptrLLR[x] = vvvScoreForFusions[n][nn][0];
					ptrCorr[x] = vvvScoreForFusions[n][nn][1];
					ptrKS[x] = vvvScoreForFusions[n][nn][2];
					if(vvRecogClasses[n][nn] == nHandModelIndex)
					{
							uptr[x] = 255;
					}
					else
					{
							uptr[x] = 0;
					}
			}
		}

	}


	static void AutoLabelHTKMAP_BS0_Fusion_with_Scores(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq, vector<cv::Mat>& vScore0, vector<cv::Mat>& vScore1, vector<cv::Mat>& vLLR, vector<cv::Mat>& vCorr)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
	#ifdef _MSC_VER
		string mkdirFolder = "mkdir " + szFEFolder;
	#elif defined(__GNUC__)
		string mkdirFolder = "mkdir -p " + szFEFolder;
	#endif
		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		vector< vector<double> > vvScore0;
		vector< vector<double> > vvScore1;
		vector< vector< vector<double> > > vvvScoreForFusions;
		//HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		//HTKGMMMAP_BS0_Fusion_With_Scores(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1);
		HTKGMMMAP_BS0_Fusion_With_Scores(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1, vvvScoreForFusions);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		vScore0.resize(N);
		vScore1.resize(N);
		vLLR.resize(N);
		vCorr.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			// Scores0
			vScore0[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Scores1
			vScore1[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Decision
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			vLLR[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			vCorr[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
					int x = nn / yrange;
					int y = nn - x * yrange;
					uchar* uptr = vMasks[n].ptr<uchar>(y);
					double* ptr0 = vScore0[n].ptr<double>(y);
					double* ptr1 = vScore1[n].ptr<double>(y);
					double* ptrLLR = vLLR[n].ptr<double>(y);
					double* ptrCorr = vCorr[n].ptr<double>(y);
					ptr0[x] = vvScore0[n][nn];
					ptr1[x] = vvScore1[n][nn];
					//cout<<"["<<n<<","<<nn<<"]:"<<vvvScoreForFusions[n][nn][0]<<","<<vvvScoreForFusions[n][nn][1]<<endl;
					ptrLLR[x] = vvvScoreForFusions[n][nn][0];
					ptrCorr[x] = vvvScoreForFusions[n][nn][1];
					if(vvRecogClasses[n][nn] == nHandModelIndex)
					{
							uptr[x] = 255;
					}
					else
					{
							uptr[x] = 0;
					}
			}
		}

	}

	static void AutoLabelHTKMAP_BS0_Fusion_with_Scores(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq, vector<cv::Mat>& vScore0, vector<cv::Mat>& vScore1)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		vector< vector<double> > vvScore0;
		vector< vector<double> > vvScore1;
		vector< vector< vector<double> > > vvvRecogFusionScores;
		//HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		HTKGMMMAP_BS0_Fusion_With_Scores(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1);
		
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		vScore0.resize(N);
		vScore1.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			// Scores0
			vScore0[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Scores1
			vScore1[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Decision
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
					int x = nn / yrange;
					int y = nn - x * yrange;
					uchar* uptr = vMasks[n].ptr<uchar>(y);
					double* ptr0 = vScore0[n].ptr<double>(y);
					double* ptr1 = vScore1[n].ptr<double>(y);
					ptr0[x] = vvScore0[n][nn];
					ptr1[x] = vvScore1[n][nn];
					if(vvRecogClasses[n][nn] == nHandModelIndex)
					{
							uptr[x] = 255;
					}
					else
					{
							uptr[x] = 0;
					}
			}
		}

	}

	static void Approx_AutoLabel_BS0(string szPath, bool*** vvvDecisionTable, cv::Mat& mask, int shrinkratio, int& nHandModelIndex, cv::Mat& ref)
	{
		vector< vector<double> > vvData;
		IplImage iplref = ref.operator IplImage();
		IplImage* pref = &iplref;
		vector< vector<double> > vvCorrs;
		vector<double> vKSs;
		//int yrange = ExtLowLevelImgFeature::ExtColorFeatures_BS0(szPath, vvData, shrinkratio, pref, vvCorrs);
		int yrange = ExtLowLevelImgFeature::ExtColorFeatures_BS0(szPath, vvData, shrinkratio, pref, vvCorrs, vKSs);

		vector<int> vrecogclasses;
		// classifying skin/nonskin for thes blocks
		vrecogclasses.resize(vvData.size());
		ApproxRecognizeByTable(vvData, vvvDecisionTable, vrecogclasses);
		size_t nResultsSize = vrecogclasses.size();
		mask = cv::Mat(yrange, nResultsSize/yrange, CV_8UC1);
		for(size_t nn = 0; nn < nResultsSize; nn++)
		{
			int x = nn / yrange;
			int y = nn - x * yrange;
			uchar* uptr = mask.ptr<uchar>(y);
			double Corr = 0.1140 * vvCorrs[nn][0] + 0.5870 * vvCorrs[nn][1] + 0.2989 * vvCorrs[nn][2];
			if(vrecogclasses[nn] == nHandModelIndex && Corr < CCThreshold)
			{
				uptr[x] = 255;
			}
			else
			{
				uptr[x] = 0;
			}
		}
		//cvReleaseImage(&pimg);
		//cvReleaseImage(&pfimg);

	}

	static void AutoLabelHTKMAPPerFrame(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

               // struct stat st;
                if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		HTKGMMAdaptPerFrame(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
				int x = nn / yrange;
				int y = nn - x * yrange;
				uchar* uptr = vMasks[n].ptr<uchar>(y);
				if(vvRecogClasses[n][nn] == nHandModelIndex)
				{
					uptr[x] = 255;
				}
				else
				{
					uptr[x] = 0;
				}
			}
		}
	}

	static void AutoLabelHTKMAPASeq(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		HTKGMMAdaptASeq(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
				int x = nn / yrange;
				int y = nn - x * yrange;
				uchar* uptr = vMasks[n].ptr<uchar>(y);
				if(vvRecogClasses[n][nn] == nHandModelIndex)
				{
					uptr[x] = 255;
				}
				else
				{
					uptr[x] = 0;
				}
			}
		}
	}
	static void AutoLabelSMAPPerFrame(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		HTKGMMSMAPPerFrame(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
				int x = nn / yrange;
				int y = nn - x * yrange;
				uchar* uptr = vMasks[n].ptr<uchar>(y);
				if(vvRecogClasses[n][nn] == nHandModelIndex)
				{
					uptr[x] = 255;
				}
				else
				{
					uptr[x] = 0;
				}
			}
		}
	}

	static void AutoLabelSMAPPerFrame_BS0(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath, cv::Mat& ref)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		//HTKGMMSMAPPerFrame_BS0(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vvvCorrs);
		HTKGMMSMAPPerFrame_BS0(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vvvCorrs, vvKSs);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
				int x = nn / yrange;
				int y = nn - x * yrange;
				uchar* uptr = vMasks[n].ptr<uchar>(y);
				if(vvRecogClasses[n][nn] == nHandModelIndex)
				{
						uptr[x] = 255;
				}
				else
				{
						uptr[x] = 0;
				}
			}
		}
	}


	static void AutoLabelSMAPASeq(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		HTKGMMSMAPASeq(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
				int x = nn / yrange;
				int y = nn - x * yrange;
				uchar* uptr = vMasks[n].ptr<uchar>(y);
				if(vvRecogClasses[n][nn] == nHandModelIndex)
				{
					uptr[x] = 255;
				}
				else
				{
					uptr[x] = 0;
				}
			}
		}
	}

	static void AutoLabelSMAPASeq_BS0(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath, cv::Mat& ref)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		//HTKGMMSMAPASeq_BS0(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vvvCorrs);
		HTKGMMSMAPASeq_BS0(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vvvCorrs, vvKSs);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
					int x = nn / yrange;
					int y = nn - x * yrange;
					uchar* uptr = vMasks[n].ptr<uchar>(y);
					if(vvRecogClasses[n][nn] == nHandModelIndex)
					{
							uptr[x] = 255;
					}
					else
					{
							uptr[x] = 0;
					}
			}
		}
	}



	// Get skin/nonskin color label, only consider RGB mean
	static void AutoLabelPart(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio, int& nHandModelIndex)
	{
		vector< vector<double> > vvData;
		int yrange = ExtLowLevelImgFeature::ExtColorFeatureRGBOnly(szPath, vvData, shrinkratio);

		//DWORD dw1 = 0, dw2 = 0;
		//IplImage* pimg;
		//pimg = cvLoadImage(szPath.c_str());
		//cv::Mat img(pimg);
		//cv::Mat fimg;
		//img.convertTo(fimg, CV_32FC3);
		//int nratio = (int) shrinkratio;
		//double dbBlkNum = (double)(img.cols * img.rows/ shrinkratio);
		//double dbLLR = 1.0;
		//int xrange = img.cols / shrinkratio;
		//int yrange = img.rows / shrinkratio;
		//mask = cv::Mat(yrange, xrange, CV_8UC1, cv::Scalar(0, 0, 0, 0));
		//DWORD dwcol = 0; 
		//DWORD dwrow = 0; 
		//DWORD dwcvt = 0;
		//DWORD dwmv = 0; 
		//DWORD dwfet = 0; 
		//DWORD dwlik = 0; 
		//vector< vector<double> > vvData;

		//for(int x = 0; x < xrange; x++)
		//{
		//	cv::Scalar fmeanxy, stddevxy;
		//	cv::Scalar fmeanxyhsv, stddevxyhsv;
		//	cv::Scalar fmeanxyx, stddevxyx;
		//	cv::Scalar fmeanxyy, stddevxyy;
		//	cv::Mat colx = fimg.colRange(x*shrinkratio, (x+1)*shrinkratio); 
		//	// Compute color feature, note that 
		//	// we only need some part of the feature
		//	// but the current code is not flexible enough to just compute them
		//	// so this may need to be changed to enhance speed
		//	for(int y = 0; y < yrange; y++)
		//	{	
		//		cv::Mat sub = colx.rowRange(y*shrinkratio, (y+1)*shrinkratio);
		//		cv::Mat sub8bits;
		//		sub.convertTo(sub8bits, CV_8UC3);
		//		cv::Mat subhsv;
		//		cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
		//		cv::meanStdDev(sub, fmeanxy, stddevxy);
		//		cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);
		//		double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
		//		double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
		//		double dbBMean = fmeanxy.val[0];
		//		double dbGMean = fmeanxy.val[1];
		//		double dbRMean = fmeanxy.val[2];
		//		double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
		//		double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
		//		double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
		//		double dbVMean = fmeanxyhsv.val[2];
		//		double dbSMean = fmeanxyhsv.val[1];
		//		double dbHMean = fmeanxyhsv.val[0];
		//		vector<double> vdbTemp;
		//		vdbTemp.push_back(dbAvgStdev);
		//		vdbTemp.push_back(dbIMean);
		//		vdbTemp.push_back(dbBMean);
		//		vdbTemp.push_back(dbGMean);
		//		vdbTemp.push_back(dbRMean);
		//		vdbTemp.push_back(dbExBMean);
		//		vdbTemp.push_back(dbExGMean);
		//		vdbTemp.push_back(dbExRMean);
		//		vdbTemp.push_back(dbVMean);
		//		vdbTemp.push_back(dbSMean);
		//		vdbTemp.push_back(dbHMean);
		//		vvData.push_back(vdbTemp);
		//	}
		//}
		vector<int> vrecogclasses;
		// classifying skin/nonskin for thes blocks
		vrecogclasses.resize(vvData.size());
		Test2ClassesGMMsPart(vvData, vGMMs, vrecogclasses);
		size_t nResultsSize = vrecogclasses.size();
		for(size_t nn = 0; nn < nResultsSize; nn++)
		{
			int x = nn / yrange;
			int y = nn - x * yrange;
			uchar* uptr = mask.ptr<uchar>(y);
			if(vrecogclasses[nn] == nHandModelIndex)
			{
				uptr[x] = 255;
			}
			else
			{
				uptr[x] = 0;
			}
		}
		//cvReleaseImage(&pimg);
	}

	// Get skin/nonskin color label with stochastic feature matching
	static void AutoLabelPartSM(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int shrinkratio, int& nHandModelIndex)
	{
		vector< vector<double> > vvData;
		int yrange = ExtLowLevelImgFeature::ExtColorFeatureRGBOnly(szPath, vvData, shrinkratio);

		//DWORD dw1 = 0, dw2 = 0;
		//IplImage* pimg;
		//pimg = cvLoadImage(szPath.c_str());
		//cv::Mat img(pimg);
		//cv::Mat fimg;
		//img.convertTo(fimg, CV_32FC3);
		//int nratio = (int) shrinkratio;
		//double dbBlkNum = (double)(img.cols * img.rows/ shrinkratio);
		//double dbLLR = 1.0;
		//int xrange = img.cols / shrinkratio;
		//int yrange = img.rows / shrinkratio;
		//mask = cv::Mat(yrange, xrange, CV_8UC1, cv::Scalar(0, 0, 0, 0));
		//DWORD dwcol = 0; 
		//DWORD dwrow = 0; 
		//DWORD dwcvt = 0;
		//DWORD dwmv = 0; 
		//DWORD dwfet = 0; 
		//DWORD dwlik = 0; 
		//vector< vector<double> > vvData;

		//for(int x = 0; x < xrange; x++)
		//{
		//	cv::Scalar fmeanxy, stddevxy;
		//	cv::Scalar fmeanxyhsv, stddevxyhsv;
		//	cv::Scalar fmeanxyx, stddevxyx;
		//	cv::Scalar fmeanxyy, stddevxyy;
		//	cv::Mat colx = fimg.colRange(x*shrinkratio, (x+1)*shrinkratio); 
		//	// Compute color feature, note that 
		//	// we only need some part of the feature
		//	// but the current code is not flexible enough to just compute them
		//	// so this may need to be changed to enhance speed
		//	for(int y = 0; y < yrange; y++)
		//	{	
		//		cv::Mat sub = colx.rowRange(y*shrinkratio, (y+1)*shrinkratio);
		//		cv::Mat sub8bits;
		//		sub.convertTo(sub8bits, CV_8UC3);
		//		cv::Mat subhsv;
		//		cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
		//		cv::meanStdDev(sub, fmeanxy, stddevxy);
		//		cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);
		//		double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
		//		double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
		//		double dbBMean = fmeanxy.val[0];
		//		double dbGMean = fmeanxy.val[1];
		//		double dbRMean = fmeanxy.val[2];
		//		double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
		//		double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
		//		double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
		//		double dbVMean = fmeanxyhsv.val[2];
		//		double dbSMean = fmeanxyhsv.val[1];
		//		double dbHMean = fmeanxyhsv.val[0];
		//		vector<double> vdbTemp;
		//		vdbTemp.push_back(dbAvgStdev);
		//		vdbTemp.push_back(dbIMean);
		//		vdbTemp.push_back(dbBMean);
		//		vdbTemp.push_back(dbGMean);
		//		vdbTemp.push_back(dbRMean);
		//		vdbTemp.push_back(dbExBMean);
		//		vdbTemp.push_back(dbExGMean);
		//		vdbTemp.push_back(dbExRMean);
		//		vdbTemp.push_back(dbVMean);
		//		vdbTemp.push_back(dbSMean);
		//		vdbTemp.push_back(dbHMean);
		//		vvData.push_back(vdbTemp);
		//	}
		//}

		vector<int> vrecogclasses0;
		vector<int> vrecogclasses;
		// classifying skin/nonskin for thes blocks
		GMMAdaptRcogAFrame(vGMMs, vvData, vrecogclasses);
		size_t nResultsSize = vrecogclasses.size();
		for(size_t nn = 0; nn < nResultsSize; nn++)
		{
			int x = nn / yrange;
			int y = nn - x * yrange;
			uchar* uptr = mask.ptr<uchar>(y);
			if(vrecogclasses[nn] == nHandModelIndex)
			{
				uptr[x] = 255;
			}
			else
			{
				uptr[x] = 0;
			}
		}
		//cvReleaseImage(&pimg);
	}
	
	static int GMMHTKMAPBySCP(string GMMPaths, string DataSCPPaths, double& Acc)
	{
		unordered_map<string, GMM> hGMMs;
			vector<GMM> GMMs;
			GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
			GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
		
		char caBuf1[4096];
		char caBuf2[4096];
		memset(caBuf1, 0, 4096);
		memset(caBuf2, 0, 4096);
		ifstream fin;
		fin.open(DataSCPPaths.c_str());
		string szTmp = DataSCPPaths + "_Tmp";
		string szTmpMLF = DataSCPPaths + "_Tmp.mlf";
		//string HERest = "/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
		//string config = "/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
		//string words = "/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
		//string GMMCopy = GMMPaths + "_Tmp";
		//string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
		string GMMCopy = GMMPaths + "_Tmp";
	#ifdef _MSC_VER
		string HERest = "C:\\debug_hand_detection\\tools\\HTK\\bin.win32\\HERest.exe"; //"C:\\CollectedDataCopy\\recognizer\\HTK\\bin.win32\\HERest.exe";
		string config = "C:\\debug_hand_detection\\config\\config"; //"C:\\Transfer\\AdjustGMMNew\\config";
		string words = "C:\\debug_hand_detection\\config\\words"; //"C:\\Transfer\\AdjustGMMNew\\words";
		string szCopyGMM = "copy " + GMMPaths + " " + GMMCopy;
	#elif defined(__GNUC__)
		string HERest = "HERest"; //"/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
		string config = "/nv/pec1/ycheng31/local/include/AdjustGMMNew/config"; //"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
		string words = "/nv/pec1/ycheng31/local/include/AdjustGMMNew/words"; //"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
		string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
	#endif
		Acc = 0.0;
		double Acc0 = 0.0;
		double N_1 = 0;
		if(fin.is_open())
		{
			while(fin.getline(caBuf2, 4096) && fin.getline(caBuf1, 4096))
			{
				
							string szHand(caBuf1);
							string szBkGnd(caBuf2);
				if(szHand[szHand.length() - 1] == '\r' || szHand[szHand.length() - 1] == '\n')
				{
					caBuf1[szHand.length() - 1] = '\0';
					szHand = caBuf1;
				}
							if(szBkGnd[szBkGnd.length() - 1] == '\r' || szBkGnd[szBkGnd.length() - 1] == '\n')
							{
									caBuf2[szBkGnd.length() - 1] = '\0';
					szBkGnd = caBuf2;
							}
				ofstream fout;
				fout.open(szTmp.c_str());
				if(!fout.is_open())
				{
					cerr<<"Cannot Open Temp File so stop it"<<endl;
					return -1;
				}
				fout<<caBuf1<<endl<<caBuf2<<endl;
				fout.close();
				ofstream fout02;
				fout02.open(szTmpMLF.c_str());
				if(!fout02.is_open())
				{
					cerr<<"Cannot open mlf file for written!"<<endl;
					return -1;
				}
				fout02<<"#!MLF!#"<<endl;
							vector< vector<double> > vvDataHand;
							vector< vector<double> > vvDataBkGnd;			
							GMMAdaptTool::ReadHTKMFCCFeature(vvDataHand, szHand);
							GMMAdaptTool::ReadHTKMFCCFeature(vvDataBkGnd, szBkGnd);
				if(vvDataHand.size())
				{
					string szOutFileName1;
					int nPos001 = szHand.find_last_of("/");
					int nPos002 = szHand.find_last_of(".");
					szOutFileName1 = szHand.substr(nPos001 + 1, nPos002 - nPos001 - 1);
					fout02<<"\"*/"<<szOutFileName1<<".lab\""<<endl;
					fout02<<"Hand\n.\n";
				}
				if(vvDataBkGnd.size())
				{
									string szOutFileName2;
									int nPos001 = szBkGnd.find_last_of("/");
									int nPos002 = szBkGnd.find_last_of(".");
									szOutFileName2 = szBkGnd.substr(nPos001 + 1, nPos002 - nPos001 - 1);
									fout02<<"\"*/"<<szOutFileName2<<".lab\""<<endl;
									fout02<<"BkGnd\n.\n";			
				}
				fout02.close();
				system(szCopyGMM.c_str());
				string sysCmd = HERest;
				sysCmd.append(" -C ");
				sysCmd.append(config);
				sysCmd.append(" -S ");
				sysCmd.append(szTmp);
				sysCmd.append(" -I ");
				sysCmd.append(szTmpMLF);
				sysCmd.append(" -H ");
				sysCmd.append(GMMCopy);
				sysCmd.append(" -u mpw ");
				sysCmd.append(words);					
				system(sysCmd.c_str());
				unordered_map<string, GMM> hGMMsTmp;
				vector<GMM> GMMsTmp;
				GMMAdaptTool::LoadGMM(GMMCopy, hGMMsTmp);
				GMMAdaptTool::GMMHash2VectBH(hGMMsTmp, GMMsTmp);
				double AccHand = 0.0;
				double AccBkGnd = 0.0;
				double AccHand0 = 0.0;
				double AccBkGnd0 = 0.0;
				if(vvDataHand.size())
				{
					vector<int> vHandClasses(vvDataHand.size(), 1);
					vector<int> vRecogHandClasses, vRecogHandClasses0;
					vector<int> confres, confres0;
					Test2ClassesGMMs(vvDataHand, vHandClasses, GMMsTmp, vRecogHandClasses, AccHand, confres);
					Test2ClassesGMMs(vvDataHand, vHandClasses, GMMs, vRecogHandClasses0, AccHand0, confres0);
				}
				if(vvDataBkGnd.size())
				{
					vector<int> vBkGndClasses(vvDataBkGnd.size(), 0);
					vector<int> vRecogBkGndClasses, vRecogBkGndClasses0;
					vector<int> confres, confres0;
					Test2ClassesGMMs(vvDataBkGnd, vBkGndClasses, GMMsTmp, vRecogBkGndClasses, AccBkGnd, confres);
					Test2ClassesGMMs(vvDataBkGnd, vBkGndClasses, GMMs, vRecogBkGndClasses0, AccBkGnd0, confres0);
				}
				double N1 = (double) vvDataHand.size();
				double N2 = (double) vvDataBkGnd.size();
				double Accn = (AccHand * N1 + AccBkGnd * N2)/(N1+N2);
				double Accn0 = (AccHand0 * N1 + AccBkGnd0 * N2)/(N1+N2);
				Acc = (Accn * (N1+N2) + Acc * N_1)/(N_1+N1+N2);
				Acc0 = (Accn0 * (N1+N2) + Acc0 * N_1)/(N_1+N1+N2);
				N_1 += (N1+N2);
			}
			fin.close();
			cout<<"Original accuracy is "<<Acc0<<endl;
			cout<<"After framewise adjustment, accuracy is "<<Acc<<endl;
		}
		return 0;
	}

	static int GMMHTKMAPBySCP2(string GMMPaths, string DataSCPPaths, double& Acc)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
		GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
		char caBuf1[4096];
		char caBuf2[4096];
		memset(caBuf1, 0, 4096);

		memset(caBuf2, 0, 4096);
		ifstream fin;
		fin.open(DataSCPPaths.c_str());
		string szTmp = DataSCPPaths + "_Tmp";
		string szTmpMLF = DataSCPPaths + "_Tmp.mlf";
		//string HERest = "/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
		//string config = "/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
		//string words = "/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
		//string GMMCopy = GMMPaths + "_Tmp";
		//string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
		string GMMCopy = GMMPaths + "_Tmp";
	#ifdef _MSC_VER
		string HERest = "C:\\debug_hand_detection\\tools\\HTK\\bin.win32\\HERest.exe"; //"C:\\CollectedDataCopy\\recognizer\\HTK\\bin.win32\\HERest.exe";
		string config = "C:\\debug_hand_detection\\config\\config"; //"C:\\Transfer\\AdjustGMMNew\\config";
		string words = "C:\\debug_hand_detection\\config\\words"; //"C:\\Transfer\\AdjustGMMNew\\words";
		string szCopyGMM = "copy " + GMMPaths + " " + GMMCopy;
	#elif defined(__GNUC__)
		string HERest = "HERest"; //"/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
		string config = "/nv/pec1/ycheng31/local/include/AdjustGMMNew/config";//"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
		string words = "/nv/pec1/ycheng31/local/include/AdjustGMMNew/words";//"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
		string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
	#endif
		Acc = 0.0;
		double N_1 = 0;
		unordered_map< string, vector<string> > hvFiles;
		vector<string> vszKeys;
		double Acc0 = 0;
		if(fin.is_open())
		{
			while(fin.getline(caBuf1, 4096))
			{
							string sz1(caBuf1);
				if(caBuf1[sz1.length() - 1] == '\r' || caBuf1[sz1.length() - 1] == '\n')
				{
					caBuf1[sz1.length() - 1] = '\0';
					sz1 = caBuf1;
				}
							int nPos01 = (int)sz1.find_last_of("_");
							int nPos00 = (int)sz1.find_last_of("/");
							string szSub1 = sz1.substr(nPos00 + 1, nPos01 - nPos00 - 1);
							int nPos02 = (int) szSub1.find_last_of("_");
							string szKey = szSub1.substr(0, nPos02);
				if(hvFiles.find(szKey) == hvFiles.end())
				{
					vszKeys.push_back(szKey);
				}
				hvFiles[szKey].push_back(sz1);
			}
			fin.close();
			for(size_t n = 0; n < vszKeys.size(); n++)
			{
				size_t nFiles = hvFiles[vszKeys[n]].size();
				ofstream fout001, fout002;
				fout001.open(szTmp.c_str());
				fout002.open(szTmpMLF.c_str());
				double Accn = 0.0;
				if(!fout001.is_open())
				{
					cerr<<"Cannot open scp file for writing"<<endl;
					return -1;
				}
				if(!fout002.is_open())
				{
					cerr<<"Cannot open mlf file for writing"<<endl;
					return -1;
				}
				fout002<<"#!MLF!#"<<endl;
				vector< vector<double> > vvData;
				vector<int> vClasses, vRecogClasses;
				for(size_t m = 0; m < nFiles; m++)
				{
					vector< vector<double> > vvDatam;
					string szName = hvFiles[vszKeys[n]][m];
					GMMAdaptTool::ReadHTKMFCCFeature(vvDatam, szName);
					fout001<<szName<<endl;
					int nPos0001 = (int)szName.find_last_of("/");
					int nPos0002 = (int)szName.find_last_of(".");
					string szNameSub = szName.substr(nPos0001 + 1, nPos0002 - nPos0001 - 1);
					int c = 0;
					fout002<<"\"*/"<<szNameSub<<".lab\""<<endl;
					if(0 <= (int)szName.find("Hand"))
					{
						c = 1;
						fout002<<"Hand\n.\n";
					}
					else
					{
						fout002<<"BkGnd\n.\n";
					}
					for(size_t i = 0; i < vvDatam.size(); i++)
					{
						vvData.push_back(vvDatam[i]);
						vClasses.push_back(c);	
					}
					
				}
				fout001.close();
				fout002.close();
				system(szCopyGMM.c_str());
							string sysCmd = HERest;
							sysCmd.append(" -C ");
							sysCmd.append(config);
							sysCmd.append(" -S ");
							sysCmd.append(szTmp);
							sysCmd.append(" -I ");
							sysCmd.append(szTmpMLF);
							sysCmd.append(" -H ");
							sysCmd.append(GMMCopy);
							sysCmd.append(" -u mpw ");
							sysCmd.append(words);
							system(sysCmd.c_str());
				unordered_map<string, GMM> hGMMsTmp;
				vector<GMM> GMMsTmp;
				GMMAdaptTool::LoadGMM(GMMCopy, hGMMsTmp);
				GMMAdaptTool::GMMHash2VectBH(hGMMsTmp, GMMsTmp);
				vector<int> vRecogClasses0;
				double Accn0 = 0.0;
				vector<int> confres0, confres;
				Test2ClassesGMMs(vvData, vClasses, GMMs, vRecogClasses0, Accn0, confres0);
				Test2ClassesGMMs(vvData, vClasses, GMMsTmp, vRecogClasses, Accn, confres);							
				double Nn = (double)vClasses.size();
				Acc0 = (Accn0 * Nn + Acc0 * N_1)/(N_1 + Nn);
				Acc = (Accn * Nn + Acc * N_1)/(N_1 + Nn);
				N_1 += Nn;
			}	
			fin.close();
			cout<<"Original accuracy is "<<Acc0<<endl;
			cout<<"After adjust through all frames, accuracy is "<<Acc<<endl;
		}
		return 0;
	}


	static int GMMHTKMAPByLABTestNonSkinAdaptOnly(string GMMPaths, string DataLABPaths, string szRefImgList, double& Acc, ofstream& fout)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
		GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
		int nPos0001 = DataLABPaths.find_last_of("/");
		int nPob0001 = DataLABPaths.find_last_of("\\");
		int nPos01 = (nPob0001>nPos0001)?nPob0001:nPos0001;
		string szCurrFolder = DataLABPaths.substr(0, nPos01);
		string szDelim = "";
		if(nPos01 == nPob0001) szDelim = "\\";
		else if(nPos01 == nPos0001) szDelim = "/";
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		//string sysMkDirFE = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string sysMkDirFE = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string sysMkDirFE = "mkdir -p " + szFEFolder;
		#endif

		struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(sysMkDirFE.c_str());
		string szTmp = DataLABPaths + "_Tmp";
		string szTmpMLF = DataLABPaths + "_Tmp.mlf";
		string GMMCopy = GMMPaths + "_Tmp";
		string GMMCopyCopy  = GMMPaths + "_TmpTmp";
	#ifdef _MSC_VER
		string HERest = "C:\\debug_hand_detection\\tools\\HTK\\bin.win32\\HERest.exe"; //"C:\\CollectedDataCopy\\recognizer\\HTK\\bin.win32\\HERest.exe";
		string config = "C:\\debug_hand_detection\\config\\config"; //"C:\\Transfer\\AdjustGMMNew\\config";
		string words = "C:\\debug_hand_detection\\config\\words"; //"C:\\Transfer\\AdjustGMMNew\\words";
		string szCopyGMM = "copy " + GMMPaths + " " + GMMCopy;
	#elif defined(__GNUC__)
		string HERest = "HERest"; //"/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
		string config = "/nv/pec1/ycheng31/local/include/AdjustGMMNew/config";//"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
		string words = "/nv/pec1/ycheng31/local/include/AdjustGMMNew/words"; //"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
		string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
	#endif
		string sysCmd = HERest;
		sysCmd.append(" -C ");
		sysCmd.append(config);
		sysCmd.append(" -S ");
		sysCmd.append(szTmp);
		sysCmd.append(" -I ");
		sysCmd.append(szTmpMLF);
		sysCmd.append(" -H ");
		string sysCmd2 = sysCmd;
		sysCmd.append(GMMCopy);
		sysCmd.append(" -u mpw ");
		sysCmd.append(words);					
		sysCmd2.append(GMMCopy);
		sysCmd2.append(" -u mpw ");
		sysCmd2.append(words);					
		Acc = 0.0;
		double Acc0 = 0.0;
		double N_1 = 0;
		int BLKSIZE = 16;
		// Now, load label
		int nCountX = 0;
		int nCountO = 0;
		unordered_map<string, int> hszJPGFiles;
		unordered_map<string, int> hszJPGFiles1;
		unordered_map<string, int> hszJPGFiles0;
		//unordered_map<string, int> hszCurrJPGFiles;
		//unordered_map<string, JPGLabelInfo> hJPGLavelInfo;
		//unordered_map< string, vector<int> > hJPGClasses;
		//vector<string> vszKeyOrder;
		string szLabelPathName = DataLABPaths;
		//unordered_map<string, unordered_map<string, int> > hhTrainFiles;

		//int nCountO1 = 0;
		//int nCountX1 = 0;
		//int nCountO10 = 0;
		//int nCountX10 = 0;
		//unordered_map<string, JPGLabelInfo> hJPGLavelInfo1;
		//unordered_map<string, JPGLabelInfo> hJPGLavelInfo0;
		LoadJPGLabelInfo lb, lb0, lb1;
		string szLabelPathName1 = DataLABPaths + "_BKGNDMAPPFHTK";
		string szLabelPathName0 = DataLABPaths + "_Orig01";
		//==============================================================================//
		// inserted, 6/7/2011 You-Chi Cheng
		//===============================================================================//
		vector<string> vszKeyOrderTmp;
		//===============================================================================//
		//GMMAdaptTool::LoadJPGInfo(BLKSIZE, 4, nCountX, nCountO,
		//							hszJPGFiles,
		//							//hszCurrJPGFiles,
		//							hJPGLavelInfo,
		//							hJPGClasses, 
		//							vszKeyOrderTmp,
		//							szLabelPathName,
		//							hhTrainFiles);
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		//==============================================================================//
		// inserted, 6/7/2011 You-Chi Cheng
		//===============================================================================//
		FixlabelInfoByRefFileList(lb.m_hJPGLavelInfo, lb.m_hJPGClasses, lb.m_vszKeyOrder, szRefImgList, BLKSIZE);
		unordered_map<string, int> hSeqStart;
		unordered_map<string, int> hSeqStop;
		unordered_map<string, vector<string> > hvSeqs;
		vector<string> vSeqs;
		TestImageStartStop(lb.m_vszKeyOrder, hSeqStart, hSeqStop, hvSeqs, vSeqs);
		// New Version:
		vector< vector<double> > vvCurrFeaturesTmp0001;
		for(size_t n = 0; n < vSeqs.size(); n++)
		{
			// Adaptation using estimated background only:
			vector< vector<double> > vvRecog1Tmp0001;
			string szCurrFE01Tmp = szFEFolder + vSeqs[n] +"_BkGnd_Not_Skin.fe";	
			for(size_t i = 0; i < hSeqStart[vSeqs[n]]; i++)
			{
				string szCurrFilePath = hvSeqs[vSeqs[n]][i];
				string szImagePath;
				string szImageName;
				string szSeqName;
				string szExt;
				GetFileNameInfo(szCurrFilePath, szImagePath, szImageName, szSeqName, szExt);
				vector<int> vTmpRecogClasses0001;
				double Accn0Tmp = 0.0;
				JPGLabelInfo* pinfo = &lb.m_hJPGLavelInfo[szCurrFilePath];
				vector< vector<double> > vvCurrFeaturesTmp0001;
				ExtColorFeature(szCurrFilePath, pinfo, vvCurrFeaturesTmp0001, BLKSIZE);
				vector<int> vTmpClasses0001(vvCurrFeaturesTmp0001.size(), 0);
				//vector<int> vTmpRecogClasses0001;
				//double Accn0Tmp = 0.0;
				vector<int> confres0Tmp;
				Test2ClassesGMMs(vvCurrFeaturesTmp0001, vTmpClasses0001, GMMs, vTmpRecogClasses0001, Accn0Tmp, confres0Tmp);
				for(size_t n000101 = 0; n000101 < vvCurrFeaturesTmp0001.size(); n000101++)
				{
					if(vTmpRecogClasses0001[n000101] != 0)
					{
						vvRecog1Tmp0001.push_back(vvCurrFeaturesTmp0001[n000101]);
					}
				}			
			}
			ofstream foutTmp;
			foutTmp.open(szTmp.c_str());
			if(!foutTmp.is_open())
			{
				cerr<<"Cannot Open Temp File so stop it"<<endl;
				return -1;
			}
			ofstream fout02Tmp;
			fout02Tmp.open(szTmpMLF.c_str());
			if(!fout02Tmp.is_open())
			{
				cerr<<"Cannot open mlf file for written!"<<endl;
				return -1;
			}
			fout02Tmp<<"#!MLF!#"<<endl;
			if(0 < vvRecog1Tmp0001.size())
			{
				WriteFeatureAsHTKMFCC(vvRecog1Tmp0001, szCurrFE01Tmp, 9);
				foutTmp<<szCurrFE01Tmp<<endl;
				fout02Tmp<<"\"*/"<<vSeqs[n]<<"_BkGnd_Not_Skin.lab\""<<endl;
				fout02Tmp<<"BkGnd\n.\n";
			}
			foutTmp.close();
			fout02Tmp.close();
			// Perform HTK:
			if(0 < vvRecog1Tmp0001.size())
			{
				system(szCopyGMM.c_str());
				system(sysCmd.c_str());
			}
			// 
			unordered_map<string, GMM> hGMMsTmp;
			vector<GMM> GMMsTmp;
			GMMAdaptTool::LoadGMM(GMMCopy, hGMMsTmp);
			GMMAdaptTool::GMMHash2VectBH(hGMMsTmp, GMMsTmp);
			// Now that the model is updated, let's go back to skin block decision:
			for(size_t i = 0; i < hvSeqs[vSeqs[n]].size(); i++)
			{
				string szCurrFilePath = hvSeqs[vSeqs[n]][i];
				string szImagePathi, szImageNamei, szSeqNamei, szExti;
				GetFileNameInfo(szCurrFilePath, szImagePathi, szImageNamei, szSeqNamei, szExti);
						string szCurrFE1 = szFEFolder + szImageNamei +"_Hand.fe";
						string szCurrFE0 = szFEFolder + szImageNamei +"_BkGnd.fe";
				JPGLabelInfo* pinfo = &lb.m_hJPGLavelInfo[szCurrFilePath];
				vector< vector<double> > vvCurrFeatures;
				ExtColorFeature(szCurrFilePath, pinfo, vvCurrFeatures, BLKSIZE);
				int nHeight = pinfo->m_nBlockHeight;
				int nWidth = pinfo->m_nBlockWidth;
				vector<int> vCurrRecogClasses0, vCurrRecogClasses;
				double Accn0 = 0.0;
				double Accn = 0.0;
				vector<int> confres0;
				Test2ClassesGMMs(vvCurrFeatures, lb.m_hJPGClasses[szCurrFilePath], GMMsTmp, vCurrRecogClasses0, Accn0, confres0);
				vector< vector<double> > vvRecog0, vvRecog1;
				// Collect those recognized as background and hand, respectively
				for(int k1 = 0; k1 < nWidth; k1++)
				{
					for(int k2 = 0; k2 < nHeight; k2++)
					{
						int nCurrIndex = k1 * nHeight + k2;
						if(vCurrRecogClasses0[nCurrIndex] == 0)
						{
							vvRecog0.push_back(vvCurrFeatures[nCurrIndex]);
						}
						else if(vCurrRecogClasses0[nCurrIndex] == 1)
						{
							vvRecog1.push_back(vvCurrFeatures[nCurrIndex]);
							lb0.MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1, nWidth, nHeight);
							if(lb0.m_hszJPGFiles.find(szCurrFilePath) == lb0.m_hszJPGFiles.end())
							{
								lb0.m_vszKeyOrder.push_back(szCurrFilePath);
							}					
							lb0.m_hszJPGFiles[szCurrFilePath] = 1;
							//GMMAdaptTool::MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1,
												//											hJPGLavelInfo0,
												//											nCountO10,
												//											nCountX10,
												//											nWidth,
												//											nHeight,
												//											BLKSIZE);
						}
					}
				}
				
				if(/*hJPGLavelInfo0*/lb0.m_hJPGLavelInfo.find(szCurrFilePath) == lb0.m_hJPGLavelInfo/*hJPGLavelInfo0*/.end())
				{
					/*hJPGLavelInfo0*/lb0.m_hJPGLavelInfo[szCurrFilePath].SetWidthHeight(nWidth*BLKSIZE, nHeight*BLKSIZE, BLKSIZE);
					/*hJPGLavelInfo0*/lb0.m_hJPGLavelInfo[szCurrFilePath].m_szPathName = szCurrFilePath;				
				}
				// Prepare mlf and scp for HTK
				ofstream fout;
				fout.open(szTmp.c_str());
				if(!fout.is_open())
				{
					cerr<<"Cannot Open Temp File so stop it"<<endl;
					return -1;
				}
				ofstream fout02;
				fout02.open(szTmpMLF.c_str());
				if(!fout02.is_open())
				{
					cerr<<"Cannot open mlf file for written!"<<endl;
					return -1;
				}
				fout02<<"#!MLF!#"<<endl;
				if(0 < vvRecog0.size())
				{
					WriteFeatureAsHTKMFCC(vvRecog0, szCurrFE0, 9);
					fout<<szCurrFE0<<endl;
					fout02<<"\"*/"<<szImageNamei<<"_BkGnd.lab\""<<endl;
					fout02<<"BkGnd\n.\n";
				}
				if(0 < vvRecog1.size())
				{
					WriteFeatureAsHTKMFCC(vvRecog1, szCurrFE1, 9);
					fout<<szCurrFE1<<endl;
					fout02<<"\"*/"<<szImageNamei<<"_Hand.lab\""<<endl;
					fout02<<"Hand\n.\n";
				}
				fout.close();
				fout02.close();
				// Run HTK 
				system(sysCmd.c_str());
				//cout<<szImageNamei<<endl;
				unordered_map<string, GMM> hGMMsTmpTmp;
				vector<GMM> GMMsTmpTmp;
				GMMAdaptTool::LoadGMM(GMMCopy, hGMMsTmpTmp);
				GMMAdaptTool::GMMHash2VectBH(hGMMsTmpTmp, GMMsTmpTmp);
				//vector<int> vCurrRecogClasses;
				vector<int> confres;
				Test2ClassesGMMs(vvCurrFeatures, lb.m_hJPGClasses[szCurrFilePath], GMMsTmpTmp, vCurrRecogClasses, Accn, confres);
				hszJPGFiles1[szCurrFilePath] = 1;
				hszJPGFiles0[szCurrFilePath] = 1;
				for(int k1 = 0; k1 < nWidth; k1++)
				{
					for(int k2 = 0; k2 < nHeight; k2++)
					{
						int nCurrIndex = k1 * nHeight + k2;
						if(vCurrRecogClasses[nCurrIndex] == 1)
						{
							lb1.MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1, nWidth, nHeight);
							if(lb1.m_hszJPGFiles.find(szCurrFilePath) == lb0.m_hszJPGFiles.end())
							{
								lb1.m_vszKeyOrder.push_back(szCurrFilePath);
							}					
							lb1.m_hszJPGFiles[szCurrFilePath] = 1;
							//GMMAdaptTool::MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1, 
							//												hJPGLavelInfo1, 	
							//												nCountO1,
							//												nCountX1,
							//												nWidth,
							//												nHeight,
							//												BLKSIZE);
						}
					}
				}
				if(/*hJPGLavelInfo1*/lb0.m_hJPGLavelInfo.find(szCurrFilePath) == /*hJPGLavelInfo1*/lb0.m_hJPGLavelInfo.end())
				{
					/*hJPGLavelInfo1*/lb0.m_hJPGLavelInfo[szCurrFilePath].SetWidthHeight(nWidth*BLKSIZE, nHeight*BLKSIZE, BLKSIZE);
					/*hJPGLavelInfo1*/lb0.m_hJPGLavelInfo[szCurrFilePath].m_szPathName = szCurrFilePath;				
				}
				double NN = (double)(nWidth * nHeight);
				Acc = (Accn * NN + Acc * N_1)/(N_1+NN);
				Acc0 = (Accn0 * NN + Acc0 * N_1)/(N_1+NN);
				N_1 += NN;	
			}
		}
		if(fout.is_open())
		{
			fout<<"Original result is "<<Acc0<<endl;
			fout<<"HTK MAP for each frame, result is "<<Acc<<endl;
		}
		else
		{
			cout<<"Original result is "<<Acc0<<endl;
			cout<<"HTK MAP for each frame, result is "<<Acc<<endl;
		}
		lb1.SaveJPGInfo(szLabelPathName1);
		lb0.SaveJPGInfo(szLabelPathName0);

		//GMMAdaptTool::SaveJPGInfo(szLabelPathName1, 
		//			hJPGLavelInfo1,
		//			vszKeyOrder, 
		//			hszJPGFiles1);
	 //   GMMAdaptTool::SaveJPGInfo(szLabelPathName0,
	 //                               hJPGLavelInfo0,
	 //                               vszKeyOrder,
	 //                               hszJPGFiles0);

		return 0;
	}

	static void FindCurrentFolderANDDelim(string& szFileName, string& szFolder, string& szDelim)
	{
		int nPos0001 = szFileName.find_last_of("/");
		int nPob0001 = szFileName.find_last_of("\\");
		int nPos01 = (nPob0001>nPos0001)?nPob0001:nPos0001;
		szDelim = (nPob0001>nPos0001)?"\\":"/";
		szFolder = (nPos01 <= 0)?"":szFileName.substr(0, nPos01);
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


	static void FindCurrentFileName(string& szPathName, string& szName)
	{
		int nPos001s = szPathName.find_last_of("/");
		int nPos001b = szPathName.find_last_of("\\");
		int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
		int nPos002 = szPathName.find_last_of(".");
		nPos002 = (nPos002 > 0)?nPos002:(szPathName.length() - 1);
		szName = szPathName.substr(nPos001 + 1, nPos002 - nPos001 - 1);
	}

	static void FindCurrentSeqNamePath(string& szPathName, string& szName, string& szPath)
	{
		int nPos001s = szPathName.find_last_of("/");
		int nPos001b = szPathName.find_last_of("\\");
		int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
		int nPos002 = szPathName.find_last_of("_");
		nPos002 = (nPos002 > 0)?nPos002:(szPathName.length() - 1);
		szName = szPathName.substr(nPos001 + 1, nPos002 - nPos001);
		szPath = szPathName.substr(0, nPos001 + 1);
	}

	static void MakeFullFilePath(string& szDir, string& szName, string& szExt, string& szFullPathName)
	{
		string szdelim;
		#ifdef _WIN32
			szdelim = "\\";
		#else
			szdelim = "/";
		#endif
		int nPos001s = szDir.find_last_of("/");
		int nPos001b = szDir.find_last_of("\\");
		int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
		if(nPos001 < szDir.length() - 1){
			szDir = szDir + szdelim;
		}
		szFullPathName = szDir + szName + szExt;
	}

	static long LoadImagesFromPaths(vector<string>& vszPaths, vector<cv::Mat>& vImgs, vector<string>& vszLoadedFiles)
	{
		long N = (long) vszPaths.size();
		if(N <= 0)
		{
			cerr<<"Cannot load images, no path specified"<<endl;
			return ERROR_INVALID_DATA;
		}
		for(long n = 0; n < N; n++)
		{
			//cv::Mat imgtmp = cv::imread(vszPaths[n]);
			// =====================================================================//
			// For data normalize
			// =====================================================================//
			//cv::Mat img;// = img.clone();
			//imgtmp.copyTo(img);
			//img.release();
			//IplImage iplimgtmp = imgtmp.operator IplImage();
			//IplImage iplimg = img.operator IplImage();
			//IplImage* piplimg = &iplimg;
			//IplImage* piplimgtmp = &iplimgtmp;
			//NormalizeImageRGB(piplimgtmp, piplimg);
			//img = cv::Mat(piplimg);
			// =====================================================================//
			cv::Mat img = cv::imread(vszPaths[n]);
			if(img.rows > 0 && img.cols > 0)
			{
				vImgs.push_back(img);
				vszLoadedFiles.push_back(vszPaths[n]);
			}
		}
		return S_OK;
	}

	static long ExtColorFeaturesFromImages(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, int BLOCKSIZE)
	{
		long N = (long) vImgs.size();
		if(N <= 0)
		{
			cerr<<"No images for feature extraction"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvvFeatures.resize(N);
		for(long n = 0; n < N; n++)
		{
			IplImage img = vImgs[n].operator IplImage();
			IplImage* pimg = &img;
			ExtLowLevelImgFeature::ExtColorFeatures(pimg, vvvFeatures[n], BLOCKSIZE);
		}
		return S_OK;
	}

	//static long ExtColorFeaturesFromImages_BS0(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, int BLOCKSIZE, cv::Mat& ref, vector< vector< vector<double> > >& vvvCorrs)
	static long ExtColorFeaturesFromImages_BS0(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, int BLOCKSIZE, cv::Mat& ref, vector< vector< vector<double> > >& vvvCorrs, vector< vector<double> >& vvKSs)
	{
		long N = (long) vImgs.size();
		if(N <= 0)
		{
			cerr<<"No images for feature extraction"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvvFeatures.resize(N);
		vvvCorrs.resize(N);
		vvKSs.resize(N);
		IplImage iplref = ref.operator IplImage();
		IplImage* pref = &iplref;
		for(long n = 0; n < N; n++)
		{
			IplImage img = vImgs[n].operator IplImage();
			IplImage* pimg = &img;
			//int yrange = ExtLowLevelImgFeature::ExtColorFeatures_BS0(pimg, vvvFeatures[n], BLOCKSIZE, pref, vvvCorrs[n]);
			int yrange = ExtLowLevelImgFeature::ExtColorFeatures_BS0(pimg, vvvFeatures[n], BLOCKSIZE, pref, vvvCorrs[n], vvKSs[n]);
		}
		return S_OK;
	}


	//static long ExtColorFeaturesFromImagesStartStop(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, int BLOCKSIZE, int nStart, int nStop, IplImage*& pRefImg, vector< vector< vector<double> > >& vvvCorrs)
	static long ExtColorFeaturesFromImagesStartStop(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, int BLOCKSIZE, int nStart, int nStop, IplImage*& pRefImg, vector< vector< vector<double> > >& vvvCorrs, vector< vector<double> >& vvKSs)
	{
		long N = (long) vImgs.size();
		if(N <= 0)
		{
			cerr<<"No images for feature extraction"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvvFeatures.resize(N);
		vvvCorrs.resize(N);
		vvKSs.resize(N);
		for(long n = nStart; n < nStop; n++)
		{
			IplImage img = vImgs[n].operator IplImage();
			IplImage* pimg = &img;
			//ExtLowLevelImgFeature::ExtColorFeatures_BS0(pimg, vvvFeatures[n], BLOCKSIZE, pRefImg, vvvCorrs[n]);
			ExtLowLevelImgFeature::ExtColorFeatures_BS0(pimg, vvvFeatures[n], BLOCKSIZE, pRefImg, vvvCorrs[n], vvKSs[n]);
			//ExtLowLevelImgFeature::ExtColorFeatures(pimg, vvvFeatures[n], BLOCKSIZE);
		}
		return S_OK;
	}

	static long ExtSLICColorFeaturesFromImagesStartStop(vector<string>& vszLoadedFiles, vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, int nStart, int nStop, IplImage*& pRefImg, vector< vector< vector<double> > >& vvvCorrs, vector< vector<double> >& vvKSs)
	{
		long N = (long) vImgs.size();
		if(N <= 0)
		{
			cerr<<"No images for feature extraction"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvvFeatures.resize(N);
		vvvCorrs.resize(N);
		vvKSs.resize(N);
		for(long n = nStart; n < nStop; n++)
		{
			IplImage img = vImgs[n].operator IplImage();
			IplImage* pimg = &img;
			string szSLICPathn;
			GuessSLICPath(vszLoadedFiles[n], szSLICPathn);
			cv::Size imsize = vImgs[n].size();
			vector< vector<cv::Point> > vvSLICLabels;
			bool SLIC_is_given = My_SLIC_Interface::LoadSLICDat(szSLICPathn, imsize.width, imsize.height, vvSLICLabels);
			ExtLowLevelImgFeature::ExtColorFeatures_SLIC_BS0(pimg, vvvFeatures[n], vvSLICLabels, 200, pRefImg, vvvCorrs[n], vvKSs[n], imsize, SLIC_is_given);
		}
		return S_OK;
	}


	static long TestAFrameFeatures(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector<int>& vRecogClasses, 
								vector< vector< vector<double> > >& vvvRecogVec , long index)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}

		long I = (long)GMMs.size();
		if(vvvRecogVec.size() != I)
		{
			vvvRecogVec.resize(I);
		}
		long n = index;
		RecognizeByGMMs(vvvFeatures[n], GMMs, vRecogClasses);
		long M = (long)vRecogClasses.size();
		for(long m = 0; m < M; m++)
		{
			int i = vRecogClasses[m];
			vvvRecogVec[i].push_back(vvvFeatures[n][m]);
		}
		return 0;
	}

	static long TestAFrameFeatures_Cheat(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector<int>& vRecogClasses_Cheat, 
								vector< vector< vector<double> > >& vvvRecogVec , long index)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}

		long I = (long)GMMs.size();
		if(vvvRecogVec.size() != I)
		{
			vvvRecogVec.resize(I);
		}
		long n = index;
		long M = (long)vRecogClasses_Cheat.size();
		for(long m = 0; m < M; m++)
		{
			int i = vRecogClasses_Cheat[m];
			vvvRecogVec[i].push_back(vvvFeatures[n][m]);
		}
		return 0;
	}

	static long TestAFrameFeatures_BS0_Cheat(vector< vector< vector<double> > >& vvvFeatures,
					vector<GMM>& GMMs,
					vector<int>& vRecogClassesCheat,
					vector< vector< vector<double> > >& vvvRecogVec , long index, 
					vector< vector< vector<double> > >& vvvCorrs)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
				cerr<<"No features for recognition"<<endl;
				return ERROR_INVALID_DATA;
		}

		long I = (long)GMMs.size();
		if(vvvRecogVec.size() != I)
		{
				vvvRecogVec.resize(I);
		}
		long n = index;
		//RecognizeByGMMs(vvvFeatures[n], GMMs, vRecogClasses);
		//vector<int> vRecogClasses(vRecogClassesCheat.size(), 0);
		for(size_t kk = 0; kk < vRecogClassesCheat.size(); kk++)
		{
			double dbCorr = 0.1140 * vvvCorrs[n][kk][0] + 0.5870 * vvvCorrs[n][kk][1] + 0.2989 * vvvCorrs[n][kk][2];
			//vvRecogClasses[n][kk] = (dbCorr < 0.5)?vvRecogClasses[n][kk]:0;
			if(vRecogClassesCheat[kk] == 1)
			{
				if(dbCorr < CCThreshold)
				{
					vvvRecogVec[1].push_back(vvvFeatures[n][kk]);
				}
			}
			else if(vRecogClassesCheat[kk] == 0)
			{
				if(dbCorr >= CCThreshold)
				{
					vvvRecogVec[0].push_back(vvvFeatures[n][kk]);
				}
			}
		}
	}


	static long TestAFrameFeatures_BS0_Scores(vector< vector< vector<double> > >& vvvFeatures,
					vector<GMM>& GMMs,
					vector< vector<double> >& vvRecogScores,
					long index, 
					vector< vector< vector<double> > >& vvvCorrs)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
				cerr<<"No features for recognition"<<endl;
				return ERROR_INVALID_DATA;
		}

		long I = (long)GMMs.size();
		long n = index;
		GMMsScores(vvvFeatures[n], GMMs, vvRecogScores);
		for(size_t kk = 0; kk < vvRecogScores.size(); kk++)
		{
			vvRecogScores[kk].push_back(vvvCorrs[n][kk][0]);
			vvRecogScores[kk].push_back(vvvCorrs[n][kk][1]);
			vvRecogScores[kk].push_back(vvvCorrs[n][kk][2]);
		}
	}

	static long TestAFrameFeatures_BS0(vector< vector< vector<double> > >& vvvFeatures,
					vector<GMM>& GMMs,
					vector<int>& vRecogClasses,
					vector< vector< vector<double> > >& vvvRecogVec , long index, 
					vector< vector< vector<double> > >& vvvCorrs)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
				cerr<<"No features for recognition"<<endl;
				return ERROR_INVALID_DATA;
		}

		long I = (long)GMMs.size();
		if(vvvRecogVec.size() != I)
		{
				vvvRecogVec.resize(I);
		}
		long n = index;
		RecognizeByGMMs(vvvFeatures[n], GMMs, vRecogClasses);
		for(size_t kk = 0; kk < vRecogClasses.size(); kk++)
		{
			double dbCorr = 0.1140 * vvvCorrs[n][kk][0] + 0.5870 * vvvCorrs[n][kk][1] + 0.2989 * vvvCorrs[n][kk][2];
			//vvRecogClasses[n][kk] = (dbCorr < 0.5)?vvRecogClasses[n][kk]:0;
			if(vRecogClasses[kk] == 1)
			{
				if(dbCorr >= CCThreshold)
				{
					vRecogClasses[kk] = 0;
				}
			}
		}
		long M = (long)vRecogClasses.size();
		for(long m = 0; m < M; m++)
		{
			int i = vRecogClasses[m];
			vvvRecogVec[i].push_back(vvvFeatures[n][m]);
		}
	}

	static long TestAFrameFeatures_BS0_V2(vector< vector< vector<double> > >& vvvFeatures,
					vector<GMM>& GMMs,
					vector<int>& vRecogClasses,
					vector< vector< vector<double> > >& vvvRecogVec , long index, 
					vector< vector< vector<double> > >& vvvCorrs)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
				cerr<<"No features for recognition"<<endl;
				return ERROR_INVALID_DATA;
		}

		long I = (long)GMMs.size();
		if(vvvRecogVec.size() != I)
		{
				vvvRecogVec.resize(I);
		}
		long n = index;
		RecognizeByGMMs(vvvFeatures[n], GMMs, vRecogClasses);
		for(size_t kk = 0; kk < vRecogClasses.size(); kk++)
		{
			double dbCorr = 0.1140 * vvvCorrs[n][kk][0] + 0.5870 * vvvCorrs[n][kk][1] + 0.2989 * vvvCorrs[n][kk][2];
			//vvRecogClasses[n][kk] = (dbCorr < 0.5)?vvRecogClasses[n][kk]:0;
			if(vRecogClasses[kk] == 1)
			{
				if(dbCorr < CCThreshold)
				{
					vvvRecogVec[1].push_back(vvvFeatures[n][kk]);
				}
				else
				{
					vRecogClasses[kk] = 0;
				}
			}
			else if(vRecogClasses[kk] == 0)
			{
				if(dbCorr >= CCThreshold)
				{
					vvvRecogVec[0].push_back(vvvFeatures[n][kk]);
				}
			}
		}
	}

	static long TestAllFeatures(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector< vector<int> >& vvRecogClasses)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvRecogClasses.resize(N);

		for(long n = 0; n < N; n++)
		{
			RecognizeByGMMs(vvvFeatures[n], GMMs, vvRecogClasses[n]);
		}
	}

	static long TestAllFeatures_Scores(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector< vector<int> >& vvRecogClasses,
								vector< vector< vector<double> > >& vvvRecogScores
								)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvRecogClasses.resize(N);
		vvvRecogScores.resize(N);
		for(long n = 0; n < N; n++)
		{
			GMMsScoresRecog(vvvFeatures[n], GMMs, vvRecogClasses[n], vvvRecogScores[n]);
			//RecognizeByGMMs(vvvFeatures[n], GMMs, vvRecogClasses[n]);
		}
	}

	static long TestAllFeatures_BS0_Corrs_Only(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector< vector<int> >& vvRecogClasses,
								vector< vector< vector<double> > >& vvvCorrs)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvRecogClasses.resize(N);

		for(long n = 0; n < N; n++)
		{
			if(vvvFeatures[n].size())
			{
				//RecognizeByGMMs(vvvFeatures[n], GMMs, vvRecogClasses[n]);
				//cout<<"n is "<<n<<endl;
				vvRecogClasses[n] = vector<int>(vvvFeatures[n].size(), 0);
				for(size_t kk = 0; kk < vvRecogClasses[n].size(); kk++)
				{
					//cerr<<"n is "<<n<<", kk is "<<kk<<endl;
					double dbCorr = 0.1140 * vvvCorrs[n][kk][0] + 0.5870 * vvvCorrs[n][kk][1] + 0.2989 * vvvCorrs[n][kk][2];
					//vvRecogClasses[n][kk] = (dbCorr < 0.5)?vvRecogClasses[n][kk]:0;
					if(dbCorr < CCThreshold)
					{
						vvRecogClasses[n][kk] = 1;
					}
				}
			}
		}
	}

	static long TestAllFeatures_BS0(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector< vector<int> >& vvRecogClasses,
								vector< vector< vector<double> > >& vvvCorrs,
								vector< vector<double> >& vvKSs = vector< vector<double> >()
								)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvRecogClasses.resize(N);

		for(long n = 0; n < N; n++)
		{
			if(vvvFeatures[n].size())
			{
				RecognizeByGMMs(vvvFeatures[n], GMMs, vvRecogClasses[n]);
				//cout<<"n is "<<n<<endl;
				for(size_t kk = 0; kk < vvRecogClasses[n].size(); kk++)
				{
					//cerr<<"n is "<<n<<", kk is "<<kk<<endl;
					double dbCorr = 0.1140 * vvvCorrs[n][kk][0] + 0.5870 * vvvCorrs[n][kk][1] + 0.2989 * vvvCorrs[n][kk][2];
					//vvRecogClasses[n][kk] = (dbCorr < 0.5)?vvRecogClasses[n][kk]:0;
					if(vvRecogClasses[n][kk] == 1)
					{
						if(dbCorr >= CCThreshold)
						{
							vvRecogClasses[n][kk] = 0;
						}
					}
				}
			}
		}
	}

	static long AllFeatures_BS0_Scores(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector< vector< vector<double> > >& vvvRecogScores,
								vector< vector< vector<double> > >& vvvCorrs,
								vector< vector<double> >& vvKSs = vector< vector<double> >())
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvvRecogScores.resize(N);

		for(long n = 0; n < N; n++)
		{
			if(vvvFeatures[n].size())
			{
				GMMsScores(vvvFeatures[n], GMMs, vvvRecogScores[n]);
				//RecognizeByGMMs(vvvFeatures[n], GMMs, vvRecogClasses[n]);
				//cout<<"n is "<<n<<endl;
				for(size_t kk = 0; kk < vvvRecogScores[n].size(); kk++)
				{
					vvvRecogScores[n][kk].push_back(vvvCorrs[n][kk][0]);
					vvvRecogScores[n][kk].push_back(vvvCorrs[n][kk][1]);
					vvvRecogScores[n][kk].push_back(vvvCorrs[n][kk][2]);
					if(vvKSs.size() == vvvFeatures[n].size())
					{
						vvvRecogScores[n][kk].push_back(vvKSs[n][kk]);
					}
				}
			}
		}
	}

	static long TestAllFeatures(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector< vector<int> >& vvRecogClasses, 
								vector< vector< vector<double> > >& vvvRecogVec )
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvRecogClasses.resize(N);

		long I = (long)GMMs.size();
		if(vvvRecogVec.size() != I)
		{
			vvvRecogVec.resize(I);
		}
		for(long n = 0; n < N; n++)
		{
			RecognizeByGMMs(vvvFeatures[n], GMMs, vvRecogClasses[n]);
			long M = (long)vvRecogClasses[n].size();
			for(long m = 0; m < M; m++)
			{
				int i = vvRecogClasses[n][m];
				vvvRecogVec[i].push_back(vvvFeatures[n][m]);
			}
		}
	}
	
	static long TestAllFeatures_BS0(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector< vector<int> >& vvRecogClasses,
								vector< vector< vector<double> > >& vvvRecogVec, 
								vector< vector< vector<double> > >& vvvCorrs)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvRecogClasses.resize(N);
                long I = (long)GMMs.size();
                if(vvvRecogVec.size() != I)
                {
                                vvvRecogVec.resize(I);
                }
		for(long n = 0; n < N; n++)
		{
			//cerr<<"n is "<<n<<endl;
			if(vvvFeatures[n].size())
			{
				RecognizeByGMMs(vvvFeatures[n], GMMs, vvRecogClasses[n]);
				for(size_t kk = 0; kk < vvRecogClasses[n].size(); kk++)
				{
					//cerr<<"n is "<<n<<", kk is "<<kk<<endl;
					double dbCorr = 0.1140 * vvvCorrs[n][kk][0] + 0.5870 * vvvCorrs[n][kk][1] + 0.2989 * vvvCorrs[n][kk][2];
					//vvRecogClasses[n][kk] = (dbCorr < 0.5)?vvRecogClasses[n][kk]:0;
					if(vvRecogClasses[n][kk] == 1)
					{
						if(dbCorr >= CCThreshold)
						{
							vvRecogClasses[n][kk] = 0;
						}
					}
				}
				long M = (long)vvRecogClasses[n].size();
				for(long m = 0; m < M; m++)
				{
					int i = vvRecogClasses[n][m];
					vvvRecogVec[i].push_back(vvvFeatures[n][m]);
				}
			}
		}
	}



	static long TestAllFeatures_BS0_V2(vector< vector< vector<double> > >& vvvFeatures, 
								vector<GMM>& GMMs, 
								vector< vector<int> >& vvRecogClasses,
								vector< vector< vector<double> > >& vvvRecogVec, 
								vector< vector< vector<double> > >& vvvCorrs)
	{
		long N = (long) vvvFeatures.size();
		if(N <= 0)
		{
			cerr<<"No features for recognition"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvRecogClasses.resize(N);
                long I = (long)GMMs.size();
                if(vvvRecogVec.size() != I)
                {
                                vvvRecogVec.resize(I);
                }
		for(long n = 0; n < N; n++)
		{
			//cerr<<"n is "<<n<<endl;
			if(vvvFeatures[n].size())
			{
				RecognizeByGMMs(vvvFeatures[n], GMMs, vvRecogClasses[n]);
				for(size_t kk = 0; kk < vvRecogClasses[n].size(); kk++)
				{
					//cerr<<"n is "<<n<<", kk is "<<kk<<endl;
					double dbCorr = 0.1140 * vvvCorrs[n][kk][0] + 0.5870 * vvvCorrs[n][kk][1] + 0.2989 * vvvCorrs[n][kk][2];
					//vvRecogClasses[n][kk] = (dbCorr < 0.5)?vvRecogClasses[n][kk]:0;


					if(vvRecogClasses[n][kk] == 1)
					{
						if(dbCorr < CCThreshold)
						{
							vvvRecogVec[1].push_back(vvvFeatures[n][kk]);
						}
						else
						{
							vvRecogClasses[n][kk] = 0;
						}
					}
					else if(vvRecogClasses[n][kk] == 0)
					{
						if(dbCorr >= CCThreshold)
						{
							vvvRecogVec[0].push_back(vvvFeatures[n][kk]);
						}
					}

				}
			}
		}
	}




	static long PrepareHTKMAPCommands(string GMMPaths, string DataPaths, string& GMMCopy, string& szTmpSCP, string& szTmpMLF, string& szCmd, int nIter = 1)
	{
		szTmpSCP = DataPaths + "_Tmp.scp";
		szTmpMLF = DataPaths + "_Tmp.mlf";
		string GMMCopyTmp; 
		char buf1[256];
		char buf2[256];
		memset(buf1, 0, 256);
		memset(buf2, 0, 256);
		if(nIter == 1)
		{
			char buf1[256];
			sprintf(buf2, "%i", nIter);
			//itoa(nIter, buf2,10);
			GMMCopyTmp = GMMPaths;
			GMMCopy = GMMCopyTmp + "_Tmp" + string(buf2);
		}
		else
		{
			int nIterLast = nIter - 1;
			sprintf(buf1, "%i", nIterLast); //itoa(nIterLast, buf1,10);
			sprintf(buf2, "%i", nIter); //itoa(nIter, buf2,10);
			GMMCopyTmp = GMMPaths + "_Tmp" + string(buf1);
			GMMCopy = GMMPaths + "_Tmp" + string(buf2);
		}
	#ifdef _MSC_VER
		string HERest = "C:\\debug_hand_detection\\tools\\HTK\\bin.win32\\HERest.exe"; //"C:\\CollectedDataCopy\\recognizer\\HTK\\bin.win32\\HERest.exe";
		string config = "C:\\debug_hand_detection\\config\\config"; //"C:\\Transfer\\AdjustGMMNew\\config";
		string words = "C:\\debug_hand_detection\\config\\words"; //"C:\\Transfer\\AdjustGMMNew\\words";
		string szCopyGMM = "copy " + GMMCopyTmp + " " + GMMCopy;
	#elif defined(__GNUC__)
		string HERest = "HERest"; //"/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
		string config = "/nv/pec1/ycheng31/local/include/AdjustGMMNew/config"; //"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
		string words = "/nv/pec1/ycheng31/local/include/AdjustGMMNew/words"; //"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
		string szCopyGMM = "cp " + GMMCopyTmp + " " + GMMCopy;
	#endif
		// Copy GMM as Temp
		system(szCopyGMM.c_str());
		szCmd = HERest;
		szCmd.append(" -C ");
		szCmd.append(config);
		szCmd.append(" -S ");
		szCmd.append(szTmpSCP);
		szCmd.append(" -I ");
		szCmd.append(szTmpMLF);
		szCmd.append(" -H ");
		szCmd.append(GMMCopy);
		szCmd.append(" -u mpw ");
		szCmd.append(words);
		return 0;
	}

	static long PrepareTempDocsForAdapt(vector< vector< vector<double> > >& vvvRecogVecm, vector<string>& vszClassLabs, string& szName, string& szPath, string& szTmpSCP, string& szTmpMLF)
	{
		ofstream fout1, fout2;
		fout1.open(szTmpSCP.c_str());
		fout2.open(szTmpMLF.c_str());
		if(!fout1.is_open())
		{
			cerr<<"Cannot open temp scp file for writing!\n"<<szTmpSCP<<endl;
			return S_FALSE;
		}
		if(!fout2.is_open())
		{
			cerr<<"Cannot open temp mlf file for writing!\n"<<szTmpSCP<<endl;
			return S_FALSE;
		}
		long I1 = (long) vvvRecogVecm.size();
		long I2 = (long) vszClassLabs.size();
		if(I1 != I2)
		{
			cerr<<"# Labels and # names didn't match!"<<endl;
			return ERROR_INVALID_DATA;
		}
		fout2<<"#!MLF!#"<<endl;
		bool bHasAdapt = false;
		for(long i = 0; i < I1; i++)
		{
			string szOutFE = szPath + szName + "_" + vszClassLabs[i] + ".fe";
			
			if(vvvRecogVecm[i].size())
			{
				WriteFeatureAsHTKMFCC(vvvRecogVecm[i], szOutFE, 9);
				fout1<<szOutFE<<endl;
				fout2<<"\"*/"<<szName<<"_"<<vszClassLabs[i]<<".lab\""<<endl;
				fout2<<vszClassLabs[i]<<"\n.\n";
				bHasAdapt = true;
			}

		}
		fout1.close();
		fout2.close();
		if(bHasAdapt)
			return S_OK;
		else
			return ERROR_EMPTY;
	}

	static long Do_Skin_HTKMAP(vector< vector< vector<double> > >& vvvRecogVecm, string GMMPaths, string szFEFolder, string szImagePathName, unordered_map<string, GMM>& hNewGMMs, int nIter = 1)
	{
		string GMMCopy;
		string szTmpSCP;
		string szTmpMLF; 
		string szCmd;
		string szDelim = "";
		string szCurrFolder = "";
		string szName;
		
		//FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		FindCurrentFileName(szImagePathName, szName);
		//string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		string szDataPaths = szFEFolder + szName;
		//string szmkdir = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string szmkdir = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string szmkdir = "mkdir -p " + szFEFolder;
		#endif

		struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(szmkdir.c_str());
		long hr = PrepareHTKMAPCommands(GMMPaths, szDataPaths, GMMCopy, szTmpSCP, szTmpMLF, szCmd, nIter);
		if(hr != S_OK)
		{
			return hr;
		}
		vector<string> vszLab(2);
		vszLab[0] = "BkGnd";
		vszLab[1] = "Hand";
		hr = PrepareTempDocsForAdapt(vvvRecogVecm, vszLab, szName, szFEFolder, szTmpSCP, szTmpMLF);
		if(hr != S_OK)
		{
			return hr;
		}
		system(szCmd.c_str());
		LoadGMM(GMMCopy, hNewGMMs);
	}

	static long PrepareSMAPCommands(string GMMPaths, string DataPaths, string& GMMCopy, string& szTmpSCP, string& szTmpMLF, string& szCmd)
	{
		szTmpSCP = DataPaths + "_Tmp.scp";
		szTmpMLF = DataPaths + "_Tmp.mlf";
		GMMCopy = GMMPaths + "_Tmp_SMAP";

#ifdef _MSC_VER
		string treeadapt_smap = "";
		cerr<<"Current version didn't support windows OS\n";
		return E_NOTIMPL;
#else
		string treeadapt_smap = "/cluster/users/ycheng31/SMAP/treeadapt_smap/treeadapt_smap";
		string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
#endif
		// Copy GMM as Temp
		//system(szCopyGMM.c_str());
		szCmd = treeadapt_smap;
		szCmd.append(" -i ");
		szCmd.append(GMMPaths);
		szCmd.append(" -o ");
		szCmd.append(GMMCopy);
		szCmd.append(" -m ");
		szCmd.append(szTmpMLF);
		szCmd.append(" -s ");
		szCmd.append(szTmpSCP);
		szCmd.append(" -d 1");
		return 0;
	}

	static long Do_Skin_SMAP(vector< vector< vector<double> > >& vvvRecogVecm, string GMMPaths, string szFEFolder, string szImagePathName, unordered_map<string, GMM>& hNewGMMs)
	{
		string GMMCopy;
		string szTmpSCP;
		string szTmpMLF; 
		string szCmd;
		string szDelim = "";
		string szCurrFolder = "";
		string szName;
		
		//FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		FindCurrentFileName(szImagePathName, szName);
		//string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		string szDataPaths = szFEFolder + szName;
		//string szmkdir = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string szmkdir = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string szmkdir = "mkdir -p " + szFEFolder;
		#endif

		struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(szmkdir.c_str());

		long hr = PrepareSMAPCommands(GMMPaths, szDataPaths, GMMCopy, szTmpSCP, szTmpMLF, szCmd);
		if(hr != S_OK)
		{
			return hr;
		}
		vector<string> vszLab(2);
		vszLab[0] = "BkGnd";
		vszLab[1] = "Hand";
		hr = PrepareTempDocsForAdapt(vvvRecogVecm, vszLab, szName, szFEFolder, szTmpSCP, szTmpMLF);
		if(hr != S_OK)
		{
			return hr;
		}
		system(szCmd.c_str());
		LoadGMM(GMMCopy, hNewGMMs);
		return S_OK;
	}

	
	static long BaselineASeq(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}

		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}

		TestAllFeatures(vvvFeatures, GMMs, vvRecogClasses) ;
		return S_OK;
	}

	static long BaselineASeq_NOBS(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}

		TestAllFeatures(vvvFeatures, GMMs, vvRecogClasses) ;
		//cvReleaseImage(&pBkGnd);
		return S_OK;
	}

	static long BaselineASeq_BS0_Corrs_Only(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}

		TestAllFeatures_BS0_Corrs_Only(vvvFeatures, GMMs, vvRecogClasses, vvvCorrs) ;
		//cvReleaseImage(&pBkGnd);
		return S_OK;
	}




	static long BaselineASeq_BS0(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}

		TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses, vvvCorrs) ;
		//cvReleaseImage(&pBkGnd);
		return S_OK;
	}


	static long HTKGMMAdaptPerFrame_Cheat(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvRecogClassesCheat)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		for(long n = 0; n < NN; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;
			//TestAFrameFeatures(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n);
			// Cheating here:
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= GMMITER; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAFrameFeatures(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n);
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			//TestAFrameFeatures_BS0(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n, vvvCorrs);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		return S_OK;
	}



	static long HTKGMMAdaptPerFrame(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		for(long n = 0; n < NN; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= GMMITER; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAFrameFeatures(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n);
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		return S_OK;
	}

	static long HTKGMMAdaptPerFrame(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvxyrange)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		for(long n = 0; n < NN; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;

			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= GMMITER; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAFrameFeatures(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n);
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static long HTKGMMAdaptASeq(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;

		vector<GMM> GMMsCopy;
		GMMsCopy = GMMs;
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		for(int nIter = 1; nIter <= GMMITER; nIter++)
		{
			vvvRecogVec0.clear();
			hNewGMMs.clear();
			newGMMs.clear();
			TestAllFeatures(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0) ;
			Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			GMMsCopy.clear();
			GMMsCopy = newGMMs;
		}

		//TestAllFeatures(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0) ;
		//unordered_map<string, GMM> hNewGMMs;
		//vector<GMM> newGMMs;
		//string szSeqName, szCurrPath;
		//FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		//string szTmpName = szCurrPath + szSeqName;
		//Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs);
		//GMMHash2VectBH(hNewGMMs, newGMMs);
		TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		return S_OK;
	}

	static long HTKGMMAdaptASeq(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvxyrange)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;

		vector<GMM> GMMsCopy;
		GMMsCopy = GMMs;
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		for(int nIter = 1; nIter <= GMMITER; nIter++)
		{
			vvvRecogVec0.clear();
			hNewGMMs.clear();
			newGMMs.clear();
			TestAllFeatures(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0) ;
			string szSeqName, szCurrPath;
			FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
			string szTmpName = szCurrPath + szSeqName;
			Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			GMMsCopy.clear();
			GMMsCopy = newGMMs;
		}
		TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static long HTKGMMSMAPPerFrame(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		for(long n = 0; n < NN; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;
			TestAFrameFeatures(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n);
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		return S_OK;
	}

	static long HTKGMMSMAPPerFrame(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvxyrange)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		for(long n = 0; n < NN; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;
			TestAFrameFeatures(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n);
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static long HTKGMMMAPPerFrame_BS0_Cheat(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq, string& GMMPath, vector< vector<int> >& vvRecogClassesCheat)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		string szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		for(long n = nStartSeq; n < nStopSeq; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;

			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= GMMITER; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAFrameFeatures_BS0_Cheat(vvvFeatures, GMMsCopy, vvRecogClassesCheat[n], vvvRecogVec0, n, vvvCorrs);
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			TestAFrameFeatures_BS0(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n, vvvCorrs);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		//TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses, vvvCorrs) ;
		//cvReleaseImage(&pBkGnd);
		return S_OK;
	}




	static long HTKGMMMAPPerFrame_BS0_Scores(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector< vector<double> > >& vvvRecogScores, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq, string& GMMPath)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		string szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		for(long n = nStartSeq; n < nStopSeq; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;
			//TestAFrameFeatures_BS0_V2(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;


			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			for(int nIter = 1; nIter <= GMMITER; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAFrameFeatures_BS0_V2(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
				//TestAFrameFeatures_BS0(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			//TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			if(GMMITER < 1)
			{
				newGMMs = GMMs;
			}
			vector< vector<double> > vvScores;
			TestAFrameFeatures_BS0_Scores(vvvFeatures, newGMMs, vvScores, n, vvvCorrs);
			vvvRecogScores.push_back(vvScores);
			//TestAFrameFeatures_BS0(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n, vvvCorrs);
			//vvRecogClasses.push_back(vRecogClasses1);
		}
		//TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses, vvvCosrrs) ;
		//cvReleaseImage(&pBkGnd);
		return S_OK;
	}


	static long HTKGMMMAPPerFrame_BS0(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq, string& GMMPath)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		string szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		for(long n = nStartSeq; n < nStopSeq; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;
			//TestAFrameFeatures_BS0_V2(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;


			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			for(int nIter = 1; nIter <= GMMITER; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAFrameFeatures_BS0_V2(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			//TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			TestAFrameFeatures_BS0(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n, vvvCorrs);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		//TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses, vvvCorrs) ;
		//cvReleaseImage(&pBkGnd);
		return S_OK;
	}

	static long HTKGMMMAPASeq_BS0_Scores(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector< vector<double> > >& vvvRecogScores, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq, string& GMMPath)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		
		vector< vector<int> > vvRecogClasses0;
		vector< vector< vector<double> > > vvvRecogVec0;
		//vector< vector<int> > vvRecogVec0;
		//TestAllFeatures_BS0_V2(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
		//cvReleaseImage(&pBkGnd);
		
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
        	string szFolder, szDelim;
	        FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
        	string szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
	        //string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
		        system(mkdirFolder.c_str());

		vector<GMM> GMMsCopy;
		GMMsCopy = GMMs;
		//unordered_map<string, GMM> hNewGMMs;
		//vector<GMM> newGMMs;
		for(int nIter = 1; nIter <= GMMITER; nIter++)
		{
			vvvRecogVec0.clear();
			hNewGMMs.clear();
			newGMMs.clear();
			TestAllFeatures_BS0_V2(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
			//TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
			Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			GMMsCopy.clear();
			GMMsCopy = newGMMs;
		}
		if(GMMITER < 1)
		{
			newGMMs = GMMs;
		}
		//TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		//TestAllFeatures_BS0(vvvFeatures, newGMMs, vvRecogClasses, vvvCorrs);
		AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs, vvKSs);

		return S_OK;
	}

	static long HTKGMMMAPASeq_SLIC_BS0_Scores(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector< vector<double> > >& vvvRecogScores, vector<string>& vszLoadadFiles, int& nStartSeq, int& nStopSeq, string& GMMPath)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtSLICColorFeaturesFromImagesStartStop(vszLoadadFiles, vImgs, vvvFeatures, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		
		vector< vector<int> > vvRecogClasses0;
		vector< vector< vector<double> > > vvvRecogVec0;
		
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		string szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector<GMM> GMMsCopy;
		GMMsCopy = GMMs;
		//unordered_map<string, GMM> hNewGMMs;
		//vector<GMM> newGMMs;
		for(int nIter = 1; nIter <= GMMITER; nIter++)
		{
			vvvRecogVec0.clear();
			hNewGMMs.clear();
			newGMMs.clear();
			TestAllFeatures_BS0_V2(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
			//TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
			Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			GMMsCopy.clear();
			GMMsCopy = newGMMs;
		}
		if(GMMITER < 1)
		{
			newGMMs = GMMs;
		}
		AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs, vvKSs);

		return S_OK;
	}






	static long HTKGMMMAPASeq_BS0(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq, string& GMMPath)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		
		vector< vector<int> > vvRecogClasses0;
		vector< vector< vector<double> > > vvvRecogVec0;
		//vector< vector<int> > vvRecogVec0;
		//TestAllFeatures_BS0_V2(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
		//cvReleaseImage(&pBkGnd);
		
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
	        string szFolder, szDelim;    
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		string szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector<GMM> GMMsCopy;
		GMMsCopy = GMMs;
		//unordered_map<string, GMM> hNewGMMs;
		//vector<GMM> newGMMs;
		for(int nIter = 1; nIter <= GMMITER; nIter++)
		{
			vvvRecogVec0.clear();
			hNewGMMs.clear();
			newGMMs.clear();
			TestAllFeatures_BS0_V2(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
			Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			GMMsCopy.clear();
			GMMsCopy = newGMMs;
		}
		//TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		TestAllFeatures_BS0(vvvFeatures, newGMMs, vvRecogClasses, vvvCorrs);

		return S_OK;
	}





	static long HTKGMMSMAPPerFrame_BS0(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq, string& GMMPath)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		string szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		for(long n = nStartSeq; n < nStopSeq; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;
			TestAFrameFeatures_BS0(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			//TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			TestAFrameFeatures_BS0(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n, vvvCorrs);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		//TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses, vvvCorrs) ;
		//cvReleaseImage(&pBkGnd);
		return S_OK;
	}


	//static long HTKGMMSMAPPerFrame_BS0(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvxyrange, cv::Mat& ref, vector< vector< vector<double> > >& vvvCorrs)
	static long HTKGMMSMAPPerFrame_BS0(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvxyrange, cv::Mat& ref, vector< vector< vector<double> > >& vvvCorrs, vector< vector<double> >& vvKSs)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		//hr = ExtColorFeaturesFromImages_BS0(vImgs, vvvFeatures, BLOCKSIZE, ref, vvvCorrs); 
		hr = ExtColorFeaturesFromImages_BS0(vImgs, vvvFeatures, BLOCKSIZE, ref, vvvCorrs, vvKSs); 
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		for(long n = 0; n < NN; n++)
		{
			vector<int> vRecogClasses0;
			vector<int> vRecogClasses1;
			vector< vector< vector<double> > > vvvRecogVec0;
			vector< vector< vector<double> > > vvvRecogVec1;
			TestAFrameFeatures_BS0(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs);
			GMMHash2VectBH(hNewGMMs, newGMMs);
			//TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
			TestAFrameFeatures_BS0(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n, vvvCorrs);
			vvRecogClasses.push_back(vRecogClasses1);
		}
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static long HTKGMMSMAPASeq(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;
		TestAllFeatures(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0) ;
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs);
		GMMHash2VectBH(hNewGMMs, newGMMs);
		TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		return S_OK;
	}

	static long HTKGMMSMAPASeq(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvxyrange)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		hr = ExtColorFeaturesFromImages(vImgs, vvvFeatures, BLOCKSIZE);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;
		TestAllFeatures(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0) ;
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs);
		GMMHash2VectBH(hNewGMMs, newGMMs);
		TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static long HTKGMMSMAPASeq_BS0(vector<GMM>& GMMs, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, int& nStartSeq, int& nStopSeq, string& GMMPath)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = 0;
		int nStop = 0;
		TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		nStartSeq = nDetStopBkGnd + 1;
		nStopSeq = nStop;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		
		vector< vector<int> > vvRecogClasses0;
		vector< vector< vector<double> > > vvvRecogVec0;
		//vector< vector<int> > vvRecogVec0;
		TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
		//cvReleaseImage(&pBkGnd);
		
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;

	        string szFolder, szDelim;
                FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
                string szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
                //string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
                	system(mkdirFolder.c_str());
	
		Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs);
		GMMHash2VectBH(hNewGMMs, newGMMs);
		//TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		TestAllFeatures_BS0(vvvFeatures, newGMMs, vvRecogClasses, vvvCorrs);
		return S_OK;
	}

	//static long HTKGMMSMAPASeq_BS0(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvxyrange, cv::Mat& ref, vector< vector< vector<double> > >& vvvCorrs)
	static long HTKGMMSMAPASeq_BS0(vector<GMM>& GMMs, string& GMMPath, string& szFEFolder, vector<string>& vszImagePathNames, vector< vector<int> >& vvRecogClasses, vector<string>& vszLoadadFiles, int BLOCKSIZE, vector< vector<int> >& vvxyrange, cv::Mat& ref, vector< vector< vector<double> > >& vvvCorrs, vector< vector<double> >& vvKSs)
	{
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		//hr = ExtColorFeaturesFromImages_BS0(vImgs, vvvFeatures, BLOCKSIZE, ref, vvvCorrs);
		hr = ExtColorFeaturesFromImages_BS0(vImgs, vvvFeatures, BLOCKSIZE, ref, vvvCorrs, vvKSs);
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;
		TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0, vvvCorrs);
	
		unordered_map<string, GMM> hNewGMMs;
		vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs);
		GMMHash2VectBH(hNewGMMs, newGMMs);
		//TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		TestAllFeatures_BS0(vvvFeatures, newGMMs, vvRecogClasses, vvvCorrs);
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}



	static long HTKGMMMAP_BS0_Fusion(
		vector<GMM>& GMMs, 
		string& GMMPath, 
		string& szFEFolder, 
		vector<string>& vszImagePathNames, 
		vector< vector<int> >& vvRecogClasses, 
		vector<string>& vszLoadadFiles, 
		int BLOCKSIZE, 
		vector< vector<int> >& vvxyrange, 
		cv::Mat& ref, 
		vector<GMM>& vFusionGMMs, 
		bool bFrameBatch, 
		int nAdaptMode, 
		int nStartSeq,
		int nStopSeq
		)
	{
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
//		hr = ExtColorFeaturesFromImages_BS0(vImgs, vvvFeatures, BLOCKSIZE, ref, vvvCorrs);
		//int nDetStopBkGnd = 0;
		//int nStop = 0;
		//TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		//int nStartSeq = nDetStopBkGnd + 1;
		//int nStopSeq = nStop;
		int nDetStopBkGnd = nStartSeq - 1;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		// background frames with null vectors!
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;
		//TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0, vvvCorrs);
	
		//unordered_map<string, GMM> hNewGMMs;
		//vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		vector< vector< vector<double> > > vvvRecogScores;
		int nNumIter = GMMITER;
		if(nAdaptMode == 1)
			nNumIter = 0;
		if(bFrameBatch)
		{
			// Framewise adaptation
			for(long n = 0; n < vvvFeatures.size(); n++)
			{
				if(n >= nStartSeq && n < nStopSeq)
				{
					vector<int> vRecogClasses0;
					vector<int> vRecogClasses1;
					vector< vector< vector<double> > > vvvRecogVec0;
					vector< vector< vector<double> > > vvvRecogVec1;
					//TestAFrameFeatures_BS0_V2(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
					unordered_map<string, GMM> hNewGMMs;
					vector<GMM> newGMMs;

					vector<GMM> GMMsCopy;
					GMMsCopy = GMMs;

					for(int nIter = 1; nIter <= nNumIter; nIter++)
					{
						vvvRecogVec0.clear();
						hNewGMMs.clear();
						newGMMs.clear();
						//TestAFrameFeatures_BS0_V2(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						TestAFrameFeatures_BS0(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
						GMMHash2VectBH(hNewGMMs, newGMMs);
						GMMsCopy.clear();
						GMMsCopy = newGMMs;
					}
					//TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
					if(nNumIter < 1)
					{
						newGMMs = GMMs;
					}
					vector< vector<double> > vvScores;
					TestAFrameFeatures_BS0_Scores(vvvFeatures, newGMMs, vvScores, n, vvvCorrs);
					vvvRecogScores.push_back(vvScores);
				}
				else
				{
					vector< vector<double> > vvScores;
					vvvRecogScores.push_back(vvScores);
				}
			}

		}
		else
		{
			// Batch adaptation
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= nNumIter; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				//TestAllFeatures_BS0_V2(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			if(nNumIter < 1)
			{
				newGMMs = GMMs;
			}
			AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs);
		}
		vvRecogClasses.resize(vvvFeatures.size());
		if(nAdaptMode == 0)
		{
			// MAP only
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvRecogClasses[i].resize(vvvRecogScores[i].size(), 0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						if(vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0] > 0.0)
						{
							vvRecogClasses[i][j] = 1;
						}
					}
				}
			}
		}
		else
		{
			// Score fusions
			vector< vector< vector<double> > > vvvScoreForFusions(vvvFeatures.size());
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvvScoreForFusions[i].resize(vvvRecogScores[i].size());
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						double LLR = vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0];
						double sLLR = 1.0/(1.0+exp(-0.5*LLR));
						double sCorr = 0.1140*vvvCorrs[i][j][0] + 0.5870 * vvvCorrs[i][j][1] + 0.2989*vvvCorrs[i][j][2];
						//cout<<LLR<<"\t"<<sLLR<<"\t"<<sCorr<<endl;
						vvvScoreForFusions[i][j].push_back(sLLR);
						vvvScoreForFusions[i][j].push_back(sCorr);
					}
				}
			}
			TestAllFeatures(vvvScoreForFusions, vFusionGMMs, vvRecogClasses);
		}

		//Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs);
		//GMMHash2VectBH(hNewGMMs, newGMMs);
		//TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		//TestAllFeatures_BS0(vvvFeatures, newGMMs, vvRecogClasses, vvvCorrs);
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static long HTKGMMMAP_BS0_Fusion_With_Scores(
		vector<GMM>& GMMs, 
		string& GMMPath, 
		string& szFEFolder, 
		vector<string>& vszImagePathNames, 
		vector< vector<int> >& vvRecogClasses, 
		vector<string>& vszLoadadFiles, 
		int BLOCKSIZE, 
		vector< vector<int> >& vvxyrange, 
		cv::Mat& ref, 
		vector<GMM>& vFusionGMMs, 
		bool bFrameBatch, 
		int nAdaptMode, 
		int nStartSeq,
		int nStopSeq,
		vector< vector<double> >& vvScore0, 
		vector< vector<double> >& vvScore1,
		vector< vector< vector<double> > >& vvvScoreForFusions
		)
	{
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
//		hr = ExtColorFeaturesFromImages_BS0(vImgs, vvvFeatures, BLOCKSIZE, ref, vvvCorrs);
		//int nDetStopBkGnd = 0;
		//int nStop = 0;
		//TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		//int nStartSeq = nDetStopBkGnd + 1;
		//int nStopSeq = nStop;
		int nDetStopBkGnd = nStartSeq - 1;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		// background frames with null vectors!
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;
		//TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0, vvvCorrs);
	
		//unordered_map<string, GMM> hNewGMMs;
		//vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		vector< vector< vector<double> > > vvvRecogScores;
		int nNumIter = GMMITER;
		if(nAdaptMode == 1)
			nNumIter = 0;
		if(bFrameBatch)
		{
			// Framewise adaptation
			for(long n = 0; n < vvvFeatures.size(); n++)
			{
				if(n >= nStartSeq && n < nStopSeq)
				{
					vector<int> vRecogClasses0;
					vector<int> vRecogClasses1;
					vector< vector< vector<double> > > vvvRecogVec0;
					vector< vector< vector<double> > > vvvRecogVec1;
					//TestAFrameFeatures_BS0_V2(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
					unordered_map<string, GMM> hNewGMMs;
					vector<GMM> newGMMs;

					vector<GMM> GMMsCopy;
					GMMsCopy = GMMs;

					for(int nIter = 1; nIter <= nNumIter; nIter++)
					{
						vvvRecogVec0.clear();
						hNewGMMs.clear();
						newGMMs.clear();
						//TestAFrameFeatures_BS0_V2(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						TestAFrameFeatures_BS0(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
						GMMHash2VectBH(hNewGMMs, newGMMs);
						GMMsCopy.clear();
						GMMsCopy = newGMMs;
					}
					//TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
					if(nNumIter < 1)
					{
						newGMMs = GMMs;
					}
					vector< vector<double> > vvScores;
					TestAFrameFeatures_BS0_Scores(vvvFeatures, newGMMs, vvScores, n, vvvCorrs);
					vvvRecogScores.push_back(vvScores);
				}
				else
				{
					vector< vector<double> > vvScores;
					vvvRecogScores.push_back(vvScores);
				}
			}

		}
		else
		{
			// Batch adaptation
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= nNumIter; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				//TestAllFeatures_BS0_V2(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			if(nNumIter < 1)
			{
				newGMMs = GMMs;
			}
			AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs);
		}
		vvRecogClasses.resize(vvvFeatures.size());
		vvScore0.resize(vvvFeatures.size());
		vvScore1.resize(vvvFeatures.size());
		if(nAdaptMode == 0)
		{
			// MAP only
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvRecogClasses[i].resize(vvvRecogScores[i].size(), 0);
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						vvScore0[i][j] = vvvRecogScores[i][j][0];
						vvScore1[i][j] = vvvRecogScores[i][j][1];
						if(vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0] > 0.0)
						{
							vvRecogClasses[i][j] = 1;
						}
					}
				}
			}
		}
		else
		{
			// Score fusions
			//vector< vector< vector<double> > > vvvScoreForFusions(vvvFeatures.size());
			vvvScoreForFusions.resize(vvvFeatures.size());
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvvScoreForFusions[i].resize(vvvRecogScores[i].size());
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						double LLR = vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0];
						double sLLR = 1.0/(1.0+exp(-0.5*LLR));
						double sCorr = 0.1140*vvvCorrs[i][j][0] + 0.5870 * vvvCorrs[i][j][1] + 0.2989*vvvCorrs[i][j][2];
						//cout<<LLR<<"\t"<<sLLR<<"\t"<<sCorr<<endl;
						vvvScoreForFusions[i][j].push_back(sLLR);
						vvvScoreForFusions[i][j].push_back(sCorr);
					}
				}
			}
			vector< vector< vector<double> > > vvvRecogFusionScores;
			TestAllFeatures_Scores(vvvScoreForFusions, 
								vFusionGMMs, 
								vvRecogClasses,
								vvvRecogFusionScores
								);
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
				{
					vvScore0[i][j] = vvvRecogFusionScores[i][j][0];
					vvScore1[i][j] = vvvRecogFusionScores[i][j][1];
				}
			}
			//TestAllFeatures(vvvScoreForFusions, vFusionGMMs, vvRecogClasses);
		}

		//Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs);
		//GMMHash2VectBH(hNewGMMs, newGMMs);
		//TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		//TestAllFeatures_BS0(vvvFeatures, newGMMs, vvRecogClasses, vvvCorrs);
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}


	static long HTKGMMMAP_BS0_Fusion_With_Scores(
		vector<GMM>& GMMs, 
		string& GMMPath, 
		string& szFEFolder, 
		vector<string>& vszImagePathNames, 
		vector< vector<int> >& vvRecogClasses, 
		vector<string>& vszLoadadFiles, 
		int BLOCKSIZE, 
		vector< vector<int> >& vvxyrange, 
		cv::Mat& ref, 
		vector<GMM>& vFusionGMMs, 
		bool bFrameBatch, 
		int nAdaptMode, 
		int nStartSeq,
		int nStopSeq,
		vector< vector<double> >& vvScore0, 
		vector< vector<double> >& vvScore1
		)
	{
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
//		hr = ExtColorFeaturesFromImages_BS0(vImgs, vvvFeatures, BLOCKSIZE, ref, vvvCorrs);
		//int nDetStopBkGnd = 0;
		//int nStop = 0;
		//TestOnlineBkgndFrames1Dot2(vImgs, nDetStopBkGnd, nStop);
		//int nStartSeq = nDetStopBkGnd + 1;
		//int nStopSeq = nStop;
		int nDetStopBkGnd = nStartSeq - 1;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs);
		// background frames with null vectors!
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;
		//TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0, vvvCorrs);
	
		//unordered_map<string, GMM> hNewGMMs;
		//vector<GMM> newGMMs;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		vector< vector< vector<double> > > vvvRecogScores;
		int nNumIter = GMMITER;
		if(nAdaptMode == 1)
			nNumIter = 0;
		if(bFrameBatch)
		{
			// Framewise adaptation
			for(long n = 0; n < vvvFeatures.size(); n++)
			{
				if(n >= nStartSeq && n < nStopSeq)
				{
					vector<int> vRecogClasses0;
					vector<int> vRecogClasses1;
					vector< vector< vector<double> > > vvvRecogVec0;
					vector< vector< vector<double> > > vvvRecogVec1;
					//TestAFrameFeatures_BS0_V2(vvvFeatures, GMMs, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
					unordered_map<string, GMM> hNewGMMs;
					vector<GMM> newGMMs;

					vector<GMM> GMMsCopy;
					GMMsCopy = GMMs;

					for(int nIter = 1; nIter <= nNumIter; nIter++)
					{
						vvvRecogVec0.clear();
						hNewGMMs.clear();
						newGMMs.clear();
						//TestAFrameFeatures_BS0_V2(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						TestAFrameFeatures_BS0(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
						GMMHash2VectBH(hNewGMMs, newGMMs);
						GMMsCopy.clear();
						GMMsCopy = newGMMs;
					}
					//TestAFrameFeatures(vvvFeatures, newGMMs, vRecogClasses1, vvvRecogVec1, n);
					if(nNumIter < 1)
					{
						newGMMs = GMMs;
					}
					vector< vector<double> > vvScores;
					TestAFrameFeatures_BS0_Scores(vvvFeatures, newGMMs, vvScores, n, vvvCorrs);
					vvvRecogScores.push_back(vvScores);
				}
				else
				{
					vector< vector<double> > vvScores;
					vvvRecogScores.push_back(vvScores);
				}
			}

		}
		else
		{
			// Batch adaptation
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= nNumIter; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				//TestAllFeatures_BS0_V2(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			if(nNumIter < 1)
			{
				newGMMs = GMMs;
			}
			AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs);
		}
		vvRecogClasses.resize(vvvFeatures.size());
		vvScore0.resize(vvvFeatures.size());
		vvScore1.resize(vvvFeatures.size());
		if(nAdaptMode == 0)
		{
			// MAP only
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvRecogClasses[i].resize(vvvRecogScores[i].size(), 0);
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						vvScore0[i][j] = vvvRecogScores[i][j][0];
						vvScore1[i][j] = vvvRecogScores[i][j][1];
						if(vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0] > 0.0)
						{
							vvRecogClasses[i][j] = 1;
						}
					}
				}
			}
		}
		else
		{
			// Score fusions
			vector< vector< vector<double> > > vvvScoreForFusions(vvvFeatures.size());
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvvScoreForFusions[i].resize(vvvRecogScores[i].size());
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						double LLR = vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0];
						double sLLR = 1.0/(1.0+exp(-0.5*LLR));
						double sCorr = 0.1140*vvvCorrs[i][j][0] + 0.5870 * vvvCorrs[i][j][1] + 0.2989*vvvCorrs[i][j][2];
						//cout<<LLR<<"\t"<<sLLR<<"\t"<<sCorr<<endl;
						vvvScoreForFusions[i][j].push_back(sLLR);
						vvvScoreForFusions[i][j].push_back(sCorr);
					}
				}
			}
			vector< vector< vector<double> > > vvvRecogFusionScores;
			TestAllFeatures_Scores(vvvScoreForFusions, 
								vFusionGMMs, 
								vvRecogClasses,
								vvvRecogFusionScores
								);
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
				{
					vvScore0[i][j] = vvvRecogFusionScores[i][j][0];
					vvScore1[i][j] = vvvRecogFusionScores[i][j][1];
				}
			}
			//TestAllFeatures(vvvScoreForFusions, vFusionGMMs, vvRecogClasses);
		}

		//Do_Skin_SMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs);
		//GMMHash2VectBH(hNewGMMs, newGMMs);
		//TestAllFeatures(vvvFeatures, newGMMs, vvRecogClasses);
		//TestAllFeatures_BS0(vvvFeatures, newGMMs, vvRecogClasses, vvvCorrs);
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}





	static int GMMByLAB(string GMMPaths, string DataLABPaths, double& Acc)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		
		int BLOCKSIZE = 16;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName0 = DataLABPaths + "_Orig01";
		LoadJPGLabelInfo lb, lb0;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			BaselineASeq(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;
			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb0.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		lb0.SaveJPGInfo(szLabelPathName0);
		return 0;
//		unordered_map<string, GMM> hGMMs;
//		vector<GMM> GMMs;
//		GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
//		GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
//
//		// Get Current Folder
//		int nPos0001 = DataLABPaths.find_last_of("/");
//		int nPob0001 = DataLABPaths.find_last_of("\\");
//		int nPos01 = (nPob0001>nPos0001)?nPob0001:nPos0001;
//		string szCurrFolder = DataLABPaths.substr(0, nPos01);
//
//		string szDelim = "";
//		if(nPos01 == nPob0001) szDelim = "\\";
//		else if(nPos01 == nPos0001) szDelim = "/";
//		
//		// Create FE store folder
//		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
//		string sysMkDirFE = "mkdir -p " + szFEFolder;
		//#ifdef _MSC_VER
		//	string sysMkDirFE = "mkdir " + szFEFolder;
		//#elif defined(__GNUC__)
		//	string sysMkDirFE = "mkdir -p " + szFEFolder;
		//#endif
//		system(sysMkDirFE.c_str());
//		string szTmp = DataLABPaths + "_Tmp";
//		string szTmpMLF = DataLABPaths + "_Tmp.mlf";
//		string GMMCopy = GMMPaths + "_Tmp";
//		
//		// Output Label
//		string szLabelPathName = DataLABPaths;
//		string szLabelPathName1 = DataLABPaths + "_MAPPFHTK";
//		string szLabelPathName0 = DataLABPaths + "_Orig01";
//
//		Acc = 0.0;
//		double Acc0 = 0.0;
//		double N_1 = 0;
//		int BLKSIZE = 16;
//		// int nCountX = 0;
//		// int nCountO = 0;
//		unordered_map<string, int> hszJPGFiles;
//		unordered_map<string, int> hszJPGFiles1;
//		unordered_map<string, int> hszJPGFiles0;
//
////#ifdef _WIN32
////		string szCopyGMM = "copy " + GMMPaths + " " + GMMCopy;
////#else
////		string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
////#endif
//
//		// Now, load label
//		LoadJPGLabelInfo lb, lb0, lb1;
//		lb.m_szLabelPathName = szLabelPathName;
//		lb.LoadJPGInfo();
//
//		// Loop through all files included in the label file:
//		for(size_t n = 0; n < lb.m_vszKeyOrder.size(); n++)
//		{
//			string szCurrFilePath = lb.m_vszKeyOrder[n];
//
//			int nPos001s = szCurrFilePath.find_last_of("/");
//			int nPos001b = szCurrFilePath.find_last_of("\\");
//			int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
//			int nPos002 = szCurrFilePath.find_last_of(".");
//			string szCurrFileName = szCurrFilePath.substr(nPos001 + 1, nPos002 - nPos001 - 1);
//			string szCurrFE1 = szFEFolder + szCurrFileName +"_Hand.fe";
//			string szCurrFE0 = szFEFolder + szCurrFileName +"_BkGnd.fe";
//
//			JPGLabelInfo* pinfo = &lb.m_hJPGLavelInfo[szCurrFilePath];
//			vector< vector<double> > vvCurrFeatures;
//			//ExtColorFeatureRGBOnly(szCurrFilePath, pinfo, vvCurrFeatures, BLKSIZE);
//			ExtLowLevelImgFeature::ExtColorFeatureRGBOnly(szCurrFilePath, vvCurrFeatures, BLKSIZE);
//			//ExtColorFeature(szCurrFilePath, pinfo, vvCurrFeatures, BLKSIZE);
//			int nHeight = pinfo->m_nBlockHeight;
//			int nWidth = pinfo->m_nBlockWidth;
//			vector<int> vCurrRecogClasses0, vCurrRecogClasses;
//			double Accn0 = 0.0;
//			double Accn = 0.0;
//			vector<int> confres0, confres;
//			Test2ClassesGMMs(vvCurrFeatures, lb.m_hJPGClasses[szCurrFilePath], GMMs, vCurrRecogClasses0, Accn0, confres0);
//			vector< vector<double> > vvRecog0, vvRecog1;
//			// Collect those recognized as background and hand, respectively
//			for(int k1 = 0; k1 < nWidth; k1++)
//			{
//				for(int k2 = 0; k2 < nHeight; k2++)
//				{
//					int nCurrIndex = k1 * nHeight + k2;
//					if(vCurrRecogClasses0[nCurrIndex] == 0)
//					{
//						vvRecog0.push_back(vvCurrFeatures[nCurrIndex]);
//					}
//					else if(vCurrRecogClasses0[nCurrIndex] == 1)
//					{
//						vvRecog1.push_back(vvCurrFeatures[nCurrIndex]);
//						lb0.MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1, nWidth, nHeight);
//						if(lb0.m_hszJPGFiles.find(szCurrFilePath) == lb0.m_hszJPGFiles.end())
//						{
//							lb0.m_vszKeyOrder.push_back(szCurrFilePath);
//						}					
//						lb0.m_hszJPGFiles[szCurrFilePath] = 1;
//					}
//				}
//			}
//			double NN = (double)(nWidth * nHeight);
//			Acc0 = (Accn0 * NN + Acc0 * N_1)/(N_1+NN);
//			N_1 += NN;
//		}
//		cout<<"Original result is "<<Acc0<<endl;
//		//lb1.SaveJPGInfo(szLabelPathName1);
//		lb0.SaveJPGInfo(szLabelPathName0);
//		return 0;
	}

	static int GMMByLAB(string GMMPaths, string DataLABPaths, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		
		int BLOCKSIZE = 16;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName0 = DataLABPaths + "_Orig01";
		LoadJPGLabelInfo lb, lb0;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			BaselineASeq(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;
			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					//cout<<"inside loop "<<mm<<endl;
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb0.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb0.SaveJPGInfo(szLabelPathName0);
		return 0;
	}


	static int GMMByLAB_NOBS(string GMMPaths, string DataLABPaths, string szRefSCP, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		
		int BLOCKSIZE = 16;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName0 = DataLABPaths + "_Orig01";
		LoadJPGLabelInfo lb, lb0;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		
		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);

		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			int nStartSeq =0;
			int nStopSeq = 0;
			BaselineASeq_NOBS(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;
			// Get differences:
			for(size_t nn = nStartSeq; nn < nStopSeq/*hvszKeyNames[vszKeys[n]].size()*/; nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					//cout<<"inside loop "<<mm<<endl;
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb0.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb0.SaveJPGInfo(szLabelPathName0);
		return 0;
	}

	static int GMMByLAB_BS0_Corrs_Only(string GMMPaths, string DataLABPaths, string szRefSCP, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		
		int BLOCKSIZE = 16;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName0 = DataLABPaths + "_Orig01";
		LoadJPGLabelInfo lb, lb0;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		
		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);

		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			int nStartSeq =0;
			int nStopSeq = 0;
			BaselineASeq_BS0_Corrs_Only(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;
			// Get differences:
			for(size_t nn = nStartSeq; nn < nStopSeq/*hvszKeyNames[vszKeys[n]].size()*/; nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					//cout<<"inside loop "<<mm<<endl;
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb0.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb0.SaveJPGInfo(szLabelPathName0);
		return 0;
	}




	static int GMMByLAB_BS0(string GMMPaths, string DataLABPaths, string szRefSCP, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		
		int BLOCKSIZE = 16;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName0 = DataLABPaths + "_Orig01";
		LoadJPGLabelInfo lb, lb0;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		
		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);

		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			int nStartSeq =0;
			int nStopSeq = 0;
			BaselineASeq_BS0(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;
			// Get differences:
			for(size_t nn = nStartSeq; nn < nStopSeq/*hvszKeyNames[vszKeys[n]].size()*/; nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					//cout<<"inside loop "<<mm<<endl;
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb0.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb0.SaveJPGInfo(szLabelPathName0);
		return 0;
	}


	static int GMMHTKMAPByLAB_Cheat(string GMMPaths, string DataLABPaths, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_MAPPF";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector< vector<int> > vvRecogClassesCheat(hvszKeyNames[vszKeys[n]].size());
			vector<string> vszLoadadFiles;
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
				size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
				size_t hhww = hh*ww;
				vvRecogClassesCheat[nn].reserve(hhww);
				vvRecogClassesCheat[nn].resize(hhww);
				for(size_t mm = 0; mm < hhww; mm++)
				{
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					vvRecogClassesCheat[nn][mm] = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
				}
			}


			HTKGMMAdaptPerFrame_Cheat(GMMs, GMMPaths, szFEFolder, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, vvRecogClassesCheat);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					//int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					int nclass = vvRecogClassesCheat[nn][mm];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;
	}



	static int GMMHTKMAPByLAB_BS0_Cheat(string GMMPaths, string DataLABPaths, string szRefSCP, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_SMAPPF";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();

		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);

		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;

			vector< vector<int> > vvRecogClassesCheat(hvszKeyNames[vszKeys[n]].size());
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
				size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
				size_t hhww = hh*ww;
				vvRecogClassesCheat[nn].reserve(hhww);
				vvRecogClassesCheat[nn].resize(hhww);
				for(size_t mm = 0; mm < hhww; mm++)
				{
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					vvRecogClassesCheat[nn][mm] = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
				}
			}


			int nStartSeq = 0;
			int nStopSeq = 0;
			HTKGMMMAPPerFrame_BS0_Cheat(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq, GMMPaths, vvRecogClassesCheat);

			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = nStartSeq; nn < nStopSeq; nn++) //hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn - nStartSeq].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn - nStartSeq][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn - nStartSeq][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;
	}

	static int GMMHTKMAPByLAB_BS0_Scores(string& GMMPaths, string& DataLABPaths, string& szRefSCP, string& szOutData)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_SMAPPF";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();

		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);

		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		FILE* pfOut = NULL;
		//string tmpFolder, tmpDelim;
		//GMMAdaptTool::FindCurrentFolderANDDelim(szOutData, tmpFolder, tmpDelim);
		//CVideoFE tmpfe;
		//tmpfe.MkdirHier(tmpFolder);
		pfOut = fopen(szOutData.c_str(), "wb");
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			vector< vector< vector<double> > > vvvRecogScores;
			int nStartSeq = 0;
			int nStopSeq = 0;
			HTKGMMMAPPerFrame_BS0_Scores(GMMs, hvszKeyNames[vszKeys[n]], vvvRecogScores, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq, GMMPaths);
			for(size_t n1 = 0; n1 < vvvRecogScores.size(); ++n1)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][n1 + nStartSeq];
				size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
				size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
				for(size_t n2 = 0; n2 < vvvRecogScores[n1].size(); ++n2)
				{
					size_t mm2 =  n2 / hh;
					size_t mm1 = n2 - (mm2 * hh);
					double dbclass = (double)lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					fwrite(&dbclass, sizeof(double), 1, pfOut);
					for(size_t n3 = 0; n3 < vvvRecogScores[n1][n2].size(); ++n3)
					{
						double db = vvvRecogScores[n1][n2][n3];
						fwrite(&db, sizeof(double), 1, pfOut);
					}
				}
			}
			cout<<"Sequence "<<n<<" is computed!"<<endl;
		}
		fclose(pfOut);

		return 0;
	}

	static int GMMHTKMAPByLAB2_BS0_Scores(string& GMMPaths, string& DataLABPaths, string& szRefSCP, string& szOutData)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_MAPAL";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();

		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);
		
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		FILE* pfOut = NULL;
		//string tmpFolder, tmpDelim;
		//GMMAdaptTool::FindCurrentFolderANDDelim(szOutData, tmpFolder, tmpDelim);
		//CVideoFE tmpfe;
		//tmpfe.MkdirHier(tmpFolder);
		pfOut = fopen(szOutData.c_str(), "wb");
		for(long n = 0; n < NN; n++)
		{
			vector< vector< vector<double> > > vvvRecogScores;
			vector<string> vszLoadadFiles;
			int nStartSeq = 0;
			int nStopSeq = 0;
			HTKGMMMAPASeq_BS0_Scores(GMMs, hvszKeyNames[vszKeys[n]], vvvRecogScores, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq, GMMPaths);
			for(size_t n1 = 0; n1 < vvvRecogScores.size(); ++n1)
			{
				if(vvvRecogScores[n1].size())
				{
					string szCurrFileName = hvszKeyNames[vszKeys[n]][n1];
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					for(size_t n2 = 0; n2 < vvvRecogScores[n1].size(); ++n2)
					{
						size_t mm2 =  n2 / hh;
						size_t mm1 = n2 - (mm2 * hh);
						double dbclass = (double)lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
						fwrite(&dbclass, sizeof(double), 1, pfOut);
						for(size_t n3 = 0; n3 < vvvRecogScores[n1][n2].size(); ++n3)
						{
							double db = vvvRecogScores[n1][n2][n3];
							fwrite(&db, sizeof(double), 1, pfOut);
						}
					}
				}
			}
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

		}
		fclose(pfOut);
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		return 0;

	}

	static int GMMHTKMAPByLAB2_SLIC_BS0_Scores(string& GMMPaths, string& DataLABPaths, string& szRefSCP, string& szOutData)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_MAPAL";

		LoadJPGSLICLabelInfo lb;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGSLICLabel();		

		// Fix JPG Label Info
		FixSLIClabelInfoByRefFileList(lb, szRefSCP);
		
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		FILE* pfOut = NULL;

		pfOut = fopen(szOutData.c_str(), "wb");
		for(long n = 0; n < NN; n++)
		{
			vector< vector< vector<double> > > vvvRecogScores;
			vector<string> vszLoadadFiles;
			int nStartSeq = 0;
			int nStopSeq = 0;
			HTKGMMMAPASeq_SLIC_BS0_Scores(GMMs, hvszKeyNames[vszKeys[n]], vvvRecogScores, vszLoadadFiles, nStartSeq, nStopSeq, GMMPaths);

			for(size_t n1 = 0; n1 < vvvRecogScores.size(); ++n1) // image index
			{
				if(vvvRecogScores[n1].size())
				{
					string szCurrFileName = hvszKeyNames[vszKeys[n]][n1];
					for(size_t n2 = 0; n2 < vvvRecogScores[n1].size(); ++n2) // feature rows
					{
						double dbclass = 0.0;
						if(lb.m_hJPGSLICLabelInfo[szCurrFileName].find(n2) != lb.m_hJPGSLICLabelInfo[szCurrFileName].end())
							dbclass = 1.0;
						fwrite(&dbclass, sizeof(double), 1, pfOut);
						for(size_t n3 = 0; n3 < vvvRecogScores[n1][n2].size(); ++n3) // feature columns, i.e., features at each row
						{
							double db = vvvRecogScores[n1][n2][n3];
							fwrite(&db, sizeof(double), 1, pfOut);
						}
					}
				}
			}
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

		}
		fclose(pfOut);
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		return 0;

	}




	static int GMMHTKMAPByLAB_BS0(string GMMPaths, string DataLABPaths, string szRefSCP, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_SMAPPF";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();

		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);

		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;

			//vector< vector<int> > vvRecogClassesCheat(hvszKeyNames[vszKeys[n]].size());
			//for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			//{
			//	string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
			//	size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
			//	size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
			//	size_t hhww = hh*ww;
			//	vvRecogClassesCheat[nn].reserve(hhww);
			//	vvRecogClassesCheat[nn].resize(hhww);
			//	for(size_t mm = 0; mm < hhww; mm++)
			//	{
			//		size_t mm2 =  mm / hh;
			//		size_t mm1 = mm - (mm2 * hh);
			//		vvRecogClassesCheat[nn][mm] = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
			//	}
			//}


			int nStartSeq = 0;
			int nStopSeq = 0;
			HTKGMMMAPPerFrame_BS0(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq, GMMPaths);

			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = nStartSeq; nn < nStopSeq; nn++) //hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn - nStartSeq].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn - nStartSeq][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn - nStartSeq][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;
	}



	static int GMMHTKMAPByLAB2_BS0(string GMMPaths, string DataLABPaths, string szRefSCP, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_MAPAL";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();

		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);
		
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			int nStartSeq = 0;
			int nStopSeq = 0;
			HTKGMMMAPASeq_BS0(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq, GMMPaths);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;

	}




	static int GMMHTKMAPByLAB(string GMMPaths, string DataLABPaths, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_MAPPF";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;

			HTKGMMAdaptPerFrame(GMMs, GMMPaths, szFEFolder, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;
	//	unordered_map<string, GMM> hGMMs;
	//	vector<GMM> GMMs;
	//	GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
	//	GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
	//	string szDelim = "";
	//	string szCurrFolder = "";

	//	FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);

	//	// Make feature directory
	//	string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
	//	string sysMkDirFE = "mkdir -p " + szFEFolder;
	//#ifdef _MSC_VER
	//	string sysMkDirFE = "mkdir " + szFEFolder;
	//#elif defined(__GNUC__)
	//	string sysMkDirFE = "mkdir -p " + szFEFolder;
	//#endif
	//	system(sysMkDirFE.c_str());

	//	// Temp MLF and SCP (for Adaptation)
	//	string szTmp = DataLABPaths + "_Tmp";
	//	string szTmpMLF = DataLABPaths + "_Tmp.mlf";

	//	//string HERest = "C:\\CollectedDataCopy\\recognizer\\HTK\\bin.win32\\HERest.exe"; //"/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
	//	//string config = "C:\\CollectedDataCopy\\recognizer\\block_config"; //"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
	//	//string words = "C:\\CollectedDataCopy\\recognizer\\block_word"; //"/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
	//	//string GMMCopy = GMMPaths + "_Tmp";
	//	//string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;

	//	// Copy GMM

	//	string GMMCopy = GMMPaths + "_Tmp";
	//#ifdef _MSC_VER
	//	string HERest = "C:\\debug_hand_detection\\tools\\HTK\\bin.win32\\HERest.exe"; //"C:\\CollectedDataCopy\\recognizer\\HTK\\bin.win32\\HERest.exe";
	//	string config = "H:\\Transfer\\AdjustGMMNew\\config";
	//	string words = "H:\\Transfer\\AdjustGMMNew\\words";
	//	string szCopyGMM = "copy " + GMMPaths + " " + GMMCopy;
	//#elif defined(__GNUC__)
	//	string HERest = "/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
	//	string config = "/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
	//	string words = "/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
	//	string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
	//#endif

	//	double Acc = 0.0;
	//	double Acc0 = 0.0;
	//	double N_1 = 0;
	//	int BLKSIZE = 16;
	//	// Now, load label
	//	int nCountX = 0;
	//	int nCountO = 0;

	//	unordered_map<string, int> hszJPGFiles;
	//	unordered_map<string, int> hszJPGFiles1;
	//	unordered_map<string, int> hszJPGFiles0;
	//	//unordered_map<string, int> hszCurrJPGFiles;
	//	//unordered_map<string, JPGLabelInfo> hJPGLavelInfo;
	//	//unordered_map< string, vector<int> > hJPGClasses;
	//	//vector<string> vszKeyOrder;
	//	string szLabelPathName = DataLABPaths;
	//	//unordered_map<string, unordered_map<string, int> > hhTrainFiles;

	//	//int nCountO1 = 0;
	//	//int nCountX1 = 0;
	//	//int nCountO10 = 0;
	//	//int nCountX10 = 0;
	//	//unordered_map<string, JPGLabelInfo> hJPGLavelInfo1;
	//	//unordered_map<string, JPGLabelInfo> hJPGLavelInfo0;
	//	string szLabelPathName1 = DataLABPaths + "_MAPPFHTK";
	//	string szLabelPathName0 = DataLABPaths + "_Orig01";
	//	LoadJPGLabelInfo lb, lb0, lb1;
	//	lb.m_szLabelPathName = szLabelPathName;

	//	lb.LoadJPGInfo();
	//	//GMMAdaptTool::LoadJPGInfo(BLKSIZE, 4, nCountX, nCountO,
	//	//							hszJPGFiles,
	//	//							//hszCurrJPGFiles,
	//	//							hJPGLavelInfo,
	//	//							hJPGClasses, 
	//	//							vszKeyOrder,
	//	//							szLabelPathName,
	//	//							hhTrainFiles);
	//	// Loop through all files included in the label file:
	//	int nTP0 = 0;
	//	int nFP0 = 0;
	//	int nFN0 = 0;
	//	int nTN0 = 0;
	//	int nTP = 0;
	//	int nFP = 0;
	//	int nFN = 0;
	//	int nTN = 0;
	//	for(size_t n = 0; n < lb.m_vszKeyOrder.size(); n++)
	//	{
	//		string szCurrFilePath = lb.m_vszKeyOrder[n];
	//		string szCurrFileName = "";
	//		FindCurrentFileName(szCurrFilePath, szCurrFileName);
	//		// Self adaptation data path name
	//		string szCurrFE1 = szFEFolder + szCurrFileName +"_Hand.fe";
	//		string szCurrFE0 = szFEFolder + szCurrFileName +"_BkGnd.fe";

	//		vector< vector<double> > vvCurrFeatures;
	//		ExtLowLevelImgFeature::ExtColorFeatureRGBOnly(szCurrFilePath, vvCurrFeatures, BLKSIZE);
	//		vector<int> vCurrRecogClasses0;
	//		Test2ClassesGMMs(vvCurrFeatures, GMMs, vCurrRecogClasses0);

	//		JPGLabelInfo* pinfo = &lb.m_hJPGLavelInfo[szCurrFilePath];
	//		//ExtColorFeature(szCurrFilePath, pinfo, vvCurrFeatures, BLKSIZE);
	//		//ExtColorFeatureRGBOnly(szCurrFilePath, pinfo, vvCurrFeatures, BLKSIZE);
	//		int nHeight = pinfo->m_nBlockHeight;
	//		int nWidth = pinfo->m_nBlockWidth;
	//		vector<int> /*vCurrRecogClasses0,*/ vCurrRecogClasses;
	//		double Accn0 = 0.0;
	//		double Accn = 0.0;
	//		vector<int> confres0, confres;
	//		Test2ClassesGMMs(vvCurrFeatures, lb.m_hJPGClasses[szCurrFilePath], GMMs, vCurrRecogClasses0, Accn0, confres0);
	//		nTP0 += confres0[0];
	//		nFP0 += confres0[1];
	//		nFN0 += confres0[2];
	//		nTN0 += confres0[3];

	//		vector< vector<double> > vvRecog0, vvRecog1;
	//		// Collect those recognized as background and hand, respectively
	//		for(int k1 = 0; k1 < nWidth; k1++)
	//		{
	//			for(int k2 = 0; k2 < nHeight; k2++)
	//			{
	//				int nCurrIndex = k1 * nHeight + k2;
	//				if(vCurrRecogClasses0[nCurrIndex] == 0)
	//				{
	//					vvRecog0.push_back(vvCurrFeatures[nCurrIndex]);
	//				}
	//				else if(vCurrRecogClasses0[nCurrIndex] == 1)
	//				{
	//					vvRecog1.push_back(vvCurrFeatures[nCurrIndex]);
	//					lb0.MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1, nWidth, nHeight);
	//					if(lb0.m_hszJPGFiles.find(szCurrFilePath) == lb0.m_hszJPGFiles.end())
	//					{
	//						lb0.m_vszKeyOrder.push_back(szCurrFilePath);
	//					}					
	//					lb0.m_hszJPGFiles[szCurrFilePath] = 1;
	//					//GMMAdaptTool::MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1,
	//										//                                                                                                hJPGLavelInfo0,
	//										//                                                                                                nCountO10,
	//										//                                                                                                nCountX10,
	//										//                                                                                                nWidth,
	//										//                                                                                                nHeight,
	//										//                                                                                                BLKSIZE);
	//				}
	//			}
	//		}
	//		// Prepare mlf and scp for HTK
	//		ofstream fout;
	//		fout.open(szTmp.c_str());
	//		if(!fout.is_open())
	//		{
	//			cerr<<"Cannot Open Temp File so stop it"<<endl;
	//			return -1;
	//		}
	//		ofstream fout02;
	//		fout02.open(szTmpMLF.c_str());
	//		if(!fout02.is_open())
	//		{
	//			cerr<<"Cannot open mlf file for written!"<<endl;
	//			return -1;
	//		}
	//		fout02<<"#!MLF!#"<<endl;
	//		if(0 < vvRecog0.size())
	//		{
	//			WriteFeatureAsHTKMFCC(vvRecog0, szCurrFE0, 9);
	//			fout<<szCurrFE0<<endl;
	//			fout02<<"\"*/"<<szCurrFileName<<"_BkGnd.lab\""<<endl;
	//			fout02<<"BkGnd\n.\n";
	//		}
	//		if(0 < vvRecog1.size())
	//		{
	//			WriteFeatureAsHTKMFCC(vvRecog1, szCurrFE1, 9);
	//			fout<<szCurrFE1<<endl;
	//			fout02<<"\"*/"<<szCurrFileName<<"_Hand.lab\""<<endl;
	//			fout02<<"Hand\n.\n";
	//		}
	//		fout.close();
	//		fout02.close();
	//		// Copy GMM as Temp
	//		system(szCopyGMM.c_str());
	//		string sysCmd = HERest;
	//		sysCmd.append(" -C ");
	//		sysCmd.append(config);
	//		sysCmd.append(" -S ");
	//		sysCmd.append(szTmp);
	//		sysCmd.append(" -I ");
	//		sysCmd.append(szTmpMLF);
	//		sysCmd.append(" -H ");
	//		sysCmd.append(GMMCopy);
	//		sysCmd.append(" -u mpw ");
	//		sysCmd.append(words);					
	//		// Use HERest to adapt the model
	//		system(sysCmd.c_str());

	//		// Load adapted model
	//		unordered_map<string, GMM> hGMMsTmp;
	//		vector<GMM> GMMsTmp;
	//		GMMAdaptTool::LoadGMM(GMMCopy, hGMMsTmp);
	//		GMMAdaptTool::GMMHash2VectBH(hGMMsTmp, GMMsTmp);
	//		//vector<int> confres;
	//		Test2ClassesGMMs(vvCurrFeatures, lb.m_hJPGClasses[szCurrFilePath], GMMsTmp, vCurrRecogClasses, Accn, confres);
	//		nTP += confres[0];
	//		nFP += confres[1];
	//		nFN += confres[2];
	//		nTN += confres[3];
	//		hszJPGFiles1[szCurrFilePath] = 1;
	//		hszJPGFiles0[szCurrFilePath] = 1;
	//		for(int k1 = 0; k1 < nWidth; k1++)
	//		{
	//			for(int k2 = 0; k2 < nHeight; k2++)
	//			{
	//				int nCurrIndex = k1 * nHeight + k2;
	//				if(vCurrRecogClasses[nCurrIndex] == 1)
	//				{
	//					lb1.MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1, nWidth, nHeight);
	//					//GMMAdaptTool::MarkSelectedPointWithoutDrawing(szCurrFilePath, k1, k2, 1, 
	//					//												hJPGLavelInfo1, 	
	//					//												nCountO1,
	//					//												nCountX1,
	//					//												nWidth,
	//					//												nHeight,
	//					//												BLKSIZE);
	//					if(lb1.m_hszJPGFiles.find(szCurrFilePath) == lb1.m_hszJPGFiles.end())
	//					{
	//						lb1.m_vszKeyOrder.push_back(szCurrFilePath);
	//					}					
	//					lb1.m_hszJPGFiles[szCurrFilePath] = 1;
	//				}
	//			}
	//		}
	//		double NN = (double)(nWidth * nHeight);
	//		Acc = (Accn * NN + Acc * N_1)/(N_1+NN);
	//		Acc0 = (Accn0 * NN + Acc0 * N_1)/(N_1+NN);
	//		N_1 += NN;
	//	}
	//	double dbPrecision0 = (double)nTP0/((double)nTP0 + (double)nFP0);
	//	double dbRecall0 = (double)nTP0/((double)nTP0 + (double)nFN0);
	//	double dbAccuracy0 = ((double)nTP0 + (double)nTN0)/((double)nTP0 + (double)nTN0 + (double)nFN0 + (double)nFP0);

	//	double dbPrecision = (double)nTP/((double)nTP + (double)nFP);
	//	double dbRecall = (double)nTP/((double)nTP + (double)nFN);
	//	double dbAccuracy = ((double)nTP + (double)nTN)/((double)nTP + (double)nTN + (double)nFN + (double)nFP);
	//	metrics.resize(14);
	//	metrics[0] = (double) nTP0;
	//	metrics[1] = (double) nFP0;
	//	metrics[2] = (double) nFN0;
	//	metrics[3] = (double) nTN0;
	//	metrics[4] = (double) nTP;
	//	metrics[5] = (double) nFP;
	//	metrics[6] = (double) nFN;
	//	metrics[7] = (double) nTN;
	//	metrics[8] = dbPrecision0;
	//	metrics[9] = dbRecall0;
	//	metrics[10] = dbAccuracy0;
	//	metrics[11] = dbPrecision;
	//	metrics[12] = dbRecall;
	//	metrics[13] = dbAccuracy;

	//	cout<<"Original result is "<<Acc0<<endl;
	//	cout<<"HTK MAP for each frame, result is "<<Acc<<endl;
	//	lb1.SaveJPGInfo(szLabelPathName1);
	//	lb0.SaveJPGInfo(szLabelPathName0);
	//	//GMMAdaptTool::SaveJPGInfo(szLabelPathName1, 
	//	//			hJPGLavelInfo1,
	//	//			vszKeyOrder, 
	//	//			hszJPGFiles1);
	// //       GMMAdaptTool::SaveJPGInfo(szLabelPathName0,
	// //                               hJPGLavelInfo0,
	// //                               vszKeyOrder,
	// //                               hszJPGFiles0);

	//	return 0;
	}

	static int GMMHTKMAPByLAB2(string GMMPaths, string DataLABPaths, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_MAPAL";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;

			HTKGMMAdaptASeq(GMMs, GMMPaths, szFEFolder, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;
	//	unordered_map<string, GMM> hGMMs;
	//	vector<GMM> GMMs;
	//	GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
	//	GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
	//	
	//	int nPos0001 = DataLABPaths.find_last_of("/");
	//	int nPob0001 = DataLABPaths.find_last_of("\\");
	//	int nPos01 = (nPob0001>nPos0001)?nPob0001:nPos0001;
	//	string szCurrFolder = DataLABPaths.substr(0, nPos01);
	//	string szDelim = "";
	//	if(nPos01 == nPob0001) szDelim = "\\";
	//	else if(nPos01 == nPos0001) szDelim = "/";
	//	string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
	//	string sysMkDirFE = "mkdir -p " + szFEFolder;
	//	system(sysMkDirFE.c_str());
	//	string szTmp = DataLABPaths + "_Tmp";
	//	string szTmpMLF = DataLABPaths + "_Tmp.mlf";
	//	string GMMCopy = GMMPaths + "_Tmp";
	//#ifdef _MSC_VER
	//	string HERest = "C:\\debug_hand_detection\\tools\\HTK\\bin.win32\\HERest.exe"; //"C:\\CollectedDataCopy\\recognizer\\HTK\\bin.win32\\HERest.exe";
	//	string config = "H:\\Transfer\\AdjustGMMNew\\config";
	//	string words = "H:\\Transfer\\AdjustGMMNew\\words";
	//	string szCopyGMM = "copy " + GMMPaths + " " + GMMCopy;
	//#elif defined(__GNUC__)
	//	string HERest = "/cluster/users/ycheng31/CollectedDataCopy/recognizer/HTK/bin/HERest";
	//	string config = "/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/config";
	//	string words = "/cluster/users/ycheng31/SkinGMMTools/AdjustGMMNew/words";
	//	string szCopyGMM = "cp " + GMMPaths + " " + GMMCopy;
	//#endif
	//	double Acc = 0.0;
	//	double Acc0 = 0.0;
	//	double N_1 = 0;
	//	int BLKSIZE = 16;
	//	// Now, load label
	//	int nCountX = 0;
	//	int nCountO = 0;
	//	unordered_map<string, int> hszJPGFiles;
	//	unordered_map<string, int> hszJPGFiles1;
	//	unordered_map<string, int> hszJPGFiles0;
	//	//unordered_map<string, int> hszCurrJPGFiles;
	//	//unordered_map<string, JPGLabelInfo> hJPGLavelInfo;
	//	//unordered_map< string, vector<int> > hJPGClasses;
	//	//vector<string> vszKeyOrder;
	//	string szLabelPathName = DataLABPaths;
	//	//unordered_map<string, unordered_map<string, int> > hhTrainFiles;

	//	int nCountO1 = 0;
	//	int nCountX1 = 0;
	//	int nCountO10 = 0;
	//	int nCountX10 = 0;
	//	//unordered_map<string, JPGLabelInfo> hJPGLavelInfo1;
	//	//unordered_map<string, JPGLabelInfo> hJPGLavelInfo0;
	//	string szLabelPathName1 = DataLABPaths + "_MAPHTKAL";
	//	string szLabelPathName0 = DataLABPaths + "_Orig02";
	//	LoadJPGLabelInfo lb, lb0, lb1;
	//	//GMMAdaptTool::LoadJPGInfo(BLKSIZE, 4, nCountX, nCountO,
	//	//							hszJPGFiles,
	//	//							//hszCurrJPGFiles,
	//	//							hJPGLavelInfo,
	//	//							hJPGClasses, 
	//	//							vszKeyOrder,
	//	//							szLabelPathName,
	//	//							hhTrainFiles);
	//	lb.m_szLabelPathName = szLabelPathName;
	//	lb.LoadJPGInfo();
	//	vector<string> vszSeqKeys;
	//	unordered_map<string, vector<string> > hszSeqNames;
	//	int nTP0 = 0;
	//	int nFP0 = 0;
	//	int nFN0 = 0;
	//	int nTN0 = 0;
	//	int nTP = 0;
	//	int nFP = 0;
	//	int nFN = 0;
	//	int nTN = 0;
	//	for(size_t n = 0; n < lb.m_vszKeyOrder.size(); n++)
	//	{
	//		string szCurrFilePath = lb.m_vszKeyOrder[n];
	//		int nPos001s = szCurrFilePath.find_last_of("/");
	//		int nPos001b = szCurrFilePath.find_last_of("\\");
	//		int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
	//		int nPos002 = szCurrFilePath.find_last_of(".");
	//		string szCurrFileName = szCurrFilePath.substr(nPos001 + 1, nPos002 - nPos001 - 1);
	//		//2010_0521_1114_50_000015
	//		int nPos003 = szCurrFileName.find_last_of("_");
	//		string szSeqKey = szCurrFileName.substr(0, nPos003);
	//		if(hszSeqNames.find(szSeqKey) == hszSeqNames.end())
	//		{
	//			vszSeqKeys.push_back(szSeqKey);
	//		}
	//		hszSeqNames[szSeqKey].push_back(szCurrFilePath);
	//	}

	//	// Loop through all files included in the label file:
	//	for(size_t n = 0; n < vszSeqKeys.size(); n++)
	//	{
	//		size_t nNumFramesn = hszSeqNames[vszSeqKeys[n]].size();
	//		// Prepare mlf and scp for HTK
	//		ofstream fout;
	//		fout.open(szTmp.c_str());
	//		if(!fout.is_open())
	//		{
	//			cerr<<"Cannot Open Temp File so stop it"<<endl;
	//			return -1;
	//		}
	//		ofstream fout02;
	//		fout02.open(szTmpMLF.c_str());
	//		if(!fout02.is_open())
	//		{
	//			cerr<<"Cannot open mlf file for written!"<<endl;
	//			return -1;
	//		}
	//		fout02<<"#!MLF!#"<<endl;
	//		vector< vector<double> > vvCurrFeaturesMFrames;
	//		vector< vector< vector<double> > > vvvCurrFeatures;
	//		double Accn0 = 0.0;
	//		double Accn = 0.0;
	//		double N_1m = 0.0;
	//		double Nmm = 0.0;
	//		double N_1m2 = 0.0;
	//		double Nmm2 = 0.0;
	//		for(size_t m = 0; m < nNumFramesn; m++)
	//		{
	//			string szCurrFilePath = hszSeqNames[vszSeqKeys[n]][m];
	//			int nPos001s = szCurrFilePath.find_last_of("/");
	//			int nPos001b = szCurrFilePath.find_last_of("\\");
	//			int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
	//			int nPos002 = szCurrFilePath.find_last_of(".");
	//			string szCurrFileName = szCurrFilePath.substr(nPos001 + 1, nPos002 - nPos001 - 1);
	//			string szCurrFE1 = szFEFolder + szCurrFileName +"_Hand.fe";
	//			string szCurrFE0 = szFEFolder + szCurrFileName +"_BkGnd.fe";
	//			JPGLabelInfo* pinfo = &lb.m_hJPGLavelInfo[szCurrFilePath];
	//			vector< vector<double> > vvCurrFeatures;
	//			ExtColorFeatureRGBOnly(szCurrFilePath, pinfo, vvCurrFeatures, BLKSIZE);
	//			vvvCurrFeatures.push_back(vvCurrFeatures);
	//			for(size_t mm = 0; mm < vvCurrFeatures.size(); mm++)
	//			{
	//				vvCurrFeaturesMFrames.push_back(vvCurrFeatures[mm]);
	//			}
	//			int nHeight = pinfo->m_nBlockHeight;
	//			int nWidth = pinfo->m_nBlockWidth;
	//			vector<int> vCurrRecogClasses0, vCurrRecogClasses;
	//			double Accn0m = 0.0;
	//			vector<int> confres0m;
	//			Test2ClassesGMMs(vvCurrFeatures, lb.m_hJPGClasses[szCurrFilePath], GMMs, vCurrRecogClasses0, Accn0m, confres0m);
	//			nTP0 += confres0m[0];
	//			nFP0 += confres0m[1];
	//			nFN0 += confres0m[2];
	//			nTN0 += confres0m[3];
	//			Nmm = (double)(vCurrRecogClasses0.size());
	//			Accn0 = (N_1m * Accn0 + Accn0m *  Nmm)/(N_1m + Nmm);
	//			N_1m += Nmm;
	//			vector< vector<double> > vvRecog0, vvRecog1;
	//			// Collect those recognized as background and hand, respectively
	//			for(int k1 = 0; k1 < nWidth; k1++)
	//			{
	//				for(int k2 = 0; k2 < nHeight; k2++)
	//				{
	//					int nCurrIndex = k1 * nHeight + k2;
	//					if(vCurrRecogClasses0[nCurrIndex] == 0)
	//					{
	//						vvRecog0.push_back(vvCurrFeatures[nCurrIndex]);
	//					}
	//					else if(vCurrRecogClasses0[nCurrIndex] == 1)
	//					{
	//						vvRecog1.push_back(vvCurrFeatures[nCurrIndex]);
	//						lb0.MarkSelectedPointWithoutDrawing(hszSeqNames[vszSeqKeys[n]][m], k1, k2, 1, nWidth, nHeight);
	//						if(lb0.m_hszJPGFiles.find(szCurrFilePath) == lb0.m_hszJPGFiles.end())
	//						{
	//							lb0.m_vszKeyOrder.push_back(szCurrFilePath);
	//						}					
	//						lb0.m_hszJPGFiles[szCurrFilePath] = 1;
	//						//GMMAdaptTool::MarkSelectedPointWithoutDrawing(hszSeqNames[vszSeqKeys[n]][m], k1, k2, 1, 
	//						//												hJPGLavelInfo0, 	
	//						//												nCountO10,
	//						//												nCountX10,
	//						//												nWidth,
	//						//												nHeight,
	//						//												BLKSIZE);
	//					}
	//				}
	//			}
	//			if(0 < vvRecog0.size())
	//			{
	//				WriteFeatureAsHTKMFCC(vvRecog0, szCurrFE0, 9);
	//				fout<<szCurrFE0<<endl;
	//				fout02<<"\"*/"<<szCurrFileName<<"_BkGnd.lab\""<<endl;
	//				fout02<<"BkGnd\n.\n";
	//			}
	//			if(0 < vvRecog1.size())
	//			{
	//				WriteFeatureAsHTKMFCC(vvRecog1, szCurrFE1, 9);
	//				fout<<szCurrFE1<<endl;
	//				fout02<<"\"*/"<<szCurrFileName<<"_Hand.lab\""<<endl;
	//				fout02<<"Hand\n.\n";
	//			}
	//		}
	//		fout.close();
	//		fout02.close();
	//		// Copy GMM as Temp
	//		system(szCopyGMM.c_str());
	//		string sysCmd = HERest;
	//		sysCmd.append(" -C ");
	//		sysCmd.append(config);
	//		sysCmd.append(" -S ");
	//		sysCmd.append(szTmp);
	//		sysCmd.append(" -I ");
	//		sysCmd.append(szTmpMLF);
	//		sysCmd.append(" -H ");
	//		sysCmd.append(GMMCopy);
	//		sysCmd.append(" -u mpw ");
	//		sysCmd.append(words);					
	//		// Use HERest to adapt the model
	//		system(sysCmd.c_str());

	//		// Load adapted model
	//		unordered_map<string, GMM> hGMMsTmp;
	//		vector<GMM> GMMsTmp;
	//		GMMAdaptTool::LoadGMM(GMMCopy, hGMMsTmp);
	//		GMMAdaptTool::GMMHash2VectBH(hGMMsTmp, GMMsTmp);
	//		double NN = 0.0;
	//		for(size_t m = 0; m < nNumFramesn; m++)
	//		{
	//			vector<int> vCurrRecogClassesm;
	//			double Accnm = 0.0;
	//			vector<int> confres0m;
	//			Test2ClassesGMMs(vvvCurrFeatures[m], lb.m_hJPGClasses[hszSeqNames[vszSeqKeys[n]][m]], GMMsTmp, vCurrRecogClassesm, Accnm, confres0m);
	//			nTP += confres0m[0];
	//			nFP += confres0m[1];
	//			nFN += confres0m[2];
	//			nTN += confres0m[3];
	//			hszJPGFiles1[hszSeqNames[vszSeqKeys[n]][m]] = 1;
	//			hszJPGFiles0[hszSeqNames[vszSeqKeys[n]][m]] = 1;
	//			Nmm2 = (double)vCurrRecogClassesm.size();
	//			Accn = (Nmm2 * Accnm + N_1m2 * Accn)/(Nmm2 + N_1m2);
	//			N_1m2 += Nmm2;
	//			NN += Nmm2;
	//			int nHeight = lb.m_hJPGLavelInfo[hszSeqNames[vszSeqKeys[n]][m]].m_nBlockHeight;
	//			int nWidth = lb.m_hJPGLavelInfo[hszSeqNames[vszSeqKeys[n]][m]].m_nBlockWidth;
	//			string szTempKeyName = hszSeqNames[vszSeqKeys[n]][m];
	//			for(int k1 = 0; k1 < nWidth; k1++)
	//			{
	//				for(int k2 = 0; k2 < nHeight; k2++)
	//				{
	//					int nCurrIndex = k1 * nHeight + k2;
	//					if(vCurrRecogClassesm[nCurrIndex] == 1)
	//					{
	//						lb1.MarkSelectedPointWithoutDrawing(szTempKeyName, k1, k2, 1, nWidth, nHeight);
	//						if(lb1.m_hszJPGFiles.find(szTempKeyName) == lb1.m_hszJPGFiles.end())
	//						{
	//							lb1.m_vszKeyOrder.push_back(szTempKeyName);
	//						}					
	//						lb1.m_hszJPGFiles[szTempKeyName] = 1;
	//						//GMMAdaptTool::MarkSelectedPointWithoutDrawing(hszSeqNames[vszSeqKeys[n]][m], k1, k2, 1, 
	//						//												hJPGLavelInfo1, 	
	//						//												nCountO1,
	//						//												nCountX1,
	//						//												nWidth,
	//						//												nHeight,
	//						//												BLKSIZE);
	//					}
	//				}
	//			}
	//		}
	//		N_1m += Nmm;
	//		Acc = (Accn * NN + Acc * N_1)/(N_1+NN);
	//		Acc0 = (Accn0 * NN + Acc0 * N_1)/(N_1+NN);
	//		N_1 += NN;
	//	}	
	//	double dbPrecision0 = (double)nTP0/((double)nTP0 + (double)nFP0);
	//	double dbRecall0 = (double)nTP0/((double)nTP0 + (double)nFN0);
	//	double dbAccuracy0 = ((double)nTP0 + (double)nTN0)/((double)nTP0 + (double)nTN0 + (double)nFN0 + (double)nFP0);

	//	double dbPrecision = (double)nTP/((double)nTP + (double)nFP);
	//	double dbRecall = (double)nTP/((double)nTP + (double)nFN);
	//	double dbAccuracy = ((double)nTP + (double)nTN)/((double)nTP + (double)nTN + (double)nFN + (double)nFP);
	//	metrics.resize(14);
	//	metrics[0] = (double) nTP0;
	//	metrics[1] = (double) nFP0;
	//	metrics[2] = (double) nFN0;
	//	metrics[3] = (double) nTN0;
	//	metrics[4] = (double) nTP;
	//	metrics[5] = (double) nFP;
	//	metrics[6] = (double) nFN;
	//	metrics[7] = (double) nTN;
	//	metrics[8] = dbPrecision0;
	//	metrics[9] = dbRecall0;
	//	metrics[10] = dbAccuracy0;
	//	metrics[11] = dbPrecision;
	//	metrics[12] = dbRecall;
	//	metrics[13] = dbAccuracy;

	//	cout<<"Original result is "<<Acc0<<endl;
	//	cout<<"HTK MAP for each sequence, result is "<<Acc<<endl;
	//	lb1.SaveJPGInfo(szLabelPathName1);
	//	lb0.SaveJPGInfo(szLabelPathName0);
	//	//GMMAdaptTool::SaveJPGInfo(szLabelPathName1, 
	//	//			hJPGLavelInfo1,
	//	//			vszKeyOrder, 
	//	//			hszJPGFiles1);	
	//	//GMMAdaptTool::SaveJPGInfo(szLabelPathName0, 
	//	//			hJPGLavelInfo0,
	//	//			vszKeyOrder, 
	//	//			hszJPGFiles0);
	//	return 0;
	}




	static int GMMSMAPByLAB(string GMMPaths, string DataLABPaths, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_SMAPPF";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;

			HTKGMMSMAPPerFrame(GMMs, GMMPaths, szFEFolder, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;
	}

	static int GMMSMAPByLAB2(string GMMPaths, string DataLABPaths, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_MAPAL";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;

			HTKGMMSMAPASeq(GMMs, GMMPaths, szFEFolder, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;

	}



	static int GMMSMAPByLAB_BS0(string GMMPaths, string DataLABPaths, string szRefSCP, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_SMAPPF";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();

		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);

		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			int nStartSeq = 0;
			int nStopSeq = 0;
			HTKGMMSMAPPerFrame_BS0(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq, GMMPaths);

			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = nStartSeq; nn < nStopSeq; nn++) //hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn - nStartSeq].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn - nStartSeq][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn - nStartSeq][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;
	}

	static int GMMSMAPByLAB2_BS0(string GMMPaths, string DataLABPaths, string szRefSCP, vector<double>& metrics)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		LoadGMM(GMMPaths, hGMMs);
		GMMHash2VectBH(hGMMs, GMMs);
		int BLOCKSIZE = 16;
		string szCurrFolder, szDelim;
		FindCurrentFolderANDDelim(DataLABPaths, szCurrFolder, szDelim);
		string szFEFolder = szCurrFolder + szDelim + "MyFE" + szDelim;
		
		string szLabelPathName = DataLABPaths;
		string szLabelPathName1 = DataLABPaths + "_MAPAL";
		LoadJPGLabelInfo lb, lb1;
		lb.m_szLabelPathName = szLabelPathName;
		lb.LoadJPGInfo();

		// Fix JPG Label Info
		FixlabelInfoByRefFileList(lb, szRefSCP, BLOCKSIZE);
		
		long N = (long) lb.m_vszKeyOrder.size();
		unordered_map< string, vector<string> > hvszKeyNames;
		vector<string> vszKeys;
		for(long n = 0; n < N; n++)
		{
			string szName, szPath;
			FindCurrentSeqNamePath(lb.m_vszKeyOrder[n], szName, szPath);
			if(hvszKeyNames.find(szName) == hvszKeyNames.end())
			{
				vszKeys.push_back(szName);
			}
			hvszKeyNames[szName].push_back(lb.m_vszKeyOrder[n]);
		}
		long NN = (long) vszKeys.size();
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		for(long n = 0; n < NN; n++)
		{
			vector< vector<int> > vvRecogClasses;
			vector<string> vszLoadadFiles;
			int nStartSeq = 0;
			int nStopSeq = 0;
			HTKGMMSMAPASeq_BS0(GMMs, hvszKeyNames[vszKeys[n]], vvRecogClasses, vszLoadadFiles, BLOCKSIZE, nStartSeq, nStopSeq, GMMPaths);
			cout<<"Sequence "<<n<<" is recognized!"<<endl;

			// Get differences:
			for(size_t nn = 0; nn < hvszKeyNames[vszKeys[n]].size(); nn++)
			{
				string szCurrFileName = hvszKeyNames[vszKeys[n]][nn];
				for(size_t mm = 0; mm < vvRecogClasses[nn].size(); mm++)
				{
					size_t hh = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockHeight;
					size_t ww = (size_t)lb.m_hJPGLavelInfo[szCurrFileName].m_nBlockWidth;
					size_t mm2 =  mm / hh;
					size_t mm1 = mm - (mm2 * hh);
					int nclass = lb.m_hJPGLavelInfo[szCurrFileName].m_ppnLabel[mm1][mm2];
					if(nclass == vvRecogClasses[nn][mm])
					{
						if(nclass == 1)
							TP++;
						else
							TN++;
					}
					else
					{
						if(nclass == 0)
							FP++;
						else
							FN++;
					}
					if(1 == vvRecogClasses[nn][mm])
					{
						lb1.MarkSelectedPointWithoutDrawing(szCurrFileName, mm2, mm1, 1, ww, hh);
					}
				}
			}
		}
		//double Acc = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics.resize(7);
		double dbPrecision = (double)(TP) / (double)(TP + FP);
		double dbRecall = (double)(TP) / (double)(TP + FN);
		double dbAccuracy = (double)(TP+TN)/(double)(TP+TN+FP+FN);
		metrics[0] = (double) TP;
		metrics[1] = (double) FP;
		metrics[2] = (double) FN;
		metrics[3] = (double) TN;
		metrics[4] = dbPrecision;
		metrics[5] = dbRecall;
		metrics[6] = dbAccuracy;
		lb1.SaveJPGInfo(szLabelPathName1);
		return 0;

	}


	static int Test2ClassesGMMs(string GMMPaths, string DataPaths, double& Acc, int nFrameSize)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		vector< vector<double> > vvData;
		vector<int> vClasses, recogClasses;
		GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
		GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
		GMMAdaptTool::LoadLogData(DataPaths, vvData, vClasses);
		vector<int> confres;
		int nretval = Test2ClassesGMMs(vvData, vClasses, GMMs, recogClasses, Acc, confres);
		return nretval;
	}

	static int Test2ClassesGMMsImg(string GMMPaths, string DataPaths, double& Acc, int nFrameSize)
	{
		unordered_map<string, GMM> hGMMs;
		vector<GMM> GMMs;
		vector< vector<double> > vvData;
		vector<int> vClasses, recogClasses;
		GMMAdaptTool::LoadGMM(GMMPaths, hGMMs);
		GMMAdaptTool::GMMHash2VectBH(hGMMs, GMMs);
		GMMAdaptTool::LoadLogData(DataPaths, vvData, vClasses);
		vector<int> confres;
		int nretval = Test2ClassesGMMs(vvData, vClasses, GMMs, recogClasses, Acc, confres);
		return nretval;
	}

	static int TestBoosted2ClassesGMMs(vector<string> vGMMPaths, string DataPaths, double& Acc, vector<double>& vAlpha)
	{
		unordered_map<string, GMM> hGMMs;
		//vector<GMM> GMMs;
		vector< vector<double> > vvData;
		vector<int> vClasses, recogClasses;
		size_t T = vGMMPaths.size();
		if(T != vAlpha.size())
		{
			cerr<<"number of models is not equal to number of weights"<<endl;
			return -1;
		}
		vector< vector<GMM> > vGMMs(T);
		GMMAdaptTool::LoadLogData(DataPaths, vvData, vClasses);
		for(size_t t = 0; t < T; t++)
		{
			unordered_map<string, GMM> hTmpGMM;
			GMMAdaptTool::LoadGMM(vGMMPaths[t], hTmpGMM);
			GMMAdaptTool::GMMHash2VectBH(hTmpGMM, vGMMs[t]);
		}	
		int nretval = TestBoosted2ClassesGMMs(vvData, vClasses, vGMMs, vAlpha, recogClasses, Acc);
		
		return nretval;
	}

	static int TestBoosted2ClassesGMMs(vector< vector<double> >& vvData, vector<int>& vClasses, vector< vector<GMM> >& vGMMs, vector<double>& vAlpha, vector<int>& recogClasses, double& Acc)
	{
		int N = (int)vvData.size();
		recogClasses.resize(N, 0);
		int I = (int)vGMMs[0].size();
		for(size_t t = 1; t < vGMMs.size(); t++)
		{
			if(I != vGMMs[t].size())
			{
				cerr<<"Not all GMMs has the same number of classes"<<endl;
				return -1;
			}
		}
		int ErrCount = 0;
		for(int n = 0; n < N; n++)
		{
			//double maxni = GMMs[0].CompLikelihoodPart(vvData[n], vDimConsidered);//GMMs[0].CompLikelihood(vvData[n]);
			double maxni = 0.0;
			int argmaxni = 0;
			for(size_t t = 0; t < vGMMs.size(); t++)
			{
				maxni += vAlpha[t] * vGMMs[t][0].CompLikelihood(vvData[n]);
			}
			for(int i = 1; i < I ;i++)
			{
				//double dbi = GMMs[i].CompLikelihoodPart(vvData[n], vDimConsidered);//GMMs[i].CompLikelihood(vvData[n]);
				double dbi = 0.0;
				
				for(size_t t = 0; t < vGMMs.size(); t++)
				{			
					dbi += vAlpha[t] * vGMMs[t][i].CompLikelihood(vvData[n]);
				}
				if(dbi > maxni)
				{
					maxni = dbi;
					argmaxni = i;
				}
			}
			recogClasses[n] = argmaxni;
		}
		for(int n = 0; n < N; n++)
		{
			if(recogClasses[n] != vClasses[n])
			{
				ErrCount++;
			}	
		}
		Acc = (double)(N - ErrCount) / (double) N;
		return 1;
	}

	static void CompMean(vector< vector<double> >& vvData, vector<double>& vMean)
	{
		int N = (int)vvData.size();
		double dbN = (double)N;
		if(N <= 0)
		{
			cerr<<"Null Data"<<endl;
			return;
		}
		int D = (int)vvData[0].size();
		vMean = vector<double>(D, 0.0);
		
		for(int n = 0; n < N; n++)
		{
			for(int d = 0; d < D; d++)
			{
				vMean[d] += vvData[n][d];
			}
		}
		for(int d = 0; d < D; d++)
		{
			vMean[d] /= dbN;
		}
	}

	static void CompDiagVar(vector< vector<double> >& vvData, matrix<double>& diagVar)
	{
		int N = (int)vvData.size();
		double dbN = (double)N;
		if(N <= 0)
		{
			cerr<<"Null Data"<<endl;
			return;
		}
		int D = (int)vvData[0].size();
		diagVar.SetSize(D, D);
		diagVar.Null();
		vector<double> vMean;
		CompMean(vvData, vMean);
		for(int n = 0; n < N; n++)
		{
			for(int d = 0; d < D; d++)
			{
				diagVar(d, d) += pow(vvData[n][d] - vMean[d], 2.0);
			}
		}	
		for(int d = 0; d < D; d++)
		{
			diagVar(d, d) /= dbN;
		}

	}
	static void GetLargestConnectedComponents(vector<int>& orig, vector<int>&refined)
	{
		refined = vector<int>(orig.size(), 0);
		unsigned char* img = new unsigned char[WIDTH * HEIGHT];
		memset(img, 0, WIDTH * HEIGHT);
		for(int nn1 = 0; nn1 < WIDTH; nn1++)
		{
			for(int nn2 = 0; nn2 < HEIGHT; nn2++)
			{
				int nn12 = nn1 * HEIGHT + nn2;
				int nn122 = nn2 * WIDTH + nn1;
				if(orig[nn12] == 1)
				{
					img[nn122] = 1;
				}
			}
		}
		unsigned char *out_uc = new unsigned char[WIDTH * HEIGHT];
		memset(out_uc, 0, WIDTH * HEIGHT);
		ConnectedComponents cc(WIDTH * HEIGHT);
		cc.connected(img, out_uc, WIDTH, HEIGHT,
			std::equal_to<unsigned char>(),
			false);	
		cout<<endl;
		
		unordered_map<unsigned char, int> hCount;
		//hash_map<unsigned char, int> hCount;

		for(int r=0; r<HEIGHT; ++r) {
			for(int c=0; c<WIDTH; ++c)
			{
				unsigned char val = out_uc[r*WIDTH+c];
				if(val != 0 && img[r*WIDTH+c] != 0)
				{
					if(hCount.find(out_uc[r*WIDTH+c]) ==hCount.end())
					{
						hCount[out_uc[r*WIDTH+c]] = 1;
					}
					else
					{
						hCount[out_uc[r*WIDTH+c]]++;
					}
				}
			}
		}		
		unordered_map<unsigned char, int>::iterator itr;
		//hash_map<unsigned char, int>::iterator itr;
		int nmaxcount = 0;
		unsigned char nkey = 0;
		for(itr = hCount.begin();itr != hCount.end(); itr++)
		{
			if(itr->second > nmaxcount)
			{
				nmaxcount = itr->second;
				nkey = itr->first;
			}
		}
		for(int r=0; r<HEIGHT; ++r) 
		{
			for(int c=0; c<WIDTH; ++c)
			{
				if(out_uc[r*WIDTH+c] == nkey)
				{
					refined[c * HEIGHT + r] = 1;
				}
				else
				{
					refined[c * HEIGHT + r] = 0;
				}
			}
		}
		delete [] img;
		delete [] out_uc;
	}
		
	static void ExtColorFeature(string szFileName, JPGLabelInfo* pinfo, vector< vector<double> >& vvFeatures, int BLKSIZE)
	{
		ExtLowLevelImgFeature::ExtColorFeatures(szFileName, vvFeatures, BLKSIZE);
	//	IplImage* pimg;
	//#ifdef LPF
	//	CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
	//	cvSetReal2D(hlpf, 0, 0, 0.1);
	//	cvSetReal2D(hlpf, 0, 1, 0.1);
	//	cvSetReal2D(hlpf, 0, 2, 0.1);
	//	cvSetReal2D(hlpf, 1, 0, 0.1);
	//	cvSetReal2D(hlpf, 1, 1, 0.8);
	//	cvSetReal2D(hlpf, 1, 2, 0.1);
	//	cvSetReal2D(hlpf, 2, 0, 0.1);
	//	cvSetReal2D(hlpf, 2, 1, 0.1);
	//	cvSetReal2D(hlpf, 2, 2, 0.1);
	//	IplImage* pimg0 = cvLoadImage(szFileName.c_str());
	//	IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
	//	IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
	//	pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
	//	cvConvert(pimg0, pimg1);
	//	cvFilter2D(pimg1, pimg2, hlpf);
	//	cvConvert(pimg2, pimg);
	//	cvReleaseImage(&pimg0);
	//	cvReleaseImage(&pimg1);
	//	cvReleaseImage(&pimg2);
	//	cvReleaseMat(&hlpf);
	//#else
	//	pimg = cvLoadImage(szFileName.c_str());
	//#endif
	//	cv::Mat img(pimg);
	//	cv::Mat fimg;
	//	img.convertTo(fimg, CV_32FC3);

	//	for(int x = 0; x < img.cols / BLKSIZE; x++)
	//	{
	//		cv::Scalar fmeanxy, stddevxy;
	//		cv::Scalar fmeanxyhsv, stddevxyhsv;
	//		cv::Scalar fmeanxyx, stddevxyx;
	//		cv::Scalar fmeanxyy, stddevxyy;
	//		cv::Mat colx = fimg/*fres*/.colRange(x*16, (x+1)*16/*min((x+1)*8, img.cols - 1)*/); 
	//		for(int y = 0; y < img.rows / BLKSIZE; y++)
	//		{	
	//			cv::Mat sub = colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
	//			cv::Mat sub8bits;
	//			sub.convertTo(sub8bits, CV_8UC3);
	//			cv::Mat subhsv;
	//			cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
	//			cv::meanStdDev(sub, fmeanxy, stddevxy);
	//			cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);
	//			//JPGLabelInfo* pinfo = &m_hJPGLavelInfo[m_vszKeyOrder[n]];
	//			//int nLabel = pinfo->m_ppnLabel[y][x];
	//			double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
	//			double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
	//			double dbBMean = fmeanxy.val[0];
	//			double dbGMean = fmeanxy.val[1];
	//			double dbRMean = fmeanxy.val[2];
	//			double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
	//			double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
	//			double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
	//			double dbVMean = fmeanxyhsv.val[2];
	//			double dbSMean = fmeanxyhsv.val[1];
	//			double dbHMean = fmeanxyhsv.val[0];
	//			vector<double> vdbTemp;
	//			//vdbTemp.push_back(dbAvgStdev);
	//			vdbTemp.push_back(dbIMean);
	//			vdbTemp.push_back(dbBMean);
	//			vdbTemp.push_back(dbGMean);
	//			vdbTemp.push_back(dbRMean);
	//			vdbTemp.push_back(dbExBMean);
	//			vdbTemp.push_back(dbExGMean);
	//			vdbTemp.push_back(dbExRMean);
	//			//vdbTemp.push_back(dbVMean);
	//			//vdbTemp.push_back(dbSMean);
	//			//vdbTemp.push_back(dbHMean);
	//			vvFeatures.push_back(vdbTemp);
	//		}
	//	}
	//	img.release();
	//	cvReleaseImage(&pimg);
	}

	static void ExtColorFeatureRGBOnly(string szFileName, JPGLabelInfo* pinfo, vector< vector<double> >& vvFeaturesSkin, vector< vector<double> >& vvFeaturesNonSkin)
	{
		if(pinfo) {
			vector< vector<double> > vvFeatures;
			ExtLowLevelImgFeature::ExtColorFeatureRGBOnly(szFileName, vvFeatures, pinfo->BLOCKSIZE);
			int nidx = 0;
			for(int x = 0; x < pinfo->m_nBlockWidth; x++)
			{
				for(int y = 0; y < pinfo->m_nBlockHeight; y++)
				{
					if(0 == pinfo->m_ppnLabel[y][x])
					{
						vvFeaturesNonSkin.push_back(vvFeatures[nidx]);
					}
					else if(1 == pinfo->m_ppnLabel[y][x])
					{
						vvFeaturesSkin.push_back(vvFeatures[nidx]);
					}
					nidx++;
				}
			}
		}
		//pinfo->
	//	IplImage* pimg;
	//#ifdef LPF
	//	CvMat* hlpf = cvCreateMat(3, 3, CV_32FC1);
	//	cvSetReal2D(hlpf, 0, 0, 0.1);
	//	cvSetReal2D(hlpf, 0, 1, 0.1);
	//	cvSetReal2D(hlpf, 0, 2, 0.1);
	//	cvSetReal2D(hlpf, 1, 0, 0.1);
	//	cvSetReal2D(hlpf, 1, 1, 0.8);
	//	cvSetReal2D(hlpf, 1, 2, 0.1);
	//	cvSetReal2D(hlpf, 2, 0, 0.1);
	//	cvSetReal2D(hlpf, 2, 1, 0.1);
	//	cvSetReal2D(hlpf, 2, 2, 0.1);
	//	IplImage* pimg0 = cvLoadImage(szFileName.c_str());
	//	IplImage* pimg1 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
	//	IplImage* pimg2 = cvCreateImage(cvSize(pimg0->width,pimg0->height), IPL_DEPTH_32F, pimg0->nChannels);
	//	pimg = cvCreateImage(cvSize(pimg0->width,pimg0->height), pimg0->depth, pimg0->nChannels);
	//	cvConvert(pimg0, pimg1);
	//	cvFilter2D(pimg1, pimg2, hlpf);
	//	cvConvert(pimg2, pimg);
	//	cvReleaseImage(&pimg0);
	//	cvReleaseImage(&pimg1);
	//	cvReleaseImage(&pimg2);
	//	cvReleaseMat(&hlpf);
	//#else
	//	pimg = cvLoadImage(szFileName.c_str());
	//#endif
	//	cv::Mat img(pimg);
	//	cv::Mat fimg;
	//	img.convertTo(fimg, CV_32FC3);

	//	for(int x = 0; x < img.cols / BLKSIZE; x++)
	//	{
	//		cv::Scalar fmeanxy, stddevxy;
	//		cv::Scalar fmeanxyhsv, stddevxyhsv;
	//		cv::Scalar fmeanxyx, stddevxyx;
	//		cv::Scalar fmeanxyy, stddevxyy;
	//		cv::Mat colx = fimg/*fres*/.colRange(x*16, (x+1)*16/*min((x+1)*8, img.cols - 1)*/); 
	//		for(int y = 0; y < img.rows / BLKSIZE; y++)
	//		{	
	//			cv::Mat sub = colx.rowRange(y*BLKSIZE, (y+1)*BLKSIZE/*min((y+1)*8, img.rows - 1)*/);
	//			cv::Mat sub8bits;
	//			sub.convertTo(sub8bits, CV_8UC3);
	//			cv::Mat subhsv;
	//			cv::cvtColor(sub8bits, subhsv, CV_BGR2HSV);
	//			cv::meanStdDev(sub, fmeanxy, stddevxy);
	//			cv::meanStdDev(subhsv, fmeanxyhsv, stddevxyhsv);
	//			//JPGLabelInfo* pinfo = &m_hJPGLavelInfo[m_vszKeyOrder[n]];
	//			//int nLabel = pinfo->m_ppnLabel[y][x];
	//			//double dbAvgStdev = (stddevxy.val[0] + stddevxy.val[1] + stddevxy.val[2])/3.0;
	//			//double dbIMean = (fmeanxy.val[0] + fmeanxy.val[1] + fmeanxy.val[2]) / 3.0;
	//			double dbBMean = fmeanxy.val[0];
	//			double dbGMean = fmeanxy.val[1];
	//			double dbRMean = fmeanxy.val[2];
	//			//double dbExBMean = 2 * fmeanxy.val[0] - (fmeanxy.val[1] + fmeanxy.val[2]);
	//			//double dbExGMean = 2 * fmeanxy.val[1] - (fmeanxy.val[0] + fmeanxy.val[2]);
	//			//double dbExRMean = 2 * fmeanxy.val[2] - (fmeanxy.val[0] + fmeanxy.val[1]);
	//			//double dbVMean = fmeanxyhsv.val[2];
	//			//double dbSMean = fmeanxyhsv.val[1];
	//			//double dbHMean = fmeanxyhsv.val[0];
	//			vector<double> vdbTemp;
	//			//vdbTemp.push_back(dbAvgStdev);
	//			//vdbTemp.push_back(dbIMean);
	//			vdbTemp.push_back(dbBMean);
	//			vdbTemp.push_back(dbGMean);
	//			vdbTemp.push_back(dbRMean);
	//			//vdbTemp.push_back(dbExBMean);
	//			//vdbTemp.push_back(dbExGMean);
	//			//vdbTemp.push_back(dbExRMean);
	//			//vdbTemp.push_back(dbVMean);
	//			//vdbTemp.push_back(dbSMean);
	//			//vdbTemp.push_back(dbHMean);
	//			vvFeatures.push_back(vdbTemp);
	//		}
	//	}
	//	img.release();
	//	cvReleaseImage(&pimg);
	}

	// SLIC version of the above
	static void ExtColorFeatureRGBOnly(string szFileName, unordered_map<int, int>& slic1, vector< vector<double> >& vvFeaturesSkin, vector< vector<double> >& vvFeaturesNonSkin)
	{
		vector< vector<double> > vvFeatures;
		string szSLICFileName;
		GMMAdaptTool::GuessSLICPath(szFileName, szSLICFileName);
			
		ExtLowLevelImgFeature::ExtColorFeatureSLICRGBOnly_LoadLabel(szFileName, vvFeatures, szSLICFileName);
			
		if(slic1.size() > 0) {
			for(size_t st = 0; st < vvFeatures.size(); ++st)
			{
				int id = (int)st;
				if(slic1.find(id) != slic1.end())
				{
					vvFeaturesSkin.push_back(vvFeatures[st]);
				}
				else
				{
					vvFeaturesNonSkin.push_back(vvFeatures[st]);
				}
			}
		}
		else {
			for(size_t st = 0; st < vvFeatures.size(); ++st)
			{
				vvFeaturesNonSkin.push_back(vvFeatures[st]);
			}
		}
	}
	//

	// Get the minimum and maximum value of matrix stored in vector< vector<double> > form at n-th dimension
	static HRESULT GetvvRangeAtnD(const vector< vector<double> >& vvData, unsigned int n, double& dbMin, double& dbMax)
	{
		size_t stvvSize = vvData.size();
		if(stvvSize == 0)
		{
			return ERROR_EMPTY;
		}
		if(vvData[0].size() < n + 1)
		{
			return ERROR_INVALID_INDEX;
		}
		double l_dbMin = vvData[0][n];
		double l_dbMax = vvData[0][n];
		for(size_t t = 1; t < stvvSize; t++)
		{
			if(vvData[t].size() < n + 1)
			{
				return ERROR_INVALID_INDEX;
			}
			l_dbMin = (vvData[t][n] < l_dbMin)?vvData[t][n]:l_dbMin;
			l_dbMax = (vvData[t][n] > l_dbMax)?vvData[t][n]:l_dbMax;
		}
		dbMin = l_dbMin;
		dbMax = l_dbMax;
		return S_OK;
	}
	
	// normalize so that the sequence are in [0,1]
	static HRESULT Makein01(const vector< vector<double> >& vvData, vector< vector<double> >& vvNewData)
	{
		HRESULT hr = S_OK;
		vector<double> vMean;
		size_t nSize = vvData.size();
		size_t nDim = vvData[0].size();
		vvNewData.resize(nSize);
		if(nDim >= 2)
		{
			for(unsigned int n = 0; n < nDim; n++)
			{
				double dbMin = 0.0;
				double dbMax = 0.0;
				if(S_OK == (hr = GetvvRangeAtnD(vvData, n, dbMin, dbMax)) )
				{
					double dbhalfRange = dbMax - dbMin;

					for(unsigned int t = 0; t < nSize; t++)
					{
						if(vvNewData[t].size() != nDim)
						{
							vvNewData[t].resize(nDim);
						}
						if(dbhalfRange > 0.0)
						{
							vvNewData[t][n] = (vvData[t][n] - dbMin)/dbhalfRange;
						}
						else
						{
							vvNewData[t][n] = 0.0;
							hr = ERROR_INVALID_DATA;
						}
						//if(n < 2)
						//{
						//	vvNewData[t][n] = (vvData[t][n] - dbMin)/dbhalfRange;
						//}
						//else
						//{
						//	vvNewData[t][n] = vvData[t][n];
						//}
					}
				}
				else
				{
					return hr;
				}
			}
			return hr;
		}
		return S_FALSE;
	}

	// Make vector< vector<double> > data zero mean and range from [-1, 1]
	static HRESULT MakeZeroMeanAndUnitRange(const vector< vector<double> >& vvData, vector< vector<double> >& vvNewData)
	{
		HRESULT hr = S_OK;
		vector<double> vMean;
		size_t nSize = vvData.size();
		if(nSize == 0)
		{
			return ERROR_INVALID_DATA;
		}
		size_t nDim = vvData[0].size();
		if(nDim == 0)
		{
			return ERROR_INVALID_DATA;
		}
		if(S_OK == (hr = CalcvvMean(vvData, vMean)) )
		{
			vvNewData.resize(nSize);
			for(size_t n = 0; n < nDim; n++)
			{
				double dbMin = 0.0;
				double dbMax = 0.0;
				if(S_OK == (hr = GetvvRangeAtnD(vvData, n, dbMin, dbMax)) )
				{
					double dbhalfRange = max(fabs(dbMax - vMean[n]), fabs(vMean[n] - dbMin));
					for(size_t t = 0; t < nSize; t++)
					{
						if(vvNewData[t].size() != nDim)
						{
							vvNewData[t].resize(nDim);
						}
						vvNewData[t][n] = (vvData[t][n] - vMean[n])/dbhalfRange;
					}
				}
				else
				{
					return hr;
				}
			}
		}
		return hr;
	}

	// Calculate the mean vector of a vector< vector<double> > data matrix by trating each vector<double> as a data point
	static HRESULT CalcvvMean(const vector< vector<double> >& vvData, vector<double>& vMean)
	{
		size_t stvvSize = vvData.size();
		if(stvvSize == 0)
		{
			return ERROR_EMPTY;
		}
		size_t nDim = vvData[0].size();
		if(nDim == 0)
		{
			return ERROR_INVALID_DATA;
		}
		double dbSize = (double)stvvSize;
		vMean.resize(nDim, 0.0);
		for(size_t t = 1; t < stvvSize; t++)
		{
			if(vvData[t].size() != nDim)
			{
				return ERROR_INVALID_DATA;
			}
			for(size_t n = 0; n < nDim; n++)
			{
				vMean[n] += vvData[t][n]/dbSize;
			}
		}
 
		return S_OK;
	}

	static long ReadHTKMFCCFeature(vector< vector<double> >& vvData, string szInFileName )
	{
		FILE* fin = NULL;
		long nSamples = 0;
		long sampPeriod = 0;
		short sampSize = 0;
		short parmKind = 0;
		fin = fopen(szInFileName.c_str(), "rb");
		if(fin == NULL)
		{
			return 2L;
		}
		fread(&nSamples, sizeof(long), 1, fin);
		fread(&sampPeriod, sizeof(long), 1, fin);
		fread(&sampSize, sizeof(short), 1, fin);
		fread(&parmKind, sizeof(short), 1, fin);
		long nSize = sampSize / sizeof(float);
		float fTmp = 0.0;
		vvData.resize(nSamples);
		for(long n = 0; n < nSamples; n++)
		{
			for(long m = 0; m < nSize; m++)
			{
				fread(&fTmp, sizeof(float), 1, fin);
				vvData[n].push_back((double)fTmp);
			}
		}
		fclose(fin);
		return 0;
	}

	//static long WriteFeatureAsHTKMFCC(const vector< vector<double> >& vvData, string szOutFileName, short sampKind)
	//{
	//	vector< vector<double> > vvNewData;
	//	size_t nSize = vvData.size();
	//	//for(size_t t = 0; t < vvNewData.size(); t++)
	//	//{
	//	//	for(size_t d = 0; d < vvNewData[t].size(); d++)
	//	//	{
	//	//	}
	//	//}
	//	short sSize = (short)vvData[0].size();
	//	int nSamples = (int)vvData.size();
	//	int sampPeriod = 666666;//10000000/15;
	//	short sampSize = sSize * sizeof(float);
	//	short parmKind = sampKind; //9; //2054; // MFCC_Z//9; // USER //2054; // MFCC_Z

	//	FILE* hFout = 0;
	//	hFout = fopen(szOutFileName.c_str(), "wb");
	//	if(hFout)
	//	{
	//		fwrite(&nSamples, sizeof(int), 1, hFout);
	//		fwrite(&sampPeriod, sizeof(int), 1, hFout);
	//		fwrite(&sampSize, sizeof(short), 1, hFout);
	//		fwrite(&parmKind, sizeof(short), 1, hFout);
	//		for(unsigned int n = 0; n < nSize; n++)
	//		{
	//			float fTmp = 0.0;
	//			for(short s = 0; s < sSize; s++)
	//			{
	//				fTmp = (float)vvData[n][s];
	//				fwrite(&fTmp, sizeof(float), 1, hFout);
	//			}
	//		}
	//	}
	//	else
	//	{
	//		return ERROR_OPEN_FAILED;
	//	}
	//	if(hFout)
	//	{
	//		fclose(hFout);
	//	}
	//	return 0;
	//}

	// Write the feature vector as HTK format MFCC feature file
	static HRESULT WriteFeatureAsHTKMFCC(const vector< vector<double> >& vvData, string szOutFileName, short parmKind=2054, int sampPeriod = 666666)
	{
		vector< vector<double> > vvNewData;
		size_t nSize = vvData.size();
		//cout<<"before makein01\n";
		HRESULT hr = Makein01(vvData, vvNewData);
		if(hr == S_OK)
		{
			//cout<<"ok, start writing htk feature\n";
			short sSize = (short)vvData[0].size();
			int nSamples = (int)vvData.size();
			//int sampPeriod = 666666;//10000000/15;
			short sampSize = sSize * sizeof(float);
			//short parmKind = 2054; // MFCC_Z//9; // USER //2054; // MFCC_Z

			FILE* hFout = 0;
			hFout = fopen(szOutFileName.c_str(), "wb");
			if(hFout)
			{
				fwrite(&nSamples, sizeof(int), 1, hFout);
				fwrite(&sampPeriod, sizeof(int), 1, hFout);
				fwrite(&sampSize, sizeof(short), 1, hFout);
				fwrite(&parmKind, sizeof(short), 1, hFout);
				for(unsigned int n = 0; n < nSize; n++)
				{
					float fTmp = 0.0;
					for(short s = 0; s < sSize; s++)
					{
						//cout<<"before writing "<<n<<" th feature at "<<s<<" th dimension\n";
						fTmp = (float)vvNewData[n][s];
						fwrite(&fTmp, sizeof(float), 1, hFout);
						//cout<<"after writing "<<n<<" th feature at "<<s<<" th dimension\n";

					}
				}
			}
			else
			{
				return ERROR_OPEN_FAILED;
			}
			if(hFout)
			{
				fclose(hFout);
			}
		}
		return hr;
	}

	// Write the feature vector as HTK format MFCC feature file and its lower and upper halves

	static HRESULT WriteFeatureAsHTKMFCCAndItsHalves(const vector< vector<double> >& vvData, string szOutFileName, string szOutFileName_uhalf, string szOutFileName_lhalf)
	{
		vector< vector<double> > vvNewData;
		size_t nSize = vvData.size();
		HRESULT hr = Makein01(vvData, vvNewData);
		if(hr == S_OK)
		{
			short sSize = (short)vvData[0].size();
			int nSamples = (int)vvData.size();
			int sampPeriod = 666666;//10000000/15;
			short sampSize = sSize * sizeof(float);
			short sampSizeU = sSize/2 * sizeof(float);
			short sampSizeL = (sSize - sSize/2) * sizeof(float);
			short parmKind = 2054; // MFCC_Z  //9; // USER //2054; // MFCC_Z
			
			FILE* hFout = 0;
			FILE* hFoutU = 0;
			FILE* hFoutL = 0;
			hFout = fopen(szOutFileName.c_str(), "wb");
			hFoutU = fopen(szOutFileName_uhalf.c_str(), "wb");
			hFoutL = fopen(szOutFileName_lhalf.c_str(), "wb");
			if(hFout)
			{
				fwrite(&nSamples, sizeof(int), 1, hFout);
				fwrite(&sampPeriod, sizeof(int), 1, hFout);
				fwrite(&sampSize, sizeof(short), 1, hFout);
				fwrite(&parmKind, sizeof(short), 1, hFout);
				for(unsigned int n = 0; n < nSize; n++)
				{
					float fTmp = 0.0;
					for(short s = 0; s < sSize; s++)
					{
						fTmp = (float)vvNewData[n][s];
						fwrite(&fTmp, sizeof(float), 1, hFout);
					}
				}
				fclose(hFout);
			}
			else 
			{
				hr = ERROR_OPEN_FAILED;
			}
			if(hFoutU)
			{
				fwrite(&nSamples, sizeof(int), 1, hFoutU);
				fwrite(&sampPeriod, sizeof(int), 1, hFoutU);
				fwrite(&sampSizeU, sizeof(short), 1, hFoutU);
				fwrite(&parmKind, sizeof(short), 1, hFoutU);
				for(unsigned int n = 0; n < nSize; n++)
				{
					float fTmp = 0.0;
					for(short s = 0; s < sSize/2; s++)
					{
						fTmp = (float)vvNewData[n][s];
						fwrite(&fTmp, sizeof(float), 1, hFoutU);
					}
				}
				fclose(hFoutU);
			}
			else 
			{
				hr = ERROR_OPEN_FAILED;
			}
			if(hFoutL)
			{
				fwrite(&nSamples, sizeof(int), 1, hFoutL);
				fwrite(&sampPeriod, sizeof(int), 1, hFoutL);
				fwrite(&sampSizeL, sizeof(short), 1, hFoutL);
				fwrite(&parmKind, sizeof(short), 1, hFoutL);
				for(unsigned int n = 0; n < nSize; n++)
				{
					float fTmp = 0.0;
					for(short s = sSize/2; s < sSize; s++)
					{
						fTmp = (float)vvNewData[n][s];
						fwrite(&fTmp, sizeof(float), 1, hFoutL);
					}
				}
				fclose(hFoutL);
			}
			else 
			{
				hr = ERROR_OPEN_FAILED;
			}
		}
		else if(hr == ERROR_INVALID_DATA)
		{
		}
		return hr;
	}

	static HRESULT WriteFeatureAsDEBUGTRAJ(const vector< vector<double> >& vvData, string szOutFileName1)
	{
		vector< vector<double> > vvNewData;
		size_t nSize = vvData.size();
		HRESULT hr = Makein01(vvData, vvNewData);
		if(hr == S_OK)
		{
			IplImage* pImg1 = cvCreateImage(cvSize(500, 500), IPL_DEPTH_8U, 3);
			cvSet(pImg1, CV_RGB(255, 255, 255));
			vector<int> vPrevData(2, 0);
			for(unsigned int n = 0; n < nSize; n++)
			{
				double dx1 = floor(vvNewData[n][0] * 500.0 + 0.5);
				double dy1 = floor(vvNewData[n][1] * 500.0 + 0.5);
				int nx1 = (int)dx1;
				int ny1 = (int)dy1;
				if(n == 0)
				{
					cvCircle(pImg1, cvPoint(nx1, ny1), 3, CV_RGB(0, 255, 0), 2);
				}
				else if(n == nSize - 1)
				{
					cvCircle(pImg1, cvPoint(nx1, ny1), 3, CV_RGB(255, 0, 0), 2);
				}
				else
				{
					cvCircle(pImg1, cvPoint(nx1, ny1), 1, CV_RGB(0, 0, 0), 1);
				}
				if(n > 0)
				{
					cvLine(pImg1, cvPoint(vPrevData[0], vPrevData[1]), cvPoint(nx1, ny1), CV_RGB(0, 0, 0), 1);
				}
				vPrevData[0] = nx1;
				vPrevData[1] = ny1;
			}
			cvSaveImage(szOutFileName1.c_str(), pImg1);
			cvReleaseImage(&pImg1);
		}
		return hr;
	}

	static HRESULT WriteFeatureAsDEBUGTRAJ(const vector< vector<double> >& vvData, string szOutFileName1, string szOutFileName2)
	{
		vector< vector<double> > vvNewData;
		size_t nSize = vvData.size();
		HRESULT hr = Makein01(vvData, vvNewData);
		if(hr == S_OK)
		{
			IplImage* pImg1 = cvCreateImage(cvSize(500, 500), IPL_DEPTH_8U, 3);
			IplImage* pImg2 = cvCreateImage(cvSize(500, 500), IPL_DEPTH_8U, 3);
			cvSet(pImg1, CV_RGB(255, 255, 255));
			cvSet(pImg2, CV_RGB(255, 255, 255));
			vector<int> vPrevData(4, 0);
			for(unsigned int n = 0; n < nSize; n++)
			{
				double dx1 = floor(vvNewData[n][0] * 500.0 + 0.5);
				double dy1 = floor(vvNewData[n][1] * 500.0 + 0.5);
				double dx2 = floor(vvNewData[n][2] * 500.0 + 0.5);
				double dy2 = floor(vvNewData[n][3] * 500.0 + 0.5);
				int nx1 = (int)dx1;
				int ny1 = (int)dy1;
				int nx2 = (int)dx2;
				int ny2 = (int)dy2;
				if(n == 0)
				{
					cvCircle(pImg1, cvPoint(nx1, ny1), 3, CV_RGB(0, 255, 0), 2);
					cvCircle(pImg2, cvPoint(nx2, ny2), 3, CV_RGB(0, 255, 0), 2);
				}
				else if(n == nSize - 1)
				{
					cvCircle(pImg1, cvPoint(nx1, ny1), 3, CV_RGB(255, 0, 0), 2);
					cvCircle(pImg2, cvPoint(nx2, ny2), 3, CV_RGB(255, 0, 0), 2);
				}
				else
				{
					cvCircle(pImg1, cvPoint(nx1, ny1), 1, CV_RGB(0, 0, 0), 1);
					cvCircle(pImg2, cvPoint(nx2, ny2), 1, CV_RGB(0, 0, 0), 1);
				}
				if(n > 0)
				{
					cvLine(pImg1, cvPoint(vPrevData[0], vPrevData[1]), cvPoint(nx1, ny1), CV_RGB(0, 0, 0), 1);
					cvLine(pImg2, cvPoint(vPrevData[2], vPrevData[3]), cvPoint(nx2, ny2), CV_RGB(0, 0, 0), 1);
				}
				vPrevData[0] = nx1;
				vPrevData[1] = ny1;
				vPrevData[2] = nx2;
				vPrevData[3] = ny2;
			}
			cvSaveImage(szOutFileName1.c_str(), pImg1);
			cvSaveImage(szOutFileName2.c_str(), pImg2);
			cvReleaseImage(&pImg1);
			cvReleaseImage(&pImg2);
		}
		return hr;
	}
	
	
	static HRESULT AccBackground(const CvArr* src, CvArr*& dst)
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
		cvAcc(src, dst);
		return S_OK;
	}

	// Divide the number of accumulated background by its number of frames
	static HRESULT AvgBackground(const CvArr* src, CvArr*& dst, unsigned int nnum)
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

	// This function use simply variance change magnitude to detect hand onset (start) and offset (end)
	static HRESULT TestOnlineBkgndFrames1Dot2(vector<IplImage*>& vpImgs, int& nDetStopBkGnd, int& nStop)
	{
		IplImage* pBkGnd = 0;
		IplImage* pcurr = 0;
		IplImage* pcurr2 = 0;
		IplImage* pdst = 0;
		//IplImage* pImg = 0;
		int nStart = 0;
		vector<double> vdbStdev;
		vector<double> vdbDifStdev;
		
		nStop = (int)vpImgs.size() - 1;
		// Accumulate background information
		double dbMostNegDif = 0.0;
		double dbMostPosDif = 0.0;
		double dbMaxStdev = 0.0;
		double dbMinStdev = DBL_MAX; //1.79769e+308;
		int nNewStop = 0;
		for(int n = nStart; n <= nStop; n++)
		{
			IplImage* pcurr = vpImgs[n];
			// OK, let's assume the first image is background, 
			// this is a big assumption, but there is nothing else we can do.
			if(n == nStart)			
			{
				pBkGnd = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
				CvScalar s0;
				s0.val[0] = s0.val[1] = s0.val[2] = s0.val[3] = 0;
				cvSet(pBkGnd, s0);
				GMMAdaptTool::AccBackground((const CvArr*)pcurr, (CvArr*&)pBkGnd);
			}
			if(n > nStart)
			{
				pcurr2 = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
				cvConvert(pcurr, pcurr2);
				pdst = cvCreateImage(cvSize(pcurr->width, pcurr->height), IPL_DEPTH_32F, 1);
				cvSub(pcurr2, pBkGnd, pdst);
				CvScalar std_dev, mean1;
				cvAvgSdv(pdst, &mean1, &std_dev);
				cvReleaseImage(&pdst);
				cvReleaseImage(&pcurr2);
				vdbStdev.push_back(std_dev.val[0]);
				if(n - nStart > 1)
				{
					double dbDif = vdbStdev[n - nStart - 1] - vdbStdev[n - nStart - 2];
					if(dbDif > 0 && dbDif > dbMostPosDif)
					{
						dbMostPosDif = dbDif;
					}
					if(dbDif < 0 && dbDif < dbMostNegDif)
					{
						dbMostNegDif = dbDif;
					}
					vdbDifStdev.push_back(dbDif);
				}
				if(std_dev.val[0] < dbMinStdev)
					dbMinStdev = std_dev.val[0];
				if(std_dev.val[0] > dbMaxStdev)
					dbMaxStdev = std_dev.val[0];
				//cvReleaseImage(&pdst);
			}
		}
		cvReleaseImage(&pBkGnd);
		nNewStop = nStop;
		for(int t = 0; t < (int)vdbStdev.size(); t++)
		{
			double db1 = (vdbStdev[t] - dbMinStdev) / (dbMaxStdev - dbMinStdev);
			if(db1 > 0.5)
			{
				nDetStopBkGnd = nStart + t - 1;
				if(nDetStopBkGnd < 1)
					nDetStopBkGnd = 1;
				break;
			}
		}
		for(int t = (int)vdbStdev.size() - 1; t >= 0; t--)
		{
			double db1 = (vdbStdev[t] - dbMinStdev) / (dbMaxStdev - dbMinStdev);
			if(db1 > 0.25)
			{
				nNewStop = nStart + t + 1;
				if(nNewStop > nStop)
					nNewStop = nStop;
				break;
			}
		}
		//for(int n = nStop; n > nNewStop; n--)
		//{
		//	pImg = vpImgs[n];
		//	vpImgs.pop_back();
		//	cvReleaseImage(&pImg);
		//}
		nStop = nNewStop;

		if(nDetStopBkGnd == 0)
		{
			return ERROR_EMPTY;
		}
		vdbStdev.clear();
		vdbDifStdev.clear();
		return S_OK;
	}	
	
	static HRESULT TestOnlineBkgndFrames1Dot2(vector<cv::Mat>& vImgs, int& nDetStopBkGnd, int& nStop)
	{
		cv::Mat BkGnd;
		cv::Mat curr2;
		cv::Mat dst;
		//IplImage* pImg = 0;
		int nStart = 0;
		vector<double> vdbStdev;
		vector<double> vdbDifStdev;
		
		nStop = (int)vImgs.size() - 1;
		// Accumulate background information
		double dbMostNegDif = 0.0;
		double dbMostPosDif = 0.0;
		double dbMaxStdev = 0.0;
		double dbMinStdev = DBL_MAX; //1.79769e+308;
		int nNewStop = 0;
		for(int n = nStart; n <= nStop; n++)
		{
			//curr = vImgs[n];
			// OK, let's assume the first image is background, 
			// this is a big assumption, but there is nothing else we can do.
			if(n == nStart)			
			{
				vImgs[n].convertTo(BkGnd, CV_64F);
			}
			if(n > nStart)
			{
				//curr2 = cv::Mat(
				
				vImgs[n].convertTo(curr2, CV_64F);
				dst = curr2 - BkGnd;
				cv::Scalar std_dev, mean1;
				cv::meanStdDev(dst, mean1, std_dev);
				
				vdbStdev.push_back(std_dev.val[0]);
				if(n - nStart > 1)
				{
					double dbDif = vdbStdev[n - nStart - 1] - vdbStdev[n - nStart - 2];
					if(dbDif > 0 && dbDif > dbMostPosDif)
					{
						dbMostPosDif = dbDif;
					}
					if(dbDif < 0 && dbDif < dbMostNegDif)
					{
						dbMostNegDif = dbDif;
					}
					vdbDifStdev.push_back(dbDif);
				}
				if(std_dev.val[0] < dbMinStdev)
					dbMinStdev = std_dev.val[0];
				if(std_dev.val[0] > dbMaxStdev)
					dbMaxStdev = std_dev.val[0];
				//cvReleaseImage(&pdst);
			}
		}
		nNewStop = nStop;
		for(int t = 0; t < (int)vdbStdev.size(); t++)
		{
			double db1 = (vdbStdev[t] - dbMinStdev) / (dbMaxStdev - dbMinStdev);
			if(db1 > 0.5)
			{
				nDetStopBkGnd = nStart + t - 1;
				if(nDetStopBkGnd < 1)
					nDetStopBkGnd = 1;
				break;
			}
		}
		for(int t = (int)vdbStdev.size() - 1; t >= 0; t--)
		{
			double db1 = (vdbStdev[t] - dbMinStdev) / (dbMaxStdev - dbMinStdev);
			if(db1 > 0.25)
			{
				nNewStop = nStart + t + 1;
				if(nNewStop > nStop)
					nNewStop = nStop;
				break;
			}
		}
		nStop = nNewStop;

		if(nDetStopBkGnd == 0)
		{
			return ERROR_EMPTY;
		}
		vdbStdev.clear();
		vdbDifStdev.clear();
		return S_OK;
	}	


	//========================================================================================//
	// Assume all images in each sequence are included in vImageNames.
	// vImageNames: Input, path names of all images
	// hSeqStart: Output, estimated image index where hand occur
	// hSeqStop: Output, estimated image index where hand disappear
	// hSeqs: Output, hash table for each name of sequence
	// vSeqs: Output, key value for image sequence in order
	//========================================================================================//
	static void TestImageStartStop(vector<string>& vImageNames, unordered_map<string, int>& hSeqStart, unordered_map<string, int>& hSeqStop, 
		unordered_map<string, vector<string> >& hvSeqs,
		vector<string>& vSeqs)
	{
		size_t nNoImgs = vImageNames.size();
		size_t nNoSeqs = 0;
		for(size_t n = 0; n < nNoImgs; n++)
		{
			//================================================================================//
			// Get the name of sequence
			//================================================================================//
			string szCurrFilePath = vImageNames[n];
			int nPos001s = szCurrFilePath.find_last_of("/");
			int nPos001b = szCurrFilePath.find_last_of("\\");
			int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
			int nPos002 = szCurrFilePath.find_last_of(".");
			string szCurrFileName = szCurrFilePath.substr(nPos001 + 1, nPos002 - nPos001 - 1);
			int nPos0001seq = szCurrFileName.find_last_of('_');
			string szCurrSeqKey = szCurrFileName.substr(0, nPos0001seq + 1);
			//================================================================================//
			if(hvSeqs.find(szCurrSeqKey) == hvSeqs.end())
			{
				vSeqs.push_back(szCurrSeqKey);
			}
			hvSeqs[szCurrSeqKey].push_back(szCurrFilePath);
		}
		nNoSeqs = vSeqs.size();
		for(size_t n = 0; n < nNoSeqs; n++)
		{
			vector<IplImage*> vpImgs;
			size_t nImgINSeq = hvSeqs[vSeqs[n]].size();
			for(size_t i = 0; i < nImgINSeq; i++)
			{
				// Load and store the pointer of gray scale version of images
				string szni = hvSeqs[vSeqs[n]][i];
				IplImage* pImg = cvLoadImage(szni.c_str(), 0);
				IplImage* pImgq = cvCreateImage(cvSize(64, 48), pImg->depth, pImg->nChannels);
				//IplImage* pImgq2 = cvCreateImage(cvSize(pImg->width, pImg->height), pImg->depth, pImg->nChannels);
				int nsizeW = (int)floor((double)pImg->width / 64.0 + 0.5);
				int nsizeH = (int)floor((double)pImg->height / 48.0 + 0.5);
				//================================================================================//
				// NO LPF
				//================================================================================//
				//CvMat* pFilter = cvCreateMat(nsizeH,nsizeW,CV_32FC1);
				//for(int n1 = 0; n1 < nsizeH; n1++)
				//	for(int n2 = 0; n2 < nsizeW; n2++) 
				//		cvSetReal2D(pFilter, n1, n2, 1.0/(nsizeH*nsizeW));
				//================================================================================//
				//cvFilter2D(pImg, pImgq2, pFilter);
				//cvResize(pImgq2, pImgq, CV_INTER_AREA);
				cvResize(pImg, pImgq, CV_INTER_AREA);
				vpImgs.push_back(pImgq);
				//cvReleaseMat(&pFilter);
				cvReleaseImage(&pImg);
				//cvReleaseImage(&pImgq2);
			}
			int nStopBkGnd = 0;
			int nStop = (int)nImgINSeq - 1;
			//================================================================================//
			// Actually test where the background start and stop
			//================================================================================//
			TestOnlineBkgndFrames1Dot2(vpImgs, nStopBkGnd, nStop);
			hSeqStart[vSeqs[n]] = nStopBkGnd + 1;
			hSeqStop[vSeqs[n]] = nStop;
			//================================================================================//
			for(size_t i = 0; i < vpImgs.size(); i++)
			{
				IplImage* pImg = vpImgs[i];
				cvReleaseImage(&pImg);
				vpImgs[i] = NULL;
			}
		}	
	}

	//========================================================================================//
	// Assume all images in each sequence are included in vImageNames.
	// vImageNames: Input, path names of all images
	// hJPGLabels: Input, information for each image sequences
	// hSeqStart: Output, exact image index where hand occur
	// hSeqStop: Output, exact image index where hand disappear
	// hSeqs: Output, hash table for each name of sequence
	// vSeqs: Output, key value for image sequence in order
	//========================================================================================//
	static void GetImageStartStop(vector<string>& vImageNames, unordered_map<string, JPGLabelInfo>& hJPGLabels, unordered_map<string, int>& hSeqStart, unordered_map<string, int>& hSeqStop,
		unordered_map<string, vector<string> >& hvSeqs,
		vector<string>& vSeqs)
	{
		size_t nNoImgs = vImageNames.size();
		size_t nNoSeqs = 0;
		for(size_t n = 0; n < nNoImgs; n++)
		{
			//================================================================================//
			// Get the name of sequence
			//================================================================================//
			string szCurrFilePath = vImageNames[n];
			int nPos001s = szCurrFilePath.find_last_of("/");
			int nPos001b = szCurrFilePath.find_last_of("\\");
			int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
			int nPos002 = szCurrFilePath.find_last_of(".");
			string szCurrFileName = szCurrFilePath.substr(nPos001 + 1, nPos002 - nPos001 - 1);
			int nPos0001seq = szCurrFileName.find_last_of('_');
			string szCurrSeqKey = szCurrFileName.substr(0, nPos0001seq + 1);
			//================================================================================//
			if(hvSeqs.find(szCurrSeqKey) == hvSeqs.end())
			{
				vSeqs.push_back(szCurrSeqKey);
			}
			hvSeqs[szCurrSeqKey].push_back(szCurrFilePath);
		}
		nNoSeqs = vSeqs.size();
		for(size_t n = 0; n < nNoSeqs; n++)
		{
			size_t nImgINSeq = hvSeqs[vSeqs[n]].size();
			bool bFoundStart = false;
			for(size_t i = 0; i < nImgINSeq; i++)
			{
				string szCurrFilePath = hvSeqs[vSeqs[n]][i];
				
				if(!bFoundStart && hJPGLabels[szCurrFilePath].m_vvOnes.size())
				{
					bFoundStart = true;
					hSeqStart[vSeqs[n]] = i;
				}
				else if(bFoundStart && 0 == hJPGLabels[szCurrFilePath].m_vvOnes.size())
				{
					//bFoundStop = true;
					hSeqStop[vSeqs[n]] = i - 1;
				}
			}
		}
	}

	static void GetFileNameInfo(string& szImagePathName, string& szImagePath, string& szImageName, string& szSeqName, string& szExt)
	{
		//================================================================================//
		// Get the name of sequence
		//================================================================================//
		int nPos001s = szImagePathName.find_last_of("/");
		int nPos001b = szImagePathName.find_last_of("\\");
		int nPos001 = (nPos001s>nPos001b)?nPos001s:nPos001b;
		int nPos002 = szImagePathName.find_last_of(".");
		szImagePath = szImagePathName.substr(0, nPos001 + 1);
		szImageName = szImagePathName.substr(nPos001 + 1, nPos002 - nPos001 - 1);
		int nPos0001seq = szImageName.find_last_of('_');
		szSeqName = szImageName.substr(0, nPos0001seq + 1);
		szExt = szImagePathName.substr(nPos002);
		//================================================================================//
	}

	//========================================================================================//
	// Assign empty label for those not listed in hJPGLavelInfo but was listed in 
	// the file szRefFilePath
	// hJPGLavelInfo: Input/Output, hash table with label information on list of files
	// hJPGClasses: Input/Output, hash table to store the ground truth classes
	// vKeysInOrder: Output, key(file paths) in the order of listed file
	// szRefFilePath: Input, path name of the file list with list of image paths
	// BLKSIZE: Input, block size
	//========================================================================================//
	static void FixlabelInfoByRefFileList(LoadJPGLabelInfo& lb, string szRefFilePath, int BLKSIZE)
	{
		return FixlabelInfoByRefFileList(lb.m_hJPGLavelInfo, lb.m_hJPGClasses, lb.m_vszKeyOrder, szRefFilePath, BLKSIZE);
	}
	static void FixlabelInfoByRefFileList(unordered_map<string, JPGLabelInfo>& hJPGLavelInfo, unordered_map< string, vector<int> >& hJPGClasses, vector<string>& vKeysInOrder, string szRefFilePath, int BLKSIZE)
	{
		// First get the list of sequence:
		unordered_map<string, int> hSeqsCount;
		unordered_map<string, JPGLabelInfo>::iterator itr;
		for(itr = hJPGLavelInfo.begin(); itr != hJPGLavelInfo.end(); itr++)
		{
			string szCurrPathName = itr->second.m_szPathName;
			string szImagePath;
			string szImageName;
			string szSeqName;
			string szExt;
			GetFileNameInfo(szCurrPathName, szImagePath, szImageName, szSeqName, szExt);
			hSeqsCount[szSeqName]++;
		}
		ifstream fin(szRefFilePath.c_str());
		vector<string> vszRefFileLists;
		string szBuf;
		if(!fin.is_open())
		{
			cerr<<"cannot open reference file list: "<<szRefFilePath<<endl;
			return;
		}
		vKeysInOrder.clear();
		while(getline(fin, szBuf))
		{
			if(szBuf[szBuf.length() - 1] == '\r')
			{
				szBuf = szBuf.substr(0, szBuf.length()-1);
			}
			string szImagePath;
			string szImageName;
			string szSeqName;
			string szExt;
			GetFileNameInfo(szBuf, szImagePath, szImageName, szSeqName, szExt);
			if(hSeqsCount.find(szSeqName) != hSeqsCount.end())
			{
				vKeysInOrder.push_back(szBuf);
				if(hJPGLavelInfo.find(szBuf) == hJPGLavelInfo.end())
				{
					IplImage* pImg = cvLoadImage(szBuf.c_str());
					hJPGLavelInfo[szBuf].SetWidthHeight(pImg->width, pImg->height, BLKSIZE);
					hJPGLavelInfo[szBuf].m_szPathName = szBuf;
					int nBlockWidth = (int)floor((double)pImg->width/(double)BLKSIZE);
					int nBlockHeight = (int)floor((double)pImg->height/(double)BLKSIZE);
					hJPGClasses[szBuf] = vector<int>(nBlockWidth*nBlockHeight, 0);
					cvReleaseImage(&pImg);
				}
			}
		}
		fin.close();
	}

	static void FixSLIClabelInfoByRefFileList(LoadJPGSLICLabelInfo& lb, string szRefFilePath)
	{
		return FixSLIClabelInfoByRefFileList(lb.m_hJPGSLICLabelInfo, lb.m_vszKeyOrder, szRefFilePath);
	}
	static void FixSLIClabelInfoByRefFileList(unordered_map<string, unordered_map<int, int> >& hJPGSLICLabelInfo, vector<string>& vKeysInOrder, string szRefFilePath)
	{
		// First get the list of sequence:
		unordered_map<string, int> hSeqsCount;
		unordered_map<string, unordered_map<int, int> >::iterator itr;
		for(itr = hJPGSLICLabelInfo.begin(); itr != hJPGSLICLabelInfo.end(); itr++)
		{
			string szCurrPathName = itr->first;
			string szImagePath;
			string szImageName;
			string szSeqName;
			string szExt;
			GetFileNameInfo(szCurrPathName, szImagePath, szImageName, szSeqName, szExt);
			hSeqsCount[szSeqName]++;
		}
		ifstream fin(szRefFilePath.c_str());
		vector<string> vszRefFileLists;
		string szBuf;
		if(!fin.is_open())
		{
			cerr<<"cannot open reference file list: "<<szRefFilePath<<endl;
			return;
		}
		vKeysInOrder.clear();
		while(getline(fin, szBuf))
		{
			if(szBuf[szBuf.length() - 1] == '\r')
			{
				szBuf = szBuf.substr(0, szBuf.length()-1);
			}
			string szImagePath;
			string szImageName;
			string szSeqName;
			string szExt;
			GetFileNameInfo(szBuf, szImagePath, szImageName, szSeqName, szExt);
			if(hSeqsCount.find(szSeqName) != hSeqsCount.end())
			{
				vKeysInOrder.push_back(szBuf);
				if(hJPGSLICLabelInfo.find(szBuf) == hJPGSLICLabelInfo.end())
				{
					hJPGSLICLabelInfo[szBuf] = unordered_map<int, int>();
				}
			}
		}
		fin.close();
	}


	static void AccBackGroundbyIndex(vector<IplImage*>& vpImgs, int nStart, int nStop, IplImage*& pAccImg)
	{
		if(0 == vpImgs.size())
		{
			cerr<<"Empty vector of images"<<endl;
			return;
		}
		
		if(nStop < nStart || nStop >= vpImgs.size() || nStart < 0)
		{
			cerr<<"Invalid Indices"<<endl;
			return;
		}
		int nwidth = vpImgs[nStart]->width;
		int nheight = vpImgs[nStart]->height;
		int ndepth = vpImgs[nStart]->depth;
		int nchannels = vpImgs[nStart]->nChannels;
		for(int n = nStart; n <= nStop; n++)
		{
			if(nwidth != vpImgs[n]->width || nheight != vpImgs[n]->height || 
				ndepth != vpImgs[n]->depth || nchannels != vpImgs[n]->nChannels)
			{
				cerr<<"Size mismatch at index = "<<n<<" of current stream"<<endl;
				return;
			}
		}
		if(pAccImg)
		{
			cvReleaseImage(&pAccImg);
		}
		
		IplImage* pAccImgTmp = cvCreateImage(cvSize(nwidth, nheight), IPL_DEPTH_32F, nchannels);
		IplImage* pAccImgTmpTmp = cvCreateImage(cvSize(nwidth, nheight), IPL_DEPTH_32F, nchannels);
		cvSetZero(pAccImgTmp);
		pAccImg = cvCreateImage(cvSize(nwidth, nheight), ndepth, nchannels);
		for(int n = nStart; n <= nStop; n++)
		{
			cvConvert(vpImgs[n], pAccImgTmpTmp);
			GMMAdaptTool::AccBackground((CvArr*)pAccImgTmpTmp, (CvArr*&)pAccImgTmp);
		}
		cvConvertScale(pAccImgTmp, pAccImg, 1.0);
		cvReleaseImage(&pAccImgTmp);
		cvReleaseImage(&pAccImgTmpTmp);
		//cvConvertScale(pAccImgTmp, pAccImg, (double>(1.0)/(double)(nStop - nStart + 1)));
	}

	static void AccBackGroundbyIndex(vector<cv::Mat>& vImgs, int nStart, int nStop, cv::Mat& AccImg)
	{
		if(0 == vImgs.size())
		{
			cerr<<"Empty vector of images"<<endl;
			return;
		}
		
		if(nStop < nStart || nStop >= vImgs.size() || nStart < 0)
		{
			cerr<<"Invalid Indices"<<endl;
			return;
		}
		
		int nwidth = vImgs[nStart].cols;
		int nheight = vImgs[nStart].rows;
		int ndepth = vImgs[nStart].depth();
		int nchannels = vImgs[nStart].channels();
		for(int n = nStart; n <= nStop; n++)
		{
			if(nwidth != vImgs[n].cols || nheight != vImgs[n].rows || 
				ndepth != vImgs[n].depth() || nchannels != vImgs[n].channels())
			{
				cerr<<"Size mismatch at index = "<<n<<" of current stream"<<endl;
				return;
			}
		}
		
		IplImage* pAccImgTmp = cvCreateImage(cvSize(nwidth, nheight), IPL_DEPTH_32F, nchannels);
		IplImage* pAccImgTmpTmp = cvCreateImage(cvSize(nwidth, nheight), IPL_DEPTH_32F, nchannels);
		cvSetZero(pAccImgTmp);
		//IplImage* pAccImg = cvCreateImage(cvSize(nwidth, nheight), ndepth, nchannels);
		for(int n = nStart; n <= nStop; n++)
		{
			IplImage img = vImgs[n].operator IplImage();
			cvConvert(&img, pAccImgTmpTmp);
			GMMAdaptTool::AccBackground((CvArr*)pAccImgTmpTmp, (CvArr*&)pAccImgTmp);
		}
		cv::Mat tmp(pAccImgTmp);
		AccImg = tmp * ((double)(1.0)/(double)(nStop - nStart + 1));
		//cvConvertScale(pAccImgTmp, pAccImg, 1.0);

		cvReleaseImage(&pAccImgTmp);
		cvReleaseImage(&pAccImgTmpTmp);
		//cvConvertScale(pAccImgTmp, pAccImg, (double>(1.0)/(double)(nStop - nStart + 1)));
	}

	// Normalize image color by RGB histogram
	static void NormalizeImageRGB(IplImage* pImg, IplImage*& pNewImg, double ptr = 0.01)
	{
		int height = pImg->height;
		int width  = pImg->width;
		int nC = pImg->nChannels;
		int depth = pImg->depth;
		double dPts = (double)(width * height);
		double dptrnum = dPts * ptr;
		if(pNewImg != NULL)
		{
			if(!(pNewImg->height == height && pNewImg->width == width && 
				pNewImg->nChannels == nC && pNewImg->depth == depth))
			{
				cvReleaseImage(&pNewImg);
				pNewImg = cvCreateImage(cvSize(width, height), depth, nC);
			}
		}
		else
		{
			pNewImg = cvCreateImage(cvSize(width, height), depth, nC);
		}
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
		}
	}

	static int NormalizeImageSeq(vector<cv::Mat>& vImgs, vector<cv::Mat>& vNewImgs, double dbtr = 0.01)
	{
		vector< vector<double> > vHist;
		double dblight = 0.0;
		if(0<GetImgesHistRGB(vImgs, vHist, dblight, dbtr))
		{
			vNewImgs.reserve(vImgs.size());
			vNewImgs.resize(vImgs.size());
			for(size_t n = 0; n < vImgs.size(); ++n)
			{
				vNewImgs[n] = cv::Mat(vImgs[n].rows, vImgs[n].cols, CV_8UC3);
				for(int l = 0; l < vImgs[n].rows; ++l)
				{
					for(int k = 0; k < vImgs[n].cols; ++k)
					{
						for(int c = 0; c < vImgs[n].channels(); ++c)
						{
							double dbval = (double)vImgs[n].at<Vec3b>(l,k)[c];
							if(dbval >= dblight)
								dbval = 255.0;
							else
								dbval = dbval / dblight * 255;
							vNewImgs[n].at<Vec3b>(l,k)[c] = (uchar) floor(dbval+0.5);
						}
					}
				}
			}
		}
		return 0;
	}


	static int GetImgesHistRGB(vector<cv::Mat>& vImgs, vector< vector<double> >& vHist, double& dblight, double dbtr = 0.01)
	{
		if(0 >= (int)vImgs.size())
		{
			cerr<<"No image to be counting\n";
			return -1;
		}
		int nC = vImgs[0].channels();
		for(size_t n = 1; n < vImgs.size(); ++n)
		{
			int nCn = vImgs[n].channels();
			if(nCn != nC)
			{
				cerr<<"Images have different number of channels\n";
				return -1;
			}
		}
		vHist.resize(nC);
		for(int c = 0; c < nC; ++c)
		{
			vHist[c] = vector<double>(256,0.0);
		}
		double dbptrnum = 0.0;
		for(size_t n = 0; n < vImgs.size(); ++n)
		{
			int width = vImgs[n].cols;
			int height = vImgs[n].rows;
			dbptrnum += (double)(width * height);
			//IplImage ipltmp = vImgs[n].operator IplImage();
			//			//IplImage* pipltmp = &ipltmp;
			for(int l = 0; l < height; ++l)
			{
				for(int k = 0; k < width; ++k)
				{
					//CvScalar slk = cvGet2D(pipltmp, l, k);
					for(int c = 0; c < nC; ++c)
					{
						uchar ucVal = vImgs[n].at<Vec3b>(l,k)[c];
						//unsigned int unVal = (unsigned int) ucVal;
						double dbVal = (double) ucVal;
						vHist[c][(int)floor(dbVal+0.5)]++;
					}
				}

			}
		}
		dbptrnum *= dbtr;
		vector<double> vtr(nC, 0.0);
		dblight = 0.0;
		for(int c = 0; c < nC; ++c)
		{
			double dbtmpcb = 0.0;
			double dbtmptr = 255.0;
			for(int b = 255; b >=0; b--)
			{
				dbtmpcb += vHist[c][b];
				if(dbtmpcb >= dbptrnum)
				{
					dbtmptr = (double)b;
					break;
				}
			}
			vtr[c] = dbtmptr;
			if(vtr[c] > dblight)
				dblight = vtr[c];
		}
		return 1;
	}


	//static void LoadMyBinMat(string szPathName, matrix<double>& mat)
	//{
	//	FILE* f = NULL;
	//	int n1 = 0;
	//	int n2 = 0;
	//	double dbBuf = 0.0;
	//	f = fopen(szPathName.c_str(), "rb");
	//	if(!f)
	//	{
	//		cerr<<"Cannot open "<<szPathName<<endl;
	//		return;
	//	}
	//	fread(&n1, sizeof(int), 1, f);
	//	fread(&n2, sizeof(int), 1, f);
	//	mat.SetSize(n1, n2);
	//	for(int i1 = 0; i1 < n1; ++i1)
	//	{
	//		for(int i2 = 0; i2 < n2; ++i2)
	//		{
	//			fread(&dbBuf, sizeof(double), 1, f);
	//			mat(i1, i2) = dbBuf;
	//		}
	//	}
	//	fclose(f);
	//}
	//
	//static void WriteMyBinMat(string szPathName, matrix<double>& mat)
	//{
	//	FILE* f = NULL;
	//	int n1 = (int)mat.RowNo();
	//	int n2 = (int)mat.ColNo();
	//	double dbBuf = 0.0;
	//	f = fopen(szPathName.c_str(), "wb");
	//	if(!f)
	//	{
	//		cerr<<"Cannot open "<<szPathName<<endl;
	//		return;
	//	}
	//	fwrite(&n1, sizeof(int), 1, f);
	//	fwrite(&n2, sizeof(int), 1, f);
	//	for(int i1 = 0; i1 < n1; ++i1)
	//	{
	//		for(int i2 = 0; i2 < n2; ++i2)
	//		{
	//			dbBuf = mat(i1, i2);
	//			fwrite(&dbBuf, sizeof(double), 1, f);
	//		}
	//	}
	//	fclose(f);
	//}
	
	// static int GenEpsilon(matrix<double>& mat, int nDims, double eps);
	// static int GMMBiasAdapt(vector< vector<double> >& vvData, matrix<double>& newmub, vector<GMM>& oldGMMs, vector<GMM>& newGMMs);
	// static int GMMBiasAdapt2(vector< vector<double> >& vvData, matrix<double>& Sigmab, matrix<double>& newmub, matrix<double>& newSigmab, vector<GMM>& oldGMMs, vector<GMM>& newGMMs, int nIter);
	// static int GMMHash2VectBH(unordered_map<string, GMM>& hGMMs, vector<GMM>& vGMMs);
	// static int Test2ClassesGMMs(vector< vector<double> >& vvData, vector<int>& classes, vector<GMM>& GMMs, vector<int>& recogClasses, double& Acc, vector<int>& confres);
	// static int Test2ClassesGMMs(vector< vector<double> >& vvData, vector<GMM>& GMMs, vector<int>& recogClasses);
	// static int Test2ClassesGMMs(string GMMPaths, string DataPaths, double& Acc, int nFrameSize);
	// static int Test2ClassesGMMsImg(string GMMPaths, string DataPaths, double& Acc, int nFrameSize);
	// static int TestBoosted2ClassesGMMs(vector<string> vGMMPaths, string DataPaths, double& Acc, vector<double>& vAlpha);
	// static int TestBoosted2ClassesGMMs(vector< vector<double> >& vvData, vector<int>& vClasses, vector< vector<GMM> >& vGMMs, vector<double>& vAlpha, vector<int>& recogClasses, double& Acc);
	// static int GMMAdaptRcogbyFrame(vector<GMM>& oldGMMs, vector< vector<double> >& vvData, vector<int>& vclasses, double& Acc, int nFrameSize);
	// static int GMMAdaptRcogbyFrame2(vector<GMM>& oldGMMs, vector< vector<double> >& vvData, vector<int>& vclasses, double& Acc, int nFrameSize);
	// static int GMMAdaptRcogbyFrame(string GMMPaths, string DataPaths, double& Acc, int nFrameSize);
	// static int GMMAdaptRcogbyFrame2(string GMMPaths, string DataPaths, double& Acc, int nFrameSize);
	// static int GMMHTKMAPBySCP(string GMMPaths, string DataSCPPaths, double& Acc);
	// static int GMMHTKMAPBySCP2(string GMMPaths, string DataSCPPaths, double& Acc);
	// static int GMMByLAB(string GMMPaths, string DataLABPaths, double& Acc);
	// static int GMMHTKMAPByLAB(string GMMPaths, string DataLABPaths, vector<double>& metrics);
	// static int GMMHTKMAPByLAB2(string GMMPaths, string DataLABPaths, vector<double>& metrics);
	// static int GMMHTKMAPByLABTestNonSkinAdaptOnly(string GMMPaths, string DataLABPaths, string szRefImgList, double& Acc, ofstream& fout);
	// static void CompMean(vector< vector<double> >& vvData, vector<double>& vMean);
	// static void CompDiagVar(vector< vector<double> >& vvData, matrix<double>& diagVar);
	// static void GetLargestConnectedComponents(vector<int>& orig, vector<int>&refined);	
	// static long ReadHTKMFCCFeature(vector< vector<double> >& vvData, string szInFileName );

	//void GMMAdaptTool::MarkSelectedPointWithoutDrawing(string szCurrFile, int x, int y, int val, 
	//	unordered_map<string, JPGLabelInfo>& hJPGLavelInfo, 
	//	int& nCountO,
	//	int& nCountX,
	//	int nWidth,
	//	int nHeight,
	//	int BLKSIZE)
	//{
	//	if(hJPGLavelInfo.find(szCurrFile) == hJPGLavelInfo.end())
	//	{
	//		hJPGLavelInfo[szCurrFile].SetWidthHeight(nWidth*BLKSIZE, nHeight*BLKSIZE, BLKSIZE);
	//		hJPGLavelInfo[szCurrFile].m_szPathName = szCurrFile;
	//	}
	//	if(val == 0)
	//	{
	//		hJPGLavelInfo[szCurrFile].SetAsNonSkinPixel(x, y);
	//	}
	//	else if(val == 1)
	//	{
	//		if(hJPGLavelInfo[szCurrFile].m_ppnLabel[y][x] == 1)
	//		{
	//			hJPGLavelInfo[szCurrFile].SetAsNonSkinPixel(x, y);
	//			nCountO--;
	//		}
	//		else
	//		{
	//			if(hJPGLavelInfo[szCurrFile].m_ppnLabel[y][x] == -1)
	//			{
	//				nCountX--;
	//			}
	//			hJPGLavelInfo[szCurrFile].SetAsSkinPixel(x, y);
	//			nCountO++;
	//		}
	//	}
	//	else if(val == -1)
	//	{
	//		if(hJPGLavelInfo[szCurrFile].m_ppnLabel[y][x] == -1)
	//		{
	//			hJPGLavelInfo[szCurrFile].SetAsNonSkinPixel(x, y);
	//			nCountX--;
	//		}
	//		else
	//		{
	//			if(hJPGLavelInfo[szCurrFile].m_ppnLabel[y][x] == 1)
	//			{
	//				nCountO--;
	//			}
	//			hJPGLavelInfo[szCurrFile].SetAsNotSurePixel(x, y);
	//			nCountX++;
	//		}
	//	}
	//}
	//
	//void GMMAdaptTool::SaveJPGInfo(string szLabelPathName, 
	//	unordered_map<string, JPGLabelInfo>& hJPGLavelInfo,
	//	vector<string>& vszJPGFiles,
	//	//unordered_map<string, int>& hszCurrJPGFiles,
	//	unordered_map<string, int>& hszJPGFiles)
	//{
	//	ofstream fout;
	//	fout.open(szLabelPathName.c_str());
	//	size_t stFiles = hJPGLavelInfo.size();
	//	size_t stCurrFileLists = vszJPGFiles.size();
	//	vector<string> vszOthers;
	//	int nIndex = 0;
	//	// Compare list of files:
	//	//for(hash_map<string, JPGLabelInfo>::iterator itr = hJPGLavelInfo.begin();itr!=hJPGLavelInfo.end(); itr++)
	//	//{
	//	//	string szKey = itr->first;
	//	//	if(hszCurrJPGFiles.find(szKey) == hszCurrJPGFiles.end())
	//	//	{
	//	//		vszOthers.push_back(szKey);
	//	//	}
	//	//}
	//	for(int nnn = 0; nnn < (int)vszJPGFiles.size(); nnn++)
	//	{
	//		if(hJPGLavelInfo.find(vszJPGFiles[nnn]) != hJPGLavelInfo.end())
	//		{
	//			nIndex = nnn;
	//		}
	//	}
	//	//for(hash_map<string, JPGLabelInfo>::iterator itr = m_hJPGLavelInfo.begin();itr!=m_hJPGLavelInfo.end(); itr++)
	//	for(size_t st = 0; st < (size_t)nIndex; st++)
	//	{
	//		//string szKey = itr->first;
	//		string szKey = vszJPGFiles[st];
	//		fout<<hJPGLavelInfo[szKey].m_szPathName<<","<<hJPGLavelInfo[szKey].m_nBlockWidth
	//			<<","<<hJPGLavelInfo[szKey].m_nBlockHeight<<","<<hszJPGFiles[szKey]<<endl;
	//		for(int mm = 0; mm < hJPGLavelInfo[szKey].m_nBlockWidth; mm++)
	//		{
	//			for(int nn = 0; nn < hJPGLavelInfo[szKey].m_nBlockHeight; nn++)
	//			{
	//				if(hJPGLavelInfo[szKey].m_ppnLabel[nn][mm] == 1)
	//				{
	//					fout<<mm<<","<<nn<<","<<1<<endl;
	//				}
	//			}
	//		}
	//		for(int mm = 0; mm < hJPGLavelInfo[szKey].m_nBlockWidth; mm++)
	//		{
	//			for(int nn = 0; nn < hJPGLavelInfo[szKey].m_nBlockHeight; nn++)
	//			{
	//				if(hJPGLavelInfo[szKey].m_ppnLabel[nn][mm] == -1)
	//				{
	//					fout<<mm<<","<<nn<<","<<-1<<endl;
	//				}
	//			}
	//		}
	//		fout<<"."<<endl;
	//	}
	//	// Others
	//	//for(size_t st = 0; st < vszOthers.size(); st++)
	//	//{
	//	//	string szKey = vszOthers[st];
	//	//	fout<<hJPGLavelInfo[szKey].m_szPathName<<","<<hJPGLavelInfo[szKey].m_nBlockWidth
	//	//		<<","<<hJPGLavelInfo[szKey].m_nBlockHeight<<","<<hszJPGFiles[szKey]<<endl;
	//	//	for(int mm = 0; mm < hJPGLavelInfo[szKey].m_nBlockWidth; mm++)
	//	//	{
	//	//		for(int nn = 0; nn < hJPGLavelInfo[szKey].m_nBlockHeight; nn++)
	//	//		{
	//	//			if(hJPGLavelInfo[szKey].m_ppnLabel[nn][mm] == 1)
	//	//			{
	//	//				fout<<mm<<","<<nn<<","<<1<<endl;
	//	//			}
	//	//		}
	//	//	}
	//	//	for(int mm = 0; mm < hJPGLavelInfo[szKey].m_nBlockWidth; mm++)
	//	//	{
	//	//		for(int nn = 0; nn < hJPGLavelInfo[szKey].m_nBlockHeight; nn++)
	//	//		{
	//	//			if(hJPGLavelInfo[szKey].m_ppnLabel[nn][mm] == -1)
	//	//			{
	//	//				fout<<mm<<","<<nn<<","<<-1<<endl;
	//	//			}
	//	//		}
	//	//	}
	//	//	fout<<"."<<endl;
	//	//}
	//
	//	fout.close();
	//}
	//	
	/*
	void GMMAdaptTool::LoadJPGInfo(int BLOCKSIZE, int NUMFRAMETR, int& m_nCountX, int& m_nCountO,
							unordered_map<string, int>& m_hszJPGFiles,
							//unordered_map<string, int>& m_hszCurrJPGFiles,
							unordered_map<string, JPGLabelInfo>& m_hJPGLavelInfo,
							unordered_map< string, vector<int> >& m_hJPGClasses,
							vector<string>& m_vszKeyOrder,
							string& m_szLabelPathName,
							unordered_map<string, unordered_map<string, int> >& m_hhTrainFiles)
	{
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
	*/
		
	//static void LoadJPGInfo(int BLOCKSIZE, int NUMFRAMETRi, int& m_nCountX, int& m_nCountO, 
	//		unordered_map<string, int>& m_hszJPGFiles,
	//		//unordered_map<string, int>& m_hszCurrJPGFiles,
	//		unordered_map<string, JPGLabelInfo>& m_hJPGLavelInfo,
	//		unordered_map< string, vector<int> >& m_hJPGClasses, 
	//		vector<string>& m_vszKeyOrder,
	//		string& m_szLabelPathName,
	//		unordered_map<string, unordered_map<string, int> >& m_hhTrainFiles);
	// static void ExtColorFeature(string szFileName, JPGLabelInfo* pinfo, vector< vector<double> >& vvFeatures, int BLKSIZE);
	// static void ExtColorFeatureRGBOnly(string szFileName, JPGLabelInfo* pinfo, vector< vector<double> >& vvFeatures, int BLKSIZE);

	//static void MarkSelectedPointWithoutDrawing(string szCurrFile, int x, int y, int val, 
	//	unordered_map<string, JPGLabelInfo>& hJPGLavelInfo, 
	//	int& nCountO,
	//	int& nCountX,
	//	int nWidth,
	//	int nHeight,
	//	int BLKSIZE);
	//static void SaveJPGInfo(string szLabelPathName, 
	//	unordered_map<string, JPGLabelInfo>& hJPGLavelInfo,
	//	vector<string>& vszJPGFiles,
	//	//unordered_map<string, int>& hszCurrJPGFiles,
	//	unordered_map<string, int>& hszJPGFiles);
	// static long WriteFeatureAsHTKMFCC(const vector< vector<double> >& vvData, string szOutFileName, short sampKind);
	// static HRESULT GMMAdaptTool::AccBackground(const CvArr* src, CvArr*& dst);
	// static HRESULT TestOnlineBkgndFrames1Dot2(vector<IplImage*> vpImgs, int& nDetStopBkGnd, int& nStop); 
	// static void TestImageStartStop(vector<string>& vImageNames, unordered_map<string, int>& hSeqStart, unordered_map<string, int>& hSeqStop, 
	// unordered_map<string, vector<string> >& hvSeqs,
	// vector<string>& vSeqs);
	// static void GetImageStartStop(vector<string>& vImageNames, unordered_map<string, JPGLabelInfo>& hJPGLabels, unordered_map<string, int>& hSeqStart, unordered_map<string, int>& hSeqStop,
		// unordered_map<string, vector<string> >& hvSeqs,
		// vector<string>& vSeqs);
	// static void GetFileNameInfo(string& szImagePathName, string& szImagePath, string& szImageName, string& szSeqName, string& szExt);
	// static void FixlabelInfoByRefFileList(unordered_map<string, JPGLabelInfo>& hJPGLavelInfo, unordered_map< string, vector<int> >& hJPGClasses, vector<string>& vKeysInOrder, string szRefFilePath, int BLKSIZE);



	//int m_nCountX;
	//int m_nCountO;
	//unordered_map<string, int> m_hszJPGFiles;
	//unordered_map<string, int> m_hszCurrJPGFiles;
	//unordered_map<string, JPGLabelInfo> m_hJPGLavelInfo;
	//vector<string> m_vszKeyOrder;
	//string m_szLabelPathName;
	//unordered_map<string, unordered_map<string, int> > m_hhTrainFiles;
	static void AutoLabelSLICRGBOnly(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int nsp, int& nHandModelIndex)
	{
		vector< vector<double> > vvData;
		vector< vector<cv::Point> > vvSLICLabels;
		// 1. extract image blocks corresponding to each block
		cv::Size imsize;
		int nsps = ExtLowLevelImgFeature::ExtColorFeatureSLICRGBOnly(szPath, vvData, vvSLICLabels, nsp, imsize);

		vector<int> vrecogclasses;
		// classifying skin/nonskin for the blocks
		vrecogclasses.resize(vvData.size());
		Test2ClassesGMMsFull(vvData, vGMMs, vrecogclasses);
		size_t nResultsSize = vrecogclasses.size();
		mask = cv::Mat(imsize, CV_8UC1);
		int n_my_nn = 0;
		for(size_t nn = 0; nn < nsp; ++nn)
		{
			size_t n_sp_size_nn = vvSLICLabels[nn].size();
			if(n_sp_size_nn > 0)
			{
				for(size_t mm = 0; mm < vvSLICLabels[nn].size(); ++mm)
				{
					if(vrecogclasses[n_my_nn] == nHandModelIndex)
						mask.at<uchar>(vvSLICLabels[nn][mm].x, vvSLICLabels[nn][mm].y) = 255;
					else
						mask.at<uchar>(vvSLICLabels[nn][mm].x, vvSLICLabels[nn][mm].y) = 0;
				}
				n_my_nn++;
			}
		}
		//cvReleaseImage(&pimg);
		//cvReleaseImage(&pfimg);
	}

	static void AutoLabel_SLIC_BS0(string szPath, vector<GMM>& vGMMs, cv::Mat& mask, int nsp, int& nHandModelIndex, cv::Mat& ref)
	{
		vector< vector<double> > vvData;
		vector< vector<cv::Point> > vvSLICLabels;
		IplImage iplref = ref.operator IplImage();
		IplImage* pref = &iplref;
		cv::Size imsize;
		vector< vector<double> > vvCorrs;
		vector<double> vKSs;
		//int nsps = ExtLowLevelImgFeature::ExtColorFeatures_SLIC_BS0(szPath, vvData, vvSLICLabels, nsp, pref, vvCorrs, imsize);
		int nsps = ExtLowLevelImgFeature::ExtColorFeatures_SLIC_BS0(szPath, vvData, vvSLICLabels, nsp, pref, vvCorrs, vKSs, imsize);
		vector<int> vrecogclasses;
		// classifying skin/nonskin for thes blocks
		vrecogclasses.resize(vvData.size());
		Test2ClassesGMMsFull(vvData, vGMMs, vrecogclasses);
		size_t nResultsSize = vrecogclasses.size();
		mask = cv::Mat(imsize, CV_8UC1);
		int n_my_nn = 0;
		for(size_t nn = 0; nn < nsp; ++nn)
		{
			size_t n_sp_size_nn = vvSLICLabels[nn].size();
			
			if(n_sp_size_nn > 0)
			{
				double Corr = 0.1140 * vvCorrs[nn][0] + 0.5870 * vvCorrs[nn][1] + 0.2989 * vvCorrs[nn][2];
				for(size_t mm = 0; mm < vvSLICLabels[nn].size(); ++mm)
				{
					if(vrecogclasses[n_my_nn] == nHandModelIndex && Corr < CCThreshold)
						mask.at<uchar>(vvSLICLabels[nn][mm].x, vvSLICLabels[nn][mm].y) = 255;
					else
						mask.at<uchar>(vvSLICLabels[nn][mm].x, vvSLICLabels[nn][mm].y) = 0;
				}
				n_my_nn++;
			}
		}
	}
	
	//static long ExtColorFeaturesFromImagesStartStop_SLIC(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, vector< vector< vector<cv::Point> > >& vvvSLICLabels, int nsp, int nStart, int nStop, IplImage*& pRefImg, vector< vector< vector<double> > >& vvvCorrs, cv::Size& imsize)
	static long ExtColorFeaturesFromImagesStartStop_SLIC(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, vector< vector< vector<cv::Point> > >& vvvSLICLabels, int nsp, int nStart, int nStop, IplImage*& pRefImg, vector< vector< vector<double> > >& vvvCorrs, vector< vector<double> >& vvKSs, cv::Size& imsize, bool SLIC_is_given = false)
	{
		long N = (long) vImgs.size();
		if(N <= 0)
		{
			cerr<<"No images for feature extraction"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvvFeatures.resize(N);
		vvvCorrs.resize(N);
		vvKSs.resize(N);
		vvvSLICLabels.resize(N);
		for(long n = nStart; n < nStop; n++)
		{
			IplImage img = vImgs[n].operator IplImage();
			IplImage* pimg = &img;
			//ExtLowLevelImgFeature::ExtColorFeatures_SLIC_BS0(pimg, vvvFeatures[n], vvvSLICLabels[n], nsp, pRefImg, vvvCorrs[n], imsize);
			ExtLowLevelImgFeature::ExtColorFeatures_SLIC_BS0(pimg, vvvFeatures[n], vvvSLICLabels[n], nsp, pRefImg, vvvCorrs[n], vvKSs[n], imsize, SLIC_is_given);
		}
		return S_OK;
	}	

	static long HTKGMMMAP_BS0_Fusion_With_Scores_SLIC(
		vector<GMM>& GMMs, 
		string& GMMPath, 
		string& szFEFolder, 
		vector<string>& vszImagePathNames, 
		vector< vector<int> >& vvRecogClasses, 
		vector<string>& vszLoadadFiles, 
		vector< vector< vector<cv::Point> > >& vvvSLICLabels, 
		int nsp, 
		vector< vector<int> >& vvxyrange, 
		cv::Mat& ref, 
		vector<GMM>& vFusionGMMs, 
		bool bFrameBatch, 
		int nAdaptMode, 
		int nStartSeq,
		int nStopSeq,
		vector< vector<double> >& vvScore0, 
		vector< vector<double> >& vvScore1
		)
	{
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		cv::Size imsize;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = nStartSeq - 1;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		hr = ExtColorFeaturesFromImagesStartStop_SLIC(vImgs, vvvFeatures, vvvSLICLabels, nsp, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs, imsize);
		// background frames with null vectors!
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;

		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		vector< vector< vector<double> > > vvvRecogScores;
		int nNumIter = GMMITER;
		if(nAdaptMode == 1)
			nNumIter = 0;
		if(bFrameBatch)
		{
			// Framewise adaptation
			for(long n = 0; n < vvvFeatures.size(); n++)
			{
				if(n >= nStartSeq && n < nStopSeq)
				{
					vector<int> vRecogClasses0;
					vector<int> vRecogClasses1;
					vector< vector< vector<double> > > vvvRecogVec0;
					vector< vector< vector<double> > > vvvRecogVec1;
					unordered_map<string, GMM> hNewGMMs;
					vector<GMM> newGMMs;

					vector<GMM> GMMsCopy;
					GMMsCopy = GMMs;

					for(int nIter = 1; nIter <= nNumIter; nIter++)
					{
						vvvRecogVec0.clear();
						hNewGMMs.clear();
						newGMMs.clear();
						TestAFrameFeatures_BS0(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
						GMMHash2VectBH(hNewGMMs, newGMMs);
						GMMsCopy.clear();
						GMMsCopy = newGMMs;
					}
					if(nNumIter < 1)
					{
						newGMMs = GMMs;
					}
					vector< vector<double> > vvScores;
					TestAFrameFeatures_BS0_Scores(vvvFeatures, newGMMs, vvScores, n, vvvCorrs);
					vvvRecogScores.push_back(vvScores);
				}
				else
				{
					vector< vector<double> > vvScores;
					vvvRecogScores.push_back(vvScores);
				}
			}

		}
		else
		{
			// Batch adaptation
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= nNumIter; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			if(nNumIter < 1)
			{
				newGMMs = GMMs;
			}
			AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs);
		}
		vvRecogClasses.resize(vvvFeatures.size());
		vvScore0.resize(vvvFeatures.size());
		vvScore1.resize(vvvFeatures.size());
		if(nAdaptMode == 0)
		{
			// MAP only
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvRecogClasses[i].resize(vvvRecogScores[i].size(), 0);
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						vvScore0[i][j] = vvvRecogScores[i][j][0];
						vvScore1[i][j] = vvvRecogScores[i][j][1];
						if(vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0] > 0.0)
						{
							vvRecogClasses[i][j] = 1;
						}
					}
				}
			}
		}
		else
		{
			// Score fusions
			vector< vector< vector<double> > > vvvScoreForFusions(vvvFeatures.size());
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvvScoreForFusions[i].resize(vvvRecogScores[i].size());
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						double LLR = vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0];
						double sLLR = 1.0/(1.0+exp(-0.5*LLR));
						double sCorr = 0.1140*vvvCorrs[i][j][0] + 0.5870 * vvvCorrs[i][j][1] + 0.2989*vvvCorrs[i][j][2];
						vvvScoreForFusions[i][j].push_back(sLLR);
						vvvScoreForFusions[i][j].push_back(sCorr);
					}
				}
			}
			vector< vector< vector<double> > > vvvRecogFusionScores;
			TestAllFeatures_Scores(vvvScoreForFusions, 
								vFusionGMMs, 
								vvRecogClasses,
								vvvRecogFusionScores
								);
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
				{
					vvScore0[i][j] = vvvRecogFusionScores[i][j][0];
					vvScore1[i][j] = vvvRecogFusionScores[i][j][1];
				}
			}
		}

		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static long HTKGMMMAP_BS0_Fusion_With_Scores_SLIC(
		vector<GMM>& GMMs, 
		string& GMMPath, 
		string& szFEFolder, 
		vector<string>& vszImagePathNames, 
		vector< vector<int> >& vvRecogClasses, 
		vector<string>& vszLoadadFiles, 
		vector< vector< vector<cv::Point> > >& vvvSLICLabels,
		int nsp, 
		vector< vector<int> >& vvxyrange, 
		cv::Mat& ref, 
		vector<GMM>& vFusionGMMs, 
		bool bFrameBatch, 
		int nAdaptMode, 
		int nStartSeq,
		int nStopSeq,
		vector< vector<double> >& vvScore0, 
		vector< vector<double> >& vvScore1,
		vector< vector< vector<double> > >& vvvScoreForFusions
		)
	{
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		cv::Size imsize;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = nStartSeq - 1;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		// Try to load SLIC data here
		vvvSLICLabels.resize(vszLoadadFiles.size());
		bool SLIC_is_given = true;
		for(size_t nimgidx = 0; nimgidx < vszLoadadFiles.size(); ++nimgidx)
		{
			string szFoldertmp;
			string szNametmp;
			string szDelimtmp;
			string szExttmp;
			string szSLICdirtmp;
			string szSLICdattmp;
			//string sztmp = "E:\\CollectedDataCopy\\Modifications\\CollectedDataCopy_blind_heq\\SLIC\\Train\\";
			GMMAdaptTool::ParseFileName(vszLoadadFiles[nimgidx], szFoldertmp, szNametmp, szExttmp, szDelimtmp);
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
			if(0 != szSLICdattmp.compare("")) {
				My_SLIC_Interface::LoadSLICDat(szSLICdattmp, vImgs[nimgidx].cols, vImgs[nimgidx].rows, vvvSLICLabels[nimgidx]); 
				if(0 == (int)vvvSLICLabels[nimgidx].size()) {
					SLIC_is_given = false;
				}
			}
		}
		//hr = ExtColorFeaturesFromImagesStartStop(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs);
		//hr = ExtColorFeaturesFromImagesStartStop_SLIC(vImgs, vvvFeatures, vvvSLICLabels, nsp, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, imsize);
		hr = ExtColorFeaturesFromImagesStartStop_SLIC(vImgs, vvvFeatures, vvvSLICLabels, nsp, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs, imsize, SLIC_is_given);
		// background frames with null vectors!
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		vector< vector< vector<double> > > vvvRecogScores;
		int nNumIter = GMMITER;
		if(nAdaptMode == 1)
			nNumIter = 0;
		if(bFrameBatch)
		{
			// Framewise adaptation
			for(long n = 0; n < vvvFeatures.size(); n++)
			{
				if(n >= nStartSeq && n < nStopSeq)
				{
					vector<int> vRecogClasses0;
					vector<int> vRecogClasses1;
					vector< vector< vector<double> > > vvvRecogVec0;
					vector< vector< vector<double> > > vvvRecogVec1;
					unordered_map<string, GMM> hNewGMMs;
					vector<GMM> newGMMs;

					vector<GMM> GMMsCopy;
					GMMsCopy = GMMs;

					for(int nIter = 1; nIter <= nNumIter; nIter++)
					{
						vvvRecogVec0.clear();
						hNewGMMs.clear();
						newGMMs.clear();
						TestAFrameFeatures_BS0(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
						GMMHash2VectBH(hNewGMMs, newGMMs);
						GMMsCopy.clear();
						GMMsCopy = newGMMs;
					}
					if(nNumIter < 1)
					{
						newGMMs = GMMs;
					}
					vector< vector<double> > vvScores;
					TestAFrameFeatures_BS0_Scores(vvvFeatures, newGMMs, vvScores, n, vvvCorrs);
					vvvRecogScores.push_back(vvScores);
				}
				else
				{
					vector< vector<double> > vvScores;
					vvvRecogScores.push_back(vvScores);
				}
			}

		}
		else
		{
			// Batch adaptation
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= nNumIter; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			if(nNumIter < 1)
			{
				newGMMs = GMMs;
			}
			AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs);
		}
		vvRecogClasses.resize(vvvFeatures.size());
		vvScore0.resize(vvvFeatures.size());
		vvScore1.resize(vvvFeatures.size());
		if(nAdaptMode == 0)
		{
			// MAP only
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvRecogClasses[i].resize(vvvRecogScores[i].size(), 0);
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						vvScore0[i][j] = vvvRecogScores[i][j][0];
						vvScore1[i][j] = vvvRecogScores[i][j][1];
						if(vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0] > 0.0)
						{
							vvRecogClasses[i][j] = 1;
						}
					}
				}
			}
		}
		else
		{
			// Score fusions
			//vector< vector< vector<double> > > vvvScoreForFusions(vvvFeatures.size());
			vvvScoreForFusions.resize(vvvFeatures.size());
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvvScoreForFusions[i].resize(vvvRecogScores[i].size());
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						double LLR = vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0];
						double sLLR = 1.0/(1.0+exp(-0.5*LLR));
						double sCorr = 0.1140*vvvCorrs[i][j][0] + 0.5870 * vvvCorrs[i][j][1] + 0.2989*vvvCorrs[i][j][2];
						vvvScoreForFusions[i][j].push_back(sLLR);
						vvvScoreForFusions[i][j].push_back(sCorr);
					}
				}
			}
			vector< vector< vector<double> > > vvvRecogFusionScores;
			TestAllFeatures_Scores(vvvScoreForFusions, 
								vFusionGMMs, 
								vvRecogClasses,
								vvvRecogFusionScores
								);
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
				{
					vvScore0[i][j] = vvvRecogFusionScores[i][j][0];
					vvScore1[i][j] = vvvRecogFusionScores[i][j][1];
				}
			}
		}
		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static void AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int nsp, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq, vector<cv::Mat>& vScore0, vector<cv::Mat>& vScore1)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector< vector<cv::Point> > > vvvSLICLabels;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		vector< vector<double> > vvScore0;
		vector< vector<double> > vvScore1;
		vector< vector< vector<double> > > vvvRecogFusionScores;
		//HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		HTKGMMMAP_BS0_Fusion_With_Scores_SLIC(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, vvvSLICLabels, nsp, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1);
		
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		vScore0.resize(N);
		vScore1.resize(N);
		
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0]);
			int yrange = (vvxyrange[n][1]);
			// Scores0
			vScore0[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Scores1
			vScore1[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Decision
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			int n_my_nn = 0;
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
				size_t n_sp_size_nn = vvvSLICLabels[n][nn].size();
				if(n_sp_size_nn > 0)
				{
					for(size_t mm = 0; mm < vvvSLICLabels[n][nn].size(); ++mm)
					{
						vScore0[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvScore0[n][nn];
						vScore1[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvScore1[n][nn];
						if(vvRecogClasses[n][n_my_nn] == nHandModelIndex)
							vMasks[n].at<uchar>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = 255;
						else
							vMasks[n].at<uchar>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = 0;
					}
					n_my_nn++;
				}				
			}
		}
	}

	static void AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int nsp, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq, vector<cv::Mat>& vScore0, vector<cv::Mat>& vScore1, vector<cv::Mat>& vLLR, vector<cv::Mat>& vCorr)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector< vector<cv::Point> > > vvvSLICLabels;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		vector< vector<double> > vvScore0;
		vector< vector<double> > vvScore1;
		vector< vector< vector<double> > > vvvScoreForFusions;
		vector< vector< vector<double> > > vvvRecogFusionScores;
		//HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		HTKGMMMAP_BS0_Fusion_With_Scores_SLIC(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, vvvSLICLabels, nsp, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1, vvvScoreForFusions);
		
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		vScore0.resize(N);
		vScore1.resize(N);
		vLLR.resize(N);
		vCorr.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0]);
			int yrange = (vvxyrange[n][1]);
			// Scores0
			vScore0[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Scores1
			vScore1[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Decision
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			vLLR[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			vCorr[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			int n_my_nn = 0;
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
				size_t n_sp_size_nn = vvvSLICLabels[n][nn].size();
				if(n_sp_size_nn > 0)
				{
					for(size_t mm = 0; mm < vvvSLICLabels[n][nn].size(); ++mm)
					{
						vScore0[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvScore0[n][nn];
						vScore1[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvScore1[n][nn];
						vLLR[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvvScoreForFusions[n][nn][0];
						vCorr[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvvScoreForFusions[n][nn][1];

						if(vvRecogClasses[n][n_my_nn] == nHandModelIndex)
							vMasks[n].at<uchar>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = 255;
						else
							vMasks[n].at<uchar>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = 0;
					}
					n_my_nn++;
				}				
			}
		}
	}

	static void AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int nsp, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq, vector<cv::Mat>& vScore0, vector<cv::Mat>& vScore1, vector<cv::Mat>& vLLR, vector<cv::Mat>& vCorr, vector<cv::Mat>& vKS)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector< vector<cv::Point> > > vvvSLICLabels;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		vector< vector<double> > vvScore0;
		vector< vector<double> > vvScore1;
		vector< vector< vector<double> > > vvvScoreForFusions;
		vector< vector< vector<double> > > vvvRecogFusionScores;
		//HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		HTKGMMMAP_BS0_Fusion_With_Scores_SLIC(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, vvvSLICLabels, nsp, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1, vvvScoreForFusions);
		
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		vScore0.resize(N);
		vScore1.resize(N);
		vLLR.resize(N);
		vCorr.resize(N);
		vKS.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0]);
			int yrange = (vvxyrange[n][1]);
			// Scores0
			vScore0[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Scores1
			vScore1[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Decision
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			vLLR[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			vCorr[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			vKS[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			int n_my_nn = 0;
			for(size_t nn = 0; nn < nResultsSize; nn++)
			{
				size_t n_sp_size_nn = vvvSLICLabels[n][nn].size();
				if(n_sp_size_nn > 0)
				{
					for(size_t mm = 0; mm < vvvSLICLabels[n][nn].size(); ++mm)
					{
						vScore0[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvScore0[n][nn];
						vScore1[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvScore1[n][nn];
						vLLR[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvvScoreForFusions[n][nn][0];
						vCorr[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvvScoreForFusions[n][nn][1];
						vKS[n].at<double>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = vvvScoreForFusions[n][nn][2];

						if(vvRecogClasses[n][n_my_nn] == nHandModelIndex)
							vMasks[n].at<uchar>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = 255;
						else
							vMasks[n].at<uchar>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = 0;
					}
					n_my_nn++;
				}				
			}
		}
	}


	//static HRESULT FindSLICCorrs(
	//	vector<string>& vszImagePathNames, 
	//	vector<string>& vszLoadadFiles, 
	//	vector< vector< vector<double> > >& vvvFeatures,
	//	vector< vector< vector<cv::Point> > >& vvvSLICLabels, 
	//	vector< vector< vector<double> > >& vvvCorrs,
	//	vector<cv::Mat>& vSLICmasks,
	//	cv::Mat& ref,
	//	int nStart,
	//	int nStop,
	//	int nsp
	//)
	static HRESULT FindSLICCorrs(
		vector<string>& vszImagePathNames, 
		vector<string>& vszLoadadFiles, 
		vector< vector< vector<double> > >& vvvFeatures,
		vector< vector< vector<cv::Point> > >& vvvSLICLabels, 
		vector< vector< vector<double> > >& vvvCorrs,
		vector< vector<double> >& vvKSs,
		vector<cv::Mat>& vSLICmasks,
		cv::Mat& ref,
		int nStart,
		int nStop,
		int nsp
	)
	{
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		// Try to load SLIC data here
		vvvSLICLabels.resize(vszLoadadFiles.size());
		bool SLIC_is_given = true;
		for(size_t nimgidx = 0; nimgidx < vszLoadadFiles.size(); ++nimgidx)
		{
			string szFoldertmp;
			string szNametmp;
			string szDelimtmp;
			string szExttmp;
			string szSLICdirtmp;
			string szSLICdattmp;
			//string sztmp = "E:\\CollectedDataCopy\\Modifications\\CollectedDataCopy_blind_heq\\SLIC\\Train\\";
			GMMAdaptTool::ParseFileName(vszLoadadFiles[nimgidx], szFoldertmp, szNametmp, szExttmp, szDelimtmp);
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
			if(0 != szSLICdattmp.compare("")) {
				My_SLIC_Interface::LoadSLICDat(szSLICdattmp, vImgs[nimgidx].cols, vImgs[nimgidx].rows, vvvSLICLabels[nimgidx]); 
				if(0 == (int)vvvSLICLabels[nimgidx].size()) {
					SLIC_is_given = false;
				}
			}
		}
		if(hr == S_OK)
		{
			size_t N_imgs = vszImagePathNames.size();
			vvvSLICLabels.resize(N_imgs);
			vvvCorrs.resize(N_imgs);
			vvKSs.resize(N_imgs);
			vvvFeatures.resize(N_imgs);
			vSLICmasks.resize(N_imgs);
			IplImage iplref = ref.operator IplImage();
			IplImage* pref = &iplref;			
			for(int n = nStart; n <= nStop; ++n)
			{
				cv::Size imsize;
				cv::Mat img = vImgs[n];
				vSLICmasks[n] = cv::Mat(img.size(), CV_8UC1);
				IplImage iplimg = img.operator IplImage();
				IplImage* pimg = &iplimg;
				if(!SLIC_is_given)
				{
					My_SLIC_Interface slic_interface;
					slic_interface.Eval_SLICO_Single_Frame(img, nsp, vvvSLICLabels[n]);
				}
				//int nsps = ExtLowLevelImgFeature::ExtColorFeatures_SLIC_BS0(pimg, vvvFeatures[n], vvvSLICLabels[n], nsp, pref, vvvCorrs[n], imsize);
				int nsps = ExtLowLevelImgFeature::ExtColorFeatures_SLIC_BS0(pimg, vvvFeatures[n], vvvSLICLabels[n], nsp, pref, vvvCorrs[n], vvKSs[n], imsize, SLIC_is_given);
				for(size_t nn = 0; nn < vvvCorrs[n].size(); ++nn)
				{
					double Corr = 0.1140 * vvvCorrs[n][nn][0] + 0.5870 * vvvCorrs[n][nn][1] + 0.2989 * vvvCorrs[n][nn][2];
					for(size_t mm = 0; mm < vvvSLICLabels[n][nn].size(); ++mm)
					{
						if(Corr < CCThreshold)
							vSLICmasks[n].at<uchar>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = 255;
						else
							vSLICmasks[n].at<uchar>(vvvSLICLabels[n][nn][mm].x, vvvSLICLabels[n][nn][mm].y) = 0;
					}
				}
			}
		}
		return hr;
	}

	//static long ExtColorFeaturesFromImagesStartStop_SLIC_MASK(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, int BLOCKSIZE, int nStart, int nStop, IplImage*& pRefImg, vector< vector< vector<double> > >& vvvCorrs, vector< unordered_map<int,int> >& vhValids, vector< cv::Mat >& vSLICMasks)
	static long ExtColorFeaturesFromImagesStartStop_SLIC_MASK(vector<cv::Mat>& vImgs, vector< vector< vector<double> > >& vvvFeatures, int BLOCKSIZE, int nStart, int nStop, IplImage*& pRefImg, vector< vector< vector<double> > >& vvvCorrs, vector< vector<double> >& vvKSs, vector< unordered_map<int,int> >& vhValids, vector< cv::Mat >& vSLICMasks)
	{
		long N = (long) vImgs.size();
		if(N <= 0)
		{
			cerr<<"No images for feature extraction"<<endl;
			return ERROR_INVALID_DATA;
		}
		vvvFeatures.resize(N);
		vvvCorrs.resize(N);
		vvKSs.resize(N);
		vhValids.resize(N);
		for(long n = nStart; n < nStop; n++)
		{
			IplImage img = vImgs[n].operator IplImage();
			IplImage* pimg = &img;
			//ExtLowLevelImgFeature::ExtColorFeatures_BS0_SLIC_MASK(pimg, vvvFeatures[n], BLOCKSIZE, pRefImg, vvvCorrs[n], vhValids[n], vSLICMasks[n]);
			ExtLowLevelImgFeature::ExtColorFeatures_BS0_SLIC_MASK(pimg, vvvFeatures[n], BLOCKSIZE, pRefImg, vvvCorrs[n], vvKSs[n], vhValids[n], vSLICMasks[n]);
			//ExtLowLevelImgFeature::ExtColorFeatures(pimg, vvvFeatures[n], BLOCKSIZE);
		}
		return S_OK;
	}
	
	static long HTKGMMMAP_BS0_Fusion_With_Scores_SLIC_MASK(
		vector<GMM>& GMMs, 
		string& GMMPath, 
		string& szFEFolder, 
		vector<string>& vszImagePathNames, 
		vector< vector<int> >& vvRecogClasses, 
		vector<string>& vszLoadadFiles, 
		int BLOCKSIZE, 
		vector< vector<int> >& vvxyrange, 
		cv::Mat& ref, 
		vector<GMM>& vFusionGMMs, 
		bool bFrameBatch, 
		int nAdaptMode, 
		int nStartSeq,
		int nStopSeq,
		vector< vector<double> >& vvScore0, 
		vector< vector<double> >& vvScore1,
		vector< cv::Mat >& vSLICMasks, 
		vector< unordered_map<int, int> >& vhValids, 
		vector< vector< vector<double> > >& vvvScoreForFusions
		)
	{
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = nStartSeq - 1;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop_SLIC_MASK(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vhValids, vSLICMasks);
		hr = ExtColorFeaturesFromImagesStartStop_SLIC_MASK(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs, vhValids, vSLICMasks);
		// background frames with null vectors!
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;
		//TestAllFeatures_BS0(vvvFeatures, GMMs, vvRecogClasses0, vvvRecogVec0, vvvCorrs);
	
		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		vector< vector< vector<double> > > vvvRecogScores;
		int nNumIter = GMMITER;
		if(nAdaptMode == 1)
			nNumIter = 0;
		if(bFrameBatch)
		{
			// Framewise adaptation
			for(long n = 0; n < vvvFeatures.size(); n++)
			{
				if(n >= nStartSeq && n < nStopSeq)
				{
					vector<int> vRecogClasses0;
					vector<int> vRecogClasses1;
					vector< vector< vector<double> > > vvvRecogVec0;
					vector< vector< vector<double> > > vvvRecogVec1;
					unordered_map<string, GMM> hNewGMMs;
					vector<GMM> newGMMs;

					vector<GMM> GMMsCopy;
					GMMsCopy = GMMs;

					for(int nIter = 1; nIter <= nNumIter; nIter++)
					{
						vvvRecogVec0.clear();
						hNewGMMs.clear();
						newGMMs.clear();
						TestAFrameFeatures_BS0(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
						GMMHash2VectBH(hNewGMMs, newGMMs);
						GMMsCopy.clear();
						GMMsCopy = newGMMs;
					}
					if(nNumIter < 1)
					{
						newGMMs = GMMs;
					}
					vector< vector<double> > vvScores;
					TestAFrameFeatures_BS0_Scores(vvvFeatures, newGMMs, vvScores, n, vvvCorrs);
					vvvRecogScores.push_back(vvScores);
				}
				else
				{
					vector< vector<double> > vvScores;
					vvvRecogScores.push_back(vvScores);
				}
			}

		}
		else
		{
			// Batch adaptation
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= nNumIter; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			if(nNumIter < 1)
			{
				newGMMs = GMMs;
			}
			AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs);
		}
		vvRecogClasses.resize(vvvFeatures.size());
		vvScore0.resize(vvvFeatures.size());
		vvScore1.resize(vvvFeatures.size());
		if(nAdaptMode == 0)
		{
			// MAP only
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvRecogClasses[i].resize(vvvRecogScores[i].size(), 0);
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						vvScore0[i][j] = vvvRecogScores[i][j][0];
						vvScore1[i][j] = vvvRecogScores[i][j][1];
						if(vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0] > 0.0)
						{
							vvRecogClasses[i][j] = 1;
						}
					}
				}
			}
		}
		else
		{
			// Score fusions
			vvvScoreForFusions.resize(vvvFeatures.size());
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvvScoreForFusions[i].resize(vvvRecogScores[i].size());
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						double LLR = vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0];
						double sLLR = 1.0/(1.0+exp(-0.5*LLR));
						double sCorr = 0.1140*vvvCorrs[i][j][0] + 0.5870 * vvvCorrs[i][j][1] + 0.2989*vvvCorrs[i][j][2];
						vvvScoreForFusions[i][j].push_back(sLLR);
						vvvScoreForFusions[i][j].push_back(sCorr);
					}
				}
			}
			vector< vector< vector<double> > > vvvRecogFusionScores;
			TestAllFeatures_Scores(vvvScoreForFusions, 
								vFusionGMMs, 
								vvRecogClasses,
								vvvRecogFusionScores
								);
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
				{
					vvScore0[i][j] = vvvRecogFusionScores[i][j][0];
					vvScore1[i][j] = vvvRecogFusionScores[i][j][1];
				}
			}
		}

		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}

	static long HTKGMMMAP_BS0_Fusion_With_Scores_SLIC_MASK(
		vector<GMM>& GMMs, 
		string& GMMPath, 
		string& szFEFolder, 
		vector<string>& vszImagePathNames, 
		vector< vector<int> >& vvRecogClasses, 
		vector<string>& vszLoadadFiles, 
		int BLOCKSIZE, 
		vector< vector<int> >& vvxyrange, 
		cv::Mat& ref, 
		vector<GMM>& vFusionGMMs, 
		bool bFrameBatch, 
		int nAdaptMode, 
		int nStartSeq,
		int nStopSeq,
		vector< vector<double> >& vvScore0, 
		vector< vector<double> >& vvScore1,
		vector< unordered_map<int, int> >& vhValids, 
		vector< cv::Mat >& vSLICMasks
		)
	{
		vector< vector< vector<double> > > vvvCorrs;
		vector< vector<double> > vvKSs;
		long N = (long) vszImagePathNames.size();
		long I = (long) GMMs.size();
		vector< vector< vector<double> > > vvvFeatures;
		if(N <= 0)
		{
			cerr<<"The image paths are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		if(I <= 0)
		{
			cerr<<"The GMMs are empty"<<endl;
			return ERROR_INVALID_DATA;
		}
		vector<cv::Mat> vImgs;
		HRESULT hr = LoadImagesFromPaths(vszImagePathNames, vImgs, vszLoadadFiles);
		long NN = (long)vImgs.size();
		if(NN <= 0 || hr != S_OK)
		{
			cerr<<"No images are loaded"<<endl;
			return ERROR_INVALID_DATA;
		}
		int nDetStopBkGnd = nStartSeq - 1;
		IplImage* pBkGnd = NULL;
		cv::Mat BkGnd;
		AccBackGroundbyIndex(vImgs, 0, nDetStopBkGnd, BkGnd);
		IplImage BkGndipl = BkGnd.operator IplImage();
		pBkGnd = &BkGndipl;
		//hr = ExtColorFeaturesFromImagesStartStop_SLIC_MASK(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vhValids, vSLICMasks);
		hr = ExtColorFeaturesFromImagesStartStop_SLIC_MASK(vImgs, vvvFeatures, BLOCKSIZE, nStartSeq, nStopSeq, pBkGnd, vvvCorrs, vvKSs, vhValids, vSLICMasks);		
		// background frames with null vectors!
		if(hr != S_OK)
		{
			cerr<<"Error while extracting color features"<<endl;
			return hr;
		}
		string szFolder, szDelim;
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());

		vector< vector<int> > vvRecogClasses0;

		vector< vector< vector<double> > > vvvRecogVec0;
		vector< vector< vector<double> > > vvvRecogVec1;

		string szSeqName, szCurrPath;
		FindCurrentSeqNamePath(vszLoadadFiles[0], szSeqName, szCurrPath);
		string szTmpName = szCurrPath + szSeqName;
		vector< vector< vector<double> > > vvvRecogScores;
		int nNumIter = GMMITER;
		if(nAdaptMode == 1)
			nNumIter = 0;
		if(bFrameBatch)
		{
			// Framewise adaptation
			for(long n = 0; n < vvvFeatures.size(); n++)
			{
				if(n >= nStartSeq && n < nStopSeq)
				{
					vector<int> vRecogClasses0;
					vector<int> vRecogClasses1;
					vector< vector< vector<double> > > vvvRecogVec0;
					vector< vector< vector<double> > > vvvRecogVec1;
					unordered_map<string, GMM> hNewGMMs;
					vector<GMM> newGMMs;

					vector<GMM> GMMsCopy;
					GMMsCopy = GMMs;

					for(int nIter = 1; nIter <= nNumIter; nIter++)
					{
						vvvRecogVec0.clear();
						hNewGMMs.clear();
						newGMMs.clear();
						TestAFrameFeatures_BS0(vvvFeatures, GMMsCopy, vRecogClasses0, vvvRecogVec0, n, vvvCorrs);
						Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, vszLoadadFiles[n], hNewGMMs, nIter);
						GMMHash2VectBH(hNewGMMs, newGMMs);
						GMMsCopy.clear();
						GMMsCopy = newGMMs;
					}
					if(nNumIter < 1)
					{
						newGMMs = GMMs;
					}
					vector< vector<double> > vvScores;
					TestAFrameFeatures_BS0_Scores(vvvFeatures, newGMMs, vvScores, n, vvvCorrs);
					vvvRecogScores.push_back(vvScores);
				}
				else
				{
					vector< vector<double> > vvScores;
					vvvRecogScores.push_back(vvScores);
				}
			}

		}
		else
		{
			// Batch adaptation
			vector<GMM> GMMsCopy;
			GMMsCopy = GMMs;
			unordered_map<string, GMM> hNewGMMs;
			vector<GMM> newGMMs;
			for(int nIter = 1; nIter <= nNumIter; nIter++)
			{
				vvvRecogVec0.clear();
				hNewGMMs.clear();
				newGMMs.clear();
				//TestAllFeatures_BS0_V2(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				TestAllFeatures_BS0(vvvFeatures, GMMsCopy, vvRecogClasses0, vvvRecogVec0, vvvCorrs) ;
				Do_Skin_HTKMAP(vvvRecogVec0, GMMPath, szFEFolder, szTmpName, hNewGMMs, nIter);
				GMMHash2VectBH(hNewGMMs, newGMMs);
				GMMsCopy.clear();
				GMMsCopy = newGMMs;
			}
			if(nNumIter < 1)
			{
				newGMMs = GMMs;
			}
			AllFeatures_BS0_Scores(vvvFeatures, newGMMs, vvvRecogScores, vvvCorrs);
		}
		vvRecogClasses.resize(vvvFeatures.size());
		vvScore0.resize(vvvFeatures.size());
		vvScore1.resize(vvvFeatures.size());
		if(nAdaptMode == 0)
		{
			// MAP only
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvRecogClasses[i].resize(vvvRecogScores[i].size(), 0);
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						vvScore0[i][j] = vvvRecogScores[i][j][0];
						vvScore1[i][j] = vvvRecogScores[i][j][1];
						if(vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0] > 0.0)
						{
							vvRecogClasses[i][j] = 1;
						}
					}
				}
			}
		}
		else
		{
			// Score fusions
			vector< vector< vector<double> > > vvvScoreForFusions(vvvFeatures.size());
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				if(i >= nStartSeq && i < nStopSeq)
				{
					vvvScoreForFusions[i].resize(vvvRecogScores[i].size());
					vvScore0[i].resize(vvvRecogScores[i].size(), 0.0);
					vvScore1[i].resize(vvvRecogScores[i].size(), 0.0);
					for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
					{
						double LLR = vvvRecogScores[i][j][1] - vvvRecogScores[i][j][0];
						double sLLR = 1.0/(1.0+exp(-0.5*LLR));
						double sCorr = 0.1140*vvvCorrs[i][j][0] + 0.5870 * vvvCorrs[i][j][1] + 0.2989*vvvCorrs[i][j][2];
						vvvScoreForFusions[i][j].push_back(sLLR);
						vvvScoreForFusions[i][j].push_back(sCorr);
					}
				}
			}
			vector< vector< vector<double> > > vvvRecogFusionScores;
			TestAllFeatures_Scores(vvvScoreForFusions, 
								vFusionGMMs, 
								vvRecogClasses,
								vvvRecogFusionScores
								);
			for(int i = 0; i < vvvFeatures.size(); ++i)
			{
				for(size_t j = 0; j < vvvRecogScores[i].size(); ++j)
				{
					vvScore0[i][j] = vvvRecogFusionScores[i][j][0];
					vvScore1[i][j] = vvvRecogFusionScores[i][j][1];
				}
			}
		}

		vvxyrange.resize(NN);
		for(long n = 0; n < NN; n++)
		{
			vvxyrange[n].push_back(vImgs[n].cols);
			vvxyrange[n].push_back(vImgs[n].rows);
		}
		return S_OK;
	}	
	
	static void AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC_MASK(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq, vector<cv::Mat>& vScore0, vector<cv::Mat>& vScore1, vector<cv::Mat>& vSLICMasks, vector<cv::Mat>& vLLR, vector<cv::Mat>& vCorr)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
	#ifdef _MSC_VER
		string mkdirFolder = "mkdir " + szFEFolder;
	#elif defined(__GNUC__)
		string mkdirFolder = "mkdir -p " + szFEFolder;
	#endif
		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		vector< vector<double> > vvScore0;
		vector< vector<double> > vvScore1;
		vector< vector< vector<double> > > vvvScoreForFusions;
		//HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		//HTKGMMMAP_BS0_Fusion_With_Scores(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1);
		vector< unordered_map<int, int> > vhValids;
		HTKGMMMAP_BS0_Fusion_With_Scores_SLIC_MASK(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1, vSLICMasks, vhValids, vvvScoreForFusions);
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		vScore0.resize(N);
		vScore1.resize(N);
		vLLR.resize(N);
		vCorr.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			// Scores0
			vScore0[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Scores1
			vScore1[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Decision
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			vLLR[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			vCorr[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			size_t nResultSizeTmp = (size_t)(xrange * yrange);
			for(size_t nn = 0; nn < nResultSizeTmp; nn++)
			{
					int x = nn / yrange;
					int y = nn - x * yrange;
					uchar* uptr = vMasks[n].ptr<uchar>(y);
					double* ptr0 = vScore0[n].ptr<double>(y);
					double* ptr1 = vScore1[n].ptr<double>(y);
					double* ptrLLR = vLLR[n].ptr<double>(y);
					double* ptrCorr = vCorr[n].ptr<double>(y);
					if(vhValids[n].find(nn) != vhValids[n].end())
					{
						int nntmp = vhValids[n][nn];
						ptr0[x] = vvScore0[n][nntmp];
						ptr1[x] = vvScore1[n][nntmp];
						//cout<<"["<<n<<","<<nn<<"]:"<<vvvScoreForFusions[n][nn][0]<<","<<vvvScoreForFusions[n][nn][1]<<endl;
						ptrLLR[x] = vvvScoreForFusions[n][nntmp][0];
						ptrCorr[x] = vvvScoreForFusions[n][nntmp][1];
						if(vvRecogClasses[n][nntmp] == nHandModelIndex)
						{
								uptr[x] = 255;
						}
						else
						{
								uptr[x] = 0;
						}

					}
					else
					{
						ptr0[x] = -LZERO;
						ptr1[x] = LZERO;
						//cout<<"["<<n<<","<<nn<<"]:"<<vvvScoreForFusions[n][nn][0]<<","<<vvvScoreForFusions[n][nn][1]<<endl;
						ptrLLR[x] = -LZERO;
						ptrCorr[x] = 1.0;
						uptr[x] = 0;
					}
			}
		}

	}

	static void AutoLabelHTKMAP_BS0_Fusion_with_Scores_SLIC_MASK(vector<string>& vszImagePathNames, vector<GMM>& vGMMs, vector<cv::Mat>& vMasks, int shrinratio, int& nHandModelIndex, string& GMMPath, cv::Mat& ref, vector<GMM>& vFusionGMMs, bool bFrameBatch, int nAdaptMode, int nStartSeq, int nStopSeq, vector<cv::Mat>& vScore0, vector<cv::Mat>& vScore1, vector<cv::Mat>& vSLICMasks)
	{
		string szFEFolder;
		vector< vector<int> > vvRecogClasses;
		vector<string> vszLoadadFiles;
		vector< vector<int> > vvxyrange;
		vector< vector< vector<double> > > vvvCorrs;
		string szFolder = "";
		string szDelim = "";
		FindCurrentFolderANDDelim(GMMPath, szFolder, szDelim);
		szFEFolder = szFolder + szDelim + "MyFE" + szDelim;
		//string mkdirFolder = "mkdir -p " + szFEFolder;
		#ifdef _MSC_VER
			string mkdirFolder = "mkdir " + szFEFolder;
		#elif defined(__GNUC__)
			string mkdirFolder = "mkdir -p " + szFEFolder;
		#endif

		//struct stat st;
		if(!file_exist(szFEFolder)) // if(stat(szFEFolder.c_str(),&st) != 0)
			system(mkdirFolder.c_str());
		vector< vector<double> > vvScore0;
		vector< vector<double> > vvScore1;
		vector< vector< vector<double> > > vvvRecogFusionScores;
		//HTKGMMMAP_BS0_Fusion(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq);
		vector< unordered_map<int, int> > vhValids;
		HTKGMMMAP_BS0_Fusion_With_Scores_SLIC_MASK(vGMMs, GMMPath, szFEFolder, vszImagePathNames, vvRecogClasses, vszLoadadFiles, shrinratio, vvxyrange, ref, vFusionGMMs, bFrameBatch, nAdaptMode, nStartSeq, nStopSeq, vvScore0, vvScore1, vhValids, vSLICMasks);
		
		size_t N = vvRecogClasses.size();
		vMasks.resize(N);
		vScore0.resize(N);
		vScore1.resize(N);
		for(size_t n = 0; n < N; n++)
		{
			size_t nResultsSize = vvRecogClasses[n].size();
			int xrange = (vvxyrange[n][0] / shrinratio);
			int yrange = (vvxyrange[n][1] / shrinratio);
			// Scores0
			vScore0[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Scores1
			vScore1[n] = cv::Mat(cv::Size(xrange, yrange), CV_64FC1);
			// Decision
			vMasks[n] = cv::Mat(cv::Size(xrange, yrange), CV_8UC1);
			size_t nResultSizeTmp = (size_t)(xrange * yrange);
			for(size_t nn = 0; nn < nResultSizeTmp; nn++)
			{
				int x = nn / yrange;
				int y = nn - x * yrange;
				uchar* uptr = vMasks[n].ptr<uchar>(y);
				double* ptr0 = vScore0[n].ptr<double>(y);
				double* ptr1 = vScore1[n].ptr<double>(y);
				if(vhValids[n].find(nn) != vhValids[n].end())
				{
					int nntmp = vhValids[n][nn];
					ptr0[x] = vvScore0[n][nntmp];
					ptr1[x] = vvScore1[n][nntmp];
					if(vvRecogClasses[n][nntmp] == nHandModelIndex)
					{
							uptr[x] = 255;
					}
					else
					{
							uptr[x] = 0;
					}

				}
				else
				{
					ptr0[x] = -LZERO;
					ptr1[x] = LZERO;
					uptr[x] = 0;
				}
			}
		}

	}
	

	static void GuessSLICPath(string& szPathName, string& szSLICPathName)
	{
		string szFoldertmp;
		string szNametmp;
		string szDelimtmp;
		string szExttmp;
		string szSLICdirtmp;
		string szSLICdattmp;
		GMMAdaptTool::ParseFileName(szPathName, szFoldertmp, szNametmp, szExttmp, szDelimtmp);
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


};
