//---------------------------------------------------------------------------

#pragma hdrstop
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <dos.h>
#include <conio.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <sstream>
#include <dir.h>
#include <vector>

#include "UtilGraph.h"
#include "TSD.h"

#define Z_PREFIX
#include "GZlib\zlib.h"

const double LenStripX=10.; //26.
const double LenStepX=10.1; //26.4
const double LenStripY=1000.; //3460.
const double LenStripZ=7.; //

//������� ����������� ������������
int UseNBM=1;
int Nstrip=32;
std::vector<std::vector<std::vector<unsigned long> > > MatrCrossHitStripPlane;
std::vector<unsigned long> nMatrCrossHitStripPlane;
//������� ����������� ������������ �� LAM
std::vector<std::vector<std::vector<unsigned long> > > MatrCrossLamStripPlane;
std::vector<unsigned long> nMatrCrossLamStripPlane;
//������� ����������� ������������ ��� ������ �� LAM
std::vector<std::vector<std::vector<unsigned long> > > MatrCrossHitLamStripPlane;
std::vector<unsigned long> nMatrCrossHitLamStripPlane;


//��� �������� ������� ������������ �������
short ChanStrip[NCHANADC];//������� ������������ ����� ��� -> ����� ������
int idChanStrip=0; //������� ������������� ������� ������������ ChanStrip

int ChanProc=-1; //�������������� ����� ������������: -1 ��� ��������� ������
int porog=10; //����� � �������

static unsigned long nlenbuf=1000; //������ ������ ������ � TD
T_TSD TD; //����������������� ������ ������ ������� ��� ������ �� ����

PAR_OBR Par; //��������� ���������

const int LenBufTD=1000;
int nBufTD=0; //������� ������� � ������� BufTD
T_TSD BufTD[LenBufTD]; //������ �������
char FileData[256]; //��� ����� ��� ������ �������
int iBufTD=0;

static std::string DIRALL="STATALL"; //��� ���������� ��� ������ � ������� �� ���� �����
static char NameDirStart[256]; //������ ��� ���������� ������� ���������
static char NameDirWrite[256]; //������ ��� ���������� ���� ������������
static char NameDirSrc[256]; //������ ��� �������� ���������� c ������� �����
static char Obrabdir[256]; //�������������� ���������� ������ ����
static char NameRun[256]; //��� ���������� ����
static char file_read[256];
static char buf[512];
static char nameout[256];
static short currentrun; //����� �������� ������

double APMatrStrip[NCHANADC][NCHANADC]; //������� ���� �������� ���� ������
double sAPMatrStrip[NCHANADC][NCHANADC]; //������� ���� ��������� �������� ���� ������
unsigned long NPMatrStrip[NCHANADC][NCHANADC]; //������� ����� �������� ���� ������
unsigned long LamMatrStrip[NCHANADC][NCHANADC]; //������� ����� LAM
double AtotalStrip[NCHANADC]; //������� ���� �������� �� ������ ����� ����� ���� �������� ����
double sAtotalStrip[NCHANADC]; //������� ���� ��������� �������� �� ������ ����� ����� ���� �������� ����
double NtotalStrip[NCHANADC]; //������� ����� ������� ����� ����� ���� �������� ����
#define NLEN 50 //3460/(34.6*2)
double MatrLenAStrip[NLEN][NCHANADC]; //����� �������� �� ���������� �� ���
double MatrLenAPStrip[NLEN][NCHANADC]; //����� �������� � ������� �� ���������� �� ���
double MatrLenSAStrip[NLEN][NCHANADC]; //����� ��������� �������� �� ���������� �� ���
double MatrLenSAPStrip[NLEN][NCHANADC]; //����� ��������� �������� � ������� �� ���������� �� ���
unsigned long MatrLenNStrip[NLEN][NCHANADC]; //����������� ����� ������� �� ���������� �� ���
unsigned long MatrLenNPStrip[NLEN][NCHANADC]; //����������� ����� ������� � ������� �� ���������
unsigned long MatrLenLamStrip[NLEN][NCHANADC]; //����������� ����� ������� � LAM
									   // �� ���������� �� ���
unsigned long MNLamStrip[NCHANADC]; //����������� ����� ������� � LAM

const int NplMatrix=5;
/*
#define NM_X 181
#define NM_Y 81
#define D_X 5.
#define D_Y 20.
*/
#define NM_X 39
#define NM_Y 55
#define D_X 5.
#define D_Y 20.
double hMatrix[NplMatrix]={0,5,10,15,20};
double MatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long nMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long LnMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long noBMnMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long tMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long ind_nMatrixXY[NplMatrix][NM_Y][NM_X][NCHANADC]; //������� ����������� ������ �������������� ����������
double RMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long RnMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long RLnMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long RnoBMnMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long RtMatrixXY[NplMatrix][NM_Y][NM_X]; //������� ����������� ������ �������������� ����������
                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ���������
unsigned long Rind_nMatrixXY[NplMatrix][NM_Y][NM_X][NCHANADC]; //������� ����������� ������ �������������� ����������

                                     //�� ������ 0, 5, 10, 15,20 �� �� ��������� ��������� ��� ������ ������� ���
double AmpDet[9]; //������� �������� ����������� ������������� ��� ������ ���� ����� � ��������� (��� 10 ��)
double SigDet[9]; //����� ����������� ������������� ��� ������ ���� ����� � ���������(��� 10 ��)
double QDet[9]; //������� �������� ��������� ������������� ��� ������ ���� ����� � ��������� (��� 10 ��)
double SigQDet[9]; //����� ��������� ������������� ��� ������ ���� ����� � ���������(��� 10 ��)
double t0Det[9]; //������� �������� ������� ������� ������������� ��� ������ ���� ����� � ��������� (��� 10 ��)
double Sigt0Det[9]; //����� ������� ������� ������������� ��� ������ ���� ����� � ���������(��� 10 ��)
double TDet[9]; //������� �������� ������ �������� ������������� ��� ������ ���� ����� � ��������� (��� 10 ��)
double SigTDet[9]; //����� ������ �������� ������������� ��� ������ ���� ����� � ���������(��� 10 ��)
unsigned long nDet[9];//���������� ������� � ��������
unsigned long nDetAll[9];//���������� ������� ����������� ����� ��������
double SredLen[9]; //������� ����� ����� � ���������
double SigLen[9]; //����� ����� ����� � ���������
double SredLenAll[9]; //������� ����� ����� � ���������
double SigLenAll[9]; //����� ����� ����� � ���������
long LamCurrStrip=0;

double NTrigURG=0.; //����� ����� ��������� ������
double NEventURG=0.; //����� ����� ������� ������ �� ��� ������� ������� � ������
double NDataBM=0.; //����� ����� ������� � ��������� ��
double NDataBMBad=0.; //��������� ������ �� ���� ����� 15 ��� 16
double NTrackAll=0.; //����� ����� ������ ��������� ����� ������
double NAmpPorogAll=0.; //����� ����� ������� � ������ �������, ����� �� ����� ������ ���� ���������� ������/
//��� Theta 0-15
double NTrackStrip=0.; //����� ����� ������ ��������� ����� ����� ������ ������ �������� ������
double NTrackStripFull=0.; //����� ����� ����������� ������ ��������� ����� ����� ������ ������ �������� ������
double NAmpPorogStrip=0.; //����� ����� ������� � ������� ����� ����� ������ ������, ����� �� ����� ������ ���� ���������� ������/
double NLamStrip=0.; //����� ����� ������� � ������� ����� ����� ������ ������, ����� �� ����� ������ ��� LAM/
double NLamCurrStrip=0.; //���������� LAM �� ������� ����� ����� (+-5 ��) ������� ������ ���� (0-15 �� ������)/
double NAmpCurrStrip=0.;
double NLamAmpPorogCurrStrip=0.; //���������� LAM �� ����������� ������� ����� ����� (+-5 ��) ������� ������ ���� (0-15 �� ������)/
double NTrackBox=0.; //����� ����� ������ ��������� ����� ������� �������� ������
double NAmpPorogBox=0.; //����� ����� ������� � �������, ����� �� ����� ������ ���� ���������� ������/
double NDataBox=0.; //����� ����� ������� c �������, ����� ���� ������ ��
double NLamBox=0.; //����� ����� ������� � �������, ����� �� ����� ������ ��� LAM
double NLamAmpPorogBox=0.;	//����� ����� LAM ��� ������������ �� �� ������� � ������� (0-15 ��. �� ������) ���������� ����� ��
// ��������� ������ �� ���������� DirRead -> ��������� � ���������� DirWrite
short Obrab_data(char *DirRead,char *DirWrite);

double BMSmesX=-8.,BMSmesY=0.,BMSmesZ=0., BMSmesD=0.; //�������� ��������� ����� ������������ ��
//---------------------------------------------------------------------------
//���������� ����� ����������� ������ �������������� �����
//���������� ����� �������� ����� ������� ��������������� ������
static short SetPointBox(double* p0,double* p1, TPoint3D boxmin,TPoint3D boxmax)
{
// ��������� ���������� �������������� �������� ����� ������
//     ��� ���������� ����� ����������� ������
static double ABCD[6][4]=  {  0.,1., 0., 0.,  //Y0=1075.8-275
			 0., 1., 0., 420.,   //Y1=Y0+420-275.(1495.8)
			 1., 0., 0.,  0.,   //X0=2372.4
			 1., 0., 0., 3108.75,   //X1=X0+3119.55
			 0., 0., 1.,  0.,   //Z0=5307.7
			 0., 0., 1., 2699.25    //Z1=Z0-2715.95
			  };
static double Rectmin[3]={ 0.,0.,0.};
static double Rectmax[3]={ 3108.75,420.,2699.25};
	double fx,fa,pt[2][3],t[6],a[3],p[3];
	short i,j,k,l;
        ABCD[0][3]=boxmin.y;
		ABCD[1][3]=boxmax.y;
        ABCD[2][3]=boxmin.x;
        ABCD[3][3]=boxmax.x;
        ABCD[4][3]=boxmin.z;
		ABCD[5][3]=boxmax.z;
        Rectmin[0]=boxmin.x; Rectmin[1]=boxmin.y; Rectmin[2]=boxmin.z;
        Rectmax[0]=boxmax.x; Rectmax[1]=boxmax.y; Rectmax[2]=boxmax.z;
	k=0;
	for(i = 0; i < 3; i++)
	{
		a[i] = p1[i] - p0[i];
	}

	for(i = 0; i < 6; i++)
	{
	  fa=ABCD[i][0]*a[0]+ABCD[i][1]*a[1]+ABCD[i][2]*a[2] ;
	  if(fa == 0.)continue;
	  fx = ABCD[i][0]*p0[0]+ABCD[i][1]*p0[1]+ABCD[i][2]*p0[2]-ABCD[i][3] ;
	  t[k]=-fx/fa;
	  for(j = 0; j < 3; j++)
	  {
		 p[j]=p0[j]+a[j]*t[k];
		 if(p[j] < Rectmin[j]-0.01)goto m1;
		 if(p[j] > Rectmax[j]+0.01)goto m1;
	  }
	  for(j = 0; j < 3; j++)pt[k][j]=p[j];
	  k++ ;
	  if(k > 1)
	  {  l=0;
		  for(j=0; j < 3; j++)
		  {
			 if( fabs(pt[0][j]-pt[1][j]) < 1.)l++ ;
		  }
		  if(l < 3)break;
		  k-- ;
	  }
m1:	continue;
	}
	if(k < 2)return 0;
	fa=0;
	for(i = 0; i < 3; i++)fa += (pt[1][i] - pt[0][i])*a[i];
	if(fa > 0.)
	{
	for(i = 0; i < 3; i++)p0[i] = pt[0][i];
	for(i = 0; i < 3; i++)p1[i] = pt[1][i];
	}
	else
	{
	for(i = 0; i < 3; i++)p0[i] = pt[1][i];
	for(i = 0; i < 3; i++)p1[i] = pt[0][i];
	}
	return k;
}
//---------------------------------------------------------------------------
//������� ����� � �����������
TPoint3D AxisTurn(const TPoint3D Point, const double Angle);
TPoint3D AxisTurn(const TPoint3D Point, const double Angle)
{
    TPoint3D NewPoint;
    NewPoint.x=Point.x*cos(Angle)+Point.y*sin(Angle);
	NewPoint.y=-Point.x*sin(Angle)+Point.y*cos(Angle);
	NewPoint.z=Point.z;
    return NewPoint;
}
void ReadBufTD(char *namefile)
{
//  char buf[80];
/*
  int handle=open(namefile,O_RDONLY|O_BINARY, S_IREAD);
  if(handle==-1)return;
  nBufTD=filelength(handle)/sizeof(T_TSD);
  read(handle,BufTD, nBufTD*sizeof(T_TSD));
  close(handle);
*/
  gzFile gzf;
  gzf=gzopen(namefile,"rb");
  int i=0;
  while(1)
  {
    if(gzread(gzf,&BufTD[i],sizeof(T_TSD))!=sizeof(T_TSD))break;
    i++;
  }
  gzclose(gzf);
  nBufTD=i;
  iBufTD=0;
}

//����� ������ ��������� ��������
TPoint3D Pdet0,Ce1;
//TPoint3D ChtMin,ChtMax; //������������� ���� �������������� ������� �������������
TPoint3D ChtMin[NCHANADC],ChtMax[NCHANADC]; //������������� ���� �������������� ������� ������
TPoint3D BoxMin,BoxMax; //���� �������������� ���� ������� ������
//������� ������������ ����� ������ - ��������� ����������
// yx - ��� "Y" , yz-��� "X"
struct TStripCoord
{
  double x;
  double y;
  double z;
};

extern short SetConfig(char *namecnf);

extern struct CONFIG_DECOR *CNF; // ��������� ������������ ����� ������
#define NCNTR 2  // ���-�� ������������ � ����� ���
#define NCHAN 4  // ���-�� ������� � ����� �����������
#define NGROUP8 76 // ����������� ��������� ���-�� ����� �� 8 ���
#define NGROUP16 38  // � ������ ������������ 38 ����� �� 16 ���
#define NGROUP64 10 // ����������� ��������� ���-�� ����� �� 64 ����
#define NPOINTMAX 320 // ����������� ��������� ����� � X-�� ��� � Y-��
#define MAXXSTRIP 320
#define MAXYSTRIP 288
#define NPOINTMAXX 320 // ����������� ��������� ����� � X-�� ��� � Y-��
#define NPOINTMAXY 288 // ����������� ��������� ����� � X-�� ��� � Y-��

extern short StripXmin,StripXmax,StripYmin,StripYmax,PLmin,PLmax;

//struct TStripCoord oneTablXY[8][NPOINTMAX];//[Sm][IPL][Strip]
extern struct TStripCoord *TablXY[16];
// ������� ����� ������� �� ��������� ���������
// ���������� ������, ��� �� ��� ������ ������� ��������.
extern double Pdcr7[16][3]; //[Sm][x,y,z], ��������� � SetNewCoordDcr()
//������ ����������� ��� X � ����������� ������� ��������� �����������
extern double VdcrX[16][3]; //[Sm][x,y,z], ��������� � SetNewCoordDcr()
//������ ����������� ��� Y � ����������� ������� ��������� �����������
extern double VdcrY[16][3]; //[Sm][x,y,z], ��������� � SetNewCoordDcr()
//������ ����������� ��� Z � ����������� ������� ��������� �����������
extern double VdcrZ[16][3]; //[Sm][x,y,z], ��������� � SetNewCoordDcr()

long SMconfig[4]={9,-1,-1,-1}; //����� �㯥ମ�㫥� ���������� � PCI ���⠬
long SMwork[4]={9,-1,-1,-1}; //����� �㯥ମ�㫥� ���������� � PCI ���⠬

//������� ������������ ������ ������ ���(c 1) ������ ��� (c 0)
/*
static int ADCtoPMT[NCHANADC]=
{  1,60, 2,59, 3,58, 4,57,  //1-8
   9,52,10,51,11,50,12,49,  //9-16
  17,44,18,43,19,42,20,41,  //17-24
  25,36,26,35,27,34,28,33,  //25-32
  32,37,31,38,30,39,29,40,  //33-40
  24,45,23,46,22,47,21,48,  //41-48
  16,53,15,54,14,55,13,56,  //49-56
   8,61, 7,62, 6,63, 5,64}; //57-64
*/
/*
static int ADCtoPMT[NCHANADC]=
{  1,2, 3,4, 5,6, 7,8,  //1-8
   9,10,11,12,13,14,15,16,  //9-16
  17,18,19,20,21,22,23,24,  //17-24
  25,26,27,28,29,30,31,32  //25-32
};
*/
static int ADCtoPMT[NCHANADC]=
{  32,31, 30,29, 28,27, 26,25,  //1-8
   24,23,22,21,20,19,18,17,  //9-16
  16,15,14,13,12,11,10,9,  //17-24
  8,7,6,5,4,3,2,1  //25-32
};
int PMTtoADC[NCHANADC]; //�������� � ADCtoPMT �������

#pragma argsused
int main(int argc, char* argv[])
{
  int done,donedir,handle,Nrun;
  int i,j,n;

  struct ffblk ffile; // ��� ������ ������
  struct ffblk fdir; // ��� ������ ��� ����������
//������� ��� �������� ���������� � ������
   if(argc < 2)
   {
m0:  printf("Usage: ObrabTMD.exe name_dir chan porog\r\nExample:\r\n"
			"  ObrabTMD.exe f:\\ 10 10\r\n   ObrabTMD.exe G:\\ChDet\\ -1 10 BM_1.txt");
	 printf("\r\nPress SPACE for exit...");
	 getch();
     exit(1);
   }
   if(argc > 2)
   {
    ChanProc=atoi(argv[2]); //�������������� �����
    ChanProc -=1;
   }
   if(ChanProc > NCHANADC-1)ChanProc=-1;
   if(ChanProc < -1)ChanProc=-1;
   if(argc > 3)porog=atoi(argv[3]); //����� � �������
   if(argc > 4)
   {
     std::ifstream source;
     source.open(argv[4]);
     if(!source)
	 {
	   printf("Error open ChanStrip file: %s",argv[4]); exit(1);
     }
     std::string srcline;
     std::istringstream iss;
	 //�����������
	 std::getline(source,srcline,'\n');
     std::getline(source,srcline,'\n');
     //������� ����������� ����� ��� -> ����� ������
     //����� �����
      for(int i=0; i < NCHANADC; i++)
      {
        std::getline(source,srcline,'\n');
        iss.clear();
		iss.str(srcline);
        int nomchan=0;
        int val=0;
        iss >> nomchan >> val;
        nomchan--; val--;
        if(nomchan < 0 || nomchan > NCHANADC-1)continue;
        ChanStrip[nomchan] = val;
      }
	 source.close();
     idChanStrip=1;
   }
   std::string srcline,svar;
   std::ifstream src;
   src.open("SMESBM.");
   if(src)
   {
     std::getline(src,srcline,'\n'); //������ �����������
     src >> svar; src >> svar; //DX
     if(svar.length() > 0)BMSmesX=atof(svar.c_str());
     src >> svar; src >> svar; //DY
     if(svar.length() > 0)BMSmesY=atof(svar.c_str());
     src >> svar; src >> svar; //DZ
     if(svar.length() > 0)BMSmesZ=atof(svar.c_str());
     src >> svar; src >> svar; //DD
	 if(svar.length() > 0)BMSmesD=atof(svar.c_str());
     src.close();
   }
   else
   {
     std::ofstream out;
     out.open("SMESBM");
     out << "�������� ��� ��������� ��������� ObrabTMD.exe\n";
     out << "DX\t"<<BMSmesX << "\n";
     out << "DY\t"<<BMSmesY << "\n";
     out << "DZ\t"<<BMSmesZ << "\n";
     out << "DD\t"<<BMSmesD << "\n";
     out.close();
   }
   strcpy(NameDirSrc,argv[1]);
   if(NameDirSrc[strlen(NameDirSrc)-1]!='\\')goto m0;
  //������� ��� ���������� ������ ��������� ���������
  strcpy(NameDirStart,"X:\\");
  NameDirStart[0]='A'+getdisk();
  i=getcurdir(0, NameDirStart+3);
  if(i){printf("Error getcurdir\n"); exit(1); }
  //�������� ���������� ��� �����������
  sprintf(buf,"STAT%02ld",porog);
  if(ChanProc >=0)
  {
    sprintf(buf+strlen(buf),"_%02ld",ChanProc+1);
  }
  else
  {
	strcat(buf,"ALL");
  }
  DIRALL=buf;
  done=findfirst(DIRALL.c_str(),&fdir,FA_DIREC); //���� �� ����� ��� ?
  findclose(&fdir);
  if(done) //����� ���
  {
    mkdir(DIRALL.c_str()); //������
    done=findfirst(DIRALL.c_str(),&fdir,FA_DIREC); //���������
	findclose(&fdir);
    if(done) // ���������� �� �������
    {
      printf("No make subdir %s",DIRALL.c_str()); exit(1);
    }
  }
         //���������� ��� ���������� ��� ������ ���� ������
         strcpy(NameDirWrite,NameDirStart);
         strcat(NameDirWrite,"\\");
         strcat(NameDirWrite,DIRALL.c_str());
         chdir(NameDirWrite);

   memset(APMatrStrip,0,sizeof(APMatrStrip)); //������� ���� �������� ���� ������
   memset(sAPMatrStrip,0,sizeof(sAPMatrStrip)); //������� ���� �������� ���� ������
   memset(NPMatrStrip,0,sizeof(NPMatrStrip)); //������� ����� �������� ���� ������
   memset(LamMatrStrip,0,sizeof(LamMatrStrip)); //������� ����� LAM
   memset(AtotalStrip,0,sizeof(AtotalStrip)); //������� ���� �������� �� ������ ����� ����� ���� �������� ����
   memset(sAtotalStrip,0,sizeof(sAtotalStrip)); //������� ���� �������� �� ������ ����� ����� ���� �������� ����
   memset(NtotalStrip,0,sizeof(NtotalStrip)); //������� ����� ������� ����� ����� ���� �������� ����
//#define NLEN 50 //3460/(34.6*2)
   memset(MatrLenAStrip,0,sizeof(MatrLenAStrip)); //����� �������� �� ���������� �� ���
   memset(MatrLenAPStrip,0,sizeof(MatrLenAPStrip)); //����� �������� �� ���������� �� ���
   memset(MatrLenSAStrip,0,sizeof(MatrLenSAStrip)); //����� �������� �� ���������� �� ���
   memset(MatrLenSAPStrip,0,sizeof(MatrLenSAPStrip)); //����� �������� �� ���������� �� ���
   memset(MatrLenNStrip,0,sizeof(MatrLenNStrip)); //����������� ����� ������� �� ���������� �� ���
   memset(MatrLenNPStrip,0,sizeof(MatrLenNPStrip)); //����������� ����� ������� � ������� �� ���������
   memset(MatrLenLamStrip,0,sizeof(MatrLenLamStrip)); //����������� ����� ������� � LAM
									   // �� ���������� �� ���
   memset(MNLamStrip,0,sizeof(MNLamStrip)); //����������� ����� ������� � LAM


   memset(MatrixXY,0,sizeof(MatrixXY));
   memset(nMatrixXY,0,sizeof(nMatrixXY));
   memset(LnMatrixXY,0,sizeof(LnMatrixXY));
   memset(noBMnMatrixXY,0,sizeof(noBMnMatrixXY));
   memset(tMatrixXY,0,sizeof(tMatrixXY));
   memset(ind_nMatrixXY,0,sizeof(ind_nMatrixXY));
   memset(RMatrixXY,0,sizeof(RMatrixXY));
   memset(RnMatrixXY,0,sizeof(RnMatrixXY));
   memset(RLnMatrixXY,0,sizeof(RLnMatrixXY));
   memset(RnoBMnMatrixXY,0,sizeof(RnoBMnMatrixXY));
   memset(RtMatrixXY,0,sizeof(RtMatrixXY));
   memset(Rind_nMatrixXY,0,sizeof(Rind_nMatrixXY));
   memset(AmpDet,0,sizeof(AmpDet));
   memset(SigDet,0,sizeof(SigDet));
   memset(QDet,0,sizeof(QDet));
   memset(SigQDet,0,sizeof(SigQDet));
   memset(t0Det,0,sizeof(t0Det));
   memset(Sigt0Det,0,sizeof(Sigt0Det));
   memset(TDet,0,sizeof(TDet));
   memset(SigTDet,0,sizeof(SigTDet));
   memset(SredLen,0,sizeof(SredLen));
   memset(SigLen,0,sizeof(SigLen));
   memset(SredLenAll,0,sizeof(SredLenAll));
   memset(SigLenAll,0,sizeof(SigLenAll));
   memset(nDet,0,sizeof(nDet));
   memset(nDetAll,0,sizeof(nDetAll));

 //������� ����������� ������������
 nMatrCrossHitStripPlane.resize(UseNBM);
 memset(&nMatrCrossHitStripPlane[0],0,sizeof(nMatrCrossHitStripPlane[0])*UseNBM);
 MatrCrossHitStripPlane.resize(UseNBM);
   for(int iBM=0; iBM < UseNBM; iBM++)
   {
    MatrCrossHitStripPlane[iBM].resize(Nstrip);
    for(int iStrip=0; iStrip < Nstrip; iStrip++)
    {
      MatrCrossHitStripPlane[iBM][iStrip].resize(Nstrip);
      memset(&MatrCrossHitStripPlane[iBM][iStrip][0],0,sizeof(MatrCrossHitStripPlane[iBM][iStrip][0])*Nstrip);
    }
   }
 //������� ����������� ������������ �� LAM
 nMatrCrossLamStripPlane.resize(UseNBM);
 memset(&nMatrCrossLamStripPlane[0],0,sizeof(nMatrCrossLamStripPlane[0])*UseNBM);
 MatrCrossLamStripPlane.resize(UseNBM);
   for(int iBM=0; iBM < UseNBM; iBM++)
   {
    MatrCrossLamStripPlane[iBM].resize(Nstrip);
    for(int iStrip=0; iStrip < Nstrip; iStrip++)
    {
      MatrCrossLamStripPlane[iBM][iStrip].resize(Nstrip);
      memset(&MatrCrossLamStripPlane[iBM][iStrip][0],0,sizeof(MatrCrossLamStripPlane[iBM][iStrip][0])*Nstrip);
    }
   }
 //������� ����������� ������������ ��� ������ �� LAM
 nMatrCrossHitLamStripPlane.resize(UseNBM);
 memset(&nMatrCrossHitLamStripPlane[0],0,sizeof(nMatrCrossHitLamStripPlane[0])*UseNBM);
 MatrCrossHitLamStripPlane.resize(UseNBM);
   for(int iBM=0; iBM < UseNBM; iBM++)
   {
    MatrCrossHitLamStripPlane[iBM].resize(Nstrip);
    for(int iStrip=0; iStrip < Nstrip; iStrip++)
    {
      MatrCrossHitLamStripPlane[iBM][iStrip].resize(Nstrip);
      memset(&MatrCrossHitLamStripPlane[iBM][iStrip][0],0,sizeof(MatrCrossHitLamStripPlane[iBM][iStrip][0])*Nstrip);
    }
   }


   Par.inverse=1;
   Par.strob[0]=100;
   Par.strob[1]=600;
   Par.R=50.;
   Par.porogImpulse[0]=-1000;
   Par.porogImpulse[1]=porog;
   Par.porogPDS[0]=0; //�� ������������
   Par.porogPDS[1]=0; //�� ������������
   for(int chan=0; chan < NCHANADC; chan++)
   {
    int ipmt=ADCtoPMT[chan]-1;
    PMTtoADC[ipmt]=chan;
   }

  StripXmin=0; StripXmax=319;
  StripYmin=6; StripYmax=282;
  PLmin=0; PLmax=7;
  strcpy(buf,NameDirStart);
  strcat(buf,"\\Decor.cnf");
  SetConfig(buf);
  //��������� ���������� �������
 double xcs[NCHANADC];
  for(int is=0; is < NCHANADC; is++)
 {
  if(is < 16)
   xcs[is]=LenStepX*is;
  else
   xcs[is]=LenStepX*(is-16)-5.;
  //������� ������� ������ 23x3000x7
  double lx=LenStripX, ly=LenStripY, lz=LenStripZ;
  Ce1.x=xcs[is];
  Ce1.y=0;//0;//Ce0.y;
  if(is < 16)
   Ce1.z=0;//Ce0.z;
  else
   Ce1.z=7.;//Ce0.z;
  //��������� ������� ��������������� ��������
  TPoint3D p0,p1,p2,p3,p4;
  p0.x=Ce1.x;
  p0.y=Ce1.y;
  p0.z=Ce1.z;
  p1=p0;
  p1.z += lz;
  p2=p1;
  p2.y += ly;
  p3=p2;
  p3.x += lx;
  p4=p3;
  p4.y -= ly;
  ChtMin[is]=p0;
  ChtMax[is]=p3;
 }
  double lx=LenStripX, ly=LenStripY, lz=LenStripZ;
  //�������� ����� � ����� ��� ������
  Ce1.x=NCHANADC/4.*LenStepX;
  Ce1.y=ly/2.;
  Ce1.z=5.;
//  BoxMin=TPoint3D(-10.,-5.,-5.);
//  BoxMax=TPoint3D(LenStepX*NCHANADC/2+5.,ly+5.,2.*lz+5.);
  BoxMin=TPoint3D(ChtMin[16].x-5,ChtMin[16].y-5.,ChtMin[0].z-5.);
  BoxMax=TPoint3D(ChtMax[15].x+5,ChtMax[15].y+5.,ChtMax[16].z+5.);
  //������� ���������� ������� � ��������������� ����� � �����

//����� ������������� ���������� NameDirSrc
//strcpy(NameDirSrc,"D:\\Programming\\C++\\ScMH_SM\\ObrabBMURG_N\\DATA_test\\");
 strcpy(file_read,NameDirSrc);
 strcat(file_read,"\\*.*");
// ������ ���� ����� �������������
  donedir=findfirst(file_read,&fdir,FA_DIREC);
// ���� �� ��������������
  NEventURG=0;
  while (!donedir)
  {
	 if( (fdir.ff_attrib&FA_DIREC)==0)goto m1;
	 if(strcmp(fdir.ff_name,".")==0 ||strcmp(fdir.ff_name,"..")==0)
	 {
		goto m1;
	 }
	 //���������� ��� ���������� ������ ����� ������ ������
	 strcpy(Obrabdir,NameDirSrc);
	 strcat(Obrabdir,"\\");
	 strcat(Obrabdir,fdir.ff_name);
         strcpy(NameRun,fdir.ff_name);
	 // ������� ����� ��������
	 strupr(NameRun);
	 strncpy(buf,NameRun,4);
         buf[4]=0;
	 currentrun=atol(buf);
         if(currentrun==0)goto m1; //����� �����
	 if(Obrab_data(Obrabdir,NameDirWrite))goto m1; //������
         NEventURG +=TD.dtsd.Nevent;
	 printf("\t%s\n",fdir.ff_name);
		 Nrun++;
m1: donedir = findnext(&fdir); //���� ��������� �������������
  }
  findclose(&fdir);
  printf("End of runs\n");
//  getch();
  //����� ����������� ������ �������������� ���������
      for(int ipl=0; ipl < NplMatrix; ipl++)
      {
       int handle=-1;
       int handleN=-1;
       int handleE=-1;
       int handleEL=-1;
       int handleEnoBM=-1;
       sprintf(buf,"XYA_%02ld.dat",int(hMatrix[ipl]));
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"XYN_%02ld.dat",int(hMatrix[ipl]));
       if((handleN = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"XYE_%02ld.dat",int(hMatrix[ipl]));
       if((handleE = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"XYEL_%02ld.dat",int(hMatrix[ipl]));
       if((handleEL = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"XYEnoBM_%02ld.dat",int(hMatrix[ipl]));
       if((handleEnoBM = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }

       if(lseek(handle,0,SEEK_END)==0)
       {
         sprintf(buf,"!Xmin=%lf,Xmax=%lf,Ymin=%lf,Ymax=%lf\n",
                   -NM_X/2*D_X,NM_X/2*D_X,-NM_Y/2*D_Y,NM_Y/2*D_Y);
		 write(handle,buf,strlen(buf));
		 write(handleN,buf,strlen(buf));
         write(handleE,buf,strlen(buf));
         write(handleEL,buf,strlen(buf));
         write(handleEnoBM,buf,strlen(buf));
       }
	   for(int iy=0; iy < NM_Y; iy++)
       {
		 for(int ix=0; ix < NM_X; ix++)
         {
          if(nMatrixXY[ipl][iy][ix])
          {
           MatrixXY[ipl][iy][ix] /=nMatrixXY[ipl][iy][ix];
          }
          if(ix !=NM_X-1)
            sprintf(buf,"%.3lf\t",MatrixXY[ipl][iy][ix]);
          else
			sprintf(buf,"%.3lf\n",MatrixXY[ipl][iy][ix]);
          write(handle,buf,strlen(buf));
          if(ix !=NM_X-1)
			sprintf(buf,"%ld\t",nMatrixXY[ipl][iy][ix]);
		  else
            sprintf(buf,"%ld\n",nMatrixXY[ipl][iy][ix]);
          write(handleN,buf,strlen(buf));

          double eff=0.;
          if(tMatrixXY[ipl][iy][ix])
            eff=double(nMatrixXY[ipl][iy][ix])/tMatrixXY[ipl][iy][ix];
          if(ix !=NM_X-1)
            sprintf(buf,"%lf\t",eff);
          else
            sprintf(buf,"%lf\n",eff);
          write(handleE,buf,strlen(buf));

          eff=0.;
          if(tMatrixXY[ipl][iy][ix])
            eff=double(LnMatrixXY[ipl][iy][ix])/tMatrixXY[ipl][iy][ix];
          if(ix !=NM_X-1)
            sprintf(buf,"%lf\t",eff);
          else
            sprintf(buf,"%lf\n",eff);
          write(handleEL,buf,strlen(buf));
          eff=0.;
          if(tMatrixXY[ipl][iy][ix])
            eff=double(noBMnMatrixXY[ipl][iy][ix])/tMatrixXY[ipl][iy][ix];
          if(ix !=NM_X-1)
            sprintf(buf,"%lf\t",eff);
          else
            sprintf(buf,"%lf\n",eff);
          write(handleEnoBM,buf,strlen(buf));
         }
       }
       close(handle);
       close(handleN);
       close(handleE);
       close(handleEL);
       close(handleEnoBM);
      }
   //��� ������ ������������ ������� ���������
  //����� ����������� ������ �������������� ���������
      for(int ipl=0; ipl < NplMatrix; ipl++)
      {
       int handle=-1;
       int handleN=-1;
       int handleE=-1;
       int handleEL=-1;
       int handleEnoBM=-1;
       sprintf(buf,"RXYA_%02ld.dat",int(hMatrix[ipl]));
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"RXYN_%02ld.dat",int(hMatrix[ipl]));
       if((handleN = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"RXYE_%02ld.dat",int(hMatrix[ipl]));
       if((handleE = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"RXYEL_%02ld.dat",int(hMatrix[ipl]));
       if((handleEL = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"RXYEnoBM_%02ld.dat",int(hMatrix[ipl]));
       if((handleEnoBM = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }

       if(lseek(handle,0,SEEK_END)==0)
       {
         sprintf(buf,"!Xmin=%lf,Xmax=%lf,Ymin=%lf,Ymax=%lf\n",
                   -NM_X/2*D_X,NM_X/2*D_X,-NM_Y/2*D_Y,NM_Y/2*D_Y);
		 write(handle,buf,strlen(buf));
		 write(handleN,buf,strlen(buf));
         write(handleE,buf,strlen(buf));
         write(handleEL,buf,strlen(buf));
         write(handleEnoBM,buf,strlen(buf));
       }
	   for(int iy=0; iy < NM_Y; iy++)
       {
		 for(int ix=0; ix < NM_X; ix++)
         {
          if(RnMatrixXY[ipl][iy][ix])
          {
           RMatrixXY[ipl][iy][ix] /=RnMatrixXY[ipl][iy][ix];
          }
          if(ix !=NM_X-1)
            sprintf(buf,"%.3lf\t",RMatrixXY[ipl][iy][ix]);
          else
			sprintf(buf,"%.3lf\n",RMatrixXY[ipl][iy][ix]);
          write(handle,buf,strlen(buf));
          if(ix !=NM_X-1)
			sprintf(buf,"%ld\t",RnMatrixXY[ipl][iy][ix]);
		  else
            sprintf(buf,"%ld\n",RnMatrixXY[ipl][iy][ix]);
          write(handleN,buf,strlen(buf));

          double eff=0.;
          if(RtMatrixXY[ipl][iy][ix])
            eff=double(RnMatrixXY[ipl][iy][ix])/RtMatrixXY[ipl][iy][ix];
          if(ix !=NM_X-1)
            sprintf(buf,"%lf\t",eff);
          else
            sprintf(buf,"%lf\n",eff);
          write(handleE,buf,strlen(buf));
          eff=0.;
          if(RtMatrixXY[ipl][iy][ix])
            eff=double(RLnMatrixXY[ipl][iy][ix])/RtMatrixXY[ipl][iy][ix];
          if(ix !=NM_X-1)
            sprintf(buf,"%lf\t",eff);
          else
            sprintf(buf,"%lf\n",eff);
          write(handleEL,buf,strlen(buf));
          eff=0.;
          if(RtMatrixXY[ipl][iy][ix])
            eff=double(RnoBMnMatrixXY[ipl][iy][ix])/RtMatrixXY[ipl][iy][ix];
          if(ix !=NM_X-1)
            sprintf(buf,"%lf\t",eff);
          else
            sprintf(buf,"%lf\n",eff);
          write(handleEnoBM,buf,strlen(buf));
         }
       }
       close(handle);
       close(handleN);
       close(handleE);
       close(handleEL);
       close(handleEnoBM);
      }

      //����������� � ����� ���� �������
	  double DD[NplMatrix];
      for(int ipl=0; ipl < NplMatrix; ipl++)
      {
        DD[ipl]=0.;
        //����� ������� �� �������
        double Xaver=0, Y2aver=0, Yaver=0, XYaver=0;
        int a=NM_Y, b=NM_X;
        long YNotEmptyLines=0;
        for(int iy=0; iy < NM_Y; iy++)
		{
			double CoodPoint=0;
            double Noise=0;
            //������� ������� �����
            for(int ix=0; ix < NM_X; ix++)
                if(nMatrixXY[ipl][iy][ix]>0)
				{
                    CoodPoint++;
                    Noise+=nMatrixXY[ipl][iy][ix];
                }
            if(CoodPoint)
                Noise/=CoodPoint;
            else
                Noise=0;
            //��������� ��������� ��������� �� ��������� ����������� ���
			double XPosAver=0;
			long XCount=0;
            for(int ix=0; ix < NM_X; ix++)
				if(nMatrixXY[ipl][iy][ix]>Noise)
                {
                    XPosAver+=nMatrixXY[ipl][iy][ix]*ix;
                    XCount+=nMatrixXY[ipl][iy][ix];
                }
            if(XCount)
            {
                YNotEmptyLines++;
                XPosAver/=XCount;
                Xaver+=XPosAver;
                XYaver+=XPosAver*iy;
                Yaver+=iy;
				Y2aver+=pow(iy,2);
            }
        }
        if(YNotEmptyLines)
        {
            double KK=YNotEmptyLines*Y2aver-pow(Yaver,2);
            double K=1.;
            if(KK)
              K=(YNotEmptyLines*XYaver-Yaver*Xaver)/KK;  //� ����������� �����
            K*=double(D_X)/D_Y;  //� ����������� ��
            DD[ipl]=LenStripY*sin(atan(K));   //�������� � ��
        }
      }
      int handleZim=-1;
      sprintf(buf,"Shift_Angle.dat");
       if((handleZim = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
	  handleZim = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE);
	  sprintf(buf,"Npl, mm\tShift,mm\n");
      write(handleZim,buf,strlen(buf));
      double DDaver=0;
      for(int ipl=0; ipl < NplMatrix; ipl++)
      {
        sprintf(buf,"%02ld\t%lf\n",int(hMatrix[ipl]),DD[ipl]);
        write(handleZim,buf,strlen(buf));
		DDaver+=DD[ipl];
      }
      DDaver/=NplMatrix;
      sprintf(buf,"Average\t%f\n",DDaver);
      write(handleZim,buf,strlen(buf));
      close(handleZim);
  //����������� �������������
  {
	   int handle=-1;
	   sprintf(buf,"AmpLen.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
		 {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"L\terrL\tLall\terrLall\tEff\tA\terrA\tT\terrT\n");
       write(handle,buf,strlen(buf));
       for(int i=0; i < 9; i++)
       {
		 if(nDet[i]==0)continue;
         double nall=nDetAll[i];
         double n=nDet[i];
         double sred=SredLen[i]/n;
         double dis=SigLen[i]/n;
         double err=sqrt(fabs(dis-sred*sred)/n);
         sprintf(buf,"%.1lf\t%.2lf",sred,err);
         write(handle,buf,strlen(buf));
         sred=SredLenAll[i]/nall;
         dis=SigLenAll[i]/nall;
         err=sqrt(fabs(dis-sred*sred)/nall);
         sprintf(buf,"\t%.1lf\t%.2lf",sred,err);
         write(handle,buf,strlen(buf));
         sprintf(buf,"\t%.4lf",n/nall);
         write(handle,buf,strlen(buf));
         sred=AmpDet[i]/n;
         dis=SigDet[i]/n;
         err=sqrt(fabs(dis-sred*sred)/n);
         sprintf(buf,"\t%.1lf\t%.2lf",sred,err);
		 write(handle,buf,strlen(buf));
         sred=TDet[i]/n;
         dis=SigTDet[i]/n;
		 err=sqrt(fabs(dis-sred*sred)/n);
         sprintf(buf,"\t%.1lf\t%.2lf\n",sred,err);
         write(handle,buf,strlen(buf));
       }
       close(handle);
  }
  //������� �������������� � �������
  {
//double APMatrStrip[NCHANADC][NCHANADC]; //������� ���� �������� ���� ������
//double sAPMatrStrip[NCHANADC][NCHANADC]; //������� ���� ��������� �������� ���� ������
//unsigned long NPMatrStrip[NCHANADC][NCHANADC]; //������� ����� �������� ���� ������
//double AtotalStrip[NCHANADC]; //������� ���� �������� �� ������ ����� ����� ���� �������� ����
//double sAtotalStrip[NCHANADC]; //������� ���� ��������� �������� �� ������ ����� ����� ���� �������� ����
//double NtotalStrip[NCHANADC]; //������� ����� ������� ����� ����� ���� �������� ����
//#define NLEN 50 //3460/(34.6*2)
//double MatrLenAStrip[NLEN][NCHANADC]; //����� �������� �� ���������� �� ���
//double MatrLenAPStrip[NLEN][NCHANADC]; //����� �������� � ������� �� ���������� �� ���
//double MatrLenSAStrip[NLEN][NCHANADC]; //����� ��������� �������� �� ���������� �� ���
//double MatrLenSAPStrip[NLEN][NCHANADC]; //����� ��������� �������� � ������� �� ���������� �� ���
//unsigned long MatrLenNStrip[NLEN][NCHANADC]; //����������� ����� ������� �� ���������� �� ���
//unsigned long MatrLenNPStrip[NLEN][NCHANADC]; //����������� ����� ������� � ������� �� ���������
                                                // �� ���������� �� ���
    //������� ������� �������� ���� ������
       int handle=-1;
       sprintf(buf,"M_AP.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"Strip");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
         sprintf(buf,"\tS%02d\terrS%02d",is+1,is+1);
         write(handle,buf,strlen(buf));
       }
       write(handle,"\n",1);
       for(int is=0; is < NCHANADC; is++)
       {
         sprintf(buf,"%d",is+1);
         write(handle,buf,strlen(buf));
         for(int istrip=0; istrip < NCHANADC; istrip++)
		 {
           double sred=0;
           double err=0;
           double n=NPMatrStrip[is][istrip];
           if(n)
		   {
              sred=APMatrStrip[is][istrip]/n;
              err=sAPMatrStrip[is][istrip]/n;
              err=sqrt(fabs(err-sred*sred)/n);
           }
           sprintf(buf,"\t%.2lf\t%.3lf",sred,err);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
	   }
       close(handle);
    //������� ������� �������� ���� ������
       sprintf(buf,"M_M_AP.dat");
	   if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"!Xmin=1.0,Xmax=32.0,Ymin=1.0,Ymax=32.0\n");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
	   {
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           double sred=0;
           double n=NPMatrStrip[is][istrip];
           if(n)
           {
              sred=APMatrStrip[is][istrip]/n;
           }
           if(istrip)
            sprintf(buf,"\t%.2lf",sred);
           else
            sprintf(buf,"%.2lf",sred);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //������� ������������ ������������ � �������
       sprintf(buf,"M_EI.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
	   sprintf(buf,"Strip");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
		 sprintf(buf,"\tS%02d\terrS%02d",is+1,is+1);
         write(handle,buf,strlen(buf));
	   }
       write(handle,"\n",1);
       for(int is=0; is < NCHANADC; is++)
       {
         sprintf(buf,"%d",is+1);
         write(handle,buf,strlen(buf));
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
		   double ver=0;
           double err=0;
		   double n=NtotalStrip[istrip];
           if(n)
           {
              ver=NPMatrStrip[is][istrip]/n;
              err=sqrt(ver*(1.-ver)/n);
           }
		   sprintf(buf,"\t%lf\t%lf",ver,err);
		   write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //������� ������������ ������������ LAM
       sprintf(buf,"M_EILam.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
	   sprintf(buf,"Strip");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
		 sprintf(buf,"\tS%02d\terrS%02d",is+1,is+1);
         write(handle,buf,strlen(buf));
	   }
       write(handle,"\n",1);
       for(int is=0; is < NCHANADC; is++)
       {
         sprintf(buf,"%d",is+1);
         write(handle,buf,strlen(buf));
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
		   double ver=0;
           double err=0;
		   double n=NtotalStrip[istrip];
           if(n)
           {
              ver=LamMatrStrip[is][istrip]/n;
              err=sqrt(ver*(1.-ver)/n);
           }
		   sprintf(buf,"\t%lf\t%lf",ver,err);
		   write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //������� ������������ ��������� LAM � ������������ ������� �� ������
       sprintf(buf,"M_EffStripLamPorog.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
	 sprintf(buf,"Strip\tEffLAM\terrEff\tEffP\terrEffP\n");
         write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
         sprintf(buf,"%d",is+1);
         write(handle,buf,strlen(buf));
	   double ver=0;
           double err=0;
	   double verLAM=0;
           double errLAM=0;
		   double n=NtotalStrip[is];
           if(n)
           {
              ver=NPMatrStrip[is][is]/n;
              err=sqrt(ver*(1.-ver)/n);
              verLAM=MNLamStrip[is]/n;
              errLAM=sqrt(verLAM*(1.-verLAM)/n);
           }
		   sprintf(buf,"\t%lf\t%lf\t%lf\t%lf",verLAM,errLAM,ver,err);
		   write(handle,buf,strlen(buf));
         write(handle,"\n",1);
       }
       close(handle);

       //������� ������������ ������������ � �������
       sprintf(buf,"M_M_EI.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"!Xmin=1.0,Xmax=32.0,Ymin=1.0,Ymax=32.0\n");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           double ver=0;
           double n=NtotalStrip[istrip];
           if(n)
           {
              ver=NPMatrStrip[is][istrip]/n;
		   }
           if(istrip)
			sprintf(buf,"\t%lf",ver);
           else
            sprintf(buf,"%lf",ver);
		   write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //������� ������������ ������������ LAM
       sprintf(buf,"M_M_EILam.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"!Xmin=1.0,Xmax=32.0,Ymin=1.0,Ymax=32.0\n");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           double ver=0;
           double n=NtotalStrip[istrip];
           if(n)
           {
              ver=LamMatrStrip[is][istrip]/n;
		   }
           if(istrip)
			sprintf(buf,"\t%lf",ver);
           else
            sprintf(buf,"%lf",ver);
		   write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //����������� �������� �� �����
       sprintf(buf,"M_ALEN.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"Len");
	   write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
         sprintf(buf,"\tS%02d\terrS%02d",is+1,is+1);
         write(handle,buf,strlen(buf));
       }
       write(handle,"\n",1);
	   for(int in=0; in < NLEN; in++)
       {
		 sprintf(buf,"%lf",(in+0.5)*ChtMax[0].y/NLEN);
         write(handle,buf,strlen(buf));
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           double sred=0;
           double err=0;
           double n=MatrLenNStrip[in][istrip];
           if(n)
           {
              sred=MatrLenAStrip[in][istrip]/n;
              err=MatrLenSAStrip[in][istrip]/n;
              err=sqrt(fabs(err-sred*sred)/n);
           }
		   sprintf(buf,"\t%.2lf\t%.3lf",sred,err);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //������� ����������� �������� �� �����
	   sprintf(buf,"M_M_ALEN.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }

	   sprintf(buf,"!Xmin=1,0,Xmax=32.0,Ymin=0.,Ymax=%lf\n",ChtMax[0].y);
       write(handle,buf,strlen(buf));

       for(int in=0; in < NLEN; in++)
       {
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           double sred=0;
           double n=MatrLenNStrip[in][istrip];
           if(n)
		   {
              sred=MatrLenAStrip[in][istrip]/n;
		   }
           if(istrip)
            sprintf(buf,"\t%.2lf",sred);
           else
            sprintf(buf,"%.2lf",sred);
           write(handle,buf,strlen(buf));
		 }
         write(handle,"\n",1);
       }
	   close(handle);

       //����������� �������� ���� ������ �� �����
       sprintf(buf,"M_ALENP.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
		 {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"Len");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
         sprintf(buf,"\tS%02d\terrS%02d",is+1,is+1);
         write(handle,buf,strlen(buf));
       }
       write(handle,"\n",1);
       for(int in=0; in < NLEN; in++)
       {
         sprintf(buf,"%lf",(in+0.5)*ChtMax[0].y/NLEN);
         write(handle,buf,strlen(buf));
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           double sred=0;
           double err=0;
           double n=MatrLenNPStrip[in][istrip];
           if(n)
		   {
              sred=MatrLenAPStrip[in][istrip]/n;
              err=MatrLenSAPStrip[in][istrip]/n;
              err=sqrt(fabs(err-sred*sred)/n);
           }
		   sprintf(buf,"\t%.2lf\t%.2lf",sred,err);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
	   }
       close(handle);

       //������� ����������� �������� ���� ������ �� �����
	   sprintf(buf,"M_M_ALENP.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"!Xmin=1,0,Xmax=32.0,Ymin=0.,Ymax=%lf\n",ChtMax[0].y);
       write(handle,buf,strlen(buf));
	   for(int in=0; in < NLEN; in++)
	   {
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
		   double sred=0;
           double n=MatrLenNPStrip[in][istrip];
           if(n)
           {
              sred=MatrLenAPStrip[in][istrip]/n;
           }
           if(istrip)
            sprintf(buf,"\t%.2lf",sred);
           else
            sprintf(buf,"%.2lf",sred);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //����������� ����������� ������������ ���� ������ �� �����
       sprintf(buf,"M_ELENP.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"Len");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
	   {
         sprintf(buf,"\tS%02d\terrS%02d",is+1,is+1);
         write(handle,buf,strlen(buf));
       }
       write(handle,"\n",1);
       for(int in=0; in < NLEN; in++)
       {
		 sprintf(buf,"%lf",(in+0.5)*ChtMax[0].y/NLEN);
		 write(handle,buf,strlen(buf));
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           double ver=0;
		   double err=0;
           double n=MatrLenNStrip[in][istrip];
		   if(n)
           {
              ver=MatrLenNPStrip[in][istrip]/n;
              err=sqrt(fabs(ver*(1.-ver)/n));
		   }
           sprintf(buf,"\t%lf\t%lf",ver,err);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
	   }
       close(handle);

       //������� ����������� ����������� ������������ ���� ������ �� �����
       sprintf(buf,"M_M_ELENP.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"!Xmin=1,0,Xmax=32.0,Ymin=0.,Ymax=%lf\n",ChtMax[0].y);
       write(handle,buf,strlen(buf));
       for(int in=0; in < NLEN; in++)
       {
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           double ver=0;
           double n=MatrLenNStrip[in][istrip];
           if(n)
           {
              ver=MatrLenNPStrip[in][istrip]/n;
           }
           if(istrip)
             sprintf(buf,"\t%lf",ver);
           else
			sprintf(buf,"%lf",ver);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
	   close(handle);


	   //������� ����������� ����������� ������������ LAM �� �����
	   sprintf(buf,"M_M_ELENLam.dat");
	   if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
		 {printf("\nError Opening File %s\n",buf); return 0; }
	   sprintf(buf,"!Xmin=1,0,Xmax=32.0,Ymin=0.,Ymax=%lf\n",ChtMax[0].y);
	   write(handle,buf,strlen(buf));
	   for(int in=0; in < NLEN; in++)
	   {
		 for(int istrip=0; istrip < NCHANADC; istrip++)
		 {
		   double ver=0;
		   double n=MatrLenNStrip[in][istrip];
		   if(n)
		   {
			ver=MatrLenLamStrip[in][istrip]/n;
           }
		   if(istrip)
			sprintf(buf,"\t%lf",ver);
		   else
			sprintf(buf,"%lf",ver);
		   write(handle,buf,strlen(buf));
		 }
		 write(handle,"\n",1);
	   }
	   close(handle);


	   //����������� ����� ������� ���� ������ �� �����
	   sprintf(buf,"M_NLENP.dat");
	   if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
		 {printf("\nError Opening File %s\n",buf); return 0; }
	   sprintf(buf,"Len");
	   write(handle,buf,strlen(buf));
	   for(int is=0; is < NCHANADC; is++)
	   {
         sprintf(buf,"\tS%02d",is+1);
         write(handle,buf,strlen(buf));
       }
	   write(handle,"\n",1);
       for(int in=0; in < NLEN; in++)
       {
         sprintf(buf,"%lf",(in+0.5)*ChtMax[0].y/NLEN);
         write(handle,buf,strlen(buf));
		 for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           sprintf(buf,"\t%ld",MatrLenNPStrip[in][istrip]);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //������� ����� ������� ���� ������ �� �����
       sprintf(buf,"M_M_NLENP.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
	   sprintf(buf,"!Xmin=1,0,Xmax=32.0,Ymin=0.,Ymax=%lf\n",ChtMax[0].y);
       write(handle,buf,strlen(buf));
       for(int in=0; in < NLEN; in++)
       {
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           if(istrip)
			 sprintf(buf,"\t%ld",MatrLenNPStrip[in][istrip]);
           else
			sprintf(buf,"%ld",MatrLenNPStrip[in][istrip]);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //����������� ����� ������� �� �����
	   sprintf(buf,"M_NLEN.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"Len");
       write(handle,buf,strlen(buf));
       for(int is=0; is < NCHANADC; is++)
       {
         sprintf(buf,"\tS%02d",is+1);
         write(handle,buf,strlen(buf));
       }
	   write(handle,"\n",1);
       for(int in=0; in < NLEN; in++)
       {
         sprintf(buf,"%lf",(in+0.5)*ChtMax[0].y/NLEN);
		 write(handle,buf,strlen(buf));
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           sprintf(buf,"\t%ld",MatrLenNStrip[in][istrip]);
           write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);

       //������� ����� ������� �� �����
       sprintf(buf,"M_M_NLEN.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
         {printf("\nError Opening File %s\n",buf); return 0; }
	   sprintf(buf,"!Xmin=1,0,Xmax=32.0,Ymin=0.,Ymax=%lf\n",ChtMax[0].y);
	   write(handle,buf,strlen(buf));
       for(int in=0; in < NLEN; in++)
       {
         for(int istrip=0; istrip < NCHANADC; istrip++)
         {
           if(istrip)
             sprintf(buf,"\t%ld",MatrLenNStrip[in][istrip]);
           else
			sprintf(buf,"%ld",MatrLenNStrip[in][istrip]);
		   write(handle,buf,strlen(buf));
         }
         write(handle,"\n",1);
       }
       close(handle);
/*
double NTrackAll=0.; //����� ����� ������ ��������� ����� ������
double NAmpPorogAll=0.; //����� ����� ������� � ������ �������, ����� �� ����� ������ ���� ���������� ������/
//��� Theta 0-15
double NTrackStrip=0.; //����� ����� ������ ��������� ����� ����� ������ ������ �������� ������
double NAmpPorogStrip=0.; //����� ����� ������� � ������� ����� ����� ������ ������, ����� �� ����� ������ ���� ���������� ������/
double NTrackBox=0.; //����� ����� ������ ��������� ����� ������� �������� ������
double NAmpPorogBox=0.; //����� ����� ������� � �������, ����� �� ����� ������ ���� ���������� ������/
*/

//���������� ������������� ������������ � BMeff.dat
	   sprintf(buf,"BMeff.dat");
       if((handle = open(buf,O_RDWR | O_TEXT | O_CREAT|O_TRUNC, S_IWRITE)) == -1)
		 {printf("\nError Opening File %s\n",buf); return 0; }
       sprintf(buf,"������������ �������� ��������� ������ ������\n");
	   write(handle,buf,strlen(buf));
       sprintf(buf,"DX= %lg , DY= %lg , DZ= %lg , DD= %lg\n",BMSmesX,BMSmesY,BMSmesZ,BMSmesD);
       write(handle,buf,strlen(buf));
       sprintf(buf,"����� ����� ��������� ������ �� ��� �������� �������\n");
       write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NEventURG);
       write(handle,buf,strlen(buf));
//       sprintf(buf,"����� ����� ��������� ������ � ����������� ������ ������� ��\n");
//       write(handle,buf,strlen(buf));
//	   sprintf(buf,"%lu\n",(unsigned long)NTrigURG);
//       write(handle,buf,strlen(buf));
       sprintf(buf,"���������� ������� � ������� �� ��\n");
       write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NDataBM);
       write(handle,buf,strlen(buf));
       sprintf(buf,"���������� ������� � ���������� ������� �� ��\n");
       write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NDataBMBad);
       write(handle,buf,strlen(buf));
       sprintf(buf,"������������� ������� ������ �� �� ��� �������� ������\n");
       write(handle,buf,strlen(buf));
       double p=0,errp=0;
       if(NEventURG)
       {
//        p=NDataBM/NTrigURG;
//		errp=sqrt(p*(1-p)/NTrigURG);
        p=NDataBM/NEventURG;
		errp=sqrt(p*(1-p)/NEventURG);
       }
       sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
       write(handle,buf,strlen(buf));

	   sprintf(buf,"\n���������� ��� ������������ ������ ������ � ������� %lf ������� ���\n", Par.porogImpulse[1]);
       write(handle,buf,strlen(buf));
       sprintf(buf,"\n����� ����� ������ ��������� ����� ������\n");
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NTrackAll);
       write(handle,buf,strlen(buf));
       sprintf(buf,"���������� ������� �� �������������� ��\n");
	   write(handle,buf,strlen(buf));
       sprintf(buf,"%lu\n",(unsigned long)NAmpPorogAll);
       write(handle,buf,strlen(buf));
       sprintf(buf,"������������� ������������ �� ���� ������\n");
       write(handle,buf,strlen(buf));
       p=0;errp=0;
       if(NTrackAll)
       {
		p=NAmpPorogAll/NTrackAll;
        errp=sqrt(p*(1-p)/NTrackAll);
       }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));


	   sprintf(buf,"\n\n    ��� ������ � �������� ����� 0-15 ��������\n");
       write(handle,buf,strlen(buf));
	   sprintf(buf,"\n����� ������ ��������� ����� ������� ��\n");
       write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NTrackBox);
       write(handle,buf,strlen(buf));
//������� ������ �� �� ������ � �������� ����� 0-15 �������� ��������� ����� ������� ��
	   sprintf(buf,"���������� ������� � ������� �� �� ������, ��������� ����� ������� ��\n");
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NDataBox);
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"������������� ������� ������ �� �� ������, ��������� ����� ������� ��\n");
	   write(handle,buf,strlen(buf));
	   p=0;errp=0;
	   if(NTrackBox)
	   {
			p = NDataBox/NTrackBox;
			errp=sqrt(p*(1-p)/NTrackBox);
	   }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));

       sprintf(buf,"���������� ������� �� �������������� ��\n");
       write(handle,buf,strlen(buf));
       sprintf(buf,"%lu\n",(unsigned long)NAmpPorogBox);
       write(handle,buf,strlen(buf));
       sprintf(buf,"������������� ������������ �� �� ������ ����� ������� ��\n");
       write(handle,buf,strlen(buf));
       p=0;errp=0;
	   if(NTrackBox)
       {
		p=NAmpPorogBox/NTrackBox;
        errp=sqrt(p*(1-p)/NTrackBox);
       }
       sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
       write(handle,buf,strlen(buf));

//����� LAM �� ������ � �������� ����� 0-15 �������� ��������� ����� ������� ��
	   sprintf(buf,"���������� ����� LAM �� ������ ��������� ����� ������� ��\n");
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NLamBox);
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"������������� ����� LAM �� ������ ��������� ����� ������� ��\n");
	   write(handle,buf,strlen(buf));
	   p=0;errp=0;
	   if(NTrackBox)
	   {
			p = NLamBox/NTrackBox;
			errp=sqrt(p*(1-p)/NTrackBox);
	   }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));

//����� LAM ��� ������������� �� �� ������ � �������� ����� 0-15 �������� ��������� ����� ������� ��
	   sprintf(buf,"���������� LAM �� ������������� �� �� ������ ��������� ����� ������� ��\n");
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NLamAmpPorogBox);
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"������������� LAM ��� ������������� �� �� ������ ��������� ����� ������� ��\n");
	   write(handle,buf,strlen(buf));
	   p=0;errp=0;
	   if(NAmpPorogBox)
	   {
			p = NLamAmpPorogBox/NAmpPorogBox;
			errp=sqrt(p*(1-p)/NAmpPorogBox);
	   }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));


       sprintf(buf,"\n����� ������ ��������� ����� ����� ������ ������\n");
       write(handle,buf,strlen(buf));
       sprintf(buf,"%lu\n",(unsigned long)NTrackStrip);
       write(handle,buf,strlen(buf));
	   sprintf(buf,"���������� ������� �� �������������� ��\n");
       write(handle,buf,strlen(buf));
       sprintf(buf,"%lu\n",(unsigned long)NAmpPorogStrip);
       write(handle,buf,strlen(buf));
	   sprintf(buf,"������������� ������������ �� �� ������ ����� ����� ������ ������\n");
       write(handle,buf,strlen(buf));
       p=0;errp=0;
	   if(NTrackStrip)
	   {
		p=NAmpPorogStrip/NTrackStrip;
		errp=sqrt(p*(1-p)/NTrackStrip);
	   }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));

//LAM �� ������ � �������� ����� 0-15 �������� ��������� ����� ����� (+-4 ��) ������ ������
	   sprintf(buf,"���������� ����� LAM �� ������ ��������� ����� ����� (+-4 ��) ������ ������\n");
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NLamStrip);
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"������������� ��������� LAM �� ������ ��������� ����� ����� (+-4 ��) ������ ������\n");
	   write(handle,buf,strlen(buf));
	   p=0;errp=0;
	   if(NTrackStrip)
	   {
			p = NLamStrip/NTrackStrip;
			errp=sqrt(p*(1-p)/NTrackStrip);
	   }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));

//LAM �� ������� ����� ����� (+-4 ��) ������� ������ ���� (0-15 �������� �� ������)
       sprintf(buf,"\n����� ����������� ������� (+-4 ��) ������� (0-15 �������� �� ������)\n");
       write(handle,buf,strlen(buf));
       sprintf(buf,"%lu\n",(unsigned long)NTrackStripFull);
       write(handle,buf,strlen(buf));
	   sprintf(buf,"���������� LAM �� ������� ����� ����� (+-4 ��) ������� ������ ����\n");
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NLamCurrStrip);
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"������������� ��������� LAM �� ������� ����� ����� (+-4 ��) ������� ������ ����\n");
	   write(handle,buf,strlen(buf));
	   p=0;errp=0;
	   if(NTrackStripFull)
	   {
			p = NLamCurrStrip/NTrackStripFull;
			errp=sqrt(p*(1-p)/NTrackStripFull);
	   }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));

//������������ ������� ����� ����� (+-4 ��) ������� ������ ���� (0-15 �������� �� ������)
	   sprintf(buf,"���������� ����������� ������� ����� ����� (+-4 ��) ������� ������ ����\n");
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NAmpCurrStrip);
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"������������� ������������ ������� ����� ����� (+-4 ��) ������� ������ ����\n");
	   write(handle,buf,strlen(buf));
	   p=0;errp=0;
	   if(NTrackStrip)
	   {
			p = NAmpCurrStrip/NTrackStripFull;
			errp=sqrt(p*(1-p)/NTrackStripFull);
	   }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));

//LAM �� ����������� ������� ����� ����� (+-4 ��) ������� ������ ���� (0-15 �������� �� ������)
	   sprintf(buf,"���������� LAM �� ����������� ������� ����� ����� (+-4 ��) ������� ������ ����\n");
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"%lu\n",(unsigned long)NLamAmpPorogCurrStrip);
	   write(handle,buf,strlen(buf));
	   sprintf(buf,"������������� ��������� LAM �� ����������� ������� ����� ����� (+-4 ��) ������� ������ ����\n");
	   write(handle,buf,strlen(buf));
	   p=0;errp=0;
	   if(NAmpCurrStrip)
	   {
			p = NLamAmpPorogCurrStrip/NAmpCurrStrip;
			errp=sqrt(p*(1-p)/NAmpCurrStrip);
	   }
	   sprintf(buf,"%.2lf +- %.3lf %%\n",p*100.,errp*100.);
	   write(handle,buf,strlen(buf));


	   close(handle);

  }
 if(nMatrCrossHitStripPlane[0])
 {
  std::ofstream out;
    for(size_t iBM=0; iBM < MatrCrossHitStripPlane.size(); iBM++)
    {
      sprintf(buf,"MatrCrossStrip_%02ld.dat",iBM);
      out.open(buf);
      sprintf(buf,"!Xmin=1,Xmax=%lg,Ymin=1,Ymax=%lg Nrec=",double(Nstrip),
                            double(Nstrip));
      out << buf << nMatrCrossHitStripPlane[iBM] <<"\n";

      for(int i=0; i < Nstrip; i++)
      {
       for(int j=0; j < Nstrip; j++)
       {
        unsigned long n=MatrCrossHitStripPlane[iBM][i][j];
        if(j)
         sprintf(buf,"\t%ld",n);
        else
         sprintf(buf,"%ld",n);
        out << buf;
       }
       out << "\n";
      }
      out.close();
    }
 }
 if(nMatrCrossLamStripPlane[0])
 {
  std::ofstream out;
    for(size_t iBM=0; iBM < MatrCrossLamStripPlane.size(); iBM++)
    {
      sprintf(buf,"MatrCrossStripLam_%02ld.dat",iBM);
      out.open(buf);
      sprintf(buf,"!Xmin=1,Xmax=%lg,Ymin=1,Ymax=%lg Nrec=",double(Nstrip),
                            double(Nstrip));
      out << buf << nMatrCrossLamStripPlane[iBM] <<"\n";

      for(int i=0; i < Nstrip; i++)
      {
       for(int j=0; j < Nstrip; j++)
       {
        unsigned long n=MatrCrossLamStripPlane[iBM][i][j];
        if(j)
         sprintf(buf,"\t%ld",n);
        else
         sprintf(buf,"%ld",n);
        out << buf;
       }
       out << "\n";
      }
      out.close();
    }
 }
 if(nMatrCrossHitLamStripPlane[0])
 {
  std::ofstream out;
    for(size_t iBM=0; iBM < MatrCrossHitLamStripPlane.size(); iBM++)
    {
      sprintf(buf,"MatrCrossStripHitLam_%02ld.dat",iBM);
      out.open(buf);
      sprintf(buf,"!Xmin=1,Xmax=%lg,Ymin=1,Ymax=%lg Nrec=",double(Nstrip),
                            double(Nstrip));
      out << buf << nMatrCrossHitLamStripPlane[iBM] <<"\n";

      for(int i=0; i < Nstrip; i++)
      {
       for(int j=0; j < Nstrip; j++)
       {
        unsigned long n=MatrCrossHitLamStripPlane[iBM][i][j];
        if(j)
         sprintf(buf,"\t%ld",n);
        else
         sprintf(buf,"%ld",n);
        out << buf;
       }
       out << "\n";
      }
      out.close();
    }
 }

  chdir(NameDirStart);
   return 0;
}
//---------------------------------------------------------------------------
// ��������� ������ �� ���������� DirRead -> ��������� � ���������� DirWrite
short Obrab_data(char *DirRead,char *DirWrite)
{
//��������� �� ���� ����� ����� �������� �
// �������� ������ �� ������������� DirRead, ������ � ��� ������
// � ������ � ��������� ���� ������������� �� Phi ���������� �������
// ��������������� �� ����� ����� �����, ��� ��������� ��������� Theta
 char NameFileSearch[256];
 char NameFileRead[256];
 char NameFileWrite[256];
 int hand_out[4],done;
 struct ffblk ffile; // ��� ������ ������
 struct ffblk fdir; // ��� ������ ��� ����������
 double torad=atan(1.)/45.; //������� � �������
    done=findfirst(DirWrite,&fdir,FA_DIREC); //���������
    findclose(&fdir);
    if(done) // ���������� �� �������
    {
	  printf("No make subdir %s",DirWrite); return 1;
    }

    done=findfirst(DirRead,&fdir,FA_DIREC); //���������
    findclose(&fdir);
    if(done) // ���������� �� �������
	{
      printf("No make subdir %s",DirRead); return 1;
	}

 strcpy(NameFileSearch,DirRead);
 strcat(NameFileSearch,"\\*.uad");
 done=findfirst(NameFileSearch,&ffile,FA_ARCH);
 if(done)return 1;
 //��� �����������
   strcpy(NameFileWrite,DirWrite);
   strcat(NameFileWrite,"\\");
   sprintf(NameFileWrite+strlen(NameFileWrite),"All_StatUAD.dat");
   hand_out[0]=open(NameFileWrite,O_RDWR|O_CREAT|O_TEXT,S_IWRITE);
   if(hand_out[0]==-1)
   {
     printf("\r\n   No create %s",NameFileWrite); return 1;
   }
   strcpy(NameFileWrite,DirWrite);
   strcat(NameFileWrite,"\\");
   sprintf(NameFileWrite+strlen(NameFileWrite),"T_StatUAD.dat");
   hand_out[1]=open(NameFileWrite,O_RDWR|O_CREAT|O_TEXT,S_IWRITE);
   if(hand_out[1]==-1)
   {
     close(hand_out[0]);
     printf("\r\n   No create %s",NameFileWrite); return 1;
   }
   strcpy(NameFileWrite,DirWrite);
   strcat(NameFileWrite,"\\");
   sprintf(NameFileWrite+strlen(NameFileWrite),"TR_StatUAD.dat");
   hand_out[2]=open(NameFileWrite,O_RDWR|O_CREAT|O_TEXT,S_IWRITE);
   if(hand_out[2]==-1)
   {
     close(hand_out[0]);
     close(hand_out[1]);
     printf("\r\n   No create %s",NameFileWrite); return 1;
   }
   strcpy(NameFileWrite,DirWrite);
   strcat(NameFileWrite,"\\");
   sprintf(NameFileWrite+strlen(NameFileWrite),"TRD_StatUAD.dat");
   hand_out[3]=open(NameFileWrite,O_RDWR|O_CREAT|O_TEXT,S_IWRITE);
   if(hand_out[3]==-1)
   {
     close(hand_out[0]);
     close(hand_out[1]);
     close(hand_out[2]);
	 printf("\r\n   No create %s",NameFileWrite); return 1;
   }
   lseek(hand_out[1],0,SEEK_END);
   lseek(hand_out[2],0,SEEK_END);
   lseek(hand_out[3],0,SEEK_END);
   if(lseek(hand_out[0],0,SEEK_END)==0)
   {
     sprintf(buf,"Nrun\tNevent\tidTrack\tPhi\tTheta\tLrad\tPds\tsigPds\tAmp\tT\n");
     write(hand_out[0],buf,strlen(buf));
	 write(hand_out[1],buf,strlen(buf));
	 write(hand_out[2],buf,strlen(buf));
     write(hand_out[3],buf,strlen(buf));
     sprintf(buf,"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
     write(hand_out[0],buf,strlen(buf));
     write(hand_out[1],buf,strlen(buf));
     write(hand_out[2],buf,strlen(buf));
     write(hand_out[3],buf,strlen(buf));
     sprintf(buf,"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
     write(hand_out[0],buf,strlen(buf));
     write(hand_out[1],buf,strlen(buf));
     write(hand_out[2],buf,strlen(buf));
	 write(hand_out[3],buf,strlen(buf));
     sprintf(buf,"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
     write(hand_out[0],buf,strlen(buf));
     write(hand_out[1],buf,strlen(buf));
     write(hand_out[2],buf,strlen(buf));
     write(hand_out[3],buf,strlen(buf));
     sprintf(buf,"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
     write(hand_out[0],buf,strlen(buf));
     write(hand_out[1],buf,strlen(buf));
     write(hand_out[2],buf,strlen(buf));
     write(hand_out[3],buf,strlen(buf));
   }
   int hamp[3]={-1,-1,-1};
   char *comamp[3]={"1020","4050","7080"};
   for(int i=0; i < 3; i++)
   {
    strcpy(NameFileWrite,DirWrite);
    strcat(NameFileWrite,"\\");
	sprintf(NameFileWrite+strlen(NameFileWrite),"AmpSpc_%s.dat",comamp[i]);
    hamp[i]=open(NameFileWrite,O_RDWR|O_CREAT|O_TEXT,S_IWRITE);
	if(lseek(hamp[i],0,SEEK_END)==0)
	{
	 sprintf(buf,"Nrun\tNevent\tidTrack\tPhi\tTheta\tLrad\tPds\tsigPds\tAmp\tT\n");
	 write(hamp[i],buf,strlen(buf));
	 sprintf(buf,"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
	 write(hamp[i],buf,strlen(buf));
	 sprintf(buf,"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
	 write(hamp[i],buf,strlen(buf));
	 sprintf(buf,"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
	 write(hamp[i],buf,strlen(buf));
	 sprintf(buf,"0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
	 write(hamp[i],buf,strlen(buf));
	}
   }
 int idBM=0; //���� ������� ������ � ��
 do
 {
   strcpy(NameFileRead,DirRead);
   strcat(NameFileRead,"\\");
   strcat(NameFileRead,ffile.ff_name);
   ReadBufTD(NameFileRead);
//   nBufTD = 1000; //?????
//   printf("\nEvents in buffer: %d\n", nBufTD);
   //������������ �� ��������
   for(int i=0; i < nBufTD; i++)
   {
     memcpy(&TD,&BufTD[i],sizeof(TD));
     if(TD.dtsd.start[0]==0)continue;
     if(idChanStrip)
      memcpy(TD.ChanStrip,ChanStrip,sizeof(ChanStrip));
     NTrigURG++;
     //�������� �� ������ ASIC'�
     idBM=1; //�� ��������� ������� �� � �������
     {
       int MASKAD=1;
       short d=TD.Data[0];
       for(int i=1; i < NCHANADC; i++)
       {
         if(TD.Data[i] == d)MASKAD++;
       }
       if(MASKAD==32)idBM=0; //�� ��������� ������ �� ASIC'�
       if(TD.Data[0] == 15 || TD.Data[0] == 16)NDataBMBad++;
     }

     if(TD.Data[0]==0)idBM=0; //�� ��������� ������ �� ASIC'�
     NDataBM +=idBM;
     //���������, ��� ������ ����

     //��������� ������ ������������
	 // chan - ����� ������
     // TSD - �������� ������ ������������
     // Par - ��������� ���������
     // �� ������:
     // ������������ ��������: 1- ��������� ������ ������� ������� ������, 0 - �������
     // Res - ���������� ���������
     RES_OBRABTSD Res;
     RES_OBRABTSD ResTemp;
	 memset(&Res,0,sizeof(Res));
	 memset(&ResTemp,0,sizeof(ResTemp));
	 int idimp=0;
     int chanmin=0;
     int chanmax=NCHANADC;
     int stripmin=0,stripmax=NCHANADC;
      //�������� ������� ����������� ������������
     if(idBM)
      for(int iBM=0; iBM < UseNBM; iBM++)
      {
       for(int chan=0; chan < Nstrip; chan++)
       {
	   int ipmt=ADCtoPMT[chan]-1; //����� ������ ���
	   int istrip=TD.ChanStrip[ipmt];
	   int ifdimp=0;
           memset(&ResTemp,0,sizeof(ResTemp));
           ifdimp=Obrab_TSD(chan, TD, Par, ResTemp);
         if(ifdimp==0)continue;
         int iX=istrip; //����� ������������ ������
         for(int chanC=0; chanC < Nstrip; chanC++)
         {
	   int ipmtC=ADCtoPMT[chanC]-1; //����� ������ ���
	   int istripC=TD.ChanStrip[ipmtC];
	   int ifdimpC=0;
           memset(&ResTemp,0,sizeof(ResTemp));
           ifdimpC=Obrab_TSD(chanC, TD, Par, ResTemp);
           if(ifdimpC==0)continue;
           int iY=istripC;
           MatrCrossHitStripPlane[iBM][iY][iX]++;
           nMatrCrossHitStripPlane[iBM]++;
         }
       }
       //�� Lam
       for(int chan=0; chan < Nstrip; chan++)
       {
	   int ipmt=ADCtoPMT[chan]-1; //����� ������ ���
	   int istrip=TD.ChanStrip[ipmt];
         if(TD.Lam[chan]==0)continue;
         int iX=istrip; //����� ������������ ������
         for(int chanC=0; chanC < Nstrip; chanC++)
         {
	   int ipmtC=ADCtoPMT[chanC]-1; //����� ������ ���
	   int istripC=TD.ChanStrip[ipmtC];
           if(TD.Lam[chanC])
           {
            int iY=istripC;
            MatrCrossLamStripPlane[iBM][iY][iX]++;
            nMatrCrossLamStripPlane[iBM]++;
           }
         }
       }
       //��� ������ �� Lam
       for(int chan=0; chan < Nstrip; chan++)
       {
	   int ipmt=ADCtoPMT[chan]-1; //����� ������ ���
	   int istrip=TD.ChanStrip[ipmt];
	   int ifdimp=0;
           memset(&ResTemp,0,sizeof(ResTemp));
           ifdimp=Obrab_TSD(chan, TD, Par, ResTemp);
         if(ifdimp==0)continue;
         if(TD.Lam[chan]==0)continue;
         int iX=istrip; //����� ������������ ������
         for(int chanC=0; chanC < Nstrip; chanC++)
         {
	   int ipmtC=ADCtoPMT[chanC]-1; //����� ������ ���
	   int istripC=TD.ChanStrip[ipmtC];
           if(TD.Lam[chanC])
           {
            int iY=istripC;
            MatrCrossHitLamStripPlane[iBM][iY][iX]++;
            nMatrCrossHitLamStripPlane[iBM]++;
           }
         }
       }

      }

     if(ChanProc >=0)
     {
       stripmin=ChanProc; stripmax=ChanProc+1;
	 }
	 double ampstrip[NCHANADC]; //��������� �� �������
	 long lamstrip[NCHANADC]; //���� LAM ��� ������. ��������� 22.06.2012 18:36
	 memset(ampstrip,0,sizeof(ampstrip));
	 memset(lamstrip,0,sizeof(lamstrip));
	 long LamAny=0;
	 for(int chan=chanmin; chan < chanmax; chan++)
	 {
	   int ipmt=ADCtoPMT[chan]-1; //����� ������ ���
	   int istrip=TD.ChanStrip[ipmt];
	   int ifdimp=0;
           memset(&ResTemp,0,sizeof(ResTemp));
           if(idBM)
             ifdimp=Obrab_TSD(chan, TD, Par, ResTemp);
	   ampstrip[istrip]=ResTemp.amppik;
	   lamstrip[istrip] = ResTemp.lam;	//���� LAM ��� ������. ��������� 22.06.2012 18:37
	   LamAny += ResTemp.lam;
//	   LamAny += lamstrip[istrip];
	   if(istrip <stripmin || istrip >= stripmax)continue;
	   if(ifdimp)
	   {
		 Res.pds +=ResTemp.pds;
		 Res.sigmapds +=ResTemp.sigmapds*ResTemp.sigmapds;
		 Res.amppik +=ResTemp.amppik;
		 Res.timp +=ResTemp.timp;
		 Res.T +=ResTemp.T;
		 Res.lam +=	ResTemp.lam;
		 idimp++;
//		 printf("\nidimp = %d\tResTemp.lam = %d\t\tRes.lam = %d\tLamAny = %d\tRes.amppik = %f\n",idimp,ResTemp.lam,Res.lam,LamAny,Res.amppik); //???
	   }
	 }
	 if(idimp)
	 {
		 Res.pds /=idimp;
		 Res.sigmapds=sqrt(Res.sigmapds/idimp);
		 Res.amppik /=idimp;
		 Res.timp /=idimp;
		 Res.Q /=idimp;
		 Res.T /=idimp;
	 }
     //����� ����������� �� ��������������

     TPoint3D P1,P2,cross1,cross2;
     double p0[3],p1[3];
     int kpointC=0;
     int kpointCB=0;
	 double L=0.; //����� ����������� ����� � ���������
	 double ThetaG=-1.,PhiG=-1.;
	 double dx,dy,dz;
	 int iz=1;
         int TrackInBox=0;
	 if( (TD.dtsd.id_trek&0xFF) !=3)
	 {
		NTrackAll++;

		if(idimp)NAmpPorogAll++;

		TPoint3D Ce0;
		int Sm=SMwork[0];
		int stripX0=(TD.X0_BM[0]-1)*16+TD.X0_BM[1]-1;
		int stripY0=(TD.Y0_BM[0]-1)*16+TD.Y0_BM[1]-1;
		Ce0.x=Pdcr7[Sm][0]+VdcrZ[Sm][0]*TablXY[Sm][7*NPOINTMAX+stripX0].z; //plane "X"
		Ce0.y=Pdcr7[Sm][1]+VdcrX[Sm][1]*TablXY[Sm][7*NPOINTMAX+stripY0].x; //plane "Y"
		Ce0.z=Pdcr7[Sm][2]+5.;//+51.2;//+42.;
		P1.x=TD.dtsd.chp0[0]; P1.y=TD.dtsd.chp0[1]; P1.z=TD.dtsd.chp0[2];
		P2.x=TD.dtsd.chp1[0]; P2.y=TD.dtsd.chp1[1]; P2.z=TD.dtsd.chp1[2];
		//������� ���������� ����� ������������ ������� �����
		P1.x-=Ce0.x+BMSmesX; P1.y-=Ce0.y+BMSmesY; P1.z-=Ce0.z+BMSmesZ;
		P2.x-=Ce0.x+BMSmesX; P2.y-=Ce0.y+BMSmesY; P2.z-=Ce0.z+BMSmesZ;
		//�������� ������������ ��������� �� ���� BMSmesD ������������ ����� {0,0,0}
		double dPhi=atan(BMSmesD/LenStripY);
		P1=AxisTurn(P1,dPhi);
		P2=AxisTurn(P2,dPhi);
   //--------------------------------------------------------
	  dx=P2.x-P1.x;
	  dy=P2.y-P1.y;
	  dz=P2.z-P1.z;
	  double r=sqrt(dx*dx+dy*dy+dz*dz);
	  if(dz < 0)
	  {
		iz=-1;
		dx =-dx;
		dy =-dy;
		dz =-dz;
	  }
	  PhiG=0;
	  dx /=r; dy /=r; dz /=r;
	  if(dx||dy)
		 PhiG=atan2(dy,dx)/torad;
	  if(PhiG < 0.)PhiG +=360;
	  ThetaG=acos(dz)/torad;
	  p0[0]=P1.x; p0[1]=P1.y; p0[2]=P1.z;
	  p1[0]=P2.x; p1[1]=P2.y; p1[2]=P2.z;
	  if(ChanProc >=0)
	   kpointC=SetPointBox(p0, p1, ChtMin[ChanProc], ChtMax[ChanProc]);
	  else
	  if(ThetaG >= 0. && ThetaG < 15.)
	  {
	   int instrip=0;
		for(int istrip=0; istrip < NCHANADC; istrip++)
		{
		  TPoint3D cs0,cs1, csc;
		  double ps0[3],ps1[3];
		  memcpy(ps0,p0,sizeof(p0));
		  memcpy(ps1,p1,sizeof(p1));
		  //������� ����������� ����� � ���������� �������
		  cs0=ChtMin[istrip];
		  cs1=ChtMax[istrip];
                  csc.x=(cs0.x+cs1.x)/2.;
		  cs0.x =csc.x-4; cs1.x =csc.x+4;
		  int kp=SetPointBox(ps0, ps1, cs0, cs1);
		  if(kp < 2)continue; //���� �� ������ ���� �����
		  instrip++;
		  double ys=(ps0[1]+ps1[1])/2.;
		  int idzone=(ys/ChtMax[istrip].y*4);
		  if(idzone > 3)idzone =3;
		  double amp=ampstrip[istrip]; //��������� � ����� ������
		  LamCurrStrip=lamstrip[istrip]; //Lam �� ������ ������
		  //int idlen=(ChtMax[istrip].y-ys)/ChtMax[istrip].y*NLEN;
		  int idlen=ys/ChtMax[istrip].y*NLEN;
		  if(idlen >=NLEN)idlen=NLEN-1;
		  if(idlen < 0)idlen=0;
		  AtotalStrip[istrip] +=amp; //��� ���������
		  sAtotalStrip[istrip] +=amp*amp; //��� ���������
		  NtotalStrip[istrip]++;
		  MatrLenAStrip[idlen][istrip]+=amp; //����������� �� ���������� �� ���
		  MatrLenSAStrip[idlen][istrip]+=amp*amp; //����������� �� ���������� �� ���
		  MatrLenNStrip[idlen][istrip]++; //����������� �� ���������� �� ���

		  if(LamCurrStrip)	//���� ���� Lam
		  {
			MatrLenLamStrip[idlen][istrip]++;
			MNLamStrip[istrip]++;
			NLamCurrStrip++;
//			printf("NLamCurrStrip: %2.f\tistrip: %d\n",NLamCurrStrip,istrip+1);
		  }
		  if(amp >=Par.porogImpulse[1])//���� ������
		  {
			MatrLenAPStrip[idlen][istrip]+=amp; //����������� ��� ������ �� ���������� �� ���
			MatrLenSAPStrip[idlen][istrip]+=amp*amp; //����������� ��� ������ �� ���������� �� ���
			MatrLenNPStrip[idlen][istrip]++;
			NAmpCurrStrip++;
			if(LamCurrStrip)NLamAmpPorogCurrStrip++;
//			printf("NAmpCurrStrip: %2.f\tNLamAmpPorogCurrStrip: %2.f\ti: %d\n",NAmpCurrStrip,NLamAmpPorogCurrStrip,i+1);
		  }
		  for(int is=0; is < NCHANADC; is++)
		  {
                   if(lamstrip[is])
                   {
		    LamMatrStrip[is][istrip]++;
                   }
		   if(ampstrip[is] < Par.porogImpulse[1])continue;
		   APMatrStrip[is][istrip]+=ampstrip[is];
		   sAPMatrStrip[is][istrip]+=ampstrip[is]*ampstrip[is];
		   NPMatrStrip[is][istrip]++;
		  }

//	���������� �����	Astrip_
		  strcpy(NameFileWrite,DirWrite);
		  strcat(NameFileWrite,"\\");
		  sprintf(NameFileWrite+strlen(NameFileWrite),"Astrip_%02d_%01d.dat",istrip+1,idzone);
		  int hxyz=open(NameFileWrite,O_RDWR|O_CREAT|O_TEXT,S_IWRITE);
		  if(lseek(hxyz,0,SEEK_END)==0)
		  {
			  for(int is=0; is < NCHANADC; is++)
			  {
			   if(is==0)
				sprintf(buf,"A%02d\tL%02d",is+1,is+1);
			   else
				sprintf(buf,"\tA%02d\tL%02d",is+1,is+1);
			   write(hxyz,buf,strlen(buf));
			  }
			  write(hxyz,"\n",1);
		  }
		  for(int is=0; is < NCHANADC; is++)
		  {
			if(is==0)
			 sprintf(buf,"%.2lf\t%1d",ampstrip[is],lamstrip[is]);
			else
			 sprintf(buf,"\t%.2lf\t%1d",ampstrip[is],lamstrip[is]);
			write(hxyz,buf,strlen(buf));
		  }
		  write(hxyz,"\n",1);
		  close(hxyz);
		}
		if(instrip)
		{
			NTrackStrip++;
			NTrackStripFull+=instrip;
			if(LamAny)NLamStrip++;	//���� ���� LAM ��� ������� ����� ����� ����������� ����� ������ �� �������������� ����� ���� LAM. ��������� 25.06.12 17:22
//			if(lamstrip[istrip])NLamCurrStrip++;	//���� ���� LAM ��� ������� ����� ����� ����������� ����� ������ � ������������ �� �� �������������� ����� ���� LAM. ��������� 26.06.12 17:44
			if(idimp)NAmpPorogStrip++;
		}
		int kp=SetPointBox(p0, p1, BoxMin, BoxMax);
                if(ChanProc < 0)
                  kpointCB=kp;
	   if(kp >=2)
	   {
               TrackInBox=1;
		NTrackBox++;
//???		printf("LamAny: %d\tNTrackBox: %2.f\ti: %d\n",LamAny,NTrackBox,i+1);
		if(LamAny)
		{
			NLamBox++;	//���� ���� LAM ��� ������� ����� ���������� ����� ������� �� �� �������������� ����� ����� LAM. ��������� 25.06.12 13:00
//???			printf("NLamBox: %2.f\ti: %d\n",NLamBox,i+1);
		}
                NDataBox +=idBM;
		if(idimp)
		{
			NAmpPorogBox++;
			if(LamAny)NLamAmpPorogBox++;
//			else printf("NLamBox: %2.f\tLamAny: %d\tNAmpPorogBox: %2.f\ti: %d\n",NLamBox,LamAny,NAmpPorogBox,i+1);
//			if(LamAny)NLamAmpPorogBox++;	//���� ���� LAM ��� ������������ �� �� ����� ���������� ����� ������� �� �� �������������� ����� ����� LAM. ��������� 25.06.12 13:00
		}


		for(int chan=0; chan < NCHANADC; chan++)
		{
		  int ipmt=ADCtoPMT[chan]-1; //����� ������ ���
		  int istrip=TD.ChanStrip[ipmt];
		  double amp=-(TD.Data[chan]-TD.PDS[chan]);
		  if(amp < Par.porogImpulse[1])continue; //��������� ���������
		  strcpy(NameFileWrite,DirWrite);
		  strcat(NameFileWrite,"\\");
		  sprintf(NameFileWrite+strlen(NameFileWrite),"XYZstrip_%02d.dat",istrip+1);
		  int hxyz=open(NameFileWrite,O_RDWR|O_CREAT|O_TEXT,S_IWRITE);
          if(lseek(hxyz,0,SEEK_END)==0)
          {
			sprintf(buf,"X\tY\tZ\tAmp\n");
            write(hxyz,buf,strlen(buf));
          }
		  sprintf(buf,"%lf\t%lf\t%lf\t%.2lf\n",(p0[0]+p1[0])/2.,(p0[1]+p1[1])/2.,
                      (p0[2]+p1[2])/2.,amp);
          write(hxyz,buf,strlen(buf));
          close(hxyz);
        }//chan
       }
	  }//ChanProc < 0
     }//Track
     if(kpointC >= 2)
     {
       cross1.x=p0[0]; cross1.y=p0[1]; cross1.z=p0[2];
	   cross2.x=p1[0]; cross2.y=p1[1]; cross2.z=p1[2];
       L=(cross2.x-cross1.x)*(cross2.x-cross1.x)+
         (cross2.y-cross1.y)*(cross2.y-cross1.y)+
         (cross2.z-cross1.z)*(cross2.z-cross1.z);
     }
     if(L)L=sqrt(L);
     int maska=0;
	 if( (TD.dtsd.id_trek&0x0FF) < 2 && ThetaG >=0.)
     {
       maska++;
       if(L > 0.)
       {
        maska++;
        if(idimp)
           maska++;
	   }
	 }
/*
     sprintf(buf,"%ld\t%ld\t%ld",TD.dtsd.Nrun,TD.dtsd.Nevent,TD.dtsd.id_trek&0x0FF);
     for(int ii=0; ii <= maska; ii++)
	   write(hand_out[ii],buf,strlen(buf));

     if(ThetaG >=0.)
      sprintf(buf,"\t%.2lf\t%.2lf",PhiG,ThetaG);
     else
      sprintf(buf,"\t-\t-");
     for(int ii=0; ii <= maska; ii++)
       write(hand_out[ii],buf,strlen(buf));

	 sprintf(buf,"\t%.2lf\t%.2lf\t%.2lf",L,Res.pds,Res.sigmapds);
     for(int ii=0; ii <= maska; ii++)
	   write(hand_out[ii],buf,strlen(buf));
    if(idimp)
	 sprintf(buf,"\t%.2lf\t%.2lf\n",Res.amppik,Res.T);
    else
     sprintf(buf,"\t-\t-\n");
    for(int ii=0; ii <= maska; ii++)
	  write(hand_out[ii],buf,strlen(buf));
*/
   int lenrad=-1;
    if(maska >=2)
    {
      lenrad=L/10.;
      if(lenrad >=0 && lenrad < 9)
      {
		nDetAll[lenrad]++;
		SredLenAll[lenrad] +=L;
        SigLenAll[lenrad] +=L*L;
      }
    }
    if(maska==3 && ChanProc >=0)
	{
      int len=lenrad;
      if(len >=0 && len < 9)
      {
        double amp=Res.amppik;
         SredLen[len] +=L;
         SigLen[len] +=L*L;
		 AmpDet[len] +=amp;
         SigDet[len] +=amp*amp;
		 nDet[len]++;
         sprintf(buf,"%ld\t%ld\t%ld",TD.dtsd.Nrun,TD.dtsd.Nevent,TD.dtsd.id_trek&0x0FF);
		 sprintf(buf+strlen(buf),"\t%.2lf\t%.2lf",PhiG,ThetaG);
         sprintf(buf+strlen(buf),"\t%.2lf\t%.2lf\t%.2lf",L,Res.pds,Res.sigmapds);
         sprintf(buf+strlen(buf),"\t%.2lf\t%.2lf\n",Res.amppik,Res.T);
         if(len==1)
          write(hamp[0],buf,strlen(buf));
         if(len==4)
          write(hamp[1],buf,strlen(buf));
         if(len==7)
          write(hamp[2],buf,strlen(buf));
      }
    }
//    if(idimp && Res.amppik < -4.e-3 && (ThetaG >= 0. && ThetaG < 10.))
    if(idimp && (ThetaG >= 0. && ThetaG < 15.))
	{
	  //����� ����������� ������ �������������� ���������
	  for(int ipl=0; ipl < NplMatrix; ipl++)
          {
		double h=hMatrix[ipl];
                double z=Ce1.z+h;
                double t=(z-P1.z)*iz/dz;
                double x=P1.x+iz*dx*t-Ce1.x;
		double y=P1.y+iz*dy*t-Ce1.y;
		int ix=(x/D_X+NM_X/2)+0.5;
		int iy=(y/D_Y+NM_Y/2)+0.5;
		if(ix >=0 && ix < NM_X && iy >=0 && iy < NM_Y)
		{
		  nMatrixXY[ipl][iy][ix]++;
		  MatrixXY[ipl][iy][ix] +=(Res.amppik);
		}
	  }
	}
        if(Res.lam && (ThetaG >= 0. && ThetaG < 15.))
	{
	  //����� ����������� ������ �������������� ���������
	  for(int ipl=0; ipl < NplMatrix; ipl++)
          {
		double h=hMatrix[ipl];
                double z=Ce1.z+h;
                double t=(z-P1.z)*iz/dz;
                double x=P1.x+iz*dx*t-Ce1.x;
		double y=P1.y+iz*dy*t-Ce1.y;
		int ix=(x/D_X+NM_X/2)+0.5;
		int iy=(y/D_Y+NM_Y/2)+0.5;
		if(ix >=0 && ix < NM_X && iy >=0 && iy < NM_Y)
		{
		  LnMatrixXY[ipl][iy][ix]++;
		}
	  }
	}
        if(idBM==0 && (ThetaG >= 0. && ThetaG < 15.))
	{
	  //����� ����������� ������ �������������� ���������
	  for(int ipl=0; ipl < NplMatrix; ipl++)
          {
		double h=hMatrix[ipl];
                double z=Ce1.z+h;
                double t=(z-P1.z)*iz/dz;
                double x=P1.x+iz*dx*t-Ce1.x;
		double y=P1.y+iz*dy*t-Ce1.y;
		int ix=(x/D_X+NM_X/2)+0.5;
		int iy=(y/D_Y+NM_Y/2)+0.5;
		if(ix >=0 && ix < NM_X && iy >=0 && iy < NM_Y)
		{
		  noBMnMatrixXY[ipl][iy][ix]++;
		}
	  }
	}
	if((ThetaG >= 0. && ThetaG < 15.))
	{
	  //����� ����������� ������ �������������� ���������
	  for(int ipl=0; ipl < NplMatrix; ipl++)
	  {
		double h=hMatrix[ipl];
		double z=Ce1.z+h;
		double t=(z-P1.z)*iz/dz;
		double x=P1.x+iz*dx*t-Ce1.x;
		double y=P1.y+iz*dy*t-Ce1.y;
		int ix=(x/D_X+NM_X/2)+0.5;
		int iy=(y/D_Y+NM_Y/2)+0.5;
		if(ix >=0 && ix < NM_X && iy >=0 && iy < NM_Y)
		{
		  tMatrixXY[ipl][iy][ix]++;
		}
	  }
	}
//��� ������ � ������� ���������
     if(kpointCB >=2 || kpointC >= 2)
     {
       cross1.x=p0[0]; cross1.y=p0[1]; cross1.z=p0[2];
	   cross2.x=p1[0]; cross2.y=p1[1]; cross2.z=p1[2];
       L=(cross2.x-cross1.x)*(cross2.x-cross1.x)+
         (cross2.y-cross1.y)*(cross2.y-cross1.y)+
         (cross2.z-cross1.z)*(cross2.z-cross1.z);
     }
     if(L)L=sqrt(L);
     maska=0;
	 if( (TD.dtsd.id_trek&0x0FF) < 2 && ThetaG >=0.)
     {
       maska++;
       if(L > 0.)
       {
        maska++;
        if(idimp)
           maska++;
	   }
	 }

    if(maska >=3 && idimp && (ThetaG >= 0. && ThetaG < 15.))
	{
	  //����� ����������� ������ �������������� ���������
	  for(int ipl=0; ipl < NplMatrix; ipl++)
          {
		double h=hMatrix[ipl];
                double z=Ce1.z+h;
                double t=(z-P1.z)*iz/dz;
                double x=P1.x+iz*dx*t-Ce1.x;
		double y=P1.y+iz*dy*t-Ce1.y;
		int ix=(x/D_X+NM_X/2)+0.5;
		int iy=(y/D_Y+NM_Y/2)+0.5;
		if(ix >=0 && ix < NM_X && iy >=0 && iy < NM_Y)
		{
		  RnMatrixXY[ipl][iy][ix]++;
		  RMatrixXY[ipl][iy][ix] +=(Res.amppik);
		}
	  }
	}
    if(maska >=2 && Res.lam && (ThetaG >= 0. && ThetaG < 15.))
	{
	  //����� ����������� ������ �������������� ���������
	  for(int ipl=0; ipl < NplMatrix; ipl++)
          {
		double h=hMatrix[ipl];
                double z=Ce1.z+h;
                double t=(z-P1.z)*iz/dz;
                double x=P1.x+iz*dx*t-Ce1.x;
		double y=P1.y+iz*dy*t-Ce1.y;
		int ix=(x/D_X+NM_X/2)+0.5;
		int iy=(y/D_Y+NM_Y/2)+0.5;
		if(ix >=0 && ix < NM_X && iy >=0 && iy < NM_Y)
		{
		  RLnMatrixXY[ipl][iy][ix]++;
		}
	  }
	}
    if(maska >=2 && idBM==0 && (ThetaG >= 0. && ThetaG < 15.))
	{
	  //����� ����������� ������ �������������� ���������
	  for(int ipl=0; ipl < NplMatrix; ipl++)
          {
		double h=hMatrix[ipl];
                double z=Ce1.z+h;
                double t=(z-P1.z)*iz/dz;
                double x=P1.x+iz*dx*t-Ce1.x;
		double y=P1.y+iz*dy*t-Ce1.y;
		int ix=(x/D_X+NM_X/2)+0.5;
		int iy=(y/D_Y+NM_Y/2)+0.5;
		if(ix >=0 && ix < NM_X && iy >=0 && iy < NM_Y)
		{
		  RnoBMnMatrixXY[ipl][iy][ix]++;
		}
	  }
	}

	if(maska >=2 && (ThetaG >= 0. && ThetaG < 15.))
	{
	  //����� ����������� ������ �������������� ���������
	  for(int ipl=0; ipl < NplMatrix; ipl++)
	  {
		double h=hMatrix[ipl];
		double z=Ce1.z+h;
		double t=(z-P1.z)*iz/dz;
		double x=P1.x+iz*dx*t-Ce1.x;
		double y=P1.y+iz*dy*t-Ce1.y;
		int ix=(x/D_X+NM_X/2)+0.5;
		int iy=(y/D_Y+NM_Y/2)+0.5;
		if(ix >=0 && ix < NM_X && iy >=0 && iy < NM_Y)
		{
		  RtMatrixXY[ipl][iy][ix]++;
		}
	  }
	}

   }
 }while( (done=findnext(&ffile)) ==0 );
 for(int ii=0; ii < 4; ii++)
   close(hand_out[ii]);
 for(int ii=0; ii < 3; ii++)
   close(hamp[ii]);
 findclose(&ffile);
 return 0;

}