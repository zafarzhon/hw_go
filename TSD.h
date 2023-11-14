#ifndef TSD_H
#define TSD_H

#pragma pack(push,1)
struct TEVENT  {
 unsigned long st;    //����� �������
 unsigned char tr[2]; // ������ ����[local,extern], ����� ����������� ����������
 unsigned short ti[2]; // ����� ����� ����� �������� ���������� �������
 unsigned long dt[64];// �������������� ������
};
#pragma pack(pop)
#pragma pack(push,1)
#ifndef dTDateTimeKadr
#define dTDateTimeKadr
//��������� ��� �������� ���� � ������� � ���� ������� ���
//��������� �� ��������� ������, �����
struct TDateTimeKadr
{
 union{
  struct{
  unsigned char hsecond; //����� �������
  unsigned char second;  //�������
  unsigned char minute;  //������
  unsigned char hour; //���
  }tm;
  unsigned long time;
 };
 union{
  struct{
  unsigned char day; //����
  unsigned char month; //�����
  unsigned short year; //���
  }dt;
  unsigned long date;
 };
};
#endif
#pragma pack(pop)

#pragma pack(push,1)
struct DataTSD //��������� ����� ������� ��� TSD
{
      char start[6];  // �������� ����� ������ ������
      short Tip;  // ��� ������:0-Config,1-�������,2-Experement event,3-Noise
      unsigned long Nrun; // ����� �������� ����
      unsigned long Nevent; // ����� �������� �������
      TDateTimeKadr dt; //����� �������
      short id_trek; //������� ������� ����� � URAGAN�
      //���������� ����� � ����������� � ������� ������
      double chp0[3];
      double chp1[3];
      TEVENT ev; //������������ ������ �������
};
#pragma pack(pop)

#define NCHANADC 32
#define NCHANOSC 4

#pragma pack(push,1)
struct CMarocSetting {
	char SerialNumber[8]; // �������� ����� ������
	short ID; // ������������� ������������� ����������
	unsigned char Word[12]; // ������ 12 ���������
	unsigned char Gain[NCHANADC]; // 64 ���������
	unsigned char Dac[12]; // 12 ��������� ���
        unsigned char Word93,Word94;
	short MessureType; // 0 - ���������
	short DataType; // 0 - ���� ������� ��������, 1-���� ���
	short ReadTO; // ������� ������
	short WriteTO; // ������� ������
};
#pragma pack(pop)

#pragma pack(push,1)
struct T_TSD
{
  DataTSD dtsd;//��������� ����� ������� ��� TSD �� URAGAN�
  CMarocSetting MS;//������������ ASIC
  short timepds;//�������� ������� ����� ����������� ����������� � �������
  short timenoise;//�������� ������� ����� ����������� ����� � �������
  short nzappds; //���������� �������� ��� ����������
  char X0_BM[2]; //����� X �������� ��������� BM � ������ [�����, ���]
  char Y0_BM[2]; //����� Y �������� ��������� BM � ������ [�����, ���]
  short ChanStrip[NCHANADC];//������� ������������ ����� ��� -> ����� ������
  short nchan; //�o�������� ������� ��� = NCHANADC
  char flagmon[2]; //[pds, noise] =1 - ����� ��������
  float noisetime;// ����� ��������� �����, ���
  float PDS[NCHANADC]; //����������, ������ ���
  float sigPDS[NCHANADC]; //����� �����������, ������ ���
  long Noise[NCHANADC]; //���� ������� �� ����� noisetime
  unsigned char Lam[NCHANADC];// LAM
  short Data[NCHANADC];//������ �� ������� ������ ���
  short Tempra[2]; //����������� �� ����������� ������� PETIROC
};
#pragma pack(pop)

//��������� ���������
struct PAR_OBR
{
  int inverse; //1- ������������� ������, 0 - ������������� ������
  double strob[2]; //������� ������, ������� ������ ��������� �������
  double R; //����������� �������������, ��
  double porogImpulse[2]; //������������� ����� ��������� ���������
                          // ��������
  double porogPDS[2]; //������������� ����� ����������� ����������
};

struct RES_OBRABTSD
{
  double pds; //��������� � ��
  double sigmapds; //����� ���������� � ��
  double amppik; //������� ��������� � �� �� ������� ����������
  double timp; //����� ��������� �������� (���������� ������ 3sigma(pds) ), ����
  double Q; //����� �������� ������ ������ �� ������� ����������, �����������
  double T; //������������ �������� (���������� ������ 3sigma(pds) ), ����
  double Noise;//������� �����, ���^-1
  long lam; //=1 -������� ������� LAM
};

extern CMarocSetting *Device;

extern bool BM_Open(); // ����� ����������, ����������� � ��������� ��������� �� ����� settings/default.cfg
extern bool BM_Close(); // ���������� ����������
extern bool BM_Upload(const char *FileName); // �������� � MAROC2 �������� �� �����
extern bool BM_ExtractData(short *Data, char *Lam); // ������ ������.
// ���������� true � ������ �������� ������ ������ � false � ��������� �������� � ��������� ������ .
// ����� Data ����� �� ������ ������ 64 ���������.
extern bool BM_Noise(long *Mean); // ����������� ������� �����. Mean[64] - �������� ���� �� �������
extern bool BM_Pedestal(float *Mean, float *Sigma, long Stat); // ������� ����������. Mean[64] ������� ��������� ��� ������� ������.
// Stat[64] ����� ��� ������� ������. Stat - ���������� �������� ���
extern bool BM_FIFOClear(); // ������� FIFO

//��������� ������ ������������
// chan - ����� ������
// TSD - �������� ������ ������������
// Par - ��������� ���������
// �� ������:
// ������������ ��������: 1- ��������� ������ �������, 0 - �������
// Res - ���������� ���������
extern int Obrab_TSD(int chan, T_TSD &TSD, PAR_OBR &Par, RES_OBRABTSD &Res);

#endif
