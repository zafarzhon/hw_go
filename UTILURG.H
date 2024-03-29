#ifndef UTILURG_H
#define UTILURG_H
struct TIMEOUT_S
 {
	 unsigned short maska_timeout; // ��᪠ �ᯮ��㥬�� �������
	 unsigned short flag_timeout; // ��᪠ 䫠��� ���=1:⠩���� ��⥪
	 unsigned long  count[8]; // ����� ⠩���⮢ � ⨪��
	 unsigned long  ticker[8]; // ���稪 ⠩���⮢ � ⨪��
 };
extern struct TIMEOUT_S *timeout;
extern unsigned short *flag_timeout; // ��� UTILURG.CPP
extern short USE_TIMEOUT; // ��᪠ ����祭��� ���稪�� ⠩���⮢
extern unsigned short kb; // ��� ����⮩ ������
extern unsigned long Nrun; // ����� ࠭�
extern unsigned long Nevent; // ����� ᮡ���
extern unsigned long Nsob; // ������⢮ ᮡ�⨩ � ࠭�

#define NCNTR 2  // ���-�� ����஫��஢ � ����� ���
#define NCHAN 4  // ���-�� ������� � ����� ����஫���
#define NGROUP8 76 // ���ᨬ��쭮 ��������� ���-�� ��㯯 �� 8 ���
#define NGROUP16 38  // � ������ ���䨣��樨 38 ��㯯 �� 16 ���
#define NGROUP64 10 // ���ᨬ��쭮 ��������� ���-�� ��㯯 �� 64 ���
#define NPOINTMAX 320 // ���ᨬ��쭮 ��������� ��⮢ � X-�� ��� � Y-��
#define MAXXSTRIP 320
#define MAXYSTRIP 288

// ������� ���䨣��樨 ������ ������
struct CONFIG_ONE
{
	 short include; // 0- �� � ࠡ�� , 1- � ࠡ��
	 short Nx; // ���-�� ��� � ������ X
	 short Ny; // ���-�� ��� � ������ Y
	 float Xx,Yx,Zx; // ���न���� ��ࢮ�� ��� ������ � ������ X
	 float Xy,Yy,Zy; // ���न���� ��ࢮ�� ��� ������ � ������ Y
	 float VXx,VYx,VZx; // ���ࠢ���騩 ����� � ������ X
	 float VXy,VYy,VZy; // ���ࠢ���騩 ����� � ������ Y
};
// ������� ���䨣��樨 �ᥣ� ������
struct CONFIG_DECOR
{
	  struct CONFIG_ONE cnf[128];
};
extern struct CONFIG_DECOR *CNF; // ������� ���䨣��樨 �ᥣ� ������

// ����⨪� �ࠡ��뢠���
#pragma pack(push,1)
struct STATISTIC {
						 unsigned long nlam[NCHAN];
						 unsigned long nhit[NCHAN][NGROUP8*8];
					  };
struct STATKADR {
					 short Sm;//supermodul
					 unsigned long Nrun; //Run
					 unsigned long Nkadr; //Kadr
  unsigned long timestart; //�६� ��砫� ���������� ����
  unsigned long datestart;
  unsigned long timestop; //�६� ���� ���������� ����
  unsigned long datestop;
					 STATISTIC Hit[NCNTR];
	        };
#pragma pack(pop)
// ����⨪� �ࠡ��뢠��� ������� ��� ࠭�
extern struct STATISTIC *Hitstat[NCNTR];
// ����⨪� �ࠡ��뢠��� ������� �⤥�쭮 �� ����� ������ ����
extern struct STATKADR *StatKadr;

#pragma pack(push,1)
struct EVENT  {
					 unsigned long st;    //����� ᮡ���
					 unsigned char tr[2]; // ����� ����[local,extern], ��᪠ �ࠡ�⠢�� ���᪮�⥩
					 unsigned short ti[2]; // ����� �६� �६� �������� ��᫥����� ᮡ���
					 unsigned long dt[64];// ������஢���� �����
				  };
#pragma pack(pop)

extern void Init_Stat();
extern void StartWork();
extern void ExitWork();
//---- ������஢���� -----
extern int Monit();
//--- ����७�� �㬮� ���᪮�⥩ -----
extern int Noise();
//------����ࠫ�� ०�� -----------
extern void Pusto();
//------- ����� ��ᯮ��樨 ------
extern int Expoz();

extern void Write_logfile(char *reg,char *mes); // ������ � log-file

#endif
