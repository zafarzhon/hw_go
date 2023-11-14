#ifndef UTILURG_H
#define UTILURG_H
struct TIMEOUT_S
 {
	 unsigned short maska_timeout; // маска используемых каналов
	 unsigned short flag_timeout; // маска флагов бит=1:таймаут истек
	 unsigned long  count[8]; // конец таймаутов в тиках
	 unsigned long  ticker[8]; // счетчик таймаутов в тиках
 };
extern struct TIMEOUT_S *timeout;
extern unsigned short *flag_timeout; // для UTILURG.CPP
extern short USE_TIMEOUT; // маска включенных счетчиков таймаутов
extern unsigned short kb; // код нажатой клавиши
extern unsigned long Nrun; // номер рана
extern unsigned long Nevent; // номер события
extern unsigned long Nsob; // Количество событий в ране

#define NCNTR 2  // кол-во контроллеров в одной ЭВМ
#define NCHAN 4  // кол-во каналов в одном контроллере
#define NGROUP8 76 // максимально возможное кол-во групп по 8 бит
#define NGROUP16 38  // в полной конфигурации 38 групп по 16 бит
#define NGROUP64 10 // максимально возможное кол-во групп по 64 бита
#define NPOINTMAX 320 // максимально возможное битов в X-ах или в Y-ах
#define MAXXSTRIP 320
#define MAXYSTRIP 288

// Структура конфигурации одного канала
struct CONFIG_ONE
{
	 short include; // 0- не в работе , 1- в работе
	 short Nx; // кол-во бит в камерах X
	 short Ny; // кол-во бит в камерах Y
	 float Xx,Yx,Zx; // координаты первого бита данных в камерах X
	 float Xy,Yy,Zy; // координаты первого бита данных в камерах Y
	 float VXx,VYx,VZx; // Направляющий вектор в камерах X
	 float VXy,VYy,VZy; // Направляющий вектор в камерах Y
};
// Структура конфигурации всего ДЕКОРа
struct CONFIG_DECOR
{
	  struct CONFIG_ONE cnf[128];
};
extern struct CONFIG_DECOR *CNF; // Структура конфигурации всего ДЕКОРа

// Статистика срабатываний
#pragma pack(push,1)
struct STATISTIC {
						 unsigned long nlam[NCHAN];
						 unsigned long nhit[NCHAN][NGROUP8*8];
					  };
struct STATKADR {
					 short Sm;//supermodul
					 unsigned long Nrun; //Run
					 unsigned long Nkadr; //Kadr
  unsigned long timestart; //Время начала заполнения кадра
  unsigned long datestart;
  unsigned long timestop; //Время конца заполнения кадра
  unsigned long datestop;
					 STATISTIC Hit[NCNTR];
	        };
#pragma pack(pop)
// статистика срабатываний каналов для рана
extern struct STATISTIC *Hitstat[NCNTR];
// статистика срабатываний каналов отдельно при наборе одного кадра
extern struct STATKADR *StatKadr;

#pragma pack(push,1)
struct EVENT  {
					 unsigned long st;    //номер события
					 unsigned char tr[2]; // статус биты[local,extern], маска сработавших плоскостей
					 unsigned short ti[2]; // живое время время ожидания последнего события
					 unsigned long dt[64];// закодированные данные
				  };
#pragma pack(pop)

extern void Init_Stat();
extern void StartWork();
extern void ExitWork();
//---- Мониторирование -----
extern int Monit();
//--- Измерение шумов плоскостей -----
extern int Noise();
//------Нейтральный режим -----------
extern void Pusto();
//------- Режим экспозиции ------
extern int Expoz();

extern void Write_logfile(char *reg,char *mes); // запись в log-file

#endif
