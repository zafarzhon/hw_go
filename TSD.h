#ifndef TSD_H
#define TSD_H

#pragma pack(push,1)
struct TEVENT  {
 unsigned long st;    //номер события
 unsigned char tr[2]; // статус биты[local,extern], маска сработавших плоскостей
 unsigned short ti[2]; // живое время время ожидания последнего события
 unsigned long dt[64];// закодированные данные
};
#pragma pack(pop)
#pragma pack(push,1)
#ifndef dTDateTimeKadr
#define dTDateTimeKadr
//структура для хранения даты и времени в виде удобном для
//сравнения по критериям раньше, позже
struct TDateTimeKadr
{
 union{
  struct{
  unsigned char hsecond; //сотые секунды
  unsigned char second;  //секунда
  unsigned char minute;  //минута
  unsigned char hour; //час
  }tm;
  unsigned long time;
 };
 union{
  struct{
  unsigned char day; //день
  unsigned char month; //месяц
  unsigned short year; //год
  }dt;
  unsigned long date;
 };
};
#endif
#pragma pack(pop)

#pragma pack(push,1)
struct DataTSD //Структура метки события для TSD
{
      char start[6];  // Ключевое слово начала записи
      short Tip;  // Тип записи:0-Config,1-монитор,2-Experement event,3-Noise
      unsigned long Nrun; // Номер текущего рана
      unsigned long Nevent; // Номер текущего события
      TDateTimeKadr dt; //Время события
      short id_trek; //Признак наличия трека в URAGANе
      //Координаты трека в миллиметрах в системе НЕВОДа
      double chp0[3];
      double chp1[3];
      TEVENT ev; //Кодированные данные события
};
#pragma pack(pop)

#define NCHANADC 32
#define NCHANOSC 4

#pragma pack(push,1)
struct CMarocSetting {
	char SerialNumber[8]; // Серийный номер Айзека
	short ID; // Идентефикатор подключенного устройства
	unsigned char Word[12]; // Первые 12 регистров
	unsigned char Gain[NCHANADC]; // 64 усилителя
	unsigned char Dac[12]; // 12 регистров ЦАП
        unsigned char Word93,Word94;
	short MessureType; // 0 - подсветка
	short DataType; // 0 - если смотрим триггеры, 1-если АЦП
	short ReadTO; // Таймаут чтения
	short WriteTO; // Таймаут записи
};
#pragma pack(pop)

#pragma pack(push,1)
struct T_TSD
{
  DataTSD dtsd;//Структура метки события для TSD от URAGANа
  CMarocSetting MS;//Конфигурация ASIC
  short timepds;//Интервал мремени между измерениями пьедесталов в минутах
  short timenoise;//Интервал мремени между измерениями шумов в минутах
  short nzappds; //Количество запусков для пьедестала
  char X0_BM[2]; //Точка X привязки координат BM и УРАГАН [карта, бит]
  char Y0_BM[2]; //Точка Y привязки координат BM и УРАГАН [карта, бит]
  short ChanStrip[NCHANADC];//Таблица соответствия канал АЦП -> номер стрипа
  short nchan; //Кoличество каналов АЦП = NCHANADC
  char flagmon[2]; //[pds, noise] =1 - новое значение
  float noisetime;// время измерения шумов, сек
  float PDS[NCHANADC]; //Пьедесталы, каналы АЦП
  float sigPDS[NCHANADC]; //Сигма пьедесталов, каналы АЦП
  long Noise[NCHANADC]; //Счёт каналов за время noisetime
  unsigned char Lam[NCHANADC];// LAM
  short Data[NCHANADC];//Отклик на событие каналы АЦП
  short Tempra[2]; //Температура по внутреннему датчику PETIROC
};
#pragma pack(pop)

//Параметры обработки
struct PAR_OBR
{
  int inverse; //1- отрицательный сигнал, 0 - положительный сигнал
  double strob[2]; //Границы строба, области поиска полезного сигнала
  double R; //Нагрузочное сопротивление, Ом
  double porogImpulse[2]; //Двухсторонний порог выделения полезного
                          // импульса
  double porogPDS[2]; //Двухсторонний порог определения пьедестала
};

struct RES_OBRABTSD
{
  double pds; //Пьедестал в мВ
  double sigmapds; //Сигма пьедестала в мВ
  double amppik; //Пиковая амплитуда в мВ за вычетом пьедестала
  double timp; //Время появления импульса (превышение порога 3sigma(pds) ), нсек
  double Q; //Заряд импульса внутри строба за вычетом пьедестала, пикоКулонах
  double T; //длительность импульса (превышение порога 3sigma(pds) ), нсек
  double Noise;//Частота шумов, сек^-1
  long lam; //=1 -Наличие сигнала LAM
};

extern CMarocSetting *Device;

extern bool BM_Open(); // Поиск устройства, подключение и начальные установки из файла settings/default.cfg
extern bool BM_Close(); // Отключение устройства
extern bool BM_Upload(const char *FileName); // Загрузка в MAROC2 настроек из файла
extern bool BM_ExtractData(short *Data, char *Lam); // Чтение данных.
// Возвращает true в случае удачного чтения данных и false с очищенным массивом в противном случае .
// Объем Data лучше не делать меньше 64 элементов.
extern bool BM_Noise(long *Mean); // Поканальный подсчет шумов. Mean[64] - выходной счет по каналам
extern bool BM_Pedestal(float *Mean, float *Sigma, long Stat); // Подсчет пъедестала. Mean[64] средняя амплитуда для каждого канала.
// Stat[64] сигма для каждого канала. Stat - количество запусков АЦП
extern bool BM_FIFOClear(); // Очистка FIFO

//Обработка данных осциллографа
// chan - номер канала
// TSD - исходные данные осциллографа
// Par - параметры обработки
// На выходе:
// Возвращаемое значение: 1- обработка прошла успешно, 0 - неудача
// Res - результаты обработки
extern int Obrab_TSD(int chan, T_TSD &TSD, PAR_OBR &Par, RES_OBRABTSD &Res);

#endif
