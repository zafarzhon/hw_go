#ifndef UTILGRAPH_H
#define UTILGRAPH_H
#include <windows.h>
#include <math.h>
//Точка в трёхмерном пространстве
#ifndef TPoint3D_D
#define TPoint3D_D
struct TPoint3D
{
  double x;
  double y;
  double z;
  TPoint3D()
  {
   x=0; y=0; z=0;
  }
  TPoint3D(double a, double b, double c)
  {
   x=a; y=b; z=c;
  }
};
#endif

struct TEllipse3D {
 //Начальные параметры
  TPoint3D center; //центр в исходной системе координат
  double radius; //Радиус
  double nphi,ntheta; //Направление вектора нормали в градусах
  //Исходные контрольные точки
  TPoint3D p[13];
  //Контрольные точки в системе координат отображения
  TPoint3D vp[13];
  //Окончательные точки для вызова PolyBezier
  POINT P[13];
  //Поворот плоскости точек p, перпендикулярно нормали
  void Rotate();
  //Преобразование начальных параметров
  // в исходные контрольные точки
  void ToBezier();
  //Конструктор
  TEllipse3D()
  {
    center=TPoint3D(0.,0.,0.);
    radius=0.;
    nphi=0.;
    ntheta=0.;
    ToBezier();
  }
  //Формирующий конструктор
  TEllipse3D(TPoint3D c,double r,double phi,double theta)
  {
    center=c;
    radius=r;
    nphi=phi;
    ntheta=theta;
    ToBezier();
  }
  void SetP()
  {
    for(int i=0; i < 13; i++)
    {
      P[i].x=vp[i].x;
      P[i].y=vp[i].y;
    }
  }
  //Для ортогональной проекции
  void SetP(double sx,double cx,double sy,double cy)
  {
    for(int i=0; i < 13; i++)
    {
      P[i].x=vp[i].x*sx+cx;
      P[i].y=vp[i].y*sy+cy;
    }
  }
  //Для перспективной проекции
  void SetP(double h,double xv,double yv,double sx,double sy)
  {
    for(int i=0; i < 13; i++)
    {
      P[i].x=(vp[i].x*h/vp[i].z+xv)*sx;
      P[i].y=(vp[i].y*h/vp[i].z+yv)*sy;
    }
  }
};
//Векторное произведение двух векторов
//Возвращает длину полученного вектора
extern double Vmul(double *v1,double *v2,double *v);//v=[v1 x v2]
//Нормализация вектора
extern double Vnorm(double *v);

//Нахождение максимального сближения двух отрезков
// На входе:
// p0[2] - первая и вторая точки первого отрезка
// p1[2] - первая и вторая точки второго отрезка
//На выходе:
// pp[2] - точки максимального сближения на первом и втором отрезке
// Возвращает:
//        минимальное расстояние между отрезками
extern double GetSectionDistance(TPoint3D *p0,TPoint3D *p1,TPoint3D *pp);
//Поиск точек пересечения прямой с поверхностью конуса
// Конус задаётся:
// точками центров окружностей оснований Ce1 и Ce2
// радиусами Cp1.x и Cp2.x
// нормалями в градусах phi: Cp1.y и Cp2.y
//                    theta: Cp1.z и Cp2.z
// Прямая задаётся двумя точками P1 и P2
// bStoreT = true - сохранение значений функционала при поиске в файле
// Функция возвращает:
//    0 - точек не найдено
//    1 - найдена одна точка cross1
//    2 - найдено две точки  cross1 и cross2
extern int CroosLineCone(TPoint3D Ce1, TPoint3D Cp1, TPoint3D Ce2, TPoint3D Cp2,
                  TPoint3D P1, TPoint3D P2, TPoint3D &cross1, TPoint3D &cross2,
                  bool bStoreT);
#endif