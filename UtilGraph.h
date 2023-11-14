#ifndef UTILGRAPH_H
#define UTILGRAPH_H
#include <windows.h>
#include <math.h>
//����� � ��������� ������������
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
 //��������� ���������
  TPoint3D center; //����� � �������� ������� ���������
  double radius; //������
  double nphi,ntheta; //����������� ������� ������� � ��������
  //�������� ����������� �����
  TPoint3D p[13];
  //����������� ����� � ������� ��������� �����������
  TPoint3D vp[13];
  //������������� ����� ��� ������ PolyBezier
  POINT P[13];
  //������� ��������� ����� p, ��������������� �������
  void Rotate();
  //�������������� ��������� ����������
  // � �������� ����������� �����
  void ToBezier();
  //�����������
  TEllipse3D()
  {
    center=TPoint3D(0.,0.,0.);
    radius=0.;
    nphi=0.;
    ntheta=0.;
    ToBezier();
  }
  //����������� �����������
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
  //��� ������������� ��������
  void SetP(double sx,double cx,double sy,double cy)
  {
    for(int i=0; i < 13; i++)
    {
      P[i].x=vp[i].x*sx+cx;
      P[i].y=vp[i].y*sy+cy;
    }
  }
  //��� ������������� ��������
  void SetP(double h,double xv,double yv,double sx,double sy)
  {
    for(int i=0; i < 13; i++)
    {
      P[i].x=(vp[i].x*h/vp[i].z+xv)*sx;
      P[i].y=(vp[i].y*h/vp[i].z+yv)*sy;
    }
  }
};
//��������� ������������ ���� ��������
//���������� ����� ����������� �������
extern double Vmul(double *v1,double *v2,double *v);//v=[v1 x v2]
//������������ �������
extern double Vnorm(double *v);

//���������� ������������� ��������� ���� ��������
// �� �����:
// p0[2] - ������ � ������ ����� ������� �������
// p1[2] - ������ � ������ ����� ������� �������
//�� ������:
// pp[2] - ����� ������������� ��������� �� ������ � ������ �������
// ����������:
//        ����������� ���������� ����� ���������
extern double GetSectionDistance(TPoint3D *p0,TPoint3D *p1,TPoint3D *pp);
//����� ����� ����������� ������ � ������������ ������
// ����� �������:
// ������� ������� ����������� ��������� Ce1 � Ce2
// ��������� Cp1.x � Cp2.x
// ��������� � �������� phi: Cp1.y � Cp2.y
//                    theta: Cp1.z � Cp2.z
// ������ ������� ����� ������� P1 � P2
// bStoreT = true - ���������� �������� ����������� ��� ������ � �����
// ������� ����������:
//    0 - ����� �� �������
//    1 - ������� ���� ����� cross1
//    2 - ������� ��� �����  cross1 � cross2
extern int CroosLineCone(TPoint3D Ce1, TPoint3D Cp1, TPoint3D Ce2, TPoint3D Cp2,
                  TPoint3D P1, TPoint3D P2, TPoint3D &cross1, TPoint3D &cross2,
                  bool bStoreT);
#endif