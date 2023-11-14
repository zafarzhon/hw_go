#include "UtilGraph.h"
#include <stdio.h>
#include <fcntl.h>
#include <io.h>
#include <sys\stat.h>

//��������� ������������ ���� ��������
//���������� ����� ����������� �������
double Vmul(double *v1,double *v2,double *v)//v=[v1 x v2]
{
  double r=0;
  v[0]=v1[1]*v2[2]-v1[2]*v2[1];
  v[1]=v1[2]*v2[0]-v1[0]*v2[2];
  v[2]=v1[0]*v2[1]-v1[1]*v2[0];
  r=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  return r;
}
double Vnorm(double *v)//������������ �������
{
  int i;
  double d=0.;
  for(i=0; i < 3; i++)d+=v[i]*v[i];
  if(d==0.)return d;
  d=sqrt(d);
  for(i=0; i < 3; i++)v[i]/=d;
  return d;
}
//������� ��������� ����� p, ��������������� �������
void TEllipse3D::Rotate()
{
  TPoint3D temp[13]; //��������� ������
  double v11,v21,v31,v12,v22,v32,v13,v23,v33;
  double v43;
  double vX[3],vY[3],vZ[3];
  double pn,tn,ph[2],t[2],cphi,sphi,cteta,steta;

   if(nphi >= 360.) nphi -=360.;
   if(nphi <= -360.) nphi +=360.;

   if(ntheta >= 360.) ntheta -=360.;
   if(ntheta <= -360.)ntheta +=360.;

  //��������� � �������
  pn =nphi*atan(1.)/45.; tn=(180-ntheta)*atan(1.)/45.;

  ph[0]=0.;
  ph[1]=pn;
  t[0]=tn;
  t[1]=0.;
 for(int k=0; k < 2; k++)
 {
  cphi=cos(ph[k]); sphi=sin(ph[k]);
  cteta=cos(t[k]); steta=sin(t[k]);
  //������ ����������� ����� ����� ������ Z
  vZ[0]=-cphi*steta;
  vZ[1]=-sphi*steta;
  vZ[2]=-cteta;
  //������ ����� Y
  vY[0]=-cphi*cteta;
  vY[1]=-sphi*cteta;
  vY[2]=steta;
  //������ vZ x vY
  vX[0]=vZ[1]*vY[2]-vZ[2]*vY[1];
  vX[1]=vZ[2]*vY[0]-vZ[0]*vY[2];
  vX[2]=vZ[0]*vY[1]-vZ[1]*vY[0];
  //������� �������������� ������ ��������� �� ������� � ������� �����������
  v11=vX[0]; v21=vX[1]; v31=vX[2];
  v12=vY[0]; v22=vY[1]; v32=vY[2];
  v13=vZ[0]; v23=vZ[1]; v33=vZ[2];
  v43=0;//rho;
  for(int i=0; i < 13; i++)
  {
    TPoint3D ps=p[i];
    TPoint3D to=temp[i];
    to.x=(ps.x*v11+ps.y*v21+ps.z*v31);
    to.y=(-ps.x*v12-ps.y*v22-ps.z*v32);
    to.z=(ps.x*v13+ps.y*v23+ps.z*v33+v43);
    p[i]=to;
  }
 }
}
// ���������� ����� ��� ��������� �����
// � �������������� ������ �����
void TEllipse3D::ToBezier()
{
  const double M=0.55228474983;
  double dx=radius*(1.-M);
  double dy=radius*(1.-M);
  p[0].x=radius;
  p[0].y=0.;
  p[0].z=0.;
  p[1].x=radius;
  p[1].y=radius-dy;
  p[1].z=0.;
  p[2].x=radius-dx;
  p[2].y=radius;
  p[2].z=0.;
  p[3].x=0.;
  p[3].y=radius;
  p[3].z=0.;

  p[4].x=-radius+dx;
  p[4].y=radius;
  p[4].z=0.;
  p[5].x=-radius;
  p[5].y=radius-dy;
  p[5].z=0.;
  p[6].x=-radius;
  p[6].y=0.;
  p[6].z=0.;

  p[7].x=-radius;
  p[7].y=-radius+dy;
  p[7].z=0.;
  p[8].x=-radius+dx;
  p[8].y=-radius;
  p[8].z=0.;
  p[9].x=0.;
  p[9].y=-radius;
  p[9].z=0.;

  p[10].x=radius-dx;
  p[10].y=-radius;
  p[10].z=0.;
  p[11].x=radius;
  p[11].y=-radius+dy;
  p[11].z=0.;
  p[12].x=radius;
  p[12].y=0.;
  p[12].z=0.;
  //������������, ��� �� ����������� ������� ������� ��������� � �������� normal
  Rotate();
  //��������� �����
  for(int i=0; i < 13; i++)
  {
    p[i].x+=center.x;
    p[i].y+=center.y;
    p[i].z+=center.z;
  }
}
//------------------------------------------------------------------
//��������������� ������� ���������� ������������ ������� 3x3
// a[0] b[0] c[0]
// a[1] b[1] c[1]
// a[2] b[2] c[2]
static double DetMatrix(double *a, double *b, double *c)
{
  double d;
  d  = a[0]*(b[1]*c[2]-b[2]*c[1]);
  d += b[0]*(a[2]*c[1]-a[1]*c[2]);
  d += c[0]*(a[1]*b[2]-a[2]*b[1]);
  return d;
}
//���������� ������������� ��������� ���� ��������
// �� �����:
// p0[2] - ������ � ������ ����� ������� �������
// p1[2] - ������ � ������ ����� ������� �������
//�� ������:
// pp[2] - ����� ������������� ��������� �� ������ � ������ �������
// ����������:
//        ����������� ���������� ����� ���������
double GetSectionDistance(TPoint3D *p0,TPoint3D *p1,TPoint3D *pp)
{
   double d0[3],d1[3]; //������� ������� � ������� ��������
   double vt0[3],vt1[3]; //������������ �������� ��������
   double vv[3]; //��������� ������������
   double det; //������������ �������
   double f;//��������������� ����������
   double t0,t1; //��������� ������
   double pp0[3]; //������ �� ������ ������� � ������ ������� ��������
   double r0,r1,r; //����� �������� � ���������� ����� ����;
   //������ ������ ����� ��� �������� �� ���������
   pp[0].x=p0[0].x; pp[0].y=p0[0].y; pp[0].z=p0[0].z;
   pp[1].x=p1[0].x; pp[1].y=p1[0].y; pp[1].z=p1[0].z;
  //������� ��������
   d0[0]=p0[1].x-p0[0].x; d0[1]=p0[1].y-p0[0].y; d0[2]=p0[1].z-p0[0].z;
   d1[0]=p1[1].x-p1[0].x; d1[1]=p1[1].y-p1[0].y; d1[2]=p1[1].z-p1[0].z;
   //��������� ����� ��������
   r0=sqrt(d0[0]*d0[0]+d0[1]*d0[1]+d0[2]*d0[2]);
   r1=sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2]);
   //��������� ������������ ������� ��������
   if(r0)
   {
     vt0[0]=d0[0]/r0;
     vt0[1]=d0[1]/r0;
     vt0[2]=d0[2]/r0;
   }
   if(r1)
   {
     vt1[0]=d1[0]/r1;
     vt1[1]=d1[1]/r1;
     vt1[2]=d1[2]/r1;
   }
  if(r0==0. || r1==0.)
  {
   if(r0==0. && r1==0.)//���� ������ � ������ ������� ������� �����
   {
     //��������� ����� �� ��������
     pp[0].x=p0[0].x; pp[0].y=p0[0].y; pp[0].z=p0[0].z;
     pp[1].x=p1[0].x; pp[1].y=p1[0].y; pp[1].z=p1[0].z;
   }
   else
    if(r0==0.)//���� ������ ������� ������� �����
    {
     //������� ����� �� ������ �������
     pp[0].x=p0[0].x; pp[0].y=p0[0].y; pp[0].z=p0[0].z;
     //��������� � ��� ����� �� ������ �������
     t1=vt1[0]*(pp[0].x-p1[0].x)+vt1[1]*(pp[0].y-p1[0].y)+vt1[2]*(pp[0].z-p1[0].z);
     if(t1 < 0.)t1=0;
     if(t1 > r1)t1=r1;
     pp[1].x=p1[0].x+vt1[0]*t1; pp[1].y=p1[0].y+vt1[1]*t1; pp[1].z=p1[0].z+vt1[2]*t1;
    }
    else  //���� ������ ������� ������� �����
     {
      //������� ����� �� ������ �������
      pp[1].x=p1[0].x; pp[1].y=p1[0].y; pp[1].z=p1[0].z;
      //��������� � ��� ����� �� ������ �������
      t0=vt0[0]*(pp[1].x-p0[0].x)+vt0[1]*(pp[1].y-p0[0].y)+vt0[2]*(pp[1].z-p0[0].z);
      if(t0 < 0.)t0=0;
      if(t0 > r0)t0=r0;
      pp[0].x=p0[0].x+vt0[0]*t0; pp[0].y=p0[0].y+vt0[1]*t0; pp[0].z=p0[0].z+vt0[2]*t0;
     }
  }//r0==0 || r1==0
  else //������� �� ������� �����
  {
   //��������� ������������ �������� vt0 � vt1
   vv[0]= vt0[1]*vt1[2]-vt0[2]*vt1[1];
   vv[1]=-vt0[0]*vt1[2]+vt0[2]*vt1[0];
   vv[2]= vt0[0]*vt1[1]-vt0[1]*vt1[0];
   //����������� �� ������ ������
   f =vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2];
   if(fabs(f)<1.e-6) //������� �����������
   {
     //������� ����� �� ������ �������
     t0=0.;
     pp[0]=p0[0];
     //��������� � ��� ����� �� ������ �������
     t1=vt1[0]*(pp[0].x-p1[0].x)+vt1[1]*(pp[0].y-p1[0].y)+vt1[2]*(pp[0].z-p1[0].z);
     t1 /=r1;
     if(t1 < 0.)
     {
      t1=0.;
      //��������� � ��� ����� �� ������ �������
      t0=vt0[0]*(p1[0].x-p0[0].x)+vt0[1]*(p1[0].y-p0[0].y)+vt0[2]*(p1[0].z-p0[0].z);
      t0 /=r0;
      if(t0 < 0.)
          t0=0.;
      else
        if(t0 > 1.)
           t0=1.;
     }
     else
      if(t1 > 1.)
          t1=1.;
   }
   else //�� �����������
   {
    f=sqrt(f);
    vv[0] /=f; vv[1] /=f; vv[2] /=f;
    //������������ ������� ���������
    double a[3],b[3],c[3],d[3];
    a[0]=vt0[0]; a[1]=vt0[1]; a[2]=vt0[2];
    b[0]= -vt1[0]; b[1]=-vt1[1]; b[2]=-vt1[2];
    c[0]=vv[0]; c[1]=vv[1]; c[2]=vv[2];
    d[0]=p1[0].x-p0[0].x; d[1]=p1[0].y-p0[0].y; d[2]=p1[0].z-p0[0].z;
    //������ ������������ �������
    det=DetMatrix(a,b,c);
    if(det==0.)return 0.; //��� �� ������ ����
    //������ t0
    f=DetMatrix(d,b,c);
    t0=f/det/r0;
    //������ t1
    f=DetMatrix(a,d,c);
    t1=f/det/r1;
    //��� ���� ������� ���� �� ��������� �����
    //����� ��������� � ������ ��� � ����� ������-�� �� ��������
    if(t0 < 0. || t0 > 1. || t1 < 0. || t1 > 1.)
    {
      //������� �������� ��� ��������� ����������
      // � ����� ����������� ���������
      //1. ����� �������� ��������
      t0=0; t1=0;
      r=pow(p0[0].x-p1[0].x,2)+pow(p0[0].y-p1[0].y,2)+pow(p0[0].z-p1[0].z,2);
      double rmin=r; //������ �����������
      //2. ����� ������� ��������
      r=pow(p0[1].x-p1[1].x,2)+pow(p0[1].y-p1[1].y,2)+pow(p0[1].z-p1[1].z,2);
      if(r < rmin)
      {
        t0=1.; t1=1.;
      }
      //3. ����� ������� ������� � ������ �������
      r=pow(p0[0].x-p1[1].x,2)+pow(p0[0].y-p1[1].y,2)+pow(p0[0].z-p1[1].z,2);
      if(r < rmin)
      {
        t0=0.; t1=1.;
      }
      //4. ����� ������ ������� � ������� �������
      r=pow(p0[1].x-p1[0].x,2)+pow(p0[1].y-p1[0].y,2)+pow(p0[1].z-p1[0].z,2);
      if(r < rmin)
      {
        t0=1.; t1=0.;
      }
     //5. �������� ������������� �� ������ ������� ������� �� ������
      //��������� � ��� ����� �� ������ �������
      double t=vt1[0]*(p0[0].x-p1[0].x)+vt1[1]*(p0[0].y-p1[0].y)+vt1[2]*(p0[0].z-p1[0].z);
      t /=r1;
      if(t+1.e-6 >= 0. && t <= 1.+1.e-6)
      {
       r=pow(p0[0].x-p1[0].x-d1[0]*t,2)+pow(p0[0].y-p1[0].y-d1[1]*t,2)+
         pow(p0[0].z-p1[0].z-d1[2]*t,2);
       if(r < rmin)
       {
        rmin=r;
        t0=0.; t1=t;
       }
      }
      //6. �������� ������������� �� ����� ������� ������� �� ������
       t=vt1[0]*(p0[1].x-p1[0].x)+vt1[1]*(p0[1].y-p1[0].y)+vt1[2]*(p0[1].z-p1[0].z);
       t /=r1;
       if(t+1.e-6 >= 0. && t <= 1.+1.e-6)
       {
        r=pow(p0[0].x-p1[0].x-d1[0]*t,2)+pow(p0[0].y-p1[0].y-d1[1]*t,2)+
          pow(p0[0].z-p1[0].z-d1[2]*t,2);
        if(r < rmin)
        {
         rmin=r;
         t0=1.; t1=t;
        }
       }
       //7. �������� ������������� �� ������ ������� ������� �� ������
       //��������� � ��� ����� �� ������ �������
        t=vt0[0]*(p1[0].x-p0[0].x)+vt0[1]*(p1[0].y-p0[0].y)+vt0[2]*(p1[0].z-p0[0].z);
        t /=r0;
       if(t+1.e-6 >= 0. && t <= 1.+1.e-6)
       {
        r=pow(p0[0].x+d0[0]*t-p1[0].x,2)+pow(p0[0].y+d0[1]*t-p1[0].y,2)+
          pow(p0[0].z+d0[2]*t-p1[0].z,2);
        if(r < rmin)
        {
         rmin=r;
         t0=t; t1=0.;
        }
       }
       //8. �������� ������������� �� ����� ������� ������� �� ������
        //��������� � ��� ����� �� ������ �������
        t=vt0[0]*(p1[1].x-p0[0].x)+vt0[1]*(p1[1].y-p0[0].y)+vt0[2]*(p1[1].z-p0[0].z);
        t /=r0;
       if(t+1.e-6 >= 0. && t <= 1.+1.e-6)
       {
        r=pow(p0[0].x+d0[0]*t-p1[1].x,2)+pow(p0[0].y+d0[1]*t-p1[1].y,2)+
          pow(p0[0].z+d0[2]*t-p1[1].z,2);
        if(r < rmin)
        {
         rmin=r;
         t0=t; t1=1.;
        }
       }
    }//�� �������� ��������
   }//�� �����������
    //��������� ����� ��������
   if(t0==0.)
     pp[0]=p0[0];
   else
    if(fabs(t0-1.)<1.e-6)
     pp[0]=p0[1];
    else
    {
      pp[0].x=p0[0].x+d0[0]*t0; pp[0].y=p0[0].y+d0[1]*t0; pp[0].z=p0[0].z+d0[2]*t0;
    }
   if(t1==0.)
     pp[1]=p1[0];
   else
    if(fabs(t1-1.)<1.e-6)
     pp[1]=p1[1];
    else
    {
     pp[1].x=p1[0].x+d1[0]*t1; pp[1].y=p1[0].y+d1[1]*t1; pp[1].z=p1[0].z+d1[2]*t1;
    }
  }//������� �� ������� �����
   //��������� ���������� ����� ���������
   r  =pow(pp[1].x-pp[0].x,2);
   r +=pow(pp[1].y-pp[0].y,2);
   r +=pow(pp[1].z-pp[0].z,2);
   r=sqrt(r);
   return r;
}
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
int CroosLineCone(TPoint3D Ce1, TPoint3D Cp1, TPoint3D Ce2, TPoint3D Cp2,
                  TPoint3D P1, TPoint3D P2, TPoint3D &cross1, TPoint3D &cross2,
                  bool bStoreT)
{
   double xp,yp,zp ; //������ ����� ������
   double xpv,ypv,zpv ; //������ ����� ������ � ��������� �������
   double ax, ay,az; //������������ ������ ������
   double axv, ayv,azv; //������������ ������ ������ � ��������� �������
   double bx0,by0,bz0; //��������� ������, ������ �� ������ ������ ����������
                    // � ������ ������
   double bx,by,bz; //��������� ������, ������ �� ������ ������ ����������
                    // � ������ ������ � ���������� �������
   double d; //����� ������� (C1,C2)
   double R1, R2; //������� �����������
   double dR; // ������� �������� R2-R1
   double r; //��������������� ����������
   double td; // ���������� �� ������ ������ ���������� �� ����������,
              // � ������� ���� ����������� �����
   double t[2]; // ����, � ������� ����� ���������� ����������

   double torad=atan(1.)/45. ; //������� �� �������� � �������
   TPoint3D Ce10, Ce20;
   TPoint3D P10, P20;
   //�������� �������� ����� ������
   P10=P1; P20=P2;
   //�������� �������� ����� ������� �����������
   Ce10=Ce1;
   Ce20=Ce2;
   xp=P1.x; yp=P1.y; zp=P1.z;
   ax=P2.x-P1.x; ay=P2.y-P1.y; az=P2.z-P1.z;
   r=sqrt(ax*ax+ay*ay+az*az);
   if(r < 1.e-3)return 0; //����� ����� ����������
   ax /=r; ay /=r; az /=r;
   bx0=Ce2.x-Ce1.x; by0=Ce2.y-Ce1.y; bz0=Ce2.z-Ce1.z;
   d=sqrt(bx0*bx0+by0*by0+bz0*bz0);
   if(d < 1.e-3)return 0; //����� ����� ���������
   bx0 /=d; by0 /=d; bz0 /=d;
   //������ ������� ������� � �������� ������ ��������� ���,
   // ����� ��� Z ���� ���������� ����� ������� ����������� ������
   // ����� ����� ����������
   // ������� ����� phi � theta ��� ��������
   double phi=Cp1.y*torad;
   double theta=Cp1.z*torad;
   //�������� ������� ��������
   double v11,v21,v31;
   double v12,v22,v32;
   double v13,v23,v33;
   v11=-sin(phi);            v21= cos(phi);             v31=0.;
   v12=-cos(phi)*cos(theta); v22=-sin(phi)*cos(theta); v32= sin(theta);
   v13= cos(phi)*sin(theta); v23= sin(phi)*sin(theta); v33= cos(theta);
   double pn[2],tn[2];
   pn[0]=Cp1.y*torad; tn[0]=Cp1.z*torad;
   pn[1]=Cp2.y*torad; tn[1]=Cp2.z*torad;
   //��������������� ������������ ����������� ������� �����������
   double k1=-cos(pn[0])*cos(tn[0]);
   double k2=-sin(pn[0]);
   double k3=-sin(pn[0])*cos(tn[0]);
   double k4=cos(pn[0]);
   double k5=sin(tn[0]);
/*
   double T1=0;//k1*v11+k3*v21+k5*v31;
   double T2=1;//k2*v11+k4*v21;
   double T3=1;//k1*v12+k3*v22+k5*v32;
   double T4=0;//k2*v12+k4*v22;
   double T5=0;//k1*v13+k3*v23+k5*v33;
   double T6=0;//k2*v13+k4*v23;
*/   
   //������� ����� ����� ������� �����������
   double x0=Ce10.x;
   double y0=Ce10.y;
   double z0=Ce10.z;
   Ce1.x=0.;//(Ce10.x-x0)*v11+(Ce10.y-y0)*v21+(Ce10.z-z0)*v31;
   Ce1.y=0.;//(Ce10.x-x0)*v12+(Ce10.y-y0)*v22+(Ce10.z-z0)*v32;
   Ce1.z=0.;//(Ce10.x-x0)*v13+(Ce10.y-y0)*v23+(Ce10.z-z0)*v33;
   Ce2.x=(Ce20.x-x0)*v11+(Ce20.y-y0)*v21+(Ce20.z-z0)*v31;
   Ce2.y=(Ce20.x-x0)*v12+(Ce20.y-y0)*v22+(Ce20.z-z0)*v32;
   Ce2.z=(Ce20.x-x0)*v13+(Ce20.y-y0)*v23+(Ce20.z-z0)*v33;
   bx=Ce2.x-Ce1.x; by=Ce2.y-Ce1.y; bz=Ce2.z-Ce1.z;
   d=sqrt(bx*bx+by*by+bz*bz);
   if(d < 1.e-3)return 0; //����� ����� ���������
   bx /=d; by /=d; bz /=d;
   R1=Cp1.x; R2=Cp2.x;
   dR=R2-R1; //����� ���� ���� ��� ��������
   //������� �����
   P1.x=(P10.x-x0)*v11+(P10.y-y0)*v21+(P10.z-z0)*v31;
   P1.y=(P10.x-x0)*v12+(P10.y-y0)*v22+(P10.z-z0)*v32;
   P1.z=(P10.x-x0)*v13+(P10.y-y0)*v23+(P10.z-z0)*v33;
   P2.x=(P20.x-x0)*v11+(P20.y-y0)*v21+(P20.z-z0)*v31;
   P2.y=(P20.x-x0)*v12+(P20.y-y0)*v22+(P20.z-z0)*v32;
   P2.z=(P20.x-x0)*v13+(P20.y-y0)*v23+(P20.z-z0)*v33;
   xpv=P1.x; ypv=P1.y; zpv=P1.z;
   axv=P2.x-P1.x; ayv=P2.y-P1.y; azv=P2.z-P1.z;
   r=sqrt(axv*axv+ayv*ayv+azv*azv);
   if(r < 1.e-3)return 0; //����� ����� ����������
   axv /=r; ayv /=r; azv /=r;
   //��� �������� ��� �����������
   if((fabs(dR) < 1.e-4 || (xpv==0 && ypv==0)) && fabs(axv) < 1.e-4 && fabs(ayv) < 1.e-4)
   {
    //��������� �� ������ ����������
     int idPoint=0;
     TPoint3D to[2];
        double tp=(Ce1.z-zpv)/azv;
        double x=xpv+axv*tp;
        double y=ypv+ayv*tp;
//        double z=zpv+azv*tp;
        if(x*x+y*y <= R1*R1)
        {
          to[0].x=xp+ax*tp;
          to[0].y=yp+ay*tp;
          to[0].z=zp+az*tp;
          idPoint |=1;
        }
        tp=(Ce2.z-zpv)/azv;
        x=xpv+axv*tp;
        y=ypv+ayv*tp;
//        z=zpv+azv*tp;
        if(x*x+y*y <= R2*R2)
        {
          to[1].x=xp+ax*tp;
          to[1].y=yp+ay*tp;
          to[1].z=zp+az*tp;
          idPoint |=2;
        }
        if(idPoint == 0)return 0;
        if(idPoint == 1)
        {
         cross1=to[0];
         return 1;
        }
        if(idPoint == 2)
        {
         cross1=to[1];
         return 1;
        }
        cross1=to[0];
        cross2=to[1];
        return 2;
   }
   //��������� ������ ��� azv �� ������ ����
   if(fabs(azv) < 1.e-4)
   {
     td=zpv/bz;
     if(td > d || td < 0.)return 0;
     double A=axv*(dR/d*td+R1);
     double B=-ayv*(dR/d*td+R1);
     double C=xpv*ayv-ypv*axv;
     double AB=sqrt(A*A+B*B);
     if(AB < 1.e-4)return 0; //��������� �� ������ ����������
     A /=AB; B /=AB;
     double D=-C/AB;
     if(fabs(D) - 1. > 0. )return 0;
     double alpha=0;
     if(A || B)alpha=atan2(B,A);
     t[0]=alpha+acos(D);
     t[1]=alpha-acos(D);
     TPoint3D to[2];
       int idPoint=0;
       for(int i=0; i < 2; i++)
       {
           double ct=cos(t[i]);
           double st=sin(t[i]);
           idPoint |=(1 << i);
           double radius=dR/d*td+R1;
           to[i].x=radius*ct*k1+radius*st*k2+Ce10.x+bx0*td;
           to[i].y=radius*ct*k3+radius*st*k4+Ce10.y+by0*td;
           to[i].z=radius*ct*k5+Ce10.z+bz0*td;
       }
       cross1=to[0];
       cross2=to[1];
       return 2;
   }
   //��������� ������ ��� azv �� ������ ����
   //������� ���� td=[A*cos(t)+B*sin(t)+C] / [D*cos(t)+E*sin(t)+F]
   double A=dR/d*(xpv-axv/azv*zpv)-R1*bz*axv/azv;
   double B=-dR/d*(ypv-ayv/azv*zpv)+R1*bz*ayv/azv;
   double C=-bz*ayv/azv*(xpv-axv/azv*zpv)+bz*axv/azv*(ypv-ayv/azv*zpv);
   double AB=sqrt(A*A+B*B);
   if(AB < 1.e-4)return 0;
   A /=AB; B /=AB;
   double D=-C/AB;
   if(fabs(D) - 1. > 0. )return 0;
   double alpha=0;
   if(A || B)alpha=atan2(B,A);
   t[0]=alpha+acos(D);
   t[1]=alpha-acos(D);
       TPoint3D to[2];
       int idPoint=0;
       for(int i=0; i < 2; i++)
       {
           double ct=cos(t[i]);
           double st=sin(t[i]);
           if(fabs(ayv) >= 1.e-4)
           {
            td=(ypv-ayv/azv*zpv-R1*ct)/(dR/d*ct-bz*ayv/azv);
           }
           else
           if(fabs(axv) >= 1.e-4)
           {
            td=(xpv-axv/azv*zpv-R1*st)/(dR/d*st-bz*axv/azv);
           }
           else
           {
             if(fabs(dR) > 1.e-4)
               td=(sqrt(xpv*xpv+ypv*ypv)-R1)*d/dR;
             else
              continue;
           }
           if(td > d || td < 0.)continue;
           idPoint |=(1 << i);
           double radius=dR/d*td+R1;
           to[i].x=radius*ct*k1+radius*st*k2+Ce10.x+bx0*td;
           to[i].y=radius*ct*k3+radius*st*k4+Ce10.y+by0*td;
           to[i].z=radius*ct*k5+Ce10.z+bz0*td;

           if(fabs(ayv) < 1.e-4 && fabs(axv) < 1.e-4)
           {
             if(fabs(to[i].x - xp) >= 1.e-4 && fabs(to[i].y - yp) >= 1.e-4)
             {
               idPoint ^=(1 << i);
             }

           }
       }
       if(idPoint == 0)
       {
        //�������� ����������� � ����������� ������
        double tp=(Ce1.z-zpv)/azv;
        double x=xpv+axv*tp;
        double y=ypv+ayv*tp;
//        double z=zpv+azv*tp;
        if(x*x+y*y <= R1*R1)
        {
          to[0].x=xp+ax*tp;
          to[0].y=yp+ay*tp;
          to[0].z=zp+az*tp;
          idPoint |=1;
        }
        tp=(Ce2.z-zpv)/azv;
        x=xpv+axv*tp;
        y=ypv+ayv*tp;
//        z=zpv+azv*tp;
        if(x*x+y*y <= R2*R2)
        {
          to[1].x=xp+ax*tp;
          to[1].y=yp+ay*tp;
          to[1].z=zp+az*tp;
          idPoint |=2;
        }
        if(idPoint == 0)return 0;
       }
       if(idPoint == 1)
       {
         cross1=to[0];
        //�������� ����������� � ����������� ������
        double tp=(Ce1.z-zpv)/azv;
        double x=xpv+axv*tp;
        double y=ypv+ayv*tp;
//        double z=zpv+azv*tp;
        if(x*x+y*y <= R1*R1)
        {
          to[1].x=xp+ax*tp;
          to[1].y=yp+ay*tp;
          to[1].z=zp+az*tp;
          idPoint |=2;
        }
        tp=(Ce2.z-zpv)/azv;
        x=xpv+axv*tp;
        y=ypv+ayv*tp;
//        z=zpv+azv*tp;
        if(x*x+y*y <= R2*R2)
        {
          to[1].x=xp+ax*tp;
          to[1].y=yp+ay*tp;
          to[1].z=zp+az*tp;
          idPoint |=2;
        }
        if(idPoint == 1)return 1;
       }
       if(idPoint == 2)
       {
         cross1=to[1];
        //�������� ����������� � ����������� ������
        double tp=(Ce1.z-zpv)/azv;
        double x=xpv+axv*tp;
        double y=ypv+ayv*tp;
//        double z=zpv+azv*tp;
        if(x*x+y*y <= R1*R1)
        {
          to[0].x=xp+ax*tp;
          to[0].y=yp+ay*tp;
          to[0].z=zp+az*tp;
          idPoint |=1;
        }
        tp=(Ce2.z-zpv)/azv;
        x=xpv+axv*tp;
        y=ypv+ayv*tp;
//        z=zpv+azv*tp;
        if(x*x+y*y <= R2*R2)
        {
          to[0].x=xp+ax*tp;
          to[0].y=yp+ay*tp;
          to[0].z=zp+az*tp;
          idPoint |=1;
        }
        if(idPoint == 2)return 1;
       }
       cross1=to[0];
       cross2=to[1];
       return 2;
}
