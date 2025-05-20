//+------------------------------------------------------------------+
//|                                               TestInterfaces.mqh |
//|            Copyright 2003-2022 Sergey Bochkanov (ALGLIB project) |
//|                             Copyright 2012-2023, MetaQuotes Ltd. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
//| Implementation of ALGLIB library in MetaQuotes Language 5        |
//|                                                                  |
//| The features of the library include:                             |
//| - Linear algebra (direct algorithms, EVD, SVD)                   |
//| - Solving systems of linear and non-linear equations             |
//| - Interpolation                                                  |
//| - Optimization                                                   |
//| - FFT (Fast Fourier Transform)                                   |
//| - Numerical integration                                          |
//| - Linear and nonlinear least-squares fitting                     |
//| - Ordinary differential equations                                |
//| - Computation of special functions                               |
//| - Descriptive statistics and hypothesis testing                  |
//| - Data analysis - classification, regression                     |
//| - Implementing linear algebra algorithms, interpolation, etc.    |
//|   in high-precision arithmetic (using MPFR)                      |
//|                                                                  |
//| This file is free software; you can redistribute it and/or       |
//| modify it under the terms of the GNU General Public License as   |
//| published by the Free Software Foundation (www.fsf.org); either  |
//| version 2 of the License, or (at your option) any later version. |
//|                                                                  |
//| This program is distributed in the hope that it will be useful,  |
//| but WITHOUT ANY WARRANTY; without even the implied warranty of   |
//| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     |
//| GNU General Public License for more details.                     |
//+------------------------------------------------------------------+
#include <Math\Alglib\alglib.mqh>

//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Func                            |
//+------------------------------------------------------------------+
class CNDimensional_Func1 : public CNDimensional_Func
  {
public:
   //--- constructor, destructor
                     CNDimensional_Func1(void) {}
                    ~CNDimensional_Func1(void) {}

   virtual void      Func(double &x[],double &func,CObject &obj);
   virtual void      Func(CRowDouble &x,double &func,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=100*(x0+3)^4 + (x1-3)^4        |
//+------------------------------------------------------------------+
void CNDimensional_Func1::Func(double &x[],double &func,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
  }
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=100*(x0+3)^4 + (x1-3)^4        |
//+------------------------------------------------------------------+
void CNDimensional_Func1::Func(CRowDouble &x,double &func,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Func                            |
//+------------------------------------------------------------------+
class CNDimensional_Func2 : public CNDimensional_Func
  {
public:
   //--- constructor, destructor
                     CNDimensional_Func2(void) {}
                    ~CNDimensional_Func2(void) {}

   virtual void      Func(double &x[],double &func,CObject &obj);
   virtual void      Func(CRowDouble &x,double &func,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=(x0^2+1)^2 + (x1-1)^2          |
//+------------------------------------------------------------------+
void CNDimensional_Func2::Func(double &x[],double &func,CObject &obj)
  {
   func=MathPow(x[0]*x[0]+1,2)+MathPow(x[1]-1,2);
  }
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=(x0^2+1)^2 + (x1-1)^2          |
//+------------------------------------------------------------------+
void CNDimensional_Func2::Func(CRowDouble &x,double &func,CObject &obj)
  {
   func=MathPow(x[0]*x[0]+1,2)+MathPow(x[1]-1,2);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Func                            |
//+------------------------------------------------------------------+
class CNDimensional_Bad_Func : public CNDimensional_Func
  {
public:
   //--- constructor, destructor
                     CNDimensional_Bad_Func(void) {}
                    ~CNDimensional_Bad_Func(void) {}

   virtual void      Func(double &x[],double &func,CObject &obj);
   virtual void      Func(CRowDouble &x,double &func,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function, i.e. function with      |
//| incorrectly calculated derivatives                               |
//+------------------------------------------------------------------+
void CNDimensional_Bad_Func::Func(double &x[],double &func,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
  }
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function, i.e. function with      |
//| incorrectly calculated derivatives                               |
//+------------------------------------------------------------------+
void CNDimensional_Bad_Func::Func(CRowDouble &x,double &func,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Grad                            |
//+------------------------------------------------------------------+
class CNDimensional_Grad1 : public CNDimensional_Grad
  {
public:
   //--- constructor, destructor
                     CNDimensional_Grad1(void) {}
                    ~CNDimensional_Grad1(void) {}

   virtual void      Grad(double &x[],double &func,double &grad[],CObject &obj);
   virtual void      Grad(CRowDouble &x,double &func,CRowDouble &grad,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=100*(x0+3)^4 + (x1-3)^4 and its|
//| derivatives df/d0 and df/dx1                                     |
//+------------------------------------------------------------------+
void CNDimensional_Grad1::Grad(double &x[],double &func,
                               double &grad[],CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
   grad[0]=400*MathPow(x[0]+3,3);
   grad[1]=4*MathPow(x[1]-3,3);
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_Grad1::Grad(CRowDouble &x,double &func,
                               CRowDouble &grad,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
   grad.Set(0,400*MathPow(x[0]+3,3));
   grad.Set(1,4*MathPow(x[1]-3,3));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Grad                            |
//+------------------------------------------------------------------+
class CNDimensional_Grad2 : public CNDimensional_Grad
  {
public:
   //--- constructor, destructor
                     CNDimensional_Grad2(void) {}
                    ~CNDimensional_Grad2(void) {}

   virtual void      Grad(double &x[],double &func,double &grad[],CObject &obj);
   virtual void      Grad(CRowDouble &x,double &func,CRowDouble &grad,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=(x0^2+1)^2 + (x1-1)^2 and its  |
//| derivatives df/d0 and df/dx1                                     |
//+------------------------------------------------------------------+
void CNDimensional_Grad2::Grad(double &x[],double &func,
                               double &grad[],CObject &obj)
  {
   func=MathPow(x[0]*x[0]+1,2)+MathPow(x[1]-1,2);
   grad[0]=4*(x[0]*x[0]+1)*x[0];
   grad[1]=2*(x[1]-1);
  }
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=(x0^2+1)^2 + (x1-1)^2 and its  |
//| derivatives df/d0 and df/dx1                                     |
//+------------------------------------------------------------------+
void CNDimensional_Grad2::Grad(CRowDouble &x,double &func,
                               CRowDouble &grad,CObject &obj)
  {
   func=MathPow(x[0]*x[0]+1,2)+MathPow(x[1]-1,2);
   grad.Set(0,4*(x[0]*x[0]+1)*x[0]);
   grad.Set(1,2*(x[1]-1));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Grad                            |
//+------------------------------------------------------------------+
class CNDimensional_Bad_Grad : public CNDimensional_Grad
  {
public:
   //--- constructor, destructor
                     CNDimensional_Bad_Grad(void) {}
                    ~CNDimensional_Bad_Grad(void) {}

   virtual void      Grad(double &x[],double &func,double &grad[],CObject &obj);
   virtual void      Grad(CRowDouble &x,double &func,CRowDouble &grad,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function, i.e. function with      |
//| incorrectly calculated derivatives                               |
//+------------------------------------------------------------------+
void CNDimensional_Bad_Grad::Grad(double &x[],double &func,double &grad[],CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
   grad[0]=40*MathPow(x[0]+3,3);
   grad[1]=40*MathPow(x[1]-3,3);
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_Bad_Grad::Grad(CRowDouble &x,double &func,CRowDouble &grad,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
   grad.Set(0,40*MathPow(x[0]+3,3));
   grad.Set(1,40*MathPow(x[1]-3,3));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Grad                            |
//+------------------------------------------------------------------+
class CNDimensional_S1_Grad : public CNDimensional_Grad
  {
public:
   //--- constructor, destructor
                     CNDimensional_S1_Grad(void) {}
                    ~CNDimensional_S1_Grad(void) {}

   virtual void      Grad(double &x[],double &func,double &grad[],CObject &obj);
   virtual void      Grad(CRowDouble &x,double &func,CRowDouble &grad,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(x)=(1+x)^(-0.2)+(1-x)^(-0.3)+1000*x   |
//| and its gradient. function is trimmed when we calculate it near  |
//| the singular points or outside of the [-1,+1]. Note that we do   |
//| NOT calculate gradient in this case.                             |
//+------------------------------------------------------------------+
void CNDimensional_S1_Grad::Grad(double &x[],double &func,double &grad[],CObject &obj)
  {
   if((x[0]<=-0.999999999999) || (x[0]>=+0.999999999999))
     {
      func=1.0E+300;
      return;
     }
   func=MathPow(1+x[0],-0.2)+MathPow(1-x[0],-0.3)+1000*x[0];
   grad[0]=-0.2*MathPow(1+x[0],-1.2)+0.3*MathPow(1-x[0],-1.3)+1000;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_S1_Grad::Grad(CRowDouble &x,double &func,CRowDouble &grad,CObject &obj)
  {
   if((x[0]<=-0.999999999999) || (x[0]>=+0.999999999999))
     {
      func=1.0E+300;
      return;
     }
   func=MathPow(1+x[0],-0.2)+MathPow(1-x[0],-0.3)+1000*x[0];
   grad.Set(0,-0.2*MathPow(1+x[0],-1.2)+0.3*MathPow(1-x[0],-1.3)+1000);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Hess                            |
//+------------------------------------------------------------------+
class CNDimensional_Hess1 : public CNDimensional_Hess
  {
public:
   //--- constructor, destructor
                     CNDimensional_Hess1(void) {}
                    ~CNDimensional_Hess1(void) {}

   virtual void      Hess(double &x[],double &func,double &grad[],CMatrixDouble &hess,CObject &obj);
   virtual void      Hess(CRowDouble &x,double &func,CRowDouble &grad,CMatrixDouble &hess,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=100*(x0+3)^4 + (x1-3)^4        |
//| its derivatives df/d0 and df/dx1                                 |
//| and its Hessian.                                                 |
//+------------------------------------------------------------------+
void CNDimensional_Hess1::Hess(double &x[],double &func,double &grad[],
                               CMatrixDouble &hess,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
   grad[0]=400*MathPow(x[0]+3,3);
   grad[1]=4*MathPow(x[1]-3,3);
   hess.Set(0,0,1200*MathPow(x[0]+3,2));
   hess.Set(0,1,0);
   hess.Set(1,0,0);
   hess.Set(1,1,12*MathPow(x[1]-3,2));
  }
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=100*(x0+3)^4 + (x1-3)^4        |
//| its derivatives df/d0 and df/dx1                                 |
//| and its Hessian.                                                 |
//+------------------------------------------------------------------+
void CNDimensional_Hess1::Hess(CRowDouble &x,double &func,CRowDouble &grad,
                               CMatrixDouble &hess,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
   grad.Set(0,400*MathPow(x[0]+3,3));
   grad.Set(1,4*MathPow(x[1]-3,3));
   hess.Set(0,0,1200*MathPow(x[0]+3,2));
   hess.Set(0,1,0);
   hess.Set(1,0,0);
   hess.Set(1,1,12*MathPow(x[1]-3,2));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Hess                            |
//+------------------------------------------------------------------+
class CNDimensional_Hess2 : public CNDimensional_Hess
  {
public:
   //--- constructor, destructor
                     CNDimensional_Hess2(void) {}
                    ~CNDimensional_Hess2(void) {}

   virtual void      Hess(double &x[],double &func,double &grad[],CMatrixDouble &hess,CObject &obj);
   virtual void      Hess(CRowDouble &x,double &func,CRowDouble &grad,CMatrixDouble &hess,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=(x0^2+1)^2 + (x1-1)^2          |
//| its gradient and Hessian                                         |
//+------------------------------------------------------------------+
void CNDimensional_Hess2::Hess(double &x[],double &func,double &grad[],
                               CMatrixDouble &hess,CObject &obj)
  {
   func=MathPow(x[0]*x[0]+1,2)+MathPow(x[1]-1,2);
   grad[0]=4*(x[0]*x[0]+1)*x[0];
   grad[1]=2*(x[1]-1);
   hess.Set(0,0,12*x[0]*x[0]+4);
   hess.Set(0,1,0);
   hess.Set(1,0,0);
   hess.Set(1,1,2);
  }
//+------------------------------------------------------------------+
//| This callback calculates f(x0,x1)=(x0^2+1)^2 + (x1-1)^2          |
//| its gradient and Hessian                                         |
//+------------------------------------------------------------------+
void CNDimensional_Hess2::Hess(CRowDouble &x,double &func,CRowDouble &grad,
                               CMatrixDouble &hess,CObject &obj)
  {
   func=MathPow(x[0]*x[0]+1,2)+MathPow(x[1]-1,2);
   grad.Set(0,4*(x[0]*x[0]+1)*x[0]);
   grad.Set(1,2*(x[1]-1));
   hess.Set(0,0,12*x[0]*x[0]+4);
   hess.Set(0,1,0);
   hess.Set(1,0,0);
   hess.Set(1,1,2);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Hess                            |
//+------------------------------------------------------------------+
class CNDimensional_Bad_Hess : public CNDimensional_Hess
  {
public:
   //--- constructor, destructor
                     CNDimensional_Bad_Hess(void) {}
                    ~CNDimensional_Bad_Hess(void) {}

   virtual void      Hess(double &x[],double &func,double &grad[],CMatrixDouble &hess,CObject &obj);
   virtual void      Hess(CRowDouble &x,double &func,CRowDouble &grad,CMatrixDouble &hess,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function, i.e. function with      |
//| incorrectly calculated derivatives                               |
//+------------------------------------------------------------------+
void CNDimensional_Bad_Hess::Hess(double &x[],double &func,double &grad[],
                                  CMatrixDouble &hess,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
   grad[0]=40*MathPow(x[0]+3,3);
   grad[1]=40*MathPow(x[1]-3,3);
   hess.Set(0,0,120*MathPow(x[0]+3,2));
   hess.Set(0,1,1);
   hess.Set(1,0,1);
   hess.Set(1,1,120*MathPow(x[1]-3,2));
  }
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function, i.e. function with      |
//| incorrectly calculated derivatives                               |
//+------------------------------------------------------------------+
void CNDimensional_Bad_Hess::Hess(CRowDouble &x,double &func,CRowDouble &grad,
                                  CMatrixDouble &hess,CObject &obj)
  {
   func=100*MathPow(x[0]+3,4)+MathPow(x[1]-3,4);
   grad.Set(0,40*MathPow(x[0]+3,3));
   grad.Set(1,40*MathPow(x[1]-3,3));
   hess.Set(0,0,120*MathPow(x[0]+3,2));
   hess.Set(0,1,1);
   hess.Set(1,0,1);
   hess.Set(1,1,120*MathPow(x[1]-3,2));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_FVec                            |
//+------------------------------------------------------------------+
class CNDimensional_FVec1 : public CNDimensional_FVec
  {
public:
   //--- constructor, destructor
                     CNDimensional_FVec1(void) {}
                    ~CNDimensional_FVec1(void) {}

   virtual void      FVec(double &x[],double &fi[],CObject &obj);
   virtual void      FVec(CRowDouble &x,CRowDouble &fi,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f0(x0,x1)=100*(x0+3)^4,                                          |
//| f1(x0,x1)=(x1-3)^4                                               |
//+------------------------------------------------------------------+
void CNDimensional_FVec1::FVec(double &x[],double &fi[],CObject &obj)
  {
   fi[0]=10*MathPow(x[0]+3,2);
   fi[1]=MathPow(x[1]-3,2);
  }
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f0(x0,x1)=100*(x0+3)^4,                                          |
//| f1(x0,x1)=(x1-3)^4                                               |
//+------------------------------------------------------------------+
void CNDimensional_FVec1::FVec(CRowDouble &x,CRowDouble &fi,CObject &obj)
  {
   fi.Set(0,10*MathPow(x[0]+3,2));
   fi.Set(1,MathPow(x[1]-3,2));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_FVec                            |
//+------------------------------------------------------------------+
class CNDimensional_FVec2 : public CNDimensional_FVec
  {
public:
   //--- constructor, destructor
                     CNDimensional_FVec2(void) {}
                    ~CNDimensional_FVec2(void) {}

   virtual void      FVec(double &x[],double &fi[],CObject &obj);
   virtual void      FVec(CRowDouble &x,CRowDouble &fi,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f0(x0,x1)=100*(x0+3)^4,                                          |
//| f1(x0,x1)=(x1-3)^4                                               |
//+------------------------------------------------------------------+
void CNDimensional_FVec2::FVec(double &x[],double &fi[],CObject &obj)
  {
   fi[0]=x[0]*x[0]+1;
   fi[1]=x[1]-1;
  }
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f0(x0,x1)=100*(x0+3)^4,                                          |
//| f1(x0,x1)=(x1-3)^4                                               |
//+------------------------------------------------------------------+
void CNDimensional_FVec2::FVec(CRowDouble &x,CRowDouble &fi,CObject &obj)
  {
   fi.Set(0,x[0]*x[0]+1);
   fi.Set(1,x[1]-1);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_FVec                            |
//+------------------------------------------------------------------+
class CNDimensional_Bad_FVec : public CNDimensional_FVec
  {
public:
   //--- constructor, destructor
                     CNDimensional_Bad_FVec(void) {}
                    ~CNDimensional_Bad_FVec(void) {}

   virtual void      FVec(double &x[],double &fi[],CObject &obj);
   virtual void      FVec(CRowDouble &x,CRowDouble &fi,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function, i.e. function with      |
//| incorrectly calculated derivatives                               |
//+------------------------------------------------------------------+
void CNDimensional_Bad_FVec::FVec(double &x[],double &fi[],CObject &obj)
  {
   fi[0]=10*MathPow(x[0]+3,2);
   fi[1]=MathPow(x[1]-3,2);
  }
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function, i.e. function with      |
//| incorrectly calculated derivatives                               |
//+------------------------------------------------------------------+
void CNDimensional_Bad_FVec::FVec(CRowDouble &x,CRowDouble &fi,CObject &obj)
  {
   fi.Set(0,10*MathPow(x[0]+3,2));
   fi.Set(1,MathPow(x[1]-3,2));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_FVec                            |
//+------------------------------------------------------------------+
class CNDimensional_NSFunc1_FVec : public CNDimensional_FVec
  {
public:
   //--- constructor, destructor
                     CNDimensional_NSFunc1_FVec(void) {}
                    ~CNDimensional_NSFunc1_FVec(void) {}

   virtual void      FVec(CRowDouble &x,CRowDouble &fi,CObject &obj);
  };
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_NSFunc1_FVec::FVec(CRowDouble &x,CRowDouble &fi,CObject &obj)
  {
//--- this callback calculates
//---     f0(x0,x1) = 2*|x0|+|x1|
   fi.Set(0,2*MathAbs(x[0])+MathAbs(x[1]));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Jac                             |
//+------------------------------------------------------------------+
class CNDimensional_Jac1 : public CNDimensional_Jac
  {
public:
   //--- constructor, destructor
                     CNDimensional_Jac1(void) {}
                    ~CNDimensional_Jac1(void) {}

   virtual void      Jac(double &x[],double &fi[],CMatrixDouble &jac,CObject &obj);
   virtual void      Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f0(x0,x1)=100*(x0+3)^4,                                          |
//| f1(x0,x1)=(x1-3)^4                                               |
//| and Jacobian matrix J=[dfi/dxj]                                  |
//+------------------------------------------------------------------+
void CNDimensional_Jac1::Jac(double &x[],double &fi[],CMatrixDouble &jac,
                             CObject &obj)
  {
   fi[0]=10*MathPow(x[0]+3,2);
   fi[1]=MathPow(x[1]-3,2);
   jac.Set(0,0,20*(x[0]+3));
   jac.Set(0,1,0);
   jac.Set(1,0,0);
   jac.Set(1,1,2*(x[1]-3));
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_Jac1::Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,
                             CObject &obj)
  {
   fi.Set(0,10*MathPow(x[0]+3,2));
   fi.Set(1,MathPow(x[1]-3,2));
   jac.Set(0,0,20*(x[0]+3));
   jac.Set(0,1,0);
   jac.Set(1,0,0);
   jac.Set(1,1,2*(x[1]-3));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Jac                             |
//+------------------------------------------------------------------+
class CNDimensional_Jac2 : public CNDimensional_Jac
  {
public:
   //--- constructor, destructor
                     CNDimensional_Jac2(void) {}
                    ~CNDimensional_Jac2(void) {}

   virtual void      Jac(double &x[],double &fi[],CMatrixDouble &jac,CObject &obj);
   virtual void      Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f0(x0,x1)=x0^2+1                                                 |
//| f1(x0,x1)=x1-1                                                   |
//| and Jacobian matrix J=[dfi/dxj]                                  |
//+------------------------------------------------------------------+
void CNDimensional_Jac2::Jac(double &x[],double &fi[],CMatrixDouble &jac,
                             CObject &obj)
  {
   fi[0]=x[0]*x[0]+1;
   fi[1]=x[1]-1;
   jac.Set(0,0,2*x[0]);
   jac.Set(0,1,0);
   jac.Set(1,0,0);
   jac.Set(1,1,1);
  }
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f0(x0,x1)=x0^2+1                                                 |
//| f1(x0,x1)=x1-1                                                   |
//| and Jacobian matrix J=[dfi/dxj]                                  |
//+------------------------------------------------------------------+
void CNDimensional_Jac2::Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,
                             CObject &obj)
  {
   fi.Set(0,x[0]*x[0]+1);
   fi.Set(1,x[1]-1);
   jac.Set(0,0,2*x[0]);
   jac.Set(0,1,0);
   jac.Set(1,0,0);
   jac.Set(1,1,1);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_Jac                             |
//+------------------------------------------------------------------+
class CNDimensional_Bad_Jac : public CNDimensional_Jac
  {
public:
   //--- constructor, destructor
                     CNDimensional_Bad_Jac(void) {}
                    ~CNDimensional_Bad_Jac(void) {}

   virtual void      Jac(double &x[],double &fi[],CMatrixDouble &jac,CObject &obj);
   virtual void      Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function,                         |
//| i.e. function with incorrectly calculated derivatives            |
//+------------------------------------------------------------------+
void CNDimensional_Bad_Jac::Jac(double &x[],double &fi[],CMatrixDouble &jac,
                                CObject &obj)
  {
   fi[0]=10*MathPow(x[0]+3,2);
   fi[1]=MathPow(x[1]-3,2);
   jac.Set(0,0,20*(x[0]+3));
   jac.Set(0,1,0);
   jac.Set(1,0,1);
   jac.Set(1,1,20*(x[1]-3));
  }
//+------------------------------------------------------------------+
//| This callback calculates 'bad' function,                         |
//| i.e. function with incorrectly calculated derivatives            |
//+------------------------------------------------------------------+
void CNDimensional_Bad_Jac::Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,
                                CObject &obj)
  {
   fi.Set(0,10*MathPow(x[0]+3,2));
   fi.Set(1,MathPow(x[1]-3,2));
   jac.Set(0,0,20*(x[0]+3));
   jac.Set(0,1,0);
   jac.Set(1,0,1);
   jac.Set(1,1,20*(x[1]-3));
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class CNDimensional_NLCFunc1_Jac : public CNDimensional_Jac
  {
public:
   //--- constructor, destructor
                     CNDimensional_NLCFunc1_Jac(void) {}
                    ~CNDimensional_NLCFunc1_Jac(void) {}

   virtual void      Jac(double &x[],double &fi[],CMatrixDouble &jac,CObject &obj);
   virtual void      Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj);
  };
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_NLCFunc1_Jac::Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj)
  {
//--- this callback calculates
//---     f0(x0,x1) = -x0+x1
//---     f1(x0,x1) = x0^2+x1^2-1
//--- and Jacobian matrix J = [dfi/dxj]
   fi.Set(0,-x[0]+x[1]);
   fi.Set(1,x[0]*x[0]+x[1]*x[1]-1.0);
   jac.Set(0,0,-1.0);
   jac.Set(0,1,+1.0);
   jac.Row(1,x*2+0);
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_NLCFunc1_Jac::Jac(double &x[],double &fi[],CMatrixDouble &jac,CObject &obj)
  {
//--- this callback calculates
//---     f0(x0,x1) = -x0+x1
//---     f1(x0,x1) = x0^2+x1^2-1
//--- and Jacobian matrix J = [dfi/dxj]
   fi[0]=-x[0]+x[1];
   fi[1]=x[0]*x[0]+x[1]*x[1]-1.0;
   jac.Set(0,0,-1.0);
   jac.Set(0,1,+1.0);
   jac.Set(1,0,x[0]*2);
   jac.Set(1,1,x[1]*2);
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class CNDimensional_NLCFunc2_Jac : public CNDimensional_Jac
  {
public:
   //--- constructor, destructor
                     CNDimensional_NLCFunc2_Jac(void) {}
                    ~CNDimensional_NLCFunc2_Jac(void) {}

   virtual void      Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj);
  };
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_NLCFunc2_Jac::Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj)
  {
//--- this callback calculates
//---     f0(x0,x1,x2) = x0+x1
//---     f1(x0,x1,x2) = x2-exp(x0)
//---     f2(x0,x1,x2) = x0^2+x1^2-1
//--- and Jacobian matrix J = [dfi/dxj]
   fi.Set(0,x[0]+x[1]);
   fi.Set(1,x[2]-MathExp(x[0]));
   fi.Set(2,x[0]*x[0]+x[1]*x[1]-1.0);
   jac.Set(0,0,1.0);
   jac.Set(0,1,1.0);
   jac.Set(0,2,0.0);
   jac.Set(1,0,-MathExp(x[0]));
   jac.Set(1,1,0.0);
   jac.Set(1,2,1.0);
   jac.Set(2,0,2*x[0]);
   jac.Set(2,1,2*x[1]);
   jac.Set(2,2,0.0);
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class CNDimensional_NSFunc1_Jac : public CNDimensional_Jac
  {
public:
   //--- constructor, destructor
                     CNDimensional_NSFunc1_Jac(void) {}
                    ~CNDimensional_NSFunc1_Jac(void) {}

   virtual void      Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj);
  };
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_NSFunc1_Jac::Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj)
  {
//--- this callback calculates
//---     f0(x0,x1) = 2*|x0|+|x1|
//--- and Jacobian matrix J = [ df0/dx0, df0/dx1 ]
   fi.Set(0,2*MathAbs(x[0])+MathAbs(x[1]));
   jac.Set(0,0,2*MathSign(x[0]));
   jac.Set(0,1,MathSign(x[1]));
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class CNDimensional_NSFunc2_Jac : public CNDimensional_Jac
  {
public:
   //--- constructor, destructor
                     CNDimensional_NSFunc2_Jac(void) {}
                    ~CNDimensional_NSFunc2_Jac(void) {}

   virtual void      Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj);
  };
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNDimensional_NSFunc2_Jac::Jac(CRowDouble &x,CRowDouble &fi,CMatrixDouble &jac,CObject &obj)
  {
//--- this callback calculates function vector
//---     f0(x0,x1) = 2*|x0|+x1
//---     f1(x0,x1) = x0-1
//---     f2(x0,x1) = -x1-1
//--- and Jacobian matrix J
//---         [ df0/dx0   df0/dx1 ]
//---     J = [ df1/dx0   df1/dx1 ]
//---         [ df2/dx0   df2/dx1 ]
   fi.Set(0,2*MathAbs(x[0])+MathAbs(x[1]));
   jac.Set(0,0,2*MathSign(x[0]));
   jac.Set(0,1,MathSign(x[1]));
   fi.Set(1,x[0]-1);
   jac.Set(1,0,1);
   jac.Set(1,1,0);
   fi.Set(2,-x[1]-1);
   jac.Set(2,0,0);
   jac.Set(2,1,-1);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_PFunc                           |
//+------------------------------------------------------------------+
class CNDimensional_CX_1_Func : public CNDimensional_PFunc
  {
public:
   //--- constructor, destructor
                     CNDimensional_CX_1_Func(void) {}
                    ~CNDimensional_CX_1_Func(void) {}

   virtual void      PFunc(double &c[],double &x[],double &func,CObject &obj);
   virtual void      PFunc(CRowDouble &c,CRowDouble &x,double &func,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(c,x)=exp(-c0*sqr(x0)) where x is a    |
//| position on X-axis and c is adjustable parameter                 |
//+------------------------------------------------------------------+
void CNDimensional_CX_1_Func::PFunc(CRowDouble &c,CRowDouble &x,double &func,
                                    CObject &obj)
  {
   func=MathExp(-c[0]*MathPow(x[0],2));
  }
//+------------------------------------------------------------------+
//| This callback calculates f(c,x)=exp(-c0*sqr(x0)) where x is a    |
//| position on X-axis and c is adjustable parameter                 |
//+------------------------------------------------------------------+
void CNDimensional_CX_1_Func::PFunc(double &c[],double &x[],double &func,
                                    CObject &obj)
  {
   func=MathExp(-c[0]*MathPow(x[0],2));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_PFunc                           |
//+------------------------------------------------------------------+
class CNDimensional_Debt_Func : public CNDimensional_PFunc
  {
public:
   //--- constructor, destructor
                     CNDimensional_Debt_Func(void) {}
                    ~CNDimensional_Debt_Func(void) {}

   virtual void      PFunc(double &c[],double &x[],double &func,CObject &obj);
   virtual void      PFunc(CRowDouble &c,CRowDouble &x,double &func,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f(c,x)=c[0]*(1+c[1]*(pow(x[0]-1999,c[2])-1))                     |
//+------------------------------------------------------------------+
void CNDimensional_Debt_Func::PFunc(double &c[],double &x[],double &func,
                                    CObject &obj)
  {
   func=c[0]*(1+c[1]*(MathPow(x[0]-1999,c[2])-1));
  }
//+------------------------------------------------------------------+
//| This callback calculates                                         |
//| f(c,x)=c[0]*(1+c[1]*(pow(x[0]-1999,c[2])-1))                     |
//+------------------------------------------------------------------+
void CNDimensional_Debt_Func::PFunc(CRowDouble &c,CRowDouble &x,double &func,
                                    CObject &obj)
  {
   func=c[0]*(1+c[1]*(MathPow(x[0]-1999,c[2])-1));
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_PGrad                           |
//+------------------------------------------------------------------+
class CNDimensional_CX_1_Grad : public CNDimensional_PGrad
  {
public:
   //--- constructor, destructor
                     CNDimensional_CX_1_Grad(void) {}
                    ~CNDimensional_CX_1_Grad(void) {}

   virtual void      PGrad(double &c[],double &x[],double &func,double &grad[],CObject &obj);
   virtual void      PGrad(CRowDouble &c,CRowDouble &x,double &func,CRowDouble &grad,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(c,x)=exp(-c0*sqr(x0)) and gradient    |
//| G={df/dc[i]} where x is a position on X-axis and c is adjustable |
//| parameter.                                                       |
//| IMPORTANT: gradient is calculated with respect to C, not to X    |
//+------------------------------------------------------------------+
void CNDimensional_CX_1_Grad::PGrad(double &c[],double &x[],double &func,
                                    double &grad[],CObject &obj)
  {
   func=MathExp(-c[0]*MathPow(x[0],2));
   grad[0]=-MathPow(x[0],2)*func;
  }
//+------------------------------------------------------------------+
//| This callback calculates f(c,x)=exp(-c0*sqr(x0)) and gradient    |
//| G={df/dc[i]} where x is a position on X-axis and c is adjustable |
//| parameter.                                                       |
//| IMPORTANT: gradient is calculated with respect to C, not to X    |
//+------------------------------------------------------------------+
void CNDimensional_CX_1_Grad::PGrad(CRowDouble &c,CRowDouble &x,double &func,
                                    CRowDouble &grad,CObject &obj)
  {
   func=MathExp(-c[0]*MathPow(x[0],2));
   grad.Set(0,-MathPow(x[0],2)*func);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_PHess                           |
//+------------------------------------------------------------------+
class CNDimensional_CX_1_Hess : public CNDimensional_PHess
  {
public:
   //--- constructor, destructor
                     CNDimensional_CX_1_Hess(void) {}
                    ~CNDimensional_CX_1_Hess(void) {}

   virtual void      PHess(double &c[],double &x[],double &func,double &grad[],CMatrixDouble &hess,CObject &obj);
   virtual void      PHess(CRowDouble &c,CRowDouble &x,double &func,CRowDouble &grad,CMatrixDouble &hess,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(c,x)=exp(-c0*sqr(x0)), gradient       |
//| G={df/dc[i]} and Hessian H={d2f/(dc[i]*dc[j])} where x is a      |
//| position on X-axis and c is adjustable parameter.                |
//| IMPORTANT: gradient/Hessian are calculated with respect to C,    |
//| not to X                                                         |
//+------------------------------------------------------------------+
void CNDimensional_CX_1_Hess::PHess(double &c[],double &x[],double &func,
                                    double &grad[],CMatrixDouble &hess,
                                    CObject &obj)
  {
   func=MathExp(-c[0]*MathPow(x[0],2));
   grad[0]=-MathPow(x[0],2)*func;
   hess.Set(0,0,MathPow(x[0],4)*func);
  }
//+------------------------------------------------------------------+
//| This callback calculates f(c,x)=exp(-c0*sqr(x0)), gradient       |
//| G={df/dc[i]} and Hessian H={d2f/(dc[i]*dc[j])} where x is a      |
//| position on X-axis and c is adjustable parameter.                |
//| IMPORTANT: gradient/Hessian are calculated with respect to C,    |
//| not to X                                                         |
//+------------------------------------------------------------------+
void CNDimensional_CX_1_Hess::PHess(CRowDouble &c,CRowDouble &x,double &func,
                                    CRowDouble &grad,CMatrixDouble &hess,
                                    CObject &obj)
  {
   func=MathExp(-c[0]*MathPow(x[0],2));
   grad.Set(0,-MathPow(x[0],2)*func);
   hess.Set(0,0,MathPow(x[0],4)*func);
  }
//+------------------------------------------------------------------+
//| Derived class from CNDimensional_ODE_RP                          |
//+------------------------------------------------------------------+
class CNDimensional_ODE_Function_1_Dif : public CNDimensional_ODE_RP
  {
public:
   //--- constructor, destructor
                     CNDimensional_ODE_Function_1_Dif(void) {}
                    ~CNDimensional_ODE_Function_1_Dif(void) {}

   virtual void      ODE_RP(double &y[],double x,double &dy[],CObject &obj);
   virtual void      ODE_RP(CRowDouble &y,double x,CRowDouble &dy,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(y[],x)=-y[0]                          |
//+------------------------------------------------------------------+
void CNDimensional_ODE_Function_1_Dif::ODE_RP(double &y[],double x,double &dy[],CObject &obj)
  {
   dy[0]=-y[0];
  }
//+------------------------------------------------------------------+
//| This callback calculates f(y[],x)=-y[0]                          |
//+------------------------------------------------------------------+
void CNDimensional_ODE_Function_1_Dif::ODE_RP(CRowDouble &y,double x,CRowDouble &dy,CObject &obj)
  {
   dy.Set(0,-y[0]);
  }
//+------------------------------------------------------------------+
//| Derived class from CIntegrator1_Func                             |
//+------------------------------------------------------------------+
class CInt_Function_1_Func : public CIntegrator1_Func
  {
public:
   //--- constructor, destructor
                     CInt_Function_1_Func(void) {}
                    ~CInt_Function_1_Func(void) {}

   virtual void      Int_Func(double x,double xminusa,double bminusx,double &y,CObject &obj);
  };
//+------------------------------------------------------------------+
//| This callback calculates f(x)=exp(x)                             |
//+------------------------------------------------------------------+
void CInt_Function_1_Func::Int_Func(double x,double xminusa,double bminusx,
                                    double &y,CObject &obj)
  {
   y=MathExp(x);
  }
//+------------------------------------------------------------------+
//| A comparison of the two numbers                                  |
//+------------------------------------------------------------------+
bool Doc_Test_Int(int val,int test_val)
  {
   return(val==test_val);
  }
//+------------------------------------------------------------------+
//| A comparison of the two numbers                                  |
//+------------------------------------------------------------------+
bool Doc_Test_Bool(bool val,bool test_val)
  {
   return(val==test_val);
  }
//+------------------------------------------------------------------+
//| A comparison of two numbers with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Real(double val,double test_val,double _threshold)
  {
//--- calculation
   double s=_threshold>=0?1.0:MathAbs(test_val);
   double threshold=MathAbs(_threshold);
//--- return result
   return(MathAbs(val-test_val)/s<=threshold);
  }
//+------------------------------------------------------------------+
//| A comparison of two numbers with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Complex(complex &val,complex &test_val,double _threshold)
  {
//--- calculation
   double s=_threshold>=0?1.0:CMath::AbsComplex(test_val);
   double threshold=MathAbs(_threshold);
//--- return result
   return(CMath::AbsComplex(val-test_val)/s<=threshold);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors                                      |
//+------------------------------------------------------------------+
bool Doc_Test_Int_Vector(int &val[],int &test_val[])
  {
//--- check
   if(CAp::Len(val)!=CAp::Len(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Len(val); i++)
      if(val[i]!=test_val[i])
         return(false);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors                                      |
//+------------------------------------------------------------------+
bool Doc_Test_Int_Vector(CRowInt &val,CRowInt &test_val)
  {
//--- check
   if(CAp::Len(val)!=CAp::Len(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Len(val); i++)
      if(val[i]!=test_val[i])
         return(false);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors                                      |
//+------------------------------------------------------------------+
bool Doc_Test_Int_Vector(CRowInt &val,int &test_val[])
  {
//--- check
   if(CAp::Len(val)!=CAp::Len(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Len(val); i++)
      if(val[i]!=test_val[i])
         return(false);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Comparison of the two matrices                                   |
//+------------------------------------------------------------------+
bool Doc_Test_Int_Matrix(CMatrixInt &val,CMatrixInt &test_val)
  {
//--- check
   if(CAp::Rows(val)!=CAp::Rows(test_val))
      return(false);
//--- check
   if(CAp::Cols(val)!=CAp::Cols(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Rows(val); i++)
      for(int j=0; j<CAp::Cols(val); j++)
         if(val.Get(i,j)!=test_val.Get(i,j))
            return(false);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Real_Vector(double &val[],double &test_val[],double _threshold)
  {
//--- check
   if(CAp::Len(val)!=CAp::Len(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Len(val); i++)
     {
      double s=_threshold>=0?1.0:MathAbs(test_val[i]);
      double threshold=MathAbs(_threshold);
      //--- check
      if(MathAbs(val[i]-test_val[i])/s>threshold)
         return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Real_Vector(CRowDouble &val,CRowDouble &test_val,double _threshold)
  {
//--- check
   if(CAp::Len(val)!=CAp::Len(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Len(val); i++)
     {
      double s=_threshold>=0?1.0:MathAbs(test_val[i]);
      double threshold=MathAbs(_threshold);
      //--- check
      if(MathAbs(val[i]-test_val[i])/s>threshold)
         return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Real_Vector(CRowDouble &val,double &test_val[],double _threshold)
  {
//--- check
   if(CAp::Len(val)!=CAp::Len(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Len(val); i++)
     {
      double s=_threshold>=0?1.0:MathAbs(test_val[i]);
      double threshold=MathAbs(_threshold);
      //--- check
      if(MathAbs(val[i]-test_val[i])/s>threshold)
         return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Real_Vector(CRowDouble &val,vector<double> &test_val,double _threshold)
  {
//--- check
   if(CAp::Len(val)!=(int)test_val.Size())
      return(false);
//--- comparison
   for(int i=0; i<CAp::Len(val); i++)
     {
      double s=_threshold>=0?1.0:MathAbs(test_val[i]);
      double threshold=MathAbs(_threshold);
      //--- check
      if(MathAbs(val[i]-test_val[i])/s>threshold)
         return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Real_Matrix(CMatrixDouble &val,CMatrixDouble &test_val,double _threshold)
  {
//--- check
   if(CAp::Rows(val)!=CAp::Rows(test_val))
      return(false);
//--- check
   if(CAp::Cols(val)!=CAp::Cols(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Rows(val); i++)
      for(int j=0; j<CAp::Cols(val); j++)
        {
         double s=_threshold>=0?1.0:MathAbs(test_val.Get(i,j));
         double threshold=MathAbs(_threshold);
         //--- check
         if(MathAbs(val.Get(i,j)-test_val.Get(i,j))/s>threshold)
            return(false);
        }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Real_Matrix(CMatrixDouble &val,matrix<double> &test_val,double _threshold)
  {
//--- check
   if(CAp::Rows(val)!=test_val.Rows())
      return(false);
//--- check
   if(CAp::Cols(val)!=test_val.Cols())
      return(false);
//--- comparison
   for(int i=0; i<CAp::Rows(val); i++)
      for(int j=0; j<CAp::Cols(val); j++)
        {
         double s=_threshold>=0?1.0:MathAbs(test_val[i,j]);
         double threshold=MathAbs(_threshold);
         //--- check
         if(MathAbs(val.Get(i,j)-test_val[i,j])/s>threshold)
            return(false);
        }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two vectors with an accuracy                     |
//+------------------------------------------------------------------+
bool Doc_Test_Complex_Vector(complex    &val[],complex    &test_val[],double _threshold)
  {
//--- check
   if(CAp::Len(val)!=CAp::Len(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Len(val); i++)
     {
      double s=_threshold>=0?1.0:CMath::AbsComplex(test_val[i]);
      double threshold=MathAbs(_threshold);
      //--- check
      if(CMath::AbsComplex(val[i]-test_val[i])/s>threshold)
         return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| A comparison of two matrices with an accuracy                    |
//+------------------------------------------------------------------+
bool Doc_Test_Complex_Matrix(CMatrixComplex &val,CMatrixComplex &test_val,double _threshold)
  {
//--- check
   if(CAp::Rows(val)!=CAp::Rows(test_val))
      return(false);
//--- check
   if(CAp::Cols(val)!=CAp::Cols(test_val))
      return(false);
//--- comparison
   for(int i=0; i<CAp::Rows(val); i++)
      for(int j=0; j<CAp::Cols(val); j++)
        {
         double s=_threshold>=0?1.0:CMath::AbsComplex(test_val.Get(i,j));
         double threshold=MathAbs(_threshold);
         //--- check
         if(CMath::AbsComplex(val.Get(i,j)-test_val.Get(i,j))/s>threshold)
            return(false);
        }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Add to the vector random element                                 |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Adding_Element(int &x[])
  {
//--- size calculation
   int n=ArraySize(x);
//--- increasing the length of the vector
   ArrayResize(x,n+1);
//--- set value
   x[n]=MathRand();
  }
//+------------------------------------------------------------------+
//| Add to the vector random element                                 |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Adding_Element(double &x[])
  {
//--- size calculation
   int n=ArraySize(x);
//--- increasing the length of the vector
   ArrayResize(x,n+1);
//--- set value
   x[n]=CMath::RandomReal();
  }
//+------------------------------------------------------------------+
//| Add to the vector random element                                 |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Adding_Element(CRowDouble &x)
  {
//--- size calculation
   int n=x.Size();
//--- increasing the length of the vector
   x.Resize(n+1);
//--- set value
   x.Set(n,CMath::RandomReal());
  }
//+------------------------------------------------------------------+
//| Add to the vector random element                                 |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Adding_Element(complex    &x[])
  {
//--- size calculation
   int n=ArraySize(x);
//--- increasing the length of the vector
   ArrayResize(x,n+1);
//--- set value
   x[n].real=CMath::RandomReal();
   x[n].imag=CMath::RandomReal();
  }
//+------------------------------------------------------------------+
//| Removing the number of vector                                    |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Deleting_Element(int &x[])
  {
//--- size calculation
   int n=ArraySize(x);
//--- reduction length of the vector
   ArrayResize(x,n-1);
  }
//+------------------------------------------------------------------+
//| Removing the number of vector                                    |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Deleting_Element(double &x[])
  {
//--- size calculation
   int n=ArraySize(x);
//--- reduction length of the vector
   ResetLastError();
   CAp::Assert(ArrayResize(x,n-1)==n-1,StringFormat("%s: error of resize array: %6d",__FUNCTION__,GetLastError()));
  }
//+------------------------------------------------------------------+
//| Removing the number of vector                                    |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Deleting_Element(CRowDouble &x)
  {
//--- size calculation
   int n=x.Size();
//--- reduction length of the vector
   x.Resize(n-1);
  }
//+------------------------------------------------------------------+
//| Removing the number of vector                                    |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Deleting_Element(CRowInt &x)
  {
//--- size calculation
   int n=x.Size();
//--- reduction length of the vector
   x.Resize(n-1);
  }
//+------------------------------------------------------------------+
//| Removing the number of vector                                    |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Deleting_Element(complex    &x[])
  {
//--- size calculation
   int n=ArraySize(x);
//--- reduction length of the vector
   ArrayResize(x,n-1);
  }
//+------------------------------------------------------------------+
//| Add a row in the matrix                                          |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Adding_Row(CMatrixInt &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- increase the dimension
   x.Resize(n+1,m);
//--- set values
   for(int i=0; i<m; i++)
      x.Set(n,i,MathRand());
  }
//+------------------------------------------------------------------+
//| Add a row in the matrix                                          |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Adding_Row(CMatrixDouble &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- increase the dimension
   x.Resize(n+1,m);
//--- set values
   for(int i=0; i<m; i++)
      x.Set(n,i,CMath::RandomReal());
  }
//+------------------------------------------------------------------+
//| Add a row in the matrix                                          |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Adding_Row(CMatrixComplex &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- increase the dimension
   x.Resize(n+1,m);
//--- set values
   for(int i=0; i<m; i++)
     {
      x.SetRe(n,i,CMath::RandomReal());
      x.SetIm(n,i,CMath::RandomReal());
     }
  }
//+------------------------------------------------------------------+
//| Deleting a row from the matrix                                   |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Deleting_Row(CMatrixInt &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- reduction of dimension
   x.Resize(n-1,m);
  }
//+------------------------------------------------------------------+
//| Deleting a row from the matrix                                   |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Deleting_Row(CMatrixDouble &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- reduction of dimension
   x.Resize(n-1,m);
  }
//+------------------------------------------------------------------+
//| Deleting a row from the matrix                                   |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Deleting_Row(CMatrixComplex &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- reduction of dimension
   x.Resize(n-1,m);
  }
//+------------------------------------------------------------------+
//| Add a col in the matrix                                          |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Adding_Col(CMatrixInt &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- increase the dimension
   x.Resize(n,m+1);
//--- set values
   for(int i=0; i<m; i++)
      x.Set(i,m,MathRand());
  }
//+------------------------------------------------------------------+
//| Add a col in the matrix                                          |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Adding_Col(CMatrixDouble &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- increase the dimension
   x.Resize(n,m+1);
//--- set values
   for(int i=0; i<m; i++)
      x.Set(i,m,CMath::RandomReal());
  }
//+------------------------------------------------------------------+
//| Add a col in the matrix                                          |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Adding_Col(CMatrixComplex &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- increase the dimension
   x.Resize(n,m+1);
//--- set values
   for(int i=0; i<m; i++)
     {
      x.SetRe(i,m,CMath::RandomReal());
      x.SetIm(i,m,CMath::RandomReal());
     }
  }
//+------------------------------------------------------------------+
//| Deleting a col from the matrix                                   |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Deleting_Col(CMatrixInt &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- reduction of dimension
   x.Resize(n,m-1);
  }
//+------------------------------------------------------------------+
//| Deleting a col from the matrix                                   |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Deleting_Col(CMatrixDouble &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- reduction of dimension
   x.Resize(n,m-1);
  }
//+------------------------------------------------------------------+
//| Deleting a col from the matrix                                   |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Deleting_Col(CMatrixComplex &x)
  {
//--- cols and rows calculation
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- reduction of dimension
   x.Resize(n,m-1);
  }
//+------------------------------------------------------------------+
//| Set number to a random position                                  |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Value(int &x[],int &val)
  {
//--- size calculation
   int n=ArraySize(x);
//--- set value
   if(n!=0)
      x[CMath::RandomInteger(n)]=val;
  }
//+------------------------------------------------------------------+
//| Set number to a random position                                  |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Value(double &x[],double val)
  {
//--- size calculation
   int n=ArraySize(x);
//--- set value
   if(n!=0)
      x[CMath::RandomInteger(n)]=val;
  }
//+------------------------------------------------------------------+
//| Set number to a random position                                  |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Value(CRowDouble &x,double val)
  {
//--- size calculation
   int n=x.Size();
//--- set value
   if(n!=0)
      x.Set(CMath::RandomInteger(n),val);
  }
//+------------------------------------------------------------------+
//| Set number to a random position                                  |
//+------------------------------------------------------------------+
void Spoil_Vector_By_Value(complex    &x[],complex    &val)
  {
//--- size calculation
   int n=ArraySize(x);
//--- set value
   if(n!=0)
      x[CMath::RandomInteger(n)]=val;
  }
//+------------------------------------------------------------------+
//| Set number to a random position                                  |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Value(CMatrixInt &x,int val)
  {
//--- get cols and rows
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- set value
   if(n!=0 && m!=0)
      x.Set(CMath::RandomInteger(n),CMath::RandomInteger(m),val);
  }
//+------------------------------------------------------------------+
//| Set number to a random position                                  |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Value(CMatrixDouble &x,double val)
  {
//--- get cols and rows
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- set value
   if(n!=0 && m!=0)
      x.Set(CMath::RandomInteger(n),CMath::RandomInteger(m),val);
  }
//+------------------------------------------------------------------+
//| Set number to a random position                                  |
//+------------------------------------------------------------------+
void Spoil_Matrix_By_Value(CMatrixComplex &x,complex    &val)
  {
//--- get cols and rows
   int n=CAp::Rows(x);
   int m=CAp::Cols(x);
//--- set value
   if(n!=0 && m!=0)
      x.Set(CMath::RandomInteger(n),CMath::RandomInteger(m),val);
  }
//+------------------------------------------------------------------+
//| Function test and exception handling                             |
//+------------------------------------------------------------------+
bool Func_spoil_scenario(int _spoil_scenario,bool &_TestResult)
  {
   if(CAp::exception_happened==true)
     {
      //--- check
      if(_spoil_scenario==-1)
         _TestResult=false;
      //--- reset exception
      CAp::exception_happened=false;
      //--- return result
      return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Nearest neighbor search, KNN queries                             |
//+------------------------------------------------------------------+
void TEST_NNeighbor_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
   int           nx;
   int           ny;
   int           normtype;
   CKDTreeShell  kdt;
   double        x[];
   CMatrixDouble r;
   int           k;
   CMatrixDouble tempmatrix;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- allocation
      a.Resize(4,2);
      //--- initialization
      a.Set(0,0,0);
      a.Set(0,1,0);
      a.Set(1,0,0);
      a.Set(1,1,1);
      a.Set(2,0,1);
      a.Set(2,1,0);
      a.Set(3,0,1);
      a.Set(3,1,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(a,CInfOrNaN::NegativeInfinity());
      //--- change values
      nx=2;
      ny=0;
      normtype=2;
      //--- function call
      CAlglib::KDTreeBuild(a,nx,ny,normtype,kdt);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x,2);
      //--- initialization
      x[0]=-1;
      x[1]=0;
      //--- function call
      k=CAlglib::KDTreeQueryKNN(kdt,x,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(k,1);
      //--- function call
      CAlglib::KDTreeQueryResultsX(kdt,r);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(1,2);
      //--- initialization
      tempmatrix.Set(0,0,0);
      tempmatrix.Set(0,1,0);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(r,tempmatrix,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Subsequent queries; buffered functions must use previously       |
//| allocated storage (if large enough), so buffer may contain some  |
//| info from previous call                                          |
//+------------------------------------------------------------------+
void TEST_NNeighbor_T_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
   int           nx;
   int           ny;
   int           normtype;
   CKDTreeShell  kdt;
   double        x[];
   CMatrixDouble rx;
   int           k;
   CMatrixDouble tempmatrix;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- allocation
      a.Resize(4,2);
      //--- initialization
      a.Set(0,0,0);
      a.Set(0,1,0);
      a.Set(1,0,0);
      a.Set(1,1,1);
      a.Set(2,0,1);
      a.Set(2,1,0);
      a.Set(3,0,1);
      a.Set(3,1,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(a,CInfOrNaN::NegativeInfinity());
      //--- change values
      nx=2;
      ny=0;
      normtype=2;
      //--- allocation
      rx.Resize(0,0);
      //--- function call
      CAlglib::KDTreeBuild(a,nx,ny,normtype,kdt);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x,2);
      //--- initialization
      x[0]=2;
      x[1]=0;
      //--- function call
      k=CAlglib::KDTreeQueryKNN(kdt,x,2,true);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(k,2);
      //--- function call
      CAlglib::KDTreeQueryResultsX(kdt,rx);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(2,2);
      //--- initialization
      tempmatrix.Set(0,0,1);
      tempmatrix.Set(0,1,0);
      tempmatrix.Set(1,0,1);
      tempmatrix.Set(1,1,1);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(rx,tempmatrix,0.05);
      //--- allocation
      ArrayResize(x,2);
      //--- initialization
      x[0]=-2;
      x[1]=0;
      //--- function call
      k=CAlglib::KDTreeQueryKNN(kdt,x,1,true);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(k,1);
      //--- function call
      CAlglib::KDTreeQueryResultsX(kdt,rx);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(2,2);
      //--- initialization
      tempmatrix.Set(0,0,0);
      tempmatrix.Set(0,1,0);
      tempmatrix.Set(1,0,1);
      tempmatrix.Set(1,1,1);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(rx,tempmatrix,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Serialization of KD-trees                                        |
//+------------------------------------------------------------------+
void TEST_NNeighbor_D_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
   int           nx;
   int           ny;
   int           normtype;
   CKDTreeShell  kdt0;
   CKDTreeShell  kdt1;
   string        s;
   double        x[];
   CMatrixDouble r0;
   CMatrixDouble r1;
   CMatrixDouble tempmatrix;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- allocation
      a.Resize(4,2);
      //--- initialization
      a.Set(0,0,0);
      a.Set(0,1,0);
      a.Set(1,0,0);
      a.Set(1,1,1);
      a.Set(2,0,1);
      a.Set(2,1,0);
      a.Set(3,0,1);
      a.Set(3,1,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(a,CInfOrNaN::NegativeInfinity());
      //--- change values
      nx=2;
      ny=0;
      normtype=2;
      //--- allocation
      r0.Resize(0,0);
      r1.Resize(0,0);
      //--- Build tree and serialize it
      CAlglib::KDTreeBuild(a,nx,ny,normtype,kdt0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::KDTreeSerialize(kdt0,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::KDTreeUnserialize(s,kdt1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Compare results from KNN queries
      ArrayResize(x,2);
      //--- initialization
      x[0]=-1;
      x[1]=0;
      //--- function call
      CAlglib::KDTreeQueryKNN(kdt0,x,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::KDTreeQueryResultsX(kdt0,r0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::KDTreeQueryKNN(kdt1,x,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::KDTreeQueryResultsX(kdt1,r1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(1,2);
      //--- initialization
      tempmatrix.Set(0,0,0);
      tempmatrix.Set(0,1,0);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(r0,tempmatrix,0.05);
      //--- allocation
      tempmatrix.Resize(1,2);
      //--- initialization
      tempmatrix.Set(0,0,0);
      tempmatrix.Set(0,1,0);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(r1,tempmatrix,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Basic functionality (moments,adev,median,percentile)             |
//+------------------------------------------------------------------+
void TEST_BaseStat_D_Base(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double mean;
   double variance;
   double skewness;
   double kurtosis;
   double adev;
   double p;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(x,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x[i]=i*i;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- Here we demonstrate calculation of sample moments
      //--- (mean,variance,skewness,kurtosis)
      CAlglib::SampleMoments(x,mean,variance,skewness,kurtosis);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(mean,28.5,0.01);
      _TestResult=_TestResult && Doc_Test_Real(variance,801.1667,0.01);
      _TestResult=_TestResult && Doc_Test_Real(skewness,0.5751,0.01);
      _TestResult=_TestResult && Doc_Test_Real(kurtosis,-1.2666,0.01);
      //--- Average deviation
      CAlglib::SampleAdev(x,adev);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(adev,23.2,0.01);
      //--- Median and percentile
      CAlglib::SampleMedian(x,v);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,20.5,0.01);
      p=0.5;
      //--- check
      if(_spoil_scenario==3)
         p=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         p=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         p=CInfOrNaN::NegativeInfinity();
      //--- function call
      CAlglib::SamplePercentile(x,p,v);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,20.5,0.01);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Correlation (covariance) between two random variables            |
//+------------------------------------------------------------------+
void TEST_BaseStat_D_C2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
     {
      //--- We have two samples - x and y,and want to measure dependency between them
      ArrayResize(x,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x[i]=i*i;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,10);
      //--- initialization
      for(int i=0; i<10; i++)
         y[i]=i;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- Three dependency measures are calculated:
      //--- * covariation
      //--- * Pearson correlation
      //--- * Spearman rank correlation
      v=CAlglib::Cov2(x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,82.5,0.001);
      //--- function call
      v=CAlglib::PearsonCorr2(x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,0.9627,0.001);
      //--- function call
      v=CAlglib::SpearmanCorr2(x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,1.000,0.001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Correlation (covariance) between components of random vector     |
//+------------------------------------------------------------------+
void TEST_BaseStat_D_CM(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble x;
   CMatrixDouble c;
   CMatrixDouble tempmatrix;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- X is a sample matrix:
      //--- * I-th row corresponds to I-th observation
      //--- * J-th column corresponds to J-th variable
      x.Resize(5,3);
      //--- initialization
      x.Set(0,0,1);
      x.Set(0,1,0);
      x.Set(0,2,1);
      x.Set(1,0,1);
      x.Set(1,1,1);
      x.Set(1,2,0);
      x.Set(2,0,-1);
      x.Set(2,1,1);
      x.Set(2,2,0);
      x.Set(3,0,-2);
      x.Set(3,1,-1);
      x.Set(3,2,1);
      x.Set(4,0,-1);
      x.Set(4,1,0);
      x.Set(4,2,9);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- Three dependency measures are calculated:
      //--- * covariation
      //--- * Pearson correlation
      //--- * Spearman rank correlation
      //--- Result is stored into C,with C[i,j] equal to correlation
      //--- (covariance) between I-th and J-th variables of X.
      CAlglib::CovM(x,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(3,3);
      //--- initialization
      tempmatrix.Set(0,0,1.8);
      tempmatrix.Set(0,1,0.6);
      tempmatrix.Set(0,2,-1.4);
      tempmatrix.Set(1,0,0.6);
      tempmatrix.Set(1,1,0.7);
      tempmatrix.Set(1,2,-0.8);
      tempmatrix.Set(2,0,-1.4);
      tempmatrix.Set(2,1,-0.8);
      tempmatrix.Set(2,2,14.7);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(c,tempmatrix,0.01);
      //--- function call
      CAlglib::PearsonCorrM(x,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(3,3);
      //--- initialization
      tempmatrix.Set(0,0,1);
      tempmatrix.Set(0,1,0.535);
      tempmatrix.Set(0,2,-0.272);
      tempmatrix.Set(1,0,0.535);
      tempmatrix.Set(1,1,1);
      tempmatrix.Set(1,2,-0.249);
      tempmatrix.Set(2,0,-0.272);
      tempmatrix.Set(2,1,-0.249);
      tempmatrix.Set(2,2,1);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(c,tempmatrix,0.01);
      //--- function call
      CAlglib::SpearmanCorrM(x,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(3,3);
      //--- initialization
      tempmatrix.Set(0,0,1);
      tempmatrix.Set(0,1,0.556);
      tempmatrix.Set(0,2,-0.306);
      tempmatrix.Set(1,0,0.556);
      tempmatrix.Set(1,1,1);
      tempmatrix.Set(1,2,-0.75);
      tempmatrix.Set(2,0,-0.306);
      tempmatrix.Set(2,1,-0.75);
      tempmatrix.Set(2,2,1);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(c,tempmatrix,0.01);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Correlation (covariance) between two random vectors              |
//+------------------------------------------------------------------+
void TEST_BaseStat_D_CM2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble x;
   CMatrixDouble y;
   CMatrixDouble tempmatrix;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- X and Y are sample matrices:
      //--- * I-th row corresponds to I-th observation
      //--- * J-th column corresponds to J-th variable
      x.Resize(5,3);
      //--- initialization
      x.Set(0,0,1);
      x.Set(0,1,0);
      x.Set(0,2,1);
      x.Set(1,0,1);
      x.Set(1,1,1);
      x.Set(1,2,0);
      x.Set(2,0,-1);
      x.Set(2,1,1);
      x.Set(2,2,0);
      x.Set(3,0,-2);
      x.Set(3,1,-1);
      x.Set(3,2,1);
      x.Set(4,0,-1);
      x.Set(4,1,0);
      x.Set(4,2,9);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- allocation
      y.Resize(5,2);
      //--- initialization
      y.Set(0,0,2);
      y.Set(0,1,3);
      y.Set(1,0,2);
      y.Set(1,1,1);
      y.Set(2,0,-1);
      y.Set(2,1,6);
      y.Set(3,0,-9);
      y.Set(3,1,9);
      y.Set(4,0,7);
      y.Set(4,1,1);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Value(y,CInfOrNaN::NegativeInfinity());
      CMatrixDouble c;
      //--- Three dependency measures are calculated:
      //--- * covariation
      //--- * Pearson correlation
      //--- * Spearman rank correlation
      //--- Result is stored into C,with C[i,j] equal to correlation
      //--- (covariance) between I-th variable of X and J-th variable of Y.
      CAlglib::CovM2(x,y,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(3,2);
      //--- initialization
      tempmatrix.Set(0,0,4.1);
      tempmatrix.Set(0,1,-3.25);
      tempmatrix.Set(1,0,2.45);
      tempmatrix.Set(1,1,-1.5);
      tempmatrix.Set(2,0,13.45);
      tempmatrix.Set(2,1,-5.75);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(c,tempmatrix,0.01);
      //--- function call
      CAlglib::PearsonCorrM2(x,y,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(3,2);
      //--- initialization
      tempmatrix.Set(0,0,0.519);
      tempmatrix.Set(0,1,-0.699);
      tempmatrix.Set(1,0,0.497);
      tempmatrix.Set(1,1,-0.518);
      tempmatrix.Set(2,0,0.596);
      tempmatrix.Set(2,1,-0.433);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(c,tempmatrix,0.01);
      //--- function call
      CAlglib::SpearmanCorrM2(x,y,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(3,2);
      //--- initialization
      tempmatrix.Set(0,0,0.541);
      tempmatrix.Set(0,1,-0.649);
      tempmatrix.Set(1,0,0.216);
      tempmatrix.Set(1,1,-0.433);
      tempmatrix.Set(2,0,0.433);
      tempmatrix.Set(2,1,-0.135);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(c,tempmatrix,0.01);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Tests ability to detect errors in inputs                         |
//+------------------------------------------------------------------+
void TEST_BaseStat_T_Base(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double mean;
   double variance;
   double skewness;
   double kurtosis;
   double adev;
   double p;
   double v;
   double x1[];
   double x2[];
   double x3[];
   double x4[];
   double x5[];
   double x6[];
   double x7[];
   double x8[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<34; _spoil_scenario++)
     {
      //--- first,we test short form of functions
      ArrayResize(x1,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x1[i]=i*i;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x1,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x1,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x1,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::SampleMoments(x1,mean,variance,skewness,kurtosis);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x2,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x2[i]=i*i;
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(x2,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(x2,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(x2,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::SampleAdev(x2,adev);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x3,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x3[i]=i*i;
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(x3,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(x3,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(x3,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::SampleMedian(x3,v);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x4,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x4[i]=i*i;
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(x4,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(x4,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(x4,CInfOrNaN::NegativeInfinity());
      //--- change value
      p=0.5;
      //--- check
      if(_spoil_scenario==12)
         p=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==13)
         p=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==14)
         p=CInfOrNaN::NegativeInfinity();
      //--- function call
      CAlglib::SamplePercentile(x4,p,v);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- and then we test full form
      ArrayResize(x5,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x5[i]=i*i;
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(x5,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(x5,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(x5,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Deleting_Element(x5);
      //--- function call
      CAlglib::SampleMoments(x5,10,mean,variance,skewness,kurtosis);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x6,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x6[i]=i*i;
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Value(x6,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==20)
         Spoil_Vector_By_Value(x6,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==21)
         Spoil_Vector_By_Value(x6,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==22)
         Spoil_Vector_By_Deleting_Element(x6);
      //--- function call
      CAlglib::SampleAdev(x6,10,adev);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x7,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x7[i]=i*i;
      //--- check
      if(_spoil_scenario==23)
         Spoil_Vector_By_Value(x7,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==24)
         Spoil_Vector_By_Value(x7,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==25)
         Spoil_Vector_By_Value(x7,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==26)
         Spoil_Vector_By_Deleting_Element(x7);
      //--- function call
      CAlglib::SampleMedian(x7,10,v);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x8,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x8[i]=i*i;
      //--- check
      if(_spoil_scenario==27)
         Spoil_Vector_By_Value(x8,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==28)
         Spoil_Vector_By_Value(x8,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==29)
         Spoil_Vector_By_Value(x8,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==30)
         Spoil_Vector_By_Deleting_Element(x8);
      //--- change value
      p=0.5;
      //--- check
      if(_spoil_scenario==31)
         p=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==32)
         p=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==33)
         p=CInfOrNaN::NegativeInfinity();
      //--- function call
      CAlglib::SamplePercentile(x8,10,p,v);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Tests ability to detect errors in inputs                         |
//+------------------------------------------------------------------+
void TEST_BaseStat_T_CovCorr(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double        v;
   CMatrixDouble c;
   double        x1[];
   double        x2[];
   double        x3[];
   double        y1[];
   double        y2[];
   double        y3[];
   double        x1a[];
   double        x2a[];
   double        x3a[];
   double        y1a[];
   double        y2a[];
   double        y3a[];
   CMatrixDouble x4;
   CMatrixDouble x5;
   CMatrixDouble x6;
   CMatrixDouble x7;
   CMatrixDouble x8;
   CMatrixDouble x9;
   CMatrixDouble x10;
   CMatrixDouble x11;
   CMatrixDouble x12;
   CMatrixDouble x13;
   CMatrixDouble x14;
   CMatrixDouble x15;
   CMatrixDouble y10;
   CMatrixDouble y11;
   CMatrixDouble y12;
   CMatrixDouble y13;
   CMatrixDouble y14;
   CMatrixDouble y15;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<126; _spoil_scenario++)
     {
      //--- 2-sample short-form cov/corr are tested
      ArrayResize(x1,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x1[i]=i*i;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x1,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x1,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x1,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x1);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x1);
      //--- allocation
      ArrayResize(y1,10);
      //--- initialization
      for(int i=0; i<10; i++)
         y1[i]=i;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y1,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y1,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y1,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y1);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y1);
      //--- function call
      v=CAlglib::Cov2(x1,y1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x2,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x2[i]=i*i;
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(x2,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(x2,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(x2,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==13)
         Spoil_Vector_By_Adding_Element(x2);
      //--- check
      if(_spoil_scenario==14)
         Spoil_Vector_By_Deleting_Element(x2);
      //--- allocation
      ArrayResize(y2,10);
      //--- initialization
      for(int i=0; i<10; i++)
         y2[i]=i;
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(y2,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(y2,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(y2,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Adding_Element(y2);
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Deleting_Element(y2);
      //--- function call
      v=CAlglib::PearsonCorr2(x2,y2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x3,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x3[i]=i*i;
      //--- check
      if(_spoil_scenario==20)
         Spoil_Vector_By_Value(x3,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==21)
         Spoil_Vector_By_Value(x3,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==22)
         Spoil_Vector_By_Value(x3,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==23)
         Spoil_Vector_By_Adding_Element(x3);
      //--- check
      if(_spoil_scenario==24)
         Spoil_Vector_By_Deleting_Element(x3);
      //--- allocation
      ArrayResize(y3,10);
      //--- initialization
      for(int i=0; i<10; i++)
         y3[i]=i;
      //--- check
      if(_spoil_scenario==25)
         Spoil_Vector_By_Value(y3,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==26)
         Spoil_Vector_By_Value(y3,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==27)
         Spoil_Vector_By_Value(y3,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==28)
         Spoil_Vector_By_Adding_Element(y3);
      //--- check
      if(_spoil_scenario==29)
         Spoil_Vector_By_Deleting_Element(y3);
      //--- function call
      v=CAlglib::SpearmanCorr2(x3,y3);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- 2-sample full-form cov/corr are tested
      ArrayResize(x1a,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x1a[i]=i*i;
      //--- check
      if(_spoil_scenario==30)
         Spoil_Vector_By_Value(x1a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==31)
         Spoil_Vector_By_Value(x1a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==32)
         Spoil_Vector_By_Value(x1a,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==33)
         Spoil_Vector_By_Deleting_Element(x1a);
      //--- allocation
      ArrayResize(y1a,10);
      //--- initialization
      for(int i=0; i<10; i++)
         y1a[i]=i;
      //--- check
      if(_spoil_scenario==34)
         Spoil_Vector_By_Value(y1a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==35)
         Spoil_Vector_By_Value(y1a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==36)
         Spoil_Vector_By_Value(y1a,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==37)
         Spoil_Vector_By_Deleting_Element(y1a);
      //--- function call
      v=CAlglib::Cov2(x1a,y1a,10);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x2a,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x2a[i]=i*i;
      //--- check
      if(_spoil_scenario==38)
         Spoil_Vector_By_Value(x2a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==39)
         Spoil_Vector_By_Value(x2a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==40)
         Spoil_Vector_By_Value(x2a,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==41)
         Spoil_Vector_By_Deleting_Element(x2a);
      //--- allocation
      ArrayResize(y2a,10);
      //--- initialization
      for(int i=0; i<10; i++)
         y2a[i]=i;
      //--- check
      if(_spoil_scenario==42)
         Spoil_Vector_By_Value(y2a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==43)
         Spoil_Vector_By_Value(y2a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==44)
         Spoil_Vector_By_Value(y2a,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==45)
         Spoil_Vector_By_Deleting_Element(y2a);
      //--- function call
      v=CAlglib::PearsonCorr2(x2a,y2a,10);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(x3a,10);
      //--- initialization
      for(int i=0; i<10; i++)
         x3a[i]=i*i;
      //--- check
      if(_spoil_scenario==46)
         Spoil_Vector_By_Value(x3a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==47)
         Spoil_Vector_By_Value(x3a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==48)
         Spoil_Vector_By_Value(x3a,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==49)
         Spoil_Vector_By_Deleting_Element(x3a);
      //--- allocation
      ArrayResize(y3a,10);
      //--- initialization
      for(int i=0; i<10; i++)
         y3a[i]=i;
      //--- check
      if(_spoil_scenario==50)
         Spoil_Vector_By_Value(y3a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==51)
         Spoil_Vector_By_Value(y3a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==52)
         Spoil_Vector_By_Value(y3a,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==53)
         Spoil_Vector_By_Deleting_Element(y3a);
      //--- function call
      v=CAlglib::SpearmanCorr2(x3a,y3a,10);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- vector short-form cov/corr are tested.
      x4.Resize(5,3);
      //--- initialization
      x4.Set(0,0,1);
      x4.Set(0,1,0);
      x4.Set(0,2,1);
      x4.Set(1,0,1);
      x4.Set(1,1,1);
      x4.Set(1,2,0);
      x4.Set(2,0,-1);
      x4.Set(2,1,1);
      x4.Set(2,2,0);
      x4.Set(3,0,-2);
      x4.Set(3,1,-1);
      x4.Set(3,2,1);
      x4.Set(4,0,-1);
      x4.Set(4,1,0);
      x4.Set(4,2,9);
      //--- check
      if(_spoil_scenario==54)
         Spoil_Matrix_By_Value(x4,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==55)
         Spoil_Matrix_By_Value(x4,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==56)
         Spoil_Matrix_By_Value(x4,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::CovM(x4,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      x5.Resize(5,3);
      //--- initialization
      x5.Set(0,0,1);
      x5.Set(0,1,0);
      x5.Set(0,2,1);
      x5.Set(1,0,1);
      x5.Set(1,1,1);
      x5.Set(1,2,0);
      x5.Set(2,0,-1);
      x5.Set(2,1,1);
      x5.Set(2,2,0);
      x5.Set(3,0,-2);
      x5.Set(3,1,-1);
      x5.Set(3,2,1);
      x5.Set(4,0,-1);
      x5.Set(4,1,0);
      x5.Set(4,2,9);
      //--- check
      if(_spoil_scenario==57)
         Spoil_Matrix_By_Value(x5,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==58)
         Spoil_Matrix_By_Value(x5,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==59)
         Spoil_Matrix_By_Value(x5,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::PearsonCorrM(x5,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      x6.Resize(5,3);
      //--- initialization
      x6.Set(0,0,1);
      x6.Set(0,1,0);
      x6.Set(0,2,1);
      x6.Set(1,0,1);
      x6.Set(1,1,1);
      x6.Set(1,2,0);
      x6.Set(2,0,-1);
      x6.Set(2,1,1);
      x6.Set(2,2,0);
      x6.Set(3,0,-2);
      x6.Set(3,1,-1);
      x6.Set(3,2,1);
      x6.Set(4,0,-1);
      x6.Set(4,1,0);
      x6.Set(4,2,9);
      //--- check
      if(_spoil_scenario==60)
         Spoil_Matrix_By_Value(x6,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==61)
         Spoil_Matrix_By_Value(x6,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==62)
         Spoil_Matrix_By_Value(x6,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::SpearmanCorrM(x6,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- vector full-form cov/corr are tested.
      x7.Resize(5,3);
      //--- initialization
      x7.Set(0,0,1);
      x7.Set(0,1,0);
      x7.Set(0,2,1);
      x7.Set(1,0,1);
      x7.Set(1,1,1);
      x7.Set(1,2,0);
      x7.Set(2,0,-1);
      x7.Set(2,1,1);
      x7.Set(2,2,0);
      x7.Set(3,0,-2);
      x7.Set(3,1,-1);
      x7.Set(3,2,1);
      x7.Set(4,0,-1);
      x7.Set(4,1,0);
      x7.Set(4,2,9);
      //--- check
      if(_spoil_scenario==63)
         Spoil_Matrix_By_Value(x7,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==64)
         Spoil_Matrix_By_Value(x7,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==65)
         Spoil_Matrix_By_Value(x7,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==66)
         Spoil_Matrix_By_Deleting_Row(x7);
      //--- check
      if(_spoil_scenario==67)
         Spoil_Matrix_By_Deleting_Col(x7);
      //--- function call
      CAlglib::CovM(x7,5,3,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      x8.Resize(5,3);
      //--- initialization
      x8.Set(0,0,1);
      x8.Set(0,1,0);
      x8.Set(0,2,1);
      x8.Set(1,0,1);
      x8.Set(1,1,1);
      x8.Set(1,2,0);
      x8.Set(2,0,-1);
      x8.Set(2,1,1);
      x8.Set(2,2,0);
      x8.Set(3,0,-2);
      x8.Set(3,1,-1);
      x8.Set(3,2,1);
      x8.Set(4,0,-1);
      x8.Set(4,1,0);
      x8.Set(4,2,9);
      //--- check
      if(_spoil_scenario==68)
         Spoil_Matrix_By_Value(x8,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==69)
         Spoil_Matrix_By_Value(x8,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==70)
         Spoil_Matrix_By_Value(x8,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==71)
         Spoil_Matrix_By_Deleting_Row(x8);
      //--- check
      if(_spoil_scenario==72)
         Spoil_Matrix_By_Deleting_Col(x8);
      //--- function call
      CAlglib::PearsonCorrM(x8,5,3,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      x9.Resize(5,3);
      //--- initialization
      x9.Set(0,0,1);
      x9.Set(0,1,0);
      x9.Set(0,2,1);
      x9.Set(1,0,1);
      x9.Set(1,1,1);
      x9.Set(1,2,0);
      x9.Set(2,0,-1);
      x9.Set(2,1,1);
      x9.Set(2,2,0);
      x9.Set(3,0,-2);
      x9.Set(3,1,-1);
      x9.Set(3,2,1);
      x9.Set(4,0,-1);
      x9.Set(4,1,0);
      x9.Set(4,2,9);
      //--- check
      if(_spoil_scenario==73)
         Spoil_Matrix_By_Value(x9,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==74)
         Spoil_Matrix_By_Value(x9,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==75)
         Spoil_Matrix_By_Value(x9,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==76)
         Spoil_Matrix_By_Deleting_Row(x9);
      //--- check
      if(_spoil_scenario==77)
         Spoil_Matrix_By_Deleting_Col(x9);
      //--- function call
      CAlglib::SpearmanCorrM(x9,5,3,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- cross-vector short-form cov/corr are tested.
      x10.Resize(5,3);
      //--- initialization
      x10.Set(0,0,1);
      x10.Set(0,1,0);
      x10.Set(0,2,1);
      x10.Set(1,0,1);
      x10.Set(1,1,1);
      x10.Set(1,2,0);
      x10.Set(2,0,-1);
      x10.Set(2,1,1);
      x10.Set(2,2,0);
      x10.Set(3,0,-2);
      x10.Set(3,1,-1);
      x10.Set(3,2,1);
      x10.Set(4,0,-1);
      x10.Set(4,1,0);
      x10.Set(4,2,9);
      //--- check
      if(_spoil_scenario==78)
         Spoil_Matrix_By_Value(x10,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==79)
         Spoil_Matrix_By_Value(x10,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==80)
         Spoil_Matrix_By_Value(x10,CInfOrNaN::NegativeInfinity());
      //--- allocation
      y10.Resize(5,2);
      //--- initialization
      y10.Set(0,0,2);
      y10.Set(0,1,3);
      y10.Set(1,0,2);
      y10.Set(1,1,1);
      y10.Set(2,0,-1);
      y10.Set(2,1,6);
      y10.Set(3,0,-9);
      y10.Set(3,1,9);
      y10.Set(4,0,7);
      y10.Set(4,1,1);
      //--- check
      if(_spoil_scenario==81)
         Spoil_Matrix_By_Value(y10,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==82)
         Spoil_Matrix_By_Value(y10,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==83)
         Spoil_Matrix_By_Value(y10,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::CovM2(x10,y10,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      x11.Resize(5,3);
      //--- initialization
      x11.Set(0,0,1);
      x11.Set(0,1,0);
      x11.Set(0,2,1);
      x11.Set(1,0,1);
      x11.Set(1,1,1);
      x11.Set(1,2,0);
      x11.Set(2,0,-1);
      x11.Set(2,1,1);
      x11.Set(2,2,0);
      x11.Set(3,0,-2);
      x11.Set(3,1,-1);
      x11.Set(3,2,1);
      x11.Set(4,0,-1);
      x11.Set(4,1,0);
      x11.Set(4,2,9);
      //--- check
      if(_spoil_scenario==84)
         Spoil_Matrix_By_Value(x11,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==85)
         Spoil_Matrix_By_Value(x11,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==86)
         Spoil_Matrix_By_Value(x11,CInfOrNaN::NegativeInfinity());
      //--- allocation
      y11.Resize(5,2);
      //--- initialization
      y11.Set(0,0,2);
      y11.Set(0,1,3);
      y11.Set(1,0,2);
      y11.Set(1,1,1);
      y11.Set(2,0,-1);
      y11.Set(2,1,6);
      y11.Set(3,0,-9);
      y11.Set(3,1,9);
      y11.Set(4,0,7);
      y11.Set(4,1,1);
      //--- check
      if(_spoil_scenario==87)
         Spoil_Matrix_By_Value(y11,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==88)
         Spoil_Matrix_By_Value(y11,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==89)
         Spoil_Matrix_By_Value(y11,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::PearsonCorrM2(x11,y11,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      x12.Resize(5,3);
      //--- initialization
      x12.Set(0,0,1);
      x12.Set(0,1,0);
      x12.Set(0,2,1);
      x12.Set(1,0,1);
      x12.Set(1,1,1);
      x12.Set(1,2,0);
      x12.Set(2,0,-1);
      x12.Set(2,1,1);
      x12.Set(2,2,0);
      x12.Set(3,0,-2);
      x12.Set(3,1,-1);
      x12.Set(3,2,1);
      x12.Set(4,0,-1);
      x12.Set(4,1,0);
      x12.Set(4,2,9);
      //--- check
      if(_spoil_scenario==90)
         Spoil_Matrix_By_Value(x12,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==91)
         Spoil_Matrix_By_Value(x12,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==92)
         Spoil_Matrix_By_Value(x12,CInfOrNaN::NegativeInfinity());
      //--- allocation
      y12.Resize(5,2);
      //--- initialization
      y12.Set(0,0,2);
      y12.Set(0,1,3);
      y12.Set(1,0,2);
      y12.Set(1,1,1);
      y12.Set(2,0,-1);
      y12.Set(2,1,6);
      y12.Set(3,0,-9);
      y12.Set(3,1,9);
      y12.Set(4,0,7);
      y12.Set(4,1,1);
      //--- check
      if(_spoil_scenario==93)
         Spoil_Matrix_By_Value(y12,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==94)
         Spoil_Matrix_By_Value(y12,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==95)
         Spoil_Matrix_By_Value(y12,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::SpearmanCorrM2(x12,y12,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- cross-vector full-form cov/corr are tested.
      x13.Resize(5,3);
      //--- initialization
      x13.Set(0,0,1);
      x13.Set(0,1,0);
      x13.Set(0,2,1);
      x13.Set(1,0,1);
      x13.Set(1,1,1);
      x13.Set(1,2,0);
      x13.Set(2,0,-1);
      x13.Set(2,1,1);
      x13.Set(2,2,0);
      x13.Set(3,0,-2);
      x13.Set(3,1,-1);
      x13.Set(3,2,1);
      x13.Set(4,0,-1);
      x13.Set(4,1,0);
      x13.Set(4,2,9);
      //--- check
      if(_spoil_scenario==96)
         Spoil_Matrix_By_Value(x13,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==97)
         Spoil_Matrix_By_Value(x13,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==98)
         Spoil_Matrix_By_Value(x13,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==99)
         Spoil_Matrix_By_Deleting_Row(x13);
      //--- check
      if(_spoil_scenario==100)
         Spoil_Matrix_By_Deleting_Col(x13);
      //--- allocation
      y13.Resize(5,2);
      //--- initialization
      y13.Set(0,0,2);
      y13.Set(0,1,3);
      y13.Set(1,0,2);
      y13.Set(1,1,1);
      y13.Set(2,0,-1);
      y13.Set(2,1,6);
      y13.Set(3,0,-9);
      y13.Set(3,1,9);
      y13.Set(4,0,7);
      y13.Set(4,1,1);
      //--- check
      if(_spoil_scenario==101)
         Spoil_Matrix_By_Value(y13,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==102)
         Spoil_Matrix_By_Value(y13,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==103)
         Spoil_Matrix_By_Value(y13,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==104)
         Spoil_Matrix_By_Deleting_Row(y13);
      //--- check
      if(_spoil_scenario==105)
         Spoil_Matrix_By_Deleting_Col(y13);
      //--- function call
      CAlglib::CovM2(x13,y13,5,3,2,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      x14.Resize(5,3);
      //--- initialization
      x14.Set(0,0,1);
      x14.Set(0,1,0);
      x14.Set(0,2,1);
      x14.Set(1,0,1);
      x14.Set(1,1,1);
      x14.Set(1,2,0);
      x14.Set(2,0,-1);
      x14.Set(2,1,1);
      x14.Set(2,2,0);
      x14.Set(3,0,-2);
      x14.Set(3,1,-1);
      x14.Set(3,2,1);
      x14.Set(4,0,-1);
      x14.Set(4,1,0);
      x14.Set(4,2,9);
      //--- check
      if(_spoil_scenario==106)
         Spoil_Matrix_By_Value(x14,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==107)
         Spoil_Matrix_By_Value(x14,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==108)
         Spoil_Matrix_By_Value(x14,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==109)
         Spoil_Matrix_By_Deleting_Row(x14);
      //--- check
      if(_spoil_scenario==110)
         Spoil_Matrix_By_Deleting_Col(x14);
      //--- allocation
      y14.Resize(5,2);
      //--- initialization
      y14.Set(0,0,2);
      y14.Set(0,1,3);
      y14.Set(1,0,2);
      y14.Set(1,1,1);
      y14.Set(2,0,-1);
      y14.Set(2,1,6);
      y14.Set(3,0,-9);
      y14.Set(3,1,9);
      y14.Set(4,0,7);
      y14.Set(4,1,1);
      //--- check
      if(_spoil_scenario==111)
         Spoil_Matrix_By_Value(y14,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==112)
         Spoil_Matrix_By_Value(y14,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==113)
         Spoil_Matrix_By_Value(y14,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==114)
         Spoil_Matrix_By_Deleting_Row(y14);
      //--- check
      if(_spoil_scenario==115)
         Spoil_Matrix_By_Deleting_Col(y14);
      //--- function call
      CAlglib::PearsonCorrM2(x14,y14,5,3,2,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      x15.Resize(5,3);
      //--- initialization
      x15.Set(0,0,1);
      x15.Set(0,1,0);
      x15.Set(0,2,1);
      x15.Set(1,0,1);
      x15.Set(1,1,1);
      x15.Set(1,2,0);
      x15.Set(2,0,-1);
      x15.Set(2,1,1);
      x15.Set(2,2,0);
      x15.Set(3,0,-2);
      x15.Set(3,1,-1);
      x15.Set(3,2,1);
      x15.Set(4,0,-1);
      x15.Set(4,1,0);
      x15.Set(4,2,9);
      //--- check
      if(_spoil_scenario==116)
         Spoil_Matrix_By_Value(x15,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==117)
         Spoil_Matrix_By_Value(x15,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==118)
         Spoil_Matrix_By_Value(x15,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==119)
         Spoil_Matrix_By_Deleting_Row(x15);
      //--- check
      if(_spoil_scenario==120)
         Spoil_Matrix_By_Deleting_Col(x15);
      //--- allocation
      y15.Resize(5,2);
      //--- initialization
      y15.Set(0,0,2);
      y15.Set(0,1,3);
      y15.Set(1,0,2);
      y15.Set(1,1,1);
      y15.Set(2,0,-1);
      y15.Set(2,1,6);
      y15.Set(3,0,-9);
      y15.Set(3,1,9);
      y15.Set(4,0,7);
      y15.Set(4,1,1);
      //--- check
      if(_spoil_scenario==121)
         Spoil_Matrix_By_Value(y15,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==122)
         Spoil_Matrix_By_Value(y15,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==123)
         Spoil_Matrix_By_Value(y15,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==124)
         Spoil_Matrix_By_Deleting_Row(y15);
      //--- check
      if(_spoil_scenario==125)
         Spoil_Matrix_By_Deleting_Col(y15);
      //--- function call
      CAlglib::SpearmanCorrM2(x15,y15,5,3,2,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Real matrix inverse                                              |
//+------------------------------------------------------------------+
void TEST_MatInv_D_R1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
   CMatrixDouble tempmatrix;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
     {
      //--- allocation
      a.Resize(2,2);
      //--- initialization
      a.Set(0,0,1);
      a.Set(0,1,-1);
      a.Set(1,0,1);
      a.Set(1,1,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(a,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(a);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(a);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(a);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(a);
      //--- create variables
      int info;
      CMatInvReportShell rep;
      //--- function call
      CAlglib::RMatrixInverse(a,info,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(2,2);
      //--- initialization
      tempmatrix.Set(0,0,0.5);
      tempmatrix.Set(0,1,0.5);
      tempmatrix.Set(1,0,-0.5);
      tempmatrix.Set(1,1,0.5);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && Doc_Test_Real_Matrix(a,tempmatrix,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(rep.GetR1(),0.5,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(rep.GetRInf(),0.5,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Complex matrix inverse                                           |
//+------------------------------------------------------------------+
void TEST_MatInv_D_C1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixComplex a;
   CMatrixComplex tempmatrix;
   complex        tempcomplex1;
   complex        tempcomplex2;
   complex        cnan;
   complex        cpositiveinfinity;
   complex        cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
     {
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- initialization
      tempcomplex1.real=0;
      tempcomplex1.imag=1;
      tempcomplex2.real=0;
      tempcomplex2.imag=-0.5;
      //--- allocation
      a.Resize(2,2);
      //--- initialization
      a.Set(0,0,tempcomplex1);
      a.Set(0,1,-1);
      a.Set(1,0,tempcomplex1);
      a.Set(1,1,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(a,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(a,cnegativeinfinity);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(a);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(a);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(a);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(a);
      //--- create variables
      int info;
      CMatInvReportShell rep;
      //--- function call
      CAlglib::CMatrixInverse(a,info,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(2,2);
      //--- initialization
      tempmatrix.Set(0,0,tempcomplex2);
      tempmatrix.Set(0,1,tempcomplex2);
      tempmatrix.Set(1,0,-0.5);
      tempmatrix.Set(1,1,0.5);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && Doc_Test_Complex_Matrix(a,tempmatrix,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(rep.GetR1(),0.5,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(rep.GetRInf(),0.5,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| SPD matrix inverse                                               |
//+------------------------------------------------------------------+
void TEST_MatInv_D_SPD1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
   CMatrixDouble tempmatrix;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
     {
      //--- allocation
      a.Resize(2,2);
      //--- initialization
      a.Set(0,0,2);
      a.Set(0,1,1);
      a.Set(1,0,1);
      a.Set(1,1,2);
      //--- check
      do
        {
         if(_spoil_scenario==0)
            Spoil_Matrix_By_Value(a,CInfOrNaN::NaN());
         //--- check
         if(_spoil_scenario==1)
            Spoil_Matrix_By_Value(a,CInfOrNaN::PositiveInfinity());
         //--- check
         if(_spoil_scenario==2)
            Spoil_Matrix_By_Value(a,CInfOrNaN::NegativeInfinity());
        }
      while(_spoil_scenario>=0 && _spoil_scenario<=2 && CApServ::IsFiniteRTrMatrix(a,2,false));
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(a);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(a);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(a);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(a);
      //--- create variables
      int info;
      CMatInvReportShell rep;
      //--- function call
      CAlglib::SPDMatrixInverse(a,info,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(2,2);
      //--- initialization
      tempmatrix.Set(0,0,0.666666);
      tempmatrix.Set(0,1,-0.333333);
      tempmatrix.Set(1,0,-0.333333);
      tempmatrix.Set(1,1,0.666666);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && Doc_Test_Real_Matrix(a,tempmatrix,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| HPD matrix inverse                                               |
//+------------------------------------------------------------------+
void TEST_MatInv_D_HPD1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixComplex a;
   CMatrixComplex tempmatrix;
   complex    cnan;
   complex    cpositiveinfinity;
   complex    cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
     {
      //--- allocation
      a.Resize(2,2);
      //--- initialization
      a.Set(0,0,2);
      a.Set(0,1,1);
      a.Set(1,0,1);
      a.Set(1,1,2);
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(a,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(a,cnegativeinfinity);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(a);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(a);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(a);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(a);
      //--- create variables
      int info;
      CMatInvReportShell rep;
      //--- function call
      CAlglib::HPDMatrixInverse(a,info,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      tempmatrix.Resize(2,2);
      //--- initialization
      tempmatrix.Set(0,0,0.666666);
      tempmatrix.Set(0,1,-0.333333);
      tempmatrix.Set(1,0,-0.333333);
      tempmatrix.Set(1,1,0.666666);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && Doc_Test_Complex_Matrix(a,tempmatrix,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Real matrix inverse: singular matrix                             |
//+------------------------------------------------------------------+
void TEST_MatInv_T_R1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
//--- allocation
   a.Resize(2,2);
//--- initialization
   a.Set(0,0,1);
   a.Set(0,1,-1);
   a.Set(1,0,-2);
   a.Set(1,1,2);
//--- create variables
   int info;
   CMatInvReportShell rep;
//--- function call
   CAlglib::RMatrixInverse(a,info,rep);
//--- handling exceptions
   Func_spoil_scenario(_spoil_scenario,_TestResult);
//--- check result
   _TestResult=_TestResult && Doc_Test_Int(info,-3);
   _TestResult=_TestResult && Doc_Test_Real(rep.GetR1(),0.0,0.00005);
   _TestResult=_TestResult && Doc_Test_Real(rep.GetRInf(),0.0,0.00005);
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Complex matrix inverse: singular matrix                          |
//+------------------------------------------------------------------+
void TEST_MatInv_T_C1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixComplex a;
   complex        tempcomplex1;
   complex        tempcomplex2;

   tempcomplex1.real=0;
   tempcomplex1.imag=1;
   tempcomplex2.real=0;
   tempcomplex2.imag=-1;
//--- allocation
   a.Resize(2,2);
//--- initialization
   a.Set(0,0,tempcomplex1);
   a.Set(0,1,tempcomplex2);
   a.Set(1,0,-2);
   a.Set(1,1,2);
//--- create variables
   int info;
   CMatInvReportShell rep;
//--- function call
   CAlglib::CMatrixInverse(a,info,rep);
//--- handling exceptions
   Func_spoil_scenario(_spoil_scenario,_TestResult);
//--- check result
   _TestResult=_TestResult && Doc_Test_Int(info,-3);
   _TestResult=_TestResult && Doc_Test_Real(rep.GetR1(),0.0,0.00005);
   _TestResult=_TestResult && Doc_Test_Real(rep.GetRInf(),0.0,0.00005);
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Attempt to use SPD function on nonsymmetrix matrix               |
//+------------------------------------------------------------------+
void TEST_MatInv_E_SPD1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
//--- allocation
   a.Resize(2,2);
//--- initialization
   a.Set(0,0,1);
   a.Set(0,1,0);
   a.Set(1,0,1);
   a.Set(1,1,1);
//--- create variables
   int info;
   CMatInvReportShell rep;
//--- function call
   CAlglib::SPDMatrixInverse(a,info,rep);
//--- handling exceptions
   if(Func_spoil_scenario(_spoil_scenario,_TestResult))
      _TestResult=false;
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Attempt to use SPD function on nonsymmetrix matrix               |
//+------------------------------------------------------------------+
void TEST_MatInv_E_HPD1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixComplex a;
//--- allocation
   a.Resize(2,2);
//--- initialization
   a.Set(0,0,1);
   a.Set(0,1,0);
   a.Set(1,0,1);
   a.Set(1,1,1);
//--- create variables
   int info;
   CMatInvReportShell rep;
//--- function call
   CAlglib::HPDMatrixInverse(a,info,rep);
//--- handling exceptions
   if(Func_spoil_scenario(_spoil_scenario,_TestResult))
      _TestResult=false;
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization by CG                                     |
//+------------------------------------------------------------------+
void TEST_MinCG_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   CObject             obj;
   CNDimensional_Grad1 fgrad;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- with nonlinear conjugate gradient method.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- create a variable
      double epsg=0.0000000001;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      //--- create a variable
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      //--- create a variable
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      int maxits=0;
      //--- create variables
      CMinCGStateShell state;
      CMinCGReportShell rep;
      //--- function call
      CAlglib::MinCGCreate(x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGOptimize(state,fgrad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization with additional settings and restarts     |
//+------------------------------------------------------------------+
void TEST_MinCG_D_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   CObject             obj;
   CNDimensional_Grad1 fgrad;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<18; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- with nonlinear conjugate gradient method.
      //--- Several advanced techniques are demonstrated:
      //--- * upper limit on step size
      //--- * restart from new point
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=0.0000000001;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      double stpmax=0.1;
      //--- check
      if(_spoil_scenario==12)
         stpmax=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==13)
         stpmax=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==14)
         stpmax=CInfOrNaN::NegativeInfinity();
      int maxits=0;
      //--- create variables
      CMinCGStateShell state;
      CMinCGReportShell rep;
      //--- first run
      CAlglib::MinCGCreate(x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGSetStpMax(state,stpmax);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGOptimize(state,fgrad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      //--- second run - algorithm is restarted with mincgrestartfrom()
      ArrayResize(x,2);
      //--- initialization
      x[0]=10;
      x[1]=10;
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::MinCGRestartFrom(state,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGOptimize(state,fgrad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization by CG with numerical differentiation      |
//+------------------------------------------------------------------+
void TEST_MinCG_NumDiff(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   CObject             obj;
   CNDimensional_Func1 ffunc;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- using numerical differentiation to calculate gradient.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=0.0000000001;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      double diffstep=1.0e-6;
      //--- check
      if(_spoil_scenario==12)
         diffstep=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==13)
         diffstep=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==14)
         diffstep=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinCGStateShell state;
      CMinCGReportShell rep;
      //--- function call
      CAlglib::MinCGCreateF(x,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGOptimize(state,ffunc,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization by CG, function with singularities        |
//+------------------------------------------------------------------+
void TEST_MinCG_FTRIM(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double                x[];
   double                temparray[];
   CObject               obj;
   CNDimensional_S1_Grad fs1grad;
   CNDimensional_Rep     frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x)=(1+x)^(-0.2) + (1-x)^(-0.3) + 1000*x.
      //--- This function has singularities at the boundary of the [-1,+1],but technique called
      //--- "function trimming" allows us to solve this optimization problem.
      //--- See http://www.CAlglib::net/optimization/tipsandtricks.php#ftrimming for more information
      //--- on this subject.
      ArrayResize(x,1);
      //--- initialization
      x[0]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=1.0e-6;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinCGStateShell state;
      CMinCGReportShell rep;
      //--- function call
      CAlglib::MinCGCreate(x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGOptimize(state,fs1grad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinCGResults(state,x,rep);
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=-0.99917305;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.000005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization with bound constraints                    |
//+------------------------------------------------------------------+
void TEST_MinBLEIC_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   double              bndl[];
   double              bndu[];
   CObject             obj;
   CNDimensional_Grad1 fgrad;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<22; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- subject to bound constraints -1<=x<=+1,-1<=y<=+1,using BLEIC optimizer.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- allocation
      ArrayResize(bndl,2);
      //--- initialization
      bndl[0]=-1;
      bndl[1]=-1;
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(bndl,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(bndl);
      //--- allocation
      ArrayResize(bndu,2);
      //--- initialization
      bndu[0]=1;
      bndu[1]=1;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(bndu,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Deleting_Element(bndu);
      CMinBLEICStateShell state;
      CMinBLEICReportShell rep;
      //--- These variables define stopping conditions for the underlying CG algorithm.
      //--- They should be stringent enough in order to guarantee overall stability
      //--- of the outer iterations.
      //--- We use very simple condition - |g|<=epsg
      double epsg=0.000001;
      //--- check
      if(_spoil_scenario==7)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==8)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==9)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==10)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==11)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==12)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==13)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==14)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==15)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- These variables define stopping conditions for the outer iterations:
      //--- * epso controls convergence of outer iterations;algorithm will stop
      //---   when difference between solutions of subsequent unconstrained problems
      //---   will be less than 0.0001
      //--- * epsi controls amount of infeasibility allowed in the final solution
      double epso=0.00001;
      //--- check
      if(_spoil_scenario==16)
         epso=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==17)
         epso=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==18)
         epso=CInfOrNaN::NegativeInfinity();
      double epsi=0.00001;
      //--- check
      if(_spoil_scenario==19)
         epsi=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==20)
         epsi=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==21)
         epsi=CInfOrNaN::NegativeInfinity();
      //--- Now we are ready to actually optimize something:
      //--- * first we create optimizer
      //--- * we add boundary constraints
      //--- * we tune stopping conditions
      //--- * and,finally,optimize and obtain results...
      CAlglib::MinBLEICCreate(x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetInnerCond(state,epsg,epsf,epsx);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetOuterCond(state,epso,epsi);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICOptimize(state,fgrad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- ...and evaluate these results
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-1;
      temparray[1]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization with linear inequality constraints        |
//+------------------------------------------------------------------+
void TEST_MinBLEIC_D_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   int                 ct[];
   CMatrixDouble       c;
   CObject             obj;
   CNDimensional_Grad1 fgrad;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<24; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- subject to inequality constraints:
      //--- * x>=2 (posed as general linear constraint),
      //--- * x+y>=6
      //--- using BLEIC optimizer.
      ArrayResize(x,2);
      //--- initialization
      x[0]=5;
      x[1]=5;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- allocation
      c.Resize(2,3);
      //--- initialization
      c.Set(0,0,1);
      c.Set(0,1,0);
      c.Set(0,2,2);
      c.Set(1,0,1);
      c.Set(1,1,1);
      c.Set(1,2,6);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Value(c,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Value(c,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Value(c,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Row(c);
      //--- check
      if(_spoil_scenario==7)
         Spoil_Matrix_By_Deleting_Col(c);
      //--- allocation
      ArrayResize(ct,2);
      //--- initialization
      ct[0]=1;
      ct[1]=1;
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(ct);
      //--- create variables
      CMinBLEICStateShell state;
      CMinBLEICReportShell rep;
      //--- These variables define stopping conditions for the underlying CG algorithm.
      //--- They should be stringent enough in order to guarantee overall stability
      //--- of the outer iterations.
      //--- We use very simple condition - |g|<=epsg
      double epsg=0.000001;
      //--- check
      if(_spoil_scenario==9)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==12)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==13)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==14)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==15)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==16)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==17)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- These variables define stopping conditions for the outer iterations:
      //--- * epso controls convergence of outer iterations;algorithm will stop
      //---   when difference between solutions of subsequent unconstrained problems
      //---   will be less than 0.0001
      //--- * epsi controls amount of infeasibility allowed in the final solution
      double epso=0.00001;
      //--- check
      if(_spoil_scenario==18)
         epso=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==19)
         epso=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==20)
         epso=CInfOrNaN::NegativeInfinity();
      double epsi=0.00001;
      //--- check
      if(_spoil_scenario==21)
         epsi=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==22)
         epsi=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==23)
         epsi=CInfOrNaN::NegativeInfinity();
      //--- Now we are ready to actually optimize something:
      //--- * first we create optimizer
      //--- * we add linear constraints
      //--- * we tune stopping conditions
      //--- * and,finally,optimize and obtain results...
      CAlglib::MinBLEICCreate(x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetLC(state,c,ct);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetInnerCond(state,epsg,epsf,epsx);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetOuterCond(state,epso,epsi);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICOptimize(state,fgrad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- ...and evaluate these results
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=2;
      temparray[1]=4;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization with bound constraints and numerical      |
//| differentiation                                                  |
//+------------------------------------------------------------------+
void TEST_MinBLEIC_NumDiff(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   double              bndl[];
   double              bndu[];
   CObject             obj;
   CNDimensional_Func1 ffunc;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<25; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- subject to bound constraints -1<=x<=+1,-1<=y<=+1,using BLEIC optimizer.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- allocation
      ArrayResize(bndl,2);
      //--- initialization
      bndl[0]=-1;
      bndl[1]=-1;
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(bndl,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(bndl);
      //--- allocation
      ArrayResize(bndu,2);
      //--- initialization
      bndu[0]=1;
      bndu[1]=1;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(bndu,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Deleting_Element(bndu);
      CMinBLEICStateShell state;
      CMinBLEICReportShell rep;
      //--- These variables define stopping conditions for the underlying CG algorithm.
      //--- They should be stringent enough in order to guarantee overall stability
      //--- of the outer iterations.
      //--- We use very simple condition - |g|<=epsg
      double epsg=0.000001;
      //--- check
      if(_spoil_scenario==7)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==8)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==9)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==10)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==11)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==12)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==13)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==14)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==15)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- These variables define stopping conditions for the outer iterations:
      //--- * epso controls convergence of outer iterations;algorithm will stop
      //---   when difference between solutions of subsequent unconstrained problems
      //---   will be less than 0.0001
      //--- * epsi controls amount of infeasibility allowed in the final solution
      double epso=0.00001;
      //--- check
      if(_spoil_scenario==16)
         epso=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==17)
         epso=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==18)
         epso=CInfOrNaN::NegativeInfinity();
      double epsi=0.00001;
      //--- check
      if(_spoil_scenario==19)
         epsi=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==20)
         epsi=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==21)
         epsi=CInfOrNaN::NegativeInfinity();
      //--- This variable contains differentiation step
      double diffstep=1.0e-6;
      //--- check
      if(_spoil_scenario==22)
         diffstep=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==23)
         diffstep=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==24)
         diffstep=CInfOrNaN::NegativeInfinity();
      //--- Now we are ready to actually optimize something:
      //--- * first we create optimizer
      //--- * we add boundary constraints
      //--- * we tune stopping conditions
      //--- * and,finally,optimize and obtain results...
      CAlglib::MinBLEICCreateF(x,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetInnerCond(state,epsg,epsf,epsx);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetOuterCond(state,epso,epsi);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICOptimize(state,ffunc,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- ...and evaluate these results
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-1;
      temparray[1]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization by BLEIC, function with singularities     |
//+------------------------------------------------------------------+
void TEST_MinBLEIC_FTRIM(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double                x[];
   double                temparray[];
   CObject               obj;
   CNDimensional_S1_Grad fs1grad;
   CNDimensional_Rep     frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<18; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x)=(1+x)^(-0.2) + (1-x)^(-0.3) + 1000*x.
      //--- This function is undefined outside of (-1,+1) and has singularities at x=-1 and x=+1.
      //--- Special technique called "function trimming" allows us to solve this optimization problem
      //--- - withusing boundary constraints!
      //--- See http://www.CAlglib::net/optimization/tipsandtricks.php#ftrimming for more information
      //--- on this subject.
      ArrayResize(x,1);
      //--- initialization
      x[0]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=1.0e-6;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      double epso=1.0e-6;
      //--- check
      if(_spoil_scenario==12)
         epso=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==13)
         epso=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==14)
         epso=CInfOrNaN::NegativeInfinity();
      double epsi=1.0e-6;
      //--- check
      if(_spoil_scenario==15)
         epsi=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==16)
         epsi=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==17)
         epsi=CInfOrNaN::NegativeInfinity();
      //--- create variables
      CMinBLEICStateShell state;
      CMinBLEICReportShell rep;
      //--- function call
      CAlglib::MinBLEICCreate(x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetInnerCond(state,epsg,epsf,epsx);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICSetOuterCond(state,epso,epsi);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICOptimize(state,fs1grad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinBLEICResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=-0.99917305;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.000005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple unconstrained MCPD model (no entry/exit states)           |
//+------------------------------------------------------------------+
void TEST_MCPD_Simple1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble track0;
   CMatrixDouble track1;
   CMatrixDouble tempmatrix;
   CMatrixDouble p;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- The very simple MCPD example
      //--- We have a loan portfolio. Our loans can be in one of two states:
      //---*normal loans ("good" ones)
      //---*past due loans ("bad" ones)
      //--- We assume that:
      //--- * loans can transition from any state to any other state. In
      //---   particular,past due loan can become "good" one at any moment
      //---   with same (fixed) probability. Not realistic,but it is toy example :)
      //--- * portfolio size does not change over time
      //
      //--- Thus,we have following model
      //---     state_new=P*state_old
      //--- where
      //---         ( p00  p01 )
      //---     P = (          )
      //---         ( p10  p11 )
      //--- We want to model transitions between these two states using MCPD
      //--- approach (Markov Chains for Proportional/Population Data),i.e.
      //--- to restore hidden transition matrix P using actual portfolio data.
      //--- We have:
      //--- * poportional data,i.e. proportion of loans in the normal and past
      //---   due states (not portfolio size measured in some currency,although
      //---   it is possible to work with population data too)
      //--- * two tracks,i.e. two sequences which describe portfolio
      //---   evolution from two different starting states: [1,0] (all loans
      //---   are "good") and [0.8,0.2] (only 80% of portfolio is in the "good"
      //---   state)
      CMCPDStateShell s;
      CMCPDReportShell rep;
      //--- allocation
      track0.Resize(5,2);
      //--- initialization
      track0.Set(0,0,1);
      track0.Set(0,1,0);
      track0.Set(1,0,0.95);
      track0.Set(1,1,0.05);
      track0.Set(2,0,0.9275);
      track0.Set(2,1,0.0725);
      track0.Set(3,0,0.91738);
      track0.Set(3,1,0.08263);
      track0.Set(4,0,0.91282);
      track0.Set(4,1,0.08718);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(track0,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(track0,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(track0,CInfOrNaN::NegativeInfinity());
      //--- allocation
      track1.Resize(4,2);
      //--- initialization
      track1.Set(0,0,0.8);
      track1.Set(0,1,0.2);
      track1.Set(1,0,0.86);
      track1.Set(1,1,0.14);
      track1.Set(2,0,0.887);
      track1.Set(2,1,0.113);
      track1.Set(3,0,0.89915);
      track1.Set(3,1,0.10085);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Value(track1,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Value(track1,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Value(track1,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::MCPDCreate(2,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDAddTrack(s,track0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDAddTrack(s,track1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDSolve(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDResults(s,p,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Hidden matrix P is equal to
      //---         ( 0.95  0.50 )
      //---         (            )
      //---         ( 0.05  0.50 )
      //--- which means that "good" loans can become "bad" with 5% probability,
      //--- while "bad" loans will return to good state with 50% probability.
      tempmatrix.Resize(2,2);
      //--- initialization
      tempmatrix.Set(0,0,0.95);
      tempmatrix.Set(0,1,0.5);
      tempmatrix.Set(1,0,0.05);
      tempmatrix.Set(1,1,0.5);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(p,tempmatrix,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple MCPD model (no entry/exit states) with equality           |
//| constraints                                                      |
//+------------------------------------------------------------------+
void TEST_MCPD_Simple2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble track0;
   CMatrixDouble track1;
   CMatrixDouble tempmatrix;
   CMatrixDouble p;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- Simple MCPD example
      //--- We have a loan portfolio. Our loans can be in one of three states:
      //--- * normal loans
      //--- * past due loans
      //--- * charged off loans
      //--- We assume that:
      //--- * normal loan can stay normal or become past due (but not charged off)
      //--- * past due loan can stay past due,become normal or charged off
      //--- * charged off loan will stay charged off for the rest of eternity
      //--- * portfolio size does not change over time
      //--- Not realistic,but it is toy example :)
      //--- Thus,we have following model
      //---     state_new=P*state_old
      //--- where
      //---         ( p00  p01    )
      //---     P = ( p10  p11    )
      //---         (      p21  1 )
      //--- i.e. four elements of P are known a priori.
      //--- Although it is possible (given enough data) to In order to enforce
      //--- this property we set equality constraints on these elements.
      //--- We want to model transitions between these two states using MCPD
      //--- approach (Markov Chains for Proportional/Population Data),i.e.
      //--- to restore hidden transition matrix P using actual portfolio data.
      //--- We have:
      //--- * poportional data,i.e. proportion of loans in the current and past
      //---   due states (not portfolio size measured in some currency,although
      //---   it is possible to work with population data too)
      //--- * two tracks,i.e. two sequences which describe portfolio
      //---   evolution from two different starting states: [1,0,0] (all loans
      //---   are "good") and [0.8,0.2,0.0] (only 80% of portfolio is in the "good"
      //---   state)
      CMCPDStateShell s;
      CMCPDReportShell rep;
      //--- allocation
      track0.Resize(5,3);
      //--- initialization
      track0.Set(0,0,1);
      track0.Set(0,1,0);
      track0.Set(0,2,0);
      track0.Set(1,0,0.95);
      track0.Set(1,1,0.05);
      track0.Set(1,2,0);
      track0.Set(2,0,0.9275);
      track0.Set(2,1,0.06);
      track0.Set(2,2,0.0125);
      track0.Set(3,0,0.911125);
      track0.Set(3,1,0.061375);
      track0.Set(3,2,0.0275);
      track0.Set(4,0,0.896256);
      track0.Set(4,1,0.0609);
      track0.Set(4,2,0.042844);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(track0,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(track0,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(track0,CInfOrNaN::NegativeInfinity());
      //--- allocation
      track1.Resize(5,3);
      //--- initialization
      track1.Set(0,0,0.8);
      track1.Set(0,1,0.2);
      track1.Set(0,2,0);
      track1.Set(1,0,0.86);
      track1.Set(1,1,0.09);
      track1.Set(1,2,0.05);
      track1.Set(2,0,0.862);
      track1.Set(2,1,0.0655);
      track1.Set(2,2,0.0725);
      track1.Set(3,0,0.85165);
      track1.Set(3,1,0.059475);
      track1.Set(3,2,0.088875);
      track1.Set(4,0,0.838805);
      track1.Set(4,1,0.057451);
      track1.Set(4,2,0.103744);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Value(track1,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Value(track1,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Value(track1,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::MCPDCreate(3,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDAddTrack(s,track0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDAddTrack(s,track1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDAddEC(s,0,2,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDAddEC(s,1,2,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDAddEC(s,2,2,1.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDAddEC(s,2,0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDSolve(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MCPDResults(s,p,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Hidden matrix P is equal to
      //---         ( 0.95 0.50      )
      //---         ( 0.05 0.25      )
      //---         (      0.25 1.00 )
      //--- which means that "good" loans can become past due with 5% probability,
      //--- while past due loans will become charged off with 25% probability or
      //--- return back to normal state with 50% probability.
      tempmatrix.Resize(3,3);
      //--- initialization
      tempmatrix.Set(0,0,0.95);
      tempmatrix.Set(0,1,0.5);
      tempmatrix.Set(0,2,0);
      tempmatrix.Set(1,0,0.05);
      tempmatrix.Set(1,1,0.25);
      tempmatrix.Set(1,2,0);
      tempmatrix.Set(2,0,0);
      tempmatrix.Set(2,1,0.25);
      tempmatrix.Set(2,2,1);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Matrix(p,tempmatrix,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization by L-BFGS                                 |
//+------------------------------------------------------------------+
void TEST_MinLBFGS_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   CObject             obj;
   CNDimensional_Grad1 fgrad;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- using LBFGS method.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=0.0000000001;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinLBFGSStateShell state;
      CMinLBFGSReportShell rep;
      //--- function call
      CAlglib::MinLBFGSCreate(1,x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSOptimize(state,fgrad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization with additional settings and restarts     |
//+------------------------------------------------------------------+
void TEST_MinLBFGS_D_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   CObject             obj;
   CNDimensional_Grad1 fgrad;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<18; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- using LBFGS method.
      //--- Several advanced techniques are demonstrated:
      //--- * upper limit on step size
      //--- * restart from new point
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=0.0000000001;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      double stpmax=0.1;
      //--- check
      if(_spoil_scenario==12)
         stpmax=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==13)
         stpmax=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==14)
         stpmax=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinLBFGSStateShell state;
      CMinLBFGSReportShell rep;
      //--- first run
      CAlglib::MinLBFGSCreate(1,x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSSetStpMax(state,stpmax);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSOptimize(state,fgrad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      //--- second run - algorithm is restarted
      ArrayResize(x,2);
      //--- initialization
      x[0]=10;
      x[1]=10;
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::MinLBFGSRestartFrom(state,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSOptimize(state,fgrad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization by L-BFGS with numerical differentiation  |
//+------------------------------------------------------------------+
void TEST_MinLBFGS_NumDiff(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              temparray[];
   CObject             obj;
   CNDimensional_Func1 ffunc;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x,y)=100*(x+3)^4+(y-3)^4
      //--- using numerical differentiation to calculate gradient.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=0.0000000001;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      double diffstep=1.0e-6;
      //--- check
      if(_spoil_scenario==12)
         diffstep=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==13)
         diffstep=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==14)
         diffstep=CInfOrNaN::NegativeInfinity();
      int maxits=0;
      CMinLBFGSStateShell state;
      CMinLBFGSReportShell rep;
      //--- function call
      CAlglib::MinLBFGSCreateF(1,x,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSOptimize(state,ffunc,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization by LBFGS, function with singularities     |
//+------------------------------------------------------------------+
void TEST_MinLBFGS_FTRIM(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double                x[];
   double                temparray[];
   CObject               obj;
   CNDimensional_S1_Grad fs1grad;
   CNDimensional_Rep     frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of f(x)=(1+x)^(-0.2) + (1-x)^(-0.3) + 1000*x.
      //--- This function has singularities at the boundary of the [-1,+1],but technique called
      //--- "function trimming" allows us to solve this optimization problem.
      //--- See http://www.CAlglib::net/optimization/tipsandtricks.php#ftrimming for more information
      //--- on this subject.
      ArrayResize(x,1);
      //--- initialization
      x[0]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=1.0e-6;
      //--- check
      if(_spoil_scenario==3)
         epsg=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsg=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsg=CInfOrNaN::NegativeInfinity();
      double epsf=0;
      //--- check
      if(_spoil_scenario==6)
         epsf=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsf=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsf=CInfOrNaN::NegativeInfinity();
      double epsx=0;
      //--- check
      if(_spoil_scenario==9)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinLBFGSStateShell state;
      CMinLBFGSReportShell rep;
      //--- function call
      CAlglib::MinLBFGSCreate(1,x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSOptimize(state,fs1grad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLBFGSResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=-0.99917305;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.000005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Solving y'=-y with ODE solver                                    |
//+------------------------------------------------------------------+
void TEST_ODESolver_D1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double        y[];
   double        x[];
   int           m;
   double        xtbl[];
   double        temparray[];
   CMatrixDouble ytbl;
   CMatrixDouble tempmatrix;
   CObject       obj;
   CNDimensional_ODE_Function_1_Dif fode;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,1);
      //--- initialization
      y[0]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(x,4);
      //--- initialization
      x[0]=0;
      x[1]=1;
      x[2]=2;
      x[3]=3;
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double eps=0.00001;
      //--- check
      if(_spoil_scenario==7)
         eps=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==8)
         eps=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==9)
         eps=CInfOrNaN::NegativeInfinity();
      double h=0;
      //--- check
      if(_spoil_scenario==10)
         h=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==11)
         h=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==12)
         h=CInfOrNaN::NegativeInfinity();
      //--- create variables
      CODESolverStateShell s;
      CODESolverReportShell rep;
      //--- function call
      CAlglib::ODESolverRKCK(y,x,eps,h,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::ODESolverSolve(s,fode,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::ODESolverResults(s,m,xtbl,ytbl,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=0;
      temparray[1]=1;
      temparray[2]=2;
      temparray[3]=3;
      //--- allocation
      tempmatrix.Resize(4,1);
      //--- initialization
      tempmatrix.Set(0,0,1);
      tempmatrix.Set(1,0,0.367);
      tempmatrix.Set(2,0,0.135);
      tempmatrix.Set(3,0,0.05);
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(m,4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(xtbl,temparray,0.005);
      _TestResult=_TestResult && Doc_Test_Real_Matrix(ytbl,tempmatrix,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Complex FFT: simple example                                      |
//+------------------------------------------------------------------+
void TEST_FFT_Complex_D1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex z[];
   complex temparray[];
   complex tempcomplex1;
   complex tempcomplex2;
   complex cnan;
   complex cpositiveinfinity;
   complex cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- first we demonstrate forward FFT:
      //--- [1i,1i,1i,1i] is converted to [4i,0,0,0]
      tempcomplex1.real=0;
      tempcomplex1.imag=1;
      tempcomplex2.real=0;
      tempcomplex2.imag=4;
      //--- allocation
      ArrayResize(z,4);
      //--- initialization
      z[0]=tempcomplex1;
      z[1]=tempcomplex1;
      z[2]=tempcomplex1;
      z[3]=tempcomplex1;
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(z,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(z,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(z,cnegativeinfinity);
      //--- function call
      CAlglib::FFTC1D(z);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=tempcomplex2;
      temparray[1]=0;
      temparray[2]=0;
      temparray[3]=0;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex_Vector(z,temparray,0.0001);
      //--- now we convert [4i,0,0,0] back to [1i,1i,1i,1i]
      //--- with backward FFT
      CAlglib::FFTC1DInv(z);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=tempcomplex1;
      temparray[1]=tempcomplex1;
      temparray[2]=tempcomplex1;
      temparray[3]=tempcomplex1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex_Vector(z,temparray,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Complex FFT: advanced example                                    |
//+------------------------------------------------------------------+
void TEST_FFT_Complex_D2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex z[];
   complex temparray[];
   complex tempcomplex1;
   complex tempcomplex2;
   complex tempcomplex3;
   complex cnan;
   complex cpositiveinfinity;
   complex cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- first we demonstrate forward FFT:
      //--- [0,1,0,1i] is converted to [1+1i,-1-1i,-1-1i,1+1i]
      tempcomplex1.real=0;
      tempcomplex1.imag=1;
      //--- allocation
      ArrayResize(z,4);
      //--- initialization
      z[0]=0;
      z[1]=1;
      z[2]=0;
      z[3]=tempcomplex1;
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(z,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(z,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(z,cnegativeinfinity);
      //--- function call
      CAlglib::FFTC1D(z);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex2.real=1;
      tempcomplex2.imag=1;
      tempcomplex3.real=-1;
      tempcomplex3.imag=-1;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=tempcomplex2;
      temparray[1]=tempcomplex3;
      temparray[2]=tempcomplex3;
      temparray[3]=tempcomplex2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex_Vector(z,temparray,0.0001);
      //--- now we convert result back with backward FFT
      CAlglib::FFTC1DInv(z);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=0;
      temparray[1]=1;
      temparray[2]=0;
      temparray[3]=tempcomplex1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex_Vector(z,temparray,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Real FFT: simple example                                         |
//+------------------------------------------------------------------+
void TEST_FFT_Real_D1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double  x[];
   complex tempcarray[];
   double  temparray[];
   complex f[];
   double  x2[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- first we demonstrate forward FFT:
      //--- [1,1,1,1] is converted to [4,0,0,0]
      ArrayResize(x,4);
      //--- initialization
      x[0]=1;
      x[1]=1;
      x[2]=1;
      x[3]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::FFTR1D(x,f);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(tempcarray,4);
      //--- initialization
      tempcarray[0]=4;
      tempcarray[1]=0;
      tempcarray[2]=0;
      tempcarray[3]=0;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex_Vector(f,tempcarray,0.0001);
      //--- now we convert [4,0,0,0] back to [1,1,1,1]
      //--- with backward FFT
      CAlglib::FFTR1DInv(f,x2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=1;
      temparray[1]=1;
      temparray[2]=1;
      temparray[3]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x2,temparray,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Real FFT: advanced example                                       |
//+------------------------------------------------------------------+
void TEST_FFT_Real_D2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double  x[];
   double  temparray[];
   complex tempcarray[];
   complex f[];
   double  x2[];
   complex tempcomplex1;
   complex tempcomplex2;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- first we demonstrate forward FFT:
      //--- [1,2,3,4] is converted to [10,-2+2i,-2,-2-2i]
      //--- note that output array is self-adjoint:
      //--- * f[0]=conj(f[0])
      //--- * f[1]=conj(f[3])
      //--- * f[2]=conj(f[2])
      ArrayResize(x,4);
      //--- initialization
      x[0]=1;
      x[1]=2;
      x[2]=3;
      x[3]=4;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::FFTR1D(x,f);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex1.real=-2;
      tempcomplex1.imag=2;
      tempcomplex2.real=-2;
      tempcomplex2.imag=-2;
      //--- allocation
      ArrayResize(tempcarray,4);
      //--- initialization
      tempcarray[0]=10;
      tempcarray[1]=tempcomplex1;
      tempcarray[2]=-2;
      tempcarray[3]=tempcomplex2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex_Vector(f,tempcarray,0.0001);
      //--- now we convert [10,-2+2i,-2,-2-2i] back to [1,2,3,4]
      CAlglib::FFTR1DInv(f,x2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=1;
      temparray[1]=2;
      temparray[2]=3;
      temparray[3]=4;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x2,temparray,0.0001);
      //--- remember that F is self-adjoint? It means that we can pass just half
      //--- (slightly larger than half) of F to inverse real FFT and still get our result.
      //--- I.e. instead [10,-2+2i,-2,-2-2i] we pass just [10,-2+2i,-2] and everything works!
      //--- NOTE: in this case we should explicitly pass array length (which is 4) to ALGLIB;
      //--- if not,it will automatically use array length to determine FFT size and
      //--- will erroneously make half-length FFT.
      ArrayResize(f,3);
      //--- initialization
      f[0]=10;
      f[1]=tempcomplex1;
      f[2]=-2;
      //--- function call
      CAlglib::FFTR1DInv(f,4,x2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=1;
      temparray[1]=2;
      temparray[2]=3;
      temparray[3]=4;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x2,temparray,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| error detection in backward FFT                                  |
//+------------------------------------------------------------------+
void TEST_FFT_Complex_E1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex tempcomplex1;
   complex tempcomplex2;
   complex z[];
   complex temparray[];
   complex cnan;
   complex cpositiveinfinity;
   complex cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(z,4);
      //--- initialization
      z[0]=0;
      z[1]=2;
      z[2]=0;
      z[3]=-2;
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(z,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(z,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(z,cnegativeinfinity);
      //--- function call
      CAlglib::FFTC1DInv(z);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex1.real=0;
      tempcomplex1.imag=1;
      tempcomplex2.real=0;
      tempcomplex2.imag=-1;
      //--- allocation
      ArrayResize(temparray,4);
      //--- initialization
      temparray[0]=0;
      temparray[1]=tempcomplex1;
      temparray[2]=0;
      temparray[3]=tempcomplex2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex_Vector(z,temparray,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Integrating f=exp(x) by adaptive integrator                      |
//+------------------------------------------------------------------+
void TEST_AutoGK_D1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CObject              obj;
   CInt_Function_1_Func fint;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- This example demonstrates integration of f=exp(x) on [0,1]:
      //--- * first,autogkstate is initialized
      //--- * then we call integration function
      //--- * and finally we obtain results with autogkresults() call
      double a=0;
      //--- check
      if(_spoil_scenario==0)
         a=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==1)
         a=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==2)
         a=CInfOrNaN::NegativeInfinity();
      double b=1;
      //--- check
      if(_spoil_scenario==3)
         b=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         b=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         b=CInfOrNaN::NegativeInfinity();
      //--- create variables
      CAutoGKStateShell s;
      double v;
      CAutoGKReportShell rep;
      //--- function call
      CAlglib::AutoGKSmooth(a,b,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::AutoGKIntegrate(s,fint,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::AutoGKResults(s,v,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,1.7182,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Interpolation and differentiation using barycentric              |
//| representation                                                   |
//+------------------------------------------------------------------+
void TEST_PolInt_D_CalcDiff(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- Here we demonstrate polynomial interpolation and differentiation
      //--- of y=x^2-x sampled at [0,1,2]. Barycentric representation of polynomial is used.
      ArrayResize(x,3);
      //--- initialization
      x[0]=0;
      x[1]=1;
      x[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      double t=-1;
      //--- check
      if(_spoil_scenario==10)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         t=CInfOrNaN::NegativeInfinity();
      //--- create variables
      double v;
      double dv;
      double d2v;
      CBarycentricInterpolantShell p;
      //--- barycentric model is created
      CAlglib::PolynomialBuild(x,y,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- barycentric interpolation is demonstrated
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      //--- barycentric differentation is demonstrated
      CAlglib::BarycentricDiff1(p,t,v,dv);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(dv,-3.0,0.00005);
      //--- second derivatives with barycentric representation
      CAlglib::BarycentricDiff1(p,t,v,dv);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(dv,-3.0,0.00005);
      //--- function call
      CAlglib::BarycentricDiff2(p,t,v,dv,d2v);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(dv,-3.0,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(d2v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Conversion between power basis and barycentric representation    |
//+------------------------------------------------------------------+
void TEST_PolInt_D_Conv(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double a[];
   double temparray[];
   double a2[];
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
     {
      //--- Here we demonstrate conversion of y=x^2-x
      //--- between power basis and barycentric representation.
      ArrayResize(a,3);
      //--- initialization
      a[0]=0;
      a[1]=-1;
      a[2]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(a,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(a,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(a,CInfOrNaN::NegativeInfinity());
      double t=2;
      //--- check
      if(_spoil_scenario==3)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::NegativeInfinity();
      //--- create a variable
      CBarycentricInterpolantShell p;
      //--- a=[0,-1,+1] is decomposition of y=x^2-x in the power basis:
      //---     y=0 - 1*x + 1*x^2
      //--- We convert it to the barycentric form.
      CAlglib::PolynomialPow2Bar(a,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- now we have barycentric interpolation;we can use it for interpolation
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.005);
      //--- we can also convert back from barycentric representation to power basis
      CAlglib::PolynomialBar2Pow(p,a2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,3);
      //--- initialization
      temparray[0]=0;
      temparray[1]=-1;
      temparray[2]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(a2,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation on special grids (equidistant,          |
//| Chebyshev I/II)                                                  |
//+------------------------------------------------------------------+
void TEST_PolInt_D_Spec(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y_eqdist[];
   double y_cheb1[];
   double y_cheb2[];
   double a_eqdist[];
   double a_cheb1[];
   double a_cheb2[];
   double temparray[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
     {
      //--- Temporaries:
      //--- * values of y=x^2-x sampled at three special grids:
      //---   * equdistant grid spanning [0,2],  x[i]=2*i/(N-1),i=0..N-1
      //---   * Chebyshev-I grid spanning [-1,+1],x[i]=1 + Cos(PI*(2*i+1)/(2*n)),i=0..N-1
      //---   * Chebyshev-II grid spanning [-1,+1],x[i]=1 + Cos(PI*i/(n-1)),i=0..N-1
      //--- * barycentric interpolants for these three grids
      //--- * vectors to store coefficients of quadratic representation
      ArrayResize(y_eqdist,3);
      //--- initialization
      y_eqdist[0]=0;
      y_eqdist[1]=0;
      y_eqdist[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y_eqdist,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y_eqdist,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y_eqdist,CInfOrNaN::NegativeInfinity());
      //--- allocation
      ArrayResize(y_cheb1,3);
      //--- initialization
      y_cheb1[0]=-0.116025;
      y_cheb1[1]=0;
      y_cheb1[2]=1.616025;
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(y_cheb1,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y_cheb1,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y_cheb1,CInfOrNaN::NegativeInfinity());
      //--- allocation
      ArrayResize(y_cheb2,3);
      //--- initialization
      y_cheb2[0]=0;
      y_cheb2[1]=0;
      y_cheb2[2]=2;
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y_cheb2,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y_cheb2,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(y_cheb2,CInfOrNaN::NegativeInfinity());
      //--- create variables
      CBarycentricInterpolantShell p_eqdist;
      CBarycentricInterpolantShell p_cheb1;
      CBarycentricInterpolantShell p_cheb2;
      //--- First,we demonstrate construction of barycentric interpolants on
      //--- special grids. We unpack power representation to ensure that
      //--- interpolant was built correctly.
      //--- In all three cases we should get same quadratic function.
      CAlglib::PolynomialBuildEqDist(0.0,2.0,y_eqdist,p_eqdist);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::PolynomialBar2Pow(p_eqdist,a_eqdist);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,3);
      //--- initialization
      temparray[0]=0;
      temparray[1]=-1;
      temparray[2]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(a_eqdist,temparray,0.00005);
      //--- function call
      CAlglib::PolynomialBuildCheb1(-1,+1,y_cheb1,p_cheb1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::PolynomialBar2Pow(p_cheb1,a_cheb1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,3);
      //--- initialization
      temparray[0]=0;
      temparray[1]=-1;
      temparray[2]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(a_cheb1,temparray,0.00005);
      //--- function call
      CAlglib::PolynomialBuildCheb2(-1,+1,y_cheb2,p_cheb2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::PolynomialBar2Pow(p_cheb2,a_cheb2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,3);
      //--- initialization
      temparray[0]=0;
      temparray[1]=-1;
      temparray[2]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(a_cheb2,temparray,0.00005);
      //--- Now we demonstrate polynomial interpolation withconstruction
      //--- of the barycentricinterpolant structure.
      //--- We calculate interpolant value at x=-2.
      //--- In all three cases we should get same f=6
      double t=-2;
      //--- check
      if(_spoil_scenario==9)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==10)
         t=CInfOrNaN::NegativeInfinity();
      double v;
      //--- function call
      v=CAlglib::PolynomialCalcEqDist(0.0,2.0,y_eqdist,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,6.0,0.00005);
      //--- function call
      v=CAlglib::PolynomialCalcCheb1(-1,+1,y_cheb1,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,6.0,0.00005);
      //--- function call
      v=CAlglib::PolynomialCalcCheb2(-1,+1,y_cheb2,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,6.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation,full list of parameters.                |
//+------------------------------------------------------------------+
void TEST_PolInt_T_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(x,3);
      //--- initialization
      x[0]=0;
      x[1]=1;
      x[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      double t=-1;
      //--- check
      if(_spoil_scenario==8)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==9)
         t=CInfOrNaN::NegativeInfinity();
      //--- create variables
      CBarycentricInterpolantShell p;
      double v;
      //--- function call
      CAlglib::PolynomialBuild(x,y,3,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation,full list of parameters.                |
//+------------------------------------------------------------------+
void TEST_PolInt_T_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(y);
      t=-1;
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         t=CInfOrNaN::NegativeInfinity();
      CBarycentricInterpolantShell p;
      //--- function call
      CAlglib::PolynomialBuildEqDist(0.0,2.0,y,3,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation,full list of parameters.                |
//+------------------------------------------------------------------+
void TEST_PolInt_T_3(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=-0.116025;
      y[1]=0;
      y[2]=1.616025;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(y);
      t=-1;
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         t=CInfOrNaN::NegativeInfinity();
      CBarycentricInterpolantShell p;
      //--- function call
      CAlglib::PolynomialBuildCheb1(-1.0,+1.0,y,3,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- create a variable
      double v;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation,full list of parameters.                |
//+------------------------------------------------------------------+
void TEST_PolInt_T_4(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double a;
   double b;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(y);
      t=-2;
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         t=CInfOrNaN::NegativeInfinity();
      a=-1;
      //--- check
      if(_spoil_scenario==6)
         a=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         a=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         a=CInfOrNaN::NegativeInfinity();
      b=1;
      //--- check
      if(_spoil_scenario==9)
         b=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         b=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         b=CInfOrNaN::NegativeInfinity();
      //--- create a variable
      CBarycentricInterpolantShell p;
      //--- function call
      CAlglib::PolynomialBuildCheb2(a,b,y,3,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,6.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation,full list of parameters.                |
//+------------------------------------------------------------------+
void TEST_PolInt_T_5(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(y);
      t=-1;
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         t=CInfOrNaN::NegativeInfinity();
      //--- function call
      v=CAlglib::PolynomialCalcEqDist(0.0,2.0,y,3,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation,full list of parameters.                |
//+------------------------------------------------------------------+
void TEST_PolInt_T_6(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double a;
   double b;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=-0.116025;
      y[1]=0;
      y[2]=1.616025;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(y);
      t=-1;
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         t=CInfOrNaN::NegativeInfinity();
      a=-1;
      //--- check
      if(_spoil_scenario==6)
         a=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         a=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         a=CInfOrNaN::NegativeInfinity();
      b=1;
      //--- check
      if(_spoil_scenario==9)
         b=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         b=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         b=CInfOrNaN::NegativeInfinity();
      //--- function call
      v=CAlglib::PolynomialCalcCheb1(a,b,y,3,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation,full list of parameters.                |
//+------------------------------------------------------------------+
void TEST_PolInt_T_7(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double a;
   double b;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(y);
      t=-2;
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         t=CInfOrNaN::NegativeInfinity();
      a=-1;
      //--- check
      if(_spoil_scenario==6)
         a=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         a=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         a=CInfOrNaN::NegativeInfinity();
      b=1;
      //--- check
      if(_spoil_scenario==9)
         b=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==10)
         b=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         b=CInfOrNaN::NegativeInfinity();
      //--- function call
      v=CAlglib::PolynomialCalcCheb2(a,b,y,3,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,6.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation: y=x^2-x,equidistant grid, barycentric  |
//| form                                                             |
//+------------------------------------------------------------------+
void TEST_PolInt_T_8(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      t=-1;
      //--- check
      if(_spoil_scenario==3)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::NegativeInfinity();
      CBarycentricInterpolantShell p;
      //--- function call
      CAlglib::PolynomialBuildEqDist(0.0,2.0,y,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind),  |
//| barycentric form                                                 |
//+------------------------------------------------------------------+
void TEST_PolInt_T_9(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double a;
   double b;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=-0.116025;
      y[1]=0;
      y[2]=1.616025;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      t=-1;
      //--- check
      if(_spoil_scenario==3)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::NegativeInfinity();
      a=-1;
      //--- check
      if(_spoil_scenario==5)
         a=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==6)
         a=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==7)
         a=CInfOrNaN::NegativeInfinity();
      b=1;
      //--- check
      if(_spoil_scenario==8)
         b=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==9)
         b=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==10)
         b=CInfOrNaN::NegativeInfinity();
      //--- create a variable
      CBarycentricInterpolantShell p;
      //--- function call
      CAlglib::PolynomialBuildCheb1(a,b,y,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind), |
//| barycentric form                                                 |
//+------------------------------------------------------------------+
void TEST_PolInt_T_10(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double a;
   double b;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      t=-2;
      //--- check
      if(_spoil_scenario==3)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::NegativeInfinity();
      a=-1;
      //--- check
      if(_spoil_scenario==5)
         a=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==6)
         a=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==7)
         a=CInfOrNaN::NegativeInfinity();
      b=1;
      //--- check
      if(_spoil_scenario==8)
         b=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==9)
         b=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==10)
         b=CInfOrNaN::NegativeInfinity();
      //--- create a variable
      CBarycentricInterpolantShell p;
      //--- function call
      CAlglib::PolynomialBuildCheb2(a,b,y,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,6.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation: y=x^2-x,equidistant grid               |
//+------------------------------------------------------------------+
void TEST_PolInt_T_11(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      t=-1;
      //--- check
      if(_spoil_scenario==3)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::NegativeInfinity();
      //--- function call
      v=CAlglib::PolynomialCalcEqDist(0.0,2.0,y,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation: y=x^2-x,Chebyshev grid (first kind)    |
//+------------------------------------------------------------------+
void TEST_PolInt_T_12(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
   double t;
   double a;
   double b;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=-0.116025;
      y[1]=0;
      y[2]=1.616025;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      t=-1;
      //--- check
      if(_spoil_scenario==3)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::NegativeInfinity();
      a=-1;
      //--- check
      if(_spoil_scenario==5)
         a=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==6)
         a=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==7)
         a=CInfOrNaN::NegativeInfinity();
      b=+1;
      //--- check
      if(_spoil_scenario==8)
         b=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==9)
         b=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==10)
         b=CInfOrNaN::NegativeInfinity();
      //--- function call
      v=CAlglib::PolynomialCalcCheb1(a,b,y,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial interpolation: y=x^2-x,Chebyshev grid (second kind)   |
//+------------------------------------------------------------------+
void TEST_PolInt_T_13(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double y[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(y,3);
      //--- initialization
      y[0]=0;
      y[1]=0;
      y[2]=2;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      double t=-2;
      //--- check
      if(_spoil_scenario==3)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==4)
         t=CInfOrNaN::NegativeInfinity();
      double a=-1;
      //--- check
      if(_spoil_scenario==5)
         a=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==6)
         a=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==7)
         a=CInfOrNaN::NegativeInfinity();
      double b=+1;
      //--- check
      if(_spoil_scenario==8)
         b=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==9)
         b=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==10)
         b=CInfOrNaN::NegativeInfinity();
      double v;
      //--- function call
      v=CAlglib::PolynomialCalcCheb2(a,b,y,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,6.0,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Piecewise linear spline interpolation                            |
//+------------------------------------------------------------------+
void TEST_Spline1D_D_Linear(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   double t;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- We use piecewise linear spline to interpolate f(x)=x^2 sampled
      //--- at 5 equidistant nodes on [-1,+1].
      ArrayResize(x,5);
      //--- initialization
      x[0]=-1;
      x[1]=-0.5;
      x[2]=0;
      x[3]=0.5;
      x[4]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,5);
      //--- initialization
      y[0]=1;
      y[1]=0.25;
      y[2]=0;
      y[3]=0.25;
      y[4]=1;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      t=0.25;
      //--- check
      if(_spoil_scenario==10)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         t=CInfOrNaN::NegativeInfinity();
      //--- create a variable
      CSpline1DInterpolantShell s;
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- build spline
      CAlglib::Spline1DBuildLinear(x,y,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- calculate S(0.25) - it is quite different from 0.25^2=0.0625
      v=CAlglib::Spline1DCalc(s,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,0.125,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Cubic spline interpolation                                       |
//+------------------------------------------------------------------+
void TEST_Spline1D_D_Cubic(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
     {
      //--- We use cubic spline to interpolate f(x)=x^2 sampled
      //--- at 5 equidistant nodes on [-1,+1].
      //--- First,we use default boundary conditions ("parabolically terminated
      //--- spline") because cubic spline built with such boundary conditions
      //--- will exactly reproduce any quadratic f(x).
      //--- Then we try to use natural boundary conditions
      //---     d2S(-1)/dx^2=0.0
      //---     d2S(+1)/dx^2=0.0
      //--- and see that such spline interpolated f(x) with small error.
      ArrayResize(x,5);
      //--- initialization
      x[0]=-1;
      x[1]=-0.5;
      x[2]=0;
      x[3]=0.5;
      x[4]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,5);
      //--- initialization
      y[0]=1;
      y[1]=0.25;
      y[2]=0;
      y[3]=0.25;
      y[4]=1;
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      double t=0.25;
      //--- check
      if(_spoil_scenario==8)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==9)
         t=CInfOrNaN::NegativeInfinity();
      //--- create variables
      double v;
      CSpline1DInterpolantShell s;
      int natural_bound_type=2;
      //--- Test exact boundary conditions: build S(x),calculare S(0.25)
      //--- (almost same as original function)
      CAlglib::Spline1DBuildCubic(x,y,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::Spline1DCalc(s,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,0.0625,0.00001);
      //--- Test natural boundary conditions: build S(x),calculare S(0.25)
      //--- (small interpolation error)
      CAlglib::Spline1DBuildCubic(x,y,5,natural_bound_type,0.0,natural_bound_type,0.0,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::Spline1DCalc(s,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,0.0580,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Differentiation on the grid using cubic splines                  |
//+------------------------------------------------------------------+
void TEST_Spline1D_D_GridDiff(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   double d1[];
   double d2[];
   double temparray1[];
   double temparray2[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
     {
      //--- We use cubic spline to do grid differentiation,i.e. having
      //--- values of f(x)=x^2 sampled at 5 equidistant nodes on [-1,+1]
      //--- we calculate derivatives of cubic spline at nodes WITHOUT
      //--- CONSTRUCTION OF SPLINE OBJECT.
      //--- There are efficient functions spline1dgriddiffcubic() and
      //--- spline1dgriddiff2cubic() for such calculations.
      //--- We use default boundary conditions ("parabolically terminated
      //--- spline") because cubic spline built with such boundary conditions
      //--- will exactly reproduce any quadratic f(x).
      //--- Actually,we could use natural conditions,but we feel that
      //--- spline which exactly reproduces f() will show us more
      //--- understandable results.
      ArrayResize(x,5);
      //--- initialization
      x[0]=-1;
      x[1]=-0.5;
      x[2]=0;
      x[3]=0.5;
      x[4]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,5);
      //--- initialization
      y[0]=1;
      y[1]=0.25;
      y[2]=0;
      y[3]=0.25;
      y[4]=1;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- We calculate first derivatives: they must be equal to 2*x
      CAlglib::Spline1DGridDiffCubic(x,y,d1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray1,5);
      //--- initialization
      temparray1[0]=-2;
      temparray1[1]=-1;
      temparray1[2]=0;
      temparray1[3]=1;
      temparray1[4]=2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(d1,temparray1,0.0001);
      //--- Now test griddiff2,which returns first AND second derivatives.
      //--- First derivative is 2*x,second is equal to 2.0
      CAlglib::Spline1DGridDiff2Cubic(x,y,d1,d2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray2,5);
      //--- initialization
      temparray2[0]=2;
      temparray2[1]=2;
      temparray2[2]=2;
      temparray2[3]=2;
      temparray2[4]=2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(d1,temparray1,0.0001);
      _TestResult=_TestResult && Doc_Test_Real_Vector(d2,temparray2,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Resampling using cubic splines                                   |
//+------------------------------------------------------------------+
void TEST_Spline1D_D_ConvDiff(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x_old[];
   double y_old[];
   double x_new[];
   double y_new[];
   double d1_new[];
   double d2_new[];
   double temparray1[];
   double temparray2[];
   double temparray3[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
     {
      //--- We use cubic spline to do resampling,i.e. having
      //--- values of f(x)=x^2 sampled at 5 equidistant nodes on [-1,+1]
      //--- we calculate values/derivatives of cubic spline on
      //--- another grid (equidistant with 9 nodes on [-1,+1])
      //--- WITHOUT CONSTRUCTION OF SPLINE OBJECT.
      //--- There are efficient functions spline1dconvcubic(),
      //--- spline1dconvdiffcubic() and spline1dconvdiff2cubic()
      //--- for such calculations.
      //--- We use default boundary conditions ("parabolically terminated
      //--- spline") because cubic spline built with such boundary conditions
      //--- will exactly reproduce any quadratic f(x).
      //--- Actually,we could use natural conditions,but we feel that
      //--- spline which exactly reproduces f() will show us more
      //--- understandable results.
      ArrayResize(x_old,5);
      //--- initialization
      x_old[0]=-1;
      x_old[1]=-0.5;
      x_old[2]=0;
      x_old[3]=0.5;
      x_old[4]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x_old,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x_old,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x_old,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x_old);
      //--- allocation
      ArrayResize(y_old,5);
      //--- initialization
      y_old[0]=1;
      y_old[1]=0.25;
      y_old[2]=0;
      y_old[3]=0.25;
      y_old[4]=1;
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y_old,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y_old,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y_old,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y_old);
      //--- allocation
      ArrayResize(x_new,9);
      //--- initialization
      x_new[0]=-1;
      x_new[1]=-0.75;
      x_new[2]=-0.5;
      x_new[3]=-0.25;
      x_new[4]=0;
      x_new[5]=0.25;
      x_new[6]=0.5;
      x_new[7]=0.75;
      x_new[8]=1;
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(x_new,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(x_new,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(x_new,CInfOrNaN::NegativeInfinity());
      //--- First,conversion withdifferentiation.
      CAlglib::Spline1DConvCubic(x_old,y_old,x_new,y_new);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray1,9);
      //--- initialization
      temparray1[0]=1;
      temparray1[1]=0.5625;
      temparray1[2]=0.25;
      temparray1[3]=0.0625;
      temparray1[4]=0;
      temparray1[5]=0.0625;
      temparray1[6]=0.25;
      temparray1[7]=0.5625;
      temparray1[8]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(y_new,temparray1,0.0001);
      //--- Then,conversion with differentiation (first derivatives only)
      CAlglib::Spline1DConvDiffCubic(x_old,y_old,x_new,y_new,d1_new);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray2,9);
      //--- initialization
      temparray2[0]=-2;
      temparray2[1]=-1.5;
      temparray2[2]=-1;
      temparray2[3]=-0.5;
      temparray2[4]=0;
      temparray2[5]=0.5;
      temparray2[6]=1.0;
      temparray2[7]=1.5;
      temparray2[8]=2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(y_new,temparray1,0.0001);
      _TestResult=_TestResult && Doc_Test_Real_Vector(d1_new,temparray2,0.0001);
      //--- Finally,conversion with first and second derivatives
      CAlglib::Spline1DConvDiff2Cubic(x_old,y_old,x_new,y_new,d1_new,d2_new);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray3,9);
      //--- initialization
      temparray3[0]=2;
      temparray3[1]=2;
      temparray3[2]=2;
      temparray3[3]=2;
      temparray3[4]=2;
      temparray3[5]=2;
      temparray3[6]=2;
      temparray3[7]=2;
      temparray3[8]=2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(y_new,temparray1,0.0001);
      _TestResult=_TestResult && Doc_Test_Real_Vector(d1_new,temparray2,0.0001);
      _TestResult=_TestResult && Doc_Test_Real_Vector(d2_new,temparray3,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Unconstrained dense quadratic programming                        |
//+------------------------------------------------------------------+
void TEST_MinQP_D_U1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
   double b[];
   double x0[];
   double x[];
   double temparray[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1)=x0^2 + x1^2 -6*x0 - 4*x1
      //--- Exact solution is [x0,x1]=[3,2]
      //--- We provide algorithm with starting point,although in this case
      //--- (dense matrix,no constraints) it can work withsuch information.
      //--- IMPORTANT: this solver minimizes  following  function:
      //---     f(x)=0.5*x'*A*x + b'*x.
      //--- Note that quadratic term has 0.5 before it. So if you want to minimize
      //--- quadratic function,you should rewrite it in such way that quadratic term
      //--- is multiplied by 0.5 too.
      //--- For example,our function is f(x)=x0^2+x1^2+...,but we rewrite it as
      //---     f(x)=0.5*(2*x0^2+2*x1^2) + ....
      //--- and pass diag(2,2) as quadratic term - NOT diag(1,1)!
      a.Resize(2,2);
      //--- initialization
      a.Set(0,0,2);
      a.Set(0,1,0);
      a.Set(1,0,0);
      a.Set(1,1,2);
      //--- check
      if(_spoil_scenario>=0 && _spoil_scenario<=2)
         do
           {
            if(_spoil_scenario==0)
               Spoil_Matrix_By_Value(a,CInfOrNaN::NaN());
            //--- check
            if(_spoil_scenario==1)
               Spoil_Matrix_By_Value(a,CInfOrNaN::PositiveInfinity());
            //--- check
            if(_spoil_scenario==2)
               Spoil_Matrix_By_Value(a,CInfOrNaN::NegativeInfinity());
           }
         while(CApServ::IsFiniteRTrMatrix(a,2,false));
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(a);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(a);
      //--- allocation
      ArrayResize(b,2);
      //--- initialization
      b[0]=-6;
      b[1]=-4;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(b,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(b,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(b,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(b);
      //--- allocation
      ArrayResize(x0,2);
      //--- initialization
      x0[0]=0;
      x0[1]=1;
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(x0,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(x0,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(x0,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Deleting_Element(x0);
      //--- create variables
      CMinQPStateShell state;
      CMinQPReportShell rep;
      //--- function call
      CAlglib::MinQPCreate(2,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPSetQuadraticTerm(state,a);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPSetLinearTerm(state,b);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPSetStartingPoint(state,x0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=3;
      temparray[1]=2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Constrained dense quadratic programming                          |
//+------------------------------------------------------------------+
void TEST_MinQP_D_BC1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble a;
   double b[];
   double x0[];
   double bndl[];
   double bndu[];
   double temparray[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<17; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1)=x0^2 + x1^2 -6*x0 - 4*x1
      //--- subject to bound constraints 0<=x0<=2.5,0<=x1<=2.5
      //--- Exact solution is [x0,x1]=[2.5,2]
      //--- We provide algorithm with starting point. With such small problem good starting
      //--- point is not really necessary,but with high-dimensional problem it can save us
      //--- a lot of time.
      //--- IMPORTANT: this solver minimizes  following  function:
      //---     f(x)=0.5*x'*A*x + b'*x.
      //--- Note that quadratic term has 0.5 before it. So if you want to minimize
      //--- quadratic function,you should rewrite it in such way that quadratic term
      //--- is multiplied by 0.5 too.
      //--- For example,our function is f(x)=x0^2+x1^2+...,but we rewrite it as
      //---     f(x)=0.5*(2*x0^2+2*x1^2) + ....
      //--- and pass diag(2,2) as quadratic term - NOT diag(1,1)!
      a.Resize(2,2);
      //--- initialization
      a.Set(0,0,2);
      a.Set(0,1,0);
      a.Set(1,0,0);
      a.Set(1,1,2);
      //--- check
      //--- check
      if(_spoil_scenario>=0 && _spoil_scenario<=2)
         do
           {
            if(_spoil_scenario==0)
               Spoil_Matrix_By_Value(a,CInfOrNaN::NaN());
            //--- check
            if(_spoil_scenario==1)
               Spoil_Matrix_By_Value(a,CInfOrNaN::PositiveInfinity());
            //--- check
            if(_spoil_scenario==2)
               Spoil_Matrix_By_Value(a,CInfOrNaN::NegativeInfinity());
           }
         while(CApServ::IsFiniteRTrMatrix(a,2,false));
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(a);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(a);
      //--- allocation
      ArrayResize(b,2);
      //--- initialization
      b[0]=-6;
      b[1]=-4;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(b,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(b,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(b,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(b);
      //--- allocation
      ArrayResize(x0,2);
      //--- initialization
      x0[0]=0;
      x0[1]=1;
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(x0,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(x0,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(x0,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Deleting_Element(x0);
      //--- allocation
      ArrayResize(bndl,2);
      //--- initialization
      bndl[0]=0;
      bndl[1]=0;
      //--- check
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(bndl,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==14)
         Spoil_Vector_By_Deleting_Element(bndl);
      //--- allocation
      ArrayResize(bndu,2);
      //--- initialization
      bndu[0]=2.5;
      bndu[1]=2.5;
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(bndu,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Deleting_Element(bndu);
      double x[];
      CMinQPStateShell state;
      CMinQPReportShell rep;
      //--- function call
      CAlglib::MinQPCreate(2,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPSetQuadraticTerm(state,a);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPSetLinearTerm(state,b);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPSetStartingPoint(state,x0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=2.5;
      temparray[1]=2;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear least squares optimization using function vector only  |
//+------------------------------------------------------------------+
void TEST_MinLM_D_V(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              s[];
   double              temparray[];
   CObject             obj;
   CNDimensional_FVec1 ffvec;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1)=f0^2+f1^2,where
      //---     f0(x0,x1)=10*(x0+3)^2
      //---     f1(x0,x1)=(x1-3)^2
      //--- using "V" mode of the Levenberg-Marquardt optimizer.
      //--- Optimization algorithm uses:
      //--- * function vector f[]={f1,f2}
      //--- No other information (Jacobian,gradient,etc.) is needed.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      ArrayResize(s,2);
      ArrayInitialize(s,1);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,CInfOrNaN::NegativeInfinity());
      double epsx=0.0000000001;
      //--- check
      if(_spoil_scenario==6)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsx=CInfOrNaN::NegativeInfinity();
      int maxits=0;
      //--- objects of classes
      CMinLMStateShell  state;
      CMinLMReportShell rep;
      //--- function call
      CAlglib::MinLMCreateV(2,x,0.0001,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMOptimize(state,ffvec,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear least squares optimization using function vector and   |
//| Jacobian                                                         |
//+------------------------------------------------------------------+
void TEST_MinLM_D_VJ(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
//--- TEST minlm_d_vj
//---      Nonlinear least squares optimization using function vector and Jacobian
   _TestResult=true;
//--- create variables
   double              x[];
   double              s[];
   double              temparray[];
   CObject             obj;
   CNDimensional_FVec1 ffvec;
   CNDimensional_Jac1  fjac;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1)=f0^2+f1^2,where
      //---     f0(x0,x1)=10*(x0+3)^2
      //---     f1(x0,x1)=(x1-3)^2
      //--- using "VJ" mode of the Levenberg-Marquardt optimizer.
      //--- Optimization algorithm uses:
      //--- * function vector f[]={f1,f2}
      //--- * Jacobian matrix J={dfi/dxj}.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      ArrayResize(s,2);
      ArrayInitialize(s,1);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,CInfOrNaN::NegativeInfinity());
      double epsx=0.0000000001;
      //--- check
      if(_spoil_scenario==6)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinLMStateShell state;
      CMinLMReportShell rep;
      //--- function call
      CAlglib::MinLMCreateVJ(2,x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMOptimize(state,ffvec,fjac,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      //_TestResult = _TestResult  &&  Doc_Test_Int(rep.GetTerminationType(), 4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear Hessian-based optimization for general functions       |
//+------------------------------------------------------------------+
void TEST_MinLM_D_FGH(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              s[];
   double              temparray[];
   CObject             obj;
   CNDimensional_Func1 ffunc;
   CNDimensional_Grad1 fgrad;
   CNDimensional_Hess1 fhess;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1)=100*(x0+3)^4+(x1-3)^4
      //--- using "FGH" mode of the Levenberg-Marquardt optimizer.
      //--- F is treated like a monolitic function withinternal structure,
      //--- i.e. we do NOT represent it as a sum of squares.
      //--- Optimization algorithm uses:
      //--- * function value F(x0,x1)
      //--- * gradient G={dF/dxi}
      //--- * Hessian H={d2F/(dxi*dxj)}
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      ArrayResize(s,2);
      ArrayInitialize(s,1);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,CInfOrNaN::NegativeInfinity());
      double epsx=0.0000000001;
      //--- check
      if(_spoil_scenario==6)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinLMStateShell state;
      CMinLMReportShell rep;
      //--- function call
      CAlglib::MinLMCreateFGH(x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMOptimize(state,ffunc,fgrad,fhess,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      //_TestResult = _TestResult  &&  Doc_Test_Int(rep.GetTerminationType(), 4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Bound constrained nonlinear least squares optimization           |
//+------------------------------------------------------------------+
void TEST_MinLM_D_VB(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              s[];
   double              bndl[];
   double              bndu[];
   double              temparray[];
   CObject             obj;
   CNDimensional_FVec1 ffvec;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1)=f0^2+f1^2,where
      //---     f0(x0,x1)=10*(x0+3)^2
      //---     f1(x0,x1)=(x1-3)^2
      //--- with boundary constraints
      //---     -1 <= x0 <= +1
      //---     -1 <= x1 <= +1
      //--- using "V" mode of the Levenberg-Marquardt optimizer.
      //--- Optimization algorithm uses:
      //--- * function vector f[]={f1,f2}
      //--- No other information (Jacobian,gradient,etc.) is needed.
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- allocation
      ArrayResize(bndl,2);
      //--- initialization
      bndl[0]=-1;
      bndl[1]=-1;
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(bndl,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(bndl);
      //--- allocation
      ArrayResize(bndu,2);
      //--- initialization
      bndu[0]=1;
      bndu[1]=1;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(bndu,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Deleting_Element(bndu);
      ArrayResize(s,2);
      ArrayInitialize(s,1);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(s,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(s,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(s,CInfOrNaN::NegativeInfinity());
      double epsx=0.0000000001;
      //--- check
      if(_spoil_scenario==10)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==11)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==12)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinLMStateShell state;
      CMinLMReportShell rep;
      //--- function call
      CAlglib::MinLMCreateV(2,x,0.0001,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMOptimize(state,ffvec,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-1;
      temparray[1]=1;
      //--- check result
      //_TestResult = _TestResult  &&  Doc_Test_Int(rep.GetTerminationType(), 4);
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Efficient restarts of LM optimizer                               |
//+------------------------------------------------------------------+
void TEST_MinLM_D_Restarts(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              s[];
   double              temparray[];
   CObject             obj;
   CNDimensional_FVec1 ffvec1;
   CNDimensional_FVec2 ffvec2;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1)=f0^2+f1^2,where
      //---     f0(x0,x1)=10*(x0+3)^2
      //---     f1(x0,x1)=(x1-3)^2
      //--- using several starting points and efficient restarts.
      ArrayResize(s,2);
      ArrayInitialize(s,1);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(s,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(s,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(s,CInfOrNaN::NegativeInfinity());
      double epsx=0.0000000001;
      //--- check
      if(_spoil_scenario==3)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==4)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==5)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinLMStateShell state;
      CMinLMReportShell rep;
      //--- create optimizer using minlmcreatev()
      ArrayResize(x,2);
      //--- initialization
      x[0]=10;
      x[1]=10;
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::MinLMCreateV(2,x,0.0001,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMOptimize(state,ffvec1,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      //--- restart optimizer using minlmrestartfrom()
      //--- we can use different starting point,different function,
      //--- different stopping conditions,but problem size
      //--- must remain unchanged.
      ArrayResize(x,2);
      //--- initialization
      x[0]=4;
      x[1]=4;
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==14)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- function call
      CAlglib::MinLMRestartFrom(state,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMOptimize(state,ffvec2,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=0;
      temparray[1]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear least squares optimization, FJ scheme (obsolete, but   |
//| supported)                                                       |
//+------------------------------------------------------------------+
void TEST_MinLM_T_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double              x[];
   double              s[];
   double              temparray[];
   CObject             obj;
   CNDimensional_Func1 ffunc;
   CNDimensional_Jac1  fjac;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      double epsg=0.0000000001;
      ArrayResize(s,2);
      ArrayInitialize(s,1);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,CInfOrNaN::NegativeInfinity());
      double epsx=0.0000000001;
      //--- check
      if(_spoil_scenario==6)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsx=CInfOrNaN::NegativeInfinity();
      int maxits=0;
      CMinLMStateShell state;
      CMinLMReportShell rep;
      //--- function call
      CAlglib::MinLMCreateFJ(2,x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMOptimize(state,ffunc,fjac,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear least squares optimization, FGJ scheme (obsolete, but  |
//| supported)                                                       |
//+------------------------------------------------------------------+
void TEST_MinLM_T_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double  x[];
   double  s[];
   double  temparray[];
   CObject obj;
   CNDimensional_Func1 ffunc;
   CNDimensional_Grad1 fgrad;
   CNDimensional_Jac1  fjac;
   CNDimensional_Rep   frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(x,2);
      //--- initialization
      x[0]=0;
      x[1]=0;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      ArrayResize(s,2);
      ArrayInitialize(s,1);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,CInfOrNaN::NegativeInfinity());
      double epsx=0.0000000001;
      //--- check
      if(_spoil_scenario==6)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==7)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==8)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      CMinLMStateShell state;
      CMinLMReportShell rep;
      //--- function call
      CAlglib::MinLMCreateFGJ(2,x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMOptimize(state,ffunc,fgrad,fjac,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::MinLMResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=-3;
      temparray[1]=3;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,temparray,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear fitting using function value only                      |
//+------------------------------------------------------------------+
void TEST_LSFit_D_NLF(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble           x;
   double                  y[];
   double                  c[];
   double                  temparray[];
   double                  w[];
   CObject                 obj;
   CNDimensional_CX_1_Func fcx1func;
   CNDimensional_Rep       frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<19; _spoil_scenario++)
     {
      //--- In this example we demonstrate exponential fitting
      //--- by f(x)=exp(-c*x^2)
      //--- using function value only.
      //--- Gradient is estimated using combination of numerical differences
      //--- and secant updates. diffstep variable stores differentiation step
      //--- (we have to tell algorithm what step to use).
      x.Resize(11,1);
      //--- initialization
      x.Set(0,0,-1);
      x.Set(1,0,-0.8);
      x.Set(2,0,-0.6);
      x.Set(3,0,-0.4);
      x.Set(4,0,-0.2);
      x.Set(5,0,0);
      x.Set(6,0,0.2);
      x.Set(7,0,0.4);
      x.Set(8,0,0.6);
      x.Set(9,0,0.8);
      x.Set(10,0,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(x);
      //--- allocation
      ArrayResize(y,11);
      //--- initialization
      y[0]=0.22313;
      y[1]=0.382893;
      y[2]=0.582748;
      y[3]=0.786628;
      y[4]=0.941765;
      y[5]=1;
      y[6]=0.941765;
      y[7]=0.786628;
      y[8]=0.582748;
      y[9]=0.382893;
      y[10]=0.22313;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(c,1);
      //--- initialization
      c[0]=0.3;
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(c,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(c,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(c,CInfOrNaN::NegativeInfinity());
      double epsx=0.000001;
      //--- check
      if(_spoil_scenario==13)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==14)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==15)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      int info;
      CLSFitStateShell state;
      CLSFitReportShell rep;
      double diffstep=0.0001;
      //--- check
      if(_spoil_scenario==16)
         diffstep=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==17)
         diffstep=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==18)
         diffstep=CInfOrNaN::NegativeInfinity();
      //--- Fitting withweights
      CAlglib::LSFitCreateF(x,y,c,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitFit(state,fcx1func,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitResults(state,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1.5;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,2);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.05);
      //--- Fitting with weights
      //--- (you can change weights and see how it changes result)
      ArrayResize(w,11);
      //--- initialization
      w[0]=1;
      w[1]=1;
      w[2]=1;
      w[3]=1;
      w[4]=1;
      w[5]=1;
      w[6]=1;
      w[7]=1;
      w[8]=1;
      w[9]=1;
      w[10]=1;
      //--- check
      if(_spoil_scenario==22)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==23)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==24)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==25)
         Spoil_Vector_By_Adding_Element(w);
      //--- check
      if(_spoil_scenario==26)
         Spoil_Vector_By_Deleting_Element(w);
      //--- function call
      CAlglib::LSFitCreateWF(x,y,w,c,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitFit(state,fcx1func,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitResults(state,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1.5;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,2);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear fitting using gradient                                 |
//+------------------------------------------------------------------+
void TEST_LSFit_D_NLFG(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble           x;
   double                  y[];
   double                  c[];
   double                  temparray[];
   double                  w[];
   CObject                 obj;
   CNDimensional_CX_1_Func fcx1func;
   CNDimensional_CX_1_Grad fcx1grad;
   CNDimensional_Rep       frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
     {
      //--- In this example we demonstrate exponential fitting
      //--- by f(x)=exp(-c*x^2)
      //--- using function value and gradient (with respect to c).
      x.Resize(11,1);
      //--- initialization
      x.Set(0,0,-1);
      x.Set(1,0,-0.8);
      x.Set(2,0,-0.6);
      x.Set(3,0,-0.4);
      x.Set(4,0,-0.2);
      x.Set(5,0,0);
      x.Set(6,0,0.2);
      x.Set(7,0,0.4);
      x.Set(8,0,0.6);
      x.Set(9,0,0.8);
      x.Set(10,0,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(x);
      //--- allocation
      ArrayResize(y,11);
      //--- initialization
      y[0]=0.22313;
      y[1]=0.382893;
      y[2]=0.582748;
      y[3]=0.786628;
      y[4]=0.941765;
      y[5]=1;
      y[6]=0.941765;
      y[7]=0.786628;
      y[8]=0.582748;
      y[9]=0.382893;
      y[10]=0.22313;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(c,1);
      //--- initialization
      c[0]=0.3;
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(c,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(c,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(c,CInfOrNaN::NegativeInfinity());
      double epsx=0.000001;
      //--- check
      if(_spoil_scenario==13)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==14)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==15)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      int info;
      CLSFitStateShell state;
      CLSFitReportShell rep;
      //--- Fitting withweights
      CAlglib::LSFitCreateFG(x,y,c,true,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitFit(state,fcx1func,fcx1grad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitResults(state,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1.5;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,2);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.05);
      //--- Fitting with weights
      //--- (you can change weights and see how it changes result)
      ArrayResize(w,11);
      //--- initialization
      w[0]=1;
      w[1]=1;
      w[2]=1;
      w[3]=1;
      w[4]=1;
      w[5]=1;
      w[6]=1;
      w[7]=1;
      w[8]=1;
      w[9]=1;
      w[10]=1;
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Adding_Element(w);
      //--- check
      if(_spoil_scenario==20)
         Spoil_Vector_By_Deleting_Element(w);
      //--- function call
      CAlglib::LSFitCreateWFG(x,y,w,c,true,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitFit(state,fcx1func,fcx1grad,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitResults(state,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1.5;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,2);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear fitting using gradient and Hessian                     |
//+------------------------------------------------------------------+
void TEST_LSFit_D_NLFGH(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble           x;
   double                  y[];
   double                  c[];
   double                  temparray[];
   double                  w[];
   CObject                 obj;
   CNDimensional_CX_1_Func fcx1func;
   CNDimensional_CX_1_Grad fcx1grad;
   CNDimensional_CX_1_Hess fcx1hess;
   CNDimensional_Rep       frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
     {
      //--- In this example we demonstrate exponential fitting
      //--- by f(x)=exp(-c*x^2)
      //--- using function value,gradient and Hessian (with respect to c)
      x.Resize(11,1);
      //--- initialization
      x.Set(0,0,-1);
      x.Set(1,0,-0.8);
      x.Set(2,0,-0.6);
      x.Set(3,0,-0.4);
      x.Set(4,0,-0.2);
      x.Set(5,0,0);
      x.Set(6,0,0.2);
      x.Set(7,0,0.4);
      x.Set(8,0,0.6);
      x.Set(9,0,0.8);
      x.Set(10,0,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(x);
      //--- allocation
      ArrayResize(y,11);
      //--- initialization
      y[0]=0.22313;
      y[1]=0.382893;
      y[2]=0.582748;
      y[3]=0.786628;
      y[4]=0.941765;
      y[5]=1;
      y[6]=0.941765;
      y[7]=0.786628;
      y[8]=0.582748;
      y[9]=0.382893;
      y[10]=0.22313;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(c,1);
      //--- initialization
      c[0]=0.3;
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(c,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(c,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(c,CInfOrNaN::NegativeInfinity());
      double epsx=0.000001;
      //--- check
      if(_spoil_scenario==13)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==14)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==15)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      int info;
      CLSFitStateShell state;
      CLSFitReportShell rep;
      //--- Fitting withweights
      CAlglib::LSFitCreateFGH(x,y,c,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitFit(state,fcx1func,fcx1grad,fcx1hess,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitResults(state,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1.5;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,2);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.05);
      //--- Fitting with weights
      //--- (you can change weights and see how it changes result)
      ArrayResize(w,11);
      //--- initialization
      w[0]=1;
      w[1]=1;
      w[2]=1;
      w[3]=1;
      w[4]=1;
      w[5]=1;
      w[6]=1;
      w[7]=1;
      w[8]=1;
      w[9]=1;
      w[10]=1;
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Adding_Element(w);
      //--- check
      if(_spoil_scenario==20)
         Spoil_Vector_By_Deleting_Element(w);
      //--- function call
      CAlglib::LSFitCreateWFGH(x,y,w,c,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitFit(state,fcx1func,fcx1grad,fcx1hess,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitResults(state,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1.5;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,2);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Bound contstrained nonlinear fitting using function value only   |
//+------------------------------------------------------------------+
void TEST_LSFit_D_NLFB(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble           x;
   double                  y[];
   double                  c[];
   double                  temparray[];
   double                  w[];
   double                  bndl[];
   double                  bndu[];
   CObject                 obj;
   CNDimensional_CX_1_Func fcx1func;
   CNDimensional_Rep       frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<23; _spoil_scenario++)
     {
      //--- In this example we demonstrate exponential fitting by
      //---     f(x)=exp(-c*x^2)
      //--- subject to bound constraints
      //---     0.0 <= c <= 1.0
      //--- using function value only.
      //--- Gradient is estimated using combination of numerical differences
      //--- and secant updates. diffstep variable stores differentiation step
      //--- (we have to tell algorithm what step to use).
      //--- Unconstrained solution is c=1.5,but because of constraints we should
      //--- get c=1.0 (at the boundary).
      x.Resize(11,1);
      //--- initialization
      x.Set(0,0,-1);
      x.Set(1,0,-0.8);
      x.Set(2,0,-0.6);
      x.Set(3,0,-0.4);
      x.Set(4,0,-0.2);
      x.Set(5,0,0);
      x.Set(6,0,0.2);
      x.Set(7,0,0.4);
      x.Set(8,0,0.6);
      x.Set(9,0,0.8);
      x.Set(10,0,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(x);
      //--- allocation
      ArrayResize(y,11);
      //--- initialization
      y[0]=0.22313;
      y[1]=0.382893;
      y[2]=0.582748;
      y[3]=0.786628;
      y[4]=0.941765;
      y[5]=1;
      y[6]=0.941765;
      y[7]=0.786628;
      y[8]=0.582748;
      y[9]=0.382893;
      y[10]=0.22313;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(c,1);
      //--- initialization
      c[0]=0.3;
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(c,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(c,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(c,CInfOrNaN::NegativeInfinity());
      //--- allocation
      ArrayResize(bndl,1);
      //--- initialization
      bndl[0]=0;
      //--- check
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(bndl,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==14)
         Spoil_Vector_By_Deleting_Element(bndl);
      //--- allocation
      ArrayResize(bndu,1);
      //--- initialization
      bndu[0]=1;
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(bndu,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Deleting_Element(bndu);
      double epsx=0.000001;
      //--- check
      if(_spoil_scenario==17)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==18)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==19)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int maxits=0;
      int info;
      CLSFitStateShell state;
      CLSFitReportShell rep;
      double diffstep=0.0001;
      //--- check
      if(_spoil_scenario==20)
         diffstep=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==21)
         diffstep=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==22)
         diffstep=CInfOrNaN::NegativeInfinity();
      //--- function call
      CAlglib::LSFitCreateF(x,y,c,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitFit(state,fcx1func,frep,0,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitResults(state,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear fitting with custom scaling and bound constraints      |
//+------------------------------------------------------------------+
void TEST_LSFit_D_NLScale(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble           x;
   double                  y[];
   double                  c[];
   double                  bndl[];
   double                  bndu[];
   double                  s[];
   double                  temparray[];
   CObject                 obj;
   CNDimensional_Debt_Func fdebtfunc;
   CNDimensional_Rep       frep;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<27; _spoil_scenario++)
     {
      //--- In this example we demonstrate fitting by
      //---     f(x)=c[0]*(1+c[1]*((x-1999)^c[2]-1))
      //--- subject to bound constraints
      //---     -INF  < c[0] < +INF
      //---      -10 <= c[1] <= +10
      //---      0.1 <= c[2] <= 2.0
      //--- Data we want to fit are time series of Japan national debt
      //--- collected from 2000 to 2008 measured in USD (dollars,not
      //--- millions of dollars).
      //--- Our variables are:
      //---     c[0] - debt value at initial moment (2000),
      //---     c[1] - direction coefficient (growth or decrease),
      //---     c[2] - curvature coefficient.
      //--- You may see that our variables are badly scaled - first one
      //--- is order of 10^12,and next two are somewhere ab1 in
      //--- magnitude. Such problem is difficult to solve withsome
      //--- kind of scaling.
      //--- That is exactly where lsfitsetscale() function can be used.
      //--- We set scale of our variables to [1.0E12,1,1],which allows
      //--- us to easily solve this problem.
      //--- You can try commenting lsfitsetscale() call - and you will
      //--- see that algorithm will fail to converge.
      x.Resize(9,1);
      //--- initialization
      x.Set(0,0,2000);
      x.Set(1,0,2001);
      x.Set(2,0,2002);
      x.Set(3,0,2003);
      x.Set(4,0,2004);
      x.Set(5,0,2005);
      x.Set(6,0,2006);
      x.Set(7,0,2007);
      x.Set(8,0,2008);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(x);
      //--- allocation
      ArrayResize(y,9);
      //--- initialization
      y[0]=4323239600000.0;
      y[1]=4560913100000.0;
      y[2]=5564091500000.0;
      y[3]=6743189300000.0;
      y[4]=7284064600000.0;
      y[5]=7050129600000.0;
      y[6]=7092221500000.0;
      y[7]=8483907600000.0;
      y[8]=8625804400000.0;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(c,3);
      //--- initialization
      c[0]=1.0e+13;
      c[1]=1;
      c[2]=1;
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(c,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(c,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(c,CInfOrNaN::NegativeInfinity());
      double epsx=1.0e-5;
      //--- check
      if(_spoil_scenario==13)
         epsx=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==14)
         epsx=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==15)
         epsx=CInfOrNaN::NegativeInfinity();
      //--- allocation
      ArrayResize(bndl,3);
      //--- initialization
      bndl[0]=-CInfOrNaN::PositiveInfinity();
      bndl[1]=-10;
      bndl[2]=0.1;
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(bndl,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Deleting_Element(bndl);
      //--- allocation
      ArrayResize(bndu,3);
      //--- initialization
      bndu[0]=CInfOrNaN::PositiveInfinity();
      bndu[1]=10;
      bndu[2]=2;
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Value(bndu,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Deleting_Element(bndu);
      //--- allocation
      ArrayResize(s,3);
      //--- initialization
      s[0]=1.0e+12;
      s[1]=1;
      s[2]=1;
      //--- check
      if(_spoil_scenario==20)
         Spoil_Vector_By_Value(s,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==21)
         Spoil_Vector_By_Value(s,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==22)
         Spoil_Vector_By_Value(s,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==23)
         Spoil_Vector_By_Deleting_Element(s);
      //--- create variables
      int maxits=0;
      int info;
      CLSFitStateShell state;
      CLSFitReportShell rep;
      double diffstep=1.0e-5;
      //--- check
      if(_spoil_scenario==24)
         diffstep=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==25)
         diffstep=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==26)
         diffstep=CInfOrNaN::NegativeInfinity();
      //--- function call
      CAlglib::LSFitCreateF(x,y,c,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitFit(state,fdebtfunc,frep,false,obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      CAlglib::LSFitResults(state,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,3);
      //--- initialization
      temparray[0]=4.142560e+12;
      temparray[1]=0.43424;
      temparray[2]=0.565376;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,2);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,-0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Unconstrained (general) linear least squares fitting with and    |
//| withweights                                                      |
//+------------------------------------------------------------------+
void TEST_LSFit_D_Lin(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble fmatrix;
   int           info;
   double        c[];
   double        y[];
   double        temparray[];
   double        w[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
     {
      //--- In this example we demonstrate linear fitting by f(x|a)=a*exp(0.5*x).
      //--- We have:
      //--- * y - vector of experimental data
      //--- * fmatrix -  matrix of basis functions calculated at sample points
      //---              Actually,we have only one basis function F0=exp(0.5*x).
      fmatrix.Resize(11,1);
      //--- initialization
      fmatrix.Set(0,0,0.606531);
      fmatrix.Set(1,0,0.67032);
      fmatrix.Set(2,0,0.740818);
      fmatrix.Set(3,0,0.818731);
      fmatrix.Set(4,0,0.904837);
      fmatrix.Set(5,0,1);
      fmatrix.Set(6,0,1.105171);
      fmatrix.Set(7,0,1.221403);
      fmatrix.Set(8,0,1.349859);
      fmatrix.Set(9,0,1.491825);
      fmatrix.Set(10,0,1.648721);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(fmatrix,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(fmatrix,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(fmatrix,CInfOrNaN::NegativeInfinity());
      //--- allocation
      ArrayResize(y,11);
      //--- initialization
      y[0]=1.133719;
      y[1]=1.306522;
      y[2]=1.504604;
      y[3]=1.554663;
      y[4]=1.884638;
      y[5]=2.072436;
      y[6]=2.257285;
      y[7]=2.534068;
      y[8]=2.622017;
      y[9]=2.897713;
      y[10]=3.219371;
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      //--- create a variable
      CLSFitReportShell rep;
      //--- Linear fitting withweights
      CAlglib::LSFitLinear(y,fmatrix,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1.9865;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.00005);
      //--- Linear fitting with individual weights.
      //--- Slightly different result is returned.
      ArrayResize(w,11);
      //--- initialization
      w[0]=1.414213;
      w[1]=1;
      w[2]=1;
      w[3]=1;
      w[4]=1;
      w[5]=1;
      w[6]=1;
      w[7]=1;
      w[8]=1;
      w[9]=1;
      w[10]=1;
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Adding_Element(w);
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Deleting_Element(w);
      //--- function call
      CAlglib::LSFitLinearW(y,w,fmatrix,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,1);
      //--- initialization
      temparray[0]=1.983354;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Constrained (general) linear least squares fitting with and      |
//| withweights                                                      |
//+------------------------------------------------------------------+
void TEST_LSFit_D_Linc(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   CMatrixDouble fmatrix;
   CMatrixDouble cmatrix;
   double        y[];
   double        c[];
   double        temparray[];
   double        w[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
     {
      //--- In this example we demonstrate linear fitting by f(x|a,b)=a*x+b
      //--- with simple constraint f(0)=0.
      //--- We have:
      //--- * y - vector of experimental data
      //--- * fmatrix -  matrix of basis functions sampled at [0,1] with step 0.2:
      //---                  [ 1.0   0.0 ]
      //---                  [ 1.0   0.2 ]
      //---                  [ 1.0   0.4 ]
      //---                  [ 1.0   0.6 ]
      //---                  [ 1.0   0.8 ]
      //---                  [ 1.0   1.0 ]
      //---              first column contains value of first basis function (constant term)
      //---              second column contains second basis function (linear term)
      //--- * cmatrix -  matrix of linear constraints:
      //---                  [ 1.0  0.0  0.0 ]
      //---              first two columns contain coefficients before basis functions,
      //---              last column contains desired value of their sum.
      //---              So [1,0,0] means "1*constant_term + 0*linear_term=0"
      ArrayResize(y,6);
      //--- initialization
      y[0]=0.072436;
      y[1]=0.246944;
      y[2]=0.491263;
      y[3]=0.5223;
      y[4]=0.714064;
      y[5]=0.921929;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      fmatrix.Resize(6,2);
      //--- initialization
      fmatrix.Set(0,0,1);
      fmatrix.Set(0,1,0);
      fmatrix.Set(1,0,1);
      fmatrix.Set(1,1,0.2);
      fmatrix.Set(2,0,1);
      fmatrix.Set(2,1,0.4);
      fmatrix.Set(3,0,1);
      fmatrix.Set(3,1,0.6);
      fmatrix.Set(4,0,1);
      fmatrix.Set(4,1,0.8);
      fmatrix.Set(5,0,1);
      fmatrix.Set(5,1,1);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Value(fmatrix,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Value(fmatrix,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Matrix_By_Value(fmatrix,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Matrix_By_Adding_Row(fmatrix);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Matrix_By_Adding_Col(fmatrix);
      //--- check
      if(_spoil_scenario==10)
         Spoil_Matrix_By_Deleting_Row(fmatrix);
      //--- check
      if(_spoil_scenario==11)
         Spoil_Matrix_By_Deleting_Col(fmatrix);
      //--- allocation
      cmatrix.Resize(1,3);
      //--- initialization
      cmatrix.Set(0,0,1);
      cmatrix.Set(0,1,0);
      cmatrix.Set(0,2,0);
      //--- check
      if(_spoil_scenario==12)
         Spoil_Matrix_By_Value(cmatrix,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==13)
         Spoil_Matrix_By_Value(cmatrix,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==14)
         Spoil_Matrix_By_Value(cmatrix,CInfOrNaN::NegativeInfinity());
      //--- create variables
      int info;
      CLSFitReportShell rep;
      //--- Constrained fitting withweights
      CAlglib::LSFitLinearC(y,fmatrix,cmatrix,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=0;
      temparray[1]=0.932933;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.0005);
      //--- Constrained fitting with individual weights
      ArrayResize(w,6);
      //--- initialization
      w[0]=1;
      w[1]=1.414213;
      w[2]=1;
      w[3]=1;
      w[4]=1;
      w[5]=1;
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Adding_Element(w);
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Deleting_Element(w);
      //--- function call
      CAlglib::LSFitLinearWC(y,w,fmatrix,cmatrix,info,c,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- allocation
      ArrayResize(temparray,2);
      //--- initialization
      temparray[0]=0;
      temparray[1]=0.938322;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,temparray,0.0005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Unconstrained polynomial fitting                                 |
//+------------------------------------------------------------------+
void TEST_LSFit_D_Pol(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   double w[];
   double xc[];
   double yc[];
   int    dc[];
   int    m;
   double t;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
     {
      //--- This example demonstrates polynomial fitting.
      //--- Fitting is done by two (M=2) functions from polynomial basis:
      //---     f0=1
      //---     f1=x
      //--- Basically,it just a linear fit;more complex polynomials may be used
      //--- (e.g. parabolas with M=3,cubic with M=4),but even such simple fit allows
      //--- us to demonstrate polynomialfit() function in action.
      //--- We have:
      //--- * x      set of abscissas
      //--- * y      experimental data
      //--- Additionally we demonstrate weighted fitting,where second point has
      //--- more weight than other ones.
      ArrayResize(x,11);
      //--- initialization
      x[0]=0;
      x[1]=0.1;
      x[2]=0.2;
      x[3]=0.3;
      x[4]=0.4;
      x[5]=0.5;
      x[6]=0.6;
      x[7]=0.7;
      x[8]=0.8;
      x[9]=0.9;
      x[10]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,11);
      //--- initialization
      y[0]=0;
      y[1]=0.05;
      y[2]=0.26;
      y[3]=0.32;
      y[4]=0.33;
      y[5]=0.43;
      y[6]=0.6;
      y[7]=0.6;
      y[8]=0.77;
      y[9]=0.98;
      y[10]=1.02;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      m=2;
      t=2;
      //--- check
      if(_spoil_scenario==10)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==11)
         t=CInfOrNaN::NegativeInfinity();
      //--- create variables
      int info;
      CBarycentricInterpolantShell p;
      CPolynomialFitReportShell rep;
      //--- Fitting withindividual weights
      //--- NOTE: result is returned as barycentricinterpolant structure.
      //---       if you want to get representation in the power basis,
      //---       you can use barycentricbar2pow() function to convert
      //---       from barycentric to power representation (see docs for
      //---       POLINT subpackage for more info).
      CAlglib::PolynomialFit(x,y,m,info,p,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.011,0.002);
      //--- Fitting with individual weights
      //--- NOTE: slightly different result is returned
      ArrayResize(w,11);
      //--- initialization
      w[0]=1;
      w[1]=1.414213562;
      w[2]=1;
      w[3]=1;
      w[4]=1;
      w[5]=1;
      w[6]=1;
      w[7]=1;
      w[8]=1;
      w[9]=1;
      w[10]=1;
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==14)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Adding_Element(w);
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Deleting_Element(w);
      //--- allocation
      ArrayResize(xc,0);
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Adding_Element(xc);
      //--- allocation
      ArrayResize(yc,0);
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Adding_Element(yc);
      //--- allocation
      ArrayResize(dc,0);
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Adding_Element(dc);
      //--- function call
      CAlglib::PolynomialFitWC(x,y,w,xc,yc,dc,m,info,p,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.023,0.002);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Constrained polynomial fitting                                   |
//+------------------------------------------------------------------+
void TEST_LSFit_D_Polc(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   double w[];
   double xc[];
   double yc[];
   int    dc[];
   int    m;
   int    info;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<29; _spoil_scenario++)
     {
      //--- This example demonstrates polynomial fitting.
      //--- Fitting is done by two (M=2) functions from polynomial basis:
      //---     f0=1
      //---     f1=x
      //--- with simple constraint on function value
      //---     f(0)=0
      //--- Basically,it just a linear fit;more complex polynomials may be used
      //--- (e.g. parabolas with M=3,cubic with M=4),but even such simple fit allows
      //--- us to demonstrate polynomialfit() function in action.
      //--- We have:
      //--- * x      set of abscissas
      //--- * y      experimental data
      //--- * xc     points where constraints are placed
      //--- * yc     constraints on derivatives
      //--- * dc     derivative indices
      //---          (0 means function itself,1 means first derivative)
      ArrayResize(x,2);
      //--- initialization
      x[0]=1;
      x[1]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,2);
      //--- initialization
      y[0]=0.9;
      y[1]=1.1;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(w,2);
      //--- initialization
      w[0]=1;
      w[1]=1;
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==13)
         Spoil_Vector_By_Adding_Element(w);
      //--- check
      if(_spoil_scenario==14)
         Spoil_Vector_By_Deleting_Element(w);
      //--- allocation
      ArrayResize(xc,1);
      //--- initialization
      xc[0]=0;
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(xc,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(xc,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(xc,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Adding_Element(xc);
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Deleting_Element(xc);
      //--- allocation
      ArrayResize(yc,1);
      //--- initialization
      yc[0]=0;
      //--- check
      if(_spoil_scenario==20)
         Spoil_Vector_By_Value(yc,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==21)
         Spoil_Vector_By_Value(yc,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==22)
         Spoil_Vector_By_Value(yc,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==23)
         Spoil_Vector_By_Adding_Element(yc);
      //--- check
      if(_spoil_scenario==24)
         Spoil_Vector_By_Deleting_Element(yc);
      //--- allocation
      ArrayResize(dc,1);
      //--- initialization
      dc[0]=0;
      //--- check
      if(_spoil_scenario==25)
         Spoil_Vector_By_Adding_Element(dc);
      //--- check
      if(_spoil_scenario==26)
         Spoil_Vector_By_Deleting_Element(dc);
      double t=2;
      //--- check
      if(_spoil_scenario==27)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==28)
         t=CInfOrNaN::NegativeInfinity();
      m=2;
      //--- create variables
      CBarycentricInterpolantShell p;
      CPolynomialFitReportShell rep;
      //--- function call
      CAlglib::PolynomialFitWC(x,y,w,xc,yc,dc,m,info,p,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.000,0.001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Unconstrained fitting by penalized regression spline             |
//+------------------------------------------------------------------+
void TEST_LSFit_D_Spline(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   int    info;
   double v;
   double rho;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<19; _spoil_scenario++)
     {
      //--- In this example we demonstrate penalized spline fitting of noisy data
      //--- We have:
      //--- * x - abscissas
      //--- * y - vector of experimental data,straight line with small noise
      ArrayResize(x,10);
      //--- initialization
      x[0]=0;
      x[1]=0.1;
      x[2]=0.2;
      x[3]=0.3;
      x[4]=0.4;
      x[5]=0.5;
      x[6]=0.6;
      x[7]=0.7;
      x[8]=0.8;
      x[9]=0.9;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,10);
      //--- initialization
      y[0]=0.1;
      y[1]=0;
      y[2]=0.3;
      y[3]=0.4;
      y[4]=0.3;
      y[5]=0.4;
      y[6]=0.62;
      y[7]=0.68;
      y[8]=0.75;
      y[9]=0.95;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      //--- create variables
      CSpline1DInterpolantShell s;
      CSpline1DFitReportShell rep;
      //--- Fit with VERY small amount of smoothing (rho=-5.0)
      //--- and large number of basis functions (M=50).
      //--- With such small regularization penalized spline almost fully reproduces function values
      rho=-5.0;
      //--- check
      if(_spoil_scenario==10)
         rho=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==11)
         rho=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==12)
         rho=CInfOrNaN::NegativeInfinity();
      //--- function call
      CAlglib::Spline1DFitPenalized(x,y,50,rho,info,s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      //--- function call
      v=CAlglib::Spline1DCalc(s,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,0.10,0.01);
      //--- Fit with VERY large amount of smoothing (rho=10.0)
      //--- and large number of basis functions (M=50).
      //--- With such regularization our spline should become close to the straight line fit.
      //--- We will compare its value in x=1.0 with results obtained from such fit.
      rho=+10.0;
      //--- check
      if(_spoil_scenario==13)
         rho=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==14)
         rho=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==15)
         rho=CInfOrNaN::NegativeInfinity();
      //--- function call
      CAlglib::Spline1DFitPenalized(x,y,50,rho,info,s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      //--- function call
      v=CAlglib::Spline1DCalc(s,1.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,0.969,0.001);
      //--- In real life applications you may need some moderate degree of fitting,
      //--- so we try to fit once more with rho=3.0.
      rho=+3.0;
      //--- check
      if(_spoil_scenario==16)
         rho=CInfOrNaN::NaN();
      //--- check
      if(_spoil_scenario==17)
         rho=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==18)
         rho=CInfOrNaN::NegativeInfinity();
      //--- function call
      CAlglib::Spline1DFitPenalized(x,y,50,rho,info,s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Int(info,1);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial fitting, full list of parameters.                     |
//+------------------------------------------------------------------+
void TEST_LSFit_T_PolFit_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   int    info;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(x,11);
      //--- initialization
      x[0]=0;
      x[1]=0.1;
      x[2]=0.2;
      x[3]=0.3;
      x[4]=0.4;
      x[5]=0.5;
      x[6]=0.6;
      x[7]=0.7;
      x[8]=0.8;
      x[9]=0.9;
      x[10]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,11);
      //--- initialization
      y[0]=0;
      y[1]=0.05;
      y[2]=0.26;
      y[3]=0.32;
      y[4]=0.33;
      y[5]=0.43;
      y[6]=0.6;
      y[7]=0.6;
      y[8]=0.77;
      y[9]=0.98;
      y[10]=1.02;
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      int m=2;
      double t=2;
      //--- check
      if(_spoil_scenario==8)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==9)
         t=CInfOrNaN::NegativeInfinity();
      //--- create variables
      CBarycentricInterpolantShell p;
      CPolynomialFitReportShell rep;
      //--- function call
      CAlglib::PolynomialFit(x,y,11,m,info,p,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.011,0.002);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial fitting, full list of parameters.                     |
//+------------------------------------------------------------------+
void TEST_LSFit_T_PolFit_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   double w[];
   double xc[];
   double yc[];
   int    dc[];
   int    m;
   double t;
   int    info;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<14; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(x,11);
      //--- initialization
      x[0]=0;
      x[1]=0.1;
      x[2]=0.2;
      x[3]=0.3;
      x[4]=0.4;
      x[5]=0.5;
      x[6]=0.6;
      x[7]=0.7;
      x[8]=0.8;
      x[9]=0.9;
      x[10]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,11);
      //--- initialization
      y[0]=0;
      y[1]=0.05;
      y[2]=0.26;
      y[3]=0.32;
      y[4]=0.33;
      y[5]=0.43;
      y[6]=0.6;
      y[7]=0.6;
      y[8]=0.77;
      y[9]=0.98;
      y[10]=1.02;
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(w,11);
      //--- initialization
      w[0]=1;
      w[1]=1.414213562;
      w[2]=1;
      w[3]=1;
      w[4]=1;
      w[5]=1;
      w[6]=1;
      w[7]=1;
      w[8]=1;
      w[9]=1;
      w[10]=1;
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(w);
      //--- allocation
      ArrayResize(xc,0);
      ArrayResize(yc,0);
      ArrayResize(dc,0);
      //--- initialization
      m=2;
      t=2;
      //--- check
      if(_spoil_scenario==12)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==13)
         t=CInfOrNaN::NegativeInfinity();
      //--- create variables
      CBarycentricInterpolantShell p;
      CPolynomialFitReportShell rep;
      //--- function call
      CAlglib::PolynomialFitWC(x,y,w,11,xc,yc,dc,0,m,info,p,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.023,0.002);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Polynomial fitting,full list of parameters.                      |
//+------------------------------------------------------------------+
void TEST_LSFit_T_PolFit_3(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double x[];
   double y[];
   double w[];
   double xc[];
   double yc[];
   int    dc[];
   int    m;
   double t;
   int    info;
   double v;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<23; _spoil_scenario++)
     {
      //--- allocation
      ArrayResize(x,2);
      //--- initialization
      x[0]=1;
      x[1]=1;
      //--- check
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      //--- allocation
      ArrayResize(y,2);
      //--- initialization
      y[0]=0.9;
      y[1]=1.1;
      //--- check
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      //--- allocation
      ArrayResize(w,2);
      //--- initialization
      w[0]=1;
      w[1]=1;
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(w,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(w,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(w,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(w);
      //--- allocation
      ArrayResize(xc,1);
      //--- initialization
      xc[0]=0;
      //--- check
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(xc,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(xc,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==14)
         Spoil_Vector_By_Value(xc,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==15)
         Spoil_Vector_By_Deleting_Element(xc);
      //--- allocation
      ArrayResize(yc,1);
      //--- initialization
      yc[0]=0;
      //--- check
      if(_spoil_scenario==16)
         Spoil_Vector_By_Value(yc,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(yc,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==18)
         Spoil_Vector_By_Value(yc,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==19)
         Spoil_Vector_By_Deleting_Element(yc);
      //--- allocation
      ArrayResize(dc,1);
      //--- initialization
      dc[0]=0;
      //--- check
      if(_spoil_scenario==20)
         Spoil_Vector_By_Deleting_Element(dc);
      m=2;
      t=2;
      //--- check
      if(_spoil_scenario==21)
         t=CInfOrNaN::PositiveInfinity();
      //--- check
      if(_spoil_scenario==22)
         t=CInfOrNaN::NegativeInfinity();
      //--- create variables
      CBarycentricInterpolantShell p;
      CPolynomialFitReportShell rep;
      //--- function call
      CAlglib::PolynomialFitWC(x,y,w,2,xc,yc,dc,1,m,info,p,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- function call
      v=CAlglib::BarycentricCalc(p,t);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(v,2.000,0.001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, real matrix, short form                 |
//+------------------------------------------------------------------+
void TEST_MatDet_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double        a;
   CMatrixDouble b;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
     {
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,1);
      b.Set(0,1,2);
      b.Set(1,0,2);
      b.Set(1,1,1);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(b);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- function call
      a=CAlglib::RMatrixDet(b);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(a,-3,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, real matrix, full form                  |
//+------------------------------------------------------------------+
void TEST_MatDet_D_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double        a;
   CMatrixDouble b;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
     {
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,5);
      b.Set(0,1,4);
      b.Set(1,0,4);
      b.Set(1,1,5);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- function call
      a=CAlglib::RMatrixDet(b,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(a,9,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, complex matrix, short form              |
//+------------------------------------------------------------------+
void TEST_MatDet_D_3(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex        a;
   complex        tempcomplex1;
   complex        tempcomplex2;
   complex        tempcomplex3;
   CMatrixComplex b;
   complex        cnan;
   complex        cpositiveinfinity;
   complex        cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
     {
      //--- initialization
      tempcomplex1.real=1;
      tempcomplex1.imag=1;
      tempcomplex2.real=1;
      tempcomplex2.imag=-1;
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,tempcomplex1);
      b.Set(0,1,2);
      b.Set(1,0,2);
      b.Set(1,1,tempcomplex2);
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,cnegativeinfinity);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(b);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- function call
      a=CAlglib::CMatrixDet(b);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex3.real=-2;
      tempcomplex3.imag=0;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex(a,tempcomplex3,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, complex matrix, full form               |
//+------------------------------------------------------------------+
void TEST_MatDet_D_4(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex        a;
   complex        tempcomplex1;
   complex        tempcomplex2;
   complex        tempcomplex3;
   CMatrixComplex b;
   complex        cnan;
   complex        cpositiveinfinity;
   complex        cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
     {
      //--- initialization
      tempcomplex1.real=0;
      tempcomplex1.imag=5;
      tempcomplex2.real=0;
      tempcomplex2.imag=4;
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,tempcomplex1);
      b.Set(0,1,4);
      b.Set(1,0,tempcomplex2);
      b.Set(1,1,5);
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,cnegativeinfinity);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- function call
      a=CAlglib::CMatrixDet(b,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex3.real=0;
      tempcomplex3.imag=9;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex(a,tempcomplex3,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, complex matrix with zero imaginary part,|
//| short form                                                       |
//+------------------------------------------------------------------+
void TEST_MatDet_D_5(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex        a;
   complex        tempcomplex;
   CMatrixComplex b;
   complex        cnan;
   complex        cpositiveinfinity;
   complex        cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
     {
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,9);
      b.Set(0,1,1);
      b.Set(1,0,2);
      b.Set(1,1,1);
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,cnegativeinfinity);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(b);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- function call
      a=CAlglib::CMatrixDet(b);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex.real=7;
      tempcomplex.imag=0;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex(a,tempcomplex,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, real matrix, full form                  |
//+------------------------------------------------------------------+
void TEST_MatDet_T_0(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double        a;
   CMatrixDouble b;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
     {
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,3);
      b.Set(0,1,4);
      b.Set(1,0,-4);
      b.Set(1,1,3);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- function call
      a=CAlglib::RMatrixDet(b,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(a,25,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, real matrix, LU, short form             |
//+------------------------------------------------------------------+
void TEST_MatDet_T_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double        a;
   CMatrixDouble b;
   int           p[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,1);
      b.Set(0,1,2);
      b.Set(1,0,2);
      b.Set(1,1,5);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(b);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- allocation
      ArrayResize(p,2);
      //--- initialization
      p[0]=1;
      p[1]=1;
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Adding_Element(p);
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(p);
      //--- function call
      a=CAlglib::RMatrixLUDet(b,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(a,-5,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, real matrix, LU, full form              |
//+------------------------------------------------------------------+
void TEST_MatDet_T_2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   double        a;
   int           p[];
   CMatrixDouble b;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,5);
      b.Set(0,1,4);
      b.Set(1,0,4);
      b.Set(1,1,5);
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NaN());
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,CInfOrNaN::PositiveInfinity());
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,CInfOrNaN::NegativeInfinity());
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- allocation
      ArrayResize(p,2);
      //--- initialization
      p[0]=0;
      p[1]=1;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Deleting_Element(p);
      //--- function call
      a=CAlglib::RMatrixLUDet(b,p,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Real(a,25,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, complex matrix, full form               |
//+------------------------------------------------------------------+
void TEST_MatDet_T_3(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex        a;
   CMatrixComplex b;
   complex        tempcomplex1;
   complex        tempcomplex2;
   complex        cnan;
   complex        cpositiveinfinity;
   complex        cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
     {
      //--- initialization
      tempcomplex1.real=0;
      tempcomplex1.imag=5;
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,tempcomplex1);
      b.Set(0,1,4);
      b.Set(1,0,-4);
      b.Set(1,1,tempcomplex1);
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,cnegativeinfinity);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- function call
      a=CAlglib::CMatrixDet(b,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex2.real=-9;
      tempcomplex2.imag=0;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex(a,tempcomplex2,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, complex matrix, LU, short form          |
//+------------------------------------------------------------------+
void TEST_MatDet_T_4(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex        a;
   int            p[];
   complex        tempcomplex1;
   complex        tempcomplex2;
   CMatrixComplex b;
   complex        cnan;
   complex        cpositiveinfinity;
   complex        cnegativeinfinity;
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- initialization
      tempcomplex1.real=0;
      tempcomplex1.imag=5;
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,1);
      b.Set(0,1,2);
      b.Set(1,0,2);
      b.Set(1,1,tempcomplex1);
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,cnegativeinfinity);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Adding_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Adding_Col(b);
      //--- check
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==6)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- allocation
      ArrayResize(p,2);
      //--- initialization
      p[0]=1;
      p[1]=1;
      //--- check
      if(_spoil_scenario==7)
         Spoil_Vector_By_Adding_Element(p);
      //--- check
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(p);
      //--- function call
      a=CAlglib::CMatrixLUDet(b,p);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex2.real=0;
      tempcomplex2.imag=-5;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex(a,tempcomplex2,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Determinant calculation, complex matrix, LU, full form           |
//+------------------------------------------------------------------+
void TEST_MatDet_T_5(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- create variables
   complex        a;
   complex        tempcomplex1;
   complex        tempcomplex2;
   CMatrixComplex b;
   complex        cnan;
   complex        cpositiveinfinity;
   complex        cnegativeinfinity;
   int            p[];
//--- testing
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- initialization
      tempcomplex1.real=0;
      tempcomplex1.imag=4;
      //--- allocation
      b.Resize(2,2);
      //--- initialization
      b.Set(0,0,5);
      b.Set(0,1,tempcomplex1);
      b.Set(1,0,4);
      b.Set(1,1,5);
      //--- initialization
      cnan.real=CInfOrNaN::NaN();
      cnan.imag=CInfOrNaN::NaN();
      cpositiveinfinity.real=CInfOrNaN::PositiveInfinity();
      cpositiveinfinity.imag=CInfOrNaN::PositiveInfinity();
      cnegativeinfinity.real=CInfOrNaN::NegativeInfinity();
      cnegativeinfinity.imag=CInfOrNaN::NegativeInfinity();
      //--- check
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(b,cnan);
      //--- check
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(b,cpositiveinfinity);
      //--- check
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(b,cnegativeinfinity);
      //--- check
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(b);
      //--- check
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(b);
      //--- allocation
      ArrayResize(p,2);
      //--- initialization
      p[0]=0;
      p[1]=1;
      //--- check
      if(_spoil_scenario==5)
         Spoil_Vector_By_Deleting_Element(p);
      //--- function call
      a=CAlglib::CMatrixLUDet(b,p,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- initialization
      tempcomplex2.real=25;
      tempcomplex2.imag=0;
      //--- check result
      _TestResult=_TestResult && Doc_Test_Complex(a,tempcomplex2,0.0001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
//--- check
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
//--- change total result
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Matrix multiplication (single-threaded)                          |
//+------------------------------------------------------------------+
void TEST_Ablas_D_Gemm(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   _spoil_scenario=-1;
   matrix<double> A={{2,1},{1,3}};
   CMatrixDouble a=A;
   matrix<double> B={{2,1},{0,1}};
   CMatrixDouble b=B;
   CMatrixDouble c=matrix<double>::Zeros(2,2);
//--- rmatrixgemm() function allows us to calculate matrix product C:=A*B or
//--- to perform more general operation, C:=alpha*op1(A)*op2(B)+beta*C,
//--- where A, B, C are rectangular matrices, op(X) can be X or X^T,
//--- alpha and beta are scalars.
//--- This function:
//--- * can apply transposition and/or multiplication by scalar to operands
//--- * can use arbitrary part of matrices A/B (given by submatrix offset)
//--- * can store result into arbitrary part of C
//--- * for performance reasons requires C to be preallocated
//--- Parameters of this function are:
//--- * M, N, K            -   sizes of op1(A) (which is MxK), op2(B) (which
//---                          is KxN) and C (which is MxN)
//--- * Alpha              -   coefficient before A*B
//--- * A, IA, JA          -   matrix A and offset of the submatrix
//--- * OpTypeA            -   transformation type:
//---                          0 - no transformation
//---                          1 - transposition
//--- * B, IB, JB          -   matrix B and offset of the submatrix
//--- * OpTypeB            -   transformation type:
//---                          0 - no transformation
//---                          1 - transposition
//--- * Beta               -   coefficient before C
//--- * C, IC, JC          -   preallocated matrix C and offset of the submatrix
//--- Below we perform simple product C:=A*B (alpha=1, beta=0)
//--- IMPORTANT: this function works with preallocated C, which must be large
//---            enough to store multiplication result.
   int    m=2;
   int    n=2;
   int    k=2;
   double alpha=1.0;
   int    ia=0;
   int    ja=0;
   int    optypea=0;
   int    ib=0;
   int    jb=0;
   int    optypeb=0;
   double beta=0.0;
   int    ic=0;
   int    jc=0;

   CAlglib::RMatrixGemm(m,n,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
     {
      matrix<double> check={{4.0,3.0},{2.0,4.0}};
      CMatrixDouble Check=check;
      if(Func_spoil_scenario(_spoil_scenario,_TestResult))
         _TestResult=_TestResult && Doc_Test_Real_Matrix(c,Check,0.0001);
     }
//--- Now we try to apply some simple transformation to operands: C:=A*B^T
   optypeb=1;
   CAlglib::RMatrixGemm(m,n,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
     {
      matrix<double> check={{5,1},{5,3}};
      CMatrixDouble Check=check;
      if(Func_spoil_scenario(_spoil_scenario,_TestResult))
         _TestResult=_TestResult && Doc_Test_Real_Matrix(c,Check,0.0001);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Symmetric rank-K update (single-threaded)                        |
//+------------------------------------------------------------------+
void TEST_Ablas_D_Syrk(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   _spoil_scenario=-1;
//--- rmatrixsyrk() function allows us to calculate symmetric rank-K update
//--- C := beta*C + alpha*A'*A, where C is square N*N matrix, A is square K*N
//--- matrix, alpha and beta are scalars. It is also possible to update by
//--- adding A*A' instead of A'*A.
//--- Parameters of this function are:
//--- * N, K       -   matrix size
//--- * Alpha      -   coefficient before A
//--- * A, IA, JA  -   matrix and submatrix offsets
//--- * OpTypeA    -   multiplication type:
//---                  * 0 - A*A^T is calculated
//---                  * 2 - A^T*A is calculated
//--- * Beta       -   coefficient before C
//--- * C, IC, JC  -   preallocated input/output matrix and submatrix offsets
//--- * IsUpper    -   whether upper or lower triangle of C is updated;
//---                  this function updates only one half of C, leaving
//---                  other half unchanged (not referenced at all).
//--- Below we will show how to calculate simple product C:=A'*A
//--- NOTE: beta=0 and we do not use previous value of C, but still it
//---       MUST be preallocated.
   int    n=2;
   int    k=1;
   double alpha=1.0;
   int    ia=0;
   int    ja=0;
   int    optypea=2;
   double beta=0.0;
   int    ic=0;
   int    jc=0;
   bool   isupper=true;
   CMatrixDouble a;
     {
      matrix<double> init={{1,2}};
      a=init;
     }
//--- preallocate space to store result
   CMatrixDouble c=matrix<double>::Zeros(2,2);
//--- calculate product, store result into upper part of c
   CAlglib::RMatrixSyrk(n,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isupper);
   if(Func_spoil_scenario(_spoil_scenario,_TestResult))
     {
      //--- output result.
      //--- IMPORTANT: lower triangle of C was NOT updated!
      matrix<double> check={{1,2},{0,4}};
      CMatrixDouble Check=check;
      _TestResult=_TestResult && Doc_Test_Real_Matrix(c,Check,0.0001);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Basis test for complex matrix functions (correctness and presence|
//| of SMP support)                                                  |
//+------------------------------------------------------------------+
void TEST_Ablas_T_Complex(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   _spoil_scenario=-1;
   CMatrixComplex a;
   CMatrixComplex b;
   CMatrixComplex c;
//--- test cmatrixgemm()
     {
      matrix<complex> init={{(0+2i),(0+1i)},{(1+0i),(3+0i)}};
      a=init;
     }
     {
      matrix<complex> init={{(2+0i),(1+0i)},{(0+0i),(1+0i)}};
      b=init;
     }
   c=matrix<complex>::Full(2,2,0);
   CAlglib::CMatrixGemm(2,2,2,(complex)(1.0),a,0,0,0,b,0,0,0,(complex)0.0,c,0,0);
   if(Func_spoil_scenario(_spoil_scenario,_TestResult))
     {
      matrix<complex> check={{(0+4i),(0+3i)},{(2+0i),(4+0i)}};
      CMatrixComplex Check=check;
      _TestResult=_TestResult && Doc_Test_Complex_Matrix(c,Check,0.0001);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Basic operations with sparse matrices                            |
//+------------------------------------------------------------------+
void TEST_Sparse_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<1; _spoil_scenario++)
     {
      //--- This example demonstrates creation/initialization of the sparse matrix
      //--- and matrix-vector multiplication.
      //--- First, we have to create matrix and initialize it. Matrix is initially created
      //--- in the Hash-Table format, which allows convenient initialization. We can modify
      //--- Hash-Table matrix with SparseSet() and sparseadd() functions.
      //--- NOTE: Unlike CRS format, Hash-Table representation allows you to initialize
      //--- elements in the arbitrary order. You may see that we initialize a[0][0] first,
      //--- then move to the second row, and then move back to the first row.
      CSparseMatrix s;
      CAlglib::SparseCreate(2,2,s);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,0,0,2.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,1,1,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,0,1,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseAdd(s,1,1,4.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Now S is equal to
      //---   [ 2 1 ]
      //---   [   5 ]
      //--- Lets check it by reading matrix contents with SparseGet().
      //--- You may see that with SparseGet() you may read both non-zero
      //--- and zero elements.
      double v;
      v=CAlglib::SparseGet(s,0,0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,2.0000,0.005);
      v=CAlglib::SparseGet(s,0,1);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,1.0000,0.005);
      v=CAlglib::SparseGet(s,1,0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,0.0000,0.005);
      v=CAlglib::SparseGet(s,1,1);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,5.0000,0.005);
      //--- After successful creation we can use our matrix for linear operations.
      //--- However, there is one more thing we MUST do before using S in linear
      //--- operations: we have to convert it from HashTable representation (used for
      //--- initialization and dynamic operations) to CRS format with sparseconverttocrs()
      //--- call. If you omit this call, ALGLIB will generate exception on the first
      //--- attempt to use S in linear operations.
      CAlglib::SparseConvertToCRS(s);
      //--- Now S is in the CRS format and we are ready to do linear operations.
      //--- Lets calculate A*x for some x.
      CRowDouble x;
        {
         vector<double> init={1,-1};
         x=init;
        }
      if(_spoil_scenario==0)
         Spoil_Vector_By_Deleting_Element(x);
      CRowDouble y;
      CAlglib::SparseMV(s,x,y);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={1.000,-5.000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Advanced topic: creation in the CRS format.                      |
//+------------------------------------------------------------------+
void TEST_Sparse_D_CRS(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<2; _spoil_scenario++)
     {
      //--- This example demonstrates creation/initialization of the sparse matrix in the
      //--- CRS format.
      //--- Hash-Table format used by default is very convenient (it allows easy
      //--- insertion of elements, automatic memory reallocation), but has
      //--- significant memory and performance overhead. Insertion of one element
      //--- costs hundreds of CPU cycles, and memory consumption is several times
      //--- higher than that of CRS.
      //--- When you work with really large matrices and when you can tell in
      //--- advance how many elements EXACTLY you need, it can be beneficial to
      //--- create matrix in the CRS format from the very beginning.
      //--- If you want to create matrix in the CRS format, you should:
      //--- * use sparsecreatecrs() function
      //--- * know row sizes in advance (number of non-zero entries in the each row)
      //--- * initialize matrix with SparseSet() - another function, sparseadd(), is not allowed
      //--- * initialize elements from left to right, from top to bottom, each
      //---   element is initialized only once.
      CSparseMatrix s;
      CRowInt row_sizes;
        {
         int init[]={2,2,2,1};
         row_sizes=init;
        }
      if(_spoil_scenario==0)
         Spoil_Vector_By_Deleting_Element(row_sizes);
      CAlglib::SparseCreateCRS(4,4,row_sizes,s);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,0,0,2.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,0,1,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,1,1,4.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,1,2,2.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,2,2,3.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,2,3,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(s,3,3,9.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Now S is equal to
      //---   [ 2 1     ]
      //---   [   4 2   ]
      //---   [     3 1 ]
      //---   [       9 ]
      //--- We should point that we have initialized S elements from left to right,
      //--- from top to bottom. CRS representation does NOT allow you to do so in
      //--- the different order. Try to change order of the SparseSet() calls above,
      //--- and you will see that your program generates exception.
      //--- We can check it by reading matrix contents with SparseGet().
      //--- However, you should remember that SparseGet() is inefficient on
      //--- CRS matrices (it may have to pass through all elements of the row
      //--- until it finds element you need).
      double v=CAlglib::SparseGet(s,0,0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,2.0000,0.005);
      v=CAlglib::SparseGet(s,2,3);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,1.0000,0.005);
      //--- you may see that you can read zero elements (which are not stored) with SparseGet()
      v=CAlglib::SparseGet(s,3,2);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,0.0000,0.005);
      //--- After successful creation we can use our matrix for linear operations.
      //--- Lets calculate A*x for some x.
      CRowDouble x;
        {
         double init[]={1,-1,1,-1};
         x=init;
        }
      if(_spoil_scenario==1)
         Spoil_Vector_By_Deleting_Element(x);
      CRowDouble y;
      CAlglib::SparseMV(s,x,y);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={1.000,-2.000,2.000,-9};
      _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Solving positive definite sparse system using Skyline(SKS) solver|                                                                |
//+------------------------------------------------------------------+
void TEST_SolveSKS_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<4; _spoil_scenario++)
     {
      //--- This example demonstrates creation/initialization of the sparse matrix
      //--- in the SKS (Skyline) storage format and solution using SKS-based direct
      //--- solver.
      //--- First, we have to create matrix and initialize it. Matrix is created
      //--- in the SKS format, using fixed bandwidth initialization function.
      //--- Several points should be noted:
      //--- 1. SKS sparse storage format also allows variable bandwidth matrices;
      //---    we just do not want to overcomplicate this example.
      //--- 2. SKS format requires you to specify matrix geometry prior to
      //---    initialization of its elements with SparseSet(). If you specified
      //---    bandwidth=1, you can not change your mind afterwards and call
      //---    SparseSet() for non-existent elements.
      //--- 3. Because SKS solver need just one triangle of SPD matrix, we can
      //---    omit initialization of the lower triangle of our matrix.
      int n=4;
      int bandwidth=1;
      CSparseMatrix s;
      CAlglib::SparseCreateSKSBand(n,n,bandwidth,s);
      CAlglib::SparseSet(s,0,0,2.0);
      CAlglib::SparseSet(s,0,1,1.0);
      CAlglib::SparseSet(s,1,1,3.0);
      CAlglib::SparseSet(s,1,2,1.0);
      CAlglib::SparseSet(s,2,2,3.0);
      CAlglib::SparseSet(s,2,3,1.0);
      CAlglib::SparseSet(s,3,3,2.0);
      //--- Now we have symmetric positive definite 4x4 system width bandwidth=1:
      //---     [ 2 1     ]   [ x0]]   [  4 ]
      //---     [ 1 3 1   ]   [ x1 ]   [ 10 ]
      //---     [   1 3 1 ] * [ x2 ] = [ 15 ]
      //---     [     1 2 ]   [ x3 ]   [ 11 ]
      //--- After successful creation we can call SKS solver.
      CRowDouble b;
        {
         double init[]={4,10,15,11};
         b=init;
        }
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(b,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(b,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(b,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(b);
      CSparseSolverReport rep;
      CRowDouble x;
      bool isuppertriangle=true;
      CAlglib::SparseSPDSolveSKS(s,isuppertriangle,b,x,rep);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={1.0000,2.0000,3.0000,4.0000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Solution of sparse linear systems with CG                        |
//+------------------------------------------------------------------+
void TEST_LinCG_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<4; _spoil_scenario++)
     {
      //--- This example illustrates solution of sparse linear systems with
      //--- conjugate gradient method.
      //--- Suppose that we have linear system A*x=b with sparse symmetric
      //--- positive definite A (represented by SparseMatrix object)
      //---         [ 5 1       ]
      //---         [ 1 7 2     ]
      //---     A = [   2 8 1   ]
      //---         [     1 4 1 ]
      //---         [       1 4 ]
      //--- and right part b
      //---     [  7 ]
      //---     [ 17 ]
      //--- b = [ 14 ]
      //---     [ 10 ]
      //---     [  6 ]
      //--- and we want to solve this system using sparse linear CG. In order
      //--- to do so, we have to create left part (SparseMatrix object) and
      //--- right part (dense array).
      //--- Initially, sparse matrix is created in the Hash-Table format,
      //--- which allows easy initialization, but do not allow matrix to be
      //--- used in the linear solvers. So after construction you should convert
      //--- sparse matrix to CRS format (one suited for linear operations).
      //--- It is important to note that in our example we initialize full
      //--- matrix A, both lower and upper triangles. However, it is symmetric
      //--- and sparse solver needs just one half of the matrix. So you may
      //--- save about half of the space by filling only one of the triangles.
      CSparseMatrix a;
      CAlglib::SparseCreate(5,5,a);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,0,0,5.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,0,1,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,1,0,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,1,1,7.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,1,2,2.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,2,1,2.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,2,2,8.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,2,3,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,3,2,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,3,3,4.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,3,4,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,4,3,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,4,4,4.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Now our matrix is fully initialized, but we have to do one more
      //--- step - convert it from Hash-Table format to CRS format (see
      //--- documentation on sparse matrices for more information about these
      //--- formats).
      //--- If you omit this call, ALGLIB will generate exception on the first
      //--- attempt to use A in linear operations.
      CAlglib::SparseConvertToCRS(a);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Initialization of the right part
      CRowDouble b;
        {
         double init[]={7,17,14,10,6};
         b=init;
        }
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(b,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(b,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(b,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(b);
      //--- Now we have to create linear solver object and to use it for the
      //--- solution of the linear system.
      //--- NOTE: lincgsolvesparse() accepts additional parameter which tells
      //---       what triangle of the symmetric matrix should be used - upper
      //---       or lower. Because we've filled both parts of the matrix, we
      //---       can use any part - upper or lower.
      CLinCGState s;
      CLinCGReport rep;
      CRowDouble x;
      CAlglib::LinCGCreate(5,s);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::LinCGSolveSparse(s,a,true,b);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::LinCGResult(s,x,rep);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,1);
      double check[]={1.000,2.000,1.000,2.000,1.000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Solution of sparse linear systems with CG                        |
//+------------------------------------------------------------------+
void TEST_LinLSQR_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<4; _spoil_scenario++)
     {
      //--- This example illustrates solution of sparse linear least squares problem
      //--- with LSQR algorithm.
      //--- Suppose that we have least squares problem min|A*x-b| with sparse A
      //--- represented by SparseMatrix object
      //---         [ 1 1 ]
      //---         [ 1 1 ]
      //---     A = [ 2 1 ]
      //---         [ 1   ]
      //---         [   1 ]
      //--- and right part b
      //---     [ 4 ]
      //---     [ 2 ]
      //--- b = [ 4 ]
      //---     [ 1 ]
      //---     [ 2 ]
      //--- and we want to solve this system in the least squares sense using
      //--- LSQR algorithm. In order to do so, we have to create left part
      //--- (SparseMatrix object) and right part (dense array).
      //--- Initially, sparse matrix is created in the Hash-Table format,
      //--- which allows easy initialization, but do not allow matrix to be
      //--- used in the linear solvers. So after construction you should convert
      //--- sparse matrix to CRS format (one suited for linear operations).
      CSparseMatrix a;
      CAlglib::SparseCreate(5,2,a);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,0,0,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,0,1,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,1,0,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,1,1,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,2,0,2.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,2,1,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,3,0,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,4,1,1.0);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Now our matrix is fully initialized, but we have to do one more
      //--- step - convert it from Hash-Table format to CRS format (see
      //--- documentation on sparse matrices for more information about these
      //--- formats).
      //--- If you omit this call, ALGLIB will generate exception on the first
      //--- attempt to use A in linear operations.
      CAlglib::SparseConvertToCRS(a);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Initialization of the right part
      CRowDouble b;
        {
         double init[]={4,2,4,1,2};
         b=init;
        }
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(b,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(b,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(b,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(b);
      //--- Now we have to create linear solver object and to use it for the
      //--- solution of the linear system.
      CLinLSQRState s;
      CLinLSQRReport rep;
      CRowDouble x;
      CAlglib::LinLSQRCreate(5,2,s);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::LinLSQRSolveSparse(s,a,b);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::LinLSQRResults(s,x,rep);
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,4);
      double check[]={1.000,2.000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Linearly constrained dense quadratic programming                 |
//+------------------------------------------------------------------+
void TEST_MinQP_D_LC1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<16; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1
      //--- subject to linear constraint x0+x1<=2
      //--- Exact solution is [x0,x1] = [1.5,0.5]
      //--- IMPORTANT: this solver minimizes  following  function:
      //---     f(x) = 0.5*x'*A*x + b'*x.
      //--- Note that quadratic term has 0.5 before it. So if you want to minimize
      //--- quadratic function, you should rewrite it in such way that quadratic term
      //--- is multiplied by 0.5 too.
      //--- For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as
      //---     f(x) = 0.5*(2*x0^2+2*x1^2) + ....
      //--- and pass diag(2,2) as quadratic term - NOT diag(1,1)!
      matrix<double> A={{2,0},{0,2}};
      CMatrixDouble a=A;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(a,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(a,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(a);
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(a);
      double b[];
        {
         double init[]={-6,-4};
         ArrayCopy(b,init);
        }
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(b,AL_NaN);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(b,AL_POSINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(b,AL_NEGINF);
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(b);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      if(_spoil_scenario==12)
         Spoil_Vector_By_Deleting_Element(s);
      matrix<double> C={{1.0,1.0,2.0}};
      CMatrixDouble c=C;
      if(_spoil_scenario==13)
         Spoil_Matrix_By_Value(c,AL_NaN);
      if(_spoil_scenario==14)
         Spoil_Matrix_By_Value(c,AL_POSINF);
      if(_spoil_scenario==15)
         Spoil_Matrix_By_Value(c,AL_NEGINF);
      int Ct[]={-1};
      CRowInt ct=Ct;
      double x[];
      CMinQPStateShell state;
      CMinQPReportShell rep;
      //--- create solver, set quadratic/linear terms
      CAlglib::MinQPCreate(2,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetQuadraticTerm(state,a);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetLinearTerm(state,b);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetLC(state,c,ct);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set scale of the parameters.
      //--- It is strongly recommended that you set scale of your variables.
      //--- Knowing their scales is essential for evaluation of stopping criteria
      //--- and for preconditioning of the algorithm steps.
      //--- You can find more information on scaling at http://www.CAlglib::net/optimization/scaling.php
      //--- NOTE: for convex problems you may try using minqpsetscaleautodiag()
      //---       which automatically determines variable scales.
      CAlglib::MinQPSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Solve problem with BLEIC-based QP solver.
      //--- This solver is intended for problems with moderate (up to 50) number
      //--- of general linear constraints and unlimited number of box constraints.
      //--- Default stopping criteria are used.
      CAlglib::MinQPSetAlgoBLEIC(state,0.0,0.0,0.0,0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={1.500,0.500};
         _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.05);
        }
      //--- Solve problem with DENSE-AUL solver.
      //--- This solver is optimized for problems with up to several thousands of
      //--- variables and large amount of general linear constraints. Problems with
      //--- less than 50 general linear constraints can be efficiently solved with
      //--- BLEIC, problems with box-only constraints can be solved with QuickQP.
      //--- However, DENSE-AUL will work in any (including unconstrained) case.
      //--- Default stopping criteria are used.
      CAlglib::MinQPSetAlgoDenseAUL(state,1.0e-9,1.0e+4,5);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={1.500,0.500};
         _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.05);
        }
      //--- Solve problem with QuickQP solver.
      //--- This solver is intended for medium and large-scale problems with box
      //--- constraints, and...
      //--- ...Oops! It does not support general linear constraints, -5 returned as completion code!
      CAlglib::MinQPSetAlgoQuickQP(state,0.0,0.0,0.0,0,true);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),-5);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Unconstrained sparse quadratic programming                       |
//+------------------------------------------------------------------+
void TEST_MinQP_D_U2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1,
      //--- with quadratic term given by sparse matrix structure.
      //--- Exact solution is [x0,x1] = [3,2]
      //--- We provide algorithm with starting point, although in this case
      //--- (dense matrix, no constraints) it can work without such information.
      //--- IMPORTANT: this solver minimizes  following  function:
      //---     f(x) = 0.5*x'*A*x + b'*x.
      //--- Note that quadratic term has 0.5 before it. So if you want to minimize
      //--- quadratic function, you should rewrite it in such way that quadratic term
      //--- is multiplied by 0.5 too.
      //--- For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as
      //---     f(x) = 0.5*(2*x0^2+2*x1^2) + ....
      //--- and pass diag(2,2) as quadratic term - NOT diag(1,1)!
      CSparseMatrix a;
      double b[];
        {
         double init[]={-6,-4};
         ArrayCopy(b,init);
        }
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(b,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(b,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(b,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(b);
      double x0[];
        {
         double init[]={0,1};
         ArrayCopy(x0,init);
        }
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(x0);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(s);
      double x[];
      CMinQPStateShell state;
      CMinQPReportShell rep;
      //--- initialize sparsematrix structure
      CAlglib::SparseCreate(2,2,0,a);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,0,0,2.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SparseSet(a,1,1,2.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- create solver, set quadratic/linear terms
      CAlglib::MinQPCreate(2,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetQuadraticTermSparse(state,a,true);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetLinearTerm(state,b);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetStartingPoint(state,x0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set scale of the parameters.
      //--- It is strongly recommended that you set scale of your variables.
      //--- Knowing their scales is essential for evaluation of stopping criteria
      //--- and for preconditioning of the algorithm steps.
      //--- You can find more information on scaling at http://www.CAlglib::net/optimization/scaling.php
      //--- NOTE: for convex problems you may try using minqpsetscaleautodiag()
      //---       which automatically determines variable scales.
      CAlglib::MinQPSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Solve problem with BLEIC-based QP solver.
      //--- This solver is intended for problems with moderate (up to 50) number
      //--- of general linear constraints and unlimited number of box constraints.
      //--- It also supports sparse problems.
      //--- Default stopping criteria are used.
      CAlglib::MinQPSetAlgoBLEIC(state,0.0,0.0,0.0,0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={3,2};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonconvex quadratic programming                                  |
//+------------------------------------------------------------------+
void TEST_MinQP_D_NonConvex(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of nonconvex function
      //---     F(x0,x1) = -(x0^2+x1^2)
      //--- subject to constraints x0,x1 in [1.0,2.0]
      //--- Exact solution is [x0,x1] = [2,2].
      //--- Non-convex problems are harded to solve than convex ones, and they
      //--- may have more than one local minimum. However, ALGLIB solves may deal
      //--- with such problems (altough they do not guarantee convergence to
      //--- global minimum).
      //--- IMPORTANT: this solver minimizes  following  function:
      //---     f(x) = 0.5*x'*A*x + b'*x.
      //--- Note that quadratic term has 0.5 before it. So if you want to minimize
      //--- quadratic function, you should rewrite it in such way that quadratic term
      //--- is multiplied by 0.5 too.
      //--- For example, our function is f(x)=-(x0^2+x1^2), but we rewrite it as
      //---     f(x) = 0.5*(-2*x0^2-2*x1^2)
      //--- and pass diag(-2,-2) as quadratic term - NOT diag(-1,-1)!
      matrix<double> A={{-2,0},{0,-2}};
      CMatrixDouble a=A;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(a,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(a,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(a);
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(a);
      double x0[];
        {
         double init[]={1,1};
         ArrayCopy(x0,init);
        }
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(x0);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      if(_spoil_scenario==12)
         Spoil_Vector_By_Deleting_Element(s);
      double bndl[];
        {
         double init[]={1.0,1.0};
         ArrayCopy(bndl,init);
        }
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(bndl,AL_NaN);
      if(_spoil_scenario==14)
         Spoil_Vector_By_Deleting_Element(bndl);
      double bndu[];
        {
         double init[]={2.0,2.0};
         ArrayCopy(bndu,init);
        }
      if(_spoil_scenario==15)
         Spoil_Vector_By_Value(bndu,AL_NaN);
      if(_spoil_scenario==16)
         Spoil_Vector_By_Deleting_Element(bndu);
      double x[];
      CMinQPStateShell state;
      CMinQPReportShell rep;
      //--- create solver, set quadratic/linear terms, constraints
      CAlglib::MinQPCreate(2,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetQuadraticTerm(state,a);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetStartingPoint(state,x0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set scale of the parameters.
      //--- It is strongly recommended that you set scale of your variables.
      //--- Knowing their scales is essential for evaluation of stopping criteria
      //--- and for preconditioning of the algorithm steps.
      //--- You can find more information on scaling at http://www.CAlglib::net/optimization/scaling.php
      //--- NOTE: there also exists minqpsetscaleautodiag() function
      //---       which automatically determines variable scales; however,
      //---       it does NOT work for non-convex problems.
      CAlglib::MinQPSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Solve problem with BLEIC-based QP solver.
      //--- This solver is intended for problems with moderate (up to 50) number
      //--- of general linear constraints and unlimited number of box constraints.
      //--- It may solve non-convex problems as long as they are bounded from
      //--- below under constraints.
      //--- Default stopping criteria are used.
      CAlglib::MinQPSetAlgoBLEIC(state,0.0,0.0,0.0,0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={2,2};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.005);
      //--- Solve problem with DENSE-AUL solver.
      //--- This solver is optimized for problems with up to several thousands of
      //--- variables and large amount of general linear constraints. Problems with
      //--- less than 50 general linear constraints can be efficiently solved with
      //--- BLEIC, problems with box-only constraints can be solved with QuickQP.
      //--- However, DENSE-AUL will work in any (including unconstrained) case.
      //--- Algorithm convergence is guaranteed only for convex case, but you may
      //--- expect that it will work for non-convex problems too (because near the
      //--- solution they are locally convex).
      //--- Default stopping criteria are used.
      CAlglib::MinQPSetAlgoDenseAUL(state,1.0e-9,1.0e+4,5);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.005);
      //--- Hmm... this problem is bounded from below (has solution) only under constraints.
      //--- What it we remove them?
      //--- You may see that BLEIC algorithm detects unboundedness of the problem,
      //--- -4 is returned as completion code. However, DENSE-AUL is unable to detect
      //--- such situation and it will cycle forever (we do not test it here).
      double nobndl[];
        {
         double init[]={-AL_POSINF,-AL_POSINF};
         ArrayCopy(nobndl,init);
        }
      if(_spoil_scenario==17)
         Spoil_Vector_By_Value(nobndl,AL_NaN);
      if(_spoil_scenario==18)
         Spoil_Vector_By_Deleting_Element(nobndl);
      double nobndu[];
        {
         double init[]={AL_POSINF,+AL_POSINF};
         ArrayCopy(nobndu,init);
        }
      if(_spoil_scenario==19)
         Spoil_Vector_By_Value(nobndu,AL_NaN);
      if(_spoil_scenario==20)
         Spoil_Vector_By_Deleting_Element(nobndu);
      CAlglib::MinQPSetBC(state,nobndl,nobndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPSetAlgoBLEIC(state,0.0,0.0,0.0,0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinQPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.GetTerminationType(),-4);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Basic linear programming example                                 |
//+------------------------------------------------------------------+
void TEST_MinLP_Basic(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
     {
      //--- This example demonstrates how to minimize
      //---     F(x0,x1) = -0.1*x0 - x1
      //--- subject to box constraints
      //---     -1 <= x0,x1 <= +1
      //--- and general linear constraints
      //---     x0 - x1 >= -1
      //---     x0 + x1 <=  1
      //--- We use dual simplex solver provided by ALGLIB for this task. Box
      //--- constraints are specified by means of constraint vectors bndl and
      //--- bndu (we have bndl<=x<=bndu). General linear constraints are
      //--- specified as AL<=A*x<=AU, with AL/AU being 2x1 vectors and A being
      //--- 2x2 matrix.
      //--- NOTE: some/all components of AL/AU can be +-INF, same applies to
      //---       bndl/bndu. You can also have AL[I]=AU[i] (as well as
      //---       BndL[i]=BndU[i]).
      matrix<double> A={{1,-1},{1,+1}};
      CMatrixDouble a=A;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(a,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Deleting_Row(a);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Deleting_Col(a);
      double Al[]={-1,-AL_POSINF};
      CRowDouble al=Al;
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(al,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(al);
      double Au[]={AL_POSINF,+1};
      CRowDouble au=Au;
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(au,AL_NaN);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Deleting_Element(au);
      double C[]={-0.1,-1};
      CRowDouble c=C;
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(c,AL_NaN);
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(c);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Deleting_Element(s);
      CRowDouble bndl=vector<double>::Full(2,-1);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Value(bndl,AL_NaN);
      if(_spoil_scenario==12)
         Spoil_Vector_By_Deleting_Element(bndl);
      CRowDouble bndu=vector<double>::Ones(2);
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(bndu,AL_NaN);
      if(_spoil_scenario==14)
         Spoil_Vector_By_Deleting_Element(bndu);
      CRowDouble x;
      CMinLPState state;
      CMinLPReport rep;
      CAlglib::MinLPCreate(2,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set cost vector, box constraints, general linear constraints.
      //--- Box constraints can be set in one call to minlpsetbc() or minlpsetbcall()
      //--- (latter sets same constraints for all variables and accepts two scalars
      //--- instead of two vectors).
      //--- General linear constraints can be specified in several ways:
      //--- * minlpsetlc2dense() - accepts dense 2D array as input; sometimes this
      //---   approach is more convenient, although less memory-efficient.
      //--- * minlpsetlc2() - accepts sparse matrix as input
      //--- * minlpaddlc2dense() - appends one row to the current set of constraints;
      //---   row being appended is specified as dense vector
      //--- * minlpaddlc2() - appends one row to the current set of constraints;
      //---   row being appended is specified as sparse set of elements
      //--- Independently from specific function being used, LP solver uses sparse
      //--- storage format for internal representation of constraints.
      CAlglib::MinLPSetCost(state,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinLPSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinLPSetLC2Dense(state,a,al,au,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set scale of the parameters.
      //--- It is strongly recommended that you set scale of your variables.
      //--- Knowing their scales is essential for evaluation of stopping criteria
      //--- and for preconditioning of the algorithm steps.
      //--- You can find more information on scaling at http://www.CAlglib::net/optimization/scaling.php
      CAlglib::MinLPSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Solve
      CAlglib::MinLPOptimize(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinLPResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={0,1};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.0005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinearly constrained optimization (inequality constraints)    |
//+------------------------------------------------------------------+
void TEST_MinNLC_D_Inequality(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_NLCFunc1_Jac Jac;
   CNDimensional_Rep Rep;
   CObject Obj;
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x0,x1) = -x0+x1
      //--- subject to box constraints
      //---    x0>=0, x1>=0
      //--- and nonlinear inequality constraint
      //---    x0^2 + x1^2 - 1 <= 0
      CRowDouble x0=vector<double>::Zeros(2);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      double epsx=0.000001;
      if(_spoil_scenario==6)
         epsx=AL_NaN;
      if(_spoil_scenario==7)
         epsx=AL_POSINF;
      if(_spoil_scenario==8)
         epsx=AL_NEGINF;
      int maxits=0;
      CRowDouble bndl=vector<double>::Zeros(2);
      CRowDouble bndu=vector<double>::Full(2,AL_POSINF);
      CMinNLCState state;
      //--- Create optimizer object and tune its settings:
      //--- * epsx=0.000001  stopping condition for inner iterations
      //--- * s=[1,1]        all variables have unit scale; it is important to
      //---                  tell optimizer about scales of your variables - it
      //---                  greatly accelerates convergence and helps to perform
      //---                  some important integrity checks.
      CAlglib::MinNLCCreate(2,x0,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Choose one of the nonlinear programming solvers supported by minnlc
      //--- optimizer:
      //--- * SQP - sequential quadratic programming NLP solver
      //--- * AUL - augmented Lagrangian NLP solver
      //--- * SLP - successive linear programming NLP solver
      //--- Different solvers have different properties:
      //--- * SQP needs less function evaluations than any other solver, but it
      //---   has much higher iteration cost than other solvers (a QP subproblem
      //---   has to be solved during each step)
      //--- * AUL solver has cheaper iterations, but needs more target function
      //---   evaluations
      //--- * SLP is the most robust solver provided by ALGLIB, but it performs
      //---   order of magnitude more iterations than SQP.
      //--- In the code below we set solver to be AUL but then override it with SLP,
      //--- and then with SQP, so the effective choice is to use SLP. We recommend
      //--- you to use SQP at least for early prototyping stages, and then switch
      //--- to AUL if possible.
      double rho=1000.0;
      int outerits=5;
      CAlglib::MinNLCSetAlgoAUL(state,rho,outerits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetAlgoSLP(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetAlgoSQP(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set constraints:
      //--- 1. boundary constraints are passed with minnlcsetbc() call
      //--- 2. nonlinear constraints are more tricky-you can not "pack" general
      //---    nonlinear function into double precision array. That's why
      //---    MinNLCSetNLC() does not accept constraints itself - only constraint
      //---    counts are passed: first parameter is number of equality constraints,
      //---    second one is number of inequality constraints.
      //---    As for constraining functions - these functions are passed as part
      //---    of problem Jacobian (see below).
      //--- NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
      //---       linear and general nonlinear constraints. This example does not
      //---       show how to work with general linear constraints, but you can
      //---       easily find it in documentation on minnlcsetlc() function.
      CAlglib::MinNLCSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetNLC(state,0,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Activate OptGuard integrity checking.
      //--- OptGuard monitor helps to catch common coding and problem statement
      //--- issues, like:
      //--- * discontinuity of the target/constraints (C0 continuity violation)
      //--- * nonsmoothness of the target/constraints (C1 continuity violation)
      //--- * erroneous analytic Jacobian, i.e. one inconsistent with actual
      //---   change in the target/constraints
      //--- OptGuard is essential for early prototyping stages because such
      //--- problems often result in premature termination of the optimizer
      //--- which is really hard to distinguish from the correct termination.
      //--- IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
      //---            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
      //---            Other OptGuard checks add moderate overhead, but anyway
      //---            it is better to turn them off when they are not needed.
      CAlglib::MinNLCOptGuardSmoothness(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCOptGuardGradient(state,0.001);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and test results.
      //--- Optimizer object accepts vector function and its Jacobian, with first
      //--- component (Jacobian row) being target function, and next components
      //--- (Jacobian rows) being nonlinear equality and inequality constraints.
      //--- So, our vector function has form
      //---     {f0,f1} = { -x0+x1 , x0^2+x1^2-1 }
      //--- with Jacobian
      //---         [  -1    +1  ]
      //---     J = [            ]
      //---         [ 2*x0  2*x1 ]
      //--- with f0 being target function, f1 being constraining function. Number
      //--- of equality/inequality constraints is specified by minnlcsetnlc(),
      //--- with equality ones always being first, inequality ones being last.
      CMinNLCReport rep;
      CRowDouble x1;
      CAlglib::MinNLCOptimize(state,Jac,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCResults(state,x1,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={1.0000,0.0000};
         _TestResult=_TestResult && Doc_Test_Real_Vector(x1,check,0.005);
        }
      //--- Check that OptGuard did not report errors
      //--- NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
      //---       1.0 to some of its components.
      COptGuardReport ogrep;
      CAlglib::MinNLCOptGuardResults(state,ogrep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_badgradsuspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc0suspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc1suspected,false);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinearly constrained optimization (equality constraints)      |
//+------------------------------------------------------------------+
void TEST_MinNLC_D_Equality(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_NLCFunc1_Jac Jac;
   CNDimensional_Rep Rep;
   CObject Obj;
//CAp::TraceFile("IPM.DETAILED,SLP.DETAILED,SLP.PROBING,SQP.DETAILED,SQP.PROBING,AGS.DETAILED.SAMPLE,OPTGUARD.ALWAYS,OPTIMIZERS.X,DSS.DETAILED,RBF.DETAILED,SCHOLESKY.SS", "alglib_interfaces.trs");
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x0,x1) = -x0+x1
      //--- subject to nonlinear equality constraint
      //---    x0^2 + x1^2 - 1 = 0
      CRowDouble x0=vector<double>::Zeros(2);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      double epsx=0.000001;
      if(_spoil_scenario==6)
         epsx=AL_NaN;
      if(_spoil_scenario==7)
         epsx=AL_POSINF;
      if(_spoil_scenario==8)
         epsx=AL_NEGINF;
      int maxits=0;
      CMinNLCState state;
      //--- Create optimizer object and tune its settings:
      //--- * epsx=0.000001  stopping condition for inner iterations
      //--- * s=[1,1]        all variables have unit scale
      CAlglib::MinNLCCreate(2,x0,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Choose one of the nonlinear programming solvers supported by minnlc
      //--- optimizer:
      //--- * SLP - successive linear programming NLP solver
      //--- * AUL - augmented Lagrangian NLP solver
      //--- Different solvers have different properties:
      //--- * SLP is the most robust solver provided by ALGLIB: it can solve both
      //---   convex and nonconvex optimization problems, it respects box and
      //---   linear constraints (after you find feasible point it won't move away
      //---   from the feasible area) and tries to respect nonlinear constraints
      //---   as much as possible. It also usually needs less function evaluations
      //---   to converge than AUL.
      //---   However, it solves LP subproblems at each iterations which adds
      //---   significant overhead to its running time. Sometimes it can be as much
      //---   as 7x times slower than AUL.
      //--- * AUL solver is less robust than SLP - it can violate box and linear
      //---   constraints at any moment, and it is intended for convex optimization
      //---   problems (although in many cases it can deal with nonconvex ones too).
      //---   Also, unlike SLP it needs some tuning (penalty factor and number of
      //---   outer iterations).
      //---   However, it is often much faster than the current version of SLP.
      //--- In the code below we set solver to be AUL but then override it with SLP,
      //--- so the effective choice is to use SLP. We recommend you to use SLP at
      //--- least for early prototyping stages.
      //--- You can comment line with SLP if you want to solve your problem with
      //--- AUL solver.
      double rho=1000.0;
      int outerits=5;
      CAlglib::MinNLCSetAlgoAUL(state,rho,outerits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetAlgoSLP(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set constraints:
      //--- Nonlinear constraints are tricky-you can not "pack" general
      //--- nonlinear function into double precision array. That's why
      //--- MinNLCSetNLC() does not accept constraints itself - only constraint
      //--- counts are passed: first parameter is number of equality constraints,
      //--- second one is number of inequality constraints.
      //--- As for constraining functions - these functions are passed as part
      //--- of problem Jacobian (see below).
      //--- NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
      //---       linear and general nonlinear constraints. This example does not
      //---       show how to work with general linear constraints, but you can
      //---       easily find it in documentation on minnlcsetbc() and
      //---       minnlcsetlc() functions.
      CAlglib::MinNLCSetNLC(state,1,0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Activate OptGuard integrity checking.
      //--- OptGuard monitor helps to catch common coding and problem statement
      //--- issues, like:
      //--- * discontinuity of the target/constraints (C0 continuity violation)
      //--- * nonsmoothness of the target/constraints (C1 continuity violation)
      //--- * erroneous analytic Jacobian, i.e. one inconsistent with actual
      //---   change in the target/constraints
      //--- OptGuard is essential for early prototyping stages because such
      //--- problems often result in premature termination of the optimizer
      //--- which is really hard to distinguish from the correct termination.
      //--- IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
      //---            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
      //---            Other OptGuard checks add moderate overhead, but anyway
      //---            it is better to turn them off when they are not needed.
      CAlglib::MinNLCOptGuardSmoothness(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCOptGuardGradient(state,0.001);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and test results.
      //--- Optimizer object accepts vector function and its Jacobian, with first
      //--- component (Jacobian row) being target function, and next components
      //--- (Jacobian rows) being nonlinear equality and inequality constraints.
      //--- So, our vector function has form
      //---     {f0,f1} = { -x0+x1 , x0^2+x1^2-1 }
      //--- with Jacobian
      //---         [  -1    +1  ]
      //---     J = [            ]
      //---         [ 2*x0  2*x1 ]
      //--- with f0 being target function, f1 being constraining function. Number
      //--- of equality/inequality constraints is specified by minnlcsetnlc(),
      //--- with equality ones always being first, inequality ones being last.
      CMinNLCReport rep;
      CRowDouble x1;
      CAlglib::MinNLCOptimize(state,Jac,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCResults(state,x1,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={0.70710,-0.70710};
         _TestResult=_TestResult && Doc_Test_Real_Vector(x1,check,0.005);
        }
      //--- Check that OptGuard did not report errors
      //--- NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
      //---       1.0 to some of its components.
      COptGuardReport ogrep;
      CAlglib::MinNLCOptGuardResults(state,ogrep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_badgradsuspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc0suspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc1suspected,false);
      _TestResult=_TestResult && (_spoil_scenario==-1);
      //CAp::TraceDisable();
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
//CAp::TraceDisable();
  }
//+------------------------------------------------------------------+
//| Nonlinearly constrained optimization with mixed equality/        |
//| inequality constraints                                           |
//+------------------------------------------------------------------+
void TEST_MinNLC_D_Mixed(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_NLCFunc2_Jac Jac;
   CNDimensional_Rep Rep;
   CObject Obj;
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x0,x1) = x0+x1
      //--- subject to nonlinear inequality constraint
      //---    x0^2 + x1^2 - 1 <= 0
      //--- and nonlinear equality constraint
      //---    x2-exp(x0) = 0
      CRowDouble x0=vector<double>::Zeros(3);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(3);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      double epsx=0.000001;
      if(_spoil_scenario==6)
         epsx=AL_NaN;
      if(_spoil_scenario==7)
         epsx=AL_POSINF;
      if(_spoil_scenario==8)
         epsx=AL_NEGINF;
      int maxits=0;
      CMinNLCState state;
      CMinNLCReport rep;
      CRowDouble x1;
      //--- Create optimizer object and tune its settings:
      //--- * epsx=0.000001  stopping condition for inner iterations
      //--- * s=[1,1]        all variables have unit scale
      //--- * upper limit on step length is specified (to avoid probing locations where exp() is large)
      CAlglib::MinNLCCreate(3,x0,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetSTPMax(state,10.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Choose one of the nonlinear programming solvers supported by minnlc
      //--- optimizer:
      //--- * SLP - successive linear programming NLP solver
      //--- * AUL - augmented Lagrangian NLP solver
      //--- Different solvers have different properties:
      //--- * SLP is the most robust solver provided by ALGLIB: it can solve both
      //---   convex and nonconvex optimization problems, it respects box and
      //---   linear constraints (after you find feasible point it won't move away
      //---   from the feasible area) and tries to respect nonlinear constraints
      //---   as much as possible. It also usually needs less function evaluations
      //---   to converge than AUL.
      //---   However, it solves LP subproblems at each iterations which adds
      //---   significant overhead to its running time. Sometimes it can be as much
      //---   as 7x times slower than AUL.
      //--- * AUL solver is less robust than SLP - it can violate box and linear
      //---   constraints at any moment, and it is intended for convex optimization
      //---   problems (although in many cases it can deal with nonconvex ones too).
      //---   Also, unlike SLP it needs some tuning (penalty factor and number of
      //---   outer iterations).
      //---   However, it is often much faster than the current version of SLP.
      //--- In the code below we set solver to be AUL but then override it with SLP,
      //--- so the effective choice is to use SLP. We recommend you to use SLP at
      //--- least for early prototyping stages.
      //--- You can comment line with SLP if you want to solve your problem with
      //--- AUL solver.
      double rho=1000.0;
      int outerits=5;
      CAlglib::MinNLCSetAlgoAUL(state,rho,outerits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCSetAlgoSLP(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set constraints:
      //--- Nonlinear constraints are tricky-you can not "pack" general
      //--- nonlinear function into double precision array. That's why
      //--- minnlcsetnlc() does not accept constraints itself - only constraint
      //--- counts are passed: first parameter is number of equality constraints,
      //--- second one is number of inequality constraints.
      //--- As for constraining functions - these functions are passed as part
      //--- of problem Jacobian (see below).
      //--- NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
      //---       linear and general nonlinear constraints. This example does not
      //---       show how to work with boundary or general linear constraints, but you
      //---       can easily find it in documentation on minnlcsetbc() and
      //---       minnlcsetlc() functions.
      CAlglib::MinNLCSetNLC(state,1,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Activate OptGuard integrity checking.
      //--- OptGuard monitor helps to catch common coding and problem statement
      //--- issues, like:
      //--- * discontinuity of the target/constraints (C0 continuity violation)
      //--- * nonsmoothness of the target/constraints (C1 continuity violation)
      //--- * erroneous analytic Jacobian, i.e. one inconsistent with actual
      //---   change in the target/constraints
      //--- OptGuard is essential for early prototyping stages because such
      //--- problems often result in premature termination of the optimizer
      //--- which is really hard to distinguish from the correct termination.
      //--- IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
      //---            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
      //---            Other OptGuard checks add moderate overhead, but anyway
      //---            it is better to turn them off when they are not needed.
      CAlglib::MinNLCOptGuardSmoothness(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCOptGuardGradient(state,0.001);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and test results.
      //--- Optimizer object accepts vector function and its Jacobian, with first
      //--- component (Jacobian row) being target function, and next components
      //--- (Jacobian rows) being nonlinear equality and inequality constraints.
      //--- So, our vector function has form
      //---     {f0,f1,f2} = { x0+x1 , x2-exp(x0) , x0^2+x1^2-1 }
      //--- with Jacobian
      //---         [  +1      +1       0 ]
      //---     J = [-exp(x0)  0        1 ]
      //---         [ 2*x0    2*x1      0 ]
      //--- with f0 being target function, f1 being equality constraint "f1=0",
      //--- f2 being inequality constraint "f2<=0". Number of equality/inequality
      //--- constraints is specified by minnlcsetnlc(), with equality ones always
      //--- being first, inequality ones being last.
      CAlglib::MinNLCOptimize(state,Jac,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNLCResults(state,x1,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={-0.70710,-0.70710,0.49306};
         _TestResult=_TestResult && Doc_Test_Real_Vector(x1,check,0.005);
        }
      //--- Check that OptGuard did not report errors
      //--- NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
      //---       1.0 to some of its components.
      COptGuardReport ogrep;
      CAlglib::MinNLCOptGuardResults(state,ogrep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_badgradsuspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc0suspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc1suspected,false);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonsmooth unconstrained optimization                             |
//+------------------------------------------------------------------+
void TEST_MinNS_D_Unconstrained(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_NSFunc1_Jac Jac;
   CNDimensional_Rep Rep;
   CObject Obj;
   for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x0,x1) = 2*|x0|+|x1|
      //--- using nonsmooth nonlinear optimizer.
      CRowDouble x0=vector<double>::Ones(2);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      double epsx=0.00001;
      if(_spoil_scenario==6)
         epsx=AL_NaN;
      if(_spoil_scenario==7)
         epsx=AL_POSINF;
      if(_spoil_scenario==8)
         epsx=AL_NEGINF;
      double radius=0.1;
      if(_spoil_scenario==9)
         radius=AL_NaN;
      if(_spoil_scenario==10)
         radius=AL_POSINF;
      if(_spoil_scenario==11)
         radius=AL_NEGINF;
      double rho=0.0;
      if(_spoil_scenario==12)
         rho=AL_NaN;
      if(_spoil_scenario==13)
         rho=AL_POSINF;
      if(_spoil_scenario==14)
         rho=AL_NEGINF;
      int maxits=0;
      CMinNSState state;
      CMinNSReport rep;
      CRowDouble x1;
      //--- Create optimizer object, choose AGS algorithm and tune its settings:
      //--- * radius=0.1     good initial value; will be automatically decreased later.
      //--- * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
      //---                  because we do not have such constraints
      //--- * epsx=0.000001  stopping conditions
      //--- * s=[1,1]        all variables have unit scale
      CAlglib::MinNSCreate(2,x0,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetAlgoAGS(state,radius,rho);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and test results.
      //--- Optimizer object accepts vector function and its Jacobian, with first
      //--- component (Jacobian row) being target function, and next components
      //--- (Jacobian rows) being nonlinear equality and inequality constraints
      //--- (box/linear ones are passed separately by means of minnssetbc() and
      //--- minnssetlc() calls).
      //--- If you do not have nonlinear constraints (exactly our situation), then
      //--- you will have one-component function vector and 1xN Jacobian matrix.
      //--- So, our vector function has form
      //---     {f0} = { 2*|x0|+|x1| }
      //--- with Jacobian
      //---         [                       ]
      //---     J = [ 2*sign(x0)   sign(x1) ]
      //---         [                       ]
      //--- NOTE: nonsmooth optimizer requires considerably more function
      //---       evaluations than smooth solver - about 2N times more. Using
      //---       numerical differentiation introduces additional (multiplicative)
      //---       2N speedup.
      //---       It means that if smooth optimizer WITH user-supplied gradient
      //---       needs 100 function evaluations to solve 50-dimensional problem,
      //---       then AGS solver with user-supplied gradient will need about 10.000
      //---       function evaluations, and with numerical gradient about 1.000.000
      //---       function evaluations will be performed.
      //--- NOTE: AGS solver used by us can handle nonsmooth and nonconvex
      //---       optimization problems. It has convergence guarantees, i.e. it will
      //---       converge to stationary point of the function after running for some
      //---       time.
      //---       However, it is important to remember that "stationary point" is not
      //---       equal to "solution". If your problem is convex, everything is OK.
      //---       But nonconvex optimization problems may have "flat spots" - large
      //---       areas where gradient is exactly zero, but function value is far away
      //---       from optimal. Such areas are stationary points too, and optimizer
      //---       may be trapped here.
      //---       "Flat spots" are nonsmooth equivalent of the saddle points, but with
      //---       orders of magnitude worse properties - they may be quite large and
      //---       hard to avoid. All nonsmooth optimizers are prone to this kind of the
      //---       problem, because it is impossible to automatically distinguish "flat
      //---       spot" from true solution.
      //---       This note is here to warn you that you should be very careful when
      //---       you solve nonsmooth optimization problems. Visual inspection of
      //---       results is essential.
      CAlglib::MinNSOptimize(state,Jac,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSResults(state,x1,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real_Vector(x1,vector<double>::Zeros(2),0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonsmooth unconstrained optimization with numerical              |
//| differentiation                                                  |
//+------------------------------------------------------------------+
void TEST_MinNS_D_Diff(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_NSFunc1_FVec FVec;
   CNDimensional_Rep Rep;
   CObject Obj;
   for(_spoil_scenario=-1; _spoil_scenario<18; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x0,x1) = 2*|x0|+|x1|
      //--- using nonsmooth nonlinear optimizer with numerical
      //--- differentiation provided by ALGLIB.
      //--- NOTE: nonsmooth optimizer requires considerably more function
      //---       evaluations than smooth solver - about 2N times more. Using
      //---       numerical differentiation introduces additional (multiplicative)
      //---       2N speedup.
      //---       It means that if smooth optimizer WITH user-supplied gradient
      //---       needs 100 function evaluations to solve 50-dimensional problem,
      //---       then AGS solver with user-supplied gradient will need about 10.000
      //---       function evaluations, and with numerical gradient about 1.000.000
      //---       function evaluations will be performed.
      CRowDouble x0=vector<double>::Ones(2);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      double epsx=0.00001;
      if(_spoil_scenario==6)
         epsx=AL_NaN;
      if(_spoil_scenario==7)
         epsx=AL_POSINF;
      if(_spoil_scenario==8)
         epsx=AL_NEGINF;
      double diffstep=0.000001;
      if(_spoil_scenario==9)
         diffstep=AL_NaN;
      if(_spoil_scenario==10)
         diffstep=AL_POSINF;
      if(_spoil_scenario==11)
         diffstep=AL_NEGINF;
      double radius=0.1;
      if(_spoil_scenario==12)
         radius=AL_NaN;
      if(_spoil_scenario==13)
         radius=AL_POSINF;
      if(_spoil_scenario==14)
         radius=AL_NEGINF;
      double rho=0.0;
      if(_spoil_scenario==15)
         rho=AL_NaN;
      if(_spoil_scenario==16)
         rho=AL_POSINF;
      if(_spoil_scenario==17)
         rho=AL_NEGINF;
      int maxits=0;
      CMinNSState state;
      CMinNSReport rep;
      CRowDouble x1;
      //--- Create optimizer object, choose AGS algorithm and tune its settings:
      //--- * radius=0.1     good initial value; will be automatically decreased later.
      //--- * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
      //---                  because we do not have such constraints
      //--- * epsx=0.000001  stopping conditions
      //--- * s=[1,1]        all variables have unit scale
      CAlglib::MinNSCreateF(2,x0,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetAlgoAGS(state,radius,rho);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and test results.
      //--- Optimizer object accepts vector function, with first component
      //--- being target function, and next components being nonlinear equality
      //--- and inequality constraints (box/linear ones are passed separately
      //--- by means of minnssetbc() and minnssetlc() calls).
      //--- If you do not have nonlinear constraints (exactly our situation), then
      //--- you will have one-component function vector.
      //--- So, our vector function has form
      //---     {f0} = { 2*|x0|+|x1| }
      CAlglib::MinNSOptimize(state,FVec,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSResults(state,x1,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real_Vector(x1,vector<double>::Zeros(2),0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void TEST_MinNS_D_BC(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_NSFunc1_Jac Jac;
   CNDimensional_Rep Rep;
   CObject Obj;
   for(_spoil_scenario=-1; _spoil_scenario<17; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x0,x1) = 2*|x0|+|x1|
      //--- subject to box constraints
      //---        1 <= x0 < +INF
      //---     -INF <= x1 < +INF
      //--- using nonsmooth nonlinear optimizer.
      CRowDouble x0=vector<double>::Ones(2);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      vector<double> Bndl={1,-AL_POSINF};
      CRowDouble bndl=Bndl;
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(bndl,AL_NaN);
      CRowDouble bndu=vector<double>::Full(2,AL_POSINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(bndu,AL_NaN);
      double epsx=0.00001;
      if(_spoil_scenario==8)
         epsx=AL_NaN;
      if(_spoil_scenario==9)
         epsx=AL_POSINF;
      if(_spoil_scenario==10)
         epsx=AL_NEGINF;
      double radius=0.1;
      if(_spoil_scenario==11)
         radius=AL_NaN;
      if(_spoil_scenario==12)
         radius=AL_POSINF;
      if(_spoil_scenario==13)
         radius=AL_NEGINF;
      double rho=0.0;
      if(_spoil_scenario==14)
         rho=AL_NaN;
      if(_spoil_scenario==15)
         rho=AL_POSINF;
      if(_spoil_scenario==16)
         rho=AL_NEGINF;
      int maxits=0;
      CMinNSState state;
      CMinNSReport rep;
      CRowDouble x1;
      //--- Create optimizer object, choose AGS algorithm and tune its settings:
      //--- * radius=0.1     good initial value; will be automatically decreased later.
      //--- * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
      //---                  because we do not have such constraints
      //--- * epsx=0.000001  stopping conditions
      //--- * s=[1,1]        all variables have unit scale
      CAlglib::MinNSCreate(2,x0,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetAlgoAGS(state,radius,rho);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set box constraints.
      //--- General linear constraints are set in similar way (see comments on
      //--- minnssetlc() function for more information).
      //--- You may combine box, linear and nonlinear constraints in one optimization
      //--- problem.
      CAlglib::MinNSSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and test results.
      //--- Optimizer object accepts vector function and its Jacobian, with first
      //--- component (Jacobian row) being target function, and next components
      //--- (Jacobian rows) being nonlinear equality and inequality constraints
      //--- (box/linear ones are passed separately by means of minnssetbc() and
      //--- minnssetlc() calls).
      //--- If you do not have nonlinear constraints (exactly our situation), then
      //--- you will have one-component function vector and 1xN Jacobian matrix.
      //--- So, our vector function has form
      //---     {f0} = { 2*|x0|+|x1| }
      //--- with Jacobian
      //---         [                       ]
      //---     J = [ 2*sign(x0)   sign(x1) ]
      //---         [                       ]
      //--- NOTE: nonsmooth optimizer requires considerably more function
      //---       evaluations than smooth solver - about 2N times more. Using
      //---       numerical differentiation introduces additional (multiplicative)
      //---       2N speedup.
      //---       It means that if smooth optimizer WITH user-supplied gradient
      //---       needs 100 function evaluations to solve 50-dimensional problem,
      //---       then AGS solver with user-supplied gradient will need about 10.000
      //---       function evaluations, and with numerical gradient about 1.000.000
      //---       function evaluations will be performed.
      //--- NOTE: AGS solver used by us can handle nonsmooth and nonconvex
      //---       optimization problems. It has convergence guarantees, i.e. it will
      //---       converge to stationary point of the function after running for some
      //---       time.
      //---       However, it is important to remember that "stationary point" is not
      //---       equal to "solution". If your problem is convex, everything is OK.
      //---       But nonconvex optimization problems may have "flat spots" - large
      //---       areas where gradient is exactly zero, but function value is far away
      //---       from optimal. Such areas are stationary points too, and optimizer
      //---       may be trapped here.
      //---       "Flat spots" are nonsmooth equivalent of the saddle points, but with
      //---       orders of magnitude worse properties - they may be quite large and
      //---       hard to avoid. All nonsmooth optimizers are prone to this kind of the
      //---       problem, because it is impossible to automatically distinguish "flat
      //---       spot" from true solution.
      //---       This note is here to warn you that you should be very careful when
      //---       you solve nonsmooth optimization problems. Visual inspection of
      //---       results is essential.
      CAlglib::MinNSOptimize(state,Jac,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSResults(state,x1,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={1.0000,0.0000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x1,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonsmooth nonlinearly constrained optimization                   |
//+------------------------------------------------------------------+
void TEST_MinNS_D_NLC(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_NSFunc2_Jac Jac;
   CNDimensional_Rep Rep;
   CObject Obj;
   for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x0,x1) = 2*|x0|+|x1|
      //--- subject to combination of equality and inequality constraints
      //---      x0  =  1
      //---      x1 >= -1
      //--- using nonsmooth nonlinear optimizer. Although these constraints
      //--- are linear, we treat them as general nonlinear ones in order to
      //--- demonstrate nonlinearly constrained optimization setup.
      CRowDouble x0=vector<double>::Ones(2);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      double epsx=0.00001;
      if(_spoil_scenario==6)
         epsx=AL_NaN;
      if(_spoil_scenario==7)
         epsx=AL_POSINF;
      if(_spoil_scenario==8)
         epsx=AL_NEGINF;
      double radius=0.1;
      if(_spoil_scenario==9)
         radius=AL_NaN;
      if(_spoil_scenario==10)
         radius=AL_POSINF;
      if(_spoil_scenario==11)
         radius=AL_NEGINF;
      double rho=50.0;
      if(_spoil_scenario==12)
         rho=AL_NaN;
      if(_spoil_scenario==13)
         rho=AL_POSINF;
      if(_spoil_scenario==14)
         rho=AL_NEGINF;
      int maxits=0;
      CMinNSState state;
      CMinNSReport rep;
      CRowDouble x1;
      //--- Create optimizer object, choose AGS algorithm and tune its settings:
      //--- * radius=0.1     good initial value; will be automatically decreased later.
      //--- * rho=50.0       penalty coefficient for nonlinear constraints. It is your
      //---                  responsibility to choose good one - large enough that it
      //---                  enforces constraints, but small enough in order to avoid
      //---                  extreme slowdown due to ill-conditioning.
      //--- * epsx=0.000001  stopping conditions
      //--- * s=[1,1]        all variables have unit scale
      CAlglib::MinNSCreate(2,x0,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetAlgoAGS(state,radius,rho);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetCond(state,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Set general nonlinear constraints.
      //--- This part is more tricky than working with box/linear constraints - you
      //--- can not "pack" general nonlinear function into double precision array.
      //--- That's why minnssetnlc() does not accept constraints itself - only
      //--- constraint COUNTS are passed: first parameter is number of equality
      //--- constraints, second one is number of inequality constraints.
      //--- As for constraining functions - these functions are passed as part
      //--- of problem Jacobian (see below).
      //--- NOTE: MinNS optimizer supports arbitrary combination of boundary, general
      //---       linear and general nonlinear constraints. This example does not
      //---       show how to work with general linear constraints, but you can
      //---       easily find it in documentation on minnlcsetlc() function.
      CAlglib::MinNSSetNLC(state,1,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and test results.
      //--- Optimizer object accepts vector function and its Jacobian, with first
      //--- component (Jacobian row) being target function, and next components
      //--- (Jacobian rows) being nonlinear equality and inequality constraints
      //--- (box/linear ones are passed separately by means of minnssetbc() and
      //--- minnssetlc() calls).
      //--- Nonlinear equality constraints have form Gi(x)=0, inequality ones
      //--- have form Hi(x)<=0, so we may have to "normalize" constraints prior
      //--- to passing them to optimizer (right side is zero, constraints are
      //--- sorted, multiplied by -1 when needed).
      //--- So, our vector function has form
      //---     {f0,f1,f2} = { 2*|x0|+|x1|,  x0-1, -x1-1 }
      //--- with Jacobian
      //---         [ 2*sign(x0)   sign(x1) ]
      //---     J = [     1           0     ]
      //---         [     0          -1     ]
      //--- which means that we have optimization problem
      //---     min{f0} subject to f1=0, f2<=0
      //--- which is essentially same as
      //---     min { 2*|x0|+|x1| } subject to x0=1, x1>=-1
      //--- NOTE: AGS solver used by us can handle nonsmooth and nonconvex
      //---       optimization problems. It has convergence guarantees, i.e. it will
      //---       converge to stationary point of the function after running for some
      //---       time.
      //---       However, it is important to remember that "stationary point" is not
      //---       equal to "solution". If your problem is convex, everything is OK.
      //---       But nonconvex optimization problems may have "flat spots" - large
      //---       areas where gradient is exactly zero, but function value is far away
      //---       from optimal. Such areas are stationary points too, and optimizer
      //---       may be trapped here.
      //---       "Flat spots" are nonsmooth equivalent of the saddle points, but with
      //---       orders of magnitude worse properties - they may be quite large and
      //---       hard to avoid. All nonsmooth optimizers are prone to this kind of the
      //---       problem, because it is impossible to automatically distinguish "flat
      //---       spot" from true solution.
      //---       This note is here to warn you that you should be very careful when
      //---       you solve nonsmooth optimization problems. Visual inspection of
      //---       results is essential.
      CAlglib::MinNSOptimize(state,Jac,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinNSResults(state,x1,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={1.0000,0.0000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x1,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Nonlinear optimization with box constraints                      |
//+------------------------------------------------------------------+
void TEST_MinBC_D_1(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_Grad1 Grad;
   CNDimensional_Rep Rep;
   CObject Obj;
   for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x,y) = 100*(x+3)^4+(y-3)^4
      //--- subject to box constraints
      //---     -1<=x<=+1, -1<=y<=+1
      //--- using MinBC optimizer with:
      //--- * initial point x=[0,0]
      //--- * unit scale being set for all variables (see minbcsetscale for more info)
      //---*stopping criteria set to "terminate after short enough step"
      //--- * OptGuard integrity check being used to check problem statement
      //---   for some common errors like nonsmoothness or bad analytic gradient
      //--- First, we create optimizer object and tune its properties:
      //--- * set box constraints
      //--- * set variable scales
      //--- * set stopping criteria
      CRowDouble x=vector<double>::Zeros(2);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Deleting_Element(s);
      CRowDouble bndl=vector<double>::Full(2,-1);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(bndl,AL_NaN);
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(bndl);
      CRowDouble bndu=vector<double>::Ones(2);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(bndu,AL_NaN);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Deleting_Element(bndu);
      CMinBCState state;
      double epsg=0;
      if(_spoil_scenario==11)
         epsg=AL_NaN;
      if(_spoil_scenario==12)
         epsg=AL_POSINF;
      if(_spoil_scenario==13)
         epsg=AL_NEGINF;
      double epsf=0;
      if(_spoil_scenario==14)
         epsf=AL_NaN;
      if(_spoil_scenario==15)
         epsf=AL_POSINF;
      if(_spoil_scenario==16)
         epsf=AL_NEGINF;
      double epsx=0.000001;
      if(_spoil_scenario==17)
         epsx=AL_NaN;
      if(_spoil_scenario==18)
         epsx=AL_POSINF;
      if(_spoil_scenario==19)
         epsx=AL_NEGINF;
      int maxits=0;
      CAlglib::MinBCCreate(x,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Then we activate OptGuard integrity checking.
      //--- OptGuard monitor helps to catch common coding and problem statement
      //--- issues, like:
      //--- * discontinuity of the target function (C0 continuity violation)
      //--- * nonsmoothness of the target function (C1 continuity violation)
      //--- * erroneous analytic gradient, i.e. one inconsistent with actual
      //---   change in the target/constraints
      //--- OptGuard is essential for early prototyping stages because such
      //--- problems often result in premature termination of the optimizer
      //--- which is really hard to distinguish from the correct termination.
      //--- IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
      //---            DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
      //---            Other OptGuard checks add moderate overhead, but anyway
      //---            it is better to turn them off when they are not needed.
      CAlglib::MinBCOptGuardSmoothness(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCOptGuardGradient(state,0.001);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and evaluate results
      CMinBCReport rep;
      CAlglib::MinBCOptimize(state,Grad,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={-1,1};
         _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.005);
        }
      //--- Check that OptGuard did not report errors
      //--- NOTE: want to test OptGuard? Try breaking the gradient - say, add
      //---       1.0 to some of its components.
      COptGuardReport ogrep;
      CAlglib::MinBCOptGuardResults(state,ogrep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_badgradsuspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc0suspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc1suspected,false);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void TEST_MinBC_NumDif(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   CNDimensional_Func1 Func;
   CNDimensional_Rep Rep;
   CObject Obj;

   for(_spoil_scenario=-1; _spoil_scenario<23; _spoil_scenario++)
     {
      //--- This example demonstrates minimization of
      //---     f(x,y) = 100*(x+3)^4+(y-3)^4
      //--- subject to box constraints
      //---    -1<=x<=+1, -1<=y<=+1
      //--- using MinBC optimizer with:
      //--- * numerical differentiation being used
      //--- * initial point x=[0,0]
      //--- * unit scale being set for all variables (see minbcsetscale for more info)
      //---*stopping criteria set to "terminate after short enough step"
      //--- * OptGuard integrity check being used to check problem statement
      //---   for some common errors like nonsmoothness or bad analytic gradient
      CRowDouble x=vector<double>::Zeros(2);
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      CRowDouble s=vector<double>::Ones(2);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Value(s,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(s,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(s,AL_NEGINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Deleting_Element(s);
      CRowDouble bndl=vector<double>::Full(2,-1);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(bndl,AL_NaN);
      if(_spoil_scenario==8)
         Spoil_Vector_By_Deleting_Element(bndl);
      CRowDouble bndu=vector<double>::Ones(2);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(bndu,AL_NaN);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Deleting_Element(bndu);
      CMinBCState state;
      double epsg=0;
      if(_spoil_scenario==11)
         epsg=AL_NaN;
      if(_spoil_scenario==12)
         epsg=AL_POSINF;
      if(_spoil_scenario==13)
         epsg=AL_NEGINF;
      double epsf=0;
      if(_spoil_scenario==14)
         epsf=AL_NaN;
      if(_spoil_scenario==15)
         epsf=AL_POSINF;
      if(_spoil_scenario==16)
         epsf=AL_NEGINF;
      double epsx=0.000001;
      if(_spoil_scenario==17)
         epsx=AL_NaN;
      if(_spoil_scenario==18)
         epsx=AL_POSINF;
      if(_spoil_scenario==19)
         epsx=AL_NEGINF;
      int maxits=0;
      double diffstep=1.0e-6;
      if(_spoil_scenario==20)
         diffstep=AL_NaN;
      if(_spoil_scenario==21)
         diffstep=AL_POSINF;
      if(_spoil_scenario==22)
         diffstep=AL_NEGINF;
      //--- Now we are ready to actually optimize something:
      //--- * first we create optimizer
      //--- * we add boundary constraints
      //--- * we tune stopping conditions
      //--- * and, finally, optimize and obtain results...
      CAlglib::MinBCCreateF(x,diffstep,state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCSetBC(state,bndl,bndu);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCSetScale(state,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCSetCond(state,epsg,epsf,epsx,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Then we activate OptGuard integrity checking.
      //--- Numerical differentiation always produces "correct" gradient
      //--- (with some truncation error, but unbiased). Thus, we just have
      //--- to check smoothness properties of the target: C0 and C1 continuity.
      //--- Sometimes user accidentally tries to solve nonsmooth problems
      //--- with smooth optimizer. OptGuard helps to detect such situations
      //--- early, at the prototyping stage.
      CAlglib::MinBCOptGuardSmoothness(state);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Optimize and evaluate results
      CMinBCReport rep;
      CAlglib::MinBCOptimize(state,Func,Rep,Obj);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCResults(state,x,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={-1,1};
         _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.005);
        }
      //--- Check that OptGuard did not report errors
      //--- Want to challenge OptGuard? Try to make your problem
      //--- nonsmooth by replacing 100*(x+3)^4 by 100*|x+3| and
      //--- real-run optimizer.
      COptGuardReport ogrep;
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MinBCOptGuardResults(state,ogrep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc0suspected,false);
      _TestResult=_TestResult && Doc_Test_Bool(ogrep.m_nonc1suspected,false);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple model built with IDW-MSTAB algorithm                      |
//+------------------------------------------------------------------+
void TEST_IDW_D_MSTAB(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- This example illustrates basic concepts of the IDW models:
      //--- creation and evaluation.
      //--- Suppose that we have set of 2-dimensional points with associated
      //--- scalar function values, and we want to build an IDW model using
      //--- our data.
      //--- NOTE: we can work with N-dimensional models and vector-valued functions too :)
      //--- Typical sequence of steps is given below:
      //--- 1. we create IDW builder object
      //--- 2. we attach our dataset to the IDW builder and tune algorithm settings
      //--- 3. we generate IDW model
      //--- 4. we use IDW model instance (evaluate, serialize, etc.)
      double v;
      //--- Step 1: IDW builder creation.
      //--- We have to specify dimensionality of the space (2 or 3) and
      //--- dimensionality of the function (scalar or vector).
      //--- New builder object is empty - it has not dataset and uses
      //--- default model construction settings
      CIDWBuilder builder;
      CAlglib::IDWBuilderCreate(2,1,builder);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Step 2: dataset addition
      //--- XY contains two points - x0=(-1,0) and x1=(+1,0) -
      //--- and two function values f(x0)=2, f(x1)=3.
      matrix<double> XY={{-1,0,2},{+1,0,3}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::IDWBuilderSetPoints(builder,xy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Step 3: choose IDW algorithm and generate model
      //--- We use modified stabilized IDW algorithm with following parameters:
      //--- * SRad - set to 5.0 (search radius must be large enough)
      //--- IDW-MSTAB algorithm is a state-of-the-art implementation of IDW which
      //--- is competitive with RBFs and bicubic splines. See comments on the
      //--- idwbuildersetalgomstab() function for more information.
      CIDWModelShell model;
      CIDWReport rep;
      CAlglib::IDWBuilderSetAlgoMSTAB(builder,5.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::IDWFit(builder,model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Step 4: model was built, evaluate its value
      v=CAlglib::IDWCalc2(model,1.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,3.000,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| IDW model serialization/unserialization                          |
//+------------------------------------------------------------------+
void TEST_IDW_D_Serialize(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- This example shows how to serialize and unserialize IDW model.
      //--- Suppose that we have set of 2-dimensional points with associated
      //--- scalar function values, and we have built an IDW model using
      //--- our data.
      //--- This model can be serialized to string or stream. ALGLIB supports
      //--- flexible (un)serialization, i.e. you can move serialized model
      //--- representation between different machines (32-bit or 64-bit),
      //--- different CPU architectures (x86/64, ARM) or even different
      //--- programming languages supported by ALGLIB (C#, C++, ...).
      //--- Our first step is to build model, evaluate it at point (1,0),
      //--- and serialize it to string.
      string s;
      double v;
      matrix<double> XY={{-1,0,2},{+1,0,3}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CIDWBuilder builder;
      CIDWModelShell model;
      CIDWModelShell model2;
      CIDWReport rep;
      CAlglib::IDWBuilderCreate(2,1,builder);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::IDWBuilderSetPoints(builder,xy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::IDWBuilderSetAlgoMSTAB(builder,5.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::IDWFit(builder,model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      v=CAlglib::IDWCalc2(model,1.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,3.000,0.005);
      //--- Serialization + unserialization to a different instance
      //--- of the model class.
      CAlglib::IDWSerialize(model,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::IDWUnserialize(s,model2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Evaluate unserialized model at the same point
      v=CAlglib::IDWCalc2(model2,1.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,3.000,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Parametric Ramer-Douglas-Peucker approximation                   |
//+------------------------------------------------------------------+
void TEST_Parametric_RDP(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
     {
      //--- We use RDP algorithm to approximate parametric 2D curve given by
      //--- locations in t=0,1,2,3 (see below), which form piecewise linear
      //--- trajectory through D-dimensional space (2-dimensional in our example).
      //---     |
      //---     |
      //---     -     *     *     X2................X3
      //---     |                .
      //---     |               .
      //---     -     *     *  .  *     *     *     *
      //---     |             .
      //---     |            .
      //---     -     *     X1    *     *     *     *
      //---     |      .....
      //---     |  ....
      //---     X0----|-----|-----|-----|-----|-----|---
      int npoints=4;
      int ndimensions=2;
      matrix<double> X={{0,0},{2,1},{3,3},{6,3}};
      CMatrixDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(x);
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(x);
      //--- Approximation of parametric curve is performed by another parametric curve
      //--- with lesser amount of points. It allows to work with "compressed"
      //--- representation, which needs smaller amount of memory. Say, in our example
      //--- (we allow points with error smaller than 0.8) approximation will have
      //--- just two sequential sections connecting X0 with X2, and X2 with X3.
      //---     |
      //---     |
      //---     -     *     *     X2................X3
      //---     |               .
      //---     |             .
      //---     -     *     .     *     *     *     *
      //---     |         .
      //---     |       .
      //---     -     .     X1    *     *     *     *
      //---     |   .
      //---     | .
      //---     X0----|-----|-----|-----|-----|-----|---
      CMatrixDouble y;
      int idxy[];
      int nsections;
      int limitcnt=0;
      double limiteps=0.8;
      if(_spoil_scenario==5)
         limiteps=AL_POSINF;
      if(_spoil_scenario==6)
         limiteps=AL_NEGINF;
      CAlglib::ParametricRDPFixed(x,npoints,ndimensions,limitcnt,limiteps,y,idxy,nsections);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      int check[]={0,2,3};
      _TestResult=_TestResult && Doc_Test_Int(nsections,2);
      _TestResult=_TestResult && Doc_Test_Int_Vector(idxy,check);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Bilinear spline interpolation                                    |
//+------------------------------------------------------------------+
void TEST_Spline2D_Bilinear(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<16; _spoil_scenario++)
     {
      //--- We use bilinear spline to interpolate f(x,y)=x^2+2*y^2 sampled
      //--- at (x,y) from [0.0, 0.5, 1.0] X [0.0, 1.0].
      double X[]={0.0,0.5,1.0};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      double Y[]={0.0,1.0};
      CRowDouble y=Y;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      double F[]={0.00,0.25,1.00,2.00,2.25,3.00};
      CRowDouble f=F;
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(f,AL_NaN);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(f,AL_POSINF);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(f,AL_NEGINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(f);
      double vx=0.25;
      if(_spoil_scenario==12)
         vx=AL_POSINF;
      if(_spoil_scenario==13)
         vx=AL_NEGINF;
      double vy=0.50;
      if(_spoil_scenario==14)
         vy=AL_POSINF;
      if(_spoil_scenario==15)
         vy=AL_NEGINF;
      double v;
      CSpline2DInterpolantShell s;
      //--- build spline
      CAlglib::Spline2DBuildBilinearV(x,3,y,2,f,1,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- calculate S(0.25,0.50)
      v=CAlglib::Spline2DCalc(s,vx,vy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,1.1250,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Bilinear spline interpolation                                    |
//+------------------------------------------------------------------+
void TEST_Spline2D_Bicubic(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<16; _spoil_scenario++)
     {
      //--- We use bilinear spline to interpolate f(x,y)=x^2+2*y^2 sampled
      //--- at (x,y) from [0.0, 0.5, 1.0] X [0.0, 1.0].
      double X[]={0.0,0.5,1.0};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      double Y[]={0.0,1.0};
      CRowDouble y=Y;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      double F[]={0.00,0.25,1.00,2.00,2.25,3.00};
      CRowDouble f=F;
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(f,AL_NaN);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(f,AL_POSINF);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(f,AL_NEGINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(f);
      double vx=0.25;
      if(_spoil_scenario==12)
         vx=AL_POSINF;
      if(_spoil_scenario==13)
         vx=AL_NEGINF;
      double vy=0.50;
      if(_spoil_scenario==14)
         vy=AL_POSINF;
      if(_spoil_scenario==15)
         vy=AL_NEGINF;
      double v;
      double dx;
      double dy;
      double dxy;
      CSpline2DInterpolantShell s;
      //--- build spline
      CAlglib::Spline2DBuildBicubicV(x,3,y,2,f,1,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- calculate S(0.25,0.50)
      v=CAlglib::Spline2DCalc(s,vx,vy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,1.0625,0.00005);
      //--- calculate derivatives
      CAlglib::Spline2DDiff(s,vx,vy,v,dx,dy,dxy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,1.0625,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(dx,0.5000,0.00005);
      _TestResult=_TestResult && Doc_Test_Real(dy,2.0000,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Fitting bicubic spline to irregular data                         |
//+------------------------------------------------------------------+
void TEST_Spline2D_Fit_Blocklls(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
     {
      //--- We use bicubic spline to reproduce f(x,y)=1/(1+x^2+2*y^2) sampled
      //--- at irregular points (x,y) from [-1,+1]*[-1,+1]
      //--- We have 5 such points, located approximately at corners of the area
      //--- and its center -  but not exactly at the grid. Thus, we have to FIT
      //--- the spline, i.e. to solve least squares problem
      matrix<double> XY={{-0.987,-0.902,0.359},{0.948,-0.992,0.347},{-1.000,1.000,0.333},{1.000,0.973,0.339},{0.017,0.180,0.968}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Deleting_Row(xy);
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Deleting_Col(xy);
      //--- First step is to create spline2dbuilder object and set its properties:
      //--- * d=1 means that we create vector-valued spline with 1 component
      //--- * we specify dataset xy
      //--- * we rely on automatic selection of interpolation area
      //--- * we tell builder that we want to use 5x5 grid for an underlying spline
      //--- * we choose least squares solver named BlockLLS and configure it by
      //---   telling that we want to apply zero nonlinearity penalty.
      //--- NOTE: you can specify non-zero lambdav if you want to make your spline
      //---       more "rigid", i.e. to penalize nonlinearity.
      //--- NOTE: ALGLIB has two solvers which fit bicubic splines to irregular data,
      //---       one of them is BlockLLS and another one is FastDDM. Former is
      //---       intended for moderately sized grids (up to 512x512 nodes, although
      //---       it may take up to few minutes); it is the most easy to use and
      //---       control spline fitting function in the library. Latter, FastDDM,
      //---       is intended for efficient solution of large-scale problems
      //---       (up to 100.000.000 nodes). Both solvers can be parallelized, but
      //---       FastDDM is much more efficient. See comments for more information.
      CSpline2DBuilder builder;
      int d=1;
      double lambdav=0.000;
      CAlglib::Spline2DBuilderCreate(d,builder);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DBuilderSetPoints(builder,xy,5);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DBuilderSetGrid(builder,5,5);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DBuilderSetAlgoBlockLLS(builder,lambdav);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Now we are ready to fit and evaluate our results
      CSpline2DInterpolantShell s;
      CSpline2DFitReport rep;
      CAlglib::Spline2DFit(builder,s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- evaluate results - function value at the grid is reproduced exactly
      double v;
      v=CAlglib::Spline2DCalc(s,-1,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,0.333000,0.005);
      //--- check maximum error - it must be nearly zero
      _TestResult=_TestResult && Doc_Test_Real(rep.m_maxerror,0.000,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Unpacking bilinear spline                                        |
//+------------------------------------------------------------------+
void TEST_Spline2D_Unpack(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- We build bilinear spline for f(x,y)=x+2*y+3*xy for (x,y) in [0,1].
      //--- Then we demonstrate how to unpack it.
      double X[]={0.0,1.0};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      double Y[]={0.0,1.0};
      CRowDouble y=Y;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      double F[]={0.00,1.00,2.00,6.00};
      CRowDouble f=F;
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(f,AL_NaN);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(f,AL_POSINF);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(f,AL_NEGINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(f);
      CMatrixDouble c;
      int m;
      int n;
      int d;
      CSpline2DInterpolantShell s;
      //--- build spline
      CAlglib::Spline2DBuildBilinearV(x,2,y,2,f,1,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- unpack and test
      CAlglib::Spline2DUnpackV(s,m,n,d,c);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      matrix<double> check={{0,1,0,1,0,2,0,0,1,3,0,0,0,0,0,0,0,0,0,0}};
      _TestResult=_TestResult && Doc_Test_Real_Matrix(c,check,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Copy and transform                                               |
//+------------------------------------------------------------------+
void TEST_Spline2D_CopyTrans(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<16; _spoil_scenario++)
     {
      //--- We build bilinear spline for f(x,y)=x+2*y for (x,y) in [0,1].
      //--- Then we apply several transformations to this spline.
      double X[]={0.0,1.0};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      double Y[]={0.0,1.0};
      CRowDouble y=Y;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      double F[]={0.00,1.00,2.00,3.00};
      CRowDouble f=F;
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(f,AL_NaN);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(f,AL_POSINF);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(f,AL_NEGINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(f);
      CSpline2DInterpolantShell s;
      CSpline2DInterpolantShell snew;
      double v;
      CAlglib::Spline2DBuildBilinearV(x,2,y,2,f,1,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- copy spline, apply transformation x:=2*xnew, y:=4*ynew
      //--- evaluate at (xnew,ynew) = (0.25,0.25) - should be same as (x,y)=(0.5,1.0)
      CAlglib::Spline2DCopy(s,snew);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DLinTransXY(snew,2.0,0.0,4.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      v=CAlglib::Spline2DCalc(snew,0.25,0.25);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,2.500,0.00005);
      //--- copy spline, apply transformation SNew:=2*S+3
      CAlglib::Spline2DCopy(s,snew);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DLinTransF(snew,2.0,3.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      v=CAlglib::Spline2DCalc(snew,0.5,1.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,8.000,0.00005);
      //--- Same example, but for vector spline (f0,f1) = {x+2*y, 2*x+y}
      double F2[]={0.00,0.00,1.00,2.00,2.00,1.00,3.00,3.00};
      CRowDouble f2=F2;
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(f2,AL_NaN);
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(f2,AL_POSINF);
      if(_spoil_scenario==14)
         Spoil_Vector_By_Value(f2,AL_NEGINF);
      if(_spoil_scenario==15)
         Spoil_Vector_By_Deleting_Element(f2);
      CRowDouble vr;
      CAlglib::Spline2DBuildBilinearV(x,2,y,2,f2,2,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- copy spline, apply transformation x:=2*xnew, y:=4*ynew
      CAlglib::Spline2DCopy(s,snew);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DLinTransXY(snew,2.0,0.0,4.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DCalcV(snew,0.25,0.25,vr);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={2.500,2.000};
         _TestResult=_TestResult && Doc_Test_Real_Vector(vr,check,0.00005);
        }
      //--- copy spline, apply transformation SNew:=2*S+3
      CAlglib::Spline2DCopy(s,snew);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DLinTransF(snew,2.0,3.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DCalcV(snew,0.5,1.0,vr);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={8.000,7.000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(vr,check,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Copy and transform                                               |
//+------------------------------------------------------------------+
void TEST_Spline2d_Vector(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
     {
      //--- We build bilinear vector-valued spline (f0,f1) = {x+2*y, 2*x+y}
      //--- Spline is built using function values at 2x2 grid: (x,y)=[0,1]*[0,1]
      //--- Then we perform evaluation at (x,y)=(0.1,0.3)
      double X[]={0.0,1.0};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      double Y[]={0.0,1.0};
      CRowDouble y=Y;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      double F[]={0.00,0.00,1.00,2.00,2.00,1.00,3.00,3.00};
      CRowDouble f=F;
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(f,AL_NaN);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(f,AL_POSINF);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(f,AL_NEGINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(f);
      CSpline2DInterpolantShell s;
      CRowDouble vr;
      CAlglib::Spline2DBuildBilinearV(x,2,y,2,f,2,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::Spline2DCalcV(s,0.1,0.3,vr);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={0.700,0.500};
      _TestResult=_TestResult && Doc_Test_Real_Vector(vr,check,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Trilinear spline interpolation                                   |
//+------------------------------------------------------------------+
void TEST_Spline3D_Trilinear(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<22; _spoil_scenario++)
     {
      //--- We use trilinear spline to interpolate f(x,y,z)=x+xy+z sampled
      //--- at (x,y,z) from [0.0, 1.0] X [0.0, 1.0] X [0.0, 1.0].
      //--- We store x, y and z-values at local arrays with same names.
      //--- Function values are stored in the array F as follows:
      //---     f[0]     (x,y,z) = (0,0,0)
      //---     f[1]     (x,y,z) = (1,0,0)
      //---     f[2]     (x,y,z) = (0,1,0)
      //---     f[3]     (x,y,z) = (1,1,0)
      //---     f[4]     (x,y,z) = (0,0,1)
      //---     f[5]     (x,y,z) = (1,0,1)
      //---     f[6]     (x,y,z) = (0,1,1)
      //---     f[7]     (x,y,z) = (1,1,1)
      double X[]={0.0,1.0};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      CRowDouble y=X;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      CRowDouble z=X;
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(z,AL_NaN);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(z,AL_POSINF);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(z,AL_NEGINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(z);
      double F[]={0,1,0,2,1,2,1,3};
      CRowDouble f=F;
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(f,AL_NaN);
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(f,AL_POSINF);
      if(_spoil_scenario==14)
         Spoil_Vector_By_Value(f,AL_NEGINF);
      if(_spoil_scenario==15)
         Spoil_Vector_By_Deleting_Element(f);
      double vx=0.50;
      if(_spoil_scenario==16)
         vx=AL_POSINF;
      if(_spoil_scenario==17)
         vx=AL_NEGINF;
      double vy=0.50;
      if(_spoil_scenario==18)
         vy=AL_POSINF;
      if(_spoil_scenario==19)
         vy=AL_NEGINF;
      double vz=0.50;
      if(_spoil_scenario==20)
         vz=AL_POSINF;
      if(_spoil_scenario==21)
         vz=AL_NEGINF;
      double v;
      CSpline3DInterpolant s;
      //--- build spline
      CAlglib::Spline3DBuildTrilinearV(x,2,y,2,z,2,f,1,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- calculate S(0.5,0.5,0.5)
      v=CAlglib::Spline3DCalc(s,vx,vy,vz);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,1.2500,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Vector-valued trilinear spline interpolation                     |
//+------------------------------------------------------------------+
void TEST_Spline3D_Vector(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<22; _spoil_scenario++)
     {
      //--- We use trilinear vector-valued spline to interpolate {f0,f1}={x+xy+z,x+xy+yz+z}
      //--- sampled at (x,y,z) from [0.0, 1.0] X [0.0, 1.0] X [0.0, 1.0].
      //--- We store x, y and z-values at local arrays with same names.
      //--- Function values are stored in the array F as follows:
      //---     f[0]     f0, (x,y,z) = (0,0,0)
      //---     f[1]     f1, (x,y,z) = (0,0,0)
      //---     f[2]     f0, (x,y,z) = (1,0,0)
      //---     f[3]     f1, (x,y,z) = (1,0,0)
      //---     f[4]     f0, (x,y,z) = (0,1,0)
      //---     f[5]     f1, (x,y,z) = (0,1,0)
      //---     f[6]     f0, (x,y,z) = (1,1,0)
      //---     f[7]     f1, (x,y,z) = (1,1,0)
      //---     f[8]     f0, (x,y,z) = (0,0,1)
      //---     f[9]     f1, (x,y,z) = (0,0,1)
      //---     f[10]    f0, (x,y,z) = (1,0,1)
      //---     f[11]    f1, (x,y,z) = (1,0,1)
      //---     f[12]    f0, (x,y,z) = (0,1,1)
      //---     f[13]    f1, (x,y,z) = (0,1,1)
      //---     f[14]    f0, (x,y,z) = (1,1,1)
      //---     f[15]    f1, (x,y,z) = (1,1,1)
      double X[]={0.0,1.0};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      CRowDouble y=X;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      CRowDouble z=X;
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(z,AL_NaN);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Value(z,AL_POSINF);
      if(_spoil_scenario==10)
         Spoil_Vector_By_Value(z,AL_NEGINF);
      if(_spoil_scenario==11)
         Spoil_Vector_By_Deleting_Element(z);
      double F[]={0,0,1,1,0,0,2,2,1,1,2,2,1,2,3,4};
      CRowDouble f=F;
      if(_spoil_scenario==12)
         Spoil_Vector_By_Value(f,AL_NaN);
      if(_spoil_scenario==13)
         Spoil_Vector_By_Value(f,AL_POSINF);
      if(_spoil_scenario==14)
         Spoil_Vector_By_Value(f,AL_NEGINF);
      if(_spoil_scenario==15)
         Spoil_Vector_By_Deleting_Element(f);
      double vx=0.50;
      if(_spoil_scenario==16)
         vx=AL_POSINF;
      if(_spoil_scenario==17)
         vx=AL_NEGINF;
      double vy=0.50;
      if(_spoil_scenario==18)
         vy=AL_POSINF;
      if(_spoil_scenario==19)
         vy=AL_NEGINF;
      double vz=0.50;
      if(_spoil_scenario==20)
         vz=AL_POSINF;
      if(_spoil_scenario==21)
         vz=AL_NEGINF;
      CSpline3DInterpolant s;
      //--- build spline
      CAlglib::Spline3DBuildTrilinearV(x,2,y,2,z,2,f,2,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- calculate S(0.5,0.5,0.5) - we have vector of values instead of single value
      CRowDouble v;
      CAlglib::Spline3DCalcV(s,vx,vy,vz,v);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={1.2500,1.5000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(v,check,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void TEST_RBF_D_HRBF(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- This example illustrates basic concepts of the RBF models: creation, modification,
      //--- evaluation.
      //--- Suppose that we have set of 2-dimensional points with associated
      //--- scalar function values, and we want to build a RBF model using
      //--- our data.
      //--- NOTE: we can work with 3D models too :)
      //--- Typical sequence of steps is given below:
      //--- 1. we create RBF model object
      //--- 2. we attach our dataset to the RBF model and tune algorithm settings
      //--- 3. we rebuild RBF model using QNN algorithm on new data
      //--- 4. we use RBF model (evaluate, serialize, etc.)
      double v;
      //--- Step 1: RBF model creation.
      //--- We have to specify dimensionality of the space (2 or 3) and
      //--- dimensionality of the function (scalar or vector).
      //--- New model is empty - it can be evaluated,
      //--- but we just get zero value at any point.
      CRBFModel model;
      CAlglib::RBFCreate(2,1,model);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      v=CAlglib::RBFCalc2(model,0.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,0.000,0.005);
      //--- Step 2: we add dataset.
      //--- XY contains two points - x0=(-1,0) and x1=(+1,0) -
      //--- and two function values f(x0)=2, f(x1)=3.
      //--- We added points, but model was not rebuild yet.
      //--- If we call RBFCalc2(), we still will get 0.0 as result.
      matrix<double> XY={{-1,0,2},{+1,0,3}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::RBFSetPoints(model,xy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      v=CAlglib::RBFCalc2(model,0.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,0.000,0.005);
      //--- Step 3: rebuild model
      //--- After we've configured model, we should rebuild it -
      //--- it will change coefficients stored internally in the
      //--- rbfmodel structure.
      //--- We use hierarchical RBF algorithm with following parameters:
      //--- * RBase - set to 1.0
      //--- * NLayers - three layers are used (although such simple problem
      //---   does not need more than 1 layer)
      //--- * LambdaReg - is set to zero value, no smoothing is required
      CRBFReport rep;
      CAlglib::RBFSetAlgoHierarchical(model,1.0,3,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFBuildModel(model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,1);
      //--- Step 4: model was built
      //--- After call of RBFBuildModel(), RBFCalc2() will return
      //--- value of the new model.
      v=CAlglib::RBFCalc2(model,0.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,2.500,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Working with vector functions                                    |
//+------------------------------------------------------------------+
void TEST_RBF_D_Vector(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- Suppose that we have set of 2-dimensional points with associated VECTOR
      //--- function values, and we want to build a RBF model using our data.
      //--- Typical sequence of steps is given below:
      //--- 1. we create RBF model object
      //--- 2. we attach our dataset to the RBF model and tune algorithm settings
      //--- 3. we rebuild RBF model using new data
      //--- 4. we use RBF model (evaluate, serialize, etc.)
      CRowDouble x;
      CRowDouble y;
      //--- Step 1: RBF model creation.
      //--- We have to specify dimensionality of the space (equal to 2) and
      //--- dimensionality of the function (2-dimensional vector function).
      //--- New model is empty - it can be evaluated,
      //--- but we just get zero value at any point.
      CRBFModel model;
      CAlglib::RBFCreate(2,2,model);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double X[]={+1,+1};
      x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      CAlglib::RBFCalc(model,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={0.000,0.000};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.005);
        }
      //--- Step 2: we add dataset.
      //--- XY arrays containt four points:
      //--- * (x0,y0) = (+1,+1), f(x0,y0)=(0,-1)
      //--- * (x1,y1) = (+1,-1), f(x1,y1)=(-1,0)
      //--- * (x2,y2) = (-1,-1), f(x2,y2)=(0,+1)
      //--- * (x3,y3) = (-1,+1), f(x3,y3)=(+1,0)
      matrix<double> XY={{+1,+1,0,-1},{+1,-1,-1,0},{-1,-1,0,+1},{-1,+1,+1,0}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==3)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==4)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==5)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::RBFSetPoints(model,xy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- We added points, but model was not rebuild yet.
      //--- If we call RBFCalc(), we still will get 0.0 as result.
      CAlglib::RBFCalc(model,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real_Vector(y,vector<double>::Zeros(2),0.005);
      //--- Step 3: rebuild model
      //--- We use hierarchical RBF algorithm with following parameters:
      //--- * RBase - set to 1.0
      //--- * NLayers - three layers are used (although such simple problem
      //---   does not need more than 1 layer)
      //--- * LambdaReg - is set to zero value, no smoothing is required
      //--- After we've configured model, we should rebuild it -
      //--- it will change coefficients stored internally in the
      //--- rbfmodel structure.
      CRBFReport rep;
      CAlglib::RBFSetAlgoHierarchical(model,1.0,3,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFBuildModel(model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,1);
      //--- Step 4: model was built
      //--- After call of RBFBuildModel(), RBFCalc() will return
      //--- value of the new model.
      CAlglib::RBFCalc(model,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={0.000,-1.000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//|RBF models - working with polynomial term                         |
//+------------------------------------------------------------------+
void TEST_RBF_D_PolTerm(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- This example show how to work with polynomial term
      //--- Suppose that we have set of 2-dimensional points with associated
      //--- scalar function values, and we want to build a RBF model using
      //--- our data.
      //--- We use hierarchical RBF algorithm with following parameters:
      //--- * RBase - set to 1.0
      //--- * NLayers - three layers are used (although such simple problem
      //---   does not need more than 1 layer)
      //--- * LambdaReg - is set to zero value, no smoothing is required
      double v;
      CRBFModel model;
      matrix<double> XY={{-1,0,2},{+1,0,3}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CRBFReport rep;
      CAlglib::RBFCreate(2,1,model);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFSetPoints(model,xy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFSetAlgoHierarchical(model,1.0,3,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- By default, RBF model uses linear term. It means that model
      //--- looks like
      //---     f(x,y) = SUM(RBF[i]) + a*x + b*y + c
      //--- where RBF[i] is I-th radial basis function and a*x+by+c is a
      //--- linear term. Having linear terms in a model gives us:
      //--- (1) improved extrapolation properties
      //--- (2) linearity of the model when data can be perfectly fitted
      //---     by the linear function
      //--- (3) linear asymptotic behavior
      //--- Our simple dataset can be modelled by the linear function
      //---     f(x,y) = 0.5*x + 2.5
      //--- and RBFBuildModel() with default settings should preserve this
      //--- linearity.
      int nx;
      int ny;
      int nc;
      int modelversion;
      CMatrixDouble xwr;
      CMatrixDouble c;
      CAlglib::RBFBuildModel(model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,1);
      CAlglib::RBFUnpack(model,nx,ny,xwr,nc,c,modelversion);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         matrix<double> check={{0.500,0.000,2.500}};
         _TestResult=_TestResult && Doc_Test_Real_Matrix(c,check,0.005);
        }
      //--- asymptotic behavior of our function is linear
      v=CAlglib::RBFCalc2(model,1000.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,502.50,0.05);
      //--- Instead of linear term we can use constant term. In this case
      //--- we will get model which has form
      //---     f(x,y) = SUM(RBF[i]) + c
      //--- where RBF[i] is I-th radial basis function and c is a constant,
      //--- which is equal to the average function value on the dataset.
      //--- Because we've already attached dataset to the model the only
      //--- thing we have to do is to call rbfsetconstterm() and then
      //--- rebuild model with RBFBuildModel().
      CAlglib::RBFSetConstTerm(model);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFBuildModel(model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,1);
      CAlglib::RBFUnpack(model,nx,ny,xwr,nc,c,modelversion);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         matrix<double> check={{0.000,0.000,2.500}};
         _TestResult=_TestResult && Doc_Test_Real_Matrix(c,check,0.005);
        }     //--- asymptotic behavior of our function is constant
      v=CAlglib::RBFCalc2(model,1000.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,2.500,0.005);
      //--- Finally, we can use zero term. Just plain RBF without polynomial
      //--- part:
      //---     f(x,y) = SUM(RBF[i])
      //--- where RBF[i] is I-th radial basis function.
      CAlglib::RBFSetZeroTerm(model);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFBuildModel(model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,1);
      CAlglib::RBFUnpack(model,nx,ny,xwr,nc,c,modelversion);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         matrix<double> check=matrix<double>::Zeros(1,3);
         _TestResult=_TestResult && Doc_Test_Real_Matrix(c,check,0.005);
        }
      //--- asymptotic behavior of our function is just zero constant
      v=CAlglib::RBFCalc2(model,1000.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,0.000,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Serialization/unserialization                                    |
//+------------------------------------------------------------------+
void TEST_RBF_D_Serialize(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- This example show how to serialize and unserialize RBF model
      //--- Suppose that we have set of 2-dimensional points with associated
      //--- scalar function values, and we want to build a RBF model using
      //--- our data. Then we want to serialize it to string and to unserialize
      //--- from string, loading to another instance of RBF model.
      //--- Here we assume that you already know how to create RBF models.
      string s;
      double v;
      CRBFModel model0;
      CRBFModel model1;
      matrix<double> XY={{-1,0,2},{+1,0,3}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CRBFReport rep;
      //--- model initialization
      CAlglib::RBFCreate(2,1,model0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFSetPoints(model0,xy);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFSetAlgoHierarchical(model0,1.0,3,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFBuildModel(model0,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,1);
      //--- Serialization - it looks easy,
      //--- but you should carefully read next section.
      CAlglib::RBFSerialize(model0,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::RBFUnserialize(s,model1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- both models return same value
      v=CAlglib::RBFCalc2(model0,0.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,2.500,0.005);
      v=CAlglib::RBFCalc2(model1,0.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,2.500,0.005);
      //--- Previous section shows that model state is saved/restored during
      //--- serialization. However, some properties are NOT serialized.
      //--- Serialization saves/restores RBF model, but it does NOT saves/restores
      //--- settings which were used to build current model. In particular, dataset
      //--- which was used to build model, is not preserved.
      //--- What does it mean in for us?
      //--- Do you remember this sequence: RBFCreate-RBFSetPoints-RBFBuildModel?
      //--- First step creates model, second step adds dataset and tunes model
      //--- settings, third step builds model using current dataset and model
      //--- construction settings.
      //--- If you call RBFBuildModel() without calling RBFSetPoints() first, you
      //--- will get empty (zero) RBF model. In our example, model0 contains
      //--- dataset which was added by RBFSetPoints() call. However, model1 does
      //--- NOT contain dataset - because dataset is NOT serialized.
      //--- This, if we call RBFBuildModel(model0,rep), we will get same model,
      //--- which returns 2.5 at (x,y)=(0,0). However, after same call model1 will
      //--- return zero - because it contains RBF model (coefficients), but does NOT
      //--- contain dataset which was used to build this model.
      //--- Basically, it means that:
      //--- * serialization of the RBF model preserves anything related to the model
      //---   EVALUATION
      //--- * but it does NOT creates perfect copy of the original object.
      CAlglib::RBFBuildModel(model0,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      v=CAlglib::RBFCalc2(model0,0.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,2.500,0.005);
      CAlglib::RBFBuildModel(model1,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      v=CAlglib::RBFCalc2(model1,0.0,0.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,0.000,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple hierarchical clusterization with Euclidean distance       |
//| function                                                         |
//+------------------------------------------------------------------+
void TEST_Clst_AHC(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- The very simple clusterization example
      //--- We have a set of points in 2D space:
      //---     (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
      //---  |
      //---  |     P3
      //---  |
      //---  | P1
      //---  |             P4
      //---  | P0          P2
      //---  |-------------------------
      //--- We want to perform Agglomerative Hierarchic Clusterization (AHC),
      //--- using complete linkage (default algorithm) and Euclidean distance
      //--- (default metric).
      //--- In order to do that, we:
      //--- * create clusterizer with ClusterizerCreate()
      //--- * set points XY and metric (2=Euclidean) with ClusterizerSetPoints()
      //--- * run AHC algorithm with ClusterizerRunAHC
      //--- You may see that clusterization itself is a minor part of the example,
      //--- most of which is dominated by comments :)
      CClusterizerState s;
      CAHCReport rep;
      matrix<double> XY={{1,1},{1,2},{4,1},{2,3},{4,1.5}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::ClusterizerCreate(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerSetPoints(s,xy,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunAHC(s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Now we've built our clusterization tree. Rep.z contains information which
      //--- is required to build dendrogram. I-th row of rep.z represents one merge
      //--- operation, with first cluster to merge having index rep.z[I,0] and second
      //--- one having index rep.z[I,1]. Merge result has index NPoints+I.
      //--- Clusters with indexes less than NPoints are single-point initial clusters,
      //--- while ones with indexes from NPoints to 2*NPoints-2 are multi-point
      //--- clusters created during merges.
      //--- In our example, Z=[[2,4], [0,1], [3,6], [5,7]]
      //--- It means that:
      //--- * first, we merge C2=(P2) and C4=(P4),    and create C5=(P2,P4)
      //--- * then, we merge  C2=(P0) and C1=(P1),    and create C6=(P0,P1)
      //--- * then, we merge  C3=(P3) and C6=(P0,P1), and create C7=(P0,P1,P3)
      //--- * finally, we merge C5 and C7 and create C8=(P0,P1,P2,P3,P4)
      //--- Thus, we have following dendrogram:
      //---      ------8-----
      //---      |          |
      //---      |      ----7----
      //---      |      |       |
      //---   ---5---   |    ---6---
      //---   |     |   |    |     |
      //---   P2   P4   P3   P0   P1
      CMatrixInt check;
      check.Resize(4,2);
      check.Set(0,0,2);
      check.Set(0,1,4);
      check.Set(1,0,0);
      check.Set(1,1,1);
      check.Set(2,0,3);
      check.Set(2,1,6);
      check.Set(3,0,5);
      check.Set(3,1,7);
      _TestResult=_TestResult && Doc_Test_Int_Matrix(rep.m_z,check);
      //--- We've built dendrogram above by reordering our dataset.
      //--- Without such reordering it would be impossible to build dendrogram without
      //--- intersections. Luckily, ahcreport structure contains two additional fields
      //--- which help to build dendrogram from your data:
      //--- * rep.p, which contains permutation applied to dataset
      //--- * rep.pm, which contains another representation of merges
      //--- In our example we have:
      //--- * P=[3,4,0,2,1]
      //--- * PZ=[[0,0,1,1,0,0],[3,3,4,4,0,0],[2,2,3,4,0,1],[0,1,2,4,1,2]]
      //--- Permutation array P tells us that P0 should be moved to position 3,
      //--- P1 moved to position 4, P2 moved to position 0 and so on:
      //---   (P0 P1 P2 P3 P4) => (P2 P4 P3 P0 P1)
      //--- Merges array PZ tells us how to perform merges on the sorted dataset.
      //--- One row of PZ corresponds to one merge operations, with first pair of
      //--- elements denoting first of the clusters to merge (start index, end
      //--- index) and next pair of elements denoting second of the clusters to
      //--- merge. Clusters being merged are always adjacent, with first one on
      //--- the left and second one on the right.
      //--- For example, first row of PZ tells us that clusters [0,0] and [1,1] are
      //--- merged (single-point clusters, with first one containing P2 and second
      //--- one containing P4). Third row of PZ tells us that we merge one single-
      //--- point cluster [2,2] with one two-point cluster [3,4].
      //--- There are two more elements in each row of PZ. These are the helper
      //--- elements, which denote HEIGHT (not size) of left and right subdendrograms.
      //--- For example, according to PZ, first two merges are performed on clusterization
      //--- trees of height 0, while next two merges are performed on 0-1 and 1-2
      //--- pairs of trees correspondingly.
        {
         int check1[]={3,4,0,2,1};
         _TestResult=_TestResult && Doc_Test_Int_Vector(rep.m_p,check1);
        }
      check.Resize(4,6);
      check.Set(0,0,0);
      check.Set(0,1,0);
      check.Set(0,2,1);
      check.Set(0,3,1);
      check.Set(0,4,0);
      check.Set(0,5,0);
      check.Set(1,0,3);
      check.Set(1,1,3);
      check.Set(1,2,4);
      check.Set(1,3,4);
      check.Set(1,4,0);
      check.Set(1,5,0);
      check.Set(2,0,2);
      check.Set(2,1,2);
      check.Set(2,2,3);
      check.Set(2,3,4);
      check.Set(2,4,0);
      check.Set(2,5,1);
      check.Set(3,0,0);
      check.Set(3,1,1);
      check.Set(3,2,2);
      check.Set(3,3,4);
      check.Set(3,4,1);
      check.Set(3,5,2);
      _TestResult=_TestResult && Doc_Test_Int_Matrix(rep.m_pm,check);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple k-means clusterization                                    |
//+------------------------------------------------------------------+
void TEST_Clst_KMeans(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- The very simple clusterization example
      //--- We have a set of points in 2D space:
      //---     (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
      //---  |
      //---  |     P3
      //---  |
      //---  | P1
      //---  |             P4
      //---  | P0          P2
      //---  |-------------------------
      //--- We want to perform k-means++ clustering with K=2.
      //--- In order to do that, we:
      //--- * create clusterizer with ClusterizerCreate()
      //--- * set points XY and metric (must be Euclidean, distype=2) with ClusterizerSetPoints()
      //--- * (optional) set number of restarts from random positions to 5
      //--- * run k-means algorithm with clusterizerrunkmeans()
      //--- You may see that clusterization itself is a minor part of the example,
      //--- most of which is dominated by comments :)
      CClusterizerState s;
      CKmeansReport rep;
      matrix<double> XY={{1,1},{1,2},{4,1},{2,3},{4,1.5}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::ClusterizerCreate(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerSetPoints(s,xy,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerSetKMeansLimits(s,5,0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunKMeans(s,2,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- We've performed clusterization, and it succeeded (completion code is +1).
      //--- Now first center is stored in the first row of rep.c, second one is stored
      //--- in the second row. rep.cidx can be used to determine which center is
      //--- closest to some specific point of the dataset.
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,1);
      //--- We called ClusterizerSetPoints() with disttype=2 because k-means++
      //--- algorithm does NOT support metrics other than Euclidean. But what if we
      //--- try to use some other metric?
      //--- We change metric type by calling ClusterizerSetPoints() one more time,
      //--- and try to run k-means algo again. It fails.
      CAlglib::ClusterizerSetPoints(s,xy,0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunKMeans(s,2,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(rep.m_terminationtype,-5);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Clusterization with different linkage types                      |
//+------------------------------------------------------------------+
void TEST_Clst_Linkage(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- We have a set of points in 1D space:
      //---     (P0,P1,P2,P3,P4) = (1, 3, 10, 16, 20)
      //--- We want to perform Agglomerative Hierarchic Clusterization (AHC),
      //--- using either complete or single linkage and Euclidean distance
      //--- (default metric).
      //--- First two steps merge P0/P1 and P3/P4 independently of the linkage type.
      //--- However, third step depends on linkage type being used:
      //--- * in case of complete linkage P2=10 is merged with [P0,P1]
      //--- * in case of single linkage P2=10 is merged with [P3,P4]
      CClusterizerState s;
      CAHCReport rep;
      matrix<double> XY={{1},{3},{10},{16},{20}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CRowInt cidx;
      CRowInt cz;
      CAlglib::ClusterizerCreate(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerSetPoints(s,xy,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- use complete linkage, reduce set down to 2 clusters.
      //--- print clusterization with ClusterizerGetKClusters(2).
      //--- P2 must belong to [P0,P1]
      CAlglib::ClusterizerSetAHCAlgo(s,0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunAHC(s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerGetKClusters(rep,2,cidx,cz);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         int check[]={1,1,1,0,0};
         _TestResult=_TestResult && Doc_Test_Int_Vector(cidx,check);
        }
      //--- use single linkage, reduce set down to 2 clusters.
      //--- print clusterization with ClusterizerGetKClusters(2).
      //--- P2 must belong to [P2,P3]
      CAlglib::ClusterizerSetAHCAlgo(s,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunAHC(s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerGetKClusters(rep,2,cidx,cz);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      int check[]={0,0,1,1,1};
      _TestResult=_TestResult && Doc_Test_Int_Vector(cidx,check);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Clusterization with different metric types                       |
//+------------------------------------------------------------------+
void TEST_Clst_Distance(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- We have three points in 4D space:
      //---     (P0,P1,P2) = ((1, 2, 1, 2), (6, 7, 6, 7), (7, 6, 7, 6))
      //--- We want to try clustering them with different distance functions.
      //--- Distance function is chosen when we add dataset to the clusterizer.
      //--- We can choose several distance types - Euclidean, city block, Chebyshev,
      //--- several correlation measures or user-supplied distance matrix.
      //--- Here we'll try three distances: Euclidean, Pearson correlation,
      //--- user-supplied distance matrix. Different distance functions lead
      //--- to different choices being made by algorithm during clustering.
      CClusterizerState s;
      CAHCReport rep;
      int disttype;
      matrix<double> XY={{1,2,1,2},{6,7,6,7},{7,6,7,6}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::ClusterizerCreate(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- With Euclidean distance function (disttype=2) two closest points
      //--- are P1 and P2, thus:
      //--- * first, we merge P1 and P2 to form C3=[P1,P2]
      //--- * second, we merge P0 and C3 to form C4=[P0,P1,P2]
      disttype=2;
      CAlglib::ClusterizerSetPoints(s,xy,disttype);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunAHC(s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CMatrixInt check;
      check.Resize(2,2);
      check.Set(0,0,1);
      check.Set(0,1,2);
      check.Set(1,0,0);
      check.Set(1,1,3);
      _TestResult=_TestResult && Doc_Test_Int_Matrix(rep.m_z,check);
      //--- With Pearson correlation distance function (disttype=10) situation
      //--- is different - distance between P0 and P1 is zero, thus:
      //--- * first, we merge P0 and P1 to form C3=[P0,P1]
      //--- * second, we merge P2 and C3 to form C4=[P0,P1,P2]
      disttype=10;
      CAlglib::ClusterizerSetPoints(s,xy,disttype);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunAHC(s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      check.Set(0,0,0);
      check.Set(0,1,1);
      check.Set(1,0,2);
      check.Set(1,1,3);
      _TestResult=_TestResult && Doc_Test_Int_Matrix(rep.m_z,check);
      //--- Finally, we try clustering with user-supplied distance matrix:
      //---     [ 0 3 1 ]
      //--- P = [ 3 0 3 ], where P[i,j] = dist(Pi,Pj)
      //---     [ 1 3 0 ]
      //--- * first, we merge P0 and P2 to form C3=[P0,P2]
      //--- * second, we merge P1 and C3 to form C4=[P0,P1,P2]
      CMatrixDouble d;
      matrix<double> D={{0,3,1},{3,0,3},{1,3,0}};
      d=D;
      CAlglib::ClusterizerSetDistances(s,d,true);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunAHC(s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      check.Set(0,0,0);
      check.Set(0,1,2);
      check.Set(1,0,1);
      check.Set(1,1,3);
      _TestResult=_TestResult && Doc_Test_Int_Matrix(rep.m_z,check);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Obtaining K top clusters from clusterization tree                |
//+------------------------------------------------------------------+
void TEST_Clst_KClusters(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- We have a set of points in 2D space:
      //---     (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
      //---  |
      //---  |     P3
      //---  |
      //---  | P1
      //---  |             P4
      //---  | P0          P2
      //---  |-------------------------
      //--- We perform Agglomerative Hierarchic Clusterization (AHC) and we want
      //--- to get top K clusters from clusterization tree for different K.
      CClusterizerState s;
      CAHCReport rep;
      matrix<double> XY={{1,1},{1,2},{4,1},{2,3},{4,1.5}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CRowInt cidx;
      CRowInt cz;
      CAlglib::ClusterizerCreate(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerSetPoints(s,xy,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::ClusterizerRunAHC(s,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- with K=5, every points is assigned to its own cluster:
      //--- C0=P0, C1=P1 and so on...
      CAlglib::ClusterizerGetKClusters(rep,5,cidx,cz);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         int check[]={0,1,2,3,4};
         _TestResult=_TestResult && Doc_Test_Int_Vector(cidx,check);
        }
      //--- with K=1 we have one large cluster C0=[P0,P1,P2,P3,P4,P5]
      CAlglib::ClusterizerGetKClusters(rep,1,cidx,cz);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         int check[]={0,0,0,0,0};
         _TestResult=_TestResult && Doc_Test_Int_Vector(cidx,check);
        }
      //--- with K=3 we have three clusters C0=[P3], C1=[P2,P4], C2=[P0,P1]
      CAlglib::ClusterizerGetKClusters(rep,3,cidx,cz);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      int check[]={2,2,1,0,1};
      _TestResult=_TestResult && Doc_Test_Int_Vector(cidx,check);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple classification with random forests                        |
//+------------------------------------------------------------------+
void TEST_RandomForest_CLS(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- The very simple classification example: classify points (x,y) in 2D space
      //--- as ones with x>=0 and ones with x<0 (y is ignored, but our classifier
      //--- has to find it).
      //--- First, we have to create decision forest builder object, load dataset and
      //--- specify training settings. Our dataset is specified as matrix, which has
      //--- following format:
      //---     x0 y0 class0
      //---     x1 y1 class1
      //---     x2 y2 class2
      //---     ....
      //--- Here xi and yi can be any values (and in fact you can have any number of
      //--- independent variables), and classi MUST be integer number in [0,NClasses)
      //--- range. In our example we denote points with x>=0 as class #0, and
      //--- ones with negative xi as class #1.
      //--- NOTE: if you want to solve regression problem, specify NClasses=1. In
      //---       this case last column of xy can be any numeric value.
      //--- For the sake of simplicity, our example includes only 4-point dataset.
      //--- However, random forests are able to cope with extremely large datasets
      //--- having millions of examples.
      CDecisionForestBuilder builder;
      int nvars=2;
      int nclasses=2;
      int npoints=4;
      matrix<double> XY={{1,1,0},{1,-1,0},{-1,1,1},{-1,-1,1}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::DFBuilderCreate(builder);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::DFBuilderSetDataset(builder,xy,npoints,nvars,nclasses);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- in our example we train decision forest using full sample - it allows us
      //--- to get zero classification error. However, in practical applications smaller
      //--- values are used: 50%, 25%, 5% or even less.
      CAlglib::DFBuilderSetSubsampleRatio(builder,1.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- we train random forest with just one tree; again, in real life situations
      //--- you typically need from 50 to 500 trees.
      int ntrees=1;
      CDecisionForestShell forest;
      CDFReportShell rep;
      CAlglib::DFBuilderBuildRandomForest(builder,ntrees,forest,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- with such settings (100% of the training set is used) you can expect
      //--- zero classification error. Beautiful results, but remember - in real life
      //--- you do not need zero TRAINING SET error, you need good generalization.
      _TestResult=_TestResult && Doc_Test_Real(rep.GetRelClsError(),0.0000,0.00005);
      //--- now, let's perform some simple processing with dfprocess()
      double x[]={+1,0};
      double y[];
      CAlglib::DFProcess(forest,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={+1,0};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.0005);
        }
      //--- another option is to use dfprocess0() which returns just first component
      //--- of the output vector y. ideal for regression problems and binary classifiers.
      double y0=CAlglib::DFProcess0(forest,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(y0,1.000,0.0005);
      //--- finally, you can use dfclassify() which returns most probable class index (i.e. argmax y[i]).
      int i=CAlglib::DFClassify(forest,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(i,0);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple regression with decision forest                           |
//+------------------------------------------------------------------+
void TEST_RandomForest_Reg(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- The very simple regression example: model f(x,y)=x+y
      //--- First, we have to create DF builder object, load dataset and specify
      //--- training settings. Our dataset is specified as matrix, which has following
      //--- format:
      //---     x0 y0 f0
      //---     x1 y1 f1
      //---     x2 y2 f2
      //---     ....
      //--- Here xi and yi can be any values, and fi is a dependent function value.
      //--- NOTE: you can also solve classification problems with DF models, see
      //---       another example for this unit.
      CDecisionForestBuilder builder;
      int nvars=2;
      int nclasses=1;
      int npoints=4;
      matrix<double> XY={{1,1,+2},{1,-1,0},{-1,1,0},{-1,-1,-2}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::DFBuilderCreate(builder);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::DFBuilderSetDataset(builder,xy,npoints,nvars,nclasses);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- in our example we train decision forest using full sample - it allows us
      //--- to get zero classification error. However, in practical applications smaller
      //--- values are used: 50%, 25%, 5% or even less.
      CAlglib::DFBuilderSetSubsampleRatio(builder,1.0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- we train random forest with just one tree; again, in real life situations
      //--- you typically need from 50 to 500 trees.
      int ntrees=1;
      CDecisionForestShell model;
      CDFReportShell rep;
      CAlglib::DFBuilderBuildRandomForest(builder,ntrees,model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- with such settings (full sample is used) you can expect zero RMS error on the
      //--- training set. Beautiful results, but remember - in real life you do not
      //--- need zero TRAINING SET error, you need good generalization.
      _TestResult=_TestResult && Doc_Test_Real(rep.GetRMSError(),0.0000,0.00005);
      //--- now, let's perform some simple processing with dfprocess()
      double x[]={+1,+1};
      double y[];
      CAlglib::DFProcess(model,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={+2};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.0005);
        }
      //--- another option is to use dfprocess0() which returns just first component
      //--- of the output vector y. ideal for regression problems and binary classifiers.
      double y0=CAlglib::DFProcess0(model,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(y0,2.000,0.0005);
      //--- there also exist another convenience function, dfclassify(),
      //--- but it does not work for regression problems - it always returns -1.
      int i=CAlglib::DFClassify(model,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(i,-1);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Linear regression used to build the very basic model and unpack  |
//| coefficients                                                     |
//+------------------------------------------------------------------+
void TEST_LinReg_D_Basic(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
//--- In this example we demonstrate linear fitting by f(x|a) = a*exp(0.5*x).
//--- We have:
//--- * xy - matrix of basic function values (exp(0.5*x)) and expected values
   matrix<double> XY={{0.606531,1.133719},{0.670320,1.306522},{0.740818,1.504604},{0.818731,1.554663},{0.904837,1.884638},{1.000000,2.072436},{1.105171,2.257285},{1.221403,2.534068},{1.349859,2.622017},{1.491825,2.897713},{1.648721,3.219371}};
   CMatrixDouble xy=XY;
   int info;
   int nvars;
   CLinearModelShell model;
   CLRReportShell rep;
   CRowDouble c;
   CAlglib::LRBuildZ(xy,11,1,info,model,rep);
//--- handling exceptions
   if(Func_spoil_scenario(-1,_TestResult))
      _TestResult=_TestResult && Doc_Test_Int(info,1);
   CAlglib::LRUnpack(model,c,nvars);
//--- handling exceptions
   if(Func_spoil_scenario(-1,_TestResult))
     {
      double check[]={1.98650,0.00000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(c,check,0.00005);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| SMA(k) filter                                                    |
//+------------------------------------------------------------------+
void TEST_Filters_D_SMA(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- Here we demonstrate SMA(k) filtering for time series.
      double X[]={5,6,7,8};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      //--- Apply filter.
      //--- We should get [5, 5.5, 6.5, 7.5] as result
      CAlglib::FilterSMA(x,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={5,5.5,6.5,7.5};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| EMA(alpha) filter                                                |
//+------------------------------------------------------------------+
void TEST_Filters_D_EMA(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- Here we demonstrate EMA(0.5) filtering for time series.
      double X[]={5,6,7,8};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      //--- Apply filter.
      //--- We should get [5, 5.5, 6.25, 7.125] as result
      CAlglib::FilterEMA(x,0.5);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={5,5.5,6.25,7.125};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| LRMA(k) filter                                                   |
//+------------------------------------------------------------------+
void TEST_Filters_D_LRMA(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- Here we demonstrate LRMA(3) filtering for time series.
      double X[]={7,8,8,9,12,12};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      //--- Apply filter.
      //--- We should get [7.0000, 8.0000, 8.1667, 8.8333, 11.6667, 12.5000] as result
      CAlglib::FilterLRMA(x,3);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={7.0000,8.0000,8.1667,8.8333,11.6667,12.5000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(x,check,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple SSA analysis demo                                         |
//+------------------------------------------------------------------+
void TEST_SSA_D_Basic(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- Here we demonstrate SSA trend/noise separation for some toy problem:
      //--- small monotonically growing series X are analyzed with 3-tick window
      //--- and "top-K" version of SSA, which selects K largest singular vectors
      //--- for analysis, with K=1.
      CSSAModel s;
      double X[]={0,0.5,1,1,1.5,2};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      //--- First, we create SSA model, set its properties and add dataset.
      //--- We use window with width=3 and configure model to use direct SSA
      //--- algorithm - one which runs exact O(N*W^2) analysis - to extract
      //--- one top singular vector. Well, it is toy problem :)
      //--- NOTE: SSA model may store and analyze more than one sequence
      //---       (say, different sequences may correspond to data collected
      //---       from different devices)
      CAlglib::SSACreate(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSASetWindow(s,3);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAddSequence(s,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSASetAlgoTopKDirect(s,1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Now we begin analysis. Internally SSA model stores everything it needs:
      //--- data, settings, solvers and so on. Right after first call to analysis-
      //--- related function it will analyze dataset, build basis and perform analysis.
      //--- Subsequent calls to analysis functions will reuse previously computed
      //--- basis, unless you invalidate it by changing model settings (or dataset).
      CRowDouble trend;
      CRowDouble noise;
      CAlglib::SSAAnalyzeSequence(s,x,trend,noise);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={0.3815,0.5582,0.7810,1.0794,1.5041,2.0105};
      _TestResult=_TestResult && Doc_Test_Real_Vector(trend,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple SSA forecasting demo                                      |
//+------------------------------------------------------------------+
void TEST_SSA_D_Forecast(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- Here we demonstrate SSA forecasting on some toy problem with clearly
      //--- visible linear trend and small amount of noise.
      CSSAModel s;
      double X[]={0.05,0.96,2.04,3.11,3.97,5.03,5.98,7.02,8.02};
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      //--- First, we create SSA model, set its properties and add dataset.
      //--- We use window with width=3 and configure model to use direct SSA
      //--- algorithm - one which runs exact O(N*W^2) analysis - to extract
      //--- two top singular vectors. Well, it is toy problem :)
      //--- NOTE: SSA model may store and analyze more than one sequence
      //---       (say, different sequences may correspond to data collected
      //---       from different devices)
      CAlglib::SSACreate(s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSASetWindow(s,3);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAddSequence(s,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSASetAlgoTopKDirect(s,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Now we begin analysis. Internally SSA model stores everything it needs:
      //--- data, settings, solvers and so on. Right after first call to analysis-
      //--- related function it will analyze dataset, build basis and perform analysis.
      //--- Subsequent calls to analysis functions will reuse previously computed
      //--- basis, unless you invalidate it by changing model settings (or dataset).
      //--- In this example we show how to use SSAForecastLast() function, which
      //--- predicts changed in the last sequence of the dataset. If you want to
      //--- perform prediction for some other sequence, use ssaforecastsequence().
      CRowDouble trend;
      CAlglib::SSAForecastLast(s,3,trend);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Well, we expected it to be [9,10,11]. There exists some difference,
      //--- which can be explained by the artificial noise in the dataset.
      double check[]={9.0005,9.9322,10.8051};
      _TestResult=_TestResult && Doc_Test_Real_Vector(trend,check,0.005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Real-time SSA algorithm with fast incremental updates            |
//+------------------------------------------------------------------+
void TEST_SSA_D_Realtime(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
     {
      //--- Suppose that you have a constant stream of incoming data, and you want
      //--- to regularly perform singular spectral analysis of this stream.
      //--- One full run of direct algorithm costs O(N*Width^2) operations, so
      //--- the more points you have, the more it costs to rebuild basis from
      //--- scratch.
      //--- Luckily we have incremental SSA algorithm which can perform quick
      //--- updates of already computed basis in O(K*Width^2) ops, where K
      //--- is a number of singular vectors extracted. Usually it is orders of
      //--- magnitude faster than full update of the basis.
      //--- In this example we start from some initial dataset x0. Then we
      //--- start appending elements one by one to the end of the last sequence.
      //--- NOTE: direct algorithm also supports incremental updates, but
      //---       with O(Width^3) cost. Typically K<<Width, so specialized
      //---       incremental algorithm is still faster.
      CSSAModel s1;
      CMatrixDouble a1;
      CRowDouble sv1;
      int w;
      int k;
      double X0[]={0.009,0.976,1.999,2.984,3.977,5.002};
      CRowDouble x0=X0;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x0,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x0,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x0,AL_NEGINF);
      CAlglib::SSACreate(s1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSASetWindow(s1,3);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAddSequence(s1,x0);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- set algorithm to the real-time version of top-K, K=2
      CAlglib::SSASetAlgoTopKRealtime(s1,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- one more interesting feature of the incremental algorithm is "power-up" cycle.
      //--- even with incremental algorithm initial basis calculation costs O(N*Width^2) ops.
      //--- if such startup cost is too high for your real-time app, then you may divide
      //--- initial basis calculation across several model updates. It results in better
      //--- latency at the price of somewhat lesser precision during first few updates.
      CAlglib::SSASetPowerUpLength(s1,3);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- now, after we prepared everything, start to add incoming points one by one;
      //--- in the real life, of course, we will perform some work between subsequent update
      //--- (analyze something, predict, and so on).
      //--- After each append we perform one iteration of the real-time solver. Usually
      //--- one iteration is more than enough to update basis. If you have REALLY tight
      //--- performance constraints, you may specify fractional amount of iterations,
      //--- which means that iteration is performed with required probability.
      double updateits=1.0;
      if(_spoil_scenario==3)
         updateits=AL_NaN;
      if(_spoil_scenario==4)
         updateits=AL_POSINF;
      if(_spoil_scenario==5)
         updateits=AL_NEGINF;
      CAlglib::SSAAppendPointAndUpdate(s1,5.951,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAppendPointAndUpdate(s1,7.074,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAppendPointAndUpdate(s1,7.925,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAppendPointAndUpdate(s1,8.992,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAppendPointAndUpdate(s1,9.942,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAppendPointAndUpdate(s1,11.051,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAppendPointAndUpdate(s1,11.965,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAppendPointAndUpdate(s1,13.047,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAppendPointAndUpdate(s1,13.970,updateits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s1,a1,sv1,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Ok, we have our basis in a1[] and singular values at sv1[].
      //--- But is it good enough? Let's print it.
        {
         matrix<double> check={{0.510607,0.753611},{0.575201,0.058445},{0.639081,-0.654717}};
         _TestResult=_TestResult && Doc_Test_Real_Matrix(a1,check,0.0005);
        }
      //--- Ok, two vectors with 3 components each.
      //--- But how to understand that is it really good basis?
      //--- Let's compare it with direct SSA algorithm on the entire sequence.
      CSSAModel s2;
      CMatrixDouble a2;
      CRowDouble sv2;
      double X2[]={0.009,0.976,1.999,2.984,3.977,5.002,5.951,7.074,7.925,8.992,9.942,11.051,11.965,13.047,13.970};
      CRowDouble x2=X2;
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(x2,AL_NaN);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(x2,AL_POSINF);
      if(_spoil_scenario==8)
         Spoil_Vector_By_Value(x2,AL_NEGINF);
      //---
      CAlglib::SSACreate(s2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSASetWindow(s2,3);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAAddSequence(s2,x2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSASetAlgoTopKDirect(s2,2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::SSAGetBasis(s2,a2,sv2,w,k);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- it is exactly the same as one calculated with incremental approach!
        {
         matrix<double> check={{0.510607,0.753611},{0.575201,0.058445},{0.639081,-0.654717}};
         _TestResult=_TestResult && Doc_Test_Real_Matrix(a2,check,0.0005);
        }
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple classification with KNN model                             |
//+------------------------------------------------------------------+
void TEST_KNN_Cls(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- The very simple classification example: classify points (x,y) in 2D space
      //--- as ones with x>=0 and ones with x<0 (y is ignored, but our classifier
      //--- has to find it).
      //--- First, we have to create KNN builder object, load dataset and specify
      //--- training settings. Our dataset is specified as matrix, which has following
      //--- format:
      //---     x0 y0 class0
      //---     x1 y1 class1
      //---     x2 y2 class2
      //---     ....
      //--- Here xi and yi can be any values (and in fact you can have any number of
      //--- independent variables), and classi MUST be integer number in [0,NClasses)
      //--- range. In our example we denote points with x>=0 as class #0, and
      //--- ones with negative xi as class #1.
      //--- NOTE: if you want to solve regression problem, specify dataset in similar
      //---       format, but with dependent variable(s) instead of class labels. You
      //---       can have dataset with multiple dependent variables, by the way!
      //--- For the sake of simplicity, our example includes only 4-point dataset and
      //--- really simple K=1 nearest neighbor search. Industrial problems typically
      //--- need larger values of K.
      CKNNBuilder builder;
      int nvars=2;
      int nclasses=2;
      int npoints=4;
      matrix<double> XY={{1,1,0},{1,-1,0},{-1,1,1},{-1,-1,1}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::KNNBuilderCreate(builder);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::KNNBuilderSetDatasetCLS(builder,xy,npoints,nvars,nclasses);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- we build KNN model with k=1 and eps=0 (exact k-nn search is performed)
      int k=1;
      double eps=0;
      CKNNModel model;
      CKNNReport rep;
      CAlglib::KNNBuilderBuildKNNModel(builder,k,eps,model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- with such settings (k=1 is used) you can expect zero classification
      //--- error on training set. Beautiful results, but remember - in real life
      //--- you do not need zero TRAINING SET error, you need good generalization.
      _TestResult=_TestResult && Doc_Test_Real(rep.m_RelCLSError,0.0000,0.00005);
      //--- now, let's perform some simple processing with knnprocess()
      double X[]={+1,0};
      CRowDouble x=X;
      CRowDouble y;
      CAlglib::KNNProcess(model,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={+1,0};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.0005);
        }
      //--- another option is to use knnprocess0() which returns just first component
      //--- of the output vector y. ideal for regression problems and binary classifiers.
      double y0=CAlglib::KNNProcess0(model,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(y0,1.000,0.0005);
      //--- finally, you can use knnclassify() which returns most probable class index (i.e. argmax y[i]).
      int i=CAlglib::KNNClassify(model,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(i,0);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Simple classification with KNN model                             |
//+------------------------------------------------------------------+
void TEST_KNN_Reg(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- The very simple regression example: model f(x,y)=x+y
      //--- First, we have to create KNN builder object, load dataset and specify
      //--- training settings. Our dataset is specified as matrix, which has following
      //--- format:
      //---     x0 y0 f0
      //---     x1 y1 f1
      //---     x2 y2 f2
      //---     ....
      //--- Here xi and yi can be any values, and fi is a dependent function value.
      //--- By the way, with KNN algorithm you can even model functions with multiple
      //--- dependent variables!
      //--- NOTE: you can also solve classification problems with KNN models, see
      //---       another example for this unit.
      //--- For the sake of simplicity, our example includes only 4-point dataset and
      //--- really simple K=1 nearest neighbor search. Industrial problems typically
      //--- need larger values of K.
      CKNNBuilder builder;
      int nvars=2;
      int nout=1;
      int npoints=4;
      matrix<double> XY={{1,1,+2},{1,-1,0},{-1,1,0},{-1,-1,-2}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::KNNBuilderCreate(builder);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::KNNBuilderSetDatasetReg(builder,xy,npoints,nvars,nout);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- we build KNN model with k=1 and eps=0 (exact k-nn search is performed)
      int k=1;
      double eps=0;
      CKNNModel model;
      CKNNReport rep;
      CAlglib::KNNBuilderBuildKNNModel(builder,k,eps,model,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- with such settings (k=1 is used) you can expect zero RMS error on the
      //--- training set. Beautiful results, but remember - in real life you do not
      //--- need zero TRAINING SET error, you need good generalization.
      _TestResult=_TestResult && Doc_Test_Real(rep.m_RMSError,0.0000,0.00005);
      //--- now, let's perform some simple processing with knnprocess()
      double X[]={+1,+1};
      CRowDouble x=X;
      CRowDouble y;
      CAlglib::KNNProcess(model,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={+2};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.0005);
        }
      //--- another option is to use knnprocess0() which returns just first component
      //--- of the output vector y. ideal for regression problems and binary classifiers.
      double y0=CAlglib::KNNProcess0(model,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(y0,2.000,0.0005);
      //--- there also exist another convenience function, knnclassify(),
      //--- but it does not work for regression problems - it always returns -1.
      int i=CAlglib::KNNClassify(model,x);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Int(i,-1);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Regression problem with one output (2=>1)                        |
//+------------------------------------------------------------------+
void TEST_NN_Regr(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- The very simple example on neural network: network is trained to reproduce
      //--- small 2x2 multiplication table.
      //--- NOTE: we use network with excessive amount of neurons, which guarantees
      //---       almost exact reproduction of the training set. Generalization ability
      //---       of such network is rather low, but we are not concerned with such
      //---       questions in this basic demo.
      CMLPTrainer trn;
      CMultilayerPerceptronShell network;
      CMLPReportShell rep;
      //--- Training set:
      //--- * one row corresponds to one record A*B=C in the multiplication table
      //--- * first two columns store A and B, last column stores C
      //--- [1 * 1 = 1]
      //--- [1 * 2 = 2]
      //--- [2 * 1 = 2]
      //--- [2 * 2 = 4]
      matrix<double> XY={{1,1,1},{1,2,2},{2,1,2},{2,2,4}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      //--- Network is created.
      //--- Trainer object is created.
      //--- Dataset is attached to trainer object.
      CAlglib::MLPCreateTrainer(2,1,trn);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPCreate1(2,5,1,network);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPSetDataset(trn,xy,4);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Network is trained with 5 restarts from random positions
      CAlglib::MLPTrainNetwork(trn,network,5,rep);
      //--- 2*2=?
      double x[]={2,2};
      double y[]={0};
      CAlglib::MLPProcess(network,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={4.000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Regression problem with multiple outputs (2=>2)                  |
//+------------------------------------------------------------------+
void TEST_NN_Regr_N(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- Network with 2 inputs and 2 outputs is trained to reproduce vector function:
      //---     (x0,x1) => (x0+x1, x0*x1)
      //--- Informally speaking, we want neural network to simultaneously calculate
      //--- both sum of two numbers and their product.
      //--- NOTE: we use network with excessive amount of neurons, which guarantees
      //---       almost exact reproduction of the training set. Generalization ability
      //---       of such network is rather low, but we are not concerned with such
      //---       questions in this basic demo.
      CMLPTrainer trn;
      CMultilayerPerceptronShell network;
      CMLPReportShell rep;
      //--- Training set. One row corresponds to one record [A,B,A+B,A*B].
      //--- [ 1   1  1+1  1*1 ]
      //--- [ 1   2  1+2  1*2 ]
      //--- [ 2   1  2+1  2*1 ]
      //--- [ 2   2  2+2  2*2 ]
      matrix<double> XY={{1,1,2,1},{1,2,3,2},{2,1,3,2},{2,2,4,4}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      //--- Network is created.
      //--- Trainer object is created.
      //--- Dataset is attached to trainer object.
      CAlglib::MLPCreateTrainer(2,2,trn);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPCreate1(2,5,2,network);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPSetDataset(trn,xy,4);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Network is trained with 5 restarts from random positions
      CAlglib::MLPTrainNetwork(trn,network,5,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- 2+1=?
      //--- 2*1=?
      double x[]={2,1};
      double y[]={0,0};
      CAlglib::MLPProcess(network,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={3.000,2.000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Binary classification problem                                    |
//+------------------------------------------------------------------+
void TEST_NN_Cls2(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- Suppose that we want to classify numbers as positive (class 0) and negative
      //--- (class 1). We have training set which includes several strictly positive
      //--- or negative numbers - and zero.
      //--- The problem is that we are not sure how to classify zero, so from time to
      //--- time we mark it as positive or negative (with equal probability). Other
      //--- numbers are marked in pure deterministic setting. How will neural network
      //--- cope with such classification task?
      //--- NOTE: we use network with excessive amount of neurons, which guarantees
      //---       almost exact reproduction of the training set. Generalization ability
      //---       of such network is rather low, but we are not concerned with such
      //---       questions in this basic demo.
      CMLPTrainer trn;
      CMultilayerPerceptronShell network;
      CMLPReportShell rep;
      double x[]={0};
      double y[]={0,0};
      //--- Training set. One row corresponds to one record [A => class(A)].
      //--- Classes are denoted by numbers from 0 to 1, where 0 corresponds to positive
      //--- numbers and 1 to negative numbers.
      //--- [ +1  0]
      //--- [ +2  0]
      //--- [ -1  1]
      //--- [ -2  1]
      //--- [  0  0]   !! sometimes we classify 0 as positive, sometimes as negative
      //--- [  0  1]   !!
      matrix<double> XY={{+1,0},{+2,0},{-1,1},{-2,1},{0,0},{0,1}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      //--- When we solve classification problems, everything is slightly different from
      //--- the regression ones:
      //--- 1. Network is created. Because we solve classification problem, we use
      //---    mlpcreatec1() function instead of mlpcreate1(). This function creates
      //---    classifier network with SOFTMAX-normalized outputs. This network returns
      //---    vector of class membership probabilities which are normalized to be
      //---    non-negative and sum to 1.0
      //--- 2. We use mlpcreatetrainercls() function instead of mlpcreatetrainer() to
      //---    create trainer object. Trainer object process dataset and neural network
      //---    slightly differently to account for specifics of the classification
      //---    problems.
      //--- 3. Dataset is attached to trainer object. Note that dataset format is slightly
      //---    different from one used for regression.
      CAlglib::MLPCreateTrainerCls(1,2,trn);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPCreateC1(1,5,2,network);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPSetDataset(trn,xy,6);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Network is trained with 5 restarts from random positions
      CAlglib::MLPTrainNetwork(trn,network,5,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Test our neural network on strictly positive and strictly negative numbers.
      //--- IMPORTANT! Classifier network returns class membership probabilities instead
      //--- of class indexes. Network returns two values (probabilities) instead of one
      //--- (class index).
      //--- Thus, for +1 we expect to get [P0,P1] = [1,0], where P0 is probability that
      //--- number is positive (belongs to class 0), and P1 is probability that number
      //--- is negative (belongs to class 1).
      //--- For -1 we expect to get [P0,P1] = [0,1]
      //--- Following properties are guaranteed by network architecture:
      //--- * P0>=0, P1>=0   non-negativity
      //--- * P0+P1=1        normalization
      x[0]=1;
      CAlglib::MLPProcess(network,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={1.000,0.000};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.05);
        }
      x[0]=-1;
      CAlglib::MLPProcess(network,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={0.000,1.000};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.05);
        }
      //--- But what our network will return for 0, which is between classes 0 and 1?
      //--- In our dataset it has two different marks assigned (class 0 AND class 1).
      //--- So network will return something average between class 0 and class 1:
      //---     0 => [0.5, 0.5]
      x[0]=0;
      CAlglib::MLPProcess(network,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={0.500,0.500};
      _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Multiclass classification problem                                |
//+------------------------------------------------------------------+
void TEST_NN_Cls3(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- Suppose that we want to classify numbers as positive (class 0) and negative
      //--- (class 1). We also have one more class for zero (class 2).
      //--- NOTE: we use network with excessive amount of neurons, which guarantees
      //---       almost exact reproduction of the training set. Generalization ability
      //---       of such network is rather low, but we are not concerned with such
      //---       questions in this basic demo.
      CMLPTrainer trn;
      CMultilayerPerceptronShell network;
      CMLPReportShell rep;
      double x[]={0};
      double y[]={0,0,0};
      //--- Training set. One row corresponds to one record [A => class(A)].
      //--- Classes are denoted by numbers from 0 to 2, where 0 corresponds to positive
      //--- numbers, 1 to negative numbers, 2 to zero
      //--- [ +1  0]
      //--- [ +2  0]
      //--- [ -1  1]
      //--- [ -2  1]
      //--- [  0  2]
      matrix<double> XY={{+1,0},{+2,0},{-1,1},{-2,1},{0,2}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      //--- When we solve classification problems, everything is slightly different from
      //--- the regression ones:
      //--- 1. Network is created. Because we solve classification problem, we use
      //---    mlpcreatec1() function instead of mlpcreate1(). This function creates
      //---    classifier network with SOFTMAX-normalized outputs. This network returns
      //---    vector of class membership probabilities which are normalized to be
      //---    non-negative and sum to 1.0
      //--- 2. We use mlpcreatetrainercls() function instead of mlpcreatetrainer() to
      //---    create trainer object. Trainer object process dataset and neural network
      //---    slightly differently to account for specifics of the classification
      //---    problems.
      //--- 3. Dataset is attached to trainer object. Note that dataset format is slightly
      //---    different from one used for regression.
      CAlglib::MLPCreateTrainerCls(1,3,trn);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPCreateC1(1,5,3,network);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPSetDataset(trn,xy,5);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Network is trained with 5 restarts from random positions
      CAlglib::MLPTrainNetwork(trn,network,5,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Test our neural network on strictly positive and strictly negative numbers.
      //--- IMPORTANT! Classifier network returns class membership probabilities instead
      //--- of class indexes. Network returns three values (probabilities) instead of one
      //--- (class index).
      //--- Thus, for +1 we expect to get [P0,P1,P2] = [1,0,0],
      //--- for -1 we expect to get [P0,P1,P2] = [0,1,0],
      //--- and for 0 we will get [P0,P1,P2] = [0,0,1].
      //--- Following properties are guaranteed by network architecture:
      //--- * P0>=0, P1>=0, P2>=0    non-negativity
      //--- * P0+P1+P2=1             normalization
      x[0]=1;
      CAlglib::MLPProcess(network,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={1.000,0.000,0.000};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.05);
        }
      x[0]=-1;
      CAlglib::MLPProcess(network,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
        {
         double check[]={0.000,1.000,0.000};
         _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.05);
        }
      x[0]=0;
      CAlglib::MLPProcess(network,x,y);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      double check[]={0.000,0.000,1.000};
      _TestResult=_TestResult && Doc_Test_Real_Vector(y,check,0.05);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Advanced example on trainer object                               |
//+------------------------------------------------------------------+
void TEST_NN_TrainerObject(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
     {
      //--- Trainer object is used to train network. It stores dataset, training settings,
      //--- and other information which is NOT part of neural network. You should use
      //--- trainer object as follows:
      //--- (1) you create trainer object and specify task type (classification/regression)
      //---     and number of inputs/outputs
      //--- (2) you add dataset to the trainer object
      //--- (3) you may change training settings (stopping criteria or weight decay)
      //--- (4) finally, you may train one or more networks
      //--- You may interleave stages 2...4 and repeat them many times. Trainer object
      //--- remembers its internal state and can be used several times after its creation
      //--- and initialization.
      CMLPTrainer trn;
      //--- Stage 1: object creation.
      //--- We have to specify number of inputs and outputs. Trainer object can be used
      //--- only for problems with same number of inputs/outputs as was specified during
      //--- its creation.
      //--- In case you want to train SOFTMAX-normalized network which solves classification
      //--- problems,  you  must  use  another  function  to  create  trainer  object:
      //--- mlpcreatetrainercls().
      //--- Below we create trainer object which can be used to train regression networks
      //--- with 2 inputs and 1 output.
      CAlglib::MLPCreateTrainer(2,1,trn);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Stage 2: specification of the training set
      //--- By default trainer object stores empty dataset. So to solve your non-empty problem
      //--- you have to set dataset by passing to trainer dense or sparse matrix.
      //--- One row of the matrix corresponds to one record A*B=C in the multiplication table.
      //--- First two columns store A and B, last column stores C
      //---     [1 * 1 = 1]   [ 1 1 1 ]
      //---     [1 * 2 = 2]   [ 1 2 2 ]
      //---     [2 * 1 = 2] = [ 2 1 2 ]
      //---     [2 * 2 = 4]   [ 2 2 4 ]
      matrix<double> XY={{1,1,1},{1,2,2},{2,1,2},{2,2,4}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      CAlglib::MLPSetDataset(trn,xy,4);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Stage 3: modification of the training parameters.
      //--- You may modify parameters like weights decay or stopping criteria:
      //--- * we set moderate weight decay
      //--- * we choose iterations limit as stopping condition (another condition - step size -
      //---   is zero, which means than this condition is not active)
      double wstep=0.000;
      if(_spoil_scenario==3)
         wstep=AL_NaN;
      if(_spoil_scenario==4)
         wstep=AL_POSINF;
      if(_spoil_scenario==5)
         wstep=AL_NEGINF;
      int maxits=100;
      CAlglib::MLPSetDecay(trn,0.01);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPSetCond(trn,wstep,maxits);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Stage 4: training.
      //--- We will train several networks with different architecture using same trainer object.
      //--- We may change training parameters or even dataset, so different networks are trained
      //--- differently. But in this simple example we will train all networks with same settings.
      //--- We create and train three networks:
      //--- * network 1 has 2x1 architecture     (2 inputs, no hidden neurons, 1 output)
      //--- * network 2 has 2x5x1 architecture   (2 inputs, 5 hidden neurons, 1 output)
      //--- * network 3 has 2x5x5x1 architecture (2 inputs, two hidden layers, 1 output)
      //--- NOTE: these networks solve regression problems. For classification problems you
      //---       should use mlpcreatec0/c1/c2 to create neural networks which have SOFTMAX-
      //---       normalized outputs.
      CMultilayerPerceptronShell net1;
      CMultilayerPerceptronShell net2;
      CMultilayerPerceptronShell net3;
      CMLPReportShell rep;
      CAlglib::MLPCreate0(2,1,net1);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPCreate1(2,5,1,net2);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPCreate2(2,5,5,1,net3);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPTrainNetwork(trn,net1,5,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPTrainNetwork(trn,net2,5,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPTrainNetwork(trn,net3,5,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Early stopping ensembles                                         |
//+------------------------------------------------------------------+
void TEST_NN_Ensembles_ES(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
     {
      //--- This example shows how to train early stopping ensebles.
      CMLPTrainer trn;
      CMLPEnsembleShell ensemble;
      CMLPReportShell rep;
      //--- Training set: f(x)=1/(x^2+1)
      //--- One row corresponds to one record [x,f(x)]
      matrix<double> XY={{-2.0,0.2},{-1.6,0.3},{-1.3,0.4},{-1,0.5},{-0.6,0.7},{-0.3,0.9},{0,1},{2.0,0.2},{1.6,0.3},{1.3,0.4},{1,0.5},{0.6,0.7},{0.3,0.9}};
      CMatrixDouble xy=XY;
      if(_spoil_scenario==0)
         Spoil_Matrix_By_Value(xy,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Matrix_By_Value(xy,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Matrix_By_Value(xy,AL_NEGINF);
      //--- Trainer object is created.
      //--- Dataset is attached to trainer object.
      //--- NOTE: it is not good idea to use early stopping ensemble on sample
      //---       as small as ours (13 examples). It is done for demonstration
      //---       purposes only. Ensemble training algorithm won't find good
      //---       solution on such small sample.
      CAlglib::MLPCreateTrainer(1,1,trn);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPSetDataset(trn,xy,13);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- Ensemble is created and trained. Each of 50 network is trained
      //--- with 5 restarts.
      CAlglib::MLPECreate1(1,4,1,50,ensemble);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      CAlglib::MLPTrainEnsembleES(trn,ensemble,5,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| Monotone interpolation                                           |
//+------------------------------------------------------------------+
void TEST_Spline1D_D_Monotone(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
     {
      //--- Spline built witn spline1dbuildcubic() can be non-monotone even when
      //--- Y-values form monotone sequence. Say, for x=[0,1,2] and y=[0,1,1]
      //--- cubic spline will monotonically grow until x=1.5 and then start
      //--- decreasing.
      //--- That's why ALGLIB provides special spline construction function
      //--- which builds spline which preserves monotonicity of the original
      //--- dataset.
      //--- NOTE: in case original dataset is non-monotonic, ALGLIB splits it
      //--- into monotone subsequences and builds piecewise monotonic spline.
      //
      vector<double> X={ 0,1,2 };
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Adding_Element(x);
      if(_spoil_scenario==4)
         Spoil_Vector_By_Deleting_Element(x);
      vector<double> Y={ 0,1,1 };
      CRowDouble y=Y;
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==8)
         Spoil_Vector_By_Adding_Element(y);
      if(_spoil_scenario==9)
         Spoil_Vector_By_Deleting_Element(y);
      CSpline1DInterpolantShell s;
      //--- build spline
      CAlglib::Spline1DBuildMonotone(x,y,s);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      //--- calculate S at x = [-0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
      //--- you may see that spline is really monotonic
      double v=CAlglib::Spline1DCalc(s,-0.5);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,0.0000,0.00005);
      v=CAlglib::Spline1DCalc(s,0.0);
      _TestResult=_TestResult && Doc_Test_Real(v,0.0000,0.00005);
      v=CAlglib::Spline1DCalc(s,+0.5);
      _TestResult=_TestResult && Doc_Test_Real(v,0.5000,0.00005);
      v=CAlglib::Spline1DCalc(s,1.0);
      _TestResult=_TestResult && Doc_Test_Real(v,1.0000,0.00005);
      v=CAlglib::Spline1DCalc(s,1.5);
      _TestResult=_TestResult && Doc_Test_Real(v,1.0000,0.00005);
      v=CAlglib::Spline1DCalc(s,2.0);
      _TestResult=_TestResult && Doc_Test_Real(v,1.0000,0.00005);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| 4-parameter logistic fitting                                     |
//+------------------------------------------------------------------+
void TEST_LSFit_T_4pl(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<8; _spoil_scenario++)
     {
      vector<double> X={ 1,2,3,4,5,6,7,8 };
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      vector<double> Y={ 0.06313223,0.44552624,0.61838364,0.71385108,0.77345838,0.81383140,0.84280033,0.86449822 };
      CRowDouble y=Y;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      int n=8;
      double a;
      double b;
      double c;
      double d;
      CLSFitReportShell rep;
      //--- Test logisticfit4() on carefully designed data with a priori known answer.
      CAlglib::LogisticFit4(x,y,n,a,b,c,d,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(a,-1.000,0.01);
      _TestResult=_TestResult && Doc_Test_Real(b,1.200,0.01);
      _TestResult=_TestResult && Doc_Test_Real(c,0.900,0.01);
      _TestResult=_TestResult && Doc_Test_Real(d,1.000,0.01);
      //--- Evaluate model at point x=0.5
      double v;
      v=CAlglib::LogisticCalc4(0.5,a,b,c,d);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(v,-0.33874308,0.001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
//| 5-parameter logistic fitting                                     |
//+------------------------------------------------------------------+
void TEST_LSFit_T_5pl(int &_spoil_scenario,bool &_TestResult,bool &_TotalResult)
  {
   _TestResult=true;
   for(_spoil_scenario=-1; _spoil_scenario<8; _spoil_scenario++)
     {
      vector<double> X={ 1,2,3,4,5,6,7,8 };
      CRowDouble x=X;
      if(_spoil_scenario==0)
         Spoil_Vector_By_Value(x,AL_NaN);
      if(_spoil_scenario==1)
         Spoil_Vector_By_Value(x,AL_POSINF);
      if(_spoil_scenario==2)
         Spoil_Vector_By_Value(x,AL_NEGINF);
      if(_spoil_scenario==3)
         Spoil_Vector_By_Deleting_Element(x);
      vector<double> Y={ 0.1949776139,0.5710060208,0.726002637,0.8060434158,0.8534547965,0.8842071579,0.9054773317,0.9209088299 };
      CRowDouble y=Y;
      if(_spoil_scenario==4)
         Spoil_Vector_By_Value(y,AL_NaN);
      if(_spoil_scenario==5)
         Spoil_Vector_By_Value(y,AL_POSINF);
      if(_spoil_scenario==6)
         Spoil_Vector_By_Value(y,AL_NEGINF);
      if(_spoil_scenario==7)
         Spoil_Vector_By_Deleting_Element(y);
      int n=8;
      double a;
      double b;
      double c;
      double d;
      double g;
      CLSFitReportShell rep;
      //--- Test logisticfit5() on carefully designed data with a priori known answer.
      CAlglib::LogisticFit5(x,y,n,a,b,c,d,g,rep);
      //--- handling exceptions
      if(!Func_spoil_scenario(_spoil_scenario,_TestResult))
         continue;
      _TestResult=_TestResult && Doc_Test_Real(a,-1.000,0.01);
      _TestResult=_TestResult && Doc_Test_Real(b,1.200,0.01);
      _TestResult=_TestResult && Doc_Test_Real(c,0.900,0.01);
      _TestResult=_TestResult && Doc_Test_Real(d,1.000,0.01);
      _TestResult=_TestResult && Doc_Test_Real(g,1.200,0.01);
      //--- Evaluate model at point x=0.5
      double v;
      v=CAlglib::LogisticCalc5(0.5,a,b,c,d,g);
      _TestResult=_TestResult && Doc_Test_Real(v,-0.2354656824,0.001);
      _TestResult=_TestResult && (_spoil_scenario==-1);
     }
   if(!_TestResult)
      PrintFormat("%-32s FAILED",__FUNCTION__);
   _TotalResult=_TotalResult && _TestResult;
  }
//+------------------------------------------------------------------+
