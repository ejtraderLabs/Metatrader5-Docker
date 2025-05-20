//+------------------------------------------------------------------+
//|                                               TestInterfaces.mq5 |
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
//| This program is free software; you can redistribute it and/or    |
//| modify it under the terms of the GNU General Public License as   |
//| published by the Free Software Foundation (www.fsf.org); either  |
//| version 2 of the License, or (at your option) any later version. |
//|                                                                  |
//| This program is distributed in the hope that it will be useful,  |
//| but WITHOUT ANY WARRANTY; without even the implied warranty of   |
//| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     |
//| GNU General Public License for more details.                     |
//+------------------------------------------------------------------+
#include "TestInterfaces.mqh"

#define CHECK_RESULT if(!_TestResult) PrintFormat("Test %d failed",test); test++;
//+------------------------------------------------------------------+
//| Testing script                                                   |
//+------------------------------------------------------------------+
void OnStart()
  {
//--- total result
   bool _TotalResult=true;
//--- test result
   bool _TestResult;
   int  test=1;
//--- spoil scenario
   int  _spoil_scenario;
   Print("MQL5 interface tests. Please wait...");
   Print("0/152");
//--- testing
   TEST_Ablas_D_Gemm(_spoil_scenario,_TestResult,_TotalResult);          //1
   CHECK_RESULT
   TEST_Ablas_D_Syrk(_spoil_scenario,_TestResult,_TotalResult);          //2
   CHECK_RESULT
   TEST_Ablas_T_Complex(_spoil_scenario,_TestResult,_TotalResult);       //3
   CHECK_RESULT
   TEST_Sparse_D_1(_spoil_scenario,_TestResult,_TotalResult);            //4
   CHECK_RESULT
   TEST_Sparse_D_CRS(_spoil_scenario,_TestResult,_TotalResult);          //5
   CHECK_RESULT
   TEST_SolveSKS_D_1(_spoil_scenario,_TestResult,_TotalResult);          //6
   CHECK_RESULT
   TEST_LinCG_D_1(_spoil_scenario,_TestResult,_TotalResult);             //7
   CHECK_RESULT
   TEST_LinLSQR_D_1(_spoil_scenario,_TestResult,_TotalResult);           //8
   CHECK_RESULT
   TEST_NNeighbor_D_1(_spoil_scenario,_TestResult,_TotalResult);         //9
   CHECK_RESULT
   TEST_NNeighbor_T_2(_spoil_scenario,_TestResult,_TotalResult);         //10
   CHECK_RESULT
   TEST_NNeighbor_D_2(_spoil_scenario,_TestResult,_TotalResult);         //11
   CHECK_RESULT
   TEST_BaseStat_D_Base(_spoil_scenario,_TestResult,_TotalResult);       //12
   CHECK_RESULT
   TEST_BaseStat_D_C2(_spoil_scenario,_TestResult,_TotalResult);         //13
   CHECK_RESULT
   TEST_BaseStat_D_CM(_spoil_scenario,_TestResult,_TotalResult);         //14
   CHECK_RESULT
   TEST_BaseStat_D_CM2(_spoil_scenario,_TestResult,_TotalResult);        //15
   CHECK_RESULT
   TEST_BaseStat_T_Base(_spoil_scenario,_TestResult,_TotalResult);       //16
   CHECK_RESULT
   TEST_BaseStat_T_CovCorr(_spoil_scenario,_TestResult,_TotalResult);    //17
   CHECK_RESULT
   TEST_MatInv_D_R1(_spoil_scenario,_TestResult,_TotalResult);           //18
   CHECK_RESULT
   TEST_MatInv_D_C1(_spoil_scenario,_TestResult,_TotalResult);           //19
   CHECK_RESULT
   TEST_MatInv_D_SPD1(_spoil_scenario,_TestResult,_TotalResult);         //20
   CHECK_RESULT
   TEST_MatInv_D_HPD1(_spoil_scenario,_TestResult,_TotalResult);         //21
   CHECK_RESULT
   TEST_MatInv_T_R1(_spoil_scenario,_TestResult,_TotalResult);           //22
   CHECK_RESULT
   TEST_MatInv_T_C1(_spoil_scenario,_TestResult,_TotalResult);           //23
   CHECK_RESULT
   TEST_MatInv_E_SPD1(_spoil_scenario,_TestResult,_TotalResult);         //24
   CHECK_RESULT
   TEST_MatInv_E_HPD1(_spoil_scenario,_TestResult,_TotalResult);         //25
   CHECK_RESULT
   TEST_MinCG_D_1(_spoil_scenario,_TestResult,_TotalResult);             //26
   CHECK_RESULT
   TEST_MinCG_D_2(_spoil_scenario,_TestResult,_TotalResult);             //27
   CHECK_RESULT
   TEST_MinCG_NumDiff(_spoil_scenario,_TestResult,_TotalResult);         //28
   CHECK_RESULT
   TEST_MinCG_FTRIM(_spoil_scenario,_TestResult,_TotalResult);           //29
   CHECK_RESULT
   TEST_MinBLEIC_D_1(_spoil_scenario,_TestResult,_TotalResult);          //30
   CHECK_RESULT
   TEST_MinBLEIC_D_2(_spoil_scenario,_TestResult,_TotalResult);          //31
   CHECK_RESULT
   TEST_MinBLEIC_NumDiff(_spoil_scenario,_TestResult,_TotalResult);      //32
   CHECK_RESULT
   TEST_MinBLEIC_FTRIM(_spoil_scenario,_TestResult,_TotalResult);        //33
   CHECK_RESULT
   TEST_MCPD_Simple1(_spoil_scenario,_TestResult,_TotalResult);          //34
   CHECK_RESULT
   TEST_MCPD_Simple2(_spoil_scenario,_TestResult,_TotalResult);          //35
   CHECK_RESULT
   TEST_MinLBFGS_D_1(_spoil_scenario,_TestResult,_TotalResult);          //36
   CHECK_RESULT
   TEST_MinLBFGS_D_2(_spoil_scenario,_TestResult,_TotalResult);          //37
   CHECK_RESULT
   TEST_MinLBFGS_NumDiff(_spoil_scenario,_TestResult,_TotalResult);      //38
   CHECK_RESULT
   TEST_MinLBFGS_FTRIM(_spoil_scenario,_TestResult,_TotalResult);        //39
   CHECK_RESULT
   TEST_ODESolver_D1(_spoil_scenario,_TestResult,_TotalResult);          //40
   CHECK_RESULT
   TEST_FFT_Complex_D1(_spoil_scenario,_TestResult,_TotalResult);        //41
   CHECK_RESULT
   TEST_FFT_Complex_D2(_spoil_scenario,_TestResult,_TotalResult);        //42
   CHECK_RESULT
   TEST_FFT_Real_D1(_spoil_scenario,_TestResult,_TotalResult);           //43
   CHECK_RESULT
   TEST_FFT_Real_D2(_spoil_scenario,_TestResult,_TotalResult);           //44
   CHECK_RESULT
   TEST_FFT_Complex_E1(_spoil_scenario,_TestResult,_TotalResult);        //45
   CHECK_RESULT
   TEST_AutoGK_D1(_spoil_scenario,_TestResult,_TotalResult);             //46
   CHECK_RESULT
   TEST_PolInt_D_CalcDiff(_spoil_scenario,_TestResult,_TotalResult);     //47
   CHECK_RESULT
   TEST_PolInt_D_Conv(_spoil_scenario,_TestResult,_TotalResult);         //48
   CHECK_RESULT
   TEST_PolInt_D_Spec(_spoil_scenario,_TestResult,_TotalResult);         //49
   CHECK_RESULT
   TEST_PolInt_T_1(_spoil_scenario,_TestResult,_TotalResult);            //50
   CHECK_RESULT
//--- 50 blocks were successful
   Print("50/152");
   TEST_PolInt_T_2(_spoil_scenario,_TestResult,_TotalResult);            //51
   CHECK_RESULT
   TEST_PolInt_T_3(_spoil_scenario,_TestResult,_TotalResult);            //52
   CHECK_RESULT
   TEST_PolInt_T_4(_spoil_scenario,_TestResult,_TotalResult);            //53
   CHECK_RESULT
   TEST_PolInt_T_5(_spoil_scenario,_TestResult,_TotalResult);            //54
   CHECK_RESULT
   TEST_PolInt_T_6(_spoil_scenario,_TestResult,_TotalResult);            //55
   CHECK_RESULT
   TEST_PolInt_T_7(_spoil_scenario,_TestResult,_TotalResult);            //56
   CHECK_RESULT
   TEST_PolInt_T_8(_spoil_scenario,_TestResult,_TotalResult);            //57
   CHECK_RESULT
   TEST_PolInt_T_9(_spoil_scenario,_TestResult,_TotalResult);            //58
   CHECK_RESULT
   TEST_PolInt_T_10(_spoil_scenario,_TestResult,_TotalResult);           //59
   CHECK_RESULT
   TEST_PolInt_T_11(_spoil_scenario,_TestResult,_TotalResult);           //60
   CHECK_RESULT
   TEST_PolInt_T_12(_spoil_scenario,_TestResult,_TotalResult);           //61
   CHECK_RESULT
   TEST_PolInt_T_13(_spoil_scenario,_TestResult,_TotalResult);           //62
   CHECK_RESULT
   TEST_Spline1D_D_Linear(_spoil_scenario,_TestResult,_TotalResult);     //63
   CHECK_RESULT
   TEST_Spline1D_D_Cubic(_spoil_scenario,_TestResult,_TotalResult);      //64
   CHECK_RESULT
   TEST_Spline1D_D_GridDiff(_spoil_scenario,_TestResult,_TotalResult);   //65
   CHECK_RESULT
   TEST_Spline1D_D_ConvDiff(_spoil_scenario,_TestResult,_TotalResult);   //66
   CHECK_RESULT
   TEST_MinQP_D_U1(_spoil_scenario,_TestResult,_TotalResult);            //67
   CHECK_RESULT
   TEST_MinQP_D_BC1(_spoil_scenario,_TestResult,_TotalResult);           //68
   CHECK_RESULT
   TEST_MinLM_D_V(_spoil_scenario,_TestResult,_TotalResult);             //69
   CHECK_RESULT
   TEST_MinLM_D_VJ(_spoil_scenario,_TestResult,_TotalResult);            //70
   CHECK_RESULT
   TEST_MinLM_D_FGH(_spoil_scenario,_TestResult,_TotalResult);           //71
   CHECK_RESULT
   TEST_MinLM_D_VB(_spoil_scenario,_TestResult,_TotalResult);            //72
   CHECK_RESULT
   TEST_MinLM_D_Restarts(_spoil_scenario,_TestResult,_TotalResult);      //73
   CHECK_RESULT
   TEST_MinLM_T_1(_spoil_scenario,_TestResult,_TotalResult);             //74
   CHECK_RESULT
   TEST_MinLM_T_2(_spoil_scenario,_TestResult,_TotalResult);             //75
   CHECK_RESULT
   TEST_LSFit_D_NLF(_spoil_scenario,_TestResult,_TotalResult);           //76
   CHECK_RESULT
   TEST_LSFit_D_NLFG(_spoil_scenario,_TestResult,_TotalResult);          //77
   CHECK_RESULT
   TEST_LSFit_D_NLFGH(_spoil_scenario,_TestResult,_TotalResult);         //78
   CHECK_RESULT
   TEST_LSFit_D_NLFB(_spoil_scenario,_TestResult,_TotalResult);          //79
   CHECK_RESULT
   TEST_LSFit_D_NLScale(_spoil_scenario,_TestResult,_TotalResult);       //80
   CHECK_RESULT
   TEST_LSFit_D_Lin(_spoil_scenario,_TestResult,_TotalResult);           //81
   CHECK_RESULT
   TEST_LSFit_D_Linc(_spoil_scenario,_TestResult,_TotalResult);          //82
   CHECK_RESULT
   TEST_LSFit_D_Pol(_spoil_scenario,_TestResult,_TotalResult);           //83
   CHECK_RESULT
   TEST_LSFit_D_Polc(_spoil_scenario,_TestResult,_TotalResult);          //84
   CHECK_RESULT
   TEST_LSFit_D_Spline(_spoil_scenario,_TestResult,_TotalResult);        //85
   CHECK_RESULT
   TEST_LSFit_T_PolFit_1(_spoil_scenario,_TestResult,_TotalResult);      //86
   CHECK_RESULT
   TEST_LSFit_T_PolFit_2(_spoil_scenario,_TestResult,_TotalResult);      //87
   CHECK_RESULT
   TEST_LSFit_T_PolFit_3(_spoil_scenario,_TestResult,_TotalResult);      //88
   CHECK_RESULT
   TEST_MatDet_D_1(_spoil_scenario,_TestResult,_TotalResult);            //89
   CHECK_RESULT
   TEST_MatDet_D_2(_spoil_scenario,_TestResult,_TotalResult);            //90
   CHECK_RESULT
   TEST_MatDet_D_3(_spoil_scenario,_TestResult,_TotalResult);            //91
   CHECK_RESULT
   TEST_MatDet_D_4(_spoil_scenario,_TestResult,_TotalResult);            //92
   CHECK_RESULT
   TEST_MatDet_D_5(_spoil_scenario,_TestResult,_TotalResult);            //93
   CHECK_RESULT
   TEST_MatDet_T_0(_spoil_scenario,_TestResult,_TotalResult);            //94
   CHECK_RESULT
   TEST_MatDet_T_1(_spoil_scenario,_TestResult,_TotalResult);            //95
   CHECK_RESULT
   TEST_MatDet_T_2(_spoil_scenario,_TestResult,_TotalResult);            //96
   CHECK_RESULT
   TEST_MatDet_T_3(_spoil_scenario,_TestResult,_TotalResult);            //97
   CHECK_RESULT
   TEST_MatDet_T_4(_spoil_scenario,_TestResult,_TotalResult);            //98
   CHECK_RESULT
   TEST_MatDet_T_5(_spoil_scenario,_TestResult,_TotalResult);            //99
   CHECK_RESULT
   TEST_MinQP_D_LC1(_spoil_scenario,_TestResult,_TotalResult);           //100
   CHECK_RESULT
//--- 100 blocks were successful
   Print("100/152");
   TEST_MinQP_D_U2(_spoil_scenario,_TestResult,_TotalResult);            //101
   CHECK_RESULT
   TEST_MinQP_D_NonConvex(_spoil_scenario,_TestResult,_TotalResult);     //102
   CHECK_RESULT
   TEST_MinLP_Basic(_spoil_scenario,_TestResult,_TotalResult);           //103
   CHECK_RESULT
   TEST_MinNLC_D_Inequality(_spoil_scenario,_TestResult,_TotalResult);   //104
   CHECK_RESULT
   TEST_MinBC_D_1(_spoil_scenario,_TestResult,_TotalResult);             //105
   CHECK_RESULT
   TEST_MinBC_NumDif(_spoil_scenario,_TestResult,_TotalResult);          //106
   CHECK_RESULT
   TEST_IDW_D_MSTAB(_spoil_scenario,_TestResult,_TotalResult);           //107
   CHECK_RESULT
   TEST_IDW_D_Serialize(_spoil_scenario,_TestResult,_TotalResult);       //108
   CHECK_RESULT
   TEST_Parametric_RDP(_spoil_scenario,_TestResult,_TotalResult);        //109
   CHECK_RESULT
   TEST_Spline2D_Bilinear(_spoil_scenario,_TestResult,_TotalResult);     //110
   CHECK_RESULT
   TEST_Spline2D_Bicubic(_spoil_scenario,_TestResult,_TotalResult);      //111
   CHECK_RESULT
   TEST_Spline2D_Unpack(_spoil_scenario,_TestResult,_TotalResult);       //112
   CHECK_RESULT
   TEST_Spline2D_CopyTrans(_spoil_scenario,_TestResult,_TotalResult);    //113
   CHECK_RESULT
   TEST_Spline3D_Trilinear(_spoil_scenario,_TestResult,_TotalResult);    //114
   CHECK_RESULT
   TEST_Spline3D_Vector(_spoil_scenario,_TestResult,_TotalResult);       //115
   CHECK_RESULT
   TEST_Clst_AHC(_spoil_scenario,_TestResult,_TotalResult);              //116
   CHECK_RESULT
   TEST_Clst_Linkage(_spoil_scenario,_TestResult,_TotalResult);          //117
   CHECK_RESULT
   TEST_Clst_Distance(_spoil_scenario,_TestResult,_TotalResult);         //118
   CHECK_RESULT
   TEST_Clst_KClusters(_spoil_scenario,_TestResult,_TotalResult);        //119
   CHECK_RESULT
   TEST_RandomForest_CLS(_spoil_scenario,_TestResult,_TotalResult);      //120
   CHECK_RESULT
   TEST_RandomForest_Reg(_spoil_scenario,_TestResult,_TotalResult);      //121
   CHECK_RESULT
   TEST_Filters_D_SMA(_spoil_scenario,_TestResult,_TotalResult);         //122
   CHECK_RESULT
   TEST_Filters_D_EMA(_spoil_scenario,_TestResult,_TotalResult);         //123
   CHECK_RESULT
   TEST_Filters_D_LRMA(_spoil_scenario,_TestResult,_TotalResult);        //124
   CHECK_RESULT
   TEST_SSA_D_Basic(_spoil_scenario,_TestResult,_TotalResult);           //125
   CHECK_RESULT
   TEST_SSA_D_Forecast(_spoil_scenario,_TestResult,_TotalResult);        //126
   CHECK_RESULT
   TEST_KNN_Reg(_spoil_scenario,_TestResult,_TotalResult);               //127
   CHECK_RESULT
   TEST_NN_Regr(_spoil_scenario,_TestResult,_TotalResult);               //128
   CHECK_RESULT
   TEST_NN_Regr_N(_spoil_scenario,_TestResult,_TotalResult);             //129
   CHECK_RESULT
   TEST_NN_Cls2(_spoil_scenario,_TestResult,_TotalResult);               //130
   CHECK_RESULT
   TEST_NN_Cls3(_spoil_scenario,_TestResult,_TotalResult);               //131
   CHECK_RESULT
   TEST_NN_TrainerObject(_spoil_scenario,_TestResult,_TotalResult);      //132
   CHECK_RESULT
   TEST_NN_Ensembles_ES(_spoil_scenario,_TestResult,_TotalResult);       //133
   CHECK_RESULT
   TEST_MinNLC_D_Equality(_spoil_scenario,_TestResult,_TotalResult);     //134
   CHECK_RESULT
   TEST_MinNLC_D_Mixed(_spoil_scenario,_TestResult,_TotalResult);        //135
   CHECK_RESULT
   TEST_MinNS_D_Unconstrained(_spoil_scenario,_TestResult,_TotalResult); //136
   CHECK_RESULT
   TEST_MinNS_D_Diff(_spoil_scenario,_TestResult,_TotalResult);          //137
   CHECK_RESULT
   TEST_MinNS_D_BC(_spoil_scenario,_TestResult,_TotalResult);            //138
   CHECK_RESULT
   TEST_MinNS_D_NLC(_spoil_scenario,_TestResult,_TotalResult);           //139
   CHECK_RESULT
   TEST_Spline2D_Fit_Blocklls(_spoil_scenario,_TestResult,_TotalResult); //140
   CHECK_RESULT
   TEST_Spline2d_Vector(_spoil_scenario,_TestResult,_TotalResult);       //141
   CHECK_RESULT
   TEST_RBF_D_HRBF(_spoil_scenario,_TestResult,_TotalResult);            //142
   CHECK_RESULT
   TEST_RBF_D_Vector(_spoil_scenario,_TestResult,_TotalResult);          //143
   CHECK_RESULT
   TEST_RBF_D_PolTerm(_spoil_scenario,_TestResult,_TotalResult);         //144
   CHECK_RESULT
   TEST_RBF_D_Serialize(_spoil_scenario,_TestResult,_TotalResult);       //145
   CHECK_RESULT
   TEST_Clst_KMeans(_spoil_scenario,_TestResult,_TotalResult);           //146
   CHECK_RESULT
   TEST_LinReg_D_Basic(_spoil_scenario,_TestResult,_TotalResult);        //147
   CHECK_RESULT
   TEST_SSA_D_Realtime(_spoil_scenario,_TestResult,_TotalResult);        //148
   CHECK_RESULT
   TEST_KNN_Cls(_spoil_scenario,_TestResult,_TotalResult);               //149
   CHECK_RESULT
   TEST_Spline1D_D_Monotone(_spoil_scenario,_TestResult,_TotalResult);   //150
   CHECK_RESULT
   TEST_LSFit_T_4pl(_spoil_scenario,_TestResult,_TotalResult);           //151
   CHECK_RESULT
   TEST_LSFit_T_5pl(_spoil_scenario,_TestResult,_TotalResult);           //152
   CHECK_RESULT
//--- all blocks were successful
   Print("152/152");
//--- print total result
   Print("Result = ",_TotalResult);
  }
//+------------------------------------------------------------------+
