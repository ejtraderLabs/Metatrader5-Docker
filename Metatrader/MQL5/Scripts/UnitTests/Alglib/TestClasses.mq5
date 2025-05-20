//+------------------------------------------------------------------+
//|                                                  TestClasses.mq5 |
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
#define PrintElapsed(function,check,elapsed) PrintFormat("%-40s: %10s in %14s",function,(check ? "PASSED  " : "- FAILED -"),elapsed);

#include "TestClasses.mqh"

#ifndef _DEBUG
#property script_show_inputs
#endif

input bool InpSilent=true;   // Do not show extended log
input uint InpSeed=UINT_MAX; // Random seed
//+------------------------------------------------------------------+
//| Testing script                                                   |
//+------------------------------------------------------------------+
void OnStart()
  {
   bool  check;
   ulong start_mcs;
   ulong stop_mcs;
//--- initialization
   bool silent=InpSilent;
   uint seed=GetTickCount();
   if(InpSeed!=UINT_MAX)
      seed=InpSeed;
//--- seed
   PrintFormat("RandomSeed = %u",seed);
//--- start time
   datetime start_time=TimeLocal();

//--- check class CHighQualityRand                                            // 1
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestHQRndUnit::TestHQRnd(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CHighQualityRand",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CTSort                                                      // 2
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestTSortUnit::TestTSort(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CTSort",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CNearestNeighbor                                            // 3
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestNearestNeighborUnit::TestNearestNeighbor(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CNearestNeighbor",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CAblas                                                      // 4
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestAblasUnit::TestAblas(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CAblas",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CBaseStat                                                   // 5
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestBaseStatUnit::TestBaseStat(seed);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CBaseStat",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CBdSS                                                       // 6
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestBdSSUnit::TestBdSS(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CBdSS",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CDForest                                                    // 7
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestDForestUnit::TestDForest(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CDForest",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CBlas                                                       // 8
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestBlasUnit::TestBlas(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CBlas",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CKMeans                                                     // 9
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestKMeansUnit::TestKMeans(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CKMeans",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CHblas                                                      // 10
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestHblasUnit::TestHblas(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CHblas",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CReflections                                                // 11
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestReflectionsUnit::TestReflections(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CReflections",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CComplexReflections                                         // 12
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestCReflectionsUnit::TestCReflections(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CComplexReflections",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSblas                                                      // 13
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSblasUnit::TestSblas(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSblas",check,GetElapsed(stop_mcs-start_mcs));

//--- check class COrtFac                                                     // 14
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestOrtFacUnit::TestOrtFac(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("COrtFac",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CEigenVDetect                                               // 15
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestEVDUnit::TestEVD(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CEigenVDetect",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMatGen                                                     // 16
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMatGenUnit::TestMatGen(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMatGen",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CTrFac                                                      // 17
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestTrFacUnit::TestTrFac(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CTrFac",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CTrLinSolve                                                 // 18
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestTrLinSolveUnit::TestTrLinSolve(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CTrLinSolve",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSafeSolve                                                  // 19
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSafeSolveUnit::TestSafeSolve(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSafeSolve",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CRCond                                                      // 20
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestRCondUnit::TestRCond(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CRCond",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMatInv                                                     // 21
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMatInvUnit::TestMatInv(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMatInv",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CLDA                                                        // 22
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestLDAUnit::TestLDA(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CLDA",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CGammaFunc                                                  // 23
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestGammaFuncUnit::TestGammaFunc(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CGammaFunc",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CBdSingValueDecompose                                       // 24
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestBdSVDUnit::TestBdSVD(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CBdSingValueDecompose",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSingValueDecompose                                         // 25
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSVDUnit::TestSVD(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSingValueDecompose",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CLinReg                                                     // 26
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestLinRegUnit::TestLinReg(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CLinReg",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CXblas                                                      // 27
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestXBlasUnit::TestXBlas(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CXblas",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CDenseSolver                                                // 28
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestDenseSolverUnit::TestDenseSolver(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CDenseSolver",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinCG                                                      // 30
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinCGUnit::TestMinCG(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinCG",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinBLEIC                                                   // 31
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinBLEICUnit::TestMinBLEIC(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinBLEIC",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMarkovCPD                                                  // 32
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMCPDUnit::TestMCPD(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMarkovCPD",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CFbls                                                       // 33
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestFblsUnit::TestFbls(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CFbls",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinLBFGS                                                   // 34
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinLBFGSUnit::TestMinLBFGS(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinLBFGS",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMLPTrain                                                   // 35
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMLPTrainUnit::TestMLPTrain(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMLPTrain",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMLPE                                                       // 36
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMLPEUnit::TestMLPE(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMLPE",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CPCAnalysis                                                 // 37
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestPCAUnit::TestPCA(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CPCAnalysis",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CODESolver                                                  // 38
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestODESolverUnit::TestODESolver(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CODESolver",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CFastFourierTransform                                       // 39
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestFFTUnit::TestFFT(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CFastFourierTransform",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CConv                                                       // 40
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestConvUnit::TestConv(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CConv",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CCorr                                                       // 41
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestCorrUnit::TestCorr(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CCorr",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CFastHartleyTransform                                       // 42
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestFHTUnit::TestFHT(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CFastHartleyTransform",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CGaussQ                                                     // 43
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestGQUnit::TestGQ(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CGaussQ",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CGaussKronrodQ                                              // 44
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestGKQUnit::TestGKQ(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CGaussKronrodQ",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CAutoGK                                                     // 45
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestAutoGKUnit::TestAutoGK(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CAutoGK",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CIDWInt                                                     // 46
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestIDWIntUnit::TestIDWInt(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CIDWInt",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CRatInt                                                     // 47
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestRatIntUnit::TestRatInt(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CRatInt",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CPolInt                                                     // 48
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestPolIntUnit::TestPolInt(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CPolInt",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSpline1D                                                   // 49
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSpline1DUnit::TestSpline1D(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSpline1D",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinLM                                                      // 50
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinLMUnit::TestMinLM(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinLM",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CLSFit                                                      // 51
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestLSFitUnit::TestLSFit(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CLSFit",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CPSpline                                                    // 52
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestPSplineUnit::TestPSpline(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CPSpline",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSpline2D                                                   // 53
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSpline2DUnit::TestSpline2D(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSpline2D",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSpdGEVD                                                    // 54
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSpdGEVDUnit::TestSpdGEVD(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSpdGEVD",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CInverseUpdate                                              // 55
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestInverseUpdateUnit::TestInverseUpdate(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CInverseUpdate",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSchur                                                      // 56
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSchurUnit::TestSchur(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSchur",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CNlEq                                                       // 57
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestNlEqUnit::TestNlEq(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CNlEq",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CChebyshev                                                  // 58
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestChebyshevUnit::TestChebyshev(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CChebyshev",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CHermite                                                    // 59
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestHermiteUnit::TestHermite(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CHermite",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CLaguerre                                                   // 60
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestLaguerreUnit::TestLaguerre(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CLaguerre",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CLegendre                                                   // 61
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestLegendreUnit::TestLegendre(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CLegendre",check,GetElapsed(stop_mcs-start_mcs));

//--- check class AlglibBasics                                                // 62
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestAlglibBasicsUnit::TestAlglibBasics(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("AlglibBasics",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSparse                                                     // 63
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSparseUnit::TestSparse(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSparse",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CAblasF                                                     // 64
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestAblasFUnit::TestAblasF(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CAblasF",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CPolynomialSolver                                           // 65
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestPolynomialSolverUnit::TestPolynomialSolver(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CPolynomialSolver",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CDirectSparseSolvers                                        // 66
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestDirectSparseSolversUnit::TestDirectSparseSolvers(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CDirectSparseSolvers",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CIterativeSparse                                            // 67
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestIterativeSparseUnit::TestIterativeSparse(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CIterativeSparse",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CLinCG                                                      // 68
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestLinCGUnit::TestLinCG(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CLinCG",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CNormEstimator                                              // 69
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestNormEstimatorUnit::TestNormEstimator(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CNormEstimator",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CLinLSQR                                                    // 70
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestLinLSQRUnit::TestLinLSQR(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CLinLSQR",check,GetElapsed(stop_mcs-start_mcs));

//--- check class COptServ                                                    // 71
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestOptservUnit::TestOptserv(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("COptServ",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CCQModels                                                   // 72
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestCQModelsUnit::TestCQModels(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CCQModels",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSNNLS                                                      // 73
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSNNLSUnit::TestSNNLS(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSNNLS",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSActiveSets                                                // 74
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSActiveSetsUnit::TestSActiveSets(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSActiveSets",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinQP                                                      // 75
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinQPUnit::TestMinQP(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinQP",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinLP                                                      // 76
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinLPUnit::TestMinLP(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinLP",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinNLC                                                     // 77
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinNLCUnit::TestMinNLC(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinNLC",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinNS                                                      // 78
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinNSUnit::TestMinNS(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinNS",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMinBC                                                      // 79
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMinBCUnit::TestMinBC(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMinBC",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CNormalDistr                                                // 80
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestNormalDistrUnit::TestNormalDistr(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CNormalDistr",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CWilcoxonSignedRank                                         // 81
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestWSRUnit::TestWSR(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CWilcoxonSignedRank",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMannWhitneyU                                               // 82
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMannWhitneyUUnit::TestMannWhitneyU(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMannWhitneyU",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSignTest                                                   // 83
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSTestUnit::TestSTest(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSignTest",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CStudentTests                                               // 84
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestStudentTestsUnit::TestStudentTests(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CStudentTests",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CFitSphere                                                  // 85
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestFitSphereUnit::TestFitSphere(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CFitSphere",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSpline3D                                                   // 86
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSpline3DUnit::TestSpline3D(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSpline3D",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CRBF                                                        // 87
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestRBFUnit::TestRBF(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CRBF",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CMLPBase                                                    // 88
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestMLPBaseUnit::TestMLPBase(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CMLPBase",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CFilters                                                    // 89
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestFiltersUnit::TestFilters(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CFilters",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CSSA                                                        // 90
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestSSAUnit::TestSSA(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CSSA",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CKNN                                                        // 91
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestKNNUnit::TestKNN(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CKNN",check,GetElapsed(stop_mcs-start_mcs));

//--- check class CClustering                                                 // 92
   _RandomSeed=seed;
   start_mcs=GetMicrosecondCount();
   check=CTestClusteringUnit::TestClustering(silent);
   stop_mcs=GetMicrosecondCount();
   PrintElapsed("CClustering",check,GetElapsed(stop_mcs-start_mcs));

//--- finish time
   ulong elapsed_time=TimeLocal()-start_time;
   int   seconds=(int)elapsed_time%60;
   int   minutes=(int)elapsed_time/60;
   PrintFormat("Test passed in %d min %02d sec",minutes,seconds);
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
string GetElapsed(ulong microseconds)
  {
   int    sec=int(microseconds/1000000);
   int    mcs=int(microseconds%1000000);
   string elapsed=StringFormat("%d.%06d sec",sec,mcs);

   return(elapsed);
  }
//+------------------------------------------------------------------+
