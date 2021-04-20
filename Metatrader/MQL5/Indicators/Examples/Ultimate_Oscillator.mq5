//+------------------------------------------------------------------+
//|                                          Ultimate_Oscillator.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "2009-2020, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
#include <MovingAverages.mqh>
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 5
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  DodgerBlue
//--- input parameters
input int InpFastPeriod=7;     // Fast ATR period
input int InpMiddlePeriod=14;  // Middle ATR period
input int InpSlowPeriod=28;    // Slow ATR period
input int InpFastK=4;          // Fast K
input int InpMiddleK=2;        // Middle K
input int InpSlowK=1;          // Slow K
//--- indicator buffers
double    ExtUOBuffer[];
double    ExtBPBuffer[];
double    ExtFastATRBuffer[];
double    ExtMiddleATRBuffer[];
double    ExtSlowATRBuffer[];
//--- indicator handles
int       ExtFastATRhandle;
int       ExtMiddleATRhandle;
int       ExtSlowATRhandle;

double    ExtDivider;
int       ExtMaxPeriod;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtUOBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtBPBuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(2,ExtFastATRBuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(3,ExtMiddleATRBuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(4,ExtSlowATRBuffer,INDICATOR_CALCULATIONS);
//--- set accuracy
   IndicatorSetInteger(INDICATOR_DIGITS,2);
//--- set levels
   IndicatorSetInteger(INDICATOR_LEVELS,2);
   IndicatorSetDouble(INDICATOR_LEVELVALUE,0,30);
   IndicatorSetDouble(INDICATOR_LEVELVALUE,1,70);
//--- set maximum and minimum for subwindow
   IndicatorSetDouble(INDICATOR_MINIMUM,0);
   IndicatorSetDouble(INDICATOR_MAXIMUM,100);
//--- set first bar from which index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,InpSlowPeriod-1);
//--- name for DataWindow and indicator subwindow label
   string short_name=StringFormat("UOS(%d,%d,%d)",InpFastPeriod,InpMiddlePeriod,InpSlowPeriod);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
//--- get handles
   ExtFastATRhandle=iATR(Symbol(),0,InpFastPeriod);
   ExtMiddleATRhandle=iATR(Symbol(),0,InpMiddlePeriod);
   ExtSlowATRhandle=iATR(Symbol(),0,InpSlowPeriod);
//---
   ExtDivider=InpFastK+InpMiddleK+InpSlowK;
   ExtMaxPeriod=InpSlowPeriod;
   if(ExtMaxPeriod<InpMiddlePeriod)
      ExtMaxPeriod=InpMiddlePeriod;
   if(ExtMaxPeriod<InpFastPeriod)
      ExtMaxPeriod=InpFastPeriod;
  }
//+------------------------------------------------------------------+
//| Ultimate Oscillator                                              |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,
                const int prev_calculated,
                const datetime &time[],
                const double &open[],
                const double &high[],
                const double &low[],
                const double &close[],
                const long &tick_volume[],
                const long &volume[],
                const int &spread[])
  {
   if(rates_total<ExtMaxPeriod)
      return(0);
//--- not all data may be calculated
   int calculated=BarsCalculated(ExtFastATRhandle);
   if(calculated<rates_total)
     {
      Print("Not all data of ExtFastATRhandle is calculated (",calculated," bars). Error ",GetLastError());
      return(0);
     }
   calculated=BarsCalculated(ExtMiddleATRhandle);
   if(calculated<rates_total)
     {
      Print("Not all data of ExtMiddleATRhandle is calculated (",calculated," bars). Error ",GetLastError());
      return(0);
     }
   calculated=BarsCalculated(ExtSlowATRhandle);
   if(calculated<rates_total)
     {
      Print("Not all data of ExtSlowATRhandle is calculated (",calculated," bars). Error ",GetLastError());
      return(0);
     }
//--- we can copy not all data
   int to_copy;
   if(prev_calculated>rates_total || prev_calculated<0)
      to_copy=rates_total;
   else
     {
      to_copy=rates_total-prev_calculated;
      if(prev_calculated>0)
         to_copy++;
     }
//--- get ATR buffers
   if(IsStopped()) // checking for stop flag
      return(0);
   if(CopyBuffer(ExtFastATRhandle,0,0,to_copy,ExtFastATRBuffer)<=0)
     {
      Print("getting ExtFastATRhandle is failed! Error ",GetLastError());
      return(0);
     }
   if(IsStopped()) // checking for stop flag
      return(0);
   if(CopyBuffer(ExtMiddleATRhandle,0,0,to_copy,ExtMiddleATRBuffer)<=0)
     {
      Print("getting ExtMiddleATRhandle is failed! Error ",GetLastError());
      return(0);
     }
   if(IsStopped()) // checking for stop flag
      return(0);
   if(CopyBuffer(ExtSlowATRhandle,0,0,to_copy,ExtSlowATRBuffer)<=0)
     {
      Print("getting ExtSlowATRhandle is failed! Error ",GetLastError());
      return(0);
     }
//--- preliminary calculations
   int i,start;
   if(prev_calculated==0)
     {
      ExtBPBuffer[0]=0.0;
      ExtUOBuffer[0]=0.0;
      //--- set value for first InpSlowPeriod bars
      for(i=1; i<=InpSlowPeriod; i++)
        {
         ExtUOBuffer[i]=0.0;
         double true_low=MathMin(low[i],close[i-1]);
         ExtBPBuffer[i]=close[i]-true_low;
        }
      //--- now we are going to calculate from start index in main loop
      start=InpSlowPeriod+1;
     }
   else
      start=prev_calculated-1;
//--- the main loop of calculations
   for(i=start; i<rates_total && !IsStopped(); i++)
     {
      double true_low=MathMin(low[i],close[i-1]);
      ExtBPBuffer[i]=close[i]-true_low;           // buying pressure

      if(ExtFastATRBuffer[i]!=0.0 &&
         ExtMiddleATRBuffer[i]!=0.0 &&
         ExtSlowATRBuffer[i]!=0.0)
        {
         double raw_uo=InpFastK*SimpleMA(i,InpFastPeriod,ExtBPBuffer)/ExtFastATRBuffer[i]+
                       InpMiddleK*SimpleMA(i,InpMiddlePeriod,ExtBPBuffer)/ExtMiddleATRBuffer[i]+
                       InpSlowK*SimpleMA(i,InpSlowPeriod,ExtBPBuffer)/ExtSlowATRBuffer[i];
         ExtUOBuffer[i]=raw_uo/ExtDivider*100;
        }
      else
         ExtUOBuffer[i]=ExtUOBuffer[i-1]; // set current Ultimate value as previous Ultimate value
     }
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
