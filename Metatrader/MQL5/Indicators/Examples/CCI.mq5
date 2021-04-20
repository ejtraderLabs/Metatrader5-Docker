//+------------------------------------------------------------------+
//|                                                          CCI.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2020, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Commodity Channel Index"
#include <MovingAverages.mqh>
//---
#property indicator_separate_window
#property indicator_buffers       4
#property indicator_plots         1
#property indicator_type1         DRAW_LINE
#property indicator_color1        LightSeaGreen
#property indicator_level1       -100.0
#property indicator_level2        100.0
#property indicator_applied_price PRICE_TYPICAL
//--- input parametrs
input int  InpCCIPeriod=14; // Period
//--- indicator buffers
double     ExtSPBuffer[];
double     ExtDBuffer[];
double     ExtMBuffer[];
double     ExtCCIBuffer[];

int        ExtCCIPeriod;
double     ExtMultiplyer;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- check for input value of period
   if(InpCCIPeriod<=0)
     {
      ExtCCIPeriod=14;
      PrintFormat("Incorrect value for input variable InpCCIPeriod=%d. Indicator will use value=%d for calculations.",InpCCIPeriod,ExtCCIPeriod);
     }
   else
      ExtCCIPeriod=InpCCIPeriod;
   ExtMultiplyer=0.015/ExtCCIPeriod;
//--- define buffers
   SetIndexBuffer(0,ExtCCIBuffer);
   SetIndexBuffer(1,ExtDBuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(2,ExtMBuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(3,ExtSPBuffer,INDICATOR_CALCULATIONS);
//--- indicator name
   string short_name=StringFormat("CCI(%d)",ExtCCIPeriod);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
//--- indexes draw begin settings
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,ExtCCIPeriod-1);
//--- number of digits of indicator value
   IndicatorSetInteger(INDICATOR_DIGITS,2);
  }
//+------------------------------------------------------------------+
//| Custom indicator iteration function                              |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,
                const int prev_calculated,
                const int begin,
                const double &price[])
  {
   int start=(ExtCCIPeriod-1)+begin;
   if(rates_total<start)
      return(0);
//--- correct draw begin
   if(begin>0)
      PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,start+(ExtCCIPeriod-1));
//--- calculate position
   int pos=prev_calculated-1;
   if(pos<start)
      pos=start;
//--- main cycle
   for(int i=pos; i<rates_total && !IsStopped(); i++)
     {
      ExtSPBuffer[i]=SimpleMA(i,ExtCCIPeriod,price);
      //--- calculate D
      double tmp_d=0.0;
      for(int j=0; j<ExtCCIPeriod; j++)
         tmp_d+=MathAbs(price[i-j]-ExtSPBuffer[i]);
      ExtDBuffer[i]=tmp_d*ExtMultiplyer;
      //--- calculate M
      ExtMBuffer[i]=price[i]-ExtSPBuffer[i];
      //--- calculate CCI
      if(ExtDBuffer[i]!=0.0)
         ExtCCIBuffer[i]=ExtMBuffer[i]/ExtDBuffer[i];
      else
         ExtCCIBuffer[i]=0.0;
     }
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
