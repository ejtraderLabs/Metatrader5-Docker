//+------------------------------------------------------------------+
//|                                                          CHV.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2020, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Chaikin Volatility"
#include <MovingAverages.mqh>
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 3
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  DodgerBlue
//--- enum
enum SmoothMethod
  {
   SMA=0,// Simple MA
   EMA=1 // Exponential MA
  };
//--- input parameters
input int          InpSmoothPeriod=10;  // Smoothing period
input int          InpCHVPeriod=10;     // CHV period
input SmoothMethod InpSmoothType=EMA;   // Smoothing method
//--- indicator buffers
double ExtCHVBuffer[];
double ExtHLBuffer[];
double ExtSHLBuffer[];

int    ExtSmoothPeriod,ExtCHVPeriod;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
   string ma_name=EnumToString(InpSmoothType);
//--- check inputs
   if(InpSmoothPeriod<=0)
     {
      ExtSmoothPeriod=10;
      PrintFormat("Incorrect value for input variable InpSmoothPeriod=%d. Indicator will use value=%d for calculations.",InpSmoothPeriod,ExtSmoothPeriod);
     }
   else
      ExtSmoothPeriod=InpSmoothPeriod;
   if(InpCHVPeriod<=0)
     {
      ExtCHVPeriod=10;
      PrintFormat("Incorrect value for input variable InpCHVPeriod=%d. Indicator will use value=%d for calculations.",InpCHVPeriod,ExtCHVPeriod);
     }
   else
      ExtCHVPeriod=InpCHVPeriod;
//--- define buffers
   SetIndexBuffer(0,ExtCHVBuffer);
   SetIndexBuffer(1,ExtHLBuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(2,ExtSHLBuffer,INDICATOR_CALCULATIONS);
//--- set draw begin
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,ExtSmoothPeriod+ExtCHVPeriod-1);
//--- set name and index label
   string params=StringFormat("(%d,%s)",ExtSmoothPeriod,ma_name);
   IndicatorSetString(INDICATOR_SHORTNAME,"Chaikin Volatility"+params);
   PlotIndexSetString(0,PLOT_LABEL,"CHV"+params);
//--- round settings
   IndicatorSetInteger(INDICATOR_DIGITS,1);
  }
//+------------------------------------------------------------------+
//| Custom indicator iteration function                              |
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
   int i,pos,pos_chv;
//--- check for rates total
   pos_chv=ExtCHVPeriod+ExtSmoothPeriod-2;
   if(rates_total<pos_chv)
      return(0);
//--- start working
   if(prev_calculated<1)
      pos=0;
   else
      pos=prev_calculated-1;
//--- fill H-L(i) buffer
   for(i=pos; i<rates_total && !IsStopped(); i++)
      ExtHLBuffer[i]=high[i]-low[i];
//--- calculate smoothed H-L(i) buffer
   if(pos<ExtSmoothPeriod-1)
     {
      pos=ExtSmoothPeriod-1;
      for(i=0; i<pos; i++)
         ExtSHLBuffer[i]=0.0;
     }
   if(InpSmoothType==SMA)
      SimpleMAOnBuffer(rates_total,prev_calculated,0,ExtSmoothPeriod,ExtHLBuffer,ExtSHLBuffer);
   else
      ExponentialMAOnBuffer(rates_total,prev_calculated,0,ExtSmoothPeriod,ExtHLBuffer,ExtSHLBuffer);
//--- correct calc position
   if(pos<pos_chv)
      pos=pos_chv;
//--- calculate CHV buffer
   for(i=pos; i<rates_total && !IsStopped(); i++)
     {
      if(ExtSHLBuffer[i-ExtCHVPeriod]!=0.0)
         ExtCHVBuffer[i]=100.0*(ExtSHLBuffer[i]-ExtSHLBuffer[i-ExtCHVPeriod])/ExtSHLBuffer[i-ExtCHVPeriod];
      else
         ExtCHVBuffer[i]=0.0;
     }
//---
   return(rates_total);
  }
//+------------------------------------------------------------------+
