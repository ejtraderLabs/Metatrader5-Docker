//+------------------------------------------------------------------+
//|                                                          CHO.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2020, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Chaikin Oscillator"
#include <MovingAverages.mqh>
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 4
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  LightSeaGreen
//--- input parameters
input int                 InpFastMA=3;                // Fast MA period
input int                 InpSlowMA=10;               // Slow MA period
input ENUM_MA_METHOD      InpSmoothMethod=MODE_EMA;   // MA method
input ENUM_APPLIED_VOLUME InpVolumeType=VOLUME_TICK;  // Volumes
//--- indicator buffers
double ExtCHOBuffer[];
double ExtFastEMABuffer[];
double ExtSlowEMABuffer[];
double ExtADBuffer[];
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtCHOBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtFastEMABuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(2,ExtSlowEMABuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(3,ExtADBuffer,INDICATOR_CALCULATIONS);
//--- set accuracy
   IndicatorSetInteger(INDICATOR_DIGITS,0);
//--- sets first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,InpSlowMA);
//--- name for DataWindow and indicator subwindow label
   string short_name=StringFormat("CHO(%d,%d)",InpSlowMA,InpFastMA);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
  }
//+------------------------------------------------------------------+
//| Chaikin Oscillator                                               |
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
   if(rates_total<InpSlowMA)
      return(0);
//--- preliminary calculations
   int i,start;
   if(prev_calculated<2)
      start=0;
   else
      start=prev_calculated-2;
//--- calculate AD buffer
   if(InpVolumeType==VOLUME_TICK)
     {
      for(i=start; i<rates_total && !IsStopped(); i++)
        {
         ExtADBuffer[i]=AD(high[i],low[i],close[i],tick_volume[i]);
         if(i>0)
            ExtADBuffer[i]+=ExtADBuffer[i-1];
        }
     }
   else
     {
      for(i=start; i<rates_total && !IsStopped(); i++)
        {
         ExtADBuffer[i]=AD(high[i],low[i],close[i],volume[i]);
         if(i>0)
            ExtADBuffer[i]+=ExtADBuffer[i-1];
        }
     }
//--- calculate EMA on array ExtADBuffer
   AverageOnArray(InpSmoothMethod,rates_total,prev_calculated,0,InpFastMA,ExtADBuffer,ExtFastEMABuffer);
   AverageOnArray(InpSmoothMethod,rates_total,prev_calculated,0,InpSlowMA,ExtADBuffer,ExtSlowEMABuffer);
//--- calculate chaikin oscillator
   for(i=start; i<rates_total && !IsStopped(); i++)
      ExtCHOBuffer[i]=ExtFastEMABuffer[i]-ExtSlowEMABuffer[i];
//--- return value of prev_calculated for next call
   return(rates_total);
  }
//+------------------------------------------------------------------+
//| calculate AD                                                     |
//+------------------------------------------------------------------+
double AD(double high,double low,double close,long volume)
  {
   double res=0.0;
//---
   double sum=(close-low)-(high-close);
   if(sum!=0.0)
     {
      if(high!=low)
         res=(sum/(high-low))*volume;
     }
//---
   return(res);
  }
//+------------------------------------------------------------------+
//| calculate average on array                                       |
//+------------------------------------------------------------------+
void AverageOnArray(const int mode,const int rates_total,const int prev_calculated,const int begin,
                    const int period,const double& source[],double& destination[])
  {
   switch(mode)
     {
      case MODE_EMA:
         ExponentialMAOnBuffer(rates_total,prev_calculated,begin,period,source,destination);
         break;
      case MODE_SMMA:
         SmoothedMAOnBuffer(rates_total,prev_calculated,begin,period,source,destination);
         break;
      case MODE_LWMA:
         LinearWeightedMAOnBuffer(rates_total,prev_calculated,begin,period,source,destination);
         break;
      default:
         SimpleMAOnBuffer(rates_total,prev_calculated,begin,period,source,destination);
     }
  }
//+------------------------------------------------------------------+
