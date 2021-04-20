//+------------------------------------------------------------------+
//|                                                          DPO.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2020, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Detrended Price Oscillator"
#include <MovingAverages.mqh>
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 2
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  DodgerBlue
//--- input parameters
input int InpDetrendPeriod=12; // Period
//--- indicator buffers
double    ExtDPOBuffer[];
double    ExtMABuffer[];

int       ExtMAPeriod;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- get length of cycle for smoothing
   ExtMAPeriod=InpDetrendPeriod/2+1;
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtDPOBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtMABuffer,INDICATOR_CALCULATIONS);
//--- set accuracy
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits+1);
//--- set first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,ExtMAPeriod-1);
//--- name for DataWindow and indicator subwindow label
   string short_name=StringFormat("DPO(%d)",InpDetrendPeriod);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
  }
//+------------------------------------------------------------------+
//| Detrended Price Oscillator                                       |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,
                const int prev_calculated,
                const int begin,
                const double &price[])
  {
   int start;
   int first_index=begin+ExtMAPeriod-1;
//--- preliminary filling
   if(prev_calculated<first_index)
     {
      ArrayInitialize(ExtDPOBuffer,0.0);
      start=first_index;
      if(begin>0)
         PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,first_index);
     }
   else
      start=prev_calculated-1;
//--- calculate simple moving average
   SimpleMAOnBuffer(rates_total,prev_calculated,begin,ExtMAPeriod,price,ExtMABuffer);
//--- the main loop of calculations
   for(int i=start; i<rates_total && !IsStopped(); i++)
      ExtDPOBuffer[i]=price[i]-ExtMABuffer[i];
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
