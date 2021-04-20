//+------------------------------------------------------------------+
//|                                                         TRIX.mq5 |
//|                   Copyright 2009-2017, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2017, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Triple Exponential Average"
#include <MovingAverages.mqh>
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 4
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  Red
#property indicator_width1  1
#property indicator_label1  "TRIX"
#property indicator_applied_price PRICE_CLOSE
//--- input parameters
input int InpPeriodEMA=14;               // EMA period
//--- indicator buffers
double    TRIX_Buffer[];
double    EMA[];
double    SecondEMA[];
double    ThirdEMA[];
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,TRIX_Buffer,INDICATOR_DATA);
   SetIndexBuffer(1,EMA,INDICATOR_CALCULATIONS);
   SetIndexBuffer(2,SecondEMA,INDICATOR_CALCULATIONS);
   SetIndexBuffer(3,ThirdEMA,INDICATOR_CALCULATIONS);
//--- sets first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,3*InpPeriodEMA-3);
//--- name for index label
   string short_name=StringFormat("TRIX(%d)",InpPeriodEMA);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
   PlotIndexSetString(0,PLOT_LABEL,short_name);
//--- indicator digits
   IndicatorSetInteger(INDICATOR_DIGITS,5);
  }
//+------------------------------------------------------------------+
//| Triple Exponential Average                                       |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,
                const int prev_calculated,
                const int begin,
                const double &price[])
  {
   if(rates_total<3*InpPeriodEMA-3)
      return(0);
//---
   int start;
   if(prev_calculated==0)
     {
      start=3*(InpPeriodEMA-1);
      for(int i=0; i<start; i++)
         TRIX_Buffer[i]=EMPTY_VALUE;
     }
   else
      start=prev_calculated-1;
//--- calculate EMA
   ExponentialMAOnBuffer(rates_total,prev_calculated,0,InpPeriodEMA,price,EMA);
//--- calculate EMA on EMA array
   ExponentialMAOnBuffer(rates_total,prev_calculated,InpPeriodEMA-1,InpPeriodEMA,EMA,SecondEMA);
//--- calculate EMA on EMA array on EMA array
   ExponentialMAOnBuffer(rates_total,prev_calculated,2*InpPeriodEMA-2,InpPeriodEMA,SecondEMA,ThirdEMA);
//--- calculate TRIX
   for(int i=start; i<rates_total && !IsStopped(); i++)
     {
      if(ThirdEMA[i-1]!=0.0)
         TRIX_Buffer[i]=(ThirdEMA[i]-ThirdEMA[i-1])/ThirdEMA[i-1];
      else
         TRIX_Buffer[i]=0.0;
     }
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
