//+------------------------------------------------------------------+
//|                                                        FrAMA.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2020, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Fractal Adaptive Moving Average"
//--- indicator settings
#property indicator_chart_window
#property indicator_buffers 1
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  DarkBlue
#property indicator_width1  1
#property indicator_label1  "FrAMA"
#property indicator_applied_price PRICE_CLOSE
//--- input parameters
input int InpPeriodFrAMA=14;            // FrAMA period
input int InpShift=0;                   // Indicator's shift
//--- indicator buffer
double    FrAmaBuffer[];
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,FrAmaBuffer,INDICATOR_DATA);
//--- sets first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,2*InpPeriodFrAMA-1);
//--- sets indicator shift
   PlotIndexSetInteger(0,PLOT_SHIFT,InpShift);
//--- name for labels
   string short_name=StringFormat("FrAMA(%d)",InpPeriodFrAMA);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
   PlotIndexSetString(0,PLOT_LABEL,short_name);
  }
//+------------------------------------------------------------------+
//| Fractal Adaptive Moving Average                                  |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,
                const int prev_calculated,
                const int begin,
                const double &price[])
  {
   if(rates_total<2*InpPeriodFrAMA)
      return(0);

   int start;
//--- start calculations
   if(prev_calculated==0)
     {
      start=2*InpPeriodFrAMA-1;
      for(int i=0; i<=start; i++)
         FrAmaBuffer[i]=price[i];
     }
   else
      start=prev_calculated-1;
//--- main cycle
   double math_log_2=MathLog(2.0);
   for(int i=start; i<rates_total && !IsStopped(); i++)
     {
      double hi1=iHigh(_Symbol,_Period,iHighest(_Symbol,_Period,MODE_HIGH,InpPeriodFrAMA,rates_total-i-1));
      double lo1=iLow(_Symbol,_Period,iLowest(_Symbol,_Period,MODE_LOW,InpPeriodFrAMA,rates_total-i-1));
      double hi2=iHigh(_Symbol,_Period,iHighest(_Symbol,_Period,MODE_HIGH,InpPeriodFrAMA,rates_total-i+InpPeriodFrAMA-1));
      double lo2=iLow(_Symbol,_Period,iLowest(_Symbol,_Period,MODE_LOW,InpPeriodFrAMA,rates_total-i+InpPeriodFrAMA-1));
      double hi3=iHigh(_Symbol,_Period,iHighest(_Symbol,_Period,MODE_HIGH,2*InpPeriodFrAMA,rates_total-i-1));
      double lo3=iLow(_Symbol,_Period,iLowest(_Symbol,_Period,MODE_LOW,2*InpPeriodFrAMA,rates_total-i-1));
      double n1=(hi1-lo1)/InpPeriodFrAMA;
      double n2=(hi2-lo2)/InpPeriodFrAMA;
      double n3=(hi3-lo3)/(2*InpPeriodFrAMA);
      double d=(MathLog(n1+n2)-MathLog(n3))/math_log_2;
      double alfa=MathExp(-4.6*(d-1.0));
      //---
      FrAmaBuffer[i]=alfa*price[i]+(1-alfa)*FrAmaBuffer[i-1];
     }
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
