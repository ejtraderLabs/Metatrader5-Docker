//+------------------------------------------------------------------+
//|                                                         VROC.mq5 |
//|                   Copyright 2009-2017, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2017, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Volume Rate of Change"
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 1
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  Green
#property indicator_style1  0
#property indicator_width1  1
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
//--- input parametrs
input int                 InpPeriodVROC=25;           // Period
input ENUM_APPLIED_VOLUME InpVolumeType=VOLUME_TICK;  // Volumes
//--- indicator buffer
double ExtVROCBuffer[];

int    ExtPeriodVROC;
//+------------------------------------------------------------------+
//| VROC initialization function                                     |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- check for input value
   if(InpPeriodVROC<=1)
     {
      ExtPeriodVROC=25;
      PrintFormat("Incorrect value for input variable InpPeriodVROC=%d. Indicator will use value=%d for calculations.",
                  InpPeriodVROC,ExtPeriodVROC);
     }
   else
      ExtPeriodVROC=InpPeriodVROC;
//--- define index buffer
   SetIndexBuffer(0,ExtVROCBuffer);
//--- set indicator short name
   IndicatorSetString(INDICATOR_SHORTNAME,"VROC("+string(ExtPeriodVROC)+")");
//--- set indicator digits
   IndicatorSetInteger(INDICATOR_DIGITS,2);
//--- set draw begin
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,ExtPeriodVROC-1);
  }
//+------------------------------------------------------------------+
//| Volume Rate of Change                                            |
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
   if(rates_total<ExtPeriodVROC)
      return(0);
//--- starting work
   int pos=prev_calculated-1;
   if(pos<ExtPeriodVROC-1)
     {
      pos=ExtPeriodVROC-1;
      for(int i=0; i<pos; i++)
         ExtVROCBuffer[i]=0.0;
     }
//--- main cycle by volume type
   if(InpVolumeType==VOLUME_TICK)
      CalculateVROC(pos,rates_total,tick_volume);
   else
      CalculateVROC(pos,rates_total,volume);
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
//| Calculate VROC by volume argument                                |
//+------------------------------------------------------------------+
void CalculateVROC(const int pos,const int rates_total,const long& volume[])
  {
   for(int i=pos; i<rates_total && !IsStopped(); i++)
     {
      double prev_volume=(double)(volume[i-(ExtPeriodVROC-1)]);
      double curr_volume=(double)volume[i];
      //--- calculate VROC
      if(prev_volume!=0.0)
         ExtVROCBuffer[i]=100.0*(curr_volume-prev_volume)/prev_volume;
      else
         ExtVROCBuffer[i]=ExtVROCBuffer[i-1];
     }
  }
//+------------------------------------------------------------------+
