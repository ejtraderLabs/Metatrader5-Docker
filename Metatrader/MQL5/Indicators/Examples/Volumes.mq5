//+------------------------------------------------------------------+
//|                                                      Volumes.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "2009-2020, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 2
#property indicator_plots   1
#property indicator_type1   DRAW_COLOR_HISTOGRAM
#property indicator_color1  Green,Red
#property indicator_style1  0
#property indicator_width1  1
#property indicator_minimum 0.0
//--- input data
input ENUM_APPLIED_VOLUME InpVolumeType=VOLUME_TICK; // Volumes
//--- indicator buffers
double ExtVolumesBuffer[];
double ExtColorsBuffer[];
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- buffers
   SetIndexBuffer(0,ExtVolumesBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtColorsBuffer,INDICATOR_COLOR_INDEX);
//--- name for DataWindow and indicator subwindow label
   IndicatorSetString(INDICATOR_SHORTNAME,"Volumes");
//--- indicator digits
   IndicatorSetInteger(INDICATOR_DIGITS,0);
  }
//+------------------------------------------------------------------+
//|  Volumes                                                         |
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
   if(rates_total<2)
      return(0);
//--- starting work
   int pos=prev_calculated-1;
//--- correct position
   if(pos<1)
     {
      ExtVolumesBuffer[0]=0;
      pos=1;
     }
//--- main cycle
   if(InpVolumeType==VOLUME_TICK)
      CalculateVolume(pos,rates_total,tick_volume);
   else
      CalculateVolume(pos,rates_total,volume);
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CalculateVolume(const int pos,const int rates_total,const long& volume[])
  {
   ExtVolumesBuffer[0]=(double)volume[0];
   ExtColorsBuffer[0]=0.0;
//---
   for(int i=pos; i<rates_total && !IsStopped(); i++)
     {
      double curr_volume=(double)volume[i];
      double prev_volume=(double)volume[i-1];
      //--- calculate indicator
      ExtVolumesBuffer[i]=curr_volume;
      if(curr_volume>prev_volume)
         ExtColorsBuffer[i]=0.0;
      else
         ExtColorsBuffer[i]=1.0;
     }
//---
  }
//+------------------------------------------------------------------+
