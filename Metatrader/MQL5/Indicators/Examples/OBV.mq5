//+------------------------------------------------------------------+
//|                                                          OBV.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2020, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "On Balance vol"
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 1
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  DodgerBlue
#property indicator_label1  "OBV"
//--- input parametrs
input ENUM_APPLIED_VOLUME InpVolumeType=VOLUME_TICK; // Volumes
//--- indicator buffer
double ExtOBVBuffer[];
//+------------------------------------------------------------------+
//| On Balance vol initialization function                        |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- define indicator buffer
   SetIndexBuffer(0,ExtOBVBuffer);
//--- set indicator digits
   IndicatorSetInteger(INDICATOR_DIGITS,0);
  }
//+------------------------------------------------------------------+
//| On Balance vol                                                |
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
//--- starting calculation
   int pos=prev_calculated-1;
//--- correct position, when it's first iteration
   if(pos<1)
     {
      pos=1;
      if(InpVolumeType==VOLUME_TICK)
         ExtOBVBuffer[0]=(double)tick_volume[0];
      else
         ExtOBVBuffer[0]=(double)volume[0];
     }
//--- main cycle
   if(InpVolumeType==VOLUME_TICK)
      CalculateOBV(pos,rates_total,close,tick_volume);
   else
      CalculateOBV(pos,rates_total,close,volume);
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
//| Calculate OBV by volume argument                                 |
//+------------------------------------------------------------------+
void CalculateOBV(int start_pos,
                  int rates_total,
                  const double& close[],
                  const long& volume[])
  {
   for(int i=start_pos; i<rates_total && !IsStopped(); i++)
     {
      double vol=(double)volume[i];
      double prev_close=close[i-1];
      double curr_close=close[i];
      //--- fill ExtOBVBuffer
      if(curr_close<prev_close)
         ExtOBVBuffer[i]=ExtOBVBuffer[i-1]-vol;
      else
        {
         if(curr_close>prev_close)
            ExtOBVBuffer[i]=ExtOBVBuffer[i-1]+vol;
         else
            ExtOBVBuffer[i]=ExtOBVBuffer[i-1];
        }
     }
  }
//+------------------------------------------------------------------+
