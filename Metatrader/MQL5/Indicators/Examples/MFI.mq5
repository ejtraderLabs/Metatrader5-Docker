//+------------------------------------------------------------------+
//|                                                          MFI.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2020, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Money Flow Index"
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers    1
#property indicator_plots      1
#property indicator_type1      DRAW_LINE
#property indicator_color1     DodgerBlue
#property indicator_maximum    100.0
#property indicator_minimum    0.0
#property indicator_level1     20.0
#property indicator_level2     80.0
#property indicator_levelcolor Silver
#property indicator_levelstyle 2
#property indicator_levelwidth 1
//--- input parameters
input int                 InpMFIPeriod=14;            // Period
input ENUM_APPLIED_VOLUME InpVolumeType=VOLUME_TICK;  // Volumes
//--- indicator buffer
double ExtMFIBuffer[];

int    ExtMFIPeriod;
//+------------------------------------------------------------------+
//| Money Flow Index initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- check for input value
   if(InpMFIPeriod<=0)
     {
      ExtMFIPeriod=14;
      Print("Parameter InpMFIPeriod has wrong value. Indicator will use value ",ExtMFIPeriod);
     }
   else
      ExtMFIPeriod=InpMFIPeriod;
//--- indicator buffer
   SetIndexBuffer(0,ExtMFIBuffer);
//--- name for DataWindow and indicator subwindow label
   IndicatorSetString(INDICATOR_SHORTNAME,"MFI("+string(ExtMFIPeriod)+")");
//--- set draw begin
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,ExtMFIPeriod);
//--- set indicator digits
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits);
  }
//+------------------------------------------------------------------+
//| Money Flow Index                                                 |
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
   if(rates_total<ExtMFIPeriod)
      return(0);

   int start_position;
//--- start working
   if(prev_calculated<ExtMFIPeriod)
      start_position=ExtMFIPeriod;
   else
      start_position=prev_calculated-1;
//--- calculate MFI by volume
   if(InpVolumeType==VOLUME_TICK)
      CalculateMFI(start_position,rates_total,high,low,close,tick_volume);
   else
      CalculateMFI(start_position,rates_total,high,low,close,volume);
//--- OnCalculate done. Return new prev_calculated
   return(rates_total);
  }
//+------------------------------------------------------------------+
//| Calculate MFI by volume from argument                            |
//+------------------------------------------------------------------+
void CalculateMFI(const int start_position,
                  const int rates_total,
                  const double &high[],
                  const double &low[],
                  const double &close[],
                  const long &volume[])
  {
   for(int i=start_position; i<rates_total && !IsStopped(); i++)
     {
      double positive=0.0;
      double negative=0.0;
      double current_tp=TypicalPrice(high[i],low[i],close[i]);
      for(int j=1; j<=ExtMFIPeriod; j++)
        {
         int    index=i-j;
         double previous_tp=TypicalPrice(high[index],low[index],close[index]);
         if(current_tp>previous_tp)
            positive+=volume[index+1]*current_tp;
         if(current_tp<previous_tp)
            negative+=volume[index+1]*current_tp;
         current_tp=previous_tp;
        }
      if(negative!=0.0)
         ExtMFIBuffer[i]=100.0-100.0/(1+positive/negative);
      else
         ExtMFIBuffer[i]=100.0;
     }
  }
//+------------------------------------------------------------------+
//| Calculate typical price                                          |
//+------------------------------------------------------------------+
double TypicalPrice(const double high_price,const double low_price,const double close_price)
  {
   return((high_price+low_price+close_price)/3.0);
  }
//+------------------------------------------------------------------+
