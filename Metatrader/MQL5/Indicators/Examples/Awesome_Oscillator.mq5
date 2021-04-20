//+------------------------------------------------------------------+
//|                                           Awesome_Oscillator.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "2009-2020, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"

//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 4
#property indicator_plots   1
#property indicator_type1   DRAW_COLOR_HISTOGRAM
#property indicator_color1  Green,Red
#property indicator_width1  1
#property indicator_label1  "AO"
//--- indicator buffers
double ExtAOBuffer[];
double ExtColorBuffer[];
double ExtFastBuffer[];
double ExtSlowBuffer[];
//--- handles for MAs
int    ExtFastSMAHandle;
int    ExtSlowSMAHandle;

#define DATA_LIMIT  33
#define FAST_PERIOD 5
#define SLOW_PERIOD 34
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtAOBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtColorBuffer,INDICATOR_COLOR_INDEX);
   SetIndexBuffer(2,ExtFastBuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(3,ExtSlowBuffer,INDICATOR_CALCULATIONS);
//--- set accuracy
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits+1);
//--- sets first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,33);
//--- name for DataWindow
   IndicatorSetString(INDICATOR_SHORTNAME,"AO");
//--- get handles
   ExtFastSMAHandle=iMA(NULL,0,FAST_PERIOD,0,MODE_SMA,PRICE_MEDIAN);
   ExtSlowSMAHandle=iMA(NULL,0,SLOW_PERIOD,0,MODE_SMA,PRICE_MEDIAN);
  }
//+------------------------------------------------------------------+
//|  Awesome Oscillator                                              |
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
   if(rates_total<=DATA_LIMIT)
      return(0);
//--- not all data may be calculated
   int calculated=BarsCalculated(ExtFastSMAHandle);
   if(calculated<rates_total)
     {
      Print("Not all data of ExtFastSMAHandle is calculated (",calculated," bars). Error ",GetLastError());
      return(0);
     }
   calculated=BarsCalculated(ExtSlowSMAHandle);
   if(calculated<rates_total)
     {
      Print("Not all data of ExtSlowSMAHandle is calculated (",calculated," bars). Error ",GetLastError());
      return(0);
     }
//--- we can copy not all data
   int to_copy;
   if(prev_calculated>rates_total || prev_calculated<0)
      to_copy=rates_total;
   else
     {
      to_copy=rates_total-prev_calculated;
      if(prev_calculated>0)
         to_copy++;
     }
//--- get FastSMA buffer
   if(IsStopped()) // checking for stop flag
      return(0);
   if(CopyBuffer(ExtFastSMAHandle,0,0,to_copy,ExtFastBuffer)<=0)
     {
      Print("Getting fast SMA is failed! Error ",GetLastError());
      return(0);
     }
//--- get SlowSMA buffer
   if(IsStopped()) // checking for stop flag
      return(0);
   if(CopyBuffer(ExtSlowSMAHandle,0,0,to_copy,ExtSlowBuffer)<=0)
     {
      Print("Getting slow SMA is failed! Error ",GetLastError());
      return(0);
     }
//--- first calculation or number of bars was changed
   int i,start;
   if(prev_calculated<=DATA_LIMIT)
     {
      for(i=0; i<DATA_LIMIT; i++)
         ExtAOBuffer[i]=0.0;
      start=DATA_LIMIT;
     }
   else
      start=prev_calculated-1;
//--- main loop of calculations
   for(i=start; i<rates_total && !IsStopped(); i++)
     {
      ExtAOBuffer[i]=ExtFastBuffer[i]-ExtSlowBuffer[i];
      if(ExtAOBuffer[i]>ExtAOBuffer[i-1])
         ExtColorBuffer[i]=0.0; // set color Green
      else
         ExtColorBuffer[i]=1.0; // set color Red
     }
//--- return value of prev_calculated for next call
   return(rates_total);
  }

//+------------------------------------------------------------------+
