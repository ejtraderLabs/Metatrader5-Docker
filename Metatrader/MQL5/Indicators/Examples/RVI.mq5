//+------------------------------------------------------------------+
//|                                                          RVI.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright   "2009-2020, MetaQuotes Software Corp."
#property link        "http://www.mql5.com"
#property description "Relative Vigor Index"
//--- indicator settings
#property indicator_separate_window
#property indicator_buffers 2
#property indicator_plots   2
#property indicator_type1   DRAW_LINE
#property indicator_type2   DRAW_LINE
#property indicator_color1  Green
#property indicator_color2  Red
#property indicator_label1  "RVI"
#property indicator_label2  "Signal"
//--- input parameters
input int InpRVIPeriod=10; // Period
//--- indicator buffers
double    ExtRVIBuffer[];
double    ExtSignalBuffer[];

#define TRIANGLE_PERIOD  3
#define AVERAGE_PERIOD   (TRIANGLE_PERIOD*2)
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtRVIBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtSignalBuffer,INDICATOR_DATA);
   IndicatorSetInteger(INDICATOR_DIGITS,3);
//--- sets first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,(InpRVIPeriod-1)+TRIANGLE_PERIOD);
   PlotIndexSetInteger(1,PLOT_DRAW_BEGIN,(InpRVIPeriod-1)+AVERAGE_PERIOD);
//--- name for DataWindow and indicator subwindow label
   string short_name=StringFormat("RVI(%d)",InpRVIPeriod);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
   PlotIndexSetString(0,PLOT_LABEL,short_name);
   PlotIndexSetString(1,PLOT_LABEL,"Signal("+string(InpRVIPeriod)+")");
  }
//+------------------------------------------------------------------+
//| Relative Vigor Index                                             |
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
   if(rates_total<=InpRVIPeriod+AVERAGE_PERIOD+2)
      return(0);

   int i,start;
//--- last counted bar will be recounted
   start=InpRVIPeriod+2;
   if(prev_calculated>InpRVIPeriod+TRIANGLE_PERIOD+2)
      start=prev_calculated-1;
//--- set empty value for uncalculated bars
   if(prev_calculated==0)
     {
      for(i=0; i<InpRVIPeriod+TRIANGLE_PERIOD; i++)
         ExtRVIBuffer[i]=0.0;
      for(i=0; i<InpRVIPeriod+AVERAGE_PERIOD; i++)
         ExtSignalBuffer[i]=0.0;
     }
//--- RVI counted in the 1-st buffer
   for(i=start; i<rates_total && !IsStopped(); i++)
     {
      double sum_up=0.0;
      double sum_down=0.0;
      for(int j=i; j>i-InpRVIPeriod; j--)
        {
         double value_up=close[j]-open[j]+2*(close[j-1]-open[j-1])+2*(close[j-2]-open[j-2])+close[j-3]-open[j-3];
         double value_down=high[j]-low[j]+2*(high[j-1]-low[j-1])+2*(high[j-2]-low[j-2])+high[j-3]-low[j-3];
         sum_up+=value_up;
         sum_down+=value_down;
        }
      if(sum_down!=0.0)
         ExtRVIBuffer[i]=sum_up/sum_down;
      else
         ExtRVIBuffer[i]=sum_up;
     }
//--- signal line counted in the 2-nd buffer
   start=InpRVIPeriod+TRIANGLE_PERIOD+2;
   if(prev_calculated>InpRVIPeriod+AVERAGE_PERIOD+2)
      start=prev_calculated-1;
   for(i=start; i<rates_total && !IsStopped(); i++)
      ExtSignalBuffer[i]=(ExtRVIBuffer[i]+2*ExtRVIBuffer[i-1]+2*ExtRVIBuffer[i-2]+ExtRVIBuffer[i-3])/AVERAGE_PERIOD;

//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
