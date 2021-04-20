//+------------------------------------------------------------------+
//|                                               Price_Channell.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "2009-2020, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
//--- indicator settings
#property indicator_chart_window
#property indicator_buffers 3
#property indicator_plots   2
#property indicator_type1   DRAW_FILLING
#property indicator_type2   DRAW_LINE
#property indicator_color1  DodgerBlue,Gray
#property indicator_color2  Blue
#property indicator_label1  "Channel upper;Channel lower"
#property indicator_label2  "Channel median"
//--- input parameters
input int InpChannelPeriod=22; // Period
//--- indicator buffers
double    ExtHighBuffer[];
double    ExtLowBuffer[];
double    ExtMiddBuffer[];
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtHighBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtLowBuffer,INDICATOR_DATA);
   SetIndexBuffer(2,ExtMiddBuffer,INDICATOR_DATA);
//--- set accuracy
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits);
//--- set first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,InpChannelPeriod);
//--- line shifts when drawing
   PlotIndexSetInteger(0,PLOT_SHIFT,1);
   PlotIndexSetInteger(1,PLOT_SHIFT,1);
//--- name for DataWindow and indicator label
   string short_name=StringFormat("Price Channel(%d)",InpChannelPeriod);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
   short_name=StringFormat("Channel(%d) upper;Channel(%d) lower",InpChannelPeriod,InpChannelPeriod);
   PlotIndexSetString(0,PLOT_LABEL,short_name);
   short_name=StringFormat("Median(%d)",InpChannelPeriod);
   PlotIndexSetString(1,PLOT_LABEL,short_name);
//--- set drawing line empty value
   PlotIndexSetDouble(0,PLOT_EMPTY_VALUE,0.0);
   PlotIndexSetDouble(1,PLOT_EMPTY_VALUE,0.0);
  }
//+------------------------------------------------------------------+
//| Price Channell                                                   |
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
   if(rates_total<InpChannelPeriod)
      return(0);

   int start;
//--- preliminary calculations
   if(prev_calculated==0)
      start=InpChannelPeriod;
   else
      start=prev_calculated-1;
//--- the main loop of calculations
   for(int i=start; i<rates_total && !IsStopped(); i++)
     {
      ExtHighBuffer[i]=Highest(high,InpChannelPeriod,i);
      ExtLowBuffer[i]=Lowest(low,InpChannelPeriod,i);
      ExtMiddBuffer[i]=(ExtHighBuffer[i]+ExtLowBuffer[i])/2.0;;
     }
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
//| get highest value for range                                      |
//+------------------------------------------------------------------+
double Highest(const double& array[],const int range,const int from_index)
  {
   double res=array[from_index];
   for(int i=from_index-1; i>from_index-range && i>=0; i--)
      if(res<array[i])
         res=array[i];
   return(res);
  }
//+------------------------------------------------------------------+
//| get lowest value for range                                       |
//+------------------------------------------------------------------+
double Lowest(const double& array[],const int range,const int from_index)
  {
   double res=array[from_index];
   for(int i=from_index-1; i>from_index-range && i>=0; i--)
      if(res>array[i])
         res=array[i];
   return(res);
  }
//+------------------------------------------------------------------+
