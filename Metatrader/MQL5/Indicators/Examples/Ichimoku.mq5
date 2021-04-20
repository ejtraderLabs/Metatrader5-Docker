//+------------------------------------------------------------------+
//|                                                     Ichimoku.mq5 |
//|                   Copyright 2009-2017, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "2009-2017, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
#property description "Ichimoku Kinko Hyo"
//--- indicator settings
#property indicator_chart_window
#property indicator_buffers 5
#property indicator_plots   4
#property indicator_type1   DRAW_LINE
#property indicator_type2   DRAW_LINE
#property indicator_type3   DRAW_FILLING
#property indicator_type4   DRAW_LINE
#property indicator_color1  Red
#property indicator_color2  Blue
#property indicator_color3  SandyBrown,Thistle
#property indicator_color4  Lime
#property indicator_label1  "Tenkan-sen"
#property indicator_label2  "Kijun-sen"
#property indicator_label3  "Senkou Span A;Senkou Span B"
#property indicator_label4  "Chikou Span"
//--- input parameters
input int InpTenkan=9;     // Tenkan-sen
input int InpKijun=26;     // Kijun-sen
input int InpSenkou=52;    // Senkou Span B
//--- indicator buffers
double    ExtTenkanBuffer[];
double    ExtKijunBuffer[];
double    ExtSpanABuffer[];
double    ExtSpanBBuffer[];
double    ExtChikouBuffer[];
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtTenkanBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtKijunBuffer,INDICATOR_DATA);
   SetIndexBuffer(2,ExtSpanABuffer,INDICATOR_DATA);
   SetIndexBuffer(3,ExtSpanBBuffer,INDICATOR_DATA);
   SetIndexBuffer(4,ExtChikouBuffer,INDICATOR_DATA);
//---
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits+1);
//--- sets first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,InpTenkan);
   PlotIndexSetInteger(1,PLOT_DRAW_BEGIN,InpKijun);
   PlotIndexSetInteger(2,PLOT_DRAW_BEGIN,InpSenkou-1);
//--- lines shifts when drawing
   PlotIndexSetInteger(2,PLOT_SHIFT,InpKijun);
   PlotIndexSetInteger(3,PLOT_SHIFT,-InpKijun);
//--- change labels for DataWindow
   PlotIndexSetString(0,PLOT_LABEL,"Tenkan-sen("+string(InpTenkan)+")");
   PlotIndexSetString(1,PLOT_LABEL,"Kijun-sen("+string(InpKijun)+")");
   PlotIndexSetString(2,PLOT_LABEL,"Senkou Span A;Senkou Span B("+string(InpSenkou)+")");
  }
//+------------------------------------------------------------------+
//| Ichimoku Kinko Hyo                                               |
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
   int start;
//---
   if(prev_calculated==0)
      start=0;
   else
      start=prev_calculated-1;
//--- main loop
   for(int i=start; i<rates_total && !IsStopped(); i++)
     {
      ExtChikouBuffer[i]=close[i];
      //--- tenkan sen
      double price_max=Highest(high,InpTenkan,i);
      double price_min=Lowest(low,InpTenkan,i);
      ExtTenkanBuffer[i]=(price_max+price_min)/2.0;
      //--- kijun sen
      price_max=Highest(high,InpKijun,i);
      price_min=Lowest(low,InpKijun,i);
      ExtKijunBuffer[i]=(price_max+price_min)/2.0;
      //--- senkou span a
      ExtSpanABuffer[i]=(ExtTenkanBuffer[i]+ExtKijunBuffer[i])/2.0;
      //--- senkou span b
      price_max=Highest(high,InpSenkou,i);
      price_min=Lowest(low,InpSenkou,i);
      ExtSpanBBuffer[i]=(price_max+price_min)/2.0;
     }
//---
   return(rates_total);
  }
//+------------------------------------------------------------------+
//| get price_max value for range                                      |
//+------------------------------------------------------------------+
double Highest(const double& array[],const int range,int from_index)
  {
   double res=0;
//---
   res=array[from_index];
   for(int i=from_index; i>from_index-range && i>=0; i--)
      if(res<array[i])
         res=array[i];
//---
   return(res);
  }
//+------------------------------------------------------------------+
//| get price_min value for range                                       |
//+------------------------------------------------------------------+
double Lowest(const double& array[],const int range,int from_index)
  {
   double res=0;
//---
   res=array[from_index];
   for(int i=from_index; i>from_index-range && i>=0; i--)
      if(res>array[i])
         res=array[i];
//---
   return(res);
  }
//+------------------------------------------------------------------+
