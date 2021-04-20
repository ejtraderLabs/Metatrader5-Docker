//+------------------------------------------------------------------+
//|                                        Custom Moving Average.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "2009-2020, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"

//--- indicator settings
#property indicator_chart_window
#property indicator_buffers 1
#property indicator_plots   1
#property indicator_type1   DRAW_LINE
#property indicator_color1  Red
//--- input parameters
input int            InpMAPeriod=13;         // Period
input int            InpMAShift=0;           // Shift
input ENUM_MA_METHOD InpMAMethod=MODE_SMMA;  // Method
//--- indicator buffer
double ExtLineBuffer[];
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtLineBuffer,INDICATOR_DATA);
//--- set accuracy
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits+1);
//--- set first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,InpMAPeriod);
//--- line shifts when drawing
   PlotIndexSetInteger(0,PLOT_SHIFT,InpMAShift);
//--- name for DataWindow
   string short_name;
   switch(InpMAMethod)
     {
      case MODE_EMA :
         short_name="EMA";
         break;
      case MODE_LWMA :
         short_name="LWMA";
         break;
      case MODE_SMA :
         short_name="SMA";
         break;
      case MODE_SMMA :
         short_name="SMMA";
         break;
      default :
         short_name="unknown ma";
     }
   IndicatorSetString(INDICATOR_SHORTNAME,short_name+"("+string(InpMAPeriod)+")");
//--- set drawing line empty value
   PlotIndexSetDouble(0,PLOT_EMPTY_VALUE,0.0);
  }
//+------------------------------------------------------------------+
//|  Moving Average                                                  |
//+------------------------------------------------------------------+
int OnCalculate(const int rates_total,
                const int prev_calculated,
                const int begin,
                const double &price[])
  {
   if(rates_total<InpMAPeriod-1+begin)
      return(0);
//--- first calculation or number of bars was changed
   if(prev_calculated==0)
     {
      ArrayInitialize(ExtLineBuffer,0);
      PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,InpMAPeriod-1+begin);
     }
//--- calculation
   switch(InpMAMethod)
     {
      case MODE_EMA:
         CalculateEMA(rates_total,prev_calculated,begin,price);
         break;
      case MODE_LWMA:
         CalculateLWMA(rates_total,prev_calculated,begin,price);
         break;
      case MODE_SMMA:
         CalculateSmoothedMA(rates_total,prev_calculated,begin,price);
         break;
      case MODE_SMA:
         CalculateSimpleMA(rates_total,prev_calculated,begin,price);
         break;
     }
//--- return value of prev_calculated for next call
   return(rates_total);
  }
//+------------------------------------------------------------------+
//|   simple moving average                                          |
//+------------------------------------------------------------------+
void CalculateSimpleMA(int rates_total,int prev_calculated,int begin,const double &price[])
  {
   int i,start;
//--- first calculation or number of bars was changed
   if(prev_calculated==0)
     {
      start=InpMAPeriod+begin;
      //--- set empty value for first start bars
      for(i=0; i<start-1; i++)
         ExtLineBuffer[i]=0.0;
      //--- calculate first visible value
      double first_value=0;
      for(i=begin; i<start; i++)
         first_value+=price[i];
      first_value/=InpMAPeriod;
      ExtLineBuffer[start-1]=first_value;
     }
   else
      start=prev_calculated-1;
//--- main loop
   for(i=start; i<rates_total && !IsStopped(); i++)
      ExtLineBuffer[i]=ExtLineBuffer[i-1]+(price[i]-price[i-InpMAPeriod])/InpMAPeriod;
  }
//+------------------------------------------------------------------+
//|  exponential moving average                                      |
//+------------------------------------------------------------------+
void CalculateEMA(int rates_total,int prev_calculated,int begin,const double &price[])
  {
   int    i,start;
   double SmoothFactor=2.0/(1.0+InpMAPeriod);
//--- first calculation or number of bars was changed
   if(prev_calculated==0)
     {
      start=InpMAPeriod+begin;
      ExtLineBuffer[begin]=price[begin];
      for(i=begin+1; i<start; i++)
         ExtLineBuffer[i]=price[i]*SmoothFactor+ExtLineBuffer[i-1]*(1.0-SmoothFactor);
     }
   else
      start=prev_calculated-1;
//--- main loop
   for(i=start; i<rates_total && !IsStopped(); i++)
      ExtLineBuffer[i]=price[i]*SmoothFactor+ExtLineBuffer[i-1]*(1.0-SmoothFactor);
  }
//+------------------------------------------------------------------+
//|  linear weighted moving average                                  |
//+------------------------------------------------------------------+
void CalculateLWMA(int rates_total,int prev_calculated,int begin,const double &price[])
  {
   int    weight=0;
   int    i,l,start;
   double sum=0.0,lsum=0.0;
//--- first calculation or number of bars was changed
   if(prev_calculated<=InpMAPeriod+begin+2)
     {
      start=InpMAPeriod+begin;
      //--- set empty value for first start bars
      for(i=0; i<start; i++)
         ExtLineBuffer[i]=0.0;
     }
   else
      start=prev_calculated-1;

   for(i=start-InpMAPeriod,l=1; i<start; i++,l++)
     {
      sum   +=price[i]*l;
      lsum  +=price[i];
      weight+=l;
     }
   ExtLineBuffer[start-1]=sum/weight;
//--- main loop
   for(i=start; i<rates_total && !IsStopped(); i++)
     {
      sum             =sum-lsum+price[i]*InpMAPeriod;
      lsum            =lsum-price[i-InpMAPeriod]+price[i];
      ExtLineBuffer[i]=sum/weight;
     }
  }
//+------------------------------------------------------------------+
//|  smoothed moving average                                         |
//+------------------------------------------------------------------+
void CalculateSmoothedMA(int rates_total,int prev_calculated,int begin,const double &price[])
  {
   int i,start;
//--- first calculation or number of bars was changed
   if(prev_calculated==0)
     {
      start=InpMAPeriod+begin;
      //--- set empty value for first start bars
      for(i=0; i<start-1; i++)
         ExtLineBuffer[i]=0.0;
      //--- calculate first visible value
      double first_value=0;
      for(i=begin; i<start; i++)
         first_value+=price[i];
      first_value/=InpMAPeriod;
      ExtLineBuffer[start-1]=first_value;
     }
   else
      start=prev_calculated-1;
//--- main loop
   for(i=start; i<rates_total && !IsStopped(); i++)
      ExtLineBuffer[i]=(ExtLineBuffer[i-1]*(InpMAPeriod-1)+price[i])/InpMAPeriod;
  }
//+------------------------------------------------------------------+
