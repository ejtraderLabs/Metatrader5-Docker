//+------------------------------------------------------------------+
//|                                                    Envelopes.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "2009-2020, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
//--- indicator settings
#property indicator_chart_window
#property indicator_buffers 3
#property indicator_plots   2
#property indicator_type1   DRAW_LINE
#property indicator_type2   DRAW_LINE
#property indicator_color1  Blue
#property indicator_color2  Red
#property indicator_label1  "Upper band"
#property indicator_label2  "Lower band"
//--- input parameters
input int                InpMAPeriod=14;              // Period
input int                InpMAShift=0;                // Shift
input ENUM_MA_METHOD     InpMAMethod=MODE_SMA;        // Method
input ENUM_APPLIED_PRICE InpAppliedPrice=PRICE_CLOSE; // Applied price
input double             InpDeviation=0.1;            // Deviation
//--- indicator buffers
double                   ExtUpBuffer[];
double                   ExtDownBuffer[];
double                   ExtMABuffer[];

int                      ExtMAHandle;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- indicator buffers mapping
   SetIndexBuffer(0,ExtUpBuffer,INDICATOR_DATA);
   SetIndexBuffer(1,ExtDownBuffer,INDICATOR_DATA);
   SetIndexBuffer(2,ExtMABuffer,INDICATOR_CALCULATIONS);
//---
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits+1);
//--- sets first bar from what index will be drawn
   PlotIndexSetInteger(0,PLOT_DRAW_BEGIN,InpMAPeriod-1);
//--- name for DataWindow
   string short_name=StringFormat("Env(%d)",InpMAPeriod);
   IndicatorSetString(INDICATOR_SHORTNAME,short_name);
   PlotIndexSetString(0,PLOT_LABEL,short_name+" Upper");
   PlotIndexSetString(1,PLOT_LABEL,short_name+" Lower");
//--- line shifts when drawing
   PlotIndexSetInteger(0,PLOT_SHIFT,InpMAShift);
   PlotIndexSetInteger(1,PLOT_SHIFT,InpMAShift);
//---
   ExtMAHandle=iMA(NULL,0,InpMAPeriod,0,InpMAMethod,InpAppliedPrice);
  }
//+------------------------------------------------------------------+
//| Envelopes                                                        |
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
   if(rates_total<InpMAPeriod)
      return(0);

   int calculated=BarsCalculated(ExtMAHandle);
   if(calculated<rates_total)
     {
      Print("Not all data of ExtMAHandle is calculated (",calculated," bars). Error ",GetLastError());
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
//--- get ma buffer
   if(IsStopped()) // checking for stop flag
      return(0);
   if(CopyBuffer(ExtMAHandle,0,0,to_copy,ExtMABuffer)<=0)
     {
      Print("Getting MA data is failed! Error ",GetLastError());
      return(0);
     }
//--- preliminary calculations
   int start=prev_calculated-1;
   if(start<InpMAPeriod)
      start=InpMAPeriod;
//--- the main loop of calculations
   for(int i=start; i<rates_total && !IsStopped(); i++)
     {
      ExtUpBuffer[i]=(1+InpDeviation/100.0)*ExtMABuffer[i];
      ExtDownBuffer[i]=(1-InpDeviation/100.0)*ExtMABuffer[i];
     }
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
