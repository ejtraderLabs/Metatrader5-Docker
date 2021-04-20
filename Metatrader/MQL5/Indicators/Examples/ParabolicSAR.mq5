//+------------------------------------------------------------------+
//|                                                 ParabolicSAR.mq5 |
//|                   Copyright 2009-2020, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "2009-2020, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
//--- indicator settings
#property indicator_chart_window
#property indicator_buffers 3
#property indicator_plots   1
#property indicator_type1   DRAW_ARROW
#property indicator_color1  DodgerBlue
//--- input parametrs
input double InpSARStep=0.02;    // Step
input double InpSARMaximum=0.2;  // Maximum
//--- indicator buffers
double ExtSARBuffer[];
double ExtEPBuffer[];
double ExtAFBuffer[];

int    ExtLastRevPos;
bool   ExtDirectionLong;
double ExtSarStep;
double ExtSarMaximum;
//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
void OnInit()
  {
//--- checking input data
   if(InpSARStep<0.0)
     {
      ExtSarStep=0.02;
      PrintFormat("Input parametr InpSARStep has incorrect value. Indicator will use value %d for calculations.",
                  ExtSarStep);
     }
   else
      ExtSarStep=InpSARStep;
   if(InpSARMaximum<0.0)
     {
      ExtSarMaximum=0.2;
      PrintFormat("Input parametr InpSARMaximum has incorrect value. Indicator will use value %d for calculations.",
                  ExtSarMaximum);
     }
   else
      ExtSarMaximum=InpSARMaximum;
//--- indicator buffers
   SetIndexBuffer(0,ExtSARBuffer);
   SetIndexBuffer(1,ExtEPBuffer,INDICATOR_CALCULATIONS);
   SetIndexBuffer(2,ExtAFBuffer,INDICATOR_CALCULATIONS);
//--- set arrow symbol
   PlotIndexSetInteger(0,PLOT_ARROW,159);
//--- set indicator digits
   IndicatorSetInteger(INDICATOR_DIGITS,_Digits);
//--- set label name
   string short_name=StringFormat("SAR(%.2f,%.2f)",ExtSarStep,ExtSarMaximum);
   PlotIndexSetString(0,PLOT_LABEL,short_name);

   ExtLastRevPos=0;
   ExtDirectionLong=false;
  }
//+------------------------------------------------------------------+
//| Custom indicator iteration function                              |
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
   if(rates_total<3)
      return(0);
//--- detect current position
   int pos=prev_calculated-1;
//--- correct position
   if(pos<1)
     {
      //--- first pass, set as SHORT
      pos=1;
      ExtAFBuffer[0]=ExtSarStep;
      ExtAFBuffer[1]=ExtSarStep;
      ExtSARBuffer[0]=high[0];
      ExtLastRevPos=0;
      ExtDirectionLong=false;
      ExtSARBuffer[1]=GetHigh(pos,ExtLastRevPos,high);
      ExtEPBuffer[0]=low[pos];
      ExtEPBuffer[1]=low[pos];
     }
//---main cycle
   for(int i=pos; i<rates_total-1 && !IsStopped(); i++)
     {
      //--- check for reverse
      if(ExtDirectionLong)
        {
         if(ExtSARBuffer[i]>low[i])
           {
            //--- switch to SHORT
            ExtDirectionLong=false;
            ExtSARBuffer[i]=GetHigh(i,ExtLastRevPos,high);
            ExtEPBuffer[i]=low[i];
            ExtLastRevPos=i;
            ExtAFBuffer[i]=ExtSarStep;
           }
        }
      else
        {
         if(ExtSARBuffer[i]<high[i])
           {
            //--- switch to LONG
            ExtDirectionLong=true;
            ExtSARBuffer[i]=GetLow(i,ExtLastRevPos,low);
            ExtEPBuffer[i]=high[i];
            ExtLastRevPos=i;
            ExtAFBuffer[i]=ExtSarStep;
           }
        }
      //--- continue calculations
      if(ExtDirectionLong)
        {
         //--- check for new High
         if(high[i]>ExtEPBuffer[i-1] && i!=ExtLastRevPos)
           {
            ExtEPBuffer[i]=high[i];
            ExtAFBuffer[i]=ExtAFBuffer[i-1]+ExtSarStep;
            if(ExtAFBuffer[i]>ExtSarMaximum)
               ExtAFBuffer[i]=ExtSarMaximum;
           }
         else
           {
            //--- when we haven't reversed
            if(i!=ExtLastRevPos)
              {
               ExtAFBuffer[i]=ExtAFBuffer[i-1];
               ExtEPBuffer[i]=ExtEPBuffer[i-1];
              }
           }
         //--- calculate SAR for tomorrow
         ExtSARBuffer[i+1]=ExtSARBuffer[i]+ExtAFBuffer[i]*(ExtEPBuffer[i]-ExtSARBuffer[i]);
         //--- check for SAR
         if(ExtSARBuffer[i+1]>low[i] || ExtSARBuffer[i+1]>low[i-1])
            ExtSARBuffer[i+1]=MathMin(low[i],low[i-1]);
        }
      else
        {
         //--- check for new Low
         if(low[i]<ExtEPBuffer[i-1] && i!=ExtLastRevPos)
           {
            ExtEPBuffer[i]=low[i];
            ExtAFBuffer[i]=ExtAFBuffer[i-1]+ExtSarStep;
            if(ExtAFBuffer[i]>ExtSarMaximum)
               ExtAFBuffer[i]=ExtSarMaximum;
           }
         else
           {
            //--- when we haven't reversed
            if(i!=ExtLastRevPos)
              {
               ExtAFBuffer[i]=ExtAFBuffer[i-1];
               ExtEPBuffer[i]=ExtEPBuffer[i-1];
              }
           }
         //--- calculate SAR for tomorrow
         ExtSARBuffer[i+1]=ExtSARBuffer[i]+ExtAFBuffer[i]*(ExtEPBuffer[i]-ExtSARBuffer[i]);
         //--- check for SAR
         if(ExtSARBuffer[i+1]<high[i] || ExtSARBuffer[i+1]<high[i-1])
            ExtSARBuffer[i+1]=MathMax(high[i],high[i-1]);
        }
     }
//--- OnCalculate done. Return new prev_calculated.
   return(rates_total);
  }
//+------------------------------------------------------------------+
//| Find highest price from start to current position                |
//+------------------------------------------------------------------+
double GetHigh(int curr_pos,int start,const double& high[])
  {
   double result=high[start];
//---
   for(int i=start+1; i<=curr_pos; i++)
      if(result<high[i])
         result=high[i];
//---
   return(result);
  }
//+------------------------------------------------------------------+
//| Find lowest price from start to current position                 |
//+------------------------------------------------------------------+
double GetLow(int curr_pos,int start,const double& low[])
  {
   double result=low[start];
//---
   for(int i=start+1; i<=curr_pos; i++)
      if(result>low[i])
         result=low[i];
//---
   return(result);
  }
//+------------------------------------------------------------------+
