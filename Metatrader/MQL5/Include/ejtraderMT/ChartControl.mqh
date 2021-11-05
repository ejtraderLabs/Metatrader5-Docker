

#property copyright "ejtrader"
#property link      "https://github.com/ejtraderLabs/MQL5-ejtraderMT"

#define CHART_CONTROL true

int CHART_DATA_PORT=15560;
int CHART_INDICATOR_DATA_PORT=15562;

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
Socket chartDataSocket(context,ZMQ_PULL);
Socket chartIndicatorDataSocket(context,ZMQ_PUB);

// Variables for controlling chart
struct ChartWindow
  {
   long              id; // Internal id
   string            chartId; // UUID
  };
  
ChartWindow chartWindows[];
int chartWindowCount = 0;

struct ChartWindowIndicator
  {
   long              id; // Internal id
   string            indicatorId; // UUID
   int               indicatorHandle; // Internal id/handle
  };

ChartWindowIndicator chartWindowIndicators[];
int chartWindowIndicatorCount = 0;

//+----------------------------------------------------------------------+
//| Function return index of chart window array by chart window id string|
//+----------------------------------------------------------------------+
int GetChartWindowIdxByChartWindowId(string chartWindowId)
  {
   for(int i=0; i<chartWindowCount; i++)
     {
      if(chartWindows[i].chartId == chartWindowId)
        {
         return i;
        }
     }
   return -1;
  }

//+------------------------------------------------------------------+
//| Open new chart or add indicator to chart                         |
//+------------------------------------------------------------------+
void ChartControl(CJAVal &dataObject)
  {

   string actionType=dataObject["actionType"].ToStr();

   if(actionType=="ADDINDICATOR")
     {
      AddChartIndicator(dataObject);
     }
   else
      if(actionType=="OPEN")
        {
         OpenChart(dataObject);
        }
  }

//+------------------------------------------------------------------+
//| Open new chart                                                   |
//+------------------------------------------------------------------+
void OpenChart(CJAVal &dataObject)
  {

   string chartId=dataObject["chartId"].ToStr();
   string symbol=dataObject["symbol"].ToStr();
   string chartTF=dataObject["chartTF"].ToStr();

   chartWindowCount++;
   ArrayResize(chartWindows,chartWindowCount);

   int idx = chartWindowCount-1;

   chartWindows[idx].chartId = chartId;

   ENUM_TIMEFRAMES period = GetTimeframe(chartTF);
   chartWindows[idx].id = ChartOpen(symbol, period);
   ChartSetInteger(chartWindows[idx].id, CHART_AUTOSCROLL, false);

   CJAVal message;
   message["error"]=(bool) false;
   message["chartId"] = (string) chartId;
   message["mtChartId"] = (string) chartWindows[idx].id;

   string t=message.Serialize();
   if(debug)
      Print(t);
   InformClientSocket(dataSocket,t);
  }

//+------------------------------------------------------------------+
//| Add JsonAPIIndicator indicator to chart                          |
//+------------------------------------------------------------------+
void AddChartIndicator(CJAVal &dataObject)
  {

   string chartIdStr=dataObject["chartId"].ToStr();
   string chartIndicatorId=dataObject["chartIndicatorId"].ToStr();
   int chartIndicatorSubWindow=dataObject["chartIndicatorSubWindow"].ToInt();
   string shortName = dataObject["shortName"].ToStr();

   int chartIdx = GetChartWindowIdxByChartWindowId(chartIdStr);
   long chartId = chartWindows[chartIdx].id;

   double chartIndicatorHandle = iCustom(ChartSymbol(chartId),ChartPeriod(chartId),"ejtraderMTIndicator",chartIndicatorId,shortName); //linelabel,colorstyle,linetype,linestyle,linewidth);

   if(ChartIndicatorAdd(chartId, chartIndicatorSubWindow, chartIndicatorHandle))
     {
      chartWindowIndicatorCount++;
      ArrayResize(chartWindowIndicators,chartWindowIndicatorCount);
      int indicatorIdx = chartWindowIndicatorCount-1;
      chartWindowIndicators[indicatorIdx].indicatorId = chartIndicatorId;
      chartWindowIndicators[indicatorIdx].indicatorHandle = chartIndicatorHandle;
     }
   if(!CheckError(__FUNCTION__))
     {
      CJAVal message;
      message["error"]=(bool) false;
      message["chartId"] = (string) chartIdStr;

      string t=message.Serialize();
      if(debug)
         Print(t);
      InformClientSocket(dataSocket,t);
     }
  }
//+------------------------------------------------------------------+
