
#property copyright "2021 ejtrader"
#property link      "https://github.com/ejtraderLabs"
#property version   "1.00"

#include <StringToEnumInt.mqh>
#include <Zmq/Zmq.mqh>
#include <Json.mqh>

// Set ports and host for ZeroMQ
string HOST="localhost";
int CHART_SUB_PORT=15562;

// ZeroMQ Cnnections
Context context("EJTRADERMT");
Socket chartSubscriptionSocket(context,ZMQ_SUB);

//--- input parameters
#property indicator_buffers 31
#property indicator_plots   30

input string            IndicatorId="";
input string            ShortName="ejtraderMTIndicator";

//--- indicator settings
double                  B0[], B1[], B2[], B3[], B4[], B5[], B6[], B7[], B8[], B9[], B10[];
double                  B11[], B12[], B13[], B14[], B15[], B16[], B17[], B18[], B19[], B20[];
double                  B21[], B22[], B23[], B24[], B25[], B26[], B27[], B28[], B29[];
bool                    debug = false;
int                     activeBufferCount = 0;
long                    mtChartId = 0;
bool                    setFormingCandleBlank = true;

//+------------------------------------------------------------------+
//| Custom indicator initialization function                         |
//+------------------------------------------------------------------+
int OnInit()

  {
// TODO subscribe only to own IndicatorId topic
// Subscribe to all topics
   chartSubscriptionSocket.setSubscribe("");
   chartSubscriptionSocket.setLinger(1000);
// https://www.randonomicon.com/zmq/2018/09/19/zmq-bind-vs-connect.html
// Number of messages to buffer in RAM.
// https://dzone.com/articles/zeromq-flow-control-and-other
   chartSubscriptionSocket.setReceiveHighWaterMark(1000);
   bool result = chartSubscriptionSocket.connect(StringFormat("tcp://%s:%d", HOST, CHART_SUB_PORT));
   if(result == false)
     {
      Print("Failed to subscrbe on port ", CHART_SUB_PORT);
     }
   else
     {
      Print("Accepting Chart Indicator data on port ", CHART_SUB_PORT);
     }

//--- indicator buffers mapping;
   ArraySetAsSeries(B0,true);
   ArraySetAsSeries(B1,true);
   ArraySetAsSeries(B2,true);
   ArraySetAsSeries(B3,true);
   ArraySetAsSeries(B4,true);
   ArraySetAsSeries(B5,true);
   ArraySetAsSeries(B6,true);
   ArraySetAsSeries(B7,true);
   ArraySetAsSeries(B8,true);
   ArraySetAsSeries(B9,true);
   ArraySetAsSeries(B10,true);
   ArraySetAsSeries(B11,true);
   ArraySetAsSeries(B12,true);
   ArraySetAsSeries(B13,true);
   ArraySetAsSeries(B14,true);
   ArraySetAsSeries(B15,true);
   ArraySetAsSeries(B16,true);
   ArraySetAsSeries(B17,true);
   ArraySetAsSeries(B18,true);
   ArraySetAsSeries(B19,true);
   ArraySetAsSeries(B20,true);
   ArraySetAsSeries(B21,true);
   ArraySetAsSeries(B22,true);
   ArraySetAsSeries(B23,true);
   ArraySetAsSeries(B24,true);
   ArraySetAsSeries(B25,true);
   ArraySetAsSeries(B26,true);
   ArraySetAsSeries(B27,true);
   ArraySetAsSeries(B28,true);
   ArraySetAsSeries(B29,true);

   SetIndexBuffer(0,B0,INDICATOR_CALCULATIONS);
   SetIndexBuffer(1,B1,INDICATOR_CALCULATIONS);
   SetIndexBuffer(2,B2,INDICATOR_CALCULATIONS);
   SetIndexBuffer(3,B3,INDICATOR_CALCULATIONS);
   SetIndexBuffer(4,B4,INDICATOR_CALCULATIONS);
   SetIndexBuffer(5,B5,INDICATOR_CALCULATIONS);
   SetIndexBuffer(6,B6,INDICATOR_CALCULATIONS);
   SetIndexBuffer(7,B7,INDICATOR_CALCULATIONS);
   SetIndexBuffer(8,B8,INDICATOR_CALCULATIONS);
   SetIndexBuffer(9,B9,INDICATOR_CALCULATIONS);
   SetIndexBuffer(10,B10,INDICATOR_CALCULATIONS);
   SetIndexBuffer(11,B11,INDICATOR_CALCULATIONS);
   SetIndexBuffer(12,B12,INDICATOR_CALCULATIONS);
   SetIndexBuffer(13,B13,INDICATOR_CALCULATIONS);
   SetIndexBuffer(14,B14,INDICATOR_CALCULATIONS);
   SetIndexBuffer(15,B15,INDICATOR_CALCULATIONS);
   SetIndexBuffer(16,B16,INDICATOR_CALCULATIONS);
   SetIndexBuffer(17,B17,INDICATOR_CALCULATIONS);
   SetIndexBuffer(18,B18,INDICATOR_CALCULATIONS);
   SetIndexBuffer(19,B19,INDICATOR_CALCULATIONS);
   SetIndexBuffer(20,B20,INDICATOR_CALCULATIONS);
   SetIndexBuffer(21,B21,INDICATOR_CALCULATIONS);
   SetIndexBuffer(22,B22,INDICATOR_CALCULATIONS);
   SetIndexBuffer(23,B23,INDICATOR_CALCULATIONS);
   SetIndexBuffer(24,B24,INDICATOR_CALCULATIONS);
   SetIndexBuffer(25,B25,INDICATOR_CALCULATIONS);
   SetIndexBuffer(26,B26,INDICATOR_CALCULATIONS);
   SetIndexBuffer(27,B27,INDICATOR_CALCULATIONS);
   SetIndexBuffer(28,B28,INDICATOR_CALCULATIONS);
   SetIndexBuffer(29,B29,INDICATOR_CALCULATIONS);

//---
   IndicatorSetString(INDICATOR_SHORTNAME,ShortName);

   return(INIT_SUCCEEDED);
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void OnDeinit(const int reason)
  {
// Print("INDI DEINIT ",reason);
  }


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void SetStyle(int bufferIdx, string linelabel, color colorstyle, int linetype, int linestyle, int linewidth)
  {
   PlotIndexSetString(bufferIdx,PLOT_LABEL,linelabel);
   PlotIndexSetInteger(bufferIdx,PLOT_LINE_COLOR,0,colorstyle);
   PlotIndexSetInteger(bufferIdx,PLOT_DRAW_TYPE,linetype);
   PlotIndexSetInteger(bufferIdx,PLOT_LINE_STYLE,linestyle);
   PlotIndexSetInteger(bufferIdx,PLOT_LINE_WIDTH,linewidth);
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

// While a new candle is forming, set the current value to be empty


   if(rates_total>prev_calculated && setFormingCandleBlank)
     {
      B0[0] = EMPTY_VALUE;
      B1[0] = EMPTY_VALUE;
      B2[0] = EMPTY_VALUE;
      B3[0] = EMPTY_VALUE;
      B4[0] = EMPTY_VALUE;
      B5[0] = EMPTY_VALUE;
      B6[0] = EMPTY_VALUE;
      B7[0] = EMPTY_VALUE;
      B8[0] = EMPTY_VALUE;
      B9[0] = EMPTY_VALUE;
      B10[0] = EMPTY_VALUE;
      B11[0] = EMPTY_VALUE;
      B12[0] = EMPTY_VALUE;
      B13[0] = EMPTY_VALUE;
      B14[0] = EMPTY_VALUE;
      B15[0] = EMPTY_VALUE;
      B16[0] = EMPTY_VALUE;
      B17[0] = EMPTY_VALUE;
      B18[0] = EMPTY_VALUE;
      B19[0] = EMPTY_VALUE;
      B20[0] = EMPTY_VALUE;
      B21[0] = EMPTY_VALUE;
      B22[0] = EMPTY_VALUE;
      B23[0] = EMPTY_VALUE;
      B24[0] = EMPTY_VALUE;
      B25[0] = EMPTY_VALUE;
      B26[0] = EMPTY_VALUE;
      B27[0] = EMPTY_VALUE;
      B28[0] = EMPTY_VALUE;
      B29[0] = EMPTY_VALUE;
     }

//--- return value of prev_calculated for next call
   return(rates_total);
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void SubscriptionHandler(ZmqMsg &chartMsg)
  {
   CJAVal message;
// Get data from request
   string msg=chartMsg.getData();
   if(debug)
      Print("Processing:"+msg);
// Deserialize msg to CJAVal array
   if(!message.Deserialize(msg))
     {
      Alert("Deserialization Error");
      ExpertRemove();
     }
   if(message["chartIndicatorId"]==IndicatorId)
     {

      if(message["action"]=="PLOT" && message["actionType"]=="DATA")
        {
         int bufferIdx = message["indicatorBufferId"].ToInt();
         if(bufferIdx == 0)
           {
            WriteToBuffer(message, B0);
            SetIndexBuffer(0,B0,INDICATOR_DATA);
           }
         if(bufferIdx == 1)
           {
            WriteToBuffer(message, B1);
            SetIndexBuffer(1,B1,INDICATOR_DATA);
           }
         if(bufferIdx == 2)
           {
            WriteToBuffer(message, B2);
            SetIndexBuffer(2,B2,INDICATOR_DATA);
           }
         if(bufferIdx == 3)
           {
            WriteToBuffer(message, B3);
            SetIndexBuffer(3,B3,INDICATOR_DATA);
           }
         if(bufferIdx == 4)
           {
            WriteToBuffer(message, B4);
            SetIndexBuffer(4,B4,INDICATOR_DATA);
           }
         if(bufferIdx == 5)
           {
            WriteToBuffer(message, B5);
            SetIndexBuffer(5,B5,INDICATOR_DATA);
           }
         if(bufferIdx == 6)
           {
            WriteToBuffer(message, B6);
            SetIndexBuffer(6,B6,INDICATOR_DATA);
           }
         if(bufferIdx == 7)
           {
            WriteToBuffer(message, B7);
            SetIndexBuffer(7,B7,INDICATOR_DATA);
           }
         if(bufferIdx == 8)
           {
            WriteToBuffer(message, B8);
            SetIndexBuffer(8,B8,INDICATOR_DATA);
           }
         if(bufferIdx == 9)
           {
            WriteToBuffer(message, B9);
            SetIndexBuffer(9,B9,INDICATOR_DATA);
           }
         if(bufferIdx == 10)
           {
            WriteToBuffer(message, B10);
            SetIndexBuffer(10,B10,INDICATOR_DATA);
           }
         if(bufferIdx == 11)
           {
            WriteToBuffer(message, B11);
            SetIndexBuffer(11,B11,INDICATOR_DATA);
           }
         if(bufferIdx == 12)
           {
            WriteToBuffer(message, B12);
            SetIndexBuffer(12,B12,INDICATOR_DATA);
           }
         if(bufferIdx == 13)
           {
            WriteToBuffer(message, B13);
            SetIndexBuffer(13,B13,INDICATOR_DATA);
           }
         if(bufferIdx == 14)
           {
            WriteToBuffer(message, B14);
            SetIndexBuffer(14,B14,INDICATOR_DATA);
           }
         if(bufferIdx == 15)
           {
            WriteToBuffer(message, B15);
            SetIndexBuffer(15,B15,INDICATOR_DATA);
           }
         if(bufferIdx == 16)
           {
            WriteToBuffer(message, B16);
            SetIndexBuffer(16,B16,INDICATOR_DATA);
           }
         if(bufferIdx == 17)
           {
            WriteToBuffer(message, B17);
            SetIndexBuffer(17,B17,INDICATOR_DATA);
           }
         if(bufferIdx == 18)
           {
            WriteToBuffer(message, B18);
            SetIndexBuffer(18,B18,INDICATOR_DATA);
           }
         if(bufferIdx == 19)
           {
            WriteToBuffer(message, B19);
            SetIndexBuffer(19,B19,INDICATOR_DATA);
           }
         if(bufferIdx == 20)
           {
            WriteToBuffer(message, B20);
            SetIndexBuffer(20,B20,INDICATOR_DATA);
           }
         if(bufferIdx == 21)
           {
            WriteToBuffer(message, B21);
            SetIndexBuffer(21,B21,INDICATOR_DATA);
           }
         if(bufferIdx == 22)
           {
            WriteToBuffer(message, B22);
            SetIndexBuffer(22,B22,INDICATOR_DATA);
           }
         if(bufferIdx == 23)
           {
            WriteToBuffer(message, B23);
            SetIndexBuffer(23,B23,INDICATOR_DATA);
           }
         if(bufferIdx == 24)
           {
            WriteToBuffer(message, B24);
            SetIndexBuffer(24,B24,INDICATOR_DATA);
           }
         if(bufferIdx == 25)
           {
            WriteToBuffer(message, B25);
            SetIndexBuffer(25,B25,INDICATOR_DATA);
           }
         if(bufferIdx == 26)
           {
            WriteToBuffer(message, B26);
            SetIndexBuffer(26,B26,INDICATOR_DATA);
           }
         if(bufferIdx == 27)
           {
            WriteToBuffer(message, B27);
            SetIndexBuffer(27,B27,INDICATOR_DATA);
           }
         if(bufferIdx == 28)
           {
            WriteToBuffer(message, B28);
            SetIndexBuffer(28,B28,INDICATOR_DATA);
           }
         if(bufferIdx == 29)
           {
            WriteToBuffer(message, B29);
            SetIndexBuffer(29,B29,INDICATOR_DATA);
           }
         ChartRedraw(mtChartId);
        }
      else
         if(message["action"]=="PLOT" && message["actionType"]=="ADDBUFFER")
           {
            string linelabel = message["style"]["linelabel"].ToStr();
            string colorstyleStr = message["style"]["color"].ToStr();
            string linetypeStr = message["style"]["linetype"].ToStr();
            string linestyleStr = message["style"]["linestyle"].ToStr();
            int linewidth = message["style"]["linewidth"].ToInt();
            setFormingCandleBlank = message["style"]["blankforming"].ToBool();

            color colorstyle = StringToColor(colorstyleStr);
            int linetype = StringToEnumInt(linetypeStr);
            int linestyle = StringToEnumInt(linestyleStr);

            SetStyle(activeBufferCount, linelabel, colorstyle, linetype, linestyle, linewidth);
            activeBufferCount = activeBufferCount + 1;

            ClearBuffer(activeBufferCount-1);
           }
     }
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void Clear(double &buffer[])
  {
   int bufferSize = ArraySize(buffer);
   for(int i=0; i<bufferSize; i++)
     {
      buffer[i] = EMPTY_VALUE;
     }
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void ClearBuffer(int bufferIdx)
  {
   switch(bufferIdx)
     {
      case  0:
        {
         Clear(B0);
        }
      case  1:
        {
         Clear(B1);
        }
      case  2:
        {
         Clear(B2);
        }
      case  3:
        {
         Clear(B3);
        }
      case  4:
        {
         Clear(B4);
        }
      case  5:
        {
         Clear(B5);
        }
      case  6:
        {
         Clear(B6);
        }
      case  7:
        {
         Clear(B7);
        }
      case  8:
        {
         Clear(B8);
        }
      case  9:
        {
         Clear(B9);
        }
      case  10:
        {
         Clear(B10);
        }
      case  11:
        {
         Clear(B11);
        }
      case  12:
        {
         Clear(B12);
        }
      case  13:
        {
         Clear(B13);
        }
      case  14:
        {
         Clear(B14);
        }
      case  15:
        {
         Clear(B15);
        }
      case  16:
        {
         Clear(B16);
        }
      case  17:
        {
         Clear(B17);
        }
      case  18:
        {
         Clear(B18);
        }
      case  19:
        {
         Clear(B19);
        }
      case  20:
        {
         Clear(B20);
        }
      case  21:
        {
         Clear(B21);
        }
      case  22:
        {
         Clear(B22);
        }
      case  23:
        {
         Clear(B23);
        }
      case  24:
        {
         Clear(B24);
        }
      case  25:
        {
         Clear(B25);
        }
      case  26:
        {
         Clear(B26);
        }
      case  27:
        {
         Clear(B27);
        }
      case  28:
        {
         Clear(B28);
        }
      case  29:
        {
         Clear(B29);
        }
      break;
      default:
        {} break;
     }
  }

//+------------------------------------------------------------------+
//| Update indicator buffer function                                 |
//+------------------------------------------------------------------+
void WriteToBuffer(CJAVal &message, double &buffer[])
  {

   int bufferSize = ArraySize(buffer);
   int messageDataSize = message["data"].Size();

// calculate the buffer offset
   MqlRates r[];
   mtChartId =(datetime)message["mtChartId"].ToInt();
   datetime fromDate=(datetime)message["fromDate"].ToInt();
   datetime toDate=TimeCurrent();
   ENUM_TIMEFRAMES period = ChartPeriod(mtChartId);
   string symbol = ChartSymbol(mtChartId);
   int rateCount;

   rateCount = CopyRates(symbol, period, fromDate, toDate, r);
   int offset = rateCount - 1;

// write to buffer
   for(int i=0; i<messageDataSize; i++)
     {
      // don't add more elements than the automatically sized buffer array can
      if(i+offset<bufferSize)
        {
         double val = message["data"][i].ToDbl();
         if(val >= EMPTY_VALUE)
            val = EMPTY_VALUE;
         buffer[i+offset] = val;
        }
     }
  }


//+------------------------------------------------------------------+
//| Check for new indicator data function                            |
//+------------------------------------------------------------------+
void CheckMessages()
  {
// This is a workaround for Timer(). It is needed, because OnTimer() works if the indicator is manually added to a chart, but not with ChartIndicatorAdd()

   ZmqMsg chartMsg;

// Recieve chart instructions stream from client via live Chart socket.
   chartSubscriptionSocket.recv(chartMsg,true);

// Request recieved
   if(chartMsg.size()>0)
     {
      // Handle subscription SubscriptionHandler()
      SubscriptionHandler(chartMsg);
      ChartRedraw(ChartID());
     }
  }

//+------------------------------------------------------------------+
//| OnTimer() workaround function                                    |
//+------------------------------------------------------------------+
// Gets triggered by the OnTimer() function of the JsonAPI Expert script
void OnChartEvent(const int id,
                  const long &lparam,
                  const double &dparam,
                  const string &sparam)
  {
   if(id==CHARTEVENT_CUSTOM+222)
      CheckMessages();
  }
//+----------------------------------------------------
