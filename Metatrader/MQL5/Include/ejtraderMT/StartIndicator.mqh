#property copyright "ejtrader"
#property link      "https://github.com/ejtraderLabs/MQL5-ejtraderMT"

#define START_INDICATOR true

int INDICATOR_DATA_PORT=15559;

Socket indicatorDataSocket(context,ZMQ_PUSH);

struct Indicator
  {
   long              id; // Internal id
   string            indicatorId; // UUID
   int               indicatorHandle; // Internal id/handle
   int               indicatorParamCount; // Number of parameters to be passed to the indicator
   int               indicatorBufferCount; // Numnber of buffers to be returned bythe indicator
  };

Indicator indicators[];
int indicatorCount = 0;

//+------------------------------------------------------------------+
//| Check if string is a representation of a number                  |
//+------------------------------------------------------------------+
bool IsNumberAsString(string str)
  {
// MQL5 seems to return true if the values are the same, no matter the data type, in this case comparing str and dbl/int.
// (str "2.1" == double 2.1) will return true.
   double dbl = StringToDouble(str);
   int integer = StringToInteger(str);
// Compaing to both int and double to cover both cases
   if(str==dbl || str==integer)
      return true;
   else
      return false;
  }

//+------------------------------------------------------------------+
//| Get index of indicator handler array by indicator id string      |
//+------------------------------------------------------------------+
int GetIndicatorIdxByIndicatorId(string indicatorId)
  {
   for(int i=0; i<indicatorCount; i++)
     {
      if(indicators[i].indicatorId == indicatorId)
        {
         return i;
        }
     }
   return -1;
  }

//+------------------------------------------------------------------+
//| Start new indicator or request indicator data                    |
//+------------------------------------------------------------------+
void IndicatorControl(CJAVal &dataObject)
  {

   string actionType=dataObject["actionType"].ToStr();

   if(actionType=="REQUEST")
     {
      GetIndicatorResult(dataObject);
     }
   else
      if(actionType=="ATTACH")
        {
         StartIndicator(dataObject);
        }
  }

//+------------------------------------------------------------------+
//| Start new indicator instance                                     |
//+------------------------------------------------------------------+
void StartIndicator(CJAVal &dataObject)
  {

// TODO map Indicators Constants https://www.mql5.com/en/docs/constants/indicatorconstants

   string symbol=dataObject["symbol"].ToStr();
   string chartTF=dataObject["chartTF"].ToStr();
   string id=dataObject["id"].ToStr();
   string indicatorName=dataObject["name"].ToStr();

   indicatorCount++;
   ArrayResize(indicators,indicatorCount);

   int idx = indicatorCount-1;

   indicators[idx].indicatorId = id;
   indicators[idx].indicatorBufferCount = dataObject["linecount"].ToInt();

   double params[];
   indicators[idx].indicatorParamCount = dataObject["params"].Size();
   for(int i=0; i<indicators[idx].indicatorParamCount; i++)
     {
      // TODO test it. Is it ok to pass EnumInts as Doubles for params?
      ArrayResize(params, i+1);
      string paramStr = dataObject["params"][i].ToStr();
      if(IsNumberAsString(paramStr))
         params[i] = StringToDouble(paramStr);
      else
        {
         params[i] = StringToEnumInt(paramStr);
         mControl.mResetLastError(); // TODO find where the Error 4003 is craeted in StringToEnumInt
        }
     }

   ENUM_TIMEFRAMES period = GetTimeframe(chartTF);

// Case construct for passing variable parameter count to the iCustom function is used, because MQL5 does not seem to support expanding an array to a function parameter list
   switch(indicators[idx].indicatorParamCount)
     {
      case 0:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName);
         break;
      case 1:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0]);
         break;
      case 2:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1]);
         break;
      case 3:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1], params[2]);
         break;
      case 4:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1], params[2], params[3]);
         break;
      case 5:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1], params[2], params[3], params[4]);
         break;
      case 6:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1], params[2], params[3], params[4], params[5]);
         break;
      case 7:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
         break;
      case 8:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7]);
         break;
      case 9:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8]);
         break;
      case 10:
         indicators[idx].indicatorHandle = iCustom(symbol,period,indicatorName, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9]);
         break;
      default:
         // TODO error handling
         break;
     }

   CJAVal message;

   if(mControl.mGetLastError())
     {
      int lastError = mControl.mGetLastError();
      string desc = mControl.mGetDesc();
      mControl.Check();

      message["error"]=(bool) true;
      message["lastError"]=(string) lastError;
      message["description"]=desc;
      message["function"]=(string) __FUNCTION__;
      string t=message.Serialize();
      if(debug)
         Print(t);
      InformClientSocket(indicatorDataSocket,t);
     }
   else
     {
      message["error"]=(bool) false;
      message["id"] = (string) id;

      string t=message.Serialize();
      if(debug)
         Print(t);
      InformClientSocket(indicatorDataSocket,t);
     }

  }

//+------------------------------------------------------------------+
//| Get indicator results                                            |
//+------------------------------------------------------------------+
void GetIndicatorResult(CJAVal &dataObject)
  {

   datetime fromDate=dataObject["fromDate"].ToInt();
   string id=dataObject["id"].ToStr();
   string indicatorName=dataObject["indicatorName"].ToStr();

   int idx = GetIndicatorIdxByIndicatorId(id);

   double values[2];

   CJAVal results;
// Cycle through all avaliable buffer positions
   for(int i=0; i<indicators[idx].indicatorBufferCount; i++)
     {
      values[0] = 0.0;
      values[1] = 0.0;
      results[i] = 0.0;
      if(idx >= 0)
        {
         if(CopyBuffer(indicators[idx].indicatorHandle, i, fromDate, 1, values) < 0)
           {
            if(mControl.mGetLastError())
              {
               CJAVal message;
               int lastError = mControl.mGetLastError();
               string desc = mControl.mGetDesc();
               mControl.Check();

               message["error"]=(bool) true;
               message["lastError"]=(string) lastError;
               message["description"]=desc;
               message["function"]=(string) __FUNCTION__;
               string t=message.Serialize();
               if(debug)
                  Print(t);
               InformClientSocket(indicatorDataSocket,t);
              }
           }
         results[i] = DoubleToString(values[0]);
        }
     }

   CJAVal message;
   message["error"]=(bool) false;
   message["id"] = (string) id;
   message["data"].Set(results);

   string t=message.Serialize();
   if(debug)
      Print(t);
   InformClientSocket(indicatorDataSocket,t);

  }
//+------------------------------------------------------------------+
