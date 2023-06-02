//+------------------------------------------------------------------+
//|                                                      ProjectName |
//|                                      Copyright 2020, CompanyName |
//|                                       http://www.companyname.net |
//+------------------------------------------------------------------+

#property copyright "Copyright 2022, ejtrader."
#property link "https://github.com/ejtraderLabs"
#property version "3.16"
#property description "ejtraderMT"
#property description "See github link for documentation"

#include <Trade/AccountInfo.mqh>
#include <Trade/DealInfo.mqh>
#include <Trade/Trade.mqh>
#include <Zmq/Zmq.mqh>
#include <Json.mqh>
#include <StringToEnumInt.mqh>
#include <ControlErrors.mqh>

// Set ports and host for ZeroMQ
string HOST = "192.168.1.154";
int SYS_PORT = 15557;

// ZeroMQ Connections
Context context("EJTRADERMT");
Socket sysSocket(context, ZMQ_REP);

// Load ejtraderMTincludes
// Required:
#include <ejtraderMT/HistoryInfo.mqh>
#include <ejtraderMT/Broker.mqh>
#include <ejtraderMT/Calendar.mqh>

// Global variables \\
bool debug = true;
bool liveStream = false;
bool connectedFlag = true;
int deInitReason = -1;

// Variables for handling price data stream
struct SymbolSubscription
{
   string symbol;
   string chartTf;
   datetime lastBar;
};
SymbolSubscription symbolSubscriptions[];
int symbolSubscriptionCount = 0;

// Error handling
ControlErrors mControl;
datetime tm;
//+------------------------------------------------------------------+
//| Bind ZMQ sockets to ports                                        |
//+------------------------------------------------------------------+
bool BindSockets()
{
   sysSocket.setLinger(1000);

   bool result = false;
   result = sysSocket.connect(StringFormat("tcp://%s:%d", HOST, SYS_PORT));
   if (result == false)
   {
      return result;
   }
   else
   {
      Print("Bound 'System' socket on port ", SYS_PORT);
   }

   return result;
}

//+------------------------------------------------------------------+
//| Expert initialization function                                   |
//+------------------------------------------------------------------+
int OnInit()
{

   // Setting up error reporting
   mControl.SetAlert(true);
   mControl.SetSound(false);
   mControl.SetWriteFlag(false);

   /* Bindinig ZMQ ports on init */
   // Skip reloading of the EA script when the reason to reload is a chart timeframe change
   if (deInitReason != REASON_CHARTCHANGE)
   {

      EventSetMillisecondTimer(1);

      int bindSocketsDelay = 65; // Seconds to wait if binding of sockets fails.
      int bindAttemtps = 2;      // Number of binding attemtps

      Print("Binding sockets...");

      for (int i = 0; i < bindAttemtps; i++)
      {
         if (BindSockets())
            return (INIT_SUCCEEDED);
         else
         {
            Print("Binding sockets failed. Waiting ", bindSocketsDelay, " seconds to try again...");
            Sleep(bindSocketsDelay * 1000);
         }
      }

      Print("Binding of sockets failed permanently.");
      return (INIT_FAILED);
   }

   return (INIT_SUCCEEDED);
}

//+------------------------------------------------------------------+
//| Expert deinitialization function                                 |
//+------------------------------------------------------------------+
void OnDeinit(const int reason)
{

   deInitReason = reason;

   // Skip reloading of the EA script when the reason to reload is a chart timeframe change
   if (reason != REASON_CHARTCHANGE)
   {
      Print(__FUNCTION__, " Deinitialization reason: ", getUninitReasonText(reason));

      Print("Unbinding 'System' socket on port ", SYS_PORT, "..");
      sysSocket.disconnect(StringFormat("tcp://%s:%d", HOST, SYS_PORT));

      // Shutdown ZeroMQ Context
      context.shutdown();
      context.destroy(0);

      EventKillTimer();
   }
}

//+------------------------------------------------------------------+
//| Check if subscribed to symbol and timeframe combination          |
//+------------------------------------------------------------------+
bool HasChartSymbol(string symbol, string chartTF)
{
   for (int i = 0; i < ArraySize(symbolSubscriptions); i++)
   {
      if (symbolSubscriptions[i].symbol == symbol && symbolSubscriptions[i].chartTf == chartTF)
      {
         return true;
      }
   }
   return false;
}

//+------------------------------------------------------------------+
//| Stream live price data                                           |
//+------------------------------------------------------------------+
void StreamPriceData()
{
   // If liveStream == true, push last candle to sysSocket.

   if (liveStream)
   {
      CJAVal last;
      if (TerminalInfoInteger(TERMINAL_CONNECTED))
      {
         connectedFlag = true;
         for (int i = 0; i < symbolSubscriptionCount; i++)
         {
            string symbol = symbolSubscriptions[i].symbol;
            string chartTF = symbolSubscriptions[i].chartTf;
            datetime lastBar = symbolSubscriptions[i].lastBar;
            CJAVal Data;
            ENUM_TIMEFRAMES period = GetTimeframe(chartTF);

            datetime thisBar = 0;
            float price;
            MqlTick tick;
            MqlRates rates[1];
            int spread[1];

            if (chartTF == "TICK")
            {
               if (SymbolInfoTick(symbol, tick) != true)
               { /*mControl.Check();*/
               }
               thisBar = (datetime)tick.time_msc;
            }
            else
            {
               if (CopyRates(symbol, period, 1, 1, rates) != 1)
               { /*mControl.Check();*/
               }
               if (CopySpread(symbol, period, 1, 1, spread) != 1)
               { /*mControl.Check();*/
                  ;
               }
               thisBar = (datetime)rates[0].time;
            }

            if (lastBar != thisBar)
            {
               if (lastBar != 0) // skip first price data after startup/reset
               {
                  if (chartTF == "TICK")
                  {
                     Data[0] = (long)tick.time_msc;
                     Data[1] = (double)tick.bid;
                     Data[2] = (double)tick.ask;
                  }
                  else
                  {
                     Data[0] = (long)rates[0].time;
                     Data[1] = (double)rates[0].open;
                     Data[2] = (double)rates[0].high;
                     Data[3] = (double)rates[0].low;
                     Data[4] = (double)rates[0].close;
                     Data[5] = (double)rates[0].tick_volume;
                  }
                  last["status"] = (string) "CONNECTED";
                  last["symbol"] = (string)symbol;
                  last["timeframe"] = (string)chartTF;
                  last["data"].Set(Data);

                  string t = last.Serialize();
                  if (debug)
                     Print(t);
                  InformClientSocket(sysSocket, t);
                  symbolSubscriptions[i].lastBar = thisBar;
               }
               else
                  symbolSubscriptions[i].lastBar = thisBar;
            }
         }
      }
      else
      {
         // send disconnect message only once
         if (connectedFlag)
         {
            last["status"] = (string) "DISCONNECTED";
            string t = last.Serialize();
            if (debug)
               Print(t);
            InformClientSocket(sysSocket, t);
            connectedFlag = false;
         }
      }
   }
}

//+------------------------------------------------------------------+
//| Expert timer function                                            |
//+------------------------------------------------------------------+
void OnTimer()
{
   tm = TimeTradeServer();
   // Stream live price data
   StreamPriceData();

   ZmqMsg request;

   // Get request from client via System socket.
   sysSocket.recv(request, true);

   // Request recived
   if (request.size() > 0)
   {
      // Pull request to RequestHandler().
      RequestHandler(request);
   }
}
//+------------------------------------------------------------------+
//| Request handler                                                  |
//+------------------------------------------------------------------+
void RequestHandler(ZmqMsg &request)
{

   CJAVal incomingMessage;

   ResetLastError();
   // Get data from reguest
   string msg = request.getData();

   if (debug)
      Print("Processing:" + msg);

   if (!incomingMessage.Deserialize(msg))
   {
      mControl.mSetUserError(65537, GetErrorID(65537));
      CheckError(__FUNCTION__);
   }

   // Send response to System socket that request was received
   // Some historical data requests can take a lot of time

   // Process action command
   string action = incomingMessage["action"].ToStr();

   if (action == "CONFIG")
      ScriptConfiguration(incomingMessage);
   else if (action == "ACCOUNT")
      GetAccountInfo();
   else if (action == "LISTSYMBOLS")
      GetSymbolsList();
   else if (action == "LISTSYMBOLSDETAILS")
      GetSymbolsListDetails();
   else if (action == "BALANCE")
      GetBalanceInfo();
   else if (action == "HISTORY")
      HistoryInfo(incomingMessage);
   else if (action == "TRADE")
      TradingModule(incomingMessage);
   else if (action == "POSITIONS")
      GetPositions(incomingMessage);
   else if (action == "ORDERS")
      GetOrders(incomingMessage);
   else if (action == "CALENDAR")
      GetEconomicCalendar(incomingMessage);
   else
   {
      mControl.mSetUserError(65538, GetErrorID(65538));
      CheckError(__FUNCTION__);
   }
}

//+------------------------------------------------------------------+
//| Reconfigure the script params                                    |
//+------------------------------------------------------------------+
void ScriptConfiguration(CJAVal &dataObject)
{

   string symbol = dataObject["symbol"].ToStr();
   string chartTF = dataObject["chartTF"].ToStr();

   ArrayResize(symbolSubscriptions, symbolSubscriptionCount + 1);
   symbolSubscriptions[symbolSubscriptionCount].symbol = symbol;
   symbolSubscriptions[symbolSubscriptionCount].chartTf = chartTF;
   // to initialze with value 0 skips the first price
   symbolSubscriptions[symbolSubscriptionCount].lastBar = 0;
   symbolSubscriptionCount++;

   mControl.mResetLastError();
   SymbolInfoString(symbol, SYMBOL_DESCRIPTION);
   if (!CheckError(__FUNCTION__))
      ActionDoneOrError(ERR_SUCCESS, __FUNCTION__, "ERR_SUCCESS");
}

//+------------------------------------------------------------------+
//| Account information                                              |
//+------------------------------------------------------------------+
void GetAccountInfo()
{

   CJAVal info;

   info["error"] = false;
   info["broker"] = AccountInfoString(ACCOUNT_COMPANY);
   info["currency"] = AccountInfoString(ACCOUNT_CURRENCY);
   info["server"] = AccountInfoString(ACCOUNT_SERVER);
   info["trading_allowed"] = TerminalInfoInteger(TERMINAL_TRADE_ALLOWED);
   info["bot_trading"] = AccountInfoInteger(ACCOUNT_TRADE_EXPERT);
   info["balance"] = AccountInfoDouble(ACCOUNT_BALANCE);
   info["equity"] = AccountInfoDouble(ACCOUNT_EQUITY);
   info["margin"] = AccountInfoDouble(ACCOUNT_MARGIN);
   info["margin_free"] = AccountInfoDouble(ACCOUNT_MARGIN_FREE);
   info["margin_level"] = AccountInfoDouble(ACCOUNT_MARGIN_LEVEL);
   info["time"] = string(tm); // sendig time to ejtraderMT for localtime dataframe

   string t = info.Serialize();
   if (debug)
      Print(t);
   InformClientSocket(sysSocket, t);
}

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void GetEconomicCalendar(CJAVal &dataObject)
{
   string actionType = dataObject["actionType"].ToStr();

   string symbol = dataObject["symbol"].ToStr();
   if (actionType == "DATA")
   {

      CJAVal data, d;

      datetime fromDate = (datetime)dataObject["fromDate"].ToInt();
      datetime toDate = TimeCurrent();
      if (dataObject["toDate"].ToInt() != NULL)
         toDate = (datetime)dataObject["toDate"].ToInt();

      CALENDAR Calendar;
      string Currencies[2];
      int Size;

      if (symbol == NULL)
      {
         Size = Calendar.Set(NULL, CALENDAR_IMPORTANCE_NONE, TimeToString(fromDate, TIME_DATE), TimeToString(toDate, TIME_DATE));
      }
      else
      {
         Currencies[0] = ::SymbolInfoString(symbol, SYMBOL_CURRENCY_BASE);
         Currencies[1] = ::SymbolInfoString(symbol, SYMBOL_CURRENCY_PROFIT);
         // Tomou eventos para todas as moedas (NULL) começando com a menor (NENHUM).
         Size = Calendar.Set(Currencies, CALENDAR_IMPORTANCE_NONE, TimeToString(fromDate, TIME_DATE), TimeToString(toDate, TIME_DATE));
      }

      if (Size)
      {
         for (int i = 0; i < Size; i++)
         {
            string string_split = Calendar[i].ToString();

            string sep = "|"; // A separator as a character
            ushort u_sep;     // The code of the separator character
            string result[];
            u_sep = StringGetCharacter(sep, 0);
            //--- Split the string to substrings
            int k = StringSplit(string_split, u_sep, result);
            for (int b = 0; b < k; b++)
            {

               data[i][b] = result[b];
            }
         }
         d["data"].Set(data);
      }
      else
      {
         d["data"].Add(data);
      }
      Print("Finished preparing Calender data");

      string t = d.Serialize();
      if (debug)
         Print(t);
      InformClientSocket(sysSocket, t);
   }
}
//+------------------------------------------------------------------+
//| Balance information                                              |
//+------------------------------------------------------------------+
void GetBalanceInfo()
{

   CJAVal info;
   info["balance"] = AccountInfoDouble(ACCOUNT_BALANCE);
   info["equity"] = AccountInfoDouble(ACCOUNT_EQUITY);
   info["margin"] = AccountInfoDouble(ACCOUNT_MARGIN);
   info["margin_free"] = AccountInfoDouble(ACCOUNT_MARGIN_FREE);

   string t = info.Serialize();
   if (debug)
      Print(t);
   InformClientSocket(sysSocket, t);
}

//+------------------------------------------------------------------+
//| All symbol list                                                  |
//+------------------------------------------------------------------+
void GetSymbolsList()
{
   CJAVal Data;
   CJAVal info;
   int total = SymbolsTotal(false);
   for (int i = 0; i < total; i++)
   {
      CJAVal symbolData;
      string symbol = SymbolName(i, false);
      symbolData = (string)symbol;
      Data.Add(symbolData);
   }

   info["status"] = (string) "CONNECTED";
   info["data"].Set(Data);

   string t = info.Serialize();
   if (debug)
      Print(t);
   InformClientSocket(sysSocket, t);
}

//+------------------------------------------------------------------+
//| All symbol list and Details                                      |
//+------------------------------------------------------------------+
void GetSymbolsListDetails()
{
   CJAVal Data;
   CJAVal info;
   int total = SymbolsTotal(false);
   for (int i = 0; i < total; i++)
   {
      CJAVal symbolData;
      string symbol = SymbolName(i, false);
      string assetType = GetSymbolType(symbol);
      datetime expiration;
      string expiration_date = "No Expiration";
      if (SymbolInfoInteger(symbol, SYMBOL_EXPIRATION_TIME, expiration) && StringLen(string(expiration)) != 0)
      {
         expiration_date = TimeToString(expiration, TIME_DATE | TIME_MINUTES);
      }
      symbolData[0] = (string)symbol;
      symbolData[1] = (string)expiration_date;
      symbolData[2] = (string)assetType;

      Data.Add(symbolData);
   }

   info["status"] = (string) "CONNECTED";
   info["data"].Set(Data);

   string t = info.Serialize();
   if (debug)
      Print(t);
   InformClientSocket(sysSocket, t);
}

//+------------------------------------------------------------------+
//| Push historical data to ZMQ socket                               |
//+------------------------------------------------------------------+
bool PushHistoricalData(CJAVal &data)
{
   string t = data.Serialize();
   if (debug)
      Print(t);
   InformClientSocket(sysSocket, t);
   return true;
}

//+------------------------------------------------------------------+
//| Convert chart timeframe from string to enum                      |
//+------------------------------------------------------------------+
ENUM_TIMEFRAMES GetTimeframe(string chartTF)
{

   ENUM_TIMEFRAMES tf;
   tf = NULL;

   if (chartTF == "TICK")
      tf = PERIOD_CURRENT;

   if (chartTF == "M1")
      tf = PERIOD_M1;

   if (chartTF == "M5")
      tf = PERIOD_M5;

   if (chartTF == "M15")
      tf = PERIOD_M15;

   if (chartTF == "M30")
      tf = PERIOD_M30;

   if (chartTF == "H1")
      tf = PERIOD_H1;

   if (chartTF == "H2")
      tf = PERIOD_H2;

   if (chartTF == "H3")
      tf = PERIOD_H3;

   if (chartTF == "H4")
      tf = PERIOD_H4;

   if (chartTF == "H6")
      tf = PERIOD_H6;

   if (chartTF == "H8")
      tf = PERIOD_H8;

   if (chartTF == "H12")
      tf = PERIOD_H12;

   if (chartTF == "D1")
      tf = PERIOD_D1;

   if (chartTF == "W1")
      tf = PERIOD_W1;

   if (chartTF == "MN1")
      tf = PERIOD_MN1;

   // if tf == NULL an error will be raised in config function
   return (tf);
}

//+------------------------------------------------------------------+
//| Trade confirmation                                               |
//+------------------------------------------------------------------+
void OrderDoneOrError(bool error, string funcName, CTrade &trade)
{

   CJAVal conf;

   conf["error"] = (bool)error;
   conf["retcode"] = (int)trade.ResultRetcode();
   conf["desription"] = (string)GetRetcodeID(trade.ResultRetcode());
   // conf["deal"]=(int) trade.ResultDeal();
   conf["order"] = (int)trade.ResultOrder();
   conf["volume"] = (double)trade.ResultVolume();
   conf["price"] = (double)trade.ResultPrice();
   conf["bid"] = (double)trade.ResultBid();
   conf["ask"] = (double)trade.ResultAsk();
   conf["function"] = (string)funcName;

   string t = conf.Serialize();
   if (debug)
      Print(t);
   InformClientSocket(sysSocket, t);
}

//+------------------------------------------------------------------+
//| Error reporting                                                  |
//+------------------------------------------------------------------+
bool CheckError(string funcName)
{
   int lastError = mControl.mGetLastError();
   if (lastError)
   {
      string desc = mControl.mGetDesc();
      if (debug)
         Print("Error handling source: ", funcName, " description: ", desc);
      Print("Error handling source: ", funcName, " description: ", desc);
      mControl.Check();
      ActionDoneOrError(lastError, funcName, desc);
      return true;
   }
   else
      return false;
}

//+------------------------------------------------------------------+
//| Action confirmation                                              |
//+------------------------------------------------------------------+
void ActionDoneOrError(int lastError, string funcName, string desc)
{

   CJAVal conf;

   conf["error"] = (bool)true;
   if (lastError == 0)
      conf["error"] = (bool)false;

   conf["lastError"] = (string)lastError;
   conf["description"] = (string)desc;
   conf["function"] = (string)funcName;

   string t = conf.Serialize();
   if (debug)
      Print(t);
   InformClientSocket(sysSocket, t);
}

//+------------------------------------------------------------------+
//| Inform Client via socket                                         |
//+------------------------------------------------------------------+
void InformClientSocket(Socket &workingSocket, string replyMessage)
{

   // non-blocking
   workingSocket.send(replyMessage, true);
   // TODO: Array out of range error
   mControl.mResetLastError();
   // mControl.Check();
}

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
string GetSymbolType(string symbol)
{
   long type;
   if (SymbolInfoInteger(symbol, SYMBOL_TRADE_CALC_MODE, type))
   {
      switch (int(type))
      {

      case SYMBOL_CALC_MODE_CFD:
         return ("CFD");
         break;

      case SYMBOL_CALC_MODE_CFDINDEX:
         return ("INDEX");
         break;

      case SYMBOL_CALC_MODE_CFDLEVERAGE:
         return ("CFD LEVERAGE");
         break;
      case SYMBOL_CALC_MODE_FOREX:
         return ("FOREX");
         break;

      case SYMBOL_CALC_MODE_FOREX_NO_LEVERAGE:
         return ("FOREX_NO_LEVERAGE");
         break;

      case SYMBOL_CALC_MODE_SERV_COLLATERAL:
         return ("SERV_COLLATERAL");
         break;

      case SYMBOL_CALC_MODE_EXCH_BONDS:
         return ("BONDS");
         break;

      case SYMBOL_CALC_MODE_EXCH_BONDS_MOEX:
         return ("BONDS_MOEX");
         break;

      case SYMBOL_CALC_MODE_EXCH_FUTURES:
         return ("FUTURES");
         break;

      case SYMBOL_CALC_MODE_EXCH_FUTURES_FORTS:
         return ("FUTURES_FORTS");
         break;

      case SYMBOL_CALC_MODE_EXCH_OPTIONS_MARGIN:
         return ("OPTIONS_MARGIN");
         break;

      case SYMBOL_CALC_MODE_EXCH_STOCKS:
         return ("STOCKS");
         break;

      case SYMBOL_CALC_MODE_EXCH_STOCKS_MOEX:
         return ("STOCKS_MOEX");
         break;

      default:
         return ("OPTIONS");
      }
   }
   else
   {
      return ("Unknown");
   }
}
//+------------------------------------------------------------------+
//| Get retcode message by retcode id                                |
//+------------------------------------------------------------------+
string GetRetcodeID(int retcode)
{

   switch (retcode)
   {
   case 10004:
      return ("TRADE_RETCODE_REQUOTE");
      break;
   case 10006:
      return ("TRADE_RETCODE_REJECT");
      break;
   case 10007:
      return ("TRADE_RETCODE_CANCEL");
      break;
   case 10008:
      return ("TRADE_RETCODE_PLACED");
      break;
   case 10009:
      return ("TRADE_RETCODE_DONE");
      break;
   case 10010:
      return ("TRADE_RETCODE_DONE_PARTIAL");
      break;
   case 10011:
      return ("TRADE_RETCODE_ERROR");
      break;
   case 10012:
      return ("TRADE_RETCODE_TIMEOUT");
      break;
   case 10013:
      return ("TRADE_RETCODE_INVALID");
      break;
   case 10014:
      return ("TRADE_RETCODE_INVALID_VOLUME");
      break;
   case 10015:
      return ("TRADE_RETCODE_INVALID_PRICE");
      break;
   case 10016:
      return ("TRADE_RETCODE_INVALID_STOPS");
      break;
   case 10017:
      return ("TRADE_RETCODE_TRADE_DISABLED");
      break;
   case 10018:
      return ("TRADE_RETCODE_MARKET_CLOSED");
      break;
   case 10019:
      return ("TRADE_RETCODE_NO_MONEY");
      break;
   case 10020:
      return ("TRADE_RETCODE_PRICE_CHANGED");
      break;
   case 10021:
      return ("TRADE_RETCODE_PRICE_OFF");
      break;
   case 10022:
      return ("TRADE_RETCODE_INVALID_EXPIRATION");
      break;
   case 10023:
      return ("TRADE_RETCODE_ORDER_CHANGED");
      break;
   case 10024:
      return ("TRADE_RETCODE_TOO_MANY_REQUESTS");
      break;
   case 10025:
      return ("TRADE_RETCODE_NO_CHANGES");
      break;
   case 10026:
      return ("TRADE_RETCODE_SERVER_DISABLES_AT");
      break;
   case 10027:
      return ("TRADE_RETCODE_CLIENT_DISABLES_AT");
      break;
   case 10028:
      return ("TRADE_RETCODE_LOCKED");
      break;
   case 10029:
      return ("TRADE_RETCODE_FROZEN");
      break;
   case 10030:
      return ("TRADE_RETCODE_INVALID_FILL");
      break;
   case 10031:
      return ("TRADE_RETCODE_CONNECTION");
      break;
   case 10032:
      return ("TRADE_RETCODE_ONLY_REAL");
      break;
   case 10033:
      return ("TRADE_RETCODE_LIMIT_ORDERS");
      break;
   case 10034:
      return ("TRADE_RETCODE_LIMIT_VOLUME");
      break;
   case 10035:
      return ("TRADE_RETCODE_INVALID_ORDER");
      break;
   case 10036:
      return ("TRADE_RETCODE_POSITION_CLOSED");
      break;
   case 10038:
      return ("TRADE_RETCODE_INVALID_CLOSE_VOLUME");
      break;
   case 10039:
      return ("TRADE_RETCODE_CLOSE_ORDER_EXIST");
      break;
   case 10040:
      return ("TRADE_RETCODE_LIMIT_POSITIONS");
      break;
   case 10041:
      return ("TRADE_RETCODE_REJECT_CANCEL");
      break;
   case 10042:
      return ("TRADE_RETCODE_LONG_ONLY");
      break;
   case 10043:
      return ("TRADE_RETCODE_SHORT_ONLY");
      break;
   case 10044:
      return ("TRADE_RETCODE_CLOSE_ONLY");
      break;

   default:
      return ("TRADE_RETCODE_UNKNOWN=" + IntegerToString(retcode));
      break;
   }
}

//+------------------------------------------------------------------+
//| Get error message by error id                                    |
//+------------------------------------------------------------------+
string GetErrorID(int error)
{

   switch (error)
   {
   // Custom errors
   case 65537:
      return ("ERR_DESERIALIZATION");
      break;
   case 65538:
      return ("ERR_WRONG_ACTION");
      break;
   case 65539:
      return ("ERR_WRONG_ACTION_TYPE");
      break;
   case 65540:
      return ("ERR_CLEAR_SUBSCRIPTIONS_FAILED");
      break;
   case 65541:
      return ("ERR_RETRIEVE_DATA_FAILED");
      break;
   case 65542:
      return ("ERR_CVS_FILE_CREATION_FAILED");
      break;

   default:
      return ("ERR_CODE_UNKNOWN=" + IntegerToString(error));
      break;
   }
}

//+------------------------------------------------------------------+
//| Return a textual description of the deinitialization reason code |
//+------------------------------------------------------------------+
string getUninitReasonText(int reasonCode)
{
   string text = "";
   //---
   switch (reasonCode)
   {
   case REASON_ACCOUNT:
      text = "Account was changed";
      break;
   case REASON_CHARTCHANGE:
      text = "Symbol or timeframe was changed";
      break;
   case REASON_CHARTCLOSE:
      text = "Chart was closed";
      break;
   case REASON_PARAMETERS:
      text = "Input-parameter was changed";
      break;
   case REASON_RECOMPILE:
      text = "Program " + __FILE__ + " was recompiled";
      break;
   case REASON_REMOVE:
      text = "Program " + __FILE__ + " was removed from chart";
      break;
   case REASON_TEMPLATE:
      text = "New template was applied to chart";
      break;
   default:
      text = "Another reason";
   }
   //---
   return text;
}
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
