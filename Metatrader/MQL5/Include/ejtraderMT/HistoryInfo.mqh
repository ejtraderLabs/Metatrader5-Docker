
//+------------------------------------------------------------------+
#property copyright "ejtrader"
#property link      "https://github.com/ejtraderLabs/MQL5-ejtraderMT"
//+------------------------------------------------------------------+
//| Get historical data                                              |
//+------------------------------------------------------------------+
void HistoryInfo(CJAVal &dataObject)
  {

   string actionType = dataObject["actionType"].ToStr();
   string chartTF = dataObject["chartTF"].ToStr();
   string symbol=dataObject["symbol"].ToStr();
// bool correctTickHistory=dataObject["correctTickHistory"].ToBool();

// Write CVS fle to local directory
   if(actionType=="WRITE" && chartTF=="TICK")
     {

      CJAVal data, d, msg;
      MqlTick tickArray[];
      string fileName=symbol + "-" + chartTF + ".csv";  // file name
      string directoryName="Data"; // directory name
      string outputFile=directoryName+"\\"+fileName;

      ENUM_TIMEFRAMES period=GetTimeframe(chartTF);
      datetime fromDate=(datetime)dataObject["fromDate"].ToInt();
      datetime toDate=TimeCurrent();
      if(dataObject["toDate"].ToInt()!=NULL)
         toDate=(datetime)dataObject["toDate"].ToInt();

      Print("Fetching HISTORY");
      Print("1) Symbol: "+symbol);
      Print("2) Timeframe: Ticks");
      Print("3) Date from: "+TimeToString(fromDate));
      if(dataObject["toDate"].ToInt()!=NULL)
         Print("4) Date to:"+TimeToString(toDate));

      int tickCount = 0;
      ulong fromDateM = StringToTime(fromDate);
      ulong toDateM = StringToTime(toDate);

      tickCount=CopyTicksRange(symbol,tickArray,COPY_TICKS_ALL,1000*(ulong)fromDateM,1000*(ulong)toDateM);
      if(tickCount < 0)
        {
         mControl.mSetUserError(65541, GetErrorID(65541));
        }
      CheckError(__FUNCTION__);

      Print("Preparing data of ", tickCount, " ticks for ", symbol);
      /*
      if(correctTickHistory)
         CorrectTicks(symbol,tickArray);
      */
      int file_handle=FileOpen(outputFile, FILE_WRITE | FILE_CSV);
      if(file_handle!=INVALID_HANDLE)
        {
         msg["status"] = (string) "CONNECTED";
         msg["type"] = (string) "NORMAL";
         msg["data"] = (string) StringFormat("Writing to: %s\\%s", TerminalInfoString(TERMINAL_DATA_PATH), outputFile);
         if(liveStream)
            InformClientSocket(liveSocket, msg.Serialize());
         ActionDoneOrError(ERR_SUCCESS, __FUNCTION__, "ERR_SUCCESS");
         //ActionDoneOrError(ERR_SUCCESS  , __FUNCTION__, dataSocket);
         // Inform client that file is avalable for writing

         PrintFormat("%s file is available for writing",fileName);
         PrintFormat("File path: %s\\Files\\",TerminalInfoString(TERMINAL_DATA_PATH));
         //--- write the time and values of signals to the file
         for(int i=0; i<tickCount; i++)
           {
            FileWrite(file_handle,tickArray[i].time_msc, ",", tickArray[i].bid, ",", tickArray[i].ask);
            msg["status"] = (string) "CONNECTED";
            msg["type"] = (string) "FLUSH";
            msg["data"] = (string) tickArray[i].time_msc;
            if(liveStream)
               InformClientSocket(liveSocket, msg.Serialize());
           }
         //--- close the file
         FileClose(file_handle);
         PrintFormat("Data is written, %s file is closed",fileName);
         msg["status"] = (string) "DISCONNECTED";
         msg["type"] = (string) "NORMAL";
         msg["data"] = (string) StringFormat("Writing to: %s\\%s", outputFile, " is finished");
         if(liveStream)
            InformClientSocket(liveSocket, msg.Serialize());

        }
      else
        {
         // File is not available for writing
         mControl.mSetUserError(65542, GetErrorID(65542));
         CheckError(__FUNCTION__);
        }
      connectedFlag=false;
     }

// Write CVS fle to local directory
   else
      if(actionType=="WRITE" && chartTF!="TICK")
        {

         CJAVal c, d;
         MqlRates r[];
         int spread[];
         string fileName=symbol + "-" + chartTF + ".csv";  // file name
         string directoryName="Data"; // directory name
         string outputFile=directoryName+"//"+fileName;

         int barCount;
         ENUM_TIMEFRAMES period=GetTimeframe(chartTF);
         datetime fromDate=(datetime)dataObject["fromDate"].ToInt();
         datetime toDate=TimeCurrent();
         if(dataObject["toDate"].ToInt()!=NULL)
            toDate=(datetime)dataObject["toDate"].ToInt();

         Print("Fetching HISTORY");
         Print("1) Symbol :"+symbol);
         Print("2) Timeframe :"+EnumToString(period));
         Print("3) Date from :"+TimeToString(fromDate));
         if(dataObject["toDate"].ToInt()!=NULL)
            Print("4) Date to:"+TimeToString(toDate));

         barCount=CopyRates(symbol,period,fromDate,toDate,r);
         if(CopySpread(symbol,period, fromDate, toDate, spread)!=1)
           {
            mControl.mSetUserError(65541, GetErrorID(65541));
           }

         Print("Preparing tick data of ", barCount, " ticks for ", symbol);
         int file_handle=FileOpen(outputFile, FILE_WRITE | FILE_CSV);
         if(file_handle!=INVALID_HANDLE)
           {
            ActionDoneOrError(ERR_SUCCESS, __FUNCTION__, "ERR_SUCCESS");;
            PrintFormat("%s file is available for writing",outputFile);
            PrintFormat("File path: %s\\Files\\",TerminalInfoString(TERMINAL_DATA_PATH));
            //--- write the time and values of signals to the file
            for(int i=0; i<barCount; i++)
               FileWrite(file_handle,r[i].time, ",", r[i].open, ",", r[i].high, ",", r[i].low, ",", r[i].close, ",", r[i].tick_volume, spread[i]);
            //--- close the file
            FileClose(file_handle);
            PrintFormat("Data is written, %s file is closed", outputFile);
           }
         else
           {
            mControl.mSetUserError(65542, GetErrorID(65542));
            CheckError(__FUNCTION__);
           }
        }

      else
         if(actionType=="DATA" && chartTF=="TICK")
           {

            CJAVal data, d;
            MqlTick tickArray[];

            ENUM_TIMEFRAMES period=GetTimeframe(chartTF);
            datetime fromDate=(datetime)dataObject["fromDate"].ToInt();
            datetime toDate=TimeCurrent();
            if(dataObject["toDate"].ToInt()!=NULL)
               toDate=(datetime)dataObject["toDate"].ToInt();

            if(debug)
              {
               Print("Fetching HISTORY");
               Print("1) Symbol: "+symbol);
               Print("2) Timeframe: Ticks");
               Print("3) Date from: "+TimeToString(fromDate));
               if(dataObject["toDate"].ToInt()!=NULL)
                  Print("4) Date to:"+TimeToString(toDate));
              }

            int tickCount = 0;
            ulong fromDateM = StringToTime(fromDate);
            ulong toDateM = StringToTime(toDate);

            tickCount=CopyTicksRange(symbol,tickArray, COPY_TICKS_ALL, 1000*(ulong)fromDateM, 1000*(ulong)toDateM);
            Print("Preparing tick data of ", tickCount, " ticks for ", symbol);
            /*
            if(correctTickHistory)
               CorrectTicks(symbol,tickArray);
            */
            if(tickCount)
              {
               for(int i=0; i<tickCount; i++)
                 {
                  data[i][0]=(long)   tickArray[i].time_msc;
                  data[i][1]=(double) tickArray[i].bid;
                  data[i][2]=(double) tickArray[i].ask;
                 }
               d["data"].Set(data);
              }
            else
              {
               d["data"].Add(data);
              }
            Print("Finished preparing tick data");

            d["symbol"]=symbol;
            d["timeframe"]=chartTF;

            PushHistoricalData(d);
           }

         else
            if(actionType=="DATA" && chartTF!="TICK")
              {

               CJAVal c, d;
               MqlRates r[];
               int spread[];
               int barCount=0;
               ENUM_TIMEFRAMES period=GetTimeframe(chartTF);
               datetime fromDate=(datetime)dataObject["fromDate"].ToInt();
               datetime toDate=TimeCurrent();
               if(dataObject["toDate"].ToInt()!=NULL)
                  toDate=(datetime)dataObject["toDate"].ToInt();

               if(debug)
                 {
                  Print("Fetching HISTORY");
                  Print("1) Symbol :"+symbol);
                  Print("2) Timeframe :"+EnumToString(period));
                  Print("3) Date from :"+TimeToString(fromDate));
                  if(dataObject["toDate"].ToInt()!=NULL)
                     Print("4) Date to:"+TimeToString(toDate));
                 }

               barCount=CopyRates(symbol, period, fromDate, toDate, r);
               if(CopySpread(symbol,period, fromDate, toDate, spread)!=1) { /*mControl.Check();*/ }

               if(barCount)
                 {
                  for(int i=0; i<barCount; i++)
                    {
                     c[i][0]=(long)   r[i].time;
                     c[i][1]=(double) r[i].open;
                     c[i][2]=(double) r[i].high;
                     c[i][3]=(double) r[i].low;
                     c[i][4]=(double) r[i].close;
                     c[i][5]=(double) r[i].tick_volume;
                     c[i][6]=(int) spread[i];
                    }
                  d["data"].Set(c);
                 }
               else
                 {
                  d["data"].Add(c);
                 }

               d["symbol"]=symbol;
               d["timeframe"]=chartTF;

               PushHistoricalData(d);
              }

            else
               if(actionType=="TRADES")
                 {
                  CDealInfo tradeInfo;
                  CJAVal trades, data;

                  if(HistorySelect(0,TimeCurrent()))
                    {
                     // Get total deals in history
                     int total = HistoryDealsTotal();
                     ulong ticket; // deal ticket

                     for(int i=0; i<total; i++)
                       {
                        if((ticket=HistoryDealGetTicket(i))>0)
                          {
                           tradeInfo.Ticket(ticket);
                           data["ticket"]=(long) tradeInfo.Ticket();
                           data["time"]=(long) tradeInfo.Time();
                           data["price"]=(double) tradeInfo.Price();
                           data["volume"]=(double) tradeInfo.Volume();
                           data["symbol"]=(string) tradeInfo.Symbol();
                           data["type"]=(string) tradeInfo.TypeDescription();
                           data["entry"]=(long) tradeInfo.Entry();
                           data["profit"]=(double) tradeInfo.Profit();

                           trades["trades"].Add(data);
                          }
                       }
                    }
                  else
                    {
                     trades["trades"].Add(data);
                    }

                  string t=trades.Serialize();
                  if(debug)
                     Print(t);
                  InformClientSocket(dataSocket,t);
                 }
               else
                 {
                  mControl.mSetUserError(65538, GetErrorID(65538));
                  CheckError(__FUNCTION__);
                 }
  }
//+------------------------------------------------------------------+
