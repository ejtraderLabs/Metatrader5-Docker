#property copyright "ejtrader"
#property link      "https://github.com/ejtraderLabs/MQL5-ejtraderMT"



//+------------------------------------------------------------------+
//| Fetch positions information                                      |
//+------------------------------------------------------------------+
void GetPositions(CJAVal &dataObject)
{
   CPositionInfo myposition;
   CJAVal data;
   
   // Get positions
   int positionsTotal=PositionsTotal();
   
   // Create empty array if no positions
   if(!positionsTotal)
   {
      CJAVal emptyPosition;
      data["positions"].Add(emptyPosition);
   }
      
   // Go through positions in a loop
   for(int i=0; i<positionsTotal; i++)
   {
      mControl.mResetLastError();
      
      // Get the position ticket directly by index
      ulong ticket = PositionGetTicket(i);
      
      if(ticket && myposition.SelectByTicket(ticket))
      {
         // Create a new position object for each position
         CJAVal position;
         
         position["id"]=PositionGetInteger(POSITION_IDENTIFIER);
         position["magic"]=PositionGetInteger(POSITION_MAGIC);
         position["symbol"]=PositionGetString(POSITION_SYMBOL);
         position["type"]=EnumToString(ENUM_POSITION_TYPE(PositionGetInteger(POSITION_TYPE)));
         position["time_setup"]=PositionGetInteger(POSITION_TIME);
         position["open"]=PositionGetDouble(POSITION_PRICE_OPEN);
         position["stoploss"]=PositionGetDouble(POSITION_SL);
         position["takeprofit"]=PositionGetDouble(POSITION_TP);
         position["volume"]=PositionGetDouble(POSITION_VOLUME);

         data["positions"].Add(position);
      }
      CheckError(__FUNCTION__);
   }
   
   data["error"]=(bool) false;
   
   string t=data.Serialize();
   if(debug)
      Print(t);
   InformClientSocket(sysSocket,t);
}

//+------------------------------------------------------------------+
//| Fetch orders information                                         |
//+------------------------------------------------------------------+
void GetOrders(CJAVal &dataObject)
  {
   mControl.mResetLastError();

   COrderInfo myorder;
   CJAVal data, order;

// Get orders
   if(HistorySelect(0,TimeCurrent()))
     {
      int ordersTotal = OrdersTotal();
      // Create empty array if no orders
      if(!ordersTotal)
        {
         data["error"]=(bool) false;
         data["orders"].Add(order);
        }

      for(int i=0; i<ordersTotal; i++)
        {
         if(myorder.Select(OrderGetTicket(i)))
           {
            order["id"]=(string) myorder.Ticket();
            order["magic"]=OrderGetInteger(ORDER_MAGIC);
            order["symbol"]=OrderGetString(ORDER_SYMBOL);
            order["type"]=EnumToString(ENUM_ORDER_TYPE(OrderGetInteger(ORDER_TYPE)));
            order["time_setup"]=OrderGetInteger(ORDER_TIME_SETUP);
            order["open"]=OrderGetDouble(ORDER_PRICE_OPEN);
            order["stoploss"]=OrderGetDouble(ORDER_SL);
            order["takeprofit"]=OrderGetDouble(ORDER_TP);
            order["volume"]=OrderGetDouble(ORDER_VOLUME_INITIAL);

            data["error"]=(bool) false;
            data["orders"].Add(order);
           }
         // Error handling
         CheckError(__FUNCTION__);
        }
     }

   string t=data.Serialize();
   if(debug)
      Print(t);
   InformClientSocket(sysSocket,t);
  }

//+------------------------------------------------------------------+
//| Trading module                                                   |
//+------------------------------------------------------------------+
void TradingModule(CJAVal &dataObject)
  {
   mControl.mResetLastError();
   CTrade trade;

   string   actionType = dataObject["actionType"].ToStr();
   string   symbol=dataObject["symbol"].ToStr();
   SymbolInfoString(symbol, SYMBOL_DESCRIPTION);
   CheckError(__FUNCTION__);

   int      idNimber=dataObject["id"].ToInt();
   double   volume=dataObject["volume"].ToDbl();
   double   SL=dataObject["stoploss"].ToDbl();
   double   TP=dataObject["takeprofit"].ToDbl();
   double   price=NormalizeDouble(dataObject["price"].ToDbl(),_Digits);
   double   deviation=dataObject["deviation"].ToDbl();
   string   comment=dataObject["comment"].ToStr();

// Order expiration section
   ENUM_ORDER_TYPE_TIME exp_type = ORDER_TIME_GTC;
   datetime expiration = 0;
   if(dataObject["expiration"].ToInt() != 0)
     {
      exp_type = ORDER_TIME_SPECIFIED;
      expiration=dataObject["expiration"].ToInt();
     }

// Market orders
   if(actionType=="ORDER_TYPE_BUY" || actionType=="ORDER_TYPE_SELL")
     {
      ENUM_ORDER_TYPE orderType=ORDER_TYPE_BUY;
      price = SymbolInfoDouble(symbol,SYMBOL_ASK);
      if(actionType=="ORDER_TYPE_SELL")
        {
         orderType=ORDER_TYPE_SELL;
         price=SymbolInfoDouble(symbol,SYMBOL_BID);
        }

      if(trade.PositionOpen(symbol,orderType,volume,price,SL,TP,comment))
        {
         OrderDoneOrError(false, __FUNCTION__, trade);
         return;
        }
     }

// Pending orders
   else
      if(actionType=="ORDER_TYPE_BUY_LIMIT" || actionType=="ORDER_TYPE_SELL_LIMIT" || actionType=="ORDER_TYPE_BUY_STOP" || actionType=="ORDER_TYPE_SELL_STOP")
        {
         if(actionType=="ORDER_TYPE_BUY_LIMIT")
           {
            if(trade.BuyLimit(volume,price,symbol,SL,TP,ORDER_TIME_GTC,expiration,comment))
              {
               OrderDoneOrError(false, __FUNCTION__, trade);
               return;
              }
           }
         else
            if(actionType=="ORDER_TYPE_SELL_LIMIT")
              {
               if(trade.SellLimit(volume,price,symbol,SL,TP,ORDER_TIME_GTC,expiration,comment))
                 {
                  OrderDoneOrError(false, __FUNCTION__, trade);
                  return;
                 }
              }
            else
               if(actionType=="ORDER_TYPE_BUY_STOP")
                 {
                  if(trade.BuyStop(volume,price,symbol,SL,TP,ORDER_TIME_GTC,expiration,comment))
                    {
                     OrderDoneOrError(false, __FUNCTION__, trade);
                     return;
                    }
                 }
               else
                  if(actionType=="ORDER_TYPE_SELL_STOP")
                    {
                     if(trade.SellStop(volume,price,symbol,SL,TP,ORDER_TIME_GTC,expiration,comment))
                       {
                        OrderDoneOrError(false, __FUNCTION__, trade);
                        return;
                       }
                    }
        }
      // Position modify
      else
         if(actionType=="POSITION_MODIFY")
           {
            if(trade.PositionModify(idNimber,SL,TP))
              {
               OrderDoneOrError(false, __FUNCTION__, trade);
               return;
              }
           }
         // Position close partial
         else
            if(actionType=="POSITION_PARTIAL")
              {
               if(trade.PositionClosePartial(idNimber,volume))
                 {
                  OrderDoneOrError(false, __FUNCTION__, trade);
                  return;
                 }
              }
            // Position close by id
            else
               if(actionType=="POSITION_CLOSE_ID")
                 {
                  if(trade.PositionClose(idNimber))
                    {
                     OrderDoneOrError(false, __FUNCTION__, trade);
                     return;
                    }
                 }
               // Position close by symbol
               else
                  if(actionType=="POSITION_CLOSE_SYMBOL")
                    {
                     if(trade.PositionClose(symbol))
                       {
                        OrderDoneOrError(false, __FUNCTION__, trade);
                        return;
                       }
                    }
                  // Modify pending order
                  else
                     if(actionType=="ORDER_MODIFY")
                       {
                        if(trade.OrderModify(idNimber,price,SL,TP,ORDER_TIME_GTC,expiration))
                          {
                           OrderDoneOrError(false, __FUNCTION__, trade);
                           return;
                          }
                       }
                     // Cancel pending order
                     else
                        if(actionType=="ORDER_CANCEL")
                          {
                           if(trade.OrderDelete(idNimber))
                             {
                              OrderDoneOrError(false, __FUNCTION__, trade);
                              return;
                             }
                          }
                        // Action type dosen't exist
                        else
                          {
                           mControl.mSetUserError(65538, GetErrorID(65538));
                           CheckError(__FUNCTION__);
                          }

// This part of the code runs if order was not completed
   OrderDoneOrError(true, __FUNCTION__, trade);
  }

//+------------------------------------------------------------------+
//| TradeTransaction function                                        |
//+------------------------------------------------------------------+
void OnTradeTransaction(const MqlTradeTransaction &trans,
                        const MqlTradeRequest &request,
                        const MqlTradeResult &result)
  {

   ENUM_TRADE_TRANSACTION_TYPE  trans_type=trans.type;
   switch(trans.type)
     {
      case  TRADE_TRANSACTION_REQUEST:
        {
         CJAVal data, req, res;

         req["action"]=EnumToString(request.action);
         req["order"]=(int) request.order;
         req["symbol"]=(string) request.symbol;
         req["volume"]=(double) request.volume;
         req["price"]=(double) request.price;
         req["stoplimit"]=(double) request.stoplimit;
         req["sl"]=(double) request.sl;
         req["tp"]=(double) request.tp;
         req["deviation"]=(int) request.deviation;
         req["type"]=EnumToString(request.type);
         req["type_filling"]=EnumToString(request.type_filling);
         req["type_time"]=EnumToString(request.type_time);
         req["expiration"]=(int) request.expiration;
         req["comment"]=(string) request.comment;
         req["position"]=(int) request.position;
         req["position_by"]=(int) request.position_by;

         res["retcode"]=(int) result.retcode;
         res["result"]=(string) GetRetcodeID(result.retcode);
         res["deal"]=(int) result.order;
         res["order"]=(int) result.order;
         res["volume"]=(double) result.volume;
         res["price"]=(double) result.price;
         res["comment"]=(string) result.comment;
         res["request_id"]=(int) result.request_id;
         res["retcode_external"]=(int) result.retcode_external;

         data["request"].Set(req);
         data["result"].Set(res);

         string t=data.Serialize();
         if(debug)
            Print(t);
         InformClientSocket(sysSocket,t);
        }
      break;
      default:
        {} break;
     }
  }
//+------------------------------------------------------------------+
