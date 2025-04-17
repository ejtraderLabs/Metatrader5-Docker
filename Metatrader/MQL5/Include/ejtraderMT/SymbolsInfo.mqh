#property copyright "ejtrader"
#property link      "https://github.com/ejtraderLabs/MQL5-ejtraderMT"

void GetSymbolInfo(CJAVal &dataObject)
{
   mControl.mResetLastError();
   
   string symbol = dataObject["symbol"].ToStr();
   
   // Validate symbol exists
   if(!SymbolSelect(symbol, true)) {
      ActionDoneOrError(ERR_MARKET_UNKNOWN_SYMBOL, __FUNCTION__, "Symbol not found");
      return;
   }
   
   CJAVal symbolData;
   
   // General symbol info
   symbolData["name"] = symbol;
   symbolData["description"] = SymbolInfoString(symbol, SYMBOL_DESCRIPTION);
   symbolData["currency_base"] = SymbolInfoString(symbol, SYMBOL_CURRENCY_BASE);
   symbolData["currency_profit"] = SymbolInfoString(symbol, SYMBOL_CURRENCY_PROFIT);
   symbolData["currency_margin"] = SymbolInfoString(symbol, SYMBOL_CURRENCY_MARGIN);
   symbolData["type"] = GetSymbolType(symbol);
   
   // Prices
   symbolData["bid"] = SymbolInfoDouble(symbol, SYMBOL_BID);
   symbolData["ask"] = SymbolInfoDouble(symbol, SYMBOL_ASK);
   symbolData["last"] = SymbolInfoDouble(symbol, SYMBOL_LAST);
   symbolData["bidhigh"] = SymbolInfoDouble(symbol, SYMBOL_BIDHIGH);
   symbolData["bidlow"] = SymbolInfoDouble(symbol, SYMBOL_BIDLOW);
   symbolData["askhigh"] = SymbolInfoDouble(symbol, SYMBOL_ASKHIGH);
   symbolData["asklow"] = SymbolInfoDouble(symbol, SYMBOL_ASKLOW);
   symbolData["lasthigh"] = SymbolInfoDouble(symbol, SYMBOL_LASTHIGH);
   symbolData["lastlow"] = SymbolInfoDouble(symbol, SYMBOL_LASTLOW);
   
   // Trade parameters
   symbolData["point"] = SymbolInfoDouble(symbol, SYMBOL_POINT);
   symbolData["digits"] = SymbolInfoInteger(symbol, SYMBOL_DIGITS);
   symbolData["spread"] = SymbolInfoInteger(symbol, SYMBOL_SPREAD);
   symbolData["spread_float"] = (bool)SymbolInfoInteger(symbol, SYMBOL_SPREAD_FLOAT);
   
   // Volume specifications
   symbolData["volume_min"] = SymbolInfoDouble(symbol, SYMBOL_VOLUME_MIN);
   symbolData["volume_max"] = SymbolInfoDouble(symbol, SYMBOL_VOLUME_MAX);
   symbolData["volume_step"] = SymbolInfoDouble(symbol, SYMBOL_VOLUME_STEP);
   symbolData["volume_limit"] = SymbolInfoDouble(symbol, SYMBOL_VOLUME_LIMIT);
   
   // Contract and tick info
   symbolData["trade_contract_size"] = SymbolInfoDouble(symbol, SYMBOL_TRADE_CONTRACT_SIZE);
   symbolData["trade_tick_size"] = SymbolInfoDouble(symbol, SYMBOL_TRADE_TICK_SIZE);
   symbolData["trade_tick_value"] = SymbolInfoDouble(symbol, SYMBOL_TRADE_TICK_VALUE);
   symbolData["trade_tick_value_profit"] = SymbolInfoDouble(symbol, SYMBOL_TRADE_TICK_VALUE_PROFIT);
   symbolData["trade_tick_value_loss"] = SymbolInfoDouble(symbol, SYMBOL_TRADE_TICK_VALUE_LOSS);
   
   // Trade modes and restrictions
   symbolData["trade_mode"] = EnumToString((ENUM_SYMBOL_TRADE_MODE)SymbolInfoInteger(symbol, SYMBOL_TRADE_MODE));
   symbolData["trade_calc_mode"] = EnumToString((ENUM_SYMBOL_CALC_MODE)SymbolInfoInteger(symbol, SYMBOL_TRADE_CALC_MODE));
   symbolData["trade_stops_level"] = SymbolInfoInteger(symbol, SYMBOL_TRADE_STOPS_LEVEL);
   symbolData["trade_freeze_level"] = SymbolInfoInteger(symbol, SYMBOL_TRADE_FREEZE_LEVEL);
   
   // Swap information
   symbolData["swap_mode"] = EnumToString((ENUM_SYMBOL_SWAP_MODE)SymbolInfoInteger(symbol, SYMBOL_SWAP_MODE));
   symbolData["swap_long"] = SymbolInfoDouble(symbol, SYMBOL_SWAP_LONG);
   symbolData["swap_short"] = SymbolInfoDouble(symbol, SYMBOL_SWAP_SHORT);
   symbolData["swap_rollover3days"] = EnumToString((ENUM_DAY_OF_WEEK)SymbolInfoInteger(symbol, SYMBOL_SWAP_ROLLOVER3DAYS));
   
   // Session information
   symbolData["session_volume"] = SymbolInfoDouble(symbol, SYMBOL_SESSION_VOLUME);
   symbolData["session_deals"] = SymbolInfoInteger(symbol, SYMBOL_SESSION_DEALS);
   symbolData["session_buy_orders"] = SymbolInfoInteger(symbol, SYMBOL_SESSION_BUY_ORDERS);
   symbolData["session_sell_orders"] = SymbolInfoInteger(symbol, SYMBOL_SESSION_SELL_ORDERS);
   
   // Expiration information (for futures/options)
   datetime expiration;
   if(SymbolInfoInteger(symbol, SYMBOL_EXPIRATION_TIME, expiration))
      symbolData["expiration_time"] = (long)expiration;
   
   symbolData["margin_hedged_use_leg"] = (bool)SymbolInfoInteger(symbol, SYMBOL_MARGIN_HEDGED_USE_LEG);
   symbolData["order_mode"] = SymbolInfoInteger(symbol, SYMBOL_ORDER_MODE);
   symbolData["filling_mode"] = SymbolInfoInteger(symbol, SYMBOL_FILLING_MODE);
   
   CJAVal data;
   data["error"] = (bool)false;
   data["data"].Set(symbolData);
   
   string t = data.Serialize();
   if(debug)
      Print(t);
   InformClientSocket(sysSocket, t);
}
