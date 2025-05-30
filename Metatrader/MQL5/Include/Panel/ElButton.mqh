//+------------------------------------------------------------------+
//|                                                   MBookPanel.mqh |
//|                        Copyright 2015, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Copyright 2015, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"
#include "Node.mqh"
#include "ElChart.mqh"
//+------------------------------------------------------------------+
//| Button state                                                     |
//+------------------------------------------------------------------+
enum ENUM_BTN_PUSH_STATE
{
   PUSH_OFF,             // The button is unpressed (off)
   PUSH_ON               // The button is pressed (on)
};

//+------------------------------------------------------------------+
//| Button                                                           |
//+------------------------------------------------------------------+
class CElButton : public CElChart
{
public:
   CElButton(void);
   ENUM_BTN_PUSH_STATE State(void);
   bool State(ENUM_BTN_PUSH_STATE state);
   virtual void OnPushButton(ENUM_BTN_PUSH_STATE state);
};
//+------------------------------------------------------------------+
//| By default the button displays the same text both when           |
//| pressed and when unpressed.                                      |
//+------------------------------------------------------------------+
CElButton::CElButton() : CElChart(OBJ_BUTTON)
{
   BorderType(BORDER_RAISED);
   BackgroundColor(clrWhiteSmoke);
   
}
//+------------------------------------------------------------------+
//| Returns PUSH_OFF when the button is unpressed or not displayed,  |
//| Returns PUSH_ON when the button is pressed.                           |
//+------------------------------------------------------------------+
ENUM_BTN_PUSH_STATE CElButton::State(void)
{
   if(!IsShowed())
      return PUSH_OFF;
   if(ObjectGetInteger(ChartID(), m_name, OBJPROP_STATE) > 0)
      return PUSH_ON;
   return PUSH_OFF;
}
//+------------------------------------------------------------------+
//| Sets button state                                                |
//+------------------------------------------------------------------+
bool CElButton::State(ENUM_BTN_PUSH_STATE state)
{
   bool isPush = state == PUSH_ON ? true : false;
   bool res = false;
   if(IsShowed())
   {
      if(ObjectSetInteger(ChartID(), m_name, OBJPROP_STATE, isPush))
      {
         OnPushButton(state);
         return true;
      }
      return false;
   }
   return false;
}
//+------------------------------------------------------------------+
//| The event performed once a button is pressed                     |
//+------------------------------------------------------------------+
void CElButton::OnPushButton(ENUM_BTN_PUSH_STATE state)
{
}