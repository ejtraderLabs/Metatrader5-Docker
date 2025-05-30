//+------------------------------------------------------------------+
//|                                                         Node.mqh |
//|                        Copyright 2015, MetaQuotes Software Corp. |
//|                                              http://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Copyright 2015, MetaQuotes Software Corp."
#property link      "http://www.mql5.com"

#include <Object.mqh>
#include <Arrays\ArrayObj.mqh>
#include "\Events\Event.mqh"
#include "\Events\EventChartObjClick.mqh"
#include "\Events\EventRefresh.mqh"
#define NAME_SIZE 8
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class CNode : public CObject
  {
private:
   static int        m_count;
protected:
   string            m_name;
   bool              m_showed;        // Object visibility flag
                     CNode();
   CArrayObj         m_elements;
   virtual void      OnShow();
   virtual void      OnHide();         // Hides child elements
   color             FontColor;
   virtual void      OnClick(void);
   virtual void      OnRefresh(CEventRefresh *event);
public:
   string            Name(void);
   bool              IsShowed(void);
   virtual void      Show();
   virtual void      Hide();           // Hides the current element
   virtual void      Event(CEvent *event);
   void              Tooltip(string message);
   void              AddElement(CNode* element);
                    ~CNode();
  };
static int CNode::m_count=0;
//+------------------------------------------------------------------+
//|  Generate random uniq name                                       |
//+------------------------------------------------------------------+
CNode::CNode()
  {
   m_count++;
   uchar name[NAME_SIZE];
   for(int i=0; i<NAME_SIZE; i++)
      name[i]=(uchar)(65+(MathRand()%25));
   m_name=CharArrayToString(name);
   m_name= CharArrayToString(name)+"_"+((string)m_count);
   m_showed = false;
   FontColor=C'255,255,0';
  }
//+------------------------------------------------------------------+
//| Add new element                                                  |
//+------------------------------------------------------------------+
void CNode::AddElement(CNode *element)
{
   m_elements.Add(element);
}
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNode::Show(void)
  {
   OnShow();
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
CNode::~CNode(void)
  {
   Hide();
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNode::OnShow(void)
  {
   for(int i=0; i<m_elements.Total(); i++)
     {
      CNode *node=m_elements.At(i);
      node.Show();
     }
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNode::Hide(void)
  {
   OnHide();
   ObjectDelete(ChartID(),m_name);
   m_showed = false;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNode::OnHide(void)
  {
   for(int i=0; i<m_elements.Total(); i++)
     {
      CObject *obj=m_elements.At(i);
      //if(CheckPointer(obj)!=POINTER_DYNAMIC)continue;
      if(CheckPointer(obj)==POINTER_INVALID)continue;
      CNode *node=m_elements.At(i);
      node.Hide();
     }

  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNode::Event(CEvent *event)
  {
   if(event.EventType()==EVENT_CHART_OBJECT_CLICK)
     {
      CEventChartObjClick *eventObjClick=event;
      if(eventObjClick.ObjectName()==m_name)
         OnClick(/*eventObjClick*/);
     }
   else if(event.EventType()==EVENT_FREFRESH)
     {
      CEventRefresh *refresh=event;
      OnRefresh(refresh);
     }
   for(int i=0; i<m_elements.Total(); i++)
     {
      CNode *node=m_elements.At(i);
      node.Event(event);
     }
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
bool CNode::IsShowed(void)
  {
   //return m_showed;
   if(ObjectFind(ChartID(),m_name)<0)
      return false;
   return true;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNode::OnClick(void)
  {
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNode::OnRefresh(CEventRefresh *event)
  {
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CNode::Tooltip(string message)
  {
   ObjectSetString(ChartID(),m_name,OBJPROP_TOOLTIP,message);
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
string CNode::Name(void)
  {
   return m_name;
  }
//+------------------------------------------------------------------+
