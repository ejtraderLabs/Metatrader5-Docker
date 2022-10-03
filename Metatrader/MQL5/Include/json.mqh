//+------------------------------------------------------------------ß
//|                                                            JAson |
//|    This software is licensed under the MIT https://goo.gl/eyJgHe |
//+------------------------------------------------------------------+
#property copyright "Copyright © 2006-2017"
#property version "1.12"
#property strict
//------------------------------------------------------------------	enum enJAType
enum enJAType
{
  jtUNDEF, jtNULL, jtBOOL, jtINT, jtDBL, jtSTR, jtARRAY, jtOBJ
};

//------------------------------------------------------------------	class CJAVal
class CJAVal
{
  public:
  virtual void Clear (enJAType jt = jtUNDEF, bool savekey = false)
  {
    m_parent = NULL;

    if (!savekey)
      groupName = "";

    m_type = jt;
    m_bv   = false;
    m_iv   = 0;
    m_dv   = 0;
    m_prec = 5;
    m_sv   = "";
    ArrayResize (group, 0, 100);
  }
  virtual bool Copy (const CJAVal &a)
  {
    groupName = a.groupName;
    CopyData (a);
    return true;
  }
  virtual void CopyData (const CJAVal &a)
  {
    m_type = a.m_type;
    m_bv   = a.m_bv;
    m_iv   = a.m_iv;
    m_dv   = a.m_dv;
    m_prec = a.m_prec;
    m_sv   = a.m_sv;
    CopyArr (a);
  }
  virtual void CopyArr (const CJAVal &a)
  {
    int n = ArrayResize (group, ArraySize (a.group));

    for (int i = 0; i < n; i++)
    {
      group [i] = a.group [i];
      group [i].m_parent = GetPointer (this);
    }
  }

  public:
  CJAVal            group [];   // group
  string            groupName;  // group name
  string            m_lkey;
  CJAVal           *m_parent;
  enJAType          m_type;
  bool              m_bv;
  long              m_iv;
  double            m_dv;
  int               m_prec;
  string            m_sv;
  static int        code_page;

  public:
  CJAVal ()
  {
    Clear ();
  }
  CJAVal (CJAVal *aparent, enJAType atype)
  {
    Clear ();
    m_type = atype;
    m_parent = aparent;
  }
  CJAVal (enJAType t, string a)
  {
    Clear ();
    FromStr (t, a);
  }
  CJAVal (const int a)
  {
    Clear ();
    m_type = jtINT;
    m_iv   = a;
    m_dv   = (double)m_iv;
    m_sv   = IntegerToString (m_iv);
    m_bv   = m_iv != 0;
  }
  CJAVal (const long a)
  {
    Clear ();
    m_type = jtINT;
    m_iv   = a;
    m_dv   = (double)m_iv;
    m_sv   = IntegerToString (m_iv);
    m_bv   = m_iv != 0;
  }
  CJAVal (const double a, int aprec = -100)
  {
    Clear ();
    m_type = jtDBL;
    m_dv = a;

    if (aprec > -100)
      m_prec = aprec;

    m_iv = (long)m_dv;
    m_sv = DoubleToString (m_dv, m_prec);
    m_bv = m_iv != 0;
  }
  CJAVal (const bool a)
  {
    Clear ();
    m_type = jtBOOL;
    m_bv   = a;
    m_iv   = m_bv;
    m_dv   = m_bv;
    m_sv   = IntegerToString (m_iv);
  }
  CJAVal (const CJAVal &a)
  {
    Clear ();
    Copy (a);
  }
  ~CJAVal ()
  {
    Clear ();
  }

  public:
  int Size ()
  {
    return ArraySize (group);
  }
  virtual bool IsNumeric ()
  {
    return m_type == jtDBL || m_type == jtINT;
  }
  virtual CJAVal *FindKey (string akey)
  {
    for (int i = Size () - 1; i >= 0; --i)
      if (group [i].groupName == akey)
        return GetPointer (group [i]);
    return NULL;
  }
  virtual CJAVal *HasKey (string akey, enJAType atype = jtUNDEF)
  {
    CJAVal *e = FindKey (akey);
    if (CheckPointer (e) != POINTER_INVALID)
    {
      if (atype == jtUNDEF || atype == e.m_type)
        return GetPointer (e);
    }
    return NULL;
  }
  virtual CJAVal *operator [] (string akey);
  virtual CJAVal *operator [] (int i);
  void operator= (const CJAVal &a)
  {
    Copy (a);
  }
  void operator= (const int a)
  {
    m_type = jtINT;
    m_iv   = a;
    m_dv   = (double)m_iv;
    m_bv   = m_iv != 0;
  }
  void operator= (const long a)
  {
    m_type = jtINT;
    m_iv   = a;
    m_dv   = (double)m_iv;
    m_bv   = m_iv != 0;
  }
  void operator= (const double a)
  {
    m_type = jtDBL;
    m_dv   = a;
    m_iv   = (long)m_dv;
    m_bv   = m_iv != 0;
  }
  void operator= (const bool a)
  {
    m_type = jtBOOL;
    m_bv   = a;
    m_iv   = (long)m_bv;
    m_dv   = (double)m_bv;
  }
  void operator= (string a)
  {
    m_type = (a != NULL) ? jtSTR : jtNULL;
    m_sv   = a;
    int err = GetLastError();
    m_iv   = StringToInteger (m_sv);
    m_dv   = StringToDouble (m_sv);
    if ((err == 0) && (GetLastError() != 0)) { //Ignore conversion errors if they occur
   ResetLastError();
}
    m_bv   = a != NULL;
  }

  bool operator== (const int a)
  {
    return m_iv == a;
  }
  bool operator== (const long a)
  {
    return m_iv == a;
  }
  bool operator== (const double a)
  {
    return m_dv == a;
  }
  bool operator== (const bool a)
  {
    return m_bv == a;
  }
  bool operator== (string a)
  {
    return m_sv == a;
  }

  bool operator!= (const int a)
  {
    return m_iv != a;
  }
  bool operator!= (const long a)
  {
    return m_iv != a;
  }
  bool operator!= (const double a)
  {
    return m_dv != a;
  }
  bool operator!= (const bool a)
  {
    return m_bv != a;
  }
  bool operator!= (string a)
  {
    return m_sv != a;
  }

  long ToInt () const
  {
    return m_iv;
  }
  double ToDbl () const
  {
    return m_dv;
  }
  bool ToBool () const
  {
    return m_bv;
  }
  string ToStr ()
  {
    return m_sv;
  }

  virtual void FromStr (enJAType t, string a)
  {
    m_type = t;
    switch (m_type)
    {
      case jtBOOL:
        m_bv = (StringToInteger (a) != 0);
        m_iv = (long)m_bv;
        m_dv = (double)m_bv;
        m_sv = a;
        break;
      case jtINT:
        m_iv = StringToInteger (a);
        m_dv = (double)m_iv;
        m_sv = a;
        m_bv = m_iv != 0;
        break;
      case jtDBL:
        m_dv = StringToDouble (a);
        m_iv = (long)m_dv;
        m_sv = a;
        m_bv = m_iv != 0;
        break;
      case jtSTR:
        m_sv   = Unescape (a);
        m_type = (m_sv != NULL) ? jtSTR : jtNULL;
        m_iv   = StringToInteger (m_sv);
        m_dv   = StringToDouble (m_sv);
        m_bv   = m_sv != NULL;
        break;
    }
  }
  virtual string GetStr (char &js [], int i, int slen)
  {
    if (slen == 0)
      return "";
    char cc [];
    ArrayCopy (cc, js, 0, i, slen);
    return CharArrayToString (cc, 0, WHOLE_ARRAY, CJAVal::code_page);
  }

  virtual void Set (const CJAVal &a)
  {
    if (m_type == jtUNDEF)
      m_type = jtOBJ;
    CopyData (a);
  }
  virtual void      Set (const CJAVal &list []);
  virtual CJAVal *Add (const CJAVal &item)
  {
    if (m_type == jtUNDEF)
      m_type = jtARRAY;
    /*ASSERT(m_type==jtOBJ || m_type==jtARRAY);*/ return AddBase (item);
  } // добавление
  virtual CJAVal *Add (const int a)
  {
    CJAVal item (a);
    return Add (item);
  }
  virtual CJAVal *Add (const long a)
  {
    CJAVal item (a);
    return Add (item);
  }
  virtual CJAVal *Add (const double a, int aprec = -2)
  {
    CJAVal item (a, aprec);
    return Add (item);
  }
  virtual CJAVal *Add (const bool a)
  {
    CJAVal item (a);
    return Add (item);
  }
  virtual CJAVal *Add (string a)
  {
    CJAVal item (jtSTR, a);
    return Add (item);
  }
  virtual CJAVal *AddBase (const CJAVal &item) // добавление
  {
    int c = Size ();
    ArrayResize (group, c + 1, 100);
    group [c] = item;
    group [c].m_parent = GetPointer (this);
    return GetPointer (group [c]);
  }
  virtual CJAVal *New ()
  {
    if (m_type == jtUNDEF)
      m_type = jtARRAY;
    /*ASSERT(m_type==jtOBJ || m_type==jtARRAY);*/ return NewBase ();
  } // добавление
  virtual CJAVal *NewBase () // добавление
  {
    int c = Size ();
    ArrayResize (group, c + 1, 100);
    return GetPointer (group [c]);
  }

  virtual string    Escape (string a);
  virtual string    Unescape (string a);
  public:
  virtual void      Serialize (string &js, bool bf = false, bool bcoma = false);
  virtual string Serialize ()
  {
    string js;
    Serialize (js);
    return js;
  }
  virtual bool      Deserialize (char &js [], int slen, int &i);
  virtual bool      ExtrStr (char &js [], int slen, int &i);
  virtual bool      Deserialize (string js, int acp = CP_ACP)
  {
    int i = 0;
    Clear ();
    CJAVal::code_page = acp;
    char arr [];
    int slen = StringToCharArray (js, arr, 0, WHOLE_ARRAY, CJAVal::code_page);
    return Deserialize (arr, slen, i);
  }
  virtual bool Deserialize (char &js [], int acp = CP_ACP)
  {
    int i = 0;
    Clear ();
    CJAVal::code_page = acp;
    return Deserialize (js, ArraySize (js), i);
  }
};

int CJAVal::code_page = CP_ACP;

//------------------------------------------------------------------	operator[]
CJAVal *CJAVal::operator[] (string akey)
{
  if (m_type == jtUNDEF)
    m_type = jtOBJ;
  CJAVal *v = FindKey (akey);
  if (v)
    return v;
  CJAVal b (GetPointer (this), jtUNDEF);
  b.groupName = akey;
  v = Add (b);
  return v;
}
//------------------------------------------------------------------	operator[]
CJAVal *CJAVal::operator[] (int i)
{
  if (m_type == jtUNDEF)
    m_type = jtARRAY;
  while (i >= Size ())
  {
    CJAVal b (GetPointer (this), jtUNDEF);
    if (CheckPointer (Add (b)) == POINTER_INVALID)
      return NULL;
  }
  return GetPointer (group [i]);
}
//------------------------------------------------------------------	Set
void CJAVal::Set (const CJAVal &list [])
{
  if (m_type == jtUNDEF)
    m_type = jtARRAY;
  int n = ArrayResize (group, ArraySize (list), 100);
  for (int i = 0; i < n; ++i)
  {
    group [i] = list [i];
    group [i].m_parent = GetPointer (this);
  }
}
//------------------------------------------------------------------	Serialize
void CJAVal::Serialize (string &js, bool bkey/*=false*/, bool coma/*=false*/)
{
  if (m_type == jtUNDEF)
    return;
  if (coma)
    js += ",";
  if (bkey)
    js += StringFormat ("\"%s\":", groupName);
  int _n = Size ();
  switch (m_type)
  {
    case jtNULL:
      js += "null";
      break;
    case jtBOOL:
      js += (m_bv ? "true" : "false");
      break;
    case jtINT:
      js += IntegerToString (m_iv);
      break;
    case jtDBL:
      js += DoubleToString (m_dv, 5);
      break;
    case jtSTR:
    {
      string ss = Escape (m_sv);
      if (StringLen (ss) > 0)
        js += StringFormat ("\"%s\"", ss);
      else
        js += "null";
    }
      break;
    case jtARRAY:
      js += "[";
      for (int i = 0; i < _n; i++)
        group [i].Serialize (js, false, i > 0);
      js += "]";
      break;
    case jtOBJ:
      js += "{";
      for (int i = 0; i < _n; i++)
        group [i].Serialize (js, true, i > 0);
      js += "}";
      break;
  }
}

//------------------------------------------------------------------	Deserialize
bool CJAVal::Deserialize (char &js [], int slen, int &i)
{
  string num = "0123456789+-.eE";
  int i0 = i;
  for (; i < slen; i++)
  {
    char c = js [i];
    if (c == 0)
      break;
    switch (c)
    {
      case '\t':
      case '\r':
      case '\n':
      case ' ': // пропускаем из имени пробелы
        i0 = i + 1;
        break;

      case '[': // начало массива. создаём объекты и забираем из js
      {
        i0 = i + 1;
        if (m_type != jtUNDEF) // если значение уже имеет тип, то это ошибка
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        m_type = jtARRAY; // задали тип значения
        i++;
        CJAVal val (GetPointer (this), jtUNDEF);
        while (val.Deserialize (js, slen, i))
        {
          if (val.m_type != jtUNDEF)
            Add (val);
          if (val.m_type == jtINT || val.m_type == jtDBL || val.m_type == jtARRAY)
            i++;
          val.Clear ();
          val.m_parent = GetPointer (this);
          if (js [i] == ']')
            break;
          i++;
          if (i >= slen)
          {
            Print (groupName + " " + string (__LINE__));
            return false;
          }
        }
        if (i >= slen)
          return false;
        return js [i] == ']' || js [i] == 0;
      }
        break;
      case ']':
        if (!m_parent)
          return false;
        return m_parent.m_type == jtARRAY; // конец массива, текущее значение должны быть массивом

      case ':':
      {
        if (m_lkey == "")
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        CJAVal val (GetPointer (this), jtUNDEF);
        CJAVal *oc = Add (val); // тип объекта пока не определён
        oc.groupName = m_lkey;
        m_lkey = ""; // задали имя ключа
        i++;
        if (!oc.Deserialize (js, slen, i))
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        break;
      }
      case ',': // разделитель значений // тип значения уже должен быть определён
        i0 = i + 1;
        if (!m_parent && m_type != jtOBJ)
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        else
          if (m_parent)
          {
            if (m_parent.m_type != jtARRAY && m_parent.m_type != jtOBJ)
            {
              Print (groupName + " " + string (__LINE__));
              return false;
            }
            if (m_parent.m_type == jtARRAY && m_type == jtUNDEF)
              return true;
          }
        break;

        // примитивы могут быть ТОЛЬКО в массиве / либо самостоятельно
      case '{': // начало объекта. создаем объект и забираем его из js
        i0 = i + 1;
        if (m_type != jtUNDEF) // ошибка типа
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        m_type = jtOBJ; // задали тип значения
        i++;
        if (!Deserialize (js, slen, i)) // вытягиваем его
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        if (i >= slen)
          return false;
        return js [i] == '}' || js [i] == 0;
        break;
      case '}':
        return m_type == jtOBJ; // конец объекта, текущее значение должно быть объектом

      case 't':
      case 'T': // начало true
      case 'f':
      case 'F': // начало false
        if (m_type != jtUNDEF) // ошибка типа
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        m_type = jtBOOL; // задали тип значения
        if (i + 3 < slen)
        {
          if (StringCompare (GetStr (js, i, 4), "true", false) == 0)
          {
            m_bv = true;
            i += 3;
            return true;
          }
        }
        if (i + 4 < slen)
        {
          if (StringCompare (GetStr (js, i, 5), "false", false) == 0)
          {
            m_bv = false;
            i += 4;
            return true;
          }
        }
        Print (groupName + " " + string (__LINE__));
        return false; // не тот тип или конец строки
        break;
      case 'n':
      case 'N': // начало null
        if (m_type != jtUNDEF) // ошибка типа
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        m_type = jtNULL; // задали тип значения
        if (i + 3 < slen)
          if (StringCompare (GetStr (js, i, 4), "null", false) == 0)
          {
            i += 3;
            return true;
          }
        Print (groupName + " " + string (__LINE__));
        return false; // не NULL или конец строки
        break;

      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
      case '-':
      case '+':
      case '.': // начало числа
      {
        if (m_type != jtUNDEF) // ошибка типа
        {
          Print (groupName + " " + string (__LINE__));
          return false;
        }
        bool dbl = false; // задали тип значения
        int is = i;
        while (js [i] != 0 && i < slen)
        {
          i++;
          if (StringFind (num, GetStr (js, i, 1)) < 0)
            break;
          if (!dbl)
            dbl = (js [i] == '.' || js [i] == 'e' || js [i] == 'E');
        }
        m_sv = GetStr (js, is, i - is);
        if (dbl)
        {
          m_type = jtDBL;
          m_dv = StringToDouble (m_sv);
          m_iv = (long)m_dv;
          m_bv = m_iv != 0;
        }
        else // уточнии тип значения
        {
          m_type = jtINT;
          m_iv = StringToInteger (m_sv);
          m_dv = (double)m_iv;
          m_bv = m_iv != 0;
        }
        i--;
        return true; // отодвинулись на 1 символ назад и вышли
        break;
      }
      case '\"': // начало или конец строки
        if (m_type == jtOBJ) // если тип еще неопределён и ключ не задан
        {
          i++;
          int is = i;
          if (!ExtrStr (js, slen, i)) // это ключ, идём до конца строки
          {
            Print (groupName + " " + string (__LINE__));
            return false;
          }
          m_lkey = GetStr (js, is, i - is);
        }
        else
        {
          if (m_type != jtUNDEF) // ошибка типа
          {
            Print (groupName + " " + string (__LINE__));
            return false;
          }
          m_type = jtSTR; // задали тип значения
          i++;
          int is = i;
          if (!ExtrStr (js, slen, i))
          {
            Print (groupName + " " + string (__LINE__));
            return false;
          }
          FromStr (jtSTR, GetStr (js, is, i - is));
          return true;
        }
        break;
    }
  }
  return true;
}

//------------------------------------------------------------------	ExtrStr
bool CJAVal::ExtrStr (char &js [], int slen, int &i)
{
  for (; js [i] != 0 && i < slen; i++)
  {
    char c = js [i];
    if (c == '\"')
      break; // конец строки
    if (c == '\\' && i + 1 < slen)
    {
      i++;
      c = js [i];
      switch (c)
      {
        case '/':
        case '\\':
        case '\"':
        case 'b':
        case 'f':
        case 'r':
        case 'n':
        case 't':
          break; // это разрешенные
        case 'u': // \uXXXX
        {
          i++;
          for (int j = 0; j < 4 && i < slen && js [i] != 0; j++, i++)
          {
            if (!((js [i] >= '0' && js [i] <= '9') || (js [i] >= 'A' && js [i] <= 'F') || (js [i] >= 'a' && js [i] <= 'f'))) // не hex
            {
              Print (groupName + " " + CharToString (js [i]) + " " + string (__LINE__));
              return false;
            }
          }
          i--;
          break;
        }
        default:
          break; /*{ return false; } // неразрешенный символ с экранированием */
      }
    }
  }
  return true;
}
//------------------------------------------------------------------	Escape
string CJAVal::Escape (string a)
{
  ushort as [], s [];
  int n = StringToShortArray (a, as);
  if (ArrayResize (s, 2 * n) != 2 * n)
    return NULL;
  int j = 0;
  for (int i = 0; i < n; i++)
  {
    switch (as [i])
    {
      case '\\':
        s [j] = '\\';
        j++;
        s [j] = '\\';
        j++;
        break;
      case '"':
        s [j] = '\\';
        j++;
        s [j] = '"';
        j++;
        break;
      case '/':
        s [j] = '\\';
        j++;
        s [j] = '/';
        j++;
        break;
      case 8:
        s [j] = '\\';
        j++;
        s [j] = 'b';
        j++;
        break;
      case 12:
        s [j] = '\\';
        j++;
        s [j] = 'f';
        j++;
        break;
      case '\n':
        s [j] = '\\';
        j++;
        s [j] = 'n';
        j++;
        break;
      case '\r':
        s [j] = '\\';
        j++;
        s [j] = 'r';
        j++;
        break;
      case '\t':
        s [j] = '\\';
        j++;
        s [j] = 't';
        j++;
        break;
      default:
        s [j] = as [i];
        j++;
        break;
    }
  }
  a = ShortArrayToString (s, 0, j);
  return a;
}
//------------------------------------------------------------------	Unescape
string CJAVal::Unescape (string a)
{
  ushort as [], s [];
  int n = StringToShortArray (a, as);
  if (ArrayResize (s, n) != n)
    return NULL;
  int j = 0, i = 0;
  while (i < n)
  {
    ushort c = as [i];
    if (c == '\\' && i < n - 1)
    {
      switch (as [i + 1])
      {
        case '\\':
          c = '\\';
          i++;
          break;
        case '"':
          c = '"';
          i++;
          break;
        case '/':
          c = '/';
          i++;
          break;
        case 'b':
          c = 8;
          /*08='\b'*/;
          i++;
          break;
        case 'f':
          c = 12;
          /*0c=\f*/ i++;
          break;
        case 'n':
          c = '\n';
          i++;
          break;
        case 'r':
          c = '\r';
          i++;
          break;
        case 't':
          c = '\t';
          i++;
          break;
        case 'u': // \uXXXX
        {
          i += 2;
          ushort k = 0;
          for (int jj = 0; jj < 4 && i < n; jj++, i++)
          {
            c = as [i];
            ushort h = 0;
            if (c >= '0' && c <= '9')
              h = c - '0';
            else
              if (c >= 'A' && c <= 'F')
                h = c - 'A' + 10;
              else
                if (c >= 'a' && c <= 'f')
                  h = c - 'a' + 10;
                else
                  break; // не hex
            k += h * (ushort)pow (16, (3 - jj));
          }
          i--;
          c = k;
          break;
        }
      }
    }
    s [j] = c;
    j++;
    i++;
  }
  a = ShortArrayToString (s, 0, j);
  return a;
}
//+------------------------------------------------------------------+
