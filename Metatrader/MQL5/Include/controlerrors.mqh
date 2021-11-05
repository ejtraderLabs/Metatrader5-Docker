//+------------------------------------------------------------------+
//|                                                ControlErrors.mqh |
//|                                             Copyright KlimMalgin |
//| The library should be located in directory:                      |
//| MetaTrader 5/MQL5/Include/                                       |
//| https://www.mql5.com/en/articles/70                              |
//+------------------------------------------------------------------+
#property copyright "KlimMalgin"
#property link      ""



class ControlErrors
{
private:

   // Flags that define what types of reports need to be enabled
   bool _PlaySound;    // Play or don't play a sound when an error occurs.
   bool _PrintInfo;    // Add error details to the journal of Expert Advisors
   bool _AlertInfo;    // Generate Alert with error details
   bool _WriteFile;    // Record reports on errors into a file or not
   
   // A structure for storing error data elements that use this structure
   struct Code
   {
      int code;      // Error code
      string desc;   // Description of the error code
   };
   Code Errors[];    // Array that contains error codes and their descriptions
   Code _UserError;  // Stores information about a custome error
   Code _Error;      // Stores information about the last error of any type
   
   // Different service properties
   short  _CountErrors;     // Number of errors stored in array Errors[]
   string _PlaySoundFile;   // File that will be played for an alert sound
   string _DataPath;        // Path to the log storing directory

   
public:
   // Constructor
   ControlErrors(void);
   
   // Methods for setting flags
   void SetSound(bool value);          // Play or don't play a sound when an error occurs
   void SetPrint(bool value);          // Enter error data the the journal of Expert Advisors or not
   void SetAlert(bool value);          // Generate an Alert message or not
   void SetWriteFlag(bool flag);       // Set the writing flag. true - keep logs, false - do not keep
   
   // Methods for working with errors
   int  mGetLastError();            // Returns contents of the system variable _LastError
   int  mGetError();                // Returns code of the last obtained error
   int  mGetTypeError();            // Returns error type (Custom = 1 ore predefined = 0)
   void mResetLastError();          // Resets the contents of the system variable _LastError
   void mSetUserError(ushort value, string desc = "");   // Sets the custom error
   void mResetUserError();          // Resets class fields that contain information about the custom error
   void mResetError();              // Resets the structure that contains information about the last error
   string mGetDesc(int nErr = 0);   // Returns error description by the number, or that of the current error of no number
   int Check(string st = "");       // Method to check the current system state for errors
   
   // Alert methods (Alert, Print, Sound)
   void mAlert(string message = "");
   void mPrint(string message = "");
   void mSound();
      
   // Various service methods
   void SetPlaySoundFile(string file); // Method sets the file name to play an sound
   void SetWritePath(string path);     // Set the path to store logs   
   int mFileWrite(string message = "");// Record into a file the available information about the last error
};

void ControlErrors::ControlErrors(void)
{
   SetAlert(false);
   SetPrint(false);
   SetSound(false);
   SetWriteFlag(false);
   SetPlaySoundFile("alert.wav");
   SetWritePath("LogErrors.txt"); 
   
   _CountErrors = 150;

   ArrayResize(Errors, _CountErrors);
   // Return codes of a trade server
   Errors[0].code = 10004;Errors[0].desc = "Requote";
   Errors[1].code = 10006;Errors[1].desc = "Request rejected";
   Errors[2].code = 10007;Errors[2].desc = "Request canceled by trader";
   Errors[3].code = 10008;Errors[3].desc = "Order placed";
   Errors[4].code = 10009;Errors[4].desc = "Request is completed";
   Errors[5].code = 10010;Errors[5].desc = "Request is partially completed";
   Errors[6].code = 10011;Errors[6].desc = "Request processing error";
   Errors[7].code = 10012;Errors[7].desc = "Request canceled by timeout";
   Errors[8].code = 10013;Errors[8].desc = "Invalid request";
   Errors[9].code = 10014;Errors[9].desc = "Invalid volume in the request";
   Errors[10].code = 10015;Errors[10].desc = "Invalid price in the request";
   Errors[11].code = 10016;Errors[11].desc = "Invalid stops in the request";
   Errors[12].code = 10017;Errors[12].desc = "Trade is disabled";
   Errors[13].code = 10018;Errors[13].desc = "Market is closed";
   Errors[14].code = 10019;Errors[14].desc = "There is not enough money to fulfill the request";
   Errors[15].code = 10020;Errors[15].desc = "Prices changed";
   Errors[16].code = 10021;Errors[16].desc = "There are no quotes to process the request";
   Errors[17].code = 10022;Errors[17].desc = "Invalid order expiration date in the request";
   Errors[18].code = 10023;Errors[18].desc = "Order state changed";
   Errors[19].code = 10024;Errors[19].desc = "Too frequent requests";
   Errors[20].code = 10025;Errors[20].desc = "No changes in request";
   Errors[21].code = 10026;Errors[21].desc = "Autotrading disabled by server";
   Errors[22].code = 10027;Errors[22].desc = "Autotrading disabled by client terminal";
   Errors[23].code = 10028;Errors[23].desc = "Request locked for processing";
   Errors[24].code = 10029;Errors[24].desc = "Order or position frozen";
   Errors[25].code = 10030;Errors[25].desc = "Invalid order filling type";

   // Common Errors
   Errors[26].code = 4001;Errors[26].desc = "Unexpected internal error";
   Errors[27].code = 4002;Errors[27].desc = "Wrong parameter in the inner call of the client terminal function";
   Errors[28].code = 4003;Errors[28].desc = "Wrong parameter when calling the system function";
   Errors[29].code = 4004;Errors[29].desc = "Not enough memory to perform the system function";
   Errors[30].code = 4005;Errors[30].desc = "The structure contains objects of strings and/or dynamic arrays and/or structure of such objects and/or classes";
   Errors[31].code = 4006;Errors[31].desc = "Array of a wrong type, wrong size, or a damaged object of a dynamic array";
   Errors[32].code = 4007;Errors[32].desc = "Not enough memory for the relocation of an array, or an attempt to change the size of a static array";
   Errors[33].code = 4008;Errors[33].desc = "Not enough memory for the relocation of string";
   Errors[34].code = 4009;Errors[34].desc = "Not initialized string";
   Errors[35].code = 4010;Errors[35].desc = "Invalid date and/or time";
   Errors[36].code = 4011;Errors[36].desc = "Requested array size exceeds 2 GB";
   Errors[37].code = 4012;Errors[37].desc = "Wrong pointer";
   Errors[38].code = 4013;Errors[38].desc = "Wrong type of pointer";
   Errors[39].code = 4014;Errors[39].desc = "System function is not allowed to call";

   // Charts
   Errors[40].code = 4101;Errors[40].desc = "Wrong chart ID";
   Errors[41].code = 4102;Errors[41].desc = "Chart does not respond";
   Errors[42].code = 4103;Errors[42].desc = "Chart not found";
   Errors[43].code = 4104;Errors[43].desc = "No Expert Advisor in the chart that could handle the event";
   Errors[44].code = 4105;Errors[44].desc = "Chart opening error";
   Errors[45].code = 4106;Errors[45].desc = "Failed to change chart symbol and period";
   Errors[46].code = 4107;Errors[46].desc = "Wrong parameter for timer";
   Errors[47].code = 4108;Errors[47].desc = "Failed to create timer";
   Errors[48].code = 4109;Errors[48].desc = "Wrong chart property ID";
   Errors[49].code = 4110;Errors[49].desc = "Error creating screenshots";
   Errors[50].code = 4111;Errors[50].desc = "Error navigating through chart";
   Errors[51].code = 4112;Errors[51].desc = "Error applying template";
   Errors[52].code = 4113;Errors[52].desc = "Subwindow containing the indicator was not found";

   // Graphical Objects
   Errors[53].code = 4201;Errors[53].desc = "Error working with a graphical object";
   Errors[54].code = 4202;Errors[54].desc = "Graphical object was not found";
   Errors[55].code = 4203;Errors[55].desc = "Wrong ID of a graphical object property";
   Errors[56].code = 4204;Errors[56].desc = "Unable to get date corresponding to the value";
   Errors[57].code = 4205;Errors[57].desc = "Unable to get value corresponding to the date";
   
   // MarketInfo   
   Errors[58].code = 4301;Errors[58].desc = "Unknown symbol";
   Errors[59].code = 4302;Errors[59].desc = "Symbol is not selected in MarketWatch";
   Errors[60].code = 4303;Errors[60].desc = "Wrong identifier of a symbol property";
   Errors[61].code = 4304;Errors[61].desc = "Time of the last tick is not known (no ticks)";
   
   // History Access
   Errors[62].code = 4401;Errors[62].desc = "Requested history not found";
   Errors[63].code = 4402;Errors[63].desc = "Wrong ID of the history property";
   
   // Global_Variables
   Errors[64].code = 4501;Errors[64].desc = "Global variable of the client terminal is not found";
   Errors[65].code = 4502;Errors[65].desc = "Global variable of the client terminal with the same name already exists";
   Errors[66].code = 4510;Errors[66].desc = "Email sending failed";
   Errors[67].code = 4511;Errors[67].desc = "Sound playing failed";
   Errors[68].code = 4512;Errors[68].desc = "Wrong identifier of the program property";
   Errors[69].code = 4513;Errors[69].desc = "Wrong identifier of the terminal property";
   Errors[70].code = 4514;Errors[70].desc = "File sending via ftp failed";
   
   // Custom Indicator Buffers
   Errors[71].code = 4601;Errors[71].desc = "Not enough memory for the distribution of indicator buffers";
   Errors[72].code = 4602;Errors[72].desc = "Wrong indicator buffer index";
   
   // Custom Indicator Properties
   Errors[73].code = 4603;Errors[73].desc = "Wrong ID of the custom indicator property";
   
   // Account
   Errors[74].code = 4701;Errors[74].desc = "Wrong account property ID";
   Errors[75].code = 4751;Errors[75].desc = "Wrong trade property ID;";
   Errors[76].code = 4752;Errors[76].desc = "Trading by Expert Advisors prohibited";
   Errors[77].code = 4753;Errors[77].desc = "Position not found";
   Errors[78].code = 4754;Errors[78].desc = "Order not found";
   Errors[79].code = 4755;Errors[79].desc = "Deal not found";
   Errors[80].code = 4756;Errors[80].desc = "Trade request sending failed";
   Errors[81].code = 4757;Errors[81].desc = "Timeout exceeded when selecting (searching) specified data";
   
   // Indicators
   Errors[82].code = 4801;Errors[82].desc = "Unknown symbol";
   Errors[83].code = 4802;Errors[83].desc = "Indicator cannot be created";
   Errors[84].code = 4803;Errors[84].desc = "Not enough memory to add the indicator";
   Errors[85].code = 4804;Errors[85].desc = "The indicator cannot be applied to another indicator";
   Errors[86].code = 4805;Errors[86].desc = "Error applying an indicator to chart";
   Errors[87].code = 4806;Errors[87].desc = "Requested data not found";
   Errors[88].code = 4807;Errors[88].desc = "Wrong index of the requested indicator buffer";
   Errors[89].code = 4808;Errors[89].desc = "Wrong number of parameters when creating an indicator";
   Errors[90].code = 4809;Errors[90].desc = "No parameters when creating an indicator";
   Errors[91].code = 4810;Errors[91].desc = "The first parameter in the array must be the name of the custom indicator";
   Errors[92].code = 4811;Errors[92].desc = "Invalid parameter type in the array when creating an indicator";
   
   // Depth of Market
   Errors[93].code = 4901;Errors[93].desc = "Depth Of Market can not be added";
   Errors[94].code = 4902;Errors[94].desc = "Depth Of Market can not be removed";
   Errors[95].code = 4903;Errors[95].desc = "The data from Depth Of Market can not be obtained";
   Errors[96].code = 4904;Errors[96].desc = "Error in subscribing to receive new data from Depth Of Market";
   
   // File Operations
   Errors[97].code = 5001;Errors[97].desc = "More than 64 files cannot be opened at the same time";
   Errors[98].code = 5002;Errors[98].desc = "Invalid file name";
   Errors[99].code = 5003;Errors[99].desc = "Too long file name";
   Errors[100].code = 5004;Errors[100].desc = "File opening error";
   Errors[101].code = 5005;Errors[101].desc = "Not enough memory for cache to read";
   Errors[102].code = 5006;Errors[102].desc = "File deleting error";
   Errors[103].code = 5007;Errors[103].desc = "A file with this handle was closed, or was not opened at all";
   Errors[104].code = 5008;Errors[104].desc = "Wrong file handle";
   Errors[105].code = 5009;Errors[105].desc = "The file must be opened for writing";
   Errors[106].code = 5010;Errors[106].desc = "The file must be opened for reading";
   Errors[107].code = 5011;Errors[107].desc = "The file must be opened as a binary one";
   Errors[108].code = 5012;Errors[108].desc = "The file must be opened as a text";
   Errors[109].code = 5013;Errors[109].desc = "The file must be opened as a text or CSV";
   Errors[110].code = 5014;Errors[110].desc = "The file must be opened as CSV";
   Errors[111].code = 5015;Errors[111].desc = "File reading error";
   Errors[112].code = 5016;Errors[112].desc = "String size must be specified, because the file is opened as binary";
   Errors[113].code = 5017;Errors[113].desc = "A text file must be for string arrays, for other arrays - binary";
   Errors[114].code = 5018;Errors[114].desc = "This is not a file, this is a directory";
   Errors[115].code = 5019;Errors[115].desc = "File does not exist";
   Errors[116].code = 5020;Errors[116].desc = "File can not be rewritten";
   Errors[117].code = 5021;Errors[117].desc = "Wrong directory name";
   Errors[118].code = 5022;Errors[118].desc = "Directory does not exist";
   Errors[119].code = 5023;Errors[119].desc = "This is a file, not a directory";
   Errors[120].code = 5024;Errors[120].desc = "The directory cannot be removed";
   
   // String Casting
   Errors[121].code = 5030;Errors[121].desc = "No date in the string";
   Errors[122].code = 5031;Errors[122].desc = "Wrong date in the string";
   Errors[123].code = 5032;Errors[123].desc = "Wrong time in the string";
   Errors[124].code = 5033;Errors[124].desc = "Error converting string to date";
   Errors[125].code = 5034;Errors[125].desc = "Not enough memory for the string";
   Errors[126].code = 5035;Errors[126].desc = "The string length is less than expected";
   Errors[127].code = 5036;Errors[127].desc = "Too large number, more than ULONG_MAX";
   Errors[128].code = 5037;Errors[128].desc = "Invalid format string";
   Errors[129].code = 5038;Errors[129].desc = "Amount of format specifiers more than the parameters";
   Errors[130].code = 5039;Errors[130].desc = "Amount of parameters more than the format specifiers";
   Errors[131].code = 5040;Errors[131].desc = "Damaged parameter of string type";
   Errors[132].code = 5041;Errors[132].desc = "Position outside the string";
   Errors[133].code = 5042;Errors[133].desc = "0 added to the string end, a useless operation";
   Errors[134].code = 5043;Errors[134].desc = "Unknown data type when converting to a string";
   Errors[135].code = 5044;Errors[135].desc = "Damaged string object";
   
   // Operations with Arrays
   Errors[136].code = 5050;Errors[136].desc = "Copying incompatible arrays. String array can be copied only to a string array, and a numeric array - in numeric array only";
   Errors[137].code = 5051;Errors[137].desc = "The receiving array is declared as AS_SERIES, and it is of insufficient size";
   Errors[138].code = 5052;Errors[138].desc = "Too small array, the starting position is outside the array";
   Errors[139].code = 5053;Errors[139].desc = "An array of zero length";
   Errors[140].code = 5054;Errors[140].desc = "Must be a numeric array";
   Errors[141].code = 5055;Errors[141].desc = "Must be a one-dimensional array";
   Errors[142].code = 5056;Errors[142].desc = "Timeseries cannot be used";
   Errors[143].code = 5057;Errors[143].desc = "Must be an array of type double";
   Errors[144].code = 5058;Errors[144].desc = "Must be an array of type float";
   Errors[145].code = 5059;Errors[145].desc = "Must be an array of type long";
   Errors[146].code = 5060;Errors[146].desc = "Must be an array of type int";
   Errors[147].code = 5061;Errors[147].desc = "Must be an array of type short";
   Errors[148].code = 5062;Errors[148].desc = "Must be an array of type char";
}

void ControlErrors::SetAlert(bool value)
{
   _AlertInfo = value;
}

void ControlErrors::SetPrint(bool value)
{
   _PrintInfo = value;
}

void ControlErrors::SetSound(bool value)
{
   _PlaySound = value;
}

void ControlErrors::SetWriteFlag(bool flag)
{
   _WriteFile = flag;
}

void ControlErrors::SetWritePath(string path)
{
   _DataPath = path;
}

void ControlErrors::SetPlaySoundFile(string file)
{
   _PlaySoundFile = file;
}

int ControlErrors::mGetLastError(void)
{
   _Error.code = GetLastError();
   _Error.desc = mGetDesc(_Error.code);
   return _Error.code;
}

void ControlErrors::mResetLastError(void)
{
   ResetLastError();
}

int ControlErrors::mGetError(void)
{
   return _Error.code;
}

void ControlErrors::mResetError(void)
{
   _Error.code = 0;
   _Error.desc = "";
}

int ControlErrors::mGetTypeError(void)
{
   if (mGetError() < ERR_USER_ERROR_FIRST)
   {
      return 0;
   }
   else if (mGetError() >= ERR_USER_ERROR_FIRST)
   {
      return 1;
   }
   return -1;
}

void ControlErrors::mSetUserError(ushort value, string desc = "")
{
   SetUserError(value);
   _UserError.code = value;
   _UserError.desc = desc;
}

void ControlErrors::mResetUserError(void)
{
   _UserError.code = 0;
   _UserError.desc = "";
}

string ControlErrors::mGetDesc(int nErr=0)
{
   int ErrorNumber = 0;
   string ReturnDesc = "";
   
   ErrorNumber = (mGetError()>0)?mGetError():ErrorNumber;
   ErrorNumber = (nErr>0)?nErr:ErrorNumber;
   
   if ((ErrorNumber > 0) && (ErrorNumber < ERR_USER_ERROR_FIRST))
   {
      for (int i = 0;i<_CountErrors;i++)
      {
         if (Errors[i].code == ErrorNumber)
         {
            ReturnDesc = Errors[i].desc;
            break;
         }
      }
   }
   else if (ErrorNumber > ERR_USER_ERROR_FIRST)
   {
      ReturnDesc = (_UserError.desc=="")?"Custom error":_UserError.desc;
   }
      
   if (ReturnDesc == ""){return "Unknown error code: "+(string)ErrorNumber;}
   return ReturnDesc;
}

void ControlErrors::mAlert(string message="")
{
   if (_AlertInfo == true)
   {
      if (message == "")
      {
         if (mGetError() > 0)
         {
            Alert("Error ",mGetError()," - ",mGetDesc());
         }
      }
      else
      {
         Alert(message);
      }   
   }
}

void ControlErrors::mPrint(string message="")
{
   if (_PrintInfo == true)
   {
      if (message == "")
      {
         if (mGetError() > 0)
         {
            Print("Error ",mGetError()," - ",mGetDesc());
         }
      }
      else
      {
         Print(message);
      }
   }
}

void ControlErrors::mSound(void)
{
   if (_PlaySound == true)
   {
      PlaySound(_PlaySoundFile);
   }
}

int ControlErrors::Check(string st="")
{
   int errNum = 0;
   errNum = mGetLastError();
   mFileWrite();
   mAlert(st);
   mPrint(st);
   mSound();
   mResetError();
   mResetLastError();
   mResetUserError();
   return errNum;
}

int ControlErrors::mFileWrite(string message = "")
{
   int      handle  = 0,
            _return = 0;
   datetime time    = TimeCurrent();
   string   text    = (message != "")?message:time+" - Error "+mGetError()+" "+mGetDesc();
   
   if (_WriteFile == true)
   {
      handle = FileOpen(_DataPath,FILE_READ|FILE_WRITE|FILE_TXT|FILE_ANSI);
      if (handle != INVALID_HANDLE)
      {
         ulong size = FileSize(handle);
         FileSeek(handle,size,SEEK_SET);
         _return = FileWrite(handle,text);
         FileClose(handle);
      }
   }
   return _return;
}


