#property copyright     "Copyright © 2021 Artem Maltsev (Vivazzi)"
#property link          "https://vivazzi.pro"
#property version		"1.00"
#property description   "Requests is a simple HTTP library for mql4, built for human beings."
#property library

#define HINTERNET int
#define BOOL int
#define INTERNET_PORT int
#define DWORD int
#define DWORD_PTR int
#define LPDWORD int&
#define LPVOID uchar& 
#define LPCTSTR string&

string nill = "";

#import	"Kernel32.dll"
	DWORD GetLastError(int);
#import

#import "wininet.dll"
	DWORD InternetAttemptConnect(DWORD dwReserved);
	HINTERNET InternetOpenW(LPCTSTR lpszAgent, DWORD dwAccessType, LPCTSTR lpszProxyName, LPCTSTR lpszProxyBypass, DWORD dwFlags);
	HINTERNET InternetConnectW(HINTERNET hInternet, LPCTSTR lpszServerName, INTERNET_PORT nServerPort, LPCTSTR lpszUsername, LPCTSTR lpszPassword, DWORD dwService, DWORD dwFlags, DWORD_PTR dwContext);
	HINTERNET HttpOpenRequestW(HINTERNET hConnect, LPCTSTR lpszVerb, LPCTSTR lpszObjectName, LPCTSTR lpszVersion, LPCTSTR lpszReferer, int /*LPCTSTR* */ lplpszAcceptTypes, uint/*DWORD*/ dwFlags, DWORD_PTR dwContext);
	BOOL HttpSendRequestW(HINTERNET hRequest, LPCTSTR lpszHeaders, DWORD dwHeadersLength, LPVOID lpOptional[], DWORD dwOptionalLength);
	HINTERNET InternetOpenUrlW(HINTERNET hInternet, LPCTSTR lpszUrl, LPCTSTR lpszHeaders, DWORD dwHeadersLength, uint/*DWORD*/ dwFlags, DWORD_PTR dwContext);
	BOOL InternetReadFile(HINTERNET hFile, LPVOID lpBuffer[], DWORD dwNumberOfBytesToRead, LPDWORD lpdwNumberOfBytesRead);
	BOOL InternetCloseHandle(HINTERNET hInternet);
	BOOL InternetSetOptionW(HINTERNET hInternet, DWORD dwOption, LPDWORD lpBuffer, DWORD dwBufferLength);
	BOOL InternetQueryOptionW(HINTERNET hInternet, DWORD dwOption, LPDWORD lpBuffer, LPDWORD lpdwBufferLength);
#import

#define OPEN_TYPE_PRECONFIG               0
#define INTERNET_SERVICE_HTTP             3

#define INTERNET_FLAG_PRAGMA_NOCACHE      0x00000100
#define INTERNET_FLAG_KEEP_CONNECTION     0x00400000
#define INTERNET_FLAG_SECURE              0x00800000
#define INTERNET_FLAG_RELOAD              0x80000000
#define INTERNET_OPTION_SECURITY_FLAGS    31

#define ERROR_INTERNET_INVALID_CA         12045
#define SECURITY_FLAG_IGNORE_UNKNOWN_CA   0x00000100

#include "classes/request_data.mqh"
#include "classes/response.mqh"
#include "classes/_url_parts.mqh"


class Requests {
public:
	int port;
    string host;
	string path;
	string parameters;
	string user;
	string password;

	int h_session;
	int h_connect;

    Response response;

    Requests() {
        response = Response();
        h_session = -1; h_connect = -1; user = ""; password = "";
    }

    ~Requests() {
        close();
    }

    bool check_dll(string &error) {
        error = "";

        return(true);
    }

    bool is_same_host(string url) {
        _UrlParts _url_parts;
        _url_parts.split(url);

        if (_url_parts.host != host) return(false);

        return(true);
    }

    bool open(string url) {
        if (h_session > 0 || h_connect > 0) close();  // if connection was opend, close session and connect

        if (InternetAttemptConnect(0) != 0) { Print("requests: ERROR InternetAttemptConnect"); return(false); }

        string user_agent = "Microsoft Internet Explorer";
        int service = INTERNET_SERVICE_HTTP;

        h_session = InternetOpenW(user_agent, OPEN_TYPE_PRECONFIG, nill, nill, 0);
        if (h_session <= 0) { Print("requests: ERROR create session, InternetOpenW()"); close(); return(false); }

        h_connect = InternetConnectW(h_session, host, port, user, password, service, 0, 0);
	    if (h_connect <= 0) { Print("requests: ERROR create connect, InternetConnectW()"); close(); return(false); }

        return(true);
    }

    void close() {
        if (h_session>0) {InternetCloseHandle(h_session); h_session = -1;#ifdef DEBUG_REQUESTS Print("requests: close session");#endif}
	    if (h_connect>0) {InternetCloseHandle(h_connect); h_connect = -1;#ifdef DEBUG_REQUESTS Print("requests: close connect");#endif}
    }

    void update_url_parts(string url) {
        _UrlParts _url_parts;
        _url_parts.split(url);

        port = _url_parts.port;
        host = _url_parts.host;
        path = _url_parts.path;
        parameters = _url_parts.parameters;
    }

    void read(int h_request, string &out) {
        out = "";

        uchar ch[100];
        int dwBytes, h = -1;

        while(InternetReadFile(h_request, ch, 100, dwBytes)) {
            if (dwBytes <= 0) break; out += CharArrayToString(ch, 0, dwBytes);
        }
    }

    // --- GET ---
    Response get(string url, string& arr_data[][]) {
        return get(url, RequestData::to_str(arr_data));
    }

    Response get(string url, RequestData& _request_data) {
        return get(url, _request_data.to_str());
    }

    Response get(string url, string _str_data="") {
        string error;
        bool is_success;

        is_success = check_dll(error);
        if (!is_success) {
            response.error = error;
            return(response);
        };

        // append _str_data to end of url
        if (_str_data != "" && _str_data != NULL) {
            int index = StringFind(url, "?");
            if (index >=0) {
                if (index != StringLen(url) - 1) url += "&";
            } else {
                url += "?";
            }
            url += _str_data;
        }

        update_url_parts(url);

        if (!is_same_host(url) || h_session <= 0 || h_connect <= 0) open(url);

        response.url = url;
        response.parameters = parameters;

        int h_url = InternetOpenUrlW(h_session, url, nill, 0, INTERNET_FLAG_RELOAD|INTERNET_FLAG_PRAGMA_NOCACHE, 0);
        if (h_url <= 0) {
            error = "requests: ERROR InternetOpenUrlW()";
	        Print(error);
	        response.error = error;
            return(response);
        }

        read(h_url, response.text);

	    InternetCloseHandle(h_url);
        return response;
    }

    // --- POST ---
    Response post(string url, string& arr_data[][]) {
        return post(url, RequestData::to_str(arr_data));
    }

    Response post(string url, RequestData& _request_data) {
        return post(url, _request_data.to_str());
    }

    Response post(string url, string _str_data) {
        string error;
        bool is_success;

        is_success = check_dll(error);
        if (!is_success) {
            response.error = error;
            return(response);
        };

        update_url_parts(url);

        if (!is_same_host(url) || h_session <= 0 || h_connect <= 0) open(url);

        response.url = url;
        response.parameters = parameters;

        uchar data[];
        int h_request, h_send = 0;
        string method = "POST";
        string http_version = "HTTP/1.1";

        StringToCharArray(_str_data, data);

        uint flags = INTERNET_FLAG_KEEP_CONNECTION|INTERNET_FLAG_RELOAD|INTERNET_FLAG_PRAGMA_NOCACHE;
	    if (port == 443) flags |= INTERNET_FLAG_SECURE;

	    h_request = HttpOpenRequestW(h_connect, method, path, http_version, nill, 0, flags, 0);
	    if (h_request <= 0) {
	        error = "requests: ERROR HttpOpenRequestW";
	        Print(error);
	        response.error = error;
	        InternetCloseHandle(h_connect);
	        return(response);
	    }

	    string headers = "Content-Type: application/x-www-form-urlencoded";

	    int trying = 0;
        while (trying < 3) {
            trying++;
            h_send = HttpSendRequestW(h_request, headers, StringLen(headers), data, StringLen(_str_data));
            if (h_send <= 0)  {
                int err = 0; err = GetLastError(); Print("requests: ERROR HttpSendRequestW = " + (string)err + " (" + (string)trying + " trying)");
                if (err != ERROR_INTERNET_INVALID_CA) {
                    int dwFlags;
                    int dwBuffLen = sizeof(dwFlags);
                    InternetQueryOptionW(h_request, INTERNET_OPTION_SECURITY_FLAGS, dwFlags, dwBuffLen);
                    dwFlags |= SECURITY_FLAG_IGNORE_UNKNOWN_CA;
                    int res = InternetSetOptionW(h_request, INTERNET_OPTION_SECURITY_FLAGS, dwFlags, sizeof(dwFlags));
                    if (!res) { Print("requests: ERROR InternetSetOptionW = ", GetLastError()); break; }
                }
                else break;
            }
            else break;
        }

        if (h_send > 0) read(h_request, response.text);

        InternetCloseHandle(h_request);
        InternetCloseHandle(h_send);

        return response;
    }

    // --- SEND (GET or POST) ---
    Response send(string method, string url, string& arr_data[][]) {
        return send(method, url, RequestData::to_str(arr_data));
    }

    Response send(string method, string url, RequestData& _request_data) {
        return send(method, url, _request_data.to_str());
    }

    Response send(string method, string url, string _str_data) {
        StringToUpper(method);
        if (method == "POST") return post(url, _str_data);
        if (method == "GET") return get(url, _str_data);
        Print("requests: ERROR " + method + ": method is wrong. Available methods: POST, GET");

        return response;
    }

};