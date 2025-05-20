#property copyright     "Copyright © 2021 Artem Maltsev (Vivazzi)"
#property link          "https://vivazzi.pro"
#property description   "_UrlParts class for requests.mqh"
#property library


class _UrlParts {
public:
    int port;
    string host;
    string path;
    string parameters;

    bool split(string url) {
        int index = StringFind(url, "://");
        if (index == -1) {
            Print("requests: ERROR " + url + ": url is wrong");
            return(false);
        }

        string protocol = StringSubstr(url, 0, index);

        if (protocol == "http") port = 80;
        else if (protocol == "https") port = 443;
        else {
            Print("requests: ERROR " + protocol + ": protocol is wrong. Available protocols: http, https");
            return(false);
        }

        string full_path = StringSubstr(url, index + 3, StringLen(url));

        index = StringFind(full_path, "/");

        if (index == -1) {
            host = full_path;
            path = "/";
        } else {
            host = StringSubstr(full_path, 0, index);
            path = StringSubstr(full_path, index, StringLen(full_path));
        }

        index = StringFind(full_path, "?");
        if (index == -1) parameters = "";
        else parameters = StringSubstr(full_path, index + 1, StringLen(full_path));

        return(true);
    }
};