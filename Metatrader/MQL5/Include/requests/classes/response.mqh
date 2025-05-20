#property copyright     "Copyright © 2021 Artem Maltsev (Vivazzi)"
#property link          "https://vivazzi.pro"
#property description   "Response class for requests.mqh"
#property library


class Response {
public:
    string text;
    string status_code;
    string error;
    string url;
    string parameters;

    Response() {};
    ~Response() {};

    Response::Response(const Response &r) {
        text = r.text;
        status_code = r.status_code;
        error = r.error;
        url = r.url;
        parameters = r.parameters;
    }
};