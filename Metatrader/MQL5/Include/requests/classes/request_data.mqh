#property copyright     "Copyright © 2021 Artem Maltsev (Vivazzi)"
#property link          "https://vivazzi.pro"
#property description   "RequestData class for requests.mqh"
#property library

/*
    Helper class for simple create request data

    USAGE:

    RequestData request_data;
    request_data.add("par_1", "foo");
    request_data.add("par_2", "bar");

    Requests requests;
    Response response = requests.get(url, request_data);
    Print("Response = " + response.text);

    ---

    You can replace value of pair using the same name in request_data.add():

    RequestData request_data;
    request_data.add("par_1", "foo");
    request_data.add("par_2", "bar");

    Print(request_data.to_str());
    // "par_1=foo&par_2=bar"

    request_data.add("par_2", "super_bar");
    Print(request_data.to_str());
    // "par_1=foo&par_2=super_bar"

    ---

    Use request_data.remove() for clear data and fill request_data new data:

    RequestData request_data;
    request_data.add("par_1", "foo");
    request_data.add("par_2", "bar");
    request_data.add("par_3", "baz");
    Print(request_data.to_str());
    // "par_1=foo&par_2=bar&par_3=baz"

    request_data.remove("par_2");  // removes data pair with specific name
    Print(request_data.to_str());
    // "par_1=foo&par_3=baz"

    request_data.remove();  // removes all data pairs
    Print(request_data.to_str());
    // ""

    ---

    Use static method to_str(string& _data[][]) if you have array of pairs:

    string array_data[2][2];
    array_data[0][0] = "par_1"; array_data[0][1] = "foo";
    array_data[1][0] = "par_2"; array_data[1][1] = "bar";
    Print(RequestData::to_str(array_data));
    // "par_1=foo&par_2=bar&par_3=baz"
*/
class RequestData {
private:
    struct RequestDataPair {
        string name;
        string value;
    };

    RequestDataPair pairs[];

public:
    void add(string name, string value) {
        int pairs_len = ArraySize(pairs);
        bool updated = false;

        for (int i=0; i < pairs_len; i++) {
            if (pairs[i].name == name) {
                pairs[i].value = value;
                updated = true;
                break;
            }
        }
        if (!updated) {
            ArrayResize(pairs, ++pairs_len);
            pairs[pairs_len-1].name = name;
            pairs[pairs_len-1].value = value;
        }
    }

    void remove() {
        ArrayFree(pairs);
    }

    void remove(string name) {
        int pairs_len = ArraySize(pairs);
        bool need_remove = false;

        for (int i=0; i < pairs_len; i++) {
            if (pairs[i].name == name) need_remove = true;

            if (need_remove) {
                if (i+1 == pairs_len) ArrayResize(pairs, --pairs_len);
                else pairs[i] = pairs[i+1];
            }
        }
    }

    string to_str() {
        string res = "";
        int pairs_len = ArraySize(pairs);

        for (int i=0; i < pairs_len; i++) {
            res += pairs[i].name + "=" + pairs[i].value;
            if (i != pairs_len - 1) res += "&";
        }

        return res;
    }

    static string to_str(string& _data[][]) {
        string res = "";
        int _data_len = ArrayRange(_data, 0);

        for (int i=0; i < _data_len; i++)
            if (_data[i][0] != NULL) {
                res += _data[i][0] + "=" + _data[i][1];
                if (i != _data_len - 1) res += "&";
            }

        return res;
    }
};