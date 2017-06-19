#include <iostream>
#include <tuple>
#include <vector>
#include <map>

using namespace std;

template<typename RType>
struct Dataset {
    typedef RType record_type;
    typedef pair<typename RType::first_type,
            vector<typename RType::second_type>> samekeyrecords_type;

    std::vector<RType> record_list;
    Dataset() { index = 0; }
    int groupbykey() {
        for (auto record : record_list) {
            indices[record.first].push_back(record.second);
        }
        indices_iter = indices.begin();
        map_index = 0;
        for (auto iter = indices.begin(); iter != indices.end(); iter++) {
            indices_vect.push_back(*iter);
        }
        return 0;
    }
    int groupbyvalue();
    // Whether has more elements
    bool hasmore() { return index < record_list.size(); }
    bool hasmore_samekeyrecords() { return map_index < indices_vect.size(); }
    // Return element one by one
    RType& next() { return record_list[index++]; }
    // return a record list with same key, groupbykey() is necessary in advance
    samekeyrecords_type next_samekeyrecords() {
        return indices_vect[map_index++];
        // return *(indices_iter++);
    }
    // Return a record list with same value, groupbyvalue() is necessary in adv
    // std::vector<RType&> next_samevaluerecords();

    int index;
    int map_index;
    map<typename RType::first_type,vector<typename RType::second_type>> indices;
    typename map<typename RType::first_type, 
        vector<typename RType::second_type>>::iterator indices_iter;
    vector<samekeyrecords_type> indices_vect;
};

template<typename T>
std::tuple<typename T::record_type> getrecord(T& t) {
    return std::make_tuple(t.next());
}

template<typename P1, typename... Param>
tuple<typename P1::record_type, typename Param::record_type...>
getrecord(P1& p1, Param&... param) {
    if (p1.hasmore())
        return tuple_cat(make_tuple(p1.next()), getrecord(param...));
    else {
        typename P1::record_type tmp_record;
        return tuple_cat(make_tuple(tmp_record), getrecord(param...));
    }
}

template<typename T>
std::tuple<typename T::samekeyrecords_type> get_samekeyrecords(T& t) {
    return std::make_tuple(t.next_samekeyrecords());
}

template<typename P1, typename... Param>
tuple<typename P1::samekeyrecords_type, typename Param::samekeyrecords_type...>
get_samekeyrecords(P1& p1, Param&... param) {
    if (p1.hasmore_samekeyrecords())
        return tuple_cat(make_tuple(p1.next_samekeyrecords()),
                get_samekeyrecords(param...));
    else {
        typename P1::samekeyrecords_type tmp_samekeyrecords;
        return tuple_cat(make_tuple(tmp_samekeyrecords),
                get_samekeyrecords(param...));
    }
}

template<typename... Ts>
bool mapfun(tuple<Ts...> args) {
    if (get<0>(args).first == "" && get<1>(args).first == "") return false;
    if (get<0>(args).first != "") {
        cout << "(" << get<0>(args).first << ", " << get<0>(args).second << ")";
    }
    if (get<1>(args).first != "") {
        cout << "(" << get<1>(args).first << ", " << get<1>(args).second << ")";
    }
    return true;
}

template<typename T, typename... Ts>
int assembler(T& t, Ts&... args) {
    bool finish = false;
    tuple<typename T::record_type, typename Ts::record_type...> records;
    do {
        records = getrecord(t, args...);
        finish = mapfun(records);
    } while (finish);
    return 0;
}

template<typename... Ts>
bool reduce(tuple<Ts...> args) {
    if (get<0>(args).first == "" && get<1>(args).first == "") return false;
    if (get<0>(args).first != "") {
        cout << "(" << get<0>(args).first << ", " << get<0>(args).second.size()
            << ")";
    }
    if (get<1>(args).first != "") {
        cout << "(" << get<1>(args).first << ", " << get<1>(args).second.size()
            << ")";
    }
    return true;
}

template<typename T, typename... Ts>
int assembler_recordlist(T& t, Ts&... args) {
    bool finish = false;
    tuple<typename T::samekeyrecords_type, typename Ts::samekeyrecords_type...>
        records;
    do {
        records = get_samekeyrecords(t, args...);
        finish = reduce(records);
    } while (finish);
    return 0;
}

int main() {
    typedef pair<string, int> WC_Type;
    Dataset<WC_Type> inds1, inds2;
    inds1.record_list = {{"the", 1}, {"focus", 1}, {"find", 1}, {"the", 1}};
    inds2.record_list = {{"an", 1}, {"focus", 1}, {"person", 1}, {"good", 1}};

    // assembler<Dataset<WC_Type>, Dataset<WC_Type>>(inds1, inds2);
    inds1.groupbykey();
    inds2.groupbykey();
    assembler_recordlist<Dataset<WC_Type>, Dataset<WC_Type>>(inds1, inds2);
    return 0;
}
