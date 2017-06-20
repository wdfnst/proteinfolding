#include <cstdint>
#include <vector>
#include <iostream>
#include <string>
#include <typeinfo>

typedef std::string KeyType;
typedef int ValueType;
typedef std::pair<std::vector<KeyType>, ValueType> Record;

template<typename RType>
struct Dataset {
	std::vector<RType> record_list;
	int groupbykey();
	int groupbyvalue();
	// Return element one by one
	RType next();
	// Pre-groupbykey() is necessary
	std::vector<RType&> next_samekeylist();
	// Pre-groupbyvalue() is necessary
	std::vector<RType&> next_samevaluelist();

    typedef RType value_type;
};

template <class... Ts> struct Assembler {};
template <class T, class... Ts>
struct Assembler<T, Ts...> : Assembler<Ts...> {
    Assembler(T t, Ts... ts) : Assembler<Ts...>(ts...), tail(t) {}
    T tail;
    size_t input_size;
};
template <size_t, class> struct elem_type_holder;

template <class T, class... Ts>
struct elem_type_holder<0, Assembler<T, Ts...>> {
  typedef T type;
};

template <size_t k, class T, class... Ts>
struct elem_type_holder<k, Assembler<T, Ts...>> {
  typedef typename elem_type_holder<k - 1, Assembler<Ts...>>::type type;
};

template <size_t k, class... Ts>
typename std::enable_if<
    k == 0, typename elem_type_holder<0, Assembler<Ts...>>::type&>::type
get(Assembler<Ts...>& t) {
  return t.tail;
}

template <size_t k, class T, class... Ts>
typename std::enable_if<
    k != 0, typename elem_type_holder<k, Assembler<T, Ts...>>::type&>::type
get(Assembler<T, Ts...>& t) {
  Assembler<Ts...>& base = t;
  return get<k - 1>(base);
}

template<class... Ts>
void reduce(std::tuple<Ts...> records_2dlist) {
//     typename elem_type_holder<0, decltype(records_2dlist)>::type
//         list0 = get<0>(records_2dlist);
    // do somthing on list0
}

template<typename T>
std::tuple<T> get2dlist(const T &t) {
    return std::make_tuple(t.next_samekeylist());
}

template<typename P1, typename... Param>
std::tuple<P1, Param...> get2dlist(P1 p1, Param... param) {
    return std::tuple_cat(std::make_tuple(p1.next_samekeylist()),
            get2dlist(param...));
}

template<class... Ts>
int run(Assembler<Ts...> as) {
    std::tuple<> records_2dlist, old_records_2dlist;
    do {
        // make records_2dlist
        records_2dlist = get2dlist(as);
        for (int i = 0; i < sizeof...(Ts) + 1; ++i) {
            typedef typename elem_type_holder<0, decltype(as)>::type RType;
            std::tuple_cat(records_2dlist, std::make_tuple('a'));
        }
        // reduce(records_2dlist);
    } while (true);
    return 0;
}

int main(int argc, char** argv) {
    typedef Dataset<Record> DS_Type;
    DS_Type ds1, ds2;
    Assembler<DS_Type, DS_Type> as(ds1, ds2);
    run<DS_Type, DS_Type>(as);
    return 0;
}
