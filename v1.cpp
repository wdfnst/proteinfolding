/*************************************************************************
 > File Name: v1.cpp
 > Author: Weidong, ZHANG
 > Mail: zhangwd@pku.edu.cn
 > Created Time: Sun May 14 19:12:04 2017
 ************************************************************************/

#include <iostream>

using namespace std;

struct S {

    typedef int value_type;
};

int main()
{
    S::value_type  value = 9;
    cout << "value:" << value << endl;
    return 0;
}
