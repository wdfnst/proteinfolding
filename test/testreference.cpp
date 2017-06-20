/*************************************************************************
 > File Name: testreference.cpp
 > Author: Weidong, ZHANG
 > Mail: zhangwd@pku.edu.cn
 > Created Time: Sat May 13 16:52:04 2017
 ************************************************************************/

#include <iostream>

using namespace std;

struct A {
    string name;
    int age;
    int classno[3];
};

class B{
private:
    A &a;
public:
    B(A &a_) : a(a_) { }

    int display() {
        cout << a.name << endl;
        cout << a.age << endl;
        for (int i = 0; i < 3; ++i) {
            cout << a.classno[i] << " ";
        }
        cout << endl;
        return 0;
    }
};

int main()
{
    A a{"zhangsan" , 12, {55, 66, 77}};

    B b(a);
    b.display();

    cout << "After change:\n";
    a.classno[0] = 100;
    a.age = 28;
    b.display();
    return 0;
}
