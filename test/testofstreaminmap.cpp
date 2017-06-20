/*************************************************************************
 > File Name: testofstreaminmap.cpp
 > Author: Weidong, ZHANG
 > Mail: zhangwd@pku.edu.cn
 > Created Time: Mon Jun 19 11:05:03 2017
 ************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <memory>

using namespace std;

class Foo {
private:
    vector<ofstream> ofsvector;
public:
    int fill(string filename) {
        ofstream ofs;
        ofs.open(filename);
        ofsvector.push_back(move(ofs));

        return 0;
    }

    int test_write(string msg) {
        vector<ofstream>::iterator iter = ofsvector.begin();
        for (; iter != ofsvector.end(); iter++) {
            iter->write(msg.c_str(), msg.size());
        }

        return 0;
    }
    ~Foo() {
        vector<ofstream>::iterator iter = ofsvector.begin();
        for (; iter != ofsvector.end(); iter++) {
            if (iter->is_open()) {
                iter->close();
            }
        }
    }
};

int main()
{
    Foo foo;
    string filename = "testfile";
    foo.fill(filename);
    foo.test_write("This is a test file msg.\n");
    return 0;
}
